"""POST /api/v1/sources/{vmh,rhea}/fetch (US3 / FR-007). Upstreams mocked."""
from unittest.mock import patch

from django.test import TestCase

from reactions.tests.helpers import (
    NO_VMH_MATCH,
    found_all,
    make_key,
    make_user,
    post_json,
)
from reactions.utils.to_mol import smiles_to_mol

VMH_URL = '/api/v1/sources/vmh/fetch'
RHEA_URL = '/api/v1/sources/rhea/fetch'


def water_mol(*args, **kwargs):
    return smiles_to_mol('O')[0], None, 'Water'


@patch('reactions.services.reaction_processing.search_metabolites_vmh', found_all)
@patch('reactions.services.reaction_processing.check_reaction_vmh',
       return_value=NO_VMH_MATCH)
class VmhFetchTests(TestCase):
    def setUp(self):
        self.user = make_user()
        self.key = make_key(self.user)

    @patch('reactions.utils.to_mol.vmh_to_mol', side_effect=water_mol)
    @patch('reactions.services.external_sources.find_reaction_by_abbreviation',
           return_value={'formula': '1 h2o[c] -> 1 h2o[e]',
                         'subsystem': 'Transport, extracellular'})
    def test_fetch_vmh_reaction_is_processed(self, _find, _mol, _check):
        response = post_json(self.client, VMH_URL,
                             {'abbreviation': 'H2Ot'}, self.key)
        self.assertEqual(response.status_code, 201)
        body = response.json()
        self.assertTrue(body['balanced_count'])
        self.assertEqual(body['direction'], 'forward')
        self.assertEqual(body['subsystem'], 'Transport, extracellular')

    @patch('reactions.services.external_sources.vmh_api_key_configured',
           return_value=True)
    @patch('reactions.services.external_sources.find_reaction_by_abbreviation',
           return_value=None)
    def test_unknown_abbreviation(self, _find, _configured, _check):
        response = post_json(self.client, VMH_URL,
                             {'abbreviation': 'NOPE'}, self.key)
        self.assertEqual(response.status_code, 422)
        self.assertEqual(response.json()['error']['code'], 'unknown_identifier')

    @patch('reactions.services.external_sources.find_reaction_by_abbreviation',
           return_value=None)
    def test_missing_vmh_credentials_explains_setup(self, _find, _check):
        # No VMH_API_KEY in the test environment: a failed lookup must say the
        # server lacks credentials, not pretend the reaction does not exist.
        response = post_json(self.client, VMH_URL,
                             {'abbreviation': 'HEX1'}, self.key)
        self.assertEqual(response.status_code, 502)
        error = response.json()['error']
        self.assertEqual(error['code'], 'upstream_unavailable')
        self.assertIn('VMH_API_KEY', error['message'])
        self.assertIn('Rhea', error['message'])

    @patch('reactions.services.external_sources.find_reaction_by_abbreviation',
           return_value={'formula': ''})
    def test_reaction_without_formula(self, _find, _check):
        response = post_json(self.client, VMH_URL,
                             {'abbreviation': 'EMPTY'}, self.key)
        self.assertEqual(response.status_code, 422)

    @patch('reactions.services.external_sources.find_reaction_by_abbreviation',
           return_value={'formula': '1.5 h2o[c] -> 1 h2o[e]'})
    def test_fractional_stoichiometry_rejected_not_truncated(self, _find, _check):
        response = post_json(self.client, VMH_URL,
                             {'abbreviation': 'FRAC'}, self.key)
        self.assertEqual(response.status_code, 422)
        error = response.json()['error']
        self.assertEqual(error['code'], 'validation_error')
        self.assertIn('fractional', error['message'])

    @patch('reactions.services.external_sources.find_reaction_by_abbreviation',
           return_value={'formula': '500 h2o[c] -> 500 h2o[e]'})
    def test_excessive_imported_stoichiometry_rejected(self, _find, _check):
        response = post_json(self.client, VMH_URL,
                             {'abbreviation': 'BIG'}, self.key)
        self.assertEqual(response.status_code, 422)
        self.assertIn('100', response.json()['error']['message'])


@patch('reactions.services.reaction_processing.search_metabolites_vmh', found_all)
@patch('reactions.services.reaction_processing.check_reaction_vmh',
       return_value=NO_VMH_MATCH)
class RheaFetchTests(TestCase):
    def setUp(self):
        self.user = make_user()
        self.key = make_key(self.user)

    @patch('reactions.utils.to_mol.chebi_to_mol', side_effect=water_mol)
    @patch('reactions.services.external_sources.get_from_rhea',
           return_value={
               'direction': 'reversible',
               'subs_sch': [1], 'prod_sch': [1],
               'subs_types': ['ChEBI ID'], 'prods_types': ['ChEBI ID'],
               'substrates': ['CHEBI:15377'], 'products': ['CHEBI:15377'],
               'subs_comps': ['-'], 'prods_comps': ['-'],
               'subsystem': '',
           })
    def test_fetch_rhea_reaction_is_processed(self, _rhea, _mol, _check):
        response = post_json(self.client, RHEA_URL,
                             {'rhea_id': '10031'}, self.key)
        self.assertEqual(response.status_code, 201)
        body = response.json()
        self.assertTrue(body['balanced_count'])
        self.assertEqual(body['direction'], 'bidirectional')

    @patch('reactions.services.external_sources.get_from_rhea',
           return_value={'error': 'No reaction found with this abbreviation'})
    def test_unknown_rhea_id(self, _rhea, _check):
        response = post_json(self.client, RHEA_URL,
                             {'rhea_id': 'RHEA:0'}, self.key)
        self.assertEqual(response.status_code, 422)
        self.assertEqual(response.json()['error']['code'], 'unknown_identifier')

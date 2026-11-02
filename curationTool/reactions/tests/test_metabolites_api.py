"""POST /api/v1/metabolites/resolve (US2 / FR-004)."""
from unittest.mock import patch

import requests
from django.test import TestCase
from rdkit import Chem

from reactions.models import SavedMetabolite
from reactions.tests.helpers import make_key, make_user, post_json
from reactions.utils.to_mol import smiles_to_mol

URL = '/api/v1/metabolites/resolve'
NOT_IN_VMH = (False, '', '', '')


def make_saved_water(owner, name='private water'):
    mol = Chem.MolFromSmiles('O')
    return SavedMetabolite.objects.create(
        owner=owner, name=name, inchi=Chem.MolToInchi(mol),
        inchi_key=Chem.MolToInchiKey(mol), smiles='O', mol_file='')


class ResolveMetaboliteTests(TestCase):
    def setUp(self):
        self.user = make_user()
        self.key = make_key(self.user)

    def test_requires_auth(self):
        response = post_json(self.client, URL, {'identifier': 'O', 'type': 'smiles'})
        self.assertEqual(response.status_code, 401)

    @patch('reactions.services.metabolite_resolution.search_vmh',
           return_value=NOT_IN_VMH)
    def test_resolve_smiles_offline(self, _):
        response = post_json(self.client, URL,
                             {'identifier': 'O', 'type': 'smiles'}, self.key)
        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertEqual(body['formula'], 'H2O')
        self.assertEqual(body['charge'], 0)
        self.assertTrue(body['inchi_key'].startswith('XLYOFNOQVPJJNP'))
        self.assertFalse(body['vmh']['found'])

    @patch('reactions.services.metabolite_resolution.search_vmh',
           return_value=(True, 'https://vmh.example/h2o', 'h2o', 'Water'))
    @patch('reactions.utils.to_mol.vmh_to_mol',
           side_effect=lambda abbr: (smiles_to_mol('O')[0], None, 'Water'))
    def test_resolve_vmh_abbreviation(self, _mol, _search):
        response = post_json(self.client, URL,
                             {'identifier': 'h2o', 'type': 'vmh'}, self.key)
        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertEqual(body['name'], 'Water')
        self.assertEqual(body['vmh'],
                         {'found': True, 'abbreviation': 'h2o',
                          'miriam': 'https://vmh.example/h2o'})

    def test_unknown_identifier(self):
        response = post_json(self.client, URL,
                             {'identifier': 'not(((smiles', 'type': 'smiles'},
                             self.key)
        self.assertEqual(response.status_code, 422)
        self.assertEqual(response.json()['error']['code'], 'unknown_identifier')

    def test_unsupported_type_names_supported_ones(self):
        response = post_json(self.client, URL,
                             {'identifier': 'x', 'type': 'kegg'}, self.key)
        self.assertEqual(response.status_code, 422)
        error = response.json()['error']
        self.assertEqual(error['code'], 'validation_error')
        self.assertIn('smiles', error['message'])
        self.assertIn('chebi', error['message'])

    @patch('reactions.services.metabolite_resolution.any_to_mol',
           side_effect=requests.ConnectionError('boom'))
    def test_upstream_unavailable(self, _):
        response = post_json(self.client, URL,
                             {'identifier': 'CHEBI:15377', 'type': 'chebi'},
                             self.key)
        self.assertEqual(response.status_code, 502)
        self.assertEqual(response.json()['error']['code'], 'upstream_unavailable')

    def test_blank_identifier_rejected(self):
        for identifier in ('', '   '):
            response = post_json(self.client, URL,
                                 {'identifier': identifier, 'type': 'smiles'},
                                 self.key)
            self.assertEqual(response.status_code, 422)
            self.assertEqual(response.json()['error']['code'], 'validation_error')

    def test_empty_smiles_is_invalid_in_util(self):
        # RDKit parses '' as a valid zero-atom molecule; the util must not.
        mol, error, _ = smiles_to_mol('')
        self.assertIsNone(mol)
        self.assertIn('Invalid SMILES', error)

    @patch('reactions.services.metabolite_resolution.search_vmh',
           return_value=NOT_IN_VMH)
    def test_own_saved_metabolite_resolves(self, _):
        saved = make_saved_water(self.user)
        response = post_json(self.client, URL,
                             {'identifier': str(saved.pk), 'type': 'saved'},
                             self.key)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['formula'], 'H2O')

    def test_foreign_saved_metabolite_is_unknown(self):
        other = make_user(name='other')
        saved = make_saved_water(other)
        response = post_json(self.client, URL,
                             {'identifier': str(saved.pk), 'type': 'saved'},
                             self.key)
        self.assertEqual(response.status_code, 422)
        error = response.json()['error']
        self.assertEqual(error['code'], 'unknown_identifier')
        # the private name must not leak into the error
        self.assertNotIn('private water', error['message'])

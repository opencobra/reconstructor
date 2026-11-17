"""Rhea import parsing (reactions/utils/get_from_rhea.py); HTTP is mocked."""
from unittest.mock import patch

import requests
from django.test import TestCase

from reactions.utils.get_from_rhea import get_from_rhea

HEADER = 'Reaction identifier\tEquation\tChEBI identifier\n'
ROWS = {
    'RHEA:25290': 'RHEA:25290\tethanol + NAD(+) = acetaldehyde + NADH + H(+)\t'
                  'CHEBI:16236;CHEBI:57540;CHEBI:15343;CHEBI:57945;CHEBI:15378\n',
    'RHEA:20309': 'RHEA:20309\t2 H2O2 = O2 + 2 H2O\tCHEBI:16240;CHEBI:15379;CHEBI:15377\n',
    # Searching the bare number of a directional ID finds its master reaction.
    '10031': 'RHEA:10028\tD-glutamate + O2 + H2O = 2-oxoglutarate + H2O2 + NH4(+)\t'
             'CHEBI:29986;CHEBI:15379;CHEBI:15377;CHEBI:16810;CHEBI:16240;CHEBI:28938\n',
    '424242': 'RHEA:16505\tchorismate = 4-hydroxybenzoate + pyruvate\t'
              'CHEBI:29748;CHEBI:17879;CHEBI:15361\n',
    'RHEA:77777': 'RHEA:77777\tA + B = C\tCHEBI:1;CHEBI:2\n',
}


def rhea(url, params=None, timeout=None):
    assert timeout, 'Rhea requests need a timeout'
    resp = requests.Response()
    resp.status_code = 200
    # Unknown queries: header plus a blank line, as Rhea answers today.
    resp._content = (HEADER + ROWS.get(params['query'], '\n')).encode('utf-8')
    return resp


@patch('reactions.utils.get_from_rhea.requests.get', side_effect=rhea)
class RheaImportTests(TestCase):
    def test_names_containing_plus_signs(self, _):
        # 'NAD(+)' must not be split into extra components (that used to turn
        # acetaldehyde into a substrate).
        result = get_from_rhea('RHEA:25290')
        self.assertEqual(result['substrates'], ['CHEBI:16236', 'CHEBI:57540'])
        self.assertEqual(result['products'], ['CHEBI:15343', 'CHEBI:57945', 'CHEBI:15378'])
        self.assertEqual(result['subs_types'], ['ChEBI ID', 'ChEBI ID'])
        self.assertEqual(result['direction'], 'forward')

    def test_stoichiometry_is_kept(self, _):
        result = get_from_rhea('RHEA:20309')
        self.assertEqual((result['subs_sch'], result['prod_sch']), ([2], [1, 2]))

    def test_directional_id_resolves_to_master(self, _):
        result = get_from_rhea('rhea:10031')
        self.assertEqual(result['rhea_id'], 'RHEA:10028')
        self.assertEqual(len(result['substrates']), 3)

    def test_unrelated_master_is_not_accepted(self, _):
        self.assertEqual(get_from_rhea('RHEA:424242'),
                         {'error': 'No reaction found with this abbreviation'})

    def test_unknown_and_malformed(self, _):
        self.assertEqual(get_from_rhea('RHEA:1'),
                         {'error': 'No reaction found with this abbreviation'})
        self.assertIn('error', get_from_rhea('16505'))
        self.assertIn('Could not match', get_from_rhea('RHEA:77777')['error'])

    def test_rhea_outage_raises(self, _):
        failing = requests.Response()
        failing.status_code = 500
        with patch('reactions.utils.get_from_rhea.requests.get', return_value=failing):
            with self.assertRaises(requests.RequestException):
                get_from_rhea('RHEA:25290')

    def test_web_endpoint_reports_outage(self, _):
        with patch('reactions.utils.get_from_rhea.requests.get',
                   side_effect=requests.ConnectionError('down')):
            body = self.client.post('/fetch_rhea_rxn', '{"reactionAbbreviation": "RHEA:25290"}',
                                    content_type='application/json').json()
        self.assertEqual(body, {'error': 'Could not reach Rhea. Please try again shortly.'})

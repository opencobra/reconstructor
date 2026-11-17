"""ChEBI 2.0 REST client (reactions/utils/chebi_api.py); HTTP is mocked."""
import json
from unittest.mock import patch

import requests
from django.core.cache import cache
from django.test import TestCase
from rdkit import Chem

from reactions.utils import chebi_api
from reactions.utils.to_mol import chebi_to_mol
from reactions.utils.to_smiles import chebi_to_smiles

DIOXYGEN = {
    'id': 15379, 'chebi_accession': 'CHEBI:15379', 'name': 'dioxygen', 'ascii_name': 'dioxygen',
    'default_structure': {'smiles': 'O=O', 'standard_inchi': 'InChI=1S/O2/c1-2',
                          'standard_inchi_key': 'MYMOFIZGZYHOMD-UHFFFAOYSA-N'},
}
GLUCOPYRANOSE = {
    'id': 4167, 'chebi_accession': 'CHEBI:4167', 'name': '<small>D</small>-glucopyranose',
    'ascii_name': 'D-glucopyranose',
    'default_structure': {'smiles': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',
                          'standard_inchi': 'InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1',
                          'standard_inchi_key': 'WQZGKKKJIJFFOK-GASJEMHNSA-N'},
}


def http(status=200, body=None):
    resp = requests.Response()
    resp.status_code = status
    resp._content = json.dumps(body if body is not None else {}).encode('utf-8')
    resp.headers['Content-Type'] = 'application/json'
    resp.url = 'https://www.ebi.ac.uk/chebi/backend/api/public/'
    return resp


def hits(*items):
    return http(200, {'total': len(items), 'number_pages': 1,
                      'results': [{'_source': item} for item in items]})


class FakeChebi:
    def __init__(self, compounds=(), searches=None):
        self.compounds = {c['chebi_accession']: c for c in compounds}
        # ChEBI's text search is case-insensitive
        self.searches = {term.lower(): result for term, result in (searches or {}).items()}
        self.calls = []

    def __call__(self, url, params=None, timeout=None, headers=None):
        self.calls.append((url, params))
        if url.endswith('/compounds/'):
            requested = params['chebi_ids'].split(',')
            return http(200, {
                accession: ({'standardized_chebi_id': accession, 'exists': True,
                             'data': self.compounds[accession]}
                            if accession in self.compounds else
                            {'standardized_chebi_id': accession, 'exists': False, 'data': None})
                for accession in requested})
        if '/compound/' in url:
            accession = url.rstrip('/').rsplit('/', 1)[1]
            compound = self.compounds.get(accession)
            return http(200, compound) if compound else http(404, {'detail': 'Not found.'})
        if url.endswith('/es_search/'):
            return self.searches.get(params['term'].lower(), http(404, {'detail': 'Not found.'}))
        raise AssertionError(url)


class ChebiClientTests(TestCase):
    def setUp(self):
        cache.clear()

    def fake(self, **kwargs):
        fake = FakeChebi(**kwargs)
        patcher = patch('reactions.utils.chebi_api.requests.get', side_effect=fake)
        patcher.start()
        self.addCleanup(patcher.stop)
        return fake

    def test_id_forms(self):
        for value in ('CHEBI:15379', 'chebi:15379', 'CHEBI15379', ' 15379 '):
            self.assertEqual(chebi_api.normalize_chebi_id(value), 'CHEBI:15379')

    def test_compound(self):
        fake = self.fake(compounds=[DIOXYGEN])
        compound = chebi_api.get_compound('15379')
        self.assertEqual(compound, {'chebi_accession': 'CHEBI:15379', 'name': 'dioxygen',
                                    'inchi': 'InChI=1S/O2/c1-2',
                                    'inchi_key': 'MYMOFIZGZYHOMD-UHFFFAOYSA-N',
                                    'smiles': 'O=O'})
        chebi_api.get_compound('CHEBI:15379')
        self.assertEqual(len(fake.calls), 1)  # cached
        self.assertTrue(fake.calls[0][0].endswith('/compound/CHEBI:15379/'))

    def test_unknown_and_invalid(self):
        self.fake()
        self.assertIsNone(chebi_api.get_compound('CHEBI:999999999'))
        self.assertIsNone(chebi_api.get_compound('not an id'))
        self.assertIsNone(chebi_api.find_compound_by_name('zzqqnotachemical'))

    def test_outage_raises_request_exception(self):
        with patch('reactions.utils.chebi_api.requests.get', return_value=http(500)):
            with self.assertRaises(requests.RequestException):
                chebi_api.get_compound('CHEBI:15379')

    def test_exact_name_beats_ranking(self):
        # ChEBI ranks 'dioxygen(2+)' above 'dioxygen' itself.
        self.fake(searches={'dioxygen': hits(
            {'chebi_accession': 'CHEBI:29371', 'ascii_name': 'dioxygen(2+)', 'stars': 3},
            {'chebi_accession': 'CHEBI:29372', 'ascii_name': 'dioxygen(.1+)', 'stars': 3},
            {'chebi_accession': 'CHEBI:15379', 'ascii_name': 'dioxygen', 'stars': 3})})
        self.assertEqual(chebi_api.find_compound_by_name('Dioxygen')['chebi_accession'],
                         'CHEBI:15379')

    def test_synonym_and_top_hit(self):
        self.fake(searches={
            'D-glucose': hits({'chebi_accession': 'CHEBI:4167', 'ascii_name': 'D-glucopyranose',
                               'name': '<small>D</small>-glucopyranose', 'stars': 3}),
            'gluco': hits({'chebi_accession': 'CHEBI:4167', 'ascii_name': 'D-glucopyranose'},
                          {'chebi_accession': 'CHEBI:17234', 'ascii_name': 'glucose'})})
        self.assertEqual(chebi_api.find_compound_by_name('D-glucose')['chebi_accession'],
                         'CHEBI:4167')
        self.assertEqual(chebi_api.find_compound_by_name('gluco')['chebi_accession'],
                         'CHEBI:4167')

    def test_chebi_to_mol_by_id_and_name(self):
        self.fake(compounds=[DIOXYGEN, GLUCOPYRANOSE], searches={'D-glucose': hits(
            {'chebi_accession': 'CHEBI:4167', 'ascii_name': 'D-glucopyranose'})})
        mol, error, name = chebi_to_mol('CHEBI:15379', name=False)
        self.assertIsNone(error)
        self.assertEqual((Chem.MolToInchiKey(mol), name), ('MYMOFIZGZYHOMD-UHFFFAOYSA-N', 'dioxygen'))
        mol, error, name = chebi_to_mol('D-glucose', name=True)
        self.assertIsNone(error)
        self.assertEqual((Chem.MolToInchiKey(mol), name),
                         ('WQZGKKKJIJFFOK-GASJEMHNSA-N', 'D-glucopyranose'))
        smiles, error = chebi_to_smiles('CHEBI:15379', name=False)
        self.assertIsNone(error)
        self.assertEqual(Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)), 'MYMOFIZGZYHOMD-UHFFFAOYSA-N')

    def test_prefetch_uses_one_request_then_the_cache(self):
        fake = self.fake(compounds=[DIOXYGEN, GLUCOPYRANOSE])
        chebi_api.prefetch_compounds(['CHEBI:15379', '4167', 'CHEBI:999999999', 'CHEBI:15379'])
        self.assertEqual(len(fake.calls), 1)
        self.assertEqual(fake.calls[0][1]['chebi_ids'], 'CHEBI:15379,CHEBI:4167,CHEBI:999999999')
        self.assertEqual(chebi_api.get_compound('CHEBI:4167')['name'], 'D-glucopyranose')
        self.assertIsNone(chebi_api.get_compound('CHEBI:999999999'))  # cached "no such entry"
        chebi_api.prefetch_compounds(['CHEBI:15379'])  # already cached: no request
        self.assertEqual(len(fake.calls), 1)

    def test_failed_prefetch_falls_back_to_single_lookups(self):
        single = FakeChebi(compounds=[DIOXYGEN])

        def flaky(url, **kwargs):
            if url.endswith('/compounds/'):
                return http(504)
            return single(url, **kwargs)

        with patch('reactions.utils.chebi_api.requests.get', side_effect=flaky):
            with self.assertLogs('reactions.utils.chebi_api', 'WARNING'):
                chebi_api.prefetch_compounds(['CHEBI:15379'])
            self.assertEqual(chebi_api.get_compound('CHEBI:15379')['name'], 'dioxygen')

    def test_any_to_mol_batches_chebi_ids(self):
        from reactions.utils.to_mol import any_to_mol
        fake = self.fake(compounds=[DIOXYGEN, GLUCOPYRANOSE])
        mols, errors, names = any_to_mol(['CHEBI:15379', 'CHEBI:4167'], ['ChEBI ID', 'ChEBI ID'],
                                         None, side='substrates')
        self.assertEqual(errors, [None, None])
        self.assertEqual(names, ['dioxygen', 'D-glucopyranose'])
        self.assertEqual([url.rsplit('/', 2)[-2] for url, _ in fake.calls], ['compounds'])

    def test_chebi_to_mol_unknown(self):
        self.fake()
        mol, error, _ = chebi_to_mol('CHEBI:999999999', name=False)
        self.assertIsNone(mol)
        self.assertIn('does not exist in ChEBI', error)

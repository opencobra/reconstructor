"""VMH public API client (reactions/utils/vmh_api.py) and the views that use it.

VMH is faked at the single network seam (vmh_api._send), so request building,
response parsing, caching, rate-limit handling and fallbacks run for real.
"""
import json
import os
import re
import tempfile
import time
from unittest.mock import patch

import requests
from django.core.cache import caches
from django.test import TestCase, override_settings
from rdkit import Chem

from reactions.utils import vmh_api
from reactions.utils.gen_vmh_abbrs import check_met_abbr_exists, check_reaction_abbr_exists
from reactions.utils.search_vmh import check_reaction_vmh, decode_formula, search_vmh
from reactions.utils.to_mol import vmh_to_mol
from reactions.utils.vmh_api import VMHUnavailableError

WATER_KEY = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
O2_KEY = 'MYMOFIZGZYHOMD-UHFFFAOYSA-N'

O2 = {'met_id': 9, 'abbreviation': 'o2', 'fullName': 'O2', 'smile': 'O=O',
      'inchiString': 'InChI=1S/O2/c1-2', 'inchiKey': O2_KEY, 'charge': 0,
      'chargedFormula': 'O2', 'pubChemId': '977', 'synonyms': 'Dioxygen***O2'}
H2O = {'abbreviation': 'h2o', 'fullName': 'Water', 'smile': '[H]O[H]',
       'inchiString': 'InChI=1S/H2O/h1H2', 'inchiKey': WATER_KEY, 'charge': 0}
CAFFEINE_RELATIVE = {'abbreviation': '13dmxnthn', 'fullName': '1,3-Dimethylxanthine',
                     'smile': 'Cn1c(=O)c2[nH]cnc2n(C)c1=O', 'inchiString': '', 'charge': 0}
HEX1 = {'rxn_id': 12548, 'abbreviation': 'HEX1', 'description': 'Hexokinase (D-Glucose:ATP)',
        'formula': 'atp[c] + glc_D[c] -> adp[c] + g6p[c] + h[c] ',
        'subsystem': 'Glycolysis/gluconeogenesis'}


def response(status=200, body=None, *, text=None, headers=None):
    resp = requests.Response()
    resp.status_code = status
    if body is not None:
        resp._content = json.dumps(body).encode('utf-8')
        resp.headers['Content-Type'] = 'application/json; charset=utf-8'
    else:
        resp._content = (text or '').encode('utf-8')
        resp.headers['Content-Type'] = 'text/html; charset=utf-8'
    resp.encoding = 'utf-8'
    for name, value in (headers or {}).items():
        resp.headers[name] = value
    return resp


def met_handler(*records, query='', former=()):
    return response(200, {'metData': list(records), 'formerAbbrs': list(former),
                          'originalQuery': query, 'canonicalAbbr': query})


def sparql_rows(*rows):
    return response(200, {
        'head': {'vars': sorted({k for row in rows for k in row})},
        'results': {'bindings': [
            {k: {'type': 'literal', 'value': v} for k, v in row.items()} for row in rows]},
    })


def reactions(*rows):
    return response(200, {'totalRecords': len(rows), 'page': 1, 'pageLength': 500,
                          'reactions': list(rows)})


NOT_FOUND = response(404, {'error': 'No metabolite found for abbreviation "x"'})
NEXT_JS_404 = response(404, text='<!DOCTYPE html><title>Page not found - VMH</title>')


class FakeVMH:
    """Stands in for vmh_api._send: routes calls to canned responses and
    records them."""

    def __init__(self):
        self.routes = []
        self.calls = []

    def on(self, method, path, result, when=None):
        self.routes.append((method, path, when, result))

    def __call__(self, method, url, *, params=None, data=None, headers=None, timeout=None):
        assert '/api/public' in url, url
        call = {'method': method, 'url': url, 'path': url.split('/api/public', 1)[1],
                'params': dict(params or {}), 'data': dict(data or {}),
                'headers': dict(headers or {}), 'timeout': timeout}
        self.calls.append(call)
        for route_method, path, when, result in self.routes:
            if route_method == method and path == call['path'] and (when is None or when(call)):
                if callable(result):
                    result = result(call)
                if isinstance(result, BaseException):
                    raise result
                return result
        raise AssertionError(f"unexpected VMH call {method} {call['path']} "
                             f"{call['params']} {call['data']}")

    def count(self, path):
        return sum(1 for call in self.calls if call['path'] == path)


def query_has(text):
    return lambda call: text in call['data'].get('query', '')


IS_EXPORT = query_has('FILTER(!BOUND(?key))')


class VMHTestCase(TestCase):
    """Fresh caches and a FakeVMH for every test."""

    def setUp(self):
        caches['vmh'].clear()
        caches['default'].clear()
        vmh_api._index_memory.update(by_key=None, built_at=0.0, base=None)
        vmh_api._retired_host_warned.clear()
        self.vmh = FakeVMH()
        patcher = patch('reactions.utils.vmh_api._send', side_effect=self.vmh)
        patcher.start()
        self.addCleanup(patcher.stop)


class ConfigurationTests(VMHTestCase):
    @override_settings(VMH_API_BASE_URL='', VMH_BASE_URL='')
    def test_defaults_point_at_vmh_life(self):
        self.assertEqual(vmh_api.vmh_api_base_url(), 'https://www.vmh.life/api/public')
        self.assertEqual(vmh_api.vmh_metabolite_url('glc_D'),
                         'https://www.vmh.life/metabolite/glc_D')
        self.assertEqual(vmh_api.vmh_reaction_url('HEX1'),
                         'https://www.vmh.life/reaction/HEX1')

    @override_settings(VMH_API_BASE_URL='https://vmh.chatimd.org/api/public',
                       VMH_BASE_URL='https://vmh.chatimd.org')
    def test_retired_host_is_replaced(self):
        with self.assertLogs('reactions.utils.vmh_api', 'WARNING'):
            self.assertEqual(vmh_api.vmh_api_base_url(), 'https://www.vmh.life/api/public')
            self.assertEqual(vmh_api.vmh_public_base_url(), 'https://www.vmh.life')

    @override_settings(VMH_API_BASE_URL='')
    def test_requests_go_to_vmh_life_api(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'))
        vmh_api.lookup_metabolite('o2')
        self.assertEqual(self.vmh.calls[0]['url'],
                         'https://www.vmh.life/api/public/getMetHandler')
        self.assertEqual(self.vmh.calls[0]['params'], {'id': 'o2'})


class MetaboliteLookupTests(VMHTestCase):
    def test_found(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'))
        record = vmh_api.lookup_metabolite('o2')
        self.assertEqual(record['abbreviation'], 'o2')
        self.assertEqual(record['fullName'], 'O2')
        self.assertEqual(record['smile'], 'O=O')
        self.assertEqual(record['inchiString'], 'InChI=1S/O2/c1-2')
        self.assertFalse(record['is_former'])
        self.assertEqual(vmh_api.find_metabolite_by_abbreviation('o2')['abbreviation'], 'o2')

    def test_case_insensitive(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='O2'))
        self.assertEqual(vmh_api.find_metabolite_by_abbreviation('O2')['abbreviation'], 'o2')

    def test_former_abbreviation_resolves_but_is_flagged(self):
        self.vmh.on('GET', '/getMetHandler',
                    met_handler(CAFFEINE_RELATIVE, query='theopln', former=['theopln']))
        record = vmh_api.lookup_metabolite('theopln')
        self.assertTrue(record['is_former'])
        self.assertEqual(record['abbreviation'], '13dmxnthn')
        # historical helper: current abbreviations only
        self.assertIsNone(vmh_api.find_metabolite_by_abbreviation('theopln'))
        # but a former abbreviation is not free for new metabolites
        self.assertTrue(check_met_abbr_exists('theopln'))

    def test_unrelated_record_is_not_a_hit(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='9'))
        self.assertIsNone(vmh_api.lookup_metabolite('9'))

    def test_not_found_is_cached(self):
        self.vmh.on('GET', '/getMetHandler', NOT_FOUND)
        self.assertIsNone(vmh_api.lookup_metabolite('zz_nope'))
        self.assertIsNone(vmh_api.lookup_metabolite('zz_nope'))
        self.assertEqual(self.vmh.count('/getMetHandler'), 1)
        self.assertFalse(check_met_abbr_exists('zz_nope'))

    def test_answers_are_cached(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'))
        for _ in range(3):
            vmh_api.lookup_metabolite('o2')
        self.assertEqual(self.vmh.count('/getMetHandler'), 1)

    def test_non_json_404_means_the_api_moved(self):
        # What the retired /_api/ paths return today: a Next.js HTML page.
        self.vmh.on('GET', '/getMetHandler', NEXT_JS_404)
        with self.assertLogs('reactions.utils.vmh_api', 'ERROR'):
            with self.assertRaises(VMHUnavailableError) as ctx:
                vmh_api.lookup_metabolite('o2')
        self.assertIn('endpoint was not found', str(ctx.exception))
        with self.assertLogs('reactions.utils.vmh_api', 'ERROR'):
            with self.assertRaises(VMHUnavailableError):
                vmh_api.lookup_metabolite('o2')
        self.assertEqual(self.vmh.count('/getMetHandler'), 2)  # failures are not cached

    def test_server_error_and_network_failures_are_unavailable(self):
        for failure in (response(503, {'error': 'down'}),
                        requests.ConnectionError('refused'),
                        requests.Timeout('slow'),
                        response(200, text='<html>maintenance</html>')):
            caches['vmh'].clear()
            self.vmh.routes = []
            self.vmh.on('GET', '/getMetHandler', failure)
            with self.assertLogs('reactions.utils.vmh_api'):
                with self.assertRaises(VMHUnavailableError):
                    vmh_api.lookup_metabolite('o2')
                # never-raise helpers degrade to "not found" (and log)
                self.assertIsNone(vmh_api.find_metabolite_by_abbreviation('o2'))

    def test_full_name_lookup_is_exact_and_case_insensitive(self):
        self.vmh.on('GET', '/advancedSearch/metabolite', response(200, {
            'totalRecords': 2, 'page': 1, 'pageLength': 10000, 'metabolites': [
                {'abbreviation': 'glc_D6p', 'fullName': 'D-glucose 6-phosphate'},
                {'abbreviation': 'glc_D', 'fullName': 'D-glucose'}]}),
            when=lambda call: call['params'].get('searchTerm') == 'D-Glucose')
        self.assertEqual(vmh_api.lookup_metabolite_by_name('D-Glucose')['abbreviation'], 'glc_D')
        self.assertTrue(vmh_api.metabolite_name_exists('D-Glucose'))


class RateLimitAndKeyTests(VMHTestCase):
    def test_429_pauses_all_requests_until_reset(self):
        reset = int(time.time()) + 600
        self.vmh.on('GET', '/getMetHandler', response(
            429, {'error': 'Too many requests'},
            headers={'x-ratelimit-limit': '100', 'x-ratelimit-remaining': '0',
                     'x-ratelimit-reset': str(reset)}))
        with self.assertLogs('reactions.utils.vmh_api', 'WARNING'):
            with self.assertRaises(VMHUnavailableError) as ctx:
                vmh_api.lookup_metabolite('o2')
        exc = ctx.exception
        self.assertTrue(exc.rate_limited)
        self.assertTrue(590 <= exc.retry_after <= 601, exc.retry_after)
        self.assertIn('10 minutes', exc.user_message())
        self.assertIn('VMH_API_KEY', exc.user_message())  # no key configured in tests
        # Any further call fails fast without contacting VMH.
        with self.assertRaises(VMHUnavailableError):
            vmh_api.lookup_reaction('HEX1')
        self.assertEqual(len(self.vmh.calls), 1)

    @override_settings(VMH_API_KEY='secret-key')
    def test_key_is_sent_when_configured(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'))
        vmh_api.lookup_metabolite('o2')
        self.assertEqual(self.vmh.calls[0]['headers'].get('x-api-key'), 'secret-key')

    def test_no_key_header_without_key(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'))
        vmh_api.lookup_metabolite('o2')
        self.assertNotIn('x-api-key', self.vmh.calls[0]['headers'])

    @override_settings(VMH_API_KEY='revoked-key')
    def test_rejected_key_is_dropped_and_request_retried(self):
        self.vmh.on('GET', '/getMetHandler', response(401, {'error': 'bad key'}),
                    when=lambda call: 'x-api-key' in call['headers'])
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'))
        with self.assertLogs('reactions.utils.vmh_api', 'ERROR'):
            self.assertEqual(vmh_api.lookup_metabolite('o2')['abbreviation'], 'o2')
        self.vmh.on('GET', '/getRxnHandler', response(200, {'rxnData': [HEX1]}))
        vmh_api.lookup_reaction('HEX1')
        self.assertNotIn('x-api-key', self.vmh.calls[-1]['headers'])

    @override_settings(VMH_API_KEY='unknown-key')
    def test_key_treated_as_anonymous_is_reported(self):
        # VMH answers an unrecognised key as an anonymous (rate-limited) client.
        self.vmh.on('GET', '/getMetHandler', response(
            200, {'metData': [O2]}, headers={'x-ratelimit-limit': '100',
                                             'x-ratelimit-remaining': '95'}))
        with self.assertLogs('reactions.utils.vmh_api', 'WARNING') as logs:
            vmh_api.lookup_metabolite('o2')
        self.assertIn('anonymous', '\n'.join(logs.output))


class StructureLookupTests(VMHTestCase):
    def test_inchikey_via_sparql_post(self):
        self.vmh.on('POST', '/sparql', sparql_rows({'abbr': 'o2', 'name': 'O2'}),
                    when=query_has(f'vmh:inchiKey "{O2_KEY}"'))
        record = vmh_api.lookup_metabolite_by_inchikey(O2_KEY.lower())
        self.assertEqual((record['abbreviation'], record['fullName']), ('o2', 'O2'))
        self.assertEqual(self.vmh.calls[0]['data']['format'], 'json')

    def test_several_matches_pick_the_first_abbreviation(self):
        self.vmh.on('POST', '/sparql', sparql_rows({'abbr': 'h2o_b'}, {'abbr': 'h2o'},
                                                   {'abbr': 'h2o'}))
        self.assertEqual(vmh_api.lookup_metabolite_by_inchikey(WATER_KEY)['abbreviation'],
                         'h2o')

    def test_sparql_failure_falls_back_to_regex_search(self):
        self.vmh.on('POST', '/sparql', response(500, {'error': 'Internal server error.'}))
        self.vmh.on('GET', '/advancedSearch/metabolite', response(200, {
            'metabolites': [{'abbreviation': 'o2', 'fullName': 'O2', 'inchiKey': O2_KEY},
                            {'abbreviation': 'o2x', 'fullName': 'other',
                             'inchiKey': 'MYMOFIZGZYHOMD-UHFFFAOYSA-X'}]}),
            when=lambda call: call['params'].get('inchiKeyRegex') == f'^{O2_KEY}$')
        with self.assertLogs('reactions.utils.vmh_api', 'WARNING'):
            self.assertEqual(vmh_api.lookup_metabolite_by_inchikey(O2_KEY)['abbreviation'],
                             'o2')
        fallback = [c for c in self.vmh.calls if c['path'] == '/advancedSearch/metabolite'][0]
        self.assertGreaterEqual(fallback['timeout'], 60)  # the regex filter is slow

    def test_rate_limited_sparql_does_not_fall_back(self):
        self.vmh.on('POST', '/sparql', response(429, {}, headers={'Retry-After': '30'}))
        with self.assertLogs('reactions.utils.vmh_api', 'WARNING'):
            with self.assertRaises(VMHUnavailableError):
                vmh_api.lookup_metabolite_by_inchikey(O2_KEY)
        self.assertEqual(self.vmh.count('/advancedSearch/metabolite'), 0)

    def test_malformed_inchikey_makes_no_request(self):
        self.assertIsNone(vmh_api.lookup_metabolite_by_inchikey('not-a-key'))
        self.assertEqual(self.vmh.calls, [])

    def test_smiles_literal_is_escaped(self):
        smiles = r'C\C=C/C'
        self.vmh.on('POST', '/sparql', sparql_rows({'abbr': 'but2en'}))
        self.assertEqual(vmh_api.lookup_metabolite_by_smiles(smiles)['abbreviation'], 'but2en')
        self.assertIn(r'vmh:smiles "C\\C=C/C"', self.vmh.calls[0]['data']['query'])

    def test_inchikey_computed_from_inchi_when_missing(self):
        self.vmh.on('POST', '/sparql', sparql_rows({'abbr': 'h2o'}),
                    when=query_has(WATER_KEY))
        row = vmh_api.match_metabolite_structure(inchi_string='InChI=1S/H2O/h1H2')
        self.assertEqual(row['abbreviation'], 'h2o')

    def test_search_vmh_finds_molecule(self):
        self.vmh.on('POST', '/sparql', sparql_rows({'abbr': 'h2o', 'name': 'Water'}),
                    when=query_has(WATER_KEY))
        found, miriam, abbr, name = search_vmh(Chem.MolFromSmiles('O'),
                                               return_abbr=True, return_name=True)
        self.assertEqual((found, abbr, name), (True, 'h2o', 'Water'))
        self.assertEqual(miriam, 'https://www.vmh.life/metabolite/h2o')

    def test_search_vmh_strict_raises_but_default_degrades(self):
        self.vmh.on('POST', '/sparql', response(503, {}))
        self.vmh.on('GET', '/advancedSearch/metabolite', response(503, {}))
        mol = Chem.MolFromSmiles('O')
        with self.assertLogs('reactions.utils.vmh_api'):
            self.assertEqual(search_vmh(mol), (False, ''))
            with self.assertRaises(VMHUnavailableError):
                search_vmh(mol, strict=True)


@override_settings(VMH_DERIVED_INDEX=True, VMH_API_KEY='test-key')
class DerivedIndexTests(VMHTestCase):
    """VMH stores many metabolites with a SMILES but no InChIKey."""

    # VMH's export of metabolites without an InChIKey (water's key is derived
    # from its SMILES; the R-group structure yields none).
    EXPORT = sparql_rows({'abbr': 'h2o', 'smiles': '[H]O[H]'},
                         {'abbr': 'h2o', 'smiles': 'O'},
                         {'abbr': 'rgroup', 'smiles': '[*]C(=O)O'})

    def setUp(self):
        super().setUp()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        override = override_settings(VMH_CACHE_DIR=self.tmp.name)
        override.enable()
        self.addCleanup(override.disable)
        # Direct InChIKey lookups find nothing: VMH has no key for these.
        self.vmh.on('POST', '/sparql', sparql_rows(), when=query_has('vmh:inchiKey "'))
        self.vmh.on('POST', '/sparql', sparql_rows(), when=query_has('vmh:smiles "'))
        self.vmh.on('GET', '/getMetHandler', met_handler(H2O, query='h2o'))

    def test_keyless_metabolite_found_via_derived_key(self):
        self.vmh.on('POST', '/sparql', self.EXPORT,
                    when=IS_EXPORT)
        record = vmh_api.lookup_metabolite_by_inchikey(WATER_KEY)
        self.assertEqual(record['abbreviation'], 'h2o')
        # the index is built once and reused (from memory, then from disk)
        vmh_api.lookup_metabolite_by_inchikey('AAAAAAAAAAAAAA-BBBBBBBBBB-N')
        vmh_api._index_memory.update(by_key=None, built_at=0.0, base=None)
        vmh_api.lookup_metabolite_by_inchikey('CCCCCCCCCCCCCC-DDDDDDDDDD-N')
        exports = [c for c in self.vmh.calls if IS_EXPORT(c)]
        self.assertEqual(len(exports), 1)
        self.assertTrue(os.path.exists(vmh_api._index_path()))

    @override_settings(VMH_API_KEY='')
    def test_not_built_without_api_key(self):
        # Anonymous SPARQL answers are truncated: a partial index would look complete.
        self.assertIsNone(vmh_api.lookup_metabolite_by_inchikey(WATER_KEY))
        self.assertFalse(any(IS_EXPORT(call) for call in self.vmh.calls))

    def test_unavailable_index_does_not_cache_not_found(self):
        self.vmh.on('POST', '/sparql', response(503, {}),
                    when=IS_EXPORT)
        with self.assertLogs('reactions.utils.vmh_api', 'WARNING'):
            self.assertIsNone(vmh_api.lookup_metabolite_by_inchikey(WATER_KEY))
        # later, with the index available, the same key is found
        caches['vmh'].delete(vmh_api._state_key('index-build-failed'))
        self.vmh.routes = [r for r in self.vmh.routes if r[3].status_code != 503]
        self.vmh.on('POST', '/sparql', self.EXPORT,
                    when=IS_EXPORT)
        self.assertEqual(vmh_api.lookup_metabolite_by_inchikey(WATER_KEY)['abbreviation'],
                         'h2o')

    def test_failed_build_is_not_retried_immediately(self):
        self.vmh.on('POST', '/sparql', response(503, {}),
                    when=IS_EXPORT)
        with self.assertLogs('reactions.utils.vmh_api', 'WARNING'):
            vmh_api.lookup_metabolite_by_inchikey(WATER_KEY)
            vmh_api.lookup_metabolite_by_inchikey(O2_KEY)
        exports = [c for c in self.vmh.calls if IS_EXPORT(c)]
        self.assertEqual(len(exports), 1)


class ReactionLookupTests(VMHTestCase):
    def test_found(self):
        self.vmh.on('GET', '/getRxnHandler', response(200, {'rxnData': [HEX1]}))
        record = vmh_api.lookup_reaction('HEX1')
        self.assertEqual(record['formula'].strip(), 'atp[c] + glc_D[c] -> adp[c] + g6p[c] + h[c]')
        self.assertFalse(record['is_former'])

    def test_other_letter_case(self):
        # getRxnHandler is case-sensitive; the retired API was not.
        self.vmh.on('GET', '/getRxnHandler', response(200, {'rxnData': [HEX1]}),
                    when=lambda call: call['params'] == {'id': 'HEX1'})
        self.vmh.on('GET', '/getRxnHandler', response(404, {'error': 'No matching reaction found.'}))
        self.assertEqual(vmh_api.find_reaction_by_abbreviation('hex1')['abbreviation'], 'HEX1')
        self.assertTrue(check_reaction_abbr_exists('hex1'))
        self.assertEqual([c['params']['id'] for c in self.vmh.calls], ['hex1', 'HEX1'])

    def test_former_abbreviation(self):
        self.vmh.on('GET', '/getRxnHandler', response(200, {
            'rxnData': [HEX1], 'formerRxnAbbrs': ['HEX1old'], 'redirectAbbr': 'HEX1'}))
        self.assertIsNone(vmh_api.find_reaction_by_abbreviation('HEX1old'))
        self.assertEqual(vmh_api.find_reaction_by_abbreviation(
            'HEX1old', allow_former=True)['abbreviation'], 'HEX1')

    def test_unknown(self):
        self.vmh.on('GET', '/getRxnHandler', response(404, {'error': 'No matching reaction found.'}))
        self.assertIsNone(vmh_api.lookup_reaction('Nope1'))
        self.assertFalse(check_reaction_abbr_exists('Nope1'))
        # as given, upper- and lower-case; then cached
        self.assertEqual([c['params']['id'] for c in self.vmh.calls], ['Nope1', 'NOPE1', 'nope1'])

    def test_description_exists(self):
        self.vmh.on('GET', '/advancedSearch/reaction', reactions(HEX1))
        self.assertTrue(vmh_api.reaction_description_exists('hexokinase (d-glucose:atp)'))
        self.assertFalse(vmh_api.reaction_description_exists('Hexokinase'))

    def test_existence_checks_raise_when_unavailable(self):
        # They gate submissions to VMH: "unknown" must never read as "free".
        self.vmh.on('GET', '/getRxnHandler', requests.ConnectionError('down'))
        self.vmh.on('GET', '/getMetHandler', requests.ConnectionError('down'))
        with self.assertLogs('reactions.utils.vmh_api'):
            with self.assertRaises(VMHUnavailableError):
                check_reaction_abbr_exists('NEWRXN')
            with self.assertRaises(VMHUnavailableError):
                check_met_abbr_exists('newmet')


class FormulaRegexTests(TestCase):
    CASES = [
        (['atp', 'glc_D'], ['adp', 'g6p', 'h'], 'atp[c] + glc_D[c] -> adp[c] + g6p[c] + h[c] ', True),
        (['glc_D', 'atp'], ['h', 'g6p', 'adp'], 'atp[c] + glc_D[c] -> h[c] + adp[c] + g6p[c]', True),
        (['atp', 'glc_D'], ['adp', 'g6p', 'h'], 'adp[c] + g6p[c] + h[c] -> atp[c] + glc_D[c]', False),
        # ABC transporters contain atp/h2o/adp/pi/h plus a cargo: wrong term count
        (['atp', 'h2o'], ['adp', 'pi', 'h'],
         'atp[c] + h2o[c] + chol[e] -> adp[c] + pi[c] + h[c] + chol[c]', False),
        (['h', 'r_o2'], ['h2o2', 'o2'], '2.0 h[c] + 2.0 r_o2[c] -> h2o2[c] + o2[c] ', True),
        (['h2o'], ['h2o'], 'h2o[e]  <=> h2o[c]', True),
        (['h2o'], [], 'h2o[e] <=>', True),
        (['h2o'], [], 'h2o[e] <=> h2o[c]', False),
        (['h'], ['o2'], 'h2o2[c] -> o2[c]', False),  # 'h' is not 'h2o2'
    ]

    def test_sides_and_term_counts(self):
        for substrates, products, formula, expected in self.CASES:
            pattern = vmh_api.build_formula_regex(substrates, products)
            self.assertEqual(bool(re.search(pattern, formula)), expected,
                             f'{substrates} -> {products} vs {formula!r}')

    def test_special_characters_are_escaped(self):
        pattern = vmh_api.build_formula_regex(['a.b'], ['c(d)'])
        self.assertIn(r'a\.b\[', pattern)
        self.assertIn(r'c\(d\)\[', pattern)
        self.assertIsNone(re.search(pattern, 'axb[c] -> c(d)[c]'))


class DecodeFormulaTests(TestCase):
    def test_stoichiometry_compartments_direction(self):
        self.assertEqual(
            decode_formula('2.0 h[c] + 2.0 r_o2[c] -> h2o2[c] + o2[c] '),
            (['h', 'r_o2'], ['h2o2', 'o2'], [2.0, 2.0], [1.0, 1.0],
             ['c', 'c'], ['c', 'c'], 'forward'))

    def test_exchange_and_irregular_spacing(self):
        self.assertEqual(decode_formula('h2o[e] <=>')[1], [])
        self.assertEqual(decode_formula('h2o[e]  <=> h2o[p]')[:2], (['h2o'], ['h2o']))
        self.assertEqual(decode_formula('h2o[e]  <=> h2o[p]')[6], 'bidirectional')

    def test_no_arrow(self):
        with self.assertRaises(ValueError):
            decode_formula('h2o[e] + o2[c]')


class CheckReactionVmhTests(VMHTestCase):
    def check(self, subsystem='Glycolysis/gluconeogenesis', subs_comps=None,
              subs_sch=(1, 1), types=('VMH', 'VMH')):
        mol = Chem.MolFromSmiles('O')
        return check_reaction_vmh(
            ['atp', 'glc_D'], ['adp', 'g6p', 'h'], list(subs_sch), [1, 1, 1],
            list(types), ['VMH', 'VMH', 'VMH'], [mol, mol], [mol, mol, mol],
            'forward', subsystem, subs_comps or ['c', 'c'], ['c', 'c', 'c'])

    def test_exact_match_uses_one_request(self):
        self.vmh.on('GET', '/advancedSearch/reaction', reactions(HEX1))
        result = self.check()
        self.assertEqual(result, {'found': True, 'similar': False,
                                  'url': 'https://www.vmh.life/reaction/HEX1',
                                  'formula': HEX1['formula']})
        self.assertEqual(len(self.vmh.calls), 1)
        self.assertIn('formulaRegex', self.vmh.calls[0]['params'])

    def test_different_subsystem_or_compartment_is_similar(self):
        self.vmh.on('GET', '/advancedSearch/reaction', reactions(HEX1))
        self.assertTrue(self.check(subsystem='Other')['similar'])
        self.assertTrue(self.check(subs_comps=['m', 'c'])['similar'])

    def test_exact_match_preferred_over_earlier_similar_one(self):
        mito = dict(HEX1, abbreviation='HEX1m', formula='atp[m] + glc_D[m] -> adp[m] + g6p[m] + h[m]')
        self.vmh.on('GET', '/advancedSearch/reaction', reactions(mito, HEX1))
        result = self.check()
        self.assertFalse(result['similar'])
        self.assertEqual(result['url'], 'https://www.vmh.life/reaction/HEX1')

    def test_stoichiometry_mismatch_is_not_found(self):
        self.vmh.on('GET', '/advancedSearch/reaction', reactions(HEX1))
        self.assertFalse(self.check(subs_sch=(2, 1))['found'])

    def test_metabolite_not_in_vmh_skips_reaction_search(self):
        self.vmh.on('POST', '/sparql', sparql_rows())  # structure not in VMH
        result = self.check(types=('SMILES', 'VMH'))
        self.assertFalse(result['found'])
        self.assertEqual(self.vmh.count('/advancedSearch/reaction'), 0)

    def test_vmh_unavailable_reports_error_not_crash(self):
        self.vmh.on('GET', '/advancedSearch/reaction', response(502, {}))
        with self.assertLogs('reactions.utils.vmh_api'):
            result = self.check()
        self.assertFalse(result['found'])
        self.assertIn('Could not check the reaction against VMH', result['error'])


class GeneLookupTests(VMHTestCase):
    GENE = {'geneData': [{'symbol': 'GCK', 'name': 'glucokinase', 'entrez_id': '2645'}],
            'rxnData': [HEX1]}

    def test_symbol(self):
        self.vmh.on('GET', '/getGeneHandler', response(200, self.GENE),
                    when=lambda call: call['params'] == {'id': 'GCK'})
        rows = vmh_api.gene_data_new('GCK')
        self.assertEqual((rows[0]['symbol'], rows[0]['gene_number']), ('GCK', '2645'))
        self.assertEqual(vmh_api.gene_reaction_formulas('GCK'), [HEX1['formula']])

    def test_entrez_id_is_mapped_through_sparql(self):
        self.vmh.on('POST', '/sparql', sparql_rows({'symbol': 'GCK'}),
                    when=query_has('vmh:entrezId "2645"'))
        self.vmh.on('GET', '/getGeneHandler', response(200, self.GENE),
                    when=lambda call: call['params'] == {'id': 'GCK'})
        self.assertEqual(vmh_api.gene_data_new('2645.1')[0]['symbol'], 'GCK')

    def test_unknown_gene(self):
        self.vmh.on('GET', '/getGeneHandler', response(404, {'error': 'Gene not found.'}))
        self.assertEqual(vmh_api.gene_data_new('NOTAGENE'), [])
        self.assertFalse(vmh_api.gene_exists_new('NOTAGENE'))


class VmhToMolTests(VMHTestCase):
    def test_vmh_structure(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'))
        mol, error, name = vmh_to_mol('o2')
        self.assertIsNone(error)
        self.assertEqual(Chem.MolToInchiKey(mol), O2_KEY)
        self.assertEqual(name, 'O2')

    def test_unknown(self):
        self.vmh.on('GET', '/getMetHandler', NOT_FOUND)
        mol, error, _ = vmh_to_mol('zz_nope')
        self.assertIsNone(mol)
        self.assertIn('does not exist in VMH', error)

    def test_renamed(self):
        self.vmh.on('GET', '/getMetHandler',
                    met_handler(CAFFEINE_RELATIVE, query='theopln', former=['theopln']))
        _, error, _ = vmh_to_mol('theopln')
        self.assertIn('use `13dmxnthn`', error)

    def test_unavailable_raises_unless_a_local_mol_file_exists(self):
        self.vmh.on('GET', '/getMetHandler', requests.ConnectionError('down'))
        with self.assertLogs('reactions.utils.vmh_api'):
            with self.assertRaises(VMHUnavailableError):
                vmh_to_mol('o2')
        with tempfile.TemporaryDirectory() as directory:
            with open(os.path.join(directory, 'o2.mol'), 'w', encoding='utf-8') as handle:
                handle.write(Chem.MolToMolBlock(Chem.MolFromSmiles('O=O')))
            with override_settings(MOL_FILE_PATH=directory), \
                    self.assertLogs('reactions.utils.vmh_api'):
                mol, error, _ = vmh_to_mol('o2')
        self.assertIsNone(error)
        self.assertEqual(Chem.MolToInchiKey(mol), O2_KEY)


class VerifyMetaboliteViewTests(VMHTestCase):
    URL = '/verify_metabolite/'

    def verify(self, metabolite, mtype='VMH', user_id='null'):
        return self.client.post(self.URL, {'metabolite': metabolite, 'type': mtype,
                                           'userID': user_id})

    def test_vmh_metabolite_found(self):
        # The reported bug: every metabolite came back "not found in VMH".
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'))
        response_ = self.verify('o2')
        self.assertEqual(response_.status_code, 200)
        body = response_.json()
        self.assertTrue(body['found'])
        self.assertEqual(body['abbr'], 'o2')
        self.assertEqual(body['miriam'], 'https://www.vmh.life/metabolite/o2')
        self.assertEqual(body['atom_counts'], {'O': 2})
        self.assertEqual(body['charge'], 0)

    def test_input_is_trimmed(self):
        self.vmh.on('GET', '/getMetHandler', met_handler(O2, query='o2'),
                    when=lambda call: call['params'] == {'id': 'o2'})
        self.assertTrue(self.verify('  o2 ').json()['found'])

    def test_former_abbreviation_switches_to_current(self):
        self.vmh.on('GET', '/getMetHandler',
                    met_handler(CAFFEINE_RELATIVE, query='theopln', former=['theopln']))
        body = self.verify('theopln').json()
        self.assertTrue(body['found'])
        self.assertEqual(body['abbr'], '13dmxnthn')
        self.assertIn('former VMH abbreviation', body['notice'])

    def test_unknown_metabolite(self):
        self.vmh.on('GET', '/getMetHandler', NOT_FOUND)
        response_ = self.verify('zz_nope')
        self.assertEqual(response_.status_code, 404)
        self.assertEqual(response_.json()['message'], 'Metabolite `zz_nope` not found in VMH.')

    def test_vmh_unavailable_is_not_reported_as_not_found(self):
        self.vmh.on('GET', '/getMetHandler', response(429, {}, headers={'Retry-After': '120'}))
        with self.assertLogs('reactions.utils.vmh_api'):
            response_ = self.verify('o2')
        self.assertEqual(response_.status_code, 503)
        body = response_.json()
        self.assertTrue(body['error'])
        self.assertTrue(body['vmh_unavailable'])
        self.assertNotIn('not found', body['message'])
        self.assertIn('2 minutes', body['message'])

    def test_structure_found_in_vmh_when_signed_out(self):
        # Signed-out browsers post userID="null"; that used to raise a 500.
        self.vmh.on('POST', '/sparql', sparql_rows({'abbr': 'h2o', 'name': 'Water'}),
                    when=query_has(WATER_KEY))
        body = self.verify('O', mtype='SMILES').json()
        self.assertTrue(body['found'])
        self.assertEqual((body['abbr'], body['name']), ('h2o', 'Water'))
        self.assertNotIn('warning', body)

    def test_structure_check_unavailable_warns_but_resolves(self):
        self.vmh.on('POST', '/sparql', response(503, {}))
        self.vmh.on('GET', '/advancedSearch/metabolite', response(503, {}))
        with self.assertLogs('reactions.utils.vmh_api'):
            response_ = self.verify('O', mtype='SMILES')
        self.assertEqual(response_.status_code, 200)
        body = response_.json()
        self.assertFalse(body['found'])
        self.assertEqual(body['atom_counts'], {'H': 2, 'O': 1})
        self.assertTrue(body['vmh_unavailable'])
        self.assertIn('treated as a new metabolite', body['warning'])

    def test_get_and_empty_requests(self):
        self.assertEqual(self.client.get(self.URL).status_code, 405)
        self.assertTrue(self.client.post(self.URL, {'metabolite': ' ', 'type': 'VMH'}).json()['error'])
        self.assertTrue(self.client.post(self.URL, {'metabolite': 'o2'}).json()['error'])
        self.assertEqual(self.vmh.calls, [])


class GetFromVmhViewTests(VMHTestCase):
    URL = '/get_from_vmh/'

    def fetch(self, abbreviation):
        return self.client.post(self.URL, json.dumps({'reactionAbbreviation': abbreviation}),
                                content_type='application/json').json()

    def test_import(self):
        self.vmh.on('GET', '/getRxnHandler', response(200, {'rxnData': [HEX1]}))
        body = self.fetch('HEX1')
        self.assertEqual(body['substrates'], ['atp', 'glc_D'])
        self.assertEqual(body['products'], ['adp', 'g6p', 'h'])
        self.assertEqual(body['subs_types'], ['VMH', 'VMH'])
        self.assertEqual(body['direction'], 'forward')
        self.assertEqual(body['subsystem'], 'Glycolysis/gluconeogenesis')

    def test_unknown(self):
        self.vmh.on('GET', '/getRxnHandler', response(404, {'error': 'No matching reaction found.'}))
        self.assertEqual(self.fetch('NOPE')['error'], 'No reaction found with this abbreviation')

    def test_unavailable(self):
        self.vmh.on('GET', '/getRxnHandler', requests.Timeout('slow'))
        with self.assertLogs('reactions.utils.vmh_api'):
            error = self.fetch('HEX1')['error']
        self.assertIn('Could not fetch `HEX1` from VMH', error)

    def test_bad_requests(self):
        self.assertEqual(self.client.get(self.URL).status_code, 400)
        self.assertEqual(self.client.post(self.URL, 'not json',
                                          content_type='application/json').status_code, 400)


class VmhAvailabilityViewTests(VMHTestCase):
    URL = '/check_vmh_availability/'

    def check(self, **payload):
        return self.client.post(self.URL, json.dumps(payload), content_type='application/json')

    def test_reaction_name_and_abbreviation_taken(self):
        self.vmh.on('GET', '/getRxnHandler', response(200, {'rxnData': [HEX1]}))
        self.vmh.on('GET', '/advancedSearch/reaction', reactions(HEX1))
        body = self.check(type='reaction', name='Hexokinase (D-Glucose:ATP)', abbr='HEX1').json()
        self.assertTrue(body['name_in_vmh'])
        self.assertTrue(body['abbr_in_vmh'])

    def test_metabolite_available(self):
        self.vmh.on('GET', '/getMetHandler', NOT_FOUND)
        self.vmh.on('GET', '/advancedSearch/metabolite', response(200, {'metabolites': []}))
        body = self.check(type='metabolite', name='My new metabolite', abbr='mynewmet').json()
        self.assertTrue(body['name_ok'])
        self.assertTrue(body['abbr_ok'])

    def test_unavailable(self):
        self.vmh.on('GET', '/getMetHandler', response(500, {}))
        with self.assertLogs('reactions.utils.vmh_api'):
            response_ = self.check(type='metabolite', name='x', abbr='y')
        self.assertEqual(response_.status_code, 503)
        self.assertEqual(response_.json()['status'], 'error')

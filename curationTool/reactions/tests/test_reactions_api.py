"""POST/GET /api/v1/reactions (US1 + US3 list/get; FR-005/FR-006/SC-002)."""
import json
from unittest.mock import patch

from django.test import TestCase

from reactions.models import Reaction, SavedMetabolite
from reactions.tests.helpers import (
    NO_VMH_MATCH,
    found_all,
    found_none,
    make_key,
    make_user,
    post_json,
    smiles_met,
)

URL = '/api/v1/reactions'


def water_transport():
    return {
        'substrates': [smiles_met('O', 'water', compartment='c')],
        'products': [smiles_met('O', 'water', compartment='e')],
        'direction': 'forward',
    }


@patch('reactions.services.reaction_processing.search_metabolites_vmh', found_all)
@patch('reactions.services.reaction_processing.check_reaction_vmh',
       return_value=NO_VMH_MATCH)
class CreateReactionTests(TestCase):
    def setUp(self):
        self.user = make_user()
        self.key = make_key(self.user)

    def test_balanced_water_transport(self, _):
        response = post_json(self.client, URL, water_transport(), self.key)
        self.assertEqual(response.status_code, 201)
        body = response.json()
        self.assertTrue(body['balanced_count'])
        self.assertTrue(body['balanced_charge'])
        self.assertEqual(body['subs_atoms'], {'O': 1, 'H': 2})
        self.assertEqual(body['metabolite_formulas'], ['H2O', 'H2O'])
        self.assertFalse(body['vmh']['found'])
        self.assertIsInstance(body['reaction_id'], int)

    def test_verdict_matrix(self, _):
        """SC-002: representative reactions get the expected verdicts (offline)."""
        cases = [
            # (substrates, products, mass_balanced, charge_balanced)
            ([('O', 1)], [('O', 1)], True, True),                      # transport
            ([('O', 1)], [('OO', 1)], False, True),                    # H2O -> H2O2
            ([('O', 1)], [('[OH-]', 1)], False, False),                # charge loss
            ([('O', 2)], [('OO', 1), ('[H][H]', 1)], True, True),      # 2 H2O -> H2O2 + H2
            ([('[Na+]', 1), ('[Cl-]', 1)], [('[Na+]', 1), ('[Cl-]', 1)], True, True),
        ]
        for i, (subs, prods, mass_ok, charge_ok) in enumerate(cases):
            payload = {
                'substrates': [smiles_met(s, f'sub{i}{j}', stoichiometry=n)
                               for j, (s, n) in enumerate(subs)],
                'products': [smiles_met(p, f'prod{i}{j}', stoichiometry=n)
                             for j, (p, n) in enumerate(prods)],
                'direction': 'forward',
            }
            body = post_json(self.client, URL, payload, self.key).json()
            self.assertEqual(bool(body['balanced_count']), mass_ok,
                             f'mass verdict case {i}: {body["balanced_count"]}')
            self.assertEqual(bool(body['balanced_charge']), charge_ok,
                             f'charge verdict case {i}: {body["balanced_charge"]}')

    def test_unknown_identifier_names_side_and_index(self, _):
        payload = water_transport()
        payload['substrates'][0]['identifier'] = 'bad((('
        response = post_json(self.client, URL, payload, self.key)
        self.assertEqual(response.status_code, 422)
        error = response.json()['error']
        self.assertEqual(error['code'], 'unknown_identifier')
        self.assertIn('substrates[0]', error['message'])

    def test_validation_rejects_empty_reaction(self, _):
        response = post_json(self.client, URL,
                             {'substrates': [], 'products': []}, self.key)
        self.assertEqual(response.status_code, 422)

    def test_validation_rejects_bad_compartment(self, _):
        payload = water_transport()
        payload['substrates'][0]['compartment'] = 'zz'
        response = post_json(self.client, URL, payload, self.key)
        self.assertEqual(response.status_code, 422)
        self.assertIn('compartment', response.json()['error']['message'])

    def test_validation_rejects_excessive_stoichiometry(self, _):
        payload = water_transport()
        payload['substrates'][0]['stoichiometry'] = 101
        response = post_json(self.client, URL, payload, self.key)
        self.assertEqual(response.status_code, 422)
        self.assertIn('100', response.json()['error']['message'])

    def test_foreign_saved_metabolite_rejected(self, _):
        other = make_user(name='other')
        theirs = SavedMetabolite.objects.create(
            owner=other, name='their secret', inchi='InChI=1S/H2O/h1H2',
            inchi_key='XLYOFNOQVPJJNP-UHFFFAOYSA-N', smiles='O', mol_file='')
        payload = water_transport()
        payload['substrates'][0] = {'identifier': str(theirs.pk), 'type': 'saved'}
        response = post_json(self.client, URL, payload, self.key)
        self.assertEqual(response.status_code, 422)
        error = response.json()['error']
        self.assertEqual(error['code'], 'unknown_identifier')
        self.assertIn('substrates[0]', error['message'])
        self.assertNotIn('their secret', error['message'])
        # and nothing was copied into the caller's collection
        self.assertEqual(SavedMetabolite.objects.filter(owner=self.user).count(), 0)

    def test_vmh_check_error_string_degrades_to_not_found(self, mock_check):
        # check_reaction_vmh returns a json-encoded error STRING when VMH is
        # unreachable; creation must still succeed with vmh.found == False.
        mock_check.return_value = json.dumps(
            {'error': 'Failed to fetch data from API. Status code: no_response'})
        response = post_json(self.client, URL, water_transport(), self.key)
        self.assertEqual(response.status_code, 201)
        body = response.json()
        self.assertTrue(body['balanced_count'])
        self.assertEqual(body['vmh'],
                         {'found': False, 'similar': False,
                          'url': None, 'formula': None})

    def test_row_encoding_matches_web_format(self, _):
        body = post_json(self.client, URL, water_transport(), self.key).json()
        reaction = Reaction.objects.get(pk=body['reaction_id'])
        self.assertEqual(json.loads(reaction.balanced_count), [True])
        self.assertEqual(json.loads(reaction.substrates_types), ['SMILES'])
        self.assertEqual(json.loads(reaction.subs_comps), ['c'])
        self.assertEqual(json.loads(reaction.substrates_names), ['water'])
        self.assertEqual(json.loads(reaction.subs_sch), [1])
        self.assertEqual(reaction.direction, 'forward')
        self.assertTrue(reaction.reaction_signature)


@patch('reactions.services.reaction_processing.check_reaction_vmh',
       return_value=NO_VMH_MATCH)
class SavedMetabolitePersistenceTests(TestCase):
    def setUp(self):
        self.user = make_user()
        self.key = make_key(self.user)

    @patch('reactions.services.reaction_processing.search_metabolites_vmh', found_none)
    def test_non_vmh_metabolites_become_saved(self, _):
        body = post_json(self.client, URL, water_transport(), self.key).json()
        self.assertEqual(
            SavedMetabolite.objects.filter(owner=self.user, name='water').count(), 1)
        reaction = Reaction.objects.get(pk=body['reaction_id'])
        self.assertEqual(json.loads(reaction.substrates_types), ['Saved'])

    @patch('reactions.services.reaction_processing.search_metabolites_vmh', found_none)
    def test_duplicate_metabolite_name_rejected(self, _):
        post_json(self.client, URL, water_transport(), self.key)
        payload = {
            'substrates': [smiles_met('OCC', 'water')],  # different molecule, same name
            'products': [smiles_met('OCC', 'ethanol', compartment='e')],
            'direction': 'forward',
        }
        response = post_json(self.client, URL, payload, self.key)
        self.assertEqual(response.status_code, 409)
        self.assertEqual(response.json()['error']['code'], 'duplicate_name')
        # failed processing must not leave an orphan row behind
        self.assertEqual(
            Reaction.objects.exclude(
                pk__in=self.user.created_reactions.values('reaction')).count(), 0)

    @patch('reactions.services.reaction_processing.search_metabolites_vmh', found_none)
    def test_failed_create_rolls_back_saved_metabolites(self, _):
        post_json(self.client, URL, water_transport(), self.key)  # saves 'water'
        payload = {
            # new metabolite that would persist before the failure...
            'substrates': [smiles_met('CCO', 'ethanol')],
            # ...then a name clash with the existing 'water' aborts the request
            'products': [smiles_met('CC', 'water', compartment='e')],
            'direction': 'forward',
        }
        response = post_json(self.client, URL, payload, self.key)
        self.assertEqual(response.status_code, 409)
        self.assertFalse(
            SavedMetabolite.objects.filter(owner=self.user, name='ethanol').exists())
        self.assertEqual(
            Reaction.objects.exclude(
                pk__in=self.user.created_reactions.values('reaction')).count(), 0)


@patch('reactions.services.reaction_processing.search_metabolites_vmh', found_all)
@patch('reactions.services.reaction_processing.check_reaction_vmh',
       return_value=NO_VMH_MATCH)
class SaveListGetTests(TestCase):
    def setUp(self):
        self.user = make_user()
        self.key = make_key(self.user)
        self.other = make_user(name='other')
        self.other_key = make_key(self.other)

    def _create(self, key=None):
        return post_json(self.client, URL, water_transport(),
                         key or self.key).json()['reaction_id']

    def test_save_list_get_roundtrip(self, _):
        reaction_id = self._create()
        response = post_json(self.client, f'{URL}/{reaction_id}/save',
                             {'name': 'water transport'}, self.key)
        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.json()['saved'])
        self.assertTrue(
            self.user.saved_reactions.filter(pk=reaction_id).exists())

        listing = self.client.get(URL, HTTP_AUTHORIZATION=f'Bearer {self.key}').json()
        self.assertEqual(listing['count'], 1)
        item = listing['results'][0]
        self.assertEqual(item['reaction_id'], reaction_id)
        self.assertEqual(item['name'], 'water transport')
        self.assertTrue(item['balanced_mass'])

        detail = self.client.get(f'{URL}/{reaction_id}',
                                 HTTP_AUTHORIZATION=f'Bearer {self.key}')
        self.assertEqual(detail.status_code, 200)
        self.assertTrue(detail.json()['balanced_count'])

    def test_duplicate_reaction_name_rejected(self, _):
        first = self._create()
        second = self._create()
        post_json(self.client, f'{URL}/{first}/save', {'name': 'rxn'}, self.key)
        response = post_json(self.client, f'{URL}/{second}/save',
                             {'name': 'rxn'}, self.key)
        self.assertEqual(response.status_code, 409)

    def test_foreign_reaction_is_not_found(self, _):
        reaction_id = self._create()
        post_json(self.client, f'{URL}/{reaction_id}/save',
                  {'name': 'mine'}, self.key)
        for url in (f'{URL}/{reaction_id}', f'{URL}/999999'):
            response = self.client.get(
                url, HTTP_AUTHORIZATION=f'Bearer {self.other_key}')
            self.assertEqual(response.status_code, 404)
            self.assertEqual(response.json()['error']['code'], 'not_found')
        # saving someone else's reaction is also not possible
        response = post_json(self.client, f'{URL}/{reaction_id}/save',
                             {'name': 'stolen'}, self.other_key)
        self.assertEqual(response.status_code, 404)

    def test_list_only_shows_own(self, _):
        mine = self._create()
        post_json(self.client, f'{URL}/{mine}/save', {'name': 'mine'}, self.key)
        theirs = self._create(self.other_key)
        post_json(self.client, f'{URL}/{theirs}/save', {'name': 'theirs'},
                  self.other_key)
        listing = self.client.get(URL, HTTP_AUTHORIZATION=f'Bearer {self.key}').json()
        self.assertEqual([r['reaction_id'] for r in listing['results']], [mine])

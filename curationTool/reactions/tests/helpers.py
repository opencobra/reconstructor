"""Shared fixtures for the API v1 test suites. External sources (VMH, ChEBI,
Rhea, PubChem) are always mocked; chemistry is computed offline via RDKit."""
import json

from reactions.models import User
from reactions.services import api_keys

# check_reaction_vmh result meaning "not in VMH" (network call mocked out).
NO_VMH_MATCH = {'found': False, 'similar': False, 'url': None, 'formula': None}


def make_user(name='tester', password='pw'):
    return User.objects.create(name=name, password=password)


def make_key(user):
    return api_keys.generate_key(user)


def auth(key):
    return {'HTTP_AUTHORIZATION': f'Bearer {key}'}


def post_json(client, url, payload, key=None):
    extra = auth(key) if key else {}
    return client.post(url, data=json.dumps(payload),
                       content_type='application/json', **extra)


def smiles_met(smiles, name, stoichiometry=1, compartment='c'):
    return {'identifier': smiles, 'type': 'smiles',
            'stoichiometry': stoichiometry, 'compartment': compartment,
            'name': name}


def found_all(mols, types, request, side='substrates', **kwargs):
    """search_metabolites_vmh stub: everything already in VMH."""
    return [True] * len(mols), ['https://vmh.example/m'] * len(mols)


def found_none(mols, types, request, side='substrates', **kwargs):
    """search_metabolites_vmh stub: nothing in VMH (exercises SavedMetabolite path)."""
    return [False] * len(mols), [''] * len(mols)

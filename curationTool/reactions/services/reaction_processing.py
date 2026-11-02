"""Reaction processing for the v1 API (FR-005/FR-006).

Runs the same domain pipeline as the web view ``input_reaction``
(reactions/views/reaction_views.py) — any_to_mol → get_mol_info →
construct_reaction_rxnfile → get_reaction_info → check_reaction_vmh →
search_metabolites_vmh → SavedMetabolite persistence — and stores Reaction rows
with the identical JSON field encoding, so API-created reactions render in the
web UI. Web-only concerns (file uploads, drawings, RDT atom mapping, template
rendering) are deliberately absent; atom mapping is always skipped here, like
the web flow's skipAtomMapping path.
"""
import json
import logging

import requests
from django.db.models import Q
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from reactions.models import CreatedReaction, Reaction, SavedMetabolite
from reactions.reaction_info import get_reaction_info
from reactions.services.constants import API_TO_INTERNAL_TYPE, COMPARTMENTS
from reactions.services.errors import (
    DuplicateNameError,
    NotFoundError,
    ProcessingError,
    ServiceValidationError,
    UnknownIdentifierError,
    UnsupportedTypeError,
    UpstreamServiceError,
)
from reactions.utils.get_mol_info import get_mol_info
from reactions.utils.process_strings import construct_reaction_rxnfile
from reactions.utils.search_vmh import check_reaction_vmh, search_metabolites_vmh
from reactions.utils.to_mol import any_to_mol, smiles_with_explicit_hydrogens
from reactions.utils.utils import safe_json_loads

SKIP_ATOM_MAPPING_IMAGE = '/images/atom_mapping_skip.png'

logger = logging.getLogger(__name__)


def _normalize_vmh_check(raw):
    """check_reaction_vmh has three return shapes: a full dict, a bare
    {'found': False}, and a json-encoded error string when VMH is unreachable.
    Degrade the error case to 'not found' — the same policy vmh_api applies to
    metabolite lookups — and fill missing keys."""
    if isinstance(raw, str):
        try:
            parsed = json.loads(raw)
        except ValueError:
            parsed = {'error': raw}
        logger.warning('VMH reaction check unavailable: %s',
                       parsed.get('error', raw))
        raw = {}
    return {
        'found': raw.get('found', False),
        'similar': raw.get('similar', False),
        'url': raw.get('url'),
        'formula': raw.get('formula'),
    }


def _normalize_side(items, side):
    """Validate one side's metabolite dicts and split them into parallel lists."""
    identifiers, types, stoichs, comps, names = [], [], [], [], []
    for index, item in enumerate(items):
        api_type = item.get('type')
        internal = API_TO_INTERNAL_TYPE.get(api_type)
        if internal is None:
            raise UnsupportedTypeError(
                f"Unsupported type '{api_type}' for {side}[{index}].",
                detail={'side': side, 'index': index})
        try:
            stoich = int(item.get('stoichiometry', 1))
        except (TypeError, ValueError):
            raise ServiceValidationError(
                f'Stoichiometry values must be whole numbers ({side}[{index}]).',
                detail={'side': side, 'index': index})
        if stoich < 1:
            raise ServiceValidationError(
                f'Stoichiometry must be >= 1 ({side}[{index}]).',
                detail={'side': side, 'index': index})
        comp = item.get('compartment', 'c')
        if comp not in COMPARTMENTS:
            raise ServiceValidationError(
                f"Unknown compartment '{comp}' for {side}[{index}]. "
                f"Accepted compartments: {', '.join(COMPARTMENTS)}.",
                detail={'side': side, 'index': index})
        identifiers.append(str(item.get('identifier', '')).strip())
        types.append(internal)
        stoichs.append(stoich)
        comps.append(comp)
        names.append((item.get('name') or '').strip())
    return identifiers, types, stoichs, comps, names


def _resolve_side(identifiers, types, side):
    try:
        mols, errors, auto_names = any_to_mol(identifiers, types, None, side=side)
    except requests.RequestException as exc:
        raise UpstreamServiceError(
            f'An upstream metabolite source is unavailable: {exc}')
    except Exception as exc:
        raise ProcessingError(f'Failed to process metabolites: {exc}')
    problems = []
    for index, error in enumerate(errors):
        if error is not None or mols[index] is None:
            problems.append({
                'side': side,
                'index': index,
                'identifier': identifiers[index],
                'message': error or f"Could not resolve '{identifiers[index]}'",
            })
    if problems:
        raise UnknownIdentifierError(
            '; '.join(f"{p['side']}[{p['index']}] '{p['identifier']}': {p['message']}"
                      for p in problems),
            detail={'problems': problems})
    return mols, auto_names


def _persist_saved_metabolites(mols, found_list, names, types, identifiers, user):
    """Mirror of input_reaction's save_metabolites: store non-VMH metabolites as
    the user's SavedMetabolites and rewrite type/identifier to 'Saved'/pk."""
    for i, (mol, found) in enumerate(zip(mols, found_list)):
        if found:
            continue
        mol_smiles = Chem.MolToSmiles(mol, allHsExplicit=True)
        mol_smiles = smiles_with_explicit_hydrogens(mol_smiles)
        mol = Chem.MolFromSmiles(mol_smiles)
        mol = Chem.AddHs(mol, explicitOnly=True)
        inchi_key = MolToInchiKey(mol)
        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol)
        mol_w = MolWt(mol)
        mol_formula = CalcMolFormula(mol)
        try:
            if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                raise ValueError('Embedding failed for molecule.')
            AllChem.UFFOptimizeMolecule(mol)
            mol_file = Chem.MolToMolBlock(mol)
        except Exception:
            mol_file = Chem.MolToMolBlock(mol)
        existing = SavedMetabolite.objects.filter(
            Q(owner=user) | Q(shared_with=user), inchi_key=inchi_key).first()
        if existing is None:
            if SavedMetabolite.objects.filter(
                    Q(owner=user) | Q(shared_with=user),
                    name=names[i].strip()).exists():
                raise DuplicateNameError(
                    f'Metabolite with name {names[i]} already exists in your '
                    'saved metabolites. Choose a different name or use the '
                    'existing metabolite.')
            existing = SavedMetabolite.objects.create(
                owner=user,
                name=names[i].strip(),
                inchi_key=inchi_key,
                inchi=inchi,
                smiles=smiles,
                mol_w=mol_w,
                mol_formula=mol_formula,
                mol_file=mol_file,
                source_type=types[i],
                original_identifier=types[i],
            )
        types[i] = 'Saved'
        identifiers[i] = existing.pk


def process_reaction(substrates, products, direction='forward', subsystem=None,
                     organs=None, user=None):
    """Process and persist a reaction submission. Returns (reaction, response dict).

    ``substrates``/``products`` are lists of dicts with keys identifier, type
    (API values), stoichiometry, compartment, name.
    """
    organs = organs or []
    if not substrates and not products:
        raise ServiceValidationError(
            'Enter at least one substrate or one product before creating a reaction.')
    if user is None and any(
            (m.get('type') == 'saved') for m in list(substrates) + list(products)):
        raise ServiceValidationError(
            'Cannot use saved metabolites without signing in.')

    subs_ids, subs_types, subs_sch, subs_comp, subs_names = _normalize_side(
        substrates, 'substrates')
    prod_ids, prod_types, prod_sch, prod_comp, prod_names = _normalize_side(
        products, 'products')

    # Same storage values the web flow uses ('forward' / 'bidirectional').
    internal_direction = 'forward' if direction == 'forward' else 'bidirectional'

    subs_mols, subs_auto = _resolve_side(subs_ids, subs_types, 'substrates')
    prod_mols, prod_auto = _resolve_side(prod_ids, prod_types, 'products')

    subs_names = [n or a or str(i) for n, a, i in zip(subs_names, subs_auto, subs_ids)]
    prod_names = [n or a or str(i) for n, a, i in zip(prod_names, prod_auto, prod_ids)]

    mol_data = get_mol_info(subs_mols + prod_mols)
    if len(mol_data['smiles']) != len(subs_mols) + len(prod_mols):
        raise ProcessingError(
            'Failed to compute molecular information for one or more metabolites.')

    rxn_file = construct_reaction_rxnfile(
        subs_mols, subs_sch, prod_mols, prod_sch, subs_names, prod_names)

    subs_inchi_keys = sorted(
        (MolToInchiKey(mol), int(stoich)) for mol, stoich in zip(subs_mols, subs_sch))
    prods_inchi_keys = sorted(
        (MolToInchiKey(mol), int(stoich)) for mol, stoich in zip(prod_mols, prod_sch))
    reaction_signature = json.dumps(
        {'substrates': subs_inchi_keys, 'products': prods_inchi_keys}, sort_keys=True)

    reaction = Reaction()
    reaction.save()
    try:
        reaction_info = get_reaction_info(rxn_file, internal_direction)

        try:
            vmh_found = _normalize_vmh_check(check_reaction_vmh(
                subs_ids, prod_ids, subs_sch, prod_sch, subs_types, prod_types,
                subs_mols, prod_mols, internal_direction, subsystem,
                subs_comp, prod_comp))
        except requests.RequestException as exc:
            raise UpstreamServiceError(f'VMH is unavailable: {exc}')

        try:
            subs_found, subs_miriams = search_metabolites_vmh(
                subs_ids, subs_types, None, side='substrates')
            prod_found, prod_miriams = search_metabolites_vmh(
                prod_ids, prod_types, None, side='products')
        except requests.RequestException as exc:
            raise UpstreamServiceError(f'VMH is unavailable: {exc}')

        if user:
            _persist_saved_metabolites(
                subs_mols, subs_found, subs_names, subs_types, subs_ids, user)
            _persist_saved_metabolites(
                prod_mols, prod_found, prod_names, prod_types, prod_ids, user)

        reaction.Organs = json.dumps(organs)
        reaction.subs_sch = json.dumps(subs_sch)
        reaction.prods_sch = json.dumps(prod_sch)
        reaction.substrates_types = json.dumps(subs_types)
        reaction.products_types = json.dumps(prod_types)
        reaction.subs_comps = json.dumps(subs_comp)
        reaction.prods_comps = json.dumps(prod_comp)
        reaction.substrates_names = json.dumps(subs_names)
        reaction.products_names = json.dumps(prod_names)
        reaction.substrates = json.dumps(subs_ids)
        reaction.products = json.dumps(prod_ids)
        reaction.direction = internal_direction
        reaction.subsystem = subsystem
        reaction.visualization = json.dumps([SKIP_ATOM_MAPPING_IMAGE])
        reaction.molc_formula = json.dumps([reaction_info['molc_formula']])
        reaction.balanced_count = json.dumps([reaction_info['balanced_count']])
        reaction.balanced_charge = json.dumps([reaction_info['balanced_charge']])
        reaction.subs_atoms = json.dumps([reaction_info['subs_atoms']])
        reaction.prods_atoms = json.dumps([reaction_info['prods_atoms']])
        reaction.subs_charge = json.dumps([reaction_info['subs_charge']])
        reaction.prods_charge = json.dumps([reaction_info['prods_charge']])
        reaction.symb_to_name = json.dumps([reaction_info['symb_to_name']])
        reaction.subs_found = json.dumps(subs_found)
        reaction.subs_miriams = json.dumps(subs_miriams)
        reaction.prod_found = json.dumps(prod_found)
        reaction.prod_miriams = json.dumps(prod_miriams)
        reaction.vmh_found = vmh_found['found']
        reaction.vmh_found_similar = vmh_found['similar']
        reaction.vmh_url = json.dumps(
            vmh_found['url']) if vmh_found['found'] else None
        reaction.vmh_formula = json.dumps(
            vmh_found['formula']) if vmh_found['found'] else None
        reaction.metabolite_names = json.dumps(subs_names + prod_names)
        reaction.metabolite_formulas = json.dumps(mol_data['metabolite_formulas'])
        reaction.metabolite_smiles = json.dumps(mol_data['smiles'])
        reaction.metabolite_inchis = json.dumps(mol_data['inchis'])
        reaction.metabolite_inchi_keys = json.dumps(mol_data['inchi_keys'])
        reaction.metabolite_mol_weights = json.dumps(mol_data['mol_weights'])
        reaction.metabolite_charges = json.dumps(mol_data['metabolite_charges'])
        reaction.metabolite_mol_file_strings = json.dumps(
            mol_data['metabolite_mol_file_strings'])
        reaction.stereo_counts = json.dumps(mol_data['stereo_counts'])
        reaction.stereo_locations_list = json.dumps(
            mol_data['stereo_locations_list'])
        reaction.reaction_signature = reaction_signature
        reaction.save()
    except Exception:
        reaction.delete()
        raise

    if user:
        # Same bookkeeping the web UI does via its create_reaction endpoint, and
        # what makes the new reaction visible to this user's GET/save calls.
        CreatedReaction.objects.create(user=user, reaction=reaction)

    return reaction, serialize_reaction(reaction)


def serialize_reaction(reaction):
    """Shape a stored Reaction into the ProcessedReaction response dict."""
    def first(field):
        values = safe_json_loads(field)
        if isinstance(values, list) and values:
            return values[0]
        return values

    return {
        'reaction_id': reaction.id,
        'direction': reaction.direction,
        'subsystem': reaction.subsystem,
        'molc_formula': first(reaction.molc_formula),
        'balanced_count': first(reaction.balanced_count),
        'balanced_charge': first(reaction.balanced_charge),
        'subs_atoms': first(reaction.subs_atoms),
        'prods_atoms': first(reaction.prods_atoms),
        'subs_charge': first(reaction.subs_charge),
        'prods_charge': first(reaction.prods_charge),
        'symb_to_name': first(reaction.symb_to_name),
        'metabolite_names': safe_json_loads(reaction.metabolite_names) or [],
        'metabolite_formulas': safe_json_loads(reaction.metabolite_formulas) or [],
        'metabolite_charges': safe_json_loads(reaction.metabolite_charges) or [],
        'metabolite_smiles': safe_json_loads(reaction.metabolite_smiles) or [],
        'metabolite_inchi_keys': safe_json_loads(reaction.metabolite_inchi_keys) or [],
        'metabolite_mol_weights': safe_json_loads(reaction.metabolite_mol_weights) or [],
        'subs_found': safe_json_loads(reaction.subs_found) or [],
        'subs_miriams': safe_json_loads(reaction.subs_miriams) or [],
        'prod_found': safe_json_loads(reaction.prod_found) or [],
        'prod_miriams': safe_json_loads(reaction.prod_miriams) or [],
        'vmh': {
            'found': reaction.vmh_found,
            'similar': reaction.vmh_found_similar,
            'url': safe_json_loads(reaction.vmh_url) if reaction.vmh_url else None,
            'formula': safe_json_loads(reaction.vmh_formula) if reaction.vmh_formula else None,
        },
    }


def get_user_reaction(user, reaction_id):
    """Fetch a reaction visible to the user (saved or created), else NotFoundError."""
    reaction = Reaction.objects.filter(pk=reaction_id).filter(
        Q(saved_by_users=user) | Q(created_by_users__user=user)).distinct().first()
    if reaction is None:
        raise NotFoundError('Reaction not found.')
    return reaction


def save_reaction_for_user(user, reaction_id, name, description=''):
    """Mirror of the web save flow: set name/description and add to the user's
    saved set. Only reactions visible to the user (created or already saved by
    them) can be saved; duplicate names (per user) are rejected."""
    reaction = get_user_reaction(user, reaction_id)
    if user.saved_reactions.filter(short_name=name).exclude(pk=reaction.pk).exists():
        raise DuplicateNameError(
            f"You already have a saved reaction named '{name}'. "
            'Choose a different name.')
    reaction.short_name = name
    reaction.description = description or ''
    reaction.save()
    user.saved_reactions.add(reaction)
    return reaction

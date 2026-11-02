"""Resolve a single metabolite identifier to structure + VMH presence (FR-004)."""
import requests
from rdkit import Chem

from reactions.services.constants import (
    API_TO_INTERNAL_TYPE,
    SUPPORTED_TYPES,
    VMH_NOT_CONFIGURED_HINT,
)
from reactions.services.errors import (
    ServiceValidationError,
    UnknownIdentifierError,
    UnsupportedTypeError,
    UpstreamServiceError,
)
from reactions.services.saved_metabolites import visible_to_user
from reactions.utils.get_mol_info import get_mol_info
from reactions.utils.search_vmh import search_vmh
from reactions.utils.to_mol import any_to_mol
from reactions.utils.vmh_api import vmh_api_key_configured


def resolve_metabolite(identifier, id_type, user=None):
    """Return a dict matching the MetaboliteResolution schema.

    Raises UnsupportedTypeError / UnknownIdentifierError / UpstreamServiceError.
    """
    internal_type = API_TO_INTERNAL_TYPE.get(id_type)
    if internal_type is None:
        raise UnsupportedTypeError(
            f"Unsupported identifier type '{id_type}'. "
            f"Supported types: {', '.join(SUPPORTED_TYPES)}.")

    identifier = str(identifier).strip()
    if not identifier:
        raise ServiceValidationError('identifier must not be empty.')

    # Saved-metabolite pks are global in the resolver utils; only ids the
    # caller owns (or was shared) may be dereferenced here.
    if internal_type == 'Saved' and not visible_to_user(identifier, user):
        raise UnknownIdentifierError(
            f"No saved metabolite with id '{identifier}' in your collection.",
            detail={'identifier': identifier, 'type': id_type})

    try:
        mols, errors, names = any_to_mol(
            [identifier], [internal_type], None, side='substrates')
    except requests.RequestException as exc:
        raise UpstreamServiceError(
            f'Upstream source for type {id_type} is unavailable: {exc}')
    except Exception as exc:
        raise UnknownIdentifierError(
            f"Could not resolve identifier '{identifier}' (type {id_type}): {exc}")

    if errors[0] or mols[0] is None:
        if internal_type == 'VMH' and not vmh_api_key_configured():
            raise UpstreamServiceError(
                f"Cannot resolve VMH identifier '{identifier}': "
                f'{VMH_NOT_CONFIGURED_HINT}')
        raise UnknownIdentifierError(
            errors[0] or f"Could not resolve identifier '{identifier}' (type {id_type}).",
            detail={'identifier': identifier, 'type': id_type})

    mol = mols[0]
    if mol.GetNumAtoms() == 0:
        raise UnknownIdentifierError(
            f"Identifier '{identifier}' resolved to an empty structure.",
            detail={'identifier': identifier, 'type': id_type})
    info = get_mol_info([mol])

    try:
        found, miriam, abbr, vmh_name = search_vmh(
            mol, return_abbr=True, return_name=True)
    except requests.RequestException as exc:
        raise UpstreamServiceError(f'VMH is unavailable: {exc}')

    charge = info['metabolite_charges'][0]
    return {
        'name': names[0] or vmh_name or None,
        'formula': info['metabolite_formulas'][0],
        'charge': int(charge) if charge is not None else None,
        'smiles': info['smiles'][0],
        'inchi': info['inchis'][0],
        'inchi_key': info['inchi_keys'][0],
        'mol_weight': info['mol_weights'][0],
        'vmh': {
            'found': bool(found),
            'abbreviation': abbr or None,
            'miriam': miriam or None,
        },
    }

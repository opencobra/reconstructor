"""Import reactions from VMH / Rhea and run them through the processing
pipeline (FR-007). Reuses the same helpers as the web fetch endpoints."""
import requests

from reactions.services.constants import VMH_NOT_CONFIGURED_HINT
from reactions.services.errors import (
    ServiceValidationError,
    UnknownIdentifierError,
    UpstreamServiceError,
)
from reactions.services.reaction_processing import process_reaction
from reactions.utils.get_from_rhea import get_from_rhea
from reactions.utils.search_vmh import decode_formula
from reactions.utils.vmh_api import (
    find_reaction_by_abbreviation,
    vmh_api_key_configured,
)

# Internal type strings -> API type values (inverse of API_TO_INTERNAL_TYPE).
_INTERNAL_TO_API = {
    'VMH': 'vmh',
    'ChEBI ID': 'chebi',
    'ChEBI Name': 'chebi_name',
    'PubChem ID': 'pubchem',
    'SwissLipids': 'swisslipids',
    'SMILES': 'smiles',
    'Saved': 'saved',
}


def _coefficient(stoich, identifier):
    """Whole-number coefficient, or a clear error: silently truncating a
    fractional upstream coefficient (e.g. '1.5 h2o[c]') would change the
    imported reaction and its balance verdicts."""
    value = float(stoich)
    if not value.is_integer():
        raise ServiceValidationError(
            f"Stoichiometry {value} for '{identifier}' is fractional; only "
            'whole-number coefficients can be imported.')
    return int(value)


def _side(identifiers, types, stoichs, comps):
    return [
        {
            'identifier': identifier,
            'type': _INTERNAL_TO_API.get(mtype, mtype),
            'stoichiometry': _coefficient(stoich, identifier),
            'compartment': comp,
            'name': None,
        }
        for identifier, mtype, stoich, comp in zip(identifiers, types, stoichs, comps)
    ]


def fetch_vmh_reaction(abbreviation, user):
    """Fetch a VMH reaction by abbreviation and process it for the user."""
    try:
        reaction_data = find_reaction_by_abbreviation(abbreviation)
    except requests.RequestException as exc:
        raise UpstreamServiceError(f'VMH is unavailable: {exc}')
    if not reaction_data:
        if not vmh_api_key_configured():
            # Without credentials every lookup comes back empty; saying
            # "not found" would be misleading.
            raise UpstreamServiceError(
                f"Could not fetch '{abbreviation}' from VMH: "
                f'{VMH_NOT_CONFIGURED_HINT}')
        raise UnknownIdentifierError(
            f"No VMH reaction found with abbreviation '{abbreviation}'.")
    formula = reaction_data.get('formula', '')
    if not formula:
        raise UnknownIdentifierError(
            f"Reaction '{abbreviation}' has no formula in VMH.")
    (substrates, products, subs_sch, prod_sch,
     subs_comps, prod_comps, direction) = decode_formula(formula)
    api_direction = 'forward' if direction == 'forward' else 'reversible'
    return process_reaction(
        substrates=_side(substrates, ['VMH'] * len(substrates), subs_sch, subs_comps),
        products=_side(products, ['VMH'] * len(products), prod_sch, prod_comps),
        direction=api_direction,
        subsystem=reaction_data.get('subsystem') or None,
        user=user,
    )


def fetch_rhea_reaction(rhea_id, user):
    """Fetch a Rhea reaction by id (with or without the RHEA: prefix) and
    process it for the user."""
    rhea_id = str(rhea_id).strip()
    if not rhea_id.upper().startswith('RHEA:'):
        rhea_id = f'RHEA:{rhea_id}'
    try:
        result = get_from_rhea(rhea_id)
    except requests.RequestException as exc:
        raise UpstreamServiceError(f'Rhea is unavailable: {exc}')
    if 'error' in result:
        raise UnknownIdentifierError(f"Rhea: {result['error']}")
    api_direction = ('forward' if result.get('direction', 'forward') == 'forward'
                     else 'reversible')
    return process_reaction(
        substrates=_side(result['substrates'], result['subs_types'],
                         result['subs_sch'], result['subs_comps']),
        products=_side(result['products'], result['prods_types'],
                       result['prod_sch'], result['prods_comps']),
        direction=api_direction,
        subsystem=result.get('subsystem') or None,
        user=user,
    )

import re

import requests

RHEA_URL = 'https://www.rhea-db.org/rhea/'
TIMEOUT = 30
# Rhea gives each reaction four consecutive identifiers: the master (undefined
# direction) and its left-to-right, right-to-left and bidirectional variants.
_GROUP_SIZE = 4
_COMPONENT_SEPARATOR = re.compile(r'\s+\+\s+')  # names contain '+', e.g. 'NAD(+)'
_ARROW = re.compile(r'\s+(<=>|=>|<=|=)\s+')
_STOICHIOMETRY = re.compile(r'^(\d+)\s+(.*)$')


def _query_rhea(query):
    """Rows (dicts keyed by column title) of a Rhea TSV search."""
    response = requests.get(RHEA_URL, params={
        "query": query,
        "columns": "rhea-id,equation,chebi-id",
        "format": 'tsv',
        "limit": 10,
    }, timeout=TIMEOUT)
    if response.status_code != 200:
        raise requests.HTTPError(
            f"Failed to fetch data from Rhea. Status code: {response.status_code}",
            response=response)
    lines = [line for line in response.text.splitlines() if line.strip()]
    if not lines:
        return []
    headers = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        values = line.split("\t")
        if len(values) == len(headers):
            rows.append(dict(zip(headers, values)))
    return rows


def _rhea_number(rhea_id):
    match = re.match(r'^RHEA:(\d+)$', (rhea_id or '').strip(), re.IGNORECASE)
    return int(match.group(1)) if match else None


def _find_reaction(rhea_id):
    """The Rhea row for an identifier. Directional and bidirectional
    identifiers (e.g. RHEA:10031) resolve to their master reaction."""
    rows = _query_rhea(rhea_id)
    if rows:
        return rows[0]
    number = _rhea_number(rhea_id)
    if number is None:
        return None
    # "RHEA:10031" finds nothing, but the bare number finds the master
    # reaction of its group; only accept a master from the same group.
    for row in _query_rhea(str(number)):
        master = _rhea_number(row.get("Reaction identifier"))
        if master is not None and master <= number < master + _GROUP_SIZE:
            return row
    return None


def _parse_side(side):
    components = []
    for component in _COMPONENT_SEPARATOR.split(side.strip()):
        component = component.strip()
        if not component:
            continue
        match = _STOICHIOMETRY.match(component)
        if match:
            components.append((int(match.group(1)), match.group(2)))
        else:
            components.append((1, component))
    return components


def get_from_rhea(reaction_abbreviation):
    reaction_abbreviation = (reaction_abbreviation or '').strip()
    if 'RHEA:' != reaction_abbreviation[:5].upper():
        return {"error": "RHEA ID should start with 'RHEA:'"}
    reaction_abbreviation = 'RHEA:' + reaction_abbreviation[5:].strip()

    reaction_info = _find_reaction(reaction_abbreviation)
    if not reaction_info:
        return {"error": "No reaction found with this abbreviation"}

    equation = reaction_info.get("Equation", "")
    chebi_ids = [chebi.strip() for chebi in reaction_info.get("ChEBI identifier", "").split(';')
                 if chebi.strip()]
    if not equation or not chebi_ids:
        return {"error": "Equation or ChEBI IDs not found"}

    arrow = _ARROW.search(equation)
    if not arrow:
        return {"error": "Unknown equation format"}
    substrates = _parse_side(equation[:arrow.start()])
    products = _parse_side(equation[arrow.end():])
    if len(substrates) + len(products) != len(chebi_ids):
        return {"error": f"Could not match the Rhea equation to its ChEBI identifiers: {equation}"}
    if arrow.group(1) == '<=>':
        direction = "reversible"
    elif arrow.group(1) == '<=':
        direction = "reverse"
    else:
        direction = "forward"  # '=' (undefined direction) and '=>'

    substrate_chebis = chebi_ids[:len(substrates)]
    product_chebis = chebi_ids[len(substrates):]
    return {
        "rhea_id": reaction_info.get("Reaction identifier", reaction_abbreviation),
        "equation": equation,
        "direction": direction,
        "subs_sch": [stoichiometry for stoichiometry, _ in substrates],
        "prod_sch": [stoichiometry for stoichiometry, _ in products],
        "subs_types": ["ChEBI ID"] * len(substrate_chebis),
        "prods_types": ["ChEBI ID"] * len(product_chebis),
        "substrates": substrate_chebis,
        "products": product_chebis,
        "subs_comps": ["-"] * len(substrate_chebis),
        "prods_comps": ["-"] * len(product_chebis),
        "subsystem": ""
    }

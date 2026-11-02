"""Constructor MCP server (feature 002-curation-api-mcp).

A thin stdio MCP wrapper over the Constructor API v1. It talks ONLY HTTP to
/api/v1/ — it never imports the Django app or touches the database — and it
exposes no tool that writes to VMH.

Configuration (environment):
    CONSTRUCTOR_BASE_URL  Base URL of a running Constructor instance
                          (default: http://localhost:8000)
    CONSTRUCTOR_API_KEY   Your personal API key (generate it at
                          <base-url>/api-key/ while signed in)

Run:  python server.py            # stdio transport, for MCP clients
      python server.py --smoke    # list tools and exit (no key needed)
"""
import json
import os
import sys

import httpx
from mcp.server.mcpserver import MCPServer

BASE_URL = os.environ.get('CONSTRUCTOR_BASE_URL', 'http://localhost:8000').rstrip('/')
API_KEY = os.environ.get('CONSTRUCTOR_API_KEY', '')

ID_TYPES = 'vmh, chebi, chebi_name, pubchem, swisslipids, smiles, saved'
COMPARTMENTS = "c, e, g, m, l, n, r, x, '-'"

mcp = MCPServer(
    'constructor',
    instructions=(
        'Curation tools for Constructor (metabolic reaction curation for '
        'genome-scale reconstructions, VMH ecosystem). Reactions are created '
        "and balance-checked in the user's private workspace; nothing is ever "
        'submitted to the shared VMH database.'),
)


def _call(method, path, payload=None):
    """Call the Constructor API and return a JSON string, or an
    'Error [code]: message' string the assistant can relay."""
    if not API_KEY:
        return ('Error [config]: CONSTRUCTOR_API_KEY is not set. Generate a '
                f'key at {BASE_URL}/api-key/ and set the environment variable.')
    try:
        response = httpx.request(
            method, f'{BASE_URL}/api/v1{path}',
            json=payload,
            headers={'Authorization': f'Bearer {API_KEY}'},
            timeout=120,
        )
    except httpx.HTTPError as exc:
        return f'Error [connection]: could not reach {BASE_URL}: {exc}'
    try:
        body = response.json()
    except ValueError:
        return f'Error [protocol]: non-JSON response (HTTP {response.status_code})'
    if response.status_code >= 400:
        error = body.get('error', {})
        return (f"Error [{error.get('code', response.status_code)}]: "
                f"{error.get('message', 'request failed')}")
    return json.dumps(body, indent=2)


@mcp.tool()
def resolve_metabolite(identifier: str, id_type: str) -> str:
    """Resolve a metabolite identifier to name, molecular formula, charge,
    structure (SMILES/InChI) and whether it already exists in VMH.

    id_type must be one of: vmh, chebi, chebi_name, pubchem, swisslipids,
    smiles, saved. Examples: ('h2o', 'vmh'), ('CHEBI:15377', 'chebi'),
    ('O', 'smiles')."""
    return _call('POST', '/metabolites/resolve',
                 {'identifier': identifier, 'type': id_type})


@mcp.tool()
def create_reaction(substrates: list, products: list,
                    direction: str = 'forward',
                    subsystem: str | None = None) -> str:
    """Create and balance-check a reaction. Returns mass-balance and
    charge-balance verdicts, per-side atom/charge tallies, resolved metabolite
    info and whether an identical/similar reaction exists in VMH. The reaction
    is stored in the user's private workspace (NOT submitted to VMH); use
    save_reaction to add it to their saved collection.

    Each substrate/product is a dict:
      {"identifier": "h2o", "type": "vmh", "stoichiometry": 1,
       "compartment": "c", "name": null}
    type: one of vmh, chebi, chebi_name, pubchem, swisslipids, smiles, saved
    (name is required for type 'smiles'). compartment: one of c, e, g, m, l,
    n, r, x, '-'. direction: 'forward' or 'reversible'."""
    return _call('POST', '/reactions', {
        'substrates': substrates,
        'products': products,
        'direction': direction,
        'subsystem': subsystem,
    })


@mcp.tool()
def save_reaction(reaction_id: int, name: str, description: str = '') -> str:
    """Save a previously created/fetched reaction to the user's collection
    under a name (names are unique per user)."""
    return _call('POST', f'/reactions/{reaction_id}/save',
                 {'name': name, 'description': description})


@mcp.tool()
def list_reactions(page: int = 1, page_size: int = 25) -> str:
    """List the user's saved reactions (id, name, formula, balance status)."""
    return _call('GET', f'/reactions?page={page}&page_size={page_size}')


@mcp.tool()
def get_reaction(reaction_id: int) -> str:
    """Get full details (balance verdicts, metabolites, VMH match) of one of
    the user's reactions."""
    return _call('GET', f'/reactions/{reaction_id}')


@mcp.tool()
def fetch_vmh_reaction(abbreviation: str) -> str:
    """Import a reaction from VMH by its abbreviation (e.g. 'HEX1') and run it
    through processing: returns it balance-checked, ready to review or save.
    Read-only with respect to VMH."""
    return _call('POST', '/sources/vmh/fetch', {'abbreviation': abbreviation})


@mcp.tool()
def fetch_rhea_reaction(rhea_id: str) -> str:
    """Import a reaction from Rhea by id (e.g. '10031' or 'RHEA:10031') and
    run it through processing: returns it balance-checked, ready to review or
    save."""
    return _call('POST', '/sources/rhea/fetch', {'rhea_id': rhea_id})


if __name__ == '__main__':
    if '--smoke' in sys.argv:
        import asyncio
        tools = asyncio.run(mcp.list_tools())
        print(f'constructor MCP server: {len(tools)} tools')
        for tool in tools:
            print(f'  - {tool.name}')
        sys.exit(0)
    mcp.run()

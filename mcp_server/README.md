# Constructor MCP Server

A thin [MCP](https://modelcontextprotocol.io) (stdio) wrapper over the Constructor
API v1, so AI assistants (Claude Desktop, Claude Code, any MCP client) can resolve
metabolites, build and balance-check reactions, import from VMH/Rhea, and manage your
saved reactions — conversationally.

It talks **only HTTP** to `/api/v1/` of a Constructor instance and exposes **no tool
that writes to VMH**.

## Setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Get an API key: sign in to Constructor, open `<base-url>/api-key/`, generate your key
(it is shown only once).

Quick check:

```bash
CONSTRUCTOR_API_KEY=ctor_... python server.py --smoke   # lists the 7 tools
```

## Register with an MCP client

```json
{
  "mcpServers": {
    "constructor": {
      "command": "/path/to/reconstructor/mcp_server/.venv/bin/python",
      "args": ["/path/to/reconstructor/mcp_server/server.py"],
      "env": {
        "CONSTRUCTOR_BASE_URL": "https://constructor.humanmetabolism.org",
        "CONSTRUCTOR_API_KEY": "ctor_..."
      }
    }
  }
}
```

For a local instance use `"CONSTRUCTOR_BASE_URL": "http://localhost:8000"`.

## Tools

| Tool | What it does |
|------|--------------|
| `resolve_metabolite` | Identifier (VMH/ChEBI/PubChem/SwissLipids/SMILES/saved) → name, formula, charge, structure, VMH presence |
| `create_reaction` | Build a reaction and get mass/charge balance verdicts + VMH match |
| `save_reaction` | Save a created/fetched reaction to your collection |
| `list_reactions` | List your saved reactions |
| `get_reaction` | Full details of one of your reactions |
| `fetch_vmh_reaction` | Import a VMH reaction by abbreviation, balance-checked |
| `fetch_rhea_reaction` | Import a Rhea reaction by id, balance-checked |

Example prompt: *“Fetch VMH reaction HEX1 and tell me whether it is mass- and
charge-balanced.”*

## Errors

API errors surface as `Error [<code>]: <message>` (codes: `unauthorized`,
`rate_limited`, `validation_error`, `unknown_identifier`, `upstream_unavailable`,
`not_found`, `duplicate_name`, …) so the assistant can react; connection problems
appear as `Error [connection]: …`.

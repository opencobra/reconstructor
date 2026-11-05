# Constructor MCP Server

A thin [MCP](https://modelcontextprotocol.io) (stdio) wrapper over the Constructor
API v1, so AI assistants (Claude Desktop, Claude Code, any MCP client) can resolve
metabolites, build and balance-check reactions, import from VMH/Rhea, and manage your
saved reactions — conversationally.

It talks **only HTTP** to `/api/v1/` of a Constructor instance and exposes **no tool
that writes to VMH**.

## Get an API key

Sign in to Constructor, open `<base-url>/api-key/`, and generate your key
(`ctor_...` — it is shown only once).

## Install

Pick whichever suits you; all three expose the same 7 tools.

### Option A — `uvx` one-liner (recommended; needs [uv](https://docs.astral.sh/uv/))

No clone, no venv. Add to your MCP client config (Claude Desktop:
`claude_desktop_config.json`; Claude Code: `.mcp.json` or `claude mcp add`):

```json
{
  "mcpServers": {
    "constructor": {
      "command": "uvx",
      "args": ["constructor-mcp"],
      "env": {
        "CONSTRUCTOR_BASE_URL": "https://constructor.humanmetabolism.org",
        "CONSTRUCTOR_API_KEY": "ctor_..."
      }
    }
  }
}
```

Until the package is on PyPI, point `uvx` at this repository instead:

```json
      "args": ["--from",
               "git+https://github.com/opencobra/reconstructor#subdirectory=mcp_server",
               "constructor-mcp"],
```

### Option B — one-click Claude Desktop bundle (`.mcpb`)

Download `constructor-mcp.mcpb` (from a release, or build it yourself with
`./build_mcpb.sh`), open it with Claude Desktop, and fill in the two settings
it asks for — Constructor URL and your API key (stored in the OS keychain).
No terminal, no config files. Requires a system `python3` (≥ 3.10).

### Option C — from source

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
CONSTRUCTOR_API_KEY=ctor_... python server.py --smoke   # lists the 7 tools
```

then register it with absolute paths:

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

Note: VMH-backed tools (`fetch_vmh_reaction`, resolving `vmh` identifiers)
additionally require the **Constructor server** to have `VMH_API_KEY` set in its
environment. If it does not, those tools return an `upstream_unavailable` error
that says so and explains the fix; SMILES/ChEBI/PubChem resolution and Rhea
imports work without it.

# Contract: Constructor MCP Server

Transport: stdio (official `mcp` Python SDK, FastMCP). Location: `mcp_server/server.py`.

**Config (env)**: `CONSTRUCTOR_BASE_URL` (default `http://localhost:8000`),
`CONSTRUCTOR_API_KEY` (required). The server is a pure HTTP client of `/api/v1/` —
it imports nothing from the Django app (FR-009).

**Error passthrough**: an API error becomes a tool result of the form
`Error [<code>]: <message>` so the assistant can relay/act on it (US4 scenario 2);
transport/connection failures are reported as `Error [connection]: ...`.

## Tools (7) — 1:1 with API v1

| Tool | Wraps | Arguments |
|------|-------|-----------|
| `resolve_metabolite` | POST /metabolites/resolve | `identifier: str`, `id_type: str` |
| `create_reaction` | POST /reactions | `substrates: list[dict]`, `products: list[dict]`, `direction: str = "forward"`, `subsystem: str \| None` |
| `save_reaction` | POST /reactions/{id}/save | `reaction_id: int`, `name: str`, `description: str = ""` |
| `list_reactions` | GET /reactions | `page: int = 1`, `page_size: int = 25` |
| `get_reaction` | GET /reactions/{id} | `reaction_id: int` |
| `fetch_vmh_reaction` | POST /sources/vmh/fetch | `abbreviation: str` |
| `fetch_rhea_reaction` | POST /sources/rhea/fetch | `rhea_id: str` |

Each tool's docstring states: supported identifier types, compartment codes, that
`create_reaction` returns mass/charge balance verdicts, and that nothing is ever
submitted to VMH (US4 scenario 3: no VMH-write tool exists).

## Client registration example (README)

```json
{
  "mcpServers": {
    "constructor": {
      "command": "python",
      "args": ["/path/to/reconstructor/mcp_server/server.py"],
      "env": {
        "CONSTRUCTOR_BASE_URL": "https://constructor.humanmetabolism.org",
        "CONSTRUCTOR_API_KEY": "ctor_..."
      }
    }
  }
}
```

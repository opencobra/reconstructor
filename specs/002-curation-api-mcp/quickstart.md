# Quickstart: Curation API & MCP Server

Validation guide for feature 002. Contracts: [api-v1.md](./contracts/api-v1.md),
[mcp-tools.md](./contracts/mcp-tools.md).

## Prerequisites

```bash
cd /path/to/reconstructor
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt          # now includes django-ninja
```

## Run the automated tests (SQLite, upstreams mocked)

```bash
cd curationTool
python manage.py test reactions --settings=reactions_project.test_settings
```

Expected: all tests pass; suites cover api-key lifecycle, metabolite resolve, reaction
create/balance/save/list/get, owner isolation, error envelope, VMH/Rhea fetch.

## Smoke-run the server locally (no Postgres needed)

```bash
cd curationTool
python manage.py migrate --settings=reactions_project.test_settings
python manage.py runserver --settings=reactions_project.test_settings
```

1. Register a user in the web UI (`http://localhost:8000`), then open
   `http://localhost:8000/api-key/` and generate a key (shown once).
2. Health + auth:

   ```bash
   curl http://localhost:8000/api/v1/health
   curl -H "Authorization: Bearer $KEY" http://localhost:8000/api/v1/me
   ```

3. Balance check a known-balanced reaction offline (SMILES, no upstream calls):

   ```bash
   curl -X POST http://localhost:8000/api/v1/reactions \
     -H "Authorization: Bearer $KEY" -H "Content-Type: application/json" \
     -d '{"substrates":[{"identifier":"O","type":"smiles","stoichiometry":1,"compartment":"c","name":"water"}],
          "products":[{"identifier":"O","type":"smiles","stoichiometry":1,"compartment":"e","name":"water"}],
          "direction":"forward"}'
   ```

   Expect `201` with `balanced_count`/`balanced_charge` true (transport reaction).
4. Save, list, fetch:

   ```bash
   curl -X POST .../api/v1/reactions/<id>/save -d '{"name":"water transport"}' ...
   curl -H "Authorization: Bearer $KEY" .../api/v1/reactions
   ```

5. OpenAPI docs: open `http://localhost:8000/api/v1/docs`.

## MCP server

```bash
python3 -m venv mcp_server/.venv && source mcp_server/.venv/bin/activate
pip install -r mcp_server/requirements.txt
CONSTRUCTOR_BASE_URL=http://localhost:8000 CONSTRUCTOR_API_KEY=$KEY \
  python mcp_server/server.py --smoke   # lists tools + calls resolve_metabolite once
```

Register in an MCP client with the snippet in `mcp_server/README.md`, then ask the
assistant e.g. *"Fetch VMH reaction HEX1 and tell me if it is balanced"* (SC-004).

## Acceptance sweep (maps to spec)

| Check | Spec |
|-------|------|
| Balanced + unbalanced SMILES reactions return correct verdicts | US1, SC-002 |
| Unknown identifier → 422 naming side+index | US1-2, FR-008 |
| Saved reaction visible in web UI saved list | US1-3 |
| Foreign reaction id → 404, no leak | US3-2, FR-006 |
| Old key fails after rotation | US5, FR-003 |
| No VMH-write endpoint/tool anywhere | FR-010, SC-005 |
| Existing web test-free flows unchanged (manual spot check: create+save in UI) | FR-011, SC-006 |

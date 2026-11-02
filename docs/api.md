# Constructor API v1

Constructor exposes its core curation capabilities as a versioned JSON API under
`/api/v1/`, so scripts and AI agents can resolve metabolites, build and
balance-check reactions, and manage saved reactions without the web UI.

- Interactive docs (Swagger UI): `<base-url>/api/v1/docs`
- OpenAPI schema: `<base-url>/api/v1/openapi.json`
- Full design artifacts: [`specs/002-curation-api-mcp/`](../specs/002-curation-api-mcp/)
- MCP server (for AI assistants): [`mcp_server/README.md`](../mcp_server/README.md)

Submitting reactions or metabolites **to VMH is deliberately not possible** through
the API or the MCP server. That remains a credentialed, web-UI-only action.

## Authentication

Every endpoint except `GET /health` requires a personal API key:

1. Sign in to Constructor and open `<base-url>/api-key/`.
2. Click **Generate** — the full key (`ctor_...`) is shown **only once**; store it
   somewhere safe. Only a SHA-256 hash is kept on the server.
3. Generating again **rotates** the key: the old key stops working immediately.

Send the key either way:

```bash
curl -H "Authorization: Bearer ctor_..." <base-url>/api/v1/me
curl -H "X-API-Key: ctor_..."            <base-url>/api/v1/me
```

## Rate limiting

Requests are throttled per key (default **60/minute**, configurable with the
`API_THROTTLE_RATE` environment variable, e.g. `100/m`). Exceeding it returns
`429` with code `rate_limited`.

## Errors

All non-2xx responses share one envelope:

```json
{"error": {"code": "unknown_identifier", "message": "...", "detail": null}}
```

| HTTP | Code | Meaning |
|------|------|---------|
| 401 | `unauthorized` | Missing, malformed, or rotated key |
| 404 | `not_found` | Reaction does not exist **or is not yours** |
| 409 | `duplicate_name` | Reaction/metabolite name already used by you |
| 422 | `validation_error` | Malformed body (message names the field) |
| 422 | `unsupported_type` | Unknown metabolite `type` (message lists supported ones) |
| 422 | `unknown_identifier` | Identifier could not be resolved (for reactions, the message names the side and index, e.g. `substrates[0]`) |
| 429 | `rate_limited` | Throttle exceeded |
| 500 | `processing_error` | Unexpected failure in the curation pipeline |
| 502 | `upstream_unavailable` | VMH / ChEBI / PubChem / SwissLipids / Rhea unreachable |

## Metabolite inputs

Reactions and resolution take metabolites in this shape:

```json
{"identifier": "h2o", "type": "vmh", "stoichiometry": 1, "compartment": "c", "name": null}
```

- `type`: `vmh`, `chebi`, `chebi_name`, `pubchem`, `swisslipids`, `smiles`, `saved`
  (file upload and structure drawing are web-UI-only). `saved` identifiers must
  reference metabolites you own or that were shared with you.
- `stoichiometry`: whole number between 1 and 100 (default 1).
- `compartment`: one of `c`, `e`, `g`, `m`, `l`, `n`, `r`, `x`, `-` (default `c`).
- `name`: **required for `type: "smiles"`**; optional otherwise.

Metabolites that are not found in VMH are automatically stored in your personal
saved-metabolite collection (reusable later with `type: "saved"` and the metabolite
name as identifier). A name clash with a *different* structure returns
409 `duplicate_name`.

## Endpoints

### `GET /api/v1/health` — liveness (no auth)

```json
{"status": "ok", "version": "v1"}
```

### `GET /api/v1/me` — identify the caller

```json
{"user_id": 1, "name": "farid", "full_name": "Farid Zare",
 "key_prefix": "ctor_AbC1", "key_created_at": "2026-11-02T12:00:00+00:00"}
```

### `POST /api/v1/metabolites/resolve`

```bash
curl -X POST <base-url>/api/v1/metabolites/resolve \
  -H "Authorization: Bearer ctor_..." -H "Content-Type: application/json" \
  -d '{"identifier": "CHEBI:15377", "type": "chebi"}'
```

Returns name, formula, charge, SMILES, InChI, InChI key, molecular weight, and
whether VMH already has it (`vmh.found`, `vmh.abbreviation`, `vmh.miriam`).

### `POST /api/v1/reactions` — create + balance-check

```bash
curl -X POST <base-url>/api/v1/reactions \
  -H "Authorization: Bearer ctor_..." -H "Content-Type: application/json" \
  -d '{
    "substrates": [{"identifier": "h2o", "type": "vmh", "compartment": "c"}],
    "products":   [{"identifier": "h2o", "type": "vmh", "compartment": "e"}],
    "direction": "forward"
  }'
```

- `direction`: `forward` or `reversible`.
- `subsystem` (string) and `organs` (list of strings) are optional.
- At least one substrate or product is required.

Returns `201` with a `ProcessedReaction`: `reaction_id`, mass verdict
(`balanced_count`), charge verdict (`balanced_charge`), per-side atom counts
(`subs_atoms`, `prods_atoms`) and charges, per-metabolite names/formulas/charges/
SMILES/InChI keys/molecular weights, VMH presence per metabolite (`subs_found`,
`prod_found`, miriams), and whether the whole reaction already exists in VMH
(`vmh.found` / `vmh.similar` / `vmh.url` / `vmh.formula`).

The reaction is persisted (it appears in the web UI) but is **not yet in your
saved list** — mirroring the web flow, where creation and saving are separate steps.

### `POST /api/v1/reactions/{id}/save`

```bash
curl -X POST <base-url>/api/v1/reactions/42/save \
  -H "Authorization: Bearer ctor_..." -H "Content-Type: application/json" \
  -d '{"name": "water transport", "description": ""}'
```

Names the reaction and adds it to your saved collection. `409 duplicate_name` if
you already have a reaction with that name; `404` if the reaction is not yours.

### `GET /api/v1/reactions?page=1&page_size=25`

Lists **your saved reactions only**:

```json
{"count": 1, "page": 1, "page_size": 25,
 "results": [{"reaction_id": 42, "name": "water transport",
              "rxn_formula": "...", "balanced_mass": true,
              "balanced_charge": true, "vmh_found": false}]}
```

### `GET /api/v1/reactions/{id}`

Full `ProcessedReaction` view of a reaction you created or saved. Anything else —
including other users' reactions — is `404 not_found` (existence is not leaked).

### `POST /api/v1/sources/vmh/fetch`

```bash
curl -X POST <base-url>/api/v1/sources/vmh/fetch \
  -H "Authorization: Bearer ctor_..." -H "Content-Type: application/json" \
  -d '{"abbreviation": "HEX1"}'
```

Imports the VMH reaction, runs it through the same balance-check pipeline, and
returns `201 ProcessedReaction`. Unknown abbreviation → `422 unknown_identifier`;
VMH down → `502 upstream_unavailable`. Reactions with fractional coefficients
(e.g. `0.5 o2[c]`) are rejected with `422 validation_error` rather than silently
rounded.

### `POST /api/v1/sources/rhea/fetch`

Same as above with `{"rhea_id": "10031"}` (the `RHEA:` prefix is optional).

## Worked example

```bash
BASE=http://localhost:8000
KEY=ctor_...

# 1. Is water already in VMH?
curl -s -X POST $BASE/api/v1/metabolites/resolve \
  -H "Authorization: Bearer $KEY" -H "Content-Type: application/json" \
  -d '{"identifier": "O", "type": "smiles"}'

# 2. Create 2 H2O -> H2O2 + H2 and check balance
RXN=$(curl -s -X POST $BASE/api/v1/reactions \
  -H "Authorization: Bearer $KEY" -H "Content-Type: application/json" \
  -d '{"substrates": [{"identifier": "O", "type": "smiles", "name": "water", "stoichiometry": 2}],
       "products": [{"identifier": "OO", "type": "smiles", "name": "hydrogen peroxide"},
                    {"identifier": "[H][H]", "type": "smiles", "name": "dihydrogen"}],
       "direction": "forward"}')
ID=$(echo "$RXN" | python3 -c 'import json,sys; print(json.load(sys.stdin)["reaction_id"])')

# 3. Save it
curl -s -X POST $BASE/api/v1/reactions/$ID/save \
  -H "Authorization: Bearer $KEY" -H "Content-Type: application/json" \
  -d '{"name": "water disproportionation"}'

# 4. List your reactions
curl -s $BASE/api/v1/reactions -H "Authorization: Bearer $KEY"
```

## Using it from an AI assistant (MCP)

The repository ships a Model Context Protocol server
([`mcp_server/`](../mcp_server/)) that wraps this API as 7 conversational tools
for Claude Desktop, Claude Code, or any MCP client. Setup, client configuration,
and the tool list are in [`mcp_server/README.md`](../mcp_server/README.md).

# Implementation Plan: Curation API & MCP Server

**Branch**: `002-curation-api-mcp` | **Date**: 2026-11-02 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `specs/002-curation-api-mcp/spec.md`

## Summary

Expose Constructor's core curation capabilities (metabolite resolution, reaction
creation with mass/charge balance verdicts, saved-reaction management, VMH/Rhea import)
as a versioned, API-key-authenticated JSON API mounted at `/api/v1/` inside the existing
Django monolith, implemented with Django Ninja over a new `reactions/services/`
orchestration layer that reuses the existing utils chain (the domain logic the web UI
already runs). A thin stdio MCP server (`mcp_server/`, official MCP Python SDK) wraps
the API 1:1 so AI assistants can curate conversationally. Push-to-VMH is excluded from
both surfaces. Constitution amended (v1.1.0, user-approved 2026-11-02) to sanction the
API tier.

## Technical Context

**Language/Version**: Python 3.12 (server), Django 5.0.1

**Primary Dependencies**: django-ninja ~=1.4 (new), pydantic 2.11 (existing), RDKit
2023.9 (existing), requests/httpx (existing); MCP server: `mcp` SDK ≥1.2 + httpx

**Storage**: Existing PostgreSQL via Django ORM; one new table (`ApiKey`); SQLite via
`test_settings.py` for tests/local smoke runs

**Testing**: Django `TestCase` under `reactions/tests/` (new package), external sources
mocked; MCP stdio smoke test; manual acceptance per spec scenarios

**Target Platform**: Linux server (existing Docker/gunicorn deployment unchanged); MCP
server runs wherever the curator's MCP client runs (stdio)

**Project Type**: Web application (Django monolith) + small sidecar CLI package (MCP)

**Performance Goals**: Reaction processing ≤10 s typical (≤6 metabolites, no upstream
fetch); metabolite resolve ≤5 s for VMH/SMILES (SC-003); atom mapping (RDT) always
skipped on the API path

**Constraints**: Zero behaviour change to existing pages/endpoints (FR-011, SC-006);
no VMH writes from API/MCP (FR-010); additive migration only; one key per user;
rate limit 60/min/key default (FR-012)

**Scale/Scope**: ~8 API endpoints + 1 model + 1 migration + services package (~4
modules) + key-management page + 7 MCP tools; research-tool traffic (tens of users)

## Constitution Check

*GATE: evaluated against constitution v1.1.0 (amended this feature, see below).*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Server-Rendered Django Core | PASS (amended) | Web UI untouched; API tier is the sanctioned exception added in v1.1.0 (user approved 2026-11-02). API renders no pages; key-management page is a normal Django view + template. |
| II. Modular Vanilla JS | PASS | Only JS change: a small `apiKey.js` module for the key page, file-per-concern, Semantic UI idioms, no build step. |
| III. UX Consistency | PASS | Key page uses existing Semantic UI card/button/modal patterns and header layout. |
| IV. Data & Integration Integrity | PASS | Reaction schema unchanged; `ApiKey` is additive; service persists reactions with identical field encoding to the web flow; balance semantics computed by the same shared functions; VMH/Rhea contracts untouched; no VMH writes. |
| V. Spec Traceability | PASS | Every task in tasks.md maps to FR-001…FR-013 / SC-001…SC-006. |
| VI. Minimal Footprint | PASS | One new dependency in the app (django-ninja); MCP deps isolated in `mcp_server/`; shared utils extended additively (one `SMILES` branch) instead of parallel converters. |

**Constitution amendment (performed with this plan)**: v1.0.0 → v1.1.0 (MINOR).
Principle I gains: *"Sanctioned exception (2026-11-02, feature 002): a versioned,
token-authenticated programmatic API under `/api/v1/` (Django Ninja) with a thin MCP
wrapper may expose curation capabilities, provided the web UI remains server-rendered,
domain logic is shared via `reactions/services/` and existing utils, and no programmatic
surface writes to VMH."* Technology Constraints updated accordingly.

## Project Structure

### Documentation (this feature)

```text
specs/002-curation-api-mcp/
├── plan.md              # This file
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
├── quickstart.md        # Phase 1 output
├── contracts/
│   ├── api-v1.md        # REST contract (endpoints, schemas, errors)
│   └── mcp-tools.md     # MCP tool surface contract
├── checklists/requirements.md
└── tasks.md             # Phase 2 output (/speckit-tasks)
```

### Source Code (repository root)

```text
curationTool/
├── reactions/
│   ├── models.py                    # + ApiKey model (additive)
│   ├── migrations/000X_apikey.py    # new, additive
│   ├── services/                    # NEW shared orchestration layer
│   │   ├── __init__.py
│   │   ├── errors.py                # typed service errors → API error codes
│   │   ├── api_keys.py              # generate/rotate/authenticate
│   │   ├── metabolite_resolution.py # resolve_metabolite(identifier, type)
│   │   ├── reaction_processing.py   # process_reaction(submission, user)
│   │   └── external_sources.py      # fetch_vmh_reaction / fetch_rhea_reaction
│   ├── api/                         # NEW Django Ninja v1 API
│   │   ├── __init__.py              # NinjaAPI instance, auth, throttle, error handlers
│   │   ├── auth.py                  # ApiKeyAuth (Bearer / X-API-Key)
│   │   ├── schemas.py               # pydantic request/response models
│   │   └── routes/
│   │       ├── metabolites.py
│   │       ├── reactions.py
│   │       └── sources.py           # VMH / Rhea fetch
│   ├── views/user_views.py          # + api_key_page / api_key_generate (session auth)
│   ├── templates/reactions/api_key.html   # NEW minimal Semantic UI page
│   ├── static/reactions/js/apiKey.js      # NEW page module
│   ├── utils/to_mol.py              # + 'SMILES' branch (additive)
│   ├── tests/                       # NEW package (replaces empty tests.py)
│   │   ├── __init__.py
│   │   ├── test_api_keys.py
│   │   ├── test_metabolites_api.py
│   │   ├── test_reactions_api.py
│   │   └── test_sources_api.py
│   └── urls.py                      # + api key page routes
├── reactions_project/
│   ├── urls.py                      # + mount api at /api/v1/
│   ├── settings.py                  # + API_THROTTLE_RATE (env), nothing else
│   └── test_settings.py             # NEW: SQLite override for tests/smoke
mcp_server/                          # NEW sidecar (not a Django app)
├── server.py                        # FastMCP stdio server, 7 tools
├── requirements.txt                 # mcp, httpx
└── README.md                        # setup + client config snippet
requirements.txt                     # + django-ninja
```

**Structure Decision**: Django monolith extended in place (`reactions/api/`,
`reactions/services/`) so the API shares the ORM and utils; the MCP server is a sidecar
directory outside `curationTool/`, which structurally prevents ORM access (FR-009) and
keeps its dependencies out of the web app image.

## Design Decisions (from research.md)

- **R1** Django Ninja at `/api/v1/`, OpenAPI at `/api/v1/openapi.json`.
- **R2** `ApiKey`: one-to-one with custom `User`; SHA-256 digest stored, shown once,
  rotate = overwrite; `Authorization: Bearer` or `X-API-Key`.
- **R3** Service layer orchestrates existing utils; `input_reaction` view untouched
  (see Complexity Tracking); persisted rows byte-compatible with web-created ones.
- **R4** Additive `SMILES` type in `any_to_mol`; API accepts
  `vmh|chebi|chebi_name|pubchem|swisslipids|smiles|saved`; file/draw rejected.
- **R5** Ninja `AuthRateThrottle`, default 60/min/key, env-tunable; 429 + error envelope.
- **R6** `mcp_server/server.py` (FastMCP, stdio), env-configured, HTTP-only to the API.
- **R7** SQLite `test_settings.py`; mocked upstreams; offline SMILES balance tests;
  MCP stdio smoke test.
- **R8** Constitution v1.1.0 amendment (Principle I exception + tech constraints).

**Error envelope (FR-013)**: every non-2xx response is
`{"error": {"code": "<machine_code>", "message": "<human text>", "detail": {...}}}`
with codes: `unauthorized`, `rate_limited`, `validation_error`, `unknown_identifier`,
`unsupported_type`, `upstream_unavailable`, `not_found`, `duplicate_name`,
`processing_error`.

## Complexity Tracking

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| Orchestration logic exists in both `input_reaction` (web) and `services/reaction_processing.py` (API) rather than one callsite | `input_reaction` is the app's untested critical path and carries web-only concerns (uploads, drawings, RDT visualisation, template error rendering); rewiring it now directly contradicts FR-011/SC-006 (zero web regressions) | "Refactor the view through the service in the same feature" — rejected: no regression harness exists for the web flow; the actual balance/chemistry logic is already a single implementation (shared utils), so verdict divergence is impossible. Follow-up noted for a dedicated refactor feature once API tests establish the harness. |
| New top-level `mcp_server/` directory (outside the Django tree) | FR-009 requires MCP to use only the public API; placing it outside `curationTool/` enforces that structurally and keeps `mcp`/`httpx` out of the web image | "Put it in `curationTool/` as an app" — rejected: it is not a Django app, would tempt ORM imports, and `mcp/` would shadow the SDK package name. |

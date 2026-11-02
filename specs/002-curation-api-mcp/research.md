# Research: Curation API & MCP Server

**Feature**: 002-curation-api-mcp | **Date**: 2026-11-02

## R1. API framework

- **Decision**: Django Ninja (`django-ninja ~= 1.4`), mounted at `/api/v1/` from the
  existing `reactions_project/urls.py`.
- **Rationale**: Runs inside the existing Django monolith (no second service, satisfies
  Minimal Footprint); request/response schemas are pydantic models — `pydantic 2.11` is
  already in `requirements.txt`; OpenAPI document generated automatically at
  `/api/v1/openapi.json` (needed for docs and for keeping the MCP tool surface honest);
  built-in throttling module covers FR-012 without a new dependency.
- **Alternatives considered**: Django REST Framework — heavier (serializers, routers,
  renderer stack), no pydantic, OpenAPI needs extra packages; plain Django JSON views —
  no schema/validation/OpenAPI, would hand-roll what Ninja provides; FastAPI — a second
  ASGI service beside Django, violates the monolith constraint.

## R2. API key storage & authentication

- **Decision**: New `ApiKey` model, one-to-one with the existing custom
  `reactions.models.User`. The key is `ctor_` + `secrets.token_urlsafe(32)`, shown once;
  only its SHA-256 hex digest is stored (`key_hash`, unique indexed) plus a display
  `prefix` (first 12 chars). Authentication hashes the presented key and looks up by
  digest. Regeneration overwrites the hash (old key dies instantly). Header:
  `Authorization: Bearer <key>` (also accepts `X-API-Key`).
- **Rationale**: The project uses a custom User (not `django.contrib.auth`), so DRF/Ninja
  bundled token models don't apply. Hash-only storage means a DB leak does not leak keys;
  digest lookup is index-based and does not compare secrets byte-by-byte in Python, so
  timing attacks reveal nothing about the key material.
- **Alternatives considered**: Plaintext key column (simpler, but the existing User model
  already stores plaintext passwords — not a pattern to extend); signed JWTs (no server-
  side revocation without a blocklist; overkill for one key per user).

## R3. Service-layer seam (sharing curation logic with the web UI)

- **Decision**: The shared domain layer is the existing utils chain —
  `any_to_mol` → `get_mol_info` → `construct_reaction_rxnfile` → `get_reaction_info`
  → `check_reaction_vmh` → `search_metabolites_vmh` — plus a new orchestration package
  `curationTool/reactions/services/` (`reaction_processing.py`,
  `metabolite_resolution.py`, `external_sources.py`, `api_keys.py`) that the API calls.
  The service persists `Reaction` rows with **the same field encoding** the web view
  uses (JSON-dumped lists, types mutated to `Saved` for newly stored metabolites), so
  API-created reactions render correctly in the web UI. `input_reaction` (the web view)
  is **not refactored in this feature**: it is the app's critical path, has no test
  harness, and handles web-only concerns (file uploads, ChemDoodle drawings, RDT atom
  mapping, template error rendering). Balance verdicts cannot diverge because both
  surfaces call the same `get_reaction_info`/`get_mol_info` functions.
- **Rationale**: FR-011 requires shared curation logic and zero web regressions; the
  utils *are* the curation logic. Rewiring the 400-line view through the new service in
  the same change would risk exactly the regressions FR-011/SC-006 forbid. Recorded in
  plan Complexity Tracking as deliberate, with a follow-up note.
- **Alternatives considered**: Full refactor of `input_reaction` to call the service
  (rejected for risk, see above); calling the view functions from the API via request
  faking (rejected: couples API to HTML/form semantics and `request.FILES`).
- **Verified**: `request` is used inside `any_to_mol`/`search_metabolites_vmh`/
  `get_fields` only for `request.FILES` (file/draw types). All identifier types work
  with `request=None` — `check_reaction_vmh` already calls them that way internally.

## R4. SMILES input type

- **Decision**: Add a `'SMILES'` branch to `any_to_mol` (and a `smiles_to_mol` helper):
  `Chem.MolFromSmiles` + explicit-H normalisation identical to the other converters.
  API type values: `vmh`, `chebi`, `chebi_name`, `pubchem`, `swisslipids`, `smiles`,
  `saved` — mapped to the internal type strings. `MDL Mol file` and `Draw` are rejected
  by API validation (out of scope per spec Assumptions).
- **Rationale**: Spec FR-004 includes SMILES; the converter set has no raw-SMILES entry
  today. An additive branch in the shared util keeps one dispatch table (Principle VI)
  and gives tests a fully offline input type (no upstream fetch).

## R5. Rate limiting

- **Decision**: Django Ninja's `AuthRateThrottle` applied API-wide, default `60/min`
  per key, configurable via env (`API_THROTTLE_RATE`). Over-limit returns HTTP 429 with
  the standard error envelope (FR-012, FR-013).
- **Alternatives considered**: `django-ratelimit` (extra dependency), reverse-proxy
  limiting (not portable to local/dev runs).

## R6. MCP server

- **Decision**: Separate top-level directory `mcp_server/` containing a single-file
  FastMCP server (`server.py`, official `mcp` Python SDK, stdio transport), its own
  `requirements.txt` (`mcp`, `httpx`) and README. It is configured by env vars
  `CONSTRUCTOR_BASE_URL` and `CONSTRUCTOR_API_KEY` and talks **only** HTTP to
  `/api/v1/` (FR-009). Tools (7): `resolve_metabolite`, `create_reaction`,
  `save_reaction`, `list_reactions`, `get_reaction`, `fetch_vmh_reaction`,
  `fetch_rhea_reaction`. API errors are returned as MCP tool errors with the API's
  machine-readable code + message passed through.
- **Rationale**: stdio is the transport every MCP client supports today and needs no
  hosting change; a directory named `mcp/` would shadow the `mcp` package — hence
  `mcp_server/`. Keeping it outside `curationTool/` makes "never touches the ORM"
  structurally true (it cannot import Django app code without a path hack).
- **Alternatives considered**: Hosted streamable-HTTP MCP endpoint inside Django
  (deployment/auth complexity; can be added later on top of the same API); generating
  the server from the OpenAPI doc (less readable tool descriptions for assistants).

## R7. Testing strategy

- **Decision**: Convert `reactions/tests.py` (empty, 3 lines) into a `tests/` package.
  Django `TestCase` tests for: api-key generate/rotate/auth, metabolite resolve, reaction
  create/save/list/get, permission isolation, error envelope, VMH/Rhea fetch. External
  calls (VMH, ChEBI, Rhea) are mocked with `unittest.mock.patch`; balance-check tests use
  `smiles` inputs so RDKit computes everything offline. A new
  `reactions_project/test_settings.py` switches the DB to on-disk SQLite so tests and a
  local smoke server run without the Postgres container. MCP server gets an in-process
  smoke test using the SDK's stdio client (list tools, call one against a mocked/live
  API).
- **Rationale**: Postgres isn't guaranteed on dev machines (none running locally); the
  schema uses nothing Postgres-specific (TextFields + JSONField work on SQLite).
  Deterministic offline tests satisfy SC-002 comparisons without flaky upstreams.

## R8. Constitution gate

- **Decision**: Amend the constitution to v1.1.0 (MINOR): Principle I gains an explicit
  sanctioned exception for a versioned, token-authenticated programmatic API under
  `/api/v1/` plus a thin MCP wrapper, with the conditions that the web UI stays
  server-rendered, the API shares domain logic through `reactions/services/` + existing
  utils, and no API/MCP surface writes to VMH. User approved the API tier on 2026-11-02.
- **Rationale**: Principle I currently reads "not as a separate API/SPA tier"; building
  without amending would leave a standing violation that every future plan must
  re-justify. Governance section permits amendment with Sync Impact Report + version
  bump.

# Tasks: Curation API & MCP Server

**Input**: Design documents from `specs/002-curation-api-mcp/`
**Prerequisites**: plan.md, research.md, data-model.md, contracts/api-v1.md, contracts/mcp-tools.md, quickstart.md

**Tests**: Included — the spec's acceptance scenarios demand verdict-equivalence and
auth/isolation guarantees, and the constitution says unit-testable backend logic SHOULD
be tested. Upstreams (VMH/ChEBI/Rhea/PubChem) are mocked; balance tests use offline
SMILES inputs.

**Organization**: Tasks grouped by user story. US2 (resolve) is implemented before US1
(create) — both are P1, and US1 reuses US2's type-mapping groundwork; each remains
independently testable.

## Phase 1: Setup

- [X] T001 Add `django-ninja` to `requirements.txt`; create project venv (`.venv/`) and install requirements
- [X] T002 [P] Create `curationTool/reactions_project/test_settings.py` — imports `settings`, overrides DATABASES to on-disk SQLite, keeps everything else; verify `python manage.py migrate --settings=reactions_project.test_settings` runs
- [X] T003 [P] Replace empty `curationTool/reactions/tests.py` with package `curationTool/reactions/tests/__init__.py`

## Phase 2: Foundational (blocking all user stories)

- [X] T004 [P] Create `curationTool/reactions/services/__init__.py` and `curationTool/reactions/services/errors.py` — typed errors (`ServiceValidationError`, `UnknownIdentifierError`, `UnsupportedTypeError`, `DuplicateNameError`, `NotFoundError`, `UpstreamServiceError`, `ProcessingError`) carrying `code` + `detail` per data-model.md
- [X] T005 Add `ApiKey` model to `curationTool/reactions/models.py` (one-to-one User, `key_hash` unique, `prefix`, `created_at`, `last_used_at`) and generate additive migration in `curationTool/reactions/migrations/`
- [X] T006 Create `curationTool/reactions/services/api_keys.py` — `generate_key(user)` (returns raw key once, stores SHA-256, rotates if exists), `authenticate(raw_key) -> User | None` (digest lookup, update `last_used_at`)
- [X] T007 [P] Create `curationTool/reactions/api/schemas.py` — pydantic models: `MetaboliteInput` (type enum `vmh|chebi|chebi_name|pubchem|swisslipids|smiles|saved`, stoichiometry ≥1, compartment validated), `ReactionSubmission`, `ProcessedReaction`, `MetaboliteResolution`, `SaveReactionIn`, list/pagination models, `ErrorEnvelope`
- [X] T008 Create `curationTool/reactions/api/auth.py` — `ApiKeyAuth` accepting `Authorization: Bearer` and `X-API-Key`, delegating to `services.api_keys.authenticate`
- [X] T009 Create `curationTool/reactions/api/__init__.py` — `NinjaAPI(version="v1")` with `AuthRateThrottle` (rate from new `API_THROTTLE_RATE` setting read from env in `curationTool/reactions_project/settings.py`, default `60/m`), exception handlers mapping service errors → error envelope (per data-model.md table), `GET /health` (no auth), `GET /me`; mount at `/api/v1/` in `curationTool/reactions_project/urls.py`
- [X] T010 Foundational tests in `curationTool/reactions/tests/test_api_keys.py` — generate/rotate/authenticate unit tests; `/health` open; `/me` 401 without key, 200 with key, 401 after rotation with old key

**Checkpoint**: API skeleton auths and errors correctly — user story phases can begin.

## Phase 3: User Story 2 — Resolve a metabolite from an identifier (P1)

**Goal**: `POST /api/v1/metabolites/resolve` returns name/formula/charge/structure/VMH
presence for any supported identifier type.

**Independent test**: resolve a SMILES offline and a mocked ChEBI/VMH id; nonsense id →
422 naming supported types.

- [X] T011 [US2] Add `'SMILES'` branch + `smiles_to_mol` helper to `curationTool/reactions/utils/to_mol.py` (additive; same AddHs/property-cache normalisation as sibling converters)
- [X] T012 [US2] Create `curationTool/reactions/services/metabolite_resolution.py` — API-type → internal-type mapping, `resolve_metabolite(identifier, id_type)` using `any_to_mol(request=None)` + `get_mol_info` + `search_vmh`/`get_vmh_miriam`; raises typed errors (unknown id, unsupported type, upstream down)
- [X] T013 [US2] Create `curationTool/reactions/api/routes/metabolites.py` — `POST /metabolites/resolve` per contracts/api-v1.md; register router in `curationTool/reactions/api/__init__.py`
- [X] T014 [US2] Tests in `curationTool/reactions/tests/test_metabolites_api.py` — offline SMILES resolve (water: formula H2O, charge 0); mocked VMH abbreviation resolve; unknown identifier → 422 `unknown_identifier`; bad type → 422 listing supported types; upstream failure mock → 502 `upstream_unavailable`

**Checkpoint**: US2 fully functional and testable on its own.

## Phase 4: User Story 1 — Programmatic reaction creation with balance check (P1)

**Goal**: `POST /api/v1/reactions` processes + persists a reaction with web-identical
verdicts; `POST /api/v1/reactions/{id}/save` adds it to the caller's collection.

**Independent test**: balanced & unbalanced SMILES reactions return correct verdicts
offline; save appears in user's saved set (and web UI).

- [X] T015 [US1] Create `curationTool/reactions/services/reaction_processing.py` — `process_reaction(submission, user) -> (Reaction, ProcessedReaction-dict)`: side validation, `any_to_mol(request=None)`, per-metabolite error mapping (side+index), `get_mol_info`, `construct_reaction_rxnfile`, InChIKey signature, `get_reaction_info` (atom mapping always skipped, placeholder visualization like web skip-path), `check_reaction_vmh`, `search_metabolites_vmh(request=None)`, SavedMetabolite persistence for non-VMH metabolites (dedup by inchi_key, name-collision → `DuplicateNameError`, mutate types to `'Saved'`), assign + save all Reaction fields with the exact JSON-encoding `input_reaction` uses
- [X] T016 [US1] Add `save_reaction_for_user(user, reaction_id, name, description)` to `curationTool/reactions/services/reaction_processing.py` — ownership-safe fetch, per-user duplicate short_name check, set fields, add to `user.saved_reactions`
- [X] T017 [US1] Create `curationTool/reactions/api/routes/reactions.py` — `POST /reactions` (201) and `POST /reactions/{id}/save` per contract; register router
- [X] T018 [US1] Tests in `curationTool/reactions/tests/test_reactions_api.py` (VMH checks mocked) — balanced SMILES reaction → `balanced_count`/`balanced_charge` true; deliberately unbalanced → false with per-side atom tallies; unknown identifier → 422 naming substrate index; save → 200 and reaction in `user.saved_reactions`; duplicate name → 409; persisted row fields are valid JSON matching web encoding (spot-check `substrates_types`, `balanced_count`)

**Checkpoint**: MVP complete (US2 + US1): resolve → create/check → save, all via API.

## Phase 5: User Story 3 — Saved reactions & VMH/Rhea sources (P2)

**Goal**: list/get own reactions; import from VMH/Rhea in processed form.

**Independent test**: list/get round-trip; foreign id → 404; mocked VMH/Rhea fetch
returns processed reaction.

- [X] T019 [US3] Add `GET /reactions` (paginated, only caller's saved set) and `GET /reactions/{id}` (saved-or-created by caller; else `NotFoundError`) to `curationTool/reactions/api/routes/reactions.py`, with serializer helper reading the JSON-encoded Reaction fields
- [X] T020 [US3] Create `curationTool/reactions/services/external_sources.py` — `fetch_vmh_reaction(abbreviation)` (reuse `find_reaction_by_abbreviation` + `decode_formula`) and `fetch_rhea_reaction(rhea_id)` (reuse `get_from_rhea`), both building a `ReactionSubmission` and delegating to `process_reaction`
- [X] T021 [US3] Create `curationTool/reactions/api/routes/sources.py` — `POST /sources/vmh/fetch`, `POST /sources/rhea/fetch`; register router
- [X] T022 [US3] Tests in `curationTool/reactions/tests/test_sources_api.py` — mocked VMH fetch → 201 processed; unknown abbr → 422; Rhea mocked likewise; plus in `test_reactions_api.py`: list shows only own reactions, foreign `reaction_id` → 404 `not_found`

**Checkpoint**: Full API surface complete.

## Phase 6: User Story 4 — AI assistant curates via MCP (P2)

**Goal**: stdio MCP server exposing the 7 tools, HTTP-only against `/api/v1/`.

**Independent test**: SDK stdio client lists 7 tools; `resolve_metabolite` round-trips
against a live local server; API error surfaces as `Error [code]: message`.

- [X] T023 [US4] Create `mcp_server/server.py` — FastMCP server `constructor`; 7 tools per contracts/mcp-tools.md with assistant-grade docstrings (identifier types, compartments, "never writes to VMH"); env config `CONSTRUCTOR_BASE_URL`/`CONSTRUCTOR_API_KEY`; shared httpx client; error passthrough; `--smoke` flag that lists tools and exits
- [X] T024 [P] [US4] Create `mcp_server/requirements.txt` (`mcp`, `httpx`) and `mcp_server/README.md` (setup, env vars, client-config JSON snippet, tool table)
- [X] T025 [US4] Create `mcp_server/test_smoke.py` — spawns the server over stdio via the MCP SDK client: asserts 7 tools listed, asserts no tool name mentions VMH submission, calls `resolve_metabolite` against a running local API (skips gracefully if server not up); run it against the quickstart smoke server

**Checkpoint**: Assistant-ready; US4 acceptance scenarios executable.

## Phase 7: User Story 5 — Curator manages their access key (P3)

**Goal**: signed-in curator generates/rotates key in the web app; shown once.

**Independent test**: with session, page renders and generate returns key; old key then
401s; without session → redirected/denied.

- [X] T026 [US5] Add `api_key_page` (GET, renders template with prefix/created if key exists) and `api_key_generate` (POST, requires `request.session['userID']`, returns full key JSON once) to `curationTool/reactions/views/user_views.py`; routes `api-key/` and `api-key/generate/` in `curationTool/reactions/urls.py`
- [X] T027 [US5] Create `curationTool/reactions/templates/reactions/api_key.html` (extends existing layout, Semantic UI card/button/copy field, shown-once warning) and `curationTool/reactions/static/reactions/js/apiKey.js` (fetch + display + copy; file-per-concern)
- [X] T028 [US5] Extend `curationTool/reactions/tests/test_api_keys.py` — session-less POST rejected; with session: generate returns `ctor_` key, page shows prefix; rotation invalidates old key end-to-end through `/api/v1/me`

**Checkpoint**: All user stories complete.

## Phase 8: Polish & Cross-Cutting

- [X] T029 [P] Document the API + MCP in `README.md` (new section: base URL, auth, endpoint table, OpenAPI link, MCP setup pointer to `mcp_server/README.md`)
- [X] T030 Run full suite `python manage.py test reactions --settings=reactions_project.test_settings`; then quickstart live smoke sweep (migrate, runserver, register user, key page, curl health/me/reactions, OpenAPI docs page, MCP smoke) and fix anything found
- [X] T031 [P] Audit sweep: grep `reactions/api/` + `mcp_server/` for any `add_to_vmh`/VMH-write imports (must be none, FR-010/SC-005); confirm no existing file behaviour changed except additive edits (`git diff` review vs. plan structure)

## Dependencies & Execution Order

- **Phase 1 → Phase 2**: setup before foundation. T001 blocks everything; T002/T003 parallel.
- **Phase 2 before all stories**: T004→(T006,T007); T005→T006→T008→T009→T010.
- **US2 (Phase 3)**: independent after foundation. T011→T012→T013→T014.
- **US1 (Phase 4)**: after foundation; reuses US2's T011/T012 type mapping (soft dependency). T015→T016→T017→T018.
- **US3 (Phase 5)**: needs US1's T015 (processing) for fetch endpoints. T019 ∥ T020; then T021→T022.
- **US4 (Phase 6)**: needs the API surface (US1–US3) deployed locally. T023 ∥ T024; then T025.
- **US5 (Phase 7)**: only needs foundation (T006). Can run parallel to Phases 3–6.
- **Phase 8 last.**

## Parallel Example

```text
After T010 (foundation done):
  Agent A: T011–T014 (US2)      Agent B: T026–T028 (US5)
Then:
  Agent A: T015–T018 (US1)      Agent B: T023/T024 scaffolding (US4, mocked)
```

## Implementation Strategy

**MVP** = Phases 1–4 (foundation + US2 + US1): a key-authed API that resolves
metabolites and creates balance-checked, saveable reactions. Each later phase ships an
independently testable increment (US3 round-trips/imports, US4 MCP, US5 self-serve
keys). Stop-and-validate checkpoints after every phase; full acceptance sweep in T030.

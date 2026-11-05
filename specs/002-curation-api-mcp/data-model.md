# Data Model: Curation API & MCP Server

**Feature**: 002-curation-api-mcp | **Date**: 2026-11-02

## New entities

### ApiKey (new table, additive migration)

| Field | Type | Constraints | Notes |
|-------|------|-------------|-------|
| `user` | OneToOneField → `reactions.User` | `on_delete=CASCADE`, related_name `api_key` | One active key per curator (spec Assumption). Cascade satisfies the deleted-user edge case: key dies with the account. |
| `key_hash` | CharField(64) | unique, indexed | SHA-256 hex digest of the full key. Raw key never stored. |
| `prefix` | CharField(16) | — | First 12 chars of the key (`ctor_XXXXXXX`) for display/identification. |
| `created_at` | DateTimeField | auto_now_add | |
| `last_used_at` | DateTimeField | null=True | Updated (best effort) on successful auth. |

**Key format**: `ctor_` + `secrets.token_urlsafe(32)` (~48 chars). Shown exactly once at
generation/rotation (FR-003).

**State transitions**: none → generated → rotated (hash replaced; old key invalid
immediately) → deleted with user.

## Existing entities touched (no schema change)

- **User** (`reactions.models.User`): gains reverse `api_key`. No field changes.
- **Reaction**: created/read by the API with the *same field encoding* as the web flow
  (JSON-dumped lists in TextFields; `substrates_types` entries mutated to `'Saved'` when
  a non-VMH metabolite is persisted for the user — mirrors `input_reaction`). No schema
  change.
- **SavedMetabolite**: the processing service creates these for signed-in users exactly
  as the web flow does (same dedup by `inchi_key`, same name-collision error).
- **User.saved_reactions** (M2M): "save reaction" adds to this set after setting
  `short_name`/`description` — mirrors `save_user_reaction` semantics, including the
  duplicate-name check (FR-006, US1 scenario 4).

## Service-layer value objects (not persisted; pydantic schemas in `api/schemas.py`)

### MetaboliteInput
`identifier: str`, `type: vmh|chebi|chebi_name|pubchem|swisslipids|smiles|saved`,
`stoichiometry: int ≥ 1`, `compartment: str` (validated against the compartment list
used by the web UI), `name: str | None` (required for non-VMH types that create
SavedMetabolites, else derived).

### ReactionSubmission
`substrates: [MetaboliteInput]`, `products: [MetaboliteInput]` (at least one side
non-empty), `direction: forward|reversible`, `subsystem: str | None`,
`organs: [str] | None`. List-length consistency is inherent (per-item objects), so
FR-008's mismatch case is impossible by construction at the API (still validated in the
service for other callers).

### ProcessedReaction (response)
`reaction_id`, `molc_formula`, `balanced_count: {balanced: bool, detail}`,
`balanced_charge: {balanced: bool, subs_charge, prods_charge}`, `subs_atoms`,
`prods_atoms`, `symb_to_name`, per-metabolite: `names, formulas, charges, smiles,
inchis, inchi_keys, mol_weights, found_in_vmh, miriams`, `vmh: {found, similar, url,
formula}`, `direction`, `subsystem`.

### MetaboliteResolution (response)
`name, formula, charge, smiles, inchi, inchi_key, mol_weight, vmh: {found, abbr,
miriam}`.

### Error envelope
`{"error": {"code", "message", "detail?"}}` — codes enumerated in plan.md (FR-013).

## Typed service errors → HTTP

| Service error | HTTP | code |
|---------------|------|------|
| `ValidationError` (unknown type, bad compartment, empty reaction, bad stoich) | 422 | `validation_error` / `unsupported_type` |
| `UnknownIdentifierError` (names identifier + side + index) | 422 | `unknown_identifier` |
| `DuplicateNameError` (metabolite or reaction name) | 409 | `duplicate_name` |
| `NotFoundError` / foreign owner | 404 | `not_found` (no existence leak, FR-006) |
| `UpstreamServiceError` (VMH/ChEBI/Rhea/PubChem down) | 502 | `upstream_unavailable` |
| `ProcessingError` (RDKit/unexpected) | 500 | `processing_error` |
| Bad/missing key | 401 | `unauthorized` |
| Throttled | 429 | `rate_limited` |

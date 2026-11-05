# Contract: Constructor API v1

Base path: `/api/v1/` (Django Ninja; OpenAPI at `/api/v1/openapi.json`, interactive docs
at `/api/v1/docs`).

**Auth**: every endpoint except `GET /health` requires `Authorization: Bearer <key>`
(or `X-API-Key: <key>`). 401 `unauthorized` otherwise. Per-key throttle 60/min →
429 `rate_limited`.

**Errors**: non-2xx bodies are `{"error": {"code", "message", "detail?"}}`
(codes per data-model.md).

## Endpoints

### GET /api/v1/health
No auth. → `200 {"status": "ok", "version": "v1"}`

### GET /api/v1/me
→ `200 {"user_id", "name", "full_name", "key_prefix", "key_created_at"}`

### POST /api/v1/metabolites/resolve
Body: `{"identifier": "CHEBI:15377", "type": "chebi"}`
→ `200 MetaboliteResolution` (name, formula, charge, smiles, inchi, inchi_key,
mol_weight, vmh{found, abbr, miriam})
Errors: 422 `unsupported_type` (lists supported types), 422 `unknown_identifier`,
502 `upstream_unavailable`.

### POST /api/v1/reactions
Body: `ReactionSubmission`
```json
{
  "substrates": [{"identifier": "h2o", "type": "vmh", "stoichiometry": 1,
                   "compartment": "c", "name": null}],
  "products":   [{"identifier": "O", "type": "smiles", "stoichiometry": 1,
                   "compartment": "c", "name": "water (custom)"}],
  "direction": "forward",
  "subsystem": null,
  "organs": []
}
```
Processes (balance check + VMH match) and persists the reaction (not yet in the user's
saved list — mirrors web flow).
→ `201 ProcessedReaction` (includes `reaction_id`)
Errors: 422 `validation_error` / `unknown_identifier` (with side+index),
409 `duplicate_name` (saved-metabolite name clash), 502 `upstream_unavailable`.
Atom-mapping visualisation is never produced on this path (spec Assumption; balance
verdicts do not depend on it).

### POST /api/v1/reactions/{reaction_id}/save
Body: `{"name": "my hydration rxn", "description": ""}`
Sets short name/description and adds the reaction to the caller's saved list
(semantics of the web `save_reaction` flow).
→ `200 {"reaction_id", "name", "saved": true}`
Errors: 404 `not_found` (nonexistent), 409 `duplicate_name` (name already used by this
user).

### GET /api/v1/reactions?page=1&page_size=25
Lists only the caller's saved reactions.
→ `200 {"count", "page", "page_size", "results": [{"reaction_id", "name",
"rxn_formula", "balanced_mass": bool, "balanced_charge": bool, "vmh_found"}]}`

### GET /api/v1/reactions/{reaction_id}
Full `ProcessedReaction` view of one of the caller's reactions (saved or created by the
caller). Foreign/unknown id → 404 `not_found` (no existence leak).

### POST /api/v1/sources/vmh/fetch
Body: `{"abbreviation": "HEX1"}` → decodes the VMH formula and runs the same processing
pipeline. → `201 ProcessedReaction`.
Errors: 422 `unknown_identifier` (abbr not in VMH / no formula), 502
`upstream_unavailable`.

### POST /api/v1/sources/rhea/fetch
Body: `{"rhea_id": "10031"}` → same, from Rhea. → `201 ProcessedReaction`.
Same errors.

## Explicitly absent (FR-010)

No endpoint submits reactions or metabolites to VMH. The web-only VMH submission flow
is untouched and remains session-authenticated outside `/api/v1/`.

## Key management (web, session-authenticated — NOT under /api/v1/)

- `GET /api-key/` — Semantic UI page; shows prefix + created date if a key exists.
- `POST /api-key/generate/` — session user required (`request.session['userID']`);
  creates or rotates; returns the full key once (JSON), page shows it with a copy
  button and a "shown only once" warning.

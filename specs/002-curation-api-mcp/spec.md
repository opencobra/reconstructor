# Feature Specification: Curation API & MCP Server

**Feature Branch**: `002-curation-api-mcp`

**Created**: 2026-11-02

**Status**: Draft

**Input**: User description: "REST API (v1) and MCP server for Constructor reaction curation. Expose Constructor's core curation capabilities as a versioned, token-authenticated JSON API (/api/v1/) with a thin MCP server wrapper so AI agents and scripts can: resolve metabolites from identifiers (VMH, ChEBI, SwissLipids, PubChem), create and balance-check reactions (mass/charge), get reaction details, list/save a user's reactions, and search/fetch reactions from VMH and Rhea. Auth via per-user API keys. Push-to-VMH is explicitly out of scope for the MCP surface. Service layer extracted so template views and API share domain logic."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Programmatic reaction creation with balance check (Priority: P1)

A curator (or a script acting for them) submits a reaction — substrates and products given
as identifiers (VMH abbreviation, ChEBI ID, SwissLipids ID, PubChem CID, or SMILES) with
stoichiometry and compartments — and receives the fully processed reaction: resolved
metabolite names and formulas, mass-balance and charge-balance verdicts, and whether an
identical or similar reaction already exists in VMH. They can then save the reaction to
their personal collection.

**Why this priority**: Reaction creation + balance checking is Constructor's core value.
Making it scriptable unlocks batch curation (today each reaction requires manual entry in
the web UI) and is the foundation every other story builds on.

**Independent Test**: With a valid access key, POST a known-balanced reaction (e.g.
VMH `h2o` + `co2` → `h2co3`) and verify the response reports mass and charge balanced;
POST a deliberately unbalanced one and verify it reports the imbalance. Save one and see
it appear in the web UI's saved reactions.

**Acceptance Scenarios**:

1. **Given** a valid access key, **When** the curator submits a reaction whose metabolites
   are all valid VMH abbreviations, **Then** the response contains resolved names,
   molecular formulas, charge per metabolite, and mass/charge balance verdicts identical
   to what the web UI would show for the same input.
2. **Given** a reaction submission with an unknown metabolite identifier, **When** it is
   processed, **Then** the response is an error that names the offending identifier and
   its position (substrate/product index) and no reaction is stored.
3. **Given** a processed reaction, **When** the curator saves it with a name, **Then** it
   appears in their saved reactions both via the API and in the web UI.
4. **Given** a save request with a name the curator already used, **When** it is
   processed, **Then** the response flags the duplicate name and nothing is overwritten.

---

### User Story 2 - Resolve a metabolite from an identifier (Priority: P1)

A curator or agent looks up a single metabolite by identifier (VMH abbreviation, ChEBI,
SwissLipids, PubChem CID, or SMILES) and gets back its resolved name, molecular formula,
charge, structure (SMILES/InChI), and whether it is already present in VMH.

**Why this priority**: Metabolite resolution is the first step of every curation session
and is useful standalone (e.g. "is this ChEBI compound already in VMH?"). It exercises
the identifier-conversion machinery end to end with minimal surface.

**Independent Test**: Request resolution of a known ChEBI ID and verify name/formula/
charge match the web UI's resolution of the same input; request a nonsense ID and verify
a clear error.

**Acceptance Scenarios**:

1. **Given** a valid access key, **When** a known identifier is resolved, **Then** the
   response includes name, formula, charge, canonical structure, and a found-in-VMH flag
   (with the VMH abbreviation when found).
2. **Given** an identifier of an unsupported or mistyped type, **When** resolution is
   attempted, **Then** the response is an error stating the supported identifier types.

---

### User Story 3 - Work with saved reactions and VMH/Rhea sources (Priority: P2)

A curator lists their saved reactions, fetches full details of one, and imports a
reaction from VMH or Rhea by its identifier as a starting point for curation.

**Why this priority**: Round-tripping with the existing collection and external databases
makes the API usable for real workflows, but it depends on Story 1's processing for
imports to be useful.

**Independent Test**: Save a reaction via Story 1, then list and fetch it; fetch a known
VMH reaction abbreviation and a known Rhea ID and verify the returned reactions are
processed (balanced verdicts present).

**Acceptance Scenarios**:

1. **Given** a curator with saved reactions, **When** they list reactions, **Then** they
   receive their reactions (id, name, formula, balance status) and no other user's.
2. **Given** a saved reaction id belonging to another user, **When** details are
   requested, **Then** access is denied without revealing the reaction's existence.
3. **Given** a valid VMH reaction abbreviation or Rhea ID, **When** it is fetched,
   **Then** the response contains the processed reaction ready for balance review.

---

### User Story 4 - AI assistant curates via MCP (Priority: P2)

A curator connects an AI assistant (e.g. Claude) to Constructor via an MCP server
configured with their access key. The assistant can resolve metabolites, build and
balance-check reactions, fetch from VMH/Rhea, and list/save the curator's reactions —
conversationally, without the curator touching the web UI.

**Why this priority**: This is the headline capability of the feature, but it is a thin
wrapper: it only works once the API (Stories 1–3) exists, and it must never offer more
than the API does.

**Independent Test**: Configure the MCP server against a running instance with a test
key; from an MCP client, invoke each tool once and verify results match direct API calls.

**Acceptance Scenarios**:

1. **Given** a configured MCP connection, **When** the assistant lists available tools,
   **Then** each v1 capability (resolve metabolite, create/balance-check reaction, get
   reaction, list/save reactions, fetch VMH, fetch Rhea) appears with a description
   sufficient for the assistant to use it unaided.
2. **Given** a tool invocation that the API rejects (bad key, bad input), **When** the
   MCP tool returns, **Then** the error is passed through in a form the assistant can
   relay and act on (not a crash or silent failure).
3. **Given** the full MCP tool list, **When** inspected, **Then** no tool submits
   reactions or metabolites to VMH (the shared production database).

---

### User Story 5 - Curator manages their access key (Priority: P3)

A signed-in curator generates a personal access key from the web application, sees it
once, and can revoke/regenerate it at any time, immediately invalidating the old key.

**Why this priority**: Required for the other stories' auth, but minimal: one key per
user, generate/revoke, no scopes or expiry in v1.

**Independent Test**: Generate a key while signed in, call the API with it successfully,
regenerate, and verify the old key now fails.

**Acceptance Scenarios**:

1. **Given** a signed-in curator, **When** they request an access key, **Then** a key is
   shown to them and works for API calls attributed to their account.
2. **Given** a regenerated key, **When** the old key is used, **Then** the request is
   rejected as unauthorized.
3. **Given** no key or an invalid key, **When** any API operation is attempted, **Then**
   the request is rejected as unauthorized with no data revealed.

---

### Edge Cases

- A metabolite identifier is valid but the upstream source (ChEBI/PubChem/SwissLipids/
  VMH/Rhea) is unreachable or times out → the error distinguishes "bad input" from
  "upstream unavailable" so callers can retry.
- A reaction submission mixes identifier types across substrates/products → supported;
  each metabolite is resolved according to its own stated type.
- Mismatched list lengths (e.g. 3 substrates but 2 stoichiometries) → rejected with a
  message naming the mismatch, before any processing.
- An unknown compartment label → rejected with the list of accepted compartments.
- Reaction processing succeeds but the structure drawing/visualisation step fails → the
  balance verdicts are still returned; the response notes the missing visualisation.
- A very large batch (hundreds of sequential calls from a script) → per-key rate limiting
  protects the shared server; limited callers receive an explicit "slow down" signal, not
  a hang.
- An access key belonging to a since-deleted user → rejected as unauthorized.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST expose a versioned programmatic interface (v1) for curation
  operations, separate from the web pages, whose responses are structured data (not
  HTML).
- **FR-002**: System MUST authenticate every v1 operation via a per-user access key and
  attribute all resulting data (saved reactions, created metabolites) to that user.
- **FR-003**: A signed-in curator MUST be able to generate, view once, and regenerate
  their access key from the web application; regeneration invalidates the previous key.
- **FR-004**: System MUST resolve a single metabolite from any of: VMH abbreviation,
  ChEBI ID, SwissLipids ID, PubChem CID, or SMILES string — returning name, molecular
  formula, charge, canonical structure, and VMH presence.
- **FR-005**: System MUST process a submitted reaction (metabolites with types,
  stoichiometry, compartments, direction) and return mass-balance and charge-balance
  verdicts, per-side atom/charge tallies, resolved metabolite info, and VMH
  identical/similar matches — with results identical to the web UI for the same input.
- **FR-006**: System MUST let a curator save a processed reaction under a name, list
  their saved reactions, and fetch full details of one; curators MUST NOT be able to
  read or modify other users' reactions through v1.
- **FR-007**: System MUST fetch a reaction from VMH by abbreviation and from Rhea by ID,
  returning it in processed form (as FR-005).
- **FR-008**: System MUST reject malformed submissions (unknown identifier types,
  length-mismatched lists, unknown compartments, missing required fields) with errors
  that name the specific problem and offending field, before persisting anything.
- **FR-009**: The MCP server MUST expose exactly the v1 capabilities as tools (resolve
  metabolite, create/balance-check reaction, save reaction, list reactions, get
  reaction, fetch VMH, fetch Rhea) and MUST obtain all data through the v1 interface —
  never by direct database access.
- **FR-010**: Neither the v1 interface nor the MCP tools may submit reactions or
  metabolites to VMH; the existing web-only VMH submission flow is unchanged.
- **FR-011**: The web application's existing behaviour MUST be unchanged: the pages,
  their endpoints, and reaction semantics (stoichiometry, balance, compartments) remain
  as they are; v1 and the web UI share one implementation of the curation logic so
  verdicts cannot diverge.
- **FR-012**: v1 MUST apply per-key rate limiting with an explicit over-limit response.
- **FR-013**: Every v1 error response MUST carry a machine-readable error code and a
  human-readable message.

### Key Entities

- **Access Key**: A secret credential owned by exactly one curator; grants v1 access as
  that curator; can be regenerated (old value invalidated). One active key per curator.
- **Reaction (existing)**: The curated reaction record — metabolites, stoichiometry,
  compartments, direction, balance verdicts, VMH match info. v1 reads and creates these;
  it does not change their meaning or structure.
- **Curator / User (existing)**: The account that owns saved reactions and the access
  key.
- **Metabolite Resolution Result**: Transient answer to "what is this identifier?" —
  name, formula, charge, structure, VMH presence. Not persisted by resolution alone.
- **MCP Tool**: A named, described operation exposed to AI assistants; each maps
  one-to-one onto a v1 capability.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A curator with an access key can go from identifiers to a saved,
  balance-checked reaction in at most 3 programmatic calls (resolve optional →
  create/check → save).
- **SC-002**: For a representative set of at least 5 reactions (balanced and unbalanced,
  mixed identifier types), v1 balance verdicts and formulas match the web UI's results
  exactly, 100% of cases.
- **SC-003**: Reaction processing via v1 for a typical reaction (≤6 metabolites, no
  upstream fetch) completes in under 10 seconds; metabolite resolution of a local/VMH
  identifier in under 5 seconds.
- **SC-004**: An AI assistant connected via MCP completes the workflow "fetch reaction X
  from VMH and tell me if it is balanced" using only the provided tools, with zero
  manual steps besides approving tool use.
- **SC-005**: 100% of v1 and MCP operations are impossible without a valid access key,
  and 0 operations exist on either surface that write to VMH.
- **SC-006**: Existing web UI flows (create, save, list reactions) behave identically
  before and after the feature ships (no regressions in the existing test suite and
  spot-checked manually).

## Assumptions

- The existing account system (name/password users) is reused; no new identity provider
  is introduced. API keys attach to existing user accounts.
- One access key per user is sufficient for v1 (no scopes, expiry, or multiple keys).
- v1 exposes curation of reactions only; metabolite saving, templates, flags, gene info,
  organ data, leaderboard, and the tabbed workspace remain web-only for now.
- File-upload metabolite types (MDL/MOL drawings) are out of scope for v1; identifier
  and SMILES inputs cover the programmatic use cases.
- Atom-mapping visualisations (slow, external-tool dependent) are not required in v1
  responses; balance verdicts do not depend on them.
- The MCP server runs wherever the curator chooses (locally or alongside the app) and is
  configured with the Constructor base URL and the curator's access key.
- Rate limits are modest defaults suitable for a research tool (not a hard SLA) and can
  be tuned in deployment configuration.

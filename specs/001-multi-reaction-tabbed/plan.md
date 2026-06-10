# Implementation Plan: Multi-Reaction Tabbed Workspace

**Branch**: `001-multi-reaction-tabbed` | **Date**: 2026-06-10 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `specs/001-multi-reaction-tabbed/spec.md`

## Summary

Generalise the single-reaction main page into a **tabbed, multi-reaction workspace**. Each header
tab is a per-user **Reaction Group**; a group shows several reactions together; selecting one opens
it in the existing detail panels. Curators can add reactions from their saved set, remove them from a
group (without deleting them), and clone a reaction into a new editable variant. The core technical
move is replacing the **single global `reactionData`** with **per-reaction state keyed by an open-
reaction id**, plus new Django models/endpoints for groups and a header tab strip that swaps groups
without a full page reload.

## Technical Context

**Language/Version**: Python 3 (Django), client-side ES + jQuery (no transpiler)

**Primary Dependencies**: Django; Semantic UI 2.4 (CSS), jQuery 3.6, FontAwesome (CDN); existing
ChemDoodle, Reaction Decoder Tool (atom mapping), VMH/Rhea integrations, MATLAB HTTP bridge

**Storage**: Django ORM (existing DB); new tables for reaction groups + membership ordering

**Testing**: Django test layout (`reactions/tests.py`) for backend endpoints; manual acceptance
scenarios from spec.md for front-end (thin automated JS surface, per constitution)

**Target Platform**: Web (server-rendered Django templates), desktop browsers first

**Project Type**: Web application (Django monolith: server-rendered templates + static JS)

**Performance Goals**: Tab switch swaps the visible group with no full page reload; a group of ~10–20
reactions renders without perceptible layout jank

**Constraints**: No SPA framework, no front-end build step (Constitution I, II); reuse existing panels
and `display*.js` renderers; preserve reaction chemical/mathematical semantics (Constitution IV);
do not break the existing VMH-staging "workspace" feature

**Scale/Scope**: Single primary screen (`home_page.html`) + ~3–5 new endpoints + 1–2 new models +
several new/edited JS modules and CSS; per-user data

## Constitution Check

*GATE: must pass before design; re-check after design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Server-Rendered Django Core | PASS | Tabs/groups are Django models + views + template changes; AJAX for group switch, no SPA tier. |
| II. Modular Vanilla JS | PASS | New behaviour added as new file-per-concern JS modules (`tabsController.js`, `groupView.js`, `openReactionState.js`) reusing `WorkspacePanels` and `display*.js`. |
| III. UX Consistency | PASS | Tab strip and reaction cards use Semantic UI menu/card components and existing CSS idioms. |
| IV. Data & Integration Integrity | PASS | Reaction schema unchanged; new models reference `Reaction`; clone copies fields without altering semantics; migrations additive. |
| V. Spec Traceability | PASS | Every task below maps to an FR/SC in spec.md (see tasks.md coverage). |
| VI. Minimal Footprint | NEEDS DECISION | Reuse/extend existing `Workspace` model vs add a dedicated `ReactionGroup` model — see Structure Decision + Complexity Tracking. |

No CRITICAL violations. The one open decision (grouping model) is recorded below, not deferred.

## Project Structure

### Documentation (this feature)

```text
specs/001-multi-reaction-tabbed/
├── spec.md                  # Feature spec
├── plan.md                  # This file
├── tasks.md                 # Ordered tasks (/speckit-tasks output)
└── implementation-report.md # Milestone deliverable (planning/run report)
```

### Source Code (repository root)

```text
curationTool/reactions/
├── models.py                        # + ReactionGroup (or extend Workspace) + ordered membership
├── migrations/                      # + additive migration for grouping tables
├── urls.py                          # + routes: groups CRUD, add/remove/clone in group
├── views/
│   └── group_views.py               # NEW: group + membership + clone endpoints (JSON)
├── templates/reactions/
│   └── home_page.html               # header tab strip (~.top-menu) + multi-reaction container in #workspacePanels
└── static/reactions/
    ├── js/
    │   ├── openReactionState.js     # NEW: per-open-reaction state store (replaces single global reactionData)
    │   ├── tabsController.js        # NEW: header tab strip create/rename/delete/switch (AJAX)
    │   ├── groupView.js             # NEW: render group members as cards, select → open in panels
    │   ├── cloneReaction.js         # NEW: clone-with-edits flow
    │   ├── Creatediv.js             # EDIT: WorkspacePanels scoped to the selected open reaction
    │   ├── Displayalldivs.js        # EDIT: load a *selected* reaction's data via the state store
    │   └── display*.js              # EDIT (minimal): read from state store instead of global reactionData
    └── css/
        ├── home_page.css            # tab strip + group layout
        └── (reactants/reactioninfo.css …)  # per-card adjustments as needed
```

**Structure Decision (grouping model)**: Introduce a **dedicated `ReactionGroup` model** (owner FK to
`User`, `name`, ordered M2M / through-model to `Reaction`, `position`/`is_active`) rather than
overloading the existing `Workspace` model. Rationale: `Workspace` is a single per-user
(`OneToOneField`) bag that already backs the **VMH-staging** page (`send_to_workspace` /
`remove_from_workspace` / `vmh_workspace` in `reactions/views/vmh_views.py`); curators need *many*
named groups, and reusing `Workspace` would (a) collide semantically with VMH staging and (b) require
breaking its one-per-user shape. A new model keeps both features intact (Constitution IV) at the cost
of one added model (justified in Complexity Tracking).

## Design Notes

### State refactor (FR-003, SC-006) — the central change

- Today, a saved reaction is loaded by `Displayalldivs.js` into a **single global `reactionData`**
  consumed by `reactantsdisplay.js`, `displayChemInfoDiv.js`, `displayMetaboliteInfo.js`, etc.
- Introduce `openReactionState.js`: a store keyed by `openReactionId` (the selected group member or a
  clone's temp id), each entry holding the same shape today's `reactionData` carries. Provide
  `getActiveReaction()` / `setActiveReaction(id)` accessors. The `display*.js` modules change from
  reading `reactionData` to reading `getActiveReaction()`. This is a mechanical, low-semantic-risk
  edit but touches many modules — call it out as the highest-effort task cluster.

### Header tabs (FR-005–FR-009)

- Add a Semantic UI tab/menu strip inside `.top-menu` ([home_page.html:209](../../curationTool/reactions/templates/reactions/home_page.html#L209)),
  alongside the existing logo/About nav and user dropdown.
- `tabsController.js` handles create/rename/delete/switch via AJAX to `group_views.py`; switching
  fetches the target group's membership and re-renders the group view (no page reload, FR-006).

### Group view + detail panels (FR-001, FR-002)

- Restructure `#workspacePanels` ([home_page.html:319](../../curationTool/reactions/templates/reactions/home_page.html#L319)) into a
  **master/detail**: a group list/grid (reaction cards: short name, formula, balance status) plus the
  existing detail panels for the selected reaction. Reuse `WorkspacePanels` (`Creatediv.js`) but scope
  panel visibility to the active open reaction. `displaysavedReaction.js` patterns are reused for
  rendering card summaries.

### Add / remove from saved (FR-010–FR-012)

- Reuse the existing saved-reactions retrieval; add an "add from saved" picker (Semantic UI modal,
  mirroring existing modals in `home_page.html`). Endpoints mirror the existing
  `send_to_workspace`/`remove_from_workspace` shape but target `ReactionGroup`. Duplicate-guard in the
  endpoint and the UI.

### Clone with edits (FR-013–FR-015)

- `cloneReaction.js`: deep-copy the selected reaction's fields into a new unsaved open-reaction state
  entry (temp id), add it to the group view as "new/unsaved", let the user edit via existing panels,
  and persist via the existing save path (creating a new `Reaction` row). Original row untouched
  (FR-014). Discard removes the temp entry without persisting (FR-015).

### Safety & fallback (FR-017, FR-018)

- Unsaved-edit guard on tab switch / remove: track a dirty flag per open reaction; warn-or-save/
  discard (reuse existing modal patterns).
- No logged-in user → render a single ephemeral default tab and the current single-reaction behaviour
  (SC-005); group persistence endpoints are no-ops/anonymous-safe.

### Migration path (no regression)

- The current single-reaction page becomes the **one-group / one-selected-reaction** case. Seed a
  default group per user on first load so existing users land in a working tab. The
  `python manage.py runserver` smoke path must still create a brand-new reaction with an empty group.

### Validation strategy

- Backend: Django tests for group CRUD, add/remove (no underlying delete), and clone (new row +
  original unchanged).
- Front-end: manual acceptance scenarios SC-001…SC-006 against a running server; explicitly verify
  the VMH-staging workspace page still functions (regression check for the model decision).

## Complexity Tracking

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|--------------------------------------|
| New `ReactionGroup` model (vs reuse `Workspace`) | Curators need many named, ordered groups per user; `Workspace` is `OneToOne` and already owns VMH-staging semantics | Reusing `Workspace` would break its one-per-user shape and collide with the VMH-staging feature, risking Constitution IV (data/integration integrity) |

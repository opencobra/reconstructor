---
description: "Task list for Multi-Reaction Tabbed Workspace"
---

# Tasks: Multi-Reaction Tabbed Workspace

**Input**: Design documents from `specs/001-multi-reaction-tabbed/`

**Prerequisites**: plan.md (required), spec.md (required for user stories)

**Tests**: Backend tasks include Django tests (group CRUD / add-remove / clone are persistence-
critical). Front-end is validated via the spec's manual acceptance scenarios.

**Organization**: Grouped by user story so each story is an independently testable increment.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: can run in parallel (different files, no dependency)
- **[Story]**: US1–US4, or FOUND (foundational), SETUP, POLISH
- File paths are relative to repo root.

## Phase 1: Setup (Shared)

- [ ] T001 [SETUP] Confirm feature branch `001-multi-reaction-tabbed` and spec-kit artifacts present; review `.specify/memory/constitution.md` gates before coding.
- [ ] T002 [SETUP] Snapshot current behaviour: load `home_page.html` and the VMH-staging `VMH_Workspace` page on a running server to baseline what must not regress.

## Phase 2: Foundational (Blocking Prerequisites)

**⚠️ Must complete before US1–US4. Centred on the state refactor + grouping model.**

- [ ] T003 [FOUND] Add `ReactionGroup` model (owner FK, `name`, `is_active`) + ordered through-model membership to `Reaction` (`position`) in `curationTool/reactions/models.py`. *(FR-007, FR-009; Key Entity: Reaction Group)*
- [ ] T004 [FOUND] Generate additive migration for the grouping tables in `curationTool/reactions/migrations/`; verify it does not touch `Reaction` or `Workspace`. *(Constitution IV)*
- [ ] T005 [FOUND] Seed/get-or-create a default group per user on first main-page load. *(FR-008, SC-005)*
- [ ] T006 [FOUND] Create `curationTool/reactions/static/reactions/js/openReactionState.js`: per-open-reaction state store keyed by `openReactionId`, with `getActiveReaction()`/`setActiveReaction(id)`/`setDirty(id)` accessors, mirroring the current `reactionData` shape. *(FR-003, SC-006)*
- [ ] T007 [FOUND] Migrate `Displayalldivs.js` and the `display*.js` renderers (`reactantsdisplay.js`, `displayChemInfoDiv.js`, `displayMetaboliteInfo.js`, `displayAtomMapping.js`, `status.js`, `reactioninfoformhandling.js`) to read from `getActiveReaction()` instead of the global `reactionData`. *(FR-003, SC-006)* — highest-effort cluster.
- [ ] T008 [FOUND] Scope `WorkspacePanels` in `Creatediv.js` to the active open reaction so panel show/hide acts on the selected group member. *(FR-002)*

**Checkpoint**: per-reaction state + grouping model exist; the page still works for one reaction.

## Phase 3: User Story 1 — View a group on one page (P1) 🎯 MVP

**Goal**: Render the current group's reactions together; select one to edit in the panels.

- [ ] T009 [US1] Add `GET` group-contents endpoint in `curationTool/reactions/views/group_views.py` returning the active group's ordered members. *(FR-001)*
- [ ] T010 [US1] Wire the route in `curationTool/reactions/urls.py`. *(FR-001)*
- [ ] T011 [US1] Restructure `#workspacePanels` in `templates/reactions/home_page.html` into master (group list/grid) + detail (existing panels). *(FR-001, FR-002)*
- [ ] T012 [P] [US1] Create `static/reactions/js/groupView.js`: render reaction cards (short name, formula, balance status) reusing `displaysavedReaction.js` summary patterns; selecting a card calls `setActiveReaction()` and loads the detail panels. *(FR-001, FR-002)*
- [ ] T013 [P] [US1] Add group/card + master-detail styles in `static/reactions/css/home_page.css`. *(Constitution III)*
- [ ] T014 [US1] Empty-group state that still allows Create Reaction (default behaviour preserved). *(FR-004, SC-005)*
- [ ] T015 [US1] Django test: group-contents endpoint returns ordered members for a user. *(SC-001)*

**Checkpoint**: US1 independently testable — group of ≥10 reactions renders; any one edits/saves.

## Phase 4: User Story 2 — Header tabs / groups (P2)

**Goal**: Header tab strip to create/rename/delete/switch groups without page reload.

- [ ] T016 [US2] Add group CRUD endpoints (create/rename/delete) in `group_views.py` + routes in `urls.py`; delete removes the group, not its reactions. *(FR-009)*
- [ ] T017 [US2] Add the Semantic UI tab strip inside `.top-menu` in `home_page.html`. *(FR-005)*
- [ ] T018 [US2] Create `static/reactions/js/tabsController.js`: render tabs, mark active, and on click fetch + render the target group via `groupView.js` (AJAX, no reload). *(FR-005, FR-006)*
- [ ] T019 [US2] Persist active group and membership per user; restore on load; ensure ≥1 default tab always present. *(FR-007, FR-008)*
- [ ] T020 [P] [US2] Tab-strip styles in `home_page.css`. *(Constitution III)*
- [ ] T021 [US2] Django test: create/rename/delete group; delete does not delete member reactions. *(FR-009)*

**Checkpoint**: US2 independently testable — switch tabs, membership/selection preserved, no reload (SC-002).

## Phase 5: User Story 3 — Add / remove from saved (P2)

**Goal**: Add reactions from Saved Reactions into a group; remove without deleting.

- [ ] T022 [US3] Add `add_to_group` / `remove_from_group` endpoints in `group_views.py` + routes, mirroring existing `send_to_workspace`/`remove_from_workspace` in `vmh_views.py`; duplicate-guard on add. *(FR-010, FR-011, FR-012)*
- [ ] T023 [US3] "Add from saved" picker modal in `home_page.html` (reuse existing modal/`displaysavedReaction.js` list patterns). *(FR-010)*
- [ ] T024 [US3] Wire add/remove actions in `groupView.js`; on remove, drop from group view but keep in Saved Reactions. *(FR-011)*
- [ ] T025 [US3] Django test: add persists membership across reload; remove leaves the saved reaction intact; duplicate add is prevented. *(SC-003, FR-012)*

**Checkpoint**: US3 independently testable — add/remove behaves; no underlying data loss.

## Phase 6: User Story 4 — Clone with edits (P3)

**Goal**: Clone a reaction into a new editable variant; original unaffected.

- [ ] T026 [US4] Create `static/reactions/js/cloneReaction.js`: deep-copy active reaction's state into a new temp open-reaction entry marked unsaved/new; add to group view. *(FR-013)*
- [ ] T027 [US4] Reuse the existing save path so saving a clone creates a **new** `Reaction` row; ensure original row's persisted fields are untouched. *(FR-014, FR-015)*
- [ ] T028 [US4] Discard-unsaved-clone path removes the temp entry without persisting. *(FR-015)*
- [ ] T029 [US4] Django test: clone → change one metabolite + GPR → save yields a distinct reaction; original unchanged. *(SC-004)*

**Checkpoint**: US4 independently testable — clone-with-edits produces a distinct reaction.

## Phase 7: Cross-cutting safety & polish

- [ ] T030 [POLISH] Dirty-flag guard: warn-or-save/discard on tab switch and on remove (reuse existing modal patterns). *(FR-017)*
- [ ] T031 [POLISH] Not-logged-in fallback: single ephemeral tab + current single-reaction behaviour, no errors. *(FR-018, SC-005)*
- [ ] T032 [POLISH] Large-group layout check (≥20 reactions): scroll/paging, no layout breakage. *(Edge Cases)*
- [ ] T033 [POLISH] Regression check: VMH-staging `VMH_Workspace` page still works after the model decision. *(Constitution IV)*
- [ ] T034 [POLISH] Accessibility pass: tabs and cards keyboard-operable / screen-reader-reasonable. *(Constitution III)*

## Dependencies

- Phase 2 (T003–T008) blocks all user stories.
- US1 (Phase 3) is the MVP and precedes US2/US3 UI work that renders into the group view.
- US2/US3 can be built in parallel after US1; US4 (Phase 6) depends on US1 detail-editing + state store.
- Polish (Phase 7) follows the stories it guards.

## Requirement → Task Coverage

| Requirement | Tasks |
|-------------|-------|
| FR-001 | T009, T011, T012, T015 |
| FR-002 | T008, T011, T012 |
| FR-003 | T006, T007 |
| FR-004 | T014 |
| FR-005 | T017, T018 |
| FR-006 | T018 |
| FR-007 | T003, T019 |
| FR-008 | T005, T019 |
| FR-009 | T016, T021 |
| FR-010 | T022, T023 |
| FR-011 | T022, T024 |
| FR-012 | T022, T025 |
| FR-013 | T026 |
| FR-014 | T027, T029 |
| FR-015 | T027, T028 |
| FR-016 | (whole-feature invariant; enforced via T029/T033 + Constitution IV) |
| FR-017 | T030 |
| FR-018 | T031 |
| SC-001 | T015 |
| SC-002 | T018, T019 |
| SC-003 | T025 |
| SC-004 | T029 |
| SC-005 | T005, T014, T031 |
| SC-006 | T006, T007 |

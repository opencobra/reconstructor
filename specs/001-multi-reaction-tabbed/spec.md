# Feature Specification: Multi-Reaction Tabbed Workspace

**Feature Branch**: `001-multi-reaction-tabbed`

**Created**: 2026-06-10

**Status**: Draft

**Input**: User description: "Refactor the main page so a curator can work with multiple reactions at once (a group), switch between groups via header tabs, quickly add/remove reactions from the group using their saved reactions, and clone an existing reaction with small edits (different metabolite, different GPR, etc.). Professional, user-friendly UI/UX."

## Overview

Today the Constructor main page (`home_page.html`) is a **single-reaction workspace**: the
curator creates, edits, or views exactly one reaction at a time across the side-nav panels
(Reactants, Atom Mapping, Chem Info, Metabolite Info, References, Ext. Links, Gene Info,
Comments). All front-end rendering reads one global `reactionData` object.

This feature generalises the workspace to operate on **groups of reactions** presented under
**tabs in the header**, while preserving every existing single-reaction capability. The intent
is faster batch curation: assemble a working set, pull in saved reactions, and spin up variants
by cloning-with-edits — without losing the rich per-reaction editing the tool already provides.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - View a group of reactions on one page (Priority: P1)

A curator opens the main page and sees the reactions in their current group rendered together
(e.g. a list/grid of reaction cards), rather than a single reaction. Selecting a reaction in the
group opens it in the familiar detailed panels for editing.

**Why this priority**: This is the foundational shift from single-reaction to group. Without it,
none of the other stories (tabs, add/remove, clone) have a surface to act on. It is the MVP.

**Independent Test**: With a user who has several saved reactions seeded into one group, load the
main page and confirm all group members render together and any one can be opened into the detail
panels and edited/saved exactly as today.

**Acceptance Scenarios**:

1. **Given** a logged-in curator whose current group contains 3 reactions, **When** they open the
   main page, **Then** all 3 reactions are visible together in the workspace with their short
   name / formula, and no reaction's data bleeds into another's panels.
2. **Given** the group view, **When** the curator selects one reaction, **Then** that reaction
   loads into the Reactants/Chem Info/… panels and can be edited and saved using the existing
   controls.
3. **Given** an empty group, **When** the curator opens the main page, **Then** they see an empty
   state that still allows creating a brand-new reaction (current default behaviour preserved).

---

### User Story 2 - Switch between groups via header tabs (Priority: P2)

A tab controller in the header lets the curator keep several groups open and switch between them.
Each tab corresponds to one group of reactions; switching tabs swaps the visible group without a
full page reload.

**Why this priority**: Tabs are the organising metaphor the milestone calls for and make multi-
group curation usable, but they depend on the group view (US1) existing first.

**Independent Test**: Create two groups, switch between their tabs, and confirm each tab shows its
own membership and that switching does not reload the page or lose unsaved edits in the other tab
beyond the documented rules.

**Acceptance Scenarios**:

1. **Given** two groups exist, **When** the curator clicks a tab in the header, **Then** the
   workspace shows that group's reactions and the active tab is visually indicated.
2. **Given** the curator is on Tab A, **When** they switch to Tab B and back, **Then** Tab A's
   membership and selection are preserved (persisted per user).
3. **Given** a curator with no groups, **When** the page loads, **Then** a single default tab/group
   is presented so the page is never tab-less.

---

### User Story 3 - Quick add / remove reactions from the group (Priority: P2)

From within a group the curator can quickly add reactions from their **saved reactions** and remove
reactions from the group. Removing from a group does not delete the underlying saved reaction.

**Why this priority**: This is what makes a group a *working set* the curator actively shapes. It
depends on US1 and reuses the existing saved-reactions data.

**Independent Test**: From a group, open the add-from-saved picker, add a saved reaction, confirm it
appears in the group and persists; then remove it and confirm it leaves the group but still exists
in the user's saved reactions.

**Acceptance Scenarios**:

1. **Given** the curator has saved reactions, **When** they use "add from saved" and pick one,
   **Then** it is added to the current group and the membership persists across reload.
2. **Given** a reaction in the group, **When** the curator removes it from the group, **Then** it
   disappears from the group but remains in the user's Saved Reactions list (no data loss).
3. **Given** a saved reaction already in the group, **When** the curator tries to add it again,
   **Then** the system prevents a duplicate or clearly indicates it is already present.

---

### User Story 4 - Clone a reaction with edits (Priority: P3)

The curator can clone an existing reaction in the group into a new, independent reaction and make
small changes to the clone (e.g. swap a metabolite, change a compartment, edit the GPR) before
saving it. The original is unaffected.

**Why this priority**: High-value accelerator for creating reaction variants, but it builds on the
group + detail-editing surface from US1–US3, so it comes last.

**Independent Test**: Clone a reaction in the group, change one metabolite and the GPR on the clone,
save it, and confirm the clone is a distinct saved reaction while the original is unchanged.

**Acceptance Scenarios**:

1. **Given** a reaction in the group, **When** the curator chooses "clone", **Then** a new editable
   reaction pre-filled from the original appears in the group, marked as unsaved/new.
2. **Given** a clone, **When** the curator swaps a substrate/product metabolite or edits the GPR and
   saves, **Then** the clone persists as a separate reaction and the original reaction's data is
   unchanged.
3. **Given** a clone that has not yet been edited or saved, **When** the curator discards it,
   **Then** no new reaction is persisted.

### Edge Cases

- **Large groups**: a group with many reactions must remain usable (no layout breakage); paging or
  scrolling behaviour is defined in the plan.
- **Unsaved edits on tab switch / remove**: switching tabs or removing a reaction with unsaved edits
  must warn or follow a documented save/discard rule so work is not silently lost.
- **Anonymous / not-logged-in users**: group persistence requires a user; the page must degrade to
  the current single-reaction behaviour (or a single ephemeral tab) when no user is logged in.
- **Concurrent edits**: the same saved reaction appearing in a group while also edited elsewhere —
  define which state wins on save.
- **Clone of an unbalanced / incomplete reaction**: cloning must copy state faithfully, including
  "not found" / unbalanced indicators, without silently "fixing" them.

## Requirements *(mandatory)*

### Functional Requirements

**Group workspace (US1)**

- **FR-001**: The main page MUST render the reactions belonging to the curator's current group
  together in one view, instead of only a single reaction.
- **FR-002**: Selecting a reaction in the group MUST open it in the existing detail panels
  (Reactants, Atom Mapping, Chem Info, Metabolite Info, References, Ext. Links, Gene Info, Comments)
  for editing and saving, reusing current behaviour.
- **FR-003**: The system MUST keep each reaction's state isolated so that editing one group member
  never corrupts another's data (replacing the single global `reactionData` with per-reaction
  state).
- **FR-004**: All existing single-reaction workflows (Create Reaction, Save Reaction, Reset, skip
  atom mapping) MUST remain available and behave as today for the selected reaction.

**Header tabs / groups (US2)**

- **FR-005**: The header MUST present a tab controller where each tab corresponds to one group of
  reactions, with the active tab visually indicated.
- **FR-006**: Switching tabs MUST swap the visible group without a full page reload.
- **FR-007**: Groups and their membership MUST persist per user across sessions.
- **FR-008**: The system MUST always present at least one (default) group/tab, so the workspace is
  never without a group.
- **FR-009**: Users MUST be able to create, rename, and delete groups/tabs; deleting a group MUST
  NOT delete the underlying saved reactions it referenced.

**Add / remove from saved (US3)**

- **FR-010**: Users MUST be able to add reactions from their Saved Reactions into the current group.
- **FR-011**: Users MUST be able to remove a reaction from a group, and removal MUST NOT delete the
  underlying saved reaction.
- **FR-012**: The system MUST prevent (or clearly flag) adding a reaction that is already a member of
  the current group.

**Clone with edits (US4)**

- **FR-013**: Users MUST be able to clone a reaction in the group into a new, independent, editable
  reaction pre-filled from the original.
- **FR-014**: Edits to a clone (including swapping a metabolite, changing a compartment, or editing
  the GPR) MUST NOT modify the original reaction.
- **FR-015**: A saved clone MUST be persisted as a distinct reaction; an unsaved clone that is
  discarded MUST NOT be persisted.

**Cross-cutting**

- **FR-016**: The feature MUST preserve reaction chemical/mathematical semantics (substrates,
  products, stoichiometry, charge/mass balance, compartments); grouping/tabbing/cloning are
  organisational and MUST NOT alter reaction meaning. Cone/solver conventions use `x` for primal and
  `s` for dual where relevant.
- **FR-017**: Unsaved edits MUST be protected on tab switch and on remove via a documented
  warn-or-save/discard rule (no silent loss).
- **FR-018**: When no user is logged in, the page MUST degrade to existing single-reaction behaviour
  (or a single ephemeral group) without error.

### Key Entities *(include if feature involves data)*

- **Reaction Group (Tab)**: a named, per-user, ordered collection of reactions shown under one header
  tab. Attributes: owner (user), name, ordered membership of reactions, optional "active/last-opened"
  marker. NOTE: the codebase already has a `Workspace` model (`user → reactions` M2M) that today backs
  the separate "send to VMH" staging page; the plan MUST decide whether to extend/rename that concept
  or introduce a dedicated grouping model, and MUST avoid a naming collision between this in-page
  "group/tab" and the existing VMH-staging "workspace".
- **Reaction**: existing `Reaction` model; unchanged in meaning. Group membership and clones reference
  `Reaction` rows; clones are new `Reaction` rows.
- **Saved Reactions**: existing per-user saved set, the source for "add from saved".

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A curator can view a group of at least 10 reactions on the main page together, open any
  one into the detail panels, and edit/save it, with no cross-reaction data corruption (US1).
- **SC-002**: Switching between two header tabs swaps the visible group **without a full page reload**
  and preserves each tab's membership and selection (US2).
- **SC-003**: Adding a saved reaction to a group and reloading the page shows the reaction still in the
  group; removing it from the group leaves it present in Saved Reactions (US3).
- **SC-004**: Cloning a reaction, changing one metabolite and the GPR, and saving yields a distinct
  saved reaction while the original is byte-for-byte unchanged in its persisted fields (US4).
- **SC-005**: With no logged-in user, the main page loads and a new reaction can still be created
  (no regression of default behaviour) (FR-018).
- **SC-006**: No workflow path depends on the single global `reactionData` after the refactor; each
  open reaction has isolated state (FR-003).

## Assumptions

- Group membership and tabs are **per-user** and require a logged-in user to persist; anonymous use
  falls back to the current single-reaction experience.
- "Saved reactions" refers to the user's existing saved-reaction set; no new sharing/permission model
  is introduced by this feature.
- The detail-editing panels and their `display*.js` renderers are reused as-is for a *selected*
  reaction; the refactor changes how state is held and how reactions are listed/selected, not the
  per-reaction editing semantics.
- Mobile/small-screen optimisation of the multi-reaction layout is desirable but may be staged after
  the desktop experience.
- The existing VMH-staging "workspace" feature continues to function; this feature must not break it.

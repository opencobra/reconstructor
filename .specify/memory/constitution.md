# Constructor (Reconstruction Curation Tool) Constitution

<!--
Sync Impact Report
- Version: (none) → 1.0.0  (initial ratification)
- Principles defined: I Server-Rendered Django Core; II Modular Vanilla JS;
  III UX Consistency (Semantic UI); IV Data & Integration Integrity;
  V Spec Traceability; VI Minimal Footprint
- Added sections: Technology Constraints; Development Workflow & Quality Gates
- Templates reviewed: plan-template.md (Constitution Check gate ✅ aligns),
  spec-template.md ✅, tasks-template.md ✅
- Follow-up TODOs: none
-->

## Core Principles

### I. Server-Rendered Django Core

The application is a Django, server-rendered web app. Pages are produced by Django
views rendering templates under `curationTool/reactions/templates/`; data lives in the
Django ORM (`curationTool/reactions/models.py`). New behaviour MUST be expressed as
Django views + URL routes + templates, not as a separate API/SPA tier. Introducing a
client-side framework (React, Vue, etc.) or a separate frontend service is a breaking
architectural change and MUST be justified in the plan's Complexity Tracking and
approved before work begins.

### II. Modular Vanilla JS (One Concern Per File)

Front-end interactivity is plain ES + jQuery, organised as **one module per concern**
under `curationTool/reactions/static/reactions/js/` (e.g. `Creatediv.js`,
`Displayalldivs.js`, `savereaction.js`). There is **no build step / bundler / npm
toolchain** for app code; scripts are loaded directly in templates. New front-end logic
MUST follow the existing file-per-concern convention and reuse established helpers
(`WorkspacePanels`, the `display*.js` renderers, `displaysavedReaction.js`) rather than
introducing parallel mechanisms. A build step or framework MAY be proposed only with an
explicit migration rationale in the plan.

### III. UX Consistency (Semantic UI)

The interface uses Semantic UI components, FontAwesome icons, and the established
header / side-nav / workspace-panel layout. New UI MUST match these patterns —
component styles, button/label/icon idioms, and the existing CSS files in
`static/reactions/css/` — so the product stays visually coherent and professional.
Interactions MUST remain keyboard-operable and screen-reader-reasonable. Net-new
visual paradigms require a stated UX rationale.

### IV. Data & Integration Integrity

The `Reaction` schema and the external integrations — VMH, Rhea, ChemDoodle, the
MATLAB/HTTP bridge, atom mapping (RDT) — encode domain meaning. Changes MUST NOT alter
the mathematical or chemical semantics of a reaction (substrates/products,
stoichiometry, charge/mass balance, compartments) as a side effect of a UI or workflow
change. Database migrations MUST be additive and reviewed; destructive migrations
require explicit justification and a backup/restore note.

### V. Spec Traceability (NON-NEGOTIABLE)

Every implementation decision MUST be auditable back to a functional requirement (FR)
or success criterion (SC) in the active feature `spec.md`. Plans and tasks that
introduce behaviour absent from the spec MUST trigger a spec update first (the
analyze step exists to catch this). No `/speckit-implement` run proceeds while a
CRITICAL traceability gap is open.

### VI. Minimal Footprint

Prefer extending existing utilities, models, views, and templates over adding new
files. New helper modules, models, or endpoints are introduced only when an existing
one cannot reasonably carry the behaviour, and the plan MUST say why. Keep diagnostics
and added UI concise; avoid speculative generality (YAGNI).

## Technology Constraints

- **Backend**: Python 3 / Django (app package `curationTool`, primary app
  `reactions`). ORM-backed models; views split by concern under `reactions/views/`.
- **Frontend**: HTML templates + Semantic UI (CSS) + jQuery + vanilla JS modules +
  FontAwesome, all via the template `<script>`/`<link>` includes (CDN or
  `static/`). No app-level build pipeline.
- **External systems**: VMH and Rhea (reaction/metabolite sources), ChemDoodle
  (drawing), Reaction Decoder Tool (atom mapping), and a MATLAB integration over HTTP
  (see `MATLAB_HTTP_MIGRATION.md`). Treat these as integration boundaries whose
  contracts must be preserved.
- **Deployment**: Docker / gunicorn (see `docker/`, `gunicorn.conf.py`,
  `deploy.sh`). Changes must remain runnable under the existing container setup.

## Development Workflow & Quality Gates

- **Spec-Driven flow**: features follow `constitution → specify → (clarify) → plan →
  tasks → (analyze) → implement`, with `spec.md`, `plan.md`, `tasks.md` committed under
  `specs/NNN-feature/`.
- **Constitution Check gate**: each `plan.md` MUST include a Constitution Check
  asserting compliance with Principles I–VI (or recording justified exceptions in
  Complexity Tracking).
- **Analyze before implement**: run the analyze step after tasks; resolve CRITICAL
  findings (especially Principle V traceability) before implementing.
- **Validation**: the app has a thin automated-test surface, so feature acceptance is
  primarily via the per-user-story Acceptance Scenarios in `spec.md`, exercised
  manually against a running server (`python manage.py runserver`). Backend logic that
  can be unit-tested SHOULD be, under the app's Django test layout.
- **Commits**: use the repo's conventional prefixes (`feat:`, `fix:`, `refactor:`,
  `docs:`).

## Governance

This constitution supersedes ad-hoc practice for this repository. Amendments are made by
editing this file with a Sync Impact Report and a version bump per semantic versioning:
**MAJOR** for incompatible governance/principle removals or redefinitions, **MINOR** for
a new principle or materially expanded guidance, **PATCH** for clarifications. Plans and
reviews MUST verify compliance; unavoidable deviations are recorded in the plan's
Complexity Tracking with a rationale and a simpler-alternative-rejected note.

**Version**: 1.0.0 | **Ratified**: 2026-06-10 | **Last Amended**: 2026-06-10

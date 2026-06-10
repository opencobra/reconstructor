# Spec-Driven Development (Spec-Kit)

This repository uses [GitHub Spec-Kit](https://github.com/github/spec-kit/) for spec-driven
development. The Spec-Kit skills are installed under `.claude/skills/speckit-*` and the project
scaffold lives under `.specify/`.

## Constitution

The project's governing principles are in
[`.specify/memory/constitution.md`](../.specify/memory/constitution.md) (v1.0.0).

## Features

Each feature is specified under `specs/NNN-feature-name/` with `spec.md`, `plan.md`, and `tasks.md`.
An implementation report is added after `/speckit-implement` has built the feature.

- **001 — Multi-Reaction Tabbed Workspace**:
  [spec](../specs/001-multi-reaction-tabbed/spec.md) ·
  [plan](../specs/001-multi-reaction-tabbed/plan.md) ·
  [tasks](../specs/001-multi-reaction-tabbed/tasks.md)

## Workflow

`constitution → specify → (clarify) → plan → tasks → (analyze) → implement`, run via the
`/speckit-*` skills. Artifacts are committed so every implementation decision is auditable back to a
spec requirement.

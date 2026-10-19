---
name: start-rfc
description: Start a new Mantle Bloom RFC by choosing its number and setting up its proposal worktree. Use for requests to propose or write a new RFC; use continue-rfc for existing RFCs.
---

# Start work on an RFC

Use this skill only to start a new RFC. To resume, revise, review, accept, gate, or break
down an existing RFC, use [`continue-rfc`](../continue-rfc/SKILL.md), which owns existing
RFC state and next-step routing. An RFC is the parent decision for a model change; its
implementation issues start only after acceptance and any diagnosis gates pass.

## 1. Define the new RFC

- Read `AGENTS.md`, `rfc/README.md`, and relevant RFCs, issues and analysis to avoid
  proposing a duplicate or conflicting change.
- For a new proposal, choose the next `NNNN` by checking the index and existing filenames.
  Do not reuse a number or create a second RFC for a question an existing RFC already
  covers; extend or supersede the existing RFC when that fits.
- Before proposing a model change, read `docs/architecture.md` and
  `docs/simulation-model.md`. Read `docs/debugging.md` when investigating saved worlds.
- Use [`rfc-authoring`](../rfc-authoring/SKILL.md) to investigate, write, cite, and open the
  new draft RFC PR. The project owner accepts it by merging the PR, or by explicitly
  accepting an already-merged draft.

## 2. Set up the RFC worktree

- For a new RFC, create a proposal worktree from `origin/main`; the primary checkout's
  local `main` may lag behind `origin`. For worktree discovery and safe path resolution,
  follow [`continue-rfc`](../continue-rfc/SKILL.md#find-rfc-and-evidence-worktrees).

  ```bash
  git fetch origin
  git worktree add -b <agent>/rfc-<NNNN>-<slug> ../mantle-bloom-rfc-<NNNN> origin/main
  ```

  `<agent>` is `claude` or `codex`. Keep this worktree for the RFC proposal. Once accepted,
  each implementation sub-issue gets its own isolated worktree and branch from `origin/main`
  through [`start-issue`](../start-issue/SKILL.md); never switch the proposal worktree to an
  issue branch.
- Run `./bin/setup.sh` in a new worktree. Each worktree needs its own `backend/.venv`; do
  not symlink one from another tree.
- Pick free, non-user ports only if running the app. Ports 8000 and 5173 belong to the
  user's instance. If the RFC needs before/after replay numbers, follow `replay-ab` and use
  its detached `origin/main` baseline.

When the new proposal work begins, follow
[Starting or resuming an RFC](../project-management/SKILL.md#starting-or-resuming-an-rfc)
to update its project status. Existing RFC state, gates, implementation breakdown and
withdrawal are owned by [`continue-rfc`](../continue-rfc/SKILL.md).

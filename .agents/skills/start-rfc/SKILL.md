---
name: start-rfc
description: "Start or advance a Mantle Bloom RFC workstream: set up its worktree, confirm owner acceptance, run diagnosis gates, plan implementation issues, and coordinate implementation PRs. Use rfc-authoring to choose the RFC number and write its scientific content."
---

# Start or advance an RFC workstream

Use this skill when RFC work begins or advances, including starting a proposal, accepting
a draft, running its diagnosis gates, or planning accepted implementation work. It owns
workstream setup and lifecycle coordination, not the scientific prose. Use
[`rfc-authoring`](../rfc-authoring/SKILL.md) to choose a new RFC's number and write or
revise its content. If this is a reduced-context handoff, first use
[`continue-rfc`](../continue-rfc/SKILL.md) to reconstruct state, then return here for the
next lifecycle action. An RFC is the parent decision for a model change; implementation
implementation work does not begin until the RFC is accepted and any required diagnosis
gates pass.

## 1. Establish the RFC workstream

- Read `AGENTS.md`, `rfc/README.md`, the relevant RFC, tracking issue and PR. Verify current
  GitHub and repository state rather than relying on a status line or conversation summary.
- For a new proposal, ask [`rfc-authoring`](../rfc-authoring/SKILL.md) to choose the next
  RFC number. Do not choose or duplicate the number here. Before opening another proposal,
  check for an existing RFC covering the question.
- For project-board status changes, follow the relevant how-to in
  [`project-management`](../project-management/SKILL.md).

## 2. Set up or recover the RFC worktree

- Use `git worktree list --porcelain` and the path-resolution guidance in
  [`continue-rfc`](../continue-rfc/SKILL.md#find-rfc-and-evidence-worktrees) to discover
  an accessible RFC worktree and evidence. Never switch branches or overwrite another
  task's work.
- For a new RFC, create its proposal worktree from `origin/main` after `rfc-authoring` has
  selected the number. The primary checkout's local `main` may lag behind `origin`.
- For an existing draft whose worktree is unavailable, follow the isolated continuation
  branch procedure in [`continue-rfc`](../continue-rfc/SKILL.md#find-rfc-and-evidence-worktrees).
  For an accepted RFC, use an isolated worktree from `origin/main` when code or RFC files
  must be edited.

  ```bash
  git fetch origin
  git worktree add -b <agent>/rfc-<NNNN>-<slug> ../mantle-bloom-rfc-<NNNN> origin/main
  ```

  `<agent>` is `claude` or `codex`. Keep the proposal tree for RFC work. Each
  implementation sub-issue gets its own isolated worktree and branch from `origin/main`
  through [`start-issue`](../start-issue/SKILL.md); never switch the proposal worktree to
  an issue branch.
- Run `./bin/setup.sh` in a new worktree. Each worktree needs its own `backend/.venv`; do
  not symlink one from another tree.
- Use `rfc-authoring` to write the proposal and open its RFC proposal PR. Pick free,
  non-user ports only if running the app. For before/after replay numbers, follow
  `replay-ab`.

## 3. Advance RFC lifecycle

- **Draft:** the project owner accepts an RFC by merging its RFC PR. If it has already
  been merged while still marked Draft, only the project owner's explicit confirmation
  counts. Before gates or implementation, verify acceptance from GitHub or an explicit
  owner statement. If neither is clear, ask the project owner; a generic request to start
  work is not acceptance. Any agent may ask, but only the owner decides. Set the status to
  Accepted and create or update the tracking issue immediately, listing diagnosis gates
  and planned implementation steps.
- **Accepted with diagnosis gates:** run the specified gate, verify matching prior evidence
  before an expensive rerun, and record the result in the RFC. If a result calls the gate or
  acceptance criteria into question, propose a change and get the project owner's decision
  before applying a revised criterion.
- **Accepted: implementation planning:** create dependency-ordered sub-issues as
  children of the tracking issue using [`create-issue`](../create-issue/SKILL.md). Create
  issues at acceptance when their scope is independent of gate results; otherwise wait for
  the gate to determine the approach and scope. Do not begin implementation until required
  gates pass. Start each issue with [`start-issue`](../start-issue/SKILL.md), then open an
  implementation PR for that issue. An RFC with sub-issues may therefore have multiple
  implementation PRs, generally one per sub-issue. Update the tracking checklist as gates
  and implementation steps finish.
- **RFC content needs revision:** use [`rfc-authoring`](../rfc-authoring/SKILL.md). Use
  [`rfc-review`](../rfc-review/SKILL.md) for review work and
  [`complete-rfc`](../complete-rfc/SKILL.md) to implement or withdraw an RFC.

At each start or resumption, follow
[Starting or resuming an RFC](../project-management/SKILL.md#starting-or-resuming-an-rfc)
to update its project status. `rfc-authoring` opens the RFC proposal PR; this skill
coordinates and can open implementation PRs through the issue workstreams. Do not restate
board-status procedures here.

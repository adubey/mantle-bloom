---
name: start-rfc
description: Set up a local Mantle Bloom RFC worktree and coordinate its proposal with linked implementation issues. Use when starting or continuing an RFC workstream.
---

# Start work on an RFC

An RFC is the parent decision for a model change. It may be implemented through several
issues, but those issues start only after the RFC is accepted and any diagnosis gates have
passed. The RFC worktree is the shared local workspace for its proposal and linked issues;
each issue still gets its own branch and reviewable PR.

## 1. Identify the RFC and its state

- Read `rfc/README.md`, the relevant RFCs, and any linked issues or PRs. Check whether this
  is a new proposal, a draft under review, an accepted RFC awaiting diagnosis or issues, or
  an RFC whose implementation is already underway.
- **Confirm acceptance before working on an existing RFC.** If the user asks to start work on
  an RFC whose status is Draft, and hasn't said they accept it, ask whether they accept the
  RFC (as written, or with changes) before running its gates or opening issues. An RFC PR
  merged as a draft doesn't count as acceptance. If they accept, set its status to Accepted
  in the RFC header and the `rfc/README.md` index, with the date, and open its tracking
  issue (below). If not, treat it as a draft under review (section 3).
- **Open a tracking issue when an RFC is accepted.** Title it `RFC NNNN: <title>
  (tracking)`. Its body is a checklist, in dependency order, of every diagnosis gate and
  solution the RFC plans, each marked with the gate result it depends on, then the
  acceptance runs and setting the status to Implemented. Link the tracking issue from the
  RFC's **Related** line. Before opening it, check that the RFC doesn't already have one.
- For a new proposal, choose the next `NNNN` by checking the index and existing filenames.
  Do not reuse a number or create a second RFC for a question an existing RFC already
  covers; extend or supersede the existing RFC when that fits.
- Before proposing a model change, read `docs/architecture.md` and
  `docs/simulation-model.md`. Read `docs/debugging.md` when investigating saved worlds.

## 2. Set up the RFC worktree

- Worktrees are local to a machine and its agentic framework. Check whether the RFC worktree
  is accessible in the current environment; don't assume a Claude Code worktree is also
  available to Codex, or vice versa. If it isn't available here, set up a local RFC
  worktree from the canonical state that is accessible through Git: use the RFC PR branch
  to continue a draft, and `origin/main` for a new or accepted RFC. Unpushed work in another
  environment is not available through its worktree alone.
- For an RFC PR branch available on `origin`, fetch it and create a local worktree that
  tracks it:

  ```bash
  git fetch origin
  git worktree add --track ../mantle-bloom-rfc-<NNNN> origin/<rfc-branch>
  ```

- For a new RFC, or an accepted RFC whose content is already on `origin/main`, create a
  local RFC worktree from `origin/main`; the primary checkout's local `main` may lag behind
  `origin`.

  ```bash
  git fetch origin
  git worktree add -b <agent>/rfc-<NNNN>-<slug> ../mantle-bloom-rfc-<NNNN> origin/main
  ```

  `<agent>` is `claude` or `codex`. For an RFC already in review, continue its existing
  branch rather than creating a competing one. Once accepted, linked implementation issues
  use this RFC worktree directory but create their own branches from `origin/main` with
  `start-issue`; they do not go on the RFC branch.
- Run `./bin/setup.sh` in a new worktree. Each worktree needs its own `backend/.venv`; do
  not symlink one from another tree.
- Pick free, non-user ports only if running the app. Ports 8000 and 5173 belong to the
  user's instance. If the RFC needs before/after replay numbers, follow `replay-ab` and use
  its detached `origin/main` baseline.

## 3. Choose the next action from the RFC state

- **New or draft RFC:** use `rfc-authoring` for evidence, literature, design, citations,
  acceptance criteria, and the RFC PR. Keep the RFC in `rfc/NNNN-slug.md` and update
  `rfc/README.md`. Do not open implementation issues while the RFC is a draft.
- **Draft under review:** revise the same RFC and PR; use `rfc-review` when the task is to
  review rather than author it.
- **Accepted with diagnosis gates:** run the specified measurements, record the outcome in
  the RFC, and select the design the gate supports before opening implementation issues.
- **Accepted and ready to implement:** break the accepted design into the necessary
  dependency-ordered issues. Each issue should cite the RFC, name its own deliverable and
  acceptance checks, and avoid duplicating work assigned to another issue. Add each one to
  the tracking issue's checklist, with its PR once open. Start each
  issue using `start-issue`; it will reuse this RFC worktree and create an issue branch from
  `origin/main`.
- **Implementation underway:** use the tracking issue to identify remaining work, and tick
  items off as their PRs merge or their gates report. When a gate result rules a planned
  step out, strike it there with the reason. Update the RFC status to Implemented only
  after all of its implementation issues have landed, then close the tracking issue.

Keep issue creation and RFC acceptance aligned with the lifecycle in `rfc/README.md`: the
project owner accepts an RFC by merging its RFC PR, or by saying so for an RFC already merged
as a draft, and the RFC remains the record of the decision after implementation.

---
name: start-rfc
description: Set up a new Mantle Bloom RFC worktree and coordinate its proposal with linked implementation issues. Use when starting a new RFC workstream.
---

# Start work on an RFC

An RFC is the parent decision for a model change. It may be implemented through several
issues, but those issues start only after the RFC is accepted and any diagnosis gates have
passed. Keep proposal work in its RFC worktree; each implementation sub-issue gets its own
isolated issue worktree, branch and reviewable PR. For a cold handoff or an existing RFC
workstream, use [`continue-rfc`](../continue-rfc/SKILL.md) to reconstruct state first.

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

- Worktrees are local to a machine and its agentic framework. From this repository, run
  `git worktree list --porcelain` before creating another tree. Each block's `worktree`
  field is a checkout root; `branch` and `HEAD` identify its current ref and commit. Match
  these to the RFC PR or expected run commit, then resolve a repo-relative file beneath
  that root (for example `<worktree-root>/rfc/NNNN-slug.md`). Quote paths containing
  spaces. The listing doesn't identify the agent that created a tree, prove that its files
  are accessible, or show worktrees registered in another clone or host. Check access
  read-only. If an existing draft's tree isn't reachable, follow
  [`continue-rfc`](../continue-rfc/SKILL.md) for an isolated continuation branch; for a
  new RFC or an accepted RFC on `origin/main`, create a new local branch there. Never
  check out or push to another agent's shared RFC branch. Unpushed work in another
  environment is not available through its worktree alone.

- For a new RFC, or an accepted RFC whose content is already on `origin/main`, create a
  local RFC worktree from `origin/main`; the primary checkout's local `main` may lag behind
  `origin`.

  ```bash
  git fetch origin
  git worktree add -b <agent>/rfc-<NNNN>-<slug> ../mantle-bloom-rfc-<NNNN> origin/main
  ```

  `<agent>` is `claude` or `codex`. Once accepted, create every implementation sub-issue
  through the RFC tracking issue. Each sub-issue gets its own isolated worktree and branch
  from `origin/main` with `start-issue`; never switch the RFC proposal worktree to an issue
  branch.
- Run `./bin/setup.sh` in a new worktree. Each worktree needs its own `backend/.venv`; do
  not symlink one from another tree.
- Pick free, non-user ports only if running the app. Ports 8000 and 5173 belong to the
  user's instance. If the RFC needs before/after replay numbers, follow `replay-ab` and use
  its detached `origin/main` baseline.

## 3. Make run results survive a handoff

- Before running a diagnosis gate, find the required inputs and check for existing results
  in the RFC, linked issue/PR, and any referenced `analysis/` directory. On the same
  machine, `git worktree list --porcelain` can reveal another registered worktree, and its
  files may still be accessible even when another agent created it. Inspect relevant
  worktrees read-only; don't change branches or clean up another task's files.
- Write durable findings to `analysis/issue<N>/README.md` (for issue gates) or the RFC's
  diagnosis-gate section. Include the command/script, input save name and SHA, seed if
  applicable, step schedule and controls, commit for each compared run, the measured
  values, and repo-relative paths to checked-in scripts or compact JSON. Keep only compact
  JSON in the repository; full saves and replay outputs stay out of it. Never put a
  personal username, home directory, absolute machine path, or private mount path in an
  RFC, issue, PR, or checked-in analysis report. If another same-machine agent needs a
  local output's location, put it in an ignored, local-only handoff note and verify the
  note is ignored; identify the artifact in the durable report by a neutral basename or
  run label and its checksum where useful. See `replay-ab` for report fields.
- Treat temporary output paths as machine-local clues, not as the record of a result. If a
  result is needed for a later decision, summarize the evidence and reproduction details
  in the durable report before the current run context disappears. Don't copy large saves
  or raw replay output into the repo, and don't publish personal filesystem paths.

When RFC work starts or resumes, follow [Starting or resuming an RFC](../project-management/SKILL.md#starting-or-resuming-an-rfc)
to update its project status.

## 4. Choose the next action from the RFC state

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
  issue using `start-issue`, which creates an isolated issue worktree and branch from
  `origin/main`. Follow [Creating an issue](../project-management/SKILL.md#creating-an-issue)
  to set its initial project status.
- **Implementation underway:** use the tracking issue to identify remaining work, and tick
  items off as their PRs merge or their gates report. When a gate result rules a planned
  step out, strike it there with the reason. When all of its implementation issues have
  landed, close it out as Implemented with `complete-rfc`.
- **Falsified, superseded or no longer wanted:** report the evidence to the project owner.
  Once they decide to withdraw it, use `complete-rfc`.

Keep issue creation and RFC acceptance aligned with the lifecycle in `rfc/README.md`: the
project owner accepts an RFC by merging its RFC PR, or by saying so for an RFC already merged
as a draft, and the RFC remains the record of the decision after implementation.

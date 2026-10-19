---
name: continue-rfc
description: Resume, pick up, or continue work on an existing Mantle Bloom RFC or RFC-linked issue. Reconstruct its state, worktree, evidence and next action before changing anything.
---

# Continue an RFC workstream

Use this when a chat or agent is taking over an existing RFC, its diagnosis gates, or an
issue linked to it. Rebuild the project state from durable records before relying on a
handoff summary. For a new RFC, use [`start-rfc`](../start-rfc/SKILL.md); for a scientific
revision use [`rfc-authoring`](../rfc-authoring/SKILL.md); for a review use
[`rfc-review`](../rfc-review/SKILL.md).

## 1. Identify the decision and requested work

- Read `AGENTS.md`, `rfc/README.md`, and the named RFC in full, plus its linked PR/review
  discussion and referenced analysis reports. If no RFC number is given, use the issue/PR
  links and index to identify it; don't guess from a similar title. For an assigned issue,
  follow `start-issue` for its comments, scope and implementation readiness rather than
  duplicating those issue checks here.
- Record the RFC status, branch/PR, accepted decision, diagnosis-gate results, tracking
  issue, dependent RFCs/issues and the deliverable requested in this chat. Check these
  against current GitHub state when available; stale RFC headers and handoff notes may not
  reflect a merge or gate result.

## 2. Find RFC and evidence worktrees

- This is read-only discovery for RFC state and run evidence. For an assigned
  implementation issue, `start-issue` is the authority for choosing or creating its
  worktree and branch.
- From the current repository, inspect `git rev-parse --show-toplevel`,
  `git status --short --branch`, and `git worktree list --porcelain`. The porcelain output
  has one block per registered worktree: `worktree` gives its root path, `branch` gives
  the checked-out ref when attached, and `HEAD` gives its commit. Match the branch and
  commit to the RFC/PR or reported run; the output does not say which agent owns a tree.
  If using Codex, also inspect attached worktrees/artifacts when available. If using Claude
  Code, inspect its `.claude/worktrees/` directory when present. These may expose trees
  absent from the current checkout's listing; confirm each tree's repository, branch,
  commit and dirty state before treating it as evidence.
- A path in the output is the root of that checkout, not the location of a file inside it.
  Resolve repo-relative paths from that root: for example, for `analysis/issue<N>/README.md`,
  inspect `<worktree-root>/analysis/issue<N>/README.md`, or run
  `sed -n '1,200p' "<worktree-root>/rfc/NNNN-slug.md"` and
  `git -C "<worktree-root>" status --short --branch`. Quote the root path because it may
  contain spaces. If the path is inaccessible, the worktree belongs to another clone/host,
  or the target file isn't there, treat it as unavailable. If the task requires editing the
  RFC, use the isolated-branch procedure below; otherwise continue from durable records.
  Don't assume matching relative paths mean matching files.
- If using an attached/managed worktree, treat it as another clue: verify its root,
  branch, commit and dirty state against the porcelain listing before using it.
- Worktrees are local to a machine and framework. If the likely RFC worktree is listed
  and accessible, inspect its RFC, branch and status read-only. Never switch branches in
  a dirty worktree or one another task may be using. Don't reset, clean, stash, or delete
  another task's work.
- If the task edits an existing draft RFC and its worktree is unavailable, don't check out
  or push to the existing shared RFC branch. Fetch `origin`, create a new worktree on a
  unique branch based on the RFC PR branch, and submit the revision from that branch:

  ```bash
  git fetch origin
  git worktree add -b <agent>/rfc-<NNNN>-continue-<slug> \
    ../mantle-bloom-rfc-<NNNN>-continue-<slug> origin/<rfc-branch>
  ```

  Confirm the new branch name is unused locally and remotely first. Propose the new branch
  as a PR targeting the existing RFC branch so its current draft history is retained; do
  not force-push or push commits to the other branch. For accepted RFC edits, fetch
  `origin` and create a unique branch/worktree from `origin/main`; preserve the same
  isolation. Never base work on the primary checkout's possibly stale local `main`.
- A remote PR branch does not contain uncommitted or unpushed work from another machine.
  If the RFC or run depends on such work and it cannot be located in an accessible
  worktree or artifact, state the missing dependency and continue independent reading or
  analysis; don't silently recreate or claim that work.

## 3. Recover run information before repeating work

- Follow references from the RFC and issue into `analysis/`, `bin/debug/`, and known local
  artifact paths. Search by the RFC/issue number, seed, save name or distinctive run
  label. On the same machine, another agent's worktree may still be registered or its
  temporary outputs may remain accessible; use `git worktree list --porcelain` and inspect
  only relevant paths. Do not scan broad temporary directories or assume a path exists on
  another host.
- Prefer a checked-in `analysis/issue<N>/README.md` or RFC gate result. Verify it records
  the script/command, save name and SHA, seed, step schedule and controls, commit(s),
  measured values and conclusion. Check referenced compact JSON and scripts against the
  current files. A save or raw replay output may be local and unavailable to download;
  don't treat an unverified filename or a run on a different commit as equivalent evidence.
- Keep machine-specific paths private. Use the local-only handoff procedure in
  [`AGENTS.md`](../../../AGENTS.md#workflow); use a neutral artifact basename/run label and
  checksum in repository-facing records.
- If a usable raw output is found, identify its producing commit and inputs before using
  it. If required details are missing, reconstruct them from the command history or
  neighboring worktree/report where possible. Do not rerun a costly gate until checking
  existing reports and accessible outputs. When a run is necessary, write durable
  reproduction details and findings to the RFC or `analysis/issue<N>/README.md` before
  handing off; keep saves and full replay output out of Git.
- For scientific/code questions, read `docs/architecture.md`, `docs/simulation-model.md`,
  and the relevant topic page. Read `docs/debugging.md` for saved-world diagnostics.
  Trace the current implementation before relying on an RFC's provisional touchpoints.

## 4. Return to the RFC workflow

Once state, evidence and the requested deliverable are reconstructed, use
[`start-rfc`](../start-rfc/SKILL.md) to begin or advance the RFC lifecycle. Use
[`rfc-authoring`](../rfc-authoring/SKILL.md) for RFC content and its number when creating a
new RFC, [`rfc-review`](../rfc-review/SKILL.md) for review, and
[`start-issue`](../start-issue/SKILL.md) for an assigned implementation issue. This skill
owns handoff reconstruction and recovery; those skills own the work they route to.

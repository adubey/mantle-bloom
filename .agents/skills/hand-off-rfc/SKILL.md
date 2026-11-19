---
name: hand-off-rfc
description: Hand an unfinished Mantle-Bloom RFC workstream to the next agent. Push the work into a PR (draft or ready, with peer reviews if ready), put the saves where the next agent can find them (B2 if another machine), and leave one handoff issue with exactly the context needed. Use when asked to hand off, hand over, pass on, wrap up or stop work on an RFC or an RFC-linked issue that someone else will continue.
---

# Hand off an RFC workstream

Use this when you are stopping work on an RFC, or on an issue under one, and another agent
will continue it. The next agent starts with no access to this chat, so the PR, the saves
and one issue must carry everything. The agent receiving it uses
[`pick-up-rfc`](../pick-up-rfc/SKILL.md).

This skill coordinates the handoff and links to the skills that own each step. It is not for:

- evidence that challenges the RFC's design, where another session should investigate and
  revise (use [`edit-rfc`](../edit-rfc/SKILL.md));
- finishing or withdrawing the RFC (use [`complete-rfc`](../complete-rfc/SKILL.md));
- a handoff the user hasn't asked for. Stopping because your own context is long is a good
  moment to *propose* one, but wait for the user to agree.

**Approval.** When the user invokes this skill, that authorizes pushing your own agent
branch, opening a PR, requesting peer reviews, pushing files to B2 and creating or
commenting on the handoff issue. It does not authorize merging, force-pushing, or touching
another task's branch (see [`AGENTS.md`](../../../AGENTS.md#workflow)).

## 0. Settle the facts with the user

Ask only what you can't determine yourself, and ask it in one message:

- **Will the next agent run on another machine?** This decides where the files go (section
  2) and how reviews work (section 1). If the user doesn't say, assume the same machine and
  state that assumption in your final message so they can correct it.
- **Draft or ready?** Default to a draft unless the work is complete and passes
  [`pre-push-check`](../pre-push-check/SKILL.md). A ready PR must be merged by the owner
  before the next agent can start, so say that when you ask.
- **Does the user want to name the next agent or the issue?** Otherwise find it yourself.

**Running work.** A handoff should start after your runs have finished. Check for replays,
probes, test runs, servers and queued [`compute-budget`](../compute-budget/SKILL.md)
entries that are yours. If any are running, ask the user to choose, and don't pick for them:

- **(a) Wait** for them to finish, then hand off with their outputs included.
- **(b) Cancel** them and hand off now. Stop them cleanly (`./bin/stop.sh` for servers,
  kill your own processes only), and record in the issue which runs were cancelled, with
  their command and how far they got, so the next agent can decide whether to rerun them.

Identify the RFC, its tracking issue, the issue you were working on, and your branch and
worktree. If you can't tie the work to an RFC, say so and ask.

## 1. Put the work in a PR

Skip this step only if you have no files the next agent needs (a pure investigation whose
findings go entirely in the issue). Say that you skipped it and why.

1. **Collect the files.** `git status --short --branch` in your worktree. Stage by name,
   because the worktree may hold unrelated changes, and never use `git stash`. Include
   what the next agent needs and nothing else:
   - code and tests for the change in progress, even if partial;
   - `bin/debug/` scripts used for diagnostics, with their usage docstrings;
   - compact findings in `analysis/issue<N>/README.md` (or `analysis/exploration-<slug>/`),
     with the command, save name and SHA, seed, steps, commit, measured values and
     conclusion that [`continue-rfc`](../continue-rfc/SKILL.md#3-recover-run-information-before-repeating-work)
     expects;
   - an RFC edit, if the RFC itself changed.

   Leave out saves, full replay output, logs, personal paths and scratch files. If you find
   files you didn't create and can't explain, leave them alone and mention them in the issue.
2. **Run [`pre-push-check`](../pre-push-check/SKILL.md) before the first push**, then commit
   and push your own agent branch. For a ready PR it must pass. For a draft it may fail, but
   put its outcome in the PR so the next agent knows what is already clean; a draft can't
   hide a failure. Use the commit style in `AGENTS.md`. A commit that is knowingly incomplete
   says so in its body, for example "WIP: ... fails test_x". Branch from `origin/main` per the
   [`worktree` skill](../worktree/SKILL.md) if you don't have an issue branch yet, and use the
   issue-evidence pattern in [`edit-rfc`](../edit-rfc/SKILL.md#1-split-in-the-origin-session)
   when the evidence shouldn't wait for a PR.
3. **Open the PR** with the `AGENTS.md` format: a `## Summary` saying what changed and why,
   replay numbers if the model moved, and `Closes #<N>` only if this PR really completes
   that issue (otherwise `Part of #<N>`). Add a `## Handoff status` section that states
   what is done, what is not, and what is failing or unchecked. Open it as a **draft** unless
   the user said ready.
4. **If it is ready (not a draft)**, the `backend` and `frontend` checks should be green or
   running, and it needs peer reviews before you ask the owner to merge:
   - Request **one required** reviewer and optionally a second, different agent, with
     [`peer-review`](../peer-review/SKILL.md) (`bin/peer_review.py request`, one `--reviewer`
     each). `request` can't mark a reviewer optional, so you keep the distinction yourself:
     the first is required, the second is advisory. Ask the user which agents or models they
     want if they haven't said; otherwise pick two other than yourself.
   - **You** wait for the required review, address it per `peer-review` section 4, and only
     then ask the owner to merge. `bin/peer_review.py wait --pr <N>` returns only when every
     requested slot is done, so run it in the background but don't let it gate the handoff:
     once the required reviewer's comment is on the PR, proceed. If the advisory review is
     still running, say so in the issue and the pick-up agent reads it afterward.
   - Peer-review state is **per machine** (`mantle-bloom-handoff/reviews/`), which is why the
     outgoing agent owns the required review whichever machine the next agent is on. If the
     user can't wait for it, convert the PR to a draft (`gh pr ready --undo <N>`) and say that
     reviews were not done.

Never merge. Merging is the owner's job. **The next agent starts from `origin/main`, so a
ready PR must be merged before it picks up.** Say so in the handoff (section 3) and tell the
user that the PR needs merging, and that the pick-up waits for it. A draft PR isn't merged:
the next agent takes the changes from it into its own branch and PR, then closes the draft
(see `pick-up-rfc`). Either way, the next agent never pushes to your branch.

## 2. Make the saves findable

List every save or output the next agent needs and that isn't in Git: the input saves it
should start from, and any outputs of runs you made that it shouldn't redo. Leave out
anything cheap to regenerate or no longer relevant, and say in the issue if you deliberately
dropped something.

1. **Put them in one place**, `~/mantle-bloom-captures/<rfc-or-issue>/`, per
   [`AGENTS.md`](../../../AGENTS.md#workflow), never `/tmp` (volatile) or the repo. Move or
   copy them there if they are elsewhere, and write a `SHA256SUMS` file beside them. Run it
   from inside that directory, so it records bare file names:
   `(cd ~/mantle-bloom-captures/<rfc-or-issue> && shasum -a 256 <files> > SHA256SUMS)`.
2. **Same machine:** that directory is enough. Put the directory name and any private note
   in `$handoff` (below); the issue gets file names and checksums only.
3. **Another machine:** upload with the [`b2-saves`](../b2-saves/SKILL.md) skill:

   ```bash
   cd ~/mantle-bloom-captures/<rfc-or-issue>
   cp SHA256SUMS SHA256SUMS.<rfc-or-issue>     # the bucket keeps one object per name
   bin/b2_saves.py push <saves...> SHA256SUMS.<rfc-or-issue> --note "<rfc-or-issue> handoff"
   bin/b2_saves.py push "$handoff/rfc-<NNNN>-handoff-<issue>.md" --note "<rfc-or-issue> handoff note"
   bin/b2_saves.py list      # confirm every file appears, with the expected size
   ```

   The note lives outside the captures directory, so it needs its own `push` (write it first,
   in step 4). `push` may rename a file whose name collides, or reuse an existing object with
   the same content, so read each file's **actual bucket name** from `list` and put it in the
   issue's evidence. Don't delete the local copies until the user says the next agent has
   them.
4. **Private note.** Put local-only details in
   `$handoff/rfc-<NNNN>-handoff-<issue>.md`, resolving `$handoff` as `AGENTS.md` describes:
   worktree path and branch, captures directory, ports, any log or diagnostic output that
   isn't worth uploading, and exact commands to resume. Never copy its contents or any
   personal path into the issue or PR.
5. **Confirm nothing is running** (section 0 settled this). Leave your worktree in place:
   the next agent doesn't use it (see `pick-up-rfc`), but it may hold evidence.

## 3. Find or create the handoff issue

1. **Find the RFC tracking issue** from `rfc/README.md` or the RFC's header. An RFC still in
   draft has none (the tracker is created at acceptance, per the
   [RFC lifecycle](../../../rfc/README.md#lifecycle)). Then the **proposal PR** is the home of
   the handoff: post the `## Handoff` as a comment there, create no sub-issue, and
   `pick-up-rfc` looks for it on that PR. Otherwise list the tracker's sub-issues:

   ```bash
   gh api repos/:owner/:repo/issues/<tracking>/sub_issues --jq '.[] | "\(.number) \(.state) \(.title)"'
   gh issue view <tracking>      # the checklist shows the dependency order
   ```

2. **Reuse an open sub-issue** when its deliverable is what the next agent should do, which
   is usually the issue you were working on if it isn't done, or the next one in the
   checklist. Post the handoff as a **comment** headed `## Handoff`; don't edit someone
   else's issue body.
3. **Otherwise create a new sub-issue** with [`create-issue`](../create-issue/SKILL.md),
   with `--parent <tracking>`, add it to the tracking checklist in dependency order, and put
   it on the project board per [`project-management`](../project-management/SKILL.md#creating-an-issue)
   as `Ready`, or `On Hold` if it waits on a ready PR's merge (the prerequisite rule in `project-management`). Whoever finds the PR merged moves it to `Ready`, or `In Progress` if starting. Use the handoff text as the body, starting with `## Handoff`. Follow the RFC
   lifecycle's rules about which issues are eligible; if the RFC isn't accepted or a gate
   hasn't passed, the handoff is for that gate or revision step, not for implementation.
4. If your own issue is unfinished, leave it open and say in its comments that the work
   moved to the handoff issue (or the same issue, if you reused it). Set board statuses per
   `project-management` (`On Hold` if it now depends on the handoff issue).

### What goes in the handoff

The rule: **include a line only if the next agent would act differently without it, and
couldn't get it from the RFC, the issue history, the PR diff or the code.** Link instead of
copying. Don't narrate the session, restate the RFC, paste diffs or logs, list every
command you ran, or include opinions you didn't act on. Omit a section when it has nothing
in it. Use neutral run labels and no personal paths.

```markdown
## Handoff

**Next action:** <one deliverable, and the check that shows it is done>

**State**
- RFC NNNN (<status>), tracking #<T>; this issue: #<N>
- PR #<P>, <draft: take its changes into your own branch, then close it | ready: pick up
  after it merges>, branch `<agent>/issue-<N>-<slug>` at `<short sha>`
- Reviews: <not requested | requested from X, Y on <date> | complete, link>
- Checks: <pre-push-check result; what fails or was not run, and why>

**Evidence** (for each run the next agent needs)
- `<run label>`: command `<...>`, save `<name>` sha256 `<...>`, seed `<...>`, steps `<...>`,
  commit `<sha>`, result `<number>`. Where: <B2 as `<bucket name>` (`bin/b2_saves.py pull <bucket name>`; local name `<name>`) | this
  machine, see private note `<file name>`>

**Don't repeat / ruled out:** <what failed or was null, and the number that showed it>

**Open decisions:** <what needs the owner, and what you recommend>

**Held out:** <seeds or saves that must not be used for tuning, if any>
```

If a section's content is already in the PR's `## Handoff status`, link it instead.

## 4. Check, then hand over

Before you say it is done, verify and report each item, with any failure stated plainly:

- the branch is pushed and `git status` shows nothing of yours uncommitted;
- the PR exists in the right state, and any review requests were made;
- each file named in the issue exists where the issue says, and its checksum matches;
- the issue links the PR, and the tracking checklist and board are updated;
- nothing of yours is still running, or the user chose to cancel it and the issue says what
  was cancelled.

End with a prompt for the user to paste into the next session, and a link to the issue and PR:

> Use `pick-up-rfc` for RFC NNNN. The handoff issue is #N and the private note is
> `rfc-NNNN-handoff-N.md` (on this machine | in B2).

If the PR is ready, say that it needs merging before the next agent starts. Also say whether
you assumed same-machine, draft or ready, and anything you left out on
purpose.

---
name: edit-rfc
description: "Split work into a new session when evidence challenges an existing Mantle Bloom RFC: write a durable brief, explore and revise, and return the outcome to affected issues and waiting sessions. Not for a revision already agreed within the current session (use rfc-authoring), exploration with no RFC at stake (use explore-problem), or plain resumption (use continue-rfc)."
---

# Edit an RFC from new evidence

Use this when work on an RFC, usually in a `continue-rfc` or `start-issue` session, produces
evidence that may change the RFC. The session where the problem appears is the *origin
session*. The new session that investigates it and revises the RFC is the *edit session*.
Splitting them keeps each session's context short and leaves a record that a new session can
work from.

This skill coordinates the split and the return, and links to the skills that own each step:
[`explore-problem`](../explore-problem/SKILL.md) for the investigation,
[`rfc-authoring`](../rfc-authoring/SKILL.md) for the RFC content,
[`complete-rfc`](../complete-rfc/SKILL.md) for a withdrawal, and
[`project-management`](../project-management/SKILL.md) for board status.

Follow [Communicate clearly](../explore-problem/SKILL.md#communicate-clearly) in every
message, brief, comment and note. A later session has only what these records say.

**Find your phase.** Use [section 1](#1-split-in-the-origin-session) if you are in the
session that found the problem. Use [section 2](#2-explore-and-edit-in-the-edit-session) if
you were given a tracker or brief. Use [section 3](#3-return-after-the-owners-decision) if
the owner has decided the outcome. Use [section 4](#4-resume-an-affected-issue) if an issue
you are working on has a return note.

**Approval.** When the project owner invokes this skill or asks for this split, that approves
an exploration tracker and a revision PR if the evidence leads there. It doesn't choose that
outcome. An agent that thinks a split is warranted proposes it and waits for the owner's
agreement. Withdrawing an RFC, starting a new RFC, closing issues and implementing a change
each need the owner's own decision.

## 1. Split: in the origin session

**Decide whether to split.** Propose a split when:

- a diagnosis challenges the RFC's premise, design, a gate or an acceptance criterion;
- a failed gate can't be explained, or a falsification criterion may have been met;
- answering the question needs new probes or literature outside the issue's scope.

Context size is hard to judge from inside a session, so use these events rather than waiting
for it to feel long, and split while this session can still write a complete brief. Don't
split for a clarification or implementation detail this session can settle in place.

**Check for an open revision.** Look on the RFC tracking issue, or the PR of a draft RFC, for
an open `## RFC revision brief`. If one covers this problem, add your evidence to its tracker
and your session to its private note instead of writing a second brief. If a different
revision of the same RFC is open, wait for it or build on it, and tell the user.

**Make the evidence durable.**

- Commit the evidence, staging the files by name, because the worktree may hold unrelated
  changes. Put debug scripts in their own commit, separate from any other work in progress, so
  the edit session can cherry-pick them cleanly.
- Push that commit to a new branch that is never rewritten, without switching this worktree
  off its branch. `<N>` is the origin issue:

  ```bash
  git push origin HEAD:refs/heads/<agent>/issue-<N>-evidence-<slug>
  ```

  A later rebase of the working branch then can't lose the commit. This branch never becomes a
  PR, so it doesn't need `pre-push-check`. Leave it in place after the RFC closes. Never use
  `git stash`.
- Put saves and large outputs in `~/mantle-bloom-captures/<rfc-or-issue>/`, as
  [`AGENTS.md`](../../../AGENTS.md#workflow) describes, not in a temporary directory or the
  handoff directory. If another machine needs them, use the
  [`b2-saves`](../b2-saves/SKILL.md) skill. Compare checksums after any copy.

**Create the tracker.** Use [`create-issue`](../create-issue/SKILL.md) to create one new issue
for this revision, linking the RFC and origin issue, and add the
[`explore-problem`](../explore-problem/SKILL.md#1-frame-the-problem-with-the-user) summary
section. The edit session uses this issue even if little exploration is needed.

**Write the brief.** It has a public and a private part, and together they must let a new
session start with no access to this chat.

Post the public part on the tracker as a comment headed `## RFC revision brief`, using
neutral run labels and no personal paths:

- the RFC, its branch or PR, the origin issue, and the design or criterion in question;
- what was observed, with numbers, and why it challenges the RFC;
- the current hypothesis and the alternatives;
- what has been tried or ruled out, including null results, and work not to repeat;
- the question the edit session should answer;
- the evidence branch and commit, and for each run the details
  [`continue-rfc`](../continue-rfc/SKILL.md#3-recover-run-information-before-repeating-work)
  expects;
- seeds that are held out or otherwise must not be used;
- the origin work's state, and how the answer could affect it;
- the gates, issues and PRs known so far that the answer could affect, as the first rows of
  the impact table;
- the private note's file name.

Link the brief from the origin issue and the tracking issue, so sessions working on sibling
issues can find it. For a draft RFC with no tracking issue, link it from the RFC PR.

Put the private part in `$handoff/rfc-<NNNN>-edit-<slug>.md`, resolving `$handoff` as
described in [`AGENTS.md`](../../../AGENTS.md#workflow):

- worktree paths, branches and local ports;
- locations of saves, logs and diagnostic output;
- exact commands to rerun the diagnosis;
- each waiting session: its issue, branch, worktree, and the prompt the user should give it
  to resume.

If the edit session may run on another machine, use the [`b2-saves`](../b2-saves/SKILL.md)
skill to send this note with the evidence.

**Hand over.** Give the user a prompt to paste into a new session:

> Use `edit-rfc` section 2 for RFC NNNN. The tracker is #T and the private note is
> `rfc-NNNN-edit-<slug>.md`.

**Choose what the origin session does next.** Pause only if the revision affects its work, and
ask the user if you can't tell. If it pauses, follow
[Investigating additional work](../project-management/SKILL.md#investigating-additional-work)
to put the issue `On Hold` with the tracker as its dependency. If the work is independent,
continue it. Both sessions follow [`compute-budget`](../compute-budget/SKILL.md) for heavy runs.

## 2. Explore and edit: in the edit session

**Check the brief first.** Read both parts, then confirm that each named commit, path and
checksum exists and matches. Treat anything missing as a missing dependency, as
[`explore-problem`](../explore-problem/SKILL.md#2-investigate-interactively) describes, and
don't silently recreate it. If you were asked to edit an RFC and no brief exists, create the
tracker and write the public part yourself from what the user gave you.

**Explore.** Follow [`explore-problem`](../explore-problem/SKILL.md), including preregistered
probes and the `## Exploration summary`. Run probes in a scratch worktree. To reproduce a
result that depends on unmerged code, use a detached worktree at the brief's evidence commit.
Treat the origin worktree as read-only. If the evidence is already sufficient, say so. If the
evidence turns out to affect a session that kept working, tell the user straight away.

**Keep this session resumable.** Keep the tracker's summary current and write long outputs to
files. If this session grows long, typically between exploring and revising, post a new brief
on the tracker, copy the waiting list to a new private note, and say that the new session
owns the rest of this work.

**Get the owner's choice.** Present the outcome as
[`explore-problem`](../explore-problem/SKILL.md#3-present-outcomes-and-ask-what-to-create)
describes, and wait for the owner's choice. Outcomes can be combined.

| Outcome | What to do | Return starts |
|---|---|---|
| Revise the RFC | Continue below | When the revision PR merges |
| RFC stands, or the result is inconclusive | If gate results, findings or scripts need keeping, open a records PR | When the owner decides, after any records PR merges |
| Small follow-up issues | Use [`create-issue`](../create-issue/SKILL.md) | When the owner decides |
| Withdraw the RFC, which you should recommend if a falsification criterion was met | Use [`complete-rfc`](../complete-rfc/SKILL.md), which settles the RFC, its issues and the index | When the completion PR merges |
| New RFC | Use [`rfc-authoring`](../rfc-authoring/SKILL.md) and [`start-rfc`](../start-rfc/SKILL.md) | When the owner chooses this outcome; next step is to wait for the new RFC, and the existing RFC stays in force until the new one is accepted |

A failed gate still blocks implementation under the
[RFC lifecycle](../../../rfc/README.md#lifecycle), even if the RFC stands. Its return note must
ask for a rerun or the owner's decision, not for work to continue.

**Revise.** Use the tracker's branch through
[`start-issue`](../start-issue/SKILL.md#start-work-on-an-issue) steps 2 and 3. This is a
documentation branch, so it doesn't mean implementation has begun. Follow
[`rfc-authoring`](../rfc-authoring/SKILL.md) for the RFC. Cherry-pick the debug scripts the
revision cites rather than re-creating them, so that Git drops its identical copy without a
conflict when the origin branch later rebases. Don't cherry-pick a commit that mixes in work
in progress.

**Assess the impact.** Build this table for the PR body. List every diagnosis gate, sub-issue
and open PR of the RFC, merged work that relied on the changed design or criterion, and
dependent RFCs:

| Item | Verdict | Reason |
|---|---|---|
| Gate, issue, PR or RFC | Unchanged, invalidated, must rerun, superseded or blocked | What in the revision causes this |

- **Unchanged:** the item stands as it is.
- **Invalidated:** its results or plan no longer hold and it needs re-scoping.
- **Must rerun:** its plan stands but its result must be measured again.
- **Superseded:** the revision replaces it, so no further work is needed.
- **Blocked:** it can't proceed until something named in the reason happens.

Below the table, state which seeds remain usable as held-out evidence. A seed or save used to
form or tune the new hypothesis can't also validate it, so a revision based on held-out runs
needs new held-out seeds.

**Propose criteria changes in the PR body.** For each changed gate or acceptance criterion,
give the old criterion, the new criterion and the reason. The owner's merge of a PR that sets
these out is their decision on those changes, and only those changes. Don't use a revised
criterion for any rerun until the PR merges, unless the owner approves it earlier.

**Open the PR and stop.** Follow [`pre-push-check`](../pre-push-check/SKILL.md) first. Use a
`## Summary` that explains the evidence and what changed, then the two tables. End with
`Refs #<tracker>` rather than `Closes`, because the tracker stays open through Return. Put
the index row, with the status line's wording, plus ordering prose and **Related** lines in
this PR. Add a pending-revision comment
with the PR link to the tracking issue and each item not marked Unchanged, so nobody starts
stale work. Tell the user the PR is ready and give them this prompt for after the merge, which
works in any session:

> Use `edit-rfc` section 3 for PR #X.

**For a draft RFC,** the revision PR targets the RFC branch, using the isolated branch
procedure in [`continue-rfc`](../continue-rfc/SKILL.md#2-find-rfc-and-evidence-worktrees) even
if the draft's worktree exists. Revise in place, with no revision note. Merging that PR doesn't
accept the RFC or change its status or index row. In the return note, waiting work either
rebases onto the RFC branch, which the note names as its base, or waits for acceptance.

## 3. Return: after the owner's decision

Any session can run this, using the tracker, the merged PR and the private note. Verify the
decision and any merge on GitHub first. If the owner closed a PR unmerged, ask them which
outcome now applies. If the outcome has no impact table, build one now and post it on the
tracking issue, since waiting issues still need verdicts. For a draft RFC with no tracking
issue, post it on the RFC PR. Then recheck the table against any work that started or landed
since.

**Settle the tracking issue.** Update its checklist and sub-issues. Merging approves criteria
changes only, so confirm each closure with the owner, unless they said the table's verdicts
stand. Settle each issue by what happened to it:

- Work that landed and met its criteria closes through
  [Completing an issue](../project-management/SKILL.md#completing-an-issue).
- Superseded or no longer needed work closes through
  [Choosing not to pursue an issue](../project-management/SKILL.md#choosing-not-to-pursue-an-issue).
- Blocked work stays open and `On Hold`.

Create new sub-issues only when the [RFC lifecycle](../../../rfc/README.md#lifecycle) allows
them. For an open PR, comment with the outcome and a recommendation, and don't close it or
push to it without the owner's go-ahead. After a withdrawal, `complete-rfc` has already done
this step, so post the return notes only.

Any later change to checked-in files needs a new branch from `origin/main`, because merged
branches are deleted. Don't hold return notes for it.

**Post a return note** on the tracking issue, every issue in the private note's waiting list,
and every item whose verdict isn't Unchanged. Link the brief, and use this form:

```markdown
## RFC revision return

- **Brief:** <link>
- **Outcome:** <PR and merged commit, or why no edit was made>
- **What changed:** <what this means for this issue, in plain language>
- **Impact table:** <link>
- **Next step for this issue:** <continue | rebase onto a named base and continue | rerun a
  named gate | re-scope | stop | wait for a named item>
- **Evidence:** <which earlier results and seeds stand, and which don't>
```

Use the next step that matches the verdict, then set the status through
[`project-management`](../project-management/SKILL.md#3-lifecycle-how-to):

| Verdict | Next step | Status while waiting to resume |
|---|---|---|
| Unchanged | Continue, or rebase onto the named base and continue | `Ready` |
| Must rerun | Rerun the named gate | `Ready` |
| Invalidated | Re-scope | `Ready` |
| Superseded | Stop | `Deprecated`, once the owner confirms closing it |
| Blocked | Wait for the named item | `On Hold` |

**Finish.** Close the tracker through
[Completing an issue](../project-management/SKILL.md#completing-an-issue) once the user confirms
its summary is complete. Tell the user every session that needs to act, including any still
running, with a one-line prompt that links its return note.

## 4. Resume an affected issue

Read the return note first, and ignore one this issue has already acted on.

- **Continue, or rebase and continue:** first confirm the branch belongs to this task, and use
  the **New session** case below if it doesn't. Run `git fetch origin`, then rebase the branch
  onto the base the note names. If the branch was already pushed, push it with
  `--force-with-lease`. Rebasing this task's own branch doesn't need the user's approval.
- **Rerun, re-scope, stop or wait:** follow the note. Results it marks invalid don't count as
  evidence for the revised RFC.
- **New session:** get the branch and worktree from the private note, and confirm with the
  owner that this session takes over that branch. Otherwise branch from its head under a new
  name. If no work exists yet, continue through
  [`start-issue`](../start-issue/SKILL.md#start-work-on-an-issue) step 2.

Then set the issue to `In Progress` through
[Starting an issue](../project-management/SKILL.md#starting-an-issue).

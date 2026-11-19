---
name: pick-up-rfc
description: Pick up a Mantle-Bloom RFC workstream that another agent handed off with hand-off-rfc. Read the handoff issue, check that its claims are true, fetch the saves, set up your own worktree, and start. Use when given a handoff issue, a "## Handoff" comment or a handoff prompt.
---

# Pick up an RFC handoff

Use this when you are given a handoff issue, or a prompt such as "use `pick-up-rfc` for RFC
NNNN, the handoff issue is #N". The previous agent wrote it with
[`hand-off-rfc`](../hand-off-rfc/SKILL.md). The handoff is a claim about the state of the
work at one moment, so check it before relying on it. This skill is only the pick-up. For
reconstructing an RFC with no handoff, use [`continue-rfc`](../continue-rfc/SKILL.md).

## 1. Read the handoff

- Read `AGENTS.md`, then the handoff issue in full with its comments:
  `gh issue view <N> --comments`. The handoff is the `## Handoff` body or comment (take the
  latest if there are several; the older ones are history). If the user gave no issue
  number, look through the RFC tracking issue's sub-issues
  (`gh api repos/:owner/:repo/issues/<tracking>/sub_issues`) for a `## Handoff` that nobody
  has picked up, and confirm it with the user rather than guessing.
- Read the RFC's status line and gate record (see [`start-issue`](../start-issue/SKILL.md)),
  and the handoff PR's description and `## Handoff status`. Read the rest of the RFC only
  where the next action touches it. If the RFC isn't accepted or a gate hasn't passed, the
  next action must be that step; if the handoff asks for more, raise it with the user.
- Check for an open `## RFC revision brief` and for `## RFC revision return` comments, as
  `continue-rfc` and `start-issue` describe.

## 2. Verify the claims

Check each of these against live state, and note every mismatch. Don't assume a stale
handoff is correct, and don't silently recreate anything that is missing:

- **Branch and PR:** `git fetch origin`, then the branch exists at the SHA the handoff names
  (or a descendant), and the PR's draft or ready state matches what it says.
  `gh pr view <P>` shows its state. **A ready PR must be merged before you start.** If it
  isn't, tell the user it is waiting on them to merge it, and stop; don't work from its
  branch instead. A draft PR is the exception, handled in section 4.
- **Files:** each save, capture and script the handoff names exists, and its SHA-256 matches
  the handoff. Get them (section 3) before checking.
- **Checks:** the failing and unchecked items it lists are still what they say. Don't rerun
  a whole suite to find out. Run only what the next action needs (`AGENTS.md` and
  [`compute-budget`](../compute-budget/SKILL.md) apply).

If something is missing or false, tell the user which claim failed and what it blocks, and
continue with independent work. If the whole handoff is untrustworthy (wrong branch, files
gone), stop and ask.

## 3. Get the files

- **Same machine:** the saves are in `~/mantle-bloom-captures/<rfc-or-issue>/` with a
  `SHA256SUMS` file, and the private note is
  `"$(cd "$(git rev-parse --git-common-dir)" && pwd)/mantle-bloom-handoff/<file named in the issue>"`.
  Read it for local paths and commands. Never copy its contents or any personal path into
  Git, an issue or a PR.
- **Another machine** (the issue says B2, or the paths don't exist here): use the
  [`b2-saves`](../b2-saves/SKILL.md) skill.

  ```bash
  bin/b2_saves.py list
  MANTLE_BLOOM_SAVES_DIR=~/mantle-bloom-captures/<rfc-or-issue> bin/b2_saves.py pull <name>
  (cd ~/mantle-bloom-captures/<rfc-or-issue> && shasum -a 256 -c SHA256SUMS)
  ```

  Pull `SHA256SUMS` and the private note by name the same way. If a file is not in the
  bucket, say so; the previous agent may still have it locally.
- Don't reproduce a costly run when a verified output exists. Check
  [`continue-rfc`](../continue-rfc/SKILL.md#3-recover-run-information-before-repeating-work)
  for how to judge whether evidence matches the commit and inputs.

## 4. Set up your own tree

- Create a **new** worktree and branch from `origin/main` with the
  [`worktree` skill](../worktree/SKILL.md), named with your own agent identifier, then run
  `./bin/setup.sh`. Never push to the previous agent's branch or reuse its worktree, even
  on the same machine.
- **If the handoff PR was ready and has merged,** its changes are already in `origin/main`.
  Nothing more to bring in.
- **If the handoff PR is a draft,** patch its changes into your branch, for example
  `git cherry-pick <first>^..<last>` of its commits (or `git merge --no-ff` it into your
  branch if the history is clean), and resolve conflicts against current `main`. Check
  that what arrived is what the handoff says it contains. Open **your own PR** with those
  changes and your own work, with a `## Summary` that says which changes came from the draft
  (`Part of`/`Closes #<N>` as appropriate). Then close the draft with a comment linking your
  PR (`gh pr close <P> --comment "Superseded by #<yours>"`). Close it only after your PR is
  open. Don't mark the draft ready, and don't delete its branch.
- Update the board per [`project-management`](../project-management/SKILL.md#starting-an-issue):
  the handoff issue moves to `In Progress`.
- Comment on the handoff issue that you've picked it up, with your agent, branch and anything
  that didn't match from section 2. This is the signal to everyone else that it is taken.

## 5. Handle open reviews

If the handoff says peer reviews were requested on the (now merged) PR, read what they
posted (`gh pr view <P> --comments`; `gh api repos/:owner/:repo/pulls/<P>/reviews`) and
check whether the findings were addressed before the merge. Summarize any that weren't for
the user, and fix the straightforward ones in your own PR with a `## Review fixes` section.
A reviewer slot that the previous agent opened is not yours to claim or complete.

## 6. Work, then close the loop

Continue with the handoff's next action, using [`start-issue`](../start-issue/SKILL.md),
[`start-rfc`](../start-rfc/SKILL.md) and the other skills that own that work. Follow the
`## Open decisions` and `## Held out` items as rules.

When you finish or must stop, leave the state in the issue: either complete it per
[Completing an issue](../project-management/SKILL.md#completing-an-issue), or write a new
handoff with [`hand-off-rfc`](../hand-off-rfc/SKILL.md). If the handoff was wrong or
incomplete in a way that cost you time, add a short correction comment to the issue so the
next reader isn't misled.

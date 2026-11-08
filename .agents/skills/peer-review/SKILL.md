---
name: peer-review
description: Coordinate independent multi-agent code reviews on a pull request across Antigravity, Claude Code, Codex, and Muse Code, preferring already-open windows with a headless CLI option, and notify the author agent when all reviews complete.
---

# Multi-agent peer review

Use this skill when the user asks an author agent to request independent code reviews from
other agents (`antigravity`, `claude`, `codex`, `muse`) with chosen models and reasoning
levels, or when an agent is asked to stand by for or perform one of those reviews. For an
RFC pull request, reviewers also follow [`rfc-review`](../rfc-review/SKILL.md).

State lives in `$(git rev-parse --git-common-dir)/mantle-bloom-handoff/reviews/pr-<N>.json`,
which all worktrees on the machine share and which never enters Git.

## 1. Request reviews (author agent)

1. Ensure the branch is pushed and the PR is open or updated on GitHub (`<N>` is the PR
   number).
2. Create the review request with `bin/peer_review.py request`, passing one `--reviewer`
   flag per requested reviewer in the form `AGENT[:MODEL[:REASONING]]`. Open-window mode is
   the default; pass `--spawn` only when the user asks to launch headless CLI processes
   instead of using open windows:

   ```bash
   python3 bin/peer_review.py request --pr <N> \
     --author <your-agent> \
     --reviewer <agent1>:<model1>:<reasoning1> \
     --reviewer <agent2>:<model2>:<reasoning2>
   ```

   `request` waits briefly (`--claim-wait 2.0`) to see whether standby listeners in
   already-open windows claim the slots automatically.
   - For any slot claimed by a standby window, tell the user that window has picked up the
     review.
   - For any slot still pending, give the user the printed one-line prompt to paste into
     that open window after choosing the requested model and reasoning level in its UI.

3. Start the completion watcher **in the background** so this author session stays
   available while the reviewers work, and wakes automatically when all reviews are done:

   ```bash
   python3 bin/peer_review.py wait --pr <N>
   ```

   Use your harness's background command execution (for example, `WaitMsBeforeAsync` in
   Antigravity or `run_in_background` in Claude Code) rather than blocking the foreground
   turn, unless the user asks to block.

## 2. Stand by or claim a review (reviewer agent)

### Option A: Standby listener in an open window (zero copy-paste)

When the user asks an open window to stand by as a reviewer, run:

```bash
python3 bin/peer_review.py listen --agent <your-agent> \
  [--model <active-model>] [--reasoning <active-reasoning>]
```

Run it in the background (or foreground if the window is dedicated to reviewing). As soon as
an author agent requests a review for `<your-agent>`, `listen` atomically claims the slot,
prints the review brief, and exits `0`. After completing the review (section 3), start
`listen` again if the user wants the window to remain on standby.

### Option B: Ad-hoc prompt in an open window

When the user selects a model and reasoning level in an open window and pastes a prompt such
as `Use the peer-review skill to review PR #<N> as slot '<slot>'`, claim the slot and read
its brief:

```bash
python3 bin/peer_review.py claim --pr <N> --slot <slot-or-agent> \
  [--model <active-model>] [--reasoning <active-reasoning>]
```

## 3. Perform an independent review and post to GitHub (reviewer agent)

1. **Stay independent.** Do not read or rely on other peer review comments posted on PR
   `#<N>` after the request timestamp. Form your own assessment from the code and diff:
   - Inspect `gh pr view <N>` (description, linked issue, author summary) and
     `gh pr diff <N>`, or create a temporary detached worktree at the PR's head SHA if you
     need to inspect surrounding code or run a targeted check:

     ```bash
     git fetch origin pull/<N>/head
     git worktree add --detach ../mantle-bloom-review-<N>-<slot> FETCH_HEAD
     ```

   - Remove any temporary review worktree (`git worktree remove ../mantle-bloom-review-<N>-<slot>`)
     when finished.
2. **Focus on high-signal checks.** As [`AGENTS.md`](../../../AGENTS.md#workflow) specifies,
   when the PR's author reports its tests passing, do not rerun the test suite. Check:
   - whether the change respects [`docs/design-principles.md`](../../../docs/design-principles.md),
     conservation ledgers, physical units, and timestep invariance;
   - save compatibility (`.mbworld` pickle migration) when state fields change;
   - complexity, dead/duplicate code ratchets, and edge cases;
   - targeted one-off checks the author's unit tests may not have covered.
   For an RFC PR, follow [`rfc-review`](../rfc-review/SKILL.md).
3. **Post the review on the PR.** Use `gh` (in the host execution context when sandboxed
   `gh` lacks Keychain access) to leave your review on PR `#<N>`. Begin the body with the
   exact two header lines printed by `claim` / `listen` so the author agent's watcher can
   identify it even from GitHub alone:

   ```markdown
   ## Peer Review — <agent> (model: <model>, reasoning: <reasoning>)
   <!-- mantle-bloom-peer-review: pr=<N> slot=<slot> agent=<agent> -->
   ```

   Post inline comments + a summary via `gh pr review <N> --comment --body-file <file>` (or
   `gh pr comment <N> --body-file <file>` when reviewing your own GitHub account's PR and
   `--request-changes` is disallowed). Never approve or merge the PR; approval and merging
   belong to the project owner.
4. **Notify the author agent.** Mark your slot complete in the shared handoff record:

   ```bash
   python3 bin/peer_review.py complete --pr <N> --slot <slot> \
     --url "<github-comment-or-review-url>" \
     --summary "<one-line summary of findings>"
   ```

## 4. Respond to completed reviews (author agent)

When the background `bin/peer_review.py wait --pr <N>` task finishes and wakes you:

1. Fetch the full review bodies and inline comments from GitHub:

   ```bash
   gh pr view <N> --comments
   gh api repos/:owner/:repo/pulls/<N>/reviews
   gh api repos/:owner/:repo/pulls/<N>/comments
   ```

2. Summarize both reviewers' findings for the user, noting where they agree, where they
   raised distinct issues, and whether any point is a false positive.
3. After agreeing on the fixes with the user (or applying straightforward fixes directly
   when asked), run [`pre-push-check`](../pre-push-check/SKILL.md), push the updated branch,
   add a `## Review fixes` section to the PR description or reply on the PR, and report the
   outcome.

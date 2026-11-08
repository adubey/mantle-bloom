---
name: start-issue
description: Set up an isolated branch, worktree, venv and ports for a Mantle-Bloom issue. Check its RFC context and tracking sub-issue.
---

# Start work on an issue

1. **Read the issue**, including its comments: `gh issue view <N> --comments`. Note linked
   issues and PRs, any `.mbworld` saves or replay numbers, and whether the issue cites or is
   listed under an RFC. Read the RFC's status and diagnosis gates in `rfc/` and
   `rfc/README.md`. If it is not accepted or its gates have not passed, resolve that RFC
   lifecycle step before implementing the issue. If a `## RFC revision return` comment on
   this issue or its tracking issue postdates this issue's last work, follow
   [Resume](../edit-rfc/SKILL.md#4-resume-an-affected-issue) in `edit-rfc`. If the tracking
   issue links an `## RFC revision brief` whose tracker is still open, raise it with the user
   before starting.
   When work begins, follow [Starting an issue](../project-management/SKILL.md#starting-an-issue)
   to update its project status.

2. **Choose the worktree.** Use the discovery procedure in
   [`continue-rfc`](../continue-rfc/SKILL.md#2-find-rfc-and-evidence-worktrees) to locate RFC
   evidence, not to choose an issue checkout. RFCs are large
   decisions tracked by a parent issue and sub-issues; every issue,
   including an RFC sub-issue, gets its own isolated issue worktree. Do not reuse or switch
   the RFC proposal worktree, and never switch branches in a tree another task may be using.
   Create the issue worktree (an RFC-linked issue is named for both the RFC and the issue)
   with the [`worktree` skill](../worktree/SKILL.md), which owns the path, branch name and
   creation command.

3. **Install deps** in the worktree, as the `worktree` skill describes.

4. **Pick ports** if you'll run the app. Ports 8000 and 5173 belong to the user's own
   instance. Use e.g. `8000+N` for the backend and `9000+N` for the frontend, after checking
   that they're free (`lsof -iTCP:<port> -sTCP:LISTEN`):

   ```bash
   ./bin/restart.sh --port <P>                                 # single process
   ./bin/restart.sh --dev --port <P> --frontend-port <Q>       # Vite HMR
   ./bin/stop.sh --port <P> [--frontend-port <Q>]              # same ports to stop
   ```

5. **If the work needs before/after numbers**, set up a baseline now. Use a task-specific,
   detached `origin/main` worktree made with the [`worktree` skill](../worktree/SKILL.md)
   (baseline, `--backend-only` setup); never share or reuse another task's baseline tree.

   See the `replay-ab` skill for running the comparison.

At completion, follow [Completing an issue](../project-management/SKILL.md#completing-an-issue)
to settle its project status and identify the next eligible RFC step.

---
name: start-issue
description: Set up an isolated branch, worktree, venv and ports for a Mantle Bloom issue. Check its RFC context and tracking sub-issue.
---

# Start work on an issue

1. **Read the issue**, including its comments: `gh issue view <N> --comments`. Note linked
   issues and PRs, any `.mbworld` saves or replay numbers, and whether the issue cites or is
   listed under an RFC. Read the RFC's status and diagnosis gates in `rfc/` and
   `rfc/README.md`. If it is not accepted or its gates have not passed, resolve that RFC
   lifecycle step before implementing the issue.
   When work begins, follow [Starting an issue](../project-management/SKILL.md#starting-an-issue)
   to update its project status.

2. **Choose the worktree.** Use the discovery procedure in
   [`continue-rfc`](../continue-rfc/SKILL.md#find-rfc-and-evidence-worktrees) to locate RFC
   evidence, not to choose an issue checkout. RFCs are large
   decisions tracked by a parent issue and sub-issues; every issue,
   including an RFC sub-issue, gets its own isolated issue worktree. Do not reuse or switch
   the RFC proposal worktree. Name an RFC-linked issue worktree for both the RFC and issue.
   Branch from `origin/main`, check that the branch name is unused locally and remotely,
   and never switch branches in a tree another task may be using. The primary checkout's
   local `main` may lag behind `origin`; never base work on it, and never pull or reset it.

   ```bash
   git fetch origin
   git worktree add -b <agent>/issue-<N>-<slug> ../mantle-bloom-issue-<N> origin/main
   ```

   `<agent>` is the short identifier for the agent/tool doing the work (for example,
   `claude`, `codex`, or `muse`); use your own identifier. `<slug>` is two to four words,
   e.g. `muse/issue-320-crust-transfer`. For an RFC-linked issue, use an RFC-specific path
   while creating the issue branch directly from `origin/main`:

   ```bash
   git worktree add -b <agent>/issue-<N>-<slug> \
     ../mantle-bloom-rfc-<NNNN>-issue-<N> origin/main
   ```

   In this repository, create the issue worktree with the command above. Do not rely on a
   harness-created worktree unless you have verified its base commit is `origin/main`; a
   branch rename does not change its base. If it is based on another commit, create a fresh
   issue worktree instead and leave the harness tree untouched.

3. **Install deps** in the worktree: `./bin/setup.sh`. It creates the worktree's own
   `backend/.venv`. Don't symlink another tree's venv.

4. **Pick ports** if you'll run the app. Ports 8000 and 5173 belong to the user's own
   instance. Use e.g. `8000+N` for the backend and `9000+N` for the frontend, after checking
   that they're free (`lsof -iTCP:<port> -sTCP:LISTEN`):

   ```bash
   ./bin/restart.sh --port <P>                                 # single process
   ./bin/restart.sh --dev --port <P> --frontend-port <Q>       # Vite HMR
   ./bin/stop.sh --port <P> [--frontend-port <Q>]              # same ports to stop
   ```

5. **If the work needs before/after numbers**, set up a baseline now. Use a task-specific,
   detached `origin/main` worktree; never share or reuse another task's baseline tree:

   ```bash
   git worktree add --detach ../mantle-bloom-issue-<N>-baseline origin/main
   (cd ../mantle-bloom-issue-<N>-baseline && ./bin/setup.sh --backend-only)
   ```

   If that path already exists, do not switch or reuse its worktree; choose another unique
   path for this issue.

   See the `replay-ab` skill for running the comparison.

At completion, follow [Completing an issue](../project-management/SKILL.md#completing-an-issue)
to settle its project status and identify the next eligible RFC step.

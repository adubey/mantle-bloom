---
name: start-issue
description: Set up the branch, worktree, venv and ports for a Mantle Bloom issue. Check for an associated RFC and use its worktree when present.
---

# Start work on an issue

1. **Read the issue**, including its comments: `gh issue view <N> --comments`. Note linked
   issues and PRs, any `.mbworld` saves or replay numbers, and whether the issue cites or is
   listed under an RFC. Read the RFC's status and diagnosis gates in `rfc/` and
   `rfc/README.md`. If it is not accepted or its gates have not passed, resolve that RFC
   lifecycle step before implementing the issue.

2. **Choose the worktree.** For an RFC-linked issue, use the local worktree for that RFC.
   Follow `start-rfc` to locate it or set one up for the current machine and agentic
   framework. Reuse its directory, but create a separate issue branch from `origin/main`
   there so each issue remains independently reviewable. Don't switch branches in a dirty
   worktree or one another task is using; set up another local RFC worktree for parallel
   work. For an issue without an RFC, use a new issue worktree. The primary checkout's
   local `main` may lag behind `origin`; never base work on it, and never pull or reset it.

   ```bash
   git fetch origin
   git worktree add -b <agent>/issue-<N>-<slug> ../mantle-bloom-issue-<N> origin/main
   ```

   `<agent>` is `claude` or `codex`; `<slug>` is two to four words, e.g.
   `claude/issue-320-crust-transfer`. For an RFC-linked issue, run the equivalent branch
   creation from inside the RFC worktree, without adding a separate issue worktree:

   ```bash
   git switch -c <agent>/issue-<N>-<slug> origin/main
   ```

   If your harness already created a worktree for you, confirm it is the correct RFC
   worktree when applicable and that the issue branch starts from `origin/main`.

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

5. **If the work needs before/after numbers**, set up the baseline now. It's a detached
   `origin/main` worktree, reused across issues:

   ```bash
   git worktree add --detach ../mantle-bloom-origin-main origin/main \
     || git -C ../mantle-bloom-origin-main checkout --detach origin/main
   (cd ../mantle-bloom-origin-main && ./bin/setup.sh --backend-only)
   ```

   See the `replay-ab` skill for running the comparison.

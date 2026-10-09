---
name: start-issue
description: Set up an isolated worktree, branch, venv and ports before working on a mantle-bloom GitHub issue. Use when starting on an issue or any new line of work in this repo.
---

# Start work on an issue

1. **Read the issue**, including its comments: `gh issue view <N> --comments`. Note linked
   issues and PRs, and any `.mbworld` saves or replay numbers it cites.

2. **Branch off `origin/main`, in a new worktree.** The primary checkout's local `main` is
   usually behind `origin`. Never base work on it, and never pull or reset it.

   ```bash
   git fetch origin
   git worktree add -b <agent>/issue-<N>-<slug> ../mantle-bloom-issue-<N> origin/main
   ```

   `<agent>` is `claude` or `codex`; `<slug>` is two to four words, e.g.
   `claude/issue-320-crust-transfer`. If your harness already created a worktree for you,
   use it, but check that its branch starts from `origin/main`.

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

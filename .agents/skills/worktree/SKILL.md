---
name: worktree
description: Choose the path, branch name and creation command for a Mantle-Bloom git worktree (issue, RFC proposal, RFC continuation, or detached baseline). Use whenever another skill says to create a worktree, and before running `git worktree add`.
---

# Create a worktree

Single source for where worktrees go and how they are made. `start-issue`, `start-rfc`,
`continue-rfc` and `replay-ab` link here rather than repeating it.

## 1. Pick the parent directory

| Your sandbox and file tools are... | Parent directory |
|---|---|
| Scoped to the primary checkout root (for example `gemini`/Antigravity) | `worktrees/` inside the primary checkout (gitignored via `/worktrees/`, excluded from ruff) |
| Free to reach sibling directories (for example `claude`, `codex`, `muse`) | `../mantle-bloom-<name>` next to the primary checkout |

The workspace-local `worktrees/` path avoids an approval prompt on every file edit and
command, because the sandbox and file tools only auto-approve inside the workspace root.
Run `git worktree add` from the primary checkout so `worktrees/` resolves there, not
inside another worktree. If creating the branch fails in the sandbox (`.git` ref writes),
rerun just that command outside it.

## 2. Pick the name

`<name>` below is the directory name under the parent from step 1.

| Purpose | `<name>` | Branch |
|---|---|---|
| Issue | `issue-<N>` | `<agent>/issue-<N>-<slug>` |
| RFC-linked issue | `rfc-<NNNN>-issue-<N>` | `<agent>/issue-<N>-<slug>` |
| RFC proposal (never reused for issues) | `rfc-<NNNN>` | `<agent>/rfc-<NNNN>-<slug>` |
| Continuing an RFC draft whose tree is unavailable | `rfc-<NNNN>-continue-<slug>` | `<agent>/rfc-<NNNN>-continue-<slug>` |
| Before/after baseline (detached) | `issue-<N>-baseline` | none |

`<agent>` is your own short identifier (`claude`, `codex`, `gemini`, `muse`); don't use
another agent's. `<slug>` is two to four words, e.g. `muse/issue-320-crust-transfer`.

## 3. Create it

Branch from a freshly fetched `origin/main`, never the primary checkout's local `main`,
which lags. Check the branch name is unused locally and remotely, and that the path
doesn't exist. If it does, pick another path; never switch or reuse another task's tree.

```bash
git fetch origin
git worktree add -b <agent>/issue-<N>-<slug> <parent>/<name> origin/main
git worktree add --detach <parent>/<name> origin/main          # baseline
git worktree add -b <branch> <parent>/<name> origin/<rfc-branch>  # RFC continuation
```

Don't rely on a harness-created worktree unless its base commit is `origin/main`; a
branch rename doesn't change the base. Otherwise make a fresh one and leave the harness
tree untouched.

## 4. Set it up

Run `./bin/setup.sh` in the new tree (add `--backend-only` for a baseline). Each tree has
its own `backend/.venv`; never symlink another tree's. Pick ports other than 8000 and
5173, which belong to the user's own instance. If setup fails with `ENOSPC`, free disk
before retrying; a partial `.venv` is safe to delete since you created it.

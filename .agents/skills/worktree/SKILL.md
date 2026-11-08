---
name: worktree
description: Choose the path, branch name and creation command for a Mantle-Bloom git worktree (issue, RFC proposal, RFC continuation, or detached baseline). Use whenever another skill says to create a worktree, and before running `git worktree add`.
---

# Create a worktree

Single source for where worktrees go and how they are made. `start-issue`, `start-rfc`,
`continue-rfc` and `replay-ab` link here rather than repeating it.

## 1. Pick the parent directory

Default to a sibling of the primary checkout: `<primary>-<name>` (for example
`../mantle-bloom-issue-320`). Use `<primary>/worktrees/<name>` only if your sandbox and file
tools can't leave the workspace root (for example `gemini`/Antigravity), where every
sibling path would prompt for approval on each edit and command. If unsure, use the
sibling. `/worktrees/` is gitignored and excluded from ruff.

Resolve `<primary>` as an absolute path, so a relative `worktrees/` never nests inside
whichever tree you happen to be in. This works from any worktree and on old git:

```bash
primary="$(cd "$(git rev-parse --git-common-dir)/.." && pwd)"
```

If creating the branch fails in the sandbox (`.git` ref writes), rerun just that command
outside it.

## 2. Pick the name

`<name>` below is the suffix or directory name from step 1.

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

`<path>` is `"$primary-<name>"` or `"$primary/worktrees/<name>"` from step 1.

Branch from a freshly fetched `origin/main`, never the primary checkout's local `main`,
which lags. Check the branch name is unused locally and remotely, and that the path
doesn't exist. If it does, pick another path; never switch or reuse another task's tree.

```bash
git fetch origin
git worktree add -b <agent>/issue-<N>-<slug> <path> origin/main
git worktree add --detach <path> origin/main          # baseline
git worktree add -b <branch> <path> origin/<rfc-branch>  # RFC continuation
```

Don't rely on a harness-created worktree unless its base commit is `origin/main`; a
branch rename doesn't change the base. Otherwise make a fresh one and leave the harness
tree untouched.

## 4. Set it up

Run `./bin/setup.sh` in the new tree (add `--backend-only` for a baseline). Each tree has
its own `backend/.venv`; never symlink another tree's. Pick ports other than 8000 and
5173, which belong to the user's own instance. If setup fails with `ENOSPC`, free disk
before retrying; a partial `.venv` is safe to delete since you created it.

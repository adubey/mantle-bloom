# Mantle Bloom: agent guide

An Earth-like planet simulator. The plate tectonics, climate and biosphere simulation and the
FastAPI routes are in `backend/app/` (Python, numpy/scipy/numba). The React + TypeScript
canvas viewer is in `frontend/src/`. There's one human author, and most changes come from
Claude Code and Codex agents, one GitHub issue per branch.

Before changing the model, read [docs/architecture.md](docs/architecture.md) (module map,
request flow) and [docs/simulation-model.md](docs/simulation-model.md). For offline
diagnostics of saved worlds, see [docs/debugging.md](docs/debugging.md).

## Commands

Run from the repo root (or worktree root).

| What | Command |
|---|---|
| Set up or refresh deps (safe to repeat) | `./bin/setup.sh` |
| Unit tests for what changed (dev loop) | `./bin/affected_test.sh --base origin/main` |
| Full unit suite (required before push) | `./bin/unit_test.sh` (extra args go to pytest, e.g. `-k biome`) |
| Slow full-simulation tests | `./bin/stress_test.sh` |
| Python lint | `backend/.venv/bin/ruff check` |
| Frontend lint, typecheck and build | `cd frontend && npm run lint && npm run build` |
| Run the app | `./bin/restart.sh --port <P>`, or add `--dev --frontend-port <Q>` for Vite HMR |
| Stop the app | `./bin/stop.sh --port <P> [--frontend-port <Q>]` |

- `affected_test.sh` maps changed files to tests through the import graph, and falls back
  to the full suite when that mapping can't be trusted. It's a dev-loop shortcut only.
- The full suite (~1,300 tests, about 12 minutes) and ruff are the gate before pushing. CI
  (`.github/workflows/test.yml`) runs both, plus the frontend checks, on every PR.
- The frontend has no unit tests. `npm run build` runs `tsc`.
- Backend deps are pinned in `backend/constraints.txt`. Change versions only with
  `./bin/pin_deps.sh`, never by hand. Upgrades can move simulation results.

## Workflow

- **One worktree per issue, branched off `origin/main`.** Name the branch
  `<agent>/issue-<N>-<slug>`, where `<agent>` is `claude` or `codex`. The primary checkout's
  local `main` lags `origin`: don't use it as a base or baseline, and don't pull or reset it.
  Run `./bin/setup.sh` in each new worktree. Each worktree has its own `backend/.venv`, so
  don't symlink one from another tree.
- **Ports:** never use 8000 or 5173 (the user's own instance). Pick free ports for each
  worktree so issues can run side by side.
- **Before/after comparisons** use a detached `origin/main` worktree. Never use
  `git stash`: the stash is shared by every worktree, including the user's.
- **Replays are RAM-heavy** (0.6–2.7 GB each). Run at most two at once, and never alongside
  the full unit suite.
- **Reviews:** when a PR's author reports the full suite passing, don't rerun it. Run
  targeted checks they couldn't have covered instead.
- Skills for these workflows, shared by Claude Code and Codex, are in `.agents/skills/`
  (`.claude/skills` links there): `start-issue`, `replay-ab`, `pre-push-check`.

## Conventions

- **Tests:** `backend/unit_tests/test_<module>.py` mirrors `backend/app/<module>.py`. Every
  unit test runs well under a second; anything slower belongs in `backend/stress_tests/`.
  A bug fix comes with a test that fails without it.
- **Comments say why**, and cite the issue (`#314`) that motivated a non-obvious choice, as
  the existing code does.
- **Saves:** `.mbworld` files are pickled `World`s. A change to pickled state must keep old
  saves loading. See [docs/save-compatibility.md](docs/save-compatibility.md) and
  `unit_tests/test_persistence.py`.
- **Docs:** when you add a module, route or model behaviour, update `docs/architecture.md`,
  `docs/api-reference.md` or `docs/simulation-model.md` in the same PR.
- **Investigations:** diagnostic scripts go in `bin/debug/` (checked in, with a usage
  docstring). Findings go in `analysis/issue<N>/README.md`, with only compact JSON checked
  in. Saves and full replay outputs stay out of the repo.
- **Commits:** imperative subject that says what changed in model terms (e.g. "Transfer
  consumed lower-plate crust to the upper plate"), with the issue number when useful.
- **PRs:** a `## Summary` that explains what changed and why, with replay numbers in a
  table when the model moved, then `Closes #<N>`. Add a `## Review fixes` section when
  addressing review comments.

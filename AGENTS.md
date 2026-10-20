# Mantle Bloom: agent guide

An Earth-like planet simulator. The plate tectonics, climate and biosphere simulation and the
FastAPI routes are in `backend/app/` (Python, numpy/scipy/numba). The React + TypeScript
canvas viewer is in `frontend/src/`. There's one human author, and most changes come from
Claude Code and Codex agents, one GitHub issue per branch.

Before changing the model, read [docs/architecture.md](docs/architecture.md) (module map,
request flow), [docs/simulation-model.md](docs/simulation-model.md), and
[docs/design-principles.md](docs/design-principles.md). In simulation issues and PRs, name
the principle the current behavior violates and the invariant or metric that will demonstrate
the fix. For offline diagnostics of saved worlds, see [docs/debugging.md](docs/debugging.md).

## Commands

Run from the repo root (or worktree root).

| What | Command |
|---|---|
| Set up or refresh deps (safe to repeat) | `./bin/setup.sh` |
| Unit tests for what changed (required before push) | `./bin/affected_test.sh --base origin/main` |
| Full unit suite (CI runs it on a schedule) | `./bin/unit_test.sh` (extra args go to pytest, e.g. `-k biome`) |
| Slow full-simulation tests | `./bin/stress_test.sh` |
| Python lint (incl. complexity) | `backend/.venv/bin/ruff check` |
| Frontend lint (incl. complexity), typecheck and build | `cd frontend && npm run lint && npm run build` |
| Dead code and duplicate code | `./bin/quality_check.sh` |
| Run the app | `./bin/restart.sh --port <P>`, or add `--dev --frontend-port <Q>` for Vite HMR |
| Stop the app | `./bin/stop.sh --port <P> [--frontend-port <Q>]` |

- `affected_test.sh` maps changed files to tests through the import graph, and falls back
  to the full suite on its own when that mapping can't be trusted (an `__init__.py`, a data
  file, a dependency pin).
- The affected tests, ruff and `quality_check.sh` are the gate before pushing -- the same
  checks CI (`.github/workflows/test.yml`) runs on every PR, plus the frontend checks. You
  don't need to run the full suite (~1,300 tests, about 12 minutes) locally: it runs four
  times a day and the stress tests nightly (`full-suite.yml`). A failure there opens a
  "Scheduled full test run is failing" issue, which takes priority over new work.
- **Quality gates are ratchets.** Functions over the complexity cap of 20 carry a
  `# noqa: C901` or `oxlint-disable-next-line eslint/complexity`, and pre-existing dead and
  duplicate code is baselined (`vulture_whitelist.py`, `.jscpd-baseline.json`). Don't add to
  any of these to get a check passing: split the function, delete the dead code or factor
  out the copy. Remove an entry when its code goes. Editing inside a baselined clone counts
  as a new clone; refresh with `cd frontend && npm run dupes:baseline` only after deciding
  the copy should stay.
- The frontend has no unit tests. `npm run build` runs `tsc`.
- Backend deps are pinned in `backend/constraints.txt`. Change versions only with
  `./bin/pin_deps.sh`, never by hand. Upgrades can move simulation results.

## Workflow

- **One issue branch and worktree per issue, branched off `origin/main`.** Name the branch
  `<agent>/issue-<N>-<slug>`, where `<agent>` is `claude` or `codex`. RFCs are large
  decisions tracked by a parent issue with sub-issues; each sub-issue gets its own isolated
  issue worktree. Never switch the RFC proposal worktree to an issue branch. Use
  `start-rfc` to begin or advance an RFC workstream, `continue-rfc` to reconstruct an RFC
  handoff, and `rfc-authoring` to choose new RFC numbers and write proposal content. The
  primary checkout's local `main` lags `origin`: don't use it as a base or baseline, and
  don't pull or reset it. Run `./bin/setup.sh` in each new worktree. Each worktree has its own
  `backend/.venv`, so don't symlink one from another tree.
- **Ports:** never use 8000 or 5173 (the user's own instance). Pick free ports for each
  worktree so issues can run side by side.
- **Before/after comparisons** use a task-specific detached `origin/main` worktree. Never
  share or reuse another task's baseline worktree. Never use `git stash`: the stash is
  shared by every worktree, including the user's.
- **Replays and scratch simulation probes are RAM-heavy** (0.6–2.7 GB each). Run only one
  replay or probe at a time, never alongside another such run or a unit or stress test. Follow
  this limit regardless of the skill used.
- **Reviews:** when a PR's author reports its tests passing, don't rerun them. Run targeted
  checks they couldn't have covered instead.
- **Merging is the user's job.** Open the PR and stop: never run `gh pr merge`. `main` only
  takes changes through a PR with the `backend` and `frontend` checks green, and refuses
  force pushes and deletion, for everyone. Every branch outside `claude/*` and `codex/*` is
  the user's: don't push to, rebase, force-push or delete it. Force-push (with
  `--force-with-lease`) only your own agent branch.
- **Never discard work you didn't create.** Don't run `git reset --hard`, `git checkout --`,
  `git restore`, `git clean`, `git stash`, `git branch -D`, `git worktree remove --force` or
  `rm -rf` against anything but your own worktree and scratch files. The user's checkout
  may hold uncommitted work that exists nowhere else.
- **Private local evidence handoffs:** use a `mantle-bloom-handoff/` directory inside the
  shared Git common directory for same-machine paths and other private pointers. Resolve it
  with `handoff="$(cd "$(git rev-parse --git-common-dir)" && pwd)/mantle-bloom-handoff"`,
  then create it with `mkdir -p "$handoff"`; this works with older Git versions and all
  worktrees of the clone share it. Never copy its
  contents or personal paths into a commit, issue, PR, or checked-in report. Keep durable
  run labels, checksums, commands, and findings in the appropriate analysis report or RFC.
- Skills for these workflows, shared by Claude Code and Codex, are in `.agents/skills/`
  (`.claude/skills` links there): `explore-problem`, `create-issue`, `start-issue`,
  `start-rfc`, `continue-rfc`, `project-management`, `replay-ab`,
  `pre-push-check`, `rfc-authoring`, `rfc-review`, `complete-rfc`, `documentation`.

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
  `docs/api-reference.md` or the relevant `docs/simulation-model/<topic>.md` page in the same
  PR. `docs/simulation-model.md` is navigation only; use its module map to choose a topic.
  Use the [`documentation` skill](.agents/skills/documentation/SKILL.md) for documentation edits.
- **Investigations:** diagnostic scripts go in `bin/debug/` (checked in, with a usage
  docstring). Findings go in `analysis/issue<N>/README.md`, with only compact JSON checked
  in. Saves and full replay outputs stay out of the repo.
- **RFCs:** a model change that is large, uncertain or rests on the science starts as an RFC
  in `rfc/NNNN-slug.md` ([rfc/README.md](rfc/README.md) has the index and lifecycle). Use
  `rfc-authoring` for its content and `start-rfc` to advance acceptance, gates and
  implementation PRs; review RFCs with `rfc-review`, and close them out (implemented or
  withdrawn) with `complete-rfc`. Open
  implementation issues after acceptance when their scope is independent of diagnosis
  gates; wait for gate results first when they determine the issue scope or implementation
  approach. Implementation must not begin until required gates pass.
- **Commits:** imperative subject that says what changed in model terms (e.g. "Transfer
  consumed lower-plate crust to the upper plate"), with the issue number when useful.
- **PRs:** a `## Summary` that explains what changed and why, with replay numbers in a
  table when the model moved, then `Closes #<N>`. Add a `## Review fixes` section when
  addressing review comments.
- **Links in chats:** when mentioning a PR or another chat in a chat message, include a
  clickable link to it rather than plain text.

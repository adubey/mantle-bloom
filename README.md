# Mantle Bloom



https://github.com/user-attachments/assets/0535bac6-e602-41e7-82f7-c4a096244443



An Earth-like planet simulator, including plate tectonics, climate and biosphere. The planet can have completely different terrain from Earth but is similar in other ways including size, composition and insolation. 

## Where to look

- **[docs/architecture.md](docs/architecture.md)** -- stack, request flow, how world state
  is held.
- **[docs/simulation-model.md](docs/simulation-model.md)** -- the actual model: plate-local
  frames, mantle flow, boundary evolution, line regularization, merge/split, boundary point
  reassignment, projections, and why each simplification was an acceptable line to draw.
- **[docs/api-reference.md](docs/api-reference.md)** -- the three backend routes.
- **[docs/debugging.md](docs/debugging.md)** -- offline `python -m app.*` diagnostic dumps and
  the `?deb`-gated debug map views/panels for checking a long run's health.
- **[docs/packaging.md](docs/packaging.md)** -- building the single-process desktop binary
  (`bin/package.sh`).

## Getting started

### Prerequisites

- **Python 3.10+** (backend: fastapi, uvicorn, numpy, scipy, numba, pytest, ruff --
  `backend/requirements.txt`, with exact versions pinned in `backend/constraints.txt`).
- **Node.js 20+** (frontend: React + TypeScript + Vite -- `npm ci` in `frontend/` pulls in
  everything, including TypeScript itself).

### First-time setup

```bash
git clone https://github.com/adubey/mantle-bloom/
cd mantle-bloom
./bin/setup.sh
```

Creates `backend/.venv`, installs the pinned backend deps into it, and installs the frontend's
from `package-lock.json`. Safe to re-run after pulling; run it in each new git worktree too.
`PYTHON=/path/to/python3.x ./bin/setup.sh` picks the interpreter. To upgrade the pinned
versions on purpose, run `./bin/pin_deps.sh` and commit the diff.

### Running it

```bash
./bin/restart.sh
```

Builds the frontend and starts the **single-process app** in the background on `:8000` --
FastAPI serves the frontend bundle from its own origin (`backend/app/desktop.py`), the same
setup `bin/package.sh` freezes into a binary. Waits for it to respond, then prints the URL.
Open `http://localhost:8000`, click **Generate World**, then **Step** or **Play**.

For active frontend work, `--dev` runs the **two-process** setup instead -- uvicorn plus
the Vite dev server (hot module reload) on `:5173`:

```bash
./bin/restart.sh --dev   # then open http://localhost:5173
```

`--port` / `--frontend-port` (the latter `--dev` only) move the ports; `--stay-awake` wraps
the processes in `caffeinate`. Logs land in `/tmp/mantle-bloom-<port>.log` (and
`/tmp/mantle-bloom-frontend-<frontend-port>.log` in `--dev`).

Stop everything with `./bin/stop.sh` (pass the same `--port` you started with, if any).

### As a standalone desktop binary

`./bin/package.sh` bundles the single-process app into a self-contained executable (no
Python or Node needed to run it, opens a browser on launch) under `dist/`. Run it on the OS
you want to ship for -- PyInstaller doesn't cross-compile; CI
(`.github/workflows/package.yml`) produces the macOS and Windows builds. See
[docs/packaging.md](docs/packaging.md).

### Running the tests

```bash
./bin/unit_test.sh
```

Runs the backend's fast unit tests (`backend/unit_tests/`, every test well under a second)
-- the frontend has no unit test framework set up yet. Extra arguments are forwarded to
pytest, e.g. `./bin/unit_test.sh -k biome`.

When you're iterating on a handful of files, run only the tests actually affected by your
changes (committed since `main` plus anything uncommitted) instead of the whole suite:

```bash
./bin/affected_test.sh
```

It maps changed files to test files through `app`'s own import graph -- a module's own
tests, plus every module that transitively imports it -- and falls back to the full suite
whenever that mapping can't be trusted. It's the pre-push check, the same one CI runs on
every PR; CI runs the full suite on a schedule.

The suite's slow, full-simulation tests (many-step integration/determinism checks, anywhere
from a few seconds to several minutes each) live separately in `backend/stress_tests/`,
run with:

```bash
./bin/stress_test.sh
```

## Project layout

```
mantle-bloom/
  backend/
    app/             # simulation pipeline + FastAPI routes -- see docs/architecture.md
    unit_tests/      # fast pytest suite (well under a second per test)
    stress_tests/    # slow pytest suite (full-simulation integration/determinism checks)
  frontend/
    src/             # React + TypeScript + Canvas map viewer
  docs/              # you are here
  bin/
    setup.sh         # create backend/.venv and install backend + frontend deps
    pin_deps.sh      # regenerate backend/constraints.txt (deliberate dependency upgrades)
    restart.sh       # start/restart the single-process app (--dev for the Vite HMR setup)
    stop.sh          # stop everything restart.sh started
    package.sh       # build the self-contained desktop binary (see docs/packaging.md)
    unit_test.sh     # run the backend's fast unit test suite
    affected_test.sh # run only the unit tests affected by this branch's changes
    stress_test.sh   # run the backend's slow, full-simulation test suite
```


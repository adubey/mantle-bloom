# mantle-bloom

An Earth-like planet simulator, including plate tectonics, climate and biosphere. The planet can have completely different terrain from Earth but is similar in other ways including size, composition and insolation. 

## Where to look

- **[docs/architecture.md](docs/architecture.md)** -- stack, request flow, how world state
  is held.
- **[docs/simulation-model.md](docs/simulation-model.md)** -- the actual model: plate-local
  frames, mantle flow, boundary evolution, line regularization, merge/split, boundary point
  reassignment, projections, and why each simplification was an acceptable line to draw.
- **[docs/api-reference.md](docs/api-reference.md)** -- the three backend routes.
- **[docs/debugging.md](docs/debugging.md)** -- diagnostic views, endpoints, and the
  `python -m app.plate_diagnostics` offline dump for checking a long run's health.
- **[docs/packaging.md](docs/packaging.md)** -- building the single-process desktop binary
  (`bin/package.sh`); **[docs/HOSTING.md](docs/HOSTING.md)** -- putting it on the public web.

## Getting started

### Prerequisites

- **Python 3.10+** (backend: fastapi, uvicorn, numpy, scipy, pytest, httpx --
  `backend/requirements.txt`).
- **Node.js 18+** (frontend: React + TypeScript + Vite -- `npm install` in `frontend/`
  pulls in everything, including TypeScript itself).

### First-time setup

```bash
git clone <this repo>
cd mantle-bloom

cd backend
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
deactivate
cd ..

cd frontend
npm install
cd ..
```

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
    restart.sh       # start/restart the single-process app (--dev for the Vite HMR setup)
    stop.sh          # stop everything restart.sh started
    package.sh       # build the self-contained desktop binary (see docs/packaging.md)
    unit_test.sh     # run the backend's fast unit test suite
    stress_test.sh   # run the backend's slow, full-simulation test suite
```

## Known limitations

- **v1 is elevation-only.** No climate, hydrology, erosion, or biomes yet -- see
  [docs/simulation-model.md#known-simplifications](docs/simulation-model.md#known-simplifications)
  for this and the other deliberate scoping decisions.
- **Single world, in memory.** World state lives in backend memory only, one world at a
  time -- calling `/world/generate` again (or restarting the backend without saving first)
  replaces/loses whatever was there. Use the **File...** button's Save/Load World to
  persist a world to disk and bring it back later -- see
  [docs/api-reference.md](docs/api-reference.md)'s `/world/save`/`/world/load` (no
  cross-version compatibility promise; it's a pickle of the whole in-memory state, not a
  stable interchange format).

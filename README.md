# mantle-bloom

A sphere-native plate tectonics simulator: plates are spherical polygons (not grid cells)
carrying their own elevation data, driven by Euler poles fit each step to a simple mantle
convection model. v1 is elevation-only, built around a sphere-native representation
specifically to avoid the problems a fixed lat/lon grid runs into once plates carry terrain
across it: pole topology, latitude distortion, and lossy resampling.

## Where to look

- **[docs/architecture.md](docs/architecture.md)** -- stack, request flow, how world state
  is held.
- **[docs/simulation-model.md](docs/simulation-model.md)** -- the actual model: plate-local
  frames, mantle flow, boundary evolution, line regularization, merge/split, boundary point
  reassignment, projections, and why each simplification was an acceptable line to draw.
- **[docs/api-reference.md](docs/api-reference.md)** -- the three backend routes.

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

Starts the backend (FastAPI/uvicorn on `:8000`) and the frontend (on `:5173`) in the
background and waits for both to respond. By default the frontend is built and served as a
production build (`vite build` + `vite preview`); pass `--dev` to run the Vite dev server
(hot module reload) instead:

```bash
./bin/restart.sh --dev
```

Logs land in `/tmp/mantle-bloom-backend.log` and `/tmp/mantle-bloom-frontend.log`. Open
`http://localhost:5173`, click **Generate World**, then **Step** or **Play**.

Stop everything with `./bin/stop.sh`.

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
    restart.sh       # start/restart the backend and frontend server
    stop.sh          # stop everything restart.sh started
    unit_test.sh     # run the backend's fast unit test suite
    stress_test.sh   # run the backend's slow, full-simulation test suite
```

## Known limitations

- **v1 is elevation-only.** No climate, hydrology, erosion, or biomes yet -- see
  [docs/simulation-model.md#known-simplifications](docs/simulation-model.md#known-simplifications)
  for this and the other deliberate scoping decisions.
- **No persistence, single world.** World state lives in backend memory only, one world at
  a time; restarting the backend or calling `/world/generate` again loses whatever was
  there.

# mantle-bloom

A sphere-native plate tectonics simulator: plates are spherical polygons (not grid cells)
carrying their own elevation data, driven by Euler poles fit each step to a simple mantle
convection model. v1 is elevation-only -- a successor to
[plate-sim](https://github.com/adubey/plate-sim)'s grid-based approach, built to fix the
sphere/grid problems documented in that project's own docs (pole topology, latitude
distortion, lossy resampling of plate-carried terrain).

## Where to look

- **[docs/architecture.md](docs/architecture.md)** -- stack, request flow, how world state
  is held.
- **[docs/simulation-model.md](docs/simulation-model.md)** -- the actual model: plate-local
  frames, mantle flow, boundary evolution, garbage collection, merge/split, projections, and
  why each simplification was an acceptable line to draw.
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
./restart.sh
```

Starts the backend (FastAPI/uvicorn on `:8000`) and the frontend dev server (Vite on
`:5173`) in the background and waits for both to respond. Logs land in
`/tmp/mantle-bloom-backend.log` and `/tmp/mantle-bloom-frontend.log`. Open
`http://localhost:5173`, click **Generate World**, then **Step** or **Play**.

Stop everything with `./stop.sh`.

### Running the tests

```bash
cd backend
source .venv/bin/activate
python -m pytest tests/ -q
```

## Project layout

```
mantle-bloom/
  backend/
    app/             # simulation pipeline + FastAPI routes -- see docs/architecture.md
    tests/           # pytest suite
  frontend/
    src/             # React + TypeScript + Canvas map viewer
  docs/              # you are here
  restart.sh         # start/restart the backend and frontend dev server
  stop.sh            # stop everything restart.sh started
```

## Known limitations

- **v1 is elevation-only.** No climate, hydrology, erosion, or biomes yet -- see
  [docs/simulation-model.md#known-simplifications](docs/simulation-model.md#known-simplifications)
  for this and the other deliberate scoping decisions.
- **No persistence, single world.** World state lives in backend memory only, one world at
  a time; restarting the backend or calling `/world/generate` again loses whatever was
  there.
- **No plate-boundary polygon overlay.** The rendered map is elevation-line points only; the
  rough Voronoi outline kept per plate is cosmetic and not kept in sync with boundary
  evolution (see
  [docs/simulation-model.md#known-simplifications](docs/simulation-model.md#known-simplifications)).

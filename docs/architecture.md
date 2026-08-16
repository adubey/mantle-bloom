# Architecture

## Stack

- **Backend:** Python, FastAPI + uvicorn, numpy (rotation/vector math), scipy
  (`SphericalVoronoi` for initial plate seeding, `cKDTree` for boundary-adjacency and
  merge-detection queries, `scipy.cluster.vq.kmeans2` for split clustering), pytest.
- **Frontend:** React + TypeScript via Vite, plain HTML `<canvas>` — no mapping/charting
  library, the same choice plate-sim made and the same reasoning: the data is simple enough
  (points + a hand-rolled colormap) that a library would be more ceremony than the problem
  needs.

## Request flow

```
Browser (App.tsx)
  │
  │  POST /world/generate  { seed, num_plates }
  ▼
FastAPI (main.py)
  │  world.generate_world(seed, num_plates) -- builds the plate mosaic, stores it as the
  │  single in-memory World (see below)
  ▼
  { seed, elapsed_years, num_plates }

Browser then fetches, for whichever projection is selected:
  GET /world/render?projection=behrmann|eckert4
  → every plate's elevation-line nodes, projected to 2D -- see api-reference below

Time-stepping:
  POST /world/step  { years }
  → world.step_world(world, years): refit Euler poles, rotate, evolve boundaries, apply
    topology changes, occasionally garbage-collect -- see simulation-model.md
  → browser re-fetches /world/render
```

The frontend never holds simulation state -- it's a thin client that re-fetches
`/world/render` after every `generate`/`step` call and redraws the canvas from scratch.

<a id="world-state"></a>
## World state (`backend/app/world.py`)

Unlike plate-sim (which keeps up to 5 generated worlds in memory, addressed by id),
mantle-bloom keeps exactly **one** world at a time, in a module-level dict in `main.py` --
simpler, matching the v1 "elevation view only" scope. A `World` holds:

- `plates: list[Plate]` -- see below.
- `mantle_centers` -- the convection-cell flow field (see
  [simulation-model.md#mantle-flow](simulation-model.md#mantle-flow)), fixed for the life of
  the world.
- `elapsed_years`, `steps_since_gc`, `next_plate_id` (a monotonically increasing counter so
  a plate created by a split never collides with an existing id, even after other plates
  have been removed).

Each `Plate` (`backend/app/plates.py`) is:

- `frame` -- a 3x3 rotation matrix, the plate's *entire* position and orientation state.
  Rotating a plate is one matrix multiply; nothing else about the plate needs to change.
- `omega` -- current angular velocity (Euler pole direction x rate), refit from the mantle
  flow field every step.
- `lines: list[ElevationLine]` -- the actual carried terrain, each a set of elevation
  samples at fixed plate-local longitudes along one plate-local latitude. This is the
  central data structure; see
  [simulation-model.md#plate-local-frames](simulation-model.md#plate-local-frames).
- `boundary_local` -- a rough polygon outline (from the initial spherical Voronoi cell),
  kept only for a cosmetic loop overlay. It rotates rigidly along with everything else but
  is *not* consulted by boundary evolution, garbage collection, or merge/split -- those all
  work directly off the elevation-line nodes (see
  [simulation-model.md#boundary-evolution](simulation-model.md#boundary-evolution)).

## The simulation pipeline, module by module

```
geometry.py       unit-vector <-> lat/lon conversion, Rodrigues-formula rotation matrices,
                   plate-local coordinate frames
projections.py     Behrmann and Eckert IV map projections, vectorized
noise.py           cheap smooth sphere noise (sum of sinusoids) for initial terrain texture
plates.py          Plate/ElevationLine data structures, the plate-local lattice sweep
                    shared by generation and merge, initial plate generation via
                    SphericalVoronoi
mantle.py           cubed-sphere convection-cell flow field, per-plate Euler-pole
                    least-squares fit
boundary.py         per-step boundary adjacency detection (k-d tree against every other
                     plate's current nodes), convergent/divergent/transform classification,
                     elevation deltas, line growth/shrinkage
line_regrid.py       periodic line-spacing regularization ("garbage collection")
merge_split.py       plate consumption, continental-collision merging, mantle-flow-driven
                     splitting
world.py             World/Plate orchestration: generate_world, step_world
main.py              FastAPI routes tying it all together
```

See [simulation-model.md](simulation-model.md) for what each of these actually computes and
why.

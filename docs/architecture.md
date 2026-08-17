# Architecture

## Stack

- **Backend:** Python, FastAPI + uvicorn, numpy (rotation/vector math), scipy (`cKDTree`
  for initial plate assignment, boundary-adjacency, and merge-detection queries,
  `scipy.cluster.vq.kmeans2` for split clustering), Pillow (server-side map rendering, see
  `render_image.py`), pytest.
- **Frontend:** React + TypeScript via Vite, plain HTML `<canvas>` — no mapping/charting
  library, the same choice plate-sim made and the same reasoning: the frontend's whole job
  is decoding a PNG the backend already rendered and drawing it (`ctx.drawImage`), which a
  library would be more ceremony than the problem needs.

## Request flow

```
Browser (App.tsx)
  │
  │  POST /world/generate  { seed, num_continents }
  ▼
FastAPI (main.py)
  │  world.generate_world(seed, num_continents) -- builds the plate mosaic (total plate
  │  count is chosen automatically from the seed; num_continents is the UI's continents
  │  slider -- see simulation-model.md#initial-plate-generation), stores it as the single
  │  in-memory World (see below)
  ▼
  { seed, elapsed_years, num_plates, events }

Browser then fetches, for whichever projection/map view/resolution is selected:
  GET /world/render?projection=...&view=...&width=...&height=...
  → a PNG, base64-encoded, rendered entirely server-side -- see api-reference below and
    simulation-model.md#render-image

Time-stepping:
  POST /world/step  { years }
  → world.step_world(world, years): refit Euler poles, rotate, evolve boundaries, apply
    topology changes (at most one collision merge per step, only after a sustained 50-100
    Myr collision -- see simulation-model.md#merge-and-split), occasionally garbage-collect
  → browser re-fetches /world/render, and appends any new `events` to the console
```

The frontend never holds simulation state, and no longer holds any *rendering* logic either
-- it's a thin client that re-fetches `/world/render` after every `generate`/`step` call (or
a projection/map-view change) and draws the returned PNG onto the canvas (`MapCanvas.tsx`;
the event console is the one exception: it just displays whatever `events` the last
`generate`/`step` response carried, which is already the full current log -- see
api-reference.md).

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
- `collision_progress: dict[(int, int), float]` -- sustained-collision tracking for
  merge_split.py, pair of plate ids -> accumulated convergent years (see
  [simulation-model.md#merge-and-split](simulation-model.md#merge-and-split)).
- `events: list[(float, str)]` -- the event log for the UI's console (elapsed_years,
  message), capped at `MAX_EVENT_LOG_LENGTH = 200` entries. Appended to via `World.log_event`.

Each `Plate` (`backend/app/plates.py`) is:

- `frame` -- a 3x3 rotation matrix, the plate's *entire* position and orientation state.
  Rotating a plate is one matrix multiply; nothing else about the plate needs to change.
- `omega` -- current angular velocity (Euler pole direction x rate), refit from the mantle
  flow field every step.
- `lines: list[ElevationLine]` -- the actual carried terrain, each a set of elevation
  samples at fixed plate-local longitudes along one plate-local latitude. This is the
  central data structure; see
  [simulation-model.md#plate-local-frames](simulation-model.md#plate-local-frames).

A plate has no separately-tracked boundary polygon at all -- an earlier version kept one
(`boundary_local`, frozen at generation and rotated rigidly thereafter) purely for the
"Plates" map view's outline overlay, and it visibly drifted out of sync with the real
territory after enough stepping (looking like plates overlapping, since it was never
touched by boundary evolution, merge, or split). `Plate.outline_world()` replaces it: every
render, the outline is traced live from each line's current two endpoints -- the actual
edge boundary evolution maintains -- so it can never be stale (see
[simulation-model.md#boundary-evolution](simulation-model.md#boundary-evolution)).

## The simulation pipeline, module by module

```
geometry.py       unit-vector <-> lat/lon conversion, Rodrigues-formula rotation matrices,
                   plate-local coordinate frames
projections.py     Behrmann and Eckert IV map projections, vectorized
noise.py           cheap smooth sphere noise (sum of sinusoids) for initial terrain texture
plates.py          Plate/ElevationLine data structures, the plate-local lattice sweep
                    shared by generation and merge, initial plate generation (nearest-seed
                    tiling), the live per-plate outline used by the "Plates" map view
mantle.py           cubed-sphere convection-cell flow field, per-plate Euler-pole
                    least-squares fit
boundary.py         per-step boundary adjacency detection (k-d tree against every other
                     plate's current nodes), convergent/divergent/transform classification,
                     elevation deltas, line growth/shrinkage
line_regrid.py       periodic line-spacing regularization ("garbage collection")
merge_split.py       plate consumption, sustained-collision continental merging (50-100 Myr,
                     at most one per step), mantle-flow-driven splitting, event log messages
gaps.py              periodic whole-sphere coverage sweep: absorbs gaps into a bordering
                     plate or spawns a new one where no plate dominates
world.py             World/Plate orchestration: generate_world, step_world
main.py              FastAPI routes
render_image.py      renders /world/render's requested view/resolution to a PNG server-side
                     (see simulation-model.md#render-image)
```

See [simulation-model.md](simulation-model.md) for what each of these actually computes and
why.

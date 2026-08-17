# Architecture

## Stack

- **Backend:** Python, FastAPI + uvicorn, numpy (rotation/vector math), scipy (`cKDTree`
  for initial plate assignment, boundary-adjacency, and merge-detection queries,
  `scipy.cluster.vq.kmeans2` for split clustering), Pillow (server-side map rendering, see
  `render_image.py`), pytest.
- **Frontend:** React + TypeScript via Vite, plain HTML `<canvas>` — no mapping/charting
  library, the same choice plate-sim made and the same reasoning: the frontend's whole job
  is decoding a PNG the backend already rendered and drawing it (`ctx.drawImage`), which a
  library would be more ceremony than the problem needs. Two exceptions: `rotation.ts` — a
  small, dependency-free port of just enough backend geometry/projection math to drive the
  "rotate the planet" drag gesture and preview it live client-side (see
  [simulation-model.md#rotating-the-view](simulation-model.md#rotating-the-view)) — and the
  Plate Inspector view (`PlateInspector.tsx`), which renders and drives its whole interaction
  client-side from raw JSON rather than a baked PNG (see
  [simulation-model.md#plate-inspector](simulation-model.md#plate-inspector)). Neither makes
  this a mapping library, and the real detailed rendering stays entirely server-side for every
  other view.

## Request flow

```
Browser (App.tsx)
  │
  │  POST /world/generate  { seed, continental_fraction, land_fraction }
  ▼
FastAPI (main.py)
  │  world.generate_world(seed, continental_fraction, land_fraction) -- builds the plate
  │  mosaic (total plate count is chosen automatically from the seed; the two fractions are
  │  the UI's generation sliders -- see simulation-model.md#initial-plate-generation),
  │  stores it as the single in-memory World (see below)
  ▼
  { seed, elapsed_years, num_plates, events }

Browser then fetches, for whichever projection/map view/resolution/rotation is selected:
  GET /world/render?projection=...&view=...&width=...&height=...&rotation=...
  → a PNG, base64-encoded, rendered entirely server-side -- see api-reference below and
    simulation-model.md#render-image. `rotation` (see simulation-model.md#rotating-the-view)
    is the map's current view orientation, driven by a long-press-and-drag gesture on the
    canvas -- client-local view state, not simulation state, so it's sent fresh with every
    render call rather than stored server-side.

Time-stepping:
  POST /world/step  { years }
  → world.step_world(world, years): refit Euler poles, rotate, evolve boundaries, apply
    topology changes (at most one collision merge per step, only after a sustained 50-100
    Myr collision -- see simulation-model.md#merge-and-split), occasionally garbage-collect
  → browser re-fetches /world/render, and appends any new `events` to the console

When the "Plate Inspector" map view is active, the browser instead (also on every
generate/step, but never on a projection/rotation-only change) fetches:
  GET /world/plates
  → every plate's outline + metadata as JSON, not a PNG -- see api-reference.md and
    simulation-model.md#plate-inspector. Clicking a plate sends its unprojected true
    lat/lon to GET /world/plate_at, a nearest-node lookup answering "which plate is here."
```

The frontend never holds simulation state, and holds only one small piece of *rendering*
logic -- it's a thin client that re-fetches `/world/render` after every `generate`/`step`
call (or a projection/map-view/rotation change) and draws the returned PNG onto the canvas
(`MapCanvas.tsx`). That component also owns the rotate-the-planet drag gesture -- a
long-press-and-drag interaction that previews a cheap wireframe graticule client-side (via
`rotation.ts`) while dragging, then commits the real rotation and lets the usual
`/world/render` re-fetch replace it with the actual detailed frame (see
simulation-model.md#rotating-the-view); that gesture itself was later extracted into
`rotationDrag.ts`'s `useRotationDrag` hook so `PlateInspector.tsx` could reuse it for a
completely different kind of content (translucent plate ellipses, click-to-select,
Tab/Shift+Tab). The event console is the one other exception: it just displays whatever
`events` the last `generate`/`step` response carried, which is already the full current log
-- see api-reference.md).

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
- `axial_tilt_deg` -- a fixed generation-time property like `seed`, read by `climate.py`'s
  insolation calculation on every future render (see
  [simulation-model.md#climate](simulation-model.md#climate)). The one deliberate exception
  to climate otherwise being fully stateless.
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
                   plate-local coordinate frames, local tangent bases + azimuthal-equidistant
                   projection (for the Plate Inspector's bounding-ellipse fit)
ellipse.py         minimum-volume enclosing ellipse (Khachiyan's algorithm), sphere-agnostic
                   pure 2D math -- see simulation-model.md#plate-inspector
projections.py     Behrmann and Eckert IV map projections, vectorized
noise.py           cheap smooth sphere noise (sum of sinusoids) for initial terrain texture
plates.py          Plate/ElevationLine data structures, the plate-local lattice sweep
                    shared by generation and merge, initial plate generation (nearest-seed
                    tiling), the live per-plate outline used by the "Plates" map view, the
                    Plate Inspector's bounding-ellipse fit and nearest-plate click hit-test
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
climate.py           temperature/wind/currents/humidity/precipitation, computed fresh per
                     render on their own fixed equirectangular grid (see
                     simulation-model.md#climate)
main.py              FastAPI routes
render_image.py      renders /world/render's requested view/resolution to a PNG server-side
                     (see simulation-model.md#render-image and simulation-model.md#climate)
```

See [simulation-model.md](simulation-model.md) for what each of these actually computes and
why.

# Architecture

## Stack

- **Backend:** Python, FastAPI + uvicorn, numpy (rotation/vector math), scipy (`cKDTree`
  for initial plate assignment, boundary-adjacency, and merge-detection queries,
  `scipy.cluster.vq.kmeans2` for split clustering), Pillow (server-side map rendering, see
  `render_image.py`), pytest.
- **Frontend:** React + TypeScript via Vite, plain HTML `<canvas>` — no mapping/charting
  library: the frontend's whole job is decoding a PNG the backend already rendered and
  drawing it (`ctx.drawImage`), which a library would be more ceremony than the problem
  needs. Three exceptions: `rotation.ts` — a
  small, dependency-free port of just enough backend geometry/projection math to drive the
  "rotate the planet" drag gesture and preview it live client-side (see
  [simulation-model.md#rotating-the-view](simulation-model.md#rotating-the-view)) — and the
  Plate Inspector (`PlateInspector.tsx`), River Inspector (`RiverInspector.tsx`), and Lake
  Inspector (`LakeInspector.tsx`) views, which each render and drive their whole interaction
  client-side from raw JSON rather than a baked PNG (see
  [simulation-model.md#plate-inspector](simulation-model.md#plate-inspector),
  [simulation-model.md#river-inspector](simulation-model.md#river-inspector), and
  [simulation-model.md#lake-inspector](simulation-model.md#lake-inspector)). None of this makes
  it a mapping library, and the real detailed rendering stays entirely server-side for every
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
    Myr collision -- see simulation-model.md#merge-and-split), erode elevation and route
    rivers/lakes/glaciers from the world's current climate (every step -- see
    simulation-model.md#erosion, simulation-model.md#hydrology, and
    simulation-model.md#glaciation), relax submerged
    continental crust toward a shelf-or-deep-water target (every step -- see
    simulation-model.md#bathymetry), roll each active volcano's own eruption chance (every
    step -- see simulation-model.md#volcanism), and occasionally fill gaps, detect divergent
    boundaries and spawn any new volcanic fields they warrant, and regularize line spacing,
    and -- on the steps in between -- reassign misplaced boundary points (see
    simulation-model.md#reassignment)
  → browser re-fetches /world/render, and appends any new `events` to the console

When the "Plate Inspector" map view is active, the browser instead (also on every
generate/step, but never on a projection/rotation-only change) fetches:
  GET /world/plates
  → every plate's outline + metadata as JSON, not a PNG -- see api-reference.md and
    simulation-model.md#plate-inspector. Clicking a plate sends its unprojected true
    lat/lon to GET /world/plate_at, a nearest-node lookup answering "which plate is here."

When the "River Inspector" map view is active, the browser instead fetches (same cadence):
  GET /world/rivers
  → every distinct river network's flow-edge segments + mouth metadata, plus the current
    coastline (coastline_segments -- see simulation-model.md#coastline, the same segments
    drawn server-side into the temperature/humidity/precipitation renders below), as JSON,
    not a PNG -- see api-reference.md and simulation-model.md#river-inspector. Grouped fresh
    from world.hydrology_cache/climate_cache on every call (no persistent river identity
    across steps -- see hydrology.group_rivers), empty before the first step. Clicking a
    river sends its unprojected true lat/lon to GET /world/river_at, the same nearest-node
    hit-test pattern as plate_at.

When the "Lake Inspector" map view is active, the browser instead fetches (same cadence):
  GET /world/lakes
  → every currently-visible lake's basin info (floor, outlet, current water level, inflowing
    rivers), as JSON, not a PNG -- see api-reference.md and
    simulation-model.md#lake-inspector. Regrouped fresh from world.hydrology_cache on every
    call, same non-persistent-identity contract as river_id, empty before the first step.
    Clicking anywhere on land -- not just a visible lake -- sends its unprojected true lat/lon
    to GET /world/lake_at, which always resolves to something informative: a lake, a dry but
    real basin, "drains straight to the ocean, no basin here," or open ocean.
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

mantle-bloom keeps exactly **one** world at a time, in a module-level dict in `main.py` --
simpler, matching the v1 "elevation view only" scope. A `World` holds:

- `plates: list[Plate]` -- see below.
- `mantle_centers` -- the convection-cell flow field (see
  [simulation-model.md#mantle-flow](simulation-model.md#mantle-flow)), fixed for the life of
  the world.
- `elapsed_years`, `steps_since_regularize`, `steps_since_reassign` (the two periodic-
  maintenance cadence counters -- see simulation-model.md#line-regularization and
  simulation-model.md#reassignment, deliberately never both due on the same step),
  `next_plate_id` (a monotonically increasing counter so a plate created by a split never
  collides with an existing id, even after other plates have been removed).
- `collision_progress: dict[(int, int), float]` -- sustained-collision tracking for
  merge_split.py, pair of plate ids -> accumulated convergent years (see
  [simulation-model.md#merge-and-split](simulation-model.md#merge-and-split)).
- `volcanic_field_plate_ids: set[int]` -- plate ids currently tracked as an active volcanic
  field, removed once diluted below volcanism.VOLCANO_FRACTION_DORMANT_THRESHOLD (see
  [simulation-model.md#volcanism](simulation-model.md#volcanism)).
- `axial_tilt_deg` -- a fixed generation-time property like `seed`, read by `climate.py`'s
  insolation calculation on every future render (see
  [simulation-model.md#climate](simulation-model.md#climate)).
- `sea_level_m`/`solar_multiplier` -- live-adjustable via `POST /world/controls` (the UI's
  "Controls" window, unlike every other generation-time property here), read fresh by
  `climate.py`/`hydrology.py`/`bathymetry.py`/`render_image.py` on every call rather than
  cached, and forced to an immediate `climate_cache` recompute when changed (see
  api-reference.md) so a render/stats call right after doesn't wait for the next step. These,
  together with `axial_tilt_deg` above, are the only deliberate exceptions to climate
  otherwise being fully stateless.
- `events: list[(float, str)]` -- the event log for the UI's console (elapsed_years,
  message), capped at `MAX_EVENT_LOG_LENGTH = 200` entries. Appended to via `World.log_event`.
- `climate_cache`/`hydrology_cache` -- this step's climate/flow-routing snapshot, populated
  by `erosion.py` (which needs a fresh one every step regardless) and reused by
  `/world/stats`, a climate map render, and river/lake rendering so they don't each trigger
  their own recomputation the same turn (see simulation-model.md#climate and
  simulation-model.md#hydrology). Up to one step stale by design, not a bug.

Each `Plate` (`backend/app/plates.py`) is:

- `frame` -- a 3x3 rotation matrix, the plate's *entire* position and orientation state.
  Rotating a plate is one matrix multiply; nothing else about the plate needs to change.
- `omega` -- current angular velocity (Euler pole direction x rate), refit from the mantle
  flow field every step.
- `lines: list[ElevationLine]` -- the actual carried terrain, each a set of elevation
  samples at fixed plate-local longitudes along one plate-local latitude. This is the
  central data structure; see
  [simulation-model.md#plate-local-frames](simulation-model.md#plate-local-frames). Each
  line also carries `channel_depth`/`channel_width`/`lake_depth`/`silt_depth`/`glacier_depth`
  (persistent, land-only -- see simulation-model.md#hydrology and
  simulation-model.md#glaciation) and `is_volcano`/`volcano_active_years_remaining`
  (persistent -- see simulation-model.md#volcanism) as ordinary parallel arrays right
  alongside `elevation` itself, so they rotate with the plate for free, no advection scheme
  needed. Any such field should be threaded via `dataclasses.replace(line, ...)` at a call
  site that only changes elevation/a value or two (not `theta`), rather than an explicit
  field-by-field reconstruction -- see `plates.ElevationLine`'s own docstring for a real bug
  that pattern caused (`is_volcano` silently wiped every step at two call sites that predated
  it).

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
line_regrid.py       periodic line-spacing regularization
merge_split.py       plate consumption, sustained-collision continental merging (50-100 Myr,
                     at most one per step), mantle-flow-driven splitting, event log messages
gaps.py              periodic whole-sphere coverage sweep: absorbs gaps into a bordering
                     plate or spawns a new one where no plate dominates, event log messages
                     for newly spawned plates
volcanism.py          periodic (same cadence as gaps.py) detection of divergent boundary
                     gaps and new continental "volcanic field" plate spawning, plus every-
                     step eruption/field-lifecycle bookkeeping (see
                     simulation-model.md#volcanism)
reassign.py          periodic pass (staggered against gaps.py's cadence) that hands a node
                     over to a neighboring plate once most of its nearest neighbors belong to
                     it, event log messages for each reassignment
world.py             World/Plate orchestration: generate_world, step_world
climate.py           temperature/wind/currents/humidity/precipitation, computed fresh on
                     their own fixed equirectangular grid -- every render, and now every
                     step too, to drive erosion.py (see simulation-model.md#climate)
erosion.py           every-step rain/river/weathering/glacier erosion + downstream sediment
                     deposition + glacier flattening, elevation deltas driven by climate.py's
                     current fields and hydrology.py's flow routing/glacier state (see
                     simulation-model.md#erosion) -- the weather-influences-geology half of
                     the coupling; climate.py's own elevation-reading mechanics (lapse rate,
                     mountain wind deflection, orographic rain shadow) are the other half
hydrology.py         every-step flow routing over the geology node cloud (a k-nearest-
                     neighbor graph, not a grid): basin-spill detection, steepest-
                     descent flow direction, downstream flow accumulation, glacier
                     accumulation/melt/flow -- feeds erosion.py's river/glacier erosion and
                     deposition and the rendered river/lake/glacier overlay (see
                     simulation-model.md#hydrology and simulation-model.md#glaciation); also
                     groups the flat is_river mask into distinct drainage networks
                     (group_rivers) and answers the River Inspector's click hit-test
                     (river_at), on demand rather than every step (see
                     simulation-model.md#river-inspector)
lakes.py              every-step lake growth/evaporation/merge/split/silt, an explicit n-ary
                     tree of Lake objects built from a depression-hierarchy pass over
                     hydrology.py's own k-NN graph -- called from hydrology.compute_hydrology,
                     projects back down into the same flat lake_depth array every other
                     consumer already reads (see simulation-model.md#lakes-are-an-explicit-tree);
                     also rebuilt fresh, on demand, by main.py to answer the Lake Inspector's
                     GET /world/lakes and GET /world/lake_at (see
                     simulation-model.md#lake-inspector)
bathymetry.py        every-step relaxation of submerged continental crust toward a shelf
                     (near land) or deep-water (far from land) target elevation (see
                     simulation-model.md#bathymetry)
coastline.py          traces the land/ocean and lake boundary as line segments over
                     climate.py's own grid, on demand (not every step) -- drawn into the
                     temperature/humidity/precipitation renders and sent as JSON alongside
                     GET /world/rivers and GET /world/lakes (see simulation-model.md#coastline)
main.py              FastAPI routes
render_image.py      renders /world/render's requested view/resolution to a PNG server-side
                     (see simulation-model.md#render-image and simulation-model.md#climate),
                     plus /world/animate's animated-GIF rendering (render_animation_gif)
persistence.py       whole-World save/load to a single opaque pickle file (File > Save/Load
                     World -- see api-reference.md's /world/save//world/load)
geodesic.py          geodesic-icosahedron hex/pentagon dome tiling + elevation/biome
                     sampling for File > Export Hex Grid (see docs/hex-export-format.md),
                     independent of the plate simulation's own node cloud
```

See [simulation-model.md](simulation-model.md) for what each of these actually computes and
why.

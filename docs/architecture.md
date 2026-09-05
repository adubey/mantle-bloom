# Architecture

## Stack

- **Backend:** Python, FastAPI + uvicorn, numpy (rotation/vector math), scipy (`cKDTree`
  for initial plate assignment, boundary-adjacency, merge-detection, and plate-defragment
  queries, `scipy.cluster.vq.kmeans2` for split clustering,
  `scipy.sparse.csgraph.connected_components` for plate defragmentation), Pillow
  (server-side map rendering, see `render_image.py`), pytest.
- **Frontend:** React + TypeScript via Vite, plain HTML `<canvas>` — no mapping/charting
  library: the frontend's whole job is decoding a PNG the backend already rendered and
  drawing it (`ctx.drawImage`), which a library would be more ceremony than the problem
  needs. Three exceptions: `rotation.ts` — a
  small, dependency-free port of just enough backend geometry/projection math to drive the
  "rotate the planet" drag gesture and preview it live client-side (see
  [simulation-model.md#rotating-the-view](simulation-model.md#rotating-the-view)) — and the
  Plate Inspector (`PlateInspector.tsx`), River Inspector (`RiverInspector.tsx`), Lake
  Inspector (`LakeInspector.tsx`), and Plates & Faults (`PlatesAndFaults.tsx`) views,
  which each render and drive their whole interaction
  client-side from raw JSON rather than a baked PNG (see
  [simulation-model.md#plate-inspector](simulation-model.md#plate-inspector),
  [simulation-model.md#river-inspector](simulation-model.md#river-inspector),
  [simulation-model.md#lake-inspector](simulation-model.md#lake-inspector), and
  [simulation-model.md#fault-inspector](simulation-model.md#fault-inspector)). None of this makes
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
  → world.step_world(world, years): every plate refits its Euler pole and rotates
    (`LithospherePlate.shift`), then every plate reconciles its actual footprint against the
    sphere minus every other live plate's own territory (`LithospherePlate.deform`, in a
    freshly randomized order each turn) -- collision/subduction uplift or trench elevation where a plate's
    rotated territory now overlaps a neighbor's, rift fill (or, occasionally, a fresh
    volcano) where it opens unclaimed space, transform uplift where it's merely close --
    see simulation-model.md#boundary-evolution. Then topology changes (at most one collision
    merge per step, only after a sustained 50-100
    Myr collision -- see simulation-model.md#merge-and-split), erode elevation and route
    rivers/lakes/glaciers from the world's current climate (every step -- see
    simulation-model.md#erosion, simulation-model.md#hydrology, and
    simulation-model.md#glaciation), and roll each active volcano's own eruption chance (every
    step -- see simulation-model.md#volcanism). (Submerged crust's depth is set directly by
    isostasy, not a per-step relaxation; `bathymetry.py` now only grades continent/ocean
    margins once, at generation -- see simulation-model.md#bathymetry.) Line regularization and "claim adjacent
    territory" (a plate growing toward its own pole, or reclaiming ground a subducted
    neighbor vacated) now happen inline inside every `deform()` call rather than on a
    periodic cadence -- the old separate gap-filling and boundary-point-reassignment passes
    are gone entirely (see simulation-model.md#gap-filling and
    simulation-model.md#reassignment)
  → browser re-fetches /world/render, and appends any new `events` to the console

Every POST /world/step call also advances the atmospheric wind solve -- not a separate
mode/endpoint (see simulation-model.md#mode-toggle):
  → world.py's step_world → _advance_fluid_dynamics(world, node_cloud), gated on
    World.simulate_climate_biomes the same way erosion/hydrology already are:
    atmosphere_cfd.step_atmosphere_cfd(world, state, seconds): a genuine time-integrated
    shallow-water solve (real Coriolis/pressure-gradient physics, CFL-stable substepping --
    see simulation-model.md#ocean-atmospheric-fluid-dynamics), a fixed real-time increment per
    tectonics step (one simulated day) regardless of the tectonic `years` requested, against
    terrain refreshed from the world's *current* elevation/climate each step (refresh_forcing)
  → browser re-fetches /world/render as usual (wind/temperature draw off this state; the
    oceanCurrents/humidity/precipitation/biome views are diagnostic climate.py fields -- see
    climate.py's own module docstring)

When the "Plate Inspector" or "Plates & Faults" map view is active, the browser instead
(also on every generate/step, but never on a projection/rotation-only change) fetches:
  GET /world/plates
  → every plate's outline + metadata as JSON, not a PNG -- see api-reference.md and
    simulation-model.md#plate-inspector. Clicking a plate sends its unprojected true
    lat/lon to GET /world/plate_at, a nearest-node lookup answering "which plate is here."
  "Plates & Faults" also fetches GET /world/faults, GET /world/earthquakes and
    GET /world/volcanoes (fault traces + fault systems, recent earthquake epicentres, and
    current volcano vents) and draws them over the plates; a sidebar toggle hides the
    earthquake + volcano overlay. Faults are display-only there -- plate selection is what
    click/Tab drives.

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

GET /world/faults (see faults.py / simulation-model.md#faults) returns every intraplate
fault as a true-frame trace polyline + type/motion/age metadata plus the fault systems and
the shared coastline_segments, as JSON not a PNG -- see api-reference.md. It backs the
"Plates & Faults" view (above). GET /world/fault_at is a nearest-trace hit-test (currently
unused by the frontend since faults aren't selectable there).
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
- `elapsed_years`, `next_plate_id` (a monotonically increasing counter so a plate created by
  a split never collides with an existing id, even after other plates have been removed).
  Line regularization and gap-filling no longer run on a periodic cadence (both happen
  inline inside every `LithospherePlate.deform()` call now -- see
  simulation-model.md#boundary-evolution), so the counters that used to gate them
  (`steps_since_regularize`, `steps_since_reassign`) are gone.
- `steps_taken` -- a plain `step_world` call count (not a year count -- step sizes vary and
  the thing it gates is about accumulated topology drift, not elapsed time). Drives the
  cadence of `merge_split.py`'s geometric plate defragmentation pass
  (`DEFRAG_INTERVAL_STEPS`, see
  [simulation-model.md#merge-and-split](simulation-model.md#merge-and-split)). A plain-int
  dataclass default is a class attribute, so worlds pickled before this field existed still
  load, reading 0.
- `collision_progress: dict[(int, int), float]` -- sustained-collision tracking for
  merge_split.py, pair of plate ids -> accumulated convergent years (see
  [simulation-model.md#merge-and-split](simulation-model.md#merge-and-split)).
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
- `atmosphere_cfd_state` -- the always-on atmospheric wind solver state
  (`atmosphere_cfd.AtmosphereCFDState`), populated once by `generate_world` and never `None`
  again for the rest of that world's life -- not a mode to switch into (see
  simulation-model.md#mode-toggle). Genuinely prognostic (wind + air temperature evolving
  continuously, step to step) rather than recomputed fresh like `climate_cache` above. Ocean
  currents and precipitation are diagnostic in climate.py, not CFD state -- see
  simulation-model.md#ocean-atmospheric-fluid-dynamics.

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
touched by `deform()`, merge, or split). `Plate.outline_world()` replaces it: every render
(and, now, every `deform()` call, which uses this same outline to decide what's contested vs.
open territory), the outline is traced live from each line's current two endpoints -- the
actual edge `deform()` maintains -- so it can never be stale (see
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
elevation_lines.py  ElevationLine data structure, node density/spacing (TARGET_LINE_SPACING_RAD,
                    line_spacing_rad), the plate-local lattice sweep shared by generation and
                    merge, and periodic line-spacing regularization (formerly line_regrid.py)
rtree_index.py      minimal bulk-loaded (STR-packed) R-tree over 2D points -- box/nearest-
                    neighbor queries, used by PlateWithRTree
plates.py          Plate (ABC) / PlateWithLines data structures -- identity, territory,
                    the plate-local lattice, the per-row outline / row-lookup fast path,
                    node iteration and field access -- plus initial plate generation
                    (nearest-seed tiling), the live per-plate outline used by the "Plates"
                    map view and by `deform()`'s own contested/open classification, and
                    the Plate Inspector's bounding-ellipse fit and nearest-plate click
                    hit-test. The tectonic engine itself -- `shift()`/`deform()`, per-turn
                    Euler-pole refit + rotation, polygon-containment boundary
                    classification, elevation/Hc/Hm deltas, line growth/shrinkage, over-
                    stretched-rift volcano spawning, claiming adjacent territory, and
                    inline line regularization -- lives on `LithospherePlate`
                    (lithosphere_plate.py)
mantle.py           cubed-sphere convection-cell flow field, per-plate Euler-pole
                    least-squares fit
boundary.py         `closing_rate` (used only by merge_split.py now, to confirm two
                     continental plates are actively converging, not just neighbors) and the
                     couple of threshold constants merge_split.py shares with plates.py's
                     `deform()`
merge_split.py       plate consumption, sustained-collision continental merging (50-100 Myr,
                     at most one per step), mantle-flow-driven splitting, periodic geometric
                     defragmentation (severed-lobe / stranded-node cleanup deform() can't do),
                     event log messages
gaps.py              whole-sphere coverage maintenance, periodic (world.step_world, same
                     cadence as merge_split's defragmentation): finds any region no plate's
                     lines currently reach -- e.g. ocean floor a fully-subducted plate
                     vacated with no neighbour left nearby to grow into it -- and spawns a
                     new oceanic plate to cover it (`lithosphere_plate.new_plate`).
                     deform()'s own per-step boundary growth only ever extends a line from
                     an existing node, so it structurally can't reach a region with no
                     nearby line at all; this is the periodic backstop for that gap
volcanism.py          every-step eruption lifecycle for existing volcano nodes (active-years
                     countdown, per-step eruption roll, elevation/mineral_deposit_m growth) --
                     volcanic-field *creation* now happens inline inside `deform()`'s own
                     overstretched-rift handling, not a separate periodic detection pass (see
                     simulation-model.md#volcanism)
faults.py            every-step intraplate fault-line lifecycle (stress-weighted Poisson
                     spawn near boundaries, Andersonian regime pick, sub-parallel fault
                     sets, per-regime relief within ~45 km of the trace, lock-up into
                     permanent scars) plus fault *systems* (FaultSystem: a long curving
                     master lineament + a wide strand family, one level above the lone
                     trace) -- an *additive* layer that never touches deform()'s own
                     boundary classification; re-homed across merges/splits by
                     reconcile_faults from world.step_world; backs the "Plates & Faults" map
                     view (GET /world/faults, /world/fault_at) -- see simulation-model.md#faults
worldsketch.py       parses a drawn/loaded coastline image (Generate World's "Human-made" tab)
                     into land/mountain/river masks + plate seed sites -- see
                     simulation-model.md#worldsketch. Consumed by lithosphere_plate.py's
                     generate_plates (its `sketch` param), not part of the per-step pipeline
world.py             World/Plate orchestration: generate_world, step_world
climate.py           temperature/wind/currents/humidity/precipitation, computed fresh on
                     their own fixed equirectangular grid -- every render, and now every
                     step too, to drive erosion.py (see simulation-model.md#climate)
erosion.py           every-step rain/river/weathering/glacier erosion + downstream sediment
                     deposition + glacier flattening + coastal planation/infill feedback
                     (near-sea-level wave-cut planation + sheltered-shelf silting, emergent
                     barrier islands), elevation deltas driven by climate.py's
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
stranded_basins.py   diagnostic-only: finds endorheic below-sea-level basins with no ocean
                     drainage (the "land-locked coastal pit") in lakes.py's forest, and
                     tracks how long each has persisted across steps (world.stranded_basin_
                     tracks, reconciled from world.step_world) -- backs GET /world/stranded_
                     basins and the python -m app.stranded_basins offline dump (see
                     debugging.md)
bathymetry.py        the shelf-width constant geology.py keys off, plus a one-off
                     generation-time pass (shape_initial_bathymetry) that drowns submerged
                     continental interiors toward abyssal depth by distance from land and
                     grades every continent/ocean plate margin into a slope, so the seabed
                     reads as real relief rather than a bright shelf and hard cliffs (see
                     simulation-model.md#bathymetry)
coastline.py          traces the land/ocean and lake boundary as line segments over
                     climate.py's own grid, on demand (not every step) -- drawn into the
                     temperature/humidity/precipitation renders and sent as JSON alongside
                     GET /world/rivers and GET /world/lakes (see simulation-model.md#coastline)
fluid_dynamics.py      shared numerical primitives for atmosphere_cfd.py: real physical-unit
                     gradients/Laplacians/divergence, real Coriolis parameter, CFL-stable
                     substep sizing, semi-Lagrangian advection (see
                     simulation-model.md#ocean-atmospheric-fluid-dynamics)
atmosphere_cfd.py      the atmospheric wind solver: a genuine time-integrated (not diagnostic)
                     shallow-water simulation of wind + air temperature -- prognostic state
                     (World.atmosphere_cfd_state) that persists/evolves step to step, plus a
                     sustained latitude-banded wind forcing. Ocean currents and precipitation
                     are diagnostic in climate.py, not CFD-solved -- the shallow-water ocean
                     solver was retired (see simulation-model.md#ocean-atmospheric-fluid-dynamics)
main.py              FastAPI routes
render_image.py      renders /world/render's requested view/resolution to a PNG server-side
                     (see simulation-model.md#render-image and simulation-model.md#climate),
                     plus /world/animate's streamed H.264/MP4 rendering (stream_animation_mp4)
                     and the "speckle" coastal-dither debug overlay (see docs/debugging.md)
persistence.py       whole-World save/load to a single opaque pickle file (File > Save/Load
                     World -- see api-reference.md's /world/save//world/load)
geodesic.py          geodesic-icosahedron hex/pentagon dome tiling + elevation/biome
                     sampling for File > Export Hex Grid (see docs/hex-export-format.md),
                     independent of the plate simulation's own node cloud
```

See [simulation-model.md](simulation-model.md) for what each of these actually computes and
why.

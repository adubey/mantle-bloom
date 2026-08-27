# API Reference

Routes over a single in-memory world (see
[architecture.md#world-state](architecture.md#world-state)) -- no world id, one world at a
time.

## `POST /world/generate`

Request body:

```json
{
  "seed": 1,
  "num_plates": null,
  "continental_fraction": 0.7,
  "land_fraction": 0.29,
  "axial_tilt_deg": 23.5,
  "num_mantle_centers": 8,
  "node_density": 4.0,
  "initial_soil_maturity": 0.0,
  "climate_density": 4.0,
  "fluid_density": 4.0
}
```

All fields optional. `num_plates` defaults to `null` -- omit it (or check the frontend's
"Auto" box, the default) to let the world tile itself into a plausible plate count from
`seed` alone (`plates.MIN_AUTO_PLATES` to `plates.MAX_AUTO_PLATES`, the same range the UI's
own plate-count slider covers when "Auto" is unchecked, see
[simulation-model.md#initial-plate-generation](simulation-model.md#initial-plate-generation)),
or pass an explicit count to override it. `continental_fraction` and `land_fraction` are the
UI's two generation sliders (0 to 1, defaulting in the frontend to
`plates.DEFAULT_CONTINENTAL_FRACTION = 0.70` and `plates.DEFAULT_LAND_FRACTION = 0.29`);
both optional, each falling back to its own default behavior when omitted (an independent
per-plate coin flip for crust type, a fixed elevation floor for land) -- see
[simulation-model.md#initial-plate-generation](simulation-model.md#initial-plate-generation).
`axial_tilt_deg` is the UI's third generation slider (degrees), defaulting to
`world.DEFAULT_AXIAL_TILT_DEG = 23.5` (Earth's real tilt) -- doesn't affect plate generation,
only `climate.py`'s insolation on future renders (see
[simulation-model.md#climate](simulation-model.md#climate)), which is why it's stored on
`World` rather than consumed once here. `num_mantle_centers` defaults to
`world.DEFAULT_MANTLE_CENTERS = 8`. `node_density` is the UI's "point density" choice
(`plates.NODE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0)`, defaults to
`plates.DEFAULT_NODE_DENSITY = 4.0`) -- `400` if it isn't one of those four values; `2.0` is
a lower-resolution middle ground (half the default), `0.5` the coarsest, fastest option.
`initial_soil_maturity` is the UI's
fifth generation slider (0 to 1, defaults to `0.0` -- a fully barren starting world, no soil
on any land node) -- a one-time seed for `soil_depth`/`soil_mineral_content`/
`soil_organic_content` (see
[simulation-model.md#resources-and-soil](simulation-model.md#resources-and-soil)), not
stored on `World` afterward. `climate_density` is the UI's "climate & biome resolution"
choice (`climate.CLIMATE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0)`, defaults to
`climate.DEFAULT_CLIMATE_DENSITY = 4.0`) -- `400` if it isn't one of those four values.
Scales `climate.py`'s own simulation grid (and, scaled the same way, the Biome/Combined/
Resources/Soil-Quality views' own finer render grid) in *each* dimension -- the default `4.0`
quadruples it (16x the reference cell count), for the sharpest climate/biome maps; `0.5`
halves it, for a coarser but faster grid -- stored on `World` (unlike `initial_soil_maturity`)
since every future step/render
reads it again, same reasoning `node_density`'s own storage gives (see
[simulation-model.md#climate](simulation-model.md#climate)). `fluid_density` is the UI's
"Fluid dynamics resolution" Advanced-settings choice -- same `climate.CLIMATE_DENSITY_CHOICES`
set (capped at `2.0`/"High" rather than `4.0`/"Very High", see `climate.FLUID_DENSITY_CHOICES`'s
own comment for why) and `400` validation as `climate_density`, but independent of it: sizes
only the atmospheric wind-solver grid (`atmosphere_cfd.init_atmosphere_cfd`, immediately
populated by this same call -- see below), not the climate/biome render grid or erosion's own
climate sampling. Lets a world keep a sharp climate/biome grid while running the wind solve at
a coarser (faster) resolution, or vice versa -- see
[simulation-model.md#ocean-atmospheric-fluid-dynamics](simulation-model.md#ocean-atmospheric-fluid-dynamics).
Replaces whatever world previously existed.

`generate_world` also immediately seeds this new world's permanent
`World.atmosphere_cfd_state` (`atmosphere_cfd.init_atmosphere_cfd`) -- the atmospheric wind
solve is always on, not a mode entered later, see
[simulation-model.md#mode-toggle](simulation-model.md#mode-toggle). Ocean currents and
precipitation are diagnostic in `climate.py`, not CFD-solved.

Response: a summary --

```json
{
  "seed": 1,
  "elapsed_years": 0.0,
  "num_plates": 13,
  "events": [{ "elapsed_years": 0.0, "message": "World generated with 13 plates (4 continental)." }]
}
```

`events` is the *entire* current event log (capped at `world.MAX_EVENT_LOG_LENGTH = 200`
entries, oldest dropped first), not just what changed this call -- simplest for the frontend,
which just replaces its displayed console with it on every response, and small enough not to
matter on the wire. See
[simulation-model.md#merge-and-split](simulation-model.md#merge-and-split) for what generates
an event (merges, plates disappearing, splits -- a collision merely starting is deliberately
not logged) and how each one is worded.

## `POST /world/step`

Request body:

```json
{ "years": 2000000 }
```

Advances the current world by `years` (see
[simulation-model.md](simulation-model.md) for what a step actually does) -- including, gated
on `World.simulate_climate_biomes`, advancing `World.atmosphere_cfd_state` by a fixed
real-time increment (one simulated day) regardless of `years`, see
[simulation-model.md#mode-toggle](simulation-model.md#mode-toggle). Returns the same summary
shape as `/world/generate`, with `events` reflecting anything logged up through this step.
`404` if no world has been generated yet.

## `GET /world/save`

The "File > Save World" download -- the *entire* current world (every plate/line, mantle
field, caches, event log -- see [architecture.md#world-state](architecture.md#world-state))
pickled as a single opaque `application/octet-stream` file (`Content-Disposition:
attachment`), not JSON (see `backend/app/persistence.py`). Deliberately makes no promise of
compatibility across app versions -- pickling by class identity means a later
renamed/restructured field on `World` breaks old files, an accepted trade for "just get me
back what I had," not a stable interchange format (contrast with `/world/export_hexgrid`
below). `404` if no world has been generated yet.

Loading a file back is equivalent to running arbitrary code from its bytes (a standard
pickle caveat) -- acceptable given this server is a single-user localhost dev tool already
(see `main.py`'s CORS allowlist), the same trust boundary every other route here assumes.

## `POST /world/load`

The "File > Load World" upload -- the raw bytes of a file `/world/save` previously
produced, as the request body (`Content-Type: application/octet-stream`, not JSON).
Replaces whatever world previously existed, same as `/world/generate`. Returns the same
summary shape `/world/generate` does. `400` if the bytes aren't a valid mantle-bloom world
file (pickle can raise many different exception types on malformed or foreign input, so
this is caught broadly).

## `POST /world/controls`

Request body (all fields optional, independently settable -- the "Controls" window sends
only whichever control the user touched):

```json
{ "sea_level_m": 500.0, "solar_multiplier": 1.1, "simulate_plate_movement": true, "simulate_climate_biomes": true, "wind_model": "cfd" }
```

Live-adjusts `World.sea_level_m` (default `0.0`), `World.solar_multiplier` (default `1.0`,
scales `climate.SUNLIGHT`), `World.simulate_plate_movement`,
`World.simulate_climate_biomes` (both default `true`), and/or `World.wind_model`
(`"cfd"` default, or `"diagnostic"` -- a `400` for any other value) on the *current* world
-- no regenerate needed. Unlike `axial_tilt_deg`/`node_density`, these are meant to be tweaked
mid-simulation: every `is_ocean` check in the codebase (`climate.py`, `hydrology.py`,
`bathymetry.py`) keys off `sea_level_m` instead of a bare `elevation <= 0.0`, and
`render_image.py`'s elevation-view hypsometric coloring shifts by it too, so raising sea
level visibly floods the elevation map immediately.

`simulate_plate_movement`/`simulate_climate_biomes` let the user run just one half of
`world.step_world` -- see that function's own docstring. With `simulate_plate_movement`
`false`, a step skips every plate's `shift`/`deform` (rotation, boundary classification,
elevation deltas, growth/shrink/claim, inline regularization -- see
simulation-model.md#boundary-evolution), topology changes, and volcanism. With
`simulate_climate_biomes` `false`, a
step skips erosion (and the `climate.compute_climate` call inside it), hydrology,
bathymetry, and resource formation -- by far the most expensive part of a step. Either can
be turned off independently to watch (or speed up) just the other half; `elapsed_years`
always advances regardless of both.

`wind_model` picks the wind field feeding `climate.py`: `"cfd"` runs the shallow-water
solve every step (`_advance_fluid_dynamics`), `"diagnostic"` skips it and rebuilds wind /
air temperature from `climate.py`'s closed-form ABL formulas -- much faster steps for
~85-90% of the CFD biome map. See simulation-model.md#wind-model.

Forces an immediate `climate.compute_climate` recompute (stored back onto
`world.climate_cache`) so the very next `/world/render` or `/world/stats` call reflects the
change without waiting for a step -- `climate_cache` is otherwise only refreshed once per
step by `erosion.py`, and not at all while `simulate_climate_biomes` is `false`. `404` if no
world has been generated yet.

Response echoes back the world's current values for all five:

```json
{ "sea_level_m": 500.0, "solar_multiplier": 1.1, "simulate_plate_movement": true, "simulate_climate_biomes": true, "wind_model": "cfd" }
```

## `GET /world/render?projection=behrmann|eckert4&view=elevation|plates|platesDetail|combined|temperature|wind|oceanCurrents|humidity|precipitation|biome|resources|soilQuality|oceanCfdSediment|oceanCfdDeposition&width=1100&height=611&rotation=1,0,0,0,1,0,0,0,1`

Renders the current world as a PNG, base64-encoded. All drawing (elevation fill, plate-color
fill, boundary outlines, pole markers, rotation arcs, per-node dots) happens server-side
(see [simulation-model.md#render-image](simulation-model.md#render-image)) -- the client
just decodes and paints the image, it never sees raw coordinate data. `400` for an
unrecognized projection/view name, a width/height outside `[1, main.MAX_RENDER_DIMENSION_PX]`
(4000), or a malformed `rotation`, `404` if no world has been generated yet.

```json
{
  "projection": "behrmann",
  "elapsed_years": 2000000,
  "image_base64": "iVBORw0KGgoAAAANSUhEUgAA..."
}
```

- `view` selects what gets drawn: `"elevation"` (colored by height/depth), `"plates"`
  (colored by owning plate, plus boundary outlines/pole markers/rotation arcs), or
  `"platesDetail"` (each plate's raw elevation-line nodes as dots, colored by elevation,
  plus boundary outlines) -- the frontend's Map View dropdown picks this directly. `"biome"`
  and `"combined"` (biome-colored land, hypsometric ocean, lakes/glaciers/rivers overlaid) are
  a categorical classification (see
  [simulation-model.md#biomes](simulation-model.md#biomes)). Five more views come from
  `climate.py` (see [simulation-model.md#climate](simulation-model.md#climate)):
  `"temperature"`, `"humidity"`, and `"precipitation"` are heatmaps (each also drawing the
  current coastline -- see [simulation-model.md#coastline](simulation-model.md#coastline) --
  since a color-scale view carries no land/ocean cue on its own); `"wind"` and
  `"oceanCurrents"` draw subsampled direction/magnitude arrows, and `"oceanCurrents"`
  additionally marks detected ocean swells with small circles. `"wind"` and the land side of
  `"temperature"` draw off `World.atmosphere_cfd_state` (a real, continuously time-integrated
  shallow-water wind solve, see
  [simulation-model.md#ocean-atmospheric-fluid-dynamics](simulation-model.md#ocean-atmospheric-fluid-dynamics))
  -- or, under `wind_model == "diagnostic"` (see
  [simulation-model.md#wind-model](simulation-model.md#wind-model)), off `climate.py`'s own
  closed-form wind/air-temperature instead; everything else is a `climate.py` diagnostic (the
  shallow-water ocean solver was retired) resampled nearest-cell onto the world's HEALPix grid. `"resources"` and `"soilQuality"` (see
  [simulation-model.md#resources-and-soil](simulation-model.md#resources-and-soil)) are
  node-cloud-derived like elevation/plates, not climate-grid-derived -- `"resources"` overlays
  coal/oil & gas/mineral deposit richness on a muted land/ocean backdrop, `"soilQuality"` is a
  continuous fertility heatmap (barren to rich) plus the coastline overlay. `"oceanCfdSediment"`/
  `"oceanCfdDeposition"` (from the retired ocean solver) are not valid `view` values
  (`/world/render` rejects them with `400`, same as any other unrecognized view name).
- `rotation` is the map's current view orientation (see
  [simulation-model.md#rotating-the-view](simulation-model.md#rotating-the-view)): a
  row-major 3x3 rotation matrix as 9 comma-separated floats, applied to every real-world
  position immediately before it's projected. Optional, defaults to identity (today's
  behavior, center at lat=0/lon=0). Purely a render-time transform -- it's never stored on
  `World` and has no bearing on climate simulation results, which key off the true,
  un-rotated planetary frame regardless of what orientation is currently being viewed.
- `width`/`height` are the returned image's exact pixel dimensions. The frontend requests
  more than its canvas's displayed CSS size (`RENDER_SCALE` in `App.tsx`) for a sharper,
  retina-style render at the same on-screen footprint; line widths, dot/pole radii, and
  padding all scale with the requested width (see `render_image.py`'s `pixel_scale`) so a
  higher-resolution request doesn't also make those features look thinner.
- `image_base64` decodes to a PNG. The frontend builds `data:image/png;base64,<this>` as an
  `<img>` source and draws it onto the canvas with `drawImage` -- see `MapCanvas.tsx`.

## `POST /world/animate`

The "File > Make Animation" action -- an animated GIF of `view`/`projection`'s progress, one
frame for the world's current state plus more, each `years_per_frame` further along:

```json
{
  "projection": "eckert4",
  "view": "elevation",
  "width": 1100,
  "height": 611,
  "rotation": null,
  "years_per_frame": 1000000,
  "num_frames": 20
}
```

Same fields/validation as `/world/render`'s own query params (`rotation` is the same
9-comma-separated-float string, `400` for an unrecognized `projection`/`view` or an
out-of-range `width`/`height`), plus `num_frames` bounded to `[1,
main.MAX_ANIMATION_FRAMES]` (60) -- each frame costs a full simulation step and render, so
this bounds worst-case request time. `404` if no world has been generated yet.

**This permanently advances the world** by `(num_frames - 1) * years_per_frame` years, the
same as calling `/world/step` that many times in a row -- deliberately not a
side-effect-free preview (see `render_image.render_animation_gif`'s own docstring). Response
is the same summary shape `/world/generate`/`/world/step` return, plus:

```json
{ "image_base64": "R0lGODlhAAA..." }
```

`image_base64` decodes to an animated GIF (`data:image/gif;base64,<this>`), each frame
quantized against the first frame's own color palette so static regions don't flicker
between playback frames.

## `GET /world/plates`

The "Plate Inspector" map mode's data source (see
[simulation-model.md#plate-inspector](simulation-model.md#plate-inspector)) -- unlike
`/world/render`, this returns plain JSON, not a baked PNG; `frontend/src/PlateInspector.tsx`
renders and drives the interaction itself. Un-rotated/true-frame throughout (no `rotation`
param -- the client applies its current view rotation only at draw time). `404` if no world
has been generated yet.

```json
{
  "elapsed_years": 2000000,
  "plates": [
    {
      "plate_id": 0,
      "crust_type": "continental",
      "num_rows": 85,
      "num_points": 3437,
      "outline": [[0.98, 0.12, -0.05], ["..."]],
      "points": [[0.981423, 0.117582, -0.052207], ["..."]],
      "bounding_ellipse": {
        "center_xyz": [0.97, 0.15, -0.04],
        "diameter_a_km": 13478.2,
        "diameter_b_km": 7148.9,
        "outline": [[0.91, 0.28, -0.04], ["..."]]
      }
    }
  ]
}
```

- `num_rows` is the count of `ElevationLine`s with at least one node (a plate can carry
  zero-node placeholder rows -- see `merge_split.py`'s own consumption checks -- excluded
  here to match `outline_world()`'s own filtering). `num_points` is the total node count
  across every row.
- `outline` traces the plate's live territorial boundary (`Plate.outline_world()`) as a
  closed loop of world-space unit vectors.
- `points` is every one of the plate's `num_points` node positions individually (not just the
  outline loop) -- what the Plate Inspector plots per node, bright for the selected plate and
  dim for every other one. All coordinates in this response (`outline`, `points`,
  `bounding_ellipse`'s fields) are rounded to 6 decimal places before serializing -- far finer
  than anything visible on screen, but it keeps the payload down now that `points` alone can
  run to tens of thousands of entries for a whole world.
- `bounding_ellipse` is `null` for an empty plate (`num_points == 0`), otherwise the minimum-
  area ellipse enclosing every one of the plate's nodes (see
  [simulation-model.md#plate-inspector](simulation-model.md#plate-inspector) for the fitting
  method) -- `diameter_a_km`/`diameter_b_km` are the major/minor diameters in real km,
  `outline` is ~72 sampled world-space points around its perimeter for drawing.

## `GET /world/plate_at?lat_deg=0&lon_deg=0`

The Plate Inspector's click hit-test: which plate owns the node nearest `(lat_deg, lon_deg)`
-- both in the *true* (un-rotated) frame; the client unprojects its click through whatever
view rotation is active first (see `rotation.ts`'s `unproject`/`matTranspose`). `400` for
non-finite input, `404` if no world has been generated yet.

```json
{ "plate_id": 3 }
```

`plate_id` is `null` only if every plate is empty (shouldn't happen via the API, but a
freshly-constructed `World` with no plates at all has no node to be nearest to).

## `GET /world/rivers`

The "River Inspector" map mode's data source (see
[simulation-model.md#river-inspector](simulation-model.md#river-inspector)) -- same "plain
JSON, client renders it" contract as `/world/plates`; `frontend/src/RiverInspector.tsx`
renders and drives the interaction itself. Un-rotated/true-frame throughout, same as
`/world/plates`. Regrouped fresh from `world.hydrology_cache` on every call (see
`hydrology.group_rivers`) rather than persisted, so `river_id` is only meaningful against this
same response -- not a stable identity across steps. `rivers`/`coastline_segments` are both
`[]` before the first step (`hydrology_cache`/`climate_cache` are `None` until `erosion.py`
runs once). `404` if no world has been generated yet.

```json
{
  "elapsed_years": 2000000,
  "rivers": [
    {
      "river_id": 0,
      "num_nodes": 42,
      "segments": [[[0.98, 0.12, -0.05], [0.97, 0.13, -0.05]], ["..."]],
      "mouth_xyz": [0.91, 0.28, -0.04],
      "mouth_type": "ocean",
      "flow_rate": 3821.4,
      "speed": 2.67,
      "num_tributaries": 3
    }
  ],
  "coastline_segments": [[[0.91, 0.28, -0.04], [0.90, 0.29, -0.03]], ["..."]]
}
```

- `segments` is a flat list of flow edges (a node and its own `flow_target`), each a pair of
  world-space points -- not an ordered polyline, since a river network can branch. Every
  coordinate in this response is rounded to 6 decimal places before serializing, same as
  `/world/plates`.
- `mouth_xyz` is the network's own highest-`flow_accum` node (flow only accumulates downhill,
  so nothing else in the network can out-flow it) -- what the River Inspector draws a ring
  around for the selected river.
- `mouth_type` is `"ocean"`, `"lake"` (checked first -- a river can end at a still-draining
  lake), or `"other"` (a dry interior sink with no lake standing there).
- `flow_rate` is `flow_accum` at the mouth; `speed` is `hydrology.compute_river_speed`
  evaluated there. Both are stylized, unitless quantities -- meaningful only relative to other
  rivers in the same world, not real physical units (see
  [simulation-model.md#hydrology](simulation-model.md#hydrology)).
- `num_tributaries` counts separate headwater branches feeding the network -- see
  [simulation-model.md#river-inspector](simulation-model.md#river-inspector) for the exact
  definition (an original one; no precedent to port).
- `coastline_segments` is the land/ocean/lake boundary (see
  [simulation-model.md#coastline](simulation-model.md#coastline)), same segment-pair shape as
  `segments` above -- included because this view draws no filled backdrop at all, so without
  it there's no land/ocean/lake cue in the River Inspector whatsoever. The exact same segments
  are also drawn server-side into the `temperature`/`humidity`/`precipitation` views returned
  by `GET /world/render` (no separate endpoint or query param -- it's baked directly into
  those PNGs, same as rivers/lakes are baked into the elevation view).

## `GET /world/river_at?lat_deg=0&lon_deg=0`

The River Inspector's click hit-test: which river network owns the node nearest
`(lat_deg, lon_deg)` -- same true-frame contract as `/world/plate_at`. `river_id` is an index
into the most recent `/world/rivers` response, not a persistent id. `400` for non-finite
input, `404` if no world has been generated yet.

```json
{ "river_id": 2 }
```

`river_id` is `null` if there are no rivers at all this step (before the first step, or a
world with no land steep/wet enough to route any).

## `GET /world/lakes`

The "Lake Inspector" map mode's data source (see
[simulation-model.md#lake-inspector](simulation-model.md#lake-inspector)) -- same "plain
JSON, client renders it" contract as `/world/rivers`; `frontend/src/LakeInspector.tsx` renders
and drives the interaction itself. A "lake" here is one connected component of
`hydrology.LAKE_MIN_VISIBLE_DEPTH_M`-deep nodes -- exactly what's drawn as standing water
everywhere else in this codebase -- regrouped fresh from `world.hydrology_cache` on every call,
so `lake_id` is only meaningful against this same response, not a persistent identity across
steps. `lakes`/`coastline_segments` are both `[]` before the first step. `404` if no world has
been generated yet.

```json
{
  "elapsed_years": 12000000,
  "lakes": [
    {
      "lake_id": 0,
      "is_lake": true,
      "member_count": 6,
      "member_xyz": [[0.35, -0.11, -0.93], ["..."]],
      "floor_xyz": [0.35, -0.11, -0.93],
      "floor_elevation_m": 1170.9,
      "outlet_xyz": [0.37, -0.09, -0.93],
      "outlet_elevation_m": 1902.8,
      "water_elevation_m": 1902.8,
      "is_spilling": true,
      "inflow_rivers": [{ "mouth_xyz": [0.36, -0.10, -0.93], "flow_rate": 812.4, "num_nodes": 9 }]
    }
  ],
  "coastline_segments": [[[0.91, 0.28, -0.04], [0.90, 0.29, -0.03]], ["..."]]
}
```

- `member_xyz` is every currently-*flooded* node's own position (not the whole geometric
  catchment, which routinely includes dry higher ground on the way down to the basin's floor --
  see `main._lake_basin_summary`'s own comment) -- what the Lake Inspector plots per node,
  bright for the selected lake and dim for every other one.
- `floor_xyz`/`floor_elevation_m` is this lake's own lowest point.
- `outlet_xyz`/`outlet_elevation_m` is "the lowest point of the edge of the basin" -- the
  saddle a river out of it would source from, resolved from this step's freshly-rebuilt
  depression hierarchy (`lakes.build_lake_hierarchy`) so a lake that's currently the result of
  several smaller basins merging together reports its own *current* rim, not an interior saddle
  that's already flooded. `null` for an unresolved closed/endorheic basin with no known spill
  (a legitimate state, not missing data -- see `lakes.py`'s own docstring).
- `water_elevation_m`/`is_spilling` describe this lake's current water surface and whether
  it's actively overflowing its own outlet.
- `inflow_rivers` is every `RiverInfo` (see `/world/rivers`) whose own mouth lands somewhere in
  this basin, as a small inline summary rather than a full `RiverSummary` -- the River
  Inspector's own map mode is what draws a river's full path.
- `coastline_segments` is the same land/lake-vs-ocean boundary `/world/rivers` returns, included
  here for the same reason: this view draws no filled backdrop of its own.

## `GET /world/lake_at?lat_deg=0&lon_deg=0`

The Lake Inspector's click hit-test: which *basin* owns the node nearest `(lat_deg, lon_deg)`
-- same true-frame contract as `/world/river_at`, except it always resolves to something
informative for any land click, not just a hit on a currently-visible lake, so "click a basin
that has no water in it right now" still returns real information (the user's own spec: "give
me information about the basin nonetheless"). `400` for non-finite input, `404` if no world has
been generated yet.

```json
{ "kind": "lake", "basin": { "lake_id": 0, "is_lake": true, "...": "same shape as /world/lakes" } }
```

`kind` is one of:

- `"lake"` -- the nearest node is currently flooded; `basin` is that lake's own info (`lake_id`
  set, matching the same node's `/world/lakes` entry for this same step).
- `"basin"` -- dry land, but part of a real (possibly still-unresolved/endorheic) basin;
  `basin` is that basin's own leaf catchment info, same shape as a lake entry but `is_lake:
  false`, `lake_id: null`, `member_xyz` covering the whole geometric catchment (there's no
  flooded subset to narrow it to), and `water_elevation_m`/`is_spilling` both meaningless
  (`null`/`false`).
- `"no_basin"` -- dry land whose own steepest-descent chain drains straight to the ocean
  without ever passing through a local minimum first -- an ordinary hillslope, never part of
  any basin. `basin` is `null`.
- `"ocean"` -- the nearest node is open ocean. `basin` is `null`.

## `POST /world/export_hexgrid`

The "File > Export Hex Grid" action: tiles the sphere into a geodesic-icosahedron
hex/pentagon dome (independent of the plate simulation's own node cloud -- see
`backend/app/geodesic.py`), samples the current world's elevation/ocean/biome onto each
tile, and returns the whole tiling as JSON, for use in another application.

```json
{ "frequency": 16 }
```

`frequency` (optional, defaults to `geodesic.DEFAULT_FREQUENCY = 16`) must be one of
`geodesic.FREQUENCY_CHOICES = (8, 16, 32)` -- 642 / 2,562 / 10,242 tiles respectively (exactly
`10 * frequency**2 + 2`, the standard geodesic-icosahedron vertex count) -- same "a few sane
presets, not a free-form input" reasoning `node_density`/`climate_density` already use.
`400` for any other value, `404` if no world has been generated yet.

Response is the export file directly (not wrapped in base64 -- it's already a small,
text-safe JSON payload, unlike the PNG/GIF routes above) -- see
[hex-export-format.md](hex-export-format.md) for the full shape and, since the user-facing
point of this route is a file another application will parse, pseudocode for how a client
finds a tile's neighbors from it.

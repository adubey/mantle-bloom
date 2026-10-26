# API Reference

Routes over the server's worlds (see
[architecture.md#world-state](architecture.md#world-state)). Every world-scoped admin route
lives under `/worlds/{world_id}/...`: the server holds several worlds in one process, and
the world is always addressed explicitly (#563). Registry routes (`GET /worlds`,
`POST /worlds`, `PATCH /worlds/{world_id}`, `DELETE /worlds/{world_id}`) manage the
records themselves.

Every `/worlds/*` route requires the admin: a GitHub session for the configured owner, or a
request made on the server's own machine. Anyone else gets `401`. Cross-site writes get
`403`. See [architecture.md#access-control](architecture.md#access-control).

## `GET /healthz`

`{"status": "ok"}`, open to anyone, and doesn't wait for the world lock. It's Railway's deploy
health check (#376; see [hosting.md](hosting.md)). It can't answer before
startup has restored the checkpoint. `503` with the reason in `detail` when a checkpoint
directory is configured but couldn't be created or written.

## Login routes

- `GET /login`: the sign-in page. It redirects to `/` when the request is already the admin.
  `?error=denied|state|github|forbidden|unconfigured` shows a fixed message for a failed
  attempt.
- `GET /auth/github`: starts GitHub OAuth (with PKCE). It sets a short-lived signed `state`
  cookie and redirects to GitHub.
- `GET /auth/github/callback`: checks `state`, exchanges the code and fetches the GitHub
  user. If the user ID is the configured owner, it sets the session cookie and redirects to
  `/`. Otherwise it redirects to `/login?error=...` without a session.
- `POST /logout`: clears the session cookie and redirects (`303`) to `/view`.
- `GET /auth/session`: `{"admin": bool, "local": bool, "login": str | null,
  "login_available": bool}`. `local` means the no-login local rule applied. `login` is the
  GitHub login of a real session.

## Public viewer routes

Open to anyone, for the logged-out `/view` page (#373). See
[architecture.md#public-viewer](architecture.md#public-viewer).

- `GET /public/worlds` lists every generated public world, plus any generated
  world whose run loop is on even when private, so anyone can watch a planet
  that is currently stepping (#597): `{"worlds": [{"world_id", "name",
  "running", "revision", "seed", "elapsed_years"}]}`. A private world drops back
  out of the list when its run stops.
- `GET /public/worlds/{world_id}/world`: `{"world": {"revision", "seed", "elapsed_years"},
  "running": bool, "steps_remaining": int | null, "step_limit": 30, "step_years": 100000}`.
  `steps_remaining` is null (unlimited) for the admin. `revision` is an
  opaque string that changes whenever the world does (#374). `404` for a private,
  unknown, or still-empty world -- the listing never advertises the other cases.
- `GET /public/worlds/{world_id}/world?after=<revision>`: the same response, but held open until the world's
  revision differs from `after`, or for at most 25 s (then answered with the unchanged
  world). Viewers poll this to pick up steps that other viewers or the owner take.
- `GET /public/worlds/{world_id}/render?view=elevation|biome&width=&height=`: like `/worlds/{world_id}/render`, but always
  Eckert IV at the default orientation, with width and height at most 2400. It returns
  `{"revision", "seed", "elapsed_years", "image_base64"}`. Each view and size is rendered once
  per revision and then served from a cache, even while the owner holds the world. `400` for
  any other view or size, `404` before a world exists, `503` if a render is needed and the
  owner holds the world for more than 10 s.
- `POST /public/worlds/{world_id}/step` with `{"revision": str}`, the revision the viewer is showing: advances
  the world by 100 kyr and returns `{"world": {...}, "computed": bool, "steps_remaining"}`.
  Only a request that runs the step (`computed: true`) uses up allowance, and the response then
  sets the updated `mb_guest` cookie. When `revision` is stale, or another guest's step is
  already running, the request gets the resulting world without computing. `429` when a
  computation is needed but the allowance is spent, `503` while the owner is using the world,
  `404` for a private, unknown, or still-empty world.

## World registry routes

World identity and management metadata live outside the pickled `World` (#563),
so saves carry no registry information:

- `GET /worlds` returns every record: metadata plus, when populated, the world
  summary and run state.
- `GET /worlds/suggest-name` returns `{"name"}`: an unused exoplanet-style pick
  for the Generate World dialog to offer up front (advisory only; duplicate
  display names are allowed).
- `POST /worlds` with `{"name"?, "is_public"?}` creates an empty record
  (`201`); the record holds no world until Generate or Load populates it. An
  omitted name defaults to an exoplanet-style pick (a star name plus a lowercase
  discovery letter, `c`/`d`/`e` weighted more, #597).
  World IDs are server-generated UUIDs. `409` past `MANTLE_BLOOM_MAX_WORLDS`
  records (default 4).
- `PATCH /worlds/{world_id}` renames (`{"name"}`) or changes visibility
  (`{"is_public"}`).
- `DELETE /worlds/{world_id}` stops the world and removes its record and
  checkpoint. `409` while a generate, step, or animation is in progress on it.

## `POST /worlds/{world_id}/generate`

Request body:

```json
{
  "seed": 1,
  "cratons_first": false,
  "synthetic_relief": false,
  "surface": "quad",
  "num_plates": null,
  "continental_fraction": 0.7,
  "land_fraction": 0.29,
  "axial_tilt_deg": 23.5,
  "num_mantle_centers": 8,
  "node_density": 4.0,
  "initial_soil_maturity": 0.0,
  "climate_density": 4.0,
  "fluid_density": 4.0,
  "name": null,
  "is_public": null
}
```

All fields optional. `num_plates` defaults to `null` -- omit it (or check the frontend's
"Auto" box, the default) to let the world tile itself into a plausible plate count from
`seed` alone (`plates.MIN_AUTO_PLATES` to `plates.MAX_AUTO_PLATES`, the same range the UI's
own plate-count slider covers when "Auto" is unchecked, see
[simulation-model/foundations.md#initial-plate-generation](simulation-model/foundations.md#initial-plate-generation-platespy)),
or pass an explicit count to override it. `continental_fraction` and `land_fraction` are the
UI's two generation sliders (0 to 1, defaulting in the frontend to
`plates.DEFAULT_CONTINENTAL_FRACTION = 0.70` and `plates.DEFAULT_LAND_FRACTION = 0.29`);
both optional, each falling back to its own default behavior when omitted (an independent
per-plate coin flip for crust type, a fixed elevation floor for land). `continental_fraction`
is a share of plates on today's plates-first worlds and a share of the surface on
cratons-first ones (`plates.CONTINENTAL_AREA_FRACTION = 0.5` when omitted), a layout the
selected by `cratons_first` -- see
[simulation-model/foundations.md#initial-plate-generation](simulation-model/foundations.md#initial-plate-generation-platespy).
`cratons_first` is the Random tab's optional terrain-generation choice. When true, continental
cratons are laid out before plates and relief is generated from that layout; when false or
omitted, generation uses the plate-first layout (omission follows
`continent_layout.CRATONS_FIRST_DEFAULT`). Sketches and premade worlds ignore this choice.
`synthetic_relief` applies only to a cratons-first layout: when true, continental relief and
coastlines come from a synthetic tectonic history (sutures between cratons, cordilleras on
leading margins, interior terrain types), backed by crustal thickness; when false or omitted,
the cratons-first relief is generated as before (omission follows
`synthetic_relief.SYNTHETIC_RELIEF_DEFAULT`). `name` and `is_public` are the Generate
World dialog's name field and Public checkbox: applied to the record being generated
into (a null name keeps the record's name). See
[simulation-model/foundations.md#synthetic-relief](simulation-model/foundations.md#relief-from-a-synthetic-tectonic-history-synthetic_reliefpy).
`axial_tilt_deg` is the UI's third generation slider (degrees), defaulting to
`world.DEFAULT_AXIAL_TILT_DEG = 23.5` (Earth's real tilt) -- doesn't affect plate generation,
only `climate.py`'s insolation on future renders (see
[simulation-model/climate.md#climate](simulation-model/climate.md#climate-climatepy)), which is why it's stored on
`World` rather than consumed once here. `num_mantle_centers` defaults to
`world.DEFAULT_MANTLE_CENTERS = 8`. `node_density` is the UI's "point density" choice
(`plates.NODE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0)`, defaults to
`plates.DEFAULT_NODE_DENSITY = 4.0`) -- `400` if it isn't one of those four values; `2.0` is
a lower-resolution middle ground (half the default), `0.5` the coarsest, fastest option.
`initial_soil_maturity` is the UI's
fifth generation slider (0 to 1, defaults to `0.0` -- a fully barren starting world, no soil
on any land node) -- a one-time seed for `soil_depth`/`soil_mineral_content`/
`soil_organic_content` (see
[simulation-model/surface-processes.md#resources-and-soil](simulation-model/surface-processes.md#resources-and-soil-geologypy-plus-volcanismpys-own-eruption-roll)), not
stored on `World` afterward. `climate_density` is the UI's "climate & biome resolution"
choice (`climate.CLIMATE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0)`, defaults to
`climate.DEFAULT_CLIMATE_DENSITY = 4.0`) -- `400` if it isn't one of those four values.
Scales `climate.py`'s own simulation grid (and, scaled the same way, the Biome/Combined/
Resources/Soil-Quality views' own finer render grid) in *each* dimension -- the default `4.0`
quadruples it (16x the reference cell count), for the sharpest climate/biome maps; `0.5`
halves it, for a coarser but faster grid -- stored on `World` (unlike `initial_soil_maturity`)
since every future step/render
reads it again, same reasoning `node_density`'s own storage gives (see
[simulation-model/climate.md#climate](simulation-model/climate.md#climate-climatepy)). `fluid_density` is the UI's
"Fluid dynamics resolution" Advanced-settings choice -- same `climate.CLIMATE_DENSITY_CHOICES`
set (capped at `2.0`/"High" rather than `4.0`/"Very High", see `climate.FLUID_DENSITY_CHOICES`'s
own comment for why) and `400` validation as `climate_density`, but independent of it: sizes
only the atmospheric wind-solver grid (`atmosphere_cfd.init_atmosphere_cfd`, immediately
populated by this same call -- see below), not the climate/biome render grid or erosion's own
climate sampling. Lets a world keep a sharp climate/biome grid while running the wind solve at
a coarser (faster) resolution, or vice versa -- see
[simulation-model/climate.md#ocean-atmospheric-fluid-dynamics](simulation-model/climate.md#atmospheric-wind-solver-fluid_dynamicspy-atmosphere_cfdpy).
Replaces whatever world previously existed.

`sketch` (the Generate World dialog's "Human-made" tab) is optional and omitted by the
"Random" tab entirely -- `{"image_base64": "<a PNG>"}`, the drawn or loaded coastline (see
[simulation-model/foundations.md#worldsketch](simulation-model/foundations.md#human-made-worlds-sketch-driven-generation-worldsketchpy) for the ink convention and
what it changes). When present, land/sea for continental crust comes from this image instead
of noise (`land_fraction` above is then ignored), and plate sites/boundaries are fit to its
landmasses and open ocean rather than scattered at random. `400` if the bytes aren't a
decodable image.

`generate_world` also immediately seeds this new world's permanent
`World.atmosphere_cfd_state` (`atmosphere_cfd.init_atmosphere_cfd`) -- the atmospheric wind
solve is always on, not a mode entered later, see
[simulation-model/climate.md#mode-toggle](simulation-model/climate.md#always-on-not-a-mode). Ocean currents and
precipitation are diagnostic in `climate.py`, not CFD-solved.

Response: a summary --

```json
{
  "seed": 1,
  "elapsed_years": 0.0,
  "num_plates": 13,
  "surface": "quad",
  "events": [{ "elapsed_years": 0.0, "message": "World generated with 13 plates (4 continental)." }],
  "replay_logs": []
}
```

`surface` is always `"quad"` (plates are sparse adaptive quad patches), except that a loaded
save with no plates reports `"empty"`.

The summary arrives as the final line of a newline-delimited JSON stream
(`application/x-ndjson`): `{"type": "progress", "fraction": f}` at the three phase
boundaries, an overdue-gated `{"type": "heartbeat", "fraction": f}` tick wherever a slow
sub-phase runs with no other output (#611), then the summary above as `{"type": "done", ...}`.

`events` is the *entire* current event log (capped at `world.MAX_EVENT_LOG_LENGTH = 200`
entries, oldest dropped first), not just what changed this call -- simplest for the frontend,
which just replaces its displayed console with it on every response, and small enough not to
matter on the wire. See
[simulation-model/tectonics.md#merge-and-split](simulation-model/tectonics.md#merge-and-split) for what generates
an event (merges, plates disappearing, splits, a plate fragmenting into disconnected
landmasses or shedding stranded nodes -- a collision merely starting is deliberately not
logged) and how each one is worded.

`replay_logs` lists the records this server process's replay log wrote, in the same
`{elapsed_years, message}` shape (newest 50; empty unless `MANTLE_BLOOM_REPLAY_LOG_DIR` is
set, see [`GET /worlds/{world_id}/replay-logs`](#get-worldreplay-logs)). They aren't world state, so they
stay out of `events` and out of saves.

## `POST /worlds/{world_id}/step`

Request body:

```json
{ "years": 2000000 }
```

Advances the current world by `years` (see
[simulation-model.md](simulation-model.md) for what a step actually does) -- including, gated
on `World.simulate_climate_biomes`, advancing `World.atmosphere_cfd_state` by a fixed
real-time increment (one simulated day) regardless of `years`, see
[simulation-model/climate.md#mode-toggle](simulation-model/climate.md#always-on-not-a-mode). Returns the same summary
shape as `/worlds/{world_id}/generate`, with `events` reflecting anything logged up through this step.
`400` for a `years` that isn't positive, finite, and at most `main.MAX_STEP_YEARS`
(100 Myr, #557), `404` if no world has been generated yet.

Like generate, the response streams newline-delimited JSON: `{"type": "progress",
"fraction": f}` after each per-plate shift/deform unit, an overdue-gated
`{"type": "heartbeat", "fraction": f}` tick at the slow tail's sub-boundaries (#611),
then the summary as `{"type": "done", ...}`.

## `POST /worlds/{world_id}/run`, `GET /worlds/{world_id}/run`

Play and Stop (#375, see [architecture.md#run-and-checkpoint](architecture.md#run-and-checkpoint)).

```json
{ "running": true, "step_years": 100000 }
```

`running: true` has the server step the world by `step_years` repeatedly, pausing 0.4 s
between steps, until a `running: false` request (or a failed step) stops it. Omit
`step_years` to keep the current size; changing it while running applies from the next step.
Stop takes effect before the next step; a step already under way finishes. The status
survives a restart when checkpointing is configured. `400` for a `step_years` that isn't a
positive number, `404` for `running: true` with no world, `500` if the status took effect but
couldn't be saved to the checkpoint (it's retried with the next world checkpoint; until then a
restart would go back to the previously saved status). `step_years` gets the same validation as
`/worlds/{world_id}/step`'s `years` (`400` unless positive, finite, and at most 100 Myr, #557).

Both return:

```json
{ "running": true, "step_years": 100000.0, "stepping": false, "error": null, "run_version": 14,
  "world": { "revision": "3f9a1c2e-12", "seed": 0, "elapsed_years": 1200000.0 } }
```

`stepping` is true while the loop is taking a step; it can stay true just after Stop, until
the step under way lands. `error` says why the run last stopped on its own (a step raised, or
there was no world), and clears on the next `POST`. `run_version` moves on every change to the
status or `error`, and when a step finishes after Stop. `world` is `null` before one exists.

`GET /worlds/{world_id}/run?after=<revision>&run_version=<n>` is a long poll like `/public/worlds/{world_id}/world`'s: it
answers when the revision moves past `after` or the run state moves past `run_version` (so a
tab sees Stop pressed in another one), or after 25 s. Without `run_version`, only the revision
counts.

## `GET /worlds/{world_id}/replay-logs`

The records the server's replay log has on disk, newest first (#442, see
[architecture.md#replay-log](architecture.md#replay-log)). Behind the `/worlds/{world_id}/*` guard like
every world route.

```json
{ "enabled": true, "directory": "/Users/me/.mantle-bloom/replay-logs/8442",
  "records": [{ "name": "20261021T044502.881373Z-seed42-4.00myr-craton_step_loss",
                "path": "/Users/me/.mantle-bloom/replay-logs/8442/20261021T044502.881373Z-seed42-4.00myr-craton_step_loss",
                "triggers": ["craton_step_loss"], "seed": 42, "elapsed_years_before": 4000000.0,
                "step_years": 1000000.0, "recorded_at": "2026-10-21T04:45:02.881373+00:00" }] }
```

Without `MANTLE_BLOOM_REPLAY_LOG_DIR`: `{"enabled": false, "directory": null, "records": []}`.
[debugging/replay-logs.md](debugging/replay-logs.md) covers replaying a record.

## `GET /worlds/{world_id}/save`

The "File > Save World" download -- the *entire* current world (every plate surface, mantle
field, caches, event log -- see [architecture.md#world-state](architecture.md#world-state))
pickled as a single opaque `application/octet-stream` file (`Content-Disposition:
attachment`), not JSON (see `backend/app/persistence.py`). The file is a versioned envelope
that declares the world's plate surface (`quad`, or `empty` with no plates). It is not a stable interchange
format (contrast with `/worlds/{world_id}/export_hexgrid` below): pickling by class identity ties it to
the code. [save-compatibility.md](save-compatibility.md) lists which older saves this build
reads and how line-backed saves move to sparse quads. `404` if no world has been generated
yet.

Loading a file back is equivalent to running arbitrary code from its bytes (a standard
pickle caveat) -- acceptable given this server is a single-user localhost dev tool already
(see `main.py`'s CORS allowlist), the same trust boundary every other route here assumes.

## `POST /worlds/{world_id}/load`

The "File > Load World" upload -- the raw bytes of a file `/worlds/{world_id}/save` previously
produced, as the request body (`Content-Type: application/octet-stream`, not JSON).
Replaces whatever world previously existed, same as `/worlds/{world_id}/generate`. Returns the same
summary shape `/worlds/{world_id}/generate` does; its `surface` is `"quad"` or `"empty"`.

A line-backed save (written before #251 retired that surface) is converted to sparse quads on
load, one way (see [save-compatibility.md](save-compatibility.md)). The summary's
`surface_conversion` then reports what the conversion changed; it is `null` for a world that
was never converted, and stays with the world through later saves.

`400` with `detail` `"invalid or incompatible world file: <reason>"` when the bytes aren't a
readable save (corrupt, truncated, or not a mantle-bloom file) or are in a format version
newer than this build reads (the save envelope or one of its sparse-quad plates). Nothing
falls back to a partial load.

`503` while a generation, step or animation is in progress, so a running step can't keep
changing a world that has already been replaced.

## `POST /worlds/{world_id}/controls`

Request body (all fields optional, independently settable -- the "Controls" window sends
only whichever control the user touched):

```json
{ "sea_level_m": 500.0, "solar_multiplier": 1.1, "ice_age_period_years": 300000, "simulate_plate_movement": true, "simulate_climate_biomes": true, "wind_model": "cfd", "fault_deformation_mode": "fault" }
```

Live-adjusts `World.sea_level_m` (default `0.0`), `World.solar_multiplier` (default `1.0`,
scales `climate.SUNLIGHT` and the `insolation_w_m2` diagnostic), `World.ice_age_period_years` (default `0.0` = disabled; the full
period in years of a glacial↔interglacial temperature cycle -- see
simulation-model/hydrology.md#glaciation, negative is clamped to `0`), `World.simulate_plate_movement`,
`World.simulate_climate_biomes` (both default `true`), `World.wind_model`
(`"diagnostic"` default, or `"cfd"` -- a `400` for any other value), and/or
`World.fault_deformation_mode` (`"boundary"` default, or `"fault"` / `"both"`) on the
*current* world
-- no regenerate needed. Unlike `axial_tilt_deg`/`node_density`, these are meant to be tweaked
mid-simulation: every `is_ocean` check in the codebase (`climate.py`, `hydrology.py`,
`bathymetry.py`) keys off `sea_level_m` instead of a bare `elevation <= 0.0`, and
`render_image.py`'s elevation-view hypsometric coloring shifts by it too, so raising sea
level visibly floods the elevation map immediately.

`simulate_plate_movement`/`simulate_climate_biomes` let the user run just one half of
`world.step_world` -- see that function's own docstring. With `simulate_plate_movement`
`false`, a step skips every plate's `shift`/`deform` (rotation, boundary classification,
elevation deltas, growth/shrink/claim, inline regularization -- see
simulation-model/tectonics.md#boundary-evolution), topology changes, and volcanism. With
`simulate_climate_biomes` `false`, a
step skips erosion (and the `climate.compute_climate` call inside it), hydrology,
bathymetry, and resource formation -- by far the most expensive part of a step. Either can
be turned off independently to watch (or speed up) just the other half; `elapsed_years`
always advances regardless of both.

`wind_model` picks the wind field feeding `climate.py`: `"cfd"` runs the shallow-water
solve every step (`_advance_fluid_dynamics`), `"diagnostic"` skips it and rebuilds wind /
air temperature from `climate.py`'s closed-form ABL formulas -- much faster steps; see
simulation-model/climate.md#wind-model for how closely it matches the CFD (#529).

`fault_deformation_mode` (`"boundary"` default, or `"fault"` / `"both"` -- a `400` for any
other value) picks where plate-boundary deformation is applied: `"boundary"` is the smooth
distance-band thickening at the polygon edge, `"fault"` localises it onto active fault traces
(and scales the fault-relief layer up to compensate), `"both"` runs both. See
simulation-model/tectonics.md#faults.

The body also accepts the fourteen **geomorphic-budget tuning multipliers**
(`world.TUNING_MULTIPLIER_FIELDS`): `rain_erosion_multiplier`, `river_erosion_multiplier`,
`wind_erosion_multiplier`, `ocean_erosion_multiplier`, `coastal_leveling_multiplier`,
`glacier_erosion_multiplier`, `seismic_erosion_multiplier`, `river_deposition_multiplier`,
`ocean_deposition_multiplier`, `collision_uplift_multiplier`,
`collision_uplift_reach_multiplier`, `volcanism_multiplier`, `fault_relief_multiplier`,
`plume_multiplier` (plume production rate; `0` thins every scheduled event away, so the
run is a matched plume-off twin -- see RFC 0007). Each is a dimensionless scale on
one erosion / deposition / uplift / volcanism process, default `1.0` (== untuned, bit for
bit), and a negative value is a `400`. See
[simulation-model/surface-processes.md#tuning-knobs](simulation-model/surface-processes.md#live-tuning-knobs-worldtuning_multiplier_fields) for what each one scales
and measured effect sizes.

Forces an immediate `climate.compute_climate` recompute (stored back onto
`world.climate_cache`) so the very next `/worlds/{world_id}/render` or `/worlds/{world_id}/stats` call reflects the
change without waiting for a step -- `climate_cache` is otherwise only refreshed once per
step by `erosion.py`, and not at all while `simulate_climate_biomes` is `false`. `404` if no
world has been generated yet.

Response echoes back the world's current values for all six controls plus all twelve tuning
multipliers:

```json
{ "sea_level_m": 500.0, "solar_multiplier": 1.1, "simulate_plate_movement": true, "simulate_climate_biomes": true, "wind_model": "cfd", "fault_deformation_mode": "boundary",
  "rain_erosion_multiplier": 1.0, "river_erosion_multiplier": 1.0, "wind_erosion_multiplier": 1.0, "ocean_erosion_multiplier": 1.0,
  "coastal_leveling_multiplier": 1.0, "glacier_erosion_multiplier": 1.0, "seismic_erosion_multiplier": 1.0, "river_deposition_multiplier": 1.0,
  "ocean_deposition_multiplier": 1.0, "collision_uplift_multiplier": 1.0, "collision_uplift_reach_multiplier": 1.0, "volcanism_multiplier": 1.0 }
```

## `GET /worlds/{world_id}/render?projection=behrmann|eckert4|globe&view=elevation|plates|platesDetail|speckle|combined|temperature|wind|oceanCurrents|humidity|precipitation|biome|resources|soilQuality|geomorph|elevReason|overlapAge|oceanCfdSediment|oceanCfdDeposition&width=1100&height=611&rotation=1,0,0,0,1,0,0,0,1`

Renders the current world as a PNG, base64-encoded. All drawing (elevation fill, plate-color
fill, boundary outlines, pole markers, rotation arcs, per-node dots) happens server-side
(see [simulation-model/rendering-and-inspectors.md#render-image](simulation-model/rendering-and-inspectors.md#render-image-render_imagepy)) -- the client
just decodes and paints the image, it never sees raw coordinate data. `400` for an
unrecognized projection/view name, a width/height outside `[1, main.MAX_RENDER_DIMENSION_PX]`
(4000), or a malformed `rotation`, `404` if no world has been generated yet. `globe` draws the
hemisphere facing the viewer as an orthographic disc (see
[simulation-model/rendering-and-inspectors.md#projections](simulation-model/rendering-and-inspectors.md#projections-projectionspy)).

```json
{
  "projection": "behrmann",
  "elapsed_years": 2000000,
  "image_base64": "iVBORw0KGgoAAAANSUhEUgAA..."
}
```

- `view` selects what gets drawn: `"elevation"` (colored by height/depth), `"plates"`
  (colored by owning plate, plus boundary outlines/pole markers/rotation arcs), or
  `"platesDetail"` (each plate's raw terrain nodes as dots, colored by elevation,
  plus boundary outlines) -- the frontend's Map View dropdown picks this directly.
  `"speckle"` is a debug overlay for diagnosing dithered coastlines (see
  [debugging/map-views.md#speckle-coastal-dither-overlay](debugging/map-views.md#speckle-coastal-dither-overlay)):
  a muted land/ocean backdrop with every node within 120 m of sea level drawn as a dot
  colored by its *coastal-dither fraction* -- the share of its nearest neighbours on the
  opposite side of the waterline -- and an oversized magenta marker on nodes at/above `0.75`
  (a near-isolated land speck or pond). `"biome"`
  and `"combined"` (Köppen-colored land, pelagic-province + hypsometric ocean,
  lakes/glaciers/rivers overlaid) are a categorical classification -- 31 Köppen-Geiger land
  classes + 10 pelagic ocean classes (see
  [simulation-model/climate.md#biomes](simulation-model/climate.md#biomes-biomespy)). `"combined"` alone returns an
  **RGBA** PNG whose alpha channel encodes a per-pixel class id for the client's
  click-to-identify feature: `code = 255 - alpha`, where `code` is `0` for a gap between
  cells, `biome_index + 1` (`1`..`41`, land Köppen then ocean pelagic) for a classified cell,
  `42` for a lake, or `43` for glacier cover. Alpha is a data channel, not real transparency
  (the range is only ~`212`..`255`). Five more views come from
  `climate.py` (see [simulation-model/climate.md#climate](simulation-model/climate.md#climate-climatepy)):
  `"temperature"`, `"humidity"`, and `"precipitation"` are heatmaps (each also drawing the
  current coastline -- see [simulation-model/hydrology.md#coastline](simulation-model/hydrology.md#coastline-coastlinepy) --
  since a color-scale view carries no land/ocean cue on its own); `"wind"` and
  `"oceanCurrents"` draw subsampled direction/magnitude arrows, and `"oceanCurrents"`
  additionally marks detected ocean swells with small circles. `"wind"` and the land side of
  `"temperature"` draw off `World.atmosphere_cfd_state` (a real, continuously time-integrated
  shallow-water wind solve, see
  [simulation-model/climate.md#ocean-atmospheric-fluid-dynamics](simulation-model/climate.md#atmospheric-wind-solver-fluid_dynamicspy-atmosphere_cfdpy))
  -- or, under `wind_model == "diagnostic"` (see
  [simulation-model/climate.md#wind-model](simulation-model/climate.md#wind-model-diagnostic-default-vs-cfd)), off `climate.py`'s own
  closed-form wind/air-temperature instead; everything else is a `climate.py` diagnostic (the
  shallow-water ocean solver was retired) resampled nearest-cell onto the world's HEALPix grid. `"resources"` and `"soilQuality"` (see
  [simulation-model/surface-processes.md#resources-and-soil](simulation-model/surface-processes.md#resources-and-soil-geologypy-plus-volcanismpys-own-eruption-roll)) are
  node-cloud-derived like elevation/plates, not climate-grid-derived -- `"resources"` overlays
  coal/oil & gas/mineral deposit richness on a muted land/ocean backdrop, `"soilQuality"` is a
  continuous fertility heatmap (barren to rich) plus the coastline overlay. `"geomorph"` (a
  debug view, see [debugging/map-views.md](debugging/map-views.md#erosion--deposition-geomorph)) colours every node by its net elevation change
  over the last step (`erosion.ErosionResult.net_elevation_change_m` off `World.erosion_cache`
  -- erosion minus every deposition pathway, no tectonics) on a diverging warm/cool scale plus
  the coastline overlay, so the per-step lumpiness of near-sea-level deposition is legible; a
  flat neutral field before the world has been stepped (or on a freshly loaded save, which
  doesn't persist the cache). `"elevReason"` (a debug view, see
  [debugging/map-views.md](debugging/map-views.md#elevreason-render-view-last-elevation-change)) colours every
  node by which process last moved its elevation -- a persistent per-node code
  (the `elev_change_reason` node field), so unlike `"geomorph"` it accumulates over the whole
  run: grey where terrain is untouched since generation, warm where crust is being built,
  cool/pale where it's being planed down or buried. `"overlapAge"` (a debug view, see
  [debugging/map-views.md](debugging/map-views.md#overlapage-render-view-plate-overlap-onset)) colours every node
  currently sitting on top of another plate by how long it has
  (`elapsed_years - overlap_onset_years`) over a muted land/ocean backdrop plus
  the coastline, pale where fresh and magenta where stuck for tens of Myr; all-backdrop is
  the healthy case. `"oceanCfdSediment"`/
  `"oceanCfdDeposition"` (from the retired ocean solver) are not valid `view` values
  (`/worlds/{world_id}/render` rejects them with `400`, same as any other unrecognized view name).
  `"cratonLayer"` / `"waterLayer"` (`render_image.LAYER_VIEWS`) aren't maps of their own:
  they return a transparent RGBA PNG painted only where there are cratons, or open ocean and
  visible lakes/seas, for the client-drawn "Plates & Faults" view to composite under its
  vectors. `/worlds/{world_id}/render` accepts them and `/worlds/{world_id}/animate` doesn't.
- `rotation` is the map's current view orientation (see
  [simulation-model/rendering-and-inspectors.md#rotating-the-view](simulation-model/rendering-and-inspectors.md#rotating-the-view)): a
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

## `GET /worlds/{world_id}/cross-section?plane_angle_deg=0&position_index=0&width=1100&height=611&show_boundaries=false&show_faults=false`

Renders a surface-sampled shallow mantle section for the current world. `plane_angle_deg`
selects a great-circle plane in `[-90, 90]`; `0` is the equator and either endpoint is a
pole-crossing plane. `position_index` selects one of 120 evenly spaced centers along that
circle in `[0, 119]`. `show_boundaries` and `show_faults` optionally draw surface markers
where the arc crosses plate boundaries and fault traces. `width` and `height` use the same `[1,
main.MAX_RENDER_DIMENSION_PX]` limit as `/worlds/{world_id}/render`. `404` means there is no current world;
`409` means the world has no active surface nodes. The JSON response contains `image_base64`
(PNG), `world_revision`, `elapsed_years`, the request's plane and position, `arc_length_km`,
`section_depth_km`, and `vertical_exaggeration`.

The PNG samples existing per-column surface fields along the selected 3,000 km arc at a rate
derived from `World.node_density`. Categorical crust type is nearest-node sampled, with no
cross-boundary blending. Crust and mantle-lithosphere thickness are drawn from their stored
columns; the crustal temperature profile is illustrative and melt eligibility is hatched
separately. The mantle below the stored lithosphere is schematic. Pending magma parcel origins
and stored volcano vents are shown at the surface. Tangential mantle-flow directions are
sampled from the surface field and drawn midway through the schematic asthenosphere to
distinguish them from surface features; this is not depth-resolved mantle flow or a subsurface
parcel trajectory. The color legend is rendered by the client. Section zoom is a client-side
view control and does not change the sampling footprint or server request size.

## `POST /worlds/{world_id}/animate`

The toolbar's "⏺ Record" action (opens the animation dialog, then starts this) -- an
H.264/MP4 video of `view`/`projection`'s progress, one frame for the world's current state
plus more, each `years_per_frame` further along:

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

Same fields/validation as `/worlds/{world_id}/render`'s own query params (`rotation` is the same
9-comma-separated-float string, `400` for an unrecognized `projection`/`view` or an
out-of-range `width`/`height`), plus `num_frames` bounded to `[1,
main.MAX_ANIMATION_FRAMES]` (480) -- each frame costs a full simulation step and render, so
this bounds worst-case request time. `step_years` gets the same validation as
`/worlds/{world_id}/step`'s `years` (`400` unless positive, finite, and at most 100 Myr, #557), and
`steps_per_frame` is bounded to `[1, main.MAX_STEPS_PER_FRAME]` (100, the dialog's largest
choice), so one frame holds the compute gate for at most 100 steps and a whole run for
at most the dialog's own maximum of 47,900. The animation dialog's "keep going until Stop is
pressed" mode just sends `num_frames` at that ceiling and expects the client to call
`POST /worlds/{world_id}/animate/stop` (below) rather than waiting for it to be reached on its own.
`404` if no world has been generated yet, `503` if a step or another animation is already
in progress.

**This permanently advances the world** by up to `(num_frames - 1) * years_per_frame` years
-- less if the run is ended early via `/worlds/{world_id}/animate/stop` -- the same as calling
`/worlds/{world_id}/step` that many times in a row -- deliberately not a side-effect-free preview (see
`render_image.stream_animation_mp4`'s own docstring).

The response is a stream of newline-delimited JSON objects (`application/x-ndjson`), one
per line, so the client can drive a progress bar instead of blocking on one opaque request:

```json
{"type": "progress", "frame": 1, "total": 20, "image_base64": "iVBORw0KGgo..."}
{"type": "progress", "frame": 2, "total": 20, "image_base64": "iVBORw0KGgo..."}
...
{"type": "done", "mime": "video/mp4", "video_base64": "AAAAIGZ0eXBpc29t...", "stopped_early": false, "seed": 9, "elapsed_years": 19000000, "num_plates": 6, "events": []}
```

While a frame's step batch runs, `{"type": "heartbeat", "frame": n, "total": N,
"elapsed_years": ...}` ticks flow whenever a step lands more than a minute after the
previous line (#611) -- imageless proof the run is still stepping, so the client's
idle watchdog doesn't abort a slow batch as stalled. The frontend ignores their shape;
their bytes alone reset its timer.

Each `progress` line carries that frame's own rendered PNG in `image_base64` (same encoding
as `/worlds/{world_id}/render`) -- the frontend paints it straight onto the main map so the display
doubles as the live animation step while the run holds the world lock. The final `done` line
carries the same summary shape `/worlds/{world_id}/generate`/`/worlds/{world_id}/step` return, plus `video_base64`
(decode via `data:video/mp4;base64,<this>` -- H.264, playable directly in a `<video>`
element), `mime`, and `stopped_early` (true if `/worlds/{world_id}/animate/stop` was what ended the run
short of `num_frames`, rather than reaching it on its own). If rendering raises partway
through, a `{"type": "error", "detail": "..."}` line is emitted instead -- the HTTP status
is already `200` by then.

## `POST /worlds/{world_id}/animate/stop`

The toolbar's "⏹ Stop" action while recording -- signals the in-progress `/worlds/{world_id}/animate`
stream to end cleanly after its current frame rather than the next one: the stream still
flushes the encoder and its `done` line still carries a complete (just shorter) video,
unlike the client simply dropping the connection. No request body; always `200`, including
when no animation is currently running (a no-op).

## `GET /worlds/{world_id}/plates`

The "Plate Inspector" map mode's data source (see
[simulation-model/rendering-and-inspectors.md#plate-inspector](simulation-model/rendering-and-inspectors.md#plate-inspector)) -- unlike
`/worlds/{world_id}/render`, this returns plain JSON, not a baked PNG; `frontend/src/PlateInspector.tsx`
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
      "num_points": 3437,
      "speed_cm_per_yr": 4.12,
      "at_max_rate": false,
      "euler_pole": { "lat_deg": 22.0, "lon_deg": -138.0 },
      "age_years": 4300000.0,
      "median_elevation_m": 210.0,
      "submerged_fraction": 0.22,
      "overlaps": [{ "plate_id": 6, "fraction": 0.058, "since_years": 178000000.0 }],
      "collisions": [{ "plate_id": 6, "years": 30700000.0 }],
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

- `num_points` is the plate's node (quad cell) count.
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
  [simulation-model/rendering-and-inspectors.md#plate-inspector](simulation-model/rendering-and-inspectors.md#plate-inspector) for the fitting
  method) -- `diameter_a_km`/`diameter_b_km` are the major/minor diameters in real km,
  `outline` is ~72 sampled world-space points around its perimeter for drawing.
- `speed_cm_per_yr` is `|Plate.omega|` as a surface speed at the planet's radius;
  `at_max_rate` is `true` when the plate is pinned at `mantle.MAX_PLATE_RATE` (15 cm/yr).
  `euler_pole` is the rotation axis in lat/lon degrees (`null` for a stationary plate).
  `age_years` is the simulated time since the plate was created (generation, split, or
  merge). These come straight from the dynamic torque solve (`torque.py`).
- `median_elevation_m` / `submerged_fraction` summarise the plate's own node elevations
  against `world.sea_level_m`. A continental plate reading mostly submerged is a red flag for
  over-stretching (see [GitHub issue #119](https://github.com/adubey/mantle-bloom/issues/119)).
- `overlaps` lists the other plates this plate's territory currently sits on top of --
  `fraction` is the share of *this* plate's nodes whose centres lie inside `plate_id`'s
  territory (`plates.compute_node_overlap`; plates that merely share a boundary don't
  count). `since_years` is
  the earliest `elapsed_years` any of this plate's still-overlapping nodes first went over
  another plate (the `overlap_onset_years` node field, stamped by
  `merge_split.update_overlap_tracking`) -- `null` if the save predates the field or the
  overlap only appeared this step; not partner-specific. See the `overlapAge` debug render
  view. `collisions` surfaces
  `merge_split.update_collision_progress`'s sustained-collision timers (accumulated
  convergent years) for pairs involving this plate. Both are diagnostics for the long-run
  plate-geometry degradation tracked in [GitHub issue #119](https://github.com/adubey/mantle-bloom/issues/119).

## `GET /worlds/{world_id}/plate_at?lat_deg=0&lon_deg=0`

The Plate Inspector's click hit-test: which plate owns the node nearest `(lat_deg, lon_deg)`
-- both in the *true* (un-rotated) frame; the client unprojects its click through whatever
view rotation is active first (see `rotation.ts`'s `unproject`/`matTranspose`). `400` for
non-finite input, `404` if no world has been generated yet.

```json
{ "plate_id": 3 }
```

`plate_id` is `null` only if every plate is empty (shouldn't happen via the API, but a
freshly-constructed `World` with no plates at all has no node to be nearest to).

## `GET /worlds/{world_id}/rivers`

The "River Inspector" map mode's data source (see
[simulation-model/hydrology.md#river-inspector](simulation-model/hydrology.md#river-inspector)) -- same "plain
JSON, client renders it" contract as `/worlds/{world_id}/plates`; `frontend/src/RiverInspector.tsx`
renders and drives the interaction itself. Un-rotated/true-frame throughout, same as
`/worlds/{world_id}/plates`. Regrouped fresh from `world.hydrology_cache` on every call (see
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
      "length_m": 2625000.0,
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

- `length_m` is the summed geodesic length of `segments`, in meters. River visibility
  uses this physical measure instead of a node-count cutoff, which would move with mesh density.
- `segments` is a flat list of line segments, normally between neighbouring nodes
  (`hydrology.flow_segments`), each a pair of world-space points -- not an ordered polyline,
  since a river network can branch. An ordinary member contributes one segment, to its own
  `flow_target`. A spilling lake sink whose `flow_target` lies beyond the rim contributes its
  spill path instead: several segments, stopping at visible lake water (#520; see
  [simulation-model/hydrology.md](simulation-model/hydrology.md#lake-vs-sea-tiers-a-depth-ceiling-by-size-2026-09-10)). Every
  coordinate in this response is rounded to 6 decimal places before serializing, same as
  `/worlds/{world_id}/plates`.
- `mouth_xyz` is the network's downstream-most member, the location the River Inspector marks
  for the selected river.
- `mouth_type` is `"ocean"` or `"lake"` (based on what the mouth's own `flow_target` empties
  into, not the mouth itself, and `"lake"` for a sink whose spill path reaches visible lake
  water -- see
  [simulation-model/hydrology.md#river-inspector](simulation-model/hydrology.md#river-inspector)), or `"other"`
  (a true dead-end sink, `flow_target == -1` -- the prognostic lake solver integrates its
  stored volume each step, so this flips to `"lake"` once the basin's water crosses
  `LAKE_MIN_VISIBLE_DEPTH_M`).
- `flow_rate` is physical discharge at the river mouth in m^3/s, from
  `HydrologyFields.discharge_m3_s`. `speed` is the stylized, unitless
  `hydrology.compute_river_speed` value evaluated there.
- `num_tributaries` counts separate headwater branches feeding the network -- see
  [simulation-model/hydrology.md#river-inspector](simulation-model/hydrology.md#river-inspector) for the exact
  definition (an original one; no precedent to port).
- `coastline_segments` is the land/ocean/lake boundary (see
  [simulation-model/hydrology.md#coastline](simulation-model/hydrology.md#coastline-coastlinepy)), same segment-pair shape as
  `segments` above -- included because this view draws no filled backdrop at all, so without
  it there's no land/ocean/lake cue in the River Inspector whatsoever. The exact same segments
  are also drawn server-side into the `temperature`/`humidity`/`precipitation` views returned
  by `GET /worlds/{world_id}/render` (no separate endpoint or query param -- it's baked directly into
  those PNGs, same as rivers/lakes are baked into the elevation view).

## `GET /worlds/{world_id}/river_at?lat_deg=0&lon_deg=0`

The River Inspector's click hit-test: which river network owns the node nearest
`(lat_deg, lon_deg)` -- same true-frame contract as `/worlds/{world_id}/plate_at`. `river_id` is an index
into the most recent `/worlds/{world_id}/rivers` response, not a persistent id. `400` for non-finite
input, `404` if no world has been generated yet.

```json
{ "river_id": 2 }
```

`river_id` is `null` if there are no rivers at all this step (before the first step, or a
world with no land steep/wet enough to route any).

## `GET /worlds/{world_id}/lakes`

The "Lake Inspector" map mode's data source (see
[simulation-model/hydrology.md#lake-inspector](simulation-model/hydrology.md#lake-inspector)) -- same "plain
JSON, client renders it" contract as `/worlds/{world_id}/rivers`; `frontend/src/LakeInspector.tsx` renders
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
      "inflow_rivers": [{ "mouth_xyz": [0.36, -0.10, -0.93], "flow_rate": 812.4, "num_nodes": 9, "length_m": 562500.0 }],
      "outflow_river": {
        "segments": [[[0.37, -0.09, -0.93], [0.38, -0.08, -0.92]], ["..."]],
        "mouth_xyz": [0.44, 0.02, -0.90],
        "mouth_type": "ocean",
        "flow_rate": 4108.2,
        "num_nodes": 14,
        "length_m": 875000.0
      }
    }
  ],
  "coastline_segments": [[[0.91, 0.28, -0.04], [0.90, 0.29, -0.03]], ["..."]]
}
```

- `member_xyz` is every currently-*flooded* node's own position (not the whole geometric
  catchment, which routinely includes dry higher ground on the way down to the basin's floor --
  see `main._lake_basin_summary`'s own comment) -- what the Lake Inspector plots per node,
  bright for the selected lake and dim for every other one.
- `floor_elevation_m` is this lake's lowest ground, counting each cell's sub-cell relief (the
  level at which the lake model holds no water), so it is never above `water_elevation_m`.
  `floor_xyz` is the lowest cell's centre. A lake wet only in its floor cell's low ground has
  its surface below that cell's mean elevation (#344).
- `outlet_xyz`/`outlet_elevation_m` is "the lowest point of the edge of the basin" -- the
  saddle a river out of it would source from, resolved from this step's freshly-rebuilt
  depression hierarchy (`lakes.build_lake_hierarchy`) so a lake that's currently the result of
  several smaller basins merging together reports its own *current* rim, not an interior saddle
  that's already flooded. `null` for an unresolved closed/endorheic basin with no known spill
  (a legitimate state, not missing data -- see `lakes.py`'s own docstring).
- `water_elevation_m`/`is_spilling` describe this lake's final same-step water surface and
  spill state from the prognostic hierarchy solve (`lakes.solve_lakes_prognostic`, #588). A
  capped lake reports spilling when surplus leaves its body even if its depth cap sits below
  the terrain rim.
- `inflow_rivers` is every `RiverInfo` (see `/worlds/{world_id}/rivers`) whose own mouth lands somewhere in
  this basin, as a small inline summary rather than a full `RiverSummary` -- the River
  Inspector's own map mode is what draws a river's full path. Each `flow_rate` is physical
  discharge in m^3/s; `length_m` is its summed geodesic flow-edge length.
- `outflow_river` is, unlike an inflow, found by which river network (if any) owns the node
  just *outside* the rim on the far side (`Lake.outlet_target_idx`) -- the same full shape
  `/worlds/{world_id}/rivers` returns (minus `river_id`/`speed`/`num_tributaries`, meaningless outside that
  view's own ranked list), since the Lake Inspector draws this one's whole path. Its
  `flow_rate` is physical discharge in m^3/s. `null` when
  there's no live outflow to show: a closed/endorheic basin with no known spill, a resolved rim
  that isn't currently `is_spilling`, or a spilling rim whose flow hasn't yet cleared
  `hydrology.RIVER_FLOW_PERCENTILE` to register as its own river network yet.
- `coastline_segments` is the same land/lake-vs-ocean boundary `/worlds/{world_id}/rivers` returns, included
  here for the same reason: this view draws no filled backdrop of its own.

## `GET /worlds/{world_id}/lake_at?lat_deg=0&lon_deg=0`

The Lake Inspector's click hit-test: which *basin* owns the node nearest `(lat_deg, lon_deg)`
-- same true-frame contract as `/worlds/{world_id}/river_at`, except it always resolves to something
informative for any land click, not just a hit on a currently-visible lake, so "click a basin
that has no water in it right now" still returns real information (the user's own spec: "give
me information about the basin nonetheless"). `400` for non-finite input, `404` if no world has
been generated yet.

```json
{ "kind": "lake", "basin": { "lake_id": 0, "is_lake": true, "...": "same shape as /worlds/{world_id}/lakes" } }
```

`kind` is one of:

- `"lake"` -- the nearest node is currently flooded; `basin` is that lake's own info (`lake_id`
  set, matching the same node's `/worlds/{world_id}/lakes` entry for this same step).
- `"basin"` -- dry land, but part of a real (possibly still-unresolved/endorheic) basin;
  `basin` is that basin's own leaf catchment info, same shape as a lake entry but `is_lake:
  false`, `lake_id: null`, `member_xyz` covering the whole geometric catchment (there's no
  flooded subset to narrow it to), and `water_elevation_m`/`is_spilling` both meaningless
  (`null`/`false`).
- `"no_basin"` -- dry land whose own steepest-descent chain drains straight to the ocean
  without ever passing through a local minimum first -- an ordinary hillslope, never part of
  any basin. `basin` is `null`.
- `"ocean"` -- the nearest node is open ocean. `basin` is `null`.

## `GET /worlds/{world_id}/faults`

One of the "Plates & Faults" map mode's data sources (see
[simulation-model/rendering-and-inspectors.md#fault-inspector](simulation-model/rendering-and-inspectors.md#fault-inspector)) -- same "plain
JSON, client renders it" contract as `/worlds/{world_id}/rivers`; `frontend/src/PlatesAndFaults.tsx`
renders the fault traces + systems over the plates (display-only there -- plate selection is
what interaction drives). These are **intraplate** faults (see `faults.py` /
[simulation-model/tectonics.md#faults](simulation-model/tectonics.md#faults)), distinct from plate boundaries.
Unlike `river_id` / `lake_id`, `fault_id` is a **persistent** identity for the fault's whole
lifetime (a monotonic `World.next_fault_id` counter), so a selection survives a step. `faults`
is `[]` before the first step. `404` if no world has been generated yet.

```json
{
  "elapsed_years": 12000000,
  "faults": [
    {
      "fault_id": 3,
      "plate_id": 5,
      "kind": "reverse",
      "trace": [[0.35, -0.11, -0.93], ["..."]],
      "active": true,
      "slip_rate_m_per_myr": 4200.0,
      "cumulative_offset_m": 18500.0,
      "age_myr": 4.4,
      "lifespan_myr": 12.0,
      "dip_deg": 31.0,
      "length_km": 63.0,
      "set_id": null,
      "system_id": 2,
      "birth_distance_from_boundary_km": 210.0,
      "distance_from_boundary_km": 260.0
    }
  ],
  "fault_systems": [
    {
      "system_id": 2,
      "plate_id": 5,
      "kind": "reverse",
      "trace": [[0.34, -0.10, -0.93], ["..."]],
      "active": true,
      "length_km": 2840.0,
      "age_myr": 4.4,
      "lifespan_myr": 96.0
    }
  ]
}
```

- `kind` is `"normal"`, `"reverse"`, or `"strike_slip"` (Andersonian regime, picked from the
  local closing rate at spawn).
- `trace` is the fault's polyline in **true** (un-rotated) world coordinates, an ordered list
  of unit vectors (rounded), not a branch-capable edge list.
- `active` is `false` once the fault has locked up and survives only as an inert scar.
- `cumulative_offset_m` is total slip since birth; `set_id` is non-`null` (and equal to the
  family's first member's `fault_id`) for a member of a tight sub-parallel fault family.
- `system_id` is non-`null` for a strand of a **fault system** -- the entry in `fault_systems`
  with the matching id holds that zone's long, gently curving master lineament (`trace`) and
  belt-scale metadata (see [simulation-model/tectonics.md#faults](simulation-model/tectonics.md#faults)). The
  master trace applies no relief of its own; its strands (the `faults` entries carrying the
  `system_id`) do. `fault_systems` is `[]` before the first step.
- `distance_from_boundary_km` is recomputed live each call (nearest cross-plate node), so it
  changes as the plate moves, unlike the stored `birth_distance_from_boundary_km`.

## `GET /worlds/{world_id}/earthquakes`

Recent earthquakes (see `faults.Earthquake` / [simulation-model/tectonics.md#faults](simulation-model/tectonics.md#faults))
for the fading epicentre overlay in the "Plates & Faults" view. The backend prunes anything
older than `faults.EARTHQUAKE_RETAIN_MYR` (~5 Myr) every step, so this list only ever holds
the last few Myr of activity. `[]` before the first step. `404` if no world has been generated yet.

## `GET /worlds/{world_id}/volcanoes`

Every current volcano node (see `volcanism.py` / `elevation_lines.is_volcano`) for the
"Plates & Faults" view's activity overlay -- same "plain JSON, client renders it" contract.
Un-rotated/true-frame world positions; `active` is `true` while the vent still has eruption
potential (`volcano_active_years_remaining > 0`) vs. a dormant cone. `[]` on a world with no
volcanoes yet. `404` if no world has been generated yet.

```json
{
  "elapsed_years": 12000000,
  "volcanoes": [
    { "position": [-0.73, 0.15, -0.67], "active": false }
  ]
}
```

```json
{
  "elapsed_years": 12000000,
  "earthquakes": [
    {
      "earthquake_id": 41,
      "fault_id": 3,
      "plate_id": 5,
      "kind": "reverse",
      "epicenter": [0.35, -0.11, -0.93],
      "magnitude": 7.12,
      "age_myr": 0.0,
      "birth_years": 12000000
    }
  ]
}
```

- `epicenter` is a single unit vector in **true** (un-rotated) world coordinates.
- `magnitude` is a moment magnitude (trace length + slip rate + an overlap-born bonus + a
  small per-event jitter), clamped to roughly `[4, 9]`.
- `age_myr` is `(elapsed_years - birth_years) / 1e6` -- the overlay fades a marker to zero by
  `EARTHQUAKE_RETAIN_MYR`.

## `GET /worlds/{world_id}/fault_at?lat_deg=0&lon_deg=0`

The Fault Line Inspector's click hit-test: the `fault_id` of the trace nearest
`(lat_deg, lon_deg)` within ~120 km, else `null` -- same true-frame contract as
`/worlds/{world_id}/river_at`. `400` for non-finite input, `404` if no world has been generated yet.

```json
{ "fault_id": 3 }
```

## `GET /worlds/{world_id}/stranded_basins`

Endorheic depressions whose floor sits below sea level and that have **no drainage path to
the ocean at all** -- the "land-locked coastal pit" pathology from
[GitHub issue #122](https://github.com/adubey/mantle-bloom/issues/122)'s coastal-speckle section, which the event log records but drowns in
near-sea-level transient-pond churn. Read straight off this step's already-resolved
depression hierarchy (`world.hydrology_cache.lake_forest`): a top-level basin with no known
spill (`lakes.Lake.max_depth is None`) whose floor is below `world.sea_level_m`. `stranded_basins`
is `[]` before the first step (or when no basin is stranded -- the healthy case). `404` if no
world has been generated yet. See [debugging/offline-cli.md](debugging/offline-cli.md#stranded-basins-cli) for how to read this and the
matching offline dump (`python -m app.hydroclimate.stranded_basins`).

```json
{
  "elapsed_years": 85100000,
  "sea_level_m": 0.0,
  "stranded_basins": [
    {
      "floor_elevation_m": -1771.0,
      "depth_below_sea_level_m": 1771.0,
      "catchment_node_count": 435,
      "flooded_node_count": 412,
      "water_elevation_m": -1750.3,
      "centroid_xyz": [0.62, 0.55, -0.21],
      "centroid_lat_deg": -12.3,
      "centroid_lon_deg": 45.6,
      "floor_xyz": [0.61, 0.56, -0.22],
      "first_seen_years": 72700000,
      "persisted_years": 12400000,
      "steps_seen": 124
    }
  ]
}
```

- entries are sorted deepest-first (`depth_below_sea_level_m` descending).
- `catchment_node_count` is the full geometric catchment (every node draining into the pit);
  `flooded_node_count` is how many of those currently hold visible standing water.
- `water_elevation_m` is the current standing-water surface, or `null` if the basin is dry.
- `first_seen_years`/`persisted_years`/`steps_seen` come from `world.stranded_basin_tracks`,
  a cross-step tracker `world.step_world` reconciles by centroid proximity every hydrology
  step (at most one step stale, same tolerance as `hydrology_cache` itself). A basin new this
  step reports `persisted_years: 0`, `steps_seen: 1`.

## `POST /worlds/{world_id}/export_hexgrid`

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

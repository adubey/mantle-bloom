# Architecture

## Stack

- **Backend:** Python, FastAPI + uvicorn, numpy (rotation/vector math), scipy (`cKDTree`
  for initial plate assignment, boundary-adjacency, merge-detection, and plate-defragment
  queries, `scipy.cluster.vq.kmeans2` for split clustering,
  `scipy.sparse.csgraph.connected_components` for plate defragmentation), Pillow
  (server-side map rendering, see `render_image.py`), pytest.
- **Frontend:** React + TypeScript via Vite, plain HTML `<canvas>` — no mapping/charting
  library: data-backed raster maps decode a PNG the backend already rendered and draw it
  (`ctx.drawImage`), which a library would be more ceremony than the problem needs. Other
  client-side views are `rotation.ts` — a
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
  [simulation-model.md#fault-inspector](simulation-model.md#fault-inspector)). Those four
  share `InspectorCanvas.tsx` (the focusable canvas, rotate-drag, click-to-select and
  Tab/Shift+Tab cycling), and with `MapCanvas.tsx` they share `screenProjector.ts`, which
  projects true world points through the live view rotation. The shallow
  mantle cross-section preview (`MantleCrossSection.tsx`) also draws its locator and schematic
  shell slice client-side; it has no sampled world data yet. The detailed data-backed rendering
  remains server-side.

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
    render call rather than stored server-side. The shallow mantle cross-section has its own
    `GET /world/cross-section` request path for plane angle, arc position, and output size;
    `cross_section.py` samples the existing surface columns along a node-density-scaled corridor
    and renders a shallow PNG. It draws stored crust/lithosphere thicknesses, an illustrative
    crustal geotherm, and supported surface vents and magma-source markers. Surface-sampled
    mantle-flow directions are projected into the schematic asthenosphere; deeper mantle flow
    remains unresolved. The frontend renders the shared legend over either placeholder or PNG.

Time-stepping:
  POST /world/step  { years }
  → world.step_world(world, years): every plate refits its Euler pole and rotates
    (`Plate.shift`), then every plate reconciles its actual footprint against the
    sphere minus every other live plate's own territory (`Plate.deform` -- quad_tectonics.py --
    in a freshly randomized order each turn) -- collision/subduction uplift or trench elevation where a plate's
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
    margins once, at generation -- see simulation-model.md#bathymetry.) Boundary advance into
    open ground happens inside every `deform()` call; whatever it leaves uncovered is filled
    by `gaps.py` every step (see simulation-model.md#gap-filling)
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
    hydroclimate/climate.py's own module docstring)

When the "Plate Inspector" or "Plates & Faults" map view is active, the browser instead
(also on every generate/step, but never on a projection/rotation-only change) fetches:
  GET /world/plates
  → every plate's outline + metadata as JSON, not a PNG -- see api-reference.md and
    simulation-model.md#plate-inspector. Clicking a plate sends its unprojected true
    lat/lon to GET /world/plate_at, a nearest-node lookup answering "which plate is here."
  "Plates & Faults" also fetches GET /world/faults, GET /world/earthquakes and
    GET /world/volcanoes (fault traces + fault systems, recent earthquake epicentres, and
    current volcano vents) and draws them over the plates; sidebar toggles hide the
    earthquake + volcano overlay and the faults. Its Cratons and Ocean / lake overlays are
    GET /world/render's transparent "cratonLayer" / "waterLayer" PNGs, fetched only while
    enabled (and, unlike the JSON above, refetched on projection/rotation changes), drawn
    under the vectors with the water translucent so submerged cratons show through. Faults are display-only there -- plate selection is what
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

The frontend never holds simulation state. Data-backed raster views re-fetch `/world/render`
after every `generate`/`step` call (or a projection/map-view/rotation change) and draw the
returned PNG onto the canvas (`MapCanvas.tsx`). That component also owns the rotate-the-planet drag gesture -- a
long-press-and-drag interaction that previews a cheap wireframe graticule client-side (via
`rotation.ts`) while dragging, then commits the real rotation and lets the usual
`/world/render` re-fetch replace it with the actual detailed frame (see
simulation-model.md#rotating-the-view); that gesture itself was later extracted into
`rotationDrag.ts`'s `useRotationDrag` hook so `PlateInspector.tsx` could reuse it for a
completely different kind of content (translucent plate ellipses, click-to-select,
Tab/Shift+Tab). The event console is the one other exception: it just displays whatever
`events` the last `generate`/`step` response carried, which is already the full current log
-- see api-reference.md).

<a id="access-control"></a>
## Access control (`backend/app/auth.py`)

The full app is single-owner (#372). `auth.AuthMiddleware` sits inside the CORS middleware
and checks every request before routing:

- `GET /` from someone who isn't the admin redirects to `/view`, the logged-out page (the
  frontend's `main.tsx` renders `Viewer.tsx` there and `App.tsx` everywhere else).
- Every `/world/*` route answers `401` to anyone who isn't the admin, including direct API
  calls.
- Any non-GET request whose `Origin` is another site (or `Sec-Fetch-Site: cross-site` when
  there's no `Origin`) is refused with `403`. This is the CSRF guard. Same-site means the
  request's own `Host` or the origin of `MANTLE_BLOOM_PUBLIC_URL`, so a proxy that rewrites
  `Host` still works. The Vite dev origin counts as same-site only for local requests.

The guard matches paths the way the router does, with any `root_path` prefix stripped, so
serving the app under a prefix doesn't open `/<prefix>/world/*`.

A request is the admin when it carries a valid session cookie, or when it was made on this
machine. A session cookie (`mb_session`) is an HMAC-signed `{uid, login, exp}` that
`/auth/github/callback` sets only when GitHub's numeric user ID matches the configured owner.
It's `Secure`, `HttpOnly` and `SameSite=Lax`, lasts 14 days, and isn't stored server-side.
A request counts as made on this machine when the client is loopback, `Host` is loopback,
and there are no proxy headers (`Forwarded`, `Via`, any `X-Forwarded-*`, `X-Real-IP`,
`CF-Connecting-IP` and similar). That covers `bin/restart.sh` and the desktop build, so
neither needs login. The OAuth flow uses PKCE and a signed, short-lived `state` cookie.

Configuration is read from the environment on every request:

| Variable | Meaning |
|---|---|
| `MANTLE_BLOOM_GITHUB_CLIENT_ID`, `MANTLE_BLOOM_GITHUB_CLIENT_SECRET` | The GitHub OAuth app. Its callback URL is `<public URL>/auth/github/callback`. |
| `MANTLE_BLOOM_GITHUB_OWNER_ID` | The one GitHub user ID allowed in (`gh api user -q .id`). |
| `MANTLE_BLOOM_SESSION_SECRET` | At least 32 bytes, used to sign cookies. Rotating it signs everyone out. |
| `MANTLE_BLOOM_PUBLIC_URL` | Optional. The public `https://` URL. The OAuth callback is built from it, and its origin is trusted for writes. Needed behind a TLS-terminating or `Host`-rewriting proxy. |
| `MANTLE_BLOOM_LOCAL_NO_LOGIN` | `0` turns off the local no-login rule. Set it behind a same-host proxy that rewrites `Host` to loopback and adds no forwarding headers. |

Missing any of the first four leaves login unconfigured, and the guard fails closed: remote
requests never get admin access, and `/login` says login isn't configured.

<a id="public-viewer"></a>
### The public viewer (`backend/app/public_view.py`)

Logged-out visitors land on `/view` (#373), which shows the shared world in the Eckert IV
projection, in the biome or elevation view only, and can step it by a fixed 100 kyr. It
reaches the world only through `/public/world`, `/public/render` and `/public/step` (see
api-reference.md). Those routes accept no projection, rotation, step size or other setting,
so creating, saving, loading, animating and reconfiguring a world all stay behind the
`/world/*` guard above.

Each browser may trigger at most 30 step computations. The count is kept in an HMAC-signed
`mb_guest` cookie, so it can't be edited down, though clearing the cookie resets it (a usage
limit, not abuse protection). It's signed with `MANTLE_BLOOM_SESSION_SECRET`, or, when that's
unset, a per-process secret that a restart replaces. The admin is never charged.

Viewers share one world, so they share the work of following it too (#374). Every change to
the world (generate, load, step, animation, controls, a guest step) publishes a new revision
through `WorldFeed`, a counter behind a per-process epoch, together with a summary snapshot
taken under the world lock. Viewers long-poll `/public/world?after=<revision>` and re-render
when it changes, so one viewer's step updates every page without anyone else computing.
Concurrent guest steps against the same revision run the simulation once, and the others get
its result for free. `/public/render` caches each (view, size) per revision in `RenderCache`,
so the first viewer after a change renders it and the rest are served without the world lock.
Any route that changes the world must call `main._world_changed()` (or `_set_world`) before
it releases `_world_lock`, or viewers keep seeing, and the cache keeps serving, the old world.

Only a request that runs the simulation is charged. A step request names the world revision
(seed and elapsed years) that its viewer is looking at. If the world has already moved past
that revision, the request gets the current world for free. If another guest's step is
running, the request waits for it and then gets that result for free. At most one guest step
runs at a time, and a guest never queues behind the owner: guest steps get `503` while the
owner holds the world, and guest renders get `503` after waiting 10 s.

<a id="run-and-checkpoint"></a>
## Running and checkpointing (`backend/app/run_control.py`, `backend/app/checkpoint.py`)

Play and Stop in the full app set a run status on the server (`POST /world/run`, behind the
`/world/*` guard, #375). While it's running, a background thread (`WorldRunner`) steps the
world by the chosen "Years per step", pausing 0.4 s between steps so renders and viewers get
the world lock. It keeps going with no tab open. A step that raises stops the run, and
`GET /world/run` reports why; an automatic stop never overrides a newer Play or Stop. The
full app follows a running world by long-polling `GET /world/run?after=<revision>&run_version=<n>`
and refreshing whenever the revision moves. It keeps following after Stop until the server
reports no step under way, so the last step still reaches the map, and a run-state change
(Stop from another tab) wakes the poll through `WorldFeed.wake`. `POST
/world/step` and animations still work alongside it; the run just waits its turn for the lock.

When `MANTLE_BLOOM_CHECKPOINT_DIR` (or, failing that, Railway's `RAILWAY_VOLUME_MOUNT_PATH`)
names a directory, the world and the run status are checkpointed there, so a hosted instance
that sleeps or restarts comes back where it was. Unset, as in a local run, nothing is saved.

- `world.mbworld` is an ordinary save (see [save-compatibility.md](save-compatibility.md)),
  rewritten after every *completed* change: `_world_changed()` writes it under the world lock
  for generate, load, controls, and every whole step (owner, guest or run loop). An animation
  writes it after each frame, so one that ends badly keeps the frames it finished. A
  step that raises partway still publishes what it left to viewers, but isn't checkpointed,
  so a restart returns to the last whole step. Pickling a default world takes about 30 ms and
  writes about 60 MB.
- `run_status.json` holds `running` and `step_years`, rewritten on every Play/Stop and when a
  run stops on its own. If that write fails, `POST /world/run` answers `500` and the status
  is written again with the next world checkpoint.
- Each file is written to a temporary file in the same directory, fsynced, and renamed over
  the old one, so a process killed mid-write leaves the previous checkpoint whole; leftover
  temporaries are deleted at startup. Nothing relies on a pre-sleep or shutdown hook.
- At startup (the app's lifespan) the world is restored exactly as saved, with no catch-up
  for the time spent asleep, and the run resumes if it was running. An unreadable world file
  is renamed to `world.mbworld.invalid-<time>` and the server starts with no world, stopped;
  an unreadable status file means stopped. Each case is logged. A directory that can't be
  created or written means no checkpointing at all, and `GET /healthz` answers `503`, so a
  hosted deploy with a broken volume fails its health check (#376).

Checkpoints assume one process owns the directory: keep the service at one replica while it
uses the volume.

The hosted instance (#376) is the root `Dockerfile` on Railway, configured by
`.railway/railway.ts`: one uvicorn process serving the built frontend and the API, with the
volume at `/data` and `GET /healthz` as the deploy health check. See
[HOSTING.md](HOSTING.md).

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
- `steps_taken` -- a plain `step_world` call count (not a year count -- step sizes vary and
  the thing it gates is about accumulated topology drift, not elapsed time). Drives the
  cadence of `merge_split.py`'s geometric plate defragmentation pass
  (`DEFRAG_INTERVAL_STEPS`, see
  [simulation-model.md#merge-and-split](simulation-model.md#merge-and-split)). A plain-int
  dataclass default is a class attribute, so worlds pickled before this field existed still
  load, reading 0.
- `collision_progress: dict[(int, int), float]` -- sustained-collision tracking for
  merge_split.py, pair of plate ids -> accumulated convergent (or craton-jammed, #401) years (see
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
- `stats_history: list[dict]` -- one `stats.compute_stats(self)` snapshot per real advance
  (generation, then every step), for the Stats panel's history charts. Deliberately
  uncapped, unlike `events` above. Appended to via `World.record_stats` (skipped, not forced,
  when climate wasn't computed this step -- see that method's own docstring), exposed via
  `GET /world/stats_history`; persisted in every saved `.mbworld` so a Load restores the
  charts instead of resetting them to empty, unlike the stateless `GET /world/stats` snapshot.
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
- the actual carried terrain: a sparse set of active cells on an equiangular cube-sphere
  lattice in the plate's own local frame (`PlateWithSparseQuadPatch`,
  `backend/app/sparse_quad_patch.py`), one terrain node per cell at its centre. This is the
  central data structure; see
  [simulation-model.md#plate-local-frames](simulation-model.md#plate-local-frames). Every
  per-node field in `surface_fields.SURFACE_FIELDS` -- `elevation`, crustal and
  mantle-lithosphere thickness, `channel_depth`/`channel_width`/`lake_depth`/`silt_depth`/
  `glacier_depth` (persistent, land-only -- see simulation-model.md#hydrology and
  simulation-model.md#glaciation), `is_volcano`/`volcano_active_years_remaining`
  (persistent -- see simulation-model.md#volcanism) and the provenance fields -- is one
  parallel array over the plate's cells, so it rotates with the plate for free, no advection
  scheme needed. Read fields with `Plate.collect` / `Plate.surface_nodes` and write them with
  `Plate.set_fields_on_plate`; a cell edit (insert, remove, refine, merge) carries every
  field through the cell's remap class, so a new field needs only its `SURFACE_FIELDS` entry.

A plate has no separately-tracked boundary polygon at all -- an earlier version kept one
(`boundary_local`, frozen at generation and rotated rigidly thereafter) purely for the
"Plates" map view's outline overlay, and it visibly drifted out of sync with the real
territory after enough stepping. `Plate.outline_world()` replaces it: the outline is traced
live from the plate's boundary cells, cached per topology/geometry revision, so it can never
be stale (see [simulation-model.md#boundary-evolution](simulation-model.md#boundary-evolution)).

## The simulation pipeline, module by module

```
geometry.py       unit-vector <-> lat/lon conversion, Rodrigues-formula rotation matrices,
                   plate-local coordinate frames, local tangent bases + azimuthal-equidistant
                   projection (for the Plate Inspector's bounding-ellipse fit)
ellipse.py         minimum-volume enclosing ellipse (Khachiyan's algorithm), sphere-agnostic
                   pure 2D math -- see simulation-model.md#plate-inspector
projections.py     Behrmann and Eckert IV map projections, vectorized
noise.py           cheap smooth sphere noise (sum of sinusoids) for initial terrain texture
elevation_lines.py  terrain-node constants and codes (planet radius, node density/spacing --
                    TARGET_LINE_SPACING_RAD, line_spacing_rad -- elevation-change and crust-type
                    codes) and the plate-local lattice sweep; the name predates #251, which
                    retired the line surface it used to hold
rtree_index.py      minimal bulk-loaded (STR-packed) R-tree over 2D points -- box/nearest-
                    neighbor queries, used by PlateWithRTree
plates.py          PlateSurface / Plate (ABC) -- node iteration and bulk field access, plate
                    identity and territory, `shift()` (per-turn Euler-pole refit + rotation),
                    the live per-plate outline used by the "Plates" map view, whole-world
                    overlap detection, and the Plate Inspector's bounding-ellipse fit and
                    nearest-plate click hit-test
sparse_quad_patch.py
                   PlateWithSparseQuadPatch, the cell-backed PlateSurface implementation
                   and the only one. Its row/column intervals are derived topology caches used
                   by cell-native algorithms; they do not expose row ends through the shared
                   surface API.
lithosphere_plate.py
                   plate generation (nearest-seed tiling, reference Hc/Hm, composite relief)
                   and the per-node column physics -- boundary classification context,
                   elevation/Hc/Hm deltas, decompression melting and rift volcanism
quad_tectonics.py  the tectonic engine, `deform()`: boundary retreat and advance on the cell
                   graph, rift opening, interior subduction, gap filling
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
gaps.py              whole-sphere coverage maintenance, every step (world.step_world): finds
                     any region no plate currently covers and grows the adjacent plate(s)
                     into it (`quad_tectonics.fill_gap`), or -- for ocean floor a
                     fully-subducted plate vacated with no neighbour left nearby -- spawns a
                     new plate to cover it (`lithosphere_plate.new_plate`). deform()'s own
                     boundary advance only grows a plate from its existing edge, so this is
                     the backstop for seams and vacated regions
hydroclimate/volcanism.py          every-step eruption lifecycle for existing volcano nodes (active-years
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
hydroclimate/        climate, hydrology, and surface-process package: atmospheric and
                     diagnostic climate, flow routing/lakes, erosion, geology, eustasy,
                     bathymetry, surface cover, volcanism, and their shared numerics
hydroclimate/climate.py           temperature/wind/currents/humidity/precipitation, computed fresh on
                     their own fixed equirectangular grid -- every render, and now every
                     step too, to drive erosion.py (see simulation-model.md#climate)
hydroclimate/erosion.py           every-step rain/river/weathering/glacier erosion + downstream sediment
                     deposition + glacier flattening + coastal planation/infill feedback
                     (near-sea-level wave-cut planation + sheltered-shelf silting, emergent
                     barrier islands), elevation deltas driven by climate.py's
                     current fields and hydrology.py's flow routing/glacier state (see
                     simulation-model.md#erosion) -- the weather-influences-geology half of
                     the coupling; climate.py's own elevation-reading mechanics (lapse rate,
                     mountain wind deflection, orographic rain shadow) are the other half
hydroclimate/hydrology.py         every-step flow routing over the geology node cloud (a k-nearest-
                     neighbor graph, not a grid): basin-spill detection, steepest-
                     descent flow direction, downstream flow accumulation, glacier
                     mass balance and shape (via ice_sheet.py) -- feeds erosion.py's river/glacier erosion and
                     deposition and the rendered river/lake/glacier overlay (see
                     simulation-model.md#hydrology and simulation-model.md#glaciation); also
                     groups the flat is_river mask into distinct drainage networks
                     (group_rivers) and answers the River Inspector's click hit-test
                     (river_at), on demand rather than every step (see
                     simulation-model.md#river-inspector)
hydroclimate/ice_sheet.py         every-step glacier solve, called from hydrology.compute_hydrology:
                     groups ice into bodies, conserves each body's volume with a backward-
                     Euler update of its extent, shapes it with a plastic ice-sheet profile
                     (an ordered-upwind numba sweep on the k-NN graph), books an ice ledger,
                     and diagnoses balance flux and the basal sliding speed erosion.py's
                     glacial abrasion reads (see simulation-model.md#glaciation)
hydroclimate/lakes.py              every-step lake growth/evaporation/merge/split/silt, an explicit n-ary
                     tree of Lake objects built from a depression-hierarchy pass over
                     hydroclimate/hydrology.py's own k-NN graph -- called from hydrology.compute_hydrology,
                     projects back down into the same flat lake_depth array every other
                     consumer already reads (see simulation-model.md#lakes-are-an-explicit-tree);
                     also rebuilt fresh, on demand, by main.py to answer the Lake Inspector's
                     GET /world/lakes and GET /world/lake_at (see
                     simulation-model.md#lake-inspector)
hydroclimate/stranded_basins.py   diagnostic-only: finds endorheic below-sea-level basins with no ocean
                     drainage (the "land-locked coastal pit") in lakes.py's forest, and
                     tracks how long each has persisted across steps (world.stranded_basin_
                     tracks, reconciled from world.step_world) -- backs GET /world/stranded_
                     basins and the python -m app.hydroclimate.stranded_basins offline dump (see
                     debugging.md)
hydroclimate/bathymetry.py        the shelf-width constant geology.py keys off, plus a one-off
                     generation-time pass (shape_initial_bathymetry) that drowns submerged
                     continental interiors toward abyssal depth by distance from land and
                     grades every continent/ocean plate margin into a slope, so the seabed
                     reads as real relief rather than a bright shelf and hard cliffs (see
                     simulation-model.md#bathymetry)
hydroclimate/coastline.py          traces the land/ocean and lake boundary as line segments over
                     hydroclimate/climate.py's own grid, on demand (not every step) -- drawn into the
                     temperature/humidity/precipitation renders and sent as JSON alongside
                     GET /world/rivers and GET /world/lakes (see simulation-model.md#coastline)
hydroclimate/fluid_dynamics.py      shared numerical primitives for atmosphere_cfd.py: real physical-unit
                     gradients/Laplacians/divergence, real Coriolis parameter, CFL-stable
                     substep sizing, semi-Lagrangian advection (see
                     simulation-model.md#ocean-atmospheric-fluid-dynamics)
hydroclimate/atmosphere_cfd.py      the atmospheric wind solver: a genuine time-integrated (not diagnostic)
                     shallow-water simulation of wind + air temperature -- prognostic state
                     (World.atmosphere_cfd_state) that persists/evolves step to step, plus a
                     sustained latitude-banded wind forcing. Ocean currents and precipitation
                     are diagnostic in climate.py, not CFD-solved -- the shallow-water ocean
                     solver was retired (see simulation-model.md#ocean-atmospheric-fluid-dynamics)
main.py              FastAPI routes
auth.py              single-owner GitHub login: the request guard, OAuth routes and signed
                     session cookies (see #access-control above)
run_control.py       the run/stopped status and the background thread that steps a running
                     world (see #run-and-checkpoint above)
checkpoint.py        atomic world and run-status checkpoints on persistent storage, and
                     restoring them at startup (see #run-and-checkpoint above)
public_view.py       the logged-out viewer's step allowance cookie, the coordinator
                     that shares one guest step among concurrent requests, and the
                     world revision feed and render cache viewers follow (see
                     #public-viewer above)
render_image.py      renders /world/render's requested view/resolution to a PNG server-side
                     (see simulation-model.md#render-image and simulation-model.md#climate),
                     plus /world/animate's streamed H.264/MP4 rendering (stream_animation_mp4)
                     and the "speckle" coastal-dither debug overlay (see docs/debugging.md)
persistence.py       whole-World save/load to a single versioned pickle envelope (File >
                     Save/Load World -- see api-reference.md's /world/save//world/load and
                     save-compatibility.md)
legacy_conversion.py the legacy reader for line-backed saves and their one-way conversion
                     to sparse quads on load (issues #248/#251, see save-compatibility.md)
geodesic.py          geodesic-icosahedron hex/pentagon dome tiling + elevation/biome
                     sampling for File > Export Hex Grid (see docs/hex-export-format.md),
                     independent of the plate simulation's own node cloud
```

See [simulation-model.md](simulation-model.md) for what each of these actually computes and
why.

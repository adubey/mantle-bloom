# Debugging & Diagnostic Views

There is no dedicated "debug build" or hidden flag. Instead, a handful of views, endpoints,
and one offline tool exist purely to answer *"is this world's geometry / climate / hydrology
healthy, or has a long run degraded it?"* -- they surface numbers the ordinary Elevation /
Biome / Climate renders don't. This page collects them.

Most of them were built while chasing the two long-run pathologies still tracked in
[TODO.md](TODO.md): plate geometry degrading over tens of My (pole winding, unbounded
overlap, over-stretched continental plates) and speckled low-relief coastlines. That
section's **"Diagnostic views & debug output"** heading lists what has landed and what is
still worth building.

---

## Plate diagnostics dump -- `python -m app.plate_diagnostics`

An offline, read-only dump of a saved world's plate geometry. It **never starts the server
or binds a port** -- it loads a `.mbworld` file directly (same pickle format as *File >
Load*, see [`persistence.py`](../backend/app/persistence.py)) and prints to stdout.

```bash
cd backend
source .venv/bin/activate
python -m app.plate_diagnostics ~/Downloads/mantle-bloom-seed888151728-85000000y.mbworld
python -m app.plate_diagnostics <save.mbworld> --json      # structured, for scripting
```

It reuses `main._plate_summary` / `main._plate_overlaps` -- the exact code path behind
`GET /world/plates` and the in-app Plate Inspector -- so the CLI and the UI can never
disagree.

### What it prints

```
mantle-bloom plate diagnostics
  seed:          505070493
  elapsed:       25,500,000 yr  (~255 steps @ 100 ky)
  node_density:  4.0
  sea level:     0.0 m
  plates:        18

per-plate  ( * = railed at MAX_PLATE_RATE 15.0 cm/yr;  ! = continental & >50% submerged )
   id  crust         nodes  rows   age     speed     pole lat,lon   med.elev  submrg
  ----------------------------------------------------------------------------------
    0  oceanic        3364   121   255     4.59     +21.7,   -6.9      -3997   0.99
    1  continental    6843    89   255     2.43     -24.1, +169.3       1391   0.08
    ...

territory overlaps  (share of THIS plate's nodes sitting on top of another; >= 0.5%, full list in --json)
    1 -> 11    22.3%
    1 -> 7      6.0%
    ...

sustained-collision timers  (world.collision_progress)
  (  1,  11)    25.5 My
  (  4,  13)    25.4 My
  ...

node budget
  total nodes:             137,261
  clean-tiling estimate:   130,577   (4*pi / line_spacing_rad(4.0)^2)
  ratio:                      1.05x   (+5%)
```

### How to read it

- **`speed` + the `*` flag.** `*` means the plate is pinned at `mantle.MAX_PLATE_RATE`
  (15 cm/yr). One or two railed plates is normal (genuine slab pull). *Most* plates railed --
  and especially *every oceanic* plate at exactly 15.0 -- is the pathology in TODO.md's
  plate-geometry item 1. (The stiff-basal-drag bug that pinned *all* plates from step 1 was
  fixed 2026-08-30; a mostly-ocean world can still rail its oceanic plates for real reasons.)
- **`med.elev` + `submrg` + the `!` flag.** `!` marks a **continental** plate with more than
  half its own nodes at/below sea level -- the signature of an over-stretched continental
  plate whose interior the bathymetry model has (correctly) oceanised. A healthy continental
  plate sits well above sea level with `submrg` well under 0.5; oceanic plates near 1.0 are
  expected.
- **`rows` vs `nodes`.** `rows` is `ElevationLine`s with at least one node. A plate with far
  more nodes-per-row than its neighbours (≫ ~150 at `node_density=4`) is a winding row --
  see plate-geometry item 3 / the "streaking" symptom.
- **territory overlaps.** `A -> B  X%` = X% of plate A's own nodes sit within half a target
  node spacing of a node owned by B (ordinary shared boundaries are ~one full spacing apart,
  so this only fires on genuine overlap). The text dump hides entries below 0.5%; `--json`
  has the full list. A stable double-digit overlap that is *not* also in the collision
  timers will never trigger the merge path -- TODO.md plate-geometry item 4.
- **sustained-collision timers** are `world.collision_progress` -- accumulated convergent
  years per plate pair (`merge_split.update_collision_progress`). Compared against the
  50--100 My merge threshold, these tell you which overlaps are on track to heal and which
  are stuck.
- **node budget.** `clean-tiling estimate` is how many nodes a gap-free, non-overlapping
  lattice would put on the whole sphere at this `node_density` (`4*pi / line_spacing_rad²`
  -- ~32.6k at 1x, ~130k at 4x). `ratio` is the world's actual total over that. Up to
  ~1.15x is the documented bounded-envelope / randomized-order effect; 1.5--1.75x is the
  long-run node blowup (plate-geometry item 5).

Tests: [`unit_tests/test_plate_diagnostics.py`](../backend/unit_tests/test_plate_diagnostics.py).

---

## Plate Inspector diagnostic fields (in-app)

The **Plate Inspector** map view (`frontend/src/PlateInspector.tsx`, fed by
`GET /world/plates`) plots every plate's nodes and outline client-side and, per plate,
reports the same motion / shape / overlap numbers the offline dump does:
`speed_cm_per_yr` + `at_max_rate` (shown red), `euler_pole`, `age_steps`,
`median_elevation_m` + `submerged_fraction` (red when a continental plate is >50%
submerged), `overlaps`, and `collisions`. Full field reference:
[api-reference.md](api-reference.md) (`GET /world/plates`); the panel layout and the
bounding-ellipse fit are in
[simulation-model.md#plate-inspector](simulation-model.md#plate-inspector).

Use the Inspector for the *visual* read -- concentric rings (pole winding), a plate's dots
sitting inside a neighbour's (overlap), a long straight sawtooth chord across open ocean (an
over-extended lattice) -- and the dump for the numbers behind it.

---

## `platesDetail` render view

`GET /world/render?view=platesDetail` draws each plate's raw `ElevationLine` nodes as dots
coloured by elevation (not the smoothed territory fill that `plates` uses). It is the
fastest way to *see* lattice-level damage -- winding rows, stray one-node "teeth", a
staircase plate edge -- against the actual elevation field. See
[api-reference.md](api-reference.md) (`GET /world/render`) for the full view list.

---

## `geomorph` render view (erosion & deposition rate)

`GET /world/render?view=geomorph` (Map View dropdown: **Debug > Erosion & Deposition**)
colours every node by its net elevation change over the last step --
`erosion.ErosionResult.net_elevation_change_m` (post-erosion elevation minus pre-erosion, so
erosion minus every deposition pathway plus the small flatten/lake-siltation terms; *not*
tectonic deform, isostasy, or volcanism), retained on `World.erosion_cache` purely for this
view. A diverging scale: warm brown/orange where the step net-lowered a node, cool blue where
it net-raised one, a flat neutral grey in the +-few-metre band so only the lumps stand out,
clamped past +-60 m/step. The coastline is overlaid for orientation.

What it's for: the per-step deposition in the near-sea-level band is wildly lumpy -- a
+200 m spike on one node, ~0 on its neighbour -- which is the mechanism behind the coastal
checkerboard (see [TODO.md](TODO.md), "Speckled low-relief coastlines"), but is invisible in
every other view. Step the world once with climate & biomes on, then switch to this view and
look along a drowned shelf: a clean coastal plain deposits smoothly (uniform pale colour), a
dithering one shows a salt-and-pepper mix of saturated warm and cool cells. Use it as a
before/after for any coastal-feedback change instead of an ad-hoc script.

`erosion_cache` is `None` until the first climate/erosion step (and on a freshly loaded save
-- it isn't persisted), where the view falls back to a flat neutral field plus the coastline
rather than erroring.

---

## River & Lake Inspectors

`GET /world/rivers` / `GET /world/lakes` and their map views
(`RiverInspector.tsx` / `LakeInspector.tsx`) render flow networks and lake basins from raw
JSON. They are primarily feature views, but the Lake Inspector is also where persistent
endorheic basins show up -- see
[simulation-model.md#river-inspector](simulation-model.md#river-inspector) and
[simulation-model.md#lake-inspector](simulation-model.md#lake-inspector).

---

## Event log

The `events` list on `GET /world/summary` (the UI's event console) logs lake
formation/splits and other discrete events.

### Lake-churn aggregation -- `lakes.summarize_lake_events`

On a long run over a dithering low-relief coast the lake solver produces hundreds of
near-sea-level transient merge/split transitions per My -- one pair per puddle per step (see
[TODO.md](TODO.md) "Speckled low-relief coastlines"). Left raw, these flood the console and
bury real basin/tectonic events.

`lakes.step_lakes` now returns structured `lakes.LakeEvent`s
(`kind` / `node_count` / `elevation_m` / `basin_count`, plus a `.message` property with the
same wording as before) instead of pre-formatted strings.
[`erosion.py`](../backend/app/erosion.py) runs a step's events through
`lakes.summarize_lake_events(events, world.sea_level_m)` before logging:

- A transition whose water surface is more than `lakes.NEAR_SEA_LEVEL_EVENT_BAND_M`
  (**15 m**) from the current sea level is a genuine basin event -- logged individually,
  unchanged.
- Transitions **within** that band are the coastal-pond churn. A lone one still logs
  verbatim; **two or more in one step collapse to a single line** --
  `"38 transient coastal ponds churned near sea level this step (22 merged, 16 split)."`

So a persistent deep endorheic basin (e.g. the real ~435-node lake oscillating near
-1770 m) stays visible in the console while the checkerboard shelf contributes at most one
aggregate line per step. The band is measured against `world.sea_level_m`, so it tracks a
sea-level control change. Tests:
[`unit_tests/test_lakes.py`](../backend/unit_tests/test_lakes.py)
(`test_summarize_lake_events_*`).

---

## Still worth building

From TODO.md's "Diagnostic views & debug output" section, not yet implemented:

1. **Speckle / coastal-dither overlay render mode** -- colour every near-sea-level node by
   the fraction of its neighbours on the opposite side of the waterline; flag ≥ 0.75.
   Instantly distinguishes a checkerboard coast from a clean shoreline.
2. **Per-node geomorph-rate view** -- render `ErosionResult.sediment_deposited` (or net
   `dElev`/step) as a diverging map. The lumpiness of near-sea-level deposition is invisible
   in every current view but is the whole coastal-speckle mechanism.
3. **Stranded-basin report** -- a `/world/lakes`-style listing of endorheic basins whose
   floor is below sea level and that aren't ocean-connected: node count, floor elevation,
   centroid, persistence.

(Event-log dedup for lake churn -- formerly item 4 -- landed 2026-08-31; see the
Lake-churn aggregation section above.)

See [TODO.md](TODO.md) for the full rationale on each.

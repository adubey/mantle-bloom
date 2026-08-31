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

## River & Lake Inspectors

`GET /world/rivers` / `GET /world/lakes` and their map views
(`RiverInspector.tsx` / `LakeInspector.tsx`) render flow networks and lake basins from raw
JSON. They are primarily feature views, but the Lake Inspector is also where persistent
endorheic basins show up -- see
[simulation-model.md#river-inspector](simulation-model.md#river-inspector) and
[simulation-model.md#lake-inspector](simulation-model.md#lake-inspector).

---

## Stranded-basin report -- `GET /world/stranded_basins` + `python -m app.stranded_basins`

A "stranded basin" is the **land-locked coastal pit** from TODO.md's coastal-speckle
section: an endorheic depression whose floor sits *below sea level* and that has **no
drainage path to the ocean at all**. Such a node is neither hydrology's connectivity-aware
`is_ocean` nor above sea level, so the marine sink, coastal planation, and lake infill all
skip it -- it churns (merge/split) in the event log every step and never drains or fills.
The event log has this today but drowns it in near-sea-level transient-pond spam (see below);
this surfaces the same thing as one clean list.

The criterion is read straight off this step's already-resolved depression hierarchy
(`hydrology.HydrologyFields.lake_forest`, `lakes.build_lake_hierarchy`): a **top-level** basin
whose `max_depth is None` (lakes.py's own "no known spill to the ocean" state) *and* whose
`floor_elevation` is below `world.sea_level_m`. Roots only -- an endorheic root is the maximal
"no drainage" unit and its floor is the min over every descendant, so a deep sub-basin is
already covered.

Both the endpoint and the offline dump go through `stranded_basins.find_stranded_basins` /
`enrich_with_persistence`, so they can't disagree. Persistence -- *how long* each pit has
been there -- comes from `world.stranded_basin_tracks`, a small cross-step tracker
`world.step_world` reconciles each hydrology step by matching this step's basins to last
step's by centroid proximity (the same lightweight first-seen-per-key idea
`world.collision_progress` uses for plate pairs; diagnostic only, nothing in the physics
reads it back). It's persisted in the save, so the offline dump reports real persistence
numbers as of save time.

```bash
cd backend
source .venv/bin/activate
python -m app.stranded_basins ~/Downloads/mantle-bloom-seed888151728-85000000y.mbworld
python -m app.stranded_basins <save.mbworld> --json
```

```
mantle-bloom stranded-basin diagnostics
  seed:          888151728
  elapsed:       85,100,000 yr  (~851 steps @ 100 ky)
  node_density:  4.0
  sea level:     0.0 m
  stranded basins: 2   (endorheic, floor below sea level, no ocean drainage)

     floor  depth<SL   catch  flooded    water   centroid lat,lon         persisted
  ------------------------------------------------------------------------------------
     -4560      4560     512      480    -4400     -31.4,   +88.7   18.2 My (182 steps)
     -1771      1771     435      412    -1750     -12.3,   +45.6   12.4 My (124 steps)
```

- **`floor` / `depth<SL`** -- basin floor elevation and how far below sea level that is.
- **`catch` / `flooded`** -- the full geometric catchment node count vs. how many members
  currently hold visible standing water.
- **`water`** -- current standing-water surface elevation (`--` if bone dry).
- **`persisted`** -- elapsed years (and approx 100-ky steps) since a basin first appeared at
  this centroid. A large number here is the signal: a pit that's been stranded for tens of My
  is a real drainage/infill gap, not a one-step transient.

An empty list is the healthy case -- most seeds never strand a basin. The report needs a
hydrology snapshot in the save (a world stepped at least once with climate on); a
never-stepped world reports nothing.

Test: [`unit_tests/test_stranded_basins.py`](../backend/unit_tests/test_stranded_basins.py).

---

## Event log

The `events` list on `GET /world/summary` (the UI's event console) logs lake
formation/splits and other discrete events. **Known limitation:** on a long run over a
dithering low-relief coast it floods with hundreds of near-sea-level transient
"N-node lake formed/split at ~0 m" lines per My, which buries real basin/tectonic events.
Dedup/severity for this is TODO.md's "Diagnostic views & debug output" item 4.

---

## Still worth building

From TODO.md's "Diagnostic views & debug output" section, not yet implemented:

1. **Speckle / coastal-dither overlay render mode** -- colour every near-sea-level node by
   the fraction of its neighbours on the opposite side of the waterline; flag ≥ 0.75.
   Instantly distinguishes a checkerboard coast from a clean shoreline.
2. **Per-node geomorph-rate view** -- render `ErosionResult.sediment_deposited` (or net
   `dElev`/step) as a diverging map. The lumpiness of near-sea-level deposition is invisible
   in every current view but is the whole coastal-speckle mechanism.

The **stranded-basin report** (was item 3) landed -- see the section above.

See [TODO.md](TODO.md) for the full rationale on each.

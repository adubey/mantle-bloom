# Debugging & Diagnostic Tools

Working through a degradation (plate geometry going bad on long runs, coastlines dithering
pixel-by-pixel) usually needs a number the ordinary map views don't surface. This page collects
what exists for that: a handful of **command-line tools** that run offline against a saved
world, and a cluster of **UI elements** -- extra map views, side panels, and a dialog tab --
that only appear in the frontend when the page is loaded with `?deb` in the URL.

Most of this was built while chasing two long-run pathologies: plate geometry degrading over
tens of My (pole winding, unbounded overlap, over-stretched continental plates) and speckled
low-relief coastlines (see [GitHub issue #122](https://github.com/adubey/mantle-bloom/issues/122)).

---

## Debugging command-line tools

Four offline, read-only dumps of a saved world. Each **never starts the server or binds a
port** -- three load a `.mbworld` file directly (same pickle format as *File > Load*, see
[`persistence.py`](../backend/app/persistence.py)) and print to stdout; `phase_budget_diagnostics`
can do that too, but (unlike the other three, which only ever read a save) it also *steps* the
world forward by a given number of years before reporting -- see its own section below. Each of
the first three reuses the exact code behind the equivalent API endpoint / UI panel, so the CLI
and the UI can never disagree.

```bash
cd backend
source .venv/bin/activate
python -m app.plate_diagnostics <save.mbworld>            # plate geometry
python -m app.stranded_basins <save.mbworld>               # land-locked sub-sea-level pits
python -m app.lake_hierarchy_diagnostics <save.mbworld>    # lake-merge-forest depth/size
python -m app.phase_budget_diagnostics <save.mbworld> --years 1000000  # per-phase Hc/Hm budget
```

Every tool also takes `--json` for structured, scriptable output.

### `python -m app.plate_diagnostics` -- plate geometry

Reuses `main._plate_summary` / `main._plate_overlaps` -- the exact code path behind
`GET /world/plates` and the in-app Plate Inspector.

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

**How to read it:**

- **`speed` + the `*` flag.** `*` means the plate is pinned at `mantle.MAX_PLATE_RATE`
  (15 cm/yr). One or two railed plates is normal (genuine slab pull); *most* plates railed --
  and especially *every oceanic* plate at exactly 15.0 -- is a sign the torque balance has
  broken down. (The stiff-basal-drag bug that pinned *all* plates from step 1 was fixed
  2026-08-30; a mostly-ocean world can still rail its oceanic plates for real reasons.)
- **`med.elev` + `submrg` + the `!` flag.** `!` marks a **continental** plate with more than
  half its own nodes at/below sea level -- the signature of an over-stretched continental
  plate whose interior the bathymetry model has (correctly) oceanised. A healthy continental
  plate sits well above sea level with `submrg` well under 0.5; oceanic plates near 1.0 are
  expected.
- **`rows` vs `nodes`.** `rows` is `ElevationLine`s with at least one node. A plate with far
  more nodes-per-row than its neighbours (≫ ~150 at `node_density=4`) is a winding row --
  the "streaking" symptom.
- **territory overlaps.** `A -> B  X%` = X% of plate A's own nodes sit within half a target
  node spacing of a node owned by B (ordinary shared boundaries are ~one full spacing apart,
  so this only fires on genuine overlap). The text dump hides entries below 0.5%; `--json`
  has the full list. A stable double-digit overlap that is *not* also in the collision timers
  will never trigger the merge path.
- **sustained-collision timers** are `world.collision_progress` -- accumulated convergent
  years per plate pair (`merge_split.update_collision_progress`). Compared against the
  50--100 My merge threshold, these tell you which overlaps are on track to heal and which
  are stuck.
- **node budget.** `clean-tiling estimate` is how many nodes a gap-free, non-overlapping
  lattice would put on the whole sphere at this `node_density` (`4*pi / line_spacing_rad²`
  -- ~32.6k at 1x, ~130k at 4x). `ratio` is the world's actual total over that. Up to
  ~1.15x is the documented bounded-envelope / randomized-order effect; 1.5--1.75x is a
  long-run node-count blowup worth investigating.

Tests: [`unit_tests/test_plate_diagnostics.py`](../backend/unit_tests/test_plate_diagnostics.py).

<a id="stranded-basins-cli"></a>
### `python -m app.stranded_basins` -- land-locked sub-sea-level pits

A "stranded basin" is an endorheic depression whose floor sits *below sea level* and that has
**no drainage path to the ocean at all**. Such a node is neither hydrology's connectivity-aware
`is_ocean` nor above sea level, so the marine sink, coastal planation, and lake infill all skip
it -- it churns (merge/split) in the event log every step and never drains or fills. The event
log has this today but drowns it in near-sea-level transient-pond spam (see
[Event log](#event-log) below); this surfaces the same thing as one clean list.

The criterion is read straight off this step's already-resolved depression hierarchy
(`hydrology.HydrologyFields.lake_forest`, `lakes.build_lake_hierarchy`): a **top-level** basin
whose `max_depth is None` (lakes.py's own "no known spill to the ocean" state) *and* whose
`floor_elevation` is below `world.sea_level_m`. Roots only -- an endorheic root is the maximal
"no drainage" unit and its floor is the min over every descendant, so a deep sub-basin is
already covered.

Both the CLI and `GET /world/stranded_basins` go through
`stranded_basins.find_stranded_basins` / `enrich_with_persistence`, so they can't disagree.
Persistence -- *how long* each pit has been there -- comes from `world.stranded_basin_tracks`,
a small cross-step tracker `world.step_world` reconciles each hydrology step by matching this
step's basins to last step's by centroid proximity (the same lightweight first-seen-per-key
idea `world.collision_progress` uses for plate pairs; diagnostic only, nothing in the physics
reads it back). It's persisted in the save, so the CLI reports real persistence numbers as of
save time.

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

### `python -m app.lake_hierarchy_diagnostics` -- lake-merge-forest depth/size

Measures two numbers nothing else does: the longest root-to-leaf chain in
`lakes.build_lake_hierarchy`'s merge forest, and a histogram of leaf-catchment node counts.
Worth checking before tuning `lakes.SILT_ACCUMULATION_COEFFICIENT`, or deciding whether a
depression pre-fill pass is warranted -- a deep merge forest (thousands of tiny
sub-resolution catchments each spilling into the next, rather than siltation collapsing them
into a handful of real basins) makes both of those guesswork otherwise.

Same shape as the other two dumps -- reads `world.hydrology_cache.lake_forest`, never starts
the server:

```
mantle-bloom lake-hierarchy diagnostics
  seed:          888151728
  elapsed:       5,000,000 yr  (~50 steps @ 100 ky)
  node_density:  4.0
  roots:         68
  leaf catchments: 82
  hierarchy depth: max 4   mean 1.21

leaf catchment size (node count)
         1-1: 0
         2-5: 6
        6-20: 15
       21-100: 39
      101-500: 20
        501+: 2

root-to-leaf hierarchy depth (chain length)
         1-1: 60
         2-5: 8
        6-20: 0
       21-100: 0
      101-500: 0
      501-1000: 0
       1001+: 0
```

- **`hierarchy depth`** -- `max` is the headline number; `mean` is over roots only (a forest
  of mostly-unmerged 1-level leaves still reports mean close to 1 even if one long cascade
  exists, so read `max` first).
- **Leaf catchment size vs. hierarchy depth histograms** -- leaf size is a statement about
  how many genuinely tiny sub-resolution depressions exist right now; depth is a statement
  about how long the spill *chains* between them run. A world could have many tiny leaves
  that all merge shallowly (low depth, e.g. a wide flat plain with lots of small independent
  pits), or few leaves chained very deep (a long river-like cascade of saddles) -- the two
  numbers answer different questions about the same pathology.

An empty/all-zero report is the healthy case -- most young or smooth worlds never build a
deep cascade. The report needs a hydrology snapshot in the save (a world stepped at least
once with climate on); a never-stepped world reports nothing, same convention as
`stranded_basins`.

**Issue #144 (closed).** Built for that issue's item #3; items #1 (silt coefficient's
magnitude vs. tectonic roughening) and #2 (whether a depression pre-fill pass belongs in the
pipeline) are answered together in `lakes.py`'s own comment above
`SILT_ACCUMULATION_COEFFICIENT` -- including a re-check against issue #117's exact repro seed
(23097282), regenerated from scratch and extended to 300 My, specifically to rule out issue
#161's since-fixed unbounded Hc/Hm growth (spurious deep, narrow wells that this module's
catchment detection could have mistaken for real depressions) as a confound in the original
measurement. Depth stayed bounded and self-resolving the whole run; no pre-fill pass is
warranted.

Test: [`unit_tests/test_lake_hierarchy_diagnostics.py`](../backend/unit_tests/test_lake_hierarchy_diagnostics.py).

### `python -m app.phase_budget_diagnostics` -- per-phase Hc/Hm budget

Built for [GitHub issue #216](https://github.com/adubey/mantle-bloom/issues/216)
("Investigate long-run Hc/Hm decline"). Every mechanism that touches crustal thickness (Hc) /
mantle-lithosphere thickness (Hm) -- convergent/divergent deformation, arc magmatism, oceanic
cooling relaxation, decompression melting, boundary growth/shrink, row claiming,
regularization, plate merges/cleanup/relatticing, failed rifts, erosion -- is individually
instrumented (see [`phase_budget.py`](../backend/app/phase_budget.py)'s module docstring for
the full mechanism-to-call-site mapping). This tool turns that instrumentation on, steps the
world forward, and reports the accumulated before/after node count and sum(Hc)/sum(Hm) per
phase, split by both plate type (`crust_type`) and per-node type (`crust_type_code`) -- so a
mean-thickness decline can be attributed to a specific mechanism, to changing node population,
or to reclassification, rather than treated as one undifferentiated trend.

Unlike the other three tools above, this one *replays* rather than only reads: it takes
`--years` (required) and steps the world forward that far, in ordinary 100 ky increments,
with `debug_diagnostics` turned on and the budget reset first -- so the report always
describes exactly the given interval, whether that's a short replay from a real save (the
issue's own "replay short intervals from this save to establish current contributors") or a
fresh world (`--seed`/`--node-density`, no save given, for a fully reproducible run):

```bash
python -m app.phase_budget_diagnostics <save.mbworld> --years 1000000
python -m app.phase_budget_diagnostics --seed 611937962 --node-density 4 --years 1000000
```

```
mantle-bloom Hc/Hm phase budget (GitHub issue #216)
  seed:          611937962
  node_density:  4.0
  interval:      334,100,000 -> 335,100,000 yr  (10 steps @ 100,000 yr)

  phase                          calls    d(count)       d(sum Hc) m       d(sum Hm) m
  ------------------------------------------------------------------------------------
  line_growth_shrink             ...
  line_end_stretch                ...
  line_end_arc_grow               ...
  line_end_retreat                ...
  line_end_accretion              ...
  line_interior_carve             ...
  line_regularization             ...
  convergent_deformation          ...
  ...

continental/oceanic node-type split (per phase, resolved against crust_type_code)
  phase                          cont d(sum Hc)   ocean d(sum Hc)    cont d(sum Hm)   ocean d(sum Hm)
  ---------------------------------------------------------------------------------------------------
  ...
```

**How to read it:**

- **`d(sum Hc)`/`d(sum Hm)`** are this interval's net change, summed over every node the phase
  touched (`sum_hc_after - sum_hc_before`), not a mean -- comparable directly across phases
  regardless of how many nodes each one happened to touch. A large `calls` count with a small
  `d(sum Hc)` means the phase ran often but nearly canceled out (e.g. `decompression_melting`
  resetting a column to almost the same reference thickness it already had); a small `calls`
  count with a large delta means a rare but individually large event (a merge, a relattice).
- **`d(count)`** is the touched slice's own node-count change -- nonzero only for phases that
  can add/remove nodes (`line_growth_shrink` and its five sub-phases below, `corner_notch_fill`,
  cleanup/merge/relattice/regularization); the pure Hc/Hm-mutation phases
  (`convergent_deformation`, `divergent_deformation`, `arc_magmatism`,
  `oceanic_cooling_relaxation`, `decompression_melting`) always report 0 here, so any thickness
  change they show is a real per-node change, not a population effect.
- **`line_growth_shrink`'s five sub-phases** -- `line_end_stretch` (endpoint stretch-thinning:
  `d(count)` always 0, a genuine per-node Hc/Hm change), `line_end_arc_grow` (brand-new
  arc-margin nodes: `d(count)` > 0, `d(sum)` is pure population growth at the seed column),
  `line_end_retreat` (a line end or interior run dropped: `d(count)` < 0, `d(sum)` the volume
  that left the line, whether or not it's later recovered), `line_end_accretion` (that volume,
  where flagged for suture accretion, reappearing on the surviving edge nodes: `d(count)` == 0,
  `d(sum)` the recovered share, capped -- so `line_end_retreat`'s loss plus `line_end_accretion`'s
  gain nets to the accretion cap's own overflow, not the whole retreated volume), and
  `line_interior_carve` (the oceanic-only mid-row subduction carve, never redistributed:
  `d(count)` < 0). `line_growth_shrink` itself is still recorded as the whole call's net total --
  a cross-check that these five sub-phases' own `d(count)`/`d(sum Hc)`/`d(sum Hm)` add up to its
  own (their raw before/after counts aren't directly comparable, since `line_growth_shrink`
  records the whole line on every call while each sub-phase records only the slice it touches).
- **the node-type split table** separates a genuine Hc/Hm change (both columns move) from pure
  reclassification (one column drops while the other rises by roughly the same amount, with
  little net change in the top table's combined total) -- `decompression_melting` is the
  phase most likely to show this, since it can flip a node's effective type.
- **`--json`**'s `scopes` dict on each phase row also has `continental_plate`/`oceanic_plate`
  (split by the *owning plate's* `crust_type` rather than each node's own effective type --
  the two disagree wherever a plate carries nodes of the other type, e.g. an accreted terrane)
  and `area_before_m2`/`area_after_m2` (`count * node_area_m2(node_density)`, the "nominal
  covered area" the issue's acceptance criteria ask for -- not runnable as an independent
  areal integral over the globe; see the issue's own caveat on why summed Hc isn't one either).
- This tool answers issue #216's first acceptance-criteria item (instrumentation); it does not
  by itself separate genuine geological sinks from numerical drift -- that's a reading exercise
  against the mechanism list in the issue itself, using this tool's own numbers as the evidence.

Test: [`unit_tests/test_phase_budget_diagnostics.py`](../backend/unit_tests/test_phase_budget_diagnostics.py).

---

## Debugging UI elements

### The `?deb` URL flag

Load the frontend with `?deb` anywhere in the URL query string (e.g.
`http://localhost:5173/?deb`) and a cluster of developer-only UI appears that is otherwise
**absent entirely, not just disabled** -- an ordinarily-generated world's sidebar and dialogs
stay uncluttered by default. It's a one-time flag read from `window.location.search` at module
load (`frontend/src/App.tsx`'s `DEBUG_UI` constant), not something that changes while the app
is open, and it works identically against the packaged prod server and the dev server -- it
only ever looks at the URL the browser loaded, never which backend is serving it.

`?deb` gates three things:

1. The Map View select's **"Debug >"** option group (below).
2. The **Event Console** and **corner-notch log** side panels (below).
3. The Generate World dialog's **"Debugging Worlds"** tab (below).

Everything documented under this heading requires `?deb` in the URL unless noted otherwise.

### Debug map views (Map View → "Debug >")

Like every other render view, most of these are just `GET /world/render?view=...` PNGs (see
[api-reference.md](api-reference.md)'s `/world/render`), billed the same as any render and
carrying no simulation side effects. Three of the ten entries in the group --
**Plate Inspector**, **Rivers**, and **Lake Inspector** -- are client-rendered JSON views
instead; they're primarily feature/inspection tools rather than diagnostic-metric views, but
they still live behind `?deb` because they were built alongside this tooling and haven't been
promoted to the always-visible "Maps" group.

#### Speckle (coastal-dither) overlay

**View:** `?view=speckle` · `render_image._render_speckle_view`

##### What it's for

The "speckled low-relief coastlines" problem (see
[GitHub issue #122](https://github.com/adubey/mantle-bloom/issues/122)): a marginally-submerged
flat shelf whose per-node elevation noise is larger than its own height above/below sea
level, so neighbouring nodes flip land↔ocean and the coast renders as a checkerboard instead
of a shoreline. On the Elevation / Biome / Combined views that just looks like a fuzzy,
slightly-noisy coast; there was no way to see *where* the coast is a genuine checkerboard vs.
a clean line, or to make a legible before/after for a coastal-feedback change without an
ad-hoc script. This view is that script, checked in.

##### The metric

For every elevation node in the raw plate node cloud (`plates.collect_all_points`), with
`sea_level = World.sea_level_m`:

- **class** = `elevation > sea_level` (land) vs. `<= sea_level` (ocean). Raw elevation only --
  no `hydrology` connectivity filter, deliberately, so this matches what the investigation
  scripts computed and so an enclosed sub-sea-level lake shore still shows up.
- **near** = `|elevation - sea_level| < SPECKLE_NEAR_BAND_M` (120 m). Only these nodes are
  drawn; everything else is just backdrop.
- **coastal-dither fraction** = of a near node's `SPECKLE_NEIGHBOR_K` (8) nearest neighbours,
  the share that are the *opposite* class from the node itself. `0.0` = the whole
  neighbourhood agrees (a coherent shoreline); higher = more disagreement.

`coastal_dither_fraction(points, elevation, sea_level_m) -> (fraction, near)` is a plain
module-level function -- call it directly from a probe script against a loaded `.mbworld`.

**Reading the numbers.** The metric's natural scale is set by taking *k* nearest neighbours
on an irregular 2D node cloud:

| situation | fraction |
|---|---|
| coherent shoreline (monotonic ramp across sea level) | `< ~0.35` |
| a perfect land/ocean checkerboard | `~0.5` (the 4 orthogonal neighbours flip, the 4 diagonals don't) |
| random per-node dither | `~0.5` |
| a genuinely isolated speck (one land node ringed entirely by ocean, or vice versa) | `→ 1.0` |

So `SPECKLE_FLAG_FRACTION` (0.75) flags **isolated specks** -- the single-pixel islands and
ponds -- not the mixed zone. The colour ramp puts the `~0.5` checkerboard band firmly in
"hot" (orange) territory below the flag threshold.

##### Reading the render

- **Backdrop:** muted olive land (`SPECKLE_LAND_BACKDROP_RGB`) / dark blue ocean
  (`SPECKLE_OCEAN_BACKDROP_RGB`), split at raw sea level.
- **Near-sea-level nodes:** a dot per node, coloured by fraction -- green (clean) → yellow →
  orange (`~checkerboard`) → red (approaching isolated). `speckle_colors()` /
  `_SPECKLE_STOP_F` / `_SPECKLE_STOP_RGB`.
- **Flagged nodes** (fraction ≥ 0.75): an oversized **magenta** square (`SPECKLE_FLAG_RGB`),
  so isolated specks stand out over the ramp even at a glance.

A clean coast reads as a thin green thread one node wide. A dithering drowned shelf reads as
a broad orange/red smear with magenta flecks. Inland lake shores also light up (the metric
doesn't know they aren't ocean) -- usually useful, occasionally noise.

##### Doing a before/after

```python
from pathlib import Path
from app import render_image
from app.persistence import load_world_bytes    # or world.generate_world + step_world

world = load_world_bytes(Path("~/Downloads/mantle-bloom-seed888151728-85000000y.mbworld").expanduser().read_bytes())
frac, near = render_image.coastal_dither_fraction(
    *render_image.plates.collect_all_points(world.plates)[:2], world.sea_level_m
)
print(f"near={near.sum()}  flagged={(frac >= 0.75).sum()}  mean_frac={frac[near].mean():.3f}")
open("/tmp/speckle.png", "wb").write(
    render_image.render_png(world, "behrmann", "speckle", 1400, 770)
)
```

Step the world N times each way (feedback change on vs. off) and compare `flagged`,
`frac[near].mean()`, and the two PNGs. A fix should drop `flagged` and `mean_frac` and visibly
thin the smear.

##### Constants (`render_image.py`)

| constant | default | meaning |
|---|---|---|
| `SPECKLE_NEAR_BAND_M` | 120.0 | half-width of the sea-level band the overlay draws |
| `SPECKLE_NEIGHBOR_K` | 8 | nearest neighbours averaged for the fraction |
| `SPECKLE_FLAG_FRACTION` | 0.75 | fraction at/above which a node gets the magenta marker |

#### Points (`platesDetail`)

`GET /world/render?view=platesDetail` draws each plate's raw `ElevationLine` nodes as dots
coloured by elevation (not the smoothed territory fill that the "Plates & Faults" view uses).
It is the fastest way to *see* lattice-level damage -- winding rows, stray one-node "teeth", a
staircase plate edge -- against the actual elevation field. See
[api-reference.md](api-reference.md) (`GET /world/render`) for the full view list.

#### Erosion & Deposition (`geomorph`)

`GET /world/render?view=geomorph` colours every node by its net elevation change over the
last step -- `erosion.ErosionResult.net_elevation_change_m` (post-erosion elevation minus
pre-erosion, so erosion minus every deposition pathway plus the small flatten/lake-siltation
terms; *not* tectonic deform, isostasy, or volcanism), retained on `World.erosion_cache`
purely for this view. A diverging scale: warm brown/orange where the step net-lowered a node,
cool blue where it net-raised one, a flat neutral grey in the +-few-metre band so only the
lumps stand out, clamped past +-60 m/step. The coastline is overlaid for orientation.

What it's for: the per-step deposition in the near-sea-level band is wildly lumpy -- a
+200 m spike on one node, ~0 on its neighbour -- which is the mechanism behind the coastal
checkerboard (see [GitHub issue #122](https://github.com/adubey/mantle-bloom/issues/122)),
but is invisible in every other view. Step the world once with climate & biomes on, then
switch to this view and look along a drowned shelf: a clean coastal plain deposits smoothly
(uniform pale colour), a dithering one shows a salt-and-pepper mix of saturated warm and cool
cells. Use it as a before/after for any coastal-feedback change instead of an ad-hoc script.

`erosion_cache` is `None` until the first climate/erosion step (and on a freshly loaded save
-- it isn't persisted), where the view falls back to a flat neutral field plus the coastline
rather than erroring.

#### `elevReason` render view (last elevation change)

`GET /world/render?view=elevReason` colours every node by
`ElevationLine.elev_change_reason` -- one categorical
`elevation_lines.ELEV_CHANGE_*` code per node recording *which process last moved that node's
elevation* by more than `ELEV_CHANGE_MIN_DELTA_M` (2 m in a step). Warm hues = crust being
built (collision / subduction-arc / transform / rift / new crust / volcano), cool blues =
crust being planed down or buried (erosion / deposition / coastal planation / submarine),
pale = glacial flattening, plain grey = **untouched since generation**. Coastline overlaid
for orientation. `render_image._render_elev_reason_view` / `_ELEV_REASON_RGB`; labels in
`elevation_lines.ELEV_CHANGE_LABELS`, frontend legend in `legendData.ts`'s
`ELEV_REASON_ENTRIES` (hand-synced, same precedent as the biome palette).

What it's for: "there should be more terrain features -- why is so much of this world flat?"
The `geomorph` view above shows *this step's* rate; this shows the *standing* provenance
accumulated over the run. A large grey (NONE) expanse on land means that terrain was never
tectonically built -- its only relief is the generation-time noise texture, slowly being
worn/buried away. Large erosion/deposition/coastal-leveling expanses mean it *is* being
actively flattened now. Collision / subduction-arc / rift belts are where relief is still
being made -- if those are thin or absent while the continents are large and quiescent, the
tectonic engine has stalled and nothing is replacing the relief erosion removes.

Provenance is written by `lithosphere_plate.deform` (tectonic codes, re-stamped
every step a belt stays active), `volcanism` (eruptions), and `erosion` (geomorphic codes).
A structural code is **sticky**: erosion only overrides it when this step's net geomorphic
change is itself large (`ELEV_CHANGE_STRUCTURAL_OVERRIDE_M_PER_MYR`, ~100 m/Myr), so ordinary
background wash on an actively-rising range doesn't relabel it "erosion".

Unlike `geomorph`, the field is *persistent* (rides the plate through rotation/split/merge,
survives save/load). But a save written before this field existed -- or a world never stepped
since -- reads all-NONE and fills in over the next few steps. `ELEV_CHANGE_MIN_DELTA_M` /
`ELEV_CHANGE_STRUCTURAL_OVERRIDE_M_PER_MYR` (elevation_lines.py) tune the two thresholds.

#### `overlapAge` render view (plate overlap onset)

`GET /world/render?view=overlapAge` draws a muted land/ocean backdrop (same grid as
Elevation) overlaid with one dot per node that is **currently** sitting on top of another
plate's territory, coloured by how long it has been --
`world.elapsed_years - ElevationLine.overlap_onset_years`
(`render_image._render_overlap_age_view` / `overlap_age_colors`). Pale yellow = a fresh
overlap (transient envelope slop, self-correcting); deepening through orange to
magenta-purple = stuck for tens of Myr (a real stalled collision the merge path never
resolves). Clamped at 60 Myr.

`overlap_onset_years` is a per-node `ElevationLine` field stamped every step by
`merge_split.update_overlap_tracking`, which goes through the same
`plates.compute_node_overlap` the Plate Inspector's own overlap fractions use -- so the map,
the `since_years` number in `GET /world/plates` / `App.tsx` (`#21 (15%, since 178 My)`), and
`python -m app.plate_diagnostics`'s `21 -> 0  15.2%  since 178.0 My` line are all the same
underlying node set. A save written before the field existed loads with every onset at 0
(`ElevationLine.__getattr__` backfill) and shows an all-backdrop map / `since` omitted
until it is stepped. An all-backdrop map on a stepped world is the healthy case.

The overlap itself now also *heals*: `torque.classify_boundary_nodes` tests deep-interior
nodes (inside a neighbour's bounding sphere, not just within `reach_rad` of one) for polygon
containment, so a plate that has slid well over another has those deep nodes classified
`contested` -- continental crust thickens there (`rheology.apply_convergent_deformation` ->
mountain uplift), oceanic crust subducts. Before this, only the boundary-local band was
classified and a deep overlap just sat.

#### Added/removed points (`nodeAge`)

`GET /world/render?view=nodeAge` draws the same muted land/ocean backdrop as `overlapAge`,
overlaid with two independent dot layers: warm dots (pale yellow -> amber -> burnt orange ->
dark rust) for still-live nodes created recently, and cool dots (pale blue -> sky blue -> deep
blue -> navy) for nodes no longer part of any plate that were removed recently
(`render_image._render_node_age_view` / `node_added_colors` / `node_removed_colors`). Both
ramps are clamped at 20 Myr -- older activity fades to plain backdrop rather than pinning at a
saturated top-of-ramp colour that would misleadingly read as "still happening now." Answers
"where has the lattice actually been changing, and in which direction" -- exactly the question
a persistent gap raises: is either side even *trying* to grow into it, or has activity nearby
gone quiet?

**Added** comes from `ElevationLine.node_created_years`, a permanent, write-once per-node
timestamp stamped once at
`lithosphere_plate.LithospherePlate._seed_and_erupt_new_nodes` -- the one choke point every
node-creation call site (`_claim_adjacent_territory`, `_fill_corner_notch`, the arc-margin
append branch in `_grow_or_shrink_line_for_deform`) funnels through, so this view is
comprehensive across every place brand-new crust actually originates, not just one code path.
Unlike `overlap_onset_years`, it never reverts: a node's birth date is permanent, carried
through `regularize_line`'s resample by nearest-neighbour (never interpolated -- averaging two
birth years would invent a meaningless date). The sentinel is `-1.0`, not `0.0` -- a save's
initial genesis nodes (from generation, not from a tracked creation event) read as "predates
tracking," not "created at year 0," and a save written before this field existed backfills the
same way. `plates.collect_all_node_created_years` is the accessor.

**Removed** has no live node to carry a field on, so it comes from a short-lived side buffer
instead: `World.removed_points_log`, a `(world_xyz, removed_years, plate_id)` list appended to
by `World.record_removed_points` at every node-removal site --
`_grow_or_shrink_line_for_deform`'s end-retreat and interior-subduction carve, `merge_split.
merge_plates` (the absorbed plate's own points, recorded before the fused resample), `Plate.
defragment` (stranded fragments shed below `min_fragment_nodes`), and `merge_split.
remove_defunct_plates` (a whole plate's points, when it's dropped). Capped by count
(`MAX_REMOVED_POINTS_LOG`, 20,000), not age -- subduction removes nodes essentially every step
on a full-size save, so an age-based window would grow unboundedly at high `node_density`.
Backfills to an empty list on a save written before this field existed.

All-backdrop (no dots of either colour) means nothing has been created or removed recently --
the common case between bursts of boundary activity, *not* necessarily unhealthy the way an
empty `overlapAge` is (that view's empty case specifically means "no stuck overlaps," a
positive signal; this one's empty case just means "quiet right now").

##### Click-to-inspect -- `GET /world/node_at`

Clicking a point on the `nodeAge` view calls `GET /world/node_at?lat_deg&lon_deg`, which --
unlike `/world/sample_at` (climate-grid lookup) -- reports the actual nearest live
`ElevationLine` node's own plate-local `phi`/`theta` (recovered by inverse-transforming its
world position through its owning plate's frame, which round-trips exactly to the node's own
creation-time value), elevation, and `node_created_years`, plus -- independently, since a
removed point has no live node to report instead -- the nearest `removed_points_log` entry if
one sits within `_REMOVED_POINT_MATCH_RAD` (~127 km) of the click. `plates.nearest_node_index`
is the shared hit-test helper (mirrors `nearest_plate_id`, but returns the node's own index
into the `collect_all_points` order rather than just its owning plate id, so a caller wanting
both plate ownership and per-node fields needn't build the k-d tree twice).

#### Plate Inspector, Rivers, Lake Inspector, Plates & Faults

The remaining four entries in the "Debug >" group are documented as their own feature views,
not here:

- **Plate Inspector** (`GET /world/plates`, `frontend/src/PlateInspector.tsx`) plots every
  plate's nodes and outline client-side and reports the same motion / shape / overlap numbers
  the CLI plate dump does -- see [Plate Inspector diagnostic fields](#plate-inspector-fields)
  below for the field walkthrough, and
  [simulation-model.md#plate-inspector](simulation-model.md#plate-inspector) for the panel
  layout and bounding-ellipse fit.
- **Rivers** and **Lake Inspector** (`GET /world/rivers` / `GET /world/lakes`,
  `RiverInspector.tsx` / `LakeInspector.tsx`) render flow networks and lake basins from raw
  JSON -- see
  [simulation-model.md#river-inspector](simulation-model.md#river-inspector) and
  [simulation-model.md#lake-inspector](simulation-model.md#lake-inspector).
- **Plates & Faults** (`GET /world/faults`, `GET /world/earthquakes`,
  `frontend/src/PlatesAndFaults.tsx`) merges plate outlines with the fault-trace and
  earthquake overlays -- see
  [simulation-model.md#fault-inspector](simulation-model.md#fault-inspector).

All four are primarily inspection/feature tools rather than pathology-hunting diagnostics, but
all four currently require `?deb` to reach, same as the views above.

<a id="plate-inspector-fields"></a>
##### Plate Inspector diagnostic fields

Per plate, alongside the geometry: `speed_cm_per_yr` + `at_max_rate` (railed at
`mantle.MAX_PLATE_RATE`, shown red), `euler_pole` (lat/lon), `age_steps`,
`median_elevation_m` + `submerged_fraction` (red when a continental plate is >50% under
water), `overlaps` (which other plates this one's territory sits on top of, by what fraction
of its own nodes, and `since_years` -- the earliest `elapsed_years` any still-overlapping node
first went over another plate; `main._plate_overlaps` / `ElevationLine.overlap_onset_years`,
see the `overlapAge` view above), and `collisions` (`world.collision_progress` timers
involving the plate). Full field reference: [api-reference.md](api-reference.md)
(`GET /world/plates`).

Use the Inspector for the *visual* read -- concentric rings (pole winding), a plate's dots
sitting inside a neighbour's (overlap), a long straight sawtooth chord across open ocean (an
over-extended lattice) -- and `python -m app.plate_diagnostics` for the numbers behind it.

### Event Console & corner-notch log panel

Two side panels, both rendered only when `?deb` is set (`frontend/src/App.tsx`'s
`EventConsole` / `CornerNotchLogPanel`). Without `?deb` the underlying data is still reachable
through the API -- `GET /world/summary`'s `events` field, and `GET /world/corner_notch_log`
directly -- only the in-app viewers are hidden.

<a id="event-log"></a>
#### Event log

The `events` list on `GET /world/summary`, shown in the Event Console, logs lake
formation/splits and other discrete events.

**Lake-churn aggregation -- `lakes.summarize_lake_events`.** On a long run over a dithering
low-relief coast the lake solver produces hundreds of near-sea-level transient merge/split
transitions per My -- one pair per puddle per step (see
[GitHub issue #122](https://github.com/adubey/mantle-bloom/issues/122)). Left raw, these flood
the console and bury real basin/tectonic events.

`lakes.step_lakes` returns structured `lakes.LakeEvent`s (`kind` / `node_count` /
`elevation_m` / `basin_count`, plus a `.message` property with the same wording as before)
instead of pre-formatted strings. [`erosion.py`](../backend/app/erosion.py) runs a step's
events through `lakes.summarize_lake_events(events, world.sea_level_m)` before logging:

- A transition whose water surface is more than `lakes.NEAR_SEA_LEVEL_EVENT_BAND_M`
  (**15 m**) from the current sea level is a genuine basin event -- logged individually,
  unchanged.
- Transitions **within** that band are the coastal-pond churn. A lone one still logs
  verbatim; **two or more in one step collapse to a single line** --
  `"38 transient coastal ponds churned near sea level this step (22 merged, 16 split)."`

So a persistent deep endorheic basin (e.g. a real ~435-node lake oscillating near -1770 m)
stays visible in the console while the checkerboard shelf contributes at most one aggregate
line per step. The band is measured against `world.sea_level_m`, so it tracks a sea-level
control change. This aggregation always runs (it isn't gated by `?deb`) -- only the console
that displays the result is. Tests:
[`unit_tests/test_lakes.py`](../backend/unit_tests/test_lakes.py)
(`test_summarize_lake_events_*`).

#### Corner-notch decision log (`_fill_corner_notch_frontier`)

`GET /world/corner_notch_log` exposes a verbose, structured, per-call record of what
`lithosphere_plate.LithospherePlate._fill_corner_notch_frontier` -- the triple-junction/
diagonal-residual gap-filling fallback `_stretch_end` and `_claim_adjacent_territory`
structurally can't reach -- actually decided each time it ran, and why. Deliberately **not**
part of the always-on Event Console: this can fire once per plate per step, far higher volume
than that log is meant to carry (the same reasoning behind `lakes.summarize_lake_events`'s own
aggregation, just solved here by giving the verbose detail its own separate, off-by-default
channel instead of collapsing it).

Recording is gated by `World.debug_diagnostics` (`bool`, off by default), toggled via
`POST /world/controls`'s `debug_diagnostics` field or the **Controls window → Tectonics tab →
"Log corner-notch decisions"** checkbox -- this checkbox is visible regardless of `?deb`, but
the panel that displays what it logs (below the Event Console) only renders with `?deb`, so
flipping it on without `?deb` records data with nowhere in the UI to see it (still readable via
`GET /world/corner_notch_log`). While off, `World.log_corner_notch` is a no-op and nothing is
recorded, so there's no cost on an ordinary play session. Entries land in
`World.corner_notch_log`, capped by count (`MAX_CORNER_NOTCH_LOG_LENGTH`, 2000) like
`World.events`.

##### Reading an entry

Every entry carries `plate_id`, `outcome`, `nodes_added`, `elapsed_years`, and `algorithm:
"frontier"`; most also carry enough of the call's own geometry to place it:

| `outcome` | Meaning |
|---|---|
| `no_neighbours` | Early return -- no real neighbour plate nearby at all, so "uncovered space next to my own edge" would just be the rest of the sphere (the guard that stops a lone plate from growing its entire perimeter every step). |
| `no_own_lines` | Early return -- this plate has no nodes/lines to notch-fill from. |
| `no_candidate_rows` | The scanned window (`window_rad`, `phi_lo`/`phi_hi`) found no lattice point that's both near a neighbour and not already covered by anyone -- the healthy, common case for an already-well-tiled boundary. |
| `claimed` | Ended with `nodes_added > 0` new nodes appended -- `gap_fill_frontier.fill_gap_by_growing_plates` grew this plate's own lines (or opened new ones) to cover them. |
| `no_claim` | Candidate rows existed but `fill_gap_by_growing_plates` claimed nothing (e.g. the connect-radius walk stalled immediately). |

##### Using it on a real save

Load a save, enable diagnostics, then step it -- the panel fills in per-plate-per-step, so
watching a specific known-bad junction's plate ids across several steps shows directly whether
`_fill_corner_notch_frontier` is even attempting that boundary (`no_neighbours`/`no_own_lines`
would mean it never gets that far), finding nothing to claim
(`no_candidate_rows`/`no_claim`), or claiming a window that turns out too small (`claimed` with
a low `nodes_added` relative to the gap's real size) -- each a different next step for a fix,
instead of guessing blind from the rendered map alone.

### "Debugging Worlds" tab -- scripted plate scenarios

The Generate World dialog's fourth tab (`POST /world/generate_debug`,
`backend/app/debug_worlds.py`) builds tiny, low-resolution (`node_density=0.5`, the coarsest
real choice), hand-placed plate configurations for fast iteration on the gap-filling problem,
independent of any real save's history. `GET /world/debug_scenarios` lists the current
options; picking one and pressing Generate replaces the current world exactly like the
Random/Human-made/Premade tabs, with `debug_diagnostics` already on (see the corner-notch log
above).

Each plate's motion is **pinned**, not torque-driven: `World.pinned_omegas` (`plate_id ->` a
fixed angular velocity) makes `LithospherePlate.shift` use that value verbatim every step,
bypassing `torque.shift_plate`'s own ridge-push/slab-pull/basal-drag recompute entirely for
that plate (`torque.apply_omega_and_rotate` is the shared rotate-and-report tail both the real
and pinned paths use). This is what makes a scenario reliable to iterate on: "these two plates
diverge" stays true step after step, rather than however the force balance against an
otherwise-irrelevant `mantle_centers` field happens to settle it. `pinned_omegas` is empty for
every ordinarily-generated world -- this is purely a debug-world mechanism.

A scenario is a set of seed lat/lons, each plate's crust type, and a list of pairwise
`(plate_i, plate_j, "divergent" | "convergent")` boundary relationships;
`debug_worlds._omegas_from_relationships` converts that into each plate's own single rotation
(the sum of its own boundaries' tangential "push" contributions, converted to a rotation axis
via `omega = cross(seed_xyz, tangent_direction)`) so a plate sitting at a junction of several
boundaries gets one rotation reflecting all of them at once, the same way a real plate's
motion is a net result of every boundary force acting on it -- just picked directly here
instead of emerging from `torque.py`'s force balance.

Current scenarios (`debug_worlds.DEBUG_SCENARIOS`):

| Scenario | Layout | Purpose |
|---|---|---|
| `two_plate_divergent` | 2 plates, straight rift | Baseline: ordinary end-growth (`_stretch_end`) should close this alone -- `_fill_corner_notch_frontier` should log mostly `no_claim`/`no_candidate_rows`, never a real gap. |
| `two_plate_convergent` | 2 plates, closing boundary | Subduction/collision baseline. |
| `triple_junction_mixed` | 3 plates, one junction, 2 divergent legs + 1 convergent | **The exact case `_fill_corner_notch_frontier`'s own docstring calls out** (confirmed on a real save, seed 430031492) -- confirmed (see `test_triple_junction_mixed_scenario_exercises_fill_corner_notch_frontier`) to drive real `claimed` activity within a handful of steps. Start here. |
| `four_plate_grid` | 4 plates, 2x2, all edges divergent | Four simultaneous triple-junction-like corners at once, around one shared center point. |
| `five_plate_irregular` | 5 plates, irregular ring, mixed relationships | Closest single scenario to a real save's messiness while staying small enough to iterate on quickly. |

Workflow for chasing a specific corner-notch failure mode: enable `?deb`, generate
`triple_junction_mixed`, step it a handful of times (Play/Stop, small "Years per step" for
fine-grained observation), and watch the "Added/removed points" and "Plate overlap age" views
alongside the corner-notch log panel together -- the combination this whole diagnostic suite
was built to let you read at once, rather than switching between four separate tools with no
shared time axis.

These tiny scripted scenarios don't naturally exercise
`gaps.fill_gaps_by_growing_neighbours` (every plate starts adjacent to every gap that can open
here) -- checking the whole-sphere sweep needs a real or synthetically-vacated-plate world
instead, see `unit_tests/test_gaps.py`'s `fill_gaps_by_growing_neighbours` tests for the shape
of that setup.

---

## Still worth building

Two tools were scoped out but never started:

1. **A map-view render of the stranded-basin report.** `GET /world/stranded_basins` already
   returns `centroid_xyz`/`floor_xyz` ready for one -- currently the report is CLI/JSON-only
   (see [`python -m app.stranded_basins`](#stranded-basins-cli) above), with no map view to
   see the basins in place.
2. **A `/world/sample_at` field for the `geomorph` / `elevReason` click-popup.** The popup is
   only wired for elevation/biome/combined today (see `GET /world/sample_at` in
   [api-reference.md](api-reference.md)).

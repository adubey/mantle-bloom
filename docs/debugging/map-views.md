# Debug map views

These views appear under Map View → "Debug >" when the page is loaded with
[`?deb`](debug-ui.md#the-deb-url-flag).

Like every other render view, most of these are just `GET /worlds/{world_id}/render?view=...` PNGs (see
[api-reference.md](../api-reference.md)'s `/worlds/{world_id}/render`), billed the same as any render and
carrying no simulation side effects. Three of the ten entries in the group --
**Plate Inspector**, **Rivers**, and **Lake Inspector** -- are client-rendered JSON views
instead; they're primarily feature/inspection tools rather than diagnostic-metric views, but
they still live behind `?deb` because they were built alongside this tooling and haven't been
promoted to the always-visible "Maps" group.

## Speckle (coastal-dither) overlay

**View:** `?view=speckle` · `render_image._render_speckle_view`

### What it's for

The "speckled low-relief coastlines" problem (see
[GitHub issue #122](https://github.com/adubey/mantle-bloom/issues/122)): a marginally-submerged
flat shelf whose per-node elevation noise is larger than its own height above/below sea
level, so neighbouring nodes flip land↔ocean and the coast renders as a checkerboard instead
of a shoreline. On the Elevation / Biome / Combined views that just looks like a fuzzy,
slightly-noisy coast; there was no way to see *where* the coast is a genuine checkerboard vs.
a clean line, or to make a legible before/after for a coastal-feedback change without an
ad-hoc script. This view is that script, checked in.

### The metric

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

### Reading the render

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

### Doing a before/after

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

### Constants (`render_image.py`)

| constant | default | meaning |
|---|---|---|
| `SPECKLE_NEAR_BAND_M` | 120.0 | half-width of the sea-level band the overlay draws |
| `SPECKLE_NEIGHBOR_K` | 8 | nearest neighbours averaged for the fraction |
| `SPECKLE_FLAG_FRACTION` | 0.75 | fraction at/above which a node gets the magenta marker |

## Points (`platesDetail`)

`GET /worlds/{world_id}/render?view=platesDetail` draws each plate's raw terrain nodes (quad cell centres)
as dots coloured by elevation (not the smoothed territory fill that the "Plates & Faults" view
uses). It is the fastest way to *see* lattice-level damage -- stray one-cell "teeth", ragged or
overlapping plate edges -- against the actual elevation field. Clicking a dot
(`GET /worlds/{world_id}/elevation_point_at`) reports the nearest node's plate, plate-local `phi`/`theta`
and elevation, and marks it on the map. See
[api-reference.md](../api-reference.md) (`GET /worlds/{world_id}/render`) for the full view list.

## Erosion & Deposition (`geomorph`)

`GET /worlds/{world_id}/render?view=geomorph` colours every node by its net elevation change over the
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

## `elevReason` render view (last elevation change)

`GET /worlds/{world_id}/render?view=elevReason` colours every node by its
`elev_change_reason` field -- one categorical
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

Provenance is written by the tectonic column pass (`lithosphere_plate.deform_columns`;
tectonic codes, re-stamped every step a belt stays active), `volcanism` (eruptions), and `erosion` (geomorphic codes).
A structural code is **sticky**: erosion only overrides it when this step's net geomorphic
change is itself large (`ELEV_CHANGE_STRUCTURAL_OVERRIDE_M_PER_MYR`, ~100 m/Myr), so ordinary
background wash on an actively-rising range doesn't relabel it "erosion".

Unlike `geomorph`, the field is *persistent* (rides the plate through rotation/split/merge,
survives save/load). But a save written before this field existed -- or a world never stepped
since -- reads all-NONE and fills in over the next few steps. `ELEV_CHANGE_MIN_DELTA_M` /
`ELEV_CHANGE_STRUCTURAL_OVERRIDE_M_PER_MYR` (elevation_lines.py) tune the two thresholds.

## `overlapAge` render view (plate overlap onset)

`GET /worlds/{world_id}/render?view=overlapAge` draws a muted land/ocean backdrop (same grid as
Elevation) overlaid with one dot per node that is **currently** sitting on top of another
plate's territory, coloured by how long it has been --
`world.elapsed_years - overlap_onset_years`
(`render_image._render_overlap_age_view` / `overlap_age_colors`). Pale yellow = a fresh
overlap (transient envelope slop, self-correcting); deepening through orange to
magenta-purple = stuck for tens of Myr (a real stalled collision the merge path never
resolves). Clamped at 60 Myr.

`overlap_onset_years` is a per-node surface field stamped every step by
`merge_split.update_overlap_tracking`, which goes through the same
`plates.compute_node_overlap` the Plate Inspector's own overlap fractions use -- so the map,
the `since_years` number in `GET /worlds/{world_id}/plates` / `App.tsx` (`#21 (15%, since 178 My)`), and
`python -m app.plate_diagnostics`'s `21 -> 0  15.2%  since 178.0 My` line are all the same
underlying node set. A save written before the field existed loads with every onset at its
default and shows an all-backdrop map / `since` omitted until it is stepped. An all-backdrop map on a stepped world is the healthy case.

The overlap itself now also *heals*: `torque.classify_boundary_nodes` tests deep-interior
nodes (inside a neighbour's bounding sphere, not just within `reach_rad` of one) for polygon
containment, so a plate that has slid well over another has those deep nodes classified
`contested` -- continental crust thickens there (`rheology.apply_convergent_deformation` ->
mountain uplift), oceanic crust subducts. Before this, only the boundary-local band was
classified and a deep overlap just sat.

## Added/removed points (`nodeAge`)

`GET /worlds/{world_id}/render?view=nodeAge` draws the same muted land/ocean backdrop as `overlapAge`,
overlaid with two independent dot layers: warm dots (pale yellow -> amber -> burnt orange ->
dark rust) for still-live nodes created recently, and cool dots (pale blue -> sky blue -> deep
blue -> navy) for nodes no longer part of any plate that were removed recently
(`render_image._render_node_age_view` / `node_added_colors` / `node_removed_colors`). Both
ramps are clamped at 20 Myr -- older activity fades to plain backdrop rather than pinning at a
saturated top-of-ramp colour that would misleadingly read as "still happening now." Answers
"where has the lattice actually been changing, and in which direction" -- exactly the question
a persistent gap raises: is either side even *trying* to grow into it, or has activity nearby
gone quiet?

**Added** comes from the `node_created_years` surface field, a permanent, write-once per-node
timestamp stamped when a cell is inserted (`quad_tectonics._new_cell_fields`, which boundary
advance, rift opening and gap filling all go through), so this view covers every place
brand-new crust originates. Unlike `overlap_onset_years`, it never reverts: a node's birth date
is permanent, carried through merges and remaps by nearest neighbour (never interpolated --
averaging two birth years would invent a meaningless date). The sentinel is `-1.0`, not `0.0`
-- a world's initial genesis nodes (from generation, not from a tracked creation event) read as
"predates tracking," not "created at year 0," and a save written before this field existed
backfills the same way. `plates.collect_all_node_created_years` is the accessor.

**Removed** has no live node to carry a field on, so it comes from a short-lived side buffer
instead: `World.removed_points_log`, a `(world_xyz, removed_years, plate_id)` list appended to
by `World.record_removed_points` at every node-removal site -- boundary retreat and
interior subduction (`quad_tectonics._retreat`), `merge_split.merge_plates` (the absorbed
plate's own points), `Plate.defragment` (stranded fragments cut off below
`min_fragment_nodes`, whether later accreted or dropped), and
`merge_split.remove_defunct_plates` (a whole plate's points, when it's dropped). Capped by count
(`MAX_REMOVED_POINTS_LOG`, 20,000), not age -- subduction removes nodes essentially every step
on a full-size save, so an age-based window would grow unboundedly at high `node_density`.
Backfills to an empty list on a save written before this field existed.

All-backdrop (no dots of either colour) means nothing has been created or removed recently --
the common case between bursts of boundary activity, *not* necessarily unhealthy the way an
empty `overlapAge` is (that view's empty case specifically means "no stuck overlaps," a
positive signal; this one's empty case just means "quiet right now").

### Click-to-inspect -- `GET /worlds/{world_id}/node_at`

Clicking a point on the `nodeAge` view calls `GET /worlds/{world_id}/node_at?lat_deg&lon_deg`, which --
unlike `/worlds/{world_id}/sample_at` (climate-grid lookup) -- reports the actual nearest live
node's own plate-local `phi`/`theta` (recovered by inverse-transforming its
world position through its owning plate's frame, which round-trips exactly to the node's own
creation-time value), elevation, and `node_created_years`, plus -- independently, since a
removed point has no live node to report instead -- the nearest `removed_points_log` entry if
one sits within `_REMOVED_POINT_MATCH_RAD` (~127 km) of the click. `plates.nearest_node_index`
is the shared hit-test helper (mirrors `nearest_plate_id`, but returns the node's own index
into the `collect_all_points` order rather than just its owning plate id, so a caller wanting
both plate ownership and per-node fields needn't build the k-d tree twice).

## Cratons (`craton`)

`GET /worlds/{world_id}/render?view=craton` (Debug > Cratons) colours every node by its craton state
(`render_image._render_craton_view` / `craton_colors`, see `cratons.py`). Deep red marks
cratons that predate the run, seeded at generation. Orange marks cratons that formed during
it. Both fade toward grey as a column's cratonic share thins. Ordinary continental crust is
grey, warming to amber as its formation clock (`stable_continental_myr`) approaches
`CRATON_FORMATION_MYR`. Teal marks belts, welded or rifted crust that never forms craton
(`cratons.belt_mask`, RFC 0009). Oceanic crust is dark blue.

Destruction is recorded in `World.craton_ledger`, one account per mechanism (rifting,
subduction, delamination, collisional reworking, erosion, topology removal), and
`cratons.balance_error_m3` checks it against the live volume. `audit_clipped_m3` is craton a
catch-all phase audit had to clip to its column, kept apart from the real mechanisms (#340).
`unattributed_m3` should stay near zero; anything there escaped every instrumented site.
Every stats snapshot carries the ledger (`craton_ledger_km3`), so the Stats window's *Craton
destroyed (cumulative)* series shows when cratons were lost. A step whose losses to one
mechanism exceed `CRATON_EVENT_FRACTION` of the live cratons is also written to the Event
Console.

To attribute a loss this view shows, see [Craton scripts](scripts-cratons.md): a save's
ledger history (`ledger_history.py`), the site behind each booking (`split_craton_rework.py`),
the writer behind an audit clip (`probe_craton_clips.py`) and the cause of a merge's loss
(`measure_merge_placement.py`). Losses at plate contacts are measured by the scripts in
[Collision and plate-motion scripts](scripts-collisions.md).

## Plate Inspector, Rivers, Lake Inspector, Plates & Faults

The remaining four entries in the "Debug >" group are documented as their own feature views,
not here:

- **Plate Inspector** (`GET /worlds/{world_id}/plates`, `frontend/src/PlateInspector.tsx`) plots every
  plate's nodes and outline client-side and reports the same motion / shape / overlap numbers
  the CLI plate dump does -- see [Plate Inspector diagnostic fields](#plate-inspector-fields)
  below for the field walkthrough, and
  [simulation-model/rendering-and-inspectors.md#plate-inspector](../simulation-model/rendering-and-inspectors.md#plate-inspector) for the panel
  layout and bounding-ellipse fit.
- **Rivers** and **Lake Inspector** (`GET /worlds/{world_id}/rivers` / `GET /worlds/{world_id}/lakes`,
  `RiverInspector.tsx` / `LakeInspector.tsx`) render flow networks and lake basins from raw
  JSON -- see
  [simulation-model/hydrology.md#river-inspector](../simulation-model/hydrology.md#river-inspector) and
  [simulation-model/hydrology.md#lake-inspector](../simulation-model/hydrology.md#lake-inspector).
- **Plates & Faults** (`GET /worlds/{world_id}/faults`, `GET /worlds/{world_id}/earthquakes`,
  `frontend/src/PlatesAndFaults.tsx`) merges plate outlines with the fault-trace and
  earthquake overlays. Its Overlays panel also toggles faults, cratons and a translucent
  ocean/lake layer (cratons show through it), and clicking the legend's Cratons row isolates
  the cratons the way a fault-type row isolates that regime -- see
  [simulation-model/rendering-and-inspectors.md#fault-inspector](../simulation-model/rendering-and-inspectors.md#fault-inspector).

All four are primarily inspection/feature tools rather than pathology-hunting diagnostics, but
all four currently require `?deb` to reach, same as the views above.

<a id="plate-inspector-fields"></a>
### Plate Inspector diagnostic fields

Per plate, alongside the geometry: `speed_cm_per_yr` + `at_max_rate` (railed at
`mantle.MAX_PLATE_RATE`, shown red), `euler_pole` (lat/lon), `age_years`,
`median_elevation_m` + `submerged_fraction` (red when a continental plate is >50% under
water), `overlaps` (which other plates this one's territory sits on top of, by what fraction
of its own nodes, and `since_years` -- the earliest `elapsed_years` any still-overlapping node
first went over another plate; `main._plate_overlaps` / the `overlap_onset_years` field,
see the `overlapAge` view above), and `collisions` (`world.collision_progress` timers
involving the plate). Full field reference: [api-reference.md](../api-reference.md)
(`GET /worlds/{world_id}/plates`).

Use the Inspector for the *visual* read -- concentric rings (pole winding), a plate's dots
sitting inside a neighbour's (overlap), a long straight sawtooth chord across open ocean (an
over-extended lattice) -- and `python -m app.plate_diagnostics` for the numbers behind it.

## Still worth building

Two tools were scoped out but never started:

1. **A map-view render of the stranded-basin report.** `GET /worlds/{world_id}/stranded_basins` already
   returns `centroid_xyz`/`floor_xyz` ready for one -- currently the report is CLI/JSON-only
   (see [`python -m app.hydroclimate.stranded_basins`](offline-cli.md#stranded-basins-cli)), with no map view to
   see the basins in place.
2. **A `/worlds/{world_id}/sample_at` field for the `geomorph` / `elevReason` click-popup.** The popup is
   only wired for elevation/biome/combined today (see `GET /worlds/{world_id}/sample_at` in
   [api-reference.md](../api-reference.md)).

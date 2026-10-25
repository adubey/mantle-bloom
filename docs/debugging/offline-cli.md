# Offline diagnostic CLIs

Four offline, read-only dumps of a saved world. Each **never starts the server or binds a
port** -- three load a `.mbworld` file directly (same pickle format as *File > Load*, see
[`persistence.py`](../../backend/app/persistence.py)) and print to stdout; `phase_budget_diagnostics`
can do that too, but (unlike the other three, which only ever read a save) it also *steps* the
world forward by a given number of years before reporting -- see its own section below. Each of
the first three reuses the exact code behind the equivalent API endpoint / UI panel, so the CLI
and the UI can never disagree.

```bash
cd backend
source .venv/bin/activate
python -m app.plate_diagnostics <save.mbworld>            # plate geometry
python -m app.hydroclimate.stranded_basins <save.mbworld>            # land-locked sub-sea-level pits
python -m app.hydroclimate.lake_hierarchy_diagnostics <save.mbworld> # lake-merge-forest depth/size
python -m app.phase_budget_diagnostics <save.mbworld> --years 1000000  # per-phase Hc/Hm budget
```

Every tool also takes `--json` for structured, scriptable output.

## `python -m app.plate_diagnostics` -- plate geometry

Reuses `main._plate_summary` / `main._plate_overlaps` -- the exact code path behind
`GET /worlds/{world_id}/plates` and the in-app Plate Inspector.

```
mantle-bloom plate diagnostics
  seed:          505070493
  elapsed:       25,500,000 yr  (~255 steps @ 100 ky)
  node_density:  4.0
  sea level:     0.0 m
  plates:        18

per-plate  ( * = railed at MAX_PLATE_RATE 15.0 cm/yr;  ! = continental & >50% submerged )
   id  crust         nodes   age     speed     pole lat,lon   med.elev  submrg
  ----------------------------------------------------------------------------
    0  oceanic        3364   255     4.59     +21.7,   -6.9      -3997   0.99
    1  continental    6843   255     2.43     -24.1, +169.3       1391   0.08
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
- **territory overlaps.** `A -> B  X%` = X% of plate A's own nodes have their centres inside
  B's territory (plates that merely share a boundary don't count). The text dump hides entries below 0.5%; `--json`
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

Tests: [`unit_tests/test_plate_diagnostics.py`](../../backend/unit_tests/test_plate_diagnostics.py).

<a id="stranded-basins-cli"></a>
## `python -m app.hydroclimate.stranded_basins` -- land-locked sub-sea-level pits

A "stranded basin" is an endorheic depression whose floor sits *below sea level* and that has
**no drainage path to the ocean at all**. Such a node is neither hydrology's connectivity-aware
`is_ocean` nor above sea level, so the marine sink, coastal planation, and lake infill all skip
it -- it churns (merge/split) in the event log every step and never drains or fills. The event
log has this today but drowns it in near-sea-level transient-pond spam (see
[Event log](debug-ui.md#event-log)); this surfaces the same thing as one clean list.

The criterion is read straight off this step's already-resolved depression hierarchy
(`hydrology.HydrologyFields.lake_forest`, `lakes.build_lake_hierarchy`): a **top-level** basin
whose `max_depth is None` (lakes.py's own "no known spill to the ocean" state) *and* whose
`floor_elevation` is below `world.sea_level_m`. Roots only -- an endorheic root is the maximal
"no drainage" unit and its floor is the min over every descendant, so a deep sub-basin is
already covered.

Both the CLI and `GET /worlds/{world_id}/stranded_basins` go through
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

Test: [`unit_tests/test_stranded_basins.py`](../../backend/unit_tests/test_stranded_basins.py).

## `python -m app.hydroclimate.lake_hierarchy_diagnostics` -- lake-merge-forest depth/size

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

Test: [`unit_tests/test_lake_hierarchy_diagnostics.py`](../../backend/unit_tests/test_lake_hierarchy_diagnostics.py).

## `python -m app.phase_budget_diagnostics` -- per-phase Hc/Hm budget

Built for [GitHub issue #216](https://github.com/adubey/mantle-bloom/issues/216)
("Investigate long-run Hc/Hm decline"). Every mechanism that touches crustal thickness (Hc) /
mantle-lithosphere thickness (Hm) -- convergent/divergent deformation, arc magmatism, oceanic
cooling relaxation, decompression melting, boundary advance/retreat, gap filling, plate
merges/cleanup, failed rifts, erosion, the end-of-step cap clamp -- is individually instrumented (see [`phase_budget.py`](../../backend/app/phase_budget.py)'s module docstring for
the full mechanism-to-call-site mapping). This tool turns that instrumentation on, steps the
world forward, and reports the accumulated before/after node count and sum(Hc)/sum(Hm) per
phase, split by both plate type (`crust_type`) and per-node type (`crust_type_code`) -- so a
mean-thickness decline can be attributed to a specific mechanism, to changing node population,
or to reclassification, rather than treated as one undifferentiated trend.

For issue #317 the same observations also feed a separate debug-only Hm source/sink ledger
([`hm_ledger.py`](../../backend/app/hm_ledger.py)). Hm is a thermal/mechanical thickness, so this
is an attribution identity, not rock-mass conservation: `live change = sources - sinks +
type reclassification + signed residual`. The JSON report keeps `hm_source_sink_accounts` separate from `phases` to
prevent double booking, reports closure independently for continental and oceanic nodes, and
includes actual per-cell Hm-cap entries/exits under each phase's `hm_cap_transitions`. Gap
filling is bracketed as `gap_fill`; it was the large previously-unattributed oceanic writer.

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

  phase                          calls    d(count)     d(area) km2     d(Hc vol) km3     d(Hm vol) km3
  ---------------------------------------------------------------------------------------------------
  boundary_advance               ...
  boundary_retreat               ...
  convergent_deformation         ...
  ...

continental/oceanic node-type split (per phase, resolved against crust_type_code)
  phase                          cont d(Hc) km3   ocean d(Hc) km3    cont d(Hm) km3   ocean d(Hm) km3
  ---------------------------------------------------------------------------------------------------
  ...
```

**How to read it:**

- **`d(Hc vol)`/`d(Hm vol)`** are this interval's net volume change, summed over every node
  the phase touched with each node weighted by its own accounting area
  (`Plate.accounting_areas_m2`, the exact cell areas), not a mean -- comparable directly across phases regardless of how many nodes each one happened to
  touch. Rows are sorted by `|d(Hc vol)|`. A large `calls` count with a small
  `d(Hc vol)` means the phase ran often but nearly canceled out (e.g. `decompression_melting`
  resetting a column to almost the same reference thickness it already had); a small `calls`
  count with a large delta means a rare but individually large event (a merge).
- **`d(area)`** is the touched slice's change in covered area, from the same per-node areas.
  `d(count)` and `d(area)` can disagree in sign, since a phase can swap large cells for small
  ones.
- **`d(count)`** is the touched slice's own node-count change -- nonzero only for phases that
  can add/remove nodes (`boundary_advance`/`boundary_retreat`, `gap_fill`, cleanup/merge);
  the pure Hc/Hm-mutation phases (`convergent_deformation`,
  `divergent_deformation`, `arc_magmatism`, `oceanic_cooling_relaxation`,
  `decompression_melting`) always report 0 here, so any thickness change they show is a real
  per-node change, not a population effect.
- **the node-type split table** separates a genuine Hc/Hm change (both columns move) from pure
  reclassification (one column drops while the other rises by roughly the same amount, with
  little net change in the top table's combined total) -- `decompression_melting` is the
  phase most likely to show this, since it can flip a node's effective type.
- **`--json`**'s `scopes` dict on each phase row also has `continental_plate`/`oceanic_plate`
  (split by the *owning plate's* `crust_type` rather than each node's own effective type --
  the two disagree wherever a plate carries nodes of the other type, e.g. an accreted terrane)
  and, per scope, the raw `count_*`, unweighted per-node `sum_hc_*`/`sum_hm_*`, area-weighted
  `area_*_m2` and `hc_volume_*_m3`/`hm_volume_*_m3` totals behind the table. The area fields
  are recorded per node as each phase runs (issue #257) -- on a quad world they can't be
  recovered from the counts afterwards.
- **`hm_closure`** in JSON gives the opening and closing live Hm inventories, booked source,
  sink, signed node-type reclassification totals, and residual. `hm_suture_budget` separately reports
  connected donor-front Hm, the volume placed on survivors, the unplaced remainder, and
  bilateral donor→neighbor totals. Those transfer counters are diagnostic only; the live
  change is already booked by `boundary_retreat`.
- This tool answers issue #216's first acceptance-criteria item (instrumentation); it does not
  by itself separate genuine geological sinks from numerical drift -- that's a reading exercise
  against the mechanism list in the issue itself, using this tool's own numbers as the evidence.

Test: [`unit_tests/test_phase_budget_diagnostics.py`](../../backend/unit_tests/test_phase_budget_diagnostics.py).

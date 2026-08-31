# TODO

Tracked follow-up work. Each item: what, why it was deferred, and enough of a starting
point that picking it up doesn't need a fresh investigation.

---

## Diagnostic ("ABL") wind model: close the last ~5-10% gap to the CFD

**Status:** shipped behind `World.wind_model` (Controls window; `"diagnostic"` default,
`"cfd"` opt-in for the full shallow-water solve). See `docs/simulation-model.md#wind-model`.

**Where it stands.** `"diagnostic"` skips the shallow-water solve entirely and rebuilds
wind from `climate.compute_wind` plus air temperature from
`climate.compute_air_temperature_diagnostic` (radiative-equilibrium temperature +
`AIR_TEMP_DIFFUSION_ITERATIONS` gentle Jacobi passes). Measured against a 12-15 step CFD
run at `fluid_density=2.0` across three seeds:

| | vs CFD |
|---|---|
| land-biome agreement | ~84-91% |
| precipitation rel-RMS | ~8-11% |
| air-temperature rel-RMS | ~9-11% |
| step wall-clock | ~6x faster at `fluid_density=2.0` (12 steps: 42s -> 7s); ~1.5x at 0.5, more at 4.0. CFD substep loop itself ~15x. |

The wind *field* barely matters downstream (swapping CFD wind for `compute_wind` while
keeping CFD air temperature holds land-biome agreement at ~95%). Almost the entire residual
gap is **air temperature**: `compute_air_temperature_diagnostic` reproduces the CFD's
near-radiative-equilibrium field to ~9% but not its genuine advective/diffusive structure,
nor the fact that the CFD temperature is still crawling toward equilibrium from its
maritime-moderated bootstrap for the first dozen-plus steps.

**What was tried and rejected.** Adding a single semi-Lagrangian downwind advection of the
equilibrium field along the diagnostic wind (blend 0.4-0.5, 12-18 deg) made agreement
*worse* (~68-73%) at every setting swept -- either the `1/lambda_radiative` characteristic
length is being mis-estimated or the one-shot advection over-corrects. Adding any maritime
moderation (`compute_air_temperature`'s nearest-ocean pull) also made it worse -- the CFD
does not do that.

**Options for closing it, roughly in effort order:**

1. **Quasi-steady linear correction.** Solve the steady advection-diffusion-relaxation
   balance `u.grad(T) = kappa*lap(T) + lambda*(T_eq - T)` approximately -- a handful of
   Jacobi/Gauss-Seidel sweeps of that elliptic operator on the climate grid, seeded from
   `T_eq`. Pure numpy, ~5-15 ms. This is the principled version of the failed one-shot
   advection: iterating the operator instead of applying it once should avoid the
   over-correction. Constants (`kappa`, `lambda`) can start from
   `atmosphere_cfd.TEMPERATURE_DIFFUSIVITY_M2_S` / `RADIATIVE_RELAXATION_PER_S`.

2. **Lightweight prognostic temperature.** Keep a persistent `(H, W)` air-temperature field
   on `World` (or on a slimmed diagnostic state), advected each step by `compute_wind`'s
   output + diffused + relaxed toward `T_eq` -- an explicit temperature-only integrator with
   no gravity-wave CFL, so ~10-30 substeps not ~1000. Recovers the "still relaxing from the
   bootstrap" transient the closed-form version can't. Costs a new persisted field and a
   save/load bump.

3. **Accept it.** 85-91% land biomes at 1/10th-1/25th the cost is already the intended
   trade. Document the gap in the Controls tooltip (done) and move on.

**Validation harness.** `test_diagnostic_wind_model_tracks_the_cfd_biome_map`
(`backend/stress_tests/test_climate.py`) already steps one world ~10 times each way and
asserts land-biome agreement > 0.78 and precipitation rel-RMS < 0.25 -- a regression floor,
not a tight bound. If this item is picked up, tighten those numbers as the approximation
improves, and add per-seed / per-density coverage (the agreement is config-sensitive:
~72% at `climate_density=4.0 / fluid_density=1.0`, ~84-89% at 2.0/2.0).

**Tuning note.** `AIR_TEMP_DIFFUSION_ITERATIONS` / `AIR_TEMP_DIFFUSION_WEIGHT` (climate.py)
were picked from a 3-seed sweep at 15 steps; 6 passes at weight 0.5 was the best of
{raw, 3, 6, 8, 12}. Worth re-sweeping if the land-temperature mapping
(`LAND_TEMP_MIN_C`/`RANGE`) or the biome thresholds change.

**Also:** `_advance_fluid_dynamics`'s `compute_climate(skip_moisture=True)` forcing snapshot
is still built every step in `"cfd"` mode (~20-60 ms) -- unrelated to this flag, but the
cheapest remaining CFD-path saving if that mode's per-step cost is ever revisited.

---

## Plate geometry degrades on long runs: pole winding, unbounded overlap, bad split siblings

**Status:** the three named sub-bugs are **fixed** 2026-08-30 (details below), but a
long-run re-verify (2026-08-31, `~/Downloads/mantle-bloom-seed888151728-85000000y.mbworld`,
851 steps / 85.1 My, `node_density=4`) shows the *degradation is not gone* -- it has just
shifted from "pole winding + winding-driven node blowup" to a cluster of related failures
below. This save is the same seed the coastal-speckle work used, stepped far longer, and its
geometry is visibly bad in the Plate Inspector. **Open follow-ups, most impactful first:**

1. **Every oceanic plate rails at `MAX_PLATE_RATE` (15 cm/yr).** 6 of 7 oceanic plates sit
   at *exactly* `mantle.MAX_PLATE_RATE` at 85 My (the 7th at 11.7); all 6 continental plates
   are a healthy 0.7-3.9 cm/yr. Stepping the loaded save 4 more times with current `main`
   does **not** relax them -- so this is not just a stale pre-fix save. Bug 3's backward-Euler
   basal-drag fix genuinely settled the *continental* plates but the *oceanic* ones are still
   pinned: `torque`'s explicit slab-pull (+ ridge-push + collision) alone exceeds the clamp
   ceiling for them, every step. Either slab-pull's scale (`SUBDUCTION_LSINK_M`,
   `slab_pull_torque`) is hot, or a mostly-ocean world genuinely subducts oceanic crust on
   nearly every margin and 15 cm/yr is simply the cap it deserves -- but 6/7 at the *exact*
   ceiling forever, driving fast oceanic plates into every continent, is the upstream cause
   of most of what follows. Check `integrate_omega`'s output against the clamp for an
   oceanic plate on this seed and see how far over it lands.

2. **Over-stretched continental plates that have mostly drowned.** Plate 3 (id 3):
   continental, 22,936 nodes, **71% of its own nodes at/below sea level**, median elevation
   **-3113 m**, bounding ellipse ~20,500 x 9,800 km (half the planet), spanning lon -92..+95
   / lat -67..+18. Its interior is pinned at *exactly* `ABYSSAL_REFERENCE_DEPTH_M`
   (-5222.22 m) -- `bathymetry._subside_offshore_continental_crust` correctly draws crust
   >1400 km from its own land down to abyssal depth, so once a continental plate is stretched
   to ~2x its natural size the model *correctly* oceanises the middle and you get a giant
   "continental" plate that is 80% deep ocean. Plate 17 (its age-43 split sibling) is the
   same story at 63% submerged / -87 m median. The stretch itself is the bug: `_grow_or_
   shrink_line_for_deform` keeps extending row ends into every newly-opened gap as the plate
   shears (euler pole for plate 3 sits at 24 N, 139 W -- far from the plate, so its rotation
   is nearly pure shear across the plate body), and nothing caps a plate's total area
   against its crustal volume.

3. **The "streaking" in the Plate Inspector / platesDetail view** is plate 3's (and a
   neighbour's) territorial boundary running as a long, straight, regular-sawtooth diagonal
   chord across ~140 deg of open ocean -- the staircase edge of a plate whose rows each end
   one `dtheta` step further than the row below (a thin triangular tongue grown by repeated
   row-end extension). In the raw `elevation` view the same feature is only a faint
   shallower-ocean smear; it is mostly a boundary-polyline artifact of an over-extended
   lattice, not a separate elevation bug. Node *spacing* within plate 3 is a clean 1.00x target both
   along-row and across-row -- the lattice is regular, just far too large.

4. **Stalled continental-continental overlaps that never heal.** Plate 4 has **17% of its
   nodes sitting on top of plate 2** (fresh, age-9, 29,778-node plate 2 -- likely a bad
   split/merge partition), and this overlap is *stable* across 4 stepped My. It is not even
   in `world.collision_progress` (the closing-rate test doesn't see it as converging), so
   the sustained-collision -> merge path will never fire; continental crust never retreats
   (`shrinkable_all` is all-False for a continental self-plate). Plates 3/6 overlap 6% and
   *are* tracked -- but `collision_progress[(3,6)]` is only 30.7 My against a 50-100 My merge
   threshold, so they will keep overlapping for ~20 My more. `_merge_probability` also drops
   toward `MERGE_PROBABILITY_FLOOR` (0.02) once a pair's combined node share passes
   `MERGE_SIZE_UNLIKELY_FRACTION` (0.25) -- plate 3 alone is ~16% of the world's nodes, so
   most of its collision pairs are effectively unmergeable and just overlap indefinitely.

5. **Node-count blowup persists** (the original headline symptom): ~140k nodes at 85 My for
   `node_density=4` vs a clean-tiling estimate of ~130k *at 1x* -- consistent with the
   ~15-75%-over range the 2026-08-30 table recorded, i.e. not fixed, just no longer
   dominated by the winding rows.

The **k-means split-cluster quality** noted under bug 2 / bug 3 (velocity-space clusters cut
from spatially-intermingled points -> disjoint daughters that drift back over each other,
with euler poles fit far from the daughter body) is the common root of (2), (3) and (4) and
is still entirely open.

Bug 1 (pole winding) **fixed** 2026-08-30 -- wrap guard in
`plates._grow_or_shrink_line_for_deform` + pole-cap margin in `_claim_adjacent_territory` +
`elevation_lines.regularize_line` unwinding already-wound rows on load. Verified against the
158.6 My save: 191 over-wound rows -> 0 after one step, `theta` max span 10,489 deg -> 360
deg, ~229k nodes -> ~223k.

Bug 2 (split/defrag siblings) **fixed** 2026-08-30 -- the partition, not the re-growth, was
the fault: a great-circle cut can slice a row (one small circle of local latitude) so that
one daughter is left holding *two* arcs with the other daughter's territory in the gap
between them, and `outline_world` / the `contains_batch` row-lookup fast path / `regularize_line`
all assume each line is a single contiguous arc from `theta[0]` to `theta[-1]` -- so that
daughter's envelope claimed the gap (every sibling node in it read as contested), and the
next `regularize_line` pass resampled straight across the gap, growing fresh nodes through
the sibling. Fix: `elevation_lines.largest_contiguous_run` reduces every masked row to its
longest gap-free arc; wired into `PlateWithLines.split` / `LithospherePlate.split` /
`_plates_from_node_masks` (partition time) and `regularize_line` (so worlds saved before this
self-heal instead of getting worse each pass). Dropped slivers along the cut are re-grown by
ordinary gap-fill/deform. Not addressed: the quality of the velocity-space k-means clusters
`maybe_split_plate` cuts on (poles fit from spatially-intermingled clusters can still send
two disjoint daughters drifting back over each other) -- that overlaps bug 3.

Bug 3 (every plate at `MAX_PLATE_RATE`) **fixed** 2026-08-30 -- the TODO's original guess
(v1 `fit_euler_pole` over-fitting) was stale: the running engine is the v2 torque balance,
and `torque.integrate_omega` was doing a plain *explicit* Euler step on the basal-drag term,
which is stiff enough (the asthenosphere coupling relaxes a plate toward the local
mantle-flow rate in far less than one 100-ky step) that the step overshoots the true
solution by ~19 orders of magnitude every call -- `mantle.clamp_rate` then pins the result
at `MAX_PLATE_RATE` on step 1 and forever after, for every plate, from fresh gen on. Fix:
`torque.basal_drag_coefficients` returns Eq. 10's drag as its exact affine-in-omega split
`tau(omega) = b - K @ omega` (`b`/`K` are cheap 3-vector / 3x3 sums over the plate's nodes),
and `integrate_omega` now solves the backward-Euler system `(I + g K) omega_new = I omega_old
+ g (tau_explicit + b)` -- unconditionally stable, drag implicit, every other (bounded,
geometry-driven) torque still explicit. Drag-only plates now settle at their local
mantle-flow rate (~1-5 cm/yr for seed 936513024); plates still rail at `MAX` only when slab
pull etc. genuinely demand it (~0-2 of 8 across sample seeds, vs 8/8 before). Not addressed:
the v1-era plate-motion section of `docs/simulation-model.md` still describes the retired
`fit_euler_pole`/`VELOCITY_DAMPING` model and needs a full rewrite for the torque engine
(bigger than this bug).

Diagnosed 2026-08-30 from two save files of seed 936513024
(`~/Downloads/mantle-bloom-seed936513024-90000000y.mbworld` and `...-158600000y.mbworld`,
1586 steps of 100 ky). Symptoms show up in the Plate Inspector: plates with concentric
circles (33, 36, 38), plates whose points overlap a neighbour's (14/22, and far worse
33/38, 36/39, 22/35), plates with holes (36), combs of stranded one-node "teeth" (11).

**What's measured.**

| | fresh gen | 90 My (900 steps) | 158.6 My (1586 steps) |
|---|---|---|---|
| plates | 16 | 16 | 26 |
| total nodes (4x density; clean tiling ~130k) | 32,651 @ 1x | ~149k (~15% over) | ~229k (~75% over) |
| plates whose territory reaches a local pole (\|phi\| > 88 deg) | 0 | 0 | 3 |
| elevation-lines winding past 360 deg of theta | 0 | 0 | 191 (one row spans ~10,500 deg) |
| plates pinned at exactly `MAX_PLATE_RATE` (bug 3, now fixed) | 16/16 | 16/16 | 26/26 |

The world is healthy at 90 My (the ~15% overlap is the documented bounded envelope/
randomized-order effect). Nearly all the damage happens in the 90 -> 159 My window.

**Three distinct bugs, most impactful first.**

1. **Pole winding / no periodic-theta guard (causes the concentric circles, the 36-style
   holes, and most of the node-count blowup).**
   - `PlateWithLines._claim_adjacent_territory` (`plates.py` ~L1617) adds new phi rows
     outward all the way to `max_phi_limit = pi/2 - spacing/2` -- it will march a plate
     right onto its own local pole whenever the space is open.
   - `_grow_or_shrink_line_for_deform` (`plates.py` ~L1523) then extends those near-pole
     rows with `dtheta = spacing_rad / max(np.cos(line.phi), 1e-3)`. Near the pole
     `cos(phi) -> 0`, floored at `1e-3`, so **each inserted node jumps ~9.8 rad (~561 deg)
     in theta**. Nothing checks whether the row has already closed a 2*pi loop, and the
     open-space test (`_count_open_prefix`, world-space) barely moves near the pole, so the
     ring keeps winding indefinitely. That's the 10,500-deg row and plate 38 at ~210
     nodes/row.
   - Smaller time steps do **not** help: they only drop `n_distance_cap` to 1 (one node/
     step instead of a few), but that node still jumps ~561 deg and there are ~10x more
     steps. Confirmed this is geometric + threshold, not integration-step-size.
   - **Fixed** as: (a) `_claim_adjacent_territory` keeps `POLE_CAP_MARGIN_MULT` (4) target
     spacings clear of `+-pi/2`, leaving a small unclaimed polar cap; (b) end-growth caps a
     row's span at one revolution (`_ROW_FULL_REVOLUTION_SLACK`), and the end stops once the
     loop closes; (c) `regularize_line` / `needs_regularizing` unwind an already-over-wound
     row to its outermost single revolution, so old saves self-heal. Not done: closing the
     ring's two ends into a genuine periodic loop (they just stop a `dtheta` apart), and
     pulling *generation*'s lattice back from the pole too (left as-is -- a plate that owns
     its pole at generation keeps its small rings; only growth toward the pole is capped).

2. **Split/defragmentation produces overlapping siblings.** ~~Plates 33/38, 36/39, 22/35
   overlap 46-93% by node count.~~ **Fixed 2026-08-30** -- see the Status block above. The
   partition itself was leaving a row as two arcs across the sibling's territory, which every
   consumer of an `ElevationLine` (all of which assume one contiguous arc per row) then read
   as claimed. `largest_contiguous_run` in `elevation_lines.py` keeps every partitioned row a
   single arc.

3. **Every plate railed at `MAX_PLATE_RATE` (15 cm/yr).** ~~True at *both* epochs, all
   plates.~~ **Fixed 2026-08-30** -- see the Status block above. Not `fit_euler_pole`
   (v1, retired): `torque.integrate_omega` was stepping the stiff basal-drag term with plain
   explicit Euler, overshooting by ~19 orders of magnitude every call so `clamp_rate` pinned
   every plate at the ceiling from step 1. `basal_drag_coefficients` + a backward-Euler solve
   for the drag term in `integrate_omega` fixes it; plates now spread across a realistic
   1-15 cm/yr and only genuinely slab-pull-driven ones reach `MAX`. This also takes the heat
   off the bounded-overlap mechanism that (1) and (2) were overwhelming.

**Not bugs:** plate 23's "two blobs" is a single connected node cloud (1 component even at
2.5x spacing) -- a crescent plate whose min-area enclosing ellipse balloons in the Inspector.
Plate 11's "teeth" are 8 legitimate one-to-two-node stub rows from heavy subduction
(`deform()` never deletes a line's last node); defragmentation is meant to prune them but
runs only every `DEFRAG_INTERVAL_STEPS`. Related: the severed-lobe defrag work on
`fix/plate-defragmentation` (unmerged).

---

## Speckled low-relief coastlines: a drowned flat shelf dithers pixel-by-pixel across sea level

**Status: PARTIALLY RESOLVED.** The coastal planation + infill feedback (option 2, "What
landed" below) cut the density-1 dither roughly in half; barrier islands (option 3) fall out
of it in emergent form. But a fresh look at
`~/Downloads/mantle-bloom-seed888151728-85000000y.mbworld` (2026-08-31, seed 888151728,
851 steps / 85.1 My, `node_density=4`) shows the coast is **still visibly speckled** -- a
broad checkerboard fringe along the south and east coasts of the plate-3 land mass around
5 N / -10 E. New diagnosis below; the feedback is running but is being *swamped*.

**Why it's still speckled at 85 My (measured, per 100-ky step, in an ~8x30 deg box over the
worst stretch).** The user's hypothesis was right that this should be a coastal plain -- but
the problem is *not* a sediment shortage:

- Sediment supply is **abundant**: `sediment_deposited` median **11.8 m/step**, p90 53 m,
  max **256 m**. Rain 1967 mm, flow_accum up to ~96,000.
- It is deposited **in lumps**: the top 10% of near-sea-level nodes in the box receive
  **51%** of all deposition; 82 of 768 band nodes get >50 m in a *single* step while the
  median node in the same band gets ~0. `route_downstream` drops `DEPOSITION_FRACTION` of
  everything passing through onto whichever discretised flow-graph nodes it routes through
  (`DEPOSITION_MIN_FLOW_M` is only 0.05, so ~half the band counts as a "slow depositing
  river"), and `_spread_beach_sediment` / `_spread_coastal_infill` concentrate on a few
  weight-favoured nodes rather than spreading evenly across the shelf.
- Net elevation change is then **+-50 m per step** (p10/p90 -46 / +51 m, max +237) around a
  ~0 m mean -- right where the surface sits within +-25 m of sea level. **315 of 768 band
  nodes flip land<->ocean every single step.** That is the checkerboard, in motion.
- Nothing removes the lumps: the ground is flat, so `river`/`rain` erosion is ~0 by design
  (relief-gated).
- `coastal_planation_amount` **is not firing** where it's needed: median planation in the
  band is **0.0 m** (only 31% of band nodes nonzero) even though `coastal_openness` there
  (median 0.44) is well above `PLANATION_EXPOSURE_REF`. Planation only touches land within
  `PLANATION_BAND_M` (60 m) *above* sea level and is then gated hard by the prominence
  reweight, so it ignores everything that is already a hair below the waterline (that's
  infill's job) and everything that isn't standing proud of its neighbours.

**Fix direction (not yet done).** The deposition needs to be *rate-limited or diffused* in
the near-sea-level band -- cap per-node per-step deposition at some fraction of local relief,
or spread each routed lump across its k nearest coastal neighbours (mass-conserving) instead
of dropping it on one node -- and planation + infill should engage as **one symmetric pass**
across a band that straddles sea level (grind the +30 m bumps and fill the -30 m ponds
together toward a common local datum) rather than two separately-gated passes that each miss
half the checkerboard. The transient-lake log spam and the stranded deep basins (below) are
both downstream symptoms of the same dithering shelf.

**Related symptom -- event-log flooding.** The world's event log is almost entirely
"N-node lake formed/split ... at elevation ~0 m" -- hundreds of these per My, one pair per
dithering puddle every step. `world.log_event` / `hydrology.lake_events` should either
dedupe/aggregate near-sea-level transient lakes ("~40 coastal ponds churned this step")
or suppress lakes whose floor is within a few metres of sea level and whose lifetime is one
step. As-is the log is unusable for spotting real basin events (there are genuine ones
buried in it: a persistent ~435-node lake oscillating around -1770 m).

**Related symptom -- stranded sub-sea-level basins.** The log also shows persistent
endorheic depressions well below sea level (a 435-node basin at ~-1770 m, another cluster at
~-4560 m) that merge/split every step and never drain or fill. These are the "land-locked
coastal pit" the loose-ends note already flags -- an isolated sub-sea-level node ringed by
land is neither `hydro.is_ocean` nor above sea level, so planation, the ocean sink, *and*
infill all skip it. Worth a dedicated "interior basin below sea level" infill/relaxation
term (or letting lake siltation actually keep up with them).

**Original resolution (density-1) retained below.** Investigated + validated from the same
save at the time.

Driven against that save, 25 steps (2.5 My), feedback on vs off:

| in the 7-17 N / 3-15 E box | feedback on | feedback off |
|---|---|---|
| isolated near-sea-level nodes | 21 -> 8 | 21 -> 19 |
| world-wide speckle fraction | 6.6% -> 4.6% | 6.6% -> 6.1% |

Visually: the sheltered north coast and the south spit resolve to a clean shoreline; the
broad central promontory converts to a coherent marsh/intertidal mosaic rather than a
checkerboard (arguably correct -- it *is* a drowned coastal plain). Pushing the coefficients
harder (`PLANATION_RATE` 400-600, sharper prominence/shelter) made the dither *worse* --
amplifying node-scale openness noise into +-10 m elevation swings -- so the mild committed
values are the right operating point. The old render despeckle hack was measured to be a
near-no-op on top of the feedback (before/after renders near-identical with and without it)
and was removed.

**Loose ends.** Coefficients are from-scratch starting points, worth a sweep if the
low-relief coast ever looks wrong again. A land-locked-coastal-pit infill sink (an isolated
sub-sea-level node ringed by land, which is neither `is_ocean` nor above sea level, so both
planation and the ocean-sink infill skip it -- it lingers as a stranded deep speck) was
prototyped and reverted: it didn't clearly help the metric and the planation/infill
asymmetry it exposed needs more thought. No automated test reproduces the density-4
checkerboard (only a stability regression floor landed) -- see the Tests note below.

**Symptom.** Around 10-13 deg N, 6-12 deg E the Elevation / Biome / Combined maps show
isolated single-pixel islands and single-pixel ponds strung along the coast -- a
checkerboard, not a coastline.

**What's actually there.** That whole strip is one continental plate (id 17) -- a drowned
continental shelf / coastal plain, *not* an ocean-continent boundary. In the transition
band the node elevations sit right on the waterline: median -4 m, ~55% of nodes within
+-25 m of sea level, ~86% within +-50 m, full range about -57..+111 m. The per-node
elevation *noise* (deform deltas, erosion/deposition, isostasy relaxation, generation
noise) is larger than the surface's height above/below sea level, so neighbouring nodes
flip land<->ocean. At render resolution (~0.56 deg/cell, `min(GRID_SPACING_RAD,
line_spacing_rad(node_density))`) each ~62 km node becomes a 1-2 px cell, so the dither
renders as literal single pixels. In a 5-13 deg E / 9-15 deg N sample box, 8 of 92 land
cells and 9 of 73 ocean cells were (near-)isolated. It is genuine node-cloud data (the
Plate Inspector shows the same dither), not purely a render artefact -- but the nearest-node
resample with no coastal cleanup is what turns a fuzzy zone into a checkerboard.

**Why the model makes it and never heals it.**

- `bathymetry.shape_initial_bathymetry` (the margin-grading pass) runs **once at
  generation** and is explicitly coast-guarded (`COAST_GUARD_DEPTH_M = 60 m`), so it leaves
  the shallow shelf hovering at 0 and never revisits it.
- Weathering is relief-gated (`WEATHERING_RELIEF_REFERENCE_SLOPE`, `erosion.py`) -- nearly
  zero on flat terrain, by design, to stop coastal plains drowning -- so the small positive
  land bumps are stable.
- `erosion.coastal_erosion_amount` is symmetric across the shoreline and *not* relief-gated:
  it eats the shallow-ocean nodes about as fast as the low-land nodes, and
  `_spread_marine_sediment` carries the spoil to *lower / deeper* ocean nodes. In a
  low-energy, low-relief embayment this *deepens* the gaps and moves the fill offshore --
  the opposite of what should happen.
- `_spread_beach_sediment` only redistributes what flow routing actually delivers to the
  coast; this dry, flat strip has almost no flow accumulation, so there is nothing to
  prograde with.
- Nothing anywhere looked at coastal *connectivity* ("this ocean cell is nearly landlocked",
  "this land cell is nearly surrounded by water") -- no planation, progradation, spit,
  barrier, or lagoon logic. (Fixed by "What landed" below: `_coastal_openness` is exactly
  that connectivity measure.)

Net (before the fix): a marginally-submerged flat sheet is a stable fixed point that just
dithers forever.

**What landed (option 2 + emergent option 3), 2026-08-31.** A per-step pass in
`erosion.apply_erosion`, all mass-conserving via `np.add.at` (no new `_flatten`-style
term):

- `_coastal_openness(points, is_ocean)` -- a "wave exposure" proxy in [0, 1]: the fraction
  of nodes within `COASTAL_OPENNESS_RANGE_KM` (~150 km fetch scale) that are open ocean
  (hydrology's connectivity-aware `is_ocean`, so inland-lake shores read as fully enclosed).
  Two `cKDTree.query_ball_point(..., return_length=True)` radius counts, density-independent.
- `coastal_planation_amount(...)` -- land within `PLANATION_BAND_M` of sea level is ground
  down toward a wave-cut platform sitting `PLANATION_UNDERCUT_M * exposure` *below* sea
  level (so a genuinely exposed low sheet is cut into open water, not left dithering on the
  waterline), rate `PLANATION_RATE_M_PER_MYR * exposure * proximity * prominence`.
- `_spread_coastal_infill(...)` -- the "shallow + sheltered sink with priority" the plan
  called for: planed rock + `COASTAL_INFILL_MARINE_FRACTION` of the submarine/coastal pool
  (the rest still goes to `_spread_marine_sediment`'s downhill-to-deep spread) is spread
  onto sheltered shallow ocean, weighted by `shelter * proximity * hollow * headroom`, and
  a sheltered sink fills to `INFILL_MARSH_CREST_M * shelter` *above* sea level -- a silted
  embayment emerges as marsh (`biomes.classify_wetland`).
- **Barrier islands** are emergent: a shallow sink with land within `BARRIER_LANDWARD_KM`
  that still faces open water (`openness >= BARRIER_MIN_OPENNESS`) gets a `BARRIER_PRIORITY`
  weight boost and a cap raised `BARRIER_CREST_M` above sea level, so a shore-parallel bar
  breaches; the water it then encloses loses open-ocean neighbours, so next step its own
  openness falls and the back-barrier lagoon silts up. No explicit lagoon flag, no
  longshore-direction field.
- **Prominence** (`PROMINENCE_REF_M` / `PROMINENCE_MAX`, from each node's height above its
  flow-graph neighbourhood mean): waves plane protrusions and fill hollows. Only reweights
  planation/infill; the openness field alone, at a ~150 km fetch scale, is too smooth to
  resolve a land/ocean call at ~60 km node spacing.

Tests: `unit_tests/test_erosion.py` (`_coastal_openness`, `coastal_planation_amount`,
`_spread_coastal_infill` -- direction + exact mass conservation), and
`stress_tests/test_world_stepping.py::test_coastal_feedback_stays_stable_over_many_steps`
(a regression floor only -- a fast, faithful reproduction of the density-4 checkerboard for
a stress test proved infeasible: injecting a `+-40 m` node checkerboard triggers unrelated
elevation instability, and a sea-level jump on a density-1 world just makes a rough newborn
coast whose transient roughening swamps the feedback's slow ~My effect).

**Options considered, in effort order.**

1. **Render-only cleanup (cosmetic).** A majority / morphological filter on the render's
   land-ocean field flips isolated 1-cell specks to match their surroundings -- zero physics
   risk, but the simulation still holds the dithering shelf. Shipped 2026-08-31 as
   `render_image._despeckle_coastal_elevation`, then **deleted** once option 2 landed (it
   was measured to be near-redundant on top of the feedback).

2. **Coastal planation + infill feedback (the real fix)** -- **shipped 2026-08-31**, see
   "What landed" and the Status table above. Cuts the node-level dither roughly in half and
   turns the drowned shelf into a marsh/intertidal coast instead of a checkerboard.

3. **Barrier islands (the user's framing)** -- **delivered in emergent form** (see "What
   landed"). A genuine wave-exposure / fetch field plus an explicit longshore-transport
   direction (for real spits and drift-aligned bars, not just fetch-sheltered accretion)
   remains an optional future refinement.

---

## Diagnostic views & debug output (from the 2026-08-31 seed-888151728 investigation)

Working through the plate-geometry and coastal-speckle degradation above needed several
numbers the program didn't surface. What landed, and what's still worth building:

**Landed 2026-08-31 -- Plate Inspector motion / shape / overlap fields.** `GET /world/plates`
(`main._plate_summary`) and the Plate Inspector panel (`App.tsx`) now report, per plate:
`speed_cm_per_yr` + `at_max_rate` (railed-at-`MAX_PLATE_RATE` flag, shown in red),
`euler_pole` (lat/lon), `age_steps`, `median_elevation_m` + `submerged_fraction` (red when a
continental plate is >50% submerged), `overlaps` (which other plates this one's territory
sits on top of, and by what fraction of its own nodes -- `main._plate_overlaps`, one global
cKDTree pair query), and `collisions` (`world.collision_progress` timers involving the
plate). `mantle.rad_per_yr_to_cm_per_yr` is the new unit helper. Test:
`test_plates_endpoint_reports_motion_shape_and_overlap_diagnostics`. These three numbers --
"every oceanic plate is at 15.0 cm/yr", "plate 3 is 71% submerged", "plate 4 is 17% inside
plate 2" -- are what turned a vague "the plates look wrong" into the specific follow-ups in
the plate-geometry section.

**Still worth building:**

1. **A speckle / coastal-dither overlay render mode.** Colour every node whose elevation is
   within a threshold of sea level by the fraction of its k nearest neighbours on the
   *opposite* side of the waterline (the exact metric the investigation scripts compute:
   `near = |elev - sea_level| < 120 m`; `frac = mean(neighbour is opposite class)`; flag
   `>= 0.75`). Instantly shows where the coast is a checkerboard vs a clean shoreline, and
   makes a before/after for any coastal-feedback change legible without an ad-hoc script.
   Natural home: a `render_image.py` view alongside `platesDetail`, or a boolean overlay on
   the `elevation` / `biome` views.

2. **A per-node geomorph-rate view.** Render `ErosionResult.sediment_deposited` (and/or net
   `dElev` this step) as a diverging map. The lumpiness of deposition in the near-sea-level
   band -- 256 m on one node, ~0 on its neighbour -- is invisible in every current view but
   is the whole coastal-speckle mechanism. `ErosionResult` already carries the arrays;
   this is a render path plus maybe a `/world/sample_at` field.

3. **Stranded-basin report.** A `/world/lakes`-style endpoint (or a field on it) listing
   endorheic basins whose floor is below sea level and that are not connected to the ocean:
   node count, floor elevation, centroid, how long they've persisted. The event log has this
   information today but drowned in transient-coastal-pond spam (see the coastal section).

4. **Event-log dedup / severity for lake churn.** `world.log_event` should collapse the
   hundreds of per-step "N-node lake formed/split at ~0 m" messages into one aggregate line,
   or drop one-step lakes with a near-sea-level floor entirely, so the log is usable for
   real basin/tectonic events again.

5. **A standalone `python -m app.<something> <save.mbworld>` plate-diagnostics dump.** The
   investigation used a throwaway script (scratchpad `probe*.py`) to load a `.mbworld` and
   print the per-plate table, overlap matrix, `collision_progress`, and node-count vs
   clean-tiling estimate. A checked-in version -- reusing `main._plate_overlaps` and
   `_plate_summary` -- would make the "is this save's geometry healthy?" check a one-liner
   for the long-run re-verify the plate-geometry section keeps asking for.

---

## Deferred work found in code/doc comments

A sweep of every source file (2026-08-27) turned up **no `TODO`/`FIXME` comments** -- this
codebase tracks follow-ups here, not inline. The items below are the loose ends the comments
*do* describe, collected so picking one up doesn't need a fresh grep.

### `PlateWithRTree` line regularization

**Where:** `plates.py` (`PlateWithRTree`), noted in `docs/architecture.md` (plates.py entry:
"inline line regularization (`PlateWithRTree`'s own versions are still a TODO)").

`PlateWithLines.deform()` does inline line growth/shrinkage and regularization per turn.
`PlateWithRTree` -- the R-tree-backed variant -- doesn't carry equivalent regularization
logic. If `PlateWithRTree` is meant to become a drop-in replacement, it needs its own
version of the per-turn node density / spacing upkeep (`elevation_lines.py`,
`TARGET_LINE_SPACING_RAD`), otherwise its lines drift out of spec over a long run.

### `bvh.py` tree-vs-tree traversal is built but not wired into the sim

**Where:** `torque.py` ~line 100 ("Not yet exercised at runtime; kept available and
independently tested").

`bvh.py`'s `query_nearest_cross` (tree-vs-tree) is validated against brute force
(`unit_tests/v2/test_bvh.py`) but nothing in the running simulation calls it. The per-step
per-plate nearest-neighbour query in `torque.py` deliberately uses `scipy` `cKDTree` instead
(compiled batch query beats the pure-Python BVH recursion badly at ~16k nodes/plate -- ~2
min vs a couple of seconds for an 8-step run). The intended call site for the BVH version is
a smaller, less frequent query -- `merge_split.py`'s collision-pair proximity check is
called out as "the natural next call site." Either wire it in there, or drop `bvh.py`'s
cross-traversal path if it's never going to be used.

### `World.volcanic_field_plate_ids` is dead state

**Where:** `world.py` ~line 70 ("Nothing populates this set any more ... kept for now").

Since `PlateWithLines.deform()` started spawning overstretched-rift volcanoes as new nodes
on the plate's own line (rather than as a separately tracked volcanic-field plate), nothing
adds to `volcanic_field_plate_ids`. Per-node eruption rolling in `volcanism.py` reads
`is_volcano` directly and doesn't need it. It's retained only as a place to report a field
"cooling" if per-field tracking ever comes back. Decide: revive the tracking, or remove the
field (and bump save/load).

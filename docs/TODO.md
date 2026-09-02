# TODO

Tracked follow-up work. Each item: what, why it was deferred, and enough of a starting
point that picking it up doesn't need a fresh investigation.

---

## Lake-hierarchy caterpillar trees: lean on siltation to keep them shallow

**Context.** On an old, rough world (seed 23097282 @ 79.2 My) `lakes.build_lake_hierarchy`
produced a merge forest whose deepest subtree was ~3,500 levels -- a near-linear chain, one
`Lake` per catchment in a long spill cascade. `lakes._resolve` was recursive and overran
Python's recursion limit; **that crash is fixed** (2026-09-01) by making `_resolve` an
explicit-stack post-order walk, matching `build_lake_hierarchy` / `_catchment_roots` which
were already iterative for the same reason. So this is no longer a crash, just a smell.

**Why the tree gets that deep.** `build_lake_hierarchy` unions catchments pairwise in
ascending saddle-elevation order; a chain of many small closed basins each spilling into the
next hangs one new component off the growing blob per merge, giving a depth ~= the number of
catchments in that drainage network. Thousands of tiny sub-resolution depressions is itself
the pathology -- the same dithering-shelf / stranded-basin family already tracked below.

**Direction.** Lake sediment deposition (`_water_balance`'s `SILT_ACCUMULATION_COEFFICIENT`
term, folded into `elevation` by `erosion.py`) should be filling these pits in faster than
tectonic roughening digs them, collapsing the cascade back toward a handful of real basins.
It clearly isn't keeping up. Look at: whether the silt term actually reaches the shallow
transient pits (it only deposits under standing water this step, `elevation < new_level`, so
a pit that never holds water gets nothing), the coefficient's magnitude vs. the per-step
roughening rate, and whether a cheap "fill depressions below N nodes / below M metres of
relief straight into `elevation`" pre-pass belongs in erosion or terrain relaxation. Shares
root cause with "stranded sub-sea-level basins" and the coastal-dither work below.

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

1. **Every oceanic plate rails at `MAX_PLATE_RATE` (15 cm/yr).** **FIXED 2026-09-02 (seed
   495717634 @ 254.8 My -- 31/31 oceanic plates at *exactly* the clamp).** Two causes, both
   in `torque.py`:
   - `slab_pull_torque` applied pull to the plate's *entire* near-neighbour band
     (`dist_to_neighbor <= reach_rad`), ridge and transform stretches included -- its own
     docstring said "contested" but the code never checked. Now gated by
     `subducting_boundary_mask` (oceanic self, in-band, closing rate past
     `TRANSFORM_RATE_THRESHOLD`).
   - Spec Eq. 8 has slab *pull* with no matching resistance. A real slab's descent is resisted
     mostly by viscous coupling to the surrounding mantle, not by basal drag on the trailing
     surface plate -- so `integrate_omega`'s steady state wanted 16-65 cm/yr for every oceanic
     plate and `clamp_rate` pinned them all at 15. New `slab_drag_coefficient_matrix`
     (`SLAB_MANTLE_VISCOSITY_PA_S = 1e21` over `SUBDUCTION_LSINK_M`, folded into the implicit
     `K` like basal drag) makes oceanic speed self-regulate: median ~4-5 cm/yr on this seed,
     31/31-railed -> 0-3 railed, only genuinely fast slab-pull-driven plates near the cap.
   `test_torque.py::test_subducting_mask_only_flags_a_converging_oceanic_boundary` /
   `test_slab_drag_keeps_a_subducting_plate_off_the_clamp` pin both. This removes the upstream
   driver of most of what follows in this list.

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

   **Partly addressed 2026-09-01** (seed 559394024 @ 199 My). `torque.classify_boundary_nodes`
   used to polygon-test only nodes within `reach_rad` (~3 spacings) of a neighbour node, so a
   plate slid *deep* over another had its deep-interior overlapping nodes classified as
   neither contested nor divergent -- `rheology` never touched them, so the overlap just sat
   with no forcing at all. It now also tests any node inside a neighbour's bounding sphere: a
   deep continental overlap classifies `contested` -> `apply_convergent_deformation` thickens
   Hc/Hm -> mountain uplift (the overlap crumples in place, as a real collision should), and a
   deep oceanic overlap classifies `contested` -> subduction deletion -> the overlap actually
   heals. Verified on the save: the plate-21/0 15%-overlap cleared within ~15 steps and
   plate-21's overridden nodes rose ~+400 m. **(2026-09 correction: the continental half of
   that -- `apply_convergent_deformation` thickening Hc -> uplift -- did nothing at the time,
   because the Mohr-Coulomb yield scale was mis-calibrated and no collision thickened crust
   at any speed. The "+400 m" was the oceanic subduction-deletion path. Calibration fixed --
   see "Orogenic renewal was completely dead" above -- so the continental case now genuinely
   crumples-and-thickens.)** **Still open:** the `_merge_probability` size
   floor for a genuinely huge pair (direction 3 below) -- though the split tuning (see
   "Plate count only decreases" below) now keeps a single plate from reaching ~16% of the
   world in the first place. `ElevationLine.overlap_onset_years` +
   `merge_split.update_overlap_tracking` + the `overlapAge` debug view now record *which*
   nodes have been overlapping and *since when*, so a genuinely stuck overlap is legible
   rather than a bare current fraction.

   **Further addressed 2026-09-02** (seed 656865324 @ 60 My): continental contested ends now
   *retreat* against a continental neighbour too (direction 1 above), so a stalled suture
   consumes its own territory overlap geologically rather than waiting on the forced-merge
   timer (direction 3). Plate 8's 17%-since-33-My overlap on plate 6 -- below the 30%
   forced-merge threshold, so `overlap_progress` never fired -- drains to ~2% within ~1.5 My.
   The retreated shortening feeds extra Hc thickening (`CONTINENTAL_COLLISION_SHORTENING_
   BOOST`) so the overlap crumples into an orogen.

5. **Node-count blowup persists** (the original headline symptom): ~140k nodes at 85 My for
   `node_density=4` vs a clean-tiling estimate of ~130k *at 1x* -- consistent with the
   ~15-75%-over range the 2026-08-30 table recorded, i.e. not fixed, just no longer
   dominated by the winding rows.

The **k-means split-cluster quality** noted under bug 2 / bug 3 (velocity-space clusters cut
from spatially-intermingled points -> disjoint daughters that drift back over each other,
with euler poles fit far from the daughter body) is the common root of (2), (3) and (4) and
is still entirely open.

6. **Plate count only decreases -- rifting essentially never fires.** Partly fixed
   2026-09-01 (seed 559394024). With the torque engine a large continental plate gets a
   *good* rigid-rotation fit, so `merge_split.SPLIT_RMS_RESIDUAL_THRESHOLD` (9 cm/yr) was
   never tripped; `SPLIT_SIZE_CERTAIN_RIFT_RAD = pi` only relaxed the gates for
   near-hemisphere plates; and the great-circle cut between the two k-means flow centroids
   frequently left one half below `SPLIT_MIN_NODES` (1200) so `plate.split()` rejected it.
   Instrumented over 70 steps: plates repeatedly cleared both physics gates only to fail the
   size floor, and true rift-splits basically never happened while merges + oceanic
   consumption ran unchecked -- 12 plates decaying to 9 in a 180-step reproduction, oceanic
   plates 6 -> 3. Tuning: `SPLIT_RMS_RESIDUAL_THRESHOLD` 9->6 cm/yr, `SPLIT_MIN_POLE_SEPARATION`
   6->4, `SPLIT_MIN_NODES` 1200->700, `SPLIT_SIZE_CERTAIN_RIFT_RAD` pi->2.2 rad,
   `SPLIT_MIN_AGE_STEPS` 15->20; plus `apply_topology_changes` now splits **at most one plate
   per step** (same incremental rule the merge path uses), so a freshly-generated world
   staggers its rifts instead of shattering all at once (1.8 rad was tried first and
   shattered every plate within ~30 My). A 200-step repro now oscillates 15-26 plates in a
   genuine churn (splits + merges + consumption). **Still open:** the k-means cluster-quality
   problem above -- a daughter's euler pole is still fit from spatially-intermingled velocity
   clusters, so the cut geometry can still be poor; and there is still no mechanism to spawn a
   *new* plate for a large region a consumed oceanic plate vacated (the old `gaps.py`
   fallback, never ported). Plate count held healthy in the repro without it, so a spawn net
   was not added.

### Node-count creep: continental boundaries grow but never retreat (2026-09-01 investigation)

**The dominant driver of the node-count blowup (item 5) and of the "plates overlap
neighbours" / over-stretched-continent geometry (items 2-4) is one asymmetry, and it is not
k-means split quality.** `lithosphere_plate.LithospherePlate.deform` sets
`shrinkable_all = np.zeros_like(...)` for a continental self-plate (only oceanic crust
subducts), so a continental line's *contested* (leading) end is a no-op every step while its
*uncontested / divergent* (trailing) end still grows -- each continental row ratchets
outward and never back. Measured (seed 936513024, `node_density=1`, 100-ky steps,
climate/biomes off, instrumented `_grow_or_shrink_line_for_deform` / `_claim_adjacent_
territory` / `regularize_line` / `apply_topology_changes`):

| | 0 My | 80 My | 160 My |
|---|---|---|---|
| continental nodes | 15,976 | 22,120 (+38%) | 33,097 (+107%) |
| oceanic nodes | 16,675 | 12,620 (-24%) | 4,826 (-71%, being consumed) |
| total / clean-tiling ratio | 1.00 | 1.06 | 1.16 |

Cumulative over the 1600-step run: `endgrow_continental +19,364`, `claimrow_continental
+1,321`, **`endshrink_continental` exactly 0** -- continental crust adds ~21k boundary nodes
and removes none. Oceanic is ~balanced (`endgrow +45,549 + claim +2,801` vs `endshrink
-55,610` plus `-6,723` topology). A second seed (42) tracks: continental +43% by 80 My,
`endshrink_continental` 0. The *total* ratio understates the damage because ocean
consumption partly masks it -- the geometric symptom is the unbounded continental growth and
the envelope overlap it drives, exactly items 2/3/4.

**The land-area side of this is fixed 2026-09-01, the node-count side is not.** The same
ratchet meant `_grow_or_shrink_line_for_deform` / `_claim_adjacent_territory` seeded every
new continental boundary node with `reference_thickness("continental")` -- Hc 35 km /
Hm 100 km -> `isostatic_elevation` = **+200 m dry land** -- so a continental plate growing
into vacated ocean permanently converted sea floor into land (measured land fraction climbed
0.27 -> 0.48 and mean planet elevation rose ~1.7 km over 180 Myr on seed 559394024, a
separate reported bug: "the amount of land is increasing over time"). Fix:
`lithosphere_plate.growth_seed_thickness()` -- brand-new areal crust is *always* the oceanic
reference column regardless of plate type ("any gap that opens on the sphere is floored by
sea-floor spreading"). A continental plate's new margin now lands ~-3.5 km (a drowned
passive margin / accreted terrane), so land fraction goes flat and the "giant 80%-drowned
continental plate" of item 2 is now the *expected* reading rather than a bathymetry
artefact. **This does not touch the node count** -- the rows still ratchet outward and
`endshrink_continental` is still 0; directions 1-2 below are still the fix for that. It also
left a smaller residual: land fraction drifts *down* over hundreds of Myr (the continental
node ratchet keeps tiling drowned ~-3.5 km accreted margin onto every continental plate, and
a fixed sea level can't compensate the way a real ocean would).

**Orogenic renewal was completely dead until 2026-09 -- fixed.** The residual above was
originally blamed partly on "erosion planing continents with no orogenic renewal." It turned
out there was *no* orogenic renewal at all, from any collision, ever:
`rheology.EFFECTIVE_LITHOSPHERE_VISCOSITY_PA_S_PER_M` (the closing-rate -> normal-stress
proxy scale) was ~3-4 orders of magnitude too small, so `rheology.plastic_strain_rate_per_
myr` never cleared the Mohr-Coulomb yield stress at *any* closing rate `mantle.MAX_PLATE_
RATE` (15 cm/yr) permits -- `apply_convergent_deformation` returned Hc unchanged on every
contested node. Instrumented on seed 926698457 @ 239.6 My: ~48k continental nodes/step
through the convergent path (~24k actually converging, median 1.6 cm/yr), **0 past yield,
0.0 m Hc added.** So the `overlapAge` view's multi-plate collisions stuck since 126-176 My
never built the mountains this doc and item 4 claim, continents only ever thinned (rift +
erosion-isostasy) and drowned, and land fell monotonically (~-0.003/step with climate on
that seed). Fix: `EFFECTIVE_LITHOSPHERE_VISCOSITY_PA_S_PER_M` 3e13 -> 1e17, so a sustained
~3 cm/yr collision sits a few x past yield (plastic strain ~0.016/Myr, Hc doubles over
~45 Myr -- the Himalaya/Tibet timescale the constant's docstring always claimed). A fresh
seed-926698457 world stepped 300x (node_density 2, climate off): land trajectory unchanged
(-0.019 vs -0.021 over 30 Myr), node-count creep **halved** (+3.6% vs +7.5%), and
continental orogen crust (p95 Hc) now *builds* (+682 m) where it used to slowly erode away
(-212 m). On the collapsed save the decline roughly halves and continental median Hc climbs
instead of falling. `unit_tests/test_rheology.py` pins the calibration. The node-ratchet
residual (directions 1-2) is unchanged -- that is now the dominant land-fraction driver.

**Naive fix rejected.** Forcing continental contested ends to retreat like oceanic
(`shrinkable = contested`) cut continental growth from +38% to +7.7% at 80 My, but (a) the
freed ground is immediately re-claimed by the oceanic neighbour so the *total* barely moves
(1.06 -> 1.04 at 80 My, still climbing to ~1.10 by 110 My), and (b) aggressive edge deletion
severs continental lobes -> defragmentation spawns spurious plates (plate count 18 -> 23+ by
110 My). Not viable as-is.

**Directions worth trying, roughly in effort order:**

1. **Retreat continental edges only against an oceanic neighbour, or only a deep contested
   run.** **DONE 2026-09-02, then extended to continental neighbours 2026-09-02
   (seed 656865324 @ 60 My).** `LithospherePlate.deform` now sets `shrinkable_all =
   _runs_of_at_least(contested_all, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN)` for continental
   crust (was all-False, then briefly `contested_all & inputs.neighbor_is_oceanic`), and the
   interior-subduction carve is gated to `crust_type == "oceanic"` so a continental row can
   never be carved mid-line into a spurious defrag plate. `CONTINENTAL_CONTESTED_RETREAT_MIN_
   RUN = 3` + the pre-existing `n_distance_cap` (= 1 at real continental drift rates) are the
   "gate on run length, cap at 1 node/step" this direction called for; the earlier
   ocean-only gate is dropped because a stalled continent-continent suture needs to consume
   its overlap too (item 4, below). The retreated shortening is channelled into extra plastic
   thickening at the contested nodes (`rheology.CONTINENTAL_COLLISION_SHORTENING_BOOST`, a >1
   `fault_factor` multiplier) so the consumed overlap builds real relief.
   `test_lithosphere_continental_contested_edge_retreats` /
   `test_continent_continent_suture_thickens_faster_than_the_bare_yield_rate` pin it. On
   seed 656865324 the plate-8/6 17%-since-33-My overlap drains to ~2% within ~1.5 My, plate
   count flat at 16, node count flat. **Caveat, as predicted here:** against an *oceanic*
   neighbour the freed ground is re-claimed as new oceanic crust so the *node count* barely
   moves, but it stops the *land-fraction* bleed. Direction 2 is still the complement for the
   raw node count.
2. **Cap a plate's total footprint against its crustal volume.** **DONE 2026-09-02** -- the
   volume-budget growth gate, mechanism 1 of the
   [Continental ratchet: solution design](#continental-ratchet-solution) below (which also
   carries the retreat-mechanism inventory and suture-orientation regime analysis). Halves
   continental-node growth over a 120 My run. The parallel-suture leading-row drop
   (mechanism 2 there) is the remaining piece.
3. **Make frozen continent-continent overlaps actually resolve.** **DONE 2026-09-02.**
   `merge_split.update_overlap_progress` / `World.overlap_progress` is a second sustained-timer
   -- the territory-overlap sibling of `collision_progress` -- and `pop_ready_forced_merge`
   fuses any continental pair that has interpenetrated by >= `FORCED_MERGE_OVERLAP_FRACTION`
   (0.30) for >= `FORCED_MERGE_SUSTAINED_YEARS` (30 My), bypassing the size/closing-rate roll
   entirely. This is what the "or force a merge once an overlap has been stable-and-large for
   N steps" option called for -- and it catches the case `collision_progress` structurally
   can't (item 4: a pair overlapping too completely to register any closing rate). Separately,
   `_merge_probability` now adds a `SPEED_MERGE_BOOST` term: `|omega_a - omega_b| /
   MAX_PLATE_RATE` lifts a large pair's odds toward certainty (a fast head-on collision merges
   more readily than a slow graze). The `MERGE_PROBABILITY_FLOOR` itself is unchanged.
   `test_forced_merge_*` / `test_merge_probability_speed_boost_*` pin both. The deep-overlap ->
   `contested` -> Hc-thickening path (partly-addressed note below) still applies in the
   meantime, so a pair crumples-and-thickens while its forced-merge timer runs.

<a id="continental-ratchet-solution"></a>
### Continental ratchet: solution design (2026-09-02)

Design pass on how to actually stop the ratchet -- item 5, this section, and the land-fraction
sweep further down all reduce to it. Grounded in a full read of `lithosphere_plate.deform`,
`_grow_or_shrink_line_for_deform`, `_claim_adjacent_territory`,
`rheology.apply_convergent_deformation`, `elevation_lines.regularize_line` /
`split_into_contiguous_runs`, `torque.classify_boundary_nodes`.

**What retreat can and can't do today.** A plate is a stack of `ElevationLine` rows at fixed
plate-local `phi`, each a theta-sorted node array plus parallel Hc/Hm arrays; isostasy derives
`elevation` from Hc (you never set elevation directly). Growth/shrink is:
- **end-only per row** -- `_grow_or_shrink_line_for_deform` trims `shrinkable` nodes only
  where the run reaches `theta[0]` / `theta[-1]` (`contested_run_from_end`).
- the interior carve is gated `crust_type == "oceanic"` -- carving a continental row mid-span
  severs the landmass into a spurious defrag plate.
- `_claim_adjacent_territory` only ever *adds* a row past a phi extreme. **Nothing anywhere
  removes a whole leading row.**
- the 2026-09-02 continental retreat (`_runs_of_at_least(contested_all, 3)`) does **not**
  conserve the retreated column's volume -- it drops the mass and multiplies `fault_factor`
  by `CONTINENTAL_COLLISION_SHORTENING_BOOST = 2.5` as a proxy.

**Suture orientation vs. the row grid decides which retreat op is even possible.** Three
regimes, two of them with no implementation:

| Suture vs. *that plate's* rows | Contested nodes land... | Today |
|---|---|---|
| **Transverse** (crosses many phi-rows) | bunched at one theta-end of each row | end-trim works (runs >=3, ends only) |
| **Oblique** (enters a row mid-span) | a run mid-row, live nodes both sides | nothing -- interior carve is oceanic-only; the tongue just thickens in place |
| **Parallel** (compression along phi) | the whole frontmost row(s), full theta width | nothing -- no uncontested end, no whole-row removal -> continental plate *cannot* retreat |

Gotcha: the grid is plate-local, each plate has its own `frame`, so one suture is
simultaneously *parallel* to plate A's rows and *transverse* to plate B's. You cannot punt
the parallel case hoping the neighbour handles retreat -- A's node pile ratchets regardless.

**On the "delete a node from each plate, respawn a thicker one" idea** (mass-conserving suture
consumption). Right direction, more honest than the `fault_factor` fudge, lets
`CONTINENTAL_COLLISION_SHORTENING_BOOST` be deleted. Three refinements:
1. *Not symmetric.* Each plate's `deform()` reads the other's polygon live. If both retreat
   their frontmost node and the rate doesn't track the closing rate, you either never heal
   the overlap or open a gap between two colliding continents that classifies `divergent` ->
   spurious rift. Pick an indentor (larger plate / larger `|omega|` share of the closing
   rate); the *overridden* plate loses its node, the *indentor* keeps its node and absorbs
   the deleted column's Hc (this is Tibet).
2. *Cap retreat at the closing distance* -- the same `n_distance_cap` the oceanic path uses,
   else gap/overlap oscillation.
3. *"Higher elevation" = larger Hc on the surviving node.* Sum the two columns' `Hc * area`;
   `regularize_line` re-evens spacing next pass; isostasy lifts it.
   Cleanest framing: accretion / terrane-transfer (move B's Hc onto A) rather than delete-and-
   respawn. Only really helps the transverse regime (already partly handled) plus continent-
   continent (forced-merge timer already backstops) -- so **lower priority than the volume
   cap.**

**Mechanisms, ranked.**

1. **Volume-budget growth gate (this is "direction 2" above) -- highest leverage, do first.**
   **DONE 2026-09-02 (seed 60461418, 120 My at node_density=1).**
   `lithosphere_plate.deform` computes, for a continental self-plate only:
   `n_continental = count(Hc >= CONTINENTAL_BUDGET_HC_FRACTION * REFERENCE_HC_CONTINENTAL_M)`
   (fraction 0.6) and `suppress_growth = len(own_points) > CONTINENTAL_AREA_BUDGET_MULT *
   n_continental` (mult 1.8). When set, the two `grow_end` blocks in
   `_grow_or_shrink_line_for_deform` are skipped and `_claim_adjacent_territory` is not
   called at all this step -- but retreat (`shrinkable`), divergent thinning and convergent
   thickening all keep running, so an over-budget plate thins / drowns / crumples back toward
   its crustal volume rather than merely freezing. Regime- and neighbour-independent.
   Effect over 120 My: continental-node growth **+7.4% vs +18.5%** with the gate disabled
   (roughly halved -- the same order as the rheology-calibration fix); total node count
   +1.6% vs +2.3% (ocean consumption masks most of it, as predicted); land fraction a hair
   better (0.428 vs 0.417 at 120 My -- the gate deliberately lets divergent thinning
   continue, so it doesn't stop the interior-drowning half). `plates.py` v1 engine
   untouched (not the running one). `test_lithosphere_continental_volume_budget_suppresses_
   growth` (unit) + `test_continental_volume_budget_bounds_the_boundary_ratchet` (stress)
   pin it. **Left for the leading-row drop (2 below):** the *parallel-suture* regime, where a
   plate's frontmost row is entirely contested -- the gate freezes its growth but nothing
   removes the row, so that pile-up still can't retreat.
2. **Make `_claim_adjacent_territory` reversible -- a leading-row drop.** The structural fix
   for the parallel-suture regime: the existing claim logic with the sign flipped, run at the
   same point in `deform()`. If a plate's outermost phi-row is >= ~70% contested for >= N
   sustained steps, delete the whole row. Whole-row removal keeps the plate contiguous (the
   lobe-severing hazard is specific to *mid*-plate carving), so far safer than splitting rows.
3. **Suture consumption as accretion, replacing (not stacking on)
   `CONTINENTAL_COLLISION_SHORTENING_BOOST`.** As above -- physical honesty, lower urgency.
4. **Periodic conservative continental re-lattice.** `build_lines_from_lattice` already
   rebuilds a plate's rows from an outline + ownership predicate. Every K steps, refit the
   lattice to the *current outline* and redistribute the existing total `sum(Hc * area)` onto
   the new node set -- the 2-D generalisation of `regularize_line`. `grow_into` was rejected
   for per-*step* use (its coverage radius balloons the plate); as a periodic re-fit-to-
   outline that objection may not hold. Prototype-worthy.

**Recommendation.** Volume cap (1) first as the regime-free runaway-killer; leading-row drop
(2) for the parallel-suture gap; suture accretion (3) later for honesty. The design rule the
regime table implies: for continental crust, prefer whole-row ops + volume caps over mid-row
carving/splitting.

**Progress:** mechanism (1) landed 2026-09-02 (see its entry above). (2)-(4) still open.

**Fixed here (2026-09-01): the v1 pole-winding guards were never ported to the v2 engine.**
"Bug 1" (below) added a `ring_room()` one-revolution cap in
`_grow_or_shrink_line_for_deform` and a `POLE_CAP_MARGIN_MULT` clearance in
`_claim_adjacent_territory` -- but only to `plates.PlateWithLines`. The running
`lithosphere_plate.LithospherePlate` overrides of both still used `max_phi_limit = pi/2 -
spacing/2` (marches onto the pole) and had no revolution cap at all, leaning entirely on
`regularize_line`'s after-the-fact unwind (rows over-wound past 2*pi and were unwound again
every step -- churn, wasted RNG draws, near-pole 1-3 node rings). Both guards are now ported
(`ring_room()` + `POLE_CAP_MARGIN_MULT`), with `test_lithosphere_deform_never_winds_a_row_
past_a_full_revolution` / `test_lithosphere_claim_adjacent_territory_keeps_a_margin_from_
the_local_pole` mirroring the v1 tests. Effect on the headline node count is small at the
seeds checked (pole winding was a minor contributor next to the continental ratchet -- max
row span 185 deg -> 147 deg, near-pole nodes ~halved), but it removes the per-step
regularize churn and brings the two engines back to parity.

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

## Land fraction slowly declines over a long run

**The erosion half is FIXED (2026-09-01, erosional isostatic compensation):** `apply_erosion`
books its whole per-step geomorphic change against `crustal_thickness_m` and moves `elevation`
by exactly the resulting `isostatic_elevation` delta, so an unloaded crustal root rebounds
(only ~1/6 of subaerial erosion, ~1/4 of submarine, survives as a surface drop) -- the same
delta idiom `deform()` uses for tectonic Hc/Hm changes. On seed 331015891 the node-level
land-fraction decline dropped ~7x. See `lithosphere.isostatic_elevation` (per-node rho_c),
`test_erosion.py::test_apply_erosion_thins_crust_where_it_erodes_and_isostasy_compensates`,
`test_world_smoke.py::test_erosion_books_against_hc_and_survives_deform`. The `elevReason`
view's `moved`/override gates key off the *raw* geomorphic move, since compensation shrinks
the surface expression ~5x without changing which process is shaping the column.

**The tectonic half is NOT fixed** -- see the sweep below. Land keeps falling on a long run
because plate movement drowns continental crust faster than anything lifts it, and that is
now the dominant driver (erosion contributes < 15%).

### Toggle sweep: land vs node count, seed 60461418 @ 69 My (2026-09-02, corrected)

**Reported symptom.** On `~/Downloads/mantle-bloom-seed60461418-69000000y.mbworld` (seed
60461418, 69 steps of 1 My, `node_density=4`, 19 plates, 16 continental and 12 of those
already >50% submerged): with **Climate turned OFF**, elevation-point count keeps rising and
the Stats-panel land fraction sits still. Both observations are **correct** -- see below.

**Two land-fraction numbers, and they disagree.** `node LF` = fraction of *nodes* whose live
`elevation` is above sea level (recomputed directly every sample). `grid LF` = what
`stats.compute_stats` / the Stats panel report -- the climate grid's `is_ocean` fraction,
which comes from `hydrology.sample_is_ocean` resampling `world.hydrology_cache.is_ocean`.
**`hydrology_cache` is only rebuilt inside `erosion.apply_erosion`, so with
`simulate_climate_biomes` off it is frozen at its step-69 state and `grid LF` cannot move at
all** -- 0.202 forever, no matter how far elevations actually sink. So the panel showing
"land static" with Climate off is a stale-cache artefact, not a real measurement; the live
`node LF` under it is falling the whole time.

**Long-run, everything ON vs Climate OFF (both stepped 69 -> 219 = +150 My).**

| | nodes | continental | oceanic | plates | node LF (live) | grid LF (panel) | node mean elev |
|---|---|---|---|---|---|---|---|
| **all on**, 69 My  | 131,979 | 111,268 | 20,711 | 19 | 0.244 | 0.202 | -2892 m |
| **all on**, 219 My | 132,693  (+0.5%) | 117,866  (**+5.9%**) | 14,827  (**-28%**) | 18 | **0.147** | **0.104** | -3195 m |
| **Climate OFF**, 69 My  | 131,979 | 111,268 | 20,711 | 19 | 0.244 | 0.202 | -2892 m |
| **Climate OFF**, 219 My | 134,614  (**+2.0%**) | 116,758  (**+4.9%**) | 17,856  (**-14%**) | 16 | **0.144** | **0.202** (frozen) | -3271 m |

Climate-OFF node count: 131,979 -> 132,485 (step 99) -> 133,680 (149) -> 134,614 (219),
climbing ~+2% per 100 My in the back half and still accelerating. **So the reported symptom
is real -- the elevation-point count does keep going up.** (An earlier version of this note
called it "flat", from a 40-step run that stopped at step 109 -- +0.2%, before the ramp;
retracted.)

**What's actually happening.**

1. **Continental boundaries ratchet outward ~+5-6% per 150 My, in every config**
   (`lithosphere_plate.deform` sets `shrinkable_all` to almost nothing for a continental
   self-plate -- see the "Node-count creep" item above). This is the robust, config-independent
   signal and the root cause of the elevation-point growth.

2. **Whether that shows up in the *total* depends on how fast oceanic crust is consumed to
   compensate**, and that is trajectory-sensitive (slab pull in `torque.py` reads bathymetry,
   so once elevations diverge the two runs consume different oceanic plates at different
   rates): -28% oceanic with everything on nearly cancels the continental growth (total
   +0.5%); -14% with Climate off does not (total +2.0% and rising). Either way the continental
   node pile only grows.

3. **Land loss is ~85-90% tectonic.** node LF falls 0.244 -> 0.147 (all on) vs 0.244 -> 0.144
   (Climate off) over the same 150 My -- erosion barely changes it. With erosion off the total
   node count still grows yet node LF still drops, so *existing* above-sea continental nodes
   are being pushed underwater -- `deform()` subsidence (divergent thinning of shear-stretched
   over-large plates + Airy isostasy), plus every newly grown margin node seeded at the
   drowned oceanic reference column (~-3.5 km, `lithosphere_plate.growth_seed_thickness`). Not
   instrumented which of the two dominates.

4. **Volcanism and the wind model are irrelevant** to both metrics (40-step sweep: baseline
   vs no-volcanism vs `wind_model="cfd"` differ < 0.002 in land fraction, < 15 m in mean
   elevation). **Plate movement is the whole story**: with `simulate_plate_movement` off, node
   count is exactly flat and node LF barely moves (0.244 -> 0.242 in 40 My).

**Direction.**

- **The stale-`hydrology_cache` `grid LF` freeze is its own bug** -- the Stats panel and
  `/world/stats` `land_fraction` silently stop tracking reality whenever Climate is off. Either
  refresh a cheap connectivity mask each step regardless of `simulate_climate_biomes`, or mark
  the stat stale in the response / panel. (`node LF` -- a bare `elevation > sea_level` count --
  is always right and would be a good panel addition on its own.)
- **The node-count + land-loss driver is the continental boundary ratchet** -- see
  [Continental ratchet: solution design](#continental-ratchet-solution). The volume-cap
  (direction 2 / mechanism 1 there) **landed 2026-09-02** and halves the continental-node
  growth over a long run (+7.4% vs +18.5% at 120 My); it does not fix the land-fraction
  decline, because it deliberately lets divergent thinning keep drowning an over-budget
  plate's interior. Directions 1 (contested-run retreat) and 3 (forced merge) were partly
  landed earlier and did not stop either trend on their own.

---

## Speckled low-relief coastlines: a drowned flat shelf dithers pixel-by-pixel across sea level

**Status: PARTIALLY RESOLVED (round 2).** Round 1 -- the coastal planation + infill feedback
(option 2) -- cut the density-1 dither roughly in half; barrier islands (option 3) fall out
of it in emergent form. A fresh density-4 look
(`~/Downloads/mantle-bloom-seed888151728-85000000y.mbworld`, 2026-08-31, seed 888151728,
851 steps / 85.1 My) then showed the coast **still visibly speckled** along the south/east
coasts of the plate-3 land mass around 5 N / -10 E -- the feedback was running but *swamped*
by lumpy river deposition (diagnosis below). Round 2 (2026-08-31, "What landed (round 2)")
de-clumps that deposition with a distributary spread and merges planation + infill into one
symmetric leveling pass across a band straddling sea level. Still needs a fresh density-4
measurement pass to confirm the checkerboard is gone (see Verification in the round-2 notes).

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

**What landed (round 2), 2026-08-31.** All in `erosion.py`, all mass-conserving via
`np.add.at`, stateless per-step (same idiom as the other `_spread_*` helpers), no new
persistent fields:

- **Distributary redirect.** In `apply_erosion`, after the beach spread: a near-sea-level
  land node (`|elev - sea_level| <= COASTAL_LEVELING_BAND_M`, 45 m) flagged `is_depositing`
  (the existing slow-big-river `DEPOSITION_*` test) has `DELTA_REDIRECT_FRACTION` (0.8) of the
  lump `route_downstream` piled on it pulled back out and fed to the leveling fill spread,
  which scatters it across the band's below-datum hollows -- the emergent distributary fan.
  The rest stays put as the active channel bar / natural levee.
- **`coastal_planation_amount` + `_spread_coastal_infill` -> `coastal_leveling_grind` +
  `_spread_coastal_leveling`**, one symmetric pass over a band straddling sea level. Every
  band node has a local target datum (`leveling_datum_m`: a single continuous function of
  wave exposure -- `sea_level - LEVELING_PLATFORM_UNDERCUT_M * exposure + LEVELING_MARSH_CREST_M
  * shelter`, or `sea_level + BARRIER_CREST_M` for a barrier candidate). The grind half planes
  down every node standing above its datum (a just-submerged shoal included -- the old
  planation gate ignored everything below sea level); the fill half silts up every node below
  its datum, land or ocean alike (a dry interdistributary low fills as readily as sheltered
  shallow water). `prominence` / `hollow` still reweight both at node scale.
- The fill side is a **capped iterative water-fill**: each below-datum sink has a hard
  per-step capacity (its metres of room to the datum, bounded by one `LEVELING_RATE` step),
  and `LEVELING_FILL_ITERS` passes distribute each source's remaining load across its
  still-open sinks (`LEVELING_SPREAD_NEIGHBOR_COUNT`), capping and carrying the overflow. That
  is what makes a 200 m lump genuinely spread across the flat plain instead of piling onto the
  single most-weighted node -- the failure mode round-1's `_spread_coastal_infill` still had.
- The fill pool is the ground-off rock + `COASTAL_INFILL_MARINE_FRACTION` (0.5) of the
  submarine/coastal erosion pool + the distributary redirect; the rest of the sea-side pool
  still spreads to deep water via `_spread_marine_sediment`. `LEVELING_LOCAL_FRACTION` (0.25)
  of a spread stays put only when the source is itself a below-datum node; a pure source
  spreads in full, and never back onto itself; whatever no sink had room for bounces back.

Renamed constants: `PLANATION_RATE_M_PER_MYR -> LEVELING_RATE_M_PER_MYR` (**250 -> 60**;
round-1's 250 over-planed the drowned shelf), `PLANATION_EXPOSURE_REF -> LEVELING_EXPOSURE_REF`,
`PLANATION_UNDERCUT_M -> LEVELING_PLATFORM_UNDERCUT_M`, `INFILL_MARSH_CREST_M ->
LEVELING_MARSH_CREST_M`; `PLANATION_BAND_M` + `INFILL_DEPTH_M` collapsed into
`COASTAL_LEVELING_BAND_M` (45 m); `COASTAL_INFILL_MARINE_FRACTION` kept 0.5. New:
`DELTA_REDIRECT_FRACTION`, `LEVELING_LOCAL_FRACTION`, `LEVELING_MIN_OPENNESS`,
`LEVELING_SPREAD_NEIGHBOR_COUNT`, `LEVELING_FILL_ITERS`.

Tests: `unit_tests/test_erosion.py` (`leveling_datum_m`, `coastal_leveling_grind`,
`_spread_coastal_leveling` -- direction, exact mass conservation incl. the no-sink fallback,
and a single-lump-declumps-across-neighbours check) and the renamed stability floor
`stress_tests/test_world_stepping.py::test_coastal_feedback_stays_stable_over_many_steps`.

**Measured so far (apply_erosion-only loop, 10 steps, seeds 830054688 / 505070493 / 443034896
-- the seed-888151728 density-4 save was not on disk and the older saves can't `step_world`,
so this is a proxy; feedback ON vs OFF).** In the near-sea-level band:

| | flip fraction | band \|dElev\|/step p90 | deposition top-10% share |
|---|---|---|---|
| seed 830054688 | 0.32 -> 0.25 | 37 -> 21 m | 0.35 -> 0.23 |
| seed 505070493 | 0.36 -> 0.30 | 41 -> 25 m | 0.60 -> 0.39 |
| seed 443034896 | 0.32 -> 0.23 | 40 -> 24 m | 0.72 -> 0.55 |

Every metric improves on every save -- the deposition-concentration drop is the de-clumping
working directly. Round-1's pass on the *same* proxy went the wrong way (flip 0.40 / 0.50,
p90 ~73 / ~90 m) -- its 250 m/My planation rate and one-shot infill scatter both amplify
node-scale noise here; round 2 is a clear improvement over both round 1 and no feedback.

**Still to verify (round 2).** No fresh density-4 `step_world` measurement yet (needs a
loadable density-4 save or a long from-scratch run). Re-run the 2026-08-31 investigation
against a seed-888151728-style world (~25 `step_world` steps of 100 ky), in the ~8x30 deg
box over 5 N / -10 E: band nodes flipping land<->ocean per step (was 315/768), deposition
top-10% share (was 51%), per-step |dElev| p10/p90 in the band (was +-50 m). The
transient-lake log spam and the stranded deep basins (below) are downstream symptoms of the
same dithering shelf and should ease.

**Related symptom -- event-log flooding.** ~~The world's event log is almost entirely
"N-node lake formed/split ... at elevation ~0 m" -- hundreds of these per My.~~ **Fixed
2026-08-31** by `lakes.summarize_lake_events` (see the "Diagnostic views & debug output"
section below): a step's near-sea-level transients collapse to one aggregate line, so real
basin events are legible again (e.g. the genuine persistent ~435-node lake oscillating
around -1770 m, which still logs individually). The dither this was a symptom of is still
open (see "Fix direction" above).

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

**What landed (option 2 + emergent option 3), 2026-08-31 -- round 1, superseded by round 2
above.** `coastal_planation_amount` / `_spread_coastal_infill` / the `PLANATION_*` /
`INFILL_DEPTH_M` / `INFILL_MARSH_CREST_M` constants named here were renamed and merged into
the symmetric leveling pass; `_coastal_openness`, the `BARRIER_*` / `PROMINENCE_*` fields,
and the emergent-barrier / prominence mechanics carried over unchanged. A per-step pass in
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

**Options considered (all resolved).** Render-only despeckle: shipped then deleted as
near-redundant. Coastal planation + infill feedback (the real fix): shipped 2026-08-31.
Barrier islands: delivered in emergent form. Optional future refinement: a genuine
wave-exposure / fetch field plus an explicit longshore-transport direction for real spits
and drift-aligned bars.

---

## Diagnostic views & debug output (from the 2026-08-31 seed-888151728 investigation)

Working through the plate-geometry and coastal-speckle degradation above needed several
numbers the program didn't surface. **Landed** (all 2026-08-31 unless noted; see
`docs/debugging.md`):

- **Plate Inspector motion / shape / overlap fields** -- `GET /world/plates` +
  `main._plate_summary` + the Plate Inspector panel report per plate `speed_cm_per_yr` /
  `at_max_rate` / `euler_pole` / `age_steps` / `median_elevation_m` / `submerged_fraction` /
  `overlaps` / `collisions`.
- **`GET /world/stranded_basins` + `python -m app.stranded_basins`** -- every endorheic basin
  whose floor is below sea level, deepest-first, with `world.stranded_basin_tracks`
  persistence.
- **`python -m app.plate_diagnostics <save.mbworld>`** -- offline per-plate motion/shape
  table, overlap list, collision timers, node count vs clean-tiling estimate. `--json` too.
- **Lake-churn event aggregation** -- `lakes.step_lakes` returns structured `LakeEvent`s;
  `lakes.summarize_lake_events` collapses a step's near-sea-level transients to one line.
- **`GET /world/render?view=geomorph`** (Map View: Debug > Erosion & Deposition) -- colours
  each node by `erosion.ErosionResult.net_elevation_change_m` on a diverging scale.
- **`ElevationLine.overlap_onset_years` + `GET /world/render?view=overlapAge`** (2026-09-01)
  -- when each node first started overlapping another plate; `since_years` on the plates
  endpoint / inspector / `plate_diagnostics`.

**Still worth building:**

1. **A speckle / coastal-dither overlay render mode.** Colour every node whose elevation is
   within a threshold of sea level by the fraction of its k nearest neighbours on the
   *opposite* side of the waterline (the exact metric the investigation scripts compute:
   `near = |elev - sea_level| < 120 m`; `frac = mean(neighbour is opposite class)`; flag
   `>= 0.75`). Instantly shows where the coast is a checkerboard vs a clean shoreline, and
   makes a before/after for any coastal-feedback change legible without an ad-hoc script.
   Natural home: a `render_image.py` view alongside `platesDetail`, or a boolean overlay on
   the `elevation` / `biome` views.
2. **A map-view render of the stranded-basin report** -- the endpoint already returns
   `centroid_xyz` / `floor_xyz` ready for one.
3. **A `/world/sample_at` field for the geomorph view's click-popup** -- the popup is only
   wired for elevation / biome / combined today.

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

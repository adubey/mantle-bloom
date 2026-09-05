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
   The retreated column's volume is conserved as accretion onto the plate's own leading edge
   (`_redistribute_accreted_column`, mechanism 3 below -- **done 2026-09-02**), so the overlap
   crumples into an orogen.

   **Overlap severity now also scales directly, 2026-09-04.** Beyond the mechanisms above
   (which resolve the overlap over time), a *currently* overlapping pair now feels a
   proportionally stronger force while it lasts: the near-field contested-band uplift rate and
   `torque.collision_friction_torque`'s resistive braking are both scaled up by how much of
   the plate's own near-boundary band is contested right now (`OVERLAP_UPLIFT_SEVERITY_GAIN` /
   `OVERLAP_FRICTION_SEVERITY_GAIN`), and the same severity feeds `Plate.internal_stress`
   (direction 7 further down) so a large plate stuck in a deep overlap is pushed toward
   rifting, not just toward merging or crumpling in place. See docs/simulation-model.md's
   "Accumulated breakup stress" subsection.

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
   clusters, so the cut geometry can still be poor.

   **The "no mechanism to spawn a new plate for a large vacated region" half is FIXED
   2026-09-04.** Confirmed on the 399 My save below that this does matter on long enough
   runs even though the <=200-step repro stayed healthy without it: a new `gaps.py` (the old
   fallback of the same name, ported to the current `LithospherePlate` engine) sweeps the
   whole-sphere lattice on the same cadence as `merge_split.defragment_plates`
   (`world.step_world`, every `gaps.GAP_FILL_INTERVAL_STEPS` steps), finds any connected
   region at least `gaps.MIN_GAP_NODES` (scaled by `node_density`, same reasoning as
   `SPLIT_MIN_NODES`) large that no plate's lines currently reach, and spawns a new oceanic
   plate over it via `lithosphere_plate.new_plate` (extended with an `is_owned` predicate so
   it claims only that region, not the whole sphere). Deliberately *only* spawns -- it does
   not also try to absorb a gap into a single dominant bordering plate the way the pre-
   refactor module did, to avoid handing continental plates an easy way to feed the
   continental-growth ratchet (directions 1-2 above); ordinary per-step boundary growth
   already does gradual absorption at a plate's own edge once the new plate gives it a real
   neighbour again. `unit_tests/test_gaps.py` pins it (no-op on a freshly-generated world,
   spawns near-full replacement crust for a removed plate's vacated footprint, ignores a
   gap below `MIN_GAP_NODES`). **2026-09-04 addendum:** the spawned plate is no longer
   unconditionally oceanic -- see docs/TODO.md's "`gaps.py`'s plate-spawn is a stopgap"
   section's own 2026-09-04 addendum.

7. **Large plates had no independent pressure to break apart beyond their own instantaneous
   size.** **Addressed 2026-09-04.** `SPLIT_SIZE_CERTAIN_RIFT_RAD` (direction 6 above) only
   ever reads a plate's *current* angular radius -- a plate that has stayed large for a long
   stretch, or that keeps getting shoved by a sustained territorial overlap, got no extra
   credit toward rifting beyond whatever that snapshot already gave it. `Plate.internal_stress`
   (`merge_split.accumulate_plate_stress`) is a real, persisted accumulator: it rises with a
   plate's own radius every step, plus an overlap-driven top-up weighted toward whichever plate
   in an overlapping pair is the *larger* one (`relative_largeness`), and decays without
   renewed forcing. `maybe_split_plate` relaxes the same residual/pole-separation gates by
   whichever of size-based or stress-based relaxation is further along. See
   docs/simulation-model.md's "Accumulated breakup stress" subsection for the full mechanism
   and constants; `unit_tests/test_merge_split.py`'s
   `test_accumulated_stress_alone_can_push_a_weak_split_over_the_gate` isolates it from the
   pre-existing size-based path.

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
   its overlap too (item 4, below). The retreated column's crustal volume is conserved as
   accretion onto the plate's own surviving leading edge -- the attached mantle lithosphere
   thickening in proportion -- via `_redistribute_accreted_column` /
   `SUTURE_ACCRETION_SPREAD_NODES` (mechanism 3 below, **done 2026-09-02**, replacing the
   earlier `CONTINENTAL_COLLISION_SHORTENING_BOOST` fudge), so the consumed overlap builds
   real relief in proportion to what it ate.
   `test_lithosphere_continental_contested_edge_retreats` /
   `test_continent_continent_suture_consumes_its_overlap_as_mass_conserving_accretion` /
   `test_redistribute_accreted_column_conserves_crustal_volume` pin it. On
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

**What retreat can and can't do today.** *(The "nothing removes a whole leading row" gap
below was closed 2026-09-02 by ranked-mechanism 2 -- `_retreat_contested_leading_rows`. The
rest of this design pass stands.)* A plate is a stack of `ElevationLine` rows at fixed
plate-local `phi`, each a theta-sorted node array plus parallel Hc/Hm arrays; isostasy derives
`elevation` from Hc (you never set elevation directly). Growth/shrink is:
- **end-only per row** -- `_grow_or_shrink_line_for_deform` trims `shrinkable` nodes only
  where the run reaches `theta[0]` / `theta[-1]` (`contested_run_from_end`).
- the interior carve is gated `crust_type == "oceanic"` -- carving a continental row mid-span
  severs the landmass into a spurious defrag plate.
- `_claim_adjacent_territory` only ever *adds* a row past a phi extreme. **Nothing anywhere
  removes a whole leading row.**
- the 2026-09-02 continental retreat (`_runs_of_at_least(contested_all, 3)`) **does** conserve
  the retreated column's volume as of mechanism 3 below (**done 2026-09-02**):
  `_redistribute_accreted_column` thrusts the dropped nodes' summed Hc/Hm onto the surviving
  leading edge, replacing the old `CONTINENTAL_COLLISION_SHORTENING_BOOST = 2.5` `fault_factor`
  fudge. A retreat against an *oceanic* neighbour still drops the column (real subduction).

**Suture orientation vs. the row grid decides which retreat op is even possible.** Three
regimes, two of them with no implementation:

| Suture vs. *that plate's* rows | Contested nodes land... | Today |
|---|---|---|
| **Transverse** (crosses many phi-rows) | bunched at one theta-end of each row | end-trim works (runs >=3, ends only) |
| **Oblique** (enters a row mid-span) | a run mid-row, live nodes both sides | nothing -- interior carve is oceanic-only; the tongue just thickens in place |
| **Parallel** (compression along phi) | the whole frontmost row(s), full theta width | whole-row drop after 5 My sustained >=70% override (`_retreat_contested_leading_rows`, done 2026-09-02) |

Gotcha: the grid is plate-local, each plate has its own `frame`, so one suture is
simultaneously *parallel* to plate A's rows and *transverse* to plate B's. You cannot punt
the parallel case hoping the neighbour handles retreat -- A's node pile ratchets regardless.

**On the "delete a node from each plate, respawn a thicker one" idea** (mass-conserving suture
consumption). Right direction, more honest than the `fault_factor` fudge, lets
`CONTINENTAL_COLLISION_SHORTENING_BOOST` be deleted. **Landed 2026-09-02 as the within-plate
variant** (`_redistribute_accreted_column`); the refinements below are how it was scoped:
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
2. **Make `_claim_adjacent_territory` reversible -- a leading-row drop. DONE 2026-09-02**
   (branch `feat/leading-row-drop`). The structural fix for the parallel-suture regime: the
   existing claim logic with the sign flipped, run at the same point in `deform()`.
   `LithospherePlate._retreat_contested_leading_rows` (continental crust only) drops a plate's
   outermost phi-row at either extreme once `LEADING_ROW_CONTESTED_FRACTION` (0.7) of its
   nodes have been contested for a cumulative `LEADING_ROW_RETREAT_SUSTAINED_YEARS` (5 My) of
   deform time -- a per-plate `_leading_row_retreat_years` tally keyed by extreme (survives a
   rotation since rows are plate-local; resets on merge/split/load, which only delays a drop).
   Whole-row removal keeps the plate contiguous -- the lobe-severing hazard is specific to
   *mid*-row carving -- so it is gated only by `LEADING_ROW_DROP_MIN_ROWS` (4). It does not
   plumb the dropped column's volume anywhere (parity with the 2026-09-02 end-retreat); the
   newly-exposed row is contested next step and thickens through the ordinary
   `CONTINENTAL_COLLISION_SHORTENING_BOOST` path. The drop/claim asymmetry is deliberate:
   `_claim_adjacent_territory` adds a row the instant space opens, retreat waits out 5 My of
   sustained override so a transient boundary wobble can't thrash a stable margin.
   `test_lithosphere_contested_leading_row_is_dropped_after_sustained_override` pins it.

   **Measured (seed 936513024, node_density 1, climate off, 800 x 100-ky steps -- the same
   config as the "Node-count creep" table above).** The drop *fires* readily: 82 whole-row
   drops over 80 My across ~9 continental plates (several plates dropping a row every few Myr
   as their leading edge is continuously overridden). But aggregate continental node count is
   **unchanged** vs the drop disabled: +12.8% vs +12.9% at 80 My (total nodes +2.4% vs +2.2%,
   plate count 20 vs 19). Two reasons, both expected: (a) on a young from-scratch world almost
   every overridden continental margin faces *ocean*, and the freed ground re-tiles as new
   oceanic crust within a step or two; (b) the drop removes one outermost row per extreme per
   5 My while the plate's *trailing* (divergent) edge keeps growing unbounded every step --
   only the volume cap (1) caps that. So this is a **correctness fix for the stalled
   parallel continent-continent suture** (the `overlapAge` view's frozen multi-plate
   collisions -- verified directly: a plate whose whole front row a neighbour overruns loses
   that row after 5 My instead of never), **not** a node-count-creep fix.

   **Still open:** the volume cap (1) is still the regime-free runaway-killer and remains the
   priority -- it is what stops the trailing-edge growth this mechanism can't touch.
3. **Suture consumption as accretion, replacing (not stacking on)
   `CONTINENTAL_COLLISION_SHORTENING_BOOST`.** **DONE 2026-09-02.** When a continental
   contested end retreats against a *continental* neighbour (`accrete_all = shrinkable_all &
   ~neighbor_is_oceanic` in `deform`), `_grow_or_shrink_line_for_deform` now captures the
   dropped nodes' summed Hc and `_redistribute_accreted_column` thrusts it back onto the
   `SUTURE_ACCRETION_SPREAD_NODES` (3) surviving leading-edge nodes, with the matching
   isostatic-elevation bump -- an imbricate thrust wedge, mass-conserving because node area is
   constant, up to a `SUTURE_ACCRETION_MAX_HC_M` (~2.4x reference Hc) ceiling past which the
   root delaminates (a never-healing suture would otherwise pile every consumed column onto
   the same retreating-edge nodes forever -- measured Hc ran to ~190 km and climbing over
   30 My without it; capped it plateaus at ~87 km, p95/p99/median unchanged from baseline).
   `CONTINENTAL_COLLISION_SHORTENING_BOOST` (the flat 2.5x `fault_factor` proxy) is
   deleted; the convergent path at those nodes is now just the ordinary yield-limited
   thickening. A retreat against an *oceanic* neighbour is untouched -- that column genuinely
   subducts. `test_redistribute_accreted_column_conserves_crustal_volume` (exact mass check) +
   `test_continent_continent_suture_consumes_its_overlap_as_mass_conserving_accretion` (edge
   thickens kilometres against a continental neighbour, ~nothing against an oceanic one) pin
   it. **Not done -- deferred:** the cross-plate *indentor* asymmetry of refinement 1 (larger
   plate keeps its node and absorbs the smaller's column -- true terrane transfer). It needs
   `deform()` to write a *neighbour's* lines, which the engine never does; the within-plate
   version conserves each plate's own consumed crust onto its own belt (a two-sided orogen,
   geologically fine) and the forced-merge timer still backstops a genuinely stuck pair.
4. **Periodic conservative continental re-lattice.** `build_lines_from_lattice` already
   rebuilds a plate's rows from an outline + ownership predicate. Every K steps, refit the
   lattice to the *current outline* and redistribute the existing total `sum(Hc * area)` onto
   the new node set -- the 2-D generalisation of `regularize_line`. `grow_into` was rejected
   for per-*step* use (its coverage radius balloons the plate); as a periodic re-fit-to-
   outline that objection may not hold. Prototype-worthy.

**Recommendation.** Volume cap (1) first as the regime-free runaway-killer; ~~leading-row drop
(2) for the parallel-suture gap~~ (done 2026-09-02); suture accretion (3) ~~later for honesty~~ **done
2026-09-02** (within-plate imbricate version; cross-plate terrane transfer deferred). The
design rule the regime table implies: for continental crust, prefer whole-row ops + volume
caps over mid-row carving/splitting.

**Progress:** mechanisms (1) (2) and (3) landed 2026-09-02 (see their entries above).
(4) still open.

**Fixed here (2026-09-01): the v1 pole-winding guards were never ported to the v2 engine.**
"Bug 1" (below) added a `ring_room()` one-revolution cap in
`_grow_or_shrink_line_for_deform` and a `POLE_CAP_MARGIN_MULT` clearance in
`_claim_adjacent_territory` -- but only to the v1 `plates.PlateWithLines`. The running
`lithosphere_plate.LithospherePlate` overrides of both still used `max_phi_limit = pi/2 -
spacing/2` (marches onto the pole) and had no revolution cap at all, leaning entirely on
`regularize_line`'s after-the-fact unwind (rows over-wound past 2*pi and were unwound again
every step -- churn, wasted RNG draws, near-pole 1-3 node rings). Both guards are now ported
(`ring_room()` + `POLE_CAP_MARGIN_MULT`), with `test_lithosphere_deform_never_winds_a_row_
past_a_full_revolution` / `test_lithosphere_claim_adjacent_territory_keeps_a_margin_from_
the_local_pole` mirroring the (since-removed) v1 tests. Effect on the headline node count is small at the
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
longest gap-free arc; wired into `LithospherePlate.split` /
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
   - `LithospherePlate._claim_adjacent_territory` (`lithosphere_plate.py`) adds new phi rows
     outward all the way to `max_phi_limit = pi/2 - spacing/2` -- it will march a plate
     right onto its own local pole whenever the space is open.
   - `_grow_or_shrink_line_for_deform` (`lithosphere_plate.py`) then extends those near-pole
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
because plate movement drowns continental crust faster than anything lifts it. **Erosion's
net effect on land fraction is within noise of zero** -- if anything marginally
land-*preserving* (deposition + coastal leveling). The earlier "~10% erosional" / "erosion
contributes < 15%" wording overstated it ~5x and is retracted.

### Toggle sweep: land vs node count, seed 60461418 @ 69 My (2026-09-02, re-run)

**Reported symptom.** On `~/Downloads/mantle-bloom-seed60461418-69000000y.mbworld` (seed
60461418, 69 steps of 1 My, `node_density=4`, 19 plates, 16 continental and 12 of those
already >50% submerged): the elevation-point count keeps rising and land keeps falling.
Confirmed.

**Metric.** `land fraction` throughout = **(count of nodes with live `elevation` >
`sea_level_m`) / (total node count)** -- computed straight off the plates every sample, not
the Stats panel's number. (The panel's `land_fraction` reads `world.hydrology_cache.is_ocean`
via `hydrology.sample_is_ocean`, and `hydrology_cache` is only rebuilt inside
`erosion.apply_erosion` -- so with Climate/erosion off it freezes at the toggle-step value
and stops tracking reality. That stale-cache freeze is its own bug, see Direction below.)

**Clean isolation -- each subsystem toggled alone, stepped 69 -> 219 My (+150 My).**

| config | land frac 69 My -> 219 My | delta | above-sea nodes 69 -> 219 | total nodes 69 -> 219 |
|---|---|---|---|---|
| control: plate movement OFF **and** erosion OFF | 0.2442 -> 0.2442 | **0.000** | 32,230 -> 32,230 | 131,979 (flat) |
| **erosion only** (plate movement OFF, erosion/climate ON) | 0.2442 -> 0.2424 | **-0.002** | 32,230 -> 31,989 (-241) | 131,979 (flat) |
| **tectonics only** (erosion OFF, plate movement ON) | 0.2442 -> 0.1445 | **-0.100** | 32,230 -> 19,451 (**-12,779**) | 131,979 -> 134,614 (**+2.0%**) |
| volcanism OFF (else all on) | 0.2442 -> 0.1470 | -0.097 | 32,230 -> 19,628 | -> 133,559 (+1.2%) |
| **baseline** (everything on) | 0.2442 -> 0.147 | **-0.097** | ~32,230 -> ~19,500 | -> 132,693 (+0.5%) |
| Climate OFF (`simulate_climate_biomes=False`) | 0.2442 -> 0.144 | -0.100 | -> ~19,000 | -> 134,614 (+2.0%) |
| `wind_model="cfd"` (else all on) | 0.2442 -> 0.1517 | -0.093 | 32,230 -> 20,044 | -> 132,155 (+0.1%) |

The `"cfd"` row lands within 0.005 land fraction / 40 m mean elevation of baseline across the
full 150 My -- same magnitude as the baseline-vs-volcanism-off spread, i.e. trajectory noise,
not an effect. (`erosion OFF` and `Climate OFF` also land within 0.0005 / a few hundred nodes
of each other at every step -- two independent ways of disabling erosion agree, since the
fluid solve does nothing to `elevation` without erosion.)

**What this says.**

1. **Land loss is ~98% tectonic.** "Tectonics only" reproduces the full baseline decline
   (-0.100 vs -0.097 -- the small overshoot means erosion is *slowing* the loss slightly, not
   driving it). "Erosion only" is -0.002 over 150 My and non-monotonic (dips to -0.0022,
   recovers to -0.0018) -- erosion is in rough equilibrium with its own deposition and the
   coastal-leveling pass. The 2026-09-01 isostatic-compensation fix is holding.

2. **The drowning is `deform()` pushing existing above-sea continental crust under water.**
   In the tectonics-only run 12,779 nodes that were above sea level at 69 My are below it by
   219 My, while the node cloud only grew by 2,635 -- so this is genuine subsidence of
   standing crust, not just dilution by new deep nodes. Mechanism: over-stretched continental
   plates (euler poles far from the plate body -> near-pure shear) thinning their interiors
   via `rheology.apply_divergent_deformation` + Airy isostasy, plus every newly grown margin
   node seeded at the drowned oceanic reference column (~-3.5 km,
   `lithosphere_plate.growth_seed_thickness`). Not instrumented which of the two dominates.

3. **Node-count growth is the continental boundary ratchet.** Continental nodes climb +5-6%
   per 150 My in every plate-movement-on config (`lithosphere_plate.deform` leaves
   `shrinkable_all` ~empty for a continental self-plate -- see the "Node-count creep" item
   above). Whether it shows in the *total* depends on how fast oceanic crust is subducted to
   compensate, which is trajectory-sensitive (slab pull in `torque.py` reads bathymetry, so
   once elevations diverge the runs consume different oceanic plates): +0.5% total with
   everything on, +2.0% and rising with erosion or climate off.

4. **Volcanism and the wind model are irrelevant.** Plate movement is the whole story: with
   `simulate_plate_movement` off, node count is exactly flat and land fraction moves -0.002 in
   150 My (all of it erosion).

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

- **Arc accretion at active margins -- crust-building counterweight. Landed 2026-09-02.**
  "An oceanic plate subducting under a continent should make more continent." Two halves,
  both in the torque engine, both gated by the volume-budget cap so neither can revive the
  old land-area runaway:
  - *Areal* (`lithosphere_plate`): a continental plate's *leading* edge growing into space a
    subducting oceanic slab is vacating seeds `ARC_MARGIN_SEED_HC_M`/`_HM_M` (~28/55 km, a
    shallow ~-450 m forearc that builds to land) + `ELEV_CHANGE_SUBDUCTION_ARC`, not
    `growth_seed_thickness`'s drowned oceanic column. The active-margin signal is a node in
    the arc band within `ARC_MARGIN_END_SCAN_NODES` of the growing end, or a recent
    subduction-arc provenance stamp there.
  - *Magmatic* (`rheology.apply_arc_magmatic_thickening`): juvenile Hc added across the whole
    arc *band* -- this plate's continental nodes within `reach_rad` (~500 km) of a converging
    oceanic neighbour, distance-faded from the trench -- **not** yield-gated (unlike
    `apply_convergent_deformation`, which conserves volume and only bites the few-node contact
    line). `ARC_MAGMATIC_HC_RATE_M_PER_MYR = 450`, gentle convergence-rate dependence.
  - **Measured** (seed 926698457, node_density 0.5, climate off, arc-ON vs a rate-0/oceanic-
    seed control, 120 My): land-fraction decline **-0.066 vs -0.073** (~10% slower, +0.006-
    0.008 LF, positive at every checkpoint), mean Hc **+150-250 m**, node count flat (no
    runaway). New: `test_rheology.py::test_arc_magmatism_*`,
    `test_plates.py::test_lithosphere_active_margin_grows_arc_crust_not_ocean_floor` /
    `test_lithosphere_arc_magmatism_thickens_the_continental_margin_band`.
  - Real but modest on its own -- it defends the *margins*. Paired with eustatic sea level
    (below) it stops being marginal.

- **Eustatic sea level -- the big one. Landed 2026-09-02 (`eustasy.py`).** `World.sea_level_m`
  is no longer a fixed input: it's re-solved every step from a **conserved ocean water
  volume** against the world's current hypsometry. Every node has equal area, so ocean volume
  is proportional to the summed water column `W = sum_i max(0, sea_level - z_i)`; `W` is
  snapshotted at generation (`eustasy.initialize_water_budget`, from the flat starting sea
  level) and `eustasy.update_sea_level` (called unconditionally at the end of `step_world`)
  bisects `total_water_column(h) == W` for the new `h`. Deepening a basin (spreading,
  subduction) or drowning a continent raises `total_water_column` at every `h`, so the solved
  `h` drops -- the eustatic fall that hands land back as freeboard, which a fixed sea level
  never did. The `/world/controls` slider now sets `W` to whatever floats the *current*
  hypsometry at the requested level (`eustasy.set_sea_level_via_water_budget`) -- i.e. it
  adds/removes ocean water, which is then itself conserved. Persistence: `ocean_water_column_m`
  is backfilled from an old save's own hypsometry + sea level so loading doesn't jump the
  shoreline. New `World.ocean_water_column_m`; `test_eustasy.py`.
  - **Measured** (seed 926698457, nd 0.5, climate off, 100 My): eustasy alone cuts the
    land-fraction decline from ~-0.052 to ~-0.025 vs a fixed sea level -- roughly halved. Sea
    level falls to ~-150 to -300 m as basins deepen (with visible ~150 m jumps at discrete
    plate-consumption events -- a smoothing pass is a worthwhile follow-up); trajectory goes
    from monotonic decline to slight-rise-then-slow-decline.
  - **All three on** (arc + eustasy + failed rifts) vs a full baseline (none), same seed,
    150 My: land fraction holds **+0.02 to +0.03 above baseline at every checkpoint**, starts
    by *rising* (0.266 -> 0.279 by 25 My) instead of declining immediately, and ends 0.203 vs
    0.182 -- the 150 My decline shrinks from -0.076 to -0.063 and is front-loaded with a gain.
    Plate-count churn stays healthy (14 -> ~21). Chaotic run-to-run divergence makes finer
    attribution unreliable (as the sweep above keeps noting).

- **Failed rifts. Landed 2026-09-02 (`merge_split.RIFT_SUCCESS_PROBABILITY = 0.55`).** A
  continental plate that clears every split gate (flow-fit, pole separation, size, viable
  great-circle cut) now only actually *breaks up* with this probability; otherwise the rift
  arrests (an aulacogen: North Sea, Benue Trough). `LithospherePlate.apply_failed_rift` books
  a one-off thinning (`FAILED_RIFT_THINNING_FRACTION = 0.10`) in a band
  (`FAILED_RIFT_BAND_MULT = 2.5` spacings) around the would-be cut great circle -- a sag
  basin, still thick continental crust, *not* oceanised -- and `reset_age()` puts the plate on
  the normal split cooldown. This directly cuts the rate of the sustained divergent-thinning +
  decompression-melting a *successful* rift inflicts on both daughters' margins. Also hardened
  `apply_failed_rift` against `maybe_split_plate`'s known degenerate-`cut_normal` case
  (spatially-intermingled k-means clusters -> non-unit normal -> skip, no aulacogen).
  Probability tuned to keep the healthy ~18-26 plate-count churn (docs/TODO.md "Plate count
  only decreases"). `test_merge_split.py::test_apply_failed_rift_*` /
  `test_maybe_split_plate_with_failed_outcome_*`.

### Fault mode re-verified; thin belts + weak volcanism addressed (2026-09-04)

**Trigger.** User-supplied `~/Downloads/mantle-bloom-seed52459390-188100000y.mbworld`
(node_density=4, 189 steps / 188.1 My): `land_fraction` 0.122 vs the 0.29 generation target.
Asked to check three things: are collision mountains too thin, is there enough volcanism, and
are glaciers/lakes/rivers/volcanoes building enough plains/plateaus.

**Re-verified the above sections still apply under `fault_deformation_mode`.** Every fix above
was landed/measured before intraplate faults existed (`9bdb7da`, `9fd9287`,
2026-09-03) or before boundary faults covered every plate boundary (`f4b42b9`, this morning) --
not something to just assume still holds. `LithospherePlate.deform()` now scales its
convergent thickening *and* divergent thinning by `fault_influence` in `"fault"` mode
(`lithosphere_plate.py:472`, `:514-516`). Read `faults.generate_boundary_faults`: it rebuilds a
fault family the *entire* length of every plate boundary every step (not a sparse/aging
Poisson layer), so `fault_influence` (80 km reach, 0.06 floor) stays close to 1.0 almost
everywhere the gated deform() logic actually fires -- the gating should be close to a no-op in
aggregate. Confirmed empirically: seed 52459390, node_density=1, land_fraction=0.29/
continental_fraction=0.70 (UI defaults), 150 My, `boundary` vs `fault` mode:

| config | 0 My | 50 My | 100 My | 150 My |
|---|---|---|---|---|
| `boundary` | 0.2900 | 0.2923 | 0.2404 | 0.2181 |
| `fault` (default) | 0.2900 | 0.2947 | 0.2396 | 0.2081 |

Within ~0.01 at every checkpoint -- noise-level, same magnitude as this section's own
baseline-vs-volcanism-off spread above. The 2026-09-03 default change is not an additional
regression on top of what's documented above.

**Thin mountains -- confirmed, fixed.** On the trigger save, `elev_change_reason` land-node
share: `ELEV_CHANGE_COLLISION` 2.4%, `ELEV_CHANGE_COLLISION_FAR_FIELD` 0.0% (never appears).
Root cause: `orogen_dilation_nodes` (the near-field ring beyond the bare contested/overlap
band) was `round((orogen_reach - 1.0) * COLLISION_REACH_DILATION_NODES_PER_UNIT) if orogen_reach
> 1.0 else 0` -- a hard zero at `collision_uplift_reach_multiplier`'s own untuned default
(1.0). A real orogen is a broad crumple zone (Himalaya ~500 km), not just the suture line.
Fixed by making the ring linear in the knob from 0 (`round(orogen_reach *
COLLISION_REACH_DILATION_NODES_PER_UNIT) if orogen_reach > 0.0 else 0`) -- the default now
carries a real +2-node ring each side, not none, without touching the knob's own default value
(TUNING_MULTIPLIER_FIELDS' "1.0 == untuned" contract, `test_tuning_knobs.py`, is a base-rate
change, not a default-value change). Also found, and left alone: `fault` mode's own boundary-
fault relief layer (master trace + 2 strands 40/80 km inland, `faults.py`) already gives a
second, genuinely position-distinct uplift source up to ~125 km inland that `boundary` mode
never had -- but at `BOUNDARY_FAULT_RELIEF_SCALE * REVERSE_UPLIFT_M_PER_MYR` = 66 m/Myr, ~8% of
the core `CONVERGENT_MOUNTAIN_RATE_M_PER_MYR` (800) rate, it's a minor contributor, not the
fix. Measured (seed 559394024, `test_tuning_knobs.py`): `reach=3.0` land fraction 0.2900 ->
0.2137 at 150 My vs 0.2081 unmodified `fault` default -- a real but modest +0.006 nudge, a
shape/width fix for this specific complaint, not a reversal of the decline documented above.

**Insufficient volcanism -- confirmed, more severe than expected, fixed.** Trigger save: 112
volcano nodes ever spawned across 188 My / 24 plates, 0 still active, only 17 above sea level,
mean elevation of every volcano node -4583 m. Volcano nodes only ever spawn where
`LithospherePlate.deform` finds a rift boundary stretched too thin (`volcanism.py`'s own
docstring) -- never at a convergent (arc) or intraplate (hotspot) setting, Earth's two dominant
land-building volcanic regimes; no code changes for that here, it's a materially bigger change
than this pass's scope. What *is* fixed: `ERUPTION_ELEVATION_M` 100 -> 300 m,
`ERUPTION_RATE_PER_MYR` 3.0 -> 5.0 (`elevation_lines.py`/`volcanism.py`) -- a volcano that does
erupt now tops out around 1500 m of gross relief instead of ~300 m, which used to read as
erosion noise within a few steps. Measured: `volcanism_multiplier=5.0` on top of this tracks
*identically* to the unmodified baseline through 50 My in aggregate land fraction (0.2947 both)
-- expected, given how few nodes are ever volcanic; this is a fix for visible peaks, not one
expected to move the aggregate number.

**Not enough plains/plateaus -- not an independent mechanism; a symptom of the above two.**
There's no dedicated "build a plateau" step, and this section's own erosion-half fix already
established erosion/deposition as land-fraction-neutral over a long run -- it just planes down
whatever relief tectonics/volcanism created. 52% of the trigger save's land last changed via
erosion, 18% via deposition, only 2.4%/~0% via collision/volcanism -- nothing was actively
building broad elevated terrain. The reach-multiplier fix above directly addresses this: its
`COLLISION_REACH_NEAR_FIELD_FACTOR` (0.4x rate) near-field ring *is* a lower-relief, plateau-
like apron around every belt, so no separate mechanism was added for this specifically.

**New: Elevation view "Mountains" / "Plains & Plateaus" legend toggles**
(`render_image._classify_terrain_relief`, `/world/render`'s `show_mountains`/
`show_plains_plateaus`). Local relief (elevation range within `TERRAIN_RELIEF_RADIUS_KM` = 50
km, a real-world radius so it reads the same across `node_density`), not raw elevation, decides
the split -- every land node is exactly one of the two (mountain: relief >=
`MOUNTAIN_RELIEF_THRESHOLD_M` = 400; plains/plateau: below), matching how real geomorphology
tells a rugged range from a high flat plateau. Computed only for the Elevation view, only when
at least one toggle is on (`_render_grid_arrays`' `include_terrain_relief` param, opt-in,
every other caller unaffected) -- baked server-side as a translucent wash over the ordinary
hypsometric colour, not a flat swap, so the elevation itself stays legible underneath.

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

`LithospherePlate.deform()` does inline line growth/shrinkage and regularization per turn.
`PlateWithRTree` -- the R-tree-backed variant -- doesn't carry equivalent regularization
logic. If `PlateWithRTree` is meant to become a drop-in replacement, it needs its own
version of the per-turn node density / spacing upkeep (`elevation_lines.py`,
`TARGET_LINE_SPACING_RAD`), otherwise its lines drift out of spec over a long run.

---

## Intraplate faults: follow-ups

**Status:** shipped (`faults.py`, `World.faults` / `World.fault_systems`; "Fault lines" map
view). See `docs/simulation-model.md#faults`. The layer is additive -- it never touches
`deform()`'s own boundary classification -- so it can be tuned or reverted in isolation.

**6. Faults along every plate boundary + motion-based boundary classification. DONE (2026-09).**
`faults.generate_boundary_faults` now lays a fault family along every classified stretch of
every plate boundary each step (`World.boundary_faults`, rebuilt from scratch -- they track
the moving edge, so no ageing / scars / reconcile), classified by the local closing rate:
convergent → reverse, divergent → normal, transform → strike-slip. This is what makes the
"Last elevation change" view show fault structure along the boundaries and gives
`fault_influence` something dense to gate on (see 3a -- now moot). In the same change,
`torque.classify_boundary_nodes` switched from geometric (`contested` overlap = the only
"convergent") to **motion-based** (closing-rate sign over the near-boundary band), so a
converging boundary builds an orogen before it visibly overlaps and transform boundaries are
their own class with a real (gentle) `TRANSFORM_UPLIFT_RATE_M_PER_MYR` term. Geometric
`contested` is kept only as the node-deletion / continental-retreat trigger.
Open follow-ons: boundary-fault `fault_id`s churn every step, so the "Plates & Faults" click
inspector can't reliably target one; `BOUNDARY_FAULT_RELIEF_SCALE` (0.3) and the transform
rate are eyeballed and fold into item 3's calibration sweep.

**1. Trace and fault-system lengths. DONE (2026-09).** Added `FaultSystem` as a first-class
object one level above the lone fault (`World.fault_systems`, `SYSTEM_SPAWN_FRACTION` = 0.18
of spawns): a long, gently curving master lineament (`SYSTEM_LENGTH_*`, to
`SYSTEM_LENGTH_MAX_KM` = 5500 km) warped by a couple of low-frequency lateral lobes, with a
5–16-strand sub-parallel family scattered along it -- `SYSTEM_STRAND_SPACING_KM` (65 km)
along strike, `SYSTEM_BELT_HALF_WIDTH_KM` (130 km) across -- each strand a `Fault` carrying
the `system_id`, drawn from a **widened** length distribution (`SYSTEM_STRAND_LENGTH_*`,
median 150 km, tail to 1300 km). The master trace is a scaffold: it applies no relief, the
strands do. Systems outlive their strands (`SYSTEM_LIFESPAN_MIN..MAX_MYR`, 25–140 Myr) then
become inert scar bundles (`_cull_inactive_systems`). The lone-fault / tight-set path
(`LENGTH_MAX_KM` = 200 km, `SET_*`) is unchanged -- a system is layered *on top*, not a
replacement. Rendered in the "Fault lines" view as a shaded belt + dashed master centerline
under the strand family (`FaultInspector.tsx`), with `system_id` in `/world/faults` and the
inspector panel. Not done here: the master lineament doesn't spawn *fresh* strands over its
lifetime (the family is all born at once); and item 3's calibration sweep now matters more
(systems add a lot of cumulative strand relief).

**2. The trace applies relief but doesn't shear the node field. PARTLY DONE (2026).**
`_apply_plate_fault_relief` raises/lowers `elevation` near an active trace per regime, but the
`ElevationLine` nodes on either side are never actually displaced *along* the fault -- a
strike-slip fault with tens of km of `cumulative_offset_m` leaves piercing points (a river
valley, a ridge crest) exactly where they were. Real strike-slip offset of pre-existing
features is the visually recognisable thing about a fault like the San Andreas. Doing that
properly (moving nodes tangent to the trace by the along-strike slip) is still open.

What *did* land: `World.fault_deformation_mode` (`"fault"` **default** / `"boundary"` /
`"both"`, Controls window). In `"fault"` mode `LithospherePlate.deform` gates its own
convergent/divergent thickening by `faults.fault_influence()` (distance to an active fault
trace, `FAULT_DEFORM_REACH_KM` 80 km / floor 0.06) and `_apply_plate_fault_relief` widens
its reach, so plate-boundary deformation reads as a segmented belt tracking the fault
families instead of a smooth band at the polygon edge. For that to work faults have to be
*on* the boundary: seed placement is now a separate, sharply boundary-peaked kernel
(`SPAWN_PLACE_DECAY_LEN_KM` 200 km / floor 0.004) decoupled from the broader stress weight
that sets the spawn *count* -- without it the plate's huge interior-node count dragged the
median seed ~600 km off the boundary. Faults also spawn wherever two plates' node clouds
overlap (`OVERLAP_STRESS_WEIGHT` lifts both weights on overlapped nodes). `FAULT_RELIEF_MODE_RATE_SCALE`
was dropped 3→1: with faults hugging the contact it barely gets gated there, so a >1 rate
double-counted the bands and pinned hypsometry at `MAX_ELEVATION_M`. See
`docs/simulation-model.md#faults`.

**3. Spawn rate / relief magnitudes are eyeballed, not calibrated.** `BASE_SPAWN_RATE_PER_MYR`
(3.0 systems/Myr sphere-wide at full stress) and the per-regime `*_M_PER_MYR` relief
constants were picked to look plausible and stay well under the boundary rates, not measured
against a target fault density or a hypsometry budget. Worth a sweep: fault count vs.
`elapsed_years` across seeds, and whether the cumulative fault relief moves the land-fraction
/ hypsography numbers the [land-fraction decline](#land-fraction-slowly-declines-over-a-long-run)
work cares about.

**3a. Sparse-fault worlds: `"fault"` mode ≈ `"boundary"` mode. MOOT (2026-09, item 6).**
Boundary faults (`generate_boundary_faults`) now line every boundary regardless of the
intraplate spawn rate, so `fault_influence` gates the bands onto real fault families
everywhere -- there is no longer a "sparse-fault seed" where the mode collapses to
`"boundary"`. The original text: with only ~20 active intraplate faults the boundary-hugging
traces covered the few short contested zones densely enough that `fault_influence` ≈ 1 along
them, so the gating was a near no-op.

**4. No live tuning knob.** Unlike the geomorphic budget, none of the fault *constants* are
exposed in the Controls window (still true). `World.fault_deformation_mode` -- the
`"boundary"` / `"fault"` / `"both"` selector added for item 2 -- *is* a Controls select now,
but that's a model switch, not a magnitude knob. If the fault relief / spawn constants turn
out to matter for a world's look they should still join the live knobs.

**5. Faults don't influence seismicity / hazard output. DONE (2026).** `faults._generate_earthquakes`
emits one characteristic `Earthquake` per active fault per step once the fault has slipped at
least `MIN_STEP_SLIP_FOR_QUAKE_M` (a real active fault ruptures thousands of times per Myr --
one representative event is all we keep). Magnitude is from trace length + slip rate + an
overlap bonus; epicentres are transient world-frame points in `World.earthquakes`, pruned
after `EARTHQUAKE_RETAIN_MYR` (~5 Myr). `erosion.py` reads them for a local seismic-erosion
burst (`_earthquake_erosion_multiplier`, `EARTHQUAKE_EROSION_*`); `GET /world/earthquakes`
and the fading epicentre overlay in the "Fault lines" view expose them. The largest event
each step (if `>= EARTHQUAKE_LOG_MIN_MW`) is logged to the event console.

---

## Very-long-run collapse (399 My): confirms the "still open" items above are the live bottleneck

**Context.** `~/Downloads/mantle-bloom-seed920135003-399300000y.mbworld` (132 steps, ~3 Myr/
step, `node_density=4`) is the longest-run save inspected so far (prior investigations above
top out around 254.8 My). Five symptoms were reported against it (missing elevation over
most of the world, plate 0 overlapping several small oceanic plates, plates 0/1 being huge,
diagonal-stripe artifacts on the elevation/biome maps, and lakes appearing inside the ocean).
Loading it and cross-checking against the "Plate geometry degrades on long runs" section
above shows these are the *same, already-tracked* still-open gaps, just further along:

1. **~42% of the sphere has zero elevation nodes.** Rasterizing every plate's
   `all_points_and_elevation()` onto a 1x1 degree grid: 27,576 / 64,800 cells empty. The
   `plates` render view (nearest-plate-by-frame ownership) still fills the whole sphere with
   no gaps, but `platesDetail` (actual node clouds) shows enormous dead zones -- so plates 0/1
   still nominally *own* most of the globe by Voronoi territory, they just have no populated
   rows there. This is exactly item 6's open half ("there is still no mechanism to spawn a new
   plate for a large region a consumed oceanic plate vacated (the old `gaps.py` fallback,
   never ported)") plus mechanism 4 ("periodic conservative continental re-lattice", also
   still open) -- both were left unaddressed because the churn repros tested (<=180 steps)
   stayed healthy without them. At 132 steps the oceanic fleet has been ground down to 7 tiny
   plates (361-6,850 nodes each, 10,932 nodes total, vs 83,166 continental) and nothing
   refills the vacated ocean floor. **The spawn half is fixed 2026-09-04** -- see item 6's
   update above. On this exact save, `gaps.fill_gaps` finds one dominant 40,474-lattice-point
   void (plus a handful of <=101-point residuals well under `MIN_GAP_NODES`, left alone) and
   fills it with a new 43,004-node oceanic plate 23, dropping empty 1x1-degree coverage from
   42.5% to 13.5% in one pass -- the residual is mostly ordinary sparse polar-lattice
   coverage, not further void. The re-lattice mechanism (4) remains open for the continental
   side of the ratchet.
2. **Plate 0 overlaps 6 different oceanic plates (3, 7, 11, 16, 20, 22); plate 1 overlaps 2
   more (17, 21).** `merge_split._deep_continental_overlap_fractions` /
   `update_overlap_progress` / `pop_ready_forced_merge` only ever look at *continental-
   continental* pairs (`if pid_a not in continental: continue`, both sides) -- there is no
   forced-resolution backstop for a continental-oceanic overlap that's stalled the same way,
   only the per-step retreat/subduction path. These are exactly the 7 surviving oceanic
   remnants from (1): they're wedged against the one continental landmass because that's the
   only place they still have any territory, and whatever is throttling their subduction
   retreat (per-step `n_distance_cap`, or simply not being detected as `contested` against
   such a large neighbour) hasn't finished consuming them after 132 steps. Worth an offline
   `plate_diagnostics.py` pass on this save to see whether these nodes are even reaching
   `torque.classify_boundary_nodes`'s `contested` mask, or falling through some other gap.
3. **Plates 0 (21,503 nodes) and 1 (61,706 nodes) are supercontinent-scale but *not* stuck --
   they're recent merge products.** `plate.age_steps` is 18 and 15 respectively (`reset_age()`
   fires on merge), against `SPLIT_MIN_AGE_STEPS = 20`; both already clear
   `maybe_split_plate`'s residual-fit and pole-separation gates by 3-6 orders of magnitude
   (checked directly: `_fit_residual_rms` ~5.4e-9 vs a ~1.4e-9/9.2e-10 threshold). So this
   looks less like "should have split and didn't" and more like a snapshot mid-cycle, 2-5
   steps from being eligible to rift again -- consistent with the documented 15-26-plate
   merge/split oscillation, just caught at the "just merged" end of it. The still-open
   k-means split-cluster-quality problem (poles fit from spatially-intermingled velocity
   clusters) is the likely reason two plates this size keep re-forming via merge in the first
   place, rather than cutting cleanly the first time.
4. **The diagonal-stripe artifacts on the elevation/biome/platesDetail renders** are the
   documented "streaking" (item 3 in the plate-geometry section): a triangular tongue grown by
   repeated row-end extension, its staircase edge rendered as a long regular sawtooth. Not a
   new bug -- just visually worse here because of how much small-plate clutter (2) is packed
   into the one populated corner of the map.
5. **Lakes rendered inside open ocean -- new hypothesis, not yet confirmed.** Not investigated
   deeply (deferred per request), but `hydrology.connected_ocean_mask` classifies "ocean" as
   the single largest connected below-sea-level component of the actual node k-NN graph (plus
   any second component within `SECOND_OCEAN_SIZE_FRACTION` of its size) -- anything smaller
   is treated as an ordinary endorheic basin and handed to `lakes.py`. A below-sea-level patch
   on the far side of one of (1)'s node-cloud gaps would be topologically severed from the
   main ocean in that k-NN graph even though it's geographically part of the same sea, so it
   would misclassify as a giant "lake" rather than ocean. If so, this is a downstream symptom
   of (1) rather than an independent bug in `lakes.py` itself -- worth checking first before
   touching lake code directly.

**Net read:** nothing here is a new bug distinct from what this doc already tracks; it's
confirmation that (6)'s new-plate-spawn gap and mechanism (4)'s re-lattice were the load-bearing
missing pieces once a world runs long enough, and a plausible (unconfirmed) mechanism linking
the lake-in-ocean report to the same node-cloud gaps. (6)'s spawn gap is now fixed (`gaps.py`,
2026-09-04); mechanism (4), items 2/3 (the small-oceanic-plate overlaps and the resulting
supercontinent-scale plates 0/1), and the lake-in-ocean hypothesis are all still open --
plausible next candidates now that the dominant void itself no longer masks them on a
re-run of this save.

---

## `gaps.py`'s plate-spawn is a stopgap, not the real fix

**Status: known hack, flagged at merge time (2026-09-04).** `gaps.fill_gaps` (see the
"Very-long-run collapse" section above and item 6 further up) closes the immediate
coverage hole, but the mechanism itself doesn't match how new ocean floor actually forms.

**Why it's a hack.** Real ocean floor is produced *continuously*, along a mid-ocean ridge,
as two already-existing oceanic plates pull apart and decompression melting accretes new
crust onto both of their own divergent edges -- new sea floor is always born attached to a
plate that's already there and already moving. `gaps.py` instead runs as a periodic,
whole-sphere patch: it notices a big enough hole *after the fact* (every
`GAP_FILL_INTERVAL_STEPS`), then conjures an entire fully-formed plate into existence out of
nothing, with an Euler pole fit post-hoc from the local mantle-flow field rather than
inherited from ridge-push off any real neighbour. It works (see the fixed save,
node coverage 42.5% -> 13.5% empty in one pass) but it's a bookkeeping fix for a modeling
gap, not a simulation of the actual process -- an artifact of `deform()`'s own per-step
boundary growth only ever extending a line from a node that already exists (see `gaps.py`'s
own module docstring), so a region with *no* nearby line at all has structurally no way to
ever grow back on its own.

**Direction.** The real fix is upstream of `gaps.py` entirely: give `LithospherePlate.deform()`
(or a sibling mechanism) a way to originate new oceanic crust directly at a divergent
boundary -- i.e. model mid-ocean-ridge spreading itself, rather than letting the void get big
enough that a whole-sphere sweep has to notice it later. That would mean an oceanic plate's
own divergent edge can always keep pace with its neighbour's retreat (no ordinary rift ever
outruns per-step growth in the first place), and the only remaining role for something like
`gaps.py` would be a true edge case (e.g. a divergent boundary between two plates that both
fully vanished the same step) rather than the routine, sizable voids seen on long runs today.
Whoever picks this up should start from `_grow_or_shrink_line_for_deform` /
`_claim_adjacent_territory` in `lithosphere_plate.py` -- the existing per-line growth this
would need to extend -- and treat `gaps.py` as the thing this work should let shrink back to
a rare fallback, not as the permanent mechanism.

**Partially addressed 2026-09-04 -- the type decision, not the underlying mechanism.** This
stopgap's own hardcoded "oceanic" -- regardless of where the gap actually sits -- is fixed:
`_spawn_plate_from_gap` now decides per-node type by real local context (continental only
where a gap point genuinely hugs a still-standing continental coastline, see
`GAP_LAND_ADOPTION_RADIUS_MULT`; oceanic everywhere else, unchanged from before), and the
spawned plate's own `crust_type` label is the actual majority of what it ended up with, not a
hardcoded constant (see docs/simulation-model.md#per-node-crust-type). Separately, the *local*
decompression-melting mechanism this section's own "Direction" gestures at
(`rheology.apply_divergent_deformation` + `LithospherePlate.deform`'s existing melting
handling) now also decides oceanic-vs-continental by whether the melting node was still above
sea level, rather than always resetting to the oceanic reference regardless of context -- see
docs/simulation-model.md's "Whole-sphere coverage" section. **Still open, unchanged by this
work:** `gaps.py` itself is still exactly the periodic, whole-sphere, conjure-a-fully-formed-
plate stopgap described above -- nothing here touches its cadence, its post-hoc Euler-pole
fit, or its fundamental "notices a hole after the fact" character. The real fix this section
calls for (ordinary end-growth structurally never able to outrun a large-enough rift) is still
entirely open.

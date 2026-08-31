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

**Status:** all three sub-bugs **fixed** 2026-08-30 (details below); a long-run re-verify
against a fresh >150 My save is still worth doing to confirm the node-count blowup is
actually gone rather than just slower, and the k-means split-cluster quality noted under
bug 2 / bug 3 is still open.

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

**Status:** investigated 2026-08-31 from
`~/Downloads/mantle-bloom-seed888151728-85000000y.mbworld` (seed 888151728, 851 steps /
85.1 My, `node_density=4`, `climate_density=4`). Interim render-only mitigation **applied**
(see "Interim mitigation" below).

Option 2 (coastal planation + infill feedback) **implemented + validated 2026-08-31** -- see
"What landed" below. Driven against the seed-888151728 save, 25 steps (2.5 My):

| in the 7-17 N / 3-15 E box | feedback on | feedback off |
|---|---|---|
| isolated near-sea-level nodes | 21 -> 8 | 21 -> 19 |
| world-wide speckle fraction | 6.6% -> 4.6% | 6.6% -> 6.1% |

Visually: the sheltered north coast and the south spit resolve to a clean shoreline; the
broad central promontory converts to a coherent marsh/intertidal mosaic rather than a
checkerboard (arguably correct -- it *is* a drowned coastal plain). Pushing the coefficients
harder (`PLANATION_RATE` 400-600, sharper prominence/shelter) made the dither *worse* --
amplifying node-scale openness noise into +-10 m elevation swings -- so the mild committed
values are the right operating point.

The interim render hack (option 1) is **still in place** but is now nearly a no-op after the
feedback runs (isolated specks it targets are mostly gone; the residual is contiguous marsh
fringe it can't grab). Deciding whether to delete it, keep it as belt-and-braces for the
first few steps after loading an old save, or replace it with a milder pass, is the only
open thread. A land-locked-coastal-pit infill sink (an isolated sub-sea-level node ringed by
land, which is neither `is_ocean` nor above sea level, so both planation and the ocean-sink
infill skip it) was prototyped and reverted -- it didn't clearly help the metric and the
planation/infill asymmetry it exposed needs more thought.

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
- Nothing anywhere looks at coastal *connectivity* ("this ocean cell is nearly landlocked",
  "this land cell is nearly surrounded by water"). `grep` confirms no planation,
  progradation, spit, barrier, or lagoon logic.

Net: a marginally-submerged flat sheet is a stable fixed point that just dithers forever.

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
   land-ocean field flips isolated 1-cell specks to match their surroundings. Kills the
   checkerboard in every view, zero physics risk -- but the simulation still holds a flat
   sheet balanced on the waterline, and the Plate Inspector / raw nodes still show it. This
   is the **interim mitigation** still applied (see below).

2. **Coastal planation + infill feedback (the real fix)** -- **implemented + validated at
   density 4**, see "What landed" and the Status table above. Cuts the node-level dither
   roughly in half; the render hack (option 1) is now nearly redundant on top of it.

3. **Barrier islands (the user's framing)** -- **delivered in emergent form** (see "What
   landed"). A genuine wave-exposure / fetch field plus an explicit longshore-transport
   direction (for real spits and drift-aligned bars, not just fetch-sheltered accretion)
   remains an optional future refinement.

**Interim mitigation (applied 2026-08-31 -- now nearly a no-op on top of option 2; decide
whether to keep, drop, or soften it).**
`render_image._despeckle_coastal_elevation` -- a render-only pass that snaps an isolated
near-sea-level node (within `_DESPECKLE_BAND_M` of sea level, and with at least
`_DESPECKLE_MAJORITY` of its `_DESPECKLE_NEIGHBORS` nearest near-sea-level neighbours on the
opposite side of the waterline) to its neighbour-median elevation, so it renders as part of
the surrounding land or ocean. Restricted to the near-sea-level node subset so a genuine
steep coast (whose land nodes climb out of the band immediately) is untouched. Called on
the gathered `all_elevation` in `_render_grid_arrays` (Elevation view), `_biome_fields`
(Biome / Combined), and `_resource_fields` (Resources / Soil Quality) before their
nearest-node resample; never touches `world.plates`. Covered by
`unit_tests/test_render_image.py` (`test_despeckle_*`). If it goes: delete the function, its
three constants, the three call sites (`_render_grid_arrays` / `_biome_fields` /
`_resource_fields`), and those three tests + the `_latlon_patch_points` helper. Post-option-2
its effect is marginal (before/after renders of the density-4 save are near-identical with
and without it), so this is now a judgement call, not a blocked task.

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

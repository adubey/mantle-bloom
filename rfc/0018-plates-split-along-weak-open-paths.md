# RFC 0018: Plates split along weak paths that are pulling apart

- **Status:** Accepted (revised 2026-11-14 after the #894 exploration of the D1 prototype; revised 2026-11-14 after #848 runtime evidence and #876 rift-candidate exploration; revised 2026-11-09 after spectral-weighting, runtime and failed-rift review; revised 2026-11-08 after G2 sample-size review; revised 2026-11-06 after G2 non-worsening criterion review; revised 2026-11-06 after G2 metric review; revised 2026-11-04 after G1 gate review; owner confirmed 2026-11-03)
- **Date:** 2026-11-03
- **Revisions:** dated entries, newest first, are in the [Revision note](#revision-note) near the end.
- **Related:** #487, #699, #770, #797, #802, #836, #846, #848, #850, #851, #866, #876, #888, #894; RFC [0007](0007-mantle-plumes-and-plateaus.md) (plumes); RFCs [0009](0009-craton-mosaic.md) (cratons-first worlds), [0010](0010-rifts-follow-weak-lithosphere.md) (accepted predecessor; withdrawal follows gates)

## Problem

When a plate splits today, the model divides it with a great circle. The line is chosen from
the difference between two mantle-flow patterns, but it does not follow the plate's weak
places or check whether the two new plates would move apart. It can also leave a very small
piece. A plate's interior never stretches before it splits. Extension happens only at plate
margins, where `quad_tectonics._open_rift` thins the cells within two rows of newly inserted
cells ([#876 step 3](https://github.com/adubey/mantle-bloom/issues/876#issuecomment-6874155833)).
Rifts therefore do not form, grow, fail or leave scars inside plates. They appear as
instantaneous cuts.

Two earlier fixes did not solve this. RFC 0010's minimum-cost cut did not route around cratons
reliably. Issue #770's inherited pseudo-sutures did not reduce tiny breakups: the 0.2 and 0.4
route-cost trials produced more tiny outcomes and no seed improved
([#770](https://github.com/adubey/mantle-bloom/issues/770)).

The graph-cut method that this RFC first proposed carries a serious runtime risk and does not
measure physics:

- **Runtime.** Eigensolves are about 85% of its search time, with up to 97 sparse solves per
  decision. 94.6% of mature-world opportunities have no structurally valid candidate, yet a
  perfect exact skip of those still leaves about 1.84 times the baseline step time against G3's
  2% limit. Lower refinement caps change partitions and `P_break`. The only exact shortcut found,
  stopping at the first repeated partition, saves 16% of solves. Untested output-preserving
  optimizations (warm-started solves, candidate reuse across steps) are not ruled out, but none
  tested comes near the limit
  ([`analysis/issue848/README.md`](../analysis/issue848/README.md)).
- **Validity is geometry, not load.** On 125 opportunities (1,491 candidates), neither a
  drag-based load map nor a thin-shell stress map separates valid from invalid cuts (rank AUC
  0.47–0.56). Plate size does (0.86). Breakup probability among valid cuts does follow load
  (AUC 0.89 on the thin-shell map) ([`analysis/issue876/README.md`](../analysis/issue876/README.md)).
- **One step's load is far below rock strength.** Along valid cuts the median thin-shell load is
  about 7% of column capacity. RFC 0010 Gate 2 found about 12% on its cuts
  ([`analysis/issue487/README.md`](../analysis/issue487/README.md)). Intact rock does not fail
  under a single step's load.

The offline damage experiments for RFC 0010 Gate 3b
([#699](https://github.com/adubey/mantle-bloom/issues/699)) let damage accumulate under slab pull
and mantle currents. They produced broad damage that fragmented plates into hundreds or
thousands of pieces, never two daughters. The diagnosis recorded there was that damaged rock
did not weaken, so a crack could not concentrate load; failures that were not kept also healed.
Issue [#866](https://github.com/adubey/mantle-bloom/issues/866)'s wide below-sea-level basins,
which #836 traces to margin thinning, show the same missing feedback from the other side:
extension spread too wide.

The first prototype of this RFC's damage bands tested each cell against its own strength. On
development worlds it made no split in 60 Myr and nibbled at plate edges
([#894](https://github.com/adubey/mantle-bloom/issues/894),
[`analysis/issue894/README.md`](../analysis/issue894/README.md)):

- **The load is flat and weak rock is at the edges.** The median extensional load is 0.50 TN/m
  at a plate's edge and 0.58 TN/m ten cells in, about a tenth of strength. Median strength is
  4.3 TN/m at the edge and 9.7 TN/m ten cells in. Only 3–6% of cells are loaded past their
  strength, and 66% of straining cell-steps are oceanic.
- **Softened damage spreads into patches.** Bands of 50 or more cells are a median 3 times
  longer than wide, lie a median 31° from the tension direction, and reach about 900 by 500 km.
  A cell's load does not fall when its neighbour yields, so softening widens the patch.
- **Whole sections are pulled hardest across the middle.** The straight section carrying the
  largest total pull leaves at least 20% of the plate on its smaller side in 67–80% of plates
  on five development worlds. Between snapshots 10 Myr apart the pull direction stayed within
  30° in 47 of 56 cases.
- **The kept daughter minimum excludes most large splits.** `SPLIT_MIN_NODES = 700` applies to
  each daughter. The median plate has about 1,700 cells, so only near-even splits pass, and a
  split into two large but unequal plates (for example 1,190 and 510 cells) is refused.

**Design principles:** the great-circle cut and the instantaneous graph cut are geometry
shortcuts that do not let rifting follow its own physics (principle 8). The graph cut's runtime
conflicts with principle 2. This proposal keeps extension, weakening and thinning in physical
units (principles 4 and 5) and judges the result on multiple seeds and visible maps
(principle 9).

The project owner's goal for this revision ([#876 decisions](https://github.com/adubey/mantle-bloom/issues/876#issuecomment-6874579897))
is a rifting process that runs efficiently and makes realistic-looking planets. It is not a
proof of new rift mechanics. Rifting on Earth is well documented, so the model should reproduce
its look and timing cheaply. Weakly constrained quantities are parameters, with literature
defaults where they exist and a tuning search against Earth targets.

## Proposed solutions

**What nature tells us.**

- Rifts reuse old weak zones. Of eighteen Atlantic and Indian Ocean conjugate margin pairs, all
  but the Madagascar–Seychelles–India system opened along former sutures
  [[Buiter & Torsvik, 2014]](#ref-buiter2014). Continents tend to break up parallel to ancient
  orogenic belts [[Vauchez et al., 1997]](#ref-vauchez1997).
- Rifts have a long slow phase, under 10 mm/yr full rate, then an abrupt increase in divergence
  before continental rupture. Brune et al. attribute the increase to the nonlinear decay of rift
  strength [[Brune et al., 2016]](#ref-brune2016). In model terms, weakening feeds back on
  extension.
- Rifts localize. Extension may concentrate in a narrow zone or spread over a region wider than
  the lithosphere is thick, depending on crustal thickness and heat flow [[Buck, 1991]](#ref-buck1991).
  Numerical rift models produce strong localization by letting frictional strength decrease
  with accumulated strain [[Huismans & Beaumont, 2003]](#ref-huismans2003).
- Many rifts fail. In failed rifts the force balance becomes dominated by lithospheric
  strengthening, and activity stops before rupture [[Brune et al., 2023]](#ref-brune2023).
  Geological inventories document failed systems beside rifts that reached seafloor spreading
  [[Elling et al., 2022]](#ref-elling2022).
- Thinned lithosphere cools back toward its original state on a thermal time constant of about
  62.8 Myr for 125 km lithosphere [[McKenzie, 1978]](#ref-mckenzie1978).
- Plumes help. Twelve of the eighteen margins are linked to large igneous provinces, with
  breakup at the same time as the province [[Buiter & Torsvik, 2014]](#ref-buiter2014).
  Courtillot et al. argue that a plume is a prerequisite for opening a major ocean basin
  [[Courtillot et al., 1999]](#ref-courtillot1999). Some margins broke up without one.

**Alternatives.**

1. **Keep the great-circle cut.** Simple and clean, but it ignores weakness, opening and rift
   history, and makes interiors rigid until the cut.
2. **Keep RFC 0010's detected, routed rifts.** Detailed setting and route gates, but #770 found
   that weaker pseudo-sutures did not reduce tiny breakups, and it shares the instantaneous-cut
   limitation.
3. **The graph-cut candidate search this RFC first proposed.** It chooses a well-formed cut, but
   its validity checks measure geometry, its eigensolves cost several times the step budget, and
   no exact shortcut brings it within G3 (#848).
4. **Damage bands grown cell by cell.** Every cell carries accumulated rift strain and strains
   when its own load exceeds its own softened strength; a band that happens to cross the plate
   splits it. This was the 2026-11-14 design. Its prototype strained the weakest rock, spread
   into patches and never split a plate (#894).
5. **Rift lines (recommended).** Each plate's forces already say where it is being pulled apart
   as a whole: the section carrying the largest total pull. Within a strip around that section,
   the cheapest path through the rock's strength is the plate's rift line. A rift nucleates at
   the line's weakest point and its tips grow along the line, straining, thinning and softening
   only the line's cells. The plate splits when the intact remainder of the line can no longer
   carry the pull. A line whose pull drops stops growing and leaves a scar that heals partly
   and can reopen. The rock still carries all the memory, and the cost is one path search on a
   strip.

Alternative 5 matches the documented behavior: slow then fast, along inherited weakness, often
failing, aided by plumes. It keeps the #699 lesson that damage must soften rock, and adds the
#894 lesson: failure has to be judged for a whole section, and strain has to be confined to
one line, or damage follows the weakest rock and spreads. It also addresses principle 2's
runtime constraint.

## Design

### Breakup bookkeeping kept from the accepted method

- Minimum plate age (`SPLIT_MIN_AGE_YEARS`) and split cooldown remain hard gates.
- Every applied split goes through a mask-based split helper that calls
  `_plates_from_node_masks()`, validates daughter sizes and ownership, and resets daughter ages.
- A plate must still have at least `SPLIT_MIN_NODES` cells (scaled by node density) to split.

Two production rules change (owner decision, 2026-11-14, recorded in the Revision note):

- **Daughter minimum.** The smaller daughter must hold at least `a_min = 15%` of the parent's
  area. This replaces the requirement that each daughter have `SPLIT_MIN_NODES` cells, which
  allowed only near-even splits on typical plates.
- **Stray pieces.** A wiggly line can cut off small pieces that touch neither end of the plate.
  Pieces of at most `a_stray = 1%` of the parent's area each join the daughter they share the
  most edge length with. A larger extra piece still cancels the split.
- The success path's craton-loss diagnostic (`cratons.record_rift_cut`) and belt marking use one
  trace-derived band around the split path. No great-circle band remains in the split path.
- The generated starting world is not changed.

The candidate search, its `O/R/S` breakup prior, proposal distribution, momentum score, keyed
breakup draws and great-circle failed-rift thinning are removed. Failure is now a band that stops
growing, not a sampled event.

### State in the rock

Three new per-cell surface fields move with the rock:

- `rift_strain` (dimensionless, ≥ 0): accumulated extensional thinning strain, `−ln β` for the
  cell's own thinning. It never decays and records the permanent, fault-like part of the damage.
- `rift_warm_strain` (dimensionless, 0 ≤ `rift_warm_strain` ≤ `rift_strain`): the share of that
  strain whose thermal weakening has not yet healed. It decays with time constant `τ_heal`.
- `rift_broken` (boolean): the cell has broken through (step 4). It is a separate latch, so
  `rift_strain` keeps meaning thinning strain for softening, `ε_cand` and the recorder.

**Remap.** Adaptive coarsening and plate merges must not heal damage, so area-weighted averaging
is wrong for these fields: one broken cell coarsened with three intact siblings would drop from
`ε = 2.0` to 0.5. `rift_strain` coarsens and merges by maximum, like `channel_depth` and
`breach_notch_depth_m`, which are special-cased in `sparse_quad_patch` coarsening and
`quad_merge`. `rift_warm_strain` is intensive (area-weighted), then clamped to at most the
result's `rift_strain`. `rift_broken` is `BOOLEAN_PROVENANCE` (any). Refinement copies all three
to the children.

**Lifecycle when rock is replaced.** New rock carries no inherited rift damage:

- A column reset by `lithosphere_plate._erupt_melted_nodes` (margin or column-pass eruption)
  resets all three fields to zero.
- A cell inserted at a margin by `quad_tectonics._open_rift` takes the area-weighted damage of its
  donors scaled by its `stretch_share`. Its magmatic share contributes zero, and `rift_broken` is
  false.
- An arc or other freshly accreted cell starts at zero.
- Removed rock (subduction, consumption in a merge) takes its damage with it, and nothing is
  redistributed.

Old saves load with all three at their zero defaults (`surface_fields.SURFACE_FIELDS`), so a
pre-revision world starts with no rift damage. Its belts, faults and plume heat still mark
inherited weakness.

### Each step, for every plate

The damage update (steps 1–5) runs once per step for **every** plate, in its own pass before the
split pass. It does not depend on split eligibility, plate order or whether another plate split
this step. `merge_split` splits at most one plate per step and skips the rest once one has split,
and age and cooldown gates exclude new daughters, so a damage update inside `maybe_split_plate`
would stop straining and healing on those plates. Only the breakthrough decision (step 6) is
gated by minimum age and cooldown. Daughters inherit their cells' damage fields through the
mask-based split.

Within the update, the loads, the column capacities `N_cap`, the plate's rift line and the
membership of its band (step 3) are frozen at the start of the step. The band's cells then
sub-step their strain, healing and softening together, and all advance synchronously, so the result does not depend on cell order. The line's
position and the coupling to thinning (through `N_cap`) are operator-split at the step, so
results depend on step length only within the tolerance D3 measures.

**1. Load, `N` (N/m, extensional, depth-integrated).** Each cell carries a symmetric 2×2
membrane-resultant tensor `𝐍` in its tangent plane. Tension is positive. The extensional load is
the largest principal value, `N = max(0, λ_max(𝐍))`, and its principal direction gives the band's
strike (perpendicular to the largest tension). `𝐍` is estimated by line cuts through the plate's
balanced load (proxy C, adopted after D0; see the Revision note):

- **Balanced load.** Sum the per-node forces the torque engine already computes
  (`torque.gather_per_node_force_field`: slab pull, slab drag, ridge push, basal drag and contact
  drag), project them into each cell's tangent plane, and remove their area-weighted
  rigid-rotation part, as `thin_shell.solve_thin_shell` does. What remains has zero net torque, so
  it is the part that stresses the plate rather than turning it. Slab pull and basal drag enter
  with their own signs; collision contact and ridge push push inward, so colliding margins come
  out compressive.
- **Line cuts.** In a frame tangent to the plate's centroid, for each of six directions `n̂`,
  sort cells along `n̂`. For the straight line through cell `x` normal to `n̂`, the balanced force
  on the far side, projected on `n̂`, divided by the plate's width there, is the resultant across
  that line, `N_n̂(x)` (N/m). The width is the cell area within half a line spacing of the line,
  divided by the spacing. A running sum over the sorted cells gives every line in `O(N log N)`.
- **Tensor.** Fit `n̂ᵀ𝐍n̂ = N_n̂` over the six directions by least squares, rotate the result into
  each cell's tangent basis, and multiply by `k_cut`.
- **Section pull.** The same running sums give the total pull across each straight section,
  `P_n̂(s) = k_cut × (balanced force on the far side) · n̂` (N), before division by the width. It
  is the force trying to separate the two sides of the plate along that section. Step 3 uses it
  to place the rift line.
- **Plumes:** while a cell has plume heat it adds an isotropic tension `N_plume 𝐈`, standing in
  for spreading above a plume head. The heat ends partway through a step at
  `t_exp = plume_heat_myr`. Averaging the load over the step before applying the nonlinear strain
  law would be timestep-dependent: a background load of 0.8 strength plus a plume load of 0.4
  strength for half the step averages to exactly the strength and gives no strain, though the
  heated half strains. The update therefore splits the step at `t_exp`, with full plume load and
  plume weakening before it and none after, and sub-steps each part (step 3).

D0 compared this estimate with the thin-shell resultant `N_ext` (`thin_shell.solve_thin_shell`),
which solves the same force field exactly but costs about 1.3 s per plate. On 451 development
plates it passed every D0 criterion with `k_cut = 1.33`, close to the physical value of 1,
at about 2 ms per plate. The shell solve is not on the step path; D1 runs it in shadow as a check.

**2. Strength, `N_eff` (N/m).** It is `N_eff = k_cap N_cap(kind, Hc, Hm, ΔT_moho) × M × W`, where:

- `N_cap` is the column's integrated extensional capacity to normal-fault yield, from the
  RFC 0010 Gate 2 yield envelope (`analysis/issue487/`, hydrostatic reference, `K0 = 0.8`). It
  is precomputed once on a grid of crust kind, Hc, Hm and Moho thermal lag, and interpolated.
  The 100 km continental anchor is 6.28 TN/m.
- `k_cap` is a single capacity scale. Gate 2 showed that the unconstrained background stress
  ratio `K0` moves the over-capacity share from 3% to 58% of cells, so the absolute scale cannot
  be derived and is calibrated instead.
- `M` is the product of the accepted inherited-weakness multipliers, clamped to `[0.12, 1]`: belt
  0.65, suture 0.55, volcano 0.65, plume heat 0.65 and fault 0.60. The plume multiplier applies
  only before `t_exp` in the split step. The `thinning` multiplier is dropped, because thinning
  now enters through Hc, Hm and `W`.
- `W = 1 - min(s_max, s_fault φ(ε̄) + s_warm φ(ε̄_warm))` is strain softening, with
  `φ(ε) = min(1, ε / ε_soft)`. `W` lies in `[1 - s_max, 1]`. `ε̄` and `ε̄_warm` are the two strain
  fields averaged over a disk of diameter `w_0` around the cell (area-weighted). This non-local
  average gives the band a physical width scale: softening cannot concentrate into a zone
  narrower than about `w_0`, whatever the node spacing. Local softening alone would make band
  width follow the grid. `s_fault` and `s_warm` are reductions: `s_fault = 0.40` matches the
  accepted fault multiplier 0.60. At the defaults `s_fault + s_warm = 0.80 < s_max = 0.88`, so the
  cap does not bind; it guards the D1 search once `s_warm > 0.48`.

**3. Rift line and strain.** Break-up is decided for a whole line across the plate. Testing each
cell against its own strength strains the weakest rock wherever it is (#894), so cells strain
only on the plate's rift line.

- **Pull axis.** For each of the six directions, take the section with the largest `P_n̂(s)`
  among sections that leave at least `a_min` of the plate's area on each side. Directions whose
  peak is positive and at least `c_dir` of the largest peak are candidates. A plate with no
  positive peak has no line this step.
- **Corridor.** For a candidate direction, the corridor is the contiguous run of sections around
  the peak whose pull is at least `c_corr` of the peak, again leaving `a_min` on each side. It is
  a strip across the plate, normal to the pull.
- **Line.** One Dijkstra search (`scipy.sparse.csgraph`) inside the corridor finds the cheapest
  edge-connected chain from the plate-boundary cells at one end of the strip to those at the
  other. The ends are the corridor's plate-boundary cells in the lowest and highest tenth of
  the coordinate along the strip. Intact cells cost `N_eff × ds` and broken cells cost zero, so
  the line follows belts, old rift strain, plume-heated and thin rock, and passes around
  cratons. A single shortest path cannot revisit a cell. A direction whose corridor is not
  connected from end to end is dropped.
- **Section load.** With `n̂ᵢ` the in-plane unit normal to the line at cell `i` and `ds_i` the
  line's length through the cell (half the centre-to-centre distances to its two line
  neighbours, so adaptive resolution does not bias the sum), the pull across the line uses the
  path-normal component `N⊥,i = max(0, n̂ᵢᵀ 𝐍ᵢ n̂ᵢ)`, not the largest principal value, which may
  point along the line:

  ```text
  P = Σ_line N⊥,i ds_i        S = Σ_intact N_eff,i ds_i        (both in N)
  ```

  The section load `Λ` is always finite, with an explicit branch for a line that has no intact
  cell left:

  - `P ≤ 0`: `Λ = 0`. The line neither grows nor splits, whether or not it is fully broken.
  - `P > 0` and `S > 0`: `Λ = P / S`.
  - `P > 0` and `S = 0`: the line is **open**. No ratio is formed. An open line ranks above
    every line that is not open (open lines rank among themselves by `P`), its growth ratio
    takes its cap (below), and step 6 splits it as soon as age and cooldown allow.

  The candidate with the largest `Λ` is the plate's rift line for the step. No line is stored.
  Strained cells are softer and broken cells cost nothing, so a line that already holds a rift
  has a larger `Λ` and is found again. The memory is in the rock.
- **Nucleus and tips.** Let `d_i` be the distance along the line from cell `i` to the nearest
  line cell with `rift_strain ≥ ε_cand`. If the line has none, the nucleus is the line cell with
  the largest `N⊥ / N_eff`, typically a plume-heated cell or a belt, and `d_i` is measured from
  it.
- **Rate.** Line cells strain when the section load passes the growth threshold:

  ```text
  r = min(Λ / Λ_grow, r_max)        (r = r_max on an open line)
  dε_i/dt = A · exp(-d_i / ℓ_tip),    A = ε̇_0 · max(0, r - 1)^n        [1/yr]
  ```

  The cap `r_max` bounds the rate at `ε̇_0 (r_max - 1)^n`, so no state gives an infinite or
  undefined rate. Cells inside the rift strain at the full rate `A`. Cells ahead of a tip
  strain more slowly until they reach `ε_cand` and join it.
- **Tip speed on the mesh.** Distances along the line come in steps of the cell spacing `h`, so
  the tip moves one cell at a time. At constant `A`, with `x = h / ℓ_tip`:
  - **Start-up.** The nucleus reaches `ε_cand` after `ε_cand / A`. An undamaged cell one
    spacing ahead strains at `A e^{-x}` and reaches `ε_cand` after `ε_cand e^{x} / A`.
  - **Steady advance.** A cell has strained at `A e^{-kx}` while it was `k` cells ahead, so it
    arrives one cell ahead with `A T / (e^{x} - 1)` already accumulated, where `T` is the time
    per advance. Setting that to `ε_cand` gives `T = ε_cand (e^{x} - 1) / A` and

    ```text
    v = h / T = (A ℓ_tip / ε_cand) · x / (e^{x} - 1)
    ```

    The first factor is the continuum speed; the second, the mesh factor, is 0.58 at `x = 1`
    and 0.77 at `x = 0.5`, and tends to 1 as the mesh is refined.
  - **Resolution.** Halving the node density multiplies `h` by `√2`. With `ℓ_tip = h` that
    changes the mesh factor from 0.58 to 0.45, a 22% slower tip. With `ℓ_tip = 2h` it changes
    from 0.77 to 0.69, an 11% slower tip. The starting value is therefore `ℓ_tip = 125 km`.
  - D1's propagation check measures start-up and steady speed numerically at both densities
    before the search (Diagnosis gates).
- **Width.** Cells within `w_i / 2` of line cell `i` strain at its rate (the largest rate among
  the line cells in reach), with `w_i = min(3 w_0, w_0 (1 + ε_i))`. The band starts `w_0` wide
  and widens as it stretches, up to the super-sea limit in the realism targets.
- **Failure.** When `Λ` falls below `Λ_grow`, because the pull dropped, turned, or the warm
  part healed, nothing strains. The rift is a scar that keeps its `rift_strain`.

- **Band.** The band is the line's cells plus every plate cell within `w_i / 2` of line cell
  `i`, fixed at the start of the step. A band cell with no eligible flank recipient anywhere on
  the plate (step 4) is held: it does not strain this step.

Strain increments `rift_strain` and `rift_warm_strain` equally; `rift_warm_strain` then decays as
`d(warm)/dt = -warm / τ_heal`. `rift_strain` is the extensional strain actually applied to the
column. It increases only while at least one of the crust and mantle lithosphere is above its
floor (step 4), and each cell stops at its latch strain, the strain at which it breaks: the
smaller of `ε_break` and `ln(Hc / RIFT_CRITICAL_THICKNESS_M)` at the start of the step. A cell
can therefore become broken, soft or cheap to cross only by strain that thinned it.

As the rift strains, softening lowers `N_eff` on the line, and broken cells leave the
denominator, so `Λ` and the rate rise. This is the slow-to-fast transition. The rate therefore
depends nonlinearly on the strains, so each part of the step is integrated in explicit
sub-steps:

- Each sub-step lasts `δt = min(t_left, Δε_max / max_i rate_i)`, so no band cell moves by more
  than `Δε_max`. A cell's increment is `min(rate_i δt, latch_i - ε_i)`: it reaches its latch
  strain only when its own rate over the elapsed time carries it there.
- `Λ`, `d_i` and `W` are recomputed from the current strains at the start of every sub-step,
  and the decay uses its exact exponential.
- Sub-steps continue until the step's time is used. No count limit converts remaining time into
  strain.
- The count is bounded: every rate is at most `ε̇_0 (r_max - 1)^n`, so each sub-step except the
  last lasts at least `Δε_max / (ε̇_0 (r_max - 1)^n)`. At the starting values that is at most
  15 sub-steps per 1 Myr step, and at most 54 at the top of D1's `ε̇_0` range. A part that would
  exceed `K_fail = 4096` sub-steps stops the run with an error, which can happen only through
  a defect. D3 reports the sub-step counts.

D1's integration check compares long steps with a converged fine-step reference (Diagnosis
gates).

Scale check, at constant `Λ = 2 Λ_grow` (so `A = ε̇_0 = 8e-8 /yr`), `ℓ_tip = 125 km` and the
default spacing `h = 62.5 km`:

- A rift cell strains at Brune et al.'s slow phase (under 10 mm/yr full rate) over a 125 km
  band and reaches `ε_break` in 25 Myr without softening.
- The nucleus reaches `ε_cand` after 1.25 Myr, and the first cell ahead after 2.1 Myr. Pre-strain
  then builds up, and within a few cells the advance settles at `T = 0.81 Myr` per cell, or
  77 km/Myr.
- Two tips cross a 5,000 km line in about 34 Myr at this constant rate, sooner as `Λ` rises.

These are within the tens of Myr that Earth's rifts take [[Brune et al., 2016]](#ref-brune2016).
They are first guesses at constant load, not a validated propagation law. D1 measures the
propagation and tunes the knobs.

**4. Thinning, subsidence and relief.** Pure-shear thinning over the step scales the cell's Hc
and Hm by `exp(-Δε)` [[McKenzie, 1978]](#ref-mckenzie1978), so the band subsides isostatically
into a rift valley. The model's cells cannot gain area inside a plate, so the removed volume
moves to the band's flanks. The flank-transfer helper adapts `_allocate_stretch`'s volume and
tracer rule in reverse. At margins, `_open_rift` thins pre-existing donors and gives their volume
to newly inserted cells. Here the straining cells of the band are the donors, and the recipients are cells
of the **same plate** within `RIFT_FLANK_ROWS` hops that did not strain this step, in proportion
to their area. Straining cells never receive, so volume is not shuffled within the band.

- **Volume moved.** The transferred volume is computed from the donor's thickness change after
  clamping to the floors (below), so a floored donor gives up only what it actually lost.
- **Tracers.** Hc and Hm move with continental material, mobile cover and its continental share,
  and the plume-head and arc tracers, as `_allocate_stretch` moves them. Craton crust follows
  `_open_rift`'s existing rule: it thins in place on the donor, is not moved to the flanks, and the
  lost volume is booked to the `rifted_stretch_m3` craton sink, which counts toward D3's
  craton-loss limit.
- **Recipients.** Recipients are unbroken plate cells outside the step's band. If a donor has
  none within `RIFT_FLANK_ROWS`, the reach widens one hop at a time with no upper limit (one
  multi-source search per step). The band is at most `3 w_0` wide, so a recipient is normally
  within three hops. A band cell with no recipient anywhere on the plate (every other cell
  broken or in the band) is held for the step (step 3): it does not strain, so its
  `rift_strain`, softening and latch do not change and no volume moves. Held cells and donors
  that needed more than four hops are counted for D1.

Donor and recipient elevations are then recomputed isostatically. Volumes are conserved exactly,
and the flanks rise slightly, like rift shoulders.

Interior band cells do not erupt before breakup. Thinning stops at the floor
`RIFT_CRITICAL_THICKNESS_M` for crust and `MIN_MANTLE_LITHOSPHERE_THICKNESS_M` for mantle
lithosphere. A cell is latched broken (`rift_broken = true`) when its crust reaches the floor or
its `rift_strain` reaches `ε_break`, whichever comes first. Because `rift_strain` counts only
applied strain (step 3), both latches record thinning that happened. For 35 km crust the floor comes first
(at strain 1.95); for crust thicker than about 37 km, `ε_break` does. `lithosphere_plate._erupt_melted_nodes` resets a
melting land cell to the 35 km continental reference column, which would re-thicken the band, so
it is not called here. New oceanic crust appears only after breakup, when the daughters separate
and the existing margin path (`_open_rift`) fills the gap. Normal-fault
relief comes from the existing fault machinery: the spawn weight in
`faults._maybe_spawn_faults` gains a term proportional to the local strain rate, so active
bands grow sub-parallel normal faults, striking along the band's principal direction.

**5. Physical candidates (diagnostic).** A physical candidate is a connected component of cells
with `rift_strain ≥ ε_cand`. The damage fields cannot recover a candidate's history: onset
time, onset width, time since last strain or rate history. A diagnostic recorder in the gate
scripts therefore tracks candidates, outside saved state:

- **Sampling:** at every step, record each component's cells, area, length (its longest geodesic
  extent along the principal strike), width (area divided by length), mean strain rate and the
  share of its length on inherited weakness.
- **Identity:** a component keeps the identity of the previous step's component it overlaps most
  (by area). When components merge, the merged one keeps the oldest identity and the others are
  recorded as merged. When one splits, the largest piece keeps the identity and the rest start
  new ones with the parent's onset time.
- **Onset** is the first sample of an identity. **Onset width** is its width there. A candidate
  is *active* while any of its cells strain, *failed* once none has strained for `t_fail`, and
  *successful* when its path splits the plate.
- **Censoring:** candidates still active at the end of a run are censored. They are counted, but
  excluded from duration and failure-share estimates, which use survival estimates over censored
  and completed candidates. Candidates already present when a run starts (from a loaded save)
  have unknown onset and are excluded from onset statistics.

Production code needs only the three damage fields; the recorder runs in the D1–D3 harness.

**6. Breakthrough and split.** The rift line from step 3 is already an ordered chain of cells,
each sharing an edge with the next, whose two ends are plate-boundary cells. The plate splits
when the intact part of the line, the ligament, can no longer carry the pull across it:

```text
P ≥ S        (the line is open, or Λ ≥ 1, with P > 0)
```

and at least `f_break` of the line's length is broken (step 4), so a line never splits before a
rift has grown along it. Another rift's broken cells also cost zero, so a line can link rifts.
An open line splits as soon as age and cooldown allow. While cooldown blocks it, its band keeps
straining at the capped rate, which is finite and stops at each cell's latch. A line with
`P ≤ 0` never splits, whether or not it is fully broken: tension parallel to the line never
splits it, and a broken line whose pull has stopped stays an unsplit, failed rift.

Remove the line's cells from the plate's adjacency graph. The two largest remaining components
are the daughters. Any other component is a stray piece: if it holds at most `a_stray` of the
parent's area it joins the daughter it shares the most edge length with, and otherwise the split
is not applied. Each line cell then joins the daughter with which it shares the most edge
length, with ties going to the side whose fitted rotation it matches best. The split is applied
only if both daughters are connected after that assignment, the smaller holds at least `a_min`
of the parent's area, and the ownership checks pass. Otherwise the rift keeps growing. The rule
is deterministic. Randomness enters through the model's sampled plume events and stochastic
faults, consistent with the owner's rule to threshold endogenous causes and use probability
only for exogenous ones.

### Starting values

| Knob | Start | Source | Tuning |
|---|---:|---|---|
| `ε̇_0` | `8e-8 /yr` | Brune et al. 2016's slow phase (< 10 mm/yr) over a 125 km band | D1 search over 2e-8 to 3e-7 /yr |
| `n` | `1` | Design choice | Fixed; sensitivity at 2 |
| `ε_soft` | `1.0` | Design choice; linear softening form of Huismans & Beaumont 2003 | Fixed; sensitivity at 0.5 and 2 |
| `Λ_grow` | `0.3` | Design choice. Intact interior sections carry a median 0.10–0.40 of their strength on five development worlds (#894) | D1 search over 0.1–0.6 |
| `ℓ_tip` | `125 km` | Design choice; twice the default spacing keeps the mesh factor of the tip speed within 11% between default and half density (step 3) | D1 search over 60–250 km; D1 propagation and resolution checks |
| `r_max` | `10` | Design choice; bounds the rate at 9 `ε̇_0` (with `n = 1`), above the 3× slow-to-fast target | Fixed; sensitivity at 5 and 20 |
| `a_min` | `0.15` | Owner decision; above the accepted 7% tiny-breakup limit | Fixed |
| `a_stray` | `0.01` | Design choice; #894's stray pieces were 1–27 cells | Fixed; D1 reports splits it cancels |
| `c_dir`, `c_corr` | `0.8`, `0.8` | Design choices | Fixed; D1 line check reports sensitivity at 0.6 and 0.9 |
| `f_break` | `0.5` | Design choice | Fixed; sensitivity at 0.25 and 0.75 |
| `s_fault`, `s_warm`, `s_max` | `0.40`, `0.40`, `0.88` | `s_fault` matches the accepted fault reduction; `s_max` matches the accepted 0.12 clamp | D1 search on `s_warm` |
| `τ_heal` | `62.8 Myr` | McKenzie 1978 thermal time constant | Fixed; sensitivity at 30 and 120 Myr |
| `k_cap` | `1.0` | #487 reference envelope | D1 search over 0.3–3 (the `K0` bracket). Only the ratio `k_cap Λ_grow` sets growth onset; `k_cap` alone sets break-up |
| `k_cut` | `1.33` | D0 fit to thin-shell `N_ext` on 451 development plates | Fixed after D0 |
| `N_plume` | `0.5 TN/m` | Design choice, under a tenth of the 100 km continental capacity (6.28 TN/m) | D1 search over 0–2 TN/m |
| `w_0` | `125 km` | Owner's onset width (100–150 km) | D1 resolution check |
| `ε_cand` | `0.1` | Design choice | Diagnostic only |
| `ε_break` | `2.0` | Design choice; about `β = 7.4`, which takes 35 km crust below 5 km | Fixed |
| `t_fail` | `5 Myr` | Design choice | Diagnostic only |
| `RIFT_FLANK_ROWS` | `2` | Matches `K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION` | Fixed |
| `K_fail` | `4096` | Defect guard, far above the bound of 54 sub-steps per 1 Myr step (step 3) | Fixed |
| `Δε_max` | `0.05` | Design choice | D3 timestep check |
| Inherited multipliers | accepted RFC 0018 values | Accepted design choices | Fixed |

At the default node density the cell spacing is 62.5 km, so a band at the `w_0 = 125 km` scale
is about two cells wide. That is the owner's chosen onset width
([#876 decisions](https://github.com/adubey/mantle-bloom/issues/876#issuecomment-6874579897)),
set by the strip around the line and the non-local softening average, not by the grid. The band
widens with its strain (step 3). D1 checks the width at two node densities.

### Runtime and saves

Production runs, per plate per step: the line-cut load tensor from existing per-node forces
(about 2 ms per plate in D0, six sorts plus a 3×3 rigid-rotation fit); a table interpolation for
strength; and up to six corridor path searches, usually one to three, each on a strip of the
plate. The disk average and the strain update touch only the line's band. The connected-component
check runs only when a line passes the breakthrough test. The candidate
recorder runs only in the D1–D3 harness. The hypothesis is milliseconds per plate, small against
a production step, though under plates-only stepping it could be a noticeable share. This is not
yet measured; D3 measures it on production settings. The thin-shell solve is not
on the step path. Saves gain two float surface fields and one boolean, and old saves load with zeros.

### Diagnosis gates

Development worlds are the 20 RFC 0018 development worlds
([`analysis/issue797/`](../analysis/issue797/README.md)) and the #848 development seeds. The
held-out sample is the ten-seed revised-method sample already set aside: registry seeds 12–20
plus the first frozen reserve seed `1479220500`, none of which has been run for a revised method.
No held-out seed is used for calibration or tuning.

1. **D0: calibrate the cheap load.** On development worlds, fit the load proxy's free parameters
   (originally `L_drag` and `L_pull`; for the adopted proxy C, `k_cut`) so the cheap tensor matches the thin-shell resultant under the same forces. First check
   manufactured cases with known answers: one-sided slab pull, two slabs pulling apart, a slab
   opposed by collision contact, and pure collision (compression only). The proxy must give the
   shell's sign (tension or compression) in each. A further case puts tension parallel to an
   intact path: the ligament rule must not split it. Then, on development worlds over cells where
   the shell is loaded, require:
   - **sign agreement** (extension versus none) on at least 80% of cells;
   - a **principal-direction difference** under 30° on at least 70% of extensional cells;
   - **Spearman correlation** of `N` with `N_ext` of at least 0.5;
   - a **median ratio** `N / N_ext` between 0.5 and 2 in the top decile, by crust kind.

   If any of these fails, use alternative B: run the thin-shell solve every `K` steps
   (start `K = 5 Myr`), hold `𝐍` between solves, and measure its cost in D3.

   **Result (2026-11-14):** the originally specified local proxy (signed drag differences with
   coupling length `L_drag`, plus slab, ridge and contact forces spread over a decay length
   `L_pull`) failed. Sign agreement was 74.6%, direction 52.9% and Spearman 0.49; the magnitude
   ratio passed; 3 of 4 manufactured cases matched. The line-cut proxy C passed every criterion:
   sign 82.4%, direction 71.9%, Spearman 0.69, ratios 0.57 continental and 0.54 oceanic, and 4 of 4
   manufactured cases, with `k_cut = 1.33`. The owner adopted C, so alternative B is not used
   ([`analysis/issue797/`](../analysis/issue797/README.md), [PR #887](https://github.com/adubey/mantle-bloom/pull/887)).
2. **D1: check the lines, then tune the mechanism on development worlds.**
   - **Line check (before any stepping search).** On the 20 development worlds at 0 Myr, and
     after 50 Myr of plates-only stepping with the production splitter at 1 Myr steps, compute
     every eligible plate's rift line (step 3) without straining anything. Report the smaller
     side's area share, the line's length over the straight distance between its ends, the
     craton and oceanic shares of its length, `Λ`, and whether cutting along it would give a
     valid two-daughter split. Repeat the 50 Myr worlds from 0.5 Myr steps on three seeds. Pass
     when at least 90% of lines give a valid split and the line crosses less craton than the straight chord between its ends on at least
     two thirds of plates that have cratons in the corridor. If this fails, the line placement
     is revised before any search is run.
   - **Integration check (before the search).** On a manufactured plate with a fixed force
     field and a fixed line, integrate 60 Myr with steps of 0.01 Myr (the reference), 0.5, 1, 2
     and 5 Myr. Include a case whose rate sits at the `r_max` cap, which needs the most
     sub-steps. Pass when, at every common output time, accumulated `rift_strain` on the band
     is within 0.05 of the reference, the set of broken cells differs from the reference by at
     most one cell at each tip, and the split time is within one coarse step of the reference.
     A failure means the integration or the operator split is revised before any search.
   - **Propagation check (before the search).** On a manufactured straight line at constant
     `Λ`, measure the nucleation time, the first advance and the steady tip speed at the default
     and half node density, for `ℓ_tip` at 60, 125 and 250 km. Report them against step 3's
     discrete law. If the measured steady speed differs from the law by more than 20%, the scale
     check and `ℓ_tip`'s range are corrected before the stage 1 predeclaration.
   - **Search.** Run 200 Myr with plates-only stepping (`simulate_climate_biomes = False`) at
     1 Myr steps, as G2 did, across a predeclared parameter search: a Latin hypercube over
     `ε̇_0`, `s_warm`, `k_cap`, `N_plume`, `Λ_grow` and `ℓ_tip`, as in #313. Score each setting
     against the realism targets below. The search runs in two stages (owner decision,
     2026-11-14):
     - **Stage 1:** a coarse search on 5 of the 20 development worlds, with the settings and the
       five seeds predeclared before the first run.
     - **Stage 2:** the best few stage-1 settings on all 20 development worlds. Freeze the best
       before D2.
   - **Timestep.** Repeat the frozen setting at 0.5 Myr steps on three development seeds and
     report the realism targets. D3's half-step comparison is the binding timestep check.
   - **Fragmentation.** Record the #699 failure mode for every run: the largest damaged share of
     a plate and the number of damaged components at the time of each split.
   - **Oceanic rifts.** Record successful rifts by crust kind: a rift is *oceanic* when more than
     half its line length runs through oceanic crust. If more than one third of successful rifts
     in the frozen setting are oceanic, apply the predeclared fallback: oceanic cells latch when
     both crust and mantle lithosphere reach their floors or at `ε_break` (strain still counts
     only while a layer thins), and `k_cap` gains a separate oceanic factor searched over 1–5.
     Then repeat D1.
   - **Also reported:** held band cells and far flank transfers (step 4), splits cancelled by
     the stray-piece limit, the share of rift strain left off the final line, and sub-step
     counts.
   - **Shadow shell.** On the stage-2 runs, solve the thin shell every 10 Myr in shadow and
     report D0's statistics for C on the evolving worlds. This is reporting, not a pass
     condition.
   - **Resolution.** Repeat the frozen setting on five development seeds at half the default
     node density. Onset width and median onset-to-breakup duration must stay within 25% of the
     default-density values, or band width or tip speed is following the grid.
3. **D2: realism on held-out worlds.** With the frozen knobs, compare the rift-line method with
   the production great-circle method on the ten held-out seeds for 200 Myr at 1 Myr steps, from identical
   generated worlds, settings and step schedules. The method passes when it meets every realism
   target, including the split-frequency target below, which is judged on the rift-line arm
   alone:
   - **Rift-arm count.** The rift-line arm makes at least 30 actual splits across the sample.
     If it makes fewer, add seed pairs in order from the frozen reserve list, excluding
     `1479220500`, up to 10 more. If it still has fewer than 30, D2 fails. Adding seeds cannot
     rescue a method that splits rarely, because the rate test below is computed over every
     pair run.
   - **Rift-arm rate.** Over all `n` pairs, the rift-line arm's pooled split rate (splits per
     world per 100 Myr) is at least half the production arm's, and at least 70% of worlds have
     at least one rift-line split. The factor of one half is a design choice. The production
     splitter keeps plate counts in range today, and D3's 10–40 plate limit is the independent
     check. Report the rate ratio with a 95% paired seed-bootstrap interval (10,000 resamples of
     whole seed pairs); the pass uses the point estimate.
   - **Sample sizes.** Each realism target states the events it is computed over. A target
     computed over successful rifts needs at least 20 of them, and one computed over physical
     candidates needs at least 20 with an observed onset; the failure share uses the survival
     estimate over at least 20 candidates, censored ones included. A target below its minimum is
     *not evaluable* and counts as not met. Report each target's `n` with its seed-bootstrap
     interval.
   - **Tiny breakups (paired comparison, kept).** No more than `max(1, floor(0.1 n))` pairs may
     have more tiny breakups than their baseline (`D_s = T_baseline,s - T_rift,s < 0`; smaller
     daughter under 7% of parent area). With `a_min = 15%` the rift-line arm cannot make a tiny
     breakup, so this is met by construction. It is reported for comparison with G2.
4. **D3: runtime and timestep.** Keep the accepted G3 limits: median runtime within 2% of the
   paired production baseline on three held-out seeds, from controlled serial runs with profiling
   kept separate; craton loss at most 2% per Gyr; and 10–40 plates at every output. Repeat at half
   the step size on three predeclared held-out seeds. Plate count and craton loss must stay
   within 20% of the normal-step results, as accepted. So must the tiny-breakup share; when the
   normal-step share is zero, the half-step run must also have zero tiny breakups. Breakup count
   and median rift duration are added, also within 20%.

**Realism targets** (median over physical candidates unless stated):

| Target | Pass range | Basis |
|---|---|---|
| Onset width (first time `ε ≥ ε_cand`) | 100–150 km | Owner design choice |
| Width growth | Successful bands are wider at breakup than at onset | Owner design choice |
| Onset-to-breakup duration | 10–50 Myr | Tens of Myr [[Brune et al., 2016]](#ref-brune2016); range is a design choice |
| Slow-to-fast transition | Strain rate in the last quarter before breakup at least 3 times the first quarter | Abrupt acceleration [[Brune et al., 2016]](#ref-brune2016); factor is a design choice |
| Breakups along inherited weakness | At least two thirds of successful paths have most of their length on belts, sutures, faults or earlier rift strain, where such marks exist on the plate | 17 of 18 margins [[Buiter & Torsvik, 2014]](#ref-buiter2014); threshold lowered for the model's incomplete suture signal |
| Failed rifts | Between 20% and 80% of physical candidates fail | Failure is common [[Brune et al., 2023]](#ref-brune2023) [[Elling et al., 2022]](#ref-elling2022); range is a design choice |
| No super seas | No below-sea-level basin made by interior rift thinning is wider than three onset band widths (about 375 km) across its narrow axis, outside crust thicker and hotter than the plate median | Narrow rifts dominate; wide rifts need thick, hot crust [[Buck, 1991]](#ref-buck1991); #866 |
| Localization | At each split, damaged cells form at most three components above 1% of plate area, and the damaged share stays under 30% of plate area | #699's fragmentation failure mode; thresholds are design choices |
| No microplates | Every split leaves both daughters with at least `a_min` of the parent's area. The distribution of the smaller share is reported; no ratio is targeted | Owner's aim: a split makes two plate-sized daughters, not a microplate at the edge (#894) |
| Line shape | Median line length at least 1.1 times the straight distance between its ends | Owner's aim that rifts are not straight lines (#894); threshold is a design choice |
| Plume link | Reported: share of breakups with plume heat on the path | 12 of 18 margins [[Buiter & Torsvik, 2014]](#ref-buiter2014); descriptive only |
| Splits happen | The rift-line arm alone makes at least 30 actual splits, its pooled split rate is at least half the production splitter's on the same worlds, and at least 70% of worlds have a rift-line split. In D1 the production rate comes from one production run per development world at the same settings and steps | Plates must keep splitting, not only thin (PR #897 review); the factor and counts are design choices |

**Pictures.** For five held-out seeds, show map sequences of one successful and one failed rift,
from onset to breakup or to failure: land, water, cratons, belts, rift strain and plate outlines,
in a shared projection and legend.

### Acceptance criteria

- [x] D0 records the fitted load proxy with its agreement statistics: C passed with
      `k_cut = 1.33` and was adopted (2026-11-14).
- [ ] D1 records the line check, the parameter search, the scores and the frozen knob values,
      using development seeds only.
- [ ] D2 passes on the ten held-out seeds (with the reserve-list fallback): every realism target
      is met and evaluable, the rift-line arm alone makes at least 30 splits at a pooled rate of
      at least half the production arm's with a split on at least 70% of worlds, and no more than
      `max(1, floor(0.1 n))` pairs have `D_s < 0`.
- [ ] D3 passes: median runtime within 2% of baseline, craton loss at most 2% per Gyr, 10–40
      plates at every output, and half-step plate count, craton loss, tiny-breakup share (with
      the zero-share rule), breakup count and rift duration within 20%.
- [ ] Every applied split makes exactly two connected plates, assigns each parent cell to one
      daughter, and preserves crust and mantle ledgers. Interior thinning conserves crust,
      continental material and tracer volumes.
- [ ] Old saves load with zero rift strain and continue stepping (`test_persistence.py`).
- [ ] Unit tests cover the bounds of the strength factors, softening monotonicity, timestep-robust
      strain integration, flank volume and tracer conservation, healing, the broken latch (no
      interior eruption), remap without healing (coarsen a quartet with one broken cell; it stays
      broken and its strain does not fall), damage reset on eruption and fresh crust, margin
      insertion inheriting only its stretched share, plume expiry (a crossing step matches the
      step split at expiry), flank transfer with floored donors and with recipients beyond
      `RIFT_FLANK_ROWS` (volume and thickness bounds hold), a band cell with no recipient held
      over enough steps to pass its break strain (its `rift_strain`, softening and latch never
      change and the plate never splits), `rift_strain` never exceeding the strain applied to
      the column, the sub-step bound and the `K_fail` error, the damage update on a plate after another plate splits and on a
      daughter in cooldown, path-normal loading (tension along the line never splits), directional
      load signs (compression never
      registers as extension), the section pull of manufactured force fields, line placement
      (around a craton, along a belt, and inside the `a_min` limits), the same line found again
      once it holds a rift, tip start-up and steady speed against the discrete law under a
      constant `Λ`, no strain off the band, the integrated ligament rule, the section load's
      branches (an open line with positive pull splits when the gates allow; an open or partly
      broken line with zero pull neither grows nor splits; an open line in cooldown keeps finite
      rates and strains, then splits when cooldown ends), no infinite or undefined value in any
      branch, the `f_break` rule, stray-piece attachment and its limit, the `a_min` daughter
      minimum, line validity and two-daughter splits.
- [ ] The pictures above are in the gate report.
- [ ] `docs/simulation-model/tectonics.md` documents rift damage, bands and breakup, and
      `docs/architecture.md` lists any new module.

### Falsification criteria

- **Whole-plate pull does not place good lines.** If D1's line check fails (lines are invalid or
  cross cratons as often as straight chords), the premise that the section of
  largest pull marks where a plate should break is wrong for this model's forces. Revise the
  placement or add an interior load, such as plume uplift, before any search.
- **Lines do not persist.** If rifts are abandoned because the chosen line jumps between steps
  (most rift strain ends up off the final line in most runs), memory in the rock is not enough
  to hold a line. Store the line, or reject the method.
- **Damage still spreads.** If the localization target fails in most D1 runs, confining strain
  to the line's band has not fixed #699's failure mode. Diagnose before D2.
- **Tips do not propagate as a front.** If D1's propagation check finds tip speed set mainly by
  the node spacing (the half-density speed differs by more than 25% at every tested `ℓ_tip`),
  the tip law is a grid artifact. Replace it with a continuous tip position before the search.
- **The cheap load stops standing in for stress.** C passed D0 on development snapshots. If D1's
  shadow shell solves on evolving worlds fall below D0's thresholds, revise the load model before
  D2; alternative B (periodic shell solves) remains the fallback.
- **No knob setting looks right.** If no D1 setting meets most realism targets on development
  worlds, the mechanism is missing a process. Diagnose which target fails before adding knobs.
- **Rifts ignore weakness.** If successful paths follow inherited weakness no more often than
  matched random edge-to-edge paths, the strength multipliers are not doing their job. Revisit
  them.
- **Timestep dependence.** If D3's half-step comparison fails, make the strain integration
  stable or reject the design.
- **Craton safety or bookkeeping fails.** If craton loss exceeds 2% per Gyr or a ledger does not
  close, stop and diagnose before changing a bound.

## Detailed design (provisional touchpoints)

- `backend/app/rift_damage.py` (new): load proxy and section pull, strength table, rift-line
  placement, strain update, healing, flank transfer, physical-candidate diagnostics and the
  breakthrough rule.
- `backend/app/surface_fields.py`: add `rift_strain`, `rift_warm_strain` (float, zero default)
  and `rift_broken` (`BOOLEAN_PROVENANCE`) to `SURFACE_FIELDS`.
- `backend/app/sparse_quad_patch.py` coarsening and `backend/app/quad_merge.py`: add
  `rift_strain` to the max-preserving special cases beside `channel_depth`; clamp
  `rift_warm_strain` to `rift_strain` after remap.
- `backend/app/lithosphere_plate.py` (`_erupt_melted_nodes`) and `backend/app/quad_tectonics.py`
  (`_open_rift`, cell insertion): reset or inherit the damage fields as in the lifecycle rules.
- `backend/app/merge_split.py` (or the step pipeline that calls it): run the damage update for
  every plate in a pass of its own before the split pass.
- `backend/app/merge_split.py`: `maybe_split_plate` keeps its age and cooldown gates and calls
  `rift_damage` for the breakthrough decision and path (step 6 only), replacing the great-circle cut, its flow
  gates and the 0.55 success roll. Remove its `apply_failed_rift` call: failed rifts are bands
  that stop growing.
- `backend/app/plates.py` and `backend/app/sparse_quad_patch.py`: the mask-based split helper
  around `_plates_from_node_masks()`, as in the accepted design, with the `a_min` daughter
  minimum and stray-piece attachment. The great-circle `plate.split` and its per-daughter
  `SPLIT_MIN_NODES` check are no longer on the split path.
- `backend/app/quad_tectonics.py`: reuse the `_allocate_stretch` volume rule for flank transfer.
  Interior band cells do not call `lithosphere_plate._erupt_melted_nodes` (see step 4).
- `backend/app/torque.py`: reuse `gather_per_node_force_field` for the balanced load; no change
  to plate-motion forces. The line-cut estimate starts from `bin/debug/rift_damage_prototype.cut_tensor`.
- `backend/app/faults.py`: add a strain-rate term to the spawn weight in `_maybe_spawn_faults`.
- `backend/app/cratons.py`: change `record_rift_cut(world, plate, points, cut_normal, band_sin,
  kind)` to take the path's band mask instead of a cut normal and band width, so no great-circle
  geometry remains at its call site in `maybe_split_plate`. Belt marking
  (`mark_belt`) already accepts a mask.
- Strength table: precomputed from the Gate 2 envelope (`bin/debug/diagnose_gate2_column_envelope.py`)
  as a small checked-in data file, or computed at import; decide in implementation.
- `bin/debug/` (new scripts): D0 calibration against `thin_shell`, D1 line check and search
  driver, realism-target scorer and map sequences. `bin/debug/compare_rift_load_maps.py` is the
  starting point for D0. The #888 prototype (`rift_damage_mechanism.py`, `rift_damage_step.py`,
  PR #890) keeps its kernels for strength, softening, healing, flank transfer and the recorder;
  its per-cell strain test and band-first path search are replaced.
- Tests: `test_rift_damage.py` (new), `test_merge_split.py`, `test_sparse_quad_patch.py`,
  `test_cratons.py`, `test_faults.py`, `test_persistence.py`.
- Docs: `docs/simulation-model/tectonics.md`, `docs/architecture.md`.

## Assumptions and uncertainty

- **The load proxy is a calibrated stand-in.** Drag and pull are scalar estimates of the largest
  extensional resultant, not a stress tensor. Their orientation is implicit. D0 tests them
  against the thin-shell diagnostic, which itself assumes linear elasticity about a zero-load
  reference (#487).
- **Capacity scale.** `k_cap` absorbs the unconstrained background stress, pore pressure and
  rheology (#487 Gate 2), so the absolute threshold is a calibration, not a prediction.
- **Softening.** The linear softening form follows numerical rift models; its amount and strain
  scale are design choices tuned for appearance and timing.
- **Fixed-grid extension.** Cells inside a plate cannot gain area, so moving band volume to its
  flanks stands in for the band widening. It keeps volumes exact but places thickening on the
  shoulders. An alternative is to hold the volume until breakup and release it at the new margins.
- **Healing.** One thermal time constant stands in for lithospheric cooling and strengthening.
  The permanent part stands in for fault damage.
- **Plume load.** `N_plume` is a design choice; plume weakening already exists (RFC 0007).
- **Sutures.** Generated pseudo-sutures are still #770 work. Until a per-cell suture field lands,
  belts carry inherited weakness, so the inherited-weakness target is lowered.
- **Oceanic columns.** The model's reference oceanic crust (`REFERENCE_HC_OCEANIC_M`, 11 km)
  reaches the 5 km floor after a strain of about 0.79, against about 1.95 for 35 km continental
  crust. In the #487 reference envelope its capacity is also lower: mature oceanic 4.20 TN/m,
  against 15.82 TN/m continental (`analysis/issue487/gate2-column-envelope.json`). Oceanic
  interiors therefore strain more easily and latch earlier than continental ones. D1 records it and
  has a predeclared fallback.
- **Line-cut proxy.** C gives each straight line's average resultant, not the local stress: it
  smooths how stress turns near plate edges and corners, straight lines can cross non-convex plates
  more than once, and the centroid-tangent frame distorts very large plates. D0 found it good enough
  on development snapshots (direction within 30° on 71.9% of cells, just above the 70% bar), and
  D1's shadow shell solves check it on evolving worlds.
- **Width scale.** The non-local softening disk `w_0` sets band width directly, which is a
  regularization choice, not a prediction. D1 checks that width does not follow the node spacing.
- **One line per plate.** A plate grows one rift at a time, along straight-section pull in six
  directions. Triple junctions, branching rifts and curved pull on very large plates are not
  represented. In #894, 20–33% of plates had their largest pull at an edge; with `a_min` those
  plates use their best interior section, which may carry little pull and never rift.
- **Whole-section failure.** Using one `Λ` for the line assumes the plate transmits the pull
  along the whole section. It replaces a stress solve around a growing crack with a ratio and a
  tip decay length `ℓ_tip`, which is a regularization choice like `w_0`.
- **Small daughters.** With `a_min = 15%`, a plate at the `SPLIT_MIN_NODES` floor can shed a
  daughter of about 105 cells. D3's 10–40 plate limit and the merge rules are the check on
  small-plate build-up.
- **No uplift load.** The model has no buoyancy force from plume uplift. Plumes act through
  `N_plume` and weakening only.
- **Not geological rates.** Realism targets judge appearance and timing on model worlds. They are
  not calibrated rift-frequency estimates.

## Revision note

**2026-11-14 — rift lines replace cell-by-cell damage bands (revision after #894).** The
previous version strained every cell whose own load exceeded its own softened strength, and
split a plate when a broken band's cheapest path reached the boundary at both ends. The #888
prototype of it ([PR #890](https://github.com/adubey/mantle-bloom/pull/890)) made no split in
60 Myr on development seed 106533960. The exploration in
[#894](https://github.com/adubey/mantle-bloom/issues/894)
([`analysis/issue894/README.md`](../analysis/issue894/README.md)) found:

- **Load is flat; strength is low at edges.** Median load 0.50 TN/m at the edge against
  0.58 TN/m ten cells in; median strength 4.3 against 9.7 TN/m. The cells that pass a per-cell
  test are the weakest rock: margins and oceanic crust (66% of straining cell-steps).
- **Damage spreads.** Large bands were patches a median 3 times longer than wide, not oriented
  across the pull.
- **Path and split rules rejected the usable cuts.** 34% of path searches were dropped for
  overlapping legs. 42 ligament-passing paths that would have made two plate-sized daughters
  (between 67/33 and 80/20) were rejected for stray pieces of 1–27 cells and for the 700-cell daughter minimum.
- **Total pull peaks across plate interiors** in 67–80% of plates on five development worlds,
  and its direction is stable over 10 Myr in 47 of 56 cases.

Review of the revision ([PR #897](https://github.com/adubey/mantle-bloom/pull/897)) added: an
explicit, finite branch for a fully broken line; a capped growth ratio `r_max`; the discrete tip
law with `ℓ_tip = 125 km`; strain integration without a forced latch; strain counted only when
applied; and a D2 count and rate test on the rift-line arm alone.

The owner chose to revise to rift lines and to change two kept production rules: the daughter
minimum becomes 15% of the parent's area, and stray pieces of at most 1% of the parent's area
join the adjacent daughter. Step 3 now places one line per plate from the section pull and
strains only that line's band, with tips that grow along it. Step 6 splits along the same line.
The damage fields, load proxy C, strength, softening, healing, thinning, the broken latch, the
recorder and the path-normal ligament rule are kept. D0's result stands, because the section
pull comes from the same sums as proxy C.

| Gate | Previous criterion | Revised criterion | Reason |
|---|---|---|---|
| D1 | Two-stage Latin-hypercube search over `ε̇_0`, `n`, `ε_soft`, `s_warm`, `k_cap` and `N_plume` | An integration check (long steps against a 0.01 Myr reference, including a rate at the cap) and a propagation check (tip speed against the discrete law at two densities) on manufactured plates; a line check on the 20 development worlds first (valid splits on at least 90% of lines, less craton than the straight chord on two thirds of plates; smaller-side share reported). Then the same two-stage search over `ε̇_0`, `s_warm`, `k_cap`, `N_plume`, `Λ_grow` and `ℓ_tip`, with `n` and `ε_soft` fixed | The line placement is the new uncertain step and can be checked without stepping. Two new knobs replace two that only shaped the per-cell law, keeping six searched knobs for 15 settings |
| D1 half-density check | Onset width and duration within 25% | Unchanged; now also guards tip speed | `ℓ_tip` is a length scale that must not follow the grid |
| Integration | Sub-steps capped at `K_sub = 64`; a cell needing more was latched broken | No count cap; sub-steps are bounded by the capped rate (at most 54 per 1 Myr step in D1's range); a defect guard stops the run at 4096 | A count cap that forces a latch can create a rupture from the step length alone (PR #897 review) |
| Withheld thinning | A donor with no recipient within four hops kept accumulating `rift_strain` without thinning | Recipients are searched without a hop limit; a band cell with none anywhere is held and does not strain; `rift_strain` counts only applied strain | A cell could break without ever thinning (PR #897 review) |
| D1 oceanic fallback | Oceanic cells latch only on `ε_break` | Oceanic cells latch when both layers reach their floors or at `ε_break` | Strain past the crust floor must still thin a layer to count |
| Realism targets | Ten targets | Two added: no daughter under `a_min` of the parent (smaller share reported, no ratio targeted), and median line length at least 1.1 times the straight distance | The owner's aims for this revision: two plate-sized daughters, not a microplate, along a line that is not straight |
| D2 | At most one of ten pairs with more tiny breakups than baseline; 50 combined splits across both arms; every realism target | The rift-line arm alone makes at least 30 splits (reserve fallback up to 10 more pairs, else fail), its pooled split rate is at least half the production arm's, and at least 70% of worlds have a rift-line split; each realism target needs at least 20 events or counts as not met; every realism target; the paired tiny-breakup rule is kept but is met by construction with `a_min = 15%` | The combined count could be met by the baseline alone while the new method rarely split (PR #897 review); the daughter minimum changed |
| D3 | Runtime, craton loss, plate count and half-step limits | Unchanged | |

| Item | Verdict | Reason |
|---|---|---|
| D0 | Stands | Proxy C and `k_cut` are unchanged |
| [#888](https://github.com/adubey/mantle-bloom/issues/888) D1 | Re-scope | Replace the prototype's strain test and path search with steps 3 and 6 here, add the line check, then post the stage 1 predeclaration |
| [PR #890](https://github.com/adubey/mantle-bloom/pull/890), [PR #896](https://github.com/adubey/mantle-bloom/pull/896) | Keep as the base | Their kernels, recorder and rejection log are reused |
| [#894](https://github.com/adubey/mantle-bloom/issues/894) | Complete when this revision is accepted | Its summary became this revision |
| [#846](https://github.com/adubey/mantle-bloom/issues/846) production implementation | Blocked | Waits for D1 |
| [#797](https://github.com/adubey/mantle-bloom/issues/797) tracking | Update | D1 gains the line check |
| Production great-circle splitter | Unchanged | It remains the D2 baseline with its own daughter rule |

**2026-11-14 — load proxy C adopted after D0.** D0 ([PR #887](https://github.com/adubey/mantle-bloom/pull/887))
showed that the specified local load proxy failed three of its four statistical criteria and one
manufactured case. A line-cut proxy through the balanced load passed all of them, at a tenth of the
local proxy's cost. The owner chose to adopt it instead of the predeclared fallback of periodic
thin-shell solves. Step 1 of the Design now specifies C, with `k_cut = 1.33` replacing `L_drag`
and `L_pull`. D1 adds shadow shell solves as a check. D0's criteria are unchanged, and no other gate
changes. [#888](https://github.com/adubey/mantle-bloom/issues/888) runs D1.

**2026-11-14 — accepted.** The owner merged [PR #880](https://github.com/adubey/mantle-bloom/pull/880) and chose to run D1 in two stages: a coarse search on 5 development worlds, then the best settings on all 20. D1 and D2 state their 1 Myr step size, and D1 adds a 0.5 Myr check of the frozen setting, as the RFC lifecycle requires.

**2026-11-14 — damage bands replace the eigensolve candidate search (revision proposed after #848 and #876).**
The previous version chose a split by a weighted-Laplacian graph cut: up to twelve candidate cuts per
plate, each refined by repeated sparse eigensolves, filtered by structural checks, scored by an
`O/R/S` breakup prior and drawn as one keyed stochastic event. Failed attempts thinned a band along
the drawn path.

The evidence against keeping it:

- **Runtime ([#848](https://github.com/adubey/mantle-bloom/issues/848), [`analysis/issue848/`](../analysis/issue848/README.md)).**
  Eigensolves are about 85% of search time. A perfect exact skip of opportunities with no valid
  candidate still leaves about 1.84 times the baseline step time. Lower refinement caps change
  partitions and `P_break`, and the one exact shortcut saves 16% of solves.
- **Validity is geometric ([#876](https://github.com/adubey/mantle-bloom/issues/876), [`analysis/issue876/`](../analysis/issue876/README.md)).**
  Structural validity is unrelated to extensional load relative to strength (AUC 0.47–0.56), so
  valid cuts are bookkeeping candidates, not physical ones. Valid cuts carry about 7% of column
  capacity under one step's load.
- **Damage needs softening ([#699](https://github.com/adubey/mantle-bloom/issues/699), [`analysis/issue487/`](../analysis/issue487/README.md)).**
  Damage that does not weaken rock fragments plates. The absolute strength threshold depends on an
  unconstrained background stress, so it must be calibrated.

The owner then set the goal and rules for this revision in [#876](https://github.com/adubey/mantle-bloom/issues/876#issuecomment-6874579897):
efficient, realistic-looking rifting rather than new science; parameters with literature defaults
and a tuning search; loads from mantle drag, slab pull and plumes; bands that start about
100–150 km wide and grow; physical candidates that can start anywhere; and gates on realism and
runtime. This revision replaces the candidate search, prior and draws with rift strain stored in
each cell. Strain grows where load exceeds a softened strength and thins the band in place, and a
split follows the cheapest path through the remaining strength once the intact ligament fails. It
keeps the minimum age, cooldown, mask-based split helper, trace-band craton diagnostics, tiny-breakup
definition, held-out sample and the G3 limits.

| Gate | Previous criterion | Proposed criterion | Reason |
|---|---|---|---|
| G0 | Candidate geometry, opening, availability and solve counts on 20 development worlds for direct-cost weighting | **D0:** calibrate the cheap directional load tensor against the thin-shell diagnostic: manufactured sign cases, then sign agreement, principal direction, rank correlation and magnitude ratio on development worlds; otherwise switch to periodic shell solves | The candidate search is removed; the load proxy is the new uncertain input |
| G1 | `O/R/S` prior monotonicity and path-distribution checks | **D1:** predeclared parameter search on development worlds against the realism targets, including the #699 fragmentation check | The stochastic prior is removed; tuning replaces probability calibration, as the owner directed |
| G2 | Ten revised held-out pairs (registry 12–20 plus reserve `1479220500`); at most one `D_s < 0`; at least 50 combined actual splits, with the frozen reserve-list fallback | **D2:** same ten seeds, same tiny-breakup non-worsening limit, same combined 50-split count and reserve fallback, plus every realism target | Owner gates on realism; the held-out sample is still unrun for any revised method |
| G3 | Runtime within 2%; craton loss at most 2% per Gyr; 10–40 plates; half-step plate count, craton loss and tiny-breakup share within 20% (zero-share rule) | **D3:** all accepted limits unchanged, with breakup count and rift duration added to the half-step comparison | Owner gates on runtime; limits unchanged |

| Item | Verdict | Reason |
|---|---|---|
| G0–G2 results for either weighting | Historical | They evaluate the removed candidate method |
| G3 | Re-scoped as D3 | Limits unchanged; run after implementation |
| [#846](https://github.com/adubey/mantle-bloom/issues/846) production implementation | Blocked, re-scope | Implement the damage-band method after D0 and D1 |
| [#797](https://github.com/adubey/mantle-bloom/issues/797) tracking | Re-scope | Replace the G0–G3 checklist with D0–D3 |
| [#848](https://github.com/adubey/mantle-bloom/issues/848) | Close after a return comment | Every measured implementation and tested exact shortcut stays far above G3's limit; untested equivalent optimizations are not ruled out, but the method they would speed up is replaced. Before closing, post a consolidated return on #848: the runtime findings, the failed-rift finding (the experimental callback never applied failed-rift thinning; failed rifts are now bands left in the rock) and the G0–G3 to D0–D3 gate change |
| [#836](https://github.com/adubey/mantle-bloom/issues/836), [#866](https://github.com/adubey/mantle-bloom/issues/866) | Addressed in part | Interior localization is the proposed fix for wide subsidence; margin `_open_rift` is unchanged |
| [#699](https://github.com/adubey/mantle-bloom/issues/699) | Evidence | Its fragmentation result motivates softening; no further gate depends on it |
| RFC 0010 | Unchanged | Remains accepted until this RFC's gates and lifecycle decision settle |

**2026-11-09 — weak-path weighting, failed-rift outcomes and runtime evidence (revision proposed in #848).**
The accepted method used `conductance = 1 / edge_cost` in its weighted graph Laplacian even though
`edge_cost` is the resistance to cutting. Since the Laplacian quadratic form is
`xᵀLx = Σ wᵢⱼ(xᵢ - xⱼ)²`, this mapping penalizes separating the low-resistance edges that the
proposal intends to cut. In the controlled comparison in
[#850](https://github.com/adubey/mantle-bloom/issues/850), the inverse mapping crossed both
planted weak links in 0/100 cycles; direct cost crossed both in 100/100. Across 68 frozen G0
opportunities, direct cost shared no candidate partition fingerprints with the inverse mapping,
changed the selected partition in all 68, and selected the same anchor in 9/68. Direct cost
therefore matches the weak-link direction better, but it is not a validated final method: its
selected development cuts had median net opening `−0.059 cm/yr`, so geometry, availability and
the accepted `O/R/S` prior must be rechecked before held-out runs. The direct-cost snapshot
comparison made 6,072 sparse eigensolves in 259.6 s versus 6,481 in 311.6 s for inverse cost;
these are diagnostic timings, not full-simulation runtime evidence. [PR #851](https://github.com/adubey/mantle-bloom/pull/851)
records the reproducible comparison.

The first ten G2 pairs had median wall times of 8,363.5 s for the graph callback and 505.2 s
for baseline (16.25×), but the run record lacks sufficient isolation metadata to call this a
controlled G3 failure. The graph callback also returned on a failed breakup without applying
the RFC-required thinning and craton diagnostics along a selected path. This means G2 remains
evidence for the experimental callback and the old candidate method only; it does not establish
runtime or full failed-rift behavior for the revised method. The accepted G3 runtime limit
remains unchanged until a controlled paired measurement exists. No output-preserving runtime
optimization has been demonstrated, so the candidate method's ability to meet the 2% limit
remains unknown.

The inverse weighting conflicts with principle 8 because it steers the partition away from
the weak paths the breakup process is intended to use. The solve cost raises principle 2's
runtime constraint, and principle 9 requires controlled evidence before calling G3 failed or
passed.

This revision proposes `edge_cost` as the positive, dimensionless Laplacian adjacency weight,
preserving the cost meaning: lower-cost edges carry a smaller separation penalty. It also makes
the event outcome explicit. For candidate `i`, let `qᵢ` be its proposal probability and `pᵢ` its
frozen breakup prior; let `q₀` be the probability of proposing no structurally valid path. Draw
breakup with `P_break = Σᵢ qᵢpᵢ`. On breakup, retain the conditional path distribution
`qᵢpᵢ / P_break`. On no breakup, choose a failed path with conditional mass
`qᵢ(1−pᵢ) / (1−P_break)`, or choose no path with mass `q₀ / (1−P_break)`. The same path then
drives the split or failed-rift thinning and craton diagnostics. The no-path outcome makes no
geological mark. This is the candidate-first proposal-and-outcome model written in event-first
form; it preserves the successful path distribution and gives a failed attempt a path without
applying a second independent success roll. Event and path draws use separate keyed streams.

| Item | Verdict | Reason |
|---|---|---|
| G0 | Must rerun | The Laplacian weighting changes candidate partitions and selected traces; remeasure geometry, candidate availability, opening, score and runtime on development worlds. |
| G1 | Must rerun | Revalidate the frozen `O/R/S` prior and structural checks against the revised candidate distribution, and test failed-path selection and no-path outcomes. |
| G2 | Must rerun | The existing pass is for the inverse-weight experimental callback and omits failed-rift effects. Use only still-unrun held-out seeds for the revised method after the revised sample is frozen. |
| G3 | Blocked | Run after the method and failed-rift semantics are implemented; use controlled paired full simulations. Keep the 2% runtime limit pending valid measurements. |
| #846 | Blocked | Keep production implementation paused until the revision and required G0/G1 gates are accepted and pass. |
| #847 | Unchanged as historical evidence | It records the old method's G0–G2 results; it is not validation for this revision. |
| #850 / #851 | Unchanged as completed diagnostic evidence | The controlled weighting result stands; it motivates a revised candidate but does not validate it. |
| #797 | Re-scope after revision decision | Update the gate checklist and implementation dependency after the owner decides on this revision. |
| RFC 0010 | Blocked | It remains accepted until the revised RFC 0018 gates and lifecycle decision are settled. |

The first ten G2 seeds remain excluded from tuning and revised-method validation. Of registry
seeds 11–20, nine remain wholly unrun. Seed `205240614` has a completed baseline-only run and
should not be described as unseen; this revision excludes it from a clean revised-method
sample. To retain ten pairs, the proposed new sample is the nine untouched registry seeds
12–20 plus the first frozen reserve seed, `1479220500`, subject to the owner's approval of this
revised gate. The remaining 29 reserve seeds stay available under the preregistered reserve
rule. No seed in the proposed sample is used for development or tuning.

| Gate | Previously accepted evidence or criterion | Proposed revision | Reason |
|---|---|---|---|
| G0 | Passed on 20 development worlds for inverse weighting | Repeat on 20 development worlds with direct-cost weighting; compare with uniform control and report candidate geometry, opening, availability and solve cost | Direct cost changed the selected partition in all 68 measured opportunities, so the old geometry and runtime evidence does not apply. |
| G1 | Passed controlled probability and structural checks for the old candidate flow | Repeat the `O/R/S` checks and validate success-path, failed-path and no-path distributions | The candidate paths change, and the old callback had no failed-path outcome. |
| G2 | The amended gate evaluated the first ten registry seeds; at most one of ten could have `D_s < 0` and the pair set needed at least 50 actual splits | For the revised method, evaluate seeds 12–20 plus frozen reserve seed `1479220500`; retain the one-of-ten and 50-split limits | The first ten have already evaluated the old method. Seed 11 (`205240614`) has a baseline-only run, leaving nine wholly unseen registry seeds. This sample change needs owner approval before use. |
| G3 | Median runtime within 2% of paired production baseline; existing craton, plate-count and timestep limits | Keep all limits. Use paired, serialized, same-host and same-settings full runs with profiling separate | Existing wall times are not controlled enough to justify changing the 2% limit. |

**2026-11-08 — G0–G2 diagnosis results recorded.** Revised G0 and the coupled-probability
G1 diagnostics pass on the 20 preregistered development worlds. G2 passes under the
owner-directed 10-pair amendment merged in [PR #845](https://github.com/adubey/mantle-bloom/pull/845):
all first 10 frozen pairs completed, zero had `D_s < 0`, and they produced 255 combined actual
splits. The G2 amendment was made after those ten outcomes were observed, so it is post-result
evidence with lower precision and supports non-worsening on the evaluated pairs, not positive
efficacy. Production implementation is tracked in [#846](https://github.com/adubey/mantle-bloom/issues/846). Gate details and compact evidence are in
[`analysis/issue797/README.md`](../analysis/issue797/README.md). Production implementation,
implementation-level ledger checks, G3 and the remaining acceptance criteria are still open.

**2026-11-08 — G2 sample-size amendment.** The owner requested reducing the G2 held-out
sample from 20 pairs to 10 because the graph runs are taking a long time. The first ten
seeds in the frozen registry are the evaluation set; all ten paired runs are complete. Keep
the existing non-worsening tolerance at 10% of the sample, rounded to one seed: at most one
of ten may have `D_s < 0`. Keep the combined 50-actual-split requirement and all other
simulation settings and outcome definitions unchanged. The ten pairs produced 255 actual
splits (157 baseline, 98 graph), zero tiny splits in either arm, and `D_s = 0` for every pair.
Thus the amended criteria pass on these ten pairs.

This is a post-result, owner-directed amendment: the first ten outcomes were observed before
the sample size was reduced. It is not preregistered or confirmatory, and the smaller sample
has less precision. The result supports non-worsening on these ten pairs; it does not show
positive efficacy. The ten-pair results and amendment are recorded in
[`analysis/issue797/g2-sample-size-amendment.json`](../analysis/issue797/g2-sample-size-amendment.json).

| Gate | Previous criterion | Revised criterion | Reason |
|---|---|---|---|
| G2 sample and non-worsening threshold | Complete all 20 frozen held-out pairs; at most two may have `D_s < 0`; at least 50 combined actual split events | Complete the first 10 seeds in the frozen registry; at most one may have `D_s < 0`; at least 50 combined actual split events | The owner requested a shorter evaluation because graph runs take a long time; preserving the 10% tolerance gives a limit of one worsened pair out of ten. The change was made after the first ten outcomes were observed. |

**2026-11-06 — G2 non-worsening criterion amendment.** This revision is proposed in
[PR #824](https://github.com/adubey/mantle-bloom/pull/824). The owner approved replacing the
positive-effect requirements after the first four held-out pairs had `D_s = 0` in every pair:
neither method produced a tiny breakup (66 baseline splits and 38 graph splits, zero tiny
events in either arm). The previous amendment required a positive mean
`D_s`, a paired bootstrap interval wholly above zero, and improvement on at least two seeds.
At that stage, the criterion passed when no more than two of the 20 held-out seeds had `D_s < 0`;
`D_s = 0` counts as non-worsening. This permits a result with no improvement if the graph
method also avoids adding tiny breakups. This is an owner-selected tolerance, not a geological
rate threshold. The seed list and simulation settings remain fixed; all four observed pairs
remain in the analysis. Because the criterion changed after four outcomes were observed,
label the result as a second amended G2 analysis, not as preregistered. Continue to report the
mean `D_s`, its 95% paired seed-bootstrap interval, and the number of seeds with `D_s > 0` as
descriptive evidence, but they are no longer pass conditions.

| Gate | Previous criterion | Revised criterion | Reason |
|---|---|---|---|
| G2 outcome threshold | Mean `D_s > 0`, two-sided 95% paired seed-bootstrap interval above zero, and `D_s > 0` on at least two seeds | No more than two of the 20 held-out seeds have `D_s < 0`; `D_s = 0` is non-worsening. Complete all 20 pairs and meet the existing 50-event requirement. | The first four pairs showed zero tiny events in both methods. Requiring positive improvement and a positive interval would reject ties even though the owner's goal is to avoid making tiny breakups worse. |

**2026-11-06 — G2 tiny-breakup endpoint amendment.** This revision merged in
[PR #823](https://github.com/adubey/mantle-bloom/pull/823). The owner approved this change after
the first two held-out seed pairs exposed a zero-baseline failure: both methods had zero tiny
breakups, so the relative reduction divided by a zero baseline share and was undefined. The
relative share also conditions on the number of splits, which can differ between methods
under the same 200 Myr exposure. The two pairs had 42 baseline splits and 19 graph splits,
with zero tiny breakups in either arm (`S_baseline = 0`, `S_graph = 0`, so the share
difference is 0 percentage points); their paired tiny-count differences are both zero.
These are preliminary outcomes from the G2 diagnostic run associated with [#797](https://github.com/adubey/mantle-bloom/issues/797). The seed
list and simulation protocol remain fixed. Because the endpoint changed after two pairs were
observed, report the revised result as an amended G2 analysis, retain the originally
preregistered share-based result as a secondary analysis, and do not describe the revised
endpoint as preregistered. Keep both completed pairs in all summaries and uncertainty
calculations.

| Gate | Earlier criterion | Revised criterion | Reason |
|---|---|---|---|
| G2 primary outcome | At least one-third relative reduction in the share of actual splits with a smaller daughter under 7% of parent area, and improvement on at least two seeds | For each seed, `D_s = T_baseline,s - T_graph,s`, the difference in tiny-breakup counts over 200 Myr. Pass when the mean paired difference is positive, its two-sided 95% paired seed-bootstrap interval excludes zero in the beneficial direction, and `D_s > 0` for at least two seeds. Use 10,000 resamples of whole seed pairs. | Relative reduction is undefined at a zero baseline share. Paired count differences are defined even when an arm has no splits or no tiny splits, and measure event burden at the fixed world exposure. |
| G2 supporting measures | Tiny-breakup share and daughter mass distribution | Report the original relative share reduction, percentage-point difference `S_baseline - S_graph`, and smaller-daughter area and mass distributions as secondary measures. | Preserve comparison with #770 while showing conditional split quality and daughter sizes. |
| G3 timestep comparison | Tiny-breakup share within 20% of the normal-step result | Keep the 20% relative comparison when the normal-step share is nonzero; if it is zero, require zero tiny breakups at half step. Report paired counts and shares. | Avoid dividing by a zero normal-step share. |

The current outcome threshold in the later 2026-11-06 amendment, further amended on
2026-11-08 for sample size, supersedes the intermediate positive-effect criterion in this
table. The intermediate endpoint and criterion were both
changed transparently after held-out outcomes had been observed.

The G2 endpoint is an amended analysis because it was selected after two seed-pair outcomes
were observed. The original relative-share result remains reported for transparency. The
seed list and the minimum of 50 actual split events across the seed set, with its frozen
reserve rule, remain unchanged.

**2026-11-04 — G1 gate review.** This revision proposed [PR #810](https://github.com/adubey/mantle-bloom/pull/810). The previous G1 required at least 10 observed opportunities
and a usable candidate in at least two thirds of trigger-positive opportunities. The 2/3 rate
was a project-selected usability cutoff, not a geological estimate. Geological rift systems
are not independent draws from the simulation's per-trigger opportunities, and observations
show both successful and failed systems whose outcomes depend on evolving force balance,
weakening and competition. The literature does not establish a universal success fraction for
this denominator; the mechanisms separating successful and failed rifts remain debated
[[Brune et al., 2023]](#ref-brune2023), and mapped failed systems illustrate the variety of
incomplete rift histories [[Elling et al., 2022]](#ref-elling2022).

The old G1 rate criterion is removed. G1 now checks structural invariants and the stated
probability model's behavior on controlled cases, while development-world outcome counts are
reported descriptively with uncertainty. No aggregate per-opportunity fraction is called a
geological pass threshold. G1 does not claim geological rate calibration; that would need a
separate mapping between simulated opportunity exposure and geological rift episodes.

## References

- <a id="ref-brune2016"></a>Brune, S., Williams, S. E., Butterworth, N. P., & Müller, R. D.
  (2016). Abrupt plate accelerations shape rifted continental margins. *Nature*, 536, 201–204.
  [doi:10.1038/nature18319](https://doi.org/10.1038/nature18319)
- <a id="ref-brune2023"></a>Brune, S., Kolawole, F., Olive, J.-A., et al. (2023).
  Geodynamics of continental rift initiation and evolution. *Nature Reviews Earth & Environment*,
  4, 235–253. [doi:10.1038/s43017-023-00391-3](https://doi.org/10.1038/s43017-023-00391-3)
- <a id="ref-buck1991"></a>Buck, W. R. (1991). Modes of continental lithospheric extension.
  *Journal of Geophysical Research*, 96(B12), 20161–20178.
  [doi:10.1029/91JB01485](https://doi.org/10.1029/91JB01485)
- <a id="ref-buiter2014"></a>Buiter, S. J. H., & Torsvik, T. H. (2014). A review of Wilson Cycle
  plate margins: A role for mantle plumes in continental break-up along sutures? *Gondwana
  Research*, 26, 627–653. [doi:10.1016/j.gr.2014.02.007](https://doi.org/10.1016/j.gr.2014.02.007)
- <a id="ref-courtillot1999"></a>Courtillot, V., Jaupart, C., Manighetti, I., Tapponnier, P., &
  Besse, J. (1999). On causal links between flood basalts and continental breakup. *Earth and
  Planetary Science Letters*, 166, 177–195.
  [doi:10.1016/S0012-821X(98)00282-9](https://doi.org/10.1016/S0012-821X(98)00282-9)
- <a id="ref-elling2022"></a>Elling, R., Stein, S., Stein, C. A., & Gefeke, K. (2022).
  Three major failed rifts in central North America: Similarities and differences. *GSA Today*,
  32(6), 4–11. [doi:10.1130/GSATG518A.1](https://doi.org/10.1130/GSATG518A.1)
- <a id="ref-huismans2003"></a>Huismans, R. S., & Beaumont, C. (2003). Symmetric and asymmetric
  lithospheric extension: Relative effects of frictional-plastic and viscous strain softening.
  *Journal of Geophysical Research*, 108(B10), 2496.
  [doi:10.1029/2002JB002026](https://doi.org/10.1029/2002JB002026)
- <a id="ref-mckenzie1978"></a>McKenzie, D. (1978). Some remarks on the development of
  sedimentary basins. *Earth and Planetary Science Letters*, 40, 25–32.
  [doi:10.1016/0012-821X(78)90071-7](https://doi.org/10.1016/0012-821X(78)90071-7)
- <a id="ref-vauchez1997"></a>Vauchez, A., Barruol, G., & Tommasi, A. (1997). Why do
  continents break-up parallel to ancient orogenic belts? *Terra Nova*, 9, 62–66.
  [doi:10.1111/j.1365-3121.1997.tb00003.x](https://doi.org/10.1111/j.1365-3121.1997.tb00003.x)

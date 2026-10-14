# RFC 0008: Cratons jam collisions: braking and fusion

- **Status:** Draft
- **Date:** 2026-10-14
- **Related:** #274 (cratons), #340 (craton partition rules, the craton brake), #318 and #315
  (collision polarity), #331 (dock terranes as cells); RFCs
  [0001](0001-merge-cap-craton-loss.md) (merges must not destroy cratons),
  [0002](0002-craton-trench-attrition.md) (cratons resist sinking),
  [0003](0003-rifting-around-cratons.md)

## Problem

When a continent carrying a craton reaches a plate boundary, the model gives the craton almost
no say in what happens next. Three behaviours follow from the code as it stands.

**Plates don't slow down.** Plate motion is a torque balance (`torque.shift_plate`): slab
pull, ridge push, basal drag and collision friction. A craton enters it in one place: the
collision brake `torque.CRATON_COLLISION_FRICTION_GAIN` scales friction at contested
continent-continent nodes by the craton strength of the *weaker* side
(`torque.contact_craton_m`, the lesser of the two sides' `craton_crust_m`). So:
- only craton against craton brakes. A craton pushing into ordinary crust adds nothing; the
  weak side retreats and the relief term (`COLLISION_RELIEF_FRICTION_GAIN`) carries the
  resistance;
- nothing brakes before contact, however close the craton is;
- a craton reaching a trench, rather than another continent, never brakes. The collision
  mask needs a continental plate with a continental neighbour;
- **the brake is inert anyway.** `collision_drag_coefficients` puts
  `COLLISION_FRICTION_REFERENCE_PA` (5 × 10⁷) where `basal_drag_coefficients` puts
  μ/d_s = 5 × 10¹⁴ Pa·s/m: the coefficient multiplies the physical relative velocity in m/s.
  At 4 cm/yr (1.27 × 10⁻⁹ m/s) that is 0.063 Pa of traction, and the overlap, relief and
  craton gains raise it by at most 3 × 4 × 4 = 48×. Basal drag at the same slip rate is
  6 × 10⁵ Pa. Collision friction, craton gain included, is about 10⁷ times too weak per unit
  area to matter.

**Cratons pull plates into trenches.** Slab pull (`slab_pull_torque`) acts on an oceanic
plate's converging boundary nodes and scales with each node's own Hm, at the mantle
lithosphere's excess density (3,300 − 3,250 kg/m³). It ignores the node's crust. An oceanic
plate carrying a continental terrane or a craton therefore pulls that block into the trench at
full strength, and a cratonic keel, being thicker, pulls *harder* than ocean floor. On Earth
it is the other way round: continental crust thicker than ~15 km makes the lithosphere buoyant
enough to jam a subduction zone [[Cloos, 1993]](#ref-cloos1993). Cratonic mantle is
chemically depleted, which offsets its cold, dense thermal state
[[Jordan, 1978]](#ref-jordan1978).

**Jams don't end in fusion.** Since #340, boundary retreat never consumes a craton cell
(`cratons.retreat_allowed`). Where both sides of a collision front are cratonic, the front
jams. Where one side is, the front falls back to the other side's retreat (the polarity
fallback in `lithosphere_plate.boundary_context`). Fusion, meanwhile, ignores cratons:
- a pair counts as colliding only while it is close *and* closing faster than
  `TRANSFORM_RATE_THRESHOLD` (1 cm/yr) at 4 or more nodes
  (`merge_split.find_continental_collision_pairs`);
- its accumulated collision time is **dropped** the first step it stops closing
  (`update_collision_progress`). A craton jam that brakes the plates below 1 cm/yr resets
  the clock that would have fused them;
- the merge roll (`_merge_probability`) depends on the pair's combined size and relative
  rotation rate only;
- the backstop forces a fusion only after 30% node overlap for 30 Myr
  (`FORCED_MERGE_OVERLAP_FRACTION`, `FORCED_MERGE_SUSTAINED_YEARS`). By then cratons can be
  stacked on each other, which is what RFC 0001 found destroying a quarter of the world's
  cratons in one merge;
- merges are between continental plates only. A cratonic terrane on an oceanic plate can't
  retreat, and attrition skips terrane cells (the `trench` mask in `quad_tectonics._retreat`),
  so it sits jammed at the trench with nothing to resolve it but issue #331's docking.

Nothing records any of this. The craton ledger books losses by sink when they happen. It
can't say whether a loss happened at a contact that jammed and didn't fuse, and
`World.collision_progress` holds only the pairs colliding now.

**Evidence.** Seeds 7, 11 and 2024 ran to 1 Gyr (`bin/debug/compare_cratons.py`, origin/main
at 1b4f686, 5 Myr steps; `analysis/rfc0008/`). Sinks are given as shares of the starting
craton:

| Seed | Destroyed, % of start craton | Delaminated | Rifted | Eroded | Subducted | Collision reworked | Merges / splits |
|---|---|---|---|---|---|---|---|
| 7 | 84% | 54% | 13% | 10% | 5.9% | 1.2% | 55 / 115 |
| 11 | 65% | 18% | 9.0% | 21% | 14% | 2.1% | 27 / 109 |
| 2024 | 67% | 31% | 16% | 10% | 9.3% | 1.3% | 36 / 104 |

- Cratons lose two thirds or more of their volume in 1 Gyr.
- The largest sink is `delaminated_m3`. It arrives in lumps: up to 15% of the starting
  craton in one 100 Myr window on seed 7, and 11% on seed 2024. That is the signature RFC 0001
  found for the merge cap, which books its clips to `delaminated_m3`, but the ledger can't
  separate the cap from root delamination until 0001's gate splits it.
- These runs include #358, which made the split cooldown 2 Myr at any step length. At 5 Myr
  steps that tripled splits and merges against earlier tuning runs, so the numbers aren't
  comparable with #378, #385 or #388.

What the numbers don't show is how much of that loss happened at jammed contacts, how often
jams reset collision time, or how much slab pull comes from buoyant columns. That is the
first gate (Design).

**Owner's decision.** The project owner has decided that both effects are to be modelled: a
craton approaching a plate boundary slows the plates, and promotes their fusion. This RFC
designs how, and its gate measures the baseline and chooses between variants. It doesn't
decide *whether*.

## Proposed solutions

**What happens in nature.**
- **Collisions slow plates.** India's convergence with Asia dropped from ~15 to ~4 cm/yr
  between 50 and 35 Ma, coincident with the onset of collision. A force balance attributes
  the slowdown mainly to the resistance of the thickening Himalaya–Tibet, ~5–6 × 10¹² N/m
  [[Copley et al., 2010]](#ref-copley2010). Others add the loss of slab pull as continental
  lithosphere entered the trench [[van Hinsbergen et al., 2011]](#ref-vanhinsbergen2011).
  Convergence slowed again, by more than 40%, between 20 and 10 Ma
  [[Molnar & Stock, 2009]](#ref-molnar2009).
- **Slab pull matters.** Attached upper-mantle slabs supply about half of the driving force
  on plates [[Conrad & Lithgow-Bertelloni, 2002]](#ref-conrad2002), so losing it at a trench
  is a first-order change in a plate's force balance.
- **Buoyant lithosphere jams subduction.** 100 km lithosphere with 30 km of continental crust
  is ~0.09 g/cm³ less dense than the asthenosphere and resists subduction. Only crust thicker
  than ~15 km makes the lithosphere buoyant enough to jam a subduction zone. A jam rearranges
  plate motions, generally with a new subduction zone elsewhere [[Cloos, 1993]](#ref-cloos1993).
- **After a jam, the slab detaches.** Continental collision is followed by attempted
  subduction of buoyant continental lithosphere, and the dense oceanic slab breaks off
  [[Davies & von Blanckenburg, 1995]](#ref-davies1995). Subduction at that contact ends, and
  the two continents are joined: the suture of a Wilson cycle
  [[Wilson, 1966]](#ref-wilson1966).
- **Cratons are the strong blocks.** Cratonic lithosphere is thick, buoyant and stiff
  [[Jordan, 1978]](#ref-jordan1978), and survives because weaker mobile belts around it
  take the strain [[Lenardic et al., 2003]](#ref-lenardic2003).

The literature ties slowing to buoyant continental lithosphere at the trench and to the
resistance of the collision zone, not to cratons specifically. A craton is the extreme case:
the thickest, most buoyant column, and one that never subducts.

**Slowing.**
0. **Make contact resistance dimensionally real (required).** Replace the collision
   brake's reference stress with a drag coefficient in Pa·s/m, sized so that a fully jammed
   contact at the post-collision India–Asia convergence rate resists with
   [[Copley et al., 2010]](#ref-copley2010)'s Himalaya–Tibet force, ~5 × 10¹² N/m of trench.
   Without this, neither solution below can slow anything. It also switches on the existing
   relief term for every continental collision, which today does nothing.
1. **Buoyant columns resist subduction (recommended).** Slab pull per node scales with the
   column's net negative buoyancy, not its Hm alone. A cratonic keel is taken as neutrally
   buoyant [[Jordan, 1978]](#ref-jordan1978), so it adds no pull. Continental crust thicker
   than the jam threshold [[Cloos, 1993]](#ref-cloos1993) stops pulling, and adds a
   resistance through solution 0's contact drag. The plate slows as its
   buoyant margin, and then its craton, reaches the trench. This is the mechanism the
   literature supports, and it fixes the cratonic-keel pull, which is backwards today.
2. **A proximity brake (recommended, behind the gate).** At collision contacts, the craton
   factor of the brake reads craton within `CRATON_BRAKE_REACH_KM` of the contact on
   *either* side, weighted by distance. That replaces the weaker side's craton at the contact
   node. A craton approaching a front then slows it before it arrives, and a craton
   indenting weak crust brakes too.
   - It departs from the reasoning in `torque.py` that collision resistance is set by the
     weak side's deformation [[Copley et al., 2010]](#ref-copley2010), which the relief term
     already carries. Its basis is a design choice: a strong block close to the front leaves
     a narrower weak belt to take the shortening.
   - It risks double counting the relief term. It is built only if the gate shows contacts
     where a craton indents weak crust without slowing.
3. **Raise `CRATON_COLLISION_FRICTION_GAIN`.** It multiplies an inert coefficient, so on
   its own it does nothing. Even after solution 0, it acts only craton on craton, at
   contact, and misses trenches and approach.
4. **Basal drag under keels.** A thick keel adds drag on the asthenosphere. Contested:
   tectosphere sped plates up in some configurations (the `torque.py` comment on
   `CRATON_COLLISION_FRICTION_GAIN`). It also slows a cratonic plate everywhere, not near a
   boundary. Not proposed.

**Fusion.**
5. **A jammed craton contact keeps its collision time (recommended).** A pair stays tracked
   while it is in contact and craton-jammed, even when the brake has stopped it closing.
   Only separation or divergence drops it.
6. **Jammed contacts lock (recommended).** A craton-jammed pair fuses at a hazard rate set
   by how much of its contact is jammed, with mean time `CRATON_LOCK_MYR` at a fully jammed
   contact. The supercontinent size gate doesn't apply to it: a locked boundary is locked
   whatever the plates' size. Slab breakoff and suturing on Earth take a few to tens of Myr
   after the buoyant block arrives, and India and Asia still converge 50 Myr after first
   contact, so this is a calibration knob, not a measured rate.
7. **Raise the merge probability with craton contact.** A craton term in
   `_merge_probability`. Simpler, but it keeps today's step-keyed roll, which depends on step
   length, and it doesn't stop the jam resetting collision time. Fallback only.
8. **Cratonic terranes on oceanic plates.** These can't merge, since the carrier is oceanic.
   They are resolved by #331's docking (the terrane's cells pass to the overriding plate,
   craton intact), which this RFC depends on rather than duplicates. Solution 1 slows the
   carrier as the terrane arrives.

**Recommendation:** a measurement gate, then 0, 1, 5 and 6, with 2 if the gate shows
unbraked craton indentation. Fusion at jam onset is also what RFC 0001's solution 3 asked for ("stop
craton-on-craton overlap from building"). This RFC takes that over. RFC 0001's placement fix
and honest booking must land first or with it, since more fusions will carry cratons.

## Design

**Units.** Forces per unit trench length in N/m; stresses in Pa; thicknesses in m; densities
in kg/m³; rates in cm/yr (`mantle.cm_per_yr_to_rad_per_yr`); times in Myr; hazard rates in
1/Myr.

**Phase 1: measurement gate.** A new probe, `bin/debug/probe_craton_contacts.py`, replays
seeds 7, 11 and 2024 to 1 Gyr at 5 Myr steps, one at a time. Each step, for every pair of
plates in contact (any crust type), it records:
- **contact:** length, normal closing rate, and craton area within `CRATON_BRAKE_REACH_KM`
  of the contact on each side (the "approach" zone);
- **jam:** a contact node is jammed when it is a retreat candidate that holds craton, so
  `retreat_allowed` blocks it, or when both sides at the node are cratonic. A retreat
  candidate is read *before* the no-craton gate:
  - on an active polarized front, the frozen lower-plate eligibility after the run-size gate
    (`collision_polarity` computes it as `eligible`, then stores only
    `retreatable = eligible & retreat_allowed`). The prepass keeps it as a new mask field,
    `run_eligible`. `BoundaryContext.shrinkable` can't be used here: on a front where some
    lower cells can retreat, it is replaced by the craton-filtered `retreatable`, so a
    cratonic lower cell facing ordinary crust would satisfy neither clause;
  - elsewhere, `BoundaryContext.shrinkable`, which `_retreat` filters by
    `retreat_allowed` only afterwards;
- **merge state:** `collision_progress` before and after the step; whether it was dropped
  while the pair was still in contact and jammed; merge rolls and outcomes; overlap fraction
  and the craton-on-craton overlap volume;
- **craton loss near contacts:** the per-node drop in `craton_crust_m` over the step on cells
  within `CRATON_BRAKE_REACH_KM` of each contact. It is booked to the pair, with the sinks
  booked that step (the per-node booking method of `bin/debug/attribute_rift_boundary_loss.py`);
- **buoyant pull:** per oceanic plate, the share of slab-pull torque magnitude from
  continental-coded or cratonic subducting nodes;
- **contact force:** at each jammed contact, the resisting force per unit trench length
  the brake would apply, in N/m, to check against Copley's scale.

It reports, per seed:
- J, the share of all craton destruction within reach of a jammed contact whose pair didn't
  fuse within 50 Myr;
- how often collision time was dropped on a jammed pair;
- the buoyant-pull share;
- the mean closing rate at contacts with and without craton in reach. Equal rates are the
  "unbraked indentation" that triggers solution 2.

These are the baseline for the acceptance criteria and select the variants:

| Gate outcome | Build |
|---|---|
| Any | 1, 5 and 6 (owner's decision) |
| Closing rate at craton contacts within 10% of craton-free contacts after 1 | 2 as well |
| Buoyant-pull share under 1% on every seed | 1 still, for the keel fix; its resistance is not tuned |
| J under 5% on every seed | Build as above; record that fusion isn't a craton-survival lever and drop its loss acceptance criterion |

**Gate result (2026-10-14, #395).** Seeds 7, 11 and 2024 to 1 Gyr on origin/main 3b9d9c1
([analysis/rfc0008](../analysis/rfc0008/README.md#phase-1-gate-issue-395)):

| | Seed 7 | Seed 11 | Seed 2024 |
|---|---|---|---|
| J | 22% | 57% | 40% |
| Jammed share of continental collision length | 44% | 56% | 37% |
| Collision-time drops on jammed pairs (of all drops) | 28 of 106 | 22 of 54 | 53 of 111 |
| Buoyant share of slab pull | 6.4% | 4.5% | 4.0% |
| Closing rate at converging continental contacts, craton in reach / clear (cm/yr) | 1.70 / 1.59 | 1.84 / 1.12 | 1.48 / 1.61 |
| Craton lost in fusions, % of start craton | 54% | 18% | 31% |

- Build 0, 1, 5 and 6. Solution 1's resistance is tuned: the buoyant-pull share is above 1%
  on every seed. The loss criterion stays: J is above 5% on every seed, and jams cover far
  more than 5% of collision length, so the "jams are rare" falsification isn't met.
- Solution 2 is decided after 0 and 1, as the table says. On the baseline, craton contacts
  don't close more slowly (the brake is inert: a median 2.5–7.8 × 10⁴ N/m at jammed contacts).
- Craton lost in fusions equals each seed's `delaminated_m3` share. So the merge cap is
  almost the whole delamination sink, which answers the question RFC 0001's gate asks. 91
  of 101 fusions were of a pair jammed at the time, and the craton-on-craton overlap they
  fused sums to 53–90% of the starting craton. That is the baseline for the 75% overlap
  criterion. It also means solution 6 must not land before RFC 0001's placement fix.
- At the contacts' own closing rates (about 1.5 cm/yr) the contact drag gives a median of
  0.8–0.9 × 10¹² N/m at jammed contacts, just under the 10¹²–10¹³ acceptance band (p90
  1.8–2.2 × 10¹²). The force criterion is read after 0 and 1, at the slowed rates.

**Slowing: buoyant columns (solution 1).** For each node in `subducting_boundary_mask`, a
net buoyancy per unit area:

`B = g × [(ρ_m − ρ_a) × Hm × (1 − s) − (ρ_a − ρ_cc) × Hc_cont]`

- s = `cratons.strength(craton_crust_m)` ∈ [0, 1]: a full-strength craton's keel is neutral.
- Hc_cont is the node's Hc where it is continental-coded (`effective_is_continental_from_codes`),
  else 0. Oceanic crust is left out as today, since it turns to eclogite on subduction
  [[Cloos, 1993]](#ref-cloos1993).
- ρ_cc = `RHO_CONTINENTAL_CRUST` (2,700), ρ_m = `RHO_LITHOSPHERE_MANTLE` (3,300),
  ρ_a = `RHO_ASTHENOSPHERE` (3,250).

Slab pull per node becomes `max(B, 0) × SUBDUCTION_LSINK_M × ds`, the same form as today. For
ocean floor (Hc_cont = 0, s = 0) it is unchanged. A buoyant node (B < 0) pulls nothing.

*Jam mask.* Its own test, independent of the slab-pull activation threshold. A node is in
the jam mask when it is continental-coded, on the side being forced down, in the boundary
band of a neighbour (`dist_to_neighbor ≤ reach_rad`), and the contact isn't diverging:
normal closing rate ≥ −`TRANSFORM_RATE_THRESHOLD`. "Forced down" means a node of an
oceanic plate (a terrane or craton riding ocean floor), or a lower-plate node of a polarized
continental front (`collision_polarity`'s `lower` mask). The upper plate is never in the
mask; it feels the equal and opposite force (below). `subducting_boundary_mask` needs closing above
+1 cm/yr, so a jam mask derived from it would drop the resistance (and slab drag) as soon
as the brake slowed the contact below 1 cm/yr. Driving forces would then speed it back up,
and the contact would flicker. A band that runs through zero keeps the resistance while the
contact stalls, and the drag form below vanishes with the relative velocity, so it can't
hold a contact that is pulling apart. The jam factor is

`j = clip((Hc_cont − JAM_MIN_CRUST_M) / (JAM_FULL_CRUST_M − JAM_MIN_CRUST_M), 0, 1)`

so a column at 15 km of continental crust starts to jam and one at 30 km jams fully. Both
plates at the contact feel it: the neighbour's band nodes read j at their nearest node on
the jammed side (through `neighbor_node_index`, as `contact_craton_m` does), so the force
is equal and opposite.

**Contact drag (solution 0).** Collision friction and the jam resistance share one
coefficient in `collision_drag_coefficients`' implicit `b − K ω` form, now in Pa·s/m:

`η = CONTACT_RESISTANCE_N_PER_M / (W × v_ref)`, per node `c_i = η × g_i × A_i × R² / SECONDS_PER_YEAR`

- `CONTACT_RESISTANCE_N_PER_M` = 5 × 10¹² N/m, the Himalaya–Tibet resistance
  [[Copley et al., 2010]](#ref-copley2010).
- v_ref = 4 cm/yr (1.27 × 10⁻⁹ m/s), India–Asia convergence after the slowdown
  [[Copley et al., 2010]](#ref-copley2010).
- W = `BOUNDARY_FORCE_REACH_MULTIPLIER` × spacing × R, the band width the term acts over:
  375 km at density 1. Computed from the world's spacing, so the force per trench length
  doesn't change with node density.
- g_i = clip(max(j, relief_factor, c), 0, 1) per node, with c the craton factor (today's
  weaker-side strength, or solution 2's proximity factor). A contact with no buoyant
  lower-side crust, no relief and no craton resists with nothing until relief builds. That
  matches Copley, whose resistance comes from mountain building. One bounded factor replaces the
  three multiplicative gains, which were tuned against an inert coefficient. The overlap
  severity term is dropped: the torque already sums over more nodes for a wider contact.

At density 1, η ≈ 1.05 × 10¹⁶ Pa·s/m, about 21 times basal drag's μ/d_s, but over a 375 km band
rather than the whole plate. Integrated across the band, a fully jammed contact (g = 1)
resists with `η × W × v` = 5 × 10¹² N/m at 4 cm/yr, twice the ridge push Copley estimates
(~2.5 × 10¹² N/m), and in proportion at other rates. Being linear in velocity, it slows a
driven plate toward the speed where driving and resistance balance. It doesn't reverse it,
which the implicit solve guarantees.

*Scale check.* For the reference continental column (Hc = 35 km, Hm = 100 km, s = 0),
B = 9.8 × (50 × 10⁵ − 550 × 3.5 × 10⁴) ≈ −1.4 × 10⁸ Pa: buoyant, as Cloos finds. For ocean
floor at reference (Hm = 40 km), B ≈ +2 × 10⁷ Pa, today's pull. The gate probe reports the
contact force per trench length at jammed contacts, which must sit at 10¹²–10¹³ N/m at
convergence rates of a few cm/yr.

*Polarity.* The jam mask applies on the plate in `subducting_boundary_mask`, which today is
only oceanic plates. A continental plate's lower-plate side at a polarized front
(`collision_polarity`) already gets the collision brake. Extending slab pull to continental
plates' oceanic margins is out of scope (Assumptions).

**Slowing: proximity brake (solution 2, if the gate selects it).** For each contested node in
the collision mask, and each jam-mask node:

`c = max(own, other)`, where `own = max over own cratonic cells within R of s × exp(−d / R)`

and `other` is the same on the neighbour plate, with R = `CRATON_BRAKE_REACH_KM`. c ∈ [0, 1]
replaces `cratons.strength(contact_craton_m)` in the brake. It costs one KD-tree query per
contact node, against each plate's cratonic cells (`cKDTree` on the craton subset, built once
per plate per step).

**Fusion: keeping collision time (solution 5).** A pair is *jammed* when at least
`MERGE_MIN_CONTACT_NODES` of its contact nodes are jammed (the gate's definition, including
the `run_eligible` candidates on polarized fronts), and it isn't diverging: mean normal
closing rate ≥ −`TRANSFORM_RATE_THRESHOLD`. `update_collision_progress`
keeps (and keeps accumulating) a jammed pair's progress, even when fewer than 4 nodes close
faster than the threshold. Progress is dropped only when the pair is neither colliding nor
jammed.

**Fusion: locking (solution 6).** For a jammed pair, `f` is the jammed share of the contact
length, in [0, 1]. The pair fuses as a Poisson process with rate `λ = f / CRATON_LOCK_MYR`.
Three things keep it step-robust.

- **Episodes, not steps, key the draw.** A jam *episode* starts when a pair becomes jammed and
  ends only after it has been unjammed for `CRATON_JAM_GAP_MYR` (5 Myr), so brief
  flickers don't start new episodes. Each pair counts its episodes, and episode k draws its
  threshold E ~ Exp(1) from `np.random.default_rng((seed, pair, k, LOCK_SEED_TAG))`. Nothing
  in the key depends on when the onset was detected, so a jam first seen at 2 Myr under one
  schedule and 5 Myr under another draws the same E. Keying by `round(onset_years)` would
  give unrelated thresholds.
- **Hazard only over the jammed interval.** Within a step from t₀ to t₁, f is taken as linear
  between its values at the step's start and end (f₀ = 0 if the pair wasn't jammed at t₀). So
  the step adds `ΔH = (f₀ + f₁) / 2 × years / CRATON_LOCK_MYR`. An onset inside the step is
  credited only for the part of the step after it, as linear interpolation places it,
  rather than the whole step.
- **Fusion times, and fusion in time order.** A pair whose H crosses E during the step has
  a fusion time t* in [t₀, t₁], from inverting the linear hazard. Every pair that locks in
  the step is fused in order of t*, through `_merge_plates_with_craton_audit`, before the
  ordinary roll. Locks are *not* limited to one per step. Before each fusion, the pair is
  revalidated: plate ids that an earlier fusion absorbed map to the survivor, and a pair
  that has become one plate is dropped. The ordinary roll and the forced merge keep their
  one-per-step limit.

With these, a fusion time agrees across step partitions to within one step: E ×
`CRATON_LOCK_MYR` / f for constant f, with independent locks resolving in the same order. What
can still differ is everything downstream. Positions, contacts and later jams evolve at step
resolution, so lock *counts* over a long run agree only statistically. The acceptance
tolerance is set for that.

The ordinary roll (`update_collision_progress`) is keyed by `round(elapsed_years)` and fires
once per step, so its probability per Myr depends on step length. That is out of scope here,
and noted for #313.

**Conservation.** Nothing here moves crust, craton or material. Slowing changes ω only.
Fusion uses the existing merge path, whose losses RFC 0001 fixes and books. Ledgers are
unchanged.

**Runtime.**
- Solution 1: per subducting node, one extra thickness read and the buoyancy arithmetic,
  inside `shift_plate`'s existing pass. O(boundary nodes).
- Solution 2: a KD-tree over each plate's cratonic cells (≤ 10⁵ at density 4, built once
  per step) and one query per contested node.
- Solutions 5 and 6: per jammed pair, a count over contact nodes the collision search
  already finds.

Under 2% of step time at density 4, to be profiled.

**Save compatibility.** One new `World` field, `craton_jams`: `{pair: (episode_index, hazard,
threshold, f_last, unjammed_since_years)}`. It is backfilled as `{}` on load, per `docs/save-compatibility.md`. An old
save's jammed pairs start their hazard on the first step after loading. No new per-node field.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `JAM_MIN_CRUST_M` | 15,000 m | Crust > ~15 km jams subduction [[Cloos, 1993]](#ref-cloos1993) | Fixed |
| `JAM_FULL_CRUST_M` | 30,000 m | 30 km continental crust on 100 km lithosphere resists subduction [[Cloos, 1993]](#ref-cloos1993) | Sensitivity at 25 and 35 km |
| Keel neutrality (factor `1 − s`) | full at s = 1 | Isopycnic tectosphere [[Jordan, 1978]](#ref-jordan1978) | Sensitivity with `1 − 0.5 s` |
| `CRATON_BRAKE_REACH_KM` (solution 2) | 500 km | Design choice: `CRATON_FORMATION_MARGIN_KM`, how far inside a margin cratons form | Sweep 250–1,000 km |
| `CONTACT_RESISTANCE_N_PER_M` | 5 × 10¹² N/m | Himalaya–Tibet resisting force [[Copley et al., 2010]](#ref-copley2010) | Sensitivity at 2.5 and 10 × 10¹² |
| v_ref | 4 cm/yr | Post-collision India–Asia convergence [[Copley et al., 2010]](#ref-copley2010) | Fixed (it trades off with the force) |
| `CRATON_JAM_GAP_MYR` | 5 Myr | Design choice: one tuning step, so a jam isn't split by a single unjammed step | Sensitivity at 2 and 10 Myr |
| `CRATON_LOCK_MYR` | 20 Myr | Calibration assumption: breakoff and suturing within tens of Myr of a buoyant block's arrival [[Davies & von Blanckenburg, 1995]](#ref-davies1995) | Sweep 10–50 Myr |
| `COLLISION_FRICTION_REFERENCE_PA`, the three friction gains | Retired | Replaced by `CONTACT_RESISTANCE_N_PER_M` and g | — |

Tune with the gate probe and `bin/debug/compare_cratons.py` on seeds 7, 11 and 2024 to 1 Gyr,
plus a regenerated seed 456862899 to 3 Gyr, one at a time. Watch:
- J and the craton sinks per Gyr;
- merges per Gyr, plate count, and the largest plate's share of the surface;
- closing rates at craton contacts, and mean plate speed.

## Detailed design (provisional touchpoints)

- `backend/app/torque.py`
  - `slab_pull_torque`: per-node buoyancy B from Hc, Hm, crust codes and craton (needs
    `craton_crust_m` in `BoundaryForceInputs`, new).
  - New `jam_mask(plate, inputs, subducting)` and `jam_factor(...)`: buoyant subducting
    nodes and j.
  - `collision_drag_coefficients`: the Pa·s/m coefficient η and the bounded factor g
    (solution 0); accept the jam mask and j; for solution 2, take c from new
    `craton_proximity(plate, neighbours, inputs)` in place of `contact_craton_m`.
  - New jam test independent of `subducting_boundary_mask`; the neighbour-side read of j.
  - `shift_plate`: wire both.
  - New constants `JAM_MIN_CRUST_M`, `JAM_FULL_CRUST_M`, `CRATON_BRAKE_REACH_KM`,
    `CONTACT_RESISTANCE_N_PER_M`, `CONTACT_RESISTANCE_REF_CM_PER_YR`, with comments citing
    this RFC; retire `COLLISION_FRICTION_REFERENCE_PA` and the three gains.
- `backend/app/collision_polarity.py`: keep `run_eligible` (the run-qualified lower-plate
  eligibility before the no-craton gate) on the frozen masks.
- `backend/app/merge_split.py`
  - `find_continental_collision_pairs`: also return jammed pairs with their jammed share f
    (new jam test, reusing the contact query).
  - `update_collision_progress`: keep progress for jammed pairs.
  - New `update_craton_jams(world, years, jammed)`: episodes, interpolated hazard,
    thresholds; returns the pairs that lock this step with their t*.
  - `apply_topology_changes`: fuse locking pairs in t* order, revalidating each, before the
    ordinary roll.
  - New constants `CRATON_LOCK_MYR`, `CRATON_JAM_GAP_MYR`, `_LOCK_SEED_TAG`.
- `backend/app/world.py`: `craton_jams` field; `backend/app/persistence.py` backfill if needed.
- `bin/debug/probe_craton_contacts.py` (new): the gate probe.
- `bin/debug/compare_cratons.py`: report merges by kind (roll, lock, forced) and closing rates
  at craton contacts.
- Tests:
  - `backend/unit_tests/test_torque.py`:
    - a reference ocean-floor node pulls exactly as today;
    - a continental node at 35 km pulls nothing and brakes; at 10 km it pulls (reduced) and
      doesn't brake;
    - a full-strength cratonic keel adds no pull;
    - j, c and g stay in [0, 1] over a grid of inputs;
    - the jam brake never reverses the relative velocity in one implicit step;
    - a fully jammed straight contact at 4 cm/yr resists with 5 × 10¹² N/m per unit trench
      length, within 1%, at densities 1 and 4;
    - a driven contact slowing through 1 cm/yr keeps its jam resistance: across the
      threshold, the resisting torque is continuous and the contact doesn't speed back up.
  - `backend/unit_tests/test_merge_split.py`:
    - a jammed pair keeps its progress when closing stalls; a separating pair loses it;
    - a mixed polarized front, where a cratonic lower-plate cell faces ordinary crust and
      other lower cells can retreat, counts that cell as jammed;
    - with f = 1, fusion happens at E × `CRATON_LOCK_MYR` to within one step, for one 10 Myr
      step, ten 1 Myr steps, and 5 Myr steps;
    - an onset inside a step draws the same E, and credits only the jammed part of the step,
      under 1, 2 and 5 Myr steps;
    - three disjoint pairs that all lock within one 10 Myr step all fuse in that step, in t*
      order, and a pair made stale by an earlier fusion is revalidated or dropped;
    - a lock ignores the size gate; a non-jammed pair still rolls as today.
  - `backend/unit_tests/test_collision_polarity.py`: `run_eligible` holds the cratonic
    lower-plate candidates that `retreatable` excludes.
  - `backend/unit_tests/test_persistence.py`: a save without `craton_jams` loads.
- Docs: `docs/simulation-model.md` (plate motion forces, "Merge and split", Cratons);
  `docs/debugging.md` (the probe).

## Assumptions and uncertainty

- **Cratons stand in for buoyancy and strength.** The literature ties slowing to buoyant
  continental lithosphere and to collision-zone resistance in general. This RFC uses
  continental crust thickness for buoyancy and craton strength for keel neutrality and
  proximity. Non-cratonic thick continental crust jams too, which is intended.
- **Model densities aren't Cloos's.** With the model's 2,700/3,300/3,250 kg/m³, the reference
  continental column is more buoyant than Cloos's ~0.09 g/cm³. The jam thresholds come from
  his crust thicknesses, not his densities.
- **Continental plates have no slab pull.** A continental plate with an oceanic leading edge
  (India before collision) isn't pulled by its slab in this model, so it can't *lose* slab
  pull when its continent arrives. Its slowing comes only from the brake. Fixing that changes
  every plate's driving forces and needs its own RFC.
- **The lock time is a guess.** Nothing published gives a mean time from craton contact to
  plate fusion, and the model's "merge" (two plates moving as one) is coarser than suturing.
  The sweep bounds it.
- **Turning the brake on changes every collision.** Solution 0 makes collision friction
  real for all continental collisions, not only cratonic ones. Plate speeds near collisions
  will drop on every seed, so the acceptance runs bound mean plate speed. Copley's force is
  one modern collision, and v_ref trades off against it.
- **The proximity brake is a design choice**, not a published mechanism, and it may double
  count the relief term. It's gated, and falsifiable (below).
- **Fusion raises plate-size variance.** Locks bypass the supercontinent size gate, so large
  plates may grow faster. The acceptance runs watch it.

## Acceptance criteria

- [ ] **Gate:** the probe's baseline (J, dropped collision time on jammed pairs, buoyant-pull
      share, closing rates with and without craton) is reported per seed in `analysis/`, and
      this RFC is updated with the variants it selects, before any solution is implemented.
- [ ] On seeds 7, 11 and 2024 to 1 Gyr, and seed 456862899 to 3 Gyr:
  - [ ] no jammed pair has its collision time dropped while still in contact and jammed;
  - [ ] the mean closing rate at contacts with craton in reach falls by at least 25% against
        `origin/main`, and no plate's speed drops below 20% of its pre-contact speed from
        braking alone (Copley's ~15 to ~4 cm/yr is about 27%);
  - [ ] if J ≥ 5% at the gate, craton destruction within reach of jammed contacts falls by at
        least half;
  - [ ] the craton-on-craton overlap volume at fusion falls by at least 75% (RFC 0001's
        upstream cause);
  - [ ] the time-averaged largest plate stays under 40% of the surface, and mean plate count
        changes by less than 25%.
- [ ] Ledgers balance (craton within 1e-9 of live volume; continental material closes).
- [ ] Fusion times are step-robust (the unit tests above). On a 1 Gyr replay of one seed at
      2 Myr and 5 Myr steps, lock counts agree within 20%; downstream topology diverges at
      step resolution, so exact agreement isn't expected.
- [ ] The resisting force at jammed contacts, from the probe, is 10¹²–10¹³ N/m of trench at
      the contacts' convergence rates.
- [ ] Step time at density 4 changes by less than 2%.
- [ ] Unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green; saves
      from before the change load.

## Falsification criteria

- **Jams are rare.** If the gate finds jammed contacts on under 5% of collision contact
  length on every seed, cratons rarely reach boundaries in the model, and the effects can't
  matter for craton survival. Build them as decided, but drop the loss criterion, and look
  for the loss elsewhere (RFCs 0002, 0003).
- **Slowing stalls the world.** If mean plate speed drops by more than 30%, or plates freeze
  at most contacts, even at the knobs' gentlest settings, the brake double counts. Drop
  solution 2 first; then lower j's reference stress.
- **Fusion builds a supercontinent.** If the largest plate exceeds 40% of the surface at any
  `CRATON_LOCK_MYR` in the sweep, locking is too strong. Require a minimum jammed contact
  length, or restore a size dependence.
- **Losses move.** If total craton destruction per Gyr falls by less than half of the drop
  near jammed contacts, cratons survive the jam and are lost after fusion (merge stacking:
  RFC 0001) or elsewhere. Revise toward that sink.

## References

- <a id="ref-cloos1993"></a>Cloos, M. (1993). Lithospheric buoyancy and collisional
  orogenesis: Subduction of oceanic plateaus, continental margins, island arcs, spreading
  ridges, and seamounts. *Geological Society of America Bulletin*, 105, 715–737.
  [doi:10.1130/0016-7606(1993)105<0715:LBACOS>2.3.CO;2](https://doi.org/10.1130/0016-7606(1993)105%3C0715:LBACOS%3E2.3.CO;2)
- <a id="ref-conrad2002"></a>Conrad, C. P., & Lithgow-Bertelloni, C. (2002). How mantle slabs
  drive plate tectonics. *Science*, 298, 207–209.
  [doi:10.1126/science.1074161](https://doi.org/10.1126/science.1074161)
- <a id="ref-copley2010"></a>Copley, A., Avouac, J.-P., & Royer, J.-Y. (2010). India-Asia
  collision and the Cenozoic slowdown of the Indian plate: Implications for the forces
  driving plate motions. *Journal of Geophysical Research: Solid Earth*, 115, B03410.
  [doi:10.1029/2009JB006634](https://doi.org/10.1029/2009JB006634)
- <a id="ref-davies1995"></a>Davies, J. H., & von Blanckenburg, F. (1995). Slab breakoff: A
  model of lithosphere detachment and its test in the magmatism and deformation of
  collisional orogens. *Earth and Planetary Science Letters*, 129, 85–102.
  [doi:10.1016/0012-821X(94)00237-S](https://doi.org/10.1016/0012-821X(94)00237-S)
- <a id="ref-jordan1978"></a>Jordan, T. H. (1978). Composition and development of the
  continental tectosphere. *Nature*, 274, 544–548.
  [doi:10.1038/274544a0](https://doi.org/10.1038/274544a0)
- <a id="ref-lenardic2003"></a>Lenardic, A., Moresi, L.-N., & Mühlhaus, H. (2003). Longevity
  and stability of cratonic lithosphere: Insights from numerical simulations of coupled
  mantle convection and continental tectonics. *Journal of Geophysical Research: Solid
  Earth*, 108(B6), 2303. [doi:10.1029/2002JB001859](https://doi.org/10.1029/2002JB001859)
- <a id="ref-molnar2009"></a>Molnar, P., & Stock, J. M. (2009). Slowing of India's
  convergence with Eurasia since 20 Ma and its implications for Tibetan mantle dynamics.
  *Tectonics*, 28, TC3001. [doi:10.1029/2008TC002271](https://doi.org/10.1029/2008TC002271)
- <a id="ref-vanhinsbergen2011"></a>van Hinsbergen, D. J. J., Steinberger, B., Doubrovine,
  P. V., & Gassmöller, R. (2011). Acceleration and deceleration of India-Asia convergence
  since the Cretaceous: Roles of mantle plumes and continental collision. *Journal of
  Geophysical Research: Solid Earth*, 116, B06101.
  [doi:10.1029/2010JB008051](https://doi.org/10.1029/2010JB008051)
- <a id="ref-wilson1966"></a>Wilson, J. T. (1966). Did the Atlantic close and then re-open?
  *Nature*, 211, 676–681. [doi:10.1038/211676a0](https://doi.org/10.1038/211676a0)

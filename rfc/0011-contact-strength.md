# RFC 0011: Collisions are held by lithospheric strength

- **Status:** Draft
- **Date:** 2026-10-16
- **Related:** #449 (exploration), #452 (measurement gate; PRs #453, #457); #443 (RFC 0001
  tracking), #396 (RFC 0008 tracking), #402 (lock and fuse, held for this RFC); #314
  (shortening cascade), #340 (craton partition rules), #315 (slab pull), #310 (Hm cap);
  RFCs [0001](0001-merge-cap-craton-loss.md) (merges must not destroy cratons),
  [0008](0008-craton-collision-jamming.md) (cratons jam collisions; this RFC supersedes its
  solution 6), [0009](0009-craton-mosaic.md) (cratons first)

## Problem

When two continents collide, the model lets them push into each other with nothing to stop
them where the rock can't give way. Each plate keeps its full column, so the same ground ends
up under two columns of crust. That double-counted volume builds for tens of Myr, and a plate
merge eventually has to fit both columns under the 84 km Hc cap
(`SUTURE_ACCRETION_MAX_HC_M`). What doesn't fit is clipped. RFC 0001's gate (#448) put 97.5%
of merge craton loss on this, and found 13 of 15 merges preceded by craton-on-craton overlap.
The H4 carry fix (#450) and the per-node crust rules (#451) cut the clipped craton by an
order of magnitude, but the overlap still builds. The merge isn't the cause: it is where the
overlap meets the cap.

**Why it builds: four things in the code.**

1. **Plate motion doesn't respond to whether a contact can deform.** `torque.shift_plate`
   integrates basal drag, slab pull and drag, ridge push and the contact drag, then rotates
   the plate rigidly. The contact drag (`torque.contact_drag_pa_s_per_m`, RFC 0008's #399) is
   η × g × v, a damper proportional to the closing speed. Against a steady drive a damper
   lowers the steady speed but never reaches zero, and nothing grows with how far one plate
   has pushed into the other.
2. **Deformation is vetoed where there is craton.** Boundary retreat never consumes a craton
   cell (`cratons.retreat_allowed`, #340). The shortening cascade (`shortening.py`, #314) gives
   a full craton nine times the yield stress of ordinary crust
   (`CASCADE_CRATON_YIELD_GAIN` = 8), so it passes shortening on rather than thickening. What
   the cascade can't place "stays in the boundary overlap, where boundary retreat and suture
   accretion consume it" (`shortening.py`'s docstring). At craton, retreat doesn't. So the
   model already has a yield law for *where* shortening goes inside a plate. But shortening
   that nothing can take up feeds nothing back to the force balance: the plates keep closing.
3. **The contact force isn't a pair force.** `shift_plate` solves each plate with its
   neighbour's ω held fixed, and each side's coefficient comes from its own nodes. At jammed
   contacts the two sides' drag forces differ by a median 25% (#452).
4. **A merge doesn't recognise a finished collision.** A pair merges after a random
   0–100 Myr timer and a size-dependent roll (`merge_split.update_collision_progress`,
   `_collision_threshold_years`, `_merge_probability`), or after 30 Myr above 30% overlap
   (`FORCED_MERGE_OVERLAP_FRACTION`, `pop_ready_forced_merge`). Neither checks whether the
   plates still move relative to each other.

**Evidence: the #452 gate.** `bin/debug/probe_collision_budget.py` replays a seed and taps
the step without feeding back: every tapped run ends with exactly the craton volume of an
untapped run. Each step, for each contact pair (continental on both sides, or craton within
500 km), it measures:
- the overlap and its double-counted volume D (both columns' Hc over the overlap ground);
- which overlapping cells are **blocked**, by craton (RFC 0008's jam test) or by the Hc cap;
- the **stall force**: the normal force per metre of contact that would hold the pair's mean
  closing at zero. It re-solves both plates' own shift systems jointly, which reproduces the
  model's closing rate to a median 0.002 cm/yr.

The runs used seeds 326494039, 748935664, 580014012 and 111065103, from `generate_world`
with defaults, to 1 Gyr, on `6579830` (after #450 and #451). The 5 Myr step check failed, so
the gate decided on 2 Myr steps. Pooled:

| Measure | 2 Myr steps | 5 Myr steps |
|---|---|---|
| Blocked share of D's net growth, all pairs | 73% | 56% |
| … continental pairs | 84% | 79% |
| … in the 30 Myr before merges | 95% | 87% |
| Of that blocked share, the share blocked by the cap | ≈ 0 | ≈ 0 |
| Stall force at craton-blocked contacts, p10 / p50 / p90 | 0.4 / 2.0 / 7.8 × 10¹² N/m | 0.5 / 2.8 / 9.8 × 10¹² N/m |
| Craton-blocked contact-steps with stall ≤ 5 × 10¹² / 10¹³ / 2 × 10¹³ N/m | 81% / 93% / 98% | 71% / 90% / 97% |
| Craton-blocked contacts' median closing rate | 0.22 cm/yr | 0.26 cm/yr |
| Stall force at deformable contacts, p50 | 4.9 × 10¹² N/m | 6.5 × 10¹² N/m |
| Slab pull's share of that stall force, p50 | 81% | 85% |
| `merge_capped_m3`, four seeds (M km³) | 394 | 166 |

The gate also found:
- **Merges happen while the plates still move apart.** Every merge that lost craton to the
  cap was jammed at the merge, for a median 53 Myr. The pair was still rotating apart at a
  median 0.13°/Myr (about 1.4 cm/yr at the surface).
- **#402's ceiling.** D at the start of each such jam was 29% of D at the merge, so merging
  at jam onset would remove at most about 70% of what those merges reconcile.
- **The push from thickened crust is weak in this model.** The model's mantle lithosphere
  (3,300 kg/m³) is denser than its asthenosphere (3,250 kg/m³), so the push from a contact
  band's own column vanishes once Hm passes about 130–200 km. At Tibet-like bands
  (Hc 65–75 km) it is a median 1.7 × 10¹² N/m at 2 Myr, and −1.6 × 10¹² N/m at 5 Myr, where
  the bands' roots are thicker. The crust-only Airy figure is 8.1 × 10¹² N/m.

[analysis/issue452/README.md](../analysis/issue452/README.md) has the per-seed tables and the
method.

**What the evidence doesn't show.**
- **The stall force is static.** It is computed from one step's systems, with the force
  spread evenly over the contact's edge and its mean closing held at zero. It says how much
  force a contact would need, not how the plates evolve once the force is applied. Stopping
  one contact passes load to the plates' other boundaries.
- **Why craton blocks is still lumped.** The gate's craton class covers two different
  contacts: a cratonic retreat candidate facing weaker crust (an indenter), and craton on
  both sides. Gate G-R below separates them.
- **Whether jammed pairs converge toward matching motion** wasn't measured over time. G-R
  measures it.
- **The numbers depend on step size.** The 5 Myr step coarsens the overlap to about one cell
  and thickens the bands' Hm.

## Proposed solutions

**What happens in nature.**
- **Collisions resist with forces of about 10¹² to 10¹³ N/m, and convergence continues.**
  India–Asia convergence dropped from ~15 to ~4 cm/yr between 50 and 35 Ma. A force balance
  explains that with a resistance from Himalaya–Tibet mountain building of ~5–6 × 10¹² N/m
  [[Copley et al., 2010]](#ref-copley2010). The vertically integrated stress available in
  most tectonic settings has an upper limit near 10 × 10¹² N/m. Ridge push is 2–3 × 10¹² N/m,
  and the integrated stress reaches ~7 × 10¹² N/m near the Tibetan plateau
  [[Porth, 2000]](#ref-porth2000). India keeps moving because mantle flow drags its thick
  cratonic root, with a maximum drag comparable to slab pull (10¹³ N/m). The resisting force
  at the plate boundary is close to the one from Tibet's gravitational potential energy
  [[Li et al., 2024]](#ref-li2024).
- **Cratons survive because they are strong, and the belts around them are weak.** In coupled
  convection and tectonics models, a high brittle yield stress for cratonic lithosphere as a
  whole is the effective and robust way to keep cratons stable. The yield-stress contrast
  needed is smaller when cratons are bordered by weak mobile belts
  [[Lenardic et al., 2003]](#ref-lenardic2003). Elastic thickness, a proxy for integrated
  strength, is high for old cratons and lower for orogenic belts and rifts
  [[Burov & Watts, 2006]](#ref-burov2006).
- **Cratons can yield and thicken.** Depleted mantle lithosphere can grow into a > 200 km
  cratonic root by lateral compression followed by gravitational self-thickening
  [[Wang et al., 2018]](#ref-wang2018). A craton under enough compression is shortened, not
  consumed.

The literature supports a resisting force bounded by the lithosphere's integrated strength,
of order 10¹² to 10¹³ N/m. It doesn't support a hard wall, or a push from thickened crust
alone that stops collisions.

**Solutions.**
1. **Hold blocked contacts with a yield-limited contact force (recommended).** One normal
   force per contact pair, equal and opposite on the two plates, solved jointly with their
   motion. Where the contact is blocked (craton, or the Hc cap), the force holds the blocked
   edge from closing, up to the contact's integrated strength. Past that, the contact yields:
   it keeps closing at the rate the excess drive allows, and solution 4 absorbs that
   convergence without losing volume. This answers the gate's G1 and G3: the overlap builds
   at blocked contacts, and the force to hold them is mostly within lithospheric strength.
   It also turns the cascade's yield law (#314) into a force-balance term. The existing
   contact drag stays as the damper.
2. **A push from thickened crust.** The gravitational potential energy contrast of the
   contact band, in the force balance at every continental contact (#449's original
   proposal). **Not proposed.** The gate found the overlap builds where crust can't thicken,
   so there is no push there. In this model's columns the push is also small or negative
   wherever Hm is thick. It could be added later as part of the resistance at deformable
   contacts, once Hm thickening is physical.
3. **A hard non-penetration rule.** Blocked contacts never close, whatever the drive.
   **Not proposed.** Cratons do deform under enough compression
   [[Wang et al., 2018]](#ref-wang2018), and an unbounded force would freeze any pair whose
   drive exceeds real strength: 2% of craton-blocked contact-steps at 2 × 10¹³ N/m.
4. **Absorb the overlap as it forms, without changing plate motion.** Let retreat consume
   cratonic cells, and carry their craton onto the receiving plate, so no overlap is left to
   double-count. **Not proposed on its own.** It conserves volume, but it consumes craton at
   every collision, where nature holds them. It is the yield branch of solution 1.
5. **Lock jammed contacts and fuse them** (RFC 0008 solution 6, #402): a hazard rate with mean
   time `CRATON_LOCK_MYR`. **Superseded.** It stops the overlap building only by merging; at
   best it removes about 70% of what merges reconcile (#452's ceiling). It merges plates
   still rotating apart, and its rate has no measured basis. Kept as the fallback if this
   RFC's merge criterion is falsified.

**Merging.** Under solution 1, a held contact stops closing, but the plates may still slide or
rotate relative to each other. A pair merges when its contact is held and the two plates'
rotations match: relative rotation |ω_a − ω_b| below a threshold for a sustained time. That
replaces the timer, the size roll and the 30% backstop for continental pairs.

**Recommendation:** gate G-R, then solution 1 with solution 4 as its yield branch, and merging
on matched rotation. This supersedes RFC 0008's solution 6 (#402), and takes over RFC 0001's
solution 3 ("stop craton-on-craton overlap from building"). RFC 0001's solutions 1 and 2,
already disqualified by its gate, are left to that RFC to withdraw.

## Design

**Units.** Forces per unit contact length in N/m; torques in N·m; angular velocities in
rad/yr (`Plate.omega`), with relative rotation reported in °/Myr; closing rates in cm/yr;
lengths in m; times in Myr.

### Gate G-R: what blocks, and do jammed pairs converge?

`probe_collision_budget.py` (extended on this RFC's branch) replays the four seeds at 2 Myr
steps. It splits the craton-blocked class into:
- **craton veto:** a cratonic retreat candidate facing non-cratonic crust, an indenter;
- **craton on both sides.**

It also records each pair's relative rotation |ω_a − ω_b| every step, and summarises jammed
episodes of at least 3 steps by their relative rotation at the start against at the end.

| Outcome | Design consequence |
|---|---|
| Craton on both sides makes most of the blocked growth | The yield force uses the stronger side's craton (both are strong). |
| Craton veto makes most of it | The veto contact should yield at the weak side's strength (Copley's resistance is the weak side's deformation): `F_y` from the weaker side's craton strength, and on yield the weak upper plate takes the shortening. |
| Relative rotation falls by half or more in most jams | The contact force and the existing drag bring rotations together; the merge criterion needs nothing more. |
| Relative rotation doesn't fall in most jams | Add a tangential (sliding) resistance, bounded by friction times the normal force, to the contact solve, so held pairs can come to matched motion. |

**Result:** *pending (runs in progress).*

### Phase 1: one contact force per pair

Restructure the shift so plates in contact are solved together.

- **Assemble every plate's system before rotating any of them.** `torque.shift_plate`
  becomes three steps: assemble each plate's `(lhs, rhs)` as `integrate_omega` builds them,
  solve, then rotate (`apply_omega_and_rotate`, unchanged).
- **One drag per contact.** The contact drag's coefficient is computed once per contact
  from both sides' nodes, and applied with opposite signs. The two plates' rows are coupled
  through the drag's dependence on the other plate's ω (`b = K ω_other`, today held fixed).
  The global system has 3N unknowns for N plates, with a 3 × 3 block per plate and an
  off-diagonal block per contact pair. N is tens to a hundred, so a dense solve is cheap.
- **No behaviour change beyond symmetry.** With no contacts this is today's solve. With
  contacts it is today's solve with the neighbour's ω no longer lagged. The #452 probe already
  solves pairs this way and reproduces the model's closing rates to 0.002 cm/yr.

### Phase 2: blocked contacts are held up to their strength

For each contact pair *c* with blocked edge length `L_b,c` > 0:
- **Blocked edge.** The edge nodes (`collision_polarity.CONTACT_EDGE_SPACINGS`) whose cell is
  blocked, read from the prepass's frozen frame as the gate does:
  - craton: a retreat candidate holding craton (`collision_polarity.retreat_candidates`), or
    craton on both sides;
  - cap: either column at the Hc cap.
- **Contact force.** `f_c ≥ 0` N/m, normal to the contact, pushing the plates apart, spread
  over the blocked edge nodes of both sides. It enters each plate's rhs as a torque, as the
  probe's lever does (`_edge_lever`).
- **Constraint.** The blocked edge's mean normal closing at the end of the step,
  `v_c(f) = v_c,0 + Σ_d W_cd f_d`, satisfies `v_c ≤ 0` while `f_c < F_y,c`. The matrix `W` is
  each contact's closing response to each force, from the phase 1 system.
- **Yield.** `f_c ≤ F_y,c`. A contact at `f_c = F_y,c` is yielding and keeps closing at
  `v_c > 0`.

This is a box-constrained linear complementarity problem in the contact forces. Projected
Gauss–Seidel over the contacts converges in tens of sweeps for tens of contacts, and the
forces then go back into the phase 1 solve. It is solved with the implicit (end-of-step) ω,
as `integrate_omega` treats drag, so a held contact stops closing whatever the step length.

**Yield strength.** `F_y = F_0 + s (F_c − F_0)`, with `s` the contact's craton strength
(`cratons.strength`) per G-R's outcome:
- `F_0` = 5 × 10¹² N/m: the resistance of a deforming continental collision zone, Copley's
  Himalaya–Tibet force.
- `F_c` = 1.5 × 10¹³ N/m: above the ~10¹³ N/m upper limit of available tectonic stress
  [[Porth, 2000]](#ref-porth2000). A full craton yields only under unusual drive, as
  [[Lenardic et al., 2003]](#ref-lenardic2003) require for craton longevity.
- **Scale check against the gate:** 81% of craton-blocked contact-steps need ≤ 5 × 10¹² N/m
  and 98% ≤ 2 × 10¹³ N/m. With these values, fewer than 10% of craton-blocked contacts
  should yield.

**What the force does not do.**
- It doesn't act at deformable contacts. There retreat and the cascade take up the
  convergence, and the existing drag damps it.
- It doesn't act before contact.
- It doesn't replace slab pull's loss at a jam (RFC 0008 solution 1).

### Phase 3: a yielding contact absorbs the overlap without losing volume

At a yielding contact, the blocked edge closes by `v_c × Δt` this step. That convergence is
taken up in the deform pass:
- **Cap-blocked cells:** the existing staged placement and ceiling overflow
  (`quad_tectonics._place_ceiling_overflow`, `_accrete_onto_survivors`).
- **Craton-blocked cells.** `retreat_allowed` admits a cratonic cell on a yielding contact,
  up to the yield's area budget `v_c Δt L_b,c`.
  - The retreated column goes through suture accretion like any other.
  - Its craton is carried with its Hc onto cratonic receivers, as root thickening by lateral
    compression [[Wang et al., 2018]](#ref-wang2018). This is a transfer, not a sink, with
    the H4 rule (#450): craton ≤ Hc on every receiver.
  - Onto non-cratonic receivers it becomes orogenic crust, booked to `collision_reworked_m3`
    as today.
  - The budget is in area per year, so it doesn't depend on step length.

### Phase 4: merge on matched rotation

A continental pair merges when:
- its contact is **held**: the constraint is active and not yielding over at least half its
  blocked length;
- its relative rotation `|ω_a − ω_b|` is below `MERGE_MATCHED_ROTATION_DEG_PER_MYR`;
- and both have held continuously for `MERGE_HOLD_MYR`, accumulated in years
  (`World.contact_hold_years`, new).

**What it replaces.** The criterion replaces `_collision_threshold_years`' timer and
`_merge_probability`'s roll for continental pairs. The overlap backstop
(`pop_ready_forced_merge`) stays as a logged safety net and should fire rarely (acceptance
below). `merge_omega` already conserves angular momentum, so a matched pair merges with
little change in motion.

`MERGE_HOLD_MYR` starts at 20 Myr, the value RFC 0008 chose for `CRATON_LOCK_MYR`. It is a
design choice: India and Asia still converge 50 Myr after contact, so the literature gives no
single suturing time.

### Timestep robustness

- **Phase 2:** the constraint is on the end-of-step closing rate, so a held contact gains no
  overlap at any step length. The lattice still quantizes the overlap to about a cell.
- **Phase 3:** yield absorbs an area set by rate × years.
- **Phase 4:** the hold counts years.

Acceptance compares 2 and 5 Myr runs of the same seeds.

### Conservation and ownership

- **Momentum.** The contact force is internal to each pair, equal and opposite, so it adds
  no net torque to the two plates together.
- **Volume.** Phase 3 moves volume and books what it can't place, as suture accretion does
  today. Craton carried onto cratonic receivers is a transfer between cells, not a ledger
  sink, so the craton ledger balances without a new account.
- **One owner per resistance.** At blocked contacts the contact force owns the resistance to
  closing, and the existing drag is kept for its damping. `CONTACT_RESISTANCE_N_PER_M`
  (5 × 10¹² at 4 cm/yr) and `F_y` are separate terms: the drag goes to zero as closing
  stops, while the constraint force holds.

### Cost and saves

- **Runtime.** Phase 1's global solve is a dense 3N × 3N system once per step, about 300 ×
  300 at most. Phase 2's complementarity solve is tens of sweeps over tens of contacts. The
  boundary searches already happen in the shift and the prepass, and can share
  `torque.BoundarySearchCache`. Expected cost: well under 5% of a step; acceptance sets a
  10% budget.
- **Saves.**
  - `World.contact_hold_years` (new, a dict of plate pair to years) defaults to empty, so old
    saves load with every pair unheld.
  - `collision_progress` and `overlap_progress` stay in the save format. Phase 4 stops
    reading `collision_progress` for continental pairs.

### Knobs

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `CONTACT_YIELD_N_PER_M` (`F_0`) | 5 × 10¹² N/m | [[Copley et al., 2010]](#ref-copley2010) | Literature-constrained; sensitivity at 3 × 10¹² and 10¹³ on the four seeds, 2 Myr. |
| `CRATON_CONTACT_YIELD_N_PER_M` (`F_c`) | 1.5 × 10¹³ N/m | Above [[Porth, 2000]](#ref-porth2000)'s 10¹³ N/m limit; design choice | Sensitivity at 10¹³ and 3 × 10¹³. Metric: share of craton-blocked contacts yielding, mean plate speed. |
| `MERGE_MATCHED_ROTATION_DEG_PER_MYR` | 0.05 °/Myr (~0.5 cm/yr at the surface) | Design choice, below the gate's 0.13 °/Myr at capped merges | Calibrated on the four seeds: merge count within the acceptance band. |
| `MERGE_HOLD_MYR` | 20 Myr | Design choice (RFC 0008's `CRATON_LOCK_MYR`) | Sensitivity at 10 and 40 Myr. |
| Tangential friction (if G-R calls for it) | coefficient 0.1 | Design choice | Only if G-R says so; tuned on relative rotation at merges. |

## Detailed design (provisional touchpoints)

- **`backend/app/torque.py`.**
  - Split `shift_plate` into assemble (new `assemble_shift_system`), solve and rotate
    (`apply_omega_and_rotate`).
  - New `solve_contacts`, holding phase 1's joint solve and phase 2's complementarity solve.
  - `collision_drag_coefficients` and `contact_drag_pa_s_per_m` move to one coefficient per
    pair.
  - New constants `CONTACT_YIELD_N_PER_M` and `CRATON_CONTACT_YIELD_N_PER_M`.
- **`backend/app/world.py`.**
  - `step_world_progress`'s shift loop calls the assemble → solve → rotate sequence. The
    progress units still count per plate.
  - New `World.contact_hold_years`.
- **`backend/app/collision_polarity.py`.** The blocked-edge classification reuses
  `retreat_candidates` and `PairContact`. It is read from the last frozen frame, as
  `torque.forced_down_mask` already does.
- **`backend/app/cratons.py`** and **`backend/app/quad_tectonics.py`.**
  - `retreat_allowed` takes the yielding contacts' cells (phase 3).
  - `_retreat` honours the yield's area budget.
  - `_accrete_onto_survivors` carries craton onto cratonic receivers.
- **`backend/app/merge_split.py`.**
  - `update_collision_progress` (or a new `update_contact_holds`) runs the phase 4 criterion.
  - `_collision_threshold_years` and `_merge_probability` retire for continental pairs.
  - `pop_ready_forced_merge` stays as a logged backstop.
- **`backend/app/shortening.py`.** Unchanged. Its `returned_unrouted_m2` at held contacts
  should fall, since the held contact stops demanding shortening; its `phase_budget` booking
  is the check.
- **Tests:**
  - `backend/unit_tests/test_torque.py`: pair-force symmetry; a blocked pair held under a
    drive below `F_y`; yield above it; step-length independence of a held contact.
  - `test_merge_split.py`: the merge criterion.
  - `test_quad_tectonics.py` and `test_cratons.py`: yield absorption conserves Hc and craton.
  - `test_persistence.py`: old saves load with `contact_hold_years` empty.
- **Probe and docs.**
  - `bin/debug/probe_collision_budget.py` measures acceptance.
  - `docs/simulation-model/tectonics.md` (plate motion, merge and split) and
    `docs/architecture.md` (step flow) get the changes.

## Assumptions and uncertainty

- **Blocked cells stand in for "can't deform".** The probe's classes are rules: retreat veto,
  craton on both sides, the cap. They mark where today's model can't take up convergence,
  not a measured strength. Phase 2 uses the same classes, so it is consistent with the gate,
  but not more physical than the rules it reads.
- **The strengths are weakly constrained.** `F_0` rests on one well-studied collision.
  `F_c` is a design choice above the available-stress limit. Integrated strength varies with
  geotherm, composition and strain rate. The sensitivity runs bracket both.
- **The force is spread evenly and acts on the edge's mean closing.** A long contact with a
  short blocked segment is held where the segment is. The rest of the contact still closes
  where it can deform. Whether a mean constraint stops the overlap at every blocked node is
  checked by D's growth at held contacts (acceptance).
- **Slab pull supplies most of the drive.** Holding a contact puts the load on the plates'
  other boundaries and on basal drag. Plates may slow, splits may change (26% of splits were
  of jammed plates at 2 Myr), and new contacts may jam. The acceptance bounds speed and
  plate count; the falsification criteria say what a collapse would mean.
- **Matching rotation may not be reachable without tangential resistance.** G-R decides
  before phase 4 is built.
- **Step size.** The gate's numbers move with step size. Acceptance runs at 2 Myr, with 5 Myr
  checked for the same verdicts.

## Acceptance criteria

All runs use the four seeds 326494039, 748935664, 580014012 and 111065103, `generate_world`
defaults, 1 Gyr, at 2 Myr steps unless stated, measured with `probe_collision_budget.py`
against pinned `6579830`.

- [ ] **Gate G-R** run and recorded in this RFC, with the phase 2 strength rule and the phase 4
  tangential term chosen from it.
- [ ] **Phase 1.** The two sides' contact forces differ by < 1% at every contact. Mean plate
  speed and plate count within 10% of baseline. Closing rates at deformable contacts within
  10%.
- [ ] **Phase 2.**
  - D's net growth at blocked contacts down ≥ 80% against baseline, pooled.
  - Craton-blocked contacts closing at a median ≤ 0.02 cm/yr (baseline 0.22).
  - At most 15% of craton-blocked contact-steps yielding.
  - Mean plate speed down no more than 25%, and no plate held below 20% of its pre-contact
    speed for more than 50 Myr, unless it is fully surrounded by held contacts.
- [ ] **Phase 3.** The craton and continental-material ledgers close to rounding. Hc volume
  is conserved through yield absorption, except bookings to `no_outlet_subducted_m3`. Craton
  ≤ Hc on every cell.
- [ ] **Phase 4.**
  - Every continental merge has `|ω_a − ω_b|` ≤ `MERGE_MATCHED_ROTATION_DEG_PER_MYR` at the
    merge, except backstop merges, which number < 10% of merges.
  - D at merge down ≥ 75% against baseline.
  - `merge_capped_m3` down ≥ 90% against baseline (394 M km³).
  - Merge count within ±30% of baseline.
- [ ] **Step size.** At 5 Myr, the same seeds meet the phase 2 and 4 criteria, and the blocked
  share of D's net growth is within 10 points of 2 Myr.
- [ ] **Cost and compatibility.** Step runtime +10% at most. Old saves load and step.
  `ruff`, `quality_check.sh` and the affected unit tests pass.

## Falsification criteria

- **Holding blocked contacts stalls the planet.** Mean plate speed falls by more than half, or
  most plates end below 20% of their pre-contact speed, with `F_y` inside the literature
  range. Then the model's drive, not the contact resistance, is wrong: slab pull supplies
  ~80% of the stall force (#315). Revise toward the driving forces, and treat #402 as the
  practical fix meanwhile.
- **Overlap still builds at held contacts.** D grows at contacts whose mean blocked closing
  is ≤ 0, so the constraint misses where the overlap comes from: rotation about a nearby pole,
  or a short blocked segment on a long contact. Revise the constraint to act per segment
  before dropping the idea.
- **Most blocked contacts yield at literature strengths.** More than half of craton-blocked
  contact-steps sit at `F_y` with `F_c` ≤ 3 × 10¹³ N/m. Then this model's drive at contacts
  is above anything lithosphere can resist, and solution 4 (absorb the overlap) becomes the
  main path, not the yield branch.
- **Held pairs never come to matched rotation,** even with G-R's tangential term. Then the
  phase 4 criterion is unreachable: fall back to RFC 0008's #402 (lock and fuse) for merging,
  and keep phases 1–3.

## References

- <a id="ref-burov2006"></a>Burov, E. B., & Watts, A. B. (2006). The long-term strength of
  continental lithosphere: "jelly sandwich" or "crème brûlée"? *GSA Today*, 16(1), 4–10.
  [doi:10.1130/1052-5173(2006)016<4:TLTSOC>2.0.CO;2](https://doi.org/10.1130/1052-5173(2006)016%3C4:TLTSOC%3E2.0.CO;2)
- <a id="ref-copley2010"></a>Copley, A., Avouac, J.-P., & Royer, J.-Y. (2010). India-Asia
  collision and the Cenozoic slowdown of the Indian plate: Implications for the forces
  driving plate motions. *Journal of Geophysical Research: Solid Earth*, 115, B03410.
  [doi:10.1029/2009JB006634](https://doi.org/10.1029/2009JB006634)
- <a id="ref-lenardic2003"></a>Lenardic, A., Moresi, L.-N., & Mühlhaus, H. (2003). Longevity
  and stability of cratonic lithosphere: Insights from numerical simulations of coupled
  mantle convection and continental tectonics. *Journal of Geophysical Research: Solid
  Earth*, 108(B6), 2303. [doi:10.1029/2002JB001859](https://doi.org/10.1029/2002JB001859)
- <a id="ref-li2024"></a>Li, Y., Liu, L., Li, S., Peng, D., Cao, Z., & Li, X. (2024). Cenozoic
  India-Asia collision driven by mantle dragging the cratonic root. *Nature Communications*,
  15, 6674. [doi:10.1038/s41467-024-51107-0](https://doi.org/10.1038/s41467-024-51107-0)
- <a id="ref-porth2000"></a>Porth, R. (2000). A strain-rate-dependent force model of
  lithospheric strength. *Geophysical Journal International*, 141(3), 647–660.
  [doi:10.1046/j.1365-246x.2000.00115.x](https://doi.org/10.1046/j.1365-246x.2000.00115.x)
- <a id="ref-wang2018"></a>Wang, H., van Hunen, J., & Pearson, D. G. (2018). Making Archean
  cratonic roots by lateral compression: A two-stage thickening and stabilization model.
  *Tectonophysics*, 746, 562–571.
  [doi:10.1016/j.tecto.2016.12.001](https://doi.org/10.1016/j.tecto.2016.12.001)

# RFC 0021: Plates are driven by their slabs: correct and calibrate the plate forces and the rifting knobs together

- **Status:** Draft
- **Date:** 2026-11-17
- **Related:** #942 (diagnosis, [`analysis/issue942/`](../analysis/issue942/README.md)), #888 (RFC 0018 D1,
  where the slow plates were found), #797 (RFC 0018 tracking), #396 (RFC 0008 tracking), #399, #423, #760 (the
  contact brake and its tuners), #468 (boundary-load normalization), #783 (plate speed in the Stats panel);
  #846 and PR #950 (RFC 0018's rift-line method in the model, behind `World.split_method`); RFCs
  [0008](0008-craton-collision-jamming.md) (contact brake), [0012](0012-physics-first-simulation.md)
  (quasi-static motion solve), [0018](0018-plates-split-along-weak-open-paths.md) (rifting; its tunable
  knobs are in scope here)

## Problem

Plates move several times slower than Earth's, and they are moved by the wrong force.

**Speed.** On the 20 RFC 0018 development worlds, brake on, plates-only stepping, 200 Myr at 1 Myr steps, the
area-weighted crust speed is **0.83 cm/yr** (continental 0.75, oceanic 0.92). That is 4.8–6.0 times below the
working Earth figure of about 4–5 cm/yr, which G0 checks against a published plate model. 70% of plate area moves
under 1 cm/yr, the fastest plate-sample in 9,442 is 4.3 cm/yr, and worlds range from 0.26 to 1.90 cm/yr. Speed is
flat over the 200 Myr ([`analysis/issue942/`](../analysis/issue942/README.md), "Speed census").

**What sets it.** A per-plate breakdown of the quasi-static motion solve (`torque.solve_plates_quasistatic`) takes
each torque's component along the plate's solved rotation axis (#942):

| Pooled share | 20 worlds, brake on | 2 worlds, brake on | 2 worlds, brake off |
|---|---:|---:|---:|
| **Driving:** prescribed mantle flow (`mantle.flow_at` through basal coupling) | 74% | 80% | 56% |
| slab pull | 20% | 15% | 42% |
| ridge push | 0.1% | 0.1% | 0.6% |
| **Resisting:** contact brake (RFC 0008) | 58% | 59% | 0 |
| basal drag | 26% | 26% | 44% |
| sunk-slab drag | 9% | 7% | 50% |
| Crust speed (cm/yr) | 0.83 | 1.05 | 2.31 |

Four findings follow.

1. **The prescribed flow is the driver and the ceiling.** Mantle flow alone would carry plates at 2.78 cm/yr on
   the 20 worlds; a world's speed follows that figure with a rank correlation of 0.85. The flow comes from 8
   random, fixed convection centres of 2–6 cm/yr (`mantle.MANTLE_FLOW_REFERENCE_RATE` = 4 cm/yr), a constant set
   on 2026-08-15 for the first motion model, which fitted each plate's rotation directly to the flow. It predates
   the torque engine, sunk-slab drag and the brake.
2. **The contact brake removes 55–60% of plate speed.** On #942's two worlds speed goes from 1.05 to 2.31 cm/yr
   with it off; on #888's five worlds from 0.63 / 0.77 to 1.55 / 2.12 cm/yr (continental / oceanic). Its craton
   component alone changes nothing measurable.
3. **Slabs are a net brake.** With the brake off, sunk-slab drag dissipates more work than slab pull does, pooled
   −22% of the pull's work. Of the drag's dissipation, 40% is on trench-parallel motion that pull cannot offset;
   34% of the subducting length is dragged but not pulled, because its column is buoyant; and on a quarter of the
   length the trench-normal speed is already past the pull-to-drag terminal speed. Mean pull is 2.9×10¹² N/m,
   59% of the 4.9×10¹² N/m for the 40 km oceanic reference column.
4. **Ridge push is missing.** `ridge_push_force_vectors_n` uses `½ g (ρ_a − ρ_w) e²`, the water-load term alone,
   with the ridge 711 m above the abyssal plain: 5.6×10⁹ N/m. Copley, Avouac & Royer estimate about 2.5×10¹² N/m
   for the Indian plate [[Copley et al., 2010]](#ref-copley2010). The missing part is the lithosphere's own
   density deficit at the ridge.

Behind findings 3 and 4 is one column. Oceanic mantle lithosphere relaxes toward
`lithosphere.REFERENCE_HM_OCEANIC_M` = 40 km whatever its age (`rheology.relax_oceanic_mantle_lithosphere`); on
one development world at 60 Myr the median oceanic cell is 102 Myr old and 39.9 km thick (#888). So old ocean
floor is no denser, thicker or higher-standing than young, and neither slab pull nor ridge push can grow with it.

**On Earth** slabs drive plates. Forsyth & Uyeda found that the forces on the downgoing slab "control the velocity
of the oceanic plates" and are "an order of magnitude stronger than any other force", with slab-attached plates at
a terminal velocity where the slab's pull is "nearly balanced with the resistance acting on the slab"
[[Forsyth & Uyeda, 1975]](#ref-forsyth1975). Slab pull and slab suction account for about 60% and 40% of the forces
on plates today, and subducting plates move about 4 times faster than non-subducting ones
[[Conrad & Lithgow-Bertelloni, 2004]](#ref-conrad2004). In the model, plates with a boundary past the subduction
gate are 1.8 times faster than the rest, but about 1.3 of that factor is faster flow under them, and the share of a
plate's boundary that is down-going and closing has no rank correlation with its speed (−0.05).

**Why it matters beyond speed.** Slow plates build little pull across continents, which holds back rifting and
spreading: with the brake off #888's rift-line prototype split 10.9 plates per world per 100 Myr against 4.4 with
it on. The rift load in RFC 0018 is built from the same per-node forces (`torque.gather_per_node_force_field`), so
any change to the forces moves RFC 0018's calibration, and the two have to be tuned together.

**Principles.** The current behaviour violates principle 8, *let each process follow its own physics*: plate
motion inherits the pattern and speed of a prescribed flow field, not the forces of the plates' own slabs, ridges
and contacts. A retune of the flow rate or the brake to restore the mean speed would also violate principle 9:
*do not tune a parameter to compensate* for a missing or misformed force. The invariants that will show the fix
are the force shares above, the subducting/non-subducting speed ratio, and the census speed distribution, with the
breakdown's terms summing to zero at every sample.

## Proposed solutions

### What happens in nature

- **Slab pull is the largest driver, and it is resisted mainly on the slab.** Slab forces set the velocity of
  slab-attached plates, which reach a terminal velocity [[Forsyth & Uyeda, 1975]](#ref-forsyth1975). Observed
  motions are best predicted if upper-mantle slabs are attached and their pull is about half of the total driving
  force [[Conrad & Lithgow-Bertelloni, 2002]](#ref-conrad2002); predicted and observed motions are reconciled most
  easily if at least half of the excess upper-mantle slab weight contributes to pull [[van Summeren et al.,
  2012]](#ref-vansummeren2012).
- **Little of a slab's weight is spent shearing the surrounding mantle.** In laboratory subduction the net pull on
  the trailing plate is up to 8–12% of the slab's negative buoyancy, up to twice ridge push; the rest drives
  rollback-induced flow (~70%), bends the plate at the trench (~15–30%) and overcomes slab–mantle shear (0–8%)
  [[Schellart, 2004]](#ref-schellart2004). Net slab pull has an upper limit of about twice ridge push on the
  Nazca plate [[Richardson & Cox, 1984]](#ref-richardson1984).
- **Ridge push is of order 10¹² N/m.** About 2.5×10¹² N/m for India, with a net subduction force of similar size
  and a collision resistance of 5–6×10¹² N/m [[Copley et al., 2010]](#ref-copley2010).
- **Basal traction is under about 1 MPa and mostly resists.** Copley et al. need basal tractions linear in plate
  velocity with an asthenospheric viscosity of about 2–5×10¹⁹ Pa·s over 100–150 km; Bird, Liu & Rucker find mean
  basal shear tractions of no more than 1 MPa under the six largest slabless plates [[Bird et al.,
  2008]](#ref-bird2008). Flow in the low-viscosity zone may be driven by the plates (Couette) or by pressure from
  within the mantle (Poiseuille) [[Richards & Lenardic, 2018]](#ref-richards2018), so basal traction can drive
  locally, but on Earth the flow under plates is itself largely set up by slabs ("slab suction").
- **Collision slows plates by a few times, not to a halt.** India slowed from about 15 to about 4 cm/yr of
  convergence as Tibet thickened [[Copley et al., 2010]](#ref-copley2010), about 27%.

### Alternatives

**A. Retune the existing constants only.** Search today's constants (the flow rate, the brake, the viscosities,
`L_sink`, the gate threshold) and the rifting knobs with no force formula changed. Cheap and immediately reversible.
The breakdown predicts that it can reach T1's mean speed but not the force balance, because each defect is in the
form of a term rather than its value:

- Ridge push is `½ g (ρ_a − ρ_w) e²`. Its only inputs are densities and the ridge's height: reaching 10¹² N/m
  needs a ridge about 9.5 km high, and Earth's own 3 km gives only 10¹¹, because the form leaves out the
  lithosphere's thickening.
- Sunk-slab drag resists all tangential motion. Lowering `SLAB_MANTLE_VISCOSITY_PA_S` cuts the trench-normal drag
  that sets a slab's terminal speed by the same factor as the trench-parallel drag that should not be there.
- Pull comes from an Hm that relaxes to 40 km at any age, so `L_sink` and `f_pull` scale every slab alike and
  cannot make old slabs pull harder.
- Loosening the gate adds slabs that, under today's drag, brake: in #942's brake-off run the gated share went
  from 11% (20 worlds, brake on) to 66% (two worlds, brake off), and slabs still netted −22% of their pull's work.
- That leaves the flow rate. About 8 cm/yr would probably give Earth's mean, with plates carried by the random
  centres: basal traction at about Bird's 1 MPa bound, world speed following the flow draw (T7), and slab-attached
  plates no faster (T2, T3).

This is an inference from one breakdown, not a measured result, so A is not rejected on argument: G2 runs it as a
**control arm** with every correction off. If the control arm meets as many targets as the corrected search, A
wins (see falsification).

**B. Correct the force forms the diagnosis found, then calibrate forces, brake and rifting knobs together
(recommended).** Four form corrections, each separately switchable so a gate can adopt or reject it:

1. Ridge push from the columns' gravitational potential energy, not the water load alone.
2. Sunk-slab drag on the slab's own sinking motion only: the trench-normal part, on the length that is pulled.
3. Slab pull and ridge push from an age-based "force column", a cooling-plate thickness from each cell's age,
   leaving the plate's own Hm, isostasy and bathymetry untouched.
4. A subduction test that does not depend on the previous step already closing faster than 1 cm/yr.

Then one joint calibration of the force constants, the RFC 0008 brake and RFC 0018's rifting knobs on development
worlds, judged on held-out worlds. It stays within the current quasi-static solve and the current cost, and every
change is a force the diagnosis measured. It still prescribes the flow pattern.

**C. Replace the convection centres with flow driven by the slabs** (slab suction from subduction history, as in
[[Lithgow-Bertelloni & Richards, 1998]](#ref-lithgowbertelloni1998)). The most physical: the flow under plates
would follow the trenches, and the 0.85 dependence of world speed on a random draw would go. It needs a mantle-flow
solve and a slab history the model does not keep. Deferred to a later RFC; B makes the flow's weight a knob, so
C can follow without redoing B.

**D. Make oceanic Hm itself follow age.** The physical form of B.3: the column thickens as it cools, so
bathymetry, ocean volume, land fraction and strength all change. It moves every surface ledger and the RFC 0018
strength map at once. Out of scope; B.3's force column is the cheap stand-in, and G1 reports what D would change.

**Recommendation: B.** It addresses each finding with the force that is missing or misformed, keeps the motion
solve, and puts every knob that sets speed and rifting through one calibration, so they do not drift against each
other. A and C stay available: A as G2's control arm, C as follow-up work.

## Design

### Force corrections (each behind a switch, adopted or rejected at G1)

**B.1 Ridge push from potential energy.** For a divergent node, the push per metre of ridge is the difference in
vertically integrated lithostatic pressure between the node's column and the reference column it will become,
down to a compensation depth below both:

    F_rp = ∫₀^{d_c} [p_node(z) − p_ref(z)] dz        (N/m)

Each column is water to its isostatic elevation (`lithosphere.isostatic_elevation`), then crust, mantle
lithosphere and asthenosphere at the model's densities. The force acts along the existing push direction
(`−direction_to_neighbor`) with the existing shared front length (#468). Scale check, hand calculation from the
model's constants: with today's columns (ridge Hm 8 km, reference 40 km) it gives 5.3×10¹¹ N/m, 94 times the
current term; with B.3's force column at 50–100 Myr it gives 2.0–2.1×10¹² N/m, close to Copley's 2.5×10¹² N/m.

**B.2 Sunk-slab drag on sinking only.** `slab_drag_coefficient_matrix` uses `c_i (I − p_i p_iᵀ)`, resisting the
node's whole tangential velocity. The slab resists its own descent, which a surface node feels as motion toward
the trench. The corrected form projects on the trench normal `n_i` (`direction_to_neighbor`):

    K_slab = Σ_i c_i (p_i × n_i)(p_i × n_i)ᵀ,   over nodes with positive pull

so the drag acts only on the pulled length and only on trench-normal motion. Trench-parallel motion and unpulled
buoyant length are left to basal drag, which already covers every node. `K_slab` stays symmetric positive
semidefinite, so the implicit solve is unchanged. In #942's brake-off budget the trench-parallel part alone was
46% of the pull's work; without it the pooled slab work would have netted about +24% in place of −22%.

**B.3 Age-based force column.** For slab pull and ridge push only, an oceanic cell's mantle lithosphere is taken
as `H_f = min(h₁ √(age / 1 Myr), H_max) − Hc`, with `age` from `node_created_years`, `h₁` = 13 km and
`H_max` = 95 km. Continental cells keep their Hm. Pull per node is then `max(B(H_f), 0) × L_sink × ds` as now
(`column_buoyancy_pa`). Scale check, from the model's constants: 54, 81 and 84 km at 25, 50 and 100 Myr give 6.6,
9.9 and 10.3×10¹² N/m of pull, against 4.9 at 40 km today. The plate's own Hm and its elevation are not changed;
whether to make them follow age is alternative D.

This is the same cooling law, with the same two values, that #950's rift-line method uses for ocean-floor strength
(`oceanic_age_strength`, `rift_strength.oceanic_strength_hm`, with `ocean_cooling_m_per_sqrt_myr` and
`ocean_plate_thickness_m` in `RiftLineParams`). B.3 must share one law and one pair of values with it, so the column
that resists rifting is the column that pulls and pushes. The implementation issue decides whether the law moves to
`lithosphere.py` with both callers reading it, or B.3 reads it from the world's rift parameters; either way a
setting cannot give the two different values.

**B.4 Subduction test.** Today a node pulls only if it is down-going crust within reach of a neighbour and the
previous step's motion closes it faster than 1 cm/yr (`subducting_boundary_mask`, `TRANSFORM_RATE_THRESHOLD`).
Only 11% of the closing down-going length passes, and passing depends on speed already gained. G1 compares two
replacements: (a) the same test against a ramp, pull scaled by `clip(closing / v_gate, 0, 1)` from zero closing;
(b) a geometric test, down-going crust whose neighbour is not diverging from it and whose collision polarity
(`collision_polarity`) puts it under. Whichever is adopted, the pull and drag masks stay the same set, as now.

### Mantle flow

The prescribed flow stays, with its strength a knob: `MANTLE_FLOW_REFERENCE_RATE` and the number of centres. With
B.1–B.3, the flow is meant to fall to a secondary driver; its target share is set in the calibration targets below,
not by a fixed value here.

### Knobs

Development worlds only for every search; held-out worlds only for acceptance.

| Knob | Starting value | Source | Tuning |
|---|---|---|---|
| `MANTLE_FLOW_REFERENCE_RATE` | 4 cm/yr | Design choice (2026-08-15, v1 motion) | G2 over 1–8 cm/yr |
| `ASTHENOSPHERE_VISCOSITY_PA_S` | 5×10¹⁹ Pa·s | [[Copley et al., 2010]](#ref-copley2010): 2–5×10¹⁹ over 100–150 km | G2 over 2×10¹⁹–1×10²⁰ |
| `SHEAR_ZONE_THICKNESS_M` | 100 km | Copley et al.: 100–150 km | Fixed; only μ/d_s enters |
| `SLAB_MANTLE_VISCOSITY_PA_S` | 10²¹ Pa·s | Design choice | G2 over 10²⁰–10²¹ |
| `SUBDUCTION_LSINK_M` | 250 km | Design choice | G2 over 150–660 km |
| `f_pull` (new): share of slab buoyancy applied as pull | 1.0 | [[van Summeren et al., 2012]](#ref-vansummeren2012): about 20% (anchored roots) to at least half | G2 over 0.2–1.0 |
| `h₁`, `H_max` (B.3; #950's `ocean_cooling_m_per_sqrt_myr`, `ocean_plate_thickness_m`) | 13 km, 95 km | #888's ocean-strength law; a cooling-plate design choice | Fixed, one pair for force and strength; sensitivity at 11 km and 125 km |
| `v_gate` (new, B.4a) | 1 cm/yr | Current threshold | G1 picks the test; G2 over 0.25–2 cm/yr if (a) |
| `contact_drag_multiplier` | 1.0 | RFC 0008, #760 | G2 over 0–1 |
| `CONTACT_RESISTANCE_N_PER_M` | 5×10¹² N/m | [[Copley et al., 2010]](#ref-copley2010): 5–6×10¹² | Fixed; the multiplier scales it |
| `craton_drag_multiplier` | 1.0 | #760 | Fixed; no measurable effect (#888) |
| `COLLISION_RELIEF_FRICTION_SCALE_M` | 4 km | RFC 0008 | Fixed; sensitivity at 2 and 8 km |
| RFC 0018's searched rift knobs, in `World.rift_line_params` (#950): `eps_dot0`, `s_warm`, `k_cap`, `n_plume`, `lam_grow`, `l_tip_m` | The setting RFC 0018's D1 re-run freezes on the #950 code | RFC 0018 | G2, over D1's ranges, recentred on the re-run's setting |
| RFC 0018's strength-map knobs: `mark_power`, `craton_root_m` | The re-run's values (#950 defaults: 1, 0) | #888's strength-map arms | G2, over the ranges those arms tried |
| `route_mark_power` | The re-run's value (#950 default: 1) | #888 | G2 only if the re-run keeps it away from 1; route-only marks gave 90 Myr median rifts in #888 |
| Every other `RiftLineParams` field (`a_min`, `c_dir`, `c_corr`, `f_break`, healing, sub-step limits) | #950's values | RFC 0018 | Fixed: set by RFC 0018 or numerical. `a_min` (15%) is RFC 0018's guard against microplates and stays its decision |

`f_pull`, the viscosities and `L_sink` trade off against each other (the terminal speed goes as
`f_pull B / (μ_slab/d_s)`), so G2 reports the ridge of equally good settings, not one point.

The rift knobs change how many plates there are and how large, not the forces that drive them; they reach speed
only through plate size and the ridges new rifts open. They are tuned here because their effect depends on the
forces: the same `k_cap` splits about 2.5 times as often with the brake off as with it on (#888), and B.1–B.3
raise the pull across plates further.

### Calibration targets

Pooled over development worlds, area-weighted, brake as calibrated:

| # | Target | Value | Source |
|---|---|---|---|
| T1 | Mean crust speed | 3–6 cm/yr | Working Earth figure 4–5 cm/yr, checked in G0 |
| T2 | Speed ratio, slab-attached to other plates | at least 2.5 | Earth about 4 [[Conrad & Lithgow-Bertelloni, 2004]](#ref-conrad2004); 3–4 [[van Summeren et al., 2012]](#ref-vansummeren2012) |
| T3 | Slab pull share of the driving torque | at least 50% | [[Conrad & Lithgow-Bertelloni, 2002]](#ref-conrad2002); [[Conrad & Lithgow-Bertelloni, 2004]](#ref-conrad2004) |
| T4 | Mean basal traction under plates without slabs | at most 1 MPa | [[Bird et al., 2008]](#ref-bird2008) |
| T5 | Ridge push at ridges older than 50 Myr on the flank | 1–4×10¹² N/m | [[Copley et al., 2010]](#ref-copley2010) |
| T6 | Plate speed after continental collision | at least 20% of pre-contact | RFC 0008's criterion, from Copley's ~27% |
| T7 | Spread of world speeds, fastest to slowest | at most 3 (today 7.3) | Design choice: a world's speed should not be set by the random flow draw |
| T8 | RFC 0008's collision measures | not worse than `origin/main` (see acceptance) | RFC 0008 |
| T9 | RFC 0018's realism targets | as many as RFC 0018 defines, scored by its own rules | RFC 0018 |

T2 needs a definition that does not depend on motion, because #942 found the current one circular: a plate is
"slab-attached" when at least a set share of its boundary length is down-going crust under B.4's test at the
reference geometry, with the share fixed in G0 before any G2 run.

### Diagnosis gates

Development worlds are RFC 0018's 20 (#888); every gate runs plates-only (`simulate_climate_biomes = False`) at
1 Myr steps for 200 Myr unless stated, and reports at 0.5 Myr on three of them. From G0 on, every run splits plates
with `World.split_method = "rift_lines"` (#950) at the rifting baseline below, unless stated.

1. **G0: targets and harness.**
   - Take the rifting baseline from RFC 0018's D1 re-run on the #950 code (#888), and record its setting, its
     brake setting (RFC 0018 runs with `contact_drag_multiplier = 0` for now) and the commit it ran on. D1's
     earlier tables are not a baseline: #950's reviews changed plume heat at the split test, flank transfer under
     the crust cap and capped donors, so the production method does not reproduce them. If the re-run has not
     frozen a setting when G0 starts, G0 waits for it rather than choosing one.
   - Check T1's Earth figure against a published plate model (area-weighted mean of `|ω × r|` from published
     Euler poles and plate areas), and record the result. If it is outside 4–5 cm/yr, T1's band moves with it,
     with the owner's agreement.
   - Fix T2's slab-attached definition and threshold on development worlds' generation snapshots, before any
     probe.
   - Freeze a held-out set for this RFC by the RFC 0018 scheme (SHA-256 of
     `mantle-bloom/rfc0021/v1/heldout/{index:02d}`, first four bytes big-endian masked with `0x7fffffff`,
     excluding every seed already in `analysis/issue797/g2-seed-preregistration.json` and the development set),
     20 seeds, committed before any G2 run.
   - Extend #942's brake-off breakdown from 2 to the 20 development worlds, on the rifting baseline, so G1's
     baseline is not two worlds and not the great-circle split. Measure the rift-line method's step cost there.
   - Build the scorer: `plate_torque_breakdown.py`'s measures plus RFC 0008's collision measures plus RFC 0018's
     realism score, one JSON per setting.
2. **G1: one correction at a time.** Each of B.1, B.2, B.3 and B.4 (a and b) alone, then all adopted ones
   together, brake on and brake off, on five development worlds, as causal probes with predictions posted before
   each run. A correction is adopted if it moves its own measure the predicted way (B.1: ridge push into T5's
   band; B.2: slab work positive; B.3: pull per metre up; B.4: gated share of closing down-going length up) and
   breaks no invariant (the breakdown sums to zero; the solver converges; ledgers close). G1 reports each one's
   effect on speed, T2 and the split rate, and what alternative D would change in bathymetry and land fraction.
   If no combination raises the T2 ratio above 1.5 at the current constants, the premise is in doubt (see
   falsification).
3. **G2: joint calibration.** A scrambled Latin hypercube over the G2 knobs above, with the adopted corrections.
   Stage 1: 24 settings on 5 development worlds; stage 2: the best 3 on all 20, scored on T1–T9 (targets met,
   ties broken by total violation as in #888). Reports the trade-off ridge between `f_pull`, viscosities and
   `L_sink`.
   - **Control arm (alternative A).** The same design, run the same way, over the knobs that exist today
     (everything in the table except `f_pull` and `v_gate`) with every correction off. It tests the claim that
     constants alone cannot meet T2–T7; its stage 2 best is scored beside the corrected one.
4. **G3: held-out acceptance** of the frozen setting, at 1 Myr and at 0.5 Myr on three seeds, plus one held-out
   world for 50 Myr with full stepping (climate and surface processes on) to check that plates-only stepping did
   not mislead.

Rough cost, from #942's measured 5.7 min per 200 Myr world plates-only (0.54 GB peak) on the great-circle split:
G0's 20-world brake-off breakdown 2 h; G1 about 12 runs × 5 worlds, 6 h; G2 stage 1, 24 × 5, 11 h, and the same
again for the control arm; stage 2, 3 × 20 for each arm, 12 h; G3 about 3 h plus the full-stepping run (about
385 s/Myr, 5 h). The rift-line method costs more per step; G0 measures how much, and these figures scale with it.
All in series under the compute budget.

### Runtime and saves

B.1 integrates two columns per divergent node, B.2 changes a 3×3 sum, B.3 reads one existing field: each is a
vectorized pass over boundary nodes, expected under 1% of step time; G1 measures it. New knobs that are world
tuners become `World` fields with defaults equal to today's behaviour where they are switches, added to the tuner
list in `world.py` so old saves load with them; the constants stay module constants. Saves made before the change
load and then follow the new forces, which changes their future, as any force change does. No pickled field is
removed.

## Detailed design (provisional touchpoints)

- `backend/app/torque.py`: `ridge_push_force_vectors_n` (B.1), `slab_drag_coefficient_matrix` and
  `slab_drag_force_vectors_n` (B.2), `slab_pull_force_n` and `column_buoyancy_pa` callers (B.3, `f_pull`),
  `subducting_boundary_mask` (B.4); the diagnostics in `solve_plates_quasistatic`.
- `backend/app/lithosphere.py`: the force column from age (B.3) beside `REFERENCE_HM_OCEANIC_M`, if the cooling
  law moves here from `rift_strength` so both callers share it.
- `backend/app/mantle.py`: `MANTLE_FLOW_REFERENCE_RATE`, `generate_convection_centers`.
- `backend/app/world.py`: new switches and tuners beside `contact_drag_multiplier` and `craton_drag_multiplier`.
- RFC 0018's rift-line method (#950): `rift_strength.RiftLineParams` and `World.rift_line_params` for the rift
  knobs, `rift_strength.oceanic_strength_hm` for the cooling law B.3 shares, and `rift_damage.force_inputs`, which
  reads the same per-node forces B.1–B.4 change.
- Tests: `backend/unit_tests/test_torque.py` (each corrected form against a manufactured case: ridge push equals
  the GPE difference of two known columns; slab drag has no trench-parallel component; pull grows with age),
  `test_quasistatic.py` (the solve still converges and stays symmetric), `test_world.py` and
  `test_persistence.py` (old saves load with the new fields).
- Scripts: `bin/debug/plate_torque_breakdown.py`, `summarize_torque_breakdown.py`, `summarize_speed_census.py`
  (#942); a new `bin/debug/calibrate_plate_forces.py` (new) for G2 settings and scoring.
- Docs: `docs/simulation-model/tectonics.md` and `docs/simulation-model/foundations.md` (the mantle flow and
  torque engine sections), `docs/api-reference.md` if new tuners reach the API.

## Assumptions and uncertainty

- **The force column is a stand-in.** B.3 gives slab pull and ridge push an age-dependent thickness that the
  plate's own column does not have, so the forces describe a deeper, older ocean than the one drawn. G1 reports
  the gap; closing it is alternative D.
- **The density contrast is fixed.** The model's mantle lithosphere is 50 kg/m³ denser than asthenosphere at
  every age and depth; Earth's thermal contrast is larger near the surface and depends on temperature. `f_pull`
  and `L_sink` absorb the difference, which is why their values are calibration assumptions, not measurements.
- **The flow stays prescribed.** Slab suction, the other half of Earth's slab drive, has no counterpart except the
  random flow. A good calibration under B still has a flow pattern unrelated to the trenches.
- **Earth's mean speed is a working figure.** T1's 4–5 cm/yr is from the issue and is checked in G0.
- **Plates-only stepping.** Every gate except G3's one run leaves out climate and surface processes, which change
  relief and crust and so the brake and buoyancy. G3 checks this on one world.
- **Shared knobs with RFC 0018.** RFC 0018's rifting baseline is being re-run on the #950 code with the brake off,
  against today's forces. Calibrating its knobs here assumes RFC 0018's owner-approved targets stay as they are;
  the owner may waive RFC 0018 failures judged fixable by tuning here, and those waivers are recorded against T9.
- **The case against alternative A is argued, not measured.** It rests on one breakdown of the force forms. G2's
  control arm measures it.
- **Two-world evidence.** The brake-off breakdown and the slab work budget are on two worlds; G0 extends the
  first to 20.

## Acceptance criteria

- [ ] **G0** recorded in this RFC: the rifting baseline and the commit it ran on; T1's Earth figure checked; T2's
      definition fixed; the held-out set frozen and committed; the 20-world brake-off breakdown reported.
- [ ] **G1** recorded: each correction's predeclared probe and outcome, and the corrections adopted.
- [ ] **G2** recorded: stage 1 and 2 scores for the corrected search and the control arm, the frozen setting and
      its trade-off ridge.
- [ ] On the 20 held-out worlds at 1 Myr steps for 200 Myr, plates-only, with the frozen setting:
  - [ ] T1: mean crust speed 3–6 cm/yr;
  - [ ] T2: slab-attached plates at least 2.5 times faster than the others;
  - [ ] T3: slab pull at least 50% of the driving torque;
  - [ ] T4: mean basal traction under slabless plates at most 1 MPa;
  - [ ] T5: ridge push on old flanks 1–4×10¹² N/m;
  - [ ] T7: fastest-to-slowest world speed ratio at most 3;
  - [ ] T9: RFC 0018's realism targets, scored by its rules, at least as many met as the rifting baseline (the
        D1 re-run's setting on today's forces), apart from any the owner waives.
- [ ] T6 and T8, RFC 0008's measures, against `origin/main` on the same held-out worlds: no jammed pair has its
      collision time dropped while still in contact; no plate's speed falls below 20% of its pre-contact speed
      from braking alone; the craton-on-craton overlap volume at fusion is not more than 25% above `origin/main`'s;
      the time-averaged largest plate stays under 40% of the surface; mean plate count within 25%.
- [ ] Step check: on three held-out worlds, 0.5 Myr steps give T1 within 15% and T2 within 20% of 1 Myr steps.
- [ ] Full stepping: on one held-out world for 50 Myr, mean crust speed within 25% of its plates-only run.
- [ ] The breakdown's terms sum to zero at every sample (blocked-edge residual reported); ledgers close within
      1e-9; the solver converges at every step.
- [ ] Step time at density 4 changes by less than 2%.
- [ ] Unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green; saves from before the change
      load.

## Falsification criteria

- **Slab forces cannot drive these plates.** If G1's combined corrections, brake off, leave slab pull under a
  third of the driving torque or the T2 ratio under 1.5, the premise that the missing drive is in the force forms
  is wrong: the cause is in the column (alternative D) or the flow (alternative C). Revise toward D or C.
- **Constants alone are enough.** If G2's control arm meets at least as many of T1–T7 as the corrected search,
  with T2 and T3 among them, the corrections are not needed: adopt the control arm's setting (alternative A) and
  drop B.1–B.4.
- **The corrections only trade one knob for another.** If G2's best settings meet T1 only with
  `MANTLE_FLOW_REFERENCE_RATE` above 6 cm/yr or `contact_drag_multiplier` under 0.2, the speed still comes from the
  flow or from removing resistance, not from slabs. Revise; do not accept on T1 alone.
- **Rifting and motion cannot both be met.** If no G2 setting meets T1–T3 and RFC 0018's realism targets together,
  the two are not compatible under one force model; report which targets conflict and ask the owner which to keep.
- **The flow draw still decides a world's speed.** If T7 fails with every adopted correction, the prescribed flow
  is the dominant control whatever the slabs do, and alternative C is needed.
- **Timestep dependence.** If the step check fails, the corrections interact with the step (B.4's test reads the
  previous step), and B.4 is revised before acceptance.

## References

- <a id="ref-bird2008"></a>Bird, P., Liu, Z., & Rucker, W. K. (2008). Stresses that drive the plates from below:
  Definitions, computational path, model optimization, and error analysis. *Journal of Geophysical Research:
  Solid Earth*, 113. [doi:10.1029/2007JB005460](https://doi.org/10.1029/2007JB005460)
- <a id="ref-conrad2002"></a>Conrad, C. P., & Lithgow-Bertelloni, C. (2002). How mantle slabs drive plate
  tectonics. *Science*, 298. [doi:10.1126/science.1074161](https://doi.org/10.1126/science.1074161)
- <a id="ref-conrad2004"></a>Conrad, C. P., & Lithgow-Bertelloni, C. (2004). The temporal evolution of plate
  driving forces: Importance of "slab suction" versus "slab pull" during the Cenozoic. *Journal of Geophysical
  Research: Solid Earth*, 109. [doi:10.1029/2004JB002991](https://doi.org/10.1029/2004JB002991)
- <a id="ref-copley2010"></a>Copley, A., Avouac, J.-P., & Royer, J.-Y. (2010). India-Asia collision and the
  Cenozoic slowdown of the Indian plate: Implications for the forces driving plate motions. *Journal of
  Geophysical Research: Solid Earth*, 115. [doi:10.1029/2009JB006634](https://doi.org/10.1029/2009JB006634)
- <a id="ref-forsyth1975"></a>Forsyth, D., & Uyeda, S. (1975). On the relative importance of the driving forces of
  plate motion. *Geophysical Journal International*, 43(1), 163–200.
  [doi:10.1111/j.1365-246X.1975.tb00631.x](https://doi.org/10.1111/j.1365-246X.1975.tb00631.x)
- <a id="ref-lithgowbertelloni1998"></a>Lithgow-Bertelloni, C., & Richards, M. A. (1998). The dynamics of Cenozoic
  and Mesozoic plate motions. *Reviews of Geophysics*, 36. [doi:10.1029/97RG02282](https://doi.org/10.1029/97RG02282)
- <a id="ref-richards2018"></a>Richards, M. A., & Lenardic, A. (2018). The Cathles parameter (Ct): A geodynamic
  definition of the asthenosphere and implications for the nature of plate tectonics. *Geochemistry, Geophysics,
  Geosystems*, 19. [doi:10.1029/2018GC007664](https://doi.org/10.1029/2018GC007664)
- <a id="ref-richardson1984"></a>Richardson, R. M., & Cox, B. L. (1984). Evolution of oceanic lithosphere: A
  driving force study of the Nazca Plate. *Journal of Geophysical Research*, 89(B12).
  [doi:10.1029/JB089iB12p10043](https://doi.org/10.1029/JB089iB12p10043)
- <a id="ref-schellart2004"></a>Schellart, W. P. (2004). Quantifying the net slab pull force as a driving
  mechanism for plate tectonics. *Geophysical Research Letters*, 31.
  [doi:10.1029/2004GL019528](https://doi.org/10.1029/2004GL019528)
- <a id="ref-vansummeren2012"></a>van Summeren, J., Conrad, C. P., & Lithgow-Bertelloni, C. (2012). The importance
  of slab pull and a global asthenosphere to plate motions. *Geochemistry, Geophysics, Geosystems*, 13.
  [doi:10.1029/2011GC003873](https://doi.org/10.1029/2011GC003873)

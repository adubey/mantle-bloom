# RFC 0021: Plates are driven by their slabs: correct and calibrate the plate forces and the rifting knobs together

- **Status:** Accepted (owner merged PR #949 and confirmed acceptance on 2026-11-17)
- **Tracking issue:** [#957](https://github.com/adubey/mantle-bloom/issues/957)
- **Date:** 2026-11-17
- **Related:** #942 (diagnosis, [`analysis/issue942/`](../analysis/issue942/README.md)), #888 (RFC 0018 D1,
  where the slow plates were found), #797 (RFC 0018 tracking), #396 (RFC 0008 tracking), #399, #423, #760 (the
  contact brake and its tuners), #468 (boundary-load normalization), #783 (plate speed in the Stats panel);
  #846 and PR #950 (RFC 0018's rift-line method in the model, behind `World.split_method`), PR #955 (its
  baseline scores); RFCs
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
solver's own residual torque a negligible share of its driving torque at every sample (today at most 5×10⁻¹⁵).

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
the trench. Let `p_i` be the node's unit position and `n̂_i` the trench normal: `direction_to_neighbor` projected
into the tangent plane and normalized, `n̂_i = (n_i − (n_i·p_i) p_i) / |n_i − (n_i·p_i) p_i|`. The node's velocity is
`v_i = ω × p_i` (as in `slab_drag_force_vectors_n`), so its trench-normal speed is
`v_i · n̂_i = (ω × p_i) · n̂_i = ω · (p_i × n̂_i)`. A drag `F_i = −c_i (v_i · n̂_i) n̂_i` then has torque
`p_i × F_i = −c_i (p_i × n̂_i)(p_i × n̂_i)ᵀ ω`, so the coefficient matrix, in the same rotation space as today's
`c_i (I − p_i p_iᵀ)`, is

    K_slab = Σ_i c_i (p_i × n̂_i)(p_i × n̂_i)ᵀ,   over nodes with positive pull

`p_i × n̂_i` lies along the trench as a surface direction, but as a rotation axis it is the one that moves `p_i`
along `n̂_i`; the velocity-space form of the same drag would be `n̂_i n̂_iᵀ`. The drag acts only on the pulled
length and only on trench-normal motion. Trench-parallel motion and unpulled
buoyant length are left to basal drag, which already covers every node. `K_slab` stays symmetric positive
semidefinite, so the implicit solve is unchanged. In #942's brake-off budget the trench-parallel part alone was
46% of the pull's work; without it the pooled slab work would have netted about +24% in place of −22%.

**B.3 Age-based force column.** For slab pull and ridge push only, an oceanic cell's mantle lithosphere is taken
as `H_f = max(Hm, min(h₁ √(age / 1 Myr), H_max) − Hc)`, with `age` from `node_created_years`, `h₁` = 13 km,
`H_max` = 95 km and `Hc` the cell's crust (the oceanic reference is `REFERENCE_HC_OCEANIC_M` = 11 km). The floor at
the cell's own `Hm` matters near the ridge: without it `H_f` would be negative under about 0.7 Myr with 11 km of
crust, and B.1's column integral would be undefined there. With it a young ridge keeps the model's
`YOUNG_RIDGE_HM_M` = 8 km, and a cell without a known age keeps its `Hm`. This is the form #950's
`rift_strength.oceanic_strength_hm` already uses, which also thins the column by the rift strain applied to it. Continental cells keep their Hm. Pull per node is then `max(B(H_f), 0) × L_sink × ds` as now
(`column_buoyancy_pa`). Scale check, from the model's constants with `Hc` = 11 km: 54, 81 and 84 km at 25, 50 and
100 Myr give 6.6,
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
| RFC 0018's searched rift knobs, in `World.rift_line_params` (#950): `eps_dot0`, `s_warm`, `k_cap`, `n_plume`, `lam_grow`, `l_tip_m` | The rifting baseline: 4.6×10⁻⁸ /yr, 0.159, 0.14, 6.2×10¹¹ N/m, 0.544, 143 km | D1 stage 1 setting 0 with `k_cap` 0.14 (#888, #955) | G2, over D1's ranges |
| `mark_power` | 2.8 | #955's deeper-marks arm | G2 over 1–2.8 |
| `route_mark_power` | 1 | #955's route-marks arm used 2.8 | G2 over 1–2.8 |
| `craton_root_m` | 0 | #955: 80 km lost a target and a quarter of the splits | Fixed at 0 |
| `oceanic_age_strength` | on | #888, #955 | Fixed; shares B.3's law |
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
| T1 | Mean crust speed, no-net-rotation frame | 3–5 cm/yr | PB2002's 3.76 cm/yr in the same frame [[Bird, 2003]](#ref-bird2003) (G0; was 3–6 cm/yr against a working 4–5) |
| T2 | Speed ratio, slab-attached to other plates | at least 2.5 | Earth about 4 [[Conrad & Lithgow-Bertelloni, 2004]](#ref-conrad2004); 3–4 [[van Summeren et al., 2012]](#ref-vansummeren2012) |
| T3 | Slab pull share of the driving torque | at least 50% | [[Conrad & Lithgow-Bertelloni, 2002]](#ref-conrad2002); [[Conrad & Lithgow-Bertelloni, 2004]](#ref-conrad2004) |
| T4 | Mean basal traction under plates without slabs | at most 1 MPa | [[Bird et al., 2008]](#ref-bird2008) |
| T5 | Ridge push at ridges older than 50 Myr on the flank | 1–4×10¹² N/m | [[Copley et al., 2010]](#ref-copley2010) |
| T6 | Plate speed after continental collision | at least 20% of pre-contact | RFC 0008's criterion, from Copley's ~27% |
| T7 | Spread of world speeds, fastest to slowest | at most 3 (today 7.3) | Design choice: a world's speed should not be set by the random flow draw |
| T8 | RFC 0008's collision measures | not worse than `origin/main` (see acceptance) | RFC 0008 |
| T9 | RFC 0018's realism targets | as many as RFC 0018 defines, scored by its own rules | RFC 0018 |

T2 needs a definition that does not depend on the plate's current speed, because #942 found the current one
circular. G0 found the geometric test this paragraph first proposed (B.4b at the generation snapshot) useless: at
generation every plate has zero rotation, and the boundaries cut through continuously dated ocean floor, so 98% of
plate area counted as slab-attached at a 20% share (see G0 below). The owner chose instead (2026-11-18): a plate is
**slab-attached when it has a slab of its own**, that is, when the oceanic floor it lost under its neighbours in the
last 10 Myr (`collision_polarity`'s consumption evidence) is at least **10%** of its area. The owner fixed the share on
G0's 20-world baseline run before any G1 probe (2026-11-18): at 10% the slab-attached plates hold 32% of plate area,
close to the share of Earth's surface on the plates Forsyth and Uyeda found slab-attached, and the baseline's ratio is
1.5–1.7 at every share from 5% to 50%, so the choice does not decide the verdict. The test is built by past convergence, not the plate's current speed.

### Diagnosis gates

Development worlds are RFC 0018's 20 (#888); every gate runs plates-only (`simulate_climate_biomes = False`) at
1 Myr steps for 200 Myr unless stated, and reports at 0.5 Myr on three of them. From G0 on, every run splits plates
with `World.split_method = "rift_lines"` (#950) at the rifting baseline below, unless stated.

**Rifting baseline.** RFC 0018's D1 strength-map arms were rerun on the production method (#955, on #950 at
`6bda6906`; `analysis/issue888/rift-lines-baseline/`), brake off, on the five stage 1 worlds. D1's prototype tables
are not a baseline: #950's reviews changed plume heat at the split test, flank transfer under the crust cap and
capped donors.

| Brake off, `k_cap` 0.14 | Base | Craton root | Deeper marks | Route marks |
|---|---:|---:|---:|---:|
| RFC 0018 targets met (of 10) | 5 | 4 | 6 | 6 |
| Splits per world per 100 Myr (need 3.85) | 4.5 | 3.5 | 7.8 | 3.1 |
| Splits mostly on weak rock (need 67%) | 53% | 49% | 79% | 71% |
| Completed-rift median (Myr, 10–50) | 63 | 70 | 64 | 52 |
| Failed share (20–80%) | 93% | 90% | 88% | 93% |
| Transition ratio (need 3) | 2.75 | 1.62 | 1.50 | 1.04 |
| Localization violations (need 0) | 2 | 1 | 9 | 0 |

The **rifting baseline** for G0 and G1 is the deeper-marks arm (`mark_power` 2.8): it ties route marks on targets
met and has the most splits, so it shows the split rate's response to the forces most clearly. Route marks enter
G2 through `route_mark_power`. The failed share, the completed-rift duration and the transition ratio fail on every
arm. They are what T9 most needs from the joint calibration, and the ones the owner may waive. Each plates-only
world took 5–6 minutes and about 0.5 GB on the rift-line method, as on the great-circle split.

1. **G0: targets and harness.**
   - Record the rifting baseline above in this RFC, with its run labels and commit.
   - Check T1's Earth figure against a published plate model (area-weighted mean of `|ω × r|` from published
     Euler poles and plate areas), and record the result. If it is outside 4–5 cm/yr, T1's band moves with it,
     with the owner's agreement.
   - Fix T2's slab-attached definition and threshold before any probe (revised by the owner after the generation
     census failed: see T2 above).
   - Freeze a held-out set for this RFC by the RFC 0018 scheme (SHA-256 of
     `mantle-bloom/rfc0021/v1/heldout/{index:02d}`, first four bytes big-endian masked with `0x7fffffff`,
     taking indices from 0 upward and skipping any value that repeats an earlier one after the mask or is already
     in `analysis/issue797/g2-seed-preregistration.json` or the development set, until 20 seeds are taken), with
     an assertion in the freezing script that the set and the development set are disjoint, committed before any
     G2 run.
   - Extend #942's brake-off breakdown from 2 to the 20 development worlds, on the rifting baseline, so G1's
     baseline is not two worlds and not the great-circle split.
   - Build the scorer: `plate_torque_breakdown.py`'s measures plus RFC 0008's collision measures plus RFC 0018's
     realism score from #955's `bin/debug/rift_line_recorder.py` and `score_rift_lines.py`, one JSON per setting.
     The breakdown does not measure two targets yet, and the scorer adds them:
     - **T4:** per node, `|basal_drag_force_vectors_n|` over the node's accounting area (Pa), area-weighted over
       plates with no node in the subducting mask, pooled over samples.
     - **T5:** per divergent node, the ridge-push force over its shared front length (N/m), length-weighted over
       ridges whose plate has oceanic floor at least 50 Myr old, the age at which B.3's column reaches 81 km of
       its 84.
     - **Closure:** the solver's `residual_torque_norm_nm / driving_torque_norm_nm` per sample, which
       `summarize_torque_breakdown.py` now reports. The breakdown's own terms sum to zero by construction,
       because the blocked-edge term is taken as the remainder, so they check nothing.
2. **G1: one correction at a time.** Each of B.1, B.2, B.3 and B.4 (a and b) alone, then all adopted ones
   together, brake on and brake off, on the five stage 1 worlds #955 scored, as causal probes with predictions posted before
   each run. A correction is adopted if it moves its own measure the predicted way (B.1: ridge push into T5's
   band; B.2: slab work positive; B.3: pull per metre up; B.4: gated share of closing down-going length up) and
   breaks no invariant (the solver's residual stays under 10⁻⁹ of its driving torque; the solver converges; ledgers
   close). G1 reports each one's
   effect on speed, T2 and the split rate, and what alternative D would change in bathymetry and land fraction.
   - **Best-case check.** G1 then runs the combined adopted corrections, brake off, at the four most
     slab-favourable corners of G2's ranges: `ASTHENOSPHERE_VISCOSITY_PA_S` 2×10¹⁹ and `f_pull` 1.0, with
     `SLAB_MANTLE_VISCOSITY_PA_S` 10²⁰ or 10²¹ and `SUBDUCTION_LSINK_M` 150 or 660 km. If none of them reaches a
     T2 ratio of 2.0, the premise fails (see falsification) before G2 is spent. Untuned constants can't explain
     a miss at the ends of the ranges G2 would search. The check costs about 4 settings × 5 worlds.
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

**G0 result (#958; evidence in [`analysis/issue958/`](../analysis/issue958/README.md)).**

- *Rifting baseline:* #955's `deep-marks` arm of `analysis/issue888/rift-lines-baseline/arms.json`, run on PR #950 at
  `6bda6906` with `run_rift_lines_arms.py` and `score_rift_lines.py` at `bd7b2c3b`; run labels
  `issue888/rift-lines-baseline/deep-marks/<seed>` for the five stage 1 worlds. On `main` it is the
  `rifting-baseline` setting (brake on) and `rifting-baseline-brake-off` setting of `analysis/issue958/settings.json`.
- *T1:* PB2002's 52 plates [[Bird, 2003]](#ref-bird2003), with their net rotation removed, move at an area-weighted
  mean of 3.76 cm/yr (RMS 4.39; 6.06 Pacific-fixed), below the working 4–5 cm/yr. The owner moved T1 to 3–5 cm/yr in
  the no-net-rotation frame (2026-11-18). `bin/debug/earth_plate_speed.py`, `earth-plate-speed.json`.
- *T2:* on the 20 generation snapshots (282 plates) B.4b's test marks a median 63% of a plate's boundary band
  slab-attached, 98% of plate area clears a 20% share, and continental-majority plates (median 69%) are not
  separated. Generation leaves every plate at zero rotation, so "not diverging" passes everywhere, and ocean floor
  ages across an ocean–ocean contact differ by a median 3.3 Myr, so the older side is a coin toss. After 10 Myr of
  stepping the test still marks about half of every band, with a rank correlation of −0.25 to +0.28 with speed.
  The owner replaced the definition (see T2 above). `t2-census.json`.
- *Held-out set:* frozen in `analysis/issue958/heldout-seeds.json` by `bin/debug/freeze_rfc0021_heldout.py`
  (20 seeds, none skipped, disjoint from the development set and RFC 0018's preregistration).
- *Harness:* `bin/debug/calibrate_plate_forces.py` runs one setting per world, recording the torque breakdown with
  T2's, T4's and T5's measures, no-net-rotation speeds and the solver closure, RFC 0008's contact probe
  (`probe_craton_contacts.py`, whose tap now passes #760's brake gains) and RFC 0018's recorder, and scores a setting
  on T1–T5, T7, T8 and T9.
- *20-world brake-off breakdown* (`rifting-baseline-brake-off`, 200 Myr, 1 Myr steps, sampled every 10 Myr;
  `g0-baseline-brake-off.json`):

  | Target | Baseline | Band | |
  |---|---:|---|---|
  | T1 mean speed, no net rotation | 1.88 cm/yr | 3–5 | fails |
  | T2 speed ratio, slab-attached (32% of area) to the rest | 1.54 | ≥ 2.5 | fails |
  | T3 slab pull share of driving torque | 55% | ≥ 50% | passes |
  | T4 basal traction under plates without slabs | 0.22 MPa | ≤ 1 | passes |
  | T5 ridge push, oceanic floor ≥ 50 Myr | 4.5×10¹⁰ N/m | 1–4×10¹² | fails, 20–90× low |
  | T7 fastest to slowest world | 2.0 | ≤ 3 | passes |
  | T9 RFC 0018 targets met | 5 | | split rate unscored (needs the great-circle rate on these worlds) |

  Resisting torque is sunk-slab drag 51% and basal drag 45%; ridge push is 0.6% of the driving torque. The solver
  closure is at most 6.5×10⁻¹⁵ and every sample converged. The harness reproduces #942's two brake-off worlds to
  four decimals. T6 is recorded per step but not yet scored. The baseline fails T1, T2 and T5 and passes T3, T4 and
  T7: slab pull already leads the drive, but sunk-slab and basal drag cancel it, as #942 found, so G1's corrections
  are the next step.
- *0.5 Myr step check* (worlds 106533960, 1097624328, 1397817354; `g0-step0.5-brake-off.json`): T1 +5.1% and T2
  −4.7% against the same worlds at 1 Myr, inside the 15% and 20% of the step-check criterion. T3, T5 and T7 move
  under 4%; T4 halves (0.91 to 0.47 MPa) on this small sample. The rift-line splitter is not step-converged: 73
  splits against 48 (+52%; 38 against 23 on 106533960), 19.2 plates against 17.5. RFC 0018's split rate therefore
  depends on step length, and G2 and G3 must score it at 1 Myr steps throughout.

Rough cost, from 5–6 min per 200 Myr world plates-only (about 0.5 GB) on either split method (#942, #955):
G0's 20-world brake-off breakdown 2 h; G1 about 12 runs × 5 worlds, 6 h, and its best-case check 4 × 5, 2 h; G2 stage 1, 24 × 5, 11 h, and the same
again for the control arm; stage 2, 3 × 20 for each arm, 12 h; G3 about 3 h plus the full-stepping run (about
385 s/Myr, 5 h). All in series under the compute budget.

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
  the GPE difference of two known columns; slab drag resists a purely trench-normal rotation and gives zero torque
  for a purely trench-parallel one, including with a `direction_to_neighbor` that is not exactly tangent; pull
  grows with age; the force column never goes below the cell's `Hm`),
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
- **Earth's mean speed is frame-dependent.** T1's figure is PB2002's in the no-net-rotation frame (G0); in the Pacific-fixed
  frame the same model gives 6.06 cm/yr. The model's prescribed flow defines a mantle frame, so T1 removes the model's own net
  rotation before comparing.
- **Plates-only stepping.** Every gate except G3's one run leaves out climate and surface processes, which change
  relief and crust and so the brake and buoyancy. G3 checks this on one world.
- **Shared knobs with RFC 0018.** The rifting baseline (#955) was scored with the brake off against today's forces,
  on five worlds. Calibrating its knobs here assumes RFC 0018's owner-approved targets stay as they are; the owner
  may waive RFC 0018 failures judged fixable by tuning here, and those waivers are recorded against T9. The base
  arm's transition ratio rose from 1.02 in D1 to 2.75 on the production method, and which of #950's changes did
  that has not been checked.
- **The case against alternative A is argued, not measured.** It rests on one breakdown of the force forms. G2's
  control arm measures it.
- **Two-world evidence.** The brake-off breakdown and the slab work budget are on two worlds; G0 extends the
  first to 20.

## Acceptance criteria

- [x] **G0** recorded in this RFC: the rifting baseline's run labels and commit; T1's Earth figure checked; T2's
      definition fixed; the held-out set frozen and committed; the 20-world brake-off breakdown reported.
- [ ] **G1** recorded: each correction's predeclared probe and outcome, and the corrections adopted.
- [ ] **G2** recorded: stage 1 and 2 scores for the corrected search and the control arm, the frozen setting and
      its trade-off ridge.
- [ ] On the 20 held-out worlds at 1 Myr steps for 200 Myr, plates-only, with the frozen setting:
  - [ ] T1: mean crust speed 3–5 cm/yr in the no-net-rotation frame;
  - [ ] T2: slab-attached plates at least 2.5 times faster than the others;
  - [ ] T3: slab pull at least 50% of the driving torque;
  - [ ] T4: mean basal traction under slabless plates at most 1 MPa;
  - [ ] T5: ridge push on old flanks 1–4×10¹² N/m;
  - [ ] T7: fastest-to-slowest world speed ratio at most 3;
  - [ ] T9: RFC 0018's realism targets, scored by its rules (`bin/debug/score_rift_lines.py`, #955), at least as
        many met as the rifting baseline's 6 of 10, apart from any the owner waives.
- [ ] T6 and T8, RFC 0008's measures, against `origin/main` on the same held-out worlds: no jammed pair has its
      collision time dropped while still in contact; no plate's speed falls below 20% of its pre-contact speed
      from braking alone; the craton-on-craton overlap volume at fusion is not more than 25% above `origin/main`'s;
      the time-averaged largest plate stays under 40% of the surface; mean plate count within 25%.
- [ ] Step check: on three held-out worlds, 0.5 Myr steps give T1 within 15% and T2 within 20% of 1 Myr steps.
- [ ] Full stepping: on one held-out world for 50 Myr, mean crust speed within 25% of its plates-only run.
- [ ] At every sample the solver's residual torque is under 10⁻⁹ of its driving torque, and the blocked-edge
      share (the breakdown's remainder) is reported; the crust and mantle-lithosphere ledgers of rift thinning
      (`rift_runtime.ledger_gap`, relative, worst over steps) close within 1e-9; the solver converges at every
      step.
- [ ] Step time at density 4 changes by less than 2%.
- [ ] Unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green; saves from before the change
      load.

## Falsification criteria

- **Slab forces cannot drive these plates.** If no setting of G1's best-case check reaches a T2 ratio of 2.0,
  the premise that the force forms cancel or omit the drive is wrong: the cause is in the column (alternative D) or
  the flow (alternative C). Revise toward D or C.
  - This replaces a bar of slab pull under a third of the driving torque or a T2 ratio under 1.5 after G1's
    combined corrections. G0's uncorrected baseline already clears both (55%, 1.54), so that bar could not fail.
  - The bar sits at the end of G1, at the slab-favourable ends of G2's ranges, so that a miss can't be blamed on
    untuned constants. 2.0 is about halfway from the baseline to T2's 2.5.
  - The owner agreed this (2026-11-18).
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

- <a id="ref-bird2003"></a>Bird, P. (2003). An updated digital model of plate boundaries. *Geochemistry, Geophysics,
  Geosystems*, 4(3). [doi:10.1029/2001GC000252](https://doi.org/10.1029/2001GC000252)
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

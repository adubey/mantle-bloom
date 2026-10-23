# RFC 0012 gate G1: force-proxy inventory and rule audit

Gate G1 of [RFC 0012](../../rfc/0012-physics-first-simulation.md) (issue #465, tracking #464).
Diagnosis only: no simulation behaviour changed. Commit measured: `dd3c620` (origin/main at the
start of the gate).

| Deliverable | Script | Output |
|---|---|---|
| Load normalization | [`bin/debug/probe_boundary_load_normalization.py`](../../bin/debug/probe_boundary_load_normalization.py) | [`load_normalization.json`](load_normalization.json) |
| Which rules bind, order dependence, band rows | [`bin/debug/audit_force_rules.py`](../../bin/debug/audit_force_rules.py) | `rules_<seed>.json` |
| Step cost by phase | `audit_force_rules.py --profile` | `profile_d1.json`, `profile_d4.json` |

## Summary

- **Load normalization fails.** Slab pull, slab drag and ridge push are summed per band
  node, so their load per boundary metre is the band's row count times the stated value:
  2.5–2.9× in synthetic cases, doubling with repeated rows, about 2.8× in generated worlds.
  Contact drag divides by the nominal band width and still varies by 20–100%. Under the
  outcomes fixed in #465 this becomes a normalization issue that G2 depends on.
- **Motion is already quasi-static and memoryless.** Relaxation times are ≤ 2×10⁻⁶ s; the
  inertia term and the ω assigned at merges and splits have no effect on the next step.
- **Shift order moves ω by a median 9–18%**, at worst several times over. A simultaneous
  solve from one snapshot is justified on order dependence alone.
- **The rate clamp binds:** the cap on 0.5–2.7% of plate-steps (≥ 1% on three seeds), the floor
  on 0.6–2.6%. G2 must report the force balance the cap hides.
- **The split "stress proxy" decides** which plates may split in 75–88% of attempts past the
  residual gate; it contains no force. The 0.55 roll fails 39–47% of attempts that clear every gate.
- **Merges happen at slow but moving contacts** (median |closing| ≈ 0.5–0.9 cm/yr, some
  opening), so G3's "locked for τ" trigger needs ε near 1 cm/yr to keep today's fusions.
- **Yield is a 1.49 cm/yr closing-rate threshold** that 74–87% of convergent nodes pass; the
  rate law never reaches saturation. Caps bind on a few percent of continental nodes.
- **Cost:** the solves are free; gathering is the motion phase's cost.

## 1. Verified inventory

### Corrections to the RFC's Problem table

| Row | What the code does | Correction |
|---|---|---|
| Rate clamp | `mantle.clamp_rate` (`mantle.py:134`) clips \|ω\| to [0.5, 15] cm/yr. Callers: `torque.integrate_omega` (every plate, every step), `merge_split.maybe_split_plate` (both daughters' fitted poles), `quad_merge.merge` (the merged plate), and `torque.merge_omega`. | `torque.merge_omega` is called only by `test_torque.py`. Production merges blend angular momentum in `quad_merge.merge` (`quad_merge.py:90`), less the momentum of crust lost in the transfer. |
| Merge timer | `update_collision_progress` accumulates convergent years per pair against a per-pair uniform 0–100 Myr threshold (`_collision_threshold_years`, seeded by the pair), then rolls `_merge_probability` each step: 1 for two small plates, falling to 0.02 as their combined node share reaches 0.25, raised by up to 0.5 × the remainder with \|ω_a − ω_b\| / 15 cm/yr, all × 0.85. A jammed pair keeps accumulating (RFC 0008, #401). Backstop: a continental pair overlapping ≥ 30% for 30 Myr fuses (`pop_ready_forced_merge`). At most one fusion per step. | As stated; the speed boost uses relative ω, not the contact's closing rate. |
| Split | `maybe_split_plate` gates, in order: ≥ 700 × density nodes, ≥ 2 Myr age, RMS misfit of the *prescribed mantle flow* (`mantle.flow_at`) to the plate's rigid rotation above 6 cm/yr × (1 − relax), k-means into two flow regimes, cut refinement, fitted-pole separation above 4 cm/yr × (1 − relax), both halves large enough, then a 0.55 roll (else a failed rift thins a 2.5-spacing band by 10%). `relax` is the larger of plate radius / 2.2 rad and `internal_stress` / 3, where `internal_stress` accumulates from plate radius and overlap and decays 2%/Myr. | The "stress proxy" is not a stress: it is an accumulator of plate size and overlap, with no force in it. The velocity being fit is the mantle flow, not plate or boundary motion. At most one split per step. |
| Convergent yield | `rheology.plastic_strain_rate_per_myr`: σ_n = 10¹⁷ Pa·s/m × closing rate; Mohr–Coulomb Y = 20 MPa + σ_n tan 30°; rate = 0.06/Myr × (σ_n − Y)/Y. | Y grows with σ_n, so yield is a closing-rate threshold: σ_n (1 − tan 30°) > 20 MPa, i.e. closing > 1.49 cm/yr. Above it the rate rises toward 0.06 × (1 − tan 30°)/tan 30° = 0.044/Myr, independent of plate strength, column, or load. On quad plates the band strain is carried inland by `shortening.ring_cascade` (#314), whose drive stress is the same proxy and whose relief pushback (`CASCADE_RELIEF_PUSHBACK` = 0.5 × ρ_c g relief) is already a gravitational resistance term. |
| Collision drag | `contact_drag_pa_s_per_m`: η = 5×10¹² N/m / (W v_ref), W = 3 spacings, v_ref = 4 cm/yr, times one bounded factor g = max(jam, relief / 4 km, craton strength) ∈ [0, 1], folded implicitly into the plate's own (b, K). | The relief, craton and jam factors are one `max`, not multipliers (RFC 0008). It is a velocity-proportional drag: zero at zero closing, so it cannot hold a contact closed. |
| Thickness caps | `MAX_CRUSTAL_THICKNESS_M` = 2.4 × 35 km = 84 km, `MAX_MANTLE_LITHOSPHERE_THICKNESS_M` = 2.4 × 100 km = 240 km (`lithosphere.py:89`); suture accretion shares the Hc cap; the cascade tapers uptake from 75% of either cap. | As stated. |
| Oceanic Hm relaxation | `relax_oceanic_mantle_lithosphere`: Hm → 40 km at 1 − e^(−0.1/Myr × dt), every oceanic node, no age. | As stated. |
| Step order | `step_world_progress` shifts plates sequentially in `world.plates` order (`world.py:900`); `deform()` runs in an order shuffled by (seed, elapsed years) (`world.py:917`). | As stated; measured below. |
| Merge ω | See the rate-clamp row: `quad_merge.merge`, then `clamp_rate`. | Location corrected. |
| Faults | `faults.py`: spawn rate 3/Myr decaying over 500 km from boundaries (interior floor 0.03), lognormal lengths (median 45 km, σ 0.6, 12–200 km), slip 150–30,000 m/Myr scaled by local closing rate, lifespans 2–25 Myr, 18% as long fault systems. Relief only. | As stated; no load returns to the motion solve. |

### Process inventory

| Process | Inputs | Output and units | Spatial support | State written | Feeds plate motion? | Material owner |
|---|---|---|---|---|---|---|
| Basal drag | own nodes, ω, prescribed mantle flow (`mantle.flow_at`) | affine torque (b, K), μ/d_s = 5×10¹⁹ / 10⁵ Pa·s/m over node area (N·m) | every node, area-weighted | none | yes, implicit; it sets the drive | none |
| Slab pull | band nodes (≤ 3 spacings) of down-going crust converging > 1 cm/yr; Hc, Hm, craton | per node B × 250 km × ds (N), B = column negative buoyancy (Pa) | boundary band, each band node × one spacing | none | yes, explicit | none |
| Slab drag | same nodes | K, 10²¹ / 10⁵ Pa·s/m × 250 km × ds per node | boundary band | none | yes, implicit | none |
| Ridge push | band nodes classified divergent; isostatic elevation | per node ½ g (ρ_a − ρ_w) e_r² ds (N), e_r = elevation above the abyssal reference | boundary band | none | yes, explicit | none |
| Collision / jam drag | contested continental band, jam factor, relief, craton | K and b from η g × node area (Pa·s/m × m²) | boundary band, divided by nominal band width | none | yes, implicit | none |
| Convergence (deform) | closing rate per band node (m/s) | fractional thickening per Myr, rate law above; ring cascade inland | band, then continental interior rings to 1,500 km | Hc, Hm, cover, craton, restite, thermal lag; ledgers | no: next step only, through Hc/Hm | per-plate `deform()`, shuffled order |
| Subduction | contested oceanic nodes, polarity prepass | node deletion, arc magmatism (450 m/Myr at 5 cm/yr, ×3 cap) | band | node set, Hc; ledgers | no (slab pull is computed separately at shift) | `deform()` |
| Rifting (boundary) | closing < −1 cm/yr | thinning, rift magmatism (150 m/Myr), new crust below 5 km | band | Hc, Hm | no | `deform()` |
| Transform | \|closing\| ≤ 1 cm/yr | classification only; faults take strike-slip regime | band | none | no shear resistance | none |
| Merge | timer, roll, overlap backstop | topology; ω from momentum blend | plate pair | plates, ω | ω reset (momentum blend + clamp) | `quad_merge` |
| Split | mantle-flow misfit, size, stress accumulator, roll | topology; daughters' ω from fitted mantle-flow poles | whole plate, great-circle cut | plates, ω, `internal_stress` | ω reset to fitted poles + clamp | `Plate.split` |
| Faults | distance to boundary, closing rate | relief | intraplate and boundary | fault list, elevation | no | `faults.py` |

Every load that reaches the motion solve is computed per plate from its own boundary band.
None is assembled per shared contact, so nothing makes the two sides' loads equal and
opposite: an ocean–ocean trench pulls both plates (`down_going_crust` is true for every
oceanic node), and each side's contact drag is computed from its own band.

## 2. Load normalization

Two oceanic plates split along the 0/180° meridian great circle. The west plate rotates east
at 5 cm/yr about the pole, so it closes on a trench at 0° and opens from a ridge at 180°;
only nodes within 60° of the equator count. Each load is divided by the counted boundary
length (R × 120°) and by its ideal value per metre for the uniform column (slab pull:
B × 250 km; ridge push: ½ g (ρ_a − ρ_w) e_r²; slab drag: (μ/d_s) × 250 km × v; contact drag:
`CONTACT_RESISTANCE_N_PER_M` at 4 cm/yr). Each density is measured at six offsets of the
east plate's lattice over one spacing: how the two lattices line up sets the gap between
their nearest nodes (0.05–1.8 spacings here), and the band is "within 3 spacings of a
neighbour node", so the gap changes how many rows fall in it. Values are the mean over the
six alignments; the alignment spread (max/min − 1) is in parentheses.

| Sweep | Slab pull | Ridge push | Slab drag | Contact drag | Spread across levels (pull / push / drag / contact) |
|---|---|---|---|---|---|
| density 1 | 2.45 (64%) | 2.28 (54%) | 2.40 (60%) | 0.82 (68%) | |
| density 2 | 2.89 (20%) | 2.13 (28%) | 2.80 (20%) | 0.98 (21%) | |
| density 4 | 2.82 (65%) | 2.11 (31%) | 2.73 (61%) | 0.94 (69%) | 18% / 8% / 16% / 20% |
| reach 1.5 spacings | 0.90 | 0.76 | 0.87 | 0.60 | |
| reach 3 (today) | 2.45 | 2.28 | 2.40 | 0.82 | |
| reach 6 | 5.41 | 5.31 | 5.28 | 0.91 | 502% / 600% / 503% / 52% |
| rows ×1 | 2.45 | 2.28 | 2.40 | 0.82 | |
| rows ×2 | 4.90 | 4.55 | 4.80 | 1.64 | 100% / 100% / 100% / 100% |

**Fails the ≤ 10% criterion for every load.** Slab pull, ridge push and slab drag add one
spacing of boundary (`ds`) for every band node, so their load per metre is the number of band
rows times the ideal: about (reach − gap) / spacing. Contact drag divides by the band's
nominal width (`BOUNDARY_FORCE_REACH_MULTIPLIER` spacings), which corrects the reach sweep's
first-order scaling but not the gap or repeated rows. Density alone moves the means by
8–20%, within the alignment spread, so the density dependence is the lattice gap, not the
spacing.

In generated worlds (§3, band rows), the subducting and divergent bands hold 2.6–3.0 nodes per
neighbour front node on every seed and checkpoint. So the model's slab pull and slab drag act
at about 2.8× their stated per-metre values, and ridge push at about 2.7×. Slab pull and slab
drag scale together, so a trench's own pull–drag balance does not move; both are inflated
against basal drag, which is integrated over area and normalized correctly.

## 3. Which rules bind

Seeds 326494039, 748935664, 580014012 and 111065103; `generate_world(seed)` at default
controls and density 1; 2 Myr steps to 500 Myr. Run times 11–13 min per seed, peak RSS
0.6–0.8 GB. Ranges below are over the four seeds; per-seed values are in `rules_<seed>.json`.

| Rule | Measured | Binds? |
|---|---|---|
| `clamp_rate` floor (0.5 cm/yr), at `integrate_omega` | 0.6–2.6% of 7,400–9,500 plate-steps | yes, ≥ 1% on three seeds |
| `clamp_rate` cap (15 cm/yr), at `integrate_omega` | 0.5–2.7% of plate-steps; unclipped 99th-percentile rate 8.7–26 cm/yr (median 2.9–4.3) | yes, ≥ 1% on three seeds |
| `clamp_rate` at merges (`quad_merge.merge`) | floor 0.8–3.2%, cap 0–2.4% of 121–158 calls | yes |
| `clamp_rate` at split poles | floor 0–1.2%, never the cap | rarely |
| Plate inertia in `integrate_omega` | slowest relaxation time 0.9–2.0×10⁻⁶ s (median 5–9×10⁻⁸ s); a 2 Myr step is ≥ 10¹⁹ of them; the step's ω differs from the quasi-static K⁻¹(τ + b) by at most 5×10⁻¹² | no: the solve is already quasi-static |
| Merge timer + roll | 6–14 fusions per seed; 4–11 pair-steps past the timer failed the roll | yes |
| Forced deep-overlap merge | 1 fusion in 41 (seed 326494039); 0 on the others | rarely |
| Contact at fusion | per-fusion median \|closing\| 0.2–6.9 cm/yr, typically 0.5–0.9; mean closing negative (opening) at 6 of 41; \|ω_a − ω_b\| R median 1.4–2.2 cm/yr | merges happen at slow, not locked, contacts |
| Split size gate | rejects 1,900–4,000 calls per seed | yes |
| Split mantle-flow residual gate | rejects 2,000–2,800 calls | yes |
| Split relaxation source, past the residual gate | the stress accumulator 75–88%, plate size 6–25%, certain (relax = 1) 0.3–6% | the accumulator decides |
| Split geometric gates (k-means, cut, pole separation, half sizes) | reject 1,900–2,700 calls | yes |
| Split 0.55 roll | 33–61 splits and 21–55 failed rifts per seed (39–47% failed) | yes |
| Convergent yield | 74–87% of 1.05–1.31 M convergent node-steps past yield; none within 5% of the rate law's saturation | the threshold excludes only 13–26%; the saturation branch never binds |
| Ring cascade | 84–92% of demanded shortening absorbed; 4–8% lost to decay, 4–9% unrouted | — |
| Hc / Hm caps | continental nodes at 99% of cap, mean over 10 Myr samples: Hc 0.6–4.0%, Hm 1.4–3.7%; peak Hc 3–21%, Hm 3–13% | yes |

**Order dependence.** On a copy of the world at 100, 300 and 500 Myr, one shift pass in three
shuffled plate orders moves each plate's ω, relative to the model's own order, by a median of
9–18% and at worst 45–446%. A plate's forces are gathered against neighbours that may already
have rotated this step, so the step's motion depends on list order.

**Band rows.** At the same checkpoints the subducting band holds 2.68–2.99 nodes per distinct
neighbour front node and the divergent band 2.62–2.95, on every seed: the factor §2 says slab
pull, slab drag and ridge push are scaled by.

**Motion has no memory between steps.** Because the solve is quasi-static, the ω a merge
(angular-momentum blend) or a split (fitted mantle-flow poles) assigns is overwritten at the
next shift, to within 10⁻¹². It still sets the plate's motion for the rest of that step's
topology and surface phases, and the clamp acts on it, but not on the next step's motion.

## 4. Cost

`audit_force_rules.py --profile`, seed 326494039, 2 Myr steps, five steps after one warm-up
step from generation, Python 3.14.6, Apple M4. `gather` and `classify` run inside
both `shift` and `deform`, so the rows overlap; phases are seconds per step.

| Phase | Density 1 (30,847 nodes, 20 plates) | Density 4 (123,802 nodes, 20 plates) |
|---|---|---|
| Median total step | 1.90 s | 4.97 s |
| `deform` (all plates) | 0.234 | 1.029 |
| `shift` (all plates) | 0.126 | 0.323 |
| `gather_boundary_force_inputs` | 0.152 | 0.379 |
| `classify_boundary_nodes` | 0.055 | 0.163 |
| force terms and inertia | 0.022 | 0.064 |
| topology (`apply_topology_changes`) | 0.023 | 0.082 |
| 3×3 solves (`integrate_omega`) | < 0.001 | < 0.001 |
| rigid rotation | 0.001 | 0.002 |

Plate motion and deformation are 27% of a density-4 step; the rest is climate, erosion and
the surface passes. Inside `shift`, neighbour gathering dominates and the solve is free. A
90-unknown dense solve is microseconds, so G2's 10% budget (≈ 0.5 s at density 4) goes to
segment assembly, which must reuse the shift's gather rather than add one. The profile is of
an early world (20 plates); at 500 Myr seeds carry 31–48 plates, so G2 should profile a
replayed late world too.

## Outcome

Applying the outcomes fixed in #465 before the runs:

- **Normalization spread > 10%** (every load, every sweep): a per-boundary-metre normalization
  of slab pull, slab drag, ridge push and contact drag is its own issue, a prerequisite of G2.
- **Retirable** (binds on < 1% of plate-steps or events on all four seeds, with a physical
  replacement already supplied): the plate inertia term in `integrate_omega`, and with it the
  merge angular-momentum blend and split fitted poles as motion state; the quasi-static
  solve already supplies the motion. No other rule qualifies: the forced-overlap merge fired
  once and the split-pole clamp on up to 1.2%. G2 checks each against accepted metrics
  before anything is removed.
- **Cap ≥ 1% of plate-steps** on three seeds: G2 reports what force balance it hides, and
  measures the floor's removal.
- **G2's baseline.** The pilot replaces the sequential per-plate shift with one simultaneous
  solve. It retains `contact_drag_pa_s_per_m` as a symmetric pairwise damper and adds a
  bounded normal force for blocked edges. Basal drag (already area-normalized), the
  rheology rate law and ring cascade (deformation, not load), and merge/split triggers (G3)
  remain in the comparison. Its order-dependence and clamp numbers are compared against
  this audit's.

## G2 pilot (preliminary; issue #476)

The pilot assembles plate-local basal and slab drag into one `3P` system, then adds symmetric
pair blocks for contact drag and a bounded normal-force column for each blocked edge carried
by the preceding frozen polarity frame. The force cap follows RFC 0011's candidate values:
5×10¹² N/m for deformable continental material, rising linearly to 1.5×10¹³ N/m with the
weaker side's craton strength. No per-node state is added. Equation tests cover zero load,
unbalanced torque, equal-and-opposite contact force, held compression, yield at capacity,
and the column-GPE reference integral.

**Early runtime check.** Same generated world and hardware, seed 326494039, node density 4,
1 Myr steps, one warm-up step followed by five timed steps. The baseline is merged
`origin/main` commit `2279d7e` (PR #473); the pilot is that commit plus the uncommitted G2
working tree. Median full-step runtime was 4.379 s baseline and 4.684 s pilot (**+6.9%**),
under the 10% target for this early-world sample. Timed baseline steps were 4.358, 4.379,
5.433, 4.433 and 4.349 s; pilot steps were 4.414, 4.428, 5.479, 5.456 and 4.684 s.
The command in each worktree generated the same seed and density, stepped once to warm up,
then measured five complete `step_world` calls.

**Initial force-balance sample.** `probe_physics_first_g2.py --seed 326494039 --steps 3
--step-myr 2 --node-density 1` recorded 19 plates and no blocked edges in the first 6 Myr.
The largest plate rate fell from 7.97 to 3.66 cm/yr over those steps. The relative torque
residual was below 4×10⁻¹⁵ of the driving-torque norm each step. This sample does not exercise
the bounded collision force.

**Completed in this evidence set:** four-seed 2 Myr / 500 Myr matched runs, four-seed 0.1/1 Myr
checks to 50 Myr, and geometry-defined ridge/collision ΔGPE scale checks. Still required before
recording a G2 outcome are full-horizon 0.1/1 Myr comparisons, the late-world density-4 profile,
and ledger/topology diagnosis across matched baselines. The corrected 500 Myr replay below
exercises blocked contacts and yielding, but its large geological/topological drift means the
candidate is not accepted. These remain pilot-development results, not a gate pass or a selected
general implementation plan.

## Extended G2 pilot (one seed, 500 Myr)

Seed 326494039, density 1, 2 Myr steps, generated independently from the same origin/main
commit `2279d7e` in each worktree. The baseline is the detached `origin/main` worktree; the
pilot includes the uncommitted G2 changes. Both runs advanced 250 complete world steps. Raw
JSON is under `/private/tmp/rfc0012-g2-baseline-500myr.json` and
`/private/tmp/rfc0012-g2-corrected-500myr.json`; these replay outputs are intentionally not
checked in.

The first 500 Myr pilot had 2,782 blocked edges but zero contact force. Inspection found that
the active-set solver rejected compliance smaller than machine epsilon, even though SI
contact compliance (m/s per N/m) is naturally much smaller. The sign check now accepts any
negative finite compliance, and an SI-scale equation test covers this case. The corrected
50 Myr smoke replay recorded 25 blocked edges and 18 yielded contacts.

| Metric at 500 Myr / across 250 steps | Baseline | G2 pilot |
|---|---:|---:|
| Median complete step | 2.588 s | 2.587 s |
| Peak plate count | 39 | 46 |
| Peak plate rate | 15.0 cm/yr (clamped) | 303.8 cm/yr |
| Steps with peak rate above 15 cm/yr | 43 | 26 |
| Steps containing yielded contacts | — | 233 |
| Yielded contacts summed over steps | — | 1,382 |
| Maximum contact force | — | 1.5×10¹³ N/m |
| Maximum relative torque residual | — | 7.6×10⁻¹⁴ of driving torque |
| Land fraction | 0.1436 | 0.2199 |
| Area-weighted Hc | 18.11 km | 19.51 km |
| Area-weighted Hm | 49.58 km | 55.10 km |

The pilot ends with 13 more plates, 7.64 percentage points more land, and thicker mean Hc/Hm.
Several material ledgers also move substantially: the pilot records 8.94×10¹⁷ m³ formed and
no merge-capped or subducted volume, versus 3.35×10¹⁷, 9.76×10¹⁶ m³, and 3.79×10¹⁴ m³ in
the baseline respectively. The simultaneous solve converged on all 250 steps (at most eight
contact sweeps); the 303.8 cm/yr peak occurred at 384 Myr, so it is not a solver nonconvergence
artifact. Runtime is within the 10% target for this one density-1 replay,
but the extreme rate excursion and geological/topological drift mean the contact candidate
is **not accepted**. This run proves the algebraic response engages; it does not show that
the selected capacity or the remaining reduced physics produces an acceptable long-run
history. The 50 Myr 0.1/1 Myr sensitivity check is recorded below; the full-horizon timestep
comparisons and late-world density-4 profile remain open. The matched 2 Myr replay and regional
column-GPE magnitude checks are recorded below.

**Peak-rate attribution.** An instrumented replay to 386 Myr found the 303.8 cm/yr maximum
on plate 18 at the 384 Myr motion solve. None of that step's 18 blocked contacts involved
plate 18; the next-fastest plate was at 5.03 cm/yr. Plate 18's pre-motion area-weighted Hc
was 11.85 km (maximum 15.18 km) and Hm was 29.35 km (maximum 31.40 km); no Hc nodes were at
the numerical cap. Thus the peak is not a blocked-contact yield response or an over-thickened
collision belt. It is a whole-plate torque-balance rate on a plate outside those blocked edges,
exposed by removing the baseline's 15 cm/yr clamp. At the 2 Myr step size, that angular rate
corresponds to about 6,076 km of surface arc at the maximum-speed latitude. An instrumented
replay found that plate 18 was not born as a fragment: it was one of the 19 generated plates,
initially 1,131 nodes and 18.75 million km². Subduction steadily consumed it; it had 156 nodes
at 300 Myr, 31 at 370 Myr, 9 at 380 Myr, then four at 384 Myr. At the peak solve it was exactly
a connected 2×2 block of level-0 cells on one cube-sphere face (each node had two edge-sharing
neighbors), covering 61,133 km² (0.012% of the planet). It had zero nodes inside another plate's
territory. Its two plates inside the 3-spacing force reach were 2 and 25; all four of its nodes
were in reach, two convergent and subducting and two divergent. The slab-pull torque norm was
1.67×10²⁴ N·m, ridge-push 1.58×10²¹ N·m, basal-flow 1.80×10²³ N·m, and assembled resistance
1.49×10²⁴ N·m. These are vector norms and should not be added as scalars. The specific hypothesis
that plate 18 was overlapped and therefore received no boundary load is falsified: it was not
overlapped and its boundary loads were assembled. Plate 18 disappeared at the end of the
following 386 Myr step, when the model logged that it had no land left and accreted its last
crust onto plate 2. The more concerning failure is that this nearly consumed, 2×2 remnant was
still treated as a freely moving rigid plate for one more motion solve after the 15 cm/yr clamp
was removed; it did not need to be disconnected or overlapped to evade cleanup. Neighbor count
here means plates within the solver's 3-spacing force reach, not polygon sides or a complete
topological neighbor list. This replay has no pressure field, so local convergence stress cannot
be inferred from absolute plate speed alone.

## Four-seed matched G2 replay (2 Myr, 500 Myr)

The G2 working-tree pilot was compared against the detached baseline at `2279d7e16c7d`
(PR #473) on all four G1 seeds. Each run generated its initial world with
`generate_world(seed=..., node_density=1)`, then advanced 250 complete steps of 2 Myr to
500 Myr. Checkpoints were sampled every 50 Myr. Baseline and pilot ran sequentially, paired
by seed, with a Numba cache under `/private/tmp` to keep the detached baseline read-only.
The pilot is the uncommitted G2 worktree based at the same `2279d7e16c7d` commit; exact
source hashes and run metadata are in [`g2_four_seed_2myr.json`](g2_four_seed_2myr.json).
Full per-step outputs remain in `/private/tmp/rfc0012-g2-four-seed/`.

Command shape from each target worktree (the probe script lives in the pilot worktree and
imports `app` through that run's `PYTHONPATH`):

```bash
NUMBA_CACHE_DIR=/private/tmp/rfc0012-g2-four-seed-numba \
PYTHONPATH=backend backend/.venv/bin/python \
  /Users/adubey/.codex/worktrees/rfc0012-g2/Mantle\ Bloom/bin/debug/probe_physics_first_g2.py \
  --seed <seed> --steps 250 --step-myr 2 --node-density 1 \
  --snapshot-every 25 --out /private/tmp/rfc0012-g2-four-seed/<model>-seed<seed>-2myr.json
```

| Seed | Land at start | Baseline land at 500 Myr | Pilot land at 500 Myr | Change | Mean Hc km (base → pilot) | Mean Hm km (base → pilot) | Plates (base → pilot) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 326494039 | 27.11% | 14.36% | 21.99% | +7.64 pp | 18.11 → 19.51 | 49.58 → 55.10 | 33 → 46 |
| 748935664 | 12.14% | 8.80% | 13.58% | +4.78 pp | 14.32 → 14.28 | 45.38 → 48.32 | 46 → 45 |
| 580014012 | 26.54% | 21.57% | 25.91% | +4.34 pp | 22.23 → 21.78 | 65.36 → 72.73 | 35 → 36 |
| 111065103 | 40.90% | 25.44% | 32.64% | +7.20 pp | 22.21 → 25.03 | 56.74 → 64.31 | 40 → 39 |

The land increase is positive on all four seeds, averaging **+5.99 percentage points**.
The starting land fraction is identical within each matched pair. Mean mantle-lithosphere
thickness rises in every pilot; mean crustal thickness rises on two seeds and is nearly
unchanged or modestly lower on the other two.

| Seed | Peak rate, base → pilot | Pilot plate at peak (Myr) | Nodes; area | Neighbors in force reach | Overlap nodes | Median step-time change |
|---:|---:|---:|---:|---:|---:|---:|
| 326494039 | 15.0 → 303.8 cm/yr | 18 (384) | 4; 61,133 km² | 2 | 0 | +10.0% |
| 748935664 | 15.0 → 119.2 cm/yr | 11 (494) | 7; 111,108 km² | 2 | 0 | +36.2% |
| 580014012 | 15.0 → 60.5 cm/yr | 33 (486) | 11; 162,663 km² | 4 | 0 | +57.9% |
| 111065103 | 15.0 → 11.0 cm/yr | 8 (444) | 52; 713,943 km² | 4 | 0 | −19.4% |

The three extreme pilot peaks occur on very small, non-overlapping plates, so the plate-18
rate excursion is not unique to one seed. Their recorded force-reach neighbors and loaded
boundary nodes show that the comparison did not reproduce the earlier overlap hypothesis.
The other seed's maximum remains below 15 cm/yr. These data make remnant handling a recurring
follow-up, while leaving its implementation for later as planned. The peak diagnostics do not
by themselves establish how each small plate formed.

All four pilot solves converged on all 250 steps and reached the 1.5×10¹³ N/m contact-force
capacity; yielded-contact totals were 1,382, 393, 971 and 638 in seed order. The candidate's
median step time changed by +10%, +36%, +58% and −19% at **density 1**. This is not the RFC's
required density-4 runtime measurement, and the full-world medians include the effects of
changed geology and climate as well as solver overhead.

The land result is not yet an acceptance result. Pilot `formed_m3` is 1.4–2.7× baseline, while
`subducted_m3` and `merge_capped_m3` are zero in every pilot; both were nonzero in every
baseline. That topology and ledger drift needs explanation before attributing the extra land
to a validated force response. The 50 Myr 0.1 and 1 Myr timestep checks and regional
ridge/collision ΔGPE scale check are recorded below. Full-horizon small-step comparisons and
density-4 runtime remain open G2 gates. These four-seed 500 Myr replays ran at 2 Myr only.

## G2 regional column-GPE scale check

Reproducible calculation: [`bin/debug/analyze_rfc0012_gpe.py`](../../bin/debug/analyze_rfc0012_gpe.py); compact inputs and outputs: [`g2_gpe_analysis.json`](g2_gpe_analysis.json). Run with `backend/.venv/bin/python bin/debug/analyze_rfc0012_gpe.py`. This is a first-order, geometry-defined magnitude screen, not a fitted column or a full thermal/isostatic inversion.

For a ridge, I used the column thermal-buoyancy approximation
`ΔGPE = ½ g Δρ H²`, plus the RFC's relief term `½ g (ρ_a − ρ_w) h²` (3,330 kg/m³ asthenosphere and 1,030 kg/m³ water). The selected ridge-to-old-flank depth difference is 3.32 km, consistent in scale with the age-depth trend in Parsons & Sclater (1977). Three explicit thermal scenarios span effective mantle density deficits of 35–65 kg/m³ over 90–120 km. Their force-per-length results are 1.51, 2.83, and 4.72 × 10¹² N/m. For the stated 1,000 km ridge segment, the corresponding net forces are 1.51, 2.83, and 4.72 × 10¹⁸ N; dividing by segment length returns the displayed N/m values. The middle scenario is near the Copley et al. (2010) order-of-magnitude ridge-push estimate of 2.5 × 10¹² N/m, while the low/high geometries bracket it.

For collision, the reference column is 35 km of crust (20 km upper plus 15 km lower), 0.5 km surface elevation, and mantle to a common 125 km compensation depth. Tibet columns use 60, 65, and 70 km crust; their surface elevation is solved from equal areal mass (Airy compensation), yielding 3.58–4.16 km for the density cases. Upper/lower crust and mantle densities are separately stated in the JSON. The resulting Tibet-minus-India GPE contrasts are 4.00, 5.07, and 5.51 × 10¹² N/m. Thus the density/geometric sweep overlaps the Copley et al. 5–6 × 10¹² N/m reference in its central and high cases, with the low case below it. The 60–70 km Tibet crust range is consistent with receiver-function measurements that report roughly 65–70 km across parts of the plateau; it is not a claim of uniform Tibet structure.

**Interpretation:** both checks land at the right order of magnitude without selecting a single column to force a match. This supports continuing G2 with a consistent column-GPE definition, but it does not validate the pilot's force assembly or demonstrate that the simulation's current relief proxy represents these loads. The ridge result is especially sensitive to effective thermal density deficit and anomaly depth (`H²`); the collision result depends on crustal thickness, density layering, compensation depth, and the Airy assumption. These remain explicit model assumptions for the next physical implementation/validation step.

Primary references: Copley, Avouac & Royer (2010), [India–Asia collision and the Cenozoic slowdown of the Indian plate](https://doi.org/10.1029/2009JB006634); Parsons & Sclater (1977), [ocean-floor bathymetry and heat flow with age](https://doi.org/10.1029/JB082i005p00803); Molnar & Lyon-Caen (1988), [simple physical aspects of mountain-belt support and evolution](https://doi.org/10.1130/SPE218-p179); receiver-function crustal constraints, [Yue et al. (2012)](https://doi.org/10.1029/2011JB008545).

## Four-seed timestep check (50 Myr)

All four RFC seeds were generated at density 1 from the matched origin/main baseline commit `2279d7e16c7d` and run to 50 Myr at 1 Myr and 0.1 Myr steps. Each baseline/pilot pair used the same seed; the 16 replays ran in series, with one simulation process active at a time. All completed successfully. Per-step records and the run manifest are in `/private/tmp/rfc0012-g2-timestep-50myr/`; compact results are in [`g2_timestep_50myr.json`](g2_timestep_50myr.json). The pilot was the same uncommitted G2 working tree used for the 2 Myr / 500 Myr comparison.

| Seed | Baseline land % (1 → 0.1 Myr) | Pilot land % (1 → 0.1 Myr) | Pilot − baseline land, pp (1 / 0.1 Myr) | Pilot peak rate cm/yr (1 → 0.1 Myr) | Pilot median step s (1 → 0.1 Myr) |
|---:|---:|---:|---:|---:|---:|
| 326494039 | 25.07 → 25.20 | 26.36 → 26.98 | +1.29 / +1.78 | 7.97 → 7.97 | 2.13 → 2.23 |
| 748935664 | 12.40 → 12.27 | 12.55 → 11.98 | +0.16 / −0.29 | 7.47 → 7.61 | 2.04 → 2.12 |
| 580014012 | 25.61 → 24.74 | 25.55 → 24.89 | −0.06 / +0.15 | 4.99 → 6.38 | 1.96 → 1.97 |
| 111065103 | 36.10 → 36.13 | 37.88 → 37.83 | +1.78 / +1.69 | 9.37 → 9.37 | 2.31 → 2.35 |

**Result:** at this 50 Myr horizon, baseline mean land changes by −0.21 percentage points from 1 to 0.1 Myr, and pilot mean land changes by −0.17 points. The pilot's mean land gain over baseline is similar at the two step sizes: +0.79 points at 1 Myr and +0.83 points at 0.1 Myr. Peak pilot rates remain 4.99–9.37 cm/yr, with no pilot step exceeding 15 cm/yr; the baseline's 15 cm/yr cap binds on seed 111065103 at both step sizes. Pilot solvers converged on every step, and their maximum relative torque residuals remain in the compact JSON.

The median-of-seed median step time is 2.08 s for the pilot at 1 Myr and 2.17 s at 0.1 Myr; the corresponding baseline values are 1.95 s and 2.06 s. That is +7.0% and +5.5% for the pilot at density 1, below the 10% target in this early sample; it does not replace the required density-4 profile. Contact activity is accumulated over 50 steps per 1 Myr run and 500 per 0.1 Myr run, so compare rates per step rather than raw totals. Across seeds, pilot active blocked contacts average 1.23 per step at 1 Myr and 1.51 at 0.1 Myr; yielded contacts average 0.69 and 0.80 per step respectively. The maximum relative torque residual across pilot steps is 7.42×10⁻¹⁵; no pilot solve failed to converge, and the contact force reached but did not exceed the configured 1.5×10¹³ N/m cap.

This is a useful early-horizon timestep check, not the full 500 Myr timestep-sensitivity gate. It shows stable land fraction and peak rates over 50 Myr, but it does not test the late remnant events seen in the 2 Myr / 500 Myr replay, nor the RFC's density-4 runtime target. The previously observed long-run topology and material-ledger drift also remains unresolved.

## G2 matched glacial-erosion sensitivity replay (250 Myr)

The 1 Myr G2 run was repeated with `glacier_erosion_multiplier` increased from 1 to 2,
holding seed, G2 code, timestep, and craton-formation rule fixed. The multiplier scales
glacial abrasion and the glacier-driven flattening pass; it does **not** change ice
accumulation or melt rates, so this is an erosion-strength sensitivity test rather than a
glacier mass-balance test. Four seeds were paired and run serially. The 200 Myr `formed_m3`
increment captures the initial cohort pulse; `eligible_host` is measured just before it, at
190 Myr. Detailed checkpoint traces are under `/private/tmp/g2-glacial-sensitivity/`, and
compact matched values are in [`g2_glacial_erosion_sensitivity.json`](g2_glacial_erosion_sensitivity.json).

| Seed | Applied glacier erosion by 200 Myr, 1× → 2× (M km³) | Craton formed at 200 Myr, 1× → 2× (M km³) | Craton formed by 250 Myr, 1× → 2× (M km³) |
|---:|---:|---:|---:|
| 326494039 | 8.47 → 17.79 | 737.44 → 736.81 | 791.25 → 772.85 |
| 748935664 | 9.87 → 11.41 | 274.98 → 266.67 | 314.61 → 284.95 |
| 580014012 | 54.03 → 117.77 | 644.21 → 698.98 | 702.63 → 738.76 |
| 111065103 | 11.47 → 20.65 | 853.43 → 927.68 | 915.97 → 994.31 |
| **Sum** | **83.84 → 167.63** | **2,510.06 → 2,630.14** | **2,724.46 → 2,790.87** |

Doubling the control approximately doubled cumulative applied glacial erosion by 200 Myr
across the four seeds. The formation response is mixed: it falls 0–3% on two seeds and rises
about 9% on the other two. Across all four, the pulse is 4.8% larger at 2×, and cumulative
formed craton by 250 Myr is 2.4% larger. The strongest erosion seed (580014012) does not form
less craton at 2×: its pulse increases by 8.5%. In two seeds, the eligible host pool at 190 Myr
also grows under 2×, while eligible-host mean Hc drops modestly in all four. Thus stronger
glacial erosion affects the evolving host and topography, but it is not a consistent brake on
the 200 Myr cohort pulse and does not explain the extra G2 craton with one directional effect.

The local-relief proxy (neighbor elevation range) fell on three seeds and rose on one; the
glaciated area above 2 km and the cumulative erosion carried by currently identified mountain
nodes also changed by seed. These are cross-sectional summaries, not same-cell causal tracks;
the optional lineage follows stable node IDs and cannot preserve history for nodes replaced by
topology remaps. Seed 580014012 also repeatedly hit the 5 km ice-depth cap, with warnings
reaching over 1,000 capped nodes and one reported shed event around 449,000 km³ in its 2× run.
That cap behavior limits interpretation of the high-ice result.

**Conclusion:** this replay does not support “glaciers are too weak” as the sole explanation
for G2's extra early craton. It only tests erosion strength, not ice mass balance, and the
formation rule still grants craton after a fixed 200 Myr in the geometric/Hc gate. A separate
formation-gate/timer sensitivity remains necessary to distinguish the synchronized
generation-time cohort from later geological evolution. The reproduction script is
[`g2_glacial_erosion_sensitivity.py`](../../bin/debug/g2_glacial_erosion_sensitivity.py); the
per-node erosion trace used a temporary instrumentation overlay, so its full historical
traces remain in `/private/tmp` rather than the repository.

## G2 craton-formation timer sensitivity (400 Myr)

This matched replay varies only `CRATON_FORMATION_MYR` across 100, 200, and 300 Myr. The
same four seeds, 1 Myr step, and G2 source were used; each run was serial. The existing host
gate stayed fixed (genuine continental crust, Hc ≤ 50 km, at least 500 km inside the plate
edge). Thus this tests how long a qualifying host must wait, not whether “quiet” should be
defined using physical strain or deformation. The runtime timer override is in
[`g2_formation_gate_sensitivity.py`](../../bin/debug/g2_formation_gate_sensitivity.py), and
compact results are in [`g2_formation_gate_sensitivity.json`](g2_formation_gate_sensitivity.json).
Full traces are in `/private/tmp/g2-formation-gate-sensitivity/`.

| Seed | Cumulative formed at 400 Myr (M km³), 100 / 200 / 300 Myr timer | Live craton at 400 Myr (M km³), 100 / 200 / 300 Myr timer | Eligible uncratonized host area at 400 Myr (M km²), 100 / 200 / 300 Myr timer |
|---:|---:|---:|---:|
| 326494039 | 1,176 / 886 / 662 | 3,768 / 3,504 / 3,289 | 1.09 / 1.62 / 3.48 |
| 748935664 | 472 / 362 / 280 | 2,004 / 1,900 / 1,828 | 1.92 / 2.97 / 6.16 |
| 580014012 | 1,241 / 854 / 592 | 4,348 / 3,967 / 3,722 | 2.57 / 5.85 / 7.37 |
| 111065103 | 1,675 / 1,131 / 796 | 4,742 / 4,219 / 3,903 | 4.67 / 9.06 / 12.82 |
| **Sum** | **4,564 / 3,233 / 2,330** | **14,862 / 13,590 / 12,741** | **10.25 / 19.51 / 29.84** |

All four seeds show a monotonic decrease in cumulative formation as the timer lengthens.
Summed cumulative formation at 400 Myr falls 29% with a 200 Myr timer and 49% with a 300 Myr
timer relative to 100 Myr. The summed live-craton stock falls less (9% and 14%), because the
stock includes the starting craton and subsequent losses, while the cumulative `formed_m3`
ledger records all formation events. Longer timers leave more qualifying but uncratonized
area at the endpoint, as expected. In every run the first formation pulse occurs exactly at
the selected timer threshold, directly confirming that the formation timer controls the
cohort timing.

**Interpretation:** the fixed 200 Myr waiting period is a substantial control on the amount
of craton formed during the first 400 Myr. This supports the hypothesis that early
timer-triggered formation contributes to the excess-craton result; it does not establish
that 100 or 300 Myr is geologically preferable, nor that the current host gate represents
low deformation. The next diagnostic should record the qualifying hosts' deformation/strain
history and test a physically defined quiet gate against the current clock. Seed 580014012
repeatedly hit the 5 km glacier-depth cap, with hundreds to over 1,000 affected nodes on some
steps and reported ice shedding/clipping; this limits interpretation of its absolute result
and remains an independent glacial-model concern.


## Four-seed full-horizon timestep check (1 Myr, 500 Myr)

The four G1 seeds were run to 500 Myr at 1 Myr steps using density 1, with each baseline/pilot pair run serially. The baseline is the detached `2279d7e16c7d` tree; the pilot is the same uncommitted G2 working tree used above. All eight replays completed successfully. Per-step outputs and run manifest are in `/private/tmp/rfc0012-g2-full-1myr/`; the compact matched data are in [`g2_timestep_500myr_1myr.json`](g2_timestep_500myr_1myr.json).

| Seed | Baseline land at 500 Myr | Pilot land at 500 Myr | Pilot − baseline | Pilot peak rate | Pilot subducted (million km³) | Pilot formed (million km³) | Pilot merge-capped (million km³) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 326494039 | 15.08% | 21.36% | +6.28 pp | 66.07 cm/yr | 0.0000 | 901.2 | 0.0000 |
| 748935664 | 10.14% | 13.18% | +3.04 pp | 152.59 cm/yr | 0.4667 | 383.8 | 0.0000 |
| 580014012 | 23.30% | 23.22% | -0.08 pp | 4.99 cm/yr | 0.0000 | 924.0 | 0.0000 |
| 111065103 | 25.43% | 32.19% | +6.76 pp | 9.37 cm/yr | 0.0000 | 1255.6 | 0.0000 |

Pilot land is higher in three seeds and effectively unchanged in one, with mean final change **+4.00 percentage points**. Compared with the 50 Myr check, the land response grows during the full horizon. This is not yet a successful G2 acceptance result: pilot peak rates still reach 152.6 cm/yr on seed 748935664 and 66.1 cm/yr on seed 326494039, while the baseline is clipped at 15 cm/yr. Pilot subduction is zero on three seeds (nonzero only on 748935664), and pilot merge-capped volume is zero on all four; baseline merge-capped volume is nonzero on three. Pilot formed volume is 1.2–1.9× baseline. These material-ledger/topology differences remain the main unresolved concern.

The pilot peak rate varies sharply by seed: 326494039: 66.1, 748935664: 152.6, 580014012: 5.0, 111065103: 9.4 cm/yr. The smaller 1 Myr step avoids the 303.8 cm/yr event seen at 2 Myr for seed 326494039, but does not establish timestep convergence: seed 748935664 reaches 152.6 cm/yr. Full-horizon 0.1 Myr runs remain the next convergence check if needed.

Seed 580014012 produced extensive 5 km glacier-cap and ice-shedding/clipping warnings in both baseline and pilot; the pilot warnings affected hundreds to over 1,000 nodes on some steps. Treat its climate and land-fraction result with caution. Density-4 runtime remains unmeasured; these density-1 runs had mean median step times of 2.763 s (baseline) and 2.790 s (pilot).

### Full-horizon timestep comparison (pre-registered)

**Question:** do the 500 Myr baseline and G2 outcomes remain comparable when the step shrinks
from 1 Myr to 0.1 Myr, beyond the existing 50 Myr sensitivity check?

**Expected direction and scale:** based on the 50 Myr check, endpoint land fraction should change
by less than 2 percentage points per seed and model, and G2's mean land-fraction advantage
should remain positive and within 2 percentage points of its 1 Myr value. Long-horizon topology
may amplify small timestep differences, so any larger change will be treated as a sensitivity
finding, not averaged away.

**Primary measures and screen:** per-seed endpoint land fraction for baseline and G2, and the
four-seed mean G2-minus-baseline land difference at each timestep. The land response passes this
screen only if each model/seed's 0.1-versus-1 Myr endpoint shift is ≤2 percentage points and
the mean G2 effect has the same sign at both steps with a difference ≤2 points. Report peak
rate, plate count, Hc/Hm, step runtime, and craton/subduction/merge ledgers as secondary outcomes;
these are descriptive, with no preselected convergence threshold.

**Controls:** generated worlds at node density 1.0; seeds 326494039, 748935664, 580014012, and
111065103; 500 Myr endpoint; baseline code `2279d7e16c7d`; G2 model code `8618185b59c5` (no
backend changes since that commit). Reuse the existing full-horizon 1 Myr traces after verifying
their inputs, and run only the eight missing 0.1 Myr cases. Run baseline then G2 within each
seed, with only one simulation process active at a time. Save full outputs and the manifest
outside the repository; check in only a compact comparison JSON and the result summary.


## Four-seed G2 ledger and quiet-gate traces (1 Myr, 500 Myr)

To determine whether zero subduction and merge-cap entries represented missing ledger events or a real change in model behavior, the same four seeds were replayed serially on baseline `2279d7e16c7d` and G2 pilot `512895dc2f3a`. Each run used `generate_world(seed=..., node_density=1)`, 500 steps of 1 Myr. The runtime wrappers in [`trace_rfc0012_ledger_gates.py`](../../bin/debug/trace_rfc0012_ledger_gates.py) record subduction-front candidates and ledger deltas, boundary deformation, strain-only lithosphere changes, suture transfers, merge inventories and caps, formation gate clocks/reset categories, erosion and ice unloading, and per-step material balance. [`summarize_rfc0012_ledger_gates.py`](../../bin/debug/summarize_rfc0012_ledger_gates.py) checks event totals against endpoint ledgers and writes the compact matched summary in [`g2_ledger_gate_traces.json`](g2_ledger_gate_traces.json). The manifest and full step-by-step JSONL traces remain in `/private/tmp/rfc0012-ledger-gates/`.

Volumes below are million km³. The subduction columns show total deeply subducted continental material and the craton subset. “Merge cap” is craton recorded as capped during plate merges.

| Seed | Merges, baseline → pilot | Deep subduction, baseline → pilot | Craton subducted, baseline → pilot | Merge cap, baseline → pilot | Craton formed, baseline → pilot | Land fraction, baseline → pilot |
|---:|---:|---:|---:|---:|---:|---:|
| 326494039 | 21 → 0 | 2,040.3 → 639.6 | 0.306 → 0.000 | 18.55 → 0.00 | 472.3 → 856.0 | 15.08% → 21.18% |
| 748935664 | 13 → 1 | 966.2 → 424.5 | 0.130 → 0.010 | 8.03 → 0.00 | 212.9 → 353.9 | 10.14% → 14.16% |
| 580014012 | 5 → 1 | 1,945.5 → 546.3 | 0.733 → 0.006 | 0.00 → 0.00 | 771.6 → 1,030.2 | 23.30% → 24.31% |
| 111065103 | 15 → 0 | 2,764.6 → 829.5 | 0.520 → 0.000 | 72.66 → 0.00 | 812.8 → 1,251.0 | 25.43% → 31.46% |

For all eight runs, deep-subduction events, merge-cap events, and formation events each sum exactly to their corresponding endpoint ledger. The residual continental-material and craton balance errors are at floating-point scale (at most 14,336 m³ and 1,536 m³ respectively, against inventories measured in roughly 10²⁰–10²² m³). Thus the zero pilot entries are not unbooked losses: these runs produced much less deep subduction and almost no merges. Pilot merge-cap volume is zero because the pilot made only zero or one merge per seed, compared with five to 21 in baseline. This supports the proposed explanation that slower collision dynamics reduce merge-trigger exposure; the traces establish the event-count and ledger changes, but do not isolate collision braking as the sole cause. The pilot also formed more craton in every seed, consistent with its lower recorded craton subduction and fewer opportunities for merge-cap loss. Land fraction rose in all four traced runs.

The quiet-gate trace records gate eligibility and reset categories at every plate step, plus aggregate absolute Hc/Hm changes in hosts that qualified as quiet. Across the four pilot runs, 425–1,499 of 12,944–17,158 plate-step samples had a local Hc or Hm change while the quiet-host mask applied. This indicates the current geometric/clock gate does not guarantee lithosphere stability. These are aggregate masked changes, not persistent cell identities through topology remaps; a future test of an endogenous strain-based quiet gate needs lineage-aware host histories. The trace also records glacier unloading area and maximum ice-deflection change. Seed 580014012 repeatedly hit the 5 km ice cap in both baseline and pilot, so its erosion/unloading diagnostics remain affected by the known glacier clipping limitation.

Reproduce the runs with the manifest's recorded commands and commit roots. Keep the eight runs serial to limit memory use; use `/private/tmp/rfc0012-ledger-gates/` as the output directory. This trace pass supplies the recorded evidence for follow-up issues #493 (early craton/quiet-gate diagnosis) and #494 (high-rate stress, heating, and melt response); it does not implement either physical model change.

## Four-seed G2 merge-opportunity diagnostic (1 Myr, 500 Myr)

To explain the lower merge count in the ledger comparison, the same four seeds were traced
serially on baseline `2279d7e16c7d` and the G2 pilot (`329fd2a5f03b`, containing model change
`9e7f9e1`). The detector records each continental pair's filter stage, close-contact count,
contact closing rates, timer progress and probability result; forced-overlap progress and the
selected merge path are also recorded. The instrumented candidate detector was checked against
the production detector on generated worlds, and the timer replay matched the production ready
list on all eight runs. Per-step traces are in `/private/tmp/rfc0012-merge-opportunities/`;
compact results are in [`g2_merge_opportunity_diagnostic.json`](g2_merge_opportunity_diagnostic.json).
Reproduce with [`trace_rfc0012_merge_opportunities.py`](../../bin/debug/trace_rfc0012_merge_opportunities.py)
and [`summarize_rfc0012_merge_opportunities.py`](../../bin/debug/summarize_rfc0012_merge_opportunities.py).

| Measure across four seeds | Baseline | G2 pilot |
|---|---:|---:|
| Continental pair checks | 196,272 | 238,071 |
| Checks with at least 3 close contact nodes | 16,805 (8.6%) | 12,665 (5.3%) |
| Convergent collision candidates | 1,001 (0.51% of checks) | 78 (0.033% of checks) |
| Collision pair-steps tracked by the persistence timer | 3,835 | 161 |
| Timer-ready pair-steps | 77 | 3 |
| Sustained-collision merges | 48 | 1 |
| Forced deep-overlap merges | 6 | 1 |
| Total merges | 54 | 2 |

The loss is not explained by the broad-phase search pruning away more pairs: the pilot checks
21% more continental pairs, and its bounding-sphere survivor fraction is similar (29.8% versus
28.4%). The decline appears at close-node qualification and especially convergence. Among
close-contact pair-steps with at least three close nodes, the share with at least three nodes closing faster than
the 1 cm/yr transform threshold falls from 6.0% to 0.6%. For rejected close pairs, median
maximum local closing speed falls from 0.58 to 0.26 cm/yr, and median relative plate rotation
falls from 1.20 to 0.43 cm/yr. The rare pilot candidate pair-steps have a similar median maximum
contact closing speed to baseline candidates (1.66 versus 1.46 cm/yr). The persistence timer
therefore sees far fewer pair-steps and rarely reaches its threshold; it is not dropping a
large pool of already-qualified pairs. Timer probability is a downstream effect too: only one
pilot pair-step passed the roll, compared with 48 baseline passes.

This supports the hypothesis that G2 collision resistance reduces sustained merge exposure,
which is consistent with fewer merges and less merge-capped craton. It does not isolate the
force law as the sole cause: the pilot's evolving topology and contact geometry also differ.
Candidate qualification requires enough spatially close nodes **and** at least three of those
nodes converging faster than 1 cm/yr, so a global plate peak speed is not itself a merge
opportunity. The diagnostic found no timer-ready decision mismatches in any run. Seed
580014012 again emitted the known 5 km glacier-cap warnings; they are unrelated to the pair
detector but remain a caveat for geological outcomes.

## Four-seed craton-ledger divergence and attribution

This analysis uses the existing eight 1 Myr / 500 Myr event traces; it does not rerun the
simulations. [`analyze_rfc0012_craton_ledger_divergence.py`](../../bin/debug/analyze_rfc0012_craton_ledger_divergence.py)
reconstructs live craton at each step as initial inventory + cumulative formation − all declared
craton sinks + the recorded balance error. The compact output has 25 Myr checkpoints, exact
1 Myr threshold-crossing times, and a final attribution for every ledger account:
[`g2_craton_ledger_divergence.json`](g2_craton_ledger_divergence.json). All values below are
million km³.

Across the four seeds, reconstructed pilot live craton exceeds baseline by **1,536.0 million
km³** at 500 Myr. The largest source of that difference is **1,221.4 million km³ of additional
formation** (about 80% of the net gap). Reduced rifting accounts for another **259.6 million
km³** of retained craton, and reduced merge-cap loss for **99.2 million km³**. Smaller terms
include 10.4 million km³ less erosion and 1.7 million km³ less craton subduction. These gains
are partly offset by 42.4 million km³ more craton reworked into ordinary crust at collision
zones. The residual after accounting for all sources and sinks is at floating-point scale.

| Seed | Initial inventory Δ | Live craton Δ at 200 Myr | Live craton Δ at 500 Myr | Additional formation | Avoided rifting | Avoided merge cap |
|---:|---:|---:|---:|---:|---:|---:|
| 326494039 | +0.7 | +307.6 | +479.2 | +383.6 | +105.9 | +18.6 |
| 748935664 | +0.7 | +79.1 | +186.8 | +141.0 | +47.3 | +8.0 |
| 580014012 | −1.7 | +324.2 | +291.3 | +258.6 | +24.8 | 0.0 |
| 111065103 | −13.6 | +348.7 | +578.7 | +438.2 | +81.6 | +72.7 |
| **Sum** | **−13.9** | **+1,059.6** | **+1,536.0** | **+1,221.4** | **+259.6** | **+99.2** |

The aggregate live-stock difference is only 17.3 million km³ at 175 Myr, then jumps to 1,059.6
million km³ at 200 Myr. All four seeds record their first formation at exactly 200 Myr, and the
pilot-minus-baseline formation difference at that checkpoint is already **1,033.5 million
km³**. It grows to 1,221.4 million km³ by 500 Myr, while avoided rift loss and merge capping
continue to widen the retained-stock difference. The 200 Myr pulse therefore explains most
of when the trajectories separate; later source and sink differences explain how the gap
continues to evolve.

The reduction in `rifted_m3` is mainly the stretch path: pilot-minus-baseline rift-path ledger
deltas sum to −11.3 million km³ for divergent thinning, −215.5 for stretch donation, and −32.7
for failed rifts. This is a new material contributor alongside merge capping and needs physical
interpretation; the merge-opportunity result alone cannot explain it. Importantly, the large
drop in **total continental deep subduction** is not the same as craton subduction: craton
subduction accounts for only 1.7 million km³ of the four-seed retained-stock difference.

The runs are seed-matched, but their recorded starting `initial_m3` inventories are not
identical: pilot minus baseline ranges from −13.6 to +0.7 million km³ (−13.9 summed). The
attribution includes this term; it is under 1% of the aggregate final gap, though it is another
reason to avoid treating these independently generated worlds as bit-for-bit paired controls.

### Event-level checks for rifting and merge-cap losses

The event audit reuses the same traces. For each merge, it compares the merge event's cap
amount with the step ledger and checks that live craton falls by that amount. Across all runs,
the sum of per-merge cap records matches the step-ledger `merge_capped_m3` exactly. Nine merges
have material cap losses (above 0.001 million km³); they account for the full 99.2 million km³
baseline-versus-pilot difference. The pilot has no positive merge-cap events. Per-event changes
in live stock plus cap are at most 358 m³, rounding error against inventories around 10²¹ m³.
The zero pilot merge-cap account is therefore a real absence of cap-producing merges, not a
missing ledger write.

The three rift-path accounts also reconcile to `rifted_m3` in every run. Across the four seeds,
divergent thinning falls from 14.22 to 2.94 million km³, stretch donation from 288.63 to 73.09,
and failed-rift thinning from 63.81 to 31.08. New event hooks attribute each loss to its
operation and source plate, with affected-cell visits and a cellwise inventory cross-check.
Event totals match each step ledger exactly; the maximum absolute cellwise residual is 32 m³.
Across the runs, the hooks recorded 2,784 baseline versus 364 pilot divergent-thinning
operations, 2,417 versus 1,059 stretch-donation operations, and 64 versus 49 failed-rift
operations. These are operation records, not unique geological rift episodes; cell counts are
visits across events, not unique cells.

The per-seed attribution shows the aggregate reduction is not uniform across mechanisms.
Failed-rift loss is lower in three pilot seeds, but increases from 1.42 to 5.56 million km³ in
seed 580014012 despite fewer affected-cell visits overall across all runs. In that seed,
failed-rift operations rise from 4 to 12, and the share of operations cutting through craton
rises from 4/4 to 12/12. That is evidence that where cuts land matters as well as how often the
operation runs. Stretch donation remains the largest change: the pilot has 73.09 versus 288.63
million km³ across the four runs, with fewer operation records and affected-cell visits.
Baseline plate 4 is the largest stretch-loss contributor in seeds 326494039 and 748935664, but
plate IDs cannot be compared between independently evolved topologies or treated as persistent
identities.

This trace establishes where ledgered loss was applied and validates the amount, but not the
causal chain from G2's force law to changed rift geometry or donor inventory. The targeted
seed-580 donor-cell trace below records that pre-operation inventory and the plate/cut state.
The result does not indicate a ledger anomaly. The compact event audit is part of
[`g2_craton_ledger_divergence.json`](g2_craton_ledger_divergence.json); full JSONL traces remain
outside the repository under `/private/tmp/rfc0012-rift-attribution/` and
`/private/tmp/rfc0012-rift-donor-detail/`.

## Seed 580014012 rift donor-cell and cut-state diagnostic

To explain why failed-rift craton loss increased in this seed even though the number of failed
cuts did not, the seed was replayed to 500 Myr at 1 Myr steps on baseline `2279d7e16c7d` and
G2 pilot `d2d7d79c1f4a`. The pilot's simulation code is unchanged from the previously traced
candidate; this revision adds instrumentation only. The trace records donor-cell identities,
local coordinates, pre/post craton thickness and area for stretch donation and failed-rift
thinning. Failed cuts also record the cut band inventory, cut normal, plate topology revision,
plate size/age/stress/frame/rotation, neighboring plates, and observed split-flow residual and
fitted-pole separation against their gate thresholds. The runtime tracer is
[`trace_rfc0012_ledger_gates.py`](../../bin/debug/trace_rfc0012_ledger_gates.py); the summary
tool is [`analyze_rfc0012_rift_donor_detail.py`](../../bin/debug/analyze_rfc0012_rift_donor_detail.py).
The compact summary is [`g2_rift_donor_detail.json`](g2_rift_donor_detail.json), while full
traces remain in `/private/tmp/rfc0012-rift-donor-detail/`.

| Seed 580014012 measure | Baseline | G2 pilot |
|---|---:|---:|
| Failed-rift cuts | 16 | 16 |
| Cuts entering a craton-bearing band | 4 (25%) | 12 (75%) |
| Craton area fraction across failed-rift bands | 5.85% | 21.65% |
| Craton volume in the craton-bearing cut bands | 114.61 million km³ | 446.37 million km³ |
| Failed-rift craton loss | 1.421 million km³ | 5.562 million km³ |
| Loss / craton volume in those cut bands | 1.240% | 1.246% |
| Stretch-donation operations / donor-cell visits | 462 / 2,060 | 113 / 358 |
| Stretch-donation craton loss | 32.887 million km³ | 5.263 million km³ |

The failed-rift loss increase is explained by **greater craton exposure to the cut bands**, not
by stronger thinning per unit of exposed craton. Across the four craton-bearing baseline cuts,
the cut bands held 114.61 million km³ of craton; the twelve pilot cuts held 446.37 million km³
(3.9× as much). Yet both runs removed about 1.24% of that band inventory. The area ledger tells
the same story: the aggregate cut-band area is only 12% larger in the pilot, while its
craton-bearing area is 4.1× larger. The total cut count is identical. G2's trajectory therefore
placed a larger share of failed-rift cuts through craton-bearing plate geometry; the per-cut
thinning response is nearly unchanged.

The step-6 cuts are a useful matched point: both runs fail a cut on plate 10 at similar plate
sizes (4,279 baseline and 4,266 pilot nodes) and lose similar craton volumes (0.575 and 0.562
million km³). Later, the two top pilot losses occur at 43 Myr on plate 3 (2.053 million km³)
and 60 Myr on plate 8 (1.353 million km³), where the cut bands contain 161.7 and 107.1 million
km³ of craton. The pilot's twelve recorded cuts all passed the same split criteria; their
flow-residual ratios range from 1.00 to 2.33 times threshold, and the fitted-pole separation
ratios from 1.09 to 3.17 times threshold. These observations locate the difference in the
selected cut exposure under diverged plate geometry and flow, but do not isolate whether it
comes from the force law changing plate geometry, the split gate, or their interaction.

Stretch donation moved in the opposite direction: the pilot had 75.5% fewer traced stretch
operations and 82.6% fewer donor-cell visits, with an 84.0% lower craton loss. Mean loss per
donor-cell visit is similar (about 0.0160 versus 0.0147 million km³), again pointing to how many
donor cells participate rather than a changed per-cell loss rule. Event losses reconcile to
the rift ledgers; the largest cellwise residual was 2 m³ in baseline and 0.5 m³ in pilot, and
the endpoint craton balance error was −512 m³ in both.

**Interpretation:** the extra failed-rift loss is not a ledger anomaly or a stronger local
thinning response. It is associated with more of the failed cut band's area and craton
inventory overlapping cratonic lithosphere in the G2 trajectory. A follow-up that separates
the force-law effect from the split selector would trace all split cuts (including successful
splits) and compare their cut-band craton exposure and gate metrics on each model trajectory.
Seed 580014012 retains the known glacier-depth-cap warnings, which affect geological
interpretation but not this closed rift ledger.

### Successful split-cut exposure follow-up (pre-registered)

**Question:** does the greater craton exposure in G2's failed-rift cuts also appear among
successful split cuts, or is it specific to the failed-rift outcome? The trace will record each
successful cut's pre-split band area, craton-bearing area and craton volume, plate inventory,
split-gate margins, and daughter inventories. This is a trajectory diagnostic; it will not
change the model or by itself accept/reject G2.

**Expected direction:** unknown. If changed force-law trajectories alter cut selection or
geometry generally, successful-cut exposure should also differ between baseline and G2. If the
failed-rift exposure gap is mostly an outcome-specific sampling effect, successful cuts may not
show the same gap.

**Primary measures:** craton-bearing area fraction and craton volume in the successful-cut
bands; report event count and per-cut distribution alongside aggregate means. A ratio of at least
1.5× in either primary measure is treated as a material exposure difference for this diagnostic;
within 20% is treated as broadly similar for this single-seed screen. These are descriptive
screening cutoffs, not model acceptance criteria. Intermediate or mixed results remain
inconclusive.

**Controls:** seed `580014012`, node density 1.0, 1 Myr steps, 500 Myr horizon, same deterministic
initialization, baseline commit `2279d7e16c7debe839bd04ab18ce0a8f895f8654`, and G2 commit
`8618185b59c5a8b08e3fe85f230373b2184b3552`. Runs will be serial. The prior matched traces
provide the failed-rift comparison; both arms will be replayed with the expanded successful-cut
trace so capture semantics are identical.

**Result:** the 500 Myr serial traces completed in 23.9 minutes baseline and 25.9 minutes G2.
The G2 trace commit is `31778b15bbf54da62b801a447396a205ad92c5fb`; it only adds diagnostics and
analysis, with no backend model changes since G2 commit `8618185b59c5a8b08e3fe85f230373b2184b3552`.
Both arms closed the craton ledger to −512 m³ and the continental-material ledger to 0 m³.
The successful-cut measurements are:

| Measure | Baseline | G2 | G2 / baseline |
|---|---:|---:|---:|
| Successful split cuts | 32 | 26 | — |
| Cuts with craton-bearing area | 14 (43.8%) | 15 (57.7%) | — |
| Craton-bearing area / total cut-band area | 14.21% | 20.02% | 1.41× |
| Craton volume in successful cut bands | 661.77 million km³ | 692.77 million km³ | 1.05× total |
| Craton volume per successful cut | 20.68 million km³ | 26.65 million km³ | 1.29× |

G2 successful cuts show a moderate increase in exposure, but neither primary measure reaches
the preregistered 1.5× materiality screen; both fall in its intermediate range. The shift is
smaller than the failed-cut exposure difference (5.85% versus 21.65% craton-bearing area, or
3.70×). Thus craton exposure is somewhat higher across successful cuts too, while the strong
exposure gap is concentrated in failed cuts. This supports an outcome-conditioned cut-exposure
effect, but it still cannot say whether G2's force law changed plate trajectories, the split
selector chose different cuts, or both. Successful-cut gate margins are broadly similar (mean
flow-fit / threshold 1.11× baseline versus 1.21× G2; mean pole-separation / threshold 1.95× versus
1.92×), so these observed events do not point to a large gate-margin shift as the explanation.

The split-cut trace ledger totals match the successful event count and band areas in both runs.
The compact event data are in [`g2_successful_split_exposure.json`](g2_successful_split_exposure.json);
the full JSONL replays remain outside the repo under `/private/tmp/rfc0012-successful-split/`.

### Full-horizon timestep comparison (pre-registered) — results

**Question:** does replaying the same four matched seeds for 500 Myr at a 0.1 Myr step materially change the 1 Myr-step G2 comparison?

**Controls:** node density 1.0; seeds `326494039`, `748935664`, `580014012`, and `111065103`; baseline model `2279d7e16c7d`; G2 model implementation `8618185b59c5` (the 0.1 Myr trace was recorded at `fa381dadafc1`, with no model-code changes since the G2 implementation). The existing 1 Myr traces were reused. All eight missing 0.1 Myr runs completed serially, with baseline and G2 paired within each seed. Full traces remain outside the repo in `/private/tmp/rfc0012-g2-full-timestep-500/`.

**Pre-registered gates:** (1) for every seed and model, the absolute endpoint land-fraction difference between 0.1 and 1 Myr is at most 2 percentage points; (2) the mean G2-minus-baseline land effect has the same sign at both steps and differs by at most 2 percentage points.

| Seed | Baseline land Δ (0.1 − 1 Myr, pp) | G2 land Δ (0.1 − 1 Myr, pp) | G2 advantage at 1 Myr (pp) | G2 advantage at 0.1 Myr (pp) | Effect Δ (pp) |
|---|---:|---:|---:|---:|---:|
| 326494039 | −1.205 | +0.067 | +6.279 | +7.550 | +1.271 |
| 748935664 | −0.387 | −0.527 | +3.041 | +2.901 | −0.140 |
| 580014012 | −0.278 | +0.197 | −0.085 | +0.390 | +0.475 |
| 111065103 | −0.288 | +0.564 | +6.758 | +7.610 | +0.852 |
| **Mean effect** | — | — | **+3.999** | **+4.613** | **+0.614** |

**Result: both preregistered gates pass.** The largest per-run land shift across timestep choices is 1.205 pp. The mean G2 land advantage remains positive and changes by only 0.614 pp (from +3.999 to +4.613 pp). Seed 580014012 changes from a very small negative G2 effect at 1 Myr (−0.085 pp) to a small positive one at 0.1 Myr (+0.390 pp); the mean-effect gate passes, while this individual seed remains sensitive in sign at a near-zero effect.

Peak rates and ledger values were secondary descriptive measures, not timestep acceptance gates. The fastest G2 outlier remains in seed 748935664 (124.6 cm/yr at 0.1 Myr; 152.6 cm/yr at 1 Myr), so it is not eliminated by the finer step and remains covered by the high-rate follow-up. Formation, subduction and merge ledgers also vary by seed and timestep; this comparison establishes land-outcome timestep robustness, not the causal explanation for G2's ledger changes or validation of the formation/quiet gate.

Compact matched results are in [`g2_timestep_500myr_0.1_vs_1.json`](g2_timestep_500myr_0.1_vs_1.json); the reused 1 Myr results are in [`g2_timestep_500myr_1myr.json`](g2_timestep_500myr_1myr.json).

### Density-4 late-world runtime profile

**Question:** does the G2 coupled motion path keep median total step time within the RFC limit of
10% over baseline at node density 4, after the model has developed a late-world plate count?

**Controls:** both arms were generated from seed `326494039` at density 4, advanced serially to
500 Myr using 2 Myr steps, then warmed by one additional step and measured for five steps.
Baseline is `2279d7e16c7debe839bd04ab18ce0a8f895f8654` (the `origin/main` baseline worktree);
G2 is `fa381dadafc174e82a5a740f6ae8723f10551a63`. Both used Python 3.14.6 and the same host.
The sample is at 502–512 Myr. The run command was:

```sh
MB_SOURCE_ROOT=<checkout> NUMBA_CACHE_DIR=/private/tmp/rfc0012-numba-cache \
  <checkout>/backend/.venv/bin/python /private/tmp/rfc0012_density4_late_profile.py 326494039 500
```

| Measure | Baseline | G2 | G2 vs baseline |
|---|---:|---:|---:|
| Median total step time | 9.008 s | 17.545 s | +94.8% |
| Timed step range (5 steps) | 7.887–10.017 s | 16.091–18.476 s | — |
| Nodes at sample | 134,533 | 129,371 | −3.8% |
| Plates after sample | 41 | 43 | — |
| Peak RSS | 1,407 MiB | 1,246 MiB | −11.5% |

**Initial result: apparent failure; later profile methods conflicted.** The first five-step G2
sample was 1.95× baseline despite fewer nodes. Its phase taps showed gather at 2.008 s/step,
classification at 0.479 s/step, and deformation at 2.565 s/step. A later G2 replay at the same
seed and horizon produced a different late-world cProfile sample: 9.34 s/step averaged over
three profiled steps. cProfile is not an acceptance timing, but this large discrepancy means the
first sample cannot yet establish a stable 94.8% regression. A matched baseline cProfile replay
was then run with the same seed, horizon, warmup and three-step profile window; its paired result
is summarized below. That cProfile pair is within the ≤10% target. The project owner later
closed G2 without a same-method unprofiled repeat; the conflicting timing methods and small
profile sample are retained as an explicit runtime caveat.

The function profile places `shift_plates_quasistatic` at 1.988 s cumulative over three steps
(0.663 s/step); `solve_quasistatic_contacts` itself takes 0.008 s total (about 0.003 s/step).
The coupled linear/complementarity solver is therefore not the hotspot. Other large cumulative
costs over the same sample include erosion (2.60 s/step), fault updates (1.89 s/step), hydrology
(1.38 s/step), and gap filling (1.26 s/step); these overlap with thread waits and are not
additive. The pstats report and raw output are under
`/private/tmp/rfc0012-density4-function-profile/`.

**Matched baseline cProfile rerun:** baseline median was 8.634 s/step; G2 median was 9.163 s/step,
an increase of 6.1% in this identically profiled sample. The baseline sample had 134,251 nodes and
41 plates; G2 had 129,305 nodes and 43 plates. Motion assembly was 0.568 s/step baseline and
0.663 s/step G2, while the G2 contact solver remained about 0.003 s/step. Thus this paired
cProfile sample is under the 10% limit. Its three samples per model are small, and the earlier
phase-tapped sample conflicts. On 2026-10-22, the project owner accepted this as a provisional
pass and closed G2 without a clean unprofiled repeat. Do not describe the runtime result as a
confirmed unprofiled ≤10% pass.
Baseline pstats and reports are under `/private/tmp/rfc0012-density4-baseline-function-profile/`.

The initial phase-tapped density-4 evolutions took 1,958 s baseline and 2,184 s G2. The matched
cProfile evolutions took about 1,846 s baseline and 1,944 s G2. Compact timing data are in
[`g2_density4_500myr_runtime.json`](g2_density4_500myr_runtime.json); the phase-tapped harness is
in `/private/tmp/rfc0012_density4_late_profile.py`.

## G2 closeout and caveats (2026-10-22)

The project owner closed G2 as a supported reduced-interface pilot, allowing RFC 0012 to
proceed to G3. This disposition does not claim that every original numerical acceptance bound
passed. The density-4 runtime evidence is a provisional pass: matched cProfile medians were
+6.1%, but only three steps per model were profiled and the earlier phase-tapped G2 sample was
+94.8%. No clean unprofiled repeat was run.

The four-seed 500 Myr timestep comparison passes the preregistered endpoint land-fraction
gates only. It does not establish that every material ledger or discrete topology event is
robust to timestep. The ledger traces reconcile zero subduction and merge-cap totals as real
outcomes; the merge-opportunity diagnostic finds far fewer convergent candidates and merges
in G2, but does not prove contact braking alone caused the geological differences. The higher
land fraction is therefore diagnostic, not model validation.

Follow-up issues #493 and #494 retain the early-craton/quiet-gate question and high-rate
remnant response. The 5 km glacier-depth cap also limits interpretation of erosion and ice
unloading for seed 580014012. None of these follow-ups blocks the G2 gate disposition; they
remain explicit model limitations for later work.

## Review-fix replay preregistration (four seeds, 2 Myr, 500 Myr)

**Question:** Do the PR #491 review fixes preserve the four-seed G2 endpoint behavior?

**Controls:** generate each world from its seed at node density 1.0; advance 250 steps of 2 Myr
to 500 Myr; compare current `origin/main` with the reviewed G2 branch at the exact commits
recorded in the run data. Pair baseline and G2 by seed, run every replay serially, and keep
full JSON outputs and logs under `/private/tmp`.

**Primary measure:** endpoint land fraction. **Secondary measures:** plate count, peak rate,
mean Hc/Hm, craton ledger, and median step time. The previous G2 result against the older
`2279d7e16c7d` baseline showed a positive land difference on all four seeds; that is context,
not a prediction for the current `origin/main`, which has since advanced. No numerical pass
threshold is preregistered for this review regression check. Report any seed whose land
difference changes sign and any change in high-rate outliers or ledger outcomes; do not
reinterpret these descriptive results as a new acceptance gate.

**Expected direction:** no directional prediction for the current-baseline comparison. The
efficiency update is algebraically equivalent; the pair-side averaging correction changes
drag when only one side is sampled, and pinned motion is corrected for scripted worlds. Those
behavior changes can alter long-run trajectories.

**Result:** all eight runs completed serially: baseline then candidate for each seed. Baseline
is current `origin/main` `f3a1e4603cbab4bbbebfae0ff87e95e48c678f8c`; the candidate is review-fix
commit `f2add3e1f5a2d2e03d2ab22945cb1da54a0a7b0c`. Full per-step JSON and logs remain in
`/private/tmp/rfc0012-review-replay/`; the compact data are in
[`g2_review_fix_replay.json`](g2_review_fix_replay.json).

| Seed | Baseline land | G2 land | Δ percentage points | G2 peak rate | Peak plate size | Formed volume, baseline → G2 (million km³) |
|---:|---:|---:|---:|---:|---:|---:|
| 326494039 | 13.26% | 22.16% | +8.90 | 130.1 cm/yr | 6 nodes | 384 → 876 |
| 748935664 | 9.35% | 13.69% | +4.34 | 161.5 cm/yr | 4 nodes | 191 → 359 |
| 580014012 | 21.51% | 26.92% | +5.42 | 8.8 cm/yr | 6 nodes | 549 → 920 |
| 111065103 | 25.73% | 32.98% | +7.24 | 9.4 cm/yr | 1,251 nodes | 716 → 1,250 |

Land remains higher in all four candidate runs, averaging **+6.47 percentage points**.
The two high-rate outliers occur on the six-node plate 18 and four-node plate 46, matching the
small-fragment pattern already assigned to #494. They do not justify adding a global rate cap
to G2. Subducted volume totals 2.12 million km³ on baseline and 0.055 million km³ on G2;
merge-capped volume totals 1,402 million km³ on baseline and zero on G2. These ledger shifts
remain real model outcomes, not missing trace entries. Mean per-step medians were 2.85 s
baseline and 2.74 s G2; treat this runtime difference as descriptive only, not a performance
gate. This rerun confirms the land direction and the small-plate speed caveat on current main;
it does not explain the formation or merge/subduction mechanisms.

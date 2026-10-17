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
- **G2's baseline.** The pilot replaces, for its contact class, `contact_drag_pa_s_per_m`'s
  velocity-proportional drag and the sequential per-plate shift. It keeps basal drag (already
  area-normalized), the rheology rate law and ring cascade (deformation, not load), and the
  merge/split triggers (G3). Its order-dependence and clamp numbers are compared against
  this audit's.

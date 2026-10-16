# Collision-resolution gate (issue #452)

The measurement gate for the collision-resolution exploration (#449). It answers five
questions: where overlap builds at continental contacts (G1); whether a push from thickened
crust would stop the deformable ones (G2); what strength would hold the blocked ones (G3);
what the 84 km cap costs before any merge (G4); and what makes the third-plate overflow (G5).
**No model change.** Compact numbers are in [gate.json](gate.json): one summary per run, and
`pooled[0]` (5 Myr steps) and `pooled[1]` (2 Myr steps) over the four seeds.

**Result: the decision table's strength row.** Most overlap builds at contacts that craton
blocks, and the force that would stop those contacts lies mostly within published
lithospheric strength. This holds at both step sizes and on both model commits measured.

The 5 Myr step check fails, so the decision rests on the 2 Myr runs, as the table asks.
Most of what fails is G2 and G5, not the deciding rows:
- At 2 Myr, H\* is below the cap at 88% of deformable contacts, against 66% at 5 Myr.
- At 5 Myr, a quarter of the third-plate overflow is continent over ocean. At 2 Myr it is
  1%: it comes from the step size, not from subduction failing.

## Method

**Model.** The runs use commit `de9d647`: origin/main `6579830`, with the H4 carry fix
(#450) and RFC 0009's per-node crust rules (#451), plus the probe, which changes no model
code. A first pass on `840c5ad`, before both, merged in PR #453 (this file and `gate.json`
at `8915070`); "first pass" below compares against it.

Each world is `world.generate_world(seed)` with every other argument at its default, as in
#403: density 1, climate on, 1 Gyr. Each seed ran at 5 Myr steps, tapped and untapped, and
then at 2 Myr steps, tapped. Runs went one at a time:

    backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed <S> \
        --out <runs>/seed<S>.jsonl --capture-dir <captures>
    backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed <S> --untapped \
        --out <runs>/untapped-<S>.jsonl
    backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed <S> --step 2e6 \
        --out <runs>/step2-seed<S>.jsonl
    backend/.venv/bin/python bin/debug/probe_collision_budget.py --summarize <every run> \
        --pool <the 5 Myr runs> --pool <the 2 Myr runs> --json analysis/issue452/gate.json

A 5 Myr run takes 11–16 minutes and a 2 Myr run about 25.

**Seeds.** The issue's screen asks for at least 3 merges that book craton to
`merge_capped_m3` in 1 Gyr. In the first pass, only 326494039 of the issued seeds passed:
281892019, 713688258 and 885010827 had 0, 0 and 1 such merges. Random replacements
748935664, 580014012 and 111065103 passed, and 96733506 didn't.

#450 and #451 cut those merges sharply. The four seeds' `merge_capped_m3` fell from 2,072
to 166 M km³ at 5 Myr steps. On the new code only 326494039 still passes, with 5. The owner
chose to keep the same four seeds for a before/after comparison rather than screen again.
G1–G4 measure contacts, which don't depend on cap-losing merges; G5 and the merge checks
have fewer events.

| Seed | Merges, 5 / 2 Myr | Booking `merge_capped_m3`, 5 / 2 Myr | `merge_capped_m3`, 5 / 2 Myr (M km³) | Craton start → end at 5 Myr (M km³) | First pass `merge_capped_m3` (M km³) |
|---|---|---|---|---|---|
| 326494039 | 47 / 29 | 5 / 4 | 58 / 128 | 2,763 → 2,484 | 646 |
| 748935664 | 21 / 24 | 1 / 0 | 6 / 0 | 1,639 → 1,372 | 301 |
| 580014012 | 34 / 18 | 2 / 3 | 97 / 241 | 3,268 → 3,025 | 797 |
| 111065103 | 47 / 26 | 1 / 5 | 5 / 25 | 3,232 → 3,077 | 328 |

**Captures.** `--capture-dir` saved the 2 merges on the new code that carry RFC 0001's
signature, at about 73 MB each. They are kept outside the repo, in
`~/mantle-bloom-captures/issue452/captures` on the author's machine. The first pass's 20
captures are in `840c5ad/captures` beside them.

**The probe doesn't perturb the run.** Every tapped 5 Myr run ends with exactly the craton
volume (all 16 digits) of its untapped run, on both commits: 12 seed runs in all.

### What the probe measures

`bin/debug/probe_collision_budget.py`'s docstring has the details. In short, per contact pair
(continental on both sides, or craton within 500 km) and step:

- **Overlap.** The overlap ground is each side's `contested` nodes facing the other plate.
  This is the containment test that `plates.compute_node_overlap` shares. **D**, the
  double-counted volume, is both columns' Hc over that ground. Each overlapping node is
  classed as:
  - **craton-blocked**: RFC 0008's jam test, that is, a retreat candidate holding craton,
    which `cratons.retreat_allowed` vetoes, or craton on both sides;
  - **cap-blocked**: either column within 50 m of `SUTURE_ACCRETION_MAX_HC_M`;
  - **deformable**: otherwise.

  D's change from one prepass to the next is booked node by node, by (plate, cell key), to
  the node's class.
- **Stall force.** The probe keeps each plate's shift solve, `(I + gK) ω = I ω_old + g(τ + b)`,
  with its contact drag split out per neighbour. It re-solves the pair jointly, with the
  contact drag coupling the two plates. A force *f* per metre of trench, normal to the
  contact and pushing the plates apart, is spread over both sides' edge nodes. The stall
  force is the *f* that brings the edge's mean normal closing to zero.

  This replaces the issue's plan to share each plate's torque out over its boundaries by
  hand. The joint solve already holds basal drag, slab pull and drag, ridge push and the
  plates' other contacts, as the model has them. The check: the joint re-solve without *f*
  reproduces the model's closing rate, with a median error of 0.002–0.003 cm/yr and a p90
  under 0.03 cm/yr. `mantle.clamp_rate` is left out.
- **Push.** The gravitational potential energy contrast between the contact band's
  continental columns (Hc and Hm at both sides' edge nodes, unclipped Airy columns from
  `isostatic_elevation`) and the reference continental column, integrated down to 400 km
  below sea level. **H\*** is the uniform band Hc, at the band's mean Hm, at which the push
  equals the stall force. The closing rate with the push added keeps only the normal
  component of the contact drag.
- **Budget.** The probe measures D's growth over the shift (this prepass against the end of
  the last deform pass), the D in cells that the deform pass consumed (retreat), and band
  thickening. It also books the deform pass's sinks: `no_outlet_subducted_m3`,
  `ceiling_overflow_no_outlet_m3` and `collision_subducted_m3`. Each plate's sinks are split
  over its contacts by contact nodes. The issue's "crust driven in" figure (closing × dt ×
  length × Hc) comes out at about 3× D's measured growth over the shift, so the shares
  below use the measured growth.

## Step check

Seed 326494039 failed the issue's 10-point test, so all four seeds were rerun at 2 Myr steps.
Pooled over the four seeds:

| Share | 5 Myr | 2 Myr | Moved |
|---|---|---|---|
| G1, blocked share of net growth | 56% | 73% | +17 |
| G1, the same on continental pairs | 79% | 84% | +5 |
| G1, blocked share in the 30 Myr before merges | 87% | 95% | +8 |
| G2, H\* below the cap | 66% | 88% | +22 |
| G2, push already ≥ stall | 13% | 21% | +8 |
| G3, craton-blocked stall ≤ 10¹³ N/m | 90% | 93% | +3 |
| G3, in 10¹²–10¹³ N/m | 70% | 68% | −2 |
| G5, ocean share of third-plate overflow | 26% | 1% | −25 |

- **G1 and G3 point the same way at both step sizes.** The finer step makes the blocked
  share larger.
- **G2 moves because Hm does.** The contact bands' median Hm is 93 km at 5 Myr and 49 km at
  2 Myr. Thinner roots push harder (G2 below).
- **G5's continent over ocean nearly vanishes at 2 Myr.** At 5 Myr a plate moves tens of km
  per step, and retreat doesn't keep up.
- **The lattice limits the depth figures.** Penetration depth falls from 1.15 to 1.02 cells
  (median), against a floor of about one cell.

**The tables below give 2 Myr numbers, with 5 Myr in brackets where it matters.**

## G1. Where does overlap grow?

The table gives the share of D's growth at **blocked** contacts (craton plus cap), at 2 Myr
steps:

| Seed | All pairs, net | All pairs, gross | All pairs, 30 Myr before merges | Continental pairs, net | Continental pairs, 30 Myr before merges |
|---|---|---|---|---|---|
| 326494039 | 75% | 43% | 100% | 84% | 100% |
| 748935664 | 77% | 46% | 91% | 84% | 91% |
| 580014012 | 74% | 35% | 100% | 86% | 100% |
| 111065103 | 67% | 29% | 91% | 83% | 92% |
| **Pooled** | **73%** (5 Myr 56%) | **37%** (41%) | **95%** (87%) | **84%** (79%) | **95%** (89%) |

- **Overlap accumulates where craton blocks.** "Net" sums each step's change, losses
  included; "gross" counts growth only. Deformable contacts open up more overlap step by
  step: blocked contacts make only 37% of gross growth. But retreat consumes most of the
  deformable contacts' overlap again. What stays, and what a merge eventually meets, is the
  craton-blocked overlap.
- **Almost all of the blocked share is craton.** The cap blocks about 1% of gross growth
  and nothing net.
- **Penetration is shallow:** the mean overlapping node sits 1.0 cell (128 km) behind its
  leading edge (median), and the deepest node per pair-step 1.3 cells (p90 3.2).

**Decides:** blocked contacts dominate, so the push row is out.

## G2. At deformable contacts, does the push stop the plates, or only slow them?

Closing pair-steps whose overlap (or, without overlap, whose edge) is mostly deformable, at
2 Myr steps:

| Seed | Pair-steps | Stall force, p50 (10¹² N/m) | H\* p50 (km) | H\* below the cap | Push already ≥ stall |
|---|---|---|---|---|---|
| 326494039 | 8,267 | 5.2 | 51.9 | 91% | 24% |
| 748935664 | 4,768 | 5.4 | 57.6 | 84% | 14% |
| 580014012 | 12,969 | 3.8 | 46.7 | 87% | 20% |
| 111065103 | 15,181 | 5.3 | 49.6 | 90% | 22% |
| **Pooled** | **41,185** | **4.9** (p10 0.7, p90 14) | **50.0** | **88%** (5 Myr 66%) | **21%** (13%) |

- **The review's estimate holds.** It put stalling at 0.6–3 × 10¹² N/m, reached at
  H ≈ 45–55 km. Measured, the median stall force is 4.9 × 10¹² N/m and the median H\* is
  50 km.
- **At today's Hc the push would stop only 21% of these contacts.** Elsewhere it would slow
  them, from a median 0.86 to 0.62 cm/yr. It would stop them only once the band thickens to
  H\*.
- **Hm sets the push's sign.** The model's mantle lithosphere (3,300 kg/m³) is denser than
  its asthenosphere (3,250 kg/m³), so a thick root pulls the column down.
  - At Hc = 70 km, the push is 6.8 × 10¹² N/m at Hm = 100 km, falls to zero at
    Hm ≈ 170 km, and is −9 × 10¹² N/m at Hm's 240 km cap. The zero crossing is at
    Hm ≈ 130 km for Hc = 50 km and Hm ≈ 200 km at the Hc cap.
  - The bands' Hm runs 12 / 49 / 122 km (p10 / p50 / p90) at 2 Myr, and 27 / 93 / 185 km at
    5 Myr. Their thickening depends on step size: the step check's main finding.
  - H\* is computed at the band's current Hm. Thickening that raises Hm along with Hc would
    push H\* up.
- **Sanity check (Tibet-like bands, Hc 65–75 km).** The model's own column gives a median
  push of 1.7 × 10¹² N/m at 2 Myr (p10 −2.1, p90 5.8 × 10¹²), and −1.6 × 10¹² at 5 Myr. The
  crust-only Airy figure gives 8.1 × 10¹² N/m. The model's thick roots remove most of the
  push, so the 3–10 × 10¹² N/m sanity band is reached only where Hm is near its reference.
- **Elevation clipping is minor:** 1% of bands at 2 Myr, 3% at 5 Myr.
- **Slab pull supplies most of the stall force (#315).** It is the drop in stall force with
  both plates' slab-pull torques removed: a median 81% (85% at 5 Myr). More than 10% of
  contacts would open without slab pull.

## G3. At blocked contacts, what strength would hold them?

Closing pair-steps whose overlap is mostly craton-blocked, at 2 Myr steps:

| Seed | Pair-steps | Stall force, p50 (10¹² N/m) | ≤ 10¹³ N/m | In 10¹²–10¹³ | Below 10¹² |
|---|---|---|---|---|---|
| 326494039 | 1,853 | 2.5 | 92% | 75% | 17% |
| 748935664 | 1,255 | 1.7 | 93% | 63% | 31% |
| 580014012 | 1,826 | 1.4 | 93% | 59% | 34% |
| 111065103 | 1,338 | 2.9 | 93% | 73% | 20% |
| **Pooled** | **6,272** | **2.0** (p10 0.43, p90 7.8) | **93%** (5 Myr 90%) | **68%** (70%) | **25%** (20%) |

- **A strength of 10¹³ N/m would hold 93% of craton-blocked contact-steps.** The integrated
  strength of continental lithosphere is 10¹²–10¹³ N/m. A quarter of these contacts would be
  held by less than 10¹² N/m.
- **Craton-blocked contacts need less force to stop than deformable ones** (median 2.0
  against 4.9 × 10¹² N/m). Yet nothing in the model resists them in proportion to how far
  they have pushed in. They keep closing, and D grows at a median 1.0 M km³/Myr per pair.
- Cap-blocked contacts are too few to judge on their own (92 closing pair-steps, median
  1.4 × 10¹² N/m, 91% ≤ 10¹³).

**Decides:** a yield law with physical strength values would hold these contacts, so the
strength row applies.

## G4. What does the cap cost before any merge?

All volumes are in M km³ over 1 Gyr, at 2 Myr steps.

| Seed | D's growth over the shift | Consumed by retreat | Left as overlap | Band thickening | Hc to `no_outlet` before merges | Hc lost in merges | Craton lost near contacts, excluding fusions | `merge_capped_m3` |
|---|---|---|---|---|---|---|---|---|
| 326494039 | 45,820 | 50% | 20% | 8,171 | 0 | 220 | 593 | 128 |
| 748935664 | 26,382 | 48% | 24% | 4,821 | 0 | 0 | 298 | 0 |
| 580014012 | 56,483 | 52% | 14% | 11,874 | 0 | 217 | 591 | 241 |
| 111065103 | 69,550 | 53% | 10% | 12,479 | 0 | 246 | 406 | 25 |

- D counts both columns over the overlap ground, so its volumes are about twice the crust
  involved.
- **On the new code the cap costs nothing before merges.** Nothing was booked to
  `no_outlet_subducted_m3` or `ceiling_overflow_no_outlet_m3` at either step size. In the
  first pass that was 0–3% of the merges' Hc loss. `collision_subducted_m3` is booked
  entirely in merges.
- **Retreat consumes about half** of the overlap the shift creates, and 10–24% is left.
- **Craton lost within 500 km of contacts, outside fusions,** exceeds `merge_capped_m3`
  several times over. That isn't the cap: it is rifting, erosion and subduction attrition
  near the margins (see each run's `craton_ledger_km3`).

**Decides:** the cap's losses before merges aren't comparable to the merges', so the RFC
needn't cover suture placement.

## G5. What causes the third-plate overflow?

At each merge, the probe measures the absorbed plate's Hc over third plates. It splits it by
whether the third plate's crust under it is continental or oceanic. All merges are included,
in M km³:

| Seed | 5 Myr: on continent / on ocean (ocean share) | 2 Myr: on continent / on ocean (ocean share) |
|---|---|---|
| 326494039 | 1,577 / 393 (20%) | 2,050 / 17 (1%) |
| 748935664 | 744 / 478 (39%) | 833 / 15 (2%) |
| 580014012 | 1,391 / 205 (13%) | 1,529 / 14 (1%) |
| 111065103 | 1,114 / 637 (36%) | 900 / 20 (2%) |
| **Pooled** | **4,825 / 1,713 (26%)** | **5,313 / 66 (1%)** |

- **Continent over ocean at merges comes from the 5 Myr step.** At 5 Myr steps it is a
  quarter of the overflow, as the review found on #448's captures (23%). At 2 Myr steps it
  is 1%. So the 5 Myr default leaves continent over ocean that retreat can't consume
  within the step; subduction isn't failing.
- **The overflow left is continent over a third continent.** This is a jam of more than two
  plates, the case the RFC should cover. At the 2 Myr merges that lost craton to the cap, it
  had existed for a median 25 Myr.

**Decides:** at 2 Myr steps, continent over ocean isn't a large share, so it needs no separate
physics issue. Whether the 5 Myr default step is too coarse for collisions is the owner's
call.

## Checks

At 2 Myr steps, pooled, with 5 Myr in brackets:

- **Pair symmetry.** At jammed contacts, the two plates' model contact-drag normal forces
  differ by a median 25% of the larger (29%; first pass 51%), p90 81%. The model's contact
  force isn't equal and opposite.
- **Splits.** 86 of 335 splits (26%) were of a plate in a jammed pair (47%), and 211 (63%)
  of a plate with overlap. Some jams still resolve by splitting.
- **Merges.** All 12 merges that lost craton to the cap were jammed at the merge, for a
  median 53 Myr (p10 9, p90 159).
  - At the merge the plates were still rotating apart at a median 0.13°/Myr (all merges
    0.16°/Myr), about 1.4 cm/yr at the surface. Their rotations didn't match.
  - **#402's ceiling.** D at the start of the jam episode was 29% of D at the merge, and
    craton-on-craton volume 20% (32% and 22% at 5 Myr). Merging at jam onset would remove
    at most about 70% of what those merges reconcile.

## Against the first pass (`840c5ad`, 5 Myr, same four seeds)

| | First pass | New code, 5 Myr | New code, 2 Myr |
|---|---|---|---|
| `merge_capped_m3`, total (M km³) | 2,072 | 166 | 394 |
| G1, blocked share of net growth | 67% | 56% | 73% |
| G1, continental pairs, net | 86% | 79% | 84% |
| G2, H\* below the cap | 75% | 66% | 88% |
| G3, craton-blocked stall p50 (10¹² N/m) | 1.8 | 2.8 | 2.0 |
| G3, ≤ 10¹³ N/m | 89% | 90% | 93% |
| G5, ocean share | 25% | 26% | 1% |

#450 and #451 cut merge losses by an order of magnitude. They don't change where overlap
builds or what would stop it.

## Decision

Rows judged on the 2 Myr runs, as the table asks when the step check fails.

| Outcome | Measured | Selected |
|---|---|---|
| Deformable contacts make most of the overlap growth (G1), and H\* < 84 km at most of them (G2) | Blocked contacts make 73% of net growth (84% on continental pairs, 95% in the 30 Myr before merges). | No |
| Blocked contacts make most of the overlap growth (G1), and the force to stop them is within 10¹²–10¹³ N/m (G3) | Blocked: as above, nearly all craton. Stall force at craton-blocked contacts: median 2.0 × 10¹² N/m, 93% ≤ 10¹³. | **Yes: RFC centred on strength** |
| The push only slows (H\* > 84 km at most contacts), or the force to stop blocked contacts is beyond any strength | H\* > cap at 12% of deformable contacts. Above 10¹³ N/m at 7% of craton-blocked contacts. | No |
| Continent on ocean is a large share of third-plate overflow (G5) | 1% at 2 Myr. The 26% at 5 Myr comes from the step size. | No |
| The cap's losses before merges are comparable to merge losses (G4) | None booked before merges. | No |
| The step check fails | It does (G1 +17, G2 +22, G5 −25 points). | Decided at 2 Myr, as above |

The checks support the items that every branch of the RFC keeps: one contact force per pair,
equal and opposite (today's sides differ by a median 25%); merging when rotations match
(today's merges rotate apart at 0.13°/Myr); and absorbing overlap without losing volume.

### For the RFC

- **Strength in place of the craton veto.** The contacts that keep overlap are the ones where
  `retreat_allowed` vetoes retreat, or where craton meets craton. They need modest forces to
  stop, mostly within lithospheric strength. A yield law that resists up to an integrated
  strength, and lets the contact deform beyond it, has room to work.
- **A push needs Hm.** With the model's densities, thickening Hm along with Hc turns the push
  negative above Hm ≈ 130–200 km. Band Hm depends on step size. A push term from the model's
  own column would act only where the root stays thin, or once it delaminates. The crust-only
  Airy formula overstates it about 5× at Tibet-like bands.
- **Slab pull drives the collisions.** At deformable contacts, about 80% of the stall force
  comes from slab pull elsewhere on the two plates (#315). Any resistance law will be tested
  mostly against that term.
- **Step size.** The collision numbers here move with the step. The RFC's acceptance runs
  should use 2 Myr steps, or show they don't depend on it.

## Caveats

- On the new code the four seeds have few cap-losing merges: 9 at 5 Myr and 12 at 2 Myr. G5's
  merge figures and the jam-onset ceiling rest on those few.
- The stall force holds the edge's mean closing at zero, with the force spread evenly along
  the contact. The force needed to stop every node would be larger.
- The pre-merge sinks are booked per plate and split over its contacts by contact nodes.
  The world totals are exact; the split per contact is approximate.
- The 2 Myr runs weren't paired with untapped runs. The probe's non-perturbation was shown at
  5 Myr, and it doesn't depend on the step.

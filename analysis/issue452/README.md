# Collision-resolution gate (issue #452)

The measurement gate for the collision-resolution exploration (#449). It answers five
questions: where overlap builds at continental contacts (G1); whether a push from thickened
crust would stop the deformable ones (G2); what strength would hold the blocked ones (G3);
what the 84 km cap costs before any merge (G4); and what makes the third-plate overflow (G5).
**No model change.** Compact numbers are in [gate.json](gate.json): one summary per run, and
`pooled` over the four qualifying seeds.

**Result: the decision table's strength row.** Most overlap builds at contacts that craton
blocks, and the force that would stop those contacts lies mostly within published
lithospheric strength. Continent over ocean is a large enough share of the third-plate
overflow to need its own issue.

## Method

The runs use model commit `840c5ad` (origin/main after #445 and #448; the H4 carry fix #450
hadn't landed). Each world is `world.generate_world(seed)` with every other argument at its
default, as in #403: density 1, climate on, 5 Myr steps, 1 Gyr. Runs went one at a time:

    backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed <S> \
        --out <runs>/seed<S>.jsonl --capture-dir <captures>
    backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed <S> --untapped \
        --out <runs>/untapped-<S>.jsonl
    backend/.venv/bin/python bin/debug/probe_collision_budget.py --summarize <runs>/seed*.jsonl \
        --pool <the four qualifying runs> --json analysis/issue452/gate.json

A tapped run takes 11–16 minutes and an untapped one 10–12.

**Seeds and screening.** A seed qualifies if its 1 Gyr run has at least 3 merges that book
craton to `merge_capped_m3`. Of the four seeds the issue drew, only 326494039 qualified.
Replacements came from `secrets.randbelow(10**9)`, in order:

| Seed | Source | Merges | Booking `merge_capped_m3` | Craton, start → end (M km³) | `merge_capped_m3` (M km³) | Qualifies |
|---|---|---|---|---|---|---|
| 326494039 | issue | 40 | 14 | 2,763 → 1,522 | 646 | yes |
| 281892019 | issue | 0 | 0 | 450 → 212 | 0 | no |
| 713688258 | issue | 6 | 0 | 1,528 → 1,406 | 0 | no |
| 885010827 | issue | 19 | 1 | 1,844 → 1,408 | 0.03 | no |
| 748935664 | draw 1 | 21 | 5 | 1,639 → 666 | 301 | yes |
| 96733506 | draw 2 | 8 | 0 | 99 → 107 | 0 | no |
| 580014012 | draw 3 | 30 | 6 | 3,268 → 1,742 | 797 | yes |
| 111065103 | draw 4 | 41 | 8 | 3,232 → 2,279 | 328 | yes |

Draws 5 and 6 (505783322, 36973970) weren't needed and weren't run. The headline numbers
pool the four qualifying seeds. The screened-out seeds are in `gate.json` too; they agree
on G2 and G3, and on G1's 30 Myr before each merge. Every row below covers one seed.

**Captures.** `--capture-dir` saved 20 merges with RFC 0001's signature, at about 73 MB each:
6, 3, 6 and 5 for the four qualifying seeds (6 is the cap per seed). They are kept outside the repo, in
`~/mantle-bloom-captures/issue452/captures` on the author's machine, so
`measure_merge_placement.py` can replay them.

**The probe doesn't perturb the run.** Every tapped run ends with exactly the craton volume
(all 16 digits) of its untapped run.

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
  reproduces the model's closing rate, with a median error of 0.000–0.001 cm/yr and a p90
  of 0.02 cm/yr. `mantle.clamp_rate` is left out.
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
  length × Hc) comes out at 2.8–3.0× D's measured growth over the shift, so the shares
  below use the measured growth.

## G1. Where does overlap grow?

The table gives the share of D's growth at **blocked** contacts (craton plus cap):

| Seed | All pairs, net | All pairs, gross | All pairs, 30 Myr before merges | Continental pairs, net | Continental pairs, 30 Myr before merges |
|---|---|---|---|---|---|
| 326494039 | 75% | 66% | 85% | 87% | 88% |
| 748935664 | 54% | 49% | 80% | 80% | 80% |
| 580014012 | 64% | 50% | 81% | 85% | 83% |
| 111065103 | 68% | 53% | 86% | 87% | 87% |
| **Pooled** | **67%** | **55%** | **84%** | **86%** | **85%** |

- "Net" sums each step's change, losses included. "Gross" counts growth only.
- **Almost all of the blocked share is craton.** The cap blocks 0.3% (net) to 2% (gross) of
  the growth. Columns at the cap are rare where overlap is growing.
- **Continental pairs are where the blocking concentrates.** The all-pairs figures also
  count ocean–continent pairs with craton in reach, which are mostly deformable.
- **Penetration is shallow.** The mean overlapping node sits 1.4 cells (179 km) behind its
  plate's leading edge (median; p90 3.3 cells). The deepest node per pair-step is at 3.4
  cells median, 9.6 cells at p90. The cell lattice quantizes these: a single row of overlap
  reads as about one cell.

**Decides:** blocked contacts dominate, so the push row is out.

## G2. At deformable contacts, does the push stop the plates, or only slow them?

Closing pair-steps whose overlap (or, without overlap, whose edge) is mostly deformable:

| Seed | Pair-steps | Stall force, p10 / p50 / p90 (10¹² N/m) | H\* p50 (km) | H\* below the cap | Push already ≥ stall | Closing now, p50 (cm/yr) | Closing with the push, p50 | Slab-pull share of stall, p50 |
|---|---|---|---|---|---|---|---|---|
| 326494039 | 2,349 | 1.0 / 6.6 / 19.1 | 52.6 | 75% | 26% | 0.90 | 0.55 | 86% |
| 748935664 | 1,395 | 1.3 / 5.3 / 18.7 | 56.1 | 69% | 18% | 0.92 | 0.77 | 80% |
| 580014012 | 4,142 | 0.4 / 3.1 / 14.2 | 52.8 | 78% | 17% | 0.67 | 0.72 | 85% |
| 111065103 | 4,675 | 0.6 / 5.0 / 17.2 | 50.6 | 74% | 24% | 0.90 | 0.59 | 87% |
| **Pooled** | **12,561** | **0.5 / 4.6 / 17.2** | **52.3** | **75%** | **21%** | **0.82** | **0.65** | **85%** |

- **The review's estimate holds.** It put stalling at 0.6–3 × 10¹² N/m, reached at
  H ≈ 45–55 km. Measured, the median stall force is 4.6 × 10¹² N/m and the median H\* is
  52 km. H\* is below the cap at three quarters of deformable contacts.
- **At today's Hc the push would stop only 22% of these contacts.** Elsewhere the push
  slows them, from a median 0.82 to 0.65 cm/yr. It would stop them only once the band
  thickens to H\*.
- **Hm sets the push's sign.** The model's mantle lithosphere (3,300 kg/m³) is denser than
  its asthenosphere (3,250 kg/m³), so a thick root pulls the column down.
  - At Hc = 70 km, the push is 6.8 × 10¹² N/m at Hm = 100 km, falls to zero at
    Hm ≈ 170 km, and is −9 × 10¹² N/m at Hm's 240 km cap. The zero crossing is at
    Hm ≈ 130 km for Hc = 50 km and Hm ≈ 200 km at the Hc cap.
  - The contact bands' Hm runs 25 / 91 / 170 km at p10 / p50 / p90.
  - H\* is computed at the band's current Hm. Thickening that raises Hm along with Hc would
    push H\* up.
- **Sanity check (Tibet-like bands, Hc 65–75 km, 930 pair-steps).** The model's own column
  gives a median push of 1.3 × 10¹² N/m (p10 −6.6 × 10¹², p90 1.1 × 10¹³). The crust-only
  Airy figure gives 8.4 × 10¹² N/m. The 3–10 × 10¹² N/m sanity band is met only where Hm is
  near the reference: seed 326494039, where Hm is the thinnest, gives 4.5 × 10¹². So
  crust-only Airy is the upper bound the review expected, and the model's thick roots
  remove most of the push.
- **Elevation clipping is minor.** 5% of the bands would clip at `MAX_ELEVATION_M`.
- **Slab pull supplies most of the stall force (#315).** It is the drop in stall force with
  both plates' slab-pull torques removed. At the median deformable contact, 85% of the
  force needed to stop it is slab pull from oceanic slabs elsewhere on the two plates. At
  more than 10% of contacts, the pair would open without slab pull (share > 1).

## G3. At blocked contacts, what strength would hold them?

Closing pair-steps whose overlap is mostly craton-blocked:

| Seed | Pair-steps | Stall force, p10 / p50 / p90 (10¹² N/m) | ≤ 10¹³ N/m | In 10¹²–10¹³ | Below 10¹² | D growth, p50 (M km³/Myr) |
|---|---|---|---|---|---|---|
| 326494039 | 1,090 | 0.7 / 3.0 / 13.9 | 83% | 68% | 15% | 1.41 |
| 748935664 | 399 | 0.6 / 1.6 / 6.2 | 94% | 66% | 28% | 2.22 |
| 580014012 | 1,260 | 0.1 / 1.0 / 5.8 | 93% | 42% | 51% | 0.55 |
| 111065103 | 1,477 | 0.2 / 2.1 / 11.1 | 88% | 56% | 33% | 0.78 |
| **Pooled** | **4,226** | **0.25 / 1.8 / 10.8** | **89%** | **56%** | **33%** | **0.91** |

- **A strength of 10¹³ N/m would hold 89% of craton-blocked contact-steps.** The integrated
  strength of continental lithosphere is 10¹²–10¹³ N/m. A third of these contacts would be
  held by less than 10¹² N/m.
- **Craton-blocked contacts need less force to stop than deformable ones** (median 1.8
  against 4.6 × 10¹² N/m). Yet nothing in the model resists them in proportion to how far
  they have pushed in. They keep closing, and D grows at a median 0.9 M km³/Myr per pair.
- Cap-blocked contacts are too few to judge on their own (39 closing pair-steps pooled,
  median 1.4 × 10¹² N/m, 95% ≤ 10¹³).

**Decides:** a yield law with physical strength values would hold these contacts, so the
strength row applies.

## G4. What does the cap cost before any merge?

Of D's growth over the shift, the deform pass's retreat consumes about half. Most of the
rest is left as overlap. All volumes are in M km³ over 1 Gyr.

| Seed | D's growth over the shift | Consumed by retreat | Left as overlap | Band thickening | Hc to `no_outlet` before merges | Hc lost in merges | Before / merges | Craton lost near contacts, excluding fusions | `merge_capped_m3` |
|---|---|---|---|---|---|---|---|---|---|
| 326494039 | 53,004 | 48% | 46% | 8,933 | 29 | 1,995 | 1.5% | 766 | 646 |
| 748935664 | 39,204 | 52% | 37% | 4,721 | 0 | 687 | 0% | 707 | 301 |
| 580014012 | 71,410 | 54% | 29% | 13,534 | 40 | 1,333 | 3.0% | 1,077 | 797 |
| 111065103 | 90,882 | 52% | 31% | 15,913 | 59 | 2,158 | 2.7% | 962 | 328 |

- D counts both columns over the overlap ground, so its volumes are about twice the crust
  involved.
- **The cap costs little before merges.** Crust booked to `no_outlet_subducted_m3` during
  deform is 0–3% of the Hc the merges lose. Nothing was booked to
  `ceiling_overflow_no_outlet_m3`. Almost all of `collision_subducted_m3` is booked in
  merges: 97–100%.
- **Craton lost within 500 km of contacts, outside fusions,** is 1.2–2.9× the merges'
  `merge_capped_m3`. That isn't the cap: it is rifting, erosion and subduction attrition
  near the margins (see each seed's `craton_ledger_km3` in `gate.json`).

**Decides:** the cap's losses before merges aren't comparable to the merges', so the RFC
needn't cover suture placement.

## G5. What causes the third-plate overflow?

At each merge, the probe measures the absorbed plate's Hc over third plates. It splits it by
whether the third plate's crust under it is continental or oceanic. All merges are included:

| Seed | On third-plate continent (M km³) | On third-plate ocean (M km³) | Ocean share |
|---|---|---|---|
| 326494039 | 4,744 | 691 | 13% |
| 748935664 | 1,218 | 1,588 | 57% |
| 580014012 | 3,071 | 1,424 | 32% |
| 111065103 | 4,650 | 766 | 14% |
| **Pooled** | **13,682** | **4,469** | **25%** |

- **A quarter of the overflow is continent over ocean.** That matches the review's 23% on
  #448's captures. At the 33 merges that lost craton to the cap, the area-weighted overlap
  over ocean had existed for a median 20 Myr by the merge (p90 48 Myr). That is four or
  more steps in which subduction left it in place. The rest, continent over a third
  continent, had existed for a median 41 Myr.
- **The share depends on step size.** At 2 Myr steps, seed 326494039's ocean share drops
  from 13% to 3%. The follow-up issue should start by checking how it depends on step size.

**Decides:** a separate issue for why subduction didn't remove the continent over ocean.

## Checks

- **Pair symmetry.** At jammed contacts, the two plates' model contact-drag normal forces
  differ by a median 51% of the larger (p90 94%). `shift_plate` solves each plate with its
  neighbour's ω held fixed, and each side's coefficient comes from its own nodes. The model's
  contact force isn't equal and opposite.
- **Splits.** 200 of the qualifying seeds' 419 splits (48%) were of a plate in a jammed pair,
  and 255 (61%) of a plate with overlap. Splitting is a common way out of a jam.
- **Merges.** All 33 merges that lost craton to the cap were jammed at the merge. They had
  been jammed for a median 45 Myr (p10 20, p90 131).
  - At the merge the plates were still rotating apart at a median 0.16°/Myr (p90
    0.32°/Myr), about 1.7 cm/yr at the surface. The edges were closing at 0.17 cm/yr. Their
    rotations didn't match.
  - **#402's ceiling.** Over those merges, D at the start of the jam episode was 50% of D
    at the merge, and craton-on-craton volume 48%. Merging at jam onset would at most
    halve what the merge has to reconcile.
- **Step size.** Seed 326494039 at 2 Myr steps against 5 Myr:

  | Share | 5 Myr | 2 Myr | Moved |
  |---|---|---|---|
  | G1, blocked share of net growth | 75% | 83% | +8 |
  | G1, blocked share in the 30 Myr before merges | 85% | 93% | +7 |
  | G1, continental pairs, net | 87% | 89% | +2 |
  | G2, H\* below the cap | 75% | 74% | −1 |
  | G2, push already ≥ stall | 26% | 26% | 0 |
  | G3, craton-blocked stall ≤ 10¹³ N/m | 83% | 81% | −2 |
  | G3, in 10¹²–10¹³ N/m | 68% | 69% | +1 |

  Every share moves by less than 10 points, so the check passes. The blocked share rises at
  the finer step. Penetration depth drops from 1.41 to 1.06 cells median: the lattice
  floor.

## Decision

| Outcome | Measured | Selected |
|---|---|---|
| Deformable contacts make most of the overlap growth (G1), and H\* < 84 km at most of them (G2) | Blocked contacts make 67% of net growth (86% on continental pairs, 84% in the 30 Myr before merges). H\* < cap at 75% of deformable contacts. | No |
| Blocked contacts make most of the overlap growth (G1), and the force to stop them is within 10¹²–10¹³ N/m (G3) | Blocked: 67–86% as above, nearly all craton. Stall force at craton-blocked contacts: median 1.8 × 10¹² N/m, 89% ≤ 10¹³. | **Yes: RFC centred on strength** |
| The push only slows (H\* > 84 km at most contacts), or the force to stop blocked contacts is beyond any strength | H\* > cap at 25% of deformable contacts. Above 10¹³ N/m at 11% of craton-blocked contacts. | No |
| Continent on ocean is a large share of third-plate overflow (G5) | 25% pooled (13–57% per seed), standing for a median 20 Myr. | **Yes: separate issue** |
| The cap's losses before merges are comparable to merge losses (G4) | 0–3% of the merges' Hc loss. | No |
| The step check fails | Every share moves under 10 points. | No |

The checks support the items that every branch of the RFC keeps: one contact force per pair,
equal and opposite (today's sides differ by 51%); merging when rotations match (today's
merges rotate apart at 0.16°/Myr); and absorbing overlap without losing volume.

### For the RFC

- **Strength in place of the craton veto.** The contacts that build overlap are the ones
  where `retreat_allowed` vetoes retreat, or where craton meets craton. They need modest
  forces to stop, mostly within lithospheric strength. A yield law that resists up to an
  integrated strength, and lets the contact deform beyond it, has room to work.
- **A push needs Hm.** With the model's densities, thickening Hm along with Hc turns the
  push negative above Hm ≈ 130–200 km. A push term from the model's own column would act
  only where the root stays thin, or once it delaminates. The crust-only Airy formula
  overstates it about 6× at Tibet-like bands.
- **Slab pull drives the collisions.** At deformable contacts, 85% of the stall force comes
  from slab pull elsewhere on the two plates (#315). Any resistance law will be tested
  mostly against that term.
- **Splits.** About half of all splits come from jammed pairs. A contact that resolves by
  strength would change the split rate as well as the merges.

## Caveats

- The runs predate the H4 carry fix (#450). That fix changes merge placement, not the
  contact physics this gate measures.
- The stall force holds the edge's mean closing at zero, with the force spread evenly along
  the contact. The force needed to stop every node would be larger.
- The pre-merge sinks are booked per plate and split over its contacts by contact nodes.
  The world totals are exact; the split per contact is approximate.
- The step check uses one seed.

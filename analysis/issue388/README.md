# Issue #388: craton boundary loss by the cut that made it (RFC 0003 phase 1b)

The second gate in [RFC 0003](../../rfc/0003-rifting-around-cratons.md). Phase 1 (#378)
found that divergent thinning and stretch donation take 83–90% of craton rifting, but only
12–18% of cuts cross craton. This phase asks how much of that boundary loss happens on
boundaries a craton-crossing cut created, since only that share (plus failed rifts) is what
steering the cut (solution 1) could remove.

## Method

`bin/debug/attribute_rift_boundary_loss.py --seed <S> --years 1e9`, seeds 7, 11 and 2024,
default density, 5 Myr steps (as `compare_cratons.py`), run one at a time. Baseline:
**origin/main** at 9251382, without #385 (PR #387). Numbers are in
[boundary_attribution.json](boundary_attribution.json).

- At each successful split, the probe marks the nodes of both daughters within the cut's
  band (`FAILED_RIFT_BAND_MULT` = 2.5 spacings, the band phase 1 counted with) as "through"
  if that band held craton before the split, else "clear". A later cut overwrites the mark.
- The marks are probe-only surface fields, registered in the probe's process only, and remap
  as intensive fractions. A cell blended from marked and unmarked cells keeps the marked
  share.
- Every later craton loss in `deform_columns` (divergent) and `_open_rift` (stretch) is
  booked per node by mark: "through", "clear" or "never" (never in a cut band).
- F = (failed-rift loss + divergent and stretch loss marked "through") / `rifted_m3`.
- Sensitivity: **wide** doubles the band (5 spacings); **ever** never lets a clear cut
  overwrite a through mark.

The probe changes nothing the model reads. `rifted_m3`, its path split and the split counts
match #378's runs exactly, and the per-node booking sums to the ledger's `rifted_divergent_m3`
and `rifted_stretch_m3` to within rounding (ratio 1.000000).

## Results (0–1 Gyr)

| Seed | Rifted, % of start craton | Failed | Splits through craton | Boundary loss marked through (band / wide) | **F** (band) | F (wide) | F, main + #385 (est.) |
|---|---|---|---|---|---|---|---|
| 7 | 4.0% | 17% | 5 of 33 | 25% / 34% | **37%** | 45% | 28% |
| 11 | 4.3% | 14% | 6 of 37 | 22% / 27% | **33%** | 37% | 25% |
| 2024 | 10.7% | 10% | 4 of 32 | 9% / 16% | **18%** | 24% | 11% |
| Mean | 6.3% | | | | **29%** | 35% | 21% |

"Boundary loss" is divergent plus stretch. The mean is unweighted; weighted by `rifted_m3` it
is 29% (band) and 35% (wide). The main + #385 column assumes #385 removes 75% of failed-rift
loss and leaves the rest unchanged: F = (0.25 failed + through) / (rifted − 0.75 failed).

By path, the share marked "through" is 10% / 17% / 7% of divergent loss and 34% / 25% / 10% of
stretch loss (seeds 7 / 11 / 2024).

## Findings

- **Most boundary loss is on nodes no cut ever marked:** 75–91% of divergent and stretch
  loss. Doubling the band only moves 5–10 points of it to "through", so this isn't loss
  leaking past a narrow band. It is boundaries that never came from a split in the run:
  the plate boundaries the world starts with, or boundaries that merges, defragmentation or
  boundary motion carried onto a craton.
- **Cuts clear of craton never cost craton later** (0% marked "clear" on every seed). Their
  band held no craton, and nothing put craton there within 1 Gyr. **ever** equals **band**
  on every seed: no through-marked node was cut again cleanly.
- **Seed 2024 is the worst rifter and the least steerable.** It loses 10.7% per Gyr, and only
  9% of its boundary loss sits on through-marked nodes, so F = 18%. Its divergent path (41% of
  `rifted_m3`) is almost all on unmarked nodes.
- **F is an upper bound on what steering removes.** It assumes every craton-crossing cut can
  be steered clear and that a steered boundary loses nothing. The RFC's admissibility tests
  will reject some candidates.

## Branch selected

On the stated baseline (origin/main), F is **29%** on the three-seed mean and 33–37% on seeds 7
and 11. That is the RFC's **25% ≤ F < 50%** branch: steering can't reach the 2%/Gyr target
alone; build solution 1 only together with solution 3, and size 3 to the remainder.

Two cautions keep this close to the line:
- Seed 2024 alone is below 25% (18%), and it is the seed furthest from target.
- With #385 merged, failed-rift loss shrinks, and F falls to about 21% on the mean (11–28% per
  seed), inside the "rift cuts aren't the path" case.

Even at F = 29%, steering removes at most 1.8%/Gyr of the 6.3%/Gyr mean, against the 4.3%/Gyr
the target needs. The larger lever is the 75–91% of boundary loss no cut accounts for, which
points at solution 3 or boundary migration.

## Not measured

- Which boundaries carry the "never" loss: initial plate boundaries, merge sutures that
  later turn divergent, or migration. A follow-up probe could mark initial boundaries at
  year 0 to split this.
- Seed 456862899 to 3 Gyr (the RFC's acceptance run, not this gate).

# Issue #388: craton boundary loss by the cut that made it (RFC 0003 phase 1b)

The second gate in [RFC 0003](../../rfc/0003-rifting-around-cratons.md). Phase 1 (#378)
found that divergent thinning and stretch donation take 83–90% of craton rifting, but only
12–18% of cuts cross craton. This phase asks how much of that boundary loss happens on
boundaries a craton-crossing cut created, since only that share (plus failed rifts) is what
steering the cut (solution 1) could remove.

## Method

`bin/debug/attribute_rift_boundary_loss.py --seed <S> --years 1e9`, seeds 7, 11 and 2024,
default density, 5 Myr steps (as `compare_cratons.py`), run one at a time, on two baselines:

- **main + #385** (origin/main at 649e5ad, failed rifts resist cratons). This is the model
  solution 1 would be built on, so it decides the branch.
- **main without #385** (origin/main at 9251382), the model phase 1 measured. #385 merged
  during this work.

Numbers are in [boundary_attribution.json](boundary_attribution.json), keyed by baseline.

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

The probe changes nothing the model reads. On each baseline, `rifted_m3`, its path split and
the split counts match that baseline's own replays (#378 and #385) exactly. The per-node
booking sums to the ledger's `rifted_divergent_m3` and `rifted_stretch_m3` to within rounding
(ratio 1.000000).

## Results (0–1 Gyr)

**main + #385** (selects the branch):

| Seed | Rifted, % of start craton | Failed | Splits through craton | Boundary loss marked through (band / wide) | **F** | F (wide) |
|---|---|---|---|---|---|---|
| 7 | 3.4% | 3% | 5 of 33 | 24% / 32% | **26%** | 34% |
| 11 | 5.4% | 1% | 7 of 28 | 14% / 19% | **15%** | 20% |
| 2024 | 8.6% | 5% | 5 of 32 | 6% / 10% | **11%** | 15% |
| Mean | 5.8% | | | | **17%** | 23% |

**main without #385**:

| Seed | Rifted, % of start craton | Failed | Splits through craton | Boundary loss marked through (band / wide) | **F** | F (wide) |
|---|---|---|---|---|---|---|
| 7 | 4.0% | 17% | 5 of 33 | 25% / 34% | **37%** | 45% |
| 11 | 4.3% | 14% | 6 of 37 | 22% / 27% | **33%** | 37% |
| 2024 | 10.7% | 10% | 4 of 32 | 9% / 16% | **18%** | 24% |
| Mean | 6.3% | | | | **29%** | 35% |

"Boundary loss" is divergent plus stretch, and "Failed" is the failed-rift share of
`rifted_m3`. Means are unweighted; weighted by `rifted_m3` they are 17% (main + #385) and 29%
(without). On main + #385, the share marked "through" is 9% / 11% / 5% of divergent loss and
33% / 16% / 6% of stretch loss (seeds 7 / 11 / 2024).

## Findings

- **Most boundary loss is on nodes no cut ever marked:** 76–94% of divergent and stretch
  loss on main + #385 (75–91% without). Doubling the band moves only 4–9 points of it to
  "through", so this isn't loss leaking past a narrow band. It is boundaries that never
  came from a split in the run: the plate boundaries the world starts with, or boundaries
  that merges, defragmentation or boundary motion carried onto a craton.
- **#385 took most of F's easy part.** Failed rifts were 10–17% of `rifted_m3`, all of it
  counted in F. With #385 they are 1–5%, and F falls by 12 points on the mean.
- **Cuts clear of craton never cost craton later** (0% marked "clear" on every seed and
  baseline). Their band held no craton, and nothing put craton there within 1 Gyr. **ever**
  equals **band** throughout: no through-marked node was cut again cleanly.
- **The worst rifter is the least steerable.** Seed 2024 loses 8.6%/Gyr, and only 6% of its
  boundary loss sits on through-marked nodes (F = 11%).
- **F is an upper bound on what steering removes.** It assumes every craton-crossing cut can
  be steered clear and that a steered boundary loses nothing. The RFC's admissibility tests
  will reject some candidates.

## Branch selected

On main + #385, F is **17%** on the three-seed mean (11–26% per seed). Even the wide band
gives 23%. That is the RFC's **F < 25%** case: "rift cuts aren't the path". Don't build
solution 1; revise the RFC toward boundary migration or solution 3.

Steering would remove at most 1.0%/Gyr of the 5.8%/Gyr mean, against the 3.8%/Gyr the
2%/Gyr target needs.

Without #385, F was 29%, the "build 1 only together with 3" branch. That baseline no longer
describes main, and it was 29% only because failed rifts, which #385 has since resisted,
counted fully toward F.

## Not measured

- Which boundaries carry the "never" loss: initial plate boundaries, merge sutures that
  later turn divergent, or migration. A probe that also marks the year-0 boundaries and each
  merge suture would split it, and would show whether a boundary-migration fix has a
  target.
- Seed 456862899 to 3 Gyr (the RFC's acceptance run, not this gate).

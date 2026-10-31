# Issue #388: craton boundary loss by the cut that made it (RFC 0003 phase 1b)

The second gate in [RFC 0003](../../rfc/abandoned/0003-rifting-around-cratons.md). Phase 1 (#378)
found that divergent thinning and stretch donation take most of craton rifting, but only a
minority of cuts cross craton. This phase asks how much of that boundary loss happens on
boundaries a craton-crossing cut created, since only that share (plus failed rifts) is what
steering the cut (solution 1) could remove.

## Method

`bin/debug/attribute_rift_boundary_loss.py --seed <S> --years 1e9`, seeds 7, 11 and 2024,
default density, 5 Myr steps (as `compare_cratons.py`), run one at a time, on three baselines:

- **current main** (origin/main at 8878cd2). It includes #385 (failed rifts resist cratons)
  and #358 (the split cooldown counted in years). **This one selects the branch.**
- main + #385, before #358 (649e5ad);
- main before #385 (9251382), the model phase 1 measured.

#358 matters most. At 5 Myr steps the old cooldown was 20 steps, 100 Myr; it is now 2 Myr,
the same as at the default 100 kyr step. So these runs split about three times as often
(104–115 splits per Gyr, against 28–33), and the earlier baselines' split pace was an
artifact of long steps. Both #358 and #385 merged while this work was under way. Numbers
are in [boundary_attribution.json](boundary_attribution.json), keyed by baseline.

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
the split counts match that baseline's own replays (#378, #385, and RFC 0008's ledger run
on current main) exactly. The per-node booking sums to the ledger's `rifted_divergent_m3`
and `rifted_stretch_m3` to within rounding (ratio 1.000000).

## Results (0–1 Gyr)

**Current main** (selects the branch):

| Seed | Rifted, % of start craton | Failed | Splits through craton | Boundary loss marked through (band / wide) | **F** | F (wide) |
|---|---|---|---|---|---|---|
| 7 | 13.2% | 7% | 45 of 115 | 34% / 46% | **39%** | 50% |
| 11 | 9.0% | 7% | 26 of 109 | 38% / 55% | **42%** | 58% |
| 2024 | 15.6% | 9% | 31 of 104 | 41% / 60% | **47%** | 64% |
| Mean | 12.6% | | | | **42%** | 57% |

**Earlier baselines** (F, band / wide):

| Seed | main + #385, before #358 | main before #385 |
|---|---|---|
| 7 | 26% / 34% | 37% / 45% |
| 11 | 15% / 20% | 33% / 37% |
| 2024 | 11% / 15% | 18% / 24% |
| Mean | 17% / 23% | 29% / 35% |

"Boundary loss" is divergent plus stretch, and "Failed" is the failed-rift share of
`rifted_m3`. Means are unweighted; weighted by `rifted_m3`, current main's F is 42% (band)
and 55% (wide). On current main, the share marked "through" is 32% / 31% / 34% of divergent
loss and 35% / 41% / 44% of stretch loss (seeds 7 / 11 / 2024).

## Findings

- **With the real split pace, cuts matter much more.** Craton rifting runs at 9–16% of the
  starting craton per Gyr, against 3–9% before #358. 24–39% of splits cross craton, against
  15–25%. 34–41% of boundary loss is on through-marked nodes, against 6–24%.
- **Most boundary loss is still on nodes no cut marked** (58–65%), but now only by a
  margin. Doubling the band moves 12–19 points of it to "through", more than before, so on
  this baseline part of that "never" share is loss that spread past the narrow band.
- **Cuts clear of craton almost never cost craton later** (0.2–0.3% of boundary loss), and
  **ever** equals **band** on every seed.
- **F is an upper bound on what steering removes.** It assumes every craton-crossing cut can
  be steered clear and that a steered boundary loses nothing. The RFC's admissibility tests
  will reject some candidates.

## Branch selected

On current main, F is **42%** on the three-seed mean (39–47% per seed). That is the RFC's
**25% ≤ F < 50%** branch: steering can't reach the 2%/Gyr target alone; build solution 1
only together with solution 3, and size 3 to the remainder. The wide band puts two seeds
above 50%, so the band width is the main uncertainty in that call. It doesn't move the
branch on the mean of the narrow band the RFC defines.

Steering would remove at most 5.4%/Gyr of the 12.6%/Gyr mean, against the 10.6%/Gyr the
target needs.

This supersedes the earlier result on main + #385 before #358 (F = 17%, the "rift cuts
aren't the path" case). That baseline's 5 Myr-step runs split at a third of the real pace.

## Not measured

- Which boundaries carry the "never" loss: initial plate boundaries, merge sutures that
  later turn divergent, or migration. A probe that marks each plate's boundaries by origin
  would split it.
- Seed 456862899 to 3 Gyr (the RFC's acceptance run, not this gate).

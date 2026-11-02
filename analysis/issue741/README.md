# Issue 741: climate noise over its 40 Myr period, and arm (m3)'s ice drops

Evidence for [#741](https://github.com/adubey/mantle-bloom/issues/741), which settles the banding
item of [#575](https://github.com/adubey/mantle-bloom/issues/575). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544).

**Principle:** #3, fast processes in aggregate over long steps. The noise stands for geography,
so swapping the redraw for a drift shouldn't move the climate on average.

**Invariant:** humidity and precipitation distributions by latitude band, averaged over a drift
period and over draws, match the per-step redraw's.

## Preregistration

Written before any run in this report. The data already seen are #575's runs, including arm
(m3)'s 14 Myr series, whose per-step ledger motivated Part B.

### Part A: climate over 48 Myr

**The baseline isn't `origin/main`.** #736 merged the drifting noise, so `origin/main` is now the
treatment. The only backend change between `9d6fb5e` (before #736) and `bf58116` (`origin/main`
now) is `climate.py`, so `9d6fb5e` is the old `origin/main` exactly.

| Arm | Tree | Draws |
|---|---|---|
| Drift (D) | `bf58116` plus this branch's measurement script | offsets 0, 1, 2, 3 |
| Redraw (R) | `9d6fb5e` for draw 0; `bf58116` with `--noise-redraw` for the others | 0, 1, 2 |

**Replays,** one at a time: seed 997271774 from 87.4 Myr and seed 649658380 from 6.4 Myr, 48 Myr
at 0.25 Myr steps, a checkpoint each Myr. On seed 997271774 the run passes the keyframe at 120 Myr.

**Extension, fixed in advance.** If any primary metric's 95% half-width on seed 997271774 exceeds
2%, add D offsets 4–7 and R offset 3 there. This reads the width only, never the gap, so it can't
select for a result.

**Statistic.** Each run gives one number per metric: its mean over checkpoints after the opening.
Checkpoints within a run aren't independent, so the draw is the unit. The gap is
`(mean D − mean R) / mean R`, with a 95% Welch interval on the difference.

**Primary metrics:** global land and ocean precipitation and humidity, and the land precipitation
band at 15–30°. The rest are descriptive: every band's land and ocean means and p10, p50 and p90,
with about 40 metrics in each save, so about two 95% intervals should exclude 0 by chance.

**Decision rule,** per metric, for bands with at least 5% of land:

| Verdict | Condition |
|---|---|
| Unchanged | the interval lies inside ±2% |
| Shifted | the interval excludes 0 on both saves, with the same sign |
| Unresolved | anything else |

**Expected:** no shift. The noise has one distribution at every instant on both arms, so a real
shift would need a nonlinear path through stored state (vegetation, lakes, rivers), which the
redraw averages and a persistent pattern doesn't.

**The period itself.** Each D draw's gap to the mean of R, in windows 0–16, 16–32 and 32–48 Myr.
If the offset belongs to the texture's overlap with this land, it should weaken or change sign as
the keyframes hand over. If most draws keep a same-sign gap above 2% in all three windows, that's
persistence through stored state, and goes to the owner.

### Part B: what happens in the ice on the step of a large drop

**Tree.** G0 prototype `ea9a187` with #575's `climate.py`, arm (m3), seed 997271774, draw 1, 0.25
Myr steps, to 101.15 Myr. It must reproduce the earlier series; its ice volume at 97.4 Myr was
10.53 M km³.

**Method.** On each named step, re-march that step in-process from its real opening ice, with:

| Cell | Climate | Surface | Opening ice |
|---|---|---|---|
| `new_new` | this step's | this step's | real (must equal the real march) |
| `old_new` | last step's | this step's | real |
| `new_old` | this step's | last step's | real |
| `old_old` | last step's | last step's | real |
| `ice_plus`, `ice_minus` | this step's | this step's | real ±0.1% |

Climate is precipitation, temperature and amplitude at nodes. Surface is elevation, deflection,
mobile cover, and the ocean masks. `world` itself, for hypsometry and trapped water, stays this
step's in every cell.

**Steps:** the drop at 97.4, the step after at 97.65, the swing pair 100.4 and 100.65, and controls
96.4 and 99.15.

**Primary measures:** the cell's closing volume against `new_new`, and the march's within-step
scatter (standard deviation of substep volumes in the step's second half).

**Decision rules:**

1. Controls give the cells' noise: the largest |`old_old` − `new_new`| across control steps.
2. A group of inputs causes a drop if restoring last step's values alone removes at least half of
   it, by more than that noise.
3. If neither group does but `old_old` still drops, nothing changed in the inputs: the march is
   leaving a state of its own accord.
4. If `ice_plus` and `ice_minus` differ by as much as the drop, closing volume is sensitive to its
   start, and one step's change can't be attributed at all.
5. If a step's change is within its march's scatter, that step-end value is a sample of the
   scatter, not a change in the ice.

A climate result leads to one more run, holding the noise from 97.15 Myr, because the climate at
nodes also changes through the surface. Where the ice went, by body and region, comes from the
same captures.

### Part B amendment, after the first factorial

Written after the first factorial run and before any run below. It is post hoc: the design
changed because of a result.

Rule 4 fired on the drop step: ±0.1% opening ice gave 13.14 and 9.70 M km³, a gap larger than the
drop (13.03 to 10.53). One march can't be attributed, so the question becomes a distribution. The
capture also showed the moving ice on equatorial plateaus near 5 km, with nodes going between 0
and about 2 km in a step.

| Run | What | Reads |
|---|---|---|
| E1 | 96.4, 97.4 and 99.15: 12 uniform opening-ice scalings (±0.01–0.5%), this step's climate; on 97.4 also last step's climate | Does the drop step differ from ordinary steps, and do last step's inputs change the odds? |
| E2 | 97.4 with the noise held from 97.15 Myr, same 12 scalings | Is the noise's own change needed? |
| E3 | Scale the 97.4 opening by 1.001 and continue to 101.15 Myr | Does an avoided drop stay avoided, or does the low state arrive anyway? |
| E4 | Repeat the factorial on 97.4 and 100.4 | The node captures a filename bug lost; volumes must match the first run |

**Reading:** the drop counts as forced if at least 10 of 12 members drop (closing below 11.5 M km³)
under this step's inputs and at most 2 do at control steps. It counts as a sensitive switch if the
fraction lies between. E1 under last step's climate, and E2, are read as a change in that fraction;
with 12 members, a gap under 5 members isn't called a difference.

### Part B second amendment, after the ensembles

Written after E1–E4 and before any run below, so it is post hoc too.

E1–E3 found a change in the climate that isn't the noise, then E4's march path showed its form: on
the drop step one substep takes the volume from 13.1 to 5.9 M km³, and it rebounds to about 9.8.

| Run | What | Reads |
|---|---|---|
| E5 | Draw 1, 97.4: log every substep, capture the solve's node fields on any that moves volume over 15% | What the collapsing substep does, and where its ice goes |
| E6 | Draw 0, 97.4: the 12-member ensemble | Does a draw with no drop have the same switch? |
| E7 | Draw 2, 96.65, its largest step: the ensemble and substep log | Is draw 2's −12.3% the same event? |

**Reading.** The collapse is a defect of the solve if the ledger can't account for it as physical
melt, or if one substep removes ice no flux or melt rate at that substep's length could. It's the
same event on draw 2 if that step holds one substep over 30% on the same plateau body. E6 and E7 use
the same 5-member gap.

**Not tested:** other draws' drops beyond these. If time allows, the same factorial runs on the largest step of
draw 2. Any change to the drift period, the noise or the march is proposed separately.

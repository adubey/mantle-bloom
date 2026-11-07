# Issue #811: time-scale interfaces and a step-size sweep

Status: first pass in progress. Day, 1 yr, century, 250 yr and 100 kyr modes, plus a 20 Myr bridge
at 250 kyr and 1 Myr steps, on seed 997271774 only, diagnostic wind. Everything below is one save and one wind mode; none of it is yet an
attribution for RFC 0015 or RFC 0016.

## Summary in plain terms

The question: the model runs plates, wind, rain, rivers, lakes and ice on different clocks, and RFC 0015
and RFC 0016 are progressing slowly. Does a mismatch between those clocks cause the trouble? One save
(seed 997271774) and the two RFC 0016 windows were replayed at many step lengths, from a day to 1 Myr.
Everything below is a measurement on those, with the numbers in the sections that follow.

- **One real mechanism found, and it only bites at short steps.** A lake's overflow is booked as a
  volume for the step and then turned into a flow by dividing by the step length. A day step makes
  the same release 365 times bigger as a flow than a 1 yr step does, so day steps show river
  discharge spikes of hundreds of thousands to millions of m3/s and 1 yr steps show none. The default
  100 kyr step does not see this.
- **From 100 yr steps up, total lake volume swings a lot from step to step.** But single large lakes
  flipping between full and empty are rare there: a handful per 200 steps, not the cause of most of
  the swings.
- **RFC 0016's stuck test events (windows D and E) depend on step size.** At the default 100 kyr
  down to 10 kyr the lake alternates between about 2.3M km3 and nearly empty. At 5 kyr and 2.5 kyr
  that mostly stops, and at 1 kyr it is gone. That is a step-size effect, but #754 had already traced
  these events to a lake-ice sublimation bug, so it does not show a general clock mismatch. On the
  current code the events no longer occur at all, so the windows need replacing regardless.
- **Ice (RFC 0015) is not measured yet.** The ice march (RFC 0015 D2) is not merged, so every ice
  number here is from the old one-shot solve. The RFC 0015 agent reports that the march has its own
  timescale problems (a regime change at 250 yr and a settle test that fails, #822).
- **So, are RFC 0015 and 0016 slow because of timescale mismatch?** Not established. For RFC 0016
  the evidence points at a specific bug that is already fixed, not at the clocks in general. For RFC
  0015 the march itself is a timescale problem, but I have not measured its effect.

Driver: [`bin/debug/replay_timescale_sweep.py`](../../bin/debug/replay_timescale_sweep.py), on
`origin/main` `d02484e6`. Save: `mantle-bloom-seed997271774-87400000y.mbworld` (87.4 Myr), ice
ages off, 200 steps per mode. Raw outputs and probe scripts' data are in
`~/mantle-bloom-captures/issue811/` (with `SHA256SUMS`); the compact summary is
[`seed997271774-modes.json`](seed997271774-modes.json).

## How `step_world` hands time to each process

| Process | Time it sees per step | Where |
|---|---|---|
| Plates, erosion, volcanism | the full `years` | `world.step_world_progress` |
| Wind, `diagnostic` (default) | none; recomputed from climate formulas every call | `climate.py` |
| Wind, `cfd` | a fixed 1 simulated day, whatever `years` is | `atmosphere_cfd.SECONDS_PER_TECTONIC_STEP` |
| Ocean currents | none; diagnostic, recomputed every call | `climate.compute_ocean_currents` |
| Rain, humidity | none; diagnostic climate snapshot (RFC 0019 plans transport) | `climate.py` |
| Rivers | steady-state routing of an annual-mean climate; river-water reservoirs use `years` | `hydrology.py` |
| Lakes | volume integrated over `years` from that climate; overflow booked as a step volume | `lakes.py` |
| Glaciers | one-shot solve over `years` today; RFC 0015 D2 (PR #806, not merged) marches in substeps of `min(250 yr, years)`, so steps under 250 yr stay a single solve and only 100 kyr and longer steps march (see the RFC 0015 note below) | `ice_sheet.py` |

The model has no real seasonal cycle (`climate.py`), so sub-year steps all see the same climate.

## Findings so far

1. **Overflow is a step volume turned into a rate.** `hydrology.py` does
   `overflow_source = lake_result.overflow_m3 / max(years, 1e-30)` and routes that as steady flow.
   A release of volume V therefore shows up as mean discharge V / dt. At day steps the same kind
   of release is 365 times larger than at annual steps, and nothing limits a release by the
   outlet's capacity. Day mode shows one-step spikes after the start-up (largest 3.15 M m3/s at
   step 21, with lake volume dropping 95.6k to 89.3k km3 in the same step).
2. **A start-up transient at step 2 that does not depend on step length.** Lake volume goes
   320k to ~100k km3 at step 2 in day, 1 yr and 100 kyr modes (97.5k, 100.8k, 118k). On seed
   910211954 (day steps) it goes 5.91M to 4.94M km3 and then 4.47M. In both saves about 470-500
   nodes (1.8-2.0 million km2) become ocean at step 2, and about 70-90 at step 3. The saves
   predate #758 (enclosed ocean-floor basins as ocean), so this is probably the new rule acting on
   an old save. That is a hypothesis: the newly-ocean nodes held none of step 1's lake water, and
   the probe's transfer counters are cumulative, so the cause of the dump is not identified.
3. **No strict lake period-2 alternation once the start-up is excluded, but lag-1 is
   inconclusive at 100 kyr.** From step 11 of the no-warm-up runs, lag-1 autocorrelation of the
   step change in world lake volume is -0.01 (day), 0.11 (1 yr) and -0.04 (100 kyr); the value
   over all 200 steps (-0.62 at 100 kyr) was the start-up transient. The warm-up rerun gives -0.41
   at 100 kyr (see below), so this statistic depends on the trajectory. Single-body flicker
   (#754's definition) is not yet computed by the driver.
4. **Lake refill follows physical time.** After the collapse lakes do not refill in 200 days,
   refill from 101k to 266k km3 over 200 yr, and are back to 278k km3 within one 100 kyr step.
5. **Ledgers close in every mode** (worst normalized lake and land residual 2e-7 or smaller).
6. **Not comparable across modes:** drift of ice and sea level. The spans differ (200 days, 200
   yr, 20 Myr). The 100 kyr run's 20% ice growth and ~125 m sea-level fall are evolution, not a
   step-size effect.

## Warm-up rerun (same save, 10 day-steps first, then 200 measured steps)

Compact summary: [`seed997271774-modes-warm.json`](seed997271774-modes-warm.json).

- **The step-2 transient is out of the measured series.** Lake volume starts at 95.6k (day),
  98.4k (1 yr) and 253k km3 (100 kyr).
- **The day run is deterministic.** Its measured steps match the no-warm-up run from step 11
  on, so replays are reproducible.
- **Day mode:** discharge spikes persist (largest 3.15 M m3/s, 133x between consecutive steps),
  one with a 6.3k km3 lake drop in the same step and the others with smaller lake changes. Only one
  step drops lake volume by more than 1,000 km3.
- **1 yr mode:** lakes rise smoothly from 98k to 266k km3 and discharge from 24k to 49k m3/s
  (largest consecutive ratio 1.26). No spikes. The high lag-1 (0.85) is the trend.
- **100 kyr mode:** lake volume moves 253k to 587k km3 with large step changes (74 of 199
  steps drop by more than 1,000 km3; largest drop 145k, largest rise 124k). Discharge
  stays between 33k and 147k m3/s. Plates move 100 kyr per step here, so some of this is
  real geology; I have not separated it from solver behaviour.
- **Lag-1 is not a reliable flicker statistic at 100 kyr.** It is -0.41 with warm-up and -0.04
  without (from step 11), because the warm-up changes the trajectory. The rate at which the step
  change flips sign is 0.44 and 0.49, close to the 0.48 of a random walk. Treat the earlier
  "no alternation in any mode" as "no strict alternation, lag-1 inconclusive at 100 kyr", and use
  #754's body-level definition for flicker.
- **What the three modes support:** discharge peaks scale with 1/dt (day spikes of millions of
  m3/s, none at 1 yr and 100 kyr), and the lake series behave differently by step size (a drift
  with few events, a smooth rise, large per-step swings). Which of these differences is a
  numerical artefact and which is physical time is not yet separated.

## Longer steps and the 20 Myr bridge (same save, warm-up, `origin/main` b7885e2e driver)

Compact summary: [`seed997271774-modes-long.json`](seed997271774-modes-long.json). Same 10 day-step
warm-up, diagnostic wind. Century and 250 yr cover 20 and 50 kyr; the bridge covers 20 Myr.

| Mode | Steps | Lake km3 first to last | Steps dropping lake by >1,000 km3 | Largest drop km3 | Max land discharge m3/s | Largest consecutive ratio | Ice km3 first to last |
|---|---|---|---|---|---|---|---|
| Day | 200 | 95.6k to 90.0k | 1 | 6.3k | 23.5k to 3.15M | 133 | 5.653M to 5.648M |
| 1 yr | 200 | 98k to 264k | 2 | 2.2k | 24k to 56k | 1.26 | 5.653M to 5.649M |
| Century | 200 | 246k to 351k | 67 | 36k | 33k to 99k | 1.63 | 5.653M to 5.651M |
| 250 yr | 200 | 261k to 347k | 63 | 70k | 33k to 103k | 1.91 | 5.653M to 5.654M |
| 100 kyr | 200 | 253k to 539k | 74 | 145k | 33k to 147k | 2.22 | 5.651M to 6.743M |
| 250 kyr | 80 | 230k to 540k | 36 | 65k | 35k to 122k | 2.00 | 5.664M to 6.857M |
| 1 Myr | 20 | 204k to 315k | 7 | 146k | 34k to 110k | 1.81 | 5.715M to 6.865M |

- **Big lake swings start between 1 yr and 100 yr steps.** 1 yr has 2 drops over 1,000 km3; century
  has 67 and 250 yr 63. The century and 250 yr runs cover 20 and 50 kyr, and their ice (about
  5.65M km3) and sea level (-403 m) do not move, so the swings are not the ice sheets or sea level
  changing. That points at the lake solver, but I have not excluded plate motion over those spans.
- **Swing size grows with step length:** largest drop 36k km3 (100 yr), 70k (250 yr), 145k (100 kyr),
  146k (1 Myr).
- **Discharge peaks stay modest above 1 yr:** consecutive ratio 1.26, 1.63, 1.91, 2.22 for 1 yr,
  century, 250 yr, 100 kyr. The 1/dt spikes are a sub-year effect (day: 133x).
- **Lag-1 of the lake change is negative in every long mode** (-0.30 century, -0.15 250 yr, -0.41
  100 kyr, -0.36 250 kyr, -0.50 1 Myr, the last from 19 changes), and the sign-reversal rate is
  0.50 to 0.62, against 0.48 for a random walk. That is suggestive of weak alternation, not proof;
  #754's single-body count is the real test and is not computed yet.
- **Ice and sea level:** ice grows about 20% and sea level falls 127 to 140 m over the 20 Myr in the
  100 kyr, 250 kyr and 1 Myr runs. The ice step gap needs per-checkpoint values and stays
  provisional until RFC 0015 D2 merges.
- **Ledgers close** to 4e-11 or better in the century, 250 yr and bridge runs.
- **Lake volume at the first measured step** is 98k (1 yr) against 246k to 261k (century, 250 yr):
  the lakes refill after the warm-up on physical time (finding 4), not a step-size effect.

## Single-body flicker (#754's definition) and a new baseline

Compact summary: [`seed997271774-body-flicker.json`](seed997271774-body-flicker.json) (event lists for
1 yr, 100 kyr and the bridges). Driver: `replay_timescale_sweep.py` now counts large-body reversals
through [`lake_body_tracking.py`](../../bin/debug/lake_body_tracking.py), shared with
`measure_lake_step_budget.py`. A reversal is a matched body over 10,000 km2 whose two consecutive
volume changes have opposite signs and each exceeds half its earlier volume.

**These seven runs are on `origin/main` dc855ca3, not d02484e6.** #812 (only connected-ocean water
floats ice, #552) merged in between and moved the lake and ice numbers, so the earlier tables are
not reproduced. Pattern and conclusions below use the new runs only. Same save, warm-up and wind.

| Mode | World-lake drops >1,000 km3 | Largest drop km3 | Max land discharge m3/s | Lag-1 | Body reversals | Comparable body-steps | Rate | Reversals of bodies >100,000 km2 |
|---|---|---|---|---|---|---|---|---|
| Day | 1 | 5.1k | 23.5k to 591k | 0.06 | 6 | 3,859 | 0.16% | 0 |
| 1 yr | 7 | 24k | 24k to 56k | 0.10 | 1 | 4,597 | 0.02% | 0 |
| Century | 64 | 72k | 33k to 100k | -0.24 | 29 | 4,634 | 0.63% | 8 |
| 250 yr | 67 | 32k | 33k to 91k | -0.19 | 102 | 4,710 | 2.2% | 8 |
| 100 kyr | 85 | 87k | 29k to 130k | -0.41 | 8 | 727 | 1.1% | 5 |
| 250 kyr (20 Myr) | 32 | 101k | 35k to 122k | -0.20 | 4 | 195 | 2.1% | 2 |
| 1 Myr (20 Myr) | 10 | 56k | 34k to 138k | -0.47 | 0 | 15 | 0 | 0 |

- **Most reversals are small bodies at the 10,000 km2 threshold.** Day events are 20 to 30 km3
  bodies; century and 250 yr are mostly hundreds of km3. They are not the RFC 0016 flicker.
- **Reversals of bodies over 100,000 km2 exist at every step from 100 yr up:** 8 (century), 8
  (250 yr), 5 (100 kyr), 2 (250 kyr), none at day or 1 yr. At 100 kyr the largest is a 164k km3,
  0.82 M km2 body that dropped 108k then rose 122k (step 103); at 250 kyr a 172k km3, 0.85 M km2
  body dropped 123k then rose 104k (step 67). They are a handful per 200 steps.
- **Whole-world swings are mostly not single-lake reversals.** 85 big-drop steps against 8
  reversals at 100 kyr, 64 against 29 at century. The remaining change is lakes growing, merging,
  splitting or draining and not refilling.
- **Comparison across modes is weak.** At 1 Myr only 15 body-steps are comparable, so zero
  reversals is not evidence of absence. Large-body counts are single digits.
- **Day-mode peak discharge fell from 3.15 M to 591k m3/s** on the new baseline; the spikes remain.
- Ledgers close to 3e-8 or better in every mode.

What this supports: large-lake reversals occur from 100 yr steps upward and not at day or 1 yr
steps, and they are rare. It does not show that RFC 0016's flicker is a step-size effect; that
needs the RFC's own windows replayed at several step sizes.

## RFC 0015 note: the glacier march is a timescale problem of its own

Source: the RFC 0015 agent's comment on #811, measured on the D2 branch (PR #806, issue #822) with
`generate_world(5, node_density=0.25)` and a single 250 kyr step. I have not reproduced these
numbers.

- **The march changes regime at 250 yr.** Substeps are fixed at `min(250 yr, years)`. Hour, day, 1 yr
  and century steps do one ice solve per step (endpoint-like, no march loop), a 250 yr step does
  exactly one substep, and only steps of 100 kyr and longer march (400 substeps at 100 kyr). Ice
  metrics compared across that boundary mix the solver regime with the step size, so the ice
  comparison needs a dedicated century / 250 yr / 100 kyr arm after D2 lands, not a read of the
  200-step table.
- **The settle test fails on generated worlds (#822).** Pair-volume change is about 5.65e-4 against
  a 1e-3 tolerance, but the three-consecutive-pairs counter never fills, so a 1 Myr step runs all
  4,000 substeps (about 2 min) and a 5 Myr step all 20,000 (about 10 min). The residual
  oscillation is the same two-substep alternation the damping was built to cancel. For the
  interface table: once D2 lands, the ice handoff averages over an unsettled alternation on these
  worlds, so "does the handoff conserve" has to be answered for the averaged output, not the last
  substep.
- **Cost per marched step, unsettled worst case:** about 30 ms per substep, so about 12 s at 100 kyr,
  2 min at 1 Myr and 10 min at 5 Myr. A 200-step 100 kyr arm pays about 40 min of march alone, which
  matters for the compute budget of the post-D2 rerun.
- **Correction to my earlier reading.** The comment says my century and 250 yr observation (ice and
  sea level do not move) predates the march merge at `dc855ca3`. `dc855ca3` is the merge of my own
  PR #815; D2 (PR #806) is still open and is not in `origin/main`. So every run here used the
  one-shot ice solve, and the observation is a statement about that solver. It is still true that
  the post-D2 rerun may differ: the march re-solves sea level from the current ice debit at every
  substep, which could move sea level in exactly those modes.

Plan: the ice metrics in every table here are labelled pre-march. After D2 merges, rerun the day,
1 yr, century, 250 yr, 100 kyr and bridge modes, with the century / 250 yr / 100 kyr ice comparison
reported separately, and size the run against the compute budget first.

## RFC 0016 windows D and E at six step sizes

Data: [`rfc0016-window-de-reversals.json`](rfc0016-window-de-reversals.json) (target volumes per step,
reversal counts, and the current-`main` result). Driver: `replay_rfc0016_river_water.py`, control arm,
ice-age cycle off, 0.4 Myr from each window's registered checkpoint (D: seed 910211954 at 225.1 Myr,
SHA-256 `6a8d6f90...`; E: 224 Myr, `4045b915...`), with the registered anchor files. Raw runs are in
`~/mantle-bloom-captures/issue811/rfc0016-windows-pin/`.

**Source.** RFC 0016's pinned source `0bda513d`, in my own detached worktree. At 100 kyr it
reproduces the registered traces exactly (D: 2,293,616 / 2,768 / 268,996 / 788,971 km3; E: 106,435 /
2,268,652 / 6,422 / 2,271,028). On current `origin/main` (`dc855ca3`) the same checkpoints at 100 kyr
do not: the anchor is dry at all four steps in both windows, and the nearest lake holds 341k to 551k
km3 (D) and 603k to 949k km3 (E). So the registered events do not occur on current code.

Reversal count: #754's rule on the anchor's body (over 10,000 km2, consecutive volume changes of
opposite sign, each over half the earlier volume). The body is the lake at the anchor point per step,
not overlap-matched across steps as #754 does it.

| Window | Step | Steps | Anchor dry | Reversals / comparable | Lake volume range km3 |
|---|---|---|---|---|---|
| D | 100 kyr | 4 | 0 | 1/2 | 2.8k to 2.29M |
| D | 25 kyr | 16 | 0 | 14/14 | 2k to 2.31M |
| D | 10 kyr | 40 | 13 | 9/9 | 8k to 2.53M |
| D | 5 kyr | 80 | 0 | 4/78 | 286k to 7.07M |
| D | 2.5 kyr | 160 | 0 | 1/158 | 343k to 3.60M |
| D | 1 kyr | 400 | 0 | 0/398 | 330k to 532k |
| E | 100 kyr | 4 | 0 | 2/2 | 6k to 2.27M |
| E | 25 kyr | 16 | 5 | 3/4 | 89k to 2.28M |
| E | 10 kyr | 40 | 12 | 6/12 | 61k to 2.43M |
| E | 5 kyr | 80 | 1 | 5/75 | 97k to 6.56M |
| E | 2.5 kyr | 160 | 0 | 19/158 | 38k to 4.54M |
| E | 1 kyr | 400 | 0 | 0/398 | 206k to 338k |

- **From 100 kyr to 10 kyr the flicker is strong.** D at 25 kyr alternates between about 2.3M and 10
  to 20k km3 on every step; at 10 kyr the anchor is also dry in 13 of 40 steps (D) and 12 of 40 (E),
  so the lake appears and disappears as well.
- **At 5 kyr and 2.5 kyr the regular alternation is gone, but D and E differ.** D falls to 4/78 and
  then 1/158. E keeps scattered reversals (5/75, then 19/158 at steps 10 to 94), not a period-2
  pattern. At 1 kyr both have none in 398 steps; step changes there are 10 to 40% of volume (up to 89k
  and 128k km3), lag-1 -0.26 (D) and -0.69 (E), below the rule's 50% bar.
- **The lake itself changes with step size, not just its noise.** The anchor lake peaks at 2.3M km3
  (100 kyr to 10 kyr), 7.07M and 6.56M (5 kyr), 3.60M and 4.54M (2.5 kyr), and 0.53M and 0.34M
  (1 kyr). Each step size is a different state of the basin, so the traces are not one lake with
  different noise. I have not examined why 5 kyr grows so large.
- **Attribution (c), for the stuck window E trace and D:** the failure is strong at 100 kyr to 10
  kyr, weaker and irregular at 5 kyr and 2.5 kyr, and absent at 1 kyr. That is a step-size
  dependence. It does not say the cause is a timescale mismatch between processes.
- **#754 already identified a cause for D and E on this source:** lake-ice sublimation charged after
  the lake solver, outside its water balance, with an integrator that switches branch on a lagged
  surface temperature (#757 moves it into the solver; #758 reclassifies the basin). The step-size
  result is consistent with that artifact depending on step length. I have not tested it; the
  absence of the events on current `main` is consistent with #757 and #758 having removed it.
- **Limits.** One trajectory per step size, 0.4 Myr, ice-age cycle off, no 100 yr or 1 yr run. The
  states differ across step sizes, so "gone at 1 kyr" does not show finer steps remove the flicker
  for the same reason on a comparable state. Window A is not run.

## Method change: warm up first

Because of finding 2, every mode now runs the same `--warmup-days` day-length steps (default 10)
before measuring, and the flicker statistics use step 11 onward. The first-pass JSON here was
taken without warm-up and is kept for the record.

## Open

- Identify the step-2 body and why it disappears (cap level, tier, rim, `#758` hysteresis).
- Hour and `cfd` modes. RFC 0016 window A, and why the lake peaks at 5 kyr.
- Second save at every mode.
- Rerun the ice metrics once RFC 0015 D2 merges (PR #806; march cost and regime boundary in the RFC 0015 note).
- The #553 depth-cap warning on seed 910211954 at step 1 (6.7e4 km3 gross overflow on lake-tier
  bodies, 3.5e4 on sea-tier bodies) is another large-overflow route.

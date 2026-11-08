# Issue #811: time-scale interfaces and a step-size sweep

Status: first pass in progress. One save (seed 997271774) in diagnostic and `cfd` wind, from day to
1 Myr steps, plus the three RFC 0016 windows at six step sizes. None of it is yet a full attribution for RFC 0015
or RFC 0016, and it predates the RFC 0015 D2 glacier march (merged 2026-11-07).

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
  these events to a lake-ice sublimation bug, so it does not show a general clock mismatch. Window A, the third gate window, shows no such swings at any step size. On the
  current code the D and E events no longer occur at all, so those windows need replacing regardless.
- **Ice (RFC 0015) is not measured yet.** The ice march (RFC 0015 D2) merged on 2026-11-07 (PR #806),
  after every run here, so every ice number is from the old one-shot solve. The RFC 0015 agent reports that the march has its own
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
  PR #815; D2 (PR #806) was still open then. It merged on 2026-11-07, after every run here. So every run here used the
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
  for the same reason on a comparable state. Window A is covered separately below.

### Window A (seed 910211954, 221.6 Myr, checkpoint `9882568e...`)

Same source, driver and rule. At 100 kyr the step-2 target volume is 863,388.83 km3, matching the
registered G0 value in `analysis/rfc0016/g0-liquid-summary.json`, and the anchor is dry at step 1
as the registered account says. I did not recompute the +495,624 km3 event change itself.

| Step | Steps | Anchor dry | Reversals / comparable | Dry counted as zero | Lake volume range km3 |
|---|---|---|---|---|---|
| 100 kyr | 4 | 1 | 0/1 | 0/2 | 783k to 863k |
| 25 kyr | 16 | 3 | 0/8 | 3/14 | 750k to 894k |
| 10 kyr | 40 | 9 | 0/17 | 12/38 | 761k to 899k |
| 5 kyr | 80 | 42 | 0/36 | 0/78 | 863k to 957k |
| 2.5 kyr | 160 | 91 | 0/67 | 0/158 | 921k to 985k |
| 1 kyr | 400 | 219 | 0/179 | 0/398 | 965k to 1.01M |

- **A has no strict reversals at any step size.** Its lake stays between about 750k and 1.0M km3, with
  no 2.3M to near-zero swings. This matches the RFC's own note that A's registered event survives
  once and is not a repeating flicker.
- **A mild alternation at 25 kyr and 10 kyr** shows only when dry steps count as zero (3/14, 12/38):
  the 25 kyr trace alternates between about 750k and 890k km3. It is gone at 5 kyr and below.
- **The anchor is dry in many fine-step outputs** (42/80, 91/160, 219/400) while the lake holds
  nearly 1M km3, so this is probably the anchor point sitting at the lake's edge, not the lake
  disappearing. I have not checked.
- **The three windows differ.** D and E show large 2.3M to near-zero swings through 10 kyr; A shows
  none. That fits per-basin causes such as #754's better than one general clock mismatch.

## CFD wind against diagnostic wind

Compact data: [`seed997271774-cfd.json`](seed997271774-cfd.json). The same seven modes as the
diagnostic set, with `--wind cfd`, on `dc855ca3` in a detached baseline worktree (the diagnostic
runs' commit, because #819 has since changed bathymetry on `main`), same 10 day-step warm-up. Smoke
check: the CFD atmosphere advances 86,400 s per step.

| Mode | Lake km3 first to last (CFD / diag) | Median max land discharge m3/s (CFD / diag) | Steps dropping lake >1,000 km3 (CFD / diag) | Body reversals (CFD / diag) |
|---|---|---|---|---|
| Day | 288k to 272k / 100k to 94k | 103k / spiky | 3 / 1 | 0 / 6 |
| 1 yr | 289k to 443k / 102k to 280k | 103k / 49k | 0 / 7 | 0 / 1 |
| Century | 369k to 916k / 252k to 314k | 168k / 72k | 23 / 64 | 13 / 29 |
| 250 yr | 387k to 884k / 261k to 348k | see json | 52 / 67 | 21 / 102 |
| 100 kyr | 386k to 582k / 253k to 547k | 184k / 87k | 99 / 85 | 32 / 8 |
| 250 kyr (20 Myr) | 383k to 327k / 230k to 354k | see json | 40 / 32 | 9 / 4 |
| 1 Myr (20 Myr) | 353k to 365k / 204k to 245k | see json | 7 / 10 | 6 / 0 |

- **The CFD clock does not follow the step.** It advances one simulated day per step in every mode,
  so each 200-step run sees 200 atmosphere days; the 250 kyr and 1 Myr bridges see 80 and 20. It is
  real time only in the day mode.
- **Global means agree, the hydrology does not.** Mean precipitation is within about 0.3% of
  diagnostic (1198 against 1201 mm in day mode) and mean temperature within 0.1 C. But CFD starts
  wetter (lake volume 289k to 387k against 100k to 261k km3) and its median largest river is about
  twice diagnostic's (103k to 184k against 49k to 87k m3/s). So the difference is in where rain
  falls, not how much.
- **The 1/dt mechanism, quantified.** In CFD day mode, step 8 has a largest land discharge of 124 M
  m3/s, far beyond any real river. Times the 86,400 s step that is about 10,720 km3, and the world
  lake volume drops 11,060 km3 in the same step: a one-step lake release booked as a one-day flow.
  This is the clearest case so far of the overflow-as-rate behaviour in finding 1, 200 times the
  diagnostic day peak.
- **Large-lake reversals are mixed:** more with CFD at 100 kyr and the bridges (32/931 against
  8/727; 9/251 against 4/195; 6/38 against 0/15), fewer at century and 250 yr. Samples are small.
- **Ledgers close to 3e-7 or better in every CFD mode** (day mode loosest, 2.8e-7).
- **What this does not show.** Whether CFD wind's larger rivers are more or less right, or whether the
  step-size sensitivity differs between wind models in any way beyond the day spike. One save. The post-RFC 0019 rain arm is not run.

## Hour mode (1 hour steps, both wind models)

Compact data: [`seed997271774-hour.json`](seed997271774-hour.json). Same save, same `dc855ca3`
baseline worktree, same 10 day-step warm-up, 200 measured steps of 1.1408e-4 yr (8.3 days in
all). Compared with day mode from [`seed997271774-cfd.json`](seed997271774-cfd.json).

| | Diagnostic hour | Diagnostic day | `cfd` hour | `cfd` day |
|---|---|---|---|---|
| Lake km3, first to last | 99.6k to 92.5k | 99.6k to 94.0k | 288k to 272k | 288k to 272k |
| Steps dropping lake >1,000 km3 | 1 | 1 | 3 | 3 |
| Largest one-step lake drop, km3 | 5,136 | 5,112 | 11,052 | 11,060 |
| Max land discharge, m3/s (median / peak) | 48k / 3.5M | 26k / 591k | 138k / 2.30G | 103k / 124M |
| Mean precipitation, mm | 1199.8 | 1200.8 | 1196.7 | 1197.7 |
| Mean air temperature, C | 16.72 | 16.72 | 16.71 | 16.71 |
| Large-body reversals | 0 of 4,019 | 6 of 3,859 | 0 of 4,000 | 0 of 4,286 |
| Worst lake residual | 4e-8 | 3e-8 | 5e-6 | 3e-7 |
| CFD atmosphere days, first to last | | | 11 to 210 | 11 to 210 |

- **Hour mode behaves like day mode.** The lake, precipitation and temperature numbers match
  day mode to within a few percent, and the same lake drops occur at the same size: 5.1k km3 in
  diagnostic and 11k km3 in `cfd`. A 24x shorter step does not make the lake solver unstable or
  change the climate.
- **The `cfd` clock is 24x too fast here.** Each hourly step advances the atmosphere one
  simulated day, so 200 steps cover 8.3 days of world time but 200 atmosphere days. This is the
  mismatch the handoff table predicted.
- **The 1/dt mechanism scales with the step.** In `cfd` the early lake release (step 8, with a
  follow-on at step 9) shows up as a largest land discharge of 2.30 G m3/s at step 8 and 0.44 G
  at step 9. Times the 3,600 s step that is about 8,270 km3 and 1,580 km3, together about 9,850
  km3, against lake drops of 11,052 and 2,819 km3 (13,871 km3) in those steps. So the
  discharge spike accounts for about 70% of the drop; the rest is outside the single largest
  cell. In day mode the same release gave 124 M m3/s, which is 97% of the drop. The peak grows
  with 1/dt (2.3 G at one hour against 0.12 G at one day, 19x for a 24x shorter step), as
  `overflow / max(years, 1e-30)` predicts. It is a bookkeeping rate, not a river.
- **The two diagnostic-run features are explained in the next section.** The 3.5 M m3/s peak at
  step 110 is a 12.7 km3 lake release booked as a one-hour flow (an earlier version of this
  note said it was not a release; that was wrong), and the 5,136 km3 drop at step 13 is a lake
  that reached the ocean.

- **No large-lake reversals at all** in either hour run. Ice moves under 0.1% and sea level
  0.5 m (diagnostic) or 0.03 m (`cfd`) over the run.
  The sampled lake window is 8.3 days, so this says little about flicker, which needs
  decades.
- **Ledgers:** diagnostic closes to 4e-8. `cfd` hour closes to 5e-6, worst at step 1 and about
  1e-6 afterward, looser than day mode (3e-7). I have not looked into why.
- **Runtime:** 6.1 s per step diagnostic and 8.7 s `cfd`, 200 steps in 22 and 30 minutes.
  Peak memory 1.7 GB.
- **What this does not show.** It is one save. Hour steps do not reach the timescales where
  the problems showed up (century and longer for lakes, 250 yr and longer for the glacier
  march), so it clears the short-step end of the range and nothing more. Ice is the old
  one-shot solver at every step size.

## Why the D/E anchor lake is largest at 5 kyr, and why A's anchor goes dry

Compact data: [`rfc0016-window-de-peak.json`](rfc0016-window-de-peak.json). Same pinned source
(`0bda513d`) and checkpoints as the windows section above. Probe scripts:
`bin/debug/trace_lake_anchor.py` (single steps) and `bin/debug/probe_anchor_shore.py` (new).

### D and E: two lake states, and the step size decides which one the run is in

Measuring the anchor lake's catchment per step shows it is in one of two states. In the small
state the catchment is 3 to 8M km2 and the lake 0.3 to 1.5M km3 over 1 to 1.8M km2. In the merged
state the catchment is 19 to 25M km2 and the lake 3 to 7M km3 over 6 to 9M km2: the anchor basin
has joined its neighbours into one inland sea. The 5 kyr "peak" in the earlier table is the
merged state.

| Step | D steps in merged state | E steps in merged state | Largest world lake, M km3 (D / E) |
|---|---|---|---|
| 1 kyr | 0 of 400 | 0 of 400 | 3.85 / 3.14 |
| 2.5 kyr | 145 of 160 | 126 of 160 | 5.87 / 6.31 |
| 5 kyr | 72 of 80 | 68 of 80 | 9.03 / 8.27 |
| 10 kyr | 1 of 40 | 0 of 40 | 4.16 / 4.22 |
| 25 kyr | 0 of 16 | 0 of 16 | 4.15 / 4.04 |
| 100 kyr | 0 of 4 | 0 of 4 | 4.12 / 4.02 |

The merged state is entered at 2.5 and 5 kyr, in the first 20 to 30 kyr, and almost never at the
other sizes. At matched times up to 15 kyr all of 1, 2.5 and 5 kyr are in the small state, with the
longer steps' lake larger (D at 15 kyr: 0.49, 0.96 and 1.44M km3).

### Same start, one step: the lake after one step depends on the step length

D's checkpoint has the anchor lake at 0.02M km3 (E's is dry). One step of each size from that same
state, climate computed once from it:

| Step | D anchor lake after one step, M km3 | World lake, M km3 | Lake merges/splits in the step |
|---|---|---|---|
| 1 kyr | 0.53 | 2.28 | 81 / 72 |
| 2.5 kyr | 0.34 | 2.28 | 92 / 84 |
| 5 kyr | 0.62 | 2.47 | 93 / 86 |
| 10 kyr | 1.19 | 3.02 | 98 / 93 |
| 25 kyr | 2.29 | 4.12 | 108 / 100 |
| 100 kyr | 2.29 | 4.12 | 116 / 107 |

The 25 and 100 kyr lakes are identical, a basin filled to its cap. The ledger residual is below
1e-11 in every case, so this is not a water-balance error. In every case the hierarchy merges and
splits 64 to 116 bodies in one step.

### Control: shorten the lake solver's internal slices, nothing else

`solve_lakes_prognostic` has a `substeps` argument that `step_world` leaves at 1, so a 100 kyr
step is one implicit slice. I forced more slices for the same single step (climate unchanged):

| Step | One slice | Shorter slices |
|---|---|---|
| 5 kyr | 0.62M km3 | 0.56 (1 kyr slices), 0.29 (200 yr slices) |
| 25 kyr | 2.29M | 1.46 (1 kyr slices) |
| 100 kyr | 2.29M | 0.85 (1 kyr slices) |

The answer moves every time the slice shortens, and it has not converged at 200 yr. So the lake
state depends on the slice length of the solver itself, within one step, at fixed climate. It
fits the D-and-E picture: longer slices fill further, and a lake that is filled further can
reach a saddle and merge.

### What this does and does not show

- **Shown:** the 5 kyr peak is the merged state, which only 2.5 and 5 kyr runs enter; one step
  from identical state gives 0.34 to 2.29M km3 depending on the step; shortening the solver's
  slices at fixed climate changes the result.
- **Not shown:** why the merged state is entered at 2.5 and 5 kyr but not at 10 kyr. E at 100 kyr
  dries out (0.11M km3) where E at 10 and 25 kyr holds 1.3 to 1.4M, so the order is not
  monotone. I did not separate the slice length from the climate refresh interval: the single
  steps with shorter slices hold the climate fixed, while the 1 kyr chain refreshes it every
  kyr, and the 1 kyr chain stays near 0.5M km3 where the 1 kyr-slice controls at 25 and 100 kyr
  reach 1.46 and 0.85. I also did not test whether the merged state is stable if started there
  at 1 kyr.
- **Relation to #754:** #754 traced the registered D and E reversals to lake-ice sublimation
  charged outside the lake balance (fixed by #757 and #758). That is a separate effect. This
  one is that the size and the merge state of the lake depend on the step length even with
  the ledger closing.

### A: the dry steps are a thin rim cell, not a vanishing lake

The window A anchor is a single cell. A 5 kyr replay with the new probe, which keeps the nearest
lake body whether or not it contains the anchor, shows a body of 0.86 to 0.95M km3, 0.66 to 0.68M
km2 and 1.3 km mean depth at every step, including the dry ones. The anchor cell is on its rim:
5.7 to 7.2 m of water through step 35, then 1.74 m at step 36, a level drop of about 5.5 m
in one step, and below the 1 m visibility threshold (`LAKE_MIN_VISIBLE_DEPTH_M`) from step 40
onward, while the ground under it rises 0.1 to 0.3 m a step. So the "dry" steps in A are the
visibility threshold at the lake's edge. It does not change the earlier reading that A has no
reversals. Two things I noticed and did not chase: the body's volume alternates every step
(0.953 and 0.887M km3 from step 28 to 34, a period-2 swing of about 7%), and the 5.5 m level
step at step 36. One window, one step size.

## The hour run's two odd features, and a start-up transient keyed to the call count

Compact data: [`seed997271774-hour-features.json`](seed997271774-hour-features.json). New probe:
`bin/debug/probe_hour_features.py`. Diagnostic wind, the hour-mode setup, `dc855ca3`.

- **Step 110, the 3.5 M m3/s peak, is a lake release.** The lake ledger shows 12.73 km3 overflowing
  out of lakes and 11.78 km3 leaving storage, in a 3,896 km2 lake of 12.0 km3 that empties. The
  largest river carries 12.7 km3 over the 3,600 s step, so it is the overflow-as-rate behaviour
  (finding 1) at a small scale.
- **Step 13, the 5,136 km3 drop, is a lake reaching the ocean.** Two bodies vanish: one of
  5,024 km3 over 100,902 km2 (mean depth 50 m) and one of 106 km3, and 97% and 100% of their former
  nodes are ocean afterwards. The ledger books 5,151 km3 of overflow to the ocean and a
  storage change of -5,136 km3. The largest river carries only 5.9 km3 that step, so the water did
  not travel by a river. The sea level rises 14.3 mm, which over a 3.6e8 km2 ocean (my figure, not
  read from the world) is 5,150 km3. The step before (step 12) has a 373.6 km3 overflow to the
  ocean and 2,193 km3 of lake-to-lake transfer, with a 1,041 km3 body falling to 106 km3.
- **Both land on the same call number whatever the step length.** In diagnostic wind the first big drop
  is at measured step 13 for hour steps (5,136 km3), day steps (5,112 km3) and 1 yr steps (3,085
  km3). In `cfd` wind the drops are at steps 8, 9 and 10 for both hour and day: 11,052, 2,819 and
  2,394 km3 against 11,060, 2,821 and 2,394 km3. A 24x change in step length does not move them.
  The sea level series agrees to within 2 mm step by step between hour and day over the first
  16 steps (-402.570, -402.525, -402.502, -402.554 for steps 1 to 4 in both).
- **It is the call count since the save was loaded, not the warm-up.** I changed the warm-up
  from 10 to 0, 5 and 20 day steps (diagnostic, hour steps). The 5,151 km3 overflow to the ocean
  happens at call 23 after load in every case (measured step 13, 23, 18 and 3), the 374 km3
  overflow at call 22, and a smaller 784 km3 one at call 20 in the 5-step case. So it is a
  deterministic start-up sequence, measured in calls, that does not scale with the step.
- **What I did not establish.** What carries the sequence from call to call. A hypothesis that fits:
  a spill chain between lakes that moves one link per call, so that it takes about 23 calls to
  reach the ocean. I did not test it, and I found no call counter in `world.py`, `hydrology.py` or
  `lakes.py` (grep for step counts and modulo cadences). I also did not run it on a save made
  after #758, so it may be specific to this older save.
- **Why it matters for #811.** The earlier modes' first 10 to 30 measured steps are partly this
  transient, in every step size. At 100 kyr, 23 calls is 2.3 Myr of simulated time; at one hour
  it is 23 hours. If it is a spill chain moving one link per call, it is a real
  time-scale mismatch (water crossing a chain of lakes at a rate set by the step count, not by
  the step length), and it would act in the long-step modes as well as the short ones. That
  is a hypothesis, not a result.

## Method change: warm up first

Because of finding 2, every mode now runs the same `--warmup-days` day-length steps (default 10)
before measuring, and the flicker statistics use step 11 onward. The first-pass JSON here was
taken without warm-up and is kept for the record.

## Open

- Identify the step-2 body and why it disappears (cap level, tier, rim, `#758` hysteresis).
- Why the merged D/E lake state is entered at 2.5 and 5 kyr but not at 10 kyr, and whether the slice length or the climate refresh interval drives it; the period-2 volume swing and the 5.5 m level step in A.
- Test whether a spill chain moving one link per call explains the call-23 start-up transient, and rerun on a post-#758 save.
- Second save at every mode.
- Rerun the ice metrics once RFC 0015 D2 merges (PR #806, merged 2026-11-07; march cost and regime boundary in the RFC 0015 note).
- The #553 depth-cap warning on seed 910211954 at step 1 (6.7e4 km3 gross overflow on lake-tier
  bodies, 3.5e4 on sea-tier bodies) is another large-overflow route.

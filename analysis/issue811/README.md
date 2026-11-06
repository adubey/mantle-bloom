# Issue #811: time-scale interfaces and a step-size sweep

Status: first pass in progress. Day, 1 yr, century, 250 yr and 100 kyr modes, plus a 20 Myr bridge
at 250 kyr and 1 Myr steps, on seed 997271774 only, diagnostic wind. Everything below is one save and one wind mode; none of it is yet an
attribution for RFC 0015 or RFC 0016.

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
| Glaciers | one-shot solve over `years` today; RFC 0015 D2 would march 250 yr substeps | `ice_sheet.py` |

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

## Method change: warm up first

Because of finding 2, every mode now runs the same `--warmup-days` day-length steps (default 10)
before measuring, and the flicker statistics use step 11 onward. The first-pass JSON here was
taken without warm-up and is kept for the record.

## Open

- Identify the step-2 body and why it disappears (cap level, tier, rim, `#758` hysteresis).
- Hour and `cfd` modes.
- Second save at every mode; body-level flicker counts.
- Rerun the ice metrics once RFC 0015 D2 merges.
- The #553 depth-cap warning on seed 910211954 at step 1 (6.7e4 km3 gross overflow on lake-tier
  bodies, 3.5e4 on sea-tier bodies) is another large-overflow route.

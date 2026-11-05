# Issue #811: time-scale interfaces and a step-size sweep

Status: first pass in progress. Day, 1 yr and 100 kyr modes on seed 997271774 only, diagnostic
wind, no warm-up. Everything below is one save and one wind mode; none of it is yet an
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
3. **No lake period-2 alternation in any mode once the start-up is excluded.** From step 11,
   lag-1 autocorrelation of the step change in world lake volume is -0.01 (day), 0.11 (1 yr)
   and -0.04 (100 kyr). The values over all 200 steps (-0.62 at 100 kyr) were the start-up
   transient. This does not rule out single-body flicker (#754's definition), which the driver
   does not yet compute.
4. **Lake refill follows physical time.** After the collapse lakes do not refill in 200 days,
   refill from 101k to 266k km3 over 200 yr, and are back to 278k km3 within one 100 kyr step.
5. **Ledgers close in every mode** (worst normalized lake and land residual 2e-7 or smaller).
6. **Not comparable across modes:** drift of ice and sea level. The spans differ (200 days, 200
   yr, 20 Myr). The 100 kyr run's 20% ice growth and ~125 m sea-level fall are evolution, not a
   step-size effect.

## Method change: warm up first

Because of finding 2, every mode now runs the same `--warmup-days` day-length steps (default 10)
before measuring, and the flicker statistics use step 11 onward. The first-pass JSON here was
taken without warm-up and is kept for the record.

## Open

- Identify the step-2 body and why it disappears (cap level, tier, rim, `#758` hysteresis).
- Hour, 250 yr, century, bridge (100 kyr / 250 kyr / 1 Myr over 20 Myr) and `cfd` modes.
- Second save at every mode; body-level flicker counts.
- Rerun the ice metrics once RFC 0015 D2 merges.
- The #553 depth-cap warning on seed 910211954 at step 1 (6.7e4 km3 gross overflow on lake-tier
  bodies, 3.5e4 on sea-tier bodies) is another large-overflow route.

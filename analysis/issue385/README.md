# Issue #385: craton resistance for failed rifts (RFC 0003 solution 2)

## Method

`bin/debug/compare_cratons.py --seeds <S> --variants on --years 1e9 --every 100e6`, seeds 7,
11 and 2024, default density, run one at a time. The baseline is the phase 1 gate run
(#378), whose branch only added diagnostics, so it is origin/main's model. Numbers are in
[replay_1gyr.json](replay_1gyr.json).

## Results (0–1 Gyr)

| Seed | Craton rifted, % of start craton (main → branch) | Expected from failed-rift loss alone | Splits | Failed rifts | Mean plates |
|---|---|---|---|---|---|
| 7 | 4.05% → 3.40% | 3.54% | 33 → 33 | 29 → 27 | 17.5 → 17.5 |
| 11 | 4.35% → 5.39% | 3.89% | 37 → 28 | 31 → 29 | 17.2 → 18.1 |
| 2024 | 10.71% → 8.64% | 9.94% | 32 → 32 | 28 → 26 | 18.6 → 15.5 |
| Mean | 6.37% → 5.81% | 5.79% | | | |

"Expected" removes 75% of origin/main's `rifted_failed_m3`, the most a full-strength
craton resists. Ledgers balance on every run (largest error 0.79 km³ on seed 2024, 5e-10 of
its initial craton).

## Reading

- Failed rifts hit cratons less often from the first one, so trajectories diverge within
  the run. Per-seed numbers mix the direct effect with a different plate history. Seed 11
  shows this most: it took 28 splits instead of 37, delaminated 153 M km³ of craton instead
  of 52, and rifted more.
- Across the three seeds, rifting fell by 0.56 points per Gyr, close to the 0.58 the
  failed-rift path alone predicts. That's the expected size: this is a consistency fix on
  10–17% of the sink.
- Failed-rift counts barely move (2 fewer per seed). Splits and mean plate count hold on
  seed 7, but seed 11 lost 9 splits and seed 2024's mean plate count fell 17%. With one
  run per seed, those can't be separated from trajectory divergence.
- As RFC 0003 expects, this alone doesn't reach 2%/Gyr. Phase 1b decides whether steering
  the cut comes next.

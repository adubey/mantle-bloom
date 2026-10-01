# Issue #259: quad sea-level jitter from represented-area churn

Quad retreat (`quad_tectonics.deform -> _retreat`) removes cells every step. Gap filling,
which grows neighbouring plates back into the vacated territory, ran only every
`gaps.GAP_FILL_INTERVAL_STEPS` (4). The area the quad surface represents therefore dropped by
about 0.5% of the sphere per step and recovered on every fourth step. The eustasy solve
conserves ocean volume, so sea level jumped with it. The fix (`gaps.gap_fill_due`) runs gap
filling on every step for cell surfaces. Line worlds keep the interval.

## Seed 1, density 1, 40 × 1 Myr steps

`bin/debug/measure_quad_area_churn.py`, second half of the run:

| | std Δ represented area / sphere | std Δ sea level | area range |
|---|---:|---:|---|
| gap fill every 4 steps (before) | 0.00747 | 46.94 m/step | 0.988 – 1.002 |
| gap fill every step on quad (after) | 0.00064 | 3.32 m/step | 0.999 – 1.001 |

## Long preset, S2 `stability_sea_level_m`

`long/` holds the quad runs from
`paired --preset long --seeds 1,2,3,4,5 --jobs 10 --out ../analysis/issue259/long`. The run
was stopped at a 2-hour limit. Quad seeds 1, 2, 4 and 5 reached 400 Myr. Quad seed 3 and
all the line runs did not, so their partial output is not kept. The fix does not change line
worlds, so the line baseline is taken from `../issue249-campaign/long`.

| seed | lines | S2 limit (2 × l + 1 m) | quad before (#249) | quad after | |
|---|---:|---:|---:|---:|---|
| 1 | 10.0 | 21.1 | 58.2 | 5.5 | pass |
| 2 | 8.3 | 17.5 | 76.3 | 5.0 | pass |
| 4 | 9.9 | 20.8 | 72.0 | 5.5 | pass |
| 5 | 7.3 | 15.5 | 73.9 | 6.4 | pass |

Every other S2 metric stayed within its limit except seed 5's `air_temperature_mean_c`. It was
0.65 against a 0.60 limit, down from 2.10 before the fix. The #249 campaign already warns on
this metric (seed 3).

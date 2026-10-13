# Issue #363: lake silt from real-unit inflow

Lake silt is now sized from the lake's real inflow (`LakeClimate.inflow_m3_per_yr`, m³/yr),
the same water that sets its level. It used to be sized from the legacy model-unit
`water_deposited`, which includes the breach surge. `SILT_ACCUMULATION_COEFFICIENT` is now a
volume fraction (m³ of silt per m³ of inflow), calibrated here so that the same amount of
silt reaches lake beds.

## Replay

Seed `997271774` from the 87.4 Myr save, 100 × 100 ky steps to 97.4 Myr (save SHA-256
`08ad9461…ed762b`, as in [issue #355](../issue355/README.md)), at the default node density
(4.0). Runs, one at a time:

- **main**: `origin/main` at `93ab3bb` (#355 steps 1–2 merged).
- **4.8e-6**: this branch with the first calibration, the coefficient that matches main's
  silt budget before each lake's room clamp.
- **branch**: this branch with the final coefficient, 4.1e-6, which matches the silt that
  actually lands.

```bash
cd backend
.venv/bin/python ../bin/debug/measure_lake_silt.py \
  ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 10 --checkpoint-myr 1 \
  --backend ../../mantle-bloom-origin-main/backend --out main.json
.venv/bin/python ../bin/debug/measure_lake_silt.py \
  ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 10 --checkpoint-myr 1 --out branch.json
```

The compact per-checkpoint comparison is
[`replay-seed997271774.json`](replay-seed997271774.json). The full replay JSONs aren't
checked in.

## Calibration

On main, `measure_lake_silt.py` taps `lakes._water_balance` and records each lake's silt
budget before the room clamp, both under main's rule and as real inflow × step years. Over
the 10 Myr, the ratio of the two totals is 4.76e-6. The per-Myr ratio ranges from 3.1e-6 to
9.4e-6. That ratio, rounded to 4.8e-6, matched the budget (128.6k vs 126.2k km³) but laid down
16% more silt (109.1k vs 94.1k km³). Most silt now goes to spilling lakes, which are full,
so the room clamp rarely binds there. Scaling by 94.1/109.1 gives 4.1e-6. That's about
7 mg/L at a bulk density of 1,600 kg/m³, far below a real river's suspended load. Erosion
delivers routed river sediment into lakes separately, so this coefficient wasn't set from
real loads.

## Results, 10 Myr totals and end point

| Measure | main | 4.8e-6 | branch (4.1e-6) |
|---|---:|---:|---:|
| Silt laid on lake beds, 10 Myr | 94.1k km³ | 109.1k km³ | 93.4k km³ |
| …budget into closed lakes (before clamp) | 34.7k km³ | 14.0k km³ | 12.1k km³ |
| …budget into spilling lakes (before clamp) | 91.4k km³ | 114.8k km³ | 96.4k km³ |
| Lakes, mean of 10 checkpoints | 25.0 | 26.3 | 23.0 |
| Lake area, mean of 10 checkpoints | 1.53 M km² | 1.47 M km² | 1.56 M km² |
| Lake volume, mean of 10 checkpoints | 140k km³ | 138k km³ | 140k km³ |
| Lakes at 97.4 Myr | 23 | 22 | 21 |
| Lake area at 97.4 Myr | 1.01 M km² | 1.12 M km² | 0.86 M km² |
| Land at 97.4 Myr (% of area) | 28.99% | 29.08% | 29.09% |
| Step time (mean) | 5.46 s | 5.31 s | 6.94 s* |
| Peak RSS | 1,999 MiB | 2,461 MiB | 2,510 MiB |

\* The branch run's steps took 10.9–11.6 s for 94–96 Myr and 4.8–5.8 s everywhere else.
The 4.8e-6 run, with identical code apart from the coefficient, averaged 5.3 s. That points
to load from elsewhere on the machine, not the change.

Lakes are counted as connected bodies of land nodes over `LAKE_MIN_VISIBLE_DEPTH_M`.

**What moved.** Total silt is unchanged by construction (−0.7%), but where it lands has
shifted. Under main's rule, a spilling lake passed its water on, so it got silt mainly from
the breach surge, which `water_deposited` carried to the sinks downstream. Closed lakes got
their settled water's share. Now every lake silts from its real inflow, a spilling lake from
the flow through it. Closed lakes' share falls from 28% to 11% of the budget. Real
through-flow lakes trap most of their inflow's sediment, so this direction is physical. But
small closed catchments now get about 65% less silt. Silt in those catchments is what
collapsed #117's lake-hierarchy cascades. This replay doesn't measure hierarchy depth.
`app.hydroclimate.lake_hierarchy_diagnostics` is the tool to check it if deep chains
reappear.

Lake count, area and volume move both ways from one checkpoint to the next, with no trend
over the 10 Myr. Mean volume differs by 0.3%. Land differs by 0.1 points, inside the ±2%
run-to-run noise (#315).

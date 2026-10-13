# Issue #355: river discharge in real units

## Step 2: route physical water once

Replay of seed `997271774` from the 87.4 Myr save, 100 × 100 ky steps to 97.4 Myr (save
SHA-256 `08ad9461…ed762b`, as in [issue #353](../issue353/README.md)). Two runs, one after
the other, each alone:

- **main**: `origin/main` at `b474481`.
- **branch**: `claude/issue-355-route-water-once` at `4512bb6`. Step 1 (#361), which it
  builds on, is behaviour-neutral.

```bash
cd backend
.venv/bin/python ../bin/debug/measure_ice_sheet.py \
  ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 10 --checkpoint-myr 1 \
  --backend ../../mantle-bloom-origin-main/backend --out main.json
.venv/bin/python ../bin/debug/measure_ice_sheet.py \
  ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 10 --checkpoint-myr 1 --out branch.json
```

The compact per-checkpoint comparison is
[`replay-seed997271774.json`](replay-seed997271774.json); the full replay JSONs are not
checked in.

| Measure | main | branch |
|---|---:|---:|
| Frozen river ice, 10 Myr total | 1.30 M km³ | 0.31 M km³ |
| …of which shed at the depth cap | 0.61 M km³ | 0 |
| Snow accumulation, 10 Myr total | 1.78 M km³ | 1.54 M km³ |
| Surface melt, 10 Myr total | 1.57 M km³ | 1.05 M km³ |
| Ice volume at 97.4 Myr | 6.03 M km³ | 5.36 M km³ |
| Ice cap at 97.4 Myr (visible ice, % of land) | 5.43% | 5.15% |
| Glacial erosion, 10 Myr total (raw capacity) | 6.24 M km³ | 5.05 M km³ |
| Land at 97.4 Myr (% of area) | 29.04% | 28.99% |
| Step time (mean) | 4.95 s | 4.95 s |
| Worst ice-ledger residual (any step) | 3 m³ | 3 m³ |
| Peak RSS | 2,531 MiB | 2,568 MiB |

**What moved.** Main froze 0.27–0.47 M km³ of river ice per Myr over the first 3 Myr, and
0.61 M km³ of it overran `GLACIER_MAX_DEPTH_M` and was shed by the cap guard. The branch
froze 0.04 M km³ per Myr in those windows and never reached the cap. Main's river ice was
sized from the legacy `water_deposited`, which also carries the breach surge
(`LAKE_BREACH_EROSION_COEFFICIENT` × lake node count) and precipitation before Budyko
evaporation. A surge that reaches a frozen sink was the likely source of the burst, but this
replay didn't isolate it. After 91 Myr the two runs freeze river ice at similar rates (main
1.1–1.9× the branch), as expected when the old stacked mm/yr and the new catchment volume
over cell area are about the same size.

The ice sheet is 11% smaller at 97.4 Myr, and the gap opens in those first 3 Myr. That
means less ice to melt (surface melt −33%) and to erode with (raw glacial erosion −19%).
Snowfall accumulation also falls by 14%, because there's less high, cold ice surface. Land
area differs by 0.05 points, inside the ±2% run-to-run noise (#315). Runtime doesn't change.

**Not measured here.** Lake inflow through spilling sinks also changed: those nodes now
lose the ordinary transit share (at most about 2% per node). This script doesn't report lake
volumes. The two-density A/B of the whole change, including discharge, is #355 step 5.

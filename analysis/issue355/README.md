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

## Step 4: river visibility across node densities

Rendered the same seed (`355004`) at node densities 1 and 4, advanced each by 10 Myr, and
compared the elevation map and `/world/rivers` summaries. Broad river placement remained
readable at both densities; the finer mesh resolves more separate river networks. Replacing
the UI's 30-node visibility cutoff with a 1,875 km physical-length cutoff yields 11 and 14
eligible networks at densities 1 and 4. The cutoff is no longer directly tied to mesh spacing.
This is a visual spot check; the step 5 replay will measure model outcomes.

Counts and setup are recorded in [step4-visual.json](step4-visual.json).

## Step 5: calibrate river forcing and check glaciation

Replayed identical saves on `origin/main` at `93ab3bb` and this branch (starting from
`357eeac`, with the step-5 coefficients below), one process at a time. The synthetic young
worlds share seed `355004` and started at 0 Myr, at node densities 1 and 4. The mature
density-4 world is seed `997271774` at 87.4 Myr. Each pair ran 100 × 100 kyr steps (10 Myr).
The starting-save SHA-256 values are `f5fce4df…544ead52` (density 1),
`9ad8d1cf…c2555cd` (density 4), and `08ad9461…ed762b` (mature). Full replay JSON stayed
under `/private/tmp/issue355-step5/`; the compact endpoints are in
[step5-comparison.json](step5-comparison.json). Reproduce with
`bin/debug/measure_discharge_calibration.py --help`; the baseline backend can be selected
with `--backend`.

| 10 Myr endpoint | Density 1, main → branch | Density 4, main → branch | Mature density 4, main → branch |
|---|---:|---:|---:|
| Land area | 18.199% → 18.214% | 18.578% → 18.571% | 28.985% → 29.007% |
| Raw river erosion | 154k → 332k km³ | 199k → 224k km³ | 357k → 323k km³ |
| Channel depth p90 | 0.32 → 4.47 m | 4.22 → 9.21 m | 211 → 203 m |
| Channel width p90 | 849 → 1,710 m | 1,224 → 1,202 m | 5,000 → 5,000 m |
| Raw glacial erosion | 260 → 261 km³ | 1,478 → 4,551 km³ | 5.05 → 4.91 M km³ |
| Applied glacial erosion | 511 → 513 km³ | 2,866 → 13,258 km³ | 7.68 → 7.52 M km³ |
| Ice volume | 55.8k → 57.6k km³ | 280k → 266k km³ | 5.36 → 5.45 M km³ |
| Mean step time | 1.73 → 1.85 s | 4.77 → 4.70 s | 4.98 → 5.41 s |

`RIVER_EROSION_COEFFICIENT` moved from 100 to **190**. At density 4, leaving it at 100
made the young world's first-Myr raw river erosion 44% lower than main. A 160 trial left
the mature world's first-Myr raw erosion 24% lower; 190 left it 10% lower, and the mature
10 Myr budget remains 9.5% lower. The young density-4 10 Myr budget is 12.6% higher.
`WIDTH_GROWTH_COEFFICIENT` moved from 50 to **75**: the young density-4 p90 width is now
within 2% of main, while the mature p90 is at the existing 5 km cap in both runs. The
young p90 channel depth is still larger on the branch, since channel depth also reflects
sediment fill, glacier scour and the different terrain trajectory. The coarse-density
result is farther from its legacy baseline, which had strongly resolution-dependent
forcing; a controlled fixed-catchment test now checks that equal contributing area gives
equal outlet discharge, river erosion and width growth at both resolutions.

The river-depositing area moved from 47% to 23% of young density-4 land and from 61% to
31% of mature land. Sweeping `DEPOSITION_MIN_FLOW_M` from 0.05 down to 0.005 raises the
branch areas only to 27% and 39%; raising the unitless speed threshold from 2 to 8 at the
current flow cutoff raises them only to 25% and 32%. The main difference is runoff and
discharge replacing raw rainfall, including dry or very small streams, so matching the old
area by changing either threshold would erase the new physical-water gate. Their numerical
values, `RIVER_SPEED_COEFFICIENT` and the breach-surge coefficient remain unchanged. The
breach surge is an erosion-only lake forcing, separately tested not to enter discharge.

The mature replay is the glacial calibration reference: raw abrasion changed −2.9% and
applied glacial erosion −2.1% at 10 Myr, with the same coefficient (430). Young density-4
abrasion more than doubled but remains about 4,600 km³ against a mature 4.9 million km³;
changing the coefficient to fit that early transient would disturb the mature reference.
The ice ledger closes to 3 m³ or better over the mature runs. One branch 1 Myr diagnostic
from the mature save booked 150k km³ of snow accumulation, 37.6k km³ of frozen lake water,
35.0k km³ of frozen river ice, 396k km³ of meltwater and 11.3k km³ of calving. The snowfall
term remains in #353 model units: `0.02 m` per (mm/yr) per Myr against about `1,090 m`
physically, or 54,500× slower. That rate ratio does not describe the whole ice-volume
change, which includes those other inputs and losses. Glacier mass-balance conversion is
still #353's separate open decision; no glacier coefficient changed here.

Total erosion-ledger removed and deposited volumes stayed within 0.3% for each pair, and
land area within 0.03 percentage points. The mature run's mean step time rose 8.6% while
the young density-4 run was 1.5% faster; the replays do not isolate profiling noise from
the different terrain evolution. A/B worlds at different densities are not geometrically
identical, so these comparisons check model behaviour, not whole-world mesh convergence.

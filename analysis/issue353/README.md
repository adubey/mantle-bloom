# Issue #353: plastic ice sheets with per-body volume conservation

Replay of seed `997271774` from the 87.4 Myr save, 100 × 100 ky steps to 97.4 Myr (save
SHA-256 `08ad9461…ed762b`, as in [issue #317](../issue317/README.md)). Two runs:

- **main**: `origin/main` at `a17b6d2` (a CI-only commit ahead of this branch's base `ca3d25d`).
- **branch**: this branch at `42ff353`.

```bash
cd backend
.venv/bin/python ../bin/debug/measure_ice_sheet.py \
  ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 10 --checkpoint-myr 1 \
  --backend ../../mantle-bloom-origin-main/backend --out main.json
.venv/bin/python ../bin/debug/measure_ice_sheet.py \
  ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 10 --checkpoint-myr 1 --out branch.json
```

The compact per-checkpoint comparison is
[`replay-seed997271774.json`](replay-seed997271774.json); the full replay JSONs are not checked
in. The runs ran one after the other, each alone.

## At 97.4 Myr

| Measure | main | branch |
|---|---:|---:|
| Ice cap (visible ice, % of land) | 7.40% | 5.43% |
| Ice volume | 5.35 M km³ | 6.03 M km³ |
| Mean / max ice depth (visible ice) | 485 / 5,000 m | 735 / 2,141 m |
| Land (% of area) | 29.03% | 29.04% |
| Sea level | −397 m | −389 m |
| Glacial erosion, 10 Myr total (raw capacity) | 7.01 M km³ | 6.24 M km³ |
| Lake water frozen into ice, 10 Myr total | 15.3 M km³ | 0.10 M km³ |
| Step time (mean) | 4.41 s | 4.87 s |
| Ice solve share of step time | — | 9.1% |
| Worst ice-ledger residual (any step) | — | 3 m³ |
| Peak RSS | 2,635 MiB | 2,954 MiB |

The ice ledger closes to a few m³ on volumes of ~10¹⁵ m³ every step. Over the replay, 54–61%
of visibly glaciated area is warm-based, and the world holds 19–27 ice bodies.

## What the numbers mean

**Ice cover drops, volume holds.** The plastic profile gathers a body's volume into a dome
whose divide thickness is `sqrt(2 c L)`, so the same volume covers less ground than `main`'s
thin sheets: ice cap 5.4% against 7.4% of land, mean depth 735 m against 485 m. Volume ends
13% above `main`: `main` loses ice to its retired ad hoc sinks (the depth-squared basal melt
and the 5 km clip, which `main` sits at throughout), while the branch only melts, calves or
overflows what its ledger books. The branch's maximum depth comes off the 5 km cap after the
first few steps. Sea level ends 8 m higher on the branch even though its ice holds more
water. I haven't traced that; lake storage, which the lake-freeze change below moves, is the
likely place the difference sits. `main`'s ice-cap trajectory still rises
across the replay (6.2% → 7.4%) while the branch's settles (5.3–6.2%), so the gap is partly
`main` drifting.

**The first step reshapes the old save.** Nothing new is persisted, so the save's ice is
reshaped body by body on the first step, with its volume conserved. Ice cap falls from 6.23%
to 5.34% in the first Myr.

**Lakes no longer pump water into the ice.** Lakes balance in real units since #297, but snow
enters the ice in model units (~50,000x slower). On `main`, ~127 lakes sit on warm ground
under ice; each step they refill with runoff and freeze 1:1 into `glacier_depth`, 15.3 M km³ over
10 Myr against ~2 M km³ of snowfall. `main` throws that ice away again through the sinks this
issue retires. With volume conserved, the first branch replay grew its ice from 5.7 to
25.7 M km³ in 10 Myr and dropped sea level 36 m. Lakes now freeze only where the annual mean is
below freezing; the branch's remaining lake ice (0.1 M km³) is from genuinely cold ground.

**Glacial erosion is calibrated on raw capacity.** The glacier term's detachment capacity
before caps (`main`: coefficient × slope × ice factor; branch: coefficient × sliding speed)
is the only measure both trees compute the same way. `GLACIER_ABRASION_COEFFICIENT = 430` was
set from an earlier 10 Myr branch run to match `main`'s capacity. The final run comes in 11%
below it because the ice and terrain trajectories differ. The branch's applied erosion
(10.3 M km³) is above its capacity because sliding ice now abrades warm, sediment-covered beds
whose mobile cover comes off faster than rock (`applied_scale > 1`). `main`'s applied glacial
erosion isn't exposed, so it can't be compared directly.

**Cost.** The ice solve is ~9% of step time, inside the issue's 10% target; steps are ~10%
slower overall (4.87 s against 4.41 s).

**Noise.** Run-to-run noise on continental area is about ±2% at 43 Myr (#315); land area
matches `main` within it here. These runs are 10 Myr on one save. The ice-cap and volume
gaps are larger than that noise and consistent across checkpoints.

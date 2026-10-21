# Issue #503: rivers that end on dry land

The investigation is in [issue #503](https://github.com/adubey/mantle-bloom/issues/503): see
its Exploration summary and its preregistered probe comments. This page records the replay
comparison for the fix on branch `claude/issue-503-rivers-to-nowhere`.

## Preregistration (committed before either run)

**Save:** seed `649658380`, 6.4 Myr (`<local-save>`, SHA-256 `ebbfb5da…7b3`).

**Runs:** 100 steps of 100 kyr each, from 6.4 to 16.4 Myr, with checkpoints every 2 Myr. The
two arms run one after the other, never in parallel.

- **main:** a detached `origin/main` worktree at `f9d3d2f`.
- **fix:** this branch at `84dd278`.

```bash
cd backend
.venv/bin/python ../bin/debug/measure_river_endings.py <local-save> --myr 10 --checkpoint-myr 2 \
    --backend <baseline-worktree>/backend --out <scratch>/main.json
.venv/bin/python ../bin/debug/measure_river_endings.py <local-save> --myr 10 --checkpoint-myr 2 \
    --out <scratch>/fix.json
```

[`measure_river_endings.py`](../../bin/debug/measure_river_endings.py) records each step's
measures. A checkpoint reports the mean over the steps since the previous one, and the
headline is the mean over all 100 steps.

### Primary metric

**Water reaching no outlet** (m³/s): inflow into warm (≥ 0 °C), ice-free land sinks that hold
no visible lake, plus inflow into `flow_target` cycles. These are rivers that reach neither the
sea nor a lake.

**Prediction:** the fix's 100-step mean is **at most 50% of main's**. I expect a reduction of
about 85–95%. The issue's probes took it from about 48,000 to about 19,000 m³/s with sibling
merges only, and a 5-step smoke run of this branch measured 2,400–6,500 m³/s.

**Decision:** if the fix's mean is at most 50% of main's, the fix works. Otherwise it doesn't.

### Secondary metrics (direction and rough size)

| Metric | Prediction |
|---|---|
| Flow cycles per step | fix about 0 (≤ 0.1 mean); main 0–3 |
| Discharge into the ocean | fix ≥ main + 5% |
| Drawn river lines ending off the water (elevation map) | fix ≤ 50% of main |
| Length-visible networks ending at `other` | fix < main |
| Warm single-node ice plugs | fix < main |
| Frozen river ice | fix ≤ 70% of main |

### No-harm guards

| Guard | Bound |
|---|---|
| Glacier volume at 16.4 Myr | within ±10% of main |
| Step time | fix ≤ 1.5× main |
| Peak RSS | fix ≤ 1.2× main |
| Visible lake nodes | reported, with no prediction: lakes may grow where they now fill within the step |

### Variability

Each arm is deterministic, so I run it once. Step-to-step variability is the standard
deviation over each checkpoint window, which the script reports. A gap smaller than either
arm's standard deviation is called inconclusive. A gap near it gets the trajectory across
checkpoints in place of the single 100-step mean.

## Results

Both arms ran one after the other on the same machine, with nothing else running, from 6.4 to
16.4 Myr. The compact numbers are in [`compare-seed649658380.json`](compare-seed649658380.json).

| Mean over 100 steps | main `f9d3d2f` | fix `b6a2a46` | Prediction | Outcome |
|---|---:|---:|---|---|
| **Water reaching no outlet** | 53,819 m³/s | **4,939 m³/s (−91%)** | ≤ 50% of main | **met** |
| … into invisible sinks | 32,192 m³/s | 4,939 m³/s | | |
| … into flow cycles | 21,627 m³/s | 0 m³/s | | |
| Flow cycles per step | 1.5 | 0.0 | ≤ 0.1 | met |
| Discharge into the ocean | 629,638 m³/s | 714,072 m³/s (+13%) | ≥ +5% | met |
| Drawn river lines ending off the water | 20.6 | 4.7 (23%) | ≤ 50% | met |
| Length-visible networks ending at `other` | 5.3 | 2.6 | lower | met |
| Warm single-node ice plugs | 2.8 | 1.7 | lower | met |
| Frozen river ice per step | 4,005 km³ | 3,921 km³ (98%) | ≤ 70% | **not met** |
| Visible lake nodes | 1,806 | 1,945 (+8%) | reported | |
| **Glacier volume at 16.4 Myr** | 7.92 M km³ | **10.59 M km³ (+34%)** | within ±10% | **guard failed** |
| Step time | 4.7 s | 5.0 s (1.06×) | ≤ 1.5× | met |
| Peak RSS | 1,938 MB | 1,999 MB (1.03×) | ≤ 1.2× | met |

**The gap is well beyond the noise.** Water reaching no outlet varies by about 20,000–28,000
m³/s from step to step on main and about 1,000–4,800 m³/s on the fix. The gap of about 48,900
m³/s holds in every 2 Myr window: 64,609 → 4,616, 50,007 → 5,354, 40,643 → 4,388, 45,594 →
4,452 and 68,241 → 5,888 m³/s.

**The primary metric is met: rivers that went nowhere now reach the sea or a lake.** Flow
cycles are gone. The water that's left in invisible sinks (about 4,900 m³/s) is mostly small
endorheic basins whose lakes stay under 1 m, as a dry climate should give.

**The glacier guard failed, and the cause is a pre-existing coupling.** Glacier volume drifts
apart steadily:

| Elapsed (Myr) | 6.4 | 8.4 | 10.4 | 12.4 | 14.4 | 16.4 |
|---|---:|---:|---:|---:|---:|---:|
| main (M km³) | 3.56 | 7.83 | 7.77 | 7.95 | 8.22 | 7.92 |
| fix (M km³) | 3.56 | 8.02 | 8.49 | 8.84 | 9.51 | 10.59 |

A second pair of 3 Myr runs recorded the whole ice ledger. These are the means per step, in
km³:

| Term | main | fix |
|---|---:|---:|
| Lake water frozen into the ice | 147,479 | 177,580 (+20%) |
| Snow accumulation | 4,274 | 4,459 |
| Frozen river ice | 3,962 | 2,974 (−25%) |
| Surface melt | 10,205 | 10,790 |
| Calved | 3,202 | 5,819 |

- **Lake water freezing 1:1 into glacier ice is about 35× snowfall on both arms.** Lakes
  balance in physical units since #297, while glacier mass balance is still in model units
  (about 50,000× slower, decision B in #353, studied in #368). That conversion already makes
  main's ice double in its first 2 Myr.
- **The fix adds to it.** By filling and spilling basins within the step, it delivers more
  water to lakes on ground below freezing, and those lakes freeze into the ice the next step.
  The extra 30,000 km³ per step accounts for the drift. Warm ice plugs (H1) are not the source:
  frozen river ice *fell* 25% over the first 3 Myr.
- **Frozen river ice barely moved over 100 steps.** Most of it is rivers freezing on ground
  genuinely below 0 °C, which the fix leaves alone.

## Decision needed

The fix meets its primary metric. Its glacier-volume guard fails because of the lake-to-ice unit
mismatch tracked in #368, which the fix amplifies but doesn't introduce. The decision on how
to proceed belongs to the project owner (see the PR).

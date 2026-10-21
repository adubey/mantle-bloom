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

Pending.

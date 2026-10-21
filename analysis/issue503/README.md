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

## Amendment: frozen lakes (committed before the second run)

The first run's glacier guard failed (see *Results: first run* below), and the owner chose to fix lake
freezing in this PR. A freezing lake now turns into ice at `hydrology.WATER_TO_ICE`, the rate
frozen river water already used, instead of 1:1. That deliberately removes main's largest ice
source (about 147,000 km³ per step, around 35× snowfall), so a "within ±10% of main" guard no
longer applies.

**Runs:** both arms again, one after the other, with the same command and the script updated to
record sea level: main `f9d3d2f` and fix2 (this branch at the commit that adds this
amendment).

| Metric | Prediction |
|---|---|
| **Primary:** water reaching no outlet, mean over 100 steps | fix2 ≤ 50% of main, as before (expected about 5,000 m³/s) |
| Lake water frozen into the ice, per step | fix2 ≤ 1/1,000 of main |
| Glacier volume at 16.4 Myr | fix2 below main, and between 2.0 and 4.5 M km³ (the save starts at 3.56 M km³). Snow, river ice and melt alone come to about −5,000 km³ per step, so the ice should shrink slowly, neither growing beyond start + 25% nor collapsing below start − 45%. |
| Sea level at 16.4 Myr | fix2 above main: less water is locked in ice. Reported; no size predicted. |
| Glacier ledger | must close every step |

The secondary metrics and the step-time and RSS guards are unchanged from the original
preregistration.

## Amendment 2: review fixes (committed before the third run)

The PR #504 review led to these changes:

- **Fill test:** uses `lakes.lake_evaporation_m3_per_yr`, including the shallow-margin boost.
- **Frozen sinks:** no longer count toward filling a basin.
- **Pit chains:** basins are swept from the highest floor down, passing released water on, and hitting the round cap now logs a warning.
- **Held unions:** pool their water at the lowest sink.
- **Breach gate:** forwards surplus from each pit's own spill level, upstream first.
- **Drawing:** fills in only river nodes below the width floor.

**Run:** fix3 at this amendment's commit. The comparison is the second run's main arm, which is deterministic, with the script unchanged since then. Same command and save.

| Metric | Prediction for fix3 |
|---|---|
| **Primary:** water reaching no outlet, mean over 100 steps | ≤ 50% of main's 53,819 m³/s; expected 4,000–8,000 m³/s, near fix2's 6,104 |
| Flow cycles per step | ≤ 0.1 |
| Discharge into the ocean | ≥ main + 5% |
| Drawn river lines ending off the water | ≤ 50% of main |
| Glacier volume at 16.4 Myr | within 2.0–4.5 M km³ |
| Glacier ledger | closes |
| Round-cap warnings in the log | none expected; any are reported |
| Step time | ≤ 1.5× main |
| Peak RSS | ≤ 1.2× main |

## Results: first run (fix `b6a2a46`, lakes still freezing 1:1)

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

The fix met its primary metric, but its glacier-volume guard failed because of the lake-to-ice
unit mismatch tracked in #368, which the fix amplified but didn't introduce. The owner chose
to fix lake freezing in this PR (see the amendment above).

## Results: second run (fix2 `a1d5783`, frozen lakes at `WATER_TO_ICE`)

Main reran identically to the first run, so the arms are deterministic. Both ran one after the
other with nothing else running.

| Mean over 100 steps, unless noted | main `f9d3d2f` | fix2 `a1d5783` | Prediction | Outcome |
|---|---:|---:|---|---|
| **Water reaching no outlet** | 53,819 m³/s | **6,104 m³/s (−89%)** | ≤ 50% of main | **met** |
| Lake water frozen into the ice, per step | 67,505 km³ | 1.5 km³ | ≤ 1/1,000 of main | met |
| Glacier volume at 16.4 Myr | 7.92 M km³ | 3.82 M km³ | below main, within 2.0–4.5 M km³ | met |
| Sea level at 16.4 Myr | −183.2 m | −172.8 m | above main | met |
| Glacier-ledger residual | 0 | 0 | closes | met |
| Flow cycles per step | 1.5 | 0.0 | ≤ 0.1 | met |
| Discharge into the ocean | 629,638 m³/s | 721,200 m³/s (+15%) | ≥ +5% | met |
| Drawn river lines ending off the water | 20.6 | 4.3 (21%) | ≤ 50% | met |
| Length-visible networks ending at `other` | 5.3 | 2.3 | lower | met |
| Warm single-node ice plugs | 2.8 | 1.5 | lower | met |
| Frozen river ice per step | 4,005 km³ | 3,909 km³ (98%) | ≤ 70% | **not met** |
| Snow accumulation per step | 4,972 km³ | 2,212 km³ | | |
| Visible lake nodes | 1,806 | 1,894 | reported | |
| Step time / peak RSS | 4.6 s / 1,888 MB | 4.8 s / 1,979 MB | ≤ 1.5× / ≤ 1.2× | met |

**Water reaching no outlet, by 2 Myr window:**

| Window | main | fix2 |
|---|---:|---:|
| 1 | 64,609 m³/s | 4,716 m³/s |
| 2 | 50,007 m³/s | 4,639 m³/s |
| 3 | 40,643 m³/s | 4,971 m³/s |
| 4 | 45,594 m³/s | 4,937 m³/s |
| 5 | 68,241 m³/s | 11,259 m³/s |

The last window is fix2's highest, still 84% below main's.

**Glacier volume and sea level:**

| Elapsed (Myr) | 6.4 | 8.4 | 10.4 | 12.4 | 14.4 | 16.4 |
|---|---:|---:|---:|---:|---:|---:|
| main glacier (M km³) | 3.56 | 7.83 | 7.77 | 7.95 | 8.22 | 7.92 |
| fix2 glacier (M km³) | 3.56 | 3.61 | 3.65 | 3.69 | 3.76 | 3.82 |
| main sea level (m) | −100.5 | −131.7 | −144.8 | −160.1 | −175.9 | −183.2 |
| fix2 sea level (m) | −100.5 | −120.7 | −134.2 | −148.3 | −160.6 | −172.8 |

### What it shows

- **Rivers reach the sea or a lake.** Water reaching no outlet falls by 89%, and flow cycles are
  gone. Discharge into the ocean rises 15%, and drawn river lines ending off the water fall
  from 20.6 to 4.3 per step.
- **The ice no longer doubles from frozen lakes.** On main, lake water freezing 1:1 doubles the
  ice in the first 2 Myr. Under fix2, the ice grows 7% over 10 Myr from snow and frozen river
  water against a smaller melt.
- **The preregistered note was wrong about the ice's direction.** I described a slow shrink,
  which assumed main's melt of about 11,500 km³ per step. Over the smaller ice area, melt is
  about 3,200 km³ per step, and frozen river ice (about 3,900) plus snow (about 2,200) exceed
  it. The prediction's bounds still held.
- **Frozen river ice barely changes.** It is rivers freezing on ground genuinely below 0 °C,
  which neither fix touches. Whether a large river should freeze solid at a below-zero annual
  mean, or flow in the warm season, is left to #368.
- **Sea level falls on both arms** (−72 m for fix2, −83 m for main over 10 Myr), so the trend
  is not from this change; less water locked in ice softens it.
- **Cost:** step time +4% and peak RSS +5%.

## Results: third run (fix3 `f692507`, review fixes)

This compares against the second run's main arm, which is deterministic. The script is
unchanged. The fix3 run ran alone.

| Mean over 100 steps, unless noted | main `f9d3d2f` | fix2 `a1d5783` | fix3 `f692507` | Prediction | Outcome |
|---|---:|---:|---:|---|---|
| **Water reaching no outlet** | 53,819 m³/s | 6,104 m³/s | **5,552 m³/s (−90%)** | ≤ 50% of main | **met** |
| Flow cycles per step | 1.5 | 0.0 | 0.0 | ≤ 0.1 | met |
| Discharge into the ocean | 629,638 m³/s | 721,200 m³/s | 722,559 m³/s (+15%) | ≥ +5% | met |
| Drawn river lines ending off the water | 20.6 | 4.3 | 4.4 (21%) | ≤ 50% | met |
| Length-visible networks ending at `other` | 5.3 | 2.3 | 2.5 | | |
| Warm single-node ice plugs | 2.8 | 1.5 | 1.5 | | |
| Frozen river ice per step | 4,005 km³ | 3,909 km³ | 3,685 km³ | | |
| Visible lake nodes | 1,806 | 1,894 | 1,903 | | |
| Glacier volume at 16.4 Myr | 7.92 M km³ | 3.82 M km³ | 3.79 M km³ | 2.0–4.5 M km³ | met |
| Sea level at 16.4 Myr | −183.2 m | −172.8 m | −173.0 m | | |
| Glacier-ledger residual | 0 | 0 | 0 | closes | met |
| Round-cap warnings | — | — | 0 | none | met |
| Step time | 4.61 s | 4.81 s | 5.41 s (1.17×) | ≤ 1.5× | met |
| Peak RSS | 1,888 MB | 1,979 MB | 1,494 MB | ≤ 1.2× | met |

**Water reaching no outlet, by 2 Myr window:**

| Window | main | fix3 |
|---|---:|---:|
| 1 | 64,609 m³/s | 4,438 m³/s |
| 2 | 50,007 m³/s | 5,013 m³/s |
| 3 | 40,643 m³/s | 4,975 m³/s |
| 4 | 45,594 m³/s | 5,527 m³/s |
| 5 | 68,241 m³/s | 7,805 m³/s |

The last window's rise is smaller than fix2's (11,259 m³/s).

**Every preregistered prediction is met.**

**Step time rose more than before:** from +4% (fix2) to +17%. The likely cause is the review
fixes' per-pit Python walks: `breaching._spill_exits` runs on every routable pit, and the
resolver sweeps basins with a path walk per opened basin. Step time stays within the 1.5×
guard, and vectorising those walks is an option if the cost matters. Peak RSS was lower in
this run than on main. RSS is noisy across runs, and this one ran alone.

### The third run is superseded (note committed before the fourth run)

After the third run, the unit tests caught a bug in the review fix to the breach gate.

- **What went wrong:** a pit that hands its water to a lake that drowns it was judged on its own catchment. So one cell of a three-cell basin could breach on the whole basin's water and only its own evaporation.
- **Effect:** `test_compute_hydrology_breaches_a_shallow_pit_instead_of_ponding_it`'s semi-arid pit was breached.
- **Fix:** such a pit is now judged as the lake it joins, and pits join at equal levels too.

The third run's numbers came from the buggy gate, so they're superseded. A fourth run (fix4, at the commit adding this note) is judged against Amendment 2's predictions, which are unchanged.

## Results: fourth run (fix4 `307a90f`, the PR's final model)

This compares against the second run's main arm, which is deterministic, and is judged
against Amendment 2's predictions. The fix4 run ran alone.

| Mean over 100 steps, unless noted | main `f9d3d2f` | fix4 `307a90f` | Prediction | Outcome |
|---|---:|---:|---|---|
| **Water reaching no outlet** | 53,819 m³/s | **6,008 m³/s (−89%)** | ≤ 50% of main | **met** |
| Flow cycles per step | 1.5 | 0.0 | ≤ 0.1 | met |
| Discharge into the ocean | 629,638 m³/s | 725,598 m³/s (+15%) | ≥ +5% | met |
| Drawn river lines ending off the water | 20.6 | 4.7 (23%) | ≤ 50% | met |
| Length-visible networks ending at `other` | 5.3 | 2.5 | | |
| Warm single-node ice plugs | 2.8 | 1.5 | | |
| Frozen river ice per step | 4,005 km³ | 3,887 km³ | | |
| Visible lake nodes | 1,806 | 1,940 | | |
| Glacier volume at 16.4 Myr | 7.92 M km³ | 3.80 M km³ | 2.0–4.5 M km³ | met |
| Sea level at 16.4 Myr | −183.2 m | −174.4 m | | |
| Glacier-ledger residual | 0 | 0 | closes | met |
| Round-cap warnings | — | 0 | none | met |
| Step time | 4.61 s | 5.49 s (1.19×) | ≤ 1.5× | met |
| Peak RSS | 1,888 MB | 1,695 MB | ≤ 1.2× | met |

**Water reaching no outlet, by 2 Myr window:**

| Window | main | fix4 |
|---|---:|---:|
| 1 | 64,609 m³/s | 5,245 m³/s |
| 2 | 50,007 m³/s | 5,656 m³/s |
| 3 | 40,643 m³/s | 5,252 m³/s |
| 4 | 45,594 m³/s | 6,456 m³/s |
| 5 | 68,241 m³/s | 7,429 m³/s |

**Every preregistered prediction is met.**

**Changes from fix2:**

- **More water reaches the sea.** Ocean discharge is 725.6k m³/s against fix2's 721.2k.
- **The late rise in water reaching no outlet is smaller:** the last window is 7,429 m³/s, against fix2's 11,259.
- **The overall mean is about the same:** 6,008 m³/s against 6,104.
- **Steps are slower:** about 1.19× main, against fix2's 1.04×. The likely cause is the breach gate's per-pit Python walks (`_spill_exits` and the topological forwarding) and the resolver's sweep.

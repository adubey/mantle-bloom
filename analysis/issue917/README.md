# Issue 917: RFC 0015 D6, sea level: debit ice above flotation

Implementation of [RFC 0015 D6](../../rfc/0015-frozen-water-in-physical-units.md#d6-the-ocean)
(sea level half) for [#917](https://github.com/adubey/mantle-bloom/issues/917). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544).

**Principle:** #1, conservation. Sea level comes from a conserved water budget, so every cubic
metre booked as ice must be water that has left the ocean.

## Pre-registered comparison

- **Saves:** the [#548](../issue548/README.md) diagnosis saves (converted to `.bloom` in B2):
  - `649658380` at 6.4 Myr (`mantle-bloom-seed649658380-6400000y.bloom`, SHA-256 `06167cb4…31131f`, deciding, ice ages off)
  - `997271774` at 87.4 Myr (`mantle-bloom-seed997271774-87400000y.bloom`, SHA-256 `07971ee2…328b8`, deciding, ice ages off)
  - `374464662` at 132.6 Myr (`mantle-bloom-seed374464662-132600000y.bloom`, SHA-256 `a9eed5da…2801a`, report only, 8 Myr ice-age period)
- **Protocol:** 10 steps at 0.25 Myr (primary) and 1 Myr (secondary) with the default `trend-1k`
  glacier march, comparing `origin/main` (`011f288d`) against the branch (`011f288d` + #917 changes).
- **Primary metric & threshold:**
  - On the deciding saves, the change in sea level from `origin/main` is expected to be positive
    and small (within G4's measured over-debit of ≤ 0.08 m, plus coupled-run drift).
  - After the change, the replay's own `sea_level_model_m` matches `measure_flotation_debit.py`'s
    generalised arm `h_generalised_m` (to within eustasy's solver tolerance, ~1 m).
  - Land and glacier water ledger residuals remain closed / unchanged relative to `origin/main`.
  - Ice volume is reported alongside sea level per step and at the end.

## Runs

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_flotation_debit.py \
    <local-save> --steps 10 --myr {0.25|1} [--ice-age-period-years 8e6] --out <scratch-output>
```

- **Captures:** full per-step JSONs and `SHA256SUMS` are in `~/mantle-bloom-captures/issue-917/`.
  Compact per-step comparison is in [`ab.json`](ab.json).
- **Resources:** peak RSS was 1.09–1.14 GB (`649658380`), 1.57–1.63 GB (`997271774`), and
  1.17–1.23 GB (`374464662`). All runs ran one at a time alongside two other background
  simulation jobs on the 16 GB machine (with owner override).

## Results

### End-of-run summary

| Save | Step (Myr) | Role | Baseline SL (m) | Branch SL (m) | Δ SL (m) | Branch max \|model − G\| (m) | Baseline ice (km³) | Branch ice (km³) | Δ ice |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| `649658380` | 0.25 | deciding | −111.310 | −111.316 | −0.005 | **0.0000** | 2,443,738 | 2,444,873 | +0.05% |
| `649658380` | 1.0 | deciding | −104.587 | −104.459 | +0.128 | **0.0000** | 2,179,809 | 2,144,466 | −1.62% |
| `997271774` | 0.25 | deciding | −473.430 | −473.016 | +0.415 | **0.0000** | 14,527,817 | 14,390,171 | −0.95% |
| `997271774` | 1.0 | deciding | −536.707 | −536.048 | +0.658 | **0.0000** | 12,509,259 | 12,393,518 | −0.93% |
| `374464662` (ice ages) | 0.25 | report only | −338.968 | −338.334 | +0.634 | **0.0000** | 4,931,289 | 4,715,742 | −4.37% |

### Generalised-arm cross-check & ledgers

- **Exact agreement with `h_generalised_m`:** Across all 50 steps of all 5 branch runs,
  `max |sea_level_model_m − h_generalised_m|` is **0.0000 m** (< 0.0001 m). Evaluating `D` at
  the candidate sea level inside `_solve_sea_level_connected` reproduces G4's generalised arm.
- **Step-0 direct debit shift:** On step 0 (before coupled ice-march drift accumulates), the
  sea-level shifts (`branch − baseline`) are `+0.000 m` (`649658380` @ 0.25 Myr), `+0.001 m`
  (`649658380` @ 1 Myr), `+0.044 m` (`997271774` @ 0.25 Myr), `+0.045 m` (`997271774` @ 1 Myr),
  and `+2.002 m` (`374464662` @ 0.25 Myr), matching G4's measured over-debit (`≤ 0.08 m` on
  deciding saves, `~2.4 m` on the ice-age save). Later steps show small coupled ice-volume
  differences (−1.6% to +0.05% on deciding saves) from the march responding to the updated
  waterline.
- **Ledger residuals:** Both land and glacier water ledgers remain closed at floating-point
  precision:

| Save | Step (Myr) | Max \|land residual\| baseline (m³) | Max \|land residual\| branch (m³) | Max \|glacier residual\| baseline (m³) | Max \|glacier residual\| branch (m³) |
|---|---:|---:|---:|---:|---:|
| `649658380` | 0.25 | 9.70 × 10³ | 4.24 × 10³ | 0 | 8 |
| `649658380` | 1.0 | 1.58 × 10⁴ | 1.59 × 10⁴ | 32 | 32 |
| `997271774` | 0.25 | 5.72 × 10³ | 5.34 × 10³ | 256 | 256 |
| `997271774` | 1.0 | 2.30 × 10⁴ | 2.55 × 10⁴ | 0 | 1.02 × 10³ |
| `374464662` (ice ages) | 0.25 | 6.09 × 10³ | 6.04 × 10³ | 64 | 256 |

### Per-step sea level & ice volume (primary 0.25 Myr deciding saves)

#### `649658380` @ 0.25 Myr

| Step | Age (Myr) | Baseline SL (m) | Branch SL (m) | Δ SL (m) | Baseline ice (km³) | Branch ice (km³) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 6.65 | −101.127 | −101.127 | +0.000 | 3,104,549 | 3,104,549 |
| 1 | 6.90 | −100.732 | −100.732 | −0.000 | 2,651,138 | 2,651,138 |
| 2 | 7.15 | −100.581 | −100.582 | −0.001 | 2,661,430 | 2,661,399 |
| 3 | 7.40 | −100.734 | −100.729 | +0.005 | 2,531,282 | 2,530,921 |
| 4 | 7.65 | −102.995 | −102.992 | +0.003 | 2,572,732 | 2,572,643 |
| 5 | 7.90 | −105.003 | −105.006 | −0.003 | 2,444,514 | 2,444,724 |
| 6 | 8.15 | −106.412 | −106.361 | +0.051 | 2,452,036 | 2,451,554 |
| 7 | 8.40 | −108.597 | −108.567 | +0.030 | 2,446,773 | 2,446,404 |
| 8 | 8.65 | −108.963 | −108.958 | +0.005 | 2,523,582 | 2,521,987 |
| 9 | 8.90 | −111.310 | −111.316 | −0.005 | 2,443,738 | 2,444,873 |

#### `997271774` @ 0.25 Myr

| Step | Age (Myr) | Baseline SL (m) | Branch SL (m) | Δ SL (m) | Baseline ice (km³) | Branch ice (km³) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 87.65 | −438.443 | −438.400 | +0.044 | 15,868,020 | 15,868,021 |
| 1 | 87.90 | −440.579 | −440.515 | +0.064 | 15,716,450 | 15,715,502 |
| 2 | 88.15 | −442.152 | −442.143 | +0.008 | 14,290,302 | 14,297,931 |
| 3 | 88.40 | −447.352 | −447.359 | −0.007 | 14,415,191 | 14,425,269 |
| 4 | 88.65 | −452.148 | −452.110 | +0.038 | 14,481,730 | 14,493,966 |
| 5 | 88.90 | −456.902 | −456.561 | +0.341 | 14,599,365 | 14,475,747 |
| 6 | 89.15 | −460.464 | −460.522 | −0.058 | 14,543,481 | 14,563,262 |
| 7 | 89.40 | −463.968 | −463.411 | +0.558 | 14,513,890 | 14,305,726 |
| 8 | 89.65 | −468.432 | −468.474 | −0.041 | 14,357,377 | 14,376,112 |
| 9 | 89.90 | −473.430 | −473.016 | +0.415 | 14,527,817 | 14,390,171 |

# Issue #1046 — Track C Baseline & Step 2 Lake Kernels in Rust

## Summary

1. **Post-#1037/#1044 baseline (`origin/main` @ `2c9ebffa`):** Measured young (`seed 0`, 124,046 nodes, 17 plates) and mature (`mantle-bloom-seed47235302-67900000y.bloom`, 67.9 Myr, 125,497 nodes, 28 plates, SHA-256 `eff1f6d7abdfe46c7af5b45ab23f136f3fac591f4d95f3457017c1e7d933216b`) stepping with `bin/debug/profile_phases.py`.
2. **Sequencing confirmation:** On the mature world, non-march step time is **10.68 s / 20.43 s (52.3%)**. Within non-march, `Updating faults` is 2.36 s/step (11.6%), `hydrology.py` routing (`Ordering the flow network` + `Routing flow downstream`) is ~1.69 s/step, `lakes.py` + `hydrology.lake_components` is ~0.79 s/step (2.366 s across 3 steps in cProfile), and `climate.py` is 0.44 s/step. Keeping the #1046 order (`lakes` → `hydrology` routing → `climate` → shared 3D spatial index → `faults`) remains sound because `faults.py` depends on the shared spatial index (step 5) while `lakes` and `hydrology` are self-contained array/graph ports.
3. **Step 2 (`lakes.rs` + `_lakes_rust.py`):** Ported `lake_components_rs`, `catchment_roots_rs`, `lake_hierarchy_events_rs`, and `_integrate_body_volume`'s multi-member Newton/bisection loop (`integrate_body_volume_rs`, triggered because `_integrate_body_volume` was **21.8%** of mature lake time in `cProfile`, above the 15% threshold). Kernels take `neighbor_idx` directly rather than `MarchGraph` because lake callers do not hold a `MarchGraph`. `_settle_prognostic_slice` (`0.017 s` own time over 3 steps) and Python `Lake`/accounting orchestration (`~0.21 s/step`) remain in Python and are deferred.

## Commands

```bash
# Young world (mean of 3 steps after 1 warmup, diagnostic wind)
/usr/bin/time -l backend/.venv/bin/python bin/debug/profile_phases.py \
  --seed 0 --scenarios step --steps 3 --warmup 1 --wind-model diagnostic \
  --json <scratch-output>/young.json

# Mature world (67.9 Myr save, mean of 3 steps after 1 warmup)
/usr/bin/time -l backend/.venv/bin/python bin/debug/profile_phases.py \
  --load <local-save> --scenarios step --steps 3 --warmup 1 --wind-model diagnostic \
  --json <scratch-output>/mature.json --cprofile <scratch-output>/mature.prof

# Mature-save 1 Myr replay A/B (67.9 -> 68.9 Myr, 10 steps)
PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
  <local-save> --myr 1.0 --checkpoint-myr 1.0 --out <scratch-output>/replay.json
```

> **Note:** All timing and replay runs ran with an explicit compute-budget override on a shared 16 GB Apple Silicon host while other background simulations were active (`0` swaps in `/usr/bin/time -l`).

## Before / After Phase Profiles

### Young World (`seed 0`, 124,046 nodes, 17 plates, `wind_model=diagnostic`)

| Metric | `origin/main` (`2c9ebffa`) | Branch (`issue-1046`) | Delta / Speedup |
|---|---:|---:|---:|
| Step wall time (mean of 3) | **17.01 s** (`16.63, 17.34, 17.06`) | **15.97 s** (`16.31, 15.65, 15.95`) | **1.07x** (-1.04 s) |
| Glacier march total (`Solving march ice` + `Re-solving march sea level`) | 9.32 s | 8.95 s | 1.04x |
| Non-march step time | **7.69 s** (45.2%) | **7.02 s** (44.0%) | **1.10x** (-0.67 s) |
| — `Resolving lakes` | **0.42 s** (2.5%) | **0.31 s** (1.9%) | **1.35x** (-0.11 s) |
| — `Routing flow downstream` (all) | 0.85 s | 0.80 s | — |
| — `Ordering the flow network` (all) | 0.74 s | 0.65 s | — |
| — `Resolving climate` | 0.49 s | 0.43 s | — |
| — `Updating faults` | 1.00 s | 1.02 s | — |
| Peak RSS (`/usr/bin/time -l`) | 990.3 MB | 854.2 MB | -13.7% |

### Mature World (`mantle-bloom-seed47235302-67900000y.bloom`, 67.9 Myr, 125,497 nodes, 28 plates)

| Metric | `origin/main` (`2c9ebffa`) | Branch (`issue-1046`) | Delta / Speedup |
|---|---:|---:|---:|
| Step wall time (mean of 3) | **20.43 s** (`21.73, 20.27, 19.29`) | **19.67 s** (`19.76, 19.27, 19.99`) | **1.04x** (-0.76 s) |
| Glacier march total (`Solving march ice` + `Re-solving march sea level`) | 9.75 s | 9.66 s | 1.01x |
| Non-march step time | **10.68 s** (52.3%) | **10.01 s** (50.9%) | **1.07x** (-0.67 s) |
| — `Resolving lakes` (tracer phase) | **0.41 s** (2.0%) | **0.22 s** (1.1%) | **1.86x** (-0.19 s) |
| — `Ordering the flow network` (all) | 0.86 s | 0.83 s | — |
| — `Routing flow downstream` (all) | 0.83 s | 0.80 s | — |
| — `Resolving climate` | 0.44 s | 0.43 s | — |
| — `Updating faults` | 2.36 s (11.6%) | 2.21 s (11.3%) | — |
| Peak RSS (`/usr/bin/time -l`) | 927.7 MB | 849.2 MB | -8.5% |

## Mature cProfile Breakdown (3 Measured Steps)

| Function | Calls | `origin/main` cumtime (s) | Branch cumtime (s) | Speedup |
|---|---:|---:|---:|---:|
| `hydrology.lake_components` | 12 | 0.7748 | 0.0233 | **33.3x** |
| `lakes.build_lake_hierarchy` (incl. `_catchment_roots`) | 3 | 0.3678 | 0.0748 | **4.9x** |
| `lakes._integrate_body_volume` | 3,978 | 0.5164 | 0.0275 | **18.8x** |
| **Combined ported lake kernels** | — | **1.6590** | **0.1256** | **13.2x** |
| `lakes.solve_lakes_prognostic` (incl. `_integrate_body_volume`) | 3 | 1.2235 | 0.6605 | 1.85x |
| `lakes._settle_prognostic_slice` (own time) | 12 | 0.0165 | 0.0135 | — |

- `_integrate_body_volume` represented `0.5164 s / (1.2235 + 0.7748 + 0.3678) = 21.8%` of mature lake time on `origin/main`, meeting the ≥15% criterion for porting its multi-member Newton/bisection loop.
- Remaining time inside `solve_lakes_prognostic` (`~0.63 s` over 3 steps, `~0.21 s/step`) is Python `Lake` object / accounting / slice orchestration (`_advance_prognostic_body`, `_active_bodies`, `_ordered_spill_bodies`), deferred from this PR.

## Mature-Save 1 Myr Replay A/B (`67.9 → 68.9 Myr`, 10 Steps)

Because `integrate_body_volume_rs` sums float64 terms sequentially rather than via `np.dot` (parity within `rtol=1e-12`), a 10-step replay was run on `mantle-bloom-seed47235302-67900000y.bloom`:

| Metric at 68.9 Myr | `origin/main` (`2c9ebffa`) | Branch (`issue-1046`) | Delta |
|---|---:|---:|---:|
| `continental_area_km2` | 354,249,271.18 | 354,249,271.18 | 0.00 |
| `hm_cap_area_km2` | 542,048.26 | 542,048.26 | 0.00 |
| `hm_volume_all_km3` | 40,616,446,488.1490 | 40,616,446,488.1495 | +0.0004 km³ (+1.1e-14 rel) |
| `land_area_km2` | 171,923,348.54 | 171,884,385.35 | -38,963.20 km² (-0.023%) |
| `sea_level_m` | -61.380 | -61.200 | +0.179 m |
| `continental_material_live_km3` | 11,487,261,825.12 | 11,487,261,614.98 | -210.15 km³ (-1.8e-8 rel) |
| `craton_balance_error_km3` | 1.024e-06 | 0.0 | -1.024e-06 |
| Runtime | 156.9 s | 159.2 s | +2.3 s (noisy shared host) |
| Peak RSS | 984.6 MB | 908.3 MB | -7.8% |

## Step 3 routing comparison: preregistration

Principle #2 calls for compiled work where per-node routing limits simulation pace. This
comparison tests whether native accumulation and uphill DAG ordering reduce routing time while
preserving flow destinations, ledgers, and reservoir state. The baseline is `f1b75db6`;
the candidate is the issue branch. Both use diagnostic wind, one warmup and three measured
100,000-year steps, sequentially on the same host. The young world uses seed 0 and densities
4/4/2. The mature world is `mantle-bloom-seed47235302-67900000y.bloom`, SHA-256
`eff1f6d7abdfe46c7af5b45ab23f136f3fac591f4d95f3457017c1e7d933216b`.

The primary performance measure is the sum of `Ordering the flow network` and `Routing flow
downstream` phase times; total step time and peak RSS are secondary. Expect a reduction in
routing time, with no material increase in RSS. As in the preceding routing slice, a phase
change below 10% or total-step change below 2% is inconclusive; replicate noisy results.
Exact routing-order, flux, and ledger parity are the principal invariants, with tight numerical
parity for reservoir arithmetic. Any observed model-output change requires a saved-world
replay comparison before claiming end-to-end parity.

For an end-to-end parity check, replay the same mature save for 1 Myr in ten 100,000-year
steps on each arm. Expect identical or numerically negligible differences in land area, sea
level, crust/mantle volumes and conservation ledgers. Compare the final values and their
trajectory at each 0.1 Myr checkpoint; investigate any model-output movement before drawing a
parity conclusion. Use `attribute_hm_ratchet.py` with `--myr 1 --checkpoint-myr 0.1` and run
the arms sequentially on an otherwise idle machine.

### Step 3 phase profiles

`profile_phases.py --scenarios step --steps 3 --warmup 1 --wind-model diagnostic`, wrapped
in `/usr/bin/time -l`, ran sequentially on the same 16 GB host with no other tracked heavy
runs. Routing and ordering are exclusive phase times summed across traced call sites; RSS is
maximum resident set size in decimal MB. The candidate was `105385b6`.

| Case | Mean step, baseline → candidate | Routing + ordering, baseline → candidate | Peak RSS, baseline → candidate |
|---|---:|---:|---:|
| Young, seed 0 | 12.935 → 12.196 s | 0.975 → 0.170 s (−82.6%) | 978 → 975 MB |
| Mature, 67.9 Myr | 14.246 → 13.312 s | 1.120 → 0.209 s (−81.3%) | 987 → 1,020 MB |

Each row is one run per arm; total-step and RSS changes should not be treated as replicated
estimates. The routing result exceeds the preregistered 10% threshold on both worlds. An
initial young baseline produced 13.638 s/step and 0.992 s routing plus ordering, but the
sandbox blocked `time` from obtaining RSS; it was repeated outside the sandbox and is excluded
from the table. The report's timing figures are not model-output parity evidence.

### Step 3 mature replay A/B

The preregistered 1 Myr replay used `attribute_hm_ratchet.py <local-save> --myr 1
--checkpoint-myr 0.1 --out <capture>` on baseline `f1b75db6` and candidate `96319a09`,
sequentially with no other tracked heavy runs. The measured metrics, inventories, mantle
budgets, shortening and mantle closure were exactly equal at every checkpoint; timing fields
varied as expected. The compact checkpoint comparison is in `routing-replay-comparison.json`.

| At 68.9 Myr | Baseline | Candidate |
|---|---:|---:|
| Land area (km²) | 171,862,737.117 | 171,862,737.117 |
| Continental area (km²) | 354,249,271.179 | 354,249,271.179 |
| Sea level (m) | −60.8954356473 | −60.8954356473 |
| Mantle volume (km³) | 40,616,446,488.168 | 40,616,446,488.168 |
| Runtime (s) | 142.845 | 134.946 |
| Peak RSS (MiB) | 1,037.1 | 1,051.0 |

This confirms parity for the reported diagnostics on this save and interval, not for every
possible world. The single-run replay timings are descriptive, not a replicated performance
estimate.

The affected suite passed 2,605 tests, and both 300,000-node routing and river-reservoir stress parity tests passed.
The full stress run was interrupted after 9 passed and one failure: the glacier seed-1 parity
case measured a sea-level difference of `0.01009711276609293 m` against its `0.01 m` limit.
The same test on the unchanged baseline failed with the identical values, confirming that
failure predates this routing port. Remaining full-suite stress cases were not run.

### Gemini review-fix profiles

After the review fixes, the same profile command, worlds, warmup and three-step schedule were
run sequentially on the same host with no other tracked heavy runs. The earlier candidate is
`105385b6`; the revised measurements include native ordering for non-NaN land elevations,
including glacier graphs with infinite off-ice values, and avoid unused erosion buffers.

| Case | Mean step, previous → revised | Routing + ordering, previous → revised | Peak RSS, previous → revised |
|---|---:|---:|---:|
| Young, seed 0 | 12.196 → 12.290 s | 0.170 → 0.172 s | 975 → 904 MB |
| Mature, 67.9 Myr | 13.312 → 13.352 s | 0.209 → 0.195 s | 1,020 → 1,048 MB |

These are single runs per arm; the differences are too small or variable to attribute a
further speed or memory change to the review fixes. The large original routing improvement
remains evident relative to the `f1b75db6` baseline above.

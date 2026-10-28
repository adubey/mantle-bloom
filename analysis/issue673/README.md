# RFC 0014: final acceptance analysis

## Scope and provenance

This report completes the checks that could be evaluated against the three reserved #538
held-out saves. It does not change the accepted lake model. All current-code replays used
backend commit `6ebf8cd19ee0098822ff9eae6685a22e20f9db32`, the ice-age control `never`, and
the preregistered save bytes:

| Seed | Start age | Density | SHA-256 |
|---:|---:|---:|---|
| 373936961 | 244.4 Myr | 4 | `9980103ce9bd7733e4671bc3bec0d48359b45a8e457bc76a73dd769f8840f342` |
| 997271774 | 130.5 Myr | 4 | `a43962a149a0ec96f1dbee3e6c675be98d50422c47e0455e062aa72a7d1da4d7` |
| 456862899 | 19,191.8 Myr | 4 | `e6a6c7514ae52601f1d8d3ae541da3c1dcb727cde141a79b0bee39d3c962e865` |

Each checksum matches the held-out registration on #538. Runs were serial; no test or other
replay overlapped them. Per-step measurements are in the compact JSON files beside this
report. Saves and full replay traces remain outside the repository.

## Results

| Acceptance check | Result | Evidence and disposition |
|---|---|---|
| Listed unit-test cases | **Audited; focused gaps added** | The endpoint-vs-C steady-state case is not applicable to production: the endpoint fast path remains disabled after #589's runtime screen. Nested/asymmetric supply and ledger closure, lake-ice bounds and fixed-forcing plateau, and missing-save `lake_ice_m` default now have focused regressions. Existing tests cover growth/shrinkage budgets, merge/split/perched bodies, hierarchy spill receiver, seasonal ice-free fraction, liquid under ice, melt/sublimation accounting, cold sinks, below-sea-level lake temperature, and silt scaling. The targeted suite passes: lakes, hydrology, persistence, and breaching. |
| Held-out combined ledger residual | **Pass** | Across ten 100 kyr steps per save, maximum `|residual| / positive external supply` was `3.82e-13` (373936961), `6.99e-13` (997271774), and `1.82e-11` (456862899). Each is below 0.1%. The third save has zero-supply steps; its maximum is on a later nonzero-supply step. |
| Water deleted outside the ledger | **Pass** | Across all 30 steps, the cold-water instrument measured zero water deleted by freezing and zero unbooked lake inflow. The step-budget instrument also found zero unbooked balance surplus. This is the represented surface-water/lake/ice budget; subsurface storage remains outside RFC 0014 per merged PR #672. |
| No pinned lakes | **Not applicable to current production** | The production endpoint optimization is disabled; all bodies use prognostic C. Therefore there are no non-spilling, uncapped endpoint bodies that are not flagged to C. |
| Large-body flicker | **Fail** | Seed 373936961 peaked at 1/4 comparable bodies (25%, step 5); seed 997271774 peaked at 1/3 (33.3%, step 8). Seed 456862899 had no body above 10,000 km² and is non-informative. The 5% per-step limit is not met on the two informative saves. The smallest follow-up is to record the identities and merge/split/spill transitions of the qualifying reversals, then determine whether topology churn or storage integration causes them before changing model behavior. |
| 10 kyr vs 100 kyr frozen-forcing solve | **Pass for C-fallback bodies** | Captured forcing was held constant for both solves. Seed 373936961: 19 matched bodies above 10,000 km², all within 5%, maximum 3.95%. Seed 997271774: 9/9 within 5%, maximum 0.90%. Seed 456862899 had no body above the area threshold. No unflagged endpoint bodies exist in the current production solver. |
| Density 4 vs 6 | **Unmeasured; needs a paired held-out input** | The reserved held-out saves are all density 4. The available density-6 save is diagnosis seed 317739385 and cannot substitute for a held-out seed. No same-seed density-6 held-out save is registered, so this comparison cannot be made without generating and preregistering a new paired input. |
| Breached pits and floor-drain toggles vs baseline | **Measured** | The #588 same-save comparison against the task-specific pre-C `origin/main` baseline `a6d33da75d9d28568e9c2c05a243acc033ee3e44` reports these counts on the three diagnosis saves; see [`analysis/issue588/replay-comparison.json`](../issue588/replay-comparison.json). C vs baseline counts were 447 vs 436 and 5,231 vs 5,115 (seed 910211954); 29 vs 29 and 793 vs 797 (635056485); 422 vs 421 and 1,841 vs 1,882 (317739385). The #503/#218 hydrology and breaching tests pass in the targeted test run. |
| Density-4 median runtime vs baseline | **Pass** | Reusing the #588 ten-step same-save comparison against baseline `a6d33da75d9d28568e9c2c05a243acc033ee3e44`, per-seed medians were 17.21 vs 18.71 s and 6.19 vs 6.87 s (C vs baseline), changes of −8.0% and −9.9%. The median of the two per-seed medians changed from 12.79 s to 11.70 s (−8.5%), within the +10% limit. See [`analysis/issue588/README.md`](../issue588/README.md). |
| Save compatibility and `lake_ice_m` default | **Pass** | The targeted persistence suite passes. A new regression removes `lake_ice_m` from a saved world's plate fields and confirms loading yields zero. |
| Simulation documentation | **Pass** | [`docs/simulation-model/hydrology.md`](../../docs/simulation-model/hydrology.md) documents the prognostic lake solver, combined lake/ice ledger, ice cap, water-surface temperature, and numerical depth guards. RFC scope now states that subsurface exchange is not represented. |

### Test audit details

Existing focused tests include:

- `test_prognostic_lake_solver_grows_and_closes_single_basin_budget` and
  `test_prognostic_lake_solver_shrinks_by_evaporation_without_reset`;
- `test_prognostic_lake_solver_merges_siblings_during_the_step`,
  `test_prognostic_lake_solver_splits_below_saddle_without_resetting_water`,
  `test_perched_water_does_not_merge_parent_or_govern_its_spill_rim`, and
  `test_prognostic_spill_uses_hierarchy_receiver_independent_of_flow_channel`;
- `test_seasonal_lake_ice_helpers_match_analytic_cycle_and_are_continuous`, cold-lake
  retention and delivery tests, and `test_compute_hydrology_sublimates_covered_lake_ice_and_closes_its_ledger`;
- `test_lake_temperature_uses_lagged_water_surface_and_connected_ocean_mask`;
- `test_silt_volume_follows_real_inflow_whatever_the_cell_size`.

Three focused tests were added for asymmetric nested source accounting, the lake-ice physical
cap/local-inventory plateau, and old-save `lake_ice_m` defaulting. The inactive endpoint path
is intentionally not tested as production behavior; #589 records why it stays disabled.

## Reproduction

Run from the repository root, with the three registered held-out saves substituted for
`<save>` and one command at a time:

```bash
backend/.venv/bin/python bin/debug/measure_lake_step_budget.py <save> \
  --myr 1 --step-years 100000 --checkpoint-myr 0.1 \
  --out analysis/issue673/heldout-budget-seed<seed>.json

backend/.venv/bin/python bin/debug/measure_cold_sink_water.py <save> \
  --myr 1 --step-years 100000 \
  --out analysis/issue673/heldout-cold-seed<seed>.json

backend/.venv/bin/python bin/debug/measure_frozen_lake_step_size.py <save> \
  --out analysis/issue673/step-size-seed<seed>.json

backend/.venv/bin/pytest backend/unit_tests/test_lakes.py backend/unit_tests/test_hydrology.py \
  backend/unit_tests/test_persistence.py backend/unit_tests/test_breaching.py
```

The density-4 runtime and baseline side-effect measurements are already recorded in
`analysis/issue588/replay-comparison.json`; the baseline commit is recorded above. The held-out
resolution comparison remains open because the preregistered set lacks a density-6 pair.

## Verification

- `backend/.venv/bin/ruff check` — passed.
- `./bin/quality_check.sh` — passed; vulture and knip reported no dead code, and jscpd found
  no new duplicates.
- `./bin/affected_test.sh --base origin/main` — 176 passed across the three modified backend
  test files (46 warnings).
- `./bin/stress_test.sh` — 58 passed (73 warnings).
- Frontend checks were not applicable; no frontend files changed.

## RFC disposition

Keep RFC 0014 **Accepted and open**. The water budget, cold-sink retention, step-size response,
tests, save compatibility, runtime, and documentation checks pass or are explicitly
not-applicable. The held-out flicker threshold fails on both informative saves, and the
held-out resolution comparison lacks a paired density-6 input. Resolve those two items before
marking the RFC implemented.

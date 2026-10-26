# Issue #588: C-default in-step lake fill-spill-merge

## Preregistration

Registered before the diagnostic replay. Compare the #588 branch with a task-specific
detached `origin/main` worktree over ten 100 kyr steps on each of the three diagnosis saves
registered on #538. Use the same save bytes and the ice-age control `never`. The held-out
acceptance saves remain reserved.

| Seed | Age | Density | SHA-256 |
|---:|---:|---:|---|
| 910211954 | 220.2 Myr | 4 | `bd4874092efc805e8b174e115fe567c82d35fc2ee606686938b5b6eb96b060d2` |
| 635056485 | 452.8 Myr | 4 | `ec27d5077bdf43b3a6f0339b5726fa4d110da9eb1f143f96a3004d79f309b045` |
| 317739385 | 200.8 Myr | 6 | `6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd` |

**Expected direction:** the C path should remove split-reset storage jumps and close each
step's lake budget. It will retain finite-step storage and may change lake volume and area
relative to the balanced-level baseline; no directional prediction is registered for those
trajectory measures. The default 64 integration slices may increase runtime; the issue does
not pre-register an allowable runtime change.

**Primary invariant:** on the branch, absolute liquid-ledger residual divided by external
lake supply is at floating-point tolerance per step. Report the maximum normalized residual
and any independent unledgered loss estimate. The baseline has no lake ledger, so its
unledgered loss estimate is measured from opening storage + external supply − evaporation −
exported overflow − closing storage; positive unexplained deficit is reported as loss, and a
negative value as creation.

**Secondary measures:** global visible lake volume and area per step, large-body flicker
reversals using the RFC's definition (consecutive volume-change signs reverse, with both
changes exceeding 50% of the prior volume, for bodies above 10,000 km²), breached-pit counts,
floor-drain toggles, mean step runtime and peak RSS. These diagnose changes and are not
standalone pass/fail thresholds for this issue.

Runs will be serial, with no tests or other replay/probe running. Record the backend commits,
save checksums, commands, runtime and RSS. Full outputs stay in scratch; check in only compact
summary metrics.


## Implementation and replay results

Implementation commit: `56102ba31b57d492dbd5058fd44985d90fe7c1eb`.
Baseline: detached `origin/main`, `a6d33da75d9d28568e9c2c05a243acc033ee3e44`.
Each arm used the same registered save bytes and ran ten 100 kyr steps with the ice-age
control set to `never`. Runs were serial. The compact per-arm/per-seed measurements are in
[`replay-comparison.json`](replay-comparison.json); full outputs remain in local scratch.
The measurement entry point is `bin/debug/measure_lake_step_budget.py` (run from
`backend/`, once per arm and save, with `--myr 1 --step-years 100000`; use `--backend
<baseline>/backend` for the baseline arm).

| Seed | Max C ledger residual / supply | Median runtime C / main (s) | Peak RSS C / main (MiB) | Final visible lake volume C / main (km³) | Final lake area C / main (km²) |
|---:|---:|---:|---:|---:|---:|
| 910211954 | 0.000449% | 20.05 / 18.71 | 5,019 / 5,343 | 3,027,030 / 1,989,244 | 6,955,943 / 3,785,076 |
| 635056485 | 0.025533% | 6.59 / 6.87 | 4,192 / 4,752 | 17,831 / 18,052 | 369,190 / 366,893 |
| 317739385 | 0.000335% | 10.67 / 12.14 | 6,136 / 6,859 | 5,955,695 / 6,711,672 | 6,546,404 / 7,061,193 |

The largest normalized C-ledger residual across all 30 branch steps was `2.5533e-4`
(0.0256% of external supply), below the issue's 0.1% bound. The largest absolute residual
was `8.473e12 m³`, on seed 635056485. Branch replays recorded zero unbooked balance
surplus and zero surplus held at a depth cap. The baseline has no ledger; the diagnostic
proxy summed legacy unbooked surplus of `3.681e17 m³` across the three replays and separately
identified `1.560e18 m³` held at caps. These are estimates from the old level-balance path,
not a direct measurement of physical water loss.

The secondary trajectory measures stayed diagnostic as preregistered. Large-body flicker
reversals were 6/95 comparable transitions versus 11/96 on seed 910211954, 0/10 versus
0/10 on seed 635056485, and 0/24 versus 0/29 on seed 317739385 (C versus main). Runtime
was mixed: C median runtime ranged from 4% faster to 7% slower by save; mean runtime ranged
from 1% faster to 10% slower. Peak RSS was lower on C for all three saves. Final lake volume
and area changed in both directions across saves, with no directional criterion registered.
Large runtime and memory spikes occurred in both arms at similar late steps.

### Protocol note

The preregistration described a 64-slice default. The implementation's C-default uses one
implicit full-step lake solve with same-step hierarchy fill/spill/merge processing; explicit
substeps remain available and are exercised in unit tests. This kept the real-world replay
cost bounded while retaining an in-step storage solve. The replay measurements above are for
the actual C default, not a 64-slice integration. The three diagnosis saves are the only
replay inputs; held-out acceptance saves remain reserved.

### Verification

- `backend/.venv/bin/ruff check` — passed on the merged issue branch.
- `./bin/quality_check.sh` — passed: vulture and knip found no dead code; jscpd found no
  clones.
- `./bin/affected_test.sh --base origin/main` — 1,707 passed in 943.67 s; 408 warnings,
  no failures.
- `./bin/stress_test.sh` — 47 passed, 9 fixture setup errors. All nine errors are in
  `stress_tests/test_main.py`, whose fixture references `app.main._state`, removed by the
  newer multi-world API on `origin/main`. Reproduced the exact `AttributeError` on a fresh
  detached `origin/main` worktree for
  `stress_tests/test_main.py::test_step_advances_elapsed_years`, using the issue worktree's
  Python environment because baseline dependency installation could not reach PyPI. The
  error is present on baseline; the hydrology and world-stepping stress cases passed.
- Frontend checks are not applicable; no frontend files changed. The replay itself completed
  before the unit and stress gates; no replay ran concurrently with tests.

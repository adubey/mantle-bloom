# Issue #875: substep the lake solve so one-step volumes do not depend on slice length

Status: complete. All #875 acceptance items measured below.

Driver: [`bin/debug/measure_lake_slice_convergence.py`](../../bin/debug/measure_lake_slice_convergence.py)
(the #853 probe, unmodified). Compact data: `slice-cap-evidence.json` (this directory).
Raw outputs stay in the private captures directory for this issue.

## Change under test

Governing principle: principle 3 in
[`docs/design-principles.md`](../../docs/design-principles.md) — one long step should
agree with several shorter steps to a stated tolerance. The lake solver charged
evaporation and sublimation at the end-of-slice level for the whole slice, so a
filling lake overpaid on long slices (#853: 13% low at a 5 kyr step, 4.2% at
25 kyr, 1.1% at 100 kyr, single slice).

Fix (no solver equations, state or parameters changed): `compute_hydrology`
passes `substeps=lake_substeps(years)` to `solve_lakes_prognostic`. The slice is
the step's hundredth clamped to [50, 1000] yr (`LAKE_MIN_SLICE_YEARS` /
`LAKE_MAX_SLICE_YEARS` in `hydrology.py`): 50 yr at 5 kyr, 250 yr at 25 kyr,
1000 yr at the 100 kyr production step. Two cost controls follow the glacier
march precedent: a settle-stop inside the solve (two slices moving total volume
by at most `LAKE_SETTLE_TOLERANCE = 1e-4` solve the remainder as one slice), and
the `SUBSTEP_LAKES` test knob (unit suite off by default in `conftest.py`, same
pattern as `MARCH_GLACIERS`; production on).

## Preregistered comparison (before any probe run)

- Checkpoints: window D `seed910211954-225.1myr.bloom` SHA-256 `c5b1f5d8…` and
  window E `seed910211954-224myr.bloom` SHA-256 `95dde684…`, the #853 saves.
- Arms, all fixed-mode (one `step_world`, climate held) from each checkpoint:
  production slice (`lake_substeps` cap), converged slice (200 yr; 10 yr at the
  5 kyr total, finer than the production 50 yr), and single slice (the pre-#875
  production behaviour and built-in baseline: the probe wrapper overrides only
  the solver's `substeps`, so arms differ solely in lake slicing).
- Primary metric: anchor volume after one step, production arm against the
  converged arm; decision threshold 1% (the #875 acceptance). Secondary:
  normalized lake-ledger residual (acceptance ~1e-11), world lake volume,
  merge/split counts, and step wall time (production against single-slice).
- Expected: production within 1% of converged at 5, 25 and 100 kyr; shift from
  single-slice about +13% / +4% / +1% (the #853 single-slice errors, reversed);
  residuals at or below 1e-11; step-time increase material but bounded (the
  #853 probe measured +23% step wall time at 100 substeps without a
  settle-stop).
- Noise: deterministic replay from a fixed checkpoint; no seed variation.
  Runtime comparisons are single runs on a shared machine and are indicative
  only.

## Results

Source: this branch at `cae04802` (branch base `origin/main`) plus the change
above. All arms fixed-mode, deterministic from the stated checkpoint.

| Window | Total | Production (slice) | Converged | Diff | Single slice | Prod shift |
|---|---|---:|---:|---:|---:|---:|
| D | 5 kyr | 393,311 (50 yr) | 393,302 (10 yr) | 0.002% | 348,607 | +12.8% |
| D | 25 kyr | 409,563 (250 yr) | 409,627 (200 yr) | 0.016% | 392,527 | +4.3% |
| D | 100 kyr | 415,171 (1000 yr) | 415,152 (200 yr) | 0.005% | 410,814 | +1.1% |
| E | 100 kyr | 987,494 (1000 yr) | 987,510 (200 yr) | 0.002% | 988,357 | −0.09% |

Anchor volumes in km³, re-measured after the peer-review fixes (L1 settle
criterion, spilling union, count cap); production arms are unchanged to the
km³, converged arms moved slightly toward them. Every production arm is
within 1% of its converged arm — the largest gap is 0.016%. Merge/split
counts barely move between arms (D 5 kyr: 93/46 production against 94/48
converged; elsewhere identical), so the effect remains the loss integration,
not discrete events. The
single-slice arms reproduce #853's values exactly (348,607 / 392,527 /
410,814 km³), confirming the arms differ only in lake slicing. Shift against
single-slice reverses #853's single-slice errors (13% / 4.2% / 1.1% low) as
preregistered. Window E's anchor was already slice-converged at 100 kyr in
#853, and its shift is negligible, also as expected. No arm enters the merged
inland-sea state (anchor catchment 0.82M km² on D).

**Ledger residuals.** Production-arm normalized residuals: 4.8e-12 (D 5 kyr),
2.4e-13 (D 25 kyr), 3.7e-13 (D 100 kyr), 1.5e-13 (E 100 kyr) — at or below the
~1e-11 acceptance level. The finest arm overall (D 5 kyr at 10 yr slices)
reads 1.5e-11, within the range #853 treated as valid (≤5.1e-11).

**Step-time cost** (single runs, shared machine, indicative only): D at 5 kyr
31.8 s single-slice against 44.0 s production (+38%); at 25 kyr 78.8 s against
97.9 s (+24%); at 100 kyr 303 s against 246 s — production *faster* in that
pair, i.e. within run-to-run noise from co-running agents, so no material
100 kyr cost is demonstrated either way. The settle-stop bounds the cost once
lakes equilibrate; the unsliced remainder of a step (tectonics, climate,
routing) dominates these wall times. Because the acceptance tolerance is met
with the substeps alone at every measured step length, the time-weighted-loss
alternative named in the issue is not needed.

**RFC 0016 gate re-read.** Re-reading the window D/E one-step anchor quantity
at the production step on this source: D moves 410,814 → 415,171 km³ (+1.06%,
the expected ~1%); E moves 988,357 → 987,494 km³ (−0.09%). Neither window is in
the merged state before or after. Registered gate volumes in
`analysis/rfc0016/README.md` were measured under the earlier target-volume
definition on an older source, so the comparable statement for the tracker is
the shift, not a replacement of those numbers.

**Tests.** `test_lakes.py` + `test_hydrology.py`: 242 passed. Affected suite
(`./bin/affected_test.sh --base origin/main -n 4`): 2424 passed in 6:51 — the
first attempt at this run, before the `SUBSTEP_LAKES` suite knob, ran past
1.5 h unfinished; the knob (suite off, production on, glacier-march pattern)
is what restored suite time. Ruff and `quality_check.sh` clean.

# Full-duration losses after lake transfers

Issue [#768](https://github.com/adubey/mantle-bloom/issues/768) addresses principle 3:
transfers settled after the timed pass must participate in the receiving lake's implicit
water balance. Principle 2 requires measuring the cost of revisiting affected bodies.

## Preregistration

Recorded before any replay. Baseline: detached `origin/main`, `dea1f6b5` (PR #761 merged).
Branch: the commit containing this preregistration, on `codex/issue-768-lake-transfer-losses`.

The implementation restores a body's previously booked slice losses before re-solving its
current storage and new transfers over the same slice duration. The new loss replaces the
old charge, including losses charged to children before a merger. External supply enters only
the timed pass. Every slice settles transfers and topology before the next slice begins.
The prognostic active set also requires storage sufficient for a common saddle surface;
floor-cell levels alone can request a merger that immediately splits when redistributed.

**Equation check.** A three-cell fixture starts split, with two 0 m floors, a 4 m shelf,
unit cell areas and a 1 m joining saddle. Supply is 30 m³ over one year and combined loss
is 1 m/yr. The final fully flooded lake must have 27 m³, at level 31/3 m. On unmodified
`dea1f6b5` it holds 29 m³ at 11 m, while its ledger still closes. Evaporation, sublimation,
and their combination must each match the independent answer within 1e-8 m and m³.
Changing which sibling runs first must preserve that answer. An unchanged receiver getting
a late transfer must replace its loss rather than charge another full-duration loss.

**Windows.** Each arm runs 12 steps of 100 kyr (1.2 Myr), ice-age cycle off. E uses #757's
checkpoint and anchor. The owner selected S200.8 as a replacement for the missing original S
checkpoint before any replay; it uses #757's S anchor.

| Window | Seed | Opening age (Myr) | SHA-256 | Anchor |
|---|---:|---:|---|---|
| E | 910211954 | 224.0 | `4045b915a2d3d9e534a84089bba7277bbd488690e51d70f4e928f58c8cf674b1` | `analysis/rfc0016/target-anchor-window-e.json` |
| S200.8 | 317739385 | 200.8 | `6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd` | `analysis/issue757/anchor-seed317739385-201.6.json` |

Both input SHA-256 values match. S200.8 spans 200.8–202.0 Myr; it is a matched comparison
on the owner's replacement, not a reproduction of original S at 201.6–202.8 Myr.

**Primary metric and decision.** `--audit-losses` captures the solved water column before
the ice split. For every wet final body, compare its trace's booked evaporation and
sublimation against the prescribed rates times final wet area times the full step duration.
The sum of absolute differences, divided by combined expected loss (minimum 1 m³), must be
at most 1e-9 on the branch. Report E's anchor at every step, especially the baseline's refill
steps 2 and 10. Dry bodies are excluded: a supply-limited body can consume water while ending
with zero wet area. The hand-calculated fixture remains the independent equation-level check.

**Accounting gates.** Lake and rain-only land normalized residuals must be <= 1e-9 every
step. Liquid and ice sub-ledgers must close. All lake merge, split, spill, capped-cycle and
silt-accounting tests must pass. A bounded-settling failure is a failure, not a usable result.

**Predictions, no magnitude threshold.** E's artificial refills should shrink because they
now face full-duration losses. The target may still swing through lagged climate and rain
partitioning; neither zero swings nor identical endpoints from different opening states is
required. Report whole-world liquid, ice and combined volume, lag-1 autocorrelation of volume
changes, mean absolute step change from step 3, sea level, largest lakes, runtime and peak RSS.
Do not infer stability from autocorrelation alone: 12 steps contain few independent cycles.
These are deterministic comparisons on two seeds, not an estimate of planetary variability.

**Controls.** The baseline should reproduce #757's final E branch (model commit `e5c697d2`):
step 2 target volume about 27.87 M km³, step 10 about 28.32 M km³, and final whole-world
liquid about 14.30 M km³. Compare exact values from its compact output where available;
use a 0.01% tolerance against unrounded data. Record the SHA and commit for every run.

**Runtime.** Follow `compute-budget` before each arm and run arms sequentially. Prefer no
other heavy run on the machine; disclose overlap. E previously peaked at 4.5 GB, S at 7.1 GB;
reserve 5 and 8 GB respectively. An increase above 25% in comparable median step runtime
triggers investigation before delivery; it is a diagnostic threshold, not proof of a regression.
Before the full runs started, the owner explicitly chose to proceed whenever memory permits
with overlapping-work timings labelled. The four arms therefore require a successful memory
check, but no longer wait for an otherwise idle machine. Runtime differences from overlapping
arms cannot establish a speedup or slowdown, so the 25% clean-timing investigation threshold
does not apply to those arms.
E main completed at 4.14 GB; both full E arms started with a 5 GB reservation. S200.8
reserves 7.5 GB against the previously measured 7.1 GB S peak. This includes the small
transient audit arrays; the shared budget check retains its own machine headroom.
While waiting for a quiet machine, a two-step E branch probe may run after a successful budget
check to test settling and the loss audit at the first refill. Its timing is excluded from
the runtime comparison; the run and any failure will still be recorded.
The first two-step probe peaked at 1.25 GB; repeat probes of that same short horizon reserve
2 GB. Full E comparisons retain the 5 GB estimate from the longer #757 runs.

```bash
backend/.venv/bin/python bin/debug/replay_lake_water_swing.py <checkpoint.mbworld> \
  --anchors <anchor.json> --myr 1.2 --step-years 100000 --checkpoint-myr 0.1 \
  --audit-losses --backend <arm-tree>/backend --out <captures>/issue768-<window>-<arm>.json
```

Full outputs remain outside the repository; only compact results and findings are committed.

## Preflight verification

- On unmodified `dea1f6b5`, all three loss variants of the new full-step fixture fail:
  final level 11 m rather than 31/3 m, even though the accounting closes.
- The first affected run exposed an oscillating merge/split in
  `test_faults_accumulate_with_a_mix_of_regimes_and_lock_up_as_scars`; it was interrupted
  after 424 passes and that failure. A targeted run reproduced the bounded-settling error.
- A tentative shared merge guard fixed that failure but changed three legacy routing/event
  tests. That diagnostic run was interrupted after 100 passes. The guard now applies only
  to the prognostic active set; the legacy behavior is preserved.
- The scoped revision passes all 88 lake tests, including partial-cell balances, both sibling
  processing orders, late transfers without a merger, loss-trace totals and the saddle-volume
  regression. Ruff and the dead/duplicate-code gates pass. The full affected run on code
  revision `65d27be8` passes: **1,969 tests**, 512 warnings, 817.54 s, two workers. Its
  command was `./bin/affected_test.sh --base origin/main -n 2 -x --ff`.

## Functional probes

The first two-step E probe ran on `b715e7aa`, starting with no other heavy run. It settled,
with worst normalized lake and land residuals 3.3e-12 and 4.9e-13. The target's step-2 liquid
volume is 464,680 km³; #757's recorded main-equivalent result is 27.87 M km³. Peak RSS was
1.25 GB. The loss audit was initially evaluated after erosion, which can mutate its shared
terrain input; its mismatches are not accepted as evidence of solver error or success.
The audit now evaluates inside the solve-return wrapper, before those mutations, and the
same probe will be repeated. Both attempts' full outputs remain in the private capture set.

The second probe (`3907747f`) preserves those water results and evaluates before erosion.
E's anchor agrees with final-level demand to 2.7e-17, but the worst nearly dry bodies show
2.8e-8 relative error. The audit now derives partial wet fraction directly from stored mean depth
(`wet = sqrt(depth / relief)`, bounded by 1), avoiding that cancellation. It also reports
the within-body water-level spread and the five largest relative errors. The 1e-9 loss
threshold is unchanged; a third probe will check this numerically stable evaluation.

The third probe (`b536ebbb`) shows that direct depth inversion does not remove the remaining
raw error. The worst body holds 0.0094 m³ at a 4,165 m absolute level, against roughly
39.9 million m³ of demand; the booked loss is 1.1 m³ below final-level demand. The solve caps
loss at available water, and the absolute level cannot represent an arbitrary demand exactly.
Raw loss errors remain reported, including the original 1e-9 gate miss. An additional audit
evaluates demand over a two-ULP interval around the reconstructed level (one rounding in the
solve's absolute level, one in the diagnostic inversion), and reports any discrepancy outside
that interval. This is a numerical representation bound, not a changed model parameter or
denominator floor. Its uncovered discrepancy is checked against the same 1e-9 threshold.

The fourth probe (`d00ff3b9`) covers the worst raw error with the representation interval,
but the largest uncovered error at step 2 is 6.54e-9 in a body holding only 6.1 m³. The
existing multi-cell integrator accepts a volume-balance residual up to 1e-8 of available
water before invoking its bisection fallback. The extra all-body 1e-9 loss gate is therefore
stricter than that numerical solve and is **not met**; it remains reported as such. This
does not change the agreed 1e-9 lake/land accounting gates or the independent unit-level
equation test. Full comparisons will report raw errors, rounding intervals and the anchor's
loss balance separately, with this numerical limit explicit. No solver precision parameter
has been changed to obtain a passing audit.

## Full comparison results

### E baseline control

The 12-step baseline on `dea1f6b5` reproduces #757's final E branch. Against that report's
rounded compact data, the largest whole-world liquid difference is 0.50 km³ and the largest
target difference is 0.04 km³, within the registered 0.01% control tolerance.

- Target liquid is 27.87 M km³ at step 2 and 28.32 M km³ at step 10. The loss audit shows
  demand shortfalls of 73.7% and 56.4% at those refills despite closing global accounts.
- Final whole-world liquid is 14.298 M km³, ice 25,229 km³, sea level -773.201 m.
- Lag-1 autocorrelation is -0.286; mean absolute step change is 6.686 M km³.
- Worst normalized lake and land residuals are 5.84e-12 and 1.34e-12; liquid and ice
  residuals are 5.84e-12 and 2.09e-16. All accounting gates pass.
- Median step time is 8.36 s; peak RSS is 4.14 GB. Two other unit workers were listed at
  the start, so runtime is reported under the owner's overlapping-load caveat.

### E comparison

| Measure | Main | Branch |
|---|---:|---:|
| Mean absolute liquid step change (km³) | 6,686,041 | 151,998 |
| Lag-1 autocorrelation of liquid changes | -0.286 | +0.367 |
| Final whole-world liquid (M km³) | 14.298 | 3.960 |
| Final lake ice (km³) | 25,229 | 14,306 |
| Final sea level (m) | -773.201 | -754.885 |
| Final target area (M km²) | 11.092 | 5.644 |
| Final target liquid (M km³) | 11.651 | 1.875 |
| Step-2 target liquid (M km³) | 27.873 | 0.465 |
| Step-10 target liquid (M km³) | 28.317 | 2.028 |
| Median instrumented step time (s) | 8.36 | 7.76 |
| Peak RSS (GB) | 4.14 | 4.88 |

The liquid step-change metric falls 97.7%. The artificial refills are removed; later growth
and changes through lagged climate remain. At the baseline's refill steps, the branch's
booked losses match the final-level demands to 2.7e-17 and 7.3e-17. Across all 12 branch
steps the anchor's worst normalized loss error is 3.4e-16. Main's corresponding refill
shortfalls are 73.7% and 56.4%.

All branch accounting gates pass: worst lake and land residuals 4.72e-12 and 9.70e-13;
liquid and ice sub-ledgers 4.72e-12 and 2.11e-16. The largest branch target ends below the
3,500 m guard; another large lake remains on that guard, as the endpoint records show.
Runtimes include the loss audit and overlap different unit/stress workloads, so the median
values do not establish a speedup. The measured memory increase is retained in the table.

**Numerical audit limit.** The extra all-body 1e-9 raw loss gate is not met. The branch's
largest raw normalized error, 0.537, comes from a body storing 2.1e-18 m³ with about 0.55 m³
of calculated demand; it lies inside a 1.66 m³ two-ULP demand interval. The largest error
outside the interval is 1.52e-8 in a body holding 6.28 m³, consistent with the pre-existing
finite-iteration precision issue identified in the probes. These flags are retained, not
counted as passing. The independent equation fixture and the materially wet E anchor pass.

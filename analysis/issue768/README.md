# Full-duration losses after lake transfers

Issue [#768](https://github.com/adubey/mantle-bloom/issues/768) addresses principle 3:
transfers settled after the timed pass must participate in the receiving lake's implicit
water balance. Principle 2 requires measuring the cost of revisiting affected bodies.
The saddle-volume check follows principle 4: floor samples must represent enough finite
cell volume to support the proposed common water surface.

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

The public per-step trajectories, anchor audits, numerical extremes and endpoint lakes are
in [replay-comparison.json](replay-comparison.json). Full outputs remain in the private
capture set with SHA-256 checksums.

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

### Owner-selected S200.8 comparison

This starts from the supplied 200.8 Myr save and ends at 202.0 Myr. It is a new matched
window, not the original 201.6–202.8 Myr S replay.

| Measure | Main | Branch |
|---|---:|---:|
| Mean absolute liquid step change (km³) | 578,594 | 445,460 |
| Lag-1 autocorrelation of liquid changes | -0.057 | +0.856 |
| Final whole-world liquid (M km³) | 8.510 | 7.364 |
| Final lake ice (km³) | 11,682 | 10,795 |
| Final sea level (m) | -1113.568 | -1111.433 |
| Final target area (M km²) | 5.154 | 4.797 |
| Final target liquid (M km³) | 7.516 | 6.440 |
| Median instrumented step time (s) | 11.25 | 11.32 |
| Peak RSS (GB) | 7.28 | 6.02 |

The liquid step-change metric falls 23.0%. At step 2, main's target body charges 49.8%
less than its final-level demand; the branch matches that demand to machine precision.
Across the branch's 12 steps the anchor's largest normalized loss error is 3.5e-16.
The visible wet component under the anchor is absent at step 1, so the script's enormous
consecutive-target ratio against a zero denominator is not interpreted as a physical ratio.

All four accounting gates pass in all **48** full-run steps. Across both windows and arms,
the worst normalized residuals are 3.61e-11 for lakes and liquid, 3.13e-12 for land, and
2.11e-16 for ice, below the agreed 1e-9 bound. S main started without another listed heavy
run; another simulation was running by the branch's start. Runtime remains qualified by
the owner's timing caveat.

S also fails the extra all-body pointwise 1e-9 loss gate in essentially dry stores. Its
worst uncovered raw-rate discrepancy occurs in a body storing 5.9e-12 m³, against about
333 m³ of calculated final-level demand and availability-limited booked loss. The two-ULP
interval does not explain that entire discrepancy. This appears associated with the existing
finite-iteration solve near dryness; it has not been resolved or counted as passing here.
The measured anchor is materially wet and passes throughout. These auxiliary flags and their
underlying quantities are retained in the compact JSON.

## Review follow-up

The [posted review](https://github.com/adubey/mantle-bloom/pull/772#issuecomment-6592623888)
approved without blocking findings. Its follow-ups were full stress coverage for fail-hard
settling, explicit capped-cycle timing under `substeps > 1`, and retaining the runtime caveat.

The cycle timing is intentional: every integration slice settles and exports its closed-cycle
surplus before the next slice adds supply. The model guide and ordering helper's docstring
now state this. A permanent positive-loss regression checks `substeps=1,4,16` against the
closed-form balance for two constantly full bodies: evaporation and sublimation each charge
once, storage remains fixed, and supply minus losses is exported to the parent outlet.
All 91 lake tests pass on this revision.

A separate independent reviewer found no actionable defect in 350 new synthetic cases:
150 nested merges against independent final-level/budget calculations and 200 mixed
split/merge cases with changed processing order. It did not rerun the reported author suite.
The probes used `np.random.default_rng(772)` and `np.random.default_rng(7721)`, respectively,
with four two-cell leaves joining at 3 m and 4 m, then a root joining at 6 m. The first
constructed supply for a selected 7–20 m final level using random positive areas, relief,
evaporation and sublimation over 0.1–1,000 years; it enforced 1e-6 m level and 1e-9 normalized
ledger bounds. The second started merged at 7–10 m, varied forcing over 1–20 years, and
permuted lake IDs with equal floors; it enforced 1e-6 m depth agreement between orders.
That second probe did not separately assert ledger residuals. Measured worst errors were not
retained; these are passing assertion bounds, not reported maxima.

Current main `7c26c2cb` was merged before the added stress coverage. Its model change is the
random-world generation default; the previously replayed saved-world stepping code is
unchanged. The full stress suite is running with one budget-approved worker.
The run completed successfully: **59 passed**, 61 warnings, 1,099.34 s, using
`./bin/stress_test.sh -n 1 -x`. A fresh affected-test gate on the merged revision follows.

## Outcome and limits

The independent unit fixture fails on main and passes on the branch, with evaporation,
sublimation and both together, either sibling order, and partial-cell geometry. The replay
anchors confirm that the intended full-step losses are charged at the formerly artificial
refills. Global and phase accounting meet the acceptance bounds. The extra all-body
pointwise accuracy gate is not met at numerical near-dry states and is explicitly disclosed.

The windows retain lagged climate feedback and cover 1.2 Myr each. Largest-lake guards remain
in effect. If #758 changes these basins' ocean classification before this lands, its affected
replays need updating against that new main. No persisted fields or frontend code changed.
Verification includes the affected suite, registered saved-world comparisons, and the
59-test full stress suite added after review.

## Five largest lakes at each endpoint

Rank is by wet-component area, as in #757. Volumes below are liquid; ice and the full per-step traces are in the compact JSON.

| Window / arm | Rank | Area (M km²) | Liquid (M km³) | Level (m) | Floor (m) | At 3,500 m guard |
|---|---:|---:|---:|---:|---:|---|
| E main | 1 | 11.092 | 11.651 | -3322.8 | -5433.1 | no |
| E main | 2 | 1.504 | 0.154 | 665.3 | 523.9 | no |
| E main | 3 | 1.434 | 0.789 | 1819.9 | -644.8 | no |
| E main | 4 | 0.677 | 1.001 | -741.6 | -4089.3 | no |
| E main | 5 | 0.431 | 0.030 | 557.3 | 213.0 | no |
| E branch | 1 | 5.644 | 1.875 | -4438.4 | -5418.0 | no |
| E branch | 2 | 1.504 | 0.154 | 665.3 | 523.6 | no |
| E branch | 3 | 0.827 | 0.116 | -4439.4 | -4723.0 | no |
| E branch | 4 | 0.694 | 1.099 | -578.5 | -4078.5 | yes |
| E branch | 5 | 0.431 | 0.030 | 557.3 | 213.0 | no |
| S200.8 main | 1 | 5.154 | 7.516 | 865.1 | -1444.0 | no |
| S200.8 main | 2 | 1.061 | 0.871 | 2925.7 | 327.1 | no |
| S200.8 main | 3 | 0.314 | 0.049 | -876.3 | -1620.6 | no |
| S200.8 main | 4 | 0.154 | 0.015 | -407.5 | -717.2 | no |
| S200.8 main | 5 | 0.114 | 0.009 | 1283.2 | 1074.1 | no |
| S200.8 branch | 1 | 4.797 | 6.440 | 646.9 | -1438.3 | no |
| S200.8 branch | 2 | 0.986 | 0.810 | 2666.9 | 471.4 | no |
| S200.8 branch | 3 | 0.311 | 0.049 | -876.3 | -1620.6 | no |
| S200.8 branch | 4 | 0.139 | 0.002 | -991.5 | -1052.4 | no |
| S200.8 branch | 5 | 0.114 | 0.009 | 1283.2 | 1074.1 | no |

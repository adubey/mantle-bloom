# Issue 973: RFC 0015 D3 implementation

Implements the accepted D3 direction in [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md#d3-glacial-erosion): integrate Cook et al. (2020)'s physical basal-sliding abrasion law through the D2 ice march, update a working bed between intervals, and transfer accumulated rock once through the existing terrain and material ledgers.

**Principles:** #5 (consistent units) and #8 (each process follows its own physics). The prior closing-ice-times-world-step estimate was then clipped with the shared neighbor-drop limit, creating a 66–73% gap between 0.25 Myr and 1 Myr steps on G3's diagnosis saves. The invariant is that each march interval contributes the physical law times interval duration, the next interval receives the bed lowered by the previous interval's abrasion, no held interval omits continuing abrasion, and the final material transfer remains within crustal headroom and closes the existing sediment ledger.

## Pre-registered replay

Before replaying, expect the integrated abrasion total at 0.25 Myr and 1 Myr to differ materially less than G3's prior 66–73% gap; a gap under 20% is the working target, while a remaining gap above 50% means D3 has not removed the dominant timestep sensitivity. Expect applied/raw to rise above G3's 0.0045–0.210 range where the removed material has available crustal headroom; any remaining loss should be explainable by headroom, mobile cover and craton resistance. Runtime is a measured cost rather than a pass threshold; compare each D3 arm with the D2 baseline under the same save, horizon and step schedule, and report slowdown alongside the interval count.

Controls: Cook et al. (2020), coefficient 1e-4 m/yr, exponent 0.65, glacier_erosion_multiplier = 1, ice ages disabled, same 10 Myr horizon and deterministic save for each paired run. Use each save at 0.25 Myr × 40 steps and 1 Myr × 10 steps. D2 baseline runs use origin/main; D3 runs use the issue branch. Do not tune coefficients after seeing these results.

### Post-merge validation preregistration (2026-11-19)

The four D3 arms below will rerun on the merged implementation (`aa217ab4`, PR #976)
because the prior full D3 table predates the shared-headroom, craton-resistance and ocean-mask
corrections. The D2 controls remain the previously recorded runs at `382ddccd`; no D2 rerun is
planned. Before any arm runs, the gate script is instrumented to record each step's march
interval count as well as its existing measures. That diagnostic-only change does not alter
the simulation model.

| Save | SHA-256 | Step schedule | Arms |
|---|---|---|---|
| `mantle-bloom-seed649658380-6400000y.bloom` | `06167cb4e45fac7141e98df2cf091999ffc40b36be4c0da3bd94703b6b31131f` | 0.25 Myr × 40; 1 Myr × 10 | D3 Cook2020 |
| `mantle-bloom-seed997271774-87400000y.bloom` | `07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8` | 0.25 Myr × 40; 1 Myr × 10 | D3 Cook2020 |

Controls are `K_G = 1e-4`, exponent `0.65`, glacier and flatten multipliers `1`, ice-age
period `0`, and 10 Myr per arm. D3 uses the merged `trend-1k` march. Existing D2 control
rows use the fixed-250 yr march at `382ddccd`, the same saves, step schedules and horizon.
Each save is deterministic; there is one run per arm, with no replicate planned unless a
result lands close enough to a threshold that machine or numerical variation could change
the decision.

Primary metric: for each save, the quarter-to-annual applied-abrasion gap is
`(applied_0.25Myr - applied_1Myr) / applied_0.25Myr`. The expected direction is below the
G3 Cook-law gaps (66.5% and 70.4%). The existing working target is <20% on each save; >50%
means D3 has not removed the dominant step-size sensitivity. A 20–50% result is partial and
will be reported without claiming the working target passed. Secondary measures are raw and
applied abrasion, applied/raw, flattening volume, march-to-transfer difference, closing ice
area/volume and land volume. Runtime, peak RSS and total march intervals will be reported
against the existing D2 control rows. No coefficients or thresholds will be tuned after
results are observed.

This preregistration is committed before the replays. Full arm JSON and checksums remain in
the private capture directory; only compact per-arm summaries will be checked in.

**Runtime-only protocol amendment (2026-11-19, before the remaining arms).** The first
seed649658380 / 0.25 Myr arm started after an empty-run budget check, but a separate unit-test
process appeared during that replay and was detected immediately after it finished. Keep its
deterministic abrasion and state measures as the primary result, but flag its elapsed time as
non-comparable. Once the unit run exits and a fresh budget check lists no other runs, repeat
that one arm solely to obtain an uncontended runtime and RSS measurement. Preserve both
captures and checksums; do not include the repeated values a second time in the abrasion gap.

**Runtime-only protocol amendment 2 (2026-11-19, after all four primary arms).** The other
three arms were run under explicit owner overrides while independent simulations were active,
so all four observed D3 elapsed times are load-affected. Do not compare any of those times to
D2. If the runtime comparison is needed to close #973, repeat each of the four D3 arms once
with no other runs listed by `compute_budget.py`; use those repeats only for runtime/RSS, not
for abrasion totals, gaps or closing states. Keep the existing D2 runtimes at `382ddccd` as the
controls. This amendment is about runtime validity only and does not change the preregistered
scientific metrics or thresholds.

| Save | Age | SHA-256 |
|---|---:|---|
| mantle-bloom-seed649658380-6400000y.bloom | 6.4 Myr | 06167cb4e45fac7141e98df2cf091999ffc40b36be4c0da3bd94703b6b31131f |
| mantle-bloom-seed997271774-87400000y.bloom | 87.4 Myr | 07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8 |

### Results

**Superseded by the post-review correction below.** These full-run numbers were captured at
`ff61e09f`, before the review-driven revision of the march/transfer material accounting. The
current branch now reserves abrasion against one shared Hc pool and transfers the same
material, so those D3 numbers do not measure the current implementation and cannot close G3.
The land/ocean-mask inversion and double craton resistance were introduced in later review
fixes and corrected afterward; neither was present in the `ff61e09f` replay. The D2 control
rows remain valid.

All eight preregistered 10 Myr replays completed: D2 baseline at `origin/main` (`382ddccd`)
and D3 at implementation commit `ff61e09f`, for each save at 0.25 Myr × 40 and 1 Myr × 10.
The captures and their `SHA256SUMS` are retained outside the repository. No other compute
runs were listed during these replay arms.

| Save | Arm | Schedule | Raw abrasion (km³) | Applied (km³) | Applied/raw | Runtime | Peak RSS |
|---|---|---:|---:|---:|---:|---:|---:|
| seed649658380 | D2 | 0.25 Myr | 4,217 | 17,274 | 4.096* | 331 s | 1.09 GB |
| seed649658380 | D2 | 1 Myr | 4,384 | 17,304 | 3.947* | 82 s | 1.23 GB |
| seed649658380 | D3 | 0.25 Myr | 7,054,389 | 7,054,389 | 1.000 | 1,760 s | 1.05 GB |
| seed649658380 | D3 | 1 Myr | 5,932,272 | 5,932,272 | 1.000 | 1,139 s | 1.03 GB |
| seed997271774 | D2 | 0.25 Myr | 169,374 | 457,849 | 2.703* | 624 s | 1.50 GB |
| seed997271774 | D2 | 1 Myr | 177,328 | 349,649 | 1.972* | 259 s | 1.58 GB |
| seed997271774 | D3 | 0.25 Myr | 135,680,425 | 100,053,493 | 0.737 | 2,262 s | 1.36 GB |
| seed997271774 | D3 | 1 Myr | 99,720,997 | 65,953,273 | 0.661 | 1,599 s | 1.42 GB |

The D3 quarter-to-annual abrasion gap is 15.9% on seed649658380 and 34.1% on seed997271774.
Thus the first save meets the <20% working target and the second does not; both are below the
50% failure threshold, and both improve on G3's 66–73% Cook-law gap. D3 removes the dominant
step-size sensitivity on these cases, while the larger save still has meaningful residual
sensitivity. Applied/raw is essentially 1.0 on the first save. On the second it is 0.66–0.74,
consistent with limited crustal headroom at some locations; max-depth guard warnings occurred
in the quarter-step run and should be interpreted alongside that result.

`*` D2 used the legacy baseline law and its raw/applied measures are not comparable to the
Cook-law D3 volumes; these rows serve as same-save runtime and closing-state controls. Relative
to D2, D3 runtime was 13.9× (annual) and 5.31× (quarter) on seed649658380, and 6.17× (annual)
and 3.63× (quarter) on seed997271774. This cost reflects continued interval marching after
ice-equilibrium diagnostics so each interval's abrasion and bed feedback are represented.

| Save | D3 0.25 Myr closing land ice area (km²) | D3 1 Myr closing land ice area (km²) | D3 0.25 Myr closing land volume (km³) | D3 1 Myr closing land volume (km³) |
|---|---:|---:|---:|---:|
| seed649658380 | 2,948,871 | 2,755,594 | 88,442,666 | 89,406,672 |
| seed997271774 | 9,737,452 | 12,238,086 | 255,421,127 | 259,873,597 |

These are two deterministic diagnostic saves, not a broad stability survey. Their D3 gap and
runtime conclusions apply only to `ff61e09f`; the current D3 step-size gap and runtime remain
unknown until the four arms are rerun.

### Post-review correction

Reviews found that the land/ocean mask in `erosion.apply_erosion` was inverted, craton
resistance was applied twice at transfer, and marching treated mobile cover and substrate
headroom as additive despite both belonging to Hc. The current code in `fb0c8772` uses one
shared Hc pool, applies resistance once, and retains abrasion when the march is disabled by
integrating the closing solve without bed feedback. A one-step 0.25 Myr smoke replay on
seed649658380 at `50cb88df` reported 262,419.27 km³ raw and applied abrasion (applied/raw
1.000; 43.7 s, peak RSS 0.89 GB). That capture predates the shared-pool and no-march fixes and
is preliminary only. It does not replace the four full D3 replay arms or validate the current
head; the prior D3 replay table is invalid until the arms are rerun.

At `0e9fa3b4`, a one-step 0.25 Myr smoke replay on seed997271774 with the `trend-1k` preset
reported 4,570,363.04 km³ raw and 3,294,829.58 km³ applied abrasion (applied/raw 0.7209;
49.5 s, peak RSS 0.95 GB). The march-applied volume and final transfer agree within 0.5 m³.
The save has ice on craton-bearing nodes, so the lower applied share is expected to include
craton resistance; this one step does not isolate the share due to resistance from headroom or
cover limits. It is a targeted smoke check, not a regression test or full G3 validation.

### Post-merge replay results (2026-11-19)

All four preregistered 10 Myr D3 arms ran from `b74d88af`; the production simulation tree is
the merged PR #976 implementation at `aa217ab4` (the follow-up commits only add evidence and
diagnostic interval counting). The D2 runtime controls below remain the same-save, same-horizon
runs at `382ddccd`.

| Save | Schedule | Raw abrasion (km³) | Applied (km³) | Applied/raw | Flattening (km³) | March intervals | D3 elapsed* | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| seed649658380 | 0.25 Myr × 40 | 6,897,278 | 6,897,278 | 1.0000 | 2,035,533 | 17,437 | 1,715 s | 1.05 GB |
| seed649658380 | 1 Myr × 10 | 5,750,860 | 5,750,860 | 1.0000 | 2,021,584 | 12,553 | 1,257 s | 0.96 GB |
| seed997271774 | 0.25 Myr × 40 | 143,832,735 | 111,931,582 | 0.7782 | 14,651,958 | 17,842 | 2,170 s | 1.35 GB |
| seed997271774 | 1 Myr × 10 | 102,514,794 | 71,136,782 | 0.6939 | 11,267,660 | 14,663 | 1,647 s | 1.30 GB |

| Save | 0.25 vs 1 Myr applied-abrasion gap | G3 Cook-law gap | Working target |
|---|---:|---:|---|
| seed649658380 | 16.62% | 66.5% | Pass (<20%) |
| seed997271774 | 36.45% | 70.4% | Partial (20–50%; below the >50% failure threshold) |

Both saves improve on G3's original step-size gap, and neither reaches the preregistered
failure range. The smaller save meets the <20% working target; the larger save does not, so D3
is a partial result against that target and should not be reported as a full pass. Applied/raw
is 1.0 on the smaller save and 0.69–0.78 on the larger save, consistent with material limits
remaining active on the larger world. The maximum absolute march-to-transfer difference is
2 m³ across the four arms.

**Runtime is not comparable.** Other agents' simulations overlapped the D3 arms, so the elapsed
times above are recorded but must not be divided by the D2 controls to claim a slowdown. D2
runtime controls are 331/82 s for seed649658380 (0.25/1 Myr) and 624/259 s for seed997271774.
The D3 runtime comparison remains open until D3 arms are timed without overlapping runs. Peak
RSS and interval totals are reported descriptively. The machine also logged the existing 5 km
glacier-depth guard during the seed997271774 replays; those warnings are part of the observed
outcome. Because the 36.45% gap is in the preregistered partial band, #973 remains open for the
owner's decision on that result as well as the runtime comparison.

The per-step compact summaries are in `d3-postmerge-results.json`. Full captures remain in the
private capture directory and are checksummed in its `SHA256SUMS` file.

### Runtime-only follow-up (2026-11-19)

All four runtime-only repeats have now been captured with no other compute jobs active. These
are the valid D3/D2 timing comparisons; the two earlier overlapping seed997271774 attempts
remain in the private capture manifest but are excluded.

| Save | Schedule | D3 runtime | Peak RSS | D2 runtime | D3 / D2 |
|---|---:|---:|---:|---:|---:|
| seed649658380 | 0.25 Myr × 40 | 1,530.7 s | 1.03 GB | 331 s | 4.62× |
| seed649658380 | 1 Myr × 10 | 1,054.22 s | 1.11 GB | 82 s | 12.86× |
| seed997271774 | 0.25 Myr × 40 | 2,000.74 s | 1.65 GB | 624 s | 3.21× |
| seed997271774 | 1 Myr × 10 | 1,643.34 s | 1.53 GB | 259 s | 6.34× |

The captures and SHA-256 values are:

- `d3-postmerge-runtime-only-seed649658380-025-retry.json`: `04485eb7170eb6882f089c4848fdd851c6a988a27c578f0325e90018cf8d3aa5`
- `d3-postmerge-runtime-only-seed649658380-100-repeat.json`: `adee35d3678efe0aae995ef73c482b4cf3c2f177a620e3da4a2e734a60f93858`
- `d3-postmerge-runtime-only-seed997271774-025-uncontended.json`: `8acaa566c756ecda5e70ec7ada6c5de523b8cea5e6d1697d3c3dc48d7cde95fa`
- `d3-postmerge-runtime-only-seed997271774-100-uncontended.json`: `9f278422f51a7f74bcfc10f532f1c796174c5c74fb055bd43b1de033a2fd5feb`

These repeats are for runtime/RSS only; their deterministic science values do not enter the
preregistered abrasion results a second time. The clean runtime matrix is complete.

### Implementation checks

- `ruff check`: passed.
- `quality_check.sh`: passed (vulture, knip and jscpd; zero new duplicate clones).
- `git diff --check`: passed.
- Focused unit tests: `test_hydrology.py` 155 passed; `test_erosion.py` 64 passed.
- Post-merge GitHub backend and frontend checks passed on PR #976 / merge commit `aa217ab4`.
- The later local post-review unit command accidentally ran the whole backend suite: 2,643
  passed, 1 failed and 8 errored. All hydrology and erosion tests passed; the unrelated B2
  failures came from sandbox-blocked loopback binds and writes to the user state-cache
  directory.
- The G3 replay harness now defaults a missing fast-path march abrasion log value to zero,
  avoiding a `KeyError` on an iceless gate step. Ruff, `py_compile`, and `quality_check.sh`
  pass on the follow-up branch.
- `affected_test.sh --base origin/main -n 2`: 2,260 passed, 12 failed and 50 errors in
  725 s. The host temp volume reached 100% capacity (116 MiB free); the reported setup and
  upload errors include `OSError: [Errno 28] No space left on device`. The full affected run
  needs a clean-volume rerun before treating the remaining failures as code regressions.
- `stress_test.sh -n 2`: interrupted after 22 completed cases, all of which passed, after
  more than 90 minutes. Workers were still CPU-bound; the suite did not finish, so this is
  not a passing stress result. The host volume was at 100% capacity during the run.

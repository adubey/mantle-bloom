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

### Implementation checks

- `ruff check`: passed.
- `quality_check.sh`: passed (vulture, knip and jscpd; zero new duplicate clones).
- `git diff --check`: passed.
- Focused unit tests: `test_hydrology.py` 155 passed; `test_erosion.py` 64 passed.
- `affected_test.sh --base origin/main -n 2`: 2,260 passed, 12 failed and 50 errors in
  725 s. The host temp volume reached 100% capacity (116 MiB free); the reported setup and
  upload errors include `OSError: [Errno 28] No space left on device`. The full affected run
  needs a clean-volume rerun before treating the remaining failures as code regressions.
- `stress_test.sh -n 2`: interrupted after 22 completed cases, all of which passed, after
  more than 90 minutes. Workers were still CPU-bound; the suite did not finish, so this is
  not a passing stress result. The host volume was at 100% capacity during the run.

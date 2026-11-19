# Issue 973: RFC 0015 D3 implementation

Implements the accepted D3 direction in [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md#d3-glacial-erosion): integrate Cook et al. (2020)'s physical basal-sliding abrasion law through the D2 ice march, update a working bed between intervals, and transfer accumulated rock once through the existing terrain and material ledgers.

**Principles:** #5 (consistent units) and #8 (each process follows its own physics). The prior closing-ice-times-world-step estimate was then clipped with the shared neighbor-drop limit, creating a 66–73% gap between 0.25 Myr and 1 Myr steps on G3's diagnosis saves. The invariant is that each march interval contributes the physical law times interval duration, the next interval receives the bed lowered by the previous interval's abrasion, no held interval omits continuing abrasion, and the final material transfer remains within crustal headroom and closes the existing sediment ledger.

## Pre-registered replay

Before replaying, expect the integrated abrasion total at 0.25 Myr and 1 Myr to differ materially less than G3's prior 66–73% gap; a gap under 20% is the working target, while a remaining gap above 50% means D3 has not removed the dominant timestep sensitivity. Expect applied/raw to rise above G3's 0.0045–0.210 range where the removed material has available crustal headroom; any remaining loss should be explainable by headroom, mobile cover and craton resistance. Runtime is a measured cost rather than a pass threshold; compare each D3 arm with the D2 baseline under the same save, horizon and step schedule, and report slowdown alongside the interval count.

Controls: Cook et al. (2020), coefficient 1e-4 m/yr, exponent 0.65, glacier_erosion_multiplier = 1, ice ages disabled, same 10 Myr horizon and deterministic save for each paired run. Use each save at 0.25 Myr × 40 steps and 1 Myr × 10 steps. D2 baseline runs use origin/main; D3 runs use the issue branch. Do not tune coefficients after seeing these results.

| Save | Age | SHA-256 |
|---|---:|---|
| mantle-bloom-seed649658380-6400000y.bloom | 6.4 Myr | 06167cb4e45fac7141e98df2cf091999ffc40b36be4c0da3bd94703b6b31131f |
| mantle-bloom-seed997271774-87400000y.bloom | 87.4 Myr | 07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8 |

### Results

**Superseded by the post-review correction below.** These full-run numbers were captured
before the land/ocean mask in material transfer was corrected; the inverted mask zeroed
applied abrasion on land. The D3 rows therefore do not represent this implementation and
must not be used to close the G3 gate. The D2 control rows remain valid.

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

These are two deterministic diagnostic saves, not a broad stability survey. The G3 target is
partially met; the remaining 34.1% gap on the larger save and the runtime cost remain explicit
follow-up considerations for RFC closure.

### Post-review correction

Review found that the land/ocean mask in `erosion.apply_erosion` was inverted, so the transfer
step discarded land abrasion even though the ice march fed it back into the bed. Commit
`b94322d2` corrects the mask. A one-step 0.25 Myr smoke replay on seed649658380 then reported
262,419.27 km³ raw and applied abrasion (applied/raw 1.000; 42.9 s, peak RSS 0.98 GB). This
confirms the transfer now matches the march for that step, but it does not replace the four
full D3 replay arms. The prior D3 replay table is invalidated until those arms are rerun.

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

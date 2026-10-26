# RFC 0014 D2(a) review replay

## Preregistration

Compare the #519 mature save `seed308161711-121.5Myr` at 100 kyr steps for 10 steps (1 Myr), from the same initial bytes, on `origin/main` and PR #610 head. Track the #519 lake IDs 26 and 13 as well as global visible lake volume and area, using `bin/debug/measure_lake_split_reset.py`.

**Primary measures:** global visible lake volume and area after step 10. **Secondary measures:** the same measures for lakes 26 and 13, each step's median runtime, and peak RSS. This is a calibration-shift report, not a pass/fail gate: report the signed and relative differences. No lake-specific run-to-run noise estimate is available; the identical inputs and deterministic replay make the paired trajectory comparison useful, but it cannot establish variation across other seeds or resolutions. The expected global direction is uncertain; cold seasonal lakes are expected to lose more water where the new sublimation term exceeds the old annual-mean evaporation.

**Controls:** same save bytes and SHA, seed, step size, 10-step horizon, tracked lake IDs, replay script, and Python environment. Run each arm alone, with no tests or other replay running. Baseline is a task-specific detached worktree at `origin/main`; branch is the PR head.


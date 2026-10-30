# RFC 0014 D2(a) review replay

## Preregistration

Compare the #519 mature save `seed308161711-121.5Myr` at 100 kyr steps for 10 steps (1 Myr), from the same initial bytes, on `origin/main` and PR #610 head. Track the #519 lake IDs 26 and 13 as well as global visible lake volume and area, using `bin/debug/measure_lake_split_reset.py`.

**Primary measures:** global visible lake volume and area after step 10. **Secondary measures:** the same measures for lakes 26 and 13, each step's median runtime, and peak RSS. This is a calibration-shift report, not a pass/fail gate: report the signed and relative differences. No lake-specific run-to-run noise estimate is available; the identical inputs and deterministic replay make the paired trajectory comparison useful, but it cannot establish variation across other seeds or resolutions. The expected global direction is uncertain; cold seasonal lakes are expected to lose more water where the new sublimation term exceeds the old annual-mean evaporation.

**Controls:** same save bytes and SHA, seed, step size, 10-step horizon, tracked lake IDs, replay script, and Python environment. Run each arm alone, with no tests or other replay running. Baseline is a task-specific detached worktree at `origin/main`; branch is the PR head.

## Results

The save hash matched the #519 record. Baseline was `origin/main` at `1d1048dd354323257e86ab03628a1728b0b729a5`; the branch code was `c5fe0f4` (replay worktree HEAD `8ca6f72`, which adds only this preregistration after the code commit). Both runs used 128,340 nodes, the same 10 × 100 kyr steps, and `measure_lake_split_reset.py --track 26 13`.

| Measure at 1 Myr | `origin/main` | PR branch | Change |
|---|---:|---:|---:|
| Global visible lake volume | 1,931,177.9 km³ | 1,931,946.4 km³ | +768.5 km³ (+0.04%) |
| Global visible lake area | 5,215,762 km² | 5,488,950 km² | +273,188 km² (+5.24%) |
| Tracked lake 26 volume | 606.8 km³ | 701.4 km³ | +94.6 km³ (+15.60%) |
| Tracked lake 13 volume | 166.1 km³ | 165.1 km³ | −1.0 km³ (−0.61%) |
| Median step runtime | 13.76 s | 13.67 s | −0.67% |
| Peak RSS | 1,244.9 MB | 1,052.6 MB | −192.3 MB |

The global volume shift is small on this one trajectory; area is more sensitive and differs by 5.24%. Within-run global lake area also varied from 4.4 to 9.8 million km², so the endpoint delta is not evidence of a general climate calibration shift. The 10 paired step runtimes show no increase near the RFC's +10% target, but a single short replay is not a performance benchmark. Peak RSS is process high-water data and is reported without interpreting it as a model improvement.

Both arms emitted the existing `basin spills did not settle within 32 rounds` warning at step 10. Raw JSON outputs remain in scratch; the checked-in record is compact.

Commands (with local paths replaced):

```bash
python -I bin/debug/measure_lake_split_reset.py <local-save> --steps 10 --step-years 100000 --track 26 13 --backend <baseline>/backend --out <scratch-output>
python -I bin/debug/measure_lake_split_reset.py <local-save> --steps 10 --step-years 100000 --track 26 13 --backend <branch>/backend --out <scratch-output>
```

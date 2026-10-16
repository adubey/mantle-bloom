# H4 craton carry replay

The original seed 11, 195 Myr H4 capture was lost in a disk event. With the user's approval,
this comparison regenerates that merge from seed 11 on `origin/main`, then replays the same
pre-merge capture on `origin/main` and the issue branch. This is a fresh generated capture, not
the lost artifact.

## Method

- Generator: `world.generate_world(11)` with defaults, 5 Myr steps, through 200 Myr.
- Captured merge: 195 Myr, plates 4 and 38; the capture SHA-256 is
  `77aaed689e073bf7bcd032f4edb11b2660f1b2dbd5365a99d5cd60b6bf9c7d05`.
- Baseline: detached `origin/main`, `840c5ad3e996cf3ecc50f928eb967ad4a99fc5a0`.
- Branch: `26b1db9` (H4 fix).
- The generator took 110 seconds. Each capture replay ran alone: 1.32 seconds on baseline and
  1.31 seconds on the branch. Peak RSS was not recorded. No test or other replay overlapped
  these runs.

Commands, from the corresponding worktree roots:

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_clips.py \
  --seed 11 --to-myr 200 --step-years 5e6 \
  --capture-dir /private/tmp/mantle-bloom-issue443

PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_merge_placement.py \
  /private/tmp/mantle-bloom-issue443/seed11-195myr-4-38.pkl \
  --json /private/tmp/mantle-bloom-issue443/baseline.json

PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_merge_placement.py \
  /private/tmp/mantle-bloom-issue443/seed11-195myr-4-38.pkl \
  --json /private/tmp/mantle-bloom-issue443/branch-final.json
```

## Results

| Merge measure | `origin/main` | Issue branch |
|---|---:|---:|
| Craton above uncapped Hc after remap (M km³) | 62.774 | <0.000001 |
| Craton above Hc after placement (M km³) | 28.597 | 0 |
| Booked to `merge_capped_m3` (M km³) | 28.597 | 14.185 |
| Final craton volume (M km³) | 884.941 | 899.353 |
| Final craton area (km²) | 16,478,816 | 16,478,816 |
| Final continental / oceanic area (km²) | 33,114,141 / 13,048,487 | 33,114,141 / 13,048,487 |
| Overflow placement share | 1.000 | 1.000 |

The branch preserves the craton's local support and removes the above-column excess. The
remaining 14.185 M km³ has no Hc headroom on supported new cells and is booked as merge
capacity loss. Total merge loss for this event falls by 50.4%. The final craton area and the
merged continental/oceanic area totals are identical between the two variants for this event;
the new regression separately verifies that an oceanic-coded region without source craton
stays oceanic. The branch result is expected not to reproduce the baseline's booked volume;
the diagnostic's `reproduced` flag is false for that reason.

Compact values are in [replay.json](replay.json). Full captures and replay output remain in
scratch storage.

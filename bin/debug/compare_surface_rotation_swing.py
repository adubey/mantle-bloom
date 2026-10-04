#!/usr/bin/env python3
"""Does the quad surface (#250) change issue #177's rotation-rate land-loss swing?

Reads a run_sweep.py output holding `avg_rotation_rate` at 0.25x and 4x for both
`--surface quad` and `--surface lines` over the same seeds, and prints, per checkpoint and
outcome stat:

- each surface's 0.25x -> 4x means and Welch z (the same unpaired comparison every earlier
  #171/#177 comment reported, so the line-backed row is directly comparable with them);
- the surface x rotation-rate interaction: per seed, swing = stat(4x) - stat(0.25x), paired
  across surfaces on that seed, so a positive/negative z says quad's swing is larger/smaller
  than lines' on the same starting worlds.

Usage:
    backend/.venv/bin/python bin/debug/compare_surface_rotation_swing.py bin/debug/results/issue177_quad_vs_lines.jsonl
"""

from __future__ import annotations

import json
import math
import statistics
import sys
from collections import defaultdict

STATS = (
    ("land_fraction", "land%", 100.0),
    ("ice_cap_fraction", "ice-cap%", 100.0),
    ("plains_fraction", "plains%", 100.0),
    ("land_volume_above_sea_km3", "land vol (M km^3)", 1e-6),
)
LOW, HIGH = 0.25, 4.0


def welch_z(a: list[float], b: list[float]) -> float:
    se = math.sqrt(statistics.variance(a) / len(a) + statistics.variance(b) / len(b))
    return (statistics.fmean(b) - statistics.fmean(a)) / se if se else float("nan")


def paired_z(diffs: list[float]) -> float:
    se = statistics.stdev(diffs) / math.sqrt(len(diffs))
    return statistics.fmean(diffs) / se if se else float("nan")


def main() -> None:
    # (surface, multiplier, checkpoint) -> seed -> row
    rows: dict[tuple[str, float, int], dict[int, dict]] = defaultdict(dict)
    with open(sys.argv[1]) as f:
        for line in f:
            r = json.loads(line)
            if r["parameter"] != "avg_rotation_rate":
                continue
            rows[(r.get("surface", "lines"), r["multiplier"], int(r["checkpoint_years"]))][r["seed"]] = r

    checkpoints = sorted({cp for (_, _, cp) in rows})
    for cp in checkpoints:
        print(f"\n=== {cp / 1e6:.0f} My ===")
        print(f"{'stat':<18} {'surface':<6} {'n':>4} {'0.25x':>9} {'4x':>9} {'swing':>9} {'z':>7}")
        for key, label, scale in STATS:
            swings: dict[str, dict[int, float]] = {}
            for surface in ("lines", "quad"):
                low, high = rows.get((surface, LOW, cp), {}), rows.get((surface, HIGH, cp), {})
                seeds = sorted(low.keys() & high.keys())
                if len(seeds) < 2:
                    continue
                a = [low[s][key] * scale for s in seeds]
                b = [high[s][key] * scale for s in seeds]
                swings[surface] = {s: (high[s][key] - low[s][key]) * scale for s in seeds}
                print(
                    f"{label:<18} {surface:<6} {len(seeds):>4} {statistics.fmean(a):>9.2f} "
                    f"{statistics.fmean(b):>9.2f} {statistics.fmean(b) - statistics.fmean(a):>+9.2f} "
                    f"{welch_z(a, b):>+7.2f}"
                )
            if len(swings) == 2:
                common = sorted(swings["lines"].keys() & swings["quad"].keys())
                diffs = [swings["quad"][s] - swings["lines"][s] for s in common]
                if len(diffs) >= 2:
                    print(
                        f"{'':<18} {'q - l':<6} {len(diffs):>4} {'':>9} {'':>9} "
                        f"{statistics.fmean(diffs):>+9.2f} {paired_z(diffs):>+7.2f}  (swing difference, paired)"
                    )


if __name__ == "__main__":
    main()

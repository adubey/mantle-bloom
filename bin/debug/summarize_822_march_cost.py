#!/usr/bin/env python3
"""#822: tabulate `probe_822_march_cost.py` runs against each world's `d2` arm.

Reads `<dir>/<world>-<preset>.json` (and the `.npz` beside each), and prints one Markdown
row per (world, preset): march and step wall time as ratios to `d2` (summed over every
step), step-0 accuracy against `d2` from the same opening state (closing ice, node-wise
L1 depth difference over `d2`'s total, sea level), the worst ice and land ledger closure
over all steps, and the substep/cap/settle counts. `--json` also writes a compact copy.

    backend/.venv/bin/python bin/debug/summarize_822_march_cost.py <scratch-dir> \\
        --worlds gen5-250k gen5-1m gen30-250k s997 s649 --json analysis/issue822/compare.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

PRESETS = ["d2", "cap-hold", "cap-endpoint", "adaptive", "cheap-substep", "cheap-substep-limited"]


def _load(directory: Path, world: str, preset: str):
    path = directory / f"{world}-{preset}.json"
    if not path.exists():
        return None, None
    return json.loads(path.read_text()), np.load(path.with_suffix(".npz"))


def _row(world: str, preset: str, run, depth, ref, ref_depth) -> dict:
    steps, ref_steps = run["per_step"], ref["per_step"]
    march_s = sum(s["march_s"] for s in steps)
    wall_s = sum(s["wall_s"] for s in steps)
    s0, r0 = steps[0], ref_steps[0]
    d0, rd0 = depth["step0"], ref_depth["step0"]
    marches = [s["march"] or {} for s in steps]
    return {
        "world": world,
        "preset": preset,
        "march_s": round(march_s, 1),
        "march_ratio": round(march_s / max(sum(s["march_s"] for s in ref_steps), 1e-9), 3),
        "step_ratio": round(wall_s / max(sum(s["wall_s"] for s in ref_steps), 1e-9), 3),
        "ice0_rel": (s0["ice_closing_km3"] - r0["ice_closing_km3"]) / max(abs(r0["ice_closing_km3"]), 1e-9),
        "depth0_l1": float(np.abs(d0 - rd0).sum() / max(np.abs(rd0).sum(), 1e-9)),
        "sea0_m": s0["sea_level_m"] - r0["sea_level_m"],
        "ice_closure_max": max(s["ice_closure"] for s in steps),
        "land_closure_max": max(s["land_closure"] for s in steps),
        "substeps": [m.get("substeps") for m in marches],
        "settled": sum(bool(m.get("settled")) for m in marches),
        "capped": sum(bool(m.get("capped")) for m in marches),
        "max_substep_years": max((m.get("max_substep_years") or 0.0) for m in marches),
        "steps": len(steps),
        "peak_rss_gb": run["peak_rss_gb"],
    }


def main() -> None:
    args = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    args.add_argument("directory", type=Path)
    args.add_argument("--worlds", nargs="+", required=True)
    args.add_argument("--json", type=Path)
    ns = args.parse_args()

    rows = []
    for world in ns.worlds:
        ref, ref_depth = _load(ns.directory, world, "d2")
        if ref is None:
            print(f"(no d2 arm for {world})")
            continue
        for preset in PRESETS:
            run, depth = _load(ns.directory, world, preset)
            if run is not None:
                rows.append(_row(world, preset, run, depth, ref, ref_depth))

    print("| World | Preset | March s | March x | Step x | Ice0 vs d2 | Depth0 L1 | Sea0 m | Ice closure | Land closure | Substeps | Settled/capped |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(
            f"| {r['world']} | {r['preset']} | {r['march_s']:.1f} | {r['march_ratio']:.2f} | {r['step_ratio']:.2f} "
            f"| {100 * r['ice0_rel']:+.2f}% | {100 * r['depth0_l1']:.2f}% | {r['sea0_m']:+.3f} "
            f"| {r['ice_closure_max']:.1e} | {r['land_closure_max']:.1e} | {'/'.join(str(s) for s in r['substeps'])} "
            f"| {r['settled']}/{r['capped']} of {r['steps']} |"
        )
    if ns.json:
        ns.json.write_text(json.dumps(rows, indent=1) + "\n")
        print(f"wrote {ns.json}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Replay #403's runs and record the cratonic area inside the genuine set (issue #409).

`cratons.diagnostics` counts every cell with craton crust, including cratonic cells thinned
below CRATON_HOST_MIN_HC_M that have left the genuine set, so genuine minus craton area
is only a lower bound on non-cratonic host, and its error can move either way between
checkpoints. This wraps `compare_cratons._snapshot` (read-only) to add
`craton_genuine_area_km2` at each checkpoint, and runs `compare_cratons.run_job` exactly as
`rebaseline_rfcs.py` does, so the genuine and total craton areas can be checked against
#403's JSON for a matching trajectory.

`--model-root` names the checkout whose model and `compare_cratons.py` to run (here a
detached `2562f63` worktree, with its own venv). Both are imported from there, never from
this script's tree: `compare_cratons` puts its own tree's backend first on `sys.path`, so a
copy from a newer tree would silently run the newer model. Seeds run one at a time.

Usage:
    $BASE/backend/.venv/bin/python bin/debug/craton_genuine_area.py --model-root $BASE \\
        --seeds 7,11,2024 --years 1e9 --out /tmp/genuine.jsonl
    $BASE/backend/.venv/bin/python bin/debug/craton_genuine_area.py --model-root $BASE \\
        --seeds 456862899 --years 3e9 --out /tmp/genuine.jsonl
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def _craton_genuine_area_km2(world, cratons) -> float:
    area = 0.0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        inside = cratons._genuine_continental(plate) & (plate.collect("craton_crust_m") > 0.0)
        area += float(cratons._areas(world, plate)[inside].sum())
    return area / 1e6


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model-root", type=Path, required=True, help="checkout whose model to run")
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--every", type=float, default=100e6)
    parser.add_argument("--step", type=float, default=5e6)
    parser.add_argument("--out", type=Path, required=True, help="JSON lines, appended one run at a time")
    args = parser.parse_args()
    root = args.model_root.resolve()
    sys.path[:0] = [str(root / "bin/debug"), str(root / "backend")]
    import compare_cratons
    from app import cratons

    if not Path(cratons.__file__).resolve().is_relative_to(root):
        raise SystemExit(f"imported {cratons.__file__}, not the model under {root}")
    commit = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    snapshot = compare_cratons._snapshot

    def tapped(world, events_seen):
        return {**snapshot(world, events_seen), "craton_genuine_area_km2": _craton_genuine_area_km2(world, cratons)}

    compare_cratons._snapshot = tapped
    for seed in (int(s) for s in args.seeds.split(",")):
        rows = compare_cratons.run_job(seed, "on", args.years, args.every, args.step, {})
        keep = ("myr", "genuine_continental_area_km2", "craton_area_km2", "craton_genuine_area_km2", "craton_volume_km3")
        with args.out.open("a") as out:
            out.write(json.dumps({"seed": seed, "model_commit": commit, "checkpoints": [{k: row[k] for k in keep} for row in rows]}) + "\n")
        print(f"seed {seed}: final genuine-restricted share {rows[-1]['craton_genuine_area_km2'] / rows[-1]['genuine_continental_area_km2']:.3f}", flush=True)


if __name__ == "__main__":
    main()

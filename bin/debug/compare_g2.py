#!/usr/bin/env python3
"""Compare two RFC 0009 gate G2 probe results (#447): origin/main against a branch.

Reads two `probe_mixed_plates.py --out` files (PR #446) and writes a compact JSON with, per
seed and run (control, mixed A, mixed B), the measures #447 asks for at the last checkpoint:
genuine continental area (% of sphere), juvenile additions and deep subduction (M km^3),
active plates, merges, splits and the worst ledger balance error. Per seed and model it adds
the control-versus-mixed gaps (mean of the mixed runs minus the control) in genuine area at
the end and over the first checkpoint to the end (100->500 Myr by default), and in juvenile
additions and deep subduction, with G2's noise floor: max(|mixed A - mixed B|, floor), the
floor 1 point for area and 0 for the ledger accounts.

Usage:
    python3 bin/debug/compare_g2.py MAIN.json BRANCH.json --out analysis/issue447/g2_ab.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

M_KM3 = 1e6  # the probe's ledger is in km^3


def _end(run: dict) -> dict:
    first, last = run["checkpoints"][1], run["checkpoints"][-1]
    ledger = last["continental_ledger_km3"]
    return {
        "genuine_pct": last["genuine_pct"],
        "genuine_change_pct": last["genuine_pct"] - first["genuine_pct"],
        "juvenile_mkm3": ledger.get("juvenile_additions_m3", 0.0) / M_KM3,
        "deeply_subducted_mkm3": ledger.get("deeply_subducted_m3", 0.0) / M_KM3,
        "active_plates": last["active_plates"],
        "merges": last["merges"],
        "splits": last["splits"],
        "mixed_plate_share": last["mixed_plate_share"],
        "max_ledger_error_rel": max(
            max(abs(c["continental_balance_error_rel"]), abs(c["craton_balance_error_rel"])) for c in run["checkpoints"]
        ),
    }


def _gaps(runs: dict[str, dict]) -> dict:
    out = {}
    for key, floor in (("genuine_pct", 1.0), ("genuine_change_pct", 1.0), ("juvenile_mkm3", 0.0), ("deeply_subducted_mkm3", 0.0)):
        control, a, b = (runs[r][key] for r in ("control", "mixed_a", "mixed_b"))
        out[key] = {"gap": (a + b) / 2 - control, "noise": max(abs(a - b), floor)}
    return out


def summarize(path: Path) -> dict:
    result = json.loads(path.read_text())
    seeds: dict[int, dict] = {}
    for run in result["runs"]:
        seeds.setdefault(run["seed"], {})[run["run"]] = _end(run)
    return {
        "model_commit": result["model_commit"],
        "myr": result["years"] / 1e6,
        "seeds": {seed: {"runs": runs, "gaps": _gaps(runs)} for seed, runs in sorted(seeds.items())},
        "g2_failing": result["summary"]["failing"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("main", type=Path)
    parser.add_argument("branch", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = {"main": summarize(args.main), "branch": summarize(args.branch)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1, default=float) + "\n")
    for name, model in out.items():
        print(name, model["model_commit"][:7], "G2 failing:", model["g2_failing"])
        for seed, data in model["seeds"].items():
            gaps = " ".join(f"{k}={v['gap']:+.1f}(noise {v['noise']:.1f})" for k, v in data["gaps"].items())
            print(f"  seed {seed}: {gaps}")


if __name__ == "__main__":
    main()

"""Pooled scoring of the stage 1 (or stage 2) runs of RFC 0018's D1 search (#888), as predeclared on the issue.

Usage (repo root, this worktree's venv):
    backend/.venv/bin/python bin/debug/score_stage1.py --runs ~/mantle-bloom-captures/issue888/stage1 \
        --baseline ~/mantle-bloom-captures/issue888/production-baseline [--out analysis/issue888/stage1-scores.json]

For each ``setting*`` directory, pools the ``recorder.pkl`` of its worlds, scores the realism targets with
``rift_damage_recorder.score`` against the production splitter's rate on the same worlds (from the baseline
directory), and ranks the settings: most targets met, then the smaller total violation, then more splits. Settings
with a missing or failed run, a ledger gap above 1e-9 or no split are listed as excluded. Targets with fewer than
``--min-events`` pooled events count as not met.
"""

import argparse
import json
import pickle
import sys
from pathlib import Path


sys.path[:0] = [str(Path(__file__).resolve().parent)]

import rift_damage_mechanism as mech  # noqa: E402
import rift_damage_recorder as rec  # noqa: E402

SCORED = ("onset_width_km", "width_growth", "duration_myr", "transition_ratio", "inherited_weakness", "failed_rifts",
          "localization", "no_microplates", "line_shape", "splits_happen")
RANGES = {"onset_width_km": (100.0, 150.0), "duration_myr": (10.0, 50.0), "failed_rifts": (0.2, 0.8)}


def violation(name: str, entry: dict) -> float:
    """How far a target is from passing, as a fraction of its range (0 when met)."""
    if entry.get("ok"):
        return 0.0
    value = entry.get("value")
    if name in RANGES and value is not None:
        lo, hi = RANGES[name]
        return max(lo - value, value - hi, 0.0) / (hi - lo)
    if name == "transition_ratio" and value is not None:
        return max(3.0 - value, 0.0) / 3.0
    if name == "inherited_weakness" and value is not None:
        return max(2 / 3 - value, 0.0) / (2 / 3)
    if name == "line_shape" and entry.get("median_length_ratio") is not None:
        return max(1.1 - entry["median_length_ratio"], 0.0) / 1.1
    if name == "no_microplates" and entry.get("min_small_share") is not None:
        return max(0.15 - entry["min_small_share"], 0.0) / 0.15
    if name == "splits_happen" and entry.get("rate_per_world_100myr") is not None and entry.get("production_rate"):
        return max(0.5 * entry["production_rate"] - entry["rate_per_world_100myr"], 0.0) / max(0.5 * entry["production_rate"], 1e-9)
    return 1.0


def events_behind(name: str, result: dict) -> int:
    stats = result["candidates"]
    if name in ("onset_width_km", "width_growth", "transition_ratio"):
        return stats["successful"] if name != "onset_width_km" else stats["counted"]
    if name == "duration_myr":
        return stats["successful"]
    if name == "failed_rifts":
        return stats["counted"]
    return result["splits"]["value"] if name in ("inherited_weakness", "no_microplates", "line_shape", "localization") else 10**9


def pool(directory: Path, production_rate: float | None, min_events: int, drop_fragments: bool = True) -> dict:
    cands, events, super_max, super_viol, myr, ledger = [], [], 0.0, 0, 0.0, 0.0
    errors = []
    for run in sorted(p for p in directory.iterdir() if p.is_dir()):
        pickled, summary = run / "recorder.pkl", run / "run.json"
        if not (pickled.exists() and summary.exists()):
            errors.append(f"{run.name}: no result")
            continue
        data = pickle.load(open(pickled, "rb"))
        offset = len(cands)
        for c in data["candidates"]:
            c.ident += offset
            cands.append(c)
        events += data["events"]
        super_max, super_viol, myr = max(super_max, data["super_sea_max_m"]), super_viol + data["super_sea_violations"], myr + data["years_myr"]
        ledger = max(ledger, max(json.load(open(summary))["ledger_relative_gap_max"].values(), default=0.0))
    result = rec.score(cands, events, 200e6, super_max, super_viol, mech.Params().a_min, production_rate, myr, drop_fragments)
    for name in SCORED:
        if events_behind(name, result) < min_events:
            result[name]["ok"], result[name]["not_evaluable"] = False, True
    result["targets_passed"] = int(sum(bool(result[n]["ok"]) for n in SCORED))
    result["violation_total"] = float(sum(violation(n, result[n]) for n in SCORED))
    result["ledger_gap"], result["errors"] = ledger, errors
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--min-events", type=int, default=5)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--keep-fragments", action="store_true", help="the original rule: count split children seen once")
    args = parser.parse_args()
    seeds = sorted(p.name for p in next(args.runs.glob("setting*")).iterdir() if p.is_dir())
    base = [json.load(open(args.baseline / f"{s}.json")) for s in seeds]
    rate = sum(b["n_splits"] for b in base) / len(base) / base[0]["years"] * 100.0
    ranked, excluded = [], []
    for directory in sorted(args.runs.glob("setting*")):
        result = pool(directory, rate, args.min_events, not args.keep_fragments)
        row = {"setting": directory.name, "targets_passed": result["targets_passed"], "violation_total": result["violation_total"],
               "splits": result["splits"]["value"], "ledger_gap": result["ledger_gap"], "errors": result["errors"],
               "met": [n for n in SCORED if result[n]["ok"]], "not_evaluable": [n for n in SCORED if result[n].get("not_evaluable")]}
        if result["errors"] or result["ledger_gap"] > 1e-9 or not row["splits"]:
            excluded.append(row)
        else:
            ranked.append(row)
    ranked.sort(key=lambda r: (-r["targets_passed"], r["violation_total"], -r["splits"]))
    out = {"production_rate_per_world_100myr": rate, "ranked": ranked, "excluded": excluded}
    print(json.dumps(out, indent=1))
    if args.out:
        args.out.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

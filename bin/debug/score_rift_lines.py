"""Pooled realism-target scores of the rift-line arms of a `run_rift_lines_arms.py` queue (RFC 0018 D1, #888).

Usage (repo root, this worktree's venv):
    backend/.venv/bin/python bin/debug/score_rift_lines.py --runs ~/mantle-bloom-captures/issue888/rift-lines-baseline \
        [--baseline-arm great-circle] [--min-events 5] [--keep-fragments] [--out analysis/issue888/rift-lines-baseline/scores.json]

The great-circle arm gives the production splitter's rate in splits per world per 100 Myr, over the same seeds,
which the "splits happen" target compares each rift-line arm with. Every other arm directory is pooled over its worlds
by `rift_line_recorder.pool` and scored on the ten targets of the #888 predeclaration, then ranked: most targets met,
then the smaller total violation, then more splits. An arm with a missing or failed run, a ledger gap above 1e-9 or
no split is listed as excluded. Each row also has the numbers behind the D1 strength-arm table.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path[:0] = [str(Path(__file__).resolve().parent), str(Path(__file__).resolve().parents[2] / "backend")]

import rift_line_recorder as rec  # noqa: E402
from app.rift_strength import RiftLineParams  # noqa: E402


def baseline_rate(directory: Path, seeds: list[str]) -> float:
    runs = [json.loads((directory / seed / "run.json").read_text()) for seed in seeds]
    return sum(r["n_splits"] for r in runs) / sum(r["years_myr"] for r in runs) * 100.0


def row(name: str, result: dict) -> dict:
    """An arm's summary: the ranking keys and the D1 strength-arm table's lines."""
    weak = result["inherited_weakness"]
    return {
        "arm": name, "targets_passed": result["targets_passed"], "violation_total": result["violation_total"],
        "splits": result["splits"]["value"], "splits_per_world_100myr": result["splits_happen"]["rate_per_world_100myr"],
        "continent_majority_splits": None if result["oceanic_share_of_successful"]["value"] is None
        else round((1.0 - result["oceanic_share_of_successful"]["value"]) * result["splits"]["value"]),
        "splits_mostly_weak": weak["value"], "weak_share_lines": weak.get("mean_path_weak_share"), "weak_share_random": weak.get("random_path_weak_share"),
        "completed_median_myr": result["duration_myr"]["value"], "failed_share": result["failed_rifts"]["value"],
        "transition_ratio": result["transition_ratio"]["value"], "localization_violations": result["localization"]["violations"],
        "onset_width_km": result["onset_width_km"]["value"], "ledger_gap": result["ledger_gap"], "errors": result["errors"],
        "met": [n for n in rec.SCORED if result[n]["ok"]], "not_evaluable": [n for n in rec.SCORED if result[n].get("not_evaluable")],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--baseline-arm", default="great-circle")
    parser.add_argument("--min-events", type=int, default=5)
    parser.add_argument("--keep-fragments", action="store_true", help="the original rule: count split children seen once")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    base_dir = args.runs / args.baseline_arm
    seeds = sorted(p.name for p in base_dir.iterdir() if p.is_dir())
    rate = baseline_rate(base_dir, seeds)
    a_min = RiftLineParams().a_min
    ranked, excluded = [], []
    for directory in sorted(p for p in args.runs.iterdir() if p.is_dir() and p.name != args.baseline_arm):
        result = rec.pool(directory, rate, a_min, args.min_events, not args.keep_fragments)
        summary = row(directory.name, result)
        missing = sorted(set(seeds) - {p.name for p in directory.iterdir() if p.is_dir()})
        summary["errors"] += [f"{seed}: not run" for seed in missing]
        (excluded if summary["errors"] or summary["ledger_gap"] > 1e-9 or not summary["splits"] else ranked).append(summary)
    ranked.sort(key=lambda r: (-r["targets_passed"], r["violation_total"], -r["splits"]))
    out = {"seeds": [int(s) for s in seeds], "production_rate_per_world_100myr": rate, "ranked": ranked, "excluded": excluded}
    text = json.dumps(out, indent=1) + "\n"
    print(text, end="")
    if args.out:
        args.out.write_text(text)


if __name__ == "__main__":
    main()

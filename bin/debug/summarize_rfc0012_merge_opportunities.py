#!/usr/bin/env python3
"""Summarize matched merge-opportunity JSONL traces into a compact comparison JSON."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import statistics


def _read_run(path: Path) -> dict:
    steps = []
    final = None
    for line in path.read_text().splitlines():
        row = json.loads(line)
        if row.get("kind") == "summary":
            final = row
        elif row.get("kind") == "step":
            steps.append(row)
    if final is None or not steps:
        raise ValueError(f"incomplete trace: {path}")

    stages = Counter()
    decision_counts = Counter()
    source_counts = Counter()
    merge_paths = Counter()
    reset_count = 0
    overlap_pair_frames = 0
    detail_values = defaultdict(lambda: defaultdict(list))
    deferred_step_count = 0
    collision_ready_deferred = 0

    for row in steps:
        search = row["collision_search"]
        stages.update(search.get("stages", {}))
        for detail in search.get("pair_details", []):
            samples = detail_values[detail["stage"]]
            for source, target in (
                ("close_contact_points", "close_contact_points"),
                ("convergent_contacts", "convergent_contacts"),
                ("relative_rotation_cm_per_yr", "relative_rotation_cm_per_yr"),
            ):
                if detail.get(source) is not None:
                    samples[target].append(detail[source])
            closing = detail.get("closing_cm_per_yr")
            if closing:
                samples["max_contact_closing_cm_per_yr"].append(closing["max"])
        progress = row["collision_progress"]
        for item in progress.get("pairs", []):
            decision_counts["tracked_pair_steps"] += 1
            source_counts[item["source"]] += 1
            if item["timer_ready"]:
                decision_counts["timer_ready_steps"] += 1
            if item["probability_pass"] is True:
                decision_counts["probability_pass_steps"] += 1
            elif item["probability_pass"] is False:
                decision_counts["probability_fail_steps"] += 1
        reset_count += len(progress.get("expired_pairs", []))
        decision_counts["returned_collision_ready_pairs"] += len(progress.get("ready_pairs", []))
        decision_counts["ready_pairs_deferred_by_one_merge_limit"] += len(progress.get("other_ready_deferred", []))
        overlap_pair_frames += len(row["overlap_progress_update"].get("pairs", []))
        if row["forced_merge"].get("deferred_by_collision_merge"):
            collision_ready_deferred += 1
        if row["forced_merge"].get("selected_pair"):
            decision_counts["forced_overlap_selected"] += 1
        actual = row.get("actual_merge")
        if actual:
            merge_paths[actual["path"]] += 1
        if row["forced_merge"].get("ready_pairs"):
            deferred_step_count += 1

    return {
        "model": path.name.split("-seed")[0],
        "seed": final["seed"],
        "commit": final["commit"],
        "steps": final["steps"],
        "step_myr": final["step_myr"],
        "runtime_seconds": final["runtime_seconds"],
        "pair_filter_stage_totals": dict(stages),
        "pair_detail_statistics": {
            stage: {
                name: {
                    "count": len(values),
                    "median": round(statistics.median(values), 4),
                    "max": round(max(values), 4),
                }
                for name, values in sorted(fields.items())
                if values
            }
            for stage, fields in sorted(detail_values.items())
        },
        "collision_timer_and_rolls": dict(decision_counts),
        "tracked_pair_step_sources": dict(source_counts),
        "collision_progress_reset_count": reset_count,
        "overlap_pair_step_records": overlap_pair_frames,
        "steps_with_ready_forced_overlap": deferred_step_count,
        "steps_forced_merge_deferred_by_collision_merge": collision_ready_deferred,
        "merge_events_by_path": dict(merge_paths),
        "collision_ready_decision_mismatches": final["counters"].get(
            "collision_ready_decision_mismatches", 0
        ),
        "final_collision_progress_pairs": final["final_collision_progress_pairs"],
        "final_overlap_progress_pairs": final["final_overlap_progress_pairs"],
    }


def summarize(directory: Path, output: Path) -> dict:
    paths = sorted(directory.glob("baseline-seed*.jsonl")) + sorted(directory.glob("pilot-seed*.jsonl"))
    runs = [_read_run(path) for path in paths]
    by_model_seed = {(row["model"], row["seed"]): row for row in runs}
    seeds = sorted({row["seed"] for row in runs})
    matched = {}
    for seed in seeds:
        baseline = by_model_seed.get(("baseline", seed))
        pilot = by_model_seed.get(("pilot", seed))
        if baseline is None or pilot is None:
            continue
        stage_keys = set(baseline["pair_filter_stage_totals"]) | set(pilot["pair_filter_stage_totals"])
        matched[str(seed)] = {
            "baseline_commit": baseline["commit"],
            "pilot_commit": pilot["commit"],
            "baseline_pair_filter_stage_totals": baseline["pair_filter_stage_totals"],
            "pilot_pair_filter_stage_totals": pilot["pair_filter_stage_totals"],
            "pilot_minus_baseline_pair_filter_stage_totals": {
                key: pilot["pair_filter_stage_totals"].get(key, 0) - baseline["pair_filter_stage_totals"].get(key, 0)
                for key in sorted(stage_keys)
            },
            "baseline_collision_timer_and_rolls": baseline["collision_timer_and_rolls"],
            "pilot_collision_timer_and_rolls": pilot["collision_timer_and_rolls"],
            "baseline_merge_events_by_path": baseline["merge_events_by_path"],
            "pilot_merge_events_by_path": pilot["merge_events_by_path"],
        }
    result = {
        "diagnostic": "RFC 0012 matched merge-opportunity tracing",
        "runs": runs,
        "matched": matched,
        "interpretation_note": (
            "Filter-stage totals count continental pair checks over all steps. A pair can appear in multiple stages only once per step; "
            "the pair details give its deepest reached stage. Bounding-sphere rejections are aggregated, while pairs that survive that "
            "prune retain their IDs and contact-rate details. Collision timers and probability rolls are reconstructed with the same "
            "deterministic inputs and cross-checked against the production function's returned ready list."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace_dir", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.trace_dir, args.out)
    print(json.dumps({"runs": len(result["runs"]), "matched_seeds": len(result["matched"]), "out": str(args.out)}))


if __name__ == "__main__":
    main()

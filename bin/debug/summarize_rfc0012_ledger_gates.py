#!/usr/bin/env python3
"""Summarize matched JSONL traces from trace_rfc0012_ledger_gates.py.

Example:

    backend/.venv/bin/python bin/debug/summarize_rfc0012_ledger_gates.py \
      /private/tmp/rfc0012-ledger-gates --out analysis/rfc0012/g2_ledger_gate_traces.json

The summary retains only matched endpoint ledgers, event counts, physical trace aggregates,
and event-to-ledger closure checks. Full JSONL traces and generated worlds stay out of Git.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path


MKM3 = 1.0e15


def _read(path: Path) -> tuple[list[dict], dict]:
    rows = [json.loads(line) for line in path.open()]
    return [row for row in rows if row["kind"] == "step"], rows[-1]


def _sum_ledger(events: list[dict], event_kind: str, ledger: str, account: str) -> float:
    return sum(
        float(event.get("ledger_delta", {}).get(ledger, {}).get(account, 0.0))
        for event in events if event["kind"] == event_kind
    )


def _run_summary(path: Path) -> dict:
    steps, final = _read(path)
    events = [event for step in steps for event in step.get("events", [])]
    counts = Counter(event["kind"] for event in events)
    continental = final["ledgers"]["continental_material_ledger"]
    craton = final["ledgers"]["craton_ledger"]
    quiet = [event for event in events if event["kind"] == "quiet_gate"]
    deformation = [event for event in events if event["kind"] == "local_lithosphere_deformation"]
    subduction = [event for event in events if event["kind"] == "subduction_front"]
    merges = [event for event in events if event["kind"] == "plate_merge"]
    sutures = [event for event in events if event["kind"] == "suture_transfer"]
    formed = sum(float(event.get("formed_this_step_m3", 0.0)) for event in events if event["kind"] == "formation_total")
    deep_subducted = _sum_ledger(events, "subduction_front", "continental_material_ledger", "deeply_subducted_m3")
    merge_capped = sum(float(event.get("merge_capped_delta_m3", 0.0)) for event in merges)
    sub_craton = _sum_ledger(events, "subduction_front", "craton_ledger", "subducted_m3")
    local_nonzero = sum(
        event["quiet_host_cells"] > 0 and (event["quiet_host_hc_delta_m3"] != 0.0 or event["quiet_host_hm_delta_m3"] != 0.0)
        for event in deformation
    )
    reason_area_steps: dict[str, float] = {}
    for event in quiet:
        for reason, values in event["reset_reason"].items():
            reason_area_steps[reason] = reason_area_steps.get(reason, 0.0) + float(values["area_m2"])
    cap_area = [
        float(event.get("erosion_unloading", {}).get("ice_unloading_area_m2", 0.0))
        for event in quiet if event.get("erosion_unloading")
    ]
    max_abs_rebound = max((
        max(
            abs(float(event.get("erosion_unloading", {}).get("ice_deflection_change_m_min", 0.0))),
            abs(float(event.get("erosion_unloading", {}).get("ice_deflection_change_m_max", 0.0))),
        )
        for event in quiet if event.get("erosion_unloading")
    ), default=0.0)
    return {
        "model": path.name.split("-seed", 1)[0],
        "seed": final["seed"],
        "commit": final["commit"],
        "steps": final["steps"],
        "step_myr": final["step_myr"],
        "runtime_min": final["runtime_seconds"] / 60.0,
        "endpoint": final["final_geology"],
        "event_counts": dict(sorted(counts.items())),
        "ledger_mkm3": {
            "deeply_subducted": continental.get("deeply_subducted_m3", 0.0) / MKM3,
            "collision_subducted": continental.get("collision_subducted_m3", 0.0) / MKM3,
            "formed_craton": craton.get("formed_m3", 0.0) / MKM3,
            "subducted_craton": craton.get("subducted_m3", 0.0) / MKM3,
            "merge_capped_craton": craton.get("merge_capped_m3", 0.0) / MKM3,
        },
        "reconciliation_m3": {
            "deep_subduction_event_minus_endpoint": deep_subducted - continental.get("deeply_subducted_m3", 0.0),
            "craton_subduction_event_minus_endpoint": sub_craton - craton.get("subducted_m3", 0.0),
            "merge_cap_event_minus_endpoint": merge_capped - craton.get("merge_capped_m3", 0.0),
            "formation_event_minus_endpoint": formed - craton.get("formed_m3", 0.0),
            "craton_balance_error": final["craton_balance_error_m3"],
            "continental_material_balance_error": final["continental_material_balance_error_m3"],
        },
        "subduction_trace": {
            "front_events": len(subduction),
            "candidate_cells": sum(event["candidate_cells"] for event in subduction),
            "removed_cells": sum(event["removed_cells"] for event in subduction),
            "removed_deep_subduction_cells": sum(event["removed_deep_subduction_cells"] for event in subduction),
            "eligible_donor_cells_not_removed": sum(event["eligible_donor_cells"] - event["removed_donor_cells"] for event in subduction),
            "budget_blocked_cells": sum(event["budget_blocked_cells"] for event in subduction),
            "suture_transfer_events": len(sutures),
            "suture_unplaced_hc_mkm3": sum(event["unplaced_m3"] for event in sutures) / MKM3,
            "suture_lost_lower_crust_mkm3": sum(event["lost_m3"] for event in sutures) / MKM3,
        },
        "merge_trace": {
            "merge_events": len(merges),
            "merge_pairs": [[event["keep_plate_id"], event["absorbed_plate_id"]] for event in merges],
            "merge_capped_event_mkm3": merge_capped / MKM3,
            "premerge_receiver_headroom_mkm3": sum(
                event["pair_before"][0]["merge_hc_headroom_m3"] for event in merges
            ) / MKM3,
        },
        "quiet_gate_trace": {
            "plate_step_samples": len(quiet),
            "formation_events": counts.get("formation_total", 0),
            "quiet_host_area_step_mkm2": sum(event["quiet_area_m2"] for event in quiet) / 1e12,
            "eligible_uncratonized_area_step_mkm2": sum(event["eligible_uncratonized_area_m2"] for event in quiet) / 1e12,
            "reset_reason_area_step_mkm2": {k: v / 1e12 for k, v in sorted(reason_area_steps.items())},
            "quiet_host_plate_steps_with_local_hc_or_hm_change": local_nonzero,
            "quiet_host_abs_hc_change_mkm3": sum(abs(event["quiet_host_hc_delta_m3"]) for event in deformation) / MKM3,
            "quiet_host_abs_hm_change_mkm3": sum(abs(event["quiet_host_hm_delta_m3"]) for event in deformation) / MKM3,
        },
        "erosion_unloading_trace": {
            "plate_step_samples": len(cap_area),
            "sum_ice_unloaded_area_step_mkm2": sum(cap_area) / 1e12,
            "max_plate_step_ice_unloaded_area_km2": max(cap_area, default=0.0) / 1e6,
            "max_abs_ice_deflection_change_m": max_abs_rebound,
        },
        "step_balance_max_abs_m3": {
            "craton": max(abs(float(step["craton_balance_error_m3"])) for step in steps),
            "continental_material": max(abs(float(step["continental_material_balance_error_m3"])) for step in steps),
        },
    }


def summarize(root: Path) -> dict:
    manifest = json.loads((root / "manifest.json").read_text())
    runs = [_run_summary(path) for path in sorted(root.glob("*-seed*-1myr-500myr.jsonl"))]
    by_pair = {}
    for seed in manifest["seeds"]:
        pair = {run["model"]: run for run in runs if run["seed"] == seed}
        if set(pair) == {"baseline", "pilot"}:
            keys = ("deeply_subducted", "subducted_craton", "merge_capped_craton", "formed_craton")
            by_pair[str(seed)] = {
                "baseline_merge_count": pair["baseline"]["merge_trace"]["merge_events"],
                "pilot_merge_count": pair["pilot"]["merge_trace"]["merge_events"],
                "baseline_land_fraction": pair["baseline"]["endpoint"]["land_fraction"],
                "pilot_land_fraction": pair["pilot"]["endpoint"]["land_fraction"],
                "pilot_minus_baseline_ledger_mkm3": {
                    key: pair["pilot"]["ledger_mkm3"][key] - pair["baseline"]["ledger_mkm3"][key]
                    for key in keys
                },
            }
    return {"manifest": manifest, "runs": runs, "matched": by_pair}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace_dir", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summarize(args.trace_dir), indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

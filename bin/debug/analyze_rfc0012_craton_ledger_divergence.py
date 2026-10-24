#!/usr/bin/env python3
"""Locate and decompose matched G2 craton-ledger divergence from JSONL traces.

The input traces are produced by ``trace_rfc0012_ledger_gates.py``. This reads completed
baseline/pilot runs, reconstructs live craton volume from initial + formed - booked sinks,
and writes 25 Myr checkpoints plus a final per-account attribution.

    backend/.venv/bin/python bin/debug/analyze_rfc0012_craton_ledger_divergence.py \
      /private/tmp/rfc0012-ledger-gates \
      --out analysis/rfc0012/g2_craton_ledger_divergence.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

MKM3 = 1.0e15
CRATON_SOURCES = ("formed_m3",)
CRATON_SINKS = (
    "rifted_m3",
    "subducted_m3",
    "delaminated_m3",
    "merge_capped_m3",
    "collision_reworked_m3",
    "eroded_m3",
    "topology_removed_m3",
    "audit_clipped_m3",
    "unattributed_m3",
)
CRATON_RIFT_PATHS = ("rifted_divergent_m3", "rifted_stretch_m3", "rifted_failed_m3")
CHECKPOINT_MYR = 25


def _read_run(path: Path) -> dict:
    steps = []
    final = None
    for line in path.open():
        row = json.loads(line)
        if row["kind"] == "step":
            steps.append(row)
        elif row["kind"] == "summary":
            final = row
    if final is None or len(steps) != final["steps"]:
        raise ValueError(f"incomplete trace: {path}")

    final_ledger = final["ledgers"]["craton_ledger"]
    initial = float(final_ledger.get("initial_m3", 0.0))
    cumulative = {name: 0.0 for name in CRATON_SOURCES + CRATON_SINKS + CRATON_RIFT_PATHS}
    checkpoints = [{"elapsed_myr": 0.0, "initial_mkm3": initial / MKM3, "formed_mkm3": 0.0,
                    "sinks_mkm3": {name: 0.0 for name in CRATON_SINKS},
                    "rift_paths_mkm3": {name: 0.0 for name in CRATON_RIFT_PATHS},
                    "live_mkm3": initial / MKM3}]
    live_series = [(0.0, initial)]
    formation_onset = None
    merge_events = []
    merge_event_cap_m3 = 0.0
    merge_event_cap_count = 0
    merge_live_cap_residual_m3 = 0.0
    rift_activity = {
        name: {"nonzero_steps": 0, "events_m3": 0.0}
        for name in CRATON_RIFT_PATHS
    }
    failed_rift = {
        "cuts": 0,
        "cuts_through_craton": 0,
        "band_area_m2": 0.0,
        "band_craton_area_m2": 0.0,
        "cut_steps": [],
    }
    rift_event_paths = {
        "divergent_thinning": {"ledger_account": "rifted_divergent_m3", "events": 0, "volume_m3": 0.0,
                               "cellwise_volume_m3": 0.0, "affected_cell_visits": 0, "by_plate": {}},
        "stretch_donation": {"ledger_account": "rifted_stretch_m3", "events": 0, "volume_m3": 0.0,
                             "cellwise_volume_m3": 0.0, "affected_cell_visits": 0, "by_plate": {}},
        "failed_rift_thinning": {"ledger_account": "rifted_failed_m3", "events": 0, "volume_m3": 0.0,
                                 "cellwise_volume_m3": 0.0, "affected_cell_visits": 0, "by_plate": {}},
    }
    for step in steps:
        delta = step["ledger_delta"].get("craton_ledger", {})
        for account in cumulative:
            cumulative[account] += float(delta.get(account, 0.0))
        for account in CRATON_RIFT_PATHS:
            amount = float(delta.get(account, 0.0))
            rift_activity[account]["events_m3"] += amount
            rift_activity[account]["nonzero_steps"] += int(amount > 0.0)
        failed_cuts = int(delta.get("failed_rift_cuts", 0.0))
        failed_rift["cuts"] += failed_cuts
        failed_rift["cuts_through_craton"] += int(delta.get("failed_rift_cuts_through_craton", 0.0))
        failed_rift["band_area_m2"] += float(delta.get("failed_rift_band_area_m2", 0.0))
        failed_rift["band_craton_area_m2"] += float(delta.get("failed_rift_band_craton_area_m2", 0.0))
        if failed_cuts:
            failed_rift["cut_steps"].append({
                "step": int(step["step"]),
                "elapsed_myr": float(step["elapsed_myr"]),
                "cuts": failed_cuts,
                "cuts_through_craton": int(delta.get("failed_rift_cuts_through_craton", 0.0)),
                "rifted_failed_m3": float(delta.get("rifted_failed_m3", 0.0)),
            })
        for event in step.get("events", []):
            if event.get("kind") == "rift_craton":
                rift_path = event["path"]
                path_audit = rift_event_paths[rift_path]
                plate_id = str(int(event["plate_id"]))
                plate_audit = path_audit["by_plate"].setdefault(plate_id, {
                    "events": 0, "volume_m3": 0.0, "affected_cell_visits": 0,
                })
                volume = float(event["craton_loss_m3"])
                cells = int(event["affected_cells"])
                path_audit["events"] += 1
                path_audit["volume_m3"] += volume
                path_audit["cellwise_volume_m3"] += float(event["cellwise_loss_m3"])
                path_audit["affected_cell_visits"] += cells
                plate_audit["events"] += 1
                plate_audit["volume_m3"] += volume
                plate_audit["affected_cell_visits"] += cells
            if event.get("kind") != "plate_merge":
                continue
            cap = float(event.get("merge_capped_delta_m3", 0.0))
            live_delta = float(event.get("craton_live_delta_m3", 0.0))
            merge_event_cap_m3 += cap
            merge_event_cap_count += int(cap > 0.0)
            merge_live_cap_residual_m3 += live_delta + cap
            if cap > 0.0:
                merge_events.append({
                    "step": int(event["step"]),
                    "elapsed_myr": float(event["elapsed_myr"]),
                    "keep_plate_id": int(event["keep_plate_id"]),
                    "absorbed_plate_id": int(event["absorbed_plate_id"]),
                    "merge_capped_m3": cap,
                    "craton_live_before_m3": float(event["craton_live_before_m3"]),
                    "craton_live_after_m3": float(event["craton_live_after_m3"]),
                    "craton_live_delta_m3": live_delta,
                })
        if formation_onset is None and cumulative["formed_m3"] > 0.0:
            formation_onset = float(step["elapsed_myr"])
        live = initial + cumulative["formed_m3"] - sum(cumulative[name] for name in CRATON_SINKS)
        live += float(step["craton_balance_error_m3"])
        elapsed = float(step["elapsed_myr"])
        live_series.append((elapsed, live))
        step_number = int(step["step"])
        if step_number % CHECKPOINT_MYR == 0 or step_number == final["steps"]:
            checkpoints.append({
                "elapsed_myr": elapsed,
                "initial_mkm3": initial / MKM3,
                "formed_mkm3": cumulative["formed_m3"] / MKM3,
                "sinks_mkm3": {name: cumulative[name] / MKM3 for name in CRATON_SINKS},
                "rift_paths_mkm3": {name: cumulative[name] / MKM3 for name in CRATON_RIFT_PATHS},
                "live_mkm3": live / MKM3,
                "balance_error_m3": float(step["craton_balance_error_m3"]),
            })

    return {
        "model": path.name.split("-seed", 1)[0],
        "seed": int(final["seed"]),
        "commit": final["commit"],
        "steps": int(final["steps"]),
        "step_myr": float(final["step_myr"]),
        "initial_mkm3": initial / MKM3,
        "formed_mkm3": cumulative["formed_m3"] / MKM3,
        "sinks_mkm3": {name: cumulative[name] / MKM3 for name in CRATON_SINKS},
        "rift_paths_mkm3": {name: cumulative[name] / MKM3 for name in CRATON_RIFT_PATHS},
        "rift_path_reconciliation_mkm3": sum(cumulative[name] for name in CRATON_RIFT_PATHS) / MKM3 - cumulative["rifted_m3"] / MKM3,
        "live_mkm3": live_series[-1][1] / MKM3,
        "balance_error_m3": float(final["craton_balance_error_m3"]),
        "formation_onset_myr": formation_onset,
        "event_trace_audit": {
            "merge_event_count": sum(1 for step in steps for event in step.get("events", []) if event.get("kind") == "plate_merge"),
            "merge_capped_event_count": merge_event_cap_count,
            "merge_event_cap_mkm3": merge_event_cap_m3 / MKM3,
            "merge_step_ledger_cap_mkm3": cumulative["merge_capped_m3"] / MKM3,
            "merge_event_to_ledger_residual_m3": merge_event_cap_m3 - cumulative["merge_capped_m3"],
            "merge_live_delta_plus_cap_residual_m3": merge_live_cap_residual_m3,
            "merge_capped_events": merge_events,
            "rift_activity": {
                name: {"nonzero_steps": value["nonzero_steps"], "mkm3": value["events_m3"] / MKM3}
                for name, value in rift_activity.items()
            },
            "failed_rift": {
                **{key: value for key, value in failed_rift.items() if key != "cut_steps"},
                "band_craton_area_fraction": (
                    failed_rift["band_craton_area_m2"] / failed_rift["band_area_m2"]
                    if failed_rift["band_area_m2"] else 0.0
                ),
                "cut_steps": failed_rift["cut_steps"],
            },
            "rift_event_paths": {
                path: {
                    "events": audit["events"],
                    "event_volume_mkm3": audit["volume_m3"] / MKM3,
                    "step_ledger_volume_mkm3": cumulative[audit["ledger_account"]] / MKM3,
                    "event_to_ledger_residual_m3": audit["volume_m3"] - cumulative[audit["ledger_account"]],
                    "cellwise_volume_mkm3": audit["cellwise_volume_m3"] / MKM3,
                    "ledger_vs_cellwise_residual_m3": audit["volume_m3"] - audit["cellwise_volume_m3"],
                    "affected_cell_visits": audit["affected_cell_visits"],
                    "by_plate": [
                        {"plate_id": int(plate_id), **values,
                         "volume_mkm3": values["volume_m3"] / MKM3}
                        for plate_id, values in sorted(
                            audit["by_plate"].items(),
                            key=lambda item: item[1]["volume_m3"], reverse=True,
                        )
                    ],
                }
                for path, audit in rift_event_paths.items()
            },
        },
        "checkpoint_series": checkpoints,
        "live_series": live_series,
    }


def _first_crossing(series: list[tuple[float, float]], threshold: float) -> float | None:
    for elapsed, value in series:
        if value >= threshold:
            return elapsed
    return None


def _first_upward_zero_crossing(series: list[tuple[float, float]]) -> float | None:
    was_negative = False
    for elapsed, value in series:
        if value < 0.0:
            was_negative = True
        elif was_negative:
            return elapsed
    return None


def _matched(baseline: dict, pilot: dict) -> dict:
    if baseline["seed"] != pilot["seed"]:
        raise ValueError("baseline/pilot seeds do not match")
    base_checkpoints = baseline["checkpoint_series"]
    pilot_checkpoints = pilot["checkpoint_series"]
    if [row["elapsed_myr"] for row in base_checkpoints] != [row["elapsed_myr"] for row in pilot_checkpoints]:
        raise ValueError(f"checkpoint mismatch for seed {baseline['seed']}")

    initial_delta = pilot["initial_mkm3"] - baseline["initial_mkm3"]
    formation_delta = pilot["formed_mkm3"] - baseline["formed_mkm3"]
    sink_delta = {
        name: pilot["sinks_mkm3"][name] - baseline["sinks_mkm3"][name]
        for name in CRATON_SINKS
    }
    rift_path_delta = {
        name: pilot["rift_paths_mkm3"][name] - baseline["rift_paths_mkm3"][name]
        for name in CRATON_RIFT_PATHS
    }
    live_delta = pilot["live_mkm3"] - baseline["live_mkm3"]
    attribution = {"initial_inventory": initial_delta, "additional_formation": formation_delta}
    attribution.update({f"avoided_{name.removesuffix('_m3')}": -delta for name, delta in sink_delta.items()})
    residual = live_delta - sum(attribution.values())
    exact_gap_series = [
        (elapsed, pilot_value / MKM3 - base_value / MKM3)
        for (elapsed, base_value), (_, pilot_value) in zip(baseline["live_series"], pilot["live_series"], strict=True)
    ]
    final_gap = exact_gap_series[-1][1]
    crossings = {
        f"{fraction:g}_of_final_gap_myr": _first_crossing(exact_gap_series, fraction * final_gap)
        for fraction in (0.1, 0.5, 0.9)
    }
    crossings["first_upward_zero_crossing_myr"] = _first_upward_zero_crossing(exact_gap_series)

    return {
        "seed": baseline["seed"],
        "baseline_commit": baseline["commit"],
        "pilot_commit": pilot["commit"],
        "initial_inventory_delta_mkm3": initial_delta,
        "baseline_formation_onset_myr": baseline["formation_onset_myr"],
        "pilot_formation_onset_myr": pilot["formation_onset_myr"],
        "baseline_event_trace_audit": baseline["event_trace_audit"],
        "pilot_event_trace_audit": pilot["event_trace_audit"],
        "final_live_craton_delta_mkm3": live_delta,
        "final_attribution_mkm3": attribution,
        "pilot_minus_baseline_rift_path_delta_mkm3": rift_path_delta,
        "baseline_rift_path_reconciliation_mkm3": baseline["rift_path_reconciliation_mkm3"],
        "pilot_rift_path_reconciliation_mkm3": pilot["rift_path_reconciliation_mkm3"],
        "reconciliation_residual_mkm3": residual,
        "threshold_crossings": crossings,
        "checkpoint_deltas": [
            {
                "elapsed_myr": base_row["elapsed_myr"],
                "live_craton_delta_mkm3": round(pilot_row["live_mkm3"] - base_row["live_mkm3"], 6),
                "formation_delta_mkm3": round(pilot_row["formed_mkm3"] - base_row["formed_mkm3"], 6),
                "sink_deltas_mkm3": {
                    name: round(pilot_row["sinks_mkm3"][name] - base_row["sinks_mkm3"][name], 6)
                    for name in CRATON_SINKS
                },
                "rift_path_deltas_mkm3": {
                    name: round(pilot_row["rift_paths_mkm3"][name] - base_row["rift_paths_mkm3"][name], 6)
                    for name in CRATON_RIFT_PATHS
                },
            }
            for base_row, pilot_row in zip(base_checkpoints, pilot_checkpoints, strict=True)
        ],
        "baseline_final": {"initial_mkm3": baseline["initial_mkm3"], "formed_mkm3": baseline["formed_mkm3"], "sinks_mkm3": baseline["sinks_mkm3"], "live_mkm3": baseline["live_mkm3"]},
        "pilot_final": {"initial_mkm3": pilot["initial_mkm3"], "formed_mkm3": pilot["formed_mkm3"], "sinks_mkm3": pilot["sinks_mkm3"], "live_mkm3": pilot["live_mkm3"]},
    }


def analyze(trace_dir: Path) -> dict:
    paths = sorted(trace_dir.glob("baseline-seed*-1myr-500myr.jsonl")) + sorted(
        trace_dir.glob("pilot-seed*-1myr-500myr.jsonl")
    )
    runs = [_read_run(path) for path in paths]
    by_model_seed = {(run["model"], run["seed"]): run for run in runs}
    seeds = sorted({run["seed"] for run in runs})
    matched = [
        _matched(by_model_seed[("baseline", seed)], by_model_seed[("pilot", seed)])
        for seed in seeds
    ]
    aggregate = {
        "initial_inventory_delta_mkm3": sum(run["initial_inventory_delta_mkm3"] for run in matched),
        "final_live_craton_delta_mkm3": sum(run["final_live_craton_delta_mkm3"] for run in matched),
        "final_attribution_mkm3": {
            name: sum(run["final_attribution_mkm3"][name] for run in matched)
            for name in matched[0]["final_attribution_mkm3"]
        },
        "pilot_minus_baseline_rift_path_delta_mkm3": {
            name: sum(run["pilot_minus_baseline_rift_path_delta_mkm3"][name] for run in matched)
            for name in CRATON_RIFT_PATHS
        },
        "checkpoint_deltas": [],
    }
    for index, base_row in enumerate(matched[0]["checkpoint_deltas"]):
        aggregate["checkpoint_deltas"].append({
            "elapsed_myr": base_row["elapsed_myr"],
            "live_craton_delta_mkm3": sum(run["checkpoint_deltas"][index]["live_craton_delta_mkm3"] for run in matched),
            "formation_delta_mkm3": sum(run["checkpoint_deltas"][index]["formation_delta_mkm3"] for run in matched),
            "sink_deltas_mkm3": {
                name: sum(run["checkpoint_deltas"][index]["sink_deltas_mkm3"][name] for run in matched)
                for name in CRATON_SINKS
            },
        })
    return {
        "diagnostic": "RFC 0012 matched craton-ledger divergence and attribution",
        "trace_dir": str(trace_dir),
        "checkpoint_interval_myr": CHECKPOINT_MYR,
        "accounting": "live = initial + formed - all declared craton sinks + recorded balance error",
        "runs": len(runs),
        "matched_seed_count": len(matched),
        "matched": matched,
        "aggregate": aggregate,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace_dir", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.trace_dir)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"runs": result["runs"], "matched_seeds": result["matched_seed_count"], "out": str(args.out)}))


if __name__ == "__main__":
    main()

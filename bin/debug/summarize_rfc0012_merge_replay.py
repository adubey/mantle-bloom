#!/usr/bin/env python3
"""Compare the preregistered RFC 0012 #547 merge-only replay matrix.

The input directory contains one JSON output per baseline/candidate arm, seed, and step
size, written by `replay_rfc0012_merge_rate.py`. The script verifies the full 4-seed ×
3-step × 2-arm matrix and writes a compact comparison JSON; raw replay outputs stay in
scratch storage.

    backend/.venv/bin/python bin/debug/summarize_rfc0012_merge_replay.py \
        --input-dir <scratch-dir> --out analysis/issue547/merge_only_replay_summary.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

SEEDS = (326494039, 748935664, 580014012, 111065103)
STEPS_MYR = (0.1, 1.0, 2.0)
ARMS = ("baseline", "candidate")
MOTION_METRICS = (
    "accumulated_driving_torque_nm_year",
    "accumulated_residual_torque_nm_year",
    "mean_solved_rate_cm_per_yr",
    "peak_solved_rate_cm_per_yr",
    "mean_absolute_closing_m_per_s",
    "peak_absolute_closing_m_per_s",
    "mean_contact_force_n_per_m",
)
QUANTILES = ("p50", "p90", "p99")
QUANTILE_FIELDS = ("hc_area_weighted_km", "hm_area_weighted_km")
M_PER_S_TO_CM_PER_YR = 3_155_760_000.0


def _percent_delta(baseline: float, candidate: float) -> float | None:
    if abs(baseline) < 1e-30:
        return 0.0 if abs(candidate) < 1e-30 else None
    return (candidate - baseline) / abs(baseline) * 100.0


def _load_run(input_dir: Path, seed: int, step_myr: float, arm: str) -> dict[str, Any]:
    path = input_dir / f"seed{seed}-step{step_myr:g}-{arm}.json"
    with path.open(encoding="utf-8") as handle:
        run = json.load(handle)
    expected = (seed, step_myr, 50.0, 4.0)
    observed = (run["seed"], run["step_myr"], run["horizon_myr"], run["node_density"])
    if observed != expected:
        raise ValueError(f"{path.name}: expected seed/step/horizon/density {expected}, got {observed}")
    if len(run["checkpoints"]) != 6 or run["checkpoints"][-1]["myr"] != 50.0:
        raise ValueError(f"{path.name}: expected checkpoints at 0, 10, 20, 30, 40, and 50 Myr")
    return run


def _topology_summary(baseline: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    base_events = baseline["topology_events"]
    candidate_events = candidate["topology_events"]
    common_prefix = 0
    while (
        common_prefix < min(len(base_events), len(candidate_events))
        and base_events[common_prefix] == candidate_events[common_prefix]
    ):
        common_prefix += 1

    def fusion_count(events: list[list[Any]]) -> int:
        return sum(
            any(term in message.lower() for term in ("merged into", "fused into", "accreted onto"))
            for _, message in events
        )

    first_difference = None
    if base_events != candidate_events:
        first_difference = {
            "baseline": base_events[common_prefix] if common_prefix < len(base_events) else None,
            "candidate": candidate_events[common_prefix] if common_prefix < len(candidate_events) else None,
        }
        for arm in ("baseline", "candidate"):
            event = first_difference[arm]
            if event is not None:
                first_difference[arm] = {"time_myr": event[0] / 1e6, "event": event[1]}

    return {
        "sequences_identical": base_events == candidate_events,
        "baseline_event_count": len(base_events),
        "candidate_event_count": len(candidate_events),
        "baseline_fusion_or_terrane_accretion_count": fusion_count(base_events),
        "candidate_fusion_or_terrane_accretion_count": fusion_count(candidate_events),
        "common_prefix_event_count": common_prefix,
        "first_difference": first_difference,
    }


def _checkpoint_trajectory(baseline: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    times = []
    rows = {
        "plate_count_baseline": [],
        "plate_count_candidate": [],
        "land_fraction_baseline": [],
        "land_fraction_candidate": [],
        "max_abs_hc_hm_quantile_delta_percent": [],
    }
    for base, changed in zip(baseline["checkpoints"], candidate["checkpoints"], strict=True):
        quantile_delta = {
            field: {
                quantile: _percent_delta(base[field][quantile], changed[field][quantile])
                for quantile in QUANTILES
            }
            for field in QUANTILE_FIELDS
        }
        times.append(base["myr"])
        rows["plate_count_baseline"].append(base["plate_count"])
        rows["plate_count_candidate"].append(changed["plate_count"])
        rows["land_fraction_baseline"].append(base["land_fraction"])
        rows["land_fraction_candidate"].append(changed["land_fraction"])
        rows["max_abs_hc_hm_quantile_delta_percent"].append(
            max(
                abs(value)
                for family in quantile_delta.values()
                for value in family.values()
                if value is not None
            )
        )
    return {"myr": times, **rows}


def _endpoint_summary(checkpoint: dict[str, Any]) -> dict[str, Any]:
    return {
        "myr": checkpoint["myr"],
        "plate_count": checkpoint["plate_count"],
        "land_fraction": checkpoint["land_fraction"],
        "hc_area_weighted_km": checkpoint["hc_area_weighted_km"],
        "hm_area_weighted_km": checkpoint["hm_area_weighted_km"],
        "ledger_balance_error_m3": checkpoint["ledger_balance_error_m3"],
        "hm_ledger_closure_error_m3": checkpoint["hm_ledger_closure_error_m3"],
        "hm_inventory_m3": checkpoint["hm_inventory_m3"],
    }


def _ledger_maxima(runs: list[dict[str, Any]]) -> dict[str, Any]:
    max_abs_balance = {name: 0.0 for name in ("continental", "craton", "mobile_cover")}
    max_abs_hm_closure = {name: 0.0 for name in ("all", "continental_node", "oceanic_node")}
    max_relative_hm_closure = dict.fromkeys(max_abs_hm_closure, 0.0)
    for run in runs:
        for checkpoint in run["checkpoints"]:
            for name, value in checkpoint["ledger_balance_error_m3"].items():
                if name in max_abs_balance:
                    max_abs_balance[name] = max(max_abs_balance[name], abs(value))
            for name, value in checkpoint.get("hm_ledger_closure_error_m3", {}).items():
                if name in max_abs_hm_closure:
                    max_abs_hm_closure[name] = max(max_abs_hm_closure[name], abs(value))
                    inventory = checkpoint["hm_inventory_m3"][name]
                    relative = abs(value) / inventory if inventory else 0.0
                    max_relative_hm_closure[name] = max(max_relative_hm_closure[name], relative)
    return {
        "max_absolute_balance_error_m3": max_abs_balance,
        "max_absolute_hm_closure_error_m3": max_abs_hm_closure,
        "max_relative_hm_closure_fraction": max_relative_hm_closure,
    }


def _make_summary(input_dir: Path) -> dict[str, Any]:
    runs = {
        (seed, step, arm): _load_run(input_dir, seed, step, arm)
        for seed in SEEDS
        for step in STEPS_MYR
        for arm in ARMS
    }
    baseline_commits = {run["commit"] for (seed, step, arm), run in runs.items() if arm == "baseline"}
    candidate_commits = {run["commit"] for (seed, step, arm), run in runs.items() if arm == "candidate"}
    if len(baseline_commits) != 1 or len(candidate_commits) != 1:
        raise ValueError(f"commit mismatch in matrix: baseline={baseline_commits}, candidate={candidate_commits}")

    pairs = []
    motion_maxima: dict[str, dict[str, Any]] = {}
    topology_differences = 0
    plate_count_differences = 0
    max_land_delta_pp = 0.0
    max_quantile_delta = {field: {quantile: 0.0 for quantile in QUANTILES} for field in QUANTILE_FIELDS}
    for seed in SEEDS:
        for step in STEPS_MYR:
            baseline = runs[(seed, step, "baseline")]
            candidate = runs[(seed, step, "candidate")]
            base_end = baseline["checkpoints"][-1]
            candidate_end = candidate["checkpoints"][-1]
            motion_delta = {
                metric: _percent_delta(baseline["motion"][metric], candidate["motion"][metric])
                for metric in MOTION_METRICS
            }
            for metric, delta in motion_delta.items():
                if delta is not None and (
                    metric not in motion_maxima or abs(delta) > abs(motion_maxima[metric]["delta_percent"])
                ):
                    motion_maxima[metric] = {
                        "seed": seed,
                        "step_myr": step,
                        "delta_percent": delta,
                    }
            topology = _topology_summary(baseline, candidate)
            topology_differences += not topology["sequences_identical"]
            plate_count_differences += base_end["plate_count"] != candidate_end["plate_count"]
            land_delta_pp = (candidate_end["land_fraction"] - base_end["land_fraction"]) * 100.0
            max_land_delta_pp = max(max_land_delta_pp, abs(land_delta_pp))
            quantile_delta = {
                field: {
                    quantile: _percent_delta(base_end[field][quantile], candidate_end[field][quantile])
                    for quantile in QUANTILES
                }
                for field in QUANTILE_FIELDS
            }
            for field, family in quantile_delta.items():
                for quantile, delta in family.items():
                    if delta is not None:
                        max_quantile_delta[field][quantile] = max(
                            max_quantile_delta[field][quantile], abs(delta)
                        )

            deformation_delta = []
            for phase in sorted(set(baseline["deformation"]) | set(candidate["deformation"])):
                for measure in ("hc_m3", "hm_m3"):
                    base_value = baseline["deformation"].get(phase, {}).get(measure, 0.0)
                    candidate_value = candidate["deformation"].get(phase, {}).get(measure, 0.0)
                    deformation_delta.append(
                        {
                            "phase": phase,
                            "measure": measure,
                            "baseline_m3": base_value,
                            "candidate_m3": candidate_value,
                            "delta_m3": candidate_value - base_value,
                            "delta_percent": _percent_delta(base_value, candidate_value),
                        }
                    )

            largest_deformation_deltas = sorted(
                deformation_delta, key=lambda row: abs(row["delta_m3"]), reverse=True
            )[:3]
            pairs.append(
                {
                    "seed": seed,
                    "step_myr": step,
                    "baseline_commit": baseline["commit"],
                    "candidate_commit": candidate["commit"],
                    "motion": {
                        metric: {
                            "baseline": baseline["motion"][metric],
                            "candidate": candidate["motion"][metric],
                            "delta_percent": motion_delta[metric],
                        }
                        for metric in MOTION_METRICS
                    },
                    "endpoint": {
                        "baseline": _endpoint_summary(base_end),
                        "candidate": _endpoint_summary(candidate_end),
                        "plate_count_delta": candidate_end["plate_count"] - base_end["plate_count"],
                        "land_fraction_delta_percentage_points": land_delta_pp,
                        "hc_hm_quantile_delta_percent": quantile_delta,
                    },
                    "topology": topology,
                    "checkpoint_trajectory": _checkpoint_trajectory(baseline, candidate),
                    "runtime": {
                        "baseline_seconds": baseline["runtime_seconds"],
                        "candidate_seconds": candidate["runtime_seconds"],
                        "delta_percent": _percent_delta(
                            baseline["runtime_seconds"], candidate["runtime_seconds"]
                        ),
                    },
                    "peak_rss_mib": {
                        "baseline": baseline["peak_rss_mib"],
                        "candidate": candidate["peak_rss_mib"],
                        "delta": candidate["peak_rss_mib"] - baseline["peak_rss_mib"],
                    },
                    "largest_absolute_deformation_deltas": largest_deformation_deltas,
                }
            )

    timestep_spread = []
    for seed in SEEDS:
        for arm in ARMS:
            reference = runs[(seed, 0.1, arm)]
            for step in (1.0, 2.0):
                comparison = runs[(seed, step, arm)]
                deltas = {
                    metric: _percent_delta(reference["motion"][metric], comparison["motion"][metric])
                    for metric in MOTION_METRICS
                }
                ref_end = reference["checkpoints"][-1]
                step_end = comparison["checkpoints"][-1]
                deltas.update(
                    {
                        "land_fraction_percent": _percent_delta(
                            ref_end["land_fraction"], step_end["land_fraction"]
                        ),
                        "plate_count_difference": step_end["plate_count"] - ref_end["plate_count"],
                    }
                )
                timestep_spread.append({"seed": seed, "arm": arm, "step_myr": step, "delta_vs_0_1": deltas})

    all_runs = list(runs.values())
    return {
        "schema_version": 1,
        "scope": "merge-only RFC 0012 #547 matched metric replay",
        "controls": {
            "seeds": list(SEEDS),
            "step_myr": list(STEPS_MYR),
            "horizon_myr": 50.0,
            "checkpoint_myr": 10.0,
            "node_density": 4.0,
            "pythonhashseed": 0,
            "baseline_commit": next(iter(baseline_commits)),
            "candidate_commit": next(iter(candidate_commits)),
            "runs": len(all_runs),
            "matched_pairs": len(pairs),
        },
        "summary": {
            "topology_sequence_differences": topology_differences,
            "endpoint_plate_count_differences": plate_count_differences,
            "max_absolute_land_fraction_delta_percentage_points": max_land_delta_pp,
            "max_absolute_endpoint_hc_hm_quantile_delta_percent": max(
                value for field in max_quantile_delta.values() for value in field.values()
            ),
            "max_absolute_motion_metric_delta_percent_by_metric": motion_maxima,
            "ledger_maxima": _ledger_maxima(all_runs),
        },
        "pairs": pairs,
        "timestep_spread_vs_0_1_myr": timestep_spread,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, help="Output path; defaults to stdout.")
    args = parser.parse_args()
    summary = _make_summary(args.input_dir)
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.out is None:
        print(rendered, end="")
        return
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()

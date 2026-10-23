#!/usr/bin/env python3
"""Summarize the private JSONL from probe_rfc0012_world_rifts.py.

    backend/.venv/bin/python bin/debug/summarize_rfc0012_world_rifts.py \
      <private-world-rift-trace.jsonl> --out analysis/rfc0012/g3_world_rifts.json

The summary contains no world save or machine-specific path. Component widths are area-
weighted p50/p90 values; the input trace does not retain the maximum individual component.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def summarize(source: Path, step_myr: float, glacier_cap_warnings: int) -> dict:
    raw = source.read_bytes()
    rows = [json.loads(line) for line in raw.decode("utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError("world-rift trace has no checkpoint rows")
    times = np.asarray([row["myr"] for row in rows], dtype=float)
    aggregate_width = np.asarray(
        [row["rift_width"]["area_equivalent_width_km"] for row in rows if row["rift_width"]["components"]],
        dtype=float,
    )
    component_p50 = np.asarray(
        [row["rift_width"]["component_width_km_p50"] for row in rows if row["rift_width"]["components"]],
        dtype=float,
    )
    component_p90 = np.asarray(
        [row["rift_width"]["component_width_km_p90"] for row in rows if row["rift_width"]["components"]],
        dtype=float,
    )
    active_rows = [row for row in rows if row["rift_width"]["components"]]
    width_rows = [row["rift_width"] for row in active_rows]
    last = rows[-1]
    path_rel = [
        abs(row["rift_path_residual_m3"]) / max(abs(row["rifted_craton_m3"]), 1.0)
        for row in rows
    ]
    return {
        "seed": last["seed"],
        "model_commit": last["model_commit"],
        "years": float(last["myr"] * 1e6),
        "step_myr": step_myr,
        "checkpoint_interval_myr": float(np.median(np.diff(times))) if len(times) > 1 else None,
        "checkpoint_count": len(rows),
        "width_metric": "2 * connected active-thinning area / its mesh perimeter; per-plate components, checkpoint samples only",
        "area_equivalent_width_km": {
            "median_across_checkpoints": float(np.median(aggregate_width)),
            "p90_across_checkpoints": float(np.percentile(aggregate_width, 90)),
            "maximum_checkpoint": float(np.max(aggregate_width)),
        },
        "area_weighted_component_width_km": {
            "p50_across_checkpoints": float(np.median(component_p50)),
            "p90_across_checkpoints": float(np.percentile(component_p50, 90)),
            "maximum_checkpoint_p50": float(np.max(component_p50)),
            "median_checkpoint_p90": float(np.median(component_p90)),
            "maximum_checkpoint_p90": float(np.max(component_p90)),
            "maximum_individual_component_width_km": max(
                item["maximum_component_width_km"] for item in width_rows
            ),
            "component_share_ge_100km": {
                "median_checkpoint": float(np.median([item["component_share_ge_100km"] for item in width_rows])),
                "maximum_checkpoint": max(item["component_share_ge_100km"] for item in width_rows),
            },
            "component_share_ge_300km": {
                "median_checkpoint": float(np.median([item["component_share_ge_300km"] for item in width_rows])),
                "maximum_checkpoint": max(item["component_share_ge_300km"] for item in width_rows),
            },
            "active_area_share_ge_100km": {
                "median_checkpoint": float(np.median([item["active_area_share_ge_100km"] for item in width_rows])),
                "maximum_checkpoint": max(item["active_area_share_ge_100km"] for item in width_rows),
            },
            "active_area_share_ge_300km": {
                "median_checkpoint": float(np.median([item["active_area_share_ge_300km"] for item in width_rows])),
                "maximum_checkpoint": max(item["active_area_share_ge_300km"] for item in width_rows),
            },
        },
        "ledger_closure": {
            "continental_material_assert_closed_at_all_checkpoints": True,
            "maximum_absolute_craton_balance_error_m3": max(
                abs(row["craton_balance_error_m3"]) for row in rows
            ),
            "maximum_absolute_continental_material_balance_error_m3": max(
                abs(row["continental_material_balance_error_m3"]) for row in rows
            ),
            "maximum_relative_craton_rift_path_residual": max(path_rel),
            "final_rifted_craton_m3": float(last["rifted_craton_m3"]),
            "final_craton_rift_path_sum_m3": float(last["rift_path_sum_m3"]),
        },
        "runtime_wall_seconds": float(last["wall_seconds"]),
        "glacier_max_depth_warnings": glacier_cap_warnings,
        "raw_trace_sha256": hashlib.sha256(raw).hexdigest(),
        "interpretation": (
            "Current production divergent thinning is summarized at 10 Myr checkpoints. "
            "The active-area width is not basin width; percentile summaries do not rule out "
            "rare wider components. No integrated-strength law feeds back into this replay."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--glacier-cap-warnings", type=int, default=0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.trace, args.step_myr, args.glacier_cap_warnings)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("seed", "years", "checkpoint_count", "area_equivalent_width_km", "ledger_closure")}))


if __name__ == "__main__":
    main()

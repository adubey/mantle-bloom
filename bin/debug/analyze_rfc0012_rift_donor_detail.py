#!/usr/bin/env python3
"""Summarize the targeted seed-580 rift donor-cell traces.

Reads the detailed baseline and pilot JSONL files produced by
``trace_rfc0012_ledger_gates.py`` and writes event-level cut/donor context plus compact
per-path and per-cell summaries. Full step traces stay outside the repository.

    backend/.venv/bin/python bin/debug/analyze_rfc0012_rift_donor_detail.py \
      /private/tmp/rfc0012-rift-donor-detail \
      --out analysis/rfc0012/g2_rift_donor_detail.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

MKM3 = 1.0e15
RIFT_PATHS = {"stretch_donation", "failed_rift_thinning"}


def _read_trace(path: Path) -> tuple[dict, list[dict], list[dict]]:
    summary = None
    events = []
    cut_steps = []
    for line in path.open():
        row = json.loads(line)
        if row["kind"] == "summary":
            summary = row
        elif row["kind"] == "step":
            ledger = row.get("ledger_delta", {}).get("craton_ledger", {})
            cuts = int(ledger.get("failed_rift_cuts", 0.0))
            if cuts:
                cut_steps.append({
                    "step": int(row["step"]),
                    "elapsed_myr": float(row["elapsed_myr"]),
                    "cuts": cuts,
                    "cuts_through_craton": int(ledger.get("failed_rift_cuts_through_craton", 0.0)),
                    "craton_loss_mkm3": float(ledger.get("rifted_failed_m3", 0.0)) / MKM3,
                    "band_area_m2": float(ledger.get("failed_rift_band_area_m2", 0.0)),
                    "band_craton_area_m2": float(ledger.get("failed_rift_band_craton_area_m2", 0.0)),
                })
            events.extend(
                event
                for event in row.get("events", [])
                if event.get("kind") == "rift_craton"
                and event.get("path") in RIFT_PATHS
                and event.get("rift_context")
            )
    if summary is None:
        raise ValueError(f"missing summary row: {path}")
    return summary, events, cut_steps


def _summarize_model(path: Path) -> dict:
    summary, events, cut_steps = _read_trace(path)
    by_path = {}
    donor_cells = {}
    failed_rift_events = []
    stretch_events = []
    for event in events:
        rift_path = event["path"]
        aggregate = by_path.setdefault(rift_path, {
            "operations": 0, "craton_loss_mkm3": 0.0, "donor_cell_visits": 0,
        })
        donors = event["rift_context"].get("donor_cells", [])
        volume = float(event["craton_loss_m3"])
        aggregate["operations"] += 1
        aggregate["craton_loss_mkm3"] += volume / MKM3
        aggregate["donor_cell_visits"] += len(donors)
        for cell in donors:
            key = (
                rift_path,
                int(event["plate_id"]),
                int(event["rift_context"].get("plate_state_before", {}).get("topology_revision", -1)),
                tuple(cell["node_id"]),
            )
            row = donor_cells.setdefault(key, {
                "path": rift_path,
                "plate_id": int(event["plate_id"]),
                "topology_revision": key[2],
                "node_id": list(key[3]),
                "local_xyz": cell["local_xyz"],
                "operations": 0,
                "craton_loss_mkm3": 0.0,
                "max_single_operation_loss_mkm3": 0.0,
            })
            cell_loss = float(cell["loss_m3"])
            row["operations"] += 1
            row["craton_loss_mkm3"] += cell_loss / MKM3
            row["max_single_operation_loss_mkm3"] = max(
                row["max_single_operation_loss_mkm3"], cell_loss / MKM3
            )
        context = event["rift_context"]
        pre_state = context.get("plate_state_before", {})
        cut = context.get("cut", {})
        donor_inventory = sum(float(cell["craton_before_m"]) * float(cell["area_m2"]) for cell in donors)
        row = {
            "step": int(event["step"]),
            "elapsed_myr": float(event["elapsed_myr"]),
            "path": rift_path,
            "plate_id": int(event["plate_id"]),
            "craton_loss_mkm3": volume / MKM3,
            "affected_cells": int(event["affected_cells"]),
            "ledger_vs_cellwise_residual_m3": float(event["ledger_vs_cellwise_residual_m3"]),
            "donor_cell_visits": len(donors),
            "donor_inventory_mkm3": donor_inventory / MKM3,
            "plate_state": {
                key: pre_state.get(key)
                for key in (
                    "plate_id", "node_count", "age_myr", "internal_stress", "craton_inventory_m3",
                    "omega_rad_per_year", "frame", "topology_revision", "geometry_revision",
                    "split_gates", "neighbour_plates",
                )
                if key in pre_state
            },
        }
        if rift_path == "failed_rift_thinning":
            gate = pre_state.get("split_gates", {})
            residual = context.get("flow_fit_rms_rad_per_year")
            pole_gap = context.get("fitted_pole_separation_rad_per_year")
            band_inventory = float(cut.get("band_craton_inventory_m3", 0.0))
            row["failed_cut"] = {
                "cut_normal_world": cut.get("cut_normal_world"),
                "band_cell_count": cut.get("band_cell_count", 0),
                "band_craton_inventory_mkm3": band_inventory / MKM3,
                "loss_fraction_of_band_inventory": volume / band_inventory if band_inventory else 0.0,
                "flow_fit_rms_rad_per_year": residual,
                "flow_residual_threshold_rad_per_year": gate.get("flow_residual_threshold_rad_per_year"),
                "flow_residual_over_threshold": (
                    residual / gate["flow_residual_threshold_rad_per_year"]
                    if residual is not None and gate.get("flow_residual_threshold_rad_per_year") else None
                ),
                "pole_separation_rad_per_year": pole_gap,
                "pole_separation_threshold_rad_per_year": gate.get("pole_separation_threshold_rad_per_year"),
                "pole_separation_over_threshold": (
                    pole_gap / gate["pole_separation_threshold_rad_per_year"]
                    if pole_gap is not None and gate.get("pole_separation_threshold_rad_per_year") else None
                ),
            }
        else:
            row["stretch"] = {
                "rifted_new_cell_count": len(context.get("rifted_new_cells", [])),
                "stretch_share_area_weighted_mean": float(event.get("stretch_share_area_weighted_mean", 0.0)),
            }
        if rift_path == "failed_rift_thinning":
            failed_rift_events.append(row)
        else:
            stretch_events.append({
                "step": row["step"],
                "elapsed_myr": row["elapsed_myr"],
                "plate_id": row["plate_id"],
                "craton_loss_mkm3": row["craton_loss_mkm3"],
                "affected_cells": row["affected_cells"],
                "donor_cell_visits": row["donor_cell_visits"],
                "donor_inventory_mkm3": row["donor_inventory_mkm3"],
                "stretch": row["stretch"],
            })
    for aggregate in by_path.values():
        aggregate["craton_loss_mkm3"] = round(aggregate["craton_loss_mkm3"], 9)
    return {
        "model": path.name.split("-seed", 1)[0],
        "commit": summary["commit"],
        "failed_rift_ledger": {
            key: float(summary["ledgers"]["craton_ledger"].get(key, 0.0))
            for key in (
                "failed_rift_cuts", "failed_rift_cuts_through_craton", "rifted_failed_m3",
                "failed_rift_band_area_m2", "failed_rift_band_craton_area_m2",
            )
        },
        "failed_rift_cut_steps": cut_steps,
        "by_path": by_path,
        "top_donor_cells": sorted(
            donor_cells.values(), key=lambda row: row["craton_loss_mkm3"], reverse=True
        )[:30],
        "failed_rift_events": failed_rift_events,
        "largest_stretch_events": sorted(
            stretch_events, key=lambda row: row["craton_loss_mkm3"], reverse=True
        )[:10],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace_dir", type=Path)
    parser.add_argument("--seed", type=int, default=580014012)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    models = [
        _summarize_model(args.trace_dir / f"{name}-seed{args.seed}-1myr-500myr.jsonl")
        for name in ("baseline", "pilot")
    ]
    output = {
        "diagnostic": "RFC 0012 targeted rift donor cell and split-state trace",
        "seed": args.seed,
        "step_myr": 1.0,
        "horizon_myr": 500.0,
        "models": models,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    cut_count = sum(len(model["failed_rift_events"]) for model in models)
    print(f"wrote {args.out} ({cut_count} craton-bearing failed-rift events)")


if __name__ == "__main__":
    main()

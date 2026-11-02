#!/usr/bin/env python3
"""Issue #757: replay a checkpoint and record whole-world lake water, sea level and the largest lakes.

Steps a saved world on the selected backend with the ice-age cycle off, as #754's windows do, and
records per step:

- whole-world liquid lake volume, the same reading as #754's `world_lake_km3`, and lake ice;
- the wet component under each registered anchor, read as #754's `target_km3` and
  `target_area_km2`;
- sea level, and the lake ledger's evaporation, sublimation, phase exchange and overflow;
- the largest lakes by area: liquid and ice volume, water level and floor, and whether the
  sea-tier depth guard (`lakes.SEA_MAX_DEPTH_M`) holds the level;
- normalized lake and land residuals, and step runtime.

With `--audit-losses`, transient solve traces compare each wet final body's booked losses
against evaporation and sublimation evaluated at its final level over the full step (#768).

The summary gives the lag-1 autocorrelation of the step change in whole-world lake volume (about
-1 is period-2 alternation), the mean absolute step change from step 3, and the largest target
ratio between consecutive steps, all computed as #754 did. Nothing here changes the model, so a run
on #754's commit must match its control.

Example:

    backend/.venv/bin/python bin/debug/replay_lake_water_swing.py <checkpoint.mbworld> \\
        --anchors analysis/rfc0016/target-anchor-window-e.json --myr 1.2 \\
        --backend <baseline worktree>/backend --out <scratch>/issue757-e-main.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from replay_lib import load_replay, parse_replay_args, run_replay

M3_PER_KM3 = 1e9
M2_PER_KM2 = 1e6
LARGEST_LAKES = 5


def add_arguments(ap) -> None:
    ap.add_argument("--anchors", type=Path, help="JSON file containing anchor_points")
    ap.add_argument("--audit-losses", action="store_true", help="Audit wet final bodies against full-step loss rates (#768)")


def audit_lake_losses(forest, elevation, years, climate, hyps, solved):
    """Compare booked losses with independent final-level demands, excluding dry bodies.

    A dry body may consume its entire supply while ending with zero wet area, so the
    final-area equation is not an equality there. The audit is for the default single slice.
    """
    from app.hydroclimate import lakes

    records = []
    for body in lakes._active_bodies(
        forest, elevation, solved.lake_depth, hyps,
    ):
        members = body.members
        wet_members = members[solved.lake_depth[members] > 0.0]
        if not len(wet_members):
            continue
        level = float(hyps.level_from_depth(solved.lake_depth[wet_members], wet_members).max())
        wet = hyps.wet_fraction(level, members)
        expected_evap = float(np.dot(
            hyps.area_m2[members] * wet,
            climate.evaporation_m_per_yr[members] * (1.0 + lakes.SHALLOW_EVAPORATION_BOOST * (1.0 - wet)),
        )) * years
        expected_sub = (
            float(np.dot(hyps.area_m2[members] * wet, climate.sublimation_m_per_yr[members])) * years
            if climate.sublimation_m_per_yr is not None else 0.0
        )
        booked_evap = float(solved.trace.evaporation_m3[members].sum())
        booked_sub = float(solved.trace.liquid_sublimation_m3[members].sum())
        scale = max(expected_evap + expected_sub, 1.0)
        records.append({
            "lake_id": body.lake_id, "nodes": len(members), "level_m": level,
            "storage_km3": float(np.dot(solved.lake_depth[members], hyps.area_m2[members])) / M3_PER_KM3,
            "member_nodes": members,
            "expected_evaporation_km3": expected_evap / M3_PER_KM3,
            "expected_sublimation_km3": expected_sub / M3_PER_KM3,
            "booked_evaporation_km3": booked_evap / M3_PER_KM3,
            "booked_sublimation_km3": booked_sub / M3_PER_KM3,
            "normalized_loss_error": (abs(booked_evap - expected_evap) + abs(booked_sub - expected_sub)) / scale,
        })
    return records


def summarize_loss_audit(records, anchor_nodes):
    """Attach end-of-step anchor indices to loss rates captured before terrain changes."""
    for record in records:
        members = record.pop("member_nodes")
        record["anchor_indices"] = [i for i, node in enumerate(anchor_nodes) if node in members]
    return {
        "worst_normalized_loss_error": max((r["normalized_loss_error"] for r in records), default=0.0),
        "largest_bodies": sorted(records, key=lambda r: -r["storage_km3"])[:LARGEST_LAKES],
        "anchor_bodies": [r for r in records if r["anchor_indices"]],
    }


def lake_record(members, hydro, areas, sea_guard_m) -> dict:
    """One wet component: its size, water and level, with the level read through the same
    hypsometry the lake solve uses."""
    depth = np.asarray(hydro.lake_depth)
    levels = hydro.hypsometry.level_from_depth(depth[members], members)
    floor = float(np.min(hydro.elevation[members]))
    level = float(np.max(levels))
    centroid = hydro.points[members].mean(axis=0)
    return {
        "nodes": len(members),
        "area_km2": float(areas[members].sum()) / M2_PER_KM2,
        "liquid_km3": float(np.dot(depth[members], areas[members])) / M3_PER_KM3,
        "ice_km3": float(np.dot(hydro.lake_ice_m[members], areas[members])) / M3_PER_KM3,
        "floor_m": floor,
        "level_m": level,
        "is_sea_tier": bool(np.any(hydro.is_sea[members])),
        "at_sea_guard": level - floor >= sea_guard_m - 1.0,
        "centroid": (centroid / np.linalg.norm(centroid)).tolist(),
    }


def summarize(steps: list[dict]) -> dict:
    """#754's swing metrics over the replayed steps."""
    out = {}
    for key in ("world_lake_km3", "world_lake_water_km3"):
        volume = np.array([step[key] for step in steps])
        change = np.diff(volume)
        out[key] = {
            "lag1_autocorrelation": (
                float(np.corrcoef(change[:-1], change[1:])[0, 1]) if len(change) > 2 else None
            ),
            "mean_abs_change_from_step3_km3": float(np.abs(change[2:]).mean()) if len(change) > 2 else None,
        }
    target = np.array([step["target_km3"] for step in steps])
    ratios = [
        max(a, b) / max(min(a, b), 1e-9) for a, b in zip(target[:-1], target[1:]) if max(a, b) > 0.0
    ]
    out["largest_consecutive_target_ratio"] = float(max(ratios)) if ratios else None
    out["worst_normalized_lake_residual"] = max(step["normalized_lake_residual"] for step in steps)
    out["worst_normalized_land_residual"] = max(step["normalized_land_residual"] for step in steps)
    return out


def main() -> None:
    args = parse_replay_args(__doc__, add_arguments)
    world, result = load_replay(args, issue=757)
    world.ice_age_period_years = 0.0
    from app import plates as plates_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.hydroclimate import hydrology, lakes  # noqa: E402

    anchors = (
        np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)
        if args.anchors is not None else np.zeros((0, 3))
    )
    result.update(steps=[], anchors=anchors.tolist())
    captured = None
    if args.audit_losses:
        original_solve = lakes.solve_lakes_prognostic

        def capture_solve(forest, elevation, previous, years, climate, hypsometry, **kwargs):
            nonlocal captured
            kwargs["diagnostic_trace"] = True
            solved = original_solve(forest, elevation, previous, years, climate, hypsometry, **kwargs)
            # Erosion updates terrain after hydrology. Evaluate now, before its shared
            # elevation arrays can change, then attach geographic anchors after the step.
            captured = audit_lake_losses(forest, elevation, years, climate, hypsometry, solved)
            return solved

        lakes.solve_lakes_prognostic = capture_solve

    def snapshot() -> dict:
        return {"elapsed_myr": world.elapsed_years / 1e6, "sea_level_m": world.sea_level_m}

    def after_step(seconds: float) -> None:
        hydro = world.hydrology_cache
        areas = plates_mod.collect_all_accounting_areas_m2(
            hydro.plates_in_order, line_spacing_rad(world.node_density)
        )
        depth = np.asarray(hydro.lake_depth)
        ice = np.asarray(hydro.lake_ice_m)
        wet = depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
        components = sorted(
            (np.asarray(nodes) for nodes in hydrology.lake_components(wet, hydro.neighbor_idx)),
            key=lambda nodes: -float(areas[nodes].sum()),
        )
        component_of = np.full(len(depth), -1, dtype=np.int64)
        for i, nodes in enumerate(components):
            component_of[nodes] = i
        targets = []
        for anchor in anchors:
            node = int(np.argmin(np.linalg.norm(hydro.points - anchor, axis=1)))
            index = int(component_of[node])
            targets.append(
                lake_record(components[index], hydro, areas, lakes.SEA_MAX_DEPTH_M) | {"rank": index}
                if index >= 0 else None
            )
        first = targets[0] if targets else None
        lake, land = hydro.lake_water_ledger, hydro.land_water_ledger
        result["steps"].append({
            "step": len(result["steps"]) + 1,
            "elapsed_myr": world.elapsed_years / 1e6,
            "runtime_s": seconds,
            "sea_level_m": world.sea_level_m,
            "world_lake_km3": float(np.dot(depth, areas)) / M3_PER_KM3,
            "world_lake_ice_km3": float(np.dot(ice, areas)) / M3_PER_KM3,
            "world_lake_water_km3": float(np.dot(depth + ice, areas)) / M3_PER_KM3,
            "lake_area_km2": float(areas[wet].sum()) / M2_PER_KM2,
            "target_km3": first["liquid_km3"] if first else 0.0,
            "target_area_km2": first["area_km2"] if first else 0.0,
            "targets": targets,
            "largest_lakes": [
                lake_record(nodes, hydro, areas, lakes.SEA_MAX_DEPTH_M) for nodes in components[:LARGEST_LAKES]
            ],
            "lake_ledger_km3": {
                key: getattr(lake, key) / M3_PER_KM3
                for key in (
                    "external_supply_m3", "evaporation_m3", "sublimation_m3", "overflow_out_m3",
                    "storage_change_m3", "ice_storage_change_m3", "freezing_m3", "melting_m3",
                )
            },
            "normalized_lake_residual": abs(lake.residual_m3) / max(lake.external_supply_m3, 1.0),
            "normalized_land_residual": abs(land.residual_m3) / (land.rain_m3 + land.meltwater_m3),
        })
        if captured is not None:
            anchor_nodes = [int(np.argmin(np.linalg.norm(hydro.points - anchor, axis=1))) for anchor in anchors]
            result["steps"][-1]["loss_audit"] = summarize_loss_audit(captured, anchor_nodes)
        liquid_residual = (
            lake.external_supply_m3 - lake.evaporation_m3 - lake.sublimation_m3
            - lake.overflow_out_m3 - lake.storage_change_m3 - lake.freezing_m3 + lake.melting_m3
        )
        ice_residual = lake.freezing_m3 - lake.melting_m3 - lake.ice_storage_change_m3
        result["steps"][-1].update(
            normalized_liquid_residual=abs(liquid_residual) / max(lake.external_supply_m3, 1.0),
            normalized_ice_residual=abs(ice_residual) / max(lake.freezing_m3 + lake.melting_m3, 1.0),
        )

    def close_window(_steps: int) -> dict:
        result["summary"] = summarize(result["steps"])
        return {"completed_steps": len(result["steps"])}

    run_replay(world, args, result, snapshot, after_step, close_window)


if __name__ == "__main__":
    main()

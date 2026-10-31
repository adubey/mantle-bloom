#!/usr/bin/env python3
"""Follow bilateral slab pull plus a native-speed geometric midline upwelling to breakup.

Usage:
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
        backend/.venv/bin/python bin/debug/probe_gate3b_slab_pull_to_breakup.py \
        --out analysis/issue699/slab-pull-to-breakup-seed2024-density4.json

The world is regenerated from seed 2024 and the 2 Myr startup used in the preceding
bilateral screen. A central line of convection centers is scaled to the starting native
mantle-flow p95 and remains fixed in the mantle frame. The diagnostic measures the
incremental basal-drag force from those centers. It does not add a dynamic 3D slab model;
the production slab-pull law uses a fixed effective sunk length.
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components

from app import cratons, geometry, lithosphere, mantle, orogeny, torque, world as world_mod
from app.sparse_quad_patch import parent_cell_keys
import probe_gate3b_component_shell as shell_audit
import probe_gate3b_midplate_current as current_audit
import probe_gate3b_slab_pull as slab_screen
import replay_gate2_capacity as capacity_audit

SEED = 2024
DENSITY = 4.0
STARTUP_YEARS = 2_000_000.0
STEP_YEARS = 2_000_000.0
MAX_HISTORY_YEARS = 1_000_000_000.0
SLICE_COUNT = 10
SHELL_CAP = 512
SHELL_INTERVAL_YEARS = 10_000_000.0
RHO_MANTLE = 3300.0


def _lineage_indices(plate, initial_keys):
    """All current material nodes descended from the starting target plate cells."""
    current = plate.surface_nodes().node_ids[:, 0].view(np.int64)
    selected = {int(key) for key in initial_keys}
    selected_ancestors = set(selected)
    for key in tuple(selected):
        ancestor = key
        while True:
            parent = int(parent_cell_keys(np.array([ancestor], dtype=np.int64))[0])
            if parent < 0:
                break
            selected_ancestors.add(parent)
            ancestor = parent
    belongs = np.zeros(len(current), dtype=bool)
    for i, raw_key in enumerate(current):
        key = int(raw_key)
        if key in selected_ancestors:
            belongs[i] = True
            continue
        ancestor = key
        while True:
            parent = int(parent_cell_keys(np.array([ancestor], dtype=np.int64))[0])
            if parent < 0:
                break
            if parent in selected:
                belongs[i] = True
                break
            ancestor = parent
    return np.flatnonzero(belongs)


def _bins_for_nodes(plate, target_indices, key_to_bin, initial_keys, initial_points, initial_bins):
    keys = plate.surface_nodes().node_ids[:, 0].view(np.int64)
    result = np.full(len(target_indices), -1, dtype=np.int64)
    initial_key_set = set(key_to_bin)
    for row, idx in enumerate(target_indices):
        key = int(keys[idx])
        if key in key_to_bin:
            result[row] = key_to_bin[key]
            continue
        ancestor = key
        while True:
            parent = int(parent_cell_keys(np.array([ancestor], dtype=np.int64))[0])
            if parent < 0:
                break
            if parent in initial_key_set:
                result[row] = key_to_bin[parent]
                break
            ancestor = parent
    unresolved = result < 0
    if unresolved.any():
        points = geometry.normalize(plate.surface_nodes().world_xyz[target_indices[unresolved]])
        # Only new/coarsened cells use this geometric fallback; stable cells follow exact lineage.
        nearest = cKDTree(initial_points).query(points, k=1)[1]
        result[unresolved] = initial_bins[nearest]
    return result


def _force_components(world, plate, source_centers):
    field, inputs, subducting = slab_screen._field_and_inputs(world, plate)
    points = geometry.normalize(field.points_xyz)
    source_velocity_m_s = (
        mantle.flow_at(points, source_centers)
        * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    )
    midline_force = (
        torque.ASTHENOSPHERE_VISCOSITY_PA_S / torque.SHEAR_ZONE_THICKNESS_M
        * field.areas_m2[:, None] * source_velocity_m_s
    )
    native_basal = field.basal_drag_n - midline_force
    return field, inputs, subducting, midline_force, native_basal


def _weighted_stats(values, weights):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not keep.any():
        return {"mean": None, "median": None, "p90": None, "min": None, "max": None}
    vals = values[keep]
    w = weights[keep]
    order = np.argsort(vals)
    sorted_vals = vals[order]
    cumulative = np.cumsum(w[order]) / w.sum()
    return {
        "mean": float(np.average(vals, weights=w)),
        "median": float(sorted_vals[np.searchsorted(cumulative, 0.5)]),
        "p90": float(sorted_vals[np.searchsorted(cumulative, 0.9)]),
        "min": float(vals.min()), "max": float(vals.max()),
    }


def _slice_snapshot(world, plate, target_indices, target_bins, source_centers, native_centers):
    field, inputs, subducting, midline_force, native_basal = _force_components(world, plate, source_centers)
    spacing = shell_audit.line_spacing_rad(world.node_density)
    front_length = torque.boundary_node_length_m(inputs, spacing)
    hc = plate.collect("crustal_thickness_m")[target_indices]
    hm = plate.collect("mantle_lithosphere_thickness_m")[target_indices]
    areas = plate.node_areas_m2()[target_indices]
    crust_density = lithosphere.node_crust_density(
        plate.collect("crust_type_code")[target_indices], plate.crust_type
    )
    mass = areas * (hc * crust_density + hm * RHO_MANTLE)
    lag = plate.collect("moho_thermal_lag_c")[target_indices]
    temp = orogeny.moho_temperature_c(hc, hm, lag)
    subducting_target = subducting[target_indices]
    thermal = np.asarray(temp)
    forces = {
        "slab_pull": field.slab_pull_n[target_indices],
        "slab_drag": field.slab_drag_n[target_indices],
        "midline_upwelling_basal_drag": midline_force[target_indices],
        "native_basal_drag": native_basal[target_indices],
        "other_production_forces": (field.ridge_push_n + field.contact_drag_n)[target_indices],
        "combined_production": field.total_n[target_indices],
    }
    slices = []
    for bin_id in range(SLICE_COUNT):
        local = target_bins == bin_id
        if not local.any():
            slices.append({"slice": bin_id, "node_count": 0, "mass_kg": 0.0})
            continue
        local_area = areas[local]
        local_hc = hc[local]
        local_hm = hm[local]
        local_mass = mass[local]
        local_temp = thermal[local]
        local_sub = subducting_target[local]
        current_indices = target_indices[local]
        column = {
            "slice": bin_id,
            "node_count": int(local.sum()),
            "area_km2": float(local_area.sum() / 1e6),
            "mass_kg": float(local_mass.sum()),
            "crust_volume_km3": float(np.sum(local_area * local_hc) / 1e9),
            "mantle_lithosphere_volume_km3": float(np.sum(local_area * local_hm) / 1e9),
            "crust_thickness_km": _weighted_stats(local_hc / 1000.0, local_area),
            "mantle_lithosphere_thickness_km": _weighted_stats(local_hm / 1000.0, local_area),
            "total_lithosphere_thickness_km": _weighted_stats((local_hc + local_hm) / 1000.0, local_area),
            "moho_temperature_c": _weighted_stats(local_temp, local_area),
            "active_subducting_nodes": int(local_sub.sum()),
            "active_slab_front_length_km": float(front_length[current_indices[local_sub]].sum() / 1000.0),
            "subducting_column_thickness_km": _weighted_stats(
                (local_hc[local_sub] + local_hm[local_sub]) / 1000.0,
                local_area[local_sub],
            ),
            "loads": {},
        }
        for name, vectors in forces.items():
            selected = np.asarray(vectors)[local]
            resultant = selected.sum(axis=0)
            mass_total = float(local_mass.sum())
            column["loads"][name] = {
                "resultant_n": resultant.tolist(),
                "resultant_magnitude_n": float(np.linalg.norm(resultant)),
                "force_per_mass_n_per_kg": float(np.linalg.norm(resultant) / mass_total) if mass_total else None,
            }
        slices.append(column)
    return {
        "elapsed_myr_after_startup": float((world.elapsed_years - STARTUP_YEARS) / 1e6),
        "plate_id": int(plate.plate_id),
        "target_material_area_km2": float(areas.sum() / 1e6),
        "full_containing_plate_area_km2": float(plate.node_areas_m2().sum() / 1e6),
        "active_slab_front_length_km": float(front_length[target_indices[subducting_target]].sum() / 1000.0),
        "effective_slab_sink_length_km": float(torque.SUBDUCTION_LSINK_M / 1000.0),
        "target_native_mantle_flow_p95_cm_per_yr": float(mantle.rad_per_yr_to_cm_per_yr(
            np.percentile(np.linalg.norm(mantle.flow_at(geometry.normalize(field.points_xyz[target_indices]), native_centers), axis=1), 95.0)
        )),
        "target_midline_source_flow_p95_cm_per_yr": float(mantle.rad_per_yr_to_cm_per_yr(
            np.percentile(np.linalg.norm(mantle.flow_at(geometry.normalize(field.points_xyz[target_indices]), source_centers), axis=1), 95.0)
        )),
        "slices": slices,
    }, field, inputs, subducting, midline_force, native_basal


def _shell_snapshot(world, plate, target_indices, field):
    result = shell_audit._plate_shell(
        plate, target_indices, field.points_xyz, field.total_n, SHELL_CAP, include_demand=True
    )
    demand = np.asarray(result.pop("component_demand_n_per_m"), dtype=float)
    graph = current_audit._component_graph(plate, target_indices)
    area = plate.node_areas_m2()[target_indices]
    capacities = capacity_audit._physical_capacities()["continental"]
    rows = {}
    for name, cap in zip(("p10", "median", "p90"), capacities):
        over = demand >= cap
        if over.any():
            count, labels = connected_components(graph[over][:, over], directed=False)
            component_areas = np.bincount(labels, weights=area[over], minlength=count)
            max_fraction = float(component_areas.max() / max(float(area.sum()), 1.0))
        else:
            max_fraction = 0.0
        rows[name] = {
            "capacity_n_per_m": float(cap),
            "max_connected_overcapacity_target_area_fraction": max_fraction,
        }
    result["capacity_cases"] = rows
    finite = demand[np.isfinite(demand)]
    result["target_demand_quantiles_n_per_m"] = np.quantile(finite, [0.1, 0.5, 0.9]).tolist()
    return result


def _write(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run(out_path: Path):
    started = time.perf_counter()
    world = world_mod.generate_world(SEED, node_density=DENSITY, cratons_first=True, synthetic_relief=False)
    world.simulate_climate_biomes = False
    cratons.ensure_initialized(world)
    world_mod.step_world(world, STARTUP_YEARS)
    target_plate, candidates = slab_screen._select(world)
    if target_plate is None or target_plate.plate_id != 0:
        raise RuntimeError("seed-2024 bilateral candidate did not reproduce as plate 0")
    initial = target_plate.surface_nodes()
    initial_keys = initial.node_ids[:, 0].view(np.int64).copy()
    initial_points, _areas, _center, _axis = slab_screen._principal_axis(target_plate)
    candidate = slab_screen._bilateral_candidate(world, target_plate)
    coord = candidate["coord"]
    slice_edges = np.linspace(float(coord.min()), float(coord.max()), SLICE_COUNT + 1)
    initial_bins = np.clip(np.digitize(coord, slice_edges[1:-1]), 0, SLICE_COUNT - 1)
    key_to_bin = {int(key): int(bin_id) for key, bin_id in zip(initial_keys, initial_bins)}

    native_centers = list(world.mantle_centers)
    native_flow = mantle.flow_at(candidate["points"], native_centers)
    native_p95_cm_yr = float(mantle.rad_per_yr_to_cm_per_yr(
        np.percentile(np.linalg.norm(native_flow, axis=1), 95.0)
    ))
    source_geometry = current_audit._source_geometry(target_plate, shell_audit.line_spacing_rad(DENSITY))
    source_centers, source_scale = current_audit._scaled_source_centers(
        source_geometry["source_centers_unit_strength"], candidate["points"], native_p95_cm_yr
    )
    world.mantle_centers = native_centers + source_centers
    geometry_record = {key: value for key, value in source_geometry.items() if key != "source_centers_unit_strength"}

    report = {
        "title": "RFC 0010 Gate 3b slab pull plus native-speed midline upwelling to breakup",
        "seed": SEED, "node_density": DENSITY,
        "startup_myr": STARTUP_YEARS / 1e6,
        "step_myr": STEP_YEARS / 1e6,
        "max_history_myr_after_startup": MAX_HISTORY_YEARS / 1e6,
        "selected_initial_plate_id": int(target_plate.plate_id),
        "selected_initial_target_area_km2": float(candidate["area_m2"] / 1e6),
        "bilateral_precondition": candidate["side_groups"],
        "selection_candidates": candidates,
        "midline_geometry": geometry_record,
        "midline_source_target_p95_cm_per_yr": native_p95_cm_yr,
        "midline_source_strength_scale": float(source_scale),
        "fixed_slab_sink_length_km": float(torque.SUBDUCTION_LSINK_M / 1000.0),
        "slice_count": SLICE_COUNT,
        "slice_projection_edges": slice_edges.tolist(),
        "shell_interval_myr": SHELL_INTERVAL_YEARS / 1e6,
        "breakup_definition": "The selected starting material lineage occupies more than one production plate after a completed plate split. Failed rifts do not count.",
        "checkpoints": [],
        "event_history": [],
        "breakup_time_myr_after_startup": None,
        "breakup_event": None,
        "stop_reason": None,
        "elapsed_wall_seconds": None,
    }
    elapsed_steps = int(round(MAX_HISTORY_YEARS / STEP_YEARS))
    last_event_elapsed_years = world.elapsed_years
    last_shell_step = -1
    target_plate_id = int(target_plate.plate_id)
    for step in range(1, elapsed_steps + 1):
        world_mod.step_world(world, STEP_YEARS)
        elapsed_myr = float((world.elapsed_years - STARTUP_YEARS) / 1e6)
        new_events = [
            (event_year, message) for event_year, message in world.events
            if last_event_elapsed_years < event_year <= world.elapsed_years
        ]
        last_event_elapsed_years = world.elapsed_years
        for event_year, message in new_events:
            report["event_history"].append({
                "elapsed_myr_after_startup": float((event_year - STARTUP_YEARS) / 1e6),
                "message": message,
            })
        split_for_target = None
        merged_target_this_step = False
        for _event_year, message in new_events:
            split = re.match(r"Plate (\d+) split into plates (\d+) and (\d+)\.", message)
            if split and int(split.group(1)) == target_plate_id:
                split_for_target = (int(split.group(2)), int(split.group(3)), message)
                break
            merge = re.match(
                r"Plates (\d+) and (\d+) (?:collided and merged|had overlapped so long they fused) into plate (\d+)",
                message,
            )
            if merge and target_plate_id in (int(merge.group(1)), int(merge.group(2))):
                target_plate_id = int(merge.group(3))
                merged_target_this_step = True
                report.setdefault("target_topology_transitions", []).append({
                    "elapsed_myr_after_startup": elapsed_myr,
                    "new_tracked_plate_id": target_plate_id,
                    "message": message,
                })
        if split_for_target is not None:
            report["breakup_time_myr_after_startup"] = elapsed_myr
            daughter_rows = []
            daughter_areas = {}
            for daughter_id in split_for_target[:2]:
                daughter = next((p for p in world.plates if p.plate_id == daughter_id), None)
                if daughter is None:
                    continue
                daughter_indices = _lineage_indices(daughter, initial_keys)
                if not len(daughter_indices):
                    daughter_indices = np.arange(daughter.node_count(), dtype=np.int64)
                daughter_bins = _bins_for_nodes(
                    daughter, daughter_indices, key_to_bin, initial_keys, initial_points, initial_bins
                )
                daughter_row, daughter_field, _i, _s, _m, _n = _slice_snapshot(
                    world, daughter, daughter_indices, daughter_bins, source_centers, native_centers
                )
                try:
                    daughter_row["shell"] = _shell_snapshot(world, daughter, daughter_indices, daughter_field)
                except Exception as exc:
                    daughter_row["shell_error"] = f"{type(exc).__name__}: {exc}"
                daughter_rows.append(daughter_row)
                daughter_areas[str(daughter_id)] = float(daughter.node_areas_m2()[daughter_indices].sum() / 1e6)
            report["breakup_event"] = {
                "elapsed_myr_after_startup": elapsed_myr,
                "daughter_plate_ids": list(split_for_target[:2]),
                "target_lineage_area_km2_by_daughter": daughter_areas,
                "daughter_load_and_slice_snapshots": daughter_rows,
                "logged_split_event": split_for_target[2],
            }
            report["stop_reason"] = "tracked plate completed a production split"
            break
        plate = next((p for p in world.plates if p.plate_id == target_plate_id), None)
        if plate is None:
            report["stop_reason"] = "selected material lineage disappeared before a completed split"
            break
        target_indices = _lineage_indices(plate, initial_keys)
        if not len(target_indices) and merged_target_this_step:
            # Merges can resample cell IDs. Continue on the resulting plate system and
            # record the geometry reset instead of silently losing the long history.
            current_nodes = plate.surface_nodes()
            initial_keys = current_nodes.node_ids[:, 0].view(np.int64).copy()
            initial_points, _area, _center, axis = slab_screen._principal_axis(plate)
            coord = initial_points @ axis
            slice_edges = np.linspace(float(coord.min()), float(coord.max()), SLICE_COUNT + 1)
            initial_bins = np.clip(np.digitize(coord, slice_edges[1:-1]), 0, SLICE_COUNT - 1)
            key_to_bin = {int(key): int(bin_id) for key, bin_id in zip(initial_keys, initial_bins)}
            target_indices = np.arange(plate.node_count(), dtype=np.int64)
            report.setdefault("slice_geometry_rebases", []).append({
                "elapsed_myr_after_startup": elapsed_myr,
                "plate_id": target_plate_id,
                "reason": "tracked cell lineage unavailable after topology resampling",
                "area_km2": float(plate.node_areas_m2().sum() / 1e6),
            })
        elif not len(target_indices):
            report["stop_reason"] = "selected material lineage disappeared before a completed split"
            break
        target_bins = _bins_for_nodes(plate, target_indices, key_to_bin, initial_keys, initial_points, initial_bins)
        row, field, _inputs, _subducting, _midline, _native = _slice_snapshot(
            world, plate, target_indices, target_bins, source_centers, native_centers
        )
        row["elapsed_myr_after_startup"] = elapsed_myr
        for event_year, message in new_events:
            if "failed" in message.lower() or "aulacogen" in message.lower():
                row.setdefault("failed_rift_events_this_step", []).append(message)
        if step % max(1, int(round(SHELL_INTERVAL_YEARS / STEP_YEARS))) == 0:
            try:
                row["shell"] = _shell_snapshot(world, plate, target_indices, field)
                last_shell_step = step
            except Exception as exc:
                row["shell_error"] = f"{type(exc).__name__}: {exc}"
                last_shell_step = step
        report["checkpoints"].append(row)
        report["elapsed_wall_seconds"] = time.perf_counter() - started
        if step % 5 == 0 or step == 1:
            _write(out_path, report)
        if step % 25 == 0:
            print(json.dumps({
                "elapsed_myr_after_startup": elapsed_myr,
                "steps": step, "plate_id": int(plate.plate_id),
                "target_area_km2": row["target_material_area_km2"],
                "active_slab_front_length_km": row["active_slab_front_length_km"],
                "midline_p95_cm_per_yr": row["target_native_mantle_flow_p95_cm_per_yr"],
                "wall_seconds": report["elapsed_wall_seconds"],
            }, sort_keys=True), flush=True)
    else:
        report["stop_reason"] = "no completed split or target-lineage loss by 1 Gyr"
        report["right_censored_at_myr_after_startup"] = MAX_HISTORY_YEARS / 1e6

    report["elapsed_wall_seconds"] = time.perf_counter() - started
    report["last_shell_step"] = last_shell_step
    _write(out_path, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print(json.dumps({
        "out": str(args.out),
        "stop_reason": report["stop_reason"],
        "breakup_time_myr_after_startup": report["breakup_time_myr_after_startup"],
        "right_censored_at_myr_after_startup": report.get("right_censored_at_myr_after_startup"),
        "checkpoints": len(report["checkpoints"]),
        "elapsed_wall_seconds": report["elapsed_wall_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

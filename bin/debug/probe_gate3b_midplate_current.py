#!/usr/bin/env python3
"""Sweep a symmetric midplate convection source through Gate 3b's shell diagnostic.

Usage:
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
        backend/.venv/bin/python bin/debug/probe_gate3b_midplate_current.py \
        --out analysis/issue699/midplate-current-seed11-density4.json

The source is a diagnostic line of upwelling convection centers with paired, symmetric
downwelling rows. It is superposed on native mantle centers, then applied through production
basal drag. Each amplitude starts from the same generated world. The reduced shell field is
diagnostic only; this script does not mutate crust, accumulate production cracks, or split a
plate.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from app import cratons, elevation_lines, geometry, lithosphere, mantle, surface_fields, world as world_mod
from app.sparse_quad_patch import parent_cell_keys
import probe_gate3b_component_shell as shell_audit
import replay_gate2_capacity as capacity_audit

SEED = 11
DENSITY = 4.0
STEP_YEARS = 2_000_000.0
DURATION_YEARS = 30_000_000.0
SOURCE_SPEEDS_CM_PER_YR = (0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 40.0, 48.0, 64.0, 96.0, 128.0)
HISTORY_CAP = 512
EVENT_AREA_FRACTION = 0.25
EVENT_MIN_DURATION_YEARS = 10_000_000.0


def _source_geometry(plate, node_spacing_rad: float) -> dict:
    points = geometry.normalize(plate.surface_nodes().world_xyz)
    areas = plate.node_areas_m2()
    center = geometry.normalize(np.sum(points * areas[:, None], axis=0))
    tangent = points - (points @ center)[:, None] * center
    covariance = np.einsum("n,ni,nj->ij", areas, tangent, tangent) / max(float(areas.sum()), 1.0)
    _eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    axis = eigenvectors[:, -1] - center * np.dot(eigenvectors[:, -1], center)
    axis = geometry.normalize(axis)
    pivot = int(np.argmax(np.abs(axis)))
    if axis[pivot] < 0.0:
        axis = -axis
    normal = geometry.normalize(np.cross(center, axis))

    equivalent_radius_m = math.sqrt(float(areas.sum()) / math.pi)
    equivalent_radius_rad = equivalent_radius_m / lithosphere.PLANET_RADIUS_M
    falloff_rad = 0.25 * equivalent_radius_rad
    sink_offset_rad = 0.5 * equivalent_radius_rad
    max_half_span = min(1.45, max(2.0 * equivalent_radius_rad, 4.0 * node_spacing_rad))
    sample_count = int(math.ceil(2.0 * max_half_span / max(node_spacing_rad / 2.0, 1e-8)))
    offsets = (np.arange(-sample_count, sample_count + 1) * (node_spacing_rad / 2.0))
    line_samples = np.cos(offsets)[:, None] * center + np.sin(offsets)[:, None] * axis
    inside = np.asarray(plate.contains_batch(line_samples), dtype=bool)
    zero = sample_count
    if not inside[zero]:
        raise ValueError("area-weighted plate center is outside the plate; no central source chord")
    left = zero
    while left > 0 and inside[left - 1]:
        left -= 1
    right = zero
    while right + 1 < len(inside) and inside[right + 1]:
        right += 1
    if right - left < 4:
        raise ValueError("geometry-only midplate chord is too short for a line source")
    chord_offsets = offsets[left:right + 1]
    source_count = min(32, max(5, int(math.ceil((chord_offsets[-1] - chord_offsets[0]) / max(falloff_rad / 2.0, 1e-8))) + 1))
    source_offsets = np.linspace(chord_offsets[0], chord_offsets[-1], source_count)
    source_points = np.cos(source_offsets)[:, None] * center + np.sin(source_offsets)[:, None] * axis
    sink_points = []
    for sign in (-1.0, 1.0):
        sink_points.append(
            math.cos(sink_offset_rad) * source_points
            + sign * math.sin(sink_offset_rad) * normal[None, :]
        )
    sink_points = np.concatenate(sink_points, axis=0)
    base_centers = [
        mantle.ConvectionCenter(position=point, strength=1.0, falloff=falloff_rad)
        for point in source_points
    ]
    base_centers.extend(
        mantle.ConvectionCenter(position=point, strength=-0.5, falloff=falloff_rad)
        for point in sink_points
    )
    return {
        "center_xyz": center.tolist(),
        "principal_axis_xyz": axis.tolist(),
        "cross_axis_xyz": normal.tolist(),
        "chord_start_offset_rad": float(chord_offsets[0]),
        "chord_end_offset_rad": float(chord_offsets[-1]),
        "chord_length_km": float((chord_offsets[-1] - chord_offsets[0]) * lithosphere.PLANET_RADIUS_M / 1000.0),
        "equivalent_radius_km": equivalent_radius_m / 1000.0,
        "source_falloff_km": falloff_rad * lithosphere.PLANET_RADIUS_M / 1000.0,
        "sink_offset_km": sink_offset_rad * lithosphere.PLANET_RADIUS_M / 1000.0,
        "upwelling_center_count": len(source_points),
        "downwelling_center_count": len(sink_points),
        "source_centers_unit_strength": base_centers,
    }


def _scaled_source_centers(base_centers, points, target_speed_cm_per_yr: float):
    if target_speed_cm_per_yr == 0.0:
        return [], 0.0
    unit_speed = np.linalg.norm(mantle.flow_at(points, base_centers), axis=1)
    p95 = float(np.percentile(unit_speed, 95.0))
    if p95 <= 0.0:
        raise ValueError("midplate source field has zero p95 speed")
    target = mantle.cm_per_yr_to_rad_per_yr(target_speed_cm_per_yr)
    scale = target / p95
    centers = [
        mantle.ConvectionCenter(c.position.copy(), c.strength * scale, c.falloff)
        for c in base_centers
    ]
    return centers, scale


def _component_graph(plate, indices):
    count = plate.node_count()
    adjacency = plate.adjacency()
    rows = np.repeat(np.arange(count, dtype=np.int64), np.diff(adjacency.offsets))
    cols = adjacency.neighbours.astype(np.int64, copy=False)
    included = np.zeros(count, dtype=bool)
    included[indices] = True
    keep = included[rows] & included[cols]
    local = np.full(count, -1, dtype=np.int64)
    local[indices] = np.arange(len(indices), dtype=np.int64)
    graph = csr_matrix(
        (np.ones(int(keep.sum())), (local[rows[keep]], local[cols[keep]])),
        shape=(len(indices), len(indices)),
    )
    return graph


def _component_indices_by_lineage(plate, initial_component_keys):
    """Follow selected initial quad cells through refinement and coarsening."""
    nodes = plate.surface_nodes()
    current = nodes.node_ids[:, 0].view(np.int64)
    selected = set(int(key) for key in initial_component_keys)
    selected_ancestors = set(selected)
    for key in tuple(selected):
        ancestor = int(key)
        while True:
            parent = int(parent_cell_keys(np.array([ancestor], dtype=np.int64))[0])
            if parent < 0:
                break
            selected_ancestors.add(parent)
            ancestor = parent
    belongs = np.zeros(len(current), dtype=bool)
    for index, raw_key in enumerate(current):
        key = int(raw_key)
        if key in selected_ancestors:
            belongs[index] = True
            continue
        ancestor = key
        while True:
            parent = int(parent_cell_keys(np.array([ancestor], dtype=np.int64))[0])
            if parent < 0:
                break
            if parent in selected:
                belongs[index] = True
                break
            ancestor = parent
    continental = elevation_lines.effective_is_continental_from_codes(
        plate.collect("crust_type_code"), plate.crust_type == "continental"
    )
    return np.flatnonzero(belongs & continental)


def _compact_shell_metrics(row: dict) -> dict:
    keys = (
        "target_representatives", "representative_count", "shell_node_count", "mesh_triangles",
        "production_plate_area_coverage", "force_transfer_residual_fraction",
        "torque_transfer_residual_fraction", "equilibrium_residual_fraction",
        "shell_energy_residual_fraction", "accounting_checks", "accounting_pass",
        "continental_component_demand_quantiles_n_per_m", "solver_stop_code",
    )
    return {key: row[key] for key in keys if key in row}


def _compact_checkpoint(point: dict) -> dict:
    if "capacity_cases" not in point:
        result = {key: value for key, value in point.items() if key not in ("accounting_1024", "accounting_512")}
        if "accounting_1024" in point:
            result["initial_mesh_accounting"] = {
                cap: _compact_shell_metrics(point[f"accounting_{cap}"]["mesh_accounting"])
                for cap in (1024, 512)
            }
            result["initial_basal_drag_energy_residual_fraction"] = {
                cap: point[f"accounting_{cap}"]["basal_drag_energy_residual_fraction"]
                for cap in (1024, 512)
            }
        return result
    keys = (
        "elapsed_myr", "source_speed_p95_cm_per_yr", "source_direction_coherence_across_line_area_weighted",
        "source_force_resultant_n", "source_torque_resultant_nm", "production_force_resultant_n",
        "production_torque_resultant_nm", "basal_drag_energy_residual_fraction", "basal_drag_energy_pass",
        "accounting_pass", "history_stopped_for_accounting", "tracked_plate_missing", "tracked_component_missing",
    )
    result = {key: point[key] for key in keys if key in point}
    result["mesh_accounting"] = _compact_shell_metrics(point["mesh_accounting"])
    result["capacity_cases"] = point["capacity_cases"]
    return result


def _compact_arm(arm: dict) -> dict:
    result = {key: value for key, value in arm.items() if key != "checkpoints"}
    result["checkpoints"] = [_compact_checkpoint(point) for point in arm.get("checkpoints", [])]
    return result


def _max_overcapacity_component(demand, capacity, areas, graph):
    over = np.asarray(demand) >= float(capacity)
    if not np.any(over):
        return 0.0, over
    subgraph = graph[over][:, over]
    count, labels = connected_components(subgraph, directed=False)
    component_areas = np.bincount(labels, weights=areas[over], minlength=count)
    return float(component_areas.max() / max(float(areas.sum()), 1.0)), over


def _evaluate_snapshot(world, plate, component_indices, graph, source_centers, time_years, cap_count):
    field = shell_audit._production_force_field(world, plate)
    expected = geometry.normalize(plate.surface_nodes().world_xyz)
    if len(field.points_xyz) != len(expected) or np.max(np.linalg.norm(expected - geometry.normalize(field.points_xyz), axis=1)) > 1e-10:
        raise ValueError("production force ordering does not match parent plate topology")
    source_only = np.zeros_like(field.total_n)
    if source_centers:
        points = geometry.normalize(field.points_xyz)
        source_velocity = mantle.flow_at(points, source_centers) * lithosphere.PLANET_RADIUS_M / shell_audit.torque.SECONDS_PER_YEAR
        coefficient = shell_audit.torque.ASTHENOSPHERE_VISCOSITY_PA_S / shell_audit.torque.SHEAR_ZONE_THICKNESS_M
        source_only = coefficient * source_velocity * field.areas_m2[:, None]
        # Source attribution is the added-current traction at the current plate velocity.
        # The total production basal field still includes native drag and plate response.
    shell_row = shell_audit._plate_shell(
        plate, component_indices, field.points_xyz, field.total_n, cap_count, include_demand=True
    )
    demand = np.asarray(shell_row.pop("component_demand_n_per_m"), dtype=float)
    areas = plate.node_areas_m2()[component_indices]
    capacities = capacity_audit._physical_capacities()["continental"]
    capacity_cases = {}
    for name, value in zip(("p10", "median", "p90"), capacities):
        max_fraction, over = _max_overcapacity_component(demand, value, areas, graph)
        craton = plate.collect("craton_crust_m")[component_indices] > 0.0
        belt = plate.collect("belt_formed_years")[component_indices] != surface_fields.BELT_UNMARKED_YEARS
        area = max(float(areas[over & craton].sum()), 0.0)
        capacity_cases[name] = {
            "capacity_n_per_m": float(value),
            "max_connected_overcapacity_area_fraction": max_fraction,
            "screen_pass_at_checkpoint": bool(max_fraction >= EVENT_AREA_FRACTION),
            "overcapacity_area_fraction_craton": area / max(float(areas.sum()), 1.0),
            "overcapacity_area_fraction_belt": float(areas[over & belt].sum() / max(float(areas.sum()), 1.0)),
            "overcapacity_nodes": int(over.sum()),
        }
    source_velocity = mantle.flow_at(geometry.normalize(field.points_xyz), source_centers) if source_centers else np.zeros_like(field.points_xyz)
    cross_axis = np.asarray(_CURRENT_GEOMETRY["cross_axis_xyz"], dtype=float)
    cross_tangent = cross_axis[None, :] - (geometry.normalize(field.points_xyz) @ cross_axis)[:, None] * geometry.normalize(field.points_xyz)
    cross_tangent = geometry.normalize(cross_tangent)
    velocity_norm = np.linalg.norm(source_velocity, axis=1)
    coherence = np.divide(
        np.abs(np.sum(source_velocity * cross_tangent, axis=1)), velocity_norm,
        out=np.zeros_like(velocity_norm), where=velocity_norm > 0.0,
    )
    basal_account = shell_audit._drag_energy(world, plate, np.arange(plate.node_count()), field)
    total_force = np.asarray(field.total_n)
    source_force = source_only
    return {
        "elapsed_myr": float(time_years / 1e6),
        "mesh_accounting": shell_row,
        "capacity_cases": capacity_cases,
        "source_speed_p95_cm_per_yr": float(mantle.rad_per_yr_to_cm_per_yr(np.percentile(velocity_norm, 95.0))),
        "source_force_resultant_n": source_force.sum(axis=0).tolist(),
        "source_torque_resultant_nm": np.cross(field.points_xyz * lithosphere.PLANET_RADIUS_M, source_force).sum(axis=0).tolist(),
        "production_force_resultant_n": total_force.sum(axis=0).tolist(),
        "production_torque_resultant_nm": np.cross(field.points_xyz * lithosphere.PLANET_RADIUS_M, total_force).sum(axis=0).tolist(),
        "basal_drag_energy_residual_fraction": basal_account["energy_residual_fraction"],
        "basal_drag_energy_pass": bool(basal_account["energy_residual_fraction"] <= 0.05),
        "source_direction_coherence_across_line_area_weighted": float(np.average(coherence, weights=field.areas_m2)),
    }


def _run_arm(speed_cm_per_yr: float, step_years: float, line_geometry: dict, initial_component_keys) -> dict:
    world = world_mod.generate_world(SEED, node_density=DENSITY, cratons_first=True, synthetic_relief=False)
    world.simulate_climate_biomes = False
    cratons.ensure_initialized(world)
    _component_area, initial_plate, component_indices = shell_audit._largest_continental_component(world)
    initial_keys = initial_plate.surface_nodes().node_ids[component_indices, 0].view(np.int64)
    if initial_plate.plate_id != _INITIAL_PLATE_ID or not np.array_equal(initial_keys, initial_component_keys):
        raise ValueError("generated initial plate/component did not reproduce the frozen selection")
    points = geometry.normalize(initial_plate.surface_nodes().world_xyz)
    source_centers, source_scale = _scaled_source_centers(
        line_geometry["source_centers_unit_strength"], points, speed_cm_per_yr
    )
    native_centers = list(world.mantle_centers)
    world.mantle_centers = native_centers + source_centers
    graph = _component_graph(initial_plate, component_indices)
    duration_by_capacity = {name: 0.0 for name in ("p10", "median", "p90")}
    prior_pass = {name: False for name in duration_by_capacity}
    event_onset = {name: None for name in duration_by_capacity}
    persistent_event_onset = {name: None for name in duration_by_capacity}
    checkpoints = []
    initial = _evaluate_snapshot(world, initial_plate, component_indices, graph, source_centers, 0.0, 1024)
    initial_512 = _evaluate_snapshot(world, initial_plate, component_indices, graph, source_centers, 0.0, HISTORY_CAP)
    initial_accounting = bool(
        initial["mesh_accounting"]["accounting_pass"]
        and initial_512["mesh_accounting"]["accounting_pass"]
        and initial["basal_drag_energy_pass"]
        and initial_512["basal_drag_energy_pass"]
    )
    checkpoints.append({"checkpoint": "initial", "accounting_1024": initial, "accounting_512": initial_512})
    if not initial_accounting:
        return {
            "source_speed_cm_per_yr": speed_cm_per_yr,
            "step_myr": step_years / 1e6,
            "eligible": False,
            "reason": "initial parent-shell or basal-drag accounting failed",
            "source_center_strength_scale": source_scale,
            "checkpoints": checkpoints,
            "event_onset_myr_by_capacity": event_onset,
        }

    steps = int(round(DURATION_YEARS / step_years))
    for _ in range(steps):
        world_mod.step_world(world, step_years)
        plate = next((item for item in world.plates if item.plate_id == _INITIAL_PLATE_ID), None)
        if plate is None:
            checkpoints.append({"elapsed_myr": float(world.elapsed_years / 1e6), "tracked_plate_missing": True})
            break
        component_indices = _component_indices_by_lineage(plate, initial_component_keys)
        if not len(component_indices):
            checkpoints.append({"elapsed_myr": float(world.elapsed_years / 1e6), "tracked_component_missing": True})
            break
        graph = _component_graph(plate, component_indices)
        try:
            snap = _evaluate_snapshot(world, plate, component_indices, graph, source_centers, world.elapsed_years, HISTORY_CAP)
        except Exception as exc:
            return {
                "source_speed_cm_per_yr": speed_cm_per_yr,
                "step_myr": step_years / 1e6,
                "eligible": False,
                "initial_accounting_pass": initial_accounting,
                "diagnostic_error": f"{type(exc).__name__}: {exc}",
                "failure_elapsed_myr": float(world.elapsed_years / 1e6),
                "persistent_event_onset_myr_by_capacity": persistent_event_onset,
                "persistent_event_by_capacity": {name: False for name in duration_by_capacity},
                "checkpoints": checkpoints,
            }
        snap["accounting_pass"] = bool(
            snap["mesh_accounting"]["accounting_pass"] and snap["basal_drag_energy_pass"]
        )
        if not snap["accounting_pass"]:
            snap["history_stopped_for_accounting"] = True
            checkpoints.append(snap)
            break
        for name, case in snap["capacity_cases"].items():
            now_passes = case["screen_pass_at_checkpoint"]
            if now_passes:
                if prior_pass[name]:
                    duration_by_capacity[name] += step_years
                else:
                    event_onset[name] = snap["elapsed_myr"]
                if duration_by_capacity[name] >= EVENT_MIN_DURATION_YEARS and event_onset[name] is not None:
                    case["persistent_event_seen"] = True
                    if persistent_event_onset[name] is None:
                        persistent_event_onset[name] = event_onset[name]
            else:
                duration_by_capacity[name] = 0.0
                event_onset[name] = None
            prior_pass[name] = now_passes
            case["current_consecutive_duration_myr"] = duration_by_capacity[name] / 1e6
        checkpoints.append(snap)
    return {
        "source_speed_cm_per_yr": speed_cm_per_yr,
        "step_myr": step_years / 1e6,
        "eligible": initial_accounting,
        "source_center_strength_scale": source_scale,
        "initial_accounting_pass": initial_accounting,
        "persistent_event_onset_myr_by_capacity": persistent_event_onset,
        "final_qualifying_streak_onset_myr_by_capacity": event_onset,
        "persistent_event_by_capacity": {
            name: any(
                point.get("capacity_cases", {}).get(name, {}).get("persistent_event_seen", False)
                for point in checkpoints
            ) for name in duration_by_capacity
        },
        "checkpoints": checkpoints,
    }


def run(out_path: Path) -> dict:
    global _CURRENT_GEOMETRY, _INITIAL_PLATE_ID
    started = time.perf_counter()
    reference = world_mod.generate_world(SEED, node_density=DENSITY, cratons_first=True, synthetic_relief=False)
    reference.simulate_climate_biomes = False
    cratons.ensure_initialized(reference)
    component_area, plate, component_indices = shell_audit._largest_continental_component(reference)
    _INITIAL_PLATE_ID = plate.plate_id
    initial_component_keys = plate.surface_nodes().node_ids[component_indices, 0].view(np.int64)
    _CURRENT_GEOMETRY = _source_geometry(plate, shell_audit.line_spacing_rad(DENSITY))
    geometry_json = {key: value for key, value in _CURRENT_GEOMETRY.items() if key != "source_centers_unit_strength"}
    arms = []
    if out_path.exists():
        try:
            previous = json.loads(out_path.read_text(encoding="utf-8"))
            arms = [
                row for row in previous.get("arms_2myr", [])
                if row.get("step_myr") == STEP_YEARS / 1e6
            ]
        except (OSError, json.JSONDecodeError):
            arms = []

    def write_partial():
        result = {
            "schema_version": 1,
            "title": "RFC 0010 Gate 3b midplate-current sweep and persistence screen",
            "hypothesis": "A symmetric midplate upwelling source can impose broad extensional load; increasing source speed may create a persistent over-capacity footprint.",
            "seed": SEED,
            "node_density": DENSITY,
            "initial_elapsed_myr": float(reference.elapsed_years / 1e6),
            "plate_id": int(_INITIAL_PLATE_ID),
            "selected_component_nodes": int(len(component_indices)),
            "selected_component_area_km2": float(component_area / 1e6),
            "source_geometry": geometry_json,
            "speed_sweep_cm_per_yr": list(SOURCE_SPEEDS_CM_PER_YR),
            "event_rule": {
                "overcapacity_ratio": "Gate 2 N_ext/N_cap >= 1",
                "connected_continental_area_fraction_min": EVENT_AREA_FRACTION,
                "duration_myr_min": EVENT_MIN_DURATION_YEARS / 1e6,
            },
            "arms_2myr": [_compact_arm(row) for row in arms if row.get("step_myr") == STEP_YEARS / 1e6],
            "arms_1myr": [_compact_arm(row) for row in arms if row.get("step_myr") == 1.0],
            "elapsed_wall_seconds": time.perf_counter() - started,
        }
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    completed_speeds = {row.get("source_speed_cm_per_yr") for row in arms}
    for speed in SOURCE_SPEEDS_CM_PER_YR:
        if speed in completed_speeds:
            continue
        try:
            arm = _run_arm(speed, STEP_YEARS, _CURRENT_GEOMETRY, initial_component_keys)
        except Exception as exc:
            arm = {
                "source_speed_cm_per_yr": speed,
                "step_myr": STEP_YEARS / 1e6,
                "eligible": False,
                "diagnostic_error": f"{type(exc).__name__}: {exc}",
                "persistent_event_by_capacity": {name: False for name in ("p10", "median", "p90")},
                "checkpoints": [],
            }
        arms.append(arm)
        write_partial()
    first_passing = next((
        index for index, arm in enumerate(arms)
        if arm.get("eligible") and any(arm.get("persistent_event_by_capacity", {}).values())
    ), None)
    if first_passing is not None:
        selected_indices = [first_passing]
        if first_passing > 0:
            selected_indices.insert(0, first_passing - 1)
        for index in selected_indices:
            arm = _run_arm(SOURCE_SPEEDS_CM_PER_YR[index], 1_000_000.0, _CURRENT_GEOMETRY, initial_component_keys)
            arms.append(arm)
            write_partial()
    write_partial()
    report = json.loads(out_path.read_text(encoding="utf-8"))
    report["interpretation"] = (
        "Diagnostic source/load and persistent capacity-screen only. Over-capacity is not a local crack, "
        "and no crack accumulation or plate split is implemented by this probe."
    )
    report["elapsed_wall_seconds"] = time.perf_counter() - started
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print(json.dumps({
        "out": str(args.out),
        "arms": len(report["arms_2myr"]),
        "persistent_events": {
            str(arm.get("source_speed_cm_per_yr")): arm.get("persistent_event_by_capacity")
            for arm in report["arms_2myr"]
        },
        "elapsed_wall_seconds": report["elapsed_wall_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Measure bilateral slab-pull loading and breakup time for an eligible generated plate.

Usage:
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
        backend/.venv/bin/python bin/debug/probe_gate3b_slab_pull.py \
        --out analysis/issue699/slab-pull-bilateral-screen.json

Uses native generated mantle flow and existing production force terms. It does not add
convection centers, modify the model, or infer that production plate splitting is caused
by shell failure.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components

from app import cratons, geometry, lithosphere, torque, world as world_mod
import probe_gate3b_component_shell as shell_audit
import probe_gate3b_midplate_current as current_audit
import replay_gate2_capacity as capacity_audit
import gate3b_diagnostic_common as shared

SEEDS = (7, 11, 2024)
DENSITY = 4.0
STARTUP_YEARS = 2_000_000.0
STEP_YEARS = 2_000_000.0
DURATION_YEARS = 30_000_000.0
SHELL_CAP = 512
SLICE_COUNT = 10
OUTER_FRACTION = 0.20
EVENT_AREA_FRACTION = 0.25
EVENT_DURATION_YEARS = 10_000_000.0
RHO_MANTLE = 3300.0


def _column_mass_kg_per_m2(plate):
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    rho_c = lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type)
    return hc * rho_c + hm * RHO_MANTLE


def _principal_axis(plate):
    points = geometry.normalize(plate.surface_nodes().world_xyz)
    areas = plate.node_areas_m2()
    center = geometry.normalize(np.sum(points * areas[:, None], axis=0))
    tangent = points - (points @ center)[:, None] * center
    cov = np.einsum("n,ni,nj->ij", areas, tangent, tangent) / max(float(areas.sum()), 1.0)
    _values, vectors = np.linalg.eigh(cov)
    axis = vectors[:, -1] - center * np.dot(vectors[:, -1], center)
    axis = geometry.normalize(axis)
    pivot = int(np.argmax(np.abs(axis)))
    if axis[pivot] < 0:
        axis = -axis
    return points, areas, center, axis


def _field_and_inputs(world, plate):
    spacing = shell_audit.line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    cache = getattr(world, "boundary_search_cache", None)
    neighbors = plate.get_neighbours(world.plates)
    inputs = torque.gather_boundary_force_inputs(plate, neighbors, spacing, reach, cache)
    _conv, divergent, _transform, contested = torque.classify_boundary_nodes(
        plate, neighbors, inputs, reach, cache
    )
    collision = torque.continental_contact_mask(plate, inputs, contested)
    subducting = torque.subducting_boundary_mask(plate, inputs, reach)
    craton_contact = torque.contact_craton_m(plate, neighbors, inputs)
    jam = torque.contact_jam_factor(plate, world, neighbors, inputs, reach)
    field = torque.gather_per_node_force_field(
        plate, world, inputs, spacing, subducting, divergent, collision, craton_contact, jam
    )
    return field, inputs, subducting


def _bilateral_candidate(world, plate):
    points, areas, center, axis = _principal_axis(plate)
    coord = points @ axis
    lo, hi = np.quantile(coord, [OUTER_FRACTION, 1.0 - OUTER_FRACTION])
    field, inputs, subducting = _field_and_inputs(world, plate)
    slab = np.asarray(field.slab_pull_n)
    active = subducting & (np.linalg.norm(slab, axis=1) > 0.0)
    sides = (coord <= lo, coord >= hi)
    side_rows = []
    own_sigma = _column_mass_kg_per_m2(plate)
    by_id = {int(item.plate_id): item for item in world.plates}
    for side_index, mask in enumerate(sides):
        indices = np.flatnonzero(active & mask)
        if not len(indices):
            return None
        axial = slab[indices] @ axis
        pull_projection = float(axial.sum())
        outward_sign = -1.0 if side_index == 0 else 1.0
        if outward_sign * pull_projection <= 0.0:
            return None
        neighbor_sigma = []
        neighbor_weights = []
        for idx in indices:
            other = by_id.get(int(inputs.neighbor_plate_id[idx]))
            other_idx = int(inputs.neighbor_node_index[idx])
            if other is None or other_idx < 0 or other_idx >= other.node_count():
                continue
            neighbor_sigma.append(float(_column_mass_kg_per_m2(other)[other_idx]))
            neighbor_weights.append(float(field.areas_m2[idx]))
        if not neighbor_sigma:
            return None
        neighbor_mean = float(np.average(neighbor_sigma, weights=neighbor_weights))
        own_mean = float(np.average(own_sigma[indices], weights=field.areas_m2[indices]))
        if neighbor_mean >= own_mean:
            return None
        side_rows.append({
            "active_nodes": int(len(indices)),
            "own_column_mass_kg_m2": own_mean,
            "neighbor_column_mass_kg_m2": neighbor_mean,
            "slab_pull_resultant_n": slab[indices].sum(axis=0).tolist(),
            "axial_outward_projection_n": float(outward_sign * pull_projection),
        })
    return {
        "plate_id": int(plate.plate_id),
        "area_m2": float(areas.sum()),
        "center_xyz": center.tolist(),
        "axis_xyz": axis.tolist(),
        "side_groups": side_rows,
        "points": points,
        "areas": areas,
        "coord": coord,
        "field": field,
    }


def _select(world):
    eligible = []
    for plate in world.plates:
        candidate = _bilateral_candidate(world, plate)
        if candidate is not None:
            eligible.append((candidate["area_m2"], plate, candidate))
    if not eligible:
        return None, []
    eligible.sort(key=lambda row: (-row[0], row[1].plate_id))
    return eligible[0][1], [
        {"plate_id": int(plate.plate_id), "area_m2": float(area), "side_groups": cand["side_groups"]}
        for area, plate, cand in eligible
    ]


def _mass_slices(plate, key_to_bin, edges, initial_points, initial_bins):
    nodes = plate.surface_nodes()
    keys = nodes.node_ids[:, 0].view(np.int64)
    points = geometry.normalize(nodes.world_xyz)
    initial_tree = cKDTree(initial_points)
    bins = np.empty(len(keys), dtype=np.int64)
    for i, key in enumerate(keys):
        known = key_to_bin.get(int(key))
        bins[i] = known if known is not None else int(initial_bins[initial_tree.query(points[i], k=1)[1]])
    areas = plate.node_areas_m2()
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    rho_c = lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type)
    mass = areas * (hc * rho_c + hm * RHO_MANTLE)
    crust_volume = areas * hc
    mantle_volume = areas * hm
    rows = []
    for s in range(SLICE_COUNT):
        mask = bins == s
        rows.append({
            "slice": s,
            "node_count": int(mask.sum()),
            "area_km2": float(areas[mask].sum() / 1e6),
            "crust_volume_km3": float(crust_volume[mask].sum() / 1e9),
            "mantle_lithosphere_volume_km3": float(mantle_volume[mask].sum() / 1e9),
            "mass_kg": float(mass[mask].sum()),
        })
    return bins, rows


def _shell_snapshot(world, plate, field):
    all_indices = np.arange(plate.node_count(), dtype=np.int64)
    result = shell_audit._plate_shell(
        plate, all_indices, field.points_xyz, field.total_n, SHELL_CAP, include_demand=True
    )
    demand = np.asarray(result.pop("component_demand_n_per_m"), dtype=float)
    areas = plate.node_areas_m2()
    graph = current_audit._component_graph(plate, all_indices)
    capacities = capacity_audit._physical_capacities()["continental"]
    cases = {}
    for name, cap in zip(("p10", "median", "p90"), capacities):
        over = demand >= cap
        count, labels = connected_components(graph[over][:, over], directed=False) if over.any() else (0, np.array([], dtype=int))
        if over.any():
            component_area = np.bincount(labels, weights=areas[over], minlength=count)
            max_fraction = float(component_area.max() / areas.sum())
        else:
            max_fraction = 0.0
        cases[name] = {
            "capacity_n_per_m": float(cap),
            "max_connected_overcapacity_footprint_fraction": max_fraction,
            "qualifies_checkpoint": bool(max_fraction >= EVENT_AREA_FRACTION),
        }
    result["capacity_cases"] = cases
    result["demand_quantiles_n_per_m"] = np.quantile(demand[np.isfinite(demand)], [0.1, 0.5, 0.9]).tolist()
    return result, demand


def _checkpoint(world, plate, key_to_bin, edges, initial_points, initial_bins):
    field, _inputs, _subducting = _field_and_inputs(world, plate)
    bins, slices = _mass_slices(plate, key_to_bin, edges, initial_points, initial_bins)
    terms = {
        "slab_pull": field.slab_pull_n,
        "slab_drag": field.slab_drag_n,
        "basal_drag": field.basal_drag_n,
        "combined": field.total_n,
    }
    for row in slices:
        mask = bins == row["slice"]
        mass = row["mass_kg"]
        row["loads"] = {}
        for name, vectors in terms.items():
            summed = np.asarray(vectors)[mask].sum(axis=0)
            row["loads"][name] = {
                "resultant_n": summed.tolist(),
                "resultant_magnitude_n": float(np.linalg.norm(summed)),
                "force_per_mass_n_per_kg": float(np.linalg.norm(summed) / mass) if mass > 0 else 0.0,
            }
    shell, demand = _shell_snapshot(world, plate, field)
    areas = plate.node_areas_m2()
    capacities = shell["capacity_cases"]
    for row in slices:
        mask = bins == row["slice"]
        local = demand[mask]
        local_area = areas[mask]
        finite = np.isfinite(local)
        if not finite.any():
            row["shell_demand"] = {"valid_nodes": 0}
            continue
        row["shell_demand"] = {
            "valid_nodes": int(finite.sum()),
            "median_n_per_m": float(np.median(local[finite])),
            "p90_n_per_m": float(np.quantile(local[finite], 0.9)),
            "area_over_p10_capacity_fraction": float(
                local_area[finite & (local >= capacities["p10"]["capacity_n_per_m"])].sum()
                / max(float(local_area[finite].sum()), 1.0)
            ),
            "area_over_p90_capacity_fraction": float(
                local_area[finite & (local >= capacities["p90"]["capacity_n_per_m"])].sum()
                / max(float(local_area[finite].sum()), 1.0)
            ),
        }
    return {
        "elapsed_myr": float(world.elapsed_years / 1e6),
        "plate_id": int(plate.plate_id),
        "plate_area_km2": float(plate.node_areas_m2().sum() / 1e6),
        "shell": shell,
        "slices": slices,
    }


def run(out_path: Path):
    started = time.perf_counter()
    chosen = None
    selection_log = []
    for seed in SEEDS:
        try:
            world = world_mod.generate_world(seed, node_density=DENSITY, cratons_first=True, synthetic_relief=False)
        except Exception as exc:
            selection_log.append({
                "seed": seed,
                "setup_error": f"{type(exc).__name__}: {exc}",
            })
            continue
        world.simulate_climate_biomes = False
        cratons.ensure_initialized(world)
        world_mod.step_world(world, STARTUP_YEARS)
        plate, candidates = _select(world)
        selection_log.append({
            "seed": seed,
            "elapsed_after_startup_myr": float(world.elapsed_years / 1e6),
            "eligible_candidates": candidates,
        })
        if plate is not None:
            chosen = (seed, world, plate)
            break
    if chosen is None:
        report = {
            "title": "RFC 0010 Gate 3b bilateral slab-pull screen",
            "eligible": False,
            "reason": "no plate met the preregistered bilateral down-going-boundary and lighter-neighbor rule",
            "seed_order": list(SEEDS), "node_density": DENSITY,
            "selection": selection_log,
            "elapsed_wall_seconds": time.perf_counter() - started,
        }
        shared.write_json_report(out_path, report)
        return report

    seed, world, plate = chosen
    candidate = _bilateral_candidate(world, plate)
    initial_nodes = plate.surface_nodes()
    initial_points = geometry.normalize(initial_nodes.world_xyz)
    keys = initial_nodes.node_ids[:, 0].view(np.int64)
    coord = candidate["coord"]
    edges = np.linspace(float(coord.min()), float(coord.max()), SLICE_COUNT + 1)
    initial_bins = np.clip(np.digitize(coord, edges[1:-1]), 0, SLICE_COUNT - 1)
    key_to_bin = {int(key): int(b) for key, b in zip(keys, initial_bins)}
    checkpoints = []
    persistent = {"p10": None, "median": None, "p90": None}
    streak = {key: 0.0 for key in persistent}
    onset = {key: None for key in persistent}
    first_split_myr = None
    split_message = None
    error = None
    initial_snapshot = _checkpoint(world, plate, key_to_bin, edges, initial_points, initial_bins)
    initial_snapshot["accounting_pass"] = bool(initial_snapshot["shell"].get("accounting_pass", False))
    checkpoints.append(initial_snapshot)
    if not initial_snapshot["accounting_pass"]:
        error = "initial full-load shell accounting failed"
    for _ in range(0 if error else int(round(DURATION_YEARS / STEP_YEARS))):
        old_event_count = len(world.events)
        world_mod.step_world(world, STEP_YEARS)
        new_events = world.events[old_event_count:]
        split_prefix = f"Plate {plate.plate_id} split into"
        found = next((message for _time, message in new_events if message.startswith(split_prefix)), None)
        if found is not None:
            first_split_myr = float(world.elapsed_years / 1e6)
            split_message = found
            break
        current = next((p for p in world.plates if p.plate_id == plate.plate_id), None)
        if current is None:
            error = "tracked plate id disappeared before a logged production split"
            break
        try:
            snap = _checkpoint(world, current, key_to_bin, edges, initial_points, initial_bins)
        except Exception as exc:  # diagnostic failures remain explicit and ineligible
            error = f"{type(exc).__name__}: {exc}"
            break
        accounting = snap["shell"]
        snap["accounting_pass"] = bool(accounting.get("accounting_pass", False))
        if not snap["accounting_pass"]:
            snap["history_stopped_for_accounting"] = True
            checkpoints.append(snap)
            break
        for name, case in accounting["capacity_cases"].items():
            if case["qualifies_checkpoint"]:
                if onset[name] is None:
                    onset[name] = snap["elapsed_myr"]
                streak[name] += STEP_YEARS
                if streak[name] >= EVENT_DURATION_YEARS and persistent[name] is None:
                    persistent[name] = onset[name]
            else:
                onset[name] = None
                streak[name] = 0.0
            case["current_consecutive_duration_myr"] = streak[name] / 1e6
        checkpoints.append(snap)
    report = {
        "title": "RFC 0010 Gate 3b bilateral slab-pull screen",
        "hypothesis": "Opposing slab-pull loads at two heavy-plate ends create sustained central extension under native mantle flow without added convection centers.",
        "eligible": error is None,
        "seed": seed, "node_density": DENSITY,
        "startup_myr": STARTUP_YEARS / 1e6,
        "step_myr": STEP_YEARS / 1e6,
        "planned_duration_myr": DURATION_YEARS / 1e6,
        "selected_plate_id": int(plate.plate_id),
        "selected_plate_area_km2": float(candidate["area_m2"] / 1e6),
        "axis_xyz": candidate["axis_xyz"],
        "bilateral_precondition": candidate["side_groups"],
        "selection": selection_log,
        "slice_count": SLICE_COUNT,
        "slice_edges_projection": edges.tolist(),
        "persistent_event_onset_myr_by_capacity": persistent,
        "first_production_breakup_myr": first_split_myr,
        "breakup_right_censored_at_myr": None if first_split_myr is not None else float(world.elapsed_years / 1e6),
        "production_split_message": split_message,
        "diagnostic_error": error,
        "checkpoints": checkpoints,
        "elapsed_wall_seconds": time.perf_counter() - started,
        "interpretation_limit": "The production split gate clusters mantle-flow disagreement and is not triggered by the diagnostic shell demand; split time is secondary and does not establish slab-pull causation.",
    }
    shared.write_json_report(out_path, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print(json.dumps({
        "out": str(args.out), "eligible": report.get("eligible"),
        "seed": report.get("seed"), "plate_id": report.get("selected_plate_id"),
        "persistent_event_onset_myr_by_capacity": report.get("persistent_event_onset_myr_by_capacity"),
        "first_production_breakup_myr": report.get("first_production_breakup_myr"),
        "censored_at_myr": report.get("breakup_right_censored_at_myr"),
        "diagnostic_error": report.get("diagnostic_error"),
        "elapsed_wall_seconds": report.get("elapsed_wall_seconds"),
    }, sort_keys=True))


if __name__ == "__main__":
    main()

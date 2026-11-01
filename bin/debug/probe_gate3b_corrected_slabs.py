#!/usr/bin/env python3
"""Run the corrected, controlled evolving-slab screen for RFC 0010 Gate 3b (#699).

Replaces the historical hinge-free screen in `probe_gate3b_evolving_slabs.py`. Each load
history advances the ribbons, solves the shell and accumulates the p10, median and p90
damage cases together, since damage does not feed back on the shell. The arms and cutoffs
are the ones frozen in the issue before this ran; see `analysis/issue699/README.md`.

Usage (one process; RAM-heavy shell mechanics):
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
      backend/.venv/bin/python bin/debug/probe_gate3b_corrected_slabs.py \
      --save <registered.mbworld> \
      --out analysis/issue699/corrected-slabs-seed447276134-14.json
"""
from __future__ import annotations

import argparse
import copy
import gc
import time
from pathlib import Path

import numpy as np
from app import geometry, torque
from app import slab_body as slab_model
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import gate3b_diagnostic_common as common
import probe_gate3b_crack_coupled as damage_model
import probe_gate3b_evolving_slabs as evolving
import probe_gate3b_m1_saved_plates as saved_screen

DURATION_MYR = 20
DETACH_MYR = 10.0
EVENT_AREA_FRACTION = 0.25
DAUGHTER_AREA_FRACTION = 0.10
CAPACITY_CASES = ("p10", "median", "p90")
SHELL_CUTOFF = 1e-3
HINGE_WORK_CUTOFF = 1e-3
MATERIAL_CUTOFF = 1e-10
PATH_CUTOFF = 1e-3


def _case(name, *, mass="column", sides=(True, True), current=True, drag=1.0, step=2.0, detach=None):
    return {"name": name, "mass": mass, "sides": sides, "current": current,
            "axial_drag_factor": drag, "step_myr": step, "detach_myr": detach}


def registered_cases():
    """Frozen load histories: drag bracket, flank, current, detachment and timestep arms."""
    cases = [_case("R_registered_mass_both_flanks", mass="registered"),
             _case("R_registered_mass_no_slab_force", mass="registered", sides=(False, False))]
    for step in (2.0, 1.0):
        tag = f"{step:g}myr"
        cases.append(_case(f"A_no_slab_force_{tag}", sides=(False, False), step=step))
        for drag, label in ((1.0, "full_drag"), (0.5, "half_drag"), (0.0, "gross")):
            cases.append(_case(f"B_stronger_flank_{label}_{tag}", sides=None, drag=drag, step=step))
            cases.append(_case(f"C_both_flanks_{label}_{tag}", drag=drag, step=step))
            cases.append(_case(f"N_both_flanks_no_current_{label}_{tag}", current=False, drag=drag, step=step))
            cases.append(_case(f"D_detach_10myr_{label}_{tag}", drag=drag, step=step, detach=DETACH_MYR))
    return cases


def baseline_components(plate):
    """Label the intact plate graph; damage is judged only inside its largest component."""
    adjacency = plate.adjacency()
    count = plate.node_count()
    graph = csr_matrix(
        (np.ones(len(adjacency.neighbours)), adjacency.neighbours, adjacency.offsets),
        shape=(count, count),
    )
    _count, labels = connected_components(graph, directed=False)
    areas = plate.node_areas_m2()
    return labels == int(np.argmax(np.bincount(labels, weights=areas)))


def baseline_adjusted_split(labels, main, areas, initial_area):
    """Count new pieces of the main component so an inherited island cannot pass or block."""
    pieces = np.unique(labels[main])
    piece_area = np.array([areas[main & (labels == piece)].sum() for piece in pieces])
    order = np.sort(piece_area)[::-1] / initial_area
    return {
        "new_components": int(len(pieces) - 1),
        "two_largest_area_fractions": order[:2].tolist(),
        "viable_split": bool(len(pieces) == 2 and order[1] >= DAUGHTER_AREA_FRACTION),
    }


def _connected_fraction(plate, over, areas):
    """Largest edge-connected over-capacity area as a fraction of plate area."""
    if not np.any(over):
        return 0.0
    adjacency = plate.adjacency()
    rows = np.repeat(np.arange(plate.node_count()), np.diff(adjacency.offsets))
    keep = over[rows] & over[adjacency.neighbours]
    graph = csr_matrix(
        (np.ones(int(keep.sum())), (rows[keep], adjacency.neighbours[keep])),
        shape=(plate.node_count(), plate.node_count()),
    )
    _count, labels = connected_components(graph, directed=False)
    return float(np.bincount(labels[over], weights=areas[over]).max() / areas.sum())


def event_window_myr(rows):
    """Longest run of outputs with the connected footprint at or above the frozen cutoff."""
    best = run = 0.0
    previous = 0.0
    for row in rows:
        run = run + row["elapsed_myr"] - previous if row["footprint"] >= EVENT_AREA_FRACTION else 0.0
        previous = row["elapsed_myr"]
        best = max(best, run)
    return best


def _affine_share(points, displacement):
    """Fraction of nonrigid displacement variance the plate-wide strain fit keeps."""
    mean = displacement - displacement.mean(axis=0)
    design = np.column_stack((points - points.mean(axis=0), np.ones(len(points))))
    fitted = design @ np.linalg.lstsq(design, mean, rcond=None)[0]
    return float(np.sum(fitted**2) / max(np.sum(mean**2), 1.0))


def _flags(case, stronger):
    if case["sides"] is None:
        return tuple(index == stronger for index in range(2))
    return tuple(case["sides"])


def _advance(body, speed, step_myr, drag_factor):
    return slab_model.advance_slab_body(
        body, step_myr * 1e6,
        mantle_viscosity_pa_s=evolving.MANTLE_VISCOSITY_PA_S,
        shear_zone_thickness_m=evolving.SHEAR_ZONE_M,
        convergence_speed_m_per_s=speed,
        axial_drag_factor=drag_factor,
    )


def _load_rows(bodies, budgets, enabled):
    return [{
        "enabled": bool(active and body.attached),
        "in_plane_hinge_force_n": float(np.linalg.norm(budget["attachment_force_local_n"])),
        "tension_n": budget["hinge_tension_n"],
        "bending_shear_n": budget["hinge_bending_shear_n"],
        "buoyancy_n": budget["axial_buoyancy_force_n"],
        "axial_drag_n": budget["axial_mantle_drag_force_n"],
        "plate_driver_n": budget["plate_driver_force_n"],
        "beam_work_residual": budget["beam_work_residual_fraction"],
        "axial_work_residual": budget["axial_work_residual_fraction"],
        "geometry": body.metrics(),
    } for body, budget, active in zip(bodies, budgets, enabled)]


def _shell_residuals(response):
    force_scale = max(response["source_force_norm_sum_n"], 1.0)
    torque_scale = max(response["source_torque_norm_sum_nm"], 1.0)
    applied = max(float(np.linalg.norm(response["applied_resultant_force_n"])), 1.0)
    return {
        "force": float(np.linalg.norm(response["force_transfer_error_n"]) / force_scale),
        "torque": float(np.linalg.norm(response["torque_transfer_error_nm"]) / torque_scale),
        "equilibrium": float(response["equilibrium_residual_n"] / applied),
        "work": float(response["work_residual_fraction"]),
    }


def _damage_rows(plate, main, damage, stress, direction, yields, capacities, thickness, step_s, elapsed):
    areas = np.asarray(plate.node_areas_m2(), dtype=float)
    rows = {}
    for name in CAPACITY_CASES:
        damage[name] = damage_model._advance_damage(damage[name], stress, yields[name], step_s)
        cracked = damage[name] >= damage_model.EPSILON_FAIL
        broken, _components, labels, _sizes = saved_screen._components(
            plate, cracked, direction, stress, yields[name]
        )
        over = (stress * thickness) >= capacities[name]
        rows[name] = {
            "elapsed_myr": elapsed,
            "footprint": _connected_fraction(plate, over, areas),
            "failed_cells": int(cracked.sum()),
            "failed_area_fraction": float(areas[cracked].sum() / areas.sum()),
            "failed_links": int(len(broken)),
            **baseline_adjusted_split(labels, main, areas, float(areas.sum())),
        }
    return rows


def _run_history(initial, template, bodies_template, records, case, capacities, thickness, main):
    plate = copy.deepcopy(template)
    bodies = copy.deepcopy(bodies_template)
    world = initial["world"]
    index = next(i for i, p in enumerate(world.plates) if int(p.plate_id) == evolving.PLATE_ID)
    world.plates[index] = plate
    enabled = _flags(case, initial["stronger_side"])
    source = np.unique(np.concatenate([np.asarray(row["source_indices"]) for row in records]))
    capacities = {name: np.where(np.isin(np.arange(len(field)), source), 1e30, field)
                  for name, field in capacities.items()}
    yields = {name: field / np.maximum(thickness, 1.0) for name, field in capacities.items()}
    damage = {name: np.zeros(plate.node_count()) for name in CAPACITY_CASES}
    step_myr = case["step_myr"]
    history = {name: [] for name in CAPACITY_CASES}
    checkpoints = []
    surface_initial = evolving._surface_inventory(plate)
    volume_initial = [body.consumed_volume_m3 for body in bodies]
    length_initial = [body.metrics()["length_m"] for body in bodies]
    fed_m = [0.0, 0.0]
    for step in range(1, int(round(DURATION_MYR / step_myr)) + 1):
        elapsed = step * step_myr
        if case["detach_myr"] is not None and elapsed > case["detach_myr"]:
            for body in bodies:
                body.attached = False
        speeds = evolving._closing_speeds(world, plate, initial["active_masks"])
        budgets = [_advance(body, speeds[i] if body.attached else 0.0, step_myr, case["axial_drag_factor"])
                   for i, body in enumerate(bodies)]
        for side, budget in enumerate(budgets):
            fed_m[side] += budget["advected_distance_m"]
        loaded = [body for body, active in zip(bodies, enabled) if active and body.attached]
        current = (evolving._midline_force(plate, initial["midline_source_centers"])
                   if case["current"] else np.zeros((plate.node_count(), 3)))
        if not loaded and not case["current"]:
            response = None
        else:
            response = slab_model.apply_slab_load_to_plate_geometry(
                plate, loaded, step_myr * 1e6,
                effective_extensional_viscosity_pa_s=evolving.EXTENSIONAL_VISCOSITY_PA_S,
                load_spread_degrees=evolving.LOAD_SPREAD_DEG,
                extensional_capacity_resultant_n_per_m=capacities["median"],
                shell_support_thickness_m=thickness,
                additional_force_local_n=geometry.to_local(plate.frame, current),
            )
        if response is None:
            demand = np.zeros(plate.node_count())
            direction = np.zeros((plate.node_count(), 3))
            direction[:, 0] = 1.0
            residuals = {"force": 0.0, "torque": 0.0, "equilibrium": 0.0, "work": 0.0}
            strain, coverage, affine = 0.0, 1.0, 0.0
        else:
            demand = np.asarray(response["max_tensile_resultant_n_per_m"], dtype=float)
            direction = geometry.normalize(np.asarray(response["principal_extension_xyz"], dtype=float))
            residuals = _shell_residuals(response)
            strain = float(response["max_log_strain"])
            coverage = float(response["shell_area_coverage_fraction"])
            affine = _affine_share(np.asarray(plate.surface_nodes().local_xyz),
                                   np.asarray(response["shell_displacement_m"]))
        stress = demand / np.maximum(thickness, 1.0)
        rows = _damage_rows(plate, main, damage, stress, direction, yields, capacities, thickness,
                            step_myr * 1e6 * torque.SECONDS_PER_YEAR, elapsed)
        for name, row in rows.items():
            history[name].append(row)
        surface = evolving._surface_inventory(plate)
        surface_residual = max(abs(surface[key] - value) / max(abs(value), 1.0)
                               for key, value in surface_initial.items())
        body_residual = max(
            abs(body.consumed_volume_m3 - start - body.kinematic_inflow_volume_m3) / start
            for body, start in zip(bodies, volume_initial)
        )
        hinge_work = max(max(b["beam_work_residual_fraction"], b["axial_work_residual_fraction"]) for b in budgets)
        path_residual = max(
            abs(body.metrics()["length_m"] - start - fed) / (start + fed)
            for body, start, fed in zip(bodies, length_initial, fed_m)
        )
        valid = bool(max(residuals.values()) <= SHELL_CUTOFF and hinge_work <= HINGE_WORK_CUTOFF
                     and path_residual <= PATH_CUTOFF
                     and max(surface_residual, body_residual) <= MATERIAL_CUTOFF)
        checkpoints.append({
            "elapsed_myr": elapsed, "fronts": _load_rows(bodies, budgets, enabled),
            "convergence_m_per_s": speeds, "shell_residuals": residuals,
            "shell_area_coverage": coverage, "max_log_strain": strain,
            "affine_displacement_share": affine,
            "demand_n_per_m_p50_p90_p99": np.percentile(demand, [50, 90, 99]).tolist(),
            "surface_material_residual": surface_residual, "slab_volume_residual": body_residual,
            "hinge_work_residual": hinge_work, "path_length_residual": path_residual, "accounting_valid": valid,
        })
        if not valid:
            break
    outcomes = {}
    for name in CAPACITY_CASES:
        rows = history[name]
        outcomes[name] = {
            "event_window_myr": event_window_myr(rows),
            "peak_footprint": max(row["footprint"] for row in rows),
            "first_failed_myr": next((r["elapsed_myr"] for r in rows if r["failed_cells"]), None),
            "first_viable_split_myr": next((r["elapsed_myr"] for r in rows if r["viable_split"]), None),
            "final": rows[-1],
            "history": rows[:: max(1, int(round(2.0 / step_myr)))],
        }
    return {**case, "sides": list(enabled), "outcomes": outcomes, "checkpoints": checkpoints[-1:],
            "first_checkpoint": checkpoints[0], "completed_myr": checkpoints[-1]["elapsed_myr"],
            "accounting_valid": all(row["accounting_valid"] for row in checkpoints)}


def run(save_path: Path, out_path: Path) -> dict:
    started = time.perf_counter()
    initial = evolving._load_target(save_path)
    base = initial["plate"]
    thickness = (base.collect("crustal_thickness_m")
                 + base.collect("mantle_lithosphere_thickness_m")).astype(float)
    capacities, references = damage_model._capacity_fields(base.collect("continental_material_m") > 0.0)
    main = baseline_components(base)
    templates = {}
    for mass, column in (("registered", False), ("column", True)):
        plate = copy.deepcopy(base)
        bodies, records = evolving._create_bodies(plate, initial, column_thickness=column)
        templates[mass] = (plate, bodies, records)
    results = []
    for case in registered_cases():
        plate, bodies, records = templates[case["mass"]]
        results.append(_run_history(initial, plate, bodies, records, case, capacities, thickness, main))
        print(f"{case['name']}: {time.perf_counter() - started:.0f} s", flush=True)
        gc.collect()
    report = {
        "schema_version": 1,
        "scope": "Offline mechanics screen with an affine plate-geometry update; not Gate 3b closure.",
        "seed": int(initial["world"].seed),
        "plate_id": evolving.PLATE_ID,
        "stronger_flank": "AB"[initial["stronger_side"]],
        "inherited_islands_area_fraction": float(
            base.node_areas_m2()[~main].sum() / base.node_areas_m2().sum()
        ),
        "capacity_references_n_per_m": references,
        "fronts": {mass: [{"side": record["side"], "thickness_m": body.thickness_m, "width_m": body.width_m,
                           "source_volume_m3": record["source_volume_m3"],
                           "antecedent_volume_m3": body.antecedent_volume_m3}
                          for body, record in zip(item[1], item[2])]
                   for mass, item in templates.items()},
        "results": results,
        "runtime_seconds": time.perf_counter() - started,
    }
    common.write_json_report(out_path, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", type=Path, required=True, help="registered .mbworld input")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.save, args.out)
    print(f"{len(report['results'])} histories, {report['runtime_seconds']:.0f} s")


if __name__ == "__main__":
    main()

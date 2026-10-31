#!/usr/bin/env python3
"""Run component-bound Gate 3b shell/load-accounting preflight.

Usage:
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
        backend/.venv/bin/python bin/debug/probe_gate3b_component_shell.py \
        --out analysis/issue699/rift-potential-component-accounting-seed11-density4.json

Regenerates the preregistered seed-11 world, takes the 2 Myr startup step, selects the
largest edge-connected continental component by represented area, then checks the shell
mesh and source/work accounting on that component. It does not score capacity or persistence.
"""
from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from app import cratons, elevation_lines, geometry, lithosphere, mantle, thin_shell, torque, world as world_mod
from app.elevation_lines import line_spacing_rad
import validate_thin_shell_plate_meshes as mesh_validation

SEED = 11
DENSITY = 4.0
STARTUP_STEP_YEARS = 2_000_000.0
MESH_CAPS = (512, 1024)
EFFECTIVE_MODULUS_PA = 70e9
POISSON_RATIO = 0.25
PLACEHOLDER_CAPACITY_N_PER_M = 1e10


def _largest_continental_component(world):
    components = []
    for plate in world.plates:
        continental = elevation_lines.effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        indices = np.flatnonzero(continental)
        if not len(indices):
            continue
        adjacency = plate.adjacency()
        rows = np.repeat(np.arange(plate.node_count(), dtype=np.int64), np.diff(adjacency.offsets))
        cols = adjacency.neighbours.astype(np.int64, copy=False)
        local_index = np.full(plate.node_count(), -1, dtype=np.int64)
        local_index[indices] = np.arange(len(indices), dtype=np.int64)
        keep = continental[rows] & continental[cols]
        graph = csr_matrix(
            (np.ones(int(keep.sum())), (local_index[rows[keep]], local_index[cols[keep]])),
            shape=(len(indices), len(indices)),
        )
        count, labels = connected_components(graph, directed=False)
        areas = plate.node_areas_m2()[indices]
        for label in range(count):
            selected = indices[labels == label]
            area = float(areas[labels == label].sum())
            components.append((area, plate, selected))
    if not components:
        raise ValueError("generated world has no edge-connected continental component")
    return max(components, key=lambda row: row[0])


def _production_force_field(world, plate):
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    cache = getattr(world, "boundary_search_cache", None)
    neighbours = plate.get_neighbours(world.plates)
    inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
    _convergent, divergent, _transform, contested = torque.classify_boundary_nodes(
        plate, neighbours, inputs, reach, cache
    )
    collision = torque.continental_contact_mask(plate, inputs, contested)
    subducting = torque.subducting_boundary_mask(plate, inputs, reach)
    craton_contact = torque.contact_craton_m(plate, neighbours, inputs)
    jam = torque.contact_jam_factor(plate, world, neighbours, inputs, reach)
    return torque.gather_per_node_force_field(
        plate, world, inputs, spacing, subducting, divergent, collision, craton_contact, jam
    )


def _representatives(points, areas, count):
    selected = np.empty(min(count, len(points)), dtype=np.int64)
    selected[0] = int(np.argmax(areas))
    nearest = np.full(len(points), np.inf)
    for i in range(1, len(selected)):
        delta = points - points[selected[i - 1]]
        nearest = np.minimum(nearest, np.einsum("ij,ij->i", delta, delta))
        nearest[selected[:i]] = -1.0
        selected[i] = int(np.argmax(nearest))
    return selected


def _plate_shell(plate, component_indices, force_points, source_force, cap, *, include_demand=False):
    """Solve on the complete production plate; retain the continent as an evaluation mask."""
    source_points = geometry.normalize(np.asarray(force_points, dtype=float))
    source_areas = plate.node_areas_m2()
    representatives = _representatives(source_points, source_areas, cap)
    rep_points = source_points[representatives]
    assignment = cKDTree(rep_points).query(source_points, k=1)[1]
    rep_areas = np.bincount(assignment, weights=source_areas, minlength=len(rep_points))
    forces = np.asarray(source_force, dtype=float)
    rep_force = np.zeros((len(rep_points), 3))
    np.add.at(rep_force, assignment, forces)
    thickness = np.maximum(
        plate.collect("crustal_thickness_m")
        + plate.collect("mantle_lithosphere_thickness_m"), 1.0
    )
    thick_sum = np.bincount(assignment, weights=source_areas * thickness, minlength=len(rep_points))
    rep_thickness = thick_sum / np.maximum(rep_areas, np.finfo(float).tiny)

    triangles, rejected, candidates = mesh_validation._candidate_triangles(plate, rep_points)
    labels = thin_shell._mesh_component_labels(triangles)
    sizes = np.bincount(labels)
    largest = int(np.argmax(sizes))
    triangles = triangles[labels == largest]
    used = np.unique(triangles)
    remap = np.full(len(rep_points), -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    local_triangles = remap[triangles]
    if len(thin_shell._mesh_component_sizes(len(used), local_triangles)) != 1:
        raise ValueError("component-bound shell is disconnected after largest-component filter")

    # Keep raw mesh coverage visible, but transfer force and column area from medoids omitted
    # by clipping to their nearest retained shell node instead of dropping those source bins.
    coverage = float(rep_areas[used].sum() / source_areas.sum())
    orphan = np.setdiff1d(np.arange(len(rep_points)), used, assume_unique=True)
    orphan_targets = cKDTree(rep_points[used]).query(rep_points[orphan], k=1)[1] if len(orphan) else np.empty(0, dtype=int)
    shell_areas = rep_areas[used].copy()
    shell_force = rep_force[used].copy()
    shell_thickness_sum = thick_sum[used].copy()
    if len(orphan):
        np.add.at(shell_areas, orphan_targets, rep_areas[orphan])
        np.add.at(shell_force, orphan_targets, rep_force[orphan])
        np.add.at(shell_thickness_sum, orphan_targets, thick_sum[orphan])
    shell_thickness = shell_thickness_sum / np.maximum(shell_areas, np.finfo(float).tiny)
    source_f = forces.sum(axis=0)
    transfer_f = shell_force.sum(axis=0)
    source_t = np.cross(source_points * lithosphere.PLANET_RADIUS_M, forces)
    transfer_t = np.cross(rep_points[used] * lithosphere.PLANET_RADIUS_M, shell_force)
    f_scale = max(float(np.linalg.norm(forces, axis=1).sum()), 1.0)
    t_scale = max(float(np.linalg.norm(source_t, axis=1).sum()), 1.0)

    solver_stats = {}
    original_lsmr = thin_shell.lsmr

    def record_lsmr(*args, **kwargs):
        result = original_lsmr(*args, **kwargs)
        solver_stats.update({
            "stop_code": int(result[1]), "iterations": int(result[2]),
            "scaled_residual_norm": float(result[3]), "scaled_normal_residual_norm": float(result[4]),
            "scaled_operator_norm": float(result[5]), "condition_estimate": float(result[6]),
            "scaled_solution_norm": float(result[7]),
        })
        return result

    thin_shell.lsmr = record_lsmr
    try:
        shell = thin_shell.solve_thin_shell(
        points_xyz=rep_points[used], areas_m2=shell_areas, applied_force_n=shell_force,
        thickness_m=shell_thickness, young_modulus_pa=EFFECTIVE_MODULUS_PA,
        extensional_capacity_resultant_n_per_m=np.full(len(used), PLACEHOLDER_CAPACITY_N_PER_M),
        poisson_ratio=POISSON_RATIO, mesh_triangles=local_triangles,
        )
    finally:
        thin_shell.lsmr = original_lsmr
    equilibrium = shell.equilibrium_residual_n / max(float(np.linalg.norm(shell.equilibrated_force_n)), 1.0)

    local_shell_points = geometry.to_local(plate.frame, rep_points[used])
    local_shell_displacement = geometry.to_local(plate.frame, shell.displacement_m)
    strain_increment = thin_shell.fit_homogeneous_strain_tensor(local_shell_points, local_shell_displacement)
    geometry_response = copy.deepcopy(plate)
    old_points = plate.surface_nodes().local_xyz
    old_areas = plate.node_areas_m2()
    old_hc = plate.collect("crustal_thickness_m")
    old_hm = plate.collect("mantle_lithosphere_thickness_m")
    old_continental = plate.collect("continental_material_m")
    geometry_response.apply_surface_strain(strain_increment)
    new_points = geometry_response.surface_nodes().local_xyz
    new_areas = geometry_response.node_areas_m2()
    angular_change = np.arccos(np.clip(np.einsum("ij,ij->i", old_points, new_points), -1.0, 1.0))

    # Independently integrate element strain energy from the solved displacement field.
    first, second = thin_shell._tangent_basis(rep_points[used])
    displacement_dofs = np.column_stack((np.sum(shell.displacement_m * first, axis=1),
                                         np.sum(shell.displacement_m * second, axis=1))).ravel()
    elastic_energy = 0.0
    modulus = np.full(len(used), EFFECTIVE_MODULUS_PA)
    for tri in local_triangles:
        _stiffness, strain_matrix, constitutive, element_area, _axis_x, _axis_y = thin_shell._triangle_element(
            rep_points[used], tri, first, second, rep_thickness[used], modulus,
            POISSON_RATIO, lithosphere.PLANET_RADIUS_M,
        )
        element_dofs = np.array([dof for node in tri for dof in (2 * node, 2 * node + 1)])
        strain = strain_matrix @ displacement_dofs[element_dofs]
        elastic_energy += 0.5 * float(strain @ constitutive @ strain) * element_area
    # The solver removes rigid-rotation load before its non-rigid membrane solve. Keep
    # rigid plate work in the separate basal-drag account above.
    external_work = 0.5 * float(np.sum(shell.equilibrated_force_n * shell.displacement_m))
    energy_residual = abs(external_work - elastic_energy) / max(abs(external_work), abs(elastic_energy), 1.0)

    component_nodes = geometry.normalize(np.asarray(force_points, dtype=float)[component_indices])
    component_mesh_map = cKDTree(rep_points[used]).query(component_nodes, k=1)[1]
    component_demand = shell.max_incremental_extensional_resultant_n_per_m[component_mesh_map]
    finite = component_demand[np.isfinite(component_demand)]
    checks = {
        "area_coverage_pass": bool(coverage >= 0.95),
        "force_transfer_pass": bool(np.linalg.norm(source_f - transfer_f) / f_scale <= 0.01),
        "torque_transfer_pass": bool(np.linalg.norm(source_t.sum(axis=0) - transfer_t.sum(axis=0)) / t_scale <= 0.01),
        "equilibrium_pass": bool(equilibrium <= 0.01),
        "shell_energy_pass": bool(energy_residual <= 0.05),
    }
    result = {
        "target_representatives": cap,
        "representative_count": int(len(rep_points)),
        "shell_node_count": int(len(used)),
        "mesh_triangles": int(len(triangles)),
        "triangle_candidates": int(candidates),
        "rejected_triangles": int(rejected),
        "triangle_component_sizes_before_filter": sizes.tolist(),
        "discarded_triangles": int(sum(sizes) - sizes[largest]),
        "orphan_representatives_transferred": int(len(orphan)),
        "production_plate_area_coverage": coverage,
        "continental_component_evaluation_nodes": int(len(component_indices)),
        "source_force_n": source_f.tolist(),
        "transferred_force_n": transfer_f.tolist(),
        "force_transfer_residual_fraction": float(np.linalg.norm(source_f - transfer_f) / f_scale),
        "source_torque_nm": source_t.sum(axis=0).tolist(),
        "transferred_torque_nm": transfer_t.sum(axis=0).tolist(),
        "torque_transfer_residual_fraction": float(np.linalg.norm(source_t.sum(axis=0) - transfer_t.sum(axis=0)) / t_scale),
        "solver_stop_code": int(shell.solver_stop_code),
        "solver_stats": solver_stats,
        "equilibrated_force_resultant_n": shell.equilibrated_force_n.sum(axis=0).tolist(),
        "equilibrium_residual_fraction": float(equilibrium),
        "incremental_external_work_j": external_work,
        "stored_elastic_energy_j": float(elastic_energy),
        "shell_energy_residual_fraction": float(energy_residual),
        "affine_geometry_response": {
            "method": "fit homogeneous strain to static thin-shell displacement, then apply once to a copied sparse-quad plate",
            "strain_increment": strain_increment.tolist(),
            "max_surface_shift_km": float(np.max(angular_change, initial=0.0) * lithosphere.PLANET_RADIUS_M / 1000.0),
            "rms_surface_shift_km": float(np.sqrt(np.mean(angular_change**2)) * lithosphere.PLANET_RADIUS_M / 1000.0),
            "area_ratio": float(new_areas.sum() / old_areas.sum()),
            "crust_volume_residual_fraction": float(
                abs(np.dot(geometry_response.collect("crustal_thickness_m"), new_areas) - np.dot(old_hc, old_areas))
                / max(abs(np.dot(old_hc, old_areas)), 1.0)
            ),
            "mantle_lithosphere_volume_residual_fraction": float(
                abs(np.dot(geometry_response.collect("mantle_lithosphere_thickness_m"), new_areas) - np.dot(old_hm, old_areas))
                / max(abs(np.dot(old_hm, old_areas)), 1.0)
            ),
            "continental_material_volume_residual_fraction": float(
                abs(np.dot(geometry_response.collect("continental_material_m"), new_areas) - np.dot(old_continental, old_areas))
                / max(abs(np.dot(old_continental, old_areas)), 1.0)
            ),
            "time_integrated": False,
            "slab_body_included": False,
        },
        "accounting_checks": checks,
        "accounting_pass": bool(all(checks.values())),
        "continental_component_demand_quantiles_n_per_m": np.quantile(finite, [0.5, 0.9, 0.95, 0.99]).tolist(),
    }
    if include_demand:
        result["component_demand_n_per_m"] = component_demand.tolist()
        stress = shell.stress_resultant_n_per_m
        matrices = np.empty((len(stress), 2, 2), dtype=float)
        matrices[:, 0, 0] = stress[:, 0]
        matrices[:, 1, 1] = stress[:, 1]
        matrices[:, 0, 1] = matrices[:, 1, 0] = stress[:, 2]
        _principal_values, principal_vectors = np.linalg.eigh(matrices)
        first, second = thin_shell._tangent_basis(rep_points[used])
        principal_xyz = (
            first * principal_vectors[:, 0, 1, None]
            + second * principal_vectors[:, 1, 1, None]
        )
        result["component_principal_extension_xyz"] = principal_xyz[component_mesh_map].tolist()
    return result


def _drag_energy(world, plate, indices, force):
    points = geometry.normalize(force.points_xyz)[indices]
    basal_force = force.basal_drag_n[indices]
    own_v = np.cross(plate.omega, points) * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    mantle_v = mantle.flow_at(points, world.mantle_centers) * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    coefficient = torque.ASTHENOSPHERE_VISCOSITY_PA_S / torque.SHEAR_ZONE_THICKNESS_M
    plate_work_rate = float(np.sum(basal_force * own_v))
    mantle_reaction_work_rate = float(np.sum((-basal_force) * mantle_v))
    dissipation_rate = float(coefficient * np.sum(np.sum((mantle_v - own_v) ** 2, axis=1) * force.areas_m2[indices]))
    residual = abs(plate_work_rate + mantle_reaction_work_rate + dissipation_rate) / max(
        abs(plate_work_rate), abs(mantle_reaction_work_rate), dissipation_rate, 1.0
    )
    return {
        "plate_work_rate_w": plate_work_rate,
        "mantle_reaction_work_rate_w": mantle_reaction_work_rate,
        "viscous_dissipation_rate_w": dissipation_rate,
        "energy_residual_fraction": float(residual),
        "reaction_force_n": (-basal_force).sum(axis=0).tolist(),
        "reaction_torque_nm": np.cross(points * lithosphere.PLANET_RADIUS_M, -basal_force).sum(axis=0).tolist(),
    }


def run(seed=SEED, density=DENSITY):
    started = time.perf_counter()
    world = world_mod.generate_world(seed, node_density=density, cratons_first=True, synthetic_relief=False)
    world.simulate_climate_biomes = False
    cratons.ensure_initialized(world)
    world_mod.step_world(world, STARTUP_STEP_YEARS)
    area, plate, indices = _largest_continental_component(world)
    field = _production_force_field(world, plate)
    expected_points = geometry.normalize(plate.surface_nodes().world_xyz)
    if len(field.points_xyz) != len(expected_points) or np.max(np.linalg.norm(expected_points - geometry.normalize(field.points_xyz), axis=1)) > 1e-10:
        raise ValueError("production force order does not match selected parent plate topology")

    component_area = float(area)
    parent_area = float(plate.node_areas_m2().sum())
    sources = {
        "basal_drag_only": field.basal_drag_n,
        "full_production_force": field.total_n,
        "production_without_basal_drag": field.total_n - field.basal_drag_n,
    }
    cases = {}
    for name, source_force in sources.items():
        cases[name] = [_plate_shell(plate, indices, field.points_xyz, source_force, cap) for cap in MESH_CAPS]
    all_accounting_pass = all(row["accounting_pass"] for rows in cases.values() for row in rows)
    drag_accounting_pass = _drag_energy(world, plate, indices, field)["energy_residual_fraction"] <= 0.05
    return {
        "schema_version": 1,
        "title": "RFC 0010 Gate 3b parent-shell accounting preflight with continental-component evaluation",
        "selection_rule": "Largest edge-connected continental component by represented area after the frozen 2 Myr startup step, before force inspection.",
        "seed": seed,
        "node_density": density,
        "elapsed_years": float(world.elapsed_years),
        "plate_id": int(plate.plate_id),
        "selected_component_nodes": int(len(indices)),
        "production_plate_nodes": int(plate.node_count()),
        "selected_component_area_km2": component_area / 1e6,
        "shell_domain": "entire parent production plate; continental component is the demand-evaluation mask",
        "production_plate_area_km2": parent_area / 1e6,
        "component_area_fraction_of_parent": component_area / parent_area,
        "basal_drag_work_account": _drag_energy(world, plate, np.arange(plate.node_count()), field),
        "basal_drag_energy_pass": bool(drag_accounting_pass),
        "all_source_mesh_accounting_pass": bool(all_accounting_pass),
        "accounting_gate_pass": bool(all_accounting_pass and drag_accounting_pass),
        "history_eligible_by_accounting": bool(all_accounting_pass and drag_accounting_pass),
        "persistence_history_run": False,
        "parent_plate_source_resultants": {
            name: {
                "force_n": np.asarray(force).sum(axis=0).tolist(),
                "torque_nm": np.cross(field.points_xyz * lithosphere.PLANET_RADIUS_M, force).sum(axis=0).tolist(),
            } for name, force in sources.items()
        },
        "cases": cases,
        "frozen_cutoffs": {
            "production_plate_area_coverage_min": 0.95,
            "force_transfer_residual_max": 0.01,
            "torque_transfer_residual_max": 0.01,
            "basal_drag_energy_residual_max": 0.05,
            "shell_equilibrium_residual_max": 0.01,
            "shell_energy_residual_max": 0.05,
        },
        "interpretation": "Initial-state accounting only; accounting gate may be passed, but no Gate 2 capacity score, rift-potential score, or persistence history is included.",
        "wall_time_seconds": time.perf_counter() - started,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--density", type=float, default=DENSITY)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.seed, args.density)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    compact = {key: value for key, value in report.items() if key != "cases"}
    compact["cases"] = {name: [{k: v for k, v in row.items() if k not in ("source_force_n", "transferred_force_n", "source_torque_nm", "transferred_torque_nm")} for row in rows] for name, rows in report["cases"].items()}
    print(json.dumps({"out": str(args.out), **compact}, sort_keys=True))


if __name__ == "__main__":
    main()

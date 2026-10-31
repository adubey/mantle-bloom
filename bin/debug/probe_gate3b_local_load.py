"""Probe a coarse offline shell-response field for the Gate 3b plate plot.

Usage:
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
        bin/debug/probe_gate3b_local_load.py \
        --out analysis/issue699/coarse-shell-response-seed11-density4.json

The deterministic coarse meshes aggregate the existing per-node plate force field onto
spherical representative cells and pass it through the reduced diagnostic membrane solver.
This does not mutate a world or implement production stress, damage, or rifting.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from app import cratons, geometry, lithosphere, thin_shell, world as world_mod
import gate3b_diagnostic_common as shared
import validate_thin_shell_plate_meshes as mesh_validation


SEED = 11
NODE_DENSITY = 4.0
MESH_CAPS = (512, 1024)
EFFECTIVE_MODULUS_PA = 70e9
CAPACITY_PLACEHOLDER_N_PER_M = 1e10


def _source_force_field(world, plate):
    return shared.production_force_field(world, plate)


def _farthest_representatives(points: np.ndarray, areas: np.ndarray, count: int) -> np.ndarray:
    """Select reproducible spherical cell medoids, starting at the largest-area cell."""
    count = min(count, len(points))
    selected = np.empty(count, dtype=np.int64)
    selected[0] = int(np.argmax(areas))
    nearest_distance = np.full(len(points), np.inf)
    for index in range(1, count):
        delta = points - points[selected[index - 1]]
        nearest_distance = np.minimum(nearest_distance, np.einsum("ij,ij->i", delta, delta))
        nearest_distance[selected[:index]] = -1.0
        selected[index] = int(np.argmax(nearest_distance))
    return selected


def _aggregate_to_mesh(plate, field, target_count: int) -> dict:
    points = geometry.normalize(np.asarray(field.points_xyz, dtype=float))
    areas = np.asarray(field.areas_m2, dtype=float)
    representatives = _farthest_representatives(points, areas, target_count)
    representative_points = points[representatives]
    assignments = cKDTree(representative_points).query(points, k=1)[1]
    representative_areas = np.bincount(assignments, weights=areas, minlength=len(representatives))
    total_force = np.asarray(field.total_n, dtype=float)
    aggregated_force = np.zeros((len(representatives), 3))
    np.add.at(aggregated_force, assignments, total_force)
    thickness_source = np.maximum(
        plate.collect("crustal_thickness_m") + plate.collect("mantle_lithosphere_thickness_m"), 1.0
    )
    thickness_sum = np.bincount(assignments, weights=areas * thickness_source, minlength=len(representatives))
    representative_thickness = thickness_sum / np.maximum(representative_areas, np.finfo(float).tiny)

    triangles, rejected, candidates = mesh_validation._candidate_triangles(plate, representative_points)
    component_labels = thin_shell._mesh_component_labels(triangles)
    component_sizes = np.bincount(component_labels)
    largest_component = int(np.argmax(component_sizes))
    discarded_triangles = int(len(triangles) - component_sizes[largest_component])
    triangles = triangles[component_labels == largest_component]
    used = np.unique(triangles)
    remap = np.full(len(representatives), -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    local_triangles = remap[triangles]
    remaining_components = thin_shell._mesh_component_sizes(len(used), local_triangles)
    if len(remaining_components) != 1 or remaining_components[0] == 0:
        raise ValueError(f"largest coarse shell component is disconnected: {remaining_components.tolist()}")

    source_force = total_force.sum(axis=0)
    transferred_force = aggregated_force[used].sum(axis=0)
    source_torque_nodes = np.cross(points * lithosphere.PLANET_RADIUS_M, total_force)
    transferred_torque_nodes = np.cross(
        representative_points[used] * lithosphere.PLANET_RADIUS_M, aggregated_force[used]
    )
    force_denominator = max(float(np.linalg.norm(total_force, axis=1).sum()), 1.0)
    torque_denominator = max(float(np.linalg.norm(source_torque_nodes, axis=1).sum()), 1.0)

    result = thin_shell.solve_thin_shell(
        points_xyz=representative_points[used],
        areas_m2=representative_areas[used],
        applied_force_n=aggregated_force[used],
        thickness_m=representative_thickness[used],
        young_modulus_pa=EFFECTIVE_MODULUS_PA,
        extensional_capacity_resultant_n_per_m=np.full(len(used), CAPACITY_PLACEHOLDER_N_PER_M),
        mesh_triangles=local_triangles,
        tolerance=1e-6,
        max_iterations=5000,
    )
    residual_fraction = result.equilibrium_residual_n / max(
        float(np.linalg.norm(result.equilibrated_force_n)), 1.0
    )
    demand = result.max_incremental_extensional_resultant_n_per_m
    # Map the solved shell nodes onto every source cell by nearest spherical representative.
    mapped = cKDTree(representative_points[used]).query(points, k=1)[1]
    full_demand = demand[mapped]
    lat_lon = np.column_stack((
        np.degrees(np.arcsin(np.clip(points[:, 2], -1.0, 1.0))),
        np.degrees(np.arctan2(points[:, 1], points[:, 0])),
    ))
    finite = demand[np.isfinite(demand)]
    quantiles = np.quantile(finite, [0.5, 0.9, 0.95, 0.99])
    return {
        "target_representatives": int(target_count),
        "representatives": int(len(representatives)),
        "mesh_nodes": int(len(used)),
        "candidate_triangles": int(candidates),
        "rejected_triangles": int(rejected),
        "accepted_triangles": int(len(triangles)),
        "triangle_component_sizes_before_largest_component_filter": component_sizes.tolist(),
        "discarded_triangles": discarded_triangles,
        "shell_components": 1,
        "area_coverage": float(representative_areas[used].sum() / areas.sum()),
        "source_force_n": source_force.tolist(),
        "transferred_force_n": transferred_force.tolist(),
        "force_transfer_residual_fraction": float(
            np.linalg.norm(source_force - transferred_force) / force_denominator
        ),
        "source_torque_nm": source_torque_nodes.sum(axis=0).tolist(),
        "transferred_torque_nm": transferred_torque_nodes.sum(axis=0).tolist(),
        "torque_transfer_residual_fraction": float(
            np.linalg.norm(source_torque_nodes.sum(axis=0) - transferred_torque_nodes.sum(axis=0))
            / torque_denominator
        ),
        "solver_stop_code": int(result.solver_stop_code),
        "equilibrium_residual_n": float(result.equilibrium_residual_n),
        "equilibrium_residual_fraction": float(residual_fraction),
        "demand_quantiles_n_per_m": {
            name: float(value) for name, value in zip(("p50", "p90", "p95", "p99"), quantiles)
        },
        "mesh_lat_lon_degrees": lat_lon[representatives[used]].tolist(),
        "mesh_demand_n_per_m": demand.tolist(),
        "_plate_demand_n_per_m": full_demand,
        "_plate_cell_area_m2": areas,
        "wall_time_seconds": {},
    }


def _rank_correlation(first: np.ndarray, second: np.ndarray) -> float:
    from scipy.stats import rankdata, spearmanr

    return float(spearmanr(rankdata(first), rankdata(second)).statistic)


def run(seed: int = SEED, node_density: float = NODE_DENSITY) -> dict:
    started = time.perf_counter()
    world = world_mod.generate_world(
        seed, node_density=node_density, cratons_first=True, synthetic_relief=False
    )
    cratons.ensure_initialized(world)
    candidates = [plate for plate in world.plates if plate.crust_type == "continental"]
    plate = max(candidates, key=lambda item: float(item.node_areas_m2().sum()))
    field = _source_force_field(world, plate)
    nodes = plate.surface_nodes()
    points = geometry.normalize(np.asarray(field.points_xyz, dtype=float))
    alignment = float(np.max(np.linalg.norm(points - geometry.normalize(nodes.world_xyz), axis=1)))
    if alignment > 1e-10 or len(points) != plate.node_count():
        raise ValueError("per-node force order does not match production surface-node order")

    results = []
    for cap in MESH_CAPS:
        mesh_started = time.perf_counter()
        results.append(_aggregate_to_mesh(plate, field, cap))
        results[-1]["wall_time_seconds"]["coarse_mesh_and_solve"] = time.perf_counter() - mesh_started
    common_first = results[0].pop("_plate_demand_n_per_m")
    common_second = results[1].pop("_plate_demand_n_per_m")
    areas = results[0].pop("_plate_cell_area_m2")
    results[1].pop("_plate_cell_area_m2")
    first_cut, second_cut = np.quantile(common_first, 0.9), np.quantile(common_second, 0.9)
    top_first, top_second = common_first >= first_cut, common_second >= second_cut
    intersection = float(areas[top_first & top_second].sum())
    union = float(areas[top_first | top_second].sum())
    return {
        "title": "RFC 0010 Gate 3b coarse shell-response diagnostic",
        "interpretation": (
            "Reduced linear static membrane response to the existing per-node force field. "
            "The mapped principal tensile resultant is a diagnostic ranking field, not strain "
            "rate, failure, damage, or a production rift law."
        ),
        "seed": seed,
        "node_density": node_density,
        "elapsed_years": float(world.elapsed_years),
        "plate_id": int(plate.plate_id),
        "plate_node_count": int(plate.node_count()),
        "plate_area_km2": float(plate.node_areas_m2().sum() / 1e6),
        "plate_neighbors": [int(item.plate_id) for item in plate.get_neighbours(world.plates)],
        "representative_rule": "deterministic farthest-point medoids; nearest-medoid force aggregation",
        "effective_modulus_pa": EFFECTIVE_MODULUS_PA,
        "poisson_ratio": 0.25,
        "source_force_gather": "production torque.gather_per_node_force_field; no world step",
        "maps": results,
        "resolution_comparison": {
            "spearman_rank_correlation": _rank_correlation(common_first, common_second),
            "top_decile_area_iou": intersection / max(union, 1.0),
        },
        "total_wall_time_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--density", type=float, default=NODE_DENSITY)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    output = run(args.seed, args.density)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "out": str(args.out),
        "resolution_comparison": output["resolution_comparison"],
        "maps": [{key: value for key, value in item.items() if key not in (
            "mesh_lat_lon_degrees", "mesh_demand_n_per_m", "source_force_n",
            "transferred_force_n", "source_torque_nm", "transferred_torque_nm"
        )} for item in output["maps"]],
        "total_wall_time_seconds": output["total_wall_time_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

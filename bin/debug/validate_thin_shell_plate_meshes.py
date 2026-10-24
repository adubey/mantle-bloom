"""Validate explicit thin-shell meshes against generated plate-cell topology.

Usage (from the repository root):
    PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
      bin/debug/validate_thin_shell_plate_meshes.py \
      --seed 7 --num-plates 12 --node-density 0.5 \
      --out analysis/issue487/representative-mesh-summary.json

The script selects the largest generated plate with an interior boundary loop and the
largest one-loop plate. It clips a Delaunay triangulation of cell centres against exact plate
containment probes. It tests a manufactured load on both the original curved mesh and a shallow
cap mapped from the same topology, sweeps localized and paired clamp placements, and checks an
independent curved-sphere solution with uniform isotropic membrane tension. A final response
sweep evaluates a parameterized phenomenological capacity; its reference threshold and weakness
mask are design inputs, not calibrated or physical strength values.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import Delaunay

from app import geometry, lithosphere, lithosphere_plate, thin_shell

_PROBE_WEIGHTS = np.array(
    [
        [1 / 3, 1 / 3, 1 / 3],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
        [0.75, 0.25, 0.0],
        [0.75, 0.0, 0.25],
        [0.0, 0.75, 0.25],
        [0.0, 0.25, 0.75],
    ]
)
_TARGET_RESULTANT_N_PER_M = 1e10


def _candidate_triangles(plate, points: np.ndarray) -> tuple[np.ndarray, int, int]:
    """Delaunay-triangulate cell centres and reject elements whose interior probes leave the plate."""
    center = geometry.normalize(points.mean(axis=0))
    denominator = points @ center
    if np.any(denominator <= 0.0):
        raise ValueError(f"plate {plate.plate_id} does not fit in one tangent hemisphere")
    first, second = thin_shell._tangent_basis(center[None, :])
    xy = np.column_stack((points @ first[0] / denominator, points @ second[0] / denominator))
    candidates = Delaunay(xy).simplices
    candidate_xyz = points[candidates]
    probes = np.sum(candidate_xyz[:, None, :, :] * _PROBE_WEIGHTS[None, :, :, None], axis=2)
    probes = geometry.normalize(probes.reshape(-1, 3)).reshape(len(candidates), -1, 3)
    accepted = plate.contains_batch(probes.reshape(-1, 3)).reshape(len(candidates), -1).all(axis=1)
    return candidates[accepted], int(np.count_nonzero(~accepted)), int(len(candidates))


def _boundary_edges(triangles: np.ndarray) -> tuple[list[tuple[int, int]], int, int]:
    """Return directed exterior edges, count non-manifold edges, and count boundary loops."""
    counts: Counter[tuple[int, int]] = Counter()
    directions: dict[tuple[int, int], tuple[int, int]] = {}
    for triangle in triangles:
        for first, second in ((triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0])):
            key = tuple(sorted((int(first), int(second))))
            counts[key] += 1
            directions[key] = (int(first), int(second))
    exterior = [directions[edge] for edge, count in counts.items() if count == 1]
    nonmanifold = sum(count > 2 for count in counts.values())
    if not exterior:
        return exterior, nonmanifold, 0
    vertices = np.unique(np.asarray(exterior, dtype=np.int64))
    local = np.full(int(triangles.max()) + 1, -1, dtype=np.int64)
    local[vertices] = np.arange(len(vertices))
    edges = local[np.asarray(exterior, dtype=np.int64)]
    graph = coo_matrix(
        (np.ones(2 * len(edges)), (np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]])),
        shape=(len(vertices), len(vertices)),
    ).tocsr()
    loops, _ = connected_components(graph, directed=False)
    return exterior, int(nonmanifold), int(loops)


def _shallow_mesh(points: np.ndarray, areas_m2: np.ndarray, triangles: np.ndarray):
    """Map representative projected node distribution to a shallow cap, retaining its triangles."""
    center = geometry.normalize(points.mean(axis=0))
    first, second = thin_shell._tangent_basis(center[None, :])
    denominator = points @ center
    xy = np.column_stack((points @ first[0] / denominator, points @ second[0] / denominator))
    xy -= xy.mean(axis=0)
    span = float(np.ptp(xy, axis=0).max())
    scale = 0.1 / span
    xy *= scale
    cap_points = np.column_stack((xy, np.sqrt(1.0 - np.sum(xy**2, axis=1))))
    first_edge = xy[triangles[:, 1]] - xy[triangles[:, 0]]
    second_edge = xy[triangles[:, 2]] - xy[triangles[:, 0]]
    signed_area = first_edge[:, 0] * second_edge[:, 1] - first_edge[:, 1] * second_edge[:, 0]
    triangles = triangles.copy()
    clockwise = signed_area < 0.0
    triangles[clockwise, 1], triangles[clockwise, 2] = triangles[clockwise, 2].copy(), triangles[clockwise, 1].copy()
    return cap_points, areas_m2 * scale**2, xy, triangles


def _manufactured_result(points: np.ndarray, areas_m2: np.ndarray, xy: np.ndarray, triangles: np.ndarray) -> dict:
    """Load every mesh boundary with an analytic constant uniaxial membrane resultant."""
    exterior, nonmanifold, loops = _boundary_edges(triangles)
    if nonmanifold:
        raise ValueError(f"explicit mesh has {nonmanifold} non-manifold edges")
    if loops == 0:
        raise ValueError("explicit mesh has no boundary loop")
    force = np.zeros_like(points)
    xy_m = xy * lithosphere.PLANET_RADIUS_M
    for first, second in exterior:
        edge = xy_m[second] - xy_m[first]
        length = float(np.linalg.norm(edge))
        outward_x = edge[1] / length  # CCW element winding: the domain lies to the left.
        nodal_load = _TARGET_RESULTANT_N_PER_M * outward_x * np.array([1.0, 0.0, 0.0]) * length * 0.5
        force[first] += nodal_load
        force[second] += nodal_load

    result = thin_shell.solve_thin_shell(
        points,
        areas_m2,
        force,
        thickness_m=100_000.0,
        young_modulus_pa=70e9,
        extensional_capacity_resultant_n_per_m=np.full(len(points), 2.0 * _TARGET_RESULTANT_N_PER_M),
        mesh_triangles=triangles,
    )
    mean = float(np.mean(result.max_incremental_extensional_resultant_n_per_m))
    return {
        "mesh_boundary_loops": loops,
        "boundary_edges": len(exterior),
        "mean_resultant_n_per_m": mean,
        "relative_error": mean / _TARGET_RESULTANT_N_PER_M - 1.0,
        "equilibrium_residual_over_load": result.equilibrium_residual_n / float(np.linalg.norm(result.equilibrated_force_n)),
        "solver_stop_code": result.solver_stop_code,
    }


def _spherical_uniform_tension_case(theta_count: int, phi_count: int) -> dict:
    """Check a spherical membrane with uniform isotropic tension and analytic edge traction."""
    theta = np.linspace(0.12, 1.25, theta_count)
    phi = np.linspace(-0.72, 0.72, phi_count)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    points = np.column_stack(
        (
            np.sin(theta_grid.ravel()) * np.cos(phi_grid.ravel()),
            np.sin(theta_grid.ravel()) * np.sin(phi_grid.ravel()),
            np.cos(theta_grid.ravel()),
        )
    )
    def node(i: int, j: int) -> int:
        return i * phi_count + j
    triangles = []
    for i in range(theta_count - 1):
        for j in range(phi_count - 1):
            triangles.extend(
                ((node(i, j), node(i + 1, j), node(i + 1, j + 1)),
                 (node(i, j), node(i + 1, j + 1), node(i, j + 1)))
            )
    triangles = np.asarray(triangles, dtype=np.int64)
    for index, triangle in enumerate(triangles):
        normal = np.cross(points[triangle[1]] - points[triangle[0]], points[triangle[2]] - points[triangle[0]])
        if np.dot(normal, points[triangle].mean(axis=0)) < 0.0:
            triangles[index, 1], triangles[index, 2] = triangles[index, 2], triangles[index, 1]

    edge_count: Counter[tuple[int, int]] = Counter()
    edge_direction: dict[tuple[int, int], tuple[int, int]] = {}
    for triangle in triangles:
        for first, second in ((triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0])):
            edge = tuple(sorted((int(first), int(second))))
            edge_count[edge] += 1
            edge_direction[edge] = (int(first), int(second))

    force = np.zeros_like(points)
    for edge, count in edge_count.items():
        if count != 1:
            continue
        first, second = edge_direction[edge]
        midpoint = geometry.normalize(points[first] + points[second])
        edge_vector = points[second] - points[first]
        tangent = geometry.normalize(edge_vector - np.dot(edge_vector, midpoint) * midpoint)
        outward_conormal = geometry.normalize(np.cross(tangent, midpoint))
        edge_angle = float(np.arccos(np.clip(np.dot(points[first], points[second]), -1.0, 1.0)))
        edge_length = edge_angle * lithosphere.PLANET_RADIUS_M
        nodal_force = _TARGET_RESULTANT_N_PER_M * outward_conormal * edge_length * 0.5
        force[first] += nodal_force
        force[second] += nodal_force

    theta_width = np.full(theta_count, theta[1] - theta[0])
    phi_width = np.full(phi_count, phi[1] - phi[0])
    theta_width[[0, -1]] *= 0.5
    phi_width[[0, -1]] *= 0.5
    node_areas = (
        lithosphere.PLANET_RADIUS_M**2
        * np.sin(theta_grid)
        * theta_width[:, None]
        * phi_width[None, :]
    ).ravel()
    result = thin_shell.solve_thin_shell(
        points,
        node_areas,
        force,
        thickness_m=100_000.0,
        young_modulus_pa=70e9,
        extensional_capacity_resultant_n_per_m=np.full(
            len(points), 2.0 * _TARGET_RESULTANT_N_PER_M
        ),
        mesh_triangles=triangles,
        tolerance=1e-10,
    )
    normalized_error = (
        result.max_incremental_extensional_resultant_n_per_m / _TARGET_RESULTANT_N_PER_M - 1.0
    )
    return {
        "theta_nodes": theta_count,
        "phi_nodes": phi_count,
        "triangle_count": len(triangles),
        "analytic_isotropic_resultant_n_per_m": _TARGET_RESULTANT_N_PER_M,
        "mean_resultant_n_per_m": float(
            np.mean(result.max_incremental_extensional_resultant_n_per_m)
        ),
        "mean_relative_error": float(np.mean(normalized_error)),
        "p95_absolute_relative_error": float(np.quantile(np.abs(normalized_error), 0.95)),
        "max_absolute_relative_error": float(np.max(np.abs(normalized_error))),
        "equilibrium_residual_over_load": result.equilibrium_residual_n
        / float(np.linalg.norm(result.equilibrated_force_n)),
        "solver_stop_code": result.solver_stop_code,
    }


def _curved_sphere_manufactured_benchmark() -> dict:
    """Refine a broad spherical patch against the uniform-tension membrane solution."""
    cases = [
        _spherical_uniform_tension_case(theta_count, phi_count)
        for theta_count, phi_count in ((5, 9), (9, 17), (17, 33))
    ]
    if any(case["equilibrium_residual_over_load"] >= 1e-5 for case in cases):
        raise AssertionError("curved-shell manufactured solve residual exceeds 1e-5 of load")
    if any(
        later["p95_absolute_relative_error"] >= earlier["p95_absolute_relative_error"]
        for earlier, later in zip(cases, cases[1:])
    ):
        raise AssertionError("curved-shell stress error must decrease under mesh refinement")
    if cases[-1]["p95_absolute_relative_error"] >= 0.03:
        raise AssertionError("fine curved-shell p95 stress error must be below 3%")
    return {
        "analytic_solution": "uniform isotropic membrane resultant N on a sphere; normal pressure 2N/R is supported and omitted from in-plane solve",
        "edge_condition": "free displacement with prescribed analytic conormal traction",
        "mesh_refinement": cases,
    }


def _phenomenological_capacity_sweep(resultant: np.ndarray, xy: np.ndarray) -> list[dict]:
    """Check threshold and prescribed weak-band behavior for an explicit capacity rule."""
    x = xy[:, 0]
    weak_band = np.abs(x - np.median(x)) <= 0.08 * float(np.ptp(x))
    weakening_fraction = 0.5
    cases = []
    for reference_capacity in (1e10, 2e10, 4e10):
        background = np.full(len(resultant), reference_capacity)
        capacity = background * np.where(weak_band, 1.0 - weakening_fraction, 1.0)
        failed = resultant >= capacity
        weak_failed = int((failed & weak_band).sum())
        strong_failed = int((failed & ~weak_band).sum())
        cases.append(
            {
                "reference_capacity_n_per_m": reference_capacity,
                "weakening_fraction": weakening_fraction,
                "weak_band_fraction_of_nodes": float(weak_band.mean()),
                "failed_nodes_total": weak_failed + strong_failed,
                "failed_nodes_in_weak_band": weak_failed,
                "failed_fraction_in_weak_band": weak_failed / int(weak_band.sum()),
                "failed_nodes_outside_weak_band": strong_failed,
                "failed_fraction_outside_weak_band": strong_failed / int((~weak_band).sum()),
            }
        )
    failed_counts = [case["failed_nodes_total"] for case in cases]
    if any(larger > smaller for smaller, larger in zip(failed_counts, failed_counts[1:])):
        raise AssertionError("increasing reference capacity must not increase failed-node count")
    middle = cases[1]
    if middle["failed_fraction_in_weak_band"] <= middle["failed_fraction_outside_weak_band"]:
        raise AssertionError("the prescribed weak band must fail more often than the background")
    return cases


def _curved_mesh_sensitivity(
    points: np.ndarray, areas_m2: np.ndarray, triangles: np.ndarray
) -> dict:
    """Sweep free edges and four localized clamp bands on the original spherical mesh."""
    center = geometry.normalize(points.mean(axis=0))
    first, second = thin_shell._tangent_basis(center[None, :])
    denominator = points @ center
    if np.any(denominator <= 0.0):
        raise ValueError("curved plate mesh does not fit in one tangent hemisphere")
    xy = np.column_stack((points @ first[0] / denominator, points @ second[0] / denominator))
    signed_area = (
        (xy[triangles[:, 1], 0] - xy[triangles[:, 0], 0])
        * (xy[triangles[:, 2], 1] - xy[triangles[:, 0], 1])
        - (xy[triangles[:, 1], 1] - xy[triangles[:, 0], 1])
        * (xy[triangles[:, 2], 0] - xy[triangles[:, 0], 0])
    )
    triangles = triangles.copy()
    clockwise = signed_area < 0.0
    triangles[clockwise, 1], triangles[clockwise, 2] = (
        triangles[clockwise, 2].copy(),
        triangles[clockwise, 1].copy(),
    )
    exterior, nonmanifold, loops = _boundary_edges(triangles)
    if nonmanifold or loops == 0:
        raise ValueError(f"curved explicit mesh has {nonmanifold} non-manifold edges and {loops} loops")

    force = np.zeros_like(points)
    boundary_nodes = np.unique(np.asarray(exterior, dtype=np.int64))
    for start, end in exterior:
        edge = (xy[end] - xy[start]) * lithosphere.PLANET_RADIUS_M
        length = float(np.linalg.norm(edge))
        projected_x_traction = edge[1] / length
        midpoint = geometry.normalize(points[start] + points[end])
        tangent_x = first[0] - np.dot(first[0], midpoint) * midpoint
        tangent_x = geometry.normalize(tangent_x)
        nodal_force = _TARGET_RESULTANT_N_PER_M * projected_x_traction * tangent_x * length * 0.5
        force[start] += nodal_force
        force[end] += nodal_force

    span = np.ptp(xy, axis=0)
    band_width = 0.12 * float(span.max())
    cases = []
    free_result = None
    edge_conditions = (
        "free", "west", "east", "south", "north", "west+east", "south+north", "all_edges"
    )
    for direction in edge_conditions:
        fixed = np.zeros(len(points), dtype=bool)
        for edge in direction.split("+"):
            if edge == "west":
                fixed |= xy[:, 0] < xy[:, 0].min() + band_width
            elif edge == "east":
                fixed |= xy[:, 0] > xy[:, 0].max() - band_width
            elif edge == "south":
                fixed |= xy[:, 1] < xy[:, 1].min() + band_width
            elif edge == "north":
                fixed |= xy[:, 1] > xy[:, 1].max() - band_width
            elif edge == "all_edges":
                fixed |= np.isin(np.arange(len(points)), boundary_nodes)
        if direction != "free":
            fixed &= np.isin(np.arange(len(points)), boundary_nodes)
        result = thin_shell.solve_thin_shell(
            points,
            areas_m2,
            force,
            thickness_m=100_000.0,
            young_modulus_pa=70e9,
            extensional_capacity_resultant_n_per_m=np.full(
                len(points), 2.0 * _TARGET_RESULTANT_N_PER_M
            ),
            mesh_triangles=triangles,
            fixed_node_mask=fixed,
        )
        if direction == "free":
            free_result = result
        resultant = result.max_incremental_extensional_resultant_n_per_m
        if free_result is None:
            raise RuntimeError("free-edge reference result was not computed")
        free_mean = float(np.mean(free_result.max_incremental_extensional_resultant_n_per_m))
        cases.append(
            {
                "edge_condition": direction,
                "clamped_boundary_nodes": int(fixed.sum()),
                "mean_resultant_n_per_m": float(np.mean(resultant)),
                "mean_change_from_free_fraction": float(np.mean(resultant) / free_mean - 1.0),
                "p95_resultant_n_per_m": float(np.quantile(resultant, 0.95)),
                "max_resultant_n_per_m": float(np.max(resultant)),
                "equilibrium_residual_over_load": result.equilibrium_residual_n
                / float(np.linalg.norm(result.equilibrated_force_n)),
                "solver_stop_code": result.solver_stop_code,
            }
        )
    return {
        "boundary_condition_policy": "Free-edge Neumann boundary is the diagnostic baseline because the single-plate input has no displacement anchor; clamps are sensitivity probes, not inferred physical plate boundaries. Coupled neighboring plates would require a multi-plate solve.",
        "projected_span_radians": span.tolist(),
        "mesh_boundary_loops": loops,
        "boundary_edges": len(exterior),
        "manufactured_resultant_n_per_m": _TARGET_RESULTANT_N_PER_M,
        "capacity_sensitivity": _phenomenological_capacity_sweep(
            free_result.max_incremental_extensional_resultant_n_per_m, xy
        ),
        "edge_cases": cases,
    }


def _validate_plate(plate) -> dict:
    """Measure hull spill and validate a shallow manufactured solve on the clipped mesh."""
    nodes = plate.surface_nodes()
    source_points = nodes.world_xyz
    hull = thin_shell._spherical_surface_triangles(source_points)
    hull_centroids = geometry.normalize(source_points[hull].mean(axis=1))
    hull_outside = ~plate.contains_batch(hull_centroids)
    candidates, rejected, total_candidates = _candidate_triangles(plate, source_points)
    used = np.unique(candidates)
    remap = np.full(len(source_points), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    triangles = remap[candidates]
    curved = _curved_mesh_sensitivity(source_points[used], nodes.area_m2[used], triangles)
    points, areas, xy, triangles = _shallow_mesh(source_points[used], nodes.area_m2[used], triangles)
    manufactured = _manufactured_result(points, areas, xy, triangles)
    return {
        "plate_id": plate.plate_id,
        "crust_type": plate.crust_type,
        "node_count": plate.node_count(),
        "actual_boundary_loops": len(plate.boundary_loops_world()),
        "convex_hull_triangles": len(hull),
        "convex_hull_centroids_outside_plate": int(hull_outside.sum()),
        "convex_hull_outside_fraction": float(hull_outside.mean()),
        "delaunay_candidate_triangles": total_candidates,
        "triangles_rejected_by_territory_probes": rejected,
        "explicit_triangles": len(triangles),
        "explicit_mesh_nodes": len(points),
        "original_curved_mesh_sensitivity": curved,
        **manufactured,
    }


def run(seed: int, num_plates: int, node_density: float) -> dict:
    """Generate a deterministic plate set and validate its largest hole-bearing and one-loop plates."""
    plates = lithosphere_plate.generate_plates(
        seed, num_plates=num_plates, node_density=node_density, cratons_first=True
    )
    hole_plates = [plate for plate in plates if len(plate.boundary_loops_world()) > 1]
    single_loop_plates = [plate for plate in plates if len(plate.boundary_loops_world()) == 1]
    if not hole_plates or not single_loop_plates:
        raise ValueError("generated sample needs at least one plate with a hole and one single-loop plate")
    selected = [
        max(hole_plates, key=lambda plate: plate.node_count()),
        max(single_loop_plates, key=lambda plate: plate.node_count()),
    ]
    return {
        "seed": seed,
        "num_plates": num_plates,
        "node_density": node_density,
        "mapping": "Generated cell-centre topology mapped to a shallow cap for manufactured-load testing",
        "target_resultant_n_per_m": _TARGET_RESULTANT_N_PER_M,
        "curved_shell_manufactured_benchmark": _curved_sphere_manufactured_benchmark(),
        "boundary_condition_selection": "Free in-plane displacement with loads applied as nodal forces is the single-plate diagnostic baseline. Fixed-edge cases are sensitivity probes because no neighboring-plate displacement anchor exists in the input. Physical coupling requires solving neighboring plates together.",
        "phenomenological_capacity_model": {
            "formula": "T_cap = T_ref * (1 - weakening_fraction * weak_zone_weight)",
            "weakening_fraction": 0.5,
            "reference_capacity_sweep_n_per_m": [1e10, 2e10, 4e10],
            "weak_zone_weight": "1 inside |projected_x - median(projected_x)| <= 0.08 * projected_x_span; otherwise 0",
            "status": "synthetic response check; parameters are not calibrated or physical",
        },
        "plates": [_validate_plate(plate) for plate in selected],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--num-plates", type=int, default=12)
    parser.add_argument("--node-density", type=float, default=0.5)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = run(args.seed, args.num_plates, args.node_density)
    rendered = json.dumps(report, indent=2) + "\n"
    if args.out is None:
        print(rendered, end="")
    else:
        args.out.write_text(rendered)


if __name__ == "__main__":
    main()

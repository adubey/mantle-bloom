"""Validate explicit thin-shell meshes against generated plate-cell topology.

Usage (from the repository root):
    PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
      bin/debug/validate_thin_shell_plate_meshes.py \
      --seed 7 --num-plates 12 --node-density 0.5 \
      --out analysis/issue487/representative-mesh-summary.json

The script selects the largest generated plate with an interior boundary loop and the
largest one-loop plate. It clips a Delaunay triangulation of cell centres against exact plate
containment probes, then maps each mesh to a shallow spherical patch with the same topology.
Manufactured edge tractions test mesh transfer without treating a broad spherical patch's
global direction field as an analytic uniform-stress solution. This is a mesh diagnostic,
not a replay of production stresses or a material-strength calibration.
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
        tensile_strength_resultant_n_per_m=np.full(len(points), 2.0 * _TARGET_RESULTANT_N_PER_M),
        mesh_triangles=triangles,
    )
    mean = float(np.mean(result.max_tensile_resultant_n_per_m))
    return {
        "mesh_boundary_loops": loops,
        "boundary_edges": len(exterior),
        "mean_resultant_n_per_m": mean,
        "relative_error": mean / _TARGET_RESULTANT_N_PER_M - 1.0,
        "equilibrium_residual_over_load": result.equilibrium_residual_n / float(np.linalg.norm(result.equilibrated_force_n)),
        "solver_stop_code": result.solver_stop_code,
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

#!/usr/bin/env python3
"""Benchmark Gate 2 membrane loads across a shared two-plate interface.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python \
      bin/debug/benchmark_gate2_coupled_shell.py \
      --output analysis/issue487/gate2-coupled-shell-benchmark.json

This is a manufactured diagnostic, not production plate coupling. A curved two-plate cap
has a shared interface and equal-and-opposite interface tractions in the isolated runs.
The joined run assembles both plate regions on shared interface nodes and applies only the
remote outer-edge tractions. It checks global balance, interface reaction continuity, and
recovered extensional resultant against the isolated-plate solves.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from app import geometry, lithosphere, thin_shell

TARGET_N = 1.0e10
RADIUS = lithosphere.PLANET_RADIUS_M


def _mesh(nx_per_plate: int = 8, ny: int = 6):
    """Create a mildly curved, two-region spherical cap with shared interface nodes."""
    x = np.linspace(-0.10, 0.10, 2 * nx_per_plate + 1)
    y = np.linspace(-0.06, 0.06, ny + 1)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    points = geometry.normalize(np.column_stack((xx.ravel(), yy.ravel(), np.ones(xx.size))))

    def node(i: int, j: int) -> int:
        return i * (ny + 1) + j

    triangles = []
    regions = []
    for i in range(2 * nx_per_plate):
        region = 0 if i < nx_per_plate else 1
        for j in range(ny):
            triangles.extend(((node(i, j), node(i + 1, j), node(i + 1, j + 1)),
                              (node(i, j), node(i + 1, j + 1), node(i, j + 1))))
            regions.extend((region, region))
    triangles = np.asarray(triangles, dtype=np.int64)
    regions = np.asarray(regions, dtype=np.int8)
    for k, tri in enumerate(triangles):
        if np.dot(np.cross(points[tri[1]] - points[tri[0]], points[tri[2]] - points[tri[0]]), points[tri].mean(axis=0)) < 0:
            triangles[k, 1], triangles[k, 2] = triangles[k, 2], triangles[k, 1]
    return points, triangles, regions, xx.ravel(), nx_per_plate, ny


def _boundary_edges(triangles: np.ndarray):
    counts: Counter[tuple[int, int]] = Counter()
    directions: dict[tuple[int, int], tuple[int, int]] = {}
    for tri in triangles:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = tuple(sorted((int(a), int(b))))
            counts[key] += 1
            directions[key] = (int(a), int(b))
    return [directions[key] for key, count in counts.items() if count == 1]


def _edge_traction(points: np.ndarray, edges, target_n: float):
    """Apply outward conormal traction for a uniform isotropic membrane resultant."""
    force = np.zeros_like(points)
    for a, b in edges:
        midpoint = geometry.normalize(points[a] + points[b])
        chord = points[b] - points[a]
        tangent = geometry.normalize(chord - np.dot(chord, midpoint) * midpoint)
        conormal = geometry.normalize(np.cross(tangent, midpoint))
        length = float(np.arccos(np.clip(np.dot(points[a], points[b]), -1, 1)) * RADIUS)
        nodal = target_n * conormal * length * 0.5
        force[a] += nodal
        force[b] += nodal
    return force

def _solve(points, triangles, force, thickness, modulus):
    return thin_shell.solve_thin_shell(
        points, np.full(len(points), 1e10), force, thickness_m=thickness,
        young_modulus_pa=modulus,
        extensional_capacity_resultant_n_per_m=np.full(len(points), 2 * TARGET_N),
        mesh_triangles=triangles, tolerance=1e-9,
    )


def _regional_internal_forces(points, triangles, regions, result, thickness, modulus):
    first, second = thin_shell._tangent_basis(points)
    displacement = result.displacement_m
    dofs = np.column_stack((np.sum(displacement * first, axis=1),
                            np.sum(displacement * second, axis=1))).ravel()
    regional = np.zeros((2, len(points), 3))
    for tri, region in zip(triangles, regions, strict=True):
        stiffness, *_ = thin_shell._triangle_element(
            points, tri, first, second, np.full(len(points), thickness),
            np.full(len(points), modulus), 0.25, RADIUS,
        )
        element_dofs = np.array([d for node in tri for d in (2 * node, 2 * node + 1)])
        local_force = stiffness @ dofs[element_dofs]
        local_xyz = local_force.reshape(3, 2)[:, :1] * first[tri] + local_force.reshape(3, 2)[:, 1:] * second[tri]
        np.add.at(regional[int(region)], tri, local_xyz)
    return regional


def run_benchmark() -> dict:
    points, triangles, regions, xcoord, nx_half, ny = _mesh()
    interface = np.zeros(len(points), dtype=bool)
    interface[nx_half * (ny + 1):(nx_half + 1) * (ny + 1)] = True
    interface_inner = interface.copy()
    interface_inner[nx_half * (ny + 1)] = False
    interface_inner[(nx_half + 1) * (ny + 1) - 1] = False
    left_nodes = xcoord <= 1e-12
    right_nodes = xcoord >= -1e-12

    thickness = 100_000.0
    modulus = 70e9
    # Joined solve: only the combined outer boundary receives tractions; the interface is internal.
    outer_edges = _boundary_edges(triangles)
    outer_force = _edge_traction(points, outer_edges, TARGET_N)
    coupled = _solve(points, triangles, outer_force, thickness, modulus)

    isolated = []
    isolated_node_maps = []
    isolated_forces = []
    for region, mask in ((0, left_nodes), (1, right_nodes)):
        global_nodes = np.flatnonzero(mask)
        local_map = np.full(len(points), -1, dtype=np.int64)
        local_map[global_nodes] = np.arange(len(global_nodes))
        region_triangles = triangles[regions == region]
        local_triangles = local_map[region_triangles]
        local_points = points[global_nodes]
        local_force = _edge_traction(local_points, _boundary_edges(local_triangles), TARGET_N)
        isolated.append(_solve(local_points, local_triangles, local_force, thickness, modulus))
        isolated_forces.append(local_force)
        isolated_node_maps.append(global_nodes)

    # The analytic sphere carries the normal Laplace pressure p=2N/R. It is omitted from
    # the in-plane solve, but include its integrated reaction in the global force audit.
    boundary_force = outer_force.sum(axis=0)
    node_area = np.zeros(len(points))
    first, second = thin_shell._tangent_basis(points)
    for tri in triangles:
        _, _, _, area, _, _ = thin_shell._triangle_element(
            points, tri, first, second, np.full(len(points), thickness),
            np.full(len(points), modulus), 0.25, RADIUS,
        )
        node_area[tri] += area / 3.0
    pressure_reaction = (2.0 * TARGET_N / RADIUS) * node_area[:, None] * points
    force_balance = boundary_force + pressure_reaction.sum(axis=0)
    torque = np.sum(np.cross(points * RADIUS, outer_force), axis=0)
    force_scale = float(np.sum(np.linalg.norm(outer_force, axis=1)))
    torque_scale = force_scale * RADIUS

    # Isolated cut-face tractions are equal and opposite away from the outer corners.
    isolated_interface_force = np.zeros((2, len(points), 3))
    for region, (nodes, local_force) in enumerate(zip(isolated_node_maps, isolated_forces, strict=True)):
        isolated_interface_force[region, nodes] = local_force
    isolated_interface_imbalance = isolated_interface_force[0, interface_inner] + isolated_interface_force[1, interface_inner]
    interface_force_scale = TARGET_N * float(np.mean(np.diff(np.linspace(-0.06, 0.06, ny + 1))) * RADIUS)
    isolated_traction_imbalance = float(np.max(np.linalg.norm(isolated_interface_imbalance, axis=1)) / interface_force_scale)
    regional_force = _regional_internal_forces(points, triangles, regions, coupled, thickness, modulus)
    interface_imbalance = regional_force[0, interface_inner] + regional_force[1, interface_inner]
    # N times each node's interface tributary length is the force scale for local continuity.
    interface_indices = np.flatnonzero(interface_inner)
    y = np.linspace(-0.06, 0.06, ny + 1)
    dy = np.diff(y) * RADIUS
    tributary = np.r_[dy[0] * 0.5, (dy[:-1] + dy[1:]) * 0.5, dy[-1] * 0.5]
    continuity_ratio = np.linalg.norm(interface_imbalance, axis=1) / (TARGET_N * tributary[1:-1])

    isolated_means = []
    coupled_means = []
    relative_differences = []
    for region, nodes in enumerate(isolated_node_maps):
        coupled_mean = float(np.mean(coupled.max_tensile_resultant_n_per_m[nodes]))
        isolated_mean = float(np.mean(isolated[region].max_tensile_resultant_n_per_m))
        coupled_means.append(coupled_mean)
        isolated_means.append(isolated_mean)
        relative_differences.append(abs(coupled_mean - isolated_mean) / TARGET_N)

    force_scale = float(np.sum(np.linalg.norm(outer_force, axis=1)))
    torque_scale = force_scale * RADIUS
    report = {
        "benchmark": "mildly curved two-plate membrane cap with a shared, displacement-continuous interface",
        "status": "manufactured verification only; not a geologic validation or production coupling model",
        "mesh": {"nodes": len(points), "triangles": len(triangles), "triangles_per_plate": [int((regions == i).sum()) for i in (0, 1)], "interface_nodes": int(interface.sum()), "angular_half_width_deg": 0.10 * 180 / np.pi, "angular_half_height_deg": 0.06 * 180 / np.pi},
        "inputs": {"analytic_isotropic_membrane_resultant_n_per_m": TARGET_N, "thickness_m": thickness, "young_modulus_pa": modulus, "poisson_ratio": 0.25, "outer_edges_loaded": True, "interface_external_load_in_joined_solve": 0, "isolated_interface_tractions_equal_and_opposite": True},
        "balance": {"outer_boundary_force_resultant_n": boundary_force.tolist(), "analytic_shell_normal_pressure_reaction_n": pressure_reaction.sum(axis=0).tolist(), "global_force_balance_over_total_boundary_load": float(np.linalg.norm(force_balance) / force_scale), "in_plane_equilibrium_residual_over_load": coupled.equilibrium_residual_n / max(float(np.linalg.norm(coupled.equilibrated_force_n)), 1.0), "outer_boundary_torque_nm": torque.tolist(), "torque_over_load_radius": float(np.linalg.norm(torque) / torque_scale), "isolated_equilibrium_residual_over_load": [s.equilibrium_residual_n / max(float(np.linalg.norm(s.equilibrated_force_n)), 1.0) for s in isolated]},
        "interface": {"shared_displacement_dofs": True, "max_traction_reaction_imbalance_over_target_resultant": float(continuity_ratio.max()), "p95_traction_reaction_imbalance_over_target_resultant": float(np.quantile(continuity_ratio, 0.95)), "interface_nodes": len(interface_indices), "isolated_cut_face_equal_opposite_traction_error": isolated_traction_imbalance},
        "resultant": {"isolated_plate_mean_n_per_m": isolated_means, "joined_plate_region_mean_n_per_m": coupled_means, "relative_difference_from_isolated": relative_differences, "maximum_relative_difference": max(relative_differences), "target_n_per_m": TARGET_N},
        "solver_stop_codes": {"joined": coupled.solver_stop_code, "isolated": [s.solver_stop_code for s in isolated]},
        "interpretation": "A successful result verifies force transfer through this shared-node membrane discretization for the manufactured case. It does not establish that the tectonic force fields or constitutive assumptions are geologically correct, and it does not model faults, slip, damage, or evolving plate boundaries.",
    }
    if report["balance"]["global_force_balance_over_total_boundary_load"] >= 1e-3:
        raise AssertionError("analytic global force balance exceeds 0.1%")
    if report["balance"]["in_plane_equilibrium_residual_over_load"] >= 1e-5:
        raise AssertionError("joined in-plane equilibrium residual/load exceeds 1e-5")
    if report["balance"]["torque_over_load_radius"] >= 1e-8:
        raise AssertionError("outer boundary torque/load-radius exceeds 1e-8")
    if report["interface"]["max_traction_reaction_imbalance_over_target_resultant"] >= 1e-5:
        raise AssertionError("interface reaction mismatch/resultant exceeds 1e-5")
    if report["resultant"]["maximum_relative_difference"] >= 1.0:
        raise AssertionError("coupled and isolated recovered resultants differ by 3% or more")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run_benchmark()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

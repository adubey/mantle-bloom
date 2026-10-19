from collections import Counter
from types import SimpleNamespace

import numpy as np

from app import lithosphere, thin_shell


def _octahedral_shell():
    points = np.array(
        [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
        dtype=float,
    )
    force = np.array(
        [[0, 1e12, 0], [0, -1e12, 0], [1e12, 0, 0], [-1e12, 0, 0], [0, 0, 1e12], [0, 0, -1e12]],
        dtype=float,
    )
    return points, force


def _shallow_cap_uniaxial_case(grid_size, membrane_resultant_n_per_m):
    from app.thin_shell import _spherical_surface_triangles

    x, y = np.meshgrid(np.linspace(-0.05, 0.05, grid_size), np.linspace(-0.05, 0.05, grid_size))
    points = np.column_stack((x.ravel(), y.ravel(), np.sqrt(1.0 - x.ravel() ** 2 - y.ravel() ** 2)))
    triangles = _spherical_surface_triangles(points)
    edges = Counter(
        tuple(sorted(edge))
        for triangle in triangles
        for edge in ((triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0]))
    )
    boundary_edges = [edge for edge, count in edges.items() if count == 1]
    xy_m = points[:, :2] * 6.371e6
    center_m = xy_m.mean(axis=0)
    force_n = np.zeros_like(points)
    for first, second in boundary_edges:
        edge = xy_m[second] - xy_m[first]
        length_m = np.linalg.norm(edge)
        outward = np.array([edge[1], -edge[0]]) / length_m
        midpoint = (xy_m[first] + xy_m[second]) * 0.5
        if np.dot(outward, midpoint - center_m) < 0.0:
            outward = -outward
        nodal_force = membrane_resultant_n_per_m * outward[0] * np.array([1.0, 0.0, 0.0]) * length_m * 0.5
        force_n[first] += nodal_force
        force_n[second] += nodal_force
    return points, triangles, force_n


def _irregular_cap_case(membrane_resultant_n_per_m):
    polygon = np.array(
        [[-0.05, -0.025], [0.005, -0.04], [0.045, -0.02], [0.052, 0.012], [0.025, 0.045], [-0.02, 0.04], [-0.048, 0.015]]
    )
    x_grid, y_grid = np.meshgrid(np.linspace(-0.05, 0.05, 9), np.linspace(-0.04, 0.045, 9))
    coordinates = []
    for x, y in zip(x_grid.ravel(), y_grid.ravel()):
        inside = False
        for i in range(len(polygon)):
            x1, y1 = polygon[i]
            x2, y2 = polygon[(i + 1) % len(polygon)]
            if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
                inside = not inside
        if inside:
            coordinates.append((x, y))
    for i in range(len(polygon)):
        for fraction in (1 / 3, 2 / 3):
            coordinates.append(polygon[i] * (1 - fraction) + polygon[(i + 1) % len(polygon)] * fraction)
    coordinates.extend(polygon.tolist())
    xy = np.unique(np.round(np.asarray(coordinates), 8), axis=0)
    points = np.column_stack((xy, np.sqrt(1.0 - np.sum(xy**2, axis=1))))
    from app.thin_shell import _spherical_surface_triangles

    triangles = _spherical_surface_triangles(points)
    edge_counts = Counter(
        tuple(sorted(edge))
        for triangle in triangles
        for edge in ((triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0]))
    )
    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    xy_m = points[:, :2] * lithosphere.PLANET_RADIUS_M
    center_m = xy_m.mean(axis=0)
    force = np.zeros_like(points)
    clamped = np.zeros(len(points), dtype=bool)
    for first, second in boundary_edges:
        edge = xy_m[second] - xy_m[first]
        length = np.linalg.norm(edge)
        outward = np.array([edge[1], -edge[0]]) / length
        midpoint = (xy_m[first] + xy_m[second]) * 0.5
        if np.dot(outward, midpoint - center_m) < 0.0:
            outward = -outward
        nodal_force = membrane_resultant_n_per_m * outward[0] * np.array([1.0, 0.0, 0.0]) * length * 0.5
        force[first] += nodal_force
        force[second] += nodal_force
        if midpoint[0] < -0.04 * lithosphere.PLANET_RADIUS_M:
            clamped[first] = True
            clamped[second] = True
    return points, triangles, force, clamped


def test_thin_shell_transfers_tangential_load_to_stress_and_strength_ratio():
    points, force = _octahedral_shell()
    result = thin_shell.solve_force_field(
        SimpleNamespace(points_xyz=points, areas_m2=np.full(len(points), 1e12), total_n=force),
        100_000.0,
        70e9,
        np.full(len(points), 1e12),
    )

    assert result.solver_stop_code == 1
    assert result.equilibrium_residual_n < 1e-9 * np.linalg.norm(result.equilibrated_force_n)
    residual_torque = np.cross(result.equilibrated_force_n, points * 6.371e6).sum(axis=0)
    assert np.linalg.norm(residual_torque) < 1e-9 * np.linalg.norm(np.cross(force, points * 6.371e6).sum(axis=0)) + 1e-5
    assert np.max(result.max_tensile_resultant_n_per_m) > 0.0
    assert np.allclose(
        result.failure_ratio,
        result.max_tensile_resultant_n_per_m / result.tensile_strength_resultant_n_per_m,
    )
    assert np.allclose(np.sum(result.displacement_m * points, axis=1), 0.0, atol=1e-9)
    doubled = thin_shell.solve_thin_shell(
        points,
        np.full(len(points), 1e12),
        2.0 * force,
        thickness_m=100_000.0,
        young_modulus_pa=70e9,
        tensile_strength_resultant_n_per_m=np.full(len(points), 1e12),
    )
    assert np.allclose(doubled.stress_resultant_n_per_m, 2.0 * result.stress_resultant_n_per_m, rtol=1e-7, atol=1e-8)
    assert np.allclose(doubled.failure_ratio, 2.0 * result.failure_ratio, rtol=1e-7, atol=1e-12)


def test_thin_shell_removes_plate_scale_rigid_rotation_load():
    points, _ = _octahedral_shell()
    force = np.cross(np.broadcast_to([2e-4, -1e-4, 3e-4], points.shape), points * 6.371e6)
    result = thin_shell.solve_thin_shell(
        points,
        np.full(len(points), 1e12),
        force,
        thickness_m=100_000.0,
        young_modulus_pa=70e9,
        tensile_strength_resultant_n_per_m=np.full(len(points), 1e12),
    )

    assert np.linalg.norm(result.equilibrated_force_n) < 1e-12 * np.linalg.norm(force)
    assert np.max(result.max_tensile_resultant_n_per_m) < 1e-10
    assert np.allclose(result.rigid_rotation_force_n, force)


def test_thin_shell_recovers_uniaxial_resultant_across_shallow_cap_meshes():
    expected_resultant = 1e10
    measured = []
    for grid_size, expected_triangle_count in ((3, 8), (5, 32), (9, 128)):
        points, triangles, force = _shallow_cap_uniaxial_case(grid_size, expected_resultant)
        assert len(triangles) == expected_triangle_count
        result = thin_shell.solve_thin_shell(
            points,
            np.full(len(points), 1e10),
            force,
            thickness_m=100_000.0,
            young_modulus_pa=70e9,
            tensile_strength_resultant_n_per_m=np.full(len(points), 2.0 * expected_resultant),
        )
        mean_resultant = np.mean(result.max_tensile_resultant_n_per_m)
        measured.append(mean_resultant)
        assert np.isclose(mean_resultant, expected_resultant, rtol=2e-3)
        assert result.equilibrium_residual_n < 1e-6 * np.linalg.norm(result.equilibrated_force_n)
        assert np.isclose(np.mean(result.failure_ratio), 0.5 * mean_resultant / expected_resultant, rtol=1e-6)
    assert (max(measured) - min(measured)) / expected_resultant < 1e-4


def test_thin_shell_handles_irregular_convex_outline_and_clamped_edge():
    expected_resultant = 1e10
    points, triangles, force, clamped_nodes = _irregular_cap_case(expected_resultant)
    assert len(points) == 67
    assert len(triangles) == 111
    assert clamped_nodes.sum() == 6
    common = dict(
        points_xyz=points,
        areas_m2=np.full(len(points), 1e10),
        applied_force_n=force,
        thickness_m=100_000.0,
        young_modulus_pa=70e9,
        tensile_strength_resultant_n_per_m=np.full(len(points), 2.0 * expected_resultant),
    )
    free = thin_shell.solve_thin_shell(**common)
    clamped = thin_shell.solve_thin_shell(**common, fixed_node_mask=clamped_nodes)
    free_mean = np.mean(free.max_tensile_resultant_n_per_m)
    clamped_mean = np.mean(clamped.max_tensile_resultant_n_per_m)

    assert np.isclose(free_mean, expected_resultant, rtol=2e-3)
    assert abs(clamped_mean - free_mean) > 0.02 * expected_resultant
    assert np.all(clamped.displacement_m[clamped_nodes] == 0.0)
    assert free.equilibrium_residual_n < 1e-5 * np.linalg.norm(free.equilibrated_force_n)
    assert clamped.equilibrium_residual_n < 1e-5 * np.linalg.norm(clamped.equilibrated_force_n)

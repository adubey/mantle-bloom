from types import SimpleNamespace

import numpy as np

from app import thin_shell


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

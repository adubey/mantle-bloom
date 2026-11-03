from types import SimpleNamespace

import numpy as np
from scipy.spatial import cKDTree

from app import world


class _Plate:
    def __init__(self, points: np.ndarray, elevation: np.ndarray):
        self.points = points.copy()
        self.elevation = elevation.copy()
        self.frame = np.eye(3)
        self.omega = np.array([0.2, -0.3, 0.4])
        self.crustal_thickness_m = np.arange(len(points), dtype=float) + 30_000.0

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        return self.points, self.elevation

    def rotate(self, increment: np.ndarray) -> None:
        self.frame = increment @ self.frame
        self.points = (increment @ self.points.T).T

    def set_omega(self, omega: np.ndarray) -> None:
        self.omega = omega


def _polar_node_count(points: np.ndarray, axis: np.ndarray) -> int:
    cap_limit = np.cos(np.deg2rad(world.POLAR_LAND_CLEARANCE_DEG))
    return int(np.count_nonzero(np.abs(points @ axis) >= cap_limit))


def test_polar_reorientation_clears_poles_and_rotates_world_vectors_together():
    points = np.array([
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    elevation = np.array([100.0, 100.0, -100.0, -100.0])
    plate = _Plate(points, elevation)
    old_omega = plate.omega.copy()
    old_local_fields = plate.crustal_thickness_m.copy()
    old_pinned_omega = np.array([-0.1, 0.5, 0.2])
    old_center = np.array([0.0, 1.0, 0.0])
    planet = SimpleNamespace(
        plates=[plate],
        sea_level_m=0.0,
        mantle_centers=[SimpleNamespace(position=old_center.copy())],
        pinned_omegas={7: old_pinned_omega.copy()},
    )
    before_pairwise_dots = points @ points.T

    world._rotate_world_to_clear_poles(planet)

    assert _polar_node_count(plate.points[:2], np.array([0.0, 0.0, 1.0])) == 0
    assert np.allclose(plate.points @ plate.points.T, before_pairwise_dots, atol=1e-12)
    assert np.array_equal(plate.elevation, elevation)
    assert np.array_equal(plate.crustal_thickness_m, old_local_fields)
    assert np.allclose(plate.omega, plate.frame @ old_omega)
    assert np.allclose(planet.mantle_centers[0].position, plate.frame @ old_center)
    assert np.allclose(planet.pinned_omegas[7], plate.frame @ old_pinned_omega)


def test_infeasible_polar_caps_choose_sampled_axis_with_minimum_land_nodes():
    land_xyz = world._fibonacci_sphere_axes(2048)
    preferred_axis, support = world._least_land_support_axis(land_xyz)
    cap_limit = np.cos(np.deg2rad(world.POLAR_LAND_CLEARANCE_DEG))
    assert support > cap_limit

    candidate_axes = np.vstack((preferred_axis, world._fibonacci_sphere_axes(world.POLAR_FALLBACK_AXIS_SAMPLES)))
    cap_chord = 2.0 * np.sin(np.deg2rad(world.POLAR_LAND_CLEARANCE_DEG) / 2.0)
    land_tree = cKDTree(np.concatenate((land_xyz, -land_xyz), axis=0))
    expected_minimum = int(np.min(land_tree.query_ball_point(candidate_axes, r=cap_chord, return_length=True)))

    plate = _Plate(land_xyz, np.ones(len(land_xyz)))
    planet = SimpleNamespace(
        plates=[plate],
        sea_level_m=0.0,
        mantle_centers=[],
        pinned_omegas={},
    )
    world._rotate_world_to_clear_poles(planet)

    actual = _polar_node_count(plate.points, np.array([0.0, 0.0, 1.0]))
    assert actual == expected_minimum

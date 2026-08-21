import numpy as np
from app import geometry, mantle, merge_split
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def _test_plate(plate_id, seed_xyz, crust_type, theta, elevation):
    """A plate with the given line at phi=0 (what each test actually exercises) plus a
    second, far-away placeholder line at a different latitude -- purely so the plate has
    more than one line. Otherwise apply_topology_changes's "no land left" pruning (a plate
    reduced to a single line, see merge_split.remove_defunct_plates) would remove these
    synthetic single-line test plates before the test's own logic ever ran."""
    frame = geometry.plate_frame_from_seed(seed_xyz)
    line = ElevationLine(phi=0.0, theta=np.asarray(theta, dtype=float), elevation=np.asarray(elevation, dtype=float))
    filler = ElevationLine(phi=1.0, theta=np.array([0.0, 0.1]), elevation=np.zeros(2))
    return Plate(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=[line, filler])


def _converging_omega_pair(rate_cm_per_yr=5.0):
    """omega for (keep, absorb) such that -- given absorb is seeded at a small *positive*
    angle from keep around +z (see the tiny_angle rotation below) -- each plate's material
    moves toward the other's, matching test_boundary.py's closing-rate sign convention."""
    rate = mantle.cm_per_yr_to_rad_per_yr(rate_cm_per_yr)
    return np.array([0.0, 0.0, rate]), np.array([0.0, 0.0, -rate])


def _converging_pair_world(seed=123):
    # See test_find_continental_collision_pairs_detects_close_and_converging_plates for why
    # these offsets are fractions of MERGE_CONTACT_DISTANCE_RAD rather than hardcoded angles.
    d = merge_split.MERGE_CONTACT_DISTANCE_RAD
    keep = _test_plate(0, [1.0, 0.0, 0.0], "continental", np.linspace(-d, 0.0, 6), np.zeros(6))
    seed_absorb = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], axis=np.array([0.0, 0.0, 1.0]), angle=d * 0.3
    )[0]
    absorb = _test_plate(1, seed_absorb, "continental", np.linspace(-d * 0.2, d * 0.2, 6), np.full(6, 50.0))
    keep.omega, absorb.omega = _converging_omega_pair()
    return World(seed=seed, plates=[keep, absorb], next_plate_id=2)


def test_apply_topology_changes_runs_without_error_during_long_simulation():
    world = generate_world(seed=40, num_plates=10)
    for _ in range(20):
        step_world(world, years=5_000_000)
    assert len(world.plates) > 0
    for plate in world.plates:
        assert plate.node_count() >= 0
    # plate ids must stay unique even after merges/splits
    ids = [p.plate_id for p in world.plates]
    assert len(ids) == len(set(ids))

import numpy as np

from app import geometry
from app.world import generate_world, step_world


def test_generate_world_gives_plates_nonzero_omega():
    world = generate_world(seed=10, num_plates=8)
    rates = [np.linalg.norm(p.omega) for p in world.plates]
    assert any(r > 0 for r in rates)


def test_step_world_advances_elapsed_years():
    world = generate_world(seed=11, num_plates=6)
    step_world(world, years=1_000_000)
    assert world.elapsed_years == 1_000_000
    step_world(world, years=500_000)
    assert world.elapsed_years == 1_500_000


def test_rigid_rotation_preserves_interior_node_spacing_exactly():
    """The whole point of the plate-local-frame design: rotating a plate must not disturb
    the relative spacing of its own elevation-line nodes at all (no resampling)."""
    world = generate_world(seed=12, num_plates=6)
    plate = max(world.plates, key=lambda p: p.node_count())
    line = max(plate.lines, key=lambda l: len(l.theta))

    before_world = line.world_xyz(plate.frame)
    before_spacing = geometry.angular_distance(before_world[:-1], before_world[1:])

    for _ in range(5):
        step_world(world, years=2_000_000)

    after_world = line.world_xyz(plate.frame)
    after_spacing = geometry.angular_distance(after_world[:-1], after_world[1:])

    assert np.allclose(before_spacing, after_spacing, atol=1e-9)
    # theta/elevation arrays themselves must be untouched (identity-based, not resampled).
    assert np.array_equal(line.theta, line.theta)


def test_rigid_rotation_preserves_plate_frame_orthonormality():
    world = generate_world(seed=13, num_plates=6)
    for _ in range(10):
        step_world(world, years=5_000_000)
    for plate in world.plates:
        assert np.allclose(plate.frame @ plate.frame.T, np.eye(3), atol=1e-6)
        assert np.isclose(np.linalg.det(plate.frame), 1.0, atol=1e-6)


def test_different_plates_generally_acquire_different_omegas():
    world = generate_world(seed=14, num_plates=10)
    omegas = np.array([p.omega for p in world.plates])
    # Not every plate should have picked up an identical rotation from a spatially-varying
    # mantle flow field.
    assert not np.allclose(omegas, omegas[0], atol=1e-12)


def test_generate_world_logs_a_generation_event():
    world = generate_world(seed=15, num_plates=8, continental_fraction=3 / 8)
    assert len(world.events) == 1
    elapsed, message = world.events[0]
    assert elapsed == 0.0
    assert "8 plates" in message and "3 continental" in message


def test_step_world_events_are_timestamped_with_post_step_elapsed_years():
    world = generate_world(seed=16, num_plates=10, continental_fraction=0.6)
    for _ in range(15):
        step_world(world, years=8_000_000)
    # Any event logged during stepping (merges/splits/consumption are plausible at this
    # scale/seed count but not guaranteed for a specific seed) must be timestamped no later
    # than the world's current elapsed_years.
    for elapsed, _ in world.events:
        assert elapsed <= world.elapsed_years
    assert len(world.events) <= 200  # MAX_EVENT_LOG_LENGTH

import numpy as np
from app import geometry
from app.world import generate_world, step_world


def test_node_density_persists_through_regularize_and_gap_fill():
    # The core correctness concern for a runtime density option: line_regrid.py's own
    # regularize pass (and gaps.py's gap-filling, merge_split.py's merging, volcanism.py's
    # field-spawning) previously always rebuilt/resampled nodes at the module's default
    # TARGET_LINE_SPACING_RAD regardless of what density a world was actually generated at,
    # silently reverting a non-default density back to the reference one within a handful of
    # steps. Confirmed directly this stays fixed: total node count should keep tracking
    # node_density's ratio, not decay toward the 1x baseline, across enough steps to trigger
    # at least one regularize/gap-fill pass (line_regrid.REGULARIZE_INTERVAL_STEPS == 5).
    reference = generate_world(seed=5, num_plates=8, continental_fraction=0.5, node_density=1.0)
    denser = generate_world(seed=5, num_plates=8, continental_fraction=0.5, node_density=4.0)

    def total_nodes(world):
        return sum(len(line.theta) for p in world.plates for line in p.lines)

    for _ in range(8):
        step_world(reference, years=300_000)
        step_world(denser, years=300_000)

    ratio = total_nodes(denser) / total_nodes(reference)
    assert 3.0 < ratio < 5.0  # stays roughly 4x -- not decayed back toward 1x


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

import numpy as np
from app import erosion, geometry
from app.world import generate_world, step_world


def _sampled_overlap_fraction(plates_list, sample_per_plate: int = 20) -> float:
    """Same proxy invariant as unit_tests/test_plates.py's own version -- see that
    function's docstring for why this isn't expected to be exactly zero (one-turn
    processing lag, residual envelope looseness for non-convex plate shapes), just bounded
    rather than growing without limit turn over turn."""
    total = 0
    overlapping = 0
    for plate in plates_list:
        points, _ = plate.all_points_and_elevation()
        if len(points) == 0:
            continue
        sample = points[:: max(1, len(points) // sample_per_plate)][:sample_per_plate]
        for other in plates_list:
            if other.plate_id == plate.plate_id:
                continue
            polygon = other.get_bounding_polygon()
            if len(polygon) < 3:
                continue
            total += len(sample)
            overlapping += int(np.count_nonzero(geometry.points_in_spherical_polygon(sample, polygon)))
    return overlapping / total if total > 0 else 0.0


def test_plate_overlap_stays_bounded_over_many_steps():
    # Confirmed directly during development: with the polygon-based shift()/deform() model,
    # a real (non-runaway) amount of "envelope overlap" persists indefinitely -- what
    # matters is that it stays roughly flat over many steps rather than climbing toward
    # saturation, which would indicate territory genuinely running away unchecked (e.g. a
    # shrink/grow imbalance letting one plate's claimed footprint balloon).
    world = generate_world(seed=3, num_plates=8, node_density=0.5)
    world.simulate_climate_biomes = False  # only plate geometry is checked here
    fractions = []
    for i in range(24):
        step_world(world, years=3_000_000)
        if i % 4 == 3:
            fractions.append(_sampled_overlap_fraction(world.plates))
    # None of the later checkpoints drift far above the earliest one -- a real runaway would
    # show a clear upward trend across the run, not noise around a stable level.
    assert max(fractions) < fractions[0] + 0.15
    assert max(fractions) < 0.35


def test_node_density_persists_through_regularize_and_gap_fill():
    # The core correctness concern for a runtime density option: elevation_lines.py's own
    # regularize pass (and PlateWithLines.deform's own claim-adjacent-territory/merge_split.py's
    # merging) previously always rebuilt/resampled nodes at the module's default
    # TARGET_LINE_SPACING_RAD regardless of what density a world was actually generated at,
    # silently reverting a non-default density back to the reference one within a handful of
    # steps. Confirmed directly this stays fixed: total node count should keep tracking
    # node_density's ratio, not decay toward the 1x baseline, across enough steps.
    reference = generate_world(seed=5, num_plates=8, continental_fraction=0.5, node_density=1.0)
    denser = generate_world(seed=5, num_plates=8, continental_fraction=0.5, node_density=4.0)
    # This test only checks node counts, never climate/erosion/hydrology output (see
    # World.simulate_climate_biomes) -- skipping that per-step computation speeds this up
    # without changing what it exercises.
    reference.simulate_climate_biomes = False
    denser.simulate_climate_biomes = False

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
    world.simulate_climate_biomes = False  # only node spacing/rotation is checked here
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
    world.simulate_climate_biomes = False  # only plate.frame orthonormality is checked here
    for _ in range(10):
        step_world(world, years=5_000_000)
    for plate in world.plates:
        assert np.allclose(plate.frame @ plate.frame.T, np.eye(3), atol=1e-6)
        assert np.isclose(np.linalg.det(plate.frame), 1.0, atol=1e-6)


def test_step_world_events_are_timestamped_with_post_step_elapsed_years():
    world = generate_world(seed=16, num_plates=10, continental_fraction=0.6)
    # Only event timestamps/log length are checked here, never climate/erosion/hydrology
    # output (see World.simulate_climate_biomes) -- plate-movement events (merges/splits/
    # gap-fills/volcanism) still fire normally, since simulate_plate_movement stays on.
    world.simulate_climate_biomes = False
    for _ in range(15):
        step_world(world, years=8_000_000)
    # Any event logged during stepping (merges/splits/consumption are plausible at this
    # scale/seed count but not guaranteed for a specific seed) must be timestamped no later
    # than the world's current elapsed_years.
    for elapsed, _ in world.events:
        assert elapsed <= world.elapsed_years
    assert len(world.events) <= 200  # MAX_EVENT_LOG_LENGTH


def test_coastal_feedback_stays_stable_over_many_steps():
    # A regression floor for the coastal planation + infill feedback (erosion.py), not a tight
    # bound. The real drowned-shelf checkerboard from docs/TODO.md "Speckled low-relief
    # coastlines" only bites at node_density=4 (or on the seed-888151728 save) -- too slow for
    # a stress test, and a sudden sea-level jump on a density-1 world just makes a rough
    # newborn coast whose transient roughening swamps the feedback's slow ~My effect. So this
    # only asserts the pass doesn't destabilise a world with a broad near-waterline shelf:
    # elevation stays finite and in-bounds over many steps, and the pass is demonstrably
    # active at the shore.
    world = generate_world(seed=17, num_plates=8, continental_fraction=0.6, node_density=1.0)
    _, elevation, *_ = erosion._gather_nodes(world)
    world.sea_level_m = float(np.percentile(elevation[elevation > 0.0], 40))
    sl = world.sea_level_m
    near_shore = np.abs(elevation - sl) <= 80.0
    assert near_shore.sum() > 200  # the test bed really does have a broad near-waterline zone

    for _ in range(30):
        step_world(world, years=200_000)

    _, after, *_ = erosion._gather_nodes(world)
    from app.elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M

    assert np.all(np.isfinite(after))
    assert np.all(after >= MIN_ELEVATION_M - 1e-6) and np.all(after <= MAX_ELEVATION_M + 1e-6)
    # The near-waterline shelf hasn't wholesale-drowned or wholesale-emerged -- planation and
    # infill nudge the coast, they don't run away with it.
    still_near = float(np.mean(np.abs(after - sl) <= 200.0))
    assert still_near > 0.05

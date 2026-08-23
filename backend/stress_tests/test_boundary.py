import numpy as np
from app import boundary, geometry, mantle
from app.plates import TARGET_LINE_SPACING_RAD, ElevationLine, PlateWithLines
from app.world import World, generate_world, step_world


def _plate(plate_id, crust_type, theta, omega, base_elevation):
    """A single-line plate at a fixed seed ([1,0,0], identity frame) so a theta offset
    directly converts to a real distance (theta_rad * PLANET_RADIUS_KM, small-angle exact
    for the scale these tests use) -- lets a test construct nodes at precise, known
    distances from a neighboring plate's own cluster."""
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(len(theta), base_elevation))
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type=crust_type, omega=np.asarray(omega, dtype=float), lines=[line])


def test_stepping_with_boundary_evolution_keeps_lines_sorted_and_elevation_bounded():
    # node_density=0.5 (the coarsest choice, an eighth of the default 4.0, see
    # plates.NODE_DENSITY_CHOICES) -- this test only checks sortedness/elevation bounds, not
    # exact node positions.
    world = generate_world(seed=21, num_plates=8, node_density=0.5)
    # This test only checks boundary.py's own output (line sorting, elevation bounds), which
    # doesn't depend on climate/erosion/hydrology at all (see World.simulate_climate_biomes) --
    # skipping that per-step computation cuts this test's runtime substantially without
    # changing what it exercises.
    world.simulate_climate_biomes = False
    for _ in range(15):
        step_world(world, years=3_000_000)

    for plate in world.plates:
        for line in plate.lines:
            assert np.all(np.diff(line.theta) > 0), "line thetas must stay strictly ascending"
            assert np.all(line.elevation >= boundary.MIN_ELEVATION_M)
            assert np.all(line.elevation <= boundary.MAX_ELEVATION_M)


def test_boundary_evolution_changes_node_counts_over_time():
    world = generate_world(seed=22, num_plates=10, node_density=0.5)  # see the test above for why this is safe
    world.simulate_climate_biomes = False  # see the test above for why this is safe here
    before = sum(p.node_count() for p in world.plates)
    for _ in range(20):
        step_world(world, years=4_000_000)
    after = sum(p.node_count() for p in world.plates)
    assert after != before

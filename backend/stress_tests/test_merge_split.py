import numpy as np
from app import geometry, mantle, merge_split
from app.elevation_lines import line_spacing_rad
from app.plates import ElevationLine, PlateWithLines, node_components
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
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=[line, filler])


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
    keep_omega, absorb_omega = _converging_omega_pair()
    keep.set_omega(keep_omega)
    absorb.set_omega(absorb_omega)
    return World(seed=seed, plates=[keep, absorb], next_plate_id=2)


def test_apply_topology_changes_runs_without_error_during_long_simulation():
    # node_density=0.5 (the coarsest choice, an eighth of the default 4.0, see
    # plates.NODE_DENSITY_CHOICES) -- this test only checks topology invariants (plate count,
    # unique ids), not exact node positions.
    world = generate_world(seed=40, num_plates=10, node_density=0.5)
    # merge_split.py only reacts to plate/node geometry, never to climate/erosion/hydrology
    # output (see World.simulate_climate_biomes) -- skipping that per-step computation cuts
    # this long-running test's cost substantially without changing what it exercises.
    world.simulate_climate_biomes = False
    for _ in range(20):
        step_world(world, years=5_000_000)
    assert len(world.plates) > 0
    for plate in world.plates:
        assert plate.node_count() >= 0
    # plate ids must stay unique even after merges/splits
    ids = [p.plate_id for p in world.plates]
    assert len(ids) == len(set(ids))


def test_long_simulation_leaves_no_plate_severed_into_disconnected_landmasses():
    # Regression guard for the failure the defragmentation pass exists to fix: subduction/
    # transform carving one Plate's node cloud into two disconnected lobes, or stranding a
    # comb of one-node rows, with maybe_split_plate (mantle-flow-only) never noticing. With
    # defragment_plates running on its cadence, every surviving plate should stay a single
    # connected component (barring a fragment younger than one defrag interval).
    world = generate_world(seed=40, num_plates=10, node_density=0.5)
    world.simulate_climate_biomes = False
    connect_radius_rad = merge_split.DEFRAG_CONNECT_RADIUS_MULT * line_spacing_rad(world.node_density)
    min_fragment_nodes = max(1, round(merge_split.DEFRAG_FRAGMENT_MIN_NODES * world.node_density))

    for step in range(40):
        step_world(world, years=4_000_000)
        # Check only on steps where defrag has just run (see DEFRAG_INTERVAL_STEPS), so a
        # plate transiently severed mid-interval isn't counted against us.
        if world.steps_taken % merge_split.DEFRAG_INTERVAL_STEPS != 0:
            continue
        for plate in world.plates:
            points, _ = plate.all_points_and_elevation()
            if len(points) < 2:
                continue
            labels = node_components(points, connect_radius_rad)
            _, counts = np.unique(labels, return_counts=True)
            big = counts[counts >= min_fragment_nodes]
            assert len(big) <= 1, (
                f"plate {plate.plate_id} has {len(big)} components >= {min_fragment_nodes} "
                f"nodes after step {world.steps_taken}: {sorted(counts.tolist(), reverse=True)[:6]}"
            )

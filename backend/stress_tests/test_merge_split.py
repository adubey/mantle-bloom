import numpy as np
from app import merge_split
from app.elevation_lines import line_spacing_rad
from app.plates import node_components
from app.world import generate_world, step_world


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

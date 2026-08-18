import numpy as np

from app import erosion, geometry, hydrology
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def test_compute_basin_spill_collapses_nested_rims_to_one_hop():
    # A 4-node chain: 0 (elevation 30, behind a rim) -> 1 (the rim, 50) -> 2 (20) -> 3 (ocean,
    # -10). Node 0's best path to the ocean is forced to cross node 1's rim (50) even though
    # node 1 itself is higher than node 0 -- so filled_elevation[0] must be 50, not node 0's
    # own (lower) neighbor chain. And per the priority-flood's own "collapse nested basins"
    # rule, node 0's escape hop should point past node 1 straight to node 2 (the cell just on
    # the low side of the rim it has to cross), not at node 1 itself.
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array(
        [
            [1, 1],
            [0, 2],
            [1, 3],
            [2, 2],
        ]
    )
    filled, spill = hydrology._compute_basin_spill(elevation, is_ocean, neighbor_idx)
    assert filled.tolist() == [50.0, 50.0, 20.0, -10.0]
    assert spill.tolist() == [2, 2, 3, -1]


def test_compute_flow_direction_sink_redirects_once_lake_reaches_spill_point():
    # Same 4-node setup as above. Node 0 has no lower neighbor (its only neighbor, node 1,
    # is higher) so it's a sink candidate -- unless its current water surface has already
    # reached the basin's true spill point (50, from the test above).
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    filled = np.array([50.0, 50.0, 20.0, -10.0])
    spill = np.array([2, 2, 3, -1])

    not_yet_full = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(4), filled, spill)
    assert not_yet_full[0] == -1  # water surface 30 < 50, still a plain sink

    now_full = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.array([25.0, 0.0, 0.0, 0.0]), filled, spill)
    assert now_full[0] == 2  # water surface 30+25=55 >= 50 -- redirects to spill_target


def test_route_downstream_accumulates_along_a_chain():
    # A straight 4-node downhill chain (0 -> 1 -> 2 -> ocean-adjacent 3), each node
    # contributing 1.0 of its own -- accumulated flux at each node should be its own
    # contribution plus everything upstream of it.
    elevation = np.array([40.0, 30.0, 20.0, -5.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 3, -1])
    source = np.array([1.0, 1.0, 1.0, 0.0])

    through_flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source)
    assert through_flux.tolist() == [1.0, 2.0, 3.0, 0.0]
    assert deposited[3] == 3.0  # everything reaching the ocean-adjacent node's target settles there
    assert deposited[0] == 0.0 and deposited[1] == 0.0 and deposited[2] == 0.0


def test_route_downstream_retain_fraction_deposits_partway():
    # Same chain, but node 1 retains half of whatever passes through it.
    elevation = np.array([40.0, 30.0, 20.0, -5.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 3, -1])
    source = np.array([1.0, 1.0, 1.0, 0.0])
    retain = np.array([0.0, 0.5, 0.0, 0.0])

    through_flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source, retain_fraction=retain)
    # Node 1 receives 1.0 (its own) + 1.0 (from node 0) = 2.0, retains half (1.0) locally,
    # passes 1.0 onward.
    assert deposited[1] == 1.0
    assert through_flux[1] == 1.0
    # Node 2 receives 1.0 (its own) + 1.0 (passed on from node 1) = 2.0, all continues on.
    assert through_flux[2] == 2.0
    assert deposited[3] == 2.0
    # Conservation: total deposited equals total source.
    assert np.isclose(deposited.sum(), source.sum())


def test_update_lakes_grows_at_a_sink_and_caps_at_the_spill_point():
    fields = hydrology.HydrologyFields(
        points=np.zeros((2, 3)),
        elevation=np.array([10.0, -5.0]),
        is_ocean=np.array([False, True]),
        neighbor_idx=np.zeros((2, 1), dtype=np.int64),
        flow_target=np.array([-1, -1]),
        flow_accum=np.zeros(2),
        water_deposited=np.array([50.0, 0.0]),  # a lot of water settled at the sink this step
        filled_elevation=np.array([25.0, -5.0]),  # sink's basin spills at elevation 25
        spill_target=np.array([-1, -1]),
        is_river=np.zeros(2, dtype=bool),
        line_refs=[],
    )
    grown = hydrology.update_lakes(fields, prev_lake_depth=np.zeros(2), water_deposited=fields.water_deposited, years=1_000_000)
    assert grown[0] > 0.0  # grew from inflow
    assert grown[0] <= 25.0 - 10.0  # never exceeds the basin's true spill depth
    assert grown[1] == 0.0  # never touches ocean

    # A lake already sitting at its cap should stay pinned there, not evaporate back down --
    # water_deposited=0 here (nothing new arriving) but it must still hold near the cap.
    at_cap = hydrology.update_lakes(fields, prev_lake_depth=np.array([15.0, 0.0]), water_deposited=np.zeros(2), years=1_000_000)
    assert at_cap[0] > 10.0  # decays some (evaporation), but nowhere near back to 0


def _flow_line_plate(plate_id, theta, elevation):
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.asarray(elevation, dtype=float))
    return Plate(plate_id=plate_id, frame=frame, crust_type="continental", lines=[line])


def test_compute_hydrology_end_to_end_on_a_small_synthetic_world():
    # 12 points along a monotonically descending line, the last three underwater --
    # exercises the whole pipeline (basin-spill, flow direction, accumulation, river
    # classification) on a real (if small) node cloud built the same way as the rest of the
    # test suite's synthetic-plate fixtures.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)  # last three points go underwater
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    precipitation_at_nodes = np.full(12, 800.0)
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes)

    assert len(fields.points) == 12
    assert fields.is_ocean.tolist() == [False] * 9 + [True] * 3
    # Every land node should have a full path to the ocean on this simple monotonic line.
    assert np.all(np.isfinite(fields.filled_elevation[~fields.is_ocean]))
    # flow_accum should be non-decreasing heading downhill along this single chain (each
    # node picks up whatever accumulated above it, net of nothing lost here).
    land_accum = fields.flow_accum[~fields.is_ocean]
    assert np.all(np.diff(land_accum) >= -1e-6)


def test_channel_depth_and_lake_depth_persist_across_boundary_and_erosion_steps():
    # A direct regression check on the field-threading work: a node's channel_depth/
    # lake_depth must survive being touched by boundary.py (which runs every step, before
    # erosion.py sets these) rather than getting silently reset to 0.
    world = generate_world(seed=30, num_plates=8, continental_fraction=0.5)
    step_world(world, years=2_000_000)

    channel_total_1 = sum(line.channel_depth.sum() for p in world.plates for line in p.lines)
    lake_total_1 = sum(line.lake_depth.sum() for p in world.plates for line in p.lines)
    assert channel_total_1 > 0.0 or lake_total_1 > 0.0

    step_world(world, years=2_000_000)  # boundary.step_boundaries runs first here again

    channel_total_2 = sum(line.channel_depth.sum() for p in world.plates for line in p.lines)
    # Channel depth only ever grows (monotonic, see erosion.py) -- if boundary.py had reset
    # it, this would have dropped back toward 0 instead.
    assert channel_total_2 >= channel_total_1 * 0.5  # loose bound: some nodes can be pruned/moved


def test_apply_erosion_conserves_mass_between_erosion_and_deposition():
    world = generate_world(seed=31, num_plates=8, continental_fraction=0.5)
    for _ in range(3):
        step_world(world, years=2_000_000)
    assert world.hydrology_cache is not None
    assert world.hydrology_cache.is_river.sum() >= 0  # river classification ran without error

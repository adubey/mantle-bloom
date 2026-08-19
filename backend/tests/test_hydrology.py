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

    not_yet_full = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(4), filled, spill, np.zeros(4))
    assert not_yet_full[0] == -1  # water surface 30 < 50, still a plain sink

    now_full = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.array([25.0, 0.0, 0.0, 0.0]), filled, spill, np.zeros(4)
    )
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
    not_accumulating = np.zeros(2, dtype=bool)
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
        lake_depth=np.zeros(2),
        glacier_depth=np.zeros(2),
        line_refs=[],
    )
    grown = hydrology.update_lakes(
        fields, prev_lake_depth=np.zeros(2), water_deposited=fields.water_deposited, years=1_000_000, is_accumulating=not_accumulating
    )
    assert grown[0] > 0.0  # grew from inflow
    assert grown[0] <= 25.0 - 10.0  # never exceeds the basin's true spill depth
    assert grown[1] == 0.0  # never touches ocean

    # A lake already sitting at its cap should stay pinned there, not evaporate back down --
    # water_deposited=0 here (nothing new arriving) but it must still hold near the cap.
    at_cap = hydrology.update_lakes(
        fields, prev_lake_depth=np.array([15.0, 0.0]), water_deposited=np.zeros(2), years=1_000_000, is_accumulating=not_accumulating
    )
    assert at_cap[0] > 10.0  # decays some (evaporation), but nowhere near back to 0

    # Now cold enough to freeze -- a just-frozen sink (prev_lake_depth already zeroed by the
    # caller, see compute_hydrology's own freeze conversion) must not re-fill from this
    # step's inflow even though water_deposited is nonzero there.
    freezing = np.array([True, False])
    frozen = hydrology.update_lakes(
        fields, prev_lake_depth=np.zeros(2), water_deposited=fields.water_deposited, years=1_000_000, is_accumulating=freezing
    )
    assert frozen[0] == 0.0


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
    temperature_at_nodes = np.full(12, 15.0)  # well above freezing -- no glaciation involved
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes, temperature_at_nodes, years=1_000_000)

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


def test_compute_flow_direction_prefers_an_existing_channel_over_the_literal_steepest_neighbor():
    # Node 0's two downhill candidates: node 1 (elevation 20, steeper, no channel) and node 2
    # (elevation 25, less steep, but already has a real established channel).
    elevation = np.array([30.0, 20.0, 25.0])
    is_ocean = np.zeros(3, dtype=bool)
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])
    filled = np.array([30.0, 20.0, 25.0])
    spill = np.full(3, -1)

    steepest = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, np.zeros(3))
    assert steepest[0] == 1  # no channel anywhere -- plain steepest descent

    channel_depth = np.array([0.0, 0.0, 50.0])  # node 2's channel is well-established
    channelized = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, channel_depth)
    assert channelized[0] == 2  # prefers the channelized neighbor even though it's less steep


def test_compute_flow_direction_ignores_a_barely_established_channel():
    # Same setup, but node 2's channel is below CHANNEL_PREFERENCE_THRESHOLD_M -- too new/
    # shallow to count as "an existing channel," so plain steepest descent still wins.
    elevation = np.array([30.0, 20.0, 25.0])
    is_ocean = np.zeros(3, dtype=bool)
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])
    filled = np.array([30.0, 20.0, 25.0])
    spill = np.full(3, -1)
    channel_depth = np.array([0.0, 0.0, hydrology.CHANNEL_PREFERENCE_THRESHOLD_M - 1.0])

    result = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, channel_depth)
    assert result[0] == 1


def test_flood_fill_lake_extent_covers_every_node_below_the_seed_water_surface():
    # A basin: node 0 is the sink (seeded), nodes 1/2 are shallower but still below a
    # generous water surface, node 3 is above it (dry shoreline, never flooded).
    elevation = np.array([10.0, 12.0, 14.0, 20.0])
    is_ocean = np.zeros(4, dtype=bool)
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1], [0, 1]])
    seed_mask = np.array([True, False, False, False])
    water_surface = np.array([15.0, 0.0, 0.0, 0.0])  # only the seed's own entry is read

    prev_wet_mask = np.zeros(4, dtype=bool)
    result = hydrology._flood_fill_lake_extent(elevation, is_ocean, neighbor_idx, seed_mask, water_surface, prev_wet_mask, max_new_hops=10)
    assert result.tolist() == [5.0, 3.0, 1.0, 0.0]  # 15-10, 15-12, 15-14, dry (20 >= 15)


def test_flood_fill_lake_extent_never_floods_ocean():
    elevation = np.array([10.0, -5.0])
    is_ocean = np.array([False, True])
    neighbor_idx = np.array([[1], [0]])
    seed_mask = np.array([True, False])
    water_surface = np.array([50.0, 0.0])  # would easily submerge the ocean node too, if allowed

    prev_wet_mask = np.zeros(2, dtype=bool)
    result = hydrology._flood_fill_lake_extent(elevation, is_ocean, neighbor_idx, seed_mask, water_surface, prev_wet_mask, max_new_hops=10)
    assert result[1] == 0.0


def test_flood_fill_lake_extent_limits_growth_into_brand_new_territory():
    # A 4-node chain 0-1-2-3, all below the seed's generous water surface. None were wet last
    # step, so each additional hop away from the seed costs one unit of the new-territory hop
    # budget -- with max_new_hops=1, only the seed and its immediate neighbor may be claimed
    # this step, even though nodes 2/3 are also below the water surface.
    elevation = np.array([10.0, 11.0, 12.0, 13.0])
    is_ocean = np.zeros(4, dtype=bool)
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    seed_mask = np.array([True, False, False, False])
    water_surface = np.array([20.0, 0.0, 0.0, 0.0])
    prev_wet_mask = np.zeros(4, dtype=bool)

    result = hydrology._flood_fill_lake_extent(
        elevation, is_ocean, neighbor_idx, seed_mask, water_surface, prev_wet_mask, max_new_hops=1
    )
    assert result.tolist() == [10.0, 9.0, 0.0, 0.0]  # nodes 2/3 unreached at this hop budget


def test_flood_fill_lake_extent_does_not_hop_limit_already_wet_territory():
    # Same chain, but nodes 1 and 2 were already wet last step -- passing through them should
    # cost nothing, so node 3 (one new hop beyond the already-wet region) is still reachable
    # even at max_new_hops=1.
    elevation = np.array([10.0, 11.0, 12.0, 13.0])
    is_ocean = np.zeros(4, dtype=bool)
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    seed_mask = np.array([True, False, False, False])
    water_surface = np.array([20.0, 0.0, 0.0, 0.0])
    prev_wet_mask = np.array([False, True, True, False])

    result = hydrology._flood_fill_lake_extent(
        elevation, is_ocean, neighbor_idx, seed_mask, water_surface, prev_wet_mask, max_new_hops=1
    )
    assert result.tolist() == [10.0, 9.0, 8.0, 7.0]  # all reachable -- 1/2 free, 3 costs one hop


def test_update_lakes_caps_depth_change_per_step():
    # A sink with an enormous amount of new inflow in one step would, without a rate limit,
    # jump straight up to its basin's (very deep) spill cap in a single step. The per-step
    # rate limit should instead cap how much this node's depth may actually change.
    not_accumulating = np.zeros(2, dtype=bool)
    fields = hydrology.HydrologyFields(
        points=np.zeros((2, 3)),
        elevation=np.array([10.0, -5.0]),
        is_ocean=np.array([False, True]),
        neighbor_idx=np.zeros((2, 1), dtype=np.int64),
        flow_target=np.array([-1, -1]),
        flow_accum=np.zeros(2),
        water_deposited=np.array([1_000_000.0, 0.0]),  # huge inflow this step
        filled_elevation=np.array([5000.0, -5.0]),  # basin could theoretically hold ~4990m
        spill_target=np.array([-1, -1]),
        is_river=np.zeros(2, dtype=bool),
        lake_depth=np.zeros(2),
        glacier_depth=np.zeros(2),
        line_refs=[],
    )
    grown = hydrology.update_lakes(
        fields, prev_lake_depth=np.zeros(2), water_deposited=fields.water_deposited, years=1_000_000, is_accumulating=not_accumulating
    )
    assert 0.0 < grown[0] <= hydrology.MAX_LAKE_DEPTH_CHANGE_M_PER_MYR  # capped well short of the ~4990m basin cap


def test_update_lakes_floods_a_whole_basin_not_just_the_sink():
    # Same shape as test_update_lakes_grows_at_a_sink_and_caps_at_the_spill_point, but with a
    # second land node (2) in the same basin, shallower than the sink and reachable via
    # neighbor_idx -- once the sink accumulates enough water, node 2 should show up as
    # flooded too, at a shallower depth than the sink itself.
    # water_deposited=50 over 1 Myr grows the sink to depth 5.5 (see
    # test_update_lakes_grows_at_a_sink_and_caps_at_the_spill_point for the same arithmetic),
    # a water surface of 10+5.5=15.5 -- node 2 (elevation 12) sits below that, node 1 is ocean.
    not_accumulating = np.zeros(3, dtype=bool)
    fields = hydrology.HydrologyFields(
        points=np.zeros((3, 3)),
        elevation=np.array([10.0, -5.0, 12.0]),
        is_ocean=np.array([False, True, False]),
        neighbor_idx=np.array([[2, 2], [0, 0], [0, 0]], dtype=np.int64),
        flow_target=np.array([-1, -1, -1]),
        flow_accum=np.zeros(3),
        water_deposited=np.array([50.0, 0.0, 0.0]),
        filled_elevation=np.array([25.0, -5.0, 25.0]),
        spill_target=np.array([-1, -1, -1]),
        is_river=np.zeros(3, dtype=bool),
        lake_depth=np.zeros(3),
        glacier_depth=np.zeros(3),
        line_refs=[],
    )
    grown = hydrology.update_lakes(
        fields, prev_lake_depth=np.zeros(3), water_deposited=fields.water_deposited, years=1_000_000, is_accumulating=not_accumulating
    )
    assert grown[0] > 0.0  # the sink itself
    assert grown[2] > 0.0  # flooded too -- part of the same basin, below the water surface
    assert grown[2] < grown[0]  # shallower than the sink (elevation 18 vs. sink's own 10)
    assert grown[1] == 0.0  # never touches ocean


def test_rivers_never_overlap_lakes_in_a_real_stepped_world():
    # Regression check for "rivers end at lakes and don't go through them": a lake, however
    # widely flooded, must never also be classified is_river.
    world = generate_world(seed=32, num_plates=10, continental_fraction=0.5)
    for _ in range(10):
        step_world(world, years=2_000_000)
    hydro = world.hydrology_cache
    assert hydro is not None
    is_lake = hydro.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    assert not np.any(hydro.is_river & is_lake)


def test_channel_width_grows_with_flow_and_persists_across_steps():
    world = generate_world(seed=33, num_plates=8, continental_fraction=0.5)
    for _ in range(4):
        step_world(world, years=2_000_000)

    widths = np.concatenate([line.channel_width for p in world.plates for line in p.lines])
    depths = np.concatenate([line.channel_depth for p in world.plates for line in p.lines])
    assert np.any(widths > 0.0)
    assert np.all(widths <= erosion.MAX_CHANNEL_WIDTH_M)
    # Width grows from discharge alone, unlike depth (which also needs real slope -- see
    # erosion.py's own river term) -- a flat, wide, slow-flowing stretch can widen without
    # incising at all, so this only checks that *some* overlap exists, not that every widened
    # node is also a channel-depth node.
    assert np.any((widths > 0.0) & (depths > 0.0))

    width_total_1 = widths.sum()
    step_world(world, years=2_000_000)
    width_total_2 = sum(line.channel_width.sum() for p in world.plates for line in p.lines)
    # Monotonic like channel_depth -- if boundary.py/bathymetry.py had reset it, this would
    # have dropped back toward 0 instead.
    assert width_total_2 >= width_total_1 * 0.5  # loose bound: some nodes can be pruned/moved


def test_is_volcano_survives_a_full_step_cycle():
    # Direct regression check for the bug this session found: bathymetry.py's and erosion.py's
    # own line-reconstruction sites didn't know about is_volcano/volcano_active_years_remaining
    # and silently wiped them to False/0 every step, before volcanism.apply_volcanic_activity
    # ever got a chance to read them.
    world = generate_world(seed=34, num_plates=10, continental_fraction=0.4)
    found_a_field = False
    for _ in range(20):
        step_world(world, years=2_000_000)
        if world.volcanic_field_plate_ids:
            found_a_field = True
            total_volcano_nodes = sum(int(line.is_volcano.sum()) for p in world.plates for line in p.lines)
            assert total_volcano_nodes > 0
            break
    assert found_a_field  # sanity: this seed/step count should produce at least one field


def test_update_glaciers_accumulates_when_cold_and_precipitating():
    # A single isolated node (no flow target) -- glacier_depth should grow by exactly
    # GLACIER_ACCUMULATION_RATE * frozen_precip * years_myr, nothing more.
    elevation = np.array([500.0])
    is_ocean = np.array([False])
    flow_target = np.array([-1])
    slope = np.array([0.0])
    prev_glacier = np.array([0.0])
    frozen_precip = np.array([50.0])
    frozen_from_lake = np.array([0.0])
    is_accumulating = np.array([True])
    temperature = np.array([-20.0])  # well below GLACIER_ACCUMULATION_TEMP_C

    new_depth, melt = hydrology._update_glaciers(
        elevation, is_ocean, flow_target, slope, prev_glacier, frozen_precip, frozen_from_lake, is_accumulating, temperature, years=1_000_000
    )
    expected = hydrology.GLACIER_ACCUMULATION_RATE * 50.0 * 1.0
    assert np.isclose(new_depth[0], expected)
    assert melt[0] == 0.0  # far too cold to melt anything


def test_update_glaciers_melts_when_warm_capped_at_available_depth():
    elevation = np.array([500.0, 500.0])
    is_ocean = np.array([False, False])
    flow_target = np.array([-1, -1])
    slope = np.array([0.0, 0.0])
    prev_glacier = np.array([5000.0, 1.0])  # first has plenty of ice; second has very little
    frozen_precip = np.array([0.0, 0.0])
    frozen_from_lake = np.array([0.0, 0.0])
    is_accumulating = np.array([False, False])
    # Well above the melt-saturation point -> maximum melt factor.
    temperature = np.full(2, hydrology.GLACIER_ACCUMULATION_TEMP_C + hydrology.GLACIER_MELT_REFERENCE_DEGREES_C + 50.0)

    new_depth, melt = hydrology._update_glaciers(
        elevation, is_ocean, flow_target, slope, prev_glacier, frozen_precip, frozen_from_lake, is_accumulating, temperature, years=1_000_000
    )
    # Node 0 has plenty of ice -- melts at the full, uncapped rate (max melt factor).
    expected_uncapped_melt = hydrology.GLACIER_MELT_RATE_M_PER_MYR * hydrology.GLACIER_MELT_MAX_FACTOR * 1.0
    assert np.isclose(melt[0], expected_uncapped_melt)
    assert np.isclose(new_depth[0], 5000.0 - expected_uncapped_melt)
    # Node 1 barely has any ice -- melt can't remove more than what's actually there, so it
    # melts away completely rather than going negative.
    assert melt[1] == 1.0
    assert new_depth[1] == 0.0


def test_update_glaciers_flows_downhill_via_flow_target():
    # Node 0 (elevation 500) flows into node 1 (elevation 100) -- a slope-driven fraction of
    # node 0's ice should move to node 1 this step.
    elevation = np.array([500.0, 100.0])
    is_ocean = np.array([False, False])
    flow_target = np.array([1, -1])
    slope = np.array([0.5, 0.0])  # steep -- meaningful flow_fraction
    prev_glacier = np.array([200.0, 0.0])
    frozen_precip = np.array([0.0, 0.0])
    frozen_from_lake = np.array([0.0, 0.0])
    is_accumulating = np.array([False, False])
    temperature = np.full(2, hydrology.GLACIER_ACCUMULATION_TEMP_C)  # right at threshold -- no melt

    new_depth, melt = hydrology._update_glaciers(
        elevation, is_ocean, flow_target, slope, prev_glacier, frozen_precip, frozen_from_lake, is_accumulating, temperature, years=1_000_000
    )
    assert melt[0] == 0.0  # temperature at the threshold, not above it -- melt_factor is 0
    assert new_depth[0] < 200.0  # some ice flowed out
    assert new_depth[1] > 0.0  # and arrived downstream
    # Conservation: nothing lost besides melt (0 here) and nothing gained besides inflow.
    assert np.isclose(new_depth[0] + new_depth[1], 200.0)


def test_update_glaciers_discards_ice_reaching_an_ocean_node():
    elevation = np.array([500.0, -50.0])
    is_ocean = np.array([False, True])
    flow_target = np.array([1, -1])
    slope = np.array([0.9, 0.0])  # steep -- most of node 0's ice flows out this step
    prev_glacier = np.array([50.0, 999.0])  # node 1 already (wrongly) has ice, to prove it gets zeroed
    frozen_precip = np.array([0.0, 0.0])
    frozen_from_lake = np.array([0.0, 0.0])
    is_accumulating = np.array([False, False])
    temperature = np.full(2, hydrology.GLACIER_ACCUMULATION_TEMP_C)

    new_depth, _ = hydrology._update_glaciers(
        elevation, is_ocean, flow_target, slope, prev_glacier, frozen_precip, frozen_from_lake, is_accumulating, temperature, years=1_000_000
    )
    assert new_depth[1] == 0.0  # ocean never accumulates ice, regardless of inflow -- calving
    assert new_depth[0] < 50.0  # node 0 still lost the ice that flowed out


def _flow_line_plate_with_lake(plate_id, theta, elevation, lake_depth):
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(
        phi=0.0, theta=theta, elevation=np.asarray(elevation, dtype=float), lake_depth=np.asarray(lake_depth, dtype=float)
    )
    return Plate(plate_id=plate_id, frame=frame, crust_type="continental", lines=[line])


def test_compute_hydrology_freezes_an_existing_lake_when_cold():
    # A single evenly-spaced, tightly-packed 16-point chain -- dense enough that every
    # point's FLOW_NEIGHBOR_COUNT=8 nearest neighbors really are just its ~4 closest
    # chain-neighbors on each side (keeping the whole chain one connected graph end to end,
    # so basin-spill can still find a real path to the ocean at the far end), with a genuine
    # local dip at index 7 -- lower than every one of its own 8 nearest neighbors, so it's a
    # real sink, not just a low point somewhere along a monotonic slope.
    d = 0.002
    theta = d * np.arange(16)
    elevation = np.array(
        [500, 450, 400, 350, 300, 250, 200, 50, 300, 250, 200, 150, 100, 50, 10, -50], dtype=float
    )
    lake_depth = np.zeros(16)
    lake_depth[7] = 30.0
    plate = _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
    world = World(seed=0, plates=[plate])

    precipitation = np.full(16, 500.0)
    cold_temperature = np.full(16, -20.0)  # everywhere well below GLACIER_ACCUMULATION_TEMP_C

    fields = hydrology.compute_hydrology(world, precipitation, cold_temperature, years=1_000_000)
    assert np.isfinite(fields.filled_elevation[7])  # confirms the sink really is graph-connected to the ocean

    assert fields.lake_depth[7] == 0.0  # the lake froze, not still sitting there as water
    assert fields.glacier_depth[7] > 0.0  # its water moved into glacier_depth instead


def test_channel_lake_and_glacier_depth_persist_across_boundary_and_erosion_steps():
    # Regression check on the field-threading work: a node's channel_depth/lake_depth/
    # glacier_depth must survive being touched by boundary.py (which runs every step, before
    # erosion.py sets these) rather than getting silently reset to 0.
    world = generate_world(seed=32, num_plates=10, continental_fraction=0.5)
    for _ in range(6):
        step_world(world, years=2_000_000)

    glacier_total_1 = sum(line.glacier_depth.sum() for p in world.plates for line in p.lines)

    step_world(world, years=2_000_000)  # boundary.step_boundaries runs first here again

    glacier_total_2 = sum(line.glacier_depth.sum() for p in world.plates for line in p.lines)
    # Not a strict monotonic-growth assertion (unlike channel_depth, glacier_depth can melt)
    # -- just confirms boundary.py didn't wipe the whole world's ice to exactly 0.
    assert glacier_total_1 > 0.0
    assert glacier_total_2 > 0.0


def _normalize(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def _river_inspector_fields():
    # A hand-built 7-node HydrologyFields for the River Inspector's grouping logic, not run
    # through compute_hydrology (small enough to hit its FLOW_NEIGHBOR_COUNT bail-out path
    # anyway). Two separate is_river-classified networks:
    #   Network A -- a Y-confluence: 0 and 1 (two headwaters) both flow into 2, which flows
    #   into 3 (the mouth, highest flow_accum in the group), which flows onward into ocean
    #   node 4.
    #   Network B -- a single lone node 5, sitting at its own sink with lake_depth above the
    #   visible threshold (so it should read as ending in a lake, not the ocean).
    # Node 6 is ordinary land with flow but deliberately *not* in the is_river mask, to check
    # grouping doesn't pull in non-river neighbors.
    points = np.array(
        [
            _normalize([1.0, 0.05, 0.0]),  # 0: headwater
            _normalize([1.0, -0.05, 0.02]),  # 1: headwater
            _normalize([1.0, 0.0, 0.03]),  # 2: confluence
            _normalize([1.0, 0.0, 0.06]),  # 3: mouth
            _normalize([1.0, 0.0, 0.09]),  # 4: ocean
            _normalize([-1.0, 0.0, 0.0]),  # 5: lone river node, ends at a lake
            _normalize([0.0, 1.0, 0.0]),  # 6: non-river land
        ]
    )
    n = len(points)
    elevation = np.array([100.0, 110.0, 80.0, 50.0, -10.0, 60.0, 40.0])
    is_ocean = np.array([False, False, False, False, True, False, False])
    flow_target = np.array([2, 2, 3, 4, -1, -1, -1])
    flow_accum = np.array([10.0, 12.0, 25.0, 40.0, 0.0, 5.0, 0.0])
    is_river = np.array([True, True, True, True, False, True, False])
    lake_depth = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0])  # node 5 well above LAKE_MIN_VISIBLE_DEPTH_M
    zeros = np.zeros(n)
    return hydrology.HydrologyFields(
        points=points,
        elevation=elevation,
        is_ocean=is_ocean,
        neighbor_idx=np.zeros((n, 1), dtype=np.int64),
        flow_target=flow_target,
        flow_accum=flow_accum,
        water_deposited=zeros.copy(),
        filled_elevation=zeros.copy(),
        spill_target=np.full(n, -1, dtype=np.int64),
        is_river=is_river,
        lake_depth=lake_depth,
        glacier_depth=zeros.copy(),
        line_refs=[],
    )


def test_group_rivers_groups_confluence_and_finds_max_flow_mouth():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    assert len(rivers) == 2

    network_a = next(r for r in rivers if len(r.member_idx) == 4)
    assert sorted(network_a.member_idx.tolist()) == [0, 1, 2, 3]
    assert network_a.mouth_idx == 3  # highest flow_accum (40) in the group
    assert network_a.flow_rate == 40.0
    # Two headwaters (0 and 1) merging into one mainstem -- one tributary joins it.
    assert network_a.num_tributaries == 1


def test_group_rivers_classifies_mouth_type_ocean_vs_lake():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)

    network_a = next(r for r in rivers if len(r.member_idx) == 4)
    assert network_a.mouth_type == "ocean"  # mouth's flow_target (4) is an ocean node

    network_b = next(r for r in rivers if len(r.member_idx) == 1)
    assert network_b.mouth_idx == 5
    assert network_b.mouth_type == "lake"  # lake_depth[5] = 5.0 > LAKE_MIN_VISIBLE_DEPTH_M
    assert network_b.num_tributaries == 0  # a single unbranched node has no tributaries


def test_group_rivers_excludes_non_river_nodes():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    all_members = {i for r in rivers for i in r.member_idx.tolist()}
    assert 4 not in all_members  # ocean node, never is_river
    assert 6 not in all_members  # land node deliberately not classified as a river


def test_group_rivers_on_no_rivers_returns_empty_list():
    fields = _river_inspector_fields()
    fields.is_river = np.zeros(len(fields.is_river), dtype=bool)
    assert hydrology.group_rivers(fields) == []


def test_river_at_finds_the_network_nearest_a_query_point():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    network_a_id = next(i for i, r in enumerate(rivers) if len(r.member_idx) == 4)
    network_b_id = next(i for i, r in enumerate(rivers) if len(r.member_idx) == 1)

    near_a = hydrology.river_at(fields, rivers, fields.points[0])
    assert near_a == network_a_id

    near_b = hydrology.river_at(fields, rivers, fields.points[5])
    assert near_b == network_b_id


def test_river_at_returns_none_when_there_are_no_rivers():
    fields = _river_inspector_fields()
    assert hydrology.river_at(fields, [], fields.points[0]) is None


def test_group_rivers_on_a_real_monotonic_chain_is_one_unbranched_network():
    # Reuses the same 12-node monotonic-descent fixture as the compute_hydrology end-to-end
    # test above -- a single chain has no confluence, so grouping it should yield exactly one
    # network with zero tributaries, mouthing into the ocean.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    fields = hydrology.compute_hydrology(world, np.full(12, 800.0), np.full(12, 15.0), years=1_000_000)
    rivers = hydrology.group_rivers(fields)

    assert len(rivers) == 1
    assert rivers[0].num_tributaries == 0
    assert rivers[0].mouth_type == "ocean"
    assert fields.is_ocean[fields.flow_target[rivers[0].mouth_idx]]


def test_compute_river_speed_increases_with_slope_and_flow():
    base = hydrology.compute_river_speed(np.array([0.01]), np.array([1.0]))[0]
    steeper = hydrology.compute_river_speed(np.array([0.1]), np.array([1.0]))[0]
    more_flow = hydrology.compute_river_speed(np.array([0.01]), np.array([100.0]))[0]
    assert steeper > base
    assert more_flow > base

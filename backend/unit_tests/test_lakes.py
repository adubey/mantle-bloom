import pytest
import numpy as np

from app.hydroclimate import hydrology, lakes


def _leaves(roots):
    return [lake for lake in lakes.iter_all_lakes(roots) if not lake.children]


def test_build_lake_hierarchy_on_a_monotonic_slope_creates_no_lakes():
    # Every node's own pure steepest-descent chain reaches the ocean without ever passing
    # through a land local minimum -- no real depression anywhere, so no Lake at all.
    elevation = np.array([40.0, 30.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])

    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert roots == []


def test_build_lake_hierarchy_single_depression():
    # 4-node fixture: node 0 sits behind a rim (node 1, elevation 50), with node 2 (20) and the
    # ocean (node 3, -10) beyond it. Node 0 has no lower neighbor at all -- it's the only real
    # local minimum -- so it's the only leaf lake; nodes 1 and 2 both drain straight to the
    # ocean via steepest descent and never become part of any lake.
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])

    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(roots) == 1
    lake = roots[0]
    assert lake.children == []
    assert lake.members.tolist() == [0]
    assert lake.floor_elevation == 30.0
    assert lake.max_depth == 50.0
    assert 1 not in lake.members.tolist() and 2 not in lake.members.tolist() and 3 not in lake.members.tolist()

    # sink_node_idx/outlet_target_idx: this leaf's own local minimum (0), and the far side
    # (node 1) of the boundary edge that set max_depth -- the actual escape route once full.
    assert lake.sink_node_idx == 0
    assert lake.outlet_node_idx == 0
    assert lake.outlet_target_idx == 1


def test_interface_pass_spills_to_a_higher_center_at_the_shared_saddle():
    # The depression contains nodes 0 and 1.  Its downstream neighbour (2) has a 5 m centre,
    # but an explicitly supplied interface geometry crosses a 4 m saddle.  The lake must spill
    # at 4 m, rather than waiting for either boundary centre (10 m and 5 m) to be submerged.
    elevation = np.array([0.0, 10.0, 5.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    interface_pass = np.maximum(elevation[:, None], elevation[neighbor_idx])
    interface_pass[1, 1] = 4.0
    interface_pass[2, 0] = 4.0

    forest = lakes.build_lake_hierarchy(
        elevation, is_ocean, neighbor_idx, interface_pass_elevation=interface_pass
    )
    lake = forest[0]
    assert lake.max_depth == 4.0
    assert lake.outlet_node_idx == 1
    assert lake.outlet_target_idx == 2

    filled, spill, _ = lakes.compute_spill_routing(forest, elevation, np.zeros(4))
    assert filled[0] == 4.0
    assert spill[0] == 2


def test_interface_pass_does_not_cross_a_genuinely_higher_divide():
    elevation = np.array([0.0, 10.0, 5.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    interface_pass = np.maximum(elevation[:, None], elevation[neighbor_idx])
    interface_pass[1, 1] = 12.0
    interface_pass[2, 0] = 12.0

    forest = lakes.build_lake_hierarchy(
        elevation, is_ocean, neighbor_idx, interface_pass_elevation=interface_pass
    )
    filled, spill, _ = lakes.compute_spill_routing(forest, elevation, np.zeros(4))
    flow, spilling = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.array([10.0, 0.0, 0.0, 0.0]), filled, spill, np.zeros(4),
    )
    assert filled[0] == 12.0
    assert flow[0] == -1
    assert not spilling[0]


def test_build_lake_hierarchy_never_includes_ocean_nodes():
    # Reuse the single-depression fixture above but scan the *whole* forest (not just the
    # root) -- ocean nodes (and nodes that drain straight to the ocean) must never appear as a
    # member anywhere in the tree, even once real lakes exist.
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])

    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(roots) == 1  # sanity: this fixture does produce a real lake
    for lake in lakes.iter_all_lakes(roots):
        assert 3 not in lake.members.tolist()


def test_build_lake_hierarchy_merges_two_basins_before_either_reaches_ocean():
    # Two depressions -- node 0's catchment {0, 1} (floor 2.0) and node 3's catchment {2, 3}
    # (floor 4.0) -- meet at a saddle (the boundary edge between nodes 1 and 2, elevation 12.0)
    # that's lower than the rim leading onward to the ocean (node 4, elevation 30.0, which
    # itself drains straight to the ocean via node 5 and so never becomes a lake member). Because
    # the saddle between the two depressions is lower than the rim to the sea, they must merge
    # with each other first.
    elevation = np.array([2.0, 12.0, 9.0, 4.0, 30.0, -10.0])
    is_ocean = np.array([False, False, False, False, False, True])
    neighbor_idx = np.array(
        [
            [1, 1],  # 0 -> 1
            [0, 2],  # 1 -> 0, 2
            [1, 3],  # 2 -> 1, 3
            [2, 4],  # 3 -> 2, 4
            [3, 5],  # 4 -> 3, 5 (ocean)
            [4, 4],  # 5 (ocean)
        ]
    )

    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(roots) == 1
    root = roots[0]

    assert sorted(root.members.tolist()) == [0, 1, 2, 3]
    assert 4 not in root.members.tolist()  # the rim node itself drains straight to the ocean
    assert root.floor_elevation == 2.0
    assert root.max_depth == 30.0  # the rim to the ocean, via node 4

    # The root is a real merge: two children, each still an independently-reachable descendant
    # whose own min_depth records the saddle (12.0) where the two original basins first
    # touched, strictly below the final rim to the ocean.
    assert len(root.children) == 2
    assert root.min_depth == 12.0
    all_lakes = list(lakes.iter_all_lakes(roots))
    node0_leaf = next(lake for lake in all_lakes if lake.members.tolist() == [0, 1])
    assert node0_leaf.floor_elevation == 2.0
    node3_leaf = next(lake for lake in all_lakes if lake.members.tolist() == [2, 3])
    assert node3_leaf.floor_elevation == 4.0


def test_build_lake_hierarchy_closed_basin_has_no_max_depth():
    # A small disconnected world with no ocean node at all -- every land node drains, via pure
    # steepest descent, to the same single interior local minimum (node 0), which has no path
    # to any sea anywhere in this graph. A legitimate endorheic basin: max_depth stays None even
    # though it was never involved in any merge.
    elevation = np.array([5.0, 15.0, 8.0])
    is_ocean = np.array([False, False, False])
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])

    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(roots) == 1
    assert roots[0].max_depth is None
    assert roots[0].children == []
    assert sorted(roots[0].members.tolist()) == [0, 1, 2]
    assert roots[0].floor_elevation == 5.0


def test_build_lake_hierarchy_on_empty_world_returns_no_lakes():
    assert lakes.build_lake_hierarchy(np.zeros(0), np.zeros(0, dtype=bool), np.zeros((0, 8), dtype=np.int64)) == []


def test_build_lake_hierarchy_flat_terrain_produces_one_bounded_lake():
    # The actual reported bug this module fixes: a large, perfectly flat catchment (six nodes,
    # all at the same elevation -- no neighbor is *strictly* lower than any other, so each node
    # starts out as its own trivial singleton catchment) behind a real rim leading to the ocean.
    # The old flat-array flood-fill grew such a region incrementally, hop-limited per step, with
    # no notion of when it was actually "done" -- this algorithm instead computes the whole
    # basin's true, bounded extent in a single pass: all six flat nodes merge into one lake
    # (every merge between them happens at their own shared floor elevation, i.e. zero depth),
    # bounded by the real rim, not an ever-expanding front.
    elevation = np.array([5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, False, False, False, False, True])
    neighbor_idx = np.array(
        [
            [1, 1],
            [0, 2],
            [1, 3],
            [2, 4],
            [3, 5],
            [4, 6],
            [5, 7],
            [6, 6],
        ]
    )

    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(roots) == 1
    root = roots[0]
    assert sorted(root.members.tolist()) == [0, 1, 2, 3, 4, 5]
    assert root.floor_elevation == 5.0
    assert root.max_depth == 20.0  # bounded by the real rim, not left open-ended


def test_iter_all_lakes_walks_every_descendant():
    elevation = np.array([2.0, 12.0, 9.0, 4.0, 30.0, -10.0])
    is_ocean = np.array([False, False, False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 5], [4, 4]])

    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    all_lakes = list(lakes.iter_all_lakes(roots))
    assert len(all_lakes) == 3  # two leaves plus the one merged root
    leaf_members = sorted(m for leaf in _leaves(roots) for m in leaf.members.tolist())
    assert leaf_members == [0, 1, 2, 3]


# -- step_lakes ---------------------------------------------------------------------------

# A single real depression: node 0 (floor 10.0), rim at node 1 (25.0), ocean at node 2 (-5.0).
# Reused by several step_lakes tests below.
_SINK_ELEVATION = np.array([10.0, 25.0, -5.0])
_SINK_IS_OCEAN = np.array([False, False, True])
_SINK_NEIGHBORS = np.array([[1, 1], [0, 2], [1, 1]])

# Two depressions -- {0, 1} (floor 2.0) and {2, 3} (floor 4.0) -- meeting at a saddle
# (min_depth 12.0), eventually spilling to the ocean (node 5) via a rim (node 4, max_depth
# 30.0). Same fixture as the build_lake_hierarchy merge test above.
_MERGE_ELEVATION = np.array([2.0, 12.0, 9.0, 4.0, 30.0, -10.0])
_MERGE_IS_OCEAN = np.array([False, False, False, False, False, True])
_MERGE_NEIGHBORS = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 5], [4, 4]])

# More water than any of these whole-cell fixtures can evaporate at its rim (a node stands for
# 1/N of the sphere, ~1e14 m^2), so a lake fed this fills to its cap.
_FLOOD_M3_PER_YR = 1.0e20


def _climate(inflow_m3_per_yr, evaporation_m_per_yr=1.0):
    inflow = np.asarray(inflow_m3_per_yr, dtype=float)
    return lakes.LakeClimate(inflow, np.full(len(inflow), evaporation_m_per_yr))


def test_prognostic_lake_solver_grows_and_closes_single_basin_budget():
    forest = lakes.build_lake_hierarchy(_SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS)
    hyps = lakes.CellHypsometry.whole_cells(_SINK_ELEVATION)
    climate = lakes.LakeClimate(
        inflow_m3_per_yr=np.zeros(3),
        evaporation_m_per_yr=np.zeros(3),
        external_supply_m3_per_yr=np.array([1.0e13, 0.0, 0.0]),
        wet_precip_m3_per_yr=np.array([1.0e13, 0.0, 0.0]),
    )

    result = lakes.solve_lakes_prognostic(
        forest, _SINK_ELEVATION, np.zeros(3), 1.0, climate, hyps, substeps=8
    )

    assert result.lake_depth[0] > 0.0
    assert np.all(result.lake_depth >= 0.0)
    assert result.ledger.external_supply_m3 == pytest.approx(1.0e13)
    assert result.ledger.wet_precip_m3 == pytest.approx(result.ledger.external_supply_m3)
    assert result.ledger.dry_runoff_m3 == 0.0
    assert result.ledger.meltwater_m3 == 0.0
    assert result.ledger.storage_change_m3 == pytest.approx(1.0e13, rel=1e-10)
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e3)
    deposited_silt_m3 = float((result.silt_deposited * hyps.area_m2).sum())
    assert 0.0 < deposited_silt_m3 <= lakes.SILT_ACCUMULATION_COEFFICIENT * 1.0e13


def test_prognostic_lake_solver_shrinks_by_evaporation_without_reset():
    forest = lakes.build_lake_hierarchy(_SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS)
    hyps = lakes.CellHypsometry.whole_cells(_SINK_ELEVATION)
    previous = np.zeros(3)
    previous[0] = 2.0
    climate = lakes.LakeClimate(np.zeros(3), np.full(3, 0.5), np.zeros(3))

    result = lakes.solve_lakes_prognostic(
        forest, _SINK_ELEVATION, previous, 1.0, climate, hyps, substeps=16
    )

    assert 0.0 < result.lake_depth[0] < previous[0]
    assert result.ledger.evaporation_m3 > 0.0
    assert result.ledger.storage_change_m3 < 0.0
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e3)


def test_prognostic_lake_solver_merges_siblings_during_the_step():
    forest = lakes.build_lake_hierarchy(_MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS)
    hyps = lakes.CellHypsometry.whole_cells(_MERGE_ELEVATION)
    supply = np.array([1.0e16, 0.0, 0.0, 1.0e16, 0.0, 0.0])
    climate = lakes.LakeClimate(np.zeros(6), np.zeros(6), supply)

    result = lakes.solve_lakes_prognostic(
        forest, _MERGE_ELEVATION, np.zeros(6), 1.0, climate, hyps, substeps=32,
        is_ocean=_MERGE_IS_OCEAN,
    )

    assert any(event.kind == "merge" for event in result.events)
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e5)
    assert np.all(result.lake_depth >= 0.0)


def test_prognostic_lake_solver_splits_below_saddle_without_resetting_water():
    forest = lakes.build_lake_hierarchy(_MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS)
    hyps = lakes.CellHypsometry.whole_cells(_MERGE_ELEVATION)
    previous = np.zeros(6)
    for member in forest[0].members:
        previous[member] = hyps.mean_depth(20.0, np.array([member]))[0]
    climate = lakes.LakeClimate(np.zeros(6), np.full(6, 0.2), np.zeros(6))

    result = lakes.solve_lakes_prognostic(
        forest, _MERGE_ELEVATION, previous, 50.0, climate, hyps, substeps=32,
        is_ocean=_MERGE_IS_OCEAN,
    )

    assert any(event.kind == "split" for event in result.events)
    assert np.any(result.lake_depth > 0.0)
    assert np.all(result.lake_depth >= 0.0)
    assert result.ledger.storage_change_m3 < 0.0
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e6)


def test_prognostic_spill_uses_hierarchy_receiver_independent_of_flow_channel():
    forest = lakes.build_lake_hierarchy(_MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS)
    hyps = lakes.CellHypsometry.whole_cells(_MERGE_ELEVATION)
    upstream = next(leaf for leaf in _leaves(forest) if leaf.sink_node_idx == 0)
    receiver_idx = upstream.outlet_target_idx
    receiver_leaf = next(leaf for leaf in _leaves(forest) if receiver_idx in leaf.members)
    supply = np.zeros(6)
    supply[0] = 1.0e15
    climate = lakes.LakeClimate(np.zeros(6), np.zeros(6), supply)

    result = lakes.solve_lakes_prognostic(
        forest, _MERGE_ELEVATION, np.zeros(6), 1.0, climate, hyps, substeps=32,
        is_ocean=_MERGE_IS_OCEAN,
    )

    assert receiver_leaf.sink_node_idx != upstream.sink_node_idx
    assert result.ledger.internal_transfer_m3 > 0.0
    assert np.any(result.lake_depth[receiver_leaf.members] > 0.0)
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e5)


def test_prognostic_capped_lake_books_only_surplus_after_evaporation():
    forest = lakes.build_lake_hierarchy(_SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS)
    hyps = lakes.CellHypsometry.whole_cells(_SINK_ELEVATION)
    supply = np.array([1.0e17, 0.0, 0.0])
    climate = lakes.LakeClimate(np.zeros(3), np.full(3, 0.5), supply)

    result = lakes.solve_lakes_prognostic(
        forest, _SINK_ELEVATION, np.zeros(3), 1.0, climate, hyps,
        is_ocean=_SINK_IS_OCEAN,
    )

    lake = forest[0]
    assert result.should_spill[lake.sink_node_idx]
    assert result.overflow_m3[lake.outlet_target_idx] > 0.0
    assert result.ledger.overflow_out_m3 > result.ledger.overflow_to_ocean_m3
    assert result.ledger.evaporation_m3 > 0.0
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e6)


def test_prognostic_solver_exports_prior_storage_outside_the_lake_hierarchy():
    forest = lakes.build_lake_hierarchy(_SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS)
    hyps = lakes.CellHypsometry.whole_cells(_SINK_ELEVATION)
    previous = np.zeros(3)
    previous[2] = 1.5  # ocean-connected in this step, but stored as lake water previously

    result = lakes.solve_lakes_prognostic(
        forest, _SINK_ELEVATION, previous, 1.0, _climate(np.zeros(3)), hyps,
        is_ocean=_SINK_IS_OCEAN,
    )

    volume = previous[2] * hyps.area_m2[2]
    assert result.lake_depth[2] == 0.0
    assert result.overflow_m3[2] == pytest.approx(volume)
    assert result.ledger.overflow_out_m3 == pytest.approx(volume)
    assert result.ledger.overflow_to_ocean_m3 == pytest.approx(volume)
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e3)


def test_prognostic_solver_exports_prior_storage_when_hierarchy_is_empty():
    hyps = lakes.CellHypsometry.whole_cells(_SINK_ELEVATION)
    previous = np.zeros(3)
    previous[0] = 2.0

    result = lakes.solve_lakes_prognostic(
        [], _SINK_ELEVATION, previous, 1.0, _climate(np.zeros(3)), hyps,
        is_ocean=_SINK_IS_OCEAN,
    )

    volume = previous[0] * hyps.area_m2[0]
    assert not np.any(result.lake_depth)
    assert result.overflow_m3[0] == pytest.approx(volume)
    assert result.ledger.overflow_out_m3 == pytest.approx(volume)
    assert result.ledger.overflow_to_ocean_m3 == 0.0
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e3)


def test_prognostic_spill_status_tracks_the_final_substep():
    forest = lakes.build_lake_hierarchy(_SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS)
    hyps = lakes.CellHypsometry.whole_cells(_SINK_ELEVATION)
    lake = forest[0]
    previous = np.zeros(3)
    previous[0] = 40.0  # starts above the lake's 25 m rim
    climate = lakes.LakeClimate(np.zeros(3), np.full(3, 10.0), np.zeros(3))

    result = lakes.solve_lakes_prognostic(
        forest, _SINK_ELEVATION, previous, 1.0, climate, hyps, substeps=2,
    )

    assert result.ledger.overflow_out_m3 > 0.0  # first slice spills
    assert result.lake_depth[0] < lake.max_depth - _SINK_ELEVATION[0]
    assert not result.should_spill[lake.sink_node_idx]
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e3)


def test_prognostic_depth_capped_siblings_break_a_closed_spill_cycle():
    elevation = np.array([0.0, 0.0])
    left = lakes.Lake(
        lake_id=0, members=np.array([0]), floor_elevation=0.0, min_depth=0.0,
        max_depth=3_000.0, outlet_node_idx=0, sink_node_idx=0, outlet_target_idx=1,
    )
    right = lakes.Lake(
        lake_id=1, members=np.array([1]), floor_elevation=0.0, min_depth=0.0,
        max_depth=3_000.0, outlet_node_idx=1, sink_node_idx=1, outlet_target_idx=0,
    )
    parent = lakes.Lake(
        lake_id=2, members=np.array([0, 1]), floor_elevation=0.0, min_depth=3_000.0,
        max_depth=None, children=[left, right],
    )
    supply = np.full(2, 1.0e20)

    result = lakes.solve_lakes_prognostic(
        [parent], elevation, np.zeros(2), 1.0,
        lakes.LakeClimate(np.zeros(2), np.zeros(2), supply),
        lakes.CellHypsometry.whole_cells(elevation), substeps=1,
    )

    assert np.all(result.should_spill)
    assert result.ledger.overflow_out_m3 > 0.0
    assert result.ledger.residual_m3 == pytest.approx(0.0, abs=1.0e7)


# -- compute_spill_routing -----------------------------------------------------------------


def test_compute_spill_routing_single_lake_matches_its_own_max_depth():
    # Dry lake (prev_lake_depth all zero): a leaf is trivially "already merged" (min_depth ==
    # floor_elevation), so its own max_depth/outlet_target_idx govern immediately -- no inflow
    # history needed to know a single, un-nested basin's rim.
    forest = lakes.build_lake_hierarchy(_SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS)
    lake = forest[0]
    assert lake.sink_node_idx == 0 and lake.max_depth == 25.0

    filled, spill, rim = lakes.compute_spill_routing(forest, _SINK_ELEVATION, np.zeros(3))
    assert filled[0] == lake.max_depth == 25.0
    assert spill[0] == lake.outlet_target_idx == 1
    assert rim[0] == lake.outlet_node_idx == 0  # a single-node catchment: its own sink is its own rim
    # An ordinary, never-a-sink node keeps its own bare elevation/-1 -- nothing downhill of it
    # ever needs a rim to escape past (see this function's own docstring).
    assert filled[1] == _SINK_ELEVATION[1] and spill[1] == -1 and rim[1] == -1
    assert filled[2] == _SINK_ELEVATION[2] and spill[2] == -1 and rim[2] == -1


def test_compute_spill_routing_routes_both_original_sinks_to_the_parents_rim_once_merged():
    # The core case this function exists for: {0, 1} and {2, 3} already merged as of
    # prev_lake_depth (node 0's water surface at 12.0 == the saddle/min_depth), so *both*
    # original leaf sink nodes (0 and 3) must now resolve to the shared parent's own rim
    # (30.0, via node 4) -- not each leaf's own already-surpassed 12.0/individual target.
    forest = lakes.build_lake_hierarchy(_MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS)
    root = forest[0]
    assert root.min_depth == 12.0 and root.max_depth == 30.0 and root.outlet_target_idx == 4
    assert root.outlet_node_idx == 3  # the merged body's own rim member, on the far side from node 4

    already_merged_prev_depth = np.array([10.0, 0.0, 0.0, 8.0, 0.0, 0.0])  # both sinks at surface 12
    filled, spill, rim = lakes.compute_spill_routing(forest, _MERGE_ELEVATION, already_merged_prev_depth)
    assert filled[0] == filled[3] == 30.0
    assert spill[0] == spill[3] == 4
    # Both original sinks share the parent's own rim (node 3), not each leaf's individual one.
    assert rim[0] == rim[3] == 3

    # Not yet merged (everything dry): each leaf's own, shallower saddle governs instead.
    filled_dry, spill_dry, rim_dry = lakes.compute_spill_routing(forest, _MERGE_ELEVATION, np.zeros(6))
    assert filled_dry[0] == filled_dry[3] == 12.0
    assert spill_dry[0] == 2  # leaf {0,1}'s own outlet_target_idx, toward node 2
    assert spill_dry[3] == 1  # leaf {2,3}'s own outlet_target_idx, toward node 1
    assert rim_dry[0] == 1  # leaf {0,1}'s own outlet_node_idx
    assert rim_dry[3] == 2  # leaf {2,3}'s own outlet_node_idx


def test_perched_water_does_not_merge_parent_or_govern_its_spill_rim():
    # Node 4 is a perched sub-basin in the left child's upstream members. Its water surface is
    # above the parent saddle, but the main lake in each child remains below it.
    elevation = np.array([2.0, 12.0, 9.0, 4.0, 8.0, -10.0])
    left = lakes.Lake(
        lake_id=0, members=np.array([0, 1, 4]), floor_elevation=2.0, min_depth=2.0,
        max_depth=12.0, outlet_node_idx=1, sink_node_idx=0,
    )
    right = lakes.Lake(
        lake_id=1, members=np.array([2, 3]), floor_elevation=4.0, min_depth=4.0,
        max_depth=12.0, outlet_node_idx=2, sink_node_idx=3,
    )
    parent = lakes.Lake(
        lake_id=2, members=np.array([0, 1, 2, 3, 4]), floor_elevation=2.0, min_depth=12.0,
        max_depth=30.0, children=[left, right], outlet_node_idx=3,
    )
    prev_lake_depth = np.array([1.0, 0.0, 0.0, 1.0, 6.0, 0.0])

    assert elevation[4] + prev_lake_depth[4] == 14.0
    assert sorted(body.lake_id for body in lakes.governing_bodies([parent], elevation, prev_lake_depth)) == [0, 1]

    depth, _, _, events = lakes.resolve_lakes(
        [parent], elevation, prev_lake_depth, years=100_000,
        climate=_climate(np.zeros(len(elevation)), 1.0),
    )
    assert events == []
    assert not depth.any()


def test_equal_floor_merge_stays_connected_when_child_depths_are_zero():
    # Equal-height sinks merge at their shared floor. Their stored mean depths are zero, but
    # their dry levels still reach the saddle, so the persisted parent must not emit a new
    # merge event or route the sinks through the child rims on every step (#520).
    elevation = np.array([12.0, 12.0, 20.0, -10.0])
    left = lakes.Lake(
        lake_id=0, members=np.array([0]), floor_elevation=12.0, min_depth=12.0,
        max_depth=12.0, outlet_node_idx=0, sink_node_idx=0, outlet_target_idx=1,
    )
    right = lakes.Lake(
        lake_id=1, members=np.array([1]), floor_elevation=12.0, min_depth=12.0,
        max_depth=12.0, outlet_node_idx=1, sink_node_idx=1, outlet_target_idx=0,
    )
    parent = lakes.Lake(
        lake_id=2, members=np.array([0, 1]), floor_elevation=12.0, min_depth=12.0,
        max_depth=20.0, children=[left, right], outlet_node_idx=1, outlet_target_idx=2,
    )
    prev_lake_depth = np.zeros(len(elevation))

    assert lakes.governing_bodies([parent], elevation, prev_lake_depth) == [parent]
    filled, spill, rim = lakes.compute_spill_routing([parent], elevation, prev_lake_depth)
    assert filled[0] == filled[1] == 20.0
    assert spill[0] == spill[1] == 2
    assert rim[0] == rim[1] == 1

    depth, _, _, events = lakes.resolve_lakes(
        [parent], elevation, prev_lake_depth, years=100_000,
        climate=_climate(np.zeros(len(elevation))),
    )
    assert not depth.any()
    assert events == []
    assert parent.current_water_elevation == 12.0


def test_compute_spill_routing_endorheic_basin_never_spills():
    elevation = np.array([5.0, 15.0, 8.0])
    is_ocean = np.array([False, False, False])
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])

    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    filled, spill, rim = lakes.compute_spill_routing(forest, elevation, np.zeros(3))
    assert np.isinf(filled[0])
    assert spill[0] == -1
    assert rim[0] == -1


def test_step_lakes_fills_a_flooded_sink_to_its_spill_point_and_silts_only_its_wet_floor():
    prev_lake_depth = np.zeros(3)
    flooded = _climate([_FLOOD_M3_PER_YR, 0.0, 0.0])

    depth, silt_deposited, _is_sea, forest, events = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, years=1_000_000,
        climate=flooded,
    )
    assert depth[0] == 25.0 - 10.0  # more inflow than it can evaporate: fills to its own max_depth
    assert depth[1] == 0.0 and depth[2] == 0.0  # rim and ocean never hold water
    assert events == []
    assert len(forest) == 1 and forest[0].max_depth == 25.0 and forest[0].is_spilling
    assert silt_deposited[0] > 0.0  # some inflow settled as silt on the wet floor this step
    assert silt_deposited[1] == 0.0 and silt_deposited[2] == 0.0  # only the underwater node


def test_step_lakes_dries_a_lake_with_no_inflow():
    prev_lake_depth = np.array([5.0, 0.0, 0.0])

    depth, _, _, _, events = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, years=100_000,
        climate=_climate(np.zeros(3)),
    )
    assert depth[0] == 0.0
    assert events == []


def test_step_lakes_merges_two_basins_once_one_reaches_the_saddle():
    # Node 0's lake was below the 12.0 saddle last step (prev depth 9.0, surface 11.0); node
    # 3's lake was dry. This step's inflow fills node 0's own lake to its own cap (12.0, the
    # saddle height) -- at that point it must merge with node 3's basin rather than spilling
    # past its own cap, since 12.0 is a real internal saddle, not the true rim (30.0).
    prev_lake_depth = np.array([9.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    depth, _, _is_sea, forest, events = lakes.step_lakes(
        _MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS, prev_lake_depth, years=1_000_000,
        climate=_climate([_FLOOD_M3_PER_YR, 0.0, 0.0, 0.0, 0.0, 0.0]),
    )
    assert len(events) == 1 and events[0].kind == "merge" and "merged" in events[0].message
    assert depth[0] == 10.0  # 12.0 (the saddle) - 2.0 (floor)
    assert depth[1] == 0.0  # 12.0 - 12.0 (node 1 sits exactly at the saddle)
    assert depth[2] == 3.0  # 12.0 - 9.0
    assert depth[3] == 8.0  # 12.0 - 4.0
    assert depth[4] == 0.0 and depth[5] == 0.0  # never part of this lake at all

    assert len(forest) == 1
    assert forest[0].current_water_elevation == 12.0


def test_step_lakes_splits_a_merged_lake_once_it_recedes_below_the_saddle():
    # Start from exactly the merged state test_step_lakes_merges_... just produced (both
    # sub-basins recorded at the 12.0 saddle), then cut off the inflow, so the balanced level
    # drops well below the saddle. It must split back into two independent basins, each at its
    # own balance -- dry, with no inflow -- not refilled to the saddle with water nothing
    # supplied (#521).
    prev_lake_depth = np.array([10.0, 0.0, 3.0, 8.0, 0.0, 0.0])

    depth, _, _is_sea, forest, events = lakes.step_lakes(
        _MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS, prev_lake_depth, years=100_000,
        climate=_climate(np.zeros(6)),
    )
    assert len(events) == 1 and events[0].kind == "split" and "split" in events[0].message
    assert not depth.any()

    assert len(forest) == 1  # the tree shape itself is unchanged -- only the *active* level split
    assert forest[0].children[0].current_water_elevation == 2.0  # each child's own dry floor
    assert forest[0].children[1].current_water_elevation == 4.0

    # Next step reads the pair as separate, so it doesn't split again.
    _, _, _, _, events2 = lakes.step_lakes(
        _MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS, depth, years=100_000,
        climate=_climate(np.zeros(6)),
    )
    assert events2 == []


def test_a_merged_lake_split_by_a_falling_supply_rebalances_each_child_at_its_own_level():
    # #521: a merged body whose supply can no longer hold it at the 12.0 saddle splits. Each
    # child must then be resolved at its own balance: the fed child keeps its own water, the
    # unfed one stays dry -- not both filled to the saddle, which created ~5 M km^3 of water per
    # step on a mature world (#519). 1.5 cells' worth of evaporation holds node 0's own basin
    # full, but the merged body balances at 4.0, where node 3 floods too.
    forest = lakes.build_lake_hierarchy(_MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS)
    root = forest[0]
    prev_lake_depth = np.array([10.0, 0.0, 3.0, 8.0, 0.0, 0.0])  # merged at the saddle
    cell_m2 = 4.0 * np.pi * (lakes.PLANET_RADIUS_KM * 1000.0) ** 2 / 6
    climate = _climate([1.5 * cell_m2, 0.0, 0.0, 0.0, 0.0, 0.0])
    out_lake_depth = np.zeros(6)
    out_silt = np.zeros(6)
    out_is_sea = np.zeros(6, dtype=bool)
    events: list[lakes.LakeEvent] = []

    lakes._resolve(
        root, _MERGE_ELEVATION, prev_lake_depth, years_myr=0.1,
        out_lake_depth=out_lake_depth, out_silt_deposited=out_silt, out_lake_is_sea=out_is_sea,
        node_area_km2=100.0, events=events, climate=climate,
    )
    assert [e.kind for e in events] == ["split"]

    # Each child balanced on its own, with fresh outputs, as a never-merged pair would be.
    hyps = lakes.CellHypsometry.whole_cells(_MERGE_ELEVATION)
    expected_depth = np.zeros(6)
    expected_silt = np.zeros(6)
    for child in root.children:
        _, tier = lakes._classify_tier(child.members, prev_lake_depth, 100.0)
        level = lakes._water_balance(
            child, _MERGE_ELEVATION, 0.1, expected_silt, tier, False, np.zeros(6, dtype=bool), climate,
        )
        assert child.current_water_elevation == level
        expected_depth[child.members] = hyps.mean_depth(level, child.members)
    assert root.children[1].current_water_elevation == 4.0  # the unfed child is dry
    np.testing.assert_allclose(out_lake_depth, expected_depth)
    assert out_lake_depth[0] > 0.0  # the fed child keeps its water
    # The discarded merged resolution's silt isn't kept on top of the children's.
    np.testing.assert_allclose(out_silt, expected_silt)


def test_step_lakes_silt_raises_the_floor_and_eventually_fills_a_small_lake_in():
    # A small, low-capacity lake (max_depth - floor = 15.0) with heavy, sustained sediment
    # inflow every step. `step_lakes` returns a per-step silt increment that erosion.py folds
    # straight into real terrain elevation (mimicked here by adding it back into `elevation`).
    # Over many steps the floor should rise monotonically and the reported water depth shrink,
    # eventually to 0 once the floor has silted all the way up to the original rim (25.0).
    elevation = _SINK_ELEVATION.astype(float).copy()
    prev_lake_depth = np.zeros(3)
    flooded = _climate([_FLOOD_M3_PER_YR, 0.0, 0.0])

    floor_over_time = []
    depths_over_time = []
    for _ in range(400):
        depth, silt_deposited, _, _, _ = lakes.step_lakes(
            elevation, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, years=1_000_000,
            climate=flooded,
        )
        elevation = elevation + silt_deposited  # erosion.py: new_elevation += hydro.silt_deposited
        prev_lake_depth = depth
        floor_over_time.append(elevation[0])
        depths_over_time.append(depth[0])

    assert elevation[0] > _SINK_ELEVATION[0]  # the floor silted upward
    assert floor_over_time == sorted(floor_over_time)  # monotonic rise, never erodes back
    assert elevation[0] <= 25.0 + 1e-6  # never silts past the basin rim
    assert depths_over_time[-1] == 0.0  # lake fully silted in despite inflow never stopping
    assert depths_over_time[0] > 0.0


def _ev(kind, elevation_m, node_count=10, basin_count=2):
    return lakes.LakeEvent(kind=kind, node_count=node_count, elevation_m=elevation_m, basin_count=basin_count)


def test_summarize_lake_events_passes_real_basin_events_through_individually():
    events = [_ev("merge", -1770.0), _ev("split", -4560.0, node_count=435)]
    lines = lakes.summarize_lake_events(events, sea_level_m=0.0)
    assert lines == [events[0].message, events[1].message]


def test_summarize_lake_events_passes_a_lone_near_sea_level_event_through_verbatim():
    events = [_ev("merge", 3.0)]
    lines = lakes.summarize_lake_events(events, sea_level_m=0.0)
    assert lines == [events[0].message]


def test_summarize_lake_events_collapses_a_flood_of_near_sea_level_churn():
    events = [_ev("merge", 2.0) for _ in range(22)] + [_ev("split", -4.0) for _ in range(16)]
    lines = lakes.summarize_lake_events(events, sea_level_m=0.0)
    assert lines == ["38 transient coastal ponds churned near sea level this step (22 merged, 16 split)."]


def test_summarize_lake_events_keeps_real_events_when_it_also_aggregates_coastal_churn():
    events = [_ev("merge", -1770.0)] + [_ev("split", 5.0) for _ in range(4)]
    lines = lakes.summarize_lake_events(events, sea_level_m=0.0)
    assert lines[0] == events[0].message
    assert lines[1] == "4 transient coastal ponds churned near sea level this step (0 merged, 4 split)."


def test_summarize_lake_events_band_is_relative_to_the_current_sea_level():
    # Same +8 m surface reads as "coastal churn" only when sea level is near it.
    events = [_ev("merge", 108.0) for _ in range(5)]
    assert len(lakes.summarize_lake_events(events, sea_level_m=100.0)) == 1  # aggregated
    assert len(lakes.summarize_lake_events(events, sea_level_m=0.0)) == 5    # all "real"


# -- lake-vs-sea tiers ----------------------------------------------------------------------
#
# `_classify_tier` is tested directly with an explicit `node_area_km2` rather than through a
# realistic-sized world: `resolve_lakes` derives that value from `4*pi*R^2 / len(elevation)`,
# so any tiny synthetic fixture (this file's 3-6 node arrays) would represent an enormous
# real-world area per node and trivially cross SEA_MIN_FLOODED_AREA_KM2 the moment anything is
# wet -- calling `_classify_tier`/`_water_balance`/`_resolve` directly with a chosen
# `node_area_km2` decouples "does the tier logic work" from "how many nodes does this fixture
# happen to have".


def test_classify_tier_promotes_once_prior_flooded_extent_crosses_the_area_threshold():
    members = np.array([0, 1, 2])
    prev_lake_depth = np.array([5.0, 0.0, 3.0])  # 2 wet members (0 and 2)

    below = (lakes.SEA_MIN_FLOODED_AREA_KM2 / 2.0) - 1.0
    is_sea, tier_max_depth = lakes._classify_tier(members, prev_lake_depth, below)
    assert is_sea is False
    assert tier_max_depth == lakes.LAKE_MAX_DEPTH_M

    above = (lakes.SEA_MIN_FLOODED_AREA_KM2 / 2.0) + 1.0
    is_sea, tier_max_depth = lakes._classify_tier(members, prev_lake_depth, above)
    assert is_sea is True
    assert tier_max_depth == lakes.SEA_MAX_DEPTH_M


def test_classify_tier_a_never_yet_flooded_lake_defaults_to_lake_tier():
    # Bone dry last step -- even an enormous node_area_km2 can't promote it; it earns sea tier
    # only after it's actually held this much water.
    members = np.array([0, 1, 2])
    is_sea, tier_max_depth = lakes._classify_tier(members, np.zeros(3), node_area_km2=1e12)
    assert is_sea is False
    assert tier_max_depth == lakes.LAKE_MAX_DEPTH_M


def test_water_balance_caps_depth_at_the_lake_tier_ceiling_even_when_the_real_rim_allows_more():
    lake = lakes._make_leaf(0, [0], [0.0], sink_node_idx=0)
    lake.max_depth = None  # unbounded real rim -- the tier cap is the only thing that can bind
    elevation = np.array([0.0])
    flooded = _climate([_FLOOD_M3_PER_YR])  # huge inflow -- would blow past any real cap otherwise
    out_silt = np.zeros(1)
    out_is_sea = np.zeros(1, dtype=bool)

    new_level = lakes._water_balance(
        lake, elevation, years_myr=1.0,
        out_silt_deposited=out_silt, tier_max_depth=lakes.LAKE_MAX_DEPTH_M, is_sea=False,
        out_lake_is_sea=out_is_sea, climate=flooded,
    )
    assert new_level == lakes.LAKE_MAX_DEPTH_M
    assert not out_is_sea[0]


def test_water_balance_sea_tier_gets_a_higher_ceiling_and_marks_is_sea():
    lake = lakes._make_leaf(0, [0], [0.0], sink_node_idx=0)
    lake.max_depth = None
    elevation = np.array([0.0])
    flooded = _climate([_FLOOD_M3_PER_YR])
    out_silt = np.zeros(1)
    out_is_sea = np.zeros(1, dtype=bool)

    new_level = lakes._water_balance(
        lake, elevation, years_myr=1.0,
        out_silt_deposited=out_silt, tier_max_depth=lakes.SEA_MAX_DEPTH_M, is_sea=True,
        out_lake_is_sea=out_is_sea, climate=flooded,
    )
    assert new_level == lakes.SEA_MAX_DEPTH_M
    assert lakes.SEA_MAX_DEPTH_M > lakes.LAKE_MAX_DEPTH_M
    assert out_is_sea[0]


# Scaled-up version of _MERGE_ELEVATION (module-level fixture above): same shape/topology
# (two catchments {0,1}/{2,3} meeting at a saddle before either reaches the ocean via node 4),
# but with the saddle (2,500 m above the floor) deliberately placed *above* LAKE_MAX_DEPTH_M
# (1,800 m) while the real rim to the ocean (5,000 m) stays well clear of both caps.
_BIG_MERGE_ELEVATION = np.array([0.0, 2500.0, 2000.0, 100.0, 5000.0, -10.0])
_BIG_MERGE_IS_OCEAN = np.array([False, False, False, False, False, True])
_BIG_MERGE_NEIGHBORS = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 5], [4, 4]])


def test_a_composite_lake_held_below_its_own_saddle_by_the_tier_cap_splits():
    # Real-world motivation: this is the mechanism that doubles as the "smaller lakes should
    # split up" fix (see lakes.py's own module docstring / _water_balance's docstring) -- a
    # depth ceiling low enough to matter sits below some basins' own internal merge saddles, so
    # a composite lake that only stayed merged because nothing capped its rise falls back below
    # its own min_depth and splits, using the ordinary split mechanism build_lake_hierarchy
    # already provides.
    forest = lakes.build_lake_hierarchy(_BIG_MERGE_ELEVATION, _BIG_MERGE_IS_OCEAN, _BIG_MERGE_NEIGHBORS)
    assert len(forest) == 1
    root = forest[0]
    assert root.min_depth == 2500.0 and root.max_depth == 5000.0  # sanity: real rim is nowhere near either cap

    # Both original sub-basins already merged as of last step: node 0's surface at 0+2500=2500,
    # node 3's at 100+2400=2500 -- both exactly at the shared saddle.
    prev_lake_depth = np.zeros(6)
    prev_lake_depth[0] = 2500.0
    prev_lake_depth[3] = 2400.0
    # Enough inflow to hold the merged lake at its tier cap, below the saddle.
    flooded = _climate([_FLOOD_M3_PER_YR, 0.0, 0.0, 0.0, 0.0, 0.0])
    out_lake_depth = np.zeros(6)
    out_silt_deposited = np.zeros(6)
    out_lake_is_sea = np.zeros(6, dtype=bool)
    events: list[lakes.LakeEvent] = []

    lakes._resolve(
        root, _BIG_MERGE_ELEVATION, prev_lake_depth, years_myr=1.0,
        out_lake_depth=out_lake_depth,
        out_silt_deposited=out_silt_deposited, out_lake_is_sea=out_lake_is_sea,
        node_area_km2=100.0,  # small -- keeps this scenario at lake tier, not sea
        events=events, climate=flooded,
    )

    assert len(events) == 1 and events[0].kind == "split"
    # Each child lands at its own balance (#521): the fed child at its own lake-tier cap
    # (floor 0 + LAKE_MAX_DEPTH_M 1800 = 1800), the unfed one dry at its floor (100) -- neither
    # at the stale saddle (2500), which would read as still merged next step and re-run the
    # same capped computation forever.
    assert root.children[0].current_water_elevation == 1800.0
    assert root.children[1].current_water_elevation == 100.0
    assert out_lake_depth[0] == 1800.0  # node 0's own floor is 0
    assert out_lake_depth[3] == 0.0
    # The cap held it at lake tier throughout -- nothing here was ever big enough to promote.
    assert not out_lake_is_sea.any()

    # Confirm the fix actually breaks the oscillation: resolving *again* from this step's own
    # output (as the next real step would), now with no inflow, must NOT immediately re-merge
    # into the same losing computation -- each child now resolves independently and drains,
    # not bouncing back up to the old saddle.
    forest2 = lakes.build_lake_hierarchy(_BIG_MERGE_ELEVATION, _BIG_MERGE_IS_OCEAN, _BIG_MERGE_NEIGHBORS)
    events2: list[lakes.LakeEvent] = []
    lakes._resolve(
        forest2[0], _BIG_MERGE_ELEVATION, out_lake_depth, years_myr=1.0,
        out_lake_depth=np.zeros(6),
        out_silt_deposited=np.zeros(6), out_lake_is_sea=np.zeros(6, dtype=bool),
        node_area_km2=100.0, events=events2, climate=_climate(np.zeros(6)),
    )
    assert not any(e.kind == "split" for e in events2)  # no repeat split -- already independent
    assert forest2[0].current_water_elevation < 1800.0  # genuinely draining, not stuck at the cap


def test_a_child_fed_to_the_saddle_churns_without_creating_water():
    # _BIG_MERGE with node 3's floor raised to 1000, so that child's own tier cap (2800) clears
    # the 2500 saddle while the merged body's (0 + 1800) doesn't. Fed, the child fills to its
    # rim, the saddle, so it creates a merge after the initial split. The merged body is then
    # capped below the saddle and becomes separate again. This is real merge/split churn; the
    # false merge from a perched member is covered separately by the #520 regression above.
    elevation = _BIG_MERGE_ELEVATION.copy()
    elevation[3] = 1000.0
    climate = _climate([0.0, 0.0, 0.0, _FLOOD_M3_PER_YR, 0.0, 0.0])
    prev_lake_depth = np.zeros(6)
    prev_lake_depth[0] = 2500.0
    prev_lake_depth[3] = 1500.0  # merged at the saddle
    for step in range(3):
        forest = lakes.build_lake_hierarchy(elevation, _BIG_MERGE_IS_OCEAN, _BIG_MERGE_NEIGHBORS)
        root = forest[0]
        out_lake_depth = np.zeros(6)
        events: list[lakes.LakeEvent] = []
        lakes._resolve(
            root, elevation, prev_lake_depth, years_myr=1.0,
            out_lake_depth=out_lake_depth,
            out_silt_deposited=np.zeros(6), out_lake_is_sea=np.zeros(6, dtype=bool),
            node_area_km2=100.0, events=events, climate=climate,
        )
        # The first step splits the stale merged state. After that, only the fed child reaches
        # the saddle, so it creates a fresh merge event rather than a false split of the parent.
        assert [e.kind for e in events] == (["split"] if step == 0 else ["merge"])
        assert out_lake_depth[0] == (0.0 if step == 0 else 1800.0)
        assert root.children[1].current_water_elevation == 2500.0  # fed child: at its own rim
        assert out_lake_depth[3] == (1500.0 if step == 0 else 800.0)
        prev_lake_depth = out_lake_depth


def _hyps(elevation, relief, area=1.0e10):
    elevation = np.asarray(elevation, dtype=float)
    return lakes.CellHypsometry(elevation, np.full(len(elevation), area), np.asarray(relief, dtype=float), area)


def test_cell_hypsometry_floods_a_cell_gradually_and_continuously():
    # One cell, mean 100 m, ground spread over 80-120 m.
    hyps = _hyps([100.0], [20.0])
    idx = np.array([0])
    assert hyps.mean_depth(80.0, idx)[0] == 0.0
    assert hyps.wet_fraction(80.0, idx)[0] == 0.0
    assert hyps.wet_fraction(100.0, idx)[0] == pytest.approx(0.5)
    assert hyps.mean_depth(100.0, idx)[0] == pytest.approx(5.0)  # 20^2 / (4 * 20)
    # Continuous where the cell becomes fully flooded, and the mean depth there is h - z.
    assert hyps.mean_depth(120.0 - 1e-9, idx)[0] == pytest.approx(20.0)
    assert hyps.mean_depth(130.0, idx)[0] == pytest.approx(30.0)
    # Mean depth inverts back to the level, partly or fully flooded.
    for level in (85.0, 100.0, 119.0, 150.0):
        depth = hyps.mean_depth(level, idx)
        assert hyps.level_from_depth(depth, idx)[0] == pytest.approx(level)


def test_cell_hypsometry_volume_and_area_have_no_whole_cell_jumps():
    elevation = np.array([0.0, 10.0, 25.0, 40.0])
    hyps = _hyps(elevation, [8.0, 8.0, 8.0, 8.0])
    idx = np.arange(4)
    levels = np.linspace(-10.0, 60.0, 7001)
    area = np.array([hyps.wet_area(h, idx) for h in levels])
    volume = np.array([hyps.volume(h, idx) for h in levels])
    assert np.all(np.diff(area) >= 0.0) and np.all(np.diff(volume) >= 0.0)
    # A 1 cm rise never adds more than a sliver of one cell's area.
    assert np.max(np.diff(area)) < 0.01 * 1.0e10
    # Whole cells, by contrast, jump by a full cell area at each centre.
    whole = _hyps(elevation, np.zeros(4))
    whole_area = np.array([whole.wet_area(h, idx) for h in levels])
    assert np.max(np.diff(whole_area)) == pytest.approx(1.0e10)


def test_level_for_volume_inverts_volume():
    elevation = np.array([0.0, 10.0, 25.0, 40.0])
    hyps = _hyps(elevation, [8.0, 3.0, 0.0, 12.0])
    idx = np.arange(4)
    for level in (-5.0, 4.0, 12.0, 26.0, 55.0):
        target = hyps.volume(level, idx)
        assert hyps.level_for_volume(target, idx, hyps.dry_level(idx), 100.0) == pytest.approx(level, abs=1e-6)


def _bowl_lake():
    # A bowl: floor cell 0, rising rings of cells. All one catchment, rim far above.
    elevation = np.array([0.0, 20.0, 20.0, 40.0, 40.0, 40.0, 40.0])
    lake = lakes._make_leaf(0, list(range(7)), elevation.tolist(), sink_node_idx=0)
    lake.max_depth = 200.0
    return lake, elevation


def _bowl_climate(inflow_m3_per_yr, evaporation_m_per_yr=1.2):
    inflow = np.zeros(7)
    inflow[0] = inflow_m3_per_yr
    return lakes.LakeClimate(inflow, np.full(7, evaporation_m_per_yr))


def test_silt_settles_in_a_lake_that_only_partly_floods_its_cells():
    # The water stays below every cell centre, so whole-cell flooding would see no wet cell.
    lake, elevation = _bowl_lake()
    hyps = _hyps(elevation, np.full(7, 10.0))
    silt = np.zeros(7)
    new_level = lakes._water_balance(
        lake, elevation, years_myr=0.1,
        out_silt_deposited=silt, tier_max_depth=lakes.LAKE_MAX_DEPTH_M, is_sea=False,
        out_lake_is_sea=np.zeros(7, dtype=bool), climate=_bowl_climate(1.0e8), hypsometry=hyps,
    )
    assert new_level < 0.0
    assert silt[0] > 0.0
    assert silt[0] <= hyps.mean_depth(new_level, np.array([0]))[0]


def test_silt_volume_follows_real_inflow_whatever_the_cell_size():
    # The same bowl at two node densities: each coarse cell split into four quarter-area cells
    # with the same ground. The lake's real inflow, and so its silt volume, is the same at
    # both (issue #363). Silt sized from model-unit inflow scaled with the cell area instead.
    lake, elevation = _bowl_lake()
    coarse = _hyps(elevation, np.full(7, 10.0))
    fine_elevation = np.repeat(elevation, 4)
    fine = _hyps(fine_elevation, np.full(28, 10.0), area=coarse.area_m2[0] / 4.0)
    fine_lake = lakes._make_leaf(0, list(range(28)), fine_elevation.tolist(), sink_node_idx=0)
    fine_lake.max_depth = lake.max_depth
    fine_inflow = np.zeros(28)
    fine_inflow[:4] = 2.0e9 / 4.0

    def silt_m3(lake_, hyps, climate):
        silt = np.zeros(len(hyps.elevation))
        lakes._water_balance(
            lake_, hyps.elevation, years_myr=0.1, out_silt_deposited=silt,
            tier_max_depth=lakes.LAKE_MAX_DEPTH_M, is_sea=False,
            out_lake_is_sea=np.zeros(len(silt), dtype=bool), climate=climate, hypsometry=hyps,
        )
        return float(np.dot(silt, hyps.area_m2))

    coarse_silt = silt_m3(lake, coarse, _bowl_climate(2.0e9))
    fine_silt = silt_m3(fine_lake, fine, lakes.LakeClimate(fine_inflow, np.full(28, 1.2)))
    assert coarse_silt == pytest.approx(lakes.SILT_ACCUMULATION_COEFFICIENT * 2.0e9 * 100_000)
    assert fine_silt == pytest.approx(coarse_silt)


def test_balanced_level_evaporates_exactly_the_inflow():
    lake, elevation = _bowl_lake()
    hyps = _hyps(elevation, np.full(7, 10.0))
    climate = _bowl_climate(2.0e10)
    level = lakes.balanced_level(lake, hyps, climate, hyps.dry_level(lake.members), 200.0)
    assert lakes.lake_evaporation_m3_per_yr(level, lake.members, hyps, climate) == pytest.approx(2.0e10, rel=1e-6)


def test_subfreezing_lake_keeps_liquid_and_balances_supply_against_sublimation():
    lake, elevation = _bowl_lake()
    hyps = _hyps(elevation, np.full(7, 10.0))
    climate = lakes.LakeClimate(
        np.full(7, 2.0e10 / 7.0),
        np.zeros(7),
        np.full(7, hydrology.LAKE_ICE_SUBLIMATION_M_PER_YR),
    )
    floor = hyps.dry_level(lake.members)
    level = lakes.balanced_level(lake, hyps, climate, floor, 200.0)
    assert level > floor
    assert lakes.lake_evaporation_m3_per_yr(level, lake.members, hyps, climate) == pytest.approx(
        2.0e10, rel=1e-6
    )


def test_balanced_level_rises_with_inflow_fills_to_the_rim_and_dries_without_water():
    lake, elevation = _bowl_lake()
    hyps = _hyps(elevation, np.full(7, 10.0))
    dry = hyps.dry_level(lake.members)
    levels = [lakes.balanced_level(lake, hyps, _bowl_climate(q), dry, 200.0) for q in (0.0, 5.0e9, 2.0e10, 1.0e13)]
    assert levels[0] == dry
    assert levels[0] < levels[1] < levels[2] < levels[3]
    assert levels[3] == 200.0  # more than a full lake can evaporate: it fills and spills


def test_a_hotter_climate_holds_a_smaller_lake():
    lake, elevation = _bowl_lake()
    hyps = _hyps(elevation, np.full(7, 10.0))
    dry = hyps.dry_level(lake.members)
    cool = lakes.balanced_level(lake, hyps, _bowl_climate(2.0e10, 0.6), dry, 200.0)
    hot = lakes.balanced_level(lake, hyps, _bowl_climate(2.0e10, 1.8), dry, 200.0)
    assert hot < cool


def test_shallow_margins_evaporate_faster_but_total_evaporation_never_falls_as_the_level_rises():
    lake, elevation = _bowl_lake()
    hyps = _hyps(elevation, np.full(7, 10.0))
    climate = _bowl_climate(0.0, 1.0)
    floor = np.array([0])
    # Half-flooded floor cell: its wet half evaporates faster than open water would.
    half = lakes.lake_evaporation_m3_per_yr(0.0, floor, hyps, climate)
    assert half == pytest.approx(0.5 * 1.0e10 * (1.0 + 0.5 * lakes.SHALLOW_EVAPORATION_BOOST))
    ice_climate = lakes.LakeClimate(np.zeros(1), np.zeros(1), np.full(1, 0.6))
    sublimation = lakes.lake_evaporation_m3_per_yr(100.0, np.array([0]), _hyps([100.0], [20.0]), ice_climate)
    assert sublimation == pytest.approx(0.5 * 1.0e10 * 0.6)
    levels = np.linspace(-15.0, 60.0, 3001)
    evaporation = [lakes.lake_evaporation_m3_per_yr(h, lake.members, hyps, climate) for h in levels]
    assert np.all(np.diff(evaporation) >= -1e-6)


def test_water_balance_with_climate_uses_the_balanced_level():
    lake, elevation = _bowl_lake()
    hyps = _hyps(elevation, np.full(7, 10.0))
    climate = _bowl_climate(2.0e10)
    new_level = lakes._water_balance(
        lake, elevation, years_myr=0.1,
        out_silt_deposited=np.zeros(7), tier_max_depth=lakes.LAKE_MAX_DEPTH_M, is_sea=False,
        out_lake_is_sea=np.zeros(7, dtype=bool), climate=climate, hypsometry=hyps,
    )
    assert new_level == pytest.approx(lakes.balanced_level(lake, hyps, climate, hyps.dry_level(lake.members), 200.0))


def test_a_pit_breached_to_its_floor_drains_its_low_ground_but_equal_siblings_still_pool():
    # Breached: pit 1 merges toward the ocean at its own floor, so it holds no water, even in
    # the low ground of its partly flooded cell.
    elevation = np.array([30.0, 0.0, 0.0, -50.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    passes = np.maximum(elevation[:, None], elevation[neighbor_idx])
    passes[1, 1] = passes[2, 0] = 0.0
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx, interface_pass_elevation=passes)
    pit = [lake for lake in lakes.iter_all_lakes(forest) if 1 in lake.members.tolist()]
    assert all(lake.drains_at_floor for lake in pit)

    # Closed: two equal-floor sinks behind a 50 m rim merge at their shared floor, but the
    # merged basin can't drain, so neither leaf drains at its floor.
    elevation = np.array([50.0, 0.0, 0.0, 50.0])
    is_ocean = np.zeros(4, dtype=bool)
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert not any(lake.drains_at_floor for lake in lakes.iter_all_lakes(forest))

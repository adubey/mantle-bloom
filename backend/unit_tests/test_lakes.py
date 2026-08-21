import numpy as np

from app import hydrology, lakes


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


def test_build_lake_hierarchy_single_depression_matches_basin_spill():
    # Same 4-node fixture as hydrology.py's own
    # test_compute_basin_spill_collapses_nested_rims_to_one_hop: node 0 sits behind a rim (node
    # 1, elevation 50), with node 2 (20) and the ocean (node 3, -10) beyond it. Node 0 has no
    # lower neighbor at all -- it's the only real local minimum -- so it's the only leaf lake;
    # nodes 1 and 2 both drain straight to the ocean via steepest descent and never become part
    # of any lake.
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])

    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(roots) == 1
    lake = roots[0]
    assert lake.children == []
    assert lake.members.tolist() == [0]
    assert lake.floor_elevation == 30.0
    assert lake.max_depth == 50.0  # matches filled_elevation[0] from the basin-spill test
    assert 1 not in lake.members.tolist() and 2 not in lake.members.tolist() and 3 not in lake.members.tolist()

    filled, _ = hydrology._compute_basin_spill(elevation, is_ocean, neighbor_idx)
    assert lake.max_depth == filled[0]


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

    filled, _ = hydrology._compute_basin_spill(elevation, is_ocean, neighbor_idx)
    assert root.max_depth == filled[0]

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

    filled, _ = hydrology._compute_basin_spill(elevation, is_ocean, neighbor_idx)
    assert root.max_depth == filled[0]


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


def test_step_lakes_grows_at_a_sink_and_caps_at_the_spill_point():
    prev_lake_depth = np.zeros(3)
    prev_silt_depth = np.zeros(3)
    water_deposited = np.array([50.0, 0.0, 0.0])
    is_accumulating = np.zeros(3, dtype=bool)

    depth, silt, forest, events = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, prev_silt_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
    )
    assert depth[0] > 0.0  # grew from inflow
    assert depth[0] <= 25.0 - 10.0  # never exceeds the basin's true spill depth (its own max_depth)
    assert depth[1] == 0.0 and depth[2] == 0.0  # rim and ocean never hold water
    assert events == []
    assert len(forest) == 1 and forest[0].max_depth == 25.0
    assert silt[0] > 0.0  # some inflow settled as silt too

    # A lake already sitting at its cap should stay pinned there, not evaporate back down.
    at_cap = np.array([15.0, 0.0, 0.0])
    depth_at_cap, _, _, _ = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, at_cap, prev_silt_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
    )
    assert depth_at_cap[0] == 15.0


def test_step_lakes_evaporates_a_dry_spell_to_nothing():
    prev_lake_depth = np.array([5.0, 0.0, 0.0])
    prev_silt_depth = np.zeros(3)
    water_deposited = np.zeros(3)
    is_accumulating = np.zeros(3, dtype=bool)

    depth, _, _, events = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, prev_silt_depth, water_deposited, years=100_000_000, is_frozen=is_accumulating
    )
    assert depth[0] == 0.0
    assert events == []


def test_step_lakes_freezes_a_lake_to_its_dry_floor_regardless_of_inflow():
    prev_lake_depth = np.array([5.0, 0.0, 0.0])
    prev_silt_depth = np.zeros(3)
    water_deposited = np.array([500.0, 0.0, 0.0])  # would otherwise grow it a lot
    is_accumulating = np.array([True, False, False])

    depth, silt, _, _ = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, prev_silt_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
    )
    assert depth[0] == 0.0
    assert silt[0] == 0.0  # frozen -- no liquid water to carry sediment either


def test_step_lakes_merges_two_basins_once_one_reaches_the_saddle():
    # Node 0's lake is already close to full (prev depth 9.0, surface 11.0 < the 12.0 saddle);
    # node 3's lake is still dry. This step's inflow pushes node 0's own lake over its own cap
    # (12.0, the saddle height) -- at that point it must merge with node 3's basin rather than
    # spilling past its own cap, since 12.0 is a real internal saddle, not the true rim (30.0).
    prev_lake_depth = np.array([9.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    prev_silt_depth = np.zeros(6)
    water_deposited = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    is_accumulating = np.zeros(6, dtype=bool)

    depth, _, forest, events = lakes.step_lakes(
        _MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS, prev_lake_depth, prev_silt_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
    )
    assert len(events) == 1 and "merged" in events[0]
    assert depth[0] == 10.0  # 12.0 (the saddle) - 2.0 (floor)
    assert depth[1] == 0.0  # 12.0 - 12.0 (node 1 sits exactly at the saddle)
    assert depth[2] == 3.0  # 12.0 - 9.0
    assert depth[3] == 8.0  # 12.0 - 4.0
    assert depth[4] == 0.0 and depth[5] == 0.0  # never part of this lake at all

    assert len(forest) == 1
    assert forest[0].current_water_elevation == 12.0


def test_step_lakes_splits_a_merged_lake_once_it_recedes_below_the_saddle():
    # Start from exactly the merged state test_step_lakes_merges_... just produced (both
    # sub-basins recorded at the 12.0 saddle), then apply a long, dry step -- strong enough
    # evaporation to drop well below the saddle. It must split back into two independent
    # basins, both starting exactly at the saddle height (continuity), not something lower.
    prev_lake_depth = np.array([10.0, 0.0, 3.0, 8.0, 0.0, 0.0])
    prev_silt_depth = np.zeros(6)
    water_deposited = np.zeros(6)
    is_accumulating = np.zeros(6, dtype=bool)

    depth, _, forest, events = lakes.step_lakes(
        _MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS, prev_lake_depth, prev_silt_depth, water_deposited, years=5_000_000, is_frozen=is_accumulating
    )
    assert len(events) == 1 and "split" in events[0]
    assert depth[0] == 10.0  # 12.0 (the saddle) - 2.0 (floor) -- both children land exactly there
    assert depth[1] == 0.0
    assert depth[2] == 3.0  # 12.0 - 9.0
    assert depth[3] == 8.0  # 12.0 - 4.0

    assert len(forest) == 1  # the tree shape itself is unchanged -- only the *active* level split
    assert forest[0].children[0].current_water_elevation == 12.0
    assert forest[0].children[1].current_water_elevation == 12.0


def test_step_lakes_silt_accumulates_and_eventually_fills_a_small_lake_in():
    # A small, low-capacity lake (max_depth - floor = 15.0) with heavy, sustained sediment
    # inflow every step -- repeated over many steps, silt should build up at its own floor and
    # the reported water depth should shrink even though nothing about the terrain's real
    # elevation or the inflow itself ever changes, eventually reaching 0 once the lake has
    # silted all the way up to its own rim.
    prev_lake_depth = np.zeros(3)
    prev_silt_depth = np.zeros(3)
    water_deposited = np.array([500.0, 0.0, 0.0])
    is_accumulating = np.zeros(3, dtype=bool)

    depths_over_time = []
    for _ in range(400):
        depth, prev_silt_depth, _, _ = lakes.step_lakes(
            _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, prev_silt_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
        )
        prev_lake_depth = depth
        depths_over_time.append(depth[0])

    assert prev_silt_depth[0] > 0.0
    # Once the floor has silted up to the original rim (25.0), the lake is gone even though
    # water_deposited never stopped arriving.
    assert depths_over_time[-1] == 0.0
    assert prev_silt_depth[0] >= 25.0 - 10.0

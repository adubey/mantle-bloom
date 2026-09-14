import numpy as np

from app import lakes


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

    already_merged_prev_depth = np.array([10.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # node 0 surface: 2+10=12
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


def test_compute_spill_routing_endorheic_basin_never_spills():
    elevation = np.array([5.0, 15.0, 8.0])
    is_ocean = np.array([False, False, False])
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])

    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    filled, spill, rim = lakes.compute_spill_routing(forest, elevation, np.zeros(3))
    assert np.isinf(filled[0])
    assert spill[0] == -1
    assert rim[0] == -1


def test_step_lakes_grows_at_a_sink_and_caps_at_the_spill_point():
    prev_lake_depth = np.zeros(3)
    water_deposited = np.array([50.0, 0.0, 0.0])
    is_accumulating = np.zeros(3, dtype=bool)

    depth, silt_deposited, _is_sea, forest, events = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
    )
    assert depth[0] > 0.0  # grew from inflow
    assert depth[0] <= 25.0 - 10.0  # never exceeds the basin's true spill depth (its own max_depth)
    assert depth[1] == 0.0 and depth[2] == 0.0  # rim and ocean never hold water
    assert events == []
    assert len(forest) == 1 and forest[0].max_depth == 25.0
    assert silt_deposited[0] > 0.0  # some inflow settled as silt on the wet floor this step
    assert silt_deposited[1] == 0.0 and silt_deposited[2] == 0.0  # only the underwater node

    # A lake already sitting at its cap should stay pinned there, not evaporate back down.
    at_cap = np.array([15.0, 0.0, 0.0])
    depth_at_cap, _, _, _, _ = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, at_cap, water_deposited, years=1_000_000, is_frozen=is_accumulating
    )
    assert depth_at_cap[0] == 15.0


def test_step_lakes_evaporates_a_dry_spell_to_nothing():
    prev_lake_depth = np.array([5.0, 0.0, 0.0])
    water_deposited = np.zeros(3)
    is_accumulating = np.zeros(3, dtype=bool)

    depth, _, _, _, events = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, water_deposited, years=100_000_000, is_frozen=is_accumulating
    )
    assert depth[0] == 0.0
    assert events == []


def test_step_lakes_freezes_a_lake_to_its_dry_floor_but_still_silts_it_in():
    prev_lake_depth = np.array([5.0, 0.0, 0.0])
    water_deposited = np.array([500.0, 0.0, 0.0])  # would otherwise grow it a lot
    is_accumulating = np.array([True, False, False])

    depth, silt_deposited, _, _, _ = lakes.step_lakes(
        _SINK_ELEVATION, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
    )
    assert depth[0] == 0.0  # frozen -- no visible standing water reported
    # ...but real inflow still carries a sediment load that settles onto the bed even under
    # ice (GitHub issue #117: without this a chronically-frozen pit never silts in at all).
    assert silt_deposited[0] > 0.0


def test_step_lakes_merges_two_basins_once_one_reaches_the_saddle():
    # Node 0's lake is already close to full (prev depth 9.0, surface 11.0 < the 12.0 saddle);
    # node 3's lake is still dry. This step's inflow pushes node 0's own lake over its own cap
    # (12.0, the saddle height) -- at that point it must merge with node 3's basin rather than
    # spilling past its own cap, since 12.0 is a real internal saddle, not the true rim (30.0).
    prev_lake_depth = np.array([9.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    water_deposited = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    is_accumulating = np.zeros(6, dtype=bool)

    depth, _, _is_sea, forest, events = lakes.step_lakes(
        _MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS, prev_lake_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
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
    # sub-basins recorded at the 12.0 saddle), then apply a long, dry step -- strong enough
    # evaporation to drop well below the saddle. It must split back into two independent
    # basins, both starting exactly at the saddle height (continuity), not something lower.
    prev_lake_depth = np.array([10.0, 0.0, 3.0, 8.0, 0.0, 0.0])
    water_deposited = np.zeros(6)
    is_accumulating = np.zeros(6, dtype=bool)

    depth, _, _is_sea, forest, events = lakes.step_lakes(
        _MERGE_ELEVATION, _MERGE_IS_OCEAN, _MERGE_NEIGHBORS, prev_lake_depth, water_deposited, years=5_000_000, is_frozen=is_accumulating
    )
    assert len(events) == 1 and events[0].kind == "split" and "split" in events[0].message
    assert depth[0] == 10.0  # 12.0 (the saddle) - 2.0 (floor) -- both children land exactly there
    assert depth[1] == 0.0
    assert depth[2] == 3.0  # 12.0 - 9.0
    assert depth[3] == 8.0  # 12.0 - 4.0

    assert len(forest) == 1  # the tree shape itself is unchanged -- only the *active* level split
    assert forest[0].children[0].current_water_elevation == 12.0
    assert forest[0].children[1].current_water_elevation == 12.0


def test_step_lakes_silt_raises_the_floor_and_eventually_fills_a_small_lake_in():
    # A small, low-capacity lake (max_depth - floor = 15.0) with heavy, sustained sediment
    # inflow every step. `step_lakes` returns a per-step silt increment that erosion.py folds
    # straight into real terrain elevation (mimicked here by adding it back into `elevation`).
    # Over many steps the floor should rise monotonically and the reported water depth shrink,
    # eventually to 0 once the floor has silted all the way up to the original rim (25.0).
    elevation = _SINK_ELEVATION.astype(float).copy()
    prev_lake_depth = np.zeros(3)
    water_deposited = np.array([500.0, 0.0, 0.0])
    is_accumulating = np.zeros(3, dtype=bool)

    floor_over_time = []
    depths_over_time = []
    for _ in range(400):
        depth, silt_deposited, _, _, _ = lakes.step_lakes(
            elevation, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, water_deposited, years=1_000_000, is_frozen=is_accumulating
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


def test_water_balance_frozen_lake_reports_no_water_but_still_deposits_silt():
    # GitHub issue #117: a chronically-frozen catchment used to get zero silt forever
    # regardless of inflow, since the old code returned before computing anything. It should
    # still report a dry floor (no visible standing water while frozen)...
    lake = lakes._make_leaf(0, [0], [0.0], sink_node_idx=0)
    elevation = np.array([0.0])
    water_deposited = np.array([500.0])
    out_silt = np.zeros(1)
    out_is_sea = np.zeros(1, dtype=bool)

    new_level = lakes._water_balance(
        lake, 0.0, elevation, water_deposited, years_myr=1.0, is_frozen=True,
        out_silt_deposited=out_silt, tier_max_depth=lakes.LAKE_MAX_DEPTH_M, is_sea=True,
        out_lake_is_sea=out_is_sea,
    )
    # ...but still silts in from this step's inflow, and never gets marked as open sea while
    # frozen even if it's big enough to otherwise qualify.
    assert new_level == lake.floor_elevation
    assert out_silt[0] > 0.0
    assert not out_is_sea[0]


def test_step_lakes_a_chronically_frozen_pit_still_silts_in_over_time():
    # Same shape as test_step_lakes_silt_raises_the_floor_and_eventually_fills_a_small_lake_in,
    # but frozen every single step -- the exact "never holds water" caterpillar-tree case from
    # GitHub issue #117. Before the fix this lake's floor never moved at all.
    elevation = _SINK_ELEVATION.astype(float).copy()
    prev_lake_depth = np.zeros(3)
    water_deposited = np.array([500.0, 0.0, 0.0])
    is_frozen = np.ones(3, dtype=bool)

    for _ in range(400):
        depth, silt_deposited, _, _, _ = lakes.step_lakes(
            elevation, _SINK_IS_OCEAN, _SINK_NEIGHBORS, prev_lake_depth, water_deposited, years=1_000_000, is_frozen=is_frozen
        )
        elevation = elevation + silt_deposited
        prev_lake_depth = depth

    assert elevation[0] > _SINK_ELEVATION[0]  # the floor silted upward despite being frozen throughout
    assert (depth == 0.0).all()  # never reports standing water while frozen


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
    water_deposited = np.array([1e9])  # huge inflow -- would blow past any real cap otherwise
    out_silt = np.zeros(1)
    out_is_sea = np.zeros(1, dtype=bool)

    new_level = lakes._water_balance(
        lake, 0.0, elevation, water_deposited, years_myr=1.0, is_frozen=False,
        out_silt_deposited=out_silt, tier_max_depth=lakes.LAKE_MAX_DEPTH_M, is_sea=False,
        out_lake_is_sea=out_is_sea,
    )
    assert new_level == lakes.LAKE_MAX_DEPTH_M
    assert not out_is_sea[0]


def test_water_balance_sea_tier_gets_a_higher_ceiling_and_marks_is_sea():
    lake = lakes._make_leaf(0, [0], [0.0], sink_node_idx=0)
    lake.max_depth = None
    elevation = np.array([0.0])
    water_deposited = np.array([1e9])
    out_silt = np.zeros(1)
    out_is_sea = np.zeros(1, dtype=bool)

    new_level = lakes._water_balance(
        lake, 0.0, elevation, water_deposited, years_myr=1.0, is_frozen=False,
        out_silt_deposited=out_silt, tier_max_depth=lakes.SEA_MAX_DEPTH_M, is_sea=True,
        out_lake_is_sea=out_is_sea,
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
    water_deposited = np.zeros(6)
    out_lake_depth = np.zeros(6)
    out_silt_deposited = np.zeros(6)
    out_lake_is_sea = np.zeros(6, dtype=bool)
    events: list[lakes.LakeEvent] = []

    lakes._resolve(
        root, _BIG_MERGE_ELEVATION, prev_lake_depth, water_deposited, years_myr=1.0,
        is_frozen=np.zeros(6, dtype=bool), out_lake_depth=out_lake_depth,
        out_silt_deposited=out_silt_deposited, out_lake_is_sea=out_lake_is_sea,
        node_area_km2=100.0,  # small -- keeps this scenario at lake tier, not sea
        events=events,
    )

    assert len(events) == 1 and events[0].kind == "split"
    # Both children land at the tier-capped level (floor 0 + LAKE_MAX_DEPTH_M 1800 = 1800), not
    # the stale old saddle (2500) -- see _resolve's own comment on why resetting to min_depth
    # here (the *ordinary* organic-dip case's continuity approximation) would completely undo
    # the cap and oscillate forever: the very next step's _prev_level check would read exactly
    # min_depth again and re-merge into the same failing computation.
    assert root.children[0].current_water_elevation == 1800.0
    assert root.children[1].current_water_elevation == 1800.0
    assert out_lake_depth[0] == 1800.0  # node 0's own floor is 0
    assert out_lake_depth[3] == 1700.0  # node 3's own floor is 100 -> 1800 - 100
    # The cap held it at lake tier throughout -- nothing here was ever big enough to promote.
    assert not out_lake_is_sea.any()

    # Confirm the fix actually breaks the oscillation: resolving *again* from this step's own
    # output (as the next real step would) must NOT immediately re-merge into the same losing
    # computation -- each child now resolves independently and keeps decaying, not bouncing
    # back up to the old saddle.
    forest2 = lakes.build_lake_hierarchy(_BIG_MERGE_ELEVATION, _BIG_MERGE_IS_OCEAN, _BIG_MERGE_NEIGHBORS)
    events2: list[lakes.LakeEvent] = []
    lakes._resolve(
        forest2[0], _BIG_MERGE_ELEVATION, out_lake_depth, water_deposited, years_myr=1.0,
        is_frozen=np.zeros(6, dtype=bool), out_lake_depth=np.zeros(6),
        out_silt_deposited=np.zeros(6), out_lake_is_sea=np.zeros(6, dtype=bool),
        node_area_km2=100.0, events=events2,
    )
    assert not any(e.kind == "split" for e in events2)  # no repeat split -- already independent
    assert forest2[0].current_water_elevation < 1800.0  # genuinely decaying, not stuck at the cap

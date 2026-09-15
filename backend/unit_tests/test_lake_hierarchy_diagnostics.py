import json

import numpy as np
from app import lakes, persistence
from app.lake_hierarchy_diagnostics import (
    build_report,
    format_report,
    hierarchy_depths,
    leaf_catchment_sizes,
    main,
)
from app.world import generate_world, step_world


def test_leaf_catchment_sizes_counts_leaves_only_not_merged_parents():
    # Same two-basin-merges-before-the-ocean fixture as
    # test_lakes.test_build_lake_hierarchy_merges_two_basins_before_either_reaches_ocean: two
    # depressions -- {0, 1} (floor 2.0) and {2, 3} (floor 4.0) -- merge with each other at
    # saddle 12.0 before either reaches the ocean rim (node 4, elevation 30.0). The forest is
    # one root (the merged parent, `members` = the union of both leaves) with two 2-node leaf
    # children; only the two leaf sizes should be counted, not the 4-node parent's.
    elevation = np.array([2.0, 12.0, 9.0, 4.0, 30.0, -10.0])
    is_ocean = np.array([False, False, False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 5], [4, 4]])
    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(roots) == 1 and len(roots[0].children) == 2  # sanity: this fixture does merge

    sizes = leaf_catchment_sizes(roots)
    assert sorted(sizes) == [2, 2]
    assert hierarchy_depths(roots) == [2]


def test_hierarchy_depths_is_one_for_a_bare_leaf_root():
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert roots[0].children == []  # sanity: a bare leaf, no merge

    assert hierarchy_depths(roots) == [1]


def test_hierarchy_depths_counts_chain_length_through_nested_merges():
    # Extends the two-basin fixture above with a third catchment merging in afterward, so the
    # forest has one root of depth 3: {0,1} (floor 2.0) and {2,3,4} (floor 4.0) merge at saddle
    # 12.0 into parent P1 (depth 2); P1 then merges with leaf {5} (floor 6.0) at a higher
    # saddle (edge 4-5, weight 30.0) into root P2, which finally reaches the ocean (node 7) via
    # the rim at node 6 (elevation 50.0). Longest root-to-leaf chain: P2 -> P1 -> either
    # 2-node/3-node leaf -> 3 levels.
    elevation = np.array([2.0, 12.0, 9.0, 4.0, 30.0, 6.0, 50.0, -10.0])
    is_ocean = np.array([False, False, False, False, False, False, False, True])
    neighbor_idx = np.array(
        [
            [1, 1],  # 0
            [0, 2],  # 1
            [1, 3],  # 2
            [2, 4],  # 3
            [3, 5],  # 4 -- now bridges to the third catchment instead of the ocean directly
            [4, 6],  # 5 -- third catchment's own local minimum
            [5, 7],  # 6 -- rim to the ocean
            [6, 6],  # 7 (ocean)
        ]
    )
    roots = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(roots) == 1
    root = roots[0]
    assert len(root.children) == 2  # sanity: P1 and the third leaf, not a 3-way flat merge
    assert sorted(leaf_catchment_sizes(roots)) == [1, 2, 3]

    depths = hierarchy_depths(roots)
    assert depths == [3]


def test_build_report_json_serializable_and_formats(capsys):
    world = generate_world(seed=7, num_plates=6)
    step_world(world, 3_000_000)
    report = build_report(world)
    json.dumps(report)  # no numpy scalars left in
    assert report["seed"] == 7
    assert report["have_hydrology_snapshot"] is True
    assert report["max_hierarchy_depth"] >= 1 or report["num_roots"] == 0
    text = format_report(report)
    assert "mantle-bloom lake-hierarchy diagnostics" in text


def test_build_report_handles_a_world_never_stepped():
    world = generate_world(seed=7, num_plates=6)
    report = build_report(world)
    assert report["have_hydrology_snapshot"] is False
    assert report["num_roots"] == 0
    assert "never stepped" in format_report(report)


def test_main_round_trips_a_saved_file(tmp_path, capsys):
    world = generate_world(seed=7, num_plates=6)
    step_world(world, 2_000_000)
    save = tmp_path / "seed7.mbworld"
    save.write_bytes(persistence.save_world_bytes(world))

    assert main([str(save)]) == 0
    assert "lake-hierarchy diagnostics" in capsys.readouterr().out
    assert main([str(save), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["seed"] == 7

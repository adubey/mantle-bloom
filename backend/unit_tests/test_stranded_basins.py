import json

import numpy as np
import pytest
from app import lakes, persistence, stranded_basins
from app.stranded_basins import (
    StrandedBasinTrack,
    build_report,
    enrich_with_persistence,
    find_stranded_basins,
    format_report,
    main,
    reconcile_world_tracks,
)
from app.world import generate_world, step_world


def _unit_points(n: int) -> np.ndarray:
    """n distinct points spread along a short arc near (1, 0, 0) -- centroid/lat-lon math just
    needs unit vectors, not a real node cloud."""
    ang = np.linspace(0.0, 0.3, n)
    return np.stack([np.cos(ang), np.sin(ang), np.zeros(n)], axis=1)


def test_find_stranded_basins_flags_an_enclosed_below_sea_level_pit():
    # No ocean anywhere in the graph: every node drains by pure steepest descent to node 0
    # (floor -1800 m), which has no path to any sea -> an endorheic root with max_depth None.
    elevation = np.array([-1800.0, 40.0, 25.0])
    is_ocean = np.array([False, False, False])
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert len(forest) == 1 and forest[0].max_depth is None  # sanity

    lake_depth = np.array([1500.0, 0.0, 0.0])  # node 0 is flooded
    basins = find_stranded_basins(forest, elevation, _unit_points(3), lake_depth, sea_level_m=0.0)

    assert len(basins) == 1
    b = basins[0]
    assert b.floor_elevation_m == -1800.0
    assert b.depth_below_sea_level_m == 1800.0
    assert b.catchment_node_count == 3
    assert b.flooded_node_count == 1
    assert b.water_elevation_m == pytest.approx(forest[0].current_water_elevation)
    assert -90.0 <= b.centroid_lat_deg <= 90.0 and -180.0 <= b.centroid_lon_deg <= 180.0
    # straight out of find_*, persistence is unresolved
    assert b.first_seen_years is None and b.persisted_years is None


def test_find_stranded_basins_ignores_ocean_connected_and_above_sea_level_basins():
    # A real depression (node 0, floor -500) but with a rim (node 1) leading on to the ocean
    # (node 2) -- max_depth gets set, so it drains eventually and is NOT stranded.
    elevation = np.array([-500.0, 30.0, -10.0])
    is_ocean = np.array([False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 1]])
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    assert forest[0].max_depth == 30.0  # sanity: reached the ocean
    assert find_stranded_basins(forest, elevation, _unit_points(3), np.zeros(3), sea_level_m=0.0) == []

    # Enclosed (max_depth None) but its floor is above sea level -- an ordinary closed basin,
    # not the land-locked coastal pit this report is about.
    high = np.array([120.0, 200.0, 150.0])
    high_forest = lakes.build_lake_hierarchy(high, np.zeros(3, dtype=bool), np.array([[1, 2], [0, 2], [0, 1]]))
    assert high_forest[0].max_depth is None
    assert find_stranded_basins(high_forest, high, _unit_points(3), np.zeros(3), sea_level_m=0.0) == []


def test_find_stranded_basins_sorts_deepest_first():
    forest = [
        lakes.Lake(lake_id=0, members=np.array([0]), floor_elevation=-200.0, min_depth=-200.0, max_depth=None),
        lakes.Lake(lake_id=1, members=np.array([1]), floor_elevation=-4000.0, min_depth=-4000.0, max_depth=None),
        lakes.Lake(lake_id=2, members=np.array([2]), floor_elevation=-1500.0, min_depth=-1500.0, max_depth=None),
    ]
    basins = find_stranded_basins(forest, np.array([-200.0, -4000.0, -1500.0]), _unit_points(3), np.zeros(3), 0.0)
    assert [b.floor_elevation_m for b in basins] == [-4000.0, -1500.0, -200.0]


def _one_basin(centroid, floor=-1000.0):
    return stranded_basins.StrandedBasin(
        floor_elevation_m=floor,
        depth_below_sea_level_m=-floor,
        catchment_node_count=10,
        flooded_node_count=5,
        water_elevation_m=floor + 100.0,
        centroid_xyz=tuple(float(c) for c in centroid),
        centroid_lat_deg=0.0,
        centroid_lon_deg=0.0,
        floor_xyz=(1.0, 0.0, 0.0),
    )


def test_enrich_with_persistence_matches_by_centroid_and_carries_first_seen():
    here = np.array([1.0, 0.0, 0.0])
    nudged = np.array([np.cos(0.01), np.sin(0.01), 0.0])  # ~0.01 rad away, well inside the gate
    far = np.array([0.0, 1.0, 0.0])

    tracks = [StrandedBasinTrack(centroid_xyz=here, first_seen_years=1_000_000.0, last_seen_years=5_000_000.0, steps_seen=40)]
    basins = [_one_basin(nudged), _one_basin(far)]
    enrich_with_persistence(basins, tracks, elapsed_years=6_000_000.0)

    matched, fresh = basins
    assert matched.first_seen_years == 1_000_000.0
    assert matched.persisted_years == 5_000_000.0
    assert matched.steps_seen == 41
    assert fresh.first_seen_years == 6_000_000.0
    assert fresh.persisted_years == 0.0
    assert fresh.steps_seen == 1


def test_enrich_with_persistence_does_not_reuse_one_track_for_two_basins():
    here = np.array([1.0, 0.0, 0.0])
    tracks = [StrandedBasinTrack(centroid_xyz=here, first_seen_years=0.0, last_seen_years=1.0, steps_seen=1)]
    a = _one_basin(here)
    b = _one_basin(np.array([np.cos(0.005), np.sin(0.005), 0.0]))
    enrich_with_persistence([a, b], tracks, elapsed_years=2_000_000.0)
    # exactly one of them matched the single track; the other is treated as new
    assert sorted([a.steps_seen, b.steps_seen]) == [1, 2]


def test_reconcile_world_tracks_accumulates_persistence_across_steps():
    world = generate_world(seed=7, num_plates=6)
    step_world(world, 2_000_000)
    first = reconcile_world_tracks(world)  # step_world already called it once; call again is idempotent-ish
    n_tracks = len(world.stranded_basin_tracks)
    assert n_tracks == len(first)
    for basin in first:
        assert basin.steps_seen >= 1
        assert basin.persisted_years is not None


def test_build_report_json_serializable_and_formats(capsys):
    world = generate_world(seed=7, num_plates=6)
    step_world(world, 3_000_000)
    report = build_report(world)
    json.dumps(report)  # no numpy scalars left in
    assert report["seed"] == 7
    assert report["have_hydrology_snapshot"] is True
    text = format_report(report)
    assert "mantle-bloom stranded-basin diagnostics" in text


def test_build_report_handles_a_world_never_stepped():
    world = generate_world(seed=7, num_plates=6)
    report = build_report(world)
    assert report["have_hydrology_snapshot"] is False
    assert report["stranded_basins"] == []
    assert "never stepped" in format_report(report)


def test_main_round_trips_a_saved_file(tmp_path, capsys):
    world = generate_world(seed=7, num_plates=6)
    step_world(world, 2_000_000)
    save = tmp_path / "seed7.mbworld"
    save.write_bytes(persistence.save_world_bytes(world))

    assert main([str(save)]) == 0
    assert "stranded-basin diagnostics" in capsys.readouterr().out
    assert main([str(save), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["seed"] == 7


def test_old_save_without_the_field_still_loads(tmp_path):
    world = generate_world(seed=7, num_plates=6)
    del world.stranded_basin_tracks  # simulate a pickle written before the field existed
    save = tmp_path / "old.mbworld"
    save.write_bytes(persistence.save_world_bytes(world))
    loaded = persistence.load_world_bytes(save.read_bytes())
    assert loaded.stranded_basin_tracks == []

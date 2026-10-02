import pickle

import numpy as np
import pytest
from app import persistence
from app.world import generate_world, step_world


def test_round_trip_preserves_a_freshly_generated_world():
    world = generate_world(seed=42, num_plates=6, surface="lines")
    data = persistence.save_world_bytes(world)
    loaded = persistence.load_world_bytes(data)

    assert loaded is not world
    assert loaded.seed == world.seed
    assert loaded.elapsed_years == world.elapsed_years
    assert len(loaded.plates) == len(world.plates)
    for original_plate, loaded_plate in zip(world.plates, loaded.plates):
        assert loaded_plate.plate_id == original_plate.plate_id
        assert loaded_plate.crust_type == original_plate.crust_type
        assert np.allclose(loaded_plate.frame, original_plate.frame)
        assert len(loaded_plate.lines) == len(original_plate.lines)
        for original_line, loaded_line in zip(original_plate.lines, loaded_plate.lines):
            assert np.allclose(loaded_line.theta, original_line.theta)
            assert np.allclose(loaded_line.elevation, original_line.elevation)


def test_round_trip_preserves_state_only_a_step_would_populate():
    # collision_progress/climate_cache/hydrology_cache/events are all empty/None on a freshly
    # generated world -- step it first so the round trip has to actually carry a dict and the
    # two cache dataclasses, not just empty defaults.
    world = generate_world(seed=7, num_plates=6)
    step_world(world, 5_000_000)
    step_world(world, 5_000_000)

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))

    assert loaded.elapsed_years == world.elapsed_years
    assert loaded.steps_taken == world.steps_taken == 2
    assert loaded.collision_progress == world.collision_progress
    assert loaded.events == world.events
    assert (loaded.climate_cache is None) == (world.climate_cache is None)
    assert (loaded.hydrology_cache is None) == (world.hydrology_cache is None)
    if world.climate_cache is not None:
        assert np.allclose(loaded.climate_cache.elevation_m, world.climate_cache.elevation_m)
    if world.hydrology_cache is not None:
        assert np.allclose(loaded.hydrology_cache.elevation, world.hydrology_cache.elevation)


def test_loading_a_world_pickled_before_steps_taken_existed_defaults_to_zero():
    # World.steps_taken is a new field; a pickle from before it existed carries no such key.
    # A plain-int dataclass default is a class attribute, so the load still succeeds and
    # reads 0 -- see World.steps_taken's own comment.
    world = generate_world(seed=3, num_plates=4)
    del world.__dict__["steps_taken"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.steps_taken == 0


def test_loading_a_world_pickled_before_removed_points_log_existed_defaults_to_empty():
    # World.removed_points_log is a `default_factory=list` field -- unlike a plain-default
    # field (steps_taken), pickle restores __dict__ directly and never calls __init__, so an
    # old save's __dict__ has no such key at all until persistence backfills it.
    world = generate_world(seed=3, num_plates=4)
    del world.__dict__["removed_points_log"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.removed_points_log == []


def test_round_trip_preserves_stats_history():
    # World.stats_history is what makes a save/load round trip keep the Stats panel's history
    # charts (see World.record_stats/GET /world/stats_history) -- generation records one
    # baseline entry, each step should add (at most) one more.
    world = generate_world(seed=11, num_plates=6)
    assert len(world.stats_history) == 1
    step_world(world, 5_000_000)
    step_world(world, 5_000_000)
    assert len(world.stats_history) == 3
    assert [s["elapsed_years"] for s in world.stats_history] == [0.0, 5_000_000.0, 10_000_000.0]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.stats_history == world.stats_history


def test_loading_a_world_pickled_before_stats_history_existed_defaults_to_empty():
    world = generate_world(seed=3, num_plates=4)
    del world.__dict__["stats_history"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.stats_history == []


def test_loading_a_world_pickled_before_gap_tracks_existed_defaults_to_empty():
    world = generate_world(seed=3, num_plates=4)
    del world.__dict__["gap_tracks"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.gap_tracks == []


def test_loading_a_world_pickled_before_corner_notch_log_existed_defaults_to_empty():
    world = generate_world(seed=3, num_plates=4)
    del world.__dict__["corner_notch_log"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.corner_notch_log == []


def test_loading_a_world_pickled_before_pending_magma_parcels_existed_defaults_to_empty():
    world = generate_world(seed=3, num_plates=4)
    del world.__dict__["pending_magma_parcels"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.pending_magma_parcels == []


def test_loading_a_save_with_the_old_water_column_budget_rebuilds_it_in_m3():
    # Issue #257: saves before the budget was area-weighted carry `ocean_water_column_m`
    # (summed metres over nominal-area nodes). Loading drops it and re-snapshots the budget
    # as a volume from the save's own hypsometry, so sea level stays put.
    world = generate_world(seed=3, num_plates=4)
    budget = world.ocean_water_volume_m3
    del world.__dict__["ocean_water_volume_m3"]
    world.__dict__["ocean_water_column_m"] = 1.0e7

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert "ocean_water_column_m" not in loaded.__dict__
    assert loaded.ocean_water_volume_m3 == pytest.approx(budget, rel=1e-9)


def test_loading_a_world_whose_lines_predate_elev_change_reason_still_steps():
    # An ElevationLine pickled before the elev_change_reason OPTIONAL_FIELD existed has no
    # _elev_change_reason backing attr (pickle restores __dict__, never calls __init__).
    # ElevationLine.__getattr__ backfills it lazily as zeros so load + step still work.
    from app.world import step_world

    world = generate_world(seed=4, num_plates=6, surface="lines")
    for plate in world.plates:
        for line in plate.lines:
            line.__dict__.pop("_elev_change_reason", None)

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    for plate in loaded.plates:
        for line in plate.lines:
            assert np.all(line.elev_change_reason == 0.0)
    step_world(loaded, years=1_000_000)  # must not raise


def test_loading_garbage_bytes_raises_a_corrupt_save_error():
    with pytest.raises(persistence.CorruptSaveError, match="not a readable mantle-bloom save"):
        persistence.load_world_bytes(b"not a pickle at all")


def test_loading_a_truncated_save_raises_a_corrupt_save_error():
    data = persistence.save_world_bytes(generate_world(seed=3, num_plates=4))
    with pytest.raises(persistence.CorruptSaveError):
        persistence.load_world_bytes(data[: len(data) // 2])


def test_loading_a_pickle_of_the_wrong_type_raises():
    with pytest.raises(persistence.CorruptSaveError, match="expected a World"):
        persistence.load_world_bytes(pickle.dumps(42))


def test_loading_another_programs_envelope_raises():
    with pytest.raises(persistence.CorruptSaveError, match="not a mantle-bloom save"):
        persistence.load_world_bytes(pickle.dumps({"format": "something-else", "version": 1, "world": None}))


def test_saves_carry_a_format_version_envelope():
    world = generate_world(seed=3, num_plates=4)
    envelope = pickle.loads(persistence.save_world_bytes(world))

    assert envelope["format"] == persistence.SAVE_FORMAT
    assert envelope["version"] == persistence.SAVE_FORMAT_VERSION


def test_loading_a_version_1_bare_world_pickle_still_works():
    world = generate_world(seed=3, num_plates=4)

    loaded = persistence.load_world_bytes(pickle.dumps(world))
    assert loaded.seed == world.seed
    assert len(loaded.plates) == len(world.plates)


def test_loading_a_save_from_a_newer_format_version_raises():
    world = generate_world(seed=3, num_plates=4)
    data = pickle.dumps({"format": persistence.SAVE_FORMAT, "version": persistence.SAVE_FORMAT_VERSION + 1, "world": world})

    with pytest.raises(ValueError, match="format version"):
        persistence.load_world_bytes(data)


def test_loading_a_version_2_envelope_without_a_surface_still_works():
    world = generate_world(seed=3, num_plates=4)
    loaded = persistence.load_world_bytes(pickle.dumps({"format": persistence.SAVE_FORMAT, "version": 2, "world": world}))
    assert len(loaded.plates) == len(world.plates)


@pytest.mark.parametrize("version", [0, 1, True, "3", None, 3.0])
def test_loading_an_envelope_with_an_invalid_version_raises(version):
    world = generate_world(seed=3, num_plates=4)
    data = pickle.dumps({"format": persistence.SAVE_FORMAT, "version": version, "surface": "lines", "world": world})
    with pytest.raises(persistence.UnsupportedSaveVersionError):
        persistence.load_world_bytes(data)


def test_saves_declare_their_surface_and_the_loader_checks_it():
    lines = generate_world(seed=3, num_plates=4, surface="lines")
    quad = generate_world(seed=3, num_plates=4, surface="quad")
    assert pickle.loads(persistence.save_world_bytes(lines))["surface"] == "lines"
    assert pickle.loads(persistence.save_world_bytes(quad))["surface"] == "quad"

    lying = pickle.dumps({"format": persistence.SAVE_FORMAT, "version": 3, "surface": "quad", "world": lines})
    with pytest.raises(persistence.CorruptSaveError, match="declares a 'quad' world"):
        persistence.load_world_bytes(lying)
    unknown = pickle.dumps({"format": persistence.SAVE_FORMAT, "version": 3, "surface": "hexes", "world": lines})
    with pytest.raises(persistence.CorruptSaveError, match="unknown surface"):
        persistence.load_world_bytes(unknown)


def test_a_world_mixing_line_and_quad_plates_is_refused():
    lines = generate_world(seed=3, num_plates=4, surface="lines")
    quad = generate_world(seed=3, num_plates=4, surface="quad")
    lines.plates.append(quad.plates[0])
    with pytest.raises(persistence.CorruptSaveError, match="mixes"):
        persistence.save_world_bytes(lines)
    with pytest.raises(persistence.CorruptSaveError, match="mixes"):
        persistence.load_world_bytes(pickle.dumps({"format": persistence.SAVE_FORMAT, "version": 2, "world": lines}))


def test_a_quad_plate_from_a_newer_build_makes_the_save_unsupported(monkeypatch):
    from app import sparse_quad_patch

    world = generate_world(seed=3, num_plates=4, surface="quad")
    monkeypatch.setattr(sparse_quad_patch, "QUAD_SURFACE_FORMAT_VERSION", sparse_quad_patch.QUAD_SURFACE_FORMAT_VERSION + 1)
    data = persistence.save_world_bytes(world)
    monkeypatch.undo()
    with pytest.raises(persistence.UnsupportedSaveVersionError, match="sparse quad surface format version"):
        persistence.load_world_bytes(data)


def test_a_line_save_loads_as_lines_unless_conversion_is_asked_for():
    world = generate_world(seed=3, num_plates=4, surface="lines")
    data = pickle.dumps(world)  # a version 1 save, like every save written before #228

    legacy = persistence.load_world_bytes(data)
    assert persistence.world_surface(legacy) == "lines"
    assert legacy.surface_conversion is None

    converted = persistence.load_world_bytes(data, convert_lines=True)
    assert persistence.world_surface(converted) == "quad"
    assert converted.surface_conversion["from"] == "lines"
    assert converted.ocean_water_volume_m3 is not None
    assert converted.sea_level_m == legacy.sea_level_m
    # Asking to convert a world that is already quad is a no-op.
    again = persistence.load_world_bytes(persistence.save_world_bytes(converted), convert_lines=True)
    assert again.surface_conversion == converted.surface_conversion

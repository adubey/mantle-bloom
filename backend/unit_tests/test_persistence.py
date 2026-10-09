import pickle

import numpy as np
import pytest
from app import persistence
from app.elevation_lines import line_spacing_rad
from app.world import generate_world, step_world

from .legacy_lines import LegacyRowLookup, line_world_like, retired_class_names


def _line_world(seed=3, num_plates=4):
    """A line-backed world, as a save from before #251 holds."""
    world = generate_world(seed=seed, num_plates=num_plates, node_density=1.0)
    return line_world_like(world, line_spacing_rad(1.0))


def _line_pickle(payload) -> bytes:
    """`payload` pickled the way a build that still had the line classes pickled it."""
    with retired_class_names():
        return pickle.dumps(payload)


def test_round_trip_preserves_a_freshly_generated_world():
    world = generate_world(seed=42, num_plates=6)
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
        assert np.array_equal(loaded_plate.cell_keys, original_plate.cell_keys)
        assert np.array_equal(loaded_plate.collect("elevation"), original_plate.collect("elevation"))


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


def test_loading_a_world_pickled_before_pending_magma_parcels_existed_defaults_to_empty():
    world = generate_world(seed=3, num_plates=4)
    del world.__dict__["pending_magma_parcels"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.pending_magma_parcels == []


def test_loading_a_save_drops_the_retired_line_engine_state():
    # #251: saves from before the line surface was retired pickled the line engine's
    # corner-notch log and gap-fill choice; a loaded world, and so its next save, carries no
    # line state.
    world = generate_world(seed=3, num_plates=4)
    world.__dict__["corner_notch_log"] = [{"plate_id": 0, "outcome": "claimed"}]
    world.__dict__["gap_fill_algorithm"] = "frontier"

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert "corner_notch_log" not in loaded.__dict__
    assert "gap_fill_algorithm" not in loaded.__dict__


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


def test_loading_a_line_save_whose_lines_predate_elev_change_reason_still_steps():
    # A line pickled before the elev_change_reason field existed has no _elev_change_reason
    # backing attr (pickle restores __dict__, never calls __init__). Conversion reads it as
    # its registry default, so load + step still work.
    world = _line_world(seed=4, num_plates=6)
    for plate in world.plates:
        for line in plate.lines:
            line.__dict__.pop("_elev_change_reason", None)

    loaded = persistence.load_world_bytes(_line_pickle(world))
    for plate in loaded.plates:
        assert np.all(plate.collect("elev_change_reason") == 0.0)
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


def test_old_hydroclimate_pickle_module_paths_resolve_to_moved_classes():
    from app import legacy_conversion
    from app.hydroclimate.climate import ClimateFields

    unpickler = legacy_conversion.legacy_unpickler(b"")
    assert unpickler.find_class("app.climate", "ClimateFields") is ClimateFields


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
    lines = _line_world()
    quad = generate_world(seed=3, num_plates=4)
    assert pickle.loads(persistence.save_world_bytes(quad))["surface"] == "quad"

    lying = _line_pickle({"format": persistence.SAVE_FORMAT, "version": 3, "surface": "quad", "world": lines})
    with pytest.raises(persistence.CorruptSaveError, match="declares a 'quad' world"):
        persistence.load_world_bytes(lying)
    unknown = _line_pickle({"format": persistence.SAVE_FORMAT, "version": 3, "surface": "hexes", "world": lines})
    with pytest.raises(persistence.CorruptSaveError, match="unknown surface"):
        persistence.load_world_bytes(unknown)


def test_a_world_mixing_line_and_quad_plates_is_refused():
    lines = _line_world()
    quad = generate_world(seed=3, num_plates=4)
    lines.plates.append(quad.plates[0])
    with pytest.raises(persistence.CorruptSaveError, match="mixes"):
        persistence.save_world_bytes(lines)
    with pytest.raises(persistence.CorruptSaveError, match="mixes"):
        persistence.load_world_bytes(_line_pickle({"format": persistence.SAVE_FORMAT, "version": 2, "world": lines}))


def test_a_quad_plate_from_a_newer_build_makes_the_save_unsupported(monkeypatch):
    from app import sparse_quad_patch

    world = generate_world(seed=3, num_plates=4)
    monkeypatch.setattr(sparse_quad_patch, "QUAD_SURFACE_FORMAT_VERSION", sparse_quad_patch.QUAD_SURFACE_FORMAT_VERSION + 1)
    data = persistence.save_world_bytes(world)
    monkeypatch.undo()
    with pytest.raises(persistence.UnsupportedSaveVersionError, match="sparse quad surface format version"):
        persistence.load_world_bytes(data)


def test_a_line_save_converts_to_quads_on_load():
    world = _line_world()
    data = _line_pickle(world)  # a version 1 save, like every save written before #228

    converted = persistence.load_world_bytes(data)
    assert persistence.world_surface(converted) == "quad"
    assert converted.surface_conversion["from"] == "lines"
    assert converted.ocean_water_volume_m3 is not None
    assert converted.sea_level_m == world.sea_level_m
    # The conversion summary rides along through later saves; a quad save never converts.
    again = persistence.load_world_bytes(persistence.save_world_bytes(converted))
    assert again.surface_conversion == converted.surface_conversion
    assert persistence.load_world_bytes(persistence.save_world_bytes(generate_world(seed=3, num_plates=4))).surface_conversion is None


def test_a_line_save_holding_a_cached_row_lookup_converts_on_load():
    # A line plate whose containment fast path had run pickled its `_RowLookup` cache, a
    # class #251 deleted; 20 of the 51 real saves #248 inventoried hold one.
    world = _line_world()
    for plate in world.plates:
        plate._row_lookup_cache = LegacyRowLookup([line.phi for line in plate.lines])

    converted = persistence.load_world_bytes(_line_pickle(world))
    assert persistence.world_surface(converted) == "quad"
    assert len(converted.plates) == len(world.plates)


def test_a_line_save_the_converter_refuses_is_a_corrupt_save():
    world = _line_world()
    world.plates[0].lines[0].__dict__["_future_field"] = np.zeros(len(world.plates[0].lines[0]))
    with pytest.raises(persistence.CorruptSaveError, match="_future_field"):
        persistence.load_world_bytes(_line_pickle(world))

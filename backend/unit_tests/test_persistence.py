import pickle

import numpy as np
import pytest
from app import persistence
from app.world import generate_world, step_world


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
        assert len(loaded_plate.lines) == len(original_plate.lines)
        for original_line, loaded_line in zip(original_plate.lines, loaded_plate.lines):
            assert np.allclose(loaded_line.theta, original_line.theta)
            assert np.allclose(loaded_line.elevation, original_line.elevation)


def test_round_trip_preserves_state_only_a_step_would_populate():
    # collision_progress/volcanic_field_plate_ids/climate_cache/hydrology_cache/events are
    # all empty/None on a freshly generated world -- step it first so the round trip has to
    # actually carry a dict, a set, and the two cache dataclasses, not just empty defaults.
    world = generate_world(seed=7, num_plates=6)
    step_world(world, 5_000_000)
    step_world(world, 5_000_000)

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))

    assert loaded.elapsed_years == world.elapsed_years
    assert loaded.steps_taken == world.steps_taken == 2
    assert loaded.collision_progress == world.collision_progress
    assert loaded.volcanic_field_plate_ids == world.volcanic_field_plate_ids
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


def test_loading_garbage_bytes_raises():
    with pytest.raises(Exception):
        persistence.load_world_bytes(b"not a pickle at all")


def test_loading_a_pickle_of_the_wrong_type_raises():
    with pytest.raises(TypeError):
        persistence.load_world_bytes(pickle.dumps(42))

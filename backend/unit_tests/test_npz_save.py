import io
import json
import zipfile

import numpy as np
import pytest
from app import npz_save, persistence
from app.replay_log import fingerprint
from app.surface_parity import state_hash
from app.world import generate_world, step_world


@pytest.fixture(scope="module")
def stepped_world():
    """A small stepped world: caches populated, cheap enough to save once per test."""
    world = generate_world(seed=11, num_plates=2, node_density=0.5)
    step_world(world, 1_000_000)
    return world


@pytest.fixture(scope="module")
def stepped_bytes(stepped_world):
    return persistence.save_world_bytes(stepped_world)


def _manifest(data: bytes) -> tuple[dict, zipfile.ZipFile]:
    zf = zipfile.ZipFile(io.BytesIO(data))
    return json.loads(zf.read("manifest.json")), zf


def test_format_constants_match():
    assert npz_save.FORMAT == persistence.SAVE_FORMAT
    assert npz_save.VERSION == persistence.SAVE_FORMAT_VERSION == 5


def test_v5_save_is_a_zip_with_a_manifest(stepped_bytes):
    assert stepped_bytes[:4] == b"PK\x03\x04"
    manifest, _ = _manifest(stepped_bytes)
    assert manifest["format"] == persistence.SAVE_FORMAT
    assert manifest["version"] == 5
    assert manifest["surface"] == "quad"


def test_round_trip_keeps_authoritative_state(stepped_world, stepped_bytes):
    loaded = persistence.load_world_bytes(stepped_bytes)

    assert loaded is not stepped_world
    assert state_hash(loaded) == state_hash(stepped_world)
    assert loaded.seed == stepped_world.seed
    assert loaded.elapsed_years == stepped_world.elapsed_years
    assert loaded.steps_taken == stepped_world.steps_taken
    assert loaded.stats_history == stepped_world.stats_history
    assert loaded.events == stepped_world.events
    assert loaded.collision_progress == stepped_world.collision_progress


def test_save_keeps_step_read_caches_and_drops_erosion_and_scratch(stepped_world, stepped_bytes):
    assert stepped_world.climate_cache is not None and stepped_world.hydrology_cache is not None

    loaded = persistence.load_world_bytes(stepped_bytes)

    # The H9 gate requires a post-load step to land on identical authoritative state,
    # so every cache the next step reads before refreshing it survives the save; only
    # the erosion cache (overwritten whole by every erosion call) and the derived
    # scratch go.
    assert loaded.hydrology_cache_step == stepped_world.hydrology_cache_step
    np.testing.assert_array_equal(loaded.climate_cache.precipitation_mm, stepped_world.climate_cache.precipitation_mm)
    np.testing.assert_array_equal(loaded.hydrology_cache.is_ocean, stepped_world.hydrology_cache.is_ocean)
    assert loaded.atmosphere_cfd_state is not None
    assert loaded.erosion_cache is None
    assert loaded.collision_polarity_frame is not None
    assert loaded.land_kdtree_cache is None
    assert loaded.node_position_tree_cache is None


def test_post_load_step_lands_on_identical_state():
    # The H9 gate at unit scale: stepping the loaded copy matches stepping the world.
    world = generate_world(seed=11, num_plates=2, node_density=0.5)
    step_world(world, 1_000_000)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    step_world(world, 100_000)
    step_world(loaded, 100_000)

    assert state_hash(loaded) == state_hash(world)


def test_fresh_world_keeps_the_empty_cache_signal():
    world = generate_world(seed=3, num_plates=2, node_density=0.5)

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))

    assert loaded.climate_cache is None
    assert loaded.hydrology_cache is None
    assert state_hash(loaded) == state_hash(world)


def test_exact_profile_keeps_all_data_caches(stepped_world):
    data = persistence.save_exact_bytes(stepped_world)
    loaded = persistence.load_world_bytes(data)

    assert loaded.climate_cache is not None
    assert loaded.hydrology_cache is not None
    assert loaded.atmosphere_cfd_state is not None
    np.testing.assert_array_equal(loaded.climate_cache.precipitation_mm, stepped_world.climate_cache.precipitation_mm)
    np.testing.assert_array_equal(loaded.hydrology_cache.is_ocean, stepped_world.hydrology_cache.is_ocean)
    # Plate references re-link to the loaded world's own plates by id (#596).
    assert [p.plate_id for p in loaded.hydrology_cache.plates_in_order] == [p.plate_id for p in stepped_world.hydrology_cache.plates_in_order]
    assert all(any(p is q for q in loaded.plates) for p in loaded.hydrology_cache.plates_in_order)
    # Derived plate caches are populated in memory but reset on load (as pickles do
    # via `__getstate__`); their content is not authoritative state, so plates compare
    # by `state_hash` while every other kept attribute compares by fingerprint.
    assert state_hash(loaded) == state_hash(stepped_world)
    skip = set(npz_save._DERIVED_SCRATCH) | {"hydrology_cache", "plates"}
    kept = {k: v for k, v in fingerprint(stepped_world).items() if k not in skip}
    got = {k: v for k, v in fingerprint(loaded).items() if k not in skip}
    assert got == kept


def test_entry_count_stays_batched(stepped_bytes):
    _, zf = _manifest(stepped_bytes)
    # One entry per plate field plus one packed entry per dataclass-list column --
    # not one entry per array (see #594: ~9,000 entries unbatched).
    assert len(zf.namelist()) < 500


def test_v5_save_is_smaller_than_the_pickle_it_replaces(stepped_world):
    import pickle as _pickle

    legacy = _pickle.dumps(
        {
            "format": persistence.SAVE_FORMAT,
            "version": 4,
            "surface": persistence.world_surface(stepped_world),
            "world": stepped_world,
        }
    )
    assert len(persistence.save_world_bytes(stepped_world)) < len(legacy)


def test_unencodable_attribute_fails_closed():
    world = generate_world(seed=3, num_plates=2, node_density=0.5)
    world.__dict__["_probe_weird"] = object()

    with pytest.raises(TypeError):
        persistence.save_world_bytes(world)


def test_loading_a_truncated_v5_save_raises(stepped_bytes):
    with pytest.raises(persistence.CorruptSaveError):
        persistence.load_world_bytes(stepped_bytes[: len(stepped_bytes) // 2])


def _remanifested(stepped_bytes, mutate):
    manifest, zf = _manifest(stepped_bytes)
    mutate(manifest)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_STORED) as out:
        out.writestr("manifest.json", json.dumps(manifest))
        for name in zf.namelist():
            if name != "manifest.json":
                out.writestr(name, zf.read(name))
    return buf.getvalue()


def test_loading_a_save_missing_its_climate_cache_rebuilds_it(stepped_world, stepped_bytes):
    # Safety net for saves that predate the kept cache: the loader recomputes a missing
    # climate cache for a stepped world instead of leaving the zeros-transpiration path.
    def drop_climate(manifest):
        manifest["world"].pop("climate_cache")

    data = _remanifested(stepped_bytes, drop_climate)
    loaded = persistence.load_world_bytes(data)

    assert loaded.climate_cache is not None
    assert np.any(loaded.climate_cache.precipitation_mm > 0)


def test_loading_a_newer_v5_raises():
    world = generate_world(seed=3, num_plates=2, node_density=0.5)
    data = persistence.save_world_bytes(world)

    def bump(manifest):
        manifest["version"] = 99

    with pytest.raises(persistence.UnsupportedSaveVersionError, match="format version"):
        persistence.load_world_bytes(_remanifested(data, bump))


def test_loading_a_v5_with_a_missing_array_raises(stepped_bytes):
    manifest, zf = _manifest(stepped_bytes)
    victim = next(n for n in zf.namelist() if n.endswith(".npy"))
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_STORED) as out:
        out.writestr("manifest.json", json.dumps(manifest))
        for name in zf.namelist():
            if name not in ("manifest.json", victim):
                out.writestr(name, zf.read(name))

    with pytest.raises(persistence.CorruptSaveError, match="missing array"):
        persistence.load_world_bytes(buf.getvalue())


def test_loading_a_v5_claiming_an_unknown_surface_raises(stepped_bytes):
    def lie(manifest):
        manifest["surface"] = "hexes"

    with pytest.raises(persistence.CorruptSaveError, match="unknown surface"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, lie))


def test_a_quad_plate_from_a_newer_build_makes_the_save_unsupported(monkeypatch):
    from app import sparse_quad_patch

    world = generate_world(seed=3, num_plates=2, node_density=0.5)
    monkeypatch.setattr(sparse_quad_patch, "QUAD_SURFACE_FORMAT_VERSION", sparse_quad_patch.QUAD_SURFACE_FORMAT_VERSION + 1)
    data = persistence.save_world_bytes(world)
    monkeypatch.undo()
    with pytest.raises(persistence.UnsupportedSaveVersionError, match="sparse quad surface format version"):
        persistence.load_world_bytes(data)

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


def _packed_column_bytes():
    """A packed (empty, missing, present) column plus zip bytes holding its entry."""
    arrays = {}
    node = npz_save._pack_arrays("f", "p", [np.zeros((0, 3)), None, np.zeros((2, 2))], arrays)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_STORED) as out:
        for name, arr in arrays.items():
            entry = io.BytesIO()
            np.save(entry, arr, allow_pickle=False)
            out.writestr(f"a/{name}.npy", entry.getvalue())
    return node, buf.getvalue()


def _thaw_packed(node, data):
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        return npz_save._thaw_packed(node, zf)


def test_packed_column_keeps_empty_arrays_distinct_from_missing():
    # An empty array and a missing one both have length 0: only the null bitmap
    # tells them apart (#596 review).
    node, data = _packed_column_bytes()

    back = _thaw_packed(node, data)

    assert isinstance(back[0], np.ndarray) and back[0].shape == (0, 3)
    assert back[1] is None
    assert back[2].shape == (2, 2)


def test_thaw_packed_rejects_a_retagged_dtype():
    node, data = _packed_column_bytes()
    node = dict(node, dtype="int32")

    with pytest.raises(persistence.CorruptSaveError, match="changed dtype"):
        _thaw_packed(node, data)


def test_thaw_packed_rejects_inconsistent_offsets():
    # The last item is shrunk with its shape, so every reshape fits but the cursor
    # stops short of the entry.
    node, data = _packed_column_bytes()
    node = dict(node, lengths=[0, 0, 2], shapes=[[0, 3], [], [1, 2]])

    with pytest.raises(persistence.CorruptSaveError, match="inconsistent offsets"):
        _thaw_packed(node, data)


def test_thaw_packed_rejects_an_item_that_does_not_fit_its_shape():
    node, data = _packed_column_bytes()
    node = dict(node, lengths=[1, 0, 4])

    with pytest.raises(persistence.CorruptSaveError, match="does not fit shape"):
        _thaw_packed(node, data)


def _first_with_key(node, key):
    if isinstance(node, dict):
        if key in node:
            return node
        for v in node.values():
            hit = _first_with_key(v, key)
            if hit is not None:
                return hit
    elif isinstance(node, list):
        for v in node:
            hit = _first_with_key(v, key)
            if hit is not None:
                return hit
    return None


def test_loading_a_v5_with_an_unknown_tag_raises(stepped_bytes):
    def retag(manifest):
        node = _first_with_key(manifest, "__npy__")
        node["__unwritten_future_tag__"] = node.pop("__npy__")

    with pytest.raises(persistence.CorruptSaveError, match="unknown tag"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, retag))


def test_loading_a_v5_with_a_bad_healpix_nside_raises(stepped_bytes):
    def spoil(manifest):
        _first_with_key(manifest, "__healpix_grid__")["__healpix_grid__"] = -4

    with pytest.raises(persistence.CorruptSaveError, match="bad HEALPix nside"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, spoil))


def test_loading_a_v5_with_duplicate_plate_ids_raises(stepped_bytes):
    def duplicate(manifest):
        manifest["plates"].append(manifest["plates"][0])

    with pytest.raises(persistence.CorruptSaveError, match="two plates with one id"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, duplicate))


def test_bytes_attributes_fail_with_a_path(stepped_world):
    stepped_world.__dict__["_bytes_probe"] = b"nope"
    try:
        with pytest.raises(TypeError, match="_bytes_probe"):
            persistence.save_world_bytes(stepped_world)
    finally:
        del stepped_world.__dict__["_bytes_probe"]


def test_numpy_bytes_scalars_fail_with_a_path(stepped_world):
    stepped_world.__dict__["_bytes_probe"] = np.bytes_(b"nope")
    try:
        with pytest.raises(TypeError, match="_bytes_probe"):
            persistence.save_world_bytes(stepped_world)
    finally:
        del stepped_world.__dict__["_bytes_probe"]


def test_exact_profile_step_lands_on_identical_state():
    # The H9 gate at unit scale for the exact profile: stepping the loaded copy
    # matches stepping the world (the approximate profile has its own pin above).
    world = generate_world(seed=11, num_plates=2, node_density=0.5)
    step_world(world, 1_000_000)
    loaded = persistence.load_world_bytes(persistence.save_exact_bytes(world))
    step_world(world, 100_000)
    step_world(loaded, 100_000)

    assert state_hash(loaded) == state_hash(world)


def test_v5_load_backfills_a_field_added_after_the_save(stepped_bytes):
    manifest, _ = _manifest(stepped_bytes)
    assert "faults" in manifest["world"]

    data = _remanifested(stepped_bytes, lambda m: m["world"].pop("faults"))

    assert persistence.load_world_bytes(data).faults == []


def test_compact_profile_drops_step_caches_and_rebuilds_climate_on_load(stepped_world):
    assert stepped_world.atmosphere_cfd_state is not None
    assert stepped_world.collision_polarity_frame is not None

    loaded = persistence.load_world_bytes(persistence.save_compact_bytes(stepped_world))

    assert loaded.erosion_cache is None
    assert loaded.hydrology_cache is None
    assert loaded.atmosphere_cfd_state is None
    assert loaded.collision_polarity_frame is None
    # The loader rebuilds the climate cache for a stepped world (the safety net in
    # `_rebuild_climate`); the other three rebuild on the first post-load step.
    assert loaded.climate_cache is not None
    assert np.any(loaded.climate_cache.precipitation_mm > 0)


def _summarize(world):
    hydro = world.hydrology_cache
    return {
        "land_frac": 1.0 - float(hydro.is_ocean.mean()),
        "sea_level_m": float(world.sea_level_m),
    }


def test_compact_profile_step_lands_within_tolerance():
    # The approximate continuation gate: a compact load rebuilds its step-read caches
    # as it goes, landing close to but not on the uninterrupted state. Tolerances are
    # ~5x the land-fraction and ~4x the sea-level divergence measured on the reference
    # world (seed 42, density 1.0, 1 Myr + 100 kyr: 2.1e-4, 0.26 m) -- a tripwire for
    # the envelope, not an exactness proof (that is H9, on the default/exact profiles).
    world = generate_world(seed=42, node_density=1.0)
    step_world(world, 1_000_000)
    loaded = persistence.load_world_bytes(persistence.save_compact_bytes(world))
    step_world(world, 100_000)
    step_world(loaded, 100_000)

    assert state_hash(loaded) != state_hash(world)
    ref, got = _summarize(world), _summarize(loaded)
    assert abs(got["land_frac"] - ref["land_frac"]) <= 1e-3
    assert abs(got["sea_level_m"] - ref["sea_level_m"]) <= 1.0


def test_compact_save_is_smaller_than_default(stepped_world, stepped_bytes):
    assert len(persistence.save_compact_bytes(stepped_world)) < len(stepped_bytes)


def test_compact_cfd_mode_load_steps_without_crashing():
    # A compact save drops the CFD state; in cfd wind mode the first post-load step
    # must cold-start it, not dereference None in `refresh_forcing` (#596 review).
    # Crash-regression only (the envelope pin above is diagnostic-mode): completes
    # and holds sane summary state.
    world = generate_world(seed=3, num_plates=2, node_density=0.5)
    world.wind_model = "cfd"
    step_world(world, 100_000)
    loaded = persistence.load_world_bytes(persistence.save_compact_bytes(world))

    assert loaded.atmosphere_cfd_state is None
    step_world(loaded, 100_000)

    assert loaded.atmosphere_cfd_state is not None
    assert 0.0 <= 1.0 - float(loaded.hydrology_cache.is_ocean.mean()) <= 1.0
    assert np.isfinite(loaded.sea_level_m)


def test_resolve_rejects_non_dataclass():
    with pytest.raises(persistence.CorruptSaveError, match="non-dataclass"):
        npz_save._resolve("builtins.dict")
    with pytest.raises(persistence.CorruptSaveError, match="non-dataclass"):
        npz_save._resolve("os.makedirs")


def test_loading_a_v5_naming_a_non_dataclass_raises(stepped_bytes):
    def retag(manifest):
        _first_with_key(manifest, "__dc__")["__dc__"] = "builtins.dict"

    with pytest.raises(persistence.CorruptSaveError, match="non-dataclass"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, retag))


def test_loading_a_v5_with_inconsistent_arrays_raises(stepped_bytes, monkeypatch):
    from app.hydroclimate import climate

    def drop_climate(manifest):
        manifest["world"]["climate_cache"] = None

    def boom(*args, **kwargs):
        raise ValueError("inconsistent test grid")

    monkeypatch.setattr(climate, "compute_climate", boom)

    with pytest.raises(persistence.CorruptSaveError, match="inconsistent arrays"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, drop_climate))


def test_loading_a_v5_with_a_wrong_typed_climate_cache_raises(stepped_bytes):
    def spoil(manifest):
        manifest["world"]["climate_cache"] = [1, 2]

    with pytest.raises(persistence.CorruptSaveError, match="not climate fields"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, spoil))


def test_loading_a_v5_plate_missing_its_surface_version_raises(stepped_bytes):
    def strip(manifest):
        del manifest["plates"][0]["_plate_surface_version"]

    with pytest.raises(persistence.CorruptSaveError, match="surface version"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, strip))


def test_loading_a_v5_missing_world_scalars_raises(stepped_bytes):
    def strip(manifest):
        manifest["world"].pop("seed")

    with pytest.raises(persistence.CorruptSaveError, match="missing world fields"):
        persistence.load_world_bytes(_remanifested(stepped_bytes, strip))

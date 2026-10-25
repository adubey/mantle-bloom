import json
import os

import pytest
from app import checkpoint, persistence
from app.checkpoint import CheckpointStore, RunStatus
from app.world import generate_world


@pytest.fixture(scope="module")
def world():
    return generate_world(seed=5, num_plates=2, node_density=0.5)


def test_world_and_status_round_trip(tmp_path, world):
    store = CheckpointStore(tmp_path)
    store.save_world(world)
    store.save_status(RunStatus(running=True, step_years=250_000.0))

    loaded, status = CheckpointStore(tmp_path).load()

    assert loaded is not world
    assert (loaded.seed, loaded.elapsed_years, len(loaded.plates)) == (world.seed, world.elapsed_years, len(world.plates))
    assert status == RunStatus(running=True, step_years=250_000.0)
    # An ordinary save, so the owner can download it and open it like any other.
    assert persistence.load_world_bytes(store.world_path.read_bytes()).seed == world.seed


def test_no_checkpoint_starts_empty_and_stopped(tmp_path, caplog):
    assert CheckpointStore(tmp_path / "new").load() == (None, RunStatus())
    assert "no saved world" in caplog.text


def test_unreadable_world_is_set_aside_and_a_running_status_stays_stopped(tmp_path, caplog):
    store = CheckpointStore(tmp_path)
    store.world_path.write_bytes(b"not a save")
    store.save_status(RunStatus(running=True, step_years=1_000.0))

    loaded, status = store.load()

    assert loaded is None
    assert status == RunStatus(running=False, step_years=1_000.0)
    assert store.load()[1] == status
    assert not store.world_path.exists()
    assert [p.read_bytes() for p in tmp_path.glob("world.mbworld.invalid-*")] == [b"not a save"]
    assert "can't load" in caplog.text and "staying stopped" in caplog.text


@pytest.mark.parametrize(
    "contents",
    ["{", json.dumps({"version": 1, "running": "yes", "step_years": 1.0}),
     json.dumps({"version": 1, "running": True, "step_years": -5}),
     json.dumps({"version": 99, "running": True, "step_years": 1.0})],
)
def test_unreadable_status_stays_stopped(tmp_path, world, caplog, contents):
    store = CheckpointStore(tmp_path)
    store.save_world(world)
    store.status_path.write_text(contents)

    loaded, status = store.load()

    assert loaded is not None
    assert status == RunStatus()
    assert "can't read" in caplog.text


def test_a_write_cut_off_partway_leaves_the_previous_checkpoint(tmp_path, world, monkeypatch):
    store = CheckpointStore(tmp_path)
    store.save_world(world)
    before = store.world_path.read_bytes()
    world.elapsed_years += 1.0

    def dies(fd):
        raise OSError("disk full")

    monkeypatch.setattr(os, "fsync", dies)
    with pytest.raises(OSError):
        store.save_world(world)
    monkeypatch.undo()
    world.elapsed_years -= 1.0

    assert store.world_path.read_bytes() == before
    assert sorted(p.name for p in tmp_path.iterdir()) == ["world.mbworld"]


def test_writes_left_by_a_killed_process_are_removed_at_startup(tmp_path, world, caplog):
    store = CheckpointStore(tmp_path)
    store.save_world(world)
    (tmp_path / ".checkpoint-abc123").write_bytes(b"half a world")

    loaded, _ = store.load()

    assert loaded.seed == world.seed
    assert sorted(p.name for p in tmp_path.iterdir()) == ["world.mbworld"]
    assert "cut off" in caplog.text


@pytest.mark.skipif(os.geteuid() == 0, reason="root can write to a read-only directory")
def test_an_unwritable_directory_fails_at_startup(tmp_path):
    # Rather than at the first checkpoint, after the deploy's health check has passed (#376).
    tmp_path.chmod(0o500)
    try:
        with pytest.raises(PermissionError):
            CheckpointStore(tmp_path).load()
    finally:
        tmp_path.chmod(0o700)


def test_registry_index_round_trip_and_per_world_stores(tmp_path, world):
    root = CheckpointStore(tmp_path)
    assert root.load_index() is None  # no index yet: a legacy directory
    records = [
        {"world_id": "a" * 32, "name": "A", "is_public": True, "created_at": 1.0, "updated_at": 2.0},
        {"world_id": "b" * 32, "name": "B", "is_public": False, "created_at": 3.0, "updated_at": 4.0},
    ]
    root.save_index(records)
    assert root.load_index() == records
    # Each world's checkpoint lives in its own subdirectory.
    root.for_world("a" * 32).save_world(world)
    root.for_world("a" * 32).save_status(RunStatus(running=True, step_years=10.0))
    loaded, status = root.for_world("a" * 32).load()
    assert (loaded.seed, status.running) == (world.seed, True)
    assert root.for_world("b" * 32).load() == (None, RunStatus())


def test_unreadable_registry_index_raises(tmp_path):
    (tmp_path / checkpoint.WORLDS_INDEX_FILE).write_text("not json")
    with pytest.raises(checkpoint.RegistryIndexError):
        CheckpointStore(tmp_path).load_index()


def test_suspicious_world_ids_fail_closed(tmp_path):
    # A tampered worlds.json must not turn for_world/rmtree into a path
    # traversal: only uuid4-hex IDs reach per-world directories, and a bad ID in
    # the index fails closed like any other corrupt index.
    root = CheckpointStore(tmp_path)
    with pytest.raises(checkpoint.RegistryIndexError):
        root.for_world("../../escape")
    with pytest.raises(checkpoint.RegistryIndexError):
        root.for_world("A" * 32)
    root.save_index([
        {"world_id": "../tampered", "name": "X", "is_public": False, "created_at": 0.0, "updated_at": 0.0}
    ])
    with pytest.raises(checkpoint.RegistryIndexError):
        root.load_index()


@pytest.mark.parametrize("field,value", [
    ("name", 42),
    ("is_public", "yes"),
    ("created_at", "garbage"),
    ("updated_at", True),
])
def test_garbage_metadata_types_fail_closed(tmp_path, field, value):
    # Timestamps are read with float() at restore outside any try: a garbage one
    # that passed validation would fail startup for every world. Reject it with
    # the index instead, like any other corruption.
    root = CheckpointStore(tmp_path)
    record = {"world_id": "a" * 32, "name": "A", "is_public": True, "created_at": 1.0, "updated_at": 2.0}
    record[field] = value
    root.save_index([record])
    with pytest.raises(checkpoint.RegistryIndexError):
        root.load_index()


def test_directory_comes_from_the_env_or_the_railway_volume(monkeypatch, tmp_path):
    monkeypatch.delenv(checkpoint.CHECKPOINT_DIR_ENV, raising=False)
    monkeypatch.delenv(checkpoint.RAILWAY_VOLUME_ENV, raising=False)
    assert CheckpointStore.from_env() is None
    monkeypatch.setenv(checkpoint.RAILWAY_VOLUME_ENV, str(tmp_path / "volume"))
    assert CheckpointStore.from_env().directory == tmp_path / "volume"
    monkeypatch.setenv(checkpoint.CHECKPOINT_DIR_ENV, str(tmp_path / "explicit"))
    assert CheckpointStore.from_env().directory == tmp_path / "explicit"

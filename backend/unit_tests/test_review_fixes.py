"""Regression tests for the #563 review findings (S1, S3-S5).

Each test targets one fixed hazard directly instead of re-driving the whole
suite: the guest-step lock must never park while holding the compute gate, a
corrupt per-world checkpoint must surface as `load_error`, a corrupt registry
index must detach checkpoints instead of being overwritten, delete must refuse
while a guest step is in flight, and streams queued at the gate must bail out
when their world is deleted before their turn arrives.
"""

import threading

import pytest
from fastapi.testclient import TestClient

from app import main
from app.main import app


@pytest.fixture
def client(tmp_path, monkeypatch):
    """A fresh registry with checkpointing pointed at a scratch directory.

    Snapshots the process-global checkpoint state so these tests don't leak it
    into (or inherit it from) their neighbors.
    """
    monkeypatch.setenv("MANTLE_BLOOM_CHECKPOINT_DIR", str(tmp_path))
    old_root, old_error = main._checkpoint_root, main._checkpoint_error
    for session in main._registry.list():
        if session.runner is not None:
            session.runner.close()
        main._registry.remove(session.world_id)
    main._restore_checkpoint()
    yield TestClient(app, base_url="http://localhost", client=("127.0.0.1", 50000))
    for session in main._registry.list():
        if session.runner is not None:
            session.runner.close()
        main._registry.remove(session.world_id)
    main._checkpoint_root, main._checkpoint_error = old_root, old_error


def _tracked_session(world_id: str):
    session = main._registry.find(world_id)
    if session.runner is not None:
        session.runner.close()
        session.runner = None
    return session


def test_guest_lock_never_parks_while_holding_the_gate(client):
    """S1: acquiring the session lock through _SessionComputeLock must be
    non-blocking -- parking there would hold the compute gate (stalling every
    other world) inside the coordinator's mutex."""
    record = client.post("/worlds", json={}).json()
    session = _tracked_session(record["world_id"])
    lock = main._SessionComputeLock(session)
    session.lock.acquire()
    try:
        assert lock.acquire(blocking=False) is False
        # The failed acquire released the gate: a fresh take succeeds.
        assert main._compute_gate.acquire(timeout=0)
        main._compute_gate.release()
    finally:
        session.lock.release()


def test_delete_refuses_while_a_guest_step_is_in_flight(client):
    """S5: a guest step sets no operation flag, so delete must check the
    coordinator instead of detaching the session from underneath it."""
    record = client.post("/worlds", json={}).json()
    session = _tracked_session(record["world_id"])
    session.guest_steps._inflight = threading.Event()
    try:
        resp = client.delete(f"/worlds/{record['world_id']}")
        assert resp.status_code == 409
    finally:
        session.guest_steps._inflight = None
    assert client.delete(f"/worlds/{record['world_id']}").status_code == 200


def test_generate_debug_refuses_while_an_operation_is_in_flight(client):
    """Review follow-up: generate_debug used to skip the operation flag, so a
    debug generation could swap session.world out from under a step queued at
    the compute gate -- the step would mutate a detached world whose work never
    reached the published state. Now it 503s like every other replacer."""
    record = client.post("/worlds", json={}).json()
    session = _tracked_session(record["world_id"])
    scenario = client.get("/debug/scenarios").json()["scenarios"][0]["name"]
    assert session.begin_operation()
    try:
        resp = client.post(f"/worlds/{record['world_id']}/generate_debug", json={"scenario": scenario, "seed": 1})
        assert resp.status_code == 503
    finally:
        session.end_operation()
    assert client.post(f"/worlds/{record['world_id']}/generate_debug", json={"scenario": scenario, "seed": 1}).status_code == 200


def test_generate_on_deleted_world_bails_without_publishing(client):
    """S5: a stream that queued at the gate before its world was deleted must
    bail out instead of generating into the detached session."""
    record = client.post("/worlds", json={}).json()
    session = _tracked_session(record["world_id"])
    session.deleted = True
    try:
        resp = client.post(f"/worlds/{record['world_id']}/generate", json={})
        assert resp.status_code == 200
        assert resp.text == ""
        assert session.world is None
        assert session.feed.summary() is None
    finally:
        session.deleted = False
        client.delete(f"/worlds/{record['world_id']}")


def test_corrupt_world_checkpoint_surfaces_as_load_error(client, tmp_path):
    """S3: a world file that can't be loaded must show up as `load_error` on
    the record, not as a silent empty world."""
    record = client.post("/worlds", json={"name": "doomed"}).json()
    session = _tracked_session(record["world_id"])
    session.checkpoint_store.directory.mkdir(parents=True, exist_ok=True)
    session.checkpoint_store.world_path.write_bytes(b"not a save")
    main._restore_checkpoint()
    assert session.load_error is not None
    body = client.get(f"/worlds/{record['world_id']}").json()
    assert body["load_error"] == session.load_error
    assert list(tmp_path.glob(f"worlds/{record['world_id']}/*.invalid-*")) != []


def test_corrupt_index_detaches_checkpoints_and_survives(client, tmp_path):
    """S4: a corrupt worlds.json must fail closed (healthz reports it) without
    later creates overwriting the evidence."""
    (tmp_path / "worlds.json").write_text("{")
    main._restore_checkpoint()
    assert main._checkpoint_root is None
    assert main._checkpoint_error
    assert client.post("/worlds", json={}).status_code == 201
    assert (tmp_path / "worlds.json").read_text() == "{"

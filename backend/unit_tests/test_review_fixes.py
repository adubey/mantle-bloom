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
        # F1: the early bail must still release the gate, the session lock and
        # the operation flag -- a pre-try return skipped the finally and wedged
        # every world at the leaked gate.
        assert main._compute_gate.acquire(timeout=0)
        main._compute_gate.release()
        assert session.lock.acquire(blocking=False)
        session.lock.release()
        assert not session.operation_busy()
    finally:
        session.deleted = False
        client.delete(f"/worlds/{record['world_id']}")


def test_animate_on_deleted_world_bails_without_leaking(client):
    """F1/F2: same early-bail guarantee for the animation stream, which also
    releases the gate between frames."""
    record = client.post("/worlds", json={}).json()
    scenario = client.get("/debug/scenarios").json()["scenarios"][0]["name"]
    assert client.post(
        f"/worlds/{record['world_id']}/generate_debug", json={"scenario": scenario, "seed": 1}
    ).status_code == 200
    session = _tracked_session(record["world_id"])
    session.deleted = True
    try:
        resp = client.post(
            f"/worlds/{record['world_id']}/animate",
            json={"step_years": 100000, "steps_per_frame": 1, "num_frames": 1},
        )
        assert resp.status_code == 200
        assert main._compute_gate.acquire(timeout=0)
        main._compute_gate.release()
        assert session.lock.acquire(blocking=False)
        session.lock.release()
        assert not session.operation_busy()
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


def test_load_refuses_while_the_session_lock_is_held(client):
    """F3: load is async, so a held session lock must 503 instead of freezing
    the event loop for the whole step holding it."""
    record = client.post("/worlds", json={}).json()
    scenario = client.get("/debug/scenarios").json()["scenarios"][0]["name"]
    assert client.post(
        f"/worlds/{record['world_id']}/generate_debug", json={"scenario": scenario, "seed": 1}
    ).status_code == 200
    saved = client.get(f"/worlds/{record['world_id']}/save").content
    session = _tracked_session(record["world_id"])
    session.lock.acquire()
    try:
        resp = client.post(f"/worlds/{record['world_id']}/load", content=saved)
        assert resp.status_code == 503
    finally:
        session.lock.release()
    assert client.post(f"/worlds/{record['world_id']}/load", content=saved).status_code == 200


def test_world_change_persists_updated_at(client, tmp_path):
    """F7: the frontend auto-selects the most recently updated world, so a
    change must persist updated_at -- not just touch it in memory."""
    import json

    def updated_at_of(world_id: str) -> float:
        index = json.loads((tmp_path / "worlds.json").read_text())
        return next(w["updated_at"] for w in index["worlds"] if w["world_id"] == world_id)

    record = client.post("/worlds", json={}).json()
    before = updated_at_of(record["world_id"])
    scenario = client.get("/debug/scenarios").json()["scenarios"][0]["name"]
    assert client.post(
        f"/worlds/{record['world_id']}/generate_debug", json={"scenario": scenario, "seed": 1}
    ).status_code == 200
    assert updated_at_of(record["world_id"]) > before


def test_failing_index_save_is_logged_not_raised(client, tmp_path):
    """Advisory 2: _save_index runs before the lock/gate release in stream
    finally blocks, so a failing write must never raise and wedge every world.
    Record ops keep working (201) while the volume is unwritable."""
    import os

    from app import checkpoint

    store = checkpoint.CheckpointStore(tmp_path / "ro")
    store.directory.mkdir(parents=True)
    os.chmod(store.directory, 0o555)
    # Plain assignment, not monkeypatch.setattr: the client fixture depends on
    # monkeypatch, so monkeypatch undoes AFTER the fixture teardown and would
    # resurrect the store for later tests. A test-level finally runs first.
    old_root = main._checkpoint_root
    main._checkpoint_root = store
    try:
        main._save_index()
        assert client.post("/worlds", json={}).status_code == 201
    finally:
        main._checkpoint_root = old_root
        os.chmod(store.directory, 0o755)


def test_concurrent_index_saves_converge_on_memory(client, tmp_path):
    """Round-4 Finding 1: concurrent saves from different worlds must not lose
    each other's updates. Every mutation is followed by its own save and saves
    are serialized, so the last completed save necessarily read after every
    completed mutation: the file must equal live memory once all threads join."""
    import json
    import threading

    sessions = [main._registry.create(name=f"w{i}") for i in range(3)]
    for session in sessions:
        main._wire_session(session)

    def hammer(session, tag: str) -> None:
        for i in range(20):
            with session.lock:
                session.name = f"{tag}-{i}"
                session.touch()
            main._save_index()

    threads = [threading.Thread(target=hammer, args=(s, f"t{n}")) for n, s in enumerate(sessions)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=60)
    assert not any(thread.is_alive() for thread in threads)
    index = {
        record["world_id"]: record
        for record in json.loads((tmp_path / "worlds.json").read_text())["worlds"]
    }
    for session in sessions:
        assert index[session.world_id] == session.record()


def test_delete_wins_over_concurrent_index_saves(client, tmp_path):
    """Round-5 P2-index: a delete's index write must not be resurrected by a
    concurrent step's stale snapshot. Saves are serialized under _index_lock, so
    any save that read before the removal completes before the delete's save --
    the final file never names the deleted world."""
    import json
    import threading

    record = client.post("/worlds", json={"name": "doomed"}).json()
    wid = record["world_id"]
    stop = threading.Event()

    def saver() -> None:
        while not stop.is_set():
            main._save_index()

    threads = [threading.Thread(target=saver) for _ in range(4)]
    for thread in threads:
        thread.start()
    try:
        assert client.delete(f"/worlds/{wid}").status_code == 200
    finally:
        stop.set()
    for thread in threads:
        thread.join(timeout=60)
    assert not any(thread.is_alive() for thread in threads)
    index = json.loads((tmp_path / "worlds.json").read_text())["worlds"]
    assert all(w["world_id"] != wid for w in index)


def test_guest_step_waiters_are_bounded(client):
    """Round-5 P1-guest: only GUEST_MAX_WAITERS requests may queue behind a
    running guest step; the rest fail fast instead of pinning worker threads."""
    import threading
    import time

    from app import public_view

    record = client.post("/worlds", json={}).json()
    session = _tracked_session(record["world_id"])
    stuck = threading.Event()
    session.guest_steps._inflight = stuck  # a step is computing; never finishes on its own

    def hold_a_slot() -> None:
        session.guest_steps.request_step(
            "rev", True, lambda: "rev", threading.Lock(), lambda: None, inflight_wait_s=60.0
        )

    holders = [threading.Thread(target=hold_a_slot, daemon=True) for _ in range(public_view.GUEST_MAX_WAITERS)]
    for holder in holders:
        holder.start()
    try:
        # Wait until every slot is occupied (poll until a non-blocking take fails,
        # releasing takes that succeed).
        deadline = time.monotonic() + 10
        while session.guest_steps._waiter_slots.acquire(blocking=False):
            session.guest_steps._waiter_slots.release()
            assert time.monotonic() < deadline, "waiters never took all slots"
            time.sleep(0.01)
        outcome = session.guest_steps.request_step(
            "rev", True, lambda: "rev", threading.Lock(), lambda: None, inflight_wait_s=0
        )
        assert outcome.status == "busy"
    finally:
        # Wake the holders and clear the inflight marker so they become computers
        # one by one, finish, and exit cleanly instead of spinning.
        stuck.set()
        session.guest_steps._inflight = None
        for holder in holders:
            holder.join(timeout=10)
            assert not holder.is_alive()


def test_restore_skips_worlds_over_the_cap(client, tmp_path):
    """Round-5 P1-memory: an index holding more worlds than the cap allows must
    not load them all. The most recently updated up to the cap restore; the rest
    are adopted as visible, deletable records with load_error, files preserved."""
    from app import checkpoint as checkpoint_mod
    from app.world_registry import WorldRegistry

    checkpoint_mod.CheckpointStore(tmp_path).save_index([
        {"world_id": "a" * 32, "name": "Old", "is_public": False, "created_at": 1.0, "updated_at": 1.0},
        {"world_id": "b" * 32, "name": "New", "is_public": False, "created_at": 2.0, "updated_at": 2.0},
    ])
    real_registry = main._registry
    main._registry = WorldRegistry(max_worlds=1)
    try:
        main._restore_checkpoint()
        sessions = {s.world_id: s for s in main._registry.list()}
        assert set(sessions) == {"a" * 32, "b" * 32}
        assert sessions["b" * 32].load_error is None
        assert "MANTLE_BLOOM_MAX_WORLDS" in (sessions["a" * 32].load_error or "")
    finally:
        main._registry = real_registry


def test_animation_rejects_steps_per_frame_over_dialog_max(client):
    """Round-5 P2-anim: one frame may hold the gate for at most the dialog's
    largest choice (100 steps)."""
    record = client.post("/worlds", json={}).json()
    scenario = client.get("/debug/scenarios").json()["scenarios"][0]["name"]
    client.post(f"/worlds/{record['world_id']}/generate_debug", json={"scenario": scenario, "seed": 1})
    wid = record["world_id"]
    assert client.post(
        f"/worlds/{wid}/animate", json={"step_years": 100000, "steps_per_frame": 101, "num_frames": 2}
    ).status_code == 400
    assert client.post(
        f"/worlds/{wid}/animate", json={"step_years": 100000, "steps_per_frame": 100, "num_frames": 2}
    ).status_code == 200


def test_blank_or_huge_world_names_are_rejected(client):
    """F9: an empty name renders as a blank option (and restores differently);
    a huge one is rewritten to worlds.json on every save."""
    assert client.post("/worlds", json={"name": ""}).status_code == 422
    assert client.post("/worlds", json={"name": "x" * 101}).status_code == 422
    record = client.post("/worlds", json={"name": "ok"}).json()
    assert client.patch(f"/worlds/{record['world_id']}", json={"name": ""}).status_code == 422
    assert client.patch(f"/worlds/{record['world_id']}", json={"name": "y" * 101}).status_code == 422
    assert client.patch(f"/worlds/{record['world_id']}", json={"name": "fine"}).json()["name"] == "fine"

"""Exercise the snapshot API and idle lifecycle without starting a real simulation."""

import asyncio
import base64
import json
import sys
from types import SimpleNamespace
from unittest.mock import AsyncMock

import httpx
import pytest
from fastapi import Request
from fastapi.testclient import TestClient

from app import checkpoint, hosted, public_view
from app.hosted_worker import HostedWorker

WORLD_ID = "a" * 32
REVISION = "aabbccdd-5"
PATH = f"/public/worlds/{WORLD_ID}"
PNG = b"\x89PNG\r\n\x1a\npreview"


def _checkpoint(tmp_path, monkeypatch):
    monkeypatch.setenv(checkpoint.CHECKPOINT_DIR_ENV, str(tmp_path))
    root = checkpoint.CheckpointStore(tmp_path)
    record = {
        "world_id": WORLD_ID, "name": "Sleep test", "is_public": True, "is_default": True,
        "created_at": 1.0, "updated_at": 1.0,
    }
    root.save_index([record])
    store = root.for_world(WORLD_ID)
    store.directory.mkdir(parents=True)
    store.world_path.write_bytes(b"checkpoint")
    world = {"revision": REVISION, "seed": 7, "elapsed_years": 100_000.0}
    store.save_hosted_summary(world)
    store.save_viewer_image({
        "version": checkpoint.VIEWER_IMAGE_VERSION,
        "generation": store.world_generation(),
        "key": {
            "seed": 7, "elapsed_years": 100_000.0, "view": public_view.GUEST_DEFAULT_VIEW,
            "width": public_view.GUEST_DEFAULT_WIDTH, "height": public_view.GUEST_DEFAULT_HEIGHT,
            "rotation": list(public_view.GUEST_IDENTITY_ROTATION),
        },
        "image_base64": base64.b64encode(PNG).decode(),
    })
    return root, record, store


def _preview(client, revision=REVISION):
    return client.get(f"{PATH}/preview.png", params={
        "revision": revision, "view": public_view.GUEST_DEFAULT_VIEW,
        "rotation": ",".join(str(value) for value in public_view.GUEST_IDENTITY_ROTATION),
    })


def test_saved_map_and_state_serve_before_waking(tmp_path, monkeypatch):
    _checkpoint(tmp_path, monkeypatch)
    warm = AsyncMock()
    monkeypatch.setattr(HostedWorker, "warm", warm)
    app = hosted.create_app()
    with TestClient(app) as client:
        assert client.get("/healthz").status_code == 200
        assert client.get("/public/worlds").json()["worlds"][0]["revision"] == REVISION
        state = client.get(f"{PATH}/world").json()
        assert state["snapshot_first"] and state["world"]["revision"] == REVISION
        assert _preview(client).content == PNG
        assert _preview(client, "old").status_code == 404
        warm.assert_not_called()
        assert client.post(f"{PATH}/warm").status_code == 202
        assert client.post(f"{PATH}/warm").status_code == 202
    warm.assert_awaited_once_with(WORLD_ID)


def test_snapshot_refuses_private_and_unmatched_checkpoint(tmp_path, monkeypatch):
    root, record, store = _checkpoint(tmp_path, monkeypatch)
    with TestClient(hosted.create_app()) as client:
        replacement = store.world_path.with_suffix(".tmp")
        replacement.write_bytes(b"replacement")
        replacement.replace(store.world_path)
        assert _preview(client).status_code == 404
        record["is_public"] = False
        root.save_index([record])
        assert client.get("/public/worlds").json() == {"worlds": []}
        assert client.get(f"{PATH}/world").status_code == 404
        assert client.post(f"{PATH}/warm").status_code == 404
        assert client.post(f"{PATH}/step", json={"revision": REVISION}).status_code == 404


def test_stale_step_does_not_start_worker(tmp_path, monkeypatch):
    _checkpoint(tmp_path, monkeypatch)
    with TestClient(hosted.create_app()) as client:
        response = client.post(f"{PATH}/step", json={"revision": "old"})
        assert response.json()["computed"] is False
        assert response.json()["world"]["revision"] == REVISION
        assert not client.app.state.worker.running


def test_saved_private_run_remains_visible_during_worker_recovery(tmp_path, monkeypatch):
    root, record, store = _checkpoint(tmp_path, monkeypatch)
    record["is_public"] = False
    root.save_index([record])
    store.save_status(checkpoint.RunStatus(running=True, step_years=100_000))
    monkeypatch.setattr(HostedWorker, "prepare", AsyncMock())
    with TestClient(hosted.create_app()) as client:
        assert client.get("/public/worlds").json()["worlds"][0]["running"]
        assert client.get(f"{PATH}/world").json()["running"]
        assert client.post(f"{PATH}/step", json={"revision": "old"}).json()["world"]["revision"] == REVISION
        assert _preview(client).status_code == 404


def test_stopped_checkpoint_restores_on_demand_with_the_same_revision(tmp_path, monkeypatch):
    from app import main
    from app.world import World
    from app.world_registry import WorldSession

    _, _, store = _checkpoint(tmp_path, monkeypatch)
    store.save_world(World(seed=7, elapsed_years=100_000))
    store.save_hosted_summary({"revision": REVISION, "seed": 7, "elapsed_years": 100_000})
    monkeypatch.setattr(main, "_hosted_worker", True)
    session = main._wire_session(WorldSession(WORLD_ID, "Sleep test", is_public=True))
    session.checkpoint_store = store
    main._restore_checkpoint_session(session)
    assert session.world is None
    assert session.feed.revision() == REVISION
    restore = session.runner.set_status

    def forbid_status_reload(status, persist=True):
        assert persist, "lazy hydration overwrote the live runner status"
        return restore(status, persist=persist)

    monkeypatch.setattr(session.runner, "set_status", forbid_status_reload)
    main._hydrate(session)
    assert session.world.elapsed_years == 100_000
    assert session.feed.revision() == REVISION


def test_failed_checkpoint_keeps_the_worker_dirty(tmp_path, monkeypatch):
    from app import main
    from app.world import World
    from app.world_registry import WorldSession

    _, _, store = _checkpoint(tmp_path, monkeypatch)
    monkeypatch.setattr(main, "_hosted_worker", True)
    session = main._wire_session(WorldSession(WORLD_ID, "Sleep test", is_public=True))
    session.world = World(seed=7)
    session.checkpoint_store = store

    def fail(_world):
        raise OSError("disk full")

    with monkeypatch.context() as patch:
        patch.setattr(store, "save_world", fail)
        assert not main._checkpoint_world(session)
        assert session.checkpoint_dirty
    session.feed.publish(session.world)
    assert main._checkpoint_world(session)
    assert not session.checkpoint_dirty
    assert store.load_hosted_summary()["revision"] == session.feed.revision()


def test_auth_guard_runs_before_proxy_and_worker_routes_are_private(tmp_path, monkeypatch):
    _checkpoint(tmp_path, monkeypatch)
    with TestClient(hosted.create_app(), follow_redirects=False) as client:
        assert client.get("/").headers["location"] == "/view"
        assert client.get("/worlds").status_code == 401
        assert client.get(f"/worlds/{WORLD_ID}/controls").status_code == 401
        assert client.get("/_worker/idle").status_code == 404
        assert client.post(f"{PATH}/warm", headers={"Origin": "https://other.example"}).status_code == 403
        assert not client.app.state.worker.running
    with TestClient(hosted.create_app(), base_url="http://localhost", client=("127.0.0.1", 50000)) as client:
        assert client.get("/worlds").status_code == 401
        assert not client.get("/auth/session").json()["admin"]
        assert client.get("/login", follow_redirects=False).status_code == 200


def test_idle_retirement_respects_active_requests_and_worker_veto(tmp_path, monkeypatch):
    async def run():
        worker = HostedWorker(tmp_path / "worker.sock", 0.01)
        worker._ready = True
        monkeypatch.setattr(HostedWorker, "running", property(lambda _worker: True))
        stop = AsyncMock()
        monkeypatch.setattr(worker, "_stop", stop)
        results = [
            httpx.Response(200, json={"idle": False}),
            httpx.Response(200, json={"idle": True}),
        ]
        probe = AsyncMock(side_effect=results)
        monkeypatch.setattr(worker._client, "get", probe)
        try:
            await worker.retire_if_idle()
            stop.assert_not_called()
            async with worker.lease(wake=False, activity=False):
                await worker.retire_if_idle()
                stop.assert_not_called()
            await worker.retire_if_idle()
            stop.assert_awaited_once()
        finally:
            await worker._client.aclose()
    asyncio.run(run())


def test_hosted_entry_point_does_not_import_simulation():
    import subprocess

    subprocess.run(
        [sys.executable, "-c", "import sys, app.hosted; assert not {'app.main', 'numpy', 'scipy', 'numba'} & sys.modules.keys()"],
        check=True,
    )


def test_first_publish_after_wake_cannot_reuse_an_uncommitted_revision(monkeypatch):
    tokens = iter(["ab" * 16, "cd" * 16])
    monkeypatch.setattr(public_view.secrets, "token_hex", lambda _size: next(tokens))
    world = SimpleNamespace(seed=7, elapsed_years=0)
    original = public_view.WorldFeed()
    original.publish(world)
    saved = original.summary()
    world.elapsed_years = 100
    original.publish(world)
    already_seen = original.revision()

    restarted = public_view.WorldFeed()
    restarted.restore(saved, "Sleep test")
    assert restarted.revision() == saved["revision"]
    world.elapsed_years = 200
    restarted.publish(world)
    assert restarted.revision() not in (already_seen, saved["revision"])


def test_failed_operation_rolls_back_and_allows_retirement(tmp_path, monkeypatch):
    from app import main
    from app.world import World
    from app.world_registry import WorldRegistry

    root, _, store = _checkpoint(tmp_path, monkeypatch)
    registry = WorldRegistry(max_worlds=1)
    monkeypatch.setattr(main, "_registry", registry)
    monkeypatch.setattr(main, "_checkpoint_root", root)
    monkeypatch.setattr(main, "_hosted_worker", True)
    monkeypatch.setattr(main, "_index_dirty", False)
    session = main._wire_session(registry.create("Sleep test", is_public=True))
    session.checkpoint_store = store
    session.world = World(seed=7)
    session.feed.publish(session.world)
    assert main._checkpoint_world(session)
    committed = session.feed.revision()
    session.world.elapsed_years = 100
    with main._compute_gate.hold(), session.lock:
        main._world_changed(session, completed=False)
    assert session.world is None
    assert session.feed.revision() == committed
    assert main._worker_idle() == {"idle": True}
    main._hydrate(session)
    assert session.world.elapsed_years == 0


def test_idle_probe_retries_failed_index_write(tmp_path, monkeypatch):
    from app import main
    from app.world_registry import WorldRegistry

    root, _, _ = _checkpoint(tmp_path, monkeypatch)
    registry = WorldRegistry(max_worlds=1)
    monkeypatch.setattr(main, "_registry", registry)
    monkeypatch.setattr(main, "_checkpoint_root", root)
    monkeypatch.setattr(main, "_index_dirty", True)
    session = main._wire_session(registry.create("Sleep test"))
    assert main._worker_idle() == {"idle": True}
    assert root.load_index() == [session.record()]


def test_cold_public_requests_do_not_wait_for_unbounded_hydration(tmp_path, monkeypatch):
    from app import main
    from app.world_registry import WorldRegistry, WorldSession

    _, _, store = _checkpoint(tmp_path, monkeypatch)
    registry = WorldRegistry(max_worlds=1)
    monkeypatch.setattr(main, "_registry", registry)
    monkeypatch.setattr(main, "_hosted_worker", True)
    monkeypatch.setattr(public_view, "GUEST_LOCK_WAIT_S", 0)
    monkeypatch.setattr(public_view, "GUEST_OWNER_STEP_WAIT_S", 0)
    session = main._wire_session(registry.adopt(WorldSession(WORLD_ID, is_public=True)))
    session.checkpoint_store = store
    main._restore_checkpoint_session(session)
    client = TestClient(main.app)
    with main._compute_gate.hold():
        assert client.get(f"{PATH}/render").status_code == 503
        assert client.post(f"{PATH}/step", json={"revision": REVISION}).status_code == 503
    assert session.world is None


def test_done_line_commits_partial_animation_batch_before_disconnect(tmp_path, monkeypatch):
    from app import main
    from app.world import World
    from app.world_registry import WorldRegistry

    root = checkpoint.CheckpointStore(tmp_path)
    registry = WorldRegistry(max_worlds=1)
    monkeypatch.setattr(main, "_hosted_worker", True)
    monkeypatch.setattr(main, "_checkpoint_root", root)
    monkeypatch.setattr(main, "_index_dirty", False)
    monkeypatch.setattr(main, "_registry", registry)
    session = main._wire_session(registry.create("Sleep test"))
    session.world = World(seed=7)

    def encode(world, *_args, **_kwargs):
        yield ("progress", 1, 2, b"png", 0.0)
        world.elapsed_years = 100
        yield ("done", b"video", True)

    monkeypatch.setattr(main.render_image, "stream_animation_mp4", encode)
    response = main.animate(session.world_id, main.AnimateRequest(
        projection="eckert4", view="elevation", width=64, height=64, num_frames=2, step_years=100, steps_per_frame=2,
    ))

    async def disconnect():
        await anext(response.body_iterator)
        line = json.loads(await anext(response.body_iterator))
        assert line["elapsed_years"] == 100
        await response.body_iterator.aclose()

    asyncio.run(disconnect())
    store = session.checkpoint_store
    assert store.load_hosted_summary()["elapsed_years"] == 100
    assert store.load()[0].elapsed_years == 100


def test_proxy_does_not_keep_cookies_or_hide_a_snapshot_when_the_worker_exits(tmp_path, monkeypatch):
    async def run():
        worker = HostedWorker(tmp_path / "worker.sock", 90)
        headers = []

        def handle(request):
            headers.append(request.headers.get("cookie"))
            if len(headers) == 3:
                worker.process = None
                raise httpx.ConnectError("worker exited", request=request)
            return httpx.Response(200, headers={"Set-Cookie": "mb_guest=someone-elses-allowance; Path=/"})

        cookies = worker._client.cookies.jar
        await worker._client.aclose()
        async with httpx.AsyncClient(
            base_url="http://simulation", transport=httpx.MockTransport(handle), cookies=cookies,
        ) as client:
            with monkeypatch.context() as patch:
                patch.setattr(worker, "_client", client)
                worker._ready = True
                worker.process = SimpleNamespace(returncode=None)
                first = Request({"type": "http", "headers": [(b"cookie", b"mb_guest=first")]})
                second = Request({"type": "http", "headers": []})
                assert (await worker.read("/first", first)).status_code == 200
                assert (await worker.read("/second", second)).status_code == 200
                assert headers == ["mb_guest=first", None]
                assert await worker.read("/second", second) is None
    asyncio.run(run())


def test_supervisor_retries_a_saved_run_and_recovers_after_a_child_exit(tmp_path, monkeypatch):
    root, _, store = _checkpoint(tmp_path, monkeypatch)
    store.save_status(checkpoint.RunStatus(running=True))
    monkeypatch.setattr(hosted, "_RECOVERY_RETRY_SECONDS", 0)
    calls = []
    worker = SimpleNamespace(ready=False, idle_seconds=1)

    async def prepare():
        calls.append("start")
        if len(calls) == 1:
            raise OSError("startup transient")
        worker.ready = True

    async def retire():
        calls.append("exit")
        worker.ready = False

    worker.prepare = prepare
    worker.retire_if_idle = retire

    async def run():
        task = asyncio.create_task(hosted._supervise_worker(root, worker))
        for _ in range(5):
            await asyncio.sleep(0)
        task.cancel()
        with pytest.raises(asyncio.CancelledError):
            await task

    asyncio.run(run())
    assert calls[:4] == ["start", "start", "exit", "start"]


def test_openapi_includes_worker_admin_routes_and_hosted_routes(tmp_path, monkeypatch):
    _checkpoint(tmp_path, monkeypatch)
    schema = AsyncMock(return_value={"paths": {"/worlds/{world_id}/generate": {"post": {}}}})
    monkeypatch.setattr(HostedWorker, "openapi", schema)
    with TestClient(hosted.create_app()) as client:
        assert client.get("/docs").status_code == 200
        paths = client.get("/openapi.json").json()["paths"]
        assert "/worlds/{world_id}/generate" in paths
        assert f"{hosted._WORLD_PATH}/preview.png" in paths
        assert client.get("/openapi.json").status_code == 200
    schema.assert_awaited_once()


def test_invalid_idle_timeout_names_the_environment_variable(monkeypatch):
    monkeypatch.setenv("MANTLE_BLOOM_IDLE_SECONDS", "invalid")
    with pytest.raises(ValueError, match="MANTLE_BLOOM_IDLE_SECONDS"):
        hosted._idle_seconds()

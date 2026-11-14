"""Exercise the snapshot API and idle lifecycle without starting a real simulation."""

import asyncio
import base64
import sys
from unittest.mock import AsyncMock

import httpx
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
            await worker._retire_if_idle()
            stop.assert_not_called()
            async with worker.lease(wake=False, activity=False):
                await worker._retire_if_idle()
                stop.assert_not_called()
            await worker._retire_if_idle()
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

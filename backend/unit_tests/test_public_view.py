import threading
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

from app import main, public_view
from app.main import app


@pytest.fixture(autouse=True)
def _env(monkeypatch):
    # No OAuth app: remote visitors are plain guests, and the allowance cookie is signed with
    # the per-process fallback secret.
    for name in ("MANTLE_BLOOM_GITHUB_CLIENT_ID", "MANTLE_BLOOM_SESSION_SECRET", "MANTLE_BLOOM_LOCAL_NO_LOGIN"):
        monkeypatch.delenv(name, raising=False)


@pytest.fixture
def steps(monkeypatch):
    """A stand-in world and step: each step advances elapsed time, and the list counts how many
    really ran. The routes only need seed/elapsed_years, so no real world is generated."""
    world = SimpleNamespace(seed=7, elapsed_years=0.0)
    ran = []

    def fake_step(w, years):
        ran.append(years)
        w.elapsed_years += years
        yield 1.0

    monkeypatch.setattr(main, "step_world_progress", fake_step)
    main._state["world"] = world
    yield ran
    main._state["world"] = None
    assert not main._world_lock.locked()


@pytest.fixture
def guest():
    return TestClient(app, base_url="https://testserver")


def _local() -> TestClient:
    return TestClient(app, base_url="http://localhost", client=("127.0.0.1", 50000))


def _revision(client) -> str:
    return client.get("/public/world").json()["world"]["revision"]


def _step(client, revision=None):
    return client.post("/public/step", json={"revision": revision or _revision(client)})


def test_world_is_null_before_one_exists(guest):
    main._state["world"] = None
    body = guest.get("/public/world").json()
    assert body["world"] is None
    assert body["steps_remaining"] == public_view.GUEST_STEP_LIMIT


@pytest.mark.parametrize(
    "method,path",
    [
        ("post", "/world/generate"),
        ("get", "/world/save"),
        ("post", "/world/load"),
        ("post", "/world/step"),
        ("post", "/world/animate"),
        ("post", "/world/animate/stop"),
        ("post", "/world/controls"),
    ],
)
def test_guest_cannot_create_save_load_or_run_the_world(guest, steps, method, path):
    assert getattr(guest, method)(path).status_code == 401
    assert steps == []


def test_guest_step_computes_and_charges_one(guest, steps):
    before = _revision(guest)
    resp = _step(guest, before)
    assert resp.status_code == 200
    body = resp.json()
    assert body["computed"] is True
    assert body["steps_remaining"] == public_view.GUEST_STEP_LIMIT - 1
    assert body["world"]["revision"] != before
    assert steps == [public_view.GUEST_STEP_YEARS]
    assert guest.get("/public/world").json()["steps_remaining"] == public_view.GUEST_STEP_LIMIT - 1


def test_stale_revision_is_served_without_computing_or_charging(guest, steps):
    other = TestClient(app, base_url="https://testserver")
    stale = _revision(guest)
    _step(other, stale)
    resp = _step(guest, stale)
    assert resp.json()["computed"] is False
    assert resp.json()["world"]["revision"] == _revision(other)
    assert "set-cookie" not in resp.headers
    assert guest.get("/public/world").json()["steps_remaining"] == public_view.GUEST_STEP_LIMIT
    assert len(steps) == 1


def test_allowance_runs_out_after_the_limit(guest, steps):
    for _ in range(public_view.GUEST_STEP_LIMIT):
        assert _step(guest).json()["computed"] is True
    current = _revision(guest)
    assert _step(guest, current).status_code == 429
    assert len(steps) == public_view.GUEST_STEP_LIMIT
    assert guest.get("/public/world").json()["steps_remaining"] == 0
    # A world someone else advanced is still free to catch up to.
    _local().post("/public/step", json={"revision": current})
    assert _step(guest, current).json()["computed"] is False


def test_forged_allowance_cookie_is_ignored(guest, steps):
    guest.cookies.set(public_view.GUEST_COOKIE, "eyJ1c2VkIjotNX0.forged")
    assert guest.get("/public/world").json()["steps_remaining"] == public_view.GUEST_STEP_LIMIT


def test_owner_steps_from_the_viewer_without_an_allowance(steps):
    client = _local()
    resp = _step(client)
    assert resp.json()["computed"] is True
    assert resp.json()["steps_remaining"] is None
    assert "set-cookie" not in resp.headers


def test_step_while_owner_holds_the_world_is_busy(guest, steps):
    revision = _revision(guest)
    main._world_lock.acquire()
    try:
        assert _step(guest, revision).status_code == 503
    finally:
        main._world_lock.release()
    assert steps == []
    assert guest.get("/public/world").json()["steps_remaining"] == public_view.GUEST_STEP_LIMIT


def test_cross_site_step_is_refused(guest, steps):
    resp = guest.post("/public/step", json={"revision": _revision(guest)}, headers={"Origin": "https://evil.example"})
    assert resp.status_code == 403
    assert steps == []


def test_render_is_fixed_to_eckert4_and_the_viewer_views(guest, steps, monkeypatch):
    calls = []
    monkeypatch.setattr(main.render_image, "render_png_base64", lambda *args: calls.append(args) or "png")
    resp = guest.get("/public/render", params={"view": "biome", "projection": "behrmann", "rotation": "0,1,0,1,0,0,0,0,1"})
    assert resp.status_code == 200
    assert resp.json()["image_base64"] == "png"
    _, projection, view, *_rest, rotation, _m, _p = calls[0]
    assert (projection, view) == ("eckert4", "biome")
    assert rotation.tolist() == main._parse_view_rotation(None).tolist()
    assert guest.get("/public/render", params={"view": "temperature"}).status_code == 400
    assert guest.get("/public/render", params={"width": 4000}).status_code == 400


def test_render_gives_up_while_owner_holds_the_world(guest, steps, monkeypatch):
    monkeypatch.setattr(public_view, "GUEST_LOCK_WAIT_S", 0.01)
    main._world_lock.acquire()
    try:
        assert guest.get("/public/render").status_code == 503
    finally:
        main._world_lock.release()


def test_concurrent_requests_share_one_computation():
    world = SimpleNamespace(seed=1, elapsed_years=0.0)
    lock = threading.Lock()
    coordinator = public_view.GuestStepCoordinator()
    started, release = threading.Event(), threading.Event()
    computed = []

    def compute(w):
        computed.append(1)
        started.set()
        release.wait(5)
        w.elapsed_years += 1

    def request(results):
        results.append(coordinator.request_step(
            public_view.world_revision(world), True, lambda: world, lock, compute
        ).status)

    first, second = [], []
    t1 = threading.Thread(target=request, args=(first,))
    t1.start()
    assert started.wait(5)
    t2 = threading.Thread(target=request, args=(second,))
    t2.start()
    release.set()
    t1.join(5)
    t2.join(5)
    assert (first, second, computed) == (["computed"], ["shared"], [1])
    assert not lock.locked()


def test_waiting_on_a_stuck_step_gives_up_as_busy():
    world = SimpleNamespace(seed=1, elapsed_years=0.0)
    coordinator = public_view.GuestStepCoordinator()
    coordinator._inflight = threading.Event()
    outcome = coordinator.request_step(
        public_view.world_revision(world), True, lambda: world, threading.Lock(), lambda w: None, inflight_wait_s=0.01
    )
    assert outcome.status == "busy"


def test_failed_step_releases_the_world_and_charges_nothing():
    world = SimpleNamespace(seed=1, elapsed_years=0.0)
    lock = threading.Lock()
    coordinator = public_view.GuestStepCoordinator()

    def boom(w):
        raise RuntimeError("degenerate geometry")

    with pytest.raises(RuntimeError):
        coordinator.request_step(public_view.world_revision(world), True, lambda: world, lock, boom)
    assert not lock.locked()
    assert coordinator._inflight is None

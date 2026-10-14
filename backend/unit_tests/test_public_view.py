import asyncio
import threading
import time
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
    main._set_world(world)
    yield ran
    main._set_world(None)
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
    main._set_world(None)
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
        ("get", "/world/run"),
        ("post", "/world/run"),
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


def _bumping(world, feed, extra=None):
    """A guest compute that advances `world` and publishes it, as main._run_guest_step does."""
    def compute():
        if extra is not None:
            extra()
        world.elapsed_years += 1
        feed.publish(world)
    return compute


def test_concurrent_requests_share_one_computation():
    world = SimpleNamespace(seed=1, elapsed_years=0.0)
    feed = public_view.WorldFeed()
    feed.publish(world)
    lock = threading.Lock()
    coordinator = public_view.GuestStepCoordinator()
    started, release = threading.Event(), threading.Event()
    computed = []

    def hold():
        computed.append(1)
        started.set()
        release.wait(5)

    seen = feed.revision()

    def request(results):
        results.append(coordinator.request_step(seen, True, feed.revision, lock, _bumping(world, feed, hold)).status)

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
    feed = public_view.WorldFeed()
    coordinator = public_view.GuestStepCoordinator()
    coordinator._inflight = threading.Event()
    outcome = coordinator.request_step(
        feed.revision(), True, feed.revision, threading.Lock(), lambda: None, inflight_wait_s=0.01
    )
    assert outcome.status == "busy"


def test_failed_step_releases_the_world_and_charges_nothing():
    feed = public_view.WorldFeed()
    lock = threading.Lock()
    coordinator = public_view.GuestStepCoordinator()

    def boom():
        raise RuntimeError("degenerate geometry")

    with pytest.raises(RuntimeError):
        coordinator.request_step(feed.revision(), True, feed.revision, lock, boom)
    assert not lock.locked()
    assert coordinator._inflight is None


def test_world_changed_after_the_revision_check_is_not_charged():
    # The owner changes the world between the coordinator's revision check and its taking the
    # world lock: the guest asked to step a world they haven't seen, so nothing is computed.
    world = SimpleNamespace(seed=1, elapsed_years=0.0)
    feed = public_view.WorldFeed()
    feed.publish(world)
    seen = feed.revision()

    class OwnerSlipsIn:
        def acquire(self, *_args, **_kwargs):
            feed.publish(world)
            return True

        def release(self):
            self.released = True

    lock = OwnerSlipsIn()
    outcome = public_view.GuestStepCoordinator().request_step(seen, True, feed.revision, lock, pytest.fail)
    assert outcome.status == "shared"
    assert lock.released


# --- Revisions, viewer updates and the render cache (#374) -----------------------------------


def test_revision_moves_on_every_publish_and_differs_across_processes():
    feed, other = public_view.WorldFeed(), public_view.WorldFeed()
    world = SimpleNamespace(seed=3, elapsed_years=5.0)
    feed.publish(world)
    other.publish(world)
    first = feed.revision()
    assert first != other.revision()
    assert feed.summary() == {"revision": first, "seed": 3, "elapsed_years": 5.0}
    # Same seed and age (e.g. reloading the same save) is still a new revision.
    feed.publish(world)
    assert feed.revision() != first
    feed.publish(None)
    assert feed.summary() is None


def test_wait_past_returns_when_another_thread_publishes():
    feed = public_view.WorldFeed()
    seen = feed.revision()

    async def wait():
        loop = asyncio.get_running_loop()
        loop.call_later(0.05, lambda: threading.Thread(target=feed.publish, args=(None,)).start())
        start = time.monotonic()
        await feed.wait_past(seen, 5.0)
        return time.monotonic() - start

    assert asyncio.run(wait()) < 1.0
    assert feed.revision() != seen
    assert feed._waiters == []


def test_wait_past_times_out_and_returns_at_once_for_a_stale_revision():
    feed = public_view.WorldFeed()
    asyncio.run(feed.wait_past(feed.revision(), 0.01))
    assert feed._waiters == []
    stale = feed.revision()
    feed.publish(None)
    asyncio.run(asyncio.wait_for(feed.wait_past(stale, 60.0), 1.0))


def test_wake_ends_a_wait_whose_condition_changed():
    feed = public_view.WorldFeed()
    seen = feed.revision()
    flag = []

    async def scenario():
        waiter = asyncio.ensure_future(feed.wait_past(seen, 5.0, changed=lambda: bool(flag)))
        await asyncio.sleep(0.01)
        assert not waiter.done()
        flag.append(True)
        threading.Thread(target=feed.wake).start()
        await asyncio.wait_for(waiter, 2.0)
        # Already changed: returns without waiting.
        await asyncio.wait_for(feed.wait_past(seen, 5.0, changed=lambda: bool(flag)), 0.5)

    asyncio.run(scenario())
    assert feed.revision() == seen


def test_render_cache_serves_only_the_current_revision():
    cache = public_view.RenderCache(max_entries=2)
    cache.put("r1", ("biome",), {"image_base64": "a"})
    assert cache.get("r1", ("biome",)) == {"image_base64": "a"}
    assert cache.get("r2", ("biome",)) is None
    cache.put("r2", ("elevation",), {"image_base64": "b"})
    assert cache.get("r1", ("biome",)) is None
    cache.put("r2", ("biome",), {"image_base64": "c"})
    cache.put("r2", ("biome", 2), {"image_base64": "d"})
    assert cache.get("r2", ("elevation",)) is None  # oldest evicted past the cap
    cache.clear()
    assert cache.get("r2", ("biome",)) is None


def test_viewers_share_one_render_per_revision(guest, steps, monkeypatch):
    calls = []
    monkeypatch.setattr(main.render_image, "render_png_base64", lambda *args: calls.append(args) or f"png{len(calls)}")
    first = guest.get("/public/render", params={"view": "biome"}).json()
    other = TestClient(app, base_url="https://testserver").get("/public/render", params={"view": "biome"}).json()
    assert other == first
    assert len(calls) == 1
    assert guest.get("/public/render", params={"view": "elevation"}).json()["image_base64"] == "png2"
    _step(guest)
    after = guest.get("/public/render", params={"view": "biome"}).json()
    assert after["image_base64"] == "png3"
    assert after["revision"] == _revision(guest) != first["revision"]
    assert after["elapsed_years"] == public_view.GUEST_STEP_YEARS


def test_cached_render_is_served_while_the_owner_holds_the_world(guest, steps, monkeypatch):
    monkeypatch.setattr(main.render_image, "render_png_base64", lambda *args: "png")
    monkeypatch.setattr(public_view, "GUEST_LOCK_WAIT_S", 0.01)
    guest.get("/public/render")
    main._world_lock.acquire()
    try:
        assert guest.get("/public/render").json()["image_base64"] == "png"
    finally:
        main._world_lock.release()


class _AnyWorld(SimpleNamespace):
    """Answers 0 for any attribute the owner's routes read that a test didn't set."""

    def __getattr__(self, name):
        return 0


@pytest.mark.parametrize("change", ["controls", "load"])
def test_owner_changes_invalidate_viewer_renders(steps, monkeypatch, change):
    monkeypatch.setattr(main.render_image, "render_png_base64", lambda *args: "png")
    monkeypatch.setattr(main, "_summary", lambda world: {})
    main._set_world(_AnyWorld(seed=7, elapsed_years=0.0))
    client = _local()
    before = client.get("/public/render").json()["revision"]
    if change == "controls":
        monkeypatch.setattr(main.climate, "compute_climate", lambda *args: None)
        assert client.post("/world/controls", json={"solar_multiplier": 1.1}).status_code == 200
    else:
        monkeypatch.setattr(main.persistence, "load_world_bytes", lambda body: _AnyWorld(
            seed=7, elapsed_years=0.0, hydrology_cache=None))
        assert client.post("/world/load", content=b"save").status_code == 200
    assert client.get("/public/render").json()["revision"] != before


def test_load_while_the_world_is_busy_is_refused(steps, monkeypatch):
    monkeypatch.setattr(main.persistence, "load_world_bytes", lambda body: _AnyWorld(
        seed=8, elapsed_years=0.0, hydrology_cache=None))
    before = _revision(_local())
    main._world_lock.acquire()
    try:
        assert _local().post("/world/load", content=b"save").status_code == 503
    finally:
        main._world_lock.release()
    assert _revision(_local()) == before


def test_viewer_long_poll_wakes_on_another_viewers_step(guest, steps):
    revision = _revision(guest)
    result = {}
    poller = threading.Thread(target=lambda: result.update(guest.get("/public/world", params={"after": revision}).json()))
    poller.start()
    time.sleep(0.05)
    assert poller.is_alive()
    stepped = _step(TestClient(app, base_url="https://testserver"), revision).json()
    poller.join(5)
    assert result["world"] == stepped["world"]
    assert result["steps_remaining"] == public_view.GUEST_STEP_LIMIT
    assert len(steps) == 1


def test_viewer_long_poll_answers_the_current_world_on_timeout(guest, steps, monkeypatch):
    monkeypatch.setattr(public_view, "GUEST_POLL_WAIT_S", 0.01)
    revision = _revision(guest)
    assert guest.get("/public/world", params={"after": revision}).json()["world"]["revision"] == revision


def test_duplicate_concurrent_step_requests_run_the_step_once(guest, steps, monkeypatch):
    started, release = threading.Event(), threading.Event()

    def slow_step(w, years):
        steps.append(years)
        started.set()
        release.wait(5)
        w.elapsed_years += years
        yield 1.0

    monkeypatch.setattr(main, "step_world_progress", slow_step)
    revision = _revision(guest)
    results = []

    def request():
        results.append(_step(TestClient(app, base_url="https://testserver"), revision).json())

    threads = [threading.Thread(target=request) for _ in range(3)]
    threads[0].start()
    assert started.wait(5)
    for t in threads[1:]:
        t.start()
    time.sleep(0.05)
    release.set()
    for t in threads:
        t.join(5)
    assert len(steps) == 1
    assert sorted(r["computed"] for r in results) == [False, False, True]
    assert len({r["world"]["revision"] for r in results}) == 1
    # A later request against the new revision is a new step.
    assert _step(guest, results[0]["world"]["revision"]).json()["computed"] is True
    assert len(steps) == 2

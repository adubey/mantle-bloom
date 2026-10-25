import asyncio
import base64
import gc
import io
import json
import math
import threading
import time
import av
import numpy as np
import pytest
from fastapi.testclient import TestClient
from PIL import Image
from app.lithosphere_plate import generate_plates
from app import main, persistence
from app.hydroclimate import hydrology
from app.main import app
from app.plates import MAX_AUTO_PLATES, MIN_AUTO_PLATES

from .river_fixtures import spill_hop_fields


_default_session = None


def _reset_app_state():
    """Drop every registered world and create one fresh empty default session
    (#563). Tests address it through the `_WorldClient` wrapper below, which
    maps the old unqualified `/world/*` paths onto its `/worlds/{id}/*` paths,
    so test bodies keep the URLs they were written against."""
    global _default_session
    for session in main._registry.list():
        if session.runner is not None:
            session.runner.close()
        main._registry.remove(session.world_id)
    _default_session = main._registry.create(name="test")
    main._wire_session(_default_session)
    _default_session.is_public = True
    _default_session.world = None
    _default_session.animation_stop_event.clear()


class _WorldClient:
    """Maps `/world/*` (and `/world/debug_scenarios`) onto the default test
    world's `/worlds/{id}/*` (and `/debug/scenarios`) paths."""

    def __init__(self, inner):
        self._inner = inner

    def _map(self, url: str) -> str:
        if url.startswith("/world/debug_scenarios"):
            return "/debug/scenarios" + url[len("/world/debug_scenarios"):]
        if url == "/world" or url.startswith("/world/") or url.startswith("/world?"):
            return f"/worlds/{_default_session.world_id}" + url[len("/world"):]
        return url

    def get(self, url, **kwargs):
        return self._inner.get(self._map(url), **kwargs)

    def post(self, url, **kwargs):
        return self._inner.post(self._map(url), **kwargs)


@pytest.fixture
def client():
    """A TestClient over a fresh app.main: one empty default world, no pending
    animation stop. Without this reset each test saw whatever world the
    previous one left behind. A test that needs a populated world generates
    its own. Teardown also fails a test that leaves the session lock held
    (releasing it first, so the leak doesn't cascade into every later test
    as a 503 or a hang)."""
    _reset_app_state()
    # A loopback client on a loopback host: the local no-login path, so these tests exercise
    # the routes rather than the auth guard (test_auth.py covers that, #372).
    yield _WorldClient(TestClient(app, base_url="http://localhost", client=("127.0.0.1", 50000)))
    _reset_app_state()
    if _default_session.lock.locked():
        _default_session.lock.release()
        pytest.fail("test left the default session lock held")


def _decode_image(body: dict) -> Image.Image:
    return Image.open(io.BytesIO(base64.b64decode(body["image_base64"])))


class _NdjsonResponse:
    """Wraps a TestClient response from /world/generate or /world/step -- both now stream
    NDJSON progress (see backend app/main.py, issue #195) instead of returning one JSON body
    -- so the many existing `resp.json()` assertions below still work unmodified: `.json()`
    returns the stream's final `{"type": "done", ...}` message's fields, same shape the old
    single-JSON-body response used to have. A non-200 response (400/404/503, raised before
    the stream ever starts) is still ordinary single-line JSON, so `.json()` falls back to
    the real response's own method for those. Every other attribute (`.status_code`, `.text`,
    `.headers`, ...) passes straight through to the wrapped response, same as the real
    `httpx.Response` would give -- see test_animate_* below for the lower-level
    `resp.text.splitlines()` parsing this wraps, used directly where a test needs the
    intermediate `"progress"` lines rather than just the final result."""

    def __init__(self, resp):
        self._resp = resp

    def __getattr__(self, name):
        return getattr(self._resp, name)

    def json(self):
        if self._resp.status_code != 200:
            return self._resp.json()
        messages = [json.loads(line) for line in self._resp.text.splitlines() if line.strip()]
        done = messages[-1]
        assert done["type"] == "done", f"expected a final done message, got {done}"
        return done


def _post_generate(client, **kwargs):
    return _NdjsonResponse(client.post("/world/generate", **kwargs))


def _post_step(client, **kwargs):
    return _NdjsonResponse(client.post("/world/step", **kwargs))


def _post_streaming_then_disconnect(path: str, body: dict, disconnect_after_chunks: int) -> int:
    """Drives app.main's ASGI app directly for one streaming POST, with a hand-rolled
    receive/send pair (the same technique Starlette's own test suite uses to simulate a
    disconnect) that reports the client as gone exactly once `send()` has recorded
    `disconnect_after_chunks` response body chunks -- TestClient's httpx-based streaming
    doesn't expose enough control over response timing to cut a stream short at an exact
    chunk deterministically, which is why this bypasses it. Used to prove /world/generate's
    and /world/step's own GeneratorExit handling (see their comments on it) really does
    finish the underlying world mutation when a client walks away mid-stream, rather than
    leaving `_state["world"]` half-built/half-stepped forever. Returns the number of body
    chunks actually delivered before the simulated disconnect (i.e. `disconnect_after_chunks`
    itself, confirming the scenario really did cut the stream short rather than racing it to
    a normal finish)."""

    async def _drive() -> int:
        body_bytes = json.dumps(body).encode()
        # Loopback client and host, like the `client` fixture, so the auth guard lets it through (#372).
        scope = {
            "type": "http",
            "asgi": {"version": "3.0", "spec_version": "2.3"},
            "http_version": "1.1",
            "method": "POST",
            "scheme": "http",
            "path": path,
            "raw_path": path.encode(),
            "query_string": b"",
            "root_path": "",
            "headers": [
                (b"host", b"localhost"),
                (b"content-type", b"application/json"),
                (b"content-length", str(len(body_bytes)).encode()),
            ],
            "client": ("127.0.0.1", 123),
            "server": ("localhost", 80),
        }

        request_body_sent = False
        disconnect_now = asyncio.Event()
        chunk_count = 0

        async def receive():
            nonlocal request_body_sent
            if not request_body_sent:
                request_body_sent = True
                return {"type": "http.request", "body": body_bytes, "more_body": False}
            # Starlette's listen_for_disconnect polls receive() in a loop once the response
            # starts streaming -- block until send() below has recorded the target chunk
            # count, then report the client as gone.
            await disconnect_now.wait()
            return {"type": "http.disconnect"}

        async def send(message):
            nonlocal chunk_count
            if message["type"] == "http.response.body":
                chunk_count += 1
                if chunk_count == disconnect_after_chunks:
                    disconnect_now.set()

        await app(scope, receive, send)

        # The GeneratorExit-triggered drain only runs once nothing references the endpoint's
        # closed-over `_stream()` generator any more -- CPython's refcounting GC calls
        # `.close()` on it, raising GeneratorExit at its suspended yield, as soon as the last
        # reference drops, which can be slightly after `await app(...)` itself returns.
        # Forcing collection here makes that prompt; the caller still polls for its own
        # expected end state afterward (see _wait_until below) rather than this function
        # guessing how long the drain's actual plate-geometry/world-generation work takes.
        gc.collect()
        return chunk_count

    return asyncio.run(_drive())


def _wait_until(predicate, timeout: float = 5.0, interval: float = 0.05) -> bool:
    """Polls `predicate()` until it returns truthy or `timeout` seconds pass -- used below to
    wait for _post_streaming_then_disconnect's GeneratorExit-triggered drain to actually
    finish applying its world mutation, without hardcoding how long that background work
    should take (a fixed sleep long enough for a 10-plate world locally could still be too
    short under CI load). Returns the last (falsy) result on timeout rather than raising, so
    the caller's own `assert` gives a clearer failure message than a bare timeout would."""
    deadline = time.monotonic() + timeout
    result = predicate()
    while not result and time.monotonic() < deadline:
        time.sleep(interval)
        result = predicate()
    return result


def test_client_fixture_discards_state_left_by_an_earlier_test(client):
    # Stand-ins for what a previous test could leave behind: its generated world and an
    # animate stop that nothing consumed. A fresh setup must clear both, or the 404 tests
    # below depend on their position in the file (#334).
    _default_session.world = object()
    _default_session.animation_stop_event.set()
    _reset_app_state()
    assert _default_session.world is None
    assert not _default_session.animation_stop_event.is_set()
    assert client.get("/world/stats").status_code == 404


def test_generate_forwards_cratons_first_choice(client, monkeypatch):
    captured = {}

    def fake_generate_world_progress(seed, **kwargs):
        captured.update(kwargs)
        yield ("world", object())

    monkeypatch.setattr(main, "generate_world_progress", fake_generate_world_progress)
    monkeypatch.setattr(main, "_set_world", lambda session, world: None)
    monkeypatch.setattr(main, "_summary", lambda session, world: {"seed": 7})

    response = _post_generate(client, json={"seed": 7, "cratons_first": True})

    assert response.json()["seed"] == 7
    assert captured["cratons_first"] is True


def test_river_summary_draws_a_spill_hop_along_neighbour_edges():
    # Issue #520: the River Inspector's segments follow sink 2's spill path over the dry basin,
    # not one straight segment from the sink to node 6 beyond the rim.
    fields = spill_hop_fields(lake_m=0.0)
    (river,) = hydrology.group_rivers(fields)

    summary = main._river_summary(fields, river, 0)

    step = float(np.linalg.norm(fields.points[1] - fields.points[0]))
    lengths = [np.linalg.norm(np.subtract(b, a)) for a, b in summary["segments"]]
    assert len(lengths) == 7
    assert max(lengths) == pytest.approx(step, rel=1e-3)


def test_render_before_generate_returns_404(client):
    resp = client.get("/world/render")
    assert resp.status_code == 404


def test_cross_section_route_renders_image_and_keeps_request_validation(client, monkeypatch):
    assert client.get("/world/cross-section").status_code == 404

    _default_session.world = object()
    monkeypatch.setattr(
        main.cross_section_view,
        "render_png_base64",
        lambda *args, **kwargs: ("cG5n", {"image_base64": "cG5n", "world_revision": "revision"}),
    )
    response = client.get("/world/cross-section")
    assert response.status_code == 200
    assert response.json()["image_base64"] == "cG5n"
    assert response.json()["world_revision"] == "revision"
    assert client.get("/world/cross-section", params={"position_index": 120}).status_code == 400
    assert client.get("/world/cross-section", params={"width": main.MAX_RENDER_DIMENSION_PX + 1}).status_code == 400


def test_step_before_generate_returns_404(client):
    resp = _post_step(client, json={"years": 1_000_000})
    assert resp.status_code == 404


def test_step_rejects_invalid_years(client):
    # #557: years must be positive, finite, and under the per-request ceiling;
    # each step costs seconds of CPU regardless of its size. Non-finite values
    # go as raw bodies: httpx's json= refuses to serialize them.
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    before = _default_session.world.elapsed_years
    for years in (0, -1_000_000, main.MAX_STEP_YEARS + 1):
        resp = _post_step(client, json={"years": years})
        assert resp.status_code == 400, years
    for raw in (b'{"years": Infinity}', b'{"years": NaN}'):
        resp = _post_step(
            client, content=raw, headers={"Content-Type": "application/json"}
        )
        assert resp.status_code == 400, raw
    # Rejected before any stepping: the world has not advanced.
    assert _default_session.world.elapsed_years == before


def test_step_accepts_ui_step_sizes_and_run_accepts_the_max(client):
    # #557 verification: ordinary UI step sizes still work, and exactly
    # MAX_STEP_YEARS is allowed (the cap is inclusive). The max itself is
    # exercised through /world/run, which validates without stepping.
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = _post_step(client, json={"years": 10_000_000})
    assert resp.status_code == 200
    assert _default_session.world.elapsed_years == 10_000_000.0
    resp = client.post(
        "/world/run", json={"running": False, "step_years": main.MAX_STEP_YEARS}
    )
    assert resp.status_code == 200


def test_animate_rejects_invalid_step_years(client):
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    for step_years in (0, -1_000_000, main.MAX_STEP_YEARS + 1):
        resp = client.post(
            "/world/animate", json={"step_years": step_years, "num_frames": 2}
        )
        assert resp.status_code == 400, step_years
    for raw in (b'{"step_years": Infinity, "num_frames": 2}', b'{"step_years": NaN, "num_frames": 2}'):
        resp = client.post(
            "/world/animate", content=raw, headers={"Content-Type": "application/json"}
        )
        assert resp.status_code == 400, raw


def test_animate_rejects_out_of_range_steps_per_frame(client):
    # #557: steps_per_frame is capped as well as floored, so a request's total
    # step count (num_frames x steps_per_frame) is bounded through both factors.
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = client.post(
        "/world/animate",
        json={"step_years": 1_000_000, "steps_per_frame": main.MAX_STEPS_PER_FRAME + 1, "num_frames": 2},
    )
    assert resp.status_code == 400


def test_run_rejects_step_years_over_the_max(client):
    # #557: /world/run gets the same MAX_STEP_YEARS cap as /world/step.
    assert (
        client.post(
            "/world/run", json={"running": False, "step_years": main.MAX_STEP_YEARS + 1}
        ).status_code
        == 400
    )


def test_plates_before_generate_returns_404(client):
    assert client.get("/world/plates").status_code == 404


def test_plate_at_before_generate_returns_404(client):
    assert client.get("/world/plate_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_rivers_before_generate_returns_404(client):
    assert client.get("/world/rivers").status_code == 404


def test_river_at_before_generate_returns_404(client):
    assert client.get("/world/river_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_sample_at_before_generate_returns_404(client):
    assert client.get("/world/sample_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_node_at_before_generate_returns_404(client):
    assert client.get("/world/node_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_elevation_point_at_before_generate_returns_404(client):
    assert client.get("/world/elevation_point_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_lakes_before_generate_returns_404(client):
    assert client.get("/world/lakes").status_code == 404


def test_faults_before_generate_returns_404(client):
    assert client.get("/world/faults").status_code == 404
    assert client.get("/world/earthquakes").status_code == 404


def test_fault_at_before_generate_returns_404(client):
    assert client.get("/world/fault_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_lake_at_before_generate_returns_404(client):
    assert client.get("/world/lake_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_stranded_basins_before_generate_returns_404(client):
    assert client.get("/world/stranded_basins").status_code == 404


def test_volcanoes_before_generate_returns_404(client):
    assert client.get("/world/volcanoes").status_code == 404


def test_save_before_generate_returns_404(client):
    assert client.get("/world/save").status_code == 404


def test_animate_before_generate_returns_404(client):
    resp = client.post("/world/animate", json={"step_years": 1_000_000, "num_frames": 2})
    assert resp.status_code == 404


def test_export_hexgrid_before_generate_returns_404(client):
    assert client.post("/world/export_hexgrid", json={}).status_code == 404


def test_stats_before_generate_returns_404(client):
    assert client.get("/world/stats").status_code == 404


def test_stats_history_before_generate_returns_404(client):
    assert client.get("/world/stats_history").status_code == 404


def test_step_disconnect_mid_stream_still_finishes_the_step(client):
    # Regression test for a code-review finding on issue #195's streaming /world/step:
    # step_world_progress mutates `world` in place across many per-plate yield points, so if
    # main.py's step() just let a client disconnect abandon the generator wherever it
    # happened to be suspended, the world would be left permanently half-stepped -- some
    # plates shifted but not deformed, elapsed_years never advanced, no topology/climate/
    # erosion/sea-level pass. step()'s own `except GeneratorExit` drains step_world_progress
    # to completion first, so the step still fully applies -- proven here by disconnecting
    # after only 3 of a 10-plate step's 22 total NDJSON lines (21 progress fractions + done).
    _post_generate(client, json={"seed": 1, "num_plates": 10})
    chunk_count = _post_streaming_then_disconnect(f"/worlds/{_default_session.world_id}/step", {"years": 1_000_000}, disconnect_after_chunks=3)
    assert chunk_count == 3  # confirms the scenario really did cut the stream short

    _wait_until(lambda: client.get("/world/summary").json()["elapsed_years"] == 1_000_000.0)
    assert client.get("/world/summary").json()["elapsed_years"] == 1_000_000.0


def test_generate_disconnect_mid_stream_still_finishes_generating(client):
    # Same regression test as test_step_disconnect_mid_stream_still_finishes_the_step above,
    # for /world/generate -- disconnecting after just the first of generate_world_progress's
    # three phase-boundary progress lines (plate/site generation) should still leave a
    # complete, usable world in `_state` once mantle-center fitting and the
    # finish_generation bootstrap finish server-side.
    chunk_count = _post_streaming_then_disconnect(
        f"/worlds/{_default_session.world_id}/generate", {"seed": 1, "num_plates": 10}, disconnect_after_chunks=1,
    )
    assert chunk_count == 1  # confirms the scenario really did cut the stream short

    _wait_until(lambda: client.get("/world/summary").status_code == 200)
    summary = client.get("/world/summary").json()
    assert summary["seed"] == 1
    assert summary["num_plates"] == 10


def test_overlapping_step_returns_503(client, monkeypatch):
    _post_generate(client, json={"seed": 1, "num_plates": 6})

    # Makes the first /world/step's own critical section deterministically overlap the
    # second's: the blocking replacement only runs (and lets app._default_session.lock be released)
    # once entered_step_world confirms the first request is already holding the lock, so the
    # second request's 503 isn't a race on request scheduling order.
    entered_step_world = threading.Event()
    release_step_world = threading.Event()

    def blocking_step_world_progress(world, years):
        entered_step_world.set()
        release_step_world.wait(timeout=5)
        yield 1.0

    import app.main as main_module

    monkeypatch.setattr(main_module, "step_world_progress", blocking_step_world_progress)

    results: list[int] = []
    t1 = threading.Thread(target=lambda: results.append(_post_step(client, json={"years": 1_000_000}).status_code))
    t1.start()
    assert entered_step_world.wait(timeout=5)

    t2_response = _post_step(client, json={"years": 1_000_000})
    release_step_world.set()
    t1.join()

    assert t2_response.status_code == 503
    assert results == [200]


def test_overlapping_generate_returns_503(client, monkeypatch):
    # Same coverage as test_overlapping_step_returns_503 above, for /world/generate's own
    # (now inline, since the shared _reject_if_busy helper this used to go through was
    # removed -- see issue #195) non-blocking `_world_lock.acquire` + 503 path.
    entered_generate = threading.Event()
    release_generate = threading.Event()

    from app.world import generate_world_progress as real_generate_world_progress

    def blocking_generate_world_progress(seed, **kwargs):
        entered_generate.set()
        release_generate.wait(timeout=5)
        yield from real_generate_world_progress(seed, **kwargs)

    import app.main as main_module

    monkeypatch.setattr(main_module, "generate_world_progress", blocking_generate_world_progress)

    results: list[int] = []
    t1 = threading.Thread(
        target=lambda: results.append(_post_generate(client, json={"seed": 1, "num_plates": 6}).status_code)
    )
    t1.start()
    assert entered_generate.wait(timeout=5)

    t2_response = _post_generate(client, json={"seed": 2, "num_plates": 6})
    release_generate.set()
    t1.join()

    assert t2_response.status_code == 503
    assert results == [200]


def test_render_waits_for_an_in_progress_step_instead_of_racing_it(client, monkeypatch):
    # Regression test: GET /world/render used to run with no lock at all, so it could read
    # plates.collect_all_points/collect_all_lake_depth/etc. (each independently walking
    # world.plates) while an in-progress step's own deform() was mid-mutation, desyncing
    # those separately-gathered arrays and crashing with an IndexError -- confirmed directly
    # by deliberately racing a step against repeated renders before this fix existed. Same
    # blocking_step_world technique as test_overlapping_step_returns_503 above, but this
    # time proving render *waits out* the lock (and then succeeds) rather than either racing
    # past it or being rejected with a 503 the way the write endpoints are -- a render is a
    # read, called far more often, so it should wait a moment, not fail.
    _post_generate(client, json={"seed": 1, "num_plates": 6})

    entered_step_world = threading.Event()
    release_step_world = threading.Event()

    def blocking_step_world_progress(world, years):
        entered_step_world.set()
        release_step_world.wait(timeout=5)
        yield 1.0

    import app.main as main_module

    monkeypatch.setattr(main_module, "step_world_progress", blocking_step_world_progress)

    step_results: list[int] = []
    t1 = threading.Thread(target=lambda: step_results.append(_post_step(client, json={"years": 1_000_000}).status_code))
    t1.start()
    assert entered_step_world.wait(timeout=5)

    render_results: list[int] = []
    t2 = threading.Thread(target=lambda: render_results.append(client.get("/world/render").status_code))
    t2.start()
    # The render thread must actually be blocked on _world_lock right now, not racing through
    # -- confirmed by it still not having finished shortly after starting, while the step
    # still holds the lock.
    t2.join(timeout=0.5)
    assert t2.is_alive()
    assert render_results == []

    release_step_world.set()
    t1.join(timeout=5)
    t2.join(timeout=5)

    assert step_results == [200]
    assert render_results == [200]


def test_step_compute_routes_do_not_race_an_in_progress_step(client, monkeypatch):
    # Regression test: /world/animate, /world/controls, /world/stats and /world/generate all
    # ran with no lock, so each could launch a Numba parallel=True kernel (fluid_dynamics.py
    # etc.) concurrently with an in-flight step's own kernel -- the default `workqueue`
    # threading layer is not safe against that and crashes the interpreter. They now share
    # `_world_lock` with step/render: the long writes (animate) 503 like an overlapping step,
    # the quick/read holders (controls, stats, generate) block and wait it out.
    _post_generate(client, json={"seed": 1, "num_plates": 6})

    entered_step_world = threading.Event()
    release_step_world = threading.Event()

    def blocking_step_world_progress(world, years):
        entered_step_world.set()
        release_step_world.wait(timeout=5)
        yield 1.0

    import app.main as main_module

    monkeypatch.setattr(main_module, "step_world_progress", blocking_step_world_progress)

    step_results: list[int] = []
    t1 = threading.Thread(target=lambda: step_results.append(_post_step(client, json={"years": 1_000_000}).status_code))
    t1.start()
    assert entered_step_world.wait(timeout=5)

    # animate is a long-running write -- rejected outright while the step holds the lock.
    animate_resp = client.post("/world/animate", json={"step_years": 1_000_000, "num_frames": 2})
    assert animate_resp.status_code == 503

    # stats is a read -- it must block, not race straight through into a climate recompute.
    stats_results: list[int] = []
    t2 = threading.Thread(target=lambda: stats_results.append(client.get("/world/stats").status_code))
    t2.start()
    t2.join(timeout=0.5)
    assert t2.is_alive()
    assert stats_results == []

    release_step_world.set()
    t1.join(timeout=5)
    t2.join(timeout=5)

    assert step_results == [200]
    assert stats_results == [200]


def test_generate_returns_summary(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6})
    assert resp.status_code == 200
    body = resp.json()
    assert body["num_plates"] == 6
    assert body["elapsed_years"] == 0.0
    assert body["seed"] == 1
    assert body["surface"] == "quad"


def test_generate_defaults_to_quad_surfaces(client):
    from app.sparse_quad_patch import PlateWithSparseQuadPatch

    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "node_density": 0.5})

    assert resp.status_code == 200
    assert all(isinstance(plate, PlateWithSparseQuadPatch) for plate in _default_session.world.plates)


def test_generate_returns_a_generation_event(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "continental_fraction": 0.5})
    body = resp.json()
    assert len(body["events"]) == 1
    assert "6 plates" in body["events"][0]["message"]
    assert "3 continental" in body["events"][0]["message"]
    assert body["events"][0]["elapsed_years"] == 0.0


def test_debug_scenarios_lists_every_scenario_with_a_label(client):
    from app import debug_worlds

    resp = client.get("/world/debug_scenarios")
    assert resp.status_code == 200
    scenarios = resp.json()["scenarios"]
    assert {s["name"] for s in scenarios} == set(debug_worlds.DEBUG_SCENARIOS)
    assert all(isinstance(s["label"], str) and s["label"] for s in scenarios)


def test_generate_debug_returns_summary_with_diagnostics_on(client):

    resp = client.post("/world/generate_debug", json={"scenario": "two_plate_divergent", "seed": 3})
    assert resp.status_code == 200
    body = resp.json()
    assert body["num_plates"] == 2
    assert body["elapsed_years"] == 0.0
    assert body["seed"] == 3
    world = _default_session.world
    assert world.debug_diagnostics is True
    assert len(world.pinned_omegas) == 2


def test_generate_debug_unknown_scenario_returns_400(client):
    resp = client.post("/world/generate_debug", json={"scenario": "not_a_real_scenario"})
    assert resp.status_code == 400


def test_generate_debug_replaces_previous_world(client):
    _post_generate(client, json={"seed": 1, "num_plates": 6})
    resp = client.post("/world/generate_debug", json={"scenario": "four_plate_grid"})
    assert resp.status_code == 200
    assert resp.json()["num_plates"] == 4


def test_generate_with_continental_fraction_gives_exact_count(client):
    # crust_type per plate isn't part of the render response (see render_image.py) -- the
    # generation event log is the documented way to confirm the exact count (also covered by
    # test_generate_returns_a_generation_event above; this locks in a second seed/count pair).
    resp = _post_generate(client, json={"seed": 1, "num_plates": 10, "continental_fraction": 0.4})
    assert "4 continental" in resp.json()["events"][0]["message"]


def test_generate_with_land_fraction(client):
    resp = _post_generate(
        client,
        json={"seed": 1, "num_plates": 10, "continental_fraction": 0.7, "land_fraction": 0.29},
    )
    assert resp.status_code == 200


def test_generate_with_node_density(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "node_density": 4.0})
    assert resp.status_code == 200


def test_generate_with_the_coarsest_node_density(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "node_density": 0.5})
    assert resp.status_code == 200


def test_generate_with_unknown_node_density_returns_400(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "node_density": 2.5})
    assert resp.status_code == 400


def test_generate_with_the_max_allowed_voronoi_points(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "voronoi_points": 2000})
    assert resp.status_code == 200


def test_generate_with_voronoi_points_over_the_max_returns_400(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "voronoi_points": 2001})
    assert resp.status_code == 400


def test_generate_with_zero_voronoi_points_returns_400(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "voronoi_points": 0})
    assert resp.status_code == 400


def test_generate_with_climate_density(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "climate_density": 2.0})
    assert resp.status_code == 200


def test_generate_with_the_finest_climate_density(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "climate_density": 4.0})
    assert resp.status_code == 200


def test_generate_with_unknown_climate_density_returns_400(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "climate_density": 3.0})
    assert resp.status_code == 400


def test_render_defaults_to_elevation_view_at_1100x611(client):
    _post_generate(client, json={"seed": 3, "num_plates": 6})
    resp = client.get("/world/render")
    assert resp.status_code == 200
    assert _decode_image(resp.json()).size == (1100, 611)


def test_render_different_views_produce_different_images(client):
    """Not a pixel-exact check (that would just re-derive render_image.py's own math) --
    just confirms the `view` param actually changes what's drawn, at the HTTP layer."""
    _post_generate(client, json={"seed": 3, "num_plates": 8, "continental_fraction": 0.5})
    bodies = {
        view: client.get("/world/render", params={"view": view, "width": 300, "height": 200}).json()["image_base64"]
        for view in ("elevation", "plates", "platesDetail")
    }
    assert len(set(bodies.values())) == 3


def test_render_different_resolutions_produce_different_size_images(client):
    _post_generate(client, json={"seed": 3, "num_plates": 6})
    small = client.get("/world/render", params={"width": 200, "height": 100}).json()
    large = client.get("/world/render", params={"width": 800, "height": 400}).json()
    assert _decode_image(small).size == (200, 100)
    assert _decode_image(large).size == (800, 400)


def test_render_unknown_projection_returns_400(client):
    _post_generate(client, json={"seed": 4, "num_plates": 6})
    resp = client.get("/world/render", params={"projection": "mercator"})
    assert resp.status_code == 400


def test_render_unknown_view_returns_400(client):
    _post_generate(client, json={"seed": 4, "num_plates": 6})
    resp = client.get("/world/render", params={"view": "topographic"})
    assert resp.status_code == 400


def test_render_out_of_range_dimensions_return_400(client):
    _post_generate(client, json={"seed": 4, "num_plates": 6})
    assert client.get("/world/render", params={"width": 0, "height": 100}).status_code == 400
    assert client.get("/world/render", params={"width": 100, "height": 0}).status_code == 400
    assert client.get("/world/render", params={"width": 100_000, "height": 100}).status_code == 400


def test_generate_replaces_previous_world(client):
    _post_generate(client, json={"seed": 5, "num_plates": 6})
    resp = _post_generate(client, json={"seed": 6, "num_plates": 8})
    assert resp.json()["num_plates"] == 8


def test_generate_without_num_plates_picks_a_plausible_count(client):
    resp = _post_generate(client, json={"seed": 7})
    assert resp.status_code == 200
    assert MIN_AUTO_PLATES <= resp.json()["num_plates"] <= MAX_AUTO_PLATES


def test_render_omitted_rotation_matches_explicit_identity(client):
    _post_generate(client, json={"seed": 8, "num_plates": 8})
    identity = "1,0,0,0,1,0,0,0,1"
    default_body = client.get("/world/render", params={"width": 300, "height": 200}).json()
    explicit_body = client.get("/world/render", params={"width": 300, "height": 200, "rotation": identity}).json()
    assert default_body["image_base64"] == explicit_body["image_base64"]


def test_render_nontrivial_rotation_changes_the_image(client):
    _post_generate(client, json={"seed": 8, "num_plates": 8})
    identity = "1,0,0,0,1,0,0,0,1"
    rotated = "0,0,1,0,1,0,-1,0,0"  # a valid 90-degree rotation matrix
    identity_body = client.get("/world/render", params={"width": 300, "height": 200, "rotation": identity}).json()
    rotated_body = client.get("/world/render", params={"width": 300, "height": 200, "rotation": rotated}).json()
    assert identity_body["image_base64"] != rotated_body["image_base64"]


def test_render_malformed_rotation_returns_400(client):
    _post_generate(client, json={"seed": 8, "num_plates": 8})
    assert client.get("/world/render", params={"rotation": "1,0,0,0,1,0,0,0"}).status_code == 400  # only 8 values
    assert client.get("/world/render", params={"rotation": "not,a,valid,rotation,matrix,at,all,here,either"}).status_code == 400
    assert client.get("/world/render", params={"rotation": "1,0,0,0,1,0,0,0,nan"}).status_code == 400


def test_plates_endpoint_matches_directly_generated_plate_state(client):
    resp = _post_generate(client, json={"seed": 11, "num_plates": 9, "continental_fraction": 0.5})
    assert resp.status_code == 200

    ground_truth = {p.plate_id: p for p in generate_plates(seed=11, num_plates=9, continental_fraction=0.5)}
    body = client.get("/world/plates").json()
    assert len(body["plates"]) == len(ground_truth)

    for entry in body["plates"]:
        truth = ground_truth[entry["plate_id"]]
        assert entry["crust_type"] == truth.crust_type
        assert entry["num_points"] == truth.node_count()
        assert len(entry["outline"]) == len(truth.outline_world())
        assert len(entry["points"]) == truth.node_count()
        if truth.node_count() > 0:
            assert entry["bounding_ellipse"] is not None
            assert entry["bounding_ellipse"]["diameter_a_km"] >= entry["bounding_ellipse"]["diameter_b_km"] >= 0.0
        else:
            assert entry["bounding_ellipse"] is None


def test_plates_endpoint_reports_motion_shape_and_overlap_diagnostics(client):
    from app import mantle

    _post_generate(client, json={"seed": 11, "num_plates": 9, "continental_fraction": 0.5})
    plates = client.get("/world/plates").json()["plates"]

    for entry in plates:
        assert entry["speed_cm_per_yr"] >= 0.0
        # A freshly generated world is clamped into [MIN, MAX]; `at_max_rate` must agree with
        # the reported speed.
        assert entry["at_max_rate"] == (entry["speed_cm_per_yr"] >= mantle.rad_per_yr_to_cm_per_yr(mantle.MAX_PLATE_RATE) - 1e-9)
        assert 0.0 <= entry["submerged_fraction"] <= 1.0
        assert entry["age_years"] >= 0.0
        # A plate never overlaps itself. (Overlap is read by containment of cell centres, so it
        # can be one-sided where two lattices meet: A's centre inside B, but no centre of B's
        # inside A.)
        for over in entry["overlaps"]:
            assert over["plate_id"] != entry["plate_id"]
            assert 0.0 < over["fraction"] <= 1.0



def test_plate_at_returns_the_owning_plate_id(client):
    _post_generate(client, json={"seed": 12, "num_plates": 8})
    plates = client.get("/world/plates").json()["plates"]
    target = next(p for p in plates if p["num_points"] > 0)
    x, y, z = target["points"][0]
    lat_deg = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    lon_deg = math.degrees(math.atan2(y, x))

    resp = client.get("/world/plate_at", params={"lat_deg": lat_deg, "lon_deg": lon_deg})
    assert resp.status_code == 200
    assert resp.json()["plate_id"] == target["plate_id"]


def test_plate_at_rejects_non_finite_query(client):
    _post_generate(client, json={"seed": 12, "num_plates": 8})
    assert client.get("/world/plate_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/plate_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_sample_at_returns_a_full_point_report(client):
    _post_generate(client, json={"seed": 12, "num_plates": 8})
    plates = client.get("/world/plates").json()["plates"]
    target = next(p for p in plates if p["num_points"] > 0)
    x, y, z = target["points"][0]
    lat_deg = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    lon_deg = math.degrees(math.atan2(y, x))

    resp = client.get("/world/sample_at", params={"lat_deg": lat_deg, "lon_deg": lon_deg})
    assert resp.status_code == 200
    body = resp.json()
    assert set(body) == {
        "lat_deg", "lon_deg", "elevation_m", "is_ocean", "biome_id", "biome",
        "temperature_c", "precipitation_mm", "plate_id",
    }
    assert body["plate_id"] == target["plate_id"]
    assert isinstance(body["biome"], str) and body["biome"]
    assert math.isfinite(body["elevation_m"])
    assert math.isfinite(body["temperature_c"])
    assert body["precipitation_mm"] >= 0


def test_sample_at_rejects_non_finite_query(client):
    _post_generate(client, json={"seed": 12, "num_plates": 8})
    assert client.get("/world/sample_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/sample_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_node_at_returns_the_owning_node_phi_theta_and_creation_year(client):
    import numpy as np

    _post_generate(client, json={"seed": 12, "num_plates": 8})
    plates = client.get("/world/plates").json()["plates"]
    target = next(p for p in plates if p["num_points"] > 0)
    x, y, z = target["points"][0]
    lat_deg = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    lon_deg = math.degrees(math.atan2(y, x))

    resp = client.get("/world/node_at", params={"lat_deg": lat_deg, "lon_deg": lon_deg})
    assert resp.status_code == 200
    body = resp.json()
    assert set(body) == {"lat_deg", "lon_deg", "node", "removed"}
    assert body["removed"] is None  # nothing has been removed yet on a freshly generated world

    node = body["node"]
    assert node["plate_id"] == target["plate_id"]
    assert math.isfinite(node["phi"]) and math.isfinite(node["theta"])
    # A freshly generated world's nodes all predate creation-time tracking (see
    # the node_created_years field's own -1.0 sentinel).
    assert node["node_created_years"] == -1.0

    # Round-trips through the owning plate's own frame back to (approximately) the query point.
    world_xyz = np.array([
        math.cos(node["phi"]) * math.cos(node["theta"]),
        math.cos(node["phi"]) * math.sin(node["theta"]),
        math.sin(node["phi"]),
    ])
    assert abs(np.linalg.norm(world_xyz) - 1.0) < 1e-9


def test_node_at_rejects_non_finite_query(client):
    _post_generate(client, json={"seed": 12, "num_plates": 8})
    assert client.get("/world/node_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/node_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_node_at_reports_a_nearby_removed_point(client):
    import numpy as np
    from app import geometry

    _post_generate(client, json={"seed": 12, "num_plates": 8})
    world = _default_session.world
    query_lat_deg, query_lon_deg = 10.0, 20.0
    query_xyz = geometry.latlon_to_xyz(math.radians(query_lat_deg), math.radians(query_lon_deg))
    world.elapsed_years = 7_000_000.0
    world.record_removed_points(np.array([query_xyz]), plate_id=99)

    resp = client.get("/world/node_at", params={"lat_deg": query_lat_deg, "lon_deg": query_lon_deg})
    assert resp.status_code == 200
    removed = resp.json()["removed"]
    assert removed is not None
    assert removed["plate_id"] == 99
    assert removed["removed_years"] == 7_000_000.0
    assert removed["distance_rad"] < 1e-6

    # Far from any removed point, no false match.
    far_resp = client.get("/world/node_at", params={"lat_deg": -query_lat_deg, "lon_deg": query_lon_deg + 180})
    assert far_resp.json()["removed"] is None


def test_elevation_point_at_rejects_non_finite_query(client):
    _post_generate(client, json={"seed": 12, "num_plates": 8})
    assert client.get("/world/elevation_point_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/elevation_point_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_elevation_point_at_returns_the_nearest_node(client):
    _post_generate(client, json={"seed": 12, "num_plates": 8})
    plates_data = client.get("/world/plates").json()["plates"]
    target = next(p for p in plates_data if p["num_points"] > 0)
    x, y, z = target["points"][0]
    lat_deg = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    lon_deg = math.degrees(math.atan2(y, x))

    resp = client.get("/world/elevation_point_at", params={"lat_deg": lat_deg, "lon_deg": lon_deg})

    assert resp.status_code == 200
    body = resp.json()
    assert set(body) == {"plate_id", "point", "point_xyz"}
    assert body["plate_id"] == target["plate_id"]
    px, py, pz = body["point_xyz"]
    assert math.isclose(px * px + py * py + pz * pz, 1.0, abs_tol=1e-6)
    assert math.isfinite(body["point"]["phi"])
    assert math.isfinite(body["point"]["theta"])


def test_surface_node_lookup_uses_storage_neutral_node_id(client):
    _post_generate(client, json={"seed": 12, "num_plates": 8})
    world = _default_session.world
    plate = next(plate for plate in world.plates if plate.node_count() > 0)
    node_id = plate.surface_nodes().node_ids[0]

    response = client.get(
        "/world/surface_node",
        params={"plate_id": plate.plate_id, "node_id_hi": str(int(node_id[0])), "node_id_lo": str(int(node_id[1]))},
    )

    assert response.status_code == 200
    body = response.json()
    assert body["plate_id"] == plate.plate_id
    assert body["node_id"] == [str(int(node_id[0])), str(int(node_id[1]))]
    assert len(body["world_xyz"]) == len(body["local_xyz"]) == 3
    assert body["area_m2"] > 0.0
    assert body["area_is_exact"] is True

    invalid = client.get(
        "/world/surface_node",
        params={"plate_id": plate.plate_id, "node_id_hi": str(2**64), "node_id_lo": "0"},
    )
    assert invalid.status_code == 400


def test_rivers_and_river_at_are_empty_before_the_first_step(client):
    # hydrology_cache is None until erosion.py runs once (see World.hydrology_cache) --
    # /world/rivers should degrade to an empty list rather than erroring, same spirit as
    # /world/plates always having *something* right after generate.
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    resp = client.get("/world/rivers")
    assert resp.status_code == 200
    assert resp.json()["rivers"] == []
    assert client.get("/world/river_at", params={"lat_deg": 0, "lon_deg": 0}).json()["river_id"] is None


def test_river_at_rejects_non_finite_query(client):
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    assert client.get("/world/river_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/river_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_lakes_and_lake_at_are_empty_before_the_first_step(client):
    # hydrology_cache is None until erosion.py runs once, same degrade-to-empty contract
    # /world/rivers already has -- see test_rivers_and_river_at_are_empty_before_the_first_step.
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    resp = client.get("/world/lakes")
    assert resp.status_code == 200
    assert resp.json()["lakes"] == []
    at = client.get("/world/lake_at", params={"lat_deg": 0, "lon_deg": 0}).json()
    assert at == {"kind": "no_basin", "basin": None}


def test_lake_summary_floor_is_the_lowest_ground_under_sub_cell_relief():
    # Issue #344: a one-cell lake wet only in its cell's low ground (relief 20 m about a
    # 180.17 m mean) has its surface below the cell mean. The reported floor is that low
    # ground, the level the lake model holds no water at, so it's never above the water.
    import types

    import numpy as np

    from app.hydroclimate import hydrology, lakes

    elevation = np.array([180.17, 185.0, 186.0])
    hyps = lakes.CellHypsometry(elevation, np.full(3, 1.0e9), np.array([20.0, 0.0, 0.0]), 1.0e9)
    members = np.array([0], dtype=np.int64)
    water = 179.32
    lake_depth = np.zeros(3)
    lake_depth[0] = hyps.mean_depth(water, members)[0]
    assert lake_depth[0] > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    lake = lakes.Lake(
        lake_id=0, members=members, floor_elevation=180.17, min_depth=180.17, max_depth=180.34,
        outlet_node_idx=0, current_water_elevation=water, sink_node_idx=0,
    )
    fields = types.SimpleNamespace(elevation=elevation, lake_depth=lake_depth, points=np.eye(3), hypsometry=hyps)

    out = main._lake_basin_summary(fields, lake, 0, [], {}, is_lake=True)

    assert out["water_elevation_m"] == pytest.approx(water)
    assert out["floor_elevation_m"] == pytest.approx(160.17)
    assert out["floor_elevation_m"] <= out["water_elevation_m"] <= out["outlet_elevation_m"]


def test_lake_at_rejects_non_finite_query(client):
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    assert client.get("/world/lake_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/lake_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_faults_is_empty_right_after_generate(client):
    # No fault spawns before the first step, same degrade-to-empty contract /world/rivers has.
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    resp = client.get("/world/faults")
    assert resp.status_code == 200
    assert resp.json()["faults"] == []
    assert resp.json()["fault_systems"] == []
    assert client.get("/world/fault_at", params={"lat_deg": 0, "lon_deg": 0}).json()["fault_id"] is None


def test_fault_at_rejects_non_finite_query(client):
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    assert client.get("/world/fault_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/fault_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_faults_returns_well_formed_entries_after_stepping(client):
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    for _ in range(20):
        _post_step(client, json={"years": 1_000_000})
    body = client.get("/world/faults").json()
    assert body["faults"], "expected some faults after 20 Myr"
    f = body["faults"][0]
    assert f["kind"] in {"normal", "reverse", "strike_slip"}
    assert len(f["trace"]) >= 2 and all(len(p) == 3 for p in f["trace"])
    assert f["length_km"] > 0.0
    assert "system_id" in f

    # Boundary faults (faults.generate_boundary_faults) come through the same list, flagged
    # `boundary` and with a null lifespan (never aged / retired).
    boundary_faults = [x for x in body["faults"] if x["boundary"]]
    assert len(boundary_faults) > 20, "boundary faults should line the plate boundaries"
    assert all(x["lifespan_myr"] is None for x in boundary_faults)
    assert all(not x["boundary"] or x["kind"] in {"normal", "reverse", "strike_slip"} for x in body["faults"])

    # Fault systems: well-formed, and every strand's system_id points at a real system.
    systems = body["fault_systems"]
    system_ids = {s["system_id"] for s in systems}
    for s in systems:
        assert s["kind"] in {"normal", "reverse", "strike_slip"}
        assert len(s["trace"]) >= 2 and all(len(p) == 3 for p in s["trace"])
        assert s["length_km"] > 0.0
    for strand in body["faults"]:
        if strand["system_id"] is not None:
            assert strand["system_id"] in system_ids
    # a click on the first fault's own midpoint hit-tests back to it
    mid = f["trace"][len(f["trace"]) // 2]
    lat = math.degrees(math.asin(max(-1.0, min(1.0, mid[2]))))
    lon = math.degrees(math.atan2(mid[1], mid[0]))
    assert client.get("/world/fault_at", params={"lat_deg": lat, "lon_deg": lon}).json()["fault_id"] == f["fault_id"]


def test_volcanoes_returns_well_formed_entries(client):
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    # Empty (or at least well-formed) right after generate; degrade-to-empty like /world/faults.
    body = client.get("/world/volcanoes").json()
    assert body["volcanoes"] == []
    for _ in range(30):
        _post_step(client, json={"years": 5_000_000})
    body = client.get("/world/volcanoes").json()
    assert body["volcanoes"], "expected some volcano nodes after 150 Myr of subduction"
    for v in body["volcanoes"]:
        assert len(v["position"]) == 3 and all(isinstance(c, float) for c in v["position"])
        assert isinstance(v["active"], bool)


def test_stranded_basins_is_empty_before_the_first_step(client):
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    resp = client.get("/world/stranded_basins")
    assert resp.status_code == 200
    body = resp.json()
    assert body["stranded_basins"] == []
    assert body["sea_level_m"] == 0.0


def test_stranded_basins_returns_well_formed_entries_after_stepping(client):
    _post_generate(client, json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    _post_step(client, json={"years": 3_000_000})
    body = client.get("/world/stranded_basins").json()
    for basin in body["stranded_basins"]:
        assert basin["depth_below_sea_level_m"] > 0.0
        assert basin["floor_elevation_m"] < body["sea_level_m"]
        assert basin["catchment_node_count"] >= 1
        assert basin["steps_seen"] >= 1
        assert basin["persisted_years"] >= 0.0


def test_stranded_basin_summary_wire_shape():
    from app.main import _stranded_basin_summary
    from app.hydroclimate.stranded_basins import StrandedBasin

    basin = StrandedBasin(
        floor_elevation_m=-1770.58,
        depth_below_sea_level_m=1770.58,
        catchment_node_count=435,
        flooded_node_count=412,
        water_elevation_m=-1750.2,
        centroid_xyz=(0.6234567, 0.55, -0.21),
        centroid_lat_deg=-12.3,
        centroid_lon_deg=45.6,
        floor_xyz=(0.61, 0.56, -0.22),
        first_seen_years=72_700_000.0,
        persisted_years=12_400_000.0,
        steps_seen=124,
    )
    out = _stranded_basin_summary(basin)
    assert out["floor_elevation_m"] == -1770.6  # rounded to 1 dp
    assert out["catchment_node_count"] == 435 and out["flooded_node_count"] == 412
    assert len(out["centroid_xyz"]) == 3 and out["centroid_xyz"][0] == 0.623457
    assert out["persisted_years"] == 12_400_000.0 and out["steps_seen"] == 124


def test_stats_returns_expected_shape(client):
    _post_generate(client, json={"seed": 13, "num_plates": 8, "continental_fraction": 0.5})
    resp = client.get("/world/stats")
    assert resp.status_code == 200
    body = resp.json()

    assert body["elapsed_years"] == 0.0
    assert 0.0 <= body["land_fraction"] <= 1.0
    assert 0.0 <= body["ocean_fraction"] <= 1.0
    assert math.isclose(body["land_fraction"] + body["ocean_fraction"], 1.0, abs_tol=1e-9)
    assert body["elevation_min_m"] <= body["elevation_mean_m"] <= body["elevation_max_m"]
    # A seed/continental_fraction generating both crust types should have both land and
    # ocean grid cells, so none of these should fall back to their None (empty-mask) case.
    for key in (
        "land_temperature_min_c", "land_temperature_max_c", "land_temperature_mean_c",
        "air_temperature_min_c", "air_temperature_max_c", "air_temperature_mean_c",
        "ocean_temperature_min_c", "ocean_temperature_max_c", "ocean_temperature_mean_c",
    ):
        assert body[key] is not None
    assert body["precipitation_min_mm"] <= body["precipitation_mean_mm"] <= body["precipitation_max_mm"]
    # Both biome breakdowns are present and each sums to ~1 over its own domain (land / ocean).
    assert math.isclose(sum(body["biome_land_fraction"].values()), 1.0, abs_tol=1e-6)
    assert math.isclose(sum(body["biome_ocean_fraction"].values()), 1.0, abs_tol=1e-6)


def test_stats_history_records_generate_and_each_step_and_survives_save_load(client):
    # The whole point of World.stats_history (see world.py/persistence.py): a save/load round
    # trip should keep the Stats panel's history charts, not silently drop them.
    _post_generate(client, json={"seed": 13, "num_plates": 8, "continental_fraction": 0.5})
    history = client.get("/world/stats_history").json()
    assert history["format"] == "columnar-v1"
    assert history["length"] == 1
    assert history["columns"]["elapsed_years"] == [0.0]
    land_biomes = history["columns"]["biome_land_fraction"]
    assert land_biomes["states"] == [2]
    assert land_biomes["keys"]
    assert all(len(values) == 1 for values in land_biomes["columns"].values())

    _post_step(client, json={"years": 1_000_000})
    history = client.get("/world/stats_history").json()
    assert history["length"] == 2
    assert history["columns"]["elapsed_years"] == [0.0, 1_000_000.0]

    saved = client.get("/world/save").content
    _post_generate(client, json={"seed": 99, "num_plates": 4})  # replace with something else
    assert client.get("/world/stats_history").json()["columns"]["elapsed_years"] == [0.0]

    client.post("/world/load", content=saved)
    reloaded_history = client.get("/world/stats_history").json()
    assert reloaded_history["columns"]["elapsed_years"] == [0.0, 1_000_000.0]
    assert reloaded_history == history


def test_save_then_load_round_trips_the_exact_world_state(client):
    _post_generate(client, json={"seed": 5, "num_plates": 8, "continental_fraction": 0.5})
    _post_step(client, json={"years": 2_000_000})
    generated_summary = client.get("/world/plates").json()

    saved = client.get("/world/save")
    assert saved.status_code == 200
    assert saved.headers["content-type"] == "application/octet-stream"
    assert "attachment" in saved.headers["content-disposition"]

    # Advance further so the live world visibly differs from the snapshot just saved.
    _post_step(client, json={"years": 2_000_000})
    assert client.get("/world/render").json()["elapsed_years"] == 4_000_000.0

    loaded = client.post("/world/load", content=saved.content)
    assert loaded.status_code == 200
    body = loaded.json()
    assert body["seed"] == 5
    assert body["elapsed_years"] == 2_000_000.0
    assert client.get("/world/plates").json() == generated_summary


def test_load_with_malformed_bytes_returns_400(client):
    resp = client.post("/world/load", content=b"definitely not a saved world")
    assert resp.status_code == 400
    assert "not a readable mantle-bloom save" in resp.json()["detail"]


def test_load_reports_a_save_from_a_newer_build_as_unsupported(client):
    import pickle

    from app import persistence

    newer = pickle.dumps({"format": persistence.SAVE_FORMAT, "version": persistence.SAVE_FORMAT_VERSION + 1, "world": None})
    resp = client.post("/world/load", content=newer)
    assert resp.status_code == 400
    assert "unsupported save format version" in resp.json()["detail"]


def test_load_converts_a_line_save_to_quads(client):
    # docs/save-compatibility.md: a line-backed save converts to sparse quads on load, one way.
    from app.elevation_lines import line_spacing_rad
    from app.world import generate_world

    from .legacy_lines import line_world_like, retired_class_names

    lines = line_world_like(generate_world(seed=5, num_plates=6, node_density=1.0), line_spacing_rad(1.0))
    with retired_class_names():
        saved = persistence.save_world_bytes(lines)

    converted = client.post("/world/load", content=saved).json()
    assert converted["surface"] == "quad"
    assert converted["surface_conversion"]["from"] == "lines"
    assert converted["num_plates"] == 6
    assert client.get("/world/render").status_code == 200
    resaved = client.get("/world/save").content
    assert client.post("/world/load", content=resaved).json()["surface"] == "quad"


def test_animate_advances_the_world_and_streams_progress_then_an_mp4(client):
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = client.post(
        "/world/animate",
        json={"projection": "eckert4", "view": "elevation", "width": 200, "height": 110, "step_years": 1_000_000, "num_frames": 3},
    )
    assert resp.status_code == 200

    messages = [json.loads(line) for line in resp.text.splitlines() if line.strip()]
    progress = [m for m in messages if m["type"] == "progress"]
    assert [m["frame"] for m in progress] == [1, 2, 3]
    assert all(m["total"] == 3 for m in progress)
    # Each progress tick carries that frame's own PNG so the client can paint it live.
    assert all(Image.open(io.BytesIO(base64.b64decode(m["image_base64"]))).size == (200, 110) for m in progress)
    # ...and its own elapsed_years, so the client can drive a live counter off the stream
    # instead of it sitting frozen until the final "done" message.
    assert [m["elapsed_years"] for m in progress] == [0.0, 1_000_000.0, 2_000_000.0]

    done = messages[-1]
    assert done["type"] == "done"
    assert done["mime"] == "video/mp4"
    assert done["elapsed_years"] == 2_000_000.0  # (num_frames - 1) * steps_per_frame * step_years

    video = base64.b64decode(done["video_base64"])
    assert video[4:12] == b"ftypisom"  # an MP4 container
    with av.open(io.BytesIO(video)) as container:
        assert sum(1 for _ in container.decode(video=0)) == 3


def test_animate_disconnect_mid_stream_still_publishes_completed_frames(client):
    # Round-4 Finding 3: animate used to defer its only _world_changed to normal
    # completion, so a disconnect left the feed revision, render cache and index
    # describing the pre-animation world even though the completed frames were
    # already checkpointed. Disconnect after the first frame's progress line;
    # the feed must still advance past the pre-animation revision.
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    before = _default_session.feed.revision()
    chunk_count = _post_streaming_then_disconnect(
        f"/worlds/{_default_session.world_id}/animate",
        {"projection": "eckert4", "view": "elevation", "width": 64, "height": 64,
         "step_years": 100_000, "num_frames": 2},
        disconnect_after_chunks=1,
    )
    assert chunk_count == 1  # confirms the scenario really did cut the stream short
    assert _wait_until(lambda: _default_session.feed.revision() != before)
    assert _default_session.feed.revision() != before


def test_animate_progress_carries_a_live_stats_snapshot(client):
    # GitHub issue #155: the Stats panel used to sit frozen for the whole recording, since
    # GET /world/stats blocks on the world lock that /world/animate holds throughout. Each
    # progress line should instead carry the latest World.stats_history snapshot so the
    # frontend can update it live off the stream, same as elapsed_years/frame already do.
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = client.post(
        "/world/animate",
        json={"projection": "eckert4", "view": "elevation", "width": 64, "height": 64, "step_years": 100_000, "num_frames": 3},
    )
    assert resp.status_code == 200

    # Steps well inside merge_split.SPLIT_MIN_AGE_YEARS, so no plate can rift mid-recording.
    progress = [json.loads(line) for line in resp.text.splitlines() if line.strip() and json.loads(line)["type"] == "progress"]
    assert [m["stats"]["elapsed_years"] for m in progress] == [0.0, 100_000.0, 200_000.0]
    assert all(m["stats"]["plate_count"] == 6 for m in progress)


def test_animate_steps_per_frame_runs_every_step_but_renders_only_the_last(client):
    # 3 frames * 4 steps_per_frame = 1 (frame 0, unstepped) + 2*4 = 8 real step_world calls
    # total, even though only 3 frames are ever rendered/encoded -- see
    # render_image.stream_animation_mp4's own docstring for why this isn't the same as one
    # bigger step_world(world, steps_per_frame * step_years) call per frame.

    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = client.post(
        "/world/animate",
        json={"width": 64, "height": 64, "step_years": 250_000, "steps_per_frame": 4, "num_frames": 3},
    )
    assert resp.status_code == 200

    messages = [json.loads(line) for line in resp.text.splitlines() if line.strip()]
    progress = [m for m in messages if m["type"] == "progress"]
    # Each rendered frame is still steps_per_frame * step_years further along than the last.
    assert [m["elapsed_years"] for m in progress] == [0.0, 1_000_000.0, 2_000_000.0]

    world = _default_session.world
    assert world.steps_taken == 8  # (num_frames - 1) * steps_per_frame real step_world calls
    assert world.elapsed_years == 2_000_000.0

    done = messages[-1]
    with av.open(io.BytesIO(base64.b64decode(done["video_base64"]))) as container:
        assert sum(1 for _ in container.decode(video=0)) == 3  # still only 3 frames encoded


def test_animate_rejects_zero_or_negative_steps_per_frame(client):
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = client.post("/world/animate", json={"step_years": 1_000_000, "steps_per_frame": 0, "num_frames": 2})
    assert resp.status_code == 400
    resp = client.post("/world/animate", json={"step_years": 1_000_000, "steps_per_frame": -1, "num_frames": 2})
    assert resp.status_code == 400


def test_animate_rejects_out_of_range_frame_counts(client):
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = client.post("/world/animate", json={"step_years": 1_000_000, "num_frames": 0})
    assert resp.status_code == 400
    resp = client.post("/world/animate", json={"step_years": 1_000_000, "num_frames": 10_000})
    assert resp.status_code == 400


def test_animate_rejects_unknown_view(client):
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = client.post("/world/animate", json={"view": "not-a-real-view", "step_years": 1_000_000, "num_frames": 2})
    assert resp.status_code == 400


def test_animate_stop_is_a_no_op_when_nothing_is_running(client):
    # No animation in flight (not even a world yet) -- the "Stop" toolbar button should never
    # itself be the thing that errors.
    resp = client.post("/world/animate/stop")
    assert resp.status_code == 200

    _post_generate(client, json={"seed": 9, "num_plates": 6})
    resp = client.post("/world/animate/stop")
    assert resp.status_code == 200


def test_animate_clears_a_stale_stop_signal_from_a_previous_call(client):
    # A stop request that arrives after its run already finished shouldn't truncate the
    # *next* run -- animate() clears the signal itself at the start of each call, so this
    # doesn't depend on ordering against a previous test.
    _post_generate(client, json={"seed": 9, "num_plates": 6})
    client.post("/world/animate/stop")  # simulate a stray/late stop signal
    resp = client.post(
        "/world/animate",
        json={"width": 64, "height": 64, "step_years": 1_000_000, "num_frames": 3},
    )
    messages = [json.loads(line) for line in resp.text.splitlines() if line.strip()]
    progress = [m for m in messages if m["type"] == "progress"]
    assert [m["frame"] for m in progress] == [1, 2, 3]  # ran to completion, not truncated
    assert messages[-1]["stopped_early"] is False


def test_export_hexgrid_returns_the_requested_tile_count(client):
    _post_generate(client, json={"seed": 11, "num_plates": 8})
    resp = client.post("/world/export_hexgrid", json={"frequency": 8})
    assert resp.status_code == 200
    body = resp.json()
    assert body["num_tiles"] == 642 == len(body["tiles"])
    assert sum(1 for t in body["tiles"] if t["is_pentagon"]) == 12


def test_export_hexgrid_rejects_unknown_frequency(client):
    _post_generate(client, json={"seed": 11, "num_plates": 8})
    resp = client.post("/world/export_hexgrid", json={"frequency": 5})
    assert resp.status_code == 400


def test_generate_with_unknown_fluid_density_returns_400(client):
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "fluid_density": 3.0})
    assert resp.status_code == 400


def test_generate_with_the_finest_allowed_fluid_density(client):
    # fluid_density is capped at 2.0 ("High") -- unlike climate_density, which allows 4.0 --
    # since Ocean/Atmospheric Fluid Dynamics now runs every step, not just opt-in.
    resp = _post_generate(client, json={"seed": 1, "num_plates": 6, "fluid_density": 2.0})
    assert resp.status_code == 200


def test_retired_cfd_sediment_views_are_not_available(client):
    # The ocean CFD (and its sediment concentration/deposition views) was retired -- see
    # render_image.py's VIEWS, which no longer lists them at all.
    _post_generate(client, json={"seed": 12, "num_plates": 6, "fluid_density": 0.5})
    for view in ["oceanCfdSediment", "oceanCfdDeposition"]:
        assert client.get("/world/render", params={"view": view}).status_code == 400


def test_climate_views_render_natively_off_the_healpix_grid(client):
    # Every CLIMATE_VIEWS member renders directly off the world's own HEALPix grid (see
    # render_image._render_climate_view) rather than resampling down to an equirectangular
    # grid first.
    from app import render_image

    _post_generate(client, json={"seed": 21, "num_plates": 6, "node_density": 0.5, "climate_density": 0.5, "fluid_density": 0.5})
    _post_step(client, json={"years": 2_000_000})
    for view in render_image.CLIMATE_VIEWS:
        resp = client.get("/world/render", params={"view": view, "width": 100, "height": 50})
        assert resp.status_code == 200
        assert len(resp.json()["image_base64"]) > 0


def test_wind_and_ocean_currents_views_are_also_always_renderable(client):
    # `wind` draws the CFD-sourced wind; `oceanCurrents` draws climate.py's diagnostic
    # currents (see climate.py's own module docstring) -- both always renderable regardless of
    # what else the world is doing.
    _post_generate(client, json={"seed": 12, "num_plates": 6, "fluid_density": 0.5})
    assert client.get("/world/render", params={"view": "wind"}).status_code == 200
    assert client.get("/world/render", params={"view": "oceanCurrents"}).status_code == 200


def test_step_advances_atmosphere_cfd_by_its_own_fixed_seconds(client):
    # Regardless of the tectonic `years` requested, the atmospheric wind solve always advances
    # by its own fixed real-time increment per step (atmosphere_cfd.SECONDS_PER_TECTONIC_STEP)
    # -- reaches into main.py's internal `_state` since this is not surfaced through any HTTP
    # response.
    from app.hydroclimate import atmosphere_cfd

    _post_generate(client, json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    client.post("/world/controls", json={"wind_model": "cfd"})  # the CFD solve only advances under "cfd" (default is "diagnostic")
    world = _default_session.world
    assert world.atmosphere_cfd_state.elapsed_seconds == 0.0

    resp = _post_step(client, json={"years": 1_000_000})
    assert resp.status_code == 200
    assert world.atmosphere_cfd_state.elapsed_seconds == pytest.approx(atmosphere_cfd.SECONDS_PER_TECTONIC_STEP)


def test_step_does_not_advance_fluid_dynamics_when_climate_biomes_paused(client):

    _post_generate(client, json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    client.post("/world/controls", json={"simulate_climate_biomes": False})
    world = _default_session.world

    resp = _post_step(client, json={"years": 1_000_000})
    assert resp.status_code == 200
    assert world.atmosphere_cfd_state.elapsed_seconds == 0.0


def test_controls_ice_age_frequency_round_trip(client):

    _post_generate(client, json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    world = _default_session.world
    assert world.ice_age_period_years == 0.0
    assert client.post("/world/controls", json={}).json()["ice_age_period_years"] == 0.0

    resp = client.post("/world/controls", json={"ice_age_period_years": 300_000})
    assert resp.status_code == 200
    assert resp.json()["ice_age_period_years"] == 300_000
    assert world.ice_age_period_years == 300_000

    # A negative period is clamped to 0 (disabled), not an error.
    resp = client.post("/world/controls", json={"ice_age_period_years": -5})
    assert resp.status_code == 200
    assert world.ice_age_period_years == 0.0


def test_controls_debug_diagnostics_round_trip(client):

    _post_generate(client, json={"seed": 12, "num_plates": 6})
    world = _default_session.world
    assert world.debug_diagnostics is False
    assert client.post("/world/controls", json={}).json()["debug_diagnostics"] is False

    resp = client.post("/world/controls", json={"debug_diagnostics": True})
    assert resp.status_code == 200
    assert resp.json()["debug_diagnostics"] is True
    assert world.debug_diagnostics is True


def test_controls_wind_model_toggle_and_validation(client):

    _post_generate(client, json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    world = _default_session.world
    assert world.wind_model == "diagnostic"

    resp = client.post("/world/controls", json={"wind_model": "cfd"})
    assert resp.status_code == 200
    assert resp.json()["wind_model"] == "cfd"
    assert world.wind_model == "cfd"

    resp = client.post("/world/controls", json={"wind_model": "diagnostic"})
    assert resp.status_code == 200
    assert resp.json()["wind_model"] == "diagnostic"
    assert world.wind_model == "diagnostic"

    # Diagnostic mode leaves the CFD solve unrun on a step.
    _post_step(client, json={"years": 1_000_000})
    assert world.atmosphere_cfd_state.elapsed_seconds == 0.0

    # A climate render still succeeds against the diagnostic fields.
    assert client.get("/world/render", params={"view": "wind"}).status_code == 200
    assert client.get("/world/render", params={"view": "temperature"}).status_code == 200

    assert client.post("/world/controls", json={"wind_model": "nonsense"}).status_code == 400


def test_controls_fault_deformation_mode_toggle_and_validation(client):

    _post_generate(client, json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    world = _default_session.world
    assert world.fault_deformation_mode == "fault"

    for mode in ("boundary", "both", "fault"):
        resp = client.post("/world/controls", json={"fault_deformation_mode": mode})
        assert resp.status_code == 200
        assert resp.json()["fault_deformation_mode"] == mode
        assert world.fault_deformation_mode == mode

    assert client.post("/world/controls", json={"fault_deformation_mode": "nonsense"}).status_code == 400



def test_controls_node_cloud_resample_mode_toggle_and_validation(client):

    _post_generate(client, json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    world = _default_session.world
    assert world.node_cloud_resample_mode == "kdtree"

    for mode in ("healpix", "kdtree"):
        resp = client.post("/world/controls", json={"node_cloud_resample_mode": mode})
        assert resp.status_code == 200
        assert resp.json()["node_cloud_resample_mode"] == mode
        assert world.node_cloud_resample_mode == mode

    assert client.post("/world/controls", json={"node_cloud_resample_mode": "nonsense"}).status_code == 400


def test_earthquakes_endpoint(client):
    _post_generate(client, json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    assert client.get("/world/earthquakes").json()["earthquakes"] == []  # none before a step
    for _ in range(20):
        _post_step(client, json={"years": 1_000_000})
    body = client.get("/world/earthquakes").json()
    for q in body["earthquakes"]:
        assert len(q["epicenter"]) == 3
        assert 3.5 <= q["magnitude"] <= 9.5
        assert q["age_myr"] >= 0.0
        assert q["trigger"] in ("tectonic", "ice_unloading")


def test_controls_tuning_multipliers_round_trip_and_validate(client):
    from app.world import TUNING_MULTIPLIER_FIELDS

    _post_generate(client, json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    world = _default_session.world

    # Default: every knob reads 1.0, and an untouched POST echoes them all.
    body = client.post("/world/controls", json={}).json()
    for name in TUNING_MULTIPLIER_FIELDS:
        assert body[name] == 1.0
        assert getattr(world, name) == 1.0

    resp = client.post("/world/controls", json={"rain_erosion_multiplier": 0.0, "collision_uplift_multiplier": 2.5})
    assert resp.status_code == 200
    assert resp.json()["rain_erosion_multiplier"] == 0.0
    assert resp.json()["collision_uplift_multiplier"] == 2.5
    assert world.rain_erosion_multiplier == 0.0
    assert world.collision_uplift_multiplier == 2.5
    assert world.river_erosion_multiplier == 1.0  # untouched knob unchanged

    assert client.post("/world/controls", json={"volcanism_multiplier": -0.5}).status_code == 400
    assert world.volcanism_multiplier == 1.0  # rejected write did not land


# --- Run status and checkpoint (#375) -------------------------------------------------------


@pytest.fixture(scope="module")
def saved_world() -> bytes:
    from app.world import generate_world

    return persistence.save_world_bytes(generate_world(seed=5, num_plates=2, node_density=0.5))


@pytest.fixture
def store(client, tmp_path, monkeypatch):
    """Checkpointing into tmp_path, and a fake step that only advances the clock. Teardown
    stops the run loop and waits out any step it has under way."""
    from app import checkpoint

    monkeypatch.setenv(checkpoint.CHECKPOINT_DIR_ENV, str(tmp_path))
    monkeypatch.setattr(main, "_checkpoint_root", checkpoint.CheckpointStore(tmp_path))
    monkeypatch.setattr(main, "_checkpoint_error", None)
    _default_session.checkpoint_store = main._checkpoint_root.for_world(_default_session.world_id)
    main._save_index()

    def fake_step(world, years):
        world.elapsed_years += years
        yield 1.0

    monkeypatch.setattr(main, "step_world_progress", fake_step)
    yield _default_session.checkpoint_store
    _default_session.runner.set_status(checkpoint.RunStatus(), persist=False)
    with _default_session.lock:
        pass


def _checkpoint(store, saved_world, running, step_years=1_000.0):
    from app.checkpoint import RunStatus

    world = persistence.load_world_bytes(saved_world)
    world.elapsed_years = 123_000.0
    store.save_world(world)
    store.save_status(RunStatus(running=running, step_years=step_years))


def _saved_elapsed(store) -> float:
    return persistence.load_world_bytes(store.world_path.read_bytes()).elapsed_years


def test_restart_restores_a_stopped_world_at_the_same_point_and_stays_stopped(client, store, saved_world):
    _checkpoint(store, saved_world, running=False)
    main._restore_checkpoint()

    assert client.get("/world/summary").json()["elapsed_years"] == 123_000.0
    assert client.get("/world/run").json()["running"] is False
    time.sleep(0.05)
    assert _default_session.world.elapsed_years == 123_000.0


def test_restart_resumes_a_running_world_from_the_same_point(client, store, saved_world):
    _checkpoint(store, saved_world, running=True, step_years=500.0)
    main._restore_checkpoint()

    run = client.get("/world/run").json()
    assert (run["running"], run["step_years"]) == (True, 500.0)
    # Resumes from the saved point: no catch-up for time spent asleep, then 500 years a step.
    assert _wait_until(lambda: _default_session.world.elapsed_years >= 123_500.0)
    assert (_default_session.world.elapsed_years - 123_000.0) % 500.0 == 0
    # Each completed step is checkpointed.
    assert _wait_until(lambda: _saved_elapsed(store) >= 123_500.0)


def test_play_and_stop_save_the_run_status(client, store, saved_world):
    main._restore_checkpoint()  # nothing saved yet: no world, stopped
    assert client.post("/world/run", json={"running": True, "step_years": 10.0}).status_code == 404
    client.post("/world/load", content=saved_world)

    assert client.post("/world/run", json={"running": True, "step_years": 10.0}).json()["running"] is True
    assert json.loads(store.status_path.read_text())["running"] is True
    start = persistence.load_world_bytes(saved_world).elapsed_years
    assert _wait_until(lambda: _default_session.world.elapsed_years > start)

    body = client.post("/world/run", json={"running": False}).json()
    assert (body["running"], body["step_years"]) == (False, 10.0)
    assert json.loads(store.status_path.read_text()) == {"version": 1, "running": False, "step_years": 10.0}


def test_run_rejects_a_bad_step_size(client, store):
    assert client.post("/world/run", json={"running": False, "step_years": 0}).status_code == 400


def test_a_failed_step_does_not_overwrite_the_checkpoint(client, store, saved_world, monkeypatch):
    _checkpoint(store, saved_world, running=False)
    main._restore_checkpoint()

    def failing_step(world, years):
        world.elapsed_years += years
        yield 0.5
        raise RuntimeError("degenerate geometry")

    monkeypatch.setattr(main, "step_world_progress", failing_step)
    lines = client.post("/world/step", json={"years": 1_000.0}).text.splitlines()
    assert json.loads(lines[-1])["type"] == "error"
    assert _default_session.world.elapsed_years == 124_000.0
    assert _saved_elapsed(store) == 123_000.0


def test_run_poll_wakes_when_another_tab_stops_the_run(client, store, saved_world):
    _checkpoint(store, saved_world, running=False)
    main._restore_checkpoint()
    client.post("/world/run", json={"running": True, "step_years": 10.0})
    result = {}

    def follow():
        # Like the full app: poll past whatever was last seen until the run has stopped.
        run = client.get("/world/run").json()
        while run["running"]:
            run = client.get("/world/run", params={"after": run["world"]["revision"], "run_version": run["run_version"]}).json()
        result.update(run)

    poller = threading.Thread(target=follow)
    started = time.monotonic()
    poller.start()
    time.sleep(0.05)
    client.post("/world/run", json={"running": False})
    poller.join(10)

    # Without the wake, the last poll would sit out the 25 s timeout: a stopped run publishes
    # no new revision.
    assert time.monotonic() - started < 5
    assert result["running"] is False


def test_a_status_that_cannot_be_saved_is_reported(client, store, saved_world, monkeypatch):
    _checkpoint(store, saved_world, running=False)
    main._restore_checkpoint()

    def fails(status):
        raise OSError("volume gone")

    monkeypatch.setattr(_default_session.checkpoint_store, "save_status", fails)
    resp = client.post("/world/run", json={"running": False, "step_years": 20.0})
    assert resp.status_code == 500
    assert client.get("/world/run").json()["step_years"] == 20.0
    monkeypatch.undo()
    _default_session.checkpoint_store = None


def test_health_check_fails_when_the_checkpoint_directory_is_unusable(client, store, tmp_path, monkeypatch):
    # A deploy whose volume can't hold a checkpoint must not pass as healthy (#376).
    from app import checkpoint

    blocker = tmp_path / "not-a-directory"
    blocker.write_text("")
    monkeypatch.setenv(checkpoint.CHECKPOINT_DIR_ENV, str(blocker))
    main._restore_checkpoint()

    resp = client.get("/healthz")
    assert resp.status_code == 503
    assert "not-a-directory" in resp.json()["detail"]


def test_an_animation_that_fails_keeps_the_frames_it_finished(client, store, saved_world, monkeypatch):
    _checkpoint(store, saved_world, running=False)
    main._restore_checkpoint()

    def frames(world, *args, **kwargs):
        world.elapsed_years += 1_000.0
        yield ("progress", 1, 3, b"png", world.elapsed_years)
        world.elapsed_years += 500.0  # a step that fails halfway
        raise RuntimeError("encoder died")

    monkeypatch.setattr(main.render_image, "stream_animation_mp4", frames)
    lines = client.post("/world/animate", json={"step_years": 1_000.0, "num_frames": 3}).text.splitlines()

    assert json.loads(lines[-1])["type"] == "error"
    assert _saved_elapsed(store) == 124_000.0


# --- Replay log (#442) --------------------------------------------------------------------------


@pytest.fixture
def replay_logged(client, saved_world, tmp_path, monkeypatch):
    """A world, and a replay log whose one trigger fires on every step of a clock-only step."""
    from app import replay_log

    def fake_step(world, years):
        world.elapsed_years += years
        yield 1.0

    log = replay_log.ReplayLog(tmp_path, min_steps_between=1, triggers={"always": lambda outcome: {}}, step_progress=fake_step)
    _default_session.replay_log = log
    _default_session.world = persistence.load_world_bytes(saved_world)
    return log


def test_replay_log_is_off_without_a_directory(client, monkeypatch):
    from app import replay_log

    monkeypatch.delenv(replay_log.REPLAY_LOG_DIR_ENV, raising=False)
    _default_session.replay_log = main._make_replay_log(_default_session.world_id)
    assert _default_session.replay_log is None
    assert client.get("/world/replay-logs").json() == {"enabled": False, "directory": None, "records": []}


def test_replay_log_starts_from_the_env(client, tmp_path, monkeypatch):
    from app import replay_log

    monkeypatch.setenv(replay_log.REPLAY_LOG_DIR_ENV, str(tmp_path))
    _default_session.replay_log = main._make_replay_log(_default_session.world_id)
    assert client.get("/world/replay-logs").json() == {
        "enabled": True,
        "directory": str(tmp_path / _default_session.world_id),
        "records": [],
    }


def test_a_logged_step_is_listed_and_shown_in_the_summary(client, replay_logged):
    done = _post_step(client, json={"years": 1_000.0}).json()
    replay_logged.flush()

    [record] = client.get("/world/replay-logs").json()["records"]
    assert record["triggers"] == ["always"]
    assert [e["message"] for e in done["replay_logs"]] == [f"Replay log: always -> {record['name']}"]
    # Not world state: the world's own event log doesn't carry it.
    assert all("Replay log" not in e["message"] for e in done["events"])


def test_a_new_world_clears_the_console_s_replay_log_messages(client, replay_logged):
    _post_step(client, json={"years": 1_000.0})
    replay_logged.flush()
    assert client.get("/world/summary").json()["replay_logs"]
    done = _post_generate(client, json={"seed": 5, "num_plates": 2, "node_density": 0.5}).json()
    assert done["replay_logs"] == []
    # The record itself stays on disk.
    assert len(client.get("/world/replay-logs").json()["records"]) == 1


def test_run_loop_and_animation_steps_are_logged(client, replay_logged):
    assert main._run_step_once(_default_session, 1_000.0, lambda: True)
    main._step_world(_default_session, _default_session.world, 1_000.0)
    replay_logged.flush()
    assert len(client.get("/world/replay-logs").json()["records"]) == 2


@pytest.fixture
def mounted_frontend(tmp_path):
    """A stand-in `vite build` output mounted on the real app, as the Railway image does
    (#376). mount_frontend adds a catch-all route, so the app's routes are put back after."""
    (tmp_path / "assets").mkdir()
    (tmp_path / "index.html").write_text("<!doctype html><title>bundle</title>")
    (tmp_path / "assets" / "index-abc.js").write_text("console.log(1)")
    routes = list(app.router.routes)
    main.mount_frontend(app, tmp_path)
    yield TestClient(app, base_url="https://mantle.example.com", follow_redirects=False)
    app.router.routes[:] = routes


def test_mounted_frontend_serves_spa_routes_and_keeps_the_backend_ones(mounted_frontend):
    for path in ("/view", "/view/some/deep/link"):
        resp = mounted_frontend.get(path)
        assert resp.status_code == 200 and "<title>bundle</title>" in resp.text, path
    assert mounted_frontend.get("/assets/index-abc.js").text == "console.log(1)"
    # Registered before the catch-all, so the bundle doesn't shadow them.
    assert "<title>bundle</title>" not in mounted_frontend.get("/login").text
    assert mounted_frontend.get("/worlds/test-world/summary").status_code == 401
    assert mounted_frontend.get("/").headers["location"].endswith("/view")


# --- Multiple worlds (#563) ------------------------------------------------------------


def _raw_client():
    from fastapi.testclient import TestClient

    return TestClient(app, base_url="http://localhost", client=("127.0.0.1", 50000))


def test_two_populated_worlds_do_not_share_state(client, saved_world):
    raw = _raw_client()
    a = raw.post("/worlds", json={"name": "A"}).json()["world_id"]
    b = raw.post("/worlds", json={"name": "B"}).json()["world_id"]
    assert a != b
    assert raw.post(f"/worlds/{a}/load", content=saved_world).status_code == 200
    assert raw.post(f"/worlds/{b}/load", content=saved_world).status_code == 200
    lines = raw.post(f"/worlds/{a}/step", json={"years": 1_000.0}).text.splitlines()
    assert json.loads(lines[-1])["type"] == "done"
    assert raw.get(f"/worlds/{a}/summary").json()["elapsed_years"] == 1_000.0
    assert raw.get(f"/worlds/{b}/summary").json()["elapsed_years"] == 0.0


def test_generate_debug_applies_dialog_name_and_public_flag(client):
    # The Generate World dialog's name + Public choice lands on the record being
    # generated into (#597).
    raw = _raw_client()
    wid = raw.post("/worlds", json={}).json()["world_id"]
    resp = raw.post(f"/worlds/{wid}/generate_debug", json={"scenario": "two_plate_divergent", "name": "Demo", "is_public": True})
    assert resp.status_code == 200
    (record,) = [w for w in raw.get("/worlds").json()["worlds"] if w["world_id"] == wid]
    assert record["name"] == "Demo"
    assert record["is_public"] is True
    listed = raw.get("/public/worlds").json()["worlds"]
    assert [w["world_id"] for w in listed] == [wid]
    assert raw.delete(f"/worlds/{wid}").status_code == 200


def test_registry_crud_public_listing_and_removed_old_routes(client, saved_world):
    raw = _raw_client()
    a = raw.post("/worlds", json={"name": "A"}).json()
    b = raw.post("/worlds", json={"name": "B"}).json()
    assert a["world_id"] != b["world_id"]
    names = [w["name"] for w in raw.get("/worlds").json()["worlds"]]
    assert {"A", "B"} <= set(names)
    assert raw.post(f"/worlds/{a['world_id']}/load", content=saved_world).status_code == 200
    raw.patch(f"/worlds/{a['world_id']}", json={"name": "Renamed", "is_public": True})
    # Only populated public worlds are listed or addressable.
    public = raw.get("/public/worlds").json()["worlds"]
    assert [w["name"] for w in public] == ["Renamed"]
    assert raw.get(f"/public/worlds/{b['world_id']}/world").status_code == 404
    assert raw.get("/world/summary").status_code == 404
    assert raw.post("/world/step", json={"years": 1}).status_code == 404
    assert raw.delete(f"/worlds/{a['world_id']}").status_code == 200
    assert raw.delete(f"/worlds/{b['world_id']}").status_code == 200
    assert {w["name"] for w in raw.get("/worlds").json()["worlds"]} == {"test"}


def test_legacy_checkpoint_imports_as_one_public_world(client, store, saved_world, tmp_path):
    from app import checkpoint

    root = checkpoint.CheckpointStore(tmp_path)
    (tmp_path / checkpoint.WORLDS_INDEX_FILE).unlink()
    saved = persistence.load_world_bytes(saved_world)
    saved.elapsed_years = 123_000.0
    root.save_world(saved)
    root.save_status(checkpoint.RunStatus(running=True, step_years=500.0))
    main._restore_checkpoint()
    restored = [s for s in main._registry.list() if s.world is not None]
    assert len(restored) == 1
    assert restored[0].is_public is True
    # The saved status was running, so the run resumes and steps on its own.
    assert restored[0].runner.status().running is True
    assert restored[0].world.elapsed_years >= 123_000.0
    assert (tmp_path / checkpoint.WORLDS_INDEX_FILE).exists()
    assert (tmp_path / checkpoint.WORLD_FILE).exists()  # legacy files left in place
    assert len(raw_list_worlds()) == 2  # the restored world plus this test's default session


def raw_list_worlds():
    return _raw_client().get("/worlds").json()["worlds"]


def test_delete_while_busy_is_refused(client):
    raw = _raw_client()
    world_id = raw.post("/worlds", json={"name": "Busy"}).json()["world_id"]
    session = main._registry.get(world_id)
    assert session.begin_operation()
    try:
        assert raw.delete(f"/worlds/{world_id}").status_code == 409
    finally:
        session.end_operation()
    assert raw.delete(f"/worlds/{world_id}").status_code == 200
    assert world_id not in {w["world_id"] for w in raw_list_worlds()}


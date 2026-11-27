"""Parity, thread-safety, cancellation, and fallback tests for Rust map rendering (#1039)."""

from __future__ import annotations

import base64
import concurrent.futures
import copy
import functools
import io
import os
import subprocess
import sys
import threading
import time
from types import SimpleNamespace

import av
import numpy as np
import pytest
from fastapi.testclient import TestClient
from PIL import Image
from scipy.ndimage import distance_transform_edt

from app import _render_rust, animation_job, checkpoint, geometry, main, persistence, plates, render_image
from app.hydroclimate import biomes, climate, hydrology, sea_ice
from app.world import World, generate_world, step_world


@functools.lru_cache(maxsize=8)
def _cached_stepped_world(seed: int = 11, steps: int = 2) -> World:
    w = generate_world(seed=seed, num_plates=6, node_density=0.5, climate_density=0.5)
    for _ in range(steps):
        step_world(w, 1_000_000.0)
    return w


def _stepped_world(seed: int = 11, steps: int = 2, *, mutate: bool = False) -> World:
    w = _cached_stepped_world(seed, steps)
    return copy.deepcopy(w) if mutate else w


def _decode_png(png_bytes: bytes) -> np.ndarray:
    return np.asarray(Image.open(io.BytesIO(png_bytes)))


def test_rust_render_extension_available():
    assert _render_rust.RUST_AVAILABLE is True
    assert _render_rust.render_pool_threads() >= 1


def test_biome_kernels_match_python(monkeypatch):
    rng = np.random.default_rng(7)
    h, w = 48, 96
    is_ocean = rng.random((h, w)) < 0.45
    temp = np.linspace(-25.0, 28.0, h)[:, None] + rng.normal(0.0, 3.0, (h, w))
    precip = np.exp(rng.uniform(3.0, 8.0, (h, w)))
    elev = np.where(is_ocean, -2000.0, rng.uniform(0.0, 4500.0, (h, w)))
    # Introduce exact ties to verify rankdata(method="average") parity.
    elev[10:20, 10:20] = 500.0
    lat = np.linspace(88.0, -88.0, h)
    has_ice = is_ocean & (temp < -2.0)
    sea_w = np.where(has_ice, 0.85, 0.15)

    cont_rs = biomes.grid_continentality(is_ocean)
    dist_rs = biomes.grid_dist_to_land_rad(is_ocean)
    w_id_rs, w_sh_rs, ru_rs = biomes.smooth_biome_field_blend(
        temp, precip, elev, np.zeros_like(elev), is_ocean,
        lat_deg=lat, axial_tilt_deg=23.5, has_sea_ice=has_ice, sea_ice_weight=sea_w,
    )
    shade_rs = biomes.biome_relative_shade_factor(w_id_rs, elev)

    monkeypatch.setattr(_render_rust, "RUST_AVAILABLE", False)
    cont_py = biomes.grid_continentality(is_ocean)
    dist_py = biomes.grid_dist_to_land_rad(is_ocean)
    w_id_py, w_sh_py, ru_py = biomes.smooth_biome_field_blend(
        temp, precip, elev, np.zeros_like(elev), is_ocean,
        lat_deg=lat, axial_tilt_deg=23.5, has_sea_ice=has_ice, sea_ice_weight=sea_w,
    )
    shade_py = biomes.biome_relative_shade_factor(w_id_py, elev)

    np.testing.assert_allclose(cont_rs, cont_py, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(dist_rs, dist_py, rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(w_id_rs, w_id_py)
    np.testing.assert_array_equal(ru_rs, ru_py)
    np.testing.assert_allclose(w_sh_rs, w_sh_py, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(shade_rs, shade_py, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("projection", ["behrmann", "eckert4", "globe"])
@pytest.mark.parametrize("view", ["combined", "biome"])
@pytest.mark.parametrize("rotated", [False, True])
def test_alpha_id_channel_bit_exact_parity(monkeypatch, projection, view, rotated):
    w = _stepped_world(seed=11, steps=2, mutate=True)
    # Inject a visible lake, an inland sea, and a glacier so all ID codes are exercised.
    if w.hydrology_cache is not None and len(w.hydrology_cache.points) > 20:
        land_nodes = np.flatnonzero(~w.hydrology_cache.is_ocean)
        if len(land_nodes) >= 6:
            plate = w.plates[0]
            n0 = min(plate.node_count(), 6)
            ld = plate._field_storage("lake_depth").copy()
            ld[:n0] = 50.0
            plate.set_fields_on_plate(lake_depth=ld)
            w.hydrology_cache.is_sea[land_nodes[0]] = True

    rot = (
        geometry.rotation_matrix(np.array([0.3, 0.7, 0.2]), 0.65)
        if rotated
        else np.eye(3)
    )
    rust_bytes = render_image.render_png(w, projection, view, 440, 244, rot)
    monkeypatch.setattr(_render_rust, "RUST_AVAILABLE", False)
    py_bytes = render_image.render_png(w, projection, view, 440, 244, rot)

    rust_rgba = _decode_png(rust_bytes)
    py_rgba = _decode_png(py_bytes)
    assert rust_rgba.shape == (244, 440, 4)
    assert py_rgba.shape == (244, 440, 4)
    np.testing.assert_array_equal(rust_rgba[:, :, 3], py_rgba[:, :, 3])


@pytest.mark.parametrize("projection", ["behrmann", "eckert4", "globe"])
@pytest.mark.parametrize("view", ["combined", "biome", "elevation"])
def test_rgb_visual_parity(monkeypatch, projection, view):
    w = _stepped_world(seed=17, steps=2)
    rot = geometry.rotation_matrix(np.array([0.0, 1.0, 0.2]), 0.4)
    kwargs = {}
    if view == "elevation":
        kwargs = {"show_mountains": True, "show_plains_plateaus": True}
    rust_img = _decode_png(render_image.render_png(w, projection, view, 440, 244, rot, **kwargs))
    monkeypatch.setattr(_render_rust, "RUST_AVAILABLE", False)
    py_img = _decode_png(render_image.render_png(w, projection, view, 440, 244, rot, **kwargs))

    rust_rgb = rust_img[:, :, :3].astype(np.int16)
    py_rgb = py_img[:, :, :3].astype(np.int16)
    diff = np.abs(rust_rgb - py_rgb)
    mean_abs = float(diff.mean())
    p99 = float(np.percentile(diff, 99.0))
    assert mean_abs <= 1.5, f"{projection}/{view} mean abs diff {mean_abs}"
    assert p99 <= 5.0, f"{projection}/{view} p99 diff {p99}"


def test_sea_ice_toggles_parity(monkeypatch):
    w = _stepped_world(seed=5, steps=2)
    for smooth, soft in [(True, False), (False, True), (False, False)]:
        rust_img = _decode_png(
            render_image.render_png(
                w, "behrmann", "combined", 320, 180, None,
                smooth_sea_ice=smooth, soft_sea_ice_edge=soft,
            )
        )
        monkeypatch.setattr(_render_rust, "RUST_AVAILABLE", False)
        py_img = _decode_png(
            render_image.render_png(
                w, "behrmann", "combined", 320, 180, None,
                smooth_sea_ice=smooth, soft_sea_ice_edge=soft,
            )
        )
        monkeypatch.setattr(_render_rust, "RUST_AVAILABLE", True)
        np.testing.assert_array_equal(rust_img[:, :, 3], py_img[:, :, 3])
        assert float(np.abs(rust_img[:, :, :3].astype(int) - py_img[:, :, :3].astype(int)).mean()) <= 1.5


def test_disable_rust_env_falls_back_cleanly():
    env = os.environ.copy()
    env["MANTLE_BLOOM_DISABLE_RUST"] = "1"
    code = (
        "from app import _render_rust, render_image; "
        "from app.world import generate_world; "
        "assert _render_rust.RUST_AVAILABLE is False; "
        "w = generate_world(seed=1, num_plates=4, node_density=0.5); "
        "b = render_image.render_png(w, 'behrmann', 'combined', 120, 68); "
        "assert b[:8] == b'\\x89PNG\\r\\n\\x1a\\n'"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        cwd=os.path.dirname(os.path.dirname(__file__)),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr


def test_concurrent_rust_render_bypasses_compute_gate_and_matches_sequential():
    w = _stepped_world(seed=3, steps=1)
    climate.compute_climate_cached(w)
    snap = _render_rust.build_render_snapshot(w, "behrmann", "combined", 240, 136, None)
    seq_png, _ = _render_rust.render_snapshot(snap, want_png=True, want_rgb=False)
    assert seq_png is not None

    # Hold _compute_gate on another thread and run step_world concurrently on a clone;
    # Rust snapshot rendering must finish immediately without waiting on _compute_gate.
    w_step = copy.deepcopy(w)
    gate_acquired = threading.Event()
    release_gate = threading.Event()
    step_done = threading.Event()
    step_err: list[BaseException] = []

    def _hold_and_step():
        with main._compute_gate.hold():
            gate_acquired.set()
            try:
                step_world(w_step, 500_000.0)
                step_done.set()
            except BaseException as exc:  # noqa: BLE001 -- re-assert cleanly on caller thread
                step_err.append(exc)
            release_gate.wait(timeout=60.0)

    t = threading.Thread(target=_hold_and_step, daemon=True)
    t.start()
    assert gate_acquired.wait(timeout=15.0)
    try:
        conc_png, _ = _render_rust.render_snapshot(snap, want_png=True, want_rgb=False)
        assert conc_png == seq_png
    finally:
        release_gate.set()
        t.join(timeout=60.0)
    assert step_err == []
    assert step_done.is_set()


def test_http_render_endpoint_bypasses_compute_gate_for_rust_views():
    w = _stepped_world(seed=3, steps=1)
    climate.compute_climate_cached(w)
    session = main._registry.create(name="rust-gate-test")
    main._wire_session(session)
    session.world = w
    gate_acquired = threading.Event()
    release_gate = threading.Event()

    def _hold_gate():
        with main._compute_gate.hold():
            gate_acquired.set()
            release_gate.wait(timeout=60.0)

    t = threading.Thread(target=_hold_gate, daemon=True)
    t.start()
    assert gate_acquired.wait(timeout=15.0)
    try:
        with TestClient(main.app, base_url="http://localhost", client=("127.0.0.1", 50000)) as client:
            resp = client.get(
                f"/worlds/{session.world_id}/render",
                params={"projection": "behrmann", "view": "combined", "width": 120, "height": 68},
            )
        assert resp.status_code == 200
        decoded = base64.b64decode(resp.json()["image_base64"])
        assert decoded[:8] == b"\x89PNG\r\n\x1a\n"
    finally:
        release_gate.set()
        t.join(timeout=60.0)
        main._registry.remove(session.world_id)


@pytest.mark.parametrize("view", ["combined", "elevation"])
def test_parallel_animation_matches_sequential_deterministically(view):
    w_seq = _stepped_world(seed=9, steps=1, mutate=True)
    w_par = copy.deepcopy(w_seq)

    seq_msgs = list(
        render_image.stream_animation_mp4(
            w_seq, "behrmann", view, 180, 100, None, 500_000.0, 1, 4,
            parallel_render=False,
        )
    )
    par_msgs = list(
        render_image.stream_animation_mp4(
            w_par, "behrmann", view, 180, 100, None, 500_000.0, 1, 4,
            parallel_render=True,
        )
    )

    seq_prog = [m for m in seq_msgs if m[0] == "progress"]
    par_prog = [m for m in par_msgs if m[0] == "progress"]
    assert [m[1] for m in seq_prog] == [1, 2, 3, 4]
    assert [m[1] for m in par_prog] == [1, 2, 3, 4]
    assert [m[4] for m in par_prog] == [m[4] for m in seq_prog]
    for ms, mp in zip(seq_prog, par_prog):
        assert mp[3] == ms[3]

    assert w_par.steps_taken == w_seq.steps_taken
    assert w_par.elapsed_years == w_seq.elapsed_years
    for ps, pp in zip(w_seq.plates, w_par.plates):
        np.testing.assert_array_equal(pp._field_storage("elevation"), ps._field_storage("elevation"))

    with av.open(io.BytesIO(par_msgs[-1][1])) as c_par, av.open(io.BytesIO(seq_msgs[-1][1])) as c_seq:
        assert sum(1 for _ in c_par.decode(video=0)) == sum(1 for _ in c_seq.decode(video=0)) == 4


def test_parallel_animation_stop_event_and_previews_leak_zero_threads():
    w = _stepped_world(seed=13, steps=1, mutate=True)
    stop_event = threading.Event()
    previews: list[tuple[bytes, float]] = []
    calls = 0

    def _step(world, years):
        nonlocal calls
        calls += 1
        time.sleep(0.03)
        step_world(world, years)
        if calls == 4:
            stop_event.set()

    def _on_preview(png_bytes, elapsed_years, _stats, _rate):
        previews.append((png_bytes, elapsed_years))

    msgs = list(
        render_image.stream_animation_mp4(
            w, "globe", "combined", 160, 90, None, 250_000.0, 3, 5,
            step_fn=_step,
            stop_event=stop_event,
            on_preview=_on_preview,
            parallel_render=True,
        )
    )
    prog = [m for m in msgs if m[0] == "progress"]
    assert [m[1] for m in prog] == [1, 2]
    assert msgs[-1][0] == "done"
    assert msgs[-1][2] is True
    assert len(previews) >= 1
    assert previews[0][0][:8] == b"\x89PNG\r\n\x1a\n"

    leaked = [
        t.name
        for t in threading.enumerate()
        if t.name.startswith("animate-render") or t.name.startswith("animate-preview")
    ]
    assert leaked == []


def test_parallel_animation_mid_batch_failure_keeps_mp4_frame_without_progress_yield():
    w = _stepped_world(seed=15, steps=1, mutate=True)
    step_calls = 0

    def _step_then_fail(world, years):
        nonlocal step_calls
        step_calls += 1
        if step_calls == 3:
            raise RuntimeError("boom")
        step_world(world, years)

    msgs = list(
        render_image.stream_animation_mp4(
            w, "behrmann", "combined", 160, 90, None, 250_000.0, 2, 4,
            step_fn=_step_then_fail,
            parallel_render=True,
        )
    )
    progress_frames = [m[1] for m in msgs if m[0] == "progress"]
    assert progress_frames == [1]

    done = msgs[-1]
    assert done[0] == "done"
    assert done[2] is False
    assert "RuntimeError: boom" in done[3]

    with av.open(io.BytesIO(done[1])) as container:
        assert sum(1 for _ in container.decode(video=0)) == 2


def test_all_land_and_all_ocean_edt_continentality_biome_parity(monkeypatch):
    land = np.ones((4, 8), dtype=bool)
    np.testing.assert_allclose(_render_rust.grid_edt(land), distance_transform_edt(land), rtol=0.0, atol=1e-12)

    all_land_ocean = np.zeros((4, 8), dtype=bool)
    all_ocean = np.ones((4, 8), dtype=bool)
    temp = np.full((4, 8), 10.0, dtype=np.float64)
    precip = np.full((4, 8), 600.0, dtype=np.float64)
    elev = np.full((4, 8), 200.0, dtype=np.float64)
    slope = np.zeros((4, 8), dtype=np.float64)
    lat_1d = np.full(4, 45.0, dtype=np.float64)

    rs_cont = biomes.grid_continentality(all_land_ocean)
    rs_dist = biomes.grid_dist_to_land_rad(all_land_ocean)
    rs_ids = biomes.smooth_biome_field(temp, precip, elev, slope, all_land_ocean, lat_deg=lat_1d)
    rs_blend = _render_rust.smooth_biome_field_blend(temp, precip, all_land_ocean, lat_1d, 23.5)

    monkeypatch.setattr(_render_rust, "RUST_AVAILABLE", False)
    py_cont = biomes.grid_continentality(all_land_ocean)
    py_dist = biomes.grid_dist_to_land_rad(all_land_ocean)
    py_ids = biomes.smooth_biome_field(temp, precip, elev, slope, all_land_ocean, lat_deg=lat_1d)

    np.testing.assert_allclose(rs_cont, py_cont, atol=1e-12)
    np.testing.assert_allclose(rs_dist, py_dist, atol=1e-12)
    np.testing.assert_array_equal(rs_ids, py_ids)
    np.testing.assert_array_equal(rs_blend[0], py_ids)
    assert np.all(rs_ids == biomes.koppen_index("Dfa"))

    # All-ocean EDT must stay 0.0 for distance-to-land and continentality.
    monkeypatch.undo()
    assert np.all(biomes.grid_continentality(all_ocean) == 0.0)
    assert np.all(biomes.grid_dist_to_land_rad(all_ocean) == 0.0)


def test_smooth_biome_field_blend_sea_ice_shape_mismatch_raises_value_error():
    temp = np.zeros((4, 8), dtype=np.float64)
    precip = np.full((4, 8), 500.0, dtype=np.float64)
    is_ocean = np.ones((4, 8), dtype=bool)
    lat = np.linspace(60.0, -60.0, 4)

    with pytest.raises(ValueError, match="has_sea_ice"):
        _render_rust.smooth_biome_field_blend(
            temp, precip, is_ocean, lat, 23.5, has_sea_ice=np.zeros((3, 8), dtype=bool)
        )
    with pytest.raises(ValueError, match="sea_ice_weight"):
        _render_rust.smooth_biome_field_blend(
            temp, precip, is_ocean, lat, 23.5, sea_ice_weight=np.zeros((3, 8), dtype=np.float64)
        )


def test_out_of_range_biome_ids_shade_factor_parity(monkeypatch):
    ids = np.array([-5, -5, 0, 0, 40, 42, 42, 99, 99, 99], dtype=np.int64)
    elev = np.array([10.0, 50.0, 100.0, 100.0, 20.0, 5.0, 15.0, 10.0, 10.0, 40.0], dtype=np.float64)
    rs_shade = biomes.biome_relative_shade_factor(ids, elev)
    monkeypatch.setattr(_render_rust, "RUST_AVAILABLE", False)
    py_shade = biomes.biome_relative_shade_factor(ids, elev)
    np.testing.assert_allclose(rs_shade, py_shade, rtol=0.0, atol=1e-12)


def test_empty_zero_node_plates_fall_back_gracefully(monkeypatch):
    w = generate_world(seed=5, num_plates=3, node_density=0.5, climate_density=0.5)
    for plate in w.plates:
        monkeypatch.setattr(
            plate,
            "all_points_and_elevation",
            lambda: (np.empty((0, 3), dtype=np.float64), np.empty(0, dtype=np.float64)),
        )
        monkeypatch.setattr(plate, "get_bounding_polygon", lambda: np.empty((0, 3), dtype=np.float64))
    assert plates.collect_all_points(w.plates) is None
    assert _render_rust.build_render_snapshot(w, "behrmann", "elevation", 120, 68, None) is None
    png_bytes = render_image.render_png(w, "behrmann", "elevation", 120, 68)
    assert png_bytes[:8] == b"\x89PNG\r\n\x1a\n"


def test_spill_path_rs_matches_python_numba_paths():
    # Diamond graph with an equal-elevation tie (source 0 -> {1, 2} -> 3), high rim 4, ocean 5.
    adj = [
        [1, 2, 4],
        [0, 3],
        [0, 3],
        [1, 2, 5],
        [0, 5],
        [3, 4],
    ]
    indptr = np.zeros(len(adj) + 1, dtype=np.int64)
    flat = []
    for i, nbrs in enumerate(adj):
        indptr[i + 1] = indptr[i] + len(nbrs)
        flat.extend(nbrs)
    indices = np.asarray(flat, dtype=np.int64)
    elev = np.array([10.0, 20.0, 20.0, 15.0, 80.0, -5.0], dtype=np.float64)
    is_ocean = np.array([False, False, False, False, False, True], dtype=np.bool_)

    for target, ceiling, max_nodes in [(5, 50.0, 64), (5, 19.0, 64), (5, 50.0, 2), (3, 50.0, 64)]:
        py_pass = hydrology._spill_pass_level(indptr, indices, elev, is_ocean, 0, target, float(ceiling), int(max_nodes))
        py_path = [] if not np.isfinite(py_pass) else list(
            hydrology._fewest_hops_path(indptr, indices, elev, is_ocean, 0, target, float(py_pass))
        )
        rs_path = _render_rust.spill_path(indptr, indices, elev, is_ocean, 0, target, ceiling, max_nodes)
        assert rs_path == py_path


def test_can_use_rust_biomes_monkeypatch_guards(monkeypatch):
    elev = np.zeros((4, 8), dtype=np.float64)
    lat = np.linspace(60.0, -60.0, 4)
    assert biomes._can_use_rust_biomes(elev, lat)[0] is True

    for attr in ("classify_koppen", "classify_koppen_soft", "classify_pelagic", "classify_pelagic_soft", "grid_continentality", "grid_dist_to_land_rad"):
        orig = getattr(biomes, attr)
        monkeypatch.setattr(biomes, attr, lambda *a, _o=orig, **kw: _o(*a, **kw))
        assert biomes._can_use_rust_biomes(elev, lat) == (False, None)
        monkeypatch.undo()


def test_pipelined_animation_boundary_checkpoint_controls_and_failure_match_sequential(tmp_path, monkeypatch):
    base_world = _stepped_world(seed=19, steps=1, mutate=True)
    world_bytes = persistence.save_world_bytes(base_world)

    def _run_mode(mode_name: str, *, disable_parallel: bool) -> tuple[float, int | None, list[tuple[int, float]]]:
        session = main._registry.create(name=f"test-{mode_name}")
        store = checkpoint.CheckpointStore(tmp_path / mode_name).for_world(session.world_id)
        session.checkpoint_store = store
        main._wire_session(session)
        session.world = persistence.load_world_bytes(world_bytes)
        job = animation_job.AnimationJob(
            total=5,
            unbounded=False,
            projection="behrmann",
            view="elevation",
            elapsed_years=session.world.elapsed_years,
        )
        assert session.begin_animation_operation(job) is None
        session.set_controls_snapshot(main._controls_state(session.world))

        step_solars: list[tuple[int, float]] = []
        calls = 0

        def _fake_step(_sess, w, years):
            nonlocal calls
            calls += 1
            step_solars.append((calls, float(w.solar_multiplier)))
            if calls == 6:
                raise RuntimeError("boom at frame 4")
            step_world(w, years)
            if calls == 2:
                session.try_queue_pending_controls({"solar_multiplier": 1.25})

        monkeypatch.setattr(main, "_step_world", _fake_step)
        if disable_parallel:
            monkeypatch.setenv("MANTLE_BLOOM_PARALLEL_RENDER", "0")
        else:
            monkeypatch.delenv("MANTLE_BLOOM_PARALLEL_RENDER", raising=False)

        req = main.AnimateRequest(
            projection="behrmann", view="elevation", width=120, height=68,
            step_years=250_000.0, steps_per_frame=2, num_frames=5,
        )
        try:
            main._run_animation_worker(session, req, None, job)
        finally:
            main._registry.remove(session.world_id)

        ckpt = persistence.load_world_bytes(store.world_path.read_bytes())
        rec = store.load_animation_progress()
        return float(ckpt.elapsed_years), None if rec is None else int(rec["completed_frames"]), step_solars

    try:
        seq_age, seq_completed, seq_solars = _run_mode("seq", disable_parallel=True)
        par_age, par_completed, par_solars = _run_mode("par", disable_parallel=False)
    finally:
        monkeypatch.delenv("MANTLE_BLOOM_PARALLEL_RENDER", raising=False)

    expected_frame3_age = base_world.elapsed_years + 4 * 250_000.0
    assert seq_completed == par_completed == 3
    assert seq_age == par_age == expected_frame3_age
    assert seq_solars == par_solars
    assert [s for _, s in par_solars] == [1.0, 1.0, 1.25, 1.25, 1.25, 1.25]


def test_parallel_animation_keeps_stream_ages_monotonic_and_never_clobbers_preview(monkeypatch):
    w = _stepped_world(seed=21, steps=1, mutate=True)
    now = [0.0]

    def _fast_clock():
        now[0] += render_image.ANIMATION_HEARTBEAT_INTERVAL_S + 1.0
        return now[0]

    monkeypatch.setattr("app.world._monotonic", _fast_clock)
    events: list[tuple[str, float]] = []

    def _on_preview(_png, elapsed_years, _stats, _rate):
        events.append(("preview", float(elapsed_years)))

    for msg in render_image.stream_animation_mp4(
        w, "behrmann", "combined", 120, 68, None, 250_000.0, 3, 5,
        on_preview=_on_preview,
        parallel_render=True,
    ):
        if msg[0] == "progress":
            events.append(("progress", float(msg[4])))
        elif msg[0] == "heartbeat":
            events.append(("heartbeat", float(msg[3])))

    ages = [age for _, age in events]
    assert ages == sorted(ages)
    assert any(kind == "preview" for kind, _ in events)
    assert any(kind == "heartbeat" for kind, _ in events)


def test_http_render_admission_bounded_by_semaphore(monkeypatch):
    w = _stepped_world(seed=23, steps=2, mutate=False)
    session = main._registry.create(name="render-limit")
    main._wire_session(session)
    session.world = w

    lock = threading.Lock()
    active = 0
    max_active = 0
    orig_render_snapshot = _render_rust.render_snapshot

    def _counted_render(*args, **kwargs):
        nonlocal active, max_active
        with lock:
            active += 1
            max_active = max(max_active, active)
        try:
            time.sleep(0.05)
            return orig_render_snapshot(*args, **kwargs)
        finally:
            with lock:
                active -= 1

    orig_sem = main._rust_render_semaphore
    main._rust_render_semaphore = threading.Semaphore(1)
    monkeypatch.setattr(_render_rust, "render_snapshot", _counted_render)
    try:
        with TestClient(main.app, base_url="http://localhost", client=("127.0.0.1", 50000)) as client:
            def _one():
                return client.get(
                    f"/worlds/{session.world_id}/render",
                    params={"projection": "behrmann", "view": "combined", "width": 120, "height": 68},
                ).status_code

            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
                statuses = [f.result() for f in [pool.submit(_one) for _ in range(4)]]
        assert statuses == [200, 200, 200, 200]
        assert max_active == 1
    finally:
        main._rust_render_semaphore = orig_sem
        main._registry.remove(session.world_id)


def test_spill_path_and_render_entrypoints_validate_inputs_without_panicking():
    rust = _render_rust._rust
    assert rust is not None

    indptr = np.array([0, 1, 2], dtype=np.int64)
    indices = np.array([1, 0], dtype=np.int64)
    elev = np.array([10.0, 5.0], dtype=np.float64)
    ocean = np.array([False, True], dtype=np.bool_)

    with pytest.raises(ValueError, match="out of bounds"):
        rust.spill_path_rs(indptr, indices, elev, ocean, 99, 1, 50.0, 16)
    with pytest.raises(ValueError, match="out of bounds"):
        rust.spill_path_rs(indptr, indices, elev, ocean, 0, 99, 50.0, 16)
    with pytest.raises(ValueError, match="shape mismatch"):
        rust.spill_path_rs(np.array([0, 1], dtype=np.int64), indices, elev, ocean, 0, 1, 50.0, 16)
    with pytest.raises(ValueError, match="neighbor index out of bounds"):
        rust.spill_path_rs(indptr, np.array([5, 0], dtype=np.int64), elev, ocean, 0, 1, 50.0, 16)

    w = _stepped_world(seed=7, steps=1)
    climate.compute_climate_cached(w)
    snap = _render_rust.build_render_snapshot(w, "behrmann", "combined", 120, 68, None)
    assert snap is not None

    with pytest.raises(ValueError, match="point array length mismatch"):
        rust.render_biome_or_combined(
            "combined", snap.projection, snap.width, snap.height, snap.view_rotation,
            snap.grid_h, snap.grid_w, snap.sea_level_m, snap.axial_tilt_deg, snap.solar_multiplier,
            snap.smooth_sea_ice, snap.soft_sea_ice_edge,
            snap.points, snap.elevation[:-1], snap.lake_depth, snap.lake_ice, snap.glacier_depth,
            snap.hillshade, snap.channel_depth, snap.channel_width,
            snap.hydro_points_same, snap.hydro_points, snap.hydro_is_ocean, snap.hydro_is_sea,
            snap.hydro_perennial_lake_ice,
            snap.clim_air_temp, snap.clim_ocean_temp, snap.clim_precip, snap.clim_lat_deg,
            snap.clim_elev_m, snap.clim_is_ocean, snap.clim_h, snap.clim_w,
            snap.river_from_xyz, snap.river_to_xyz, snap.river_width_frac,
            snap.river_color_alpha, snap.river_canyon_fill, True, False, snap.stroke_scale_override,
        )

    snap_elev = _render_rust.build_render_snapshot(w, "behrmann", "elevation", 120, 68, None)
    assert snap_elev is not None
    with pytest.raises(ValueError, match="render dimensions"):
        rust.render_elevation(
            snap_elev.projection, 0, 68, snap_elev.view_rotation,
            snap_elev.grid_spacing_rad, snap_elev.sea_level_m, snap_elev.axial_tilt_deg,
            snap_elev.solar_multiplier, snap_elev.smooth_sea_ice, snap_elev.soft_sea_ice_edge,
            snap_elev.show_mountains, snap_elev.show_plains_plateaus, snap_elev.tectonics_bbox,
            snap_elev.points, snap_elev.elevation, snap_elev.lake_depth, snap_elev.lake_ice,
            snap_elev.glacier_depth, snap_elev.hillshade, snap_elev.channel_depth,
            snap_elev.channel_width, snap_elev.terrain_relief,
            snap_elev.hydro_points_same, snap_elev.hydro_points, snap_elev.hydro_is_ocean,
            snap_elev.hydro_is_sea, snap_elev.hydro_perennial_lake_ice,
            snap_elev.clim_air_temp, snap_elev.clim_ocean_temp, snap_elev.clim_lat_deg,
            snap_elev.clim_elev_m, snap_elev.clim_is_ocean, snap_elev.clim_h, snap_elev.clim_w,
            snap_elev.river_from_xyz, snap_elev.river_to_xyz, snap_elev.river_width_frac[:-1],
            snap_elev.river_color_alpha, snap_elev.river_canyon_fill, True, False,
            snap_elev.stroke_scale_override,
        )


def test_biome_relative_shade_factor_propagates_nan_like_python(monkeypatch):
    # Bucket 5 has a single NaN (len == 1 -> 1.0); bucket 1 has a NaN in len == 3 -> all NaN;
    # bucket 2 is finite (len == 2 -> 0.75, 1.25).
    ids = np.array([1, 1, 1, 2, 2, 5], dtype=np.int64)
    elev = np.array([10.0, np.nan, 20.0, 5.0, 15.0, np.nan], dtype=np.float64)
    rs_out = biomes.biome_relative_shade_factor(ids, elev)
    monkeypatch.setattr(_render_rust, "RUST_AVAILABLE", False)
    py_out = biomes.biome_relative_shade_factor(ids, elev)
    np.testing.assert_allclose(rs_out, py_out, rtol=0.0, atol=1e-12, equal_nan=True)


def test_rust_render_constants_match_python_tables():
    rust = _render_rust._rust
    assert rust is not None
    c = rust.render_constants_rs()

    assert np.array_equal(np.asarray([list(row) for row in c["biome_colors"]], dtype=np.uint8), biomes.BIOME_COLORS)
    np.testing.assert_allclose(c["elevation_stop_e"], render_image._ELEVATION_STOP_E)
    np.testing.assert_allclose(c["elevation_stop_rgb"], render_image._ELEVATION_STOP_RGB)
    assert tuple(c["background_rgb"]) == render_image.BACKGROUND_RGB
    assert tuple(c["lake_color_rgb"]) == render_image.LAKE_COLOR_RGB
    assert tuple(c["sea_color_rgb"]) == render_image.SEA_COLOR_RGB
    assert tuple(c["glacier_color_rgb"]) == render_image.GLACIER_COLOR_RGB
    assert tuple(c["river_color_rgb"]) == render_image.RIVER_COLOR_RGB
    assert tuple(c["mountain_overlay_rgb"]) == render_image.MOUNTAIN_OVERLAY_RGB
    assert tuple(c["plains_plateau_overlay_rgb"]) == render_image.PLAINS_PLATEAU_OVERLAY_RGB
    assert c["terrain_overlay_alpha"] == pytest.approx(render_image.TERRAIN_OVERLAY_ALPHA)
    assert c["cell_overlap_factor"] == pytest.approx(render_image.CELL_OVERLAP_FACTOR)
    assert c["cell_blur_radius_px"] == pytest.approx(render_image.CELL_BLUR_RADIUS_PX)
    assert c["biome_shade_amplitude"] == pytest.approx(biomes.BIOME_SHADE_AMPLITUDE)
    assert c["continentality_saturation_rad"] == pytest.approx(biomes.CONTINENTALITY_SATURATION_RAD)
    assert c["shelf_range_rad"] == pytest.approx(biomes.SHELF_RANGE_RAD)
    assert c["seawater_freezing_c"] == pytest.approx(sea_ice.SEAWATER_FREEZING_C)
    assert c["ocean_heat_flux_w_m2"] == pytest.approx(sea_ice.OCEAN_HEAT_FLUX_W_M2)
    assert c["polar_sea_ice_fraction"] == pytest.approx(sea_ice.POLAR_SEA_ICE_FRACTION)
    assert c["perennial_cover_low"] == pytest.approx(sea_ice.PERENNIAL_COVER_LOW)
    assert c["perennial_cover_high"] == pytest.approx(sea_ice.PERENNIAL_COVER_HIGH)
    assert c["combined_lake_id_code"] == render_image.COMBINED_LAKE_ID_CODE
    assert c["combined_sea_id_code"] == render_image.COMBINED_SEA_ID_CODE


def test_stale_rust_module_missing_render_exports_falls_back_cleanly(caplog):
    assert _render_rust._resolve_rust_module(_render_rust._rust) is _render_rust._rust
    with caplog.at_level("WARNING", logger=_render_rust.__name__):
        resolved = _render_rust._resolve_rust_module(SimpleNamespace(step_ice=object()))
    assert resolved is None
    assert "bin/build_rust.py" in caplog.text





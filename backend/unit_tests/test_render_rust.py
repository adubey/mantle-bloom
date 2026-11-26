"""Parity, thread-safety, cancellation, and fallback tests for Rust map rendering (#1039)."""

from __future__ import annotations

import base64
import copy
import functools
import io
import os
import subprocess
import sys
import threading
import time

import av
import numpy as np
import pytest
from fastapi.testclient import TestClient
from PIL import Image

from app import _render_rust, geometry, main, render_image
from app.hydroclimate import biomes, climate
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

    def _hold_and_step():
        with main._compute_gate.hold():
            gate_acquired.set()
            step_world(w_step, 500_000.0)
            step_done.set()
            release_gate.wait(timeout=5.0)

    t = threading.Thread(target=_hold_and_step, daemon=True)
    t.start()
    assert gate_acquired.wait(timeout=2.0)
    try:
        conc_png, _ = _render_rust.render_snapshot(snap, want_png=True, want_rgb=False)
        assert conc_png == seq_png
    finally:
        release_gate.set()
        t.join(timeout=5.0)
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
            release_gate.wait(timeout=5.0)

    t = threading.Thread(target=_hold_gate, daemon=True)
    t.start()
    assert gate_acquired.wait(timeout=2.0)
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
        t.join(timeout=5.0)
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


"""Parity and fallback tests for `mantle_bloom_rust` lake kernels (#1046)."""

from __future__ import annotations

import os
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial import cKDTree

from app.hydroclimate import _lakes_rust, hydrology, lakes


def _random_knn_graph(rng: np.random.Generator, n: int, k: int, *, asymmetric: bool = False):
    pts = rng.normal(size=(n, 3))
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    nb = cKDTree(pts).query(pts, k + 1)[1][:, 1:].astype(np.int64)
    if asymmetric and n > k + 2:
        # Replace half the last-column entries with random nodes to force one-directional edges.
        mask = rng.random(n) < 0.5
        nb[mask, -1] = rng.integers(0, n, size=int(mask.sum()), dtype=np.int64)
    return nb


def _assert_lake_equal(rl: lakes.Lake, pl: lakes.Lake) -> None:
    assert rl.lake_id == pl.lake_id
    np.testing.assert_array_equal(rl.members, pl.members)
    assert rl.floor_elevation == pl.floor_elevation
    assert rl.min_depth == pl.min_depth
    assert rl.max_depth == pl.max_depth
    assert rl.outlet_node_idx == pl.outlet_node_idx
    assert rl.outlet_target_idx == pl.outlet_target_idx
    assert rl.sink_node_idx == pl.sink_node_idx
    assert rl.drains_at_floor == pl.drains_at_floor
    assert rl.current_water_elevation == pl.current_water_elevation
    assert len(rl.children) == len(pl.children)
    for rc, pc in zip(rl.children, pl.children):
        _assert_lake_equal(rc, pc)


def _assert_forest_equal(rust_roots: list[lakes.Lake], py_roots: list[lakes.Lake]) -> None:
    assert len(rust_roots) == len(py_roots)
    for rl, pl in zip(rust_roots, py_roots):
        _assert_lake_equal(rl, pl)


@pytest.mark.skipif(not _lakes_rust.RUST_AVAILABLE, reason="mantle_bloom_rust lake exports unavailable")
@pytest.mark.parametrize("n,k,asymmetric", [(600, 6, False), (1200, 8, True)])
def test_lake_components_rust_matches_python(monkeypatch, n: int, k: int, asymmetric: bool):
    rng = np.random.default_rng(100 + n + int(asymmetric))
    nb = _random_knn_graph(rng, n, k, asymmetric=asymmetric)
    for frac in (0.05, 0.35, 0.9, 1.0):
        is_lake = rng.random(n) < frac
        rust_groups = hydrology.lake_components(is_lake, nb)
        monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", False)
        py_groups = hydrology.lake_components(is_lake, nb)
        monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", True)

        assert len(rust_groups) == len(py_groups)
        for rg, pg in zip(rust_groups, py_groups):
            assert rg.dtype == np.int64
            np.testing.assert_array_equal(rg, pg)


@pytest.mark.skipif(not _lakes_rust.RUST_AVAILABLE, reason="mantle_bloom_rust lake exports unavailable")
def test_lake_components_rust_edge_cases(monkeypatch):
    cases = [
        (np.zeros(0, dtype=bool), np.zeros((0, 6), dtype=np.int64)),
        (np.zeros(20, dtype=bool), np.zeros((20, 6), dtype=np.int64)),
        (np.array([True, False, True, True], dtype=bool), np.zeros((4, 0), dtype=np.int64)),
        # Self-loops / isolated singletons
        (np.array([True, False, True, True], dtype=bool), np.array([[0], [1], [2], [3]], dtype=np.int64)),
    ]
    for is_lake, nb in cases:
        rust_groups = hydrology.lake_components(is_lake, nb)
        monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", False)
        py_groups = hydrology.lake_components(is_lake, nb)
        monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", True)
        assert len(rust_groups) == len(py_groups)
        for rg, pg in zip(rust_groups, py_groups):
            np.testing.assert_array_equal(rg, pg)


@pytest.mark.skipif(not _lakes_rust.RUST_AVAILABLE, reason="mantle_bloom_rust lake exports unavailable")
@pytest.mark.parametrize("integer_ties", [False, True])
@pytest.mark.parametrize("ocean_mode", ["mixed", "none", "all"])
def test_catchment_roots_rust_matches_python(monkeypatch, integer_ties: bool, ocean_mode: str):
    assert lakes._OCEAN_CATCHMENT == -2
    rng = np.random.default_rng(200 + int(integer_ties) * 10 + len(ocean_mode))
    n, k = 1000, 8
    nb = _random_knn_graph(rng, n, k, asymmetric=True)
    elev = rng.uniform(-500.0, 2000.0, size=n)
    if integer_ties:
        elev = np.round(elev / 200.0) * 200.0
    if ocean_mode == "mixed":
        is_ocean = elev < 0.0
    elif ocean_mode == "none":
        is_ocean = np.zeros(n, dtype=bool)
    else:
        is_ocean = np.ones(n, dtype=bool)

    rust_roots = lakes._catchment_roots(elev, is_ocean, nb)
    monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", False)
    py_roots = lakes._catchment_roots(elev, is_ocean, nb)
    monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", True)

    assert rust_roots == py_roots


@pytest.mark.skipif(not _lakes_rust.RUST_AVAILABLE, reason="mantle_bloom_rust lake exports unavailable")
@pytest.mark.parametrize("integer_ties", [False, True])
@pytest.mark.parametrize("use_interface_pass", [False, True])
@pytest.mark.parametrize("ocean_mode", ["mixed", "none", "all"])
def test_build_lake_hierarchy_rust_matches_python(
    monkeypatch,
    integer_ties: bool,
    use_interface_pass: bool,
    ocean_mode: str,
):
    rng = np.random.default_rng(300 + int(integer_ties) * 100 + int(use_interface_pass) * 10 + len(ocean_mode))
    n, k = 1200, 8
    nb = _random_knn_graph(rng, n, k, asymmetric=True)
    elev = rng.uniform(-400.0, 2500.0, size=n)
    if integer_ties:
        elev = np.round(elev / 150.0) * 150.0
    if ocean_mode == "mixed":
        is_ocean = elev < 0.0
    elif ocean_mode == "none":
        is_ocean = np.zeros(n, dtype=bool)
    else:
        is_ocean = np.ones(n, dtype=bool)

    passes = None
    if use_interface_pass:
        passes = np.maximum(elev[:, None], elev[nb]) - rng.uniform(0.0, 20.0, size=nb.shape)
        if integer_ties:
            passes = np.round(passes / 50.0) * 50.0

    rust_forest = lakes.build_lake_hierarchy(elev, is_ocean, nb, interface_pass_elevation=passes)
    monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", False)
    py_forest = lakes.build_lake_hierarchy(elev, is_ocean, nb, interface_pass_elevation=passes)
    monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", True)

    _assert_forest_equal(rust_forest, py_forest)


@pytest.mark.skipif(not _lakes_rust.RUST_AVAILABLE, reason="mantle_bloom_rust lake exports unavailable")
def test_build_lake_hierarchy_rust_edge_cases(monkeypatch):
    edge_cases = [
        (np.zeros(0), np.zeros(0, dtype=bool), np.zeros((0, 6), dtype=np.int64)),
        (np.array([10.0, 20.0]), np.zeros(2, dtype=bool), np.zeros((2, 0), dtype=np.int64)),
        # Self-loops / isolated single node
        (np.array([5.0]), np.array([False]), np.array([[0]], dtype=np.int64)),
        # Disconnected components with isolated self-loop nodes
        (
            np.array([10.0, 20.0, 5.0, -10.0]),
            np.array([False, False, False, True]),
            np.array([[0, 0], [1, 1], [3, 3], [2, 2]], dtype=np.int64),
        ),
    ]
    for elev, is_ocean, nb in edge_cases:
        rust_forest = lakes.build_lake_hierarchy(elev, is_ocean, nb)
        monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", False)
        py_forest = lakes.build_lake_hierarchy(elev, is_ocean, nb)
        monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", True)
        _assert_forest_equal(rust_forest, py_forest)


def test_disable_rust_env_falls_back_cleanly():
    env = os.environ.copy()
    env["MANTLE_BLOOM_DISABLE_RUST"] = "1"
    code = (
        "import numpy as np; "
        "from app.hydroclimate import _lakes_rust, hydrology, lakes; "
        "assert _lakes_rust.RUST_AVAILABLE is False; "
        "elev = np.array([30.0, 50.0, 20.0, -10.0]); "
        "is_ocean = np.array([False, False, False, True]); "
        "nb = np.array([[1, 1], [0, 2], [1, 3], [2, 2]]); "
        "forest = lakes.build_lake_hierarchy(elev, is_ocean, nb); "
        "assert len(forest) == 1; "
        "comps = hydrology.lake_components(np.array([True, False, False, False]), nb); "
        "assert len(comps) == 1 and comps[0].tolist() == [0]"
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


def test_stale_rust_module_missing_lake_exports_falls_back_cleanly(caplog):
    if _lakes_rust.RUST_AVAILABLE:
        assert _lakes_rust._resolve_rust_module(_lakes_rust._rust) is _lakes_rust._rust
    with caplog.at_level("WARNING", logger=_lakes_rust.__name__):
        resolved = _lakes_rust._resolve_rust_module(SimpleNamespace(step_ice=object()))
    assert resolved is None
    assert "bin/build_rust.py" in caplog.text


@pytest.mark.skipif(not _lakes_rust.RUST_AVAILABLE, reason="mantle_bloom_rust lake exports unavailable")
def test_native_entry_points_handle_fortran_and_strided_arrays():
    import math

    import mantle_bloom_rust

    # Fortran-order neighbor_idx and interface_pass (Codex review finding)
    nb_c = np.array([[1, 1], [0, 0], [2, 2]], dtype=np.int64)
    nb_f = np.asfortranarray(nb_c)
    is_lake = np.array([True, True, True], dtype=bool)
    ip_c, mem_c = mantle_bloom_rust.lake_components_rs(is_lake, nb_c)
    ip_f, mem_f = mantle_bloom_rust.lake_components_rs(is_lake, nb_f)
    np.testing.assert_array_equal(ip_f, ip_c)
    np.testing.assert_array_equal(mem_f, mem_c)

    elev = np.array([-0.0, 5.0, 0.0, -10.0], dtype=np.float64)
    is_ocean = np.array([False, False, False, True], dtype=bool)
    nb2_c = np.array([[1, 2], [0, 2], [1, 3], [2, 2]], dtype=np.int64)
    nb2_f = np.asfortranarray(nb2_c)
    pass_c = np.array([[5.0, 2.0], [5.0, 3.0], [3.0, 1.0], [1.0, -10.0]], dtype=np.float64)
    pass_f = np.asfortranarray(pass_c)

    np.testing.assert_array_equal(
        mantle_bloom_rust.catchment_roots_rs(elev, is_ocean, nb2_f, lakes._OCEAN_CATCHMENT),
        mantle_bloom_rust.catchment_roots_rs(elev, is_ocean, nb2_c, lakes._OCEAN_CATCHMENT),
    )

    ev_c = mantle_bloom_rust.lake_hierarchy_events_rs(elev, is_ocean, nb2_c, pass_c, lakes._OCEAN_CATCHMENT)
    ev_f = mantle_bloom_rust.lake_hierarchy_events_rs(elev, is_ocean, nb2_f, pass_f, lakes._OCEAN_CATCHMENT)
    for a, b in zip(ev_c, ev_f):
        np.testing.assert_array_equal(a, b)

    # Signed-zero endpoint maximum matches np.maximum(-0.0, +0.0) -> +0.0
    ev_z = mantle_bloom_rust.lake_hierarchy_events_rs(
        np.array([-0.0, 0.0]), np.array([False, True]), np.array([[1], [0]], dtype=np.int64), None, -2,
    )
    assert math.copysign(1.0, float(ev_z[4][0])) == 1.0


@pytest.mark.skipif(not _lakes_rust.RUST_AVAILABLE, reason="mantle_bloom_rust lake exports unavailable")
def test_lake_components_skips_dry_row_neighbors(monkeypatch):
    # Matches Python reference: dry rows' out-of-range neighbors are ignored when only wet rows matter.
    is_lake = np.array([True, False], dtype=bool)
    nb = np.array([[0], [-1]], dtype=np.int64)
    rust_groups = hydrology.lake_components(is_lake, nb)
    monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", False)
    py_groups = hydrology.lake_components(is_lake, nb)
    monkeypatch.setattr(_lakes_rust, "RUST_AVAILABLE", True)
    assert len(rust_groups) == len(py_groups) == 1
    np.testing.assert_array_equal(rust_groups[0], py_groups[0])

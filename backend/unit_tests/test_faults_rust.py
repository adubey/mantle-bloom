"""Parity for fault trace distance calculations (#1046)."""

import numpy as np
import pytest

from app import _faults_rust


def test_nearest_trace_matches_numpy_and_first_tie(monkeypatch):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "nearest_trace_rs"):
        pytest.skip("Rust extension unavailable")
    rng = np.random.default_rng(1046)
    points = np.vstack((rng.normal(size=(200, 3)), [0.0, 0.0, 0.0]))
    trace = np.vstack(([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], rng.normal(size=(12, 3)) * 3))
    native = _faults_rust.nearest_trace(points, trace)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    fallback = _faults_rust.nearest_trace(points, trace)
    np.testing.assert_array_equal(native[0], fallback[0])
    np.testing.assert_array_equal(native[1], fallback[1])
    assert native[1][-1] == 0


def test_shear_backtrace_matches_python(monkeypatch):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "shear_backtrace_rs"):
        pytest.skip("Rust extension unavailable")
    rng = np.random.default_rng(1046)
    points = rng.normal(size=(500, 3))
    points /= np.linalg.norm(points, axis=1)[:, None]
    trace = points[:8].copy()
    dip = np.array([0.0, 0.0, 1.0])
    native = _faults_rust.shear_backtrace(points, trace, dip, 0.9, 0.15, -1)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    fallback = _faults_rust.shear_backtrace(points, trace, dip, 0.9, 0.15, -1)
    np.testing.assert_array_equal(native[0], fallback[0])
    np.testing.assert_allclose(native[1], fallback[1], rtol=0, atol=3e-16)


@pytest.mark.parametrize("floor", [0.0, 0.12, 1.2, np.nan])
def test_influence_matches_numpy(monkeypatch, floor):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "influence_rs"):
        pytest.skip("Rust extension unavailable")
    distances = np.array([0.0, 0.1, 0.2, 1.0, np.inf, np.nan])
    native = _faults_rust.influence(distances, 0.2, floor)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    np.testing.assert_array_equal(native, _faults_rust.influence(distances, 0.2, floor))


def test_foreign_nearest_ties_and_empty_neighborhoods(monkeypatch):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "foreign_nearest_rs"):
        pytest.skip("Rust extension unavailable")
    points = np.array([[0., 0., 0.], [1., 0., 0.], [-1., 0., 0.], [5., 0., 0.]])
    outline = np.array([[0., 0., 0.], [3., 0., 0.], [0., 0., 0.]])
    owners = np.array([0, 1, 2, 0])
    offsets = np.array([0, 3, 4, 4])
    ids = np.array([0, 1, 2, 3])
    args = (points, outline, owners, offsets, ids, 0)
    native = _faults_rust.foreign_nearest(*args)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    fallback = _faults_rust.foreign_nearest(*args)
    for actual, expected in zip(native, fallback, strict=True):
        np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(native[1], [1, -1, -1])


def test_shear_full_plate_matches_fallback(monkeypatch):
    import copy

    from app import faults, geometry
    from app.surface_fields import SURFACE_FIELDS
    from app.world import World
    from .quad_fixtures import quad_plate
    from .test_faults import _fault

    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "shear_backtrace_rs"):
        pytest.skip("Rust extension unavailable")
    rows, columns = np.arange(-20, 21, 10), np.arange(-50, 51)
    rng = np.random.default_rng(1046)
    plate = quad_plate(0, "continental", columns=columns, rows=rows,
                       n=int(round((np.pi / 2) / 0.001)),
                       elevation=rng.uniform(100, 200, size=len(rows) * len(columns)))
    dip = geometry.normalize(geometry.local_xyz(np.array([0.]), np.array([1e-5]))[0]
                             - geometry.local_xyz(np.array([0.]), np.array([0.]))[0])
    fault = _fault(kind=faults._KIND_STRIKE_SLIP, local_phi=np.linspace(-0.03, 0.03, 7),
                   local_theta=np.zeros(7), dip_dir_local=dip, strike_sense=1,
                   slip_rate_m_per_myr=faults.SLIP_RATE_MAX_M_PER_MYR, plate_id=0)
    world = World(seed=0, plates=[plate])
    world.faults = [fault]
    world.fault_deformation_mode = "boundary"
    baseline = copy.deepcopy(world)
    faults._apply_plate_fault_shear(world, plate, 10.0)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    faults._apply_plate_fault_shear(baseline, baseline.plates[0], 10.0)
    for field in SURFACE_FIELDS:
        np.testing.assert_array_equal(plate.collect(field), baseline.plates[0].collect(field), err_msg=field)


@pytest.mark.parametrize("rates", [(-180., 55., 1.), (1., -1., 130.), (-1., 1., 130.)])
def test_relief_contribution_matches_python(monkeypatch, rates):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "relief_contrib_rs"):
        pytest.skip("Rust extension unavailable")
    rng = np.random.default_rng(1046)
    points = rng.normal(size=(500, 3))
    distances = rng.uniform(0, 2, len(points))
    mid = rng.normal(size=3)
    dip = rng.normal(size=3)
    params = (1., 0.7, 0.1, 1., 0.3, *rates)
    native = _faults_rust.relief_contrib(points, distances, mid, dip, params)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    fallback = _faults_rust.relief_contrib(points, distances, mid, dip, params)
    np.testing.assert_array_equal(native, fallback)


@pytest.mark.parametrize("kind", ["normal", "strike_slip"])
def test_relief_full_plate_matches_fallback(monkeypatch, kind):
    import copy

    from app import faults, geometry
    from app.surface_fields import SURFACE_FIELDS
    from app.world import World
    from .quad_fixtures import quad_plate
    from .test_faults import _fault

    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "relief_contrib_rs"):
        pytest.skip("Rust extension unavailable")
    rows, columns = np.arange(-20, 21, 10), np.arange(-50, 51)
    rng = np.random.default_rng(1046)
    count = len(rows) * len(columns)
    plate = quad_plate(0, "continental", columns=columns, rows=rows,
                       n=int(round((np.pi / 2) / 0.001)),
                       crustal_thickness_m=rng.uniform(32000, 48000, count),
                       continental_material_m=rng.uniform(10000, 20000, count))
    dip = geometry.normalize(geometry.local_xyz(np.array([0.]), np.array([1e-5]))[0]
                             - geometry.local_xyz(np.array([0.]), np.array([0.]))[0])
    fault = _fault(kind=kind, local_phi=np.linspace(-0.03, 0.03, 7), local_theta=np.zeros(7),
                   dip_dir_local=dip, slip_rate_m_per_myr=30000, plate_id=0)
    world = World(seed=0, plates=[plate])
    world.faults = [fault]
    baseline = copy.deepcopy(world)
    faults._apply_plate_fault_relief(world, plate, 1.0)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    faults._apply_plate_fault_relief(baseline, baseline.plates[0], 1.0)
    for field in SURFACE_FIELDS:
        actual, expected = plate.collect(field), baseline.plates[0].collect(field)
        if np.issubdtype(actual.dtype, np.floating):
            np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=1e-9, err_msg=field)
        else:
            np.testing.assert_array_equal(actual, expected, err_msg=field)


@pytest.mark.parametrize("count", [0, 1, 19, 1000])
def test_relief_transfer_parity_and_conservation(monkeypatch, count):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "relief_transfer_rs"):
        pytest.skip("Rust extension unavailable")
    rng = np.random.default_rng(1046 + count)
    hc = rng.uniform(5000, 80000, count)
    material = rng.uniform(0, 1, count) * hc
    areas = rng.uniform(1e6, 1e11, count)
    kinds = rng.integers(-1, 3, count, dtype=np.int8)
    request = rng.uniform(-40000, 40000, count)
    args = (hc, material, areas, kinds, request, 6000., 75000.)
    native = _faults_rust.relief_transfer(*args)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    fallback = _faults_rust.relief_transfer(*args)
    for a, b in zip(native, fallback, strict=True):
        np.testing.assert_allclose(a, b, rtol=3e-14, atol=1e-9)
    for kind in np.unique(kinds):
        mask = kinds == kind
        for before, after in zip((hc, material), native, strict=True):
            assert np.dot(before[mask], areas[mask]) == pytest.approx(np.dot(after[mask], areas[mask]), rel=2e-14)


@pytest.mark.parametrize("regime,owners", [([], []), ([0, 0], [1, 1]),
    ([1, 1, 1, 0, 2, 2, 3], [1, 1, 2, 2, 2, 2, 2])])
def test_boundary_runs_preserve_order_and_splits(monkeypatch, regime, owners):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "boundary_runs_rs"):
        pytest.skip("Rust extension unavailable")
    regime = np.array(regime, dtype=np.int8)
    owners = np.array(owners, dtype=np.int64)
    native = _faults_rust.boundary_runs(regime, owners)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    fallback = _faults_rust.boundary_runs(regime, owners)
    for a, b in zip(native, fallback, strict=True):
        np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("sense", [-1., 0., 1.])
def test_unloading_stress_parity(monkeypatch, sense):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "unloading_stress_rs"):
        pytest.skip("Rust extension unavailable")
    load = np.array([-1e9, -1e6, -1., 0., 1e6])
    args = (load, sense, 1 / 3, np.sin(np.pi / 6), np.cos(np.pi / 6), 0.4)
    native = _faults_rust.unloading_stress(*args)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    np.testing.assert_array_equal(native, _faults_rust.unloading_stress(*args))


@pytest.mark.parametrize("overlap", [False, True])
def test_quake_magnitude_parity(monkeypatch, overlap):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "quake_magnitude_rs"):
        pytest.skip("Rust extension unavailable")
    constants = (3.4, 1.15, 0.35, 0.4, 4., 9.3)
    args = (137., 4500., overlap, 0.19, constants)
    native = _faults_rust.quake_magnitude(*args)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    assert native == pytest.approx(_faults_rust.quake_magnitude(*args), abs=2e-15)


def test_boundary_regime_thresholds_and_nan(monkeypatch):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "boundary_regime_rs"):
        pytest.skip("Rust extension unavailable")
    distances = np.array([0., 1., 1., 1., 1.1, np.nan, 0.])
    closing = np.array([0., 1., -1., 0.5, 2., 2., np.nan])
    native = _faults_rust.boundary_regime(distances, closing, 1., 0.5)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    np.testing.assert_array_equal(native, _faults_rust.boundary_regime(distances, closing, 1., 0.5))
    np.testing.assert_array_equal(native, [3, 1, 2, 3, 0, 0, 0])


def test_moment_magnitude_parity(monkeypatch):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "moment_magnitude_rs"):
        pytest.skip("Rust extension unavailable")
    args = (153000., 0.2, 15000., 3e10, 4., 9.3)
    native = _faults_rust.moment_magnitude(*args)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    assert native == pytest.approx(_faults_rust.moment_magnitude(*args), abs=2e-15)


def test_unloading_earthquake_event_matches_fallback(monkeypatch):
    import copy
    from app import faults
    from .test_faults import _fault, _unloading_world

    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "moment_magnitude_rs"):
        pytest.skip("Rust extension unavailable")
    world, points, loads = _unloading_world(_fault(kind="reverse", dip_deg=30., fault_id=4),
                                             np.array([0., -1e6, -3e6, -1e6, 0.]))
    baseline = copy.deepcopy(world)
    native = faults.trigger_unloading_earthquakes(world, points, loads, 1e6)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    fallback = faults.trigger_unloading_earthquakes(baseline, points, loads, 1e6)
    assert len(native) == len(fallback) == 1
    assert native[0].fault_id == fallback[0].fault_id
    assert native[0].earthquake_id == fallback[0].earthquake_id
    np.testing.assert_array_equal(native[0].epicenter_world, fallback[0].epicenter_world)
    assert native[0].slip_m == pytest.approx(fallback[0].slip_m, rel=1e-15)
    assert native[0].magnitude == pytest.approx(fallback[0].magnitude, abs=2e-15)
    assert world.faults[0].cumulative_offset_m == pytest.approx(baseline.faults[0].cumulative_offset_m, rel=1e-15)


def test_tectonic_earthquake_events_match_fallback(monkeypatch):
    import copy
    from app import faults
    from .test_faults import _fault, _unloading_world

    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "quake_magnitude_rs"):
        pytest.skip("Rust extension unavailable")
    fault = _fault(kind="reverse", fault_id=4)
    world, points, _ = _unloading_world(fault, 0.)
    fault.world_polyline = points
    boundary = copy.deepcopy(fault)
    boundary.fault_id = 5
    boundary.boundary = True
    world.boundary_faults = [boundary]
    inactive = copy.deepcopy(fault)
    inactive.fault_id = 6
    inactive.active = False
    world.faults.append(inactive)
    baseline = copy.deepcopy(world)
    faults._generate_earthquakes(world, 1.)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    faults._generate_earthquakes(baseline, 1.)
    assert len(world.earthquakes) == len(baseline.earthquakes) == 2
    for a, b in zip(world.earthquakes, baseline.earthquakes, strict=True):
        assert (a.earthquake_id, a.fault_id, a.plate_id, a.kind) == (b.earthquake_id, b.fault_id, b.plate_id, b.kind)
        np.testing.assert_array_equal(a.epicenter_world, b.epicenter_world)
        assert a.magnitude == pytest.approx(b.magnitude, abs=2e-15)
        assert a.slip_m == b.slip_m
    assert world.next_earthquake_id == baseline.next_earthquake_id


def test_boundary_generation_matches_fallback(monkeypatch):
    import copy
    from app import faults
    from app.world import generate_world

    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "boundary_regime_rs"):
        pytest.skip("Rust extension unavailable")
    world = generate_world(seed=4, num_plates=6, node_density=0.5, climate_density=0.5, fluid_density=0.5)
    baseline = copy.deepcopy(world)
    faults.generate_boundary_faults(world)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    faults.generate_boundary_faults(baseline)
    assert len(world.boundary_faults) == len(baseline.boundary_faults) > 0
    for a, b in zip(world.boundary_faults, baseline.boundary_faults, strict=True):
        assert (a.fault_id, a.plate_id, a.kind, a.strike_sense) == (b.fault_id, b.plate_id, b.kind, b.strike_sense)
        assert a.slip_rate_m_per_myr == b.slip_rate_m_per_myr
        np.testing.assert_array_equal(a.world_polyline, b.world_polyline)
        np.testing.assert_array_equal(a.dip_dir_local, b.dip_dir_local)
    assert world.next_fault_id == baseline.next_fault_id


def test_boundary_resampling_parity(monkeypatch):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "resample_boundary_rs"):
        pytest.skip("Rust extension unavailable")
    rng = np.random.default_rng(1046)
    points = rng.normal(size=(99, 3))
    points /= np.linalg.norm(points, axis=1)[:, None]
    positions = np.linspace(0., len(points) - 1, 23)
    native = _faults_rust.resample_boundary(points, positions)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    np.testing.assert_array_equal(native, _faults_rust.resample_boundary(points, positions))


def test_disable_rust_override_in_fresh_process():
    import os
    import subprocess
    import sys

    env = os.environ.copy()
    env["MANTLE_BLOOM_DISABLE_RUST"] = "1"
    env["PYTHONPATH"] = "backend"
    subprocess.run([sys.executable, "-c", "from app import _faults_rust; assert _faults_rust._rust is None"],
                   env=env, check=True, capture_output=True, text=True)


def test_boundary_strands_match_fallback(monkeypatch):
    if _faults_rust._rust is None or not hasattr(_faults_rust._rust, "boundary_strands_rs"):
        pytest.skip("Rust extension unavailable")
    rng = np.random.default_rng(1046)
    master = rng.normal(size=(25, 3))
    master /= np.linalg.norm(master, axis=1)[:, None]
    centroid = np.array([1., 0., 0.])
    args = (master, centroid, 3, 0.01)
    native = _faults_rust.boundary_strands(*args)
    monkeypatch.setattr(_faults_rust, "_rust", None)
    for a, b in zip(native, _faults_rust.boundary_strands(*args), strict=True):
        np.testing.assert_array_equal(a, b)

"""Equivalence tests verifying `mantle_bloom_rust` matches Python/Numba glacier and sea-level solves."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.spatial import cKDTree

from app import lithosphere, world as world_mod
from app.elevation_lines import PLANET_RADIUS_KM
from app.hydroclimate import _glacier_rust, eustasy, hydrology, ice_sheet


def _make_sphere_patch(spacing_km: float = 120.0, half_km: float = 1200.0):
    r_m = PLANET_RADIUS_KM * 1000.0
    m = int(half_km / spacing_km) + 1
    g = np.arange(-m, m + 1, dtype=float)
    gx, gy = (a.ravel() for a in np.meshgrid(g, g))
    x = (gx + 0.5 * (gy % 2)) * spacing_km * 1e3
    y = gy * (np.sqrt(3.0) / 2.0) * spacing_km * 1e3
    pts = np.stack([np.ones_like(x), x / r_m, y / r_m], axis=1)
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    nb = cKDTree(pts).query(pts, hydrology.FLOW_NEIGHBOR_COUNT + 1)[1][:, 1:]
    area = np.full(len(pts), spacing_km**2 * (np.sqrt(3.0) / 2.0) * 1e6)
    return pts, nb, x, y, area


def _assert_close(actual: np.ndarray, expected: np.ndarray, rtol: float = 1e-12, atol: float = 1e-12):
    np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True)


def test_rust_extension_is_built_and_available():
    assert _glacier_rust.RUST_AVAILABLE is True


def test_step_ice_rust_matches_python_across_bodies_and_edge_cases(monkeypatch):
    pts, nb, x, y, area = _make_sphere_patch(spacing_km=140.0, half_km=1100.0)
    n = len(pts)
    r1 = np.hypot(x + 420e3, y)
    r2 = np.hypot(x - 520e3, y - 300e3)

    bed = np.where(r1 < 350e3, 1200.0, 200.0) - np.where(x > 750e3, 2500.0, 0.0)
    defl = np.linspace(-20.0, 50.0, n)
    elev = bed + defl
    is_ocean = bed < -500.0
    prev_depth = np.where(r1 < 220e3, 900.0, 0.0)
    # Include a 1-node isolated body (< ICE_BODY_MIN_NODES) and a floating stray node.
    prev_depth[0] = 50.0
    is_ocean[0] = False
    ocean_nodes = np.flatnonzero(is_ocean)
    if len(ocean_nodes) > 0:
        prev_depth[ocean_nodes[0]] = 50.0

    temp = np.where((r1 < 420e3) | (r2 < 260e3), -18.0, 6.0)
    temp[0] = -12.0
    amp = np.where(np.arange(n) % 7 == 0, 0.0, 12.0)
    precip = np.where(r1 < 500e3, 600.0, 200.0)
    fric = np.linspace(0.0, 0.04, n)
    dry_snow = np.where(temp < 0.0, precip * 0.8, 0.0)
    bed_frozen = np.where(np.arange(n) % 11 == 0, 25.0, 0.0)
    mobile = np.linspace(0.0, 30.0, n)
    hints = np.where(r1 < 250e3, 150_000.0, np.nan)

    climate = ice_sheet.SurfaceClimate(
        precip, temp, amp, temp < 0.0,
        frictional_heat_w_m2=fric, dry_snow_mm=dry_snow, bed_frozen_snow_mm=bed_frozen,
    )
    cost = ice_sheet.march_cost("cheap-substep-limited")

    d_rust, m_rust, ice_rust = ice_sheet.step_ice(
        pts, nb, area, elev, defl, 0.0, prev_depth, climate, mobile, is_ocean, 2_000.0,
        margin_hint_m=hints, cost=cost,
    )

    monkeypatch.setattr(_glacier_rust, "RUST_AVAILABLE", False)
    d_py, m_py, ice_py = ice_sheet.step_ice(
        pts, nb, area, elev, defl, 0.0, prev_depth, climate, mobile, is_ocean, 2_000.0,
        margin_hint_m=hints, cost=cost,
    )

    _assert_close(d_rust, d_py)
    _assert_close(m_rust, m_py)
    _assert_close(ice_rust.thickness, ice_py.thickness)
    _assert_close(ice_rust.covered, ice_py.covered)
    np.testing.assert_array_equal(ice_rust.warm_based, ice_py.warm_based)
    _assert_close(ice_rust.balance_flux_m3_per_yr, ice_py.balance_flux_m3_per_yr)
    _assert_close(ice_rust.sliding_speed_m_per_yr, ice_py.sliding_speed_m_per_yr)
    _assert_close(ice_rust.calved_m3, ice_py.calved_m3)
    assert len(ice_rust.bodies) == len(ice_py.bodies)
    for br, bp in zip(ice_rust.bodies, ice_py.bodies):
        np.testing.assert_array_equal(br.nodes, bp.nodes)
        np.testing.assert_array_equal(br.domain_nodes, bp.domain_nodes)
        assert br.margin_m == pytest.approx(bp.margin_m, rel=1e-12, abs=1e-12)
        assert br.closing_m3 == pytest.approx(bp.closing_m3, rel=1e-12, abs=1e-6)
    assert ice_rust.ledger.residual_m3 == pytest.approx(0.0, abs=1e-12 * max(ice_rust.ledger.closing_m3, 1.0))
    assert ice_rust.ledger.closing_m3 == pytest.approx(ice_py.ledger.closing_m3, rel=1e-12, abs=1e-6)
    assert ice_rust.frozen_precipitation.retained_m3 == pytest.approx(
        ice_py.frozen_precipitation.retained_m3, rel=1e-12, abs=1e-6
    )


def test_eustasy_and_ocean_mask_rust_match_python(monkeypatch):
    pts, nb, x, y, area = _make_sphere_patch(spacing_km=140.0, half_km=900.0)
    r = np.hypot(x, y)
    elev = np.where(r < 400e3, -200.0, np.where(r < 650e3, 300.0, -500.0))
    elev_nan = elev.copy()
    elev_nan[3] = np.nan
    seed = r > 700e3
    ice_water = np.where((r > 500e3) & (r < 620e3), 150.0, 0.0)
    prior_ocean = r > 680e3

    w = world_mod.World(seed=1, node_density=1.0)
    budget = float(np.dot(area, np.clip(-50.0 - elev, 0.0, None))) * 0.8

    levels_rust = eustasy._seeded_flood_levels(elev_nan, nb, seed)
    sl_rust = eustasy._solve_sea_level_connected(
        w, elev, area, budget, seed_is_ocean=seed, seed_neighbor_idx=nb,
        bracket=(-400.0, 100.0), ice_water_m=ice_water,
    )
    comp_rust = hydrology._below_sea_components(elev <= -100.0, nb)
    mask_rust = hydrology.connected_ocean_mask(
        pts, elev, -100.0, nb, areas_m2=area, prior_ocean=prior_ocean,
    )

    monkeypatch.setattr(_glacier_rust, "RUST_AVAILABLE", False)
    levels_py = eustasy._seeded_flood_levels(elev_nan, nb, seed)
    sl_py = eustasy._solve_sea_level_connected(
        w, elev, area, budget, seed_is_ocean=seed, seed_neighbor_idx=nb,
        bracket=(-400.0, 100.0), ice_water_m=ice_water,
    )
    comp_py = hydrology._below_sea_components(elev <= -100.0, nb)
    mask_py = hydrology.connected_ocean_mask(
        pts, elev, -100.0, nb, areas_m2=area, prior_ocean=prior_ocean,
    )

    np.testing.assert_array_equal(levels_rust, levels_py)
    assert sl_rust == pytest.approx(sl_py, rel=1e-12, abs=1e-12)
    np.testing.assert_array_equal(comp_rust, comp_py)
    np.testing.assert_array_equal(mask_rust, mask_py)


def test_march_helpers_rust_match_python():
    n = 64
    rng = np.random.default_rng(42)
    elev = rng.uniform(-1000.0, 3000.0, n)
    defl = rng.uniform(-200.0, 0.0, n)
    depth = rng.uniform(0.0, 1500.0, n)
    ocean = elev < 0.0
    hc = rng.choice([0.0, 11_000.0, 35_000.0], size=n)
    hm = np.full(n, 80_000.0)
    rho = np.full(n, lithosphere.RHO_CONTINENTAL_CRUST)
    temp0 = rng.uniform(-25.0, 15.0, n)

    t_rust, f_rust = _glacier_rust.march_substep_prep(elev, defl, depth, 0.0, ocean, hc, hm, rho, temp0)
    unloaded = elev - defl
    load = lithosphere.grounded_ice_load_kg_m2(depth, unloaded, 0.0, ocean)
    target = lithosphere.ice_load_deflection(hc, hm, rho, load)
    t_py = ice_sheet.air_temperature_at_elevation(temp0, defl, target)
    f_py = (t_py < hydrology.FREEZE_POINT_C) | (depth >= hydrology.GLACIER_VISIBLE_DEPTH_M)
    _assert_close(t_rust, t_py)
    np.testing.assert_array_equal(f_rust, f_py)

    sliding = rng.uniform(0.0, 40.0, n)
    cover_r, cover_p = rng.uniform(0.0, 20.0, n), None
    cover_p = cover_r.copy()
    head_r = rng.uniform(0.0, 100.0, n)
    head_p = head_r.copy()
    craton = rng.uniform(0.2, 1.0, n)
    raw_r, raw_p = np.zeros(n), np.zeros(n)
    abr_r, abr_p = np.zeros(n), np.zeros(n)
    we_r, we_p = elev.copy(), elev.copy()

    cut_r = _glacier_rust.apply_march_abrasion(sliding, 250.0, 1.5, ocean, cover_r, head_r, craton, raw_r, abr_r, we_r)
    monkeypatch_guard = _glacier_rust.RUST_AVAILABLE
    _glacier_rust.RUST_AVAILABLE = False
    try:
        cut_p = hydrology._apply_march_abrasion(sliding, 250.0, 1.5, ocean, cover_p, head_p, craton, raw_p, abr_p, we_p)
    finally:
        _glacier_rust.RUST_AVAILABLE = monkeypatch_guard

    assert cut_r == pytest.approx(cut_p, rel=1e-12, abs=1e-12)
    _assert_close(cover_r, cover_p)
    _assert_close(head_r, head_p)
    _assert_close(raw_r, raw_p)
    _assert_close(abr_r, abr_p)
    _assert_close(we_r, we_p)


def test_step_world_rust_matches_python_at_coarse_density(monkeypatch):
    w_rust = world_mod.generate_world(seed=0, node_density=0.5, climate_density=0.5, fluid_density=0.5)
    w_rust.wind_model = "diagnostic"
    for _ in range(2):
        world_mod.step_world(w_rust, 100_000.0)

    monkeypatch.setattr(_glacier_rust, "RUST_AVAILABLE", False)
    w_py = world_mod.generate_world(seed=0, node_density=0.5, climate_density=0.5, fluid_density=0.5)
    w_py.wind_model = "diagnostic"
    for _ in range(2):
        world_mod.step_world(w_py, 100_000.0)

    assert w_rust.sea_level_m == pytest.approx(w_py.sea_level_m, rel=1e-9, abs=1e-6)
    for pr, pp in zip(w_rust.plates, w_py.plates):
        _assert_close(pr.collect("elevation"), pp.collect("elevation"), rtol=1e-8, atol=1e-5)
        _assert_close(pr.collect("glacier_depth"), pp.collect("glacier_depth"), rtol=1e-8, atol=1e-5)

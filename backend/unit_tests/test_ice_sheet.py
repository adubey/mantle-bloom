import numpy as np
import pytest
from scipy.spatial import cKDTree

from app import lithosphere
from app.elevation_lines import PLANET_RADIUS_KM, line_spacing_rad
from app.hydroclimate import hydrology, ice_sheet

_R_M = PLANET_RADIUS_KM * 1000.0
_C_100KPA = 100_000.0 / (lithosphere.RHO_ICE * lithosphere.GRAVITY_M_S2)


def _lattice(spacing_km, half_km, kind="square"):
    """A flat patch of nodes near (lat 0, lon 0): positions, (x, y) in metres, cell area."""
    m = int(half_km / spacing_km) + 2
    g = np.arange(-m, m + 1, dtype=float)
    x, y = (a.ravel() for a in np.meshgrid(g, g))
    cell = spacing_km**2
    if kind == "hex":
        x = x + 0.5 * (y % 2)
        y = y * np.sqrt(3.0) / 2.0
        cell *= np.sqrt(3.0) / 2.0
    x, y = x * spacing_km * 1e3, y * spacing_km * 1e3
    points = np.stack([np.ones_like(x), x / _R_M, y / _R_M], axis=1)
    points /= np.linalg.norm(points, axis=1)[:, None]
    return points, x, y, cell * 1e6


def _neighbors(points):
    return cKDTree(points).query(points, hydrology.FLOW_NEIGHBOR_COUNT + 1)[1][:, 1:]


def _cold(n, precipitation=500.0, temperature=-30.0, amplitude=0.0):
    return ice_sheet.SurfaceClimate(
        np.broadcast_to(np.asarray(precipitation, dtype=float), (n,)).copy(),
        np.broadcast_to(np.asarray(temperature, dtype=float), (n,)).copy(),
        np.broadcast_to(np.asarray(amplitude, dtype=float), (n,)).copy(),
        np.broadcast_to(np.asarray(temperature, dtype=float), (n,)) < hydrology.FREEZE_POINT_C,
    )


def _dome_domain(points, radius_m, bed, cell_area, sea_level=-1e4):
    """Every node as one domain whose positive-balance area is a disc of `radius_m`, with
    phi the true (Euclidean) signed distance, so the profile is tested on its own."""
    n = len(points)
    graph = ice_sheet.symmetric_graph(points, _neighbors(points))
    r = np.hypot(points[:, 1], points[:, 2]) * _R_M
    return ice_sheet._build_domain(
        graph, np.arange(n), r - radius_m, bed, np.full(n, _C_100KPA), np.full(n, cell_area), sea_level,
        np.zeros(n, dtype=bool), _cold(n),
    )


@pytest.mark.parametrize("density", [0.5, 1.0, 2.0, 4.0])
@pytest.mark.parametrize("kind", ["square", "hex"])
def test_flat_bed_dome_matches_nye_profile(density, kind):
    # Nye's perfectly plastic dome on a flat bed: H = sqrt(2 c (L - r)), so the divide is
    # sqrt(2 c L) thick and the volume is (8 pi / 15) sqrt(2 c) L^2.5. The graph solve must
    # match both within 5% at every node density and on two lattice layouts.
    spacing = line_spacing_rad(density) * PLANET_RADIUS_KM
    L = 600e3
    points, _, _, cell = _lattice(spacing, L / 1e3 + 2 * spacing, kind)
    shape = ice_sheet._shape(_dome_domain(points, L, np.zeros(len(points)), cell), 0.0, -1e4)
    assert shape.thickness.max() == pytest.approx(np.sqrt(2 * _C_100KPA * L), rel=0.05)
    assert shape.volume == pytest.approx(8 * np.pi / 15 * np.sqrt(2 * _C_100KPA) * L**2.5, rel=0.05)


def test_upwind_update_is_monotone_causal_and_regular_at_the_margin():
    # ds_i/ds_up in [0, 1), s_i > s_up, and no singularity starting from H = 0.
    s_up = np.linspace(-200.0, 3000.0, 400)
    s = np.array([ice_sheet._upwind(v, 0.0, 100.0, _C_100KPA, 50e3) for v in s_up])
    assert np.all(np.diff(s) >= 0.0) and np.all(np.diff(s) < np.diff(s_up))
    assert np.all(s > s_up)
    assert ice_sheet._upwind(0.0, 0.0, 0.0, _C_100KPA, 50e3) == pytest.approx(np.sqrt(2 * _C_100KPA * 50e3))


def test_extent_family_volume_is_continuous_nondecreasing_and_never_lowers_the_surface():
    rng = np.random.default_rng(3)
    points, _, _, cell = _lattice(125.0, 900.0, "hex")
    n = len(points)
    dom = _dome_domain(points, 400e3, rng.uniform(-100.0, 100.0, n), cell)
    previous = None
    volumes = []
    for d in np.arange(-450e3, 300e3, 2e3):
        shape = ice_sheet._shape(dom, d, -1e4)
        volumes.append(shape.volume)
        if previous is not None:
            was = previous.active
            assert np.all(shape.active[was])
            assert np.all(shape.surface[was] >= previous.surface[was] - 1e-6)
        previous = shape
    volumes = np.array(volumes)
    assert np.all(np.diff(volumes) >= -1e-9 * volumes.max())
    # No jumps: a 2 km move of the margin never adds more than a percent of the sheet.
    assert np.max(np.diff(volumes)) < 0.01 * volumes[-1]


def test_a_bed_bump_inside_the_sheet_pokes_through_as_a_nunatak():
    points, x, y, cell = _lattice(62.5, 800.0)
    n = len(points)
    bed = np.zeros(n)
    bump = int(np.argmin(np.hypot(x - 200e3, y)))
    bed[bump] = 4000.0  # above the ~3 km surface a flat bed would carry there
    flat = ice_sheet._shape(_dome_domain(points, 600e3, np.zeros(n), cell), 0.0, -1e4)
    shape = ice_sheet._shape(_dome_domain(points, 600e3, bed, cell), 0.0, -1e4)
    assert np.all(np.isfinite(shape.surface[shape.active]))
    assert shape.thickness[bump] == 0.0
    ring = cKDTree(points).query(points[bump], 9)[1][1:]
    assert np.all(shape.thickness[ring] > 1000.0)
    # Ice around the bump stays; only the bump's own column (and a little shoulder) is lost.
    assert shape.volume == pytest.approx(flat.volume, rel=0.05)


def _continent(spacing_km=125.0, land_radius_km=900.0, half_km=1300.0, ocean_floor=-3000.0):
    points, x, y, cell = _lattice(spacing_km, half_km)
    r = np.hypot(x, y)
    elevation = np.where(r <= land_radius_km * 1e3, 500.0, ocean_floor)
    return points, r, elevation, np.full(len(points), cell)


def _step(points, area, elevation, depth, climate, years, sea_level=0.0, is_ocean=None):
    n = len(points)
    if is_ocean is None:
        is_ocean = elevation <= sea_level
    return ice_sheet.step_ice(
        points, _neighbors(points), area, elevation, np.zeros(n), sea_level, depth, climate, np.zeros(n), is_ocean,
        years,
    )


@pytest.mark.parametrize("min_nodes", [1, 100000], ids=["plastic", "small"])
def test_bare_snow_on_ocean_inside_glacier_domain_is_not_booked(monkeypatch, min_nodes):
    monkeypatch.setattr(ice_sheet, "ICE_BODY_MIN_NODES", min_nodes)
    points, r, elevation, area = _continent(ocean_floor=-10.0)
    ocean = elevation <= 0.0
    climate = _cold(len(points), temperature=np.where(r < 400e3, -20.0, 2.0), amplitude=10.0)
    melt_bare_snow = ice_sheet._melt_bare_snow

    def check_bare_snow(out, idx, *args):
        assert not np.any(ocean[idx])
        return melt_bare_snow(out, idx, *args)

    monkeypatch.setattr(ice_sheet, "_melt_bare_snow", check_bare_snow)
    _, _, ice = _step(points, area, elevation, np.zeros(len(points)), climate, 1000.0, is_ocean=ocean)
    if min_nodes == 1:
        domain = ice.bodies[0].domain_nodes
        assert np.any(ocean[domain])
    assert np.any(climate.rates(np.zeros(len(points)))[0][ocean] > 0.0)
    assert np.all(ice.bare_snowmelt_depth_m[ocean] == 0.0)
    assert np.dot(ice.bare_snowmelt_depth_m, area) == pytest.approx(ice.ledger.snow_melted_m3)


def _graded_climate(r, precipitation=500.0, inner_km=500.0):
    """Cold inside `inner_km`, warming outward past freezing, with a 10 C seasonal swing."""
    temperature = -25.0 + 30.0 * np.clip(r / (inner_km * 1e3), 0.0, 2.0) - 15.0
    return _cold(len(r), precipitation, temperature, 10.0)


def test_volume_contract_closes_per_body_and_across_the_world():
    points, r, elevation, area = _continent()
    climate = _graded_climate(r)
    depth = np.zeros(len(points))
    for _ in range(3):
        opening = float(np.dot(depth, area))
        depth, melt, ice = _step(points, area, elevation, depth, climate, 1e6)
        ledger = ice.ledger
        scale = max(ledger.accumulation_m3, ledger.closing_m3, 1.0)
        assert ledger.opening_m3 == pytest.approx(opening)
        assert abs(ledger.residual_m3) < 1e-9 * scale
        assert ledger.closing_m3 == pytest.approx(float(np.dot(depth, area)))
        for body in ice.bodies:
            assert abs(body.residual_m3) < 1e-9 * max(body.accumulation_m3, body.closing_m3, 1.0)
        # Ice is never created: closing <= opening + accumulation.
        assert ledger.closing_m3 <= ledger.opening_m3 + ledger.accumulation_m3 + 1e-9 * scale
        # Every booked sink reaches the water or the ocean: melt (with snow on bare ground) as
        # meltwater at its node, calving unrouted.
        assert float(np.dot(melt, area)) == pytest.approx(ledger.meltwater_m3, rel=1e-9, abs=1.0)
        bare_snowmelt = ice.bare_snowmelt_depth_m
        assert np.all((bare_snowmelt >= 0.0) & (bare_snowmelt <= melt))
        assert np.dot(bare_snowmelt, area) == pytest.approx(ledger.snow_melted_m3)
        ice_release = ice.ice_release_depth_m
        # Rust keeps each source in its own array; their sum can differ by a few ulps
        # from the combined accumulation while preserving the D8 source budget (#1019).
        np.testing.assert_allclose(melt, ice_release + bare_snowmelt, rtol=2e-15, atol=1e-12)
        assert np.dot(ice_release, area) == pytest.approx(ledger.meltwater_m3 - ledger.snow_melted_m3)
        assert np.dot(ice.surface_melt_depth_m, area) == pytest.approx(ledger.surface_melt_m3)
        assert np.all(ice.surface_melt_depth_m <= ice_release)
        assert float(ice.calved_m3.sum()) == pytest.approx(ledger.calved_m3)
    assert ledger.closing_m3 > 0.0


def _steady_state(points, area, elevation, climate, years=1e6, start=None, steps=4):
    """Settle the ice toward its balance. Physical rates (RFC 0015 D1, #789) respond in kyr,
    so Myr steps reach the steady sheet; the old 1e10 yr steps predate physical units."""
    depth = np.zeros(len(points)) if start is None else start
    for _ in range(steps):
        depth, _, ice = _step(points, area, elevation, depth, climate, years)
    return depth, ice


def test_more_snowfall_on_the_same_positive_balance_area_builds_a_bigger_sheet():
    # Snow-covered and melt-free inside 500 km, strong summer melt outside it: the positive-
    # balance area is the same disc whatever the snowfall. Twice the snow must still give a
    # larger sheet that reaches further out (accumulation magnitude, not just its sign).
    # Physical ice (RFC 0015 D1, #789) responds in kyr and pins at its margin within a few
    # of them, so the magnitude is read on one short step from bare ground, in the explicit
    # regime where the new volume is the opening balance times the step.
    points, r, elevation, area = _continent(land_radius_km=1100.0, half_km=1500.0)
    inside = r <= 500e3
    temperature = np.where(inside, -40.0, 8.0)

    def run(precipitation):
        climate = ice_sheet.SurfaceClimate(np.full(len(r), precipitation), temperature, np.full(len(r), 5.0), inside)
        depth, _, ice = _step(points, area, elevation, np.zeros(len(r)), climate, 500.0)
        return float(np.dot(depth, area)), np.count_nonzero(depth >= hydrology.GLACIER_VISIBLE_DEPTH_M), ice

    v1, n1, ice1 = run(300.0)
    v2, n2, ice2 = run(600.0)
    assert not any(b.saturated for b in ice1.bodies + ice2.bodies)
    assert v2 == pytest.approx(2.0 * v1, rel=0.02)
    assert n2 > n1


def test_backward_euler_converges_to_zero_balance_for_long_steps_and_to_explicit_for_short_ones():
    points, r, elevation, area = _continent(land_radius_km=1100.0, half_km=1500.0)
    climate = _graded_climate(r, inner_km=600.0)
    depth, ice = _steady_state(points, area, elevation, climate, years=1e6, steps=3)
    ledger = ice.ledger
    # dt * B = closing - opening, so with dt = 10 Gyr the balance is ~zero next to snowfall.
    assert not any(b.saturated for b in ice.bodies)
    assert abs(ledger.accumulation_m3 - ledger.surface_melt_m3 - ledger.basal_melt_m3) < 0.01 * ledger.accumulation_m3

    # Short steps from a growing sheet follow the explicit balance at the current extent:
    # one 200 yr step moves the opening balance (accumulation minus melt at the opening ice)
    # times the step. Physical ice (RFC 0015 D1, #789) responds in kyr, so Myr steps no
    # longer read as explicit.
    growing, _, _ = _step(points, area, elevation, np.zeros(len(points)), climate, 2e3)
    acc, sm, bm = climate.rates(growing)
    explicit_m3_per_yr = float(np.dot(acc - sm - bm, area)) / 1e6
    new, _, _ = _step(points, area, elevation, growing, climate, 200.0)
    rate_m3_per_yr = (float(np.dot(new, area)) - float(np.dot(growing, area))) / 200.0
    assert rate_m3_per_yr == pytest.approx(explicit_m3_per_yr, rel=0.05)


def test_height_mass_balance_feedback_gives_two_stable_sheets():
    # A high plateau whose bare surface melts more than it gets snow (a 1 C summer above
    # freezing), but whose surface under ~1 km of ice is cold enough to keep all its snow.
    # Bare, it stays bare; covered, it keeps a sheet (Weertman 1961).
    points, r, elevation, area = _continent(land_radius_km=1000.0, half_km=1400.0)
    elevation = np.where(elevation > 0.0, 2000.0, elevation)
    n = len(points)
    climate = ice_sheet.SurfaceClimate(np.full(n, 1000.0), np.full(n, -2.0), np.full(n, 8.0), np.ones(n, dtype=bool))
    bare, _ = _steady_state(points, area, elevation, climate, years=2e5, steps=3)
    assert bare.max() == 0.0
    start = np.where(elevation > 0.0, 1500.0, 0.0)
    covered, ice = _steady_state(points, area, elevation, climate, years=2e5, start=start, steps=3)
    assert float(np.dot(covered, area)) > 0.2 * float(np.dot(start, area))

    # Under fixed climate the chosen root moves continuously from the previous extent: a
    # short step barely moves the margin.
    margin = ice.bodies[0].margin_m
    _, _, ice = _step(points, area, elevation, covered, climate, 1e5)
    assert abs(ice.bodies[0].margin_m - margin) < 125e3


def test_grounded_ice_matches_the_ice_load_classification_and_flows_below_sea_level():
    # Land with a 150 m deep shelf around it: thick ice pushes its grounding line onto the shelf.
    points, r, elevation, area = _continent(land_radius_km=500.0)
    elevation = np.where((elevation < 0.0) & (r < 800e3), -150.0, elevation)
    start = np.where(r < 500e3, 2500.0, 0.0)
    # Physical ice (RFC 0015 D1, #789) advances the grounding line in kyr; a Myr step would
    # pile the interior past the depth guard and flatten the surface the test reads.
    depth, _, ice = _step(points, area, elevation, start, _cold(len(points)), 2e3)
    loaded = lithosphere.grounded_ice_load_kg_m2(depth, elevation, 0.0, elevation <= 0.0) > 0.0
    assert np.array_equal(depth > 0.0, loaded)
    shelf = (elevation < 0.0) & (depth > 0.0)
    assert np.any(shelf)
    neighbors = _neighbors(points)
    surface = elevation + depth
    assert np.all(np.max(np.abs(surface[neighbors[shelf]] - surface[shelf][:, None]), axis=1) > 0.0)
    assert np.any(ice.balance_flux_m3_per_yr[shelf] > 0.0)


def test_a_saturated_marine_body_calves_rather_than_melting():
    # An island in ocean too deep to ground on, too small to hold its ice: the excess calves at
    # its marine margin and is never meltwater.
    points, r, elevation, area = _continent(land_radius_km=200.0, half_km=600.0, ocean_floor=-5000.0)
    start = np.where(r <= 200e3, 4000.0, 0.0)
    depth, melt, ice = _step(points, area, elevation, start, _cold(len(points)), 1e6)
    assert ice.ledger.calved_m3 > 0.0
    assert ice.ledger.overflow_m3 == 0.0
    ledger = ice.ledger
    assert float(np.dot(melt, area)) == pytest.approx(ledger.meltwater_m3)
    # #355 step 5: the hydrology handoff is a yearly liquid-water rate. Calving leaves
    # directly for the ocean, even when it occurs in the same step as melt.
    routed_melt_m3 = hydrology.meltwater_m3_per_yr(melt, area, 1e6).sum() * 1e6
    assert routed_melt_m3 == pytest.approx(ledger.meltwater_m3 * lithosphere.RHO_ICE / lithosphere.RHO_WATER)
    assert routed_melt_m3 < ledger.meltwater_m3 + ledger.calved_m3
    assert ice.bodies[0].saturated


def test_snow_in_an_enclosed_pit_below_sea_level_stays_as_ice():
    # #552: a pit below sea level that no seawater reaches is dry land to the ice. Thin ice
    # in it used to count as afloat and calve to an ocean it never crossed the land to reach.
    # It stays, as grounded ice, and the ledger closes without calving.
    points, r, elevation, area = _continent(land_radius_km=900.0)
    pit = r < 250e3
    elevation = np.where(pit, -60.0, elevation)
    is_ocean = (elevation <= 0.0) & ~pit
    # Physical snowfall (RFC 0015 D1, #789) buries this fixture in ~4 km over 10 kyr,
    # past any holdable extent; 1 kyr keeps the pit ice thin, which is all the
    # enclosed-grounds-and-stays mechanism needs.
    depth, _, ice = _step(points, area, elevation, np.zeros(len(points)), _cold(len(points)), 1e3, is_ocean=is_ocean)
    ledger = ice.ledger
    assert ledger.calved_m3 == 0.0
    assert float(np.dot(depth[pit], area[pit])) > 0.0
    assert ledger.closing_m3 == pytest.approx(ledger.accumulation_m3 - ledger.meltwater_m3)
    assert ledger.residual_m3 == pytest.approx(0.0, abs=1e-6 * ledger.accumulation_m3)


def test_a_body_beside_an_enclosed_deep_basin_spreads_over_it_instead_of_calving():
    # #552: the island of the marine case above, but its surrounding 5 km deep water is an
    # enclosed basin the ocean doesn't reach. The ice grounds on that floor and spreads over
    # it rather than saturating and calving, and the ledger closes.
    points, r, elevation, area = _continent(land_radius_km=200.0, half_km=600.0, ocean_floor=-5000.0)
    # Physical ice (RFC 0015 D1, #789): 4000 m plus 1 Myr of snowfall saturates any
    # holdable extent; 500 m over 1 kyr keeps the same spread-over-the-enclosed-
    # basin mechanism within bounds.
    start = np.where(r <= 200e3, 500.0, 0.0)
    depth, _, ice = _step(points, area, elevation, start, _cold(len(points)), 1e3, is_ocean=np.zeros(len(points), dtype=bool))
    ledger = ice.ledger
    assert ledger.calved_m3 == 0.0
    assert not ice.bodies[0].saturated
    assert np.count_nonzero(depth[r > 200e3]) > 0
    assert ledger.residual_m3 == pytest.approx(0.0, abs=1e-6 * ledger.opening_m3)


def test_only_connected_ocean_water_floats_ice():
    bed = np.array([-100.0, -100.0, 100.0])
    ocean = np.array([True, False, False])
    thin = np.array([20.0, 20.0, 20.0])
    assert list(ice_sheet.is_grounded(thin, bed, 0.0, ocean)) == [False, True, True]
    assert list(ice_sheet.flotation_thickness_m(bed, 0.0, ocean) > 0.0) == [True, False, False]
    # lithosphere loads exactly what is_grounded keeps
    assert list(lithosphere.grounded_ice_load_kg_m2(thin, bed, 0.0, ocean) > 0.0) == [False, True, True]


def test_uniform_slab_velocity_matches_balance_velocity():
    # A chain of cells sqrt(area) apart under a uniform slab, gaining b everywhere and keeping
    # none of it: the balance flux at x is b x sqrt(area), so u = b x / H.
    n, spacing = 12, 50e3
    theta = spacing / _R_M * np.arange(n)
    points = np.stack([np.cos(theta), np.sin(theta), np.zeros(n)], axis=1)
    neighbors = np.stack([np.clip(np.arange(n) - 1, 0, n - 1), np.clip(np.arange(n) + 1, 0, n - 1)], axis=1)
    neighbors[0, 0], neighbors[-1, 1] = 1, n - 2
    area = np.full(n, spacing**2)
    thickness = np.full(n, 500.0)
    surface = 2000.0 - 10.0 * np.arange(n) + thickness
    climate = _cold(n, precipitation=400.0)
    # Physical snowfall by mass (RFC 0015 D1, #789): all 400 mm/yr falls as snow at -30 C.
    b_m_per_yr = 400.0 * 1e-3 * lithosphere.RHO_WATER / lithosphere.RHO_ICE
    for warm in (True, False):
        _, velocity, sliding = ice_sheet._diagnose_flow(
            points, neighbors, area, surface, thickness, thickness, thickness, np.ones(n), climate, np.full(n, warm), 1e6,
        )
        x = spacing * (np.arange(n) + 1)
        assert velocity == pytest.approx(b_m_per_yr * x / thickness, rel=1e-6)
        if warm:
            assert sliding == pytest.approx(ice_sheet.ICE_SLIDING_FRACTION * velocity)
        else:
            assert np.all(sliding == 0.0)  # frozen to its bed: no sliding, so no abrasion


def test_bed_thaws_under_thick_ice_and_yield_stress_softens_only_on_thawed_sediment():
    assert not ice_sheet.warm_based(np.array([600.0]), np.array([-20.0]))[0]
    assert ice_sheet.warm_based(np.array([800.0]), np.array([-20.0]))[0]
    assert ice_sheet.warm_based(np.array([1.0]), np.array([0.5]))[0]
    tau = ice_sheet.yield_stress_pa(np.array([0.0, 50.0, 50.0]), np.array([True, False, True]))
    assert tau[0] == tau[1] == ice_sheet.ICE_YIELD_STRESS_HARD_PA
    assert tau[2] == ice_sheet.ICE_YIELD_STRESS_SOFT_PA


def test_seasonal_melt_switches_on_around_the_freezing_threshold():
    # Bare-ice melt (RFC 0015 D1, #789): with no snow on the ground, degree-days above 0 C
    # ablate bare ice wherever the lapsed surface thaws. Read at 1 km of ice so the surface
    # sits 6.5 C below the given annual mean.
    def surface_melt(temperature, amplitude):
        n = len(temperature)
        clim = ice_sheet.SurfaceClimate(
            np.zeros(n), np.asarray(temperature, float), np.asarray(amplitude, float), np.ones(n, bool)
        )
        return clim.rates(np.full(n, 1000.0))[1]

    # With no swing, surface melt switches on exactly at a 0C surface mean.
    no_swing = surface_melt([6.4, 6.5, 6.6], [0.0, 0.0, 0.0])
    assert no_swing[0] == 0.0 and no_swing[1] == 0.0 and no_swing[2] > 0.0
    # A valley whose annual mean is below freezing but whose summer peaks above it melts; the
    # same mean with a swing too small to reach 0C doesn't.
    melt = surface_melt([1.5, 1.5, -8.5], [15.0, 4.0, 10.0])
    assert melt[0] > 0.0 and melt[1] == 0.0 and melt[2] == 0.0
    assert np.all(np.diff(surface_melt(np.full(4, -1.5), [10.0, 15.0, 20.0, 25.0])) > 0.0)
    # A thicker sheet's surface is colder (lapse rate), so it melts less. Snow on the
    # ground supplies this one: bare ice with no snow has nothing to melt.
    clim = ice_sheet.SurfaceClimate(np.full(2, 1000.0), np.full(2, -3.0), np.full(2, 10.0), np.ones(2, bool))
    melt = clim.rates(np.array([0.0, 1500.0]))[1]
    assert melt[1] < melt[0]


def test_seasonal_melt_prevents_permanent_ice_from_annual_mean_alone():
    # Snowy ground at a -8C annual mean: with a strong seasonal swing the summer melt outpaces
    # the snow and no ice survives; with no swing the same ground builds a glacier.
    points, r, elevation, area = _continent(land_radius_km=500.0, half_km=700.0)

    def run(amplitude):
        climate = _cold(len(points), 500.0, -8.0, amplitude)
        depth = np.zeros(len(points))
        for _ in range(5):
            depth, _, _ = _step(points, area, elevation, depth, climate, 1e6)
        return float(np.dot(depth, area))

    assert run(15.0) == 0.0
    assert run(0.0) > 0.0


def test_cold_deep_ocean_calves_legacy_floating_ice_and_accumulates_no_cap():
    points, _, elevation, area = _continent(land_radius_km=300.0, half_km=600.0, ocean_floor=-6000.0)
    ocean = elevation < 0.0
    opening = np.where(ocean, 20.0, 0.0)
    depth, _, ice = _step(points, area, elevation, opening, _cold(len(points), 100.0), 1e3)

    assert np.all(depth[ocean] == 0.0)
    assert ice.ledger.calved_m3 == pytest.approx(float(np.dot(opening, area)))
    assert ice.frozen_precipitation.ocean_retained_m3 == 0.0
    assert abs(ice.ledger.residual_m3) < 1e-9 * ice.ledger.closing_m3


def _ice_ledger(opening=0.0, surface=0.0, basal=0.0, calved=0.0, overflow=0.0, capped=0.0):
    return ice_sheet.IceLedger(
        opening_m3=opening, surface_melt_m3=surface, basal_melt_m3=basal, calved_m3=calved,
        overflow_m3=overflow, capped_m3=capped,
    )


def _retained_only(volume):
    return ice_sheet.FrozenPrecipitation(volume, 0.0, volume, 0.0)


def test_glacier_water_account_accumulates_without_any_ice_loss():
    step = ice_sheet.update_glacier_water_account(100.0, _ice_ledger(opening=5.0), _retained_only(40.0))
    assert step.closing_m3 == 140.0
    assert step.release_m3 == 0.0
    assert not step.initialized_from_ice


def test_glacier_water_account_releases_the_fraction_of_opening_ice_lost_by_destination():
    # A quarter of the opening ice goes, so a quarter of the balance does, split by what took it.
    ledger = _ice_ledger(opening=8.0, surface=0.5, basal=0.5, calved=0.5, overflow=0.25, capped=0.25)
    step = ice_sheet.update_glacier_water_account(100.0, ledger, _retained_only(20.0))
    assert step.release_m3 == pytest.approx(30.0)
    assert step.melt_release_m3 == pytest.approx(15.0)
    assert step.calving_release_m3 == pytest.approx(7.5)
    assert step.overflow_release_m3 == pytest.approx(3.75)
    assert step.guard_release_m3 == pytest.approx(3.75)
    assert step.closing_m3 == pytest.approx(90.0)
    assert abs(step.residual_m3) <= 1e-12 * step.opening_m3


def test_glacier_water_account_calving_alone_leaves_no_meltwater():
    step = ice_sheet.update_glacier_water_account(100.0, _ice_ledger(opening=4.0, calved=1.0), _retained_only(0.0))
    assert step.calving_release_m3 == pytest.approx(25.0)
    assert step.melt_release_m3 == 0.0


def test_glacier_water_account_caps_release_when_losses_exceed_the_opening_ice():
    step = ice_sheet.update_glacier_water_account(100.0, _ice_ledger(opening=1.0, surface=3.0), _retained_only(20.0))
    assert step.release_m3 == pytest.approx(120.0)
    assert step.closing_m3 == pytest.approx(0.0, abs=1e-9)


def test_glacier_water_account_releases_nothing_without_opening_ice():
    # New ice that melts in its first step has no opening fraction to follow (RFC 0016).
    step = ice_sheet.update_glacier_water_account(0.0, _ice_ledger(surface=2.0), _retained_only(20.0))
    assert step.release_m3 == 0.0
    assert step.closing_m3 == 20.0


def test_glacier_water_account_opens_a_save_from_its_opening_ice_by_mass():
    # Physical opening ice converts to water by mass (RFC 0015 D1, #789): opening at
    # WATER_TO_ICE would book ~50,000x too much water (#806 review).
    step = ice_sheet.update_glacier_water_account(None, _ice_ledger(opening=2.0), _retained_only(0.0))
    assert step.initialized_from_ice
    assert step.opening_m3 == pytest.approx(2.0 * lithosphere.RHO_ICE / lithosphere.RHO_WATER)


def test_steady_ice_returns_retained_snow_at_the_opening_conversion():
    # Physical loss balancing physical accumulation holds the account level: the ice and
    # water budgets differ by exactly the by-mass conversion the opening balance uses.
    water = 3.0e6
    ice = water * lithosphere.RHO_WATER / lithosphere.RHO_ICE
    balance = 1.0e9
    ledger = _ice_ledger(opening=balance * lithosphere.RHO_WATER / lithosphere.RHO_ICE, surface=ice)
    step = ice_sheet.update_glacier_water_account(balance, ledger, _retained_only(water))
    assert step.closing_m3 == pytest.approx(balance, rel=1e-3)


def test_frozen_precipitation_follows_the_partition_snow_shares():
    # The account books exactly the partition's dry-share snow, so rain plus frozen land
    # sums to the dry precipitation -- not the binary freeze mask, which double-counts on
    # seasonal cells and misses most snow near freezing (#806 review).
    points, r, elevation, area = _continent()
    n = len(points)
    cold = r < 500_000.0
    snow = np.where(cold, 500.0, 250.0)
    climate = ice_sheet.SurfaceClimate(
        np.full(n, 500.0), np.where(cold, -30.0, -2.0), np.where(cold, 0.0, 15.0), np.ones(n, dtype=bool),
        None, snow,
    )
    years = 100_000.0
    *_, ice = _step(points, area, elevation, np.zeros(n), climate, years)
    frozen = ice.frozen_precipitation
    land = elevation > 0.0
    assert frozen.land_m3 == pytest.approx(float(np.dot(snow * land, area)) / 1000.0 * years)
    assert frozen.retained_m3 > 0.0
    assert frozen.bare_melted_m3 > 0.0
    assert frozen.land_m3 + frozen.ocean_retained_m3 == pytest.approx(frozen.retained_m3 + frozen.bare_melted_m3)


def test_frozen_precipitation_falls_back_to_the_seasonal_share():
    # Direct fixtures predate the partition shares: frozen falls back to the ground-phase
    # seasonal share instead of the binary mask.
    points, r, elevation, area = _continent()
    n = len(points)
    temperature = np.where(r < 500_000.0, -30.0, -2.0)
    amplitude = np.where(r < 500_000.0, 0.0, 15.0)
    climate = ice_sheet.SurfaceClimate(
        np.full(n, 500.0), temperature, amplitude, np.ones(n, dtype=bool),
    )
    years = 100_000.0
    *_, ice = _step(points, area, elevation, np.zeros(n), climate, years)
    frozen = ice.frozen_precipitation
    land = elevation > 0.0
    expected = ice_sheet.snow_fraction(temperature, amplitude) * 500.0
    assert frozen.land_m3 == pytest.approx(float(np.dot(expected * land, area)) / 1000.0 * years)


def _assert_retained_snow_is_the_scaled_accumulation(ice):
    # Everywhere frozen, the physical water the ice keeps and its accumulation are the same
    # snow, by mass (RFC 0015 D1, #789): ice volume times RHO_ICE / RHO_WATER is water.
    assert ice.ledger.accumulation_m3 > 0.0
    assert ice.frozen_precipitation.retained_m3 == pytest.approx(
        ice.ledger.accumulation_m3 * lithosphere.RHO_ICE / lithosphere.RHO_WATER
    )


def test_plastic_body_keeps_snow_on_the_bare_ground_that_feeds_it():
    points, r, elevation, area = _continent()
    n = len(points)
    # Bare margin cells are cold enough to feed the body without being covered by it. The
    # step stays short (physical ice fills its whole domain within 2 kyr, leaving no bare
    # feeders, #789).
    climate = _cold(n, 500.0)
    _, _, ice = _step(points, area, elevation, np.zeros(n), climate, 500.0)
    assert ice.bodies and np.any((ice.covered < 1.0) & (elevation > 0.0))
    _assert_retained_snow_is_the_scaled_accumulation(ice)


def test_uniform_body_keeps_its_snow_and_the_ocean_around_it_keeps_nothing():
    # An island too small for a plastic profile, in water too deep for ice to ground.
    points, r, elevation, area = _continent(land_radius_km=70.0, half_km=600.0, ocean_floor=-6000.0)
    n = len(points)
    assert 0 < np.count_nonzero(elevation > 0.0) < ice_sheet.ICE_BODY_MIN_NODES
    _, _, ice = _step(points, area, elevation, np.zeros(n), _cold(n, 500.0), 100_000.0)
    _assert_retained_snow_is_the_scaled_accumulation(ice)
    assert ice.frozen_precipitation.ocean_retained_m3 == 0.0


def _still(n, precipitation=0.0, temperature=-30.0, amplitude=0.0, frictional=0.0):
    return ice_sheet.SurfaceClimate(
        np.full(n, precipitation), np.full(n, temperature), np.full(n, amplitude),
        np.full(n, temperature < 0.0), np.full(n, frictional),
    )


def test_d1_units_table_rows():
    # RFC 0015 D1's units table (#789), row by row. The DDF values and rain-snow
    # threshold are G2's (#782): DDF_snow 4, DDF_ice 8 mm water/day/K, T_rs 1.0 C.
    to_ice = lithosphere.RHO_WATER / lithosphere.RHO_ICE
    # Precipitation to ice by mass: 1000 mm/yr of water is 1000/917 m/yr of ice.
    acc, _, _ = _still(3, precipitation=1000.0).rates(np.zeros(3))
    assert acc == pytest.approx(np.full(3, 1000.0 * 1e-3 * to_ice * 1e6))
    # Surface melt ablates degree-days: 3.5 C over the whole year at DDF_ice 8.
    _, melt, _ = _still(3, temperature=10.0).rates(np.full(3, 1000.0))
    assert melt == pytest.approx(np.full(3, 8e-3 * (3.5 * 365.25) * to_ice * 1e6))
    # Degree-day melt spends the year's snow first: at a lapsed surface of -0.5 C with a
    # 1 C swing every drop is snow (the seasonal share clips: (1 - -0.5) / 1 > 1), and
    # bare ice only sees the degree-days left after the snow melts.
    snow_climate = ice_sheet.SurfaceClimate(
        np.full(2, 100.0), np.full(2, 6.0), np.ones(2), np.zeros(2, dtype=bool)
    )
    _, melt, _ = snow_climate.rates(np.full(2, 1000.0))
    snow_we = 100.0 * 1e-3
    pdd = 365.25 * hydrology.seasonal_positive_degrees(np.array([-0.5]), np.array([1.0]), 0.0)[0]
    expected = (snow_we + 8e-3 * (pdd - snow_we / 4e-3)) * to_ice * 1e6
    assert expected > snow_we * to_ice * 1e6  # the bare ice really melts too
    assert melt == pytest.approx(np.full(2, expected))
    # Heat flux to basal melt: Y / (rho_i L_f) m/yr per W/m^2.
    w_per_m_yr = ice_sheet.SECONDS_PER_YEAR / (lithosphere.RHO_ICE * ice_sheet.ICE_LATENT_HEAT_J_KG)
    assert ice_sheet.W_M2_TO_M_ICE_PER_YR == pytest.approx(w_per_m_yr)
    # Frictional heat from Pa times m/yr: 100 kPa sliding at 10 m/yr.
    assert ice_sheet.frictional_heat_w_m2(np.array([100_000.0]), np.array([10.0])) == pytest.approx(
        100_000.0 * 10.0 / ice_sheet.SECONDS_PER_YEAR
    )


def test_d1_snow_fraction_bounds_continuity_and_still_air_fallback():
    # The snow share is a share: within [0, 1]; all snow where even the summer stays
    # below the threshold, all rain where even the winter stays above it.
    cold = ice_sheet.snow_fraction(np.full(4, -30.0), np.array([0.0, 5.0, 15.0, 29.0]))
    assert np.all(cold == 1.0)
    warm = ice_sheet.snow_fraction(np.full(4, 30.0), np.array([0.0, 5.0, 15.0, 28.0]))
    assert np.all(warm == 0.0)
    # Continuous in temperature: no jumps on a fine grid, including across the seasonal
    # formula's clip edges and the still-air ramp.
    grid = np.linspace(-30.0, 30.0, 1201)
    for amplitude in (0.0, 1e-4, 5.0, 15.0):
        share = ice_sheet.snow_fraction(grid, np.full_like(grid, amplitude))
        assert np.all((share >= 0.0) & (share <= 1.0))
        # The seasonal formula has a vertical tangent at its clip edges, so a fine grid
        # still steps ~0.05 there; a true jump would step ~1.
        assert np.max(np.abs(np.diff(share))) < 0.1
    # A 0.01 C amplitude would switch over 0.02 C unclamped -- a near-step the fine grid
    # would cross in one ~1.0 jump. The amplitude clamps to half the ramp width, so the
    # switch spans the full ramp degree (#806 review): still steeper at the clip-edge
    # tangents (~0.14) than the wide-amplitude case, but continuous.
    clamped = ice_sheet.snow_fraction(grid, np.full_like(grid, 0.01))
    assert np.max(np.abs(np.diff(clamped))) < 0.2
    tiny = ice_sheet.snow_fraction(np.array([0.5, 1.0, 1.5]), np.full(3, 0.01))
    assert tiny == pytest.approx([1.0, 0.5, 0.0])
    # Where the seasonal amplitude vanishes the sinusoid has no year to spread over, so
    # the share falls back to the continuous ramp around the threshold.
    ramp = np.clip(0.5 + (ice_sheet.RAIN_SNOW_THRESHOLD_C - grid) / ice_sheet.ICE_SNOW_RAMP_C, 0.0, 1.0)
    assert ice_sheet.snow_fraction(grid, np.full_like(grid, 1e-9)) == pytest.approx(ramp)
    assert ice_sheet.snow_fraction(grid, np.zeros_like(grid)) == pytest.approx(ramp)


def test_trial_surface_repartitions_dry_precipitation_but_keeps_bed_frozen_snow():
    climate = ice_sheet.SurfaceClimate(
        precipitation_mm=np.array([800.0]),
        temperature_c=np.array([6.0]),
        seasonal_amplitude_c=np.array([10.0]),
        is_frozen=np.array([True]),
        dry_snow_mm=np.array([800.0 / 3.0]),
    )
    without_lake_snow = climate.rates(np.array([5000.0]))[0][0]
    assert without_lake_snow > climate.rates(np.array([0.0]))[0][0]

    climate.bed_frozen_snow_mm = np.array([200.0])
    climate = climate.subset(np.array([0]))
    to_snow_mm = ice_sheet.RHO_ICE / ice_sheet.RHO_WATER / 1000.0

    bare_snow_mm = climate.rates(np.array([0.0]))[0][0] * to_snow_mm
    thick_snow_mm = climate.rates(np.array([5000.0]))[0][0] * to_snow_mm

    assert bare_snow_mm == pytest.approx(200.0 + 800.0 / 3.0)
    assert thick_snow_mm == pytest.approx(1000.0)


def test_d1_basal_melt_continuous_at_the_thaw_boundary_and_bounded_as_h_to_zero():
    # At a fixed cold surface the bed thaws past a threshold thickness: conducted heat
    # G = k (-p H - T_s) / H solves to H* = -k T_s / (G + k p). The melt is zero exactly
    # there and small just past it (RFC 0015 D1, #789).
    surface = -20.0
    k = ice_sheet.ICE_THERMAL_CONDUCTIVITY_W_M_K
    p = ice_sheet.ICE_PRESSURE_MELT_C_PER_M
    g = ice_sheet.ICE_GEOTHERMAL_FLUX_W_M2
    h_star = -k * surface / (g + k * p)
    heats = ice_sheet.basal_heat_w_m2(np.array([h_star - 1.0, h_star + 1.0]), np.full(2, surface), 0.0)
    assert heats[0] < 0.0 < heats[1]
    melts = ice_sheet.W_M2_TO_M_ICE_PER_YR * np.clip(heats, 0.0, None) * 1e6
    assert melts[0] == 0.0
    deep = ice_sheet.W_M2_TO_M_ICE_PER_YR * (g - k * (-p * 5000.0 - surface) / 5000.0) * 1e6
    assert 0.0 < melts[1] < 0.05 * deep
    # No ice is no melt, and a vanishing sheet under a warm surface stays bounded: the
    # conducted term clips at zero instead of blowing up as 1/H (RFC 0015 D1, #789).
    assert ice_sheet.basal_heat_w_m2(np.array([0.0]), np.array([5.0]), 0.0)[0] < 0.0
    thin = ice_sheet.W_M2_TO_M_ICE_PER_YR * ice_sheet.basal_heat_w_m2(
        np.array([1e-6]), np.array([5.0]), 0.0
    )[0] * 1e6
    assert 0.0 < thin < 1e5


def test_d1_frictional_melt_reads_last_step_heat_inside_the_root():
    # A bed frozen under geothermal heat alone thaws with last step's frictional heat:
    # at 500 m under -20 C conduction (0.082 W/m^2) beats geothermal (0.06), and 100 kPa
    # sliding at 25 m/yr adds 0.079. The trial-thickness melt inside rates() then rises
    # with diminishing steps (RFC 0015 D1, #789).
    surface, thick = -20.0, 500.0
    cold = _still(1, frictional=0.0)
    _, _, still_basal = cold.rates(np.array([thick]))
    assert still_basal[0] == 0.0
    sliding = ice_sheet.frictional_heat_w_m2(
        ice_sheet.yield_stress_pa(np.array([0.0]), np.array([True])), np.array([25.0])
    )
    assert sliding[0] > 0.0
    warm = _still(1, temperature=surface, frictional=sliding[0])
    _, _, basal = warm.rates(np.array([thick]))
    assert basal[0] > 0.0
    assert not ice_sheet.warm_based(np.array([thick]), np.array([surface]))[0]
    assert bool(ice_sheet.warm_based(np.array([thick]), np.array([surface]), sliding)[0])
    # Thicker ice conducts less, so basal melt rises toward its geothermal-plus-friction
    # ceiling with diminishing steps (the conducted term decays as 1/H).
    near = np.array([thick, thick + 50.0, thick + 100.0])
    melts = warm.rates(near)[2]
    assert np.all(np.diff(melts) >= 0.0)
    assert melts[1] - melts[0] >= melts[2] - melts[1] - 1e-9


def test_march_cost_default_is_the_trend_1k_march():
    # #842: the production default is the trend march with a 1 kyr longest substep. The
    # fixed-250 yr "d2" march stays selectable, and is still the bare MarchCost().
    assert ice_sheet.MARCH_COST is ice_sheet.march_cost("trend-1k")
    assert ice_sheet.MARCH_COST.trend and ice_sheet.MARCH_COST.trend_max_substep_years == 1_000.0
    d2 = ice_sheet.MarchCost()
    assert ice_sheet.march_cost("d2") == d2
    assert d2.max_substeps is None and not d2.adaptive
    assert d2.sea_level_skip_m == 0.0 and not d2.local_dijkstra and d2.dijkstra_limit_m is None
    with pytest.raises(KeyError, match="cap-hold"):
        ice_sheet.march_cost("no-such-preset")
    with pytest.raises(ValueError):
        ice_sheet.MarchCost(cap_fill="extrapolate")
    with pytest.raises(ValueError):
        ice_sheet.MarchCost(max_substeps=1)


def _patchy_positive(points, seed=822):
    """A positive set of a few disjoint discs on a lattice, for the distance searches."""
    rng = np.random.default_rng(seed)
    centres = points[rng.choice(len(points), size=4, replace=False)]
    return (points @ centres.T > np.cos(250e3 / _R_M)).any(axis=1)


def test_restrict_matches_the_coo_round_trip_exactly():
    # #147: filtering the CSR rows in place must give the matrix the old COO round trip
    # built, stored order included (the distance searches break ties by it).
    from scipy.sparse import csr_matrix

    points, _, _, _ = _lattice(50.0, 1000.0, kind="hex")
    graph = ice_sheet.symmetric_graph(points, _neighbors(points))
    rng = np.random.default_rng(147)
    # The same matrix with each row's edges stored in reverse.
    flip = np.concatenate([np.arange(b - 1, a - 1, -1) for a, b in zip(graph.indptr[:-1], graph.indptr[1:])])
    unsorted = csr_matrix((graph.data[flip], graph.indices[flip], graph.indptr), shape=graph.shape)
    assert not unsorted.has_canonical_format
    for source in (graph, unsorted):
        for share in (0.0, 0.4, 1.0):
            keep = rng.random(len(points)) < share
            coo = source.tocoo()
            ok = keep[coo.row] & keep[coo.col]
            expected = csr_matrix((coo.data[ok], (coo.row[ok], coo.col[ok])), shape=source.shape)
            got = ice_sheet._restrict(source, keep)
            np.testing.assert_array_equal(got.indptr, expected.indptr)
            np.testing.assert_array_equal(got.indices, expected.indices)
            np.testing.assert_array_equal(got.data, expected.data)


def test_local_inside_distance_matches_the_full_search_exactly():
    # #822 option 3: the inside distance searched on S plus its boundary is the full
    # search exactly (a shortest path into S can be cut at its last outside node), so
    # phi and the nearest-source labels are bit-identical.
    points, _, _, _ = _lattice(50.0, 1000.0, kind="hex")
    graph = ice_sheet.symmetric_graph(points, _neighbors(points))
    positive = _patchy_positive(points)
    # Since #147 every search is the local one, so the full search is spelled out here.
    from scipy.sparse.csgraph import dijkstra

    eligible = ice_sheet._restrict(graph, np.ones(len(points), dtype=bool))
    phi_local, _ = ice_sheet._signed_distance(graph, eligible, positive)
    d_in = dijkstra(graph, indices=np.flatnonzero(~positive), min_only=True)
    d_out = dijkstra(eligible, indices=np.flatnonzero(positive), min_only=True)
    phi_full = 0.5 * (np.where(positive, 0.0, d_out) - np.where(positive, d_in, 0.0))
    assert 0 < positive.sum() < len(points)
    np.testing.assert_array_equal(phi_local, phi_full)
    everywhere = np.ones(len(points), dtype=bool)
    assert np.all(ice_sheet._signed_distance(graph, eligible, everywhere)[0] == -np.inf)


def test_limited_outside_distance_falls_back_when_it_cuts_off_a_member():
    # #822 option 3: the limited outside search is exact within the limit, and a member
    # past half the limit makes it rerun unlimited, so no body loses a node.
    points, _, _, _ = _lattice(50.0, 1000.0, kind="hex")
    graph = ice_sheet.symmetric_graph(points, _neighbors(points))
    positive = _patchy_positive(points)
    eligible = ice_sheet._restrict(graph, np.ones(len(points), dtype=bool))
    phi_full, src_full = ice_sheet._signed_distance(graph, eligible, positive)
    limited = ice_sheet.MarchCost(dijkstra_limit_m=400e3)
    near = positive.copy()
    phi_near, _ = ice_sheet._signed_distance(graph, eligible, positive, cost=limited, members=near)
    inside = np.isfinite(phi_near)
    assert np.any(~inside)
    np.testing.assert_array_equal(phi_near[inside], phi_full[inside])
    far = positive | (phi_full > 300e3)
    phi_far, src_far = ice_sheet._signed_distance(graph, eligible, positive, cost=limited, members=far)
    np.testing.assert_array_equal(phi_far, phi_full)
    np.testing.assert_array_equal(src_far, src_full)


class _CountedProfiles:
    """A stand-in for `_solve_plastic_body`'s cached shape builder over a nondecreasing
    V(d), counting the distinct extents built."""

    def __init__(self, volume_of):
        self.volume_of = volume_of
        self.built: set[float] = set()

    def __call__(self, d):
        self.built.add(d)
        return type("Profile", (), {"volume": self.volume_of(d)})()


def test_volume_extent_hint_changes_the_work_not_the_answer():
    # #147: with any hint (good, far off, outside the range, not finite) the bisection
    # returns the very same float as without one; a hint near the answer builds two or
    # three profiles where the full bisection builds over a dozen.
    def volume_of(d):
        return max(d + 300e3, 0.0) ** 2.5  # empty below -300 km, growing like a dome

    d_min, d_max, tolerance = -400e3 - 1.0, 2_750e3, 62e3
    rng = np.random.default_rng(147)
    for answer in rng.uniform(-290e3, 2_700e3, 40):
        volume = volume_of(answer)
        plain = _CountedProfiles(volume_of)
        expected = ice_sheet._volume_extent(plain, volume, d_min, d_max, tolerance)
        assert expected - tolerance < answer <= expected
        for hint in (answer, answer + 100.0, answer - 100.0, answer + 5 * tolerance, answer - 5 * tolerance,
                     d_min - 1e6, d_max + 1e6, d_min, d_max, np.nan, np.inf):
            hinted = _CountedProfiles(volume_of)
            assert ice_sheet._volume_extent(hinted, volume, d_min, d_max, tolerance, hint) == expected
            assert expected in hinted.built or expected == d_max  # the walk from d0 reuses it
        near = _CountedProfiles(volume_of)
        ice_sheet._volume_extent(near, volume, d_min, d_max, tolerance, answer + 100.0)
        assert len(near.built) <= 3 < len(plain.built)
    # Volumes no extent holds, and none at all, are the same with a hint.
    for volume in (volume_of(d_max) * 2.0, 0.0):
        expected = ice_sheet._volume_extent(_CountedProfiles(volume_of), volume, d_min, d_max, tolerance)
        assert ice_sheet._volume_extent(_CountedProfiles(volume_of), volume, d_min, d_max, tolerance, 0.0) == expected
    # A domain narrower than the tolerance, and an empty one.
    for lo, hi in ((0.0, 10e3), (5.0, 5.0)):
        expected = ice_sheet._volume_extent(_CountedProfiles(volume_of), 1.0, lo, hi, tolerance)
        assert ice_sheet._volume_extent(_CountedProfiles(volume_of), 1.0, lo, hi, tolerance, 3.0) == expected


def test_margin_hints_leave_a_marched_ice_step_bit_identical(monkeypatch):
    # #147: stepping a growing dome with the previous step's margins as hints gives the
    # same ice as stepping it without, with fewer profiles built.
    points, x, y, cell = _lattice(125.0, 1500.0, "hex")
    n = len(points)
    neighbors = _neighbors(points)
    area = np.full(n, cell)
    r = np.hypot(x, y)
    climate = _cold(n, precipitation=400.0, temperature=np.where(r < 500e3, -25.0, 8.0), amplitude=10.0)
    no = np.zeros(n, dtype=bool)
    built = []
    real_shape = ice_sheet._shape
    monkeypatch.setattr(ice_sheet, "_shape", lambda *a: built.append(1) or real_shape(*a))

    def march(hinted):
        depth, bodies, history = np.zeros(n), None, []
        for _ in range(6):
            hints = ice_sheet.margin_hints(n, bodies) if hinted and bodies is not None else None
            depth, melt, ice = ice_sheet.step_ice(
                points, neighbors, area, np.zeros(n), np.zeros(n), -1e4, depth, climate, np.zeros(n), no, 2_000.0,
                margin_hint_m=hints,
            )
            bodies = ice.bodies
            history.append((depth, melt, ice.thickness, [b.margin_m for b in bodies]))
        return history

    plain = march(False)
    plain_built = len(built)
    built.clear()
    hinted = march(True)
    assert plain[-1][0].max() > 100.0
    for (d0, m0, t0, g0), (d1, m1, t1, g1) in zip(plain, hinted):
        np.testing.assert_array_equal(d1, d0)
        np.testing.assert_array_equal(m1, m0)
        np.testing.assert_array_equal(t1, t0)
        assert g1 == g0
    assert len(built) < plain_built

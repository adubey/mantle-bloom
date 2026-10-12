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


def _step(points, area, elevation, depth, climate, years, sea_level=0.0, frozen_lake=None, allow_ocean_ice=None):
    n = len(points)
    is_ocean = elevation <= sea_level
    return ice_sheet.step_ice(
        points, _neighbors(points), area, elevation, np.zeros(n), sea_level, depth,
        np.zeros(n) if frozen_lake is None else frozen_lake, climate, np.zeros(n), is_ocean,
        np.zeros(n, dtype=bool) if allow_ocean_ice is None else allow_ocean_ice, years,
    )


def _graded_climate(r, precipitation=500.0, inner_km=500.0):
    """Cold inside `inner_km`, warming outward past freezing, with a 10 C seasonal swing."""
    temperature = -25.0 + 30.0 * np.clip(r / (inner_km * 1e3), 0.0, 2.0) - 15.0
    return _cold(len(r), precipitation, temperature, 10.0)


def test_volume_contract_closes_per_body_and_across_the_world():
    points, r, elevation, area = _continent()
    climate = _graded_climate(r)
    depth = np.zeros(len(points))
    lake = np.zeros(len(points))
    lake[np.argmin(r)] = 40.0  # a lake freezing at the centre joins the ice
    for step in range(3):
        opening = float(np.dot(depth, area))
        depth, melt, ice = _step(points, area, elevation, depth, climate, 1e6, frozen_lake=lake if step == 1 else None)
        ledger = ice.ledger
        scale = max(ledger.accumulation_m3, ledger.closing_m3, 1.0)
        assert ledger.opening_m3 == pytest.approx(opening)
        assert abs(ledger.residual_m3) < 1e-9 * scale
        assert ledger.closing_m3 == pytest.approx(float(np.dot(depth, area)))
        for body in ice.bodies:
            assert abs(body.residual_m3) < 1e-9 * max(body.accumulation_m3, body.closing_m3, 1.0)
        # Ice is never created: closing <= opening + lake ice + accumulation.
        assert ledger.closing_m3 <= ledger.opening_m3 + ledger.frozen_lake_m3 + ledger.accumulation_m3 + 1e-9 * scale
        # Every booked sink reaches the water or the ocean: melt (with snow on bare ground) as
        # meltwater at its node, calving unrouted.
        assert float(np.dot(melt, area)) == pytest.approx(ledger.meltwater_m3, rel=1e-9, abs=1.0)
        assert float(ice.calved_m3.sum()) == pytest.approx(ledger.calved_m3)
    assert ledger.closing_m3 > 0.0


def _steady_state(points, area, elevation, climate, years=1e10, start=None, steps=4):
    depth = np.zeros(len(points)) if start is None else start
    for _ in range(steps):
        depth, _, ice = _step(points, area, elevation, depth, climate, years)
    return depth, ice


def test_more_snowfall_on_the_same_positive_balance_area_builds_a_bigger_sheet():
    # Snow-covered and melt-free inside 500 km, strong summer melt outside it: the positive-
    # balance area is the same disc whatever the snowfall. Twice the snow must still give a
    # larger steady sheet that reaches further out (accumulation magnitude, not just its sign).
    points, r, elevation, area = _continent(land_radius_km=1100.0, half_km=1500.0)
    inside = r <= 500e3
    temperature = np.where(inside, -40.0, 8.0)

    def run(precipitation):
        climate = ice_sheet.SurfaceClimate(np.full(len(r), precipitation), temperature, np.full(len(r), 5.0), inside)
        depth, ice = _steady_state(points, area, elevation, climate)
        return float(np.dot(depth, area)), np.count_nonzero(depth >= hydrology.GLACIER_VISIBLE_DEPTH_M), ice

    v1, n1, ice1 = run(300.0)
    v2, n2, ice2 = run(600.0)
    assert not any(b.saturated for b in ice1.bodies + ice2.bodies)
    assert v2 > 1.1 * v1
    assert n2 > n1


def test_backward_euler_converges_to_zero_balance_for_long_steps_and_to_explicit_for_short_ones():
    points, r, elevation, area = _continent(land_radius_km=1100.0, half_km=1500.0)
    climate = _graded_climate(r, inner_km=600.0)
    depth, ice = _steady_state(points, area, elevation, climate, years=1e10, steps=3)
    ledger = ice.ledger
    # dt * B = closing - opening, so with dt = 10 Gyr the balance is ~zero next to snowfall.
    assert not any(b.saturated for b in ice.bodies)
    assert abs(ledger.accumulation_m3 - ledger.surface_melt_m3 - ledger.basal_melt_m3) < 0.01 * ledger.accumulation_m3

    # Short steps from a growing sheet: dV / dt is the explicit balance at the current extent,
    # the same whichever short step is taken.
    growing, _, _ = _step(points, area, elevation, np.zeros(len(points)), climate, 3e6)
    rates = []
    for years in (1e3, 1e4):
        new, _, ice = _step(points, area, elevation, growing, climate, years)
        rates.append((float(np.dot(new, area)) - float(np.dot(growing, area))) / years)
    assert rates[0] == pytest.approx(rates[1], rel=0.02)


def test_height_mass_balance_feedback_gives_two_stable_sheets():
    # A high plateau whose bare surface melts more than it gets snow (a 1 C summer above
    # freezing), but whose surface under ~1 km of ice is cold enough to keep all its snow.
    # Bare, it stays bare; covered, it keeps a sheet (Weertman 1961).
    points, r, elevation, area = _continent(land_radius_km=1000.0, half_km=1400.0)
    elevation = np.where(elevation > 0.0, 2000.0, elevation)
    n = len(points)
    climate = ice_sheet.SurfaceClimate(np.full(n, 1000.0), np.full(n, -2.0), np.full(n, 8.0), np.ones(n, dtype=bool))
    bare, _ = _steady_state(points, area, elevation, climate, years=2e7, steps=3)
    assert bare.max() == 0.0
    start = np.where(elevation > 0.0, 1500.0, 0.0)
    covered, ice = _steady_state(points, area, elevation, climate, years=2e7, start=start, steps=3)
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
    depth, _, ice = _step(points, area, elevation, start, _cold(len(points)), 1e6)
    loaded = lithosphere.grounded_ice_load_kg_m2(depth, elevation, 0.0) > 0.0
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
    assert ice.bodies[0].saturated


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
    b_m_per_yr = hydrology.GLACIER_ACCUMULATION_RATE * 400.0 / 1e6
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
    def surface_melt(temperature, amplitude):
        n = len(temperature)
        clim = ice_sheet.SurfaceClimate(np.zeros(n), np.asarray(temperature, float), np.asarray(amplitude, float), np.ones(n, bool))
        return clim.rates(np.zeros(n))[1]

    # With no swing, surface melt switches on exactly at a 0C annual mean.
    no_swing = surface_melt([-0.1, 0.0, 0.1], [0.0, 0.0, 0.0])
    assert no_swing[0] == 0.0 and no_swing[1] == 0.0 and no_swing[2] > 0.0
    # A valley whose annual mean is below freezing but whose summer peaks above it melts; the
    # same mean with a swing too small to reach 0C doesn't.
    melt = surface_melt([-5.0, -5.0, -15.0], [15.0, 4.0, 10.0])
    assert melt[0] > 0.0 and melt[1] == 0.0 and melt[2] == 0.0
    assert np.all(np.diff(surface_melt(np.full(4, -8.0), [10.0, 15.0, 20.0, 25.0])) > 0.0)
    # A thicker sheet's surface is colder (lapse rate), so it melts less.
    clim = ice_sheet.SurfaceClimate(np.zeros(2), np.full(2, -3.0), np.full(2, 10.0), np.ones(2, bool))
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


def test_sea_ice_cap_keeps_its_own_rule_and_other_floating_ice_calves():
    points, r, elevation, area = _continent(land_radius_km=300.0, half_km=600.0)
    n = len(points)
    ocean = elevation < 0.0
    start = np.where(ocean, 20.0, 0.0)
    cap = ocean & (np.arange(n) % 2 == 0)
    depth, _, ice = _step(points, area, elevation, start, _cold(n, 100.0), 1e6, allow_ocean_ice=cap)
    gained = hydrology.GLACIER_ACCUMULATION_RATE * 100.0
    assert depth[cap] == pytest.approx(20.0 + gained)
    assert np.all(depth[ocean & ~cap] == 0.0)
    assert ice.ledger.calved_m3 == pytest.approx(float(np.dot(start[ocean & ~cap], area[ocean & ~cap])))


def test_sea_ice_cap_is_held_to_the_depth_guard_booked_and_logged(caplog):
    # A cold spell with no summer melt: the floating cap gains 20 m on 4,990 m. The guard keeps
    # it at GLACIER_MAX_DEPTH_M, sheds the excess as meltwater, books it, and warns.
    points, r, elevation, area = _continent(land_radius_km=300.0, half_km=600.0)
    n = len(points)
    cap = elevation < 0.0
    start = np.where(cap, hydrology.GLACIER_MAX_DEPTH_M - 10.0, 0.0)
    with caplog.at_level("WARNING", logger="app.hydroclimate.ice_sheet"):
        depth, melt, ice = _step(points, area, elevation, start, _cold(n, 1000.0), 1e6, allow_ocean_ice=cap)
    ledger = ice.ledger
    assert np.all(depth[cap] == hydrology.GLACIER_MAX_DEPTH_M)
    assert ledger.capped_m3 == pytest.approx(10.0 * float(area[cap].sum()))
    assert ledger.capped_nodes >= int(cap.sum())
    assert abs(ledger.residual_m3) < 1e-9 * ledger.closing_m3
    assert float(np.dot(melt[cap], area[cap])) == pytest.approx(ledger.capped_m3)
    assert "GLACIER_MAX_DEPTH_M" in caplog.text

    caplog.clear()
    with caplog.at_level("WARNING", logger="app.hydroclimate.ice_sheet"):
        _, _, ice = _step(points, area, elevation, np.where(cap, 100.0, 0.0), _cold(n, 1000.0), 1e6, allow_ocean_ice=cap)
    assert ice.ledger.capped_m3 == 0.0 and ice.ledger.capped_nodes == 0
    assert caplog.text == ""

import numpy as np

from app import erosion, geometry, plates
from app.elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M
from app.world import World, generate_world


def test_climate_grid_indices_matches_build_grid_convention():
    # Row 0 = north pole, row increases southward; column increases eastward from lon=-180.
    world_xyz = np.array(
        [
            geometry.latlon_to_xyz(np.radians(89.0), np.radians(-179.0)),
            geometry.latlon_to_xyz(np.radians(-89.0), np.radians(179.0)),
            geometry.latlon_to_xyz(0.0, 0.0),
        ]
    )
    row, col = erosion.climate_grid_indices(world_xyz, height=90, width=180)
    assert row[0] == 0
    assert row[1] == 89
    assert col[0] == 0
    assert col[2] == 90  # lon=0 falls in the middle column


def test_compute_slope_zero_for_flat_cluster():
    # A tight cluster of same-elevation points -- nobody is lower than anybody else, so
    # both the slope and its cap (drop to lowest neighbor) must be exactly 0.
    theta = np.linspace(0.0, 0.01, 8)
    points = geometry.local_xyz(np.zeros_like(theta), theta)
    elevation = np.full(8, 500.0)
    slope, drop_m = erosion.compute_slope(points, elevation)
    assert np.allclose(slope, 0.0)
    assert np.allclose(drop_m, 0.0)


def test_compute_slope_matches_known_gradient():
    # 6 points along a line with elevation decreasing monotonically by 100m per step --
    # point 0's SLOPE_NEIGHBOR_COUNT=4 nearest neighbors are points 1-4, and since
    # elevation decreases monotonically, the *lowest* of those 4 is point 4 (farthest in
    # this small set, but lowest).
    d = 0.01
    theta = d * np.arange(6)
    points = geometry.local_xyz(np.zeros_like(theta), theta)
    elevation = 500.0 - 100.0 * np.arange(6)

    slope, drop_m = erosion.compute_slope(points, elevation)

    expected_drop = elevation[0] - elevation[4]
    expected_run_m = geometry.angular_distance(points[0], points[4]) * plates.PLANET_RADIUS_KM * 1000.0
    assert np.isclose(drop_m[0], expected_drop)
    assert np.isclose(slope[0], expected_drop / expected_run_m)


def test_weathering_relief_factor_suppresses_flat_terrain_and_saturates_on_steep_terrain():
    # Regression check for the "flat coastal plain weathers exactly as fast as a mountain"
    # bug: weathering must scale down toward 0 as slope -> 0 (same clip-and-normalize idiom
    # humidity_norm already uses), and saturate to full strength on genuinely steep terrain,
    # rather than being a flat function of wind/humidity alone.
    d = 0.01
    theta = d * np.arange(6)
    points = geometry.local_xyz(np.zeros_like(theta), theta)

    flat_elevation = np.full(6, 30.0)  # uniform -- see test_compute_slope_zero_for_flat_cluster
    flat_slope, _ = erosion.compute_slope(points, flat_elevation)
    flat_relief_factor = np.clip(flat_slope / erosion.WEATHERING_RELIEF_REFERENCE_SLOPE, 0.0, 1.0)
    assert np.allclose(flat_relief_factor, 0.0)

    steep_elevation = 30.0 + 2000.0 * np.arange(6)  # steep monotonic gradient
    steep_slope, _ = erosion.compute_slope(points, steep_elevation)
    steep_relief_factor = np.clip(steep_slope / erosion.WEATHERING_RELIEF_REFERENCE_SLOPE, 0.0, 1.0)
    # Every point except the global minimum (point 0, which has no lower neighbor at all)
    # should be steep enough to saturate weathering to full strength.
    assert np.allclose(steep_relief_factor[1:], 1.0)


def test_submarine_erosion_scales_with_slope_and_depth_only_over_ocean():
    # Flat sea floor (slope 0) erodes nothing; a steeper, deeper submerged scarp erodes more
    # than a steeper, shallower one; subaerial nodes are untouched by this term.
    is_ocean = np.array([True, True, True, False])
    elevation = np.array([-100.0, -100.0, -5000.0, 500.0])
    slope = np.array([0.0, 0.03, 0.03, 0.03])
    amt = erosion.submarine_erosion_amount(elevation, slope, is_ocean, dt_myr=5.0)
    assert amt[0] == 0.0  # flat
    assert amt[2] > amt[1] > 0.0  # deeper (more pressure) erodes faster at the same slope
    assert amt[3] == 0.0  # subaerial


def test_coastal_erosion_confined_to_band_and_peaks_near_freezing():
    elevation = np.array([0.0, 0.0, 150.0, 5000.0, -5000.0])
    freezing = np.full(5, erosion.COASTAL_FROST_PEAK_C)
    warm = np.full(5, 30.0)
    at_freezing = erosion.coastal_erosion_amount(elevation, freezing, dt_myr=1.0)
    when_warm = erosion.coastal_erosion_amount(elevation, warm, dt_myr=1.0)
    assert at_freezing[0] > when_warm[0] > 0.0  # frost adds on top of wave attack at the shore
    assert at_freezing[2] > 0.0 and at_freezing[2] < at_freezing[0]  # in-band but tapering
    assert at_freezing[3] == 0.0 and at_freezing[4] == 0.0  # far above / far below sea level


def test_apply_erosion_can_lower_a_submerged_range():
    # End to end: over a generated world, some ocean node loses net elevation once submarine +
    # coastal erosion are in play (they were entirely erosion-exempt before).
    world = generate_world(seed=20, num_plates=8)
    _, elevation_before, _, _, _, _ = erosion._gather_nodes(world)
    is_ocean = elevation_before <= 0.0
    erosion.apply_erosion(world, years=5_000_000)
    _, elevation_after, _, _, _, _ = erosion._gather_nodes(world)
    assert len(elevation_after) == len(elevation_before)
    assert np.any(elevation_after[is_ocean] < elevation_before[is_ocean] - 1e-6)


def test_apply_erosion_keeps_elevation_finite_and_changing():
    world = generate_world(seed=21, num_plates=8)
    _, elevation_before, _, _, _, _ = erosion._gather_nodes(world)
    erosion.apply_erosion(world, years=5_000_000)
    _, elevation_after, _, _, _, _ = erosion._gather_nodes(world)
    # Erosion/deposition/weathering can now raise *or* lower any given node (deposition can
    # outweigh local erosion at a floodplain or delta) -- just confirm the pass actually did
    # something and didn't produce garbage.
    assert np.all(np.isfinite(elevation_after))
    assert not np.allclose(elevation_after, elevation_before)


def test_apply_erosion_respects_elevation_bounds():
    world = generate_world(seed=22, num_plates=8)
    erosion.apply_erosion(world, years=5_000_000)
    for plate in world.plates:
        for line in plate.lines:
            assert np.all(line.elevation >= MIN_ELEVATION_M - 1e-6)
            assert np.all(line.elevation <= MAX_ELEVATION_M + 1e-6)


def test_apply_erosion_noop_for_empty_world():
    world = World(seed=0, plates=[])
    erosion.apply_erosion(world, years=1_000_000)  # must not raise
    assert world.plates == []


# --- Coastal planation + infill feedback (docs/TODO.md "Speckled low-relief coastlines") ---


def _equator_lattice(half_span_deg: float, step_deg: float) -> np.ndarray:
    """A square lattice of unit vectors straddling (lat, lon) = (0, 0), for the coastal-
    feedback helpers (which only care about relative node geometry, not which plate owns
    what)."""
    grid = np.radians(np.arange(-half_span_deg, half_span_deg + 1e-9, step_deg))
    lat, lon = np.meshgrid(grid, grid, indexing="ij")
    return geometry.latlon_to_xyz(lat.ravel(), lon.ravel())


def test_coastal_openness_tracks_open_water_fraction():
    points = _equator_lattice(half_span_deg=3.0, step_deg=0.4)
    lon_deg = np.degrees(geometry.xyz_to_latlon(points)[1])
    is_ocean = lon_deg > 0.0  # a straight N-S coastline down the prime meridian

    openness = erosion._coastal_openness(points, is_ocean)

    deep_land = np.argmin(lon_deg)  # far to the west, ringed by land
    deep_ocean = np.argmax(lon_deg)  # far to the east, ringed by ocean
    coast = np.argmin(np.abs(lon_deg) + np.abs(np.degrees(geometry.xyz_to_latlon(points)[0])))
    assert openness[deep_land] < 0.05
    assert openness[deep_ocean] > 0.95
    assert 0.3 < openness[coast] < 0.7


def test_coastal_openness_zero_when_no_open_ocean():
    points = _equator_lattice(half_span_deg=2.0, step_deg=0.4)
    assert np.all(erosion._coastal_openness(points, np.zeros(len(points), dtype=bool)) == 0.0)


def test_coastal_planation_only_near_sea_level_land_and_scales_with_exposure():
    elevation = np.array([10.0, 25.0, 100.0, -50.0])
    openness = np.full(4, 0.5)  # well above PLANATION_EXPOSURE_REF -> full exposure
    amt = erosion.coastal_planation_amount(elevation, sea_level_m=0.0, openness=openness, dt_myr=0.05)
    assert amt[0] > amt[1] > 0.0  # both in-band, shallower planes faster (bigger proximity)
    assert amt[2] == 0.0  # above the band
    assert amt[3] == 0.0  # underwater -- infill's job, not planation's

    sheltered = erosion.coastal_planation_amount(elevation, 0.0, np.full(4, 0.15), dt_myr=0.05)
    assert sheltered[0] < amt[0]  # half the exposure -> half the planation (0.15 / 0.3)

    # Never past the wave-cut platform (sea level minus PLANATION_UNDERCUT_M at full
    # exposure), however long the step -- here 12 m above + 6 m undercut = 18 m.
    huge = erosion.coastal_planation_amount(np.array([12.0]), 0.0, np.array([1.0]), dt_myr=99.0)
    assert huge[0] == 12.0 + erosion.PLANATION_UNDERCUT_M


def test_spread_coastal_infill_conserves_mass_including_fallback():
    rng = np.random.default_rng(0)
    points = _equator_lattice(half_span_deg=3.0, step_deg=0.5)
    n = len(points)
    elevation = rng.uniform(-500.0, 200.0, n)
    is_ocean = elevation <= 0.0
    openness = rng.uniform(0.0, 1.0, n)
    dist_to_land = rng.uniform(0.0, 0.02, n)
    source = np.zeros(n)
    source[rng.choice(n, 5, replace=False)] = rng.uniform(1.0, 10.0, 5)

    out = erosion._spread_coastal_infill(points, elevation, is_ocean, openness, dist_to_land, 0.0, source)
    assert np.isclose(out.sum(), source.sum())

    # No ocean at all -> every source keeps its own amount in place.
    land_only = erosion._spread_coastal_infill(
        points, np.full(n, 100.0), np.zeros(n, dtype=bool), openness, dist_to_land, 0.0, source
    )
    assert np.allclose(land_only, source)


def test_spread_coastal_infill_prefers_sheltered_shallow_and_barrier_nodes():
    # Three ocean sinks in a row just west of a land source; only geometry/openness differ.
    lat = np.radians(np.array([0.0, 0.0, 0.0, 0.0]))
    lon = np.radians(np.array([0.30, 0.10, -0.10, -0.30]))  # source, sheltered, exposed-deep, barrier
    points = geometry.latlon_to_xyz(lat, lon)
    elevation = np.array([20.0, -10.0, -2000.0, -1.0])
    is_ocean = np.array([False, True, True, True])
    openness = np.array([0.0, 0.1, 0.9, 0.4])  # sheltered bay / open abyss / barrier edge
    dist_to_land = np.array([0.0, 0.05, 0.05, 0.001])  # barrier node hugs the coast
    source = np.array([100.0, 0.0, 0.0, 0.0])

    out = erosion._spread_coastal_infill(points, elevation, is_ocean, openness, dist_to_land, 0.0, source)
    assert out[1] > 0.0  # sheltered shallow water silts up
    assert out[3] > 0.0  # barrier candidate accretes despite facing open water
    assert out[2] < out[1] and out[2] < out[3]  # exposed deep abyss gets far less
    assert np.isclose(out.sum(), 100.0)


def test_apply_erosion_coastal_feedback_keeps_a_generated_world_sane():
    world = generate_world(seed=23, num_plates=8)
    _, before, _, _, _, _ = erosion._gather_nodes(world)
    erosion.apply_erosion(world, years=2_000_000)
    _, after, _, _, _, _ = erosion._gather_nodes(world)
    assert np.all(np.isfinite(after))
    assert np.all(after >= MIN_ELEVATION_M - 1e-6) and np.all(after <= MAX_ELEVATION_M + 1e-6)
    # The feedback nudges the coast, it doesn't rewrite the whole map in one step.
    assert np.median(np.abs(after - before)) < 50.0

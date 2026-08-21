import numpy as np
from app import climate, geometry
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def _world(seed=1, num_plates=12, continental_fraction=0.7, land_fraction=0.29, steps=0, years=5_000_000):
    world = generate_world(seed, num_plates=num_plates, continental_fraction=continental_fraction, land_fraction=land_fraction)
    for _ in range(steps):
        step_world(world, years)
    return world


def test_compute_climate_handles_extreme_land_ocean_fractions():
    # All-ocean and mostly-land worlds shouldn't crash or produce NaNs -- these exercise the
    # "no land"/"almost no ocean" guard branches in swirl, coastal deflection, and humidity.
    all_ocean = climate.compute_climate(_world(seed=2, continental_fraction=0.0, land_fraction=0.0))
    assert np.all(all_ocean.is_ocean)
    assert np.all(np.isfinite(all_ocean.current_u))

    mostly_land = climate.compute_climate(_world(seed=3, continental_fraction=1.0, land_fraction=0.9))
    assert np.all(np.isfinite(mostly_land.humidity))


def test_submerged_continental_crust_is_treated_as_ocean():
    # A continental plate whose elevation is everywhere below sea level (a fully-submerged
    # continent, not just a shelf) should read as ocean everywhere -- is_ocean is derived
    # from elevation, not crust_type, so evaporation/currents/coastal effects apply to it
    # exactly like any other ocean cell.
    frame = geometry.plate_frame_from_seed(np.array([1.0, 0.0, 0.0]))
    lines = [
        ElevationLine(phi=float(phi), theta=np.linspace(-np.pi, np.pi, 30, endpoint=False), elevation=np.full(30, -500.0))
        for phi in np.linspace(-1.4, 1.4, 15)
    ]
    plate = Plate(plate_id=0, frame=frame, crust_type="continental", lines=lines)
    world = World(seed=1, plates=[plate])

    fields = climate.compute_climate(world, height=30, width=60)
    assert np.all(fields.is_ocean)
    # Evaporation should actually happen there -- an all-land reading would leave humidity at
    # 0 everywhere, with no ocean source anywhere on the grid to draw from.
    assert np.any(fields.humidity > 0.0)


def test_insolation_peaks_at_equator_and_falls_toward_poles():
    lat_deg = np.array([-90.0, -45.0, 0.0, 45.0, 90.0])
    insolation = climate.compute_insolation(lat_deg, axial_tilt_deg=0.0)
    assert insolation[2] == insolation.max()
    assert insolation[0] < insolation[2] and insolation[4] < insolation[2]
    # Symmetric with no tilt.
    assert np.isclose(insolation[0], insolation[4])
    assert np.isclose(insolation[1], insolation[3])


def test_axial_tilt_flattens_the_insolation_curve():
    lat_deg = np.linspace(-90, 90, 37)
    flat = climate.compute_insolation(lat_deg, axial_tilt_deg=0.0)
    tilted = climate.compute_insolation(lat_deg, axial_tilt_deg=23.5)
    # Tilt brings the poles' annual-mean insolation up (they catch some direct sun over the
    # sweep) without changing the equator much.
    pole_idx = 0
    equator_idx = len(lat_deg) // 2
    assert tilted[pole_idx] > flat[pole_idx]
    assert abs(tilted[equator_idx] - flat[equator_idx]) < abs(tilted[pole_idx] - flat[pole_idx])


def test_land_temperature_warmer_at_equator_than_poles_and_cooler_with_elevation():
    lat_deg = np.array([-90.0, 0.0, 90.0])
    insolation = climate.compute_insolation(lat_deg, axial_tilt_deg=0.0)
    flat_elevation = np.zeros((3, 1))
    temp_flat = climate.compute_land_temperature(insolation, flat_elevation)
    assert temp_flat[1, 0] > temp_flat[0, 0]
    assert temp_flat[1, 0] > temp_flat[2, 0]

    high_elevation = np.full((3, 1), 4000.0)
    temp_high = climate.compute_land_temperature(insolation, high_elevation)
    assert temp_high[1, 0] < temp_flat[1, 0]  # lapse-rate cooling


def test_wind_coriolis_deflection_flips_sign_by_hemisphere():
    # Coriolis parameter (sin(lat)) is antisymmetric about the equator, so the zonal
    # deflection of an otherwise-identical meridional flow should flip sign north vs south.
    north = climate.coriolis_parameter(np.array([30.0]))
    south = climate.coriolis_parameter(np.array([-30.0]))
    assert north > 0 and south < 0
    assert np.isclose(north, -south)


def test_mountain_deflection_cancels_into_slope_component():
    # A flat wind field blowing straight into a single Gaussian-bump mountain, centered in
    # the grid, should have its into-slope component reduced right at the upslope face.
    height, width = 40, 80
    elevation = np.zeros((height, width))
    cy, cx = height // 2, width // 2
    yy, xx = np.mgrid[0:height, 0:width]
    elevation += 4000.0 * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (6.0 ** 2)))
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)

    u = np.full((height, width), 5.0)  # blowing due east
    v = np.zeros((height, width))
    u2, v2 = climate._mountain_deflection(u, v, elevation, lat_deg)

    # Just upwind (west) of the peak, the eastward (into-slope) component should have
    # dropped relative to the undeflected field.
    probe = (cy, cx - 4)
    assert u2[probe] < u[probe]


def test_coastal_deflection_cancels_into_land_component():
    height, width = 30, 60
    is_ocean = np.ones((height, width), dtype=bool)
    is_ocean[:, width // 2 :] = False  # land fills the eastern half
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)

    u = np.full((height, width), 5.0)  # blowing due east, straight into the coastline
    v = np.zeros((height, width))
    u2, v2 = climate._deflect_into_coast(u, v, is_ocean, lat_deg)

    coast_col = width // 2 - 1
    probe_row = height // 2
    assert u2[probe_row, coast_col] < u[probe_row, coast_col]


def test_land_swirl_current_is_zero_with_no_land():
    height, width = 20, 40
    is_ocean = np.ones((height, width), dtype=bool)
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_grid = np.repeat(lon_deg[None, :], height, axis=0)
    from app import geometry

    world_xyz = geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))
    wind_u = np.full((height, width), 5.0)
    wind_v = np.zeros((height, width))

    swirl_u, swirl_v = climate._land_swirl_current(wind_u, wind_v, is_ocean, lat_deg, world_xyz)
    assert np.allclose(swirl_u, 0.0)
    assert np.allclose(swirl_v, 0.0)


def test_circumglobal_boost_applies_only_to_fully_open_rows():
    height, width = 10, 20
    is_ocean = np.ones((height, width), dtype=bool)
    is_ocean[3, 5] = False  # row 3 has one land cell -- not fully open
    boost = climate._circumglobal_row_boost(is_ocean)
    assert np.allclose(boost[3], 1.0)
    assert np.allclose(boost[0], climate.CIRCUMGLOBAL_SPEEDUP_FACTOR)


def test_ocean_swells_detects_convergence_between_opposing_currents():
    height, width = 20, 40
    is_ocean = np.ones((height, width), dtype=bool)
    current_u = np.zeros((height, width))
    # Flow converges on the seam at column width//2: eastward on the left, westward on the
    # right of it.
    current_u[:, : width // 2] = 3.0
    current_u[:, width // 2 :] = -3.0
    current_v = np.zeros((height, width))

    rng = np.random.default_rng(0)
    rows, cols = climate.compute_ocean_swells(current_u, current_v, is_ocean, rng)
    assert len(rows) > 0
    # Every sampled point should be near the convergence seam, not scattered randomly.
    assert np.all(np.abs(cols - width // 2) <= 2)


def test_humidity_is_higher_near_warm_ocean_and_decays_inland():
    height, width = 30, 60
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    is_ocean = np.ones((height, width), dtype=bool)
    is_ocean[:, width // 2 :] = False  # land fills the eastern half
    elevation = np.where(is_ocean, 0.0, 200.0)
    ocean_temperature = np.full((height, width), 25.0)  # warm ocean everywhere
    wind_u = np.full((height, width), 3.0)  # blows eastward -- onto the land
    wind_v = np.zeros((height, width))
    elevation_factor = np.ones((height, width))

    humidity, _ = climate.compute_humidity(is_ocean, elevation, ocean_temperature, wind_u, wind_v, elevation_factor, lat_deg)
    # The zonal sweep's transport direction is the fixed latitude-band lookup
    # (zonal_direction_for_lat), not the literal sign of wind_u -- pick a row where that band
    # direction actually runs west-to-east (ocean -> land, increasing column) so "coast" and
    # "inland" mean what they say.
    row = next(r for r in range(height) if climate.zonal_direction_for_lat(lat_deg[r : r + 1])[0] == 1.0)
    coast_col = width // 2
    inland_col = width - 1
    assert humidity[row, coast_col] > humidity[row, inland_col]


def test_precipitation_increases_with_humidity_and_orographic_lift():
    low_humidity = climate.compute_precipitation(np.array([0.1]), np.array([0.0]))
    high_humidity = climate.compute_precipitation(np.array([1.0]), np.array([0.0]))
    assert high_humidity[0] > low_humidity[0]

    no_lift = climate.compute_precipitation(np.array([0.5]), np.array([0.0]))
    with_lift = climate.compute_precipitation(np.array([0.5]), np.array([0.2]))
    assert with_lift[0] > no_lift[0]


def test_air_temperature_pulled_toward_nearby_ocean_over_far_inland():
    from app import geometry

    height, width = 30, 60
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_grid = np.repeat(lon_deg[None, :], height, axis=0)
    world_xyz = geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))

    is_ocean = np.zeros((height, width), dtype=bool)
    is_ocean[:, :3] = True  # a thin ocean strip on the west edge

    land_temperature = np.full((height, width), -40.0)  # cold land baseline everywhere
    ocean_temperature = np.full((height, width), 20.0)  # warm ocean

    air_temperature = climate.compute_air_temperature(land_temperature, ocean_temperature, is_ocean, world_xyz)
    row = height // 2
    near_coast_col = 4
    # Longitude wraps, so the literal last column (width - 1) sits right next to the ocean
    # strip at column 0 -- not actually "far inland". The antipodal column (~180 degrees of
    # longitude from the ocean strip) is the genuinely farthest point on the sphere.
    far_inland_col = width // 2 + 1
    # Land right next to the ocean should read warmer (pulled toward it) than land far away.
    assert air_temperature[row, near_coast_col] > air_temperature[row, far_inland_col]
    assert air_temperature[row, far_inland_col] < ocean_temperature[row, far_inland_col]


def test_compute_climate_cached_populates_and_reuses_cache():
    world = _world(seed=2, num_plates=8, steps=0)
    assert world.climate_cache is None

    fields = climate.compute_climate_cached(world)
    assert world.climate_cache is fields  # populated by the call itself

    # A second call must reuse the same object, not compute a fresh one -- confirmed by
    # identity (`is`), not just equal values, since a fresh computation would also produce
    # equal-looking fields for an unchanged world.
    again = climate.compute_climate_cached(world)
    assert again is fields


def test_compute_climate_cached_ignores_stale_cache_by_design():
    # Deliberately not a correctness guarantee -- see climate.py's own docstring. A manually
    # planted sentinel proves compute_climate_cached really does trust whatever's already in
    # the cache rather than checking it against current world state.
    world = _world(seed=3, num_plates=8, steps=0)
    sentinel = climate.compute_climate(world)
    world.climate_cache = sentinel
    assert climate.compute_climate_cached(world) is sentinel

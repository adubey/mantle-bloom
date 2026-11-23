from types import SimpleNamespace

import numpy as np
import pytest

from app.hydroclimate import biomes, climate, sea_ice


def test_annual_sea_ice_balance_and_half_year_threshold():
    freeze = sea_ice.SEAWATER_FREEZING_C
    temperature = np.array([freeze, freeze + 0.1, freeze, freeze - 0.1, -25, -25])
    amplitude = np.array([0, 12, 12, 12, 12, 12])
    ocean = np.array([True, True, True, True, True, False])

    ice = sea_ice.annual_balance(temperature, amplitude, ocean)

    assert ice.seasonal.tolist() == [False, True, True, True, False, False]
    assert ice.perennial.tolist() == [False, False, False, False, True, False]
    assert ice.fraction[[0, 2, 4, 5]].tolist() == [0.0, 0.5, 1.0, 0.0]
    assert (ice.fraction >= sea_ice.POLAR_SEA_ICE_FRACTION).tolist() == [False, False, True, True, True, False]


@pytest.mark.parametrize(
    ("tilt", "perennial_edge", "freezing_edge"),
    [(15.0, 64.375, 54.875), (23.5, 67.875, 51.125), (35.0, 81.875, 45.875)],
)
def test_g5_zonal_edges_match_the_gate(tilt, perennial_edge, freezing_edge):
    # Match the 0.25° cell centers in analysis/issue548/zonal_bands.json.
    latitude = np.arange(0.125, 90.0, 0.25)
    heating = climate.compute_heating_index(latitude, tilt, 1.0)
    temperature = climate.compute_land_temperature(heating, np.zeros((len(latitude), 1)))[:, 0]
    amplitude = biomes.seasonal_temp_amplitude(latitude, np.zeros_like(latitude), tilt)
    ice = sea_ice.annual_balance(temperature, amplitude, np.ones(latitude.shape, dtype=bool))

    assert latitude[ice.perennial].min() == perennial_edge
    assert latitude[ice.seasonal].min() == freezing_edge


@pytest.mark.parametrize("cooling", [False, True])
def test_current_anomaly_shifts_the_perennial_edge(cooling):
    world = SimpleNamespace(
        axial_tilt_deg=23.5, solar_multiplier=1.0,
        ice_age_period_years=1e6 if cooling else 0.0, elapsed_years=0.5e6,
    )
    latitude = np.array([55.5 if cooling else 67.9])
    heating = climate.compute_heating_index(latitude, world.axial_tilt_deg, world.solar_multiplier)
    baseline = climate.compute_ocean_temperature_baseline(heating, 1, 2) - climate.ice_age_cooling_c(world)
    fields = SimpleNamespace(
        lat_deg=latitude,
        ocean_temperature_c=baseline + np.array([[2.7, -2.7]]),
        is_ocean=np.ones((1, 2), dtype=bool),
    )
    ice = sea_ice.sample(world, fields, np.array([0, 0]), np.array([0, 1]), fields.is_ocean[0])

    assert ice.temperature_c[0] - ice.temperature_c[1] == pytest.approx(5.4)
    assert ice.perennial.tolist() == [False, True]


def test_sample_smooth_interpolates_between_coarse_climate_cells():
    world = SimpleNamespace(
        axial_tilt_deg=23.5, solar_multiplier=1.0,
        ice_age_period_years=0.0, elapsed_years=0.0,
    )
    # 6-row equirectangular grid has cell centers at +75, +45, +15, -15, -45, -75 deg.
    lat_deg = np.array([75.0, 45.0, 15.0, -15.0, -45.0, -75.0])
    sst_rows = np.array([-1.5, 10.5, 25.0, 25.0, 10.5, -1.5])
    fields = SimpleNamespace(
        lat_deg=lat_deg,
        ocean_temperature_c=np.repeat(sst_rows[:, None], 2, axis=1),
        is_ocean=np.ones((6, 2), dtype=bool),
    )
    sample_lats = np.array([75.0, 60.0, 45.0])
    sample_lons = np.array([0.0, 0.0, 0.0])
    smooth = sea_ice.sample_smooth(world, fields, sample_lats, sample_lons, np.ones(3, dtype=bool))

    # Intermediate latitude gets an interpolated temperature rather than snapping to either coarse row.
    assert smooth.temperature_c[0] < smooth.temperature_c[1] < smooth.temperature_c[2]
    assert smooth.perennial_cover[0] == pytest.approx(1.0)
    assert smooth.perennial_cover[2] == pytest.approx(0.0)


def test_perennial_cover_ramps_smoothly_across_perennial_edge():
    lat = np.array([65.0, 67.5, 72.0])
    heating = climate.compute_heating_index(lat, 23.5, 1.0)
    temperature = climate.compute_land_temperature(heating, np.zeros((len(lat), 1)))[:, 0]
    amplitude = biomes.seasonal_temp_amplitude(lat, np.zeros_like(lat), 23.5)
    ice = sea_ice.annual_balance(temperature, amplitude, np.ones_like(lat, dtype=bool))

    assert ice.perennial_cover[0] == 0.0
    assert 0.0 < ice.perennial_cover[1] < 1.0
    assert ice.perennial_cover[2] == 1.0

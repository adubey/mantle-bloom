"""Checks for geography-only eligibility and area-weighted packet summaries."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'bin' / 'debug'))
import rfc0019_g3b_calibration_packet as packet  # noqa: E402


def test_exposure_looks_upwind_and_wraps_longitude():
    ocean = np.zeros((360, 720), dtype=bool)
    ocean[180, 719] = True
    u = np.ones_like(ocean, dtype=float)
    v = np.zeros_like(u)
    lat = 90 - (np.arange(360) + 0.5) * 0.5
    result = packet.land_exposure(ocean, u, v, lat)
    assert result[180, 0] == 1
    assert result[180, 20] == 0
    u[180, 20] = 0
    assert packet.land_exposure(ocean, u, v, lat)[180, 20] == -1


def test_profiles_and_strata_use_area_and_geography_only():
    rain = np.array([[10., 20., 30., 40.]])
    area = np.array([[1., 3., 2., 4.]])
    ocean = np.array([[True, False, False, False]])
    distance = np.array([[0., 100., 1600., 1700.]])
    exposure = np.array([[-1, 1, 1, 0]])
    profiles, comparisons = packet.summaries(rain, area, ocean, distance, exposure, np.array([2.]))
    row = next(row for row in profiles if row['latitude_deg'] == [0, 5])
    assert row['all']['rain_mm_yr'] == 29.0
    assert row['land']['rain_mm_yr'] == 280 / 9
    interior = next(row for row in comparisons if row['latitude_deg'] == [0, 15]
                    and row['exposure'] == 'ocean_upwind_500km' and row['distance'] == 'interior')
    assert interior['rain_mm_yr'] == 30
    assert interior['cells'] == 1

import numpy as np
from app import boundary
from app.plates import _band_intensity


def test_closing_rate_positive_when_approaching():
    # Two points near the equator, close together in longitude.
    p_self = np.array([np.cos(0.0), np.sin(0.0), 0.0])
    p_neighbor = np.array([np.cos(0.01), np.sin(0.01), 0.0])

    # self plate spins +z (eastward motion at the equator) -> moves toward neighbor (east of it)
    approaching_omega = np.array([0.0, 0.0, 1.0])
    still_omega = np.zeros(3)
    closing = boundary.closing_rate(p_self[None, :], approaching_omega, still_omega, p_neighbor[None, :])
    assert closing[0] > 0

    # self plate spins -z (westward motion) -> moves away from neighbor
    receding_omega = np.array([0.0, 0.0, -1.0])
    closing = boundary.closing_rate(p_self[None, :], receding_omega, still_omega, p_neighbor[None, :])
    assert closing[0] < 0


def test_band_intensity_zero_outside_band_peaks_at_midpoint():
    dist = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    result = _band_intensity(dist, inner=0.1, outer=0.3)
    assert result[0] == 0.0  # below inner edge
    assert np.isclose(result[1], 0.0)  # exactly at inner edge
    assert result[2] == 1.0  # at the band's midpoint
    assert np.isclose(result[3], 0.0)  # exactly at outer edge
    assert result[4] == 0.0  # beyond outer edge
    assert np.all((result >= 0.0) & (result <= 1.0))

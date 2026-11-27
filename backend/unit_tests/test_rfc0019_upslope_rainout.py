"""Equation and water-accounting checks for the experimental upslope rule."""

from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin" / "debug"))

from rfc0019_coarse_moisture import TransportFaces, solve_bounded_budget  # noqa: E402
from rfc0019_g3_ridge_rainout import _saturation_scale_height_m, _upslope_interception  # noqa: E402

TEMP_K = np.full(4, 288.0)
LAPSE_K_PER_KM = 6.0


def _faces(first, second):
    count = len(first)
    return TransportFaces(np.array(first), np.array(second), np.ones(count),
                          np.zeros(count), np.zeros(count))


def test_upslope_rainout_direction_and_land_mask():
    grid = SimpleNamespace(is_ocean=np.array([False, False, False, True]))
    faces = _faces([0, 1, 1, 2], [1, 2, 0, 3])
    result = _upslope_interception(grid, faces, np.array([0., 250., 500., 0.]), TEMP_K, LAPSE_K_PER_KM)
    hsat = _saturation_scale_height_m(288.0, LAPSE_K_PER_KM)
    fraction = 1 - np.exp(-250 / hsat)
    np.testing.assert_allclose(result.forward, [fraction, fraction, 0, 0])
    np.testing.assert_allclose(result.reverse, [0, 0, fraction, 0])


def test_subdivided_climb_preserves_transmitted_fraction():
    grid = SimpleNamespace(is_ocean=np.zeros(3, dtype=bool))
    faces = _faces([0, 1, 0], [1, 2, 2])
    result = _upslope_interception(grid, faces, np.array([0., 125., 250.]), np.full(3, 288.), LAPSE_K_PER_KM)
    assert (1 - result.forward[0]) * (1 - result.forward[1]) == pytest.approx(1 - result.forward[2])
    flat = _upslope_interception(grid, faces, np.zeros(3), np.full(3, 288.), LAPSE_K_PER_KM)
    assert np.all(flat.forward == 0)
    with pytest.raises(ValueError):
        _upslope_interception(grid, faces, np.zeros(3), np.zeros(3), LAPSE_K_PER_KM)


def test_upslope_interception_books_water_as_rain():
    grid = SimpleNamespace(is_ocean=np.zeros(2, dtype=bool))
    faces = _faces([0], [1])
    interception = _upslope_interception(grid, faces, np.array([0., 250.]), np.full(2, 288.), LAPSE_K_PER_KM)
    budget, _, _ = solve_bounded_budget(
        np.ones(2), np.array([1., 0.]), np.ones(2), faces, np.full(2, 10.),
        interception=interception,
    )
    assert budget.intercepted_rain_kg_s[1] > 0
    assert budget.rain.sum() == pytest.approx(1.0)
    assert budget.global_relative < 1e-10

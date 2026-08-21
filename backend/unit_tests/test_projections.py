import numpy as np

from app import projections


def test_behrmann_origin():
    x, y = projections.behrmann(0.0, 0.0)
    assert np.isclose(x, 0.0)
    assert np.isclose(y, 0.0)


def test_behrmann_known_points():
    cos30 = np.cos(np.radians(30.0))
    x, y = projections.behrmann(np.radians(90.0), np.radians(180.0))
    assert np.isclose(y, 1.0 / cos30, atol=1e-9)
    assert np.isclose(x, np.pi * cos30, atol=1e-9)


def test_behrmann_is_linear_in_longitude():
    lat = np.radians(15.0)
    x1, _ = projections.behrmann(lat, np.radians(10.0))
    x2, _ = projections.behrmann(lat, np.radians(20.0))
    x3, _ = projections.behrmann(lat, np.radians(30.0))
    assert np.isclose(x2 - x1, x3 - x2, atol=1e-9)


def test_eckert4_origin():
    x, y = projections.eckert4(0.0, 0.0)
    assert np.isclose(x, 0.0, atol=1e-9)
    assert np.isclose(y, 0.0, atol=1e-9)


def test_eckert4_known_pole_value():
    c = 2.0 / np.sqrt(np.pi * (4.0 + np.pi))
    x, y = projections.eckert4(np.radians(90.0), 0.0)
    assert np.isclose(x, 0.0, atol=1e-9)
    assert np.isclose(y, c * np.pi, atol=1e-6)


def test_eckert4_antisymmetric_in_latitude():
    lon = np.radians(37.0)
    _, y_pos = projections.eckert4(np.radians(40.0), lon)
    _, y_neg = projections.eckert4(np.radians(-40.0), lon)
    assert np.isclose(y_pos, -y_neg, atol=1e-9)


def test_eckert4_vectorized_matches_scalar():
    lats = np.radians(np.array([-89.0, -45.0, 0.0, 30.0, 89.0]))
    lons = np.radians(np.array([-170.0, -10.0, 0.0, 45.0, 179.0]))
    x_vec, y_vec = projections.eckert4(lats, lons)
    for i, (lat, lon) in enumerate(zip(lats, lons)):
        x_s, y_s = projections.eckert4(lat, lon)
        assert np.isclose(x_vec[i], x_s, atol=1e-9)
        assert np.isclose(y_vec[i], y_s, atol=1e-9)


def test_project_registry_dispatch():
    x, y = projections.project("behrmann", 0.0, 0.0)
    assert np.isclose(x, 0.0) and np.isclose(y, 0.0)
    x, y = projections.project("eckert4", 0.0, 0.0)
    assert np.isclose(x, 0.0) and np.isclose(y, 0.0)


def test_project_unknown_raises():
    try:
        projections.project("mercator", 0.0, 0.0)
        assert False, "expected ValueError"
    except ValueError:
        pass

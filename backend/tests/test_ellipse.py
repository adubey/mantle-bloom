import numpy as np

from app.ellipse import min_enclosing_ellipse


def _max_containment(points_xy, fit):
    cos_a, sin_a = np.cos(fit.angle_rad), np.sin(fit.angle_rad)
    rel = points_xy - fit.center
    u = rel[:, 0] * cos_a + rel[:, 1] * sin_a
    v = -rel[:, 0] * sin_a + rel[:, 1] * cos_a
    semi_major = max(fit.semi_major, 1e-12)
    semi_minor = max(fit.semi_minor, 1e-12)
    return float((((u / semi_major) ** 2) + ((v / semi_minor) ** 2)).max())


def test_contains_every_point_on_random_clouds():
    rng = np.random.default_rng(0)
    for _ in range(20):
        n = int(rng.integers(3, 60))
        points = rng.normal(size=(n, 2)) * rng.uniform(0.1, 50, size=2)
        fit = min_enclosing_ellipse(points)
        assert _max_containment(points, fit) <= 1.0 + 1e-6


def test_matches_the_closed_form_rectangle_oracle():
    # The true MVEE of an axis-aligned rectangle's 4 corners (+-a, +-b) is analytically known:
    # semi_major = max(a,b)*sqrt(2), semi_minor = min(a,b)*sqrt(2), axis-aligned.
    a, b = 5.0, 2.0
    points = np.array([[a, b], [a, -b], [-a, b], [-a, -b]])
    fit = min_enclosing_ellipse(points)
    assert np.isclose(fit.semi_major, a * np.sqrt(2), rtol=1e-4)
    assert np.isclose(fit.semi_minor, b * np.sqrt(2), rtol=1e-4)
    assert np.allclose(fit.center, [0.0, 0.0], atol=1e-6)


def test_single_point_is_a_zero_size_ellipse():
    fit = min_enclosing_ellipse(np.array([[3.0, -2.0]]))
    assert np.allclose(fit.center, [3.0, -2.0])
    assert fit.semi_major == 0.0
    assert fit.semi_minor == 0.0


def test_two_points_degenerate_to_a_segment():
    fit = min_enclosing_ellipse(np.array([[0.0, 0.0], [4.0, 0.0]]))
    assert np.allclose(fit.center, [2.0, 0.0])
    assert np.isclose(fit.semi_major, 2.0)
    assert fit.semi_minor == 0.0
    assert np.isclose(fit.angle_rad, 0.0, atol=1e-9)


def test_collinear_points_do_not_crash_and_stay_contained():
    points = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1.5, 0.0], [0.5, 0.0]])
    fit = min_enclosing_ellipse(points)
    assert np.isfinite(fit.semi_major)
    assert np.isfinite(fit.semi_minor)
    assert _max_containment(points, fit) <= 1.0 + 1e-6


def test_all_identical_points_do_not_crash():
    points = np.full((5, 2), [3.0, 3.0])
    fit = min_enclosing_ellipse(points)
    assert np.allclose(fit.center, [3.0, 3.0])
    assert fit.semi_major == 0.0
    assert fit.semi_minor == 0.0

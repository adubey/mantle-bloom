import numpy as np
import pytest

from app import terrain_noise
from app.terrain_noise import ContinentalRelief, FractalTexture, OceanicRelief, _fibonacci_sphere

_UPLIFT_KW = dict(orogenic_units=2.4, plateau_units=1.4, plateau_relief_units=0.45)


def _points(n=2000, seed=0):
    return _fibonacci_sphere(n)


def test_continental_relief_duck_types_sphere_noise_and_is_bounded():
    r = ContinentalRelief(np.random.default_rng(1), **_UPLIFT_KW)
    pts = _points()
    s = r.sample(pts)
    u = r.uplift(pts)
    assert s.shape == (len(pts),) and np.isfinite(s).all()
    assert u.shape == (len(pts),) and np.isfinite(u).all()
    # sample() keeps roughly the old single-SphereNoise character; uplift() is non-negative.
    assert abs(np.std(s) - terrain_noise._TARGET_SAMPLE_STD) < 0.05
    assert u.min() >= 0.0
    assert u.max() < 6.0  # a few units at most -- see the _*_RELIEF budget


def test_continental_relief_shape_preserved_for_2d_input():
    r = ContinentalRelief(np.random.default_rng(2), **_UPLIFT_KW)
    grid = _fibonacci_sphere(90).reshape(9, 10, 3)
    assert r.sample(grid).shape == (9, 10)
    assert r.uplift(grid).shape == (9, 10)


@pytest.mark.parametrize("cls,kw", [(ContinentalRelief, _UPLIFT_KW), (OceanicRelief, {}), (FractalTexture, {})])
def test_relief_fields_are_deterministic_in_rng(cls, kw):
    pts = _points()
    a = cls(np.random.default_rng(7), **kw)
    b = cls(np.random.default_rng(7), **kw)
    assert np.array_equal(a.sample(pts), b.sample(pts))
    if hasattr(a, "uplift"):
        assert np.array_equal(a.uplift(pts), b.uplift(pts))


def test_uplift_has_localized_belts_not_uniform_lift():
    # A real orogenic field is concentrated: most of the sphere gets little uplift, a
    # minority carries the ranges. A uniform lift would fail this.
    r = ContinentalRelief(np.random.default_rng(3), **_UPLIFT_KW)
    u = r.uplift(_points(6000))
    assert np.median(u) < 0.5 * np.quantile(u, 0.95)
    assert np.mean(u > np.quantile(u, 0.9)) < 0.2


def test_oceanic_relief_stays_gentle():
    o = OceanicRelief(np.random.default_rng(4))
    v = o.sample(_points(6000))
    # Abyssal hills plus the odd seamount -- never a whole unit (~2 km) of relief on average.
    assert np.std(v) < 0.5
    assert np.quantile(np.abs(v), 0.99) < 1.5

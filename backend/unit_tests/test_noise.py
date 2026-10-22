import numpy as np

from app.noise import IsotropicSphereNoise, SphereNoise
from app.terrain_noise import _fibonacci_sphere


def _structure_by_axis(noise, points: np.ndarray, step: float = 0.05) -> np.ndarray:
    """Mean squared change of the field for a small step along each world axis."""
    base = noise.sample(points)
    out = []
    for axis in np.eye(3):
        moved = points + step * axis
        moved /= np.linalg.norm(moved, axis=1, keepdims=True)
        out.append(np.mean((noise.sample(moved) - base) ** 2))
    return np.asarray(out)


def test_isotropic_noise_is_deterministic_and_unit_scaled():
    points = _fibonacci_sphere(4096)
    first = IsotropicSphereNoise(np.random.default_rng(7), octaves=4, base_freq=6.0).sample(points)
    again = IsotropicSphereNoise(np.random.default_rng(7), octaves=4, base_freq=6.0).sample(points)
    other = IsotropicSphereNoise(np.random.default_rng(8), octaves=4, base_freq=6.0).sample(points)

    assert np.array_equal(first, again)
    assert not np.array_equal(first, other)
    assert 0.8 < np.std(first) < 1.2


def test_isotropic_noise_has_no_preferred_direction_unlike_single_wave_octaves():
    # #218: one plane wave per octave draws parallel stripes, so the field changes much
    # faster across the stripes than along them; many random waves per octave don't.
    points = _fibonacci_sphere(8192)
    isotropic = _structure_by_axis(IsotropicSphereNoise(np.random.default_rng(3), octaves=3, base_freq=8.0), points)
    striped = _structure_by_axis(SphereNoise(np.random.default_rng(3), octaves=1, base_freq=8.0), points)

    assert isotropic.max() / isotropic.min() < 1.5
    assert striped.max() / striped.min() > 3.0

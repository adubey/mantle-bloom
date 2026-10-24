"""Cheap smooth noise on the sphere: a small sum of sinusoids with random frequency
vectors/phases. Not a proper gradient-noise implementation, but continuous and good enough
for giving plate interiors some non-flat texture without adding an extra dependency.
"""

from __future__ import annotations

import numpy as np


class SphereNoise:
    def __init__(self, rng: np.random.Generator, octaves: int = 3, base_freq: float = 2.0):
        self._freqs = []
        self._phases = []
        self._amps = []
        freq = base_freq
        amp = 1.0
        total_amp = 0.0
        for _ in range(octaves):
            direction = rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            self._freqs.append(direction * freq)
            self._phases.append(rng.uniform(0, 2 * np.pi))
            self._amps.append(amp)
            total_amp += amp
            freq *= 2.0
            amp *= 0.5
        self._norm = 1.0 / total_amp if total_amp > 0 else 1.0

    def sample(self, xyz: np.ndarray) -> np.ndarray:
        """xyz: (..., 3) world/unit-sphere positions -> noise values in roughly [-1, 1]."""
        xyz = np.asarray(xyz, dtype=float)
        total = np.zeros(xyz.shape[:-1])
        for freq_vec, phase, amp in zip(self._freqs, self._phases, self._amps):
            total = total + amp * np.sin(xyz @ freq_vec + phase)
        return total * self._norm


class IsotropicSphereNoise:
    """Smooth noise with no preferred direction: each octave sums `waves` plane waves with
    random directions and phases, approximately a Gaussian random field. `SphereNoise` uses
    one wave per octave, so its relief runs in parallel stripes across a whole continent
    (#218). Each octave is scaled to unit variance before the `gain` falloff, and the sum is
    renormalized, so `sample` has a standard deviation close to 1 (unbounded, unlike
    `SphereNoise`'s [-1, 1])."""

    def __init__(
        self, rng: np.random.Generator, octaves: int = 3, base_freq: float = 2.0, waves: int = 32, gain: float = 0.5
    ):
        self._octaves = []
        freq, amp, total_var = base_freq, 1.0, 0.0
        for _ in range(octaves):
            directions = rng.normal(size=(waves, 3))
            directions /= np.linalg.norm(directions, axis=1, keepdims=True)
            phases = rng.uniform(0, 2 * np.pi, waves)
            # A random-phase sine has variance 1/2, so `waves` of them sum to variance waves/2.
            self._octaves.append((directions * freq, phases, amp / np.sqrt(waves / 2.0)))
            total_var += amp * amp
            freq *= 2.0
            amp *= gain
        self._norm = 1.0 / np.sqrt(total_var) if total_var > 0 else 1.0

    def sample(self, xyz: np.ndarray) -> np.ndarray:
        """xyz: (..., 3) unit-sphere positions -> noise values with standard deviation ~1."""
        xyz = np.asarray(xyz, dtype=float)
        total = np.zeros(xyz.shape[:-1])
        for freq_vecs, phases, amp in self._octaves:
            total = total + amp * np.sin(xyz @ freq_vecs.T + phases).sum(axis=-1)
        return total * self._norm

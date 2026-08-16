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

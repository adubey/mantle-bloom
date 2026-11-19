"""Physical glacial abrasion law shared by the ice march and terrain update (RFC 0015 D3)."""

from __future__ import annotations

import numpy as np

# In K * u_b**l, K has the compound units (m/yr) ** (1 - l), so the result
# remains an erosion rate in m/yr before integrating over the interval.
GLACIAL_ABRASION_COEFFICIENT = 1.0e-4
GLACIAL_ABRASION_EXPONENT = 0.65


def abrasion_depth_m(sliding_speed_m_per_yr: np.ndarray, years: float, multiplier: float) -> np.ndarray:
    """Integrate D3's basal-sliding abrasion law over one march interval, in metres."""
    sliding = np.clip(np.asarray(sliding_speed_m_per_yr, dtype=float), 0.0, None)
    multiplier = max(float(multiplier), 0.0)
    return (
        GLACIAL_ABRASION_COEFFICIENT
        * multiplier
        * np.power(sliding, GLACIAL_ABRASION_EXPONENT)
        * years
    )

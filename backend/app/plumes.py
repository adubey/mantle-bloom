"""Mantle plume head footprints (RFC 0007; issue #577 scale-check slice).

Pure geometry only: no event stream, no `World` integration, no tracers or ledger.
The head footprint derives the emplacement area from the head volume through one mean
excess thickness, so the peak excess is exactly `V / A` and typical heads reproduce
the 2-3 km relief plateaus stand above the abyssal floor
[[Kerr & Mahoney, 2007]](#ref-kerr2007). Full design and knob table live in
`rfc/0007-mantle-plumes-and-plateaus.md`; the event stream, tails, weakening,
docking, and ledger arrive with the implementation issues tracked in #573.
"""

from __future__ import annotations

import math

import numpy as np

KM3_TO_M3 = 1e9
KM2_TO_M2 = 1e6

# Mean excess crust thickness tying head volume to footprint area: with the model's
# water-loaded oceanic isostasy each km of excess raises the surface by ~155 m, so
# 15 km gives ~2.3 km of relief, inside the observed 2-3 km plateau range (#577).
PLATEAU_MEAN_EXCESS_M = 15_000.0
# LIPs cover more than 0.1 M km²; plateaus reach 2 M km² (#577).
PLATEAU_AREA_MIN_KM2 = 1e5
PLATEAU_AREA_MAX_KM2 = 2e6
# Calibration assumption: the range spans typical ≥1 M km³ LIPs to Ontong Java's
# 44-57 M km³; the distribution shape is not established (#577).
PLUME_HEAD_VOLUME_MIN_KM3 = 1e6
PLUME_HEAD_VOLUME_MAX_KM3 = 5e7


def head_footprint_area_m2(volume_m3: float) -> float:
    """Oceanic head footprint area: volume tied to area by one mean excess thickness.

    Very small heads would otherwise vanish below the grid; very large ones are
    capped at the largest observed plateaus, accepting the overshoot the RFC records
    for Ontong-Java-scale heads (issue #577).
    """
    area = volume_m3 / PLATEAU_MEAN_EXCESS_M
    return float(
        np.clip(
            area,
            PLATEAU_AREA_MIN_KM2 * KM2_TO_M2,
            PLATEAU_AREA_MAX_KM2 * KM2_TO_M2,
        )
    )


def head_kernel_radius_m(area_m2: float) -> float:
    """Flat-topped kernel radius whose own integral over the plane equals the area.

    With `w = exp(-(d/r)^4)`, the plane integral is `pi * r^2 * Gamma(3/2)`, so this
    `r` makes the kernel place exactly the head volume at peak `V / A` (#577).
    """
    return math.sqrt(area_m2 / (math.pi * math.gamma(1.5)))


def head_excess_thickness_m(distances_m: np.ndarray, volume_m3: float) -> np.ndarray:
    """Excess crust thickness over the footprint: peak `V / A` with flat-topped falloff.

    The peak is exact by construction; the integral over cells equals the head volume
    up to the grid discretization, which the emplacement pass normalizes away (#577).
    """
    distances = np.asarray(distances_m, dtype=float)
    area = head_footprint_area_m2(volume_m3)
    radius = head_kernel_radius_m(area)
    peak = volume_m3 / area
    return peak * np.exp(-((distances / radius) ** 4))

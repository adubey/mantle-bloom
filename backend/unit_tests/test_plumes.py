import math

import numpy as np
import pytest

from app import lithosphere, plumes


def _median_head_m3() -> float:
    return (
        math.sqrt(
            plumes.PLUME_HEAD_VOLUME_MIN_KM3 * plumes.PLUME_HEAD_VOLUME_MAX_KM3
        )
        * plumes.KM3_TO_M3
    )


def _peak_relief_m(volume_m3: float) -> float:
    """Relief of the footprint peak over reference oceanic lithosphere (#577)."""
    hc0 = lithosphere.REFERENCE_HC_OCEANIC_M
    hm = lithosphere.REFERENCE_HM_OCEANIC_M
    rho = lithosphere.crust_density("oceanic")
    before = float(lithosphere.isostatic_elevation(hc0, hm, rho))
    peak_excess = volume_m3 / plumes.head_footprint_area_m2(volume_m3)
    after = float(lithosphere.isostatic_elevation(hc0 + peak_excess, hm, rho))
    return after - before


def test_median_head_gives_plateau_relief():
    """Isolated-head scale check: a median head makes 2.0-2.6 km of relief (#577)."""
    volume = _median_head_m3()
    area = plumes.head_footprint_area_m2(volume)
    assert (
        plumes.PLATEAU_AREA_MIN_KM2 * plumes.KM2_TO_M2
        < area
        < plumes.PLATEAU_AREA_MAX_KM2 * plumes.KM2_TO_M2
    )
    assert 2_000.0 <= _peak_relief_m(volume) <= 2_600.0


def test_min_and_max_heads_stay_within_relief_bounds():
    """Small heads are area-clipped, Ontong-scale heads capped, both sane (#577)."""
    for volume_km3 in (
        plumes.PLUME_HEAD_VOLUME_MIN_KM3,
        plumes.PLUME_HEAD_VOLUME_MAX_KM3,
    ):
        assert 1_200.0 <= _peak_relief_m(volume_km3 * plumes.KM3_TO_M3) <= 4_500.0


def test_kernel_places_scheduled_volume():
    """The flat-topped kernel integrates to the head volume (#577)."""
    volume = _median_head_m3()
    area = plumes.head_footprint_area_m2(volume)
    radius = plumes.head_kernel_radius_m(area)
    edges = np.linspace(0.0, 4.0 * radius, 4001)
    mid = 0.5 * (edges[:-1] + edges[1:])
    ring_area = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    placed = float(np.sum(plumes.head_excess_thickness_m(mid, volume) * ring_area))
    assert placed == pytest.approx(volume, rel=5e-3)
    assert float(plumes.head_excess_thickness_m(np.zeros(1), volume)[0]) == volume / area

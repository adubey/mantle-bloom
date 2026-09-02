"""Mohr-Coulomb convergent/divergent deformation (`app.rheology`).

The one property with real teeth here is calibration: `EFFECTIVE_LITHOSPHERE_VISCOSITY_PA_S_
PER_M` has to be large enough that a *realistic* continental-collision closing rate (a few
cm/yr) actually clears the yield stress and thickens crust. It was ~3-4 orders of magnitude
too small for most of 2026 -- `apply_convergent_deformation` returned Hc unchanged at every
plate speed the sim allows, so no collision ever built a mountain and continents only ever
thinned and drowned (land fraction fell monotonically; the `overlapAge` view's stalled
multi-plate collisions never crumpled into orogens). These tests pin the calibration.
"""

import numpy as np

from app import rheology
from app.lithosphere import PLANET_RADIUS_M

SECONDS_PER_YEAR = 365.25 * 86400.0


def _closing_m_per_s(cm_per_yr: float) -> float:
    return (cm_per_yr / 100.0) / SECONDS_PER_YEAR


def test_realistic_continental_collision_thickens_crust():
    """A sustained 3 cm/yr collision -- an ordinary continent-continent convergence rate,
    well within `mantle.MAX_PLATE_RATE` (15 cm/yr) -- must thicken Hc by a geologically
    meaningful amount over a 1-Myr step (hundreds of metres on a 35-km column), not leave it
    untouched."""
    hc = np.array([35_000.0])
    hm = np.array([100_000.0])
    closing = np.array([_closing_m_per_s(3.0)])

    new_hc, new_hm = rheology.apply_convergent_deformation(
        hc, hm, closing, years_myr=1.0, fault_factor=np.array([1.0])
    )

    assert new_hc[0] - hc[0] > 200.0
    # Hc and Hm thicken together (crustal shortening drags the mantle lid along).
    assert new_hm[0] > hm[0]
    assert (new_hc[0] / hc[0]) == np.float64(new_hm[0] / hm[0])


def test_slow_graze_stays_below_yield():
    """A ~0.5 cm/yr brush is not a collision -- it stays in the elastic regime and thickens
    nothing, so the yield check still discriminates real orogeny from plates merely sliding
    past each other."""
    hc = np.array([35_000.0])
    hm = np.array([100_000.0])
    closing = np.array([_closing_m_per_s(0.5)])

    new_hc, new_hm = rheology.apply_convergent_deformation(
        hc, hm, closing, years_myr=1.0, fault_factor=np.array([1.0])
    )

    assert new_hc[0] == hc[0]
    assert new_hm[0] == hm[0]


def test_plastic_strain_rate_monotonic_and_saturates():
    """Strain rate rises with closing rate then flattens -- the linear-viscosity stand-in is
    normalized by the yield stress so it can't blow up at the fast end (`mantle.MAX_PLATE_
    RATE` head-on is ~30 cm/yr closing)."""
    rates = np.array([_closing_m_per_s(c) for c in (0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0)])
    strain = rheology.plastic_strain_rate_per_myr(rates)

    assert np.all(np.diff(strain) >= 0.0)
    assert strain[0] == 0.0  # 0.5 cm/yr: sub-yield
    assert strain[3] > 0.0  # 3 cm/yr: real orogeny
    # Saturation: the fastest closing rate is only modestly above a moderate one, not 10x.
    assert strain[-1] < 4.0 * strain[3]


def test_extension_thins_crust_and_triggers_decompression_melting():
    """The divergent branch is the same yield check with the opposite sign -- a fast rift
    thins Hc, and a column dragged below `RIFT_CRITICAL_THICKNESS_M` flags for melting."""
    hc = np.array([5_500.0, 40_000.0])
    hm = np.array([30_000.0, 100_000.0])
    opening = np.array([-_closing_m_per_s(10.0), -_closing_m_per_s(10.0)])

    new_hc, _, melting = rheology.apply_divergent_deformation(hc, hm, opening, years_myr=10.0)

    assert np.all(new_hc < hc)
    assert bool(melting[0]) and not bool(melting[1])

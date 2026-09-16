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

from app import lithosphere, rheology
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

    new_hc, new_hm, overflow_hc = rheology.apply_convergent_deformation(
        hc, hm, closing, years_myr=1.0, fault_factor=np.array([1.0])
    )

    assert new_hc[0] - hc[0] > 200.0
    # Hc and Hm thicken together (crustal shortening drags the mantle lid along).
    assert new_hm[0] > hm[0]
    assert (new_hc[0] / hc[0]) == np.float64(new_hm[0] / hm[0])
    assert overflow_hc[0] == 0.0  # nowhere near the ceiling from one ordinary step


def test_slow_graze_stays_below_yield():
    """A ~0.5 cm/yr brush is not a collision -- it stays in the elastic regime and thickens
    nothing, so the yield check still discriminates real orogeny from plates merely sliding
    past each other."""
    hc = np.array([35_000.0])
    hm = np.array([100_000.0])
    closing = np.array([_closing_m_per_s(0.5)])

    new_hc, new_hm, overflow_hc = rheology.apply_convergent_deformation(
        hc, hm, closing, years_myr=1.0, fault_factor=np.array([1.0])
    )

    assert new_hc[0] == hc[0]
    assert new_hm[0] == hm[0]
    assert overflow_hc[0] == 0.0


def test_convergent_thickening_caps_hc_and_hm_and_reports_the_overflow():
    """GitHub issue #161: unbounded, no ceiling meant a node sitting in a long-lived
    convergent regime compounded Hc/Hm exponentially forever (measured up to Hc=413,885 m /
    Hm=1,311,592 m on a real long save). A node already sitting at (or past) the ceiling must
    not grow any further, and the function must report how much growth got clipped off so the
    caller can conserve it (spread it onto the near-field foreland -- see
    lithosphere_plate.deform) rather than silently losing it."""
    hc = np.array([lithosphere.MAX_CRUSTAL_THICKNESS_M, 35_000.0])
    hm = np.array([lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M, 100_000.0])
    closing = np.full(2, _closing_m_per_s(5.0))  # a fast, well-past-yield collision

    new_hc, new_hm, overflow_hc = rheology.apply_convergent_deformation(
        hc, hm, closing, years_myr=5.0, fault_factor=np.ones(2)
    )

    assert new_hc[0] == lithosphere.MAX_CRUSTAL_THICKNESS_M  # already at the ceiling: no further growth
    assert new_hm[0] == lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M
    assert overflow_hc[0] > 0.0  # but the strain that would have thickened it further is reported...
    assert new_hc[1] < lithosphere.MAX_CRUSTAL_THICKNESS_M  # ...while an ordinary column below the ceiling is untouched
    assert overflow_hc[1] == 0.0


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


def test_arc_magmatism_adds_juvenile_crust_scaled_by_band_and_convergence():
    """A subducting oceanic slab underplates the overriding continental margin with juvenile
    melt -- extra Hc (only Hc; Hm is left for the ordinary convergent shortening), scaled by
    the caller's band `intensity` and, more gently, by how fast the slab is converging.
    Nothing at all where the boundary is not closing, or where band intensity is zero."""
    hc = np.full(4, 25_000.0)
    hm = np.full(4, 55_000.0)
    closing = np.array([_closing_m_per_s(2.0), _closing_m_per_s(6.0), -_closing_m_per_s(4.0), _closing_m_per_s(6.0)])
    intensity = np.array([1.0, 1.0, 1.0, 0.0])

    new_hc, new_hm = rheology.apply_arc_magmatic_thickening(hc, hm, closing, years_myr=1.0, intensity=intensity)

    assert np.array_equal(new_hm, hm)  # Hm untouched here
    assert new_hc[0] > hc[0]
    assert new_hc[1] > new_hc[0]  # faster convergence -> more melt flux
    assert new_hc[2] == hc[2]  # diverging: no arc
    assert new_hc[3] == hc[3]  # outside the band (intensity 0): nothing
    # A sustained ~5 cm/yr margin at *full* band intensity for tens of Myr builds a thick
    # Altiplano-scale arc welt -- the band-average node (intensity ~0.5, slower convergence)
    # sees a few km, which is the geologically ordinary case.
    hc_over_40myr, _ = rheology.apply_arc_magmatic_thickening(
        np.array([25_000.0]), np.array([55_000.0]), np.array([_closing_m_per_s(5.0)]), years_myr=40.0, intensity=np.array([1.0])
    )
    assert 8_000.0 < hc_over_40myr[0] - 25_000.0 < 22_000.0


def test_arc_magmatism_flux_saturates_at_a_fast_margin():
    """The convergence multiplier is capped so an unusually fast slab can't flux unbounded
    melt -- a 30 cm/yr margin adds only modestly more than the reference 5 cm/yr one."""
    one = np.array([1.0])
    ref, _ = rheology.apply_arc_magmatic_thickening(
        np.array([25_000.0]), np.array([55_000.0]), np.array([_closing_m_per_s(5.0)]), years_myr=1.0, intensity=one
    )
    fast, _ = rheology.apply_arc_magmatic_thickening(
        np.array([25_000.0]), np.array([55_000.0]), np.array([_closing_m_per_s(30.0)]), years_myr=1.0, intensity=one
    )
    assert fast[0] > ref[0]
    assert (fast[0] - 25_000.0) <= rheology.ARC_MAGMATIC_CONVERGENCE_CAP * (ref[0] - 25_000.0) + 1e-6


def test_delamination_melt_intrusion_only_conserves_a_bounded_fraction_of_overflow():
    """GitHub issue #145's reopened investigation: the first version of #161's overflow
    conservation thrust the *entire* clipped-off Hc onto the near-field ring instantly, fully
    conserved, every step -- which turned out to over-thicken/over-elevate the majority of a
    real save's continental land within tens of Myr. `apply_delamination_melt_intrusion`
    replaces that with a bounded melt fraction at a bounded per-Myr rate, so a large overflow
    should add real but modest growth, well under the full overflow amount, when plenty of
    time (`years_myr`) is available."""
    hc_near = np.full(10, 40_000.0)
    overflow_hc = 100_000.0  # a large one-step overflow from the core band

    new_hc = rheology.apply_delamination_melt_intrusion(hc_near, overflow_hc, years_myr=1.0)

    added_total = float(np.sum(new_hc - hc_near))
    assert added_total > 0.0
    assert added_total < overflow_hc  # not fully conserved -- most of it delaminates for good
    assert np.all(new_hc <= lithosphere.MAX_CRUSTAL_THICKNESS_M)


def test_delamination_melt_intrusion_rate_limited_not_instant():
    """However much melt fraction is theoretically available, it can't all arrive in one
    vanishingly small step -- the whole point of routing it through a bounded per-Myr rate
    instead of #161's original instant full-overflow dump. A tiny `years_myr` should add far
    less than a long one, for the same overflow."""
    hc_near = np.full(5, 40_000.0)
    overflow_hc = 50_000.0

    tiny_step = rheology.apply_delamination_melt_intrusion(hc_near, overflow_hc, years_myr=1e-4)
    long_step = rheology.apply_delamination_melt_intrusion(hc_near, overflow_hc, years_myr=5.0)

    assert float(np.sum(tiny_step - hc_near)) < float(np.sum(long_step - hc_near))


def test_delamination_melt_intrusion_no_op_without_overflow_or_receivers():
    hc_near = np.full(5, 40_000.0)
    assert np.array_equal(rheology.apply_delamination_melt_intrusion(hc_near, 0.0, years_myr=5.0), hc_near)
    empty = np.array([])
    assert rheology.apply_delamination_melt_intrusion(empty, 50_000.0, years_myr=5.0).size == 0


def test_arc_magmatism_caps_hc_at_the_same_ceiling_as_convergent_thickening():
    """Arc underplating is a second, independent (additive rather than multiplicative, so far
    slower) source of unbounded Hc growth -- a sustained multi-hundred-Myr arc with no cap at
    all would still eventually climb past what real continental crust can hold before
    delaminating. A node already at the ceiling gets no further juvenile crust added."""
    hc = np.array([lithosphere.MAX_CRUSTAL_THICKNESS_M])
    hm = np.array([100_000.0])
    closing = np.array([_closing_m_per_s(6.0)])

    new_hc, new_hm = rheology.apply_arc_magmatic_thickening(hc, hm, closing, years_myr=1.0, intensity=np.array([1.0]))

    assert new_hc[0] == lithosphere.MAX_CRUSTAL_THICKNESS_M
    assert new_hm[0] == hm[0]


def test_extension_thins_crust_and_triggers_decompression_melting():
    """The divergent branch is the same yield check with the opposite sign -- a fast rift
    thins Hc, and a column dragged below `RIFT_CRITICAL_THICKNESS_M` flags for melting."""
    hc = np.array([5_500.0, 40_000.0])
    hm = np.array([30_000.0, 100_000.0])
    opening = np.array([-_closing_m_per_s(10.0), -_closing_m_per_s(10.0)])

    new_hc, _, melting = rheology.apply_divergent_deformation(hc, hm, opening, years_myr=10.0)

    assert np.all(new_hc < hc)
    assert bool(melting[0]) and not bool(melting[1])


def test_rift_magmatism_adds_hc_only_below_onset_while_extending():
    """Nodes below RIFT_VOLCANISM_ONSET_HC_M and genuinely extending get a partial Hc
    offset; nodes still above onset, or not extending, get nothing."""
    hc = np.array([12_000.0, 25_000.0, 12_000.0])
    hm = np.array([60_000.0, 80_000.0, 60_000.0])
    closing = np.array([-_closing_m_per_s(3.0), -_closing_m_per_s(3.0), _closing_m_per_s(3.0)])

    new_hc, new_hm = rheology.apply_rift_magmatic_thickening(hc, hm, closing, years_myr=1.0)

    assert np.array_equal(new_hm, hm)  # Hm untouched
    assert new_hc[0] > hc[0]  # below onset, extending -> partial offset
    assert new_hc[1] == hc[1]  # still above onset -> nothing
    assert new_hc[2] == hc[2]  # converging, not extending -> nothing


def test_rift_magmatism_ramps_up_as_crust_approaches_full_rupture():
    """"Further rifting -> more volcanism": a column close to RIFT_CRITICAL_THICKNESS_M gets
    a bigger offset than one that just crossed onset, at the same extension rate, and never
    overshoots back above onset in one step at an ordinary rate."""
    closing = np.array([-_closing_m_per_s(3.0)])
    near_onset, _ = rheology.apply_rift_magmatic_thickening(np.array([19_000.0]), np.array([60_000.0]), closing, years_myr=1.0)
    near_rupture, _ = rheology.apply_rift_magmatic_thickening(np.array([6_000.0]), np.array([60_000.0]), closing, years_myr=1.0)

    assert (near_onset[0] - 19_000.0) < (near_rupture[0] - 6_000.0)
    assert near_onset[0] < rheology.RIFT_VOLCANISM_ONSET_HC_M


def test_rift_magmatism_flux_saturates_at_a_fast_extension_rate():
    """The extension multiplier is capped, mirroring apply_arc_magmatic_thickening's own
    convergence cap -- an unusually fast rift can't flux unbounded melt."""
    ref, _ = rheology.apply_rift_magmatic_thickening(
        np.array([12_000.0]), np.array([60_000.0]), np.array([-_closing_m_per_s(2.0)]), years_myr=1.0
    )
    fast, _ = rheology.apply_rift_magmatic_thickening(
        np.array([12_000.0]), np.array([60_000.0]), np.array([-_closing_m_per_s(20.0)]), years_myr=1.0
    )
    assert fast[0] > ref[0]
    assert (fast[0] - 12_000.0) <= rheology.RIFT_MAGMATIC_EXTENSION_CAP * (ref[0] - 12_000.0) + 1e-6

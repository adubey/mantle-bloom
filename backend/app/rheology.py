"""Elastic-viscoplastic Mohr-Coulomb deformation (spec section 2.3): replaces v1's empirical
per-Myr rate tables (`plates.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR` and friends) with a strain-
rate-driven update to the crustal/mantle-lithosphere thickness columns (`Hc`/`Hm`), gated by
a real yield check. Elevation itself is never touched directly here -- see `lithosphere.py`'s
`sync_plate_elevation`, called by `lithosphere_plate.py` after every deform() pass, which derives it
from Hc/Hm via isostasy.

Scope: a per-node scalar-stress proxy at boundary-classified nodes, not a full 2D
depth-integrated stress tensor field -- see the plan's own "Scope and fidelity calls." Every
term in Eqs. 3-5 is computed, at the same boundary-local dimensional reduction v1's own
distance-decay effects already use.
"""

from __future__ import annotations

import numpy as np

from .elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M
from . import lithosphere

# Mohr-Coulomb yield criterion Y = C + sigma_n * tan(phi), Eq. 4 -- typical crustal values.
COHESION_PA = 2e7  # C, ~20 MPa, a plausible upper-crustal cohesion
INTERNAL_FRICTION_ANGLE_RAD = np.radians(30.0)  # phi, Byerlee's-law-consistent ballpark

# Converts a boundary-normal closing/opening rate (m/s) into a normal stress proxy sigma_n
# (Pa) -- not a literal depth-integrated rheology solve, just enough of a stand-in to make
# the yield check respond to how fast two plates are actually converging/diverging, not
# merely whether they geometrically touch.
#
# Calibration (2026-09): sigma_n = this * closing_rate_m_per_s, and the closing rate for a
# real continental collision is a few cm/yr -- 3 cm/yr is ~9.5e-10 m/s. The Mohr-Coulomb
# yield stress here is COHESION_PA (2e7) plus a friction term, so sigma_n has to reach
# ~5e7 Pa before `yield_excess` is even nonzero. The original 3e13 put sigma_n at ~3 cm/yr
# around 3e4 Pa -- three orders of magnitude below yield -- so `apply_convergent_deformation`
# returned Hc unchanged at *every* plate speed mantle.MAX_PLATE_RATE (15 cm/yr) allows: the
# engine never thickened crust, never built a mountain, and continents only ever thinned
# (rifting + erosion-isostasy) and drowned, with land fraction falling monotonically and the
# stalled multi-plate overlaps the `overlapAge` view shows never crumpling into orogens.
# 1e17 puts a sustained 3 cm/yr collision a few x past yield (plastic strain ~0.016/Myr, Hc
# doubles over ~45 Myr -- the Himalaya/Tibet timescale), while a ~1 cm/yr graze stays
# sub-yield (no spurious mountains) and >=5 cm/yr saturates.
EFFECTIVE_LITHOSPHERE_VISCOSITY_PA_S_PER_M = 1e17

SECONDS_PER_YEAR = 365.25 * 86400.0

# Eq. 5's mass-conservation transport, discretized: plastic normal strain rate converts
# directly into a fractional thickness change per Myr once yielded -- this is the "rate" a
# real viscoplastic flow law would otherwise derive from the stress excess over yield;
# calibrated so a sustained, fast collision (well over yield) still builds real mountains
# over tens of Myr, matching the real Himalaya/Tibetan Plateau timescale v1's own
# CONVERGENT_MOUNTAIN_RATE_M_PER_MYR was calibrated against.
PLASTIC_THICKENING_RATE_PER_MYR_PER_YIELD_EXCESS = 0.06

# Section 2.3: below this Hc, decompression melting erupts new oceanic crust at a rift --
# the spec's own literal ~5km critical-thinning trigger, replacing v1's flat
# STRETCH_VOLCANO_PROBABILITY roll.
RIFT_CRITICAL_THICKNESS_M = 5_000.0

# Same fold-thrust-belt "not every point in a collision belt rises at the same rate" texture
# v1 modelled via REVERSE_FAULT_VALLEY_UPLIFT_FACTOR -- reused here as a multiplier on the
# plastic strain rate itself (the physically-motivated cause: a downthrown fault block
# accumulates less shortening than the thrust sheets around it), not on elevation directly.
REVERSE_FAULT_VALLEY_UPLIFT_FACTOR = 0.15


def normal_closing_rate_m_per_s(plate_omega: np.ndarray, neighbor_omega: np.ndarray, points_xyz: np.ndarray, direction_to_neighbor: np.ndarray) -> np.ndarray:
    """`boundary.closing_rate`'s own formula (relative tangential velocity projected onto the
    boundary-normal direction), reused directly -- positive means converging. Returned in
    real m/s (see torque.py's own unit-convention docstring for why the *2/SECONDS_PER_YEAR
    conversion is needed: `plate.omega` is real rad/yr, `points_xyz` are unit vectors, so
    `omega x point` is numerically an omega-equivalent "rad/yr" tangential rate that becomes a
    real m/yr velocity once scaled by the planet's actual radius)."""
    v_self = np.cross(plate_omega, points_xyz)
    v_neighbor = np.cross(neighbor_omega, points_xyz)
    closing = np.sum((v_self - v_neighbor) * direction_to_neighbor, axis=-1)
    return closing * lithosphere.PLANET_RADIUS_M / SECONDS_PER_YEAR


def yield_excess(sigma_n_pa: np.ndarray) -> np.ndarray:
    """`sigma_e - Y` (Eq. 4), clipped to >= 0 -- the *elastic* regime (sigma_e < Y) leaves
    Hc/Hm untouched entirely; only nodes at or past yield accumulate plastic strain, and the
    magnitude past yield sets how fast."""
    yield_stress = COHESION_PA + np.abs(sigma_n_pa) * np.tan(INTERNAL_FRICTION_ANGLE_RAD)
    return np.clip(np.abs(sigma_n_pa) - yield_stress, 0.0, None)


def plastic_strain_rate_per_myr(closing_rate_m_per_s: np.ndarray) -> np.ndarray:
    """The fractional-thickness-change rate (positive = thickening under convergence,
    negative = thinning under extension) a node accumulates once past yield -- Eq. 3/4's
    plastic regime, discretized into a single scalar rate rather than a full strain-rate
    tensor (see module docstring's scope note)."""
    sigma_n = closing_rate_m_per_s * EFFECTIVE_LITHOSPHERE_VISCOSITY_PA_S_PER_M
    excess = yield_excess(sigma_n)
    sign = np.sign(closing_rate_m_per_s)
    # excess grows unboundedly with |closing_rate| under this linear-viscosity stand-in;
    # normalizing by the yield stress itself keeps the rate a well-conditioned O(1)-ish
    # multiplier on PLASTIC_THICKENING_RATE_PER_MYR_PER_YIELD_EXCESS across the whole range
    # of plate speeds mantle.MIN/MAX_PLATE_RATE allow, rather than blowing up at high speed.
    yield_stress = COHESION_PA + np.abs(sigma_n) * np.tan(INTERNAL_FRICTION_ANGLE_RAD)
    normalized_excess = np.where(yield_stress > 0, excess / yield_stress, 0.0)
    return sign * normalized_excess * PLASTIC_THICKENING_RATE_PER_MYR_PER_YIELD_EXCESS


def apply_convergent_deformation(
    hc_m: np.ndarray,
    hm_m: np.ndarray,
    closing_rate_m_per_s: np.ndarray,
    years_myr: float,
    fault_factor: np.ndarray,
    strength: np.ndarray | float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Contested (convergent) nodes: mass-conserving thickening under compression. Hc and Hm
    both thicken in proportion (crustal shortening drags the attached mantle lithosphere
    along with it) -- `fault_factor` (1.0 almost everywhere, `REVERSE_FAULT_VALLEY_UPLIFT_
    FACTOR` on noise-selected downthrown blocks, same pattern v1 used) modulates how much of
    the plastic strain this particular node actually accumulates, giving the same
    discrete-thrust-sheet visual texture v1 had, now as a real strain-rate multiplier rather
    than a post-hoc elevation multiplier. `strength` is the live collision-uplift tuning knob
    (World.collision_uplift_multiplier, plus the reach knob's near-field taper -- see
    lithosphere_plate.py); a plain 1.0 default keeps every existing caller/behaviour
    unchanged."""
    rate = plastic_strain_rate_per_myr(closing_rate_m_per_s)
    rate = np.clip(rate, 0.0, None)  # convergent branch only ever thickens
    fractional_change = rate * years_myr * fault_factor * strength
    new_hc = hc_m * (1.0 + fractional_change)
    new_hm = hm_m * (1.0 + fractional_change)
    return new_hc, new_hm


def apply_divergent_deformation(hc_m: np.ndarray, hm_m: np.ndarray, closing_rate_m_per_s: np.ndarray, years_myr: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Uncontested, extensional (opening) boundary nodes: crust thins under tension. Returns
    (new_hc, new_hm, decompression_melting_mask) -- the mask marks nodes whose Hc just
    crossed below `RIFT_CRITICAL_THICKNESS_M`, the spec's own decompression-melting trigger
    (Section 2.3), for the caller to spawn fresh oceanic crust there (see lithosphere_plate.py)."""
    rate = plastic_strain_rate_per_myr(closing_rate_m_per_s)
    rate = np.clip(rate, None, 0.0)  # divergent branch only ever thins
    fractional_change = rate * years_myr
    was_above = hc_m >= RIFT_CRITICAL_THICKNESS_M
    new_hc = np.clip(hc_m * (1.0 + fractional_change), lithosphere.MIN_CRUSTAL_THICKNESS_M, None)
    new_hm = np.clip(hm_m * (1.0 + fractional_change), lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M, None)
    melting = was_above & (new_hc < RIFT_CRITICAL_THICKNESS_M)
    return new_hc, new_hm, melting


# Continental arc magmatism (Cordilleran / Andean active margins). Distinct from
# `apply_convergent_deformation`'s yield-limited plastic *shortening*: subduction dehydrates
# the down-going slab, fluxes the mantle wedge, and the melt underplates / intrudes the
# overriding continental crust -- juvenile mass added from the mantle, not conserved from the
# neighbour, and it happens whether or not the margin is at Mohr-Coulomb yield. This is the
# crust-*building* half of "an oceanic plate subducting under a continent makes more
# continent" (the areal half is `lithosphere_plate.ARC_MARGIN_SEED_*`). It acts over the
# whole arc *band* inboard of the trench (`lithosphere_plate` passes a per-node distance-
# falloff `intensity`, not just the contact line -- the contact line alone is only a few tens
# of nodes, far too narrow to matter), with a gentle extra dependence on convergence rate
# (more slab -> more flux). Calibrated so a sustained ~5 cm/yr margin at full band intensity
# adds ~9-14 km of Hc over the tens of Myr an arc is active -- the order of measured Andean
# crustal-growth rates -- without runaway (the CONTINENTAL_AREA_BUDGET_MULT volume gate still
# bounds the plate's footprint).
ARC_MAGMATIC_HC_RATE_M_PER_MYR = 450.0
ARC_REFERENCE_CONVERGENCE_M_PER_S = 0.05 / SECONDS_PER_YEAR  # 5 cm/yr
ARC_MAGMATIC_CONVERGENCE_CAP = 3.0  # a very fast margin fluxes at most 3x the reference
ARC_MIN_CONVERGENCE_M_PER_S = 0.002 / SECONDS_PER_YEAR  # 0.2 cm/yr -- below this it's a graze, no arc


def apply_arc_magmatic_thickening(
    hc_m: np.ndarray, hm_m: np.ndarray, closing_rate_m_per_s: np.ndarray, years_myr: float, intensity: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Add juvenile arc crust to the overriding continental margin. `intensity` is the
    caller's per-node band weight (1 at the trench, fading inboard). Only Hc grows -- arc
    magmatism thickens the crustal column; the attached mantle lithosphere is returned
    unchanged (the caller still runs the ordinary convergent shortening on the contested
    subset, which does drag Hm along). Nodes not actually converging
    (`closing_rate <= ARC_MIN_CONVERGENCE_M_PER_S`) get nothing."""
    active = closing_rate_m_per_s > ARC_MIN_CONVERGENCE_M_PER_S
    convergence = np.clip(closing_rate_m_per_s / ARC_REFERENCE_CONVERGENCE_M_PER_S, 0.0, ARC_MAGMATIC_CONVERGENCE_CAP)
    rate_mult = np.where(active, np.clip(0.4 + 0.6 * convergence, 0.0, ARC_MAGMATIC_CONVERGENCE_CAP), 0.0)
    new_hc = hc_m + ARC_MAGMATIC_HC_RATE_M_PER_MYR * years_myr * rate_mult * np.asarray(intensity)
    return new_hc, hm_m


def relax_young_oceanic_mantle_lithosphere(hm_m: np.ndarray, divergent_age_myr: np.ndarray, years_myr: float) -> np.ndarray:
    """Freshly-formed ridge crust starts with thin mantle lithosphere (`lithosphere.
    YOUNG_RIDGE_HM_M`) and thickens toward the reference oceanic value as it cools and ages
    -- see `lithosphere.py`'s own note on why this (age-keyed relaxation toward one fixed
    reference) stands in for real open-ended sqrt(age) thickening. Reuses the exact same
    `divergent_age_myr` field/relax-toward-target shape v1's own elevation relaxation used."""
    target = lithosphere.REFERENCE_HM_OCEANIC_M
    relax_factor = 1.0 - np.exp(-0.1 * years_myr)  # same order of magnitude as v1's DIVERGENT_RELAX_RATE_PER_MYR
    still_young = divergent_age_myr < 30.0  # Myr -- oceanic lithosphere is largely equilibrated well before this
    new_hm = np.where(still_young, hm_m + (target - hm_m) * relax_factor, hm_m)
    return new_hm


def clip_elevation_bounds(z: np.ndarray) -> np.ndarray:
    return np.clip(z, MIN_ELEVATION_M, MAX_ELEVATION_M)

"""Elastic-viscoplastic Mohr-Coulomb deformation (spec section 2.3): replaces v1's empirical
per-Myr rate tables (a flat convergent uplift rate in m/Myr and friends) with a strain-
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
# the spec's own literal ~5km critical-thinning trigger, replacing v1's flat per-event
# stretch-volcano roll.
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


def convergent_strain(
    closing_rate_m_per_s: np.ndarray, years_myr: float, fault_factor: np.ndarray, strength: np.ndarray | float = 1.0
) -> np.ndarray:
    """The fractional thickening a convergent node takes on this step, before any cap --
    `apply_convergent_deformation`'s strain, and the shortening a quad plate's band demands of
    its interior (shortening.py, issue #314)."""
    rate = np.clip(plastic_strain_rate_per_myr(closing_rate_m_per_s), 0.0, None)  # convergent branch only ever thickens
    return rate * years_myr * fault_factor * strength


def apply_convergent_deformation(
    hc_m: np.ndarray,
    hm_m: np.ndarray,
    closing_rate_m_per_s: np.ndarray,
    years_myr: float,
    fault_factor: np.ndarray,
    strength: np.ndarray | float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Contested (convergent) nodes: mass-conserving thickening under compression. Hc and Hm
    both thicken in proportion (crustal shortening drags the attached mantle lithosphere
    along with it) -- `fault_factor` (1.0 almost everywhere, `REVERSE_FAULT_VALLEY_UPLIFT_
    FACTOR` on noise-selected downthrown blocks, same pattern v1 used) modulates how much of
    the plastic strain this particular node actually accumulates, giving the same
    discrete-thrust-sheet visual texture v1 had, now as a real strain-rate multiplier rather
    than a post-hoc elevation multiplier. `strength` is the live collision-uplift tuning knob
    (World.collision_uplift_multiplier -- see lithosphere_plate.boundary_context); a plain 1.0
    default keeps every existing caller/behaviour
    unchanged.

    Hc/Hm are clipped at `lithosphere.MAX_CRUSTAL_THICKNESS_M`/`MAX_MANTLE_LITHOSPHERE_
    THICKNESS_M` -- see those constants' own comment (GitHub issue #161: with no ceiling here,
    a node sitting in a long-lived convergent regime compounds this exponentially, run after
    run, with no physical floor on how tall/thick a single column can get). The third return
    value, `overflow_hc_m`, is how much Hc growth this step actually got clipped off (>= 0,
    zero everywhere the node wasn't already at the ceiling) -- real continental crust doesn't
    just vanish at that ceiling, it spreads laterally into the surrounding foreland (a
    fold-thrust belt widening once its hinterland can't thicken any further), so the caller
    (lithosphere_plate.deform_columns) is expected to place this rather than silently dropping
    it, the same mass-conserving idiom suture accretion uses for suture retreat. Hm's own overflow is not returned/conserved -- unlike buoyant crust, an
    over-thickened mantle-lithosphere root has nowhere to spread to; it delaminates (sinks into
    the asthenosphere), a real geodynamic sink, not a modeling shortcut."""
    fractional_change = convergent_strain(closing_rate_m_per_s, years_myr, fault_factor, strength)
    uncapped_new_hc = hc_m * (1.0 + fractional_change)
    new_hc = np.clip(uncapped_new_hc, None, lithosphere.MAX_CRUSTAL_THICKNESS_M)
    new_hm = np.clip(hm_m * (1.0 + fractional_change), None, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)
    overflow_hc_m = uncapped_new_hc - new_hc
    return new_hc, new_hm, overflow_hc_m


# Lateral magma transport (GitHub issue #205, follow-up to #120's "Land fraction slowly
# declines"). Everything that thickened crust then acted right at a collision boundary (or, at
# most, a fixed ~350km ring on the *same* plate) -- #120's
# own remaining land-loss driver is that nothing carries mass from a plate being over-thickened
# by collision to a distant, over-stretched interior losing land to thinning. This is the
# *source* half of the fix (the transport/deposit half lives in the new magma_transport.py,
# called from world.step_world): a small, continuous skim off the same yield-driven strain flux
# `apply_convergent_deformation` already computes, diverted to a mobile magma parcel instead of
# thickening the node in place. It's melt, not diverted solid rock -- a fraction of the strain
# that would otherwise pile up as coherent shortened crust instead partially melts at the point
# of generation (real over-thickened collision-belt crust really does partially anatectically
# melt), and melt is exactly the thing in this model that can travel (matching every other
# magmatic-addition path here -- arc, delamination, rift underplating -- already being Hc-only,
# no Hm coupling). Because this comes *out of* `apply_convergent_deformation`'s own strain
# increment (not an independent addition on top), it correctly shrinks how much reaches that
# function's existing ceiling-overflow path (`overflow_hc_m`) rather than being additive to it.
#
# Design discussion (GitHub issue #205) went through two full review rounds before landing here
# (v3): the earlier drafts skimmed the line engine's near-field ring too, which starved that
# ring's other melt supply -- so this only ever applies to the core `convergent` mask. Earlier
# drafts also skimmed every convergent node unconditionally; a node already close enough to
# `MAX_CRUSTAL_THICKNESS_M` to be generating (or about to generate) ceiling overflow needs that
# overflow undisturbed (the caller places it -- quad_tectonics._place_ceiling_overflow), so
# `MAGMA_EXPORT_HC_CEILING_FRACTION` exempts nodes already past that
# headroom threshold -- they get back their full, undiminished strength and contribute nothing
# to the export pool.
MAGMA_EXPORT_HC_CEILING_FRACTION = 0.9
# Placeholder pending an #120-style empirical toggle sweep (against land fraction *and* #146's
# own ELEV_CHANGE_COLLISION land-share diagnostic, per the issue's own agreed validation plan) --
# no first-principles default exists for this any more than for GRANITIC_MELT_FRACTION above.
MAGMA_EXPORT_FRACTION = 0.1


def magma_export_strength_and_volume(
    hc_m: np.ndarray,
    closing_rate_m_per_s: np.ndarray,
    years_myr: float,
    fault_factor: np.ndarray,
    strength: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """For the core convergent band only (see module comment above): splits `strength` into what stays in place (`reduced_strength`, to feed
    `apply_convergent_deformation` in the caller's place of the original `strength`) and how
    much Hc-equivalent volume per node is withheld into a mobile magma parcel instead
    (`export_hc_m`).

    Computed closed-form rather than by calling `apply_convergent_deformation` twice: `strength`
    is a pure multiplier on `fractional_change`, so the volume withheld by scaling it down by
    `MAGMA_EXPORT_FRACTION` is exactly `hc_m * full_fractional_change * MAGMA_EXPORT_FRACTION`,
    where `full_fractional_change` is the *unreduced* fractional thickening
    `apply_convergent_deformation` would otherwise apply (before any ceiling clip -- the clip
    only matters to the in-place path, not to how much strain was available to skim from).

    A node already at/past `MAGMA_EXPORT_HC_CEILING_FRACTION * MAX_CRUSTAL_THICKNESS_M` is left
    completely undisturbed (`reduced_strength == strength`, `export_hc_m == 0`) -- it reverts to
    ordinary full-strength behaviour so `apply_convergent_deformation`'s own ceiling-overflow
    path is unaffected by this mechanism entirely."""
    rate = np.clip(plastic_strain_rate_per_myr(closing_rate_m_per_s), 0.0, None)
    full_fractional_change = rate * years_myr * fault_factor * strength
    has_headroom = hc_m < MAGMA_EXPORT_HC_CEILING_FRACTION * lithosphere.MAX_CRUSTAL_THICKNESS_M
    export_fraction = np.where(has_headroom, MAGMA_EXPORT_FRACTION, 0.0)
    reduced_strength = strength * (1.0 - export_fraction)
    export_hc_m = hc_m * full_fractional_change * export_fraction
    return reduced_strength, export_hc_m


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
    (`closing_rate <= ARC_MIN_CONVERGENCE_M_PER_S`) get nothing.

    Clipped at `lithosphere.MAX_CRUSTAL_THICKNESS_M` -- same ceiling `apply_convergent_
    deformation` enforces (issue #161), since this is a second, independent, unbounded-over-
    enough-Myr source of Hc growth (additive rather than multiplicative, so far slower to run
    away, but a sustained multi-hundred-Myr arc would still climb past it with no cap at all).
    No overflow to conserve here, unlike the convergent path: this mass is juvenile, added
    fresh from the mantle wedge rather than shortened out of the node's own prior column, so
    simply not adding more past the ceiling loses nothing that existed a moment ago."""
    active = closing_rate_m_per_s > ARC_MIN_CONVERGENCE_M_PER_S
    convergence = np.clip(closing_rate_m_per_s / ARC_REFERENCE_CONVERGENCE_M_PER_S, 0.0, ARC_MAGMATIC_CONVERGENCE_CAP)
    rate_mult = np.where(active, np.clip(0.4 + 0.6 * convergence, 0.0, ARC_MAGMATIC_CONVERGENCE_CAP), 0.0)
    new_hc = hc_m + ARC_MAGMATIC_HC_RATE_M_PER_MYR * years_myr * rate_mult * np.asarray(intensity)
    return np.clip(new_hc, None, lithosphere.MAX_CRUSTAL_THICKNESS_M), hm_m


# Rift magmatic underplating: the "further rifting -> more volcanism" middle stage a real
# continental rift passes through well before full rupture (RIFT_CRITICAL_THICKNESS_M's hard
# melt-through reset in `lithosphere_plate._erupt_melted_nodes`). Once extension has thinned
# crust past this onset, upwelling asthenosphere starts partially melting and the melt
# intrudes/erupts into the extending column -- a partial offset to the ongoing plastic
# thinning (GitHub issue #120, "Land fraction slowly declines", "over-stretched interiors": a
# continent's interior currently thins from full reference Hc all the way to
# RIFT_CRITICAL_THICKNESS_M with zero magmatic counterweight, unlike the convergent side's
# apply_convergent_deformation + apply_arc_magmatic_thickening pairing). Rate deliberately far
# below ARC_MAGMATIC_HC_RATE_M_PER_MYR: real rift magmatism is dominated by unseen intrusion
# rather than eruption, and globally contributes less new continental crust than arc
# accretion -- this is a partial brake on the thinning, not a replacement for genuine rifting
# (a sustained rift must still be able to reach RIFT_CRITICAL_THICKNESS_M and rupture).
RIFT_VOLCANISM_ONSET_HC_M = 20_000.0
RIFT_MAGMATIC_HC_RATE_M_PER_MYR = 150.0
RIFT_MAGMATIC_REFERENCE_EXTENSION_M_PER_S = 0.02 / SECONDS_PER_YEAR  # 2 cm/yr
RIFT_MAGMATIC_EXTENSION_CAP = 2.0
RIFT_MAGMATIC_MIN_EXTENSION_M_PER_S = 0.002 / SECONDS_PER_YEAR  # 0.2 cm/yr -- below this, no melt


def apply_rift_magmatic_thickening(
    hc_m: np.ndarray, hm_m: np.ndarray, closing_rate_m_per_s: np.ndarray, years_myr: float
) -> tuple[np.ndarray, np.ndarray]:
    """Partial Hc offset for divergent nodes already thinner than RIFT_VOLCANISM_ONSET_HC_M --
    see module constants above. Ramps from 0 intensity at onset to full intensity as Hc
    approaches RIFT_CRITICAL_THICKNESS_M ("further rifting -> more volcanism," not a step
    function), with a gentle extra dependence on extension rate (faster stretching -> more
    decompression melt) -- the same shape apply_arc_magmatic_thickening uses for convergence
    rate. Only Hc grows; Hm (mantle lithosphere) is untouched, matching that function's own
    convention (underplating adds crustal material, not mantle lithosphere). Not yield-gated,
    also matching arc magmatism: decompression melting doesn't care whether the extensional
    stress is past Mohr-Coulomb yield, only how thin the column already is and how fast it's
    still extending. Caller is expected to skip nodes already flagged `melting` this step
    (rheology.apply_divergent_deformation's own return) -- those already got the full
    melt-through reset and don't need a partial offset on top of it."""
    extending = closing_rate_m_per_s < -RIFT_MAGMATIC_MIN_EXTENSION_M_PER_S
    below_onset = hc_m < RIFT_VOLCANISM_ONSET_HC_M
    active = extending & below_onset
    onset_span = RIFT_VOLCANISM_ONSET_HC_M - RIFT_CRITICAL_THICKNESS_M
    depth_fraction = np.clip((RIFT_VOLCANISM_ONSET_HC_M - hc_m) / onset_span, 0.0, 1.0)
    extension = np.clip(-closing_rate_m_per_s / RIFT_MAGMATIC_REFERENCE_EXTENSION_M_PER_S, 0.0, RIFT_MAGMATIC_EXTENSION_CAP)
    rate_mult = np.where(active, depth_fraction * (0.4 + 0.6 * extension), 0.0)
    new_hc = hc_m + RIFT_MAGMATIC_HC_RATE_M_PER_MYR * years_myr * rate_mult
    return new_hc, hm_m


def relax_oceanic_mantle_lithosphere(hm_m: np.ndarray, is_oceanic: np.ndarray, years_myr: float) -> np.ndarray:
    """Oceanic mantle lithosphere relaxes toward the reference oceanic value every step:
    freshly-formed ridge crust starts thin (`lithosphere.YOUNG_RIDGE_HM_M`) and thickens as it
    cools, and an oceanic column thickened or thinned elsewhere settles back the same way --
    see `lithosphere.py`'s own note on why one fixed reference stands in for real open-ended
    sqrt(age) thickening. Only `is_oceanic` nodes relax: a continental terrane riding on an
    oceanic plate keeps its own keel (issue #311).

    Not age-gated. There's no seafloor-age field to gate on -- `divergent_age_myr` counts only
    *continuous* divergence and resets to 0 the step a node leaves the divergent band, so the
    old `divergent_age_myr < 30` gate passed nearly every node anyway (issue #311)."""
    target = lithosphere.REFERENCE_HM_OCEANIC_M
    relax_factor = 1.0 - np.exp(-0.1 * years_myr)  # same order of magnitude as v1's DIVERGENT_RELAX_RATE_PER_MYR
    return np.where(is_oceanic, hm_m + (target - hm_m) * relax_factor, hm_m)


def clip_elevation_bounds(z: np.ndarray) -> np.ndarray:
    return np.clip(z, MIN_ELEVATION_M, MAX_ELEVATION_M)

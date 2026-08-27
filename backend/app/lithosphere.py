"""The lithospheric column state: Airy isostasy (spec section 2), and the mass/moment-of-
inertia integrals `torque.py` needs to solve for plate motion (spec section 3.1, Eq 6).

Every `ElevationLine` node carries two new fields (`crustal_thickness_m`/
`mantle_lithosphere_thickness_m`, see elevation_lines.py) instead of an independently-set
`elevation`. `elevation` becomes a *cache*: after any mutation to Hc/Hm,
`sync_line_elevation`/`sync_plate_elevation` below recompute it via `isostatic_elevation` and
write it back onto the line, so every module downstream (render_image.py, erosion.py,
hydrology.py, stats.py, ...) keeps reading `line.elevation` exactly as before, unaware it's
derived rather than primary state.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M, PLANET_RADIUS_KM

if TYPE_CHECKING:
    from .lithosphere_plate import LithospherePlate

PLANET_RADIUS_M = PLANET_RADIUS_KM * 1000.0
GRAVITY_M_S2 = 9.81

# Densities, Section 2.1 -- kg/m^3.
RHO_CONTINENTAL_CRUST = 2700.0
RHO_OCEANIC_CRUST = 2900.0
RHO_LITHOSPHERE_MANTLE = 3300.0
RHO_ASTHENOSPHERE = 3250.0
RHO_WATER = 1000.0

# Reference column thicknesses (meters) new crust is seeded with -- Section 2.2's own
# worked examples ("Hc ~= 7km" oceanic, "Hc >= 50km" orogens) anchor the oceanic/continental
# crustal values directly; Hm (not given numerically in the spec) uses typical real
# lithospheric-mantle thicknesses away from any thermal/age gradient we don't model (see
# below). ~7+60=67km total oceanic lithosphere and ~35+100=135km total continental
# lithosphere are both realistic real-Earth orders of magnitude.
REFERENCE_HC_CONTINENTAL_M = 35_000.0
REFERENCE_HM_CONTINENTAL_M = 100_000.0
REFERENCE_HC_OCEANIC_M = 7_000.0
REFERENCE_HM_OCEANIC_M = 60_000.0
# Freshly-formed ridge crust starts thin (both crust and underlying mantle lid) and thickens
# as it ages/cools -- see rheology.py's divergent-branch relaxation of Hm toward
# REFERENCE_HM_OCEANIC_M keyed off the same `divergent_age_myr` field v1 already tracks.
# This is a deliberate scope simplification: real oceanic lithosphere thickens with the
# square root of age indefinitely; here it relaxes toward one fixed reference thickness
# instead of an open-ended age-dependent curve, which is enough to give young ridge crust a
# real (thin, weak, more buoyant) column distinct from old, thick, strongly slab-pulling
# crust approaching a trench, without carrying an extra unbounded age field.
YOUNG_RIDGE_HM_M = 8_000.0

MIN_CRUSTAL_THICKNESS_M = 500.0  # never let Hc integrate through zero
MIN_MANTLE_LITHOSPHERE_THICKNESS_M = 2_000.0


def reference_thickness(crust_type: str) -> tuple[float, float]:
    """(Hc, Hm) new crust of this type starts at."""
    if crust_type == "continental":
        return REFERENCE_HC_CONTINENTAL_M, REFERENCE_HM_CONTINENTAL_M
    return REFERENCE_HC_OCEANIC_M, YOUNG_RIDGE_HM_M


def crust_density(crust_type: str) -> float:
    return RHO_CONTINENTAL_CRUST if crust_type == "continental" else RHO_OCEANIC_CRUST


#  Eq. 1/2, taken completely literally (elevation measured from the asthenosphere's own
# natural buoyancy level, zero offset), puts realistic reference columns wildly above sea
# level: Hc=35km/Hm=100km continental crust computes to +4.4km, and Hc=7km/Hm=60km oceanic
# crust to only -0.24km (barely-submerged shelf, not deep ocean). This isn't a bug in Eq.
# 1/2 -- it's the standard, expected property of a single-column isostasy formula with no
# reference datum: real geodynamic models always calibrate elevation against a reference
# column (e.g. a mid-ocean-ridge or standard continental column), never read it off as an
# absolute value the way Eq. 1/2 alone implies. ISOSTATIC_REFERENCE_OFFSET_M is exactly that
# calibration: a single constant shift chosen so the *continental* reference column
# (REFERENCE_HC/HM_CONTINENTAL_M) lands at v1's own BASE_CONTINENTAL_M-equivalent (exactly
# 200m). The *oceanic* reference column then lands around -6289m -- deeper than v1's own
# BASE_OCEANIC_M (-3800m), but still a realistic abyssal-plain-to-trench depth, and (more
# importantly for how this actually plays) a continent-to-ocean contrast of ~6.5km, the same
# order of magnitude as real Earth's own ~4.6km average continent/ocean elevation difference.
ISOSTATIC_REFERENCE_OFFSET_M = -4184.615384615388


def isostatic_elevation(hc_m: np.ndarray, hm_m: np.ndarray, rho_c: float) -> np.ndarray:
    """Airy isostatic elevation/bathymetry, Eq. 1/2 plus `ISOSTATIC_REFERENCE_OFFSET_M`
    (see above). `rho_c` is a scalar (this plate's own crust density) since a plate's crust
    type doesn't vary node-to-node.

    The calibration offset is folded into `shifted_bracket` *before* the water-loading
    branch, not added on afterward to each branch separately -- adding it after applying
    Eq. 2's rescaling factor to the *unshifted* bracket would evaluate the two branches to
    different values at the same shifted_bracket == 0 crossing (a real discontinuity, and
    -- since the water-loaded branch's slope is steeper -- a non-monotonic kink right at sea
    level, confirmed directly as a real bug: thicker crust reading as *lower* elevation right
    across that boundary). Folding the offset in first keeps both branches equal to zero at
    the same crossing point by construction, so the whole piecewise function is continuous
    and strictly increasing in Hc everywhere, not just within either branch alone."""
    bracket = hc_m * (1.0 - rho_c / RHO_ASTHENOSPHERE) + hm_m * (1.0 - RHO_LITHOSPHERE_MANTLE / RHO_ASTHENOSPHERE)
    shifted_bracket = bracket + ISOSTATIC_REFERENCE_OFFSET_M
    water_loaded = shifted_bracket * (RHO_ASTHENOSPHERE / (RHO_ASTHENOSPHERE - RHO_WATER))
    z = np.where(shifted_bracket <= 0.0, water_loaded, shifted_bracket)
    return np.clip(z, MIN_ELEVATION_M, MAX_ELEVATION_M)


def sync_line_elevation(line, rho_c: float):
    """Recompute `line.elevation` from its current Hc/Hm columns -- call after any mutation
    to `crustal_thickness_m`/`mantle_lithosphere_thickness_m`. Returns a new `ElevationLine`
    (this module never mutates a line's arrays in place)."""
    z = isostatic_elevation(line.crustal_thickness_m, line.mantle_lithosphere_thickness_m, rho_c)
    return line.replace(elevation=z)


def sync_plate_elevation(plate: "LithospherePlate") -> None:
    """`sync_line_elevation` over every line on `plate`, written back via `set_lines` --
    the one call site that makes `elevation` track Hc/Hm after a batch of lines change."""
    rho_c = crust_density(plate.crust_type)
    plate.set_lines([sync_line_elevation(line, rho_c) for line in plate.lines])


def node_area_m2(spacing_rad: float) -> float:
    """Physical footprint area (m^2) of one lattice node at this world's line spacing.

    Nodes are deliberately placed so this is (almost exactly) constant across latitude: rows
    are equally spaced in phi (physically equidistant on a sphere regardless of latitude),
    and each row's own theta spacing is widened by 1/cos(phi) specifically so the node's
    *physical* along-row spacing is also `spacing_rad` (see elevation_lines.py's module
    docstring / docs/simulation-model.md#plate-local-frames) -- so every node, at any
    latitude, represents the same physical patch, `spacing_rad^2` steradians times R^2."""
    return (spacing_rad * PLANET_RADIUS_M) ** 2


def moment_of_inertia_tensor(points_xyz: np.ndarray, hc_m: np.ndarray, hm_m: np.ndarray, rho_c: float, spacing_rad: float) -> np.ndarray:
    """Eq. 6's mass moment of inertia tensor, discretized as a sum over nodes rather than a
    continuum integral: I_p = sum_i sigma_i * A_i * (R^2 * I_3x3 - x_i x_i^T), where
    `sigma_i = Hc_i*rho_c + Hm_i*rho_lith_mantle` is each node's own areal mass density
    (kg/m^2) and `A_i` its physical area (`node_area_m2`).

    The spec's Eq. 6 literally integrates a bare `rho_c(phi, theta)` over dA, which is
    dimensionally inconsistent (a volumetric density integrated over an area gives mass per
    unit length, not mass) -- the physically sensible reading, used here, is the full areal
    mass density of the moving lithospheric column (crust *and* the rigid mantle lid riding
    with it), not the crust alone."""
    if len(points_xyz) == 0:
        return np.zeros((3, 3))
    sigma = hc_m * rho_c + hm_m * RHO_LITHOSPHERE_MANTLE
    mass = sigma * node_area_m2(spacing_rad)
    r2 = PLANET_RADIUS_M**2
    eye = np.eye(3)
    # points_xyz are unit vectors (world xyz on the unit sphere); x_i x_i^T for the *physical*
    # position (R*point) is R^2 * (point point^T), matching r2*I3x3 - R^2*(point point^T).
    outer_unit = np.einsum("ni,nj->nij", points_xyz, points_xyz)
    per_node = mass[:, None, None] * r2 * (eye[None, :, :] - outer_unit)
    return per_node.sum(axis=0)


def angular_momentum(inertia_tensor: np.ndarray, omega: np.ndarray) -> np.ndarray:
    return inertia_tensor @ omega


def omega_from_angular_momentum(inertia_tensor: np.ndarray, angular_momentum_vec: np.ndarray) -> np.ndarray:
    """Inverts `angular_momentum` -- used both by torque integration (alpha = I^-1 tau) and
    by merge blending (conserve L, not omega, across a fusion -- see torque.merge_omega)."""
    try:
        return np.linalg.solve(inertia_tensor, angular_momentum_vec)
    except np.linalg.LinAlgError:
        # A plate with too few/degenerate nodes (near-empty, or all nodes collinear) can give
        # a singular tensor -- fall back to a damped least-squares solve rather than crashing
        # a whole step_world call over one edge-case plate.
        inertia_reg = inertia_tensor + 1e-6 * np.trace(inertia_tensor) * np.eye(3)
        return np.linalg.solve(inertia_reg, angular_momentum_vec)

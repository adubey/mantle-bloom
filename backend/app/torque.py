"""Dynamic torque-balance plate motion (spec section 3): slab-pull, ridge-push, basal drag,
and collision friction, integrated into an angular acceleration rather than fit via OLS to a
static mantle-flow field the way `mantle.fit_euler_pole` does for v1.

Unit convention, matched to the rest of this codebase (`mantle.py`'s own docstrings/
constants): `Plate.omega` is a true angular velocity in real radians/year -- radius-
independent, the same physical quantity `geometry.rotation_matrix_from_omega(omega, years)`
consumes directly, and what `mantle.flow_at` already returns (a field of tangential
"velocity" vectors evaluated at *unit*-sphere points, which are numerically identical to an
angular velocity since v = omega x r reduces to v = omega x p_unit when r = 1). The physics
in this module (slab-pull/ridge-push/basal-drag forces, moment of inertia) is expressed in
real SI units (N, kg, m, s) since that's what the force formulas in the spec are written in
-- `_real_velocity_m_per_s`/`SECONDS_PER_YEAR` are the only two places that cross between the
two unit systems.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import geometry, mantle
from .boundary import TRANSFORM_RATE_THRESHOLD
from .elevation_lines import line_spacing_rad
from . import lithosphere
from .plates import query_workers

SECONDS_PER_YEAR = 365.25 * 86400.0

# Section 3.2's own reach constants -- slab-pull/ridge-push are boundary-line integrals, so
# only nodes within a boundary-local band contribute; basal drag is a body force and applies
# to every node regardless of proximity to a boundary.
SUBDUCTION_LSINK_M = 250_000.0  # effective down-dip sunk-slab length contributing pull today
ASTHENOSPHERE_VISCOSITY_PA_S = 5e19  # within the spec's 1e19-1e21 Pa*s range
SHEAR_ZONE_THICKNESS_M = 100_000.0  # d_s, Eq. 10
# Deep-mantle viscosity the *sunk slab* shears against as it descends (`slab_drag_torque`),
# ~20x the asthenosphere value the surface plate's own base slides over -- the mantle below
# the asthenosphere is far stiffer. The spec's Eq. 8 gives slab *pull* but no matching
# resistance, so on its own it drives every subducting oceanic plate straight past
# MAX_PLATE_RATE and `clamp_rate` then pins them all there -- the "every oceanic plate railed
# at MAX" reading the overlapAge / long-run plate-geometry investigations kept seeing. A real
# slab's descent is resisted mostly by exactly this viscous coupling to the surrounding
# mantle, not by basal drag on the trailing surface plate; adding it back lets oceanic plate
# speed self-regulate (the faster a plate converges, the harder the mantle resists its slab)
# rather than everything sitting on the clamp. Within the spec's own 1e19-1e21 Pa*s bracket.
SLAB_MANTLE_VISCOSITY_PA_S = 1e21
# Typical old-ocean-floor depth ridge relief is measured against -- the isostatic elevation
# of oceanic crust already at its own reference thickness (lithosphere.REFERENCE_HC/HM_
# OCEANIC_M), so this tracks ISOSTATIC_REFERENCE_OFFSET_M automatically rather than
# hardcoding a second, driftable copy of the same number.
ABYSSAL_PLAIN_REFERENCE_ELEVATION_M = float(
    lithosphere.isostatic_elevation(
        np.array([lithosphere.REFERENCE_HC_OCEANIC_M]), np.array([lithosphere.REFERENCE_HM_OCEANIC_M]), lithosphere.RHO_OCEANIC_CRUST
    )[0]
)

# Collision friction: not a literal stress solve (see rheology.py), a reference stress scale
# applied to the [0, 1] collision intensity `rheology.py` already reports, opposing local
# relative motion at continent-continent contested nodes -- 50 MPa is a plausible
# lithospheric differential-stress magnitude.
COLLISION_FRICTION_REFERENCE_PA = 5e7

# Boundary-line integrals (slab-pull/ridge-push) treat each contributing node as owning one
# `spacing_rad * PLANET_RADIUS_M`-long stretch of the boundary -- consistent with how deform()
# already treats a line's own node spacing as the physical along-boundary resolution.


def _real_velocity_m_per_s(omega_equivalent_vectors: np.ndarray) -> np.ndarray:
    """Convert an (N,3) array of `omega`-unit ("rad/yr", radius-independent) tangential
    vectors -- e.g. `omega x unit_point` or `mantle.flow_at`'s own output -- into real
    physical velocity (m/s) at the planet's actual radius."""
    return omega_equivalent_vectors * lithosphere.PLANET_RADIUS_M / SECONDS_PER_YEAR


@dataclass
class BoundaryForceInputs:
    """One plate's own near-boundary geometry, gathered once per `compute_torque` call and
    reused across the slab-pull/ridge-push/collision terms -- avoids each force needing its
    own separate neighbour/cKDTree pass."""

    own_points: np.ndarray  # (N, 3) unit vectors, this plate's own nodes
    own_hc: np.ndarray
    own_hm: np.ndarray
    dist_to_neighbor: np.ndarray  # (N,) angular distance (rad) to nearest other-plate node, +inf if none
    direction_to_neighbor: np.ndarray  # (N, 3) unit vector from own node toward that nearest neighbour point
    neighbor_is_oceanic: np.ndarray  # (N,) bool
    neighbor_omega: np.ndarray  # (N, 3) the owning neighbour plate's own omega, for relative-velocity terms


def gather_boundary_force_inputs(plate, neighbours: list, spacing_rad: float, reach_rad: float) -> BoundaryForceInputs:
    own_points, _ = plate.all_points_and_elevation()
    own_hc = plate.collect("crustal_thickness_m")
    own_hm = plate.collect("mantle_lithosphere_thickness_m")
    n = len(own_points)
    if n == 0:
        empty3 = np.zeros((0, 3))
        return BoundaryForceInputs(empty3, np.zeros(0), np.zeros(0), np.full(0, np.inf), empty3, np.zeros(0, dtype=bool), empty3)

    no_neighbour = BoundaryForceInputs(
        own_points, own_hc, own_hm, np.full(n, np.inf), np.zeros((n, 3)), np.zeros(n, dtype=bool), np.zeros((n, 3))
    )
    if not neighbours:
        return no_neighbour

    # Query each neighbour's own cached node k-d tree (Plate.get_node_kdtree, invalidated on
    # rotate / node-set change exactly like the bounding polygon) and keep the elementwise
    # nearest, rather than concatenating every neighbour's node cloud into a fresh cKDTree on
    # every call. The same plate is a neighbour of several others and is queried in both the
    # shift and deform pass, so one cached tree per plate is built once per step instead of
    # the ~24 fresh builds this used to do. The argmin over neighbours reproduces the single
    # global-nearest the combined tree returned -- `neighbours` is in the same order the old
    # concatenation used, so distance ties break identically. scipy's compiled batch query
    # beats a pure-Python spherical-BVH recursion badly at this node count (tens of thousands
    # of points per plate at real density), so there's no accelerator worth swapping in here.
    workers = query_workers(n)
    best_dist = np.full(n, np.inf)
    best_point = np.zeros((n, 3))
    best_owner = np.full(n, -1)
    for i, neighbour in enumerate(neighbours):
        tree = neighbour.get_node_kdtree()
        if tree is None:
            continue
        dist, idx = tree.query(own_points, workers=workers)
        closer = dist < best_dist
        best_dist[closer] = dist[closer]
        best_point[closer] = tree.data[idx[closer]]
        best_owner[closer] = i

    if not np.any(best_owner >= 0):
        return no_neighbour

    direction = geometry.normalize(best_point - own_points)
    neighbor_is_oceanic = np.zeros(n, dtype=bool)
    neighbor_omega = np.zeros((n, 3))
    for i, neighbour in enumerate(neighbours):
        owned = best_owner == i
        if not np.any(owned):
            continue
        neighbor_is_oceanic[owned] = neighbour.crust_type == "oceanic"
        neighbor_omega[owned] = neighbour.omega
    dist = np.where(best_dist <= reach_rad, best_dist, np.inf)
    return BoundaryForceInputs(own_points, own_hc, own_hm, dist, direction, neighbor_is_oceanic, neighbor_omega)


def subducting_boundary_mask(plate, inputs: BoundaryForceInputs, reach_rad: float) -> np.ndarray:
    """This plate's own boundary-band nodes that are genuinely subducting: the plate is
    oceanic (only oceanic lithosphere is dense enough to descend, the same asymmetry
    `deform()`'s trench-vs-uplift branch draws), the node is within `reach_rad` of a
    neighbour, *and* the two plates are actually converging there (`closing_rate` past the
    same `TRANSFORM_RATE_THRESHOLD` `merge_split` uses). The convergence check is what this
    used to be missing: slab pull was applied to the plate's *entire* near-neighbour band,
    ridge and transform stretches included, so a plate spreading at three of its four edges
    still felt a full-perimeter pull and railed at MAX_PLATE_RATE. A ridge does not pull a
    slab down."""
    n = len(inputs.own_points)
    if plate.crust_type != "oceanic" or n == 0:
        return np.zeros(n, dtype=bool)
    band = inputs.dist_to_neighbor <= reach_rad
    v_self = np.cross(plate.omega, inputs.own_points)
    v_neighbor = np.cross(inputs.neighbor_omega, inputs.own_points)
    closing = np.sum((v_self - v_neighbor) * inputs.direction_to_neighbor, axis=-1)
    return band & (closing > TRANSFORM_RATE_THRESHOLD)


def slab_pull_torque(plate, inputs: BoundaryForceInputs, subducting_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """Eq. 8, discretized as a sum over this plate's own subducting-boundary nodes
    (`subducting_mask`, from `subducting_boundary_mask` -- passed in rather than recomputed so
    it stays in lockstep with the matching `slab_drag_torque` resistance)."""
    if not np.any(subducting_mask):
        return np.zeros(3)
    ds = spacing_rad * lithosphere.PLANET_RADIUS_M
    hm = inputs.own_hm[subducting_mask]
    force_mag = (lithosphere.RHO_LITHOSPHERE_MANTLE - lithosphere.RHO_ASTHENOSPHERE) * lithosphere.GRAVITY_M_S2 * hm * SUBDUCTION_LSINK_M * ds
    force = force_mag[:, None] * inputs.direction_to_neighbor[subducting_mask]
    r = inputs.own_points[subducting_mask] * lithosphere.PLANET_RADIUS_M
    return np.cross(r, force).sum(axis=0)


def slab_drag_coefficient_matrix(inputs: BoundaryForceInputs, subducting_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """The `K` half of Eq. 10-style drag, but for the *sunk slab* shearing against the deep
    mantle (`SLAB_MANTLE_VISCOSITY_PA_S`) over its `SUBDUCTION_LSINK_M` down-dip length,
    summed over this plate's subducting nodes. Folded into `integrate_omega`'s implicit `K`
    exactly like `basal_drag_coefficients` -- it is at least as stiff, and an explicit step on
    it would overshoot the same way. `b` is zero: the deep mantle it shears against is taken
    at rest (unlike the asthenosphere, which carries `mantle.flow_at`), so the drag is purely
    `-K @ omega`. Without this term slab pull is unopposed and every subducting oceanic plate
    rails at `MAX_PLATE_RATE` (see `SLAB_MANTLE_VISCOSITY_PA_S`)."""
    pts = inputs.own_points[subducting_mask]
    if len(pts) == 0:
        return np.zeros((3, 3))
    ds = spacing_rad * lithosphere.PLANET_RADIUS_M
    c = (
        (SLAB_MANTLE_VISCOSITY_PA_S / SHEAR_ZONE_THICKNESS_M)
        * (SUBDUCTION_LSINK_M * ds)
        * lithosphere.PLANET_RADIUS_M**2
        / SECONDS_PER_YEAR
    )
    return c * (len(pts) * np.eye(3) - np.einsum("ni,nj->nij", pts, pts).sum(axis=0))


def slab_drag_torque(plate, inputs: BoundaryForceInputs, subducting_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """`slab_drag_coefficient_matrix` evaluated at the plate's current `omega` (`-K @ omega`)
    -- a standalone entry point for tests / callers that just want the resistive torque;
    `shift_plate` itself folds the `K` into `integrate_omega` for implicit stability."""
    return -slab_drag_coefficient_matrix(inputs, subducting_mask, spacing_rad) @ plate.omega


def ridge_push_torque(plate, inputs: BoundaryForceInputs, divergent_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """Eq. 9, over this plate's own divergent/ridge boundary nodes (`divergent_mask`, the
    same geometric classification `rheology.py`'s deform() pass computes -- passed in rather
    than recomputed so the two stay in lockstep). Push direction is away from the boundary,
    into this plate's own interior -- the opposite sense from slab-pull's "toward the
    neighbour," since a ridge shoves both flanking plates apart."""
    if not np.any(divergent_mask):
        return np.zeros(3)
    ds = spacing_rad * lithosphere.PLANET_RADIUS_M
    z = lithosphere.isostatic_elevation(inputs.own_hc[divergent_mask], inputs.own_hm[divergent_mask], lithosphere.crust_density(plate.crust_type))
    e_r = np.clip(z - ABYSSAL_PLAIN_REFERENCE_ELEVATION_M, 0.0, None)
    force_mag = 0.5 * lithosphere.GRAVITY_M_S2 * (lithosphere.RHO_ASTHENOSPHERE - lithosphere.RHO_WATER) * e_r**2 * ds
    direction = -inputs.direction_to_neighbor[divergent_mask]
    # No neighbour within reach (an isolated ridge segment, direction left at zero) -- fall
    # back to the node's own outward radial tangent component is unnecessary here since
    # push-direction ambiguity without any neighbour reference just contributes no torque
    # (force magnitude still nonzero, but a zero direction vector cancels it) -- acceptable,
    # since "no neighbour nearby" also means this node isn't really acting as part of an
    # active spreading boundary the plate as a whole would feel a coherent push from.
    force = force_mag[:, None] * direction
    r = inputs.own_points[divergent_mask] * lithosphere.PLANET_RADIUS_M
    return np.cross(r, force).sum(axis=0)


def basal_drag_coefficients(plate, world, spacing_rad: float) -> tuple[np.ndarray, np.ndarray]:
    """Eq. 10's basal drag written as its exact affine-in-`omega` form, `tau(omega) = b - K @ omega`,
    over *every* node this plate owns (a body force, not a boundary integral).

    Each node's drag force is `F_i = (mu / d_s) * (v_mantle_i - v_plate_i) * A_i` with
    `v_plate_i = _real_velocity_m_per_s(omega x p_i)`, so the whole-plate drag torque
    `sum_i (R p_i) x F_i` is affine in `omega`:

        K = c * sum_i (I3 - p_i p_i^T)      -- symmetric, positive-semidefinite
        b = c * sum_i p_i x flow_i          -- `flow_i` = mantle.flow_at, this codebase's rad/yr units
        c = (mu / d_s) * A * R^2 / SECONDS_PER_YEAR

    (`c`'s `R^2 / SECONDS_PER_YEAR` folds in the two `_real_velocity_m_per_s` conversions plus
    the `R` lever arm; `sum_i p_i x (omega x p_i) == (sum_i I3 - p_i p_i^T) @ omega`.)

    Split out from the plain torque so `integrate_omega` can treat this term *implicitly*: it
    is by far the stiffest in the balance -- the asthenosphere coupling relaxes a plate toward
    the local mantle-flow rate in a small fraction of one tectonic step, so an explicit Euler
    step on it overshoots by ~19 orders of magnitude and `clamp_rate` then pins every plate at
    `MAX_PLATE_RATE` on step one, forever (the historical "every plate railed at MAX" bug)."""
    own_points, _ = plate.all_points_and_elevation()
    if len(own_points) == 0:
        return np.zeros(3), np.zeros((3, 3))
    c = (
        (ASTHENOSPHERE_VISCOSITY_PA_S / SHEAR_ZONE_THICKNESS_M)
        * lithosphere.node_area_m2(spacing_rad)
        * lithosphere.PLANET_RADIUS_M**2
        / SECONDS_PER_YEAR
    )
    flow = mantle.flow_at(own_points, world.mantle_centers)
    b = c * np.cross(own_points, flow).sum(axis=0)
    k = c * (len(own_points) * np.eye(3) - np.einsum("ni,nj->nij", own_points, own_points).sum(axis=0))
    return b, k


def basal_drag_torque(plate, world, spacing_rad: float) -> np.ndarray:
    """Eq. 10's basal drag evaluated at the plate's current `omega` -- `b - K @ omega` from
    `basal_drag_coefficients`. A standalone entry point for tests / callers that just want the
    torque; `shift_plate` itself uses the split `(b, K)` form for implicit integration."""
    b, k = basal_drag_coefficients(plate, world, spacing_rad)
    return b - k @ plate.omega


def collision_friction_torque(plate, inputs: BoundaryForceInputs, collision_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """A resistive torque at continent-continent contested nodes (`collision_mask`),
    proportional to `COLLISION_FRICTION_REFERENCE_PA` and opposing this plate's own local
    velocity relative to the colliding neighbour there -- keeps two head-on continents from
    accelerating straight through each other indefinitely."""
    if not np.any(collision_mask):
        return np.zeros(3)
    own_points = inputs.own_points[collision_mask]
    v_self = _real_velocity_m_per_s(np.cross(plate.omega, own_points))
    v_neighbor = _real_velocity_m_per_s(np.cross(inputs.neighbor_omega[collision_mask], own_points))
    relative = v_self - v_neighbor
    relative_dir = geometry.normalize(relative)
    speed = np.linalg.norm(relative, axis=-1)
    # Resistive stress scales with how fast the two plates are actually converging here (no
    # friction to overcome if they're not moving relative to each other), capped so a fast
    # collision doesn't blow up the resistive force past the reference stress scale itself.
    force = -COLLISION_FRICTION_REFERENCE_PA * lithosphere.node_area_m2(spacing_rad) * np.minimum(speed, 1.0)[:, None] * relative_dir
    r = own_points * lithosphere.PLANET_RADIUS_M
    return np.cross(r, force).sum(axis=0)


def integrate_omega(
    plate,
    explicit_torque: np.ndarray,
    drag_b: np.ndarray,
    drag_k: np.ndarray,
    inertia_tensor: np.ndarray,
    years: float,
) -> np.ndarray:
    """One angular-velocity step, implicit in the stiff basal-drag term and explicit in the
    rest:

        (I + g K) omega_new = I omega_old + g (tau_explicit + b)        g = years * SECONDS_PER_YEAR**2

    `(drag_b, drag_k)` is `basal_drag_coefficients`' affine split of Eq. 10's drag,
    `tau(omega) = b - K @ omega`; `explicit_torque` is every other (bounded, geometry-driven)
    torque. Backward Euler on the drag term is unconditionally stable for any `years` -- an
    explicit step on it does not converge at any step size this simulation uses (see
    `basal_drag_coefficients`), and the implicit drag term also dominates the system strongly
    enough to damp the explicit collision-friction term riding along in `explicit_torque`.
    Reduces to the old explicit `omega_old + g I^-1 tau` when `K = 0, b = 0`. The final
    `mantle.clamp_rate` (same physical speed bounds v1 enforced) is unchanged."""
    g = years * SECONDS_PER_YEAR**2
    lhs = inertia_tensor + g * drag_k
    rhs = inertia_tensor @ plate.omega + g * (explicit_torque + drag_b)
    try:
        new_omega = np.linalg.solve(lhs, rhs)
    except np.linalg.LinAlgError:
        # Near-empty / degenerate-footprint plate -> singular system; same damped fallback
        # `lithosphere.omega_from_angular_momentum` uses rather than crash the whole step.
        lhs_reg = lhs + 1e-6 * np.trace(lhs) * np.eye(3)
        new_omega = np.linalg.solve(lhs_reg, rhs)
    return mantle.clamp_rate(new_omega)


BOUNDARY_FORCE_REACH_MULTIPLIER = 3.0  # multiples of spacing_rad -- how boundary-local slab-pull/ridge-push/collision are


def classify_boundary_nodes(plate, neighbours: list, inputs: BoundaryForceInputs, reach_rad: float) -> tuple[np.ndarray, np.ndarray]:
    """(contested, divergent) for this plate's own nodes -- the same geometric test
    `rheology.py`'s deform() pass uses (a node is contested if it currently falls inside some
    neighbour's polygon), computed independently here since `shift()` runs *before*
    `deform()` each step and needs its own boundary read of the *pre-rotation* configuration
    to decide this step's driving forces."""
    n = len(inputs.own_points)
    if n == 0 or not neighbours:
        return np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    own_points = inputs.own_points

    # `near` (within reach_rad of a neighbour node) is the ordinary boundary-local band.
    # But a plate that has slid *deep* over a neighbour has its deep-interior overlapping
    # nodes far from any neighbour node, so they'd never get the polygon test and would sit
    # in neither `contested` nor `divergent` -- no thickening, no subduction, the overlap
    # just persists (reported: "plates 21 and 0 have quite some overlapping" -- 15% of one
    # plate's nodes on top of the other, static). Also test any node inside a neighbour's
    # bounding sphere (cheap triangle-inequality prefilter, same one
    # merge_split.find_continental_collision_pairs uses); the polygon `contains_batch` below
    # still only actually runs on nodes genuinely inside a sphere, so this costs ~one
    # angular_distance per neighbour when there's no real overlap. A deep continental
    # overlap then classifies contested -> rheology thickens Hc/Hm -> mountain uplift; a
    # deep oceanic overlap classifies contested -> subduction deletion -> the overlap heals.
    consider = inputs.dist_to_neighbor <= reach_rad
    for neighbour in neighbours:
        npts, _ = neighbour.all_points_and_elevation()
        if len(npts) == 0:
            continue
        centroid, radius = geometry.bounding_sphere(npts)
        consider |= geometry.angular_distance(own_points, centroid) <= radius

    contested = np.zeros(n, dtype=bool)
    if np.any(consider):
        consider_points = own_points[consider]
        consider_contested = np.zeros(len(consider_points), dtype=bool)
        for neighbour in neighbours:
            consider_contested |= neighbour.contains_batch(consider_points)
            if np.all(consider_contested):
                break
        contested[consider] = consider_contested
    # Only the boundary-local band (not the deep-overlap extension above) can be divergent --
    # an uncontested deep-interior node is just interior, not an opening rift.
    divergent = (inputs.dist_to_neighbor <= reach_rad) & ~contested
    return contested, divergent


def shift_plate(plate, world, other_plates: list, years: float) -> float:
    """`LithospherePlate.shift`'s real implementation: gather this step's driving/resisting
    torques from the plate's *current* (pre-rotation) boundary configuration, integrate omega,
    then rotate rigidly -- same `D` (max node displacement) contract as `PlateWithLines.shift`."""
    spacing_rad = line_spacing_rad(world.node_density)
    reach_rad = BOUNDARY_FORCE_REACH_MULTIPLIER * spacing_rad

    old_points, _ = plate.all_points_and_elevation()
    if len(old_points) == 0:
        return 0.0

    neighbours = plate.get_neighbours(other_plates)
    inputs = gather_boundary_force_inputs(plate, neighbours, spacing_rad, reach_rad)
    contested, divergent = classify_boundary_nodes(plate, neighbours, inputs, reach_rad)
    collision_mask = contested & (plate.crust_type == "continental") & ~inputs.neighbor_is_oceanic
    subducting = subducting_boundary_mask(plate, inputs, reach_rad)

    rho_c = lithosphere.crust_density(plate.crust_type)
    inertia = lithosphere.moment_of_inertia_tensor(inputs.own_points, inputs.own_hc, inputs.own_hm, rho_c, spacing_rad)

    explicit_torque = (
        slab_pull_torque(plate, inputs, subducting, spacing_rad)
        + ridge_push_torque(plate, inputs, divergent, spacing_rad)
        + collision_friction_torque(plate, inputs, collision_mask, spacing_rad)
    )
    drag_b, drag_k = basal_drag_coefficients(plate, world, spacing_rad)
    drag_k = drag_k + slab_drag_coefficient_matrix(inputs, subducting, spacing_rad)
    new_omega = integrate_omega(plate, explicit_torque, drag_b, drag_k, inertia, years)
    plate.set_omega(new_omega)

    increment = geometry.rotation_matrix_from_omega(plate.omega, years)
    plate.rotate(increment)

    new_points, _ = plate.all_points_and_elevation()
    return float(geometry.angular_distance(old_points, new_points).max())


def merge_omega(plate_a, inertia_a: np.ndarray, plate_b, inertia_b: np.ndarray) -> np.ndarray:
    """Angular-momentum-conserving blend for `LithospherePlate.merge_with`, replacing the
    base `Plate.merge_with`'s naive `(omega_a + omega_b) / 2`: L = I @ omega is additive
    across a fusion, so the combined plate's omega is `(I_a + I_b)^-1 @ (L_a + L_b)`, not a
    plain average of two possibly very differently-massed plates' rates."""
    l_total = lithosphere.angular_momentum(inertia_a, plate_a.omega) + lithosphere.angular_momentum(inertia_b, plate_b.omega)
    combined_omega = lithosphere.omega_from_angular_momentum(inertia_a + inertia_b, l_total)
    return mantle.clamp_rate(combined_omega)

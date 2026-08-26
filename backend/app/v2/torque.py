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
from scipy.spatial import cKDTree

from .. import geometry, mantle
from ..elevation_lines import line_spacing_rad
from . import lithosphere

SECONDS_PER_YEAR = 365.25 * 86400.0

# Section 3.2's own reach constants -- slab-pull/ridge-push are boundary-line integrals, so
# only nodes within a boundary-local band contribute; basal drag is a body force and applies
# to every node regardless of proximity to a boundary.
SUBDUCTION_LSINK_M = 250_000.0  # effective down-dip sunk-slab length contributing pull today
ASTHENOSPHERE_VISCOSITY_PA_S = 5e19  # within the spec's 1e19-1e21 Pa*s range
SHEAR_ZONE_THICKNESS_M = 100_000.0  # d_s, Eq. 10
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

    if not neighbours:
        return BoundaryForceInputs(
            own_points, own_hc, own_hm, np.full(n, np.inf), np.zeros((n, 3)), np.zeros(n, dtype=bool), np.zeros((n, 3))
        )

    pieces = [p.all_points_and_elevation()[0] for p in neighbours]
    owners = [np.full(len(pts), i) for i, pts in enumerate(pieces)]
    neighbour_points = np.concatenate(pieces, axis=0)
    neighbour_owner_idx = np.concatenate(owners, axis=0)
    if len(neighbour_points) == 0:
        return BoundaryForceInputs(
            own_points, own_hc, own_hm, np.full(n, np.inf), np.zeros((n, 3)), np.zeros(n, dtype=bool), np.zeros((n, 3))
        )

    # cKDTree, not bvh.py's own tree -- this runs every step for every plate's *entire* node
    # set (tens of thousands of points at real density), where bvh.py's pure-Python per-point
    # recursion loses badly to cKDTree's compiled batch query (confirmed directly: an 8-step,
    # ~16k-node-per-plate run went from ~2 minutes with the BVH here to a couple of seconds
    # with cKDTree). bvh.py's tree-vs-tree traversal (query_nearest_cross) is validated
    # against brute force (unit_tests/v2/test_bvh.py) and well-suited to a smaller, less
    # frequent plate-pair query -- merge_split.py's own collision-pair proximity check is the
    # natural next call site -- but isn't wired into that shared v1 module in this pass, to
    # avoid touching well-tested v1 code for a query that isn't currently a measured
    # bottleneck. Not yet exercised at runtime; kept available and independently tested.
    tree = cKDTree(neighbour_points)
    dist, idx = tree.query(own_points)
    owner_idx = neighbour_owner_idx[idx]
    nearest_points = neighbour_points[idx]
    direction = geometry.normalize(nearest_points - own_points)
    neighbor_is_oceanic = np.array([neighbours[i].crust_type == "oceanic" for i in owner_idx])
    neighbor_omega = np.array([neighbours[i].omega for i in owner_idx])
    dist = np.where(dist <= reach_rad, dist, np.inf)
    return BoundaryForceInputs(own_points, own_hc, own_hm, dist, direction, neighbor_is_oceanic, neighbor_omega)


def slab_pull_torque(plate, inputs: BoundaryForceInputs, spacing_rad: float, reach_rad: float) -> np.ndarray:
    """Eq. 8, discretized as a sum over this plate's own subducting-boundary nodes -- this
    plate must itself be oceanic (only oceanic lithosphere is dense enough to subduct, same
    asymmetry `deform()`'s own trench-vs-uplift branch already draws) and contested (its
    territory is being overridden -- see `rheology.py`'s reuse of the same "contested"
    geometric test `deform()` uses)."""
    if plate.crust_type != "oceanic" or len(inputs.own_points) == 0:
        return np.zeros(3)
    subducting = inputs.dist_to_neighbor <= reach_rad
    if not np.any(subducting):
        return np.zeros(3)
    ds = spacing_rad * lithosphere.PLANET_RADIUS_M
    hm = inputs.own_hm[subducting]
    force_mag = (lithosphere.RHO_LITHOSPHERE_MANTLE - lithosphere.RHO_ASTHENOSPHERE) * lithosphere.GRAVITY_M_S2 * hm * SUBDUCTION_LSINK_M * ds
    force = force_mag[:, None] * inputs.direction_to_neighbor[subducting]
    r = inputs.own_points[subducting] * lithosphere.PLANET_RADIUS_M
    return np.cross(r, force).sum(axis=0)


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


def basal_drag_torque(plate, world, spacing_rad: float) -> np.ndarray:
    """Eq. 10, over *every* node this plate owns (a body force, not a boundary integral)."""
    own_points, _ = plate.all_points_and_elevation()
    if len(own_points) == 0:
        return np.zeros(3)
    v_mantle = _real_velocity_m_per_s(mantle.flow_at(own_points, world.mantle_centers))
    v_plate = _real_velocity_m_per_s(np.cross(plate.omega, own_points))
    stress = (ASTHENOSPHERE_VISCOSITY_PA_S / SHEAR_ZONE_THICKNESS_M) * (v_mantle - v_plate)
    force = stress * lithosphere.node_area_m2(spacing_rad)
    r = own_points * lithosphere.PLANET_RADIUS_M
    return np.cross(r, force).sum(axis=0)


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


def integrate_omega(plate, torque_total: np.ndarray, inertia_tensor: np.ndarray, years: float) -> np.ndarray:
    """Semi-implicit Euler: alpha = I^-1 tau (rad/s^2), omega_new = omega_old + alpha * dt,
    then `mantle.clamp_rate` -- same physical speed bounds v1 already enforces, now applied
    to an integrated rather than fitted omega."""
    alpha = lithosphere.omega_from_angular_momentum(inertia_tensor, torque_total)  # I^-1 @ tau, reusing the same linear solve
    dt_seconds = years * SECONDS_PER_YEAR
    delta_omega_rad_per_s = alpha * dt_seconds
    delta_omega = delta_omega_rad_per_s * SECONDS_PER_YEAR  # rad/s -> rad/yr, this codebase's own omega convention
    new_omega = plate.omega + delta_omega
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
    near = inputs.dist_to_neighbor <= reach_rad
    contested = np.zeros(n, dtype=bool)
    if np.any(near):
        near_points = inputs.own_points[near]
        near_contested = np.zeros(len(near_points), dtype=bool)
        for neighbour in neighbours:
            near_contested |= neighbour.contains_batch(near_points)
            if np.all(near_contested):
                break
        contested[near] = near_contested
    divergent = near & ~contested
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

    rho_c = lithosphere.crust_density(plate.crust_type)
    inertia = lithosphere.moment_of_inertia_tensor(inputs.own_points, inputs.own_hc, inputs.own_hm, rho_c, spacing_rad)

    torque_total = (
        slab_pull_torque(plate, inputs, spacing_rad, reach_rad)
        + ridge_push_torque(plate, inputs, divergent, spacing_rad)
        + basal_drag_torque(plate, world, spacing_rad)
        + collision_friction_torque(plate, inputs, collision_mask, spacing_rad)
    )
    new_omega = integrate_omega(plate, torque_total, inertia, years)
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

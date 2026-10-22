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
-- each drag/pull term's own `c` coefficient folds in a `PLANET_RADIUS_M / SECONDS_PER_YEAR`
factor (twice, for the affine-in-omega `(b, K)` forms) at the point it crosses between the two
unit systems.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import cratons, geometry, mantle
from .boundary import TRANSFORM_RATE_THRESHOLD
from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad
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

# Contact drag (RFC 0008 solution 0, #399): collision resistance as a force per unit trench
# length. Copley, Avouac & Royer (2010) balance India's motion after the collision against
# ~5e12 N/m from the thickening Himalaya-Tibet, at ~4 cm/yr of convergence. The drag is linear
# in relative velocity, converted to per-node traction by the shared front-length weight, so a
# fully braked contact (factor 1, see `contact_drag_pa_s_per_m`) at the reference rate resists
# with CONTACT_RESISTANCE_N_PER_M per represented boundary metre. It
# replaces a 5e7 reference *stress* that multiplied the velocity in m/s, as basal drag's
# mu/d_s (5e14 Pa*s/m) does: about 1e7 times too weak to slow anything. The RFC 0008 gate
# (#395) measured a median 2.5-7.8e4 N/m at jammed contacts.
CONTACT_RESISTANCE_N_PER_M = 5e12
CONTACT_RESISTANCE_REF_CM_PER_YR = 4.0

# Reference (baseline) continental column elevation -- REFERENCE_HC/HM_CONTINENTAL_M's own
# isostatic_elevation, i.e. plates.BASE_CONTINENTAL_M (200m) recomputed from first principles
# the same way ABYSSAL_PLAIN_REFERENCE_ELEVATION_M is above, rather than importing plates.py's
# own copy of the same number.
CONTINENTAL_REFERENCE_ELEVATION_M = float(
    lithosphere.isostatic_elevation(
        np.array([lithosphere.REFERENCE_HC_CONTINENTAL_M]), np.array([lithosphere.REFERENCE_HM_CONTINENTAL_M]), lithosphere.RHO_CONTINENTAL_CRUST
    )[0]
)

# GitHub issue #176: `collision_uplift_multiplier` pumps up land volume/ice-cap coverage by
# making *existing* orogens taller rather than reclaiming new land area, because nothing in the
# model resists further shortening once a belt is already very tall. Real thickened orogenic
# crust resists further shortening more than fresh crust does (already near its own strength
# limit and gravitationally spreading under its own weight -- the same physical argument
# behind Hc's own hard `lithosphere.MAX_CRUSTAL_THICKNESS_M` ceiling), and Copley's resistance
# *is* the thickened plateau's. So a node's contact drag rises with its relief, at the velocity
# level, upstream of `apply_convergent_deformation`'s thickening rate. Full resistance at
# `COLLISION_RELIEF_FRICTION_SCALE_M` of relief above `CONTINENTAL_REFERENCE_ELEVATION_M`,
# roughly the Tibetan Plateau's mean-elevation scale (isostatic_elevation at 2x reference
# Hc/Hm computes to ~4.4km, ~4.4km relief above the 200m baseline) -- well below
# `MAX_CRUSTAL_THICKNESS_M`'s hard cap (2.4x reference), so the brake bites *before* a
# collision zone reaches the ceiling.
COLLISION_RELIEF_FRICTION_SCALE_M = 4000.0
# Where craton meets craton, no weak crust is left to crumple, so convergence has to slow or go
# elsewhere (#340). Retreat never consumes craton (cratons.retreat_allowed), so a contested
# node resists as hard as a Tibet-height plateau does, in proportion to the *weaker* side's
# craton strength (cratons.strength): its own, or the craton at the neighbour's nearest node.
# Collision resistance is one force across the contact, equal and opposite on both plates
# (Copley, Avouac & Royer 2010; Warners-Ruckstuhl et al. 2013), and it is the weak side's
# deformation that sets it: India, a cratonic indenter, is resisted by Tibet's thickening and
# the strong mantle under it (Clark 2012), not by its own strength. A craton pushing into weak
# crust therefore adds nothing here; the weak side's retreat and the relief term carry that
# resistance. The brake acts only at a contact; a keel adds no basal drag here. That leaves
# out a contested effect rather than settling it: Stoddard & Abbott (1996) found tectosphere
# sped plates up in the configurations they tested, which argues against a simple basal-drag
# slowdown but doesn't show keels never slow plates.

# Buoyant columns (RFC 0008 solution 1, #400). Continental crust thicker than ~15 km makes the
# lithosphere buoyant enough to jam a subduction zone, and 30 km of it on 100 km of lithosphere
# resists subduction (Cloos 1993). A continental-coded node being forced down at a contact
# jams in proportion between the two (`jam_factor`), and the jam brakes through the contact
# drag. The thresholds come from Cloos's crust thicknesses, not his densities: with the model's
# densities the reference continental column is more buoyant than his.
JAM_MIN_CRUST_M = 15_000.0
JAM_FULL_CRUST_M = 30_000.0

# Boundary-line integrals (slab-pull/ridge-push) share a `spacing_rad * PLANET_RADIUS_M`
# stretch among all band rows mapped to the same neighbour front node (#468).


@dataclass
class BoundaryForceInputs:
    """One plate's own near-boundary geometry, gathered once per `compute_torque` call and
    reused across the slab-pull/ridge-push/collision terms -- avoids each force needing its
    own separate neighbour/cKDTree pass."""

    own_points: np.ndarray  # (N, 3) unit vectors, this plate's own nodes
    own_hc: np.ndarray
    own_hm: np.ndarray
    own_crust_type_codes: np.ndarray  # (N,) int8, elevation_lines.CRUST_TYPE_* -- see lithosphere.node_crust_density
    dist_to_neighbor: np.ndarray  # (N,) angular distance (rad) to nearest other-plate node, +inf if none
    direction_to_neighbor: np.ndarray  # (N, 3) unit vector from own node toward that nearest neighbour point
    neighbor_is_oceanic: np.ndarray  # (N,) bool
    neighbor_omega: np.ndarray  # (N, 3) the owning neighbour plate's own omega, for relative-velocity terms
    # Which plate that nearest neighbour node belongs to, and its index in that plate's own node
    # order -- at any distance, not just within reach (a deep overlap's nearest node is still
    # its overrider's), -1 where there's no neighbour at all. The collision-polarity prepass
    # attributes evidence to the neighbour responsible (issue #318).
    neighbor_plate_id: np.ndarray | None = None  # (N,) int64
    neighbor_node_index: np.ndarray | None = None  # (N,) int64
    own_craton_m: np.ndarray | None = None  # (N,) craton_crust_m, for the keel's buoyancy (RFC 0008, #400)
    boundary_length_m: np.ndarray | None = None  # (N,) shared front-node length measure (#468)
    boundary_length_spacing_rad: float | None = None

    def __post_init__(self) -> None:
        n = len(self.own_points)
        if self.neighbor_plate_id is None:
            self.neighbor_plate_id = np.full(n, -1, dtype=np.int64)
        if self.neighbor_node_index is None:
            self.neighbor_node_index = np.full(n, -1, dtype=np.int64)
        if self.own_craton_m is None:
            self.own_craton_m = np.zeros(n)


@dataclass(frozen=True)
class PerNodeForceField:
    """Spatial external force terms on one plate's nodes, each in N.

    These are applied loads, not the internal stress that balances them. A shell or other
    load-transfer calculation must map them to stress before a rift strength comparison.
    Keeping the components separate preserves source attribution and lets that calculation
    use the exact forcing distribution instead of reconstructing it from a plate torque.
    """

    points_xyz: np.ndarray  # (N, 3) unit vectors
    areas_m2: np.ndarray  # (N,) accounting area
    slab_pull_n: np.ndarray  # (N, 3), N
    slab_drag_n: np.ndarray  # (N, 3), N
    ridge_push_n: np.ndarray  # (N, 3), N
    basal_drag_n: np.ndarray  # (N, 3), N
    contact_drag_n: np.ndarray  # (N, 3), N

    @property
    def total_n(self) -> np.ndarray:
        return self.slab_pull_n + self.slab_drag_n + self.ridge_push_n + self.basal_drag_n + self.contact_drag_n

    @property
    def total_traction_pa(self) -> np.ndarray:
        """Area-average applied traction at each node; not an in-plane stress tensor."""
        return self.total_n / self.areas_m2[:, None]

    @property
    def torque_nm(self) -> np.ndarray:
        """Moment of the spatial force field about the planet center, in N·m."""
        return np.cross(self.points_xyz * lithosphere.PLANET_RADIUS_M, self.total_n).sum(axis=0)


class BoundarySearchCache:
    """One step's per-(plate, neighbour) boundary searches, shared by the collision-polarity
    prepass (collision_polarity.py, issue #318) and `deform()`'s own `boundary_context`, so the
    prepass doesn't add a second full neighbour search to every step.

    A cached result is reused only while both plates' node sets are unchanged, judged by the
    identity of their cached node k-d trees (`Plate.get_node_kdtree`), which every geometry
    change drops. A reused result is therefore exactly what a fresh search would return, and
    `deform()` behaves as it did without the cache. A neighbour that already deformed this
    step (cells peeled or grown) is searched afresh. Lives on `World.boundary_search_cache`
    from the prepass to the end of the deform pass only."""

    def __init__(self) -> None:
        # (own id, neighbour id) -> (own tree, neighbour tree, dist, idx): nearest neighbour
        # node of every own node.
        self._nearest: dict[tuple[int, int], tuple] = {}
        # (own id, neighbour id) -> (own tree, neighbour tree, evaluated, inside): which own
        # nodes lie inside the neighbour's territory, over the nodes asked about so far.
        self._inside: dict[tuple[int, int], tuple] = {}
        # neighbour id -> (tree, centroid, radius) of its node cloud.
        self._spheres: dict[int, tuple] = {}
        # (own id, neighbour id, threshold) -> (own tree, neighbour tree, is neighbour).
        self._within: dict[tuple[int, int, float], tuple] = {}
        self.searches = 0
        self.reused = 0
        # (searches, reused) when the prepass finished, so the deform pass's share can be told
        # apart -- see collision_polarity.finish_deform_pass.
        self.prepass_mark = (0, 0)

    def neighbours(self, plate, others: list, threshold_rad: float) -> list:
        """`plate.get_neighbours(others, threshold_rad)`, pair by pair through the cache. The
        test is per pair, so the result and its order match one call over `others`."""
        own_tree = plate.get_node_kdtree()
        others = [other for other in others if other.plate_id != plate.plate_id]
        trees = [other.get_node_kdtree() for other in others]
        keys = [(plate.plate_id, other.plate_id, float(threshold_rad)) for other in others]
        stale = [
            i
            for i, key in enumerate(keys)
            if own_tree is None
            or (entry := self._within.get(key)) is None
            or entry[0] is not own_tree
            or entry[1] is not trees[i]
        ]
        if stale:
            # One call over every pair the cache can't answer: the test is per pair, and one
            # call computes this plate's own outline sphere once rather than per pair.
            found = {other.plate_id for other in plate.get_neighbours([others[i] for i in stale], threshold_rad=threshold_rad)}
            for i in stale:
                self._within[keys[i]] = (own_tree, trees[i], others[i].plate_id in found)
        return [other for other, key in zip(others, keys) if self._within[key][2]]

    def nearest(self, plate, own_tree, neighbour, tree, own_points: np.ndarray, workers: int) -> tuple[np.ndarray, np.ndarray]:
        key = (plate.plate_id, neighbour.plate_id)
        entry = self._nearest.get(key)
        if entry is not None and entry[0] is own_tree and entry[1] is tree and own_tree is not None:
            self.reused += 1
            return entry[2], entry[3]
        self.searches += 1
        dist, idx = tree.query(own_points, workers=workers)
        self._nearest[key] = (own_tree, tree, dist, idx)
        return dist, idx

    def sphere(self, neighbour) -> tuple[np.ndarray, float] | None:
        tree = neighbour.get_node_kdtree()
        if tree is None:
            return None
        entry = self._spheres.get(neighbour.plate_id)
        if entry is not None and entry[0] is tree:
            return entry[1], entry[2]
        centroid, radius = geometry.bounding_sphere(tree.data)
        self._spheres[neighbour.plate_id] = (tree, centroid, radius)
        return centroid, radius

    def inside(self, plate, neighbour, own_points: np.ndarray, mask: np.ndarray) -> np.ndarray:
        """`neighbour.contains_batch` over `own_points[mask]`, as a full-length mask (False off
        `mask`). Only the masked nodes not yet tested against this neighbour are tested."""
        own_tree = plate.get_node_kdtree()
        tree = neighbour.get_node_kdtree()
        key = (plate.plate_id, neighbour.plate_id)
        entry = self._inside.get(key)
        if entry is None or entry[0] is not own_tree or entry[1] is not tree or own_tree is None:
            entry = (own_tree, tree, np.zeros(len(own_points), dtype=bool), np.zeros(len(own_points), dtype=bool))
            self._inside[key] = entry
        evaluated, inside = entry[2], entry[3]
        need = mask & ~evaluated
        if np.any(need):
            inside[need] = neighbour.contains_batch(own_points[need])
            evaluated |= need
        return inside & mask


def gather_boundary_force_inputs(
    plate, neighbours: list, spacing_rad: float, reach_rad: float, cache: BoundarySearchCache | None = None
) -> BoundaryForceInputs:
    own_points, _ = plate.all_points_and_elevation()
    own_hc = plate.collect("crustal_thickness_m")
    own_hm = plate.collect("mantle_lithosphere_thickness_m")
    own_crust_type_codes = plate.collect("crust_type_code")
    own_craton_m = plate.collect("craton_crust_m")
    n = len(own_points)
    if n == 0:
        empty3 = np.zeros((0, 3))
        return BoundaryForceInputs(empty3, np.zeros(0), np.zeros(0), np.zeros(0, dtype=np.int8), np.full(0, np.inf), empty3, np.zeros(0, dtype=bool), empty3)

    no_neighbour = BoundaryForceInputs(
        own_points, own_hc, own_hm, own_crust_type_codes, np.full(n, np.inf), np.zeros((n, 3)), np.zeros(n, dtype=bool), np.zeros((n, 3)),
        own_craton_m=own_craton_m,
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
    best_index = np.full(n, -1, dtype=np.int64)
    own_tree = plate.get_node_kdtree() if cache is not None else None
    for i, neighbour in enumerate(neighbours):
        tree = neighbour.get_node_kdtree()
        if tree is None:
            continue
        if cache is None:
            dist, idx = tree.query(own_points, workers=workers)
        else:
            dist, idx = cache.nearest(plate, own_tree, neighbour, tree, own_points, workers)
        closer = dist < best_dist
        best_dist[closer] = dist[closer]
        best_point[closer] = tree.data[idx[closer]]
        best_owner[closer] = i
        best_index[closer] = idx[closer]

    if not np.any(best_owner >= 0):
        return no_neighbour

    direction = geometry.normalize(best_point - own_points)
    neighbor_is_oceanic = np.zeros(n, dtype=bool)
    neighbor_omega = np.zeros((n, 3))
    neighbor_plate_id = np.full(n, -1, dtype=np.int64)
    for i, neighbour in enumerate(neighbours):
        owned = best_owner == i
        if not np.any(owned):
            continue
        # The nearest neighbour node's own crust, not its plate's majority type: an ocean
        # margin on a mostly continental plate is a trench to subduct under, and a continent on
        # a mostly oceanic plate is a collision (RFC 0009, #447).
        codes = neighbour.collect("crust_type_code")[best_index[owned]]
        neighbor_is_oceanic[owned] = ~effective_is_continental_from_codes(codes, neighbour.crust_type == "continental")
        neighbor_omega[owned] = neighbour.omega
        neighbor_plate_id[owned] = neighbour.plate_id
    dist = np.where(best_dist <= reach_rad, best_dist, np.inf)
    return BoundaryForceInputs(
        own_points, own_hc, own_hm, own_crust_type_codes, dist, direction, neighbor_is_oceanic, neighbor_omega,
        neighbor_plate_id, best_index, own_craton_m,
    )


def subducting_boundary_mask(plate, inputs: BoundaryForceInputs, reach_rad: float) -> np.ndarray:
    """This plate's own boundary-band nodes that are genuinely subducting: the node's crust goes
    down (`down_going_crust`), it is within `reach_rad` of a neighbour, *and*
    the two plates are actually converging there (`closing_rate` past the same
    `TRANSFORM_RATE_THRESHOLD` `merge_split` uses). The convergence check is what this used to
    be missing: slab pull was applied to the plate's *entire* near-neighbour band, ridge and
    transform stretches included, so a plate spreading at three of its four edges still felt a
    full-perimeter pull and railed at MAX_PLATE_RATE. A ridge does not pull a slab down."""
    band = inputs.dist_to_neighbor <= reach_rad
    return band & (closing_rad_per_yr(plate, inputs) > TRANSFORM_RATE_THRESHOLD) & down_going_crust(plate, inputs, reach_rad)


def own_is_continental(plate, inputs: BoundaryForceInputs) -> np.ndarray:
    return effective_is_continental_from_codes(inputs.own_crust_type_codes, plate.crust_type == "continental")


def down_going_crust(plate, inputs: BoundaryForceInputs, reach_rad: float) -> np.ndarray:
    """Per node, whether its crust goes down when it converges on a neighbour (RFC 0009, #447):
    every oceanic-coded node, since only oceanic lithosphere is dense enough to descend (the
    same asymmetry `deform()`'s trench-vs-uplift branch draws), and the continental nodes
    `slab_attached_continent` finds. A classification of crust, not a contact mask: it holds
    interior ocean floor too, and callers add the boundary band and closing they need
    (`subducting_boundary_mask`, `jam_mask`).

    Before #447 the plate's majority type decided for every node: all of an oceanic-majority
    plate went down, at any contact, and none of a continental-majority plate did. A
    single-type plate gives the same answer as before."""
    return ~own_is_continental(plate, inputs) | slab_attached_continent(plate, inputs, reach_rad)


def slab_attached_continent(plate, inputs: BoundaryForceInputs, reach_rad: float) -> np.ndarray:
    """Per node, whether it is continental crust attached to its own plate's slab at its
    contact: some oceanic-coded node of the plate is in the band facing the same neighbour
    plate, at a contact that isn't diverging. The slab drags the continent in behind it, a
    terrane riding ocean floor or a continent behind its own ocean (India), and the
    continent's buoyancy then cuts the pull and jams the trench (`column_buoyancy_pa`,
    `jam_mask`). Not diverging rather than converging past +1 cm/yr, for the reason `jam_mask`
    gives: a jam that slows the contact must not detach the continent from its slab."""
    continental = own_is_continental(plate, inputs)
    return continental & _faces_slab_neighbour(inputs, ~continental & _jam_contact(plate, inputs, reach_rad))


def _faces_slab_neighbour(inputs: BoundaryForceInputs, slab: np.ndarray) -> np.ndarray:
    """Per node, whether its nearest neighbour plate is one of those faced by a `slab` node."""
    slab_ids = np.unique(inputs.neighbor_plate_id[slab])
    return np.isin(inputs.neighbor_plate_id, slab_ids[slab_ids >= 0])


def closing_rad_per_yr(plate, inputs: BoundaryForceInputs) -> np.ndarray:
    """Per node, the closing rate toward the nearest neighbour in rad/yr, the unit
    `TRANSFORM_RATE_THRESHOLD` is in (positive = converging)."""
    v_self = np.cross(plate.omega, inputs.own_points)
    v_neighbor = np.cross(inputs.neighbor_omega, inputs.own_points)
    return np.sum((v_self - v_neighbor) * inputs.direction_to_neighbor, axis=-1)


def jam_factor(hc_cont: np.ndarray) -> np.ndarray:
    """j in [0, 1] per node from its continental crust (`continental_crust_m`): 0 up to
    JAM_MIN_CRUST_M, 1 from JAM_FULL_CRUST_M, linear between (RFC 0008, #400)."""
    return np.clip((np.asarray(hc_cont, dtype=float) - JAM_MIN_CRUST_M) / (JAM_FULL_CRUST_M - JAM_MIN_CRUST_M), 0.0, 1.0)


def forced_down_mask(plate, world, inputs: BoundaryForceInputs, reach_rad: float) -> np.ndarray:
    """Per node of `plate`, whether it is on the side forced down at a contact: crust that goes
    down (`down_going_crust`: an oceanic node, or a continental one behind its plate's own
    slab), and the lower-plate nodes of polarized fronts (`lower_plate_mask`). Like
    `down_going_crust`, not a contact mask: `jam_mask` adds the contact."""
    return down_going_crust(plate, inputs, reach_rad) | lower_plate_mask(plate, world, len(inputs.own_points))


def lower_plate_mask(plate, world, n: int) -> np.ndarray:
    """Per node of `plate` (`n` nodes), whether it is on the lower plate of a polarized
    continental front (`collision_polarity`).

    Plates shift before the polarity prepass freezes this step's fronts, so this reads the last
    frozen frame, the previous step's, and matches its lower-plate cells to the plate's cells
    now by cell key: cells lost since then drop out, cells gained join only at the next
    prepass. A plate the frame doesn't know (new since a split or merge, or a loaded save)
    has no lower-plate nodes until then."""
    frame = getattr(world, "collision_polarity_frame", None)
    masks = frame.masks.get(plate.plate_id) if frame is not None else None
    if masks is None or not np.any(masks.lower):
        return np.zeros(n, dtype=bool)
    return np.isin(plate.cell_keys, masks.node_keys[masks.lower])


def jam_mask(plate, inputs: BoundaryForceInputs, forced_down: np.ndarray, reach_rad: float) -> np.ndarray:
    """RFC 0008's jam mask (#400): continental-coded nodes on the side being forced down
    (`forced_down_mask`), in the boundary band, at a contact that isn't diverging (closing
    >= -TRANSFORM_RATE_THRESHOLD).

    Independent of `subducting_boundary_mask`, which needs closing above +1 cm/yr: a jam read
    from it would lose its resistance as soon as the brake slowed the contact below 1 cm/yr,
    the driving forces would speed it back up, and the contact would flicker. A band through
    zero keeps the resistance while the contact stalls; the drag vanishes with the relative
    velocity, so it can't hold together a contact that is pulling apart."""
    return own_is_continental(plate, inputs) & forced_down & _jam_contact(plate, inputs, reach_rad)


def _jam_contact(plate, inputs: BoundaryForceInputs, reach_rad: float) -> np.ndarray:
    return (inputs.dist_to_neighbor <= reach_rad) & (closing_rad_per_yr(plate, inputs) >= -TRANSFORM_RATE_THRESHOLD)


def continental_crust_m(plate, inputs: BoundaryForceInputs) -> np.ndarray:
    """Per node, Hc where the node is continental-coded, else 0. Oceanic crust is left out of
    a column's buoyancy: it turns to eclogite as it subducts (Cloos 1993)."""
    return np.where(own_is_continental(plate, inputs), inputs.own_hc, 0.0)


def column_buoyancy_pa(hc_cont: np.ndarray, hm: np.ndarray, craton_m: np.ndarray) -> np.ndarray:
    """Net negative buoyancy per unit area of a column, in Pa (positive sinks), RFC 0008
    solution 1 (#400): `g [(rho_m - rho_a) Hm (1 - s) - (rho_a - rho_cc) Hc_cont]`.

    The mantle lithosphere's excess density pulls; continental crust (`hc_cont`, from
    `continental_crust_m`) holds the column up. A cratonic keel is chemically depleted, which
    offsets its cold, dense thermal state (Jordan 1978), so a full-strength craton's keel
    (`cratons.strength` s = 1) is neutral. For ocean floor (no continental crust, no craton)
    it is the mantle term alone, slab pull's buoyancy before this change."""
    s = cratons.strength(craton_m)
    return lithosphere.GRAVITY_M_S2 * (
        (lithosphere.RHO_LITHOSPHERE_MANTLE - lithosphere.RHO_ASTHENOSPHERE) * hm * (1.0 - s)
        - (lithosphere.RHO_ASTHENOSPHERE - lithosphere.RHO_CONTINENTAL_CRUST) * hc_cont
    )


def boundary_node_length_m(inputs: BoundaryForceInputs, spacing_rad: float) -> np.ndarray:
    """Length represented by each own-side node along its nearest neighbour front.

    Several rows in a boundary band can map to the same neighbour node. They share that
    front node's one along-boundary interval, rather than each contributing a full `ds`.
    This keeps line loads invariant to band width and repeated sampling (#468). Inputs
    without a valid neighbour-node identity retain the historical one-`ds` fallback for
    standalone callers and small equation tests. The result is memoized on `inputs` for
    the requested spacing; a different spacing replaces the cached array.
    """
    if inputs.boundary_length_m is not None and inputs.boundary_length_spacing_rad == spacing_rad:
        return inputs.boundary_length_m
    n = len(inputs.own_points)
    ds = spacing_rad * lithosphere.PLANET_RADIUS_M
    lengths = np.full(n, ds, dtype=float)
    valid = (inputs.neighbor_plate_id >= 0) & (inputs.neighbor_node_index >= 0) & np.isfinite(inputs.dist_to_neighbor)
    if np.any(valid):
        front_nodes = np.column_stack((inputs.neighbor_plate_id[valid], inputs.neighbor_node_index[valid]))
        _, inverse, counts = np.unique(front_nodes, axis=0, return_inverse=True, return_counts=True)
        lengths[valid] = ds / counts[inverse]
    inputs.boundary_length_m = lengths
    inputs.boundary_length_spacing_rad = spacing_rad
    return lengths


def slab_pull_force_n(plate, inputs: BoundaryForceInputs, subducting_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """Per subducting node (`subducting_mask`), the magnitude of its slab pull in N: its
    column's buoyancy (`column_buoyancy_pa`) times `SUBDUCTION_LSINK_M` and its stretch of
    boundary. A buoyant column pulls nothing (RFC 0008, #400): before, every node pulled with
    its Hm alone, so a continental terrane or a craton was dragged into the trench at full
    strength, and a thick cratonic keel pulled *harder* than ocean floor."""
    ds = boundary_node_length_m(inputs, spacing_rad)[subducting_mask]
    buoyancy = column_buoyancy_pa(
        continental_crust_m(plate, inputs)[subducting_mask], inputs.own_hm[subducting_mask], inputs.own_craton_m[subducting_mask]
    )
    return np.clip(buoyancy, 0.0, None) * SUBDUCTION_LSINK_M * ds


def slab_pull_force_vectors_n(plate, inputs: BoundaryForceInputs, subducting_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """Spatial slab-pull forces on all own nodes, in N, directed toward the neighbour."""
    forces = np.zeros_like(inputs.own_points, dtype=float)
    if np.any(subducting_mask):
        forces[subducting_mask] = (
            slab_pull_force_n(plate, inputs, subducting_mask, spacing_rad)[:, None]
            * inputs.direction_to_neighbor[subducting_mask]
        )
    return forces


def slab_pull_torque(plate, inputs: BoundaryForceInputs, subducting_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """Eq. 8, discretized as a sum over this plate's own subducting-boundary nodes
    (`subducting_mask`, from `subducting_boundary_mask` -- passed in rather than recomputed so
    it stays in lockstep with the matching `slab_drag_torque` resistance), each pulling with
    `slab_pull_force_n`."""
    force = slab_pull_force_vectors_n(plate, inputs, subducting_mask, spacing_rad)
    return np.cross(inputs.own_points * lithosphere.PLANET_RADIUS_M, force).sum(axis=0)


def slab_drag_force_vectors_n(inputs: BoundaryForceInputs, subducting_mask: np.ndarray, spacing_rad: float, omega: np.ndarray) -> np.ndarray:
    """Spatial deep-mantle slab drag on a plate's subducting nodes, in N."""
    forces = np.zeros_like(inputs.own_points, dtype=float)
    if not np.any(subducting_mask):
        return forces
    points = inputs.own_points[subducting_mask]
    lengths = boundary_node_length_m(inputs, spacing_rad)[subducting_mask]
    area_m2 = SUBDUCTION_LSINK_M * lengths
    velocity_m_s = np.cross(np.asarray(omega, dtype=float), points) * lithosphere.PLANET_RADIUS_M / SECONDS_PER_YEAR
    traction_pa = -(SLAB_MANTLE_VISCOSITY_PA_S / SHEAR_ZONE_THICKNESS_M) * velocity_m_s
    forces[subducting_mask] = traction_pa * area_m2[:, None]
    return forces


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
    ds = boundary_node_length_m(inputs, spacing_rad)[subducting_mask]
    c = (
        (SLAB_MANTLE_VISCOSITY_PA_S / SHEAR_ZONE_THICKNESS_M)
        * (SUBDUCTION_LSINK_M * ds)
        * lithosphere.PLANET_RADIUS_M**2
        / SECONDS_PER_YEAR
    )
    return c.sum() * np.eye(3) - np.einsum("n,ni,nj->ij", c, pts, pts)


def slab_drag_torque(plate, inputs: BoundaryForceInputs, subducting_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """`slab_drag_coefficient_matrix` evaluated at the plate's current `omega` (`-K @ omega`)
    -- a standalone entry point for tests / callers that just want the resistive torque;
    `shift_plate` itself folds the `K` into `integrate_omega` for implicit stability."""
    return -slab_drag_coefficient_matrix(inputs, subducting_mask, spacing_rad) @ plate.omega


def ridge_push_force_vectors_n(plate, inputs: BoundaryForceInputs, divergent_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """Spatial ridge-push force on all own nodes, in N, directed into the plate."""
    force = np.zeros_like(inputs.own_points, dtype=float)
    if not np.any(divergent_mask):
        return force
    ds = boundary_node_length_m(inputs, spacing_rad)[divergent_mask]
    rho_c = lithosphere.node_crust_density(inputs.own_crust_type_codes[divergent_mask], plate.crust_type)
    z = lithosphere.isostatic_elevation(inputs.own_hc[divergent_mask], inputs.own_hm[divergent_mask], rho_c)
    e_r = np.clip(z - ABYSSAL_PLAIN_REFERENCE_ELEVATION_M, 0.0, None)
    force_mag = 0.5 * lithosphere.GRAVITY_M_S2 * (lithosphere.RHO_ASTHENOSPHERE - lithosphere.RHO_WATER) * e_r**2 * ds
    direction = -inputs.direction_to_neighbor[divergent_mask]
    # No neighbour within reach (an isolated ridge segment, direction left at zero) -- fall
    # back to the node's own outward radial tangent component is unnecessary here since
    # push-direction ambiguity without any neighbour reference just contributes no torque
    # (force magnitude still nonzero, but a zero direction vector cancels it) -- acceptable,
    # since "no neighbour nearby" also means this node isn't really acting as part of an
    # active spreading boundary the plate as a whole would feel a coherent push from.
    force[divergent_mask] = force_mag[:, None] * direction
    return force


def ridge_push_torque(plate, inputs: BoundaryForceInputs, divergent_mask: np.ndarray, spacing_rad: float) -> np.ndarray:
    """Eq. 9; sum the moments of the spatial ridge-push force vectors."""
    force = ridge_push_force_vectors_n(plate, inputs, divergent_mask, spacing_rad)
    r = inputs.own_points * lithosphere.PLANET_RADIUS_M
    return np.cross(r, force).sum(axis=0)


def basal_drag_force_vectors_n(plate, world, spacing_rad: float) -> np.ndarray:
    """Spatial basal-drag force on each owned node, in N, before plate aggregation."""
    own_points, _ = plate.all_points_and_elevation()
    if len(own_points) == 0:
        return np.zeros((0, 3))
    own_velocity_m_s = np.cross(plate.omega, own_points) * lithosphere.PLANET_RADIUS_M / SECONDS_PER_YEAR
    mantle_velocity_m_s = mantle.flow_at(own_points, world.mantle_centers) * lithosphere.PLANET_RADIUS_M / SECONDS_PER_YEAR
    basal_coefficient_pa_s_per_m = ASTHENOSPHERE_VISCOSITY_PA_S / SHEAR_ZONE_THICKNESS_M
    return (
        basal_coefficient_pa_s_per_m
        * (mantle_velocity_m_s - own_velocity_m_s)
        * plate.accounting_areas_m2(spacing_rad)[:, None]
    )


def basal_drag_coefficients(plate, world, spacing_rad: float) -> tuple[np.ndarray, np.ndarray]:
    """Eq. 10's basal drag written as its exact affine-in-`omega` form, `tau(omega) = b - K @ omega`,
    over *every* node this plate owns (a body force, not a boundary integral).

    Each node's drag force is `F_i = (mu / d_s) * (v_mantle_i - v_plate_i) * A_i` with
    `v_plate_i = (omega x p_i) * R / SECONDS_PER_YEAR` (the real m/s velocity at the planet's
    actual radius `R`), so the whole-plate drag torque `sum_i (R p_i) x F_i` is affine in
    `omega`:

        K = sum_i c_i (I3 - p_i p_i^T)      -- symmetric, positive-semidefinite
        b = sum_i c_i p_i x flow_i          -- `flow_i` = mantle.flow_at, this codebase's rad/yr units
        c_i = (mu / d_s) * A_i * R^2 / SECONDS_PER_YEAR

    (`c`'s `R^2 / SECONDS_PER_YEAR` folds in the two rad/yr-to-m/s conversions (`v_mantle_i` and
    `v_plate_i` each contribute one power of `R / SECONDS_PER_YEAR`) plus the `R` lever arm;
    `sum_i p_i x (omega x p_i) == (sum_i I3 - p_i p_i^T) @ omega`.)

    Split out from the plain torque so `integrate_omega` can treat this term *implicitly*: it
    is by far the stiffest in the balance -- the asthenosphere coupling relaxes a plate toward
    the local mantle-flow rate in a small fraction of one tectonic step, so an explicit Euler
    step on it overshoots by ~19 orders of magnitude and `clamp_rate` then pins every plate at
    `MAX_PLATE_RATE` on step one, forever (the historical "every plate railed at MAX" bug)."""
    own_points, _ = plate.all_points_and_elevation()
    if len(own_points) == 0:
        return np.zeros(3), np.zeros((3, 3))
    # (N,) -- `c` above per node, each with its own area `A_i` (`plate.accounting_areas_m2`).
    c = (
        (ASTHENOSPHERE_VISCOSITY_PA_S / SHEAR_ZONE_THICKNESS_M)
        * plate.accounting_areas_m2(spacing_rad)
        * lithosphere.PLANET_RADIUS_M**2
        / SECONDS_PER_YEAR
    )
    flow = mantle.flow_at(own_points, world.mantle_centers)
    b = np.einsum("n,ni->i", c, np.cross(own_points, flow))
    k = c.sum() * np.eye(3) - np.einsum("n,ni,nj->ij", c, own_points, own_points)
    return b, k


def basal_drag_torque(plate, world, spacing_rad: float) -> np.ndarray:
    """Eq. 10's basal drag evaluated at the plate's current `omega` -- `b - K @ omega` from
    `basal_drag_coefficients`. A standalone entry point for tests / callers that just want the
    torque; `shift_plate` itself uses the split `(b, K)` form for implicit integration."""
    b, k = basal_drag_coefficients(plate, world, spacing_rad)
    return b - k @ plate.omega


def contact_drag_pa_s_per_m(
    plate,
    inputs: BoundaryForceInputs,
    collision_mask: np.ndarray,
    spacing_rad: float,
    craton_m: np.ndarray | None = None,
    jam: np.ndarray | None = None,
) -> np.ndarray:
    """Per node, the contact drag coefficient in Pa*s/m: traction per unit relative velocity.

    `eta_i = CONTACT_RESISTANCE_N_PER_M * dL_i / (A_i * v_ref)`, where `dL_i` is this row's
    shared front-node length and `A_i` its accounting area. Integrating `eta_i * A_i * v_ref`
    returns the stated force per represented boundary metre (RFC 0008, #468). g is the largest of the
    contact's jam factor (`jam`, from `contact_jam_factor`, #400) and, at contested nodes
    (`collision_mask`), the node's relief factor (COLLISION_RELIEF_FRICTION_SCALE_M) and the
    contact's craton strength (`craton_m`, per node the weaker side's `craton_crust_m`), so it
    stays in [0, 1]: one bounded factor in place of the multiplicative overlap, relief and
    craton gains tuned against the old inert coefficient. A contact with no buoyant crust
    forced down, no relief and no craton resists with nothing until relief builds, as
    Copley's resistance comes from mountain building."""
    n = len(inputs.own_points)
    factor = np.zeros(n) if jam is None else np.asarray(jam, dtype=float).copy()
    if np.any(collision_mask):
        rho_c = lithosphere.node_crust_density(inputs.own_crust_type_codes[collision_mask], plate.crust_type)
        elevation = lithosphere.isostatic_elevation(inputs.own_hc[collision_mask], inputs.own_hm[collision_mask], rho_c)
        relief = np.clip(elevation - CONTINENTAL_REFERENCE_ELEVATION_M, 0.0, None)
        relief_factor = np.clip(relief / COLLISION_RELIEF_FRICTION_SCALE_M, 0.0, 1.0)
        craton_factor = cratons.strength(craton_m[collision_mask]) if craton_m is not None else 0.0
        factor[collision_mask] = np.maximum(factor[collision_mask], np.maximum(relief_factor, craton_factor))
    factor = np.clip(factor, 0.0, 1.0)
    v_ref_m_per_s = CONTACT_RESISTANCE_REF_CM_PER_YR * 0.01 / SECONDS_PER_YEAR
    node_area_m2 = plate.accounting_areas_m2(spacing_rad)
    boundary_length_m = boundary_node_length_m(inputs, spacing_rad)
    return CONTACT_RESISTANCE_N_PER_M * boundary_length_m / (node_area_m2 * v_ref_m_per_s) * factor


def collision_drag_coefficients(
    plate,
    inputs: BoundaryForceInputs,
    collision_mask: np.ndarray,
    spacing_rad: float,
    craton_m: np.ndarray | None = None,
    jam: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Resistive drag at continent-continent contested nodes (`collision_mask`) and jammed
    contacts (`jam`), in `basal_drag_coefficients`' own affine `tau(omega) = b - K @ omega`
    form, with each node's coefficient from `contact_drag_pa_s_per_m`.

    GitHub issue #176: this used to be a standalone torque handed to `integrate_omega` as an
    explicit forcing term, though it was already affine in omega (`relative = v_self -
    v_neighbor`, linear in `plate.omega`). The implicit basal-drag term is stiff enough to swamp
    *any* purely explicit torque, so raising its gains barely moved a #176 sweep. In the same
    `(b, K)` form `basal_drag_coefficients`/`slab_drag_coefficient_matrix` use, `shift_plate`
    folds it into `drag_k`, where it competes with basal drag on equal footing. Being implicit
    and linear in the relative velocity, it slows a driven plate toward the speed where driving
    and resistance balance but never reverses it in one step."""
    if len(inputs.own_points) == 0:
        return np.zeros(3), np.zeros((3, 3))
    coefficient = contact_drag_pa_s_per_m(plate, inputs, collision_mask, spacing_rad, craton_m, jam)
    drag = coefficient > 0.0
    if not np.any(drag):
        return np.zeros(3), np.zeros((3, 3))
    own_points = inputs.own_points[drag]
    neighbor_omega = inputs.neighbor_omega[drag]
    area_m2 = plate.accounting_areas_m2(spacing_rad)[drag]
    c = coefficient[drag] * area_m2 * lithosphere.PLANET_RADIUS_M**2 / SECONDS_PER_YEAR
    raw_neighbor = np.cross(neighbor_omega, own_points)  # rad/yr-equivalent, not yet real m/s -- see basal_drag_coefficients
    b = np.einsum("n,ni->i", c, np.cross(own_points, raw_neighbor))
    k = np.einsum("n,nij->ij", c, np.eye(3)[None, :, :] - np.einsum("ni,nj->nij", own_points, own_points))
    return b, k


def contact_drag_force_vectors_n(
    plate,
    inputs: BoundaryForceInputs,
    collision_mask: np.ndarray,
    spacing_rad: float,
    craton_m: np.ndarray | None = None,
    jam: np.ndarray | None = None,
) -> np.ndarray:
    """Spatial contact/jam drag forces on all own nodes, in N."""
    force = np.zeros_like(inputs.own_points, dtype=float)
    coefficient = contact_drag_pa_s_per_m(plate, inputs, collision_mask, spacing_rad, craton_m, jam)
    active = coefficient > 0.0
    if np.any(active):
        relative_velocity_m_s = (
            np.cross(inputs.neighbor_omega[active] - plate.omega, inputs.own_points[active])
            * lithosphere.PLANET_RADIUS_M
            / SECONDS_PER_YEAR
        )
        force[active] = (
            coefficient[active, None]
            * plate.accounting_areas_m2(spacing_rad)[active, None]
            * relative_velocity_m_s
        )
    return force


def gather_per_node_force_field(
    plate,
    world,
    inputs: BoundaryForceInputs,
    spacing_rad: float,
    subducting_mask: np.ndarray,
    divergent_mask: np.ndarray,
    collision_mask: np.ndarray,
    craton_m: np.ndarray | None = None,
    jam: np.ndarray | None = None,
) -> PerNodeForceField:
    """Return the existing torque engine's force terms before spatial aggregation.

    The resulting node vectors preserve each source's existing units and formulation. They
    are inputs to a later internal-load calculation, not a local stress estimate themselves.
    """
    areas_m2 = np.asarray(plate.accounting_areas_m2(spacing_rad), dtype=float)
    return PerNodeForceField(
        points_xyz=np.asarray(inputs.own_points, dtype=float),
        areas_m2=areas_m2,
        slab_pull_n=slab_pull_force_vectors_n(plate, inputs, subducting_mask, spacing_rad),
        slab_drag_n=slab_drag_force_vectors_n(inputs, subducting_mask, spacing_rad, plate.omega),
        ridge_push_n=ridge_push_force_vectors_n(plate, inputs, divergent_mask, spacing_rad),
        basal_drag_n=basal_drag_force_vectors_n(plate, world, spacing_rad),
        contact_drag_n=contact_drag_force_vectors_n(plate, inputs, collision_mask, spacing_rad, craton_m, jam),
    )


def collision_friction_torque(
    plate,
    inputs: BoundaryForceInputs,
    collision_mask: np.ndarray,
    spacing_rad: float,
    craton_m: np.ndarray | None = None,
    jam: np.ndarray | None = None,
) -> np.ndarray:
    """`collision_drag_coefficients` evaluated at the plate's current `omega` (`b - K @ omega`)
    -- a standalone entry point for tests / callers that just want the resistive torque, same
    role `slab_drag_torque`/`basal_drag_torque` play for their own `(b, K)` splits; `shift_plate`
    itself uses the split form for implicit integration."""
    b, k = collision_drag_coefficients(plate, inputs, collision_mask, spacing_rad, craton_m, jam)
    return b - k @ plate.omega


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

    `(drag_b, drag_k)` is `basal_drag_coefficients`' own affine split of Eq. 10's drag (plus
    `slab_drag_coefficient_matrix`/`collision_drag_coefficients`, folded into the same `drag_k`/
    `drag_b` by `shift_plate` for the same reason -- see `collision_drag_coefficients`' own
    docstring for why an explicit torque structurally can't compete with this implicit term),
    `tau(omega) = b - K @ omega`; `explicit_torque` is every other (bounded, geometry-driven)
    torque. Backward Euler on the drag term is unconditionally stable for any `years` -- an
    explicit step on it does not converge at any step size this simulation uses (see
    `basal_drag_coefficients`). Reduces to the old explicit `omega_old + g I^-1 tau` when
    `K = 0, b = 0`. The final `mantle.clamp_rate` (same physical speed bounds v1 enforced) is
    unchanged."""
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


def boundary_force_reach_rad(spacing_rad: float) -> float:
    """Configured width of the boundary-local force integration band."""
    return BOUNDARY_FORCE_REACH_MULTIPLIER * spacing_rad

# Motion-based boundary classification (see classify_boundary_nodes). A near-boundary node
# whose closing rate is within +-this of zero is a transform (strike-slip) contact; faster
# than this and converging -> convergent, faster and opening -> divergent. 1 cm/yr in real
# m/s -- the same physical threshold `boundary.TRANSFORM_RATE_THRESHOLD` (1 cm/yr in rad/yr)
# expresses for merge_split's own closing-rate test, restated in the m/s units
# `rheology.normal_closing_rate_m_per_s` returns.
BOUNDARY_TRANSFORM_RATE_M_PER_S = 0.01 / SECONDS_PER_YEAR


def boundary_closing_rate_m_per_s(plate, inputs: BoundaryForceInputs) -> np.ndarray:
    """Per own-node relative closing rate against the nearest neighbour plate, in real m/s
    (positive = converging). Thin wrapper over `rheology.normal_closing_rate_m_per_s` so
    `classify_boundary_nodes` here and `quad_tectonics.deform` share one call/convention."""
    from . import rheology

    n = len(inputs.own_points)
    if n == 0:
        return np.zeros(0)
    return rheology.normal_closing_rate_m_per_s(
        np.asarray(plate.omega, dtype=float), inputs.neighbor_omega, inputs.own_points, inputs.direction_to_neighbor
    )


def classify_boundary_nodes(
    plate, neighbours: list, inputs: BoundaryForceInputs, reach_rad: float, cache: BoundarySearchCache | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(convergent, divergent, transform, contested) for this plate's own nodes.

    `contested` is the geometric test (a node currently falls inside some neighbour's
    polygon) -- still the trigger for node deletion / continental retreat, and still widened
    by the bounding-sphere prefilter so a plate that has slid *deep* over a neighbour is
    caught (see below). The other three partition this plate's near-boundary band
    (`dist_to_neighbor <= reach_rad`) by the *sign and size* of the relative closing rate
    (`boundary_closing_rate_m_per_s`), not by geometry: a boundary that is genuinely
    converging builds an orogen before any polygon overlap accumulates, and a transform
    contact is finally its own class rather than being lumped in with divergence. `contested`
    is folded into `convergent` so a deep overrun still thickens / subducts.

    Computed independently of `deform()` since `shift()` runs *before* `deform()` each step
    and needs its own read of the *pre-rotation* configuration to size this step's forces."""
    n = len(inputs.own_points)
    if n == 0 or not neighbours:
        z = np.zeros(n, dtype=bool)
        return z, z.copy(), z.copy(), z.copy()
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
        if cache is None:
            npts, _ = neighbour.all_points_and_elevation()
            if len(npts) == 0:
                continue
            centroid, radius = geometry.bounding_sphere(npts)
        else:
            sphere = cache.sphere(neighbour)
            if sphere is None:
                continue
            centroid, radius = sphere
        consider |= geometry.angular_distance(own_points, centroid) <= radius

    contested = np.zeros(n, dtype=bool)
    if np.any(consider):
        if cache is None:
            consider_points = own_points[consider]
            consider_contested = np.zeros(len(consider_points), dtype=bool)
            for neighbour in neighbours:
                consider_contested |= neighbour.contains_batch(consider_points)
                if np.all(consider_contested):
                    break
            contested[consider] = consider_contested
        else:
            for neighbour in neighbours:
                contested |= cache.inside(plate, neighbour, own_points, consider)
                if np.all(contested[consider]):
                    break

    # Motion partition of the near-boundary band (not the deep-overlap extension above -- an
    # uncontested deep-interior node is just interior, not an active boundary).
    band = inputs.dist_to_neighbor <= reach_rad
    closing = boundary_closing_rate_m_per_s(plate, inputs)
    thr = BOUNDARY_TRANSFORM_RATE_M_PER_S
    convergent = (band & (closing > thr)) | contested
    divergent = band & ~convergent & (closing < -thr)
    transform = band & ~convergent & ~divergent
    return convergent, divergent, transform, contested


def apply_omega_and_rotate(plate, old_points: np.ndarray, omega: np.ndarray, years: float) -> float:
    """Shared tail of `shift_plate`'s real torque-balance path and `Plate.shift`'s
    pinned-omega override (see `World.pinned_omegas`): set `omega`, rotate rigidly by it over
    `years`, and report the same `D` (max node displacement) contract every `shift`
    implementation returns. `old_points` is the plate's own pre-rotation node cloud (the
    caller already has it -- both callers need it for their own zero-node early-out too, so
    it's a parameter rather than refetched here)."""
    plate.set_omega(omega)
    increment = geometry.rotation_matrix_from_omega(plate.omega, years)
    plate.rotate(increment)
    new_points, _ = plate.all_points_and_elevation()
    return float(geometry.angular_distance(old_points, new_points).max())


def continental_contact_mask(plate, inputs: BoundaryForceInputs, contested: np.ndarray) -> np.ndarray:
    """The contested nodes that are a continent-continent contact, the ones contact drag brakes:
    the node's own crust and its nearest neighbour node's are both continental. Node by node,
    not by either plate's majority type (RFC 0009, #447)."""
    return contested & own_is_continental(plate, inputs) & ~inputs.neighbor_is_oceanic


def contact_craton_m(plate, neighbours: list, inputs: BoundaryForceInputs) -> np.ndarray:
    """Per node of `plate`, the craton of the weaker side of its contact: the lesser of its own
    `craton_crust_m` and that at its nearest neighbour node (`inputs.neighbor_node_index`),
    so both plates at one contact are braked alike (`contact_drag_pa_s_per_m`). Zero
    where there is no neighbour."""
    own = plate.collect("craton_crust_m")
    other = np.zeros(len(own))
    for neighbour in neighbours:
        mine = inputs.neighbor_plate_id == neighbour.plate_id
        if np.any(mine):
            other[mine] = neighbour.collect("craton_crust_m")[inputs.neighbor_node_index[mine]]
    return np.minimum(own, other)


def contact_jam_factor(plate, world, neighbours: list, inputs: BoundaryForceInputs, reach_rad: float) -> np.ndarray:
    """Per node of `plate`, the jam factor j of its contact (RFC 0008, #400): its own where it
    is in the jam mask (`jam_mask`), or that of the neighbour's nearest node
    (`inputs.neighbor_node_index`, as `contact_craton_m` reads) where that node is
    continental-coded and forced down, whichever is larger. Both plates at a jammed contact so
    feel one resistance, equal and opposite: the upper plate is never in the jam mask itself,
    but reads the lower plate's j across the contact.

    The neighbour's node is forced down where it is a lower-plate cell, or where the neighbour
    has a slab under this plate: an oceanic-coded nearest node of it at a contact with this
    plate that isn't diverging, the mirror of `slab_attached_continent` read from this side (#447)."""
    n = len(inputs.own_points)
    own = np.where(jam_mask(plate, inputs, forced_down_mask(plate, world, inputs, reach_rad), reach_rad), jam_factor(continental_crust_m(plate, inputs)), 0.0)
    other = np.zeros(n)
    behind_slab = _faces_slab_neighbour(inputs, inputs.neighbor_is_oceanic & _jam_contact(plate, inputs, reach_rad))
    for neighbour in neighbours:
        mine = inputs.neighbor_plate_id == neighbour.plate_id
        if not np.any(mine):
            continue
        idx = inputs.neighbor_node_index[mine]
        continental = ~inputs.neighbor_is_oceanic[mine]
        down = lower_plate_mask(neighbour, world, len(neighbour.cell_keys))[idx] | behind_slab[mine]
        hc = neighbour.collect("crustal_thickness_m")[idx]
        other[mine] = np.where(continental & down, jam_factor(hc), 0.0)
    return np.maximum(own, np.where(_jam_contact(plate, inputs, reach_rad), other, 0.0))


def predict_plate_motion(plate, world, other_plates: list, years: float) -> np.ndarray:
    """Predict this plate's first-step angular velocity without changing its state.

    Generation uses this read-only interface to build terrain from the *initial* tectonic
    regime. The normal stepping path uses the same prediction and then applies it via
    `apply_omega_and_rotate`; this function itself never sets omega or rotates geometry."""
    spacing_rad = line_spacing_rad(world.node_density)
    reach_rad = boundary_force_reach_rad(spacing_rad)

    old_points, _ = plate.all_points_and_elevation()
    if len(old_points) == 0:
        return np.asarray(plate.omega, dtype=float).copy()

    neighbours = plate.get_neighbours(other_plates)
    inputs = gather_boundary_force_inputs(plate, neighbours, spacing_rad, reach_rad)
    _convergent, divergent, _transform, contested = classify_boundary_nodes(plate, neighbours, inputs, reach_rad)
    collision_mask = continental_contact_mask(plate, inputs, contested)
    subducting = subducting_boundary_mask(plate, inputs, reach_rad)

    rho_c = lithosphere.node_crust_density(inputs.own_crust_type_codes, plate.crust_type)
    inertia = lithosphere.moment_of_inertia_tensor(
        inputs.own_points, inputs.own_hc, inputs.own_hm, rho_c, spacing_rad, area_m2=plate.accounting_areas_m2(spacing_rad)
    )

    explicit_torque = slab_pull_torque(plate, inputs, subducting, spacing_rad) + ridge_push_torque(
        plate, inputs, divergent, spacing_rad
    )
    drag_b, drag_k = basal_drag_coefficients(plate, world, spacing_rad)
    drag_k = drag_k + slab_drag_coefficient_matrix(inputs, subducting, spacing_rad)
    # Collision resistance is folded into the same implicit (b, K) drag terms as basal/slab
    # drag, not explicit_torque -- see collision_drag_coefficients' own docstring for why an
    # explicit collision torque here structurally can't compete with basal drag's stiffness.
    jam = contact_jam_factor(plate, world, neighbours, inputs, reach_rad)
    collision_b, collision_k = collision_drag_coefficients(
        plate, inputs, collision_mask, spacing_rad, contact_craton_m(plate, neighbours, inputs), jam
    )
    drag_b = drag_b + collision_b
    drag_k = drag_k + collision_k
    new_omega = integrate_omega(plate, explicit_torque, drag_b, drag_k, inertia, years)
    return np.asarray(new_omega, dtype=float)


def shift_plate(plate, world, other_plates: list, years: float) -> float:
    """`Plate.shift`'s real (non-pinned) implementation: predict omega from the plate's
    current boundary configuration, then rotate rigidly -- see `apply_omega_and_rotate` for
    the shared tail with the pinned-omega override."""
    old_points, _ = plate.all_points_and_elevation()
    if len(old_points) == 0:
        return 0.0
    new_omega = predict_plate_motion(plate, world, other_plates, years)
    return apply_omega_and_rotate(plate, old_points, new_omega, years)


def merge_omega(plate_a, inertia_a: np.ndarray, plate_b, inertia_b: np.ndarray) -> np.ndarray:
    """Angular-momentum-conserving blend for a plate merge (quad_merge.py), rather than a
    naive `(omega_a + omega_b) / 2`: L = I @ omega is additive
    across a fusion, so the combined plate's omega is `(I_a + I_b)^-1 @ (L_a + L_b)`, not a
    plain average of two possibly very differently-massed plates' rates."""
    l_total = lithosphere.angular_momentum(inertia_a, plate_a.omega) + lithosphere.angular_momentum(inertia_b, plate_b.omega)
    combined_omega = lithosphere.omega_from_angular_momentum(inertia_a + inertia_b, l_total)
    return mantle.clamp_rate(combined_omega)

"""Boundary evolution: where a plate's elevation-line nodes are close to another plate's,
classify the local relative motion (convergent/divergent/transform) and apply the
consequences -- uplift, trench deepening, ridge/rift relaxation, and, at each line's two
ends, inserting new nodes (crust created at a divergent boundary) or deleting them (crust
destroyed/folded at a convergent one).

Adjacency is *not* a maintained topological structure -- every step, each plate's nodes are
matched (via k-d tree, restricted to plates a cheap bounding-sphere check can't rule out --
see step_boundaries) against every geometrically-nearby other plate's current nodes. This is
what lets plates evolve independently (see plates.py) while boundaries still behave sensibly:
it's self-healing every step rather than requiring an always-consistent shared-edge structure.

**Convergent boundaries aren't all one effect.** Distance to the boundary alone isn't
enough -- how far an effect *reaches*, and its *shape* with distance, depends on what's
actually converging:

- **Continent-continent collision** (own plate continental, neighbor continental): a broad
  crumple zone, up to COLLISION_RANGE_RAD (400km) -- peaks right at the boundary, decays
  outward, same shape as before, just reaching much farther (e.g. the Himalaya/Tibetan
  Plateau's real deformation belt is comparably wide). On top of that, the same collision
  also drives a second, far weaker effect much farther inland: a gentle intraplate rise
  from FAR_FIELD_COLLISION_INNER_RAD to FAR_FIELD_COLLISION_OUTER_RAD (1000-3000km) -- real
  collisions transmit stress well past their own suture (the Himalayan-Tibetan collision's
  effects reach the Tien Shan/Baikal region at comparable distances; the Arabian-Eurasian
  collision's shortening is distributed across Anatolia the same way; the Laramide orogeny's
  basement uplifts sat ~1000-1500km inland of the contemporary margin). See
  FAR_FIELD_MOUNTAIN_RATE_M_PER_MYR below for why it's additive rather than replacing the
  near-field effect: at any given point at most one of the two bands is actually nonzero
  (400km and 1000km don't overlap), so "additive" only matters in the sense that both use
  `+=` against the same elevation array.
- **Oceanic-under-continental subduction** (own plate continental, neighbor oceanic): a
  volcanic arc, but *offset inland* rather than peaking at the boundary -- real arcs form
  where the subducting slab has descended deep enough to melt, not right at the trench. See
  `_band_intensity`: zero at the boundary, peaking at SUBDUCTION_ARC_INNER_RAD..
  SUBDUCTION_ARC_OUTER_RAD's midpoint (100-300km), zero again past the outer edge. The
  subducting oceanic plate's own trench (the convergent-oceanic branch below) is unaffected
  by this -- it still peaks right at the boundary, same as before.
- **Transform boundaries**: previously had no elevation effect at all. Real strike-slip
  motion can still produce local pressure-ridge relief, so this adds one -- narrower
  (TRANSFORM_RANGE_RAD, 50km) and gentler (TRANSFORM_UPLIFT_RATE_M_PER_MYR) than either
  convergent case, peaking at the boundary like collision/trench do.

**Divergent boundaries aren't uniform either.** Continental rifting stretches and thins the
crust (subsidence) over RIFT_RANGE_RAD (300km) -- much wider than oceanic ridge spreading's
unchanged FAR_THRESHOLD_RAD (~200km) reach. Both still relax exponentially toward their own
target (`_divergent_target`) -- only the reach differs by crust type, same pattern as the
convergent cases above.
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, mantle
from .plates import PLANET_RADIUS_KM, TARGET_LINE_SPACING_RAD, ElevationLine, line_spacing_rad, query_workers

if TYPE_CHECKING:
    from .world import World

FAR_THRESHOLD_RAD = 1.6 * TARGET_LINE_SPACING_RAD
EXTEND_THRESHOLD_RAD = 1.3 * TARGET_LINE_SPACING_RAD
MERGE_THRESHOLD_RAD = 0.4 * TARGET_LINE_SPACING_RAD

TRANSFORM_RATE_THRESHOLD = mantle.cm_per_yr_to_rad_per_yr(1.0)

# Continent-continent collision and a subducting oceanic plate's own trench: peaks right at
# the boundary and decays with distance (the same shape FAR_THRESHOLD_RAD already used),
# just at different reaches -- see COLLISION_RANGE_RAD below for why collision's is wider.
CONVERGENT_MOUNTAIN_RATE_M_PER_MYR = 800.0
CONVERGENT_TRENCH_RATE_M_PER_MYR = 700.0
DIVERGENT_RIDGE_TARGET_M = -1500.0  # new oceanic crust at a mid-ocean ridge
DIVERGENT_RIFT_TARGET_M = -200.0  # new continental crust in a rift valley
DIVERGENT_RELAX_RATE_PER_MYR = 0.5

# Continental rifting stretches and thins the crust over a much wider zone than a plain
# boundary's FAR_THRESHOLD_RAD reach -- the land subsides (relaxes toward
# DIVERGENT_RIFT_TARGET_M) out to RIFT_RANGE_RAD, not just right at the fault line (e.g. the
# East African Rift's actual subsidence zone is comparably wide). Continental-only: oceanic
# ridge spreading is a narrower, already-well-modeled process and keeps FAR_THRESHOLD_RAD's
# reach unchanged.
RIFT_RANGE_KM = 300.0
RIFT_RANGE_RAD = RIFT_RANGE_KM / PLANET_RADIUS_KM

# Continent-continent collision crumples a much broader belt than a plain trench/mountain
# boundary does (e.g. the Himalaya/Tibetan Plateau deformation zone) -- same decay-from-the-
# boundary shape as CONVERGENT_MOUNTAIN_RATE_M_PER_MYR already used, just a wider reach.
COLLISION_RANGE_KM = 400.0
COLLISION_RANGE_RAD = COLLISION_RANGE_KM / PLANET_RADIUS_KM

# Continent-continent collision's far-field effect: real collisions bend and reactivate
# continental interior far beyond their own suture's crumple zone (COLLISION_RANGE_RAD,
# above). Modeled as a second band, entirely separate from and much farther out than the
# near-field one -- zero out to FAR_FIELD_COLLISION_INNER_RAD (leaving the 400-1000km gap
# where the near-field band has already decayed to nothing untouched), ramping linearly from
# full intensity at the inner edge back to zero by FAR_FIELD_COLLISION_OUTER_RAD (see
# _far_field_intensity). Distances tuned against real broad collisional belts: the
# Himalayan-Tibetan collision's far-field reactivation reaches the Tien Shan/Baikal region
# roughly 1500-2500km from the suture; the Arabian-Eurasian (Zagros) collision's shortening
# is distributed across Anatolia at comparable range; the Variscan-Appalachian and Uralian
# orogenies both left deformation belts far wider than their core sutures; the Laramide
# orogeny's basement-cored uplifts sat roughly 1000-1500km inland of the margin. Deliberately
# a much smaller rate than CONVERGENT_MOUNTAIN_RATE_M_PER_MYR -- this is a gentle rise, not
# real mountain-building -- but like the near-field collision effect (and unlike the
# divergent cases' relax-toward-a-target behavior) it adds unconditionally every step it
# applies, so a long-lived collision can still accumulate it into a substantial broad rise
# (e.g. the Colorado Plateau's multi-km Laramide-and-after surface uplift).
FAR_FIELD_COLLISION_INNER_KM = 1000.0
FAR_FIELD_COLLISION_OUTER_KM = 3000.0
FAR_FIELD_COLLISION_INNER_RAD = FAR_FIELD_COLLISION_INNER_KM / PLANET_RADIUS_KM
FAR_FIELD_COLLISION_OUTER_RAD = FAR_FIELD_COLLISION_OUTER_KM / PLANET_RADIUS_KM
FAR_FIELD_MOUNTAIN_RATE_M_PER_MYR = 60.0

# Oceanic-under-continental subduction: the volcanic arc forms *inland* of the trench, not
# at it -- offset by how far the subducting slab needs to descend before it starts melting,
# not by proximity to the boundary itself. Modeled as a band (see _band_intensity), zero at
# the boundary, peaking at the band's midpoint, zero again past the outer edge -- a
# genuinely different shape from every other boundary effect here, which all peak right at
# the boundary and decay outward.
SUBDUCTION_ARC_INNER_KM = 100.0
SUBDUCTION_ARC_OUTER_KM = 300.0
SUBDUCTION_ARC_INNER_RAD = SUBDUCTION_ARC_INNER_KM / PLANET_RADIUS_KM
SUBDUCTION_ARC_OUTER_RAD = SUBDUCTION_ARC_OUTER_KM / PLANET_RADIUS_KM

# Transform (strike-slip) boundaries: real motion here produces at most local pressure-ridge
# relief, not real mountain-building -- narrower reach and a smaller peak rate than either
# convergent case. TRANSFORM_UPLIFT_RATE_M_PER_MYR has no real-world number to port the way
# the others do; picked as a clear fraction of CONVERGENT_MOUNTAIN_RATE_M_PER_MYR ("not as
# big"), a starting point rather than a derived constant.
TRANSFORM_RANGE_KM = 50.0
TRANSFORM_RANGE_RAD = TRANSFORM_RANGE_KM / PLANET_RADIUS_KM
TRANSFORM_UPLIFT_RATE_M_PER_MYR = 200.0

# Widest reach any single boundary effect needs (currently FAR_FIELD_COLLISION_OUTER_RAD) --
# sizes the candidate search (bounding-sphere prescreen and cKDTree query) so a node up to
# that far away is never excluded before the per-effect distance checks above even get to run.
MAX_BOUNDARY_EFFECT_RAD = max(
    FAR_THRESHOLD_RAD,
    COLLISION_RANGE_RAD,
    FAR_FIELD_COLLISION_OUTER_RAD,
    SUBDUCTION_ARC_OUTER_RAD,
    TRANSFORM_RANGE_RAD,
    RIFT_RANGE_RAD,
)

MIN_ELEVATION_M = -11000.0
MAX_ELEVATION_M = 9000.0

# Safety cap on how many nodes a single step can insert at one line end. Not meant to bind
# in practice -- even at MAX_PLATE_RATE with the largest step size the UI offers, the real
# gap is only ever a handful of spacing units (see the comment on _grow_or_shrink_line for
# why a fixed one-node-per-step used to fall far short of that and why it matters). A count
# along a 1D line, not an area, so this scales by ~2x (not ~4x) alongside plates.py's most
# recent resolution doubling -- closing the same physical distance now takes ~2x as many,
# half-as-far-apart nodes.
MAX_EXTEND_NODES_PER_STEP = 400

# FAR_THRESHOLD_RAD/EXTEND_THRESHOLD_RAD/MERGE_THRESHOLD_RAD/MAX_EXTEND_NODES_PER_STEP above
# are the reference (World.node_density == 1.0) values -- kept as bare module constants,
# still the default for every function below, so existing direct calls (tests, any other
# caller without a World in hand) behave exactly as before. step_boundaries itself, which
# does have a World, calls these with `plates.line_spacing_rad(world.node_density)` instead,
# so a world generated at a non-default density keeps growing/merging/reaching at a scale
# consistent with its own actual node spacing rather than the reference one. Unlike
# COLLISION_RANGE_RAD/SUBDUCTION_ARC_*_RAD/TRANSFORM_RANGE_RAD/RIFT_RANGE_RAD (genuine fixed
# physical distances, e.g. "a continent-continent suture is ~400km wide" -- unrelated to how
# many points happen to sample it), these three are explicitly defined as multiples of
# TARGET_LINE_SPACING_RAD, i.e. "roughly N node-spacings," so they're the only boundary
# thresholds that should track a change in spacing at all.
def _far_threshold_rad(spacing_rad: float) -> float:
    return 1.6 * spacing_rad


def _extend_threshold_rad(spacing_rad: float) -> float:
    return 1.3 * spacing_rad


def _merge_threshold_rad(spacing_rad: float) -> float:
    return 0.4 * spacing_rad


def _max_boundary_effect_rad(spacing_rad: float) -> float:
    return max(
        _far_threshold_rad(spacing_rad),
        COLLISION_RANGE_RAD,
        FAR_FIELD_COLLISION_OUTER_RAD,
        SUBDUCTION_ARC_OUTER_RAD,
        TRANSFORM_RANGE_RAD,
        RIFT_RANGE_RAD,
    )


def _max_extend_nodes_per_step(node_density: float) -> int:
    # A 1D count, not an area -- see MAX_EXTEND_NODES_PER_STEP's own comment for why this
    # scales by sqrt(node_density) (~2x at 4x density), not node_density itself.
    return max(1, round(MAX_EXTEND_NODES_PER_STEP * np.sqrt(node_density)))


def _divergent_target(crust_type: str) -> float:
    return DIVERGENT_RIDGE_TARGET_M if crust_type == "oceanic" else DIVERGENT_RIFT_TARGET_M


def _band_intensity(dist: np.ndarray, inner: float, outer: float) -> np.ndarray:
    """Triangular profile: 0 at and outside [inner, outer], peaking at 1.0 at the band's
    midpoint. For an effect (the subduction volcanic arc) that's strongest *offset* from the
    boundary rather than right at it -- every other boundary effect here instead decays
    monotonically from a peak at dist=0, which np.clip(1 - dist / range, 0, 1) already
    expresses directly."""
    mid = (inner + outer) / 2.0
    half_width = (outer - inner) / 2.0
    return np.clip(1.0 - np.abs(dist - mid) / half_width, 0.0, 1.0)


def _far_field_intensity(dist: np.ndarray, inner: float, outer: float) -> np.ndarray:
    """One-sided ramp: 0 below `inner` (unlike the plain decay-from-zero shape every
    near-boundary effect uses, this effect only exists *starting* at `inner`), 1.0 right at
    `inner`, decaying linearly back to 0 by `outer`. Used for the collision far-field band,
    which is deliberately offset inland rather than continuous with the boundary itself."""
    ramp = np.clip(1.0 - (dist - inner) / (outer - inner), 0.0, 1.0)
    return np.where(dist < inner, 0.0, ramp)


def closing_rate(
    points: np.ndarray, self_omega: np.ndarray, neighbor_omega: np.ndarray, neighbor_points: np.ndarray
) -> np.ndarray:
    """Positive = this plate's material is moving toward the neighbor's (convergent) at
    this point; negative = moving apart (divergent). Public because merge_split.py also
    needs it: two plates are already touching along their entire shared boundary by
    construction (plates.py's tiling has no gaps), so proximity alone can't distinguish an
    actively-colliding pair from any other pair of neighbors -- it has to check motion."""
    v_self = np.cross(self_omega, points)
    v_neighbor = np.cross(neighbor_omega, points)
    normal_dir = neighbor_points - points
    norm = np.linalg.norm(normal_dir, axis=-1, keepdims=True)
    safe_norm = np.where(norm < 1e-12, 1.0, norm)
    normal_dir = normal_dir / safe_norm
    return np.sum((v_self - v_neighbor) * normal_dir, axis=-1)


def _grow_or_shrink_line(
    line: ElevationLine,
    dist: np.ndarray,
    closing: np.ndarray,
    crust_type: str,
    spacing_rad: float = TARGET_LINE_SPACING_RAD,
    extend_threshold_rad: float = EXTEND_THRESHOLD_RAD,
    merge_threshold_rad: float = MERGE_THRESHOLD_RAD,
    max_extend_nodes: int = MAX_EXTEND_NODES_PER_STEP,
    # Only matters to a caller exercising multi-node shrink specifically (see below) --
    # every real caller (step_boundaries) always passes its own actual years, this default
    # just keeps existing single-node-shrink unit tests (which don't care about this) working
    # unchanged.
    years: float = 1_000_000.0,
) -> ElevationLine:
    """Grow or shrink a line's two ends based on this step's boundary classification there.

    Growth inserts as many nodes as it takes to actually close the gap (`dist`), not just
    one -- at a small step size the two are the same thing, but at a large `years` step (the
    UI offers up to 10 Myr/step) a fast-diverging boundary can open by many spacing units in
    a single step, and inserting only one node per step falls further and further behind.
    That leftover, perpetually-reopening gap looked to gaps.py's periodic gap-filling like a
    genuinely new, unclosable gap and kept spawning fresh micro-plates at the same busy
    boundary every interval -- fixing the actual growth rate here is what stops that at the
    source, rather than only reacting to its symptom in gaps.py.

    Shrink is the mirror case, and needed the same fix for the same reason: `closing[i] *
    years` is how far this point's material actually moves toward its cross-plate neighbor
    this step (both in radians, same units as `dist`/`spacing_rad` -- see closing_rate's own
    docstring), so a fast enough closing rate, or a large enough `years` step, can consume
    many spacing units of subducting crust in one step, not just one. Removing only one node
    per step regardless fell behind exactly the way one-node-per-step growth used to: the
    un-consumed remainder of that motion is still there next step, i.e. exactly the "nodes of
    the subducted plate overlap nodes of the remaining plate" case -- consuming the full
    `closing[i] * years` worth of nodes per step (not just one) is what keeps deletion paced
    with however fast the convergence actually is, so the overlap this guards against never
    gets the chance to accumulate in the first place.

    Deliberately *not* also widening the shrink condition itself (`dist[i] < merge_
    threshold_rad and closing[i] > TRANSFORM_RATE_THRESHOLD`) to fire on raw proximity alone,
    ignoring `closing[i]`'s sign, as a second line of defense once two nodes are already
    essentially colocated: measured directly against a freshly generated world (see
    generate_world), neighboring plates' own nearest cross-plate nodes already start as close
    as ~0.07 * spacing_rad at generation time purely from plates.py's gapless tiling, well
    inside any distance tight enough to mean "these are colocated, not just newly touching."
    A proximity-only override can't tell that freshly-divergent case apart from genuine
    subduction overlap, so it would delete nodes off a boundary that's actually supposed to
    grow. `closing[i]`'s sign is the only signal here that can, so the fix stays scoped to
    sizing *how much* gets removed once it already reads convergent, not to overriding it.

    channel_depth/channel_width/lake_depth/glacier_depth/silt_depth/is_volcano/volcano_active_
    years_remaining/soil_depth/soil_mineral_content/soil_organic_content/coal_deposit_m/
    oil_gas_deposit_m/mineral_deposit_m ride along: a surviving node keeps its own prior value
    (sliced the same way theta/elevation are), a newly-inserted node (brand new crust) starts
    at 0/False -- no history to carry, the same reasoning plates.ElevationLine's own defaulting
    already uses for a call site that doesn't pass them at all.

    Every persistent field is threaded through the same append/insert/slice as theta/
    elevation here (not dataclasses.replace -- there's no "unchanged" value to copy when the
    array itself is being reshaped), keyed in one dict rather than one local variable per
    field, so a future persistent field only needs to be added to `persistent_fields` below,
    not duplicated across every branch of this function -- see plates.ElevationLine's own
    docstring for the bug class this avoids."""
    theta = line.theta.copy()
    elevation = line.elevation.copy()
    persistent_fields = {
        "channel_depth": line.channel_depth.copy(),
        "channel_width": line.channel_width.copy(),
        "lake_depth": line.lake_depth.copy(),
        "glacier_depth": line.glacier_depth.copy(),
        "silt_depth": line.silt_depth.copy(),
        "is_volcano": line.is_volcano.copy(),
        "volcano_active_years_remaining": line.volcano_active_years_remaining.copy(),
        "soil_depth": line.soil_depth.copy(),
        "soil_mineral_content": line.soil_mineral_content.copy(),
        "soil_organic_content": line.soil_organic_content.copy(),
        "coal_deposit_m": line.coal_deposit_m.copy(),
        "oil_gas_deposit_m": line.oil_gas_deposit_m.copy(),
        "mineral_deposit_m": line.mineral_deposit_m.copy(),
    }
    if len(theta) == 0:
        return ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)

    dtheta = spacing_rad / max(np.cos(line.phi), 1e-3)
    target = _divergent_target(crust_type)

    # High end first so the low-end index (0) is unaffected by any change made here.
    if dist[-1] > extend_threshold_rad and closing[-1] < -TRANSFORM_RATE_THRESHOLD:
        n_new = min(max(int(dist[-1] / spacing_rad), 1), max_extend_nodes)
        new_theta = theta[-1] + dtheta * np.arange(1, n_new + 1)
        theta = np.append(theta, new_theta)
        elevation = np.append(elevation, np.full(n_new, target))
        for name, values in persistent_fields.items():
            fill = np.zeros(n_new, dtype=values.dtype)
            persistent_fields[name] = np.append(values, fill)
    elif dist[-1] < merge_threshold_rad and closing[-1] > TRANSFORM_RATE_THRESHOLD and len(theta) > 1:
        n_remove = min(max(int(closing[-1] * years / spacing_rad), 1), max_extend_nodes, len(theta) - 1)
        theta = theta[:-n_remove]
        elevation = elevation[:-n_remove]
        persistent_fields = {name: values[:-n_remove] for name, values in persistent_fields.items()}

    if len(theta) == 0:
        return ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)

    if dist[0] > extend_threshold_rad and closing[0] < -TRANSFORM_RATE_THRESHOLD:
        n_new = min(max(int(dist[0] / spacing_rad), 1), max_extend_nodes)
        new_theta = theta[0] - dtheta * np.arange(n_new, 0, -1)
        theta = np.insert(theta, 0, new_theta)
        elevation = np.insert(elevation, 0, np.full(n_new, target))
        for name, values in persistent_fields.items():
            fill = np.zeros(n_new, dtype=values.dtype)
            persistent_fields[name] = np.insert(values, 0, fill)
    elif dist[0] < merge_threshold_rad and closing[0] > TRANSFORM_RATE_THRESHOLD and len(theta) > 1:
        n_remove = min(max(int(closing[0] * years / spacing_rad), 1), max_extend_nodes, len(theta) - 1)
        theta = theta[n_remove:]
        elevation = elevation[n_remove:]
        persistent_fields = {name: values[n_remove:] for name, values in persistent_fields.items()}

    return ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)


def step_boundaries(world: World, years: float) -> None:
    """Mutates world.plates in place: applies boundary elevation deltas and grows/shrinks
    each line's two ends. Call after rotation, using the plates' new positions."""
    years_myr = years / 1e6

    # Scaled to this world's own node_density (see the block above _divergent_target for
    # why only these particular thresholds/cap need it) rather than reading the bare
    # reference-density module constants directly.
    spacing_rad = line_spacing_rad(world.node_density)
    far_threshold_rad = _far_threshold_rad(spacing_rad)
    extend_threshold_rad = _extend_threshold_rad(spacing_rad)
    merge_threshold_rad = _merge_threshold_rad(spacing_rad)
    max_boundary_effect_rad = _max_boundary_effect_rad(spacing_rad)
    max_extend_nodes = _max_extend_nodes_per_step(world.node_density)

    plate_by_id = {plate.plate_id: plate for plate in world.plates}

    plate_points: dict[int, np.ndarray] = {}
    for plate in world.plates:
        pieces = [line.world_xyz(plate.frame) for line in plate.lines]
        plate_points[plate.plate_id] = (
            np.concatenate(pieces, axis=0) if pieces else np.zeros((0, 3))
        )

    # Cheap bounding sphere per plate (see geometry.bounding_sphere): used below so each
    # plate's tree/query is built from only the other plates a bounding-sphere check can't
    # rule out, rather than literally every other plate concatenated -- and the query below
    # runs once per plate (all its lines' points at once) rather than once per line, which
    # together are what dominated step time once plates carry thousands of nodes each (see
    # docs/simulation-model.md's resolution note).
    spheres = {pid: geometry.bounding_sphere(pts) for pid, pts in plate_points.items() if len(pts) > 0}

    for plate in world.plates:
        own_points = plate_points[plate.plate_id]
        if not plate.lines or len(own_points) == 0 or plate.plate_id not in spheres:
            continue

        ca, ra = spheres[plate.plate_id]
        other_points_list = []
        other_owner_list = []
        for other in world.plates:
            if other.plate_id == plate.plate_id or other.plate_id not in spheres:
                continue
            cb, rb = spheres[other.plate_id]
            centroid_dist = float(geometry.angular_distance(ca, cb))
            if centroid_dist - ra - rb > max_boundary_effect_rad:
                continue  # no point of "other" can possibly be within max_boundary_effect_rad
            other_points_list.append(plate_points[other.plate_id])
            other_owner_list.append(np.full(len(plate_points[other.plate_id]), other.plate_id))
        if not other_points_list:
            continue
        other_points = np.concatenate(other_points_list, axis=0)
        other_owner = np.concatenate(other_owner_list, axis=0)
        # balanced_tree=False/compact_nodes=False trade a slightly slower query for a much
        # faster build -- the right tradeoff here since this tree is built fresh per plate,
        # per step, and queried exactly once (never reused); still an exact nearest-neighbor
        # search either way, just a different internal tree layout, so results are unchanged
        # (confirmed via the stress test suite's determinism checks). Profiled at this
        # simulation's realistic plate-count/node-count scale (~100k-node trees): roughly
        # half the combined build+query time. query_workers (see plates.py) parallelizes the
        # query only when own_points is large enough for that to pay for itself -- a plate's
        # own node count ranges from a handful (freshly spawned) to tens of thousands
        # (established), and spinning up a thread pool for a 5-point query is pure overhead.
        tree = cKDTree(other_points, balanced_tree=False, compact_nodes=False)

        dist_all, idx_all = tree.query(own_points, workers=query_workers(len(own_points)))
        neighbor_owner_all = other_owner[idx_all]
        neighbor_points_all = other_points[idx_all]
        neighbor_omega_all = np.array([plate_by_id[o].omega for o in neighbor_owner_all])
        closing_all = closing_rate(own_points, plate.omega, neighbor_omega_all, neighbor_points_all)
        # Crust type of each point's own single nearest cross-plate neighbor -- same
        # per-point python lookup pattern as neighbor_omega_all above (small plate count,
        # not a hot loop). Distinguishes oceanic-under-continental subduction (-> a volcanic
        # arc band, see subduction_all below) from continent-continent collision (-> a
        # broad crumple zone) -- both are "convergent_all", but shaped very differently.
        neighbor_crust_all = np.array([plate_by_id[o].crust_type for o in neighbor_owner_all])
        neighbor_is_oceanic_all = neighbor_crust_all == "oceanic"

        # default_intensity: the original plain decay-from-the-boundary shape, still used
        # for divergent boundaries and a subducting oceanic plate's own trench (neither
        # changed here) -- FAR_THRESHOLD_RAD's ~200km reach, unchanged.
        default_intensity_all = np.clip(1.0 - dist_all / far_threshold_rad, 0.0, 1.0)
        # collision_intensity: same shape, but reaching COLLISION_RANGE_RAD (400km) instead
        # -- a continent-continent suture crumples a much wider belt than a plain boundary.
        collision_intensity_all = np.clip(1.0 - dist_all / COLLISION_RANGE_RAD, 0.0, 1.0)
        # far_field_intensity: the collision's second, much weaker and much farther-reaching
        # band -- zero out to FAR_FIELD_COLLISION_INNER_RAD, then ramping down to zero by
        # FAR_FIELD_COLLISION_OUTER_RAD (see _far_field_intensity and the comment on the
        # constants themselves for why 1000-3000km).
        far_field_intensity_all = _far_field_intensity(
            dist_all, FAR_FIELD_COLLISION_INNER_RAD, FAR_FIELD_COLLISION_OUTER_RAD
        )
        # arc_intensity: the one non-monotonic shape here -- see _band_intensity. Zero right
        # at the boundary, peaks inland at the volcanic arc's typical offset, zero again past
        # SUBDUCTION_ARC_OUTER_RAD.
        arc_intensity_all = _band_intensity(dist_all, SUBDUCTION_ARC_INNER_RAD, SUBDUCTION_ARC_OUTER_RAD)
        # transform_intensity: plain decay again, but TRANSFORM_RANGE_RAD's much shorter
        # ~50km reach -- real strike-slip motion doesn't build real mountains.
        transform_intensity_all = np.clip(1.0 - dist_all / TRANSFORM_RANGE_RAD, 0.0, 1.0)
        # rift_intensity: same plain decay shape as default_intensity, but reaching
        # RIFT_RANGE_RAD (300km) instead -- continental rifting stretches and thins the
        # crust over a much wider zone than oceanic ridge spreading (which still uses
        # default_intensity/FAR_THRESHOLD_RAD, unchanged).
        rift_intensity_all = np.clip(1.0 - dist_all / RIFT_RANGE_RAD, 0.0, 1.0)

        # Classification is by *rate* only (not distance) -- MAX_BOUNDARY_EFFECT_RAD here
        # just bounds the candidate search; each intensity array above already zeroes itself
        # out past its own specific (narrower) reach, so no further distance masking is
        # needed below. divergent_all's own gate uses RIFT_RANGE_RAD (the wider of its two
        # cases) for the same reason.
        convergent_all = (dist_all < max_boundary_effect_rad) & (closing_all > TRANSFORM_RATE_THRESHOLD)
        divergent_all = (dist_all < RIFT_RANGE_RAD) & (closing_all < -TRANSFORM_RATE_THRESHOLD)
        transform_all = (dist_all < TRANSFORM_RANGE_RAD) & (np.abs(closing_all) <= TRANSFORM_RATE_THRESHOLD)
        # Only meaningful when this plate itself is continental (see the per-line loop
        # below) -- computed for every point regardless, matching convergent_all/
        # divergent_all's own unconditional-computation style.
        subduction_all = convergent_all & neighbor_is_oceanic_all
        collision_all = convergent_all & ~neighbor_is_oceanic_all

        target = _divergent_target(plate.crust_type)
        relax_factor = 1.0 - np.exp(-DIVERGENT_RELAX_RATE_PER_MYR * years_myr)

        new_lines = []
        offset = 0
        for line in plate.lines:
            n = len(line.theta)
            dist = dist_all[offset : offset + n]
            closing = closing_all[offset : offset + n]
            default_intensity = default_intensity_all[offset : offset + n]
            collision_intensity = collision_intensity_all[offset : offset + n]
            far_field_intensity = far_field_intensity_all[offset : offset + n]
            arc_intensity = arc_intensity_all[offset : offset + n]
            transform_intensity = transform_intensity_all[offset : offset + n]
            rift_intensity = rift_intensity_all[offset : offset + n]
            convergent = convergent_all[offset : offset + n]
            divergent = divergent_all[offset : offset + n]
            transform = transform_all[offset : offset + n]
            subduction = subduction_all[offset : offset + n]
            collision = collision_all[offset : offset + n]
            offset += n

            elevation = line.elevation.copy()
            if plate.crust_type == "continental":
                elevation[subduction] += CONVERGENT_MOUNTAIN_RATE_M_PER_MYR * years_myr * arc_intensity[subduction]
                elevation[collision] += CONVERGENT_MOUNTAIN_RATE_M_PER_MYR * years_myr * collision_intensity[collision]
                elevation[collision] += FAR_FIELD_MOUNTAIN_RATE_M_PER_MYR * years_myr * far_field_intensity[collision]
            else:
                elevation[convergent] -= CONVERGENT_TRENCH_RATE_M_PER_MYR * years_myr * default_intensity[convergent]

            elevation[transform] += TRANSFORM_UPLIFT_RATE_M_PER_MYR * years_myr * transform_intensity[transform]

            # Continental rifting reaches RIFT_RANGE_RAD (300km); oceanic ridge spreading
            # keeps the original default_intensity/FAR_THRESHOLD_RAD (~200km) reach.
            divergent_intensity = rift_intensity if plate.crust_type == "continental" else default_intensity
            elevation[divergent] += (target - elevation[divergent]) * relax_factor * divergent_intensity[divergent]

            elevation = np.clip(elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)
            # theta unchanged -- dataclasses.replace copies every other field (channel_width
            # included) from the existing line automatically.
            updated_line = dataclasses.replace(line, elevation=elevation)
            grown_line = _grow_or_shrink_line(
                updated_line,
                dist,
                closing,
                plate.crust_type,
                spacing_rad,
                extend_threshold_rad,
                merge_threshold_rad,
                max_extend_nodes,
                years=years,
            )
            if len(grown_line.theta) > 0:
                new_lines.append(grown_line)

        plate.lines = new_lines

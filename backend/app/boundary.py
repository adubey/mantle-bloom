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
  Plateau's real deformation belt is comparably wide).
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

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, mantle
from .plates import PLANET_RADIUS_KM, TARGET_LINE_SPACING_RAD, ElevationLine

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

# Widest reach any single boundary effect needs (currently COLLISION_RANGE_RAD) -- sizes the
# candidate search (bounding-sphere prescreen and cKDTree query) so a node up to that far
# away is never excluded before the per-effect distance checks above even get to run.
MAX_BOUNDARY_EFFECT_RAD = max(
    FAR_THRESHOLD_RAD, COLLISION_RANGE_RAD, SUBDUCTION_ARC_OUTER_RAD, TRANSFORM_RANGE_RAD, RIFT_RANGE_RAD
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

    channel_depth/lake_depth/glacier_depth ride along: a surviving node keeps its own prior
    value (sliced the same way theta/elevation are), a newly-inserted node (brand new crust)
    starts at 0 -- no history to carry, the same reasoning plates.ElevationLine's own
    defaulting already uses for a call site that doesn't pass them at all."""
    theta = line.theta.copy()
    elevation = line.elevation.copy()
    channel_depth = line.channel_depth.copy()
    lake_depth = line.lake_depth.copy()
    glacier_depth = line.glacier_depth.copy()
    if len(theta) == 0:
        return ElevationLine(
            phi=line.phi, theta=theta, elevation=elevation, channel_depth=channel_depth, lake_depth=lake_depth, glacier_depth=glacier_depth
        )

    dtheta = TARGET_LINE_SPACING_RAD / max(np.cos(line.phi), 1e-3)
    target = _divergent_target(crust_type)

    # High end first so the low-end index (0) is unaffected by any change made here.
    if dist[-1] > EXTEND_THRESHOLD_RAD and closing[-1] < -TRANSFORM_RATE_THRESHOLD:
        n_new = min(max(int(dist[-1] / TARGET_LINE_SPACING_RAD), 1), MAX_EXTEND_NODES_PER_STEP)
        new_theta = theta[-1] + dtheta * np.arange(1, n_new + 1)
        theta = np.append(theta, new_theta)
        elevation = np.append(elevation, np.full(n_new, target))
        channel_depth = np.append(channel_depth, np.zeros(n_new))
        lake_depth = np.append(lake_depth, np.zeros(n_new))
        glacier_depth = np.append(glacier_depth, np.zeros(n_new))
    elif dist[-1] < MERGE_THRESHOLD_RAD and closing[-1] > TRANSFORM_RATE_THRESHOLD and len(theta) > 1:
        theta = theta[:-1]
        elevation = elevation[:-1]
        channel_depth = channel_depth[:-1]
        lake_depth = lake_depth[:-1]
        glacier_depth = glacier_depth[:-1]

    if len(theta) == 0:
        return ElevationLine(
            phi=line.phi, theta=theta, elevation=elevation, channel_depth=channel_depth, lake_depth=lake_depth, glacier_depth=glacier_depth
        )

    if dist[0] > EXTEND_THRESHOLD_RAD and closing[0] < -TRANSFORM_RATE_THRESHOLD:
        n_new = min(max(int(dist[0] / TARGET_LINE_SPACING_RAD), 1), MAX_EXTEND_NODES_PER_STEP)
        new_theta = theta[0] - dtheta * np.arange(n_new, 0, -1)
        theta = np.insert(theta, 0, new_theta)
        elevation = np.insert(elevation, 0, np.full(n_new, target))
        channel_depth = np.insert(channel_depth, 0, np.zeros(n_new))
        lake_depth = np.insert(lake_depth, 0, np.zeros(n_new))
        glacier_depth = np.insert(glacier_depth, 0, np.zeros(n_new))
    elif dist[0] < MERGE_THRESHOLD_RAD and closing[0] > TRANSFORM_RATE_THRESHOLD and len(theta) > 1:
        theta = theta[1:]
        elevation = elevation[1:]
        channel_depth = channel_depth[1:]
        lake_depth = lake_depth[1:]
        glacier_depth = glacier_depth[1:]

    return ElevationLine(
        phi=line.phi, theta=theta, elevation=elevation, channel_depth=channel_depth, lake_depth=lake_depth, glacier_depth=glacier_depth
    )


def step_boundaries(world: World, years: float) -> None:
    """Mutates world.plates in place: applies boundary elevation deltas and grows/shrinks
    each line's two ends. Call after rotation, using the plates' new positions."""
    years_myr = years / 1e6

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
            if centroid_dist - ra - rb > MAX_BOUNDARY_EFFECT_RAD:
                continue  # no point of "other" can possibly be within MAX_BOUNDARY_EFFECT_RAD
            other_points_list.append(plate_points[other.plate_id])
            other_owner_list.append(np.full(len(plate_points[other.plate_id]), other.plate_id))
        if not other_points_list:
            continue
        other_points = np.concatenate(other_points_list, axis=0)
        other_owner = np.concatenate(other_owner_list, axis=0)
        tree = cKDTree(other_points)

        dist_all, idx_all = tree.query(own_points)
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
        default_intensity_all = np.clip(1.0 - dist_all / FAR_THRESHOLD_RAD, 0.0, 1.0)
        # collision_intensity: same shape, but reaching COLLISION_RANGE_RAD (400km) instead
        # -- a continent-continent suture crumples a much wider belt than a plain boundary.
        collision_intensity_all = np.clip(1.0 - dist_all / COLLISION_RANGE_RAD, 0.0, 1.0)
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
        convergent_all = (dist_all < MAX_BOUNDARY_EFFECT_RAD) & (closing_all > TRANSFORM_RATE_THRESHOLD)
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
            else:
                elevation[convergent] -= CONVERGENT_TRENCH_RATE_M_PER_MYR * years_myr * default_intensity[convergent]

            elevation[transform] += TRANSFORM_UPLIFT_RATE_M_PER_MYR * years_myr * transform_intensity[transform]

            # Continental rifting reaches RIFT_RANGE_RAD (300km); oceanic ridge spreading
            # keeps the original default_intensity/FAR_THRESHOLD_RAD (~200km) reach.
            divergent_intensity = rift_intensity if plate.crust_type == "continental" else default_intensity
            elevation[divergent] += (target - elevation[divergent]) * relax_factor * divergent_intensity[divergent]

            elevation = np.clip(elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)
            updated_line = ElevationLine(
                phi=line.phi,
                theta=line.theta,
                elevation=elevation,
                channel_depth=line.channel_depth,
                lake_depth=line.lake_depth,
                glacier_depth=line.glacier_depth,
            )
            grown_line = _grow_or_shrink_line(updated_line, dist, closing, plate.crust_type)
            if len(grown_line.theta) > 0:
                new_lines.append(grown_line)

        plate.lines = new_lines

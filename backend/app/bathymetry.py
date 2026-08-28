"""Bathymetry: the shape of the sea floor between the shelf break and the abyssal plain.

Isostasy (see lithosphere.py) sets each submerged column's own depth directly from its
crustal/mantle-lithosphere thickness, and per-plate generation seeds those thicknesses
purely from a plate's crust type -- which leaves two artefacts a real ocean floor doesn't
have, both of which `shape_initial_bathymetry` corrects **once, at generation** (nothing
here touches the ongoing simulation):

- Submerged *continental* crust reads as one uniform bright shelf everywhere, however far
  from land it sits -- there's no mechanism drowning a large continental block's interior to
  basin depth the way real thinned/drowned crust (Zealandia, submerged plateaus) subsides.
  `_subside_offshore_continental_crust` thins it toward the abyssal reference depth with
  distance from the nearest coastline.
- Every continent/ocean *plate boundary* is a vertical cliff -- full-thickness continental
  crust straight against thin oceanic crust, nothing between. Real margins ramp shelf ->
  slope -> rise -> abyssal plain across a ~100-300km transition zone.
  `_smooth_continental_margins` grades the columns across that contact.

`SHELF_RANGE_KM`/`SHELF_RANGE_RAD` (the shelf width geology.py's oil & gas formation keys
off, to tell shallow shelf water from open ocean) also lives here.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import lithosphere
from .elevation_lines import PLANET_RADIUS_KM

if TYPE_CHECKING:
    from .plates import Plate

# Real continental shelves are shallow and comparatively narrow before the "shelf break"
# drops off toward deep water -- SHELF_RANGE_KM has no exact figure to port (none was given),
# picked as a reasonable shelf width.
SHELF_RANGE_KM = 200.0
SHELF_RANGE_RAD = SHELF_RANGE_KM / PLANET_RADIUS_KM

# --- Offshore continental subsidence (_subside_offshore_continental_crust) ----------------
# Distance from the nearest coastline over which submerged continental crust grades from
# shelf depth to full basin depth. Within OFFSHORE_SHELF_KM it stays a shelf; past
# OFFSHORE_ABYSSAL_KM it's drawn all the way to the abyssal reference depth; between, a
# smoothstep. ~200 -> ~1400km spans a real shelf-plus-slope-plus-rise apron.
OFFSHORE_SHELF_KM = SHELF_RANGE_KM
OFFSHORE_ABYSSAL_KM = 1400.0
OFFSHORE_SHELF_DEPTH_M = -200.0  # target depth for still-on-the-shelf submerged crust
# Full basin depth submerged continental crust is drawn toward far from land: the same
# abyssal reference an aged oceanic column floats at, so a drowned continental interior and
# the open ocean around it read as the same deep sea floor rather than plateau-vs-basin.
ABYSSAL_REFERENCE_DEPTH_M = float(
    lithosphere.isostatic_elevation(
        np.array([lithosphere.REFERENCE_HC_OCEANIC_M]),
        np.array([lithosphere.REFERENCE_HM_OCEANIC_M]),
        lithosphere.RHO_OCEANIC_CRUST,
    )[0]
)

# --- Continent/ocean margin grading (_smooth_continental_margins) -------------------------

# How far either side of a continent/ocean plate boundary the sea floor is graded from a
# shelf-thickness column toward an abyssal one. Real continent-ocean transition zones
# (stretched, thinned continental crust plus the sediment-draped slope and rise seaward of
# the shelf break) run ~100-300km; 400 is the outer end of that, wide enough that the ramp
# reads as a gradual slope at render resolution rather than a one-cell step.
MARGIN_TRANSITION_KM = 400.0
MARGIN_TRANSITION_RAD = MARGIN_TRANSITION_KM / PLANET_RADIUS_KM

# The grading is an iterated neighbour-averaging (a restricted Laplacian smooth -- see below).
# The band is only a few nodes thick and each node's per-iteration step is capped by its
# `accept` weight, so it settles fast: this count is scaled by (band width / node spacing)
# and clamped to a modest range -- more iterations past that don't visibly change the result.
MARGIN_RELAX_ITERS_MIN = 15
MARGIN_RELAX_ITERS_MAX = 60

# A submerged node's grading weight ramps from 0 at COAST_GUARD_DEPTH_M (just below sea
# level) to 1 at MARGIN_FULL_DEPTH_M, so coastlines and the shallowest shelf hold their
# position while everything from the shelf break down grades fully into a slope.
COAST_GUARD_DEPTH_M = 60.0
MARGIN_FULL_DEPTH_M = 250.0

# A node is a neighbour of another if it's within this multiple of the local node spacing.
# Wider than one lattice step on purpose: nodes on opposite sides of a plate boundary belong
# to two independent lattices that only meet along a thin seam, and their nearest cross-seam
# links sit noticeably farther apart than in-plate neighbours -- this reach keeps those links
# in the graph, which is what lets the smoothing actually cross the contact.
NEIGHBOUR_SPACING_MULTIPLE = 2.5
_NEIGHBOUR_QUERY_K = 24  # candidate neighbours pulled per node before the distance cut


def _chord(angle_rad: float) -> float:
    """Straight-line distance between two unit vectors `angle_rad` apart -- cKDTree over
    world unit vectors measures Euclidean chord length, not arc length."""
    return 2.0 * np.sin(angle_rad / 2.0)


def _arc_km(chord: np.ndarray) -> np.ndarray:
    """Chord length between unit vectors -> great-circle distance in km."""
    return 2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0)) * PLANET_RADIUS_KM


def shape_initial_bathymetry(plates: list["Plate"]) -> None:
    """The full generation-time sea-floor fix-up, in order: first drown the interiors of
    submerged continental crust toward basin depth by distance from land, then grade the
    (now smaller) continent/ocean plate-boundary steps into slopes. Both are one-shot -- the
    ongoing simulation is untouched. See this module's docstring."""
    _subside_offshore_continental_crust(plates)
    _smooth_continental_margins(plates)


def _gather_columns(plates: list["Plate"]):
    """(points, plates_in_order, per-plate node counts, and node-aligned is_continental /
    Hc / Hm / elevation arrays) for every plate with nodes -- the common preamble both passes
    below need."""
    from . import plates as plates_mod

    points, plates_in_order = plates_mod.gather_node_positions(plates)
    counts = [p.node_count() for p in plates_in_order]
    is_continental = np.concatenate(
        [np.full(c, p.crust_type == "continental", dtype=bool) for p, c in zip(plates_in_order, counts)]
    ) if counts else np.zeros(0, dtype=bool)
    hc = np.concatenate([np.asarray(p.collect("crustal_thickness_m"), dtype=float) for p in plates_in_order]) if counts else np.zeros(0)
    hm = np.concatenate([np.asarray(p.collect("mantle_lithosphere_thickness_m"), dtype=float) for p in plates_in_order]) if counts else np.zeros(0)
    elevation = np.concatenate([np.asarray(p.collect("elevation"), dtype=float) for p in plates_in_order]) if counts else np.zeros(0)
    return points, plates_in_order, counts, is_continental, hc, hm, elevation


def _write_back(plates_in_order, counts, active, new_hc, new_hm) -> None:
    """Push per-node Hc/Hm back onto each plate (only those with a changed node) and re-sync
    its cached `elevation`."""
    offset = 0
    for plate, count in zip(plates_in_order, counts):
        sl = slice(offset, offset + count)
        offset += count
        if not np.any(active[sl]):
            continue
        plate.set_fields_on_plate(crustal_thickness_m=new_hc[sl], mantle_lithosphere_thickness_m=new_hm[sl])
        lithosphere.sync_plate_elevation(plate)


def _subside_offshore_continental_crust(plates: list["Plate"]) -> None:
    """Thin submerged continental crust toward `ABYSSAL_REFERENCE_DEPTH_M` with distance from
    the nearest coastline, so a large drowned continental block's interior reads as a deep
    basin (matching the open ocean around it) rather than a uniform bright shelf.

    - Target depth by distance to nearest land node (`elevation > 0`, any plate -- a
      geographic coastline question, not a plate-boundary one): shelf within
      `OFFSHORE_SHELF_KM`, full abyssal past `OFFSHORE_ABYSSAL_KM`, smoothstep between.
    - Applied only to submerged *continental* nodes, and only ever *downward* (this drowns
      crust; it never lifts a node the generation noise already put deep). Weighted in with
      depth over `COAST_GUARD_DEPTH_M -> MARGIN_FULL_DEPTH_M` so coastlines and the shallow
      shelf hold their position.
    - Realised by thinning `Hc` (via `lithosphere.crustal_thickness_for_submerged_elevation`)
      so `elevation` stays a faithful isostatic readout -- physically, hyper-extended /
      attenuated continental crust, which is exactly what real drowned continental interiors
      are. `Hm` is left as-is (a scope simplification: the mantle lid isn't separately
      attenuated at generation).

    Generation-time only.
    """
    points, plates_in_order, counts, is_continental, hc, hm, elevation = _gather_columns(plates)
    if len(points) < 3:
        return

    land = elevation > 0.0
    submerged_continental = is_continental & ~land
    if not np.any(land) or not np.any(submerged_continental):
        return

    from . import plates as plates_mod

    dist_to_land_chord, _ = cKDTree(points[land]).query(points, workers=plates_mod.query_workers(len(points)))
    dist_km = _arc_km(dist_to_land_chord)

    t = np.clip((dist_km - OFFSHORE_SHELF_KM) / (OFFSHORE_ABYSSAL_KM - OFFSHORE_SHELF_KM), 0.0, 1.0)
    t = t * t * (3.0 - 2.0 * t)  # smoothstep
    target_z = OFFSHORE_SHELF_DEPTH_M + t * (ABYSSAL_REFERENCE_DEPTH_M - OFFSHORE_SHELF_DEPTH_M)

    depth_w = np.clip((-elevation - COAST_GUARD_DEPTH_M) / (MARGIN_FULL_DEPTH_M - COAST_GUARD_DEPTH_M), 0.0, 1.0)
    # Downward only: weight is 0 wherever the target is not deeper than the node already is.
    deepening = submerged_continental & (target_z < elevation)
    w = np.where(deepening, depth_w, 0.0)
    if not np.any(w > 0.0):
        return

    new_z = elevation + w * (target_z - elevation)
    rho_c = lithosphere.RHO_CONTINENTAL_CRUST
    new_hc = hc.copy()
    subsided_hc = lithosphere.crustal_thickness_for_submerged_elevation(new_z, hm, rho_c)
    # Only ever thin, never thicken, and never below the floor Hc integration allows.
    new_hc = np.where(w > 0.0, np.clip(np.minimum(subsided_hc, hc), lithosphere.MIN_CRUSTAL_THICKNESS_M, None), hc)

    _write_back(plates_in_order, counts, w > 0.0, new_hc, hm)


def _smooth_continental_margins(plates: list["Plate"]) -> None:
    """Grade each node's Hc/Hm laterally across every continent/ocean plate boundary so a
    freshly generated world's sea floor ramps shelf -> slope -> abyssal plain over
    `MARGIN_TRANSITION_KM` instead of dropping vertically at the plate edge (which is what
    per-plate generation, seeding every node's column purely from its own plate's crust type,
    produces).

    Mechanism: an iterated restricted neighbour-average of the Hc/Hm columns (a Jacobi
    relaxation of the heat equation -- the classic "spread a step out into a ramp" smoother).
    Each node accepts the smoothing in proportion to `accept`:

    - 0 above sea level and grading in over the first `MARGIN_FULL_DEPTH_M` of depth, so land
      keeps all its generated relief and coastlines stay put;
    - 0 outside `MARGIN_TRANSITION_RAD` of opposite-type crust and smoothstep-graded in
      toward the boundary, so oceanic interiors -- including oceanic/oceanic boundaries, where
      a ridge or trench genuinely *is* a sharp step -- keep full abyssal depth.

    Nodes with `accept == 0` never move: they act as the fixed boundary values the relaxation
    ramps between, so a continental shelf stays a shelf until its seaward edge then dives into
    a slope, and oceanic crust near the contact is pulled up into a rise that deepens back to
    the abyssal plain as it leaves the margin.

    Generation-time only. Later tectonics re-sharpen boundaries on their own terms (fresh
    ocean floor erupting at a rift, old floor bending into a trench) -- correctly, since those
    features really are abrupt.
    """
    from . import plates as plates_mod

    points, plates_in_order, counts, is_continental, hc, hm, elevation = _gather_columns(plates)
    if len(points) < 3:
        return

    cont_pts = points[is_continental]
    ocean_pts = points[~is_continental]
    if len(cont_pts) == 0 or len(ocean_pts) == 0:
        return  # a single-crust-type world has no continent/ocean margins to smooth

    n = len(points)
    workers = plates_mod.query_workers(n)
    transition_chord = _chord(MARGIN_TRANSITION_RAD)

    # Distance from every node to the nearest node of the *opposite* crust type.
    dist_to_ocean, _ = cKDTree(ocean_pts).query(points, workers=workers)
    dist_to_cont, _ = cKDTree(cont_pts).query(points, workers=workers)
    dist_to_opposite = np.where(is_continental, dist_to_ocean, dist_to_cont)

    band = (elevation < 0.0) & (dist_to_opposite < transition_chord)
    if not np.any(band):
        return

    # How much each node accepts the relaxation: 0 (a fixed boundary value) on land, in the
    # deep ocean, and outside the transition zone; ramping to 1 with depth and with proximity
    # to the contact.
    depth_t = np.clip((-elevation - COAST_GUARD_DEPTH_M) / (MARGIN_FULL_DEPTH_M - COAST_GUARD_DEPTH_M), 0.0, 1.0)
    zone_t = np.clip(1.0 - dist_to_opposite / transition_chord, 0.0, 1.0)
    zone_w = zone_t * zone_t * (3.0 - 2.0 * zone_t)  # smoothstep -> the ramp has no kink
    accept = np.where(band, depth_t * zone_w, 0.0)

    # Immediate-neighbour graph. Only nodes that either accept smoothing or are a neighbour of
    # one need to take part -- restrict to that active set so the iteration stays cheap even
    # when the whole world is dense.
    full_tree = cKDTree(points)
    k = min(n, _NEIGHBOUR_QUERY_K)
    nbr_dist, nbr_idx = full_tree.query(points, k=k, workers=workers)
    spacing_chord = float(np.median(nbr_dist[:, 1]))
    neighbour = (nbr_dist > 0.0) & (nbr_dist <= NEIGHBOUR_SPACING_MULTIPLE * spacing_chord)

    touches_band = neighbour & accept[nbr_idx].astype(bool)  # a neighbour that accepts smoothing
    active = (accept > 0.0) | np.any(touches_band, axis=1)
    if not np.any(active):
        return

    remap = np.full(n, -1)
    remap[active] = np.arange(int(active.sum()))
    a_nbr_idx = np.where(neighbour[active], remap[nbr_idx[active]], -1)
    a_nbr_valid = a_nbr_idx >= 0
    a_nbr_idx = np.where(a_nbr_valid, a_nbr_idx, 0)

    # A margin node's neighbours are split into same-crust and opposite-crust groups, and the
    # relaxation target weights the two groups *equally* wherever both are present (rather than
    # letting whichever side happens to contribute more lattice nodes dominate). This is what
    # makes the smoothed field actually continuous across the contact: the two disjoint plate
    # lattices only touch along a thin seam with a handful of cross-links, so a plain
    # neighbour average barely couples them and the step survives.
    a_is_cont = is_continental[active]
    a_nbr_same = a_nbr_valid & (is_continental[np.where(a_nbr_valid, nbr_idx[active], 0)] == a_is_cont[:, None])
    a_nbr_opp = a_nbr_valid & ~a_nbr_same
    same_cnt = np.clip(a_nbr_same.sum(axis=1), 1, None)
    opp_cnt = np.clip(a_nbr_opp.sum(axis=1), 1, None)
    has_same = a_nbr_same.any(axis=1)
    has_opp = a_nbr_opp.any(axis=1)
    # A node with no active neighbour at all has nothing to relax toward -- freeze it.
    a_accept = np.where(has_same | has_opp, accept[active], 0.0)

    iters = int(np.clip(round(3.0 * transition_chord / spacing_chord), MARGIN_RELAX_ITERS_MIN, MARGIN_RELAX_ITERS_MAX))

    def _relax(values: np.ndarray) -> np.ndarray:
        v = values[active].copy()
        for _ in range(iters):
            same_mean = np.where(a_nbr_same, v[a_nbr_idx], 0.0).sum(axis=1) / same_cnt
            opp_mean = np.where(a_nbr_opp, v[a_nbr_idx], 0.0).sum(axis=1) / opp_cnt
            target = np.where(
                has_same & has_opp, 0.5 * (same_mean + opp_mean), np.where(has_same, same_mean, opp_mean)
            )
            v += a_accept * (target - v)
        out = values.copy()
        out[active] = v
        return out

    new_hc = np.clip(_relax(hc), lithosphere.MIN_CRUSTAL_THICKNESS_M, None)
    new_hm = np.clip(_relax(hm), lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M, None)

    _write_back(plates_in_order, counts, active, new_hc, new_hm)

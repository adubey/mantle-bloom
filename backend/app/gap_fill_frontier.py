"""Given a set of gap points and the plate(s) genuinely adjacent to them, iteratively grow
those *existing* plates to cover the gap -- one node at a time onto an existing `ElevationLine`
where one sits close enough to extend, or a brand-new single-node line where none does. Wired
at both `lithosphere_plate.py`'s `deform()` (via `LithospherePlate._fill_corner_notch_frontier`)
and `gaps.py`'s whole-sphere sweep (`gaps.fill_gaps_by_growing_neighbours`, which falls back to
spawning a brand-new plate via `gaps._spawn_plate_from_gap` only when nothing is adjacent at
all). See docs/simulation-model.md's "Frontier gap-fill" section for the full writeup.

Every claimed node is a real magma eruption, never a free area grant (see
`_erupt_melted_nodes`'s own docstring). The two ways a node can be claimed:

- **Stretch an existing line.** *Continental claimants only* (a landlocked marginal sea's
  coastline extending over its own recently-drained neighbourhood, `gaps.py`'s own
  `GAP_LAND_ADOPTION_RADIUS_MULT` case): a candidate grid-adjacent to a line this plate already
  had *before this call* draws its material down from that line's own nearest
  `K_STRETCH_SOURCE_NODES` end nodes (same row-claim mass-conservation shape
  `lithosphere_plate._claim_adjacent_territory` already uses for a whole new phi row, applied
  here at single-node granularity -- see `_stretch_extend_line`), then both the new node and
  the thinned source nodes run through the same decompression-melting eruption path
  (`_erupt_melted_nodes`) as everything else. This is the "stretches the plate it's on" case.
- **New magma-filled node.** Every *oceanic* claim (the overwhelming majority -- gaps are, per
  `gaps.py`'s own docstring, almost always open ocean a fully-subducted plate vacated), plus any
  continental candidate with no existing line to extend at all -- including every node added to
  a line *opened during this same call* (that line has no prior column of its own to draw from
  either, no matter how many nodes it's grown to by the time this call ends) -- is thin-seeded
  from scratch (`LithospherePlate._seed_and_erupt_new_nodes`, the same seeding every other
  new-crust event in `deform()` uses) with nothing thinned in exchange.

  GitHub issue #216: an oceanic claimant used to stretch here too (drawing down its own edge
  nodes exactly like the continental case), which is backwards for what this call actually
  represents -- most oceanic gap-fill is a stand-in for mid-ocean-ridge seafloor spreading,
  where the ridge continuously injects *new* magma-derived crust as the flanks separate, not a
  fixed reservoir thinning to cover more area. Mass-conserving stretch left this call's own
  creation on the same self-canceling "thin now, remelt-and-reset later" footing issue #216's
  own investigation found dominating a similarly-shaped in-`deform()` growth path
  (`LithospherePlate._stretch_end` vs. `regularize_line`) -- real seafloor doesn't work that
  way, so an oceanic claim now always erupts fresh instead, regardless of whether it's
  extending a pre-existing line.

Either way, the freshly-erupted node is typed oceanic vs. continental by whether it was above
or below sea level the instant it erupted (`_erupt_melted_nodes`'s own `melt_land`/`melt_ocean`
split) -- oceanic crust underwater, a volcano (bimodal continental-rift volcanism) above it.

Nodes are appended via `ElevationLine.with_new_nodes` -- as of this module, its first real call
site (its own docstring notes none existed before). `with_new_nodes` zero-fills every
`OPTIONAL_FIELDS` entry for the appended nodes and leaves the result theta-unsorted; both
helpers below immediately overwrite the appended tail with the real seeded values (that
docstring's own flagged pitfall: a zero-fill misreads `node_created_years`' -1.0 "unknown"
sentinel as "created at year 0") and re-sort by theta, since every other consumer of a line --
including this module's own end-adjacency checks -- assumes ascending order.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, lithosphere, rheology, terrain_noise
from .elevation_lines import CONTIGUOUS_RUN_GAP_MULT, DEFRAG_CONNECT_RADIUS_MULT, ElevationLine
from .lithosphere_plate import _TERRAIN_SEED_TAG, _erupt_melted_nodes, growth_seed_thickness

if TYPE_CHECKING:
    from .lithosphere_plate import LithospherePlate
    from .world import World

# How many of a stretched line's own existing end nodes share a newly-appended node's mass
# deficit -- same value and same reasoning as lithosphere_plate.K_NEIGHBOUR_ROWS_FOR_MASS_
# CONSERVATION (a small constant, not node-count-scaled: this is about how many nodes share a
# single stretch event, not a distance or a per-call budget).
K_STRETCH_SOURCE_NODES = 2

# One appended node is always exactly one grid step's worth of separation (unlike
# _claim_adjacent_territory's whole-row claim, whose stretch_fraction depends on how much of a
# variable-sized gap is actually behind it) -- so the mass-conservation share below always
# treats it as a full one-step stretch: thin_ratio = 1 / (K_STRETCH_SOURCE_NODES + 1). Both the
# new node's own seed and each of the K source nodes' existing column are scaled to this same
# fraction, mirroring _claim_adjacent_territory's own formula at stretch_fraction == 1.0.
_STRETCH_SHARE_COUNT = K_STRETCH_SOURCE_NODES + 1
_STRETCH_THIN_RATIO = 1.0 / _STRETCH_SHARE_COUNT

# A brand-new (non-stretch) node is seeded well below RIFT_CRITICAL_THICKNESS_M -- same
# formula and same reasoning as _fill_corner_notch_frontier's own `seed_thin_ratio`: comfortably clear
# of the +-amp texture noise _seed_and_erupt_new_nodes adds on top, so every such node is
# *guaranteed* (not just usually) to melt through and erupt via _erupt_melted_nodes, matching
# this module's own "never a silent full-thickness spawn" invariant.
_NEW_MAGMA_SEED_THIN_RATIO = min(1.0, (0.3 * rheology.RIFT_CRITICAL_THICKNESS_M) / growth_seed_thickness()[0])

# Safety-valve floors, same shape as lithosphere_plate.MIN_CORNER_FILL_NODES_PER_STEP /
# CORNER_NOTCH_MIN_WINDOW_ROWS -- the real per-call budget is scaled up from these by the
# caller (candidate count / cluster geometry); these are just the floor for a tiny call.
MIN_FRONTIER_HOPS = 1
MIN_FRONTIER_NODES = 200


def _pole_aligned_row_key(local_phi: float, spacing_rad: float) -> int:
    """Integer row index on the same pole-aligned phi grid `_fill_corner_notch_frontier`'s own
    `row_lo`/`row_hi` snap to (`iter_local_lattice`'s own bound) -- two candidate points this
    close in phi always land on the same row, regardless of which one a plate's real lines
    happen to sit exactly on."""
    max_abs_phi = np.pi / 2 - spacing_rad / 2
    return int(np.round((local_phi + max_abs_phi) / spacing_rad))


def _row_phi_value(row_key: int, spacing_rad: float) -> float:
    max_abs_phi = np.pi / 2 - spacing_rad / 2
    return -max_abs_phi + row_key * spacing_rad


def _sorted_by_theta(line: ElevationLine) -> ElevationLine:
    """`ElevationLine.with_new_nodes` appends unsorted (see its own docstring) -- every
    adjacency check in this module reads `line.theta[0]`/`line.theta[-1]` as the true low/high
    ends, so every append is immediately re-sorted. (Doesn't handle a line spanning the +-pi
    theta seam any more correctly than `_fill_corner_notch_frontier`'s own linear theta_lo/theta_hi math
    already doesn't -- a pre-existing limitation, not one this module introduces.)"""
    order = np.argsort(line.theta)
    return line.masked(order)


class _WorkingPlate:
    """This call's mutable view of one claimant plate: a plain Python copy of its lines (only
    committed back via `plate.set_lines` once, at the very end of `fill_gap_by_growing_plates`
    -- same "accumulate, apply once" shape `_fill_corner_notch_frontier` already uses), plus enough
    per-row bookkeeping to decide "stretch an existing line" vs. "erupt a new one" for each
    node as it's claimed, node by node, in this call's own theta order."""

    def __init__(self, plate: "LithospherePlate", spacing_rad: float) -> None:
        self.plate = plate
        self.lines: list[ElevationLine] = list(plate.lines)
        # Row key -> list of indices into self.lines. Multiple lines can legitimately share a
        # row (interior-subduction carve-outs elsewhere in this codebase) -- adjacency is
        # checked against every line at a row, not just the first.
        self.rows: dict[int, list[int]] = {}
        # Indices present *before* this call started -- only these draw a stretch; anything
        # opened during this call (however many nodes it grows to before this call ends) has
        # no prior column to draw from, so every node on it is a fresh eruption. See module
        # docstring.
        self.pre_existing: set[int] = set(range(len(self.lines)))
        self.own_points, _ = plate.all_points_and_elevation()
        for i, line in enumerate(self.lines):
            if len(line) == 0:
                continue
            row_key = _pole_aligned_row_key(line.phi, spacing_rad)
            self.rows.setdefault(row_key, []).append(i)

    def refresh_own_points(self) -> None:
        pts = [
            geometry.to_world(self.plate.frame, geometry.local_xyz(np.full(len(line), line.phi), line.theta))
            for line in self.lines
            if len(line) > 0
        ]
        self.own_points = np.concatenate(pts, axis=0) if pts else np.zeros((0, 3))


def _seed_new_node(
    world: "World", plate: "LithospherePlate", line_index: int, phi: float, theta: float,
    thin_ratio: float, hc0: float, hm0: float, amp: float, texture: "terrain_noise.FractalTexture",
) -> dict[str, np.ndarray]:
    world_pt = geometry.to_world(plate.frame, geometry.local_xyz(np.array([phi]), np.array([theta])))
    return plate._seed_and_erupt_new_nodes(world, line_index, world_pt, thin_ratio, hc0, hm0, amp, texture)


def _append_seeded_node(line: ElevationLine, new_theta: float, seeded: dict[str, np.ndarray]) -> ElevationLine:
    """`with_new_nodes` plus overwriting its zero-filled tail with `seeded`'s real values (see
    module docstring), re-sorted by theta."""
    extended = line.with_new_nodes(np.array([new_theta]), seeded["elevation"])
    n_old = len(line)
    overrides = {}
    for name, new_vals in seeded.items():
        if name == "elevation":
            continue
        full = getattr(extended, name).copy()
        full[n_old:] = new_vals
        overrides[name] = full
    return _sorted_by_theta(extended.replace(**overrides))


def _erupt_new_single_node(
    world: "World", plate: "LithospherePlate", line: ElevationLine, line_index: int, new_theta: float,
    hc0: float, hm0: float, amp: float, texture: "terrain_noise.FractalTexture",
) -> ElevationLine:
    """Append one brand-new (fresh-eruption, no draw-down) node to `line` -- used both to open
    a genuinely new line (a single-node `line`) and to keep extending one opened earlier this
    same call, which still has no prior column to draw from (see module docstring)."""
    seeded = _seed_new_node(world, plate, line_index, line.phi, new_theta, _NEW_MAGMA_SEED_THIN_RATIO, hc0, hm0, amp, texture)
    return _append_seeded_node(line, new_theta, seeded)


def _stretch_extend_line(
    world: "World", plate: "LithospherePlate", line: ElevationLine, line_index: int,
    at_low_end: bool, new_theta: float, hc0: float, hm0: float, amp: float, texture: "terrain_noise.FractalTexture",
) -> ElevationLine:
    """One new node, drawing its material down from `line`'s own nearest
    `K_STRETCH_SOURCE_NODES` end nodes (the end `new_theta` is extending) at
    `_STRETCH_THIN_RATIO` -- see this module's own docstring for the mass-conservation shape
    this mirrors. Returns the replacement line (new node appended, source nodes thinned)."""
    seeded = _seed_new_node(world, plate, line_index, line.phi, new_theta, _STRETCH_THIN_RATIO, hc0, hm0, amp, texture)
    extended = _append_seeded_node(line, new_theta, seeded)

    n_old = len(line)
    k = min(K_STRETCH_SOURCE_NODES, n_old)
    if k == 0:
        return extended
    # extended is theta-sorted; the new node's own index within it is wherever new_theta
    # landed, not necessarily n_old -- find the K nearest *other* nodes on the correct side
    # instead of assuming positions, since re-sorting can have moved everything.
    new_idx = int(np.searchsorted(extended.theta, new_theta))
    if at_low_end:
        source_idx = np.arange(new_idx + 1, min(new_idx + 1 + k, len(extended)))
    else:
        source_idx = np.arange(max(new_idx - k, 0), new_idx)
    if len(source_idx) == 0:
        return extended

    hc = extended.crustal_thickness_m.copy()
    hm = extended.mantle_lithosphere_thickness_m.copy()
    crust_type_code = extended.crust_type_code.copy()
    is_volcano = extended.is_volcano.copy()
    volcano_remaining = extended.volcano_active_years_remaining.copy()
    prior_elevation = extended.elevation.copy()

    old_hc = hc[source_idx].copy()
    new_hc = old_hc * _STRETCH_THIN_RATIO
    new_hm = hm[source_idx] * _STRETCH_THIN_RATIO
    melting = (old_hc >= rheology.RIFT_CRITICAL_THICKNESS_M) & (new_hc < rheology.RIFT_CRITICAL_THICKNESS_M)
    hc[source_idx] = new_hc
    hm[source_idx] = new_hm
    sub_crust_type = crust_type_code[source_idx]
    sub_is_volcano = is_volcano[source_idx]
    sub_volcano_remaining = volcano_remaining[source_idx]
    _erupt_melted_nodes(
        world, plate.plate_id, line_index,
        hc[source_idx], hm[source_idx], sub_crust_type, sub_is_volcano, sub_volcano_remaining,
        melting, prior_elevation[source_idx],
    )
    crust_type_code[source_idx] = sub_crust_type
    is_volcano[source_idx] = sub_is_volcano
    volcano_remaining[source_idx] = sub_volcano_remaining
    node_rho_c = lithosphere.node_crust_density(crust_type_code[source_idx], plate.crust_type)
    new_elevation = prior_elevation.copy()
    new_elevation[source_idx] = lithosphere.isostatic_elevation(hc[source_idx], hm[source_idx], node_rho_c)
    return extended.replace(
        crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm, crust_type_code=crust_type_code,
        is_volcano=is_volcano, volcano_active_years_remaining=volcano_remaining, elevation=new_elevation,
    )


def _claim_row_points(
    world: "World", wp: _WorkingPlate, row_key: int, phi: float, gap_tol: float,
    thetas_sorted: np.ndarray, hc0: float, hm0: float, amp: float, texture: "terrain_noise.FractalTexture",
    budget: int,
) -> tuple[int, int]:
    """Claim as many of `thetas_sorted` (ascending, raw local theta -- not snapped to any fixed
    grid, since a real plate's own lines don't reliably sit on one after regularize_line's own
    resampling) as fit within `budget` real new nodes. Deliberately does *not* require exact
    grid alignment: a candidate within `gap_tol` (`CONTIGUOUS_RUN_GAP_MULT * dtheta`, the same
    threshold `split_into_contiguous_runs` already uses to call two nodes "still one contiguous
    row") of a line's own low/high end extends that line (stretch if it pre-dates this call
    *and* this claimant is continental, fresh eruption otherwise -- see module docstring); a
    candidate that already falls *inside*
    an existing line's own theta span (within the same tolerance) is a duplicate/near-duplicate
    detection of an already-covered point and is silently absorbed, no new node. Anything with
    no nearby line opens a new one. `needs_regularizing`/`regularize_line` (run every `deform()`
    pass regardless of which gap-fill algorithm is selected) re-densifies the result to target
    spacing afterward -- the same tolerance-then-resample shape `_stretch_end` already relies
    on, see this module's own docstring.

    Mutates `wp` in place. Returns `(nodes_added, consumed)`: `nodes_added` counts only genuine
    new nodes (what the caller charges against its own budget/stats); `consumed` additionally
    counts absorbed duplicates, so the caller knows how many *original* candidates (in the same
    ascending order) this call actually resolved one way or the other -- anything beyond that
    stays unresolved for a later hop/call."""
    added = 0
    consumed = 0
    for theta in thetas_sorted:
        if added >= budget:
            break
        theta = float(theta)
        best: tuple[int, bool | None, float] | None = None
        for idx in wp.rows.get(row_key, []):
            line = wp.lines[idx]
            if len(line) == 0:
                continue
            lo, hi = float(line.theta[0]), float(line.theta[-1])
            if lo - gap_tol <= theta <= hi + gap_tol:
                if lo <= theta <= hi:
                    best = (idx, None, 0.0)
                    break
                dist = (lo - theta) if theta < lo else (theta - hi)
                if best is None or dist < best[2]:
                    best = (idx, theta < lo, dist)
        if best is not None and best[1] is None:
            consumed += 1  # duplicate/near-duplicate of an already-covered point, no-op
            continue
        if best is not None:
            idx, at_low_end, _ = best
            line = wp.lines[idx]
            # GitHub issue #216: only a *continental* claimant stretches (draws material down
            # from its own nearest end nodes) -- see this module's own docstring for why an
            # oceanic claimant always erupts fresh crust instead, regardless of whether it's
            # extending a pre-existing line.
            stretch = (idx in wp.pre_existing) and wp.plate.crust_type != "oceanic"
            wp.lines[idx] = (
                _stretch_extend_line(world, wp.plate, line, idx, at_low_end, theta, hc0, hm0, amp, texture)
                if stretch
                else _erupt_new_single_node(world, wp.plate, line, idx, theta, hc0, hm0, amp, texture)
            )
        else:
            seeded = _seed_new_node(world, wp.plate, len(wp.lines), phi, theta, _NEW_MAGMA_SEED_THIN_RATIO, hc0, hm0, amp, texture)
            new_line = ElevationLine(phi=phi, theta=np.array([theta]), **seeded)
            wp.lines.append(new_line)
            wp.rows.setdefault(row_key, []).append(len(wp.lines) - 1)
            # Deliberately never added to wp.pre_existing -- see module docstring.
        added += 1
        consumed += 1
    return added, consumed


def fill_gap_by_growing_plates(
    world: "World", gap_points: np.ndarray, claimants: list["LithospherePlate"], spacing_rad: float,
    max_hops: int | None = None, max_nodes: int | None = None,
) -> dict[int, int]:
    """Grow `claimants` to cover `gap_points` (world xyz), node by node, walking the connected
    frontier outward in `DEFRAG_CONNECT_RADIUS_MULT`-sized hops exactly like
    `_fill_corner_notch_frontier` does (so a claim never reads as a disconnected stray to
    `merge_split.defragment_plates`'s own connected-components check). Each hop: every
    still-uncovered gap point within one connect-radius of *some* claimant's current node cloud
    is assigned to whichever claimant is nearest (this is "detect adjacent plates" applied per
    point when there's more than one candidate -- a single-claimant call, the corner-notch
    site, trivially assigns everything to that one plate); assigned points are grouped by
    claimant and by pole-aligned phi row, then claimed in ascending raw theta order via
    `_claim_row_points` (see its own docstring for why no synthetic grid-snapping happens).

    Returns `{plate_id: nodes_added}` for whichever claimants actually grew (omits any that
    didn't). Never spawns a new plate itself -- an empty `claimants` list is a no-op; the
    caller (`gaps.fill_gaps_by_growing_neighbours`) falls back to spawning one of its own when
    there's truly nothing nearby to grow."""
    if len(claimants) == 0 or len(gap_points) == 0:
        return {}
    connect_radius_rad = DEFRAG_CONNECT_RADIUS_MULT * spacing_rad
    if max_hops is None:
        _centroid, radius_rad = geometry.bounding_sphere(gap_points)
        max_hops = max(MIN_FRONTIER_HOPS, int(np.ceil(2.0 * radius_rad / connect_radius_rad)))
    if max_nodes is None:
        max_nodes = max(MIN_FRONTIER_NODES, len(gap_points))

    working = [_WorkingPlate(p, spacing_rad) for p in claimants]
    hc0, hm0 = growth_seed_thickness()
    amp = hc0 * 0.1
    texture_by_plate = {
        p.plate_id: terrain_noise.FractalTexture(np.random.default_rng((world.seed, p.plate_id, _TERRAIN_SEED_TAG)))
        for p in claimants
    }

    remaining = np.ones(len(gap_points), dtype=bool)
    nodes_added: dict[int, int] = {}
    total_added = 0
    for _hop in range(max_hops):
        if not np.any(remaining) or total_added >= max_nodes:
            break
        own_points = [wp.own_points for wp in working]
        counts = [len(pts) for pts in own_points]
        if sum(counts) == 0:
            break
        tree = cKDTree(np.concatenate(own_points, axis=0))
        owner_of = np.concatenate([np.full(c, i) for i, c in enumerate(counts)])

        idx_remaining = np.flatnonzero(remaining)
        dist, nearest = tree.query(gap_points[idx_remaining])
        reachable = dist <= connect_radius_rad
        if not np.any(reachable):
            break
        owners = owner_of[nearest]

        # Aligned 1:1 with idx_remaining -- which of this hop's candidates actually got a node
        # created for them (vs. left for a later hop/call by a budget cutoff).
        covered_this_hop = np.zeros(len(idx_remaining), dtype=bool)

        for wp_index, wp in enumerate(working):
            if total_added >= max_nodes:
                break
            mine = reachable & (owners == wp_index)
            if not np.any(mine):
                continue
            mine_positions = np.flatnonzero(mine)
            local = geometry.to_local(wp.plate.frame, gap_points[idx_remaining[mine_positions]])
            local_phi, local_theta = geometry.xyz_to_latlon(local)
            row_keys = np.array([_pole_aligned_row_key(p, spacing_rad) for p in local_phi])

            for row_key in np.unique(row_keys):
                if total_added >= max_nodes:
                    break
                row_key = int(row_key)
                in_row = np.flatnonzero(row_keys == row_key)
                existing = wp.rows.get(row_key)
                ref_phi = float(wp.lines[existing[0]].phi) if existing else _row_phi_value(row_key, spacing_rad)
                dtheta = spacing_rad / max(np.cos(ref_phi), 1e-3)
                gap_tol = CONTIGUOUS_RUN_GAP_MULT * dtheta

                order = np.argsort(local_theta[in_row])
                sorted_in_row = in_row[order]  # positions into mine_positions, ascending theta
                thetas_sorted = local_theta[sorted_in_row]

                budget = max_nodes - total_added
                if budget <= 0:
                    break
                added, consumed = _claim_row_points(
                    world, wp, row_key, ref_phi, gap_tol, thetas_sorted, hc0, hm0, amp,
                    texture_by_plate[wp.plate.plate_id], budget,
                )
                total_added += added
                if added:
                    nodes_added[wp.plate.plate_id] = nodes_added.get(wp.plate.plate_id, 0) + added

                # `_claim_row_points` resolves thetas_sorted strictly in order -- the first
                # `consumed` of them (new nodes plus absorbed duplicates) are covered; anything
                # past that stays `remaining` for a later hop/call.
                covered_this_hop[mine_positions[sorted_in_row[:consumed]]] = True

        if not np.any(covered_this_hop):
            break
        remaining[idx_remaining[covered_this_hop]] = False
        for wp in working:
            wp.refresh_own_points()

    for wp in working:
        wp.plate.set_lines(wp.lines)
    return nodes_added

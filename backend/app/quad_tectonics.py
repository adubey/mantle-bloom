"""Per-step tectonics for `PlateWithSparseQuadPatch` plates (issue #228 Phase 4).

The line-backed engine (`lithosphere_plate.LithospherePlate.deform`) moves boundaries by
editing the two ends of plate-local latitude rows: end-trim for retreat, end-stretch plus
whole-row claims plus a diagonal corner-notch filler for advance, and a whole-leading-row drop
for the "parallel suture" case end-trim can't reach. Every one of those exists because a row
can only grow or shrink along its own axis. A quad surface has no preferred axis, so this
module replaces all of them with two operations on the cell graph:

- **Retreat** peels the plate's contested boundary in layers: any retreatable cell with a
  wholly exposed side (the 2D analogue of a line end) is deactivated, then the newly exposed
  layer is considered, up to this step's displacement in cells. A continental suture's
  consumed crust is thrust onto the nearest surviving cells (area-weighted, so volume is
  conserved up to the usual accretion cap), exactly as `_redistribute_accreted_column` does
  for a line end.
- **Advance** activates the empty cell across each exposed side of an eligible boundary
  cell, wherever no neighbour already covers it, again in layers. New cells are fresh crust
  through the same `seed_and_erupt_new_nodes` path the line engine's row claims use, and the
  cells behind them are stretch-thinned by the same share (`K_NEIGHBOUR_ROWS_FOR_MASS_
  CONSERVATION` layers instead of rows). An active continental margin grows juvenile arc crust
  instead (`ARC_MARGIN_SEED_*`).

Boundary classification and the per-node column physics (convergent/arc/divergent/transform
updates, melting, provenance) are shared with the line engine unchanged -- see
`lithosphere_plate.boundary_context` and `deform_columns`. Nothing here needs regularizing:
cells never drift off the lattice.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import geometry, lithosphere, phase_budget, rheology, terrain_noise
from .elevation_lines import COVERAGE_RADIUS_MULT, ELEV_CHANGE_SUBDUCTION_ARC, line_spacing_rad
from .lithosphere_plate import (
    ARC_MARGIN_END_SCAN_NODES,
    ARC_MARGIN_SEED_HC_M,
    ARC_MARGIN_SEED_HM_M,
    COLUMN_FIELDS,
    CONTINENTAL_CONTESTED_RETREAT_MIN_RUN,
    EXTEND_THRESHOLD_MULTIPLIER,
    K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION,
    MAX_CLAIM_ROWS_PER_STEP,
    MAX_EXTEND_NODES_PER_STEP,
    SUTURE_ACCRETION_MAX_HC_M,
    SUTURE_ACCRETION_SPREAD_NODES,
    _TERRAIN_SEED_TAG,
    _erupt_melted_nodes,
    boundary_context,
    deform_columns,
    growth_seed_thickness,
    seed_and_erupt_new_nodes,
)
from .plates import _contested_by_any
from .surface_fields import SURFACE_FIELDS

if TYPE_CHECKING:
    from .sparse_quad_patch import PlateWithSparseQuadPatch
    from .world import World

# How many layers of cells a boundary may advance in one step -- the same per-step ceiling the
# line engine's row claim uses (MAX_CLAIM_ROWS_PER_STEP), now applying in every direction.
MAX_ADVANCE_LAYERS_PER_STEP = MAX_CLAIM_ROWS_PER_STEP

# The eruption rng stream keys `deform_columns` / `seed_and_erupt_new_nodes` use in place of a
# line index. Distinct values keep the column pass's melting draw and each advance layer's
# new-crust and thinning draws from sharing one (seed, step, plate) stream: the advance uses
# `_ADVANCE_RNG_INDEX + 2 * layer` for new cells and the next value for the cells behind them.
_COLUMN_RNG_INDEX = 0
_ADVANCE_RNG_INDEX = 1


def _adjacency_matrix(plate: "PlateWithSparseQuadPatch") -> csr_matrix:
    graph = plate.adjacency()
    n = plate.node_count()
    rows = np.repeat(np.arange(n), np.diff(graph.offsets))
    return csr_matrix((np.ones(len(graph.neighbours), dtype=np.int8), (rows, graph.neighbours)), shape=(n, n))


def hop_distance(plate: "PlateWithSparseQuadPatch", mask: np.ndarray, width: int) -> np.ndarray:
    """Per cell, edge hops to the nearest `mask` cell, saturating at `width + 1` -- the cell
    graph's version of `lithosphere_plate._distance_to_mask_1d`."""
    dist = np.where(mask, 0, width + 1)
    if width <= 0 or not np.any(mask):
        return dist
    adjacency = _adjacency_matrix(plate)
    reached = np.asarray(mask, dtype=bool).copy()
    frontier = reached.copy()
    for hop in range(1, width + 1):
        frontier = (adjacency @ frontier.astype(np.int8) > 0) & ~reached
        if not np.any(frontier):
            break
        dist[frontier] = hop
        reached |= frontier
    return dist


def components_of_at_least(plate: "PlateWithSparseQuadPatch", mask: np.ndarray, min_size: int) -> np.ndarray:
    """`mask`, with every edge-connected component smaller than `min_size` cleared -- the cell
    graph's version of `lithosphere_plate._runs_of_at_least`."""
    mask = np.asarray(mask, dtype=bool)
    if min_size <= 1 or not np.any(mask):
        return mask.copy()
    members = np.flatnonzero(mask)
    sub = _adjacency_matrix(plate)[members][:, members]
    _, labels = connected_components(sub, directed=False)
    sizes = np.bincount(labels)
    out = np.zeros_like(mask)
    out[members[sizes[labels] >= min_size]] = True
    return out


def deform(plate: "PlateWithSparseQuadPatch", world: "World", other_plates: list, years: float, max_distance: float) -> None:
    """The quad counterpart of `LithospherePlate.deform` -- see this module's docstring."""
    if plate.node_count() == 0:
        return
    spacing_rad = line_spacing_rad(world.node_density)
    nominal_area_m2 = lithosphere.node_area_m2(spacing_rad)
    areas = plate.node_areas_m2()
    ctx = boundary_context(
        world,
        plate,
        other_plates,
        years,
        lambda contested: components_of_at_least(plate, contested, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN),
        node_weight=areas / nominal_area_m2,
    )
    near_field_dist = hop_distance(plate, ctx.convergent, ctx.orogen_dilation_nodes) if ctx.orogen_dilation_nodes > 0 else None
    columns = deform_columns(
        world,
        plate,
        ctx,
        slice(None),
        {name: plate.collect(name) for name in COLUMN_FIELDS},
        near_field_dist,
        lambda: plate.surface_nodes().local_xyz,
        areas,
        _COLUMN_RNG_INDEX,
        years,
    )
    plate.set_fields_on_plate(**columns)

    max_cells = max(1, round(MAX_EXTEND_NODES_PER_STEP * np.sqrt(world.node_density)))
    before = phase_budget.snapshot(plate) if world.debug_diagnostics else None
    survivors = _retreat(plate, world, ctx, max_distance, max_cells)
    if world.debug_diagnostics:
        after = phase_budget.snapshot(plate)
        phase_budget.record(world, plate, "boundary_retreat", *before, *after)
        before = after

    if not ctx.suppress_growth:
        # Wide enough to see every plate a MAX_ADVANCE_LAYERS_PER_STEP-deep advance could run
        # into, not just the ones inside the boundary-force reach.
        growth_neighbours = plate.get_neighbours(
            other_plates, threshold_rad=ctx.reach_rad + (MAX_ADVANCE_LAYERS_PER_STEP + 1) * spacing_rad
        )
        _advance(plate, world, ctx, survivors, growth_neighbours, spacing_rad, max_cells)
        if world.debug_diagnostics:
            phase_budget.record(world, plate, "boundary_advance", *before, *phase_budget.snapshot(plate))


def _retreat(plate: "PlateWithSparseQuadPatch", world: "World", ctx, max_distance: float, max_cells: int) -> np.ndarray:
    """Peel retreatable boundary cells, layer by layer -- see the module docstring. Returns the
    survivor mask over the pre-retreat node order."""
    n = plate.node_count()
    survivors = np.ones(n, dtype=bool)
    retreatable = ctx.shrinkable.copy()
    if not np.any(retreatable) or n <= 1:
        return survivors

    probes = plate._probe_neighbour_indices()
    has_probe = probes >= 0
    safe_probes = np.where(has_probe, probes, 0)
    weights = plate.node_areas_m2() / lithosphere.node_area_m2(ctx.spacing_rad)
    hc = plate.collect("crustal_thickness_m")
    budget_hc = ctx.oceanic_override_retreat_budget_hc
    remaining = min(max_cells, n - 1)
    removed = np.zeros(n, dtype=bool)
    for _ in range(max(1, int(max_distance / ctx.spacing_rad))):
        # Exposed now: a whole side with no leaf at all, or only leaves already peeled.
        open_half = ~has_probe | removed[safe_probes]
        exposed = np.flatnonzero(retreatable & ~removed & np.any(np.all(open_half, axis=2), axis=1))
        take = np.zeros(n, dtype=bool)
        take[exposed[:remaining]] = True
        if plate.crust_type == "continental":
            # Issue #177: a continental passive margin overridden by an oceanic neighbour may
            # only subduct as much volume as this step's arc magmatism creates -- the same
            # budget the line engine spends end by end (see `_budget_limited_removal`).
            override = np.flatnonzero(take & ~ctx.accrete)
            cost = np.cumsum(hc[override] * weights[override])
            allowed = cost <= max(float(budget_hc[0]), 0.0)
            take[override[~allowed]] = False
            if np.any(allowed):
                budget_hc[0] -= float(cost[allowed][-1])
        chosen = np.flatnonzero(take)
        if not len(chosen):
            break
        removed[chosen] = True
        remaining -= len(chosen)
        if remaining <= 0:
            break
    if not np.any(removed):
        return survivors

    own_points = ctx.own_points
    world.record_removed_points(own_points[removed], plate.plate_id)
    donors = removed & ctx.accrete
    if np.any(donors):
        _accrete_onto_survivors(plate, donors, ~removed)
    plate.remove_cells(removed)
    return ~removed


def _accrete_onto_survivors(plate: "PlateWithSparseQuadPatch", donors: np.ndarray, survivors: np.ndarray) -> None:
    """Thrust each continental suture's consumed Hc/Hm volume onto the surviving cells within
    `SUTURE_ACCRETION_SPREAD_NODES` hops behind it, as a uniform thickening -- the 2D form of
    `_redistribute_accreted_column`, which spreads a retreating line end's volume over that
    many nodes *inward along the row*. Each edge-connected run of donor cells is one suture
    front with its own band, so separate sutures on one plate don't share volume. Conserves
    volume by exact cell area up to the same caps; a front with no survivor in reach feeds its
    nearest survivor instead."""
    survivor_idx = np.flatnonzero(survivors)
    if not len(survivor_idx):
        return
    areas = plate.node_areas_m2()
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    elevation = plate.collect("elevation")
    add_hc = np.zeros(len(hc))
    add_hm = np.zeros(len(hm))

    donor_idx = np.flatnonzero(donors)
    _, labels = connected_components(_adjacency_matrix(plate)[donor_idx][:, donor_idx], directed=False)
    for label in np.unique(labels):
        front = donor_idx[labels == label]
        mask = np.zeros(len(hc), dtype=bool)
        mask[front] = True
        band = survivors & (hop_distance(plate, mask, SUTURE_ACCRETION_SPREAD_NODES) <= SUTURE_ACCRETION_SPREAD_NODES)
        if not np.any(band):
            points = plate.surface_nodes().local_xyz
            _, nearest = cKDTree(points[survivor_idx]).query(points[front].mean(axis=0))
            band[survivor_idx[nearest]] = True
        band_area = float(areas[band].sum())
        add_hc[band] += float(np.sum(hc[front] * areas[front])) / band_area
        add_hm[band] += float(np.sum(hm[front] * areas[front])) / band_area

    gained = np.flatnonzero(add_hc > 0.0)
    density = lithosphere.node_crust_density(plate.collect("crust_type_code")[gained], plate.crust_type)
    before = lithosphere.isostatic_elevation(hc[gained], hm[gained], density)
    hc[gained] = np.minimum(hc[gained] + add_hc[gained], SUTURE_ACCRETION_MAX_HC_M)
    hm[gained] = np.minimum(hm[gained] + add_hm[gained], lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)
    after = lithosphere.isostatic_elevation(hc[gained], hm[gained], density)
    elevation[gained] = rheology.clip_elevation_bounds(elevation[gained] + (after - before))
    plate.set_fields_on_plate(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm, elevation=elevation)


def _advance(
    plate: "PlateWithSparseQuadPatch",
    world: "World",
    ctx,
    survivors: np.ndarray,
    neighbours: list,
    spacing_rad: float,
    max_cells: int,
) -> None:
    """Grow the boundary into open ground, layer by layer -- see the module docstring."""
    inputs = ctx.inputs
    extend_threshold_rad = EXTEND_THRESHOLD_MULTIPLIER * spacing_rad
    # Sources for the first layer: uncontested boundary cells with open water between them
    # and the nearest neighbour -- the line engine's end-growth condition, per cell.
    eligible = (~ctx.contested & (inputs.dist_to_neighbor > extend_threshold_rad))[survivors]
    if not np.any(eligible):
        return
    arc_source = np.zeros(plate.node_count(), dtype=bool)
    if plate.crust_type == "continental":
        # Active margin: an arc-band cell, or one still stamped by a recent subduction-arc
        # step, within ARC_MARGIN_END_SCAN_NODES hops -- the per-cell version of the line
        # engine's end scan.
        signal = ctx.arc_band[survivors] | (plate.collect("elev_change_reason") == ELEV_CHANGE_SUBDUCTION_ARC)
        arc_source = hop_distance(plate, signal, ARC_MARGIN_END_SCAN_NODES - 1) < ARC_MARGIN_END_SCAN_NODES
    grow_frontier(plate, world, eligible, arc_source, neighbours, spacing_rad, MAX_ADVANCE_LAYERS_PER_STEP, max_cells)


def grow_frontier(
    plate: "PlateWithSparseQuadPatch",
    world: "World",
    eligible: np.ndarray,
    arc_source: np.ndarray,
    neighbours: list,
    spacing_rad: float,
    max_layers: int,
    max_cells: int,
    claimable=None,
) -> int:
    """Activate open empty cells across the exposed sides of `eligible` cells, then across
    the newest layer's, for up to `max_layers` layers and `max_cells` cells -- the shared
    areal-growth walk behind boundary advance and gap filling. A candidate is open when no
    plate in `neighbours` contains it and none of their nodes lies within
    `EXTEND_THRESHOLD_MULTIPLIER` spacings; `claimable(world_pts)`, when given, narrows that
    further. Cells grown from an `arc_source` cell are arc crust; the rest are fresh crust
    that stretch-thins the cells behind them. Returns how many cells were added."""
    extend_threshold_rad = EXTEND_THRESHOLD_MULTIPLIER * spacing_rad
    neighbour_points = [p.all_points_and_elevation()[0] for p in neighbours if p.node_count() > 0]
    neighbour_tree = cKDTree(np.concatenate(neighbour_points, axis=0)) if neighbour_points else None
    hc0, hm0 = growth_seed_thickness()
    amp = hc0 * 0.1
    texture = terrain_noise.FractalTexture(np.random.default_rng((world.seed, plate.plate_id, _TERRAIN_SEED_TAG)))

    # Per current node: may it grow this layer, and is it an arc margin. Kept aligned with
    # the plate's node order across each insertion.
    keys = plate.cell_keys.copy()
    state = {"eligible": np.asarray(eligible, dtype=bool), "arc": np.asarray(arc_source, dtype=bool)}
    budget = max_cells
    for layer in range(max_layers):
        if budget <= 0:
            break
        sources, candidates = plate.empty_neighbour_keys(np.flatnonzero(state["eligible"]))
        if not len(candidates):
            break
        local = plate.cell_centres_local(candidates)
        world_pts = geometry.to_world(plate.frame, local)
        open_mask = ~_contested_by_any(world_pts, neighbours)
        if claimable is not None:
            open_mask &= claimable(world_pts)
        gap_direction = np.zeros((len(candidates), 3))
        if neighbour_tree is not None:
            dist, idx = neighbour_tree.query(world_pts)
            open_mask &= dist > extend_threshold_rad
            finite = np.isfinite(dist)
            gap_direction[finite] = neighbour_tree.data[idx[finite]] - world_pts[finite]
        sources, candidates, local, world_pts, gap_direction = (
            array[open_mask][:budget] for array in (sources, candidates, local, world_pts, gap_direction)
        )
        if not len(candidates):
            break
        arc = state["arc"][sources]
        seed_rng, thin_rng = _ADVANCE_RNG_INDEX + 2 * layer, _ADVANCE_RNG_INDEX + 2 * layer + 1
        fields, thin_ratio = _new_cell_fields(plate, world, seed_rng, local, world_pts, gap_direction, sources, arc, hc0, hm0, amp, texture)
        inserted = plate.insert_cells(candidates, fields)
        if not np.any(inserted):
            break
        budget -= int(inserted.sum())

        # Re-align per-node state with the new node order: old nodes by key, new nodes
        # inherit from the source that grew them.
        new_keys = plate.cell_keys
        old_pos = np.searchsorted(keys, new_keys)
        is_old = (old_pos < len(keys)) & (keys[np.minimum(old_pos, len(keys) - 1)] == new_keys)
        new_index = plate._index_of_keys(candidates[inserted])
        grown_from = sources[inserted]
        for name, values in state.items():
            aligned = np.zeros((len(new_keys),) + values.shape[1:], dtype=values.dtype)
            aligned[is_old] = values[old_pos[is_old]]
            aligned[new_index] = values[grown_from]
            state[name] = aligned
        # Only the newest layer grows next.
        state["eligible"] = np.zeros(len(new_keys), dtype=bool)
        state["eligible"][new_index] = True
        keys = new_keys.copy()

        ordinary = ~arc[inserted]
        if np.any(ordinary):
            _thin_behind(plate, world, thin_rng, new_index[ordinary], thin_ratio[inserted][ordinary])
    return max_cells - budget


def _new_cell_fields(
    plate: "PlateWithSparseQuadPatch",
    world: "World",
    rng_index: int,
    local: np.ndarray,
    world_pts: np.ndarray,
    gap_direction: np.ndarray,
    sources: np.ndarray,
    arc: np.ndarray,
    hc0: float,
    hm0: float,
    amp: float,
    texture: "terrain_noise.FractalTexture",
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Field values for each candidate cell, and the stretch ratio its growth imposes on the
    cells behind it (1.0 for arc crust, which is fed by the slab rather than by stretching).

    The ratio is the line engine's row-claim share, `1 - f (S - 1) / S` with
    `S = K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION + 1`: `f` is how much of this cell's growth
    runs along the separation direction (the component of the outward step toward the nearest
    neighbour -- 1 with no neighbour in view, where the whole step is open ground)."""
    count = len(local)
    source_local = plate.surface_nodes().local_xyz[sources]
    outward = geometry.normalize(local - source_local)
    toward = np.zeros_like(outward)
    has_gap = np.linalg.norm(gap_direction, axis=1) > 0.0
    if np.any(has_gap):
        toward[has_gap] = geometry.normalize(geometry.to_local(plate.frame, gap_direction[has_gap]))
    alignment = np.where(has_gap, np.abs(np.sum(outward * toward, axis=1)), 1.0)
    share = K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION + 1
    thin_ratio = np.where(arc, 1.0, 1.0 - np.clip(alignment, 0.0, 1.0) * (share - 1) / share)

    fields = {name: np.full(count, SURFACE_FIELDS[name].default, dtype=SURFACE_FIELDS[name].dtype) for name in (
        "elevation", "crustal_thickness_m", "mantle_lithosphere_thickness_m", "crust_type_code", "is_volcano",
        "volcano_active_years_remaining", "elev_change_reason", "node_created_years",
    )}
    ordinary = ~arc
    if np.any(ordinary):
        seeded = seed_and_erupt_new_nodes(world, plate, rng_index, world_pts[ordinary], thin_ratio[ordinary], hc0, hm0, amp, texture)
        for name, values in seeded.items():
            fields[name][ordinary] = values
    if np.any(arc):
        density = lithosphere.crust_density(plate.crust_type)
        fields["crustal_thickness_m"][arc] = ARC_MARGIN_SEED_HC_M
        fields["mantle_lithosphere_thickness_m"][arc] = ARC_MARGIN_SEED_HM_M
        fields["elevation"][arc] = lithosphere.isostatic_elevation(np.array([ARC_MARGIN_SEED_HC_M]), np.array([ARC_MARGIN_SEED_HM_M]), density)[0]
        fields["elev_change_reason"][arc] = ELEV_CHANGE_SUBDUCTION_ARC
        fields["node_created_years"][arc] = world.elapsed_years
    return fields, thin_ratio


def _thin_behind(plate: "PlateWithSparseQuadPatch", world: "World", rng_index: int, grown: np.ndarray, ratio: np.ndarray) -> None:
    """Stretch-thin the `K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION` layers of older cells behind
    freshly grown ones by each new cell's share -- the per-cell form of the line engine's
    row-claim draw-down, including erupting any column that thins past the rift threshold.
    A cell behind several new ones takes the strongest (smallest) ratio."""
    n = plate.node_count()
    adjacency = _adjacency_matrix(plate)
    grown_mask = np.zeros(n, dtype=bool)
    grown_mask[grown] = True
    node_ratio = np.ones(n)
    node_ratio[grown] = ratio
    applied = np.ones(n)
    frontier = grown_mask.copy()
    reached = grown_mask.copy()
    for _ in range(K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION):
        rows, cols = adjacency[np.flatnonzero(frontier)].nonzero()
        sources = np.flatnonzero(frontier)[rows]
        step = np.ones(n)
        np.minimum.at(step, cols, node_ratio[sources])
        nxt = np.zeros(n, dtype=bool)
        nxt[cols] = True
        nxt &= ~reached
        if not np.any(nxt):
            break
        applied[nxt] = step[nxt]
        node_ratio[nxt] = step[nxt]
        reached |= nxt
        frontier = nxt
    thinned = applied < 1.0
    if not np.any(thinned):
        return
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    codes = plate.collect("crust_type_code")
    is_volcano = plate.collect("is_volcano")
    remaining = plate.collect("volcano_active_years_remaining")
    elevation = plate.collect("elevation")
    new_hc = hc[thinned] * applied[thinned]
    new_hm = hm[thinned] * applied[thinned]
    melting = (hc[thinned] >= rheology.RIFT_CRITICAL_THICKNESS_M) & (new_hc < rheology.RIFT_CRITICAL_THICKNESS_M)
    sub_codes, sub_volcano, sub_remaining = codes[thinned], is_volcano[thinned], remaining[thinned]
    _erupt_melted_nodes(world, plate.plate_id, rng_index, new_hc, new_hm, sub_codes, sub_volcano, sub_remaining, melting, elevation[thinned])
    hc[thinned], hm[thinned] = new_hc, new_hm
    codes[thinned], is_volcano[thinned], remaining[thinned] = sub_codes, sub_volcano, sub_remaining
    elevation[thinned] = lithosphere.isostatic_elevation(new_hc, new_hm, lithosphere.node_crust_density(sub_codes, plate.crust_type))
    plate.set_fields_on_plate(
        crustal_thickness_m=hc,
        mantle_lithosphere_thickness_m=hm,
        crust_type_code=codes,
        is_volcano=is_volcano,
        volcano_active_years_remaining=remaining,
        elevation=elevation,
    )


def fill_gap(
    world: "World", plate: "PlateWithSparseQuadPatch", gap_points: np.ndarray, others: list, spacing_rad: float, max_layers: int
) -> int:
    """Grow `plate` into the uncovered `gap_points` (world xyz) -- the quad counterpart of
    `gap_fill_frontier.fill_gap_by_growing_plates` for one claimant. The same frontier walk as
    boundary advance, from every boundary cell, restricted to cells whose centre lies within
    `COVERAGE_RADIUS_MULT` spacings of a gap point and kept clear of `others`. Returns how many
    cells were added."""
    if plate.node_count() == 0 or len(gap_points) == 0:
        return 0
    gap_tree = cKDTree(gap_points)
    coverage_radius_rad = COVERAGE_RADIUS_MULT * spacing_rad

    def near_gap(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = gap_tree.query(world_pts)
        return dist <= coverage_radius_rad

    neighbours = plate.get_neighbours(others, threshold_rad=(max_layers + 1) * spacing_rad)
    n = plate.node_count()
    return grow_frontier(
        plate, world, np.ones(n, dtype=bool), np.zeros(n, dtype=bool), neighbours, spacing_rad,
        max_layers, max(1, 2 * len(gap_points)), claimable=near_gap,
    )

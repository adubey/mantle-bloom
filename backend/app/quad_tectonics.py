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
  for a line end. An oceanic plate also carves out contested patches the peel can't reach
  from its edge (the line engine's interior-subduction carve-out, `_carve_interior`).
- **Advance** activates the empty cell across each exposed side of an eligible boundary
  cell, wherever no neighbour already covers it, again in layers. Ordinary new ground is a
  rift opening (`_open_rift`): the share of each new cell's footprint that lines up with the
  plates' separation is covered by stretching the crust within
  `K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION` cells behind it, volume-conservingly -- the 2D
  form of the line engine's `_stretch_end` -- and the rest is fresh magmatic oceanic crust.
  Columns stretched through `RIFT_CRITICAL_THICKNESS_M` erupt (breakup / ridge accretion).
  An active continental margin grows juvenile arc crust instead (`ARC_MARGIN_SEED_*`).

Boundary classification and the per-node column physics (convergent/arc/divergent/transform
updates, melting, provenance) are shared with the line engine unchanged -- see
`lithosphere_plate.boundary_context` and `deform_columns`. Nothing here needs regularizing:
cells never drift off the lattice.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components
from scipy.spatial import cKDTree

from . import geometry, lithosphere, phase_budget, rheology, terrain_noise
from .elevation_lines import (
    COVERAGE_RADIUS_MULT,
    CRUST_TYPE_CONTINENTAL,
    CRUST_TYPE_OCEANIC,
    ELEV_CHANGE_NEW_CRUST,
    ELEV_CHANGE_RIFT,
    ELEV_CHANGE_SUBDUCTION_ARC,
    ELEV_CHANGE_VOLCANO,
    effective_is_continental_from_codes,
    line_spacing_rad,
)
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
)
from .plates import _INTERIOR_SUBDUCTION_MIN_RUN, _contested_by_any
from .surface_fields import SURFACE_FIELDS

if TYPE_CHECKING:
    from .sparse_quad_patch import PlateWithSparseQuadPatch
    from .world import World

# How many layers of cells a boundary may advance in one step -- the same per-step ceiling the
# line engine's row claim uses (MAX_CLAIM_ROWS_PER_STEP), now applying in every direction.
MAX_ADVANCE_LAYERS_PER_STEP = MAX_CLAIM_ROWS_PER_STEP

# The eruption rng stream keys `deform_columns` / `_open_rift` use in place of a line index.
# Distinct values keep the column pass's melting draw and each advance layer's draws from
# sharing one (seed, step, plate) stream: the advance uses `_ADVANCE_RNG_INDEX + 2 * layer` for
# new cells and the next value for the stretched cells behind them.
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
    if plate.crust_type == "oceanic" and remaining > 0:
        removed |= _carve_interior(plate, retreatable & ~removed, ~has_probe | removed[safe_probes], remaining)
    if not np.any(removed):
        return survivors

    own_points = ctx.own_points
    world.record_removed_points(own_points[removed], plate.plate_id)
    donors = removed & ctx.accrete
    if np.any(donors):
        _accrete_onto_survivors(plate, donors, ~removed)
    plate.remove_cells(removed)
    return ~removed


def _carve_interior(plate: "PlateWithSparseQuadPatch", retreatable: np.ndarray, open_half: np.ndarray, max_cells: int) -> np.ndarray:
    """Interior subduction: the `retreatable` patches the layered peel can never reach,
    because no cell of theirs has a wholly open side (`open_half`: per cell, side and probe,
    whether that half-side borders nothing or an already-peeled cell -- the same whole-side
    test the peel uses, so a patch touching open ground only through a half-side still
    counts as unreachable) -- a neighbour overriding this plate somewhere other than its
    edge. Each such edge-connected patch of at least `_INTERIOR_SUBDUCTION_MIN_RUN` cells
    subducts, up to `max_cells` in total -- the 2D form of the line engine's mid-row
    carve-out, leaving a hole the quad surface represents directly. A patch larger than the
    remaining budget is carved partway, as a connected breadth-first prefix from one cell,
    so the hole it opens gives next step's peel an exposed edge to continue from. Oceanic
    plates only, as there: carving a continent's middle would sever it into a spurious
    defragmentation plate. Returns the mask to remove."""
    carved = np.zeros(len(retreatable), dtype=bool)
    members = np.flatnonzero(retreatable)
    if len(members) < _INTERIOR_SUBDUCTION_MIN_RUN or max_cells <= 0:
        return carved
    exposed = np.any(np.all(open_half[members], axis=2), axis=1)
    sub = _adjacency_matrix(plate)[members][:, members]
    _, labels = connected_components(sub, directed=False)
    sizes = np.bincount(labels)
    reachable = np.bincount(labels, weights=exposed) > 0
    budget = max_cells
    for label in np.flatnonzero((sizes >= _INTERIOR_SUBDUCTION_MIN_RUN) & ~reachable):
        if budget <= 0:
            break
        patch = np.flatnonzero(labels == label)
        if len(patch) > budget:
            order = breadth_first_order(sub, patch[0], directed=False, return_predecessors=False)
            patch = order[:budget]
        carved[members[patch]] = True
        budget -= len(patch)
    return carved


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
        stretch_share = _stretch_share(plate, local, sources, gap_direction, arc)
        fields = _new_cell_fields(plate, world, world_pts, arc, hc0, hm0, amp, texture)
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

        rifted = ~arc[inserted]
        if np.any(rifted):
            _open_rift(plate, world, _ADVANCE_RNG_INDEX + 2 * layer, new_index, new_index[rifted], stretch_share[inserted][rifted])
    return max_cells - budget


def _stretch_share(
    plate: "PlateWithSparseQuadPatch",
    local: np.ndarray,
    sources: np.ndarray,
    gap_direction: np.ndarray,
    arc: np.ndarray,
) -> np.ndarray:
    """Per candidate cell, the share of its footprint the plate covers by stretching its own
    crust rather than by fresh magmatic accretion: how much of the outward step from its
    source runs along the separation direction (toward the nearest neighbour node -- 1 with
    no neighbour in view, where the whole step is the plate pulling apart). The 2D form of
    the line engine's `rheology.stretch_components` split between end-stretch and row claim.
    Arc cells are fed by the slab, not by stretching: 0."""
    source_local = plate.surface_nodes().local_xyz[sources]
    outward = geometry.normalize(local - source_local)
    toward = np.zeros_like(outward)
    has_gap = np.linalg.norm(gap_direction, axis=1) > 0.0
    if np.any(has_gap):
        toward[has_gap] = geometry.normalize(geometry.to_local(plate.frame, gap_direction[has_gap]))
    alignment = np.where(has_gap, np.abs(np.sum(outward * toward, axis=1)), 1.0)
    return np.where(arc, 0.0, np.clip(alignment, 0.0, 1.0))


def _new_cell_fields(
    plate: "PlateWithSparseQuadPatch",
    world: "World",
    world_pts: np.ndarray,
    arc: np.ndarray,
    hc0: float,
    hm0: float,
    amp: float,
    texture: "terrain_noise.FractalTexture",
) -> dict[str, np.ndarray]:
    """Field values each candidate cell is inserted with. Arc cells are final: juvenile arc
    crust (`ARC_MARGIN_SEED_*`). Every other cell starts as the fresh magmatic column -- the
    oceanic growth seed plus the plate's terrain texture, explicitly typed oceanic -- which
    `_open_rift` then blends with the stretched crust it draws from behind."""
    count = len(world_pts)
    fields = {name: np.full(count, SURFACE_FIELDS[name].default, dtype=SURFACE_FIELDS[name].dtype) for name in (
        "elevation", "crustal_thickness_m", "mantle_lithosphere_thickness_m", "crust_type_code",
        "elev_change_reason", "node_created_years",
    )}
    fields["node_created_years"][:] = world.elapsed_years
    ordinary = ~arc
    if np.any(ordinary):
        hc = hc0 + amp * texture.sample(world_pts[ordinary])
        fields["crustal_thickness_m"][ordinary] = hc
        fields["mantle_lithosphere_thickness_m"][ordinary] = hm0
        fields["crust_type_code"][ordinary] = CRUST_TYPE_OCEANIC
        fields["elevation"][ordinary] = lithosphere.isostatic_elevation(hc, np.full(len(hc), hm0), lithosphere.RHO_OCEANIC_CRUST)
        fields["elev_change_reason"][ordinary] = ELEV_CHANGE_NEW_CRUST
    if np.any(arc):
        density = lithosphere.crust_density(plate.crust_type)
        fields["crustal_thickness_m"][arc] = ARC_MARGIN_SEED_HC_M
        fields["mantle_lithosphere_thickness_m"][arc] = ARC_MARGIN_SEED_HM_M
        fields["elevation"][arc] = lithosphere.isostatic_elevation(np.array([ARC_MARGIN_SEED_HC_M]), np.array([ARC_MARGIN_SEED_HM_M]), density)[0]
        fields["elev_change_reason"][arc] = ELEV_CHANGE_SUBDUCTION_ARC
    return fields


def _open_rift(
    plate: "PlateWithSparseQuadPatch",
    world: "World",
    rng_index: int,
    layer: np.ndarray,
    rifted: np.ndarray,
    stretch_share: np.ndarray,
) -> None:
    """Rift opening for this layer's `rifted` cells (node indices; `layer` is every cell this
    layer inserted). Each new cell covers `stretch_share` of its footprint by stretching the
    older cells within `K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION` hops of it, and the rest with
    the fresh magmatic column it was inserted with.

    Stretching is exactly volume-conserving (the areal form of `rheology.
    apply_stretch_thinning`, the line engine's `_stretch_end`): a donor asked to cover `D` m^2
    of new ground on top of its own area `a` thins by `a / (a + D)`, and the Hc/Hm it loses
    moves into the cells that asked, in proportion to what each asked of it. A new cell's
    demand is spread over its donor band by donor area. So a continental margin pulled apart
    thins into a widening band of stretched continental crust instead of losing it, and only
    once a column thins past `RIFT_CRITICAL_THICKNESS_M` does it erupt -- continental breakup
    or ridge accretion, through the same `_erupt_melted_nodes` path the column pass uses.
    Magmatic share is new crust from the mantle, as at any spreading ridge."""
    n = plate.node_count()
    areas = plate.node_areas_m2()
    old = np.ones(n, dtype=bool)
    old[layer] = False
    old_idx = np.flatnonzero(old)
    adjacency = _adjacency_matrix(plate).astype(float)
    reach = adjacency[rifted][:, old_idx]
    band = reach.copy()
    old_adjacency = adjacency[old_idx][:, old_idx]
    for _ in range(K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION - 1):
        reach = reach @ old_adjacency
        band = band + reach
    band = (band > 0).astype(float).multiply(areas[old_idx][None, :]).tocsr()
    band_area = np.asarray(band.sum(axis=1)).ravel()
    new_area = areas[rifted]
    demand = np.where(band_area > 0.0, stretch_share * new_area, 0.0)
    # asked[c, d]: footprint new cell c asks of donor d.
    asked = csr_matrix(band.multiply((demand / np.where(band_area > 0.0, band_area, 1.0))[:, None]))
    donor_demand = np.asarray(asked.sum(axis=0)).ravel()
    donor_area = areas[old_idx]
    ratio = donor_area / (donor_area + donor_demand)
    share_of_donor = csr_matrix(asked.multiply((1.0 / np.where(donor_demand > 0.0, donor_demand, 1.0))[None, :]))

    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    codes = plate.collect("crust_type_code")
    is_volcano = plate.collect("is_volcano")
    remaining = plate.collect("volcano_active_years_remaining")
    elevation = plate.collect("elevation")
    reason = plate.collect("elev_change_reason")
    continental = effective_is_continental_from_codes(codes, plate.crust_type == "continental")

    lost_hc = hc[old_idx] * donor_area * (1.0 - ratio)
    lost_hm = hm[old_idx] * donor_area * (1.0 - ratio)
    got_hc = share_of_donor @ lost_hc
    got_hm = share_of_donor @ lost_hm
    got_continental = share_of_donor @ (lost_hc * continental[old_idx])

    # Donors: thin in place, erupting any column that thins through the rift threshold.
    donors = old_idx[donor_demand > 0.0]
    donor_ratio = ratio[donor_demand > 0.0]
    if len(donors):
        before = lithosphere.isostatic_elevation(hc[donors], hm[donors], lithosphere.node_crust_density(codes[donors], plate.crust_type))
        new_hc, new_hm = hc[donors] * donor_ratio, hm[donors] * donor_ratio
        melting = (hc[donors] >= rheology.RIFT_CRITICAL_THICKNESS_M) & (new_hc < rheology.RIFT_CRITICAL_THICKNESS_M)
        sub_codes, sub_volcano, sub_remaining = codes[donors], is_volcano[donors], remaining[donors]
        _erupt_melted_nodes(world, plate.plate_id, rng_index + 1, new_hc, new_hm, sub_codes, sub_volcano, sub_remaining, melting, elevation[donors])
        after = lithosphere.isostatic_elevation(new_hc, new_hm, lithosphere.node_crust_density(sub_codes, plate.crust_type))
        elevation[donors] = rheology.clip_elevation_bounds(elevation[donors] + (after - before))
        hc[donors], hm[donors] = new_hc, new_hm
        codes[donors], is_volcano[donors], remaining[donors] = sub_codes, sub_volcano, sub_remaining
        reason[donors[melting]] = ELEV_CHANGE_VOLCANO

    # New cells: stretched crust from behind plus the magmatic remainder.
    magmatic = (1.0 - stretch_share) * new_area
    cell_hc = (got_hc + magmatic * hc[rifted]) / new_area
    cell_hm = (got_hm + magmatic * hm[rifted]) / new_area
    cell_codes = np.where(got_continental > 0.5 * cell_hc * new_area, CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC).astype(codes.dtype)
    cell_elevation = lithosphere.isostatic_elevation(cell_hc, cell_hm, lithosphere.node_crust_density(cell_codes, plate.crust_type))
    melting = cell_hc < rheology.RIFT_CRITICAL_THICKNESS_M
    sub_volcano, sub_remaining = is_volcano[rifted], remaining[rifted]
    _erupt_melted_nodes(world, plate.plate_id, rng_index, cell_hc, cell_hm, cell_codes, sub_volcano, sub_remaining, melting, cell_elevation)
    hc[rifted], hm[rifted], codes[rifted] = cell_hc, cell_hm, cell_codes
    is_volcano[rifted], remaining[rifted] = sub_volcano, sub_remaining
    elevation[rifted] = lithosphere.isostatic_elevation(cell_hc, cell_hm, lithosphere.node_crust_density(cell_codes, plate.crust_type))
    reason[rifted] = np.where(melting, ELEV_CHANGE_VOLCANO, np.where(stretch_share >= 0.5, ELEV_CHANGE_RIFT, ELEV_CHANGE_NEW_CRUST))
    plate.set_fields_on_plate(
        crustal_thickness_m=hc,
        mantle_lithosphere_thickness_m=hm,
        crust_type_code=codes,
        is_volcano=is_volcano,
        volcano_active_years_remaining=remaining,
        elevation=elevation,
        elev_change_reason=reason,
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

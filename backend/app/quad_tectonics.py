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
  conserved), as `_redistribute_accreted_column` does for a line end. Past the receiving
  belts it escapes along strike before any of it may delaminate (issue #290; `orogeny.py`). An oceanic plate also carves out contested patches the peel can't reach
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

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components
from scipy.spatial import cKDTree

from . import continental_ledger, cratons, geometry, lithosphere, orogeny, phase_budget, rheology, terrain_noise
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
    _ignite_early_rift_volcanoes,
    boundary_context,
    deform_columns,
    growth_seed_thickness,
)
from .plates import _INTERIOR_SUBDUCTION_MIN_RUN, _contested_by_any
from .sparse_quad_patch import lattice_points, unpack_cell_keys
from .surface_fields import CRATON_UNFORMED_YEARS, SURFACE_FIELDS

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

# Gap fill is allowed to bridge the sub-cell sliver between independently rotated quad
# lattices, but not by claiming a cell whose centre happens to be free while most of its
# footprint already belongs to a neighbour. Four samples per axis avoid making that decision
# depend on the centre point that the frontier's ordinary openness test already checks.
_GAP_FILL_SAMPLE_FRACTIONS = np.array([0.125, 0.375, 0.625, 0.875])

# A saturated suture may carry overflow through three additional belts of the same width:
# the broad orogen's own lateral spreading (issue #253).
SUTURE_ACCRETION_MAX_HOPS = 4 * SUTURE_ACCRETION_SPREAD_NODES
# Past the belts, tectonic escape (issue #290): crust squeezed out of a full orogen moves
# along strike, collision-parallel, as Indochina and Anatolia did, not straight inland. It
# reaches this many hops farther, into cells within `SUTURE_ESCAPE_MAX_ANGLE_DEG` of the
# suture's strike as seen from the front.
SUTURE_ESCAPE_EXTRA_HOPS = 4 * SUTURE_ACCRETION_SPREAD_NODES
SUTURE_ESCAPE_MAX_ANGLE_DEG = 45.0
# Only then may part of the donation delaminate, and only from the belts' eligible dense roots
# (`orogeny.delamination_capacity_m`). At most this share of any one donation, which also
# bounds cumulative suture delamination to this share of all donated crust.
SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION = 0.20


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
    continental_ledger.ensure_initialized(world)
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
    before = phase_budget.snapshot(plate, spacing_rad) if world.debug_diagnostics else None
    survivors = _retreat(plate, world, ctx, max_distance, max_cells, years)
    if world.debug_diagnostics:
        after = phase_budget.snapshot(plate, spacing_rad)
        phase_budget.record_snapshots(world, plate, "boundary_retreat", before, after)
        before = after

    if not ctx.suppress_growth:
        # Wide enough to see every plate a MAX_ADVANCE_LAYERS_PER_STEP-deep advance could run
        # into, not just the ones inside the boundary-force reach.
        growth_neighbours = plate.get_neighbours(
            other_plates, threshold_rad=ctx.reach_rad + (MAX_ADVANCE_LAYERS_PER_STEP + 1) * spacing_rad
        )
        _advance(plate, world, ctx, survivors, growth_neighbours, spacing_rad, max_cells)
        if world.debug_diagnostics:
            after = phase_budget.snapshot(plate, spacing_rad)
            phase_budget.record_snapshots(world, plate, "boundary_advance", before, after)
            before = after

    # Standing orogens collapse and flow into their surroundings (issue #290), making room
    # for next step's suture accretion before any of it may delaminate.
    orogeny.relax_orogens(plate, world, years)
    if world.debug_diagnostics:
        phase_budget.record_snapshots(world, plate, "orogenic_relief", before, phase_budget.snapshot(plate, spacing_rad))


def _retreat(
    plate: "PlateWithSparseQuadPatch", world: "World", ctx, max_distance: float, max_cells: int, years: float = 0.0
) -> np.ndarray:
    """Peel retreatable boundary cells, layer by layer -- see the module docstring. Returns the
    survivor mask over the pre-retreat node order. `years` sets how much of the receiving
    belts' dense roots may delaminate (none at 0)."""
    n = plate.node_count()
    survivors = np.ones(n, dtype=bool)
    # A craton cell holds out against consumption until it has been overlapped long enough
    # (cratons.retreat_allowed), so the younger belts around it take the shortening first.
    retreatable = ctx.shrinkable & cratons.retreat_allowed(world, plate)
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
    codes = plate.collect("crust_type_code")
    continental = effective_is_continental_from_codes(codes, plate.crust_type == "continental")
    # A continental terrane keeps being continental even when it rides an oceanic plate.
    # Treat those cells like the explicit continent-continent suture donors instead of
    # subducting them with their nominal owning plate (issue #253).
    terrane = continental if plate.crust_type == "oceanic" else np.zeros_like(continental)
    donors = removed & (ctx.accrete | terrane)
    # Everything else removed goes down the trench: its craton and its continental-derived
    # material leave the surface as deep subduction.
    subducted = removed & ~donors
    areas = plate.node_areas_m2()
    cratons.record(world, "subducted_m3", float(np.dot(plate.collect("craton_crust_m")[subducted], areas[subducted])))
    continental_ledger.record(
        world, "deeply_subducted_m3", float(np.dot(plate.collect("continental_material_m")[subducted], areas[subducted]))
    )
    if np.any(donors):
        _accrete_onto_survivors(
            plate, donors, ~removed, world, convergence_xyz=ctx.inputs.direction_to_neighbor, years=years,
            overriders=ctx.neighbours,
        )
    plate.remove_cells(removed)
    return ~removed


def _carve_interior(plate: "PlateWithSparseQuadPatch", retreatable: np.ndarray, open_half: np.ndarray, max_cells: int) -> np.ndarray:
    """Interior subduction -- the 2D form of the line engine's mid-row carve-out: remove
    contested patches the layered peel can never reach, leaving a hole the quad surface
    represents directly. Returns the mask to remove.

    Eligibility. A patch is an edge-connected component of `retreatable` cells (contested by
    a neighbour overriding this plate somewhere other than its edge) with at least
    `_INTERIOR_SUBDUCTION_MIN_RUN` cells. Oceanic plates only, as in the line engine; the
    caller never passes a continental plate, since carving a continent's middle would sever
    it into a spurious defragmentation plate.

    Reachability. A patch is left to the peel when any of its cells has a wholly open side.
    `open_half` is per cell, side and half-side probe: whether that half-side borders nothing
    or an already-peeled cell. A side is open only when both halves are, the same whole-side
    test the peel uses, so a patch touching open ground only through a half-side is
    unreachable and is carved here.

    Partial carve. At most `max_cells` cells go in total. A patch larger than what remains of
    that budget is carved partway, as a breadth-first prefix grown from one of its cells, so
    the removed part is connected and the hole it opens gives next step's peel an exposed
    edge to continue from."""
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


def _accrete_onto_survivors(
    plate: "PlateWithSparseQuadPatch",
    donors: np.ndarray,
    survivors: np.ndarray,
    world: "World | None" = None,
    convergence_xyz: np.ndarray | None = None,
    years: float = 0.0,
    overriders: list | None = None,
) -> None:
    """Thrust each continental suture's consumed Hc/Hm volume onto the surviving cells within
    `SUTURE_ACCRETION_SPREAD_NODES` hops behind it, as a uniform thickening -- the 2D form of
    `_redistribute_accreted_column`, which spreads a retreating line end's volume over that
    many nodes *inward along the row*. Each edge-connected run of donor cells is one suture
    front with its own band, so separate sutures on one plate don't share volume. Donors and
    receivers are matched by effective crust type, which keeps continental terranes on
    nominally oceanic plates in the continental reservoir. A terrane on an oceanic plate
    relocates onto oceanic footprint when its reservoir has no room.

    When the first band fills, Hc goes where `_place_suture_crust` finds room, in order
    (issue #290): three more belts (lateral spreading), cells along the suture's strike
    (tectonic escape; strike from `convergence_xyz`, the world-frame direction to the
    overriding plate per node, or the front's own long axis), the belts' eligible dense
    roots (delamination, rate-limited over `years`), the nearest capacity farther across the
    plate, then the plate's own non-continental cells (a foreland spill, retyping a cell
    continental once continental crust makes up most of it), then the overriding plate among
    `overriders` (`_hand_to_overrider`), and only what none of those can hold, a terminal
    delamination remainder.
    Hm spreads through the belts and delaminates past them, as an over-thickened mantle
    root does.

    Each front's continental-derived material (`continental_material_m`) travels with the Hc
    it placed and its delaminated share is booked to the continental ledger; its cratonic
    crust becomes ordinary orogenic crust (craton ledger: collision reworked) or delaminates.
    A relocated terrane keeps both whole. Booking needs `world`; without one only the fields
    move."""
    survivor_idx = np.flatnonzero(survivors)
    if not len(survivor_idx):
        return
    areas = plate.node_areas_m2()
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    elevation = plate.collect("elevation")
    codes = plate.collect("crust_type_code")
    codes_before = codes.copy()
    continental = effective_is_continental_from_codes(codes, plate.crust_type == "continental")
    adjacency = _adjacency_matrix(plate)
    points = plate.surface_nodes().local_xyz
    changed = np.zeros(len(hc), dtype=bool)
    elevation_before = elevation.copy()
    hc_before = hc.copy()
    hm_before = hm.copy()
    material = plate.collect("continental_material_m")
    craton = plate.collect("craton_crust_m")
    formed = plate.collect("craton_formed_years")
    # Dense-root volume each surviving continental column may shed this step, drawn down
    # as fronts use it so two fronts sharing a belt can't both spend it.
    root_capacity = np.where(
        survivors & continental, orogeny.delamination_capacity_m(hc, hm, years / 1_000_000.0) * areas, 0.0
    )
    convergence_local = (
        geometry.to_local(plate.frame, np.asarray(convergence_xyz, dtype=float)) if convergence_xyz is not None else None
    )

    for donor_type in (False, True):
        typed_donors = donors & (continental == donor_type)
        donor_idx = np.flatnonzero(typed_donors)
        if not len(donor_idx):
            continue
        _, labels = connected_components(adjacency[donor_idx][:, donor_idx], directed=False)
        for label in np.unique(labels):
            front = donor_idx[labels == label]
            hc_front_start = hc.copy()
            material_volume = float(np.dot(material[front], areas[front]))
            craton_volume = float(np.dot(craton[front], areas[front]))
            typed_survivors = survivors & (continental == donor_type)
            hc_volume = float(np.sum(hc[front] * areas[front]))
            hm_volume = float(np.sum(hm[front] * areas[front]))
            hc_room = float(
                np.sum(
                    np.maximum(SUTURE_ACCRETION_MAX_HC_M - hc[typed_survivors], 0.0)
                    * areas[typed_survivors]
                )
            )
            hm_room = float(
                np.sum(
                    np.maximum(lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M - hm[typed_survivors], 0.0)
                    * areas[typed_survivors]
                )
            )
            needs_new_footprint = hc_room < hc_volume * (1.0 - 1e-12) or hm_room < hm_volume * (1.0 - 1e-12)
            if (
                plate.crust_type == "oceanic"
                and donor_type
                and needs_new_footprint
                and np.any(survivors & ~continental)
            ):
                # The retreat consumed an isolated terrane, or every remaining continental
                # receiver is full. Move this front's column onto the nearest surviving
                # oceanic footprint, displacing that footprint's old oceanic column back into
                # the oceanic reservoir. This preserves both categorical Hc volume and total
                # Hc instead of either losing the fragment or counting its oceanic substrate
                # as newly continental.
                was_continental = continental.copy()
                relocated = _relocate_terrane_column(
                    hc, hm, codes, continental, areas, points, adjacency, front, survivors, hc_volume, hm_volume
                )
                if np.any(relocated):
                    changed |= relocated
                    # The terrane's new footprint is the cells relocation retyped; the rest of
                    # `relocated` received the oceanic columns it displaced.
                    _relocate_tracers(
                        material, craton, formed, areas, hc - hc_front_start, front, continental & ~was_continental,
                        material_volume, craton_volume,
                    )
                    continue
            if not np.any(typed_survivors):
                # Preserve the old any-type nearest-survivor fallback. Same-type placement is
                # preferred because it preserves the categorical reservoir, but the column's
                # actual Hc/Hm volume must not vanish when that representation has no receiver.
                typed_survivors = survivors
            strike = _suture_strike(points, front, convergence_local)
            shed = np.zeros(len(hc))
            filled, stages = _place_suture_crust(
                hc,
                areas,
                adjacency,
                points,
                front,
                typed_survivors,
                hc_volume,
                SUTURE_ACCRETION_MAX_HC_M,
                strike,
                root_capacity if donor_type else None,
                shed,
            )
            # The shed roots take their own columns' continental material and cratonic
            # crust with them, in proportion to the column they left.
            shed_material, shed_craton = _shed_provenance(material, craton, hc_front_start, shed)
            changed |= filled
            stuck = stages["no_outlet_delaminated_m3"]
            spill_cells = survivors & ~continental
            if donor_type and stuck > 0.0 and np.any(spill_cells):
                # Every continental receiver is full -- typically a small plate whose whole
                # continent sits at the cap, where collapse has no gradient to work with.
                # The crust thrusts out over the plate's own drowned margin or oceanic
                # foreland instead, nearest first; a cell left mostly continental crust is
                # retyped continental.
                hc_before_spill = hc.copy()
                spilled, spill_stages = _place_suture_crust(
                    hc, areas, adjacency, points, front, spill_cells, stuck, SUTURE_ACCRETION_MAX_HC_M, None, None
                )
                changed |= spilled
                left = spill_stages["no_outlet_delaminated_m3"]
                stages["foreland_spill_placed_m3"] = stuck - left
                stages["no_outlet_delaminated_m3"] = left
                retyped = spilled & (hc - hc_before_spill > 0.5 * hc)
                codes[retyped] = CRUST_TYPE_CONTINENTAL
                continental[retyped] = True
            handed_share = 0.0
            stuck = stages["no_outlet_delaminated_m3"]
            if donor_type and stuck > 0.0 and overriders and hc_volume > 0.0:
                # The plate is being consumed faster than it can hold its own crust -- a
                # microcontinent ground into a collision. What it can't keep accretes to the
                # plate overriding it, as real suture crust does.
                share = stuck / hc_volume
                handed = _hand_to_overrider(
                    world, plate, front, overriders, stuck, share * material_volume, convergence_xyz, years,
                )
                handed_share = handed / hc_volume
                stages["overrider_placed_m3"] = handed
                stages["no_outlet_delaminated_m3"] = max(stuck - handed, 0.0)
            placed_share = _carry_material(material, areas, hc - hc_front_start + shed, hc_volume, material_volume)
            lost_share = max(1.0 - placed_share - handed_share, 0.0)
            if world is not None:
                continental_ledger.record(
                    world, "delaminated_lower_crust_m3", lost_share * material_volume + float(shed_material @ areas)
                )
                cratons.record(world, "collision_reworked_m3", (1.0 - lost_share) * craton_volume)
                cratons.record(world, "delaminated_m3", lost_share * craton_volume + float(shed_craton @ areas))
                orogeny.record(world, "suture_donated_m3", hc_volume)
                for account, volume in stages.items():
                    orogeny.record(world, account, volume)
            changed |= _spread_accretion_volume(
                hm,
                areas,
                adjacency,
                front,
                typed_survivors,
                hm_volume,
                lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M,
            )

    gained = np.flatnonzero(changed)
    density_before = lithosphere.node_crust_density(codes_before[gained], plate.crust_type)
    density_after = lithosphere.node_crust_density(codes[gained], plate.crust_type)
    before = lithosphere.isostatic_elevation(hc_before[gained], hm_before[gained], density_before)
    after = lithosphere.isostatic_elevation(hc[gained], hm[gained], density_after)
    elevation[gained] = rheology.clip_elevation_bounds(elevation_before[gained] + (after - before))
    plate.set_fields_on_plate(
        crustal_thickness_m=hc,
        mantle_lithosphere_thickness_m=hm,
        crust_type_code=codes,
        elevation=elevation,
        continental_material_m=material,
        craton_crust_m=craton,
        craton_formed_years=formed,
    )


def _shed_provenance(
    material: np.ndarray, craton: np.ndarray, hc_before: np.ndarray, shed: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """In place: remove the continental material and cratonic crust that delaminated roots
    (`shed`, thickness per cell) take with them -- each cell's own fractions of the column it
    had before this front arrived. Returns the removed thicknesses (material, craton)."""
    fraction = np.divide(shed, hc_before, out=np.zeros(len(shed)), where=(shed > 0.0) & (hc_before > 0.0))
    shed_material = material * np.minimum(fraction, 1.0)
    shed_craton = craton * np.minimum(fraction, 1.0)
    material -= shed_material
    craton -= shed_craton
    return shed_material, shed_craton


def _carry_material(
    material: np.ndarray, areas: np.ndarray, hc_gain: np.ndarray, hc_volume: float, material_volume: float
) -> float:
    """In place: give each receiver the front's continental-derived material in proportion to
    the Hc volume it received (`hc_gain`, thickness). Returns the share of the front's Hc
    that was placed rather than delaminated."""
    gained = np.maximum(hc_gain, 0.0)
    placed = float(np.dot(gained, areas))
    if hc_volume <= 0.0 or placed <= 0.0:
        return 0.0
    share = min(placed / hc_volume, 1.0)
    material += share * material_volume * gained / placed
    return share


def _relocate_tracers(
    material: np.ndarray,
    craton: np.ndarray,
    formed: np.ndarray,
    areas: np.ndarray,
    hc_change: np.ndarray,
    front: np.ndarray,
    chosen: np.ndarray,
    material_volume: float,
    craton_volume: float,
) -> None:
    """In place: a relocated terrane (`_relocate_terrane_column`) carries its continental
    material, cratonic crust and craton date onto the `chosen` cells; the oceanic columns it
    displaced take their own material with them to the cells that received their Hc."""
    chosen_area = float(areas[chosen].sum())
    displaced = float(np.dot(material[chosen], areas[chosen]))
    gained = np.where(chosen, 0.0, np.maximum(hc_change, 0.0))
    placed = float(np.dot(gained, areas))
    if placed > 0.0:
        material += displaced * gained / placed
    material[chosen] = material_volume / chosen_area
    craton[chosen] = craton_volume / chosen_area
    dated = formed[front][formed[front] != CRATON_UNFORMED_YEARS]
    formed[chosen] = float(dated.min()) if len(dated) and craton_volume > 0.0 else CRATON_UNFORMED_YEARS


def _relocate_terrane_column(
    hc: np.ndarray,
    hm: np.ndarray,
    codes: np.ndarray,
    continental: np.ndarray,
    areas: np.ndarray,
    points: np.ndarray,
    adjacency: csr_matrix,
    front: np.ndarray,
    survivors: np.ndarray,
    hc_volume: float,
    hm_volume: float,
) -> np.ndarray:
    """Re-home a fully consumed continental fragment without mixing reservoir labels.

    The nearest oceanic survivor cells become the fragment and receive exactly its Hc/Hm
    volumes. Their displaced oceanic columns spread back into the remaining oceanic plate.
    """
    oceanic = survivors & ~continental
    changed = np.zeros(len(hc), dtype=bool)
    if not np.any(oceanic):
        return changed

    reached = np.zeros(len(hc), dtype=bool)
    reached[front] = True
    frontier = reached.copy()
    chosen = np.zeros(len(hc), dtype=bool)
    chosen_area = 0.0
    while (
        chosen_area * SUTURE_ACCRETION_MAX_HC_M < hc_volume
        or chosen_area * lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M < hm_volume
    ):
        next_frontier = (adjacency @ frontier.astype(np.int8) > 0) & ~reached
        if not np.any(next_frontier):
            # A plate can be briefly disconnected before the topology cleanup pass. Preserve
            # the fragment on its nearest reachable component when possible, then take only
            # as much of the nearest other component(s) as the terrane column actually needs.
            remaining = np.flatnonzero(oceanic & ~chosen)
            if not len(remaining):
                break
            # `front` and candidates are already in one plate-local frame. Chord distance is
            # monotonic with spherical distance, so sorting by it avoids building another
            # spatial index for this rare cleanup path.
            centre = points[front].mean(axis=0)
            order = np.argsort(np.linalg.norm(points[remaining] - centre, axis=1), kind="stable")
            for idx in remaining[order]:
                chosen[idx] = True
                chosen_area += float(areas[idx])
                if (
                    chosen_area * SUTURE_ACCRETION_MAX_HC_M >= hc_volume
                    and chosen_area * lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M >= hm_volume
                ):
                    break
            break
        frontier = next_frontier
        reached |= frontier
        layer = frontier & oceanic
        chosen |= layer
        chosen_area += float(areas[layer].sum())
    if not np.any(chosen):
        return changed
    if (
        chosen_area * SUTURE_ACCRETION_MAX_HC_M < hc_volume * (1.0 - 1e-12)
        or chosen_area * lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M < hm_volume * (1.0 - 1e-12)
    ):
        return changed

    displaced_hc = float(np.sum(hc[chosen] * areas[chosen]))
    displaced_hm = float(np.sum(hm[chosen] * areas[chosen]))
    oceanic_receivers = survivors & ~continental & ~chosen
    hc_room = float(
        np.sum(np.maximum(SUTURE_ACCRETION_MAX_HC_M - hc[oceanic_receivers], 0.0) * areas[oceanic_receivers])
    )
    hm_room = float(
        np.sum(
            np.maximum(lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M - hm[oceanic_receivers], 0.0)
            * areas[oceanic_receivers]
        )
    )
    if hc_room < displaced_hc * (1.0 - 1e-12) or hm_room < displaced_hm * (1.0 - 1e-12):
        return changed
    hc[chosen] = min(hc_volume / chosen_area, SUTURE_ACCRETION_MAX_HC_M)
    hm[chosen] = min(hm_volume / chosen_area, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)
    codes[chosen] = CRUST_TYPE_CONTINENTAL
    continental[chosen] = True
    changed |= chosen

    if np.any(oceanic_receivers):
        changed |= _spread_accretion_volume(
            hc,
            areas,
            adjacency,
            np.flatnonzero(chosen),
            oceanic_receivers,
            displaced_hc,
            SUTURE_ACCRETION_MAX_HC_M,
            0.0,
        )
        changed |= _spread_accretion_volume(
            hm,
            areas,
            adjacency,
            np.flatnonzero(chosen),
            oceanic_receivers,
            displaced_hm,
            lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M,
            0.0,
        )
    return changed


def _spread_accretion_volume(
    thickness: np.ndarray,
    areas: np.ndarray,
    adjacency: csr_matrix,
    front: np.ndarray,
    eligible: np.ndarray,
    volume: float,
    cap: float,
    max_delamination_fraction: float = 1.0,
) -> np.ndarray:
    """Spread one reservoir outward from a donor front, carrying cap overflow onward.

    The ordinary suture width remains the first receiving band. Only overflow reaches cells
    farther inland, one graph-hop layer at a time, through three additional suture-width belts.
    Within each accumulated band the added thickness is uniform except where a cell hits
    `cap`; water-filling then gives the rest to cells with remaining room. Anything the four
    belts cannot hold delaminates, up to `max_delamination_fraction` of this donation; larger
    overflow must find room elsewhere in the same reservoir. A suture's mantle lithosphere
    keeps the full-delamination default, and a relocated terrane's displaced oceanic columns
    use 0. A suture's crust goes through `_place_suture_crust` instead (issue #290).
    Returns the cells whose thickness changed.
    """
    changed = np.zeros(len(thickness), dtype=bool)
    if volume <= 0.0 or not np.any(eligible):
        return changed

    reached = np.zeros(len(thickness), dtype=bool)
    reached[front] = True
    frontier = reached.copy()
    receivers = np.zeros(len(thickness), dtype=bool)
    hops = 0
    remaining = volume
    placed_locally = False
    far_fallback = False
    while remaining > max(volume, 1.0) * 1e-12:
        if hops > 0:
            receivers |= frontier & eligible
        if hops >= SUTURE_ACCRETION_SPREAD_NODES and np.any(receivers):
            before_fill = remaining
            remaining, filled = _fill_capped_volume(thickness, areas, receivers, remaining, cap)
            changed |= filled
            placed_locally |= remaining < before_fill
            receivers[:] = False
        if hops >= SUTURE_ACCRETION_MAX_HOPS:
            # If the finite belt found no receiver with room, retain the old unbounded search
            # instead of treating mere distance as delamination. Once that fallback starts it
            # walks every farther graph layer until the column is placed or the graph ends.
            far_fallback |= not placed_locally
            if not far_fallback:
                must_place = max(remaining - max_delamination_fraction * volume, 0.0)
                if must_place > 0.0:
                    _, filled = _fill_capped_volume(thickness, areas, eligible, must_place, cap)
                    changed |= filled
                break
        next_frontier = (adjacency @ frontier.astype(np.int8) > 0) & ~reached
        if not np.any(next_frontier):
            # Preserve the old nearest-survivor fallback for a temporarily disconnected
            # plate, but only after every graph-reachable band has had first call.
            if remaining > max(volume, 1.0) * 1e-12:
                remaining, filled = _fill_capped_volume(thickness, areas, eligible, remaining, cap)
                changed |= filled
            break
        reached |= next_frontier
        frontier = next_frontier
        hops += 1
    return changed


def _suture_strike(points: np.ndarray, front: np.ndarray, convergence_local: np.ndarray | None) -> np.ndarray | None:
    """Unit tangent along a suture front's strike, in the plate's local frame: perpendicular to
    the front's mean direction toward the overriding plate when that is known, else the
    front's own long axis. None when neither defines one (a single cell with no direction)."""
    centre = geometry.normalize(points[front].mean(axis=0))
    if convergence_local is not None:
        toward = convergence_local[front].mean(axis=0)
        toward = toward - np.dot(toward, centre) * centre
        if np.linalg.norm(toward) > 1e-9:
            return geometry.normalize(np.cross(centre, toward))
    if len(front) < 2:
        return None
    offsets = points[front] - points[front].mean(axis=0)
    offsets = offsets - np.outer(offsets @ centre, centre)
    _, singular, axes = np.linalg.svd(offsets, full_matrices=False)
    if singular[0] <= 1e-12 or (len(singular) > 1 and singular[1] >= 0.999 * singular[0]):
        return None  # no long axis: a compact front has no strike
    return axes[0]


def _hand_to_overrider(
    world: "World | None",
    plate: "PlateWithSparseQuadPatch",
    front: np.ndarray,
    overriders: list,
    volume: float,
    material_volume: float,
    convergence_xyz: np.ndarray | None,
    years: float,
) -> float:
    """Accrete `volume` of a consumed front's Hc onto the quad plate overriding it, through the
    same staged `_place_suture_crust` placement on that plate, seeded at its cells nearest the
    front. The overrider is the candidate containing most of the front's cell centres (the
    nearest one if none does). The front's continental material goes with the Hc it placed, in
    proportion; its cratonic crust becomes ordinary orogenic crust, which the caller books. The
    overrider's own belts may shed eligible roots to make room, booked here with their own
    provenance. Whatever the overrider can't hold is left for the caller's terminal remainder.
    Returns the Hc volume placed on the overrider."""
    candidates = [p for p in overriders if p is not plate and hasattr(p, "adjacency") and p.node_count() > 0]
    if not candidates or volume <= 0.0:
        return 0.0
    front_world = geometry.to_world(plate.frame, plate.surface_nodes().local_xyz[front])
    centre = geometry.normalize(front_world.mean(axis=0))

    def rank(candidate) -> tuple[int, float]:
        inside = int(np.count_nonzero(candidate.contains_batch(front_world)))
        nearest = float(cKDTree(candidate.all_points_and_elevation()[0]).query(centre)[0])
        return -inside, nearest

    over = min(candidates, key=rank)
    over_world = over.all_points_and_elevation()[0]
    _, seed = cKDTree(over_world).query(front_world)
    seed = np.unique(seed)

    areas = over.node_areas_m2()
    hc = over.collect("crustal_thickness_m")
    hm = over.collect("mantle_lithosphere_thickness_m")
    codes = over.collect("crust_type_code")
    material = over.collect("continental_material_m")
    craton = over.collect("craton_crust_m")
    elevation = over.collect("elevation")
    continental = effective_is_continental_from_codes(codes, over.crust_type == "continental")
    eligible = continental.copy() if np.any(continental) else np.ones(len(hc), dtype=bool)
    root_capacity = np.where(continental, orogeny.delamination_capacity_m(hc, hm, years / 1_000_000.0) * areas, 0.0)
    points = over.surface_nodes().local_xyz
    convergence_local = None
    if convergence_xyz is not None:
        # The overrider converges on the front from the other side.
        toward = -np.asarray(convergence_xyz, dtype=float)[front].mean(axis=0)
        convergence_local = np.zeros_like(points)
        convergence_local[seed] = geometry.to_local(over.frame, toward[None, :])[0]
    strike = _suture_strike(points, seed, convergence_local)

    hc_before = hc.copy()
    shed = np.zeros(len(hc))
    changed, stages = _place_suture_crust(
        hc, areas, _adjacency_matrix(over), points, seed, eligible, volume, SUTURE_ACCRETION_MAX_HC_M, strike,
        root_capacity, shed,
    )
    placed = max(volume - stages["no_outlet_delaminated_m3"], 0.0)
    if placed <= 0.0:
        return 0.0
    shed_material, shed_craton = _shed_provenance(material, craton, hc_before, shed)
    _carry_material(material, areas, hc - hc_before + shed, volume, material_volume)
    if world is not None:
        continental_ledger.record(world, "delaminated_lower_crust_m3", float(shed_material @ areas))
        cratons.record(world, "delaminated_m3", float(shed_craton @ areas))
        for account in (
            "suture_belt_placed_m3", "escape_attempted_m3", "escape_placed_m3", "delamination_attempted_m3",
            "delamination_completed_m3", "far_field_placed_m3",
        ):
            orogeny.record(world, account, stages[account])
    density = lithosphere.node_crust_density(codes, over.crust_type)
    gained = np.flatnonzero(changed)
    shift = lithosphere.isostatic_elevation(hc[gained], hm[gained], density[gained]) - lithosphere.isostatic_elevation(
        hc_before[gained], hm[gained], density[gained]
    )
    elevation[gained] = rheology.clip_elevation_bounds(elevation[gained] + shift)
    over.set_fields_on_plate(crustal_thickness_m=hc, continental_material_m=material, craton_crust_m=craton, elevation=elevation)
    return placed


def _place_suture_crust(
    thickness: np.ndarray,
    areas: np.ndarray,
    adjacency: csr_matrix,
    points: np.ndarray,
    front: np.ndarray,
    eligible: np.ndarray,
    volume: float,
    cap: float,
    strike: np.ndarray | None,
    root_capacity: np.ndarray | None,
    shed: np.ndarray | None = None,
) -> tuple[np.ndarray, dict[str, float]]:
    """Place one front's donated Hc `volume` under `cap`, in the order of issue #290:

    1. **Belts** -- graph-hop bands up to `SUTURE_ACCRETION_MAX_HOPS` from the front, water-
       filled band by band once the first `SUTURE_ACCRETION_SPREAD_NODES` hops are reached.
    2. **Tectonic escape** -- the next `SUTURE_ESCAPE_EXTRA_HOPS` hops, but only cells within
       `SUTURE_ESCAPE_MAX_ANGLE_DEG` of `strike` as seen from the front. Skipped without one.
    3. **Delamination** -- the belts shed part of their own dense roots, drawn from
       `root_capacity` (volume per cell, drawn down in place) up to
       `SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION` of `volume`, and the donation takes the
       room that frees. What sinks is the receivers' old lower crust, not the incoming crust:
       each cell's shed thickness is added to `shed` for the caller to book with that cell's
       own provenance. None when `root_capacity` is None.
    4. **Far field** -- every eligible cell reached so far, then farther hops one at a time,
       then (a plate the graph doesn't connect) any eligible cell.
    5. **No outlet** -- whatever still has no room among `eligible`. The caller spills a
       continental front's remainder onto the plate's non-continental cells before booking
       any of it as terminal delamination.

    Returns the cells whose thickness changed and the volume each stage handled, keyed by
    `orogeny.BUDGET_ACCOUNTS` name."""
    changed = np.zeros(len(thickness), dtype=bool)
    stages = dict.fromkeys(
        (
            "suture_belt_placed_m3",
            "escape_attempted_m3",
            "escape_placed_m3",
            "delamination_attempted_m3",
            "delamination_completed_m3",
            "far_field_placed_m3",
            "foreland_spill_placed_m3",
            "overrider_placed_m3",
            "no_outlet_delaminated_m3",
        ),
        0.0,
    )
    if volume <= 0.0 or not np.any(eligible):
        return changed, stages
    tolerance = max(volume, 1.0) * 1e-12
    remaining = volume

    def fill(receivers: np.ndarray) -> float:
        nonlocal remaining
        before = remaining
        remaining, filled = _fill_capped_volume(thickness, areas, receivers, remaining, cap)
        changed[:] |= filled
        return before - remaining

    reached = np.zeros(len(thickness), dtype=bool)
    reached[front] = True
    frontier = reached.copy()
    hops = 0

    def step() -> bool:
        nonlocal frontier, hops
        next_frontier = (adjacency @ frontier.astype(np.int8) > 0) & ~reached
        if not np.any(next_frontier):
            return False
        reached[:] |= next_frontier
        frontier = next_frontier
        hops += 1
        return True

    # 1. Belts.
    belts = np.zeros(len(thickness), dtype=bool)
    receivers = np.zeros(len(thickness), dtype=bool)
    graph_open = True
    while remaining > tolerance and hops < SUTURE_ACCRETION_MAX_HOPS:
        graph_open = step()
        if not graph_open:
            break
        receivers |= frontier & eligible
        belts |= frontier & eligible
        if hops >= SUTURE_ACCRETION_SPREAD_NODES and np.any(receivers):
            stages["suture_belt_placed_m3"] += fill(receivers)
            receivers[:] = False
    if remaining > tolerance and np.any(receivers):
        # The graph ended inside the first band.
        stages["suture_belt_placed_m3"] += fill(receivers)

    # 2. Tectonic escape along strike.
    if remaining > tolerance and graph_open and strike is not None:
        stages["escape_attempted_m3"] = remaining
        centre = geometry.normalize(points[front].mean(axis=0))
        offsets = points - centre
        offsets -= np.outer(offsets @ centre, centre)
        length = np.linalg.norm(offsets, axis=1)
        alignment = np.abs(offsets @ strike) / np.where(length > 0.0, length, 1.0)
        along_strike = alignment >= np.cos(np.radians(SUTURE_ESCAPE_MAX_ANGLE_DEG))
        escape_end = hops + SUTURE_ESCAPE_EXTRA_HOPS
        while remaining > tolerance and hops < escape_end:
            graph_open = step()
            if not graph_open:
                break
            stages["escape_placed_m3"] += fill(frontier & eligible & along_strike)

    # 3. Delamination from the belts' eligible dense roots.
    if remaining > tolerance and root_capacity is not None:
        stages["delamination_attempted_m3"] = remaining
        available = float(root_capacity[belts].sum())
        take = min(remaining, available, SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION * volume)
        if take > 0.0:
            lost = root_capacity * np.where(belts, take / available, 0.0)
            root_capacity -= lost
            thickness -= lost / areas
            if shed is not None:
                shed += lost / areas
            changed |= lost > 0.0
            stages["delamination_completed_m3"] = take
            fill(belts)

    # 4. Far field: nearest remaining capacity first.
    if remaining > tolerance:
        stages["far_field_placed_m3"] += fill(reached & eligible)
        while remaining > tolerance and graph_open:
            graph_open = step()
            if graph_open:
                stages["far_field_placed_m3"] += fill(frontier & eligible)
        if remaining > tolerance:
            stages["far_field_placed_m3"] += fill(eligible)

    # 5. No outlet anywhere on the plate.
    if remaining > tolerance:
        stages["no_outlet_delaminated_m3"] = remaining
    return changed, stages


def _fill_capped_volume(
    thickness: np.ndarray,
    areas: np.ndarray,
    receivers: np.ndarray,
    volume: float,
    cap: float,
) -> tuple[float, np.ndarray]:
    """Water-fill `volume` across receivers, returning (unplaced volume, changed mask)."""
    changed = np.zeros(len(thickness), dtype=bool)
    idx = np.flatnonzero(receivers & (thickness < cap))
    if not len(idx) or volume <= 0.0:
        return volume, changed
    room = np.maximum(cap - thickness[idx], 0.0)
    take = min(volume, float(np.dot(room, areas[idx])))
    if take <= 0.0:
        return volume, changed

    # Find the common added thickness whose capped per-cell volumes sum to `take` despite
    # unequal quad-cell areas.
    lo, hi = 0.0, float(room.max())
    for _ in range(52):
        mid = 0.5 * (lo + hi)
        if float(np.dot(np.minimum(room, mid), areas[idx])) < take:
            lo = mid
        else:
            hi = mid
    add = np.minimum(room, hi)
    thickness[idx] += add
    changed[idx[add > 0.0]] = True
    placed = min(float(np.dot(add, areas[idx])), volume)
    return volume - placed, changed


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
    standoff: bool = True,
) -> int:
    """Activate open empty cells across the exposed sides of `eligible` cells, then across
    the newest layer's, for up to `max_layers` layers and `max_cells` cells -- the shared
    areal-growth walk behind boundary advance and gap filling. A candidate is open when no
    plate in `neighbours` contains it and, with `standoff`, none of their nodes lies within
    `EXTEND_THRESHOLD_MULTIPLIER` spacings; `claimable(keys, world_pts)`, when given, narrows
    that further. Cells grown from an `arc_source` cell are arc crust; the rest are fresh crust
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
            open_mask &= claimable(candidates, world_pts)
        gap_direction = np.zeros((len(candidates), 3))
        if neighbour_tree is not None:
            dist, idx = neighbour_tree.query(world_pts)
            if standoff:
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
        arc_added = arc[inserted]
        if np.any(arc_added):
            continental_ledger.record(
                world,
                "juvenile_additions_m3",
                float(np.dot(fields["continental_material_m"][inserted], plate.node_areas_m2()[new_index])),
            )
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
        "elev_change_reason", "node_created_years", "continental_material_m",
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
        fields["continental_material_m"][arc] = ARC_MARGIN_SEED_HC_M
    return fields


@dataclass(frozen=True)
class _StretchTransfer:
    """How a rift layer's stretched share is paid for, from `_allocate_stretch`. Donors are
    the pre-existing cells that cover some new cell's stretched footprint; every volume is in
    m * m^2 (thickness times area) and aligned to the rifted cells it is received by."""

    donor_indices: np.ndarray  # node indices of cells that thin
    donor_thinning_ratio: np.ndarray  # a / (a + D) per donor: its new thickness / its old
    received_hc_volume: np.ndarray  # Hc volume each rifted cell receives from its donors
    received_hm_volume: np.ndarray  # Hm volume each rifted cell receives
    received_continental_hc_volume: np.ndarray  # the part of received_hc_volume that was continental
    withheld_area: np.ndarray  # m^2 of each rifted cell's stretched share cratonic donors withheld
    received_continental_material_volume: np.ndarray


def _allocate_stretch(
    plate: "PlateWithSparseQuadPatch",
    inserted_indices: np.ndarray,
    rifted: np.ndarray,
    stretch_share: np.ndarray,
    hc: np.ndarray,
    hm: np.ndarray,
    continental: np.ndarray,
    continental_material: np.ndarray,
    craton: np.ndarray | None = None,
) -> _StretchTransfer:
    """Share each rifted cell's stretched footprint (`stretch_share` of its area) out over its
    donor band -- the pre-existing cells within `K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION` hops,
    weighted by donor area -- and work out what each donor loses and each rifted cell receives.

    A donor asked to cover `D` m^2 of new ground on top of its own area `a` thins by
    `a / (a + D)`; the volume that removes, `thickness * a * (1 - ratio)`, goes to the cells
    that asked, in proportion to how much each asked. Total volume is unchanged by
    construction. A rifted cell with no pre-existing cell in reach asks for nothing.

    A cratonic donor (`craton`, its cratonic thickness per node) is asked for only
    `1 - CRATON_STRETCH_RESISTANCE * strength` of its area share; what it withholds is not
    shifted onto the other donors but left to the rifted cell's magmatic share
    (`withheld_area`), so the stretch concentrates in the younger crust beside a craton."""
    areas = plate.node_areas_m2()
    preexisting_mask = np.ones(plate.node_count(), dtype=bool)
    preexisting_mask[inserted_indices] = False
    preexisting = np.flatnonzero(preexisting_mask)
    adjacency = _adjacency_matrix(plate).astype(float)
    reach = adjacency[rifted][:, preexisting]
    band = reach.copy()
    preexisting_adjacency = adjacency[preexisting][:, preexisting]
    for _ in range(K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION - 1):
        reach = reach @ preexisting_adjacency
        band = band + reach
    # band_area_by_donor[c, d]: donor d's area if it is in rifted cell c's band, else 0.
    in_band = (band > 0).astype(float)
    band_area_by_donor = in_band.multiply(areas[preexisting][None, :]).tocsr()
    band_area = np.asarray(band_area_by_donor.sum(axis=1)).ravel()
    requested_area_by_cell = np.where(band_area > 0.0, stretch_share * areas[rifted], 0.0)
    willing = np.ones(len(preexisting))
    if craton is not None:
        willing = 1.0 - cratons.CRATON_STRETCH_RESISTANCE * cratons.strength(craton[preexisting])
    # requested_area[c, d]: footprint rifted cell c asks of donor d, in m^2.
    requested_area = csr_matrix(
        in_band.multiply((areas[preexisting] * willing)[None, :]).multiply(
            (requested_area_by_cell / np.where(band_area > 0.0, band_area, 1.0))[:, None]
        )
    )
    requested_area_by_donor = np.asarray(requested_area.sum(axis=0)).ravel()
    donor_area = areas[preexisting]
    thinning_ratio = donor_area / (donor_area + requested_area_by_donor)
    # share_of_donor[c, d]: the fraction of donor d's lost volume rifted cell c receives.
    share_of_donor = csr_matrix(
        requested_area.multiply((1.0 / np.where(requested_area_by_donor > 0.0, requested_area_by_donor, 1.0))[None, :])
    )

    lost_hc_volume = hc[preexisting] * donor_area * (1.0 - thinning_ratio)
    lost_hm_volume = hm[preexisting] * donor_area * (1.0 - thinning_ratio)
    is_donor = requested_area_by_donor > 0.0
    return _StretchTransfer(
        donor_indices=preexisting[is_donor],
        donor_thinning_ratio=thinning_ratio[is_donor],
        received_hc_volume=share_of_donor @ lost_hc_volume,
        received_hm_volume=share_of_donor @ lost_hm_volume,
        received_continental_hc_volume=share_of_donor @ (lost_hc_volume * continental[preexisting]),
        withheld_area=requested_area_by_cell - np.asarray(requested_area.sum(axis=1)).ravel(),
        received_continental_material_volume=share_of_donor @ (
            continental_material[preexisting] * donor_area * (1.0 - thinning_ratio)
        ),
    )


def _open_rift(
    plate: "PlateWithSparseQuadPatch",
    world: "World",
    rng_index: int,
    inserted_indices: np.ndarray,
    rifted: np.ndarray,
    stretch_share: np.ndarray,
) -> None:
    """Rift opening for this layer's `rifted` cells (node indices; `inserted_indices` is every
    cell this layer inserted, arc cells included). Each rifted cell covers `stretch_share` of
    its footprint by stretching the pre-existing cells behind it (`_allocate_stretch`), and
    the rest with the fresh magmatic column it was inserted with. The areal form of
    `rheology.apply_stretch_thinning`, the line engine's `_stretch_end`, and exactly
    volume-conserving.

    So a continental margin pulled apart thins into a widening band of stretched continental
    crust instead of losing it, and only once a column thins past `RIFT_CRITICAL_THICKNESS_M`
    does it erupt -- continental breakup or ridge accretion, through the same
    `_erupt_melted_nodes` path the column pass uses. Magmatic share is new crust from the
    mantle, as at any spreading ridge. A cell that is mostly magmatic is a vent: it starts a
    volcano lifecycle without changing its column, as every node the line engine's row claims
    and gap filling create does (`seed_and_erupt_new_nodes`)."""
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    codes = plate.collect("crust_type_code")
    is_volcano = plate.collect("is_volcano")
    remaining = plate.collect("volcano_active_years_remaining")
    elevation = plate.collect("elevation")
    reason = plate.collect("elev_change_reason")
    craton = plate.collect("craton_crust_m")
    material = plate.collect("continental_material_m")
    continental = effective_is_continental_from_codes(codes, plate.crust_type == "continental")
    transfer = _allocate_stretch(
        plate, inserted_indices, rifted, stretch_share, hc, hm, continental, material, craton
    )

    # Donors: thin in place, erupting any column that thins through the rift threshold.
    donors, ratio = transfer.donor_indices, transfer.donor_thinning_ratio
    if len(donors):
        before = lithosphere.isostatic_elevation(hc[donors], hm[donors], lithosphere.node_crust_density(codes[donors], plate.crust_type))
        # Hm floored like `rheology.apply_stretch_thinning`: a donor too thin to melt through
        # is never reset, and repeated rifts would otherwise thin it without limit (issue #256).
        new_hc = hc[donors] * ratio
        new_hm = np.maximum(hm[donors] * ratio, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M)
        melting = (hc[donors] >= rheology.RIFT_CRITICAL_THICKNESS_M) & (new_hc < rheology.RIFT_CRITICAL_THICKNESS_M)
        sub_codes, sub_volcano, sub_remaining = codes[donors], is_volcano[donors], remaining[donors]
        material[donors] *= ratio
        hc_before_melting = new_hc.copy()
        _erupt_melted_nodes(world, plate.plate_id, rng_index + 1, new_hc, new_hm, sub_codes, sub_volcano, sub_remaining, melting, elevation[donors])
        donor_continental = effective_is_continental_from_codes(sub_codes, plate.crust_type == "continental")
        donor_juvenile = np.where(melting & donor_continental, np.maximum(new_hc - hc_before_melting, 0.0), 0.0)
        material[donors] += donor_juvenile
        continental_ledger.record(
            world,
            "juvenile_additions_m3",
            float(np.dot(donor_juvenile, plate.node_areas_m2()[donors])),
        )
        donor_oceanic = melting & ~donor_continental
        if np.any(donor_oceanic):
            removed = np.where(donor_oceanic, material[donors], 0.0)
            material[donors[donor_oceanic]] = 0.0
            continental_ledger.record(
                world, "rift_thinned_m3",
                float(np.dot(removed, plate.node_areas_m2()[donors])),
            )
        after = lithosphere.isostatic_elevation(new_hc, new_hm, lithosphere.node_crust_density(sub_codes, plate.crust_type))
        elevation[donors] = rheology.clip_elevation_bounds(elevation[donors] + (after - before))
        hc[donors], hm[donors] = new_hc, new_hm
        codes[donors], is_volcano[donors], remaining[donors] = sub_codes, sub_volcano, sub_remaining
        reason[donors[melting]] = ELEV_CHANGE_VOLCANO
        # A donor's cratonic crust thins with it, and is gone where it melted through.
        thinned = np.where(melting, 0.0, craton[donors] * ratio)
        cratons.record(world, "rifted_m3", float(np.dot(craton[donors] - thinned, plate.node_areas_m2()[donors])))
        craton[donors] = thinned

    # Rifted cells: stretched crust received from behind plus the magmatic remainder.
    cell_area = plate.node_areas_m2()[rifted]
    magmatic_area = (1.0 - stretch_share) * cell_area + np.maximum(transfer.withheld_area, 0.0)
    cell_hc = (transfer.received_hc_volume + magmatic_area * hc[rifted]) / cell_area
    cell_hm = np.maximum(
        (transfer.received_hm_volume + magmatic_area * hm[rifted]) / cell_area, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
    )
    mostly_continental = transfer.received_continental_hc_volume > 0.5 * cell_hc * cell_area
    cell_codes = np.where(mostly_continental, CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC).astype(codes.dtype)
    cell_material = transfer.received_continental_material_volume / cell_area
    magmatic_material = np.where(
        mostly_continental, magmatic_area * hc[rifted] / cell_area, 0.0
    )
    cell_material += magmatic_material
    cell_elevation = lithosphere.isostatic_elevation(cell_hc, cell_hm, lithosphere.node_crust_density(cell_codes, plate.crust_type))
    melting = cell_hc < rheology.RIFT_CRITICAL_THICKNESS_M
    sub_volcano, sub_remaining = is_volcano[rifted], remaining[rifted]
    hc_before_melting = cell_hc.copy()
    _erupt_melted_nodes(world, plate.plate_id, rng_index, cell_hc, cell_hm, cell_codes, sub_volcano, sub_remaining, melting, cell_elevation)
    cell_continental = effective_is_continental_from_codes(cell_codes, plate.crust_type == "continental")
    melt_juvenile = np.where(melting & cell_continental, np.maximum(cell_hc - hc_before_melting, 0.0), 0.0)
    cell_material += melt_juvenile
    continental_ledger.record(
        world,
        "juvenile_additions_m3",
        float(np.dot(magmatic_material + melt_juvenile, cell_area)),
    )
    cell_oceanic = melting & ~cell_continental
    if np.any(cell_oceanic):
        removed = np.where(cell_oceanic, cell_material, 0.0)
        cell_material[cell_oceanic] = 0.0
        continental_ledger.record(
            world, "rift_thinned_m3", float(np.dot(removed, cell_area))
        )
    _ignite_early_rift_volcanoes(world, plate.plate_id, rng_index, sub_volcano, sub_remaining, ~melting & (stretch_share < 0.5))
    hc[rifted], hm[rifted], codes[rifted] = cell_hc, cell_hm, cell_codes
    material[rifted] = cell_material
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
        craton_crust_m=craton,
        continental_material_m=material,
    )


def _gap_cells_mostly_uncovered(
    plate: "PlateWithSparseQuadPatch", keys: np.ndarray, others: list
) -> np.ndarray:
    """Whether more than half of each candidate cell's sampled footprint is unclaimed."""
    keys = np.asarray(keys, dtype=np.int64).reshape(-1)
    if len(keys) == 0 or not others:
        return np.ones(len(keys), dtype=bool)
    face, level, i, j = unpack_cell_keys(keys)
    scale = np.left_shift(1, level)
    u, v = np.meshgrid(_GAP_FILL_SAMPLE_FRACTIONS, _GAP_FILL_SAMPLE_FRACTIONS, indexing="xy")
    u, v = u.ravel(), v.ravel()
    a = (i[:, None] + u) / scale[:, None]
    b = (j[:, None] + v) / scale[:, None]
    local = lattice_points(
        np.broadcast_to(face[:, None], a.shape).ravel(),
        a.ravel(),
        b.ravel(),
        plate._n,
    )
    world_points = geometry.to_world(plate.frame, local)
    covered = _contested_by_any(world_points, others).reshape(len(keys), len(u))
    return np.mean(~covered, axis=1) > 0.5


def fill_gap(
    world: "World", plate: "PlateWithSparseQuadPatch", gap_points: np.ndarray, others: list, spacing_rad: float, max_layers: int
) -> int:
    """Grow `plate` into the uncovered `gap_points` (world xyz) -- the quad counterpart of
    `gap_fill_frontier.fill_gap_by_growing_plates` for one claimant. The same frontier walk as
    boundary advance, from every boundary cell, restricted to cells whose centre lies within
    `COVERAGE_RADIUS_MULT` spacings of a gap point and that no plate in `others` contains.
    A 4x4 interior sample must also show that most of the candidate footprint is uncovered;
    checking only its centre lets a rotated neighbour already own most of the cell (#268).

    Unlike boundary advance there is no standoff from the neighbours' nodes. Advance keeps a
    new cell's centre `EXTEND_THRESHOLD_MULTIPLIER` spacings from them, which is what leaves
    a seam about one cell wide between two plates that stopped short of each other; gap
    filling is what closes it. Returns how many cells were added."""
    if plate.node_count() == 0 or len(gap_points) == 0:
        return 0
    gap_tree = cKDTree(gap_points)
    coverage_radius_rad = COVERAGE_RADIUS_MULT * spacing_rad

    def claimable(keys: np.ndarray, world_pts: np.ndarray) -> np.ndarray:
        dist, _ = gap_tree.query(world_pts)
        near = dist <= coverage_radius_rad
        if np.any(near):
            near[near] &= _gap_cells_mostly_uncovered(plate, keys[near], neighbours)
        return near

    neighbours = plate.get_neighbours(others, threshold_rad=(max_layers + 1) * spacing_rad)
    n = plate.node_count()
    return grow_frontier(
        plate, world, np.ones(n, dtype=bool), np.zeros(n, dtype=bool), neighbours, spacing_rad,
        max_layers, max(1, 2 * len(gap_points)), claimable=claimable, standoff=False,
    )

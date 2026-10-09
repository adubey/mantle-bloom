"""Consumed lower-plate crust transfers to the upper plate at a polarized collision front
(issue #320, part of #315).

At a continental collision front with a frozen polarity (collision_polarity.py), only the
lower plate's cells retreat, and their mantle lithosphere goes down with the slab
(`suture_hm_subducted_m3`, issue #319). This module decides where their crust goes. It
replaces the old treatment, in which each plate thrust its own consumed crust back onto its
own survivors.

**Partition.** Each connected front of consumed lower-plate cells is one column, summed
with exact cell areas. Its mobile cover, at the top of the column, goes with the scraped
share. The crust below the cover splits three ways (`partition`):

- **Scraped** (`World.suture_scrape_fraction`): upper and middle crust, peeled off into
  thrust sheets stacked in front of the upper plate. The Himalaya is largely Indian crust
  stacked south of the suture. It goes through the staged suture placement
  (`quad_tectonics._place_on_overrider`): belts, escape along strike, the upper plate's own
  dense roots shedding to make room (booked with *their* provenance), then far field.
- **Underthrust** (`World.suture_underthrust_fraction`): lower crust pushed horizontally
  under the upper plate's front, as Indian lower crust lies under southern Tibet (Nabelek et
  al. 2009). It is added to the base of the upper plate's continental columns within
  `UNDERTHRUST_REACH_KM` of the front, water-filled under the Hc cap.
- **Lost** (`World.suture_lower_crust_loss_fraction`): lower crust that eclogitizes and
  sinks with the slab. Ultrahigh-pressure rocks show continental crust going deeper than
  100 km (Chopin 1984).

These are model knobs for a partition that varies along strike and is unresolved. India-
Asia mass balances that put the loss near half the converged crust (Ingalls et al. 2016)
depend on the reconstruction and have been challenged, so they are not a calibration
target. The defaults keep the loss small, close to issue #276's collision-loss calibration,
since about 3% of the continental inventory can pass through suture retreat every Myr
(issue #289).

**Caps and no outlet**, decided before anything moves:

1. Underthrust crust the reach can't hold under the Hc cap is thrust up with the scraped
   share instead: the orogen accommodates it at the surface rather than losing it.
2. Scraped crust the upper plate can't hold anywhere goes to the front's other overriding
   neighbours with continental crust (`quad_tectonics._hand_to_overrider`), as the old
   fallback did. Oceanic neighbours never take it: their cells aren't retyped.
3. What none of them can hold goes down with the slab: `no_outlet_subducted_m3`, and
   `collision_subducted_m3` for its continental material, as before.

**Provenance.** Continental material, restite and cratonic crust split in proportion to
each share's Hc volume, the same source-volume fractions. The mobile cover shares the fate
of the thrust-up crust it tops: `accreted_m3` where that lands, `subducted_m3` in proportion
to what had no outlet. Material and restite move with
the Hc each receiver takes. A continental lower plate's cratonic crust that reaches the upper
plate becomes ordinary orogenic crust (`collision_reworked_m3`), as suture crust always has;
a docking terrane's lands as craton (see "Terranes"). The lost share's
material goes to `collision_lower_crust_subducted_m3` in the continental-material ledger,
its craton to `subducted_m3` and its restite to `restite_subducted_m3`. Receivers' Moho is
buried under the crust they take (`orogeny.bury_moho`), and their elevation follows the
actual Hc they gained, isostatically.

**Mantle lithosphere** stays out of it: none goes onto either plate. A coherent underthrust
share of lower-plate Hm (Indian lithosphere under Tibet) is real, but adding it would
restore the Hm ratchet before #310 has a relaxation sink, so it is later scope, at zero.

**Interface.** `transfer_column` takes a `ConsumedColumn` -- volumes and a world-frame
front, nothing tied to the lower plate's node order -- and an upper plate.

**Terranes** (issue #321). A continental terrane riding an oceanic plate is consumed against
a continental overrider when its carrier subducts there. It docks onto that plate through the
same partition and placement: it jams the trench and welds onto the margin, as Wrangellia
did onto North America (Coney et al. 1980). Its upper plate is the frozen front's when the
contact is a polarized front, and otherwise the plate overriding it at its contact, where
that plate's crust is continental: a continent's margin or another terrane (`dock_targets`). Its cratonic crust stays cratonic: it moves with each share of the
crust it is part of, and each receiver keeps the older of its own craton date and the
terrane's. A craton docks the way cratons assembled, as blocks welded together without
losing their age. Its mantle lithosphere goes down with the carrier's slab, booked to
the carrier's oceanic subduction sink (`oceanic_and_deep_subduction`); none of it becomes an
upper-plate Hm stack. Unlike a continental lower plate, a terrane's crust never goes down for
lack of room: what no overrider can hold is returned to the carrier's own cells (`returned_m3`,
the column's provenance and Hm scaled to it), and the caller's own-plate path places or
relocates it (`quad_tectonics._relocate_terrane_column`). That path sinks the returned Hm
too, unless the terrane relocates and needs a mantle column under its new footprint. Docked, partial and fallback fronts
are counted in `World.suture_transfer_stats["terranes"]`.

**Roots.** Every front landing on one upper plate in a step draws down one shared allowance
of root shedding (`root_capacity`), so several fronts can't each shed a full step's worth.

**Validation.** `partition` is checked once at the start of each step's deform pass, so bad
knobs fail before any plate deforms; `summary` reports them without raising.

**Cost.** Per front, not per cell: one query on the upper plate's cached KD-tree, and one
cell graph per upper plate, shared by all its fronts' underthrust and scraped placement.
With `World.debug_diagnostics`, `World.suture_transfer_stats` counts fronts, volumes, time
per front and per step, in constant size.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse.csgraph import connected_components

from . import continental_ledger, cratons, hm_ledger, orogeny, quad_tectonics
from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from .hydroclimate import mobile_cover
from .lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M
from .surface_fields import CRATON_UNFORMED_YEARS

if TYPE_CHECKING:
    from .sparse_quad_patch import PlateWithSparseQuadPatch
    from .world import World

# Indian lower crust is imaged under southern Tibet to roughly 200 km north of the suture
# (Nabelek et al. 2009). Converted to graph hops with the world's line spacing.
UNDERTHRUST_REACH_KM = 200.0
# The three shares must sum to 1 within this.
PARTITION_SUM_TOLERANCE = 1e-9


@dataclass(frozen=True)
class CrustPartition:
    """Shares of a consumed column's crust below its mobile cover. They sum to 1."""

    scrape: float
    underthrust: float
    loss: float


def partition(world: "World") -> CrustPartition:
    """The world's configured partition, validated."""
    shares = CrustPartition(
        scrape=float(world.suture_scrape_fraction),
        underthrust=float(world.suture_underthrust_fraction),
        loss=float(world.suture_lower_crust_loss_fraction),
    )
    values = (shares.scrape, shares.underthrust, shares.loss)
    if not all(np.isfinite(v) and v >= 0.0 for v in values):
        raise ValueError(f"suture crust partition shares must be finite and non-negative, got {shares}")
    if abs(sum(values) - 1.0) > PARTITION_SUM_TOLERANCE:
        raise ValueError(f"suture crust partition shares must sum to 1, got {shares}")
    return shares


@dataclass(frozen=True)
class ConsumedColumn:
    """One consumed front's column, as volumes (m^3, exact cell areas) and a world-frame
    footprint. `approach_world` is the direction the upper plate converges on the front
    from, or None. `terrane`: a continental terrane on an oceanic carrier, which docks
    rather than losing crust for lack of room (module docstring)."""

    plate_id: int
    front_world: np.ndarray
    approach_world: np.ndarray | None
    hc_m3: float
    hm_m3: float
    mobile_cover_m3: float
    material_m3: float
    craton_m3: float
    restite_m3: float
    terrane: bool = False
    # The oldest craton date on the front (CRATON_UNFORMED_YEARS without one). A terrane's
    # craton lands with this date; a continental lower plate's is reworked.
    craton_formed_years: float = CRATON_UNFORMED_YEARS


@dataclass(frozen=True)
class TransferResult:
    """Where one column's Hc went (m^3). scraped + underthrust + lost = the column's Hc;
    scraped + underthrust overflow = upper placed + handed + unplaced + returned. Only a
    terrane returns crust, and only a continental lower plate leaves any unplaced."""

    scraped_m3: float
    underthrust_m3: float
    underthrust_placed_m3: float
    lost_m3: float
    upper_placed_m3: float
    handed_m3: float
    unplaced_m3: float
    returned_m3: float = 0.0


def column_of(
    plate: "PlateWithSparseQuadPatch", front: np.ndarray, approach_world: np.ndarray | None
) -> ConsumedColumn:
    """The column a front of `plate`'s cells (indices) carries, before they're removed. On an
    oceanic plate the front is a continental terrane."""
    areas = plate.node_areas_m2()[front]

    def volume(name: str) -> float:
        return float(np.dot(plate.collect(name)[front], areas))

    return ConsumedColumn(
        plate_id=plate.plate_id,
        front_world=plate.all_points_and_elevation()[0][front],
        approach_world=approach_world,
        hc_m3=volume("crustal_thickness_m"),
        hm_m3=volume("mantle_lithosphere_thickness_m"),
        mobile_cover_m3=volume("mobile_cover_m"),
        material_m3=volume("continental_material_m"),
        craton_m3=volume("craton_crust_m"),
        restite_m3=volume("restite_m"),
        terrane=plate.crust_type == "oceanic",
        craton_formed_years=_oldest_craton_date(plate, front),
    )


def _oldest_craton_date(plate: "PlateWithSparseQuadPatch", front: np.ndarray) -> float:
    """The oldest craton date among `front`'s cells that hold craton."""
    dated = front[plate.collect("craton_crust_m")[front] > 0.0]
    return float(plate.collect("craton_formed_years")[dated].min()) if len(dated) else CRATON_UNFORMED_YEARS


def dock_targets(
    world: "World",
    plate: "PlateWithSparseQuadPatch",
    terrane: np.ndarray,
    upper_plate_ids: np.ndarray,
    neighbour_plate_ids: np.ndarray,
    neighbour_node_index: np.ndarray,
    neighbours: list,
) -> np.ndarray:
    """Each consumed terrane cell's upper plate (issue #321): the frozen front's in
    `upper_plate_ids` where it has one, otherwise the plate overriding it at its contact --
    its nearest other plate (`neighbour_plate_ids`), when that plate's nearest node
    (`neighbour_node_index`) is continental-coded. That is a continent's margin, or another
    terrane: terranes amalgamate offshore before they reach a continent, as Wrangellia and
    the Alexander terrane did. -1 where there is none: an oceanic-coded contact leaves the
    terrane to its carrier's own-plate path. `terrane` masks the consumed continental cells
    of oceanic `plate`."""
    out = np.asarray(upper_plate_ids, dtype=np.int64).copy()
    if plate.crust_type != "oceanic":
        return out
    if world.debug_diagnostics:
        _terrane_stats(world)["polarized_cells"] += int(np.count_nonzero(terrane & (out >= 0)))
    open_cells = np.flatnonzero(terrane & (out < 0) & (neighbour_plate_ids >= 0))
    if not len(open_cells):
        return out
    by_id = {p.plate_id: p for p in neighbours if p is not plate and p.node_count() > 0}
    for q in np.unique(neighbour_plate_ids[open_cells]):
        over = by_id.get(int(q))
        if over is None:
            continue
        cells = open_cells[neighbour_plate_ids[open_cells] == q]
        nodes = neighbour_node_index[cells]
        continental = effective_is_continental_from_codes(over.collect("crust_type_code"), over.crust_type == "continental")
        valid = (nodes >= 0) & (nodes < len(continental))
        cells, nodes = cells[valid], nodes[valid]
        docking = cells[continental[nodes]]
        out[docking] = int(q)
        if world.debug_diagnostics and over.crust_type == "oceanic":
            _terrane_stats(world)["terrane_to_terrane_cells"] += len(docking)
    if world.debug_diagnostics:
        _terrane_stats(world)["contact_docking_cells"] += int(np.count_nonzero(out[open_cells] >= 0))
    return out


def transfer_fronts(
    world: "World",
    plate: "PlateWithSparseQuadPatch",
    donors: np.ndarray,
    upper_plate_ids: np.ndarray,
    convergence_xyz: np.ndarray | None,
    neighbour_plate_ids: np.ndarray | None,
    years: float,
    overriders: list | None,
    can_return: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Transfer every connected front of `donors` (this plate's consumed lower-plate cells,
    or a terrane's, each with its upper plate in `upper_plate_ids`) to its upper plate. A
    front whose upper plate is gone is left for the caller's ordinary suture accretion. Each
    front's approach direction averages `convergence_xyz` (per node, toward its nearest
    other plate) over the nodes whose nearest plate in `neighbour_plate_ids` is the upper
    plate, so a third plate at a junction doesn't set it. Returns (handled, returned): the
    donors whose column went wholly to their upper plate, and terrane donors that docked only
    in part, whose fields now hold the returned remainder for the caller's own-plate path.
    Without `can_return` (the carrier has no surviving cells to keep it on) nothing returns:
    a terrane's remainder goes down with the slab like a continental lower plate's."""
    handled = np.zeros(len(donors), dtype=bool)
    returned = np.zeros(len(donors), dtype=bool)
    candidates = donors & (upper_plate_ids >= 0)
    if not np.any(candidates):
        return handled, returned
    shares = partition(world)
    by_id = {p.plate_id: p for p in world.plates}
    adjacency = quad_tectonics._adjacency_matrix(plate)
    for upper_id in np.unique(upper_plate_ids[candidates]):
        upper = by_id.get(int(upper_id))
        if upper is None or upper is plate or upper.node_count() == 0 or not hasattr(upper, "adjacency"):
            continue
        upper_adjacency = quad_tectonics._adjacency_matrix(upper)
        cells = np.flatnonzero(candidates & (upper_plate_ids == upper_id))
        _, labels = connected_components(adjacency[cells][:, cells], directed=False)
        for label in np.unique(labels):
            front = cells[labels == label]
            facing = front
            if neighbour_plate_ids is not None and np.any(neighbour_plate_ids[front] == upper_id):
                facing = front[neighbour_plate_ids[front] == upper_id]
            column = column_of(plate, front, quad_tectonics._overrider_approach(convergence_xyz, facing))
            result = transfer_column(
                world, column, upper, overriders=overriders, years=years, shares=shares, adjacency=upper_adjacency,
                can_return=can_return,
            )
            if result.returned_m3 > 0.0:
                _keep_remainder(plate, front, column, result)
                returned[front] = True
            else:
                handled[front] = True
    return handled, returned


def _keep_remainder(plate: "PlateWithSparseQuadPatch", front: np.ndarray, column: ConsumedColumn, result: TransferResult) -> None:
    """In place: scale a partly docked terrane front's cells down to what came back. Its Hc,
    Hm, continental material, craton and restite keep the returned share of the column's Hc.
    Its returned cover (`_returned_cover_m3`) is spread over the remaining crust as one
    fraction of it, so no cell's cover can exceed its Hc."""
    share = result.returned_m3 / column.hc_m3 if column.hc_m3 > 0.0 else 0.0
    fields = {}
    for name in ("crustal_thickness_m", "mantle_lithosphere_thickness_m", "continental_material_m", "craton_crust_m", "restite_m"):
        values = plate.collect(name)
        values[front] *= share
        fields[name] = values
    # Cover rides on thrust-up crust, so it is at most that crust: the fraction is <= 1.
    fraction = min(_returned_cover_m3(column, result) / result.returned_m3, 1.0) if result.returned_m3 > 0.0 else 0.0
    cover = plate.collect("mobile_cover_m")
    cover_total = float(cover[front] @ plate.node_areas_m2()[front])
    continental = plate.collect("mobile_cover_continental_m")
    continental_share = float(continental[front] @ plate.node_areas_m2()[front]) / cover_total if cover_total > 0.0 else 0.0
    cover[front] = fraction * fields["crustal_thickness_m"][front]
    continental[front] = continental_share * cover[front]
    fields["mobile_cover_m"] = cover
    fields["mobile_cover_continental_m"] = continental
    plate.set_fields_on_plate(**fields)


def _thrust_up(result: TransferResult) -> float:
    """The crust thrust up at the front: the scraped share and the underthrust overflow."""
    return result.scraped_m3 + result.underthrust_m3 - result.underthrust_placed_m3


def _returned_cover_m3(column: ConsumedColumn, result: TransferResult) -> float:
    """The cover a terrane's returned crust carries: the returned share of the thrust-up
    crust the cover rode on. `_book` and `_keep_remainder` both use it."""
    thrust_up = _thrust_up(result)
    return max(column.mobile_cover_m3, 0.0) * min(result.returned_m3 / thrust_up, 1.0) if thrust_up > 0.0 else 0.0


def root_capacity(world: "World", upper: "PlateWithSparseQuadPatch", years: float) -> np.ndarray:
    """The dense-root volume per cell `upper` may still shed this step, shared by every front
    that lands on it (drawn down in place by the placement). Kept in the deform pass's
    step-scoped `World.suture_root_capacity`, and recomputed if `upper`'s cells have changed
    since; outside a step, a fresh allowance."""
    cache = getattr(world, "suture_root_capacity", None)
    entry = cache.get(upper.plate_id) if cache is not None else None
    if entry is not None and np.array_equal(entry[0], upper.cell_keys):
        return entry[1]
    continental = effective_is_continental_from_codes(upper.collect("crust_type_code"), upper.crust_type == "continental")
    capacity = orogeny.plate_delamination_capacity_m3(upper, continental, years)
    if cache is not None:
        cache[upper.plate_id] = (upper.cell_keys.copy(), capacity)
    return capacity


def _approach_from_seed(column: ConsumedColumn, upper: "PlateWithSparseQuadPatch", seed: np.ndarray) -> np.ndarray | None:
    """Without a known approach, the direction from `upper`'s seed cells to the front."""
    centre = column.front_world.mean(axis=0)
    toward = centre - upper.all_points_and_elevation()[0][seed].mean(axis=0)
    normal = centre / max(float(np.linalg.norm(centre)), 1e-12)
    toward -= np.dot(toward, normal) * normal
    return toward if np.linalg.norm(toward) > 1e-9 else None


def transfer_column(
    world: "World",
    column: ConsumedColumn,
    upper: "PlateWithSparseQuadPatch",
    *,
    overriders: list | None = None,
    years: float = 0.0,
    shares: CrustPartition | None = None,
    adjacency=None,
    can_return: bool = True,
) -> TransferResult:
    """Partition one consumed column and place its retained crust on `upper` -- see the
    module docstring. `overriders` are the front's other overriding neighbours, the fallback
    for crust `upper` can't hold; only those with continental crust take it. `adjacency`
    reuses `upper`'s cell graph. Books every share, the column's mobile cover included; the
    caller removes the source cells, or for a terrane that returned crust, keeps their
    remainder. Without `can_return` a terrane returns nothing."""
    started = time.perf_counter()
    shares = partition(world) if shares is None else shares
    hc = max(column.hc_m3, 0.0)
    body = max(hc - max(column.mobile_cover_m3, 0.0), 0.0)
    underthrust = shares.underthrust * body
    lost = shares.loss * body
    scraped = max(hc - underthrust - lost, 0.0)

    def share_of(volume: float) -> float:
        return volume / hc if hc > 0.0 else 0.0

    # A terrane's craton lands with its crust (and its date); a continental lower plate's is
    # reworked into orogenic crust (`_book`).
    craton = column.craton_m3 if column.terrane else 0.0

    if adjacency is None:
        adjacency = quad_tectonics._adjacency_matrix(upper)
    seed = quad_tectonics._overrider_seed(upper, column.front_world)
    approach = column.approach_world if column.approach_world is not None else _approach_from_seed(column, upper, seed)
    # The roots `upper` may shed this step, as it stood before this front's crust arrived.
    roots = root_capacity(world, upper, years)
    underthrust_placed = _underthrust(
        world, upper, seed, adjacency, underthrust,
        share_of(underthrust) * column.material_m3, share_of(underthrust) * column.restite_m3,
        share_of(underthrust) * craton, column.craton_formed_years,
    )
    # 1. What the reach can't hold under the cap is thrust up with the scraped share.
    thrust_up = scraped + (underthrust - underthrust_placed)
    upper_placed = quad_tectonics._place_on_overrider(
        world, upper, column.front_world, thrust_up, share_of(thrust_up) * column.material_m3,
        approach, years, share_of(thrust_up) * column.restite_m3, adjacency,
        seed=seed, root_capacity=roots,
        craton_volume=share_of(thrust_up) * craton, craton_formed_years=column.craton_formed_years,
    )
    # 2. Then the front's other overriding neighbours with continental crust.
    handed = 0.0
    remaining = thrust_up - upper_placed
    others = [
        p for p in (overriders or [])
        if p is not upper and p.plate_id != column.plate_id and p.node_count() > 0
        and quad_tectonics._has_continental_crust(p)
    ]
    if remaining > max(thrust_up, 1.0) * 1e-12 and others:
        handed = quad_tectonics._hand_to_overrider(
            world, upper, column.front_world, others, remaining, share_of(remaining) * column.material_m3,
            approach, years, share_of(remaining) * column.restite_m3,
            share_of(remaining) * craton, column.craton_formed_years,
        )
    # 3. The rest goes down with the slab -- or, from a terrane, back to its carrier.
    rest = max(remaining - handed, 0.0)
    returned = rest if column.terrane and can_return and rest > max(thrust_up, 1.0) * 1e-12 else 0.0
    unplaced = rest - returned

    result = TransferResult(
        scraped_m3=scraped,
        underthrust_m3=underthrust,
        underthrust_placed_m3=underthrust_placed,
        lost_m3=lost,
        upper_placed_m3=upper_placed,
        handed_m3=handed,
        unplaced_m3=unplaced,
        returned_m3=returned,
    )
    _book(world, column, upper, result, thrust_up)
    if world.debug_diagnostics:
        _record_stats(world, column, result, hc, time.perf_counter() - started)
    return result


def _underthrust(
    world: "World",
    upper: "PlateWithSparseQuadPatch",
    seed: np.ndarray,
    adjacency,
    volume: float,
    material_volume: float,
    restite_volume: float,
    craton_volume: float = 0.0,
    craton_formed_years: float = CRATON_UNFORMED_YEARS,
) -> float:
    """Add `volume` of lower crust to the base of `upper`'s continental columns within
    `UNDERTHRUST_REACH_KM` of `seed`, water-filled under the Hc cap, with its continental
    material and restite in proportion, and a carried `craton_volume` with its date
    (`quad_tectonics.settle_received_crust`). A plate
    with no continental columns takes none: continental crust isn't put under oceanic
    columns, and the caller thrusts it up instead. Returns the volume placed."""
    if volume <= 0.0:
        return 0.0
    hc = upper.collect("crustal_thickness_m")
    continental = effective_is_continental_from_codes(upper.collect("crust_type_code"), upper.crust_type == "continental")
    if not np.any(continental):
        return 0.0
    seed_mask = np.zeros(len(hc), dtype=bool)
    seed_mask[seed] = True
    hops = cratons._hops(UNDERTHRUST_REACH_KM, line_spacing_rad(world.node_density))
    reach = quad_tectonics.hop_distance(upper, seed_mask, hops, adjacency) <= hops
    hc_before = hc.copy()
    remaining, changed = quad_tectonics._fill_capped_volume(
        hc, upper.node_areas_m2(), reach & continental, volume, SUTURE_ACCRETION_MAX_HC_M
    )
    placed = volume - remaining
    if placed <= 0.0:
        return 0.0
    quad_tectonics.settle_received_crust(
        upper, hc_before, hc, hc - hc_before, changed, volume, material_volume, restite_volume,
        upper.collect("continental_material_m"), upper.collect("restite_m"),
        craton_volume=craton_volume, craton_formed_years=craton_formed_years,
    )
    return placed


def _book(
    world: "World", column: ConsumedColumn, upper: "PlateWithSparseQuadPatch", result: TransferResult, thrust_up: float
) -> None:
    """Every share in its own account. Material and restite on the receivers were moved by
    the placement; only what leaves the surface is booked here. The mobile cover rides at the
    top of the thrust-up crust, so it shares that crust's fate in proportion: metamorphosed
    into the crust it lands as, or subducted with what had no outlet. A terrane's returned
    share stays on its cells (`_keep_remainder`), unbooked: the carrier's own-plate path
    books it as it places it, so it counts as donated there rather than here."""
    hc = max(column.hc_m3, 0.0)
    lost_share = result.lost_m3 / hc if hc > 0.0 else 0.0
    unplaced_share = result.unplaced_m3 / hc if hc > 0.0 else 0.0
    returned_share = result.returned_m3 / hc if hc > 0.0 else 0.0
    kept_share = max(1.0 - lost_share - unplaced_share - returned_share, 0.0)
    continental_ledger.record(world, "collision_lower_crust_subducted_m3", lost_share * column.material_m3)
    continental_ledger.record(world, "collision_subducted_m3", unplaced_share * column.material_m3)
    if not column.terrane:
        # A terrane's kept craton landed as craton (`transfer_column`).
        cratons.record(world, "collision_reworked_m3", kept_share * column.craton_m3)
    cratons.record(world, "subducted_m3", (lost_share + unplaced_share) * column.craton_m3)
    cover = max(column.mobile_cover_m3, 0.0)
    cover_subducted = cover * min(result.unplaced_m3 / thrust_up, 1.0) if thrust_up > 0.0 else 0.0
    cover_returned = _returned_cover_m3(column, result)
    mobile_cover.record(world, "accreted_m3", max(cover - cover_subducted - cover_returned, 0.0))
    mobile_cover.record(world, "subducted_m3", cover_subducted)
    for account, volume in (
        ("suture_donated_m3", hc - result.returned_m3),
        ("terrane_returned_m3", result.returned_m3),
        ("suture_scraped_m3", result.scraped_m3),
        ("suture_underthrust_m3", result.underthrust_m3),
        ("suture_underthrust_placed_m3", result.underthrust_placed_m3),
        ("suture_lower_crust_subducted_m3", result.lost_m3),
        ("upper_plate_placed_m3", result.upper_placed_m3),
        ("overrider_placed_m3", result.handed_m3),
        ("no_outlet_subducted_m3", result.unplaced_m3),
        ("restite_subducted_m3", (lost_share + unplaced_share) * column.restite_m3),
    ):
        orogeny.record(world, account, volume)
    subducted_hm = column.hm_m3 * (1.0 - returned_share)
    if column.terrane:
        # The carrier's slab: a terrane is not a continental suture front.
        hm_ledger.record_typed_sink(world, "oceanic_and_deep_subduction", subducted_hm, continental=True)
        if world.debug_diagnostics:
            counts = _terrane_stats(world)
            counts["partial_fronts" if result.returned_m3 > 0.0 else "docked_fronts"] += 1
            counts["docked_hc_m3"] += hc - result.returned_m3
            counts["returned_hc_m3"] += result.returned_m3
            counts["docked_craton_m3"] += (1.0 - returned_share) * column.craton_m3
        return
    hm_ledger.record_suture_front(
        world,
        column.plate_id,
        [upper.plate_id],
        subducted_hm,
        0.0,
        donor_is_continental=True,
        placed_hm_continental_m3=0.0,
        subducted_hm_m3=subducted_hm,
    )


_TERRANE_COUNTERS = (
    "polarized_cells", "contact_docking_cells", "terrane_to_terrane_cells", "docked_fronts", "partial_fronts", "fallback_fronts",
    "relocated_fronts", "docked_hc_m3", "docked_craton_m3", "returned_hc_m3", "fallback_hc_m3", "relocated_hc_m3",
)


def _terrane_stats(world: "World") -> dict:
    """Terrane docking counters (issue #321), inside `World.suture_transfer_stats`."""
    stats = world.suture_transfer_stats
    if not isinstance(stats, dict):
        stats = world.suture_transfer_stats = {}
    counts = stats.setdefault("terranes", {})
    for key in _TERRANE_COUNTERS:
        counts.setdefault(key, 0)
    return counts


def record_terrane_fallback(world: "World | None", hc_m3: float, relocated: bool, returned: bool = False) -> None:
    """Count one terrane front left to its carrier's own-plate path, and whether it had to
    relocate onto new footprint. A partial dock's `returned` remainder was counted as a
    partial front already; only its relocation is counted here."""
    if world is None or not world.debug_diagnostics:
        return
    counts = _terrane_stats(world)
    if not returned:
        counts["fallback_fronts"] += 1
        counts["fallback_hc_m3"] += hc_m3
    if relocated:
        counts["relocated_fronts"] += 1
        counts["relocated_hc_m3"] += hc_m3


def _record_stats(world: "World", column: ConsumedColumn, result: TransferResult, hc: float, seconds: float) -> None:
    """Cumulative counters plus per-step cost, kept in constant size: the current step's
    fronts and seconds roll into a count of steps with fronts and the slowest step."""
    stats = world.suture_transfer_stats
    if not isinstance(stats, dict):
        stats = world.suture_transfer_stats = {}
    stats.pop("by_step", None)  # unbounded per-step history from early builds of this change
    stats["fronts"] = stats.get("fronts", 0) + 1
    stats["terrane_fronts"] = stats.get("terrane_fronts", 0) + int(column.terrane)
    stats["seconds"] = stats.get("seconds", 0.0) + seconds
    stats["max_front_seconds"] = max(stats.get("max_front_seconds", 0.0), seconds)
    stats["hc_m3"] = stats.get("hc_m3", 0.0) + hc
    for name, value in vars(result).items():
        stats[name] = stats.get(name, 0.0) + value
    if stats.get("current_step") != world.steps_taken:
        stats["current_step"] = world.steps_taken
        stats["current_step_seconds"] = 0.0
        stats["steps_with_fronts"] = stats.get("steps_with_fronts", 0) + 1
    stats["current_step_seconds"] += seconds
    stats["max_step_seconds"] = max(stats.get("max_step_seconds", 0.0), stats["current_step_seconds"])


def summary(world: "World") -> dict:
    """The configured partition (or why it's invalid) and, with debug diagnostics, the
    transfer counters. Never raises, so a bad knob can't break the stats endpoint."""
    try:
        shares = partition(world)
        report = {
            "scrape_fraction": shares.scrape,
            "underthrust_fraction": shares.underthrust,
            "lower_crust_loss_fraction": shares.loss,
        }
    except ValueError as error:
        report = {"error": str(error)}
    report["underthrust_reach_km"] = UNDERTHRUST_REACH_KM
    stats = dict(getattr(world, "suture_transfer_stats", None) or {})
    if isinstance(stats.get("terranes"), dict):
        stats["terranes"] = dict(stats["terranes"])
    return {"partition": report, "stats": stats}

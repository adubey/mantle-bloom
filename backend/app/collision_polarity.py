"""Which side of each continental collision front is the lower plate (issue #318, part of #315).

A continental collision inherits its polarity from the ocean that closed: the plate whose
ocean floor was subducting at that margin keeps its slab and goes under; the plate with the
arc stays on top. This module records that history where it happens and freezes a decision
for every active front before any plate deforms in a step. It decides; it doesn't act on the
decision. Retreat and Hm placement still treat both sides alike (#315 is the follow-up).

**Evidence.** Each step, three kinds of observations are binned onto the plates where they
happen, in each plate's own local frame, so they ride with the plate:

- `consumption`: oceanic-coded cells a plate actually lost to retreat under a neighbour
  (`record_consumption`, called from `quad_tectonics._retreat`). The plate that lost them
  gets `lower`, the neighbour that overrode them gets `upper` at the same place.
- `slab`: an oceanic plate's `torque.subducting_boundary_mask` nodes, which feel slab pull
  toward a converging neighbour. `lower` on the oceanic plate, `upper` on the neighbour's
  nearest nodes.
- `arc`: a continental plate's arc band (`lithosphere_plate.continental_arc_band`). `upper`
  on the arc's plate, `lower` on the oceanic neighbour. This is a model cue, not an
  observation: the arc band follows continental ownership, an oceanic neighbour and positive
  closing, with no test of which way a slab dips.

Each bin keeps its source, role, the neighbour at observation, when it was first and last
seen, and a count of node observations. Bins not seen for `EVIDENCE_LOOKBACK_YEARS` expire.
Because evidence sits on both sides of the contact, a third, oceanic plate that closed
between two continents can vanish and its history still survives on the continent it went
under.

**Fronts.** A front is an edge-connected stretch of contact between two plates where both
sides are continental-coded and the boundary is converging or overlapping. Fronts are matched
to the previous step's records by how much their plate-local extents overlap, one to one. A
record that matches several fronts has split (each extra front starts a child record with
the parent's polarity). Several records matching one front have grown together: those that
agree merge into the oldest, and those that disagree stay separate records, each keeping the
stretch of the front nearer its own nodes. A record keeps its polarity for life, so no
stretch of contact changes polarity while it lasts. A record expires after
`FRONT_EXPIRY_YEARS` without contact, which covers both separation and divergence, since a
diverging stretch is no longer a front; expired records are dropped before matching.

**Decision.** A new front looks up the evidence near it on both plates. The physical sources
(`consumption`, `slab`) decide when they agree; the `arc` cue decides only when no physical
evidence is present, and never overrules it. Evidence pointing both ways, or arcs on both
sides, is ambiguous. With nothing near the front, an adjacent front of the same pair lends
its polarity, then the pair's own evidence from elsewhere along their boundary decides, then
any front of the pair; only then a labelled heuristic (`_fallback`). See `_decide`.

**Freezing.** `observe_contacts` runs after every plate has shifted and before the first
`deform()`. It visits plates in id order, so its records and its per-plate lower/upper masks
(`World.collision_polarity_frame`) don't depend on `world.plates` order or deform order.

**Cost.** The prepass's boundary searches are the ones `deform()` needs, so they go through
`torque.BoundarySearchCache` and a plate that hasn't deformed yet is answered from the cache.
"""

from __future__ import annotations

import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import geometry, lithosphere, torque
from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from .lithosphere_plate import CONTINENTAL_CONTESTED_RETREAT_MIN_RUN, continental_arc_band

if TYPE_CHECKING:
    from .world import World

SOURCE_CONSUMPTION = 0
SOURCE_SLAB = 1
SOURCE_ARC = 2
SOURCE_NAMES = ("consumption", "slab", "arc")
_PHYSICAL_SOURCES = (SOURCE_CONSUMPTION, SOURCE_SLAB)
ROLE_UPPER = 1
ROLE_LOWER = -1

# Evidence bins are this many line spacings across, in the plate's own local frame. Two
# spacings keeps a long subduction zone to a few hundred bins while staying well inside the
# lookup radius below.
EVIDENCE_BIN_SPACINGS = 2.0
# A bin not observed for this long is forgotten. Long enough to bridge the last stage of an
# ocean's closure, when the arc and slab shut down before continents touch; short enough that
# a margin's polarity comes from its recent history, not from an ocean that closed long ago.
EVIDENCE_LOOKBACK_YEARS = 10_000_000.0
# Evidence counts toward a front when it lies within this many spacings of the front's nodes
# on the same plate. The arc band reaches `BOUNDARY_FORCE_REACH_MULTIPLIER` (3) spacings
# inboard of the contact, so 4 spacings covers it plus a bin's half-width.
EVIDENCE_LOOKUP_SPACINGS = 4.0
# Upper bound on evidence bins kept per plate, oldest dropped first. Ordinary worlds stay far
# below it; it only bounds a pathological case.
MAX_EVIDENCE_BINS_PER_PLATE = 20_000

# Contact nodes this close are one front: the boundary-force reach, so a stretch where the
# contact briefly stops converging for a node or two doesn't cut one front in two.
FRONT_LINK_SPACINGS = 3.0
# A front's nodes match a prior record's nodes within this distance on the same plate. A
# boundary moves at most a cell or two per step, so this follows a migrating front while
# keeping two fronts a few hundred km apart distinct.
FRONT_MATCH_SPACINGS = 3.0
# Smaller contacts are envelope fuzz, the same threshold continental retreat uses.
FRONT_MIN_NODES = CONTINENTAL_CONTESTED_RETREAT_MIN_RUN
# A new front with no evidence of its own copies a same-pair front this close (see `_decide`):
# twice the evidence lookup, comfortably past the matching tolerance.
FRONT_INHERIT_SPACINGS = 8.0
# A record with no matching front for this long expires: contact was lost, or the boundary
# stopped converging.
FRONT_EXPIRY_YEARS = 5_000_000.0

# Evidence one way is ambiguous when the opposing share exceeds this.
AMBIGUOUS_MINORITY_SHARE = 0.25
# The fallback's motion test calls two speeds tied within this, m/s (0.1 cm/yr).
FALLBACK_SPEED_TIE_M_PER_S = 0.001 / torque.SECONDS_PER_YEAR
# Topology reconciliation: a stored point stays with a descendant plate whose nearest node is
# within this many spacings of it.
REHOME_SPACINGS = 2.0


@dataclass
class EvidenceStore:
    """One plate's binned evidence, one row per (bin, source, role, neighbour)."""

    bin_rad: float
    keys: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    source: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int8))
    role: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int8))
    neighbour: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    first_seen: np.ndarray = field(default_factory=lambda: np.zeros(0))
    last_seen: np.ndarray = field(default_factory=lambda: np.zeros(0))
    count: np.ndarray = field(default_factory=lambda: np.zeros(0))

    def __len__(self) -> int:
        return len(self.keys)

    def local_points(self) -> np.ndarray:
        return _bin_centres(self.keys, self.bin_rad)


@dataclass
class CollisionFront:
    """A persistent record of one continental contact front and its frozen polarity."""

    front_id: int
    plate_ids: tuple[int, int]  # ascending
    lower_plate_id: int
    # "consumption", "slab" or "arc" when evidence decided; "inherited" when copied from
    # another front of the same pair; "fallback" otherwise.
    source: str
    # Evidence was present but contradictory (both physical directions, or arcs on both
    # sides), so the fallback decided.
    ambiguous: bool
    # "motion", "size" or "plate_id" for a fallback decision: which tier broke the tie.
    fallback_basis: str | None
    # Signed evidence weight per source at decision time, positive toward plate_ids[0] being
    # the upper plate.
    votes: dict[str, float]
    established_years: float
    last_contact_years: float
    contact_steps: int = 1
    # The record this one took its polarity from: the record it split off, or the same-pair
    # front it inherited from.
    parent_id: int | None = None
    # Where the deciding evidence came from: "front" (near this front), "pair" (elsewhere along
    # the pair's boundary), "record" (inherited), "fallback" -- see `_decide`.
    scope: str = "front"
    # plate id -> (k, 3) every node of this front on that plate at its last contact, in that
    # plate's local frame. All of them, not a sample: a topology split must see every
    # stretch of contact each descendant holds (`end_topology`).
    side_points: dict[int, np.ndarray] = field(default_factory=dict)

    @property
    def upper_plate_id(self) -> int:
        a, b = self.plate_ids
        return b if self.lower_plate_id == a else a


@dataclass
class PlateCollisionMasks:
    """One plate's frozen front membership for this step, over its node order at the prepass
    (which is its node order at the start of its `deform()`)."""

    lower: np.ndarray  # front nodes where this plate is the lower plate
    upper: np.ndarray  # front nodes where this plate is the upper plate
    front_id: np.ndarray  # record id per front node, -1 elsewhere
    retreat_eligible: np.ndarray  # lower & contested: what may be consumed
    override: np.ndarray  # upper & contested: what overrides rather than retreats


@dataclass
class PolarityFrame:
    """The prepass's frozen decisions for one step."""

    elapsed_years: float
    masks: dict[int, PlateCollisionMasks]
    # front id -> (lower plate id, upper plate id), for fronts in contact this step.
    polarity: dict[int, tuple[int, int]]


@dataclass
class _Contact:
    plate: object
    points: np.ndarray
    continental: np.ndarray
    inputs: torque.BoundaryForceInputs
    convergent: np.ndarray
    contested: np.ndarray


# --- Evidence bins ---------------------------------------------------------------------------


def _bin_half_range(bin_rad: float) -> int:
    return int(np.ceil(1.0 / bin_rad)) + 1


def _bin_keys(local_xyz: np.ndarray, bin_rad: float) -> np.ndarray:
    half = _bin_half_range(bin_rad)
    side = 2 * half + 1
    cells = np.floor(np.asarray(local_xyz, dtype=float) / bin_rad).astype(np.int64) + half
    return (cells[:, 0] * side + cells[:, 1]) * side + cells[:, 2]


def _bin_centres(keys: np.ndarray, bin_rad: float) -> np.ndarray:
    half = _bin_half_range(bin_rad)
    side = 2 * half + 1
    keys = np.asarray(keys, dtype=np.int64)
    cells = np.stack([keys // (side * side), (keys // side) % side, keys % side], axis=-1)
    return geometry.normalize((cells - half + 0.5) * bin_rad)


def _merge_rows(store: EvidenceStore, rows: dict[str, np.ndarray]) -> None:
    """Fold `rows` into `store`: one row per (bin, source, role, neighbour), keeping the
    earliest first sighting, the latest last sighting and the summed count. Sorted, so the
    result doesn't depend on the order observations arrived in."""
    keys = np.concatenate([store.keys, rows["keys"]])
    if not len(keys):
        return
    source = np.concatenate([store.source, rows["source"]])
    role = np.concatenate([store.role, rows["role"]])
    neighbour = np.concatenate([store.neighbour, rows["neighbour"]])
    first = np.concatenate([store.first_seen, rows["first_seen"]])
    last = np.concatenate([store.last_seen, rows["last_seen"]])
    count = np.concatenate([store.count, rows["count"]])
    order = np.lexsort((neighbour, role, source, keys))
    keys, source, role, neighbour = keys[order], source[order], role[order], neighbour[order]
    first, last, count = first[order], last[order], count[order]
    starts = np.flatnonzero(
        np.concatenate(
            [[True], (keys[1:] != keys[:-1]) | (source[1:] != source[:-1]) | (role[1:] != role[:-1]) | (neighbour[1:] != neighbour[:-1])]
        )
    )
    store.keys = keys[starts]
    store.source = source[starts]
    store.role = role[starts]
    store.neighbour = neighbour[starts]
    store.first_seen = np.minimum.reduceat(first, starts)
    store.last_seen = np.maximum.reduceat(last, starts)
    store.count = np.add.reduceat(count, starts)


def _keep_rows(store: EvidenceStore, keep: np.ndarray) -> None:
    store.keys = store.keys[keep]
    store.source = store.source[keep]
    store.role = store.role[keep]
    store.neighbour = store.neighbour[keep]
    store.first_seen = store.first_seen[keep]
    store.last_seen = store.last_seen[keep]
    store.count = store.count[keep]


def _store_for(world: "World", plate_id: int) -> EvidenceStore:
    store = world.collision_evidence.get(plate_id)
    if store is None:
        store = EvidenceStore(EVIDENCE_BIN_SPACINGS * line_spacing_rad(world.node_density))
        world.collision_evidence[plate_id] = store
    return store


def add_evidence(
    world: "World", plate, world_xyz: np.ndarray, source: int, role: int, neighbour_ids: np.ndarray, years: float | None = None
) -> None:
    """Bin observations at `world_xyz` (on `plate`) into its evidence store, with the
    neighbour responsible for each one."""
    world_xyz = np.asarray(world_xyz, dtype=float).reshape(-1, 3)
    if not len(world_xyz):
        return
    now = world.elapsed_years if years is None else years
    store = _store_for(world, plate.plate_id)
    n = len(world_xyz)
    _merge_rows(
        store,
        {
            "keys": _bin_keys(geometry.to_local(plate.frame, world_xyz), store.bin_rad),
            "source": np.full(n, source, dtype=np.int8),
            "role": np.full(n, role, dtype=np.int8),
            "neighbour": np.broadcast_to(np.asarray(neighbour_ids, dtype=np.int64), (n,)).copy(),
            "first_seen": np.full(n, now),
            "last_seen": np.full(n, now),
            "count": np.ones(n),
        },
    )
    _stats(world)[f"evidence_obs_{SOURCE_NAMES[source]}"] += n


def record_consumption(world: "World", plate, world_xyz: np.ndarray, is_oceanic: np.ndarray, neighbour_ids: np.ndarray, neighbours: list) -> None:
    """Oceanic-coded cells `plate` lost to retreat this step: `lower` evidence on `plate`,
    and `upper` on the neighbour that overrode each one, at the same place in its frame.
    `neighbour_ids` is each removed cell's nearest other plate, from the step's
    `BoundaryForceInputs`; `neighbours` are the plate objects those ids refer to."""
    mask = np.asarray(is_oceanic, dtype=bool) & (np.asarray(neighbour_ids) >= 0)
    if not np.any(mask):
        return
    points = np.asarray(world_xyz, dtype=float)[mask]
    ids = np.asarray(neighbour_ids, dtype=np.int64)[mask]
    add_evidence(world, plate, points, SOURCE_CONSUMPTION, ROLE_LOWER, ids)
    by_id = {p.plate_id: p for p in neighbours}
    for neighbour_id in np.unique(ids):
        neighbour = by_id.get(int(neighbour_id))
        if neighbour is not None:
            add_evidence(world, neighbour, points[ids == neighbour_id], SOURCE_CONSUMPTION, ROLE_UPPER, plate.plate_id)


def _expire_evidence(world: "World", now: float) -> None:
    for plate_id in list(world.collision_evidence):
        store = world.collision_evidence[plate_id]
        keep = store.last_seen >= now - EVIDENCE_LOOKBACK_YEARS
        if len(store) and np.count_nonzero(keep) > MAX_EVIDENCE_BINS_PER_PLATE:
            # Most recent first; ties by position so the cut is deterministic.
            order = np.lexsort((store.neighbour, store.role, store.source, store.keys, -store.last_seen))
            keep = np.zeros(len(store), dtype=bool)
            keep[order[:MAX_EVIDENCE_BINS_PER_PLATE]] = True
        _keep_rows(store, keep)
        if not len(store):
            del world.collision_evidence[plate_id]


# --- The prepass -----------------------------------------------------------------------------


_STAT_KEYS = (
    "prepass_calls",
    "prepass_seconds",
    "prepass_searches",
    "deform_searches",
    "deform_searches_reused",
    "front_steps",
    "fronts_created",
    "fronts_split",
    "fronts_merged",
    "conflicts_kept_apart",
    "fronts_expired",
    "fronts_dropped_topology",
    "fronts_split_topology",
    "fronts_fused",
    "decided_consumption",
    "decided_slab",
    "decided_arc",
    "decided_inherited",
    "decided_fallback",
    "scope_front",
    "scope_pair",
    "scope_record",
    "fallback_no_evidence",
    "fallback_ambiguous",
    "ambiguous_physical",
    "ambiguous_opposing_arcs",
    "arc_overruled",
    "fallback_motion",
    "fallback_size",
    "fallback_plate_id",
    "evidence_obs_consumption",
    "evidence_obs_slab",
    "evidence_obs_arc",
)


def _stats(world: "World") -> dict:
    stats = world.collision_polarity_stats
    for key in _STAT_KEYS:
        stats.setdefault(key, 0)
    return stats


def observe_contacts(world: "World", years: float) -> PolarityFrame:
    """The prepass: observe every plate's contacts, update evidence and front records, and
    freeze this step's polarity and masks onto `world.collision_polarity_frame`. Run after
    every plate has shifted and before any deforms."""
    start = time.perf_counter()
    stats = _stats(world)
    spacing_rad = line_spacing_rad(world.node_density)
    reach_rad = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing_rad
    cache = world.boundary_search_cache
    if cache is None:
        cache = torque.BoundarySearchCache()
    searches_before = cache.searches
    now = world.elapsed_years

    plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
    by_id = {p.plate_id: p for p in plates}
    contacts: dict[int, _Contact] = {}
    for plate in plates:
        neighbours = cache.neighbours(plate, plates, reach_rad)
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing_rad, reach_rad, cache)
        convergent, _, _, contested = torque.classify_boundary_nodes(plate, neighbours, inputs, reach_rad, cache)
        contacts[plate.plate_id] = _Contact(
            plate=plate,
            points=inputs.own_points,
            continental=effective_is_continental_from_codes(inputs.own_crust_type_codes, plate.crust_type == "continental"),
            inputs=inputs,
            convergent=convergent,
            contested=contested,
        )
        closing = torque.boundary_closing_rate_m_per_s(plate, inputs)
        _observe(world, by_id, contacts[plate.plate_id], continental_arc_band(plate, inputs, closing), SOURCE_ARC, ROLE_UPPER)
        _observe(world, by_id, contacts[plate.plate_id], torque.subducting_boundary_mask(plate, inputs, reach_rad), SOURCE_SLAB, ROLE_LOWER)
    _expire_evidence(world, now)

    fronts_by_pair = _detect_fronts(contacts, spacing_rad)
    masks = {
        pid: PlateCollisionMasks(
            lower=np.zeros(len(c.points), dtype=bool),
            upper=np.zeros(len(c.points), dtype=bool),
            front_id=np.full(len(c.points), -1, dtype=np.int64),
            retreat_eligible=np.zeros(len(c.points), dtype=bool),
            override=np.zeros(len(c.points), dtype=bool),
        )
        for pid, c in contacts.items()
    }
    # Expire first, so matching and every decision below only ever see live records.
    expired = [r for r in world.collision_fronts if r.last_contact_years < now - FRONT_EXPIRY_YEARS]
    if expired:
        stats["fronts_expired"] += len(expired)
        world.collision_fronts = [r for r in world.collision_fronts if r.last_contact_years >= now - FRONT_EXPIRY_YEARS]

    polarity: dict[int, tuple[int, int]] = {}
    for pair, fronts in fronts_by_pair.items():
        segments = _match_fronts(world, pair, fronts, contacts, spacing_rad, now)
        for segment, record in segments:
            polarity[record.front_id] = (record.lower_plate_id, record.upper_plate_id)
            for side, idx in segment.items():
                m = masks[side]
                if side == record.lower_plate_id:
                    m.lower[idx] = True
                else:
                    m.upper[idx] = True
                m.front_id[idx] = record.front_id
        stats["front_steps"] += len(segments)

    for pid, m in masks.items():
        contested = contacts[pid].contested
        m.retreat_eligible = m.lower & contested
        m.override = m.upper & contested
    frame = PolarityFrame(elapsed_years=now, masks=masks, polarity=polarity)
    world.collision_polarity_frame = frame
    stats["prepass_calls"] += 1
    stats["prepass_searches"] += cache.searches - searches_before
    stats["prepass_seconds"] += time.perf_counter() - start
    cache.prepass_mark = (cache.searches, cache.reused)
    return frame


def finish_deform_pass(world: "World") -> None:
    """Book how much of the deform pass the prepass's searches spared, and drop the step's
    search cache."""
    cache = world.boundary_search_cache
    if cache is not None:
        stats = _stats(world)
        searches, reused = cache.prepass_mark
        stats["deform_searches"] += cache.searches - searches
        stats["deform_searches_reused"] += cache.reused - reused
    world.boundary_search_cache = None


def _observe(world: "World", by_id: dict, c: _Contact, mask: np.ndarray, source: int, role: int) -> None:
    """Evidence on both sides: `role` on this plate's `mask` nodes, the opposite role on each
    responsible neighbour's nearest node. Neighbours are observed in whatever order plates
    are visited; `_merge_rows` sorts, so that order doesn't show in the store."""
    idx = np.flatnonzero(mask & (c.inputs.neighbor_plate_id >= 0))
    if not len(idx):
        return
    ids = c.inputs.neighbor_plate_id[idx]
    add_evidence(world, c.plate, c.points[idx], source, role, ids)
    for neighbour_id in np.unique(ids):
        neighbour = by_id.get(int(neighbour_id))
        if neighbour is None:
            continue
        sel = idx[ids == neighbour_id]
        points = neighbour.all_points_and_elevation()[0][c.inputs.neighbor_node_index[sel]]
        add_evidence(world, neighbour, points, source, -role, c.plate.plate_id)


def _detect_fronts(contacts: dict[int, _Contact], spacing_rad: float) -> dict[tuple[int, int], list[dict[int, np.ndarray]]]:
    """Every pair's connected continental contact fronts, each as {plate id: node indices}.
    Pairs and fronts come out in a fixed order: pairs ascending, fronts by first node."""
    sides: dict[tuple[int, int], dict[int, np.ndarray]] = defaultdict(dict)
    for pid, c in contacts.items():
        nb = c.inputs.neighbor_plate_id
        idx = np.flatnonzero(c.convergent & c.continental & (nb >= 0))
        for q in np.unique(nb[idx]):
            q = int(q)
            if q not in contacts:
                continue
            sel = idx[nb[idx] == q]
            sel = sel[contacts[q].continental[c.inputs.neighbor_node_index[sel]]]
            if len(sel):
                sides[(min(pid, q), max(pid, q))][pid] = sel

    link = FRONT_LINK_SPACINGS * spacing_rad
    out: dict[tuple[int, int], list[dict[int, np.ndarray]]] = {}
    for pair in sorted(sides):
        a, b = pair
        ia = sides[pair].get(a, np.zeros(0, dtype=np.int64))
        ib = sides[pair].get(b, np.zeros(0, dtype=np.int64))
        points = np.concatenate([contacts[a].points[ia], contacts[b].points[ib]])
        n = len(points)
        if n < FRONT_MIN_NODES:
            continue
        links = cKDTree(points).query_pairs(link, output_type="ndarray")
        graph = coo_matrix((np.ones(len(links)), (links[:, 0], links[:, 1])), shape=(n, n)) if len(links) else coo_matrix((n, n))
        _, labels = connected_components(graph, directed=False)
        _, first = np.unique(labels, return_index=True)
        fronts = []
        for label in labels[np.sort(first)]:
            members = np.flatnonzero(labels == label)
            if len(members) < FRONT_MIN_NODES:
                continue
            fronts.append({a: ia[members[members < len(ia)]], b: ib[members[members >= len(ia)] - len(ia)]})
        if fronts:
            out[pair] = fronts
    return out


def _overlap(front_points: np.ndarray, record_points: np.ndarray, tol: float) -> float:
    """Dice overlap of two point sets on one plate: the share of both sets lying within `tol`
    of the other."""
    if not len(front_points) or not len(record_points):
        return 0.0
    near_f = np.isfinite(cKDTree(record_points).query(front_points, distance_upper_bound=tol)[0])
    near_r = np.isfinite(cKDTree(front_points).query(record_points, distance_upper_bound=tol)[0])
    return float(near_f.sum() + near_r.sum()) / float(len(front_points) + len(record_points))


def _match_fronts(
    world: "World", pair: tuple[int, int], fronts: list[dict[int, np.ndarray]], contacts: dict[int, _Contact], spacing_rad: float, now: float
) -> list[tuple[dict[int, np.ndarray], CollisionFront]]:
    """This pair's fronts matched to records: (segment, record) pairs, where a segment is a
    front's nodes or, for a front two disagreeing records share, the part nearer one of them.

    Fronts take prior records one to one by overlap. A front no free record overlaps but a
    taken one does has split off it, and inherits its polarity. Several records overlapping
    one front have grown together: records that agree on polarity merge into the oldest, and
    records that disagree never merge -- each keeps the stretch nearer its own stored nodes,
    so no stretch of contact ever changes polarity while it lasts."""
    stats = _stats(world)
    tol = FRONT_MATCH_SPACINGS * spacing_rad
    records = [r for r in world.collision_fronts if r.plate_ids == pair]
    stored_world = [
        {side: geometry.to_world(contacts[side].plate.frame, r.side_points[side]) for side in pair if len(r.side_points.get(side, ()))}
        for r in records
    ]
    scores = np.zeros((len(fronts), len(records)))
    for j in range(len(records)):
        for side, points in stored_world[j].items():
            for i, front in enumerate(fronts):
                scores[i, j] += 0.5 * _overlap(contacts[side].points[front[side]], points, tol)

    assigned: dict[int, int] = {}
    taken: set[int] = set()
    candidates = sorted(
        ((-scores[i, j], records[j].front_id, i, j) for i in range(len(fronts)) for j in range(len(records)) if scores[i, j] > 0.0)
    )
    for _, _, i, j in candidates:
        if i in assigned or j in taken:
            continue
        assigned[i] = j
        taken.add(j)

    # Each matched front's group: its record plus every free record it also overlaps.
    claimed = set(taken)
    survivors: dict[int, list[int]] = {}
    absorbed: set[int] = set()
    for i in sorted(assigned):
        group = [assigned[i]] + [int(j) for j in np.flatnonzero(scores[i] > 0.0) if int(j) not in claimed]
        claimed.update(group)
        by_lower: dict[int, list[int]] = defaultdict(list)
        for j in group:
            by_lower[records[j].lower_plate_id].append(j)
        keep = []
        for lower in sorted(by_lower):
            members = sorted(by_lower[lower], key=lambda j: (records[j].established_years, records[j].front_id))
            keep.append(members[0])
            for j in members[1:]:
                absorbed.add(records[j].front_id)
                stats["fronts_merged"] += 1
        if len(keep) > 1:
            stats["conflicts_kept_apart"] += len(keep) - 1
        survivors[i] = keep
    if absorbed:
        world.collision_fronts = [r for r in world.collision_fronts if r.front_id not in absorbed]

    out: list[tuple[dict[int, np.ndarray], CollisionFront]] = []
    for i, front in enumerate(fronts):
        if i in survivors:
            keep = survivors[i]
            if len(keep) == 1:
                segments = [(front, records[keep[0]])]
            else:
                segments = _segment(front, [stored_world[j] for j in keep], contacts)
                segments = [(segment, records[keep[k]]) for k, segment in enumerate(segments) if segment is not None]
            for _, record in segments:
                record.contact_steps += 1
        else:
            overlapping = np.flatnonzero(scores[i] > 0.0)
            if len(overlapping):
                # Split: a front off a record already matched elsewhere inherits its polarity.
                parent = records[int(max(overlapping, key=lambda j: (scores[i, j], -records[j].front_id)))]
                record = _new_record(world, pair, now)
                record.lower_plate_id = parent.lower_plate_id
                record.source = parent.source
                record.ambiguous = parent.ambiguous
                record.fallback_basis = parent.fallback_basis
                record.votes = dict(parent.votes)
                record.established_years = parent.established_years
                record.parent_id = parent.front_id
                record.scope = parent.scope
                stats["fronts_split"] += 1
            else:
                record = _new_record(world, pair, now)
                _decide(world, record, front, contacts, spacing_rad)
                stats["fronts_created"] += 1
            world.collision_fronts.append(record)
            segments = [(front, record)]
        for segment, record in segments:
            record.last_contact_years = now
            for side in pair:
                if len(segment[side]):
                    record.side_points[side] = geometry.to_local(contacts[side].plate.frame, contacts[side].points[segment[side]])
            out.append((segment, record))
    return out


def _segment(front: dict[int, np.ndarray], stored: list[dict[int, np.ndarray]], contacts: dict[int, _Contact]) -> list[dict[int, np.ndarray] | None]:
    """Split one front's nodes between records that disagree: each node goes to the record
    whose stored nodes on the same plate are nearest (ties to the first). None for a record
    left with no nodes this step."""
    owner: dict[int, np.ndarray] = {}
    for side, idx in front.items():
        points = contacts[side].points[idx]
        dist = np.full((len(stored), len(idx)), np.inf)
        for k, record_points in enumerate(stored):
            if side in record_points and len(idx):
                dist[k] = cKDTree(record_points[side]).query(points)[0]
        owner[side] = np.argmin(dist, axis=0) if len(idx) else np.zeros(0, dtype=int)
    out = []
    for k in range(len(stored)):
        segment = {side: idx[owner[side] == k] for side, idx in front.items()}
        out.append(segment if any(len(v) for v in segment.values()) else None)
    return out


def _new_record(world: "World", pair: tuple[int, int], now: float) -> CollisionFront:
    record = CollisionFront(
        front_id=world.next_collision_front_id,
        plate_ids=pair,
        lower_plate_id=pair[0],
        source="fallback",
        ambiguous=False,
        fallback_basis=None,
        votes={},
        established_years=now,
        last_contact_years=now,
    )
    world.next_collision_front_id += 1
    return record


def _verdict(positive: float, negative: float) -> str | None:
    """"a" (positive wins), "b", "ambiguous", or None for no evidence."""
    total = positive + negative
    if total <= 0.0:
        return None
    if min(positive, negative) / total > AMBIGUOUS_MINORITY_SHARE:
        return "ambiguous"
    return "a" if positive > negative else "b"


def _tally(world: "World", pair: tuple[int, int], select) -> tuple[np.ndarray, np.ndarray]:
    """Evidence weight per source on both plates of `pair`, over the rows `select(side, store)`
    picks: (toward `pair[0]` upper, toward `pair[1]` upper)."""
    a, b = pair
    positive = np.zeros(len(SOURCE_NAMES))
    negative = np.zeros(len(SOURCE_NAMES))
    for side in pair:
        store = world.collision_evidence.get(side)
        if store is None or not len(store):
            continue
        rows = select(side, store)
        if not np.any(rows):
            continue
        # Upper on `a` (or lower on `b`) says `a` is upper.
        toward_a = (store.role[rows] == ROLE_UPPER) == (side == a)
        np.add.at(positive, store.source[rows][toward_a], store.count[rows][toward_a])
        np.add.at(negative, store.source[rows][~toward_a], store.count[rows][~toward_a])
    return positive, negative


def _judge(world: "World", record: CollisionFront, positive: np.ndarray, negative: np.ndarray) -> str | None:
    """Apply one tier's evidence to `record`: physical sources decide when they agree, the arc
    cue only without them. Returns "decided", "ambiguous" or None (no evidence)."""
    stats = _stats(world)
    a, b = record.plate_ids
    physical = list(_PHYSICAL_SOURCES)
    verdict = _verdict(float(positive[physical].sum()), float(negative[physical].sum()))
    arc_verdict = _verdict(float(positive[SOURCE_ARC]), float(negative[SOURCE_ARC]))
    if verdict in ("a", "b"):
        winning = positive if verdict == "a" else negative
        record.source = SOURCE_NAMES[physical[int(np.argmax(winning[physical]))]]
        record.lower_plate_id = b if verdict == "a" else a
        if arc_verdict in ("a", "b") and arc_verdict != verdict:
            stats["arc_overruled"] += 1
        return "decided"
    if verdict == "ambiguous":
        stats["ambiguous_physical"] += 1
        return "ambiguous"
    if arc_verdict in ("a", "b"):
        record.source = "arc"
        record.lower_plate_id = b if arc_verdict == "a" else a
        return "decided"
    if arc_verdict == "ambiguous":
        stats["ambiguous_opposing_arcs"] += 1
        return "ambiguous"
    return None


def _decide(world: "World", record: CollisionFront, front: dict[int, np.ndarray], contacts: dict[int, _Contact], spacing_rad: float) -> None:
    """Set a new record's polarity, from the first of these that can decide:

    1. "front": evidence on either plate within `EVIDENCE_LOOKUP_SPACINGS` of the front.
    2. "record": a live front of the same pair within `FRONT_INHERIT_SPACINGS`, whose
       polarity this one copies. A fragment beside an established front is part of the same
       collision; deciding it afresh lets it disagree, and when the two grow together one of
       them would have to give way.
    3. "pair": evidence anywhere on either plate that names the other as the neighbour --
       still this pair's own recent history, from elsewhere along their shared boundary. A
       collision often starts as a few continental cells touching along a margin whose
       subduction was recorded a few hundred km away.
    4. "record": the nearest live front of the same pair at any distance.
    5. "fallback": `_fallback`'s heuristic.

    Contradictory evidence at the front stops the search there: it is ambiguous, and evidence
    from farther away can't be trusted over it, so the fallback decides."""
    stats = _stats(world)
    a, b = record.plate_ids
    radius = EVIDENCE_LOOKUP_SPACINGS * spacing_rad

    def near_front(side: int, store: EvidenceStore) -> np.ndarray:
        other = b if side == a else a
        look = contacts[side].points[front[side]] if len(front[side]) else contacts[other].points[front[other]]
        rows = geometry.to_world(contacts[side].plate.frame, store.local_points())
        return np.isfinite(cKDTree(look).query(rows, distance_upper_bound=radius)[0])

    def naming_other(side: int, store: EvidenceStore) -> np.ndarray:
        return store.neighbour == (b if side == a else a)

    sibling, gap = _nearest_pair_record(world, record, front, contacts)

    def inherit() -> bool:
        if sibling is None:
            return False
        record.lower_plate_id = sibling.lower_plate_id
        record.source = "inherited"
        record.scope = "record"
        record.parent_id = sibling.front_id
        stats["decided_inherited"] += 1
        stats["scope_record"] += 1
        return True

    tiers = (("front", near_front), ("adjacent", None), ("pair", naming_other), ("any_record", None))
    for scope, select in tiers:
        if select is None:
            # A front of the same pair close by is the same collision: copy it. Farther away,
            # the pair's own evidence comes first.
            if (scope == "any_record" or gap <= FRONT_INHERIT_SPACINGS * spacing_rad) and inherit():
                return
            continue
        positive, negative = _tally(world, record.plate_ids, select)
        outcome = _judge(world, record, positive, negative)
        if outcome is None:
            continue
        record.votes = {name: float(positive[k] - negative[k]) for k, name in enumerate(SOURCE_NAMES)}
        record.scope = scope
        if outcome == "decided":
            stats[f"decided_{record.source}"] += 1
            stats[f"scope_{scope}"] += 1
            return
        record.ambiguous = True
        break

    record.source = "fallback"
    record.scope = "fallback"
    record.lower_plate_id, record.fallback_basis = _fallback(record.plate_ids, front, contacts, spacing_rad)
    stats["decided_fallback"] += 1
    stats["fallback_ambiguous" if record.ambiguous else "fallback_no_evidence"] += 1
    stats[f"fallback_{record.fallback_basis}"] += 1


def _nearest_pair_record(
    world: "World", record: CollisionFront, front: dict[int, np.ndarray], contacts: dict[int, _Contact]
) -> tuple[CollisionFront | None, float]:
    """The live record of the same pair whose stored nodes come closest to this front, ties to
    the lower id, and that distance; (None, inf) if the pair has no other record."""
    best: tuple[float, int] | None = None
    chosen = None
    for other in world.collision_fronts:
        if other.plate_ids != record.plate_ids or other.front_id == record.front_id:
            continue
        gap = np.inf
        for side in record.plate_ids:
            stored = other.side_points.get(side)
            if stored is None or not len(stored) or not len(front[side]):
                continue
            stored_world = geometry.to_world(contacts[side].plate.frame, stored)
            gap = min(gap, float(cKDTree(stored_world).query(contacts[side].points[front[side]])[0].min()))
        key = (gap, other.front_id)
        if best is None or key < best:
            best, chosen = key, other
    return chosen, (best[0] if best is not None else np.inf)


def _fallback(pair: tuple[int, int], front: dict[int, np.ndarray], contacts: dict[int, _Contact], spacing_rad: float) -> tuple[int, str]:
    """Heuristic polarity when the evidence can't decide. Returns (lower plate id, basis).

    1. "motion": the plate moving faster into the front, in the mantle frame (`omega` is
       absolute), goes down -- it's the one being driven into the boundary. Absolute motion
       depends on the reference frame and on trench motion, so this is a guess, not a test.
    2. "size": within `FALLBACK_SPEED_TIE_M_PER_S`, the smaller plate by area goes down.
    3. "plate_id": equal areas, the lower id goes down, so the result is still deterministic.

    Column buoyancy is deliberately not used: without a thermal density term every column in
    this model is buoyant (#315), so it can't say which side sinks."""
    a, b = pair
    speed = {}
    for side in pair:
        other = b if side == a else a
        if len(front[side]):
            points = contacts[side].points[front[side]]
            toward = contacts[side].inputs.direction_to_neighbor[front[side]]
        else:
            points = contacts[other].points[front[other]]
            toward = -contacts[other].inputs.direction_to_neighbor[front[other]]
        velocity = np.cross(np.asarray(contacts[side].plate.omega, dtype=float), points)
        speed[side] = float(np.mean(np.sum(velocity * toward, axis=-1))) * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    if abs(speed[a] - speed[b]) > FALLBACK_SPEED_TIE_M_PER_S:
        return (a if speed[a] > speed[b] else b), "motion"
    area = {side: float(np.sum(contacts[side].plate.accounting_areas_m2(spacing_rad))) for side in pair}
    if not np.isclose(area[a], area[b], rtol=1e-9, atol=0.0):
        return (a if area[a] < area[b] else b), "size"
    return a, "plate_id"


# --- Topology --------------------------------------------------------------------------------


def begin_topology(world: "World") -> dict[int, np.ndarray]:
    """Before `merge_split.apply_topology_changes`: start recording plate lineage and keep
    each plate's frame, so stored local points can still be placed in the world afterwards."""
    world.topology_lineage = []
    return {p.plate_id: p.frame.copy() for p in world.plates}


def end_topology(world: "World", frames_before: dict[int, np.ndarray]) -> None:
    """Move evidence and front records onto the plates that now carry them. A plate that
    merged into another, split, or fragmented hands each stored point to whichever of its
    descendants now holds that ground; a plate that vanished with no descendant takes its
    evidence and fronts with it. A front split between descendants keeps a record for every
    stretch of contact that survives, all with the parent's frozen decision. A stretch whose
    two sides end up on one plate is dropped: the suture is now internal."""
    lineage = world.topology_lineage or []
    world.topology_lineage = None
    live = {p.plate_id: p for p in world.plates}
    children: dict[int, list[int]] = defaultdict(list)
    for parent, child in lineage:
        children[parent].append(child)

    def changed(pid: int) -> bool:
        plate = live.get(pid)
        return plate is None or pid in children or not np.allclose(plate.frame, frames_before.get(pid, plate.frame))

    def descendants(pid: int) -> list[int]:
        seen, stack, out = {pid}, [pid], []
        if pid in live:
            out.append(pid)
        while stack:
            for child in children.get(stack.pop(), []):
                if child not in seen:
                    seen.add(child)
                    stack.append(child)
                    if child in live:
                        out.append(child)
        return sorted(set(out))

    tol = REHOME_SPACINGS * line_spacing_rad(world.node_density)

    def owners(pid: int, local: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """World positions of `local` (in `pid`'s pre-topology frame) and the descendant
        plate now holding each, -1 where none does."""
        frame = frames_before.get(pid)
        if frame is None or not len(local):
            return np.zeros((0, 3)), np.zeros(0, dtype=np.int64)
        world_xyz = geometry.to_world(frame, local)
        best = np.full(len(world_xyz), np.inf)
        owner = np.full(len(world_xyz), -1, dtype=np.int64)
        for candidate in descendants(pid):
            tree = live[candidate].get_node_kdtree()
            if tree is None:
                continue
            dist, _ = tree.query(world_xyz, distance_upper_bound=tol)
            closer = dist < best
            best[closer] = dist[closer]
            owner[closer] = candidate
        return world_xyz, owner

    stats = _stats(world)
    for pid in sorted(world.collision_evidence):
        if not changed(pid):
            continue
        store = world.collision_evidence.pop(pid)
        world_xyz, owner = owners(pid, store.local_points())
        for new_pid in np.unique(owner[owner >= 0]):
            sel = owner == new_pid
            target = _store_for(world, int(new_pid))
            _merge_rows(
                target,
                {
                    "keys": _bin_keys(geometry.to_local(live[int(new_pid)].frame, world_xyz[sel]), target.bin_rad),
                    "source": store.source[sel],
                    "role": store.role[sel],
                    "neighbour": store.neighbour[sel],
                    "first_seen": store.first_seen[sel],
                    "last_seen": store.last_seen[sel],
                    "count": store.count[sel],
                },
            )

    kept: list[CollisionFront] = []
    for record in world.collision_fronts:
        if not any(changed(side) for side in record.plate_ids):
            kept.append(record)
            continue
        # Each side's parts: (plate now holding them, their world positions), one per
        # descendant holding at least FRONT_MIN_NODES of the front's nodes.
        parts: dict[int, list[tuple[int, np.ndarray]]] = {}
        for side in record.plate_ids:
            stored = record.side_points.get(side, np.zeros((0, 3)))
            if not changed(side):
                parts[side] = [(side, geometry.to_world(live[side].frame, stored))] if len(stored) else []
                continue
            world_xyz, owner = owners(side, stored)
            parts[side] = [
                (int(pid), world_xyz[owner == pid])
                for pid in np.unique(owner[owner >= 0])
                if np.count_nonzero(owner == pid) >= FRONT_MIN_NODES
            ]
        a, b = record.plate_ids
        if not parts[a] or not parts[b]:
            stats["fronts_dropped_topology"] += 1
            continue
        # Pair each part with the nearest part across the contact, from both sides, so a split
        # on either side (or both) leaves one record per surviving stretch of contact.
        combos: set[tuple[int, int]] = set()
        for k, (_, points) in enumerate(parts[a]):
            combos.add((k, _nearest_part(points, parts[b])))
        for k, (_, points) in enumerate(parts[b]):
            combos.add((_nearest_part(points, parts[a]), k))
        pieces = []
        for ka, kb in sorted(combos):
            (pid_a, points_a), (pid_b, points_b) = parts[a][ka], parts[b][kb]
            if pid_a == pid_b:
                stats["fronts_fused"] += 1
                continue
            pieces.append((len(points_a) + len(points_b), pid_a, pid_b, points_a, points_b))
        if not pieces:
            continue
        # The largest piece keeps the record; every other surviving stretch is a child with the
        # same frozen decision.
        pieces.sort(key=lambda piece: (-piece[0], piece[1], piece[2]))
        for n, (_, pid_a, pid_b, points_a, points_b) in enumerate(pieces):
            target = record if n == 0 else _new_record(world, record.plate_ids, record.established_years)
            if n > 0:
                target.source = record.source
                target.ambiguous = record.ambiguous
                target.fallback_basis = record.fallback_basis
                target.votes = dict(record.votes)
                target.scope = record.scope
                target.last_contact_years = record.last_contact_years
                target.contact_steps = record.contact_steps
                target.parent_id = record.front_id
                stats["fronts_split_topology"] += 1
            target.lower_plate_id = pid_a if record.lower_plate_id == a else pid_b
            target.plate_ids = (min(pid_a, pid_b), max(pid_a, pid_b))
            target.side_points = {
                pid_a: geometry.to_local(live[pid_a].frame, points_a),
                pid_b: geometry.to_local(live[pid_b].frame, points_b),
            }
            kept.append(target)
    world.collision_fronts = kept


def _nearest_part(points: np.ndarray, parts: list[tuple[int, np.ndarray]]) -> int:
    """Index of the part in `parts` with a point nearest any of `points`, ties to the first."""
    gaps = [float(cKDTree(other).query(points)[0].min()) for _, other in parts]
    return int(np.argmin(gaps))


def note_lineage(world: "World", parent_id: int, child_id: int) -> None:
    """Record that `child_id` now carries ground that belonged to `parent_id` (a merge, a
    split, a fragment), for `end_topology`. A no-op outside a topology pass."""
    lineage = getattr(world, "topology_lineage", None)
    if lineage is not None:
        lineage.append((int(parent_id), int(child_id)))


# --- Reporting -------------------------------------------------------------------------------


def summary(world: "World") -> dict:
    """The replay diagnostics: decision-source shares, fallback and ambiguity rates, record
    churn, and what the prepass cost and saved."""
    stats = dict(_stats(world))
    sources = ("consumption", "slab", "arc", "inherited", "fallback")
    decided = sum(stats[f"decided_{name}"] for name in sources)
    share = (lambda n: n / decided) if decided else (lambda n: 0.0)
    observed = sum(stats[f"evidence_obs_{name}"] for name in SOURCE_NAMES)
    deform_total = stats["deform_searches"] + stats["deform_searches_reused"]
    return {
        "counters": stats,
        "fronts_active": len(world.collision_polarity_frame.polarity) if world.collision_polarity_frame is not None else 0,
        "front_records": len(world.collision_fronts),
        "evidence_bins": int(sum(len(s) for s in world.collision_evidence.values())),
        "decision_share": {name: share(stats[f"decided_{name}"]) for name in sources},
        "scope_share": {name: share(stats[f"scope_{name}"]) for name in ("front", "record", "pair")},
        "fallback_rate": share(stats["decided_fallback"]),
        "ambiguity_rate": share(stats["ambiguous_physical"] + stats["ambiguous_opposing_arcs"]),
        "evidence_obs_share": {name: (stats[f"evidence_obs_{name}"] / observed if observed else 0.0) for name in SOURCE_NAMES},
        "churn": {key: stats[key] for key in ("fronts_created", "fronts_split", "fronts_merged", "conflicts_kept_apart", "fronts_expired", "fronts_dropped_topology", "fronts_split_topology", "fronts_fused")},
        "prepass_seconds_per_step": stats["prepass_seconds"] / stats["prepass_calls"] if stats["prepass_calls"] else 0.0,
        # Searches deform() would have run without the cache, and how many the cache answered.
        "deform_search_reuse": stats["deform_searches_reused"] / deform_total if deform_total else 0.0,
        "net_extra_searches_per_step": (
            (stats["prepass_searches"] - stats["deform_searches_reused"]) / stats["prepass_calls"] if stats["prepass_calls"] else 0.0
        ),
    }

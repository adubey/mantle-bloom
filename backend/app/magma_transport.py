"""Lateral magma transport (GitHub issue #205, follow-up to #120's "Land fraction slowly
declines"). Everything that thickens crust today (`apply_convergent_deformation`, the
near-field ring, arc magmatism, delamination melt intrusion -- all in rheology.py/
lithosphere_plate.py) either acts right at a collision boundary or, at most, a fixed ~350 km
ring on the *same* plate. Nothing transports mass from a plate actively being over-thickened by
collision to a plate (or a distant part of the same plate) losing land to stretching/thinning
elsewhere -- #120's own remaining land-loss driver.

This module is the transport/deposit half of the fix. The source half -- a small, continuous
skim of `apply_convergent_deformation`'s own strain increment at core convergent nodes,
diverted into a `MagmaParcel` instead of thickening the node in place -- lives in
`rheology.magma_export_strength_and_volume` (the mechanism reasoning, including why it's melt
and not diverted solid rock, and why the near-field ring is exempt, is documented there) and is
wired in by `lithosphere_plate.LithospherePlate.deform()`, which appends parcels to
`world.pending_magma_parcels`.

Parcels are banked across steps (real magmatic transport is a slow cumulative process, in the
same 150-450 m/Myr rate range as every other magmatic-addition path in this codebase) and
processed here at a fixed cadence (see `MAGMA_TRANSPORT_INTERVAL_STEPS`, same reasoning and
cadence as `gaps.GAP_FILL_INTERVAL_STEPS`'s own whole-sphere pass) by `world.step_world`, at the
same point in the step `gaps.py`'s own whole-sphere pass runs -- after this step's topology
changes have already settled, so a deposit is never targeted at a line that subducts, splits,
or merges away in the same step it's written.

A parcel is deliberately index-free (`{origin_xyz, volume_m3, step_generated,
unplaced_cycles}`, no `(plate_id, line_index, node_index)` reference): its origin's own local
strength was already reduced in the same step it was generated, so nothing about a banked
parcel ever needs to be re-addressed on the mesh, and there is no plate/line handle that could
go stale across a merge/split while it waits to fire. `origin_xyz` is frozen at generation and
not updated for the plate motion of the (few) banked steps before this pass runs -- a small,
deliberate approximation relative to `MAGMA_TRANSPORT_RANGE_KM`'s >=1000 km scale.

Design discussion (mechanism, data structures, what existing constructs might be
simplified/removed once this lands) went through two full review rounds on GitHub issue #205
before settling on this shape (referred to as "v3" there); see that issue's own comment history
for the reasoning behind each individual choice below, most of which exists specifically to
close a gap an earlier draft's review caught:

- Destinations are continental nodes only (`elevation_lines.effective_is_continental_from_codes`
  -- the correct per-node type resolution, not a plate-level or Hc-threshold check, since a
  plate's nominal crust_type can diverge from what an individual node actually is).
- Weighted toward thinner-than-reference continental crust with an inverse-distance falloff --
  both because thin lithosphere is mechanically easier to breach and because it directly
  targets the exact nodes #120 identifies as the land-loss mechanism.
- One *global*, per-destination-node rate cap, resolved after every parcel's own candidate
  destinations for the step are gathered -- not a per-parcel/per-source cap. A per-source cap
  alone doesn't stop several unrelated collisions worldwide from each (validly, by their own
  local budget) piling onto the same globally-thinnest node, which is a different route back to
  GitHub issue #145's "unbounded reappearance concentrated on a small target set" regression
  that `apply_delamination_melt_intrusion` already had to fix once.
- The deposit uses the same before/after isostasy-*delta* idiom `LithospherePlate.deform()`
  itself uses (see that method's own comment on why an elevation overwrite would silently
  launder/erase erosion history or existing transform/far-field debt already baked into a
  line's elevation), with per-node-resolved crust density (a destination can be a
  diverged-type patch, same reason the continental filter above is per-node not per-plate).
- Hc-only at both ends, no Hm coupling: this is genuinely mantle-derived melt (see
  `rheology.magma_export_strength_and_volume`'s own docstring for why), not relocated solid
  crust dragging its own mantle-lithosphere root along -- unlike
  `lithosphere_plate._redistribute_accreted_column`, which is a different case (relocating an
  already-solid retreated column) and does couple Hm.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import lithosphere, rheology
from .elevation_lines import (
    ELEV_CHANGE_LATERAL_MAGMA,
    ELEV_CHANGE_MIN_DELTA_M,
    PLANET_RADIUS_KM,
    effective_is_continental_from_codes,
    line_spacing_rad,
)

if TYPE_CHECKING:
    from .world import World

# Cadence: a whole-sphere lattice/tree pass, cheap but not free -- same reasoning and same
# cadence as gaps.GAP_FILL_INTERVAL_STEPS, which world.step_world calls it alongside.
MAGMA_TRANSPORT_INTERVAL_STEPS = 4

# Starting sweep value only, not a physically-shared quantity: this reuses FAR_FIELD_OUTER_KM's
# order of magnitude (lithosphere_plate.py) purely because it's a plausible "how far can a
# collision's mass reach" scale to begin a sweep from -- that constant represents a different
# process (elastic stress transmission backing a bare elevation delta), so the two are free to
# diverge once this is actually tuned (see #120-style toggle-sweep note on rheology.
# MAGMA_EXPORT_FRACTION).
MAGMA_TRANSPORT_RANGE_KM = 1000.0

# Same order of magnitude as rheology.DELAMINATION_MELT_INTRUSION_RATE_M_PER_MYR (300 m/Myr) --
# melt arriving faster than this can place it in one firing is not banked for later at that
# destination, same "lost, not banked" convention that rate already uses. Applied once per
# destination *node*, not summed across the whole world -- each destination has exactly one
# ceiling, enforced once, however many parcels (from however many unrelated collisions)
# independently targeted it this firing.
MAGMA_DEPOSIT_RATE_M_PER_MYR = 300.0

# A parcel that can't find rate-cap headroom at any qualifying destination for this many
# transport-pass firings running is dropped rather than banked indefinitely -- same "arrives too
# fast/too congested to place, is lost, not banked forever" convention
# rheology.apply_delamination_melt_intrusion's own rate cap already establishes, just applied
# across firings instead of within one.
MAGMA_PARCEL_MAX_AGE_CYCLES = 3

# A parcel with less than this much unplaced volume left is considered fully placed and
# dropped -- avoids keeping a parcel alive indefinitely over a floating-point dust remainder.
_NEGLIGIBLE_VOLUME_M3 = 1.0


@dataclass
class MagmaParcel:
    """One banked unit of exported convergent-boundary melt, not yet placed -- see module
    docstring for why this is deliberately index-free (no mesh address)."""

    origin_xyz: np.ndarray  # (3,) unit vector, frozen at generation
    volume_m3: float  # remaining unplaced volume
    step_generated: int  # diagnostic only -- nothing reads this back
    unplaced_cycles: int = 0  # consecutive transport-pass firings with zero placement


class _ContinentalNodeIndex:
    """Whole-sphere snapshot of every live continental node, addressable back to its own
    `(plate_id, line_index, node_index)` -- unlike `gaps.py`'s own whole-sphere tree (which
    only ever reads), this pass has to write back into whichever specific line a destination
    node lives on."""

    def __init__(
        self,
        xyz: np.ndarray,
        plate_id: np.ndarray,
        line_index: np.ndarray,
        node_index: np.ndarray,
        hc_m: np.ndarray,
    ) -> None:
        self.tree = cKDTree(xyz) if len(xyz) > 0 else None
        self.xyz = xyz
        self.plate_id = plate_id
        self.line_index = line_index
        self.node_index = node_index
        self.hc_m = hc_m

    def __len__(self) -> int:
        return len(self.xyz)


def _build_continental_node_index(world: "World") -> _ContinentalNodeIndex:
    """Mirrors `gaps._existing_node_tree`'s whole-sphere assembly, but per-line-addressable
    (needed for the scatter-write below) and pre-filtered to continental nodes only (this
    pass's destinations are never oceanic -- an oceanic destination is just ordinary seafloor
    volcanism, not the land-fraction fix this exists for)."""
    xyz_chunks, plate_id_chunks, line_index_chunks, node_index_chunks, hc_chunks = [], [], [], [], []
    for plate in world.plates:
        own_points, _ = plate.all_points_and_elevation()
        if len(own_points) == 0:
            continue
        plate_is_continental = plate.crust_type == "continental"
        offset = 0
        for line_index, line in enumerate(plate.lines):
            n = len(line)
            if n == 0:
                continue
            sl = slice(offset, offset + n)
            offset += n
            is_continental = effective_is_continental_from_codes(line.crust_type_code, plate_is_continental)
            if not np.any(is_continental):
                continue
            xyz_chunks.append(own_points[sl][is_continental])
            plate_id_chunks.append(np.full(int(np.count_nonzero(is_continental)), plate.plate_id))
            line_index_chunks.append(np.full(int(np.count_nonzero(is_continental)), line_index))
            node_index_chunks.append(np.flatnonzero(is_continental))
            hc_chunks.append(line.crustal_thickness_m[is_continental])
    if not xyz_chunks:
        return _ContinentalNodeIndex(np.zeros((0, 3)), np.zeros(0, dtype=int), np.zeros(0, dtype=int), np.zeros(0, dtype=int), np.zeros(0))
    return _ContinentalNodeIndex(
        np.concatenate(xyz_chunks, axis=0),
        np.concatenate(plate_id_chunks),
        np.concatenate(line_index_chunks),
        np.concatenate(node_index_chunks),
        np.concatenate(hc_chunks),
    )


def _weighted_destination_pairs(
    parcel_origins_xyz: np.ndarray, dest_index: _ContinentalNodeIndex, range_rad: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Every (parcel, destination) pair within `range_rad` with nonzero weight, as three flat
    arrays `(parcel_idx, dest_idx, weight)` -- weight is thinness-relative-to-reference
    (`REFERENCE_HC_CONTINENTAL_M - hc`, floored at 0) times a linear inverse-distance falloff to
    zero at `MAGMA_TRANSPORT_RANGE_KM`, the same ramp shape
    `lithosphere_plate._far_field_intensity` already uses elsewhere in this codebase. Distance
    is the plain cKDTree chord distance over unit-sphere points, times `PLANET_RADIUS_KM` -- the
    same "chord stands in for great-circle arc length" approximation the far-field collision
    code (lithosphere_plate.py) already relies on at this same ~1000 km scale (well under 1%
    error there)."""
    origin_tree = cKDTree(parcel_origins_xyz)
    candidates = origin_tree.query_ball_tree(dest_index.tree, range_rad)

    parcel_idx_chunks, dest_idx_chunks, weight_chunks = [], [], []
    for p, dest_candidates in enumerate(candidates):
        if not dest_candidates:
            continue
        dest_candidates = np.array(dest_candidates)
        dist_km = np.linalg.norm(parcel_origins_xyz[p] - dest_index.xyz[dest_candidates], axis=1) * PLANET_RADIUS_KM
        thinness = np.clip(lithosphere.REFERENCE_HC_CONTINENTAL_M - dest_index.hc_m[dest_candidates], 0.0, None)
        falloff = np.clip(1.0 - dist_km / MAGMA_TRANSPORT_RANGE_KM, 0.0, 1.0)
        weight = thinness * falloff
        keep = weight > 0.0
        if not np.any(keep):
            continue
        parcel_idx_chunks.append(np.full(int(np.count_nonzero(keep)), p))
        dest_idx_chunks.append(dest_candidates[keep])
        weight_chunks.append(weight[keep])

    if not parcel_idx_chunks:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int), np.zeros(0)
    return np.concatenate(parcel_idx_chunks), np.concatenate(dest_idx_chunks), np.concatenate(weight_chunks)


def run_magma_transport(world: "World", banked_myr: float) -> list[str]:
    """Resolve every banked `world.pending_magma_parcels` entry against the current whole-sphere
    continental-node set, apply one global per-destination-node rate cap, scatter-write the
    realized deposits, and update `world.pending_magma_parcels` in place (partial placements
    keep their remaining volume banked; fully-placed or stale parcels are dropped -- see module
    docstring / `MAGMA_PARCEL_MAX_AGE_CYCLES`). Returns event strings for `world.log_event`;
    a no-op (returns `[]`, mutates nothing) when there are no pending parcels."""
    parcels = world.pending_magma_parcels
    if not parcels:
        return []

    dest_index = _build_continental_node_index(world)
    if len(dest_index) == 0 or dest_index.tree is None:
        # Nothing to deposit onto at all -- every parcel just ages a cycle.
        for parcel in parcels:
            parcel.unplaced_cycles += 1
        world.pending_magma_parcels = [p for p in parcels if p.unplaced_cycles < MAGMA_PARCEL_MAX_AGE_CYCLES]
        return []

    spacing_rad = line_spacing_rad(world.node_density)
    node_area_m2 = lithosphere.node_area_m2(spacing_rad)
    range_rad = MAGMA_TRANSPORT_RANGE_KM / PLANET_RADIUS_KM

    parcel_origins = np.array([p.origin_xyz for p in parcels])
    parcel_volume = np.array([p.volume_m3 for p in parcels])
    parcel_idx, dest_idx, weight = _weighted_destination_pairs(parcel_origins, dest_index, range_rad)

    placed_per_parcel = np.zeros(len(parcels))
    n_deposited = 0
    total_deposited_m = 0.0

    if len(parcel_idx) > 0:
        sum_weight_per_parcel = np.zeros(len(parcels))
        np.add.at(sum_weight_per_parcel, parcel_idx, weight)
        share = weight / sum_weight_per_parcel[parcel_idx]
        requested_hc = parcel_volume[parcel_idx] * share / node_area_m2

        total_requested_hc = np.zeros(len(dest_index))
        np.add.at(total_requested_hc, dest_idx, requested_hc)
        cap_hc = max(MAGMA_DEPOSIT_RATE_M_PER_MYR * banked_myr, 0.0)
        # Also bounded by each destination's own headroom to MAX_CRUSTAL_THICKNESS_M (the same
        # ceiling `_scatter_write_deposits` clips its actual Hc write to). Ordinary weighting
        # already keeps destinations well under this -- weight goes to 0 at
        # REFERENCE_HC_CONTINENTAL_M, far below the ceiling -- but an unusually large
        # `banked_myr` (an outsized step's own `years`) could otherwise make `cap_hc` promise
        # more than a near-ceiling destination has room for; without this, the promised amount
        # (used below for `placed_per_parcel` and the logged deposited total) would silently
        # drift from what `_scatter_write_deposits` actually writes, losing mass without
        # accounting for it -- unlike every other bounded-rate melt path in this codebase
        # (e.g. rheology.apply_delamination_melt_intrusion), where "can't be placed" is always
        # reflected in what the caller is told was placed.
        headroom_hc = np.clip(lithosphere.MAX_CRUSTAL_THICKNESS_M - dest_index.hc_m, 0.0, None)
        realized_at_dest = np.minimum(np.minimum(total_requested_hc, cap_hc), headroom_hc)
        safe_total = np.where(total_requested_hc > 0.0, total_requested_hc, 1.0)
        ratio_at_dest = np.where(total_requested_hc > 0.0, realized_at_dest / safe_total, 0.0)

        realized_hc_per_pair = requested_hc * ratio_at_dest[dest_idx]
        np.add.at(placed_per_parcel, parcel_idx, realized_hc_per_pair * node_area_m2)

        deposit_mask = realized_at_dest > 0.0
        n_deposited = int(np.count_nonzero(deposit_mask))
        total_deposited_m = float(np.sum(realized_at_dest[deposit_mask]))
        if n_deposited > 0:
            _scatter_write_deposits(world, dest_index, deposit_mask, realized_at_dest)

    new_parcels: list[MagmaParcel] = []
    for i, parcel in enumerate(parcels):
        remaining = max(parcel_volume[i] - placed_per_parcel[i], 0.0)
        parcel.unplaced_cycles = 0 if placed_per_parcel[i] > 0.0 else parcel.unplaced_cycles + 1
        if remaining <= _NEGLIGIBLE_VOLUME_M3:
            continue  # fully placed -- done
        if parcel.unplaced_cycles >= MAGMA_PARCEL_MAX_AGE_CYCLES:
            continue  # stale -- couldn't find headroom for MAGMA_PARCEL_MAX_AGE_CYCLES firings running
        parcel.volume_m3 = remaining
        new_parcels.append(parcel)
    world.pending_magma_parcels = new_parcels

    if n_deposited == 0:
        return []
    return [f"Lateral magma transport deposited {total_deposited_m:.0f} m of crust across {n_deposited} nodes."]


def _scatter_write_deposits(
    world: "World", dest_index: _ContinentalNodeIndex, deposit_mask: np.ndarray, realized_at_dest: np.ndarray
) -> None:
    """Group the capped per-destination-node deposits by `(plate_id, line_index)` and issue one
    `replace_line` per touched line -- `PlateWithLines.replace_line` (plates.py) already exists
    as the mutate-one-line-by-index primitive; this is what's new is the *dispatch* across
    plates to reach it, since every existing cKDTree use in lithosphere_plate.py is a read-only
    query into another plate's own node cloud, never a write. Reuses `deform()`'s own
    before/after isostasy-delta idiom (see module docstring) -- never a raw isostasy overwrite,
    which would silently erase erosion history or existing unbacked-relief debt already baked
    into a line's elevation."""
    plates_by_id = {plate.plate_id: plate for plate in world.plates}
    deposit_indices = np.flatnonzero(deposit_mask)

    groups: dict[tuple[int, int], list[tuple[int, float]]] = {}
    for d in deposit_indices:
        key = (int(dest_index.plate_id[d]), int(dest_index.line_index[d]))
        groups.setdefault(key, []).append((int(dest_index.node_index[d]), float(realized_at_dest[d])))

    for (plate_id, line_index), entries in groups.items():
        plate = plates_by_id.get(plate_id)
        if plate is None or line_index >= len(plate.lines):
            continue  # defensive -- shouldn't happen, this pass runs after topology settles
        line = plate.lines[line_index]
        node_idx = np.array([e[0] for e in entries])
        delta_hc = np.array([e[1] for e in entries])

        hc = line.crustal_thickness_m.copy()
        hm = line.mantle_lithosphere_thickness_m
        rho_c = lithosphere.node_crust_density(line.crust_type_code[node_idx], plate.crust_type)
        elevation_before = lithosphere.isostatic_elevation(hc[node_idx], hm[node_idx], rho_c)
        new_hc_touched = np.minimum(hc[node_idx] + delta_hc, lithosphere.MAX_CRUSTAL_THICKNESS_M)
        elevation_after = lithosphere.isostatic_elevation(new_hc_touched, hm[node_idx], rho_c)
        hc[node_idx] = new_hc_touched

        elevation = line.elevation.copy()
        new_elevation_touched = rheology.clip_elevation_bounds(elevation[node_idx] + (elevation_after - elevation_before))
        moved = np.abs(new_elevation_touched - elevation[node_idx]) >= ELEV_CHANGE_MIN_DELTA_M
        elevation[node_idx] = new_elevation_touched

        reason = line.elev_change_reason.copy()
        reason[node_idx[moved]] = ELEV_CHANGE_LATERAL_MAGMA

        plate.replace_line(line_index, line.replace(crustal_thickness_m=hc, elevation=elevation, elev_change_reason=reason))

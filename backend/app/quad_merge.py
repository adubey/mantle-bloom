"""Cross-plate merge for `PlateWithSparseQuadPatch` plates (issue #228 Phase 4).

The line engine fuses two plates by resampling both node clouds onto a fresh lattice in the
surviving plate's frame (`lithosphere_plate._merge_lines_from_resample`), carrying only Hc/Hm
by nearest neighbour and re-deriving elevation. Nearest-neighbour sampling at a constant
nominal node area can't conserve volume, and the resample discards every other field.

Quad plates keep the surviving plate's cells and IDs as they are and remap the absorbed plate
onto that plate's lattice by exact area:

- **Frame remap.** Each absorbed cell is split into `MERGE_REMAP_SUBSAMPLES`^2 sub-cells on
  its own lattice, whose exact solid angles sum to the cell's own area. Each sub-cell centre
  is carried through world space into the surviving plate's local frame and located on its
  lattice. The sub-cell areas form the overlap weights between absorbed and surviving cells.
- **Territory.** An empty surviving-lattice cell becomes active when absorbed sub-cells cover
  at least `NEW_CELL_MIN_COVERAGE` of its exact area and it isn't closer to some third
  plate's node than to either fusing plate's (the Voronoi exclusivity the line merge uses).
- **Transfer.** A new cell takes the area-weighted remap of the sub-cells inside it. For
  extensive fields that's the mean thickness, so a partly covered boundary cell isn't left
  thin. The volume of sub-cells whose cell wasn't activated, plus the coverage difference,
  is restored by one ratio per field across the new cells, so Hc/Hm and every other
  extensive field are conserved exactly by cell area. Where the absorbed plate overlapped
  the surviving one (the suture), its volume stacks onto the surviving column, with Hc
  capped at `SUTURE_ACCRETION_MAX_HC_M` as in the line merge. That cap, and Hm's own
  ceiling, are the only places volume can leave. Every other field follows its
  `surface_fields.RemapClass`, with the same rules `coarsen_cells` uses.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, lithosphere, rheology, torque
from .elevation_lines import CRUST_TYPE_CONTINENTAL, CRUST_TYPE_INHERIT, CRUST_TYPE_OCEANIC
from .lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M
from .surface_fields import SURFACE_FIELDS, RemapClass

if TYPE_CHECKING:
    from .sparse_quad_patch import PlateWithSparseQuadPatch

# Sub-cells per absorbed-cell side for the frame remap. Conservation doesn't depend on it
# (the sub-cell areas sum to each cell's exact area whatever the count); it sets how finely
# partial coverage of a surviving-lattice cell is resolved -- 1/16 of a cell at 4.
MERGE_REMAP_SUBSAMPLES = 4

# Fraction of an empty surviving-lattice cell's exact area the absorbed plate must cover for
# the cell to become part of the merged plate. A half keeps the merged footprint's area close
# to the absorbed plate's: the two lattices are rotated against each other, so boundary cells
# are partly covered, and those under half cover give their share to a neighbour instead.
NEW_CELL_MIN_COVERAGE = 0.5


def merge(keep: "PlateWithSparseQuadPatch", absorb: "PlateWithSparseQuadPatch", other_points_xyz: np.ndarray) -> None:
    """Fuse `absorb` into `keep` in place: `keep`'s frame, cells, and cell IDs survive, and
    `absorb`'s territory and fields are remapped onto `keep`'s lattice by exact area (see the
    module docstring). `omega` is the angular-momentum-conserving blend of both plates,
    with inertia from exact cell areas. `other_points_xyz` are every other live plate's
    nodes (world xyz), which `keep` doesn't grow over."""
    inertia_keep = _inertia(keep)
    inertia_absorb = _inertia(absorb)
    if absorb.node_count():
        _transfer(keep, absorb, np.asarray(other_points_xyz, dtype=float).reshape(-1, 3))
    keep.set_omega(torque.merge_omega(keep, inertia_keep, absorb, inertia_absorb))
    keep.reset_age()


def _inertia(plate: "PlateWithSparseQuadPatch") -> np.ndarray:
    return lithosphere.moment_of_inertia_tensor(
        plate.all_points_and_elevation()[0],
        plate.collect("crustal_thickness_m"),
        plate.collect("mantle_lithosphere_thickness_m"),
        lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type),
        0.0,
        area_m2=plate.node_areas_m2(),
    )


def subsample_cells(plate: "PlateWithSparseQuadPatch", per_side: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(source node, local centre, exact area in m^2) of each of `per_side`^2 sub-cells per
    active cell, on the plate's own lattice. Each cell's sub-cell areas are normalised to sum
    to its `node_areas_m2` exactly, so any partition of the sub-cells partitions the plate's
    area (and hence its volume) with no residual beyond floating-point summation."""
    from .sparse_quad_patch import PLANET_RADIUS_M, cell_areas_sr, lattice_points, unpack_cell_keys

    face, level, i, j = unpack_cell_keys(plate.cell_keys)
    count = per_side * per_side
    da, db = (offset.ravel() for offset in np.meshgrid(np.arange(per_side), np.arange(per_side), indexing="ij"))
    source = np.repeat(np.arange(len(face)), count)
    sub_face = face[source]
    sub_level = level[source]
    sub_i = i[source] * per_side + np.tile(da, len(face))
    sub_j = j[source] * per_side + np.tile(db, len(face))
    local = np.zeros((len(source), 3))
    area = np.zeros(len(source))
    for lev in np.unique(sub_level):
        at = sub_level == lev
        resolution = plate.cells_per_edge * (1 << int(lev)) * per_side
        local[at] = lattice_points(sub_face[at], sub_i[at] + 0.5, sub_j[at] + 0.5, resolution)
        area[at] = cell_areas_sr(sub_i[at], sub_j[at], resolution) * PLANET_RADIUS_M**2
    area *= (plate.node_areas_m2() / np.bincount(source, weights=area, minlength=len(face)))[source]
    return source, local, area


def _free_cell_keys(plate: "PlateWithSparseQuadPatch", local: np.ndarray) -> np.ndarray:
    """For local points in no active leaf of `plate`: the coarsest lattice cell containing
    each that overlaps no active leaf. Level 0 unless the plate is refined nearby; at the
    deepest active level it always exists, since the point isn't in any active leaf there."""
    from .sparse_quad_patch import locate_cells, pack_cell_keys, unpack_cell_keys

    keys = np.full(len(local), -1, dtype=np.int64)
    max_level = int(unpack_cell_keys(plate.cell_keys)[1].max(initial=0))
    for lev in range(max_level + 1):
        todo = np.flatnonzero(keys < 0)
        if not len(todo):
            break
        face, i, j = locate_cells(local[todo], plate.cells_per_edge * (1 << lev))
        candidate = pack_cell_keys(face, i, j, lev)
        free = ~plate._overlaps_active(candidate)
        keys[todo[free]] = candidate[free]
    return keys


def _explicit_codes(codes: np.ndarray, crust_type: str) -> np.ndarray:
    """`crust_type_code` with CRUST_TYPE_INHERIT resolved against the owning plate's type."""
    inherited = CRUST_TYPE_CONTINENTAL if crust_type == "continental" else CRUST_TYPE_OCEANIC
    return np.where(codes == CRUST_TYPE_INHERIT, inherited, codes).astype(codes.dtype)


def _transfer(keep: "PlateWithSparseQuadPatch", absorb: "PlateWithSparseQuadPatch", other_points_xyz: np.ndarray) -> None:
    source, absorb_local, sub_area = subsample_cells(absorb, MERGE_REMAP_SUBSAMPLES)
    sub_world = geometry.to_world(absorb.frame, absorb_local)
    sub_local = geometry.to_local(keep.frame, sub_world)

    # Sub-cells under an existing `keep` leaf (the suture) versus open lattice.
    in_keep = keep._leaf_key_at(sub_local)
    fresh = in_keep < 0
    candidate = _free_cell_keys(keep, sub_local[fresh])
    cell_keys, cell_of = np.unique(candidate, return_inverse=True)
    coverage = np.bincount(cell_of, weights=sub_area[fresh], minlength=len(cell_keys))
    chosen = coverage >= NEW_CELL_MIN_COVERAGE * _cell_areas_m2(keep, cell_keys)
    if len(other_points_xyz) and np.any(chosen):
        centres = geometry.to_world(keep.frame, keep.cell_centres_local(cell_keys[chosen]))
        own = np.concatenate([keep.all_points_and_elevation()[0], absorb.all_points_and_elevation()[0]])
        own_dist, _ = cKDTree(own).query(centres)
        other_dist, _ = cKDTree(other_points_xyz).query(centres)
        chosen[np.flatnonzero(chosen)[other_dist < own_dist]] = False

    old_keys = keep.cell_keys.copy()
    old_fields = {name: keep.collect(name) for name in _materialised(keep, absorb)}
    inserted = keep.insert_cells(cell_keys[chosen])
    new_keys = cell_keys[chosen][inserted]

    # Target node (in the merged order) of every sub-cell. Orphans -- sub-cells of an open
    # cell that wasn't activated -- are given the nearest cell, but only their extensive volume
    # is used (see `extensive` below); values are remapped from mapped sub-cells alone.
    target = np.full(len(source), -1, dtype=np.int64)
    target[~fresh] = keep._index_of_keys(in_keep[~fresh])
    target[fresh] = keep._index_of_keys(candidate)
    orphan = target < 0
    if np.any(orphan):
        _, nearest = cKDTree(keep.surface_nodes().local_xyz).query(sub_local[orphan])
        target[orphan] = nearest
    mapped = ~orphan

    n = keep.node_count()
    is_new = np.isin(keep.cell_keys, new_keys)
    old_index = np.searchsorted(old_keys, keep.cell_keys[~is_new])
    target_area = keep.node_areas_m2()
    mapped_area = np.bincount(target[mapped], weights=sub_area[mapped], minlength=n)

    absorb_fields = {name: absorb.collect(name) for name in old_fields}
    absorb_fields["crust_type_code"] = _relative_codes(
        _explicit_codes(absorb_fields["crust_type_code"], absorb.crust_type), keep.crust_type
    )

    def summed(values: np.ndarray, weights: np.ndarray, among: np.ndarray = mapped) -> np.ndarray:
        """Per merged node, sum of `weights * values` over the sub-cells `among` it receives."""
        return np.bincount(target[among], weights=(weights * values)[among], minlength=n)

    def weighted_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
        total = summed(np.ones(len(values)), weights)
        return np.divide(summed(values, weights), total, out=np.zeros(n), where=total > 0.0)

    def extensive(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(value on new cells, thickness added to existing cells) for an extensive field.
        A suture cell stacks the volume of the sub-cells under it. A new cell takes the
        area-weighted mean of its own sub-cells, so a boundary cell that is only partly
        covered isn't left thin; the orphans' volume and that coverage difference are then
        restored by one ratio across the new cells, which makes the total exact without
        putting stray volume in any one cell. With no new cell to take it, orphan volume
        stacks onto the nearest cell instead."""
        volume = summed(values, sub_area)
        added = np.where(is_new, 0.0, volume) / target_area
        mean = np.divide(volume, mapped_area, out=np.zeros(n), where=mapped_area > 0.0)
        placed = float(np.sum((mean * target_area)[is_new]))
        owed = float(np.sum(volume[is_new])) + float(np.sum((values * sub_area)[orphan]))
        if placed > 0.0:
            return mean * (owed / placed), added
        orphan_volume = summed(values, sub_area, orphan) / target_area
        return orphan_volume, added + np.where(is_new, 0.0, orphan_volume)

    out: dict[str, np.ndarray] = {}
    added_hc = np.zeros(n)
    for name in old_fields:
        spec = SURFACE_FIELDS[name]
        merged = np.full(n, spec.default, dtype=spec.dtype)
        merged[~is_new] = old_fields[name][old_index]
        values = absorb_fields[name][source]
        if spec.remap_class == RemapClass.EXTENSIVE:
            new_values, added = extensive(values)
            merged[is_new] = new_values[is_new]
            merged[~is_new] += added[~is_new]
            if name == "crustal_thickness_m":
                added_hc = added
        elif name == "channel_depth" or spec.remap_class == RemapClass.COUNTDOWN:
            peak = np.full(n, -np.inf)
            np.maximum.at(peak, target[mapped], values[mapped])
            merged[is_new] = peak[is_new]
        elif spec.remap_class in (RemapClass.HISTORY, RemapClass.WRITE_ONCE_HISTORY):
            valid = mapped & (values != spec.sentinel)
            earliest = np.full(n, np.inf)
            np.minimum.at(earliest, target[valid], values[valid])
            merged[is_new] = np.where(np.isfinite(earliest), earliest, spec.sentinel)[is_new]
        elif spec.remap_class == RemapClass.BOOLEAN_PROVENANCE:
            merged[is_new] = (summed(values.astype(float), np.ones(len(values))) > 0.0)[is_new]
        elif spec.remap_class == RemapClass.CATEGORICAL:
            choices = np.unique(values)
            votes = np.stack([summed((values == c).astype(float), sub_area) for c in choices])
            merged[is_new] = choices[np.argmax(votes, axis=0)][is_new]
        elif name != "elevation":
            weights = sub_area
            coupled = {"soil_mineral_content": "soil_depth", "soil_organic_content": "soil_depth", "channel_width": "channel_depth"}
            if name in coupled:
                weights = sub_area * absorb_fields[coupled[name]][source]
            merged[is_new] = weighted_mean(values, weights)[is_new]
        out[name] = merged

    # Crust: the suture cap applies to stacked columns only, never lowers a column already
    # past it, and caps Hm the way `quad_tectonics._accrete_onto_survivors` does.
    stacked = ~is_new & (added_hc > 0.0)
    old_hc = np.zeros(n)
    old_hm = np.zeros(n)
    old_hc[~is_new] = old_fields["crustal_thickness_m"][old_index]
    old_hm[~is_new] = old_fields["mantle_lithosphere_thickness_m"][old_index]
    hc, hm = out["crustal_thickness_m"], out["mantle_lithosphere_thickness_m"]
    hc[stacked] = np.maximum(old_hc[stacked], np.minimum(hc[stacked], SUTURE_ACCRETION_MAX_HC_M))
    hm[stacked] = np.maximum(old_hm[stacked], np.minimum(hm[stacked], lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M))

    # Elevation: isostasy from the new column plus the carried erosion/texture residual, as in
    # `coarsen_cells`; a stacked column keeps its own residual and moves by the isostatic change.
    density = lithosphere.node_crust_density(out["crust_type_code"], keep.crust_type)
    absorb_residual = absorb_fields["elevation"] - lithosphere.isostatic_elevation(
        absorb.collect("crustal_thickness_m"),
        absorb.collect("mantle_lithosphere_thickness_m"),
        lithosphere.node_crust_density(absorb.collect("crust_type_code"), absorb.crust_type),
    )
    residual = weighted_mean(absorb_residual[source], sub_area)
    elevation = out["elevation"]
    new_isostatic = lithosphere.isostatic_elevation(hc, hm, density)
    elevation[is_new] = new_isostatic[is_new] + residual[is_new]
    elevation[stacked] += new_isostatic[stacked] - lithosphere.isostatic_elevation(old_hc[stacked], old_hm[stacked], density[stacked])
    out["elevation"] = rheology.clip_elevation_bounds(elevation)
    keep.set_fields_on_plate(**out)


def _materialised(keep: "PlateWithSparseQuadPatch", absorb: "PlateWithSparseQuadPatch") -> list[str]:
    """Every field either plate has stored, plus the column fields the transfer always
    writes -- a field neither plate has touched stays at its registry default for free."""
    names = set(keep._fields) | set(absorb._fields)
    names |= {"elevation", "crustal_thickness_m", "mantle_lithosphere_thickness_m", "crust_type_code"}
    return sorted(names)


def _relative_codes(explicit: np.ndarray, crust_type: str) -> np.ndarray:
    """Explicit crust codes re-encoded for a plate of `crust_type`: its own type becomes
    CRUST_TYPE_INHERIT again, the other stays explicit."""
    own = CRUST_TYPE_CONTINENTAL if crust_type == "continental" else CRUST_TYPE_OCEANIC
    return np.where(explicit == own, CRUST_TYPE_INHERIT, explicit).astype(explicit.dtype)


def _cell_areas_m2(plate: "PlateWithSparseQuadPatch", keys: np.ndarray) -> np.ndarray:
    """Exact area of each (not necessarily active) cell of `plate`'s lattice."""
    from .sparse_quad_patch import PLANET_RADIUS_M, cell_areas_sr, unpack_cell_keys

    _, level, i, j = unpack_cell_keys(keys)
    areas = np.zeros(len(keys))
    for lev in np.unique(level):
        at = level == lev
        areas[at] = cell_areas_sr(i[at], j[at], plate.cells_per_edge * (1 << int(lev)))
    return areas * PLANET_RADIUS_M**2

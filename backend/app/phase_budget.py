"""Per-phase Hc/Hm budget instrumentation for GitHub issue #216 ("Investigate long-run Hc/Hm
decline: geological sinks, changing node counts, and numerical volume losses").

Gated entirely by `World.debug_diagnostics` so it costs nothing on an ordinary world. Each
call site that mutates crustal thickness (Hc) / mantle-lithosphere thickness (Hm) --
convergent/divergent deformation, arc magmatism, oceanic cooling relaxation, decompression
melting, boundary advance/retreat, gap fill, plate merges/cleanup, failed rifts, erosion, the
end-of-step cap clamp -- calls `record()` with the same node slice's Hc/Hm (and `crust_type_code`) just
before and just after that phase ran. `record()` accumulates the delta into `World.phase_budget[phase_name]`, broken out
by plate type (the owning plate's own `crust_type`) and by node type (`crust_type_code`,
resolved via `elevation_lines.effective_is_continental_from_codes` -- a node can carry a type
different from its plate's nominal one, e.g. an accreted terrane), so a run can separate a
genuine reservoir loss/gain in one mechanism from population change (more/fewer nodes) or
reclassification (nodes changing type without any Hc/Hm change at all).

`World.phase_budget` is a plain cumulative total since the world was created (or since the
last `World.reset_phase_budget()`) -- there is no separate per-step history, by design: a
caller investigating a specific interval (e.g. a short replay from a saved world) resets the
budget, steps forward, and reads it back, rather than the engine paying to retain a full
per-step time series it may never need.

Each scope carries node counts and plain per-node Hc/Hm sums (`sum_hc_*`, `sum_hm_*`), plus the
covered area and Hc/Hm volumes weighted by each node's own accounting area
(`Plate.accounting_areas_m2`: exact cell areas). Quad cells are not equal-area and a phase can
swap large cells for small ones, so only the area-weighted fields measure real crust gained or
lost (issue #257). They have to be weighted per node as each call is recorded -- the aggregate
counts can't be converted afterwards. A caller passing no areas gets the nominal area per
node.

One entry is not an Hc/Hm phase: `SHORTENING_PHASE` (issue #314) books a quad plate's
collisional shortening -- demanded by its convergent band, absorbed into its columns, or
returned to the boundary overlap -- as shortening area (m^2: metres of convergence times
metres of boundary) and as the Hc/Hm volume it builds or would have built. It carries a
`shortening` dict instead of `scopes`; the absorbed share's Hc/Hm change is already inside
`convergent_deformation`.
"""

from __future__ import annotations

import numpy as np

from typing import NamedTuple

from . import elevation_lines, hm_ledger, lithosphere

_SCOPE_FIELDS = (
    "count_before",
    "count_after",
    "sum_hc_before",
    "sum_hc_after",
    "sum_hm_before",
    "sum_hm_after",
    "area_before_m2",
    "area_after_m2",
    "hc_volume_before_m3",
    "hc_volume_after_m3",
    "hm_volume_before_m3",
    "hm_volume_after_m3",
)
# "all" = every node the call touched; "*_plate" splits by the owning plate's own `crust_type`
# (one bucket per call, the whole slice belongs to one plate); "*_node" splits by each node's
# own *effective* type (crust_type_code resolved against the plate), which can disagree with
# the plate bucket for e.g. an accreted terrane or a freshly-melted node stamped the other way.
SCOPES = ("all", "continental_plate", "oceanic_plate", "continental_node", "oceanic_node")

SHORTENING_PHASE = "convergent_shortening"
SHORTENING_FIELDS = (
    "demand_m2",
    "absorbed_m2",
    "returned_decay_m2",
    "returned_unrouted_m2",
    "absorbed_hc_m3",
    "absorbed_hm_m3",
    "returned_hc_m3",
    "returned_hm_m3",
)

CAP_SCOPES = (
    "all",
    "continental_node",
    "oceanic_node",
    "continental_craton_node",
    "continental_non_craton_node",
)
_CAP_FIELDS = (
    "newly_capped_area_m2",
    "persistently_capped_area_m2",
    "left_cap_area_m2",
)
# These snapshots can change the physical identity behind a node ID. Do not infer per-cell
# cap entries/exits from a coincidental ID match across a merge.
_UNALIGNED_CAP_TRANSITION_PHASES = {"plate_merge", "forced_plate_merge"}


class Snapshot(NamedTuple):
    """One plate's columns at an instant, in node order, with each node's accounting area."""

    hc: np.ndarray
    hm: np.ndarray
    codes: np.ndarray
    area_m2: np.ndarray
    craton_m: np.ndarray
    node_ids: np.ndarray
    owning_plate_is_continental: np.ndarray


def snapshot(plate, spacing_rad: float) -> Snapshot:
    """Hc, Hm, crust_type_code and `Plate.accounting_areas_m2`, each concatenated across the
    whole plate in node order (`Plate.collect`) -- for a phase that changes topology or works
    plate-wide (merge, cleanup removal, retreat/advance, gap fill, erosion, the cap clamp),
    where `record`'s per-node before/after arrays can't capture the whole effect."""
    nodes = plate.surface_nodes()
    return Snapshot(
        plate.collect("crustal_thickness_m"),
        plate.collect("mantle_lithosphere_thickness_m"),
        plate.collect("crust_type_code"),
        plate.accounting_areas_m2(spacing_rad),
        plate.collect("craton_crust_m"),
        nodes.node_ids,
        np.full(len(nodes.node_ids), plate.crust_type == "continental", dtype=bool),
    )


def record_snapshots(world, plate, phase: str, before: Snapshot, after: Snapshot) -> None:
    """`record` for two `snapshot`s of the same plate."""
    record(
        world, plate, phase, before.hc, before.hm, before.codes, after.hc, after.hm, after.codes,
        area_before_m2=before.area_m2, area_after_m2=after.area_m2,
        craton_before_m=before.craton_m, craton_after_m=after.craton_m,
        node_ids_before=before.node_ids, node_ids_after=after.node_ids,
        plate_is_continental_before=before.owning_plate_is_continental,
        plate_is_continental_after=after.owning_plate_is_continental,
    )


def _new_totals() -> dict:
    return {
        "calls": 0,
        "scopes": {scope: dict.fromkeys(_SCOPE_FIELDS, 0.0) for scope in SCOPES},
        "hm_cap_transitions": {
            "aligned_calls": 0,
            "unaligned_calls": 0,
            "scopes": {scope: dict.fromkeys(_CAP_FIELDS, 0.0) for scope in CAP_SCOPES},
        },
    }


def _id_rows(ids: np.ndarray) -> np.ndarray:
    ids = np.ascontiguousarray(ids)
    if ids.ndim == 1:
        ids = ids[:, None]
    return ids.view(np.dtype((np.void, ids.dtype.itemsize * ids.shape[1]))).reshape(-1)


def _record_cap_transitions(
    totals: dict,
    plate_is_continental: bool,
    hm_before: np.ndarray,
    hm_after: np.ndarray,
    codes_before: np.ndarray,
    codes_after: np.ndarray,
    area_before: np.ndarray,
    area_after: np.ndarray,
    craton_before: np.ndarray | None,
    craton_after: np.ndarray | None,
    node_ids_before: np.ndarray | None,
    node_ids_after: np.ndarray | None,
    plate_is_continental_before: bool | np.ndarray,
    plate_is_continental_after: bool | np.ndarray,
    force_unaligned: bool = False,
) -> None:
    """Book actual per-cell entries/exits at the Hm cap; never infer them from net volume."""
    transitions = totals["hm_cap_transitions"]
    if force_unaligned:
        # The phase has no trustworthy one-to-one cell correspondence. Conservatively treat
        # every capped after-cell as new and every capped before-cell as removed; persistent
        # coverage is not claimed across a topology identity reset.
        transitions["unaligned_calls"] += 1
        before_index = after_index = np.zeros(0, dtype=int)
        new_after = np.ones(len(hm_after), dtype=bool)
        removed_before = np.ones(len(hm_before), dtype=bool)
    elif node_ids_before is None or node_ids_after is None:
        if len(hm_before) != len(hm_after):
            transitions["unaligned_calls"] += 1
            return
        before_index = after_index = np.arange(len(hm_before))
        new_after = np.zeros(len(hm_after), dtype=bool)
        removed_before = np.zeros(len(hm_before), dtype=bool)
    else:
        before_ids = _id_rows(np.asarray(node_ids_before))
        after_ids = _id_rows(np.asarray(node_ids_after))
        if np.array_equal(before_ids, after_ids):
            before_index = after_index = np.arange(len(before_ids))
            new_after = np.zeros(len(hm_after), dtype=bool)
            removed_before = np.zeros(len(hm_before), dtype=bool)
        elif len(np.unique(before_ids)) != len(before_ids) or len(np.unique(after_ids)) != len(after_ids):
            transitions["unaligned_calls"] += 1
            return
        else:
            _, before_index, after_index = np.intersect1d(before_ids, after_ids, return_indices=True)
            new_after = np.ones(len(hm_after), dtype=bool)
            new_after[after_index] = False
            removed_before = np.ones(len(hm_before), dtype=bool)
            removed_before[before_index] = False

    if not force_unaligned:
        transitions["aligned_calls"] += 1
    before_capped = hm_before >= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M - 1.0
    after_capped = hm_after >= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M - 1.0
    newly_after = new_after & after_capped
    persistent_after = np.zeros(len(hm_after), dtype=bool)
    left_before = removed_before & before_capped
    if len(before_index):
        newly_after[after_index] |= ~before_capped[before_index] & after_capped[after_index]
        persistent_after[after_index] = before_capped[before_index] & after_capped[after_index]
        left_before[before_index] |= before_capped[before_index] & ~after_capped[after_index]

    cont_before = elevation_lines.effective_is_continental_from_codes(codes_before, plate_is_continental_before)
    cont_after = elevation_lines.effective_is_continental_from_codes(codes_after, plate_is_continental_after)
    craton_before_mask = np.zeros(len(hm_before), dtype=bool) if craton_before is None else np.asarray(craton_before) > 0.0
    craton_after_mask = np.zeros(len(hm_after), dtype=bool) if craton_after is None else np.asarray(craton_after) > 0.0
    masks = {
        "all": (np.ones(len(hm_before), dtype=bool), np.ones(len(hm_after), dtype=bool)),
        "continental_node": (cont_before, cont_after),
        "oceanic_node": (~cont_before, ~cont_after),
        "continental_craton_node": (cont_before & craton_before_mask, cont_after & craton_after_mask),
        "continental_non_craton_node": (cont_before & ~craton_before_mask, cont_after & ~craton_after_mask),
    }
    for scope, (before_mask, after_mask) in masks.items():
        row = transitions["scopes"][scope]
        row["newly_capped_area_m2"] += float(area_after[newly_after & after_mask].sum())
        row["persistently_capped_area_m2"] += float(area_after[persistent_after & after_mask].sum())
        row["left_cap_area_m2"] += float(area_before[left_before & before_mask].sum())


def _add(scope_totals: dict, mask_before: np.ndarray, mask_after: np.ndarray, before: tuple, after: tuple) -> None:
    """Accumulate the `mask_*`-selected nodes of `before`/`after` (each `(hc, hm, area_m2)`)."""
    for suffix, mask, (hc, hm, area) in (("before", mask_before, before), ("after", mask_after, after)):
        hc, hm, area = hc[mask], hm[mask], area[mask]
        scope_totals[f"count_{suffix}"] += int(len(hc))
        scope_totals[f"sum_hc_{suffix}"] += float(hc.sum())
        scope_totals[f"sum_hm_{suffix}"] += float(hm.sum())
        scope_totals[f"area_{suffix}_m2"] += float(area.sum())
        scope_totals[f"hc_volume_{suffix}_m3"] += float(np.dot(area, hc))
        scope_totals[f"hm_volume_{suffix}_m3"] += float(np.dot(area, hm))


def record(
    world,
    plate,
    phase: str,
    hc_before: np.ndarray,
    hm_before: np.ndarray,
    codes_before: np.ndarray,
    hc_after: np.ndarray,
    hm_after: np.ndarray,
    codes_after: np.ndarray | None = None,
    *,
    area_before_m2: np.ndarray | None = None,
    area_after_m2: np.ndarray | None = None,
    craton_before_m: np.ndarray | None = None,
    craton_after_m: np.ndarray | None = None,
    node_ids_before: np.ndarray | None = None,
    node_ids_after: np.ndarray | None = None,
    plate_is_continental_before: bool | np.ndarray | None = None,
    plate_is_continental_after: bool | np.ndarray | None = None,
) -> None:
    """Record one call's before/after Hc/Hm for `phase` against `world.phase_budget`, gated by
    `world.debug_diagnostics`. `hc_before`/`hm_before`/`codes_before` and their `_after`
    counterparts describe the same node slice at the two instants -- they need not be the same
    length (a growth/shrink phase can add or remove nodes; the "before"/"after" scope totals
    are accumulated independently so this is measured correctly either way). `codes_after`
    defaults to `codes_before` for a phase that cannot itself reclassify a node.
    `area_before_m2`/`area_after_m2` are each node's accounting area; omitted, every node gets
    the nominal `lithosphere.node_area_m2`, which a cell's exact area only approximates."""
    if not getattr(world, "debug_diagnostics", False):
        return
    if len(hc_before) == 0 and len(hc_after) == 0:
        return
    if codes_after is None:
        codes_after = codes_before
    if area_before_m2 is None or area_after_m2 is None:
        nominal_m2 = lithosphere.node_area_m2(elevation_lines.line_spacing_rad(world.node_density))
        if area_before_m2 is None:
            area_before_m2 = np.full(len(hc_before), nominal_m2)
        if area_after_m2 is None:
            area_after_m2 = np.full(len(hc_after), nominal_m2)

    plate_is_continental = plate.crust_type == "continental"
    if plate_is_continental_before is None:
        plate_is_continental_before = plate_is_continental
    if plate_is_continental_after is None:
        plate_is_continental_after = plate_is_continental
    cont_before = elevation_lines.effective_is_continental_from_codes(codes_before, plate_is_continental_before)
    cont_after = elevation_lines.effective_is_continental_from_codes(codes_after, plate_is_continental_after)

    totals = world.phase_budget.setdefault(phase, _new_totals())
    totals["calls"] += 1
    scopes = totals["scopes"]

    before = (np.asarray(hc_before), np.asarray(hm_before), np.asarray(area_before_m2))
    after = (np.asarray(hc_after), np.asarray(hm_after), np.asarray(area_after_m2))
    all_before = np.ones(len(hc_before), dtype=bool)
    all_after = np.ones(len(hc_after), dtype=bool)
    _add(scopes["all"], all_before, all_after, before, after)
    _add(scopes["continental_plate" if plate_is_continental else "oceanic_plate"], all_before, all_after, before, after)
    _add(scopes["continental_node"], cont_before, cont_after, before, after)
    _add(scopes["oceanic_node"], ~cont_before, ~cont_after, before, after)
    _record_cap_transitions(
        totals,
        plate_is_continental,
        np.asarray(hm_before),
        np.asarray(hm_after),
        np.asarray(codes_before),
        np.asarray(codes_after),
        np.asarray(area_before_m2),
        np.asarray(area_after_m2),
        craton_before_m,
        craton_after_m,
        node_ids_before,
        node_ids_after,
        plate_is_continental_before,
        plate_is_continental_after,
        force_unaligned=phase in _UNALIGNED_CAP_TRANSITION_PHASES,
    )
    hm_ledger.record_change(
        world,
        plate,
        phase,
        np.asarray(hm_before),
        np.asarray(codes_before),
        np.asarray(hm_after),
        np.asarray(codes_after),
        np.asarray(area_before_m2),
        np.asarray(area_after_m2),
        node_ids_before,
        node_ids_after,
        plate_is_continental_before,
        plate_is_continental_after,
    )


def record_shortening(world, demand_m2: np.ndarray, result, hc: np.ndarray, hm: np.ndarray) -> None:
    """Book one plate's `shortening.ShorteningResult` against `SHORTENING_PHASE`, gated by
    `world.debug_diagnostics`. Shortening area `A` taken up by (or leaving the cascade at) a
    column of thickness `H` is `A * H` of volume, so `hc`/`hm` are the columns before the
    shortening."""
    if not getattr(world, "debug_diagnostics", False):
        return
    totals = world.phase_budget.setdefault(SHORTENING_PHASE, {"calls": 0, "shortening": dict.fromkeys(SHORTENING_FIELDS, 0.0)})
    totals["calls"] += 1
    book = totals["shortening"]
    returned = result.returned_decay_m2 + result.returned_unrouted_m2
    book["demand_m2"] += float(np.sum(demand_m2))
    book["absorbed_m2"] += float(np.sum(result.absorbed_m2))
    book["returned_decay_m2"] += float(np.sum(result.returned_decay_m2))
    book["returned_unrouted_m2"] += float(np.sum(result.returned_unrouted_m2))
    book["absorbed_hc_m3"] += float(np.dot(result.absorbed_m2, hc))
    book["absorbed_hm_m3"] += float(np.dot(result.absorbed_m2, hm))
    book["returned_hc_m3"] += float(np.dot(returned, hc))
    book["returned_hm_m3"] += float(np.dot(returned, hm))

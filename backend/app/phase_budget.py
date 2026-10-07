"""Per-phase Hc/Hm budget instrumentation for GitHub issue #216 ("Investigate long-run Hc/Hm
decline: geological sinks, changing node counts, and numerical volume losses").

Gated entirely by `World.debug_diagnostics` (the same flag `World.log_corner_notch` already
uses) so it costs nothing on an ordinary world. Each call site that mutates crustal thickness
(Hc) / mantle-lithosphere thickness (Hm) -- convergent/divergent deformation, arc magmatism,
oceanic cooling relaxation, decompression melting, boundary growth/shrink, row claiming,
regularization, plate merges/cleanup/relatticing, failed rifts, erosion, the end-of-step cap
clamp -- calls `record()` with the same node slice's Hc/Hm (and `crust_type_code`) just
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
(`Plate.accounting_areas_m2`: exact cells on quad plates, nominal on line plates). Quad cells
are not equal-area and a phase can swap large cells for small ones, so on a quad world only
the area-weighted fields measure real crust gained or lost (issue #257). They have to be
weighted per node as each call is recorded -- the aggregate counts can't be converted
afterwards. A caller passing bare arrays (line-level phases) gets the nominal area per node.

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

from . import elevation_lines, lithosphere

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


class Snapshot(NamedTuple):
    """One plate's columns at an instant, in node order, with each node's accounting area."""

    hc: np.ndarray
    hm: np.ndarray
    codes: np.ndarray
    area_m2: np.ndarray


def snapshot(plate, spacing_rad: float) -> Snapshot:
    """Hc, Hm, crust_type_code and `Plate.accounting_areas_m2`, each concatenated across the
    whole plate in node order (`Plate.collect`) -- for a phase that operates plate-wide rather
    than one line at a time (row claiming, corner-notch fill, merge, relattice, cleanup
    removal, quad retreat/advance, erosion, the cap clamp), where no single line's local
    `hc`/`hm` arrays capture the whole effect."""
    return Snapshot(
        plate.collect("crustal_thickness_m"),
        plate.collect("mantle_lithosphere_thickness_m"),
        plate.collect("crust_type_code"),
        plate.accounting_areas_m2(spacing_rad),
    )


def record_snapshots(world, plate, phase: str, before: Snapshot, after: Snapshot) -> None:
    """`record` for two `snapshot`s of the same plate."""
    record(
        world, plate, phase, before.hc, before.hm, before.codes, after.hc, after.hm, after.codes,
        area_before_m2=before.area_m2, area_after_m2=after.area_m2,
    )


def _new_totals() -> dict:
    return {"calls": 0, "scopes": {scope: dict.fromkeys(_SCOPE_FIELDS, 0.0) for scope in SCOPES}}


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
) -> None:
    """Record one call's before/after Hc/Hm for `phase` against `world.phase_budget`, gated by
    `world.debug_diagnostics`. `hc_before`/`hm_before`/`codes_before` and their `_after`
    counterparts describe the same node slice at the two instants -- they need not be the same
    length (a growth/shrink phase can add or remove nodes; the "before"/"after" scope totals
    are accumulated independently so this is measured correctly either way). `codes_after`
    defaults to `codes_before` for a phase that cannot itself reclassify a node.
    `area_before_m2`/`area_after_m2` are each node's accounting area; omitted, every node gets
    the nominal `lithosphere.node_area_m2` -- right only for a line plate's nodes."""
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
    cont_before = elevation_lines.effective_is_continental_from_codes(codes_before, plate_is_continental)
    cont_after = elevation_lines.effective_is_continental_from_codes(codes_after, plate_is_continental)

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

"""Per-phase Hc/Hm budget instrumentation for GitHub issue #216 ("Investigate long-run Hc/Hm
decline: geological sinks, changing node counts, and numerical volume losses").

Gated entirely by `World.debug_diagnostics` (the same flag `World.log_corner_notch` already
uses) so it costs nothing on an ordinary world. Each call site that mutates crustal thickness
(Hc) / mantle-lithosphere thickness (Hm) -- convergent/divergent deformation, arc magmatism,
oceanic cooling relaxation, decompression melting, boundary growth/shrink, row claiming,
regularization, plate merges/cleanup/relatticing, failed rifts, erosion -- calls `record()`
with the same node slice's Hc/Hm (and `crust_type_code`) just before and just after that
phase ran. `record()` accumulates the delta into `World.phase_budget[phase_name]`, broken out
by plate type (the owning plate's own `crust_type`) and by node type (`crust_type_code`,
resolved via `elevation_lines.effective_is_continental_from_codes` -- a node can carry a type
different from its plate's nominal one, e.g. an accreted terrane), so a run can separate a
genuine reservoir loss/gain in one mechanism from population change (more/fewer nodes) or
reclassification (nodes changing type without any Hc/Hm change at all).

Boundary growth/shrink (`LithospherePlate._grow_or_shrink_line_for_deform`) is itself split
into five phases -- `line_end_stretch`, `line_end_arc_grow`, `line_end_retreat`,
`line_end_accretion`, `line_interior_carve` -- rather than one lumped `line_growth_shrink`
total, since that single call site bundles several mechanistically distinct sub-events
(endpoint stretch-thinning, brand-new arc-margin nodes, plain end/interior deletion, and
accreted-column redistribution) that the issue's own investigation found dominate the whole
Hc/Hm budget yet couldn't be told apart without their own snapshots. `line_growth_shrink`
itself is still recorded too, as the call's net total -- a cross-check that the five sub-phases'
own before-after deltas sum to its delta exactly (their raw count_before/count_after and
sum_before/sum_after aren't directly comparable to the aggregate's, since each sub-phase records
only the slice it actually touches while the aggregate records the whole line on every call).

`World.phase_budget` is a plain cumulative total since the world was created (or since the
last `World.reset_phase_budget()`) -- there is no separate per-step history, by design: a
caller investigating a specific interval (e.g. a short replay from a saved world) resets the
budget, steps forward, and reads it back, rather than the engine paying to retain a full
per-step time series it may never need. "Nominal covered area" isn't stored per phase since
it's a deterministic function of node count at fixed density: multiply any `count_*` field by
`lithosphere.node_area_m2(elevation_lines.line_spacing_rad(world.node_density))`.
"""

from __future__ import annotations

import numpy as np

from . import elevation_lines

_SCOPE_FIELDS = ("count_before", "count_after", "sum_hc_before", "sum_hc_after", "sum_hm_before", "sum_hm_after")
# "all" = every node the call touched; "*_plate" splits by the owning plate's own `crust_type`
# (one bucket per call, the whole slice belongs to one plate); "*_node" splits by each node's
# own *effective* type (crust_type_code resolved against the plate), which can disagree with
# the plate bucket for e.g. an accreted terrane or a freshly-melted node stamped the other way.
SCOPES = ("all", "continental_plate", "oceanic_plate", "continental_node", "oceanic_node")


def snapshot(plate) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(hc, hm, crust_type_code), each concatenated across every one of `plate`'s own lines in
    node order (`Plate.collect`) -- for a phase that operates plate-wide rather than one line
    at a time (row claiming, corner-notch fill, merge, relattice, cleanup removal), where no
    single line's local `hc`/`hm` arrays capture the whole effect."""
    return plate.collect("crustal_thickness_m"), plate.collect("mantle_lithosphere_thickness_m"), plate.collect("crust_type_code")


def _new_totals() -> dict:
    return {"calls": 0, "scopes": {scope: dict.fromkeys(_SCOPE_FIELDS, 0.0) for scope in SCOPES}}


def _add(scope_totals: dict, count_before: int, count_after: int,
         hc_before: np.ndarray, hc_after: np.ndarray, hm_before: np.ndarray, hm_after: np.ndarray) -> None:
    scope_totals["count_before"] += count_before
    scope_totals["count_after"] += count_after
    scope_totals["sum_hc_before"] += float(hc_before.sum())
    scope_totals["sum_hc_after"] += float(hc_after.sum())
    scope_totals["sum_hm_before"] += float(hm_before.sum())
    scope_totals["sum_hm_after"] += float(hm_after.sum())


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
) -> None:
    """Record one call's before/after Hc/Hm for `phase` against `world.phase_budget`, gated by
    `world.debug_diagnostics`. `hc_before`/`hm_before`/`codes_before` and their `_after`
    counterparts describe the same node slice at the two instants -- they need not be the same
    length (a growth/shrink phase can add or remove nodes; the "before"/"after" scope totals
    are accumulated independently so this is measured correctly either way). `codes_after`
    defaults to `codes_before` for a phase that cannot itself reclassify a node."""
    if not getattr(world, "debug_diagnostics", False):
        return
    if len(hc_before) == 0 and len(hc_after) == 0:
        return
    if codes_after is None:
        codes_after = codes_before

    plate_is_continental = plate.crust_type == "continental"
    cont_before = elevation_lines.effective_is_continental_from_codes(codes_before, plate_is_continental)
    cont_after = elevation_lines.effective_is_continental_from_codes(codes_after, plate_is_continental)

    totals = world.phase_budget.setdefault(phase, _new_totals())
    totals["calls"] += 1
    scopes = totals["scopes"]

    _add(scopes["all"], len(hc_before), len(hc_after), hc_before, hc_after, hm_before, hm_after)
    _add(
        scopes["continental_plate" if plate_is_continental else "oceanic_plate"],
        len(hc_before), len(hc_after), hc_before, hc_after, hm_before, hm_after,
    )
    _add(
        scopes["continental_node"],
        int(cont_before.sum()), int(cont_after.sum()),
        hc_before[cont_before], hc_after[cont_after], hm_before[cont_before], hm_after[cont_after],
    )
    _add(
        scopes["oceanic_node"],
        int((~cont_before).sum()), int((~cont_after).sum()),
        hc_before[~cont_before], hc_after[~cont_after], hm_before[~cont_before], hm_after[~cont_after],
    )

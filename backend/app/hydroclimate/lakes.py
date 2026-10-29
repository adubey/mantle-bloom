"""An explicit n-ary tree of `Lake` objects: each lake knows its own basin (member nodes),
floor, min/max depth, and volume, replacing the flat per-node flood-fill/rate-limit mechanism
that used to live in hydrology.py (`update_lakes`/`_flood_fill_lake_extent`, and the
MAX_LAKE_NEW_HOPS_PER_MYR/MAX_LAKE_DEPTH_CHANGE_M_PER_MYR rate-limit hacks -- both removed).
That flat-array approach bounded how *fast* a lake's flood-fill could grow per step, but never
gave a lake a real, computed boundary -- so a large, gently-sloped catchment could keep creeping
outward step after step without ever reaching a physical equilibrium. `build_lake_hierarchy`
instead computes each basin's actual floor/rim geometry in one pass, so a lake's maximum extent
is known immediately rather than discovered incrementally.

**The algorithm is two phases: catchment assignment, then a Kruskal-style merge over catchment
boundaries** -- the standard "watershed by immersion" / depression-hierarchy technique (Barnes
et al.'s fill-spill-merge family), adapted from a regular grid to this codebase's k-NN graph.

Phase 1 (`_catchment_roots`) assigns every node to the catchment of the true local minimum its
own *pure* steepest descent eventually drains to (or marks it as draining straight to the ocean,
`_OCEAN_CATCHMENT`, if its chain never passes through a land local minimum first) -- a small,
separate steepest-descent pass, deliberately **not** reusing hydrology.py's own
`_compute_flow_direction` (which additionally biases toward existing channels and redirects a
full sink toward its spill target, both river-*routing* concerns orthogonal to basin
*geometry*). Without this phase, a naive Kruskal merge over every raw graph edge would wrongly
treat every intermediate ridge/rim node passed on the way down to the ocean as its own trivial
one-node "lake" (confirmed directly while developing this: a 4-node fixture with exactly one
real depression produced a spurious extra merge level wrapping the rim node itself, before this
phase was added) -- assigning nodes to catchments first means only genuine local minima ever
become `Lake` leaves, and a catchment's own internal nodes are pre-merged as one component
before the boundary-merge phase even starts.

Phase 2 walks catchment-boundary edges (an edge whose two endpoints belong to *different*
catchments, weighted by `max(elevation[i], elevation[j])` -- the saddle/col height a path
between them is forced to cross) in ascending weight order, union-finding the catchments
together -- identical in spirit to Kruskal's minimum-spanning-tree construction, except the
"tree" being built is the *merge tree* of the elevation field itself. **A lake's own `max_depth`
is just the weight of the very next merge event it participates in**, whether that unions it
with a sibling lake (creating a parent) or with the ocean (finalizing it as a root) -- both are
"the top of the basin" in the user's own words, one spilling to another basin, the other
spilling to the sea. This falls directly out of ascending-order processing: the first merge a
catchment is ever part of is, by construction, the lowest elevation at which it has *any* path
out, i.e. exactly its rim -- so a land-land merge at weight `w` sets `max_depth = w` on *both*
merging lakes and creates a new parent (`min_depth = w`, "the point where the basins split," the
user's own spec) representing the pair going forward; a land-ocean merge at weight `w` sets
`max_depth = w` on whichever lake was on the land side, with no separate Lake object needed for
"the ocean" itself.

Every genuine catchment gets a leaf `Lake` **eagerly**, before phase 2 runs, not lazily on first
merge -- a catchment that never merges with anything (a fully enclosed, single-basin world with
no other depression and no ocean anywhere reachable, i.e. a real endorheic/closed basin) would
otherwise never get a `Lake` object at all despite being exactly the "no known spill, could fill
indefinitely" case this module needs to represent (confirmed directly: a minimal 3-node
fully-enclosed fixture produced zero lakes under a lazy-materialization version of this
function).

**No stable per-node identity is used or assumed here** (matching hydrology.py's own
"recomputed fresh, not advected" philosophy) -- `build_lake_hierarchy` is a pure function of
this step's current `elevation`/`is_ocean`/`neighbor_idx`, rebuilt from scratch every step.

**Continuity across steps needs no persistent `Lake` registry at all -- see `step_lakes`.**
Every quantity a lake needs to remember (how much water it held, whether it was already merged
with a sibling) is recoverable directly from the same persistent per-node `lake_depth` array
this module already outputs and hydrology.py already threads through boundary/regrid/
reassign/merge-split exactly like `channel_depth` -- so `step_lakes` reconciles a freshly-built
tree against *that* array (`elevation + prev_lake_depth`, the same "last step's water surface"
quantity the old per-node system itself used) rather than needing to match this step's `Lake`
objects against last step's by member-overlap. This avoids an entire class of fragility the
member-overlap approach would have: a `Lake`'s exact tree shape can change step to step purely
from ordinary terrain churn (erosion, regridding) even when nothing physical about the lake
itself changed, so matching *objects* across steps would need the same kind of "is this really
the same lake" heuristic judgment call hydrology.py's own docstring already documents as a
previously-real source of bugs for the old per-node design -- deriving continuity from the
already-persisted array sidesteps needing that judgment call at all.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

_logger = logging.getLogger(__name__)

# Silt settled on a lake's bed per cubic metre of water flowing into it (issue #363): each
# step a lake gains its real inflow (`LakeClimate.inflow_m3_per_yr`, the same water that sets
# its level, including the flow through a spilling lake's sink) times the step length times
# this, as a volume spread over its wet area. Being a volume fraction, it gives a catchment
# the same silt at any node density, and the breach surge (model units, erosion-only) never
# reaches it. The value keeps the total silt laid on lake beds (after each lake's room clamp)
# under the earlier model-unit rule (1/100 of the retired lake fill rate, times
# `water_deposited`) on the seed 997271774 replay at the default node density in
# `analysis/issue363/`, rather than being set from real sediment loads: erosion.py delivers
# routed river sediment into lakes separately (`_spread_lake_sediment`). Siltation is a slow,
# permanent process next to a lake's water level, which tracks its climate every step. A
# small, low-inflow lake plausibly silts in entirely over the tens-of-Myr span a real small
# lake survives; a large lake spreads its silt over far more wet area, so it fills in far more
# slowly, as real big lakes do. Unlike the old rule, a spilling lake now silts from the water
# passing through it, as a real through-flow lake traps sediment; before, it got almost none
# apart from the breach surge.
#
# `step_lakes` returns this step's per-node silt increment (`silt_deposited`); erosion.py folds
# it straight into real terrain `elevation` (like every other deposition pathway), so a basin
# genuinely fills in and stays filled even if it later drains or is reclassified. This module
# therefore builds every hierarchy from bare `elevation` -- last step's fill is already baked
# into it -- rather than an `elevation + prev_silt_depth` effective floor.
#
# GitHub issue #144 asked whether this coefficient's magnitude actually keeps pace with the
# tectonic roughening rate that carves the small pits it's meant to fill, since it was tuned by
# feel (the ~100x above) rather than against a measured comparison. Measured directly on a live
# 30 My run of issue #117's own repro seed (23097282): the line engine's crumpling pass (the
# actual roughening injector then, run on plate node-density increases; retired with the line
# surface in #251) shifted a line's elevation by a median ~51 m per event (p90 ~160 m, up
# to ~3,000 m) against a naive resample -- vs. this coefficient (then in model units) depositing a median ~0.02
# m/step (p90 ~0.32 m/step) on a small leaf catchment's bed -- roughly a 1,000x gap per event. Despite that gap,
# re-running issue #117's exact repro (seed 23097282 @ 79.2 My) against current code via
# `app.hydroclimate.lake_hierarchy_diagnostics` shows the pathology it named already resolved: max hierarchy
# depth stays in the single digits the whole run (peak 8, ending at 2 with only 6 leaf
# catchments, none sub-6-node) -- nowhere near the originally-reported ~3,500-level cascade.
# That's because crumple events are sparse and localized (~5/step across a whole world's worth
# of lines in that run) while siltation is continuous on every wet leaf every step; #143's fix
# (silting chronically-frozen catchments instead of skipping them) was the actual bottleneck,
# not this coefficient's magnitude. No bump is justified by this measurement -- revisit only if
# a future repro reproduces runaway depth with #143 already in place.
#
# That measurement predates issue #161's fix: `rheology.apply_convergent_deformation`'s Hc/Hm
# growth was unbounded until then, desyncing a node's rendered `elevation` from what its own
# crustal/mantle-lithosphere thickness isostatically implied by up to 16 km on a save examined
# in that issue -- exactly the kind of spurious deep, narrow well this module's catchment
# detection would mistake for a real depression. Worth re-checking whether the "no bump needed"
# call above was itself an artifact of the bug it predates, not just of the frozen-catchment fix
# it credits. Re-ran the same seed (23097282) from scratch against current code (#161's fix and
# its #145/#146/#162 follow-ups included), both at the original 79.2 My checkpoint and extended
# to 300 My. At the original checkpoint itself: max depth 3 with 34 leaf catchments (vs. the
# original's peak 8/ending 2 with only 6 leaf catchments) -- a different but not worse number,
# as expected given #161's fix changes plate-deformation behaviour from early steps on, not
# just once Hc/Hm would otherwise have blown up, so this is a genuinely different simulated
# history from the same seed rather than a literal reproduction. Extended to 300 My (12
# checkpoints, every 25 My): max hierarchy depth never exceeds 17 at any checkpoint and
# oscillates rather than trending upward (4, 17, 4, 4, 7, 15, 2, 7, 7, 2, 6, 3 across the
# 25-300 My checkpoints), with every spike resolving back down within one or two checkpoints --
# the same "transient chain forms, siltation collapses it" dynamic the original measurement
# saw, not a slow climb toward issue #117's ~3,500-level cascade. Conclusion unchanged and now
# confirmed independent of #161: still no bump justified.
#
# That same extended run also answers issue #144's item #2 (whether a cheap depression
# pre-fill pass belongs in the pipeline alongside this module's catchment+Kruskal hierarchy):
# with depth staying bounded and self-resolving out to 4x the original check's duration, a
# second geometric-fill mechanism isn't warranted -- and would cut against the precedent
# `hydrology.py`'s own module docstring documents (an earlier standalone priority-flood
# basin-spill pass, `_compute_basin_spill`, was removed because it could drift out of sync with
# this module's own hierarchy). Issue #144 is closed on that basis: item #1 measured (above),
# item #2 decided (no pre-fill pass), item #3's tooling is `app.hydroclimate.lake_hierarchy_diagnostics`.
SILT_ACCUMULATION_COEFFICIENT = 4.1e-6

# Same value elevation_lines.py/mantle.py each already redefine locally rather than import --
# following that precedent here too rather than pulling in plates.py (a much heavier module,
# and one this module has never depended on) just for one scalar.
PLANET_RADIUS_KM = 6371.0

# Lake-vs-sea tiers: a real lake and a real sea behave differently (a sea drives regional
# humidity at close to ocean strength -- see climate.py's SEA_EVAPORATION_FRACTION -- and looks
# different, hydrology.HydrologyFields.is_sea/render_image.py's SEA_COLOR), and neither should
# be allowed to grow to an arbitrary size just because its own geological rim happens to sit
# high. Without any ceiling, a closed basin fills to whatever its real rim allows regardless of
# how physically implausible that is -- confirmed directly on a 286 My save (seed 931964976):
# the basin's biggest closed lake had grown to ~4.84M km^2 and 4,626 m deep, sitting right at
# its own real spill point -- 13x the Caspian Sea's actual area and 4.5x Lake Baikal's actual
# (the real deepest lake's) depth. The ceilings below are first-pass estimates pinned to those
# two real-world reference points, not swept.
#
# Comfortably above the deepest real lake (Baikal, 1,642 m) -- an ordinary closed basin's water
# column is capped here regardless of how much deeper its own geological rim would allow.
LAKE_MAX_DEPTH_M = 1800.0
# Deeper than the Caspian Sea's real 1,025 m (a simulated rift-sea basin can plausibly run
# deeper), but still finite -- unlike LAKE_MAX_DEPTH_M/SEA_MAX_DEPTH_M this doesn't reflect an
# absence of a rim; it's a second, tier-specific ceiling on top of whatever the real rim is.
SEA_MAX_DEPTH_M = 3500.0
# Between the historical Aral Sea (~68,000 km^2, before it was largely drained) and the Caspian
# (~371,000 km^2, Earth's largest lake) -- a basin has to be genuinely sea-scale, not just "a
# big lake," to promote. A physical area, not a node count, so it's resolution-invariant
# without threading `world.node_density` into this module (which takes only bare arrays, no
# `World`, by design -- see the module docstring's own testability rationale).
SEA_MIN_FLOODED_AREA_KM2 = 200_000.0

# Partial-cell flooding (issue #297, part 1). A cell's single elevation is its mean, not its
# lowest ground. Each cell's ground is therefore spread evenly over `elevation +/- relief`,
# where relief is SUBCELL_RELIEF_FRACTION of the mean height difference to its neighbours:
# halfway to each neighbour, as a piecewise-linear surface between cell centres would be.
# A lake floods a cell gradually from `elevation - relief`, so its area and volume grow
# continuously instead of jumping by whole cells. Relief is capped so a coastal cell next to
# deep ocean doesn't get kilometres of it.
SUBCELL_RELIEF_FRACTION = 0.5
SUBCELL_RELIEF_MAX_M = 500.0


def subcell_relief(elevation: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Each cell's half-range of ground elevation about its mean (see SUBCELL_RELIEF_FRACTION)."""
    if neighbor_idx.ndim != 2 or neighbor_idx.shape[1] == 0:
        return np.zeros(len(elevation))
    spread = np.abs(elevation[neighbor_idx] - elevation[:, None]).mean(axis=1)
    return np.minimum(SUBCELL_RELIEF_FRACTION * spread, SUBCELL_RELIEF_MAX_M)


# Lake evaporation (issue #297). Water held in a partly flooded cell is a shallow margin, which
# warms faster and evaporates more than open deep water. A just-wetted cell evaporates up to
# this much faster, falling to none once the cell is fully flooded. Kept below 1 so a lake's
# total evaporation still rises with its level, which keeps the balanced level unique.
SHALLOW_EVAPORATION_BOOST = 0.3


@dataclass(frozen=True)
class LakeClimate:
    """Real-unit water balance inputs, all with shape (N,).

    `inflow_m3_per_yr` is the compatibility input used by `resolve_lakes` and direct fixtures.
    The prognostic path reads disjoint settled sources from `external_supply_m3_per_yr`, falling
    back to inflow for direct callers. Optional component arrays let the step ledger keep wet
    precipitation, dry runoff and meltwater distinct. `evaporation_m_per_yr` is open-season
    evaporation at each node; optional `sublimation_m_per_yr` is the frozen-season loss, without
    the shallow-margin boost."""

    inflow_m3_per_yr: np.ndarray
    evaporation_m_per_yr: np.ndarray
    # External sources settled in each node's catchment before any inter-lake transfer.
    # When absent, `inflow_m3_per_yr` is treated as external for direct callers and fixtures.
    external_supply_m3_per_yr: np.ndarray | None = None
    wet_precip_m3_per_yr: np.ndarray | None = None
    dry_runoff_m3_per_yr: np.ndarray | None = None
    meltwater_m3_per_yr: np.ndarray | None = None
    sublimation_m_per_yr: np.ndarray | None = None


@dataclass(frozen=True)
class LakeStepLedger:
    """Integrated lake liquid-plus-ice budget for one hydrology step, all values in m³."""

    external_supply_m3: float
    evaporation_m3: float
    overflow_to_ocean_m3: float
    overflow_out_m3: float
    storage_change_m3: float
    internal_transfer_m3: float
    residual_m3: float
    wet_precip_m3: float = 0.0
    dry_runoff_m3: float = 0.0
    meltwater_m3: float = 0.0
    sublimation_m3: float = 0.0
    ice_storage_change_m3: float = 0.0
    freezing_m3: float = 0.0
    melting_m3: float = 0.0


@dataclass(frozen=True)
class LakeSolveResult:
    lake_depth: np.ndarray
    silt_deposited: np.ndarray
    lake_is_sea: np.ndarray
    events: list[LakeEvent]
    ledger: LakeStepLedger
    overflow_m3: np.ndarray
    should_spill: np.ndarray
    filled_elevation: np.ndarray
    spill_target: np.ndarray
    rim_node_idx: np.ndarray
    trace: "LakeSolveTrace | None" = None


@dataclass(frozen=True)
class LakeTransferTrace:
    """One routed spill transaction, including its source and receiving basin IDs."""

    source_node_idx: int
    target_node_idx: int
    volume_m3: float
    internal: bool
    source_lake_id: int
    receiver_lake_id: int | None


@dataclass(frozen=True)
class LakeSolveTrace:
    """Per-node source/loss terms and spill edges for one prognostic lake solve."""

    external_supply_m3: np.ndarray
    wet_precip_m3: np.ndarray
    dry_runoff_m3: np.ndarray
    meltwater_m3: np.ndarray
    evaporation_m3: np.ndarray
    liquid_sublimation_m3: np.ndarray
    transfers: tuple[LakeTransferTrace, ...]


def _active_bodies(forest, elevation, depth, hypsometry):
    """Return the current connected water bodies, using only current-step levels."""
    return governing_bodies(forest, elevation, depth, hypsometry)


def _ordered_spill_bodies(bodies, leaf_of_node, body_for_leaf, body_by_id, node_count):
    """Topologically order active bodies by hierarchy receivers, breaking sibling cycles by
    deterministic floor order so return transfers can be deferred to the next slice."""
    active_ids = {id(body) for body in bodies}
    outgoing: dict[int, int] = {}
    indegree = {id(body): 0 for body in bodies}
    for body in bodies:
        target = body.outlet_target_idx
        leaf_id = int(leaf_of_node[target]) if 0 <= target < node_count else -1
        receiver = body_for_leaf.get(leaf_id)
        if receiver is not None and id(receiver) != id(body) and id(receiver) in active_ids:
            outgoing[id(body)] = id(receiver)
            indegree[id(receiver)] += 1
    ready = sorted(
        (body for body in bodies if indegree[id(body)] == 0),
        key=lambda body: (-body.floor_elevation, body.lake_id),
    )
    ordered = []
    while ready:
        body = ready.pop(0)
        ordered.append(body)
        receiver_id = outgoing.get(id(body))
        if receiver_id is None:
            continue
        indegree[receiver_id] -= 1
        if indegree[receiver_id] == 0:
            ready.append(body_by_id[receiver_id])
            ready.sort(key=lambda item: (-item.floor_elevation, item.lake_id))
    if len(ordered) != len(bodies):
        ordered_ids = {id(body) for body in ordered}
        ordered.extend(sorted(
            (body for body in bodies if id(body) not in ordered_ids),
            key=lambda body: (-body.floor_elevation, body.lake_id),
        ))
    return ordered


def _body_leaf_mapping(bodies, leaf_of_node):
    body_for_leaf: dict[int, Lake] = {}
    body_leaf_ids: dict[int, list[int]] = {}
    for body in bodies:
        ids = sorted({int(x) for x in leaf_of_node[body.members] if x >= 0})
        body_leaf_ids[id(body)] = ids
        for leaf_id in ids:
            body_for_leaf[leaf_id] = body
    return body_for_leaf, body_leaf_ids


def _record_body_transitions(previous_bodies, current_bodies, events, seen_events):
    previous_ids = {id(body) for body in previous_bodies}
    current_ids = {id(body) for body in current_bodies}
    for body in current_bodies:
        if id(body) not in previous_ids and body.children:
            key = (id(body), "merge")
            if key not in seen_events:
                events.append(LakeEvent("merge", len(body.members), body.min_depth, len(body.children)))
                seen_events.add(key)
    for body in previous_bodies:
        if id(body) not in current_ids and body.children:
            key = (id(body), "split")
            if key not in seen_events:
                events.append(LakeEvent("split", len(body.members), body.min_depth, len(body.children)))
                seen_events.add(key)


def _export_orphaned_storage(depth, area_m2, leaf_of_node, is_ocean, accounting):
    """Route water from prior lake cells that no longer belong to a basin this step."""
    orphaned_nodes = np.flatnonzero((leaf_of_node < 0) & (depth > 0.0))
    if not len(orphaned_nodes):
        return
    orphaned_volume = depth[orphaned_nodes] * area_m2[orphaned_nodes]
    depth[orphaned_nodes] = 0.0
    accounting.overflow_out_total += float(orphaned_volume.sum())
    accounting.overflow_by_outlet[orphaned_nodes] += orphaned_volume
    if is_ocean is not None:
        accounting.ocean_overflow_total += float(
            orphaned_volume[np.asarray(is_ocean, dtype=bool)[orphaned_nodes]].sum()
        )
    if accounting.trace_transfers is not None:
        accounting.trace_transfers.extend(
            LakeTransferTrace(int(node), -1, float(volume), False, -1, None)
            for node, volume in zip(orphaned_nodes, orphaned_volume)
        )


def _empty_lake_solve_result(elevation, depth, hypsometry, is_ocean, *, export_prior_water):
    """Build an empty result, exporting prior water only when no lake bodies remain."""
    n = len(depth)
    if export_prior_water:
        overflow = depth * hypsometry.area_m2
        total_overflow = float(overflow.sum())
        ocean_overflow = (
            float(overflow[np.asarray(is_ocean, dtype=bool)].sum()) if is_ocean is not None else 0.0
        )
        final_depth = np.zeros(n)
        storage_change = -total_overflow
    else:
        overflow = np.zeros(n)
        total_overflow = ocean_overflow = storage_change = 0.0
        final_depth = depth.copy()
    zero = np.zeros(n)
    ledger = LakeStepLedger(
        0.0, 0.0, ocean_overflow, total_overflow, storage_change, 0.0, 0.0
    )
    return LakeSolveResult(
        final_depth, zero, np.zeros(n, dtype=bool), [], ledger,
        overflow, np.zeros(n, dtype=bool), np.asarray(elevation).copy(),
        np.full(n, -1, dtype=np.int64), np.full(n, -1, dtype=np.int64),
    )


@dataclass
class _PrognosticAccounting:
    evaporation_total: float
    ocean_overflow_total: float
    overflow_out_total: float
    internal_transfer_total: float
    external_used_total: float
    wet_precip_total: float
    dry_runoff_total: float
    meltwater_total: float
    integrated_supply_at_leaf: np.ndarray
    final_spilling_ids: set[int]
    overflow_by_outlet: np.ndarray
    deferred: dict[int, float]
    trace_external_by_node: np.ndarray | None = None
    trace_wet_precip_by_node: np.ndarray | None = None
    trace_dry_runoff_by_node: np.ndarray | None = None
    trace_meltwater_by_node: np.ndarray | None = None
    trace_evaporation_by_node: np.ndarray | None = None
    trace_sublimation_by_node: np.ndarray | None = None
    trace_transfers: list[LakeTransferTrace] | None = None
    # Counts are unique active bodies within this solve; gross overflow can count water again
    # when a downstream capped body receives and passes on the same transfer (#553).
    lake_cap_body_ids: set[int] = field(default_factory=set)
    sea_cap_body_ids: set[int] = field(default_factory=set)
    lake_cap_overflow_m3: float = 0.0
    sea_cap_overflow_m3: float = 0.0


def _new_prognostic_accounting(leaf_count: int, node_count: int, diagnostic_trace: bool):
    return _PrognosticAccounting(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, np.zeros(leaf_count), set(),
        np.zeros(node_count), {},
        trace_external_by_node=np.zeros(node_count) if diagnostic_trace else None,
        trace_wet_precip_by_node=np.zeros(node_count) if diagnostic_trace else None,
        trace_dry_runoff_by_node=np.zeros(node_count) if diagnostic_trace else None,
        trace_meltwater_by_node=np.zeros(node_count) if diagnostic_trace else None,
        trace_evaporation_by_node=np.zeros(node_count) if diagnostic_trace else None,
        trace_sublimation_by_node=np.zeros(node_count) if diagnostic_trace else None,
        trace_transfers=[] if diagnostic_trace else None,
    )


def _build_solve_trace(accounting, diagnostic_trace: bool):
    if not diagnostic_trace:
        return None
    return LakeSolveTrace(
        accounting.trace_external_by_node,
        accounting.trace_wet_precip_by_node,
        accounting.trace_dry_runoff_by_node,
        accounting.trace_meltwater_by_node,
        accounting.trace_evaporation_by_node,
        accounting.trace_sublimation_by_node,
        tuple(accounting.trace_transfers),
    )


def _warn_if_depth_caps_bound(accounting: _PrognosticAccounting) -> None:
    """Warn once per solve when a lake-tier or sea-tier depth guard routes overflow (#553)."""
    lake_bodies = len(accounting.lake_cap_body_ids)
    sea_bodies = len(accounting.sea_cap_body_ids)
    if lake_bodies == 0 and sea_bodies == 0:
        return
    _logger.warning(
        "lake depth-cap guard bound: lake tier (%.0f m) bodies=%d "
        "(%.3g km^3 gross overflow), sea tier (%.0f m) bodies=%d "
        "(%.3g km^3 gross overflow); overflow is booked through lake routing",
        LAKE_MAX_DEPTH_M,
        lake_bodies,
        accounting.lake_cap_overflow_m3 / 1e9,
        SEA_MAX_DEPTH_M,
        sea_bodies,
        accounting.sea_cap_overflow_m3 / 1e9,
    )


def _record_tier_cap_overflow(accounting, body, is_sea, overflow, tier_cap_limits_body):
    """Accumulate a guard binding without changing the overflow's ordinary routing (#553)."""
    if overflow <= 1e-9 or not tier_cap_limits_body:
        return
    if is_sea:
        accounting.sea_cap_body_ids.add(id(body))
        accounting.sea_cap_overflow_m3 += overflow
    else:
        accounting.lake_cap_body_ids.add(id(body))
        accounting.lake_cap_overflow_m3 += overflow


def _source_volume(rates, members, dt_step):
    return 0.0 if rates is None or dt_step <= 0.0 else float(np.maximum(rates[members], 0.0).sum() * dt_step)


def _loss_by_node(evaporated, members, level, dt_step, evaporation, sublimation, hypsometry):
    """Allocate the integrated body loss to member nodes at the solved water level."""
    area = hypsometry.area_m2[members]
    wet = hypsometry.wet_fraction(level, members)
    evaporation_rate = np.maximum(evaporation[members], 0.0)
    sublimation_rate = (
        np.zeros_like(evaporation_rate)
        if sublimation is None else np.maximum(sublimation[members], 0.0)
    )
    evaporation_loss = (
        area * wet * evaporation_rate * (1.0 + SHALLOW_EVAPORATION_BOOST * (1.0 - wet)) * dt_step
    )
    sublimation_loss = area * wet * sublimation_rate * dt_step
    raw_total = float(evaporation_loss.sum() + sublimation_loss.sum())
    if evaporated <= 0.0:
        return evaporation_loss * 0.0, sublimation_loss * 0.0
    if raw_total > 0.0:
        scale = evaporated / raw_total
        return evaporation_loss * scale, sublimation_loss * scale
    weights = area * (evaporation_rate + sublimation_rate)
    if float(weights.sum()) <= 0.0:
        weights = area
    fallback = evaporated * weights / max(float(weights.sum()), 1.0)
    total_rate = evaporation_rate + sublimation_rate
    evap_share = np.divide(
        evaporation_rate, total_rate, out=np.ones_like(total_rate), where=total_rate > 0.0
    )
    return fallback * evap_share, fallback * (1.0 - evap_share)


def _trace_transfer(accounting, source, receiver, target, volume, *, internal, depth):
    if accounting.trace_transfers is None:
        return
    source_node = int(source.outlet_node_idx)
    if source_node < 0:
        source_node = int(source.sink_node_idx)
    target_node = int(target)
    accounting.trace_transfers.append(
        LakeTransferTrace(
            int(source_node), target_node if 0 <= target_node < len(accounting.overflow_by_outlet) else -1,
            float(volume), internal, int(source.lake_id),
            None if receiver is None else int(receiver.lake_id),
        )
    )


def _trace_deferred_export(accounting, body, sink, volume, depth):
    if accounting.trace_transfers is None:
        return
    source_node = int(body.outlet_node_idx)
    if source_node < 0:
        source_node = int(body.sink_node_idx)
    accounting.trace_transfers.append(
        LakeTransferTrace(
            int(source_node), int(sink), float(volume), False, int(body.lake_id), None,
        )
    )


def _trace_node_source(accounting, rates, climate, members, dt_step):
    if accounting.trace_external_by_node is not None and dt_step > 0.0:
        accounting.trace_external_by_node[members] += rates[members] * dt_step
        accounting.trace_wet_precip_by_node[members] += np.maximum(
            climate.wet_precip_m3_per_yr[members], 0.0
        ) * dt_step
        accounting.trace_dry_runoff_by_node[members] += np.maximum(
            climate.dry_runoff_m3_per_yr[members], 0.0
        ) * dt_step
        accounting.trace_meltwater_by_node[members] += np.maximum(
            climate.meltwater_m3_per_yr[members], 0.0
        ) * dt_step


def _trace_body_losses(accounting, evaporated, members, level, dt_step, evaporation, sublimation, hypsometry):
    if accounting.trace_evaporation_by_node is None:
        return
    evaporation_loss, sublimation_loss = _loss_by_node(
        evaporated, members, level, dt_step, evaporation, sublimation, hypsometry
    )
    accounting.trace_evaporation_by_node[members] += evaporation_loss
    accounting.trace_sublimation_by_node[members] += sublimation_loss


def _integrate_body_volume(available, members, cap_level, dry_level, dt_step, evaporation, sublimation, hypsometry):
    """Implicitly integrate storage against level-dependent lake losses for one slice."""
    if len(members) == 1:
        node = int(members[0])
        area = float(hypsometry.area_m2[node])
        z = float(hypsometry.elevation[node])
        relief = float(hypsometry.relief_m[node])
        evap = float(evaporation[node])
        sub = 0.0 if sublimation is None else float(sublimation[node])
        if relief <= 0.0:
            cap_depth = max(cap_level - z, 0.0)
            cap_volume = area * cap_depth
            cap_evaporation = area * (evap + sub) * dt_step if cap_level > z else 0.0
            surplus = available - cap_volume - cap_evaporation
            if surplus > 0.0:
                return cap_volume, cap_evaporation, surplus, cap_level
            if available <= cap_evaporation:
                return 0.0, max(available, 0.0), 0.0, z
            volume = available - cap_evaporation
            return volume, cap_evaporation, 0.0, z + volume / area

        cap_wet = float(np.clip((cap_level - (z - relief)) / (2.0 * relief), 0.0, 1.0))
        if cap_level >= z + relief:
            cap_wet = 1.0
        cap_depth = max(cap_level - z, 0.0) if cap_wet >= 1.0 else relief * cap_wet * cap_wet
        cap_volume = area * cap_depth
        cap_evaporation = area * cap_wet * (
            evap * (1.0 + SHALLOW_EVAPORATION_BOOST * (1.0 - cap_wet)) + sub
        ) * dt_step
        surplus = available - cap_volume - cap_evaporation
        if surplus > 0.0:
            return cap_volume, cap_evaporation, surplus, cap_level

        q = max(available, 0.0) / area
        full_budget = relief + (evap + sub) * dt_step
        if cap_level >= z + relief and q >= full_budget:
            level = z + q - (evap + sub) * dt_step
            volume = area * max(level - z, 0.0)
            return volume, area * (evap + sub) * dt_step, 0.0, level
        alpha = relief - evap * dt_step * SHALLOW_EVAPORATION_BOOST
        beta = dt_step * (evap * (1.0 + SHALLOW_EVAPORATION_BOOST) + sub)
        discriminant = max(beta * beta + 4.0 * alpha * q, 0.0)
        denominator = beta + float(np.sqrt(discriminant))
        wet = q / beta if abs(alpha) <= 1e-15 and beta > 0.0 else (
            2.0 * q / denominator if denominator > 0.0 else 0.0
        )
        wet = min(max(wet, 0.0), cap_wet)
        level = z - relief + 2.0 * relief * wet
        volume = area * relief * wet * wet
        evaporated = area * wet * (evap * (1.0 + SHALLOW_EVAPORATION_BOOST * (1.0 - wet)) + sub) * dt_step
        return volume, min(evaporated, available), 0.0, level

    cap_volume = hypsometry.volume(cap_level, members) if np.isfinite(cap_level) else np.inf
    z = hypsometry.elevation[members]
    relief = hypsometry.relief_m[members]
    area = hypsometry.area_m2[members]
    evap_rate = evaporation[members]
    sub_rate = np.zeros_like(evap_rate) if sublimation is None else np.maximum(sublimation[members], 0.0)

    def budget(level):
        above_low = np.maximum(level - (z - relief), 0.0)
        depth = np.divide(
            above_low * above_low, 4.0 * relief, out=np.zeros_like(z), where=relief > 0.0
        )
        depth = np.where(level >= z + relief, np.maximum(level - z, 0.0), depth)
        wet = np.where(
            level >= z + relief,
            (level > z - relief).astype(float),
            np.clip(np.divide(level - (z - relief), 2.0 * relief, out=np.zeros_like(z), where=relief > 0.0), 0.0, 1.0),
        )
        volume = float(np.dot(depth, area))
        evaporation_rate = float(np.dot(area * wet, evap_rate * (1.0 + SHALLOW_EVAPORATION_BOOST * (1.0 - wet)) + sub_rate))
        wet_area = float(np.dot(area, wet))
        partial = (relief > 0.0) & (wet > 0.0) & (wet < 1.0)
        d_wet = np.divide(1.0, 2.0 * relief, out=np.zeros_like(relief), where=partial)
        evaporation_slope = float(np.dot(
            area * (evap_rate * (1.0 + SHALLOW_EVAPORATION_BOOST - 2.0 * SHALLOW_EVAPORATION_BOOST * wet) + sub_rate),
            d_wet,
        ))
        return volume, evaporation_rate * dt_step, wet_area + evaporation_slope * dt_step

    _, cap_evaporation, _ = budget(cap_level)
    surplus = available - cap_volume - cap_evaporation
    if surplus > 0.0:
        return cap_volume, cap_evaporation, surplus, cap_level
    low, high = dry_level, cap_level
    level = float(np.clip(dry_level + available / max(float(area.sum()), 1.0), low, high))
    volume, evaporated, slope = budget(level)
    for _ in range(8):
        residual = volume + evaporated - available
        if abs(residual) <= 1e-10 * max(available, 1.0):
            break
        if residual > 0.0:
            high = level
        else:
            low = level
        candidate = level - residual / slope if slope > 0.0 else np.nan
        level = candidate if low < candidate < high else 0.5 * (low + high)
        volume, evaporated, slope = budget(level)
    if abs(volume + evaporated - available) > 1e-8 * max(available, 1.0):
        for _ in range(32):
            level = 0.5 * (low + high)
            volume, evaporated, _ = budget(level)
            if volume + evaporated > available:
                high = level
            else:
                low = level
        level = 0.5 * (low + high)
        volume, evaporated, _ = budget(level)
    return volume, min(evaporated, available), 0.0, level


def _advance_prognostic_body(
    body,
    dt_step,
    rates,
    climate,
    evaporation,
    sublimation,
    depth,
    prev_lake_depth,
    hypsometry,
    node_area_km2,
    node_count,
    is_ocean,
    accounting,
    leaves,
    body_leaf_ids,
    leaf_of_node,
    body_for_leaf,
    pending,
    processed,
):
    """Advance one active body and book its evaporation, storage and spill transfers."""
    members = body.members
    accounting.final_spilling_ids.discard(id(body))
    if len(members) == 1:
        external_rate = float(rates[members[0]]) if dt_step > 0.0 else 0.0
    else:
        external_rate = float(rates[members].sum()) if dt_step > 0.0 else 0.0
    pending_volume = pending.pop(id(body), 0.0)
    input_volume = external_rate * dt_step + pending_volume
    accounting.external_used_total += external_rate * dt_step
    _trace_node_source(accounting, rates, climate, members, dt_step)
    accounting.wet_precip_total += _source_volume(climate.wet_precip_m3_per_yr, members, dt_step)
    accounting.dry_runoff_total += _source_volume(climate.dry_runoff_m3_per_yr, members, dt_step)
    accounting.meltwater_total += _source_volume(climate.meltwater_m3_per_yr, members, dt_step)
    body_leafs = body_leaf_ids[id(body)]
    if body_leafs and external_rate > 0.0 and dt_step > 0.0:
        total_rate = float(rates[members].sum())
        if total_rate > 0.0:
            for leaf_id in body_leafs:
                leaf = leaves[leaf_id]
                accounting.integrated_supply_at_leaf[leaf_id] += (
                    external_rate * dt_step * float(rates[leaf.members].sum()) / total_rate
                )
    if input_volume <= 0.0 and not np.any(depth[members] > 0.0):
        body.current_water_elevation = hypsometry.dry_level(members)
        return
    stored = (
        float(depth[members[0]] * hypsometry.area_m2[members[0]])
        if len(members) == 1
        else float((depth[members] * hypsometry.area_m2[members]).sum())
    )
    is_sea, tier_depth = _classify_tier(members, prev_lake_depth, node_area_km2)
    rim_level = body.max_depth if body.max_depth is not None else np.inf
    tier_cap_level = body.floor_elevation + tier_depth
    cap_level = min(rim_level, tier_cap_level)
    tier_cap_limits_body = not body.drains_at_floor and tier_cap_level < rim_level - 1e-6
    dry_level = hypsometry.dry_level(members)
    if body.drains_at_floor:
        cap_level = dry_level
    available = stored + input_volume
    _, evap_volume, overflow, new_level = _integrate_body_volume(
        available, members, cap_level, dry_level, dt_step, evaporation, sublimation, hypsometry
    )
    _trace_body_losses(
        accounting, evap_volume, members, new_level, dt_step, evaporation, sublimation, hypsometry
    )
    _record_tier_cap_overflow(accounting, body, is_sea, overflow, tier_cap_limits_body)
    accounting.evaporation_total += evap_volume
    if len(members) == 1:
        node = int(members[0])
        z = float(hypsometry.elevation[node])
        relief = float(hypsometry.relief_m[node])
        above_low = max(new_level - (z - relief), 0.0)
        depth[node] = (
            max(new_level - z, 0.0)
            if new_level >= z + relief
            else above_low * above_low / (4.0 * relief) if relief > 0.0 else 0.0
        )
    else:
        depth[members] = hypsometry.mean_depth(new_level, members)
    body.current_water_elevation = new_level
    body.is_spilling = overflow > 1e-9
    if body.is_spilling:
        accounting.final_spilling_ids.add(id(body))
    if overflow <= 0.0:
        return
    target = body.outlet_target_idx
    route_source = body
    traversed = {id(body)}
    while True:
        leaf_id = int(leaf_of_node[target]) if 0 <= target < node_count else -1
        receiver = body_for_leaf.get(leaf_id)
        if receiver is None or receiver is body:
            break
        if id(receiver) in traversed:
            # Depth caps can bind below a sibling saddle. In that case two capped siblings
            # cannot merge, and their outlet edges form a closed hierarchy cycle. Export this
            # cap surplus once instead of circulating it through zero-duration flush passes.
            break
        receiver_is_capped = id(receiver) in accounting.final_spilling_ids
        receiver_rim = receiver.max_depth if receiver.max_depth is not None else np.inf
        _, receiver_tier_depth = _classify_tier(receiver.members, prev_lake_depth, node_area_km2)
        receiver_cap = min(receiver_rim, receiver.floor_elevation + receiver_tier_depth)
        receiver_dry = hypsometry.dry_level(receiver.members)
        if receiver.drains_at_floor:
            receiver_cap = receiver_dry
        artificial_cap = receiver.max_depth is None or receiver_cap < receiver_rim - 1e-6 or receiver.drains_at_floor
        if id(receiver) not in processed or not receiver_is_capped or not artificial_cap:
            if id(receiver) in processed:
                accounting.deferred[leaf_id] = accounting.deferred.get(leaf_id, 0.0) + overflow
            else:
                pending[id(receiver)] = pending.get(id(receiver), 0.0) + overflow
            _trace_transfer(accounting, route_source, receiver, target, overflow, internal=True, depth=depth)
            accounting.internal_transfer_total += overflow
            if leaf_id >= 0:
                accounting.integrated_supply_at_leaf[leaf_id] += overflow
            return
        traversed.add(id(receiver))
        _trace_transfer(accounting, route_source, receiver, target, overflow, internal=True, depth=depth)
        accounting.internal_transfer_total += overflow
        if leaf_id >= 0:
            accounting.integrated_supply_at_leaf[leaf_id] += overflow
        route_source = receiver
        target = receiver.outlet_target_idx
    accounting.overflow_out_total += overflow
    _trace_transfer(accounting, route_source, None, target, overflow, internal=False, depth=depth)
    if 0 <= target < node_count:
        accounting.overflow_by_outlet[target] += overflow
        if is_ocean is not None and bool(is_ocean[target]):
            accounting.ocean_overflow_total += overflow


def solve_lakes_prognostic(
    forest: list[Lake],
    elevation: np.ndarray,
    prev_lake_depth: np.ndarray,
    years: float,
    climate: LakeClimate,
    hypsometry: CellHypsometry,
    *,
    substeps: int = 1,
    is_ocean: np.ndarray | None = None,
    diagnostic_trace: bool = False,
) -> LakeSolveResult:
    """Integrate lake storage over this timestep, promoting and splitting bodies as water
    crosses hierarchy saddles (#588).

    Each substep advances physical volume, limits evaporation to available liquid water, and
    sends only water above a body's active cap to the leaf selected by the hierarchy outlet.
    The uniform temporal partition is deterministic. An implicit volume update keeps the fast
    evaporation response stable when one slice spans thousands of years, while spill cascades
    settle through hierarchy receivers rather than repeated terrain-routing rounds.
    """
    n = len(elevation)
    depth = np.asarray(prev_lake_depth, dtype=float).copy()
    if n == 0 or years <= 0.0 or not forest:
        return _empty_lake_solve_result(
            elevation, depth, hypsometry, is_ocean,
            export_prior_water=bool(n) and years > 0.0,
        )

    for lake in iter_all_lakes(forest):
        lake.is_spilling = False

    steps = max(1, int(substeps))
    dt = years / steps
    rates = climate.external_supply_m3_per_yr
    if rates is None:
        rates = climate.inflow_m3_per_yr
    rates = np.maximum(np.asarray(rates, dtype=float), 0.0)
    evaporation = np.maximum(np.asarray(climate.evaporation_m_per_yr, dtype=float), 0.0)
    sublimation = None if climate.sublimation_m_per_yr is None else np.maximum(
        np.asarray(climate.sublimation_m_per_yr, dtype=float), 0.0
    )
    node_area_km2 = (4.0 * np.pi * PLANET_RADIUS_KM**2) / max(n, 1)
    leaves = [lake for lake in iter_all_lakes(forest) if not lake.children]
    leaf_of_node = np.full(n, -1, dtype=np.int64)
    for i, leaf in enumerate(leaves):
        leaf_of_node[leaf.members] = i

    total_initial = float((depth * hypsometry.area_m2).sum())
    accounting = _new_prognostic_accounting(len(leaves), n, diagnostic_trace)
    _export_orphaned_storage(depth, hypsometry.area_m2, leaf_of_node, is_ocean, accounting)
    events: list[LakeEvent] = []
    seen_events: set[tuple[int, str]] = set()

    for substep in range(steps + max(1, len(leaves))):
        dt_step = dt if substep < steps else 0.0
        if dt_step > 0.0:
            accounting.final_spilling_ids.clear()
        bodies = _active_bodies(forest, elevation, depth, hypsometry)
        body_by_id = {id(body): body for body in bodies}
        body_for_leaf, body_leaf_ids = _body_leaf_mapping(bodies, leaf_of_node)
        pending: dict[int, float] = {}
        for leaf_id, volume in accounting.deferred.items():
            receiver = body_for_leaf.get(leaf_id)
            if receiver is not None:
                pending[id(receiver)] = pending.get(id(receiver), 0.0) + volume
        accounting.deferred = {}
        # Hierarchy outlets select receivers; terrain flow_target is used later for exported
        # discharge only. Sibling cycles defer return flow until the next integration slice.
        ordered_bodies = _ordered_spill_bodies(bodies, leaf_of_node, body_for_leaf, body_by_id, n)
        processed: set[int] = set()
        for body in ordered_bodies:
            if dt_step <= 0.0 and pending.get(id(body), 0.0) <= 0.0:
                processed.add(id(body))
                continue
            _advance_prognostic_body(
                body, dt_step, rates, climate, evaporation, sublimation, depth, prev_lake_depth, hypsometry,
                node_area_km2, n, is_ocean, accounting, leaves, body_leaf_ids,
                leaf_of_node, body_for_leaf, pending, processed,
            )
            processed.add(id(body))

        next_bodies = _active_bodies(forest, elevation, depth, hypsometry)
        _record_body_transitions(body_by_id.values(), next_bodies, events, seen_events)
        if substep + 1 >= steps and not accounting.deferred:
            break

    # A hierarchy cycle can survive the finite zero-duration receiver flush (for example a
    # depth-capped endorheic sibling pair). Keep its final pass-on in the physical overflow
    # field and ledger rather than dropping the pending volume at the end of the solve.
    for leaf_id, volume in accounting.deferred.items():
        sink = leaves[leaf_id].sink_node_idx
        accounting.overflow_out_total += volume
        if 0 <= sink < n:
            accounting.overflow_by_outlet[sink] += volume
            if is_ocean is not None and bool(is_ocean[sink]):
                accounting.ocean_overflow_total += volume
            _trace_deferred_export(accounting, leaves[leaf_id], sink, volume, depth)
    accounting.deferred.clear()

    lake_depth = np.zeros(n)
    lake_is_sea = np.zeros(n, dtype=bool)
    silt_deposited = np.zeros(n)
    final_bodies = _active_bodies(forest, elevation, depth, hypsometry)
    should_spill = np.zeros(n, dtype=bool)
    filled_elevation = np.asarray(elevation, dtype=float).copy()
    spill_target = np.full(n, -1, dtype=np.int64)
    rim_node_idx = np.full(n, -1, dtype=np.int64)
    final_body_leaf_ids: dict[int, list[int]] = {}
    for body in final_bodies:
        final_body_leaf_ids[id(body)] = sorted({int(x) for x in leaf_of_node[body.members] if x >= 0})
        wet_members = body.members[depth[body.members] > 0.0]
        body.current_water_elevation = (
            float(np.max(hypsometry.level_from_depth(depth[wet_members], wet_members)))
            if len(wet_members) else hypsometry.dry_level(body.members)
        )
        body.is_spilling = id(body) in accounting.final_spilling_ids
        sinks = leaf_sinks(body)
        if body.max_depth is not None:
            filled_elevation[sinks] = body.max_depth
        else:
            filled_elevation[sinks] = np.inf
        spill_target[sinks] = body.outlet_target_idx
        rim_node_idx[sinks] = body.outlet_node_idx
        if body.is_spilling:
            should_spill[sinks] = True
        lake_depth[body.members] = depth[body.members]
        is_sea, _ = _classify_tier(body.members, prev_lake_depth, node_area_km2)
        wet_fraction = hypsometry.wet_fraction(body.current_water_elevation, body.members)
        wet = body.members[wet_fraction > 0.0]
        if is_sea:
            lake_is_sea[wet] = True
        if len(wet):
            supply = sum(accounting.integrated_supply_at_leaf[i] for i in final_body_leaf_ids[id(body)])
            wet_area = float((hypsometry.area_m2[wet] * wet_fraction[wet_fraction > 0.0]).sum())
            rise = SILT_ACCUMULATION_COEFFICIENT * supply / max(wet_area, 1.0)
            silt_deposited[wet] = np.minimum(rise * wet_fraction[wet_fraction > 0.0], depth[wet])

    final_storage = float((lake_depth * hypsometry.area_m2).sum())
    storage_change = final_storage - total_initial
    # External supply and internal transfers are both included at their receiving body; the
    # global residual excludes internal transfers because each is debited and credited once.
    residual = accounting.external_used_total - accounting.evaporation_total - accounting.overflow_out_total - storage_change
    ledger = LakeStepLedger(
        external_supply_m3=accounting.external_used_total,
        evaporation_m3=accounting.evaporation_total,
        overflow_to_ocean_m3=accounting.ocean_overflow_total,
        overflow_out_m3=accounting.overflow_out_total,
        storage_change_m3=storage_change,
        internal_transfer_m3=accounting.internal_transfer_total,
        residual_m3=residual,
        wet_precip_m3=accounting.wet_precip_total,
        dry_runoff_m3=accounting.dry_runoff_total,
        meltwater_m3=accounting.meltwater_total,
    )
    _warn_if_depth_caps_bound(accounting)
    trace = _build_solve_trace(accounting, diagnostic_trace)
    return LakeSolveResult(
        lake_depth, silt_deposited, lake_is_sea, events, ledger,
        accounting.overflow_by_outlet, should_spill, filled_elevation, spill_target, rim_node_idx,
        trace,
    )


def lake_evaporation_m3_per_yr(level: float, idx: np.ndarray, hyps: "CellHypsometry", climate: LakeClimate) -> float:
    """Water loss from wet cells in m³/yr (#608); only open-water evaporation gets the
    shallow-margin boost."""
    wet = hyps.wet_fraction(level, idx)
    rate = climate.evaporation_m_per_yr[idx] * (1.0 + SHALLOW_EVAPORATION_BOOST * (1.0 - wet))
    if climate.sublimation_m_per_yr is not None:
        rate = rate + climate.sublimation_m_per_yr[idx]
    return float((hyps.area_m2[idx] * wet * rate).sum())


def balanced_level(lake: "Lake", hyps: "CellHypsometry", climate: LakeClimate, low: float, high: float) -> float:
    """The level in `[low, high]` where the lake's evaporation equals its inflow, found by
    bisection (evaporation never falls as the level rises). Returns `high` when even a full
    lake can't evaporate the inflow, and `low` when there's no inflow."""
    members = lake.members
    inflow = float(climate.inflow_m3_per_yr[members].sum())
    if inflow <= 0.0:
        return low
    if lake_evaporation_m3_per_yr(high, members, hyps, climate) <= inflow:
        return high
    for _ in range(60):
        mid = 0.5 * (low + high)
        if lake_evaporation_m3_per_yr(mid, members, hyps, climate) < inflow:
            low = mid
        else:
            high = mid
        if high - low <= 1e-6:
            break
    return 0.5 * (low + high)


@dataclass(frozen=True)
class CellHypsometry:
    """Per-cell area and sub-cell relief, used to turn a lake's water level into flooded
    area, volume, and per-cell mean water depth, and back.

    A cell with mean elevation `z` and relief `r` has ground spread evenly over `[z - r, z + r]`.
    At level `h` its wet fraction is `clip((h - (z - r)) / 2r, 0, 1)`. Its mean water depth is
    `(h - (z - r))**2 / 4r` while partly flooded, and `h - z` once fully flooded. Both are
    continuous at `h = z + r`, and the mean depth's derivative is the wet fraction. With `r = 0`
    this is whole-cell flooding: depth `max(0, h - z)`.

    `lake_depth` is the per-cell mean water depth, so `lake_depth * area` is exactly the water
    the cell holds. That's the extensive reading `surface_fields` already remaps it with.
    `reference_area_m2` converts inflow into volume. It is the mean cell area, so inflow
    keeps the same scale as before, when each cell stood for an equal share of the sphere."""

    elevation: np.ndarray
    area_m2: np.ndarray
    relief_m: np.ndarray
    reference_area_m2: float

    @classmethod
    def whole_cells(cls, elevation: np.ndarray) -> "CellHypsometry":
        """Equal-area cells with no sub-cell relief: each cell floods all at once."""
        n = max(len(elevation), 1)
        area = 4.0 * np.pi * (PLANET_RADIUS_KM * 1000.0) ** 2 / n
        return cls(elevation, np.full(len(elevation), area), np.zeros(len(elevation)), area)

    @classmethod
    def from_terrain(cls, elevation: np.ndarray, area_m2: np.ndarray, neighbor_idx: np.ndarray) -> "CellHypsometry":
        area_m2 = np.asarray(area_m2, dtype=float)
        reference = float(area_m2.mean()) if len(area_m2) else 1.0
        return cls(elevation, area_m2, subcell_relief(elevation, neighbor_idx), reference)

    def dry_level(self, idx: np.ndarray) -> float:
        """The lowest ground among `idx`: the level at which they hold no water at all."""
        return float((self.elevation[idx] - self.relief_m[idx]).min())

    def mean_depth(self, level: float, idx: np.ndarray) -> np.ndarray:
        z = self.elevation[idx]
        r = self.relief_m[idx]
        above_low = np.clip(level - (z - r), 0.0, None)
        partial = np.divide(above_low * above_low, 4.0 * r, out=np.zeros_like(z), where=r > 0.0)
        return np.where(level >= z + r, np.clip(level - z, 0.0, None), partial)

    def wet_fraction(self, level: float | np.ndarray, idx: np.ndarray) -> np.ndarray:
        z = self.elevation[idx]
        r = self.relief_m[idx]
        partial = np.divide(level - (z - r), 2.0 * r, out=np.zeros_like(z), where=r > 0.0)
        return np.where(level >= z + r, (level > z - r).astype(float), np.clip(partial, 0.0, 1.0))

    def volume(self, level: float, idx: np.ndarray) -> float:
        return float((self.mean_depth(level, idx) * self.area_m2[idx]).sum())

    def wet_area(self, level: float, idx: np.ndarray) -> float:
        return float((self.wet_fraction(level, idx) * self.area_m2[idx]).sum())

    def level_from_depth(self, depth: np.ndarray, idx: np.ndarray) -> np.ndarray:
        """The level that gives each cell its mean water depth (inverse of `mean_depth`, for
        depth > 0)."""
        z = self.elevation[idx]
        r = self.relief_m[idx]
        depth = np.clip(np.asarray(depth, dtype=float), 0.0, None)
        return np.where(depth >= r, z + depth, z - r + np.sqrt(4.0 * r * depth))

    def level_for_volume(self, target_m3: float, idx: np.ndarray, low: float, high: float) -> float:
        """The level in `[low, high]` holding `target_m3` of water. Volume is increasing and
        piecewise quadratic in the level, and its derivative is the wet area, so this uses
        Newton's method, falling back to bisection whenever a step leaves the bracket."""
        level = high
        for _ in range(60):
            excess = self.volume(level, idx) - target_m3
            if abs(excess) <= 1e-9 * max(target_m3, 1.0):
                break
            if excess > 0.0:
                high = level
            else:
                low = level
            area = self.wet_area(level, idx)
            step = level - excess / area if area > 0.0 else np.nan
            level = step if low < step < high else 0.5 * (low + high)
            if high - low <= 1e-9:
                break
        return level


def _classify_tier(
    members: np.ndarray, prev_lake_depth: np.ndarray, node_area_km2: float
) -> tuple[bool, float]:
    """Whether `members` (a Lake's own catchment, leaf or merged) counts as sea tier *this*
    step, decided from *last* step's own persisted flooded extent (`prev_lake_depth`) rather
    than the depth this step is about to compute -- the same "read the already-persisted array,
    not a value not yet computed this step" idiom `_already_merged`/`_resolve`'s own `already_merged`
    check already uses, avoiding a chicken-and-egg loop between "how big is this lake" and "how
    deep is this lake allowed to get." A lake that has never yet held any water defaults to lake
    tier (promotes one step after first crossing the threshold -- negligible). Returns
    `(is_sea, tier_max_depth)` -- the caller folds `tier_max_depth` into `_water_balance`'s own
    `cap_depth`."""
    prev_wet_count = int(np.count_nonzero(prev_lake_depth[members] > 0.0))
    is_sea = prev_wet_count * node_area_km2 >= SEA_MIN_FLOODED_AREA_KM2
    return is_sea, (SEA_MAX_DEPTH_M if is_sea else LAKE_MAX_DEPTH_M)


# Sentinel `_catchment_roots` result: this node's own pure-steepest-descent chain reaches the
# ocean without ever passing through a land local minimum first -- an ordinary hillslope node,
# never part of any Lake. Distinct from a real node index (always >= 0) and from -1 (unused
# here; hydrology.py's own flow_target uses -1 for "no target", which this module deliberately
# doesn't reuse to avoid the two meanings being confused).
_OCEAN_CATCHMENT = -2

# A lake merge/split transition whose water surface sits within this band of sea level is a
# transient coastal pond on a dithering low-relief shelf (GitHub issue #122, "Speckled low-relief
# coastlines"), not a real basin event -- on a long run over a drowned coastal plain hundreds
# of these form and split every My, one pair per puddle per step, and they bury genuine
# (deep) basin/tectonic events in the UI's event console. `summarize_lake_events` collapses a
# whole step's worth of them into a single aggregate line. Deliberately generous: the dither
# is +-50 m per step around the waterline (see the TODO), but a real lake this shallow is
# still functionally a coastal pond, so aggregating it costs nothing.
NEAR_SEA_LEVEL_EVENT_BAND_M = 15.0


@dataclass(frozen=True)
class LakeEvent:
    """One lake merge or split transition detected in `_resolve` this step. `message` is the
    same human-readable line the code used to append directly; the structured fields exist so
    `summarize_lake_events` can filter/aggregate the near-sea-level coastal-pond churn before
    it reaches `World.log_event` -- see `NEAR_SEA_LEVEL_EVENT_BAND_M`. Kept as a plain value
    object (no `World` dependency) for the same testability reason the module never calls
    `world.log_event` itself."""

    kind: str  # "merge" | "split"
    node_count: int  # nodes in the merged body (merge) / the body that just split (split)
    elevation_m: float  # water surface at the transition -- the shared saddle (`min_depth`)
    basin_count: int  # basins involved: 2 for a merge, the child count for a split

    @property
    def message(self) -> str:
        if self.kind == "merge":
            return (
                f"A {self.node_count}-node lake formed as two basins merged "
                f"at elevation {self.elevation_m:.0f}m."
            )
        return (
            f"A {self.node_count}-node lake split back into {self.basin_count} "
            f"separate basins at elevation {self.elevation_m:.0f}m."
        )


def summarize_lake_events(events: list[LakeEvent], sea_level_m: float) -> list[str]:
    """Turn this step's raw `LakeEvent`s into the lines to write to the event log. Genuine
    basin events (water surface more than `NEAR_SEA_LEVEL_EVENT_BAND_M` from sea level) pass
    through individually. Near-sea-level transients are collapsed: a lone one still passes
    through verbatim, but two or more become a single aggregate count line so a dithering
    coast can't flood the log. Real events are listed first, the aggregate (if any) last."""
    near = [e for e in events if abs(e.elevation_m - sea_level_m) <= NEAR_SEA_LEVEL_EVENT_BAND_M]
    real = [e for e in events if abs(e.elevation_m - sea_level_m) > NEAR_SEA_LEVEL_EVENT_BAND_M]
    lines = [e.message for e in real]
    if len(near) == 1:
        lines.append(near[0].message)
    elif near:
        merged = sum(1 for e in near if e.kind == "merge")
        split = len(near) - merged
        lines.append(
            f"{len(near)} transient coastal ponds churned near sea level this step "
            f"({merged} merged, {split} split)."
        )
    return lines


@dataclass
class Lake:
    """One node of the depression hierarchy. A leaf (`children == []`) is a single local
    minimum's own catchment; a parent represents two or more children merged together once
    water would rise high enough to connect them. `min_depth` is meaningless for a leaf (no
    children to split apart, so it's set equal to `floor_elevation` -- see `_make_leaf`) and is
    only ever consulted by `step_lakes`'s split check when `children` is non-empty.
    `max_depth`/`outlet_node_idx` start unset (`None`/`-1`) and are filled in the moment this
    lake first participates in a merge -- see this module's own docstring for why that's exactly
    equivalent to "reaching the top of the basin." `current_water_elevation`/`is_spilling`
    are continuity-only fields: `build_lake_hierarchy` itself always initializes
    them to a dry, unspilled default (`current_water_elevation = floor_elevation`);
    `step_lakes`'s reconciliation pass is what actually carries real values forward from the
    previous step."""

    lake_id: int
    members: np.ndarray  # (M,) int64 node indices: this lake's own catchment plus every descendant's
    floor_elevation: float  # lowest bare `elevation` among members (already reflects past silt fill)
    min_depth: float  # saddle elevation where `children` split apart; == floor_elevation if no children
    max_depth: float | None  # this lake's own rim -- elevation of the next merge event; None if unresolved
    children: list["Lake"] = field(default_factory=list)
    outlet_node_idx: int = -1  # member node where max_depth's merge happens; -1 while unresolved
    current_water_elevation: float = 0.0  # this lake's water surface; == floor_elevation means dry
    is_spilling: bool = False  # hysteresis: was this lake actively spilling last step
    # Leaf-only: this catchment's own true local minimum -- the same node `_catchment_roots`
    # already computed as this leaf's `root` key, just threaded onto the object itself so
    # `compute_spill_routing` can write per-sink-node results without recomputing it. -1 on a
    # merged parent (which, being more than one original catchment, has no single sink node of
    # its own -- see compute_spill_routing's own docstring for why a parent's *subtree* of
    # leaves is what actually matters there).
    sink_node_idx: int = -1
    # The node on the FAR side of the boundary edge that set max_depth -- i.e. where this
    # lake's overflow actually escapes to, once full. Distinct from outlet_node_idx (this
    # lake's own side of that same edge, kept unchanged for main.py's outlet_xyz display):
    # outlet_target_idx is a node in whatever lies on the other side of the rim -- a sibling
    # lake's own catchment, or the ocean itself for a land-ocean merge -- not a member of this
    # lake. -1 while unresolved (mirrors outlet_node_idx's own unresolved sentinel).
    outlet_target_idx: int = -1
    # True when this lake's outlet sits at (or below) its own floor and leads to lower ground:
    # a breached pit (breaching.py) whose notch is cut to its floor, or a lake merged at its
    # floor into a basin that does. Water in the low ground inside the floor cell drains out
    # through that same outlet, so the lake holds none. Set by `_mark_floor_drains`.
    drains_at_floor: bool = False


def _mark_floor_drains(roots: list["Lake"]) -> None:
    """Sets `drains_at_floor` top-down. A root drains at its floor if it reaches the ocean
    (`max_depth` set) at its floor. A child does if it merges at its floor and its parent
    either sits lower or itself drains at its floor. Siblings merging at a shared floor inside
    a closed basin therefore don't drain, and still pool into one lake."""
    stack = [(root, None) for root in roots]
    while stack:
        lake, parent = stack.pop()
        at_floor = lake.max_depth is not None and lake.max_depth <= lake.floor_elevation + 1e-9
        if parent is None:
            lake.drains_at_floor = at_floor
        else:
            lake.drains_at_floor = at_floor and (parent.floor_elevation < lake.floor_elevation or parent.drains_at_floor)
        stack.extend((child, lake) for child in lake.children)


def _make_leaf(lake_id: int, members: list[int], elevation_list: list[float], sink_node_idx: int) -> Lake:
    members_arr = np.array(members, dtype=np.int64)
    floor = min(elevation_list[m] for m in members)
    return Lake(
        lake_id=lake_id,
        members=members_arr,
        floor_elevation=floor,
        min_depth=floor,
        max_depth=None,
        current_water_elevation=floor,
        sink_node_idx=sink_node_idx,
    )


def _catchment_roots(elevation: np.ndarray, is_ocean: np.ndarray, neighbor_idx: np.ndarray) -> list[int]:
    """Per node: the node index of the true local minimum its own pure steepest-descent chain
    drains to, or `_OCEAN_CATCHMENT` if that chain reaches the ocean first. Steepest descent
    here is deliberately the plain "lowest neighbor, full stop" rule -- see this module's own
    docstring for why hydrology.py's own `_compute_flow_direction` (channel-preference-biased,
    spill-redirecting) isn't reused. Resolved via iterative pointer-chasing with path
    compression (the same amortized technique union-find's own `find` uses) rather than
    recursion, since a chain's length is bounded only by the graph's diameter and a recursive
    version risks Python's recursion limit on a large world."""
    n = len(elevation)
    neighbor_elev = elevation[neighbor_idx]
    own_elev = elevation[:, None]
    lower_mask = neighbor_elev < own_elev
    steepest_masked = np.where(lower_mask, neighbor_elev, np.inf)
    steepest_col = np.argmin(steepest_masked, axis=1)
    rows = np.arange(n)
    has_lower = np.isfinite(steepest_masked[rows, steepest_col])
    flow_target = np.where(has_lower, neighbor_idx[rows, steepest_col], -1).astype(np.int64)
    flow_target_list = flow_target.tolist()
    is_ocean_list = is_ocean.tolist()

    catchment_root = [-1] * n
    for start in range(n):
        if catchment_root[start] != -1:
            continue
        path = []
        x = start
        while True:
            if catchment_root[x] != -1:
                break
            path.append(x)
            if is_ocean_list[x]:
                catchment_root[x] = _OCEAN_CATCHMENT
                break
            target = flow_target_list[x]
            if target < 0:
                catchment_root[x] = x  # no lower neighbor at all -- a true local minimum
                break
            x = target
        result = catchment_root[x]
        for node in path:
            catchment_root[node] = result
    return catchment_root


def build_lake_hierarchy(  # noqa: C901
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    interface_pass_elevation: np.ndarray | None = None,
) -> list["Lake"]:
    """Returns the forest of top-level `Lake` roots: every lake that was either finalized by
    reaching the ocean (`max_depth` set) or is left as an unresolved closed/endorheic basin
    (`max_depth is None`, a legitimate "no known spill" case -- a real Caspian-Sea-style basin
    that can fill to its own rim and hold indefinitely) once every catchment-boundary edge has
    been processed. Every non-root `Lake` is still reachable by walking down through `children`.
    Ocean nodes, and any node whose own steepest descent drains straight to the ocean without
    passing through a land local minimum, never appear in any `Lake.members` -- see this
    module's own docstring for the two-phase algorithm. ``interface_pass_elevation``, when
    supplied, contains caller-validated effective passes aligned with ``neighbor_idx`` and
    replaces the default endpoint-maximum edge weight. `hydrology.compute_hydrology` supplies
    one built from `breaching.py`'s notched passage elevations (established channels plus this
    step's breach notches). Without it, the endpoint maximum of bare centres is used."""
    n = len(elevation)
    if n == 0:
        return []
    k = neighbor_idx.shape[1] if neighbor_idx.ndim == 2 else 0
    if k == 0:
        return []
    if interface_pass_elevation is not None:
        interface_pass_elevation = np.asarray(interface_pass_elevation, dtype=float)
        if interface_pass_elevation.shape != neighbor_idx.shape:
            raise ValueError("interface_pass_elevation must match neighbor_idx")
        if not np.all(np.isfinite(interface_pass_elevation)):
            raise ValueError("interface_pass_elevation must be finite")

    catchment_root = _catchment_roots(elevation, is_ocean, neighbor_idx)
    catchment_members: dict[int, list[int]] = {}
    for node, root in enumerate(catchment_root):
        if root != _OCEAN_CATCHMENT:
            catchment_members.setdefault(root, []).append(node)

    ocean_idx = n
    parent = list(range(n + 1))

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    # Pre-merge each catchment's own members into one component (so its internal edges are
    # already same-root no-ops once the boundary-merge loop below runs), and separately
    # pre-merge every ocean-draining node -- true is_ocean nodes and any node whose chain
    # reaches one -- straight into one virtual OCEAN super-root.
    elevation_list = elevation.tolist()
    component_lake: dict[int, Lake] = {}
    next_id = 0
    for root, members in catchment_members.items():
        for node in members:
            if node != root:
                parent[find(node)] = find(root)
        component_lake[root] = _make_leaf(next_id, members, elevation_list, sink_node_idx=root)
        next_id += 1
    for node, root in enumerate(catchment_root):
        if root == _OCEAN_CATCHMENT:
            parent[find(node)] = find(ocean_idx)
    ocean_root = find(ocean_idx)

    roots: list[Lake] = []

    # Every same-catchment edge is already a guaranteed no-op by this point -- both its
    # endpoints were pre-merged into one component above, unconditionally, regardless of
    # this edge's own weight (see the pre-merge comment above). Only a genuine
    # catchment-boundary edge can ever trigger a real union below, so filtering to just
    # those (a single vectorized comparison against catchment_root, already computed) skips
    # the vast majority of the raw k-NN graph's edges -- typically 98%+ of them -- without
    # changing which merges happen or their order: removing a no-op edge can't affect any
    # edge that comes after it in sorted order.
    catchment_root_arr = np.asarray(catchment_root, dtype=np.int64)
    rows_all = np.repeat(np.arange(n), k)
    cols_all = neighbor_idx.ravel()
    boundary_edge = catchment_root_arr[rows_all] != catchment_root_arr[cols_all]
    rows = rows_all[boundary_edge]
    cols = cols_all[boundary_edge]
    if interface_pass_elevation is None:
        weight = np.maximum(elevation[rows], elevation[cols])
    else:
        # Index through the unfiltered directed edge array before applying the catchment-boundary
        # mask so duplicate/reciprocal graph edges retain their exact alignment.
        # An explicitly supplied pass field is trusted geometry from the caller (see
        # breaching.interface_pass_elevation).
        weight = interface_pass_elevation.ravel()[boundary_edge]
    order = np.argsort(weight, kind="stable")
    rows_list = rows[order].tolist()
    cols_list = cols[order].tolist()
    weight_list = weight[order].tolist()

    for a, b, w in zip(rows_list, cols_list, weight_list):
        ra, rb = find(a), find(b)
        if ra == rb:
            continue
        a_ocean = ra == ocean_root
        b_ocean = rb == ocean_root
        if a_ocean and b_ocean:
            continue

        if a_ocean or b_ocean:
            land_root = rb if a_ocean else ra
            land_node = b if a_ocean else a
            ocean_node = a if a_ocean else b
            lake = component_lake.pop(land_root, None)
            if lake is not None:
                lake.max_depth = w
                lake.outlet_node_idx = land_node
                lake.outlet_target_idx = ocean_node
                roots.append(lake)
            parent[land_root] = ocean_root
            continue

        # Both sides are real catchments (or already-merged parents) -- a genuine saddle.
        # Every catchment was eagerly given a leaf above, so these should always be present;
        # the fallback below is only a defensive no-op for an unreachable degenerate case
        # (e.g. a node with no neighbors at all). `ra`/`rb` are themselves real node indices
        # (union-find representatives that trace back to a literal catchment-root node -- see
        # the pre-merge loop above), so they're valid sink_node_idx values for this fallback.
        lake_a = component_lake.pop(ra, None) or _make_leaf(next_id, [ra], elevation_list, sink_node_idx=ra)
        lake_b = component_lake.pop(rb, None) or _make_leaf(next_id + 1, [rb], elevation_list, sink_node_idx=rb)
        next_id += 2

        lake_a.max_depth = w
        lake_a.outlet_node_idx = a
        lake_a.outlet_target_idx = b
        lake_b.max_depth = w
        lake_b.outlet_node_idx = b
        lake_b.outlet_target_idx = a

        parent[ra] = rb  # simple path-compressed union, same style as hydrology._lake_component_sizes
        new_root = rb
        merged = Lake(
            lake_id=next_id,
            members=np.concatenate([lake_a.members, lake_b.members]),
            floor_elevation=min(lake_a.floor_elevation, lake_b.floor_elevation),
            min_depth=w,
            max_depth=None,
            children=[lake_a, lake_b],
            current_water_elevation=min(lake_a.floor_elevation, lake_b.floor_elevation),
        )
        next_id += 1
        component_lake[new_root] = merged

    roots.extend(component_lake.values())
    _mark_floor_drains(roots)
    return roots


def iter_all_lakes(roots: list[Lake]):
    """Every `Lake` in the forest, roots and descendants alike, depth-first -- a small shared
    traversal helper since both `step_lakes` (reconciliation) and tests need to walk the whole
    tree rather than just its top-level roots."""
    stack = list(roots)
    while stack:
        lake = stack.pop()
        yield lake
        stack.extend(lake.children)


def _water_balance(
    lake: Lake,
    elevation: np.ndarray,
    years_myr: float,
    out_silt_deposited: np.ndarray,
    tier_max_depth: float,
    is_sea: bool,
    out_lake_is_sea: np.ndarray,
    climate: LakeClimate,
    hypsometry: CellHypsometry | None = None,
) -> float:
    """This lake's new water elevation, a single scalar shared by every member: its balanced
    level (`balanced_level`, issue #297), where its evaporation over its partly flooded cells
    (`CellHypsometry`) equals its real-unit inflow from `climate`, summed over every one of its
    own members, since more than one of them can be a true sink once several originally
    separate basins are merged into one lake. Clipped to `[dry level, max_depth]` (the real
    geological rim; unbounded above for a closed/endorheic basin, `max_depth is None`), and
    *also* to `tier_max_depth` (`_classify_tier`'s lake-vs-sea depth ceiling, always finite)
    regardless of what the real rim allows -- see that function's own docstring for why an
    unbounded-by-tier basin is the actual bug being fixed here. Capping a composite
    (already-merged) lake's own rise this way is what lets it fall back below its own
    `min_depth` and split on a later `_resolve` pass, using the split mechanism that already
    exists -- a depth ceiling doubles as the fix for lakes that should never have stayed merged
    this large in the first place.

    The same real inflow also drops a small volume of sediment
    (`SILT_ACCUMULATION_COEFFICIENT` per cubic metre -- see that constant's own comment) onto the
    lake bed: written into `out_silt_deposited` as *this step's* per-node increment (erosion.py
    folds it straight into real terrain `elevation`, see that constant's comment), spread only
    over the members actually underwater this step (`elevation < new_level`) so a leaf's dry
    higher catchment nodes -- which are `members` too, see the module docstring -- never get
    raised. Clamped so one step's fill can't lift the bed past the basin rim (`max_depth`).

    A cold lake keeps its liquid water and balances like any other (RFC 0014 D2(b), #539).
    Until then a lake with any member below freezing on the annual mean was forced dry, and the
    water its inflow brought was booked nowhere: a sub-zero annual mean doesn't mean the lake
    never thaws (design principle 3), and ice cover doesn't stop inflow, storage or spill.
    Silt was still deposited on such lakes (#117), and every lake now deposits it the same way."""
    hyps = hypsometry if hypsometry is not None else CellHypsometry.whole_cells(elevation)
    members = lake.members
    # The same real inflow that sets the level sizes the silt (issue #363).
    inflow_m3_per_yr = float(climate.inflow_m3_per_yr[members].sum())

    rim_cap_depth = (lake.max_depth - lake.floor_elevation) if lake.max_depth is not None else None
    cap_depth = tier_max_depth if rim_cap_depth is None else min(rim_cap_depth, tier_max_depth)
    dry_level = hyps.dry_level(members)
    # A lake whose outlet is cut to its floor and leads lower holds no water at all, not even
    # in the low ground inside its floor cell (see Lake.drains_at_floor).
    cap_level = dry_level if lake.drains_at_floor else max(lake.floor_elevation + cap_depth, dry_level)
    new_level = balanced_level(lake, hyps, climate, dry_level, cap_level)

    wet_fraction = hyps.wet_fraction(new_level, members)
    is_wet = wet_fraction > 0.0
    wet = members[is_wet]
    if is_sea and len(wet) > 0:
        out_lake_is_sea[wet] = True
    if len(wet) > 0:
        # The silt volume settles evenly over the wet area, so a partly flooded cell gains its
        # wet fraction's share of that thickness, as a cell-mean rise.
        wet_area = float((wet_fraction[is_wet] * hyps.area_m2[wet]).sum())
        silt_m3 = SILT_ACCUMULATION_COEFFICIENT * inflow_m3_per_yr * years_myr * 1e6
        silt_rise = silt_m3 / wet_area * wet_fraction[is_wet]
        # Never lift a node's bed above this step's own water surface (sediment settles under
        # water, it doesn't pile into the air) -- which also keeps the fill at or below the
        # basin rim, since `new_level` is already clipped to `max_depth`. A cell's room is the
        # water it holds, as a mean depth.
        room = hyps.mean_depth(new_level, wet)
        out_silt_deposited[wet] += np.minimum(silt_rise, room)
    return new_level


def _already_merged(
    lake: Lake,
    elevation: np.ndarray,
    lake_depth: np.ndarray,
    hypsometry: CellHypsometry | None = None,
) -> bool:
    """Whether water in every child's main basin reached this parent's saddle.

    `members` also contains perched sub-basins in the upstream catchments. Their water level
    says nothing about whether the main bodies on both sides of this saddle are connected, so
    inspect each child's floor (its sink for a leaf) instead of the highest wet member. The
    floor cell's mean depth is inverted through the hypsometry to recover its water surface.
    At zero depth, use the cell's dry level so equal-floor basins that merge at their shared
    floor remain connected across steps (#520).

    When checking persisted water, all children must have reached the saddle. New merges during
    this step use the children’s just-resolved levels directly in `_MERGE_CHECK`.
    """
    if not lake.children:
        return True
    hyps = hypsometry if hypsometry is not None else CellHypsometry.whole_cells(elevation)
    reached_saddle: list[bool] = []
    for child in lake.children:
        floor_idx = child.sink_node_idx
        if floor_idx < 0:
            floor_idx = int(child.members[np.argmin(elevation[child.members])])
        floor_cell = np.array([floor_idx])
        depth = float(lake_depth[floor_idx])
        if depth <= 0.0:
            level = hyps.dry_level(floor_cell)
        else:
            level = float(hyps.level_from_depth(np.array([depth]), floor_cell)[0])
        reached_saddle.append(level >= lake.min_depth - 1e-6)
    return all(reached_saddle)


def compute_spill_routing(
    forest: list[Lake], elevation: np.ndarray, prev_lake_depth: np.ndarray, hypsometry: CellHypsometry | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-node `(filled_elevation, spill_target, rim_node_idx)`, the triple
    `hydrology._compute_flow_direction`/`hydrology.compute_hydrology` consult at a sink node to
    decide whether it should redirect, and where: `filled_elevation` the elevation a sink's
    water surface must reach before it escapes its basin, `spill_target` the one-hop neighbor to
    escape toward once it does, `rim_node_idx` the governing basin's own `outlet_node_idx` --
    the member on *this* lake's own side of that same boundary edge, as opposed to
    `spill_target`/`outlet_target_idx` which sits on the far side, outside the lake entirely.

    **Why callers need both sides of one edge, not just `spill_target`.** The boundary edge's
    weight (`max_depth`) is `max(elevation[outlet_node_idx], elevation[outlet_target_idx])` --
    whichever of the *two* is higher is the actual dam a lake's overflow has to cut through, and
    it's roughly a coin flip which side that turns out to be (confirmed directly: a small
    synthetic 6-node basin with the lake's own rim node higher than the far side reproduces
    exactly this). `spill_target` alone only ever exposes the far side, so a breach-erosion term
    keyed only off it silently never touches the dam on the *near* side whenever that's the
    higher one -- a real lake behind a rim like that would spill and erode its outflow channel
    forever downstream of the rim while the rim itself never wears down, so the lake could never
    actually drain. Exposing `rim_node_idx` alongside `spill_target` lets a caller erode *both*
    members of the governing edge, so the true bottleneck erodes regardless of which side it's
    on -- see hydrology.compute_hydrology's own `rim_breach_source` for where this is used.

    Formerly computed by an entirely separate,
    independent priority-flood over bare terrain (`hydrology._basin_spill_kernel`) that never
    agreed exactly with *this* module's own catchment+Kruskal hierarchy -- the two usually
    landed close but were evaluated differently and could drift apart step to step (confirmed
    directly against a real save: a merged lake's own fresh `_resolve` showed it hadn't yet
    reached its true rim, while the old independent kernel's stale, one-step-lagged threshold
    said it had already spilled). Deriving both arrays from this module's *own* tree instead
    makes `build_lake_hierarchy`'s hierarchy the single source of truth for where a basin's
    overflow goes -- see this module's own docstring for why that hierarchy, not a separate
    Dijkstra pass, is the physically-correct one.

    Every node defaults to its own bare `elevation`/`-1` (an ordinary hillslope node has no
    real rim to speak of -- any downhill neighbor already clears it via plain steepest descent,
    the same "no real rim ever blocks it" fact `_compute_flow_direction`'s own docstring relies
    on) and only a genuine sink node -- a lakes.py catchment-root leaf's own local minimum,
    `Lake.sink_node_idx` -- is ever overwritten below; see this module's own docstring for why
    a `_catchment_roots` catchment root is always exactly a `_compute_flow_direction` sink node.
    A leaf whose chain never reaches the ocean (`max_depth is None`, a genuine closed/endorheic
    basin) gets `filled_elevation = inf`/`spill_target = -1`/`rim_node_idx = -1`, matching the
    old kernel's own "unreachable" convention -- it never spills, however full it gets.

    **Which rim currently governs a sink is decided the same way `_resolve` decides "is this
    lake one connected body yet," and deliberately reuses that exact test (`_already_merged`
    against
    `min_depth`)** -- this is the one place this function depends on `prev_lake_depth`, not just
    bare terrain, and it's what fixes the multi-sink-under-one-merged-body drift: once several
    original leaf catchments have merged (as of *last* step's persisted water) into one shared
    body, every one of their original sink nodes -- each still a literal local minimum in the
    bare-terrain flow graph even while fully submerged -- must escape via the *merged* body's
    own current rim, not whatever shallower rim its own individual leaf used to have before it
    was swallowed. Walking down from each root, the first node found to already be one merged
    body (a leaf, trivially, or a parent whose children's main water bodies both reached the
    saddle last step) is
    exactly that governing rim: everything in its subtree writes that same `(max_depth,
    outlet_target_idx, outlet_node_idx)` triple. A node still short of that (a parent whose children were *not*
    yet one body last step) recurses into its children independently instead, since they may
    each still be governed by their own, still-separate, shallower rims. Because the walk always
    stops at the *first* such node encountered from the root down, and `_already_merged` checks
    each child's main water body (so a perched lake cannot mark the parent as merged), no leaf's
    `sink_node_idx` can ever be written by more than one stopping point -- the whole walk does
    `O(total tree size)` work, the same complexity class as `_resolve` itself, not `O(depth x
    leaf count)`.

    Implemented as an explicit-stack walk, matching `_catchment_roots`/`build_lake_hierarchy`/
    `_resolve`'s own iterative style for the same reason: a long spill cascade's tree can run
    thousands of levels deep, well past Python's recursion limit."""
    n = len(elevation)
    filled_elevation = elevation.copy()
    spill_target = np.full(n, -1, dtype=np.int64)
    rim_node_idx = np.full(n, -1, dtype=np.int64)
    for body in governing_bodies(forest, elevation, prev_lake_depth, hypsometry):
        sinks = leaf_sinks(body)
        filled_elevation[sinks] = body.max_depth if body.max_depth is not None else np.inf
        spill_target[sinks] = body.outlet_target_idx
        rim_node_idx[sinks] = body.outlet_node_idx
    return filled_elevation, spill_target, rim_node_idx


def governing_bodies(
    forest: list[Lake], elevation: np.ndarray, prev_lake_depth: np.ndarray, hypsometry: CellHypsometry | None = None
) -> list[Lake]:
    """The lakes whose rims currently govern their leaves' sinks, walking down from each root to
    the first lake whose children's main water bodies reached its saddle last step (see
    `compute_spill_routing`). Every leaf lies under exactly one of them."""
    bodies: list[Lake] = []
    stack: list[Lake] = list(forest)
    while stack:
        node = stack.pop()
        already_merged = _already_merged(node, elevation, prev_lake_depth, hypsometry)
        if node.children and not already_merged:
            stack.extend(node.children)
            continue
        bodies.append(node)
    return bodies


def leaf_sinks(lake: Lake) -> np.ndarray:
    """The `sink_node_idx` of every leaf under `lake` (itself, for a leaf)."""
    sinks: list[int] = []
    stack = [lake]
    while stack:
        cur = stack.pop()
        if cur.children:
            stack.extend(cur.children)
        else:
            sinks.append(cur.sink_node_idx)
    return np.asarray(sinks, dtype=np.int64)


# `_resolve`'s stack entry modes: resolve a lake, or finish a parent once its children are done.
_DESCEND, _MERGE_CHECK, _SPLIT_DONE = range(3)


def _resolve(
    lake: Lake,
    elevation: np.ndarray,
    prev_lake_depth: np.ndarray,
    years_myr: float,
    out_lake_depth: np.ndarray,
    out_silt_deposited: np.ndarray,
    out_lake_is_sea: np.ndarray,
    node_area_km2: float,
    events: list[LakeEvent],
    climate: LakeClimate,
    hypsometry: CellHypsometry | None = None,
) -> float:
    """Resolves one lake (and, for a parent, implicitly its whole subtree) into this step's
    actual water elevation, writing every member's own `out_lake_depth` exactly once, and
    detects a merge/split transition purely from whether last step's recorded depth
    (`prev_lake_depth`) showed both children's main water bodies reaching their shared
    `min_depth` -- see `_already_merged` for the connectivity check.
    per-node array rather than a persistent `Lake` registry.

    A leaf, or a parent whose children's main water bodies **already reached** the saddle last
    step, is resolved as a single connected water body over its own full `members` set.
    If that resolution's *new* level drops back below `min_depth`, it just split this step: the
    merged resolution is discarded (its silt and sea flags restored from a snapshot) and each
    child is resolved at its own balance this same step, exactly as a never-merged pair would
    be. Children used to be written full to the saddle instead (continuity), but nothing
    supplied that water, and the saddle level read back as "still merged" next step, so a lake
    whose budget couldn't hold it there sat on its saddle indefinitely (#519, #521).

    A parent that **was not** yet merged last step resolves its children independently first;
    if either child's own new level has now reached this parent's `min_depth`, they merge this
    step: continuity again pins the newly-merged body at exactly `min_depth` (any growth beyond
    the saddle in this same step is deferred to next step's now-pooled water balance) rather
    than also folding in that step's remaining growth, which would need choosing an
    otherwise-arbitrary allocation between "how much of this step's inflow arrived before vs.
    after the moment of merging."

    Returns this lake's own resolved level (informational -- only consulted by a parent's own
    merge check above; a lake that just split returns its higher child's level, since the two
    children no longer share one level).

    Implemented as an explicit-stack post-order walk of the subtree rather than plain
    recursion: the merge hierarchy `build_lake_hierarchy` produces is a near-linear chain for
    a long spill cascade (thousands of levels deep on an old, rough world), and a recursive
    resolve -- one Python frame per level -- overran the interpreter's recursion limit.
    `build_lake_hierarchy` and `_catchment_roots` are already iterative for exactly this
    reason. A parent's entry can come back after its children are done: `_MERGE_CHECK` for a
    pair that was still separate, `_SPLIT_DONE` for a merged body that just split. A resolved
    lake's level is read back off its own `current_water_elevation`.

    `out_lake_is_sea`/`node_area_km2` feed `_classify_tier` at every point this function
    writes `out_lake_depth` -- each write site classifies the *members it's actually writing*
    (a just-split child's own extent, not its former parent's; a just-merged body's own now-
    larger extent, not either child's alone), so `_classify_tier`'s own `SEA_MIN_FLOODED_AREA_
    KM2` reads on the same set of members a viewer would actually see as one lake this step."""
    hyps = hypsometry if hypsometry is not None else CellHypsometry.whole_cells(elevation)
    stack: list[tuple[Lake, int]] = [(lake, _DESCEND)]
    while stack:
        node, mode = stack.pop()

        if mode == _SPLIT_DONE:
            node.current_water_elevation = max(child.current_water_elevation for child in node.children)
            node.is_spilling = False
            continue

        if mode == _MERGE_CHECK:
            # A parent that was still split last step: its children are now resolved, so
            # decide whether this step's rise has merged them.
            child_levels = [child.current_water_elevation for child in node.children]
            if max(child_levels) >= node.min_depth:
                events.append(LakeEvent(kind="merge", node_count=len(node.members), elevation_m=node.min_depth, basin_count=2))
                is_sea, tier_max_depth = _classify_tier(node.members, prev_lake_depth, node_area_km2)
                # A child capable of independently reaching min_depth on its own (this can only
                # happen if that child's *own* tier cap is above the saddle -- e.g. a child big
                # enough to classify sea tier growing past an internal saddle lower than
                # SEA_MAX_DEPTH_M) legitimately merges -- ordinary lake physics, unrelated to
                # this section's own cap. But the *merged whole*'s own tier cap still applies to
                # the result: capping here at merge time (not just in _water_balance, which
                # only ever runs on an already-single body) is what stops the merged level from
                # sitting above the ceiling this step's caterpillar-splitting fix would then
                # immediately undo again -- min_depth alone is a real rim value that predates
                # any tier cap, so it can land above tier_max_depth just as easily as below it.
                merge_level = min(node.min_depth, node.floor_elevation + tier_max_depth)
                out_lake_depth[node.members] = hyps.mean_depth(merge_level, node.members)
                if is_sea:
                    out_lake_is_sea[node.members[elevation[node.members] < merge_level]] = True
                node.current_water_elevation = merge_level
                node.is_spilling = node.max_depth is not None and merge_level >= node.max_depth
            else:
                node.current_water_elevation = max(child_levels)
            continue

        already_merged = _already_merged(node, elevation, prev_lake_depth, hyps)

        if node.children and not already_merged:
            # Resolve children independently first, then merge-check on the way back up.
            # Pushed in reverse so the leftmost child resolves first, matching the old
            # list-comprehension order (only affects `events` ordering -- sibling subtrees
            # write disjoint node sets).
            stack.append((node, _MERGE_CHECK))
            stack.extend((child, _DESCEND) for child in reversed(node.children))
            continue

        # A leaf, or a parent already one merged body last step: single connected water body.
        if node.children:
            # If this body splits, its merged resolution is discarded, silt and sea flags too.
            saved_silt = out_silt_deposited[node.members]
            saved_sea = out_lake_is_sea[node.members]
        is_sea, tier_max_depth = _classify_tier(node.members, prev_lake_depth, node_area_km2)
        new_level = _water_balance(
            node, elevation, years_myr,
            out_silt_deposited, tier_max_depth, is_sea, out_lake_is_sea, climate, hyps,
        )
        if node.children and new_level < node.min_depth:
            events.append(LakeEvent(kind="split", node_count=len(node.members), elevation_m=node.min_depth, basin_count=len(node.children)))
            # Each child is balanced on its own water this step, under its own rim and tier
            # cap, rather than filled to the saddle (#521). Neither can rise past this saddle
            # (it is each child's own rim). A child fed enough to reach it reads as merged
            # next step and splits again: per-step churn, but no water is created (#520 part 1
            # is the false merge behind it).
            out_silt_deposited[node.members] = saved_silt
            out_lake_is_sea[node.members] = saved_sea
            stack.append((node, _SPLIT_DONE))
            stack.extend((child, _DESCEND) for child in reversed(node.children))
        else:
            out_lake_depth[node.members] = hyps.mean_depth(new_level, node.members)
            node.current_water_elevation = new_level
            node.is_spilling = node.max_depth is not None and new_level >= node.max_depth

    return lake.current_water_elevation


def resolve_lakes(
    forest: list[Lake],
    elevation: np.ndarray,
    prev_lake_depth: np.ndarray,
    years: float,
    climate: LakeClimate,
    hypsometry: CellHypsometry | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[LakeEvent]]:
    """The back half of `step_lakes` -- everything after `build_lake_hierarchy` -- pulled out
    on its own so `compute_hydrology` can build the hierarchy once, early (before flow routing,
    which now needs it too -- see `compute_spill_routing`), and pass that *same* `forest` in
    here afterward rather than paying for a second, redundant `build_lake_hierarchy` call every
    step. Resolves every top-level lake's own water balance (growth/evaporation/merge/split,
    and this step's own silt drop -- see `_resolve`/`_water_balance`). Returns `(lake_depth,
    silt_deposited, lake_is_sea, events)` -- see `step_lakes`'s own docstring for what each
    means; `forest` itself isn't returned here since the caller already has the exact object it
    passed in, now mutated in place with this step's resolved `current_water_elevation`/
    `is_spilling`.

    `node_area_km2` -- the physical area one node stands in for, `_classify_tier`'s own lake-
    vs-sea size threshold -- is derived once here from `len(elevation)` (a sphere's fixed
    surface area spread evenly over however many nodes this call sees), not passed in: every
    caller already has `elevation` at hand, and this avoids threading `world.node_density`
    into a module that otherwise takes only bare arrays -- see `SEA_MIN_FLOODED_AREA_KM2`'s own
    comment.

    `hypsometry` gives each cell's real area and sub-cell relief, for partial-cell flooding
    (see `CellHypsometry`). Without it every cell is an equal share of the sphere and floods
    all at once. `climate` resolves each lake at its real-unit balanced level (see
    `LakeClimate`), and its inflow sizes this step's silt."""
    n = len(elevation)
    lake_depth = np.zeros(n)
    silt_deposited = np.zeros(n)
    lake_is_sea = np.zeros(n, dtype=bool)
    if n == 0:
        return lake_depth, silt_deposited, lake_is_sea, []

    node_area_km2 = (4.0 * np.pi * PLANET_RADIUS_KM**2) / n
    years_myr = years / 1_000_000.0
    events: list[LakeEvent] = []
    for root in forest:
        _resolve(
            root, elevation, prev_lake_depth, years_myr,
            lake_depth, silt_deposited, lake_is_sea, node_area_km2, events, climate, hypsometry,
        )
    return lake_depth, silt_deposited, lake_is_sea, events


def step_lakes(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    prev_lake_depth: np.ndarray,
    years: float,
    climate: LakeClimate,
    hypsometry: CellHypsometry | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Lake], list[LakeEvent]]:
    """The full per-step lake update, for a caller that doesn't need to reuse the depression
    hierarchy for anything else (`compute_hydrology` itself no longer calls this directly --
    see `resolve_lakes`'s own docstring -- but every existing test, and any other future
    caller that just wants "lakes for this terrain, resolved," still can): rebuild this step's
    depression hierarchy from current terrain (`build_lake_hierarchy`, bare `elevation` --
    everything silted in through the end of last step is already folded into real `elevation`
    by erosion.py, so there's no separate effective-floor offset to add), then `resolve_lakes`
    it. Returns `(lake_depth, silt_deposited, lake_is_sea, forest, events)`: `lake_depth` an
    (N,) array aligned with `elevation` (0 wherever no lake reaches, exactly the same
    shape/meaning hydrology.py's old `update_lakes` produced, for every other consumer --
    render_image.py, coastline.py, erosion.py, group_rivers -- to keep reading unchanged);
    `silt_deposited` *this step's* per-node sediment increment for erosion.py to add straight
    into `elevation` (always >= 0, nonzero only under standing water this step); `lake_is_sea`
    an (N,) bool array, true wherever `lake_depth` is nonzero *and* that water body currently
    classifies as sea tier (`_classify_tier`) -- climate.py's own moisture source and
    render_image.py's overlay both key off it the same way they already key off `lake_depth`;
    `forest` this step's fully-resolved `Lake` tree (each lake's own `current_water_elevation`/
    `is_spilling` set to its outcome this step -- informational only, e.g. for a future
    per-lake stats view; nothing in this module reads it back next step, see the module
    docstring for why no persistent registry is needed); and `events` this step's `LakeEvent`s
    (merge/split transitions -- the caller runs them through `summarize_lake_events` and logs
    the result, mirroring merge_split.apply_topology_changes's own pattern; this module never
    calls `world.log_event` directly, keeping it testable without a `World`)."""
    forest = build_lake_hierarchy(elevation, is_ocean, neighbor_idx)
    lake_depth, silt_deposited, lake_is_sea, events = resolve_lakes(
        forest, elevation, prev_lake_depth, years, climate, hypsometry
    )
    return lake_depth, silt_deposited, lake_is_sea, forest, events

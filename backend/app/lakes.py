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

from dataclasses import dataclass, field

import numpy as np

# Lake growth/evaporation -- moved here from hydrology.py (formerly `update_lakes`'s own
# constants) unchanged in value and meaning; see hydrology.py's module docstring for the
# original tuning rationale (evaporation rescaled from a much finer-substep-suited rate down to
# this codebase's own actual step granularity). `LAKE_FILL_RATE` scales how much of a step's
# routed inflow becomes standing water; growth is spread over a lake's own member count as an
# area proxy (the same implicit "node count as area" convention `_lake_component_sizes` and
# hydrology.LAKE_BREACH_EROSION_COEFFICIENT already rely on) rather than piling the whole
# inflow onto one point, since a lake's water level is now tracked as a single scalar shared by
# every member, not a per-node depth.
LAKE_FILL_RATE = 0.12
LAKE_EVAPORATION_RATE_PER_MYR = 0.3
LAKE_EVAPORATION_BASELINE_M_PER_MYR = 0.5

# How much of a lake's own inflow settles as silt on its bed each step, in the same "spread
# over member count" units LAKE_FILL_RATE uses for water depth -- deliberately ~100x smaller:
# siltation is a slow, permanent process (silt_depth never erodes back away, see
# plates.ElevationLine's own comment) compared to a lake's water level, which rises and falls
# with ordinary season-to-season-equivalent inflow/evaporation swings every step. At this rate a
# small, low-inflow lake plausibly silts in entirely over the tens-of-Myr span a real small lake
# actually survives; a large, actively-flowing lake's own inflow is spread over far more
# members, so its own silt accumulates correspondingly slower per node, matching the real-world
# intuition that big lakes take far longer to fill in than small ones.
SILT_ACCUMULATION_COEFFICIENT = LAKE_FILL_RATE / 100.0

# Sentinel `_catchment_roots` result: this node's own pure-steepest-descent chain reaches the
# ocean without ever passing through a land local minimum first -- an ordinary hillslope node,
# never part of any Lake. Distinct from a real node index (always >= 0) and from -1 (unused
# here; hydrology.py's own flow_target uses -1 for "no target", which this module deliberately
# doesn't reuse to avoid the two meanings being confused).
_OCEAN_CATCHMENT = -2


@dataclass
class Lake:
    """One node of the depression hierarchy. A leaf (`children == []`) is a single local
    minimum's own catchment; a parent represents two or more children merged together once
    water would rise high enough to connect them. `min_depth` is meaningless for a leaf (no
    children to split apart, so it's set equal to `floor_elevation` -- see `_make_leaf`) and is
    only ever consulted by `step_lakes`'s split check when `children` is non-empty.
    `max_depth`/`outlet_node_idx` start unset (`None`/`-1`) and are filled in the moment this
    lake first participates in a merge -- see this module's own docstring for why that's exactly
    equivalent to "reaching the top of the basin." `current_water_elevation`/`is_spilling`/
    `silt_depth` are continuity-only fields: `build_lake_hierarchy` itself always initializes
    them to a dry, unspilled, silt-free default (`current_water_elevation = floor_elevation`);
    `step_lakes`'s reconciliation pass is what actually carries real values forward from the
    previous step."""

    lake_id: int
    members: np.ndarray  # (M,) int64 node indices: this lake's own catchment plus every descendant's
    floor_elevation: float  # lowest elevation among members (silt-adjusted once accumulated)
    min_depth: float  # saddle elevation where `children` split apart; == floor_elevation if no children
    max_depth: float | None  # this lake's own rim -- elevation of the next merge event; None if unresolved
    children: list["Lake"] = field(default_factory=list)
    outlet_node_idx: int = -1  # member node where max_depth's merge happens; -1 while unresolved
    current_water_elevation: float = 0.0  # this lake's water surface; == floor_elevation means dry
    is_spilling: bool = False  # hysteresis: was this lake actively spilling last step
    silt_depth: float = 0.0  # accumulated silt depth at this lake's own floor


def _make_leaf(lake_id: int, members: list[int], elevation_list: list[float]) -> Lake:
    members_arr = np.array(members, dtype=np.int64)
    floor = min(elevation_list[m] for m in members)
    return Lake(lake_id=lake_id, members=members_arr, floor_elevation=floor, min_depth=floor, max_depth=None, current_water_elevation=floor)


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


def build_lake_hierarchy(elevation: np.ndarray, is_ocean: np.ndarray, neighbor_idx: np.ndarray) -> list["Lake"]:
    """Returns the forest of top-level `Lake` roots: every lake that was either finalized by
    reaching the ocean (`max_depth` set) or is left as an unresolved closed/endorheic basin
    (`max_depth is None`, a legitimate "no known spill" case -- a real Caspian-Sea-style basin
    that can fill to its own rim and hold indefinitely) once every catchment-boundary edge has
    been processed. Every non-root `Lake` is still reachable by walking down through `children`.
    Ocean nodes, and any node whose own steepest descent drains straight to the ocean without
    passing through a land local minimum, never appear in any `Lake.members` -- see this
    module's own docstring for the two-phase algorithm."""
    n = len(elevation)
    if n == 0:
        return []
    k = neighbor_idx.shape[1] if neighbor_idx.ndim == 2 else 0
    if k == 0:
        return []

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
        component_lake[root] = _make_leaf(next_id, members, elevation_list)
        next_id += 1
    for node, root in enumerate(catchment_root):
        if root == _OCEAN_CATCHMENT:
            parent[find(node)] = find(ocean_idx)
    ocean_root = find(ocean_idx)

    roots: list[Lake] = []

    rows = np.repeat(np.arange(n), k)
    cols = neighbor_idx.ravel()
    weight = np.maximum(elevation[rows], elevation[cols])
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
            lake = component_lake.pop(land_root, None)
            if lake is not None:
                lake.max_depth = w
                lake.outlet_node_idx = land_node
                roots.append(lake)
            parent[land_root] = ocean_root
            continue

        # Both sides are real catchments (or already-merged parents) -- a genuine saddle.
        # Every catchment was eagerly given a leaf above, so these should always be present;
        # the fallback below is only a defensive no-op for an unreachable degenerate case
        # (e.g. a node with no neighbors at all).
        lake_a = component_lake.pop(ra, None) or _make_leaf(next_id, [ra], elevation_list)
        lake_b = component_lake.pop(rb, None) or _make_leaf(next_id + 1, [rb], elevation_list)
        next_id += 2

        lake_a.max_depth = w
        lake_a.outlet_node_idx = a
        lake_b.max_depth = w
        lake_b.outlet_node_idx = b

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
    lake: Lake, prev_level: float, water_deposited: np.ndarray, years_myr: float, is_freezing: bool, out_silt_depth: np.ndarray
) -> float:
    """This lake's new water elevation, generalizing the old per-node `update_lakes` formula
    (hydrology.py) to a single scalar shared by every member: evaporate `prev_level`'s depth
    (same retention/baseline constants), then grow from this step's own inflow -- the sum of
    `water_deposited` (route_downstream's per-node settled-water output) over every one of this
    lake's own members, since more than one of them can be a true sink once several originally
    separate basins are merged into one lake -- spread as a level rise over the lake's own area
    (member count, the same area proxy used elsewhere in this codebase). Clipped to
    `[floor_elevation, max_depth]` (unbounded above for a closed/endorheic basin, `max_depth is
    None`) -- note `floor_elevation`/`max_depth` are themselves already silt-adjusted by the time
    they reach here (`step_lakes` builds the hierarchy from `elevation + prev_silt_depth`, not
    bare elevation), so nothing here needs to know about silt except to also grow it: the same
    inflow that grows standing water also deposits a small amount of silt on this lake's own bed
    (`SILT_ACCUMULATION_COEFFICIENT`, ~100x smaller -- see that constant's own comment),
    written directly into `out_silt_depth` for every member. A lake with any member currently
    cold enough to accumulate ice (`hydrology.GLACIER_ACCUMULATION_TEMP_C`) freezes solid this
    step -- forced to its own dry floor with no silt growth either (no liquid water to carry
    suspended sediment), matching the old system's same all-or-nothing freeze simplification
    generalized from one node to a lake's whole connected surface."""
    if is_freezing:
        return lake.floor_elevation

    retention = np.exp(-LAKE_EVAPORATION_RATE_PER_MYR * years_myr)
    baseline_loss = LAKE_EVAPORATION_BASELINE_M_PER_MYR * years_myr
    prev_depth = max(0.0, prev_level - lake.floor_elevation)
    carried_depth = max(0.0, prev_depth * retention - baseline_loss)

    inflow = float(water_deposited[lake.members].sum())
    member_count = max(len(lake.members), 1)
    growth_depth = LAKE_FILL_RATE * inflow * years_myr / member_count
    out_silt_depth[lake.members] += SILT_ACCUMULATION_COEFFICIENT * inflow * years_myr / member_count

    cap_depth = (lake.max_depth - lake.floor_elevation) if lake.max_depth is not None else None
    new_depth = carried_depth + growth_depth
    if cap_depth is not None:
        new_depth = min(new_depth, cap_depth)
    return lake.floor_elevation + new_depth


def _prev_level(lake: Lake, elevation: np.ndarray, prev_lake_depth: np.ndarray) -> float:
    """Last step's recorded water elevation for `lake`, considering only members that were
    actually wet (`prev_lake_depth > 0`) -- **not** simply `max(elevation + prev_lake_depth)`
    over every member. A lake's own members routinely include currently-dry higher ground (a
    leaf's catchment includes every node steepest descent drains through, not just its floor;
    a parent's members include the saddle/rim node that defines its own `min_depth`/
    `max_depth`) whose *bare* elevation, with zero recorded depth, can easily sit at or above
    `min_depth` purely because that's what a saddle or rim node's elevation *is* -- taking the
    max over all members regardless of wetness would then read a completely dry lake as
    already full, or a still-separate pair of children as already merged, confirmed directly
    while writing this module's own test fixtures (a leaf's own catchment member can coincide
    with the very saddle node that defines its parent's `min_depth`). Falls back to
    `floor_elevation` (fully dry) when no member was wet at all."""
    members = lake.members
    member_depth = prev_lake_depth[members]
    wet = member_depth > 0.0
    if not wet.any():
        return lake.floor_elevation
    return max(float((elevation[members][wet] + member_depth[wet]).max()), lake.floor_elevation)


def _resolve(
    lake: Lake,
    elevation: np.ndarray,
    prev_lake_depth: np.ndarray,
    water_deposited: np.ndarray,
    years_myr: float,
    is_accumulating: np.ndarray,
    out_lake_depth: np.ndarray,
    out_silt_depth: np.ndarray,
    events: list[str],
) -> float:
    """Resolves one lake (and, for a parent, implicitly its whole subtree) into this step's
    actual water elevation, writing every member's own `out_lake_depth` exactly once, and
    detects a merge/split transition purely from whether last step's recorded depth
    (`prev_lake_depth`) already showed this lake's water reaching its own `min_depth` -- see
    `_prev_level` and this module's own docstring for why that's derived from the persistent
    per-node array rather than a persistent `Lake` registry.

    A leaf, or a parent that **was already** one merged body last step (`prev_level >=
    min_depth`), is resolved as a single connected water body over its own full `members` set.
    If that resolution's *new* level drops back below `min_depth`, it just split this step: both
    children start fresh at exactly `min_depth` (continuity -- the instant they separated, both
    sides were still at that shared height) with no further balance applied to them this same
    step, deferring their independent divergence to next step, the same "changes happen
    gradually across steps" character as everything else persistent in this codebase.

    A parent that **was not** yet merged last step resolves its children independently first;
    if either child's own new level has now reached this parent's `min_depth`, they merge this
    step: continuity again pins the newly-merged body at exactly `min_depth` (any growth beyond
    the saddle in this same step is deferred to next step's now-pooled water balance) rather
    than also folding in that step's remaining growth, which would need choosing an
    otherwise-arbitrary allocation between "how much of this step's inflow arrived before vs.
    after the moment of merging."

    Returns this lake's own resolved level (informational -- only consulted by a parent's own
    merge check above; a lake that just split returns its parent's saddle, not a per-child
    value, since the two children no longer share one level)."""
    is_freezing = bool(is_accumulating[lake.members].any())
    prev_level = _prev_level(lake, elevation, prev_lake_depth)
    already_merged = prev_level >= lake.min_depth  # trivially true for a leaf (min_depth == floor_elevation)

    if lake.children and not already_merged:
        child_levels = [
            _resolve(child, elevation, prev_lake_depth, water_deposited, years_myr, is_accumulating, out_lake_depth, out_silt_depth, events)
            for child in lake.children
        ]
        if max(child_levels) >= lake.min_depth:
            events.append(f"A {len(lake.members)}-node lake formed as two basins merged at elevation {lake.min_depth:.0f}m.")
            out_lake_depth[lake.members] = np.maximum(0.0, lake.min_depth - elevation[lake.members])
            lake.current_water_elevation = lake.min_depth
            lake.is_spilling = lake.max_depth is not None and lake.min_depth >= lake.max_depth
            return lake.min_depth
        lake.current_water_elevation = max(child_levels)
        return lake.current_water_elevation

    new_level = _water_balance(lake, prev_level, water_deposited, years_myr, is_freezing, out_silt_depth)
    if lake.children and new_level < lake.min_depth:
        events.append(f"A {len(lake.members)}-node lake split back into {len(lake.children)} separate basins at elevation {lake.min_depth:.0f}m.")
        for child in lake.children:
            out_lake_depth[child.members] = np.maximum(0.0, lake.min_depth - elevation[child.members])
            child.current_water_elevation = lake.min_depth
            child.is_spilling = False
        lake.current_water_elevation = lake.min_depth
        return lake.min_depth

    out_lake_depth[lake.members] = np.maximum(0.0, new_level - elevation[lake.members])
    lake.current_water_elevation = new_level
    lake.is_spilling = lake.max_depth is not None and new_level >= lake.max_depth
    return new_level


def step_lakes(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    prev_lake_depth: np.ndarray,
    prev_silt_depth: np.ndarray,
    water_deposited: np.ndarray,
    years: float,
    is_accumulating: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[Lake], list[str]]:
    """The full per-step lake update: rebuild this step's depression hierarchy from current
    terrain (`build_lake_hierarchy`) -- using `elevation + prev_silt_depth`, not bare
    `elevation`, so a lake's own floor/rim already reflects everything silted in through the end
    of last step (see `plates.ElevationLine.silt_depth`'s own comment and `_water_balance`'s) --
    then resolve every top-level lake's own water balance (growth/evaporation/merge/split, and
    this step's own silt growth -- see `_resolve`/`_water_balance`). Returns `(lake_depth,
    silt_depth, forest, events)`: `lake_depth` an (N,) array aligned with `elevation` (0
    wherever no lake reaches, exactly the same shape/meaning hydrology.py's old `update_lakes`
    produced, for every other consumer -- render_image.py, coastline.py, erosion.py,
    group_rivers -- to keep reading unchanged); `silt_depth` this step's updated silt array
    (monotonically >= `prev_silt_depth` -- silt never erodes back away); `forest` this step's
    fully-resolved `Lake` tree (each lake's own `current_water_elevation`/`is_spilling` set to
    its outcome this step -- informational only, e.g. for a future per-lake stats view; nothing
    in this module reads it back next step, see the module docstring for why no persistent
    registry is needed); and `events` human-readable merge/split transition strings for the
    caller to log (mirroring merge_split.apply_topology_changes's own pattern -- this module
    never calls `world.log_event` directly, keeping it testable without a `World`)."""
    n = len(elevation)
    lake_depth = np.zeros(n)
    if n == 0:
        return lake_depth, np.zeros(n), [], []

    effective_elevation = elevation + prev_silt_depth
    forest = build_lake_hierarchy(effective_elevation, is_ocean, neighbor_idx)
    years_myr = years / 1_000_000.0
    silt_depth = np.array(prev_silt_depth, dtype=float, copy=True)

    events: list[str] = []
    for root in forest:
        _resolve(root, effective_elevation, prev_lake_depth, water_deposited, years_myr, is_accumulating, lake_depth, silt_depth, events)
    return lake_depth, silt_depth, forest, events

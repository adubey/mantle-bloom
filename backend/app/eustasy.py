"""Eustatic (whole-ocean) sea level: `World.sea_level_m` derived from a conserved ocean
water *volume* against the world's current hypsometry, rather than pinned at a fixed value.

Every `is_ocean` / coastline / bathymetry check in the codebase reads `world.sea_level_m`
directly (see `World.sea_level_m`), so making sea level respond to tectonics is just
recomputing that one number each step: tectonics deepens ocean basins (sea-floor spreading,
subduction) and thins/drowns continental crust, erosion planes highlands into the sea -- all
of which change how much basin volume the same water fills, hence where the shoreline sits.
Without this, a fixed sea level means every bit of continental subsidence or new deep ocean
floor is a permanent, uncompensated loss of dry land (docs/TODO.md "Land fraction slowly
declines"): on the real Earth, opening an ocean basin drops sea level and hands that land
back as continental freeboard.

Model. Every lattice node covers the same area by construction (`lithosphere.node_area_m2`
is a function of spacing only), so the ocean's water volume is proportional to the summed
water column `W_ocean = sum_i max(0, sea_level - z_i)` over every node *actually part of the
connected ocean* -- see `_ocean_connected_mask` below for why "actually part of," not just
"below `sea_level`," matters.

`World.ocean_water_column_m` holds the world's *total* surface-water budget `W_total`,
snapshot once at generation (from the flat starting sea level, with no ice or lakes yet) and
conserved forever after. Each step `update_sea_level` splits that budget into the part locked
up on land -- `W_trapped`, the water frozen into ice caps / glaciers / mountain-top ice plus
what's standing in lakes and seas (`trapped_water_column_m`) -- and the rest, which is the
ocean, then solves the monotonic 1-D equation `ocean_water_column(h) == W_total - W_trapped`
for the new `h`. Adding deep ocean floor raises `ocean_water_column(h)` at every `h`, so the
solved `h` drops -- the eustatic fall a spreading basin produces; growing an ice age's ice
sheets, or a landlocked sea filling up, raises `W_trapped`, which drops `h` the same way --
glacio-eustasy, the ~120 m Pleistocene sea-level swing, generalized to any water taken out of
circulation. The user's sea-level slider (`POST /world/controls`) sets `W_total` to whatever
value floats the *current* hypsometry (and its *current* trapped water) at the requested
level, which then persists and is itself conserved going forward.

**Why the ocean's own volume must be connectivity-restricted, not just "below `h`."** A
landlocked lake or sea's floor routinely sits below the world's current `sea_level_m` -- a
rift basin that dropped out from under an old shoreline, or simply an ocean-floor-deep trench
a rising barrier has since sealed off from the world ocean (see hydrology.connected_ocean_mask
/lakes.py's own SEA_MAX_DEPTH_M comment for a real example: a 4.84M km^2, 4,626 m deep closed
basin on a 286 My save). That basin's own water is already counted, for real, via
`trapped_water_column_m`'s measured `lake_depth` -- but naively summing `max(0, h - z_i)` over
*every* node below `h`, connected or not, double-books it: the solve then reads that basin's
entire sub-`h` volume as "free" ocean depth requiring no water budget at all, so *less* of the
conserved budget is needed to reach the same `h` -- silently lowering every world's solved sea
level below its true value, worse the deeper and more numerous its closed basins get over a
long run (confirmed directly on a 161 My save, seed 349936951: restricting the sum to the
connected ocean raised the solved level by over 100 m). `_ocean_connected_mask` is what
excludes those nodes; `ocean_water_column_m`/`_solve_sea_level_connected` are `total_water_
column_m`/`solve_sea_level` narrowed to it, used by every `World`-level entry point below.
`total_water_column_m`/`solve_sea_level` themselves stay the bare, connectivity-oblivious
primitives (exercised directly by this module's own unit tests) -- restricting *which* nodes
they sum over is entirely the caller's job.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from . import hydrology

if TYPE_CHECKING:
    from .world import World

# Bisection tolerance for the sea-level solve, in metres of the summed water column (`W`).
# `W` is ~1e5-1e7 m·nodes for a real world, so this is a very tight relative tolerance and
# the solve still converges in ~40 iterations (each a single vectorized sum).
_SOLVE_TOLERANCE_M = 1.0
_SOLVE_MAX_ITERS = 80

# Ice is ~91.7% the density of liquid water, so a metre of `glacier_depth` (stored as a
# metre-of-ice column, see hydrology.py) is this much liquid-water column removed from the
# ocean. Lake water is already liquid, so it counts 1:1. Applied in `trapped_water_column_m`.
_ICE_WATER_EQUIVALENT = 0.917

# Lake water only counts toward the trapped budget where it's a real, visible standing body,
# not the sub-metre coastal-pond dither a low-relief shelf produces (the same threshold
# hydrology.LAKE_MIN_VISIBLE_DEPTH_M draws that line at) -- "large lakes", loosely. The
# per-node depths that clear it are summed directly; no connected-component size test, since
# a genuinely large lake is exactly a wide connected patch of cleared nodes anyway.
_LAKE_MIN_TRAPPED_DEPTH_M = 1.0


def all_elevations(world: "World") -> np.ndarray:
    """Every node's current live elevation across every plate, concatenated -- the hypsometry
    the water volume is filled against."""
    from .plates import collect_all_elevation

    return collect_all_elevation(world.plates)


def total_water_column_m(elevations: np.ndarray, sea_level_m: float) -> float:
    """`sum_i max(0, sea_level - z_i)` -- the summed depth of water standing over the whole
    lattice at `sea_level_m`. Proportional to ocean volume (every node has equal area).
    Strictly increasing in `sea_level_m`."""
    if len(elevations) == 0:
        return 0.0
    return float(np.sum(np.clip(sea_level_m - elevations, 0.0, None)))


def _seeded_connected_mask(
    elevations: np.ndarray, sea_level_m: float, neighbor_idx: np.ndarray, seed: np.ndarray
) -> np.ndarray:
    """Below-`sea_level_m` connectivity mask anchored to whichever below-`sea_level_m`
    component(s) currently contain a `seed` node -- last step's own real `is_ocean` -- rather
    than `hydrology.connected_ocean_mask`'s own "largest component wins" rule.

    Why not just call `connected_ocean_mask` directly at every candidate `h`: a bisection
    candidate can land well below the *real* ocean's typical floor while still being at or
    above some small closed pit's own (deeper) floor -- at that candidate, the pit is the
    *only* below-`h` component, or briefly the largest one, and "largest wins" would label the
    pit itself "ocean" for that one candidate. Since raising `h` further eventually submerges
    the real ocean and flips the largest-component label back, the naive largest-wins volume
    is not monotonic in `h` across the *whole* bisection bracket -- breaking the bisection this
    is used from (confirmed directly: a synthetic fixture with a pit deeper than the "ocean"
    converges to a wrong root under the naive rule). Anchoring to a fixed seed set instead is
    genuinely monotonic: raising `h` can only ever merge more nodes into the seed's own
    component (never remove one, see `connected_ocean_mask`'s own docstring for the same
    argument), so "how much water sits in the seed's component" only ever grows with `h`. This
    also transitively keeps a real second ocean recognized (`hydro.is_ocean` already includes
    one, from `connected_ocean_mask`'s own `SECOND_OCEAN_SIZE_FRACTION` handling, so the seed
    set already covers it) and correctly grows to include a landlocked lake/sea the instant
    `h` actually rises enough to reconnect it -- both sides of that same boundary edge then
    share one label, same as any other merge.

    Returns all-`False` when no seed node is below `h` at all (the seed's own component hasn't
    reached down this far at this candidate) -- correctly zero ocean volume for a candidate
    this deep, not a mislabeled non-ocean pit standing in for it."""
    below = elevations <= sea_level_m
    n = len(elevations)
    if n == 0 or not below.any():
        return below
    seed_below = seed & below
    if not seed_below.any():
        return np.zeros(n, dtype=bool)
    k = neighbor_idx.shape[1]
    rows = np.repeat(np.arange(n), k)
    cols = neighbor_idx.ravel()
    keep = below[rows] & below[cols]
    graph = coo_matrix((np.ones(int(np.count_nonzero(keep))), (rows[keep], cols[keep])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    ocean_labels = np.unique(labels[seed_below])
    return below & np.isin(labels, ocean_labels)


def _ocean_connected_mask(world: "World", elevations: np.ndarray, sea_level_m: float) -> np.ndarray:
    """Per-node bool: is this node actually part of the connected world ocean at
    `sea_level_m` -- not merely sitting below it (see this module's own docstring for why
    that distinction is what keeps a landlocked lake/sea's volume from being double-counted).

    Anchored (`_seeded_connected_mask`) to `world.hydrology_cache.is_ocean` -- last step's own
    real ocean classification -- rather than recomputing `hydrology.connected_ocean_mask`'s
    "largest wins" rule fresh at every candidate `h` (see `_seeded_connected_mask`'s own
    docstring for why that matters for the bisection this feeds). Reuses whatever k-NN graph
    the cache already built this step (`hydro.neighbor_idx`, guaranteed the same node order as
    `elevations` -- both ultimately walk `world.plates` in the same plate-major order, see
    plates.gather_node_positions) rather than paying to rebuild one from scratch.

    Falls back to `hydrology.connected_ocean_mask` itself (no seed available yet) when the
    cache has no usable `is_ocean` of its own -- e.g. a world with no ocean at all -- and to
    the bare below-`sea_level_m` mask, the old connectivity-oblivious behaviour, whenever
    there's no usable cache at all: world generation (no step has run), or a save from before
    `hydrology_cache` existed."""
    hydro = getattr(world, "hydrology_cache", None)
    if hydro is None or len(hydro.points) != len(elevations):
        return elevations <= sea_level_m
    seed = getattr(hydro, "is_ocean", None)
    if seed is None or len(seed) != len(elevations) or not seed.any():
        return hydrology.connected_ocean_mask(hydro.points, elevations, sea_level_m, hydro.neighbor_idx)
    return _seeded_connected_mask(elevations, sea_level_m, hydro.neighbor_idx, seed)


def ocean_water_column_m(world: "World", elevations: np.ndarray, sea_level_m: float) -> float:
    """`total_water_column_m`, restricted to nodes that are actually the connected ocean at
    `sea_level_m` (see `_ocean_connected_mask`) -- the quantity the eustasy solve actually
    means by "the ocean's own volume." A closed basin's own below-`sea_level_m` footprint
    contributes nothing here even when its floor plunges far deeper than `sea_level_m` itself:
    it either holds no water at all (a dry catchment, however deep), or holds real, measured
    water already accounted for by `trapped_water_column_m`."""
    mask = _ocean_connected_mask(world, elevations, sea_level_m)
    return total_water_column_m(elevations[mask], sea_level_m)


def _solve_sea_level_connected(world: "World", elevations: np.ndarray, water_column_m: float) -> float:
    """`solve_sea_level`, but bisecting against `ocean_water_column_m` (connectivity-restricted)
    rather than the bare `total_water_column_m` -- see this module's own docstring for why a
    candidate level can't be allowed to "flood" a disconnected closed basin for free. Recomputes
    the connected-ocean mask at every candidate `h` (connectivity is itself a function of `h`:
    raising it can only ever merge components, never split one -- see `connected_ocean_mask`'s
    own docstring -- so `ocean_water_column_m(h)` stays strictly increasing and bisection stays
    valid), which is why this isn't simply `total_water_column_m` over one fixed mask."""
    if len(elevations) == 0:
        return 0.0
    lo = float(np.min(elevations))
    if water_column_m <= 0.0:
        return lo
    hi = float(np.max(elevations)) + water_column_m / len(elevations) + 1.0
    for _ in range(_SOLVE_MAX_ITERS):
        mid = 0.5 * (lo + hi)
        w = ocean_water_column_m(world, elevations, mid)
        if abs(w - water_column_m) < _SOLVE_TOLERANCE_M:
            return mid
        if w < water_column_m:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def water_column_for_sea_level(world: "World", sea_level_m: float) -> float:
    return ocean_water_column_m(world, all_elevations(world), sea_level_m)


def trapped_water_column_m(world: "World") -> float:
    """Surface water currently locked up on land rather than sitting in the ocean, in the same
    "summed column over equal-area nodes" units as `total_water_column_m` -- so it can be
    subtracted straight from the conserved total budget. Two contributions, both read off the
    persisted per-node fields hydrology.py maintains: every node's `glacier_depth` (ice caps,
    valley glaciers, mountain-top ice, and any ice-age sea-ice cap -- converted from a column
    of ice to its liquid-water equivalent, `_ICE_WATER_EQUIVALENT`), and the `lake_depth` of
    every node holding a real visible lake (`_LAKE_MIN_TRAPPED_DEPTH_M`). Both are water that
    fell as precipitation and would otherwise have run back to the sea; debiting them is what
    makes sea level fall as an ice age's ice sheets grow (glacio-eustasy). 0.0 for a world
    with neither."""
    from .plates import collect_all_glacier_depth, collect_all_lake_depth

    if not world.plates:
        return 0.0
    glacier_depth = collect_all_glacier_depth(world.plates)
    lake_depth = collect_all_lake_depth(world.plates)
    ice = _ICE_WATER_EQUIVALENT * float(np.sum(np.clip(glacier_depth, 0.0, None)))
    lakes = float(np.sum(np.where(lake_depth >= _LAKE_MIN_TRAPPED_DEPTH_M, lake_depth, 0.0)))
    return ice + lakes


def solve_sea_level(elevations: np.ndarray, water_column_m: float) -> float:
    """The `h` with `total_water_column_m(elevations, h) == water_column_m`. Bisection on
    `[min z, max z + headroom]` -- `total_water_column_m` is 0 at `min z` and unbounded
    above, so a bracket always exists for any non-negative target."""
    if len(elevations) == 0:
        return 0.0
    lo = float(np.min(elevations))
    if water_column_m <= 0.0:
        return lo
    # Upper bracket: enough headroom that every node is submerged and then some.
    hi = float(np.max(elevations)) + water_column_m / len(elevations) + 1.0
    for _ in range(_SOLVE_MAX_ITERS):
        mid = 0.5 * (lo + hi)
        w = total_water_column_m(elevations, mid)
        if abs(w - water_column_m) < _SOLVE_TOLERANCE_M:
            return mid
        if w < water_column_m:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def initialize_water_budget(world: "World") -> None:
    """Snapshot the ocean water volume from the current hypsometry + current `sea_level_m` --
    called once at generation, and as the backfill for a save written before this field
    existed. Idempotent-safe: recomputing from an already-eustatic world reproduces the same
    budget (the solve is the exact inverse). Includes any water already trapped in ice/lakes
    (0 at generation, non-zero when this is the backfill for a mid-run save), so the stored
    number is always the world's *total* surface water."""
    world.ocean_water_column_m = water_column_for_sea_level(world, world.sea_level_m) + trapped_water_column_m(world)


def update_sea_level(world: "World") -> None:
    """Re-solve `world.sea_level_m` against this step's hypsometry, holding the *total*
    surface-water budget (`world.ocean_water_column_m`) fixed and subtracting the part
    currently locked up in ice caps / glaciers / lakes and seas (`trapped_water_column_m`) so
    only the ocean's own share floats the shoreline -- solved against `ocean_water_column_m`
    (connectivity-restricted, see this module's own docstring), not the bare `total_water_
    column_m`, so a landlocked basin's own below-sea-level footprint isn't also counted as
    free ocean depth on top of being debited as trapped. Cheap enough to call unconditionally
    every step. Initializes the budget on first use (older saves / a freshly constructed
    World)."""
    if getattr(world, "ocean_water_column_m", None) is None:
        initialize_water_budget(world)
        return
    elevations = all_elevations(world)
    ocean_budget = max(0.0, world.ocean_water_column_m - trapped_water_column_m(world))
    world.sea_level_m = _solve_sea_level_connected(world, elevations, ocean_budget)


def set_sea_level_via_water_budget(world: "World", sea_level_m: float) -> None:
    """The `POST /world/controls` sea-level slider: interpret "put sea level at X" as "add or
    remove ocean water until the *current* hypsometry floats at X", then let that new *total*
    water volume (ocean at X, plus whatever's currently trapped in ice/lakes) be conserved
    going forward."""
    world.sea_level_m = float(sea_level_m)
    world.ocean_water_column_m = water_column_for_sea_level(world, world.sea_level_m) + trapped_water_column_m(world)

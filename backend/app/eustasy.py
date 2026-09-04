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
water column `W_ocean = sum_i max(0, sea_level - z_i)` over every node in the world.

`World.ocean_water_column_m` holds the world's *total* surface-water budget `W_total`,
snapshot once at generation (from the flat starting sea level, with no ice or lakes yet) and
conserved forever after. Each step `update_sea_level` splits that budget into the part locked
up on land -- `W_trapped`, the water frozen into ice caps / glaciers / mountain-top ice plus
what's standing in lakes (`trapped_water_column_m`) -- and the rest, which is the ocean, then
solves the monotonic 1-D equation `total_water_column(h) == W_total - W_trapped` for the new
`h`. Adding deep ocean floor raises `total_water_column(h)` at every `h`, so the solved `h`
drops -- the eustatic fall a spreading basin produces; growing an ice age's ice sheets raises
`W_trapped`, which drops `h` the same way -- glacio-eustasy, the ~120 m Pleistocene sea-level
swing. The user's sea-level slider (`POST /world/controls`) sets `W_total` to whatever value
floats the *current* hypsometry (and its *current* trapped water) at the requested level,
which then persists and is itself conserved going forward.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

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


def water_column_for_sea_level(world: "World", sea_level_m: float) -> float:
    return total_water_column_m(all_elevations(world), sea_level_m)


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
    currently locked up in ice caps / glaciers / large lakes (`trapped_water_column_m`) so
    only the ocean's own share floats the shoreline. Cheap enough to call unconditionally
    every step. Initializes the budget on first use (older saves / a freshly constructed
    World)."""
    if getattr(world, "ocean_water_column_m", None) is None:
        initialize_water_budget(world)
        return
    elevations = all_elevations(world)
    ocean_budget = max(0.0, world.ocean_water_column_m - trapped_water_column_m(world))
    world.sea_level_m = solve_sea_level(elevations, ocean_budget)


def set_sea_level_via_water_budget(world: "World", sea_level_m: float) -> None:
    """The `POST /world/controls` sea-level slider: interpret "put sea level at X" as "add or
    remove ocean water until the *current* hypsometry floats at X", then let that new *total*
    water volume (ocean at X, plus whatever's currently trapped in ice/lakes) be conserved
    going forward."""
    world.sea_level_m = float(sea_level_m)
    world.ocean_water_column_m = water_column_for_sea_level(world, world.sea_level_m) + trapped_water_column_m(world)

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
water column `W = sum_i max(0, sea_level - z_i)` over every node in the world. `W` is fixed
at generation (from the flat starting sea level) and stored on `World.ocean_water_column_m`;
each step `update_sea_level` solves the monotonic 1-D equation `total_water_column(h) == W`
for the new `h`. Adding deep ocean floor raises `total_water_column(h)` at every `h`, so the
solved `h` drops -- the eustatic fall a spreading basin produces. The user's sea-level slider
(`POST /world/controls`) now sets `W` to whatever value floats the *current* hypsometry at
the requested level, i.e. it adds or removes ocean water (glacio-eustasy / a bigger ocean),
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
    budget (the solve is the exact inverse)."""
    world.ocean_water_column_m = water_column_for_sea_level(world, world.sea_level_m)


def update_sea_level(world: "World") -> None:
    """Re-solve `world.sea_level_m` against this step's hypsometry, holding the ocean water
    volume (`world.ocean_water_column_m`) fixed. Cheap enough to call unconditionally every
    step. Initializes the budget on first use (older saves / a freshly constructed World)."""
    if getattr(world, "ocean_water_column_m", None) is None:
        initialize_water_budget(world)
        return
    elevations = all_elevations(world)
    world.sea_level_m = solve_sea_level(elevations, world.ocean_water_column_m)


def set_sea_level_via_water_budget(world: "World", sea_level_m: float) -> None:
    """The `POST /world/controls` sea-level slider: interpret "put sea level at X" as "add or
    remove ocean water until the *current* hypsometry floats at X", then let that new water
    volume be conserved going forward."""
    world.sea_level_m = float(sea_level_m)
    world.ocean_water_column_m = water_column_for_sea_level(world, world.sea_level_m)

"""Volcanic eruption lifecycle for existing volcano nodes.

Volcanic fields themselves are now created directly by `PlateWithLines.deform` (see
plates.py) when a rift boundary has stretched too thin to keep filling with plain
ridge/rift crust -- detection/spawning/merging/isolated-growth of whole volcanic-field
*plates* used to live here as a periodic clean-up pass, but that's subsumed by deform()'s
own per-turn rift handling now (see docs/simulation-model.md and the plan this replaced).

What's left here is the *per-node* eruption lifecycle, run every step regardless of how a
volcano node came to exist: each individual volcano point has its own
`volcano_active_years_remaining` (`elevation_lines.VOLCANO_ACTIVE_MIN/MAX_YEARS`, drawn once
at creation), decremented every step. While active, it rolls a per-step eruption chance
(`1 - exp(-ERUPTION_RATE_PER_MYR * active_years_this_step / 1e6)`) and, if it erupts, adds
`elevation_lines.ERUPTION_ELEVATION_M` of new land and grows `mineral_deposit_m`.
Deterministic per `(seed, elapsed_years, plate_id, line_index)`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .elevation_lines import ELEV_CHANGE_VOLCANO, ERUPTION_ELEVATION_M, MAX_ELEVATION_M, MIN_ELEVATION_M
from .plates import PlateWithLines

if TYPE_CHECKING:
    from .world import World

# Expected number of eruption events over a volcano's full active life is
# ERUPTION_RATE_PER_MYR * (active life in Myr) -- e.g. at the low end of VOLCANO_ACTIVE
# (0.1 Myr), that's 0.3 expected events (p(>=1) ~= 26%, so most short-lived volcanoes erupt
# zero or one time); at the high end (1 Myr), 3 expected events. "Occasionally," not "every
# step" or "constantly."
ERUPTION_RATE_PER_MYR = 3.0

# Mineral deposits: real hydrothermal circulation around an active volcanic vent precipitates
# metal-rich ore (porphyry-copper/VMS-style deposits), so mineral_deposit_m is grown right
# here, at the same eruption roll that already adds ERUPTION_ELEVATION_M -- "an eruption
# deposits mineral-rich material" is exactly what that mask already means, no separate
# detection pass needed. Monotonically non-decreasing (see plates.ElevationLine), same
# self-reinforcing convention silt_depth/coal_deposit_m/oil_gas_deposit_m already use.
MINERAL_DEPOSIT_PER_ERUPTION_M = 0.5
MAX_MINERAL_DEPOSIT_M = 20.0


def apply_volcanic_activity(world: "World", years: float) -> None:
    """Every step: rolls each individual active volcano's own eruption chance, adding
    ERUPTION_ELEVATION_M wherever it erupts. Mutates world.plates in place."""
    for plate in world.plates:
        _apply_volcanic_activity_to_lines(plate, world, years)


def _apply_volcanic_activity_to_lines(plate: PlateWithLines, world: "World", years: float) -> None:
    """`PlateWithLines`' own per-line eruption roll."""
    for line_index, line in enumerate(plate.lines):
        if len(line) == 0 or not np.any(line.is_volcano):
            continue
        active_mask = line.is_volcano & (line.volcano_active_years_remaining > 0)
        if not np.any(active_mask):
            continue

        active_years_this_step = np.minimum(years, line.volcano_active_years_remaining)
        p_erupt = 1.0 - np.exp(-ERUPTION_RATE_PER_MYR * active_years_this_step / 1_000_000.0)
        rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate.plate_id, line_index))
        erupts = active_mask & (rng.random(len(line)) < p_erupt)

        new_elevation = line.elevation.copy()
        new_elevation[erupts] += ERUPTION_ELEVATION_M
        new_elevation = np.clip(new_elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)
        new_remaining = np.clip(line.volcano_active_years_remaining - years, 0.0, None)
        new_mineral_deposit = np.clip(
            line.mineral_deposit_m + np.where(erupts, MINERAL_DEPOSIT_PER_ERUPTION_M, 0.0), 0.0, MAX_MINERAL_DEPOSIT_M
        )
        # Elevation-change provenance (diagnostic only -- see elevation_lines.ELEV_CHANGE_*):
        # an eruption always adds ERUPTION_ELEVATION_M, well past the min-delta threshold.
        new_reason = np.where(erupts, ELEV_CHANGE_VOLCANO, line.elev_change_reason)

        # theta unchanged -- line.replace copies every other field (including
        # channel_width) from the existing line automatically. See plates.ElevationLine's
        # own docstring for why this pattern replaced explicit field-by-field
        # reconstruction here.
        plate.replace_line(
            line_index,
            line.replace(
                elevation=new_elevation,
                volcano_active_years_remaining=new_remaining,
                mineral_deposit_m=new_mineral_deposit,
                elev_change_reason=new_reason,
            ),
        )


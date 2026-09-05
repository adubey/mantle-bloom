"""Volcanic eruption lifecycle for existing volcano nodes.

Volcanic fields themselves are now created directly by `LithospherePlate.deform` (see
lithosphere_plate.py) when a rift boundary has stretched too thin to keep filling with plain
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
from scipy.spatial import cKDTree

from . import geometry
from .elevation_lines import (
    ELEV_CHANGE_MIN_DELTA_M,
    ELEV_CHANGE_VOLCANIC_PLAIN,
    ELEV_CHANGE_VOLCANO,
    ERUPTION_ELEVATION_M,
    MAX_ELEVATION_M,
    MIN_ELEVATION_M,
    PLANET_RADIUS_KM,
    VOLCANIC_PLAIN_ELEVATION_M,
    VOLCANIC_PLAIN_REACH_KM,
)
from .plates import PlateWithLines

if TYPE_CHECKING:
    from .world import World

# Expected number of eruption events over a volcano's full active life is
# ERUPTION_RATE_PER_MYR * (active life in Myr) -- e.g. at the low end of VOLCANO_ACTIVE
# (0.1 Myr), that's 0.5 expected events (p(>=1) ~= 39%, so most short-lived volcanoes still
# erupt zero or one time); at the high end (1 Myr), 5 expected events -> up to
# ERUPTION_ELEVATION_M * 5 = 1500 m of gross relief before dormancy. "Occasionally," not
# "every step" or "constantly." Raised alongside ERUPTION_ELEVATION_M 2026-09-04 -- see that
# constant's own comment (docs/TODO.md "Land fraction slowly declines") -- volcanism was
# contributing next to nothing to land at the old rate.
ERUPTION_RATE_PER_MYR = 5.0

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
    ERUPTION_ELEVATION_M wherever it erupts, then spreads a broader, weaker volcanic-plain
    apron around each vent that erupted this step. Mutates world.plates in place."""
    for plate in world.plates:
        erupted_points = _apply_volcanic_activity_to_lines(plate, world, years)
        if erupted_points:
            _spread_volcanic_plains(plate, world, years, erupted_points)


def _apply_volcanic_activity_to_lines(plate: PlateWithLines, world: "World", years: float) -> list[np.ndarray]:
    """`PlateWithLines`' own per-line eruption roll. Returns the world-space positions of
    every node that erupted this step (across all of the plate's lines), for
    `_spread_volcanic_plains` to spread an apron around."""
    erupted_points: list[np.ndarray] = []
    for line_index, line in enumerate(plate.lines):
        if len(line) == 0 or not np.any(line.is_volcano):
            continue
        active_mask = line.is_volcano & (line.volcano_active_years_remaining > 0)
        if not np.any(active_mask):
            continue

        active_years_this_step = np.minimum(years, line.volcano_active_years_remaining)
        # world.volcanism_multiplier (the "Controls" tuning knob, 1.0 == untuned) scales both
        # the per-step eruption probability *and* the elevation each eruption adds below, so a
        # single knob controls total volcanic land-building. 0.0 -> p_erupt == 0 everywhere.
        p_erupt = 1.0 - np.exp(
            -ERUPTION_RATE_PER_MYR * world.volcanism_multiplier * active_years_this_step / 1_000_000.0
        )
        rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate.plate_id, line_index))
        erupts = active_mask & (rng.random(len(line)) < p_erupt)

        new_elevation = line.elevation.copy()
        new_elevation[erupts] += ERUPTION_ELEVATION_M * world.volcanism_multiplier
        new_elevation = np.clip(new_elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)
        new_remaining = np.clip(line.volcano_active_years_remaining - years, 0.0, None)
        new_mineral_deposit = np.clip(
            line.mineral_deposit_m + np.where(erupts, MINERAL_DEPOSIT_PER_ERUPTION_M, 0.0), 0.0, MAX_MINERAL_DEPOSIT_M
        )
        # Elevation-change provenance (diagnostic only -- see elevation_lines.ELEV_CHANGE_*):
        # an eruption always adds ERUPTION_ELEVATION_M, well past the min-delta threshold.
        new_reason = np.where(erupts, ELEV_CHANGE_VOLCANO, line.elev_change_reason)

        if np.any(erupts):
            erupted_points.append(geometry.to_world(plate.frame, geometry.local_xyz(np.full(len(line), line.phi), line.theta))[erupts])

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
    return erupted_points


def _spread_volcanic_plains(plate: PlateWithLines, world: "World", years: float, erupted_points: list[np.ndarray]) -> None:
    """Spread a broad, low-relief apron around every vent that erupted this step -- a
    flood-basalt/shield-flank plain, distinct from the sharp point bump `_apply_volcanic_
    activity_to_lines` already applied there. Tapers linearly from
    VOLCANIC_PLAIN_ELEVATION_M at the vent to 0 at VOLCANIC_PLAIN_REACH_KM, same taper shape
    `faults._apply_plate_fault_relief` uses for a fault's own relief. Where two aprons
    overlap this step, the *larger* contribution wins (not the sum) -- a cluster of vents
    erupting the same step should read as one coalesced apron, not a runaway stack."""
    own_points, _ = plate.all_points_and_elevation()
    if len(own_points) == 0:
        return
    reach_rad = VOLCANIC_PLAIN_REACH_KM / PLANET_RADIUS_KM
    tree = cKDTree(own_points, balanced_tree=False, compact_nodes=False)
    vent_points = np.concatenate(erupted_points, axis=0)

    delta = np.zeros(len(own_points))
    for vent in vent_points:
        affected = tree.query_ball_point(vent, reach_rad)
        if not affected:
            continue
        affected = np.asarray(affected)
        d = np.linalg.norm(own_points[affected] - vent, axis=-1)
        taper = np.clip(1.0 - d / reach_rad, 0.0, 1.0)
        contrib = VOLCANIC_PLAIN_ELEVATION_M * world.volcanism_multiplier * taper
        np.maximum.at(delta, affected, contrib)

    if not np.any(delta):
        return

    new_lines = []
    offset = 0
    for line in plate.lines:
        n = len(line)
        seg_delta = delta[offset : offset + n]
        offset += n
        if not np.any(seg_delta):
            new_lines.append(line)
            continue
        new_elev = np.clip(line.elevation + seg_delta, MIN_ELEVATION_M, MAX_ELEVATION_M)
        moved = np.abs(new_elev - line.elevation) >= ELEV_CHANGE_MIN_DELTA_M
        # Don't downgrade the vent's own sharper VOLCANO stamp to the plain's -- only claim
        # nodes the point bump didn't already touch this step.
        new_reason = np.where(moved & (line.elev_change_reason != ELEV_CHANGE_VOLCANO), ELEV_CHANGE_VOLCANIC_PLAIN, line.elev_change_reason)
        new_lines.append(line.replace(elevation=new_elev, elev_change_reason=new_reason))
    plate.set_lines(new_lines)


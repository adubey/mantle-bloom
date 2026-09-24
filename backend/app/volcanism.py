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

from . import geometry, lithosphere
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
from .plates import Plate, PlateWithLines

if TYPE_CHECKING:
    from .world import World

# Expected number of eruption events over a volcano's full active life is
# ERUPTION_RATE_PER_MYR * (active life in Myr) -- e.g. at the low end of VOLCANO_ACTIVE
# (0.1 Myr), that's 0.5 expected events (p(>=1) ~= 39%, so most short-lived volcanoes still
# erupt zero or one time); at the high end (1 Myr), 5 expected events -> up to
# ERUPTION_ELEVATION_M * 5 = 1500 m of gross relief before dormancy. "Occasionally," not
# "every step" or "constantly." Raised alongside ERUPTION_ELEVATION_M 2026-09-04 -- see that
# constant's own comment (GitHub issue #120, "Land fraction slowly declines") -- volcanism was
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

        # Issue #173: an eruption's elevation gain has to come with a matching
        # crustal_thickness_m (Hc) addition, or it's "phantom" relief no crustal mass backs --
        # exactly the gap that let erosion's slope-driven terms tear down more real crust than
        # volcanism ever added. See lithosphere.back_elevation_gain's own docstring for why
        # this is anchored to the column's isostatic equilibrium rather than raw elevation.
        new_crustal_thickness, new_elevation = lithosphere.back_elevation_gain(
            line, plate, ERUPTION_ELEVATION_M * world.volcanism_multiplier, erupts
        )
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
                crustal_thickness_m=new_crustal_thickness,
                volcano_active_years_remaining=new_remaining,
                mineral_deposit_m=new_mineral_deposit,
                elev_change_reason=new_reason,
            ),
        )
    return erupted_points


def _spread_volcanic_plains(plate: Plate, world: "World", years: float, erupted_points: list[np.ndarray]) -> None:
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

    elevation = plate.collect("elevation")
    new_hc, new_elevation = lithosphere.back_elevation_gain_fields(
        elevation,
        plate.collect("crustal_thickness_m"),
        plate.collect("mantle_lithosphere_thickness_m"),
        plate.collect("crust_type_code"),
        plate.crust_type,
        delta,
        delta > 0.0,
    )
    moved = np.abs(new_elevation - elevation) >= ELEV_CHANGE_MIN_DELTA_M
    old_reason = plate.collect("elev_change_reason")
    new_reason = np.where(
        moved & (old_reason != ELEV_CHANGE_VOLCANO),
        ELEV_CHANGE_VOLCANIC_PLAIN,
        old_reason,
    )
    plate.set_fields_on_plate(
        elevation=new_elevation,
        crustal_thickness_m=new_hc,
        elev_change_reason=new_reason,
    )

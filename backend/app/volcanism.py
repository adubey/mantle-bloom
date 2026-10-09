"""Volcanic eruption lifecycle for existing volcano nodes.

Volcanic fields themselves are created directly by the tectonic engine (see
quad_tectonics.py and lithosphere_plate.py) when a rift boundary has stretched too thin to keep filling with plain
ridge/rift crust -- detection/spawning/merging/isolated-growth of whole volcanic-field
*plates* used to live here as a periodic clean-up pass, but that's subsumed by deform()'s
own per-turn rift handling now (see docs/simulation-model.md and the plan this replaced).

What's left here is the *per-node* eruption lifecycle, run every step regardless of how a
volcano node came to exist: each individual volcano point has its own
`volcano_active_years_remaining` (`elevation_lines.VOLCANO_ACTIVE_MIN/MAX_YEARS`, drawn once
at creation), decremented every step. While active, it rolls a per-step eruption chance
(`1 - exp(-ERUPTION_RATE_PER_MYR * active_years_this_step / 1e6)`) and, if it erupts, adds
`elevation_lines.ERUPTION_ELEVATION_M` of new land and grows `mineral_deposit_m`.
Deterministic per `(seed, elapsed_years, plate_id, node ID)` (`_node_uniforms`).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import continental_ledger, lithosphere, mobile_cover
from .elevation_lines import (
    ELEV_CHANGE_MIN_DELTA_M,
    ELEV_CHANGE_VOLCANIC_PLAIN,
    ELEV_CHANGE_VOLCANO,
    ERUPTION_ELEVATION_M,
    effective_is_continental_from_codes,
    PLANET_RADIUS_KM,
    VOLCANIC_PLAIN_ELEVATION_M,
    VOLCANIC_PLAIN_REACH_KM,
)
from .plates import Plate

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
# detection pass needed. Monotonically non-decreasing, same
# self-reinforcing convention coal_deposit_m/oil_gas_deposit_m already use.
MINERAL_DEPOSIT_PER_ERUPTION_M = 0.5
MAX_MINERAL_DEPOSIT_M = 20.0


def apply_volcanic_activity(world: "World", years: float) -> None:
    """Every step: rolls each individual active volcano's own eruption chance, adding
    ERUPTION_ELEVATION_M wherever it erupts, then spreads a broader, weaker volcanic-plain
    apron around each vent that erupted this step. Erupted crust buries the mobile cover it
    lands on, consolidating it in proportion to the surface lava's thickness, the elevation it
    adds (mobile_cover.VOLCANIC_SEAL_THICKNESS_M). Mutates world.plates in place."""
    continental_ledger.ensure_initialized(world)
    for plate in world.plates:
        hc_before = plate.collect("crustal_thickness_m")
        elevation_before = plate.collect("elevation")
        erupted_points = _apply_volcanic_activity_to_surface(plate, world, years)
        if erupted_points:
            _spread_volcanic_plains(plate, world, years, erupted_points)
        hc_after = plate.collect("crustal_thickness_m")
        continental = effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        continental_ledger.add_material_thickness(
            world,
            plate,
            np.maximum(hc_after - hc_before, 0.0),
            "juvenile_additions_m3",
            eligible=continental,
        )
        lava_m = np.maximum(plate.collect("elevation") - elevation_before, 0.0)
        sealed = lava_m / mobile_cover.VOLCANIC_SEAL_THICKNESS_M
        mobile_cover.end(world, plate, sealed, "volcanic_buried_m3")


def _splitmix64(x: np.ndarray) -> np.ndarray:
    """SplitMix64's finaliser, elementwise on uint64 (wrapping arithmetic)."""
    x = x + np.uint64(0x9E3779B97F4A7C15)
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def _node_uniforms(stream: tuple[int, ...], node_ids: np.ndarray) -> np.ndarray:
    """One uniform in [0, 1) per node, a pure function of `stream` and that node's own
    two-word surface ID (`SurfaceNodes.node_ids`) -- so a node's draw doesn't change when
    cells elsewhere on the plate are added or removed, which a single stream consumed in node
    order can't promise (issue #228 Phase 0a asked for exactly this of stable node IDs)."""
    salt = np.random.SeedSequence(list(stream)).generate_state(2, dtype=np.uint64)
    ids = np.asarray(node_ids, dtype=np.uint64).reshape(-1, 2)
    x = _splitmix64(_splitmix64(ids[:, 0] ^ salt[0]) ^ ids[:, 1] ^ salt[1])
    return (x >> np.uint64(11)).astype(float) * 2.0**-53


def _apply_volcanic_activity_to_surface(plate: Plate, world: "World", years: float) -> list[np.ndarray]:
    """Each active volcano's eruption roll, over the whole plate at once through the
    representation-neutral field API. Each node's draw is keyed by its stable node ID
    (`_node_uniforms`). Returns the world-space positions of every node that erupted this
    step, for `_spread_volcanic_plains` to spread an apron around."""
    is_volcano = plate.collect("is_volcano")
    if not np.any(is_volcano):
        return []
    remaining = plate.collect("volcano_active_years_remaining")
    active_mask = is_volcano & (remaining > 0)
    if not np.any(active_mask):
        return []

    active_years_this_step = np.minimum(years, remaining)
    p_erupt = 1.0 - np.exp(-ERUPTION_RATE_PER_MYR * world.volcanism_multiplier * active_years_this_step / 1_000_000.0)
    draws = _node_uniforms((world.seed, round(world.elapsed_years), plate.plate_id), plate.surface_nodes().node_ids)
    erupts = active_mask & (draws < p_erupt)

    elevation = plate.collect("elevation")
    new_crustal_thickness, new_elevation = lithosphere.back_elevation_gain_fields(
        elevation,
        plate.collect("crustal_thickness_m"),
        plate.collect("mantle_lithosphere_thickness_m"),
        plate.collect("crust_type_code"),
        plate.crust_type,
        ERUPTION_ELEVATION_M * world.volcanism_multiplier,
        erupts,
    )
    new_mineral_deposit = np.clip(
        plate.collect("mineral_deposit_m") + np.where(erupts, MINERAL_DEPOSIT_PER_ERUPTION_M, 0.0), 0.0, MAX_MINERAL_DEPOSIT_M
    )
    plate.set_fields_on_plate(
        elevation=new_elevation,
        crustal_thickness_m=new_crustal_thickness,
        volcano_active_years_remaining=np.clip(remaining - years, 0.0, None),
        mineral_deposit_m=new_mineral_deposit,
        elev_change_reason=np.where(erupts, ELEV_CHANGE_VOLCANO, plate.collect("elev_change_reason")),
    )
    world_points, _ = plate.all_points_and_elevation()
    return [world_points[erupts]] if np.any(erupts) else []


def _spread_volcanic_plains(plate: Plate, world: "World", years: float, erupted_points: list[np.ndarray]) -> None:
    """Spread a broad, low-relief apron around every vent that erupted this step -- a
    flood-basalt/shield-flank plain, distinct from the sharp point bump `_apply_volcanic_
    activity_to_surface` already applied there. Tapers linearly from
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

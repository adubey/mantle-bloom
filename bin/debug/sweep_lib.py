"""Shared plumbing for the parameter-sensitivity sweep (see run_sweep.py) -- issue: land
percent trending low / ice-cap percent trending high on long runs, suspected to trace back to
plate-collision behavior (GitHub issue #171). Kept separate from run_sweep.py so
aggregate_sweep.py can import PARAM_SPECS/CHECKPOINT_YEARS without pulling in the
multiprocessing orchestration.

Every job (one parameter, one multiplier-of-baseline, one seed) runs in a worker process from
a fresh `generate_world` -- never against one of the big `.mbworld` saves in ~/Downloads, which
are mid-history snapshots with whatever state a long play session left them in, not a
reproducible starting point a sweep can vary one knob against. `node_density=1.0` (the
production default, not debug_worlds.py's 0.5 fast-iteration density) so results reflect what a
real save actually does; a 10,000,000-year step size (one of the UI's own "years per step"
choices) keeps a 120 My run to 12 steps.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import biomes, mantle, render_image, stats  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app import lithosphere  # noqa: E402
from app.world import World, generate_world, step_world  # noqa: E402

# Five seeds pulled from real reported-problem saves under ~/Downloads (filenames encode
# seed and elapsed years), not arbitrary small integers -- so the sweep's baseline behavior is
# checked against the exact plate-generation seeds that produced the "land disappears" /
# "ice cap everywhere" reports, alongside a few other seeds that ran long without collapsing
# (for contrast). `generate_world(seed=...)` reproduces a seed's plate layout and mantle
# convection centers exactly; it does not, and cannot, reproduce a save's mid-history state --
# only ever the same starting point that produced it.
SWEEP_SEEDS = (829071382, 579428537, 896200538, 331006609, 673790085)

CHECKPOINT_YEARS = (30_000_000, 60_000_000, 90_000_000, 120_000_000)
STEP_YEARS = 10_000_000  # evenly divides every checkpoint above

NODE_DENSITY = 1.0
MULTIPLIERS = (0.25, 0.5, 1.0, 2.0, 4.0)

BASELINE_AVG_RATE_CM_YR = mantle.rad_per_yr_to_cm_per_yr(mantle.MANTLE_FLOW_REFERENCE_RATE)
BASELINE_MAX_RATE_CM_YR = mantle.rad_per_yr_to_cm_per_yr(mantle.MAX_PLATE_RATE)


@dataclass(frozen=True)
class ParamSpec:
    label: str
    unit: str
    baseline: float
    # Set for the two rotation-rate knobs (mantle.py module constants -- see
    # _apply_rotation_rate_overrides for why a plain module-attribute assignment isn't enough
    # on its own); `world_field` is set for the rest (a plain World dataclass field, see
    # world.TUNING_MULTIPLIER_FIELDS). Exactly one of the two is non-None.
    rotation_kind: str | None = None  # "avg" | "max"
    world_field: str | None = None


PARAM_SPECS: dict[str, ParamSpec] = {
    "avg_rotation_rate": ParamSpec(
        "Average continental rotation rate", "cm/yr", BASELINE_AVG_RATE_CM_YR, rotation_kind="avg"
    ),
    "max_rotation_rate": ParamSpec(
        "Maximum continental rotation rate", "cm/yr", BASELINE_MAX_RATE_CM_YR, rotation_kind="max"
    ),
    "rain_erosion": ParamSpec("Rain erosion", "x baseline", 1.0, world_field="rain_erosion_multiplier"),
    "river_erosion": ParamSpec("River erosion", "x baseline", 1.0, world_field="river_erosion_multiplier"),
    "glacier_erosion": ParamSpec("Glacier erosion", "x baseline", 1.0, world_field="glacier_erosion_multiplier"),
    "volcanism": ParamSpec("Volcanism", "x baseline", 1.0, world_field="volcanism_multiplier"),
    "collision_uplift_amount": ParamSpec(
        "Collision uplift amount", "x baseline", 1.0, world_field="collision_uplift_multiplier"
    ),
    "collision_uplift_distance": ParamSpec(
        "Collision uplift distance (reach)", "x baseline", 1.0, world_field="collision_uplift_reach_multiplier"
    ),
}

# The sentinel param name for an un-overridden run (every multiplier at 1.0) -- run once per
# seed and shared across every real parameter's own multiplier=1.0 point rather than re-run
# per parameter, since it's the same world either way. See run_sweep.build_jobs.
BASELINE_PARAM = "baseline"


def _apply_rotation_rate_overrides(avg_cm_yr: float, max_cm_yr: float) -> None:
    """Point every consumer of mantle.py's two rotation-rate constants at new values for the
    *next* generate_world/step_world call. A plain `mantle.MANTLE_FLOW_REFERENCE_RATE = x`
    module-attribute assignment alone is not enough: `generate_convection_centers`'s own
    `reference_rate` parameter and `clamp_rate`'s own `max_rate`/`min_rate` parameters default
    to these constants, and Python binds a function's default argument values once, at def
    time (module import) -- so a later attribute reassignment never reaches code that relies on
    that default. Reassigning `__defaults__` on both functions (alongside the module attributes
    themselves, which the handful of *direct* readers -- merge_split.py's collision speed
    fraction, plate_diagnostics.py/render_image.py's speed-normalized visuals -- do pick up)
    is what actually changes the rotation behavior; confirmed empirically (avg_rotation_rate up
    -> mean plate speed up; max_rotation_rate up -> more plates railing above the old cap).
    Always called unconditionally at the top of every job (even ones that don't touch rotation
    rate) so a worker process reused across jobs never leaks a previous job's override."""
    avg_rad = mantle.cm_per_yr_to_rad_per_yr(avg_cm_yr)
    max_rad = mantle.cm_per_yr_to_rad_per_yr(max_cm_yr)
    mantle.MANTLE_FLOW_REFERENCE_RATE = avg_rad
    mantle.MAX_PLATE_RATE = max_rad
    mantle.generate_convection_centers.__defaults__ = (8, avg_rad)
    mantle.clamp_rate.__defaults__ = (mantle.MIN_PLATE_RATE, max_rad)


def _land_volume_above_sea_km3(world: World) -> float:
    """Total continental+oceanic crust volume actually standing above sea level -- unlike
    stats.py's `total_continental_crust_volume_km3` (continental crustal thickness regardless of
    submergence), this is the "topographic excess" a collision pumping up a plateau adds to and
    erosion/land loss drains from, so it's the one stat here that can directly distinguish
    "land is disappearing" (this falls) from "land is fine but crust is just thicker/thinner
    underwater" (the crust-volume stat would move, this wouldn't). Same node-area-times-height
    approach as stats._total_land_area_and_continental_volume, restricted to the
    above-sea-level slice of elevation."""
    area_m2 = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    total_m3 = 0.0
    for plate in world.plates:
        _, elevation = plate.all_points_and_elevation()
        if len(elevation) == 0:
            continue
        above = elevation[elevation > world.sea_level_m] - world.sea_level_m
        if above.size:
            total_m3 += float(np.sum(above)) * area_m2
    return total_m3 / 1.0e9


def _plains_fraction_of_world(world: World) -> float | None:
    """Fraction of all elevation-line nodes (land + ocean, the same node cloud render_image's
    Elevation view classifies) whose *local relief* falls under render_image's own
    MOUNTAIN_RELIEF_THRESHOLD_M -- i.e. the exact "Plains & Plateaus" legend toggle's own
    definition, reused rather than reinvented so this stat means the same thing here as it does
    on screen. None for a degenerate empty-world edge case (never expected in practice)."""
    cloud = render_image._node_cloud_and_tree(world)
    if cloud is None:
        return None
    all_points, all_elevation, _all_owner, tree = cloud
    n = len(all_points)
    if n == 0:
        return None
    codes = render_image._classify_terrain_relief(all_points, all_elevation, tree, world.sea_level_m)
    return float(np.count_nonzero(codes == render_image._TERRAIN_PLAINS_PLATEAU)) / n


def compute_outcome_stats(world: World) -> dict:
    """The four outcome stats the sweep tracks, at one instant. `land_fraction` and
    `ice_cap_fraction` are both climate-grid-based (stats.compute_stats's own equirectangular
    resample, so they share a denominator with each other); `plains_fraction` is node-cloud-based
    (see `_plains_fraction_of_world`) since that's the only place the game already draws a
    mountain/plains distinction -- the two therefore aren't on quite the same underlying grid,
    same as the live app's own Stats panel vs. Elevation-view legend never are either."""
    snapshot = stats.compute_stats(world)
    land_fraction = snapshot["land_fraction"]
    # biomes.BIOME_NAMES holds display names ("Ice Cap"), not the Koppen letter code -- look the
    # key up via biomes.ICE_CAP rather than hardcoding "EF" (confirmed by hand this actually
    # matters: a first pass of this sweep silently read 0.0 for every one of 660 records because
    # of exactly this mismatch, even on a seed/checkpoint directly confirmed via
    # climate.compute_climate_cached to have hundreds of real Ice Cap cells).
    ice_share_of_land = snapshot["biome_land_fraction"].get(biomes.BIOME_NAMES[biomes.ICE_CAP], 0.0)
    return {
        "land_fraction": land_fraction,
        "ice_cap_fraction": land_fraction * ice_share_of_land,
        "plains_fraction": _plains_fraction_of_world(world),
        "land_volume_above_sea_km3": _land_volume_above_sea_km3(world),
        # Kept alongside for cross-checking the land-volume stat against the game's own
        # existing conservation-check metric (see stats.py's own module docstring, GitHub
        # issues #119/#120) -- not one of the four requested outcome stats itself.
        "total_continental_crust_volume_km3": snapshot["total_continental_crust_volume_km3"],
    }


def run_one_job(param_name: str, multiplier: float, seed: int) -> list[dict]:
    """Runs one (parameter, multiplier-of-baseline, seed) job from a fresh `generate_world`
    out to CHECKPOINT_YEARS' last entry, returning one outcome-stats record per checkpoint.
    `param_name == BASELINE_PARAM` runs with every knob at its untouched default (multiplier is
    still recorded, always 1.0) -- see PARAM_SPECS/BASELINE_PARAM's own docstrings for why this
    is shared across every real parameter's own multiplier=1.0 point rather than re-run per
    parameter."""
    _apply_rotation_rate_overrides(BASELINE_AVG_RATE_CM_YR, BASELINE_MAX_RATE_CM_YR)

    value = None
    spec = None if param_name == BASELINE_PARAM else PARAM_SPECS[param_name]
    if spec is not None:
        value = spec.baseline * multiplier
        if spec.rotation_kind == "avg":
            _apply_rotation_rate_overrides(value, BASELINE_MAX_RATE_CM_YR)
        elif spec.rotation_kind == "max":
            _apply_rotation_rate_overrides(BASELINE_AVG_RATE_CM_YR, value)

    world = generate_world(seed=seed, node_density=NODE_DENSITY)

    if spec is not None and spec.world_field is not None:
        setattr(world, spec.world_field, value)

    records = []
    years_done = 0.0
    for checkpoint in CHECKPOINT_YEARS:
        step_world(world, years=checkpoint - years_done)
        years_done = checkpoint
        record = {
            "parameter": param_name,
            "multiplier": multiplier,
            "value": value,
            "seed": seed,
            "checkpoint_years": checkpoint,
        }
        record.update(compute_outcome_stats(world))
        records.append(record)
    return records

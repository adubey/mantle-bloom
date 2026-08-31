"""Stranded-basin diagnostics: endorheic depressions whose floor sits below sea level and
that have no drainage path to the ocean at all.

This is the "land-locked coastal pit" pathology tracked in docs/TODO.md's coastal-speckle
section. An isolated sub-sea-level node ringed by land is neither hydrology's
connectivity-aware `is_ocean` nor above sea level, so the marine sink, coastal planation, and
lake infill all skip it -- it churns (merge/split) in the event log every step and never
drains or fills. The event log records this today but is drowned in near-sea-level
transient-pond spam (see docs/debugging.md), so this module surfaces the same information as
one clean list.

Two pieces:

- `find_stranded_basins` is a pure read of one already-resolved depression hierarchy
  (`lakes.build_lake_hierarchy`'s forest, exposed on `hydrology.HydrologyFields.lake_forest`).
  A top-level `Lake` with `max_depth is None` is exactly "no known spill to the ocean" (see
  lakes.py's own docstring); keep the ones whose floor is also below sea level.

- Persistence -- "how long has this basin been here" -- can't come from a single snapshot,
  since lakes.py deliberately keeps no cross-step `Lake` identity. `reconcile_world_tracks`
  (called once per step from `world.step_world`) matches this step's basins to last step's by
  centroid proximity and carries a first-seen timestamp forward on `world.stranded_basin_tracks`
  -- the same lightweight cross-step tracker `world.collision_progress` already is for plate
  pairs. `enrich_with_persistence` is the read-only side, for the server endpoint and the
  offline dump.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from . import geometry, hydrology, lakes

if TYPE_CHECKING:
    from .world import World

# A stranded basin drifts with its own plate at most ~MAX_PLATE_RATE (15 cm/yr) -> ~15 km per
# 100-ky step, far under this; two genuinely distinct stranded basins in the docs/TODO.md
# investigation were different depths and well separated. So a fixed ~0.08 rad (~500 km)
# centroid gate reliably re-identifies the same basin step to step without ever fusing two.
MATCH_DISTANCE_RAD = 0.08

# Convention only -- nothing in the engine hardcodes a step size -- but every run in
# docs/TODO.md and the UI's Step/Play buttons uses 100 ky, so reporting an approximate step
# count next to raw elapsed years is what makes "12.4 My (124 steps)" legible. Mirrors
# plate_diagnostics.CONVENTIONAL_YEARS_PER_STEP.
CONVENTIONAL_YEARS_PER_STEP = 100_000.0


@dataclass
class StrandedBasin:
    """One endorheic, below-sea-level, ocean-disconnected basin as of the step it's read from.
    `first_seen_years`/`persisted_years`/`steps_seen` are `None` straight out of
    `find_stranded_basins` and filled in by `enrich_with_persistence` (or
    `reconcile_world_tracks`)."""

    floor_elevation_m: float
    depth_below_sea_level_m: float  # sea_level_m - floor_elevation_m, always > 0
    catchment_node_count: int  # the full geometric catchment (Lake.members)
    flooded_node_count: int  # members currently holding visible standing water
    water_elevation_m: float | None  # current standing-water surface; None if bone dry
    centroid_xyz: tuple[float, float, float]  # unit vector, catchment centroid
    centroid_lat_deg: float
    centroid_lon_deg: float
    floor_xyz: tuple[float, float, float]  # the lowest member node's position
    first_seen_years: float | None = None
    persisted_years: float | None = None
    steps_seen: int | None = None


@dataclass
class StrandedBasinTrack:
    """One entry in `world.stranded_basin_tracks`: the minimal cross-step memory a stranded
    basin needs, reconciled by centroid proximity each step (see the module docstring)."""

    centroid_xyz: np.ndarray  # (3,) unit vector
    first_seen_years: float
    last_seen_years: float
    steps_seen: int


def find_stranded_basins(
    forest: list[lakes.Lake],
    elevation: np.ndarray,
    points: np.ndarray,
    lake_depth: np.ndarray,
    sea_level_m: float,
) -> list[StrandedBasin]:
    """Every top-level basin in `forest` that is endorheic (`max_depth is None` -- no path to
    the ocean at any fill level) *and* whose floor sits below `sea_level_m`. Roots only: an
    endorheic root is the maximal "no drainage" unit, and its floor is the min over every
    descendant, so a deep sub-basin below sea level is already covered by its root."""
    visible = hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    basins: list[StrandedBasin] = []
    for lake in forest:
        if lake.max_depth is not None:
            continue  # has a spill path to the ocean -- drains eventually, not stranded
        if lake.floor_elevation >= sea_level_m:
            continue  # endorheic but above sea level -- an ordinary closed basin, not the pit
        members = lake.members
        member_pts = points[members]
        centroid = geometry.normalize(member_pts.mean(axis=0))
        lat, lon = geometry.xyz_to_latlon(centroid)
        flooded = members[lake_depth[members] > visible]
        floor_node = int(members[np.argmin(elevation[members])])
        basins.append(
            StrandedBasin(
                floor_elevation_m=float(lake.floor_elevation),
                depth_below_sea_level_m=float(sea_level_m - lake.floor_elevation),
                catchment_node_count=int(len(members)),
                flooded_node_count=int(len(flooded)),
                water_elevation_m=float(lake.current_water_elevation) if len(flooded) > 0 else None,
                centroid_xyz=tuple(float(c) for c in centroid),
                centroid_lat_deg=round(float(np.degrees(lat)), 3),
                centroid_lon_deg=round(float(np.degrees(lon)), 3),
                floor_xyz=tuple(float(c) for c in points[floor_node]),
            )
        )
    # Deepest first -- the ones the coastal section cares about most sit at the top.
    basins.sort(key=lambda b: -b.depth_below_sea_level_m)
    return basins


def _nearest_track(
    prev_centroids: np.ndarray | None, used: np.ndarray | None, centroid: np.ndarray
) -> int | None:
    """Index of the closest not-yet-claimed previous track within `MATCH_DISTANCE_RAD` of
    `centroid`, or `None`."""
    if prev_centroids is None or len(prev_centroids) == 0:
        return None
    dist = geometry.angular_distance(prev_centroids, centroid)
    dist = np.where(used, np.inf, dist)
    j = int(np.argmin(dist))
    return j if dist[j] <= MATCH_DISTANCE_RAD else None


def enrich_with_persistence(
    basins: list[StrandedBasin], tracks: list[StrandedBasinTrack], elapsed_years: float
) -> list[StrandedBasin]:
    """Fill each basin's `first_seen_years`/`persisted_years`/`steps_seen` from whichever
    stored track its centroid matches -- read-only, no track mutation (the server endpoint and
    the offline dump both go through here against an at-most-one-step-stale
    `world.stranded_basin_tracks`, the same staleness tolerance `hydrology_cache` itself
    carries). A basin with no matching track is brand new this step: persisted 0, 1 step."""
    prev_centroids = np.array([t.centroid_xyz for t in tracks]) if tracks else None
    used = np.zeros(len(tracks), dtype=bool)
    for basin in basins:
        centroid = np.asarray(basin.centroid_xyz)
        match = _nearest_track(prev_centroids, used, centroid)
        if match is None:
            basin.first_seen_years = elapsed_years
            basin.persisted_years = 0.0
            basin.steps_seen = 1
        else:
            used[match] = True
            track = tracks[match]
            basin.first_seen_years = track.first_seen_years
            basin.persisted_years = elapsed_years - track.first_seen_years
            basin.steps_seen = track.steps_seen + 1
    return basins


def reconcile_world_tracks(world: "World") -> list[StrandedBasin]:
    """Recompute this step's stranded basins from `world.hydrology_cache` and replace
    `world.stranded_basin_tracks` with the reconciled set: matched basins keep their
    first-seen timestamp and bump `steps_seen`; unmatched basins start a fresh track; tracks
    with no basin this step are dropped. Called once per step from `world.step_world`, only on
    a step that actually recomputed hydrology. Returns the enriched basin list."""
    fields = world.hydrology_cache
    if fields is None or len(fields.points) == 0:
        world.stranded_basin_tracks = []
        return []
    basins = find_stranded_basins(
        fields.lake_forest, fields.elevation, fields.points, fields.lake_depth, world.sea_level_m
    )
    enrich_with_persistence(basins, world.stranded_basin_tracks, world.elapsed_years)
    world.stranded_basin_tracks = [
        StrandedBasinTrack(
            centroid_xyz=np.asarray(basin.centroid_xyz, dtype=float),
            first_seen_years=float(basin.first_seen_years),
            last_seen_years=world.elapsed_years,
            steps_seen=int(basin.steps_seen),
        )
        for basin in basins
    ]
    return basins


# ---------------------------------------------------------------------------
# Offline dump: python -m app.stranded_basins <save.mbworld> [--json]
# ---------------------------------------------------------------------------


def build_report(world: "World") -> dict:
    """Structured stranded-basin diagnostics for `world` -- the payload behind both the text
    dump and `--json`. Reads `world.hydrology_cache` (persisted in the save, same as the
    climate cache) and the persisted `world.stranded_basin_tracks`; a world never stepped
    with climate on has no hydrology snapshot and reports an empty list."""
    fields = getattr(world, "hydrology_cache", None)
    tracks = getattr(world, "stranded_basin_tracks", []) or []
    have_hydrology = fields is not None and len(fields.points) > 0
    basins: list[StrandedBasin] = []
    if have_hydrology:
        basins = find_stranded_basins(
            fields.lake_forest, fields.elevation, fields.points, fields.lake_depth, world.sea_level_m
        )
        enrich_with_persistence(basins, tracks, world.elapsed_years)
    return {
        "seed": world.seed,
        "elapsed_years": world.elapsed_years,
        "approx_steps": round(world.elapsed_years / CONVENTIONAL_YEARS_PER_STEP),
        "node_density": world.node_density,
        "sea_level_m": world.sea_level_m,
        "have_hydrology_snapshot": bool(have_hydrology),
        "stranded_basins": [
            {
                "floor_elevation_m": round(b.floor_elevation_m, 1),
                "depth_below_sea_level_m": round(b.depth_below_sea_level_m, 1),
                "catchment_node_count": b.catchment_node_count,
                "flooded_node_count": b.flooded_node_count,
                "water_elevation_m": None if b.water_elevation_m is None else round(b.water_elevation_m, 1),
                "centroid_lat_deg": b.centroid_lat_deg,
                "centroid_lon_deg": b.centroid_lon_deg,
                "centroid_xyz": [round(c, 6) for c in b.centroid_xyz],
                "floor_xyz": [round(c, 6) for c in b.floor_xyz],
                "first_seen_years": b.first_seen_years,
                "persisted_years": b.persisted_years,
                "steps_seen": b.steps_seen,
            }
            for b in basins
        ],
    }


def format_report(report: dict) -> str:
    lines: list[str] = []
    lines.append("mantle-bloom stranded-basin diagnostics")
    lines.append(f"  seed:          {report['seed']}")
    lines.append(
        f"  elapsed:       {report['elapsed_years']:,.0f} yr"
        f"  (~{report['approx_steps']:,} steps @ 100 ky)"
    )
    lines.append(f"  node_density:  {report['node_density']}")
    lines.append(f"  sea level:     {report['sea_level_m']:.1f} m")
    if not report["have_hydrology_snapshot"]:
        lines.append("")
        lines.append("  no hydrology snapshot in this save (never stepped with climate on) -- nothing to report")
        return "\n".join(line.rstrip() for line in lines)

    basins = report["stranded_basins"]
    lines.append(
        f"  stranded basins: {len(basins)}"
        f"   (endorheic, floor below sea level, no ocean drainage)"
    )
    lines.append("")
    if not basins:
        lines.append("  (none)")
        return "\n".join(line.rstrip() for line in lines)

    header = (
        f"  {'floor':>8} {'depth<SL':>9} {'catch':>7} {'flooded':>8} {'water':>8}"
        f"  {'centroid lat,lon':>18}  {'persisted':>20}"
    )
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))
    for b in basins:
        water = "--" if b["water_elevation_m"] is None else f"{b['water_elevation_m']:>8.0f}"
        pole = f"{b['centroid_lat_deg']:+7.1f},{b['centroid_lon_deg']:+8.1f}"
        persisted_my = (b["persisted_years"] or 0.0) / 1e6
        persisted = f"{persisted_my:.1f} My ({b['steps_seen']} steps)"
        lines.append(
            f"  {b['floor_elevation_m']:>8.0f} {b['depth_below_sea_level_m']:>9.0f}"
            f" {b['catchment_node_count']:>7} {b['flooded_node_count']:>8} {water}"
            f"  {pole:>18}  {persisted:>20}"
        )
    return "\n".join(line.rstrip() for line in lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m app.stranded_basins",
        description="List endorheic below-sea-level basins with no ocean drainage in a saved world.",
    )
    parser.add_argument("save", type=Path, help="path to a .mbworld save file")
    parser.add_argument("--json", action="store_true", help="emit the structured report as JSON")
    args = parser.parse_args(argv)

    if not args.save.is_file():
        parser.error(f"no such file: {args.save}")
    from . import persistence  # local: keeps world.py -> stranded_basins -> persistence -> world out of the import cycle

    world = persistence.load_world_bytes(args.save.read_bytes())
    report = build_report(world)

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(format_report(report))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

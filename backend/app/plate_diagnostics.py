"""Standalone plate-geometry diagnostics dump for a saved world.

    python -m app.plate_diagnostics path/to/save.mbworld
    python -m app.plate_diagnostics path/to/save.mbworld --json

Loads a `.mbworld` file (same pickle format as "File > Load", see persistence.py) and prints
the per-plate motion / shape / overlap table, the territory-overlap list, the
sustained-collision timers, and the total node count against a clean-tiling estimate -- the
exact set of numbers the plate-geometry investigation in docs/TODO.md keeps needing for its
"is this save's geometry healthy?" re-verify. It reuses `main._plate_summary` /
`main._plate_overlaps` (the same code path behind `GET /world/plates` and the Plate
Inspector), so the CLI and the UI never drift.

This is a read-only offline tool: it never starts the server or binds a port. See
docs/debugging.md for the full list of debug views and how to read this one.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from . import elevation_lines, mantle, merge_split, persistence
from .main import _plate_overlaps, _plate_summary
from .world import World

# Convention only -- nothing in the engine hardcodes a step size (see world.step_world, which
# advances by whatever `years` it's handed), but every run in docs/TODO.md and the UI's
# Step/Play buttons uses 100 ky, so reporting an approximate step count alongside the raw
# elapsed years is what makes "851 steps / 85.1 My" legible.
CONVENTIONAL_YEARS_PER_STEP = 100_000.0


def clean_tiling_node_estimate(node_density: float) -> float:
    """How many nodes a gap-free, non-overlapping lattice at this density would put on the
    whole sphere: unit-sphere area (4*pi) divided by the area one node's spacing covers
    (`line_spacing_rad**2`). The long-run node-count blowup in docs/TODO.md is measured as
    the ratio of the world's actual total to this number (~130k at node_density 4)."""
    spacing = elevation_lines.line_spacing_rad(node_density)
    return 4.0 * np.pi / (spacing * spacing)


def build_report(world: World) -> dict:
    """Structured diagnostics for `world` -- the payload behind both the text dump and
    `--json`. Plate rows are `main._plate_summary` trimmed to the scalar fields (the full
    version also carries every node coordinate, which a diagnostics dump has no use for)."""
    overlaps = _plate_overlaps(world)
    keep = (
        "plate_id", "crust_type", "num_points", "num_rows", "age_steps",
        "speed_cm_per_yr", "at_max_rate", "euler_pole", "median_elevation_m",
        "submerged_fraction", "overlaps", "collisions",
    )
    plate_rows = []
    for plate in world.plates:
        full = _plate_summary(plate, world, overlaps)
        plate_rows.append({key: full[key] for key in keep})

    total_nodes = sum(row["num_points"] for row in plate_rows)
    estimate = clean_tiling_node_estimate(world.node_density)
    collisions = sorted(
        ({"plates": [int(a), int(b)], "years": float(years)}
         for (a, b), years in world.collision_progress.items()),
        key=lambda entry: -entry["years"],
    )
    forced_merge_timers = sorted(
        ({"plates": [int(a), int(b)], "years": float(years)}
         for (a, b), years in getattr(world, "overlap_progress", {}).items()),
        key=lambda entry: -entry["years"],
    )
    return {
        "seed": world.seed,
        "elapsed_years": world.elapsed_years,
        "approx_steps": round(world.elapsed_years / CONVENTIONAL_YEARS_PER_STEP),
        "node_density": world.node_density,
        "sea_level_m": world.sea_level_m,
        "num_plates": len(world.plates),
        "max_plate_rate_cm_per_yr": mantle.rad_per_yr_to_cm_per_yr(mantle.MAX_PLATE_RATE),
        "plates": plate_rows,
        "collisions": collisions,
        "forced_merge_timers": forced_merge_timers,
        "node_budget": {
            "total_nodes": int(total_nodes),
            "clean_tiling_estimate": int(round(estimate)),
            "ratio": round(total_nodes / estimate, 3) if estimate else None,
        },
    }


def _fmt_pole(pole: dict | None) -> str:
    if pole is None:
        return "--"
    return f"{pole['lat_deg']:+6.1f},{pole['lon_deg']:+7.1f}"


# Overlap fractions below this are still real (the KDTree tolerance in main._plate_overlaps
# already excludes ordinary shared boundaries), but a fringe of a few contested nodes along
# a long boundary buries the handful of entries that matter in the text dump. The full,
# unfiltered list is always in `--json`.
_OVERLAP_DISPLAY_FLOOR = 0.005


def format_report(report: dict) -> str:
    lines: list[str] = []
    lines.append("mantle-bloom plate diagnostics")
    lines.append(f"  seed:          {report['seed']}")
    lines.append(
        f"  elapsed:       {report['elapsed_years']:,.0f} yr"
        f"  (~{report['approx_steps']:,} steps @ 100 ky)"
    )
    lines.append(f"  node_density:  {report['node_density']}")
    lines.append(f"  sea level:     {report['sea_level_m']:.1f} m")
    lines.append(f"  plates:        {report['num_plates']}")
    lines.append("")

    max_rate = report["max_plate_rate_cm_per_yr"]
    lines.append(
        f"per-plate  ( * = railed at MAX_PLATE_RATE {max_rate:.1f} cm/yr;"
        f"  ! = continental & >50% submerged )"
    )
    header = (
        f"  {'id':>3}  {'crust':<11} {'nodes':>7} {'rows':>5} {'age':>5}"
        f"  {'speed':>8}  {'pole lat,lon':>15}  {'med.elev':>9}  {'submrg':>6}"
    )
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))
    for row in sorted(report["plates"], key=lambda r: r["plate_id"]):
        railed = "*" if row["at_max_rate"] else " "
        drowned = (
            "!" if row["crust_type"] == "continental" and row["submerged_fraction"] > 0.5 else " "
        )
        rows_cell = "--" if row["num_rows"] is None else f"{row['num_rows']}"
        elev_cell = "--" if row["median_elevation_m"] is None else f"{row['median_elevation_m']:>9.0f}"
        lines.append(
            f"  {row['plate_id']:>3}  {row['crust_type']:<11} {row['num_points']:>7}"
            f" {rows_cell:>5} {row['age_steps']:>5}"
            f"  {row['speed_cm_per_yr']:>7.2f}{railed}  {_fmt_pole(row['euler_pole']):>15}"
            f"  {elev_cell}  {row['submerged_fraction']:>5.2f}{drowned}"
        )
    lines.append("")

    lines.append(
        f"territory overlaps  (share of THIS plate's nodes sitting on top of another;"
        f" >= {_OVERLAP_DISPLAY_FLOOR * 100:.1f}%, full list in --json)"
    )
    shown = 0
    for row in sorted(report["plates"], key=lambda r: r["plate_id"]):
        for over in row["overlaps"]:
            if over["fraction"] < _OVERLAP_DISPLAY_FLOOR:
                continue
            shown += 1
            since = over.get("since_years")
            since_cell = (
                f"   since {since / 1e6:.1f} My" if since is not None else ""
            )
            lines.append(
                f"  {row['plate_id']:>3} -> {over['plate_id']:<3}  {over['fraction'] * 100:5.1f}%{since_cell}"
            )
    if not shown:
        lines.append("  (none)")
    lines.append("")

    lines.append("sustained-collision timers  (world.collision_progress)")
    if report["collisions"]:
        for entry in report["collisions"]:
            a, b = entry["plates"]
            lines.append(f"  ({a:>3}, {b:>3})  {entry['years'] / 1e6:6.1f} My")
    else:
        lines.append("  (none)")
    lines.append("")

    lines.append(
        "forced-merge timers  (world.overlap_progress; a continental pair fuses at "
        f"{merge_split.FORCED_MERGE_SUSTAINED_YEARS / 1e6:.0f} My)"
    )
    if report.get("forced_merge_timers"):
        for entry in report["forced_merge_timers"]:
            a, b = entry["plates"]
            lines.append(f"  ({a:>3}, {b:>3})  {entry['years'] / 1e6:6.1f} My")
    else:
        lines.append("  (none)")
    lines.append("")

    budget = report["node_budget"]
    lines.append("node budget")
    lines.append(f"  total nodes:           {budget['total_nodes']:>9,}")
    lines.append(
        f"  clean-tiling estimate: {budget['clean_tiling_estimate']:>9,}"
        f"   (4*pi / line_spacing_rad({report['node_density']})^2)"
    )
    if budget["ratio"] is not None:
        over = (budget["ratio"] - 1.0) * 100.0
        if abs(over) < 0.5:
            over = 0.0
        lines.append(f"  ratio:                 {budget['ratio']:>9.2f}x   ({over:+.0f}%)")
    return "\n".join(line.rstrip() for line in lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m app.plate_diagnostics",
        description="Dump per-plate geometry / motion / overlap diagnostics for a saved world.",
    )
    parser.add_argument("save", type=Path, help="path to a .mbworld save file")
    parser.add_argument("--json", action="store_true", help="emit the structured report as JSON")
    args = parser.parse_args(argv)

    if not args.save.is_file():
        parser.error(f"no such file: {args.save}")
    world = persistence.load_world_bytes(args.save.read_bytes())
    report = build_report(world)

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(format_report(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())

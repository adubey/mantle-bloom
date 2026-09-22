"""Per-phase Hc/Hm budget report for GitHub issue #216 ("Investigate long-run Hc/Hm decline:
geological sinks, changing node counts, and numerical volume losses").

    python -m app.phase_budget_diagnostics --seed 611937962 --node-density 4 --years 1000000
    python -m app.phase_budget_diagnostics path/to/save.mbworld --years 1000000
    ... --json

Either replays a real `.mbworld` save forward by `--years` (this issue's own acceptance
criteria explicitly want the *current* save re-measured, not just a fresh world -- see the
issue's "Replay short intervals from this save" item) or steps a freshly generated world
(`--seed`/`--node-density`, no save given) the same way. Either way: turns on
`World.debug_diagnostics`, resets `World.phase_budget` so only this run's own interval is
measured, steps forward in ordinary `CONVENTIONAL_YEARS_PER_STEP`-sized increments (matching
plate_diagnostics.py's own convention), then prints/dumps the accumulated per-phase budget --
see phase_budget.py for exactly what each phase measures and where in the engine it's hooked
in. This is a read-only offline tool: it never starts the server or binds a port. See
docs/debugging.md for the full list of debug views and how to read this one.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import elevation_lines, lithosphere, persistence
from .phase_budget import SCOPES
from .world import World, generate_world, step_world

# Same convention plate_diagnostics.py uses -- nothing in the engine hardcodes a step size,
# but every long-run save/UI Step-Play cadence examined so far uses 100 ky.
CONVENTIONAL_YEARS_PER_STEP = 100_000.0


def build_report(world: World, total_years: float) -> dict:
    """Steps `world` forward by `total_years` (rounded to the nearest whole
    CONVENTIONAL_YEARS_PER_STEP) with instrumentation on, and returns the structured budget --
    the payload behind both the text dump and `--json`. Mutates `world` in place (advances
    it); pass an already-loaded/generated world you're fine stepping forward."""
    world.debug_diagnostics = True
    world.reset_phase_budget()
    node_area_m2 = lithosphere.node_area_m2(elevation_lines.line_spacing_rad(world.node_density))
    steps = max(1, round(total_years / CONVENTIONAL_YEARS_PER_STEP))
    start_years = world.elapsed_years

    for _ in range(steps):
        step_world(world, years=CONVENTIONAL_YEARS_PER_STEP)

    phases = []
    for phase, totals in world.phase_budget.items():
        scopes = {}
        for scope in SCOPES:
            s = totals["scopes"][scope]
            scopes[scope] = {
                **s,
                "area_before_m2": s["count_before"] * node_area_m2,
                "area_after_m2": s["count_after"] * node_area_m2,
                "delta_count": s["count_after"] - s["count_before"],
                "delta_sum_hc": s["sum_hc_after"] - s["sum_hc_before"],
                "delta_sum_hm": s["sum_hm_after"] - s["sum_hm_before"],
            }
        phases.append({"phase": phase, "calls": totals["calls"], "scopes": scopes})
    # Biggest net Hc mover first -- the ordering the issue's own "which mechanism actually
    # moves the needle" question cares about.
    phases.sort(key=lambda row: -abs(row["scopes"]["all"]["delta_sum_hc"]))

    return {
        "seed": world.seed,
        "node_density": world.node_density,
        "start_elapsed_years": start_years,
        "end_elapsed_years": world.elapsed_years,
        "steps": steps,
        "phases": phases,
    }


def format_report(report: dict) -> str:
    lines: list[str] = []
    lines.append("mantle-bloom Hc/Hm phase budget (GitHub issue #216)")
    lines.append(f"  seed:          {report['seed']}")
    lines.append(f"  node_density:  {report['node_density']}")
    lines.append(
        f"  interval:      {report['start_elapsed_years']:,.0f} -> {report['end_elapsed_years']:,.0f} yr"
        f"  ({report['steps']:,} steps @ {CONVENTIONAL_YEARS_PER_STEP:,.0f} yr)"
    )
    lines.append("")

    header = (
        f"  {'phase':<28} {'calls':>7}  {'d(count)':>10}  {'d(sum Hc) m':>16}  {'d(sum Hm) m':>16}"
    )
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))
    for row in report["phases"]:
        all_scope = row["scopes"]["all"]
        lines.append(
            f"  {row['phase']:<28} {row['calls']:>7,}  {all_scope['delta_count']:>10,}"
            f"  {all_scope['delta_sum_hc']:>16,.1f}  {all_scope['delta_sum_hm']:>16,.1f}"
        )
    lines.append("")

    lines.append("continental/oceanic node-type split (per phase, resolved against crust_type_code)")
    sub_header = f"  {'phase':<28} {'cont d(sum Hc)':>16}  {'ocean d(sum Hc)':>16}  {'cont d(sum Hm)':>16}  {'ocean d(sum Hm)':>16}"
    lines.append(sub_header)
    lines.append("  " + "-" * (len(sub_header) - 2))
    for row in report["phases"]:
        cont = row["scopes"]["continental_node"]
        ocean = row["scopes"]["oceanic_node"]
        lines.append(
            f"  {row['phase']:<28} {cont['delta_sum_hc']:>16,.1f}  {ocean['delta_sum_hc']:>16,.1f}"
            f"  {cont['delta_sum_hm']:>16,.1f}  {ocean['delta_sum_hm']:>16,.1f}"
        )
    return "\n".join(line.rstrip() for line in lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m app.phase_budget_diagnostics",
        description="Replay a world forward and report its per-phase Hc/Hm budget (GitHub issue #216).",
    )
    parser.add_argument("save", type=Path, nargs="?", help="path to a .mbworld save file (omit to generate a fresh world)")
    parser.add_argument("--seed", type=int, default=0, help="seed for a freshly generated world (ignored if `save` is given)")
    parser.add_argument("--node-density", type=float, default=4.0, help="node density for a freshly generated world (ignored if `save` is given)")
    parser.add_argument("--years", type=float, required=True, help="how many years to step forward before reporting")
    parser.add_argument("--json", action="store_true", help="emit the structured report as JSON")
    args = parser.parse_args(argv)

    if args.save is not None:
        if not args.save.is_file():
            parser.error(f"no such file: {args.save}")
        world = persistence.load_world_bytes(args.save.read_bytes())
    else:
        world = generate_world(seed=args.seed, node_density=args.node_density)

    report = build_report(world, args.years)

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(format_report(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())

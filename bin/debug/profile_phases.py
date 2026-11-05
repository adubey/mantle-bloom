"""Where does the simulation spend its wall-clock time, in user-facing phase names? (#801)

Wraps a table of functions (`PHASES` below) with a nesting-aware wall-clock tracer, then runs
three scenarios and prints a phase tree for each:

  * `generate` -- "Building world": `generate_world`, run `--generate-runs` times in one
    process. Run 1 pays numba JIT compilation, later runs don't, so the gap between them is the
    cold-start cost a freshly started server charges its first user.
  * `render`   -- the first `combined` map render after a step (every render-side cache is
    cold), then a second render (warm). The difference is the cache-fill cost.
  * `step`     -- `--steps` calls of `step_world` after `--warmup` unmeasured ones.

Unlike cProfile, whose cumulative times overlap and double-count nested calls, every tree row
here carries an *exclusive* time (own time minus time in traced children), so a parent's rows
add up to its total and the leftover is shown as `(untraced)`. A large untraced share means the
table is missing a hook, not that the time is free: drill in with `--cprofile OUT.prof`
(step scenario only) and add the hot function to `PHASES`. Hooks that never fire are listed so
a stale entry (a renamed function) gets noticed instead of silently reading as zero.

A phase wrapped inside the same label (recursion, or two hooks sharing one label) is merged
into the outer one rather than nested.

Usage (from the repository root, with this worktree's own venv):

    backend/.venv/bin/python bin/debug/profile_phases.py --seed 0 --steps 3
    backend/.venv/bin/python bin/debug/profile_phases.py --load ~/Downloads/<save>.mbworld \
        --scenarios render,step --steps 3 --json /path/outside/repo/out.json
    backend/.venv/bin/python bin/debug/profile_phases.py --scenarios step --cprofile /tmp/x.prof

Options: --seed, --node-density/--climate-density/--fluid-density (UI defaults 4/4/2),
--step-years (default 100000), --scenarios (comma list of generate,render,step), --load (use a
saved world for render/step instead of the generated one), --json (write the trees as JSON),
--min-pct (hide rows under this share of their scenario, default 1.0), --synthetic-relief
(generate with the synthetic history), --wind-model (diagnostic|cfd; the UI default is diagnostic, which skips the CFD solve entirely).

Follow the compute-budget skill first: a mature save peaks near 6 GB. Timings are wall-clock,
so run on an idle machine and compare runs only like-for-like (same machine, load, inputs).
"""

from __future__ import annotations

import argparse
import contextlib
import cProfile
import importlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import trace  # noqa: E402 -- needs the path above

# (group, label, "module:attribute") -- attribute may be dotted for a method ("Class.method").
# `group` is "<bar>/<group>", the user-visible hierarchy from analysis/issue801/traces.md. The
# tracer nests rows by the real call stack, but the by-group rollup uses each label's own
# group, so Hydrometeorology rows count toward it even though they run inside erosion. A
# label with no deeper rows is a leaf; deeper rows are profiling-only detail. Labels must be
# unique: a hook inside another hook with the same label is merged into the outer one.
LAYOUT = "Building world/Laying out the planet"
SEEDING = "Building world/Seeding the surface"
CACHE = "Building cache/Render cache"
PLATES = "Stepping/Moving plates"
HYDRO = "Stepping/Hydrometeorology (climate, water & ice)"
SURFACE = "Stepping/Reshaping the surface"
PHASES: list[tuple[str, str, str]] = [
    # --- Building world ---
    (LAYOUT, "Placing sutures and cratons", "app.lithosphere_plate:cratons_first_layout"),
    (LAYOUT, "Shaping continents", "app.continent_layout:_continents"),
    (LAYOUT, "Raising orogenic belts", "app.continent_layout:_belts"),
    (LAYOUT, "Seeding Voronoi plate partitions", "app.plate_partition:_polygon_seed_partition"),
    (LAYOUT, "Raising initial orogens", "app.initial_terrain:apply_predicted_orogens"),
    (LAYOUT, "Shaping lowland relief", "app.initial_terrain:apply_lowland_spherical_relief"),
    (SEEDING, "Seeding atmosphere and ocean state", "app.hydroclimate.atmosphere_cfd:init_atmosphere_cfd"),
    (SEEDING, "Recording stats", "app.world:World.record_stats"),
    # --- Building cache (first map render after a step) ---
    (CACHE, "k-d tree caches", "app.render_image:_node_cloud_and_tree"),
    (CACHE, "Resampling climate onto the map grid", "app.render_image:_project_climate_grid"),
    (CACHE, "Looking up nearest nodes for every map pixel", "app.render_image:_biome_fields"),
    (CACHE, "Classifying biomes", "app.render_image:_classify_biome_grid"),
    (CACHE, "Drawing rivers", "app.render_image:_draw_rivers"),
    (CACHE, "Encoding the image", "app.render_image:_encode_image"),
    # --- Stepping: Moving plates ---
    (PLATES, "Coupled plate solve", "app.torque:solve_plates_quasistatic"),
    (PLATES, "Gathering boundary forces", "app.torque:gather_boundary_force_inputs"),
    (PLATES, "Resolving blocked contacts", "app.torque:_blocked_contact_response"),
    (PLATES, "Growing plate frontiers", "app.quad_tectonics:grow_frontier"),
    (PLATES, "Filling plate gaps (per plate)", "app.quad_tectonics:fill_gap"),
    (PLATES, "Reading plate boundary context", "app.lithosphere_plate:boundary_context"),
    (PLATES, "Looking up boundary neighbours", "app.torque:BoundarySearchCache.neighbours"),
    (PLATES, "Finding neighbouring plates", "app.sparse_quad_patch:PlateWithSparseQuadPatch.get_neighbours"),
    (PLATES, "Outlining plate territory", "app.sparse_quad_patch:PlateWithSparseQuadPatch.outline_world"),
    (PLATES, "Generating boundary faults", "app.faults:generate_boundary_faults"),
    (PLATES, "Shearing plates along faults", "app.faults:_apply_plate_fault_shear"),
    (PLATES, "Applying fault relief", "app.faults:_apply_plate_fault_relief"),
    (PLATES, "Generating earthquakes", "app.faults:_generate_earthquakes"),
    (PLATES, "Triggering ice-unloading earthquakes", "app.faults:trigger_unloading_earthquakes"),
    (PLATES, "Finding gap points", "app.gaps:_find_gap_points"),
    (PLATES, "Growing neighbours into gaps", "app.gaps:_grow_adjacent_into_gap"),
    # --- Stepping: Hydrometeorology ---
    (HYDRO, "Advancing the wind and ocean solver", "app.world:_advance_fluid_dynamics"),
    (HYDRO, "Refreshing CFD forcing", "app.hydroclimate.atmosphere_cfd:refresh_forcing"),
    (HYDRO, "Finding the humidity distribution", "app.hydroclimate.climate:compute_humidity"),
    (HYDRO, "Finding the rainfall distribution", "app.hydroclimate.climate:compute_precipitation"),
    (HYDRO, "Computing ocean currents", "app.hydroclimate.climate:compute_ocean_currents"),
    (HYDRO, "Ordering the flow network", "app.hydroclimate.hydrology:_routing_order"),
    (HYDRO, "Routing flow downstream", "app.hydroclimate.hydrology:route_downstream"),
    (HYDRO, "Routing external water sources", "app.hydroclimate.hydrology:_route_external_water_sources"),
    (HYDRO, "Resolving lakes", "app.hydroclimate.lakes:solve_lakes_prognostic"),
    (HYDRO, "Solving the plastic ice body", "app.hydroclimate.ice_sheet:_plastic_body"),
    (HYDRO, "Shaping ice sheets", "app.hydroclimate.ice_sheet:_shape"),
    (HYDRO, "Tracking stranded basins", "app.hydroclimate.stranded_basins:reconcile_world_tracks"),
    # --- Stepping: Reshaping the surface ---
    (SURFACE, "Computing slopes", "app.hydroclimate.erosion:compute_slope"),
    (SURFACE, "Applying ice load", "app.hydroclimate.erosion:_apply_ice_load"),
    (SURFACE, "Routing mass wasting", "app.hydroclimate.erosion:_route_mass_wasting"),
    (SURFACE, "Depositing wind-blown sediment", "app.hydroclimate.erosion:_route_wind_deposit"),
    (SURFACE, "Spreading beach sediment", "app.hydroclimate.erosion:_spread_beach_sediment"),
    (SURFACE, "Spreading marine sediment", "app.hydroclimate.erosion:_spread_marine_sediment"),
    (SURFACE, "Spreading lake sediment", "app.hydroclimate.erosion:_spread_lake_sediment_capped"),
    (SURFACE, "Leveling coasts", "app.hydroclimate.erosion:_spread_coastal_leveling"),
    (SURFACE, "Updating cratons", "app.cratons:update"),
    (SURFACE, "Solving connected-ocean sea level", "app.hydroclimate.eustasy:_solve_sea_level_connected"),
]
LABEL_GROUP = {label: group for group, label, _spec in PHASES}

# Anything the table can't resolve (a refactor moved it) is reported, never fatal.


def _import(module_name: str):
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError:
        return None


def install() -> list[str]:
    """Patch every resolvable `PHASES` entry in place with `trace.wrap`; return the specs that
    didn't resolve. The functions the app already traces natively (the user-visible phases) are
    not in the table, so nothing is counted twice."""
    missing = []
    for group, label, spec in PHASES:
        module_name, attr = spec.split(":")
        module = _import(module_name)
        if module is None:
            missing.append(spec)
            continue
        owner = module
        *parents, name = attr.split(".")
        try:
            for parent in parents:
                owner = getattr(owner, parent)
            original = getattr(owner, name)
        except AttributeError:
            missing.append(spec)
            continue
        setattr(owner, name, trace.wrap(original, label, group))
    return missing


def _children(stats, path):
    depth = len(path)
    return sorted(
        (p for p in stats if len(p) == depth + 1 and p[:depth] == path),
        key=lambda p: -stats[p][0],
    )


def render_tree(title: str, total: float, stats, divisor: float, min_pct: float) -> list[str]:
    """Format `stats` (path -> [inclusive, exclusive, calls]) as an indented table, per-unit."""
    lines = [f"\n== {title}: {total / divisor:.2f} s ==", f"{'phase':<62}{'incl s':>8}{'own s':>8}{'%':>6}{'calls':>7}"]
    roots = [p for p in stats if len(p) == 1]
    traced_root = sum(stats[p][0] for p in roots)

    def emit(path, indent):
        incl, own, calls = stats[path]
        pct = 100.0 * incl / total if total else 0.0
        if pct >= min_pct:
            name = ("  " * indent) + path[-1]
            lines.append(f"{name:<62}{incl / divisor:>8.2f}{own / divisor:>8.2f}{pct:>6.1f}{calls / divisor:>7.1f}")
            for child in _children(stats, path):
                emit(child, indent + 1)

    for root in sorted(roots, key=lambda p: -stats[p][0]):
        emit(root, 0)
    untraced = max(total - traced_root, 0.0)
    lines.append(f"{'(untraced: not under any hook)':<62}{untraced / divisor:>8.2f}{'':>8}{100.0 * untraced / total if total else 0:>6.1f}")
    return lines


def group_rollup(total: float, stats, divisor: float, groups: dict) -> list[str]:
    """Exclusive seconds per user-visible group, by each label's own group (not its call parent)."""
    by_group: dict[str, float] = {}
    for path, (_incl, own, _calls) in stats.items():
        group = groups.get(path[-1], "?")
        by_group[group] = by_group.get(group, 0.0) + own
    traced = sum(by_group.values())
    lines = ["by group (exclusive):"]
    for group, seconds in sorted(by_group.items(), key=lambda kv: -kv[1]):
        lines.append(f"  {group:<58}{seconds / divisor:>7.2f}{100.0 * seconds / total if total else 0:>6.1f}%")
    lines.append(f"  {'(untraced)':<58}{max(total - traced, 0.0) / divisor:>7.2f}{100.0 * max(total - traced, 0.0) / total if total else 0:>6.1f}%")
    return lines


def label_totals(total: float, stats, divisor: float, limit: int = 12) -> list[str]:
    """Exclusive seconds per label summed over every parent it ran under: a function reached
    from several phases (boundary forces, frontier growth) is split across the tree above."""
    by_label: dict[str, tuple[float, float]] = {}
    for path, (_incl, own, calls) in stats.items():
        seconds, count = by_label.get(path[-1], (0.0, 0.0))
        by_label[path[-1]] = (seconds + own, count + calls)
    lines = ["top labels, own time summed over parents:"]
    for label, (seconds, count) in sorted(by_label.items(), key=lambda kv: -kv[1][0])[:limit]:
        lines.append(f"  {label:<58}{seconds / divisor:>7.2f}{100.0 * seconds / total if total else 0:>6.1f}%{count / divisor:>8.1f} calls")
    return lines


def tree_to_json(stats, divisor: float) -> list[dict]:
    return [
        {"path": list(path), "inclusive_s": row[0] / divisor, "exclusive_s": row[1] / divisor, "calls": row[2] / divisor}
        for path, row in sorted(stats.items())
    ]


def main() -> None:  # noqa: C901
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--node-density", type=float, default=4.0)
    ap.add_argument("--climate-density", type=float, default=4.0)
    ap.add_argument("--fluid-density", type=float, default=2.0)
    ap.add_argument("--step-years", type=float, default=100_000.0)
    ap.add_argument("--steps", type=int, default=3)
    ap.add_argument("--warmup", type=int, default=2)
    ap.add_argument("--generate-runs", type=int, default=2)
    ap.add_argument("--scenarios", default="generate,render,step")
    ap.add_argument("--load", type=Path)
    ap.add_argument("--cprofile", type=Path)
    ap.add_argument("--json", type=Path)
    ap.add_argument("--min-pct", type=float, default=1.0)
    ap.add_argument("--synthetic-relief", action="store_true", help="generate with the synthetic tectonic history (#218), off by default in the app")
    ap.add_argument("--wind-model", choices=("diagnostic", "cfd"), help="override World.wind_model for render/step (UI default: diagnostic)")
    args = ap.parse_args()
    scenarios = set(args.scenarios.split(","))

    from app import persistence, render_image, world as world_mod

    rec = trace.Recorder()
    recording = contextlib.ExitStack()
    recording.enter_context(trace.recording(rec))
    missing = install()
    results: dict[str, object] = {"args": {k: str(v) for k, v in vars(args).items()}}
    report: list[str] = []
    world = None

    def snapshot_scenario(name: str, total: float, divisor: float = 1.0) -> None:
        stats = rec.snapshot()
        report.extend(render_tree(name, total, stats, divisor, args.min_pct))
        report.extend(group_rollup(total, stats, divisor, {**LABEL_GROUP, **rec.groups}))
        report.extend(label_totals(total, stats, divisor))
        results[name] = {"total_s": total / divisor, "tree": tree_to_json(stats, divisor)}
        rec.reset()

    gen_kwargs = dict(node_density=args.node_density, climate_density=args.climate_density, fluid_density=args.fluid_density)
    if args.synthetic_relief:
        gen_kwargs["synthetic_relief"] = True

    if "generate" in scenarios:
        for run in range(args.generate_runs):
            rec.reset()
            start = time.perf_counter()
            world = world_mod.generate_world(args.seed, **gen_kwargs)
            kind = "cold (includes numba JIT)" if run == 0 else "warm"
            snapshot_scenario(f"Building world, run {run + 1} {kind}", time.perf_counter() - start)
    if args.load:
        with args.load.open("rb") as handle:
            world = persistence.load_world_file(handle)
    elif world is None and scenarios & {"render", "step"}:
        world = world_mod.generate_world(args.seed, **gen_kwargs)
    if world is not None and args.wind_model:
        world.wind_model = args.wind_model
    if world is not None:
        nodes = sum(len(plate.all_points_and_elevation()[1]) for plate in world.plates)
        results["world"] = {"nodes": nodes, "plates": len(world.plates), "elapsed_years": world.elapsed_years, "steps_taken": world.steps_taken}
        print(f"world: {nodes} nodes, {len(world.plates)} plates, {world.elapsed_years / 1e6:.1f} Myr, wind_model={world.wind_model}", file=sys.stderr)

    if "render" in scenarios:
        world_mod.step_world(world, args.step_years)  # caches are dropped by a step, as in the UI
        rec.reset()
        for label in ("cold (building cache)", "warm"):
            start = time.perf_counter()
            render_image.render_png(world, "eckert4", "combined", 1100, 611)
            snapshot_scenario(f"First map render after a step, {label}", time.perf_counter() - start)

    if "step" in scenarios:
        for _ in range(args.warmup):
            world_mod.step_world(world, args.step_years)
        rec.reset()
        profiler = cProfile.Profile() if args.cprofile else None
        per_step = []
        if profiler:
            profiler.enable()
        for _ in range(args.steps):
            start = time.perf_counter()
            world_mod.step_world(world, args.step_years)
            per_step.append(time.perf_counter() - start)
        if profiler:
            profiler.disable()
            profiler.dump_stats(str(args.cprofile))
        results["step_seconds"] = per_step
        snapshot_scenario(f"Stepping, mean of {args.steps} steps (per step)", sum(per_step), divisor=float(args.steps))
        report.append("per-step wall: " + ", ".join(f"{s:.2f}" for s in per_step))

    print("\n".join(report))
    unfired = sorted(
        {spec for _g, label, spec in PHASES if label not in rec.fired and spec not in missing}
    )
    if missing:
        print("\nHooks that did not resolve (renamed or moved?):\n  " + "\n  ".join(missing))
    if unfired:
        print("\nHooks that never fired in the scenarios run:\n  " + "\n  ".join(unfired))
    if args.json:
        args.json.write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()

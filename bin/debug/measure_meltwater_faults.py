#!/usr/bin/env python3
"""Issue #275 phase 5: multi-seed long runs of ice-unloading fault triggering, against a baseline.

Generates a world per seed, steps it, and every --every steps records:

- land fraction (area-weighted) and `hc_at_max_fraction` (the node-count metric `/world/stats`
  reports);
- craton survival: live craton volume and area, as a share of the step-0 values, plus the
  craton ledger's balance;
- every continental-material ledger reservoir (`continental_ledger.inventories`);
- earthquakes born since the last snapshot, split by trigger (a backend without the `trigger`
  field counts every quake as tectonic), and the largest unloading-triggered magnitude;
- the most ice any column lost in one step since the last snapshot, as a stress (MPa).

Long runs diverge chaotically once any surface field differs, so `--controlled` isolates the
coupling instead: at every --every steps of one run, it deep-copies the world and steps each copy
--controlled-steps more with and without `faults.trigger_unloading_earthquakes`, then reports the
two copies' erosion budgets and triggered-quake counts side by side (this backend only).

Seeds run in parallel. Runs against any backend, so the same script gives a `main` baseline:

    backend/.venv/bin/python bin/debug/measure_meltwater_faults.py --seeds 1 2 3 --steps 100 --out /tmp/branch.jsonl
    backend/.venv/bin/python bin/debug/measure_meltwater_faults.py --backend ../mantle-bloom-origin-main/backend \\
        --seeds 1 2 3 --steps 100 --out /tmp/main.jsonl
    backend/.venv/bin/python bin/debug/measure_meltwater_faults.py --compare /tmp/main.jsonl /tmp/branch.jsonl
    backend/.venv/bin/python bin/debug/measure_meltwater_faults.py --controlled --seeds 1 2 3 --steps 80 --every 20
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

LEDGER_KEYS = (
    "surface_continental_derived_m3",
    "continental_sediment_on_oceanic_hosts_m3",
    "delaminated_lower_crust_m3",
    "deeply_subducted_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
    "juvenile_additions_m3",
    "balance_error_m3",
)
COMPARE_KEYS = (
    "land_fraction",
    "hc_at_max_fraction",
    "craton_volume_share",
    "craton_area_share",
    "quakes_tectonic",
    "quakes_ice_unloading",
    "max_triggered_mw",
    *(f"ledger_{k}" for k in LEDGER_KEYS),
)


def _run_seed(backend: str, seed: int, steps: int, years: float, every: int) -> list[dict]:
    sys.path.insert(0, backend)
    from app import continental_ledger, cratons, erosion  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.plates import collect_all_accounting_areas_m2, collect_all_crustal_thickness  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    cap = __import__("app.lithosphere", fromlist=["x"]).MAX_CRUSTAL_THICKNESS_M
    window = {"max_unload_mpa": 0.0}
    # A reused pool worker already wrapped it for an earlier seed.
    original_apply = getattr(erosion.apply_erosion, "__wrapped__", erosion.apply_erosion)

    def wrapped(world, years, node_cloud=None):
        result = original_apply(world, years, node_cloud=node_cloud)
        change = getattr(result, "ice_load_change_pa", None) if result is not None else None
        if change is not None and len(change):
            window["max_unload_mpa"] = max(window["max_unload_mpa"], float(-np.min(change)) / 1e6)
        return result

    wrapped.__wrapped__ = original_apply
    erosion.apply_erosion = wrapped

    world = generate_world(seed=seed)
    craton0 = cratons.diagnostics(world)
    vol0 = max(craton0["craton_volume_m3"], 1.0)
    area0 = max(craton0["craton_area_m2"], 1.0)
    rows = []
    seen_ids: set[int] = set()
    quakes = {"tectonic": 0, "ice_unloading": 0}
    max_triggered_mw = 0.0
    start = time.perf_counter()
    for step in range(1, steps + 1):
        step_world(world, years)
        for q in world.earthquakes:
            if q.earthquake_id in seen_ids:
                continue
            seen_ids.add(q.earthquake_id)
            trigger = getattr(q, "trigger", "tectonic")
            quakes[trigger] = quakes.get(trigger, 0) + 1
            if trigger == "ice_unloading":
                max_triggered_mw = max(max_triggered_mw, float(q.magnitude))
        if step % every and step != steps:
            continue
        hydro = world.hydrology_cache
        order = hydro.plates_in_order
        hc = collect_all_crustal_thickness(order)
        areas = collect_all_accounting_areas_m2(order, line_spacing_rad(world.node_density))
        land = ~hydro.is_ocean
        inv = continental_ledger.inventories(world)
        craton = cratons.diagnostics(world)
        row = {
            "seed": seed,
            "step": step,
            "seconds": round(time.perf_counter() - start, 1),
            "land_fraction": float(np.sum(areas[land]) / np.sum(areas)),
            "hc_at_max_fraction": float(np.mean(hc >= cap - 1e-6)),
            "craton_volume_share": craton["craton_volume_m3"] / vol0,
            "craton_area_share": craton["craton_area_m2"] / area0,
            "craton_balance_error_m3": craton["craton_balance_error_m3"],
            "quakes_tectonic": quakes.get("tectonic", 0),
            "quakes_ice_unloading": quakes.get("ice_unloading", 0),
            "max_triggered_mw": max_triggered_mw,
            "max_unload_mpa": window["max_unload_mpa"],
            **{f"ledger_{k}": float(inv[k]) for k in LEDGER_KEYS if k in inv},
        }
        quakes = {"tectonic": 0, "ice_unloading": 0}
        max_triggered_mw = 0.0
        window["max_unload_mpa"] = 0.0
        rows.append(row)
        print(json.dumps({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    return rows


def _run_controlled(backend: str, seed: int, steps: int, years: float, every: int, extra: int) -> list[dict]:
    sys.path.insert(0, backend)
    from app import erosion, faults  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    real_trigger = faults.trigger_unloading_earthquakes
    budgets: list[dict] = []
    original_apply = getattr(erosion.apply_erosion, "__wrapped__", erosion.apply_erosion)

    def wrapped(world, years, node_cloud=None):
        result = original_apply(world, years, node_cloud=node_cloud)
        if result is not None:
            budgets.append(dict(result.budget or {}))
        return result

    wrapped.__wrapped__ = original_apply
    erosion.apply_erosion = wrapped

    def triggered_ids(world) -> set[int]:
        return {q.earthquake_id for q in world.earthquakes if getattr(q, "trigger", "tectonic") == "ice_unloading"}

    world = generate_world(seed=seed)
    rows = []
    for step in range(1, steps + 1):
        step_world(world, years)
        if step % every:
            continue
        row: dict = {"seed": seed, "step": step}
        for label, trigger in (("off", lambda *args: []), ("on", real_trigger)):
            branch = copy.deepcopy(world)
            before = triggered_ids(branch)
            faults.trigger_unloading_earthquakes = trigger
            budgets.clear()
            try:
                for _ in range(extra):
                    step_world(branch, years)
            finally:
                faults.trigger_unloading_earthquakes = real_trigger
            # New ids only; older triggered quakes may have aged out meanwhile.
            row[f"ice_quakes_{label}"] = len(triggered_ids(branch) - before)
            for key in ("mass_wasting_removed_m3", "removed_m3", "hc_cap_overflow_m3"):
                row[f"{key}_{label}"] = sum(b.get(key, 0.0) for b in budgets)
        rows.append(row)
        print(json.dumps({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    return rows


def _compare(base_path: Path, branch_path: Path) -> None:
    def load(path: Path) -> dict[tuple[int, int], dict]:
        rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        return {(r["seed"], r["step"]): r for r in rows}

    base, branch = load(base_path), load(branch_path)
    keys = sorted(set(base) & set(branch))
    seeds = sorted({s for s, _ in keys})
    for key in COMPARE_KEYS:
        print(f"\n{key}")
        for seed in seeds:
            steps = [st for s, st in keys if s == seed]
            cells = []
            for st in steps:
                a, b = base[(seed, st)].get(key, float("nan")), branch[(seed, st)].get(key, float("nan"))
                cells.append(f"{st}: {a:.4g} -> {b:.4g}")
            print(f"  seed {seed}  " + " | ".join(cells))
        # Window sums for quake counts; everything else is a snapshot, compare the windowed mean.
        late = [k for k in keys if k[1] > max(st for _, st in keys) // 2]
        if key == "max_triggered_mw":
            agg_a = max(base[k].get(key, 0.0) for k in keys)
            agg_b = max(branch[k].get(key, 0.0) for k in keys)
            print(f"  largest over all seeds/steps: {agg_a:.3g} -> {agg_b:.3g}")
        elif key.startswith("quakes_"):
            agg_a = sum(base[k].get(key, 0) for k in keys)
            agg_b = sum(branch[k].get(key, 0) for k in keys)
            print(f"  total over all seeds/steps: {agg_a} -> {agg_b}")
        elif late:
            mean_a = np.mean([base[k].get(key, np.nan) for k in late])
            mean_b = np.mean([branch[k].get(key, np.nan) for k in late])
            print(f"  mean over later half, all seeds: {mean_a:.4g} -> {mean_b:.4g}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--every", type=int, default=10, help="snapshot every Nth step")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--out", type=Path, help="JSON-lines output")
    parser.add_argument("--compare", type=Path, nargs=2, metavar=("BASE", "BRANCH"), help="compare two --out files")
    parser.add_argument("--controlled", action="store_true", help="on/off steps from identical states (see above)")
    parser.add_argument("--controlled-steps", type=int, default=3)
    args = parser.parse_args()

    if args.compare:
        _compare(*args.compare)
        return
    backend = str(args.backend.resolve())
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        if args.controlled:
            futures = [
                pool.submit(_run_controlled, backend, s, args.steps, args.years, args.every, args.controlled_steps)
                for s in args.seeds
            ]
        else:
            futures = [pool.submit(_run_seed, backend, s, args.steps, args.years, args.every) for s in args.seeds]
        rows = [row for f in futures for row in f.result()]
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text("".join(json.dumps(r) + "\n" for r in rows))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Issue #867: the river-substep assumption ladder, as a reference trajectory and branches.

Two subcommands share one measurement loop:

- `reference` generates a world (or loads a save) and steps it with the production model.
  At every checkpoint it writes a loadable `.bloom` save and two map images (combined and
  elevation, Eckert IV), so each checkpoint can be opened in the app's Lake Inspector.
- `branch` loads one checkpoint and runs each requested arm from the same bytes, one after
  the other, with the same step schedule: `off` (the baseline), `outlet-bound[:all]`,
  `uplift-bound[:all]` and `march:<substeps>` (see `app/hydroclimate/river_substeps.py`).
  Every arm starts from an identical state: a default-profile save keeps every cache the
  next step reads, and the `off` arm's first-step `state_hash` is written next to the
  reference's hash at the same age so the runner can confirm the branch is a continuation,
  not a perturbation.

Per step it records (one JSON line on stdout, the whole list in `<out-dir>/<label>-steps.json`):
inland water (visible lake area and liquid volume, lake ice, sea-tier area and volume, the
largest bodies, endorheic roots and their area), outlets (the land area whose water reaches
the ocean, and how much of last step's exporting area still exports: old-outlet retention),
channel and notch geometry, the erosion budget and its closure, the lake, land-water and
glacier ledgers, glacier ice, sea level, step runtime and peak RSS, plus the arm's own log
entry (threatened and relieved pits, or the march's slices, cost and rerouted nodes).

Examples (from the repository root; a density-4 reference run peaks near 4 GB, so follow the
compute-budget skill first):

    backend/.venv/bin/python bin/debug/river_substep_arms.py reference \\
        --seed 544349906 --node-density 4 --myr 25 --checkpoint-myr 5 \\
        --out-dir ~/mantle-bloom-captures/issue867/reference

    backend/.venv/bin/python bin/debug/river_substep_arms.py branch \\
        ~/mantle-bloom-captures/issue867/reference/reference-25.0Myr.bloom \\
        --arm off --arm outlet-bound --arm uplift-bound --arm march:10 \\
        --myr 5 --checkpoint-myr 1 --out-dir ~/mantle-bloom-captures/issue867/branch-25

Saves and images stay in `--out-dir` (outside the repository) with a `SHA256SUMS` file; only
the compact JSON belongs in `analysis/issue867/`. Nothing here changes the model: the arms
are module knobs that default to off.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, os.environ.get("MANTLE_BLOOM_BACKEND", str(Path(__file__).resolve().parents[2] / "backend")))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lake_body_tracking import components  # noqa: E402

from app import continental_ledger, persistence, plates, render_image, surface_parity  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes, river_substeps  # noqa: E402

M3_PER_KM3 = 1e9
LARGE_BODY_KM2 = 10_000.0
IMAGE_SIZE = (1600, 800)
IMAGE_VIEWS = ("combined", "elevation")


def _peak_rss_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


def _commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parent, check=True, capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _controls(world) -> dict:
    """The world's scalar fields (flags, densities, multipliers), as replay logs record them."""
    return {k: v for k, v in vars(world).items() if isinstance(v, (bool, int, float, str)) and not isinstance(v, np.generic)}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _lake_metrics(fields) -> dict:
    area = fields.hypsometry.area_m2 if fields.hypsometry is not None else lakes.CellHypsometry.whole_cells(fields.elevation).area_m2
    ice = fields.lake_ice_m if len(fields.lake_ice_m) == len(area) else np.zeros_like(area)
    water = fields.lake_depth + ice
    visible = water > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    is_sea = fields.is_sea if len(fields.is_sea) == len(area) else np.zeros(len(area), dtype=bool)
    bodies = sorted((b["area_km2"] for b in components(fields)), reverse=True)
    endorheic = [root for root in fields.lake_forest if root.max_depth is None]
    endorheic_area = sum(float(area[root.members].sum()) for root in endorheic) / 1e6
    return {
        "lake_area_km2": float(area[visible].sum() / 1e6),
        "lake_volume_km3": float(np.dot(fields.lake_depth, area) / M3_PER_KM3),
        "lake_ice_km3": float(np.dot(ice, area) / M3_PER_KM3),
        "sea_tier_area_km2": float(area[visible & is_sea].sum() / 1e6),
        "sea_tier_volume_km3": float(np.dot(np.where(is_sea, fields.lake_depth, 0.0), area) / M3_PER_KM3),
        "lake_bodies": len(bodies),
        "large_bodies": int(sum(1 for a in bodies if a >= LARGE_BODY_KM2)),
        "largest_bodies_km2": [round(a, 1) for a in bodies[:5]],
        "endorheic_roots": len(endorheic),
        "endorheic_area_km2": endorheic_area,
    }


def _outlet_metrics(fields, previous: tuple[np.ndarray, np.ndarray] | None) -> tuple[dict, tuple[np.ndarray, np.ndarray]]:
    """Exporting land area now, and the share of last step's exporting area that still exports
    (matched by nearest node, as the model resamples its own previous masks). Returns the
    metrics and this step's `(points, exported)` for the next call."""
    area = fields.hypsometry.area_m2
    exported = river_substeps.exported_to_ocean(fields.flow_target, fields.is_ocean)
    land = ~fields.is_ocean
    metrics = {
        "exporting_land_fraction": float(area[land & exported].sum() / max(area[land].sum(), 1.0)),
        "outlet_retention": None,
        "newly_stranded_area_km2": None,
    }
    if previous is not None:
        prev_points, prev_exported = previous
        from scipy.spatial import cKDTree

        _, idx = cKDTree(prev_points).query(fields.points)
        was_exporting = prev_exported[idx] & land
        still = was_exporting & exported
        metrics["outlet_retention"] = float(area[still].sum() / max(area[was_exporting].sum(), 1.0))
        metrics["newly_stranded_area_km2"] = float(area[was_exporting & ~exported].sum() / 1e6)
    return metrics, (np.asarray(fields.points), exported)


def _ledger_metrics(world) -> dict:
    fields = world.hydrology_cache
    land = fields.land_water_ledger
    glacier = fields.glacier_water_ledger
    lake = fields.lake_water_ledger
    budget = world.erosion_cache.budget if world.erosion_cache is not None else {}
    removed = budget.get("removed_m3", 0.0)
    deposited = budget.get("deposited_m3", 0.0)
    terminal = budget.get("hc_cap_overflow_terminal_m3", 0.0)
    inventory = continental_ledger.inventories(world)
    return {
        "land_water_residual_frac": (
            hydrology.full_land_water_residual_m3(land, glacier) / max(land.rain_m3 + land.meltwater_m3, 1.0)
            if land is not None and glacier is not None else None
        ),
        "lake_ledger_residual_frac": lake.residual_m3 / max(lake.external_supply_m3, 1.0) if lake is not None else None,
        "erosion_removed_km3": removed / M3_PER_KM3,
        "erosion_deposited_km3": deposited / M3_PER_KM3,
        "erosion_terminal_km3": terminal / M3_PER_KM3,
        "erosion_lake_silt_km3": budget.get("lake_silt_m3", 0.0) / M3_PER_KM3,
        "erosion_glacial_km3": budget.get("glacial_erosion_m3", 0.0) / M3_PER_KM3,
        # removed == deposited + terminal when the ocean-deposition knob is 1 (erosion.py).
        "rock_residual_frac": (removed - deposited - terminal) / max(removed, 1.0),
        "continental_tracked_km3": (
            inventory["surface_continental_derived_m3"] + sum(v for k, v in inventory.items() if k.endswith("_m3") and k != "surface_continental_derived_m3")
        ) / M3_PER_KM3,
    }


def _channel_metrics(world) -> dict:
    order = world.hydrology_cache.plates_in_order
    channel = plates.collect_all_channel_depth(order)
    notch = plates.collect_all_breach_notch_depth(order)
    glacier = plates.collect_all_glacier_depth(order)
    land = ~world.hydrology_cache.is_ocean
    area = world.hydrology_cache.hypsometry.area_m2
    return {
        "channel_depth_land_mean_m": float(channel[land].mean()) if land.any() else 0.0,
        "channel_depth_p99_m": float(np.percentile(channel[land], 99)) if land.any() else 0.0,
        "channel_depth_max_m": float(channel.max()) if len(channel) else 0.0,
        "channel_volume_km3": float(np.dot(channel, area) / M3_PER_KM3),
        "notch_sum_m": float(notch.sum()),
        "notched_nodes": int(np.count_nonzero(notch > 0.0)),
        "glacier_ice_km3": float(np.dot(glacier, area) / M3_PER_KM3),
    }


class Run:
    """One labelled trajectory: steps a world, records a row per step, writes checkpoints."""

    def __init__(self, label: str, world, args, out_dir: Path, summary: dict) -> None:
        self.label = label
        self.world = world
        self.args = args
        self.out_dir = out_dir
        self.summary = summary
        self.rows: list[dict] = []
        self.previous_outlets = None
        self.steps_path = out_dir / f"{label}-steps.json"

    def measure(self, seconds: float) -> dict:
        world = self.world
        fields = world.hydrology_cache
        outlets, self.previous_outlets = _outlet_metrics(fields, self.previous_outlets)
        arm_log = dict(river_substeps.LOG[-1]) if river_substeps.LOG else None
        row = {
            "label": self.label,
            "step": world.steps_taken,
            "elapsed_years": world.elapsed_years,
            "seconds": round(seconds, 2),
            "peak_rss_gb": round(_peak_rss_gb(), 3),
            "sea_level_m": world.sea_level_m,
            **_lake_metrics(fields),
            **outlets,
            **_channel_metrics(world),
            **_ledger_metrics(world),
            "arm_log": arm_log,
            # Branch identity: a branch's `off` arm must reproduce the reference's hash at
            # every age it shares with it.
            "state_hash": surface_parity.state_hash(world),
        }
        return row

    def checkpoint(self) -> dict:
        world = self.world
        age_myr = world.elapsed_years / 1e6
        stem = f"{self.label}-{age_myr:.1f}Myr"
        entry = {"label": self.label, "elapsed_years": world.elapsed_years, "state_hash": surface_parity.state_hash(world)}
        save_path = self.out_dir / f"{stem}.bloom"
        data = persistence.save_exact_bytes(world) if self.args.exact else persistence.save_world_bytes(world, name=stem)
        save_path.write_bytes(data)
        entry["save"] = save_path.name
        entry["save_sha256"] = hashlib.sha256(data).hexdigest()
        for view in IMAGE_VIEWS:
            image_path = self.out_dir / f"{stem}-{view}.png"
            image_path.write_bytes(render_image.render_png(world, "eckert4", view, *IMAGE_SIZE))
            entry[f"{view}_png"] = image_path.name
        # The first map render is the process's memory peak on a small world (about 2.3 GB
        # transient at 1600x800), so every later step's `peak_rss_gb` reads it; this entry
        # lets a reader attribute the peak.
        entry["peak_rss_gb_after_render"] = round(_peak_rss_gb(), 3)
        with (self.out_dir / "SHA256SUMS").open("a") as sums:
            for name in [save_path.name] + [entry[f"{view}_png"] for view in IMAGE_VIEWS]:
                sums.write(f"{_sha256(self.out_dir / name)}  {name}\n")
        return entry

    def run(self) -> None:
        args = self.args
        steps = max(1, int(round(args.myr * 1e6 / args.step_years)))
        per_checkpoint = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
        started = time.perf_counter()
        for step in range(1, steps + 1):
            t0 = time.perf_counter()
            world_mod.step_world(self.world, args.step_years)
            row = self.measure(time.perf_counter() - t0)
            self.rows.append(row)
            print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in row.items() if k != "arm_log"}), flush=True)
            if step % per_checkpoint == 0 or step == steps:
                self.summary["checkpoints"].append(self.checkpoint())
            self.steps_path.write_text(json.dumps(self.rows, indent=1))
            self.summary["runs"][self.label] = {"steps": step, "seconds": round(time.perf_counter() - started, 1), "peak_rss_gb": round(_peak_rss_gb(), 3)}
            (self.out_dir / self.summary["file"]).write_text(json.dumps(self.summary, indent=1))


def _common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--myr", type=float, required=True, help="how far to step, Myr")
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--checkpoint-myr", type=float, default=5.0, help="save and render every this many Myr")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--exact", action="store_true", help="write exact-profile saves (keeps the erosion cache; much larger)")
    parser.add_argument(
        "--glacier-march", choices=("on", "off", "keep"), default="keep",
        help="World.simulate_glacier_march for the run (default: keep the world's own flag; a generated world has it on)",
    )


def _apply_flags(world, args) -> None:
    if args.glacier_march != "keep":
        world.simulate_glacier_march = args.glacier_march == "on"


def _summary(args, kind: str, file: str) -> dict:
    return {
        "issue": 867,
        "kind": kind,
        "file": file,
        "commit": _commit(),
        "step_years": args.step_years,
        "myr": args.myr,
        "checkpoint_myr": args.checkpoint_myr,
        "save_profile": "exact" if args.exact else "default",
        "checkpoints": [],
        "runs": {},
    }


def cmd_reference(args) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = _summary(args, "reference", "summary.json")
    river_substeps.set_arm(river_substeps.RiverSubstepArm())
    t0 = time.perf_counter()
    if args.load is not None:
        data = args.load.read_bytes()
        world = persistence.load_world_bytes(data)
        summary["source"] = {"save": str(args.load), "sha256": hashlib.sha256(data).hexdigest()}
    else:
        world = world_mod.generate_world(seed=args.seed, node_density=args.node_density)
        summary["source"] = {"seed": args.seed, "node_density": args.node_density}
    _apply_flags(world, args)
    summary["generate_seconds"] = round(time.perf_counter() - t0, 1)
    summary["controls"] = _controls(world)
    summary["nodes"] = int(sum(p.node_count() for p in world.plates))
    run = Run("reference", world, args, args.out_dir, summary)
    if args.load is None:
        # The opening state too, so the branch runner can start at age zero if it wants.
        summary["checkpoints"].append(run.checkpoint())
    run.run()


def _branch_one(args, text: str) -> None:
    """One arm from the checkpoint, in this process; writes `<label>-summary.json`."""
    arm = river_substeps.parse_arm(text)
    label = arm.label.replace(":", "-")
    summary = _summary(args, "branch", f"{label}-summary.json")
    data = args.checkpoint.read_bytes()
    summary["source"] = {"save": str(args.checkpoint), "sha256": hashlib.sha256(data).hexdigest()}
    summary["arm"] = arm.label
    river_substeps.set_arm(arm)
    try:
        world = persistence.load_world_bytes(data)
        _apply_flags(world, args)
        summary["controls"] = _controls(world)
        Run(label, world, args, args.out_dir, summary).run()
    finally:
        river_substeps.set_arm(river_substeps.RiverSubstepArm())


def cmd_branch(args) -> None:
    """Each arm runs in its own process, one after the other: hydrology's per-process caches
    (`_LAST_OCEAN_TREE`, the nearest-node cache) keep the previous arm's world alive
    otherwise, and a process per arm also gives each arm its own peak RSS. `summary.json`
    collects the per-arm summaries."""
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if len(args.arm) == 1:
        _branch_one(args, args.arm[0])
        return
    for text in args.arm:
        command = [sys.executable, str(Path(__file__).resolve()), "branch", str(args.checkpoint), "--arm", text]
        command += ["--myr", str(args.myr), "--step-years", str(args.step_years), "--checkpoint-myr", str(args.checkpoint_myr)]
        command += ["--out-dir", str(args.out_dir), "--glacier-march", args.glacier_march] + (["--exact"] if args.exact else [])
        print(json.dumps({"spawn": text}), flush=True)
        subprocess.run(command, check=True)
    merged = {"issue": 867, "kind": "branch", "commit": _commit(), "arms": {}}
    for text in args.arm:
        label = river_substeps.parse_arm(text).label.replace(":", "-")
        merged["arms"][label] = json.loads((args.out_dir / f"{label}-summary.json").read_text())
    (args.out_dir / "summary.json").write_text(json.dumps(merged, indent=1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    ref = sub.add_parser("reference", help="step the production model and write checkpoints")
    ref.add_argument("--seed", type=int, default=544349906, help="the #866 seed by default")
    ref.add_argument("--node-density", type=float, default=4.0)
    ref.add_argument("--load", type=Path, default=None, help="continue a save instead of generating")
    _common_arguments(ref)
    branch = sub.add_parser("branch", help="run arms from one checkpoint")
    branch.add_argument("checkpoint", type=Path)
    branch.add_argument("--arm", action="append", required=True, help="off, outlet-bound[:all], uplift-bound[:all] or march:<substeps>; repeat")
    _common_arguments(branch)
    args = parser.parse_args()
    if args.command == "reference":
        cmd_reference(args)
    else:
        cmd_branch(args)


if __name__ == "__main__":
    main()

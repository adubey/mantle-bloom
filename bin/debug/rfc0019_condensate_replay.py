"""Replay one registered rank-16 coupled arm for #978.

Usage: backend/.venv/bin/python bin/debug/rfc0019_condensate_replay.py CHECKPOINT \
    --anchors ANCHORS --g4-runner PINNED_RUNNER.py --arm immediate|six-hour --out PRIVATE.json

Run both arms sequentially from the same inputs. See analysis/issue978/README.md.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import numpy as np
from scipy.spatial import cKDTree

from replay_lib import upstream_catchment
from rfc0019_condensate_fallout import FALLOUT_SECONDS, install_fallout, load_g4
from rfc0019_g3_compare import _sha256

ROOT = Path(__file__).resolve().parents[2]
G4_COMMIT = "c703237694a353643d7d3cad4c0f9637e157f65b"
CHECKPOINT_HASH = "5e93f0feacf0b54d4016321bade8b0d36f435760ea55c1583a1678c94f2e9b41"
ANCHOR_HASH = "b679d471f2ad0dab44d4e4208b02c7de9c9a878477d74224323ef0d2eeeeff0d"
STEP_YEARS = 10_000
STEPS = 40


def validate_inputs(args, parser, runner):
    for path, expected in ((args.save, CHECKPOINT_HASH), (args.anchors, ANCHOR_HASH)):
        if _sha256(path) != expected:
            parser.error(f"{path.name} does not match the registered rank-16 SHA-256")

    pinned_root = args.g4_runner.resolve().parents[2]
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=pinned_root, text=True).strip()
    tracked_changes = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=no", "--", "backend", "bin/debug"],
        cwd=pinned_root, text=True,
    )
    if head != G4_COMMIT or tracked_changes:
        parser.error("use a clean backend and debug tree at the #977 G4 commit")
    # Both arms import the historical backend and use the current checkout's debug helpers.
    for name in runner.FINGERPRINT_FILES:
        if name.startswith("bin/debug/") and name != "bin/debug/rfc0019_g4_coupling.py":
            if _sha256(ROOT / name) != _sha256(pinned_root / name):
                parser.error(f"shared helper {name} differs from the pinned G4 checkout")
    return pinned_root / "backend"


class BasinReadout:
    """Sample a fixed opening catchment on each step's current accounting mesh."""

    def __init__(self, hydro, anchor_points, visible_depth, replay):
        if hydro is None:
            raise ValueError("the registered checkpoint lacks its opening hydrology cache")
        self.tree = cKDTree(hydro.points, copy_data=True)
        anchor = int(self.tree.query(anchor_points, k=1)[1][0])
        wet = np.asarray(hydro.lake_depth) > visible_depth
        if not wet[anchor]:
            raise ValueError("the registered anchor must cover a lake in the opening checkpoint")
        from app.hydroclimate.hydrology import lake_components

        labels = replay.lake_bodies(hydro, np.ones(len(wet)), wet, lake_components)["labels"]
        members = np.flatnonzero(labels == labels[anchor])
        self.opening_mask = upstream_catchment(hydro.flow_target, hydro.is_ocean, members)
        self.visible_depth = visible_depth
        self.runoff = None
        self.overflow_seen = False

    def record_routing_pass(self, runoff, years, *, source_pass):
        if years != STEP_YEARS:
            raise RuntimeError(f"expected a {STEP_YEARS}-year routing pass, got {years}")
        if source_pass:
            if self.runoff is not None:
                raise RuntimeError("duplicate source routing pass before the catchment readout")
            self.runoff = np.array(runoff, copy=True)
            return
        if self.runoff is None or self.overflow_seen or np.any(runoff):
            raise RuntimeError("expected one zero-runoff overflow pass after the source pass")
        self.overflow_seen = True

    def read(self, hydro, areas, precipitation):
        if not self.overflow_seen:
            raise RuntimeError("expected both routing passes before the catchment readout")
        if self.runoff is None or self.runoff.shape != areas.shape:
            raise RuntimeError("no aligned runoff input was captured for this step")
        inside = self.opening_mask[self.tree.query(hydro.points, k=1)[1]]
        if not np.any(inside):
            raise RuntimeError("opening catchment has no overlap with the current accounting mesh")
        wet = inside & (hydro.lake_depth > self.visible_depth)
        total = {
            "nodes": int(inside.sum()),
            "area_km2": float(areas[inside].sum() / 1e6),
            "liquid_precipitation_m3": float(precipitation["liquid_m3"][inside].sum()),
            "frozen_precipitation_m3": float(precipitation["frozen_m3"][inside].sum()),
            "generated_dry_runoff_m3": float(self.runoff[inside].sum() * STEP_YEARS),
            "settled_external_lake_supply_m3": float(
                hydro.lake_external_supply_m3_per_yr[inside].sum() * STEP_YEARS
            ),
            "river_evaporation_m3": float(hydro.river_water_evaporation_m3[inside].sum()),
            "visible_lake_area_km2": float(areas[wet].sum() / 1e6),
            "visible_lake_volume_km3": float(np.dot(hydro.lake_depth[wet], areas[wet]) / 1e9),
        }
        self.runoff = None
        self.overflow_seen = False
        return total


def install_replay_readout(runner, arm):
    from app.hydroclimate import hydrology

    replay = runner.rfc0016
    load = replay.load_replay
    target_lake = replay.target_lake
    route = hydrology.route_river_reservoirs
    basin = None

    def load_with_basin(args, issue):
        nonlocal basin
        world, result = load(args, issue)
        anchors = np.asarray(json.loads(args.anchors.read_text())["anchor_points"])
        basin = BasinReadout(world.hydrology_cache, anchors, hydrology.LAKE_MIN_VISIBLE_DEPTH_M, replay)
        result["condensate_replay"] = {
            "issue": 978, "arm": arm, "fallout_seconds": FALLOUT_SECONDS if arm == "six-hour" else 0,
            "opening_catchment_nodes": int(basin.opening_mask.sum()),
            "script_sha256": _sha256(Path(__file__)),
            "fallout_script_sha256": _sha256(ROOT / "bin/debug/rfc0019_condensate_fallout.py"),
        }
        return world, result

    def measure(hydro, areas, anchors, depth, components, precipitation, previous=None):
        target, bodies = target_lake(hydro, areas, anchors, depth, components, precipitation, previous)
        target["fixed_catchment"] = basin.read(hydro, areas, precipitation)
        return target, bodies

    def capture(
        points, elevation, is_ocean, flow_target, stores_water, runoff_m3_per_yr,
        meltwater_m3_per_yr, prior_storage_m3, discharge_m3_s, years, runoff_loss, meltwater_loss, **kwargs,
    ):
        # The pinned backend passes this keyword on the source call even when its value is None.
        basin.record_routing_pass(runoff_m3_per_yr, years, source_pass="initial_discharge_m3_s" in kwargs)
        return route(
            points, elevation, is_ocean, flow_target, stores_water, runoff_m3_per_yr,
            meltwater_m3_per_yr, prior_storage_m3, discharge_m3_s, years, runoff_loss, meltwater_loss, **kwargs,
        )

    return (
        patch.object(replay, "load_replay", load_with_basin),
        patch.object(replay, "target_lake", measure),
        patch.object(hydrology, "route_river_reservoirs", capture),
    )


def validate_output(path, arm):
    steps = json.loads(path.read_text())["steps"]
    if len(steps) != STEPS:
        raise RuntimeError(f"expected {STEPS} completed steps, got {len(steps)}")
    for step in steps:
        limits = ("normalized_land_residual", "normalized_lake_residual", "normalized_river_residual")
        if not all(0 <= step[key] <= 1e-3 for key in limits) or step["normalized_stranded_overflow"] != 0:
            raise RuntimeError(f"surface ledger failed at step {step['step']}")
        if arm == "six-hour" and not 0 <= step["atmosphere_budget"]["condensate_relative_residual"] <= 1e-8:
            raise RuntimeError(f"condensate ledger failed at step {step['step']}")
    print(f"{arm}: {sum(step['target']['anchor_wet'] for step in steps)}/{STEPS} wet anchor steps")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("save", type=Path)
    parser.add_argument("--anchors", required=True, type=Path)
    parser.add_argument("--g4-runner", required=True, type=Path)
    parser.add_argument("--arm", required=True, choices=("immediate", "six-hour"))
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("--out already exists")
    runner = load_g4(args.g4_runner)
    backend = validate_inputs(args, parser, runner)
    sys.path.insert(0, str(backend))
    if args.arm == "six-hour":
        install_fallout(runner)

    runner._install_replay_wrapper()
    hooks = install_replay_readout(runner, args.arm)
    argv = [
        sys.argv[0], str(args.save), "--anchors", str(args.anchors), "--out", str(args.out),
        "--backend", str(backend), "--g4-atmosphere", "conservative", "--treatment",
        "--myr", str(STEPS * STEP_YEARS / 1e6), "--step-years", str(STEP_YEARS), "--checkpoint-myr", "0.4",
    ]
    with patch.object(sys, "argv", argv), hooks[0], hooks[1], hooks[2]:
        runner.rfc0016.main()
    validate_output(args.out, args.arm)


if __name__ == "__main__":
    main()

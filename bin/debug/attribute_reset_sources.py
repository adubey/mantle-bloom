#!/usr/bin/env python3
"""Caller-attributed juvenile accounting on a stripped current-main baseline (issue #407).

Diagnosis only for RFC 0007's source-attribution gate: it measures how much of the
booked continental juvenile flows through the existing elevation-driven melt-reset
route (`lithosphere_plate._erupt_melted_nodes` recoding by prior elevation, then
juvenile booked on the recoded host) versus every other booking caller. It shares
the stripped baseline with #406: default `world.generate_world(seed)`, cratons on,
no control overrides, 5 Myr steps, 100 Myr checkpoints. It does not implement plume
physics and changes no model code: every tap wraps a function, records, delegates,
and is restored afterwards.

Usage:
    backend/.venv/bin/python bin/debug/attribute_reset_sources.py \\
        --out analysis/issue407/attribution.json --save-dir /tmp/mantle-bloom-issue407
    backend/.venv/bin/python bin/debug/attribute_reset_sources.py --probe-only

One world runs at a time (replays are RAM-heavy). Saves stay outside Git; only the
compact JSON and analysis/issue407/README.md are reviewed.
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import platform
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_cratons  # noqa: E402
import ledger_history  # noqa: E402
from app import continental_ledger  # noqa: E402
from app import elevation_lines  # noqa: E402
from app import lithosphere_plate  # noqa: E402
from app import merge_split  # noqa: E402
from app import persistence  # noqa: E402
from app import quad_tectonics  # noqa: E402
from app import rheology  # noqa: E402
from app.elevation_lines import effective_is_continental_from_codes  # noqa: E402
from app import lithosphere  # noqa: E402

SCRIPT_NAME = Path(__file__).name
APP_ROOT = str(Path(lithosphere_plate.__file__).resolve().parent) + "/"


def _caller_tag(skip_modules: tuple[str, ...]) -> str:
    """First backend/app caller outside the skipped modules, as module.func:line."""
    for frame in inspect.stack()[2:]:
        filename = frame.filename
        if not filename.startswith(APP_ROOT):
            continue
        module = Path(filename).stem
        if module in skip_modules:
            continue
        return f"{module}.{frame.function}:{frame.lineno}"
    return "unknown"


class AttributionTaps:
    """Read-only wrappers around ledger bookings and melt resets."""

    def __init__(self) -> None:
        self.caller_volumes_m3: dict[tuple[str, str], float] = defaultdict(float)
        self.recode_calls: int = 0
        self.melt_nodes: int = 0
        self.recoded_exposed: int = 0
        self.reset_thickness_m: float = 0.0
        self.reset_volume_exact_m3: float = 0.0
        self.reset_volume_estimated_m3: float = 0.0
        self.reset_volume_estimate_basis: str = ""
        self.prior_elev_m: list[float] = []
        self.prior_code_counts: dict[str, int] = defaultdict(int)
        self._original_record = continental_ledger.record
        self._original_erupt = lithosphere_plate._erupt_melted_nodes

    def install(self) -> None:
        taps = self

        def record(world, account, volume_m3):
            taps.caller_volumes_m3[(_caller_tag(("continental_ledger",)), account)] += float(volume_m3)
            return taps._original_record(world, account, volume_m3)

        def erupt(world, plate_id, rng_index, hc, hm, codes, is_volcano, remaining, melting, prior_elevation):
            melting = np.asarray(melting, dtype=bool)
            tag = _caller_tag((SCRIPT_NAME.replace(".py", ""),))
            pre_idx = np.flatnonzero(melting) if np.any(melting) else np.array([], dtype=int)
            pre = None
            if len(pre_idx):
                plate = next((p for p in world.plates if p.plate_id == plate_id), None)
                plate_cont = bool(plate.crust_type == "continental") if plate is not None else False
                pre = {
                    "tag": tag,
                    "codes": np.asarray(codes)[melting].copy(),
                    "elev": np.asarray(prior_elevation, dtype=float)[melting].copy(),
                    "hc": np.asarray(hc, dtype=float)[melting].copy(),
                    "plate_cont": plate_cont,
                    "full_plate": plate is not None and len(np.asarray(hc)) == plate.node_count(),
                    "areas": np.asarray(plate.node_areas_m2()) if plate is not None else None,
                }
            result = taps._original_erupt(
                world, plate_id, rng_index, hc, hm, codes, is_volcano, remaining,
                melting, prior_elevation,
            )
            if pre is not None:
                taps._digest_call(
                    pre,
                    np.asarray(codes)[melting],
                    np.asarray(hc, dtype=float)[melting],
                    np.flatnonzero(melting),
                )
            return result

        continental_ledger.record = record
        lithosphere_plate._erupt_melted_nodes = erupt
        quad_tectonics._erupt_melted_nodes = erupt

    def _digest_call(self, pre, post_codes, post_hc, melt_positions) -> None:
        self.recode_calls += 1
        pre_eff = effective_is_continental_from_codes(pre["codes"], pre["plate_cont"])
        post_eff = effective_is_continental_from_codes(np.asarray(post_codes), pre["plate_cont"])
        recoded = ~pre_eff & post_eff
        self.melt_nodes += int(len(pre_eff))
        exposed = recoded & (pre["elev"] > 0.0)
        self.recoded_exposed += int(np.sum(exposed))
        thick = np.maximum(post_hc - pre["hc"], 0.0)[exposed]
        self.reset_thickness_m += float(np.sum(thick))
        for code_value, name in (
            (elevation_lines.CRUST_TYPE_INHERIT, "inherit"),
            (elevation_lines.CRUST_TYPE_OCEANIC, "oceanic"),
            (elevation_lines.CRUST_TYPE_CONTINENTAL, "continental"),
        ):
            self.prior_code_counts[name] += int(np.sum(pre["codes"][exposed] == code_value))
        self.prior_elev_m.extend(float(v) for v in pre["elev"][exposed][:50])
        if not np.any(exposed):
            return
        areas = pre["areas"]
        if pre["full_plate"] and areas is not None:
            # Full-plate call: melt-space positions align with the plate node order.
            self.reset_volume_exact_m3 += float(np.sum(thick * areas[melt_positions[exposed]]))
        else:
            mean_area = float(np.mean(areas)) if areas is not None else 0.0
            self.reset_volume_estimated_m3 += float(np.sum(thick)) * mean_area
            if not self.reset_volume_estimate_basis and areas is not None:
                spread = float(np.max(areas) - np.min(areas))
                self.reset_volume_estimate_basis = (
                    f"mean node area {mean_area:.3e} m2, max-min spread {spread:.3e} m2"
                )

    def uninstall(self) -> None:
        continental_ledger.record = self._original_record
        lithosphere_plate._erupt_melted_nodes = self._original_erupt
        quad_tectonics._erupt_melted_nodes = self._original_erupt


def run_synthetic_probe() -> dict:
    """Drive the real `_erupt_melted_nodes` with a synthetic exposed oceanic column."""
    critical = float(rheology.RIFT_CRITICAL_THICKNESS_M)
    world = SimpleNamespace(seed=7, elapsed_years=0.0)
    codes = np.array(
        [elevation_lines.CRUST_TYPE_OCEANIC, elevation_lines.CRUST_TYPE_OCEANIC,
         elevation_lines.CRUST_TYPE_CONTINENTAL]
    )
    hc = np.array([critical / 2, critical / 2, critical / 2])
    hm = np.array([60_000.0, 60_000.0, 60_000.0])
    is_volcano = np.zeros(3, dtype=bool)
    remaining = np.zeros(3)
    melting = np.array([True, True, True])
    prior_elev = np.array([400.0, -3000.0, 500.0])
    before = {"codes": codes.tolist(), "hc_m": hc.tolist(), "prior_elev_m": prior_elev.tolist()}
    lithosphere_plate._erupt_melted_nodes(
        world, 0, 0, hc, hm, codes, is_volcano, remaining, melting, prior_elev
    )
    return {
        "nodes": "[exposed oceanic, drowned oceanic, exposed continental]",
        "before": before,
        "after": {"codes": codes.tolist(), "hc_m": hc.tolist(), "is_volcano": is_volcano.tolist()},
        "constants_m": {
            "rift_critical_thickness": critical,
            "reference_hc_continental": float(lithosphere.REFERENCE_HC_CONTINENTAL_M),
            "reference_hc_oceanic": float(lithosphere.REFERENCE_HC_OCEANIC_M),
        },
        "caller_would_book_juvenile_m_per_node": [float(max(hc[0] - critical / 2, 0.0)), 0.0,
                                                  float(max(hc[2] - critical / 2, 0.0))],
    }


def run_self_test() -> dict:
    """Validate the tap arithmetic on synthetic full-plate and subset calls."""
    taps = AttributionTaps()
    taps.install()
    try:
        from unittest import mock

        results = {}
        for tag, n_nodes, full in (("lithosphere_plate.deform_columns:999", 4, True),
                                   ("quad_tectonics._open_rift:999", 4, False)):
            plate = SimpleNamespace(
                plate_id=3, crust_type="continental",
                node_count=lambda: 4,
                node_areas_m2=lambda: np.full(4, 1e12),
            )
            world = SimpleNamespace(seed=7, elapsed_years=0.0, plates=[plate])
            codes = np.array([elevation_lines.CRUST_TYPE_OCEANIC, elevation_lines.CRUST_TYPE_OCEANIC,
                              elevation_lines.CRUST_TYPE_CONTINENTAL, elevation_lines.CRUST_TYPE_OCEANIC])
            hc = np.full(4, 4_000.0)
            hm = np.full(4, 60_000.0)
            melting = np.array([True, True, True, False])
            prior_elev = np.array([400.0, -3000.0, 500.0, 100.0])
            with mock.patch.object(
                taps, "_original_erupt",
                lambda w, pid, r, h, m, c, v, rem, melt, elev: (
                    c.__setitem__((melt & (np.asarray(elev) > 0.0)),
                                      elevation_lines.CRUST_TYPE_CONTINENTAL),
                    h.__setitem__((melt & (np.asarray(elev) > 0.0)), 35_000.0),
                    h.__setitem__((melt & ~(np.asarray(elev) > 0.0)), 11_000.0),
                ),
            ):
                if not full:
                    sub = np.array([0, 1, 2])
                    codes_s, hc_s, melt_s, elev_s = codes[sub], hc[sub].copy(), melting[sub], prior_elev[sub]
                    hm_s = hm[sub].copy()
                    lithosphere_plate._erupt_melted_nodes(
                        world, 3, 0, hc_s, hm_s, codes_s, np.zeros(3, bool), np.zeros(3),
                        melt_s, elev_s,
                    )
                else:
                    lithosphere_plate._erupt_melted_nodes(
                        world, 3, 0, hc, hm, codes, np.zeros(4, bool), np.zeros(4),
                        melting, prior_elev,
                    )
            results[tag] = {
                "recoded_exposed": taps.recoded_exposed,
                "exact_m3": taps.reset_volume_exact_m3,
                "estimated_m3": taps.reset_volume_estimated_m3,
            }
            taps.recoded_exposed = 0
            taps.reset_volume_exact_m3 = 0.0
            taps.reset_volume_estimated_m3 = 0.0
        # One exposed oceanic node recoded with 31 km of reset thickness over 1e12 m2.
        assert results["lithosphere_plate.deform_columns:999"]["recoded_exposed"] == 1
        assert results["lithosphere_plate.deform_columns:999"]["exact_m3"] == 31_000.0 * 1e12
        assert results["quad_tectonics._open_rift:999"]["recoded_exposed"] == 1
        assert results["quad_tectonics._open_rift:999"]["estimated_m3"] == 31_000.0 * 1e12
        return {"self_test": "pass", "cases": results}
    finally:
        taps.uninstall()


def _measure(seed: int, years: float, every: float, step: float, save_dir: Path) -> dict:
    taps = AttributionTaps()
    original_snapshot = compare_cratons._snapshot
    captured: dict = {}

    def snapshot(world, events):
        captured["world"] = world
        base = original_snapshot(world, events)
        base["continental_material_ledger"] = continental_ledger.inventories(world)
        base["juvenile_by_caller_m3"] = {
            f"{tag} :: {account}": volume for (tag, account), volume in taps.caller_volumes_m3.items()
        }
        base["melt_recodes"] = {
            "erupt_calls": taps.recode_calls,
            "melt_nodes": taps.melt_nodes,
            "recoded_exposed_oceanic": taps.recoded_exposed,
            "reset_thickness_m": taps.reset_thickness_m,
            "reset_volume_exact_m3": taps.reset_volume_exact_m3,
            "reset_volume_estimated_m3": taps.reset_volume_estimated_m3,
        }
        return base

    compare_cratons._snapshot = snapshot
    taps.install()
    try:
        t0 = time.time()
        rows = compare_cratons.run_job(seed, "on", years, every, step, {})
        world = captured["world"]
        save_dir.mkdir(parents=True, exist_ok=True)
        save = save_dir / f"seed{seed}-stripped-407.mbworld"
        data = persistence.save_world_bytes(world)
        save.write_bytes(data)
        return {
            "seed": seed,
            "wall_s": round(time.time() - t0, 1),
            "checkpoints": rows,
            "splits_per_checkpoint": [dict(r.get("events", {})) for r in rows],
            "caller_volumes_m3": {
                f"{tag} :: {account}": volume for (tag, account), volume in taps.caller_volumes_m3.items()
            },
            "recode_summary": {
                "erupt_calls": taps.recode_calls,
                "melt_nodes": taps.melt_nodes,
                "recoded_exposed_oceanic": taps.recoded_exposed,
                "reset_thickness_m": taps.reset_thickness_m,
                "reset_volume_exact_m3": taps.reset_volume_exact_m3,
                "reset_volume_estimated_m3": taps.reset_volume_estimated_m3,
                "estimate_basis": taps.reset_volume_estimate_basis,
                "prior_code_of_reset_nodes": dict(taps.prior_code_counts),
                "prior_elev_sample_m": taps.prior_elev_m[:50],
            },
            "ledger_history": ledger_history.summarize(world.stats_history),
            "save_bytes": len(data),
            "save_sha256": hashlib.sha256(data).hexdigest(),
        }
    finally:
        taps.uninstall()
        compare_cratons._snapshot = original_snapshot


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--years", type=float, default=5e8)
    parser.add_argument("--every", type=float, default=1e8)
    parser.add_argument("--step", type=float, default=5e6)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--save-dir", type=Path, required=True)
    parser.add_argument("--probe-only", action="store_true")
    parser.add_argument("--self-test-only", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    if not args.probe_only:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
        model_diff = subprocess.check_output(
            ["git", "diff", "HEAD", "--", "backend/app", "backend/constraints.txt"], cwd=root, text=True
        )
        if model_diff:
            parser.error("baseline needs an unchanged model and constraints")
    else:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    tap_self_test = run_self_test()
    if args.self_test_only:
        print(json.dumps(tap_self_test, indent=1))
        return
    result: dict = {
        "baseline_commit": commit,
        "python": platform.python_version(),
        "constraints_sha256": hashlib.sha256((root / "backend/constraints.txt").read_bytes()).hexdigest(),
        "generation": "world.generate_world(seed), cratons on, all other arguments default; no control overrides",
        "step_years": args.step,
        "every_years": args.every,
        "years": args.years,
        "split_min_age_years": merge_split.SPLIT_MIN_AGE_YEARS,
        "tap_self_test": tap_self_test,
        "probe": run_synthetic_probe(),
        "runs": [],
    }
    if not args.probe_only:
        for seed in map(int, args.seeds.split(",")):
            result["runs"].append(_measure(seed, args.years, args.every, args.step, args.save_dir))
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(json.dumps(result, indent=1) + "\n")
    else:
        print(json.dumps(result["probe"], indent=1))


if __name__ == "__main__":
    main()

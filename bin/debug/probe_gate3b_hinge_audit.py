#!/usr/bin/env python3
"""Audit the slab-ribbon hinge load on the registered Gate 3b fronts, body only (#699).

Advances the two registered ribbons for 20 Myr with no plate coupling, shell solve or
damage, so it exposes no rift outcome. It reports the hinge load against the active
buoyant weight, its alignment with the down-dip direction and its timestep sensitivity,
with and without the conveyor drive.

Usage (one process; loads a young saved world):
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
      backend/.venv/bin/python bin/debug/probe_gate3b_hinge_audit.py \
      --save <registered.mbworld> --out analysis/issue699/hinge-audit-seed447276134-14.json
"""
from __future__ import annotations

import argparse
import copy
from pathlib import Path

import numpy as np
from app import geometry
from app import slab_body as slab_model
import gate3b_diagnostic_common as common
import probe_gate3b_evolving_slabs as evolving

STEP_CASES_MYR = (0.5, 1.0, 2.0)
DURATION_MYR = 20.0


def _weight_n(body):
    active = body.active_station_mask
    return float(np.sum(
        body.density_anomaly_kg_m3[active] * slab_model.GRAVITY_M_S2 * body.station_volume_m3[active]
    ))


def _down_dip(body):
    points = body.centerline_local_xyz_m
    segment = points[1] - points[0]
    radial = body.trench_local_xyz
    return geometry.normalize(segment - np.dot(segment, radial) * radial)


def _history(template, step_myr, speed_m_per_s):
    body = copy.deepcopy(template)
    rows = []
    for step in range(1, int(round(DURATION_MYR / step_myr)) + 1):
        budget = slab_model.advance_slab_body(
            body, step_myr * 1e6,
            mantle_viscosity_pa_s=evolving.MANTLE_VISCOSITY_PA_S,
            shear_zone_thickness_m=evolving.SHEAR_ZONE_M,
            convergence_speed_m_per_s=speed_m_per_s,
        )
        force = np.asarray(budget["attachment_force_local_n"], dtype=float)
        magnitude = float(np.linalg.norm(force))
        weight = _weight_n(body)
        rows.append({
            "elapsed_myr": step * step_myr,
            "hinge_force_n": magnitude,
            "active_buoyant_weight_n": weight,
            "hinge_over_weight": magnitude / max(weight, 1.0),
            "hinge_per_trench_length_n_per_m": magnitude / body.width_m,
            "down_dip_alignment": float(np.dot(force, _down_dip(body)) / max(magnitude, 1.0)),
            "first_segment_length_m": float(np.linalg.norm(
                body.centerline_local_xyz_m[1] - body.centerline_local_xyz_m[0]
            )),
            "body_length_m": body.metrics()["length_m"],
        })
    return rows


def run(save_path: Path, out_path: Path, column_thickness: bool = False) -> dict:
    initial = evolving._load_target(save_path)
    plate = initial["plate"]
    speeds = evolving._closing_speeds(initial["world"], plate, initial["active_masks"])
    bodies, records = evolving._create_bodies(plate, initial, column_thickness=column_thickness)
    post_transfer_speeds = evolving._closing_speeds(initial["world"], plate, initial["active_masks"])
    fronts = []
    for index, body in enumerate(bodies):
        cases = {}
        for label, speed in (("no_conveyor", 0.0), ("measured_conveyor", speeds[index])):
            cases[label] = {f"{step:g}_myr": _history(body, step, speed) for step in STEP_CASES_MYR}
        fronts.append({
            "side": records[index]["side"],
            "width_m": float(body.width_m),
            "thickness_m": float(body.thickness_m),
            "bending_rigidity_nm2": float(body.bending_rigidity_nm2),
            "convergence_speed_m_per_s": float(speeds[index]),
            "post_transfer_convergence_speed_m_per_s": float(post_transfer_speeds[index]),
            "cases": cases,
        })
    report = {
        "schema_version": 1,
        "scope": "Body-only hinge audit; no shell, damage or rift outcome.",
        "seed": int(initial["world"].seed),
        "plate_id": evolving.PLATE_ID,
        "duration_myr": DURATION_MYR,
        "fronts": fronts,
    }
    common.write_json_report(out_path, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--column-thickness", action="store_true")
    args = parser.parse_args()
    report = run(args.save, args.out, args.column_thickness)
    for front in report["fronts"]:
        print(front["side"], "h", front["thickness_m"], "w", front["width_m"],
              "v", front["convergence_speed_m_per_s"], front["post_transfer_convergence_speed_m_per_s"])
        for label, steps in front["cases"].items():
            for step, rows in steps.items():
                first, last = rows[0], rows[-1]
                print(f"  {label:18s} {step:8s} first F={first['hinge_force_n']:.3e} "
                      f"F/W={first['hinge_over_weight']:.3g} align={first['down_dip_alignment']:+.2f} | "
                      f"last F={last['hinge_force_n']:.3e} F/W={last['hinge_over_weight']:.3g} "
                      f"N/m={last['hinge_per_trench_length_n_per_m']:.3e} "
                      f"align={last['down_dip_alignment']:+.2f} seg1={last['first_segment_length_m'] / 1e3:.0f}km")


if __name__ == "__main__":
    main()

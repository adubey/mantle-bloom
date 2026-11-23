"""What the pulled trenches are made of (#995): relative-motion direction and floor age at plate-boundary nodes.

Steps one generated world plates-only under one RFC 0021 setting and, at every `--every` Myr, over every plate's
boundary-band nodes (length-weighted by `boundary_node_length_m`) records:

- the relative velocity of the plate against its nearest neighbour node, split into the part along
  `direction_to_neighbor` (the classifier's "closing rate") and the rest (tangential). Shares of band length whose
  tangential speed exceeds the normal speed, and whose |normal| is under the 1 cm/yr transform threshold;
- the pulled-trench nodes' (`subducting_boundary_mask`) own floor age (`node_created_years`) in bins, and the
  force-column Hm against the cell's own Hm.

Usage (repo root; compute-budget first; plates-only, about 1 GB):
    backend/.venv/bin/python bin/debug/boundary_geometry_probe.py --settings analysis/issue975/settings.json \\
        --setting b123-brake-off --seed 591453549 --years 200 --every 50 --out analysis/issue995/boundary-geometry-<seed>.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "backend")]

import plate_probe_common as common  # noqa: E402

AGE_BINS = (0.0, 10.0, 30.0, 60.0, 100.0, 1e9)


def measure(world) -> dict:
    from app import torque
    from app.surface_fields import SURFACE_FIELDS

    geom, motions = common.boundary_motion(world, require_neighbour=False)
    rows, per_plate = [], []
    for m in motions:
        if m.length is None:
            continue
        plate, inputs, band, length, normal = m.plate, m.inputs, m.band, m.length, m.normal
        tangential = np.linalg.norm(m.v_rel - normal[:, None] * inputs.direction_to_neighbor, axis=-1)
        sub = torque.subducting_boundary_mask(plate, inputs, geom.reach)
        created = plate.collect("node_created_years")
        undated = created == SURFACE_FIELDS["node_created_years"].sentinel
        age = np.where(undated, np.nan, (world.elapsed_years - created) / 1e6)
        v_mag = np.hypot(normal, tangential)
        per_plate.append({"band_km": float(length[band].sum() / 1e3), "area_km2": float(plate.node_areas_m2().sum() / 1e6),
                          "pulled_km": {"today": float(length[band & sub].sum() / 1e3),
                                        **{str(t): float(length[band & sub & (normal > np.sin(np.radians(t)) * v_mag)].sum() / 1e3) for t in (15, 30, 45)}}})
        rows.append({"band": band, "length": length, "normal": normal, "tangential": tangential, "sub": sub,
                     "age": age, "own_hm": inputs.own_hm, "force_hm": inputs.column_hm})
    cat = {k: np.concatenate([r[k] for r in rows]) for k in rows[0]}
    band, length = cat["band"], cat["length"]
    w = length[band]
    speed = np.hypot(cat["normal"][band], cat["tangential"][band])
    out = {"band_km": float(w.sum() / 1e3),
           "share_tangential_exceeds_normal": float(w[cat["tangential"][band] > np.abs(cat["normal"][band])].sum() / w.sum()),
           "share_abs_normal_under_1cm": float(w[np.abs(cat["normal"][band]) < 1.0].sum() / w.sum()),
           "share_closing": float(w[cat["normal"][band] > 0].sum() / w.sum()),
           "share_closing_past_1cm": float(w[cat["normal"][band] > 1.0].sum() / w.sum()),
           "median_relative_speed_cm_per_yr": float(np.median(speed))}
    # An angle test in place of the 1 cm/yr threshold: transform where |normal| < sin(theta) |v_rel|. `pulled_*` keeps
    # today's gate (closing past 1 cm/yr, down-going crust) and asks what share of it an angle test would also keep.
    normal_b, speed_b, sub_b = cat["normal"][band], speed, cat["sub"][band]
    out["angle_tests"] = {}
    for theta in (15, 30, 45):
        limit = np.sin(np.radians(theta)) * speed_b
        out["angle_tests"][str(theta)] = {
            "transform_share_of_band": float(w[np.abs(normal_b) < limit].sum() / w.sum()),
            "closing_share_of_band": float(w[normal_b > limit].sum() / w.sum()),
            "pulled_km_under_angle_and_gate": float(w[sub_b & (normal_b > limit)].sum() / 1e3),
            "share_of_todays_pulled_kept": float(w[sub_b & (normal_b > limit)].sum() / w[sub_b].sum()) if w[sub_b].sum() else None}
    out["per_plate_pulled_share_of_own_band"] = {}
    for key in ("today", "15", "30", "45"):
        shares = sorted(p["pulled_km"][key] / p["band_km"] for p in per_plate if p["band_km"] > 0)
        out["per_plate_pulled_share_of_own_band"][key] = {
            "plates": len(shares), "p10": shares[int(0.1 * (len(shares) - 1))], "median": shares[len(shares) // 2], "p90": shares[int(0.9 * (len(shares) - 1))],
            "under_5pct": sum(x < 0.05 for x in shares), "under_10pct": sum(x < 0.10 for x in shares)}
    sub = band & cat["sub"]
    ws, ages = length[sub], cat["age"][sub]
    if ws.sum():
        known = np.isfinite(ages)
        out["pulled_km"] = float(ws.sum() / 1e3)
        out["pulled_age_undated_share"] = float(ws[~known].sum() / ws.sum())
        out["pulled_age_bins_share"] = {f"{lo:g}-{hi:g}": float(ws[known & (ages >= lo) & (ages < hi)].sum() / ws.sum()) for lo, hi in zip(AGE_BINS, AGE_BINS[1:])}
        out["pulled_force_hm_over_own_hm_gt_1.05"] = float(ws[cat["force_hm"][sub] > 1.05 * cat["own_hm"][sub]].sum() / ws.sum())
        out["pulled_force_hm_mean_km"] = float(np.average(cat["force_hm"][sub], weights=ws) / 1e3)
        out["pulled_own_hm_mean_km"] = float(np.average(cat["own_hm"][sub], weights=ws) / 1e3)
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--settings", type=Path, required=True)
    p.add_argument("--setting", required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--years", type=float, default=200.0)
    p.add_argument("--step-myr", type=float, default=1.0)
    p.add_argument("--every", type=float, default=50.0)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--scratch", default="/tmp/unused", help="directory step_all may write to")
    args = p.parse_args()
    setting, world = common.open_world(args.settings, args.setting, args.seed)
    _, samples = common.step_and_sample(world, setting, args, measure, lambda sample: print(json.dumps(sample), flush=True))
    args.out.write_text(json.dumps({"seed": args.seed, "setting": args.setting, "samples": samples}, indent=1) + "\n")


if __name__ == "__main__":
    main()

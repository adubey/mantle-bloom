"""Why T2 fell at the best-case corners (#995): T2's two groups, their sizes, speeds and what drives each.

Re-reads breakdown.json of RFC 0021's G1 runs (no new run) and repeats `calibrate_plate_forces.t2_ratio`'s pairing
(plate speed T2_LEAD_MYR after the sample; slab-attached = consumed >= share of area), then reports per setting the
number of plate-samples and area share of each group, each group's area-weighted mean speed, and for each group the
area-weighted mean torque along the plate's axis (slab pull, ridge push, mantle flow, basal drag, slab drag), the
share of plates that have any subducting node, and mean plate area.

Usage: backend/.venv/bin/python bin/debug/g1_t2_groups.py --out analysis/issue995/t2-groups.json LABEL=DIR [LABEL=DIR ...]
where DIR holds <seed>/breakdown.json.
"""

import argparse
import json
from pathlib import Path

LEAD, SHARE = 10.0, 0.10
TORQUES = ("slab_pull", "ridge_push", "mantle_flow", "basal_drag", "slab_drag", "blocked_edges_residual")


def wmean(pairs):
    w = sum(b for _, b in pairs)
    return sum(a * b for a, b in pairs) / w if w else None


def groups(directory: Path) -> dict:
    out = {"attached": [], "other": []}
    for path in sorted(directory.glob("*/breakdown.json")):
        d = json.loads(path.read_text())
        by = {round(s["myr"], 6): {r["plate_id"]: r for r in s["plates"]} for s in d["samples"]}
        for myr, rows in by.items():
            later = by.get(round(myr + LEAD, 6))
            for pid, row in (rows.items() if later else ()):
                if pid in later:
                    key = "attached" if row["g0"]["consumed_10myr_km2"] >= SHARE * row["area_km2"] else "other"
                    out[key].append((row, later[pid]))
    return out


def describe(pairs: list, total_area: float) -> dict:
    rows = [r for r, _ in pairs]
    desc = {"plate_samples": len(rows), "area_share": sum(r["area_km2"] for r in rows) / total_area,
            "mean_area_Mkm2": sum(r["area_km2"] for r in rows) / len(rows) / 1e6 if rows else None,
            "speed_later_cm_per_yr": wmean([(l["speed_cm_per_yr"]["all"], r["area_km2"]) for r, l in pairs]),
            "speed_now_cm_per_yr": wmean([(r["speed_cm_per_yr"]["all"], r["area_km2"]) for r in rows]),
            "mantle_flow_only_cm_per_yr": wmean([(r["counterfactual_speed_cm_per_yr"]["mantle_flow_only"]["all"], r["area_km2"]) for r in rows]),
            "has_subducting_node": wmean([(float(r["g0"]["has_subducting_node"]), r["area_km2"]) for r in rows]),
            "continental_share": wmean([(r["continental_area_share"], r["area_km2"]) for r in rows]),
            "slab_km_per_Mkm2": wmean([(r["g0"]["geometric_slab_km"] / (r["area_km2"] / 1e6), r["area_km2"]) for r in rows])}
    for t in TORQUES:
        desc[f"torque_{t}_1e25Nm"] = wmean([(r["torque_along_axis_nm"][t] / 1e25, r["area_km2"]) for r in rows])
    return desc


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path)
    p.add_argument("settings", nargs="+")
    args = p.parse_args()
    result = {}
    for item in args.settings:
        label, _, directory = item.partition("=")
        g = groups(Path(directory).expanduser())
        total = sum(r["area_km2"] for pairs in g.values() for r, _ in pairs)
        result[label] = {k: describe(v, total) for k, v in g.items()}
        a, o = result[label]["attached"], result[label]["other"]
        print(f"{label:28s} attached {a['plate_samples']:4d} ({a['area_share']:.0%}) {a['speed_later_cm_per_yr']:.2f} cm/yr | "
              f"other {o['plate_samples']:4d} ({o['area_share']:.0%}) {o['speed_later_cm_per_yr']:.2f} cm/yr | ratio {a['speed_later_cm_per_yr'] / o['speed_later_cm_per_yr']:.2f}")
    if args.out:
        args.out.write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()

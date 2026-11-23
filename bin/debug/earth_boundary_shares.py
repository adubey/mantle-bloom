"""Earth's plate-boundary geometry from PB2002, to set against the model's (#995).

Usage (the steps file is Bird's, downloaded beforehand into a directory of its own from
http://peterbird.name/oldFTP/PB2002/; plate areas and speeds come from `analysis/issue958/earth-plate-speed.json`):
    python3 -I bin/debug/earth_boundary_shares.py --steps <dir>/PB2002_steps.dat.txt \\
        --plates analysis/issue958/earth-plate-speed.json [--out analysis/issue995/earth-boundary-shares.json]

PB2002_steps.dat gives each ~20-100 km boundary step's plate pair, length, relative speed (mm/yr) and its divergence
and left-lateral components, and a class (`SUB` subduction zone; `OSR`, `OTF`, `OCB`, `CCB`, ...). A subduction step's
pair is written `A\\B` (A subducts) or `A/B` (B subducts), checked here against Kermadec, Aleutians, Kuril and Cascadia.
Per plate: its perimeter is every step it belongs to; its pulled length is the SUB steps where it is the subducting
plate. Reports the same measures as `boundary_geometry_probe.py`: the share of perimeter that is pulled trench, the
spread across plates, and the share of boundary that slides more than it closes.
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_steps(path: Path) -> list[dict]:
    steps = []
    for line in path.read_text().splitlines():
        t = line.split()
        if len(t) < 14:
            continue
        pair = t[1].lstrip(":")
        cls = t[-1].lstrip(":").rstrip("*")
        steps.append({"pair": pair, "class": cls, "length_km": float(t[6]), "speed_mm": float(t[8]),
                      "divergence_mm": float(t[10]), "left_lateral_mm": float(t[11])})
    return steps


def split_pair(pair: str) -> tuple[str, str, str | None]:
    for sep in ("\\", "/", "-"):
        if sep in pair:
            a, b = pair.split(sep)
            return a, b, sep
    raise ValueError(pair)


def wq(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    order = np.argsort(values)
    cum = np.cumsum(weights[order]) / weights.sum()
    return float(values[order][np.searchsorted(cum, q)])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--steps", type=Path, required=True)
    p.add_argument("--plates", type=Path, required=True)
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    steps = read_steps(args.steps)
    areas = json.loads(args.plates.read_text())["plates"]
    perimeter, pulled = defaultdict(float), defaultdict(float)
    for s in steps:
        a, b, sep = split_pair(s["pair"])
        perimeter[a] += s["length_km"]
        perimeter[b] += s["length_km"]
        if s["class"] == "SUB":
            pulled[a if sep == "\\" else b] += s["length_km"]
    w = np.array([s["length_km"] for s in steps])
    div = np.array([s["divergence_mm"] for s in steps])
    ll = np.array([s["left_lateral_mm"] for s in steps])
    speed = np.array([s["speed_mm"] for s in steps])
    out = {
        "boundary_km_each_counted_once": float(w.sum()),
        "subduction_km": float(sum(pulled.values())),
        "subduction_share_of_boundary": float(sum(pulled.values()) / w.sum()),
        "share_sliding_more_than_closing_or_opening": float(w[np.abs(ll) > np.abs(div)].sum() / w.sum()),
        "share_closing_past_10mm_per_yr": float(w[-div > 10.0].sum() / w.sum()),
        "share_abs_normal_under_10mm_per_yr": float(w[np.abs(div) < 10.0].sum() / w.sum()),
        "median_relative_speed_cm_per_yr": wq(speed, w, 0.5) / 10.0,
    }
    rows = [(k, areas[k]["area_share"], pulled[k] / perimeter[k], pulled[k], areas[k]["nnr_cm_per_yr"]) for k in areas if perimeter.get(k)]
    shares = np.array([r[2] for r in rows])
    area = np.array([r[1] for r in rows])
    out["plates"] = {k: {"area_share": a, "perimeter_km": perimeter[k], "pulled_km": pulled[k], "pulled_share_of_perimeter": s, "nnr_cm_per_yr": v}
                     for k, a, s, _, v in sorted(rows, key=lambda r: -r[1])}
    out["per_plate_pulled_share"] = {
        "plates": len(rows), "unweighted_p10": float(np.quantile(shares, 0.1)), "unweighted_median": float(np.median(shares)),
        "unweighted_p90": float(np.quantile(shares, 0.9)), "area_weighted_p10": wq(shares, area, 0.1), "area_weighted_median": wq(shares, area, 0.5),
        "area_weighted_p90": wq(shares, area, 0.9), "under_5pct": int((shares < 0.05).sum()), "under_10pct": int((shares < 0.10).sum()),
        "area_share_under_5pct": float(area[shares < 0.05].sum()), "area_share_under_10pct": float(area[shares < 0.10].sum()),
        "area_share_over_30pct": float(area[shares > 0.30].sum()),
    }
    print(json.dumps({k: v for k, v in out.items() if k != "plates"}, indent=1))
    for k, v in list(out["plates"].items())[:16]:
        print(f"{k:3s} area {v['area_share']:.3f}  perimeter {v['perimeter_km']:7.0f} km  pulled {v['pulled_km']:6.0f} km  share {v['pulled_share_of_perimeter']:.2f}  speed {v['nnr_cm_per_yr']:.1f}")
    if args.out:
        args.out.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

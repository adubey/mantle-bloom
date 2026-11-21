"""Are the two plates of a contact equally fast? Earth against model samples (#998).

For each contact (plate pair) and class (converging, diverging, sliding, by the obliquity rule of
`earth_speed_reference.py`) the asymmetry is min(speed_a, speed_b) / max(speed_a, speed_b), with each plate's
area-mean no-net-rotation speed. Near 0: one plate moves, the other is nearly still (a slab-pulled plate diving
under a stationary overriding plate). Near 1: both move alike. Reported length-weighted per class, with the share
of contact length with asymmetry below 0.25 (one-sided) and above 0.5 (two-sided), and the median of the faster
and slower plate's speed. The direction-sensitive measure is `opposition` = the contact's length-weighted relative
speed over the sum of the two plates' speeds: about 1 when the plates move in opposite directions across the
contact, near |sa - sb| / (sa + sb) when they move alike. Plate speeds are per plate, so a contact's two plates are the same two numbers at every
node of the contact; only the class split is by node.

    python3 bin/debug/contact_asymmetry.py --steps <captures>/pb2002/PB2002_steps.dat.txt \\
        --plates analysis/issue958/earth-plate-speed.json \\
        --runs b123-brake-off-100myr@100 b123-100myr@100 --captures <captures> --out analysis/issue998/contact-asymmetry.json

`dir@myr` picks a sample time of every seed in `<captures>/<dir>`. numpy only.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import earth_speed_reference as earth  # noqa: E402

CLASSES = ("converging", "diverging", "sliding")
MIN_CONTACT_KM = earth.MIN_CONTACT_KM


def summarize(rows: list[tuple]) -> dict:
    """rows: (class, km, asymmetry, faster plate speed, slower plate speed, opposition)."""
    out = {}
    for cls in CLASSES:
        sel = [r for r in rows if r[0] == cls]
        if not sel:
            continue
        km = np.array([r[1] for r in sel])
        a = np.array([r[2] for r in sel])
        order = np.argsort(a)
        cum = (np.cumsum(km[order]) - 0.5 * km[order]) / km.sum()
        out[cls] = {"contacts": len(sel), "length_km": float(km.sum()),
                    "asymmetry_median": float(np.interp(0.5, cum, a[order])),
                    "share_one_sided_under_0.25": float(km[a < 0.25].sum() / km.sum()),
                    "share_two_sided_over_0.5": float(km[a > 0.5].sum() / km.sum()),
                    "faster_plate_speed_median": float(np.median([r[3] for r in sel])),
                    "slower_plate_speed_median": float(np.median([r[4] for r in sel])),
                    "opposition_median": float(np.median([r[5] for r in sel])),
                    "share_opposition_over_0.8": float(km[np.array([r[5] for r in sel]) > 0.8].sum() / km.sum())}
    return out


def earth_rows(steps: Path, plates: Path) -> list[tuple]:
    speed = {k: v["nnr_cm_per_yr"] for k, v in json.loads(plates.read_text())["plates"].items()}
    length = defaultdict(float)
    moved = defaultdict(float)
    for r in earth.read_steps(steps):
        a, b = r["pair"]
        if a in speed and b in speed:
            length[(a, b, r["cls"])] += r["km"]
            moved[(a, b, r["cls"])] += r["km"] * r["speed"]
    total = defaultdict(float)
    for (a, b, _), km in length.items():
        total[(a, b)] += km
    return [(cls, km, min(speed[a], speed[b]) / max(speed[a], speed[b]), max(speed[a], speed[b]), min(speed[a], speed[b]),
             moved[(a, b, cls)] / km / (speed[a] + speed[b]))
            for (a, b, cls), km in length.items() if total[(a, b)] >= MIN_CONTACT_KM]


def model_rows(sample: dict) -> list[tuple]:
    speed = {str(p["plate_id"]): p["speed_nnr_cm_per_yr"] for p in sample["breakdown"]["plates"]}
    rows = []
    for contact in sample["contacts"]:
        a, b = str(contact["a"]), str(contact["b"])
        nodes = np.array(contact["nodes"])
        if len(nodes) == 0 or nodes[:, 0].sum() < MIN_CONTACT_KM or a not in speed or b not in speed:
            continue
        sa, sb = speed[a], speed[b]
        if max(sa, sb) <= 0:
            continue
        cls = earth_classes(nodes[:, 1], nodes[:, 2])
        for c in CLASSES:
            km = nodes[cls == c, 0].sum()
            if km > 0:
                m = cls == c
                rows.append((c, float(km), min(sa, sb) / max(sa, sb), max(sa, sb), min(sa, sb),
                             float((nodes[m, 0] * nodes[m, 1]).sum() / km / (sa + sb))))
    return rows


def earth_classes(speed: np.ndarray, closing: np.ndarray) -> np.ndarray:
    return np.where(np.abs(closing) < earth.SLIDING_OBLIQUITY * speed, "sliding", np.where(closing > 0, "converging", "diverging"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", type=Path, required=True)
    ap.add_argument("--plates", type=Path, required=True)
    ap.add_argument("--captures", type=Path, required=True)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    result = {"Earth (PB2002)": summarize(earth_rows(args.steps, args.plates))}
    for spec in args.runs:
        name, _, myr = spec.partition("@")
        rows = []
        for path in sorted((args.captures / name).glob("*.json")):
            for sample in json.loads(path.read_text())["samples"]:
                if float(myr) == sample["myr"]:
                    rows += model_rows(sample)
        result[spec] = summarize(rows)
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    for name, block in result.items():
        print(name)
        for cls, b in block.items():
            print(f"  {cls:11s} contacts {b['contacts']:4d}  asym median {b['asymmetry_median']:.2f}  one-sided(<0.25) {b['share_one_sided_under_0.25']:.2f}  "
                  f"two-sided(>0.5) {b['share_two_sided_over_0.5']:.2f}  opposition median {b['opposition_median']:.2f} over 0.8: {b['share_opposition_over_0.8']:.2f}")


if __name__ == "__main__":
    main()

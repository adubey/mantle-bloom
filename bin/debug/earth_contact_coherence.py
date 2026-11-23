"""How alike do Earth's plates move across each contact, compared with the model's? (#998, exploration)

Earth side of `contact_flow_coherence.py`'s `actual` numbers. For every PB2002 boundary segment (steps file) the two
plates' own velocities at the segment midpoint come from Bird's Pacific-fixed Euler poles with the area-weighted net
rotation removed (no-net-rotation frame, the way `earth_plate_speed.py` does it), and the contact's opposition is
sum(km x |v_a - v_b|) / sum(km x (|v_a| + |v_b|)) over the pair's segments: about 1 = the plates move opposite ways at the
contact, near 0 = alike. Contacts of at least 500 km, as in the model. The model side is the same quantity from the
`cfc` captures. The model's own net rotation is 0.3 cm/yr against plate speeds of 4-8, so its frame is close to NNR.

    python3 -I bin/debug/earth_contact_coherence.py --pb2002 <captures>/pb2002 --cfc <captures>/cfc --out analysis/issue998/earth-contact-coherence.json

`<captures>/pb2002` holds PB2002_steps.dat.txt, PB2002_poles.dat.txt and PB2002_plates.dig.txt (see `earth_plate_speed.py`).
numpy only.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import earth_plate_speed as eps  # noqa: E402
import earth_speed_reference as ref  # noqa: E402

MIN_CONTACT_KM = ref.MIN_CONTACT_KM
GRID_POINTS = 80_000


def net_rotation(pb: Path, poles: dict) -> np.ndarray:
    outlines = eps.read_outlines(pb / "PB2002_plates.dig.txt")
    points = eps.fibonacci(GRID_POINTS)
    owner = np.full(len(points), -1)
    claims = np.zeros(len(points), dtype=int)
    names = sorted(outlines)
    for k, name in enumerate(names):
        inside = eps.winding(points, outlines[name]) == 1
        owner[inside] = k
        claims += inside
    ok = claims == 1
    omega = np.array([poles.get(name, np.zeros(3)) for name in names])[owner]
    pts, om = points[ok], omega[ok]
    projector = np.eye(3)[None] - pts[:, :, None] * pts[:, None, :]
    return np.linalg.solve(projector.sum(axis=0), np.cross(pts, np.cross(om, pts)).sum(axis=0))


def earth_contacts(pb: Path) -> list[dict]:
    poles = eps.read_poles(pb / "PB2002_poles.dat.txt")
    net = net_rotation(pb, poles)
    acc = defaultdict(lambda: defaultdict(float))
    for line in (pb / "PB2002_steps.dat.txt").read_text().splitlines():
        t = line.split()
        if len(t) < 14:
            continue
        pair = t[1].lstrip(":")
        for sep in ("\\", "/", "-"):
            if sep in pair:
                a, b = pair.split(sep)
                break
        else:
            continue
        lon = np.radians((float(t[2]) + float(t[4])) / 2.0) if abs(float(t[2]) - float(t[4])) < 180 else None
        if lon is None:  # a segment across the dateline: use the start point
            lon = np.radians(float(t[2]))
            lat = np.radians(float(t[3]))
        else:
            lat = np.radians((float(t[3]) + float(t[5])) / 2.0)
        r = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
        va = np.cross(poles.get(a, np.zeros(3)) - net, r) * eps.CM_PER_YR_PER_DEG_PER_MYR
        vb = np.cross(poles.get(b, np.zeros(3)) - net, r) * eps.CM_PER_YR_PER_DEG_PER_MYR
        km = float(t[6])
        key = tuple(sorted((a, b)))
        acc[key]["km"] += km
        acc[key]["rel"] += km * float(np.linalg.norm(va - vb))
        acc[key]["sum"] += km * float(np.linalg.norm(va) + np.linalg.norm(vb))
    return [{"pair": k, "km": v["km"], "opposition": v["rel"] / v["sum"]} for k, v in acc.items() if v["km"] >= MIN_CONTACT_KM and v["sum"] > 0]


def summarize(rows: list[tuple[float, float]]) -> dict:
    """rows: (km, opposition)."""
    km = np.array([r[0] for r in rows])
    op = np.array([r[1] for r in rows])
    order = np.argsort(op)
    cum = (np.cumsum(km[order]) - 0.5 * km[order]) / km.sum()
    return {"contacts": len(rows), "length_km": float(km.sum()), "opposition_length_weighted_mean": float((km * op).sum() / km.sum()),
            "opposition_median": float(np.interp(0.5, cum, op[order])),
            "share_alike_under_0.3": float(km[op < 0.3].sum() / km.sum()), "share_middle": float(km[(op >= 0.3) & (op <= 0.7)].sum() / km.sum()),
            "share_opposed_over_0.7": float(km[op > 0.7].sum() / km.sum())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pb2002", type=Path, required=True)
    ap.add_argument("--cfc", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    earth = earth_contacts(args.pb2002)
    result = {"Earth (PB2002, NNR)": summarize([(c["km"], c["opposition"]) for c in earth])}
    runs = [json.loads(p.read_text()) for p in sorted(args.cfc.glob("*.json"))]
    for t in (5.0, 20.0, 50.0):
        rows = [(c["km"], c["actual"]["opposition"]) for r in runs for s in r["samples"] if s["myr"] == t for c in s["contacts"] if c["actual"]["opposition"] is not None]
        result[f"model flow x4, brake 0.5, {t:g} Myr"] = summarize(rows)
    flow_rows = [(c["km"], c["flow"]["opposition"]) for r in runs for s in r["samples"] if s["myr"] == 5.0 for c in s["contacts"] if c["flow"]["opposition"] is not None]
    result["model flow alone, 5 Myr geometry"] = summarize(flow_rows)
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    for name, b in result.items():
        print(f"{name:42s} contacts {b['contacts']:4d}  mean {b['opposition_length_weighted_mean']:.2f}  median {b['opposition_median']:.2f}  "
              f"alike<0.3 {b['share_alike_under_0.3']:.2f}  middle {b['share_middle']:.2f}  opposed>0.7 {b['share_opposed_over_0.7']:.2f}")


if __name__ == "__main__":
    main()

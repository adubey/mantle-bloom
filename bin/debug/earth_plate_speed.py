"""Earth's area-weighted mean plate speed from a published plate model, for RFC 0021's T1 (#958).

Usage (any Python with numpy; the two PB2002 files are Bird's published plate model, downloaded beforehand into a
directory of their own from http://peterbird.name/oldFTP/PB2002/):
    python3 -I bin/debug/earth_plate_speed.py --plates <dir>/PB2002_plates.dig.txt --poles <dir>/PB2002_poles.dat.txt \
        [--points 200000] [--out analysis/issue958/earth-plate-speed.json]

Bird's PB2002 [Bird, 2003, doi:10.1029/2001GC000252] gives 52 plate outlines and Euler poles relative to the Pacific
plate. Each point of a near-uniform Fibonacci grid is assigned to the plate whose outline winds around it (winding
number of the outline, by the signed angles of its edges seen from the point; this is the stereographic projection
from the point's antipode, so it holds for outlines around a pole and across the dateline). The no-net-rotation
frame removes the area-weighted net rotation `(∫ (I - r rᵀ) dA)⁻¹ ∫ r × v dA`, and in that frame the script reports
the area-weighted mean of `|ω × r|`, overall and per plate, and the same in the Pacific-fixed frame for comparison.
"""

import argparse
import json
from pathlib import Path

import numpy as np

EARTH_RADIUS_CM = 6.371e8
CM_PER_YR_PER_DEG_PER_MYR = np.radians(1.0) * EARTH_RADIUS_CM / 1e6


def unit(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    return np.column_stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)))


def read_outlines(path: Path) -> dict[str, np.ndarray]:
    outlines, name, points = {}, None, []
    for line in path.read_text().splitlines():
        text = line.strip()
        if not text:
            continue
        if text.startswith("***"):
            outlines[name] = unit(np.array([p[1] for p in points]), np.array([p[0] for p in points]))
            name, points = None, []
        elif name is None:
            name = text
        else:
            lon, lat = (float(v) for v in text.split(","))
            points.append((lon, lat))
    return outlines


def read_poles(path: Path) -> dict[str, np.ndarray]:
    poles = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 4 and len(parts[0]) == 2:
            lat, lon, rate = float(parts[1]), float(parts[2]), float(parts[3])
            poles[parts[0]] = unit(np.array([lat]), np.array([lon]))[0] * rate  # deg/Myr, CCW positive
    return poles


def fibonacci(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    z = 1.0 - 2.0 * i / n
    lon = np.pi * (1.0 + 5.0**0.5) * i
    r = np.sqrt(1.0 - z * z)
    return np.column_stack((r * np.cos(lon), r * np.sin(lon), z))


def winding(points: np.ndarray, outline: np.ndarray, chunk: int = 1000) -> np.ndarray:
    """Per point, the outline's winding number around it (0 outside; see main for the sign)."""
    a, b = outline, np.roll(outline, -1, axis=0)
    out = np.zeros(len(points))
    for start in range(0, len(points), chunk):
        p = points[start:start + chunk][:, None, :]
        ua = a[None] - np.sum(a[None] * p, axis=2, keepdims=True) * p
        ub = b[None] - np.sum(b[None] * p, axis=2, keepdims=True) * p
        cross = np.sum(np.cross(ua, ub) * p, axis=2)
        dot = np.sum(ua * ub, axis=2)
        out[start:start + chunk] = np.arctan2(cross, dot).sum(axis=1) / (2.0 * np.pi)
    return np.rint(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--plates", type=Path, required=True)
    parser.add_argument("--poles", type=Path, required=True)
    parser.add_argument("--points", type=int, default=200_000)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    outlines, poles = read_outlines(args.plates), read_poles(args.poles)
    points = fibonacci(args.points)
    owner = np.full(len(points), -1)
    names = sorted(outlines)
    claims = np.zeros(len(points), dtype=int)
    for k, name in enumerate(names):
        # Seen from a point, an outline winds +1 around it when the point is inside and -1 when its antipode is
        # (the projection is from the antipode). PB2002's outlines all run the same way, inside +1, which the
        # exactly-once coverage check below confirms.
        inside = winding(points, outlines[name]) == 1
        owner[inside] = k
        claims += inside
    assigned = claims == 1
    omega = np.array([poles[name] for name in names])[owner]  # deg/Myr, Pacific-fixed
    pts, om = points[assigned], omega[assigned]
    velocity = np.cross(om, pts)
    projector = np.eye(3)[None] - pts[:, :, None] * pts[:, None, :]
    net = np.linalg.solve(projector.sum(axis=0), np.cross(pts, velocity).sum(axis=0))
    nnr = np.linalg.norm(np.cross(om - net, pts), axis=1) * CM_PER_YR_PER_DEG_PER_MYR
    pa = np.linalg.norm(velocity, axis=1) * CM_PER_YR_PER_DEG_PER_MYR
    owners = owner[assigned]
    per_plate = {names[k]: {"area_share": float(np.mean(owners == k)), "nnr_cm_per_yr": float(nnr[owners == k].mean())}
                 for k in range(len(names)) if np.any(owners == k)}
    result = {
        "model": "PB2002 (Bird 2003, doi:10.1029/2001GC000252), Pacific-fixed poles, NNR by removing the net rotation",
        "grid_points": args.points, "unassigned_share": float(np.mean(claims == 0)), "multiply_assigned_share": float(np.mean(claims > 1)),
        "net_rotation_deg_per_myr": float(np.linalg.norm(net)),
        "mean_speed_nnr_cm_per_yr": float(nnr.mean()), "rms_speed_nnr_cm_per_yr": float(np.sqrt(np.mean(nnr**2))),
        "mean_speed_pacific_fixed_cm_per_yr": float(pa.mean()),
        "plates": dict(sorted(per_plate.items(), key=lambda kv: -kv[1]["area_share"])),
    }
    text = json.dumps(result, indent=1) + "\n"
    if args.out:
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()

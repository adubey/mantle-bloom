#!/usr/bin/env python3
"""Where do a generated world's plate boundaries run, relative to its continents? (RFC 0009)

On Earth, plate boundaries mostly run through ocean; this measures how far a generated
world departs from that at t = 0. On a quarter-degree equirectangular grid (each pixel
takes its nearest node), it reports:
- plate-boundary length, and the share with land on both sides, land on one side, or none;
  for generated worlds also the share with continental crust on both sides;
- how much of the continent's edge (continental crust next to oceanic) lies within
  EDGE_NEAR_KM of a plate boundary: on Earth most continental margins are passive, inside
  a plate;
- how mixed plates are: the share of plates, and of plate area, whose continental share is
  5-95%;
- for cratons, the share of craton area within 500 and 1,000 km of a plate boundary.

`--earth` measures Earth itself: the 16 real plates of `data/plate_owner_grid.npz`
(NNR-MORVEL56 outlines) against the real coastline (`data/real_land_mask.npz`), the inputs
the premade "earth" world is built from. Not that world itself: the model's plates have one
crust type each, so it splits every real plate into land and ocean sub-plates. Land stands
in for continental crust there and leaves out the drowned shelves, and there are no craton
measures.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_plate_layout.py \
        --seeds 7,11,2024 --earth --out /tmp/layout.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_plate_layout.py \
        --seeds "" --load seed7-500myr.mbworld --out /tmp/layout-500.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from app import cratons, persistence
from app import world as world_mod
from app.elevation_lines import PLANET_RADIUS_KM, effective_is_continental_from_codes
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree

GRID_DEG = 0.25
MIXED_RANGE = (0.05, 0.95)
NEAR_KM = (500.0, 1000.0)
# An active margin's trench sits offshore of the coast, ~100-300 km at the Andes.
EDGE_NEAR_KM = 300.0


def _grid(world):
    """Per pixel: plate index, land, continental crust, craton."""
    points, plate, land, cont, craton = [], [], [], [], []
    for n, p in enumerate(q for q in world.plates if q.node_count() > 0):
        pts, elev = p.all_points_and_elevation()
        points.append(pts)
        plate.append(np.full(len(pts), n))
        land.append(elev > world.sea_level_m)
        cont.append(effective_is_continental_from_codes(p.collect("crust_type_code"), p.crust_type == "continental"))
        craton.append(p.collect("craton_crust_m") > 0.0)
    points = np.concatenate(points)
    h, w = int(180 / GRID_DEG), int(360 / GRID_DEG)
    lat = np.radians(90.0 - (np.arange(h) + 0.5) * GRID_DEG)
    lon = np.radians(-180.0 + (np.arange(w) + 0.5) * GRID_DEG)
    lon2, lat2 = np.meshgrid(lon, lat)
    xyz = np.stack([np.cos(lat2) * np.cos(lon2), np.cos(lat2) * np.sin(lon2), np.sin(lat2)], axis=-1).reshape(-1, 3)
    _, idx = cKDTree(points).query(xyz, workers=-1)
    fields = [np.concatenate(f)[idx].reshape(h, w) for f in (plate, land, cont, craton)]
    return (*fields, lat)


def _distance_km(mask: np.ndarray, lat: np.ndarray, step_km: float) -> np.ndarray:
    """Approximate distance to the nearest `mask` pixel: pixel distance scaled by the mean of
    the north-south and east-west pixel sizes at that row -- fair away from the poles."""
    dist_px = distance_transform_edt(~mask)
    return dist_px * step_km * np.sqrt(0.5 * (1.0 + np.cos(lat)[:, None] ** 2))


def _earth_grid():
    """The real plates and coastline on the same grid; unowned cells take the nearest plate."""
    from app import real_plates

    owner, _ = real_plates.load_plate_owner_grid()
    fine = real_plates.load_real_land_mask()
    h, w = owner.shape
    land = fine[: (fine.shape[0] // h) * h, : (fine.shape[1] // w) * w].reshape(h, fine.shape[0] // h, w, fine.shape[1] // w).mean(axis=(1, 3)) >= 0.5
    _, (iy, ix) = distance_transform_edt(owner < 0, return_indices=True)
    plate = owner[iy, ix].astype(int)
    lat = np.radians(90.0 - (np.arange(h) + 0.5) * 180.0 / h)
    return plate, land, land, np.zeros_like(land), lat


def _measure(fields, label: str) -> dict:
    plate, land, cont, craton, lat = fields
    h, w = plate.shape
    step_km = np.radians(GRID_DEG) * PLANET_RADIUS_KM
    # A change between east-west neighbours is a boundary segment one pixel tall; between
    # north-south neighbours, one pixel wide at that latitude.
    ew = plate != np.roll(plate, -1, axis=1)
    ns = np.zeros_like(ew)
    ns[:-1] = plate[:-1] != plate[1:]
    ew_len = np.full((h, 1), step_km) * np.ones((1, w))
    ns_len = (step_km * np.cos(lat))[:, None] * np.ones((1, w))
    out = {"label": label, "plates": int(plate.max() + 1)}
    total = float(ew_len[ew].sum() + ns_len[ns].sum())
    out["boundary_km"] = total
    for name, field in (("land", land), ("continental", cont)):
        both_ew = field & np.roll(field, -1, axis=1)
        one_ew = field ^ np.roll(field, -1, axis=1)
        both_ns = np.zeros_like(field)
        one_ns = np.zeros_like(field)
        both_ns[:-1] = field[:-1] & field[1:]
        one_ns[:-1] = field[:-1] ^ field[1:]
        both = float(ew_len[ew & both_ew].sum() + ns_len[ns & both_ns].sum())
        one = float(ew_len[ew & one_ew].sum() + ns_len[ns & one_ns].sum())
        out[f"boundary_{name}_both_share"] = both / total
        out[f"boundary_{name}_one_share"] = one / total
    boundary = ew | ns
    dist_km = _distance_km(boundary, lat, step_km)
    cont_edge_ew = cont ^ np.roll(cont, -1, axis=1)
    cont_edge_ns = np.zeros_like(cont)
    cont_edge_ns[:-1] = cont[:-1] ^ cont[1:]
    near = dist_km <= EDGE_NEAR_KM
    edge_total = float(ew_len[cont_edge_ew].sum() + ns_len[cont_edge_ns].sum())
    edge_near = float(ew_len[cont_edge_ew & near].sum() + ns_len[cont_edge_ns & near].sum())
    out["continent_edge_km"] = edge_total
    out["continent_edge_on_boundary_share"] = edge_near / edge_total if edge_total else None
    area = np.cos(lat)[:, None] * np.ones((1, w))
    shares = np.array([cont[plate == p].mean() for p in range(plate.max() + 1)])
    plate_area = np.array([area[plate == p].sum() for p in range(plate.max() + 1)])
    mixed = (shares > MIXED_RANGE[0]) & (shares < MIXED_RANGE[1])
    out["mixed_plate_share"] = float(mixed.mean())
    out["mixed_plate_area_share"] = float(plate_area[mixed].sum() / plate_area.sum())
    out["land_area_share"] = float(area[land].sum() / area.sum())
    out["continental_area_share"] = float(area[cont].sum() / area.sum())
    if craton.any():
        c_area = area[craton].sum()
        for km in NEAR_KM:
            out[f"craton_within_{km:g}km_share"] = float(area[craton & (dist_km <= km)].sum() / c_area)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--load", type=Path, action="append", default=[], help="also measure a saved .mbworld (repeatable)")
    parser.add_argument("--earth", action="store_true", help="also measure the real plates and coastline")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for seed in map(int, args.seeds.split(",")) if args.seeds else ():
        world = world_mod.generate_world(seed)
        cratons.ensure_initialized(world)
        rows.append(_measure(_grid(world), f"seed {seed}"))
        print(json.dumps(rows[-1]), flush=True)
    for path in args.load:
        world = persistence.load_world_bytes(path.read_bytes())
        rows.append(_measure(_grid(world), f"{path.name} ({world.elapsed_years / 1e6:.0f} Myr)"))
        print(json.dumps(rows[-1]), flush=True)
    if args.earth:
        rows.append(_measure(_earth_grid(), "earth (16 real plates, real coastline)"))
        print(json.dumps(rows[-1]), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()

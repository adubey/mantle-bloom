#!/usr/bin/env python3
"""Earth's plates over its cratons: a reference figure for RFC 0009.

Draws the 16 real plate outlines the model bundles (`backend/app/data/plate_owner_grid.npz`,
NNR-MORVEL56; see real_plates.py) in red over the USGS world geologic provinces map, whose
shield (orange) and platform (pink) provinces are the cratons. The USGS map is public
domain but not bundled; download it first:

    curl -sSLo provinces.jpg \
        https://upload.wikimedia.org/wikipedia/commons/a/a9/World_geologic_provinces.jpg

Its projection isn't documented. It fits a Mercator projection, x = a + b lon and
y = c + d ln tan(45 + lat / 2), to within 8 pixels at the six coastline landmarks below,
read off the 1,200 x 637 image.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_earth_cratons.py \
        --provinces provinces.jpg --out /tmp/earth_plates_cratons.jpg
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from app import real_plates
from PIL import Image, ImageDraw
from scipy.ndimage import distance_transform_edt

# (lon, lat, x, y): Cape Agulhas, Kanyakumari, Cape Farewell, Cape Horn, South East Cape
# (Tasmania), Cape Guardafui.
LANDMARKS = np.array([
    (20.0, -34.8, 671, 517), (77.5, 8.1, 859.5, 360), (-43.9, 59.8, 451, 136),
    (-67.3, -56.0, 377, 614), (146.8, -43.6, 1089, 545), (51.3, 11.8, 782, 342),
])
SCALE = 2
# Saved at this width, as JPEG when --out ends in .jpg, to keep the checked-in figure small.
OUT_WIDTH = 1600
LINE = (200, 0, 0)


def _mercator_deg(lat: np.ndarray) -> np.ndarray:
    return np.degrees(np.log(np.tan(np.radians(45.0 + np.clip(lat, -85.0, 85.0) / 2.0))))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--provinces", type=Path, required=True, help="the USGS World_geologic_provinces.jpg")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    bx, ax = np.polyfit(LANDMARKS[:, 0], LANDMARKS[:, 2], 1)
    by, ay = np.polyfit(_mercator_deg(LANDMARKS[:, 1]), LANDMARKS[:, 3], 1)
    owner, _ = real_plates.load_plate_owner_grid()
    # Cells on an outline belong to no polygon (-1): give each its nearest plate.
    _, (iy, ix) = distance_transform_edt(owner < 0, return_indices=True)
    owner = owner[iy, ix]
    h, w = owner.shape
    lat = 90.0 - (np.arange(h) + 0.5) * 180.0 / h
    lon = -180.0 + (np.arange(w) + 0.5) * 360.0 / w
    # Boundary points: midpoints between neighbouring cells of different plates.
    ew = owner != np.roll(owner, -1, axis=1)
    ns = np.zeros_like(ew)
    ns[:-1] = owner[:-1] != owner[1:]
    iy, ix = np.nonzero(ew)
    pts = [np.c_[lon[ix] + 0.5 * 360.0 / w, lat[iy]]]
    iy, ix = np.nonzero(ns)
    pts.append(np.c_[lon[ix], lat[iy] - 0.5 * 180.0 / h])
    pts = np.concatenate(pts)
    x = ax + bx * pts[:, 0]
    y = ay + by * _mercator_deg(pts[:, 1])
    image = Image.open(args.provinces).convert("RGB")
    image = image.resize((image.size[0] * SCALE, image.size[1] * SCALE), Image.LANCZOS)
    draw = ImageDraw.Draw(image)
    for shift in (0.0, -360.0 * bx, 360.0 * bx):  # the map spans a little over 360 degrees
        for px, py in zip((x + shift) * SCALE, y * SCALE):
            if 0 <= px < image.size[0] and 0 <= py < image.size[1]:
                draw.ellipse((px - 1.5, py - 1.5, px + 1.5, py + 1.5), fill=LINE)
    draw.text((10, image.size[1] - 22), "Plates: NNR-MORVEL56 (16 largest). Provinces: USGS. Shield + platform = craton.", fill=(0, 0, 0))
    image = image.resize((OUT_WIDTH, round(image.size[1] * OUT_WIDTH / image.size[0])), Image.LANCZOS)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.out, quality=85)


if __name__ == "__main__":
    main()

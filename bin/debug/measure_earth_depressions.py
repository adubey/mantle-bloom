#!/usr/bin/env python3
"""Earth's closed depressions at node scale: the reference for initial drainage (#218).

Box-averages GMT's `earth_relief_10m` grid (SRTM15+, download
https://oceania.generic-mapping-tools.org/server/earth/earth_relief/earth_relief_10m_g.grd) to
about node spacing, then priority-floods from the world ocean with 8-neighbour connectivity
(periodic in longitude). The world ocean is found at full resolution, so straits narrower than
a box (Gibraltar, the Bosporus) stay open; a box holding any world-ocean subcell floods at sea
level. Reports the area-weighted share of ice-free land (south of 60 S and Greenland excluded)
lying in a closed depression deeper than 1-1,000 m. measure_initial_relief.py's EARTH table
holds the 3x3 result.

Needs h5py, which the backend venv does not carry; use a separate venv:
    python3 -m venv /tmp/earth-venv && /tmp/earth-venv/bin/pip install h5py numpy scipy
    /tmp/earth-venv/bin/python bin/debug/measure_earth_depressions.py earth_relief_10m_g.grd 3

BOX is the box size in 10' cells: 3 is ~55 km, 5 is ~90 km.
"""

import heapq
import json
import sys

import h5py
import numpy as np
from scipy import ndimage

path, box = sys.argv[1], int(sys.argv[2])
f = h5py.File(path)
z = f["z"][:].astype(float)[:-1, :-1]  # drop duplicate edge
lat = f["lat"][:][:-1]
below = z <= 0.0
labels, n = ndimage.label(below, structure=np.ones((3, 3)))
# Wrap in longitude.
for a, b in zip(labels[:, 0], labels[:, -1]):
    if a and b and a != b:
        labels[labels == b] = a
counts = np.bincount(labels.ravel())
counts[0] = 0
ocean_fine = (labels == counts.argmax()).astype(float)
ny, nx = z.shape
ny2, nx2 = ny // box, nx // box
def boxed(a):
    return a[: ny2 * box, : nx2 * box].reshape(ny2, box, nx2, box).mean(axis=(1, 3))
z = boxed(z)
ocean_frac = boxed(ocean_fine)
lat = lat[: ny2 * box].reshape(ny2, box).mean(axis=1)
# A box with any world-ocean subcell is open water at sea level for the flood (straits stay
# open); boxes mostly world ocean are not land.
seed = ocean_frac > 0.0
ocean = ocean_frac >= 0.5
filled = np.full(z.shape, np.inf)
heap = [(min(z[i, j], 0.0), i, j) for i, j in zip(*np.nonzero(seed))]
for lv, i, j in heap:
    filled[i, j] = lv
heapq.heapify(heap)
while heap:
    level, i, j = heapq.heappop(heap)
    if level > filled[i, j]:
        continue
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == dj == 0:
                continue
            a, b = i + di, (j + dj) % nx2
            if not 0 <= a < ny2:
                continue
            lv = max(level, z[a, b])
            if lv < filled[a, b]:
                filled[a, b] = lv
                heapq.heappush(heap, (lv, a, b))

w = np.cos(np.radians(lat))[:, None] * np.ones((1, nx2))
lon = np.linspace(-180, 180, nx2, endpoint=False) + 180 / nx2
LON, LAT = np.meshgrid(lon, lat)
greenland = (LAT > 59) & (LON > -75) & (LON < -10)
land = ~ocean & (z > 0.0) & (LAT > -60) & ~greenland
depth = np.where(land, filled - z, 0.0)
W = w[land].sum()
out = {
    "box_cells": box,
    "land_area_share_of_globe": float(w[land].sum() / w.sum()),
    "depression_share_gt": {str(t): round(float(w[land & (depth > t)].sum() / W), 4) for t in (1, 10, 50, 100, 250, 500, 1000)},
    "closed_below_sea_level_share_of_land": round(float(w[~ocean & ~seed & (z <= 0) & (LAT > -60)].sum() / W), 4),
}
print(json.dumps(out, indent=1))

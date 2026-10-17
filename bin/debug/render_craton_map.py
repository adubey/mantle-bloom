#!/usr/bin/env python3
"""Draw a deliberately plain map of plates and cratons (RFC 0009, #426).

Four things only, on an equirectangular grid: oceanic crust, continental crust, craton, and
plate boundaries (black, where the majority plate changes between neighbouring pixels).
Where plates overlap, the colours are shaded darker. With --mosaic the world is seeded as RFC 0009's mosaic (probe_craton_mosaic's
seeding, formation off), and belts are drawn too. With --cratons-first the world is generated
with RFC 0009's cratons-first layout (#462). Elevation, climate and the app's
renderer are left out on purpose: the question is only where the cratons sit relative to
the plates.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py \
        --seed 7 --out /tmp/seed7.png
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py \
        --seed 7 --years 5e8 --out /tmp/seed7-500myr.png
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py \
        --seed 7 --mosaic --out /tmp/seed7-mosaic.png
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py \
        --seed 7 --cratons-first --land-fraction 0.29 --out /tmp/seed7-cratons-first.png
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py \
        --load world.mbworld --out /tmp/saved.png
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py \
        --seed 7 --years 5e8 --save-world /tmp/seed7-500myr.mbworld --out /tmp/seed7-500myr.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

OCEAN = (168, 197, 224)
CONTINENT = (233, 214, 160)
CRATON = (166, 88, 40)
BELT = (42, 157, 143)
BOUNDARY = (0, 0, 0)
# Window (pixels) for the majority plate, and the share of a window's pixels on a raw plate
# change above which it is overlap rather than one boundary (about 2 / SMOOTH_PX).
SMOOTH_PX = 9
OVERLAP_FLIPS = 0.3


def _sample(world, width: int, height: int, belt_field: str | None):
    """Per pixel: plate index, and a colour class (0 ocean, 1 continent, 2 craton, 3 belt)."""
    from app.elevation_lines import effective_is_continental_from_codes
    from scipy.spatial import cKDTree

    points, plate_idx, kind = [], [], []
    for n, plate in enumerate(p for p in world.plates if p.node_count() > 0):
        pts = plate.all_points_and_elevation()[0]
        continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        k = np.where(continental, 1, 0)
        if belt_field is not None:
            k = np.where(plate.collect(belt_field) < 1e17, 3, k)
        k = np.where(plate.collect("craton_crust_m") > 0.0, 2, k)
        points.append(pts)
        plate_idx.append(np.full(len(pts), n))
        kind.append(k)
    points, plate_idx, kind = np.concatenate(points), np.concatenate(plate_idx), np.concatenate(kind)
    lat = np.radians(90.0 - (np.arange(height) + 0.5) * 180.0 / height)
    lon = np.radians(-180.0 + (np.arange(width) + 0.5) * 360.0 / width)
    lon, lat = np.meshgrid(lon, lat)
    xyz = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1).reshape(-1, 3)
    _, idx = cKDTree(points).query(xyz, workers=-1)
    return plate_idx[idx].reshape(height, width), kind[idx].reshape(height, width)


def _edges(plate: np.ndarray) -> np.ndarray:
    """Pixels whose right or lower neighbour is on another plate."""
    edge = np.zeros_like(plate, dtype=bool)
    edge[:, :-1] |= plate[:, :-1] != plate[:, 1:]
    edge[:, -1] |= plate[:, -1] != plate[:, 0]  # the map wraps in longitude
    edge[:-1, :] |= plate[:-1, :] != plate[1:, :]
    return edge


def render(world, out: Path, width: int, belt_field: str | None, title: str) -> None:
    height = width // 2
    from scipy.ndimage import uniform_filter

    plate, kind = _sample(world, width, height, belt_field)
    # Where plates overlap, the nearest node alternates between them pixel by pixel. Take each
    # pixel's majority plate over a small window so a boundary is drawn once, and shade the
    # overlap -- where the raw plate flips in more than OVERLAP_FLIPS of the window's pixels,
    # far more than a clean boundary's one flip per row -- instead of hatching it.
    share = np.stack([uniform_filter((plate == p).astype(np.float32), SMOOTH_PX, mode="wrap") for p in range(plate.max() + 1)])
    overlap = uniform_filter(_edges(plate).astype(np.float32), SMOOTH_PX, mode="wrap") > OVERLAP_FLIPS
    plate = share.argmax(axis=0)
    rgb = np.array([OCEAN, CONTINENT, CRATON, BELT], dtype=np.uint8)[kind]
    rgb[overlap] = (rgb[overlap] * 0.75).astype(np.uint8)
    rgb[_edges(plate) & ~overlap] = BOUNDARY
    image = Image.fromarray(rgb)
    ImageDraw.Draw(image).text((8, 6), title, fill=BOUNDARY)
    out.parent.mkdir(parents=True, exist_ok=True)
    image.save(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--load", type=Path, help="render a saved .mbworld instead of generating one")
    parser.add_argument("--years", type=float, default=0.0, help="step this long (5 Myr steps) before drawing")
    parser.add_argument("--mosaic", action="store_true", help="seed RFC 0009's mosaic, formation off, and draw belts")
    parser.add_argument("--cratons-first", action="store_true", help="generate with RFC 0009's cratons-first layout")
    parser.add_argument("--land-fraction", type=float, help="the Generate World dialog's land slider (0-1)")
    parser.add_argument("--save-world", type=Path, help="also save the stepped world, to redraw it with --load")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    belt_field = None
    if args.mosaic:
        # Importing the probe registers its belt field before any app module loads.
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import probe_craton_mosaic as probe

        belt_field = probe.BELT_FIELD
    from app import cratons, persistence
    from app import world as world_mod

    if args.mosaic:
        cratons.CRATON_FORMATION_MYR = float("inf")
        cratons.seed_initial_cratons = probe.seed_mosaic
    if args.load is not None:
        world = persistence.load_world_bytes(args.load.read_bytes())
    else:
        world = world_mod.generate_world(args.seed, land_fraction=args.land_fraction, cratons_first=args.cratons_first)
    cratons.ensure_initialized(world)
    while world.elapsed_years < args.years - 1.0:
        world_mod.step_world(world, min(5e6, args.years - world.elapsed_years))
    if args.save_world is not None:
        args.save_world.write_bytes(persistence.save_world_bytes(world))
    layout = "mosaic seeding, formation off" if args.mosaic else "cratons-first layout" if args.cratons_first else "default craton rules"
    source = args.load.name if args.load is not None else f"seed {args.seed}"
    title = f"{source}, {world.elapsed_years / 1e6:.0f} Myr, {layout}. Brown craton, tan continent, blue ocean" + (", teal belt" if args.mosaic else "") + ", black plate boundaries, darker where plates overlap"
    render(world, args.out, args.width, belt_field, title)


if __name__ == "__main__":
    main()

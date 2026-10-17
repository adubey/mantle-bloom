#!/usr/bin/env python3
"""Contact sheets of cratons-first layouts at t = 0, one map per seed (RFC 0009, #470).

Each seed is generated cratons first and drawn with `render_craton_map.render` (brown
craton, tan continental crust, teal belt from the layout, blue ocean, black plate
boundaries), titled with its plate, continent and craton counts. The per-seed PNGs and
sheets of 4 x 5 maps are written to --out.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_layout_sheets.py \
        --out /tmp/sheets
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_layout_sheets.py \
        --seeds 7,11,2024 --land-fraction 0.29 --out /tmp/sheets
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_craton_map as rcm  # noqa: E402

DEFAULT_SEEDS = [7, 11, 2024] + [s for s in range(1, 40) if s not in (7, 11)]
COLUMNS, ROWS = 4, 5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default=",".join(map(str, DEFAULT_SEEDS)))
    parser.add_argument("--land-fraction", type=float, default=0.29)
    parser.add_argument("--width", type=int, default=900)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    from app import cratons, lithosphere_plate
    from app import world as world_mod

    args.out.mkdir(parents=True, exist_ok=True)
    tiles = []
    for seed in map(int, args.seeds.split(",")):
        world = world_mod.generate_world(seed, land_fraction=args.land_fraction, cratons_first=True)
        cratons.ensure_initialized(world)
        layout = lithosphere_plate.cratons_first_layout(seed, land_fraction=args.land_fraction).layout
        area = layout.lattice.area_km2
        share = area[layout.is_core].sum() / area[layout.is_continental].sum()
        n_cratons = len(np.unique(layout.core[layout.is_core]))
        title = f"seed {seed}: {len(world.plates)} plates, {len(layout.continent_area_km2)} continents, {n_cratons} cratons, craton {share:.0%} of continent"

        def belt_at(xyz: np.ndarray, _layout=layout) -> np.ndarray:
            return _layout.belt[_layout.nearest(xyz)]

        path = args.out / f"seed{seed:04d}.png"
        rcm.render(world, path, args.width, None, title, belt_at=belt_at)
        tiles.append(path)
        print(path, flush=True)
    per_sheet = COLUMNS * ROWS
    for start in range(0, len(tiles), per_sheet):
        images = [Image.open(p) for p in tiles[start : start + per_sheet]]
        w, h = images[0].size
        sheet = Image.new("RGB", (COLUMNS * w, ROWS * h + 30), "white")
        for i, image in enumerate(images):
            sheet.paste(image, ((i % COLUMNS) * w, 30 + (i // COLUMNS) * h))
        caption = f"Cratons-first layouts at t = 0, land_fraction {args.land_fraction:g}. Brown craton, tan continental crust, teal belt, blue ocean, black plate boundaries."
        ImageDraw.Draw(sheet).text((10, 8), caption, fill=(0, 0, 0))
        path = args.out / f"sheet{start // per_sheet + 1}.png"
        sheet.save(path)
        print(path, flush=True)


if __name__ == "__main__":
    main()

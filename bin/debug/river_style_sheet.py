#!/usr/bin/env python3
"""Render one saved world under several river styles, as a labelled contact sheet (#933).

The map is rendered once without rivers, then each style draws its rivers on a copy, so a
sheet of many styles costs little more than one render. Each row of the sheet is one crop
of the map and each column is one style.

Usage:
    backend/.venv/bin/python bin/debug/river_style_sheet.py SAVE OUT.png \\
        --style "A=" --style "B=alpha_max=0.5,alpha_gamma=2" \\
        --crop 0.35,0.05,0.53,0.26 [--view combined] [--projection eckert4] [--width 2200]

A style is LABEL=field=value,... over `render_image.RiverStyle`; an empty field list is the
current default. The look before #933 is
"alpha_min=0.15,alpha_max=1,water_min_percentile=0,water_ramp_gamma=0,canyon_shade=0,meander_amplitude=0,smooth_passes=0,width_max_mult=3".
A crop is x0,y0,x1,y1 as fractions of the image; with no crop the whole
map is shown.
"""

from __future__ import annotations

import argparse
import dataclasses
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from PIL import Image, ImageDraw, ImageFont  # noqa: E402

from app import persistence, render_image  # noqa: E402

LABEL_HEIGHT_PX = 34
GAP_PX = 10


def parse_style(spec: str) -> tuple[str, render_image.RiverStyle]:
    label, _, fields = spec.partition("=")
    types = {f.name: f.type for f in dataclasses.fields(render_image.RiverStyle)}
    values = {}
    for item in filter(None, fields.split(",")):
        name, _, raw = item.partition("=")
        if name not in types:
            raise SystemExit(f"unknown RiverStyle field {name!r}; choices: {sorted(types)}")
        kind = types[name]
        values[name] = raw.lower() in ("1", "true", "yes") if kind == "bool" else raw if kind == "str" else int(raw) if kind == "int" else float(raw)
    return label, render_image.RiverStyle(**values)


def render_variants(world, projection: str, view: str, width: int, styles) -> list[Image.Image]:
    """The map under each style, sharing one river-free base render."""
    captured = {}
    draw_rivers = render_image._draw_rivers

    def capture(image, *args, **kwargs):
        captured["call"] = (image.copy(), args, kwargs)  # the caller may add alpha in place
        return image

    render_image._draw_rivers = capture
    try:
        render_image.render_png(world, projection, view, width, round(width * 0.52))
    finally:
        render_image._draw_rivers = draw_rivers
    image, args, kwargs = captured["call"]
    return [draw_rivers(image.copy(), *args, **{**kwargs, "style": style}) for _, style in styles]


def contact_sheet(variants, labels, crops, zoom: float) -> Image.Image:
    font = ImageFont.load_default(size=24)
    rows = []
    for crop in crops:
        tiles = []
        for image in variants:
            w, h = image.size
            box = (round(crop[0] * w), round(crop[1] * h), round(crop[2] * w), round(crop[3] * h))
            tile = image.crop(box)
            tiles.append(tile.resize((round(tile.width * zoom), round(tile.height * zoom)), Image.LANCZOS))
        rows.append(tiles)
    tile_w = max(t.width for row in rows for t in row)
    sheet_w = len(variants) * tile_w + (len(variants) + 1) * GAP_PX
    sheet_h = LABEL_HEIGHT_PX + sum(row[0].height + GAP_PX for row in rows) + GAP_PX
    sheet = Image.new("RGB", (sheet_w, sheet_h), (24, 24, 28))
    draw = ImageDraw.Draw(sheet)
    for col, label in enumerate(labels):
        draw.text((GAP_PX + col * (tile_w + GAP_PX) + 4, 4), label, fill=(240, 240, 240), font=font)
    y = LABEL_HEIGHT_PX
    for row in rows:
        for col, tile in enumerate(row):
            sheet.paste(tile, (GAP_PX + col * (tile_w + GAP_PX), y))
        y += row[0].height + GAP_PX
    return sheet


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("save")
    parser.add_argument("out")
    parser.add_argument("--style", action="append", required=True)
    parser.add_argument("--crop", action="append", default=[])
    parser.add_argument("--view", default="combined")
    parser.add_argument("--projection", default="eckert4")
    parser.add_argument("--width", type=int, default=2200)
    parser.add_argument("--zoom", type=float, default=2.0)
    args = parser.parse_args()

    styles = [parse_style(spec) for spec in args.style]
    crops = [tuple(float(v) for v in crop.split(",")) for crop in args.crop] or [(0.0, 0.0, 1.0, 1.0)]
    with open(args.save, "rb") as f:
        world = persistence.load_world_file(f)
    variants = render_variants(world, args.projection, args.view, args.width, styles)
    sheet = contact_sheet(variants, [label for label, _ in styles], crops, args.zoom)
    buffer = io.BytesIO()
    sheet.save(buffer, format="PNG")
    Path(args.out).write_bytes(buffer.getvalue())
    print(f"wrote {args.out} ({sheet.width}x{sheet.height})")


if __name__ == "__main__":
    main()

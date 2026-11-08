#!/usr/bin/env python3
"""Draw Gate D2 routes with and without the margin limit, side by side (RFC 0010, #737). Read-only.

Replays one seed with the `probe_rift_settings.py` census and, at each listed event's onset, routes
it twice: with no margin limit, as D2's first run did, and with the frozen limit, as the amendment
does. Each PNG shows the plate in an orthographic view, colored by daughter, with craton in gold, the
plate-margin band shaded, and the route in red.

Usage, one process per seed:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/draw_d2_routes.py \
      --seed <n> --events <setting:plate:neighbour:onset_myr> [...] --out-dir <private-dir>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import probe_rift_settings as probe  # noqa: E402

PANEL = 640


def _project(points):
    centre = probe.geometry.normalize(points.mean(axis=0, keepdims=True))[0]
    east = probe.geometry.normalize(np.cross([0.0, 0.0, 1.0], centre))
    north = np.cross(centre, east)
    x, y = points @ east, points @ north
    scale = (PANEL - 60) / max(np.ptp(x), np.ptp(y))
    return 30 + (x - x.min()) * scale, PANEL - 30 - (y - y.min()) * scale, max(2, int(scale * 0.004))


def _panel(plate, state, title, scored, snapshot, cratons_mask, seed):
    points = state["points"]
    px, py, size = _project(points)
    edge = probe._edge_cells(plate)
    band = probe.ROUTE_MARGIN_BAND_KM * 1e3 / probe.lithosphere.PLANET_RADIUS_M
    near = cKDTree(points[edge]).query(points)[0] < 2.0 * np.sin(band / 2.0)
    image = Image.new("RGB", (PANEL, PANEL + 40), "white")
    canvas = ImageDraw.Draw(image)
    mask = snapshot["mask"] if snapshot else np.zeros(len(points), dtype=bool)
    for i in range(len(points)):
        fill = ("#9fbcd8" if mask[i] else "#b7d8b0") if snapshot else "#d8d8d8"
        if near[i]:
            fill = {"#9fbcd8": "#4d6f94", "#b7d8b0": "#5f8a57", "#d8d8d8": "#8a8a8a"}[fill]
        if cratons_mask[i]:
            fill = "#d9a321"
        canvas.rectangle((px[i] - size, py[i] - size, px[i] + size, py[i] + size), fill=fill)
    if snapshot:
        route = snapshot["route"]
        canvas.line([(px[i], py[i]) for i in route], fill="#c0202c", width=4)
        detail = (f"margin share {scored['margin_share']:.2f}, smaller daughter "
                  f"{scored['smaller_daughter_fraction']:.1%}, route {scored['route_km']:.0f} km")
    else:
        detail = f"no route: {scored['status']}"
    canvas.ellipse((px[seed] - 7, py[seed] - 7, px[seed] + 7, py[seed] + 7), outline="black", width=3)
    canvas.text((10, PANEL + 4), title, fill="black")
    canvas.text((10, PANEL + 20), detail, fill="black")
    return image


def draw_event(path, world, plate, state, site, now):
    cratons_mask = plate.collect("craton_crust_m") > 0.0
    belts = plate.collect("belt_formed_years") > 0.0
    limit = probe.ROUTE_MARGIN_MAX_SHARE
    panels = []
    try:
        for name, value in (("first D2 run, no margin limit", 1.0), (f"amendment, margin limit {limit}", limit)):
            probe.ROUTE_MARGIN_MAX_SHARE = value
            scored, snapshot = probe.score_site(world, plate, state, site, cratons_mask, belts)
            panels.append(_panel(plate, state, name, scored, snapshot, cratons_mask, site["seed"]))
    finally:
        probe.ROUTE_MARGIN_MAX_SHARE = limit
    sheet = Image.new("RGB", (2 * PANEL, PANEL + 70), "white")
    ImageDraw.Draw(sheet).text(
        (10, 8), f"{site['setting']}, plate {plate.plate_id}, {now:.0f} Myr: "
        f"{plate.node_areas_m2().sum() / (4 * np.pi * probe.lithosphere.PLANET_RADIUS_M ** 2):.1%} of the sphere. "
        "Darker cells lie within 150 km of the plate edge; gold is craton.", fill="black")
    for column, panel in enumerate(panels):
        sheet.paste(panel, (column * PANEL, 30))
    sheet.save(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--events", nargs="+", required=True, help="setting:plate:neighbour:onset_myr")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    wanted = {(setting, int(plate), int(neighbour), float(onset)) for setting, plate, neighbour, onset in
              (event.split(":") for event in args.events)}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    event_row = probe.Census.event_row

    def drawing_row(self, world, plate, state, site, track, now):
        if (track.setting, track.plate_id, track.neighbour, now) in wanted:
            path = args.out_dir / (f"seed{args.seed}-{track.setting}-plate{track.plate_id}-"
                                   f"neighbour{track.neighbour}-{now:.0f}myr.png")
            draw_event(path, world, plate, state, site, now)
            print(path, flush=True)
        return event_row(self, world, plate, state, site, track, now)

    probe.Census.event_row = drawing_row
    probe.run_seed(args.seed, max(onset for _s, _p, _n, onset in wanted), 1.0)


if __name__ == "__main__":
    main()

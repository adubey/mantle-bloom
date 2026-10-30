#!/usr/bin/env python3
"""Draw a full-plate view of Gate 3b crack paths and optional shell demand.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/plot_gate3b_plate_paths.py \
        --route-json analysis/issue699/interior-belt-path-seed11-density4-stepover150.json \
        --route-json analysis/issue699/interior-belt-path-seed11-density4.json \
        --out /tmp/gate3b-plate-paths.png

The plate is regenerated from the seed and density recorded in the route JSON. The full
overview shows its mesh and craton mask; the two close-ups show the route tips and local
craton context. Opposing arrows show the probe's prescribed opening axes, with their length
enlarged for visibility. This is a visualization of graph-probe inputs and results, not a
production stress field, plate-motion prediction, or simulated split.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree

EARTH_RADIUS_KM = 6_371.0
MOTION_STEP_YEARS = 2.0e6
OVERVIEW_ARROW_LENGTH_KM = 900.0
PATH_COLORS = ((213, 94, 0), (0, 114, 178))
CRATON = (197, 139, 0)


def _unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    length = math.sqrt(sum(value * value for value in vector))
    return tuple(value / length for value in vector)


def _xyz(lat_lon: list[float]) -> tuple[float, float, float]:
    lat, lon = (math.radians(value) for value in lat_lon)
    return math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)


def _lat_lon(vector: tuple[float, float, float]) -> tuple[float, float]:
    return (
        math.degrees(math.asin(max(-1.0, min(1.0, vector[2])))),
        math.degrees(math.atan2(vector[1], vector[0])),
    )


def _mollweide(lat: float, lon: float) -> tuple[float, float]:
    """Return Mollweide map coordinates for longitude centered at zero."""
    phi = math.radians(lat)
    lam = math.radians(lon)
    theta = phi
    if abs(abs(phi) - math.pi / 2) < 1e-10:
        theta = math.copysign(math.pi / 2, phi)
    else:
        for _ in range(12):
            denominator = 2 + 2 * math.cos(2 * theta)
            if abs(denominator) < 1e-12:
                break
            delta = (2 * theta + math.sin(2 * theta) - math.pi * math.sin(phi)) / denominator
            theta -= delta
            if abs(delta) < 1e-12:
                break
    return (
        2 * math.sqrt(2) / math.pi * lam * math.cos(theta),
        math.sqrt(2) * math.sin(theta),
    )


def _opening_axis(
    source_lat_lon: list[float], target_lat_lon: list[float]
) -> tuple[tuple[float, float], list[tuple[float, float]]]:
    """Recreate the prescribed tangent opening axis used by the route probe."""
    source, target = _xyz(source_lat_lon), _xyz(target_lat_lon)
    center = _unit(tuple(source[index] + target[index] for index in range(3)))
    normal = _unit((
        source[1] * target[2] - source[2] * target[1],
        source[2] * target[0] - source[0] * target[2],
        source[0] * target[1] - source[1] * target[0],
    ))
    projection = sum(normal[index] * center[index] for index in range(3))
    drive = _unit(tuple(normal[index] - projection * center[index] for index in range(3)))
    angle = OVERVIEW_ARROW_LENGTH_KM / EARTH_RADIUS_KM
    endpoints = []
    for sign in (-1, 1):
        endpoint = tuple(
            math.cos(angle) * center[index] + sign * math.sin(angle) * drive[index]
            for index in range(3)
        )
        endpoints.append(_lat_lon(endpoint))
    return _lat_lon(center), endpoints


def _font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = (
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold
        else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold
        else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    )
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


def _route_inputs(paths: list[Path]) -> list[dict[str, Any]]:
    if len(paths) != 2:
        raise ValueError("provide exactly two --route-json files (A and B)")
    results = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    shared = (results[0].get("seed"), results[0].get("node_density"), results[0].get("plate_id"))
    for result in results:
        if result.get("status") != "eligible" or not result.get("scenarios"):
            raise ValueError("each route JSON must contain an eligible probe result")
        if (result.get("seed"), result.get("node_density"), result.get("plate_id")) != shared:
            raise ValueError("both route JSONs must use the same seed, density, and plate")
        paths_by_tip = result["scenarios"][0].get("partial_tip_paths_lat_lon_degrees")
        if not paths_by_tip or len(paths_by_tip) != 2 or any(not path for path in paths_by_tip):
            raise ValueError("each route JSON must contain two partial tip paths")
    return results


def _draw_arrow(draw: ImageDraw.ImageDraw, start: tuple[float, float], end: tuple[float, float]) -> None:
    color = (41, 59, 115)
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = math.hypot(dx, dy)
    if length < 2:
        return
    ux, uy = dx / length, dy / length
    draw.line((*start, *end), fill=color, width=5)
    base_x, base_y = end[0] - ux * 11, end[1] - uy * 11
    draw.polygon(
        (end, (base_x - uy * 7, base_y + ux * 7), (base_x + uy * 7, base_y - ux * 7)),
        fill=color,
    )


def _dashed_line(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    color: tuple[int, int, int],
    width: int,
) -> None:
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = math.hypot(dx, dy)
    if not length:
        return
    ux, uy = dx / length, dy / length
    position = 0.0
    while position < length:
        stop = min(position + 9, length)
        draw.line(
            (start[0] + ux * position, start[1] + uy * position,
             start[0] + ux * stop, start[1] + uy * stop),
            fill=color,
            width=width,
        )
        position = stop + 7


def _demand_color(value: float, maximum: float) -> tuple[int, int, int]:
    """Map clipped positive shell demand from cool blue to warm red."""
    t = max(0.0, min(1.0, value / maximum)) if maximum > 0 else 0.0
    stops = ((40, 99, 170), (54, 168, 153), (246, 202, 73), (201, 65, 58))
    position = t * (len(stops) - 1)
    low = min(int(position), len(stops) - 2)
    fraction = position - low
    return tuple(round(stops[low][axis] * (1 - fraction) + stops[low + 1][axis] * fraction)
                 for axis in range(3))


def _prepare_demand_points(
    shell_data: dict[str, Any] | None, plate_lat_lon: list[tuple[float, float]]
) -> tuple[list[tuple[tuple[float, float], float]], float]:
    if shell_data is None:
        return [], 0.0
    shell_map = max(shell_data["maps"], key=lambda item: item["target_representatives"])
    if len(shell_map["mesh_lat_lon_degrees"]) != len(shell_map["mesh_demand_n_per_m"]):
        raise ValueError("shell-response map has mismatched point and demand arrays")
    mesh_demand = np.asarray(shell_map["mesh_demand_n_per_m"], dtype=float)
    demand_max = float(np.quantile(mesh_demand, 0.99))
    mesh_points = shell_map["mesh_lat_lon_degrees"]
    medoid_xyz = np.asarray([_xyz(list(point)) for point in mesh_points])
    plate_xyz = np.asarray([_xyz(list(point)) for point in plate_lat_lon])
    nearest = cKDTree(medoid_xyz).query(plate_xyz, k=1)[1]
    interpolated = mesh_demand[nearest]
    return list(zip(plate_lat_lon, interpolated.tolist())), demand_max


def _draw_shell_legend(draw: ImageDraw.ImageDraw, x: int, y: int, maximum: float) -> None:
    draw.text((x, y), "Shell demand: low", font=_font(14), fill="#39475c")
    for index in range(48):
        x0 = x + 115 + index * 4
        draw.line((x0, y + 8, x0 + 4, y + 8),
                  fill=_demand_color(maximum * index / 47, maximum), width=8)
    draw.text((x + 320, y), "high (p99)", font=_font(14), fill="#39475c")


def _draw_plate_cells(
    draw: ImageDraw.ImageDraw,
    plate_lat_lon: list[tuple[float, float]],
    demand_points: list[tuple[tuple[float, float], float]],
    demand_max: float,
    project,
) -> None:
    if demand_points:
        for point, value in demand_points:
            x, y = project(point)
            draw.ellipse((x - 5, y - 5, x + 5, y + 5), fill=_demand_color(value, demand_max))
        return
    for point in plate_lat_lon:
        x, y = project(point)
        draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill="#dce7ec")


def _draw_overview_legend(
    draw: ImageDraw.ImageDraw, ox0: int, oy1: int, shell_data: dict[str, Any] | None,
    demand_max: float,
) -> None:
    legend_y = oy1 - 83
    draw.ellipse((ox0 + 30, legend_y + 4, ox0 + 42, legend_y + 16), fill=CRATON)
    draw.text((ox0 + 50, legend_y), "cratons", font=_font(16), fill="#39475c")
    draw.line((ox0 + 144, legend_y + 10, ox0 + 182, legend_y + 10), fill=PATH_COLORS[0], width=5)
    draw.text((ox0 + 191, legend_y), "crack paths", font=_font(16), fill="#39475c")
    draw.line((ox0 + 311, legend_y + 10, ox0 + 346, legend_y + 10), fill="#293b73", width=5)
    draw.polygon(((ox0 + 352, legend_y + 10), (ox0 + 340, legend_y + 3),
                  (ox0 + 340, legend_y + 17)), fill="#293b73")
    draw.text((ox0 + 362, legend_y), "opening arrows", font=_font(16), fill="#39475c")
    draw.text((ox0 + 560, legend_y), "A / B link to close-ups →", font=_font(16, True), fill="#39475c")
    if shell_data is not None:
        _draw_shell_legend(draw, ox0 + 30, legend_y + 25, demand_max)


def _draw_local_demand(draw, demand_points, demand_max, local, local_xy, bounds) -> None:
    if not demand_points:
        return
    min_x, max_x, min_y, max_y = bounds
    for point, value in demand_points:
        x, y = local(point)
        if min_x <= x <= max_x and min_y <= y <= max_y:
            px, py = local_xy(point)
            draw.ellipse((px - 4, py - 4, px + 4, py + 4), fill=_demand_color(value, demand_max))


def render(world: Any, results: list[dict[str, Any]], out: Path, a_context_scale: float,
           shell_data: dict[str, Any] | None = None) -> None:
    from app import cratons

    cratons.ensure_initialized(world)
    plate_id = results[0]["plate_id"]
    plate = next((item for item in world.plates if item.plate_id == plate_id), None)
    if plate is None:
        raise ValueError(f"regenerated world has no plate {plate_id}")
    nodes = plate.surface_nodes()
    points = nodes.world_xyz
    craton_mask = plate.collect("craton_crust_m") > 0
    plate_lat_lon = [
        (math.degrees(math.asin(max(-1.0, min(1.0, point[2] / math.sqrt(sum(v * v for v in point)))))),
         math.degrees(math.atan2(point[1], point[0])))
        for point in points
    ]
    craton_points = [point for point, is_craton in zip(plate_lat_lon, craton_mask) if is_craton]
    demand_points, demand_max = _prepare_demand_points(shell_data, plate_lat_lon)
    route_data = [result["scenarios"][0] for result in results]
    tracks_by_case = [item["partial_tip_paths_lat_lon_degrees"] for item in route_data]

    width, height = 1900, 1270
    image = Image.new("RGB", (width, height), "#f7f8fa")
    draw = ImageDraw.Draw(image)
    title_suffix = " and shell response" if shell_data is not None else ""
    draw.text((65, 35), f"Seed {results[0]['seed']} · plate, cratons, crack paths{title_suffix}",
              font=_font(38, True), fill="#172033")
    subtitle = ("Color = nearest-medoid shell demand at generation (clipped at p99); paths use the 2 Myr step. "
                "Field is provisional; gold marks cratons." if shell_data is not None
                else "Full plate view with route locations marked; close-ups show both gaps and nearby craton context.")
    draw.text((67, 91), subtitle,
              font=_font(23), fill="#45536a")

    # Plate overview.
    ox0, oy0, ox1, oy1 = 60, 155, 1130, 1060
    draw.rounded_rectangle((ox0, oy0, ox1, oy1), radius=22, fill="white", outline="#d7dde6", width=2)
    draw.text((ox0 + 28, oy0 + 20), "Plate overview", font=_font(28, True), fill="#18243a")
    px0, py0, px1, py1 = ox0 + 48, oy0 + 86, ox1 - 40, oy1 - 112
    scale = min((px1 - px0) / (4 * math.sqrt(2)), (py1 - py0) / (2 * math.sqrt(2)))
    center_x, center_y = (px0 + px1) / 2, (py0 + py1) / 2

    def global_xy(point: tuple[float, float]) -> tuple[float, float]:
        x, y = _mollweide(*point)
        return center_x + x * scale, center_y - y * scale

    for latitude in (-60, -30, 0, 30, 60):
        draw.line([global_xy((latitude, longitude)) for longitude in range(-180, 181, 2)],
                  fill="#e6ebf0", width=1)
    for longitude in range(-150, 180, 30):
        draw.line([global_xy((latitude, longitude)) for latitude in range(-89, 90, 2)],
                  fill="#e6ebf0", width=1)
    radius_x, radius_y = 2 * math.sqrt(2) * scale, math.sqrt(2) * scale
    draw.ellipse((center_x - radius_x, center_y - radius_y,
                  center_x + radius_x, center_y + radius_y), outline="#9ba8b8", width=2)
    _draw_plate_cells(draw, plate_lat_lon, demand_points, demand_max, global_xy)
    for point in craton_points:
        x, y = global_xy(point)
        draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=CRATON)

    centers_and_ends = []
    overview_centers = []
    for tracks in tracks_by_case:
        center, endpoints = _opening_axis(tracks[0][0], tracks[1][0])
        centers_and_ends.append((center, endpoints))
        overview_centers.append(global_xy(center))
        for index, track in enumerate(tracks):
            projected = [global_xy(tuple(point)) for point in track]
            draw.line(projected, fill=PATH_COLORS[index], width=5, joint="curve")
            x, y = projected[0]
            draw.ellipse((x - 6, y - 6, x + 6, y + 6), fill=PATH_COLORS[index], outline="white", width=2)
            x, y = projected[-1]
            draw.rectangle((x - 6, y - 6, x + 6, y + 6), fill=PATH_COLORS[index], outline="white", width=2)
    for center, endpoints in centers_and_ends:
        projected_center = global_xy(center)
        for endpoint in endpoints:
            _draw_arrow(draw, projected_center, global_xy(endpoint))
    for index, (x, y) in enumerate(overview_centers):
        radius = 18
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), outline="#172033", width=3)
        draw.ellipse((x - radius + 4, y - radius + 4,
                      x + radius - 4, y + radius - 4), outline="white", width=2)
        draw.text((x + 19, y - 18), chr(65 + index), font=_font(21, True),
                  fill="#172033", stroke_width=2, stroke_fill="white")

    _draw_overview_legend(draw, ox0, oy1, shell_data, demand_max)

    # Two local windows; widen A to show the surrounding craton pattern.
    for index, (result, scenario, tracks) in enumerate(zip(results, route_data, tracks_by_case)):
        x0, y0, x1, y1 = 1160, 155 + index * 455, 1840, 590 + index * 455
        all_points = [point for track in tracks for point in track]
        lat0 = sum(point[0] for point in all_points) / len(all_points)
        lon0 = sum(point[1] for point in all_points) / len(all_points)

        def local(point: tuple[float, float]) -> tuple[float, float]:
            return ((point[1] - lon0) * 111.32 * math.cos(math.radians(lat0)),
                    (point[0] - lat0) * 111.32)

        local_points = [local(tuple(point)) for point in all_points]
        min_x, max_x = min(x for x, _ in local_points), max(x for x, _ in local_points)
        min_y, max_y = min(y for _, y in local_points), max(y for _, y in local_points)
        pad_x = max(50, (max_x - min_x) * 0.18)
        pad_y = max(55, (max_y - min_y) * 0.22)
        min_x, max_x = min_x - pad_x, max_x + pad_x
        min_y, max_y = min_y - pad_y, max_y + pad_y
        zoom = a_context_scale if index == 0 else 1.0
        middle_x, middle_y = (min_x + max_x) / 2, (min_y + max_y) / 2
        half_x, half_y = (max_x - min_x) * zoom / 2, (max_y - min_y) * zoom / 2
        min_x, max_x = middle_x - half_x, middle_x + half_x
        min_y, max_y = middle_y - half_y, middle_y + half_y

        draw.rounded_rectangle((x0, y0, x1, y1), radius=22, fill="white", outline="#d7dde6", width=2)
        draw.text((x0 + 24, y0 + 18), f"{chr(65 + index)} · {result['endpoint_gap_km']:.0f} km step-over",
                  font=_font(25, True), fill="#18243a")
        draw.text((x0 + 24, y0 + 53),
                  f"Start {result['endpoint_gap_km']:.1f} km → after arrest {scenario['remaining_tip_gap_km']:.1f} km",
                  font=_font(19), fill="#39475c")
        plot_x0, plot_y0, plot_x1, plot_y1 = x0 + 55, y0 + 95, x1 - 25, y1 - 62

        def local_xy(point: tuple[float, float]) -> tuple[float, float]:
            x, y = local(point)
            return (plot_x0 + (x - min_x) / (max_x - min_x) * (plot_x1 - plot_x0),
                    plot_y1 - (y - min_y) / (max_y - min_y) * (plot_y1 - plot_y0))

        for fraction in (0, 0.25, 0.5, 0.75, 1):
            gx = plot_x0 + fraction * (plot_x1 - plot_x0)
            gy = plot_y0 + fraction * (plot_y1 - plot_y0)
            draw.line((gx, plot_y0, gx, plot_y1), fill="#edf0f4", width=1)
            draw.line((plot_x0, gy, plot_x1, gy), fill="#edf0f4", width=1)
        _draw_local_demand(draw, demand_points, demand_max, local, local_xy,
                           (min_x, max_x, min_y, max_y))
        visible_cratons = 0
        for point in craton_points:
            x, y = local(point)
            if min_x <= x <= max_x and min_y <= y <= max_y:
                px, py = local_xy(point)
                draw.ellipse((px - 3, py - 3, px + 3, py + 3), fill=CRATON, outline=(147, 102, 0), width=1)
                visible_cratons += 1

        starts, ends = [], []
        for path_index, track in enumerate(tracks):
            projected = [local_xy(tuple(point)) for point in track]
            starts.append(projected[0])
            ends.append(projected[-1])
            draw.line(projected, fill=PATH_COLORS[path_index], width=6, joint="curve")
            x, y = projected[0]
            draw.ellipse((x - 8, y - 8, x + 8, y + 8), fill=PATH_COLORS[path_index], outline="white", width=2)
            x, y = projected[-1]
            draw.rectangle((x - 8, y - 8, x + 8, y + 8), fill=PATH_COLORS[path_index], outline="white", width=2)
        _dashed_line(draw, starts[0], starts[1], (100, 116, 139), 3)
        _dashed_line(draw, ends[0], ends[1], (39, 54, 75), 4)
        draw.text((x0 + 22, y1 - 34),
                  f"Gold = cratons ({visible_cratons} cells shown) · circles = start · squares = arrested",
                  font=_font(15), fill="#46546a")

    footer_y = 1085
    draw.rounded_rectangle((60, footer_y, 1840, 1235), radius=18,
                           fill="#eef2f7", outline="#d7dde6", width=1)
    draw.text((86, footer_y + 20),
              f"Gold = actual craton mask; A detail is widened by {a_context_scale:g}× for surrounding context.",
              font=_font(20, True), fill="#27364b")
    draw.text((86, footer_y + 58),
              "Opposing arrows show the prescribed opening axis at A and B; displayed length is enlarged for visibility.",
              font=_font(19), fill="#45536a")
    draw.text((86, footer_y + 96),
              ("Provisional static shell response: solver residual missed the preregistered 1% cutoff; "
               "no damage or split applied." if shell_data is not None else
               "Illustrative graph probe only: no production stress field, measured plate motion, or split applied."),
              font=_font(18), fill="#45536a")
    out.parent.mkdir(parents=True, exist_ok=True)
    image.save(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--route-json", type=Path, action="append", required=True,
                        help="one route JSON from probe_gate3b_interior_belt_paths.py; repeat twice, A then B")
    parser.add_argument("--shell-json", type=Path,
                        help="optional coarse shell-response output from probe_gate3b_local_load.py")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--a-context-scale", type=float, default=3.2,
                        help="linear context expansion for close-up A (default: 3.2)")
    args = parser.parse_args()
    if args.a_context_scale < 1:
        parser.error("--a-context-scale must be at least 1")
    results = _route_inputs(args.route_json)
    shell_data = json.loads(args.shell_json.read_text(encoding="utf-8")) if args.shell_json else None
    if shell_data is not None:
        route = results[0]
        if (shell_data.get("seed"), shell_data.get("node_density"), shell_data.get("plate_id")) != (
                route.get("seed"), route.get("node_density"), route.get("plate_id")):
            parser.error("--shell-json must match the route seed, density, and plate")
    from app import world as world_mod

    world = world_mod.generate_world(
        seed=int(results[0]["seed"]),
        node_density=float(results[0]["node_density"]),
        cratons_first=True,
        synthetic_relief=False,
    )
    world.simulate_climate_biomes = False
    world_mod.step_world(world, MOTION_STEP_YEARS)
    render(world, results, args.out, args.a_context_scale, shell_data)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

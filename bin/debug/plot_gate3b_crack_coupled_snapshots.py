#!/usr/bin/env python3
"""Plot the first p10 damage snapshots before and after failed links disconnect a plate.

Usage from the repository root:
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
        backend/.venv/bin/python bin/debug/plot_gate3b_crack_coupled_snapshots.py \
        --out analysis/issue699/crack-coupled-components.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from app import cratons, elevation_lines, world as world_mod
import gate3b_diagnostic_common as shared
import probe_gate3b_crack_coupled as diagnostic
import probe_gate3b_slab_pull as slab_screen
import replay_gate2_capacity as capacity_audit


def _load_fields():
    world = world_mod.generate_world(
        diagnostic.SEED, node_density=diagnostic.DENSITY,
        cratons_first=True, synthetic_relief=False,
    )
    world.simulate_climate_biomes = False
    cratons.ensure_initialized(world)
    world_mod.step_world(world, diagnostic.STARTUP_YEARS)
    plate, _candidates, candidate, _native_p95, _centers, _scale, field, midline_force, _inputs, subducting = (
        diagnostic._source_and_target(world)
    )
    source_base = np.asarray(field.total_n, dtype=float)
    demands = {}
    directions = {}
    for slab_factor in (1.0, 2.0, 4.0):
        force = source_base + (slab_factor - 1.0) * field.slab_pull_n
        _shell, demands[slab_factor], directions[slab_factor] = diagnostic._opening_capacity_demand(
            plate, force
        )
    return world, plate, candidate, field, midline_force, subducting, demands, directions


def _failed_link_components(plate, cracked, direction, stress, yield_strength):
    broken, count, labels, _sizes, _pairs = shared.local_failure_components(
        plate, cracked, direction, stress, yield_strength,
        diagnostic.RESIDUAL_STRENGTH_FRACTION,
    )
    return broken, count, labels


def _snapshots(plate, field, demands, directions):
    n = plate.node_count()
    thickness = np.maximum(
        plate.collect("crustal_thickness_m") + plate.collect("mantle_lithosphere_thickness_m"),
        1.0,
    )
    # Capacity is assigned by crust type, as in the corrected diagnostic run.
    continental = elevation_lines.effective_is_continental_from_codes(
        plate.collect("crust_type_code"), plate.crust_type == "continental"
    )
    raw = capacity_audit._physical_capacities()
    capacity = np.where(
        continental,
        float(raw["continental"][0]),
        float(raw["oceanic_mature"][0]),
    )
    yield_strength = capacity / thickness
    areas = plate.node_areas_m2()
    epsilon = np.zeros(n, dtype=float)
    seconds = diagnostic.STEP_YEARS * 365.25 * 24.0 * 3600.0
    snapshots = {}
    factors = {}
    for elapsed in range(2, 22, 2):
        slab_factor = 1.0 if elapsed < 10 else (2.0 if elapsed < 20 else 4.0)
        factors[elapsed] = slab_factor
        stress = demands[slab_factor] / thickness
        epsilon = diagnostic._advance_damage(
            epsilon, stress, yield_strength, seconds
        )
        if elapsed in (18, 20):
            cracked = epsilon >= diagnostic.EPSILON_FAIL
            broken, count, labels = _failed_link_components(
                plate, cracked, directions[slab_factor], stress, yield_strength
            )
            snapshots[elapsed] = {
                "cracked": cracked.copy(),
                "broken_pairs": broken,
                "component_count": count,
                "component_labels": labels,
                "component_cell_counts": np.bincount(labels, minlength=count),
                "cross_component_failed_link_count": int(
                    np.sum(labels[broken[:, 0]] != labels[broken[:, 1]])
                ),
                "failed_cell_count": int(cracked.sum()),
                "failed_cell_area_fraction": float(areas[cracked].sum() / areas.sum()),
                "slab_factor": slab_factor,
                "stress": stress,
            }
    return snapshots


def _map_projection(points):
    return shared.map_spherical_coordinates(points)


def _slab_vectors(field, xyz, candidate, subducting):
    forces = np.asarray(field.slab_pull_n, dtype=float)
    if forces.shape[0] != len(xyz):
        raise RuntimeError("slab-pull force field does not map one-to-one to target plate cells")
    magnitude = np.linalg.norm(forces, axis=1)
    coordinate = np.asarray(candidate["coord"], dtype=float)
    lo, hi = np.quantile(coordinate, [slab_screen.OUTER_FRACTION, 1.0 - slab_screen.OUTER_FRACTION])
    side_masks = (coordinate <= lo, coordinate >= hi)
    sample_groups = []
    for side, side_mask in enumerate(side_masks):
        valid = np.flatnonzero(np.asarray(subducting, dtype=bool) & side_mask & (magnitude > 0.0))
        if not len(valid):
            raise RuntimeError(f"bilateral slab-pull side {side + 1} has no force-bearing cells")
        # Sample each flank independently so the stronger side cannot hide the weaker one.
        sample_groups.extend((index, side) for index in valid[np.argsort(magnitude[valid])[-min(8, len(valid)):]])
    sample = np.asarray([item[0] for item in sample_groups], dtype=int)
    side_ids = np.asarray([item[1] for item in sample_groups], dtype=int)
    lat = np.arcsin(np.clip(xyz[sample, 2], -1.0, 1.0))
    lon = np.arctan2(xyz[sample, 1], xyz[sample, 0])
    east = np.column_stack((-np.sin(lon), np.cos(lon), np.zeros(len(sample))))
    north = np.column_stack((
        -np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)
    ))
    fx = np.sum(forces[sample] * east, axis=1)
    fy = np.sum(forces[sample] * north, axis=1)
    norm = np.hypot(fx, fy)
    # A fixed map-space arrow length emphasizes direction; arrow lengths do not encode force.
    arrow_length = 2.7
    u = np.divide(fx, norm, out=np.zeros_like(fx), where=norm > 0.0) * arrow_length
    u /= np.maximum(np.cos(lat), 0.25)
    v = np.divide(fy, norm, out=np.zeros_like(fy), where=norm > 0.0) * arrow_length
    return sample, u, v, side_ids


def plot(out_path: Path):
    _world, plate, candidate, field, _midline, subducting, demands, directions = _load_fields()
    snapshots = _snapshots(plate, field, demands, directions)
    points, _elevation = plate.all_points_and_elevation()
    xyz, lon, lat = _map_projection(points)
    craton_mask = plate.collect("craton_crust_m") > 0.0
    slab_indices, slab_u, slab_v, slab_side_ids = _slab_vectors(
        field, xyz, candidate, subducting
    )
    width, height = 2400, 1220
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default()
    title_font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 32)
    subtitle_font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 23)
    draw.text((width // 2, 26), "Gate 3b crack-coupled diagnostic: connected plate to first disconnection", fill="#18212b", font=title_font, anchor="mt")
    panel_top, panel_bottom, panel_width = 145, 1010, 1080
    lefts = (85, 1235)
    extent = max(float(np.ptp(lon)), float(np.ptp(lat)))
    pad = max(extent * 0.08, 2.0)
    x_min, x_max = float(lon.min() - pad), float(lon.max() + pad)
    y_min, y_max = float(lat.min() - pad), float(lat.max() + pad)
    component_colors = ["#91b4d4", "#70a98b", "#ab91c6", "#e2a86f", "#74b9bf", "#d98187", "#b4b466", "#7390bd", "#d18fb5", "#9d9d9d"]

    def map_point(x_value, y_value, left):
        x = left + 55 + (x_value - x_min) / max(x_max - x_min, 1e-9) * (panel_width - 285)
        y = panel_bottom - 35 - (y_value - y_min) / max(y_max - y_min, 1e-9) * (panel_bottom - panel_top - 75)
        return int(x), int(y)

    for left, elapsed in zip(lefts, (18, 20), strict=True):
        snap = snapshots[elapsed]
        right = left + panel_width
        draw.rounded_rectangle((left, panel_top, right, panel_bottom), radius=14, outline="#d5dce3", width=2)
        heading = f"{elapsed} Myr  |  {snap['component_count']} connected component{'s' if snap['component_count'] != 1 else ''}"
        draw.text((left + panel_width // 2, 95), heading, fill="#243447", font=subtitle_font, anchor="mm")
        stat = f"{snap['failed_cell_count']:,} locally failed cells · {snap['failed_cell_area_fraction']:.1%} of plate area · {len(snap['broken_pairs']):,} failed links"
        draw.text((left + panel_width // 2, 124), stat, fill="#586878", font=font, anchor="mm")
        for x_value in np.linspace(x_min, x_max, 5):
            x0, y0 = map_point(x_value, y_min, left)
            x1, y1 = map_point(x_value, y_max, left)
            draw.line((x0, y0, x1, y1), fill="#edf0f3", width=1)
        for y_value in np.linspace(y_min, y_max, 5):
            x0, y0 = map_point(x_min, y_value, left)
            x1, y1 = map_point(x_max, y_value, left)
            draw.line((x0, y0, x1, y1), fill="#edf0f3", width=1)
        labels = snap["component_labels"]
        for index, (x_value, y_value) in enumerate(zip(lon, lat, strict=True)):
            x, y = map_point(x_value, y_value, left)
            color = component_colors[int(labels[index]) % len(component_colors)]
            draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=color)
        for index in np.flatnonzero(craton_mask):
            x, y = map_point(lon[index], lat[index], left)
            draw.ellipse((x - 7, y - 7, x + 7, y + 7), fill="#f2c14e", outline="#372f20", width=2)
        for index in np.flatnonzero(snap["cracked"]):
            x, y = map_point(lon[index], lat[index], left)
            draw.line((x - 5, y - 5, x + 5, y + 5), fill="#b51f2e", width=2)
            draw.line((x - 5, y + 5, x + 5, y - 5), fill="#b51f2e", width=2)
        component_sizes = snap["component_cell_counts"]
        component_order = np.argsort(component_sizes)[::-1]
        component_number = np.empty(len(component_order), dtype=int)
        component_number[component_order] = np.arange(1, len(component_order) + 1)
        largest_label = int(component_order[0])
        if elapsed == 20:
            major_x, major_y = map_point(float(np.mean(lon[labels == largest_label])), float(np.mean(lat[labels == largest_label])), left)
            draw.rounded_rectangle((major_x - 18, major_y - 15, major_x + 18, major_y + 15), radius=8, fill="white", outline="#365a78", width=2)
            draw.text((major_x, major_y), str(int(component_number[largest_label])), fill="#24445f", font=subtitle_font, anchor="mm")
            callout_x = right - 210
            for rank, label_id in enumerate(component_order):
                number = rank + 1
                row_y = panel_top + 100 + rank * 49
                if label_id != largest_label:
                    cells = np.flatnonzero(labels == label_id)
                    cx = float(np.mean(lon[cells]))
                    cy = float(np.mean(lat[cells]))
                    marker_x, marker_y = map_point(cx, cy, left)
                    draw.ellipse((marker_x - 12, marker_y - 12, marker_x + 12, marker_y + 12), fill="white", outline="#501a74", width=3)
                    draw.text((marker_x, marker_y), str(number), fill="#501a74", font=font, anchor="mm")
                    draw.line((marker_x + 12, marker_y, callout_x - 8, row_y), fill="#501a74", width=2)
                draw.rounded_rectangle((callout_x, row_y - 15, right - 18, row_y + 15), radius=7, fill="white", outline="#b8a3c7", width=2)
                suffix = "cell" if component_sizes[label_id] == 1 else "cells"
                draw.text((callout_x + 10, row_y - 8), f"{number}: {component_sizes[label_id]:,} {suffix}", fill="#354455", font=font)
        for index, u, v, side in zip(slab_indices, slab_u, slab_v, slab_side_ids, strict=True):
            x0, y0 = map_point(lon[index], lat[index], left)
            x1, y1 = map_point(lon[index] + u, lat[index] + v, left)
            arrow_color = ("#152b43", "#00796b")[side]
            draw.line((x0, y0, x1, y1), fill=arrow_color, width=5)
            angle = np.arctan2(y1 - y0, x1 - x0)
            head = 16
            left_head = (x1 - head * np.cos(angle - 0.48), y1 - head * np.sin(angle - 0.48))
            right_head = (x1 - head * np.cos(angle + 0.48), y1 - head * np.sin(angle + 0.48))
            draw.polygon(((x1, y1), left_head, right_head), fill=arrow_color)
        for a, b in snap["broken_pairs"]:
            if labels[a] != labels[b]:
                p0, p1 = map_point(lon[a], lat[a], left), map_point(lon[b], lat[b], left)
                draw.line((*p0, *p1), fill="#4a1a09", width=12)
                draw.line((*p0, *p1), fill="#ff6d00", width=8)
                mid_x, mid_y = (p0[0] + p1[0]) // 2, (p0[1] + p1[1]) // 2
                draw.polygon(((mid_x, mid_y - 7), (mid_x + 7, mid_y), (mid_x, mid_y + 7), (mid_x - 7, mid_y)), fill="#ff6d00", outline="#4a1a09")
        draw.text((left + panel_width // 2, 1032), "Longitude (centered on plate)", fill="#586878", font=font, anchor="mm")
        draw.text((left + 20, panel_top + 15), "N", fill="#586878", font=font)
    legend_y = 1092
    legend_items = [("#f2c14e", "craton", "circle"), ("#b51f2e", "locally failed cell", "cross"), ("#ff6d00", "link cut between components", "line"), ("#501a74", "fragment and callout", "ring"), ("#152b43", "slab-pull side A", "arrow"), ("#00796b", "slab-pull side B", "arrow"), ("#91b4d4", "connected-component region", "dot")]
    for (color, label, kind), x in zip(legend_items, (80, 365, 690, 1065, 1400, 1690, 1985), strict=True):
        if kind == "circle":
            draw.ellipse((x, legend_y - 7, x + 15, legend_y + 8), fill=color, outline="#372f20")
        elif kind == "cross":
            draw.line((x, legend_y - 7, x + 15, legend_y + 8), fill=color, width=3)
            draw.line((x, legend_y + 8, x + 15, legend_y - 7), fill=color, width=3)
        elif kind == "line":
            draw.line((x, legend_y, x + 18, legend_y), fill=color, width=4)
        elif kind == "ring":
            draw.ellipse((x, legend_y - 8, x + 16, legend_y + 8), outline=color, width=3)
        elif kind == "arrow":
            draw.line((x, legend_y, x + 18, legend_y), fill=color, width=4)
            draw.polygon(((x + 18, legend_y), (x + 9, legend_y - 6), (x + 9, legend_y + 6)), fill=color)
        else:
            draw.ellipse((x, legend_y - 5, x + 11, legend_y + 6), fill=color)
        draw.text((x + 25, legend_y - 8), label, fill="#354455", font=font)
    note = "Arrows are sampled separately from both slab-pull flanks; equal lengths show direction, not relative force. Orange strokes are severed links. Five detached fragments contain only 1 or 2 cells."
    draw.text((width // 2, 1160), note, fill="#586878", font=font, anchor="mm")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(out_path)
    return snapshots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    snapshots = plot(args.out)
    print({
        "output": str(args.out),
        "snapshots": {
            str(time): {
                "components": row["component_count"],
                "component_cell_counts": row["component_cell_counts"].tolist(),
                "cross_component_failed_links": row["cross_component_failed_link_count"],
                "failed_cells": row["failed_cell_count"],
                "failed_links": int(len(row["broken_pairs"])),
            }
            for time, row in snapshots.items()
        },
    })


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Map the modeled drainage catchment and wet lake components through a replay.

Each requested phase selects every governing lake body whose wet members intersect a
component containing one of the fixed geographic anchors. A body's `members` are the
catchment cells assigned to that lake hierarchy body. The figure colors dry land members
orange and the selected connected wet components blue. This is a map of the model's lake
hierarchy catchments, not a claim that watershed divides are independently observed.

Run from the repository root. For an automatically selected deep lake:

    backend/.venv/bin/python bin/debug/trace_lake_drainage_phases.py <save> \\
        --anchor-step 4 --reference-volume-km3 136570 --phase before=4 --phase after=5 \\
        --myr 0.5 --step-years 100000 --out /tmp/deep-catchments.npz \\
        --figure-out /tmp/deep-catchments.png

For fixed geographic anchors, provide `--anchors` and the requested phase steps.
"""

from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree

from replay_lib import load_replay, parse_replay_args


def add_arguments(ap) -> None:
    ap.add_argument("--anchors", type=Path, help="JSON with fixed XYZ anchor_points")
    ap.add_argument("--anchor-step", type=int, help="step to select the reference component")
    ap.add_argument("--reference-volume-km3", type=float, help="target volume for automatic anchor selection")
    ap.add_argument("--anchor-count", type=int, default=8)
    ap.add_argument("--phase", action="append", required=True, help="phase label=step; repeat for each phase")
    ap.add_argument("--figure-out", type=Path, help="optional PNG figure output")


args = parse_replay_args(__doc__, add_arguments)
if args.anchors is None:
    if args.anchor_step is None or args.reference_volume_km3 is None:
        raise ValueError("provide --anchors or both --anchor-step and --reference-volume-km3")
elif args.anchor_step is not None or args.reference_volume_km3 is not None:
    raise ValueError("fixed anchors cannot be combined with automatic component selection")
if args.anchor_count < 1:
    raise ValueError("anchor count must be positive")
if args.step_years != 100_000:
    raise ValueError("RFC 0014 held-out flicker phases use the registered 100 kyr schedule")


def parse_phases(values: list[str]) -> list[tuple[str, int]]:
    phases = []
    for value in values:
        label, separator, step_text = value.partition("=")
        if not separator or not label.strip():
            raise ValueError(f"invalid phase {value!r}; use label=step")
        step = int(step_text)
        if step < 1:
            raise ValueError("phase steps must be positive")
        phases.append((label.strip(), step))
    if len({label for label, _ in phases}) != len(phases):
        raise ValueError("phase labels must be unique")
    if len({step for _, step in phases}) != len(phases):
        raise ValueError("phase steps must be unique")
    return sorted(phases, key=lambda item: item[1])


phases = parse_phases(args.phase)
last_step = max(step for _, step in phases)
required_anchor_step = args.anchor_step or 0
if required_anchor_step > last_step:
    raise ValueError("anchor step must not be later than the final requested phase")
os.environ["MANTLE_BLOOM_LAKE_LINEAGE_TRACE"] = "1"
world, header = load_replay(args, issue=675)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

world.ice_age_period_years = 0.0
visible_depth = hydrology.LAKE_MIN_VISIBLE_DEPTH_M
anchor_points = (
    np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)
    if args.anchors else None
)
if anchor_points is not None and (anchor_points.ndim != 2 or anchor_points.shape[1] != 3 or len(anchor_points) == 0):
    raise ValueError("anchor_points must be a nonempty array of XYZ unit-sphere coordinates")
phase_records: dict[str, dict] = {}
phase_arrays: dict[str, dict[str, np.ndarray]] = {}
anchor_selection = None


def component_sets(fields):
    wet = np.asarray(fields.lake_depth) > visible_depth
    groups = [np.asarray(group, dtype=np.int64) for group in hydrology.lake_components(wet, fields.neighbor_idx)]
    labels = np.full(len(wet), -1, dtype=np.int32)
    for component_id, nodes in enumerate(groups):
        labels[nodes] = component_id
    return wet, groups, labels


def choose_automatic_anchors(fields, groups):
    hyps = fields.hypsometry
    selected = min(
        groups,
        key=lambda nodes: abs(float(np.dot(fields.lake_depth[nodes], hyps.area_m2[nodes]) / 1e9)
                              - args.reference_volume_km3),
    )
    points = np.asarray(fields.points)[selected]
    center = points.mean(axis=0)
    chosen_local = [int(np.argmax(np.linalg.norm(points - center, axis=1)))]
    min_distance = np.linalg.norm(points - points[chosen_local[0]], axis=1)
    while len(chosen_local) < min(args.anchor_count, len(selected)):
        next_idx = int(np.argmax(min_distance))
        chosen_local.append(next_idx)
        min_distance = np.minimum(min_distance, np.linalg.norm(points - points[next_idx], axis=1))
    return points[np.asarray(chosen_local, dtype=int)].copy(), selected


def capture_phase(label: str, step: int, fields) -> None:
    wet, groups, labels = component_sets(fields)
    distance, nearest = cKDTree(fields.points).query(anchor_points, k=1)
    anchor_labels = labels[nearest]
    selected_ids = sorted({int(value) for value in anchor_labels if value >= 0})
    selected_wet = np.zeros(len(wet), dtype=bool)
    components = []
    for component_id in selected_ids:
        nodes = groups[component_id]
        selected_wet[nodes] = True
        area_m2 = fields.hypsometry.area_m2[nodes]
        area_fraction = fields.hypsometry.wet_fraction(
            fields.hypsometry.level_from_depth(fields.lake_depth[nodes], nodes), nodes
        )
        components.append({
            "component_id": component_id,
            "anchor_count": int(np.count_nonzero(anchor_labels == component_id)),
            "volume_km3": float(np.dot(fields.lake_depth[nodes], area_m2) / 1e9),
            "cell_footprint_area_km2": float(area_m2.sum() / 1e6),
            "hypsometric_wet_area_km2": float(np.dot(area_fraction, area_m2) / 1e6),
            "node_count": int(len(nodes)),
        })

    selected_body_ids = []
    drainage = np.zeros(len(wet), dtype=bool)
    for body in lakes.governing_bodies(
        fields.lake_forest, fields.elevation, fields.lake_depth, fields.hypsometry
    ):
        body_wet = body.members[wet[body.members]]
        if len(body_wet) and np.any(selected_wet[body_wet]):
            drainage[body.members] = True
            selected_body_ids.append(int(body.lake_id))

    area_m2 = fields.hypsometry.area_m2
    ocean = np.asarray(fields.is_ocean, dtype=bool)
    dry_catchment = drainage & ~ocean & ~selected_wet
    xyz = np.asarray(fields.points, dtype=float)
    norm = np.linalg.norm(xyz, axis=1)
    latitude = np.degrees(np.arcsin(np.clip(xyz[:, 2] / np.maximum(norm, 1e-30), -1.0, 1.0)))
    longitude = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    phase_records[label] = {
        "step": step,
        "elapsed_years": float(world.elapsed_years),
        "anchor_count": int(len(anchor_points)),
        "wet_anchor_count": int(np.count_nonzero(anchor_labels >= 0)),
        "anchor_nearest_cell_distance_km_max": float(np.max(distance) * 6371.0),
        "selected_component_count": len(components),
        "components": components,
        "governing_lake_ids": sorted(selected_body_ids),
        "drainage_member_count": int(drainage.sum()),
        "drainage_member_area_km2": float(area_m2[drainage].sum() / 1e6),
        "dry_land_catchment_area_km2": float(area_m2[dry_catchment].sum() / 1e6),
        "connected_lake_volume_km3": float(sum(item["volume_km3"] for item in components)),
        "connected_lake_cell_footprint_area_km2": float(
            sum(item["cell_footprint_area_km2"] for item in components)
        ),
        "connected_lake_hypsometric_wet_area_km2": float(
            sum(item["hypsometric_wet_area_km2"] for item in components)
        ),
    }
    phase_arrays[label] = {
        "longitude": longitude.astype(np.float32),
        "latitude": latitude.astype(np.float32),
        "node_ids": np.asarray(fields.lake_trace_node_ids, dtype=np.uint64),
        "area_m2": np.asarray(area_m2, dtype=np.float64),
        "ocean": ocean,
        "drainage": drainage,
        "dry_catchment": dry_catchment,
        "lake": selected_wet,
    }
    print(json.dumps({"phase": label, **phase_records[label]}), flush=True)


if anchor_points is None and required_anchor_step == 0:
    raise ValueError("automatic anchor selection requires --anchor-step")

for step in range(1, last_step + 1):
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    fields = world.hydrology_cache
    wet, groups, _ = component_sets(fields)
    if args.anchor_step == step:
        anchor_points, selected_nodes = choose_automatic_anchors(fields, groups)
        anchor_selection = {
            "step": step,
            "reference_volume_km3": args.reference_volume_km3,
            "selected_volume_km3": float(
                np.dot(fields.lake_depth[selected_nodes], fields.hypsometry.area_m2[selected_nodes]) / 1e9
            ),
            "anchor_points_xyz": anchor_points.tolist(),
        }
    for label, phase_step in phases:
        if step == phase_step:
            if anchor_points is None:
                raise RuntimeError("anchor points were not available at the requested phase")
            capture_phase(label, step, fields)
    print(json.dumps({"step": step, "runtime_s": time.perf_counter() - started}), flush=True)

if len(phase_records) != len(phases):
    raise RuntimeError("not all requested phases were captured")


def save_outputs() -> None:
    arrays = {f"{label}__{name}": value for label, items in phase_arrays.items() for name, value in items.items()}
    np.savez_compressed(args.out, **arrays)
    transitions = []
    for (from_label, _), (to_label, _) in zip(phases, phases[1:]):
        old = phase_arrays[from_label]
        new = phase_arrays[to_label]
        old_index = {tuple(map(int, row)): idx for idx, row in enumerate(old["node_ids"])}
        new_index = {tuple(map(int, row)): idx for idx, row in enumerate(new["node_ids"])}
        old_members = {
            tuple(map(int, old["node_ids"][idx]))
            for idx in np.flatnonzero(old["drainage"])
        }
        new_members = {
            tuple(map(int, new["node_ids"][idx]))
            for idx in np.flatnonzero(new["drainage"])
        }
        added = new_members - old_members
        removed = old_members - new_members
        persistent = old_members & new_members
        transitions.append({
            "from": from_label,
            "to": to_label,
            "persistent_member_count": len(persistent),
            "added_member_count": len(added),
            "removed_member_count": len(removed),
            "added_member_area_km2": float(sum(new["area_m2"][new_index[node]] for node in added) / 1e6),
            "removed_member_area_km2": float(sum(old["area_m2"][old_index[node]] for node in removed) / 1e6),
            "persistent_area_change_km2": float(sum(
                new["area_m2"][new_index[node]] - old["area_m2"][old_index[node]]
                for node in persistent
            ) / 1e6),
        })
    summary = {
        "issue": 675,
        "commit": header["commit"],
        "save_sha256": header["save_sha256"],
        "step_years": args.step_years,
        "anchor_selection": anchor_selection,
        "phases": phase_records,
        "catchment_membership_transitions": transitions,
        "definition": (
            "Drainage catchment is the union of governing lake hierarchy body members whose wet nodes "
            "intersect anchor-containing connected lake components. Dry land catchment and lake water "
            "are plotted separately."
        ),
    }
    args.out.with_suffix(".json").write_text(json.dumps(summary, indent=2))
    if args.figure_out is not None:
        render_figure(args.figure_out)


def render_figure(path: Path) -> None:
    labels = [label for label, _ in phases]
    scale = 2
    panel_w, panel_h = 760, 500
    margin, gap, title_h, footer_h = 34, 18, 78, 88
    width = margin * 2 + len(labels) * panel_w + (len(labels) - 1) * gap
    height = title_h + panel_h + footer_h + margin
    image = Image.new("RGB", (width * scale, height * scale), "#fbfaf7")
    draw = ImageDraw.Draw(image)
    font_path = "/System/Library/Fonts/Supplemental/Arial.ttf"
    try:
        title_font = ImageFont.truetype(font_path, 27 * scale)
        label_font = ImageFont.truetype(font_path, 19 * scale)
        small_font = ImageFont.truetype(font_path, 14 * scale)
    except OSError:
        title_font = label_font = small_font = ImageFont.load_default()
    title = "Modeled drainage catchments through the lake flicker"
    draw.text((margin * scale, 18 * scale), title, fill="#243443", font=title_font)

    all_lat, all_lon = [], []
    for label in labels:
        arr = phase_arrays[label]
        anchor = np.asarray(anchor_points)
        center = anchor.mean(axis=0)
        center /= max(float(np.linalg.norm(center)), 1e-30)
        center_lat = math.asin(float(np.clip(center[2], -1.0, 1.0)))
        center_lon = math.atan2(float(center[1]), float(center[0]))
        lon = np.radians(arr["longitude"].astype(float))
        lat = np.radians(arr["latitude"].astype(float))
        delta_lon = (lon - center_lon + math.pi) % (2.0 * math.pi) - math.pi
        x_km = delta_lon * math.cos(center_lat) * 6371.0
        y_km = (lat - center_lat) * 6371.0
        catch = arr["drainage"]
        if np.any(catch):
            all_lon.extend(x_km[catch].tolist())
            all_lat.extend(y_km[catch].tolist())
        arr["plot_x"] = x_km
        arr["plot_y"] = y_km
    x_min, x_max = min(all_lon), max(all_lon)
    y_min, y_max = min(all_lat), max(all_lat)
    x_pad = max(150.0, (x_max - x_min) * 0.14)
    y_pad = max(150.0, (y_max - y_min) * 0.14)
    x_min, x_max = x_min - x_pad, x_max + x_pad
    y_min, y_max = y_min - y_pad, y_max + y_pad

    sea_color = "#e8f0f4"
    land_color = "#d5d2c9"
    catch_color = "#e7a64e"
    lake_color = "#2679a5"
    outline_color = "#34495a"
    for panel, label in enumerate(labels):
        x0 = margin + panel * (panel_w + gap)
        y0 = title_h
        left, top = x0 * scale, y0 * scale
        right, bottom = (x0 + panel_w) * scale, (y0 + panel_h) * scale
        draw.rounded_rectangle((left, top, right, bottom), radius=13 * scale, fill=sea_color, outline="#d1dbe1", width=2 * scale)
        inset = 32 * scale
        map_left, map_top = left + inset, top + 54 * scale
        map_right, map_bottom = right - inset, bottom - 56 * scale
        panel_title = f"{label} · {phase_records[label]['elapsed_years'] / 1e6:.1f} Myr"
        draw.text((left + 20 * scale, top + 14 * scale), panel_title, fill="#243443", font=label_font)
        arr = phase_arrays[label]
        x = arr["plot_x"]
        y = arr["plot_y"]
        visible = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)
        sx = map_left + (x - x_min) / (x_max - x_min) * (map_right - map_left)
        sy = map_bottom - (y - y_min) / (y_max - y_min) * (map_bottom - map_top)
        r = 2.2 * scale
        # Draw the base terrain first, then the catchment and its connected lake water.
        for mask, color, radius in (
            (visible & ~arr["ocean"], land_color, r),
            (visible & arr["dry_catchment"], catch_color, r * 1.45),
            (visible & arr["lake"], lake_color, r * 1.8),
        ):
            for px, py in zip(sx[mask], sy[mask]):
                pxf, pyf = float(px), float(py)
                draw.ellipse((pxf - radius, pyf - radius, pxf + radius, pyf + radius), fill=color)
        meta = phase_records[label]
        detail = (
            f"lake {meta['connected_lake_volume_km3']:,.0f} km³ · "
            f"lake-cell footprint {meta['connected_lake_cell_footprint_area_km2']:,.0f} km²"
        )
        draw.text((left + 20 * scale, bottom - 38 * scale), detail, fill="#34495a", font=small_font)

    legend_y = (title_h + panel_h + 24) * scale
    lx = margin * scale
    legend_items = [(catch_color, "Dry land in selected lake catchment"), (lake_color, "Connected lake water")]
    for color, label in legend_items:
        draw.ellipse((lx, legend_y, lx + 16 * scale, legend_y + 16 * scale), fill=color, outline=outline_color)
        draw.text((lx + 24 * scale, legend_y - 2 * scale), label, fill="#34495a", font=small_font)
        lx += 350 * scale
    draw.text(
        (margin * scale, (height - 34) * scale),
        "Same local equal-distance projection and map extent within each figure. Catchments are model hierarchy members, not measured watersheds.",
        fill="#566776", font=small_font,
    )
    image.resize((width, height), Image.Resampling.LANCZOS).save(path)


args.out.parent.mkdir(parents=True, exist_ok=True)
if args.figure_out is not None:
    args.figure_out.parent.mkdir(parents=True, exist_ok=True)
save_outputs()
print(f"wrote {args.out} and {args.out.with_suffix('.json')}")
if args.figure_out is not None:
    print(f"wrote {args.figure_out}")

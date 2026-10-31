#!/usr/bin/env python3
"""Run the preregistered three-arm 20 Myr evolving-slab mechanics screen (#699).

This offline screen initializes attached 3-D ribbons from the selected plate's currently
active bilateral subduction-front columns, transfers those columns out of the surface
inventory, advances body bending/advection and shell-driven plate strain, and accumulates
the registered M1 local damage rule. It does not step the production world or alter its
production split rule.

Usage (one process; RAM-heavy shell mechanics):
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
      backend/.venv/bin/python bin/debug/probe_gate3b_evolving_slabs.py \
      --save <registered.mbworld> \
      --out analysis/issue699/evolving-slabs-seed447276134-14.json \
      --image analysis/issue699/evolving-slabs-seed447276134-14.png
"""
from __future__ import annotations

import argparse
import copy
import gc
import json
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from app import geometry, lithosphere, mantle, persistence, torque
from app import slab_body as slab_model
import probe_gate3b_crack_coupled as damage_model
import probe_gate3b_component_shell as shell_audit
import probe_gate3b_midplate_current as current_audit
import probe_gate3b_slab_pull as slab_screen
import probe_gate3b_m1_saved_plates as saved_screen

PLATE_ID = 14
STEP_YEARS = 2_000_000.0
STEPS = 10
SHELL_CAP = 512
STATION_COUNT = 8
MAX_DEPTH_M = 660_000.0
DIP_DEG = 45.0
DENSITY_ANOMALY_KG_M3 = 80.0
BENDING_MODULUS_PA = 70e9
SURFACE_TEMP_K = 273.0
GEOTHERMAL_GRADIENT_K_PER_M = 3e-4
MANTLE_VISCOSITY_PA_S = 1e21
SHEAR_ZONE_M = 20_000.0
EXTENSIONAL_VISCOSITY_PA_S = 1e24
LOAD_SPREAD_DEG = 15.0
AMBIENT_MANTLE_DENSITY_KG_M3 = 3300.0

ARMS = (("A_no_slab_force", (False, False)),
        ("B_stronger_flank_only", (True, False)),
        ("C_both_flanks", (True, True)))


def _font(size: int):
    for candidate in ("/System/Library/Fonts/Supplemental/Arial.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    return ImageFont.load_default()


def _load_target(save_path: Path):
    with save_path.open("rb") as handle:
        world = persistence.load_world_file(handle)
    plate = next((p for p in world.plates if int(p.plate_id) == PLATE_ID), None)
    if plate is None:
        raise RuntimeError(f"registered plate {PLATE_ID} is absent from the saved world")
    points, areas, center, axis = slab_screen._principal_axis(plate)
    coord = points @ axis
    lo, hi = np.quantile(coord, [slab_screen.OUTER_FRACTION, 1.0 - slab_screen.OUTER_FRACTION])
    force_field, inputs, subducting = slab_screen._field_and_inputs(world, plate)
    slab_force = np.asarray(force_field.slab_pull_n)
    active = np.asarray(subducting, dtype=bool) & (np.linalg.norm(slab_force, axis=1) > 0.0)
    side_masks = (active & (coord <= lo), active & (coord >= hi))
    side_forces = [slab_force[mask].sum(axis=0) for mask in side_masks]
    outward = [-float(np.dot(side_forces[0], axis)), float(np.dot(side_forces[1], axis))]
    if not all(np.any(mask) for mask in side_masks) or min(outward) <= 0.0:
        raise RuntimeError("registered bilateral-front precondition failed on the frozen save")
    stronger = int(np.argmax(outward))

    native_centers = list(world.mantle_centers)
    native_p95 = float(mantle.rad_per_yr_to_cm_per_yr(np.percentile(
        np.linalg.norm(mantle.flow_at(points, native_centers), axis=1), 95.0
    )))
    source_geometry = current_audit._source_geometry(
        plate, shell_audit.line_spacing_rad(float(world.node_density))
    )
    source_centers, source_scale = current_audit._scaled_source_centers(
        source_geometry["source_centers_unit_strength"], points, native_p95
    )
    current_force = _midline_force(plate, source_centers)
    initial = {
        "world": world,
        "plate": plate,
        "points": points,
        "areas": areas,
        "center": center,
        "axis": axis,
        "coord": coord,
        "force_field": force_field,
        "inputs": inputs,
        "subducting": subducting,
        "active_masks": side_masks,
        "side_forces": side_forces,
        "stronger_side": stronger,
        "native_centers": native_centers,
        "native_p95_cm_yr": native_p95,
        "midline_source_centers": source_centers,
        "midline_source_scale": source_scale,
        "midline_force": current_force,
    }
    return initial


def _midline_force(plate, source_centers):
    nodes = plate.surface_nodes()
    points = geometry.normalize(np.asarray(nodes.world_xyz, dtype=float))
    velocity = mantle.flow_at(points, source_centers) * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    coefficient = torque.ASTHENOSPHERE_VISCOSITY_PA_S / torque.SHEAR_ZONE_THICKNESS_M
    return coefficient * np.asarray(nodes.area_m2)[:, None] * velocity


def _create_bodies(plate, initial):
    points = np.asarray(plate.surface_nodes().world_xyz, dtype=float)
    areas = np.asarray(plate.node_areas_m2(), dtype=float)
    hc = np.asarray(plate.collect("crustal_thickness_m"), dtype=float)
    hm = np.asarray(plate.collect("mantle_lithosphere_thickness_m"), dtype=float)
    crust_density = lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type)
    source_thickness = hc + hm
    source_mass_density = hc * crust_density + hm * AMBIENT_MANTLE_DENSITY_KG_M3
    spacing_m = shell_audit.line_spacing_rad(float(initial["world"].node_density)) * lithosphere.PLANET_RADIUS_M
    bodies = []
    source_records = []
    keys_seen: set[int] = set()
    for side_index, mask in enumerate(initial["active_masks"]):
        indices = np.flatnonzero(mask)
        centroid = geometry.normalize(np.average(points[indices], axis=0, weights=areas[indices]))
        side_axis = (-1.0 if side_index == 0 else 1.0) * initial["axis"]
        down_dip = geometry.normalize(side_axis - np.dot(side_axis, centroid) * centroid)
        anchor_index = int(indices[np.argmax(np.linalg.norm(initial["force_field"].slab_pull_n[indices], axis=1))])
        trench = geometry.normalize(points[anchor_index])
        source_area = float(areas[indices].sum())
        source_volume = float(np.dot(areas[indices], source_thickness[indices]))
        source_mass = float(np.dot(areas[indices], source_mass_density[indices]))
        width = source_area / spacing_m
        prototype = slab_model.create_slab_body(
            body_id=side_index,
            plate_id=plate.plate_id,
            plate_frame=plate.frame,
            trench_world_xyz=trench,
            down_dip_world_xyz=down_dip,
            source_cell_keys=plate.cell_keys[indices],
            width_m=width,
            thickness_m=100_000.0,
            max_depth_m=MAX_DEPTH_M,
            station_count=STATION_COUNT,
            dip_degrees=DIP_DEG,
            density_anomaly_kg_m3=DENSITY_ANOMALY_KG_M3,
            surface_temperature_k=SURFACE_TEMP_K,
            geothermal_gradient_k_per_m=GEOTHERMAL_GRADIENT_K_PER_M,
            bending_modulus_pa=BENDING_MODULUS_PA,
        )
        thickness = source_volume / (width * prototype.metrics()["length_m"])
        body = slab_model.create_slab_body(
            body_id=side_index,
            plate_id=plate.plate_id,
            plate_frame=plate.frame,
            trench_world_xyz=trench,
            down_dip_world_xyz=down_dip,
            source_cell_keys=plate.cell_keys[indices],
            width_m=width,
            thickness_m=thickness,
            max_depth_m=MAX_DEPTH_M,
            station_count=STATION_COUNT,
            dip_degrees=DIP_DEG,
            density_anomaly_kg_m3=DENSITY_ANOMALY_KG_M3,
            surface_temperature_k=SURFACE_TEMP_K,
            geothermal_gradient_k_per_m=GEOTHERMAL_GRADIENT_K_PER_M,
            bending_modulus_pa=BENDING_MODULUS_PA,
        )
        body.station_volume_m3[:] *= source_volume / body.consumed_volume_m3
        body.bending_rigidity_nm2 = BENDING_MODULUS_PA * body.width_m * body.thickness_m**3 / 12.0
        source_keys = set(map(int, body.source_cell_keys))
        if keys_seen.intersection(source_keys):
            raise RuntimeError("bilateral source columns overlap")
        keys_seen.update(source_keys)
        bodies.append(body)
        source_records.append({
            "side": "A" if side_index == 0 else "B",
            "source_node_count": int(len(indices)),
            "source_area_km2": source_area / 1e6,
            "source_volume_km3": source_volume / 1e9,
            "source_volume_m3": source_volume,
            "source_mass_kg": source_mass,
            "source_cell_keys": list(map(int, body.source_cell_keys)),
            "source_indices": indices,
            "anchor_index": anchor_index,
            "trench_world_xyz": trench.tolist(),
            "down_dip_world_xyz": down_dip.tolist(),
            "baseline_slab_force_n": initial["side_forces"][side_index].tolist(),
            "body_scenario_initial_force_n": (geometry.to_world(plate.frame, slab_model.attachment_force_local_n(body)[0])).tolist(),
            "slice_volume_km3": (body.station_volume_m3 / 1e9).tolist(),
            "slice_mass_kg": (np.full(STATION_COUNT, source_mass / STATION_COUNT)).tolist(),
            "body_geometry": body.metrics(),
        })
    # Transfer the source columns out of the plate's extensive lithosphere and tracer fields.
    source_indices = np.concatenate([record["source_indices"] for record in source_records])
    hc[source_indices] = 0.0
    hm[source_indices] = 0.0
    continental = plate.collect("continental_material_m")
    craton = plate.collect("craton_crust_m")
    continental[source_indices] = 0.0
    craton[source_indices] = 0.0
    plate.set_fields_on_plate(
        crustal_thickness_m=hc,
        mantle_lithosphere_thickness_m=hm,
        continental_material_m=continental,
        craton_crust_m=craton,
    )
    return bodies, source_records


def _body_force_rows(plate, bodies, enabled):
    nodes = plate.surface_nodes()
    points = geometry.normalize(np.asarray(nodes.world_xyz, dtype=float))
    force = np.zeros_like(points)
    by_side = []
    for body, is_enabled in zip(bodies, enabled):
        local, _reaction = slab_model.attachment_force_local_n(body)
        world = geometry.to_world(plate.frame, local)
        trench_world = geometry.normalize(geometry.to_world(plate.frame, body.trench_local_xyz))
        index = int(np.argmax(points @ trench_world))
        if is_enabled:
            force = world
        else:
            force = np.zeros(3)
        by_side.append({"enabled": bool(is_enabled), "force_n": force.tolist(),
                        "trench_node_index": index, "trench_xyz": trench_world.tolist()})
    return by_side


def _closing_speeds(world, plate, masks):
    _field, inputs, subducting = slab_screen._field_and_inputs(world, plate)
    rate = torque.boundary_closing_rate_m_per_s(plate, inputs)
    return [max(0.0, float(np.average(rate[mask & subducting], weights=plate.node_areas_m2()[mask & subducting])))
            if np.any(mask & subducting) else 0.0 for mask in masks]


def _surface_inventory(plate):
    areas = plate.node_areas_m2()
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    return {
        "crust_volume_m3": float(np.dot(areas, hc)),
        "mantle_lithosphere_volume_m3": float(np.dot(areas, hm)),
    }


def _run_arm(initial, body_templates, source_records, arm_name, enabled, capacity_name, capacity_n_per_m, template):
    plate = copy.deepcopy(template)
    bodies = copy.deepcopy(body_templates)
    world = initial["world"]
    original_plate = next(i for i, p in enumerate(world.plates) if int(p.plate_id) == PLATE_ID)
    world.plates[original_plate] = plate
    # Initial capacity and thickness are held in the material frame; source cells have
    # transferred their lithosphere into the ribbons and are excluded from crack initiation.
    initial_thickness = np.asarray(template._gate3b_initial_thickness, dtype=float)
    capacity = np.asarray(capacity_n_per_m, dtype=float).copy()
    source_indices = np.unique(np.concatenate([np.asarray(row["source_indices"], dtype=np.int64) for row in source_records]))
    capacity[source_indices] = 1e30
    epsilon = np.zeros(plate.node_count(), dtype=float)
    initial_craton = np.asarray(template._gate3b_initial_craton, dtype=bool)
    histories = []
    last = None
    total_error = 0.0
    for step in range(1, STEPS + 1):
        speeds = _closing_speeds(world, plate, initial["active_masks"])
        bend_rows = [slab_model.advance_slab_body(
            body, STEP_YEARS,
            mantle_viscosity_pa_s=MANTLE_VISCOSITY_PA_S,
            shear_zone_thickness_m=SHEAR_ZONE_M,
            convergence_speed_m_per_s=speeds[i],
        ) for i, body in enumerate(bodies)]
        side_rows = _body_force_rows(plate, bodies, enabled)
        current_force = _midline_force(plate, initial["midline_source_centers"])
        source_force_n = current_force.sum(axis=0)
        for side in side_rows:
            source_force_n += np.asarray(side["force_n"], dtype=float)
        strain_response = slab_model.apply_slab_load_to_plate_geometry(
            plate,
            [body for body, active in zip(bodies, enabled) if active],
            STEP_YEARS,
            effective_extensional_viscosity_pa_s=EXTENSIONAL_VISCOSITY_PA_S,
            load_spread_degrees=LOAD_SPREAD_DEG,
            extensional_capacity_resultant_n_per_m=capacity,
            additional_force_local_n=geometry.to_local(plate.frame, current_force),
        )
        demand = np.asarray(strain_response["max_tensile_resultant_n_per_m"], dtype=float)
        direction = geometry.normalize(np.asarray(strain_response["principal_extension_xyz"], dtype=float))
        stress_pa = demand / np.maximum(initial_thickness, 1.0)
        yield_pa = capacity / np.maximum(initial_thickness, 1.0)
        epsilon = damage_model._advance_damage(
            epsilon, stress_pa, yield_pa, STEP_YEARS * torque.SECONDS_PER_YEAR
        )
        cracked = epsilon >= damage_model.EPSILON_FAIL
        broken, components, labels, sizes = saved_screen._components(
            plate, cracked, direction, stress_pa, yield_pa
        )
        areas = np.asarray(plate.node_areas_m2(), dtype=float)
        component_area = np.bincount(labels, weights=areas, minlength=components)
        order = np.argsort(component_area)[::-1]
        area_fraction = component_area[order] / max(float(areas.sum()), 1.0)
        viable = bool(components == 2 and len(area_fraction) == 2 and np.all(area_fraction >= 0.10))
        volume_before = _surface_inventory(plate)
        volume_after = _surface_inventory(plate)
        surface_residual = max(
            abs(volume_after[key] - value) / max(abs(value), 1.0)
            for key, value in volume_before.items()
        )
        force_scale = max(float(strain_response["source_force_norm_sum_n"]), 1.0)
        torque_scale = max(float(strain_response["source_torque_norm_sum_nm"]), 1.0)
        force_resid = float(np.linalg.norm(strain_response["force_transfer_error_n"]) / force_scale)
        torque_resid = float(np.linalg.norm(strain_response["torque_transfer_error_nm"]) / torque_scale)
        equilibrium_resid = float(
            strain_response["equilibrium_residual_n"]
            / max(float(np.linalg.norm(strain_response["applied_resultant_force_n"])), 1.0)
        )
        work_resid = float(strain_response["work_residual_fraction"])
        material_residual = max(
            abs(sum(body.consumed_volume_m3 for body in bodies) - sum(float(r["source_volume_m3"]) for r in source_records))
            / max(sum(float(r["source_volume_m3"]) for r in source_records), 1.0),
            surface_residual,
        )
        body_work_residual = max(float(row["work_residual_fraction"]) for row in bend_rows)
        valid = bool(
            len(set(k for body in bodies for k in map(int, body.source_cell_keys)))
            == sum(len(body.source_cell_keys) for body in bodies)
            and material_residual <= 1e-10
            and body_work_residual <= 0.05
            and max(force_resid, torque_resid, equilibrium_resid, work_resid) <= 1e-3
        )
        centerline = [body.metrics() for body in bodies]
        total_error = max(total_error, force_resid, torque_resid, equilibrium_resid, work_resid,
                          body_work_residual, material_residual)
        checkpoint = {
            "elapsed_myr": step * 2,
            "arm": arm_name,
            "capacity_case": capacity_name,
            "slab_force_enabled_sides": [bool(v) for v in enabled],
            "source_force_n": source_force_n.tolist(),
            "midline_force_resultant_n": current_force.sum(axis=0).tolist(),
            "slab_force_sides": side_rows,
            "body_convergence_speed_m_per_s": speeds,
            "body_motion": bend_rows,
            "bodies": centerline,
            "shell": {
                "plate_area_coverage": 1.0,
                "force_residual_fraction": force_resid,
                "torque_residual_fraction": torque_resid,
                "equilibrium_residual_fraction": equilibrium_resid,
                "work_residual_fraction": work_resid,
                "body_relaxation_work_residual_fraction": body_work_residual,
                "surface_volume_residual_fraction": surface_residual,
                "max_log_strain": strain_response["max_log_strain"],
                "strain_increment": np.asarray(strain_response["strain_increment"]).tolist(),
                "load_work_j": float(strain_response["load_work_j"]),
                "stored_elastic_energy_j": float(strain_response["stored_elastic_energy_j"]),
            },
            "material_accounting": {
                "source_volume_m3": sum(float(r["source_volume_m3"]) for r in source_records),
                "slab_volume_m3": sum(body.consumed_volume_m3 for body in bodies),
                "deep_reservoir_volume_m3": sum(body.metrics()["deep_reservoir_volume_m3"] for body in bodies),
                "relative_residual": material_residual,
            },
            "damage": {
                "failed_cells": int(cracked.sum()),
                "failed_area_fraction": float(areas[cracked].sum() / max(areas.sum(), 1.0)),
                "max_epsilon": float(epsilon.max(initial=0.0)),
                "failed_links": int(len(broken)),
                "components": int(components),
                "two_largest_area_fractions": area_fraction[:2].tolist(),
                "viable_split": viable,
            },
            "accounting_valid": valid,
        }
        histories.append(checkpoint)
        last = {"plate": plate, "craton": initial_craton, "cracked": cracked.copy(),
                "broken": broken.copy(), "labels": labels.copy(), "areas": areas.copy(),
                "force_sides": side_rows, "checkpoint": checkpoint}
        if viable:
            break
        if not valid:
            break
    return {
        "arm": arm_name,
        "capacity_case": capacity_name,
        "history": histories,
        "first_viable_split_myr": next((row["elapsed_myr"] for row in histories if row["damage"]["viable_split"]), None),
        "status": "completed" if histories and histories[-1]["accounting_valid"] else "invalid-accounting",
        "final": last,
        "max_accounting_residual_fraction": total_error,
    }


def _map_xy(points, bounds, box):
    lon = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    lat = np.degrees(np.arcsin(np.clip(points[:, 2], -1.0, 1.0)))
    x0, x1, y0, y1 = bounds
    return np.column_stack((box[0] + (lon - x0) / (x1 - x0) * (box[2] - box[0]),
                            box[3] - (lat - y0) / (y1 - y0) * (box[3] - box[1])))


def _draw_montage(path, results):
    snapshots = [result["final"] for result in results]
    all_points = [geometry.normalize(s["plate"].surface_nodes().world_xyz) for s in snapshots]
    points = np.concatenate(all_points)
    lon = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    lat = np.degrees(np.arcsin(np.clip(points[:, 2], -1.0, 1.0)))
    center = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(lon))))))
    lon = (lon - center + 180.0) % 360.0 - 180.0 + center
    bounds = (float(lon.min() - 2), float(lon.max() + 2), float(lat.min() - 2), float(lat.max() + 2))
    width, height = 2700, 2250
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font, head_font, small_font = _font(36), _font(25), _font(16)
    draw.text((width // 2, 24), "Evolving slab ribbons + endogenous plate strain · 20 Myr", fill="#18212b", font=title_font, anchor="mt")
    draw.text((width // 2, 70), "Rows: p10 / median / p90 capacity · Columns: no slab force / one flank / both flanks", fill="#586878", font=head_font, anchor="mt")
    colors = ["#8ca7c2", "#6caa8d", "#aa8cc0", "#e0a45f", "#68aeb4", "#d27d83", "#b2ae5e", "#7087b2"]
    for index, snap in enumerate(snapshots):
        row, col = divmod(index, 3)
        xoff, yoff = col * 900, 115 + row * 700
        box = (xoff + 55, yoff + 75, xoff + 865, yoff + 590)
        result = results[index]
        data = snap["checkpoint"]
        largest = data["damage"]["two_largest_area_fractions"]
        summary = f"{result['arm']} · {result['capacity_case']} · {data['elapsed_myr']} Myr · {data['damage']['components']} components · daughters {largest[0]:.2%}/{largest[1]:.2%}"
        draw.text((xoff + 450, yoff + 16), summary, fill="#263442", font=head_font, anchor="mt")
        plate = snap["plate"]
        xyz = geometry.normalize(plate.surface_nodes().world_xyz)
        xy = _map_xy(xyz, bounds, box)
        labels = snap["labels"]
        for i, (x, y) in enumerate(xy):
            draw.ellipse((x-3, y-3, x+3, y+3), fill=colors[int(labels[i]) % len(colors)])
        craton = snap["craton"]
        for i in np.flatnonzero(craton):
            x, y = xy[i]
            draw.ellipse((x-7, y-7, x+7, y+7), fill="#f2c14e", outline="#372f20", width=2)
        for i in np.flatnonzero(snap["cracked"]):
            x, y = xy[i]
            draw.line((x-5, y-5, x+5, y+5), fill="#b51f2e", width=3)
            draw.line((x-5, y+5, x+5, y-5), fill="#b51f2e", width=3)
        for a, b in snap["broken"]:
            draw.line((*xy[int(a)], *xy[int(b)]), fill="#f07818", width=4)
        for side, color in zip(snap["force_sides"], ("#152b43", "#00796b")):
            if not side["enabled"]:
                continue
            vec = np.asarray(side["force_n"], dtype=float)
            norm = np.linalg.norm(vec)
            if norm <= 0:
                continue
            anchor = geometry.normalize(np.asarray(side["trench_xyz"], dtype=float))
            lat_r = np.arcsin(np.clip(anchor[2], -1, 1))
            lon_r = np.arctan2(anchor[1], anchor[0])
            east = np.array([-np.sin(lon_r), np.cos(lon_r), 0.0])
            north = np.array([-np.sin(lat_r)*np.cos(lon_r), -np.sin(lat_r)*np.sin(lon_r), np.cos(lat_r)])
            fx, fy = float(vec @ east), float(vec @ north)
            direction = np.array([fx, -fy])
            direction /= max(np.linalg.norm(direction), 1.0)
            start = xy[side["trench_node_index"]]
            end = start + direction * 70.0
            draw.line((*start, *end), fill=color, width=6)
            angle = np.arctan2(end[1]-start[1], end[0]-start[0])
            draw.polygon((tuple(end), (end[0]-18*np.cos(angle-.5), end[1]-18*np.sin(angle-.5)),
                          (end[0]-18*np.cos(angle+.5), end[1]-18*np.sin(angle+.5))), fill=color)
        draw.rectangle(box, outline="#aab5bf", width=1)
    legend_y = 2210
    draw.text((45, legend_y), "No craton cells occur on registered plate 14. Red crosses: failed cells   Orange: removed failed links   Navy/teal arrows: slab pull sides A/B (direction; length schematic)", fill="#263442", font=small_font)
    image.save(path)


def run(save_path: Path, out_path: Path, image_path: Path):
    started = time.perf_counter()
    initial = _load_target(save_path)
    base_world = initial["world"]
    base_plate = initial["plate"]
    original_thickness = base_plate.collect("crustal_thickness_m") + base_plate.collect("mantle_lithosphere_thickness_m")
    original_craton = base_plate.collect("craton_crust_m") > 0.0
    original_continental = base_plate.collect("continental_material_m") > 0.0
    base_plate._gate3b_initial_thickness = np.asarray(original_thickness, dtype=float).copy()
    base_plate._gate3b_initial_craton = np.asarray(original_craton, dtype=bool).copy()
    bodies, source_records = _create_bodies(base_plate, initial)
    template = copy.deepcopy(base_plate)
    # Retain material-frame damage support and the original craton mask on all clones.
    template._gate3b_initial_thickness = np.asarray(original_thickness, dtype=float).copy()
    template._gate3b_initial_craton = np.asarray(original_craton, dtype=bool).copy()
    capacity_fields, capacity_reference = damage_model._capacity_fields(original_continental)
    results = []
    for capacity_name in ("p10", "median", "p90"):
        for arm_name, enabled_flags in ARMS:
            # Every capacity/force arm starts from the exact same transferred snapshot.
            result = _run_arm(initial, bodies, source_records, arm_name, enabled_flags,
                              capacity_name, capacity_fields[capacity_name], template)
            results.append(result)
            gc.collect()
    _draw_montage(image_path, results)
    clean_results = [{k: v for k, v in result.items() if k != "final"} for result in results]
    payload = {
        "schema_version": 1,
        "title": "RFC 0010 Gate 3b evolving slab ribbon mini-simulation",
        "scope": "Offline implementation/mechanics screen; not production world stepping or Gate 3b closure.",
        "seed": int(base_world.seed),
        "plate_id": PLATE_ID,
        "saved_elapsed_myr": float(base_world.elapsed_years / 1e6),
        "node_density": float(base_world.node_density),
        "native_midline_p95_cm_per_year": initial["native_p95_cm_yr"],
        "midline_source_scale": initial["midline_source_scale"],
        "stronger_flank": "A" if initial["stronger_side"] == 0 else "B",
        "initial_source_records": [
            {k: v for k, v in record.items() if k != "source_indices"}
            for record in source_records
        ],
        "initial_surface_inventory_after_transfer": _surface_inventory(base_plate),
        "capacity_references_n_per_m": capacity_reference,
        "parameters": {
            "step_myr": 2.0, "steps": STEPS, "body_dip_degrees": DIP_DEG,
            "body_max_depth_m": MAX_DEPTH_M, "body_station_count": STATION_COUNT,
            "density_anomaly_kg_m3": DENSITY_ANOMALY_KG_M3,
            "bending_modulus_pa": BENDING_MODULUS_PA,
            "surface_temperature_k": SURFACE_TEMP_K,
            "geothermal_gradient_k_per_m": GEOTHERMAL_GRADIENT_K_PER_M,
            "mantle_viscosity_pa_s": MANTLE_VISCOSITY_PA_S,
            "shear_zone_m": SHEAR_ZONE_M,
            "shell_extensional_viscosity_pa_s": EXTENSIONAL_VISCOSITY_PA_S,
            "load_spread_degrees": LOAD_SPREAD_DEG,
            "damage_law": {
                "epsilon_onset": damage_model.EPSILON_ONSET,
                "epsilon_fail": damage_model.EPSILON_FAIL,
                "viscosity_pa_s": damage_model.DAMAGE_VISCOSITY_PA_S,
                "residual_strength_fraction": damage_model.RESIDUAL_STRENGTH_FRACTION,
            },
            "midline": "M1 geometric chord source scaled so source-only mantle p95 matches the saved world's native p95",
        },
        "results": clean_results,
        "runtime_seconds": time.perf_counter() - started,
        "image": image_path.name,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", type=Path, required=True, help="registered .mbworld input")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--image", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.save, args.out, args.image)
    print(json.dumps({
        "output": str(args.out), "image": str(args.image),
        "runtime_seconds": result["runtime_seconds"],
        "results": [{
            "arm": row["arm"], "capacity": row["capacity_case"],
            "status": row["status"], "split_myr": row["first_viable_split_myr"],
            "final": row["history"][-1] if row["history"] else None,
        } for row in result["results"]],
    }, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the Gate 3b M1 damage diagnostic for one saved candidate plate.

This is a frozen-geometry diagnostic: it loads one saved world, reconstructs its
production slab-pull field, adds the same geometrically centered midline source used by
the seed-2024 M1 experiment, ramps the loads for 20 Myr, and advances the Gate 3a
damage law. It does not step the full tectonic world or invoke production rifting.

Usage (run one case per process to limit memory):
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
      backend/.venv/bin/python bin/debug/probe_gate3b_m1_saved_plates.py \
      --case 447276134-14 --out analysis/issue699/m1-saved-447276134-14.json \
      --image analysis/issue699/m1-saved-447276134-14.png
"""
from __future__ import annotations

import argparse
import gc
import json
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from app import cratons, elevation_lines, geometry, lithosphere, orogeny, persistence, torque
import gate3b_diagnostic_common as shared
import probe_gate3b_crack_coupled as diagnostic
import probe_gate3b_component_shell as shell_audit
import probe_gate3b_midplate_current as current_audit
import probe_gate3b_slab_pull as slab_screen
import validate_thin_shell_plate_meshes as mesh_validation

CASES = {
    "447276134-14": ("/Users/adubey/Downloads/mantle-bloom-seed447276134-116300000y.mbworld", 14, 116.3),
    "373936961-0": ("/Users/adubey/Downloads/mantle-bloom-seed373936961-244400000y.mbworld", 0, 244.4),
    "295569319-11": ("/Users/adubey/Downloads/mantle-bloom-seed295569319-180000000y.mbworld", 11, 180.0),
    "373936961-40": ("/Users/adubey/Downloads/mantle-bloom-seed373936961-244400000y.mbworld", 40, 244.4),
}
STEP_YEARS = 2_000_000.0
DURATION_MYR = 20
LOAD_STAGES = ((0, 1.0), (10, 2.0), (20, 4.0))
MIDLINE_MULTIPLIER = 1.0


_candidate_triangles = mesh_validation._candidate_triangles


def _candidate_triangles_without_degenerate_slivers(plate, points):
    """Drop accepted Delaunay slivers that the shell element solver cannot integrate."""
    triangles, rejected, candidates = _candidate_triangles(plate, points)
    xyz = geometry.normalize(np.asarray(points, dtype=float))
    valid = np.ones(len(triangles), dtype=bool)
    radius = lithosphere.PLANET_RADIUS_M
    for index, triangle in enumerate(triangles):
        radial = geometry.normalize(xyz[triangle].sum(axis=0))
        edge = xyz[triangle[1]] - xyz[triangle[0]]
        axis_x = geometry.normalize(edge - np.dot(edge, radial) * radial)
        axis_y = np.cross(radial, axis_x)
        xy = np.column_stack((xyz[triangle] @ axis_x, xyz[triangle] @ axis_y)) * radius
        x, y = xy[:, 0], xy[:, 1]
        twice_area = (x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0])
        valid[index] = abs(twice_area) * 0.5 > np.finfo(float).eps * radius**2
    return triangles[valid], rejected + int((~valid).sum()), candidates


# Saved cell meshes can contain a handful of accepted near-zero-area Delaunay slivers.
# Remove only elements below the exact shell solver's degeneracy threshold for this run.
mesh_validation._candidate_triangles = _candidate_triangles_without_degenerate_slivers


def _target(world, plate_id):
    plate = next((p for p in world.plates if int(p.plate_id) == plate_id), None)
    if plate is None:
        raise RuntimeError(f"saved world does not contain plate {plate_id}")
    points, areas, center, axis = slab_screen._principal_axis(plate)
    candidate = {"points": points, "areas": areas, "center_xyz": center,
                 "axis_xyz": axis, "coord": points @ axis}
    (_native_centers, native_p95, _source_geometry, _centers, scale, field, inputs,
     subducting, midline_force) = shared.add_native_midline_source(
        world, plate, candidate["points"], current_audit._source_geometry,
        current_audit._scaled_source_centers, slab_screen._field_and_inputs,
        shell_audit.line_spacing_rad(float(world.node_density)),
    )
    bilateral = slab_screen._bilateral_candidate(world, plate)
    return plate, candidate, bilateral, native_p95, scale, field, midline_force, inputs, subducting


def _stage_for(elapsed_myr):
    multiplier = LOAD_STAGES[0][1]
    for start, value in LOAD_STAGES:
        if elapsed_myr >= start:
            multiplier = value
    return multiplier


def _components(plate, cracked, direction, stress, yield_strength):
    broken, count, labels, sizes, _pairs = shared.local_failure_components(
        plate, cracked, direction, stress, yield_strength,
        diagnostic.RESIDUAL_STRENGTH_FRACTION,
    )
    return broken, count, labels, sizes


def _map_coords(points):
    return shared.map_spherical_coordinates(points)


def _draw_image(path, title, plate, candidate, field, subducting, snapshot):
    points, _ = plate.all_points_and_elevation()
    xyz, lon, lat = _map_coords(points)
    labels = snapshot["labels"]
    sizes = snapshot["sizes"]
    component_colors = ["#91b4d4", "#70a98b", "#ab91c6", "#e2a86f", "#74b9bf", "#d98187", "#b4b466", "#7390bd"]
    width, height = 1800, 1180
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 30)
    sub_font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 20)
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 17)
    draw.text((width // 2, 30), title, fill="#18212b", font=title_font, anchor="mt")
    ncomp = len(sizes)
    largest = np.sort(sizes)[::-1]
    subtitle = (f"20 Myr · p10 strength · {ncomp} components · largest {largest[0]:,} cells"
                + (f" / {largest[1]:,}" if len(largest) > 1 else " / 0")
                + f" · {snapshot['failed_cells']:,} cracked cells ({snapshot['failed_area_fraction']:.1%})")
    if ncomp == 1:
        subtitle = subtitle.replace("1 components", "1 component")
    draw.text((width // 2, 75), subtitle, fill="#586878", font=sub_font, anchor="mt")
    pad = 5.0
    x0, x1 = float(lon.min() - pad), float(lon.max() + pad)
    y0, y1 = float(lat.min() - pad), float(lat.max() + pad)
    map_box = (75, 125, 1330, 1050)
    def xy(lon_v, lat_v):
        return (int(map_box[0] + (lon_v-x0)/(x1-x0)*(map_box[2]-map_box[0])),
                int(map_box[3] - (lat_v-y0)/(y1-y0)*(map_box[3]-map_box[1])))
    for i in range(len(points)):
        x, y = xy(lon[i], lat[i])
        color = component_colors[int(labels[i]) % len(component_colors)]
        draw.ellipse((x-4, y-4, x+4, y+4), fill=color)
    craton = plate.collect("craton_crust_m") > 0
    for i in np.flatnonzero(craton):
        x, y = xy(lon[i], lat[i])
        draw.ellipse((x-7, y-7, x+7, y+7), fill="#f2c14e", outline="#372f20", width=2)
    for i in np.flatnonzero(snapshot["cracked"]):
        x, y = xy(lon[i], lat[i])
        draw.line((x-5,y-5,x+5,y+5), fill="#b51f2e", width=2)
        draw.line((x-5,y+5,x+5,y-5), fill="#b51f2e", width=2)
    # Sample force-bearing nodes independently from the two opposing plate margins.
    force = np.asarray(field.slab_pull_n)
    coord = np.asarray(candidate["coord"])
    lo, hi = np.quantile(coord, [slab_screen.OUTER_FRACTION, 1-slab_screen.OUTER_FRACTION])
    for side, mask in enumerate((coord <= lo, coord >= hi)):
        indices = np.flatnonzero(np.asarray(subducting, dtype=bool) & mask & (np.linalg.norm(force, axis=1) > 0))
        if not len(indices):
            continue
        indices = indices[np.argsort(np.linalg.norm(force[indices], axis=1))[-min(7, len(indices)):]]
        for index in indices:
            lat_r = np.arcsin(np.clip(xyz[index, 2], -1, 1))
            lon_r = np.arctan2(xyz[index,1], xyz[index,0])
            east = np.array([-np.sin(lon_r), np.cos(lon_r), 0.0])
            north = np.array([-np.sin(lat_r)*np.cos(lon_r), -np.sin(lat_r)*np.sin(lon_r), np.cos(lat_r)])
            fx, fy = float(force[index] @ east), float(force[index] @ north)
            norm = max(float(np.hypot(fx, fy)), 1.0)
            # Scale vectors to a readable fixed map length; direction is the measurement shown.
            length = 2.6
            u = fx / norm * length / max(np.cos(lat_r), .25)
            v = fy / norm * length
            xa, ya = xy(lon[index], lat[index]); xb, yb = xy(lon[index]+u, lat[index]+v)
            color = ("#152b43", "#00796b")[side]
            draw.line((xa,ya,xb,yb), fill=color, width=5)
            angle = np.arctan2(yb-ya, xb-xa)
            head = 14
            draw.polygon(((xb,yb),(xb-head*np.cos(angle-.5),yb-head*np.sin(angle-.5)),
                          (xb-head*np.cos(angle+.5),yb-head*np.sin(angle+.5))), fill=color)
    for a, b in snapshot["broken"]:
        if labels[a] != labels[b]:
            pa, pb = xy(lon[a],lat[a]), xy(lon[b],lat[b])
            draw.line((*pa,*pb), fill="#4a1a09", width=10)
            draw.line((*pa,*pb), fill="#ff6d00", width=6)
    # Component size key helps distinguish isolated cells from substantial fragments.
    draw.rounded_rectangle((1370, 135, 1765, 820), radius=12, fill="#fbfcfd", outline="#d5dce3", width=2)
    draw.text((1395,155), "Component sizes", fill="#243447", font=sub_font)
    order = np.argsort(sizes)[::-1]
    for rank, label in enumerate(order[:20], start=1):
        y = 195 + (rank-1)*29
        draw.text((1400,y), f"{rank:>2}. {sizes[label]:,} cells", fill="#354455", font=font)
    if len(order) > 20:
        draw.text((1400, 195+20*29), f"… {len(order)-20} smaller components", fill="#586878", font=font)
    draw.text((1395, 850), "Legend", fill="#243447", font=sub_font)
    legends = [("#f2c14e","craton"),("#b51f2e","locally failed cell"),("#ff6d00","failed link"),("#152b43","slab pull side A"),("#00796b","slab pull side B"),("#91b4d4","component area")]
    for i, (color,label) in enumerate(legends):
        y = 890 + i*27
        draw.ellipse((1398,y,1410,y+12), fill=color)
        draw.text((1420,y), label, fill="#354455", font=font)
    draw.text((width//2, 1120),
              "Arrows show sampled slab-pull directions (fixed display length); orange links separate components.",
              fill="#586878", font=font, anchor="mm")
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def run(case_name, out_path, image_path):
    started = time.perf_counter()
    save_path, plate_id, save_time_myr = CASES[case_name]
    with open(save_path, "rb") as handle:
        world = persistence.load_world_file(handle)
    cratons.ensure_initialized(world)
    plate, candidate, bilateral, native_p95, source_scale, field, midline_force, inputs, subducting = _target(world, plate_id)
    source_base = np.asarray(field.total_n, dtype=float)
    demand_by_factor, direction_by_factor, shell_by_factor = {}, {}, {}
    for slab_multiplier in (1.0, 2.0, 4.0):
        scaled_force = source_base + (slab_multiplier - 1.0) * field.slab_pull_n + (MIDLINE_MULTIPLIER-1.0) * midline_force
        shell, demand, direction = diagnostic._opening_capacity_demand(plate, scaled_force)
        demand_by_factor[slab_multiplier] = demand
        direction_by_factor[slab_multiplier] = direction
        shell["demand_quantiles_n_per_m"] = np.quantile(demand, [0.1, .5, .9]).tolist()
        shell_by_factor[str(slab_multiplier)] = {key: shell.get(key) for key in (
            "target_representatives", "shell_node_count", "mesh_triangles", "production_plate_area_coverage",
            "accounting_pass", "force_transfer_residual_fraction", "torque_transfer_residual_fraction",
            "equilibrium_residual_fraction", "shell_energy_residual_fraction", "demand_quantiles_n_per_m")}
    thickness = np.maximum(plate.collect("crustal_thickness_m") + plate.collect("mantle_lithosphere_thickness_m"), 1.0)
    areas = plate.node_areas_m2()
    continental = elevation_lines.effective_is_continental_from_codes(
        plate.collect("crust_type_code"), plate.crust_type == "continental"
    )
    caps, cap_refs = diagnostic._capacity_fields(continental)
    outcomes = {}
    snapshots = {}
    for case, cap in caps.items():
        epsilon = np.zeros(plate.node_count())
        seconds = STEP_YEARS * 365.25 * 24.0 * 3600.0
        for step in range(1, int(DURATION_MYR * 1e6 / STEP_YEARS) + 1):
            elapsed = step * STEP_YEARS / 1e6
            multiplier = _stage_for(elapsed)
            stress = demand_by_factor[multiplier] / thickness
            direction = direction_by_factor[multiplier]
            yield_strength = cap / thickness
            epsilon = diagnostic._advance_damage(epsilon, stress, yield_strength, seconds)
        cracked = epsilon >= diagnostic.EPSILON_FAIL
        broken, count, labels, sizes = _components(plate, cracked, direction, stress, yield_strength)
        size_order = np.argsort(sizes)[::-1]
        area_by_label = np.bincount(labels, weights=areas, minlength=count)
        area_order = np.argsort(area_by_label)[::-1]
        summary = diagnostic._connectivity_summary(plate, cracked, direction, stress, yield_strength, areas)
        summary.update({
            "component_cell_counts_largest_first": sizes[size_order].astype(int).tolist(),
            "component_area_fractions_largest_first": (area_by_label[area_order]/areas.sum()).tolist(),
            "two_largest_component_cells": [int(sizes[size_order[i]]) if i < count else 0 for i in range(2)],
            "viable_split_candidate": bool(
                count == 2
                and area_by_label[area_order[1]] >= diagnostic.MIN_DAUGHTER_AREA_FRACTION * areas.sum()
                and area_by_label[area_order[:2]].sum() >= diagnostic.MIN_INTACT_AREA_FRACTION * areas.sum()
            ),
            "failed_cell_area_km2": float(areas[cracked].sum()/1e6),
            "peak_damage_strain": float(epsilon.max()),
        })
        outcomes[case] = summary
        snapshots[case] = {"cracked": cracked, "broken": broken, "labels": labels, "sizes": sizes,
                           "failed_cells": int(cracked.sum()),
                           "failed_area_fraction": float(areas[cracked].sum()/areas.sum())}
    _draw_image(image_path, f"M1 saved-plate diagnostic · seed {case_name.split('-')[0]} · plate {plate_id}",
                plate, candidate, field, subducting, snapshots["p10"])
    magnitudes = np.linalg.norm(np.asarray(field.slab_pull_n), axis=1)
    coordinate = np.asarray(candidate["coord"])
    lo, hi = np.quantile(coordinate, [slab_screen.OUTER_FRACTION, 1-slab_screen.OUTER_FRACTION])
    side_metrics = []
    for mask in (coordinate <= lo, coordinate >= hi):
        ids = np.flatnonzero(np.asarray(subducting, dtype=bool) & mask & (magnitudes > 0))
        side_metrics.append({"active_nodes": int(len(ids)), "slab_pull_resultant_n": np.asarray(field.slab_pull_n)[ids].sum(axis=0).tolist(),
                            "slab_pull_resultant_magnitude_n": float(np.linalg.norm(np.asarray(field.slab_pull_n)[ids].sum(axis=0)))})
    side_vectors = [np.asarray(row["slab_pull_resultant_n"]) for row in side_metrics]
    side_magnitudes = [float(np.linalg.norm(vector)) for vector in side_vectors]
    slab_angle = float(np.degrees(np.arccos(np.clip(
        np.dot(side_vectors[0], side_vectors[1]) / max(side_magnitudes[0]*side_magnitudes[1], 1.0), -1.0, 1.0
    ))))
    areas = plate.node_areas_m2()
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    rho_c = lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type)
    column_mass = areas * (hc * rho_c + hm * 3300.0)
    moho_c = orogeny.moho_temperature_c(hc, hm, plate.collect("moho_thermal_lag_c"))
    front_length = torque.boundary_node_length_m(inputs, shell_audit.line_spacing_rad(float(world.node_density)))
    slice_edges = np.linspace(float(coordinate.min()), float(coordinate.max()), 11)
    slice_ids = np.clip(np.digitize(coordinate, slice_edges[1:-1]), 0, 9)
    slice_metrics = []
    for slice_id in range(10):
        mask = slice_ids == slice_id
        active = mask & np.asarray(subducting, dtype=bool)
        slice_metrics.append({
            "slice": slice_id,
            "mass_1e21_kg": float(column_mass[mask].sum()/1e21),
            "crust_volume_1e6_km3": float(np.sum(areas[mask]*hc[mask])/1e15),
            "mantle_lithosphere_volume_1e6_km3": float(np.sum(areas[mask]*hm[mask])/1e15),
            "mean_total_thickness_km": float(np.average((hc[mask]+hm[mask])/1000, weights=areas[mask])),
            "mean_active_slab_column_thickness_km": float(np.average((hc[active]+hm[active])/1000, weights=areas[active])) if active.any() else None,
            "mean_moho_temperature_c": float(np.average(moho_c[mask], weights=areas[mask])),
            "active_slab_front_length_km": float(front_length[active].sum()/1000),
            "slab_pull_resultant_1e18_n": (np.asarray(field.slab_pull_n)[mask].sum(axis=0)/1e18).tolist(),
            "midline_basal_drag_resultant_1e18_n": (np.asarray(midline_force)[mask].sum(axis=0)/1e18).tolist(),
        })
    result = {
        "case": case_name, "seed": int(case_name.split("-")[0]), "plate_id": plate_id,
        "save_time_myr": save_time_myr, "node_density": float(world.node_density),
        "diagnostic_duration_myr": DURATION_MYR, "step_myr": STEP_YEARS/1e6,
        "status": "frozen-geometry M1 diagnostic; not full tectonic-world evolution",
        "plate_area_km2": float(areas.sum()/1e6), "node_count": int(plate.node_count()),
        "continental_area_fraction": float(areas[continental].sum()/areas.sum()),
        "craton_area_fraction": float(areas[plate.collect("craton_crust_m") > 0].sum()/areas.sum()),
        "plate_mass_1e21_kg": float(column_mass.sum()/1e21),
        "crust_volume_1e6_km3": float(np.sum(areas*hc)/1e15),
        "mantle_lithosphere_volume_1e6_km3": float(np.sum(areas*hm)/1e15),
        "mean_total_thickness_km": float(np.average((hc+hm)/1000, weights=areas)),
        "mean_active_slab_column_thickness_km": float(np.average((hc[subducting]+hm[subducting])/1000, weights=areas[subducting])) if np.any(subducting) else None,
        "mean_moho_temperature_c": float(np.average(moho_c, weights=areas)),
        "native_mantle_flow_p95_cm_per_yr": native_p95, "midline_source_scale": float(source_scale),
        "slab_sink_length_km": float(torque.SUBDUCTION_LSINK_M/1000),
        "active_slab_front_length_km": float(torque.boundary_node_length_m(inputs, shell_audit.line_spacing_rad(float(world.node_density)))[subducting].sum()/1000),
        "bilateral_slab_sides": side_metrics,
        "axis_end_slab_resultant_separation_deg": slab_angle,
        "weaker_to_stronger_axis_end_slab_resultant_ratio": min(side_magnitudes)/max(side_magnitudes) if max(side_magnitudes) > 0 else 0.0,
        "ten_equal_axis_slices": slice_metrics,
        "meets_lighter_neighbor_bilateral_precondition": bilateral is not None,
        "load_schedule": [{"elapsed_myr": start, "slab_multiplier": value, "midline_multiplier": MIDLINE_MULTIPLIER} for start,value in LOAD_STAGES],
        "shell_accounting_by_slab_multiplier": shell_by_factor,
        "capacity_references_n_per_m": cap_refs,
        "capacity_cases": outcomes,
        "damage_law": {"epsilon_onset": diagnostic.EPSILON_ONSET, "epsilon_fail": diagnostic.EPSILON_FAIL,
                       "viscosity_pa_s": diagnostic.DAMAGE_VISCOSITY_PA_S,
                       "residual_strength_fraction": diagnostic.RESIDUAL_STRENGTH_FRACTION},
        "split_rule": {"required_components": 2,
                       "minimum_daughter_area_fraction": diagnostic.MIN_DAUGHTER_AREA_FRACTION,
                       "minimum_retained_two_daughter_area_fraction": diagnostic.MIN_INTACT_AREA_FRACTION},
        "force_resultants_n": {"native_slab_pull": np.asarray(field.slab_pull_n).sum(axis=0).tolist(),
                               "native_total_production": np.asarray(field.total_n).sum(axis=0).tolist(),
                               "midline_current": np.asarray(midline_force).sum(axis=0).tolist()},
        "elapsed_wall_seconds": time.perf_counter()-started,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"case": case_name, "components": {k: v["plate_components_after_failed_links"] for k,v in outcomes.items()},
                      "two_largest": {k:v["two_largest_component_cells"] for k,v in outcomes.items()},
                      "image": str(image_path), "wall_seconds": result["elapsed_wall_seconds"]}, sort_keys=True))
    del world
    gc.collect()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=sorted(CASES), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--image", type=Path, required=True)
    args = parser.parse_args()
    run(args.case, args.out, args.image)


if __name__ == "__main__":
    main()

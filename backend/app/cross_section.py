"""Surface-sampled shallow lithosphere sections (#360).

The section is a narrow cut through the existing surface-node columns. It deliberately does
not construct a volumetric mantle: only Hc/Hm columns, their illustrative thermal profile,
surface vents/flow, and the schematic asthenosphere below the modeled lithosphere are drawn.
"""

from __future__ import annotations

import base64
import io

import numpy as np
from PIL import Image, ImageDraw

from . import geometry, mantle, orogeny, plates
from .elevation_lines import (
    CRUST_TYPE_CONTINENTAL,
    CRUST_TYPE_INHERIT,
    PLANET_RADIUS_KM,
    line_spacing_rad,
)
from .hydroclimate import biomes, climate

ARC_LENGTH_KM = 3000.0
ARC_POSITION_COUNT = 120
SECTION_DEPTH_KM = 350.0
VERTICAL_EXAGGERATION = 2.0
CRUST_CONTINENTAL = (201, 143, 130)
CRUST_OCEANIC = (82, 103, 121)
MANTLE_LITHOSPHERE = (66, 90, 90)
ASTHENOSPHERE = (117, 108, 148)
MAGMA_PALETTE = ((201, 61, 61), (244, 122, 60), (255, 214, 90))
SEA_LEVEL_COLOR = (83, 159, 203)


def _world_surface(world) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    """Return one plate-major snapshot of active node centers and fields."""
    live_plates = [plate for plate in world.plates if plate.node_count()]
    if not live_plates:
        return np.zeros((0, 3)), {}, np.zeros(0, dtype=bool)
    xyz = np.concatenate([plate.all_points_and_elevation()[0] for plate in live_plates])
    fields = {
        name: np.concatenate([plate.collect(name) for plate in live_plates])
        for name in (
            "elevation",
            "crustal_thickness_m",
            "mantle_lithosphere_thickness_m",
            "crust_type_code",
            "glacier_depth",
            "is_volcano",
            "volcano_active_years_remaining",
            "moho_thermal_lag_c",
        )
    }
    plate_continental = np.concatenate(
        [np.full(plate.node_count(), plate.crust_type == "continental", dtype=bool) for plate in live_plates]
    )
    codes = fields["crust_type_code"]
    continental = np.where(codes == CRUST_TYPE_INHERIT, plate_continental, codes == CRUST_TYPE_CONTINENTAL)
    return xyz, fields, continental


def _great_circle_samples(plane_angle_deg: float, position_index: int, spacing_km: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    plane = np.radians(plane_angle_deg)
    center_angle = 2.0 * np.pi * position_index / ARC_POSITION_COUNT
    half_angle = ARC_LENGTH_KM / (2.0 * PLANET_RADIUS_KM)
    count = max(3, int(np.ceil(ARC_LENGTH_KM / max(spacing_km * 0.5, 1.0))) + 1)
    angles = np.linspace(center_angle - half_angle, center_angle + half_angle, count)
    plane_axis = np.array([np.cos(plane), 0.0, np.sin(plane)])
    along_axis = np.array([0.0, 1.0, 0.0])
    points = np.cos(angles)[:, None] * plane_axis + np.sin(angles)[:, None] * along_axis
    tangent = -np.sin(angles)[:, None] * plane_axis + np.cos(angles)[:, None] * along_axis
    distance_km = (angles - center_angle) * PLANET_RADIUS_KM
    return points, tangent, distance_km


def _sample_fields(world, points: np.ndarray) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    node_xyz, fields, continental = _world_surface(world)
    if not len(node_xyz):
        raise ValueError("world has no active surface nodes")
    _, nearest = plates.cached_node_position_tree(world, node_xyz).query(points, k=1)
    sampled = {name: values[nearest] for name, values in fields.items()}
    return sampled, continental[nearest], node_xyz


def _biome_colors(world, points: np.ndarray, elevations: np.ndarray) -> np.ndarray:
    """Sample the current climate map's natural palette, falling back to latitude/elevation."""
    try:
        climate_fields = climate.compute_climate_cached(world)
        lat_rad, lon_rad = geometry.xyz_to_latlon(points)
        grid_height, grid_width = climate_fields.biome_ids.shape
        rows = np.clip(((90.0 - np.degrees(lat_rad)) * grid_height / 180.0).astype(int), 0, grid_height - 1)
        cols = (((np.degrees(lon_rad) + 180.0) * grid_width / 360.0).astype(int)) % grid_width
        colors = biomes.BIOME_COLORS[climate_fields.biome_ids[rows, cols]].astype(float)
        return colors
    except (AttributeError, ValueError):
        lat, _ = geometry.xyz_to_latlon(points)
        cold = np.clip((np.abs(np.degrees(lat)) - 35.0) / 50.0, 0.0, 1.0)
        green = np.clip(1.0 - np.maximum(elevations, 0.0) / 5000.0, 0.35, 1.0)
        return np.stack((75 + 130 * cold, 112 + 70 * cold, 58 + 75 * cold), axis=1) * green[:, None]


def _geotherm_temperature(hc_m: np.ndarray, hm_m: np.ndarray, lag_c: np.ndarray, depth_m: np.ndarray) -> np.ndarray:
    """Illustrative conductive crust geotherm from the model's Moho thermal proxy."""
    hc = np.maximum(hc_m, 0.0)
    hm = np.maximum(hm_m, 0.0)
    k = orogeny.CRUST_CONDUCTIVITY_W_M_K
    heat = orogeny.CRUST_HEAT_PRODUCTION_W_M3
    total = np.maximum(hc + hm, 1.0)
    q0 = (k * (orogeny.LITHOSPHERE_BASE_TEMPERATURE_C - orogeny.SURFACE_TEMPERATURE_C)
          + heat * hc * (0.5 * hc + hm)) / total
    depth = np.minimum(np.maximum(depth_m, 0.0), hc)
    temp = orogeny.SURFACE_TEMPERATURE_C + q0 * depth / k - heat * depth**2 / (2.0 * k)
    moho = orogeny.moho_temperature_c(hc, hm, lag_c)
    return np.where(depth_m <= hc, temp, moho)


def _mix(color: tuple[int, int, int] | np.ndarray, other: tuple[int, int, int], amount: float) -> tuple[int, int, int]:
    return tuple(np.clip(np.asarray(color) * (1.0 - amount) + np.asarray(other) * amount, 0, 255).astype(int))


def _parcel_markers(world, sample_xyz: np.ndarray, spacing_km: float) -> list[tuple[int, float]]:
    parcels = getattr(world, "pending_magma_parcels", ()) or ()
    if not parcels:
        return []
    origins = np.asarray([parcel.origin_xyz for parcel in parcels], dtype=float)
    dots = np.clip(origins @ sample_xyz.T, -1.0, 1.0)
    nearest = np.argmax(dots, axis=1)
    distance_km = np.arccos(dots[np.arange(len(origins)), nearest]) * PLANET_RADIUS_KM
    corridor_km = max(spacing_km * 1.5, 1.0)
    volumes = np.asarray([max(float(parcel.volume_m3), 0.0) for parcel in parcels])
    max_volume = max(float(volumes.max(initial=0.0)), 1.0)
    return [(int(i), float(np.clip(v / max_volume, 0.0, 1.0))) for i, v, d in zip(nearest, volumes, distance_km) if d <= corridor_km]


def render_png_base64(world, plane_angle_deg: float, position_index: int, width: int, height: int) -> tuple[str, dict]:
    """Render the selected arc from current surface columns and return PNG plus metadata."""
    spacing_km = line_spacing_rad(world.node_density) * PLANET_RADIUS_KM
    sample_xyz, tangent, distances = _great_circle_samples(plane_angle_deg, position_index, spacing_km)
    sampled, continental, _ = _sample_fields(world, sample_xyz)
    elev = sampled["elevation"].astype(float)
    hc = np.maximum(sampled["crustal_thickness_m"].astype(float), 0.0)
    hm = np.maximum(sampled["mantle_lithosphere_thickness_m"].astype(float), 0.0)
    glacier = np.maximum(sampled["glacier_depth"].astype(float), 0.0)
    lag = sampled["moho_thermal_lag_c"].astype(float)
    volcano = sampled["is_volcano"].astype(bool)
    active_vent = volcano & (sampled["volcano_active_years_remaining"] > 0.0)
    sea_level = float(world.sea_level_m)
    land_colors = _biome_colors(world, sample_xyz, elev)
    moho_temp = orogeny.moho_temperature_c(hc, hm, lag)
    mid_crust_temp = _geotherm_temperature(hc, hm, lag, hc * 0.5)

    margin_x = max(16, round(width * 0.025))
    margin_top, margin_bottom = max(30, round(height * 0.075)), max(45, round(height * 0.11))
    # Match the step-1 schematic's geometry: the 3,000 km arc occupies 690/820 of
    # the frame, and radial depth is exaggerated twice using the same pixel/km scale.
    # This keeps server frames comparable to the placeholder at any rendered size (#360).
    arc_width = (width - 2 * margin_x) * (690.0 / 820.0)
    half_arc_angle = ARC_LENGTH_KM / (2.0 * PLANET_RADIUS_KM)
    horizontal_scale = arc_width / (2.0 * PLANET_RADIUS_KM * np.sin(half_arc_angle))
    vertical_scale = VERTICAL_EXAGGERATION * horizontal_scale
    surface_center_y = margin_top + round(height * 0.08)
    arc_angles = distances / PLANET_RADIUS_KM
    x_samples = width / 2.0 + PLANET_RADIUS_KM * np.sin(arc_angles) * horizontal_scale
    surface_y = (
        surface_center_y
        + VERTICAL_EXAGGERATION * PLANET_RADIUS_KM * (1.0 - np.cos(arc_angles)) * horizontal_scale
        - (elev - sea_level) / 1000.0 * vertical_scale
    )
    x_pixel = np.arange(width)
    nearest = np.clip(np.searchsorted(x_samples, x_pixel, side="right") - 1, 0, len(x_samples) - 1)

    image = Image.new("RGB", (width, height), (14, 20, 34))
    draw = ImageDraw.Draw(image)
    # Draw the schematic mantle only inside the curved section window; filling the whole
    # plotting rectangle made the frame read as an asthenosphere block (#360).
    points_for_parcels = sample_xyz
    marker_by_sample = dict(_parcel_markers(world, points_for_parcels, spacing_km))
    melt_eligible = continental & (moho_temp >= orogeny.MELT_ONSET_MOHO_C)
    line_samples = np.linspace(0, len(sample_xyz) - 1, min(24, len(sample_xyz)), dtype=int)
    flow = mantle.flow_at(sample_xyz[line_samples], world.mantle_centers)
    flow_along = np.sum(flow * tangent[line_samples], axis=1)
    flow_max = max(float(np.max(np.abs(flow_along), initial=0.0)), 1e-15)

    for x in range(
        max(margin_x, int(x_samples[0])),
        min(width - margin_x, int(np.ceil(x_samples[-1])) + 1),
    ):
        i = int(nearest[x])
        x0, x1 = x, min(x + 1, width - margin_x)
        top = float(surface_y[i])
        local_depth_scale = vertical_scale * np.cos(arc_angles[i])
        depth_scale = local_depth_scale / 1000.0
        bottom = top + SECTION_DEPTH_KM * local_depth_scale
        seabed = float(elev[i])
        sea_y = (
            surface_center_y
            + VERTICAL_EXAGGERATION * PLANET_RADIUS_KM * (1.0 - np.cos(arc_angles[i])) * horizontal_scale
        )
        hcm = min(float(hc[i]) / 1000.0, SECTION_DEPTH_KM)
        hml = min((float(hc[i]) + float(hm[i])) / 1000.0, SECTION_DEPTH_KM)
        crust_end = top + hcm * local_depth_scale
        lith_end = top + hml * local_depth_scale
        base_color = tuple(np.clip(land_colors[i], 0, 255).astype(int)) if continental[i] else CRUST_OCEANIC
        # Ocean fill lies above the sampled seafloor, with shallow shelves visibly lighter.
        if seabed < sea_level:
            water_depth = max(sea_level - seabed, 0.0)
            deep = np.clip(water_depth / 6000.0, 0.0, 1.0)
            water_color = _mix(SEA_LEVEL_COLOR, (8, 24, 54), float(deep) * 0.82)
            draw.rectangle((x0, max(margin_top, sea_y), x1, min(height - margin_bottom, top)), fill=water_color)
        else:
            # Subtle local slope shading makes ridges and valleys read at a glance.
            before, after = max(i - 1, 0), min(i + 1, len(elev) - 1)
            slope = abs(float(elev[after] - elev[before])) / max(abs(float(distances[after] - distances[before])) * 1000.0, 1.0)
            base_color = _mix(base_color, (212, 216, 200) if slope > 0.08 else (22, 30, 28), 0.14)
            draw.rectangle((x0, max(margin_top, round(top - glacier[i] * depth_scale)), x1, round(top)), fill=(224, 237, 242))
        draw.rectangle((x0, round(top), x1, min(round(crust_end), height - margin_bottom)), fill=base_color)
        # Lower-crust thermal shading is intentionally subtle and is not itself melt.
        if hcm > 0.0:
            lower_start = round(top + hcm * local_depth_scale * 0.5)
            lower_end = min(round(crust_end), height - margin_bottom)
            tint = float(np.clip((mid_crust_temp[i] - 250.0) / 700.0, 0.0, 1.0) * 0.16)
            if lower_end > lower_start:
                draw.rectangle((x0, lower_start, x1, lower_end), fill=_mix(base_color, (224, 93, 51), tint))
            moho_tint = float(np.clip((moho_temp[i] - 250.0) / 700.0, 0.0, 1.0) * 0.2)
            moho_band = max(round(top), lower_end - 2)
            draw.line((x0, moho_band, x1, moho_band), fill=_mix(base_color, (245, 152, 77), moho_tint))
        draw.rectangle((x0, round(crust_end), x1, min(round(lith_end), height - margin_bottom)), fill=MANTLE_LITHOSPHERE)
        draw.rectangle((x0, round(lith_end), x1, round(bottom)), fill=ASTHENOSPHERE)
        if melt_eligible[i] and hcm > 0:
            hatch_top = max(round(crust_end - 8.0 * local_depth_scale), round(top))
            hatch_bottom = min(round(crust_end), height - margin_bottom)
            if hatch_bottom > hatch_top and x % 7 == 0:
                draw.line((x0, hatch_top, x1, hatch_bottom), fill=(255, 210, 120), width=1)
        if i in marker_by_sample:
            strength = marker_by_sample[i]
            color = MAGMA_PALETTE[min(2, int(strength * 3))]
            draw.ellipse((x - 3, top - 5, x + 3, top + 1), fill=color, outline=(45, 24, 24))
        if volcano[i]:
            color = (237, 126, 60) if active_vent[i] else (183, 190, 193)
            draw.polygon(((x, top - 2), (x - 3, top - 8), (x + 3, top - 8)), fill=color, outline=(32, 29, 38))

    # Draw the surface and depth reference after fills so they stay legible.
    profile = [(int(x_samples[i]), int(surface_y[i])) for i in range(len(x_samples))]
    draw.line(profile, fill=(229, 231, 202), width=max(1, round(width / 900)))
    center_x = width // 2
    for depth_km in (100, 200, 300, 350):
        y = int(surface_y[len(surface_y) // 2] + depth_km * vertical_scale)
        draw.line((center_x - 5, y, center_x + 5, y), fill=(230, 231, 239), width=1)
        draw.text((center_x + 8, y - 6), f"{depth_km} km", fill=(222, 225, 235))
    for j, value in zip(line_samples, flow_along):
        if abs(value) < flow_max * 0.03:
            continue
        x = int(x_samples[j])
        y = int(surface_y[j] - 12)
        direction = 1 if value >= 0 else -1
        length = 9 + 13 * abs(float(value)) / flow_max
        end_x = int(x + direction * length)
        draw.line((x, y, end_x, y), fill=(71, 236, 210), width=2)
        draw.polygon(((end_x, y), (end_x - direction * 5, y - 3), (end_x - direction * 5, y + 3)), fill=(71, 236, 210))

    draw.text((margin_x, 8), "Surface-sampled lithosphere section · 3,000 km arc", fill=(235, 239, 246))
    draw.text((margin_x, 22), f"350 km below local surface · vertical relief {VERTICAL_EXAGGERATION:g}× · geotherm illustrative", fill=(168, 180, 197))
    draw.text((margin_x, 35), "Active vents retain eruption potential; dormant vents have none. Neither label means erupting now.", fill=(150, 164, 184))
    legend_y = height - max(32, round(height * 0.06))
    legend = (
        (SEA_LEVEL_COLOR, "Water"), ((224, 237, 242), "Stored glacier"),
        (CRUST_CONTINENTAL, "Continental crust"), (CRUST_OCEANIC, "Oceanic crust"),
        (MANTLE_LITHOSPHERE, "Mantle lithosphere"), (ASTHENOSPHERE, "Schematic asthenosphere"),
        (MAGMA_PALETTE[1], "Stored melt source"), ((255, 210, 120), "Melt-eligible estimate"),
        ((71, 236, 210), "Surface mantle flow"), ((237, 126, 60), "Active vent"),
        ((183, 190, 193), "Dormant vent"),
    )
    x = margin_x
    for color, label in legend:
        item_width = 12 + len(label) * 6 + 11
        if x + item_width > width - margin_x and x > margin_x:
            legend_y -= 14
            x = margin_x
        draw.rectangle((x, legend_y, x + 9, legend_y + 9), fill=color)
        draw.text((x + 12, legend_y - 1), label, fill=(204, 212, 225))
        x += item_width

    buffer = io.BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    revision = f"{id(world):x}:{world.steps_taken}:{float(world.elapsed_years).hex()}:{sea_level.hex()}:{float(world.solar_multiplier).hex()}"
    metadata = {
        "world_revision": revision,
        "elapsed_years": float(world.elapsed_years),
        "plane_angle_deg": float(plane_angle_deg),
        "position_index": int(position_index),
        "arc_length_km": ARC_LENGTH_KM,
        "section_depth_km": SECTION_DEPTH_KM,
        "vertical_exaggeration": VERTICAL_EXAGGERATION,
        "image_base64": base64.b64encode(buffer.getvalue()).decode("ascii"),
    }
    return metadata["image_base64"], metadata

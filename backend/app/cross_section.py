"""Surface-sampled shallow lithosphere sections (#360).

The section is a narrow cut through the existing surface-node columns. It deliberately does
not construct a volumetric mantle: only Hc/Hm columns, their illustrative thermal profile,
surface vents, and the schematic asthenosphere below the modeled lithosphere are drawn.
Surface-sampled mantle-flow directions are projected into that schematic asthenosphere.
"""

from __future__ import annotations

import base64
import io

import numpy as np
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree

from . import faults, geometry, mantle, orogeny, plates, plumes
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
CRUST_OCEANIC = (82, 103, 121)
MANTLE_LITHOSPHERE = (66, 90, 90)
ASTHENOSPHERE = (117, 108, 148)
MAGMA_PALETTE = ((201, 61, 61), (244, 122, 60), (255, 214, 90))
SEA_LEVEL_COLOR = (83, 159, 203)
# Plume markers reuse the magma hues (per UX request): the head band/ellipse
# shares the stored-melt orange and the tail glyph the weak-magma red. The legend's
# separate "Mantle plume" / "Plume tail" keys disambiguate them from vents.
PLUME_HEAD_FILL = (244, 122, 60)
PLUME_TAIL_COLOR = (201, 61, 61)
# Interfaces are drawn this many times larger and box-filtered down, which antialiases them
# without touching flat fills; text is added afterwards so it stays crisp (#702).
SUPERSAMPLE = 3
# Arc samples per drawn column: below 1 so every sample owns a column for its glyph.
SAMPLES_PER_COLUMN = 0.7
# Continuous fields blend over this many node spacings. Copying the nearest node drew every
# interface as a staircase one node wide (#702).
KERNEL_RADIUS_SPACINGS = 2.0
KERNEL_NEIGHBOURS = 16
CONTINUOUS_FIELDS = (
    "elevation", "crustal_thickness_m", "mantle_lithosphere_thickness_m", "glacier_depth", "moho_thermal_lag_c",
)
# Schematic: the model stores no conduit width, and PLUME_TAIL_RADIUS_KM is a surface
# footprint, which would draw a blob rather than a conduit (#702).
PLUME_TAIL_CONDUIT_RADIUS_KM = 60.0
PLUME_TAIL_FILL = (214, 92, 70)
PLUME_TAIL_RIM = (120, 40, 34)
GLYPH_OUTLINE = (45, 24, 24)


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
    """Sample columns at `points`: categorical fields from the nearest node, continuous
    ones blended, with `nearest_node` naming each sample's node.

    The Wendland C2 kernel has compact support, so a node fades out of the blend before
    it leaves the neighbourhood and the result stays continuous along the arc (#702).
    """
    node_xyz, fields, continental = _world_surface(world)
    if not len(node_xyz):
        raise ValueError("world has no active surface nodes")
    tree = plates.cached_node_position_tree(world, node_xyz)
    neighbour_ranks = list(range(1, min(KERNEL_NEIGHBOURS, len(node_xyz)) + 1))
    distances, neighbours = tree.query(points, k=neighbour_ranks)
    nearest = neighbours[:, 0]
    sampled = {name: values[nearest] for name, values in fields.items()}
    sampled["nearest_node"] = nearest
    radius = KERNEL_RADIUS_SPACINGS * line_spacing_rad(world.node_density)
    q = np.clip(distances / radius, 0.0, 1.0)
    weights = (1.0 - q) ** 4 * (4.0 * q + 1.0)
    total = weights.sum(axis=1)
    covered = total > 1e-9  # no node in reach: keep the nearest column
    for name in CONTINUOUS_FIELDS:
        blended = (weights * fields[name][neighbours].astype(float)).sum(axis=1) / np.maximum(total, 1e-9)
        sampled[name] = np.where(covered, blended, sampled[name].astype(float))
    return sampled, continental[nearest], node_xyz


def _biome_colors(world, points: np.ndarray, elevations: np.ndarray) -> np.ndarray:
    """Sample the current climate map's natural palette, falling back to latitude/elevation."""
    try:
        climate_fields = climate.compute_climate_cached(world)
        lat_rad, lon_rad = geometry.xyz_to_latlon(points)
        grid_height, grid_width = climate_fields.biome_ids.shape
        # Bilinear in display colour: biome ids are categorical, their colours are not (#702).
        v = np.clip((90.0 - np.degrees(lat_rad)) * grid_height / 180.0 - 0.5, 0.0, grid_height - 1.0)
        u = (np.degrees(lon_rad) + 180.0) * grid_width / 360.0 - 0.5
        rows, cols = v.astype(int), np.floor(u).astype(int)
        next_rows = np.minimum(rows + 1, grid_height - 1)
        row_t, col_t = (v - rows)[:, None], (u - cols)[:, None]
        grid = biomes.BIOME_COLORS[climate_fields.biome_ids].astype(float)
        west, east = cols % grid_width, (cols + 1) % grid_width
        upper = grid[rows, west] * (1.0 - col_t) + grid[rows, east] * col_t
        lower = grid[next_rows, west] * (1.0 - col_t) + grid[next_rows, east] * col_t
        return upper * (1.0 - row_t) + lower * row_t
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


def _plume_markers(world, sample_xyz: np.ndarray, spacing_km: float) -> list[tuple[int, str]]:
    """Live plume markers per section sample: (sample index, role).

    Roles are `"head"` (nearest sample to a head-stage plume), `"head_footprint"`
    (further samples inside the head's emplacement radius), and `"tail"` (nearest
    sample to a tail-stage plume). The head ellipse and tail glyphs require the
    plume centre within the same corridor rule as `_parcel_markers`, but the
    footprint band is drawn wherever the section cuts it -- a head just outside
    the corridor still paints its band, only without the ellipse. The deep tail
    conduit below the 350 km section is why the renderer leaves its body open
    at the section floor (RFC 0007, #702). Head and tail stages are disjoint by construction
    (`tail_start_years == head_end_years`), so a head-stage plume never shows
    a conduit.
    """
    live = [p for p in (getattr(world, "mantle_plumes", ()) or ()) if p is not None]
    if len(sample_xyz) == 0 or not live:
        return []
    elapsed = float(getattr(world, "elapsed_years", 0.0) or 0.0)
    positions = np.asarray([tuple(float(v) for v in p.position) for p in live], dtype=float)
    norms = np.linalg.norm(positions, axis=1)
    norms[norms <= 0.0] = 1.0
    dots = np.clip((positions / norms[:, None]) @ sample_xyz.T, -1.0, 1.0)
    nearest = np.argmax(dots, axis=1)
    distance_km = np.arccos(dots[np.arange(len(live)), nearest]) * PLANET_RADIUS_KM
    corridor_km = max(spacing_km * 1.5, 1.0)
    markers: list[tuple[int, str]] = []
    unit = positions / norms[:, None]
    for rank, (plume, index, dist) in enumerate(zip(live, nearest, distance_km)):
        if not np.isfinite(dist):
            continue
        head_live = (
            not plume.is_tail_only
            and plume.born_years <= elapsed < plume.head_end_years
        )
        tail_live = (
            plume.tail_flux_m3_per_yr > 0.0
            and plume.tail_start_years <= elapsed < plume.tail_end_years
        )
        if head_live:
            radius_m = plumes.head_kernel_radius_m(
                plumes.head_footprint_area_m2(plume.head_volume_m3)
            )
            # along-arc distance in km (radius_m is meters -- see the /1000).
            along_km = np.arccos(np.clip(sample_xyz @ unit[rank], -1.0, 1.0)) * PLANET_RADIUS_KM
            inside = np.flatnonzero(along_km <= radius_m / 1000.0)
            if dist > corridor_km:
                # Centre off-section: show where the plane cuts the footprint,
                # but no head ellipse.
                markers.extend((int(other), "head_footprint") for other in inside)
            else:
                for other in inside:
                    if int(other) != int(index):
                        markers.append((int(other), "head_footprint"))
                markers.append((int(index), "head"))
        elif tail_live and dist <= corridor_km:
            markers.append((int(index), "tail"))
    return markers


# Marker priority when plumes share a sample: heads read first, then tails,
# then footprint extent.
_PLUME_ROLE_PRIORITY = {"head": 0, "tail": 1, "head_footprint": 2}


def _merge_plume_markers(world, sample_xyz: np.ndarray, spacing_km: float) -> dict[int, str]:
    """Merge plume markers with head > tail > footprint priority, so a later
    plume's footprint never overwrites an earlier plume's head/tail marker
    at a shared sample."""
    merged: dict[int, str] = {}
    for sample_index, role in _plume_markers(world, sample_xyz, spacing_km):
        incumbent = merged.get(sample_index)
        if incumbent is None or _PLUME_ROLE_PRIORITY[role] < _PLUME_ROLE_PRIORITY[incumbent]:
            merged[sample_index] = role
    return merged


def _parcel_markers(world, sample_xyz: np.ndarray, spacing_km: float) -> list[tuple[int, float]]:
    parcels = getattr(world, "pending_magma_parcels", ()) or ()
    if not parcels:
        return []
    origins = np.asarray([parcel.origin_xyz for parcel in parcels], dtype=float)
    # A bounded tree query, not a parcels-by-samples dot matrix: a mature world holds hundreds
    # of thousands of parcels, which cost seconds per frame once the arc was sampled per drawn
    # column (#702).
    corridor_km = max(spacing_km * 1.5, 1.0)
    reach = 2.0 * np.sin(corridor_km / (2.0 * PLANET_RADIUS_KM))
    chords, nearest = cKDTree(sample_xyz).query(origins, k=1, distance_upper_bound=reach)
    volumes = np.maximum(np.fromiter((parcel.volume_m3 for parcel in parcels), dtype=float, count=len(parcels)), 0.0)
    strengths = np.clip(volumes / max(float(volumes.max(initial=0.0)), 1.0), 0.0, 1.0)
    return [(int(nearest[k]), float(strengths[k])) for k in np.flatnonzero(np.isfinite(chords))]


def _boundary_sample_indices(world, sample_xyz: np.ndarray, node_xyz: np.ndarray) -> np.ndarray:
    world_plates = getattr(world, "plates", ())
    active_plates = [plate for plate in world_plates if plate.node_count()]
    if not active_plates:
        return np.zeros(0, dtype=int)
    plate_ids = np.concatenate([
        np.full(plate.node_count(), plate.plate_id, dtype=np.int64)
        for plate in active_plates
    ])
    _, nearest_nodes = plates.cached_node_position_tree(world, node_xyz).query(sample_xyz, k=1)
    sampled_plate_ids = plate_ids[nearest_nodes]
    return np.flatnonzero(sampled_plate_ids[1:] != sampled_plate_ids[:-1]) + 1


def _fault_sample_set(world, sample_xyz: np.ndarray, spacing_km: float) -> set[int]:
    by_id = {plate.plate_id: plate for plate in getattr(world, "plates", ())}
    sample_tree = cKDTree(sample_xyz)
    samples: set[int] = set()
    for fault in (*getattr(world, "faults", ()), *getattr(world, "boundary_faults", ())):
        plate = by_id.get(fault.plate_id)
        if plate is None:
            continue
        trace = getattr(fault, "world_polyline", None)
        if trace is None or len(trace) == 0:
            trace = faults.fault_world_points(fault, plate)
        if len(trace) == 0:
            continue
        distances_to_arc, indices = sample_tree.query(trace, k=1)
        close = np.flatnonzero(distances_to_arc <= max(spacing_km * 1.5 / PLANET_RADIUS_KM, 1e-6))
        samples.update(int(indices[index]) for index in close)
    return samples


def _vent_samples(volcano: np.ndarray, node_ids: np.ndarray) -> np.ndarray:
    """The middle sample of each run owned by a volcano node: one glyph per vent, even
    though the fine arc gives a node many samples (#702)."""
    run_edges = np.r_[0, np.flatnonzero(np.diff(node_ids)) + 1, len(node_ids)]
    middles = (run_edges[:-1] + run_edges[1:] - 1) // 2
    return middles[volcano[middles]]


def _draw_conduit_column(draw, x, x0, x1, lith_end_y, bottom_y, conduits, half_width, line) -> None:
    """Paint any tail conduit covering this column: domed under the lithosphere, and down
    to this column's own asthenosphere bottom, because the conduit continues below the
    section (RFC 0007) and a closed end would read as its base (#702)."""
    for conduit_x, conduit_top in conduits:
        offset = abs(x - conduit_x)
        if offset > half_width:
            continue
        dome = conduit_top + half_width - np.sqrt(half_width**2 - offset**2)
        start = max(round(dome), lith_end_y)
        if start >= bottom_y:
            continue
        fill = PLUME_TAIL_RIM if offset > half_width - line else PLUME_TAIL_FILL
        draw.rectangle((x0, start, x1, bottom_y), fill=fill)
        draw.rectangle((x0, start, x1, start + line), fill=PLUME_TAIL_RIM)


def _draw_sample_glyphs(draw, x, top, parcel_strength, plume_role, vent_color, line) -> None:
    """One sample's surface glyphs, sized in output pixels times `line`.

    These run after every column is filled, so a later column can't clip a glyph's right
    half, and once per sample, so one glyph can't smear across the sample's width.
    """
    if parcel_strength is not None:
        color = MAGMA_PALETTE[min(2, int(parcel_strength * 3))]
        draw.ellipse((x - 3 * line, top - 5 * line, x + 3 * line, top + line), fill=color, outline=GLYPH_OUTLINE, width=line)
    if plume_role == "head":
        draw.ellipse(
            (x - 4 * line, top - 9 * line, x + 4 * line, top - line),
            fill=PLUME_HEAD_FILL, outline=GLYPH_OUTLINE, width=line,
        )
    elif plume_role == "tail":
        draw.ellipse(
            (x - 3 * line, top - 7 * line, x + 3 * line, top - line),
            fill=PLUME_TAIL_COLOR, outline=GLYPH_OUTLINE, width=line,
        )
    if vent_color is not None:
        triangle = ((x, top - 2 * line), (x - 3 * line, top - 8 * line), (x + 3 * line, top - 8 * line))
        draw.polygon(triangle, fill=vent_color, outline=(32, 29, 38))


def _draw_flow_arrows(draw, x_samples, surface_y, cos_angles, line_samples, flow_along, depth_km, asthenosphere_km, vertical_scale, line) -> None:
    flow_max = max(float(np.max(np.abs(flow_along), initial=0.0)), 1e-15)
    for j, value, depth, asthenosphere in zip(line_samples, flow_along, depth_km, asthenosphere_km):
        if abs(value) < flow_max * 0.03 or asthenosphere <= 0.0:
            continue
        x = int(x_samples[j])
        y = int(surface_y[j] + depth * vertical_scale * cos_angles[j])
        direction = 1 if value >= 0 else -1
        end_x = int(x + direction * line * (9 + 13 * abs(float(value)) / flow_max))
        draw.line((x, y, end_x, y), fill=(71, 236, 210), width=2 * line)
        head_x = end_x - direction * 5 * line
        draw.polygon(((end_x, y), (head_x, y - 3 * line), (head_x, y + 3 * line)), fill=(71, 236, 210))


def render_png_base64(
    world, plane_angle_deg: float, position_index: int, width: int, height: int,
    *, show_boundaries: bool = False, show_faults: bool = False,
    show_plumes: bool = False,
) -> tuple[str, dict]:
    """Render the selected arc from current surface columns and return PNG plus metadata.

    Geometry is drawn at SUPERSAMPLE times the output size, so `line` is one output pixel
    in drawing units.
    """
    output_width, output_height = width, height
    line = SUPERSAMPLE
    width, height = output_width * line, output_height * line
    spacing_km = line_spacing_rad(world.node_density) * PLANET_RADIUS_KM
    sample_xyz, tangent, distances = _great_circle_samples(
        plane_angle_deg, position_index, 2.0 * ARC_LENGTH_KM / (width * SAMPLES_PER_COLUMN)
    )
    sampled, continental, node_xyz = _sample_fields(world, sample_xyz)
    boundary_samples = (
        set(_boundary_sample_indices(world, sample_xyz, node_xyz).tolist()) if show_boundaries else set()
    )
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

    margin_x = line * max(16, round(output_width * 0.025))
    margin_top = line * max(30, round(output_height * 0.075))
    margin_bottom = line * max(45, round(output_height * 0.11))
    # Match the step-1 schematic's geometry: the 3,000 km arc occupies 690/820 of
    # the frame, and radial depth is exaggerated twice using the same pixel/km scale.
    # This keeps server frames comparable to the placeholder at any rendered size (#360).
    arc_width = (width - 2 * margin_x) * (690.0 / 820.0)
    half_arc_angle = ARC_LENGTH_KM / (2.0 * PLANET_RADIUS_KM)
    horizontal_scale = arc_width / (2.0 * PLANET_RADIUS_KM * np.sin(half_arc_angle))
    vertical_scale = VERTICAL_EXAGGERATION * horizontal_scale
    surface_center_y = margin_top + line * round(output_height * 0.08)
    arc_angles = distances / PLANET_RADIUS_KM
    cos_angles = np.cos(arc_angles)
    x_samples = width / 2.0 + PLANET_RADIUS_KM * np.sin(arc_angles) * horizontal_scale
    surface_y = (
        surface_center_y
        + VERTICAL_EXAGGERATION * PLANET_RADIUS_KM * (1.0 - cos_angles) * horizontal_scale
        - (elev - sea_level) / 1000.0 * vertical_scale
    )
    x_pixel = np.arange(width)
    nearest = np.clip(np.searchsorted(x_samples, x_pixel, side="right") - 1, 0, len(x_samples) - 1)
    # First pixel column owned by each sample: searchsorted(..., side="right")
    # hands column x to the sample at or left of x, so ceil(x_samples[i]) is
    # always visited and always owned by i (round() could land on i-1's column
    # when the fractional part exceeds 0.5, silently dropping that glyph).
    owner = np.ceil(x_samples).astype(int)

    image = Image.new("RGB", (width, height), (14, 20, 34))
    draw = ImageDraw.Draw(image)
    # Draw the schematic mantle only inside the curved section window; filling the whole
    # plotting rectangle made the frame read as an asthenosphere block (#360).
    marker_by_sample = dict(_parcel_markers(world, sample_xyz, spacing_km))
    plume_by_sample = _merge_plume_markers(world, sample_xyz, spacing_km) if show_plumes else {}
    fault_samples = _fault_sample_set(world, sample_xyz, spacing_km) if show_faults else set()
    melt_eligible = continental & (moho_temp >= orogeny.MELT_ONSET_MOHO_C)
    line_samples = np.linspace(0, len(sample_xyz) - 1, min(24, len(sample_xyz)), dtype=int)
    flow = mantle.flow_at(sample_xyz[line_samples], world.mantle_centers)
    flow_along = np.sum(flow * tangent[line_samples], axis=1)
    lithosphere_base_km = (hc[line_samples] + hm[line_samples]) / 1000.0
    asthenosphere_thickness_km = np.maximum(SECTION_DEPTH_KM - lithosphere_base_km, 0.0)
    # The model's mantle-flow field is tangential at the surface, not depth-resolved. Place
    # its direction markers midway through the schematic asthenosphere (#360), rather than
    # above the surface where they looked like surface-plate arrows.
    flow_depth_km = lithosphere_base_km + asthenosphere_thickness_km * 0.5
    conduit_half_width = PLUME_TAIL_CONDUIT_RADIUS_KM * horizontal_scale
    conduits = [
        (
            float(x_samples[i]),
            float(surface_y[i] + min((hc[i] + hm[i]) / 1000.0, SECTION_DEPTH_KM) * vertical_scale * cos_angles[i]),
        )
        for i, role in plume_by_sample.items() if role == "tail"
    ]
    first_column = max(margin_x, int(x_samples[0]))
    end_column = min(width - margin_x, int(np.ceil(x_samples[-1])) + 1)

    for x in range(first_column, end_column):
        i = int(nearest[x])
        x0, x1 = x, min(x + 1, width - margin_x)
        top = float(surface_y[i])
        local_depth_scale = vertical_scale * cos_angles[i]
        depth_scale = local_depth_scale / 1000.0
        bottom = top + SECTION_DEPTH_KM * local_depth_scale
        seabed = float(elev[i])
        sea_y = (
            surface_center_y
            + VERTICAL_EXAGGERATION * PLANET_RADIUS_KM * (1.0 - cos_angles[i]) * horizontal_scale
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
            moho_band = max(round(top), lower_end - 2 * line)
            draw.rectangle((x0, moho_band, x1, moho_band + line), fill=_mix(base_color, (245, 152, 77), moho_tint))
        draw.rectangle((x0, round(crust_end), x1, min(round(lith_end), height - margin_bottom)), fill=MANTLE_LITHOSPHERE)
        draw.rectangle((x0, round(lith_end), x1, round(bottom)), fill=ASTHENOSPHERE)
        _draw_conduit_column(draw, x, x0, x1, round(lith_end), round(bottom), conduits, conduit_half_width, line)
        if i in boundary_samples:
            draw.line(
                (x, round(top), x, min(round(lith_end), height - margin_bottom)),
                fill=(88, 190, 255), width=2 * line,
            )
        if i in fault_samples:
            draw.line(
                (x, round(top), x, min(round(top + 16 * local_depth_scale), height - margin_bottom)),
                fill=(255, 125, 87), width=2 * line,
            )
        if melt_eligible[i] and hcm > 0:
            hatch_top = max(round(crust_end - 8.0 * local_depth_scale), round(top))
            hatch_bottom = min(round(crust_end), height - margin_bottom)
            if hatch_bottom > hatch_top and x % (7 * line) == 0:
                draw.line((x0, hatch_top, x1, hatch_bottom), fill=(255, 210, 120), width=line)
        if plume_by_sample.get(i) in ("head", "head_footprint"):
            draw.rectangle((x0, round(top) - 4 * line, x1, round(top)), fill=PLUME_HEAD_FILL)

    vent_samples = set(_vent_samples(volcano, sampled.get("nearest_node", np.arange(len(volcano)))).tolist())
    for i in sorted({*marker_by_sample, *plume_by_sample, *vent_samples}):
        if first_column <= owner[i] < end_column:
            vent_color = None
            if i in vent_samples:
                vent_color = (237, 126, 60) if active_vent[i] else (183, 190, 193)
            _draw_sample_glyphs(
                draw, int(owner[i]), float(surface_y[i]), marker_by_sample.get(i), plume_by_sample.get(i), vent_color, line,
            )

    # Draw the surface and depth reference after fills so they stay legible.
    profile = [(int(x_samples[i]), int(surface_y[i])) for i in range(len(x_samples))]
    draw.line(profile, fill=(229, 231, 202), width=line * max(1, round(output_width / 900)))
    center_x = width // 2
    depth_ticks = []
    for depth_km in (100, 200, 300, 350):
        y = int(surface_y[len(surface_y) // 2] + depth_km * vertical_scale)
        draw.line((center_x - 5 * line, y, center_x + 5 * line, y), fill=(230, 231, 239), width=line)
        depth_ticks.append((depth_km, y // line))
    _draw_flow_arrows(
        draw, x_samples, surface_y, cos_angles, line_samples, flow_along, flow_depth_km,
        asthenosphere_thickness_km, vertical_scale, line,
    )

    image = image.resize((output_width, output_height), Image.Resampling.BOX)
    draw = ImageDraw.Draw(image)
    for depth_km, y in depth_ticks:
        draw.text((output_width // 2 + 8, y - 6), f"{depth_km} km", fill=(222, 225, 235))
    text_x = margin_x // line
    draw.text((text_x, 8), "Surface-sampled lithosphere section · 3,000 km arc", fill=(235, 239, 246))
    draw.text((text_x, 22), f"350 km below local surface · vertical relief {VERTICAL_EXAGGERATION:g}× · geotherm illustrative", fill=(168, 180, 197))
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

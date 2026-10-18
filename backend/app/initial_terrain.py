"""Generation-time relief and drainage seeded from the world's initial tectonic state.

The first plate-motion prediction is deliberately read-only: it gives generation a plausible
convergent-boundary map without advancing the tectonic clock. The resulting terrain then gets
one short climate/hydrology bootstrap so saved worlds already have a routed river graph.
"""

from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

from . import geometry, lithosphere, rheology, torque
from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from .hydroclimate import biomes, climate, hydrology
from .plates import PLANET_RADIUS_KM, collect_all_elevation, gather_node_positions
from .hydroclimate.erosion import climate_grid_indices

if TYPE_CHECKING:
    from .world import World

FIRST_MOTION_STEP_YEARS = 5_000_000.0
INITIAL_GLACIER_DEPTH_M = 150.0
INITIAL_GLACIER_MIN_ELEVATION_M = 2_500.0
INITIAL_GLACIER_MAX_TEMPERATURE_C = -3.0
OROGEN_SEARCH_RADIUS_KM = 250.0
OROGEN_SECTION_LENGTH_KM = 1_000.0
OROGEN_MAX_RELIEF_M = 3_200.0
OROGEN_WIDTH_KM = 65.0
OROGEN_OFFSET_SIGMA_KM = 90.0
OROGEN_MAX_OFFSET_KM = 200.0
OROGEN_CRATON_SCALE_M = 2_000.0
OROGEN_MIN_BOUNDARY_SAMPLES = 6
OROGEN_MAX_STRAIGHTNESS_RESIDUAL_KM = 100.0
CHANNEL_LINEAR_LIMIT_KM = 400.0
CHANNEL_TURN_LIMIT_DEG = 10.0


def apply_predicted_orogens(world: "World") -> None:
    """Seed narrow, high ridges on weak land beside initially converging contacts."""
    spacing = line_spacing_rad(world.node_density)
    reach = OROGEN_SEARCH_RADIUS_KM / PLANET_RADIUS_KM
    predicted = {
        plate.plate_id: torque.predict_plate_motion(plate, world, world.plates, FIRST_MOTION_STEP_YEARS)
        for plate in world.plates
    }
    relief_by_plate = {
        plate.plate_id: np.zeros(len(plate.cell_keys), dtype=float) for plate in world.plates
    }
    plate_by_id = {plate.plate_id: plate for plate in world.plates}
    for plate in world.plates:
        points, elevation = plate.all_points_and_elevation()
        if not len(points):
            continue
        neighbours = plate.get_neighbours(world.plates, threshold_rad=reach)
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach)
        if not neighbours:
            continue
        neighbour_omega = np.zeros_like(inputs.neighbor_omega)
        for neighbour in neighbours:
            mask = inputs.neighbor_plate_id == neighbour.plate_id
            if np.any(mask):
                neighbour_omega[mask] = predicted[neighbour.plate_id]
        closing = rheology.normal_closing_rate_m_per_s(
            predicted[plate.plate_id], neighbour_omega, inputs.own_points, inputs.direction_to_neighbor
        )
        continental = effective_is_continental_from_codes(
            inputs.own_crust_type_codes, plate.crust_type == "continental"
        )
        own_craton = plate.collect("craton_crust_m")
        candidates = (
            (inputs.dist_to_neighbor <= reach)
            & (closing > torque.BOUNDARY_TRANSFORM_RATE_M_PER_S)
            & continental
            & (elevation > world.sea_level_m)
            & (inputs.neighbor_plate_id >= 0)
            & (inputs.neighbor_node_index >= 0)
        )
        for neighbour_id in np.unique(inputs.neighbor_plate_id[candidates]):
            # Visit each contact once; the lower id supplies the boundary samples.
            if plate.plate_id >= int(neighbour_id):
                continue
            neighbour = plate_by_id[int(neighbour_id)]
            neighbour_points, neighbour_elevation = neighbour.all_points_and_elevation()
            neighbour_craton = neighbour.collect("craton_crust_m")
            neighbour_continental = effective_is_continental_from_codes(
                neighbour.collect("crust_type_code"), neighbour.crust_type == "continental"
            )
            pair_mask = candidates & (inputs.neighbor_plate_id == neighbour_id)
            own_indices = np.flatnonzero(pair_mask)
            other_indices = inputs.neighbor_node_index[own_indices]
            land_contact = (
                (neighbour_elevation[other_indices] > world.sea_level_m)
                & neighbour_continental[other_indices]
            )
            own_indices = own_indices[land_contact]
            other_indices = other_indices[land_contact]
            if len(own_indices) < OROGEN_MIN_BOUNDARY_SAMPLES:
                continue

            boundary_xyz = geometry.normalize(points[own_indices] + neighbour_points[other_indices])
            mean_craton_own = float(np.mean(own_craton[own_indices]))
            mean_craton_other = float(np.mean(neighbour_craton[other_indices]))
            use_own = mean_craton_own <= mean_craton_other
            weak_plate = plate if use_own else neighbour
            weak_points = points[own_indices] if use_own else neighbour_points[other_indices]
            weak_craton_m = mean_craton_own if use_own else mean_craton_other
            seed_key = (world.seed, plate.plate_id, neighbour.plate_id, 0x4F524F)
            rng = np.random.default_rng(seed_key)
            ridge = _initial_ridge(
                boundary_xyz, weak_points, rng, closing[own_indices], weak_craton_m,
            )
            if ridge is None:
                continue
            origin_xyz, center_xy, tangent, normal, half_length_km, peak_m = ridge
            target_points, target_elevation = weak_plate.all_points_and_elevation()
            xy = _local_xy(origin_xyz, target_points) - center_xy
            along = xy @ tangent
            across = xy @ normal
            beyond = np.maximum(np.abs(along) - half_length_km, 0.0)
            distance = np.hypot(across, beyond)
            land = target_elevation > world.sea_level_m
            gain = peak_m * np.exp(-0.5 * (distance / OROGEN_WIDTH_KM) ** 2)
            gain *= np.clip((half_length_km + OROGEN_WIDTH_KM - np.abs(along)) / OROGEN_WIDTH_KM, 0.0, 1.0)
            gain = np.where(land, gain, 0.0)
            relief_by_plate[weak_plate.plate_id] = np.maximum(relief_by_plate[weak_plate.plate_id], gain)

    for plate in world.plates:
        gain = relief_by_plate[plate.plate_id]
        if not np.any(gain > 0.0):
            continue
        points, elevation = plate.all_points_and_elevation()
        hc, elevation = lithosphere.back_elevation_gain_fields(
            elevation,
            plate.collect("crustal_thickness_m"),
            plate.collect("mantle_lithosphere_thickness_m"),
            plate.collect("crust_type_code"),
            plate.crust_type,
            gain,
            gain > 0.0,
        )
        plate.set_fields_on_plate(crustal_thickness_m=hc, elevation=elevation)


def _local_xy(origin_xyz: np.ndarray, points_xyz: np.ndarray) -> np.ndarray:
    east, north = geometry.local_tangent_basis(origin_xyz)
    return geometry.azimuthal_equidistant_forward(origin_xyz, east, north, points_xyz) * PLANET_RADIUS_KM


def _initial_ridge(boundary_xyz, weak_points, rng, closing_rates, weak_craton_m):
    """Fit a straight contact section, then offset and rotate its narrow ridge trace."""
    origin_xyz = geometry.normalize(boundary_xyz.mean(axis=0))
    boundary_xy = _local_xy(origin_xyz, boundary_xyz)
    weak_xy = _local_xy(origin_xyz, weak_points)
    centered = boundary_xy - boundary_xy.mean(axis=0)
    _, _, axes = np.linalg.svd(centered, full_matrices=False)
    principal = axes[0]
    projection = boundary_xy @ principal
    section_length = min(OROGEN_SECTION_LENGTH_KM, float(np.ptp(projection)))
    if section_length <= 0.0:
        return None
    best = None
    starts = np.linspace(float(projection.min()), float(projection.max() - section_length), 7)
    for start in starts:
        selected = (projection >= start) & (projection <= start + section_length)
        if np.count_nonzero(selected) < OROGEN_MIN_BOUNDARY_SAMPLES:
            continue
        section = boundary_xy[selected]
        center = section.mean(axis=0)
        _, _, local_axes = np.linalg.svd(section - center, full_matrices=False)
        tangent = local_axes[0]
        normal = np.array([-tangent[1], tangent[0]])
        residual = float(np.sqrt(np.mean(((section - center) @ normal) ** 2)))
        score = residual - min(np.count_nonzero(selected), 30) * 0.25
        if best is None or score < best[0]:
            best = (score, center, tangent, normal, section)
    if best is None:
        return None
    residual, center, tangent, normal, section = best
    if residual > OROGEN_MAX_STRAIGHTNESS_RESIDUAL_KM:
        return None
    half_length = min(OROGEN_SECTION_LENGTH_KM / 2.0, float(np.ptp(section @ tangent)) / 2.0)
    if half_length < 100.0:
        return None
    weak_direction = weak_xy.mean(axis=0) - boundary_xy.mean(axis=0)
    if np.dot(weak_direction, normal) < 0.0:
        normal = -normal
    offset = min(abs(float(rng.normal(0.0, OROGEN_OFFSET_SIGMA_KM))), OROGEN_MAX_OFFSET_KM)
    center = center + normal * offset
    angle = np.radians(float(rng.uniform(2.0, 10.0)) * rng.choice([-1.0, 1.0]))
    cosine, sine = np.cos(angle), np.sin(angle)
    tangent = np.array([cosine * tangent[0] - sine * tangent[1], sine * tangent[0] + cosine * tangent[1]])
    normal = np.array([-tangent[1], tangent[0]])
    if np.dot(weak_direction, normal) < 0.0:
        normal = -normal
    speed = float(np.clip(np.mean(closing_rates) / (5.0 * torque.BOUNDARY_TRANSFORM_RATE_M_PER_S), 0.25, 1.0))
    craton_factor = float(np.exp(-max(0.0, weak_craton_m) / OROGEN_CRATON_SCALE_M))
    return origin_xyz, center, tangent, normal, half_length, OROGEN_MAX_RELIEF_M * speed * craton_factor


def _turn_angle_deg(previous: np.ndarray, current: np.ndarray, following: np.ndarray) -> float:
    incoming = geometry.normalize(current - np.dot(previous, current) * current)
    outgoing = geometry.normalize(following - np.dot(following, current) * current)
    return float(np.degrees(np.arccos(np.clip(np.dot(incoming, outgoing), -1.0, 1.0))))


def _reject_long_straight_channel_nodes(fields: hydrology.HydrologyFields) -> np.ndarray:
    """Return river nodes in any route that stays nearly straight for over 400 km."""
    river = fields.is_river.copy()
    flow = fields.flow_target
    upstream_count = np.zeros(len(river), dtype=np.int32)
    sources = np.flatnonzero(river & (flow >= 0))
    np.add.at(upstream_count, flow[sources], 1)
    heads = np.flatnonzero(river & (upstream_count == 0))
    rejected = np.zeros(len(river), dtype=bool)
    for head in heads:
        path = [int(head)]
        seen = {int(head)}
        node = int(head)
        while flow[node] >= 0:
            node = int(flow[node])
            if node in seen or not river[node]:
                break
            path.append(node)
            seen.add(node)
        if len(path) < 4:
            continue
        run_start = 0
        run_distance = 0.0
        for i in range(1, len(path) - 1):
            prev, cur, nxt = (fields.points[path[j]] for j in (i - 1, i, i + 1))
            angle = _turn_angle_deg(prev, cur, nxt)
            edge_km = float(geometry.angular_distance(cur, nxt) * PLANET_RADIUS_KM)
            if angle <= CHANNEL_TURN_LIMIT_DEG:
                run_distance += edge_km
                if run_distance > CHANNEL_LINEAR_LIMIT_KM:
                    rejected[path[run_start : i + 2]] = True
            else:
                run_start = i
                run_distance = 0.0
    return rejected


def initialize_river_channels(world: "World") -> None:
    """Run one zero-duration routing bootstrap, then seed only curved river channels."""
    node_cloud = gather_node_positions(world.plates)
    climate_fields = climate.compute_climate(
        world, *climate.grid_dimensions(world.climate_density), node_cloud=node_cloud
    )
    points, plates_in_order = node_cloud
    rows, cols = climate_grid_indices(points, *climate_fields.precipitation_mm.shape)
    precipitation = climate_fields.precipitation_mm[rows, cols]
    elevation = collect_all_elevation(plates_in_order)
    temperature = np.where(
        elevation <= world.sea_level_m,
        climate_fields.ocean_temperature_c[rows, cols],
        climate_fields.air_temperature_c[rows, cols],
    )
    initial_ice = (
        (elevation >= INITIAL_GLACIER_MIN_ELEVATION_M)
        & (temperature <= INITIAL_GLACIER_MAX_TEMPERATURE_C)
        & (elevation > world.sea_level_m)
    )
    offset = 0
    for plate in plates_in_order:
        count = len(plate.cell_keys)
        plate.set_fields_on_plate(
            glacier_depth=np.where(
                initial_ice[offset : offset + count], INITIAL_GLACIER_DEPTH_M,
                plate.collect("glacier_depth"),
            )
        )
        offset += count
    seasonal = biomes.seasonal_temp_amplitude(
        climate_fields.lat_deg[rows], biomes.grid_continentality(climate_fields.is_ocean)[rows, cols],
        world.axial_tilt_deg,
    )
    fields = hydrology.compute_hydrology(
        world, precipitation, temperature, 0.0, node_cloud=node_cloud,
        seasonal_amplitude_at_nodes=seasonal,
    )
    discharge = fields.discharge_m3_s if fields.discharge_m3_s is not None else np.zeros(len(points))
    rejected_straight = _reject_long_straight_channel_nodes(fields)
    fields.is_river[rejected_straight] = False
    world.climate_cache = climate_fields
    world.hydrology_cache = fields
    world.hydrology_cache_step = world.steps_taken
    eligible = fields.is_river
    if np.any(eligible):
        river_discharge = discharge[eligible]
        rank = np.zeros(len(points))
        order = np.argsort(river_discharge, kind="stable")
        rank[np.flatnonzero(eligible)[order]] = np.linspace(0.0, 1.0, len(order))
        depth = np.where(eligible, 20.0 + 100.0 * rank**1.7, 0.0)
        width = np.where(eligible, 80.0 + 1_200.0 * rank**1.4, 0.0)
        fields.channel_width = width
        offset = 0
        for plate in plates_in_order:
            count = len(plate.cell_keys)
            plate.set_fields_on_plate(
                channel_depth=depth[offset : offset + count],
                channel_width=width[offset : offset + count],
                glacier_depth=fields.glacier_depth[offset : offset + count],
                lake_depth=fields.lake_depth[offset : offset + count],
            )
            offset += count
    else:
        offset = 0
        for plate in plates_in_order:
            count = len(plate.cell_keys)
            plate.set_fields_on_plate(
                glacier_depth=fields.glacier_depth[offset : offset + count],
                lake_depth=fields.lake_depth[offset : offset + count],
            )
            offset += count

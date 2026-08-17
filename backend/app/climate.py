"""Climate: temperature, wind, ocean currents (with swells), humidity, and precipitation.

Unlike elevation (Lagrangian, carried by each plate's own rotating frame -- see
docs/simulation-model.md#why-not-a-grid) and unlike the render grid (`render_image.py`'s
`_render_grid_arrays`, a *ragged* lat/lon sweep with variable row width, immediately
flattened to 1D), climate lives on a third, genuinely fixed-shape equirectangular array grid
-- `lat: (H,)`, `lon: (W,)`, every field `(H, W)`. This is a deliberate choice: the specific
mechanisms below (Coriolis deflection, coastal/mountain redirection, divergence-based swell
detection, wind-driven humidity advection) all lean on `np.roll` wraparound, centered-
difference gradients, and land-excluding neighbor averaging, none of which work on a ragged
lattice or an irregular point cloud. This grid exists *only* here; it is never stored on
`World` and never touches `world.plates`.

**Fully stateless.** Every field here is recomputed from scratch on every call to
`compute_climate`, from whatever the *current* plate elevation/crust_type happens to be
(sampled via the same `cKDTree` nearest-neighbor technique `_render_grid_arrays` already
uses) -- mirroring the render grid's own "recompute from scratch every call" philosophy, and
matching how plate-sim itself documents climate as something that "re-derives almost
everything downstream of elevation from scratch" every step, unlike elevation itself. No new
per-step cost: `step_world` is completely unaware this module exists. Climate is only ever
computed when a climate-related view is actually being rendered.

**Ported from plate-sim** (`~/plate-sim`, a sibling project, `docs/simulation-model.md`),
adapted from its equirectangular-array grid (`grid.py`'s `latlon_axes`) to this module's own
fixed grid, and from plate-sim's own source (`wind.py`, `ocean_currents.py`, `climate.py`)
read directly for exact formulas: latitude-banded meridional wind + Coriolis zonal
deflection, mountain deflection/Venturi/wake, Ekman-based ocean currents + coastal
deflection/smoothing/wake + land swirl + circumglobal boost, convergence-based swell
detection, semi-Lagrangian temperature advection along currents, evaporation-ceiling +
wind-driven 2D humidity advection, and orographic precipitation.

**Deliberately not ported** (mantle-bloom has no vegetation, rivers, or lakes -- these
mechanisms' *inputs* don't exist here, not a reduced-fidelity choice): humidity's
evapotranspiration term, river outflow feeding currents, lake climate influence.
**Deliberately cut** (confirmed with the user): river outflow, deep currents, and
precipitation's zonal latitude-climatology baseline (equator/mid-latitude wet bands) --
precipitation here is purely a function of humidity and orographic lift.

**Pipeline order** (breaks what would otherwise be a circular dependency -- wind needs
temperature, but the *final* ocean temperature needs currents, which need wind): insolation
-> pre-advection land/ocean baseline temperatures -> wind (from the *baseline* combined
surface temperature) -> ocean currents (from wind) -> ocean swells (from final currents) ->
final ocean temperature (baseline advected along final currents) -> air temperature (baseline
moderated toward the *final* nearest-ocean temperature) -> humidity (evaporation from *final*
ocean temperature, advected by wind) -> precipitation (from humidity + wind-over-mountains).
See `compute_climate` for the concrete call order.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, plates
from .world import World

# ---------------------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------------------

# Starting resolution -- tunable/benchmarked the same way GRID_SPACING_KM was tuned for the
# render grid. 2 degrees per cell in each direction.
GRID_HEIGHT = 90
GRID_WIDTH = 180

# A fixed reference grid width, decoupled from GRID_WIDTH, purely so the fixed-*degree*
# offset distances below (mountain/coast wake lookback, mountain tangent sampling) stay
# physically meaningful if GRID_WIDTH is retuned later -- same pattern as plate-sim's own
# module-level `_REFERENCE_CELL_DEG`.
_REFERENCE_WIDTH = 180
_REFERENCE_CELL_DEG = 360.0 / _REFERENCE_WIDTH

_EPS = 1e-9


@dataclass
class ClimateFields:
    lat_deg: np.ndarray  # (H,)
    lon_deg: np.ndarray  # (W,)
    world_xyz: np.ndarray  # (H, W, 3) unit vectors, for projecting cells onto the map
    elevation_m: np.ndarray  # (H, W)
    is_ocean: np.ndarray  # (H, W) bool
    land_temperature_c: np.ndarray  # (H, W) -- solar heating (+ lapse rate), unmoderated
    ocean_temperature_c: np.ndarray  # (H, W) -- final, current-advected
    air_temperature_c: np.ndarray  # (H, W) -- final, moderated toward nearby ocean
    wind_u: np.ndarray  # (H, W) eastward
    wind_v: np.ndarray  # (H, W) northward
    current_u: np.ndarray  # (H, W) eastward, zero over land
    current_v: np.ndarray  # (H, W) northward, zero over land
    humidity: np.ndarray  # (H, W), roughly [0, MAX_EVAPORATION_CEILING]
    precipitation_mm: np.ndarray  # (H, W)
    swell_rows: np.ndarray  # (K,) int, sampled convergence points
    swell_cols: np.ndarray  # (K,) int


# ---------------------------------------------------------------------------------------
# Grid construction + elevation/crust sampling
# ---------------------------------------------------------------------------------------


def _build_grid(height: int, width: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Row 0 = north pole, row increases southward (matches plate-sim's own convention);
    column increases eastward, wraps. Returns (lat_deg (H,), lon_deg (W,), world_xyz (H,W,3))."""
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_grid = np.repeat(lon_deg[None, :], height, axis=0)
    world_xyz = geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))
    return lat_deg, lon_deg, world_xyz


def _sample_elevation_and_crust(world: World, world_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-elevation-node resample of the *current* plate state onto the climate grid --
    same cKDTree technique render_image.py's _render_grid_arrays already uses. Returns
    (elevation_m, is_ocean), both (H, W)."""
    height, width, _ = world_xyz.shape
    points_list, elev_list, crust_list = [], [], []
    for plate in world.plates:
        pts, elev = plate.all_points_and_elevation()
        if len(pts) == 0:
            continue
        points_list.append(pts)
        elev_list.append(elev)
        crust_list.append(np.full(len(pts), plate.crust_type == "oceanic"))
    flat_xyz = world_xyz.reshape(-1, 3)
    if not points_list:
        return np.zeros((height, width)), np.ones((height, width), dtype=bool)

    all_points = np.concatenate(points_list, axis=0)
    all_elev = np.concatenate(elev_list, axis=0)
    all_is_oceanic = np.concatenate(crust_list, axis=0)
    tree = cKDTree(all_points)
    _, idx = tree.query(flat_xyz)
    elevation = all_elev[idx].reshape(height, width)
    is_ocean = all_is_oceanic[idx].reshape(height, width)
    return elevation, is_ocean


# ---------------------------------------------------------------------------------------
# Small shared array helpers (np.roll wraparound/gradients, closed-form offset sampling)
# ---------------------------------------------------------------------------------------


def _smooth_field(field: np.ndarray, iterations: int) -> np.ndarray:
    """Plain 4-neighbor Jacobi box blur, longitude-wrapping, no masking."""
    smoothed = field
    for _ in range(iterations):
        smoothed = 0.25 * (
            np.roll(smoothed, -1, axis=1) + np.roll(smoothed, 1, axis=1)
            + np.roll(smoothed, -1, axis=0) + np.roll(smoothed, 1, axis=0)
        )
    return smoothed


def _centered_gradient(field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(eastward, northward) centered-difference gradient, longitude-wrapping. Row 0 = north
    pole, so "northward" is the *negative* row direction -- matches plate-sim's own sign
    convention (`np.roll(f, 1, axis=0) - np.roll(f, -1, axis=0)`)."""
    gx = (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / 2.0
    gy = (np.roll(field, 1, axis=0) - np.roll(field, -1, axis=0)) / 2.0
    return gx, gy


def _sample_at_offset(field: np.ndarray, dir_u: np.ndarray, dir_v: np.ndarray, dist_deg: float, lat_deg: np.ndarray) -> np.ndarray:
    """Samples `field` at each cell's own position offset by `dist_deg` along (dir_u, dir_v)
    (already unit-length), nearest-cell. Same closed-form lat/lon -> row/col technique used
    throughout plate-sim for offset-position sampling (mountain tangent picking, wake
    lookback, current advection)."""
    height, width = field.shape
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    cos_lat = np.clip(np.cos(np.radians(lat_grid)), 0.15, 1.0)

    dst_lat = lat_grid + dir_v * dist_deg
    dst_lon = lon_deg[None, :] + dir_u * dist_deg / cos_lat

    dst_row = np.clip(np.round((90.0 - dst_lat) / (180.0 / height) - 0.5).astype(np.int64), 0, height - 1)
    dst_col = np.round((dst_lon + 180.0) / (360.0 / width) - 0.5).astype(np.int64) % width
    return field[dst_row, dst_col]


def _rotate_90(x: np.ndarray, y: np.ndarray, sign) -> tuple[np.ndarray, np.ndarray]:
    """Rotate (x, y) by 90 degrees (positive `sign` = counterclockwise in (east, north))."""
    return -sign * y, sign * x


def _cancel_and_redirect(
    u: np.ndarray, v: np.ndarray, nx: np.ndarray, ny: np.ndarray,
    tangent_x: np.ndarray, tangent_y: np.ndarray, has_obstacle: np.ndarray, speedup_factor: float,
) -> tuple[np.ndarray, np.ndarray]:
    """The mechanism shared by both wind-around-mountains and current-around-coast: cancel
    whatever velocity component points *into* the obstacle (along the outward normal
    `(nx, ny)`) and redirect that same magnitude into the tangential direction
    `(tangent_x, tangent_y)`, scaled up by `speedup_factor` (a real funneling/nozzle effect
    -- the same water/air squeezed along a narrower effective path has to move faster).
    Callers differ only in how `(nx, ny)` and `(tangent_x, tangent_y)` are derived -- a fixed
    hemisphere sense for ocean coastal deflection (real boundary currents have a
    Coriolis-preferred circulation direction), local topology for mountain deflection
    (mesoscale flow-splitting doesn't)."""
    normal_component = u * nx + v * ny
    into_obstacle = np.clip(-normal_component, 0.0, None)
    redirected = into_obstacle * speedup_factor
    u2 = np.where(has_obstacle, u + into_obstacle * nx + redirected * tangent_x, u)
    v2 = np.where(has_obstacle, v + into_obstacle * ny + redirected * tangent_y, v)
    return u2, v2


def _weighted_sample_without_replacement(rng: np.random.Generator, weights: np.ndarray, k: int) -> np.ndarray:
    """Indices of up to `k` points sampled without replacement, probability proportional to
    `weights` (zero-weight points never selected). Same technique plate-sim's hazards.py
    uses for earthquake/volcano placement and ocean_currents.py for swells."""
    positive = weights > 0
    n_positive = int(np.sum(positive))
    if n_positive == 0:
        return np.array([], dtype=int)
    k = min(k, n_positive)
    idx = np.flatnonzero(positive)
    p = weights[idx].astype(float)
    p /= p.sum()
    return rng.choice(idx, size=k, replace=False, p=p)


# ---------------------------------------------------------------------------------------
# Insolation
# ---------------------------------------------------------------------------------------

SUNLIGHT = 1.0
INSOLATION_FLOOR = 0.03
AXIAL_TILT_DECLINATION_SAMPLES = 24


def compute_insolation(lat_deg: np.ndarray, axial_tilt_deg: float) -> np.ndarray:
    """Annual-mean insolation, (H,) broadcastable. Flat zenith-angle cosine law with no
    tilt; with tilt, the mean of that same law (clipped at 0, sun below horizon) over
    `AXIAL_TILT_DECLINATION_SAMPLES` declinations swept between -tilt and +tilt -- the
    sub-solar latitude's annual sweep, not an actual season cycle (this model has no
    calendar; one step spans thousands to millions of years)."""
    lat_r = np.radians(lat_deg)
    if axial_tilt_deg <= 1e-6:
        return np.clip(np.cos(lat_r), INSOLATION_FLOOR, 1.0) * SUNLIGHT
    declinations = np.radians(np.linspace(-axial_tilt_deg, axial_tilt_deg, AXIAL_TILT_DECLINATION_SAMPLES))
    cos_zenith = np.cos(lat_r[:, None] - declinations[None, :])
    return np.clip(np.clip(cos_zenith, 0.0, None).mean(axis=1), INSOLATION_FLOOR, 1.0) * SUNLIGHT


# ---------------------------------------------------------------------------------------
# Temperature: land (solar heating + lapse rate), ocean baseline, air (moderated)
# ---------------------------------------------------------------------------------------

LAND_TEMP_MIN_C = -60.0
LAND_TEMP_RANGE_C = 95.0
LAPSE_RATE_C_PER_KM = 6.5
WATER_TEMP_MIN_C = -2.0
WATER_TEMP_RANGE_C = 32.0
# How far (degrees, great-circle) the ocean's moderating influence reaches onto land --
# an e-folding distance, not a hard cutoff.
MARITIME_INFLUENCE_DIST_DEG = 15.0


def compute_land_temperature(insolation_row: np.ndarray, elevation_m: np.ndarray) -> np.ndarray:
    """Solar heating only (the user's own description), plus elevation-based lapse-rate
    cooling -- kept as part of the *same* base-heating formula rather than a separate causal
    channel, matching plate-sim's own compute_temperature (mountains being cold is a
    consequence of solar heating at altitude, not an extra input)."""
    base = LAND_TEMP_MIN_C + LAND_TEMP_RANGE_C * insolation_row[:, None]
    altitude_cooling = LAPSE_RATE_C_PER_KM * np.clip(elevation_m, 0.0, None) / 1000.0
    return base - altitude_cooling


def compute_ocean_temperature_baseline(insolation_row: np.ndarray, height: int, width: int) -> np.ndarray:
    """Pre-advection zonal baseline -- narrower range and a freezing floor relative to land,
    water has far more thermal inertia (mirrors plate-sim's compute_water_temperature)."""
    base = WATER_TEMP_MIN_C + WATER_TEMP_RANGE_C * insolation_row[:, None]
    return np.repeat(base, width, axis=1) if base.shape[1] == 1 else base


def _nearest_ocean_gather(is_ocean: np.ndarray, world_xyz: np.ndarray, field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For every grid cell, the great-circle angular distance to the nearest ocean cell and
    that ocean cell's own `field` value -- a cKDTree chord-distance query over ocean-cell 3D
    positions, replacing plate-sim's lat/lon-tangent-plane BFS (`multi_source_bfs_nearest`)
    with an equivalent (true 3D distance, no pole/antimeridian special-casing needed)
    technique already idiomatic in this codebase. Returns (dist_rad, gathered_value), both
    (H, W); if there's no ocean at all, dist is +inf and gathered_value is 0 everywhere."""
    height, width = is_ocean.shape
    flat_ocean = is_ocean.reshape(-1)
    if not np.any(flat_ocean):
        return np.full((height, width), np.inf), np.zeros((height, width))

    ocean_xyz = world_xyz.reshape(-1, 3)[flat_ocean]
    ocean_values = field.reshape(-1)[flat_ocean]
    tree = cKDTree(ocean_xyz)
    chord_dist, idx = tree.query(world_xyz.reshape(-1, 3))
    # Chord distance -> great-circle angular distance (points are unit vectors).
    angular_dist = 2.0 * np.arcsin(np.clip(chord_dist / 2.0, 0.0, 1.0))
    return angular_dist.reshape(height, width), ocean_values[idx].reshape(height, width)


def compute_air_temperature(land_temperature_c: np.ndarray, ocean_temperature_c: np.ndarray, is_ocean: np.ndarray, world_xyz: np.ndarray) -> np.ndarray:
    """Land temperature's own solar-heating baseline, pulled toward the *nearest ocean
    cell's* (final, current-advected) temperature by a distance-based falloff -- "moderating
    effect of oceans," literally: over the ocean itself the influence is 1.0 (air right
    above the water is the water's own temperature), decaying inland with distance."""
    dist_rad, nearest_ocean_temp = _nearest_ocean_gather(is_ocean, world_xyz, ocean_temperature_c)
    dist_deg = np.degrees(dist_rad)
    influence = np.exp(-dist_deg / MARITIME_INFLUENCE_DIST_DEG)
    return land_temperature_c * (1.0 - influence) + nearest_ocean_temp * influence


# ---------------------------------------------------------------------------------------
# Wind
# ---------------------------------------------------------------------------------------

TRADE_WIND_MAX_LAT = 30.0
WESTERLIES_MAX_LAT = 60.0
MERIDIONAL_BASE_SPEED = 6.0
CORIOLIS_DEFLECTION_GAIN = 1.8
# Additive contribution from the real local gradient of the (pre-advection) surface
# temperature -- the user's explicit "wind is affected by temperature gradients" ask, layered
# on top of the latitude-banded structure above rather than replacing it (plate-sim's own
# docs record that deriving wind from a computed gradient *alone* produced no visible
# planetary-scale structure).
GRADIENT_WIND_COEFFICIENT = 0.4
ELEVATION_SLOWDOWN_REF_M = 4000.0
MIN_ELEVATION_SPEED_FACTOR = 0.4
MOUNTAIN_GRADIENT_SMOOTHING_ITERATIONS = 4
MOUNTAIN_OBSTACLE_ELEVATION_M = 2000.0
MOUNTAIN_OBSTACLE_RAMP_M = 800.0
MOUNTAIN_VENTURI_SPEEDUP_FACTOR = 1.4
MOUNTAIN_TANGENT_SAMPLE_STEPS = 3
MOUNTAIN_WAKE_LOOKBACK_STEPS = 5
MOUNTAIN_WAKE_MIN_SPEED_FACTOR = 0.5


def coriolis_parameter(lat_deg: np.ndarray) -> np.ndarray:
    return np.sin(np.radians(lat_deg))


def meridional_direction_for_lat(lat_deg: np.ndarray) -> np.ndarray:
    """+1 (northward) or -1 (southward): equatorward in the trade-wind and polar-easterly
    bands, poleward in the westerlies -- the near-surface branch of the real 3-cell
    circulation, kept as an empirical pattern (see module docstring)."""
    abs_lat = np.abs(lat_deg)
    hemisphere = np.where(lat_deg >= 0, 1.0, -1.0)
    toward_equator = -hemisphere
    toward_pole = hemisphere
    return np.where(abs_lat < TRADE_WIND_MAX_LAT, toward_equator, np.where(abs_lat < WESTERLIES_MAX_LAT, toward_pole, toward_equator))


def zonal_direction_for_lat(lat_deg: np.ndarray) -> np.ndarray:
    f = coriolis_parameter(lat_deg)
    meridional_dir = meridional_direction_for_lat(lat_deg)
    sign = np.sign(f * meridional_dir)
    sign = np.where(sign == 0, -1.0, sign)
    return sign


def _elevation_speed_factor(elevation_m: np.ndarray) -> np.ndarray:
    slowdown = np.clip(elevation_m / ELEVATION_SLOWDOWN_REF_M, 0.0, 1.0)
    return 1.0 - (1.0 - MIN_ELEVATION_SPEED_FACTOR) * slowdown


def _mountain_deflection(u: np.ndarray, v: np.ndarray, elevation_m: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cancels wind's into-slope component and redirects it tangentially (Venturi speedup),
    tangent side chosen by local topology -- mesoscale mountain flow-splitting has no
    Coriolis-preferred side the way basin-scale boundary currents do."""
    smoothed = _smooth_field(elevation_m, MOUNTAIN_GRADIENT_SMOOTHING_ITERATIONS)
    gx, gy = _centered_gradient(smoothed)
    magnitude = np.hypot(gx, gy)
    ux = np.divide(gx, magnitude, out=np.zeros_like(gx), where=magnitude > _EPS)
    uy = np.divide(gy, magnitude, out=np.zeros_like(gy), where=magnitude > _EPS)

    tangent1_x, tangent1_y = _rotate_90(ux, uy, 1.0)
    tangent2_x, tangent2_y = _rotate_90(ux, uy, -1.0)
    sample_dist = MOUNTAIN_TANGENT_SAMPLE_STEPS * _REFERENCE_CELL_DEG
    elev1 = _sample_at_offset(smoothed, tangent1_x, tangent1_y, sample_dist, lat_deg)
    elev2 = _sample_at_offset(smoothed, tangent2_x, tangent2_y, sample_dist, lat_deg)
    take1 = elev1 <= elev2
    tangent_x = np.where(take1, tangent1_x, tangent2_x)
    tangent_y = np.where(take1, tangent1_y, tangent2_y)

    ramp = np.clip(
        (smoothed - (MOUNTAIN_OBSTACLE_ELEVATION_M - MOUNTAIN_OBSTACLE_RAMP_M)) / (2.0 * MOUNTAIN_OBSTACLE_RAMP_M), 0.0, 1.0
    )
    has_gradient = magnitude > _EPS
    # _cancel_and_redirect expects an *outward*-from-obstacle normal (matching
    # _coastal_normal's "push away from land" convention) -- (ux, uy) points uphill, i.e.
    # into the slope, so the outward normal is its negation.
    u2, v2 = _cancel_and_redirect(u, v, -ux, -uy, tangent_x, tangent_y, has_gradient, MOUNTAIN_VENTURI_SPEEDUP_FACTOR)
    return u + ramp * (u2 - u), v + ramp * (v2 - v)


def _mountain_wake_factor(u: np.ndarray, v: np.ndarray, elevation_m: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
    """Lee-side speed factor behind mountainous terrain, walking backward along the wind's
    own (post-deflection) direction checking for upstream obstacle terrain."""
    height, width = elevation_m.shape
    is_obstacle = elevation_m > MOUNTAIN_OBSTACLE_ELEVATION_M
    speed = np.hypot(u, v)
    has_wind = speed > 1e-6
    dir_u = np.divide(u, speed, out=np.zeros_like(speed), where=has_wind)
    dir_v = np.divide(v, speed, out=np.zeros_like(speed), where=has_wind)

    nearest_step = np.full((height, width), MOUNTAIN_WAKE_LOOKBACK_STEPS + 1, dtype=np.int32)
    found = np.zeros((height, width), dtype=bool)
    for step in range(1, MOUNTAIN_WAKE_LOOKBACK_STEPS + 1):
        upstream_here = _sample_at_offset(is_obstacle, -dir_u, -dir_v, step * _REFERENCE_CELL_DEG, lat_deg)
        newly_found = upstream_here & ~found
        nearest_step[newly_found] = step
        found |= upstream_here

    in_wake = has_wind & (nearest_step <= MOUNTAIN_WAKE_LOOKBACK_STEPS)
    wake_strength = np.where(in_wake, 1.0 - (nearest_step - 1) / MOUNTAIN_WAKE_LOOKBACK_STEPS, 0.0)
    return 1.0 - wake_strength * (1.0 - MOUNTAIN_WAKE_MIN_SPEED_FACTOR)


def compute_wind(lat_deg: np.ndarray, elevation_m: np.ndarray, surface_temperature_c: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns (u, v, elevation_factor) -- elevation_factor folds together the own-cell
    terrain slowdown and the mountain wake, re-clipped to [MIN_ELEVATION_SPEED_FACTOR, 1.0],
    exposed separately from raw speed so humidity's retention formula can respond to
    terrain-driven slowdown specifically (see compute_humidity)."""
    height, width = elevation_m.shape
    meridional_dir = meridional_direction_for_lat(lat_deg)[:, None]
    f = coriolis_parameter(lat_deg)[:, None]
    own_elevation_factor = _elevation_speed_factor(elevation_m)

    v = meridional_dir * MERIDIONAL_BASE_SPEED * own_elevation_factor
    u = CORIOLIS_DEFLECTION_GAIN * f * v

    grad_x, grad_y = _centered_gradient(surface_temperature_c)
    # Surface wind flows toward warmer air (a stylized convection/pressure-gradient
    # intuition: warm air rises, drawing surface flow in to replace it) -- along the
    # temperature gradient itself, not against it.
    u = u + GRADIENT_WIND_COEFFICIENT * grad_x
    v = v + GRADIENT_WIND_COEFFICIENT * grad_y

    u, v = _mountain_deflection(u, v, elevation_m, lat_deg)
    wake_factor = _mountain_wake_factor(u, v, elevation_m, lat_deg)
    u, v = u * wake_factor, v * wake_factor

    elevation_factor = np.clip(own_elevation_factor * wake_factor, MIN_ELEVATION_SPEED_FACTOR, 1.0)
    return u, v, elevation_factor


# ---------------------------------------------------------------------------------------
# Ocean currents (+ swells)
# ---------------------------------------------------------------------------------------

CORIOLIS_DEFLECTION_DEG = 35.0
DEFLECTION_SPEEDUP_FACTOR = 1.35
COASTAL_SMOOTHING_ITERATIONS = 8
COASTAL_SMOOTHING_WEIGHT = 0.35
COASTAL_RESTORE_FRACTION = 0.12
SWIRL_PEAK_DIST_DEG = 22.0
SWIRL_DECAY_DIST_DEG = 18.0
SWIRL_MAX_SEARCH_DIST_DEG = SWIRL_PEAK_DIST_DEG + 4.0 * SWIRL_DECAY_DIST_DEG
SWIRL_SPEED_FRACTION = 0.7
CIRCUMGLOBAL_SPEEDUP_FACTOR = 1.6
WAKE_LOOKBACK_STEPS = 4
WAKE_MIN_SPEED_FACTOR = 0.45
WAKE_MIXING_ANGLE_MAX_DEG = 50.0
MAX_OCEAN_SWELLS = 25
OCEAN_SWELL_CONVERGENCE_REFERENCE = 2.5


def _ekman_current(wind_u: np.ndarray, wind_v: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hemisphere = np.where(lat_deg >= 0, 1.0, -1.0)[:, None]
    theta = np.radians(-CORIOLIS_DEFLECTION_DEG * hemisphere)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    u = wind_u * cos_t - wind_v * sin_t
    v = wind_u * sin_t + wind_v * cos_t
    return u, v


def _coastal_normal(is_ocean: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    is_land = ~is_ocean
    north_is_land = np.roll(is_land, 1, axis=0)
    south_is_land = np.roll(is_land, -1, axis=0)
    east_is_land = np.roll(is_land, -1, axis=1)
    west_is_land = np.roll(is_land, 1, axis=1)

    push_x = west_is_land.astype(float) - east_is_land.astype(float)
    push_y = south_is_land.astype(float) - north_is_land.astype(float)
    mag = np.hypot(push_x, push_y)
    has_coast = is_ocean & (mag > _EPS)
    safe_mag = np.where(mag > _EPS, mag, 1.0)
    return push_x / safe_mag, push_y / safe_mag, has_coast


def _deflect_into_coast(u: np.ndarray, v: np.ndarray, is_ocean: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nx, ny, has_coast = _coastal_normal(is_ocean)
    hemisphere = np.where(lat_deg >= 0, 1.0, -1.0)[:, None]
    tangent_x, tangent_y = _rotate_90(nx, ny, hemisphere)
    return _cancel_and_redirect(u, v, nx, ny, tangent_x, tangent_y, has_coast, DEFLECTION_SPEEDUP_FACTOR)


def _smooth_along_coast(u: np.ndarray, v: np.ndarray, is_ocean: np.ndarray, lat_deg: np.ndarray, ambient_u: np.ndarray, ambient_v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Jacobi neighbor-average (ocean-only), re-deflected each pass, leaking a small
    fraction back toward the ambient (pre-deflection) value each pass so the deflection's
    reach along a coast is real but bounded, not an unbounded diffusion of the whole basin."""
    ocean_f = is_ocean.astype(float)
    safe_count = np.where(ocean_f > 0, ocean_f, 1.0)

    def neighbor_average(field: np.ndarray) -> np.ndarray:
        masked = field * ocean_f
        total = (
            np.roll(masked, -1, axis=1) + np.roll(masked, 1, axis=1)
            + np.roll(masked, -1, axis=0) + np.roll(masked, 1, axis=0)
        )
        count = (
            np.roll(ocean_f, -1, axis=1) + np.roll(ocean_f, 1, axis=1)
            + np.roll(ocean_f, -1, axis=0) + np.roll(ocean_f, 1, axis=0)
        )
        return np.divide(total, count, out=np.zeros_like(total), where=count > 0)

    for _ in range(COASTAL_SMOOTHING_ITERATIONS):
        avg_u, avg_v = neighbor_average(u), neighbor_average(v)
        u = u + COASTAL_SMOOTHING_WEIGHT * (avg_u - u)
        v = v + COASTAL_SMOOTHING_WEIGHT * (avg_v - v)
        u = u + COASTAL_RESTORE_FRACTION * (ambient_u - u)
        v = v + COASTAL_RESTORE_FRACTION * (ambient_v - v)
        u, v = _deflect_into_coast(u, v, is_ocean, lat_deg)
    return u, v


def _land_swirl_current(wind_u: np.ndarray, wind_v: np.ndarray, is_ocean: np.ndarray, lat_deg: np.ndarray, world_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotational swirl around every landmass -- keyed off nearest *land cell* (a cKDTree
    chord-distance query, equivalent to plate-sim's own BFS-based `_land_swirl_current`,
    which is also nearest-cell, not a connected-landmass grouping), ramping up from 0 at the
    coast to full strength at SWIRL_PEAK_DIST_DEG, decaying exponentially beyond that. Real
    ocean gyres are wind-driven *and* basin-shaped; Ekman + coastal deflection alone gives
    flat latitude-banded flow with no closed loops (see module docstring) -- this is what
    actually produces current-like circulation."""
    height, width = is_ocean.shape
    is_land = ~is_ocean
    if not np.any(is_land):
        return np.zeros((height, width)), np.zeros((height, width))

    land_xyz = world_xyz.reshape(-1, 3)[is_land.reshape(-1)]
    tree = cKDTree(land_xyz)
    chord_dist, nearest_idx = tree.query(world_xyz.reshape(-1, 3))
    dist_rad = 2.0 * np.arcsin(np.clip(chord_dist / 2.0, 0.0, 1.0)).reshape(height, width)
    dist_deg = np.degrees(dist_rad)

    # Outward-from-land direction: this cell's own position minus the nearest land cell's --
    # a genuine 3D tangent-plane-ish direction rather than plate-sim's lat/lon offset (their
    # simplification for an equirectangular grid; a real vector difference works directly here).
    nearest_land_xyz = land_xyz[nearest_idx].reshape(height, width, 3)
    to_here = world_xyz - nearest_land_xyz
    # Project out the radial (through-the-sphere) component, keep only the tangent-plane
    # part, matching this module's other 2D (east, north) vector fields. `north = p x east`
    # is the standard sphere tangent-frame identity (verified: at (lat, lon) it equals
    # d/dlat of the position exactly).
    east = np.stack([-np.sin(np.radians(_lon_grid(width, height))), np.cos(np.radians(_lon_grid(width, height))), np.zeros((height, width))], axis=-1)
    north = np.cross(world_xyz, east)
    radial_east = np.sum(to_here * east, axis=-1)
    radial_north = np.sum(to_here * north, axis=-1)
    radial_mag = np.hypot(radial_east, radial_north) + _EPS
    r_east, r_north = radial_east / radial_mag, radial_north / radial_mag

    hemisphere = np.where(lat_deg >= 0, 1.0, -1.0)[:, None]
    t_east, t_north = _rotate_90(r_east, r_north, hemisphere)

    ramp_up = np.clip(dist_deg / SWIRL_PEAK_DIST_DEG, 0.0, 1.0)
    ramp_down = np.exp(-np.maximum(dist_deg - SWIRL_PEAK_DIST_DEG, 0.0) / SWIRL_DECAY_DIST_DEG)
    profile = np.where(dist_deg <= SWIRL_PEAK_DIST_DEG, ramp_up, ramp_down)
    profile = np.where(dist_deg <= SWIRL_MAX_SEARCH_DIST_DEG, profile, 0.0)
    wind_speed = np.hypot(wind_u, wind_v)
    magnitude = SWIRL_SPEED_FRACTION * wind_speed * profile
    return t_east * magnitude, t_north * magnitude


def _lon_grid(width: int, height: int) -> np.ndarray:
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    return np.repeat(lon_deg[None, :], height, axis=0)


def _circumglobal_row_boost(is_ocean: np.ndarray) -> np.ndarray:
    is_open_row = ~np.any(~is_ocean, axis=1)
    return np.where(is_open_row[:, None], CIRCUMGLOBAL_SPEEDUP_FACTOR, 1.0)


def _upstream_obstacle_step_distance(is_obstacle: np.ndarray, dir_u: np.ndarray, dir_v: np.ndarray, lat_deg: np.ndarray, steps: int) -> np.ndarray:
    height, width = is_obstacle.shape
    nearest_step = np.full((height, width), steps + 1, dtype=np.int32)
    found = np.zeros((height, width), dtype=bool)
    for step in range(1, steps + 1):
        upstream_here = _sample_at_offset(is_obstacle, -dir_u, -dir_v, step * _REFERENCE_CELL_DEG, lat_deg)
        newly_found = upstream_here & ~found
        nearest_step[newly_found] = step
        found |= upstream_here
    return nearest_step


def _apply_current_wake(u: np.ndarray, v: np.ndarray, is_ocean: np.ndarray, lat_deg: np.ndarray, mixing_noise: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    speed = np.hypot(u, v)
    has_current = speed > 1e-6
    dir_u = np.divide(u, speed, out=np.zeros_like(speed), where=has_current)
    dir_v = np.divide(v, speed, out=np.zeros_like(speed), where=has_current)

    nearest_step = _upstream_obstacle_step_distance(~is_ocean, dir_u, dir_v, lat_deg, WAKE_LOOKBACK_STEPS)
    in_wake = has_current & (nearest_step <= WAKE_LOOKBACK_STEPS)
    wake_strength = np.where(in_wake, 1.0 - (nearest_step - 1) / WAKE_LOOKBACK_STEPS, 0.0)
    speed_factor = 1.0 - wake_strength * (1.0 - WAKE_MIN_SPEED_FACTOR)

    mix_angle = np.radians((mixing_noise - 0.5) * 2.0 * WAKE_MIXING_ANGLE_MAX_DEG) * wake_strength
    cos_m, sin_m = np.cos(mix_angle), np.sin(mix_angle)
    u_mixed = u * cos_m - v * sin_m
    v_mixed = u * sin_m + v * cos_m
    return u_mixed * speed_factor, v_mixed * speed_factor


def compute_ocean_currents(
    wind_u: np.ndarray, wind_v: np.ndarray, is_ocean: np.ndarray, lat_deg: np.ndarray, world_xyz: np.ndarray, mixing_noise: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Ekman base + land swirl + circumglobal boost, then coastal deflection/smoothing/wake.
    Zero over land throughout."""
    u, v = _ekman_current(wind_u, wind_v, lat_deg)
    swirl_u, swirl_v = _land_swirl_current(wind_u, wind_v, is_ocean, lat_deg, world_xyz)
    u, v = u + swirl_u, v + swirl_v
    boost = _circumglobal_row_boost(is_ocean)
    u, v = u * boost, v * boost

    ambient_u, ambient_v = u.copy(), v.copy()
    u, v = _deflect_into_coast(u, v, is_ocean, lat_deg)
    u, v = _smooth_along_coast(u, v, is_ocean, lat_deg, ambient_u, ambient_v)
    u, v = _apply_current_wake(u, v, is_ocean, lat_deg, mixing_noise)

    u = np.where(is_ocean, u, 0.0)
    v = np.where(is_ocean, v, 0.0)
    return u, v


def compute_ocean_swells(current_u: np.ndarray, current_v: np.ndarray, is_ocean: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Convergence (negative divergence) of the current field, weighted-sampled at up to
    MAX_OCEAN_SWELLS points -- visible where two currents collide and pile water up."""
    east = np.roll(current_u, -1, axis=1)
    west = np.roll(current_u, 1, axis=1)
    north = np.roll(current_v, 1, axis=0)
    south = np.roll(current_v, -1, axis=0)
    divergence = (east - west) / 2.0 + (north - south) / 2.0
    convergence = np.clip(-divergence, 0.0, None)
    convergence = np.where(is_ocean, convergence, 0.0)
    convergence = np.clip(convergence / OCEAN_SWELL_CONVERGENCE_REFERENCE, 0.0, 1.0)

    flat = convergence.reshape(-1)
    picked = _weighted_sample_without_replacement(rng, flat, MAX_OCEAN_SWELLS)
    rows, cols = np.unravel_index(picked, convergence.shape)
    return rows, cols


# ---------------------------------------------------------------------------------------
# Ocean temperature advection
# ---------------------------------------------------------------------------------------

ADVECTION_DISTANCE_DEG = 6.0


def advect_ocean_temperature(baseline_c: np.ndarray, current_u: np.ndarray, current_v: np.ndarray, is_ocean: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
    """Semi-Lagrangian: each ocean cell samples the baseline zonal temperature at a fixed
    distance *upstream* along its own current -- the surface expression of real
    thermohaline heat transport (warm equatorial water carried poleward along boundary
    currents)."""
    height, width = baseline_c.shape
    speed = np.hypot(current_u, current_v)
    has_current = speed > 1e-6
    dir_u = np.divide(current_u, speed, out=np.zeros_like(speed), where=has_current)
    dir_v = np.divide(current_v, speed, out=np.zeros_like(speed), where=has_current)

    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    cos_lat = np.clip(np.cos(np.radians(lat_grid)), 0.15, 1.0)

    src_lat = lat_grid - dir_v * ADVECTION_DISTANCE_DEG
    src_lon = lon_deg[None, :] - dir_u * ADVECTION_DISTANCE_DEG / cos_lat
    src_row = np.clip(np.round((90.0 - src_lat) / (180.0 / height) - 0.5).astype(np.int64), 0, height - 1)
    src_col = np.round((src_lon + 180.0) / (360.0 / width) - 0.5).astype(np.int64) % width

    sampled = baseline_c[src_row, src_col]
    valid = is_ocean & has_current & is_ocean[src_row, src_col]
    return np.where(valid, sampled, baseline_c)


# ---------------------------------------------------------------------------------------
# Humidity (evaporation ceiling + wind-driven 2D advection) and precipitation
# ---------------------------------------------------------------------------------------

EVAPORATION_REFERENCE_TEMP_C = 20.0
MIN_EVAPORATION_CEILING = 0.3
MAX_EVAPORATION_CEILING = 1.4
MIN_WIND_RETENTION_FACTOR = 0.5
RETENTION_PER_CELL = 0.96
OROGRAPHIC_LIFT_SCALE_M = 600.0
OROGRAPHIC_RAIN_SHADOW_FACTOR = 0.6
PRECIP_HUMIDITY_COEFFICIENT_MM = 1500.0
OROGRAPHIC_PRECIPITATION_COEFFICIENT = 1.0


def _evaporation_ceiling(ocean_temperature_c: np.ndarray) -> np.ndarray:
    return np.clip(ocean_temperature_c / EVAPORATION_REFERENCE_TEMP_C, MIN_EVAPORATION_CEILING, MAX_EVAPORATION_CEILING)


def _retention_factor(elevation_factor_cell: np.ndarray) -> np.ndarray:
    wind_factor = MIN_WIND_RETENTION_FACTOR + (1.0 - MIN_WIND_RETENTION_FACTOR) * elevation_factor_cell
    return RETENTION_PER_CELL * wind_factor


def _orographic_retained_fraction(gain_m: np.ndarray) -> np.ndarray:
    saturation = 1.0 - np.exp(-np.clip(gain_m, 0.0, None) / OROGRAPHIC_LIFT_SCALE_M)
    retained = 1.0 - saturation * (1.0 - OROGRAPHIC_RAIN_SHADOW_FACTOR)
    return np.where(gain_m > 0, retained, 1.0)


def _humidity_zonal_sweep(is_ocean: np.ndarray, elevation_m: np.ndarray, evap_ceiling: np.ndarray, elevation_factor: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    height, width = is_ocean.shape
    zonal_dir = zonal_direction_for_lat(lat_deg)
    rows = np.arange(height)

    humidity = np.zeros((height, width))
    orographic = np.zeros((height, width))
    moisture = np.zeros(height)
    prev_elev = np.zeros(height)
    for lap in range(2):
        for i in range(width):
            cols = np.where(zonal_dir == 1, i, width - 1 - i).astype(int)
            elev_i = elevation_m[rows, cols]
            ocean_i = is_ocean[rows, cols]
            ceiling_i = evap_ceiling[rows, cols]
            ef_i = elevation_factor[rows, cols]

            land_moisture = moisture * _retention_factor(ef_i)
            after_source = np.where(ocean_i, ceiling_i, land_moisture)
            retained = _orographic_retained_fraction(elev_i - prev_elev)
            new_moisture = np.where(ocean_i, after_source, after_source * retained)
            dump = np.where(ocean_i, 0.0, after_source - new_moisture)

            if lap == 1:
                humidity[rows, cols] = new_moisture
                orographic[rows, cols] = dump
            moisture = new_moisture
            prev_elev = elev_i
    return humidity, orographic


def _meridional_bands(lat_deg: np.ndarray) -> list[np.ndarray]:
    """Contiguous latitude bands (trade/westerlies/polar per hemisphere), each an array of
    row indices ordered from the flow's diverging edge to its converging edge."""
    row_direction = -np.sign(meridional_direction_for_lat(lat_deg))  # +1 = increasing row index
    bands = []
    start = 0
    n = len(lat_deg)
    for i in range(1, n + 1):
        if i == n or row_direction[i] != row_direction[start]:
            rows = np.arange(start, i)
            if row_direction[start] < 0:
                rows = rows[::-1]
            bands.append(rows)
            start = i
    return bands


def _humidity_meridional_sweep(is_ocean: np.ndarray, elevation_m: np.ndarray, evap_ceiling: np.ndarray, elevation_factor: np.ndarray, lat_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    height, width = is_ocean.shape
    humidity = np.zeros((height, width))
    orographic = np.zeros((height, width))
    for band_rows in _meridional_bands(lat_deg):
        moisture = np.zeros(width)
        prev_elev = np.zeros(width)
        for r in band_rows:
            elev_r = elevation_m[r]
            ocean_r = is_ocean[r]
            ceiling_r = evap_ceiling[r]
            ef_r = elevation_factor[r]

            land_moisture = moisture * _retention_factor(ef_r)
            after_source = np.where(ocean_r, ceiling_r, land_moisture)
            retained = _orographic_retained_fraction(elev_r - prev_elev)
            new_moisture = np.where(ocean_r, after_source, after_source * retained)
            dump = np.where(ocean_r, 0.0, after_source - new_moisture)

            humidity[r] = new_moisture
            orographic[r] = dump
            moisture = new_moisture
            prev_elev = elev_r
    return humidity, orographic


def compute_humidity(
    is_ocean: np.ndarray, elevation_m: np.ndarray, ocean_temperature_c: np.ndarray,
    wind_u: np.ndarray, wind_v: np.ndarray, elevation_factor: np.ndarray, lat_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaporation ceiling over ocean (from local ocean temperature), wind-driven 2D
    advection onto land -- a zonal sweep and a meridional sweep, blended per-cell by each
    wind component's share of total wind magnitude (a cell whose wind is mostly north-south
    gets its humidity mostly from the meridional pass, and vice versa). No evapotranspiration
    term (needs vegetation, which doesn't exist here). Returns (humidity, orographic_dump)."""
    evap_ceiling = _evaporation_ceiling(ocean_temperature_c)
    humidity_zonal, oro_zonal = _humidity_zonal_sweep(is_ocean, elevation_m, evap_ceiling, elevation_factor, lat_deg)
    humidity_meridional, oro_meridional = _humidity_meridional_sweep(is_ocean, elevation_m, evap_ceiling, elevation_factor, lat_deg)

    abs_u, abs_v = np.abs(wind_u), np.abs(wind_v)
    total = abs_u + abs_v
    w_u = np.divide(abs_u, total, out=np.full_like(total, 0.5), where=total > 1e-9)
    w_v = 1.0 - w_u

    humidity = w_u * humidity_zonal + w_v * humidity_meridional
    orographic = w_u * oro_zonal + w_v * oro_meridional
    return humidity, orographic


def compute_precipitation(humidity: np.ndarray, orographic_dump: np.ndarray) -> np.ndarray:
    """Purely a function of humidity and orographic lift (wind carrying damp air over
    mountains) -- no zonal latitude-climatology baseline, cut per the user's own
    description."""
    return PRECIP_HUMIDITY_COEFFICIENT_MM * humidity + OROGRAPHIC_PRECIPITATION_COEFFICIENT * PRECIP_HUMIDITY_COEFFICIENT_MM * orographic_dump


# ---------------------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------------------


def compute_climate(world: World, height: int = GRID_HEIGHT, width: int = GRID_WIDTH) -> ClimateFields:
    """Runs the full climate pipeline against the world's *current* plate state. See module
    docstring for the pipeline order and why it's structured this way."""
    lat_deg, lon_deg, world_xyz = _build_grid(height, width)
    elevation_m, is_ocean = _sample_elevation_and_crust(world, world_xyz)

    insolation_row = compute_insolation(lat_deg, world.axial_tilt_deg)
    land_temperature_c = compute_land_temperature(insolation_row, elevation_m)
    ocean_baseline_c = compute_ocean_temperature_baseline(insolation_row, height, width)
    surface_temperature_c = np.where(is_ocean, ocean_baseline_c, land_temperature_c)

    wind_u, wind_v, elevation_factor = compute_wind(lat_deg, elevation_m, surface_temperature_c)

    # A cached-per-world-state noise texture standing in for turbulent current mixing --
    # deterministic in (seed, elapsed_years) so it doesn't flicker between renders of the
    # same world state, matching every other RNG use in this codebase.
    mix_rng = np.random.default_rng((world.seed, round(world.elapsed_years)))
    mixing_noise = mix_rng.random((height, width))

    current_u, current_v = compute_ocean_currents(wind_u, wind_v, is_ocean, lat_deg, world_xyz, mixing_noise)
    swell_rows, swell_cols = compute_ocean_swells(current_u, current_v, is_ocean, mix_rng)

    ocean_temperature_c = advect_ocean_temperature(ocean_baseline_c, current_u, current_v, is_ocean, lat_deg)
    air_temperature_c = compute_air_temperature(land_temperature_c, ocean_temperature_c, is_ocean, world_xyz)

    humidity, orographic_dump = compute_humidity(is_ocean, elevation_m, ocean_temperature_c, wind_u, wind_v, elevation_factor, lat_deg)
    precipitation_mm = compute_precipitation(humidity, orographic_dump)

    return ClimateFields(
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        world_xyz=world_xyz,
        elevation_m=elevation_m,
        is_ocean=is_ocean,
        land_temperature_c=land_temperature_c,
        ocean_temperature_c=ocean_temperature_c,
        air_temperature_c=air_temperature_c,
        wind_u=wind_u,
        wind_v=wind_v,
        current_u=current_u,
        current_v=current_v,
        humidity=humidity,
        precipitation_mm=precipitation_mm,
        swell_rows=swell_rows,
        swell_cols=swell_cols,
    )

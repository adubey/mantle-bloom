"""RFC 0019 offline coarse moisture budget.

Run the analytic checks with bin/debug/check_rfc0019_coarse_moisture.py.
Inputs are latitude rows ordered south to north, with winds positive east and north.
Only aligned, nested spherical grids are supported by this prototype.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from rfc0019_transport_prototypes import EARTH_RADIUS_M, sphere_grid

CLOSURE_LIMIT = 1e-8
SECONDS_PER_DAY = 86_400.0


@dataclass(frozen=True)
class Reservoirs:
    native_ids: np.ndarray
    native_area: np.ndarray
    area: np.ndarray
    is_ocean: np.ndarray
    coarse_cell: np.ndarray
    tagged_source: dict[str, np.ndarray]
    coarse_shape: tuple[int, int]
    factor: int


@dataclass(frozen=True)
class TransportFaces:
    first: np.ndarray
    second: np.ndarray
    forward: np.ndarray  # advective capacity, m²/s from first to second
    reverse: np.ndarray
    exchange: np.ndarray  # symmetric mixing conductance, m²/s


@dataclass(frozen=True)
class FaceInterception:
    """Fraction of advected vapor that rains at the receiving face endpoint."""

    forward: np.ndarray
    reverse: np.ndarray


@dataclass(frozen=True)
class Budget:
    mass: np.ndarray
    rain: np.ndarray
    face_flux: np.ndarray  # kg/s before interception, positive from first to second
    global_relative: float
    max_cell_relative: float
    summed_cell_relative: float
    clipped_negative_mass_kg: float = 0.0
    intercepted_rain_kg_s: np.ndarray | None = None


def aggregate_inputs(
    is_ocean: np.ndarray, sources: dict[str, np.ndarray], factor: int = 1,
) -> Reservoirs:
    """Group native cells by coarse footprint and surface class, summing kg/s by tag."""
    nlat, nlon = is_ocean.shape
    if factor < 1 or nlat % factor or nlon % factor or nlon // factor < 3:
        raise ValueError("expected a nested grid with at least three coarse longitude cells")
    coarse_shape = (nlat // factor, nlon // factor)
    _, _, native_area = sphere_grid(nlat, nlon)
    row = np.arange(nlat)[:, None] // factor
    col = np.arange(nlon)[None, :] // factor
    coarse_id = row * coarse_shape[1] + col
    class_id = 2 * coarse_id + is_ocean.astype(np.intp)
    classes, inverse = np.unique(class_id.ravel(), return_inverse=True)
    native_ids = inverse.reshape(is_ocean.shape)
    area = np.bincount(inverse, weights=native_area.ravel())
    tagged_source = {}
    for tag, values in sources.items():
        if values.shape != is_ocean.shape or not np.all(np.isfinite(values) & (values >= 0)):
            raise ValueError(f"invalid source grid for {tag}")
        tagged_source[tag] = np.bincount(inverse, weights=values.ravel())
    return Reservoirs(
        native_ids, native_area, area, classes % 2 == 1, classes // 2,
        tagged_source, coarse_shape, factor,
    )


def _native_faces(first, second, flow, length, spacing, external, diffusivity):
    """Keep coarse boundary segments and internal shorelines as distinct faces."""
    present = first != second
    return TransportFaces(
        first[present],
        second[present],
        np.maximum(flow, 0.0)[present],
        np.maximum(-flow, 0.0)[present],
        np.broadcast_to(np.where(external, diffusivity * length / spacing, 0.0), flow.shape)[present],
    )


def _coastal_exchange(grid: Reservoirs, coast_days: float) -> TransportFaces:
    land_ids = np.flatnonzero(~grid.is_ocean)
    ocean_ids = np.flatnonzero(grid.is_ocean)
    _, land_index, ocean_index = np.intersect1d(
        grid.coarse_cell[land_ids], grid.coarse_cell[ocean_ids], assume_unique=True, return_indices=True,
    )
    land = land_ids[land_index]
    ocean = ocean_ids[ocean_index]
    total_area = grid.area[land] + grid.area[ocean]
    exchange = grid.area[land] * (grid.area[ocean] / total_area) / (coast_days * SECONDS_PER_DAY)
    zeros = np.zeros(land.size)
    return TransportFaces(land, ocean, zeros, zeros, exchange)


def build_transport_rates(
    grid: Reservoirs, u: np.ndarray, v: np.ndarray, diffusivity: float, coast_days: float = 2.0,
    *, north_face_delta: np.ndarray | None = None,
) -> TransportFaces:
    """Build shared face capacities. Convert to donor mass rates in the solver."""
    ids = grid.native_ids
    nlat, nlon = ids.shape
    if u.shape != ids.shape or v.shape != ids.shape or not np.all(np.isfinite([u, v])):
        raise ValueError("invalid wind grid")
    if not np.isfinite(diffusivity) or diffusivity < 0 or not np.isfinite(coast_days) or coast_days <= 0:
        raise ValueError("expected nonnegative diffusivity and positive coastal exchange time")

    _, edges, _ = sphere_grid(nlat, nlon)
    _, coarse_edges, _ = sphere_grid(*grid.coarse_shape)
    coarse_lat = 0.5 * (coarse_edges[:-1] + coarse_edges[1:])
    lon_width = 2 * np.pi / nlon
    east_length = (EARTH_RADIUS_M * np.diff(edges))[:, None]
    east_spacing = (EARTH_RADIUS_M * np.cos(coarse_lat).repeat(grid.factor) * lon_width * grid.factor)[:, None]
    east = _native_faces(
        ids, np.roll(ids, -1, axis=1),
        east_length * 0.5 * (u + np.roll(u, -1, axis=1)),
        east_length, east_spacing, (np.arange(nlon) + 1) % grid.factor == 0, diffusivity,
    )

    north_length = (EARTH_RADIUS_M * np.cos(edges[1:-1]) * lon_width)[:, None]
    north_boundary = np.arange(1, nlat) % grid.factor == 0
    # Internal face spacing is unused. The coarse center spacing applies only
    # where native segments lie along a boundary between coarse cells.
    north_spacing = np.ones(nlat - 1)
    north_spacing[north_boundary] = EARTH_RADIUS_M * np.diff(coarse_lat)
    north_velocity = 0.5 * (v[:-1] + v[1:])
    if north_face_delta is not None:
        if north_face_delta.shape != (nlat - 1,) or not np.all(np.isfinite(north_face_delta)):
            raise ValueError("expected a finite velocity adjustment per internal latitude face")
        north_velocity += north_face_delta[:, None]
    north = _native_faces(
        ids[:-1], ids[1:], north_length * north_velocity,
        north_length, north_spacing[:, None], north_boundary[:, None], diffusivity,
    )
    coast = _coastal_exchange(grid, coast_days)
    return TransportFaces(
        np.concatenate((east.first, north.first, coast.first)),
        np.concatenate((east.second, north.second, coast.second)),
        np.concatenate((east.forward, north.forward, coast.forward)),
        np.concatenate((east.reverse, north.reverse, coast.reverse)),
        np.concatenate((east.exchange, north.exchange, coast.exchange)),
    )


def net_export_operator(
    area: np.ndarray, faces: TransportFaces, interception: FaceInterception | None = None,
) -> sparse.csr_matrix:
    """Map reservoir mass in kg to net export in kg/s."""
    origin = np.concatenate((faces.first, faces.second))
    destination = np.concatenate((faces.second, faces.first))
    rate = np.concatenate((faces.forward + faces.exchange, faces.reverse + faces.exchange)) / area[origin]
    count = area.size
    outgoing = np.bincount(origin, weights=rate, minlength=count)
    delivered = rate
    if interception is not None:
        _validate_interception(faces, interception)
        delivered = np.concatenate((
            faces.forward * (1.0 - interception.forward) + faces.exchange,
            faces.reverse * (1.0 - interception.reverse) + faces.exchange,
        )) / area[origin]
    incoming = sparse.csr_matrix((-delivered, (destination, origin)), shape=(count, count))
    return sparse.diags(outgoing, format="csr") + incoming


def _validate_interception(faces, interception):
    for fraction in (interception.forward, interception.reverse):
        if fraction.shape != faces.forward.shape or not np.all((fraction >= 0.0) & (fraction <= 1.0)):
            raise ValueError("expected a fraction in [0, 1] for each face and direction")


def _audit_budget(area, source, sink, faces, mass, clipped=0.0, *, overflow=None, interception=None):
    first, second = faces.first, faces.second
    count = area.size
    water = mass / area
    rain = sink * mass
    if overflow is not None:
        rain += overflow
    intercepted = None
    if interception is not None:
        intercepted = (
            np.bincount(second, weights=faces.forward * water[first] * interception.forward, minlength=count)
            + np.bincount(first, weights=faces.reverse * water[second] * interception.reverse, minlength=count)
        )
        rain += intercepted
    flux = (
        faces.forward * water[first] - faces.reverse * water[second]
        + faces.exchange * (water[first] - water[second])
    )
    incoming_net = (
        np.bincount(second, weights=flux, minlength=count)
        - np.bincount(first, weights=flux, minlength=count)
    )
    abs_local = np.abs(source + incoming_net - rain)
    total_source = float(source.sum())
    scale = total_source if total_source > 0 else 1.0
    closure = abs(float(rain.sum()) - total_source) / scale
    max_local = float(abs_local.max()) / scale
    summed_local = float(abs_local.sum()) / scale
    if not np.isfinite(summed_local) or max(closure, summed_local) > CLOSURE_LIMIT:
        raise AssertionError(f"rain and transport failed closure: global={closure:g}, local sum={summed_local:g}")
    return Budget(mass, rain, flux, closure, max_local, summed_local, clipped, intercepted)


def _stationary_operator(area, source, sink, faces, interception=None) -> sparse.csr_matrix:
    if area.ndim != 1 or not area.size or not np.all(np.isfinite(area) & (area > 0)):
        raise ValueError("expected positive finite reservoir areas")
    if source.shape != area.shape or not np.all(np.isfinite(source) & (source >= 0)):
        raise ValueError("source must be finite and nonnegative, in kg/s per reservoir")
    if sink.shape != area.shape or not np.all(np.isfinite(sink) & (sink > 0)):
        raise ValueError("a positive rain rate is required for each reservoir")
    return sparse.diags(sink, format="csr") + net_export_operator(area, faces, interception)


def solve_budget(area: np.ndarray, source: np.ndarray, sink: np.ndarray, faces: TransportFaces) -> Budget:
    """Solve for absolute column mass and audit shared face fluxes independently."""
    operator = _stationary_operator(area, source, sink, faces)
    mass = spsolve(operator.tocsc(), source)
    if not np.all(np.isfinite(mass)):
        raise AssertionError("stationary solve returned nonfinite mass")
    tolerance = 100 * np.finfo(mass.dtype).eps * np.max(np.abs(mass))
    if mass.min() < -tolerance:
        raise AssertionError("stationary solve returned negative mass beyond roundoff")
    audited = _audit_budget(area, source, sink, faces, mass)
    clipped = float(-mass[mass < 0].sum())
    if clipped:
        return _audit_budget(area, source, sink, faces, np.maximum(mass, 0.0), clipped)
    return audited


def solve_bounded_budget(area, source, sink, faces, mass_cap, *, interception=None) -> tuple[Budget, np.ndarray, int]:
    """Solve the stationary capped column and return budget, overflow rain (kg/s), and solve count."""
    operator = _stationary_operator(area, source, sink, faces, interception)
    if mass_cap.shape != area.shape or not np.all(np.isfinite(mass_cap) & (mass_cap >= 0)):
        raise ValueError("expected a nonnegative finite mass cap in kg per reservoir")
    mass = mass_cap.copy()
    residual = source - operator @ mass
    # The operator has a positive diagonal and nonpositive off-diagonals. Releasing
    # the cells that cannot sustain their cap only decreases mass and incoming flux.
    # A released cell cannot return to its cap as the remaining caps are released.
    solves = 0
    free = np.zeros(area.size, dtype=bool)
    while np.any(release := (~free & (residual < 0.0))):
        free |= release
        free_operator = operator[free][:, free].tocsc()
        # Include imports from still-saturated neighbors in the free-cell RHS.
        mass[free] = 0.0
        rhs = (source - operator @ mass)[free]
        mass[free] = spsolve(free_operator, rhs)
        residual = source - operator @ mass
        solves += 1
    tolerance = 100 * np.finfo(float).eps * mass_cap
    if not np.all(np.isfinite(mass) & (mass >= 0) & (mass <= mass_cap + tolerance)):
        raise AssertionError("bounded stationary solve violated its mass limits")
    overflow = np.where(free, 0.0, residual)
    budget = _audit_budget(area, source, sink, faces, mass, overflow=overflow, interception=interception)
    return budget, overflow, solves


def distribute_rain(grid: Reservoirs, rain: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Return native-cell rain in kg/s, preserving each reservoir's sink total."""
    if rain.shape != grid.area.shape or not np.all(np.isfinite(rain) & (rain >= 0)):
        raise ValueError("invalid reservoir rain")
    if weights.shape != grid.native_ids.shape or not np.all(np.isfinite(weights) & (weights > 0)):
        raise ValueError("rain placement weights must be finite and positive")
    weighted_area = grid.native_area * weights
    totals = np.bincount(grid.native_ids.ravel(), weights=weighted_area.ravel(), minlength=rain.size)
    ids = grid.native_ids
    return rain[ids] * (weighted_area / totals[ids])

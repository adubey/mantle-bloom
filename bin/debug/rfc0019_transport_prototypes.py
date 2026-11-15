#!/usr/bin/env python3
"""RFC 0019 G2: compare conservative finite-volume and wind-arrow prototypes.

Usage: backend/.venv/bin/python bin/debug/rfc0019_transport_prototypes.py
Runs small analytic spherical-grid cases, checks positivity and global closure,
and reports steady-solve residuals and wall time. No World or production code is used.
"""
from __future__ import annotations

import json
import argparse
import time
import warnings

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import MatrixRankWarning, spsolve

EARTH_RADIUS_M = 6_371_000.0


def sphere_grid(nlat: int, nlon: int):
    """Return equal-area latitude centers, latitude edges, and cell areas."""
    latitude_edges = np.linspace(-np.pi / 2, np.pi / 2, nlat + 1)
    sine_midpoints = (np.sin(latitude_edges[:-1]) + np.sin(latitude_edges[1:])) / 2
    latitude = np.arcsin(sine_midpoints)
    longitude_width = 2 * np.pi / nlon
    row_area = (
        EARTH_RADIUS_M**2
        * longitude_width
        * (np.sin(latitude_edges[1:]) - np.sin(latitude_edges[:-1]))
    )
    area = np.broadcast_to(row_area[:, None], (nlat, nlon)).copy()
    return latitude, latitude_edges, area


def transition_rates(u, v, latitude, latitude_edges, area, finite_volume: bool):
    """Build directed neighbor rates q[i,j] in 1/s for i -> j transfers."""
    nlat, nlon = area.shape
    longitude_width = 2 * np.pi / nlon
    rates = []
    for row in range(nlat):
        east_distance = EARTH_RADIUS_M * np.cos(latitude[row]) * longitude_width
        for col in range(nlon):
            cell = row * nlon + col
            east = row * nlon + (col + 1) % nlon
            if finite_volume:
                face_u = 0.5 * (u[row, col] + u[row, (col + 1) % nlon])
                zonal_edge_length = EARTH_RADIUS_M * (
                    latitude_edges[row + 1] - latitude_edges[row]
                )
                donor = cell if face_u > 0 else east
                zonal_rate = abs(face_u) * zonal_edge_length / area.ravel()[donor]
                if face_u > 0:
                    rates.append((cell, east, zonal_rate))
                elif face_u < 0:
                    rates.append((east, cell, zonal_rate))
            else:
                if u[row, col] > 0:
                    rates.append((cell, east, u[row, col] / east_distance))
                elif u[row, col] < 0:
                    west = row * nlon + (col - 1) % nlon
                    rates.append((cell, west, -u[row, col] / east_distance))
                for sign in (-1, 1):
                    neighbor_row = row + sign
                    if not 0 <= neighbor_row < nlat or v[row, col] * sign <= 0:
                        continue
                    neighbor = neighbor_row * nlon + col
                    center_distance = EARTH_RADIUS_M * abs(
                        latitude[neighbor_row] - latitude[row]
                    )
                    rates.append((cell, neighbor, abs(v[row, col]) / center_distance))

            if finite_volume and row < nlat - 1:
                neighbor_row = row + 1
                neighbor = neighbor_row * nlon + col
                face_v = 0.5 * (v[row, col] + v[neighbor_row, col])
                edge_latitude = latitude_edges[neighbor_row]
                meridional_edge_length = (
                    EARTH_RADIUS_M * np.cos(edge_latitude) * longitude_width
                )
                donor = cell if face_v > 0 else neighbor
                rate = abs(face_v) * meridional_edge_length / area.ravel()[donor]
                if face_v > 0:
                    rates.append((cell, neighbor, rate))
                elif face_v < 0:
                    rates.append((neighbor, cell, rate))
    return rates


def detect_sinkless_closed_class(n_cells, rates, sinks):
    """Reject a recurrent transition class that has no positive rain sink."""
    rows = [source for source, _, _ in rates]
    cols = [destination for _, destination, _ in rates]
    graph = sparse.csr_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(n_cells, n_cells)
    )
    count, labels = connected_components(graph, directed=True, connection="strong")
    outgoing = np.zeros(count, dtype=bool)
    for source, destination, _ in rates:
        if labels[source] != labels[destination]:
            outgoing[labels[source]] = True
    for component in range(count):
        member_sinks = sinks[labels == component]
        if not outgoing[component] and np.all(member_sinks <= 0):
            raise ValueError("closed transition class has no rain sink")


def solve_stationary(source, sinks, rates):
    """Solve the absolute stationary mass balance and return closure metrics."""
    n_cells = source.size
    detect_sinkless_closed_class(n_cells, rates, sinks)
    diagonal = sinks.copy()
    rows = []
    cols = []
    values = []
    for origin, destination, rate in rates:
        diagonal[origin] += rate
        rows.append(destination)
        cols.append(origin)
        values.append(-rate)
    rows.extend(range(n_cells))
    cols.extend(range(n_cells))
    values.extend(diagonal)
    matrix = sparse.csr_matrix((values, (rows, cols)), shape=(n_cells, n_cells))
    with warnings.catch_warnings():
        warnings.simplefilter("error", MatrixRankWarning)
        mass = np.asarray(spsolve(matrix, source)).reshape(source.shape)
    inflow = np.zeros(n_cells)
    outflow = np.zeros(n_cells)
    for origin, destination, rate in rates:
        transfer = rate * mass[origin]
        outflow[origin] += transfer
        inflow[destination] += transfer
    rain = sinks * mass
    residual = source + inflow - outflow - rain
    total_source = float(source.sum())
    closure = abs(float(rain.sum()) - total_source) / max(total_source, 1e-300)
    local = float(np.max(np.abs(residual))) / max(total_source, 1e-300)
    if not np.all(np.isfinite(mass)) or np.min(mass) < -1e-10:
        raise AssertionError("stationary solve produced invalid or negative water mass")
    return mass, rain, {"global_closure_relative": closure, "max_cell_residual_relative": local}


def cap_rain_to_available(opening_mass, potential_rain, delta_t):
    """Limit a finite-interval sink to opening water plus current source input."""
    available = opening_mass
    removed = min(potential_rain * delta_t, available)
    return available - removed, removed


def run_case(nlat=8, nlon=16, source_mode="ocean", wind_mode="circulation", rain_days=9.1):
    latitude, latitude_edges, area = sphere_grid(nlat, nlon)
    u = np.zeros(area.shape)
    v = np.zeros(area.shape)
    if wind_mode == "circulation":
        u[:] = 5.0
        v[:] = 0.4 * np.sin(latitude)[:, None]
    elif wind_mode == "convergent":
        u[:] = 2.0
        v[:] = -np.sign(latitude)[:, None] * 0.3
    elif wind_mode == "divergent":
        u[:] = 2.0
        v[:] = np.sign(latitude)[:, None] * 0.3
    elif wind_mode == "calm":
        pass
    elif wind_mode == "loop":
        u[:] = 5.0
    else:
        raise ValueError(f"unknown wind mode: {wind_mode}")

    land = np.broadcast_to(np.arange(nlon)[None, :] >= nlon // 2, area.shape)
    source_per_area = np.zeros(area.shape)
    if source_mode == "ocean":
        source_per_area[~land] = 3e-5
    elif source_mode == "uniform":
        source_per_area[:] = 3e-5
    elif source_mode != "zero":
        raise ValueError(f"unknown source mode: {source_mode}")
    source = (source_per_area * area).ravel()
    sinks = np.full(source.size, 1.0 / (rain_days * 86_400.0))

    outputs = {}
    for name, finite_volume in (("finite_volume", True), ("wind_arrows", False)):
        rates = transition_rates(u, v, latitude, latitude_edges, area, finite_volume)
        started = time.perf_counter()
        mass, rain, metrics = solve_stationary(source, sinks, rates)
        metrics["solve_seconds"] = time.perf_counter() - started
        metrics["minimum_column_water_kg_m2"] = float(np.min(mass / area.ravel()))
        metrics["rain_over_land_fraction"] = float(
            rain.reshape(area.shape)[land].sum() / max(rain.sum(), 1e-300)
        )
        outputs[name] = metrics
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-90x180", action="store_true",
                        help="include the climate-density 1.0 prototype-size benchmark")
    args = parser.parse_args()
    base = run_case()
    assert max(x["global_closure_relative"] for x in base.values()) < 1e-8
    assert max(x["max_cell_residual_relative"] for x in base.values()) < 1e-8
    assert min(x["minimum_column_water_kg_m2"] for x in base.values()) >= 0
    assert min(x["rain_over_land_fraction"] for x in base.values()) > 0

    zero = run_case(source_mode="zero")
    assert all(x["minimum_column_water_kg_m2"] == 0 for x in zero.values())
    calm = run_case(wind_mode="calm")
    assert all(x["global_closure_relative"] < 1e-8 for x in calm.values())
    convergence = run_case(wind_mode="convergent")
    divergence = run_case(wind_mode="divergent")
    assert all(x["global_closure_relative"] < 1e-8
               for run in (convergence, divergence) for x in run.values())
    loop = run_case(nlat=2, nlon=8, source_mode="uniform", wind_mode="loop")
    assert all(x["global_closure_relative"] < 1e-8 for x in loop.values())

    opening, rain = cap_rain_to_available(1.0, 1.0 / (9.1 * 86_400.0), 18 * 86_400.0)
    assert opening == 0.0 and rain == 1.0

    latitude, latitude_edges, area = sphere_grid(8, 16)
    land = np.broadcast_to(np.arange(area.shape[1])[None, :] >= 8, area.shape)
    ocean_source = np.zeros(area.shape)
    land_source = np.zeros(area.shape)
    ocean_source[~land] = 3e-5 * area[~land]
    land_source[land] = 1e-5 * area[land]
    sources = ocean_source.ravel() + land_source.ravel()
    winds_u = np.full(area.shape, 5.0)
    winds_v = np.zeros(area.shape)
    source_tags = {}
    orographic = {}
    for name, finite_volume in (("finite_volume", True), ("wind_arrows", False)):
        rates = transition_rates(winds_u, winds_v, latitude, latitude_edges,
                                 area, finite_volume)
        assert all(0 <= source < area.size and 0 <= destination < area.size
                   for source, destination, _ in rates)
        assert all(abs(source // area.shape[1] - destination // area.shape[1]) <= 1
                   for source, destination, _ in rates)
        if name == "finite_volume":
            assert any(source % area.shape[1] == area.shape[1] - 1
                       and destination % area.shape[1] == 0
                       for source, destination, _ in rates)
            assert all(np.isfinite(rate) and rate >= 0 for _, _, rate in rates)
        sinks = np.full(sources.size, 1.0 / (9.1 * 86_400.0))
        total_mass, total_rain, total_metrics = solve_stationary(sources, sinks, rates)
        ocean_mass, ocean_rain, _ = solve_stationary(ocean_source.ravel(), sinks, rates)
        land_mass, _, _ = solve_stationary(land_source.ravel(), sinks, rates)
        tag_error = np.max(np.abs(ocean_mass + land_mass - total_mass)) / np.max(total_mass)
        assert tag_error < 1e-8
        assert np.sum(ocean_rain.reshape(area.shape)[land]) > 0
        source_tags[name] = {
            "mass_tag_sum_relative_error": float(tag_error),
            "ocean_origin_rain_on_land_kg_s": float(
                ocean_rain.reshape(area.shape)[land].sum()
            ),
            "global_closure_relative": total_metrics["global_closure_relative"],
        }

        enhanced_sink = sinks.copy()
        target = np.flatnonzero(land.ravel())[0]
        enhanced_sink[target] *= 3.0
        _, enhanced_rain, enhanced_metrics = solve_stationary(sources, enhanced_sink, rates)
        assert enhanced_rain[target] > total_rain[target]
        assert enhanced_metrics["global_closure_relative"] < 1e-8
        orographic[name] = {
            "target_cell_rain_increase_fraction": float(
                enhanced_rain[target] / total_rain[target] - 1
            ),
            "global_closure_relative": enhanced_metrics["global_closure_relative"],
        }

    latitude, latitude_edges, area = sphere_grid(1, 4)
    rates = transition_rates(np.zeros(area.shape), np.zeros(area.shape), latitude,
                             latitude_edges, area, finite_volume=False)
    try:
        solve_stationary(np.ones(4), np.zeros(4), rates)
    except ValueError as error:
        assert "no rain sink" in str(error)
    else:
        raise AssertionError("sinkless closed class was not rejected")

    resolutions = {
        f"{nlat}x{nlon}": run_case(nlat=nlat, nlon=nlon, source_mode="uniform", wind_mode="calm")
        for nlat, nlon in ((4, 8), (8, 16), (16, 32))
    }
    scale_benchmark = (
        run_case(nlat=90, nlon=180, source_mode="uniform", wind_mode="circulation")
        if args.benchmark_90x180 else None
    )
    result = {"cases": {"ocean_source": base, "zero_source": zero, "calm": calm,
                         "convergent_wind": convergence, "divergent_wind": divergence,
                         "closed_loop": loop},
              "source_tag_linearity": source_tags,
              "enhanced_rain_sink": orographic,
              "finite_interval_rain_cap": {
                  "opening_mass_kg": 1.0,
                  "requested_rain_kg_over_18_days": 18 / 9.1,
                  "actual_rain_kg": rain,
                  "remaining_mass_kg": opening,
              },
              "uniform_source_resolution_sensitivity": resolutions,
              "climate_density_1_scale_benchmark": scale_benchmark,
              "sinkless_closed_class": "detected and rejected",
              "interpretation": "toy-grid analytic screen only; no production climate or performance claim"}
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

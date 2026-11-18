"""Run RFC 0019's focused G2 extension on synthetic grids.

Usage: backend/.venv/bin/python bin/debug/check_rfc0019_coarse_moisture.py [--out evidence.json]
No saved worlds or data downloads are used.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from rfc0019_belt_transport import BELT_CENTERS_DEG, HALF_WIDTH_DEG, broaden_convergence
from rfc0019_coarse_moisture import (
    FaceInterception,
    SECONDS_PER_DAY,
    TransportFaces,
    aggregate_inputs,
    build_transport_rates,
    distribute_rain,
    solve_bounded_budget,
    solve_budget,
)
from rfc0019_transport_prototypes import EARTH_RADIUS_M, solve_stationary, sphere_grid, transition_rates

MIXING = 30_000.0
RAIN_RATE = 1.0 / (9.1 * SECONDS_PER_DAY)


def _solve(grid, faces, source):
    return solve_budget(grid.area, source, np.full(grid.area.size, RAIN_RATE), faces)


def _relative_error(actual, reference):
    return float(np.max(np.abs(actual - reference)) / np.max(np.abs(reference)))


def _residuals(budget):
    return {
        "global_relative": budget.global_relative,
        "max_cell_relative": budget.max_cell_relative,
        "summed_cell_relative": budget.summed_cell_relative,
    }


def _zero_mixing_parity():
    shape = (90, 180)
    latitude, edges, area = sphere_grid(*shape)
    rng = np.random.default_rng(913)
    u = rng.uniform(-5, 5, shape)
    v = rng.uniform(-3, 3, shape)
    ocean = rng.random(shape) > 0.45
    source = ocean * area * rng.uniform(0, 1e-4, shape)
    grid = aggregate_inputs(ocean, {"ocean": source})
    faces = build_transport_rates(grid, u, v, diffusivity=0)
    budget = _solve(grid, faces, grid.tagged_source["ocean"])
    assert not np.any(faces.exchange)
    old_rates = transition_rates(u, v, latitude, edges, area, finite_volume=True)
    old_mass, old_rain, _ = solve_stationary(
        source.ravel(), np.full(area.size, RAIN_RATE), old_rates,
    )
    errors = {
        "mass_relative_to_max_reference": _relative_error(budget.mass, old_mass),
        "rain_relative_to_max_reference": _relative_error(budget.rain, old_rain),
    }
    assert max(errors.values()) < 1e-12
    return errors | _residuals(budget)


def _mixing_and_geometry():
    shape = (4, 8)
    calm = np.zeros(shape)
    grid = aggregate_inputs(np.ones(shape, dtype=bool), {})
    faces = build_transport_rates(grid, calm, calm, MIXING)
    first_row, first_col = np.divmod(faces.first, shape[1])
    second_row, second_col = np.divmod(faces.second, shape[1])
    east = first_row == second_row
    np.testing.assert_array_equal((second_col[east] - first_col[east]) % shape[1], 1)
    np.testing.assert_array_equal(second_col[~east], first_col[~east])
    np.testing.assert_array_equal(second_row[~east] - first_row[~east], 1)
    assert np.count_nonzero(east & (second_col == 0)) == shape[0]
    assert np.count_nonzero(~east) == (shape[0] - 1) * shape[1]

    unequal_area = grid.area.min() < grid.area.max()
    assert unequal_area
    uniform_mass = 10.0 * grid.area
    uniform = _solve(grid, faces, uniform_mass * RAIN_RATE)
    np.testing.assert_allclose(uniform.mass, uniform_mass, rtol=1e-13)
    assert np.max(np.abs(uniform.face_flux)) / uniform.rain.sum() < 1e-13

    donor = grid.native_ids[0, -1]
    receiver = grid.native_ids[0, 0]
    seam = (faces.first == donor) & (faces.second == receiver)
    source = np.zeros_like(grid.area)
    source[donor] = 1.0
    mixed = _solve(grid, faces, source)
    unmixed = _solve(grid, build_transport_rates(grid, calm, calm, 0), source)
    assert mixed.rain[donor] < unmixed.rain[donor]
    assert np.all(mixed.rain > 0)
    assert mixed.face_flux[seam].item() > 0
    zero = _solve(grid, faces, np.zeros_like(source))
    for values in (zero.mass, zero.rain, zero.face_flux):
        np.testing.assert_array_equal(values, 0.0)
    return {"seam_transfer_kg_s": float(mixed.face_flux[seam].item())} | _residuals(mixed)


def _mixing_flux_magnitude():
    shape = (6, 12)
    grid = aggregate_inputs(np.ones(shape, dtype=bool), {})
    calm = np.zeros(shape)
    faces = build_transport_rates(grid, calm, calm, MIXING)
    row, col = 1, 3  # center at 45°S, north face at 30°S
    origin = grid.native_ids[row, col]
    length_angle = np.pi / 6
    checks = {}
    for direction, receiver, conductance in (
        ("east", grid.native_ids[row, col + 1], MIXING / np.cos(np.pi / 4)),
        ("north", grid.native_ids[row + 1, col], MIXING * np.cos(np.pi / 6)),
    ):
        face = (faces.first == origin) & (faces.second == receiver)
        np.testing.assert_allclose(faces.exchange[face], [conductance], rtol=1e-13)
        areas = grid.area[[origin, receiver]]
        np.testing.assert_allclose(areas[0], EARTH_RADIUS_M**2 * length_angle
                                   * (np.sin(-np.pi / 6) - np.sin(-np.pi / 3)))
        one_face = TransportFaces(
            np.array([0]), np.array([1]), np.zeros(1), np.zeros(1), faces.exchange[face],
        )
        budget = solve_budget(areas, np.array([1.0, 0.0]), np.full(2, RAIN_RATE), one_face)
        expected = (conductance / areas[0]) / (RAIN_RATE + np.sum(conductance / areas))
        np.testing.assert_allclose(budget.face_flux, [expected], rtol=1e-13)
        checks[direction] = float(budget.face_flux[0])
    return checks


def _mirrored_wind():
    shape = (6, 12)
    rng = np.random.default_rng(19)
    u, v, native_source = rng.random((3, *shape))
    ocean = rng.random(shape) > 0.5
    grid = aggregate_inputs(ocean, {"test": native_source}, factor=2)
    mirrored = aggregate_inputs(ocean[::-1], {"test": native_source[::-1]}, factor=2)
    direct = _solve(grid, build_transport_rates(grid, u, v, MIXING), grid.tagged_source["test"])
    reflected = _solve(
        mirrored, build_transport_rates(mirrored, u[::-1], -v[::-1], MIXING),
        mirrored.tagged_source["test"],
    )
    placement = distribute_rain(grid, direct.rain, np.ones(shape))
    reversed_placement = distribute_rain(mirrored, reflected.rain, np.ones(shape))[::-1]
    np.testing.assert_allclose(placement, reversed_placement, rtol=1e-13)
    return {"rain_relative_to_max_reference": _relative_error(reversed_placement, placement)}


def _mixed_island():
    shape = (8, 16)
    ocean = np.ones(shape, dtype=bool)
    ocean[0, 0] = False
    source = np.zeros(shape)
    source[0, 1] = 1
    grid = aggregate_inputs(ocean, {"ocean": source}, factor=4)
    calm = np.zeros(shape)
    faces = build_transport_rates(grid, calm, calm, diffusivity=0)
    rain = _solve(grid, faces, grid.tagged_source["ocean"])
    ocean_id, land_id = grid.native_ids[0, 1], grid.native_ids[0, 0]
    cross = grid.is_ocean[faces.first] != grid.is_ocean[faces.second]
    ocean_to_land = np.where(grid.is_ocean[faces.first], rain.face_flux, -rain.face_flux)[cross]
    gross_import = float(np.maximum(ocean_to_land, 0).sum())
    gross_export = float(np.maximum(-ocean_to_land, 0).sum())
    land_rain = rain.rain[land_id]
    assert grid.is_ocean.size == 9 and np.count_nonzero(~grid.is_ocean) == 1
    ocean_area, land_area = grid.area[[ocean_id, land_id]]
    conductance = ocean_area * land_area / (ocean_area + land_area) / (2 * SECONDS_PER_DAY)
    a, b = conductance / ocean_area, conductance / land_area
    np.testing.assert_allclose(land_rain, a / (RAIN_RATE + a + b), rtol=1e-13)
    np.testing.assert_allclose(land_rain, gross_import - gross_export, rtol=1e-13)
    weights = np.linspace(0.5, 1.5, ocean.size).reshape(shape)
    fine = distribute_rain(grid, rain.rain, weights)
    by_reservoir = np.bincount(grid.native_ids.ravel(), weights=fine.ravel())
    np.testing.assert_allclose(by_reservoir, rain.rain, rtol=1e-13)
    np.testing.assert_allclose(fine[~ocean].sum(), land_rain, rtol=1e-13)
    np.testing.assert_array_equal(fine[grid.coarse_cell[grid.native_ids] != 0], 0.0)
    ocean_density = fine[0, 1:4] / grid.native_area[0, 1:4]
    np.testing.assert_allclose(ocean_density / ocean_density[0], weights[0, 1:4] / weights[0, 1])
    return {
        "land_rain_kg_s": float(land_rain),
        "gross_ocean_to_land_kg_s": gross_import,
        "gross_land_to_ocean_kg_s": gross_export,
        "max_native_allocation_error_kg_s": float(np.max(np.abs(by_reservoir - rain.rain))),
    } | _residuals(rain)


def _internal_shoreline_wind():
    shape = (8, 16)
    ocean = np.ones(shape, dtype=bool)
    ocean[0, 0] = False
    grid = aggregate_inputs(ocean, {}, factor=4)
    faces = build_transport_rates(grid, np.full(shape, 2.0), np.full(shape, 3.0), MIXING)
    land_id, ocean_id = grid.native_ids[0, :2]
    _, edges, _ = sphere_grid(*shape)

    # Only the east and north sides of this island are inside its footprint.
    # The west side wraps into a different coarse cell and must stay separate.
    shore = (faces.first == land_id) & (faces.second == ocean_id) & (faces.forward > 0)
    expected = [
        2.0 * EARTH_RADIUS_M * (edges[1] - edges[0]),
        3.0 * EARTH_RADIUS_M * np.cos(edges[1]) * (2 * np.pi / shape[1]),
    ]
    np.testing.assert_allclose(np.sort(faces.forward[shore]), np.sort(expected), rtol=1e-14)
    np.testing.assert_array_equal(faces.exchange[shore], 0.0)
    np.testing.assert_array_equal(faces.reverse[shore], 0.0)
    shared = (faces.first == land_id) & (faces.second == ocean_id)
    areas = grid.area[[land_id, ocean_id]]
    isolated = TransportFaces(
        np.zeros(np.count_nonzero(shared), dtype=int),
        np.ones(np.count_nonzero(shared), dtype=int),
        faces.forward[shared], faces.reverse[shared], faces.exchange[shared],
    )
    budget = solve_budget(areas, np.array([1., 0.]), np.full(2, RAIN_RATE), isolated)
    exchange = np.prod(areas) / (areas.sum() * 2 * SECONDS_PER_DAY)
    land_to_ocean = (sum(expected) + exchange) / areas[0]
    ocean_to_land = exchange / areas[1]
    expected_rain = land_to_ocean / (RAIN_RATE + land_to_ocean + ocean_to_land)
    np.testing.assert_allclose(budget.rain[1], expected_rain, rtol=1e-13)
    return {
        "internal_wind_faces": int(np.count_nonzero(shore)),
        "wind_and_exchange_rain_kg_s": float(budget.rain[1]),
    }


def _opposing_flow_and_tags():
    shape = (4, 8)
    ocean = np.indices(shape)[1] < 4
    ocean_source = ocean.astype(float)
    land_source = (~ocean).astype(float) * 0.3
    grid = aggregate_inputs(ocean, {"ocean": ocean_source, "land": land_source}, factor=2)
    u = np.where(np.indices(shape)[0] % 2 == 0, 4.0, -4.0)
    faces = build_transport_rates(grid, u, np.zeros(shape), MIXING)
    ocean_id, land_id = grid.native_ids[0, 3], grid.native_ids[0, 4]
    pair = (faces.first == ocean_id) & (faces.second == land_id)
    assert np.count_nonzero(pair) == 2
    assert faces.forward[pair].sum() > 0 and faces.reverse[pair].sum() > 0
    from_ocean = _solve(grid, faces, grid.tagged_source["ocean"])
    from_land = _solve(grid, faces, grid.tagged_source["land"])
    both = _solve(grid, faces, grid.tagged_source["ocean"] + grid.tagged_source["land"])
    np.testing.assert_allclose(from_ocean.rain + from_land.rain, both.rain, rtol=1e-13)
    np.testing.assert_allclose(from_ocean.mass + from_land.mass, both.mass, rtol=1e-13)
    cross = grid.is_ocean[faces.first] != grid.is_ocean[faces.second]
    flux = np.where(grid.is_ocean[faces.first], from_ocean.face_flux, -from_ocean.face_flux)[cross]
    gross_in = float(np.maximum(flux, 0).sum())
    gross_out = float(np.maximum(-flux, 0).sum())
    np.testing.assert_allclose(from_ocean.rain[~grid.is_ocean].sum(), gross_in - gross_out, rtol=1e-13)
    assert gross_in > 0 and gross_out > 0
    return {
        "source_tag_rain_relative_to_max_reference": _relative_error(from_ocean.rain + from_land.rain, both.rain),
        "ocean_tag_gross_import_kg_s": gross_in,
        "ocean_tag_gross_return_kg_s": gross_out,
    } | _residuals(both)


def _aligned_resolution():
    shape = (4, 8)
    ocean = np.indices(shape)[1] < 4
    rain_by_factor = {}
    for factor in (1, 4):
        mask = ocean.repeat(factor, axis=0).repeat(factor, axis=1)
        _, _, area = sphere_grid(*mask.shape)
        source = area * mask * 1e-6
        grid = aggregate_inputs(mask, {"ocean": source}, factor=factor)
        u, v = np.full(mask.shape, 3.0), np.full(mask.shape, -0.3)
        faces = build_transport_rates(grid, u, v, MIXING)
        rain_by_factor[factor] = _solve(grid, faces, grid.tagged_source["ocean"]).rain
    np.testing.assert_allclose(rain_by_factor[4], rain_by_factor[1], rtol=1e-13)
    return {"rain_relative_to_max_reference": _relative_error(rain_by_factor[4], rain_by_factor[1])}


def _targeted_belt_faces():
    shape = (90, 180)
    _, edges, _ = sphere_grid(*shape)
    centers = (edges[:-1] + edges[1:]) / 2
    v = -4 * np.sign(centers)[:, None] + np.linspace(-1, 1, shape[1])
    delta, geometry = broaden_convergence(v)
    face_latitudes = np.linspace(-90.0, 90.0, shape[0] + 1)[1:-1]
    interior = np.any(
        [np.abs(face_latitudes - center) < HALF_WIDTH_DEG for center in BELT_CENTERS_DEG],
        axis=0,
    )
    np.testing.assert_array_equal(delta[~interior], 0.0)
    grid = aggregate_inputs(np.zeros(shape, dtype=bool), {})
    u = np.zeros(shape)
    control = build_transport_rates(grid, u, v, 0.0)
    changed = build_transport_rates(grid, u, v, 0.0, north_face_delta=delta)
    nlat, nlon = shape
    np.testing.assert_array_equal(changed.first, control.first)
    np.testing.assert_array_equal(changed.second, control.second)
    capacity = changed.forward - changed.reverse - (control.forward - control.reverse)
    face_length = EARTH_RADIUS_M * np.cos(edges[1:-1]) * (2 * np.pi / nlon)
    np.testing.assert_array_equal(capacity[:nlat * nlon], 0.0)
    np.testing.assert_allclose(
        capacity[nlat * nlon:].reshape(nlat - 1, nlon),
        np.broadcast_to((face_length * delta)[:, None], (nlat - 1, nlon)),
        atol=1e-9, rtol=1e-12,
    )
    equator = geometry[1]
    assert equator["peak_convergence_after_m2_s"] < equator["peak_convergence_before_m2_s"]
    return {"changed_latitude_faces": int(np.count_nonzero(delta))}


def _bounded_columns():
    area = np.array([2.0, 1.0, 1.0])
    source = np.array([9.0, 0.0, 0.0])
    sink = np.ones_like(area)
    faces = TransportFaces(np.array([0, 1]), np.array([1, 2]), area[:2], np.zeros(2), np.zeros(2))
    # The middle column releases first. The terminal column can no longer
    # sustain its cap once the middle column exports less.
    cap = np.array([2.0, 2.0, 1.5])
    budget, overflow, solves = solve_bounded_budget(area, source, sink, faces, cap)
    np.testing.assert_allclose(budget.mass, [2.0, 1.0, 1.0])
    np.testing.assert_allclose(budget.face_flux, [2.0, 1.0])
    np.testing.assert_allclose(overflow, [5.0, 0.0, 0.0])
    np.testing.assert_allclose(budget.rain, [7.0, 1.0, 1.0])
    assert solves == 2
    for trial_source in (source, np.zeros_like(source)):
        expected = solve_budget(area, trial_source, sink, faces)
        uncapped, extra, _ = solve_bounded_budget(area, trial_source, sink, faces, cap + 10)
        np.testing.assert_allclose(uncapped.mass, expected.mass)
        np.testing.assert_allclose(uncapped.rain, expected.rain)
        np.testing.assert_array_equal(extra, 0.0)
    return {"free_solves": solves, "overflow_kg_s": float(overflow.sum())} | _residuals(budget)


def _intercepted_ridge():
    area = np.ones(3)
    source = np.array([20.0, 0.0, 0.0])
    sink = np.ones(3)
    cap = np.array([6.0, 1.0, 100.0])
    first, second = np.array([0, 1]), np.array([1, 2])
    forward, reverse, exchange = np.ones(2), np.zeros(2), np.zeros(2)
    interception = FaceInterception(np.array([0.5, 0.0]), np.zeros(2))
    faces = TransportFaces(first, second, forward, reverse, exchange)
    budget, overflow, _ = solve_bounded_budget(area, source, sink, faces, cap, interception=interception)
    np.testing.assert_allclose(budget.mass, [6.0, 1.0, 1.0])
    np.testing.assert_allclose(budget.face_flux, [6.0, 1.0])
    np.testing.assert_allclose(budget.intercepted_rain_kg_s, [0.0, 3.0, 0.0])
    np.testing.assert_allclose(overflow, [8.0, 1.0, 0.0])
    np.testing.assert_allclose(budget.rain, [14.0, 5.0, 1.0])
    mirrored = TransportFaces(2 - second, 2 - first, reverse, forward, exchange)
    mirrored_interception = FaceInterception(interception.reverse, interception.forward)
    reverse_budget, _, _ = solve_bounded_budget(
        area, source[::-1], sink, mirrored, cap[::-1], interception=mirrored_interception,
    )
    np.testing.assert_allclose(reverse_budget.rain, budget.rain[::-1])
    return _residuals(budget)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = {
        "zero_mixing_90x180": _zero_mixing_parity(),
        "mixing_and_geometry": _mixing_and_geometry(),
        "analytic_mixing_flux_kg_s": _mixing_flux_magnitude(),
        "hemisphere_mirror": _mirrored_wind(),
        "mixed_coastal_island": _mixed_island(),
        "internal_shoreline_wind": _internal_shoreline_wind(),
        "opposing_flow_and_source_tags": _opposing_flow_and_tags(),
        "aligned_resolution": _aligned_resolution(),
        "targeted_belt_faces": _targeted_belt_faces(),
        "bounded_columns": _bounded_columns(),
        "intercepted_ridge": _intercepted_ridge(),
    }
    report = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()

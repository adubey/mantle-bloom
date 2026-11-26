"""Thin Python bridge to the `mantle_bloom_rust` PyO3 acceleration extension (#1036)."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix

_DISABLE_RUST = os.environ.get("MANTLE_BLOOM_DISABLE_RUST") == "1"

if _DISABLE_RUST:
    _rust = None
else:
    try:
        import mantle_bloom_rust as _rust
    except ImportError:
        _rust = None

RUST_AVAILABLE: bool = _rust is not None


def _f64(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.float64)


def _bool(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.bool_)


def _i64(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.int64)


def build_march_graph(graph: "csr_matrix", step_m: float, neighbor_idx: np.ndarray | None = None):
    """Build a persistent Rust `MarchGraph` from a symmetric CSR graph and optional neighbor table."""
    if _rust is None:
        return None
    if not graph.has_canonical_format:
        graph = graph.copy()
        graph.sum_duplicates()
    nb = None if neighbor_idx is None else _i64(neighbor_idx)
    return _rust.MarchGraph(
        _i64(graph.indptr),
        _i64(graph.indices),
        _f64(graph.data),
        float(step_m),
        nb,
    )


def step_ice(
    rust_graph,
    area_m2: np.ndarray,
    elevation: np.ndarray,
    prior_deflection_m: np.ndarray,
    sea_level_m: float,
    prev_depth: np.ndarray,
    climate,
    mobile_cover_m: np.ndarray,
    is_ocean: np.ndarray,
    years: float,
    margin_hint_m: np.ndarray | None,
    dijkstra_limit_m: float | None,
    *,
    body_cls,
    ledger_cls,
    frozen_cls,
    ice_step_cls,
):
    """Run `mantle_bloom_rust.step_ice` and wrap outputs into Python `IceStep` dataclasses."""
    (
        depth,
        melt_m,
        thickness,
        covered,
        is_warm,
        flux,
        sliding,
        calved,
        ledger_tuple,
        frozen_tuple,
        body_tuples,
    ) = _rust.step_ice(
        rust_graph,
        _f64(area_m2),
        _f64(elevation),
        _f64(prior_deflection_m),
        float(sea_level_m),
        _f64(prev_depth),
        _f64(climate.precipitation_mm),
        _f64(climate.temperature_c),
        _f64(climate.seasonal_amplitude_c),
        _f64(mobile_cover_m),
        _bool(is_ocean),
        float(years),
        None if climate.frictional_heat_w_m2 is None else _f64(climate.frictional_heat_w_m2),
        None if climate.dry_snow_mm is None else _f64(climate.dry_snow_mm),
        None if climate.bed_frozen_snow_mm is None else _f64(climate.bed_frozen_snow_mm),
        None if margin_hint_m is None else _f64(margin_hint_m),
        None if dijkstra_limit_m is None else float(dijkstra_limit_m),
    )

    bodies = [
        body_cls(
            nodes=nodes,
            opening_m3=opening_m3,
            accumulation_m3=acc_m3,
            surface_melt_m3=sm_m3,
            basal_melt_m3=bm_m3,
            calved_m3=calved_m3,
            overflow_m3=overflow_m3,
            closing_m3=closing_m3,
            capped_m3=capped_m3,
            margin_m=margin_m,
            saturated=saturated,
            domain_nodes=domain_nodes,
        )
        for (
            nodes,
            domain_nodes,
            opening_m3,
            acc_m3,
            sm_m3,
            bm_m3,
            calved_m3,
            overflow_m3,
            closing_m3,
            capped_m3,
            margin_m,
            saturated,
        ) in body_tuples
    ]
    ledger = ledger_cls(
        opening_m3=ledger_tuple[0],
        accumulation_m3=ledger_tuple[1],
        surface_melt_m3=ledger_tuple[2],
        basal_melt_m3=ledger_tuple[3],
        calved_m3=ledger_tuple[4],
        overflow_m3=ledger_tuple[5],
        closing_m3=ledger_tuple[6],
        snow_melted_m3=ledger_tuple[7],
        capped_m3=ledger_tuple[8],
        capped_nodes=ledger_tuple[9],
    )
    frozen_precipitation = frozen_cls(
        land_m3=frozen_tuple[0],
        ocean_retained_m3=frozen_tuple[1],
        retained_m3=frozen_tuple[2],
        bare_melted_m3=frozen_tuple[3],
    )
    return depth, melt_m, ice_step_cls(
        thickness=thickness,
        covered=covered,
        warm_based=is_warm,
        balance_flux_m3_per_yr=flux,
        sliding_speed_m_per_yr=sliding,
        calved_m3=calved,
        ledger=ledger,
        bodies=bodies,
        frozen_precipitation=frozen_precipitation,
    )


def seeded_flood_levels(
    elevations: np.ndarray,
    neighbor_idx: np.ndarray,
    seed: np.ndarray,
    adjacency: tuple[np.ndarray, np.ndarray] | None = None,
    rust_graph=None,
) -> np.ndarray:
    """Compute seeded flood levels via Rust."""
    if rust_graph is not None:
        return _rust.seeded_flood_levels(_f64(elevations), _bool(seed), rust_graph=rust_graph)
    if adjacency is not None:
        return _rust.seeded_flood_levels(
            _f64(elevations),
            _bool(seed),
            indptr=_i64(adjacency[0]),
            adjacency=_i64(adjacency[1]),
        )
    return _rust.seeded_flood_levels(_f64(elevations), _bool(seed), neighbor_idx=_i64(neighbor_idx))


def solve_sea_level_connected(
    elevations: np.ndarray,
    areas_m2: np.ndarray,
    water_volume_m3: float,
    seed_is_ocean: np.ndarray,
    seed_neighbor_idx: np.ndarray,
    *,
    bracket: tuple[float, float] | None = None,
    seed_adjacency: tuple[np.ndarray, np.ndarray] | None = None,
    ice_water_m: np.ndarray | None = None,
    rust_graph=None,
) -> float:
    """Solve connected-ocean sea level via Rust fused priority-flood and bisection."""
    iw = None if ice_water_m is None else _f64(ice_water_m)
    br = None if bracket is None else (float(bracket[0]), float(bracket[1]))
    if rust_graph is not None:
        return float(
            _rust.solve_sea_level_connected(
                _f64(elevations),
                _f64(areas_m2),
                float(water_volume_m3),
                _bool(seed_is_ocean),
                ice_water_m=iw,
                bracket=br,
                rust_graph=rust_graph,
            )
        )
    if seed_adjacency is not None:
        return float(
            _rust.solve_sea_level_connected(
                _f64(elevations),
                _f64(areas_m2),
                float(water_volume_m3),
                _bool(seed_is_ocean),
                ice_water_m=iw,
                bracket=br,
                indptr=_i64(seed_adjacency[0]),
                adjacency=_i64(seed_adjacency[1]),
            )
        )
    return float(
        _rust.solve_sea_level_connected(
            _f64(elevations),
            _f64(areas_m2),
            float(water_volume_m3),
            _bool(seed_is_ocean),
            ice_water_m=iw,
            bracket=br,
            neighbor_idx=_i64(seed_neighbor_idx),
        )
    )


def below_sea_components(below: np.ndarray, neighbor_idx: np.ndarray, rust_graph=None) -> np.ndarray:
    """Label connected components of below-sea-level nodes via Rust union-find."""
    if rust_graph is not None:
        return _rust.below_sea_components(_bool(below), rust_graph=rust_graph)
    return _rust.below_sea_components(_bool(below), neighbor_idx=_i64(neighbor_idx))


def connected_ocean_mask(
    elevation: np.ndarray,
    sea_level_m: float,
    neighbor_idx: np.ndarray,
    *,
    areas_m2: np.ndarray | None = None,
    prior_ocean: np.ndarray | None = None,
    enclosed_min_area_km2: float,
    enclosed_keep_area_km2: float,
    planet_radius_km: float,
    rust_graph=None,
) -> np.ndarray:
    """Classify connected ocean + qualifying enclosed oceans via Rust union-find."""
    areas = None if areas_m2 is None else _f64(areas_m2)
    prior = None if prior_ocean is None else _bool(prior_ocean)
    if rust_graph is not None:
        return _rust.connected_ocean_mask(
            _f64(elevation),
            float(sea_level_m),
            float(enclosed_min_area_km2),
            float(enclosed_keep_area_km2),
            float(planet_radius_km),
            areas_m2=areas,
            prior_ocean=prior,
            rust_graph=rust_graph,
        )
    return _rust.connected_ocean_mask(
        _f64(elevation),
        float(sea_level_m),
        float(enclosed_min_area_km2),
        float(enclosed_keep_area_km2),
        float(planet_radius_km),
        neighbor_idx=_i64(neighbor_idx),
        areas_m2=areas,
        prior_ocean=prior,
    )


def march_substep_prep(
    working_elevation: np.ndarray,
    deflection_m: np.ndarray,
    ref_depth: np.ndarray,
    sea_level: float,
    ocean_mask: np.ndarray,
    crust_h_m: np.ndarray,
    mantle_h_m: np.ndarray,
    rho_crust: np.ndarray,
    temperature_at_nodes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute lapse-adjusted air temperature and frozen mask for a march substep in Rust."""
    return _rust.march_substep_prep(
        _f64(working_elevation),
        _f64(deflection_m),
        _f64(ref_depth),
        float(sea_level),
        _bool(ocean_mask),
        _f64(crust_h_m),
        _f64(mantle_h_m),
        _f64(rho_crust),
        _f64(temperature_at_nodes),
    )


def apply_march_abrasion(
    sliding_speed_m_per_yr: np.ndarray,
    dt: float,
    multiplier: float,
    ocean_mask: np.ndarray,
    cover_remaining: np.ndarray,
    substrate_headroom: np.ndarray,
    craton_keep: np.ndarray,
    raw_abrasion_depth: np.ndarray,
    abrasion_depth: np.ndarray,
    working_elevation: np.ndarray,
) -> float:
    """Book one interval's glacial abrasion in place in Rust; returns max bed cut (m)."""
    return float(
        _rust.apply_march_abrasion(
            _f64(sliding_speed_m_per_yr),
            float(dt),
            float(multiplier),
            _bool(ocean_mask),
            cover_remaining,
            substrate_headroom,
            _f64(craton_keep),
            raw_abrasion_depth,
            abrasion_depth,
            working_elevation,
        )
    )

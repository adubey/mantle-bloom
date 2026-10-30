#!/usr/bin/env python3
"""RFC 0014 gate G5: offline 10 Myr lake-ice-cap budget prototype (#541).

Uses the preregistered #538 diagnosis saves only. The cap recurrence is a diagnostic
implementation of D2(c), not production code: a Stefan conduction increment in physical ice
thickness, converted to water equivalent, minus the existing glacier positive-degree melt law
and a constant sublimation flux. Physical ice is capped at 2 m; the water-equivalent cap is
also limited by available lake water. Climate and bathymetry are held fixed for this screen.

Usage (one save at a time; do not run beside a replay or test):

    backend/.venv/bin/python bin/debug/gate_lake_ice_cap.py \
        ~/Downloads/mantle-bloom-seed910211954-220200000y.mbworld \
        --out analysis/issue541/seed910211954.json

The saved JSON contains all sublimation (0.25x/1x/4x) by conduction (0.5x/1x/2x)
sensitivity runs and the input save checksum. Run once for each of the three diagnosis saves.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "backend"))

RHO_ICE_KG_M3 = 917.0
RHO_WATER_KG_M3 = 1000.0
LATENT_HEAT_J_KG = 334_000.0
ICE_CONDUCTIVITY_W_M_K = 2.1
SECONDS_PER_DAY = 86400.0
SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY
YEARS = 10_000_000.0
STEP_YEARS = 100_000.0
MAX_PHYSICAL_ICE_M = 2.0
SUBLIMATION_SCALES = (0.25, 1.0, 4.0)
CONDUCTION_SCALES = (0.5, 1.0, 2.0)
DIAGNOSIS_SAVES = {
    "bd4874092efc805e8b174e115fe567c82d35fc2ee606686938b5b6eb96b060d2": 910211954,
    "ec27d5077bdf43b3a6f0339b5726fa4d110da9eb1f143f96a3004d79f309b045": 635056485,
    "6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd": 317739385,
}


def _stats(values: np.ndarray) -> dict:
    return {
        "p50_m": round(float(np.percentile(values, 50)), 6),
        "p90_m": round(float(np.percentile(values, 90)), 6),
        "p99_m": round(float(np.percentile(values, 99)), 6),
        "max_m": round(float(np.max(values, initial=0.0)), 6),
    }


def _simulate(depth_m: np.ndarray, area_m2: np.ndarray, fdd_c_days: np.ndarray,
              pdd_c_days: np.ndarray, ice_fraction: np.ndarray, base_sublimation_m_yr: float,
              melt_m_per_myr_per_20_c_days: float) -> dict:
    from app.hydroclimate import hydrology

    wet_area = float(area_m2.sum())
    conductivity = ICE_CONDUCTIVITY_W_M_K
    # d(h_i^2)/dt = 2*k*FDD/(rho_i*L_f); a fixed annual loss makes
    # d(h_i)/dt = q/h_i - loss_i, with equilibrium h_i = q/loss_i.
    q_m2_yr = conductivity * fdd_c_days * SECONDS_PER_DAY / (RHO_ICE_KG_M3 * LATENT_HEAT_J_KG)
    pdd_c_years = pdd_c_days / 365.25
    melt_i_m_yr = (melt_m_per_myr_per_20_c_days / 1_000_000.0
                   * np.clip(pdd_c_years / hydrology.GLACIER_MELT_REFERENCE_DEGREES_C,
                             0.0, hydrology.GLACIER_MELT_MAX_FACTOR))
    sublimation_we_m_yr = base_sublimation_m_yr * ice_fraction
    loss_i_m_yr = melt_i_m_yr + sublimation_we_m_yr * RHO_WATER_KG_M3 / RHO_ICE_KG_M3
    equilibrium_i_m = np.divide(q_m2_yr, loss_i_m_yr, out=np.zeros_like(q_m2_yr), where=loss_i_m_yr > 0.0)
    bed_i_m = depth_m * RHO_WATER_KG_M3 / RHO_ICE_KG_M3
    max_i_m = np.minimum(MAX_PHYSICAL_ICE_M, bed_i_m)
    shallow = depth_m <= MAX_PHYSICAL_ICE_M
    deep_area_m2 = float(area_m2[~shallow].sum())
    snapshots = []
    for t_years in (0.0, 100_000.0, 500_000.0, 1_000_000.0, 2_500_000.0,
                    5_000_000.0, 7_500_000.0, YEARS):
        # Invert the exact time integral of dh/dt=q/h-loss by vectorized bisection. This
        # avoids timestep instability when the initial sublimation flux is much larger than
        # a 100-kyr Stefan increment, while retaining the transient and its plateau.
        cap_i_m = np.zeros_like(depth_m)
        active = (equilibrium_i_m > 0.0) & (t_years > 0.0)
        low = np.zeros_like(depth_m)
        high = np.minimum(equilibrium_i_m, max_i_m)
        target = loss_i_m_yr * t_years
        for _ in range(48):
            mid = (low + high) * 0.5
            gap = np.maximum(equilibrium_i_m - mid, np.finfo(float).tiny)
            time_integral = equilibrium_i_m * np.log(
                np.maximum(equilibrium_i_m, np.finfo(float).tiny) / gap
            ) - mid
            below_target = time_integral < target
            low = np.where(active & below_target, mid, low)
            high = np.where(active & ~below_target, mid, high)
        cap_i_m = np.where(active, (low + high) * 0.5, cap_i_m)
        # The owner-directed non-permafrost rule limits surface ice to 2 m; shallow
        # water may still freeze to the bed. Keep the physical thickness in metres.
        cap_i_m = np.where(active, np.minimum(cap_i_m, max_i_m), 0.0)
        cap_we_m = cap_i_m * RHO_ICE_KG_M3 / RHO_WATER_KG_M3
        bedfast = cap_i_m >= bed_i_m - 1e-9
        volume_m3 = float(np.dot(cap_we_m, area_m2))
        snapshots.append({
            "time_myr": round(t_years / 1_000_000.0, 2),
            "cap_volume_km3": round(volume_m3 / 1e9, 6),
            "wet_area_ice_covered_pct": round(100.0 * float(area_m2[cap_we_m > 0].sum()) / wet_area, 6),
            "frozen_to_bed_wet_area_pct": round(100.0 * float(area_m2[bedfast].sum()) / wet_area, 6),
            "shallow_bedfast_wet_area_pct": round(100.0 * float(area_m2[bedfast & shallow].sum()) / wet_area, 6),
            "deep_bedfast_deep_area_pct": round(100.0 * float(area_m2[bedfast & ~shallow].sum()) / deep_area_m2, 6) if deep_area_m2 else 0.0,
        })

    wet = cap_we_m > 0.0
    volume_m3 = float(np.dot(cap_we_m, area_m2))
    sublimation_m3_yr = float(np.dot(np.where(wet, sublimation_we_m_yr, 0.0), area_m2))
    return {
        "final_thickness_we": _stats(cap_we_m),
        "final_physical_ice_thickness": _stats(cap_we_m * RHO_WATER_KG_M3 / RHO_ICE_KG_M3),
        "wet_cells": int(len(depth_m)),
        "ice_covered_cells": int(wet.sum()),
        "frozen_to_bed_wet_area_pct": round(100.0 * float(area_m2[cap_i_m >= bed_i_m - 1e-9].sum()) / wet_area, 6),
        "shallow_bedfast_wet_area_pct": round(100.0 * float(area_m2[(cap_i_m >= bed_i_m - 1e-9) & shallow].sum()) / wet_area, 6),
        "deep_bedfast_deep_area_pct": round(100.0 * float(area_m2[(cap_i_m >= bed_i_m - 1e-9) & ~shallow].sum()) / deep_area_m2, 6) if deep_area_m2 else 0.0,
        "max_physical_ice_m": MAX_PHYSICAL_ICE_M,
        "cap_volume_km3": round(volume_m3 / 1e9, 6),
        "global_sublimation_km3_yr": round(sublimation_m3_yr / 1e9, 9),
        "volume_trajectory": snapshots,
    }


def measure(save: Path) -> dict:
    from app import persistence, plates
    from app.elevation_lines import line_spacing_rad
    from app.geometry import xyz_to_latlon
    from app.hydroclimate import biomes, climate, erosion, hydrology

    raw = save.read_bytes()
    save_sha256 = hashlib.sha256(raw).hexdigest()
    if save_sha256 not in DIAGNOSIS_SAVES:
        raise ValueError("G5 accepts only the three preregistered RFC 0014 diagnosis saves")
    world = persistence.load_world_bytes(raw)
    if world.seed != DIAGNOSIS_SAVES[save_sha256]:
        raise ValueError("save checksum is registered for a different seed")
    fields = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density))
    points, plate_order = plates.gather_node_positions(world.plates)
    if len(points) == 0:
        raise ValueError("save has no live nodes")
    depth = plates.collect_all_lake_depth(plate_order)
    elevation = plates.collect_all_elevation(plate_order)
    areas = plates.collect_all_accounting_areas_m2(plate_order, line_spacing_rad(world.node_density))
    lat_rad, _ = xyz_to_latlon(points)
    row, col = erosion.climate_grid_indices(points, *fields.air_temperature_c.shape)
    ocean = fields.is_ocean[row, col]
    wet = (depth > 1e-6) & ~ocean & (areas > 0)
    if not np.any(wet):
        raise ValueError("diagnosis save contains no freshwater lake nodes")
    depth = depth[wet]
    area = areas[wet]
    row, col = row[wet], col[wet]
    lat = np.degrees(lat_rad[wet])
    surface_elevation = elevation[wet] + depth
    mean_c = fields.air_temperature_c[row, col]
    continentality = biomes.grid_continentality(fields.is_ocean)[row, col]
    amp_c = biomes.seasonal_temp_amplitude(
        lat, continentality, world.axial_tilt_deg, relief_m=surface_elevation - world.sea_level_m
    )
    fdd = hydrology.seasonal_freezing_degrees(mean_c, amp_c) * 365.25
    pdd = hydrology.seasonal_positive_degrees(mean_c, amp_c) * 365.25
    ice_fraction = 1.0 - hydrology.seasonal_fraction_above(mean_c, amp_c)

    runs = []
    for sublimation_scale in SUBLIMATION_SCALES:
        for conduction_scale in CONDUCTION_SCALES:
            # Scale Stefan's conductivity, retaining all other physical constants.
            # Scaling FDD is algebraically equivalent to scaling Stefan conductivity.
            sim = _simulate(depth, area, fdd * conduction_scale, pdd, ice_fraction,
                            hydrology.LAKE_ICE_SUBLIMATION_M_PER_YR * sublimation_scale,
                            hydrology.GLACIER_MELT_RATE_M_PER_MYR)
            sim["sublimation_scale"] = sublimation_scale
            sim["conduction_coefficient_scale"] = conduction_scale
            runs.append(sim)

    return {
        "issue": 541,
        "rfc": "0014",
        "gate": "G5",
        "save": save.name,
        "save_sha256": save_sha256,
        "seed": getattr(world, "seed", None),
        "age_myr": round(float(world.elapsed_years) / 1e6, 3),
        "axial_tilt_deg": float(world.axial_tilt_deg),
        "climate_density": float(world.climate_density),
        "lake_nodes": int(len(depth)),
        "wet_lake_area_km2": round(float(area.sum()) / 1e6, 3),
        "shallow_lake_area_km2": round(float(area[depth <= MAX_PHYSICAL_ICE_M].sum()) / 1e6, 3),
        "deep_lake_area_km2": round(float(area[depth > MAX_PHYSICAL_ICE_M].sum()) / 1e6, 3),
        "physical_ice_cap_m": MAX_PHYSICAL_ICE_M,
        "permafrost_assumption": "not represented; all modeled lakes treated as non-permafrost",
        "lake_depth_m": _stats(depth),
        "temperature_c": {"mean_p50": round(float(np.percentile(mean_c, 50)), 3),
                          "amplitude_p50": round(float(np.percentile(amp_c, 50)), 3),
                          "fdd_c_days_p50": round(float(np.percentile(fdd, 50)), 3),
                          "pdd_c_days_p50": round(float(np.percentile(pdd, 50)), 3)},
        "method": {
            "climate_and_depth": "held fixed at saved state; sampled on live nodes",
            "water_surface": "air temperature and seasonal amplitude sampled at lake nodes; amplitude depth damping off",
            "growth": "Stefan ODE d(h_i^2)/dt=2*k*FDD/(rho_i*L_f), integrated exactly against constant annual losses; physical ice capped at min(2 m, local lake depth); water-equivalent conversion 917/1000",
            "conductivity": "k=2.1 W m-1 K-1 baseline; scale factors 0.5, 1, 2",
            "melt": "existing glacier positive-degree law: 400 m/Myr at 20 annual-mean positive C, capped at 3x; converted to water equivalent",
            "sublimation": "0.6 m water equivalent per ice-covered year, multiplied by 1-f_free for each cell; scale factors 0.25, 1, 4",
            "step_years": "closed-form integration; checkpoints at 0, 0.1, 0.5, 1, 2.5, 5, 7.5, 10 Myr",
            "duration_years": YEARS,
            "limits": "No climate, lake-depth, or water-budget feedback; non-permafrost 2 m physical cap is an owner-directed empirical bound.",
        },
        "runs": runs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("save", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = measure(args.save)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: result[k] for k in ("seed", "age_myr", "lake_nodes", "wet_lake_area_km2")}))
    for run in result["runs"]:
        print(json.dumps({k: run[k] for k in ("sublimation_scale", "conduction_coefficient_scale", "cap_volume_km3", "frozen_to_bed_wet_area_pct", "shallow_bedfast_wet_area_pct", "deep_bedfast_deep_area_pct", "global_sublimation_km3_yr")}))


if __name__ == "__main__":
    main()

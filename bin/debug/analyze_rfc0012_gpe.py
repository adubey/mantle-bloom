"""Reproduce the RFC 0012 G2 column-GPE magnitude screen.

Usage: backend/.venv/bin/python bin/debug/analyze_rfc0012_gpe.py [--output PATH]

The script evaluates geometry-defined first-order ridge and India–Tibet columns,
prints the component loads, and optionally writes compact JSON. This is a scale
check, not a calibrated thermal/isostatic inversion.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

G_M_S2 = 9.81


def column_gpe(layers: list[tuple[float, float, float]], base_z_m: float) -> float:
    """Integrate rho*g*(height above the common basal datum) over layers.

    Each layer is (top_z, bottom_z, density); z is positive upward. Return J/m²,
    numerically equal to N/m. The common datum cancels for mass-balanced columns.
    """
    return sum(
        0.5 * rho * G_M_S2 * ((top - base_z_m) ** 2 - (bottom - base_z_m) ** 2)
        for top, bottom, rho in layers
    )


def continental_case(name: str, upper_rho: float, lower_rho: float, mantle_rho: float,
                     tibet_crust_km: float) -> dict[str, float | str]:
    base_depth = 125_000.0
    base_z = -base_depth
    india_surface = 500.0
    india_crust = 35_000.0
    india_upper = 20_000.0
    tibet_crust = tibet_crust_km * 1000
    tibet_upper = min(india_upper, tibet_crust)
    india_moho = india_surface - india_crust
    india_upper_base = india_surface - india_upper
    india_layers = [
        (india_surface, india_upper_base, upper_rho),
        (india_upper_base, india_moho, lower_rho),
        (india_moho, -base_depth, mantle_rho),
    ]
    # Airy mass balance to the same compensation depth determines plateau elevation.
    india_mass = upper_rho * india_upper + lower_rho * (india_crust - india_upper)
    tibet_mass = upper_rho * tibet_upper + lower_rho * (tibet_crust - tibet_upper)
    tibet_surface = india_surface + (
        india_mass - tibet_mass + mantle_rho * (tibet_crust - india_crust)
    ) / mantle_rho
    tibet_moho = tibet_surface - tibet_crust
    tibet_upper_base = tibet_surface - tibet_upper
    tibet_layers = [
        (tibet_surface, tibet_upper_base, upper_rho),
        (tibet_upper_base, tibet_moho, lower_rho),
        (tibet_moho, -base_depth, mantle_rho),
    ]
    india_gpe = column_gpe(india_layers, base_z)
    tibet_gpe = column_gpe(tibet_layers, base_z)
    india_column_mass = sum((top - bottom) * rho for top, bottom, rho in india_layers)
    tibet_column_mass = sum((top - bottom) * rho for top, bottom, rho in tibet_layers)
    return {
        "case": name,
        "upper_crust_density_kg_m3": upper_rho,
        "lower_crust_density_kg_m3": lower_rho,
        "mantle_density_kg_m3": mantle_rho,
        "india_surface_elevation_km": india_surface / 1000,
        "india_crust_thickness_km": india_crust / 1000,
        "tibet_surface_elevation_km_airy_balanced": tibet_surface / 1000,
        "tibet_crust_thickness_km": tibet_crust / 1000,
        "compensation_depth_km_below_sea_level": base_depth / 1000,
        "column_mass_mismatch_fraction": abs(tibet_column_mass - india_column_mass) / india_column_mass,
        "india_gpe_N_per_m": india_gpe,
        "tibet_gpe_N_per_m": tibet_gpe,
        "delta_gpe_N_per_m": tibet_gpe - india_gpe,
    }


def ridge_case(name: str, thermal_density_deficit: float, anomaly_depth_km: float,
               ridge_flank_depth_difference_km: float) -> dict[str, float | str]:
    depth = anomaly_depth_km * 1000
    relief = ridge_flank_depth_difference_km * 1000
    water_density = 1030.0
    asthenosphere_density = 3330.0
    mantle_water_density_contrast = asthenosphere_density - water_density
    thermal_gpe = 0.5 * G_M_S2 * thermal_density_deficit * depth**2
    bathymetry_gpe = 0.5 * G_M_S2 * mantle_water_density_contrast * relief**2
    total = thermal_gpe + bathymetry_gpe
    segment_length = 1000_000.0
    return {
        "case": name,
        "thermal_density_deficit_kg_m3": thermal_density_deficit,
        "thermal_anomaly_depth_km": anomaly_depth_km,
        "ridge_to_old_flank_depth_difference_km": ridge_flank_depth_difference_km,
        "thermal_column_gpe_N_per_m": thermal_gpe,
        "bathymetry_component_gpe_N_per_m": bathymetry_gpe,
        "delta_gpe_N_per_m": total,
        "stated_ridge_segment_length_km": segment_length / 1000,
        "net_force_on_segment_N": total * segment_length,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write compact JSON results")
    args = parser.parse_args()
    collision = [
        continental_case("low-density, 60 km Tibet crust", 2700, 2850, 3250, 60),
        continental_case("central-density, 65 km Tibet crust", 2750, 2900, 3300, 65),
        continental_case("high-density, 70 km Tibet crust", 2800, 3000, 3350, 70),
    ]
    ridge = [
        ridge_case("low thermal anomaly", 35, 90, 3.32),
        ridge_case("central thermal anomaly", 50, 105, 3.32),
        ridge_case("high thermal anomaly", 65, 120, 3.32),
    ]
    result = {
        "method": {
            "g_m_s2": G_M_S2,
            "column_integral": "GPE = integral rho*g*(z - common basal datum) dz, z positive upward",
            "units": "J/m^2 = N/m",
            "ridge_segment_length_km": 1000,
            "ridge_note": "uniform thermal density-deficit layer plus the RFC relief term; first-order bound, not a full plate-cooling solution",
            "collision_note": "20-km upper crust, remaining lower crust, mantle to 125-km common compensation depth; Tibet elevation solved for Airy equal-column-mass balance",
        },
        "ridge_force_per_length_benchmark_N_per_m": 2.5e12,
        "himalaya_tibet_force_per_length_benchmark_N_per_m": [5e12, 6e12],
        "ridge_cases": ridge,
        "collision_cases": collision,
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()

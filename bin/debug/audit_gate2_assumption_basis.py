#!/usr/bin/env python3
"""Compute literature-comparison quantities for the Gate 2 assumption audit.

This is a bookkeeping diagnostic, not a fit to legacy rift-cut outcomes. Run from the
repository root with:

    python3 bin/debug/audit_gate2_assumption_basis.py \
      --output analysis/issue487/gate2-assumption-basis.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def surface_heat_flow_mw_m2(
    crust_km: float, mantle_lid_km: float, base_temperature_c: float
) -> float:
    """Mirror the Gate 2 steady conductive profile's surface-flux equation."""
    surface_c = 10.0
    conductivity_w_m_k = 2.5
    crust_heat_w_m3 = 0.7e-6
    crust_m = crust_km * 1000.0
    mantle_m = mantle_lid_km * 1000.0
    total_m = crust_m + mantle_m
    flux_w_m2 = (
        conductivity_w_m_k * (base_temperature_c - surface_c)
        + crust_heat_w_m3 * crust_m * (0.5 * crust_m + mantle_m)
    ) / total_m
    return flux_w_m2 * 1000.0


def build_report() -> dict[str, object]:
    years_s = 365.25 * 24.0 * 60.0 * 60.0
    strain_rates = (1e-18, 1e-16, 1e-15, 1e-14)
    columns = []
    for total_km in (70.0, 90.0, 100.0, 135.0, 150.0):
        mantle_km = total_km - 35.0
        columns.append(
            {
                "total_lithosphere_km": total_km,
                "crust_km": 35.0,
                "mantle_lid_km": mantle_km,
                "surface_heat_flow_mw_m2_by_base_temperature_c": {
                    str(base_c): round(
                        surface_heat_flow_mw_m2(35.0, mantle_km, base_c), 3
                    )
                    for base_c in (1250.0, 1330.0, 1400.0)
                },
            }
        )
    return {
        "purpose": "Independent input-scale bookkeeping; no legacy event labels or loads used.",
        "thermal_profile": {
            "equation": "q_surface=[k*(T_base-T_surface)+A*Hc*(Hc/2+Hm)]/(Hc+Hm)",
            "assumptions": {
                "surface_temperature_c": 10.0,
                "conductivity_w_m_k": 2.5,
                "uniform_crust_heat_production_w_m3": 0.7e-6,
                "crust_km": 35.0,
            },
            "columns": columns,
            "limitation": "Steady 1-D conduction with fixed conductivity and no advective or magmatic heat transfer.",
        },
        "strain_rate_conversion": {
            "equation": "strain_rate_per_year = strain_rate_per_second * seconds_per_year",
            "values_per_second_to_per_year": {
                f"{rate:.0e}": rate * years_s for rate in strain_rates
            },
            "ears_geodetic_second_invariant_per_year": [0.0, 2e-8],
            "ears_geodetic_upper_equivalent_per_second": 2e-8 / years_s,
        },
        "pore_pressure_scale": {
            "hydrostatic_gradient_mpa_per_km_freshwater": 9.80665,
            "lithostatic_to_hydrostatic_gradient_ratio_for_density_kg_m3": {
                str(rho): 1000.0 / rho for rho in (2700.0, 2800.0, 2900.0, 3300.0)
            },
            "limitation": "Gradient ratios are illustrative density-only estimates; the true overburden includes changing crust and mantle densities, and connected fluids are not established throughout the mantle lithosphere.",
        },
        "k0_status": {
            "model_definition": "effective horizontal background compression / effective vertical overburden",
            "scenario_values": [0.6, 0.8, 1.0],
            "independent_constraint": "Regional stress observations constrain regime and orientation better than a spatially uniform stress magnitude ratio; K0 is not observed by this world model.",
            "interpretation": "0.8 is retained as a nominal grid midpoint, not a measured or calibrated value. Keep the full 0.6-1.0 sensitivity bracket.",
        },
        "mantle_creep": {
            "laws": ["dry dislocation-creep olivine", "wet dislocation-creep olivine"],
            "wet_water_content_h_per_1e6_si": 1000,
            "unmodeled_controls": ["water distribution", "grain size and diffusion creep", "fabric and anisotropy", "Peierls creep", "partial melt", "damage and grain-size evolution"],
            "limitation": "The dry/wet pair is a constitutive sensitivity bracket, not a mapped hydration range. Laboratory power laws are extrapolated to model strain rates.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(build_report(), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

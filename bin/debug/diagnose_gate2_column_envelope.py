#!/usr/bin/env python3
"""Explore RFC 0010 Gate 2's depth-integrated extensional strength hypothesis.

This is a standalone diagnostic, not production rheology. It combines a Mohr-Coulomb
normal-faulting increment with Arrhenius power-law creep, integrates the weaker additional
stress through representative lithosphere columns, and reports sensitivity to missing
thermal, fluid, composition, strain-rate, friction and background-stress state.

The reference columns use Mantle Bloom's Hc/Hm, densities, surface/base temperatures,
conductivity and crustal heat-production proxy. The constitutive examples are model
parameterizations from Chen & Gerya (2016), drawing on Ranalli (1995), plus a weak mafic-
granulite lower-crust alternative using Wang et al. (2012). Mantle wet/dry flow
laws follow Hirth & Kohlstedt (2003). Their lab-derived laws are extrapolated many orders
of magnitude in strain rate and assume idealized compositions; these outputs are scenarios,
not calibrated predictions.

Usage:
    backend/.venv/bin/python bin/debug/diagnose_gate2_column_envelope.py \\
        --output analysis/issue487/gate2-column-envelope.json
"""

from __future__ import annotations

import argparse
import dataclasses
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import lithosphere, orogeny  # noqa: E402

GAS_CONSTANT_J_MOL_K = 8.314462618
WATER_COH_PPM_H_SI = 1_000.0


@dataclass(frozen=True)
class FlowLaw:
    """Power-law creep parameters; stress is MPa, pressure is Pa, temperature is K."""

    label: str
    a_mpa_n_s: float
    n: float
    activation_j_mol: float
    activation_volume_m3_mol: float = 0.0
    water_exponent: float = 0.0
    water_ppm_h_si: float = 0.0


@dataclass(frozen=True)
class Layer:
    top_m: float
    bottom_m: float
    law: FlowLaw


@dataclass(frozen=True)
class Scenario:
    base_temperature_c: float
    strain_rate_s: float
    pore_pressure_ratio: float
    k0: float
    cohesion_pa: float
    friction_angle_deg: float
    lower_crust_law: str
    mantle_law: str


FLOW_LAWS = {
    # Wet quartzite parameterization used in geodynamic strength-envelope models;
    # effective water content is implicit in the selected law.
    "wet_quartzite": FlowLaw("wet quartzite", 10.0**-3.49, 2.3, 154_000.0, 8e-6),
    "plagioclase_an75": FlowLaw("plagioclase An75", 10.0**-3.48, 3.2, 238_000.0, 8e-6),
    # Wang et al. (2012) point fit: A=10^-2 MPa^-3.2 s^-1, E=244 kJ/mol.
    "mafic_granulite": FlowLaw("mafic granulite", 10.0**-2.0, 3.2, 244_000.0),
    # Hirth & Kohlstedt (2003), using C_OH in ppm H/Si.
    "dry_olivine": FlowLaw("dry olivine", 1.1e5, 3.5, 530_000.0, 15e-6),
    "wet_olivine": FlowLaw(
        "wet olivine", 20.0, 3.5, 480_000.0, 11e-6, 1.2, WATER_COH_PPM_H_SI
    ),
}


def make_layers(crust_type: str, lower_crust_law: str, mantle_law: str) -> tuple[list[Layer], float, float]:
    """Return illustrative layers and reference Hc/Hm for a broad crust class."""
    if crust_type == "continental":
        hc_m, hm_m = lithosphere.reference_thickness("continental")
    else:
        hc_m = lithosphere.REFERENCE_HC_OCEANIC_M
        hm_m = (
            lithosphere.YOUNG_RIDGE_HM_M
            if crust_type == "oceanic_young_ridge"
            else lithosphere.REFERENCE_HM_OCEANIC_M
        )
    if crust_type == "continental":
        upper_crust_m = min(15_000.0, hc_m)
        layers = [
            Layer(0.0, upper_crust_m, FLOW_LAWS["wet_quartzite"]),
            Layer(upper_crust_m, hc_m, FLOW_LAWS[lower_crust_law]),
        ]
    else:
        # The model's 11 km oceanic Hc folds basalt and sediment effects together.
        # An75 is only an effective mafic-crust rheology proxy for this column test.
        layers = [Layer(0.0, hc_m, FLOW_LAWS["plagioclase_an75"])]
    layers.append(Layer(hc_m, hc_m + hm_m, FLOW_LAWS[mantle_law]))
    return layers, hc_m, hm_m


def depth_grid(total_m: float, hc_m: float, layers: list[Layer], dz_m: float) -> np.ndarray:
    regular = np.arange(0.0, total_m, dz_m)
    boundaries = [0.0, hc_m, total_m]
    boundaries.extend(layer.top_m for layer in layers)
    boundaries.extend(layer.bottom_m for layer in layers)
    return np.unique(np.concatenate((regular, np.asarray(boundaries, dtype=float))))


def profile_fields(
    depth_m: np.ndarray,
    hc_m: float,
    hm_m: float,
    crust_type: str,
    base_temperature_c: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return steady temperature (K) and lithostatic pressure (Pa)."""
    k = orogeny.CRUST_CONDUCTIVITY_W_M_K
    heat = orogeny.CRUST_HEAT_PRODUCTION_W_M3
    total_m = max(hc_m + hm_m, 1.0)
    q_surface = (
        k * (base_temperature_c - orogeny.SURFACE_TEMPERATURE_C)
        + heat * hc_m * (0.5 * hc_m + hm_m)
    ) / total_m
    moho_c = orogeny.SURFACE_TEMPERATURE_C + q_surface * hc_m / k - heat * hc_m**2 / (2.0 * k)
    crust_depth = np.minimum(depth_m, hc_m)
    temperature_c = (
        orogeny.SURFACE_TEMPERATURE_C
        + q_surface * crust_depth / k
        - heat * crust_depth**2 / (2.0 * k)
    )
    mantle_depth = np.maximum(depth_m - hc_m, 0.0)
    mantle_heat_flux = q_surface - heat * hc_m
    temperature_c = np.where(
        depth_m <= hc_m,
        temperature_c,
        moho_c + mantle_heat_flux * mantle_depth / k,
    )
    temperature_k = np.maximum(temperature_c + 273.15, 1.0)

    rho_crust = lithosphere.crust_density(crust_type)
    rho_mantle = lithosphere.RHO_LITHOSPHERE_MANTLE
    pressure_pa = np.where(
        depth_m <= hc_m,
        rho_crust * lithosphere.GRAVITY_M_S2 * depth_m,
        lithosphere.GRAVITY_M_S2
        * (rho_crust * hc_m + rho_mantle * (depth_m - hc_m)),
    )
    return temperature_k, pressure_pa


def creep_stress_pa(law: FlowLaw, strain_rate_s: float, pressure_pa: np.ndarray, temperature_k: np.ndarray) -> np.ndarray:
    """Invert epsilon_dot=A*sigma^n*C_OH^r*exp(-(E+PV)/RT)."""
    log_a = math.log(law.a_mpa_n_s)
    log_water = law.water_exponent * math.log(max(law.water_ppm_h_si, 1.0))
    activation = law.activation_j_mol + pressure_pa * law.activation_volume_m3_mol
    log_stress_mpa = (
        math.log(strain_rate_s) - log_a - log_water
    ) / law.n + activation / (law.n * GAS_CONSTANT_J_MOL_K * temperature_k)
    return np.exp(np.clip(log_stress_mpa, -100.0, 100.0)) * 1e6


def evaluate_case(crust_type: str, scenario: Scenario, dz_m: float) -> dict[str, Any]:
    layers, hc_m, hm_m = make_layers(crust_type, scenario.lower_crust_law, scenario.mantle_law)
    total_m = hc_m + hm_m
    depth_m = depth_grid(total_m, hc_m, layers, dz_m)
    temperature_k, pressure_pa = profile_fields(
        depth_m, hc_m, hm_m, crust_type, scenario.base_temperature_c
    )
    pore_pressure_pa = scenario.pore_pressure_ratio * pressure_pa
    vertical_effective_pa = pressure_pa - pore_pressure_pa
    friction = math.radians(scenario.friction_angle_deg)
    sin_phi, cos_phi = math.sin(friction), math.cos(friction)
    horizontal_background_pa = scenario.k0 * vertical_effective_pa
    horizontal_failure_pa = (
        vertical_effective_pa * (1.0 - sin_phi) - 2.0 * scenario.cohesion_pa * cos_phi
    ) / (1.0 + sin_phi)
    brittle_increment_pa = np.maximum(horizontal_background_pa - horizontal_failure_pa, 0.0)

    ductile_stress_pa = np.empty_like(depth_m)
    law_labels = np.empty(depth_m.shape, dtype=object)
    for layer in layers:
        mask = (depth_m >= layer.top_m) & (depth_m <= layer.bottom_m)
        ductile_stress_pa[mask] = creep_stress_pa(
            layer.law, scenario.strain_rate_s, pressure_pa[mask], temperature_k[mask]
        )
        law_labels[mask] = layer.law.label
    background_differential_pa = np.abs(vertical_effective_pa - horizontal_background_pa)
    ductile_increment_pa = np.maximum(ductile_stress_pa - background_differential_pa, 0.0)

    ductile_controls = ductile_increment_pa < brittle_increment_pa
    yield_increment_pa = np.minimum(brittle_increment_pa, ductile_increment_pa)
    capacity_n_per_m = float(np.trapezoid(yield_increment_pa, depth_m))
    brittle_integral_n_per_m = float(np.trapezoid(brittle_increment_pa, depth_m))
    ductile_integral_n_per_m = float(np.trapezoid(ductile_increment_pa, depth_m))
    control_depth_fraction = float(np.trapezoid(ductile_controls.astype(float), depth_m) / total_m)
    switches = np.flatnonzero(ductile_controls[1:] != ductile_controls[:-1])
    crossover_m = float(depth_m[switches[0] + 1]) if len(switches) else None

    moho_i = int(np.argmin(np.abs(depth_m - hc_m)))
    base_i = -1
    return {
        "crust_type": crust_type,
        "scenario": {
            "base_temperature_c": scenario.base_temperature_c,
            "strain_rate_s-1": scenario.strain_rate_s,
            "pore_pressure_ratio": scenario.pore_pressure_ratio,
            "k0_horizontal_over_vertical": scenario.k0,
            "cohesion_pa": scenario.cohesion_pa,
            "friction_angle_deg": scenario.friction_angle_deg,
            "lower_crust_law": scenario.lower_crust_law,
            "mantle_law": scenario.mantle_law,
        },
        "profile": {
            "crust_thickness_m": hc_m,
            "mantle_lithosphere_thickness_m": hm_m,
            "moho_temperature_c": float(temperature_k[moho_i] - 273.15),
            "lithosphere_base_temperature_c": float(temperature_k[base_i] - 273.15),
            "lithosphere_base_pressure_gpa": float(pressure_pa[base_i] / 1e9),
        },
        "integrated_additional_strength_n_per_m": capacity_n_per_m,
        "integrated_brittle_only_n_per_m": brittle_integral_n_per_m,
        "integrated_ductile_only_n_per_m": ductile_integral_n_per_m,
        "ductile_controlled_depth_fraction": control_depth_fraction,
        "first_brittle_ductile_crossover_depth_m": crossover_m,
        "moho": {
            "ductile_flow_stress_mpa": float(ductile_stress_pa[moho_i] / 1e6),
            "controlling_branch": "ductile" if ductile_controls[moho_i] else "brittle",
            "creep_law": str(law_labels[moho_i]),
        },
    }


def scenarios(include_lower_crust: bool = True) -> list[Scenario]:
    values = itertools.product(
        (1250.0, orogeny.LITHOSPHERE_BASE_TEMPERATURE_C, 1400.0),
        (1e-16, 1e-15, 1e-14),
        (0.0, 0.42),
        # The brittle branch assumes normal faulting, with vertical compression
        # as the major principal stress. K0 > 1 needs a different failure mode.
        (0.6, 0.8, 1.0),
        (1e6, 20e6),
        (math.degrees(math.atan(0.6)), 30.0),
        (("plagioclase_an75", "mafic_granulite") if include_lower_crust else ("plagioclase_an75",)),
        ("dry_olivine", "wet_olivine"),
    )
    return [Scenario(*item) for item in values]


def percentile_summary(cases: list[dict[str, Any]]) -> dict[str, Any]:
    capacities = np.asarray([case["integrated_additional_strength_n_per_m"] for case in cases])
    ordered = sorted(cases, key=lambda case: case["integrated_additional_strength_n_per_m"])

    def extreme(case: dict[str, Any]) -> dict[str, Any]:
        return {
            "capacity_tn_per_m": case["integrated_additional_strength_n_per_m"] / 1e12,
            "scenario": case["scenario"],
            "controlling_branch": (
                "ductile"
                if case["ductile_controlled_depth_fraction"] >= 0.5
                else "brittle"
            ),
        }

    return {
        "case_count": len(cases),
        "capacity_tn_per_m_p10_median_p90": [
            float(value / 1e12) for value in np.percentile(capacities, (10, 50, 90))
        ],
        "ductile_controlled_depth_fraction_median": float(
            np.median([case["ductile_controlled_depth_fraction"] for case in cases])
        ),
        "weakest_case": extreme(ordered[0]),
        "strongest_case": extreme(ordered[-1]),
    }


def run(dz_km: float) -> dict[str, Any]:
    dz_m = dz_km * 1000.0
    reference_params = Scenario(
        base_temperature_c=orogeny.LITHOSPHERE_BASE_TEMPERATURE_C,
        strain_rate_s=1e-15,
        pore_pressure_ratio=0.42,
        k0=1.0,
        cohesion_pa=20e6,
        friction_angle_deg=30.0,
        lower_crust_law="plagioclase_an75",
        mantle_law="dry_olivine",
    )
    crust_types = ("continental", "oceanic_young_ridge", "oceanic_mature")
    reference = [evaluate_case(crust, reference_params, dz_m) for crust in crust_types]
    scenario_sets = {
        crust: scenarios(include_lower_crust=crust == "continental") for crust in crust_types
    }
    swept = {
        crust: [evaluate_case(crust, scenario, dz_m) for scenario in scenario_sets[crust]]
        for crust in crust_types
    }
    sweep_dimensions = {
        "base_temperature_c": (1250.0, 1330.0, 1400.0),
        "strain_rate_s-1": (1e-16, 1e-15, 1e-14),
        "pore_pressure_ratio": (0.0, 0.42),
        "k0_horizontal_over_vertical": (0.6, 0.8, 1.0),
        "cohesion_pa": (1e6, 20e6),
        "friction_angle_deg": (math.degrees(math.atan(0.6)), 30.0),
        "lower_crust_law": ("plagioclase_an75", "mafic_granulite"),
        "mantle_law": ("dry_olivine", "wet_olivine"),
    }

    def parameter_sensitivity(cases: list[dict[str, Any]], crust_type: str) -> dict[str, Any]:
        summaries: dict[str, Any] = {}
        dimensions = sweep_dimensions.items()
        if crust_type != "continental":
            dimensions = ((key, values) for key, values in dimensions if key != "lower_crust_law")
        for dimension, values in dimensions:
            summaries[dimension] = {
                str(value): percentile_summary(
                    [case for case in cases if case["scenario"][dimension] == value]
                )
                for value in values
            }
        return summaries

    return {
        "schema_version": 1,
        "diagnostic": "RFC 0010 Gate 2 column strength envelope; not production rheology",
        "integration_depth_step_km": dz_km,
        "material_scenarios": {
            name: dataclasses.asdict(law) for name, law in FLOW_LAWS.items()
        },
        "equations": {
            "brittle_failure_horizontal_stress_pa": "((1 - sin(phi))*sigma_v_eff - 2*C*cos(phi))/(1 + sin(phi))",
            "additional_brittle_stress_pa": "max(0, K0*sigma_v_eff - sigma_h_failure)",
            "power_law_creep": "strain_rate = A*sigma^n*COH^r*exp(-(E + P*V)/(R*T))",
            "additional_ductile_stress_pa": "max(0, sigma_creep - abs(sigma_v_eff - K0*sigma_v_eff))",
            "integrated_capacity_n_per_m": "integral(min(additional_brittle, additional_ductile), depth)",
        },
        "assumptions": [
            "Steady 1-D conduction; one conductivity and crustal heat production are used through the column proxy.",
            "One effective continental split uses 15 km quartz-rich upper crust and remaining feldspar or granulite lower crust.",
            "The entire effective oceanic crust uses one An75 plagioclase-based mafic flow law.",
            "Oceanic columns include the model's 8 km young-ridge and 40 km mature-reference lids; no ocean-age distribution is represented.",
            "Wet mantle uses 1000 ppm H/Si; the world has no mantle water field.",
            "Pore pressure is a fixed fraction of lithostatic pressure, not a solved fluid-pressure profile.",
            "The vertical stress is the maximum principal compression; K0 defines initial horizontal stress.",
            "Laboratory creep laws are extrapolated to 1e-16 to 1e-14 per second and ignore semibrittle, Peierls, diffusion-creep transitions, melt and damage feedback.",
            "The mafic-granulite central flow-law fit is from reconstituted nominally dry rock containing 0.16-0.28 wt% dissolved water; its reported parameter uncertainties are not sampled.",
        ],
        "parameter_sources": {
            "crust_and_dry_mantle_flow_laws": "Chen & Gerya (2016), doi:10.1002/2016JB012872; parameterizations based on Ranalli (1995)",
            "mafic_granulite_flow_law": "Wang et al. (2012), doi:10.1016/j.epsl.2012.08.004",
            "wet_dry_olivine_flow_laws": "Hirth & Kohlstedt (2003), doi:10.1029/138GM06",
            "brittle_friction": "Byerlee (1978), doi:10.1007/BF00876528",
            "thickness_density_and_thermal_baseline": "backend/app/lithosphere.py and backend/app/orogeny.py",
        },
        "reference_cases": reference,
        "sensitivity": {
            "sweep_values": sweep_dimensions,
            "summary_by_crust_type": {
                crust: {
                    "overall": percentile_summary(cases),
                    "marginal_by_parameter": parameter_sensitivity(cases, crust),
                }
                for crust, cases in swept.items()
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dz-km", type=float, default=0.125, help="Depth integration spacing (default: 0.125 km).")
    parser.add_argument("--output", type=Path, help="Write JSON here; stdout if omitted.")
    args = parser.parse_args()
    if not math.isfinite(args.dz_km) or args.dz_km <= 0.0:
        parser.error("--dz-km must be positive and finite")
    result = run(args.dz_km)
    rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

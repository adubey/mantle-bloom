"""Exercise RFC 0012's G2 load/motion contract for rift, transform and lock cases.

Run from the repository root:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_physics_first_g3.py

This is a compact synthetic verification probe, not a calibrated constitutive model. The
normal contact solver is reused by orienting each process response against its positive
failure/slip direction. Geometry, integrated strength and shear capacity remain inputs from
the process owning the segment.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from app.quasistatic import solve_quasistatic_contacts


def integrated_rift_strength_n_per_m(crust_thickness_m: float) -> tuple[float, dict[str, float]]:
    """Integrate a brittle/ductile yield envelope for one continental column.

    This is a literature-parameterized diagnostic end member, not the simulator's current
    constitutive law. The column has a fixed 100 km lithosphere, a piecewise-linear geotherm
    from 0°C at the surface to 650°C at the Moho and 1300°C at the LAB, Byerlee-like dry
    friction, wet-quartzite crustal creep, and dry-olivine mantle creep. The weaker of the
    brittle and ductile differential stresses is integrated over depth.
    """
    total_thickness_m = 100_000.0
    mantle_thickness_m = total_thickness_m - crust_thickness_m
    depth_m = np.linspace(0.0, total_thickness_m, 1001)
    crust_mask = depth_m <= crust_thickness_m
    temperature_k = np.where(
        crust_mask,
        273.15 + (923.15 - 273.15) * depth_m / crust_thickness_m,
        923.15 + (1573.15 - 923.15) * (depth_m - crust_thickness_m) / mantle_thickness_m,
    )
    density_kg_m3 = np.where(crust_mask, 2800.0, 3300.0)
    dz_m = depth_m[1] - depth_m[0]
    overburden_pa = np.cumsum(density_kg_m3 * 9.81 * dz_m)
    # Mohr-Coulomb friction approximation, effective pressure = full lithostatic pressure.
    brittle_strength_pa = 0.6 * overburden_pa

    # Power-law dislocation creep: strain_rate = A * sigma_MPa**n * exp(-Q / (R*T)).
    # Wet Simpson quartzite (Koch et al., 1989) and dry olivine (Hirth & Kohlstedt, 2003).
    strain_rate_per_s = 1e-15
    gas_constant_j_mol_k = 8.314462618
    quartz_a_mpa_n_s = 5.05e-6
    quartz_n = 2.61
    quartz_q_j_mol = 145e3
    olivine_a_mpa_n_s = 1.1e5
    olivine_n = 3.5
    olivine_q_j_mol = 530e3
    preexponential = np.where(crust_mask, quartz_a_mpa_n_s, olivine_a_mpa_n_s)
    stress_exponent = np.where(crust_mask, quartz_n, olivine_n)
    activation_energy_j_mol = np.where(crust_mask, quartz_q_j_mol, olivine_q_j_mol)
    ductile_strength_mpa = (
        strain_rate_per_s
        / preexponential
        * np.exp(activation_energy_j_mol / (gas_constant_j_mol_k * temperature_k))
    ) ** (1.0 / stress_exponent)
    yield_strength_pa = np.minimum(brittle_strength_pa, ductile_strength_mpa * 1e6)
    strength_n_per_m = float(np.trapezoid(yield_strength_pa, depth_m))
    return strength_n_per_m, {
        "crust_thickness_m": crust_thickness_m,
        "lithosphere_thickness_m": total_thickness_m,
        "moho_temperature_k": 923.15,
        "lithosphere_base_temperature_k": 1573.15,
        "reference_strain_rate_per_s": strain_rate_per_s,
        "peak_envelope_strength_pa": float(np.max(yield_strength_pa)),
        "integrated_strength_n_per_m": strength_n_per_m,
    }


def solve_resistive_segment(
    drive_n_per_m: float,
    capacity_n_per_m: float,
    *,
    segment_length_m: float = 100_000.0,
    drag_n_s_per_m: float = 1e27,
) -> dict[str, float | bool]:
    """Two free plates, with a shared segment at the planetary radius."""
    length = segment_length_m
    radius_m = 6_371_000.0
    seconds_per_year = 31_557_600.0
    # Independent plate drags are eta*R²/seconds-per-year in N·m per rad/year.
    # The equal-and-opposite line load integrates over segment length, then applies its
    # planetary-radius lever arm to produce torque. The rate map converts rad/year to m/s.
    drag = np.eye(2) * drag_n_s_per_m * radius_m**2 / seconds_per_year
    drive = np.array([-1.0, 1.0]) * drive_n_per_m * length * radius_m
    force_map = np.array([[1.0], [-1.0]]) * length * radius_m
    failure_rate = np.array([[-1.0, 1.0]]) * radius_m / seconds_per_year
    result = solve_quasistatic_contacts(
        drag,
        drive,
        force_map,
        failure_rate,
        np.array([capacity_n_per_m]),
        rate_tolerance=1e-13,
    )
    return {
        "drive_n_per_m": drive_n_per_m,
        "capacity_n_per_m": capacity_n_per_m,
        "resistance_n_per_m": float(result.contact_force_n_per_m[0]),
        "relative_rate_m_per_s": float(result.closing_rate[0]),
        "equal_opposite": bool(np.isclose(np.sum(result.contact_force_n_per_m[0] * force_map[:, 0]), 0.0)),
        "residual_torque_nm": float(np.max(np.abs(result.residual_torque_nm))),
        "residual_torque_fraction": float(
            np.max(np.abs(result.residual_torque_nm)) / max(np.linalg.norm(drive), 1.0)
        ),
        "converged": result.converged,
    }


def fusion_epsilon_sensitivity() -> dict[str, object]:
    """Count historical G1 fusion contacts that meet instantaneous |v| < epsilon.

    G1 records rates at fusion, not each preceding contact timestep. These counts bound
    instantaneous epsilon sensitivity only; they cannot establish a duration tau.
    """
    rows = []
    for path in sorted(Path("analysis/rfc0012").glob("rules_*.json")):
        payload = json.loads(path.read_text())
        for event in payload["merges"]["events"]:
            rows.append(
                {
                    "rate": abs(float(event["closing_abs_median_cm_yr"])),
                    "duration": float(event.get("threshold_myr", 0.0)),
                }
            )
    epsilons = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
    durations = [0.0, 10.0, 25.0, 50.0, 75.0, 100.0]
    return {
        "source": "G1 rules_*.json fusion-event contact rate summaries",
        "n_fusions": len(rows),
        "epsilon_tau_sensitivity": [
            {
                "epsilon_cm_per_yr": eps,
                "tau_myr": tau,
                "instantaneous_rate_and_threshold_duration_candidates": sum(
                    row["rate"] < eps and row["duration"] >= tau for row in rows
                ),
            }
            for eps in epsilons
            for tau in durations
        ],
        "limitation": (
            "fusion-time rate and randomized convergence threshold only; qualifying counts "
            "are upper bounds, since prior per-step rates are not recorded to prove continuous locking"
        ),
    }


def geological_scale_widening_sweep() -> dict[str, object]:
    """Measure yielded width across a prescribed 1,000 km continental cross-section.

    Each 25 km column receives the same extensional membrane force, as in a 1-D
    force-balance diagnostic. The profile is deliberately a frozen heterogeneity sweep:
    there is no lateral stress transfer, thermal evolution, or strain weakening.
    """
    x_m = np.arange(-487_500.0, 487_500.0 + 1.0, 25_000.0)
    crust_thickness_m = 15_000.0 + 15_000.0 * np.minimum(np.abs(x_m) / 250_000.0, 1.0)
    capacities = np.array(
        [integrated_rift_strength_n_per_m(float(thickness))[0] for thickness in crust_thickness_m]
    )
    forces_n_per_m = [4.3e12, 4.4e12, 4.5e12, 4.65e12, 4.8e12, 4.95e12]
    sweep = []
    for force in forces_n_per_m:
        local = [solve_resistive_segment(force, float(capacity)) for capacity in capacities]
        active = np.array([item["relative_rate_m_per_s"] > 1e-14 for item in local])
        active_x = x_m[active]
        sweep.append(
            {
                "uniform_extensional_force_n_per_m": force,
                "active_column_count": int(np.count_nonzero(active)),
                "active_width_km": float(np.count_nonzero(active) * 25.0),
                "active_x_min_km": float(active_x[0] / 1000.0) if len(active_x) else None,
                "active_x_max_km": float(active_x[-1] / 1000.0) if len(active_x) else None,
                "mean_active_extension_rate_cm_per_yr": float(
                    np.mean([item["relative_rate_m_per_s"] for item, is_active in zip(local, active) if is_active])
                    * 100.0 * 31_557_600.0
                ) if np.any(active) else 0.0,
            }
        )
    return {
        "cross_section_width_km": 1000.0,
        "column_spacing_km": 25.0,
        "column_count": int(len(x_m)),
        "crust_thickness_profile": "15 km at x=0, linear to 30 km at |x|=250 km, 30 km beyond",
        "thermal_profile_per_column": "same prescribed 0°C surface, 650°C Moho, 1300°C LAB as the strength envelope",
        "loading_assumption": "uniform extensional membrane force per along-strike metre at each column",
        "force_sweep": sweep,
        "limitation": "frozen strength heterogeneity; no lateral shear transfer, strain weakening, rift migration, or thermal evolution",
    }


def main() -> None:
    # Each candidate cut gets an explicit column yield envelope. It remains a diagnostic
    # assumption because present-day World fields do not provide a temperature profile.
    integrated_strength_n_per_m, nominal_strength_profile = integrated_rift_strength_n_per_m(20_000.0)
    cut_length_m = 100_000.0
    loaded_width_m = 100_000.0
    basal_traction_pa = 0.1e6
    basal_line_load_n_per_m = basal_traction_pa * loaded_width_m * cut_length_m / cut_length_m
    rift = []
    for total_factor in [0.5, 0.9, 1.0, 1.1, 1.5, 2.0]:
        total_drive = total_factor * integrated_strength_n_per_m
        boundary_tension = total_drive - basal_line_load_n_per_m
        result = solve_resistive_segment(total_drive, integrated_strength_n_per_m)
        result.update(
            boundary_tensile_load_n_per_m=boundary_tension,
            basal_traction_pa=basal_traction_pa,
            basal_loaded_width_m=loaded_width_m,
            basal_line_load_n_per_m=basal_line_load_n_per_m,
        )
        rift.append(result)
    transform_capacity_n_per_m = 0.5e12
    transform = [
        solve_resistive_segment(force, transform_capacity_n_per_m)
        for force in [0.25e12, 0.5e12, 0.75e12, 1.0e12]
    ]
    weak_capacity_n_per_m, weak_strength_profile = integrated_rift_strength_n_per_m(15_000.0)
    strong_capacity_n_per_m, strong_strength_profile = integrated_rift_strength_n_per_m(30_000.0)
    localization = [
        {
            "drive_n_per_m": drive,
            "weak_segment": solve_resistive_segment(drive, weak_capacity_n_per_m),
            "strong_segment": solve_resistive_segment(drive, strong_capacity_n_per_m),
        }
        for drive in np.linspace(0.5 * weak_capacity_n_per_m, 1.5 * strong_capacity_n_per_m, 7)
    ]
    output = {
        "kind": "synthetic_gate_g3",
        "units": {"load": "N/m", "rate": "m/s", "integrated_strength": "Pa*m = N/m"},
        "rift": {
            "integrated_strength_n_per_m": integrated_strength_n_per_m,
            "strength_assumption": nominal_strength_profile,
            "strength_law": {
                "brittle": "0.6 * lithostatic pressure; full lithostatic effective pressure",
                "crust_ductile": "wet Simpson quartzite power-law creep, Koch et al. 1989",
                "mantle_ductile": "dry olivine dislocation creep, Hirth & Kohlstedt 2003",
                "envelope": "min(brittle, ductile differential stress), depth-integrated",
                "thermal_profile": "0°C surface, 650°C Moho, 1300°C at 100 km LAB; piecewise linear",
                "caveat": "unvalidated rheology extrapolated to 1e-15 s^-1; no water pressure, pressure-dependent creep, damage, melt, or simulator temperature field",
            },
            "force_sweep": rift,
            "localization_profiles": [weak_strength_profile, strong_strength_profile],
            "weaker_adjacent_segment_sweep": localization,
            "interpretation": "the 15 km crust profile yields before the 30 km profile; both extend after the stronger profile's integrated strength is exceeded",
            "geological_scale_widening": geological_scale_widening_sweep(),
        },
        "transform": {
            "shear_capacity_n_per_m": transform_capacity_n_per_m,
            "force_sweep": transform,
            "interpretation": "bounded resistance suppresses slip until capacity, then excess drive produces slip",
        },
        "merge_lock": fusion_epsilon_sensitivity(),
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

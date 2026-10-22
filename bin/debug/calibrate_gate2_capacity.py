#!/usr/bin/env python3
"""Audit Gate 2 capacity scales against published lithosphere benchmarks.

This diagnostic is a conditional calibration audit, not a fitted Earth strength model.
It keeps the existing 35 km continental crust, sweeps total lithosphere thickness and
all constitutive scenarios from diagnose_gate2_column_envelope.py, and compares a
named reference state with published rift/craton scale constraints and matched-replay
loads. The published values have different stress/load definitions and are treated as
order-of-magnitude checks, not labels for fitting.

Usage:
    backend/.venv/bin/python bin/debug/calibrate_gate2_capacity.py \\
        --output analysis/issue487/gate2-capacity-calibration.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import orogeny  # noqa: E402
from diagnose_gate2_column_envelope import (  # noqa: E402
    Scenario,
    evaluate_case,
    scenarios,
)

THICKNESS_CASES_KM = (70.0, 90.0, 100.0, 135.0, 150.0)
CAPACITY_SCALE_N_PER_M = 1e12


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _capacity(case: dict[str, Any]) -> float:
    return float(case["integrated_additional_strength_n_per_m"])


def _run_column(total_thickness_km: float, dz_km: float) -> dict[str, Any]:
    crust_m = 35_000.0
    total_m = total_thickness_km * 1000.0
    mantle_lid_m = total_m - crust_m
    if mantle_lid_m <= 0.0:
        raise ValueError("total lithosphere thickness must exceed 35 km crust")

    # Explicit reference state: hydrostatic pore-pressure proxy, modest extensional
    # prestress (K0=0.8), 1330 C base, 1e-15/s, dry olivine. These are design
    # assumptions, not values uniquely inferred from literature or world state.
    anchor = Scenario(
        base_temperature_c=orogeny.LITHOSPHERE_BASE_TEMPERATURE_C,
        strain_rate_s=1e-15,
        pore_pressure_ratio=0.42,
        k0=0.8,
        cohesion_pa=20e6,
        friction_angle_deg=30.0,
        lower_crust_law="plagioclase_an75",
        mantle_law="dry_olivine",
    )
    anchor_case = evaluate_case(
        "continental", anchor, dz_km * 1000.0, continental_mantle_lithosphere_m=mantle_lid_m
    )
    coarse_anchor_case = evaluate_case(
        "continental", anchor, dz_km * 2000.0, continental_mantle_lithosphere_m=mantle_lid_m
    )
    cases = [
        evaluate_case(
            "continental", scenario, dz_km * 1000.0,
            continental_mantle_lithosphere_m=mantle_lid_m,
        )
        for scenario in scenarios(include_lower_crust=True)
    ]
    values_tn_m = np.asarray([_capacity(case) / CAPACITY_SCALE_N_PER_M for case in cases])
    q10, q25, q50, q75, q90 = np.percentile(values_tn_m, (10, 25, 50, 75, 90))
    return {
        "total_lithosphere_thickness_km": total_thickness_km,
        "crust_thickness_km": crust_m / 1000.0,
        "mantle_lithosphere_thickness_km": mantle_lid_m / 1000.0,
        "anchor_capacity_tn_per_m": _capacity(anchor_case) / CAPACITY_SCALE_N_PER_M,
        "anchor_capacity_coarse_dz_tn_per_m": _capacity(coarse_anchor_case) / CAPACITY_SCALE_N_PER_M,
        "anchor_resolution_relative_difference": abs(
            _capacity(coarse_anchor_case) - _capacity(anchor_case)
        ) / max(_capacity(anchor_case), 1.0),
        "factorial_scenario_grid": {
            "case_count": len(cases),
            "capacity_tn_per_m_quantiles": {
                "p10": float(q10),
                "p25": float(q25),
                "median": float(q50),
                "p75": float(q75),
                "p90": float(q90),
            },
            "fraction_of_grid_in_east_african_rift_3_to_8_tn_per_m_scale": float(
                np.mean((values_tn_m >= 3.0) & (values_tn_m <= 8.0))
            ),
        },
        "anchor_state": {
            "base_temperature_c": anchor.base_temperature_c,
            "strain_rate_s-1": anchor.strain_rate_s,
            "pore_pressure_ratio": anchor.pore_pressure_ratio,
            "k0_horizontal_over_vertical": anchor.k0,
            "cohesion_pa": anchor.cohesion_pa,
            "friction_angle_deg": anchor.friction_angle_deg,
            "lower_crust_law": anchor.lower_crust_law,
            "mantle_law": anchor.mantle_law,
            "status": "explicit design anchor; not inferred or uniquely calibrated",
        },
    }


def run(dz_km: float, replay_path: Path) -> dict[str, Any]:
    replay = json.loads(replay_path.read_text(encoding="utf-8"))
    load = replay["stress_measure"]
    return {
        "schema_version": 1,
        "title": "RFC 0010 Gate 2 strength/capacity calibration audit",
        "diagnostic_status": "scale-constrained scenario bracket; no unique capacity selected",
        "units": {
            "capacity": "TN/m (10^12 N/m)",
            "shell_load": "TN/m (10^12 N/m)",
        },
        "method": {
            "script": "bin/debug/calibrate_gate2_capacity.py",
            "script_sha256": _sha256(Path(__file__).resolve()),
            "column_model": "bin/debug/diagnose_gate2_column_envelope.py",
            "column_model_sha256": _sha256(Path(__file__).with_name("diagnose_gate2_column_envelope.py")),
            "matched_replay_json": str(replay_path),
            "matched_replay_json_sha256": _sha256(replay_path),
            "integration_depth_step_km": dz_km,
            "sweep": "Existing 1,728-case continental factorial scenario grid, independently repeated at each thickness.",
            "thickness_mapping": "35 km crust held fixed; remaining total thickness assigned to mantle lithosphere.",
            "quantile_interpretation": "Designed factorial grid summaries, not probabilities or posterior confidence intervals.",
        },
        "columns": [_run_column(thickness, dz_km) for thickness in THICKNESS_CASES_KM],
        "matched_replay_load_context": {
            "population": "74 source-connected craton-crossing candidates, seeds 7, 11 and 2024; 22,896 cut-band node samples.",
            "band_median_of_event_band_medians_tn_per_m": load[
                "band_median_of_event_band_medians_tn_per_m"
            ],
            "median_of_event_band_p95_tn_per_m": load["median_of_event_band_p95_tn_per_m"],
            "interpretation": "Candidate loads from the old algorithm are not observed rupture/no-rupture labels and are not used to fit capacity.",
        },
        "external_scale_checks": [
            {
                "source": "Stamps et al. (2014), doi:10.1002/2013JB010717",
                "reported_quantity": "East African Rift integrated lithosphere strength of about 30-80 MPa for a 100 km lithosphere, i.e. 3-8 TN/m; estimated tectonic stresses about 1-2.2 TN/m.",
                "use": "Regional scale check. Definitions, spatial averaging and stress state are not identical to this model's incremental normal-fault capacity.",
                "url": "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JB010717",
            },
            {
                "source": "Rajaonarison et al. (2021), doi:10.1029/2020GL090483",
                "reported_quantity": "Geodynamic model sets lithosphere at 70 km for the Eastern Rift, 90 km for the Western Rift, and 150 km for Tanzania/Congo cratons.",
                "use": "Thickness-context scenarios only; this paper's modeled thicknesses do not directly give strength values.",
                "url": "https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL090483",
            },
            {
                "source": "Duretz et al. (2025), doi:10.1029/2025GC012409",
                "reported_quantity": "Their isotropic rifting model has peak driving force about 22 TN/m; models with strong inherited anisotropy stay at or below 10 TN/m.",
                "use": "Load-scale context only; driving force is not lithosphere capacity and cannot be used as a capacity target.",
                "url": "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025GC012409",
            },
        ],
        "limitations": [
            "No rupture/no-rupture observations are available for the matched replay candidates, so fitting a capacity multiplier to their candidate loads would be circular.",
            "Replay nodes carry per-cell crustal thickness, mantle-lid thickness, crust type, craton share, belt history, Moho thermal lag, and restite; this audit fixes 35 km crust and does not evaluate that spatial field.",
            "No full depth-temperature, pore-pressure, water, layer-mineralogy, grain-size, local strain-rate, fault-fabric, or K0 field is available; this keeps the capacity choice underidentified.",
            "All thickness scenarios retain a 35 km crust; this is not representative of every rift or craton and keeps the result conditional on the current continental column proxy.",
            "The literature checks use different definitions, averaging, boundary conditions, and constitutive assumptions; agreement in TN/m is only an order-of-magnitude check.",
            "The matched replay population is conditionally selected: connected-source craton-crossing cuts only, not the 130 candidates or all tectonic states.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dz-km", type=float, default=0.125, help="Depth integration spacing (default: 0.125 km).")
    parser.add_argument(
        "--replay-json", type=Path, default=Path("analysis/issue487/gate2-matched-replay.json"),
        help="Matched replay summary JSON.",
    )
    parser.add_argument("--output", type=Path, help="Write JSON here; stdout if omitted.")
    args = parser.parse_args()
    if not math.isfinite(args.dz_km) or args.dz_km <= 0.0:
        parser.error("--dz-km must be positive and finite")
    rendered = json.dumps(run(args.dz_km, args.replay_json), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Re-score fixed Gate 2 loads with independently bounded stress/creep/fluid cases.

The screen keeps the accepted connected-source Gate 1 loads fixed, uses the upper
distributed East African Rift geodetic strain rate (2e-8/year), and crosses it with
the conditional K0 and dry/wet olivine sensitivity cases. A near-critical K0 limit is
included as a friction-law endmember, not an independent estimate. A crust-hydrostatic case
holds pore pressure at its Moho value below the crust to avoid assuming a hydrostatic
fluid column through the mantle lid. Neither pressure profile is a fluid-flow solve.

Run from the repository root:

    backend/.venv/bin/python bin/debug/replay_gate2_assumption_screen.py \
        --samples /tmp/g2-gate2-sensitivity-seed7.jsonl /tmp/g2-gate2-sensitivity-seed11.jsonl /tmp/g2-gate2-sensitivity-seed2024.jsonl \
        --replay-report /tmp/g2-gate2-sensitivity-seed7.json /tmp/g2-gate2-sensitivity-seed11.json /tmp/g2-gate2-sensitivity-seed2024.json \
        --identity-reference /tmp/issue487-g2-spatial-g1baseline-seed7.json /tmp/issue487-g2-spatial-g1baseline-seed11.json /tmp/issue487-g2-spatial-g1baseline-seed2024.json \
        --output analysis/issue487/gate2-assumption-screen.json
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from diagnose_gate2_column_envelope import Scenario  # noqa: E402
import sensitivity_gate2_capacity as sensitivity  # noqa: E402

SECONDS_PER_YEAR = 365.25 * 24.0 * 60.0 * 60.0
GEODETIC_UPPER_RATE_PER_YEAR = 2e-8
GEODETIC_UPPER_RATE_PER_SECOND = GEODETIC_UPPER_RATE_PER_YEAR / SECONDS_PER_YEAR
# Derive the limiting case from the same friction angle used by the reference column.
REFERENCE_FRICTION_ANGLE_DEG = sensitivity.REFERENCE.friction_angle_deg
REFERENCE_FRICTION_COEFFICIENT = math.tan(math.radians(REFERENCE_FRICTION_ANGLE_DEG))
SIN_REFERENCE_FRICTION_ANGLE = math.sin(math.radians(REFERENCE_FRICTION_ANGLE_DEG))
K0_CRITICAL_NORMAL_FAULT = (1.0 - SIN_REFERENCE_FRICTION_ANGLE) / (
    1.0 + SIN_REFERENCE_FRICTION_ANGLE
)
# This friction-law limit is a conditional sensitivity, not an independent stress estimate.
K0_VALUES = (K0_CRITICAL_NORMAL_FAULT, 0.6, 0.8, 1.0)
MANTLE_LAWS = ("dry_olivine", "wet_olivine")
PORE_MODES = ("hydrostatic_water", "hydrostatic_crust_moho_hold")


def _scenario(k0: float, mantle_law: str, pore_mode: str) -> Scenario:
    reference = sensitivity.REFERENCE
    return dataclasses.replace(
        reference,
        strain_rate_s=GEODETIC_UPPER_RATE_PER_SECOND,
        k0=k0,
        mantle_law=mantle_law,
        pore_pressure_mode=pore_mode,
        pore_pressure_ratio=0.0,
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _screen(
    sample_paths: list[Path], report_paths: list[Path], reference_paths: list[Path]
) -> dict[str, Any]:
    events, input_validation = sensitivity._load_capture_inputs(
        sample_paths, report_paths, reference_paths
    )
    data = sensitivity._flatten(events)
    ext = data["n_ext"].astype(float)
    cases: dict[str, dict[str, Any]] = {}
    capacity_by_case: dict[tuple[float, str, str], np.ndarray] = {}
    for k0 in K0_VALUES:
        for mantle_law in MANTLE_LAWS:
            for pore_mode in PORE_MODES:
                key = (k0, mantle_law, pore_mode)
                name = f"k0_{k0:.2f}_{mantle_law}_{pore_mode}"
                print(f"Scoring {name}", flush=True)
                capacity = sensitivity._capacity_vector(
                    data, _scenario(k0, mantle_law, pore_mode), use_lag=False
                )
                capacity_by_case[key] = capacity
                cases[name] = {
                    "assumptions": {
                        "strain_rate_s-1": GEODETIC_UPPER_RATE_PER_SECOND,
                        "strain_rate_per_year": GEODETIC_UPPER_RATE_PER_YEAR,
                        "k0_horizontal_over_vertical": k0,
                        "mantle_law": mantle_law,
                        "pore_pressure_mode": pore_mode,
                        "temperature_profile": "existing steady conductive profile, 1330 C base",
                        "cohesion_pa": 20e6,
                        "friction_coefficient": REFERENCE_FRICTION_COEFFICIENT,
                    },
                    "summary": sensitivity._summaries(data, ext, capacity),
                }

    paired_pressure_comparison = {}
    for k0 in K0_VALUES:
        for mantle_law in MANTLE_LAWS:
            full = capacity_by_case[(k0, mantle_law, "hydrostatic_water")]
            capped = capacity_by_case[(k0, mantle_law, "hydrostatic_crust_moho_hold")]
            full_over = ext >= full
            capped_over = ext >= capped
            ratio_full = ext / np.maximum(full, 1.0)
            ratio_capped = ext / np.maximum(capped, 1.0)
            paired_pressure_comparison[f"k0_{k0:.2f}_{mantle_law}"] = {
                "comparison": "crust hydrostatic then hold at Moho pressure versus hydrostatic through full lithosphere",
                "median_capacity_change_percent": float(
                    100.0 * np.median((capped - full) / np.maximum(full, 1.0))
                ),
                "changed_threshold_nodes": int(np.count_nonzero(full_over != capped_over)),
                "fraction_classification_changed": float(np.mean(full_over != capped_over)),
                "newly_over_capacity_nodes": int(np.count_nonzero(capped_over & ~full_over)),
                "no_longer_over_capacity_nodes": int(np.count_nonzero(~capped_over & full_over)),
                "spearman_capacity_rank_correlation": float(spearmanr(full, capped).statistic),
                "spearman_n_ext_over_capacity_rank_correlation": float(
                    spearmanr(ratio_full, ratio_capped).statistic
                ),
                "full_hydrostatic_over_capacity_fraction": float(np.mean(full_over)),
                "crust_hydrostatic_moho_hold_over_capacity_fraction": float(np.mean(capped_over)),
            }

    capacity_matrix = np.vstack(list(capacity_by_case.values()))
    min_capacity = np.min(capacity_matrix, axis=0)
    max_capacity = np.max(capacity_matrix, axis=0)
    return {
        "schema_version": 1,
        "title": "Gate 2 independent assumption screen on fixed connected-source loads",
        "purpose": "Scenario sensitivity only. No event identity or legacy outcome is used to fit any parameter.",
        "input_validation": input_validation,
        "provenance": {
            "gate1_baseline_commit": sensitivity.GATE1_BASELINE_COMMIT,
            "screen_script_sha256": _sha256(Path(__file__).resolve()),
            "column_model_sha256": _sha256(ROOT / "bin/debug/diagnose_gate2_column_envelope.py"),
            "sensitivity_loader_sha256": _sha256(ROOT / "bin/debug/sensitivity_gate2_capacity.py"),
        },
        "method": {
            "node_loads_fixed": True,
            "sample_count": int(len(ext)),
            "event_count": int(sum(input_validation["expected_connected_events"].values())),
            "strain_rate_source_scale": {
                "geodetic_second_invariant_per_year": [0.0, GEODETIC_UPPER_RATE_PER_YEAR],
                "tested_upper_rate_per_year": GEODETIC_UPPER_RATE_PER_YEAR,
                "tested_upper_rate_per_second": GEODETIC_UPPER_RATE_PER_SECOND,
            },
            "crossed_dimensions": {
                "k0": list(K0_VALUES),
                "mantle_law": list(MANTLE_LAWS),
                "pore_pressure_mode": list(PORE_MODES),
            },
            "crust_hydrostatic_profile": "P_f=rho_water*g*min(depth,H_crust); pressure is continuous and held at its Moho value below the crust. This is a sensitivity cap, not a fluid-flow solution or evidence for connected deep fluids.",
            "limitations": [
                "The geodetic rate is a regional surface strain invariant, not a depth-resolved local creep rate.",
                "The K0 and dry/wet endmembers remain unobserved local state and are crossed as scenarios, not assigned probabilities; the near-critical K0 case is derived from the reference friction law, not external stress data.",
                "The two pore-pressure profiles bracket a deep hydrostatic gradient assumption; neither predicts fluid connectivity.",
                "The steady conductive temperature profile remains the existing 1-D model proxy.",
                "Capacity exceedance is an offline screen, not a rupture probability or a fitted yield label.",
            ],
        },
        "scenarios": cases,
        "paired_pressure_profile_comparison": paired_pressure_comparison,
        "crossed_scenario_envelope": {
            "interpretation": "Finite crossed scenario envelope, not a probability interval. A node is robustly over capacity if N_ext exceeds the largest tested capacity; it is possibly over capacity if it exceeds the smallest.",
            "robust_over_capacity": sensitivity._summary(
                np.ones(len(ext), dtype=bool), ext, max_capacity
            ),
            "possible_over_capacity": sensitivity._summary(
                np.ones(len(ext), dtype=bool), ext, min_capacity
            ),
            "nodes_changing_classification_across_scenarios": int(
                np.count_nonzero((ext >= min_capacity) != (ext >= max_capacity))
            ),
            "fraction_changing_classification_across_scenarios": float(
                np.mean((ext >= min_capacity) != (ext >= max_capacity))
            ),
        },
    }


def main() -> None:
    args = sensitivity.capture_cli_args(__doc__)
    result = _screen(args.samples, args.replay_report, args.identity_reference)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()

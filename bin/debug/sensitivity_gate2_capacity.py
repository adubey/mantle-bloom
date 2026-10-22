#!/usr/bin/env python3
"""Re-score fixed Gate 1 shell loads over a one-factor-at-a-time capacity bracket.

The input JSONL is captured with replay_gate2_capacity.py --capture-only and contains
the exact cut-band N_ext values plus the local column state. This script does not run
the world or fit capacities to split/failed-rift outcomes.

Usage:
    backend/.venv/bin/python bin/debug/sensitivity_gate2_capacity.py \
        --samples /tmp/g2-seed7.jsonl /tmp/g2-seed11.jsonl /tmp/g2-seed2024.jsonl \
        --replay-report /tmp/g2-seed7.json /tmp/g2-seed11.json /tmp/g2-seed2024.json \
        --identity-reference /tmp/g1-seed7.json /tmp/g1-seed11.json /tmp/g1-seed2024.json \
        --output analysis/issue487/gate2-capacity-sensitivity.json
"""

from __future__ import annotations

import argparse
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
sys.path.insert(0, str(ROOT / "backend"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from app import orogeny  # noqa: E402
from diagnose_gate2_column_envelope import Scenario, evaluate_case  # noqa: E402

GATE1_BASELINE_COMMIT = "0118752d0ffe75060464204d5b94f81a60fe9a07"
EXPECTED_EVENTS = {7: 44, 11: 42, 2024: 44}
EXPECTED_CONNECTED = {7: 23, 11: 26, 2024: 25}
EXPECTED_CUTS = {7: 144, 11: 111, 2024: 125}
EXPECTED_SAMPLES = 22_896
DZ_M = 125.0

REFERENCE = Scenario(
    base_temperature_c=orogeny.LITHOSPHERE_BASE_TEMPERATURE_C,
    strain_rate_s=1e-15,
    pore_pressure_ratio=0.42,
    k0=0.8,
    cohesion_pa=20e6,
    friction_angle_deg=30.0,
    lower_crust_law="plagioclase_an75",
    mantle_law="dry_olivine",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _variants() -> list[dict[str, Any]]:
    """Use documented factorial endpoints, changing one capacity assumption at a time."""
    return [
        {"name": "reference_steady", "scenario": REFERENCE, "thermal_mode": "steady"},
        {"name": "reference_stored_moho_lag", "scenario": REFERENCE, "thermal_mode": "stored_moho_lag"},
        {"name": "base_temperature_1250_c", "scenario": dataclasses.replace(REFERENCE, base_temperature_c=1250.0), "thermal_mode": "steady"},
        {"name": "base_temperature_1400_c", "scenario": dataclasses.replace(REFERENCE, base_temperature_c=1400.0), "thermal_mode": "steady"},
        {"name": "strain_rate_1e-16_s", "scenario": dataclasses.replace(REFERENCE, strain_rate_s=1e-16), "thermal_mode": "steady"},
        {"name": "strain_rate_1e-14_s", "scenario": dataclasses.replace(REFERENCE, strain_rate_s=1e-14), "thermal_mode": "steady"},
        {"name": "zero_pore_pressure_ratio", "scenario": dataclasses.replace(REFERENCE, pore_pressure_ratio=0.0), "thermal_mode": "steady"},
        {"name": "k0_0_6", "scenario": dataclasses.replace(REFERENCE, k0=0.6), "thermal_mode": "steady"},
        {"name": "k0_1_0", "scenario": dataclasses.replace(REFERENCE, k0=1.0), "thermal_mode": "steady"},
        {"name": "cohesion_1_mpa", "scenario": dataclasses.replace(REFERENCE, cohesion_pa=1e6), "thermal_mode": "steady"},
        {"name": "friction_atan_0_6", "scenario": dataclasses.replace(REFERENCE, friction_angle_deg=math.degrees(math.atan(0.6))), "thermal_mode": "steady"},
        {"name": "mafic_granulite_lower_crust", "scenario": dataclasses.replace(REFERENCE, lower_crust_law="mafic_granulite"), "thermal_mode": "steady"},
        {"name": "wet_olivine_mantle", "scenario": dataclasses.replace(REFERENCE, mantle_law="wet_olivine"), "thermal_mode": "steady"},
    ]


def _event_key(row: dict[str, Any]) -> tuple[float, int]:
    return (float(row["myr"]), int(row["plate"]))


def _validate_seed_report(seed: int, report: dict[str, Any], reference: dict[str, Any]) -> None:
    if not report.get("capture_only") or report.get("gate1_baseline_commit") != GATE1_BASELINE_COMMIT:
        raise ValueError(f"seed {seed} capture does not identify the pinned Gate 1 source")
    if report.get("events_scored") != EXPECTED_EVENTS[seed]:
        raise ValueError(f"seed {seed} candidate count changed: {report.get('events_scored')}")
    if report.get("production_cut_count") != EXPECTED_CUTS[seed] or not report.get("gate1_cut_count_matches"):
        raise ValueError(f"seed {seed} production cut count does not match Gate 1")
    if not report.get("gate1_craton_crossing_cut_count_matches"):
        raise ValueError(f"seed {seed} craton-crossing candidate count does not match Gate 1")
    if report["diagnostic_coverage"].get("valid_shell_solves") != EXPECTED_CONNECTED[seed]:
        raise ValueError(f"seed {seed} connected-source solve count changed")
    current_keys = [_event_key(row) for row in report["events"]]
    reference_keys = [_event_key(row) for row in reference["events"]]
    if current_keys != reference_keys:
        raise ValueError(f"seed {seed} ordered candidate time/plate identities differ from Gate 1 replay")


def _read_sample_rows(paths: list[Path]) -> list[dict[str, Any]]:
    rows = []
    for path in paths:
        rows.extend(json.loads(line) for line in path.read_text().splitlines() if line.strip())
    return rows


def _accepted_events(rows: list[dict[str, Any]]) -> dict[int, list[dict[str, Any]]]:
    by_seed: dict[int, list[dict[str, Any]]] = {seed: [] for seed in EXPECTED_CONNECTED}
    for row in rows:
        seed = int(row["seed"])
        if seed not in by_seed:
            raise ValueError(f"unexpected seed in node capture: {seed}")
        if (
            row["solver_stop_code"] == 1
            and row["equilibrium_residual_over_load"] <= 0.01
            and row["mesh_band_coverage"] >= 0.95
            and row["mesh_craton_band_coverage"] >= 0.95
        ):
            by_seed[seed].append(row)
    for seed, event_rows in by_seed.items():
        if len(event_rows) != EXPECTED_CONNECTED[seed]:
            raise ValueError(f"seed {seed} has {len(event_rows)} accepted captured events, expected {EXPECTED_CONNECTED[seed]}")
    return by_seed


def _validate_sample_arrays(events: list[dict[str, Any]]) -> int:
    for event in events:
        lengths = {
            len(event[key]) for key in (
                "source_cell_index", "n_ext_n_per_m", "h_crust_m", "h_mantle_m", "moho_lag_c",
                "crust_type", "craton_crust_fraction", "inherited_belt",
            )
        }
        if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
            raise ValueError("captured per-cell arrays are misaligned or empty")
    return sum(len(row["n_ext_n_per_m"]) for row in events)


def _load_capture_inputs(
    sample_paths: list[Path], report_paths: list[Path], reference_paths: list[Path]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not (len(sample_paths) == len(report_paths) == len(reference_paths) == 3):
        raise ValueError("provide exactly three sample files, replay reports and identity references")
    reports = {int(json.loads(path.read_text())["seed"]): json.loads(path.read_text()) for path in report_paths}
    references = {int(json.loads(path.read_text())["seed"]): json.loads(path.read_text()) for path in reference_paths}
    if set(reports) != set(EXPECTED_EVENTS) or set(references) != set(EXPECTED_EVENTS):
        raise ValueError("input reports must cover exactly seeds 7, 11 and 2024")

    for seed, report in reports.items():
        _validate_seed_report(seed, report, references[seed])
    by_seed = _accepted_events(_read_sample_rows(sample_paths))
    accepted_events = [row for seed in sorted(by_seed) for row in by_seed[seed]]
    sample_count = _validate_sample_arrays(accepted_events)
    if sample_count != EXPECTED_SAMPLES:
        raise ValueError(f"captured {sample_count} accepted band nodes, expected {EXPECTED_SAMPLES}")

    metadata = {
        "baseline_commit": GATE1_BASELINE_COMMIT,
        "expected_candidate_counts": EXPECTED_EVENTS,
        "expected_connected_events": EXPECTED_CONNECTED,
        "expected_production_cuts": EXPECTED_CUTS,
        "accepted_events": sum(map(len, by_seed.values())),
        "accepted_node_samples": sample_count,
        "ordered_candidate_time_plate_keys_match_gate1": True,
        "source_files": [
            {"path": str(path), "sha256": _sha256(path)}
            for path in (*sample_paths, *report_paths, *reference_paths)
        ],
    }
    return accepted_events, metadata


def _flatten(events: list[dict[str, Any]]) -> dict[str, np.ndarray]:
    data: dict[str, list[Any]] = {
        "seed": [], "event_index": [], "kind": [], "myr": [], "plate": [], "cell": [],
        "n_ext": [], "hc": [], "hm": [], "lag": [], "crust": [], "craton_fraction": [], "belt": [],
    }
    for event in events:
        count = len(event["n_ext_n_per_m"])
        for key, value in (
            ("seed", event["seed"]), ("event_index", event["event_index"]),
            ("kind", event["kind"]), ("myr", event["myr"]), ("plate", event["plate"]),
        ):
            data[key].extend([value] * count)
        for key, field in (
            ("cell", "source_cell_index"), ("n_ext", "n_ext_n_per_m"), ("hc", "h_crust_m"),
            ("hm", "h_mantle_m"), ("lag", "moho_lag_c"), ("crust", "crust_type"),
            ("craton_fraction", "craton_crust_fraction"), ("belt", "inherited_belt"),
        ):
            data[key].extend(event[field])
    return {
        key: np.asarray(value, dtype=object if key in ("kind", "crust") else None)
        for key, value in data.items()
    }


def _capacity_vector(data: dict[str, np.ndarray], scenario: Scenario, use_lag: bool) -> np.ndarray:
    states: list[tuple[str, int, int, int]] = []
    for index in range(len(data["n_ext"])):
        states.append((
            str(data["crust"][index]),
            int(round(float(data["hc"][index]) / 100.0) * 100),
            int(round(float(data["hm"][index]) / 100.0) * 100),
            int(round(float(data["lag"][index]))) if use_lag else 0,
        ))
    unique_states = list(dict.fromkeys(states))
    capacities_by_state: dict[tuple[str, int, int, int], float] = {}
    for state in unique_states:
        crust, hc, hm, lag = state
        result = evaluate_case(
            crust, scenario, DZ_M,
            column_crust_thickness_m=float(hc),
            column_mantle_lithosphere_thickness_m=float(hm),
            moho_thermal_lag_c=float(lag),
        )
        capacities_by_state[state] = float(result["integrated_additional_strength_n_per_m"])
    print(f"  solved {len(unique_states)} unique rounded column states", flush=True)
    return np.asarray([capacities_by_state[state] for state in states], dtype=float)


def _summary(mask: np.ndarray, ext: np.ndarray, cap: np.ndarray) -> dict[str, Any] | None:
    if not np.any(mask):
        return None
    n_ext, n_cap = ext[mask], cap[mask]
    ratio = n_ext / np.maximum(n_cap, 1.0)
    over = n_ext >= n_cap
    return {
        "node_count": int(mask.sum()),
        "capacity_n_per_m_p10_median_p90": [float(v) for v in np.percentile(n_cap, (10, 50, 90))],
        "n_ext_over_capacity_median_p95": [float(v) for v in np.percentile(ratio, (50, 95))],
        "over_capacity_nodes": int(over.sum()),
        "over_capacity_fraction": float(np.mean(over)),
    }


def _summaries(data: dict[str, np.ndarray], ext: np.ndarray, cap: np.ndarray) -> dict[str, Any]:
    masks: dict[str, dict[str, np.ndarray]] = {
        "by_seed": {str(seed): data["seed"] == seed for seed in EXPECTED_CONNECTED},
        "by_crust_type": {str(label): data["crust"] == label for label in sorted(set(data["crust"]))},
        "by_event_kind": {str(kind): data["kind"] == kind for kind in sorted(set(data["kind"]))},
        "by_inherited_belt": {
            "marked": data["belt"].astype(bool), "unmarked": ~data["belt"].astype(bool),
        },
        "by_craton_share": {
            "none": data["craton_fraction"].astype(float) <= 1e-9,
            "partial": (data["craton_fraction"].astype(float) > 1e-9) & (data["craton_fraction"].astype(float) < 0.5),
            "majority": data["craton_fraction"].astype(float) >= 0.5,
        },
    }
    return {
        "all": _summary(np.ones(len(ext), dtype=bool), ext, cap),
        **{
            group: {name: _summary(mask, ext, cap) for name, mask in members.items()}
            for group, members in masks.items()
        },
    }


def run(sample_paths: list[Path], report_paths: list[Path], reference_paths: list[Path]) -> dict[str, Any]:
    events, input_validation = _load_capture_inputs(sample_paths, report_paths, reference_paths)
    data = _flatten(events)
    ext = data["n_ext"].astype(float)
    variants = _variants()
    results: dict[str, Any] = {}
    capacity_vectors: list[np.ndarray] = []
    reference_over = None
    for variant in variants:
        print(f"Scoring {variant['name']} ({variant['thermal_mode']})", flush=True)
        cap = _capacity_vector(
            data, variant["scenario"], variant["thermal_mode"] == "stored_moho_lag"
        )
        capacity_vectors.append(cap)
        over = ext >= cap
        if variant["name"] == "reference_steady":
            reference_over = over
        item = {
            "assumption_changes": dataclasses.asdict(variant["scenario"]),
            "thermal_mode": variant["thermal_mode"],
            "one_factor_at_a_time": variant["name"] != "reference_steady" and variant["name"] != "reference_stored_moho_lag",
            "summaries": _summaries(data, ext, cap),
        }
        if reference_over is not None:
            item["classification_change_vs_reference_steady"] = {
                "changed_nodes": int(np.count_nonzero(over != reference_over)),
                "changed_fraction": float(np.mean(over != reference_over)),
                "newly_over_capacity_nodes": int(np.count_nonzero(over & ~reference_over)),
                "no_longer_over_capacity_nodes": int(np.count_nonzero(~over & reference_over)),
                "spearman_rank_correlation_of_n_ext_over_n_cap": float(spearmanr(ext / cap, ext / capacity_vectors[0]).statistic),
            }
        results[variant["name"]] = item

    capacity_stack = np.vstack(capacity_vectors)
    min_cap = np.min(capacity_stack, axis=0)
    max_cap = np.max(capacity_stack, axis=0)
    report = {
        "schema_version": 1,
        "title": "RFC 0010 Gate 2 capacity assumption sensitivity on fixed matched loads",
        "purpose": "One-factor-at-a-time sensitivity on the same Gate 1 node loads; not a probability distribution, posterior interval or fitted rupture model.",
        "input_validation": input_validation,
        "method": {
            "sample_count": len(ext),
            "column_depth_step_m": DZ_M,
            "capacity_cache_rounding": {"Hc_m": 100, "Hm_m": 100, "Moho_lag_c": 1},
            "reference_scenario": dataclasses.asdict(REFERENCE),
            "one_factor_at_a_time_variants": [variant["name"] for variant in variants],
            "ranges_are_the_existing_documented_factorial_endpoints": True,
            "thermal_assumption": "Stored Moho lag is compared with steady profile at the reference rheology; other OAT changes use the steady profile.",
            "limitations": [
                "OAT isolates first-order sensitivity and does not estimate interactions among assumptions.",
                "The stored Moho-lag interpolation is a diagnostic profile assumption, not a transient thermal solution.",
                "Young-ridge classification remains the active-divergence proxy; the simulation has no oceanic crust-age field.",
                "Threshold exceedance is a screening comparison and not a rupture probability or a calibrated yield label.",
            ],
        },
        "scenarios": results,
        "one_factor_scenario_envelope": {
            "interpretation": "Minimum/maximum capacity across the finite OAT variants; not a joint uncertainty envelope.",
            "possible_over_capacity_nodes": _summary(np.ones(len(ext), dtype=bool), ext, min_cap),
            "robust_over_capacity_nodes": _summary(np.ones(len(ext), dtype=bool), ext, max_cap),
        },
    }
    report["provenance"] = {
        "script_sha256": _sha256(Path(__file__).resolve()),
        "capture_replay_script_sha256": _sha256(Path(__file__).with_name("replay_gate2_capacity.py")),
        "column_model_sha256": _sha256(Path(__file__).with_name("diagnose_gate2_column_envelope.py")),
        "shell_solver_sha256": _sha256(ROOT / "backend/app/thin_shell.py"),
        "force_field_helpers_source_sha256": _sha256(ROOT / "backend/app/torque.py"),
        "capture_execution": (
            "Pinned Gate 1 source with the current diagnostic force-field helpers and thin-shell solver overlaid; "
            "baseline production torque routines retained."
        ),
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs=3, type=Path, required=True)
    parser.add_argument("--replay-report", nargs=3, type=Path, required=True)
    parser.add_argument("--identity-reference", nargs=3, type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rendered = json.dumps(run(args.samples, args.replay_report, args.identity_reference), indent=2, allow_nan=False) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()

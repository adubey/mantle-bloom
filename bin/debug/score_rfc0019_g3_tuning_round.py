"""Aggregate compact outputs from one preregistered G3 coast-sink tuning round.

Usage: backend/.venv/bin/python bin/debug/score_rfc0019_g3_tuning_round.py \
    CAPTURE_RESULT.json ... --out PRIVATE_DIR/round-score.json
This reads only registered tuning results; it never loads a save or capture array.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path

from rfc0019_g3_mixing import TUNING_SAVES

SAVE_FAMILY = {
    "seed910211954-221.6myr.bloom": "seed910211954",
    "seed910211954-225.1myr.bloom": "seed910211954",
    "seed910211954-224myr.bloom": "seed910211954",
    "mantle-bloom-seed317739385-200800000y.bloom": "seed317739385",
}


def _candidate_summary(records):
    screens = [record["screen"] for record in records]
    share_changes = [abs(screen["land_rain_share_delta_vs_control"]) for screen in screens]
    pattern_coverage = {
        name: any(screen["geographic_support"][name]["evaluable"] for screen in screens)
        for name in ("equatorial_belt", "inland_coast", "ridge_shadow", "cold_coast")
    }
    gains_by_family = {}
    for record in records:
        gain = record["screen"]["ridge_median_contrast_gain_mm_yr"]
        if gain is not None:
            family = SAVE_FAMILY[record["save"]]
            gains_by_family.setdefault(family, []).append(gain)
    family_gain_medians = [statistics.median(gains) for gains in gains_by_family.values()]
    scores = [screen["minimax_score"] for screen in screens if screen["minimax_score"] is not None]
    return {
        "worst_case_minimax_score": max(scores) if scores else None,
        "pattern_coverage": pattern_coverage,
        "all_required_patterns_evaluated": all(pattern_coverage.values()),
        "all_tuning_safeguards_pass": all(
            all(screen["safeguards"].values()) for screen in screens
        ),
        "all_hard_gates_pass": all(screen["all_hard_gates_pass"] for screen in screens),
        "all_round_gates_pass": all(screen["eligible"] for screen in screens),
        "all_tuning_pattern_checks_pass": (
            all(pattern_coverage.values()) and all(screen["pattern_pass"] for screen in screens)
        ),
        "max_absolute_land_share_delta": max(share_changes),
        "family_balanced_median_ridge_contrast_gain_mm_yr": (
            statistics.median(family_gain_medians) if family_gain_medians else None
        ),
        "case_count": len(records),
        "independent_family_count": len(gains_by_family),
        "cases": [{"save": record["save"], "wind": record["wind"], **record["screen"]}
                  for record in records],
    }


def _validate_report_headers(reports, parser):
    rounds = {report.get("tuning_round") for report in reports}
    if len(rounds) != 1 or None in rounds:
        parser.error("all inputs must belong to one declared coast tuning round")
    round_number = next(iter(rounds))
    if round_number not in (1, 2, 3):
        parser.error("tuning_round must be 1, 2, or 3")
    if any(report.get("heldout") for report in reports):
        parser.error("held-out results cannot be used for tuning selection")
    if any(report.get("solvers") != ["finite_volume"] for report in reports):
        parser.error("tuning selection is preregistered for finite_volume only")
    if any(report.get("density") != 4 for report in reports):
        parser.error("the round score requires density-4 tuning captures")
    if len({report.get("code_commit") for report in reports}) != 1:
        parser.error("round results must use the same source commit")
    source_fingerprints = {
        (report.get("source_provenance") or {}).get("fingerprint") for report in reports
    }
    if None in source_fingerprints or len(source_fingerprints) != 1:
        parser.error("round results must use identical complete source fingerprints")
    if any(report.get("baseline_reference") != "same-run ridge_025 arm; within-run comparisons only"
           for report in reports):
        parser.error("all tuning reports must use their same-run ridge_025 baseline")
    return round_number


def _expected_arms(round_number, alphas, parser):
    if round_number == 1:
        if alphas != [1.0, 3.0]:
            parser.error("round 1 requires alpha 1 and 3")
        return {"control", "ridge_025", "control_coast_1", "coast_1",
                "control_uniform_1", "uniform_1", "control_coast_3", "coast_3",
                "control_uniform_3", "uniform_3"}
    if round_number == 2:
        if alphas not in ([1.0], [3.0]):
            parser.error("round 2 requires exactly one selected alpha, 1 or 3")
        suffix = f"{alphas[0]:g}"
        return {"control", "ridge_025", f"control_coast_{suffix}", f"coast_{suffix}",
                f"control_uniform_{suffix}", f"uniform_{suffix}"}
    if alphas != [2.0]:
        parser.error("round 3 requires alpha 2")
    return {"control", "ridge_025", "control_coast_2", "coast_2",
            "control_uniform_2", "uniform_2"}


def _collect_round_results(reports, round_number, parser):
    expected = {name: digest for name, digest in TUNING_SAVES.items()}
    seen = {}
    by_arm = {}
    settings = {}
    for report in reports:
        save = Path(report["save"]).name
        if expected.get(save) != report.get("save_digest"):
            parser.error(f"unregistered tuning save or digest: {save}")
        sink = report.get("parameters", {}).get("coast_sink")
        if not isinstance(sink, dict):
            parser.error(f"missing coast-sink settings for {save}")
        length = sink.get("length_km")
        if not isinstance(length, (int, float)) or not math.isfinite(length):
            parser.error(f"invalid coast length for {save}")
        alphas = sorted(sink.get("alphas", {}).values())
        if set(report["results"].get("finite_volume", {})) != {"diagnostic", "cfd"}:
            parser.error(f"both registered wind modes are required for {save}")
        expected_arms = _expected_arms(round_number, alphas, parser)
        key = f"{save}@L{length:g}"
        if key in seen:
            parser.error(f"duplicate tuning input and length: {key}")
        seen[key] = report["save_digest"]
        settings[key] = {"save": save, "length_km": length, "alphas": alphas}
        for wind, arms in report["results"]["finite_volume"].items():
            if set(arms) != expected_arms:
                parser.error(f"wrong arm set for round {round_number}: {save}/{wind}")
            for name, metrics in arms.items():
                if name == "control" or name.startswith("control_"):
                    continue
                if "tuning_screen" not in metrics:
                    parser.error(f"missing preregistered score for {save}/{wind}/{name}")
                if not {"all_hard_gates_pass", "eligible", "pattern_pass"}.issubset(
                    metrics["tuning_screen"],
                ):
                    parser.error(f"missing hard-gate or pattern-gate results for {save}/{wind}/{name}")
                arm_key = f"{name}@L{length:g}"
                by_arm.setdefault(arm_key, []).append({
                    "save": save,
                    "wind": wind,
                    "screen": metrics["tuning_screen"],
                })
    return seen, by_arm, settings


def _validate_round_settings(round_number, seen, settings, parser):
    if round_number not in (1, 2, 3):
        parser.error("tuning_round must be 1, 2, or 3")
    expected_names = set(TUNING_SAVES)
    if not {Path(key.split("@L")[0]).name for key in seen} == expected_names:
        parser.error("supply each registered tuning save exactly once")
    expected_lengths = {1: {750.0}, 2: {500.0, 1500.0}, 3: None}[round_number]
    actual_lengths = {float(value["length_km"]) for value in settings.values()}
    if expected_lengths is not None and actual_lengths != expected_lengths:
        parser.error(f"round {round_number} requires lengths {sorted(expected_lengths)}")
    if round_number == 3 and len(actual_lengths) != 1:
        parser.error("round 3 requires one shared coast length")
    expected_alphas = {1: {1.0, 3.0}, 2: None, 3: {2.0}}[round_number]
    actual_alpha_sets = {tuple(value["alphas"]) for value in settings.values()}
    if expected_alphas is not None and actual_alpha_sets != {tuple(sorted(expected_alphas))}:
        parser.error(f"round {round_number} has an invalid alpha set")
    if round_number == 2 and (len(actual_alpha_sets) != 1 or next(iter(actual_alpha_sets)) not in {(1.0,), (3.0,)}):
        parser.error("round 2 requires one shared alpha of 1 or 3")
    required = {f"{name}@L{length:g}" for name in expected_names for length in actual_lengths}
    if set(seen) != required:
        parser.error("supply each registered tuning save once at every required length")
    return actual_lengths


def _build_comparisons(summaries, round_number, actual_lengths):
    comparisons = {}
    for alpha in ("1", "2", "3"):
        for length in sorted(actual_lengths):
            suffix = f"@L{length:g}"
            coast = summaries.get(f"coast_{alpha}{suffix}")
            uniform = summaries.get(f"uniform_{alpha}{suffix}")
            frozen = summaries.get(f"ridge_025{suffix}")
            if coast is None or uniform is None:
                continue
            frozen_gates_pass = frozen is not None and frozen["all_round_gates_pass"]
            uniform_gates_pass = uniform["all_round_gates_pass"]
            coast_score = coast["worst_case_minimax_score"]
            uniform_score = uniform["worst_case_minimax_score"]
            frozen_score = None if frozen is None else frozen["worst_case_minimax_score"]
            coast_freeze_eligible = (
                coast["all_round_gates_pass"]
                and coast["all_required_patterns_evaluated"]
                and coast["all_tuning_pattern_checks_pass"]
                and frozen_gates_pass
                and uniform_gates_pass
                and coast_score is not None
                and frozen_score is not None
                and uniform_score is not None
                and coast_score < frozen_score
                and coast_score < uniform_score
            )
            comparison = {
                "coast_round_gates_pass": coast["all_round_gates_pass"],
                "coast_pattern_checks_pass": coast["all_tuning_pattern_checks_pass"],
                "frozen_baseline_round_gates_pass": frozen_gates_pass,
                "matched_uniform_round_gates_pass": uniform_gates_pass,
                "coast_beats_frozen_candidate": (
                    frozen_gates_pass
                    and coast_score is not None
                    and frozen_score is not None
                    and coast_score < frozen_score
                ),
                "coast_beats_matched_uniform": (
                    uniform_gates_pass and coast_score is not None and uniform_score is not None
                    and coast_score < uniform_score
                ),
                "coast_score": coast_score,
                "uniform_score": uniform_score,
                "frozen_score": frozen_score,
            }
            comparison["may_advance_to_round2"] = (
                round_number == 1
                and comparison["coast_round_gates_pass"]
                and not coast_freeze_eligible
                and comparison["coast_beats_frozen_candidate"]
                and comparison["coast_beats_matched_uniform"]
            )
            comparison["may_freeze_for_heldout_check"] = (
                coast_freeze_eligible
            )
            comparisons[f"alpha_{alpha}_L{length:g}"] = comparison
    return comparisons


def _rank_arms(summaries):
    ranked_eligible_arms = sorted(
        (name for name, summary in summaries.items() if summary["all_round_gates_pass"]),
        key=lambda name: (
            summaries[name]["worst_case_minimax_score"]
            if summaries[name]["worst_case_minimax_score"] is not None else float("inf"),
            summaries[name]["max_absolute_land_share_delta"],
            -summaries[name]["family_balanced_median_ridge_contrast_gain_mm_yr"]
            if summaries[name]["family_balanced_median_ridge_contrast_gain_mm_yr"] is not None
            else float("inf"),
        ),
    )
    excluded_arms = sorted(set(summaries) - set(ranked_eligible_arms))
    freeze_eligible_arms = sorted(
        name for name, summary in summaries.items()
        if summary["all_round_gates_pass"]
        and summary["all_required_patterns_evaluated"]
        and summary["all_tuning_pattern_checks_pass"]
    )
    return ranked_eligible_arms, excluded_arms, freeze_eligible_arms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="+", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    reports = [json.loads(path.read_text()) for path in args.results]
    round_number = _validate_report_headers(reports, parser)
    seen, by_arm, settings = _collect_round_results(reports, round_number, parser)
    actual_lengths = _validate_round_settings(round_number, seen, settings, parser)
    summaries = {name: _candidate_summary(records) for name, records in sorted(by_arm.items())}
    comparisons = _build_comparisons(summaries, round_number, actual_lengths)
    ranked_eligible_arms, excluded_arms, freeze_eligible_arms = _rank_arms(summaries)

    output = {
        "tuning_round": round_number,
        "heldout": False,
        "solver": "finite_volume",
        "density": 4,
        "protocol_id": reports[0]["source_provenance"]["protocol_id"],
        "baseline_reference": "same-run ridge_025 arm; within-run comparisons only",
        "code_commit": reports[0]["code_commit"],
        "source_provenance": reports[0]["source_provenance"],
        "save_digests": seen,
        "settings": settings,
        "candidate_summaries": summaries,
        "coast_vs_matched_uniform": comparisons,
        "ranked_eligible_arms_for_tuning": ranked_eligible_arms,
        "excluded_arms": excluded_arms,
        "freeze_eligible_arms": freeze_eligible_arms,
        "scorer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

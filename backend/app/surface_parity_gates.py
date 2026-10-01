"""Quality gates for line-vs-quad parity runs (issue #247) -- see `surface_parity` for what a run
records and docs/surface-parity.md for the tolerance table this module implements.

Every gate result is `pass`, `warn`, `fail`, `insufficient` (not enough data to judge, e.g. an
ensemble gate with fewer than `MIN_ENSEMBLE_SEEDS` seeds) or `info` (a line-baseline
finding that does not gate the quad surface). The verdict is `fail` if any quad or paired
gate fails, else `warn` if any warns, else `pass`; `insufficient` and `info` results are listed
in the report but do not change the verdict.

Comparisons are statistical, never byte-for-byte: line-index-derived random streams make
line and quad trajectories diverge from the first step (docs/plate-surface-baseline.md 4.4).
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

from .surface_parity import DEFORM_TOPOLOGY_PHASES

PASS, WARN, FAIL, INSUFFICIENT, INFO = "pass", "warn", "fail", "insufficient", "info"
_RANK = {PASS: 0, INFO: 1, INSUFFICIENT: 2, WARN: 3, FAIL: 4}

MIN_ENSEMBLE_SEEDS = 3

# Hard invariants: audit violation kinds per gate. Neutral gates apply to both surfaces;
# `quad_only` ones have no line counterpart.
HARD_GATES = (
    ("H1", "quad_topology", "Quad leaf topology valid: no overlapping or duplicate leaves, 2:1 balanced", ("quad_topology",), True),
    ("H2", "quad_neighbours", "Quad neighbours in range, symmetric, and sharing an edge", ("quad_neighbours",), True),
    ("H3", "quad_folded", "No folded quad cells", ("quad_folded",), True),
    ("H4", "fields", "Every field finite, correctly sized and within its caps; areas positive", ("fields",), False),
    ("H5", "frames", "Plate frames are proper rotations", ("frame",), False),
    ("H6", "derived_caches", "Every populated derived cache equals a fresh rebuild; atmosphere grid fixed", ("stale_plate_cache", "stale_world_cache", "atmosphere_grid"), False),
    ("H7", "revisions", "Local node-set / frame changes always bump topology / geometry revisions", ("revision",), False),
)
# Elevation/Hc/Hm caps: the engine clamps them in specific code paths only, and the line
# surface breaches them too, so a quad breach warns rather than fails.
BOUNDS_GATE = ("H11", "field_caps", "Elevation, Hc and Hm stay within their caps", ("field_bounds",))

# (id, name, metric path under checkpoint["coverage"], warn above line + this, fail above line + this)
COVERAGE_GATES = (
    ("C1", "uncovered", "uncovered", None, 0.005),
    ("C4", "nodes_inside_other_plate", "nodes_inside_other_plate", 0.0, 0.005),
)
VOID_FAIL = 0.0005

# Quad plates carry exact cells on independently rotated lattices. Where two such lattices
# meet, filling the last seam with whole cells cannot in general make both cell edges coincide:
# leaving the cell out creates uncovered ground, while inserting it creates a narrow, real
# overlap. The #249 long/stress campaigns found that overlap stabilises near 1.9% while quad
# coverage and stacking remain substantially better than lines. Keep a warning band so a
# change toward that envelope is visible, but gate the exact-cell representation against its
# own documented tolerance rather than the line outline's unlike polygon measurement (#255).
QUAD_MULTIPLY_COVERED_WARN = 0.015
QUAD_MULTIPLY_COVERED_FAIL = 0.020

CONSERVATION_WARN, CONSERVATION_FAIL = 0.05, 0.10
AREA_ACCOUNTING_WARN = 0.01

ANISOTROPIC_WARN = 0.01
ALIGNMENT_WARN = 0.3
ALIGNMENT_MIN_ANISOTROPIC = 0.001
THIN_WARN = 0.01
SHAPE_WARN = 0.001

# Per-step stability: std of first differences over the second half of the run, quad must stay
# within `STABILITY_FACTOR` x line + the series' own absolute slack.
STABILITY_FACTOR = 2.0
STABILITY_SLACK = {"air_temperature_mean_c": 0.05, "ocean_temperature_mean_c": 0.05, "precipitation_mean_mm": 5.0, "sea_level_m": 1.0, "land_fraction": 0.002}

PERFORMANCE_WARN = 1.25

HEALPIX_SAME_NODE_SLACK = 0.05
HEALPIX_DISTANCE_FACTOR = 1.25

# Ensemble parity (docs/plate-surface-baseline.md 3.4): quad mean within 2 sigma of the line
# ensemble's seed-to-seed spread.
PARITY_SIGMA = 2.0
PARITY_METRICS = (
    ("land_fraction", ("totals", "land_fraction")),
    ("continental_hc_volume_km3", ("totals", "continental_hc_volume_km3")),
    ("hc_volume_km3", ("totals", "hc_volume_km3")),
    ("plates", ("totals", "plates")),
    ("sea_level_m", ("totals", "sea_level_m")),
    ("elevation_p05", ("totals", "elevation", "p05")),
    ("elevation_p50", ("totals", "elevation", "p50")),
    ("elevation_p95", ("totals", "elevation", "p95")),
)

ISSUE147_FRAMES = Path(__file__).resolve().parents[2] / "analysis" / "issue147-profile-20260922" / "frames.csv"
# frames.csv column -> harness phase bucket(s).
ISSUE147_PHASES = {
    "plate_shift": ("shift",),
    "plate_deform": ("deform",),
    "faults": ("faults",),
    "topology_changes": ("topology",),
    "climate_erosion_hydrology": ("climate_erosion_hydrology",),
    "volcanism": ("volcanism",),
    "sea_level": ("sea_level",),
}


def _get(document: dict | None, path: tuple[str, ...]):
    for key in path:
        if document is None:
            return None
        document = document.get(key)
    return document


def _result(gate: str, name: str, status: str, **details) -> dict:
    return {"gate": gate, "name": name, "status": status, **{k: v for k, v in details.items() if v is not None}}


def _worst(statuses) -> str:
    statuses = list(statuses)
    return max(statuses, key=lambda s: _RANK[s]) if statuses else PASS


def _myr(checkpoint: dict) -> float:
    return round(float(checkpoint["elapsed_myr"]), 6)


# --- Per-run gates --------------------------------------------------------------------------


def run_gates(run: dict) -> list[dict]:
    """Hard invariants for one run. Line-surface failures are reported as `info` (baseline
    findings), never as a quad gate failure."""
    surface = run["surface"]
    by_kind = defaultdict(list)
    for violation in run["violations"]:
        by_kind[violation["kind"]].append(violation)
    failed = FAIL if surface == "quad" else INFO
    results = []
    for gate, name, description, kinds, quad_only in HARD_GATES:
        if quad_only and surface != "quad":
            continue
        found = [v for kind in kinds for v in by_kind.get(kind, [])]
        if gate == "H6":
            mismatches = [c["step"] for c in run["checkpoints"] if (_get(c, ("index_parity", "cached_kdtree_mismatches")) or 0) > 0]
            found += [{"kind": "cached_kdtree_mismatch", "step": step} for step in mismatches]
        results.append(
            _result(
                gate,
                name,
                failed if found else PASS,
                description=description,
                violations=len(found),
                first=found[0] if found else None,
                audits=run["audits"],
            )
        )
    gate, name, description, kinds = BOUNDS_GATE
    found = [v for kind in kinds for v in by_kind.get(kind, [])]
    results.append(_result(gate, name, (WARN if surface == "quad" else INFO) if found else PASS, description=description, violations=len(found), first=found[0] if found else None))
    loads = run.get("load_checks") or []
    if loads:
        bad = [c for c in loads if not c["state_identical"] or c["derived_state_after_load"]]
        results.append(_result("H8", "load_round_trip", failed if bad else PASS, description="Save/load reproduces authoritative state and drops every derived index", checks=len(loads), first=bad[0] if bad else None))
        continuation = [c["continuation_identical"] for c in loads if "continuation_identical" in c]
        if continuation:
            results.append(_result("H9", "load_continuation", PASS if all(continuation) else failed, description="A step from the loaded world matches a step from the in-memory world"))
    if surface == "quad":
        stacked = [(_myr(c), c["sample_cloud"]["stacked"]) for c in run["checkpoints"] if c["sample_cloud"]["stacked"] > 0]
        results.append(_result("H10", "quad_no_stacked_nodes", FAIL if stacked else PASS, description="No two quad nodes of one plate closer than 0.5 spacing", first=stacked[0] if stacked else None))
        for checkpoint in run["checkpoints"]:
            age = _myr(checkpoint)
            multiply = checkpoint["coverage"]["multiply_covered"]
            if multiply is None:
                results.append(_result("C3", "multiply_covered", INSUFFICIENT, age_myr=age))
                continue
            status = FAIL if multiply > QUAD_MULTIPLY_COVERED_FAIL else WARN if multiply > QUAD_MULTIPLY_COVERED_WARN else PASS
            results.append(
                _result(
                    "C3", "multiply_covered", status, age_myr=age, quad=multiply,
                    warn_above=QUAD_MULTIPLY_COVERED_WARN,
                    fail_above=QUAD_MULTIPLY_COVERED_FAIL,
                )
            )
    finite = []
    for checkpoint in run["checkpoints"]:
        climate = checkpoint["climate_hydrology"]
        values = [v for v in climate["stats"].values() if v is not None]
        if any(isinstance(v, str) or not np.isfinite(v) for v in values) or (climate.get("hydrology") and not climate["hydrology"]["finite"]):
            finite.append(_myr(checkpoint))
    results.append(_result("S1", "climate_hydrology_finite", failed if finite else PASS, description="Climate stats and hydrology fields finite at every checkpoint", first=finite[0] if finite else None))
    for result in results:
        result["surface"] = surface
        result["seed"] = run["seed"]
    return results


# --- Paired gates (one seed, quad vs lines) --------------------------------------------------


def paired_gates(lines: dict, quad: dict, line_timings: dict | None, quad_timings: dict | None) -> list[dict]:
    results = []
    line_by_age = {_myr(c): c for c in lines["checkpoints"]}
    quad_by_age = {_myr(c): c for c in quad["checkpoints"]}
    ages = sorted(set(line_by_age) & set(quad_by_age))
    line_start, quad_start = line_by_age.get(ages[0]) if ages else None, quad_by_age.get(ages[0]) if ages else None
    for age in ages:
        line, q = line_by_age[age], quad_by_age[age]
        results += _coverage_gates(line, q, age)
        if age > 0:
            results += _conservation_gates(line_start, line, quad_start, q, age)
        results += _mesh_gates(line, q, age)
        results += _healpix_gates(line, q, age)
    results += _stability_gates(lines, quad)
    if line_timings and quad_timings:
        results += _performance_gates(line_timings, quad_timings)
    for result in results:
        result["seed"] = quad["seed"]
    return results


def _coverage_gates(line: dict, quad: dict, age: float) -> list[dict]:
    results = []
    for gate, name, key, warn_margin, fail_margin in COVERAGE_GATES:
        l, q = line["coverage"][key], quad["coverage"][key]
        if l is None or q is None:
            results.append(_result(gate, name, INSUFFICIENT, age_myr=age))
            continue
        status = FAIL if q > l + fail_margin else WARN if warn_margin is not None and q > l + warn_margin else PASS
        results.append(_result(gate, name, status, age_myr=age, quad=q, lines=l, fail_above=l + fail_margin))
    void = quad["coverage"]["void"]
    results.append(_result("C2", "void", FAIL if void is not None and void > VOID_FAIL else PASS, age_myr=age, quad=void, lines=line["coverage"]["void"], fail_above=VOID_FAIL))
    return results


def _drift(start: dict, now: dict, key: str) -> float | None:
    before, after = start["totals"].get(key), now["totals"].get(key)
    if not before:
        return None
    return after / before - 1.0


def _conservation_gates(line_start: dict, line: dict, quad_start: dict, quad: dict, age: float) -> list[dict]:
    results = []
    for gate, key in (("K1", "hc_volume_km3"), ("K2", "continental_hc_volume_km3"), ("K4", "hm_volume_km3")):
        l, q = _drift(line_start, line, key), _drift(quad_start, quad, key)
        if l is None or q is None:
            results.append(_result(gate, f"{key}_drift", INSUFFICIENT, age_myr=age))
            continue
        gap = abs(q - l)
        status = FAIL if gap > CONSERVATION_FAIL else WARN if gap > CONSERVATION_WARN else PASS
        results.append(_result(gate, f"{key}_drift", status, age_myr=age, quad=q, lines=l, gap=gap, warn_above=CONSERVATION_WARN, fail_above=CONSERVATION_FAIL))
    totals, cover = quad["totals"], quad["coverage"]
    if totals.get("area_is_exact") and cover["uncovered"] is not None:
        # Exact cell areas: covered sphere + overlap must equal the summed cell area.
        expected = 1.0 - cover["uncovered"] + cover["multiply_covered"]
        error = abs(totals["area_over_sphere"] - expected)
        results.append(_result("K3", "quad_area_accounting", WARN if error > AREA_ACCOUNTING_WARN else PASS, age_myr=age, quad=totals["area_over_sphere"], expected=expected, warn_above=AREA_ACCOUNTING_WARN))
    return results


def _mesh_gates(line: dict, quad: dict, age: float) -> list[dict]:
    cloud, lattice = quad["sample_cloud"], quad["quad_lattice"] or {}
    results = [
        _result("M1", "anisotropic_fraction", WARN if cloud["anisotropic"] > ANISOTROPIC_WARN else PASS, age_myr=age, quad=cloud["anisotropic"], lines=line["sample_cloud"]["anisotropic"], warn_above=ANISOTROPIC_WARN)
    ]
    alignment = cloud["anisotropic_row_alignment"]
    if alignment is None or cloud["anisotropic"] < ALIGNMENT_MIN_ANISOTROPIC:
        results.append(_result("M2", "row_alignment", PASS, age_myr=age, quad=alignment, note="too few anisotropic nodes to have a direction"))
    else:
        results.append(_result("M2", "row_alignment", WARN if abs(alignment) > ALIGNMENT_WARN else PASS, age_myr=age, quad=alignment, lines=line["sample_cloud"]["anisotropic_row_alignment"], warn_above=ALIGNMENT_WARN))
    thin, line_thin = cloud["thin"], line["sample_cloud"]["thin"]
    results.append(_result("M3", "thin_fraction", WARN if thin > THIN_WARN or thin > line_thin else PASS, age_myr=age, quad=thin, lines=line_thin, warn_above=min(THIN_WARN, line_thin)))
    for gate, key in (("M4", "aspect_gt_4"), ("M5", "skew_gt_45")):
        value = lattice.get(key, 0.0)
        results.append(_result(gate, key, WARN if value > SHAPE_WARN else PASS, age_myr=age, quad=value, warn_above=SHAPE_WARN))
    stray = cloud["extra_components"] + cloud["isolated_nodes"]
    results.append(_result("M6", "fragments", WARN if stray else PASS, age_myr=age, quad=stray, lines=line["sample_cloud"]["extra_components"] + line["sample_cloud"]["isolated_nodes"]))
    return results


def _healpix_gates(line: dict, quad: dict, age: float) -> list[dict]:
    results = []
    line_categories = _get(line, ("index_parity", "healpix", "categories")) or {}
    quad_categories = _get(quad, ("index_parity", "healpix", "categories")) or {}
    for category, q in quad_categories.items():
        l = line_categories.get(category)
        if not l or not q["samples"] or not l["samples"]:
            results.append(_result("X1", f"healpix_{category}", INSUFFICIENT, age_myr=age, samples=q["samples"]))
            continue
        same_ok = q["same_node"] >= l["same_node"] - HEALPIX_SAME_NODE_SLACK
        distance_ok = q["distance_spacing"]["p95"] <= l["distance_spacing"]["p95"] * HEALPIX_DISTANCE_FACTOR
        results.append(
            _result(
                "X1",
                f"healpix_{category}",
                PASS if same_ok and distance_ok else WARN,
                age_myr=age,
                samples=q["samples"],
                quad_same_node=q["same_node"],
                lines_same_node=l["same_node"],
                quad_distance_p95=q["distance_spacing"]["p95"],
                lines_distance_p95=l["distance_spacing"]["p95"],
            )
        )
    return results


def _stability_gates(lines: dict, quad: dict) -> list[dict]:
    results = []
    for key, slack in STABILITY_SLACK.items():
        l, q = _jitter(lines["series"], key), _jitter(quad["series"], key)
        if l is None or q is None:
            results.append(_result("S2", f"stability_{key}", INSUFFICIENT))
            continue
        limit = STABILITY_FACTOR * l + slack
        results.append(_result("S2", f"stability_{key}", WARN if q > limit else PASS, quad=q, lines=l, warn_above=limit))
    return results


def _jitter(series: list[dict], key: str) -> float | None:
    values = np.array([row.get(key) for row in series[len(series) // 2 :] if isinstance(row.get(key), (int, float))], dtype=float)
    if len(values) < 3:
        return None
    return float(np.std(np.diff(values)))


def phase_means(timings: dict) -> dict[str, float]:
    steps = timings["steps"]
    if not steps:
        return {}
    totals = defaultdict(float)
    for step in steps:
        totals["step_total"] += step["wall_s"]
        for phase, seconds in step["phases"].items():
            totals[phase] += seconds
    totals["deform_topology"] = sum(totals[p] for p in DEFORM_TOPOLOGY_PHASES)
    return {phase: seconds / len(steps) for phase, seconds in totals.items()}


def _performance_gates(line_timings: dict, quad_timings: dict) -> list[dict]:
    line, quad = phase_means(line_timings), phase_means(quad_timings)
    results = []
    for gate, phase in (("R1", "deform_topology"), ("R2", "step_total")):
        if not line.get(phase):
            results.append(_result(gate, f"{phase}_s_per_step", INSUFFICIENT))
            continue
        ratio = quad.get(phase, 0.0) / line[phase]
        results.append(_result(gate, f"{phase}_s_per_step", WARN if ratio > PERFORMANCE_WARN else PASS, quad=quad.get(phase), lines=line[phase], ratio=ratio, warn_above=PERFORMANCE_WARN))
    return results


def issue147_phase_means(path: Path = ISSUE147_FRAMES) -> dict[str, float] | None:
    """Per-step phase means from the issue #147 profile (10 steps per frame). That run was
    under cProfile, so absolute times are inflated; compare ratios with care."""
    if not path.exists():
        return None
    with path.open() as f:
        rows = list(csv.DictReader(f))
    steps = 10 * len(rows)
    means = {column: sum(float(r[column]) for r in rows) / steps for column in ISSUE147_PHASES}
    means["step_total"] = sum(float(r["step_total"]) for r in rows) / steps
    return means


# --- Ensemble gates -------------------------------------------------------------------------


def ensemble_gates(pairs: dict[int, dict[str, dict]]) -> list[dict]:
    seeds = sorted(seed for seed, runs in pairs.items() if "lines" in runs and "quad" in runs)
    ages = None
    for seed in seeds:
        seed_ages = {_myr(c) for c in pairs[seed]["lines"]["checkpoints"]} & {_myr(c) for c in pairs[seed]["quad"]["checkpoints"]}
        ages = seed_ages if ages is None else ages & seed_ages
    results = []
    for age in sorted(ages or ()):
        if age == 0:
            continue
        for name, path in PARITY_METRICS:
            line_values, quad_values = [], []
            for seed in seeds:
                line = next(c for c in pairs[seed]["lines"]["checkpoints"] if _myr(c) == age)
                quad = next(c for c in pairs[seed]["quad"]["checkpoints"] if _myr(c) == age)
                line_values.append(_get(line, path))
                quad_values.append(_get(quad, path))
            if any(v is None for v in line_values + quad_values):
                results.append(_result("P1", name, INSUFFICIENT, age_myr=age))
                continue
            line_mean, quad_mean = float(np.mean(line_values)), float(np.mean(quad_values))
            if len(seeds) < MIN_ENSEMBLE_SEEDS:
                results.append(_result("P1", name, INSUFFICIENT, age_myr=age, seeds=len(seeds), lines_mean=line_mean, quad_mean=quad_mean, note=f"needs >= {MIN_ENSEMBLE_SEEDS} seeds"))
                continue
            spread = float(np.std(line_values, ddof=1))
            limit = PARITY_SIGMA * spread
            delta = quad_mean - line_mean
            status = FAIL if abs(delta) > limit and not np.isclose(delta, 0.0) else PASS
            results.append(_result("P1", name, status, age_myr=age, seeds=len(seeds), lines_mean=line_mean, quad_mean=quad_mean, lines_sigma=spread, delta=delta, fail_above=limit))
    return results


# --- Evaluation -----------------------------------------------------------------------------


def evaluate(runs: list[tuple[dict, dict | None]]) -> dict:
    """Judge a set of runs (as returned by `surface_parity.load_results`)."""
    pairs: dict[int, dict[str, dict]] = defaultdict(dict)
    timing_pairs: dict[int, dict[str, dict | None]] = defaultdict(dict)
    configs = {run["config"]["name"] for run, _ in runs}
    for run, timings in runs:
        pairs[run["seed"]][run["surface"]] = run
        timing_pairs[run["seed"]][run["surface"]] = timings
    results = []
    for run, _ in runs:
        results += run_gates(run)
    for seed, surfaces in sorted(pairs.items()):
        if "lines" in surfaces and "quad" in surfaces:
            results += paired_gates(surfaces["lines"], surfaces["quad"], timing_pairs[seed].get("lines"), timing_pairs[seed].get("quad"))
    results += ensemble_gates(pairs)

    gates: dict[str, dict] = {}
    for result in results:
        key = f"{result['gate']}:{result['name']}"
        summary = gates.setdefault(key, {"gate": result["gate"], "name": result["name"], "status": PASS, "results": 0, "worst": None})
        summary["results"] += 1
        if result["status"] == INFO:
            summary["line_baseline_findings"] = summary.get("line_baseline_findings", 0) + 1
        if summary["worst"] is None or _RANK[result["status"]] > _RANK[summary["worst"]["status"]]:
            summary["worst"] = result
        summary["status"] = _worst([summary["status"], result["status"]])
        if result.get("description"):
            summary["description"] = result["description"]
    verdict = _worst(r["status"] for r in results)
    return {
        "verdict": verdict if verdict in (FAIL, WARN) else PASS,
        "configs": sorted(configs),
        "seeds": sorted(pairs),
        "gates": dict(sorted(gates.items(), key=lambda item: (_gate_order(item[1]["gate"]), item[0]))),
        "line_baseline_findings": [r for r in results if r["status"] == INFO],
        "results": results,
    }


def _gate_order(gate: str) -> tuple[int, int]:
    groups = "HCKMSXRP"
    return groups.index(gate[0]) if gate[0] in groups else len(groups), int(gate[1:])


# --- Report ---------------------------------------------------------------------------------

GROUP_TITLES = {
    "H": "Hard invariants",
    "C": "Coverage and overlap",
    "K": "Conservation (actual node/cell areas)",
    "M": "Mesh quality (quad)",
    "S": "Climate and hydrology",
    "X": "Derived indexes: HEALPix vs KD-tree",
    "R": "Performance",
    "P": "Ensemble parity",
}


def _fmt(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        if value != 0 and (abs(value) < 1e-3 or abs(value) >= 1e6):
            return f"{value:.3g}"
        return f"{value:.4g}"
    return str(value)


def render_report(evaluation: dict, runs: list[tuple[dict, dict | None]], commands: list[str]) -> str:
    lines = [
        "# Line vs quad parity report",
        "",
        f"**Verdict: {evaluation['verdict'].upper()}** — configs {', '.join(evaluation['configs'])}; seeds {', '.join(map(str, evaluation['seeds']))}.",
        "",
        "Reproduce:",
        "",
        *[f"    {command}" for command in commands],
        "",
        "Tolerances and gate definitions: docs/surface-parity.md. `info` rows are line-baseline",
        "findings; they never gate the quad surface.",
        "",
    ]
    by_group = defaultdict(list)
    for key, gate in evaluation["gates"].items():
        by_group[gate["gate"][0]].append((key, gate))
    for group, title in GROUP_TITLES.items():
        if group not in by_group:
            continue
        lines += [f"## {title}", "", "| gate | status | worst case |", "|---|---|---|"]
        for key, gate in by_group[group]:
            worst = gate["worst"] or {}
            detail = ", ".join(
                f"{k}={_fmt(v)}"
                for k, v in worst.items()
                if k not in ("gate", "name", "status", "description", "first") and not isinstance(v, dict)
            )
            if worst.get("first"):
                detail += f"; first: {worst['first']}"
            status = gate["status"]
            if gate.get("line_baseline_findings"):
                status += f" (line baseline: {gate['line_baseline_findings']} finding(s))"
            lines.append(f"| {key} | {status} | {detail} |")
        lines.append("")

    lines += _checkpoint_tables(runs)
    lines += _performance_table(runs)
    return "\n".join(lines) + "\n"


_CHECKPOINT_ROWS = (
    ("plates", ("totals", "plates")),
    ("nodes", ("totals", "nodes")),
    ("land fraction", ("totals", "land_fraction")),
    ("sea level (m)", ("totals", "sea_level_m")),
    ("Hc volume (km³)", ("totals", "hc_volume_km3")),
    ("continental Hc (km³)", ("totals", "continental_hc_volume_km3")),
    ("elevation p05/p50/p95 (m)", None),
    ("uncovered", ("coverage", "uncovered")),
    ("multiply covered", ("coverage", "multiply_covered")),
    ("void", ("coverage", "void")),
    ("nodes inside other plate", ("coverage", "nodes_inside_other_plate")),
    ("stacked", ("sample_cloud", "stacked")),
    ("anisotropic", ("sample_cloud", "anisotropic")),
    ("row alignment", ("sample_cloud", "anisotropic_row_alignment")),
    ("thin", ("sample_cloud", "thin")),
    ("one-node lines", ("line_topology", "one_node_lines")),
    ("quad one-cell-thin", ("quad_lattice", "one_cell_thin")),
    ("air temperature mean (°C)", ("climate_hydrology", "stats", "air_temperature_mean_c")),
    ("precipitation mean (mm)", ("climate_hydrology", "stats", "precipitation_mean_mm")),
    ("river fraction of land", ("climate_hydrology", "hydrology", "river_fraction_of_land")),
    ("HEALPix same node", ("index_parity", "healpix", "categories", "all", "same_node")),
)


def _checkpoint_tables(runs) -> list[str]:
    pairs = defaultdict(dict)
    for run, _ in runs:
        pairs[run["seed"]][run["surface"]] = run
    out = []
    for seed, surfaces in sorted(pairs.items()):
        ordered = [s for s in ("lines", "quad") if s in surfaces]
        ages = sorted({_myr(c) for s in ordered for c in surfaces[s]["checkpoints"]})
        header = "| seed " + str(seed) + " | " + " | ".join(f"{a:g} Myr {s}" for a in ages for s in ordered) + " |"
        out += [f"## Checkpoints, seed {seed}", "", header, "|---" * (1 + len(ages) * len(ordered)) + "|"]
        for label, path in _CHECKPOINT_ROWS:
            cells = []
            for age in ages:
                for surface in ordered:
                    checkpoint = next((c for c in surfaces[surface]["checkpoints"] if _myr(c) == age), None)
                    if path is None:
                        e = _get(checkpoint, ("totals", "elevation")) or {}
                        cells.append("/".join(_fmt(e.get(k)) for k in ("p05", "p50", "p95")))
                    else:
                        cells.append(_fmt(_get(checkpoint, path)))
            out.append(f"| {label} | " + " | ".join(cells) + " |")
        out.append("")
    return out


def _performance_table(runs) -> list[str]:
    rows = [(run["seed"], run["surface"], phase_means(timings)) for run, timings in runs if timings]
    if not rows:
        return []
    reference = issue147_phase_means() if any(run["config"]["name"] == "issue147" for run, _ in runs) else None
    phases = sorted({p for _, _, means in rows for p in means}, key=lambda p: (p != "step_total", p != "deform_topology", p))
    header = "| s/step | " + " | ".join(f"seed {seed} {surface}" for seed, surface, _ in rows) + (" | #147 profile |" if reference else " |")
    out = ["## Performance (mean seconds per step)", "", header, "|---" * (1 + len(rows) + (1 if reference else 0)) + "|"]
    reverse = {bucket: column for column, buckets in ISSUE147_PHASES.items() for bucket in buckets}
    for phase in phases:
        cells = [_fmt(means.get(phase)) for _, _, means in rows]
        if reference:
            column = "step_total" if phase == "step_total" else reverse.get(phase)
            cells.append(_fmt(reference.get(column)) if column else "—")
        out.append(f"| {phase} | " + " | ".join(cells) + " |")
    out += ["", "Index builds per step (mean):", ""]
    for run, timings in runs:
        if not timings or not timings["steps"]:
            continue
        builds = defaultdict(int)
        for step in timings["steps"]:
            for name, count in step["index_builds"].items():
                builds[name] += count
        summary = ", ".join(f"{name} {count / len(timings['steps']):.1f}" for name, count in sorted(builds.items()))
        out.append(f"- seed {run['seed']} {run['surface']}: {summary}")
    if reference:
        out += ["", "The #147 profile ran under cProfile, so its absolute times are inflated."]
    out.append("")
    return out

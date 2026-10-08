"""Quality gates for plate-surface audit runs (issues #247/#249) -- see `surface_parity` for what
a run records and docs/surface-parity.md for the tolerance table this module implements.

Every gate result is `pass`, `warn`, `fail` or `insufficient` (not enough data to judge). The
verdict is `fail` if any gate fails, else `warn` if any warns, else `pass`; `insufficient`
results are listed in the report but do not change the verdict.

The line-vs-quad comparison gates these campaigns started with were retired with the line
surface (#251): every gate here has an absolute threshold.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from .surface_parity import DEFORM_TOPOLOGY_PHASES

PASS, WARN, FAIL, INSUFFICIENT = "pass", "warn", "fail", "insufficient"
_RANK = {PASS: 0, INSUFFICIENT: 1, WARN: 2, FAIL: 3}

# Hard invariants: audit violation kinds per gate.
HARD_GATES = (
    ("H1", "quad_topology", "Quad leaf topology valid: no overlapping or duplicate leaves, 2:1 balanced", ("quad_topology",)),
    ("H2", "quad_neighbours", "Quad neighbours in range, symmetric, and sharing an edge", ("quad_neighbours",)),
    ("H3", "quad_folded", "No folded quad cells", ("quad_folded",)),
    ("H4", "fields", "Every field finite, correctly sized and within its caps; areas positive", ("fields",)),
    ("H5", "frames", "Plate frames are proper rotations", ("frame",)),
    ("H6", "derived_caches", "Every populated derived cache equals a fresh rebuild; atmosphere grid fixed", ("stale_plate_cache", "stale_world_cache", "atmosphere_grid")),
    ("H7", "revisions", "Local node-set / frame changes always bump topology / geometry revisions", ("revision",)),
)
# Elevation/Hc/Hm caps: the engine clamps them in specific code paths only, so a breach warns
# rather than fails.
BOUNDS_GATE = ("H11", "field_caps", "Elevation, Hc and Hm stay within their caps", ("field_bounds",))

VOID_FAIL = 0.0005

# Quad plates carry exact cells on independently rotated lattices. Where two such lattices
# meet, filling the last seam with whole cells cannot in general make both cell edges coincide:
# leaving the cell out creates uncovered ground, while inserting it creates a narrow, real
# overlap. The #249 long/stress campaigns found that overlap stabilises near 1.9%. Keep a
# warning band so a change toward that envelope is visible (#255).
QUAD_MULTIPLY_COVERED_WARN = 0.015
QUAD_MULTIPLY_COVERED_FAIL = 0.020

AREA_ACCOUNTING_WARN = 0.01

ANISOTROPIC_WARN = 0.01
ALIGNMENT_WARN = 0.3
ALIGNMENT_MIN_ANISOTROPIC = 0.001
THIN_WARN = 0.01
SHAPE_WARN = 0.001


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
    """Every gate for one run: hard invariants across the run, then per-checkpoint coverage,
    area-accounting and mesh-quality gates."""
    by_kind = defaultdict(list)
    for violation in run["violations"]:
        by_kind[violation["kind"]].append(violation)
    results = []
    for gate, name, description, kinds in HARD_GATES:
        found = [v for kind in kinds for v in by_kind.get(kind, [])]
        if gate == "H6":
            mismatches = [c["step"] for c in run["checkpoints"] if (_get(c, ("index_parity", "cached_kdtree_mismatches")) or 0) > 0]
            found += [{"kind": "cached_kdtree_mismatch", "step": step} for step in mismatches]
        results.append(
            _result(
                gate,
                name,
                FAIL if found else PASS,
                description=description,
                violations=len(found),
                first=found[0] if found else None,
                audits=run["audits"],
            )
        )
    gate, name, description, kinds = BOUNDS_GATE
    found = [v for kind in kinds for v in by_kind.get(kind, [])]
    results.append(_result(gate, name, WARN if found else PASS, description=description, violations=len(found), first=found[0] if found else None))
    loads = run.get("load_checks") or []
    if loads:
        bad = [c for c in loads if not c["state_identical"] or c["derived_state_after_load"]]
        results.append(_result("H8", "load_round_trip", FAIL if bad else PASS, description="Save/load reproduces authoritative state and drops every derived index", checks=len(loads), first=bad[0] if bad else None))
        continuation = [c["continuation_identical"] for c in loads if "continuation_identical" in c]
        if continuation:
            results.append(_result("H9", "load_continuation", PASS if all(continuation) else FAIL, description="A step from the loaded world matches a step from the in-memory world"))
    stacked = [(_myr(c), c["sample_cloud"]["stacked"]) for c in run["checkpoints"] if c["sample_cloud"]["stacked"] > 0]
    results.append(_result("H10", "quad_no_stacked_nodes", FAIL if stacked else PASS, description="No two quad nodes of one plate closer than 0.5 spacing", first=stacked[0] if stacked else None))
    finite = []
    for checkpoint in run["checkpoints"]:
        climate = checkpoint["climate_hydrology"]
        values = [v for v in climate["stats"].values() if v is not None]
        if any(isinstance(v, str) or not np.isfinite(v) for v in values) or (climate.get("hydrology") and not climate["hydrology"]["finite"]):
            finite.append(_myr(checkpoint))
    results.append(_result("S1", "climate_hydrology_finite", FAIL if finite else PASS, description="Climate stats and hydrology fields finite at every checkpoint", first=finite[0] if finite else None))
    for checkpoint in run["checkpoints"]:
        age = _myr(checkpoint)
        results += _coverage_gates(checkpoint, age)
        results += _area_gates(checkpoint, age)
        results += _mesh_gates(checkpoint, age)
    for result in results:
        result["seed"] = run["seed"]
    return results


def _coverage_gates(checkpoint: dict, age: float) -> list[dict]:
    cover = checkpoint["coverage"]
    void = cover["void"]
    results = [_result("C2", "void", FAIL if void is not None and void > VOID_FAIL else PASS, age_myr=age, value=void, fail_above=VOID_FAIL)]
    multiply = cover["multiply_covered"]
    if multiply is None:
        results.append(_result("C3", "multiply_covered", INSUFFICIENT, age_myr=age))
    else:
        status = FAIL if multiply > QUAD_MULTIPLY_COVERED_FAIL else WARN if multiply > QUAD_MULTIPLY_COVERED_WARN else PASS
        results.append(_result("C3", "multiply_covered", status, age_myr=age, value=multiply, warn_above=QUAD_MULTIPLY_COVERED_WARN, fail_above=QUAD_MULTIPLY_COVERED_FAIL))
    return results


def _area_gates(checkpoint: dict, age: float) -> list[dict]:
    totals, cover = checkpoint["totals"], checkpoint["coverage"]
    if not totals.get("area_is_exact") or cover["uncovered"] is None:
        return []
    # Exact cell areas: covered sphere + overlap must equal the summed cell area.
    expected = 1.0 - cover["uncovered"] + cover["multiply_covered"]
    error = abs(totals["area_over_sphere"] - expected)
    return [_result("K3", "quad_area_accounting", WARN if error > AREA_ACCOUNTING_WARN else PASS, age_myr=age, value=totals["area_over_sphere"], expected=expected, warn_above=AREA_ACCOUNTING_WARN)]


def _mesh_gates(checkpoint: dict, age: float) -> list[dict]:
    cloud, lattice = checkpoint["sample_cloud"], checkpoint["quad_lattice"] or {}
    results = [_result("M1", "anisotropic_fraction", WARN if cloud["anisotropic"] > ANISOTROPIC_WARN else PASS, age_myr=age, value=cloud["anisotropic"], warn_above=ANISOTROPIC_WARN)]
    alignment = cloud["anisotropic_row_alignment"]
    if alignment is None or cloud["anisotropic"] < ALIGNMENT_MIN_ANISOTROPIC:
        results.append(_result("M2", "row_alignment", PASS, age_myr=age, value=alignment, note="too few anisotropic nodes to have a direction"))
    else:
        results.append(_result("M2", "row_alignment", WARN if abs(alignment) > ALIGNMENT_WARN else PASS, age_myr=age, value=alignment, warn_above=ALIGNMENT_WARN))
    results.append(_result("M3", "thin_fraction", WARN if cloud["thin"] > THIN_WARN else PASS, age_myr=age, value=cloud["thin"], warn_above=THIN_WARN))
    for gate, key in (("M4", "aspect_gt_4"), ("M5", "skew_gt_45")):
        value = lattice.get(key, 0.0)
        results.append(_result(gate, key, WARN if value > SHAPE_WARN else PASS, age_myr=age, value=value, warn_above=SHAPE_WARN))
    stray = cloud["extra_components"] + cloud["isolated_nodes"]
    results.append(_result("M6", "fragments", WARN if stray else PASS, age_myr=age, value=stray))
    return results


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


# --- Evaluation -----------------------------------------------------------------------------


def evaluate(runs: list[tuple[dict, dict | None]]) -> dict:
    """Judge a set of runs (as returned by `surface_parity.load_results`)."""
    configs = {run["config"]["name"] for run, _ in runs}
    results = []
    for run, _ in runs:
        results += run_gates(run)

    gates: dict[str, dict] = {}
    for result in results:
        key = f"{result['gate']}:{result['name']}"
        summary = gates.setdefault(key, {"gate": result["gate"], "name": result["name"], "status": PASS, "results": 0, "worst": None})
        summary["results"] += 1
        if summary["worst"] is None or _RANK[result["status"]] > _RANK[summary["worst"]["status"]]:
            summary["worst"] = result
        summary["status"] = _worst([summary["status"], result["status"]])
        if result.get("description"):
            summary["description"] = result["description"]
    verdict = _worst(r["status"] for r in results)
    return {
        "verdict": verdict if verdict in (FAIL, WARN) else PASS,
        "configs": sorted(configs),
        "seeds": sorted({run["seed"] for run, _ in runs}),
        "gates": dict(sorted(gates.items(), key=lambda item: (_gate_order(item[1]["gate"]), item[0]))),
        "results": results,
    }


def _gate_order(gate: str) -> tuple[int, int]:
    groups = "HCKMS"
    return groups.index(gate[0]) if gate[0] in groups else len(groups), int(gate[1:])


# --- Report ---------------------------------------------------------------------------------

GROUP_TITLES = {
    "H": "Hard invariants",
    "C": "Coverage and overlap",
    "K": "Conservation (actual node/cell areas)",
    "M": "Mesh quality (quad)",
    "S": "Climate and hydrology",
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
        "# Plate-surface audit report",
        "",
        f"**Verdict: {evaluation['verdict'].upper()}** — configs {', '.join(evaluation['configs'])}; seeds {', '.join(map(str, evaluation['seeds']))}.",
        "",
        "Reproduce:",
        "",
        *[f"    {command}" for command in commands],
        "",
        "Tolerances and gate definitions: docs/surface-parity.md.",
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
            lines.append(f"| {key} | {gate['status']} | {detail} |")
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
    ("quad one-cell-thin", ("quad_lattice", "one_cell_thin")),
    ("air temperature mean (°C)", ("climate_hydrology", "stats", "air_temperature_mean_c")),
    ("precipitation mean (mm)", ("climate_hydrology", "stats", "precipitation_mean_mm")),
    ("river fraction of land", ("climate_hydrology", "hydrology", "river_fraction_of_land")),
    ("HEALPix same node", ("index_parity", "healpix", "categories", "all", "same_node")),
)


def _checkpoint_tables(runs) -> list[str]:
    out = []
    for run, _ in sorted(runs, key=lambda pair: pair[0]["seed"]):
        ages = sorted({_myr(c) for c in run["checkpoints"]})
        header = f"| seed {run['seed']} | " + " | ".join(f"{a:g} Myr" for a in ages) + " |"
        out += [f"## Checkpoints, seed {run['seed']}", "", header, "|---" * (1 + len(ages)) + "|"]
        for label, path in _CHECKPOINT_ROWS:
            cells = []
            for age in ages:
                checkpoint = next((c for c in run["checkpoints"] if _myr(c) == age), None)
                if path is None:
                    e = _get(checkpoint, ("totals", "elevation")) or {}
                    cells.append("/".join(_fmt(e.get(k)) for k in ("p05", "p50", "p95")))
                else:
                    cells.append(_fmt(_get(checkpoint, path)))
            out.append(f"| {label} | " + " | ".join(cells) + " |")
        out.append("")
    return out


def _performance_table(runs) -> list[str]:
    rows = [(run["seed"], phase_means(timings)) for run, timings in runs if timings]
    if not rows:
        return []
    phases = sorted({p for _, means in rows for p in means}, key=lambda p: (p != "step_total", p != "deform_topology", p))
    header = "| s/step | " + " | ".join(f"seed {seed}" for seed, _ in rows) + " |"
    out = ["## Performance (mean seconds per step)", "", header, "|---" * (1 + len(rows)) + "|"]
    for phase in phases:
        out.append(f"| {phase} | " + " | ".join(_fmt(means.get(phase)) for _, means in rows) + " |")
    out += ["", "Index builds per step (mean):", ""]
    for run, timings in runs:
        if not timings or not timings["steps"]:
            continue
        builds = defaultdict(int)
        for step in timings["steps"]:
            for name, count in step["index_builds"].items():
                builds[name] += count
        summary = ", ".join(f"{name} {count / len(timings['steps']):.1f}" for name, count in sorted(builds.items()))
        out.append(f"- seed {run['seed']}: {summary}")
    out.append("")
    return out

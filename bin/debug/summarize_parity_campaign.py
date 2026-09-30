#!/usr/bin/env python3
"""Issue #249: summarize a parity campaign made of several `surface_parity.py` output
directories (one per run set) into Markdown tables -- the gate status of every set, and
seed-averaged line vs quad metrics at each shared checkpoint.

    cd backend
    .venv/bin/python ../bin/debug/summarize_parity_campaign.py ../analysis/issue249-campaign > ../analysis/issue249-campaign/summary.md

Every immediate subdirectory holding `seed*-*.json` runs is one set; each is re-judged with
`surface_parity_gates.evaluate`, so the tables always reflect the current gate code. Per-phase
timing means are written to each set's `phase-means.json` and read back from there when the
bulky `*.timings.json` files have not been kept (the performance gates then read `insufficient`
on re-judging, but the performance tables survive).
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import surface_parity, surface_parity_gates as gates  # noqa: E402

METRICS = (
    ("plates", ("totals", "plates")),
    ("land fraction", ("totals", "land_fraction")),
    ("sea level (m)", ("totals", "sea_level_m")),
    ("Hc drift since 0 Myr", None),
    ("uncovered", ("coverage", "uncovered")),
    ("multiply covered", ("coverage", "multiply_covered")),
    ("nodes inside other plate", ("coverage", "nodes_inside_other_plate")),
    ("stacked", ("sample_cloud", "stacked")),
    ("anisotropic", ("sample_cloud", "anisotropic")),
    ("thin", ("sample_cloud", "thin")),
    ("one-node lines", ("line_topology", "one_node_lines")),
    ("quad one-cell-thin", ("quad_lattice", "one_cell_thin")),
    ("quad hole loops", ("quad_lattice", "hole_loops")),
    ("air temperature (°C)", ("climate_hydrology", "stats", "air_temperature_mean_c")),
    ("precipitation (mm)", ("climate_hydrology", "stats", "precipitation_mean_mm")),
    ("river fraction of land", ("climate_hydrology", "hydrology", "river_fraction_of_land")),
    ("HEALPix same node", ("index_parity", "healpix", "categories", "all", "same_node")),
    ("HEALPix p95 distance (s)", ("index_parity", "healpix", "categories", "all", "distance_spacing", "p95")),
)


def _get(document, path):
    for key in path:
        if document is None:
            return None
        document = document.get(key)
    return document


def _value(run: dict, checkpoint: dict, label: str, path):
    if path is None:
        start = run["checkpoints"][0]["totals"]["hc_volume_km3"]
        return checkpoint["totals"]["hc_volume_km3"] / start - 1.0
    return _get(checkpoint, path)


def _fmt(value) -> str:
    return "—" if value is None else gates._fmt(float(value))


def phase_means(directory: Path, runs) -> dict:
    """{surface: {seed: {phase: mean s/step}}} from the runs' timings files, written to
    `phase-means.json` so the (large) timings files need not be kept; read back from there
    when they are gone."""
    path = directory / "phase-means.json"
    perf = defaultdict(dict)
    for run, timings in runs:
        if timings:
            perf[run["surface"]][str(run["seed"])] = gates.phase_means(timings)
    if perf:
        path.write_text(surface_parity.dumps(perf))
        return perf
    return json.loads(path.read_text()) if path.exists() else {}


def summarize(root: Path) -> str:
    sets = sorted(p for p in root.iterdir() if p.is_dir() and any(p.glob("seed*-*.json")))
    evaluations = {}
    runs_by_set = {}
    for directory in sets:
        runs = surface_parity.load_results(directory)
        runs_by_set[directory.name] = runs
        evaluations[directory.name] = gates.evaluate(runs)

    out = ["# Parity campaign summary", "", "| set | runs | verdict |", "|---|---:|---|"]
    for name, evaluation in evaluations.items():
        out.append(f"| {name} | {len(runs_by_set[name])} | {evaluation['verdict']} |")

    keys = sorted({k for e in evaluations.values() for k in e["gates"]}, key=lambda k: (gates._gate_order(k.split(":")[0]), k))
    out += ["", "## Gate status by set", "", "| gate | " + " | ".join(evaluations) + " |", "|---" * (1 + len(evaluations)) + "|"]
    for key in keys:
        cells = []
        for evaluation in evaluations.values():
            gate = evaluation["gates"].get(key)
            if gate is None:
                cells.append("—")
                continue
            cell = gate["status"]
            if gate.get("line_baseline_findings"):
                cell += f" (+{gate['line_baseline_findings']} line)"
            cells.append(cell)
        out.append(f"| {key} | " + " | ".join(cells) + " |")

    for name, runs in runs_by_set.items():
        by_surface = defaultdict(list)
        for run, _ in runs:
            by_surface[run["surface"]].append(run)
        surfaces = [s for s in ("lines", "quad") if s in by_surface]
        ages = sorted(set.intersection(*[{gates._myr(c) for r in by_surface[s] for c in r["checkpoints"]} for s in surfaces]))
        seeds = sorted({r["seed"] for r, _ in runs})
        out += ["", f"## {name}: seed means (seeds {', '.join(map(str, seeds))})", ""]
        out.append("| metric | " + " | ".join(f"{a:g} Myr {s}" for a in ages for s in surfaces) + " |")
        out.append("|---" * (1 + len(ages) * len(surfaces)) + "|")
        for label, path in METRICS:
            cells = []
            for age in ages:
                for surface in surfaces:
                    values = []
                    for run in by_surface[surface]:
                        checkpoint = next((c for c in run["checkpoints"] if gates._myr(c) == age), None)
                        value = _value(run, checkpoint, label, path) if checkpoint else None
                        if isinstance(value, (int, float)):
                            values.append(value)
                    cells.append(_fmt(np.mean(values)) if values else "—")
            out.append(f"| {label} | " + " | ".join(cells) + " |")
        perf = phase_means(root / name, runs)
        if perf.get("lines") and perf.get("quad"):
            out += ["", "| s/step (seed mean) | lines | quad | quad / lines |", "|---|---:|---:|---:|"]
            for phase in ("step_total", "deform_topology", "deform", "topology", "gap_fill", "overlap_tracking", "faults", "shift", "climate_erosion_hydrology", "magma_transport"):
                line = np.mean([m.get(phase, 0.0) for m in perf["lines"].values()])
                quad = np.mean([m.get(phase, 0.0) for m in perf["quad"].values()])
                out.append(f"| {phase} | {_fmt(line)} | {_fmt(quad)} | {_fmt(quad / line) if line else '—'} |")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    print(summarize(Path(sys.argv[1])), end="")

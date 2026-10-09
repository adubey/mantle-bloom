#!/usr/bin/env python3
"""Issue #248: convert real line-backed `.mbworld` saves to sparse quads and judge the result.

For each save this:

1. reads it through the legacy reader (`legacy_conversion.legacy_unpickler`), converts it
   (`legacy_conversion.convert_world_to_quad`) and writes the full `ConversionReport`;
2. round-trips the converted world through the current save format (which applies the usual
   load-time backfills), runs the #247 hard-invariant audit (`surface_parity.audit_world`) on
   it, and checks a second round trip reproduces its authoritative state hash;
3. continues it `--steps` steps of `--step-years`, recording totals and audits for each step;
4. checks the docs/save-compatibility.md tolerances and records each pass/fail.

`--render` writes elevation and platesDetail PNGs of the converted world.
`--write-converted` also writes each converted world, before any step, as
`<save>-quad.mbworld` -- the same file loading the save in the app and then saving produces.
Writes `<save>.json` per save plus `summary.json` and `summary.md`.

`--recheck` skips converting and re-judges the per-save JSON already in `--out` (after a
tolerance change, say).

Usage (from repo root, this repo's venv):
    backend/.venv/bin/python bin/debug/convert_legacy_saves.py SAVE.mbworld [...] \\
        --out analysis/issue248/conversions [--steps 3] [--render] [--jobs 4]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))


from app import legacy_conversion, persistence, surface_parity  # noqa: E402
from app.surface_parity_gates import QUAD_MULTIPLY_COVERED_FAIL, QUAD_MULTIPLY_COVERED_WARN  # noqa: E402
from app import world as world_mod  # noqa: E402

# docs/save-compatibility.md §5. Volume deltas are against the line world clipped to the
# Hc/Hm caps (the line engine's own end-of-step clamp removed anything past them on its next
# step), with node areas deduplicated; area fractions use those same areas.
TOLERANCES = {
    # Volume the transfer itself gained or lost: the quad total against the capped line total,
    # less what the caps removed beyond the line world's own over-cap volume.
    "crust_transfer_delta": 0.01,
    "mantle_transfer_delta": 0.015,
    # Overlap crust stacked past the caps, after spreading, and removed (suture delamination).
    "crust_delamination": 0.01,
    "mantle_delamination": 0.015,
    "plate_crust_volume_delta": 0.05,  # worst plate with >= 500 nodes that no overlap touched
    "land_fraction_delta": 0.01,  # absolute, sampled
    "continental_fraction_delta": 0.01,  # absolute, sampled
    # Rotated lattices overlap at plate edges by construction; the #247 harness's C3 fail level
    # (results between its warn and fail levels are listed as warnings).
    "multiply_covered": QUAD_MULTIPLY_COVERED_FAIL,
    "sea_level_jump_m": 1e-6,  # at conversion
}
PLATE_DELTA_MIN_NODES = 500
# A plate counts as untouched by overlaps when fewer than this share of its nodes stacked onto
# other plates or were stacked onto it.
PLATE_DELTA_MAX_OVERLAP = 0.01


def _checks(report: dict, audit: list, round_trip: bool, step_audits: list) -> dict:
    ext = report["extensive"]
    after = report["coverage"]["after"]
    untouched = [
        p
        for p in report["plates"]
        if p["line_nodes"] >= PLATE_DELTA_MIN_NODES and p["stacked_nodes"] + p["stacked_onto"] < PLATE_DELTA_MAX_OVERLAP * p["line_nodes"]
    ]
    worst_plate = max((abs(p["crust_volume_delta"]) for p in untouched), default=0.0)
    fractions = report["area_fractions"]
    delamination = {name: max(report["clamped"][name] - report["source_over_cap"][name], 0.0) for name in report["clamped"]}
    checks = {
        "crust_transfer_delta": (abs(ext["crustal_thickness_m"]["delta_vs_capped"] + delamination["crustal_thickness_m"]), TOLERANCES["crust_transfer_delta"]),
        "mantle_transfer_delta": (
            abs(ext["mantle_lithosphere_thickness_m"]["delta_vs_capped"] + delamination["mantle_lithosphere_thickness_m"]),
            TOLERANCES["mantle_transfer_delta"],
        ),
        "crust_delamination": (delamination["crustal_thickness_m"], TOLERANCES["crust_delamination"]),
        "mantle_delamination": (delamination["mantle_lithosphere_thickness_m"], TOLERANCES["mantle_delamination"]),
        "plate_crust_volume_delta": (worst_plate, TOLERANCES["plate_crust_volume_delta"]),
        "land_fraction_delta": (abs(fractions["land"]["after"] - fractions["land"]["before"]), TOLERANCES["land_fraction_delta"]),
        "continental_fraction_delta": (abs(fractions["continental"]["after"] - fractions["continental"]["before"]), TOLERANCES["continental_fraction_delta"]),
        "multiply_covered": (after["multiply_covered"], TOLERANCES["multiply_covered"]),
        "sea_level_jump_m": (abs(report["sea_level"]["after_m"] - report["sea_level"]["before_m"]), TOLERANCES["sea_level_jump_m"]),
    }
    result = {name: {"value": value, "limit": limit, "pass": bool(value <= limit)} for name, (value, limit) in checks.items()}
    exact = {
        "audit_clean": not audit,
        "round_trip_identical": round_trip,
        "steps_audit_clean": all(not a for a in step_audits),
        "provenance_kept": not any(report["preserved"]["provenance_violations"].values()),
    }
    result.update({name: {"pass": bool(ok)} for name, ok in exact.items()})
    result["multiply_covered"]["warn"] = bool(after["multiply_covered"] > QUAD_MULTIPLY_COVERED_WARN)
    return result


def _read_line_world(data: bytes):
    """The save's `World`, its plates still the legacy reader's stand-ins -- what
    `persistence.load_world_bytes` reads before it converts."""
    payload = legacy_conversion.legacy_unpickler(data).load()
    world = payload["world"] if isinstance(payload, dict) else payload
    if not legacy_conversion.is_line_world(world):
        raise ValueError("not a line-backed save")
    return world


def convert_one(path: Path, out_dir: Path, steps: int, step_years: float, render: bool, write_converted: bool) -> dict:
    started = time.perf_counter()
    lines = _read_line_world(path.read_bytes())
    # The save's own age: the converted world is stepped further below.
    elapsed_years = lines.elapsed_years
    load_s = time.perf_counter() - started
    started = time.perf_counter()
    report = legacy_conversion.convert_world_to_quad(lines).to_dict()
    convert_s = time.perf_counter() - started
    quad = persistence.load_world_bytes(persistence.save_world_bytes(lines))
    quad_totals = surface_parity.totals(quad)
    if write_converted:
        (out_dir / f"{path.stem}-quad.mbworld").write_bytes(persistence.save_world_bytes(quad))
    audit, signatures = surface_parity.audit_world(quad, {})
    round_trip = surface_parity.state_hash(persistence.load_world_bytes(persistence.save_world_bytes(quad))) == surface_parity.state_hash(quad)
    if render:
        surface_parity._render(quad, out_dir / "renders", f"{path.stem}-quad", 0, (1100, 611))

    series, step_audits = [], []
    for step in range(steps):
        world_mod.step_world(quad, step_years)
        violations, signatures = surface_parity.audit_world(quad, signatures)
        step_audits.append(violations)
        series.append({"step": step + 1, "quad": surface_parity.totals(quad)})

    document = {
        "save": path.name,
        "elapsed_years": elapsed_years,
        "seconds": {"load": load_s, "convert": convert_s},
        "report": report,
        "totals_at_conversion": {"quad": quad_totals},
        "audit": audit,
        "round_trip_identical": round_trip,
        "steps": series,
        "step_audits": step_audits,
    }
    document["checks"] = _checks(report, audit, round_trip, step_audits)
    document["pass"] = all(c["pass"] for c in document["checks"].values())
    (out_dir / f"{path.stem}.json").write_text(json.dumps(surface_parity.sig(document), indent=2, sort_keys=True) + "\n")
    return document


def _verdict(passed: bool, failed: list[str], warned: list[str]) -> str:
    if not passed:
        return "FAIL: " + ", ".join(failed)
    return "pass" + (f" (warn: {', '.join(warned)})" if warned else "")


def _summary_md(documents: list[dict]) -> str:
    rows = [
        "| save | Myr | line nodes | quad cells | Hc Δ | Hm Δ | Hc over cap in save | Hc delaminated | worst untouched plate Hc Δ | stacked nodes "
        "| land line→quad | continental line→quad | uncovered | multiply covered | sea level after steps | result |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---|",
    ]
    for d in documents:
        r, c = d["report"], d["checks"]
        after = r["coverage"]["after"]
        fractions = r["area_fractions"]
        last = d["steps"][-1] if d["steps"] else None
        sea = f"{last['quad']['sea_level_m']:.0f}" if last else "-"
        failed = [name for name, check in c.items() if not check["pass"]]
        rows.append(
            f"| {d['save'].removeprefix('mantle-bloom-').removesuffix('.mbworld')} | {d['elapsed_years'] / 1e6:.1f} | {r['source']['nodes']} | {r['result']['cells']} "
            f"| {100 * r['extensive']['crustal_thickness_m']['delta_vs_capped']:+.2f}% "
            f"| {100 * r['extensive']['mantle_lithosphere_thickness_m']['delta_vs_capped']:+.2f}% "
            f"| {100 * r['source_over_cap']['crustal_thickness_m']:.2f}% "
            f"| {100 * max(r['clamped']['crustal_thickness_m'] - r['source_over_cap']['crustal_thickness_m'], 0.0):.2f}% "
            f"| {100 * c['plate_crust_volume_delta']['value']:.2f}% "
            f"| {100 * r['targets']['stacked_fraction']:.2f}% "
            f"| {fractions['land']['before']:.3f} → {fractions['land']['after']:.3f} "
            f"| {fractions['continental']['before']:.3f} → {fractions['continental']['after']:.3f} "
            f"| {100 * after['uncovered']:.2f}% "
            f"| {100 * after['multiply_covered']:.2f}% "
            f"| {sea} | {_verdict(d['pass'], failed, [name for name, check in c.items() if check.get('warn')])} |"
        )
    return "\n".join(rows) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("saves", nargs="*", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=3)
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--render", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--write-converted", action="store_true")
    parser.add_argument("--recheck", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    documents = []
    if args.recheck:
        for path in sorted(args.out.glob("*.json")):
            if path.name == "summary.json":
                continue
            document = json.loads(path.read_text())
            document["checks"] = _checks(document["report"], document["audit"], document["round_trip_identical"], document["step_audits"])
            document["pass"] = all(c["pass"] for c in document["checks"].values())
            path.write_text(json.dumps(surface_parity.sig(document), indent=2, sort_keys=True) + "\n")
            documents.append(document)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(convert_one, path, args.out, args.steps, args.step_years, args.render, args.write_converted) for path in args.saves]
            for future in futures:
                documents.append(future.result())
    for document in documents:
        print(f"{document['save']}: {'pass' if document['pass'] else 'FAIL'}", flush=True)
    documents.sort(key=lambda d: d["elapsed_years"])
    (args.out / "summary.json").write_text(
        json.dumps(surface_parity.sig({"tolerances": TOLERANCES, "saves": {d["save"]: {"pass": d["pass"], "checks": d["checks"]} for d in documents}}), indent=2, sort_keys=True) + "\n"
    )
    (args.out / "summary.md").write_text(_summary_md(documents))
    sys.exit(0 if all(d["pass"] for d in documents) else 1)


if __name__ == "__main__":
    main()

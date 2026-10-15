#!/usr/bin/env python3
"""RFC 0009 gate G1: does a craton mosaic give rifts a way around? (#426)

Today each continent's interior is one craton, so RFC 0010's minimum cut can't route a rift
around it (analysis/rfc0010). This probe reseeds generation as RFC 0009's mosaic, in-process
only: per plate, the columns eligible for seeding with the formation margin
(CRATON_FORMATION_MARGIN_KM) are split into cores around Poisson-disc nuclei
CORE_SPACING_KM apart, by nearest nucleus in graph hops, and columns within BELT_WIDTH_KM / 2
of a partition boundary become belts instead of craton. CRATON_FORMATION_MYR is set to
infinity, so the seeded layout is the only craton for the run. Belts are tracked through
every remap by a probe-only HISTORY surface field registered before the app loads; no saved
state changes, and the probe never saves.

Every cut that goes ahead (split or failed rift) is evaluated with
prototype_mincut_split's minimum cut at a 1,000 km corridor and floor 0.05, with RFC 0010's
admissibility gates and its strength: a belt-marked column scores RIFT_BELT_STRENGTH (0)
whatever its keel. Over the flow cuts whose band crosses craton it reports:
- L, line clearance: the share with an admissible minimum cut whose line touches no craton;
- B, band exposure: craton area in the FAILED_RIFT_BAND_MULT band, minimum cut where
  admissible and flow cut otherwise, over the same sum for the flow cuts.
Per 100 Myr checkpoint it reports belt area (belt-marked genuine continental columns
without craton) as a share of genuine continental area and of its t = 0 value, and craton
cores (connected craton components per plate).

`--layout blob` is the control: today's seeding (CRATON_SEED_MARGIN_KM), also with formation
off, so the two layouts differ only in where the seeded craton is.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py \
        --out /tmp/mosaic.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py \
        --layout blob --out /tmp/blob.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py \
        --belt-width-km 600 --out /tmp/mosaic-600.json     # RFC 0009's knob sweep
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py \
        --seeds 7 --years 1e8 --out /tmp/mosaic-smoke.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py \
        --years 0 --out /tmp/mosaic-t0.json     # generation only: cores and belt share
"""

from __future__ import annotations

from app import surface_fields

# Registered before any module that reads the field registry is imported, so every plate
# carries and remaps it. HISTORY keeps the oldest contributor's date, so a merged or coarsened
# cell stays a belt when any of its contributors was one.
BELT_FIELD = "probe_belt_formed_years"
BELT_NEVER = surface_fields.CRATON_UNFORMED_YEARS
surface_fields.SURFACE_FIELDS[BELT_FIELD] = surface_fields._field(
    float, BELT_NEVER, surface_fields.RemapClass.HISTORY, sentinel=BELT_NEVER
)

import argparse  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import compare_cratons  # noqa: E402
import numpy as np  # noqa: E402
import prototype_mincut_split as mincut  # noqa: E402
from app import cratons, merge_split, quad_tectonics  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import PLANET_RADIUS_KM, line_spacing_rad  # noqa: E402
from scipy.sparse.csgraph import connected_components, dijkstra  # noqa: E402

# RFC 0009's starting knobs; --core-spacing-km and --belt-width-km sweep them.
CORE_SPACING_KM = 1500.0
BELT_WIDTH_KM = 300.0
BELT_SEED_AGE_YEARS = 1.8e9
# The gate's corridor and RFC 0010's default floor.
CORRIDOR_KM = 1000.0
FLOOR = 0.05
TAG = f"w{CORRIDOR_KM:g}_f{FLOOR:g}"
# RFC 0010's strength for a belt-marked column: the cheapest path.
RIFT_BELT_STRENGTH = 0.0
# Seeds the nuclei per (world seed, plate); separate from the world's own streams.
_NUCLEUS_SEED_TAG = 9009


def _nuclei(points: np.ndarray, rng: np.random.Generator, min_chord: float) -> list[int]:
    """Greedy Poisson-disc sample: visit `points` in random order, keeping each one at least
    `min_chord` from every point kept so far."""
    kept: list[int] = []
    for idx in rng.permutation(len(points)):
        if not kept or np.min(np.linalg.norm(points[kept] - points[idx], axis=1)) >= min_chord:
            kept.append(int(idx))
    return kept


def _partition(eligible, adjacency, points, rng, min_chord) -> tuple[np.ndarray, list[int]]:
    """Per cell, the id of its core (nearest nucleus by hops within its connected piece of
    `eligible`), -1 off it; and the nucleus count of each piece."""
    members = np.flatnonzero(eligible)
    sub = adjacency[members][:, members]
    pieces, piece = connected_components(sub, directed=False)
    labels = np.full(len(eligible), -1)
    counts = []
    next_label = 0
    for p in range(pieces):
        local = np.flatnonzero(piece == p)
        nuclei = local[_nuclei(points[members[local]], rng, min_chord)]
        counts.append(len(nuclei))
        _, _, source = dijkstra(sub, indices=nuclei, unweighted=True, min_only=True, return_predecessors=True)
        lookup = {int(n): next_label + k for k, n in enumerate(nuclei)}
        labels[members[local]] = [lookup[int(s)] for s in source[local]]
        next_label += len(nuclei)
    return labels, counts


def _belts(plate, labels, adjacency, half_hops: int) -> np.ndarray:
    """Cells of the partition within `half_hops` - 1 hops of a cell bordering another core,
    so a belt is about 2 * half_hops cells wide across a boundary."""
    graph = adjacency.tocoo()
    i, j = graph.row, graph.col
    border = (labels[i] >= 0) & (labels[j] >= 0) & (labels[i] != labels[j])
    boundary = np.zeros(len(labels), dtype=bool)
    boundary[i[border]] = True
    width = max(half_hops - 1, 0)
    return (labels >= 0) & (quad_tectonics.hop_distance(plate, boundary, width, adjacency) <= width)


SEEDING: list[dict] = []


def seed_mosaic(world) -> float:
    """RFC 0009 solution 1 in place of cratons.seed_initial_cratons; returns seeded m^3."""
    spacing_rad = line_spacing_rad(world.node_density)
    margin = cratons._hops(cratons.CRATON_FORMATION_MARGIN_KM, spacing_rad)
    half_hops = cratons._hops(0.5 * BELT_WIDTH_KM, spacing_rad)
    min_chord = 2.0 * np.sin(0.5 * CORE_SPACING_KM / PLANET_RADIUS_KM)
    seeded = 0.0
    for plate in world.plates:
        host = cratons._genuine_continental(plate)
        interior = cratons._interior_hops(plate, host, margin)
        if interior is None:
            continue
        hc = plate.collect("crustal_thickness_m")
        craton = plate.collect("craton_crust_m")
        eligible = host & (hc <= cratons.CRATON_HOST_MAX_HC_M) & (interior >= margin) & (craton <= 0.0)
        if not np.any(eligible):
            continue
        adjacency = quad_tectonics._adjacency_matrix(plate)
        points = plate.all_points_and_elevation()[0]
        rng = np.random.default_rng((world.seed, plate.plate_id, _NUCLEUS_SEED_TAG))
        labels, counts = _partition(eligible, adjacency, points, rng, min_chord)
        belt = _belts(plate, labels, adjacency, half_hops)
        core = eligible & ~belt
        areas = cratons._areas(world, plate)
        craton[core] = hc[core]
        formed = plate.collect("craton_formed_years")
        formed[core] = world.elapsed_years - cratons.CRATON_SEED_AGE_YEARS
        stable = plate.collect("stable_continental_myr")
        stable[core] = np.maximum(stable[core], cratons.CRATON_SEED_AGE_YEARS / 1e6)
        marks = plate.collect(BELT_FIELD)
        marks[belt] = world.elapsed_years - BELT_SEED_AGE_YEARS
        plate.set_fields_on_plate(craton_crust_m=craton, craton_formed_years=formed, stable_continental_myr=stable, **{BELT_FIELD: marks})
        seeded += float(np.dot(hc[core], areas[core]))
        SEEDING.append({
            "plate": plate.plate_id, "interior_mkm2": float(areas[eligible].sum() / 1e12),
            "belt_mkm2": float(areas[belt].sum() / 1e12), "nuclei_per_piece": counts,
            "cores_per_piece": _cores_per_piece(plate, eligible, core, adjacency),
        })
    return seeded


def _cores_per_piece(plate, eligible, core, adjacency) -> list[int]:
    """Connected craton components inside each connected piece of the seeded interior."""
    members = np.flatnonzero(eligible)
    _, piece = connected_components(adjacency[members][:, members], directed=False)
    core_members = np.flatnonzero(core)
    _, comp = connected_components(adjacency[core_members][:, core_members], directed=False)
    piece_of = dict(zip(members.tolist(), piece.tolist()))
    per_piece: dict[int, set[int]] = {p: set() for p in set(piece.tolist())}
    for cell, c in zip(core_members.tolist(), comp.tolist()):
        per_piece[piece_of[cell]].add(c)
    return [len(per_piece[p]) for p in sorted(per_piece)]


_KEEL_STRENGTH = mincut._strength


def _rift_strength(plate) -> np.ndarray:
    """RFC 0010's strength: craton or keel, but RIFT_BELT_STRENGTH on a belt."""
    marked = plate.collect(BELT_FIELD) < 0.5 * BELT_NEVER
    return np.where(marked, RIFT_BELT_STRENGTH, _KEEL_STRENGTH(plate))


def _layout_snapshot(world) -> dict:
    """Belt and craton areas, and craton cores, across every plate."""
    spacing = line_spacing_rad(world.node_density)
    out = {"myr": world.elapsed_years / 1e6, "genuine_mkm2": 0.0, "craton_mkm2": 0.0, "belt_mkm2": 0.0,
           "belt_marked_mkm2": 0.0, "belt_on_craton_mkm2": 0.0, "cores": 0, "cores_over_1mkm2": 0}
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        areas = plate.accounting_areas_m2(spacing)
        genuine = cratons._genuine_continental(plate)
        craton = plate.collect("craton_crust_m") > 0.0
        marked = plate.collect(BELT_FIELD) < 0.5 * BELT_NEVER
        out["genuine_mkm2"] += float(areas[genuine].sum() / 1e12)
        out["craton_mkm2"] += float(areas[craton].sum() / 1e12)
        out["belt_mkm2"] += float(areas[marked & genuine & ~craton].sum() / 1e12)
        out["belt_marked_mkm2"] += float(areas[marked].sum() / 1e12)
        out["belt_on_craton_mkm2"] += float(areas[marked & craton].sum() / 1e12)
        members = np.flatnonzero(craton)
        if len(members):
            sub = quad_tectonics._adjacency_matrix(plate)[members][:, members]
            _, comp = connected_components(sub, directed=False)
            comp_area = np.bincount(comp, weights=areas[members])
            out["cores"] += len(comp_area)
            out["cores_over_1mkm2"] += int(np.count_nonzero(comp_area >= 1e12))
    out["belt_share_of_genuine"] = out["belt_mkm2"] / out["genuine_mkm2"] if out["genuine_mkm2"] else 0.0
    return out


def _gate(rows: list[dict]) -> dict:
    """L and B over the cuts whose flow band crosses craton."""
    through = [r for r in rows if r["flow"]["through_craton"] and "skipped" not in r["variants"][TAG]]
    flow_band = sum(r["flow"]["band_craton_area_mkm2"] for r in through)
    used_band = sum(
        r["variants"][TAG]["band_craton_area_mkm2"] if r["variants"][TAG]["admissible"] else r["flow"]["band_craton_area_mkm2"]
        for r in through
    )
    clear = sum(r["variants"][TAG]["admissible"] and not r["variants"][TAG]["line_through_craton"] for r in through)
    return {
        "cuts": len(rows), "through": len(through), "line_clear_admissible": clear,
        "L": clear / len(through) if through else None,
        "flow_band_craton_mkm2": flow_band, "used_band_craton_mkm2": used_band,
        "B": used_band / flow_band if flow_band > 0.0 else None,
        "admissible": sum(r["variants"][TAG]["admissible"] for r in through),
    }


def _compact(row: dict) -> dict:
    v = row["variants"][TAG]
    keep = ("admissible", "line_through_craton", "through_craton", "band_craton_area_mkm2", "residual_ratio",
            "size_ok", "poles_ok", "residual_ok", "connected") if "skipped" not in v else ("skipped",)
    return {"kind": row["kind"], "myr": row["myr"], "plate": row["plate"], "nodes": row["nodes"],
            "plate_area_mkm2": row["plate_area_mkm2"], "plate_craton_area_mkm2": row["plate_craton_area_mkm2"],
            "flow": {k: row["flow"][k] for k in ("through_craton", "line_through_craton", "band_area_mkm2", "band_craton_area_mkm2")},
            "variants": {TAG: {k: v[k] for k in keep}}}


def _run(seed: int, layout: str, years: float, every: float, step: float) -> dict:
    for name, value in compare_cratons._DEFAULTS.items():
        setattr(cratons, name, value)
    cratons.CRATON_FORMATION_MYR = float("inf")
    mincut.CORRIDORS_KM, mincut.FLOORS = (CORRIDOR_KM,), (FLOOR,)
    mincut._strength = _rift_strength
    original_seed, original_cut = cratons.seed_initial_cratons, cratons.record_rift_cut
    rows: list[dict] = []

    def tapped(world, plate, points, cut_normal, band_sin, kind):
        rows.append(_compact(mincut._evaluate(world, plate, np.asarray(points, dtype=float), np.asarray(cut_normal, dtype=float), kind)))
        return original_cut(world, plate, points, cut_normal, band_sin, kind)

    SEEDING.clear()
    if layout == "mosaic":
        cratons.seed_initial_cratons = seed_mosaic
    cratons.record_rift_cut = tapped
    started = time.perf_counter()
    try:
        world = world_mod.generate_world(seed)
        cratons.ensure_initialized(world)
        generation_s = time.perf_counter() - started
        checkpoints = [_layout_snapshot(world)]
        next_checkpoint = every
        while world.elapsed_years < years - 1.0:
            world_mod.step_world(world, min(step, years - world.elapsed_years))
            if world.elapsed_years >= next_checkpoint - 1.0:
                checkpoints.append({**_layout_snapshot(world), "gate_so_far": _gate(rows), "wall_s": time.perf_counter() - started})
                print(json.dumps({"seed": seed, "layout": layout, **{k: checkpoints[-1][k] for k in ("myr", "belt_mkm2", "craton_mkm2", "cores")},
                                  **{k: checkpoints[-1]["gate_so_far"][k] for k in ("through", "L", "B")}}, default=float), flush=True)
                next_checkpoint += every
    finally:
        cratons.seed_initial_cratons, cratons.record_rift_cut = original_seed, original_cut
    belt0 = checkpoints[0]["belt_mkm2"]
    for c in checkpoints:
        c["belt_vs_t0"] = c["belt_mkm2"] / belt0 if belt0 else None
    return {"seed": seed, "layout": layout, "generation_s": generation_s, "wall_s": time.perf_counter() - started,
            "seeding": list(SEEDING), "gate": _gate(rows), "checkpoints": checkpoints, "cuts": rows}


def _pooled(runs: list[dict]) -> dict:
    through = sum(r["gate"]["through"] for r in runs)
    clear = sum(r["gate"]["line_clear_admissible"] for r in runs)
    flow = sum(r["gate"]["flow_band_craton_mkm2"] for r in runs)
    used = sum(r["gate"]["used_band_craton_mkm2"] for r in runs)
    return {"through": through, "L": clear / through if through else None, "B": used / flow if flow else None}


def main() -> None:
    global CORE_SPACING_KM, BELT_WIDTH_KM
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--layout", choices=("mosaic", "blob"), default="mosaic")
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--every", type=float, default=100e6)
    parser.add_argument("--step", type=float, default=5e6)
    parser.add_argument("--core-spacing-km", type=float, default=CORE_SPACING_KM)
    parser.add_argument("--belt-width-km", type=float, default=BELT_WIDTH_KM)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    CORE_SPACING_KM, BELT_WIDTH_KM = args.core_spacing_km, args.belt_width_km
    app_root = Path(merge_split.__file__).resolve().parents[2]
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=app_root, text=True).strip()
    result = {"model_commit": commit, "layout": args.layout, "years": args.years, "step": args.step,
              "corridor_km": CORRIDOR_KM, "floor": FLOOR, "core_spacing_km": CORE_SPACING_KM,
              "belt_width_km": BELT_WIDTH_KM, "runs": []}
    for seed in map(int, args.seeds.split(",")):
        result["runs"].append(_run(seed, args.layout, args.years, args.every, args.step))
        result["pooled"] = _pooled(result["runs"])
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=1, default=float) + "\n")
        print(json.dumps({"seed": seed, "gate": result["runs"][-1]["gate"], "pooled": result["pooled"]}, default=float), flush=True)


if __name__ == "__main__":
    main()

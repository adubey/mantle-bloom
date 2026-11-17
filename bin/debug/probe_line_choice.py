"""Probe of RFC 0018's rift-line choice (#888): does picking the line over the whole plate, by pull over strength,
follow inherited weakness better than the peak-pull corridor?

Usage (repo root, this worktree's venv; one heavy run, follow the compute-budget skill):
    backend/.venv/bin/python bin/debug/probe_line_choice.py --seed 106533960 --out <file.json> \
        [--setting '{"k_cap": 0.2}'] [--at 50 100 150]

Steps a generated development world plates-only with the rift-line prototype (as ``run_rift_damage_d1.py``) and, at
each ``--at`` time, compares three line choices on every plate without straining anything:

* ``corridor``: the prototype's rule, the cheapest path inside the strip around the peak-pull section;
* ``whole``: the cheapest path inside each of a set of overlapping strips that together cover the plate, for each of the
  six directions, keeping the one with the largest section load ``Λ = P / S``;
* ``whole_deep``: as ``whole``, with inherited-weakness cells a further ``--deepen`` weaker (a sensitivity arm).

It also records strength by rock class and how much of each plate's weak continental rock the corridors contain.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import rift_damage_mechanism as mech  # noqa: E402
import rift_damage_step as step_mod  # noqa: E402
import rift_line as rl  # noqa: E402
from run_rift_damage_d1 import base_parser  # noqa: E402
from app import thin_shell  # noqa: E402

WINDOW = 0.15  # area share of the plate inside one strip
STRIDE = 0.075


def windows(sections: rl.Sections, p) -> list:
    """Overlapping strips across the plate for every direction: (direction, cell mask, pull at the strip's middle)."""
    out = []
    for d in range(sections.pull.shape[1]):
        share = sections.far_share[:, d]
        centre = p.a_min + WINDOW / 2.0
        while centre <= 1.0 - p.a_min - WINDOW / 2.0 + 1e-9:
            cells = np.abs(share - centre) <= WINDOW / 2.0
            if cells.sum() >= 4:
                middle = int(np.argmin(np.abs(share - centre)))
                out.append((d, cells, float(sections.pull[middle, d])))
            centre += STRIDE
    return out


def evaluate(prep, sections, chain, direction, n_eff, p, frame) -> dict | None:
    if chain is None or len(chain) < 2:
        return None
    geometry = rl.path_geometry(prep.points, chain, *frame)
    tensor = prep.tensor + p.n_plume * prep.heated[:, None] * np.array([1.0, 1.0, 0.0])
    load = np.maximum(rl.path_normal_load(tensor[chain], geometry.normal), 0.0)
    pull = float(np.sum(load * geometry.ds))
    strength = float(np.sum((n_eff[chain] * geometry.ds)[~prep.broken[chain]]))
    lam, is_open = rl.section_load(pull, strength, p)
    ds, total = geometry.ds, max(float(geometry.ds.sum()), 1e-9)
    share = lambda flag: float(ds[flag[chain]].sum() / total)  # noqa: E731
    continental = prep.continental
    cont_len = float(ds[continental[chain]].sum())
    mask, why = rl.split_masks(prep.graph, prep.areas, chain, p, points=prep.points)
    return {"direction": direction, "lam": lam, "open": bool(is_open), "key": rl.rank_key(pull, lam, is_open), "pull_n": pull, "strength_n": strength,
            "length_km": total / 1e3, "oceanic_share": share(~continental), "weak_share": share(prep.weak),
            "weak_share_of_continental_part": float(ds[(prep.weak & continental)[chain]].sum() / cont_len) if cont_len > 0 else None,
            "valid": mask is not None, "reject": why,
            "smaller_share": None if mask is None else float(min(prep.areas[mask].sum(), prep.areas[~mask].sum()) / prep.areas.sum())}


def best_line(prep, sections, strips, n_eff, p, frame) -> dict | None:
    cost = np.where(prep.broken, 0.0, n_eff)
    best = None
    for direction, cells, _ in strips:
        chain = rl.find_line(prep.graph, cells, prep.boundary, sections.strike[:, direction], cost)
        found = evaluate(prep, sections, chain, direction, n_eff, p, frame)
        if found is not None and (best is None or found["key"] > best["key"]):
            best = found
    if best is not None:
        best.pop("key")
    return best


def strength_classes(prep, plate) -> dict:
    craton = plate.collect("craton_crust_m") > 0.0
    classes = {"oceanic_ridge": prep.kind == 1, "oceanic_mature": prep.kind == 2, "continental_weak": prep.continental & prep.weak & ~craton,
               "continental_plain": prep.continental & ~prep.weak & ~craton, "craton": prep.continental & craton}
    return {name: {"area_m2": float(prep.areas[m].sum()), "n_eff_sum": float(np.dot(prep.n_eff[m], prep.areas[m]))} for name, m in classes.items()}


def probe_plate(world, plate, p, deepen: float) -> dict | None:
    prep = step_mod.prepare_plate(world, plate, p, 1e6)
    if len(prep.points) < 40:
        return None
    field_, _ = step_mod.force_inputs(world, plate)
    sections = rl.section_profile(prep.points, prep.areas, field_.total_n, p.k_cut)
    frame = thin_shell._tangent_basis(prep.points)
    corridors = rl.candidate_corridors(sections, p)
    weak_cont = prep.weak & prep.continental
    in_corridor = np.zeros(len(prep.points), bool)
    for corridor in corridors:
        in_corridor |= corridor.cells
    started = time.perf_counter()
    corridor_line = best_line(prep, sections, [(c.direction, c.cells, c.peak_pull) for c in corridors], prep.n_eff, p, frame)
    t_corridor = time.perf_counter() - started
    strips = windows(sections, p)
    started = time.perf_counter()
    whole = best_line(prep, sections, strips, prep.n_eff, p, frame)
    t_whole = time.perf_counter() - started
    deep = best_line(prep, sections, strips, prep.n_eff * np.where(prep.weak, deepen, 1.0), p, frame)
    area = float(prep.areas.sum())
    return {"plate": int(plate.plate_id), "cells": len(prep.points), "continental_area_share": float(prep.areas[prep.continental].sum() / area),
            "weak_continental_area_share": float(prep.areas[weak_cont].sum() / area), "corridor_area_share": float(prep.areas[in_corridor].sum() / area),
            "weak_continental_inside_corridor": float(prep.areas[weak_cont & in_corridor].sum() / max(prep.areas[weak_cont].sum(), 1e-9)) if weak_cont.any() else None,
            "strips": len(strips), "seconds": {"corridor": t_corridor, "whole": t_whole}, "classes": strength_classes(prep, plate),
            "corridor": corridor_line, "whole": whole, "whole_deep": deep}


def main() -> None:
    parser = base_parser(__doc__)
    parser.add_argument("--at", type=int, nargs="+", default=[50, 100, 150], help="Myr")
    parser.add_argument("--deepen", type=float, default=0.5, help="extra strength factor on inherited-weakness cells in the sensitivity arm")
    args = parser.parse_args()
    params = mech.Params(**args.setting)
    step_mod.install(params)
    from app import world as world_module

    world = world_module.generate_world(args.seed)
    world.simulate_climate_biomes = False
    records = []
    for myr in range(1, max(args.at) + 1):
        world_module.step_world(world, 1e6)
        if myr in args.at:
            plates = [r for r in (probe_plate(world, plate, params, args.deepen) for plate in world.plates) if r is not None]
            records.append({"myr": myr, "plates": plates})
            print(json.dumps({"myr": myr, "plates": len(plates)}), flush=True)
    args.out.write_text(json.dumps({"seed": args.seed, "setting": args.setting, "window": WINDOW, "stride": STRIDE, "deepen": args.deepen,
                                    "records": records}, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()

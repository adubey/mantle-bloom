"""Does the flow's own coherence at a contact decide whether plates meet, subduct and move together? (#998, exploration)

The proposal under test: adjacent plates are usually pulled the same way by the mantle flow; where the flow converges on a
contact the plates are pushed together and subduct. For every plate-pair contact at each sample time this records
  * `flow`: the relative motion of the two plates under the flow alone, i.e. each plate's own flow-only angular velocity
    (`basal_drag_coefficients`: K omega = b, the 'mantle_flow_only' counterfactual of `plate_torque_breakdown`), evaluated
    at the contact nodes: length-weighted mean closing rate (cm/yr, + = converging), relative speed, and the opposition
    |v_rel| / (|v_a| + |v_b|) with v_a, v_b the two plates' own node velocities (about 1 = opposite, 0 = alike);
  * `actual`: the same three numbers under the solved motion (`solve_plates_quasistatic`);
  * `subducting_share`: the length share of the contact's nodes that are subducting (either side);
  * `D`: the mean downwelling potential at the contact, against the ORIGINAL centres at `--reference-scale`.
Local node velocities are used (not plate-mean speeds), so opposition here is not the number in `contact_asymmetry.py`.

    run        backend/.venv/bin/python bin/debug/contact_flow_coherence.py run --settings analysis/issue998/settings.json \\
                   --setting b123-brake-0.5-for-flow4 --seed 591453549 --at 0 5 20 50 --flow-scale 4 --out <captures>/cfc/591453549.json
    queue      same, over `--seeds-file`, one process per seed.
    aggregate  backend/.venv/bin/python bin/debug/contact_flow_coherence.py aggregate --runs <captures>/cfc --out analysis/issue998/contact-flow-coherence.json

Exploration only; no `backend/app` change. About 1 GB per seed (compute-budget first).
"""

import argparse
import copy
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "backend")]

MIN_CONTACT_KM = 500.0


def contact_rows(world, centers) -> list[dict]:
    from app import lithosphere, torque
    from app.elevation_lines import line_spacing_rad
    from downwelling_check import downwelling_potential

    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
    cache = world.boundary_search_cache
    if cache is None:
        cache = torque.BoundarySearchCache()
        world.boundary_search_cache = cache
    forms = torque.force_forms(world)
    solved = torque.solve_plates_quasistatic(world, diagnostics={})
    flow_only = {}
    for plate in plates:
        b, k = torque.basal_drag_coefficients(plate, world, spacing)
        flow_only[plate.plate_id] = np.linalg.lstsq(k, b, rcond=None)[0]
    scale = lithosphere.PLANET_RADIUS_M * 100.0
    acc = defaultdict(lambda: defaultdict(float))
    for plate in plates:
        neighbours = cache.neighbours(plate, plates, reach)
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
        torque.apply_force_forms(world, plate, neighbours, inputs, forms)
        band = (inputs.dist_to_neighbor <= reach) & (inputs.neighbor_plate_id >= 0)
        if not band.any():
            continue
        subducting = torque.subducting_boundary_mask(plate, inputs, reach)
        km = torque.boundary_node_length_m(inputs, spacing) / 2e3
        depth = downwelling_potential(inputs.own_points, centers)
        for q in np.unique(inputs.neighbor_plate_id[band]):
            m = band & (inputs.neighbor_plate_id == q)
            other = plates[[p.plate_id for p in plates].index(int(q))]
            key = tuple(sorted((plate.plate_id, int(q))))
            a = acc[key]
            w = km[m]
            for tag, omegas in (("actual", solved), ("flow", flow_only)):
                va = np.cross(omegas[plate.plate_id], inputs.own_points[m]) * scale
                vb = np.cross(omegas[other.plate_id], inputs.own_points[m]) * scale
                rel = va - vb
                a[f"{tag}_closing"] += float((w * np.sum(rel * inputs.direction_to_neighbor[m], axis=-1)).sum())
                a[f"{tag}_rel"] += float((w * np.linalg.norm(rel, axis=-1)).sum())
                a[f"{tag}_sum"] += float((w * (np.linalg.norm(va, axis=-1) + np.linalg.norm(vb, axis=-1))).sum())
            a["km"] += float(w.sum())
            a["sub_km"] += float(w[subducting[m]].sum())
            a["D"] += float((w * depth[m]).sum())
    rows = []
    for (pa, pb), a in sorted(acc.items()):
        if a["km"] < MIN_CONTACT_KM:
            continue
        row = {"a": pa, "b": pb, "km": a["km"], "subducting_share": a["sub_km"] / a["km"], "D": a["D"] / a["km"]}
        for tag in ("actual", "flow"):
            row[tag] = {"closing": a[f"{tag}_closing"] / a["km"], "rel": a[f"{tag}_rel"] / a["km"],
                        "opposition": a[f"{tag}_rel"] / a[f"{tag}_sum"] if a[f"{tag}_sum"] > 0 else None}
        rows.append(row)
    return rows


def run_one(args: argparse.Namespace) -> None:
    import calibrate_plate_forces as cal
    import run_rift_lines_d1 as d1
    from app import mantle

    setting = cal.load_setting(args.settings, args.setting)
    world = d1.make_world(SimpleNamespace(seed=args.seed, node_density=1.0, world_knobs=setting["world_knobs"],
                                          split_method=setting["split_method"], params=setting["params"]))
    cal.apply_module_knobs(setting["module_knobs"])
    centers = copy.deepcopy(world.mantle_centers)
    reference = mantle.MANTLE_FLOW_REFERENCE_RATE * args.reference_scale
    for c in centers:
        c.strength *= args.reference_scale / reference
    for c in world.mantle_centers:
        c.strength *= args.flow_scale
    samples = []

    def take(myr: float) -> None:
        samples.append({"myr": myr, "contacts": contact_rows(world, centers)})
        print(json.dumps({"myr": myr, "contacts": len(samples[-1]["contacts"])}), flush=True)

    wanted = sorted(set(args.at))
    if wanted[0] == 0:
        take(0.0)
    scratch = args.out.parent / f"{args.out.stem}-scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    stepping = SimpleNamespace(years=max(wanted), step_myr=args.step_myr, max_steps=None, snapshot_every=0.0)
    recorder = None
    if setting["split_method"] == "rift_lines":
        import rift_line_recorder as rec
        from app import rift_damage

        recorder = rec.Recorder(rift_damage.runtime(world).params.eps_cand)
    if max(wanted) > 0:
        d1.step_all(world, stepping, scratch, recorder, lambda w, myr: take(myr) if round(myr, 6) in wanted else None)
    args.out.write_text(json.dumps({"seed": args.seed, "setting": args.setting, "flow_scale": args.flow_scale, "samples": samples}) + "\n")


def queue(args: argparse.Namespace) -> None:
    seeds = json.loads(args.seeds_file.read_text())["seeds"]
    args.out.mkdir(parents=True, exist_ok=True)
    for seed in seeds:
        target = args.out / f"{seed}.json"
        if target.exists():
            continue
        command = [sys.executable, __file__, "run", "--settings", str(args.settings), "--setting", args.setting, "--seed", str(seed), "--step-myr", str(args.step_myr),
                   "--flow-scale", str(args.flow_scale), "--reference-scale", str(args.reference_scale), "--out", str(target), "--at", *map(str, args.at)]
        with open(args.out / f"{seed}.log", "w") as log:
            code = subprocess.call(command, stdout=log, stderr=subprocess.STDOUT)
        print(seed, "exit", code, time.strftime("%H:%M"), flush=True)


def weighted(rows: list[dict], key) -> float | None:
    vals = [(key(r), r["km"]) for r in rows if key(r) is not None]
    return float(sum(v * w for v, w in vals) / sum(w for _, w in vals)) if vals else None


def aggregate(args: argparse.Namespace) -> None:
    runs = [json.loads(p.read_text()) for p in sorted(args.runs.glob("*.json"))]
    result = {"seeds": len(runs), "by_time": {}}
    for t in sorted({s["myr"] for r in runs for s in r["samples"]}):
        rows = [c for r in runs for s in r["samples"] if s["myr"] == t for c in s["contacts"]]
        # Bins by the flow alone: converging (flow closing > +0.5 cm/yr), opening (< -0.5), neither (flow-only sliding or alike).
        bins = {"flow_converging": [r for r in rows if r["flow"]["closing"] > 0.5], "flow_opening": [r for r in rows if r["flow"]["closing"] < -0.5],
                "flow_neither": [r for r in rows if abs(r["flow"]["closing"]) <= 0.5]}
        # Bins by how alike the flow moves the two plates (opposition under the flow alone).
        alike = {"flow_alike_opposition_lt_0.3": [r for r in rows if r["flow"]["opposition"] is not None and r["flow"]["opposition"] < 0.3],
                 "flow_opposed_opposition_gt_0.7": [r for r in rows if r["flow"]["opposition"] is not None and r["flow"]["opposition"] > 0.7]}
        result["by_time"][str(t)] = {"contacts": len(rows), "length_km": float(sum(r["km"] for r in rows)),
                                     "flow_opposition_all": weighted(rows, lambda r: r["flow"]["opposition"]),
                                     "actual_opposition_all": weighted(rows, lambda r: r["actual"]["opposition"]),
                                     "bins": {name: {"contacts": len(b), "length_share": float(sum(r["km"] for r in b) / sum(r["km"] for r in rows)) if b else 0.0,
                                                     "subducting_share": weighted(b, lambda r: r["subducting_share"]),
                                                     "actual_opposition": weighted(b, lambda r: r["actual"]["opposition"]),
                                                     "actual_rel_speed_cm_per_yr": weighted(b, lambda r: r["actual"]["rel"]),
                                                     "mean_D": weighted(b, lambda r: r["D"])} for name, b in {**bins, **alike}.items()}}
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    for t, b in result["by_time"].items():
        print(f"t={t}: contacts {b['contacts']}  flow-only opposition {b['flow_opposition_all']}  actual {b['actual_opposition_all']}")
        for name, x in b["bins"].items():
            print(f"   {name:34s} n={x['contacts']:4d} share {x['length_share']:.2f}  subducting {x['subducting_share']}  actual opposition {x['actual_opposition']}  rel speed {x['actual_rel_speed_cm_per_yr']}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("run", run_one), ("queue", queue)):
        p = sub.add_parser(name)
        p.add_argument("--settings", type=Path, required=True)
        p.add_argument("--setting", required=True)
        p.add_argument("--at", type=float, nargs="+", default=[0, 5, 20, 50])
        p.add_argument("--step-myr", type=float, default=1.0)
        p.add_argument("--flow-scale", type=float, default=1.0)
        p.add_argument("--reference-scale", type=float, default=4.0)
        p.add_argument("--out", type=Path, required=True)
        if name == "run":
            p.add_argument("--seed", type=int, required=True)
        else:
            p.add_argument("--seeds-file", type=Path, default=Path("analysis/issue958/development-seeds.json"))
        p.set_defaults(fn=fn)
    a = sub.add_parser("aggregate")
    a.add_argument("--runs", type=Path, required=True)
    a.add_argument("--out", type=Path, required=True)
    a.set_defaults(fn=aggregate)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()

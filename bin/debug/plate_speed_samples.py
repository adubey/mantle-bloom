"""Sample the model's plate relative speeds per contact, boundary class and plate-area pair (#998, step 2).

One world, plates-only: the initial plate speeds come from the torque solve (`solve_plates_quasistatic`), recorded at
t=0 (before any step) and at the later times in `--at`. At each time it writes
  * `contacts`: every boundary node in the contact band as (length km, relative speed cm/yr, closing rate cm/yr),
    grouped by plate pair (each contact is seen from both plates; each side counts half), and the plates' areas;
  * `plates`: `plate_torque_breakdown.sample`'s per-plate rows (torque terms, slab pull, ridge push, basal drag,
    contact drag, consumed floor), for the attribution in step 4.
Classes are the Earth reference's (`earth_speed_reference.py`): sliding if |closing| < 0.5 x speed, else converging
or diverging by the sign of the closing rate; scoring is in `score_speed_samples.py`.

    run    one seed (about 1 GB, plates-only; compute-budget first):
           backend/.venv/bin/python bin/debug/plate_speed_samples.py run --settings analysis/issue975/settings.json \\
               --setting rifting-baseline --seed 591453549 --at 0 1 2 5 10 --step-myr 1 --out <captures>/591453549.json
           add `--flow-scale 2` to multiply the static mantle-flow amplitude (the #998 flow test).
    queue  seeds in series, one process each, skipping finished ones:
           backend/.venv/bin/python bin/debug/plate_speed_samples.py queue --settings ... --setting ... \\
               --seeds-file analysis/issue958/development-seeds.json --at 0 1 2 5 10 --out <captures>/<setting>
"""

import argparse
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


def record_contacts(world, omegas) -> dict:
    from app import lithosphere, torque
    from app.elevation_lines import line_spacing_rad

    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
    forms = torque.force_forms(world)
    cache = world.boundary_search_cache
    scale = lithosphere.PLANET_RADIUS_M * 100.0
    nodes = defaultdict(list)  # (a, b) -> [(km, speed, closing), ...] from either side, half weight each
    for plate in plates:
        neighbours = cache.neighbours(plate, plates, reach)
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
        torque.apply_force_forms(world, plate, neighbours, inputs, forms)
        band = (inputs.dist_to_neighbor <= reach) & (inputs.neighbor_plate_id >= 0)
        if not band.any():
            continue
        length = torque.boundary_node_length_m(inputs, spacing)
        neighbour_omega = np.zeros_like(inputs.neighbor_omega)
        for other in plates:
            neighbour_omega[inputs.neighbor_plate_id == other.plate_id] = omegas[other.plate_id]
        v_rel = (np.cross(omegas[plate.plate_id], inputs.own_points) - np.cross(neighbour_omega, inputs.own_points)) * scale
        closing = np.sum(v_rel * inputs.direction_to_neighbor, axis=-1)
        speed = np.linalg.norm(v_rel, axis=-1)
        for q in np.unique(inputs.neighbor_plate_id[band]):
            m = band & (inputs.neighbor_plate_id == q)
            key = tuple(sorted((plate.plate_id, int(q))))
            nodes[key] += [(float(a) / 2e3, round(float(b), 3), round(float(c), 3)) for a, b, c in zip(length[m], speed[m], closing[m])]
    areas = {str(p.plate_id): float(p.node_areas_m2().sum() / 1e6) for p in plates}
    return {"areas_km2": areas, "contacts": [{"a": a, "b": b, "nodes": v} for (a, b), v in sorted(nodes.items())]}


def run_one(args: argparse.Namespace) -> None:
    import calibrate_plate_forces as cal
    import plate_torque_breakdown as breakdown
    import run_rift_lines_d1 as d1
    from app import torque

    if args.seed in cal.refused_seeds():
        sys.exit(f"seed {args.seed} is held out")
    setting = cal.load_setting(args.settings, args.setting)
    world = d1.make_world(SimpleNamespace(seed=args.seed, node_density=1.0, world_knobs=setting["world_knobs"],
                                          split_method=setting["split_method"], params=setting["params"]))
    cal.apply_module_knobs(setting["module_knobs"])
    # The static mantle-flow field is drawn at generation with a fixed reference rate; scale it after generation so
    # every seed starts from the same world and only the stepped flow torque changes (exploration of #998 only).
    for center in world.mantle_centers:
        center.strength *= args.flow_scale
    started = time.perf_counter()
    samples = []

    def take(myr: float) -> None:
        omegas = torque.solve_plates_quasistatic(world, diagnostics={})
        samples.append({"myr": myr, **record_contacts(world, omegas), "breakdown": breakdown.sample(world)})
        print(json.dumps({"myr": myr, "plates": len(world.plates), "wall_s": round(time.perf_counter() - started, 1)}), flush=True)

    wanted = sorted(set(args.at))
    if wanted and wanted[0] == 0:
        take(0.0)
    final = max(wanted)
    scratch = args.out.parent / f"{args.out.stem}-scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    stepping = SimpleNamespace(years=final, step_myr=args.step_myr, max_steps=None, snapshot_every=0.0)
    recorder = None
    if setting["split_method"] == "rift_lines":
        import rift_line_recorder as rec
        from app import rift_damage

        recorder = rec.Recorder(rift_damage.runtime(world).params.eps_cand)
    if final > 0:
        d1.step_all(world, stepping, scratch, recorder, lambda w, myr: take(myr) if round(myr, 6) in wanted else None)
    args.out.write_text(json.dumps({"seed": args.seed, "setting": args.setting, "step_myr": args.step_myr, "flow_scale": args.flow_scale, "samples": samples,
                                    "wall_s": time.perf_counter() - started}, default=float) + "\n")


def queue(args: argparse.Namespace) -> None:
    seeds = args.seeds or json.loads(args.seeds_file.read_text())["seeds"]
    args.out.mkdir(parents=True, exist_ok=True)
    failed = []
    for seed in seeds:
        target = args.out / f"{seed}.json"
        if target.exists():
            continue
        command = [sys.executable, __file__, "run", "--settings", str(args.settings), "--setting", args.setting, "--seed", str(seed),
                   "--step-myr", str(args.step_myr), "--flow-scale", str(args.flow_scale), "--out", str(target), "--at", *map(str, args.at)]
        with open(args.out / f"{seed}.log", "w") as log:
            code = subprocess.call(command, stdout=log, stderr=subprocess.STDOUT)
        print(seed, "exit", code, time.strftime("%H:%M"), flush=True)
        if code:
            failed.append(seed)
    sys.exit(1 if failed else 0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("run", run_one), ("queue", queue)):
        p = sub.add_parser(name)
        p.add_argument("--settings", type=Path, required=True)
        p.add_argument("--setting", required=True)
        p.add_argument("--at", type=float, nargs="+", default=[0, 1, 2, 5, 10])
        p.add_argument("--step-myr", type=float, default=1.0)
        p.add_argument("--flow-scale", type=float, default=1.0, help="multiply every mantle convection centre's strength (default 1)")
        p.add_argument("--out", type=Path, required=True)
        if name == "run":
            p.add_argument("--seed", type=int, required=True)
        else:
            p.add_argument("--seeds", type=int, nargs="*")
            p.add_argument("--seeds-file", type=Path, default=Path("analysis/issue958/development-seeds.json"))
        p.set_defaults(fn=fn)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()

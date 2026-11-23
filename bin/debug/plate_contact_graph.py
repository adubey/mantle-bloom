"""Plate contact graphs, Earth and the model, with the relative speed on each edge (#995).

Node = a plate (sized by area); edge = two plates in contact (width by the length of contact); label = the relative
speed of the two plates along their contact, cm/yr, length-weighted; colour = the contact's mean normal rate against
its speed (red closing, blue opening, grey sliding: |normal| under 30% of speed).

    record   step one generated world plates-only and write its contact graphs (compute-budget first; about 1 GB):
             backend/.venv/bin/python bin/debug/plate_contact_graph.py record --settings analysis/issue975/settings.json \\
                 --setting b123-brake-off --seed 591453549 --years 100 --every 50 --out analysis/issue995/graph-591453549.json
    earth    the same graph from PB2002 (steps file downloaded beforehand into its own directory; plate areas from
             `analysis/issue958/earth-plate-speed.json`):
             python3 -I bin/debug/plate_contact_graph.py earth --steps <dir>/PB2002_steps.dat.txt \\
                 --plates analysis/issue958/earth-plate-speed.json --out analysis/issue995/graph-earth.json
    draw     needs matplotlib and networkx (not in the repo venv; use a scratch venv):
             python -I bin/debug/plate_contact_graph.py draw --graphs graph-earth.json graph-591453549.json@100 ... --out figure.png
             `file@myr` picks the sample of a recorded world.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
MIN_CONTACT_KM = 500.0


def record_graph(world) -> dict:
    from app import lithosphere, torque
    from app.elevation_lines import line_spacing_rad

    omegas = torque.solve_plates_quasistatic(world, diagnostics={})
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
    forms = torque.force_forms(world)
    cache = world.boundary_search_cache
    scale = lithosphere.PLANET_RADIUS_M * 100.0
    side = defaultdict(lambda: [0.0, 0.0, 0.0])  # (plate, neighbour) -> length m, length x speed, length x normal
    local: dict[int, dict] = {}  # per plate: nodes, areas, omega, for the area-weighted speeds (T1, #995)
    for plate in plates:
        neighbours = cache.neighbours(plate, plates, reach)
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
        torque.apply_force_forms(world, plate, neighbours, inputs, forms)
        local[plate.plate_id] = {"points": inputs.own_points, "areas": np.asarray(plate.accounting_areas_m2(spacing), dtype=float),
                                 "omega": omegas[plate.plate_id]}
        band = (inputs.dist_to_neighbor <= reach) & (inputs.neighbor_plate_id >= 0)
        if not band.any():
            continue
        length = torque.boundary_node_length_m(inputs, spacing)
        neighbour_omega = np.zeros_like(inputs.neighbor_omega)
        for other in plates:
            neighbour_omega[inputs.neighbor_plate_id == other.plate_id] = omegas[other.plate_id]
        v_rel = (np.cross(omegas[plate.plate_id], inputs.own_points) - np.cross(neighbour_omega, inputs.own_points)) * scale
        normal = np.sum(v_rel * inputs.direction_to_neighbor, axis=-1)
        speed = np.linalg.norm(v_rel, axis=-1)
        for q in np.unique(inputs.neighbor_plate_id[band]):
            m = band & (inputs.neighbor_plate_id == q)
            acc = side[(plate.plate_id, int(q))]
            acc[0] += length[m].sum()
            acc[1] += (length[m] * speed[m]).sum()
            acc[2] += (length[m] * normal[m]).sum()
    pairs = []
    for a, b in sorted({tuple(sorted(k)) for k in side}):
        both = [side[k] for k in ((a, b), (b, a)) if k in side]
        total = sum(s[0] for s in both)
        pairs.append({"a": a, "b": b, "contact_km": float(np.mean([s[0] for s in both]) / 1e3),
                      "speed_cm_per_yr": float(sum(s[1] for s in both) / total), "normal_cm_per_yr": float(sum(s[2] for s in both) / total)})
    areas = {p.plate_id: float(p.node_areas_m2().sum() / 1e6) for p in plates}
    from plate_torque_breakdown import _net_rotation

    net = _net_rotation(local)
    speeds = {k: float(np.dot(np.linalg.norm(np.cross(v["omega"] - net, v["points"]), axis=1) * scale, v["areas"]) / v["areas"].sum())
              for k, v in local.items()}
    total = sum(v["areas"].sum() for v in local.values())
    mean_nnr = sum(speeds[k] * v["areas"].sum() for k, v in local.items()) / total
    return {"plates": {str(k): {"area_km2": v, "speed_nnr_cm_per_yr": speeds.get(k)} for k, v in areas.items()}, "pairs": pairs,
            "mean_speed_nnr_cm_per_yr": float(mean_nnr)}


def run_record(args: argparse.Namespace) -> None:
    sys.path[:0] = [str(HERE), str(HERE.parents[1] / "backend")]
    import calibrate_plate_forces as cal
    import rift_line_recorder as rec
    import run_rift_lines_d1 as d1
    from app import rift_damage

    setting = cal.load_setting(args.settings, args.setting)
    if args.seed in cal.refused_seeds():
        sys.exit("held-out seed")
    world = d1.make_world(SimpleNamespace(seed=args.seed, node_density=1.0, world_knobs=setting["world_knobs"],
                                          split_method=setting["split_method"], params=setting["params"]))
    cal.apply_module_knobs(setting["module_knobs"])
    recorder = rec.Recorder(rift_damage.runtime(world).params.eps_cand) if setting["split_method"] == "rift_lines" else None
    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    samples = []

    def on_step(world, myr: float) -> None:
        if round(myr / args.every) * args.every == myr:
            samples.append({"myr": myr, **record_graph(world)})
            print(myr, len(samples[-1]["plates"]), "plates", len(samples[-1]["pairs"]), "pairs", flush=True)

    splits, *_ = d1.step_all(world, SimpleNamespace(years=args.years, step_myr=args.step_myr, max_steps=None, snapshot_every=0.0), scratch, recorder, on_step)
    args.out.write_text(json.dumps({"label": f"model seed {args.seed}", "seed": args.seed, "setting": args.setting, "splits": splits, "samples": samples}) + "\n")


def run_earth(args: argparse.Namespace) -> None:
    acc = defaultdict(lambda: [0.0, 0.0, 0.0])
    for line in args.steps.read_text().splitlines():
        t = line.split()
        if len(t) < 14:
            continue
        pair = t[1].lstrip(":")
        for sep in ("\\", "/", "-"):
            if sep in pair:
                a, b = pair.split(sep)
                break
        length, speed, div = float(t[6]), float(t[8]) / 10.0, float(t[10]) / 10.0
        s = acc[tuple(sorted((a, b)))]
        s[0] += length
        s[1] += length * speed
        s[2] += length * -div  # closing is positive, as in the model
    plates = json.loads(args.plates.read_text())["plates"]
    pairs = [{"a": a, "b": b, "contact_km": s[0], "speed_cm_per_yr": s[1] / s[0], "normal_cm_per_yr": s[2] / s[0]} for (a, b), s in acc.items()]
    args.out.write_text(json.dumps({"label": "Earth (PB2002)", "samples": [{"myr": 0, "plates": {k: {"area_km2": v["area_share"]} for k, v in plates.items()}, "pairs": pairs}]}) + "\n")


def run_draw(args: argparse.Namespace) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx

    items = []
    for spec in args.graphs:
        path, _, myr = spec.partition("@")
        data = json.loads(Path(path).read_text())
        samples = data["samples"]
        pick = next((s for s in samples if myr and s["myr"] == float(myr)), samples[-1])
        items.append((data["label"] + (f", {pick['myr']:g} Myr" if myr else ""), pick))
    cols = 2
    fig, axes = plt.subplots((len(items) + 1) // cols, cols, figsize=(8.5 * cols, 7.5 * ((len(items) + 1) // cols)))
    for ax, (title, sample) in zip(np.atleast_1d(axes).ravel(), items):
        earth = "Earth" in title
        areas = {k: v["area_km2"] for k, v in sample["plates"].items()}
        total = sum(areas.values())
        share = {k: v / total for k, v in areas.items()}
        keep = {k for k, v in share.items() if v >= (0.01 if earth else 0.0)}
        edges = [p for p in sample["pairs"] if str(p["a"]) in keep and str(p["b"]) in keep and p["contact_km"] >= MIN_CONTACT_KM]
        g = nx.Graph()
        for p in edges:
            g.add_edge(str(p["a"]), str(p["b"]), **p)
        g.remove_nodes_from([n for n in list(g) if g.degree(n) == 0])
        pos = nx.kamada_kawai_layout(g, weight=None)
        colours = ["#c0392b" if e["normal_cm_per_yr"] > 0.3 * e["speed_cm_per_yr"] else "#2e6fb7" if e["normal_cm_per_yr"] < -0.3 * e["speed_cm_per_yr"] else "#8a8a8a"
                   for _, _, e in g.edges(data=True)]
        widths = [0.6 + 5.0 * e["contact_km"] / max(p["contact_km"] for p in edges) for _, _, e in g.edges(data=True)]
        nx.draw_networkx_edges(g, pos, ax=ax, edge_color=colours, width=widths, alpha=0.8)
        nx.draw_networkx_nodes(g, pos, ax=ax, node_size=[160 + 5200 * share[n] for n in g], node_color="#f2e6c9", edgecolors="#555")
        nx.draw_networkx_labels(g, pos, ax=ax, font_size=8)
        nx.draw_networkx_edge_labels(g, pos, ax=ax, edge_labels={(u, v): f"{e['speed_cm_per_yr']:.1f}" for u, v, e in g.edges(data=True)}, font_size=7,
                                     bbox={"boxstyle": "round,pad=0.1", "fc": "white", "ec": "none", "alpha": 0.8})
        speeds = np.array([e["speed_cm_per_yr"] for _, _, e in g.edges(data=True)])
        weights = np.array([e["contact_km"] for _, _, e in g.edges(data=True)])
        order = np.argsort(speeds)
        median = float(speeds[order][np.searchsorted(np.cumsum(weights[order]) / weights.sum(), 0.5)])
        ax.set_title(f"{title}\n{g.number_of_nodes()} plates, {g.number_of_edges()} contacts, contact-length median {median:.1f} cm/yr", fontsize=11)
        ax.axis("off")
    for ax in np.atleast_1d(axes).ravel()[len(items):]:
        ax.axis("off")
    fig.text(0.5, 0.005, "Edge label: relative speed, cm/yr. Width: contact length. Red: closing, blue: opening, grey: sliding. Node size: area. Contacts under 500 km omitted"
             + (", Earth plates under 1% of area omitted." if any("Earth" in t for t, _ in items) else "."), ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(args.out, dpi=130)


def main() -> None:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("record")
    r.add_argument("--settings", type=Path, required=True)
    r.add_argument("--setting", required=True)
    r.add_argument("--seed", type=int, required=True)
    r.add_argument("--years", type=float, default=100.0)
    r.add_argument("--step-myr", type=float, default=1.0)
    r.add_argument("--every", type=float, default=50.0)
    r.add_argument("--scratch", default="/tmp/unused")
    r.add_argument("--out", type=Path, required=True)
    e = sub.add_parser("earth")
    e.add_argument("--steps", type=Path, required=True)
    e.add_argument("--plates", type=Path, required=True)
    e.add_argument("--out", type=Path, required=True)
    d = sub.add_parser("draw")
    d.add_argument("--graphs", nargs="+", required=True)
    d.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    {"record": run_record, "earth": run_earth, "draw": run_draw}[args.command](args)


if __name__ == "__main__":
    main()

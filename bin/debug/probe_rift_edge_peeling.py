"""Measure where the RFC 0018 damage-band prototype loads, strains and cuts, relative to plate edges (#894).

Runs the prototype (``rift_damage_step``) plates-only on a development seed and records, without changing it:

* every step, by hop distance from the plate boundary: cells, straining cells (oceanic / continental), cells at or
  above ``eps_cand``, broken cells;
* at snapshot times, per plate and hop bin: quantiles of the extensional load, the strength and their ratio, the
  inherited-weakness share, each force component's pull along the principal direction, the share of straining
  cells near a divergent boundary, the thin-shell load (optional) and the straight-section profile (total pull
  across each straight section against that section's integrated strength, with the far-side area share);
* for every breakthrough path search: band size and depth, path length and depth, and how far apart its two ends
  are compared with the plate's span.

Usage (repo root, this worktree's venv; compute-budget skill first, one heavy run at a time):
    backend/.venv/bin/python bin/debug/probe_rift_edge_peeling.py --seed 106533960 --years 60 --out OUT.json \
        [--snapshot-every 10] [--shell-at 1 30] [--setting '{"k_cap": 1.0}'] [--maps DIR --map-at 30 60]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import dijkstra

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import calibrate_rift_load_proxy as d0  # noqa: E402
import rift_damage_mechanism as mech  # noqa: E402
import rift_damage_prototype as proto  # noqa: E402
import rift_damage_step as st  # noqa: E402
import run_rift_damage_d1 as d1  # noqa: E402

from app import merge_split  # noqa: E402

HOP_EDGES = np.array([0, 1, 2, 3, 5, 10, 10**6])  # bins: 0, 1, 2, 3-4, 5-9, 10+
HOP_LABELS = ["0", "1", "2", "3-4", "5-9", "10+"]
COMPONENTS = ("slab_pull_n", "slab_drag_n", "ridge_push_n", "basal_drag_n", "contact_drag_n")
RADIUS_M = mech.RADIUS_M


def hop_bin(hops):
    return np.clip(np.searchsorted(HOP_EDGES, hops, side="right") - 1, 0, len(HOP_LABELS) - 1)


def hops_from(graph, sources):
    if not len(sources):
        return np.full(graph.shape[0], 10**6)
    out = dijkstra(graph, directed=False, indices=sources, unweighted=True, min_only=True)
    return np.where(np.isfinite(out), out, 10**6).astype(int)


def quantiles(values):
    if not len(values):
        return None
    return [float(q) for q in np.percentile(values, [50, 90, 99])]


def along(tensor, angle):
    c, s = np.cos(angle), np.sin(angle)
    return tensor[:, 0] * c * c + tensor[:, 1] * s * s + 2.0 * tensor[:, 2] * s * c


def section_profile(points, areas, force_n, n_eff, spacing_m, k_cut):
    """For straight sections in the six proxy directions: total pull across the section (N), the section's
    integrated strength (N), their ratio, and the smaller side's area share. Returns the best section overall
    and the best with at least 20% of the plate's area on the smaller side."""
    balanced = proto.balanced_tangent_force(points, areas, force_n)
    centre = points.mean(axis=0)
    centre /= np.linalg.norm(centre)
    e_u, e_v = (axis[0] for axis in st.thin_shell._tangent_basis(centre[None, :]))
    coords = RADIUS_M * np.column_stack([points @ e_u, points @ e_v])
    force_2d = np.column_stack([balanced @ e_u, balanced @ e_v])
    total_area = areas.sum()
    best = {"any": None, "interior": None, "pull_peak": None}
    for angle in proto.CUT_ANGLES:
        normal = np.array([np.cos(angle), np.sin(angle)])
        position, push = coords @ normal, force_2d @ normal
        order = np.argsort(position, kind="stable")
        ranked = position[order]
        beyond = np.empty(len(position))
        beyond[order] = np.cumsum(push[order][::-1])[::-1] - 0.5 * push[order]
        pull = k_cut * beyond  # N across the whole section
        cum_strength = np.concatenate([[0.0], np.cumsum((n_eff * areas)[order])])
        cum_area = np.concatenate([[0.0], np.cumsum(areas[order])])
        upper = np.searchsorted(ranked, position + 0.5 * spacing_m, side="right")
        lower = np.searchsorted(ranked, position - 0.5 * spacing_m, side="left")
        strength = (cum_strength[upper] - cum_strength[lower]) / spacing_m  # N
        rank = np.empty(len(position), dtype=int)
        rank[order] = np.arange(len(position))
        far = 1.0 - cum_area[rank] / total_area
        small = np.minimum(far, 1.0 - far)
        ratio = np.where(strength > 0.0, pull / np.maximum(strength, 1.0), 0.0)
        for name, mask, score in (("any", np.ones(len(ratio), bool), ratio), ("interior", small >= 0.2, ratio),
                                  ("pull_peak", np.ones(len(ratio), bool), pull)):
            if not mask.any():
                continue
            k = int(np.flatnonzero(mask)[np.argmax(score[mask])])
            row = {"ratio": float(ratio[k]), "pull_n": float(pull[k]), "strength_n": float(strength[k]),
                   "small_share": float(small[k]), "angle_deg": float(np.degrees(angle))}
            key = "pull_n" if name == "pull_peak" else "ratio"
            if best[name] is None or row[key] > best[name][key]:
                best[name] = row
    return best


class Probe:
    def __init__(self, params, snapshot_every, shell_at):
        self.params = params
        self.snapshot_every = snapshot_every
        self.shell_at = set(shell_at)
        self.per_step = []  # one row per step: counts by hop bin
        self.snapshots = []
        self.paths = []
        self.map_dir, self.map_at = None, set()
        self.hops = {}  # plate_id -> (hops from boundary, span_m)
        self.row = None
        self.myr = 0

    def start_step(self, myr):
        self.myr = myr
        z = lambda: [0] * len(HOP_LABELS)  # noqa: E731
        self.row = {"myr": myr, "cells": z(), "straining_cont": z(), "straining_ocean": z(), "cand": z(), "broken": z(),
                    "max_hop": []}

    def end_step(self):
        self.per_step.append(self.row)

    def after_damage(self, world, plate, step):
        boundary = np.flatnonzero(plate.exposed_sides().any(axis=1))
        hops = hops_from(step.graph, boundary)
        bins = hop_bin(hops)
        if len(boundary) > 1:
            sample = step.points[boundary[:: max(1, len(boundary) // 200)]]
            span = RADIUS_M * float(np.arccos(np.clip((sample @ sample.T).min(), -1.0, 1.0)))
        else:
            span = 0.0
        self.hops[plate.plate_id] = (hops, span)
        straining = step.d_eps > 0.0
        count = lambda mask: np.bincount(bins[mask], minlength=len(HOP_LABELS)).tolist()  # noqa: E731
        for key, mask in (("cells", np.ones(len(bins), bool)), ("straining_cont", straining & ~step.oceanic),
                          ("straining_ocean", straining & step.oceanic), ("cand", step.eps >= self.params.eps_cand),
                          ("broken", step.broken)):
            self.row[key] = [a + b for a, b in zip(self.row[key], count(mask))]
        self.row["max_hop"].append(int(hops[hops < 10**6].max()) if len(hops) else 0)
        if self.myr == 1 or (self.snapshot_every and self.myr % self.snapshot_every == 0):
            self.snapshot(world, plate, step, hops, bins, straining)

    def snapshot(self, world, plate, step, hops, bins, straining):
        p = self.params
        lam, _, angle = proto.principal(step.tensor)
        load = np.maximum(lam, 0.0)
        ratio = load / np.maximum(step.n_eff, 1.0)
        field, divergent = st.force_inputs(world, plate)
        spacing_m = st.line_spacing_rad(world.node_density) * RADIUS_M
        div_hops = hops_from(step.graph, np.flatnonzero(divergent))
        parts = {name: along(p.k_cut * proto.cut_tensor(step.points, step.areas, getattr(field, name), spacing_m), angle)
                 for name in COMPONENTS}
        shell = None
        if self.myr in self.shell_at:
            tensor = d0.solve_shell(plate, step.points, step.areas, field.total_n)
            if tensor is not None:
                shell = np.maximum(proto.principal(tensor)[0], 0.0)
        rows = []
        for b, label in enumerate(HOP_LABELS):
            sel = bins == b
            if not sel.any():
                continue
            row = {"hops": label, "cells": int(sel.sum()), "load_n_per_m": quantiles(load[sel]), "strength_n_per_m": quantiles(step.n_eff[sel]),
                   "ratio": quantiles(ratio[sel]), "over_strength_share": float((ratio[sel] > 1.0).mean()),
                   "oceanic_share": float(step.oceanic[sel].mean()), "weak_share": float(step.weak[sel].mean()),
                   "component_median_n_per_m": {name: float(np.median(parts[name][sel])) for name in COMPONENTS}}
            if shell is not None and np.isfinite(shell[sel]).any():
                row["shell_load_n_per_m"] = quantiles(shell[sel][np.isfinite(shell[sel])])
            rows.append(row)
        min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
        if self.map_dir is not None and self.myr in self.map_at:
            np.savez_compressed(self.map_dir / f"t{self.myr:03d}-plate{int(plate.plate_id)}.npz", points=step.points.astype(np.float32), eps=step.eps.astype(np.float32),
                                broken=step.broken, oceanic=step.oceanic, weak=step.weak, load=load.astype(np.float32), strength=step.n_eff.astype(np.float32),
                                angle=angle.astype(np.float32), hops=hops.astype(np.int32), craton=plate.collect("craton_crust_m") > 0.0,
                                elevation=step.elevation.astype(np.float32), divergent=divergent)
        self.snapshots.append({
            "myr": self.myr, "plate_id": int(plate.plate_id), "cells": int(len(hops)), "max_hop": int(hops[hops < 10**6].max()),
            "oceanic_share": float(step.oceanic.mean()), "min_nodes": min_nodes, "by_hops": rows,
            "straining": int(straining.sum()),
            "straining_within_3_of_divergent": int((straining & (div_hops <= 3)).sum()),
            "straining_oceanic": int((straining & step.oceanic).sum()),
            "cells_within_3_of_divergent": int((div_hops <= 3).sum()),
            "sections": section_profile(step.points, step.areas, field.total_n, step.n_eff, spacing_m, p.k_cut),
        })

    def record_path(self, step, boundary, band, chain):
        hops, span = self.hops.get(step.plate_id, (None, 0.0))
        row = {"myr": self.myr, "plate_id": int(step.plate_id), "plate_cells": int(len(step.points)), "band_cells": int(band.sum()),
               "band_oceanic_share": float(step.oceanic[band].mean()), "found": chain is not None}
        if hops is not None:
            row.update(band_mean_hops=float(hops[band].mean()), band_max_hops=int(hops[band].max()), plate_max_hop=int(hops[hops < 10**6].max()))
        if chain is not None:
            a, b = step.points[chain[0]], step.points[chain[-1]]
            row.update(path_cells=len(chain), end_separation_m=RADIUS_M * float(np.arccos(np.clip(a @ b, -1.0, 1.0))), plate_span_m=span,
                       path_intact_cells=int((~step.broken[chain]).sum()))
            if hops is not None:
                row["path_max_hops"] = int(hops[chain].max())
            row.update(self.path_outcome(step, boundary, band, chain))
        self.paths.append(row)

    def path_outcome(self, step, boundary, band, chain):
        """What the breakthrough rule does with this chain (mirrors ``rift_damage_step.breakthrough``), plus the band's
        shape: elongation (ratio of its two principal extents) and the angle between its long axis and the direction
        of largest tension (90 degrees = a band striking normal to the pull, as a rift should)."""
        cells = np.flatnonzero(band)
        centred = step.points[cells] - step.points[cells].mean(axis=0)
        _, sing, vt = np.linalg.svd(centred, full_matrices=False)
        first, second = st.thin_shell._tangent_basis(step.points)
        _, _, angle = proto.principal(step.tensor[cells])
        tension = np.cos(angle)[:, None] * first[cells] + np.sin(angle)[:, None] * second[cells]
        cosine = np.abs(tension @ vt[0])
        out = {"band_elongation": float(sing[0] / max(sing[1], 1e-12)), "band_axis_to_tension_deg": float(np.degrees(np.arccos(np.clip(np.median(cosine), 0, 1)))),
               "overlap": len(set(chain)) != len(chain)}
        if out["overlap"] or len(chain) < 2 or not boundary[chain[0]] or not boundary[chain[-1]] or chain[0] == chain[-1]:
            out["outcome"] = "overlap" if out["overlap"] else "bad_ends"
            return out
        path = np.asarray(chain)
        normal, ds = st._path_normals(step.points, path, first, second)
        pulling = np.maximum(0.0, proto.normal_component(step.tensor[path], normal))
        load, strength = float(np.sum(pulling * ds)), float(np.sum((step.n_eff[path] * ds)[~step.broken[path]]))
        out.update(ligament_load_n=load, ligament_strength_n=strength, path_length_m=float(ds.sum()))
        if not (load > 0.0 and load >= strength):
            out["outcome"] = "ligament_fail"
            return out
        off = np.ones(len(step.points), dtype=bool)
        off[path] = False
        count, labels = st.connected_components(step.graph[off][:, off], directed=False)
        sizes = sorted(np.bincount(labels).tolist(), reverse=True)
        out.update(outcome="ligament_pass", components=int(count), component_sizes=sizes[:4])
        return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--years", type=int, default=60, help="Myr at 1 Myr steps")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--snapshot-every", type=int, default=10)
    parser.add_argument("--shell-at", type=int, nargs="*", default=[])
    parser.add_argument("--setting", type=json.loads, default={})
    parser.add_argument("--maps", type=Path, help="directory for per-plate map arrays at the --map-at times (snapshot times only)")
    parser.add_argument("--map-at", type=int, nargs="*", default=[])
    args = parser.parse_args()
    if args.seed in d1.forbidden_seeds():
        sys.exit(f"seed {args.seed} is held-out or reserve; development seeds only")
    params = mech.Params(**args.setting)
    state = st.install(params)
    probe = Probe(params, args.snapshot_every, args.shell_at)
    if args.maps:
        args.maps.mkdir(parents=True, exist_ok=True)
        probe.map_dir, probe.map_at = args.maps, set(args.map_at)
    original_damage, original_path = st.damage_step, st.band_path

    def damage_step(world, plate, years, mechanism):
        step = original_damage(world, plate, years, mechanism)
        if step is not None:
            probe.after_damage(world, plate, step)
        return step

    def band_path(step, boundary, band):
        chain = original_path(step, boundary, band)
        probe.record_path(step, boundary, band, chain)
        return chain

    st.damage_step, st.band_path = damage_step, band_path
    from app import world as world_module

    world = world_module.generate_world(args.seed)
    world.simulate_climate_biomes = False
    started = time.perf_counter()
    for myr in range(1, args.years + 1):
        probe.start_step(myr)
        world_module.step_world(world, 1e6)
        probe.end_step()
        if myr % 10 == 0:
            print(json.dumps({"myr": myr, "plates": len(world.plates), "counters": state.counters, "wall_s": round(time.perf_counter() - started, 1)}), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "seed": args.seed, "years_myr": args.years, "params": args.setting, "hop_labels": HOP_LABELS, "counters": state.counters,
        "splits": len(state.events), "per_step": probe.per_step, "snapshots": probe.snapshots, "paths": probe.paths,
    }, default=float) + "\n")
    print("wrote", args.out)


if __name__ == "__main__":
    main()

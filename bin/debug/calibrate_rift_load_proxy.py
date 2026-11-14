"""RFC 0018 D0: calibrate the cheap rift-load tensor against the thin-shell diagnostic (#797).

Usage (repo root, this worktree's venv; follow the compute-budget skill first):
    backend/.venv/bin/python bin/debug/calibrate_rift_load_proxy.py --out analysis/issue797/d0-load-calibration.json \
        [--snapshots <rfc0018-g0-snapshots> --seeds N ...] [--saves <a.mbworld> ...]

For every plate of each world, this builds the per-node force field (``torque.gather_per_node_force_field``
through ``replay_gate2_capacity._event_force_field``), solves the thin-shell membrane diagnostic
(``thin_shell.solve_thin_shell``) for its resultant tensor, and builds the cheap tensor of
``rift_damage_prototype`` for a grid of ``L_pull`` and ``L_drag``. It picks the pair with the smallest median
|log10(N / N_ext)| over loaded extensional cells, then reports D0's frozen criteria at that pair:

* manufactured sign cases on one real plate mesh (one-sided slab pull, two slabs pulling apart, a slab
  opposed by collision, pure collision): the proxy's sign at the plate centre must match the shell's;
* a ligament case: tension parallel to an intact path must not split it, and tension across it must;
* sign agreement (extension versus none) ≥ 80% of loaded cells; principal direction within 30° on ≥ 70% of
  cells extensional in both; Spearman(N, N_ext) ≥ 0.5; top-decile median N / N_ext in [0.5, 2] for
  each crust kind.

A cell is *loaded* when its largest absolute shell principal value is at least 10% of the plate's 95th
percentile, and *extensional* when its largest principal value exceeds 5% of that scale (each field's own
scale). Plume load is not in the shell's force field, so ``N_plume`` is not calibrated here.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parent), str(HERE.parents[2] / "backend")]

import compare_rift_load_maps as maps  # noqa: E402
import replay_gate2_capacity as gate2  # noqa: E402
import rift_damage_prototype as proto  # noqa: E402

from app import elevation_lines, persistence, thin_shell, torque  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402

L_PULL_KM = (250.0, 500.0, 1000.0, 2000.0, 4000.0)
L_DRAG_M = tuple(float(v) for v in np.logspace(4, 10, 25))  # up to 1e10 m: a local drag difference stands in for drag integrated across the plate


def adjacency_pairs(plate):
    adjacency = plate.adjacency()
    counts = np.diff(adjacency.offsets)
    i = np.repeat(np.arange(len(counts)), counts)
    j = np.asarray(adjacency.neighbours)
    keep = i < j
    return i[keep], j[keep]


def solve_shell(plate, points, areas, total_n):
    """Thin-shell resultant tensor (N, 3) in each node's tangent basis; NaN off the mesh. None without a mesh."""
    mesh = maps.plate_mesh(plate, points)
    if mesh is None:
        return None
    used, triangles = mesh
    hc = plate.collect("crustal_thickness_m")[used]
    hm = plate.collect("mantle_lithosphere_thickness_m")[used]
    result = thin_shell.solve_thin_shell(points[used], areas[used], total_n[used], thickness_m=np.maximum(hc + hm, 1.0),
                                         young_modulus_pa=gate2.EFFECTIVE_MODULUS_PA,
                                         extensional_capacity_resultant_n_per_m=np.ones(len(used)), mesh_triangles=triangles)
    out = np.full((len(points), 3), np.nan)
    out[used] = result.stress_resultant_n_per_m
    return out


def plate_record(world, plate):
    points, _ = plate.all_points_and_elevation()
    field = gate2._event_force_field(world, plate)
    started = time.perf_counter()
    shell = solve_shell(plate, points, field.areas_m2, field.total_n)
    shell_s = time.perf_counter() - started
    if shell is None:
        return None
    i, j = adjacency_pairs(plate)
    started = time.perf_counter()
    parts = proto.load_parts(points, field.areas_m2, i, j, field)
    kernels = {l_km: {name: proto.kernel_tensor(points, src, l_km * 1e3) for name, src in parts.sources.items()}
               for l_km in L_PULL_KM}
    proxy_s = (time.perf_counter() - started) / len(L_PULL_KM)
    continental = np.asarray(elevation_lines.effective_is_continental_from_codes(
        plate.collect("crust_type_code"), plate.crust_type == "continental"), bool)
    return {"plate_id": int(plate.plate_id), "nodes": len(points), "shell": shell, "parts": parts, "kernels": kernels,
            "continental": continental, "shell_s": shell_s, "proxy_s": proxy_s}


def masks(record, proxy):
    """Loaded and extensional masks for one plate; NaN shell cells are excluded."""
    s1, s2, _ = proto.principal(record["shell"])
    p1, p2, _ = proto.principal(proxy)
    finite = np.isfinite(s1)
    scale_s = np.nanpercentile(np.maximum(np.abs(s1), np.abs(s2)), 95)
    scale_p = np.percentile(np.maximum(np.abs(p1), np.abs(p2)), 95)
    loaded = finite & (np.maximum(np.abs(s1), np.abs(s2)) >= 0.1 * scale_s)
    return loaded, s1 > 0.05 * scale_s, p1 > 0.05 * max(scale_p, 1e-30)


def pooled_metrics(records, l_pull_km, l_drag_m):
    rows = {"loaded": [], "ext_s": [], "ext_p": [], "s1": [], "p1": [], "ang_s": [], "ang_p": [], "cont": []}
    for record in records:
        proxy = proto.assemble(record["parts"], record["kernels"][l_pull_km], l_drag_m)
        loaded, ext_s, ext_p = masks(record, proxy)
        s1, _, ang_s = proto.principal(record["shell"])
        p1, _, ang_p = proto.principal(proxy)
        for key, value in (("loaded", loaded), ("ext_s", ext_s), ("ext_p", ext_p), ("s1", s1), ("p1", p1),
                           ("ang_s", ang_s), ("ang_p", ang_p), ("cont", record["continental"])):
            rows[key].append(value[loaded] if key != "loaded" else value[loaded])
    return {k: np.concatenate(v) for k, v in rows.items()}


def score(m) -> dict:
    both = m["ext_s"] & m["ext_p"]
    diff = np.abs(((m["ang_p"] - m["ang_s"] + np.pi / 2) % np.pi) - np.pi / 2)
    out = {"loaded_cells": int(len(m["s1"])), "sign_agreement": float(np.mean(m["ext_s"] == m["ext_p"])),
           "direction_within_30deg": float(np.mean(diff[both] < np.radians(30))) if both.any() else None,
           "spearman": float(spearmanr(np.maximum(m["p1"], 0), np.maximum(m["s1"], 0)).statistic),
           "top_decile_median_ratio": {}}
    for kind, sel in (("continental", m["cont"]), ("oceanic", ~m["cont"])):
        s1 = m["s1"][sel & m["ext_s"]]
        p1 = m["p1"][sel & m["ext_s"]]
        if len(s1) >= 10:
            top = s1 >= np.percentile(s1, 90)
            out["top_decile_median_ratio"][kind] = float(np.median(p1[top] / s1[top]))
    ext = m["ext_s"] & (m["p1"] > 0)
    out["median_abs_log10_ratio"] = float(np.median(np.abs(np.log10(m["p1"][ext] / m["s1"][ext])))) if ext.any() else np.inf
    return out


def passes(result) -> dict:
    ratios = result["top_decile_median_ratio"].values()
    return {"sign_agreement": result["sign_agreement"] >= 0.80,
            "direction": (result["direction_within_30deg"] or 0.0) >= 0.70,
            "spearman": result["spearman"] >= 0.5,
            "magnitude": bool(ratios) and all(0.5 <= r <= 2.0 for r in ratios)}


def fit(records):
    grid = []
    for l_pull in L_PULL_KM:
        for l_drag in L_DRAG_M:
            result = score(pooled_metrics(records, l_pull, l_drag))
            grid.append({"l_pull_km": l_pull, "l_drag_m": l_drag, **result})
    best = min(grid, key=lambda g: g["median_abs_log10_ratio"])
    return best, grid


def centre_category(tensor, centre) -> str:
    p1, p2, _ = proto.principal(tensor[centre])
    scale = np.nanmax(np.abs(np.concatenate([p1, p2])))
    if np.nanmedian(p1) > 0.1 * scale:
        return "tension"
    if np.nanmedian(p2) < -0.1 * scale:
        return "compression"
    return "neutral"


def manufactured_cases(world, plate, l_pull_km, l_drag_m) -> dict:
    """Boundary loads on a real plate; the proxy's centre sign must match the shell's."""
    points, _ = plate.all_points_and_elevation()
    areas = np.asarray(plate.accounting_areas_m2(line_spacing_rad(world.node_density)), float)
    centre_xyz = points.mean(axis=0) / np.linalg.norm(points.mean(axis=0))
    offsets = points - centre_xyz
    axis = np.linalg.svd(offsets - offsets.mean(axis=0), full_matrices=False)[2][0]
    along = offsets @ axis
    boundary = plate.exposed_sides().any(axis=1)
    side_a, side_b = boundary & (along < -0.3 * np.abs(along).max()), boundary & (along > 0.3 * np.abs(along).max())
    centre = np.abs(along) < 0.15 * np.abs(along).max()
    tangent_axis = axis - (points @ axis)[:, None] * points
    tangent_axis /= np.linalg.norm(tangent_axis, axis=1)[:, None]
    force = 1e12  # N per boundary node
    i, j = adjacency_pairs(plate)
    cases = {"one_sided_pull": {"slab": (side_a, -1.0)},
             "two_slabs_apart": {"slab": (side_a | side_b, None)},
             "slab_opposed_by_collision": {"slab": (side_a, -1.0), "contact": (side_b, -1.0)},
             "pure_collision": {"contact": (side_a | side_b, "inward")}}
    out = {}
    for name, loads in cases.items():
        comp = {k: np.zeros((len(points), 3)) for k in ("slab", "contact")}
        for kind, (mask, direction) in loads.items():
            outward = np.where((along < 0)[:, None], -tangent_axis, tangent_axis)
            vec = outward if direction is None else (-outward if direction == "inward" else direction * tangent_axis)
            comp[kind][mask] = force * vec[mask]
        zero = np.zeros((len(points), 3))
        field = torque.PerNodeForceField(points, areas, comp["slab"], zero, zero, zero, comp["contact"])
        shell = solve_shell(plate, points, areas, field.total_n)
        parts = proto.load_parts(points, areas, i, j, field)
        kernels = {k: proto.kernel_tensor(points, src, l_pull_km * 1e3) for k, src in parts.sources.items()}
        proxy = proto.assemble(parts, kernels, l_drag_m)
        out[name] = {"shell": centre_category(shell, centre), "proxy": centre_category(proxy, centre)}
        out[name]["match"] = out[name]["shell"] == out[name]["proxy"]
    return out


def ligament_case() -> dict:
    """Uniform uniaxial tension along the first tangent axis; a path along it must not split, one across it must."""
    n = 20
    tensor = np.tile([1e12, 0.0, 0.0], (n, 1))
    ds, intact, n_eff = np.full(n, 1e5), np.ones(n, bool), np.full(n, 5e11)
    along = proto.ligament_splits(tensor, np.tile([0.0, 1.0], (n, 1)), ds, intact, n_eff)
    across = proto.ligament_splits(tensor, np.tile([1.0, 0.0], (n, 1)), ds, intact, n_eff)
    return {"parallel_path_splits": along, "crossing_path_splits": across, "match": (not along) and across}


def load_worlds(args):
    paths = [(f"snapshot-{s}", args.snapshots.expanduser() / f"g0-snapshot-seed-{s}.bloom") for s in args.seeds]
    paths += [(f"{p.parent.name}/{p.name}", p.expanduser()) for p in args.saves]
    for label, path in paths:
        with open(path, "rb") as handle:
            yield label, persistence.load_world_file(handle)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--snapshots", type=Path)
    parser.add_argument("--seeds", type=int, nargs="*", default=[])
    parser.add_argument("--saves", type=Path, nargs="*", default=[])
    args = parser.parse_args()
    records, first_world, skipped = [], None, 0
    for label, world in load_worlds(args):
        first_world = first_world or world
        for plate in world.plates:
            record = plate_record(world, plate)
            skipped += record is None
            if record is not None:
                records.append({"world": label, **record})
        print(json.dumps({"world": label, "plates": len(records), "skipped": skipped}), flush=True)
    best, grid = fit(records)
    largest = max(first_world.plates, key=lambda p: p.node_count())
    summary = {"plates": len(records), "plates_without_mesh": skipped, "fitted": best, "passes": passes(best),
               "manufactured": manufactured_cases(first_world, largest, best["l_pull_km"], best["l_drag_m"]),
               "ligament": ligament_case(),
               "seconds": {"shell_total": float(sum(r["shell_s"] for r in records)),
                           "proxy_total_per_setting": float(sum(r["proxy_s"] for r in records))},
               "grid": [{k: g[k] for k in ("l_pull_km", "l_drag_m", "sign_agreement", "spearman",
                                           "median_abs_log10_ratio")} for g in grid]}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=1, default=float) + "\n")
    print(json.dumps({k: summary[k] for k in ("plates", "fitted", "passes", "manufactured", "ligament", "seconds")},
                     indent=1, default=float))


if __name__ == "__main__":
    main()

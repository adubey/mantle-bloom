"""Do cheap extensional-load maps line up with where RFC 0018 finds valid cuts? (#876, step 2)

Usage (repo root, this worktree's venv; follow the compute-budget skill first):
    backend/.venv/bin/python bin/debug/compare_rift_load_maps.py \
        --g0-script <measure_rfc0018_g0.py> --out analysis/issue876/load-maps.json --rows <rows.jsonl> \
        [--snapshots <rfc0018-g0-snapshots> --seeds N ...] [--saves <a.mbworld> ...]

For every split opportunity (as ``g0._is_opportunity`` defines it) this builds two maps before any eigensolve:

* **drag map** (no solve): per shared edge, the rate at which the mantle flow beneath the plate, minus the plate's
  own rigid rotation, pulls the two end nodes apart along the edge, per unit edge angle (1/yr). A node takes its
  largest positive incident rate. Basal drag only; divided by rock resistance (RFC 0018's node resistance with its
  flow-residual factor removed, so it measures rock alone: thickness, cratons, belts, plume heat, faults, thinning).
* **shell map** (one sparse least-squares solve): the incremental extensional stress resultant from the #699
  thin-shell membrane diagnostic (``thin_shell.solve_thin_shell``) under the full per-node force field (slab pull,
  slab drag, ridge push, basal drag, contact drag), divided by the column's integrated extensional capacity from
  the #487 yield envelope (``replay_gate2_capacity._column_capacity``, Moho thermal lag included). A ratio of 1
  means the column yields under the load. The capacity has no inherited-weakness term.

It then runs the rebuilt RFC 0018 generator (``rebuild_rfc0018_mass_balanced.propose`` with the exact repeat
shortcut) and records, per candidate, where its cut nodes fall in each map, and per opportunity whether any
candidate is structurally valid or has positive breakup probability. The summary reports how well plate-level and
cut-level map values separate those outcomes (rank AUC, 0.5 = no information) and the share of opportunities a
zero-false-negative threshold would skip. RFC 0018's outcomes are the comparison, not ground truth: a load map
that disagrees with them may be the better physics.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parent), str(HERE.parents[2] / "backend")]

import rebuild_rfc0018_mass_balanced as rebuild  # noqa: E402
import replay_gate2_capacity as gate2  # noqa: E402
import study_rfc0018_skip_predictors as skip_study  # noqa: E402
import validate_thin_shell_plate_meshes as mesh_validation  # noqa: E402

from app import elevation_lines, geometry, mantle, thin_shell, torque  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402


def rock_resistance(g0, world, plate, points):
    """RFC 0018 node resistance with the flow term neutral: a flow that matches one rigid rotation everywhere."""
    rigid = np.cross(plate.omega, points)
    resistance, _, _ = g0._node_resistance(world, plate, points, rigid, np.asarray([plate.omega, plate.omega]))
    return resistance


def drag_map(world, plate, points, i, j):
    """Largest positive pull-apart rate of (mantle flow - plate motion) over each node's edges, 1/yr."""
    residual = mantle.flow_at(points, world.mantle_centers) - np.cross(plate.omega, points)
    chord = points[j] - points[i]
    span = np.maximum(np.linalg.norm(chord, axis=1), 1e-15)
    rate = np.sum((residual[j] - residual[i]) * chord, axis=1) / span**2
    node = np.zeros(len(points))
    np.maximum.at(node, i, rate)
    np.maximum.at(node, j, rate)
    return node


def plate_mesh(plate, points):
    """The #699 clipped plate mesh, or None when it stays disconnected after the driver's repair."""
    try:
        triangles = mesh_validation._candidate_triangles(plate, points)[0]
    except ValueError as error:
        if "tangent hemisphere" not in str(error):
            raise
        triangles = thin_shell._spherical_surface_triangles(points)
        triangles = triangles[plate.contains_batch(geometry.normalize(points[triangles].mean(axis=1)))]
    triangles = gate2._nondegenerate_triangles(points, triangles)
    if len(triangles) and len(gate2._triangle_mesh_component_sizes(triangles)) > 1:
        triangles = gate2._repair_connected_source_mesh(plate, points, triangles)[0]
    if not len(triangles) or len(gate2._triangle_mesh_component_sizes(triangles)) > 1:
        return None
    return gate2._remap_triangle_nodes(triangles)


def column_capacity(world, plate, used):
    hc = plate.collect("crustal_thickness_m")[used]
    hm = plate.collect("mantle_lithosphere_thickness_m")[used]
    lag = plate.collect("moho_thermal_lag_c")[used]
    continental = np.asarray(elevation_lines.effective_is_continental_from_codes(
        plate.collect("crust_type_code"), plate.crust_type == "continental"), bool)[used]
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    cache = getattr(world, "boundary_search_cache", None)
    neighbours = plate.get_neighbours(world.plates)
    inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
    divergent = torque.classify_boundary_nodes(plate, neighbours, inputs, reach, cache)[1][used]
    kinds = np.where(continental, "continental", np.where(divergent, "oceanic_young_ridge", "oceanic_mature"))
    # Quantized as in the #487 driver, only to reuse column solves.
    return np.array([gate2._column_capacity(str(k), round(a / 100.0) * 100.0, round(b / 100.0) * 100.0, round(c))
                     for k, a, b, c in zip(kinds, hc, hm, lag)])


def shell_map(world, plate, points):
    """Extensional resultant over column capacity per node (NaN off the mesh), or a skip reason."""
    field = gate2._event_force_field(world, plate)
    mesh = plate_mesh(plate, points)
    if mesh is None:
        return None, "thin-shell mesh disconnected"
    used, triangles = mesh
    hc = plate.collect("crustal_thickness_m")[used]
    hm = plate.collect("mantle_lithosphere_thickness_m")[used]
    result = thin_shell.solve_thin_shell(
        points[used], field.areas_m2[used], field.total_n[used], thickness_m=np.maximum(hc + hm, 1.0),
        young_modulus_pa=gate2.EFFECTIVE_MODULUS_PA, extensional_capacity_resultant_n_per_m=np.full(len(used), 1e10),
        mesh_triangles=triangles)
    ratio = np.full(len(points), np.nan)
    capacity = column_capacity(world, plate, used)
    ratio[used] = result.max_incremental_extensional_resultant_n_per_m / np.maximum(capacity, 1.0)
    return ratio, None


def percentile_rank(values):
    """Within-plate rank in [0, 1] (NaN stays NaN)."""
    out = np.full(len(values), np.nan)
    finite = np.isfinite(values)
    order = np.argsort(np.argsort(values[finite], kind="stable"), kind="stable")
    out[finite] = order / max(finite.sum() - 1, 1)
    return out


def cut_stats(side, i, j, maps):
    cut = side[i] != side[j]
    nodes = np.unique(np.concatenate([i[cut], j[cut]]))
    stats = {}
    for name, (raw, rank) in maps.items():
        r, v = rank[nodes], raw[nodes]
        r, v = r[np.isfinite(r)], v[np.isfinite(v)]
        stats[f"{name}_rank_median"] = float(np.median(r)) if len(r) else None
        stats[f"{name}_p10"] = float(np.percentile(v, 10)) if len(v) else None  # the weakest-loaded tenth of the line
        stats[f"{name}_median"] = float(np.median(v)) if len(v) else None
    shell = maps["shell"][0][nodes]
    shell = shell[np.isfinite(shell)]
    stats["shell_yield_share"] = float(np.mean(shell >= 1.0)) if len(shell) else None
    return stats


def measure_plate(g0, world, plate):
    points, _ = plate.all_points_and_elevation()
    i, j, _ = g0._shared_edges(plate)
    started = time.perf_counter()
    drag = drag_map(world, plate, points, i, j) / rock_resistance(g0, world, plate, points)
    drag_s = time.perf_counter() - started
    started = time.perf_counter()
    shell, skip = shell_map(world, plate, points)
    shell_s = time.perf_counter() - started
    shell = np.full(len(points), np.nan) if shell is None else shell
    finite = shell[np.isfinite(shell)]
    row = {"plate_id": int(plate.plate_id), "nodes": int(len(points)), "drag_s": drag_s, "shell_s": shell_s,
           "shell_skip": skip, "shell_coverage": float(np.isfinite(shell).mean()),
           "drag_p95": float(np.percentile(drag, 95)),
           "shell_max": float(finite.max()) if len(finite) else None,
           "shell_p95": float(np.percentile(finite, 95)) if len(finite) else None,
           "shell_yield_share": float(np.mean(finite >= 1.0)) if len(finite) else None}
    maps = {"drag": (drag, percentile_rank(drag)), "shell": (shell, percentile_rank(shell))}
    started = time.perf_counter()
    try:
        candidates = rebuild.propose(g0, world, plate, cycle_shortcut=True)
        row["p_break"] = rebuild.breakup_prior(candidates)["P_break"]
    except ValueError as exc:  # the generator raises when no candidate can be proposed
        candidates, row["p_break"], row["status"] = [], 0.0, str(exc)
    row["search_s"] = time.perf_counter() - started
    row["candidates"] = [{"valid": c["structurally_valid"], "p": c.get("p", 0.0), **cut_stats(c["side_mask"], i, j, maps)}
                         for c in candidates]
    row["has_valid"] = any(c["valid"] for c in row["candidates"])
    return row


def auc(scores, labels):
    """Probability that a positive outscores a negative (ties half); None when a class is empty."""
    scores, labels = np.asarray(scores, float), np.asarray(labels, bool)
    keep = np.isfinite(scores)
    scores, labels = scores[keep], labels[keep]
    pos, neg = scores[labels], scores[~labels]
    if not len(pos) or not len(neg):
        return None
    greater = (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
    return float(greater / (len(pos) * len(neg)))


def nan_none(values):
    return np.array([np.nan if v is None else v for v in values], float)


def summarize(rows):
    has_valid = np.array([r["has_valid"] for r in rows])
    positive = np.array([r["p_break"] > 0 for r in rows])
    out = {"opportunities": len(rows), "with_valid": int(has_valid.sum()), "with_positive_p_break": int(positive.sum()),
           "shell_skipped": sum(r["shell_skip"] is not None for r in rows),
           "seconds": {k: float(sum(r[k] for r in rows)) for k in ("drag_s", "shell_s", "search_s")},
           "plate_features": {}, "cut_features": {}}
    for name in ("drag_p95", "shell_max", "shell_p95", "shell_yield_share", "nodes"):
        values = nan_none([r[name] for r in rows])
        finite = np.isfinite(values)
        out["plate_features"][name] = {
            "auc_has_valid": auc(values, has_valid), "auc_positive_p_break": auc(values, positive),
            "zero_fn_skip_has_valid": skip_study.best_zero_fn_threshold(values[finite], has_valid[finite]),
            "opportunities_scored": int(finite.sum())}
    cands = [c for r in rows for c in r["candidates"]]
    valid = np.array([c["valid"] for c in cands])
    p_pos = np.array([c["p"] > 0 for c in cands])
    out["candidates"] = {"total": len(cands), "valid": int(valid.sum()), "positive_p": int(p_pos.sum())}
    for name in ("drag_rank_median", "drag_p10", "drag_median", "shell_rank_median", "shell_p10", "shell_median",
                 "shell_yield_share"):
        values = nan_none([c[name] for c in cands])
        out["cut_features"][name] = {
            "auc_valid": auc(values, valid), "auc_positive_p": auc(values, p_pos),
            "auc_positive_p_among_valid": auc(values[valid], p_pos[valid]),
            "median_valid": float(np.nanmedian(values[valid])) if np.isfinite(values[valid]).any() else None,
            "median_invalid": float(np.nanmedian(values[~valid])) if np.isfinite(values[~valid]).any() else None}
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--g0-script", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--rows", required=True, type=Path)
    parser.add_argument("--snapshots", type=Path)
    parser.add_argument("--seeds", type=int, nargs="*", default=[])
    parser.add_argument("--saves", type=Path, nargs="*", default=[])
    args = parser.parse_args()
    g0 = rebuild.load_g0(args.g0_script.expanduser())
    sources = [(f"snapshot-{s}", args.snapshots.expanduser() / f"g0-snapshot-seed-{s}.bloom") for s in args.seeds]
    sources += [(f"{p.parent.name}/{p.name}", p.expanduser()) for p in args.saves]
    rows = []
    with open(args.rows.expanduser(), "w") as log:
        for label, path in sources:
            with open(path, "rb") as handle:
                world = rebuild.persistence.load_world_file(handle)
            for plate in world.plates:
                if g0._is_opportunity(world, plate):
                    row = {"world": label, **measure_plate(g0, world, plate)}
                    rows.append(row)
                    log.write(json.dumps(row) + "\n")
                    log.flush()
                    print(json.dumps({k: row[k] for k in ("world", "plate_id", "nodes", "has_valid", "shell_skip")}), flush=True)
    summary = summarize(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

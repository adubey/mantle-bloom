"""Summarize a ``probe_rift_edge_peeling.py`` output into the compact numbers reported in #894.

Usage:
    backend/.venv/bin/python bin/debug/summarize_rift_edge_peeling.py PROBE.json [SUMMARY.json]
"""

import json
import sys

import numpy as np


def pooled(snapshots, key, index=0):
    """Cell-weighted mean over plates of a per-hop-bin quantile (index 0 = median, 1 = p90, 2 = p99)."""
    out = {}
    for snap in snapshots:
        for row in snap["by_hops"]:
            value = row.get(key)
            if value is None:
                continue
            value = value[index] if isinstance(value, list) else value
            out.setdefault(row["hops"], []).append((value, row["cells"]))
    return {hop: float(np.average([v for v, _ in rows], weights=[w for _, w in rows])) for hop, rows in out.items()}


def summarize(data):
    labels = data["hop_labels"]
    steps = data["per_step"]
    total = lambda key: np.sum([row[key] for row in steps], axis=0)  # noqa: E731
    cells, cont, ocean = total("cells"), total("straining_cont"), total("straining_ocean")
    straining = cont + ocean
    out = {"seed": data["seed"], "years_myr": data["years_myr"], "counters": data["counters"], "splits": data["splits"], "hop_labels": labels}
    out["cell_share_by_hops"] = (cells / cells.sum()).round(4).tolist()
    out["straining_share_by_hops"] = (straining / max(straining.sum(), 1)).round(4).tolist()
    out["straining_rate_by_hops"] = (straining / np.maximum(cells, 1)).round(4).tolist()
    out["straining_oceanic_share"] = float(ocean.sum() / max(straining.sum(), 1))
    out["straining_cell_steps"] = int(straining.sum())
    out["median_plate_max_hop"] = float(np.median([h for row in steps for h in row["max_hop"]]))
    by_decade = {}
    for row in steps:
        decade = (row["myr"] - 1) // 10
        s = np.array(row["straining_cont"]) + np.array(row["straining_ocean"])
        acc = by_decade.setdefault(decade, np.zeros(len(labels)))
        acc += s
    out["straining_share_10plus_by_decade"] = {f"{10 * d + 1}-{10 * d + 10}": float(v[-1] / max(v.sum(), 1)) for d, v in by_decade.items()}
    out["straining_share_5plus_by_decade"] = {f"{10 * d + 1}-{10 * d + 10}": float(v[-2:].sum() / max(v.sum(), 1)) for d, v in by_decade.items()}
    last = steps[-1]
    out["final_cand_by_hops"], out["final_broken_by_hops"] = last["cand"], last["broken"]

    snaps = data["snapshots"]
    for key, name in (("load_n_per_m", "load"), ("strength_n_per_m", "strength"), ("ratio", "ratio"), ("shell_load_n_per_m", "shell_load")):
        out[f"{name}_median_by_hops"] = pooled(snaps, key, 0)
        out[f"{name}_p99_by_hops"] = pooled(snaps, key, 2)
    out["over_strength_share_by_hops"] = pooled(snaps, "over_strength_share")
    out["oceanic_share_by_hops"] = pooled(snaps, "oceanic_share")
    out["weak_share_by_hops"] = pooled(snaps, "weak_share")
    comps = {}
    for snap in snaps:
        for row in snap["by_hops"]:
            for name, value in row["component_median_n_per_m"].items():
                comps.setdefault(name, {}).setdefault(row["hops"], []).append((value, row["cells"]))
    out["component_median_by_hops"] = {name: {hop: float(np.average([v for v, _ in rows], weights=[w for _, w in rows])) for hop, rows in hops.items()}
                                       for name, hops in comps.items()}
    strain = sum(s["straining"] for s in snaps)
    out["snapshot_straining_within_3_of_divergent_share"] = sum(s["straining_within_3_of_divergent"] for s in snaps) / max(strain, 1)
    out["snapshot_cells_within_3_of_divergent_share"] = sum(s["cells_within_3_of_divergent"] for s in snaps) / sum(s["cells"] for s in snaps)
    out["snapshot_straining_oceanic_share"] = sum(s["straining_oceanic"] for s in snaps) / max(strain, 1)

    big = [s for s in snaps if s["cells"] >= 2 * s["min_nodes"] and s["sections"]["any"] and s["sections"]["interior"]]
    any_small = np.array([s["sections"]["any"]["small_share"] for s in big])
    any_ratio = np.array([s["sections"]["any"]["ratio"] for s in big])
    in_ratio = np.array([s["sections"]["interior"]["ratio"] for s in big])
    peak_small = np.array([s["sections"]["pull_peak"]["small_share"] for s in big])
    out["sections"] = {
        "plate_snapshots": len(big),
        "best_section_small_share_median": float(np.median(any_small)),
        "best_section_small_share_under_10pct": float((any_small < 0.1).mean()),
        "best_section_ratio_quantiles": np.percentile(any_ratio, [10, 50, 90]).tolist(),
        "interior_section_ratio_quantiles": np.percentile(in_ratio, [10, 50, 90]).tolist(),
        "interior_over_best_ratio_median": float(np.median(in_ratio / np.maximum(any_ratio, 1e-30))),
        "interior_ratio_at_least_1_share": float((in_ratio >= 1.0).mean()),
        "interior_ratio_positive_share": float((in_ratio > 0.0).mean()),
        "pull_peak_small_share_median": float(np.median(peak_small)),
        "pull_peak_small_share_at_least_20pct": float((peak_small >= 0.2).mean()),
    }

    paths = [p for p in data["paths"] if p["found"]]
    if paths:
        sep = np.array([p["end_separation_m"] / max(p["plate_span_m"], 1.0) for p in paths])
        out["paths"] = {
            "searched": len(data["paths"]), "found": len(paths),
            "band_cells_quantiles": np.percentile([p["band_cells"] for p in paths], [50, 90, 100]).tolist(),
            "band_mean_hops_quantiles": np.percentile([p["band_mean_hops"] for p in paths], [50, 90, 100]).tolist(),
            "band_max_hops_quantiles": np.percentile([p["band_max_hops"] for p in paths], [50, 90, 100]).tolist(),
            "path_cells_quantiles": np.percentile([p["path_cells"] for p in paths], [50, 90, 100]).tolist(),
            "path_max_hops_quantiles": np.percentile([p["path_max_hops"] for p in paths], [50, 90, 100]).tolist(),
            "end_separation_over_span_quantiles": np.percentile(sep, [50, 90, 100]).tolist(),
            "band_oceanic_share_mean": float(np.mean([p["band_oceanic_share"] for p in paths])),
            "plate_max_hop_median": float(np.median([p["plate_max_hop"] for p in paths])),
        }
    return out


if __name__ == "__main__":
    summary = summarize(json.load(open(sys.argv[1])))
    text = json.dumps(summary, indent=1)
    if len(sys.argv) > 2:
        open(sys.argv[2], "w").write(text + "\n")
    print(text)

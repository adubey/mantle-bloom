#!/usr/bin/env python3
"""Compare RFC 0018 resistance-to-spectral-weight mappings on frozen G0 snapshots.

Uses the candidate and resistance implementation frozen at commit a1860310, loaded
from Git history, and changes only the mapping from that fixed edge cost to Laplacian
conductance. No world is stepped or mutated. G2 held-out seeds are not inputs.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_issue850_spectral_weighting.py \
        --snapshots /path/to/rfc0018-g0-snapshots --out-json /path/to/results.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

from app import persistence


PINNED_SOURCE = "a1860310:bin/debug/measure_rfc0018_g0.py"
MAPPINGS = ("inverse_cost", "direct_cost", "uniform")


def _load_frozen_probe() -> dict:
    source = subprocess.check_output(["git", "show", PINNED_SOURCE], text=True)
    namespace: dict = {"__name__": "_frozen_rfc0018_g0_probe"}
    exec(compile(source, PINNED_SOURCE, "exec"), namespace)
    namespace["__source_sha256__"] = hashlib.sha256(source.encode("utf-8")).hexdigest()
    return namespace


def _mapping_fiedler(mapping: str, counter: dict, captured_graph: dict | None = None):
    def solve(n, i, j, edge_cost):
        counter["total"] += 1
        if counter["total"] == 1:
            counter["base"] += 1
            if captured_graph is not None:
                captured_graph["i"] = np.asarray(i).copy()
                captured_graph["j"] = np.asarray(j).copy()
                captured_graph["edge_cost"] = np.asarray(edge_cost).copy()
        else:
            counter["refinement"] += 1
        if mapping == "inverse_cost":
            conductance = 1.0 / np.maximum(edge_cost, 1.0e-10)
        elif mapping == "direct_cost":
            conductance = np.maximum(edge_cost, 1.0e-10)
        elif mapping == "uniform":
            conductance = np.ones_like(edge_cost)
        else:
            raise ValueError(mapping)
        rows = np.concatenate((i, j))
        cols = np.concatenate((j, i))
        values = np.concatenate((conductance, conductance))
        adjacency = coo_matrix((values, (rows, cols)), shape=(n, n)).tocsr()
        degree = np.asarray(adjacency.sum(axis=1)).ravel()
        inv_root = 1.0 / np.sqrt(np.maximum(degree, 1.0e-15))
        normalized = diags(inv_root) @ (diags(degree) - adjacency) @ diags(inv_root)
        indices = np.arange(1, n + 1, dtype=float)
        v0 = np.cos(indices * 0.6180339887498948) + np.sin(indices * 0.1732050807568877)
        _values, vectors = eigsh(normalized, k=2, which="SM", tol=1.0e-5, v0=v0)
        return vectors[:, 1]
    return solve


def _synthetic_graphs(count: int = 100) -> dict:
    """Compare exact dense Fiedler partitions on deterministic planted-link cycles."""
    rng = np.random.default_rng(850001)
    n = 40
    i = np.arange(n, dtype=np.int64)
    j = (i + 1) % n
    outcomes = {name: {"weak_pair_cut": 0, "any_weak_edge_cut": 0, "graphs": count,
                       "repeat_matches": 0} for name in MAPPINGS}
    for _ in range(count):
        costs = 0.9 + 0.2 * rng.random(n)
        start = int(rng.integers(0, n))
        weak_edges = np.array([start, (start + n // 2) % n])
        costs[weak_edges] = 0.05 + 0.01 * rng.random(2)
        for mapping in MAPPINGS:
            if mapping == "inverse_cost":
                conductance = 1.0 / costs
            elif mapping == "direct_cost":
                conductance = costs
            else:
                conductance = np.ones_like(costs)
            adjacency = np.zeros((n, n), dtype=float)
            adjacency[i, j] = conductance
            adjacency[j, i] = conductance
            degree = adjacency.sum(axis=1)
            inv_root = 1.0 / np.sqrt(degree)
            normalized = (inv_root[:, None] * (np.diag(degree) - adjacency)
                          * inv_root[None, :])
            _eigenvalues, eigenvectors = np.linalg.eigh(normalized)
            first = eigenvectors[:, 1]
            _repeat_values, repeat_vectors = np.linalg.eigh(normalized)
            repeated = repeat_vectors[:, 1]
            order = np.argsort(first, kind="stable")
            side = np.zeros(n, dtype=bool)
            side[order[n // 2:]] = True
            cuts = side[i] != side[j]
            pair_cut = bool(cuts[weak_edges].all())
            outcomes[mapping]["weak_pair_cut"] += int(pair_cut)
            outcomes[mapping]["any_weak_edge_cut"] += int(cuts[weak_edges].any())
            # Eigenvector sign is arbitrary. Compare the induced edge cut.
            repeat_order = np.argsort(repeated, kind="stable")
            repeat_side = np.zeros(n, dtype=bool)
            repeat_side[repeat_order[n // 2:]] = True
            outcomes[mapping]["repeat_matches"] += int(np.array_equal(cuts, repeat_side[i] != repeat_side[j]))
    return outcomes


def _snapshot_comparison(snapshot_dir: Path, manifest: dict, probe: dict) -> tuple[dict, dict]:
    rows = []
    timings = {name: {"runtime_s": 0.0, "base_eigensolves": 0,
                      "candidate_refinement_eigensolves": 0, "opportunities": 0,
                      "successful_candidates": 0, "selected_weak_anchor_costs": []}
               for name in MAPPINGS}
    for entry in manifest["worlds"]:
        if not entry.get("captured"):
            continue
        path = snapshot_dir / entry["file"]
        actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual_hash != entry["sha256"]:
            raise ValueError(f"snapshot hash mismatch: {entry['file']}")
        world = persistence.load_world_bytes(path.read_bytes())
        world_rows_start = len(rows)
        for plate in world.plates:
            if not probe["_is_opportunity"](world, plate):
                continue
            row = {"seed": int(world.seed), "plate_id": int(plate.plate_id),
                   "elapsed_myr": world.elapsed_years / 1e6, "node_count": plate.node_count(),
                   "mappings": {}}
            for mapping in MAPPINGS:
                counter = {"total": 0, "base": 0, "refinement": 0}
                captured_graph = {}
                probe["_fiedler_vector"] = _mapping_fiedler(mapping, counter, captured_graph)
                start = time.perf_counter()
                try:
                    _points, chosen, candidates, _signals = probe["choose_cut"](world, plate)
                    keys = plate.cell_keys.view(np.uint64)
                    selected_pair = tuple(int(value) for value in chosen["seed_edge"])
                    edge_lookup = {(int(keys[a]), int(keys[b])): index
                                   for index, (a, b) in enumerate(zip(captured_graph["i"],
                                                                       captured_graph["j"]))}
                    selected_index = edge_lookup.get(selected_pair)
                    initial_anchor_cost = (float(captured_graph["edge_cost"][selected_index])
                                           if selected_index is not None else None)
                    initial_cost_percentile = (float(np.mean(
                        captured_graph["edge_cost"] <= initial_anchor_cost))
                                               if initial_anchor_cost is not None else None)
                    candidate_rows = []
                    for candidate in candidates:
                        candidate_side = np.asarray(candidate["side_mask"], dtype=bool)
                        packed_side = np.packbits(candidate_side).tobytes()
                        packed_other = np.packbits(~candidate_side).tobytes()
                        mass_a = float(candidate["opening"]["mass_a_kg"])
                        mass_b = float(candidate["opening"]["mass_b_kg"])
                        if mass_a > mass_b or (mass_a == mass_b and packed_side > packed_other):
                            canonical_partition = packed_other
                        else:
                            canonical_partition = packed_side
                        pair = tuple(int(value) for value in candidate["seed_edge"])
                        candidate_edge_index = edge_lookup.get(pair)
                        candidate_anchor_cost = (
                            float(captured_graph["edge_cost"][candidate_edge_index])
                            if candidate_edge_index is not None else None)
                        candidate_rows.append([
                            pair[0], pair[1], hashlib.sha256(canonical_partition).hexdigest(),
                            candidate is chosen, candidate_anchor_cost,
                            (float(np.mean(captured_graph["edge_cost"] <= candidate_anchor_cost))
                             if candidate_anchor_cost is not None else None),
                            bool(candidate["usable_candidate"]),
                        ])
                    selected_candidate_index = next(
                        index for index, candidate in enumerate(candidates) if candidate is chosen)
                    result = {"status": "measured", "candidate_count": len(candidates),
                              "selected_anchor": chosen["seed_edge"],
                              "selected_candidate_index": selected_candidate_index,
                              "selected_anchor_initial_edge_cost": initial_anchor_cost,
                              "selected_anchor_cost_percentile": initial_cost_percentile,
                              "selected_partition_edge_cost": float(chosen["cut_cost"]),
                              "selected_usable": bool(chosen["usable_candidate"]),
                              "selected_usable_count": sum(bool(c["usable_candidate"]) for c in candidates),
                              "selected_mass_fraction": float(chosen["opening"]["smaller_mass_fraction"]),
                              "selected_opening_cm_per_yr": float(chosen["opening"]["net_opening_cm_per_yr"]),
                              "selected_trace_share": float(chosen["largest_trace_length_share"]),
                              "selected_connected": bool(chosen["admissible_connected"]),
                              "candidate_partitions": candidate_rows}
                    timings[mapping]["successful_candidates"] += len(candidates)
                    if initial_anchor_cost is not None:
                        timings[mapping]["selected_weak_anchor_costs"].append(initial_anchor_cost)
                except Exception as exc:
                    result = {"status": "error", "error": f"{type(exc).__name__}: {exc}"}
                finally:
                    elapsed = time.perf_counter() - start
                result["runtime_s"] = elapsed
                result["sparse_eigensolves"] = {"base_partition": counter["base"],
                                                "candidate_refinement": counter["refinement"],
                                                "total": counter["total"]}
                row["mappings"][mapping] = result
                timings[mapping]["runtime_s"] += elapsed
                timings[mapping]["base_eigensolves"] += counter["base"]
                timings[mapping]["candidate_refinement_eigensolves"] += counter["refinement"]
                timings[mapping]["opportunities"] += int(result["status"] == "measured")
            rows.append(row)
        print(json.dumps({"snapshot_seed": int(world.seed),
                          "eligible_opportunities": len(rows) - world_rows_start}), flush=True)
    for summary in timings.values():
        costs = summary.pop("selected_weak_anchor_costs")
        summary["selected_anchor_cost_mean"] = float(np.mean(costs)) if costs else None
        summary["selected_anchor_cost_median"] = float(np.median(costs)) if costs else None
    return {"opportunities": rows}, timings


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--snapshots", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--synthetic-only", action="store_true",
                         help="run the deterministic graph control without loading snapshots")
    args = parser.parse_args()
    manifest = json.loads((args.snapshots / "manifest.json").read_text())
    if manifest.get("protocol") != "RFC0018-G0-snapshot-v1":
        raise ValueError("unexpected snapshot manifest protocol")
    probe = _load_frozen_probe()
    original = probe["_fiedler_vector"]

    def with_mapping(n, i, j, edge_cost, mapping="inverse_cost"):
        counter = {"total": 0, "base": 0, "refinement": 0}
        solve = _mapping_fiedler(mapping, counter)
        return solve(n, i, j, edge_cost)

    probe["_fiedler_vector"] = with_mapping
    synthetic = _synthetic_graphs()
    if args.synthetic_only:
        print(json.dumps(synthetic, indent=2))
        return
    probe["_fiedler_vector"] = original
    snapshots, timings = _snapshot_comparison(args.snapshots, manifest, probe)
    report = {"issue": 850, "protocol": "RFC0018-ISSUE850-SPECTRAL-MAPPING-v1",
              "frozen_probe_source": PINNED_SOURCE,
              "frozen_probe_source_sha256": probe["__source_sha256__"],
              "snapshot_manifest_protocol": manifest["protocol"],
              "snapshot_source_commit": manifest["source_commit"],
              "mapping_definitions": {"inverse_cost": "conductance = 1 / edge_cost (current)",
                                      "direct_cost": "conductance = edge_cost",
                                      "uniform": "conductance = 1 for every shared edge"},
              "synthetic_protocol": {"seed": 850001, "graphs": 100, "nodes_per_graph": 40,
                                     "topology": "cycle with two antipodal planted weak edges",
                                     "background_edge_cost": "uniform draws from [0.9, 1.1)",
                                     "weak_edge_cost": "uniform draws from [0.05, 0.06)",
                                     "partition": "balanced 20/20 split from the second normalized-Laplacian eigenvector",
                                     "solver": "dense symmetric eigensolve; each mapping receives identical graphs",
                                     "determinism": "each graph and mapping solved twice; edge cut compared"},
              "candidate_partition_record_fields": ["anchor_cell_key_i", "anchor_cell_key_j",
                                                     "canonical_partition_sha256", "selected",
                                                     "initial_anchor_edge_cost", "initial_anchor_cost_percentile",
                                                     "usable_candidate"],
              "synthetic": synthetic, "snapshot_summary": timings,
              "snapshot_details": snapshots, "worlds_stepped": False,
              "resistance_score_prior_unchanged": True}
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    print(json.dumps({"synthetic": synthetic, "snapshot_summary": timings,
                      "opportunities": len(snapshots["opportunities"]),
                      "out_json": str(args.out_json)}, indent=2))


if __name__ == "__main__":
    main()

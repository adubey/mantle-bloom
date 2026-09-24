#!/usr/bin/env python3
"""Issue #228 Phase 0a: print a markdown comparison table from characterize_plate_surface.py
results -- one column per world, one row per headline metric. Used to build
analysis/issue228-phase0a/report.md; rerun it against a later representation's results to
compare like for like.

Usage (from repo root):
    python3 bin/debug/summarize_plate_surface.py analysis/issue228-phase0a/results/*.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# (label, section, dotted key, format)
ROWS = (
    ("Elapsed (Myr)", "world", "elapsed_myr", "{:g}"),
    ("Plates", "conservation", "plates", "{:d}"),
    ("Nodes", "conservation", "nodes", "{:,d}"),
    ("Lines", "line_topology", "lines", "{:,d}"),
    ("One-node lines", "line_topology", "one_node_line_fraction", "{:.1%}"),
    ("Rows with overlapping arcs", "line_topology", "rows_with_overlapping_arcs", "{:,d}"),
    ("Max arc stacking depth", "line_topology", "max_arc_stacking_depth", "{:d}"),
    ("Stacked nodes (< 0.5 s, same plate)", "sample_cloud", "stacked_fraction", "{:.2%}"),
    ("Anisotropic nodes (ratio > 2)", "sample_cloud", "anisotropic_fraction", "{:.2%}"),
    ("One-dimensional nodes", "sample_cloud", "one_dimensional_fraction", "{:.2%}"),
    ("Anisotropy row alignment (cos 2α)", "sample_cloud", "anisotropic_row_alignment_mean_cos2", "{:+.2f}"),
    ("Refinement jump p99", "sample_cloud", "refinement_jump.p99", "{:.2f}"),
    ("Boundary nodes", "sample_cloud", "boundary_fraction", "{:.1%}"),
    ("Thin (tendril) nodes", "sample_cloud", "thin_fraction", "{:.2%}"),
    ("Isolated nodes", "sample_cloud", "isolated_nodes", "{:d}"),
    ("Implied cells: collapsed", "implied_cells", "collapsed_fraction", "{:.2%}"),
    ("Implied cells: aspect > 4", "implied_cells", "aspect_over_4_fraction", "{:.3%}"),
    ("Implied cells: skew > 45°", "implied_cells", "skew_over_45_fraction", "{:.3%}"),
    ("Implied cells: folded", "implied_cells", "folded_fraction", "{:.4%}"),
    ("Implied cells: non-conforming", "implied_cells", "nonconforming_fraction", "{:.3%}"),
    ("Uncovered sphere", "coverage", "uncovered_sphere_fraction", "{:.2%}"),
    ("Void sphere (> 1.5 s from nodes)", "coverage", "uncovered_more_than_1_5_spacing_from_any_node_fraction", "{:.3%}"),
    ("Multiply covered sphere", "coverage", "multiply_covered_sphere_fraction", "{:.2%}"),
    ("Nodes inside another plate's outline", "coverage", "nodes_inside_another_plates_outline_fraction", "{:.2%}"),
    ("Plates with holes", "coverage", "plates_with_holes", "{:d}"),
    ("Nominal node area / sphere", "coverage", "nominal_node_area_over_sphere", "{:.3f}"),
    ("Voronoi area / nominal p05", "coverage", "voronoi_area_per_node_over_nominal.p05", "{:.2f}"),
    ("Voronoi area / nominal p95", "coverage", "voronoi_area_per_node_over_nominal.p95", "{:.2f}"),
    ("Land fraction (node)", "conservation", "land_fraction_node", "{:.1%}"),
    ("Sea level (m)", "conservation", "sea_level_m", "{:.0f}"),
    ("Σ Hc, nominal area (10⁹ km³)", "conservation", "extensive_totals.crustal_thickness_m_km3", "e9"),
    ("Σ Hc, Voronoi area (10⁹ km³)", "conservation", "extensive_totals_voronoi_weighted.crustal_thickness_m_km3", "e9"),
    ("Σ Hm, nominal area (10⁹ km³)", "conservation", "extensive_totals.mantle_lithosphere_thickness_m_km3", "e9"),
    ("Σ Hm, Voronoi area (10⁹ km³)", "conservation", "extensive_totals_voronoi_weighted.mantle_lithosphere_thickness_m_km3", "e9"),
    ("Continental Hc (effective), nominal", "conservation", "effective_continental_hc_volume_km3", "e9"),
    ("Continental Hc (effective), Voronoi", "conservation", "effective_continental_hc_volume_voronoi_km3", "e9"),
)


def lookup(doc: dict, section: str, key: str):
    value = doc.get(section, {})
    for part in key.split("."):
        if not isinstance(value, dict) or part not in value:
            return None
        value = value[part]
    return value


def fmt(value, spec: str) -> str:
    if value is None:
        return "–"
    if spec == "e9":
        return f"{value / 1e9:.2f}"
    return spec.format(value)


def main() -> None:
    paths = [Path(p) for p in sys.argv[1:] if not p.endswith(".timings.json")]
    docs = [(p.stem, json.loads(p.read_text())) for p in paths]
    print("| Metric | " + " | ".join(name for name, _ in docs) + " |")
    print("|---|" + "---:|" * len(docs))
    for label, section, key, spec in ROWS:
        print(f"| {label} | " + " | ".join(fmt(lookup(doc, section, key), spec) for _, doc in docs) + " |")


if __name__ == "__main__":
    main()

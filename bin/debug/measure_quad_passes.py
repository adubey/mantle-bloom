"""Issue #228 Phase 4: measure gap filling, volcanism, overlap tracking and lattice quality on
quad vs line worlds.

Steps one world per surface from the same seed. Every `--every` steps it reports:

- **coverage** from a fixed Fibonacci sample of the sphere tested against every plate's own
  `contains_batch`: the fraction covered by no plate and by two or more;
- **overlap tracking**: the fraction of nodes `plates.compute_node_overlap` flags (the
  node-proximity test `merge_split.update_overlap_tracking` stamps), against the fraction whose
  centre actually lies inside another plate (`contains_batch`), and the agreement of the two;
- **gap detection**: the fraction of `gaps._find_gap_points`' sweep it reports uncovered, by
  node distance and as gap filling reads this world (containment on cells);
- **volcanism**: volcano and active-volcano node counts, and eruptions since the last report
  (counted from the vents `volcanism._spread_volcanic_plains` receives);
- **quad lattice quality**: refined (level > 0) cells, cells one cell thin along a row or a
  column, cells with at most one same-plate neighbour, and hole loops.

    cd backend && .venv/bin/python ../bin/debug/measure_quad_passes.py --seed 7 --steps 60
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import gaps, lithosphere, plates as plates_mod, volcanism  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.sparse_quad_patch import PlateWithSparseQuadPatch, unpack_cell_keys  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402


def fibonacci_sphere(n: int) -> np.ndarray:
    k = np.arange(n) + 0.5
    z = 1.0 - 2.0 * k / n
    r = np.sqrt(1.0 - z * z)
    lon = np.pi * (1.0 + 5**0.5) * k
    return np.column_stack([r * np.cos(lon), r * np.sin(lon), z])


def node_areas(world, plate) -> np.ndarray:
    if hasattr(plate, "node_areas_m2"):
        return plate.node_areas_m2()
    return np.full(plate.node_count(), lithosphere.node_area_m2(line_spacing_rad(world.node_density)))


def coverage(world, sample: np.ndarray) -> dict:
    count = np.zeros(len(sample), dtype=int)
    for plate in world.plates:
        if plate.node_count():
            count += plate.contains_batch(sample)
    return {"uncovered": round(float(np.mean(count == 0)), 4), "overlapping": round(float(np.mean(count >= 2)), 4)}


def overlap_tracking(world) -> dict:
    tol = plates_mod.OVERLAP_TOLERANCE_MULT * line_spacing_rad(world.node_density)
    proximity = plates_mod.compute_node_overlap(world.plates, tol)
    flagged = inside = both = total = 0
    for plate in world.plates:
        n = plate.node_count()
        if n == 0:
            continue
        points, _ = plate.all_points_and_elevation()
        exact = np.zeros(n, dtype=bool)
        for other in plate.get_neighbours([p for p in world.plates if p is not plate]):
            exact |= other.contains_batch(points)
        mask = proximity[plate.plate_id]["overlap_mask"]
        flagged += int(mask.sum())
        inside += int(exact.sum())
        both += int((mask & exact).sum())
        total += n
    return {
        "proximity_flagged": round(flagged / total, 4),
        "inside_other": round(inside / total, 4),
        "flagged_and_inside": both,
        "precision": round(both / flagged, 3) if flagged else None,
        "recall": round(both / inside, 3) if inside else None,
    }


def gap_sweep(world) -> dict:
    """The gap-fill sweep's own uncovered fraction, by node distance (the line reading) and as
    `gaps.fill_gaps_by_growing_neighbours` now reads it for this world."""
    context = gaps._existing_node_tree(world)
    spacing = line_spacing_rad(world.node_density)
    live = [p for p in world.plates if p.node_count()]
    swept = sum(len(pts) for _, _, pts in gaps.iter_local_lattice(gaps._GLOBAL_FRAME, spacing_rad=spacing))
    by_distance = len(gaps._find_gap_points(context, spacing))
    try:
        as_filled = len(gaps._find_gap_points(context, spacing, live))
    except TypeError:  # before gap filling learned containment: it only reads node distance
        as_filled = by_distance
    return {"sweep_by_distance": round(by_distance / swept, 4), "sweep_as_filled": round(as_filled / swept, 4)}


def volcano_state(world) -> dict:
    volcano = active = 0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        is_volcano = plate.collect("is_volcano")
        volcano += int(is_volcano.sum())
        active += int((is_volcano & (plate.collect("volcano_active_years_remaining") > 0)).sum())
    return {"volcano_nodes": volcano, "active_volcanoes": active}


def lattice_quality(world) -> dict | None:
    quads = [p for p in world.plates if isinstance(p, PlateWithSparseQuadPatch) and p.node_count()]
    if not quads:
        return None
    cells = refined = thin = spur = holes = 0
    for plate in quads:
        _, level, _, _ = unpack_cell_keys(plate.cell_keys)
        exposed = plate.exposed_sides()  # -v, +u, +v, -u
        graph = plate.adjacency()
        degree = np.diff(graph.offsets)
        cells += plate.node_count()
        refined += int((level > 0).sum())
        thin += int(((exposed[:, 0] & exposed[:, 2]) | (exposed[:, 1] & exposed[:, 3])).sum())
        spur += int((degree <= 1).sum())
        holes += sum(1 for loop in plate._local_boundary_loops() if _signed_area(loop) < 0)
    return {
        "refined_cells": round(refined / cells, 4),
        "one_cell_thin": round(thin / cells, 4),
        "degree_le_1": round(spur / cells, 4),
        "hole_loops": holes,
    }


def _signed_area(loop: np.ndarray) -> float:
    """Sign of a local loop's orientation about its mean direction: + for CCW (outer)."""
    centre = loop.mean(axis=0)
    centre /= np.linalg.norm(centre)
    rel = np.roll(loop, -1, axis=0)
    return float(np.sum(np.einsum("ij,ij->i", np.cross(loop, rel), np.broadcast_to(centre, loop.shape))))


def totals(world) -> dict:
    volume = land = area = 0.0
    nodes = 0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        a = node_areas(world, plate)
        volume += float(np.sum(plate.collect("crustal_thickness_m") * a))
        land += float(np.sum(a[plate.collect("elevation") > world.sea_level_m]))
        area += float(a.sum())
        nodes += plate.node_count()
    return {"volume_1e9_km3": round(volume / 1e18, 4), "land_fraction": round(land / area, 4), "nodes": nodes, "plates": len(world.plates)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=60)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--every", type=int, default=10)
    parser.add_argument("--surfaces", default="quad,lines")
    parser.add_argument("--sample", type=int, default=200_000)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    sample = fibonacci_sphere(args.sample)
    original_spread = volcanism._spread_volcanic_plains
    report = {}
    for surface in args.surfaces.split(","):
        eruptions = [0]

        def counted_spread(plate, world, years, erupted_points, eruptions=eruptions):
            eruptions[0] += sum(len(p) for p in erupted_points)
            return original_spread(plate, world, years, erupted_points)

        volcanism._spread_volcanic_plains = counted_spread
        world = generate_world(seed=args.seed, node_density=args.density, surface=surface)
        rows = []
        step_seconds = 0.0
        for step in range(args.steps + 1):
            if step % args.every == 0:
                row = {"step": step, **totals(world), **coverage(world, sample), **overlap_tracking(world), **gap_sweep(world)}
                row.update(volcano_state(world), eruptions=eruptions[0])
                quality = lattice_quality(world)
                if quality:
                    row.update(quality)
                eruptions[0] = 0
                rows.append(row)
                print(surface, json.dumps(row), flush=True)
            if step < args.steps:
                started = time.perf_counter()
                step_world(world, args.years)
                step_seconds += time.perf_counter() - started
        report[surface] = {"rows": rows, "seconds_per_step": round(step_seconds / max(1, args.steps), 3)}
        print(surface, "seconds/step", report[surface]["seconds_per_step"], flush=True)
    volcanism._spread_volcanic_plains = original_spread

    if args.out:
        args.out.write_text(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

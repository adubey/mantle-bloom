"""Which torque sets each plate's speed: a per-plate breakdown of the quasi-static motion balance (#942).

Usage (repo root, this worktree's venv; follow the compute-budget skill, a few minutes per world):
    backend/.venv/bin/python bin/debug/plate_torque_breakdown.py --out ~/mantle-bloom-captures/issue942/breakdown \
        --seeds 106533960 1097624328 [--years 200] [--every 10] [--world-knobs '{"contact_drag_multiplier": 0.0}']

Generates each world with production defaults and steps it plates-only (``simulate_climate_biomes = False``) at
1 Myr. Every ``--every`` Myr it repeats the motion solve read-only (``torque.solve_plates_quasistatic`` changes no
geometry) and writes, per plate:

- each torque's component along the solved rotation axis, in N m. Positive drives the plate's motion and negative
  resists it; the terms sum to zero. Slab pull, ridge push and the mantle-flow half of basal drag are the loads;
  basal drag, sunk-slab drag and the contact brake (RFC 0008) are the velocity-proportional resistances; what is
  left is the blocked-edge contact force, which the solver does not report per plate.
- how much of the 1 cm/yr motion gate is open: slab pull acts only where a boundary already closes faster than
  1 cm/yr and ridge push only where it opens faster than 1 cm/yr, so the script also records the boundary length
  that moves the right way at any rate, and the slab pull that length would give.
- area-weighted mean speed of continental and oceanic crust, ``|omega| R``, area and continental share.

Writes ``<out>/<seed>.json`` (existing files are kept, so a run resumes).
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).resolve().parents[2] / "backend")]

GATE_M_PER_S_NOTE = "1 cm/yr (torque.BOUNDARY_TRANSFORM_RATE_M_PER_S, boundary.TRANSFORM_RATE_THRESHOLD)"


def _along(axis: np.ndarray, vector: np.ndarray) -> float:
    return float(np.dot(axis, vector))


def sample(world) -> dict:
    from app import lithosphere, torque
    from app.elevation_lines import line_spacing_rad

    diagnostics: dict = {}
    omegas = torque.solve_plates_quasistatic(world, diagnostics=diagnostics)
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
    cache = world.boundary_search_cache
    radius_cm = lithosphere.PLANET_RADIUS_M * 100.0

    local: dict[int, dict] = {}
    side_k: dict[tuple[int, int], dict[int, np.ndarray]] = {}
    for plate in plates:
        neighbours = cache.neighbours(plate, plates, reach)
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
        _convergent, divergent, _transform, contested = torque.classify_boundary_nodes(plate, neighbours, inputs, reach, cache)
        collision = torque.continental_contact_mask(plate, inputs, contested)
        subducting = torque.subducting_boundary_mask(plate, inputs, reach)
        basal_b, basal_k = torque.basal_drag_coefficients(plate, world, spacing)
        slab_k = torque.slab_drag_coefficient_matrix(inputs, subducting, spacing)
        jam = torque.contact_jam_factor(plate, world, neighbours, inputs, reach)
        craton_m = torque.contact_craton_m(plate, neighbours, inputs)
        for neighbour in neighbours:
            pair = (min(plate.plate_id, neighbour.plate_id), max(plate.plate_id, neighbour.plate_id))
            mask = collision & (inputs.neighbor_plate_id == neighbour.plate_id)
            _b, k = torque.collision_drag_coefficients(
                plate, inputs, mask, spacing, craton_m, jam,
                contact_gain=world.contact_drag_multiplier, craton_gain=world.craton_drag_multiplier,
            )
            side_k.setdefault(pair, {})[plate.plate_id] = k

        # The motion gate. `closing` is read from the plates' current (previous step's) omegas, as the solver reads it.
        band = inputs.dist_to_neighbor <= reach
        length = torque.boundary_node_length_m(inputs, spacing)
        down = torque.down_going_crust(plate, inputs, reach)
        closing = torque.closing_rad_per_yr(plate, inputs)
        closing_any = band & down & (closing > 0.0)
        normal = torque.boundary_closing_rate_m_per_s(plate, inputs)
        opening_any = band & ~contested & (normal < 0.0)

        points = inputs.own_points
        areas = np.asarray(plate.accounting_areas_m2(spacing), dtype=float)
        continental = torque.own_is_continental(plate, inputs)
        omega = omegas[plate.plate_id]
        speed = np.linalg.norm(np.cross(omega, points), axis=1) * radius_cm
        local[plate.plate_id] = {
            "omega": omega, "basal_b": basal_b, "basal_k": basal_k, "slab_k": slab_k,
            "slab_pull": torque.slab_pull_torque(plate, inputs, subducting, spacing),
            "slab_pull_ungated": torque.slab_pull_torque(plate, inputs, closing_any, spacing),
            "ridge_push": torque.ridge_push_torque(plate, inputs, divergent, spacing),
            "row": {
                "plate_id": int(plate.plate_id),
                "rate_cm_per_yr": float(np.linalg.norm(omega) * radius_cm),
                "area_km2": float(areas.sum() / 1e6),
                "continental_area_share": float(areas[continental].sum() / areas.sum()),
                "speed_cm_per_yr": {
                    "all": float(np.dot(speed, areas) / areas.sum()),
                    "continental": float(np.dot(speed[continental], areas[continental]) / areas[continental].sum()) if continental.any() else None,
                    "oceanic": float(np.dot(speed[~continental], areas[~continental]) / areas[~continental].sum()) if (~continental).any() else None,
                },
                "continental_area_km2": float(areas[continental].sum() / 1e6),
                "boundary_km": {
                    "band_total": float(length[band].sum() / 1e3),
                    "down_going_closing_any_rate": float(length[closing_any].sum() / 1e3),
                    "down_going_closing_past_gate": float(length[subducting].sum() / 1e3),
                    "opening_any_rate": float(length[opening_any].sum() / 1e3),
                    "opening_past_gate": float(length[divergent].sum() / 1e3),
                    "continental_collision": float(length[collision].sum() / 1e3),
                },
            },
        }

    brake = {plate.plate_id: np.zeros(3) for plate in plates}
    for (a, b), sides in side_k.items():
        if a not in local or b not in local:
            continue
        zero = np.zeros((3, 3))
        pair_k = 0.5 * (sides.get(a, zero) + sides.get(b, zero))
        relative = local[a]["omega"] - local[b]["omega"]
        brake[a] -= pair_k @ relative
        brake[b] += pair_k @ relative

    rows = []
    for plate in plates:
        item = local[plate.plate_id]
        omega = item["omega"]
        rate = float(np.linalg.norm(omega))
        row = item["row"]
        if rate > 0.0:
            axis = omega / rate
            terms = {
                "slab_pull": _along(axis, item["slab_pull"]),
                "ridge_push": _along(axis, item["ridge_push"]),
                "mantle_flow": _along(axis, item["basal_b"]),
                "basal_drag": -_along(axis, item["basal_k"] @ omega),
                "slab_drag": -_along(axis, item["slab_k"] @ omega),
                "contact_brake": _along(axis, brake[plate.plate_id]),
            }
            terms["blocked_edges_residual"] = -sum(terms.values())
            row["torque_along_axis_nm"] = terms
            row["slab_pull_ungated_along_axis_nm"] = _along(axis, item["slab_pull_ungated"])
        row["torque_norm_nm"] = {
            "slab_pull": float(np.linalg.norm(item["slab_pull"])),
            "slab_pull_ungated": float(np.linalg.norm(item["slab_pull_ungated"])),
            "ridge_push": float(np.linalg.norm(item["ridge_push"])),
            "mantle_flow": float(np.linalg.norm(item["basal_b"])),
        }
        rows.append(row)
    return {
        "myr": world.elapsed_years / 1e6,
        "plates": rows,
        "solver": {k: diagnostics[k] for k in ("solver_converged", "blocked_contacts", "yielded_contacts", "residual_torque_norm_nm", "driving_torque_norm_nm")},
    }


def run_world(seed: int, args: argparse.Namespace) -> None:
    target = args.out / f"{seed}.json"
    if target.exists():
        print(seed, "kept", flush=True)
        return
    from app import world as world_module

    started = time.perf_counter()
    world = world_module.generate_world(seed)
    world.simulate_climate_biomes = False
    knobs = json.loads(args.world_knobs)
    for name, value in knobs.items():
        if not hasattr(world, name):
            sys.exit(f"unknown world knob {name}")
        setattr(world, name, value)
    samples = []
    for step in range(1, round(args.years) + 1):
        world_module.step_world(world, 1e6)
        if step % args.every == 0:
            samples.append(sample(world))
    result = {"seed": seed, "years": args.years, "every": args.every, "world_knobs": knobs, "gate": GATE_M_PER_S_NOTE,
              "samples": samples, "wall_s": time.perf_counter() - started}
    target.write_text(json.dumps(result, indent=1) + "\n")
    print(seed, round(result["wall_s"]), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", required=True)
    parser.add_argument("--years", type=float, default=200.0)
    parser.add_argument("--every", type=int, default=10)
    parser.add_argument("--world-knobs", default="{}")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for seed in args.seeds:
        run_world(seed, args)


if __name__ == "__main__":
    main()

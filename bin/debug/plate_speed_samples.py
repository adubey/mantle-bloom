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
           add `--flow-scale 2` to multiply the static mantle-flow amplitude (the #998 flow test), or
           `--ocean-brake 0.25` to extend the contact brake to non-divergent oceanic contacts (in-process patch, #998).
    queue  seeds in series, one process each, skipping finished ones:
           backend/.venv/bin/python bin/debug/plate_speed_samples.py queue --settings ... --setting ... \\
               --seeds-file analysis/issue958/development-seeds.json --at 0 1 2 5 10 --out <captures>/<setting>
"""

import argparse
import json
import sys
import time
from collections import defaultdict

import numpy as np
from sampler_common import boundary_context, build_parser, build_world, plate_boundary, scale_flow, step_and_sample


def record_contacts(world, omegas) -> dict:
    from app import lithosphere, torque

    ctx = boundary_context(world)
    plates = ctx.plates
    scale = lithosphere.PLANET_RADIUS_M * 100.0
    nodes = defaultdict(list)  # (a, b) -> [(km, speed, closing), ...] from either side, half weight each
    for plate in plates:
        _neighbours, inputs, band = plate_boundary(world, ctx, plate)
        if not band.any():
            continue
        length = torque.boundary_node_length_m(inputs, ctx.spacing)
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


def patch_ocean_brake(gain: float) -> None:
    """Exploration of #998 only, no `backend/app` change: extend the contact brake to oceanic contacts in this process.

    The brake brakes continent-continent contested nodes only (`torque.continental_contact_mask`). With this patch
    the quasi-static solve (the only caller inside `torque.py`) also brakes every non-divergent band node
    (convergent, transform or contested, as `classify_boundary_nodes` partitions them) that is not already a
    continent-continent node, with factor 1 and `gain` x CONTACT_RESISTANCE_N_PER_M, the same linear drag on the
    relative velocity. Ridges (divergent nodes) are left free. `rift_damage`'s own use of the mask is untouched."""
    from app import torque

    if gain <= 0:
        return
    original_classify = torque.classify_boundary_nodes
    original_mask = torque.continental_contact_mask
    original_drag = torque.contact_drag_pa_s_per_m
    stash: dict[int, tuple] = {}

    def classify(plate, neighbours, inputs, reach_rad, cache=None):
        result = original_classify(plate, neighbours, inputs, reach_rad, cache)
        stash[id(inputs)] = result
        return result

    def mask(plate, inputs, contested):
        base = original_mask(plate, inputs, contested)
        if sys._getframe(1).f_globals.get("__name__") != "app.torque" or id(inputs) not in stash:
            return base
        convergent, _divergent, transform, contested_ = stash[id(inputs)]
        return base | convergent | transform | contested_

    def drag(plate, inputs, collision_mask, spacing_rad, craton_m=None, jam=None, *, contact_gain=1.0, craton_gain=1.0):
        coefficient = original_drag(plate, inputs, collision_mask, spacing_rad, craton_m, jam, contact_gain=contact_gain, craton_gain=craton_gain)
        ocean = collision_mask & ~(torque.own_is_continental(plate, inputs) & ~inputs.neighbor_is_oceanic)
        if not ocean.any():
            return coefficient
        extra = original_drag(plate, inputs, np.zeros_like(collision_mask), spacing_rad, None, ocean.astype(float), contact_gain=gain)
        return np.maximum(coefficient, extra)

    torque.classify_boundary_nodes = classify
    torque.continental_contact_mask = mask
    torque.contact_drag_pa_s_per_m = drag


def run_one(args: argparse.Namespace) -> None:
    import plate_torque_breakdown as breakdown
    from app import torque

    setting, world = build_world(args)
    # The static mantle-flow field is drawn at generation with a fixed reference rate; scale it after generation so
    # every seed starts from the same world and only the stepped flow torque changes (exploration of #998 only).
    scale_flow(world, args.flow_scale)
    patch_ocean_brake(args.ocean_brake)
    started = time.perf_counter()
    samples = []

    def take(myr: float) -> None:
        omegas = torque.solve_plates_quasistatic(world, diagnostics={})
        samples.append({"myr": myr, **record_contacts(world, omegas), "breakdown": breakdown.sample(world)})
        print(json.dumps({"myr": myr, "plates": len(world.plates), "wall_s": round(time.perf_counter() - started, 1)}), flush=True)

    step_and_sample(world, setting, args, take)
    args.out.write_text(json.dumps({"seed": args.seed, "setting": args.setting, "step_myr": args.step_myr, "flow_scale": args.flow_scale, "ocean_brake": args.ocean_brake,
                                    "samples": samples, "wall_s": time.perf_counter() - started}, default=float) + "\n")


def main() -> None:
    def extra(p: argparse.ArgumentParser) -> None:
        p.add_argument("--ocean-brake", type=float, default=0.0, help="gain (x 5e12 N/m) of a contact brake extended to non-divergent oceanic contacts (default 0, off)")

    args = build_parser(__doc__, __file__, run_one, forward=("flow_scale", "ocean_brake"), extra=extra, default_at=(0, 1, 2, 5, 10)).parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()

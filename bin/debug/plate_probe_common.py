"""Set-up shared by the #995 plate probes: open a world under an RFC 0021 setting, walk every plate's boundary band with the
relative velocity of its neighbour, and step a world while sampling it. Used by `big_plate_rift_probe.py`,
`boundary_geometry_probe.py` and `plate_contact_graph.py`; import with `bin/debug` on `sys.path` and `backend` too."""

import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np


def open_world(settings: Path, setting_name: str, seed: int):
    """The (setting, world) of one generated world; refuses held-out seeds and applies the setting's module knobs."""
    import calibrate_plate_forces as cal
    import run_rift_lines_d1 as d1

    setting = cal.load_setting(settings, setting_name)
    if seed in cal.refused_seeds():
        sys.exit("held-out seed")
    world = d1.make_world(SimpleNamespace(seed=seed, node_density=1.0, world_knobs=setting["world_knobs"],
                                          split_method=setting["split_method"], params=setting["params"]))
    cal.apply_module_knobs(setting["module_knobs"])
    return setting, world


@dataclass
class BoundaryMotion:
    """One plate's boundary-band nodes; `length`, `v_rel` (cm/yr) and `normal` (along `direction_to_neighbor`) are None when
    the plate has no band node."""

    plate: object
    inputs: object
    band: np.ndarray
    length: np.ndarray | None
    v_rel: np.ndarray | None
    normal: np.ndarray | None


def boundary_motion(world, require_neighbour: bool):
    """`(geom, motions)`: `geom` holds spacing, reach, scale (rad/yr x unit vector -> cm/yr), plates and omegas; `motions`
    yields a `BoundaryMotion` per plate. `require_neighbour` drops band nodes with no neighbour plate."""
    from app import lithosphere, torque
    from app.elevation_lines import line_spacing_rad

    omegas = torque.solve_plates_quasistatic(world, diagnostics={})
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
    geom = SimpleNamespace(spacing=spacing, reach=reach, scale=lithosphere.PLANET_RADIUS_M * 100.0, plates=plates, omegas=omegas)

    def walk():
        forms = torque.force_forms(world)
        cache = world.boundary_search_cache
        for plate in plates:
            neighbours = cache.neighbours(plate, plates, reach)
            inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
            torque.apply_force_forms(world, plate, neighbours, inputs, forms)
            band = inputs.dist_to_neighbor <= reach
            if require_neighbour:
                band = band & (inputs.neighbor_plate_id >= 0)
            if not band.any():
                yield BoundaryMotion(plate, inputs, band, None, None, None)
                continue
            neighbour_omega = np.zeros_like(inputs.neighbor_omega)
            for other in plates:
                neighbour_omega[inputs.neighbor_plate_id == other.plate_id] = omegas[other.plate_id]
            v_rel = (np.cross(omegas[plate.plate_id], inputs.own_points) - np.cross(neighbour_omega, inputs.own_points)) * geom.scale
            yield BoundaryMotion(plate, inputs, band, torque.boundary_node_length_m(inputs, spacing), v_rel,
                                 np.sum(v_rel * inputs.direction_to_neighbor, axis=-1))

    return geom, walk()


def step_and_sample(world, setting: dict, args, measure, show):
    """Step plates-only for `args.years` and, every `args.every` Myr, append `{"myr", **measure(world)}` to the samples and
    call `show(sample)`. Returns `(splits, samples)`. `args` needs years, step_myr, every and scratch."""
    import rift_line_recorder as rec
    import run_rift_lines_d1 as d1
    from app import rift_damage

    recorder = rec.Recorder(rift_damage.runtime(world).params.eps_cand) if setting["split_method"] == "rift_lines" else None
    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    samples = []

    def on_step(world, myr: float) -> None:
        if round(myr / args.every) * args.every == myr:
            samples.append({"myr": myr, **measure(world)})
            show(samples[-1])

    splits, *_ = d1.step_all(world, SimpleNamespace(years=args.years, step_myr=args.step_myr, max_steps=None, snapshot_every=0.0), scratch,
                             recorder, on_step)
    return splits, samples

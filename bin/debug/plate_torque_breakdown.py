"""Which torque sets each plate's speed: a per-plate breakdown of the quasi-static motion balance (#942).

Usage (repo root, this worktree's venv; follow the compute-budget skill, a few minutes per world):
    backend/.venv/bin/python bin/debug/plate_torque_breakdown.py --out ~/mantle-bloom-captures/issue942/breakdown \
        --seeds 106533960 1097624328 [--years 200] [--every 10] [--world-knobs '{"contact_drag_multiplier": 0.0}']

Generates each world with production defaults and steps it plates-only (``simulate_climate_biomes = False``) at
1 Myr. Every ``--every`` Myr it repeats the motion solve read-only (``torque.solve_plates_quasistatic`` changes no
geometry) and writes, per plate:

- each torque's component along the solved rotation axis, in N m. Positive drives the plate's motion and negative
  resists it. Slab pull, ridge push and the mantle-flow half of basal drag are the loads; basal drag, sunk-slab
  drag and the contact brake (RFC 0008) are the velocity-proportional resistances. The blocked-edge contact force,
  which the solver does not report per plate, is taken as what is left, so these terms sum to zero by
  construction and that sum checks nothing. The independent check is the solver's own residual over its driving
  torque (``solver.residual_torque_norm_nm / driving_torque_norm_nm``, every plate and contact at once), which
  ``summarize_torque_breakdown.py`` reports.
- how much of the 1 cm/yr motion gate is open: slab pull acts only where a boundary already closes faster than
  1 cm/yr and ridge push only where it opens faster than 1 cm/yr, so the script also records the boundary length
  that moves the right way at any rate, and the slab pull that length would give.
- area-weighted mean speed of continental and oceanic crust, ``|omega| R``, area and continental share.
- two speeds the plate would have alone, from the same snapshot: carried by mantle flow only, and under all of its
  own loads and drags with no contact; and the mean speed of the prescribed mantle flow under the plate. Both
  are least-squares solves, so a plate whose drag matrix is near-singular contributes the minimum-norm rotation,
  which can understate its speed; read pooled means of these as estimates.
- A plate with zero rotation has no axis: it enters the speed and area pools but not the torque shares.
- over its subducting nodes, the work slab pull does and the work sunk-slab drag takes, split by trench-normal and
  trench-parallel motion, how much of the pull cancels across the plate, and how much trench is dragged unpulled.
- each plate's area-weighted speed in the no-net-rotation frame, and the sample's net rotation (#958).
- RFC 0021's G0 measures (#958, `_g0_terms`): the boundary length that is slab-attached under a test that does not
  read the gate (T2), the area-weighted basal traction in Pa (T4), the ridge push per metre of front at the nodes
  that push (T5), and the oldest oceanic floor on the plate.

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


def _crust_speeds(omega: np.ndarray, points: np.ndarray, areas: np.ndarray, continental: np.ndarray, radius_cm: float) -> dict:
    """`|omega| R` and the area-weighted mean speed of all, continental and oceanic crust under `omega`, in cm/yr."""
    speed = np.linalg.norm(np.cross(omega, points), axis=1) * radius_cm

    def mean(mask: np.ndarray) -> float | None:
        return float(np.dot(speed[mask], areas[mask]) / areas[mask].sum()) if mask.any() else None

    return {"rate": float(np.linalg.norm(omega) * radius_cm), "all": mean(np.ones(len(points), dtype=bool)),
            "continental": mean(continental), "oceanic": mean(~continental)}


def _solve(k: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.linalg.lstsq(k, b, rcond=None)[0]


def _slab_terms(plate, inputs, subducting: np.ndarray, spacing: float, omega: np.ndarray) -> dict:
    """Why slab pull and sunk-slab drag net out as they do on one plate, over its subducting nodes.

    Powers are in W at the solved `omega`: force times the node's own velocity. Slab pull acts along the direction
    to the trench, so it does work only on the trench-normal part of the motion, positive where the plate moves
    toward its trench and negative where it moves away (the neighbour overrides it). Sunk-slab drag acts against
    the whole velocity, so its dissipation splits into a trench-normal and a trench-parallel part; under RFC 0021's
    B.2 (`slab_drag_form = "trench_normal"`) it acts on the normal part alone and the parallel part is zero.
    """
    from app import lithosphere, torque

    if not subducting.any():
        return {}
    length = torque.boundary_node_length_m(inputs, spacing)[subducting]
    pull = torque.slab_pull_force_n(plate, inputs, subducting, spacing)
    toward = inputs.direction_to_neighbor[subducting]
    points = inputs.own_points[subducting]
    velocity = np.cross(omega, points) * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    normal = np.einsum("ni,ni->n", velocity, toward)
    # The speed the drag resists: under B.2 along the same tangent normal the solve uses, else along the chord.
    drag_normal = np.einsum("ni,ni->n", velocity, torque.trench_normal(inputs, subducting)) if inputs.forms.slab_drag == "trench_normal" else normal
    speed_sq = np.einsum("ni,ni->n", velocity, velocity)
    drag_per_speed = torque.SLAB_MANTLE_VISCOSITY_PA_S / torque.SHEAR_ZONE_THICKNESS_M * torque.SUBDUCTION_LSINK_M * length
    drag_per_speed = drag_per_speed * torque.slab_drag_nodes(plate, inputs, subducting)[subducting]  # B.2: pulled nodes only
    pull_power = pull * normal
    moments = np.cross(points * lithosphere.PLANET_RADIUS_M, pull[:, None] * toward)
    moment_sum = float(np.linalg.norm(moments, axis=1).sum())
    return {
        "subducting_km": float(length.sum() / 1e3),
        "no_pull_km": float(length[pull <= 0.0].sum() / 1e3),  # buoyant columns: dragged but not pulled
        "moving_away_from_trench_km": float(length[normal < 0.0].sum() / 1e3),
        # Trench-normal speed already past pull / drag: there the drag exceeds the pull on the normal motion alone.
        "past_terminal_speed_km": float(length[drag_normal * drag_per_speed > pull].sum() / 1e3),
        "pull_n": float(pull.sum()),
        "hm_km_length_weighted": float(np.dot(inputs.column_hm[subducting], length) / length.sum() / 1e3),
        "pull_power_toward_trench_w": float(pull_power[normal > 0.0].sum()),
        "pull_power_away_from_trench_w": float(pull_power[normal < 0.0].sum()),
        "drag_power_trench_normal_w": float(np.dot(drag_per_speed, drag_normal**2)),
        "drag_power_trench_parallel_w": 0.0 if inputs.forms.slab_drag == "trench_normal"
        else float(np.dot(drag_per_speed, np.clip(speed_sq - drag_normal**2, 0.0, None))),
        "pull_torque_surviving_cancellation": float(np.linalg.norm(moments.sum(axis=0)) / moment_sum) if moment_sum else None,
        "pull_moment_sum_nm": moment_sum,
    }


def _pulled_km(plate, inputs, mask: np.ndarray, length: np.ndarray) -> float:
    """Length of `mask` that actually pulls, in km: only negatively buoyant nodes count, each weighted by the pull
    share the subduction test gives it (1 but under B.4a). A buoyant node is in the test's mask but pulls nothing."""
    from app import torque

    mask = torque.pulling_nodes(plate, inputs, mask)
    weight = torque.subduction_pull_weight(plate, inputs, mask)
    return float((length[mask] if weight is None else length[mask] * weight).sum() / 1e3)


def _oceanic_age_myr(plate, world) -> np.ndarray:
    """Per node, the age of oceanic floor in Myr from `node_created_years`; NaN where undated (the field's sentinel)."""
    from app.surface_fields import SURFACE_FIELDS

    created = plate.collect("node_created_years")
    undated = created == SURFACE_FIELDS["node_created_years"].sentinel
    return np.where(undated, np.nan, (world.elapsed_years - created) / 1e6)


def geometric_slab_mask(plate, inputs, reach: float, ages_by_plate: dict[int, np.ndarray], own_age: np.ndarray) -> np.ndarray:
    """RFC 0021 B.4b's subduction test, read for T2 without the 1 cm/yr gate (#958): down-going crust in the boundary
    band at a contact that isn't diverging (closing >= -1 cm/yr, `jam_mask`'s contact test), and, where both sides
    are oceanic, only the older side: denser floor goes under. Continental fronts are decided by `down_going_crust`
    (only oceanic or slab-attached continental crust goes down), so the age tie-break is only for ocean against
    ocean. An undated side counts as older, so an undated contact keeps both sides, as the gate does today."""
    from app import torque

    contact = torque._jam_contact(plate, inputs, reach) & torque.down_going_crust(plate, inputs, reach)
    if not contact.any():
        return contact
    neighbour_age = np.full(len(own_age), np.nan)
    for plate_id, ages in ages_by_plate.items():
        at = contact & (inputs.neighbor_plate_id == plate_id) & (inputs.neighbor_node_index >= 0)
        neighbour_age[at] = ages[inputs.neighbor_node_index[at]]
    younger = inputs.neighbor_is_oceanic & (own_age < neighbour_age)  # NaN compares False: undated stays in
    return contact & ~younger


def _consumed_area_km2(plate, world) -> float:
    """Oceanic floor this plate lost under a neighbour in the last `EVIDENCE_LOOKBACK_YEARS` (10 Myr), from
    `collision_polarity`'s consumption evidence: one bin per place, so this is the area of distinct bins, an upper
    bound on the consumed area at the bin's resolution. A plate with no slab has none; T2's candidate definitions
    read it as the plate having a slab of its own, built by past convergence rather than its current speed (#958)."""
    from app import collision_polarity as polarity, lithosphere

    store = (getattr(world, "collision_evidence", None) or {}).get(plate.plate_id)
    if store is None or not len(store):
        return 0.0
    recent = (store.source == polarity.SOURCE_CONSUMPTION) & (store.role == polarity.ROLE_LOWER)
    recent &= store.last_seen >= world.elapsed_years - polarity.EVIDENCE_LOOKBACK_YEARS
    return float(len(np.unique(store.keys[recent])) * (store.bin_rad * lithosphere.PLANET_RADIUS_M / 1e3) ** 2)


def _g0_terms(plate, world, inputs, masks: dict, spacing: float, omega: np.ndarray, areas: np.ndarray) -> dict:
    """RFC 0021 G0's per-plate measures (#958): T2's slab-attached length, T4's basal traction, T5's ridge push."""
    from app import lithosphere, mantle, torque

    length = torque.boundary_node_length_m(inputs, spacing)
    points = inputs.own_points
    relative = (mantle.flow_at(points, world.mantle_centers) - np.cross(omega, points)) * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    traction = torque.ASTHENOSPHERE_VISCOSITY_PA_S / torque.SHEAR_ZONE_THICKNESS_M * np.linalg.norm(relative, axis=1)
    divergent = masks["divergent"]
    push = np.linalg.norm(torque.ridge_push_force_vectors_n(plate, inputs, divergent, spacing), axis=1)[divergent]
    front = length[divergent]
    age = masks["oceanic_age_myr"]
    oceanic = ~masks["continental"] & np.isfinite(age)
    return {
        "geometric_slab_km": float(length[masks["geometric_slab"]].sum() / 1e3),
        "consumed_10myr_km2": _consumed_area_km2(plate, world),
        "basal_traction_pa": float(np.dot(traction, areas) / areas.sum()),
        "ridge_push_km": float(front.sum() / 1e3),
        "ridge_push_n_per_m": float(push.sum() / front.sum()) if front.sum() else None,
        "oceanic_age_max_myr": float(age[oceanic].max()) if oceanic.any() else None,
        "has_subducting_node": bool(masks["subducting"].any()),
    }


def _net_rotation(local: dict) -> np.ndarray:
    """The lithosphere's area-weighted net rotation, `(sum A (I - p pT))^-1 sum A (I - p pT) omega`, in rad/yr: the
    rotation that no-net-rotation speeds remove, as published plate models' NNR frame does (RFC 0021 T1, #958)."""
    k_total, moment = np.zeros((3, 3)), np.zeros(3)
    for item in local.values():
        points, areas = item["points"], item["areas"]
        k = areas.sum() * np.eye(3) - np.einsum("n,ni,nj->ij", areas, points, points)
        k_total += k
        moment += k @ item["omega"]
    return np.linalg.solve(k_total, moment) if np.linalg.det(k_total) > 0.0 else np.zeros(3)


def sample(world) -> dict:
    from app import lithosphere, mantle, torque
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
    ages = {plate.plate_id: _oceanic_age_myr(plate, world) for plate in plates}
    forms = torque.force_forms(world)
    for plate in plates:
        neighbours = cache.neighbours(plate, plates, reach)
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
        torque.apply_force_forms(world, plate, neighbours, inputs, forms)  # as the solve does (RFC 0021, #971)
        _convergent, divergent, _transform, contested = torque.classify_boundary_nodes(plate, neighbours, inputs, reach, cache)
        collision = torque.continental_contact_mask(plate, inputs, contested)
        subducting = torque.subducting_boundary_mask(plate, inputs, reach)
        basal_b, basal_k = torque.basal_drag_coefficients(plate, world, spacing)
        slab_k = torque.slab_drag_coefficient_matrix(inputs, torque.slab_drag_nodes(plate, inputs, subducting), spacing)
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
        solved = _crust_speeds(omega, points, areas, continental, radius_cm)
        g0 = _g0_terms(plate, world, inputs, {
            "divergent": divergent, "subducting": subducting, "continental": continental, "oceanic_age_myr": ages[plate.plate_id],
            "geometric_slab": geometric_slab_mask(plate, inputs, reach, ages, ages[plate.plate_id]),
        }, spacing, omega, areas)
        slab_pull = torque.slab_pull_torque(plate, inputs, subducting, spacing)
        ridge_push = torque.ridge_push_torque(plate, inputs, divergent, spacing)
        flow_speed = np.linalg.norm(mantle.flow_at(points, world.mantle_centers), axis=1) * radius_cm
        local[plate.plate_id] = {
            "points": points, "areas": areas,
            "omega": omega, "basal_b": basal_b, "basal_k": basal_k, "slab_k": slab_k,
            "slab_pull": slab_pull,
            "slab_pull_ungated": torque.slab_pull_torque(plate, inputs, closing_any, spacing),
            "ridge_push": ridge_push,
            "row": {
                "plate_id": int(plate.plate_id),
                "rate_cm_per_yr": solved["rate"],
                "area_km2": float(areas.sum() / 1e6),
                "continental_area_share": float(areas[continental].sum() / areas.sum()),
                "speed_cm_per_yr": {kind: solved[kind] for kind in ("all", "continental", "oceanic")},
                # Speeds this plate would have alone, from this snapshot's own terms: under the mantle-flow
                # coupling only (K_basal omega = b), and under all of its own loads and drags with no contact
                # (brake or blocked edge) and its neighbours' motion ignored.
                "counterfactual_speed_cm_per_yr": {
                    "mantle_flow_only": _crust_speeds(_solve(basal_k, basal_b), points, areas, continental, radius_cm),
                    "own_loads_no_contact": _crust_speeds(
                        _solve(basal_k + slab_k, basal_b + slab_pull + ridge_push), points, areas, continental, radius_cm
                    ),
                },
                "mantle_flow_under_plate_cm_per_yr": float(np.dot(flow_speed, areas) / areas.sum()),
                "slab": _slab_terms(plate, inputs, subducting, spacing, omega),
                "continental_area_km2": float(areas[continental].sum() / 1e6),
                "g0": g0,
                "boundary_km": {
                    "band_total": float(length[band].sum() / 1e3),
                    "down_going_closing_any_rate": float(length[closing_any].sum() / 1e3),
                    "down_going_closing_past_gate": float(length[subducting & closing_any].sum() / 1e3),  # inside the closing set it is divided by (B.4b also admits stationary contacts)
                    # RFC 0021 G1's B.4 measure (#971): closing down-going length that pulls, at B.4a's ramp weight.
                    "down_going_closing_pulled": float(_pulled_km(plate, inputs, subducting & closing_any, length)),
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

    net = _net_rotation(local)
    rows = []
    for plate in plates:
        item = local[plate.plate_id]
        item["row"]["speed_nnr_cm_per_yr"] = _crust_speeds(item["omega"] - net, item["points"], item["areas"],
                                                           np.ones(len(item["points"]), dtype=bool), radius_cm)["all"]
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
        "net_rotation_cm_per_yr": float(np.linalg.norm(net) * radius_cm),
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

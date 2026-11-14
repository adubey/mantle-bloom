"""Plate-level glue for the RFC 0018 damage-band prototype (D1, #888): the per-step damage pass for every plate,
the breakthrough decision and mask split, and the hooks that run them inside ``world.step_world``.

See ``rift_damage_mechanism.py`` for the kernels and the list of readings of the RFC. Not production code.

``install(params)`` replaces ``merge_split.maybe_split_plate`` (the damage pass for every plate runs once per
step on its first call, then only the breakthrough decision runs, gated by ``SPLIT_MIN_AGE_YEARS`` and
``SPLIT_MIN_NODES`` as in production) and wraps ``quad_tectonics._open_rift`` so margin insertions inherit
their stretched share of damage. ``uninstall()`` restores both.
"""

import time
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra

import replay_gate2_capacity as gate2
import rift_damage_mechanism as mech
import rift_damage_prototype as proto
from rift_damage_mechanism import Params, RADIUS_M

from app import cratons, elevation_lines, faults, lithosphere, mantle, merge_split, quad_tectonics, rheology, terranes, thin_shell  # noqa: E501
from app.elevation_lines import ELEV_CHANGE_RIFT, line_spacing_rad


@dataclass
class PlateStep:
    """What a plate's damage pass left behind: arrays for the recorder and the breakthrough decision."""

    plate_id: int
    keys: np.ndarray
    points: np.ndarray
    areas: np.ndarray
    graph: csr_matrix
    tensor: np.ndarray  # frozen load tensor (N, 3)
    n_eff: np.ndarray  # strength at the end of the step
    eps: np.ndarray
    d_eps: np.ndarray
    broken: np.ndarray
    oceanic: np.ndarray
    weak: np.ndarray  # inherited weakness marks (belt, volcano, fault) per cell
    plumed: np.ndarray  # plume heat present per cell
    elevation: np.ndarray
    hc: np.ndarray
    seconds: float = 0.0
    stats: dict = field(default_factory=dict)


@dataclass
class Mechanism:
    params: Params
    steps: dict = field(default_factory=dict)  # plate_id -> PlateStep of the latest pass
    events: list = field(default_factory=list)  # one dict per split / counter, with world time
    counters: dict = field(default_factory=lambda: {"withheld": 0, "capped": 0, "path_overlap": 0, "ligament_pass": 0,
                                                    "split_rejected": 0, "no_shell_mesh": 0})
    pass_key: float | None = None
    prev_reference: dict = field(default_factory=dict)
    seconds_damage: float = 0.0
    seconds_split: float = 0.0
    ledger_gap: dict = field(default_factory=dict)  # worst relative closure error seen per quantity


_STATE: Mechanism | None = None
_ORIGINAL: dict = {}


# --- the damage pass -----------------------------------------------------------------------------


def force_inputs(world, plate):
    """The per-node force field and the divergent-boundary mask (D0's event force field)."""
    return gate2._event_force_field(world, plate, with_divergent=True)


def fault_mask(world, plate, points, spacing_m) -> np.ndarray:
    """Cells within one cell spacing of a fault trace of this plate."""
    traces = [faults.fault_world_points(f, plate) for f in world.faults if f.plate_id == plate.plate_id]
    if not traces:
        return np.zeros(len(points), dtype=bool)
    from scipy.spatial import cKDTree

    tree = cKDTree(np.concatenate(traces))
    distance, _ = tree.query(points)
    return distance * RADIUS_M <= spacing_m


def _is_reference(hc, hm) -> np.ndarray:
    return ((hc == lithosphere.REFERENCE_HC_CONTINENTAL_M) & (hm == lithosphere.REFERENCE_HM_CONTINENTAL_M)) | (
        (hc == lithosphere.REFERENCE_HC_OCEANIC_M) & (hm == lithosphere.YOUNG_RIDGE_HM_M))


def reset_after_eruption(state: Mechanism, plate, eps, warm, broken) -> None:
    """Zero the damage of cells whose column became exactly a reference column since the last pass: that is
    what ``_erupt_melted_nodes`` leaves (new rock carries no inherited damage). In place."""
    keys, now = plate.cell_keys, _is_reference(plate.collect("crustal_thickness_m"), plate.collect("mantle_lithosphere_thickness_m"))
    previous = state.prev_reference.get(plate.plate_id)
    if previous is None:
        return
    old_keys, old_ref = previous
    at = np.clip(np.searchsorted(old_keys, keys), 0, len(old_keys) - 1)
    known = old_keys[at] == keys
    reset = now & known & ~old_ref[at]
    eps[reset], warm[reset], broken[reset] = 0.0, 0.0, False


def damage_step(world, plate, years: float, state: Mechanism) -> PlateStep | None:
    """Steps 1-4 for one plate: load, strength, strain, thinning with flank transfer, the broken latch."""
    p = state.params
    n = plate.node_count()
    if n == 0 or years <= 0.0:
        return None
    started = time.perf_counter()
    points, _ = plate.all_points_and_elevation()
    spacing_m = line_spacing_rad(world.node_density) * RADIUS_M
    eps0 = plate.collect("rift_strain")
    warm0 = plate.collect("rift_warm_strain")
    broken = plate.collect("rift_broken").astype(bool)
    reset_after_eruption(state, plate, eps0, warm0, broken)
    hc, hm = plate.collect("crustal_thickness_m"), plate.collect("mantle_lithosphere_thickness_m")
    codes, lag = plate.collect("crust_type_code"), plate.collect("moho_thermal_lag_c")
    heat = plate.collect("plume_heat_myr")
    continental = np.asarray(elevation_lines.effective_is_continental_from_codes(codes, plate.crust_type == "continental"), bool)

    field_, divergent = force_inputs(world, plate)
    areas = field_.areas_m2
    tensor = p.k_cut * proto.cut_tensor(points, areas, field_.total_n, spacing_m)
    lam, _, _ = proto.principal(tensor)
    heated = heat > 0.0
    frac_heated = np.where(heated, np.clip(heat / (years / mech.YEAR_PER_MYR), 0.0, 1.0), 0.0)
    load_bg = np.maximum(lam, 0.0)
    load_plume = np.maximum(lam + p.n_plume * heated, 0.0)

    kind = np.where(continental, 0, np.where(divergent, 1, 2))
    ncap = mech.capacity(kind, hc, hm, lag) * p.k_cap * np.where(continental, 1.0, p.oceanic_k_cap_factor)
    belt = cratons.belt_mask(plate, world.elapsed_years)
    volcano = plate.collect("is_volcano").astype(bool)
    near_fault = fault_mask(world, plate, points, spacing_m)
    base_plume = ncap * mech.inherited_multiplier(belt, volcano, heated, near_fault)
    base_quiet = ncap * mech.inherited_multiplier(belt, volcano, np.zeros(n, bool), near_fault)

    disk = mech.disk_average_matrix(points, areas, p.w0_m / 2.0)
    self_weight = np.asarray(disk.diagonal()) if p.live_self_term else np.zeros(n)
    nb_eps = disk @ eps0 - self_weight * eps0
    nb_warm = disk @ warm0 - self_weight * warm0

    floor = rheology.RIFT_CRITICAL_THICKNESS_M
    to_floor = np.log(np.maximum(hc, floor) / floor)
    latch = np.minimum(p.eps_break, eps0 + to_floor)
    if p.oceanic_latch_on_break_only:
        latch = np.where(continental, latch, p.eps_break)
    latch = np.where(broken, eps0, latch)

    w0 = mech.softening(disk @ eps0, disk @ warm0, p)
    loaded = (load_bg > base_quiet * w0) | (heated & (load_plume > base_plume * w0))
    first = mech.integrate_strain(load_plume, base_plume, eps0, warm0, nb_eps, nb_warm, self_weight, frac_heated * years, latch, p)
    second = mech.integrate_strain(load_bg, base_quiet, first.eps, first.warm, nb_eps, nb_warm, self_weight, (1.0 - frac_heated) * years, latch, p)
    eps1, warm1 = second.eps, np.minimum(second.warm, second.eps)
    d_eps = np.maximum(eps1 - eps0, 0.0)
    newly_broken = (eps1 >= latch - 1e-9) & loaded & ~broken
    state.counters["capped"] += int(first.capped.sum() + second.capped.sum())

    graph, _, _ = mech.plate_graph(plate, points)
    adjacency = graph
    columns = {
        "hc": hc, "hm": hm, "material": plate.collect("continental_material_m"), "craton": plate.collect("craton_crust_m"),
        "cover": plate.collect("mobile_cover_m"), "cover_material": plate.collect("mobile_cover_continental_m"),
    }
    plateau, arc, moment = terranes.plate_tracers(plate)
    columns.update(plateau=plateau, arc=arc, moment=moment)
    before_totals = _volumes(columns, areas)
    out = mech.thin_with_flank_transfer(areas, adjacency, d_eps, broken | newly_broken, columns, p)
    after_totals = _volumes({**columns, **{k: getattr(out, k) for k in columns}}, areas)
    _record_closure(state, before_totals, after_totals, out.craton_lost_m3, areas, columns["craton"])
    state.counters["withheld"] += int(out.withheld.sum())

    density = lithosphere.node_crust_density(codes, plate.crust_type)
    elevation = plate.collect("elevation")
    changed = (out.hc != hc) | (out.hm != hm)
    shift = lithosphere.isostatic_elevation(out.hc, out.hm, density) - lithosphere.isostatic_elevation(hc, hm, density)
    new_elevation = np.where(changed, rheology.clip_elevation_bounds(elevation + shift), elevation)
    reason = plate.collect("elev_change_reason")
    reason[(d_eps > 0.0) & ~out.withheld & changed] = ELEV_CHANGE_RIFT
    plate.set_fields_on_plate(
        crustal_thickness_m=out.hc, mantle_lithosphere_thickness_m=out.hm, continental_material_m=out.material,
        craton_crust_m=out.craton, mobile_cover_m=out.cover, mobile_cover_continental_m=np.minimum(out.cover_material, out.material),
        elevation=new_elevation, elev_change_reason=reason, rift_strain=eps1, rift_warm_strain=warm1,
        rift_broken=broken | newly_broken,
    )
    if "arc_crust_m" in mech.SURFACE_FIELDS:
        plate.set_fields_on_plate(arc_crust_m=out.arc)
    plate.set_fields_on_plate(plateau_crust_m=out.plateau, buoyant_crust_time_m_yr=out.moment)
    if out.craton_lost_m3 > 0.0:
        cratons.record(world, "rifted_stretch_m3", out.craton_lost_m3)

    eps_bar = disk @ eps1
    warm_bar = disk @ warm1
    n_eff = base_quiet * mech.softening(eps_bar, warm_bar, p)
    state.prev_reference[plate.plate_id] = (plate.cell_keys.copy(), _is_reference(out.hc, out.hm))
    step = PlateStep(plate.plate_id, plate.cell_keys.copy(), points, areas, graph, tensor, n_eff, eps1, d_eps,
                     broken | newly_broken, ~continental, belt | volcano | near_fault, heated, new_elevation, out.hc, time.perf_counter() - started)
    state.steps[plate.plate_id] = step
    return step


def _volumes(columns, areas) -> dict:
    return {name: float(np.dot(columns[name], areas)) for name in
            ("hc", "hm", "material", "cover", "cover_material", "plateau", "arc", "moment", "craton")}


def _record_closure(state, before, after, craton_lost, areas, craton_before) -> None:
    """Interior thinning moves volume; it must conserve every total. Craton is booked to its sink."""
    for name, value in before.items():
        scale = max(abs(value), 1.0)
        expected = value - (craton_lost if name == "craton" else 0.0)
        gap = abs(after[name] - expected) / scale
        state.ledger_gap[name] = max(state.ledger_gap.get(name, 0.0), gap)


# --- breakthrough --------------------------------------------------------------------------------


def _path_normals(points, path, first, second):
    """In-plane unit normals (a, b in each cell's tangent basis) and path lengths ds for an ordered chain."""
    prev_i, next_i = np.r_[path[0], path[:-1]], np.r_[path[1:], path[-1]]
    tangent = points[next_i] - points[prev_i]
    ta = np.sum(tangent * first[path], axis=1)
    tb = np.sum(tangent * second[path], axis=1)
    norm = np.maximum(np.hypot(ta, tb), 1e-15)
    normal = np.column_stack([-tb / norm, ta / norm])
    step = RADIUS_M * np.arccos(np.clip(np.sum(points[path[1:]] * points[path[:-1]], axis=1), -1.0, 1.0))
    ds = np.zeros(len(path))
    ds[:-1] += 0.5 * step
    ds[1:] += 0.5 * step
    ds[0] += 0.5 * step[0]
    ds[-1] += 0.5 * step[-1]
    return normal, ds


def _trace(pred, node):
    chain = [int(node)]
    while pred[chain[-1]] >= 0:
        chain.append(int(pred[chain[-1]]))
    return chain


def band_path(step: PlateStep, boundary: np.ndarray, band: np.ndarray):
    """Cheapest edge-connected chain from the plate boundary through ``band`` (broken cells; one component) back
    to the boundary on its far side: two multi-source Dijkstra searches from the boundary cells at each end of
    the band's long axis, joined at the band cell minimizing the sum. Broken cells cost zero and intact cells
    cost ``N_eff × length``. Returns the ordered cell chain, or None when the band has no ends to reach."""
    n = len(step.points)
    cost = np.where(step.broken, 0.0, step.n_eff)
    coo = step.graph.tocoo()
    weight = np.maximum(0.5 * (cost[coo.row] + cost[coo.col]) * coo.data, 1e-9 * coo.data)
    graph = csr_matrix((weight, (coo.row, coo.col)), shape=(n, n))
    cells = np.flatnonzero(band)
    centre = step.points[cells].mean(axis=0)
    _, _, vt = np.linalg.svd(step.points[cells] - centre, full_matrices=False)
    coord = (step.points - centre) @ vt[0]
    lo, hi = np.percentile(coord[cells], [10, 90])
    minus, plus = np.flatnonzero(boundary & (coord <= lo)), np.flatnonzero(boundary & (coord >= hi))
    if not len(minus) or not len(plus):
        return None
    d_minus, p_minus, _ = dijkstra(graph, directed=False, indices=minus, min_only=True, return_predecessors=True)
    d_plus, p_plus, _ = dijkstra(graph, directed=False, indices=plus, min_only=True, return_predecessors=True)
    total = np.where(band, d_minus + d_plus, np.inf)
    meet = int(np.argmin(total))
    if not np.isfinite(total[meet]):
        return None
    chain = _trace(p_minus, meet)[::-1] + _trace(p_plus, meet)[1:]
    return chain


def split_masks(world, plate, step: PlateStep, path: np.ndarray, min_nodes: int, why: dict | None = None):
    """Remove the path; require exactly two components; give each path cell to the component it shares the most
    edges with, ties to the side whose fitted Euler pole matches the mantle flow there better. Returns the
    boolean mask of side A (path cells assigned), or None if either daughter is disconnected or too small
    (the reason is counted in ``why``)."""
    def reject(reason):
        if why is not None:
            why[reason] = why.get(reason, 0) + 1
    n = len(step.points)
    off_path = np.ones(n, dtype=bool)
    off_path[path] = False
    sub = step.graph[off_path][:, off_path]
    count, labels = connected_components(sub, directed=False)
    if count != 2:
        reject(f"components_after_removal_{min(count, 5)}")
        return None
    side = np.full(n, -1)
    side[off_path] = labels
    pattern = (step.graph > 0).astype(np.int8).tocsr()
    for cell in path:
        row = pattern.indices[pattern.indptr[cell]:pattern.indptr[cell + 1]]
        votes = np.bincount(side[row][side[row] >= 0], minlength=2)
        if votes[0] != votes[1]:
            side[cell] = int(np.argmax(votes))
        else:
            side[cell] = _tie_side(world, step, side, cell)
    mask_a = side == 0
    for mask in (mask_a, ~mask_a):
        if int(mask.sum()) < min_nodes:
            reject("daughter_too_small")
            return None
        if connected_components(step.graph[mask][:, mask], directed=False)[0] != 1:
            reject("daughter_disconnected")
            return None
    return mask_a


def _tie_side(world, step, side, cell) -> int:
    flow = mantle.flow_at(step.points, world.mantle_centers)
    best, score = 0, -np.inf
    for k in (0, 1):
        sel = (side == k)
        pole = mantle.fit_euler_pole(step.points[sel], flow[sel])
        velocity = np.cross(pole, step.points[cell])
        s = float(velocity @ flow[cell])
        if s > score:
            best, score = k, s
    return best


def random_path_weak_share(step: PlateStep, boundary: np.ndarray, plate_id: int, elapsed_years: float, trials: int = 20) -> float:
    """Mean inherited-weakness share along matched random edge-to-edge paths: shortest paths between random pairs of
    boundary cells at least 60% of the plate's longest boundary-to-boundary span apart (the falsification test's
    baseline). Deterministic in (plate, time)."""
    rng = np.random.default_rng((plate_id, round(elapsed_years)))
    cells = np.flatnonzero(boundary)
    if len(cells) < 4:
        return float("nan")
    span = 0.6 * float(np.max(np.linalg.norm(step.points[cells][:, None] - step.points[cells][None, :], axis=2)))
    shares = []
    for _ in range(trials * 5):
        a, b = rng.choice(cells, 2, replace=False)
        if np.linalg.norm(step.points[a] - step.points[b]) < span:
            continue
        _, pred = dijkstra(step.graph, directed=False, indices=a, return_predecessors=True)
        chain = _trace(pred, b)
        shares.append(float(step.weak[chain].mean()))
        if len(shares) == trials:
            break
    return float(np.mean(shares)) if shares else float("nan")


def localization(step: PlateStep, p: Params) -> dict:
    """The #699 failure mode at a split: components of damaged cells above 1% of plate area, and damaged share."""
    damaged = step.eps >= p.eps_cand
    total = float(step.areas.sum())
    share = float(step.areas[damaged].sum() / total)
    count = 0
    if damaged.any():
        _, labels = connected_components(step.graph[damaged][:, damaged], directed=False)
        sizes = np.bincount(labels, weights=step.areas[damaged])
        count = int((sizes > 0.01 * total).sum())
    return {"damaged_share": share, "damaged_components": count}


def breakthrough(world, plate, state: Mechanism):
    """Step 6 for one plate (the caller has applied the age and size gates). Returns the daughters or None."""
    step = state.steps.get(plate.plate_id)
    if step is None or not step.broken.any():
        return None
    p = state.params
    min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
    boundary = plate.exposed_sides().any(axis=1)
    count, labels = connected_components(step.graph[step.broken][:, step.broken], directed=False)
    cells = np.flatnonzero(step.broken)
    first, second = thin_shell._tangent_basis(step.points)
    # Largest band first: it is the likeliest to span the plate.
    order = np.argsort(-np.bincount(labels, weights=step.areas[cells], minlength=count))
    for component in order:
        band = np.zeros(len(step.points), dtype=bool)
        band[cells[labels == component]] = True
        if band.sum() < 2:
            continue
        chain = band_path(step, boundary, band)
        if chain is None or len(chain) < 2 or len(set(chain)) != len(chain):
            state.counters["path_overlap"] += int(chain is not None and len(set(chain)) != len(chain))
            continue
        if not boundary[chain[0]] or not boundary[chain[-1]] or chain[0] == chain[-1]:
            continue
        path = np.asarray(chain)
        normal, ds = _path_normals(step.points, path, first, second)
        if not proto.ligament_splits(step.tensor[path], normal, ds, ~step.broken[path], step.n_eff[path]):
            continue
        state.counters["ligament_pass"] += 1
        mask_a = split_masks(world, plate, step, path, min_nodes, state.counters.setdefault("rejected_because", {}))
        if mask_a is None:
            state.counters["split_rejected"] += 1
            continue
        pieces = plate._plates_from_node_masks([mask_a, ~mask_a], [plate.plate_id, world.next_plate_id])
        if len(pieces) != 2:
            state.counters["split_rejected"] += 1
            continue
        world.next_plate_id += 1
        a, b = pieces
        a.reset_age()
        b.reset_age()
        flow = mantle.flow_at(step.points, world.mantle_centers)
        pole_a = mantle.fit_euler_pole(step.points[mask_a], flow[mask_a])
        pole_b = mantle.fit_euler_pole(step.points[~mask_a], flow[~mask_a])
        a.set_omega(mantle.clamp_rate(pole_a))
        b.set_omega(mantle.clamp_rate(pole_b))
        a.set_internal_stress(0.0)
        b.set_internal_stress(0.0)
        near_path = np.zeros(len(step.points), dtype=bool)
        near_path[path] = True
        pattern = (step.graph > 0).astype(np.int8)
        near_path |= (pattern @ near_path.astype(np.int8)) > 0
        for daughter, mask in ((a, mask_a), (b, ~mask_a)):
            cratons.mark_belt(world, daughter, near_path[mask])
        oceanic_share = float(ds[step.oceanic[path]].sum() / max(ds.sum(), 1e-9))
        state.events.append({
            "kind": "split", "elapsed_years": world.elapsed_years, "plate_id": int(plate.plate_id), "daughters": [int(a.plate_id), int(b.plate_id)],
            "path_cells": int(len(path)), "path_length_m": float(ds.sum()), "oceanic_share": oceanic_share, "oceanic": oceanic_share > 0.5,
            "intact_share": float(ds[~step.broken[path]].sum() / max(ds.sum(), 1e-9)), "weak_share": float(ds[step.weak[path]].sum() / max(ds.sum(), 1e-9)),
            "plume_share": float(ds[step.plumed[path]].sum() / max(ds.sum(), 1e-9)),
            "random_weak_share": random_path_weak_share(step, boundary, int(plate.plate_id), world.elapsed_years),
            "small_share": float(min(step.areas[mask_a].sum(), step.areas[~mask_a].sum()) / step.areas.sum()),
            "path_keys": step.keys[path].tolist(), **localization(step, p),
        })
        return a, b
    return None


# --- hooks ---------------------------------------------------------------------------------------


def _damage_pass(world, years) -> None:
    state = _STATE
    started = time.perf_counter()
    state.steps = {}
    for plate in world.plates:
        damage_step(world, plate, years, state)
    state.seconds_damage += time.perf_counter() - started


def maybe_split_plate(world, plate, step_years: float = 0.0):
    """Replacement for ``merge_split.maybe_split_plate``: the damage pass for every plate (once per step), then the
    breakthrough decision for this plate under production's age and size gates."""
    state = _STATE
    if state.pass_key != world.elapsed_years:
        state.pass_key = world.elapsed_years
        _damage_pass(world, step_years)
    min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
    if plate.node_count() < min_nodes or plate.age_years < merge_split.SPLIT_MIN_AGE_YEARS:
        return None
    started = time.perf_counter()
    try:
        return breakthrough(world, plate, state)
    finally:
        state.seconds_split += time.perf_counter() - started


def open_rift(plate, world, rng_index, inserted_indices, rifted, stretch_share, years=0.0):
    """``quad_tectonics._open_rift`` plus damage inheritance: a margin cell inserted by rift opening takes the
    area-weighted damage of its donor band times its stretch share; its magmatic share is undamaged."""
    eps, warm = plate.collect("rift_strain"), plate.collect("rift_warm_strain")
    donor_eps = np.zeros(len(rifted))
    donor_warm = np.zeros(len(rifted))
    if len(rifted) and (eps.any() or warm.any()):
        points, _ = plate.all_points_and_elevation()
        areas = plate.node_areas_m2()
        graph, _, _ = mech.plate_graph(plate, points)
        pattern = (graph > 0).astype(np.int8)
        pre = np.ones(len(points), dtype=bool)
        pre[inserted_indices] = False
        reach = csr_matrix((np.ones(len(rifted), dtype=np.int8), (np.arange(len(rifted)), rifted)), shape=(len(rifted), len(points)))
        band = (reach @ pattern).tocsr()
        for _ in range(lithosphere_rows() - 1):
            band = band + (band.multiply(pre[None, :].astype(np.int8)) @ pattern)
        band = band.multiply(pre[None, :].astype(np.int8)).tocsr()
        band.data[:] = 1
        weights = band.multiply(areas[None, :]).tocsr()
        total = np.asarray(weights.sum(axis=1)).ravel()
        safe = np.where(total > 0.0, total, 1.0)
        donor_eps, donor_warm = (weights @ eps) / safe, (weights @ warm) / safe
    _ORIGINAL["open_rift"](plate, world, rng_index, inserted_indices, rifted, stretch_share, years)
    if len(rifted) and (donor_eps.any() or donor_warm.any()):
        eps, warm, broken = plate.collect("rift_strain"), plate.collect("rift_warm_strain"), plate.collect("rift_broken")
        eps[rifted] = donor_eps * stretch_share
        warm[rifted] = np.minimum(donor_warm * stretch_share, eps[rifted])
        broken[rifted] = False
        plate.set_fields_on_plate(rift_strain=eps, rift_warm_strain=warm, rift_broken=broken)


def lithosphere_rows() -> int:
    from app import lithosphere_plate

    return lithosphere_plate.K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION


def install(params: Params | None = None) -> Mechanism:
    """Register the fields and hook the step. Call before generating or loading any world."""
    global _STATE
    mech.register_fields()
    _STATE = Mechanism(params or Params())
    if not _ORIGINAL:
        _ORIGINAL["split"], _ORIGINAL["open_rift"] = merge_split.maybe_split_plate, quad_tectonics._open_rift
    merge_split.maybe_split_plate = maybe_split_plate
    quad_tectonics._open_rift = open_rift
    return _STATE


def uninstall() -> None:
    global _STATE
    if _ORIGINAL:
        merge_split.maybe_split_plate = _ORIGINAL["split"]
        quad_tectonics._open_rift = _ORIGINAL["open_rift"]
    _STATE = None

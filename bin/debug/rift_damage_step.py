"""Plate-level glue for the RFC 0018 rift-line prototype (D1, #888): the per-step damage pass for every plate, the
breakthrough decision and mask split, and the hooks that run them inside ``world.step_world``.

See ``rift_line.py`` for the line kernels and ``rift_damage_mechanism.py`` for the strength, softening and flank
transfer kernels and the list of readings of the RFC. Not production code.

``install(params)`` replaces ``merge_split.maybe_split_plate`` (the damage pass for every plate runs once per
step on its first call, then only the breakthrough decision runs, gated by ``SPLIT_MIN_AGE_YEARS`` and
``SPLIT_MIN_NODES`` as in production) and wraps ``quad_tectonics._open_rift`` so margin insertions inherit
their stretched share of damage. ``uninstall()`` restores both. ``prepare_plate`` and ``plate_line`` compute a plate's
rift line from a world without straining anything, for the D1 line check.
"""

import time
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import replay_gate2_capacity as gate2
import rift_damage_mechanism as mech
import rift_damage_prototype as proto
import rift_line as rl
from rift_damage_mechanism import Params, RADIUS_M

from app import cratons, elevation_lines, faults, lithosphere, mantle, merge_split, quad_tectonics, rheology, terranes, thin_shell  # noqa: E501
from app.elevation_lines import ELEV_CHANGE_RIFT, line_spacing_rad
from app.surface_fields import SURFACE_FIELDS


@dataclass
class LineInfo:
    """The rift line chosen for a plate this step and what is known about it."""

    chain: np.ndarray
    geometry: rl.Geometry
    direction: int
    pull: float  # N, at the step start
    strength: float  # N
    lam: float
    is_open: bool
    peak_share: float  # smaller-side area share at the corridor's peak section
    normal_load: np.ndarray  # n̂ᵀ 𝐍 n̂ per line cell, without the plume term


@dataclass
class PlatePrep:
    """A plate's frozen inputs for one step (RFC step 3: loads, N_cap, the line and its band are frozen)."""

    points: np.ndarray
    areas: np.ndarray
    graph: csr_matrix
    boundary: np.ndarray
    tensor: np.ndarray
    eps: np.ndarray
    warm: np.ndarray
    broken: np.ndarray
    hc: np.ndarray
    hm: np.ndarray
    continental: np.ndarray
    heated: np.ndarray
    plume_left_yr: np.ndarray  # per cell
    base_quiet: np.ndarray
    base_plume: np.ndarray
    kind: np.ndarray
    lag: np.ndarray
    k_scale: np.ndarray
    mult_quiet: np.ndarray
    mult_plume: np.ndarray
    n_eff: np.ndarray  # start-of-step strength per cell (plume weakening where heated)
    weak: np.ndarray
    disk: csr_matrix
    line: LineInfo | None
    candidates: int  # corridors that gave a connected line
    hm_strength: np.ndarray | None = None  # the mantle-lithosphere thickness behind N_cap (``hm`` unless ``oceanic_age_strength``)


@dataclass
class PlateStep:
    """What a plate's damage pass left behind: arrays for the recorder and the breakthrough decision."""

    plate_id: int
    keys: np.ndarray
    points: np.ndarray
    areas: np.ndarray
    graph: csr_matrix
    boundary: np.ndarray
    eps: np.ndarray
    d_eps: np.ndarray
    broken: np.ndarray
    oceanic: np.ndarray
    weak: np.ndarray  # inherited weakness marks (belt, volcano, fault) per cell
    plumed: np.ndarray  # plume heat present per cell
    elevation: np.ndarray
    hc: np.ndarray
    line: LineInfo | None = None
    pull_end: float = 0.0  # P and S at the end of the step, for the breakthrough test
    strength_end: float = 0.0
    broken_share: float = 0.0  # length share of the line that is broken at the end of the step
    seconds: float = 0.0
    stats: dict = field(default_factory=dict)


@dataclass
class Mechanism:
    params: Params
    steps: dict = field(default_factory=dict)  # plate_id -> PlateStep of the latest pass
    events: list = field(default_factory=list)  # one dict per split
    counters: dict = field(default_factory=lambda: {"held": 0, "far_flank": 0, "line_none": 0, "line_steps": 0, "growing_steps": 0,
                                                    "off_line_strain_cell_steps": 0, "substeps_max": 0, "break_tests": 0,
                                                    "split_rejected": 0, "rejected_because": {}})
    pass_key: float | None = None
    prev_reference: dict = field(default_factory=dict)
    seconds_damage: float = 0.0
    seconds_split: float = 0.0
    ledger_gap: dict = field(default_factory=dict)  # worst relative closure error seen per quantity
    on_line_strain: float = 0.0  # strain applied on the line's own cells (summed over steps)
    lam_samples: list = field(default_factory=list)  # (Λ, open, growth ratio) at the start of each growing step
    total_strain: float = 0.0  # strain applied anywhere


_STATE: Mechanism | None = None
_ORIGINAL: dict = {}


# --- the plate's frozen inputs and its line ------------------------------------------------------


def force_inputs(world, plate):
    """The per-node force field and the divergent-boundary mask (D0's event force field)."""
    return gate2._event_force_field(world, plate, with_divergent=True)


def fault_mask(world, plate, points, spacing_m) -> np.ndarray:
    """Cells within one cell spacing of a fault trace of this plate."""
    traces = [faults.fault_world_points(f, plate) for f in world.faults if f.plate_id == plate.plate_id]
    if not traces:
        return np.zeros(len(points), dtype=bool)
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


def _line_problem(prep: PlatePrep, line: LineInfo, p: Params):
    """The frozen ``LineProblem`` of a line, its plan for flank recipients, and the latch strain of every cell."""
    chain = line.chain
    band, reach_band, reach_line = rl.band_reach(prep.points, chain, prep.eps, p)
    in_band = np.zeros(len(prep.points), dtype=bool)
    in_band[band] = True
    plan = mech.flank_plan(prep.graph, in_band, prep.broken, p.flank_rows)
    floor_c, floor_m = rheology.RIFT_CRITICAL_THICKNESS_M, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
    # The strain at which the crust reaches its floor is the cell's strain now plus ln(Hc / floor) (Hc is already thinned).
    crust_room = np.log(np.maximum(prep.hc, floor_c) / floor_c)
    latch = np.minimum(p.eps_break, prep.eps + crust_room)
    if p.oceanic_latch_on_break_only:
        mantle_room = np.log(np.maximum(prep.hm, floor_m) / floor_m)
        latch = np.where(prep.continental, latch, np.minimum(p.eps_break, prep.eps + np.maximum(crust_room, mantle_room)))
    latch = np.where(prep.broken, prep.eps, np.maximum(latch, prep.eps))
    latch[plan.band[plan.held]] = prep.eps[plan.band[plan.held]]  # a held band cell does not strain this step
    prob = rl.LineProblem(
        chain=chain, s=line.geometry.s, ds=line.geometry.ds, normal_load=line.normal_load,
        base_quiet=prep.base_quiet[chain], base_plume=prep.base_plume[chain], plume_left_yr=prep.plume_left_yr[chain],
        intact=~prep.broken[chain], disk_rows=prep.disk[chain], band=band, reach_band=reach_band, reach_line=reach_line, latch=latch,
        capacity=rl.CapacityInputs(prep.kind[chain], prep.hc[chain], prep.hm_strength[chain], prep.lag[chain], prep.k_scale[chain],
                                   prep.mult_quiet[chain], prep.mult_plume[chain], prep.eps[chain]))
    return prob, plan, latch


def prepare_plate(world, plate, p: Params, years: float, eps=None, warm=None, broken=None) -> PlatePrep:
    """Steps 1-3 for one plate without straining anything: the load, strength and softening, the section pull, and the
    rift line (the best candidate over the corridors) with its section load."""
    points, _ = plate.all_points_and_elevation()
    n = len(points)
    spacing_m = line_spacing_rad(world.node_density) * RADIUS_M
    eps = plate.collect("rift_strain") if eps is None else eps
    warm = plate.collect("rift_warm_strain") if warm is None else warm
    broken = plate.collect("rift_broken").astype(bool) if broken is None else broken
    hc, hm = plate.collect("crustal_thickness_m"), plate.collect("mantle_lithosphere_thickness_m")
    codes, lag = plate.collect("crust_type_code"), plate.collect("moho_thermal_lag_c")
    heat = plate.collect("plume_heat_myr")
    continental = np.asarray(elevation_lines.effective_is_continental_from_codes(codes, plate.crust_type == "continental"), bool)
    field_, divergent = force_inputs(world, plate)
    areas = field_.areas_m2
    tensor = p.k_cut * proto.cut_tensor(points, areas, field_.total_n, spacing_m)
    sections = rl.section_profile(points, areas, field_.total_n, p.k_cut)
    heated = heat > 0.0
    plume_left = np.where(heated, np.minimum(heat * mech.YEAR_PER_MYR, years), 0.0)
    kind = np.where(continental, 0, np.where(divergent, 1, 2))
    k_scale = p.k_cap * np.where(continental, 1.0, p.oceanic_k_cap_factor)
    hm_strength = hm
    if p.oceanic_age_strength:
        created = plate.collect("node_created_years")
        age = np.where(created == SURFACE_FIELDS["node_created_years"].sentinel, np.nan, world.elapsed_years - created)
        hm_strength = np.where(continental, hm, mech.oceanic_strength_hm(hc, hm, age, eps, p))
    if p.craton_root_m > 0.0:
        root = p.craton_root_m * cratons.strength(plate.collect("craton_crust_m")) * np.exp(-eps)
        hm_strength = np.where(continental, hm + root, hm_strength)
    ncap = mech.capacity(kind, hc, hm_strength, lag) * k_scale
    belt = cratons.belt_mask(plate, world.elapsed_years)
    volcano = plate.collect("is_volcano").astype(bool)
    near_fault = fault_mask(world, plate, points, spacing_m)
    mult_plume = mech.inherited_multiplier(belt, volcano, heated, near_fault)
    mult_quiet = mech.inherited_multiplier(belt, volcano, np.zeros(n, bool), near_fault)
    if p.mark_power != 1.0:
        mult_plume = np.clip(mult_plume ** p.mark_power, mech.MULTIPLIER_FLOOR, 1.0)
        mult_quiet = np.clip(mult_quiet ** p.mark_power, mech.MULTIPLIER_FLOOR, 1.0)
    base_plume, base_quiet = ncap * mult_plume, ncap * mult_quiet
    disk = mech.disk_average_matrix(points, areas, p.w0_m / 2.0)
    w = mech.softening(disk @ eps, disk @ warm, p)
    n_eff = np.where(heated, base_plume, base_quiet) * w
    graph, _, _ = mech.plate_graph(plate, points)
    boundary = plate.exposed_sides().any(axis=1)
    prep = PlatePrep(points, areas, graph, boundary, tensor, eps, warm, broken, hc, hm, continental, heated, plume_left, base_quiet,
                     base_plume, kind, lag, k_scale, mult_quiet, mult_plume, n_eff, belt | volcano | near_fault, disk, None, 0)
    prep.hm_strength = hm_strength
    prep.line, prep.candidates = _choose_line(prep, sections, p)
    return prep


def _choose_line(prep: PlatePrep, sections: rl.Sections, p: Params):
    first, second = thin_shell._tangent_basis(prep.points)
    cost = np.where(prep.broken, 0.0, prep.n_eff * prep.mult_quiet ** (p.route_mark_power - 1.0))
    plume_tensor = prep.tensor + p.n_plume * prep.heated[:, None] * np.array([1.0, 1.0, 0.0])
    best, best_key, count = None, None, 0
    for corridor in rl.candidate_corridors(sections, p):
        chain = rl.find_line(prep.graph, corridor.cells, prep.boundary, sections.strike[:, corridor.direction], cost)
        if chain is None or len(chain) < 2:
            continue
        count += 1
        geometry = rl.path_geometry(prep.points, chain, first, second)
        normal_load = rl.path_normal_load(prep.tensor[chain], geometry.normal)
        load = np.maximum(rl.path_normal_load(plume_tensor[chain], geometry.normal), 0.0)
        pull = float(np.sum(load * geometry.ds))
        strength = float(np.sum((prep.n_eff[chain] * geometry.ds)[~prep.broken[chain]]))
        lam, is_open = rl.section_load(pull, strength, p)
        key = rl.rank_key(pull, lam, is_open)
        if best_key is None or key > best_key:
            best_key = key
            best = LineInfo(chain, geometry, corridor.direction, pull, strength, lam, is_open, corridor.peak_share, normal_load)
    return best, count


def plate_line(world, plate, p: Params | None = None) -> dict:
    """The line check's record for one plate: its rift line without straining anything, and whether cutting along it
    gives a valid two-daughter split (RFC D1 line check)."""
    p = p or Params()
    prep = prepare_plate(world, plate, p, 1e6)
    line = prep.line
    if line is None:
        return {"plate": int(plate.plate_id), "cells": len(prep.points), "line": False, "corridors": prep.candidates}
    chain, geometry = line.chain, line.geometry
    mask, why = rl.split_masks(prep.graph, prep.areas, chain, p, points=prep.points)
    craton = plate.collect("craton_crust_m") > 0.0
    oceanic = ~prep.continental
    share = lambda flag: float(geometry.ds[flag[chain]].sum() / max(geometry.ds.sum(), 1e-9))  # noqa: E731
    return {
        "plate": int(plate.plate_id), "cells": len(prep.points), "line": True, "corridors": prep.candidates, "valid": mask is not None,
        "reject": why, "smaller_share": None if mask is None else float(min(prep.areas[mask].sum(), prep.areas[~mask].sum()) / prep.areas.sum()),
        "peak_share": line.peak_share, "length_ratio": geometry.length / max(geometry.straight, 1.0), "length_km": geometry.length / 1e3,
        "craton_share": share(craton), "chord_craton_share": chord_craton_share(prep.points, craton, chain), "oceanic_share": share(oceanic),
        "lam": line.lam, "open": line.is_open, "pull_tn_per_m_mean": float(line.pull / max(geometry.ds.sum(), 1.0) / 1e12),
        "weak_share": share(prep.weak), "cells_on_line": len(chain),
        "cratons_in_corridor": bool(craton.any()),
    }


def chord_craton_share(points, craton, chain, samples: int = 200) -> float:
    """Craton share along the great-circle chord between the line's end cells (nearest-cell lookup)."""
    a, b = points[chain[0]], points[chain[-1]]
    angle = float(np.arccos(np.clip(a @ b, -1.0, 1.0)))
    if angle < 1e-9:
        return 0.0
    t = np.linspace(0.0, 1.0, samples)
    path = (np.sin((1 - t) * angle)[:, None] * a + np.sin(t * angle)[:, None] * b) / np.sin(angle)
    _, nearest = cKDTree(points).query(path)
    return float(craton[nearest].mean())


# --- the damage pass -----------------------------------------------------------------------------


def damage_step(world, plate, years: float, state: Mechanism) -> PlateStep | None:
    """Steps 1-4 for one plate: load, strength, the rift line and its strain, thinning with flank transfer, the latch."""
    p = state.params
    n = plate.node_count()
    if n == 0 or years <= 0.0:
        return None
    started = time.perf_counter()
    eps0, warm0 = plate.collect("rift_strain"), plate.collect("rift_warm_strain")
    broken = plate.collect("rift_broken").astype(bool)
    reset_after_eruption(state, plate, eps0, warm0, broken)
    prep = prepare_plate(world, plate, p, years, eps0, warm0, broken)
    line = prep.line
    state.counters["line_steps"] += line is not None
    state.counters["line_none"] += line is None
    codes = plate.collect("crust_type_code")
    eps1, warm1, d_eps, newly_broken = eps0, warm0 * np.exp(-years / p.tau_heal_yr), np.zeros(n), np.zeros(n, bool)
    hc, hm = prep.hc, prep.hm
    out, plan, advance, prob = None, None, None, None
    if line is not None:
        prob, plan, latch = _line_problem(prep, line, p)
        advance = rl.advance_rift(prob, eps0, warm0, years, p)
        eps1, warm1 = advance.eps, np.minimum(advance.warm, advance.eps)
        d_eps = np.maximum(eps1 - eps0, 0.0)
        state.counters["growing_steps"] += float(d_eps.max()) > 0.0
        state.counters["substeps_max"] = max(state.counters["substeps_max"], advance.substeps)
        state.counters["held"] += int(plan.held.sum())
        state.counters["far_flank"] += int((plan.hops > 4).sum())
        state.total_strain += float(d_eps.sum())
        state.on_line_strain += float(d_eps[line.chain].sum())
        state.lam_samples.append((float(advance.lam_start), bool(advance.open_start), float(advance.ratio_start)))
    floor_c, floor_m = rheology.RIFT_CRITICAL_THICKNESS_M, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
    columns = {
        "hc": hc, "hm": hm, "material": plate.collect("continental_material_m"), "craton": plate.collect("craton_crust_m"),
        "cover": plate.collect("mobile_cover_m"), "cover_material": plate.collect("mobile_cover_continental_m"),
    }
    plateau, arc, moment = terranes.plate_tracers(plate)
    columns.update(plateau=plateau, arc=arc, moment=moment)
    before_totals = _volumes(columns, prep.areas)
    if plan is not None and d_eps.max() > 0.0:
        out = mech.thin_with_flank_transfer(prep.areas, plan, d_eps, columns)
        at_floor = out.hc <= floor_c + 1e-6
        if p.oceanic_latch_on_break_only:
            at_floor = np.where(prep.continental, at_floor, at_floor & (out.hm <= floor_m + 1e-6))
        newly_broken = (d_eps > 0.0) & ~broken & ((eps1 >= p.eps_break - 1e-9) | at_floor)
        after_totals = _volumes({**columns, **{k: getattr(out, k) for k in columns}}, prep.areas)
        _record_closure(state, before_totals, after_totals, out.craton_lost_m3)
    density = lithosphere.node_crust_density(codes, plate.crust_type)
    elevation = plate.collect("elevation")
    if out is not None:
        hc, hm = out.hc, out.hm
        changed = (out.hc != prep.hc) | (out.hm != prep.hm)
        shift = lithosphere.isostatic_elevation(out.hc, out.hm, density) - lithosphere.isostatic_elevation(prep.hc, prep.hm, density)
        new_elevation = np.where(changed, rheology.clip_elevation_bounds(elevation + shift), elevation)
        reason = plate.collect("elev_change_reason")
        reason[(d_eps > 0.0) & changed] = ELEV_CHANGE_RIFT
        plate.set_fields_on_plate(
            crustal_thickness_m=out.hc, mantle_lithosphere_thickness_m=out.hm, continental_material_m=out.material,
            craton_crust_m=out.craton, mobile_cover_m=out.cover, mobile_cover_continental_m=np.minimum(out.cover_material, out.material),
            elevation=new_elevation, elev_change_reason=reason)
        if "arc_crust_m" in mech.SURFACE_FIELDS:
            plate.set_fields_on_plate(arc_crust_m=out.arc)
        plate.set_fields_on_plate(plateau_crust_m=out.plateau, buoyant_crust_time_m_yr=out.moment)
        if out.craton_lost_m3 > 0.0:
            cratons.record(world, "rifted_stretch_m3", out.craton_lost_m3)
    else:
        new_elevation = elevation
    broken_now = broken | newly_broken
    plate.set_fields_on_plate(rift_strain=eps1, rift_warm_strain=warm1, rift_broken=broken_now)
    pull_end = strength_end = broken_share = 0.0
    if line is not None:
        fresh = rl.LineProblem(**{**prob.__dict__, "intact": ~broken_now[line.chain]})
        n_eff, load, _ = rl.line_state(fresh, eps1, warm1, years, p)
        pull_end = float(np.sum(load * fresh.ds))
        strength_end = float(np.sum((n_eff * fresh.ds)[fresh.intact]))
        broken_share = float(line.geometry.ds[broken_now[line.chain]].sum() / max(line.geometry.ds.sum(), 1e-9))
        outside = np.ones(n, dtype=bool)
        outside[prob.band] = False
        state.counters["off_line_strain_cell_steps"] += int(((d_eps > 0.0) & outside).sum())
    state.prev_reference[plate.plate_id] = (plate.cell_keys.copy(), _is_reference(hc, hm))
    step = PlateStep(plate.plate_id, plate.cell_keys.copy(), prep.points, prep.areas, prep.graph, prep.boundary, eps1, d_eps, broken_now,
                     ~prep.continental, prep.weak, prep.heated, new_elevation, hc, line, pull_end, strength_end, broken_share,
                     time.perf_counter() - started)
    state.steps[plate.plate_id] = step
    return step


def _volumes(columns, areas) -> dict:
    return {name: float(np.dot(columns[name], areas)) for name in
            ("hc", "hm", "material", "cover", "cover_material", "plateau", "arc", "moment", "craton")}


def _record_closure(state, before, after, craton_lost) -> None:
    """Interior thinning moves volume; it must conserve every total. Craton is booked to its sink."""
    for name, value in before.items():
        scale = max(abs(value), 1.0)
        expected = value - (craton_lost if name == "craton" else 0.0)
        gap = abs(after[name] - expected) / scale
        state.ledger_gap[name] = max(state.ledger_gap.get(name, 0.0), gap)


# --- breakthrough --------------------------------------------------------------------------------


def _tie_side_factory(world, step: PlateStep):
    flow = mantle.flow_at(step.points, world.mantle_centers)

    def tie(side, cell) -> int:
        best, score = 0, -np.inf
        for k in (0, 1):
            sel = side == k
            if not sel.any():
                continue
            pole = mantle.fit_euler_pole(step.points[sel], flow[sel])
            value = float(np.cross(pole, step.points[cell]) @ flow[cell])
            if value > score:
                best, score = k, value
        return best

    return tie


def random_path_weak_share(step: PlateStep, plate_id: int, elapsed_years: float, trials: int = 20) -> float:
    """Mean inherited-weakness share along matched random edge-to-edge paths: shortest paths between random pairs of
    boundary cells at least 60% of the plate's longest boundary-to-boundary span apart (the falsification test's
    baseline). Deterministic in (plate, time)."""
    from scipy.sparse.csgraph import dijkstra

    rng = np.random.default_rng((plate_id, round(elapsed_years)))
    cells = np.flatnonzero(step.boundary)
    if len(cells) < 4:
        return float("nan")
    span = 0.6 * float(np.max(np.linalg.norm(step.points[cells][:, None] - step.points[cells][None, :], axis=2)))
    shares = []
    for _ in range(trials * 5):
        a, b = rng.choice(cells, 2, replace=False)
        if np.linalg.norm(step.points[a] - step.points[b]) < span:
            continue
        _, pred = dijkstra(step.graph, directed=False, indices=a, return_predecessors=True)
        shares.append(float(step.weak[rl._trace(pred, b)].mean()))
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


def near_path_mask(step: PlateStep, path: np.ndarray) -> np.ndarray:
    """The path cells and their edge neighbours: the one trace-derived band for belt marking and craton accounting."""
    mask = np.zeros(len(step.points), dtype=bool)
    mask[path] = True
    return mask | ((step.graph > 0).astype(np.int8) @ mask.astype(np.int8) > 0)


def record_cut_mask(world, plate, band: np.ndarray) -> None:
    """``cratons.record_rift_cut`` for a mask band: the split's count, band area and craton area in the craton ledger."""
    cratons.ensure_ledger(world)
    areas = cratons._areas(world, plate)[band]
    craton_area = float(areas[plate.collect("craton_crust_m")[band] > 0.0].sum())
    ledger = world.craton_ledger
    ledger["split_cuts"] += 1.0
    ledger["split_cuts_through_craton"] += float(craton_area > 0.0)
    ledger["split_band_area_m2"] += float(areas.sum())
    ledger["split_band_craton_area_m2"] += craton_area


def breakthrough(world, plate, state: Mechanism):
    """Step 6 for one plate (the caller has applied the age and size gates). Returns the daughters or None."""
    step = state.steps.get(plate.plate_id)
    p = state.params
    if step is None or step.line is None or not rl.ligament_gives(step.pull_end, step.strength_end) or step.broken_share < p.f_break:
        return None
    state.counters["break_tests"] += 1
    path = step.line.chain
    mask_a, why = rl.split_masks(step.graph, step.areas, path, p, _tie_side_factory(world, step), step.points)
    if mask_a is None:
        state.counters["split_rejected"] += 1
        reasons = state.counters["rejected_because"]
        reasons[why] = reasons.get(why, 0) + 1
        return None
    geometry = step.line.geometry
    record_cut_mask(world, plate, near_path_mask(step, path))
    pieces = plate._plates_from_node_masks([mask_a, ~mask_a], [plate.plate_id, world.next_plate_id])
    if len(pieces) != 2:
        state.counters["split_rejected"] += 1
        return None
    world.next_plate_id += 1
    a, b = pieces
    reference = state.prev_reference.get(plate.plate_id)
    if reference is not None:  # daughters keep their cells' keys, so the eruption-reset memory carries over
        for daughter in (a, b):
            known = np.isin(reference[0], daughter.cell_keys)
            state.prev_reference[daughter.plate_id] = (reference[0][known], reference[1][known])
    a.reset_age()
    b.reset_age()
    flow = mantle.flow_at(step.points, world.mantle_centers)
    a.set_omega(mantle.clamp_rate(mantle.fit_euler_pole(step.points[mask_a], flow[mask_a])))
    b.set_omega(mantle.clamp_rate(mantle.fit_euler_pole(step.points[~mask_a], flow[~mask_a])))
    a.set_internal_stress(0.0)
    b.set_internal_stress(0.0)
    near_path = near_path_mask(step, path)
    for daughter, mask in ((a, mask_a), (b, ~mask_a)):
        cratons.mark_belt(world, daughter, near_path[mask])
    ds = geometry.ds
    share = lambda flag: float(ds[flag[path]].sum() / max(ds.sum(), 1e-9))  # noqa: E731
    oceanic_share = share(step.oceanic)
    state.events.append({
        "kind": "split", "elapsed_years": world.elapsed_years, "plate_id": int(plate.plate_id), "daughters": [int(a.plate_id), int(b.plate_id)],
        "path_cells": int(len(path)), "path_length_m": float(ds.sum()), "length_ratio": geometry.length / max(geometry.straight, 1.0),
        "oceanic_share": oceanic_share, "oceanic": oceanic_share > 0.5, "intact_share": float(ds[~step.broken[path]].sum() / max(ds.sum(), 1e-9)),
        "weak_share": share(step.weak), "plume_share": share(step.plumed), "open": step.line.is_open,
        "random_weak_share": random_path_weak_share(step, int(plate.plate_id), world.elapsed_years),
        "small_share": float(min(step.areas[mask_a].sum(), step.areas[~mask_a].sum()) / step.areas.sum()),
        "path_keys": step.keys[path].tolist(), **localization(step, p),
    })
    return a, b


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

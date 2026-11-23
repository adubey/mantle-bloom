"""Plate breakup along rift lines (RFC 0018, #846): the per-step damage pass and the split that follows a broken line.

This is the alternative to the great-circle split in `merge_split.maybe_split_plate`. It runs only when
`World.split_method == "rift_lines"`; the default, `"great_circle"`, never reaches this module, so existing worlds
step exactly as before. Knobs are `rift_strength.RiftLineParams`, overridden per world by `World.rift_line_params`.

Each step, `merge_split.apply_topology_changes` calls:

1. `damage_pass(world, years)`: for every plate, the load on it (`rift_strength.cut_tensor`), its strength, one rift
   line where the pull across the plate is largest and the rock cheapest to cross (`rift_lines`), the strain on that
   line's band, and the thinning that strain causes, with the removed volume moved to the band's flanks. Strain is
   stored in the rock, in the surface fields `rift_strain`, `rift_warm_strain` and `rift_broken`, so a rift that
   stalls leaves a scar the next line can reuse.
2. `breakthrough(world, plate)` for one plate at a time, under `maybe_split_plate`'s age and size gates: when at least
   `f_break` of the line is broken and the pull across it exceeds the strength left, the plate splits along it into
   two daughters.

`quad_tectonics._open_rift` calls `margin_donor_damage` and `inherit_margin_damage` so crust stretched into a new
margin cell carries its donors' strain.

`World.rift_runtime` holds the last pass's per-plate records, the split events and counters for diagnostics. It is
rebuilt on demand and never saved.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import cratons, elevation_lines, faults, lithosphere, mantle, rheology, terranes, thin_shell, torque, trace
from . import rift_lines as lines
from . import rift_strength as strength
from .elevation_lines import ELEV_CHANGE_RIFT, line_spacing_rad
from .rift_strength import RADIUS_M, RiftLineParams
from .surface_fields import SURFACE_FIELDS

SPLIT_METHODS = ("great_circle", "rift_lines")
# A flank transfer that had to look this many hops for a recipient is counted (RFC 0018 step 4 reports them).
FAR_FLANK_HOPS = 4
# Damaged patches larger than this share of a plate count as separate components at a split (#699's failure mode).
LOCALIZATION_COMPONENT_SHARE = 0.01


@dataclass
class LineInfo:
    """The rift line chosen for a plate this step and what is known about it."""

    chain: np.ndarray
    geometry: lines.Geometry
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
    hm_strength: np.ndarray  # the mantle-lithosphere thickness behind N_cap (`hm` unless a strength-map knob is on)
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
    line: LineInfo | None = None
    candidates: int = 0  # corridors that gave a connected line


@dataclass
class PlateStep:
    """What a plate's damage pass left behind, for the breakthrough decision and for diagnostics."""

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


def _new_counters() -> dict:
    return {"held": 0, "far_flank": 0, "line_none": 0, "line_steps": 0, "growing_steps": 0, "substeps_max": 0, "break_tests": 0,
            "split_rejected": 0, "rejected_because": {}, "cap_held_m3": 0.0}


@dataclass
class RiftRuntime:
    """Per-world scratch of the rift-line method. Not saved: a loaded world starts with an empty one."""

    params: RiftLineParams
    overrides: dict  # the `World.rift_line_params` these params were built from
    steps: dict = field(default_factory=dict)  # plate_id -> PlateStep of the latest pass
    events: list = field(default_factory=list)  # one dict per split
    counters: dict = field(default_factory=_new_counters)
    ledger_gap: dict = field(default_factory=dict)  # worst relative closure error of interior thinning, per quantity


def runtime(world) -> RiftRuntime:
    """The world's rift scratch, with params rebuilt whenever `World.rift_line_params` has changed."""
    state = world.rift_runtime
    if not isinstance(state, RiftRuntime):
        state = world.rift_runtime = RiftRuntime(RiftLineParams.with_overrides(world.rift_line_params), dict(world.rift_line_params))
    elif state.overrides != world.rift_line_params:
        state.params, state.overrides = RiftLineParams.with_overrides(world.rift_line_params), dict(world.rift_line_params)
    return state


# --- the plate's frozen inputs and its line ------------------------------------------------------


def force_inputs(world, plate):
    """The torque engine's per-node force field on a plate, and its divergent-boundary mask."""
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    cache = world.boundary_search_cache
    neighbours = plate.get_neighbours(world.plates)
    inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
    torque.apply_force_forms(world, plate, neighbours, inputs, torque.force_forms(world))
    _convergent, divergent, _transform, contested = torque.classify_boundary_nodes(plate, neighbours, inputs, reach, cache)
    collision = torque.continental_contact_mask(plate, inputs, contested)
    subducting = torque.subducting_boundary_mask(plate, inputs, reach)
    craton_contact = torque.contact_craton_m(plate, neighbours, inputs)
    jam = torque.contact_jam_factor(plate, world, neighbours, inputs, reach)
    forces = torque.gather_per_node_force_field(plate, world, inputs, spacing, subducting, divergent, collision, craton_contact, jam)
    return forces, np.asarray(divergent, dtype=bool)


def fault_mask(world, plate, points, spacing_m) -> np.ndarray:
    """Cells within one cell spacing of a fault trace of this plate."""
    traces = [faults.fault_world_points(f, plate) for f in world.faults if f.plate_id == plate.plate_id]
    if not traces:
        return np.zeros(len(points), dtype=bool)
    distance, _ = cKDTree(np.concatenate(traces)).query(points)
    return distance * RADIUS_M <= spacing_m


def plate_graph(plate, points) -> csr_matrix:
    """Symmetric CSR adjacency of the plate's cells with centre distances (m) as data."""
    adjacency = plate.adjacency()
    n = len(points)
    i = np.repeat(np.arange(n), np.diff(adjacency.offsets))
    j = np.asarray(adjacency.neighbours)
    length = RADIUS_M * np.arccos(np.clip(np.sum(points[i] * points[j], axis=1), -1.0, 1.0))
    return csr_matrix((np.maximum(length, 1.0), (i, j)), shape=(n, n))


def _is_reference(hc, hm) -> np.ndarray:
    return ((hc == lithosphere.REFERENCE_HC_CONTINENTAL_M) & (hm == lithosphere.REFERENCE_HM_CONTINENTAL_M)) | (
        (hc == lithosphere.REFERENCE_HC_OCEANIC_M) & (hm == lithosphere.YOUNG_RIDGE_HM_M))


def reset_erupted(hc, hm, eps, warm, broken) -> None:
    """Zero the damage of cells that erupted since the last pass (in place): new rock carries no inherited damage.
    `_erupt_melted_nodes` leaves exactly a reference column, and a strained cell cannot otherwise hold one, because
    the strain that damaged it thinned it. So a damaged cell with a reference column has erupted."""
    reset = _is_reference(hc, hm) & ((eps > 0.0) | (warm > 0.0) | broken)
    eps[reset], warm[reset], broken[reset] = 0.0, 0.0, False


def _strength_thickness(world, plate, p: RiftLineParams, hc, hm, eps, continental) -> np.ndarray:
    """The mantle-lithosphere thickness behind N_cap: the model's, or the strength-map knobs' replacement (#888)."""
    out = hm
    if p.oceanic_age_strength:
        created = plate.collect("node_created_years")
        age = np.where(created == SURFACE_FIELDS["node_created_years"].sentinel, np.nan, world.elapsed_years - created)
        out = np.where(continental, hm, strength.oceanic_strength_hm(hc, hm, age, eps, p))
    if p.craton_root_m > 0.0:
        root = p.craton_root_m * cratons.strength(plate.collect("craton_crust_m")) * np.exp(-eps)
        out = np.where(continental, hm + root, out)
    return out


def prepare_plate(world, plate, p: RiftLineParams, years: float, eps=None, warm=None, broken=None) -> PlatePrep:
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
    forces, divergent = force_inputs(world, plate)
    areas = forces.areas_m2
    balanced = strength.balanced_tangent_force(points, areas, forces.total_n)
    tensor = p.k_cut * strength.cut_tensor(points, areas, balanced, spacing_m)
    sections = lines.section_profile(points, areas, balanced, p.k_cut)
    heated = heat > 0.0
    # The plume's true remaining time, not capped at the step: the breakthrough test reads the line at t = years, and
    # capping would switch every plume off there, though the sub-steps saw it (PR #950 review). `advance_rift` only
    # sub-steps at expiries inside the step.
    plume_left = np.where(heated, heat * strength.YEAR_PER_MYR, 0.0)
    kind = np.where(continental, 0, np.where(divergent, 1, 2))
    k_scale = p.k_cap * np.where(continental, 1.0, p.oceanic_k_cap_factor)
    hm_strength = _strength_thickness(world, plate, p, hc, hm, eps, continental)
    ncap = strength.capacity(kind, hc, hm_strength, lag) * k_scale
    belt = cratons.belt_mask(plate, world.elapsed_years)
    volcano = plate.collect("is_volcano").astype(bool)
    near_fault = fault_mask(world, plate, points, spacing_m)
    mult_plume = strength.inherited_multiplier(belt, volcano, heated, near_fault, p.mark_power)
    mult_quiet = strength.inherited_multiplier(belt, volcano, np.zeros(n, bool), near_fault, p.mark_power)
    base_plume, base_quiet = ncap * mult_plume, ncap * mult_quiet
    disk = strength.disk_average_matrix(points, areas, p.w0_m / 2.0)
    w = strength.softening(disk @ eps, disk @ warm, p)
    n_eff = np.where(heated, base_plume, base_quiet) * w
    prep = PlatePrep(points, areas, plate_graph(plate, points), plate.exposed_sides().any(axis=1), tensor, eps, warm, broken, hc, hm,
                     hm_strength, continental, heated, plume_left, base_quiet, base_plume, kind, lag, k_scale, mult_quiet, mult_plume,
                     n_eff, belt | volcano | near_fault, disk)
    prep.line, prep.candidates = _choose_line(prep, sections, p)
    return prep


def _choose_line(prep: PlatePrep, sections: lines.Sections, p: RiftLineParams):
    """The corridor line with the largest section load `Λ` (open lines first). The path cost may weigh weakness marks
    more than the strength in `Λ` does (`route_mark_power`), so the route follows them while the open-or-not test
    keeps the rock's strength."""
    first, second = thin_shell._tangent_basis(prep.points)
    cost = np.where(prep.broken, 0.0, prep.n_eff * prep.mult_quiet ** (p.route_mark_power - 1.0))
    plume_tensor = prep.tensor + p.n_plume * prep.heated[:, None] * np.array([1.0, 1.0, 0.0])
    best, best_key, count = None, None, 0
    for corridor in lines.candidate_corridors(sections, p):
        chain = lines.find_line(prep.graph, corridor.cells, prep.boundary, sections.strike[:, corridor.direction], cost)
        if chain is None or len(chain) < 2:
            continue
        count += 1
        geometry = lines.path_geometry(prep.points, chain, first, second)
        normal_load = lines.path_normal_load(prep.tensor[chain], geometry.normal)
        load = np.maximum(lines.path_normal_load(plume_tensor[chain], geometry.normal), 0.0)
        pull = float(np.sum(load * geometry.ds))
        intact_strength = float(np.sum((prep.n_eff[chain] * geometry.ds)[~prep.broken[chain]]))
        lam, is_open = lines.section_load(pull, intact_strength, p)
        key = lines.rank_key(pull, lam, is_open)
        if best_key is None or key > best_key:
            best_key = key
            best = LineInfo(chain, geometry, corridor.direction, pull, intact_strength, lam, is_open, corridor.peak_share, normal_load)
    return best, count


def _line_problem(prep: PlatePrep, line: LineInfo, p: RiftLineParams):
    """The frozen `LineProblem` of a line and its plan for flank recipients."""
    chain = line.chain
    band, reach_band, reach_line = lines.band_reach(prep.points, chain, prep.eps, p)
    in_band = np.zeros(len(prep.points), dtype=bool)
    in_band[band] = True
    plan = strength.flank_plan(prep.graph, in_band, prep.broken, p.flank_rows, group=prep.continental)
    floor_c, floor_m = rheology.RIFT_CRITICAL_THICKNESS_M, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
    # The strain at which the crust reaches its floor is the cell's strain now plus ln(Hc / floor) (Hc is already thinned).
    crust_room = np.log(np.maximum(prep.hc, floor_c) / floor_c)
    latch = np.minimum(p.eps_break, prep.eps + crust_room)
    if p.oceanic_latch_on_break_only:
        mantle_room = np.log(np.maximum(prep.hm, floor_m) / floor_m)
        latch = np.where(prep.continental, latch, np.minimum(p.eps_break, prep.eps + np.maximum(crust_room, mantle_room)))
    latch = np.where(prep.broken, prep.eps, np.maximum(latch, prep.eps))
    latch[plan.band[plan.held]] = prep.eps[plan.band[plan.held]]  # a held band cell does not strain this step
    prob = lines.LineProblem(
        chain=chain, s=line.geometry.s, ds=line.geometry.ds, normal_load=line.normal_load,
        base_quiet=prep.base_quiet[chain], base_plume=prep.base_plume[chain], plume_left_yr=prep.plume_left_yr[chain],
        intact=~prep.broken[chain], disk_rows=prep.disk[chain], band=band, reach_band=reach_band, reach_line=reach_line, latch=latch,
        capacity=lines.CapacityInputs(prep.kind[chain], prep.hc[chain], prep.hm_strength[chain], prep.lag[chain], prep.k_scale[chain],
                                      prep.mult_quiet[chain], prep.mult_plume[chain], prep.eps[chain]))
    return prob, plan


# --- the damage pass -----------------------------------------------------------------------------


def _columns(plate, hc, hm) -> dict:
    plateau, arc, moment = terranes.plate_tracers(plate)
    return {"hc": hc, "hm": hm, "material": plate.collect("continental_material_m"), "craton": plate.collect("craton_crust_m"),
            "cover": plate.collect("mobile_cover_m"), "cover_material": plate.collect("mobile_cover_continental_m"),
            "plateau": plateau, "arc": arc, "moment": moment}


def _volumes(columns: dict, areas) -> dict:
    return {name: float(np.dot(columns[name], areas)) for name in strength.COLUMNS}


def _record_closure(state: RiftRuntime, before: dict, after: dict, craton_lost: float) -> None:
    """Interior thinning moves volume, so it must conserve every total. Craton is booked to its sink."""
    for name, value in before.items():
        expected = value - (craton_lost if name == "craton" else 0.0)
        gap = abs(after[name] - expected) / max(abs(value), 1.0)
        state.ledger_gap[name] = max(state.ledger_gap.get(name, 0.0), gap)


def _write_thinning(world, plate, prep: PlatePrep, out: strength.ThinningOutcome, d_eps) -> np.ndarray:
    """Write the thinned columns and the isostatic subsidence they cause. Returns the new elevation."""
    density = lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type)
    elevation = plate.collect("elevation")
    changed = (out.hc != prep.hc) | (out.hm != prep.hm)
    shift = lithosphere.isostatic_elevation(out.hc, out.hm, density) - lithosphere.isostatic_elevation(prep.hc, prep.hm, density)
    new_elevation = np.where(changed, rheology.clip_elevation_bounds(elevation + shift), elevation)
    reason = plate.collect("elev_change_reason")
    reason[(d_eps > 0.0) & changed] = ELEV_CHANGE_RIFT
    plate.set_fields_on_plate(
        crustal_thickness_m=out.hc, mantle_lithosphere_thickness_m=out.hm, continental_material_m=out.material,
        craton_crust_m=out.craton, mobile_cover_m=out.cover, mobile_cover_continental_m=np.minimum(out.cover_material, out.material),
        elevation=new_elevation, elev_change_reason=reason, plateau_crust_m=out.plateau, buoyant_crust_time_m_yr=out.moment)
    if "arc_crust_m" in SURFACE_FIELDS:
        plate.set_fields_on_plate(arc_crust_m=out.arc)
    if out.craton_lost_m3 > 0.0:
        cratons.record(world, "rifted_stretch_m3", out.craton_lost_m3)
    return new_elevation


def _strain_line(prep: PlatePrep, state: RiftRuntime, years: float):
    """Strain the line's band for the step. Returns (problem, flank plan, strain after, warm strain after)."""
    p, line = state.params, prep.line
    prob, plan = _line_problem(prep, line, p)
    advance = lines.advance_rift(prob, prep.eps, prep.warm, years, p)
    counters = state.counters
    counters["growing_steps"] += bool(np.any(advance.eps > prep.eps))
    counters["substeps_max"] = max(counters["substeps_max"], advance.substeps)
    counters["held"] += int(plan.held.sum())
    counters["far_flank"] += int((plan.hops > FAR_FLANK_HOPS).sum())
    return prob, plan, advance.eps, np.minimum(advance.warm, advance.eps)


@dataclass
class _Pending:
    """One plate's damage step, computed but not yet written to the plate."""

    prep: PlatePrep
    prob: "lines.LineProblem | None"
    eps1: np.ndarray
    warm1: np.ndarray
    d_eps: np.ndarray
    broken_now: np.ndarray
    out: "strength.ThinningOutcome | None"


def _compute_step(world, plate, years: float, state: RiftRuntime) -> _Pending | None:
    """Steps 1-4 for one plate (load, strength, the rift line and its strain, thinning with flank transfer, the
    latch), without writing anything to the plate."""
    p = state.params
    n = plate.node_count()
    if n == 0 or years <= 0.0:
        return None
    eps0, warm0 = plate.collect("rift_strain"), plate.collect("rift_warm_strain")
    broken = plate.collect("rift_broken").astype(bool)
    reset_erupted(plate.collect("crustal_thickness_m"), plate.collect("mantle_lithosphere_thickness_m"), eps0, warm0, broken)
    prep = prepare_plate(world, plate, p, years, eps0, warm0, broken)
    line = prep.line
    state.counters["line_steps" if line is not None else "line_none"] += 1
    healed = warm0 * np.exp(-years / p.tau_heal_yr)
    pending = _Pending(prep, None, eps0, healed, np.zeros(n), broken, None)
    if line is None:
        return pending
    prob, plan, eps1, warm1 = _strain_line(prep, state, years)
    d_eps = np.maximum(eps1 - eps0, 0.0)
    pending.prob, pending.eps1, pending.warm1, pending.d_eps = prob, eps1, warm1, d_eps
    if d_eps.max() <= 0.0:
        return pending
    columns = _columns(plate, prep.hc, prep.hm)
    out = strength.thin_with_flank_transfer(prep.areas, plan, d_eps, columns)
    _record_closure(state, _volumes(columns, prep.areas), _volumes({k: getattr(out, k) for k in strength.COLUMNS}, prep.areas),
                    out.craton_lost_m3)
    state.counters["cap_held_m3"] += out.cap_held_m3
    # A donor whose recipients were all at a cap gave up nothing; like a held cell it does not strain this step, so
    # its damage stays consistent with its unchanged column (PR #950 review).
    stuck = out.stuck
    eps1[stuck], warm1[stuck], d_eps[stuck] = eps0[stuck], np.minimum(healed, eps0)[stuck], 0.0
    at_floor = out.hc <= rheology.RIFT_CRITICAL_THICKNESS_M + 1e-6
    if p.oceanic_latch_on_break_only:
        at_floor = np.where(prep.continental, at_floor, at_floor & (out.hm <= lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M + 1e-6))
    # A cell breaks when the strain that thinned it took its crust to the floor, or at eps_break. There is no
    # interior eruption: a broken cell costs nothing to cross and waits for the split.
    pending.broken_now = broken | ((d_eps > 0.0) & ((eps1 >= p.eps_break - 1e-9) | at_floor))
    pending.out = out
    return pending


def _commit_step(world, plate, years: float, state: RiftRuntime, pending: _Pending) -> PlateStep:
    """Write a computed damage step to the plate and record it for the breakthrough test."""
    p, prep = state.params, pending.prep
    line, out, broken_now = prep.line, pending.out, pending.broken_now
    hc, new_elevation = prep.hc, plate.collect("elevation")
    if out is not None:
        hc, new_elevation = out.hc, _write_thinning(world, plate, prep, out, pending.d_eps)
    plate.set_fields_on_plate(rift_strain=pending.eps1, rift_warm_strain=pending.warm1, rift_broken=broken_now)
    step = PlateStep(plate.plate_id, plate.cell_keys.copy(), prep.points, prep.areas, prep.graph, prep.boundary, pending.eps1,
                     pending.d_eps, broken_now, ~prep.continental, prep.weak, prep.heated, new_elevation, hc, line)
    if line is not None:
        fresh = lines.LineProblem(**{**pending.prob.__dict__, "intact": ~broken_now[line.chain]})
        n_eff, load, _ = lines.line_state(fresh, pending.eps1, pending.warm1, years, p)
        step.pull_end = float(np.sum(load * fresh.ds))
        step.strength_end = float(np.sum((n_eff * fresh.ds)[fresh.intact]))
        step.broken_share = float(line.geometry.ds[broken_now[line.chain]].sum() / max(line.geometry.ds.sum(), 1e-9))
    state.steps[plate.plate_id] = step
    return step


@trace.traced("Finding rifts", trace.PLATES, user=True)
def damage_pass(world, years: float) -> None:
    """The damage step for every plate, once per world step, before any split is decided. Every plate is computed
    from the step-start world before any is written, so the result doesn't depend on plate order, and an error in
    any plate (a sub-step overrun, say) leaves the world untouched rather than partly strained (PR #950 review)."""
    state = runtime(world)
    pending = [(plate, _compute_step(world, plate, years, state)) for plate in world.plates]
    state.steps = {}
    for plate, step in pending:
        if step is not None:
            _commit_step(world, plate, years, state, step)


# --- breakthrough --------------------------------------------------------------------------------


def _tie_side_factory(world, step: PlateStep):
    """A line cell with as many edges to one daughter as to the other goes to the one whose fitted rotation best
    matches the mantle flow under it."""
    flow = mantle.world_flow_at(world, step.points)
    # `split_masks` updates the sides only between rounds, so every tie in a round shares them: fit the two poles once
    # per round, not once per tied cell (ties are common on a quad mesh, a sweep cost, #942).
    fitted = {"side": None, "poles": None}

    def poles(side):
        if fitted["side"] is None or not np.array_equal(fitted["side"], side):
            fitted["side"] = side.copy()
            fitted["poles"] = [mantle.fit_euler_pole(step.points[side == k], flow[side == k]) if (side == k).any() else None for k in (0, 1)]
        return fitted["poles"]

    def tie(side, cell) -> int:
        best, score = 0, -np.inf
        for k, pole in enumerate(poles(side)):
            if pole is None:
                continue
            value = float(np.cross(pole, step.points[cell]) @ flow[cell])
            if value > score:
                best, score = k, value
        return best

    return tie


def localization(step: PlateStep, p: RiftLineParams) -> dict:
    """The #699 failure mode at a split: the damaged share of the plate and its number of separate damaged patches."""
    damaged = step.eps >= p.eps_cand
    total = float(step.areas.sum())
    count = 0
    if damaged.any():
        _, labels = connected_components(step.graph[damaged][:, damaged], directed=False)
        count = int((np.bincount(labels, weights=step.areas[damaged]) > LOCALIZATION_COMPONENT_SHARE * total).sum())
    return {"damaged_share": float(step.areas[damaged].sum() / total), "damaged_components": count}


def near_path_mask(step: PlateStep, path: np.ndarray) -> np.ndarray:
    """The path cells and their edge neighbours: the band that is marked as belt and counted in the craton ledger."""
    mask = np.zeros(len(step.points), dtype=bool)
    mask[path] = True
    return mask | ((step.graph > 0).astype(np.int8) @ mask.astype(np.int8) > 0)


def _split_event(world, plate, step: PlateStep, mask_a, daughters, p: RiftLineParams) -> dict:
    path, ds = step.line.chain, step.line.geometry.ds
    share = lambda flag: float(ds[flag[path]].sum() / max(ds.sum(), 1e-9))  # noqa: E731
    geometry = step.line.geometry
    return {
        "kind": "split", "elapsed_years": world.elapsed_years, "plate_id": int(plate.plate_id),
        "daughters": [int(d.plate_id) for d in daughters], "path_cells": int(len(path)), "path_length_m": float(ds.sum()),
        "length_ratio": geometry.length / max(geometry.straight, 1.0), "oceanic_share": share(step.oceanic),
        "intact_share": share(~step.broken), "weak_share": share(step.weak), "plume_share": share(step.plumed),
        "open": step.line.is_open, "path_keys": step.keys[path].tolist(),
        "small_share": float(min(step.areas[mask_a].sum(), step.areas[~mask_a].sum()) / step.areas.sum()), **localization(step, p),
    }


def breakthrough(world, plate):
    """Step 6 for one plate: split it along its rift line once enough of the line is broken and the pull across it
    exceeds the strength left. The caller has applied the age and size gates. Returns the daughters or None."""
    state = runtime(world)
    step, p = state.steps.get(plate.plate_id), state.params
    if step is None or step.line is None or step.broken_share < p.f_break or not lines.ligament_gives(step.pull_end, step.strength_end):
        return None
    state.counters["break_tests"] += 1
    path = step.line.chain
    mask_a, why = lines.split_masks(step.graph, step.areas, path, p, _tie_side_factory(world, step), step.points)
    if mask_a is None:
        state.counters["split_rejected"] += 1
        state.counters["rejected_because"][why] = state.counters["rejected_because"].get(why, 0) + 1
        return None
    near_path = near_path_mask(step, path)
    cut_areas = cratons.band_cut_areas(world, plate, near_path)
    pieces = plate._plates_from_node_masks([mask_a, ~mask_a], [plate.plate_id, world.next_plate_id])
    if len(pieces) != 2:
        state.counters["split_rejected"] += 1
        state.counters["rejected_because"]["pieces"] = state.counters["rejected_because"].get("pieces", 0) + 1
        return None
    cratons.record_cut_areas(world, "split", *cut_areas)
    world.next_plate_id += 1
    flow = mantle.world_flow_at(world, step.points)
    for daughter, mask in zip(pieces, (mask_a, ~mask_a)):
        daughter.reset_age()
        daughter.set_omega(mantle.clamp_rate(mantle.fit_euler_pole(step.points[mask], flow[mask])))
        daughter.set_internal_stress(0.0)
        # The crust along the line was stretched to breaking: both margins are belt (RFC 0009).
        cratons.mark_belt(world, daughter, near_path[mask])
    state.events.append(_split_event(world, plate, step, mask_a, pieces, p))
    return pieces[0], pieces[1]


# --- margins -------------------------------------------------------------------------------------


def margin_donor_damage(plate, inserted_indices, rifted, rows: int):
    """The area-weighted strain of the pre-existing cells within `rows` hops of each newly inserted `rifted` margin
    cell: the cells `_allocate_stretch` thins to fill it. Returns (strain, warm strain) per rifted cell, or None when
    the plate carries no damage. Call before `_open_rift` changes the columns."""
    eps, warm = plate.collect("rift_strain"), plate.collect("rift_warm_strain")
    if not len(rifted) or not (eps.any() or warm.any()):
        return None
    points, _ = plate.all_points_and_elevation()
    areas = plate.node_areas_m2()
    pattern = (plate_graph(plate, points) > 0).astype(np.int8)
    pre = np.ones(len(points), dtype=np.int8)
    pre[inserted_indices] = 0
    reach = csr_matrix((np.ones(len(rifted), dtype=np.int8), (np.arange(len(rifted)), rifted)), shape=(len(rifted), len(points)))
    band = (reach @ pattern).tocsr()
    for _ in range(rows - 1):
        band = band + (band.multiply(pre[None, :]) @ pattern)
    band = band.multiply(pre[None, :]).tocsr()
    band.data[:] = 1
    weights = band.multiply(areas[None, :]).tocsr()
    total = np.asarray(weights.sum(axis=1)).ravel()
    safe = np.where(total > 0.0, total, 1.0)
    return (weights @ eps) / safe, (weights @ warm) / safe


def inherit_margin_damage(plate, rifted, stretch_share, donors) -> None:
    """Give each rifted margin cell its donors' strain times its stretched share: the stretched part of the cell is
    old rock and keeps its damage, the magmatic part is new and has none. A new margin cell is never broken."""
    if donors is None:
        return
    donor_eps, donor_warm = donors
    eps, warm, broken = plate.collect("rift_strain"), plate.collect("rift_warm_strain"), plate.collect("rift_broken")
    eps[rifted] = donor_eps * stretch_share
    warm[rifted] = np.minimum(donor_warm * stretch_share, eps[rifted])
    broken[rifted] = False
    plate.set_fields_on_plate(rift_strain=eps, rift_warm_strain=warm, rift_broken=broken)

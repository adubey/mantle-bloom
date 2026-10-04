"""What over-thickened continental crust does before any of it delaminates (issue #290).

Collision mostly shortens, thickens and moves continental crust around; it doesn't lose it.
Lower-crustal delamination is real, but only under specific conditions: the lower crust has
to sit deep enough to turn to eclogite, which makes it denser than the mantle beneath it,
and it has to be hot enough to detach. Reaching the numerical Hc cap is not one of those
conditions. This module holds the state the other passes need to tell those cases apart,
and the relief processes that act on standing orogens:

- **Thermal proxy.** `moho_temperature_c` is the steady-state conductive geotherm of a
  column with radiogenic crust over a non-radiogenic mantle lid. It is the only thermal
  state the model has. Thick crust heats itself, and a thick mantle lid keeps the Moho cool.
  A freshly thickened column is treated as already at steady state. That overstates how
  quickly young orogens weaken; a thermal-lag field is left to the anatexis phase.
- **Crust states.** `crust_state` sorts continental columns into cold-strong, hot-weak
  (ductile), melt-eligible, and delamination-eligible. They are diagnostics here. Of the
  thresholds behind them, only delamination eligibility and the ductile ramp change behaviour.
- **Collapse and ductile flow.** `relax_orogens` lets crust thicker than
  `COLLAPSE_ONSET_HC_M` spread down thickness gradients into neighbouring continental
  columns. Gravitational collapse and normal faulting act on any such column, hot or cold.
  Lower- and mid-crustal ductile flow adds to that as the Moho heats through the
  `DUCTILE_*_MOHO_C` ramp. Each transfer moves Hc and its continental-material provenance
  between cells. It is volume-conserving by cell area and subcycled so the result doesn't
  depend on the step size.
- **Delamination eligibility.** `delaminable_root_m` is the dense lower-crustal root a
  column could shed: crust below `ECLOGITE_DEPTH_M`, where the Moho is at least
  `DELAMINATION_MIN_MOHO_C`. `quad_tectonics._accrete_onto_survivors` lets a suture
  delaminate only from that root, at `DELAMINATION_RATE_PER_MYR`, and only after
  lateral spreading and tectonic escape have had first call.

`world.orogenic_relief_budget` books the attempted and completed volume of every process
(Hc volume, m^3). The continental-material ledger still books the terminal delamination sink.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from . import cratons, lithosphere, rheology
from .elevation_lines import effective_is_continental_from_codes

if TYPE_CHECKING:
    from .world import World

# --- Thermal proxy ----------------------------------------------------------------------
# Steady-state conduction through a crust with uniform heat production over a mantle lid
# with none, between a fixed surface temperature and a fixed temperature at the base of the
# lithosphere. 0.7 uW/m^3 is a bulk-continental-crust average (upper crust ~1-2, lower crust
# ~0.2-0.5). With these values a reference column (35 km crust, 100 km lid) has a surface heat
# flow of ~46 mW/m^2 and a Moho at ~480 C, typical of stable continent. Doubled crust on the
# same lid reaches ~950 C, and on a doubled lid ~860 C, the range inferred beneath Tibet.
SURFACE_TEMPERATURE_C = 10.0
LITHOSPHERE_BASE_TEMPERATURE_C = 1330.0
CRUST_CONDUCTIVITY_W_M_K = 2.5
CRUST_HEAT_PRODUCTION_W_M3 = 0.7e-6

# --- Crust-state thresholds -------------------------------------------------------------
# Lower crust weakens enough to flow from ~700 C (quartz- and feldspar-dominated ductile
# creep at geological strain rates); by ~900 C it is weak enough to flow freely, as a
# mid-/lower-crustal channel does.
DUCTILE_ONSET_MOHO_C = 700.0
DUCTILE_FULL_MOHO_C = 900.0
# Fluid-absent (muscovite / biotite dehydration) melting of metapelite and granitoid lower
# crust begins around 750-800 C.
MELT_ONSET_MOHO_C = 750.0
# Mafic lower crust converts to eclogite (garnet + omphacite, ~3.4-3.5 g/cm^3, denser than
# the ~3.3 mantle below it) from about 1.4 GPa, ~50 km depth. That needs the Moho hot enough
# for the reaction to run (it stalls in cold, dry granulite), and the root has to be thick
# enough to founder as a body rather than as a thin film.
ECLOGITE_DEPTH_M = 50_000.0
DELAMINATION_MIN_ROOT_M = 10_000.0
DELAMINATION_MIN_MOHO_C = 750.0
# At most this fraction of an eligible root may founder per Myr (e-folding ~20 Myr, the
# timescale of Rayleigh-Taylor removal of a thickened orogenic root).
DELAMINATION_RATE_PER_MYR = 0.05

CRUST_STATE_NOT_CONTINENTAL = 0
CRUST_STATE_COLD_STRONG = 1
CRUST_STATE_HOT_WEAK = 2
CRUST_STATE_MELT_ELIGIBLE = 3
CRUST_STATE_DELAMINATION_ELIGIBLE = 4
CRUST_STATE_NAMES = {
    CRUST_STATE_NOT_CONTINENTAL: "not_continental",
    CRUST_STATE_COLD_STRONG: "cold_strong",
    CRUST_STATE_HOT_WEAK: "hot_weak",
    CRUST_STATE_MELT_ELIGIBLE: "melt_eligible",
    CRUST_STATE_DELAMINATION_ELIGIBLE: "delamination_eligible",
}

# --- Collapse and ductile flow ----------------------------------------------------------
# Only crust above this thickness is mobile: ~2.7 km of isostatic relief on a reference lid,
# where orogenic plateaus start to collapse under their own weight (Tibet, the Altiplano and
# the Basin and Range all extend above ~3 km). Ordinary continent is never smoothed.
COLLAPSE_ONSET_HC_M = 50_000.0
# Thickness diffusivities (m^2/yr): flux per unit boundary length = K * thickness gradient.
# The ductile value is Tibet's lower-crustal channel: ~1 cm/yr through a ~15 km channel down
# a ~70->40 km thickness drop over ~500 km gives K ~ 2.5e3. Collapse and normal faulting is
# about a third of that: Tibet's ~1e-15/s E-W extension thins a 70 km plateau ~2 km/Myr.
COLLAPSE_DIFFUSIVITY_M2_PER_YR = 1_000.0
DUCTILE_DIFFUSIVITY_M2_PER_YR = 3_000.0
# Strong cratonic crust resists flowing like it resists rifting (cratons.CRATON_*_RESISTANCE).
CRATON_FLOW_RESISTANCE = 0.75
# Explicit diffusion is monotone below a Courant number of 0.5. 0.1 also holds the time-
# discretization error to well under 1% of the relief, so the outcome doesn't depend on
# how a span of time is split into steps.
FLOW_COURANT = 0.1
MAX_FLOW_SUBSTEPS = 512

BUDGET_ACCOUNTS = (
    # Suture accretion (quad_tectonics._accrete_onto_survivors), in the order it tries them.
    "suture_donated_m3",
    "suture_belt_placed_m3",
    "escape_attempted_m3",
    "escape_placed_m3",
    "delamination_attempted_m3",
    "delamination_completed_m3",
    "far_field_placed_m3",
    "foreland_spill_placed_m3",
    "overrider_placed_m3",
    "no_outlet_delaminated_m3",
    # Standing orogens (relax_orogens).
    "relief_mobile_excess_m3",
    "collapse_transferred_m3",
    "ductile_flow_transferred_m3",
)


def ensure_budget(world: "World") -> dict[str, float]:
    if not isinstance(getattr(world, "orogenic_relief_budget", None), dict):
        world.orogenic_relief_budget = {}
    for key in BUDGET_ACCOUNTS:
        world.orogenic_relief_budget.setdefault(key, 0.0)
    return world.orogenic_relief_budget


def record(world: "World | None", account: str, volume_m3: float) -> None:
    if world is None:
        return
    if account not in BUDGET_ACCOUNTS:
        raise KeyError(f"unknown orogenic relief account: {account}")
    volume_m3 = float(volume_m3)
    if not np.isfinite(volume_m3) or volume_m3 < 0.0:
        raise ValueError(f"orogenic relief volume must be finite and non-negative, got {volume_m3!r}")
    ensure_budget(world)[account] += volume_m3


def moho_temperature_c(hc_m: np.ndarray, hm_m: np.ndarray) -> np.ndarray:
    """Steady-state Moho temperature (C) of each column -- see the thermal-proxy constants.

    With heat production A in the crust (0 < z < Hc) only, T(z) = Ts + q0 z / k - A z^2 / 2k
    in the crust and linear below. Requiring T(Hc + Hm) = Tb fixes the surface heat flow:
    q0 = (k (Tb - Ts) + A Hc (Hc / 2 + Hm)) / (Hc + Hm)."""
    hc = np.maximum(np.asarray(hc_m, dtype=float), 0.0)
    hm = np.maximum(np.asarray(hm_m, dtype=float), 0.0)
    k, a = CRUST_CONDUCTIVITY_W_M_K, CRUST_HEAT_PRODUCTION_W_M3
    total = np.maximum(hc + hm, 1.0)
    q0 = (k * (LITHOSPHERE_BASE_TEMPERATURE_C - SURFACE_TEMPERATURE_C) + a * hc * (0.5 * hc + hm)) / total
    return SURFACE_TEMPERATURE_C + q0 * hc / k - a * hc**2 / (2.0 * k)


def ductile_weakness(moho_c: np.ndarray) -> np.ndarray:
    """0 (strong) .. 1 (freely flowing) lower crust, linear across the ductile ramp."""
    span = DUCTILE_FULL_MOHO_C - DUCTILE_ONSET_MOHO_C
    return np.clip((np.asarray(moho_c, dtype=float) - DUCTILE_ONSET_MOHO_C) / span, 0.0, 1.0)


def delaminable_root_m(hc_m: np.ndarray, hm_m: np.ndarray) -> np.ndarray:
    """Thickness of dense, eclogitized lower crust each column could shed: the crust below
    `ECLOGITE_DEPTH_M`, where that root is at least `DELAMINATION_MIN_ROOT_M` thick and the
    Moho at least `DELAMINATION_MIN_MOHO_C`; zero elsewhere."""
    hc = np.asarray(hc_m, dtype=float)
    root = np.maximum(hc - ECLOGITE_DEPTH_M, 0.0)
    eligible = (root >= DELAMINATION_MIN_ROOT_M) & (moho_temperature_c(hc, hm_m) >= DELAMINATION_MIN_MOHO_C)
    return np.where(eligible, root, 0.0)


def delamination_capacity_m(hc_m: np.ndarray, hm_m: np.ndarray, years_myr: float) -> np.ndarray:
    """How much of each column's eligible root may founder in `years_myr` -- the rate limit
    integrated exactly, so it doesn't depend on step size."""
    if years_myr <= 0.0:
        return np.zeros(np.shape(hc_m))
    return delaminable_root_m(hc_m, hm_m) * -np.expm1(-DELAMINATION_RATE_PER_MYR * years_myr)


def crust_state(hc_m: np.ndarray, hm_m: np.ndarray, continental: np.ndarray) -> np.ndarray:
    """Per column, one of the `CRUST_STATE_*` codes. The states nest -- a delamination-
    eligible column is also melt-eligible and hot -- so each column gets the strongest."""
    moho = moho_temperature_c(hc_m, hm_m)
    state = np.where(continental, CRUST_STATE_COLD_STRONG, CRUST_STATE_NOT_CONTINENTAL)
    state = np.where(continental & (moho >= DUCTILE_ONSET_MOHO_C), CRUST_STATE_HOT_WEAK, state)
    state = np.where(continental & (moho >= MELT_ONSET_MOHO_C), CRUST_STATE_MELT_ELIGIBLE, state)
    state = np.where(continental & (delaminable_root_m(hc_m, hm_m) > 0.0), CRUST_STATE_DELAMINATION_ELIGIBLE, state)
    return state.astype(np.int8)


def crust_state_areas_m2(world: "World") -> dict[str, float]:
    """World-wide area in each crust state (by accounting area), for diagnostics."""
    from .elevation_lines import line_spacing_rad

    spacing = line_spacing_rad(world.node_density)
    totals = {name: 0.0 for name in CRUST_STATE_NAMES.values()}
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        state = crust_state(plate.collect("crustal_thickness_m"), plate.collect("mantle_lithosphere_thickness_m"), continental)
        areas = plate.accounting_areas_m2(spacing)
        for code, name in CRUST_STATE_NAMES.items():
            totals[name] += float(areas[state == code].sum())
    return totals


def _continental_edges(plate, continental: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    graph = plate.adjacency()
    n = plate.node_count()
    rows = np.repeat(np.arange(n), np.diff(graph.offsets))
    cols = np.asarray(graph.neighbours)
    keep = (rows < cols) & continental[rows] & continental[cols]
    return rows[keep], cols[keep]


def relax_orogens(plate, world: "World", years: float) -> None:
    """Collapse and ductile flow on this quad plate's standing orogens -- see the module
    docstring. Moves Hc, `continental_material_m` in proportion, and craton crust (which
    becomes ordinary orogenic crust where it flows: craton ledger `collision_reworked_m3`)
    from thicker to thinner continental columns across shared cell edges. Hm stays: the
    crust flows over its mantle lid, decoupled from it."""
    if plate.node_count() == 0 or years <= 0.0 or not hasattr(plate, "adjacency"):
        return
    hc = plate.collect("crustal_thickness_m")
    if not np.any(hc > COLLAPSE_ONSET_HC_M):
        return
    codes = plate.collect("crust_type_code")
    continental = effective_is_continental_from_codes(codes, plate.crust_type == "continental")
    if not np.any(continental & (hc > COLLAPSE_ONSET_HC_M)):
        return
    i, j = _continental_edges(plate, continental)
    if not len(i):
        return

    hm = plate.collect("mantle_lithosphere_thickness_m")
    material = plate.collect("continental_material_m")
    craton = plate.collect("craton_crust_m")
    areas = plate.node_areas_m2()
    points = plate.surface_nodes().local_xyz
    # Flux across an edge per unit thickness difference: boundary length over centre
    # spacing. ~1 on a uniform lattice; the coarse side of a 2:1 edge shares only the fine
    # cell's edge, at three quarters of a coarse spacing.
    spacing_m = np.linalg.norm(points[i] - points[j], axis=1) * lithosphere.PLANET_RADIUS_M
    conductance = np.sqrt(np.minimum(areas[i], areas[j])) / np.maximum(spacing_m, 1.0)
    resistance = 1.0 - CRATON_FLOW_RESISTANCE * cratons.strength(craton)

    max_k = (COLLAPSE_DIFFUSIVITY_M2_PER_YR + DUCTILE_DIFFUSIVITY_M2_PER_YR) * resistance
    per_cell = np.bincount(i, conductance, len(hc)) + np.bincount(j, conductance, len(hc))
    rate = per_cell * max_k / areas
    stable_dt = FLOW_COURANT / max(float(rate.max()), 1e-300)
    substeps = int(min(MAX_FLOW_SUBSTEPS, max(1, np.ceil(years / stable_dt))))
    dt = years / substeps

    excess = np.maximum(hc - COLLAPSE_ONSET_HC_M, 0.0) * continental
    record(world, "relief_mobile_excess_m3", float(np.dot(excess, areas)))
    hc_start = hc.copy()
    craton_start = craton.copy()
    collapse_total = ductile_total = 0.0
    for _ in range(substeps):
        mobile_hc = np.maximum(hc, COLLAPSE_ONSET_HC_M)
        diff = mobile_hc[i] - mobile_hc[j]
        donor = np.where(diff > 0.0, i, j)
        receiver = np.where(diff > 0.0, j, i)
        collapse_k = COLLAPSE_DIFFUSIVITY_M2_PER_YR * resistance[donor]
        ductile_k = DUCTILE_DIFFUSIVITY_M2_PER_YR * ductile_weakness(moho_temperature_c(hc[donor], hm[donor])) * resistance[donor]
        flux = (collapse_k + ductile_k) * conductance * np.abs(diff) * dt
        # Never draw a donor below the onset in one substep, however its edges add up.
        out = np.bincount(donor, flux, len(hc))
        available = np.maximum(hc - COLLAPSE_ONSET_HC_M, 0.0) * areas
        scale = np.divide(available, out, out=np.ones(len(hc)), where=out > available)
        flux = flux * scale[donor]
        if not np.any(flux > 0.0):
            break
        thickness = flux / areas[donor]
        material_moved = flux * np.divide(material[donor], hc[donor], out=np.zeros(len(flux)), where=hc[donor] > 0.0)
        craton_moved = flux * np.divide(craton[donor], hc[donor], out=np.zeros(len(flux)), where=hc[donor] > 0.0)
        hc -= np.bincount(donor, thickness, len(hc))
        hc += np.bincount(receiver, flux, len(hc)) / areas
        material -= np.bincount(donor, material_moved, len(hc)) / areas
        material += np.bincount(receiver, material_moved, len(hc)) / areas
        craton -= np.bincount(donor, craton_moved, len(hc)) / areas
        ductile_share = np.divide(ductile_k, collapse_k + ductile_k, out=np.zeros(len(flux)), where=(collapse_k + ductile_k) > 0.0)
        ductile_total += float(np.dot(flux, ductile_share))
        collapse_total += float(np.dot(flux, 1.0 - ductile_share))
    # Each transfer carries its donor's own material and craton fractions, so neither can
    # exceed its column; these clips only absorb round-off.
    material = np.maximum(material, 0.0)
    craton = np.maximum(craton, 0.0)
    record(world, "collapse_transferred_m3", collapse_total)
    record(world, "ductile_flow_transferred_m3", ductile_total)
    cratons.record(world, "collision_reworked_m3", float(np.dot(np.maximum(craton_start - craton, 0.0), areas)))

    density = lithosphere.node_crust_density(codes, plate.crust_type)
    elevation = plate.collect("elevation")
    shift = lithosphere.isostatic_elevation(hc, hm, density) - lithosphere.isostatic_elevation(hc_start, hm, density)
    plate.set_fields_on_plate(
        crustal_thickness_m=hc,
        continental_material_m=material,
        craton_crust_m=craton,
        elevation=rheology.clip_elevation_bounds(elevation + shift),
    )

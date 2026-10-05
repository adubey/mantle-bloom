"""What over-thickened continental crust does before any of it delaminates (issue #290).

Collision mostly shortens, thickens and moves continental crust around; it doesn't lose it.
Lower-crustal delamination is real, but only under specific conditions: the lower crust has
to sit deep enough to turn to eclogite, which makes it denser than the mantle beneath it,
and it has to be hot enough to detach. Reaching the numerical Hc cap is not one of those
conditions. This module holds the state the other passes need to tell those cases apart,
and the relief processes that act on standing orogens:

- **Thermal proxy.** `moho_temperature_c` is the steady-state conductive geotherm of a
  column with radiogenic crust over a non-radiogenic mantle lid, less the column's thermal
  lag (`moho_thermal_lag_c`). Thick crust heats itself, and a thick mantle lid keeps the
  Moho cool. Thickening buries the Moho with its old temperature (`bury_moho`), so freshly
  thickened crust starts cold and the lag decays toward steady state over
  `THERMAL_RELAXATION_MYR` (`relax_thermal_lag`).
- **Crust states.** `crust_state` sorts continental columns into cold-strong, hot-weak
  (ductile), melt-eligible, and delamination-eligible.
- **Anatexis.** `anatexis` partially melts the fertile lower crust below the
  `MELT_ONSET_MOHO_C` isotherm, rate-limited. The melt rises: most of it is emplaced in the
  column's own upper crust, and `ANATEXIS_RELAMINATION_FRACTION` of it relaminates into
  thinner neighbouring columns. The refractory residue stays at the base of the column as
  `restite_m`, which can't melt again.
- **Collapse and ductile flow.** `relax_orogens` lets crust thicker than
  `COLLAPSE_ONSET_HC_M` spread down thickness gradients into neighbouring continental
  columns. Gravitational collapse and normal faulting act on any such column, hot or cold.
  Lower- and mid-crustal ductile flow adds to that as the Moho heats through the
  `DUCTILE_*_MOHO_C` ramp. Each transfer moves Hc and its continental-material provenance
  between cells. It is volume-conserving by cell area and subcycled so the result doesn't
  depend on the step size.
- **Delamination eligibility.** `delaminable_root_m` is the dense lower-crustal root a
  column could shed: the crust below `ECLOGITE_DEPTH_M` or its restite, whichever reaches
  higher, where the Moho is at least `DELAMINATION_MIN_MOHO_C`.
  `quad_tectonics._accrete_onto_survivors` lets a suture delaminate only from that root,
  restite first, at `DELAMINATION_RATE_PER_MYR`, and only after lateral spreading and
  tectonic escape have had first call.

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
# Thermal lag. Thickened crust takes ~20-40 Myr to heat to its new steady state (England &
# Thompson 1984): the diffusion time across a ~70 km orogenic crust, L^2 / (pi^2 kappa) with
# kappa = 1e-6 m^2/s, is ~16 Myr, and its own radiogenic heating sets the slower tail. One
# e-folding time stands in for both.
THERMAL_RELAXATION_MYR = 25.0

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

# --- Anatexis ---------------------------------------------------------------------------
# Fluid-absent melting of fertile lower crust yields ~10-20% melt by ~850 C and ~30-40% by
# ~950 C (Clemens & Vielzeuf 1987; Vielzeuf & Holloway 1988). The melt fraction ramps
# linearly in temperature from the onset up to `MAX_MELT_FRACTION` at `MELT_FULL_C`.
MELT_FULL_C = 950.0
MAX_MELT_FRACTION = 0.35
# Fertile crust above the solidus melts and gives up its melt at this rate per Myr (e-folding
# ~10 Myr: segregation and ascent of crustal melt take ~1-10 Myr once it is connected, and
# heating the next increment of crust past the solidus is slower).
ANATEXIS_RATE_PER_MYR = 0.1
# Most granite pools in the overlying upper crust (Himalayan leucogranites sit almost straight
# above their source), but some moves out laterally along the mid-crust before it stalls.
ANATEXIS_RELAMINATION_FRACTION = 0.3

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
    # Convergent shortening past the Hc ceiling (quad_tectonics._place_ceiling_overflow):
    # overflow = residue + melt placed + melt with no outlet.
    "ceiling_overflow_m3",
    "ceiling_overflow_residue_m3",
    "ceiling_overflow_melt_placed_m3",
    "ceiling_overflow_no_outlet_m3",
    # Anatexis (anatexis): extracted = emplaced + relaminated.
    "anatexis_melt_extracted_m3",
    "anatexis_melt_emplaced_m3",
    "anatexis_melt_relaminated_m3",
    "anatexis_residue_m3",
    # Of all delaminated suture roots, the share that was restite.
    "restite_delaminated_m3",
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


def steady_moho_temperature_c(hc_m: np.ndarray, hm_m: np.ndarray) -> np.ndarray:
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


def moho_temperature_c(hc_m: np.ndarray, hm_m: np.ndarray, lag_c: np.ndarray | None = None) -> np.ndarray:
    """Moho temperature (C) of each column: the steady state less its thermal lag `lag_c`
    (`moho_thermal_lag_c`; none means at steady state), kept between the surface and
    lithosphere-base temperatures."""
    steady = steady_moho_temperature_c(hc_m, hm_m)
    if lag_c is None:
        return steady
    return np.clip(steady - np.asarray(lag_c, dtype=float), SURFACE_TEMPERATURE_C, LITHOSPHERE_BASE_TEMPERATURE_C)


def bury_moho(
    lag_c: np.ndarray, hc_before: np.ndarray, hm_before: np.ndarray, hc_after: np.ndarray, hm_after: np.ndarray
) -> np.ndarray:
    """New thermal lag of columns whose crust or lid just changed thickness tectonically. The
    Moho keeps its temperature as shortening or stacking buries it -- rock carries its heat
    with it -- so the lag grows by however much the steady state rose. Thinning works the
    other way: an exhumed Moho is briefly hotter than its new steady state (a negative lag)."""
    steady = steady_moho_temperature_c(hc_after, hm_after)
    lag = np.asarray(lag_c, dtype=float) + steady - steady_moho_temperature_c(hc_before, hm_before)
    return np.clip(lag, steady - LITHOSPHERE_BASE_TEMPERATURE_C, steady - SURFACE_TEMPERATURE_C)


def relax_thermal_lag(lag_c: np.ndarray, years_myr: float) -> np.ndarray:
    """The lag after `years_myr` of conductive re-equilibration -- exact, so step-size free."""
    return np.asarray(lag_c, dtype=float) * np.exp(-max(years_myr, 0.0) / THERMAL_RELAXATION_MYR)


def ductile_weakness(moho_c: np.ndarray) -> np.ndarray:
    """0 (strong) .. 1 (freely flowing) lower crust, linear across the ductile ramp."""
    span = DUCTILE_FULL_MOHO_C - DUCTILE_ONSET_MOHO_C
    return np.clip((np.asarray(moho_c, dtype=float) - DUCTILE_ONSET_MOHO_C) / span, 0.0, 1.0)


def delaminable_root_m(
    hc_m: np.ndarray, hm_m: np.ndarray, lag_c: np.ndarray | None = None, restite_m: np.ndarray | None = None
) -> np.ndarray:
    """Thickness of dense lower crust each column could shed: the eclogitized crust below
    `ECLOGITE_DEPTH_M` or the column's garnet-rich restite (`restite_m`, at its base),
    whichever reaches higher, where that root is at least `DELAMINATION_MIN_ROOT_M` thick and
    the Moho at least `DELAMINATION_MIN_MOHO_C`; zero elsewhere."""
    hc = np.asarray(hc_m, dtype=float)
    root = np.maximum(hc - ECLOGITE_DEPTH_M, 0.0)
    if restite_m is not None:
        root = np.maximum(root, np.clip(np.asarray(restite_m, dtype=float), 0.0, np.maximum(hc, 0.0)))
    eligible = (root >= DELAMINATION_MIN_ROOT_M) & (moho_temperature_c(hc, hm_m, lag_c) >= DELAMINATION_MIN_MOHO_C)
    return np.where(eligible, root, 0.0)


def delamination_capacity_m(
    hc_m: np.ndarray,
    hm_m: np.ndarray,
    years_myr: float,
    lag_c: np.ndarray | None = None,
    restite_m: np.ndarray | None = None,
) -> np.ndarray:
    """How much of each column's eligible root may founder in `years_myr` -- the rate limit
    integrated exactly, so it doesn't depend on step size."""
    if years_myr <= 0.0:
        return np.zeros(np.shape(hc_m))
    return delaminable_root_m(hc_m, hm_m, lag_c, restite_m) * -np.expm1(-DELAMINATION_RATE_PER_MYR * years_myr)


def plate_delamination_capacity_m3(plate, continental: np.ndarray, years: float) -> np.ndarray:
    """Per cell of a quad plate, the root volume (m^3) it may shed over `years` (zero off
    `continental`), from its own thermal lag and restite."""
    capacity = delamination_capacity_m(
        plate.collect("crustal_thickness_m"),
        plate.collect("mantle_lithosphere_thickness_m"),
        years / 1_000_000.0,
        plate.collect("moho_thermal_lag_c"),
        plate.collect("restite_m"),
    )
    return np.where(continental, capacity * plate.node_areas_m2(), 0.0)


def crust_state(
    hc_m: np.ndarray,
    hm_m: np.ndarray,
    continental: np.ndarray,
    lag_c: np.ndarray | None = None,
    restite_m: np.ndarray | None = None,
) -> np.ndarray:
    """Per column, one of the `CRUST_STATE_*` codes. The states nest -- a delamination-
    eligible column is also melt-eligible and hot -- so each column gets the strongest."""
    moho = moho_temperature_c(hc_m, hm_m, lag_c)
    state = np.where(continental, CRUST_STATE_COLD_STRONG, CRUST_STATE_NOT_CONTINENTAL)
    state = np.where(continental & (moho >= DUCTILE_ONSET_MOHO_C), CRUST_STATE_HOT_WEAK, state)
    state = np.where(continental & (moho >= MELT_ONSET_MOHO_C), CRUST_STATE_MELT_ELIGIBLE, state)
    state = np.where(
        continental & (delaminable_root_m(hc_m, hm_m, lag_c, restite_m) > 0.0), CRUST_STATE_DELAMINATION_ELIGIBLE, state
    )
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
        state = crust_state(
            plate.collect("crustal_thickness_m"),
            plate.collect("mantle_lithosphere_thickness_m"),
            continental,
            plate.collect("moho_thermal_lag_c"),
            plate.collect("restite_m"),
        )
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


def _edge_conductance(plate, i: np.ndarray, j: np.ndarray, areas: np.ndarray) -> np.ndarray:
    """Flux across each shared edge (i, j) per unit thickness difference: boundary length over
    centre spacing. ~1 on a uniform lattice; the coarse side of a 2:1 edge shares only the fine
    cell's edge, at three quarters of a coarse spacing."""
    points = plate.surface_nodes().local_xyz
    spacing_m = np.linalg.norm(points[i] - points[j], axis=1) * lithosphere.PLANET_RADIUS_M
    return np.sqrt(np.minimum(areas[i], areas[j])) / np.maximum(spacing_m, 1.0)


def melt_fraction(temperature_c: np.ndarray) -> np.ndarray:
    """Equilibrium melt fraction of fertile crust at `temperature_c`."""
    span = MELT_FULL_C - MELT_ONSET_MOHO_C
    return MAX_MELT_FRACTION * np.clip((np.asarray(temperature_c, dtype=float) - MELT_ONSET_MOHO_C) / span, 0.0, 1.0)


def _mean_melt_fraction(t_low: np.ndarray, t_high: np.ndarray) -> np.ndarray:
    """Mean of `melt_fraction` over temperatures evenly spread from `t_low` to `t_high` -- the
    ramp's integral over the interval, divided by its width."""
    span = MELT_FULL_C - MELT_ONSET_MOHO_C

    def integral(t):  # of melt_fraction / MAX_MELT_FRACTION from the onset up to t
        x = np.clip((t - MELT_ONSET_MOHO_C) / span, 0.0, None)
        return span * np.where(x <= 1.0, 0.5 * x**2, x - 0.5)

    width = t_high - t_low
    mean = np.divide(integral(t_high) - integral(t_low), width, out=np.zeros(np.shape(width)), where=width > 1e-9)
    return MAX_MELT_FRACTION * np.where(width > 1e-9, mean, melt_fraction(t_high) / MAX_MELT_FRACTION)


def anatexis_yield_m(
    hc_m: np.ndarray, hm_m: np.ndarray, lag_c: np.ndarray, restite_m: np.ndarray, years_myr: float
) -> tuple[np.ndarray, np.ndarray]:
    """(melt, residue) thickness each column's lower crust yields over `years_myr`.

    The crust hotter than `MELT_ONSET_MOHO_C` -- below that isotherm, with temperature taken
    as linear in depth down to the (lagged) Moho -- is the melting zone. Its bottom
    `restite_m` is already depleted and can't melt again; the fertile rest above it melts at
    `ANATEXIS_RATE_PER_MYR` (integrated exactly), each increment giving up its mean
    `melt_fraction` and leaving the rest as residue."""
    hc = np.maximum(np.asarray(hc_m, dtype=float), 0.0)
    moho = moho_temperature_c(hc, hm_m, lag_c)
    hot = moho > MELT_ONSET_MOHO_C
    zone = np.where(hot, hc * (moho - MELT_ONSET_MOHO_C) / np.maximum(moho - SURFACE_TEMPERATURE_C, 1e-9), 0.0)
    fertile = np.maximum(zone - np.clip(restite_m, 0.0, None), 0.0)
    fertile_base_c = MELT_ONSET_MOHO_C + (moho - MELT_ONSET_MOHO_C) * np.divide(
        fertile, zone, out=np.zeros(len(zone)), where=zone > 0.0
    )
    processed = fertile * -np.expm1(-ANATEXIS_RATE_PER_MYR * max(years_myr, 0.0))
    melt = processed * _mean_melt_fraction(np.full(len(zone), MELT_ONSET_MOHO_C), fertile_base_c)
    return melt, processed - melt


def anatexis(plate, world: "World", years: float) -> None:
    """Partial melting of this quad plate's hot continental lower crust -- see the module
    docstring. Adds the residue to `restite_m` and moves `ANATEXIS_RELAMINATION_FRACTION` of
    the melt, with its column's continental-material fraction, into thinner continental
    neighbours across shared edges; the rest of the melt stays in its own column's upper
    crust, which changes no thickness. Relaminated melt that would push a receiver past the
    Hc cap stays home instead. Cratonic crust that melts out becomes ordinary crust (craton
    ledger `collision_reworked_m3`). The intruding melt brings its heat, so it doesn't touch
    the thermal lag."""
    if plate.node_count() == 0 or years <= 0.0 or not hasattr(plate, "adjacency"):
        return
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    lag = plate.collect("moho_thermal_lag_c")
    codes = plate.collect("crust_type_code")
    continental = effective_is_continental_from_codes(codes, plate.crust_type == "continental")
    restite = np.where(continental, np.clip(plate.collect("restite_m"), 0.0, hc), 0.0)
    melt, residue = anatexis_yield_m(hc, hm, lag, restite, years / 1_000_000.0)
    melt = np.where(continental, melt, 0.0)
    residue = np.where(continental, residue, 0.0)
    if not np.any(melt > 0.0):
        plate.set_fields_on_plate(restite_m=restite)
        return
    areas = plate.node_areas_m2()
    restite = restite + residue
    record(world, "anatexis_melt_extracted_m3", float(melt @ areas))
    record(world, "anatexis_residue_m3", float(residue @ areas))

    i, j = _continental_edges(plate, continental)
    relaminated = 0.0
    if len(i):
        material = plate.collect("continental_material_m")
        craton = plate.collect("craton_crust_m")
        donor = np.where(hc[i] >= hc[j], i, j)
        receiver = np.where(hc[i] >= hc[j], j, i)
        downhill = hc[donor] > hc[receiver]
        donor, receiver = donor[downhill], receiver[downhill]
        weight = _edge_conductance(plate, donor, receiver, areas)
        share = np.divide(weight, np.bincount(donor, weight, len(hc))[donor], out=np.zeros(len(weight)), where=weight > 0.0)
        flux = ANATEXIS_RELAMINATION_FRACTION * melt[donor] * areas[donor] * share
        room = np.maximum(lithosphere.MAX_CRUSTAL_THICKNESS_M - hc, 0.0) * areas
        incoming = np.bincount(receiver, flux, len(hc))
        flux = flux * np.divide(room, incoming, out=np.ones(len(hc)), where=incoming > room)[receiver]
        if np.any(flux > 0.0):
            hc_start = hc.copy()
            out = np.bincount(donor, flux, len(hc))
            into = np.bincount(receiver, flux, len(hc))
            material_fraction = np.divide(material, hc, out=np.zeros(len(hc)), where=hc > 0.0)
            craton_fraction = np.divide(craton, hc, out=np.zeros(len(hc)), where=hc > 0.0)
            material_moved = np.bincount(receiver, flux * material_fraction[donor], len(hc))
            craton_out = out * craton_fraction
            hc = hc + (into - out) / areas
            material = np.maximum(material + (material_moved - out * material_fraction) / areas, 0.0)
            craton = np.maximum(craton - craton_out / areas, 0.0)
            restite = np.minimum(restite, hc)
            relaminated = float(out.sum())
            cratons.record(world, "collision_reworked_m3", float(craton_out.sum()))
            density = lithosphere.node_crust_density(codes, plate.crust_type)
            shift = lithosphere.isostatic_elevation(hc, hm, density) - lithosphere.isostatic_elevation(hc_start, hm, density)
            plate.set_fields_on_plate(
                crustal_thickness_m=hc,
                continental_material_m=material,
                craton_crust_m=craton,
                elevation=rheology.clip_elevation_bounds(plate.collect("elevation") + shift),
            )
    record(world, "anatexis_melt_relaminated_m3", relaminated)
    record(world, "anatexis_melt_emplaced_m3", max(float(melt @ areas) - relaminated, 0.0))
    plate.set_fields_on_plate(restite_m=restite)


def evolve_standing_orogens(plate, world: "World", years: float) -> None:
    """One step of this quad plate's standing orogens, after its tectonics: the Moho lag
    relaxes toward steady state, hot lower crust melts (`anatexis`), and thick crust
    collapses and flows (`relax_orogens`)."""
    if plate.node_count() == 0 or years <= 0.0 or not hasattr(plate, "adjacency"):
        return
    plate.set_fields_on_plate(moho_thermal_lag_c=relax_thermal_lag(plate.collect("moho_thermal_lag_c"), years / 1_000_000.0))
    anatexis(plate, world, years)
    relax_orogens(plate, world, years)


def relax_orogens(plate, world: "World", years: float) -> None:
    """Collapse and ductile flow on this quad plate's standing orogens -- see the module
    docstring. Moves Hc, `continental_material_m` and the mobile cover (mobile_cover.py) in
    proportion, and craton crust (which becomes ordinary orogenic crust where it flows: craton
    ledger `collision_reworked_m3`) from thicker to thinner continental columns across shared
    cell edges. Hm stays: the crust flows over its mantle lid, decoupled from it."""
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
    lag = plate.collect("moho_thermal_lag_c")
    material = plate.collect("continental_material_m")
    craton = plate.collect("craton_crust_m")
    cover = plate.collect("mobile_cover_m")
    cover_material = plate.collect("mobile_cover_continental_m")
    areas = plate.node_areas_m2()
    conductance = _edge_conductance(plate, i, j, areas)
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
        ductile_k = DUCTILE_DIFFUSIVITY_M2_PER_YR * ductile_weakness(moho_temperature_c(hc[donor], hm[donor], lag[donor])) * resistance[donor]
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
        fraction = np.divide(flux, hc[donor], out=np.zeros(len(flux)), where=hc[donor] > 0.0)
        for layer in (cover, cover_material):
            moved = fraction * layer[donor]
            layer -= np.bincount(donor, moved, len(hc)) / areas
            layer += np.bincount(receiver, moved, len(hc)) / areas
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
    cover = np.maximum(cover, 0.0)
    cover_material = np.clip(cover_material, 0.0, np.minimum(cover, material))
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
        mobile_cover_m=cover,
        mobile_cover_continental_m=cover_material,
        # Refractory restite stays at the base, like Hm; only a column thinned below it
        # loses any.
        restite_m=np.minimum(plate.collect("restite_m"), hc),
        elevation=rheology.clip_elevation_bounds(elevation + shift),
    )

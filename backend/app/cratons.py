"""Persistent cratons: old, stable continental cores that resist ordinary tectonics (issue #274).

State. Three per-node surface fields carry everything (see surface_fields.py):

- ``craton_crust_m`` -- the cratonic share of the column's crust, a thickness in metres and
  never more than Hc. Extensive, so every topology change (quad merge, partition, defrag,
  legacy conversion) conserves its area-integrated volume exactly as it does Hc's, instead of
  resetting or averaging a flag away. Crust added on top later (sediment, accreted belts,
  arc magmatism) is ordinary continental crust: the cratonic volume only grows by formation.
- ``craton_formed_years`` -- provenance: the ``elapsed_years`` at which the craton's column
  stabilised (negative for cratons seeded at generation, which predate the simulation;
  surface_fields.CRATON_UNFORMED_YEARS where there is no craton). A HISTORY field, so a merged
  or coarsened cell keeps its oldest contributor's date.
- ``stable_continental_myr`` -- the formation clock: how long the column has sat as quiet
  continental interior.

Formation rule. A column is quiet when it is genuine continental crust (Hc within
[CRATON_HOST_MIN_HC_M, CRATON_HOST_MAX_HC_M] -- neither a drowned shelf nor an active orogen)
at least CRATON_FORMATION_MARGIN_KM from its plate's edge or any non-genuine cell. Once the
clock reaches CRATON_FORMATION_MYR its whole Hc becomes cratonic. World generation (and the
first load of an older save) seeds the deep interiors -- CRATON_SEED_MARGIN_KM in -- as
cratons that already exist, dated CRATON_SEED_AGE_YEARS before the start.

Resistance. Strength rises linearly with cratonic thickness up to CRATON_FULL_STRENGTH_HC_M.
It scales back erosional Hc removal (erosion.py), divergent thinning
(lithosphere_plate.deform_columns), how much of a rift's stretched footprint a craton donates
(quad_tectonics._allocate_stretch), and delays boundary consumption of a craton cell until it
has been overlapped for CRATON_RETREAT_DELAY_YEARS (quad_tectonics._retreat), so shortening
goes into the younger belts around it first.

Destruction. Every loss of cratonic volume is booked into ``World.craton_ledger`` under the
mechanism that caused it (CRATON_SINK_ACCOUNTS); ``balance_error_m3`` checks that sources,
sinks and the live inventory agree. Boundary consumption also books the continental material
it removes into the continental-material ledger (continental_ledger.py).
Line-backed plates carry the fields through every operation but neither seed nor form
cratons and get no retreat resistance: that surface is being retired (issue #251).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from . import lithosphere
from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from .mantle import PLANET_RADIUS_KM
from .surface_fields import CRATON_UNFORMED_YEARS

if TYPE_CHECKING:
    from .plates import Plate
    from .world import World

# --- Formation ---------------------------------------------------------------------------

# The Hc window a column must sit in to count as stable continental interior: at least 80%
# of the continental reference (not a drowned shelf or thinned rift flank) and below the
# thickened crust of a still-active orogen.
CRATON_HOST_MIN_HC_M = 0.8 * lithosphere.REFERENCE_HC_CONTINENTAL_M
CRATON_HOST_MAX_HC_M = 50_000.0
# How far inside its margin a column must be to keep its formation clock running.
CRATON_FORMATION_MARGIN_KM = 500.0
# How long a column must stay quiet to cratonise. Real cratonic roots stabilise over a few
# hundred Myr; long enough here that an ordinary collision belt has to stop deforming and
# relax back into the host Hc window before it joins the craton.
CRATON_FORMATION_MYR = 200.0
# Generation seeding: the deep interiors of the starting continents. 875 km (seven default
# cells) covers ~40-50% of genuine continental area on default worlds -- the order of
# Earth's shields plus platforms.
CRATON_SEED_MARGIN_KM = 875.0
CRATON_SEED_AGE_YEARS = 2.5e9

# --- Resistance --------------------------------------------------------------------------

# Cratonic thickness at which resistance is at full strength.
CRATON_FULL_STRENGTH_HC_M = 25_000.0
# Fraction of erosional Hc removal a full-strength craton column resists.
CRATON_EROSION_RESISTANCE = 0.75
# Fraction of divergent thinning a full-strength craton column resists.
CRATON_RIFT_RESISTANCE = 0.75
# Fraction of its share of a rift's stretched footprint a full-strength craton donor
# withholds -- the rest of that footprint is fresh magmatic crust instead.
CRATON_STRETCH_RESISTANCE = 0.75
# How long a full-strength craton cell must sit overlapped by another plate before boundary
# retreat may consume it. Linear in strength. Longer than merge_split's
# FORCED_MERGE_SUSTAINED_YEARS, so a collision jammed against a craton gets the chance to end
# in a suture merge rather than consuming it. Measured over 300 Myr on three seeds
# (bin/debug/compare_cratons.py): 20 Myr still let collisions rework about half the seeded
# cratons; 60 Myr cut that by 2-3x, and never consuming them kept no more continent than
# 60 Myr did while making deep craton subduction impossible.
CRATON_RETREAT_DELAY_YEARS = 60_000_000.0

# --- Ledger ------------------------------------------------------------------------------

CRATON_SOURCE_ACCOUNTS = ("initial_m3", "formed_m3")
# Each declared destruction mechanism. "collision_reworked" is cratonic crust consumed at a
# suture or reworked by faulting: the material stays on the surface as ordinary continental
# crust. "topology_removed" is craton on a plate or fragment removed outright by topology
# cleanup. "unattributed" is signed and should stay ~0 -- anything there escaped every
# instrumented site.
CRATON_SINK_ACCOUNTS = (
    "rifted_m3",
    "subducted_m3",
    "delaminated_m3",
    "collision_reworked_m3",
    "eroded_m3",
    "topology_removed_m3",
    "unattributed_m3",
)
CRATON_ACCOUNTS = CRATON_SOURCE_ACCOUNTS + CRATON_SINK_ACCOUNTS
# A single step's destruction above this fraction of the live craton volume is logged to the
# event console, naming the mechanism.
CRATON_EVENT_FRACTION = 0.01


def strength(craton_crust_m: np.ndarray) -> np.ndarray:
    """0..1 resistance per node, linear in cratonic thickness."""
    return np.clip(np.asarray(craton_crust_m, dtype=float) / CRATON_FULL_STRENGTH_HC_M, 0.0, 1.0)


def _hops(km: float, spacing_rad: float) -> int:
    return max(1, round(km / (spacing_rad * PLANET_RADIUS_KM)))


def _genuine_continental(plate: "Plate") -> np.ndarray:
    hc = plate.collect("crustal_thickness_m")
    continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
    return continental & (hc >= CRATON_HOST_MIN_HC_M)


def _interior_hops(plate: "Plate", host: np.ndarray, max_hops: int) -> np.ndarray | None:
    """Per cell, edge hops to the nearest plate edge or non-`host` cell, saturating at
    `max_hops + 1`. None on surfaces without cell adjacency (line plates)."""
    probe = getattr(plate, "_probe_neighbour_indices", None)
    if probe is None or plate.node_count() == 0:
        return None
    from .quad_tectonics import hop_distance

    edge = np.any(np.any(probe() < 0, axis=2), axis=1)
    return hop_distance(plate, edge | ~host, max_hops)


# --- Ledger bookkeeping ------------------------------------------------------------------


def ensure_ledger(world: "World") -> None:
    """Backfill the ledger's accounts without touching any plate -- safe on load."""
    if not isinstance(getattr(world, "craton_ledger", None), dict):
        world.craton_ledger = {}
    for key in CRATON_ACCOUNTS:
        if key != "initial_m3":
            world.craton_ledger.setdefault(key, 0.0)


def ensure_initialized(world: "World") -> None:
    """`ensure_ledger`, then -- the first time a world with cell-surface plates is generated
    or stepped -- seed its starting cratons. The `initial_m3` key's presence marks a seeded
    world, so a world with no eligible interior is never reseeded, while a line world stays
    unmarked so the quad world a legacy save converts into is still seeded. Loading a save
    never seeds: an older save gets its cratons on its first step."""
    ensure_ledger(world)
    if "initial_m3" not in world.craton_ledger and any(hasattr(p, "_probe_neighbour_indices") for p in world.plates):
        world.craton_ledger["initial_m3"] = seed_initial_cratons(world)


def record(world: "World", account: str, volume_m3: float) -> None:
    """Book `volume_m3` of cratonic crust into one account. Only `unattributed_m3` may go
    negative (a remap that gained volume)."""
    if account not in CRATON_ACCOUNTS:
        raise KeyError(f"unknown craton ledger account: {account}")
    volume_m3 = float(volume_m3)
    if not np.isfinite(volume_m3) or (volume_m3 < 0.0 and account != "unattributed_m3"):
        raise ValueError(f"craton ledger volume must be finite and non-negative, got {volume_m3!r}")
    if volume_m3 == 0.0:
        return
    ensure_ledger(world)
    world.craton_ledger[account] = world.craton_ledger.get(account, 0.0) + volume_m3


def _areas(world: "World", plate: "Plate") -> np.ndarray:
    return plate.accounting_areas_m2(line_spacing_rad(world.node_density))


def plate_volume_m3(world: "World", plate: "Plate") -> float:
    if plate.node_count() == 0:
        return 0.0
    return float(np.dot(plate.collect("craton_crust_m"), _areas(world, plate)))


def live_volume_m3(world: "World") -> float:
    return sum(plate_volume_m3(world, plate) for plate in world.plates)


def booked_sinks_m3(world: "World") -> float:
    ledger = getattr(world, "craton_ledger", {})
    return sum(float(ledger.get(key, 0.0)) for key in CRATON_SINK_ACCOUNTS)


def balance_error_m3(world: "World") -> float:
    """Live + destroyed - (seeded + formed); zero when every change was booked."""
    ledger = world.craton_ledger
    sources = sum(float(ledger.get(key, 0.0)) for key in CRATON_SOURCE_ACCOUNTS)
    return live_volume_m3(world) + booked_sinks_m3(world) - sources


def clip_to_column(world: "World", plate: "Plate", account: str) -> float:
    """Hold `craton_crust_m` to [0, Hc] and to continental hosts, booking what that removes
    to `account`. The clip semantics suit processes that strip crust from the top (erosion,
    fault relief): cover goes first, the craton only once Hc drops below it."""
    if plate.node_count() == 0:
        return 0.0
    craton = plate.collect("craton_crust_m")
    if not np.any(craton > 0.0):
        return 0.0
    hc = plate.collect("crustal_thickness_m")
    continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
    clipped = np.where(continental, np.clip(craton, 0.0, np.maximum(hc, 0.0)), 0.0)
    lost = float(np.dot(craton - clipped, _areas(world, plate)))
    if lost != 0.0:
        _write(plate, clipped)
        record(world, account, lost)
    return lost


def thin_with_column(world: "World", plate: "Plate", hc_before: np.ndarray, account: str) -> float:
    """After a whole-column thinning that kept node order (stretching, failed rifts), scale
    each node's cratonic crust by its Hc ratio and book the loss to `account`."""
    craton = plate.collect("craton_crust_m")
    if not np.any(craton > 0.0):
        return 0.0
    hc = plate.collect("crustal_thickness_m")
    thinned = scale_with_column(craton, hc_before, hc)
    lost = float(np.dot(craton - thinned, _areas(world, plate)))
    if lost != 0.0:
        _write(plate, thinned)
        record(world, account, lost)
    return lost


def scale_with_column(craton: np.ndarray, hc_before: np.ndarray, hc_after: np.ndarray) -> np.ndarray:
    """Cratonic crust after its column went from `hc_before` to `hc_after`: shrunk in
    proportion where the column thinned, unchanged where it thickened, never above Hc."""
    ratio = np.where((hc_after < hc_before) & (hc_before > 0.0), hc_after / np.where(hc_before > 0.0, hc_before, 1.0), 1.0)
    return np.clip(np.minimum(craton * ratio, hc_after), 0.0, None)


def _write(plate: "Plate", craton: np.ndarray) -> None:
    formed = plate.collect("craton_formed_years")
    plate.set_fields_on_plate(
        craton_crust_m=craton,
        craton_formed_years=np.where(craton > 0.0, formed, CRATON_UNFORMED_YEARS),
    )


class PhaseAudit:
    """Brackets one step phase: on `settle`, holds every plate's craton to its column (booking
    that to `clip_account`) and books whatever cratonic volume vanished without passing
    through an instrumented site -- removed nodes, remap drift -- to `residual_account`.
    A residual gain always goes to `unattributed_m3`."""

    def __init__(self, world: "World") -> None:
        self.world = world
        self.live = live_volume_m3(world)
        self.booked = booked_sinks_m3(world)

    def settle(self, clip_account: str, residual_account: str = "unattributed_m3") -> None:
        world = self.world
        for plate in world.plates:
            clip_to_column(world, plate, clip_account)
        residual = self.live - live_volume_m3(world) - (booked_sinks_m3(world) - self.booked)
        tolerance = 1e-9 * max(self.live, 1.0)
        if residual > tolerance:
            record(world, residual_account, residual)
        elif residual < -tolerance:
            record(world, "unattributed_m3", residual)


# --- Seeding, formation, and the per-step update -----------------------------------------


def seed_initial_cratons(world: "World") -> float:
    """Mark the deep interiors of the world's genuine continents as pre-existing cratons.
    Returns the seeded volume (m^3). No-op on line plates."""
    spacing_rad = line_spacing_rad(world.node_density)
    hops = _hops(CRATON_SEED_MARGIN_KM, spacing_rad)
    seeded = 0.0
    for plate in world.plates:
        host = _genuine_continental(plate)
        interior = _interior_hops(plate, host, hops)
        if interior is None:
            continue
        hc = plate.collect("crustal_thickness_m")
        chosen = host & (hc <= CRATON_HOST_MAX_HC_M) & (interior >= hops) & (plate.collect("craton_crust_m") <= 0.0)
        if not np.any(chosen):
            continue
        craton = plate.collect("craton_crust_m")
        craton[chosen] = hc[chosen]
        formed = plate.collect("craton_formed_years")
        formed[chosen] = world.elapsed_years - CRATON_SEED_AGE_YEARS
        stable = plate.collect("stable_continental_myr")
        stable[chosen] = np.maximum(stable[chosen], CRATON_SEED_AGE_YEARS / 1e6)
        plate.set_fields_on_plate(craton_crust_m=craton, craton_formed_years=formed, stable_continental_myr=stable)
        seeded += float(np.dot(hc[chosen], _areas(world, plate)[chosen]))
    return seeded


def update(world: "World", years: float) -> None:
    """End-of-step craton pass: advance every quiet column's formation clock, cratonise the
    ones that reach CRATON_FORMATION_MYR, and drop the provenance stamp of any cell whose
    craton is gone."""
    ensure_initialized(world)
    years_myr = years / 1e6
    spacing_rad = line_spacing_rad(world.node_density)
    hops = _hops(CRATON_FORMATION_MARGIN_KM, spacing_rad)
    for plate in world.plates:
        host = _genuine_continental(plate)
        interior = _interior_hops(plate, host, hops)
        if interior is None:
            continue
        hc = plate.collect("crustal_thickness_m")
        quiet = host & (hc <= CRATON_HOST_MAX_HC_M) & (interior >= hops)
        stable = np.where(quiet, plate.collect("stable_continental_myr") + years_myr, 0.0)
        craton = plate.collect("craton_crust_m")
        formed = plate.collect("craton_formed_years")
        forming = quiet & (stable >= CRATON_FORMATION_MYR) & (craton <= 0.0)
        if np.any(forming):
            craton[forming] = hc[forming]
            formed[forming] = world.elapsed_years - stable[forming] * 1e6
            record(world, "formed_m3", float(np.dot(hc[forming], _areas(world, plate)[forming])))
        plate.set_fields_on_plate(
            stable_continental_myr=stable,
            craton_crust_m=craton,
            craton_formed_years=np.where(craton > 0.0, formed, CRATON_UNFORMED_YEARS),
        )


def retreat_allowed(world: "World", plate: "Plate") -> np.ndarray:
    """Per node, whether boundary retreat may consume it this step: always for non-cratonic
    cells; for a craton cell, only once it has been overlapped by another plate for
    CRATON_RETREAT_DELAY_YEARS scaled by its strength."""
    craton = plate.collect("craton_crust_m")
    allowed = craton <= 0.0
    if np.all(allowed):
        return allowed
    onset = plate.collect("overlap_onset_years")
    overlapped_for = np.where(onset != 0.0, world.elapsed_years - onset, -np.inf)
    return allowed | (overlapped_for >= CRATON_RETREAT_DELAY_YEARS * strength(craton))


def log_destruction(world: "World", before: dict[str, float], live_before: float) -> None:
    """Log one event per mechanism whose destruction this step exceeded CRATON_EVENT_FRACTION
    of the step's starting craton volume."""
    if live_before <= 0.0:
        return
    for key in CRATON_SINK_ACCOUNTS:
        lost = world.craton_ledger.get(key, 0.0) - before.get(key, 0.0)
        if lost > CRATON_EVENT_FRACTION * live_before:
            mechanism = key.removesuffix("_m3").replace("_", " ")
            world.log_event(f"Craton destruction ({mechanism}): {lost / 1e9:,.0f} km³ ({100.0 * lost / live_before:.1f}% of cratons).")


def diagnostics(world: "World") -> dict[str, float]:
    """Live craton inventory, the persisted ledger, and its balance -- for stats and the API."""
    ensure_ledger(world)
    area = 0.0
    continental_area = 0.0
    volume = 0.0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        areas = _areas(world, plate)
        craton = plate.collect("craton_crust_m")
        area += float(areas[craton > 0.0].sum())
        continental_area += float(areas[_genuine_continental(plate)].sum())
        volume += float(np.dot(craton, areas))
    result = {f"craton_{key}": float(world.craton_ledger.get(key, 0.0)) for key in CRATON_ACCOUNTS}
    result["craton_area_m2"] = area
    result["craton_area_fraction_of_continental"] = area / continental_area if continental_area > 0.0 else 0.0
    result["craton_volume_m3"] = volume
    result["craton_balance_error_m3"] = balance_error_m3(world)
    return result

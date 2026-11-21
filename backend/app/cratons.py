"""Persistent cratons: old, stable continental cores that resist ordinary tectonics (issue #274).

State. Four per-node surface fields carry everything (see surface_fields.py):

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
- ``belt_formed_years`` -- when the column last became a belt (below), or
  surface_fields.BELT_UNMARKED_YEARS if it never has.

Formation rule. A column is quiet when it is genuine continental crust (Hc within
[CRATON_HOST_MIN_HC_M, CRATON_HOST_MAX_HC_M] -- neither a drowned shelf nor an active orogen)
at least CRATON_FORMATION_MARGIN_KM from its plate's edge or any non-genuine cell. Once the
clock reaches CRATON_FORMATION_MYR the column qualifies, and its whole Hc becomes cratonic
when the step's formation budget reaches it (#1001): CRATON_FORMATION_RATE_PER_GYR of the
genuine continental area per Gyr, longest-quiet hosts first. Earth still forms a little
craton but loses more -- ~5-10% of the craton area per Gyr -- so cratons also lose area to
plume heads, which strip the keel under them and leave ordinary continental crust
(CRATON_PLUME_DECRATONIZATION_RADIUS_KM). A present-day planet's cratonic area should fall
slowly, not grow.
A cratons-first world
(RFC 0009, continent_layout.py) starts with its generated cratonic cores; every other world
(and the first load of an older save) seeds the deep interiors -- CRATON_SEED_MARGIN_KM in.
Either way the starting cratons already exist, dated CRATON_SEED_AGE_YEARS before the start.

Belts (RFC 0009). Crust welded in a collision or stretched by a rift stays weaker than the
shields around it and is what later rifts reactivate (Vauchez et al., 1997; Heron et al.,
2016), so it never forms craton however long it stays quiet: `update` leaves belt columns out
of `quiet`. A cratons-first world starts with its generated belts (`seed_belts`, dated
CRATON_BELT_SEED_AGE_YEARS before the start). At runtime every writer goes through
`mark_belt`, which dates the continental, non-cratonic cells it is given: suture-accretion receivers (quad_tectonics._accrete_onto_survivors), crust
an overriding plate takes in (quad_tectonics.settle_received_crust, which carries
crust_transfer's upper-plate placement), stacked merge cells (quad_merge._transfer), a failed
rift's thinned band (apply_failed_rift) and the band along a successful split's cut
(merge_split.maybe_split_plate). A merged or coarsened cell is a belt only if most of its area
was (surface_fields.majority_history). Readers use `belt_mask`, so a belt could later heal
after CRATON_BELT_HEAL_MYR; that is infinite here.

Resistance. Strength rises linearly with cratonic thickness up to CRATON_FULL_STRENGTH_HC_M.
It scales back erosional Hc removal (erosion.py), divergent thinning
(lithosphere_plate.deform_columns) and failed-rift thinning (apply_failed_rift, RFC 0003)
alike (`rift_resistance`), collisional shortening (a craton passes it on to weaker
ground instead of thickening -- shortening.py), how much of a rift's stretched footprint a craton donates
(quad_tectonics._allocate_stretch), and how hard a collision front it jams brakes the plates
(torque.collision_drag_coefficients). Cratonic crust doesn't flow, melt or delaminate
(orogeny.py), and boundary retreat never consumes a craton cell whole (`retreat_allowed`): at
a collision front it jams, so the weaker crust around it takes the shortening, and at a trench
its edge wears back by subduction erosion first (`attrite`, #340).

Destruction. Every loss of cratonic volume is booked into ``World.craton_ledger`` under the
mechanism that caused it (CRATON_SINK_ACCOUNTS); ``balance_error_m3`` checks that sources,
sinks and the live inventory agree. Boundary consumption also books the continental material
it removes into the continental-material ledger (continental_ledger.py). Rifting is also split
by path (CRATON_RIFT_PATH_ACCOUNTS), and every split or failed-rift cut records how much craton
its band crosses (CRATON_RIFT_CUT_COUNTERS) -- diagnostics outside the balance (RFC 0003).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse.csgraph import connected_components

from . import lithosphere
from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from .mantle import PLANET_RADIUS_KM
from .surface_fields import BELT_UNMARKED_YEARS, CRATON_UNFORMED_YEARS

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
# How long a column must stay quiet to qualify to cratonise. Real cratonic roots stabilise
# over a few hundred Myr; long enough here that an ordinary collision belt has to stop
# deforming and relax back into the host Hc window before it joins the craton. Qualifying
# isn't forming: the formation budget below decides how much of the qualifying area does.
CRATON_FORMATION_MYR = 200.0
# Earth-like craton budget (#1001). Cratonisation peaked in the Archean and petered out by
# ~1 Ga, as the mantle cooled below the melting that made thick depleted keels; a keel
# destroyed since then (eastern North China, ~130-110 Ma) regrew thin and fertile, not
# cratonic (Liu et al. 2025, Natl Sci Rev nwaf027; Celli et al. 2020, Nat Commun 11:92).
# Meanwhile cratons lose keel area at ~5-10% per Gyr: Liu et al. (2018, JGR 123) put it at
# ~14% of the original craton area over 2 Gyr. So a present-day planet still forms a little
# craton but loses it faster. Before #1001 every host that qualified formed: 26-50% of the
# starting craton volume per Gyr on the RFC 0012 G2 seeds, against 12-20% destroyed, nearly
# all by thinning that never removed craton area (analysis/issue1001).
#
# Formation: each step earns this fraction of the genuine continental area per Gyr
# (World.craton_formation_multiplier scales it), spent on qualifying hosts. 1% of a ~45%
# cratonic share is ~2% of the craton area per Gyr.
CRATON_FORMATION_RATE_PER_GYR = 0.01
# The most credit that can bank while no host qualifies, in Myr of earnings.
CRATON_FORMATION_CREDIT_MAX_MYR = 100.0
# Destruction: a plume head strips the craton within this radius of its centre, leaving
# ordinary continental crust (southern Kaapvaal, the Angolan Shield, Tanzania: Celli et al.
# 2020). It is a lumped sink: Earth also loses keels to subduction and flat slabs (eastern
# North China, western Wyoming), which have no mechanism here, so plume heads carry the whole
# Earth rate. A stripped column's clock restarts, so it can cratonise again only by waiting
# out CRATON_FORMATION_MYR and then winning the budget behind every longer-quiet host. Heads arrive at ~0.1 per Myr and cratons cover ~20% of the sphere, so ~20 heads per
# Gyr hit craton; calibrated (analysis/issue1001) so they strip ~5-10% of the craton area
# per Gyr. World.craton_decratonization_multiplier scales the stripped area.
CRATON_PLUME_DECRATONIZATION_RADIUS_KM = 350.0
# Generation seeding without generated cores (plates-first, sketch and premade worlds): the
# deep interiors of the starting continents. 875 km (seven default cells) covers ~40-50% of
# genuine continental area on default worlds -- the order of Earth's shields plus platforms.
CRATON_SEED_MARGIN_KM = 875.0
CRATON_SEED_AGE_YEARS = 2.5e9
# Belts (RFC 0009). Generated belts are dated this long before the start: Proterozoic welds,
# younger than the cores they join, as for Laurentia (Hoffman, 1988). Provenance only.
CRATON_BELT_SEED_AGE_YEARS = 1.8e9
# How long after it was marked a belt can form craton again. Infinite: in RFC 0009 welded and
# rifted crust never heals.
CRATON_BELT_HEAL_MYR = math.inf

# --- Resistance --------------------------------------------------------------------------

# Cratonic thickness at which resistance is at full strength.
CRATON_FULL_STRENGTH_HC_M = 25_000.0
# Fraction of erosional Hc removal a full-strength craton column resists.
CRATON_EROSION_RESISTANCE = 0.75
# Fraction of divergent and failed-rift thinning a full-strength craton column resists.
CRATON_RIFT_RESISTANCE = 0.75
# Fraction of its share of a rift's stretched footprint a full-strength craton donor
# withholds -- the rest of that footprint is fresh magmatic crust instead.
CRATON_STRETCH_RESISTANCE = 0.75
# Subduction erosion and sediment subduction recycle ~2.5 km^3/yr of crust globally, ~40 km^3
# of eroded forearc per Myr per km of margin (Scholl & von Huene). Cratonic crust is lost at
# large scale only this way, by attrition at its margins (#340), so a craton edge at a trench
# wears back at this rate instead of being consumed whole. The rate is measured at
# ocean-continent margins and doesn't carry over to collision fronts, where cratons jam. It
# is the nonaccreting-margin value; accreting margins erode more slowly (~12 km^3/Myr/km),
# and the model doesn't yet tell the two apart.
CRATON_ATTRITION_KM2_PER_MYR = 40.0

# --- Ledger ------------------------------------------------------------------------------

CRATON_SOURCE_ACCOUNTS = ("initial_m3", "formed_m3")
# Each declared destruction mechanism. "collision_reworked" is cratonic crust consumed at a
# suture or reworked by faulting: the material stays on the surface as ordinary continental
# crust. "topology_removed" is craton on a plate or fragment removed outright by topology
# cleanup. "audit_clipped" is craton a catch-all phase audit (the deform phase, the end of
# the step) had to clip to its column because some writer thinned the column below it
# without booking it: kept apart from any mechanism so the real sinks can be read (#340).
# "merge_capped" is craton stacked past the suture cap in a plate merge that the overflow
# placement couldn't place (RFC 0001, #441). Merges booked it to "delaminated" before then,
# which read as root foundering (forbidden for cratons, #340); nothing books "delaminated"
# now, and it is kept so older saves' historical merge losses still balance.
# "decratonized" is craton under a plume head that lost its keel (#1001): like
# "collision_reworked", the crust stays on the surface as ordinary continental crust (Hc is
# untouched), but the column loses the strength, rift resistance and jamming the craton label
# carries. That is #340's "it is the mantle keel that goes" in the interior; the model has no
# separate keel, so the label is what goes, and with it the crust's cratonic provenance
# (#1003 separates keel state from crust).
# "unattributed" is signed and should stay ~0 -- anything there escaped every instrumented
# site.
CRATON_SINK_ACCOUNTS = (
    "rifted_m3",
    "subducted_m3",
    "delaminated_m3",
    "merge_capped_m3",
    "collision_reworked_m3",
    "eroded_m3",
    "decratonized_m3",
    "topology_removed_m3",
    "audit_clipped_m3",
    "unattributed_m3",
)
CRATON_ACCOUNTS = CRATON_SOURCE_ACCOUNTS + CRATON_SINK_ACCOUNTS
# Craton area formed and decratonized, and the formation budget earned and discarded by its
# cap (#1001), m^2: diagnostics outside the volume balance, because Earth's craton budget is
# reported in area.
CRATON_AREA_COUNTERS = (
    "formed_area_m2",
    "decratonized_area_m2",
    "formation_credit_earned_m2",
    "formation_credit_discarded_m2",
)
# rifted_m3 by path (RFC 0003 phase 1): divergent-boundary thinning
# (lithosphere_plate.deform_columns), stretch donated to new rift cells
# (quad_tectonics._allocate_stretch) and failed rifts (merge_split.maybe_split_plate). Booking
# one of these also books rifted_m3, so they are diagnostics outside the balance; their sum
# equals rifted_m3 only from when a save started tracking them.
CRATON_RIFT_PATH_ACCOUNTS = ("rifted_divergent_m3", "rifted_stretch_m3", "rifted_failed_m3")
# Every plate-split cut that cleared the split gates, by outcome ("split" broke the plate up,
# "failed_rift" became an aulacogen): how many there were, how many had craton in the cut's
# band (the band a failed rift thins), and the band's total and cratonic area -- the share of
# cuts through cratons that steering the cut (RFC 0003) would act on.
CRATON_RIFT_CUT_KINDS = ("split", "failed_rift")
CRATON_RIFT_CUT_COUNTERS = tuple(
    f"{kind}_{name}" for kind in CRATON_RIFT_CUT_KINDS
    for name in ("cuts", "cuts_through_craton", "band_area_m2", "band_craton_area_m2")
)
# A single step's destruction above this fraction of the live craton volume is logged to the
# event console, naming the mechanism.
CRATON_EVENT_FRACTION = 0.01


def strength(craton_crust_m: np.ndarray) -> np.ndarray:
    """0..1 resistance per node, linear in cratonic thickness."""
    return np.clip(np.asarray(craton_crust_m, dtype=float) / CRATON_FULL_STRENGTH_HC_M, 0.0, 1.0)


def rift_resistance(craton_crust_m: np.ndarray) -> np.ndarray:
    """0..1 share of rift thinning (divergent or failed) a column resists, clamped so a
    tuned CRATON_RIFT_RESISTANCE can't push a thinning factor out of [0, 1] (RFC 0003)."""
    return np.clip(CRATON_RIFT_RESISTANCE * strength(craton_crust_m), 0.0, 1.0)


def heated_strength(
    craton_crust_m: np.ndarray,
    plume_heat_myr: np.ndarray,
    step_years: float,
) -> np.ndarray:
    """Craton strength scaled down by plume heat (RFC 0007, #582).

    Cells carrying plume heat (`plume_heat_myr > 0`) keep only
    ``(1 - PLUME_CRATON_WEAKENING)`` of their strength; a step that extends
    beyond the remaining heat weakens only the heated fraction
    (``remaining / step``), so decay is gradual as the countdown expires. No
    heat anywhere returns exactly `strength`.
    """
    from . import plumes

    weakening = float(np.clip(plumes.PLUME_CRATON_WEAKENING, 0.0, 1.0))
    heat = np.asarray(plume_heat_myr, dtype=float)
    step_myr = max(float(step_years) / 1e6, 0.0)
    if step_myr <= 0.0:
        return strength(craton_crust_m)
    heated_fraction = np.clip(heat / step_myr, 0.0, 1.0)
    return strength(craton_crust_m) * (1.0 - weakening * heated_fraction)


def heated_rift_resistance(
    craton_crust_m: np.ndarray,
    plume_heat_myr: np.ndarray,
    step_years: float,
) -> np.ndarray:
    """Craton-rooted rift resistance weakened by plume heat (RFC 0007, #582).

    While cells carry plume heat (`plume_heat_myr > 0`), craton strength is
    scaled by ``(1 - PLUME_CRATON_WEAKENING)`` and the product is clamped to
    [0, 1] exactly like `rift_resistance`. A step that extends beyond the
    remaining heat weakens only the heated fraction (``remaining / step``),
    so decay is gradual as the countdown expires rather than an abrupt return
    to full strength. No heat anywhere returns exactly `rift_resistance`.
    Tails convention holds: tails stamp no tracer, so heat (and weakening)
    is a head-only state.
    """
    return np.clip(
        CRATON_RIFT_RESISTANCE * heated_strength(craton_crust_m, plume_heat_myr, step_years),
        0.0,
        1.0,
    )


def _hops(km: float, spacing_rad: float) -> int:
    return max(1, round(km / (spacing_rad * PLANET_RADIUS_KM)))


def _continental(plate: "Plate") -> np.ndarray:
    return effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")


def _genuine_continental(plate: "Plate") -> np.ndarray:
    return _continental(plate) & (plate.collect("crustal_thickness_m") >= CRATON_HOST_MIN_HC_M)


def _interior_hops(plate: "Plate", host: np.ndarray, max_hops: int) -> np.ndarray | None:
    """Per cell, edge hops to the nearest plate edge or non-`host` cell, saturating at
    `max_hops + 1`. None on a plate with no cells."""
    if plate.node_count() == 0:
        return None
    from .quad_tectonics import hop_distance

    edge = np.any(np.any(plate._probe_neighbour_indices() < 0, axis=2), axis=1)
    return hop_distance(plate, edge | ~host, max_hops)


# --- Ledger bookkeeping ------------------------------------------------------------------


def ensure_ledger(world: "World") -> None:
    """Backfill the ledger's accounts without touching any plate -- safe on load."""
    if not isinstance(getattr(world, "craton_ledger", None), dict):
        world.craton_ledger = {}
    for key in CRATON_ACCOUNTS + CRATON_RIFT_PATH_ACCOUNTS + CRATON_RIFT_CUT_COUNTERS + CRATON_AREA_COUNTERS:
        if key != "initial_m3":
            world.craton_ledger.setdefault(key, 0.0)


def ensure_initialized(world: "World") -> None:
    """`ensure_ledger`, then -- the first time a world with plates is generated or stepped --
    seed its starting cratons. The `initial_m3` key's presence marks a seeded world, so a
    world with no eligible interior is never reseeded. Loading a save never seeds: an older
    save (including one converted from line plates) gets its cratons on its first step."""
    ensure_ledger(world)
    if "initial_m3" not in world.craton_ledger and world.plates:
        world.craton_ledger["initial_m3"] = seed_initial_cratons(world)


def record(world: "World", account: str, volume_m3: float) -> None:
    """Book `volume_m3` of cratonic crust into one account. Only `unattributed_m3` may go
    negative (a remap that gained volume). A rift path account books rifted_m3 too."""
    if account not in CRATON_ACCOUNTS and account not in CRATON_RIFT_PATH_ACCOUNTS:
        raise KeyError(f"unknown craton ledger account: {account}")
    volume_m3 = float(volume_m3)
    if not np.isfinite(volume_m3) or (volume_m3 < 0.0 and account != "unattributed_m3"):
        raise ValueError(f"craton ledger volume must be finite and non-negative, got {volume_m3!r}")
    if volume_m3 == 0.0:
        return
    ensure_ledger(world)
    world.craton_ledger[account] = world.craton_ledger.get(account, 0.0) + volume_m3
    if account in CRATON_RIFT_PATH_ACCOUNTS:
        world.craton_ledger["rifted_m3"] += volume_m3


def record_rift_cut(world: "World", plate: "Plate", points: np.ndarray, cut_normal: np.ndarray, band_sin: float, kind: str) -> None:
    """Count one split cut of `kind` (CRATON_RIFT_CUT_KINDS) through `plate` along the great
    circle with unit normal `cut_normal`, and the area and craton area of the nodes (world
    positions `points`, in collect order) within `band_sin` of it."""
    in_band = np.abs(points @ np.asarray(cut_normal, dtype=float)) < band_sin
    record_cut_areas(world, kind, *band_cut_areas(world, plate, in_band))


def band_cut_areas(world: "World", plate: "Plate", in_band: np.ndarray) -> tuple[float, float]:
    """(area, craton area), in m2, of the nodes of `plate` in the boolean mask `in_band` (collect order). Read it
    before the plate is cut."""
    areas = _areas(world, plate)[in_band]
    return float(areas.sum()), float(areas[plate.collect("craton_crust_m")[in_band] > 0.0].sum())


def record_cut_areas(world: "World", kind: str, band_area: float, craton_area: float) -> None:
    """Count one cut of `kind` (CRATON_RIFT_CUT_KINDS) with its band's area and craton area. Both split methods book
    their cuts here, so the account names live in one place (the rift-line method's band follows its line, PR #950)."""
    if kind not in CRATON_RIFT_CUT_KINDS:
        raise KeyError(f"unknown rift cut kind: {kind}")
    ensure_ledger(world)
    ledger = world.craton_ledger
    ledger[f"{kind}_cuts"] += 1.0
    ledger[f"{kind}_cuts_through_craton"] += float(craton_area > 0.0)
    ledger[f"{kind}_band_area_m2"] += band_area
    ledger[f"{kind}_band_craton_area_m2"] += craton_area


def _areas(world: "World", plate: "Plate") -> np.ndarray:
    return plate.accounting_areas_m2(line_spacing_rad(world.node_density))


def plate_volume_m3(world: "World", plate: "Plate") -> float:
    if plate.node_count() == 0:
        return 0.0
    return float(np.dot(plate.collect("craton_crust_m"), _areas(world, plate)))


def live_volume_m3(world: "World") -> float:
    return sum(plate_volume_m3(world, plate) for plate in world.plates)


def genuine_continental_area_m2(world: "World") -> float:
    """Area of continental crust thick enough to host a craton (CRATON_HOST_MIN_HC_M)."""
    return sum(float(_areas(world, plate)[_genuine_continental(plate)].sum()) for plate in world.plates if plate.node_count())


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


def seed_cores(plate: "Plate", core: np.ndarray, elapsed_years: float = 0.0) -> None:
    """Make the whole Hc of each `core` column (a mask in node order) cratonic, dated
    CRATON_SEED_AGE_YEARS before `elapsed_years` -- generation's cratonic cores on a
    cratons-first world (RFC 0009, continent_layout.py)."""
    if not np.any(core):
        return
    hc = plate.collect("crustal_thickness_m")
    craton = plate.collect("craton_crust_m")
    craton[core] = hc[core]
    formed = plate.collect("craton_formed_years")
    formed[core] = elapsed_years - CRATON_SEED_AGE_YEARS
    stable = plate.collect("stable_continental_myr")
    stable[core] = np.maximum(stable[core], CRATON_SEED_AGE_YEARS / 1e6)
    plate.set_fields_on_plate(craton_crust_m=craton, craton_formed_years=formed, stable_continental_myr=stable)


def seed_belts(plate: "Plate", belt: np.ndarray, elapsed_years: float = 0.0) -> None:
    """Mark generation's belts (a mask in node order) on a cratons-first world (RFC 0009,
    continent_layout.py), dated CRATON_BELT_SEED_AGE_YEARS before `elapsed_years`."""
    if not np.any(belt):
        return
    formed = plate.collect("belt_formed_years")
    formed[belt] = elapsed_years - CRATON_BELT_SEED_AGE_YEARS
    plate.set_fields_on_plate(belt_formed_years=formed)


def seed_initial_cratons(world: "World") -> float:
    """The world's starting cratons; returns their volume (m^3). A cratons-first world's
    generation has already seeded its cores (`seed_cores`), so they are only booked. Every
    other world -- a sketch, a premade world, a plates-first random world or an old save's
    first step -- gets the deep interiors of its genuine continents, CRATON_SEED_MARGIN_KM in."""
    generated = live_volume_m3(world)
    if generated > 0.0:
        return generated
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


def mark_belt(world: "World", plate: "Plate", mask: np.ndarray) -> None:
    """Date the cells in `mask` (node order) as belt at `world.elapsed_years`. Only continental
    cells holding no craton are marked: belts are continental crust between the shields, and
    crust welded onto a shield is ordinary crust over it."""
    mask = np.asarray(mask, dtype=bool)
    if plate.node_count() == 0 or not np.any(mask):
        return
    mask = mask & _continental(plate) & (plate.collect("craton_crust_m") <= 0.0)
    if not np.any(mask):
        return
    formed = plate.collect("belt_formed_years")
    formed[mask] = world.elapsed_years
    plate.set_fields_on_plate(belt_formed_years=formed)


def keep_belts_off_cratons(fields: dict[str, np.ndarray]) -> None:
    """In place on a remap's output `fields`: a cell holding craton is not a belt, whatever
    the majority vote gave it. Belts and cores stay disjoint, as `mark_belt` keeps them, so a
    mark can't linger under a craton and block re-formation once that craton is lost."""
    if "belt_formed_years" in fields and "craton_crust_m" in fields:
        fields["belt_formed_years"][fields["craton_crust_m"] > 0.0] = BELT_UNMARKED_YEARS


def belt_mask(plate: "Plate", elapsed_years: float) -> np.ndarray:
    """Per node, whether the column is a belt: marked, and marked less than
    CRATON_BELT_HEAL_MYR ago."""
    formed = plate.collect("belt_formed_years")
    return (formed != BELT_UNMARKED_YEARS) & (elapsed_years - formed < CRATON_BELT_HEAL_MYR * 1e6)


def clear_oceanic_belts(plate: "Plate") -> None:
    """Unmark every belt cell whose crust is no longer continental. Writers retype crust in
    many places; the end-of-step `update` clears the marks they leave on oceanic cells, so a
    cell that becomes continental again is a belt only if a writer marks it anew."""
    if plate.node_count() == 0:
        return
    formed = plate.collect("belt_formed_years")
    stale = (formed != BELT_UNMARKED_YEARS) & ~_continental(plate)
    if np.any(stale):
        formed[stale] = BELT_UNMARKED_YEARS
        plate.set_fields_on_plate(belt_formed_years=formed)


def update(world: "World", years: float) -> None:
    """End-of-step craton pass (#1001): strip the craton under every plume head that arrived this
    step (`_decratonize_under_plume_heads`), advance every quiet column's formation clock, spend
    the step's formation budget on the qualifying hosts (`_formation_picks`), and drop the
    provenance stamp of any cell whose craton is gone. A belt column is never quiet (RFC
    0009), and a cell retyped oceanic stops being a belt: belts are continental crust
    (`clear_oceanic_belts`)."""
    ensure_initialized(world)
    # Strip first, so the step's formation budget isn't spent on a column a head strips now.
    _decratonize_under_plume_heads(world, years)
    years_myr = years / 1e6
    spacing_rad = line_spacing_rad(world.node_density)
    hops = _hops(CRATON_FORMATION_MARGIN_KM, spacing_rad)
    states = []
    continental_area = 0.0
    largest_host_cell = 0.0
    for plate in world.plates:
        clear_oceanic_belts(plate)
        host = _genuine_continental(plate)
        interior = _interior_hops(plate, host, hops)
        if interior is None:
            continue
        areas = _areas(world, plate)
        continental_area += float(areas[host].sum())
        largest_host_cell = max(largest_host_cell, float(areas[host].max(initial=0.0)))
        hc = plate.collect("crustal_thickness_m")
        quiet = host & (hc <= CRATON_HOST_MAX_HC_M) & (interior >= hops) & ~belt_mask(plate, world.elapsed_years)
        stable = np.where(quiet, plate.collect("stable_continental_myr") + years_myr, 0.0)
        craton = plate.collect("craton_crust_m")
        forming = quiet & (stable >= CRATON_FORMATION_MYR) & (craton <= 0.0)
        states.append((plate, areas, hc, stable, craton, forming))
    picks = _formation_picks(world, states, continental_area, years, largest_host_cell)
    for (plate, areas, hc, stable, craton, _), take in zip(states, picks):
        formed = plate.collect("craton_formed_years")
        if np.any(take):
            craton[take] = hc[take]
            formed[take] = world.elapsed_years - stable[take] * 1e6
            record(world, "formed_m3", float(np.dot(hc[take], areas[take])))
            world.craton_ledger["formed_area_m2"] += float(areas[take].sum())
        plate.set_fields_on_plate(
            stable_continental_myr=stable,
            craton_crust_m=craton,
            craton_formed_years=np.where(craton > 0.0, formed, CRATON_UNFORMED_YEARS),
        )


def _formation_picks(
    world: "World", states: list, continental_area_m2: float, years: float, largest_host_cell_m2: float
) -> list[np.ndarray]:
    """Per plate in `states`, the qualifying hosts that cratonise this step (#1001).

    The step earns CRATON_FORMATION_RATE_PER_GYR x `craton_formation_multiplier` of the genuine
    continental area, added to `World.craton_formation_credit_m2`. The credit is spent on whole
    cells: the longest-quiet first, and among equally quiet cells the nearest to existing
    craton, so new craton grows onto old blocks' quiet margins and as coherent blocks rather
    than speckle; remaining ties go by a seeded hash of plate and cell, so no plate's place in
    `world.plates` decides where craton forms (PR #1002 review). A cell bigger than the credit
    waits. Only then is the unspent credit capped at CRATON_FORMATION_CREDIT_MAX_MYR of
    earnings (so a long wait for hosts can't release one giant burst), floored at the largest
    host cell so a coarse grid or low multiplier can't stall formation, and the excess booked to the `formation_credit_discarded_m2` counter. Capping
    after spending keeps the budget independent of how an interval is split into steps, while
    hosts are available. A multiplier of 0 discards the credit. So, per world,
    `formation_credit_earned_m2` = `formed_area_m2` + credit + `formation_credit_discarded_m2`."""
    picks = [np.zeros(len(state[5]), dtype=bool) for state in states]
    multiplier = float(world.craton_formation_multiplier)
    rate = CRATON_FORMATION_RATE_PER_GYR * multiplier if multiplier > 0.0 and years > 0.0 else 0.0  # NaN-safe
    earned = rate * continental_area_m2 * years / 1e9
    world.craton_ledger["formation_credit_earned_m2"] += earned
    credit = float(world.craton_formation_credit_m2) + earned
    candidates = [(i, np.flatnonzero(state[5])) for i, state in enumerate(states) if np.any(state[5])]
    if credit > 0.0 and candidates:
        credit -= _spend_formation_credit(world, states, candidates, credit, picks)
    cap = max(rate * continental_area_m2 * CRATON_FORMATION_CREDIT_MAX_MYR / 1e3, largest_host_cell_m2) if rate > 0.0 else 0.0
    _bank_formation_credit(world, credit, cap)
    return picks


def _spend_formation_credit(world: "World", states: list, candidates: list, credit: float, picks: list[np.ndarray]) -> float:
    """Mark in `picks` the candidate cells `credit` buys (see `_formation_picks`); returns the
    area spent."""
    from .quad_tectonics import hop_distance

    reach = _hops(CRATON_FORMATION_MARGIN_KM, line_spacing_rad(world.node_density))
    owner, index, stable, near, area, tie = [], [], [], [], [], []
    for i, nodes in candidates:
        plate, areas, _, quiet_myr, craton, _ = states[i]
        owner.append(np.full(len(nodes), i))
        index.append(nodes)
        stable.append(quiet_myr[nodes])
        near.append(hop_distance(plate, craton > 0.0, reach)[nodes])
        area.append(areas[nodes])
        tie.append(_tie_break(world.seed, plate.plate_id, nodes))
    owner, index, stable, near, area, tie = (np.concatenate(x) for x in (owner, index, stable, near, area, tie))
    # lexsort: last key first -- longest quiet, then nearest craton, then the seeded hash.
    order = np.lexsort((index, owner, tie, near, -stable))
    spent = np.cumsum(area[order])
    count = int(np.searchsorted(spent, credit, side="right"))
    for k in order[:count]:
        picks[owner[k]][index[k]] = True
    return float(spent[count - 1]) if count else 0.0


def _tie_break(seed: int, plate_id: int, nodes: np.ndarray) -> np.ndarray:
    """A deterministic pseudo-random uint64 per (seed, plate, cell): splitmix64 of the three."""
    mask = (1 << 64) - 1
    base = ((int(seed) & mask) * 0x9E3779B97F4A7C15 + (int(plate_id) & mask) * 0xBF58476D1CE4E5B9) & mask
    with np.errstate(over="ignore"):
        x = np.uint64(base) + np.asarray(nodes, dtype=np.uint64) * np.uint64(0x94D049BB133111EB)
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return x ^ (x >> np.uint64(31))


def _bank_formation_credit(world: "World", credit: float, cap: float) -> None:
    """Carry `credit` to the next step, held to [0, cap]; book what the cap drops."""
    kept = min(max(credit, 0.0), max(cap, 0.0))
    world.craton_ledger["formation_credit_discarded_m2"] += max(credit - kept, 0.0)
    world.craton_formation_credit_m2 = kept


def _decratonize_under_plume_heads(world: "World", years: float) -> None:
    """Strip the craton within CRATON_PLUME_DECRATONIZATION_RADIUS_KM (area-scaled by
    `craton_decratonization_multiplier`) of every plume head that arrived this step (#1001).

    The keel goes, not the crust (design principle 7): the column keeps its Hc and becomes
    ordinary continental crust, so the cratonic volume is booked to `decratonized_m3` and its
    area to the `decratonized_area_m2` counter, and the column's formation clock restarts,
    as a destroyed keel regrows thin and fertile (eastern North China). Heads arrive in
    `plumes.apply` earlier in the step, with arrival times in (t - years, t]. Known gap (PR
    #1002 review): a head is found here only while it is still in `world.mantle_plumes`, so a
    step longer than a plume's whole life (~80 Myr), or the first step after the plume stream
    is unpaused, can miss heads; Myr-scale steps are unaffected."""
    multiplier = float(world.craton_decratonization_multiplier)
    if not (years > 0.0 and multiplier > 0.0):  # NaN-safe
        return
    t1 = world.elapsed_years
    t0 = t1 - years
    heads = [plume for plume in world.mantle_plumes if not plume.is_tail_only and t0 < plume.born_years <= t1]
    if not heads:
        return
    centres = np.array([plume.position for plume in heads], dtype=float)
    centres /= np.linalg.norm(centres, axis=1, keepdims=True)
    # A spherical cap's area is 2 pi R^2 (1 - cos theta), so scaling the area by the multiplier
    # scales 1 - cos; clamped at the whole sphere, so a huge multiplier can't wrap past pi.
    cos_reach = max(-1.0, 1.0 - multiplier * (1.0 - math.cos(CRATON_PLUME_DECRATONIZATION_RADIUS_KM / PLANET_RADIUS_KM)))
    stripped_area = 0.0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        craton = plate.collect("craton_crust_m")
        if not np.any(craton > 0.0):
            continue
        points = plate.all_points_and_elevation()[0]
        unit = points / np.linalg.norm(points, axis=1, keepdims=True)
        hit = (craton > 0.0) & np.any(unit @ centres.T >= cos_reach, axis=1)
        if not np.any(hit):
            continue
        areas = _areas(world, plate)
        record(world, "decratonized_m3", float(np.dot(craton[hit], areas[hit])))
        stripped_area += float(areas[hit].sum())
        craton[hit] = 0.0
        stable = plate.collect("stable_continental_myr")
        stable[hit] = 0.0
        formed = plate.collect("craton_formed_years")
        formed[hit] = CRATON_UNFORMED_YEARS
        plate.set_fields_on_plate(craton_crust_m=craton, stable_continental_myr=stable, craton_formed_years=formed)
    if stripped_area > 0.0:
        world.craton_ledger["decratonized_area_m2"] += stripped_area
        world.log_event(f"Plume head decratonized {stripped_area / 1e12:.2f} M km2 of craton (keel lost, crust kept).")


def retreat_allowed(plate: "Plate") -> np.ndarray:
    """Per node, whether boundary retreat may consume it: only cells with no craton left.
    Cratons are lost at large scale only by attrition at their margins (#340): at a collision
    front a craton cell jams, and the polarized front falls back to the other side's retreat
    (lithosphere_plate.boundary_context), so the weak crust between the blocks takes the
    shortening; craton against craton holds until the plates merge. At a trench `attrite`
    wears the cell's craton away first. Over 300 Myr on seeds 7, 11 and 2024
    (bin/debug/compare_cratons.py), this kept 87-105% of the cratons against 59-74% for the
    old 60 Myr consumption delay, and 10% more genuine continental area on two seeds (4% less
    on the third). #274 had found never consuming kept no more continent, but flow and melt
    then drained cratons anyway (analysis/issue340)."""
    return plate.collect("craton_crust_m") <= 0.0


def attrite(world: "World", plate: "Plate", margin: np.ndarray, years: float, max_m3: float | None = None) -> float:
    """Subduction erosion of the craton cells in `margin` (a trench's exposed edge) over
    `years`: each loses CRATON_ATTRITION_KM2_PER_MYR times its edge length (the square root
    of its area) of crust from its craton, which goes down the trench -- craton ledger
    `subducted_m3`, its continental material `deeply_subducted_m3`. Hc thins with it and the
    surface subsides isostatically. `max_m3` caps the total (a trench's #177 budget), scaling
    every cell's loss alike. Returns the volume lost (m^3)."""
    craton = plate.collect("craton_crust_m")
    margin = np.asarray(margin, dtype=bool) & (craton > 0.0)
    if years <= 0.0 or not np.any(margin):
        return 0.0
    from . import continental_ledger, rheology

    areas = _areas(world, plate)
    hc = plate.collect("crustal_thickness_m")
    material = plate.collect("continental_material_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    codes = plate.collect("crust_type_code")
    # km^2/Myr is m^2/yr: the volume per unit margin length, spread over the cell.
    loss = np.where(margin, np.minimum(CRATON_ATTRITION_KM2_PER_MYR * years / np.sqrt(areas), craton), 0.0)
    wanted = float(loss @ areas)
    if max_m3 is not None and wanted > max_m3:
        loss = loss * (max(max_m3, 0.0) / wanted)
    material_loss = np.minimum(loss, material)
    new_hc = hc - loss
    density = lithosphere.node_crust_density(codes, plate.crust_type)
    shift = lithosphere.isostatic_elevation(new_hc, hm, density) - lithosphere.isostatic_elevation(hc, hm, density)
    lost = float(loss @ areas)
    record(world, "subducted_m3", lost)
    continental_ledger.record(world, "deeply_subducted_m3", float(material_loss @ areas))
    plate.set_fields_on_plate(
        crustal_thickness_m=new_hc,
        continental_material_m=material - material_loss,
        restite_m=np.minimum(plate.collect("restite_m"), new_hc),
        elevation=rheology.clip_elevation_bounds(plate.collect("elevation") + shift),
    )
    _write(plate, craton - loss)
    return lost


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


def core_count(plate: "Plate") -> int:
    """Cratonic cores on `plate`: its edge-connected components of craton cells."""
    core = np.flatnonzero(plate.collect("craton_crust_m") > 0.0)
    if not len(core):
        return 0
    from .quad_tectonics import _adjacency_matrix

    return int(connected_components(_adjacency_matrix(plate)[core][:, core], directed=False)[0])


def diagnostics(world: "World") -> dict[str, float]:
    """Live craton inventory, the persisted ledger, and its balance -- for stats and the API."""
    ensure_ledger(world)
    area = 0.0
    continental_area = 0.0
    volume = 0.0
    belt_area = 0.0
    cores = 0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        areas = _areas(world, plate)
        craton = plate.collect("craton_crust_m")
        area += float(areas[craton > 0.0].sum())
        continental_area += float(areas[_genuine_continental(plate)].sum())
        volume += float(np.dot(craton, areas))
        # Only belts on continental crust: a marked cell retyped oceanic mid-step is no belt.
        belt_area += float(areas[belt_mask(plate, world.elapsed_years) & _continental(plate)].sum())
        cores += core_count(plate)
    keys = CRATON_ACCOUNTS + CRATON_RIFT_PATH_ACCOUNTS + CRATON_RIFT_CUT_COUNTERS + CRATON_AREA_COUNTERS
    result = {f"craton_{key}": float(world.craton_ledger.get(key, 0.0)) for key in keys}
    result["craton_area_m2"] = area
    result["craton_area_fraction_of_continental"] = area / continental_area if continental_area > 0.0 else 0.0
    result["craton_volume_m3"] = volume
    result["craton_belt_area_m2"] = belt_area
    result["craton_core_count"] = float(cores)
    result["craton_balance_error_m3"] = balance_error_m3(world)
    return result

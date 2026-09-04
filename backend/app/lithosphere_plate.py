"""`LithospherePlate`: the plate representation. Subclasses `PlateWithLines` (see
plates.py) rather than reinventing its plumbing -- everything not about *how a plate moves or
deforms* (outline tracing, containment, neighbour search, node iteration, the whole render/
erosion/hydrology/stats/persistence surface) is inherited unchanged, since every one of those
consumers only ever reads `Plate`'s abstract interface (`all_points_and_elevation`, `collect`,
`contains_batch`, ...), never `shift`/`deform`'s own internals. See the plan's "Why subclass
PlateWithLines" section.

Overridden here: `shift` (torque integration, torque.py), `deform` (Mohr-Coulomb/isostasy,
rheology.py + lithosphere.py), `merge_with`/`_merge_nodes_with`/`split`/`grow_into` (carrying
Hc/Hm through the representation's own resample/partition operations, which
`PlateWithLines`'s versions don't know about).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from . import geometry
from .elevation_lines import (
    ELEV_CHANGE_COLLISION,
    ELEV_CHANGE_MIN_DELTA_M,
    ELEV_CHANGE_NEW_CRUST,
    ELEV_CHANGE_RIFT,
    ELEV_CHANGE_SUBDUCTION_ARC,
    ELEV_CHANGE_TRANSFORM,
    ELEV_CHANGE_TRENCH,
    ELEV_CHANGE_VOLCANO,
    ElevationLine,
    build_lines_from_lattice,
    line_spacing_rad,
    needs_regularizing,
    regularize_line,
    split_into_contiguous_runs,
)
from .noise import SphereNoise
from .plates import (
    CONTINENTAL_FRACTION,
    MIN_AUTO_PLATES,
    MAX_AUTO_PLATES,
    MIN_OCEANIC_PLATES,
    POLE_CAP_MARGIN_MULT,
    PlateWithLines,
    _INTERIOR_SUBDUCTION_MIN_RUN,
    _ROW_FULL_REVOLUTION_SLACK,
    _contested_by_any,
    _land_noise_threshold,
    _row_median_step,
)
from . import bathymetry, lithosphere, rheology, terrain_noise, torque

EXTEND_THRESHOLD_MULTIPLIER = 1.3  # same shape as v1's plates.EXTEND_THRESHOLD_RAD
MAX_EXTEND_NODES_PER_STEP = 400

# Transform (strike-slip) boundary pressure-ridge uplift, applied as a direct elevation delta
# on the transform band in deform() (there is no net crustal shortening at a strike-slip
# contact, so it does not go through Hc/Hm). Gentler than either convergent case -- real
# transform relief is local pressure ridges / transtensional sags, not an orogen. Half v1's
# plates.TRANSFORM_UPLIFT_RATE_M_PER_MYR (200) since here it is not distance-tapered, only
# fault_influence-gated.
TRANSFORM_UPLIFT_RATE_M_PER_MYR = 100.0

# A continental line's *contested* end is allowed to retreat -- one node per step -- whether
# the overriding neighbour is oceanic (a passive margin / accretion front: the ocean slab
# descends under it and the buried continental node cedes nothing the model should keep) or
# continental (a suture whose overlapping crust is consumed into the orogen -- its volume is
# not discarded but thrust back onto the plate's own surviving leading edge, see
# `_redistribute_accreted_column` / SUTURE_ACCRETION_SPREAD_NODES, so the belt builds real
# relief in proportion to the overlap it actually eats). Left un-retreatable, a contested end
# still grows at its *other* (divergent) side every step and never back -- the continental
# node ratchet that drives the unbounded node-count creep and the slow land-fraction decline
# (docs/TODO.md "Node-count creep") -- and, for a continent-continent pile-up, a deep
# territory overlap that just sat there for tens of Myr until the forced-merge timer fused
# the pair (the `overlapAge` view's stalled multi-plate collisions). Retreat is gated:
#   - one node per step (the existing `n_distance_cap` / `max_extend_nodes` caps already do
#     this at continental drift rates), and
#   - only where the contested node is part of a run of at least this many consecutive
#     contested nodes, so a single stray node from bounding-polygon envelope fuzz can't
#     nibble a stable coastline or, worse, sever a lobe into a spurious defragmentation plate
#     (the failure the naive "retreat every continental contested node" experiment hit -- see
#     that TODO section; the interior-subduction carve below also stays oceanic-only for the
#     same lobe-severing reason).
CONTINENTAL_CONTESTED_RETREAT_MIN_RUN = 3

# Whole-row retreat -- the reverse of `_claim_adjacent_territory`, and the *only* retreat op
# available in the "parallel suture" regime: when a neighbour overrides a continental plate's
# frontmost phi-row over its full theta width there is no uncontested end for
# `_grow_or_shrink_line_for_deform` to trim, and a continental row is never carved mid-span
# (that severs the landmass into a spurious defragmentation plate), so end-trim alone leaves
# that plate physically unable to give ground -- its trailing edge still grows every step, so
# the node pile ratchets outward regardless (docs/TODO.md#continental-ratchet-solution,
# mechanism 2). Once a plate's outermost row (either phi extreme) has been at least
# LEADING_ROW_CONTESTED_FRACTION contested for a cumulative LEADING_ROW_RETREAT_SUSTAINED_YEARS
# of deform time, the whole row is dropped. Whole-row removal keeps the plate contiguous -- the
# lobe-severing hazard is specific to *mid*-row carving -- so this is safe exactly where the
# interior-subduction carve is not. Like the 2026-09-02 end-retreat this does not plumb the
# dropped column's volume anywhere; the newly-exposed frontmost row is contested next step and
# thickens through the ordinary `CONTINENTAL_COLLISION_SHORTENING_BOOST` path.
LEADING_ROW_CONTESTED_FRACTION = 0.7
LEADING_ROW_RETREAT_SUSTAINED_YEARS = 5_000_000.0
# Never drop a row that would take the plate below this many rows -- a tiny plate has no
# "leading row" worth the name and the contiguity argument gets thin.
LEADING_ROW_DROP_MIN_ROWS = 4
# Volume-budget growth gate (docs/TODO.md "Continental ratchet: solution design",
# mechanism 1). A lattice node's physical footprint is constant across the sphere by
# construction (`lithosphere.node_area_m2`), so a plate's total area is just its node count
# times that -- and its implied *mean* crustal thickness is `mean(crustal_thickness_m)`.
# The continental boundary ratchet dilutes this: `_grow_or_shrink_line_for_deform` and
# `_claim_adjacent_territory` seed every new margin node at the *oceanic* reference column
# (`growth_seed_thickness`), and nothing ever removes a whole leading row, so a
# shear-stretched continental plate tiles unbounded drowned passive-margin outward -- node
# count creeps ~+5-6% per 150 My and the plate interior isostatically oceanises into a
# "giant 80%-drowned continent" (docs/TODO.md items 2 / 5, and the land-fraction decline).
#
# The gate counts a plate's *genuine* continental nodes -- Hc at least
# `CONTINENTAL_BUDGET_HC_FRACTION` of the continental reference -- and, once the plate's
# total node count exceeds `CONTINENTAL_AREA_BUDGET_MULT` times that count, suppresses all
# areal *growth* for the step (end-growth here and whole new rows in
# `_claim_adjacent_territory`). Retreat, divergent thinning and convergent thickening keep
# running, so an over-budget plate thins / drowns / crumples back toward its crustal volume
# rather than merely freezing. A real craton sits near reference Hc across its whole area,
# nowhere near the cap; the >1 multiplier is the realistic shelf + accreted-terrane
# allowance. Regime-independent and neighbour-independent -- unlike the contested-run
# retreat it does not care how the suture sits against the row grid.
CONTINENTAL_BUDGET_HC_FRACTION = 0.6
CONTINENTAL_AREA_BUDGET_MULT = 1.8

# When a continental *suture* end retreats (a continental neighbour overrides it -- not an
# oceanic one, where the buried column genuinely subducts and is lost), the removed column's
# crustal volume is conserved: it is thrust back onto the plate's own surviving leading-edge
# nodes, spread over this many of them (an imbricate thrust wedge), the attached mantle
# lithosphere thickening in proportion. This is the mass-honest replacement for the retired
# `rheology.CONTINENTAL_COLLISION_SHORTENING_BOOST` fudge (a flat 2.5x `fault_factor`
# multiplier at continent-continent contested nodes, unrelated to how much overlap was
# actually consumed). Node area is constant per node (`lithosphere.node_area_m2`), so
# conserving volume is just moving the summed Hc of the dropped nodes onto the survivors;
# `regularize_line` re-evens the spacing next pass and isostasy lifts the thickened belt.
SUTURE_ACCRETION_SPREAD_NODES = 3

# Hard ceiling on a node's Hc after suture accretion. A suture that never heals (the
# neighbour keeps overriding) would otherwise pile every consumed column onto the same few
# retreating-edge nodes indefinitely -- Hc ran to ~190 km and climbing on a 30-My test run.
# Real orogenic crust does not stack past ~2x reference: the excess root is removed by
# lower-crustal / mantle-lithosphere delamination (and the surface by erosion). Accreted
# mass over this ceiling is dropped (delaminated), so accretion is mass-conserving only up
# to the cap -- which a normal collision, healing over ~1-2 My, never reaches.
SUTURE_ACCRETION_MAX_HC_M = 2.4 * lithosphere.REFERENCE_HC_CONTINENTAL_M

# Active-margin (Cordilleran) accretion. When a continental plate's *leading* edge grows
# into space a subducting oceanic neighbour is vacating (slab rollback / trench retreat),
# the new ground is juvenile arc + accreted-terrane crust, not abyssal sea floor -- so it is
# seeded at this intermediate column (Hc ~0.8x continental reference) rather than
# `growth_seed_thickness`'s drowned oceanic one. This is the deliberately *restricted*
# reverse of the land-area runaway that `growth_seed_thickness` documents: the runaway was
# seeding +200 m dry land on *every* growth event, including growth into open ocean far from
# any margin; seeding a thicker column *only* where the growing end abuts a genuinely
# converging oceanic slab -- and still under the `CONTINENTAL_AREA_BUDGET_MULT` volume gate
# -- is arc accretion, the dominant land-loss driver's actual physical counterweight (see
# docs/TODO.md "Land fraction slowly declines"). The seed lands as shallow forearc/shelf
# (~ -450 m) and builds to land as convergence continues via
# `rheology.apply_arc_magmatic_thickening` + ordinary convergent shortening.
ARC_MARGIN_SEED_HC_M = 28_000.0
ARC_MARGIN_SEED_HM_M = 55_000.0

# How many nodes in from a line end are scanned for an active-margin signal -- a node
# contested by an oceanic neighbour, or one still carrying a subduction-arc provenance stamp
# from a recent step -- when deciding whether that end's growth seeds arc crust or ocean
# floor. Small: the signal only has to survive the one step between the ocean's edge
# retreating and this plate's edge growing into the gap.
ARC_MARGIN_END_SCAN_NODES = 4


def growth_seed_thickness() -> tuple[float, float]:
    """(Hc, Hm) a plate seeds *brand-new areal* nodes with -- when a line grows an end into
    open water (`_grow_or_shrink_line_for_deform`) or claims a whole new phi row
    (`_claim_adjacent_territory`).

    Always the *oceanic* reference column, regardless of the growing plate's own
    `crust_type`: any gap that opens on the sphere is floored by sea-floor spreading, not by
    the neighbouring plate's crust. Seeding a continental plate's own reference column here
    (Hc 35 km / Hm 100 km -> isostatic_elevation = +200 m) was a real land-area runaway: a
    continental plate continuously grows into the space a subducting oceanic plate vacates,
    and continental crust never subducts back, so every such step converted ocean floor into
    +200 m dry land permanently -- measured land fraction climbed 0.27 -> 0.48 and mean
    planet elevation rose ~1.7 km over 180 Myr on seed 559394024. New oceanic crust on a
    continental plate lands ~-3.5 km (a drowned passive margin / accreted terrane); genuine
    continental rifting is untouched, since that thins *existing* crust
    (`rheology.apply_divergent_deformation`) rather than growing new nodes here."""
    return lithosphere.REFERENCE_HC_OCEANIC_M, lithosphere.YOUNG_RIDGE_HM_M


def _runs_of_at_least(mask: np.ndarray, min_run: int) -> np.ndarray:
    """`mask`, with every True-run shorter than `min_run` cleared to False. Used to gate
    continental-edge retreat on a genuine multi-node contested stretch rather than a
    single stray envelope-fuzz node (see CONTINENTAL_CONTESTED_RETREAT_MIN_RUN). Runs are taken
    in the plate's concatenated node order -- a run that happens to bridge two lines' worth of
    nodes is astronomically rare (line breaks sit at a plate's theta extremes) and harmless
    if it ever happens, since `_grow_or_shrink_line_for_deform` re-checks per line anyway."""
    if min_run <= 1 or not mask.any():
        return mask.copy()
    edges = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    long_enough = (ends - starts) >= min_run
    out = np.zeros_like(mask)
    for start, end in zip(starts[long_enough], ends[long_enough]):
        out[start:end] = True
    return out


def _dilate_1d(mask: np.ndarray, width: int) -> np.ndarray:
    """`mask` grown by `width` positions on each side, within the 1-D node order (used per
    line, so no wrap). `width <= 0` returns an unchanged copy. Backs the collision-uplift
    *reach* knob: a wider contested band -> the orogenic thickening spreads into a broader
    belt, the same "how far inland does a collision crumple crust" lever v1 had as
    COLLISION_RANGE_RAD."""
    if width <= 0 or not mask.any():
        return mask.copy()
    out = mask.copy()
    for shift in range(1, width + 1):
        out[shift:] |= mask[:-shift]
        out[:-shift] |= mask[shift:]
    return out


# The collision-uplift *reach* knob (World.collision_uplift_reach_multiplier) dilates the
# contested band feeding the orogenic thickening by this many nodes per unit of multiplier
# above 1.0 (so reach 3x -> +4 nodes each side of every contested stretch, at the default
# node density); below 1.0 it instead scales the thickening strength down. Reach exactly 1.0
# is a no-op either way.
COLLISION_REACH_DILATION_NODES_PER_UNIT = 2
# Near-field (dilated-but-not-contested) nodes thicken at this fraction of the contested
# rate -- a collision belt's deformation fades outward from the suture, it doesn't step.
COLLISION_REACH_NEAR_FIELD_FACTOR = 0.4


def _redistribute_accreted_column(
    persistent_fields: dict[str, np.ndarray],
    elevation: np.ndarray,
    rho_c: float,
    removed_hc: np.ndarray,
    accrete_removed: np.ndarray,
    from_high: bool,
) -> None:
    """Conserve the crustal volume of the continental-suture nodes just dropped from a line
    end (`removed_hc`, restricted to the `accrete_removed` subset) by thrusting it back onto
    the `SUTURE_ACCRETION_SPREAD_NODES` surviving nodes nearest that same end -- the attached
    mantle lithosphere thickening in proportion, and each node's elevation bumped by the
    isostatic delta. Mutates `persistent_fields`' Hc/Hm arrays and `elevation` in place.

    No-op when nothing dropped was flagged for accretion -- a passive-margin retreat against
    an *oceanic* neighbour leaves `accrete_removed` all-False, and that column is genuinely
    subducted, not preserved. Node area is constant per node, so summed Hc *is* the conserved
    volume (see SUTURE_ACCRETION_SPREAD_NODES / SUTURE_ACCRETION_MAX_HC_M)."""
    if not np.any(accrete_removed):
        return
    add_hc = float(np.sum(removed_hc[accrete_removed]))
    hc = persistent_fields["crustal_thickness_m"]
    hm = persistent_fields["mantle_lithosphere_thickness_m"]
    n = len(hc)
    if n == 0 or add_hc <= 0.0:
        return
    k = min(SUTURE_ACCRETION_SPREAD_NODES, n)
    idx = np.arange(n - k, n) if from_high else np.arange(k)
    before = lithosphere.isostatic_elevation(hc[idx], hm[idx], rho_c)
    # Crustal shortening drags the attached mantle lithosphere along in proportion (same as
    # rheology.apply_convergent_deformation). Hc is capped at SUTURE_ACCRETION_MAX_HC_M -- the
    # overflow delaminates (see the constant) -- and Hm thickens by whatever fraction Hc
    # actually grew after that cap.
    new_hc = np.minimum(hc[idx] + add_hc / k, SUTURE_ACCRETION_MAX_HC_M)
    hm[idx] *= new_hc / hc[idx]
    hc[idx] = new_hc
    after = lithosphere.isostatic_elevation(hc[idx], hm[idx], rho_c)
    elevation[idx] = rheology.clip_elevation_bounds(elevation[idx] + (after - before))


class LithospherePlate(PlateWithLines):
    """A `PlateWithLines` whose per-node state is a lithospheric column (Hc/Hm) rather than
    an independently-set elevation -- see elevation_lines.py's own note on the two new
    OPTIONAL_FIELDS this relies on."""

    def crust_density(self) -> float:
        return lithosphere.crust_density(self.crust_type)

    # -- Motion: torque.py's real implementation -----------------------------------------

    def shift(self, world: "World", years: float) -> float:  # noqa: F821 (World only for typing)
        other_plates = [p for p in world.plates if p.plate_id != self.plate_id]
        return torque.shift_plate(self, world, other_plates, years)

    # -- Deformation: rheology.py's Mohr-Coulomb/isostasy update ---------------------------

    def deform(self, world: "World", other_plates: list, years: float, max_distance: float) -> None:  # noqa: F821
        own_points, _ = self.all_points_and_elevation()
        if not self.lines or len(own_points) == 0:
            return

        # Volume-budget growth gate -- see CONTINENTAL_AREA_BUDGET_MULT. Continental crust
        # only: oceanic footprint is already bounded by subduction. Over budget -> this step
        # grows no new areal crust (end-growth below and `_claim_adjacent_territory`), but
        # still retreats / thins / thickens toward the budget.
        suppress_growth = False
        if self.crust_type == "continental":
            hc_all = self.collect("crustal_thickness_m")
            n_continental = int(np.count_nonzero(hc_all >= CONTINENTAL_BUDGET_HC_FRACTION * lithosphere.REFERENCE_HC_CONTINENTAL_M))
            suppress_growth = len(own_points) > CONTINENTAL_AREA_BUDGET_MULT * n_continental

        spacing_rad = line_spacing_rad(world.node_density)
        reach_rad = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing_rad
        extend_threshold_rad = EXTEND_THRESHOLD_MULTIPLIER * spacing_rad
        max_extend_nodes = max(1, round(MAX_EXTEND_NODES_PER_STEP * np.sqrt(world.node_density)))

        neighbours = self.get_neighbours(other_plates, threshold_rad=reach_rad)
        inputs = torque.gather_boundary_force_inputs(self, neighbours, spacing_rad, reach_rad)
        # Motion-based: `convergent_all` is the whole converging band (not just the nodes
        # that already overlap a neighbour polygon), so a boundary builds an orogen before
        # any overlap accumulates; `contested_all` (the geometric overlap subset, folded into
        # `convergent_all`) still gates node deletion / continental retreat below.
        convergent_all, divergent_all, transform_all, contested_all = torque.classify_boundary_nodes(
            self, neighbours, inputs, reach_rad
        )

        # Fault-localised deformation (World.fault_deformation_mode == "fault"): scale this
        # step's convergent thickening and divergent thinning by proximity to an active fault
        # trace, so plate-boundary transformation concentrates onto fault lines instead of a
        # smooth band at the polygon edge. `fault_influence` is all-ones (i.e. a no-op) in
        # every other mode, when the plate has no active fault, or before the first fault has
        # spawned in a fresh contested zone -- Piece-1 overlap spawning fills those in within
        # a step or two. Deliberately NOT applied to the arc band below: a volcanic arc is a
        # genuinely broad magmatic swath, not a fault-localised structure.
        if getattr(world, "fault_deformation_mode", "fault") == "fault":
            from . import faults

            fault_influence_all = faults.fault_influence(world, self, own_points)
        else:
            fault_influence_all = np.ones(len(own_points))

        neighbor_omega_all = inputs.neighbor_omega
        closing_rate_all = rheology.normal_closing_rate_m_per_s(self.omega, neighbor_omega_all, own_points, inputs.direction_to_neighbor)

        # What may retreat this step. Oceanic crust: any contested node subducts. Continental
        # crust: any contested end-node in a run of >= CONTINENTAL_CONTESTED_RETREAT_MIN_RUN
        # consecutive contested nodes -- whether the overriding neighbour is oceanic (passive
        # margin) or continental (a suture whose overlap is consumed into the orogen, the
        # retreated column's volume thrust onto the plate's own leading edge -- see
        # _redistribute_accreted_column). Envelope fuzz (a lone contested node) still can't
        # nibble a stable margin, and the interior carve below stays oceanic-only so a
        # continental row is never severed mid-line. See CONTINENTAL_CONTESTED_RETREAT_MIN_RUN
        # for the ratchet / frozen-overlap this breaks.
        if self.crust_type != "continental":
            shrinkable_all = contested_all
        else:
            shrinkable_all = _runs_of_at_least(contested_all, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN)

        # Continental suture retreat conserves the consumed column's volume by accreting it
        # onto this plate's own leading edge (_redistribute_accreted_column); a retreat where
        # the overriding neighbour is *oceanic* does not -- that column subducts and is lost.
        # Oceanic self-plates never accrete.
        if self.crust_type == "continental":
            accrete_all = shrinkable_all & ~inputs.neighbor_is_oceanic
        else:
            accrete_all = np.zeros_like(shrinkable_all)

        # Continental arc band: this plate's own nodes within `reach_rad` of a *converging
        # oceanic* neighbour -- the volcanic arc + accreted forearc / underplated wedge sits
        # inboard of the trench, a swath (~500 km at default density), not just the contact
        # line (which is only a few tens of nodes -- far too narrow to counter the land
        # decline). `arc_intensity_all` fades from 1 at the contact to ~0.3 at the band edge.
        # Feeds both the magmatic Hc thickening (below) and the arc-crust growth seed
        # (`arc_end_*` -> `_grow_or_shrink_line_for_deform`). See ARC_MARGIN_SEED_HC_M.
        arc_band_all = np.zeros(len(own_points), dtype=bool)
        arc_intensity_all = np.zeros(len(own_points))
        if self.crust_type == "continental":
            arc_band_all = (
                inputs.neighbor_is_oceanic
                & np.isfinite(inputs.dist_to_neighbor)
                & (closing_rate_all > rheology.ARC_MIN_CONVERGENCE_M_PER_S)
            )
            arc_intensity_all = np.where(
                arc_band_all, np.clip(1.0 - 0.7 * (inputs.dist_to_neighbor / reach_rad), 0.3, 1.0), 0.0
            )

        years_myr = years / 1_000_000.0
        rho_c = self.crust_density()

        # Collision-uplift tuning knobs (the "Controls" window, 1.0 == untuned -- see World).
        # `orogen_amount` scales the plastic thickening rate at contested nodes; `orogen_reach`
        # widens (>1) or narrows (<1) the belt it acts on -- see _dilate_1d /
        # COLLISION_REACH_*. Both exactly 1.0 leave apply_convergent_deformation's strength at
        # 1.0 over precisely the contested set, i.e. byte-identical to before the knobs.
        orogen_amount = world.collision_uplift_multiplier
        orogen_reach = world.collision_uplift_reach_multiplier
        orogen_contested_strength = orogen_amount * min(orogen_reach, 1.0)
        orogen_dilation_nodes = (
            round((orogen_reach - 1.0) * COLLISION_REACH_DILATION_NODES_PER_UNIT)
            if orogen_reach > 1.0 and self.crust_type == "continental"
            else 0
        )

        fault_noise = (
            SphereNoise(np.random.default_rng((world.seed, self.plate_id, 9001)), octaves=3, base_freq=9.0)
            if self.crust_type == "continental"
            else None
        )

        new_lines: list[ElevationLine] = []
        offset = 0
        for line_index, line in enumerate(self.lines):
            n = len(line)
            sl = slice(offset, offset + n)
            offset += n

            contested = contested_all[sl]
            convergent = convergent_all[sl]
            divergent = divergent_all[sl]
            transform = transform_all[sl]
            shrinkable = shrinkable_all[sl]
            accrete = accrete_all[sl]
            closing_rate = closing_rate_all[sl]
            neighbor_oceanic = inputs.neighbor_is_oceanic[sl]
            arc_band = arc_band_all[sl]
            arc_intensity = arc_intensity_all[sl]
            fault_influence = fault_influence_all[sl]  # all-ones except in "fault" mode

            # Active-margin growth seed per line end -- see ARC_MARGIN_SEED_HC_M. An end is an
            # active margin if a node within ARC_MARGIN_END_SCAN_NODES of it is in the arc
            # band, or still carries a subduction-arc provenance stamp from a recent step (the
            # ocean's edge can retreat a step before this plate's edge grows into the gap).
            arc_end_low = arc_end_high = False
            if self.crust_type == "continental" and n > 0:
                arc_signal = arc_band | (line.elev_change_reason == ELEV_CHANGE_SUBDUCTION_ARC)
                k = ARC_MARGIN_END_SCAN_NODES
                arc_end_low = bool(arc_signal[:k].any())
                arc_end_high = bool(arc_signal[-k:].any())

            hc = line.crustal_thickness_m.copy()
            hm = line.mantle_lithosphere_thickness_m.copy()
            # Isostasy-driven elevation change is applied as a *delta* on top of whatever
            # elevation already holds (elevation_before -> below), not a wholesale overwrite
            # -- erosion.py (run later this same step_world call, and every step
            # thereafter until the next deform()) mutates `elevation` directly, with no
            # notion of Hc/Hm at all. An unconditional overwrite here would silently erase
            # every step's worth of erosion the instant the *next* deform() call ran,
            # confirmed directly as a real bug (a 3-step run's own elevation stopped
            # matching isostasy(Hc, Hm) exactly the way an unconditional-overwrite design
            # would have predicted, because erosion's own contribution was still baked into
            # the *un-clipped* portion of `elevation` between tectonic uplift events -- the
            # fix is this delta, not forcing elevation back to a bare isostasy readout).
            rho_c = self.crust_density()
            elevation_before = lithosphere.isostatic_elevation(hc, hm, rho_c)

            # The band that plastically thickens: the whole converging band at
            # `orogen_contested_strength`, plus (reach knob > 1) a dilated near-field ring at
            # a faded rate. `orogen_strength` is the per-node multiplier handed to
            # apply_convergent_deformation; > 0 exactly on the nodes that thicken.
            # `apply_convergent_deformation` still gates on each node's own closing rate
            # (below yield / not actually closing -> zero strain), so a node that is
            # `convergent` only via the `contested` deep-overlap fold and is no longer
            # actively closing simply thickens at zero.
            near_field = (
                _dilate_1d(convergent, orogen_dilation_nodes) & ~convergent & ~divergent
                if orogen_dilation_nodes > 0
                else np.zeros(n, dtype=bool)
            )
            orogen_strength = np.where(convergent, orogen_contested_strength, 0.0)
            orogen_strength[near_field] = orogen_amount * COLLISION_REACH_NEAR_FIELD_FACTOR
            # "fault" mode: concentrate the shortening onto fault traces (no-op / all-ones
            # otherwise). `strength` scales apply_convergent_deformation's thickening rate.
            orogen_strength = orogen_strength * fault_influence
            thicken = orogen_strength > 0.0
            if np.any(thicken):
                fault_factor = (
                    np.where(
                        fault_noise.sample(geometry.local_xyz(np.full(n, line.phi), line.theta)) < -0.15,
                        rheology.REVERSE_FAULT_VALLEY_UPLIFT_FACTOR,
                        1.0,
                    )
                    if fault_noise is not None
                    else np.ones(n)
                )
                # The overlapping crust a continent-continent suture retreats over is not
                # lost here via a `fault_factor` boost -- its actual volume is conserved and
                # thrust onto the leading edge in `_grow_or_shrink_line_for_deform` (see
                # `_redistribute_accreted_column`). This path is just the ordinary
                # yield-limited plastic thickening.
                new_hc, new_hm = rheology.apply_convergent_deformation(
                    hc[thicken], hm[thicken], closing_rate[thicken], years_myr,
                    fault_factor[thicken], strength=orogen_strength[thicken],
                )
                hc[thicken] = new_hc
                hm[thicken] = new_hm

            # Continental arc magmatism: an oceanic slab subducting under this margin fluxes
            # the mantle wedge and underplates juvenile crust across the whole arc band --
            # extra Hc (added from the mantle, not conserved), the crust-building half of
            # "subduction under a continent makes more continent" (docs/TODO.md "Land fraction
            # slowly declines"). Separate from the contested shortening above: the band is far
            # wider than the contact line. Bounded long-term by the CONTINENTAL_AREA_BUDGET_MULT
            # volume gate.
            if np.any(arc_band):
                hc[arc_band], hm[arc_band] = rheology.apply_arc_magmatic_thickening(
                    hc[arc_band], hm[arc_band], closing_rate[arc_band], years_myr, arc_intensity[arc_band]
                )

            melting = np.zeros(n, dtype=bool)
            if np.any(divergent):
                new_hc, new_hm, melt = rheology.apply_divergent_deformation(hc[divergent], hm[divergent], closing_rate[divergent], years_myr)
                # "fault" mode: scale the thinning delta by fault proximity (all-ones
                # otherwise). Melt (decompression volcanism) still fires on the geometric
                # rift threshold -- it's a discrete event, not a rate.
                infl = fault_influence[divergent]
                hc[divergent] = hc[divergent] + infl * (new_hc - hc[divergent])
                hm[divergent] = hm[divergent] + infl * (new_hm - hm[divergent])
                melting[divergent] = melt

            prior_age = line.divergent_age_myr
            new_age = np.where(divergent, prior_age + years_myr, 0.0)
            if self.crust_type == "oceanic":
                hm = rheology.relax_young_oceanic_mantle_lithosphere(hm, new_age, years_myr)

            is_volcano = line.is_volcano.copy()
            volcano_remaining = line.volcano_active_years_remaining.copy()
            if np.any(melting):
                # Decompression melting (spec 2.3): a rift that just thinned past the
                # critical threshold erupts fresh oceanic crust in place -- same one-
                # guaranteed-eruption convention v1's stretch-volcano growth used.
                from .elevation_lines import ERUPTION_ELEVATION_M, VOLCANO_ACTIVE_MAX_YEARS, VOLCANO_ACTIVE_MIN_YEARS

                hc[melting] = lithosphere.REFERENCE_HC_OCEANIC_M
                hm[melting] = lithosphere.YOUNG_RIDGE_HM_M
                is_volcano[melting] = True
                rng = np.random.default_rng((world.seed, round(world.elapsed_years), self.plate_id, line_index))
                volcano_remaining[melting] = rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=int(melting.sum()))

            # Transform (strike-slip) pressure-ridge uplift: a modest, always-transpressional
            # bump on the transform band, kept as a direct elevation delta (like erosion's
            # own contributions) rather than an Hc change -- a strike-slip contact shoulders
            # up local relief without net crustal shortening. Gated by `fault_influence` in
            # "fault" mode so it tracks the boundary strike-slip fault families rather than
            # smearing along the whole polygon edge.
            transform_uplift = np.zeros(n)
            transform_uplift[transform] = (
                TRANSFORM_UPLIFT_RATE_M_PER_MYR * years_myr * fault_influence[transform]
            )

            elevation_after = lithosphere.isostatic_elevation(hc, hm, rho_c)
            new_elevation = rheology.clip_elevation_bounds(
                line.elevation + (elevation_after - elevation_before) + transform_uplift
            )

            # Elevation-change provenance (diagnostic only -- see elevation_lines.ELEV_CHANGE_*
            # and render_image's "elevReason" view). Stamp whichever tectonic process moved a
            # node this step, gated on ELEV_CHANGE_MIN_DELTA_M so a node barely grazed by a
            # fading boundary force keeps its older provenance. The masks partition the
            # near-boundary band by motion (convergent / divergent / transform), so a plain
            # per-mask assignment needs no priority order. `faults._apply_plate_fault_relief`
            # runs after this pass and overwrites these with a FAULT_* code wherever a
            # boundary fault of the matching regime moved the node -- that is what paints the
            # fault families along every boundary in the elevReason view.
            reason = line.elev_change_reason.copy()
            moved = np.abs(new_elevation - line.elevation) >= ELEV_CHANGE_MIN_DELTA_M
            if self.crust_type == "continental":
                # near_field (the reach knob's dilated ring) is continent-continent orogenic
                # belt too, so it carries the same COLLISION provenance as the converging core.
                reason[(convergent | near_field) & moved & ~neighbor_oceanic] = ELEV_CHANGE_COLLISION
                reason[convergent & moved & neighbor_oceanic] = ELEV_CHANGE_SUBDUCTION_ARC
                reason[arc_band & moved] = ELEV_CHANGE_SUBDUCTION_ARC
            else:
                reason[convergent & moved] = ELEV_CHANGE_TRENCH
            reason[divergent & moved] = ELEV_CHANGE_RIFT
            reason[transform & moved] = ELEV_CHANGE_TRANSFORM
            reason[melting] = ELEV_CHANGE_VOLCANO

            updated_line = line.replace(
                elevation=new_elevation,
                crustal_thickness_m=hc,
                mantle_lithosphere_thickness_m=hm,
                divergent_age_myr=new_age,
                is_volcano=is_volcano,
                volcano_active_years_remaining=volcano_remaining,
                elev_change_reason=reason,
            )
            grown_lines = self._grow_or_shrink_line_for_deform(
                updated_line,
                inputs.dist_to_neighbor[sl],
                contested,
                shrinkable,
                accrete,
                spacing_rad,
                extend_threshold_rad,
                max_extend_nodes,
                max_distance,
                world,
                line_index,
                neighbours,
                suppress_growth,
                arc_end_low,
                arc_end_high,
            )
            new_lines.extend(gl for gl in grown_lines if len(gl) > 0)

        if self.crust_type == "continental":
            new_lines = self._retreat_contested_leading_rows(new_lines, contested_all, years)

        self.set_lines(new_lines)
        if not suppress_growth:
            self._claim_adjacent_territory(world, neighbours, spacing_rad)

        for line_index, line in enumerate(self.lines):
            if needs_regularizing(line, spacing_rad):
                self.replace_line(line_index, regularize_line(line, spacing_rad))

    def _count_open_prefix(self, theta_candidates: np.ndarray, phi: float, neighbours: list) -> int:
        if len(theta_candidates) == 0 or not neighbours:
            return len(theta_candidates)
        world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.full_like(theta_candidates, phi), theta_candidates))
        contested = _contested_by_any(world_pts, neighbours)
        first_contested = np.argmax(contested) if np.any(contested) else len(contested)
        return int(first_contested)

    def _grow_or_shrink_line_for_deform(
        self,
        line: ElevationLine,
        dist: np.ndarray,
        contested: np.ndarray,
        shrinkable: np.ndarray,
        accrete: np.ndarray,
        spacing_rad: float,
        extend_threshold_rad: float,
        max_extend_nodes: int,
        max_distance: float,
        world: "World",  # noqa: F821
        line_index: int,
        neighbours: list,
        suppress_growth: bool = False,
        arc_end_low: bool = False,
        arc_end_high: bool = False,
    ) -> list[ElevationLine]:
        """Same grow/shrink shape as `PlateWithLines._grow_or_shrink_line_for_deform` (see
        that method's own docstring -- end-only growth/shrink, plus the oceanic-only
        interior-subduction carve-out that can return a row as two contiguous
        `ElevationLine`s) -- reimplemented rather than inherited only because `grow_end`
        below needs to seed fresh Hc/Hm columns instead of a flat elevation target.
        Shrinking (end and interior) is generic over every `ElevationLine.OPTIONAL_FIELDS`
        name already (Hc/Hm included, since they're threaded through `OPTIONAL_FIELDS` --
        see elevation_lines.py), so only growth needed a new Hc/Hm-aware body.

        `accrete` marks end nodes whose crustal/mantle-lithosphere volume must be conserved
        when they retreat (a continental suture -- see `_redistribute_accreted_column`);
        elsewhere retreat drops the column (oceanic subduction, or a continental passive
        margin against an oceanic slab)."""
        theta = line.theta.copy()
        elevation = line.elevation.copy()
        contested = contested.copy()
        shrinkable = shrinkable.copy()
        accrete = accrete.copy()
        dist = dist.copy()
        rho_c = self.crust_density()
        persistent_fields = {name: getattr(line, name).copy() for name in ElevationLine.OPTIONAL_FIELDS}
        if len(theta) == 0:
            return [ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)]

        dtheta = spacing_rad / max(np.cos(line.phi), 1e-3)
        n_distance_cap = max(1, int(max_distance / spacing_rad))

        # A row is a circle of local latitude -- its theta extent physically cannot exceed a
        # full revolution. Nothing here treats theta as periodic, so once end-growth has
        # closed the loop the "gap to nearest neighbour is wide open" test stays true forever
        # near a plate's own local pole (the pole cap belongs to nobody) and the row just
        # keeps winding. `ring_room()` is how many more `dtheta` nodes an end can take before
        # the row spans 2*pi; growth is capped by it, and at zero the end stops. Ported from
        # `PlateWithLines._grow_or_shrink_line_for_deform` (the v1 pole-winding fix) -- this
        # v2 override predates that fix and, without this, relied entirely on
        # `regularize_line`'s after-the-fact unwind, so rows still over-wound by up to a
        # revolution every step and were unwound the next (continuous churn, and near-pole
        # rings feeding overlap / node count).
        full_revolution_span = 2.0 * np.pi - _ROW_FULL_REVOLUTION_SLACK * dtheta

        def ring_room() -> int:
            if len(theta) < 2:
                return n_distance_cap
            return int(np.floor((full_revolution_span - (theta[-1] - theta[0])) / dtheta))

        def contested_run_from_end(mask: np.ndarray, from_high: bool) -> int:
            ordered = mask[::-1] if from_high else mask
            run = 0
            for value in ordered:
                if not value:
                    break
                run += 1
            return run

        if len(shrinkable) > 0 and shrinkable[-1]:
            n_remove = min(contested_run_from_end(shrinkable, from_high=True), n_distance_cap, max_extend_nodes, len(theta) - 1)
            if n_remove > 0:
                removed_hc = persistent_fields["crustal_thickness_m"][-n_remove:].copy()
                accrete_removed = accrete[-n_remove:].copy()
                theta, elevation = theta[:-n_remove], elevation[:-n_remove]
                contested, shrinkable, accrete, dist = contested[:-n_remove], shrinkable[:-n_remove], accrete[:-n_remove], dist[:-n_remove]
                persistent_fields = {name: values[:-n_remove] for name, values in persistent_fields.items()}
                _redistribute_accreted_column(persistent_fields, elevation, rho_c, removed_hc, accrete_removed, from_high=True)

        if len(theta) == 0:
            return [ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)]

        if shrinkable[0]:
            n_remove = min(contested_run_from_end(shrinkable, from_high=False), n_distance_cap, max_extend_nodes, len(theta) - 1)
            if n_remove > 0:
                removed_hc = persistent_fields["crustal_thickness_m"][:n_remove].copy()
                accrete_removed = accrete[:n_remove].copy()
                theta, elevation = theta[n_remove:], elevation[n_remove:]
                contested, shrinkable, accrete, dist = contested[n_remove:], shrinkable[n_remove:], accrete[n_remove:], dist[n_remove:]
                persistent_fields = {name: values[n_remove:] for name, values in persistent_fields.items()}
                _redistribute_accreted_column(persistent_fields, elevation, rho_c, removed_hc, accrete_removed, from_high=False)

        if len(theta) == 0:
            return [ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)]

        # Interior subduction: carve out each substantial mid-row `shrinkable` run the two
        # end-shrinks can't reach, splitting the row into separate arcs -- see
        # `PlateWithLines._grow_or_shrink_line_for_deform` for the full rationale. Oceanic
        # self-plate only: a continental row's oceanic-contested nodes (now `shrinkable`, see
        # deform) must only ever retreat from the ends -- carving a continental row's middle
        # would sever the landmass into a spurious defragmentation plate.
        if self.crust_type == "oceanic" and len(shrinkable) >= _INTERIOR_SUBDUCTION_MIN_RUN + 2 and shrinkable[1:-1].any():
            prev_shrink = np.concatenate([[False], shrinkable[:-1]])
            run_starts = np.nonzero(shrinkable & ~prev_shrink)[0]
            keep = np.ones(len(theta), dtype=bool)
            budget = max_extend_nodes
            for start in run_starts:
                end = start
                while end < len(shrinkable) and shrinkable[end]:
                    end += 1
                if start == 0 or end >= len(shrinkable):
                    continue
                if end - start < _INTERIOR_SUBDUCTION_MIN_RUN or budget <= 0:
                    continue
                take = min(end - start, budget)
                keep[start : start + take] = False
                budget -= take
            if not keep.all():
                theta, elevation = theta[keep], elevation[keep]
                contested, shrinkable, accrete, dist = contested[keep], shrinkable[keep], accrete[keep], dist[keep]
                persistent_fields = {name: values[keep] for name, values in persistent_fields.items()}

        # Brand-new areal crust at a growing end is normally oceanic regardless of this
        # plate's own type -- see growth_seed_thickness() for the land-area runaway that rule
        # prevents. The one exception is a continental plate's *leading* edge advancing into
        # space a subducting oceanic neighbour is vacating (`arc_end_low` / `arc_end_high`,
        # from deform's active-margin scan): that ground is juvenile arc / accreted-terrane
        # crust, seeded at the thicker ARC_MARGIN_SEED_* column and stamped as a subduction
        # arc. See ARC_MARGIN_SEED_HC_M for why this is safe against the old runaway.
        ocean_hc0, ocean_hm0 = growth_seed_thickness()

        def _end_seed(is_arc: bool) -> tuple[float, float, float, float]:
            hc_seed, hm_seed = (ARC_MARGIN_SEED_HC_M, ARC_MARGIN_SEED_HM_M) if is_arc else (ocean_hc0, ocean_hm0)
            elev_seed = float(lithosphere.isostatic_elevation(np.array([hc_seed]), np.array([hm_seed]), rho_c)[0])
            reason_seed = ELEV_CHANGE_SUBDUCTION_ARC if is_arc else ELEV_CHANGE_NEW_CRUST
            return hc_seed, hm_seed, elev_seed, reason_seed

        def _fill_new_nodes(n_new: int, hc_seed: float, hm_seed: float, reason_seed: float) -> dict[str, np.ndarray]:
            out = {}
            for name, values in persistent_fields.items():
                if name == "crustal_thickness_m":
                    fill = np.full(n_new, hc_seed)
                elif name == "mantle_lithosphere_thickness_m":
                    fill = np.full(n_new, hm_seed)
                elif name == "elev_change_reason":
                    fill = np.full(n_new, reason_seed, dtype=values.dtype)
                else:
                    fill = np.zeros(n_new, dtype=values.dtype)
                out[name] = fill
            return out

        if not suppress_growth and not contested[-1] and dist[-1] > extend_threshold_rad and ring_room() > 0:
            gap_estimate = min(dist[-1], (n_distance_cap + 1) * spacing_rad)
            n_candidates = min(max(int(gap_estimate / spacing_rad), 1), n_distance_cap, max_extend_nodes, ring_room())
            candidate_theta = theta[-1] + dtheta * np.arange(1, n_candidates + 1)
            n_new = self._count_open_prefix(candidate_theta, line.phi, neighbours)
            if n_new > 0:
                hc_seed, hm_seed, elev_seed, reason_seed = _end_seed(arc_end_high)
                new_theta = candidate_theta[:n_new]
                theta = np.append(theta, new_theta)
                elevation = np.append(elevation, np.full(n_new, elev_seed))
                for name, fill in _fill_new_nodes(n_new, hc_seed, hm_seed, reason_seed).items():
                    persistent_fields[name] = np.append(persistent_fields[name], fill)

        if not suppress_growth and not contested[0] and dist[0] > extend_threshold_rad and ring_room() > 0:
            gap_estimate = min(dist[0], (n_distance_cap + 1) * spacing_rad)
            n_candidates = min(max(int(gap_estimate / spacing_rad), 1), n_distance_cap, max_extend_nodes, ring_room())
            candidate_theta = theta[0] - dtheta * np.arange(1, n_candidates + 1)
            n_new = self._count_open_prefix(candidate_theta, line.phi, neighbours)
            if n_new > 0:
                hc_seed, hm_seed, elev_seed, reason_seed = _end_seed(arc_end_low)
                new_theta = candidate_theta[:n_new][::-1]
                theta = np.insert(theta, 0, new_theta)
                elevation = np.insert(elevation, 0, np.full(n_new, elev_seed))
                for name, fill in _fill_new_nodes(n_new, hc_seed, hm_seed, reason_seed).items():
                    persistent_fields[name] = np.insert(persistent_fields[name], 0, fill)

        result = ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)
        return split_into_contiguous_runs(result, dtheta)

    def _retreat_contested_leading_rows(
        self, new_lines: list[ElevationLine], contested_all: np.ndarray, years: float
    ) -> list[ElevationLine]:
        """Drop this plate's outermost phi-row (at either phi extreme) once a neighbour has
        overridden it -- `LEADING_ROW_CONTESTED_FRACTION` of its nodes contested -- for a
        cumulative `LEADING_ROW_RETREAT_SUSTAINED_YEARS` of deform time. The reverse of
        `_claim_adjacent_territory`, and the only retreat op the parallel-suture regime allows.
        See the constant block above.

        `contested_all` is this step's boundary classification in the concatenation order of
        the *pre-grow* `self.lines` (this runs before `set_lines(new_lines)`); the drop is
        applied to `new_lines` by matching `line.phi` (grow/shrink is theta-only, so a row's
        phi is unchanged, and a row split into two arcs by `split_into_contiguous_runs` shares
        one phi and is dropped as a unit). The sustained-time tally lives on the plate
        (`_leading_row_retreat_years`), keyed by which extreme -- it survives a rotation (rows
        are stored plate-local) but resets on merge/split/load, which only delays a drop."""
        tracker: dict[str, tuple[float, float]] = getattr(self, "_leading_row_retreat_years", None)
        if tracker is None:
            tracker = {}
            self._leading_row_retreat_years = tracker

        contested_by_phi: dict[float, list[float]] = {}
        offset = 0
        for line in self.lines:
            n = len(line)
            key = round(float(line.phi), 6)
            agg = contested_by_phi.setdefault(key, [0.0, 0.0])
            agg[0] += n
            agg[1] += float(contested_all[offset : offset + n].sum())
            offset += n

        rows_left = sorted({round(float(ln.phi), 6) for ln in new_lines if len(ln) > 0})
        if len(rows_left) < LEADING_ROW_DROP_MIN_ROWS:
            tracker.clear()
            return new_lines

        drop_phis: list[float] = []
        for extreme, phi_key in (("lo", rows_left[0]), ("hi", rows_left[-1])):
            n_nodes, n_contested = contested_by_phi.get(phi_key, (0.0, 0.0))
            fraction = n_contested / n_nodes if n_nodes else 0.0
            if fraction < LEADING_ROW_CONTESTED_FRACTION:
                tracker.pop(extreme, None)
                continue
            prev_phi, prev_years = tracker.get(extreme, (None, 0.0))
            accumulated = (prev_years if prev_phi == phi_key else 0.0) + years
            if accumulated >= LEADING_ROW_RETREAT_SUSTAINED_YEARS:
                drop_phis.append(phi_key)
                tracker.pop(extreme, None)
            else:
                tracker[extreme] = (phi_key, accumulated)

        if not drop_phis:
            return new_lines
        return [ln for ln in new_lines if round(float(ln.phi), 6) not in drop_phis]

    def _claim_adjacent_territory(self, world: "World", neighbours: list, spacing_rad: float) -> None:  # noqa: F821
        """Same shape as `PlateWithLines._claim_adjacent_territory` -- a brand-new phi row
        just past this plate's own phi extremes, where open -- seeded with fresh Hc/Hm
        (oceanic reference, see `growth_seed_thickness`) plus `terrain_noise.FractalTexture`
        on Hc (an extension of an already-shaped plate, so texture rather than a fresh
        orogen), rather than a flat elevation baseline. Keyed off `(world.seed, plate_id,
        _TERRAIN_SEED_TAG)` so the texture stays attached to this plate as it grows."""
        lines_with_nodes = [line for line in self.lines if len(line) > 0]
        if not lines_with_nodes:
            return
        ordered = sorted(lines_with_nodes, key=lambda line: line.phi)
        # Keep POLE_CAP_MARGIN_MULT target spacings clear of the local pole -- see the v1
        # POLE_CAP_MARGIN_MULT comment. Right at +-pi/2 a row's theta step (spacing / cos phi)
        # blows up and the row degenerates into a handful of sub-spacing rings that read as
        # concentric circles / holes and feed the theta-winding pathology the ring_room cap in
        # `_grow_or_shrink_line_for_deform` now guards against. This v2 override predated the
        # v1 fix and still marched a plate right onto its pole (spacing_rad / 2).
        max_phi_limit = np.pi / 2 - POLE_CAP_MARGIN_MULT * spacing_rad
        # Oceanic regardless of self.crust_type -- new sphere area is floored by sea-floor
        # spreading, not by this plate's own crust (see growth_seed_thickness()).
        hc0, hm0 = growth_seed_thickness()
        amp = hc0 * 0.1  # texture on fresh Hc, same spirit as v1's noise-on-elevation
        texture = terrain_noise.FractalTexture(
            np.random.default_rng((world.seed, self.plate_id, _TERRAIN_SEED_TAG))
        )
        new_lines: list[ElevationLine] = []

        for reference, direction in ((ordered[0], -1), (ordered[-1], 1)):
            new_phi = reference.phi + direction * spacing_rad
            if abs(new_phi) > max_phi_limit:
                continue
            dtheta = spacing_rad / max(np.cos(new_phi), 1e-3)
            span = reference.theta[-1] - reference.theta[0]
            n_cols = max(int(round(span / dtheta)) + 1, 1)
            theta_candidates = reference.theta[0] + dtheta * np.arange(n_cols)
            world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.full(n_cols, new_phi), theta_candidates))

            contested = _contested_by_any(world_pts, neighbours)
            open_mask = ~contested
            if not np.any(open_mask):
                continue

            theta_open = theta_candidates[open_mask]
            n_open = int(open_mask.sum())
            hc_open = np.full(n_open, hc0) + amp * texture.sample(world_pts[open_mask])
            hm_open = np.full(n_open, hm0)
            elevation_open = lithosphere.isostatic_elevation(hc_open, hm_open, self.crust_density())
            new_lines.append(
                ElevationLine(
                    phi=new_phi,
                    theta=theta_open,
                    elevation=elevation_open,
                    crustal_thickness_m=hc_open,
                    mantle_lithosphere_thickness_m=hm_open,
                    elev_change_reason=np.full(n_open, ELEV_CHANGE_NEW_CRUST, dtype=float),
                )
            )

        if new_lines:
            self.set_lines(list(self.lines) + new_lines)

    # -- Merge/split: carry Hc/Hm through, not just elevation -------------------------------

    def merge_with(self, other: "LithospherePlate", spacing_rad: float, coverage_radius_rad: float, other_points_xyz: np.ndarray) -> None:
        own_points, _ = self.all_points_and_elevation()
        other_points, _ = other.all_points_and_elevation()
        inertia_self = lithosphere.moment_of_inertia_tensor(
            own_points, self.collect("crustal_thickness_m"), self.collect("mantle_lithosphere_thickness_m"), self.crust_density(), spacing_rad
        )
        inertia_other = lithosphere.moment_of_inertia_tensor(
            other_points, other.collect("crustal_thickness_m"), other.collect("mantle_lithosphere_thickness_m"), other.crust_density(), spacing_rad
        )
        self._merge_nodes_with(other, spacing_rad, coverage_radius_rad, other_points_xyz)
        self.set_omega(torque.merge_omega(self, inertia_self, other, inertia_other))
        self.reset_age()

    def _merge_nodes_with(self, other: "LithospherePlate", spacing_rad: float, coverage_radius_rad: float, other_points_xyz: np.ndarray) -> None:
        keep_pts, _ = self.all_points_and_elevation()
        absorb_pts, _ = other.all_points_and_elevation()
        old_points = np.concatenate([keep_pts, absorb_pts], axis=0)
        old_hc = np.concatenate([self.collect("crustal_thickness_m"), other.collect("crustal_thickness_m")])
        old_hm = np.concatenate([self.collect("mantle_lithosphere_thickness_m"), other.collect("mantle_lithosphere_thickness_m")])
        exclude_tree = cKDTree(other_points_xyz) if len(other_points_xyz) else None
        self.set_lines(_lines_from_resample(self.frame, old_points, old_hc, old_hm, coverage_radius_rad, spacing_rad, exclude_tree))
        lithosphere.sync_plate_elevation(self)

    def grow_into(self, new_points_xyz: np.ndarray, new_elevation: np.ndarray, coverage_radius_rad: float, spacing_rad: float) -> None:
        old_points, _ = self.all_points_and_elevation()
        old_hc = self.collect("crustal_thickness_m")
        old_hm = self.collect("mantle_lithosphere_thickness_m")
        hc0, hm0 = lithosphere.reference_thickness(self.crust_type)
        combined_points = np.concatenate([old_points, new_points_xyz], axis=0)
        combined_hc = np.concatenate([old_hc, np.full(len(new_points_xyz), hc0)])
        combined_hm = np.concatenate([old_hm, np.full(len(new_points_xyz), hm0)])
        self.set_lines(_lines_from_resample(self.frame, combined_points, combined_hc, combined_hm, coverage_radius_rad, spacing_rad))
        lithosphere.sync_plate_elevation(self)

    def split(self, new_id: int, cut_normal: np.ndarray, min_nodes: int) -> tuple["LithospherePlate", "LithospherePlate"] | None:
        lines_a: list[ElevationLine] = []
        lines_b: list[ElevationLine] = []
        for line in self.lines:
            world_pts = line.world_xyz(self.frame)
            side = np.sum(world_pts * cut_normal, axis=-1) > 0
            # One contiguous ElevationLine per arc -- see PlateWithLines.split's own docstring
            # for why a great-circle cut can otherwise strand a row as two arcs, and why
            # carrying that whole makes the two daughters' envelopes overlap.
            ref = _row_median_step(line)
            if np.any(side):
                lines_a.extend(split_into_contiguous_runs(line.masked(side), ref))
            if np.any(~side):
                lines_b.extend(split_into_contiguous_runs(line.masked(~side), ref))

        if sum(len(l) for l in lines_a) < min_nodes or sum(len(l) for l in lines_b) < min_nodes:
            return None

        plate_a = LithospherePlate(plate_id=self.plate_id, frame=self.frame.copy(), crust_type=self.crust_type, lines=lines_a)
        plate_b = LithospherePlate(plate_id=new_id, frame=self.frame.copy(), crust_type=self.crust_type, lines=lines_b)
        return plate_a, plate_b

    def apply_failed_rift(self, cut_normal: np.ndarray, spacing_rad: float) -> None:
        """A rift that started but *aborted* (`merge_split.RIFT_SUCCESS_PROBABILITY`): the
        plate does not break up, but the stretched zone along the would-be cut is left as a
        thinned continental sag basin (an aulacogen -- the North Sea, the Benue Trough), not
        healed back to full thickness and not oceanised. Thins Hc/Hm by up to
        `FAILED_RIFT_THINNING_FRACTION` within `FAILED_RIFT_BAND_MULT` spacings of the cut
        great circle, tapering to zero at the band edge, and books the isostatic subsidence as
        a delta on `elevation` (the same idiom `deform` uses so erosion isn't clobbered).
        This is a one-off event -- far less crust lost than the sustained divergent thinning +
        decompression-melting a *successful* rift would inflict on both daughters' margins."""
        from .merge_split import FAILED_RIFT_BAND_MULT, FAILED_RIFT_THINNING_FRACTION

        # A well-formed cut plane's normal is a unit vector; a degenerate one (the two flow
        # clusters were spatially intermingled, so `normalize(centroid_a - centroid_b)`
        # collapsed toward zero -- a known `maybe_split_plate` failure mode, see the
        # pole-winding notes in docs/TODO.md) would put every node "next to the rift" and thin
        # the whole plate. No cut, no aulacogen -- the rift just fails silently.
        if not np.isfinite(cut_normal).all() or abs(np.linalg.norm(cut_normal) - 1.0) > 1e-3:
            return
        band_sin = float(np.sin(FAILED_RIFT_BAND_MULT * spacing_rad))
        rho_c = self.crust_density()
        new_lines: list[ElevationLine] = []
        for line in self.lines:
            if len(line) == 0:
                new_lines.append(line)
                continue
            dist_to_plane = np.abs(line.world_xyz(self.frame) @ cut_normal)
            in_band = dist_to_plane < band_sin
            if not np.any(in_band):
                new_lines.append(line)
                continue
            hc = line.crustal_thickness_m.copy()
            hm = line.mantle_lithosphere_thickness_m.copy()
            z_before = lithosphere.isostatic_elevation(hc, hm, rho_c)
            taper = np.clip(1.0 - dist_to_plane / band_sin, 0.0, 1.0)
            factor = 1.0 - FAILED_RIFT_THINNING_FRACTION * taper
            hc[in_band] = np.maximum(hc[in_band] * factor[in_band], lithosphere.MIN_CRUSTAL_THICKNESS_M)
            hm[in_band] = np.maximum(hm[in_band] * factor[in_band], lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M)
            z_after = lithosphere.isostatic_elevation(hc, hm, rho_c)
            new_elevation = rheology.clip_elevation_bounds(line.elevation + (z_after - z_before))
            reason = line.elev_change_reason.copy()
            moved = np.abs(new_elevation - line.elevation) >= ELEV_CHANGE_MIN_DELTA_M
            reason[in_band & moved] = ELEV_CHANGE_RIFT
            new_lines.append(
                line.replace(
                    elevation=new_elevation,
                    crustal_thickness_m=hc,
                    mantle_lithosphere_thickness_m=hm,
                    elev_change_reason=reason,
                )
            )
        self.set_lines(new_lines)


def _lines_from_resample(
    frame: np.ndarray,
    points: np.ndarray,
    hc: np.ndarray,
    hm: np.ndarray,
    coverage_radius_rad: float,
    spacing_rad: float,
    exclude_tree: cKDTree | None = None,
) -> list[ElevationLine]:
    """`plates._lines_from_resample`'s own algorithm, carrying Hc/Hm (nearest-point lookup)
    instead of a bare scalar elevation -- see that function's docstring for the exclusivity
    logic. `elevation` on the returned lines is a placeholder (zeros); the caller must run
    `lithosphere.sync_plate_elevation` right after to derive the real isostatic value."""
    tree = cKDTree(points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        own_dist, _ = tree.query(world_pts)
        if exclude_tree is None:
            return own_dist < coverage_radius_rad
        other_dist, _ = exclude_tree.query(world_pts)
        return (own_dist < coverage_radius_rad) & (own_dist < other_dist)

    lines: list[ElevationLine] = []
    from .elevation_lines import iter_local_lattice

    for phi, theta_candidates, world_pts in iter_local_lattice(frame, spacing_rad=spacing_rad):
        owned = is_owned(world_pts)
        if not np.any(owned):
            continue
        theta_owned = theta_candidates[owned]
        _, idx = tree.query(world_pts[owned])
        lines.append(
            ElevationLine(
                phi=phi,
                theta=theta_owned,
                elevation=np.zeros(len(theta_owned)),
                crustal_thickness_m=hc[idx],
                mantle_lithosphere_thickness_m=hm[idx],
            )
        )
    return lines


# Hc-noise amplitude that reproduces v1's own CONTINENTAL/OCEANIC_NOISE_AMPLITUDE_M elevation
# swing through the isostasy formula rather than through elevation directly: amplitude_Hc =
# amplitude_elevation / (dz/dHc) at each crust type's own reference thickness. Continental
# dz/dHc = 1 - rho_c/rho_a =~ 0.169 (dry); oceanic uses the water-loaded factor (rho_a >
# rho_w always applies at the oceanic reference depth) =~ 0.156 -- see lithosphere.py's own
# isostatic_elevation for both branches.
_HC_NOISE_AMPLITUDE_CONTINENTAL_M = 2000.0 / (1.0 - lithosphere.RHO_CONTINENTAL_CRUST / lithosphere.RHO_ASTHENOSPHERE)
_HC_NOISE_AMPLITUDE_OCEANIC_M = 900.0 / (
    (1.0 - lithosphere.RHO_OCEANIC_CRUST / lithosphere.RHO_ASTHENOSPHERE) * (lithosphere.RHO_ASTHENOSPHERE / (lithosphere.RHO_ASTHENOSPHERE - lithosphere.RHO_WATER))
)

# Amplitudes for terrain_noise.ContinentalRelief.uplift() (orogenic belts + plateaus).
# Expressed in metres of the elevation swing each contribution should produce, then divided
# by the same 2000 m the CONTINENTAL amplitude above bakes in -- so multiplying the relief
# field (in those units) by `_HC_NOISE_AMPLITUDE_CONTINENTAL_M` lands the contribution at
# the intended elevation through the isostasy formula, and continental crust keeps one
# single Hc-noise amplitude. `_OROGENIC_RELIEF_M` is a belt crest's lift over its own
# `sample()` baseline; a between-ridge basin gets ~none of it, so that is also the depth of
# the intermontane valley below the crests. Sized (with plateaus, and the occasional
# overlap of the two) to land routine peaks near 6 km and the tallest near 7-8 km, under
# MAX_ELEVATION_M = 9000 with only occasional clipping.
_OROGENIC_RELIEF_M = 5200.0
_PLATEAU_BASE_UPLIFT_M = 2800.0
_PLATEAU_INTERNAL_RELIEF_M = 900.0
_OROGENIC_RELIEF_UNITS = _OROGENIC_RELIEF_M / 2000.0
_PLATEAU_UPLIFT_UNITS = _PLATEAU_BASE_UPLIFT_M / 2000.0
_PLATEAU_RELIEF_UNITS = _PLATEAU_INTERNAL_RELIEF_M / 2000.0

# Land gate for the uplift term: it ramps from 0 to full over `sample()` values from
# `land_threshold + _UPLIFT_SEA_MARGIN` to `+ _UPLIFT_SEA_MARGIN + _UPLIFT_SEA_RAMP` (units
# of `sample()`, ~2000 m each). So orogeny/plateaus only touch crust already well above sea
# level -- they can never lift a marine node into land or (being non-negative anyway) drop a
# coastal node into the sea, keeping the land set identical to `sample()` alone.
_UPLIFT_SEA_MARGIN = 0.05
_UPLIFT_SEA_RAMP = 0.35

# RNG tag distinguishing the terrain-noise stream from every other `(seed, plate_id, ...)`
# stream a plate draws (fault noise uses 9001, etc). numpy's SeedSequence only accepts
# integers, so this is an int, not the string "terrain".
_TERRAIN_SEED_TAG = 0x7E44A1


def _continental_sealevel_noise_offset() -> float:
    """`(Hc_at_sealevel - Hc0) / amplitude` for the continental column -- the amount by
    which a node's `relief.sample()` value can sit *below* `land_threshold` and still be
    land, because the reference continental column already floats ~+200 m above sea level.
    Passed to `_land_noise_threshold` so its quantile lands on the true land/sea crossing
    (see that function). Small and negative (~-0.10)."""
    hc0, hm0 = lithosphere.reference_thickness("continental")
    hc_sealevel = float(
        lithosphere.crustal_thickness_for_submerged_elevation(
            np.array([0.0]), np.array([float(hm0)]), lithosphere.RHO_CONTINENTAL_CRUST
        )[0]
    )
    return (hc_sealevel - hc0) / _HC_NOISE_AMPLITUDE_CONTINENTAL_M


# Each plate is seeded with one "primary" site plus this many "extra" sites, and the plate's
# territory is the *union* of its own sites' Voronoi cells rather than a single cell. Merging
# a handful of adjacent cells per plate is what turns the old one-cell-per-plate tiling (every
# plate a convex-ish blob) into lumpier, more continent-like outlines. 0 recovers the exact
# old behaviour. Kept modest on purpose: the more cells a plate fuses, the more concave its
# outline can get, and `PlateWithLines`' per-row outline is only an envelope for a genuinely
# non-convex shape (see `PlateWithLines.outline_world`).
EXTRA_SITES_PER_PLATE = 2


@dataclass
class PlateTiling:
    """A merged-Voronoi partition of the sphere into `num_plates` plates. `site_xyz` is every
    Voronoi site (unit vectors); `site_plate[s]` is which plate owns site `s`'s cell. The
    first `num_plates` sites are the per-plate "primary" sites (`site_plate[:num_plates] ==
    arange(num_plates)`), the rest are extras merged into an existing plate. A world point
    belongs to whichever plate owns its nearest site."""

    site_xyz: np.ndarray
    site_plate: np.ndarray
    num_plates: int

    def primary_site(self, plate_id: int) -> np.ndarray:
        return self.site_xyz[plate_id]


def build_plate_tiling(rng: np.random.Generator, num_plates: int, extra_sites_per_plate: int = EXTRA_SITES_PER_PLATE) -> PlateTiling:
    """Place `num_plates` primary sites plus `num_plates * extra_sites_per_plate` extra sites
    uniformly on the sphere, then hand every extra site to a plate by region-growing: each
    round, the still-unassigned site closest (angularly) to any already-assigned site joins
    that site's plate. Prim-style growth keeps each plate's set of sites a compact cluster,
    so the union of their Voronoi cells stays a single lumpy blob rather than scattering
    disconnected islands across the sphere. Deterministic in `rng`."""
    num_extra = max(0, num_plates * extra_sites_per_plate)
    site_xyz = rng.normal(size=(num_plates + num_extra, 3))
    site_xyz /= np.linalg.norm(site_xyz, axis=-1, keepdims=True)

    site_plate = np.full(len(site_xyz), -1, dtype=int)
    site_plate[:num_plates] = np.arange(num_plates)

    if num_extra > 0:
        angular = np.arccos(np.clip(site_xyz @ site_xyz.T, -1.0, 1.0))
        for _ in range(num_extra):
            unassigned = np.flatnonzero(site_plate < 0)
            assigned = np.flatnonzero(site_plate >= 0)
            block = angular[np.ix_(unassigned, assigned)]
            u, a = np.unravel_index(int(np.argmin(block)), block.shape)
            site_plate[unassigned[u]] = site_plate[assigned[a]]

    return PlateTiling(site_xyz=site_xyz, site_plate=site_plate, num_plates=num_plates)


def generate_plates(
    seed: int,
    num_plates: int | None = None,
    continental_fraction: float | None = None,
    land_fraction: float | None = None,
    node_density: float = 1.0,
    extra_sites_per_plate: int = EXTRA_SITES_PER_PLATE,
) -> list[LithospherePlate]:
    """`plates.generate_plates`'s own seed-placement/Voronoi-tiling algorithm, extended so
    each plate owns the union of several adjacent Voronoi cells (see `build_plate_tiling` and
    `EXTRA_SITES_PER_PLATE`) rather than a single cell -- still deterministic per `seed`, still
    nearest-site-owns-the-node so the tiling has no gaps/overlaps by construction, just with
    lumpier, less convex plate outlines. Only the per-plate line-building step differs from
    v1: each node gets a reference Hc/Hm plus a composite relief field on Hc (see
    `terrain_noise.py` -- a low-frequency `sample()` that decides land/sea exactly as v1's
    single noise did, plus a land-gated non-negative `uplift()` carrying orogenic belts and
    plateaus; `_HC_NOISE_AMPLITUDE_*`/`_OROGENIC_*`/`_PLATEAU_*` above set the amplitudes),
    with `elevation` itself computed once via isostasy at the end."""
    rng = np.random.default_rng(seed)
    if num_plates is None:
        num_plates = int(rng.integers(MIN_AUTO_PLATES, MAX_AUTO_PLATES + 1))

    num_continents: int | None = None
    if continental_fraction is not None:
        continental_fraction = max(0.0, min(continental_fraction, 1.0))
        num_continents = round(continental_fraction * num_plates)
        num_plates = max(num_plates, num_continents + MIN_OCEANIC_PLATES)

    tiling = build_plate_tiling(rng, num_plates, extra_sites_per_plate)
    seed_xyz = tiling.site_xyz

    if num_continents is None:
        crust_types = ["continental" if rng.random() < CONTINENTAL_FRACTION else "oceanic" for _ in range(num_plates)]
    else:
        continental_indices = set(rng.choice(num_plates, size=num_continents, replace=False).tolist())
        crust_types = ["continental" if i in continental_indices else "oceanic" for i in range(num_plates)]

    # Per-site crust type (each site inherits its owning plate's) -- `_land_noise_threshold`
    # and the `is_owned` test below both index by nearest *site*, not nearest plate.
    site_crust_types = [crust_types[tiling.site_plate[s]] for s in range(len(seed_xyz))]

    owner_tree = cKDTree(seed_xyz)
    # Composite relief fields (see terrain_noise.py) -- the last consumers of `rng`, drawn in
    # a fixed order so a given seed reproduces the same terrain. `relief.sample()` stands in
    # for the old single `SphereNoise` (same std, same land/sea decision); `relief.uplift()`
    # adds the orogenic belts and plateaus, land-gated in `hc_at` below.
    relief = terrain_noise.ContinentalRelief(
        rng,
        orogenic_units=_OROGENIC_RELIEF_UNITS,
        plateau_units=_PLATEAU_UPLIFT_UNITS,
        plateau_relief_units=_PLATEAU_RELIEF_UNITS,
    )
    ocean_relief = terrain_noise.OceanicRelief(rng)

    land_threshold = None
    if land_fraction is not None:
        land_fraction = max(0.0, min(land_fraction, 1.0))
        land_threshold = _land_noise_threshold(
            owner_tree, site_crust_types, relief, land_fraction, _continental_sealevel_noise_offset()
        )

    spacing_rad = line_spacing_rad(node_density)
    plates: list[LithospherePlate] = []
    for i in range(num_plates):
        frame = geometry.plate_frame_from_seed(tiling.primary_site(i))
        crust_type = crust_types[i]
        hc0, hm0 = lithosphere.reference_thickness(crust_type)
        hc_amp = _HC_NOISE_AMPLITUDE_CONTINENTAL_M if crust_type == "continental" else _HC_NOISE_AMPLITUDE_OCEANIC_M

        def is_owned(world_pts: np.ndarray, _i: int = i) -> np.ndarray:
            _, nearest_idx = owner_tree.query(world_pts)
            return tiling.site_plate[nearest_idx] == _i

        if crust_type == "continental":
            _lt = 0.0 if land_threshold is None else land_threshold

            def hc_at(world_pts: np.ndarray, _hc0: float = hc0, _amp: float = hc_amp, _lt: float = _lt) -> np.ndarray:
                s = relief.sample(world_pts)
                gate = np.clip((s - _lt - _UPLIFT_SEA_MARGIN) / _UPLIFT_SEA_RAMP, 0.0, 1.0)
                return _hc0 + _amp * (s - _lt) + _amp * gate * relief.uplift(world_pts)
        else:

            def hc_at(world_pts: np.ndarray, _hc0: float = hc0, _amp: float = hc_amp) -> np.ndarray:
                return _hc0 + _amp * ocean_relief.sample(world_pts)

        def elevation_at(world_pts: np.ndarray) -> np.ndarray:
            return np.zeros(len(world_pts))  # placeholder; synced from Hc/Hm below

        lines = build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)
        hc_lines = []
        for line in lines:
            world_pts = line.world_xyz(frame)
            hc = np.clip(hc_at(world_pts), lithosphere.MIN_CRUSTAL_THICKNESS_M, None)
            hm = np.full(len(line), hm0)
            hc_lines.append(line.replace(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm))

        plate = LithospherePlate(plate_id=i, frame=frame, crust_type=crust_type, lines=hc_lines)
        lithosphere.sync_plate_elevation(plate)
        plates.append(plate)

    # Each plate above was seeded purely from its own crust type: submerged continental crust
    # is a uniform bright shelf however far from land it sits, and every continent/ocean
    # boundary is a vertical cliff. Drown the offshore continental interiors and grade the
    # boundary steps into slopes -- see bathymetry.shape_initial_bathymetry.
    bathymetry.shape_initial_bathymetry(plates)
    return plates


def new_plate(plate_id: int, frame: np.ndarray, crust_type: str, spacing_rad: float, seed: int) -> LithospherePlate:
    """A brand-new `LithospherePlate` covering `frame`'s entire local lattice at
    `spacing_rad`, seeded with reference Hc/Hm plus the same composite relief field
    `generate_plates` uses (see `terrain_noise.py`) -- the v2 analogue of
    `plates.generate_plates`' own per-plate initial-line construction. Keyed off
    `(seed, plate_id, _TERRAIN_SEED_TAG)` so the crust stays attached to this plate."""
    hc0, hm0 = lithosphere.reference_thickness(crust_type)
    rng = np.random.default_rng((seed, plate_id, _TERRAIN_SEED_TAG))
    if crust_type == "continental":
        relief = terrain_noise.ContinentalRelief(
            rng,
            orogenic_units=_OROGENIC_RELIEF_UNITS,
            plateau_units=_PLATEAU_UPLIFT_UNITS,
            plateau_relief_units=_PLATEAU_RELIEF_UNITS,
        )
        amp = _HC_NOISE_AMPLITUDE_CONTINENTAL_M

        def hc_at(world_pts: np.ndarray) -> np.ndarray:
            s = relief.sample(world_pts)
            gate = np.clip((s - _UPLIFT_SEA_MARGIN) / _UPLIFT_SEA_RAMP, 0.0, 1.0)
            return hc0 + amp * s + amp * gate * relief.uplift(world_pts)
    else:
        ocean_relief = terrain_noise.OceanicRelief(rng)
        amp = _HC_NOISE_AMPLITUDE_OCEANIC_M

        def hc_at(world_pts: np.ndarray) -> np.ndarray:
            return hc0 + amp * ocean_relief.sample(world_pts)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        return np.ones(len(world_pts), dtype=bool)

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        return np.zeros(len(world_pts))

    lines = build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)
    hc_lines = []
    for line in lines:
        world_pts = line.world_xyz(frame)
        hc = np.clip(hc_at(world_pts), lithosphere.MIN_CRUSTAL_THICKNESS_M, None)
        hm = np.full(len(line), hm0)
        hc_lines.append(line.replace(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm))
    plate = LithospherePlate(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=hc_lines)
    lithosphere.sync_plate_elevation(plate)
    return plate

"""Plate generation and the column physics shared by the tectonic engine.

Defined here: plate generation (`generate_plates`, `new_plate`), building
`PlateWithSparseQuadPatch` plates seeded with reference Hc/Hm and composite relief; and the
per-step column helpers (`boundary_context`, `deform_columns`, `_erupt_melted_nodes`, ...)
that quad_tectonics.py drives.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

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
    CRUST_TYPE_CONTINENTAL,
    CRUST_TYPE_OCEANIC,
    effective_is_continental_from_codes,
    line_spacing_rad,
)
from .noise import SphereNoise
from .plates import (
    CONTINENTAL_FRACTION,
    MIN_AUTO_PLATES,
    MAX_AUTO_PLATES,
    MIN_OCEANIC_PLATES,
    Plate,
    _land_noise_threshold,
)
from . import bathymetry, continental_ledger, cratons, lithosphere, magma_transport, mobile_cover, phase_budget, rheology, shortening, terrain_noise, torque, worldsketch
from .sparse_quad_patch import PlateWithSparseQuadPatch

# A boundary cell may advance into open space whose centre lies within this many node
# spacings of the plate's existing cells (quad_tectonics.py).
EXTEND_THRESHOLD_MULTIPLIER = 1.3
# Hard per-step ceiling on cells a boundary advance or retreat may touch, at the reference
# node density; scaled by sqrt(node_density) at the point of use.
MAX_EXTEND_NODES_PER_STEP = 400

# A rift's mass-conservation share: rather than seed brand-new cells at the full oceanic
# reference column "for free," the volume they would otherwise get is drawn down across this
# many layers of the plate's own existing cells behind the advancing edge *plus* the new layer
# itself, in proportion to how much genuine separation is behind the advance (see
# quad_tectonics._stretch_share). A small constant, not node-count-scaled -- this is about how
# many layers share a stretching event's mass deficit, not a distance or a per-node budget.
K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION = 2

# How many layers of cells a boundary may advance in one step. More than one, because at a
# triple junction, where three independently-oriented plates are all receding from a shared
# point, one layer per plate per step structurally cannot keep pace with three-way divergence
# (confirmed on a real save: the void between three plates widened every step). Capped so a
# single pathological step can't claim an unbounded amount.
MAX_CLAIM_ROWS_PER_STEP = 4

# A plate whose boundary is more deeply/widely overlapping a neighbour right now should
# crumple faster than one barely grazing -- on top of (not instead of) the existing
# distance-decay shape within the belt. `contested_all.mean()` (this plate's own fraction of
# near-boundary-band nodes currently classified contested) is a free, already-computed
# per-step severity signal; at severity 1.0 (the whole band contested -- a deep pile-up) the
# contested-band strength is tripled.
OVERLAP_UPLIFT_SEVERITY_GAIN = 2.0

# Transform (strike-slip) boundary pressure-ridge uplift, applied as a direct elevation delta
# on the transform band in deform() (there is no net crustal shortening at a strike-slip
# contact, so it does not go through Hc/Hm). Gentler than either convergent case -- real
# transform relief is local pressure ridges / transtensional sags, not an orogen. Half v1's
# transform uplift rate (200 m/Myr) since here it is not distance-tapered, only
# fault_influence-gated.
TRANSFORM_UPLIFT_RATE_M_PER_MYR = 100.0

# GitHub issue #189's "still open, separate from this fix" follow-up: transform_uplift is the
# last bare-elevation-delta path deform() writes -- by design (see its own comment: no net
# crustal shortening), unlike the faults.py bug #191 fixed, so it still doesn't go through
# lithosphere.back_elevation_gain. But nothing ever pulls the elevation it adds back down
# either, and unlike every other write to `elevation` here (deform()'s own tectonic Hc deltas,
# erosion.py's isostatic compensation, faults.py/volcanism.py's Hc-backed relief) this is a
# source with literally no decay or consumption mechanism. A real strike-slip pressure ridge
# isn't permanent the way a deep-rooted orogen is either -- lacking a compensating crustal
# root, it's gravitationally unstable and works itself back down over geologic time even
# without dedicated erosion attention. So any *existing* positive debt on a node decays toward
# zero every step (this step's own fresh contribution is added after, so it isn't clawed back
# before it even shows up) -- same age/relax-toward-target idiom
# `rheology.relax_oceanic_mantle_lithosphere` already uses for a different field. Picked
# so a node sitting continuously in a fully-active transform band (worst case, no relief from
# ever leaving the band) settles at a steady-state debt on the order of a few thousand meters
# -- consistent with transform relief being "modest"/"local" by design, not a second orogeny
# -- rather than drifting arbitrarily far before the hard clip.
UNBACKED_RELIEF_DECAY_PER_MYR = 0.05

# A continental plate's *contested* edge is allowed to retreat whether the overriding
# neighbour is oceanic (a passive margin / accretion front: the ocean slab descends under it and
# the buried continental column cedes nothing the model should keep) or continental (a suture
# whose overlapping crust is consumed into the orogen -- its volume is not discarded but thrust
# back onto the plate's own surviving leading edge, see SUTURE_ACCRETION_SPREAD_NODES, so the belt builds real relief in proportion to the overlap it
# actually eats). Left un-retreatable, a contested edge still grows at its *other* (divergent)
# side every step and never back -- the continental ratchet that drives unbounded area creep
# and the slow land-fraction decline (GitHub issue #119, "Node-count creep"). Retreat is gated
# to contested patches of at least this many connected cells, so a single stray cell can't
# nibble a stable coastline or, worse, sever a lobe into a spurious defragmentation plate (the
# failure the naive "retreat every continental contested node" experiment hit -- see GitHub
# issue #119; interior-subduction carving also stays oceanic-only for the same reason).
CONTINENTAL_CONTESTED_RETREAT_MIN_RUN = 3

# Volume-budget growth gate (GitHub issue #119, "Continental ratchet: solution design",
# mechanism 1). A lattice node's physical footprint is constant across the sphere by
# construction (`lithosphere.node_area_m2`), so a plate's total area is just its node count
# times that -- and its implied *mean* crustal thickness is `mean(crustal_thickness_m)`.
# The continental boundary ratchet dilutes this: boundary advance seeds every new margin cell
# at the *oceanic* reference column (`growth_seed_thickness`), so a shear-stretched
# continental plate tiles unbounded drowned passive-margin outward -- node
# count creeps ~+5-6% per 150 My and the plate interior isostatically oceanises into a
# "giant 80%-drowned continent" (GitHub issue #119 items 2 / 5, and issue #120's land-fraction decline).
#
# The gate counts a plate's *genuine* continental nodes -- Hc at least
# `CONTINENTAL_BUDGET_HC_FRACTION` of the continental reference -- and, once the plate's
# total node count exceeds `CONTINENTAL_AREA_BUDGET_MULT` times that count, suppresses all
# areal *growth* for the step (quad_tectonics.py's boundary advance). Retreat, divergent thinning and convergent thickening keep
# running, so an over-budget plate thins / drowns / crumples back toward its crustal volume
# rather than merely freezing. A real craton sits near reference Hc across its whole area,
# nowhere near the cap; the >1 multiplier is the realistic shelf + accreted-terrane
# allowance. Regime-independent and neighbour-independent.
CONTINENTAL_BUDGET_HC_FRACTION = 0.6
CONTINENTAL_AREA_BUDGET_MULT = 1.8

# When a continental *suture* edge retreats (a continental neighbour overrides it -- not an
# oceanic one, where the buried column genuinely subducts and is lost), the removed column's
# crustal volume is conserved: it is thrust back onto the plate's own surviving leading-edge
# nodes, spread over this many of them (an imbricate thrust wedge), the attached mantle
# lithosphere thickening in proportion. This is the mass-honest replacement for the retired
# `rheology.CONTINENTAL_COLLISION_SHORTENING_BOOST` fudge (a flat 2.5x `fault_factor`
# multiplier at continent-continent contested nodes, unrelated to how much overlap was
# actually consumed). Node area is constant per node (`lithosphere.node_area_m2`), so
# conserving volume is just moving the summed Hc of the dropped nodes onto the survivors, and
# isostasy lifts the thickened belt. On cells, "nodes" are hops behind the retreating edge
# (quad_tectonics.py).
SUTURE_ACCRETION_SPREAD_NODES = 3

# Hard ceiling on a node's Hc after suture accretion. A suture that never heals (the
# neighbour keeps overriding) would otherwise pile every consumed column onto the same few
# retreating-edge nodes indefinitely -- Hc ran to ~190 km and climbing on a 30-My test run.
# Real orogenic crust does not stack past ~2x reference: the excess root is removed by
# lower-crustal / mantle-lithosphere delamination (and the surface by erosion). Accreted
# mass over this ceiling is dropped (delaminated), so accretion is mass-conserving only up
# to the cap -- which a normal collision, healing over ~1-2 My, never reaches. Same ceiling
# `lithosphere.MAX_CRUSTAL_THICKNESS_M` uses for ordinary convergent thickening (issue #161)
# -- both paths agree on where continental crust actually maxes out.
SUTURE_ACCRETION_MAX_HC_M = lithosphere.MAX_CRUSTAL_THICKNESS_M

# GitHub issue #177: a continental passive margin overridden by an *oceanic* neighbour
# subducts its retreating column outright, with no rate limit of its own -- unlike the
# compensating gain on the same convergent boundary, `rheology.apply_arc_magmatic_thickening`,
# which is rate-capped at `rheology.ARC_MAGMATIC_CONVERGENCE_CAP`. The issue's own ledger
# instrumentation confirmed this uncapped-loss/capped-gain pairing is the actual driver of
# `avg_rotation_rate`'s super-linear land loss: faster rotation -> deeper/more overlap per
# step -> more uncapped subduction, while the capped creation path can't scale to match.
#
# Rather than invent a second hardcoded rate for the loss side (which would just be two
# constants racing each other, free to drift apart again as either one gets retuned), this
# plate's oceanic-override retreat is capped, this same step, by however much Hc this exact
# step's real arc-magmatic thickening is actually adding across the whole plate --
# `boundary_context` computes that total once (mirroring the calculation
# `apply_arc_magmatic_thickening` performs per node for real) into a shared budget, and
# quad_tectonics.py's retreat spends it down. Both sides are area-weighted, so summed Hc is
# directly comparable as "volume" on both sides of this cap. A continent-continent suture
# retreat is untouched by this (it already fully conserves its own volume, so it isn't the
# uncapped channel this targets), and so is an oceanic self-plate's own ordinary subduction
# (not a passive margin at all -- that loss is expected, not a bug).
#
# 1.0 is a genuinely symmetric cap (this issue's own proposal: match the two channels' rates
# directly), kept as a named multiplier rather than folded into the budget calculation itself
# so it is the one place to loosen this if a strict 1:1 turns out to starve ordinary retreat
# unrelated to the runaway this exists to fix.
OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER = 1.0


# Active-margin (Cordilleran) accretion. When a continental plate's *leading* edge grows
# into space a subducting oceanic neighbour is vacating (slab rollback / trench retreat),
# the new ground is juvenile arc + accreted-terrane crust, not abyssal sea floor -- so it is
# seeded at this intermediate column (Hc ~0.8x continental reference) rather than
# `growth_seed_thickness`'s drowned oceanic one. This is the deliberately *restricted*
# reverse of the land-area runaway that `growth_seed_thickness` documents: the runaway was
# seeding +200 m dry land on *every* growth event, including growth into open ocean far from
# any margin; seeding a thicker column *only* where the growing edge abuts a genuinely
# converging oceanic slab -- and still under the `CONTINENTAL_AREA_BUDGET_MULT` volume gate
# -- is arc accretion, the dominant land-loss driver's actual physical counterweight (see
# GitHub issue #120, "Land fraction slowly declines"). The seed lands as shallow forearc/shelf
# (~ -450 m) and builds to land as convergence continues via
# `rheology.apply_arc_magmatic_thickening` + ordinary convergent shortening.
ARC_MARGIN_SEED_HC_M = 28_000.0
ARC_MARGIN_SEED_HM_M = 55_000.0

# How many hops in from a growing edge are scanned for an active-margin signal -- a node
# contested by an oceanic neighbour, or one still carrying a subduction-arc provenance stamp
# from a recent step -- when deciding whether that edge's growth seeds arc crust or ocean
# floor. Small: the signal only has to survive the one step between the ocean's edge
# retreating and this plate's edge growing into the gap.
ARC_MARGIN_END_SCAN_NODES = 4


def growth_seed_thickness() -> tuple[float, float]:
    """(Hc, Hm) a plate seeds *brand-new areal* nodes with -- when its boundary advances into
    open water (quad_tectonics.py).

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


def _erupt_melted_nodes(
    world: "World",  # noqa: F821
    plate_id: int,
    rng_index: int,
    hc: np.ndarray,
    hm: np.ndarray,
    crust_type_code: np.ndarray,
    is_volcano: np.ndarray,
    volcano_remaining: np.ndarray,
    melting: np.ndarray,
    prior_elevation: np.ndarray,
) -> None:
    """In-place: any node flagged `melting` (its Hc just crossed below
    `RIFT_CRITICAL_THICKNESS_M` -- from divergent thinning or from rift stretching alike,
    thinning is thinning) erupts fresh crust in
    place, typed by whether it was still standing above sea level (`prior_elevation > 0`,
    this node's elevation *before* today's melting event) at the moment it melted through --
    the bimodal continental-rift-volcanism vs. ordinary mid-ocean-ridge distinction `deform()`'s
    own module-level docstring describes. Mutates hc/hm/crust_type_code/is_volcano/
    volcano_remaining in place; the caller still owns recomputing elevation from the resulting
    hc/hm afterward (isostasy needs each melted node's own new crust_type_code-driven density,
    which can now vary node to node)."""
    if not np.any(melting):
        return
    from .elevation_lines import VOLCANO_ACTIVE_MAX_YEARS, VOLCANO_ACTIVE_MIN_YEARS

    melt_land = melting & (prior_elevation > 0.0)
    melt_ocean = melting & ~melt_land
    hc[melt_land] = lithosphere.REFERENCE_HC_CONTINENTAL_M
    hm[melt_land] = lithosphere.REFERENCE_HM_CONTINENTAL_M
    crust_type_code[melt_land] = CRUST_TYPE_CONTINENTAL
    hc[melt_ocean] = lithosphere.REFERENCE_HC_OCEANIC_M
    hm[melt_ocean] = lithosphere.YOUNG_RIDGE_HM_M
    crust_type_code[melt_ocean] = CRUST_TYPE_OCEANIC
    is_volcano[melting] = True
    rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate_id, rng_index))
    volcano_remaining[melting] = rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=int(melting.sum()))


def _ignite_early_rift_volcanoes(
    world: "World",  # noqa: F821
    plate_id: int,
    rng_index: int,
    is_volcano: np.ndarray,
    volcano_remaining: np.ndarray,
    ignite: np.ndarray,
) -> None:
    """In-place: any node flagged `ignite` (crossed below `rheology.RIFT_VOLCANISM_ONSET_HC_M`
    this step without fully melting through -- see `rheology.apply_rift_magmatic_thickening`)
    starts its own point-volcano eruption lifecycle early, the same one ordinary
    decompression-melt volcanoes (`_erupt_melted_nodes`) and arc volcanoes already use
    (volcanism.py) -- a real continental rift is volcanically active well before it actually
    ruptures (GitHub issue #120, "Land fraction slowly declines", "over-stretched interiors").
    Mutates is_volcano/volcano_remaining in place; unlike `_erupt_melted_nodes` this never
    touches hc/hm/crust_type_code -- the node is still ordinary, still-thinning continental
    crust, just one that's now also erupting onto its own surface."""
    if not np.any(ignite):
        return
    from .elevation_lines import VOLCANO_ACTIVE_MAX_YEARS, VOLCANO_ACTIVE_MIN_YEARS

    is_volcano[ignite] = True
    # Distinct rng key (trailing 1) from _erupt_melted_nodes' own draw above so the two
    # event types don't share a random stream at the same (seed, step, plate, stream).
    rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate_id, rng_index, 1))
    volcano_remaining[ignite] = rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=int(ignite.sum()))


@dataclass
class BoundaryContext:
    """One plate's per-step boundary classification and the per-plate knobs derived from it --
    everything `deform()` decides before touching any node (quad_tectonics.py). Every array is per node, in the plate's own
    node order at the start of deform()."""

    own_points: np.ndarray
    spacing_rad: float
    reach_rad: float
    neighbours: list
    inputs: torque.BoundaryForceInputs
    convergent: np.ndarray
    divergent: np.ndarray
    transform: np.ndarray
    contested: np.ndarray
    shrinkable: np.ndarray
    accrete: np.ndarray
    closing_rate: np.ndarray
    arc_band: np.ndarray
    arc_intensity: np.ndarray
    fault_influence: np.ndarray
    suppress_growth: bool
    # Single-element and mutable: the retreat spends it down in place (quad_tectonics._retreat).
    oceanic_override_retreat_budget_hc: np.ndarray
    suture_hm_subduct: np.ndarray
    orogen_amount: float
    orogen_contested_strength: float
    fault_noise: SphereNoise | None


def continental_arc_band(plate: Plate, inputs: torque.BoundaryForceInputs, closing_rate: np.ndarray) -> np.ndarray:
    """A continental plate's own nodes within reach of a converging *oceanic* neighbour plate
    -- the arc band `boundary_context` thickens and grows arc crust on. A model cue, not an
    observed slab: it follows continental ownership, an oceanic neighbour and positive
    closing, with no independent test of which way the slab dips."""
    if plate.crust_type != "continental":
        return np.zeros(len(inputs.own_points), dtype=bool)
    return inputs.neighbor_is_oceanic & np.isfinite(inputs.dist_to_neighbor) & (closing_rate > rheology.ARC_MIN_CONVERGENCE_M_PER_S)


def boundary_context(
    world: "World",  # noqa: F821
    plate: Plate,
    other_plates: list,
    years: float,
    continental_retreat_runs,
    node_weight: np.ndarray | float = 1.0,
) -> BoundaryContext:
    """Classify `plate`'s boundary for this step's deform() and derive every per-plate knob
    the column and topology updates read -- see the comments below for each piece's
    rationale.

    `continental_retreat_runs(contested)` keeps the contested nodes that form a genuine
    multi-node stretch (connected components of cells, see
    `quad_tectonics.components_of_at_least`), and `node_weight` is each node's footprint in
    nominal-node units (`cell area / lithosphere.node_area_m2`; a scalar treats every node as
    one nominal footprint) so the continental area budget and the issue #177 retreat budget
    count area, not nodes."""
    own_points, _ = plate.all_points_and_elevation()
    hc_all = plate.collect("crustal_thickness_m")

    # Volume-budget growth gate -- see CONTINENTAL_AREA_BUDGET_MULT. Continental crust
    # only: oceanic footprint is already bounded by subduction. Over budget -> this step
    # grows no new areal crust, but still retreats / thins / thickens toward the budget.
    suppress_growth = False
    if plate.crust_type == "continental":
        genuine = hc_all >= CONTINENTAL_BUDGET_HC_FRACTION * lithosphere.REFERENCE_HC_CONTINENTAL_M
        if np.ndim(node_weight) == 0:
            suppress_growth = len(own_points) > CONTINENTAL_AREA_BUDGET_MULT * int(np.count_nonzero(genuine))
        else:
            suppress_growth = float(np.sum(node_weight)) > CONTINENTAL_AREA_BUDGET_MULT * float(np.sum(node_weight[genuine]))

    spacing_rad = line_spacing_rad(world.node_density)
    reach_rad = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing_rad

    # The step's shared boundary searches, when the collision-polarity prepass ran this step
    # (issue #318): a neighbour that hasn't deformed yet is answered from the prepass's own
    # search -- see torque.BoundarySearchCache.
    cache = getattr(world, "boundary_search_cache", None)
    if cache is None:
        neighbours = plate.get_neighbours(other_plates, threshold_rad=reach_rad)
    else:
        neighbours = cache.neighbours(plate, other_plates, reach_rad)
    inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing_rad, reach_rad, cache)
    # Motion-based: `convergent` is the whole converging band (not just the nodes that
    # already overlap a neighbour polygon), so a boundary builds an orogen before any overlap
    # accumulates; `contested` (the geometric overlap subset, folded into `convergent`) still
    # gates node deletion / continental retreat.
    convergent, divergent, transform, contested = torque.classify_boundary_nodes(plate, neighbours, inputs, reach_rad, cache)

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

        fault_influence = faults.fault_influence(world, plate, own_points)
    else:
        fault_influence = np.ones(len(own_points))

    closing_rate = rheology.normal_closing_rate_m_per_s(plate.omega, inputs.neighbor_omega, own_points, inputs.direction_to_neighbor)

    # What may retreat this step. Oceanic crust: any contested node subducts. Continental
    # crust: any contested boundary node in a run of >= CONTINENTAL_CONTESTED_RETREAT_MIN_RUN
    # contested nodes -- whether the overriding neighbour is oceanic (passive margin) or
    # continental (a suture whose overlap is consumed into the orogen, the retreated
    # column's volume thrust onto the plate's own leading edge -- see
    # SUTURE_ACCRETION_SPREAD_NODES). Envelope fuzz (a lone contested node) still can't
    # nibble a stable margin. See CONTINENTAL_CONTESTED_RETREAT_MIN_RUN for the ratchet /
    # frozen-overlap this breaks.
    if plate.crust_type != "continental":
        shrinkable = contested
    else:
        shrinkable = continental_retreat_runs(contested)

    # Apply frozen continental-collision fronts independently. The prepass stores each
    # side's retreatability in its original node order, before any plate can deform.
    suture_hm_subduct = np.zeros(len(own_points), dtype=bool)
    frame = getattr(world, "collision_polarity_frame", None)
    masks = getattr(frame, "masks", {}).get(plate.plate_id) if frame is not None else None
    if (
        masks is not None
        and len(masks.front_id) == len(own_points)
        and np.array_equal(masks.node_keys, plate.cell_keys)
        and len(masks.retreatable) == len(own_points)
        and len(masks.continental) == len(own_points)
    ):
        active_fronts: list[int] = []
        plate_by_id = {p.plate_id: p for p in [plate, *other_plates]}
        for front_id, pair in frame.polarity.items():
            if plate.plate_id not in pair:
                continue
            counterpart = pair[1] if pair[0] == plate.plate_id else pair[0]
            other = plate_by_id.get(counterpart)
            other_masks = getattr(frame, "masks", {}).get(counterpart)
            if (
                other is None
                or other_masks is None
                or len(other_masks.front_id) != len(other_masks.retreatable)
                or len(other_masks.front_id) != len(other_masks.continental)
                or len(other_masks.front_id) != len(other_masks.node_keys)
            ):
                continue
            own_front = masks.front_id == front_id
            other_front = other_masks.front_id == front_id
            if not np.any(own_front) or not np.any(other_front):
                continue
            own_cont = masks.continental
            other_cont = other_masks.continental
            if not np.any(own_front & own_cont) or not np.any(other_front & other_cont):
                continue
            lower_id = pair[0]
            lower_masks = masks if lower_id == plate.plate_id else other_masks
            lower_front = lower_masks.front_id == front_id
            lower_cont = own_cont if lower_id == plate.plate_id else other_cont
            lower_can_retreat = bool(np.any(lower_masks.retreatable & lower_front & lower_cont))
            # If every lower-plate candidate is too small or still cratonic, preserve the
            # ordinary bilateral retreat for this front so overlap cannot become permanent.
            if lower_can_retreat:
                active_fronts.append(front_id)
        if active_fronts:
            collision_nodes = np.isin(masks.front_id, active_fronts)
            active_lower = masks.lower & collision_nodes
            own_cont = masks.continental
            selected = masks.retreatable & active_lower & own_cont
            shrinkable = (shrinkable & ~collision_nodes) | selected
            suture_hm_subduct = selected & own_cont

    # Continental suture retreat conserves the consumed column's volume by accreting it
    # onto this plate's own leading edge; a retreat where the overriding neighbour is
    # *oceanic* does not -- that column subducts and is lost. Oceanic self-plates never
    # accrete.
    if plate.crust_type == "continental":
        accrete = shrinkable & ~inputs.neighbor_is_oceanic
    else:
        accrete = np.zeros_like(shrinkable)

    # Continental arc band: this plate's own nodes within `reach_rad` of a *converging
    # oceanic* neighbour -- the volcanic arc + accreted forearc / underplated wedge sits
    # inboard of the trench, a swath (~500 km at default density), not just the contact
    # (which is only a few tens of nodes -- far too narrow to counter the land
    # decline). `arc_intensity` fades from 1 at the contact to ~0.3 at the band edge.
    # Feeds both the magmatic Hc thickening and the arc-crust growth seed. See
    # ARC_MARGIN_SEED_HC_M.
    arc_band = continental_arc_band(plate, inputs, closing_rate)
    arc_intensity = np.zeros(len(own_points))
    if plate.crust_type == "continental":
        arc_intensity = np.where(arc_band, np.clip(1.0 - 0.7 * (inputs.dist_to_neighbor / reach_rad), 0.3, 1.0), 0.0)

    years_myr = years / 1_000_000.0

    # GitHub issue #177 direction 1 -- see OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER's own
    # comment for the full rationale. This step's real arc-magmatic creation across the
    # whole plate, computed once here (mirrors the per-node calculation in
    # `deform_columns`) into a shared budget the oceanic-override retreats spend down in place.
    oceanic_override_retreat_budget_hc = np.zeros(1)
    if plate.crust_type == "continental" and np.any(arc_band):
        hm_all = plate.collect("mantle_lithosphere_thickness_m")
        grown_hc, _ = rheology.apply_arc_magmatic_thickening(
            hc_all[arc_band], hm_all[arc_band], closing_rate[arc_band], years_myr, arc_intensity[arc_band],
        )
        grown = grown_hc - hc_all[arc_band]
        if np.ndim(node_weight) != 0:
            grown = grown * node_weight[arc_band]
        oceanic_override_retreat_budget_hc[0] = OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER * float(np.sum(grown))

    # Collision-uplift tuning knobs (the "Controls" window, 1.0 == untuned -- see World).
    # `orogen_amount` scales the plastic thickening rate at contested nodes; `orogen_reach`
    # below 1 narrows it in proportion. At 1.0 both leave apply_convergent_deformation's
    # contested-band strength at exactly 1.0. `orogen_reach` also sets how far the shortening
    # cascade carries a collision into the plate (shortening.py's `reach_scale`, issue #314).
    orogen_amount = world.collision_uplift_multiplier
    orogen_reach = world.collision_uplift_reach_multiplier
    # How much of this plate's own active margin is currently jammed in overlap right now
    # -- normalized against the near-boundary band, not the whole plate (a huge plate's
    # boundary is a small fraction of its own node count, which would dilute this to
    # near-zero for exactly the large-plate case that matters). A deeper/wider overlap
    # should crumple faster than a light graze, on top of the existing distance-decay
    # shape within the belt -- see OVERLAP_UPLIFT_SEVERITY_GAIN.
    band = inputs.dist_to_neighbor <= reach_rad
    overlap_severity = float(np.count_nonzero(contested)) / max(1, int(np.count_nonzero(band)))
    orogen_contested_strength = orogen_amount * min(orogen_reach, 1.0) * (1.0 + OVERLAP_UPLIFT_SEVERITY_GAIN * overlap_severity)

    fault_noise = (
        SphereNoise(np.random.default_rng((world.seed, plate.plate_id, 9001)), octaves=3, base_freq=9.0)
        if plate.crust_type == "continental"
        else None
    )

    return BoundaryContext(
        own_points=own_points,
        spacing_rad=spacing_rad,
        reach_rad=reach_rad,
        neighbours=neighbours,
        inputs=inputs,
        convergent=convergent,
        divergent=divergent,
        transform=transform,
        contested=contested,
        shrinkable=shrinkable,
        accrete=accrete,
        closing_rate=closing_rate,
        arc_band=arc_band,
        arc_intensity=arc_intensity,
        fault_influence=fault_influence,
        suppress_growth=suppress_growth,
        oceanic_override_retreat_budget_hc=oceanic_override_retreat_budget_hc,
        suture_hm_subduct=suture_hm_subduct,
        orogen_amount=orogen_amount,
        orogen_contested_strength=orogen_contested_strength,
        fault_noise=fault_noise,
    )


# The per-node fields `deform_columns` reads and returns.
COLUMN_FIELDS = (
    "elevation",
    "crustal_thickness_m",
    "mantle_lithosphere_thickness_m",
    "divergent_age_myr",
    "is_volcano",
    "volcano_active_years_remaining",
    "elev_change_reason",
    "crust_type_code",
    "craton_crust_m",
    "continental_material_m",
    "mobile_cover_m",
    "mobile_cover_continental_m",
)


def deform_columns(
    world: "World",  # noqa: F821
    plate: Plate,
    ctx: BoundaryContext,
    fields: dict[str, np.ndarray],
    local_xyz,
    node_area_m2: np.ndarray | float,
    rng_index: int,
    years: float,
    ceiling_overflow: np.ndarray | None = None,
    strained: dict[str, np.ndarray] | None = None,
    accommodate: Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray] | None = None,
) -> dict[str, np.ndarray]:
    """This step's in-place lithospheric column update for every node of `plate`, against
    `ctx`'s per-plate arrays -- convergent thickening, arc magmatism, divergent thinning and
    decompression melting, oceanic cooling, transform pressure ridges, and provenance
    stamping. No topology change: returns new values for every `COLUMN_FIELDS` name, same
    length as `fields`' own arrays.

    Called by the quad-surface engine once per plate. The caller supplies `local_xyz()`
    (lazily, for the fault-noise texture), `node_area_m2` (per node, or a scalar for a
    uniform footprint) and `rng_index` (the eruption rng's per-plate stream key).

    `ceiling_overflow`, when given, receives the convergent band's per-node Hc past
    `MAX_CRUSTAL_THICKNESS_M` for the caller to place (issue #290); without it that overflow
    delaminates. `strained`, when
    given, receives the columns with only this step's tectonic strain applied --
    `crustal_thickness_m` and `mantle_lithosphere_thickness_m` after convergent shortening
    and divergent thinning but before any arc or rift magma -- and `melted`, the nodes that
    melted through and were reset, for the caller's thermal bookkeeping.

    `accommodate`, when given, takes the convergent band's shortening off the band (issue
    #314): it is called with each node's demanded strain (`rheology.convergent_strain`, after
    the magma-export skim) and its pre-shortening Hc/Hm, and returns the strain each node
    actually takes up -- spread into the plate's interior and kept under the caps
    (shortening.py). The band then never overflows its ceiling."""
    convergent = ctx.convergent
    divergent = ctx.divergent
    transform = ctx.transform
    closing_rate = ctx.closing_rate
    neighbor_oceanic = ctx.inputs.neighbor_is_oceanic
    arc_band = ctx.arc_band
    arc_intensity = ctx.arc_intensity
    fault_influence = ctx.fault_influence  # all-ones except in "fault" mode
    elevation = fields["elevation"]
    n = len(elevation)
    years_myr = years / 1_000_000.0

    hc = fields["crustal_thickness_m"].copy()
    hm = fields["mantle_lithosphere_thickness_m"].copy()
    continental_material = fields["continental_material_m"].copy()
    # Mobile cover (mobile_cover.py) is the top of the column: stretching thins it with the
    # column and a column that melts through loses it; shortening and underplating leave it.
    cover = fields["mobile_cover_m"].copy()
    cover_material = fields["mobile_cover_continental_m"].copy()
    # GitHub issue #216 Hc/Hm budget checkpoints -- see phase_budget.py. `codes0` is
    # these nodes' crust_type_code, unchanged until the decompression-melting checkpoint
    # below, so every intermediate checkpoint below reuses it for both before/after.
    codes0 = fields["crust_type_code"]
    checkpoint_hc, checkpoint_hm = (hc.copy(), hm.copy()) if world.debug_diagnostics else (None, None)
    # Every phase below changes columns in place on the same nodes, so before and after share
    # each node's area (issue #257: per-cell on quad plates).
    budget_area_m2 = np.broadcast_to(node_area_m2, hc.shape)
    budget_craton_m = fields["craton_crust_m"]
    budget_areas = {
        "area_before_m2": budget_area_m2,
        "area_after_m2": budget_area_m2,
        "craton_before_m": budget_craton_m,
        "craton_after_m": budget_craton_m,
    }
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
    rho_c = lithosphere.crust_density(plate.crust_type)
    elevation_before = lithosphere.isostatic_elevation(hc, hm, rho_c)

    # The band that plastically thickens: the whole converging band at
    # `orogen_contested_strength`. `orogen_strength` is the per-node multiplier handed to
    # apply_convergent_deformation; > 0 exactly on the nodes that thicken.
    # `apply_convergent_deformation` still gates on each node's own closing rate
    # (below yield / not actually closing -> zero strain), so a node that is
    # `convergent` only via the `contested` deep-overlap fold and is no longer
    # actively closing simply thickens at zero.
    orogen_strength = np.where(convergent, ctx.orogen_contested_strength, 0.0)
    # "fault" mode: concentrate the shortening onto fault traces (no-op / all-ones
    # otherwise). `strength` scales apply_convergent_deformation's thickening rate.
    orogen_strength = orogen_strength * fault_influence
    thicken = orogen_strength > 0.0
    strained_hc, strained_hm = hc.copy(), hm.copy()
    if np.any(thicken):
        fault_factor = (
            np.where(
                ctx.fault_noise.sample(local_xyz()) < -0.15,
                rheology.REVERSE_FAULT_VALLEY_UPLIFT_FACTOR,
                1.0,
            )
            if ctx.fault_noise is not None
            else np.ones(n)
        )
        # The overlapping crust a continent-continent suture retreats over is not
        # lost here via a `fault_factor` boost -- its actual volume is conserved and
        # thrust onto the leading edge by the retreat step (see
        # SUTURE_ACCRETION_SPREAD_NODES). This path is just the ordinary
        # yield-limited plastic thickening.

        # Lateral magma export (GitHub issue #205, follow-up to #120's "Land fraction
        # slowly declines"): divert a fraction of the core convergent band's own strain
        # increment to a mobile magma parcel instead of thickening the node in place --
        # see rheology.magma_export_strength_and_volume's own docstring for the full
        # mechanism/reasoning (and on why near-ceiling nodes are exempt).
        core_idx = np.flatnonzero(thicken)[convergent[thicken]]
        used_strength = orogen_strength[thicken]
        if len(core_idx) > 0:
            reduced_strength, export_hc = rheology.magma_export_strength_and_volume(
                hc[core_idx], closing_rate[core_idx], years_myr, fault_factor[core_idx], orogen_strength[core_idx],
            )
            used_strength = used_strength.copy()
            used_strength[convergent[thicken]] = reduced_strength
            exporting = export_hc > 0.0
            if np.any(exporting):
                node_idx = core_idx[exporting]
                area = node_area_m2 if np.ndim(node_area_m2) == 0 else node_area_m2[node_idx]
                volumes_m3 = export_hc[exporting] * area
                origins = ctx.own_points[node_idx]
                world.pending_magma_parcels.extend(
                    magma_transport.MagmaParcel(
                        origin_xyz=origins[i], volume_m3=float(volumes_m3[i]), step_generated=world.steps_taken
                    )
                    for i in range(len(node_idx))
                )

        if accommodate is not None:
            demanded = np.zeros(n)
            demanded[thicken] = rheology.convergent_strain(closing_rate[thicken], years_myr, fault_factor[thicken], used_strength)
            hc, hm = shortening.apply_strain(hc, hm, accommodate(demanded, hc, hm))
            strained_hc, strained_hm = hc.copy(), hm.copy()
            overflow_hc = np.zeros(int(np.count_nonzero(thicken)))
        else:
            new_hc, new_hm, overflow_hc = rheology.apply_convergent_deformation(
                hc[thicken], hm[thicken], closing_rate[thicken], years_myr,
                fault_factor[thicken], strength=used_strength,
            )
            hc[thicken] = new_hc
            hm[thicken] = new_hm
            strained_hc[thicken] = new_hc
            strained_hm[thicken] = new_hm

        # Hc that hit MAX_CRUSTAL_THICKNESS_M this step doesn't just vanish (issue #161):
        # the caller places it (quad_tectonics._place_ceiling_overflow, issue #290).
        if ceiling_overflow is not None:
            ceiling_overflow[np.flatnonzero(thicken)[convergent[thicken]]] = overflow_hc[convergent[thicken]]

    if world.debug_diagnostics:
        phase_budget.record(world, plate, "convergent_deformation", checkpoint_hc, checkpoint_hm, codes0, hc, hm, codes0, **budget_areas)
        checkpoint_hc, checkpoint_hm = hc.copy(), hm.copy()

    # Continental arc magmatism: an oceanic slab subducting under this margin fluxes
    # the mantle wedge and underplates juvenile crust across the whole arc band --
    # extra Hc (added from the mantle, not conserved), the crust-building half of
    # "subduction under a continent makes more continent" (GitHub issue #120, "Land fraction
    # slowly declines"). Separate from the contested shortening above: the band is far
    # wider than the contact line. Bounded long-term by the CONTINENTAL_AREA_BUDGET_MULT
    # volume gate.
    hc_before_arc = hc.copy()
    if np.any(arc_band):
        hc[arc_band], hm[arc_band] = rheology.apply_arc_magmatic_thickening(
            hc[arc_band], hm[arc_band], closing_rate[arc_band], years_myr, arc_intensity[arc_band]
        )
        continental_arc = effective_is_continental_from_codes(codes0, plate.crust_type == "continental")
        juvenile = np.where(
            arc_band & continental_arc, np.maximum(hc - hc_before_arc, 0.0), 0.0
        )
        continental_material += juvenile
        continental_ledger.record(
            world, "juvenile_additions_m3", float(np.dot(juvenile, budget_area_m2))
        )

    if world.debug_diagnostics:
        phase_budget.record(world, plate, "arc_magmatism", checkpoint_hc, checkpoint_hm, codes0, hc, hm, codes0, **budget_areas)
        checkpoint_hc, checkpoint_hm = hc.copy(), hm.copy()

    prior_hc = hc.copy()
    melting = np.zeros(n, dtype=bool)
    newly_below_rift_onset = np.zeros(n, dtype=bool)
    craton = fields["craton_crust_m"]
    if np.any(divergent):
        material_before_stretch = continental_material.copy()
        new_hc, new_hm, melt = rheology.apply_divergent_deformation(hc[divergent], hm[divergent], closing_rate[divergent], years_myr)
        # "fault" mode: scale the thinning delta by fault proximity (all-ones
        # otherwise). Melt (decompression volcanism) still fires on the geometric
        # rift threshold -- it's a discrete event, not a rate.
        infl = fault_influence[divergent]
        # A craton resists the thinning itself (cratons.CRATON_RIFT_RESISTANCE), so it only
        # melts through once its own, slower thinning crosses the rift threshold.
        resisted = cratons.strength(craton[divergent]) * cratons.CRATON_RIFT_RESISTANCE
        old_hc = hc[divergent]
        hc[divergent] = old_hc + infl * (1.0 - resisted) * (new_hc - old_hc)
        hm[divergent] = hm[divergent] + infl * (1.0 - resisted) * (new_hm - hm[divergent])
        # Stretching spreads the same material over the enlarged footprint.  Keep its
        # thickness fraction locked to Hc before adding any genuinely juvenile rift magma.
        continental_material[divergent] *= np.divide(
            hc[divergent], prior_hc[divergent],
            out=np.zeros_like(hc[divergent]), where=prior_hc[divergent] > 0.0,
        )
        stretch_loss = np.maximum(material_before_stretch - continental_material, 0.0)
        continental_ledger.record(
            world, "rift_thinned_m3", float(np.dot(stretch_loss, budget_area_m2))
        )
        stretch_ratio = np.clip(
            np.divide(hc[divergent], prior_hc[divergent], out=np.zeros_like(hc[divergent]), where=prior_hc[divergent] > 0.0),
            0.0,
            1.0,
        )
        cover_before_stretch = cover.copy()
        cover[divergent] *= stretch_ratio
        cover_material[divergent] *= stretch_ratio
        mobile_cover.record(world, "rift_thinned_m3", float(np.dot(cover_before_stretch - cover, budget_area_m2)))
        melting[divergent] = np.where(
            resisted > 0.0,
            melt & (hc[divergent] < rheology.RIFT_CRITICAL_THICKNESS_M),
            melt,
        )
        strained_hc[divergent] = hc[divergent]
        strained_hm[divergent] = hm[divergent]

        # Rift magmatic underplating (see rheology.apply_rift_magmatic_thickening): a
        # partial Hc offset for nodes that thinned past RIFT_VOLCANISM_ONSET_HC_M but
        # didn't melt all the way through this step -- nodes that did melt already got
        # the full reference-column reset below and don't need this on top of it.
        magmatic_band = divergent & ~melting
        if np.any(magmatic_band):
            hc_before_rift_magmatism = hc.copy()
            new_hc_mag, new_hm_mag = rheology.apply_rift_magmatic_thickening(
                hc[magmatic_band], hm[magmatic_band], closing_rate[magmatic_band], years_myr
            )
            hc[magmatic_band] = new_hc_mag
            hm[magmatic_band] = new_hm_mag
            continental_now = effective_is_continental_from_codes(codes0, plate.crust_type == "continental")
            juvenile = np.where(
                magmatic_band & continental_now,
                np.maximum(hc - hc_before_rift_magmatism, 0.0),
                0.0,
            )
            continental_material += juvenile
            continental_ledger.record(
                world, "juvenile_additions_m3", float(np.dot(juvenile, budget_area_m2))
            )
            newly_below_rift_onset = (
                magmatic_band & (prior_hc >= rheology.RIFT_VOLCANISM_ONSET_HC_M) & (hc < rheology.RIFT_VOLCANISM_ONSET_HC_M)
            )

    if world.debug_diagnostics:
        phase_budget.record(world, plate, "divergent_deformation", checkpoint_hc, checkpoint_hm, codes0, hc, hm, codes0, **budget_areas)
        checkpoint_hc, checkpoint_hm = hc.copy(), hm.copy()

    prior_age = fields["divergent_age_myr"]
    new_age = np.where(divergent, prior_age + years_myr, 0.0)
    if plate.crust_type == "oceanic":
        oceanic_now = ~effective_is_continental_from_codes(codes0, False)
        hm = rheology.relax_oceanic_mantle_lithosphere(hm, oceanic_now, years_myr)
        if world.debug_diagnostics:
            phase_budget.record(world, plate, "oceanic_cooling_relaxation", checkpoint_hc, checkpoint_hm, codes0, hc, hm, codes0, **budget_areas)
            checkpoint_hc, checkpoint_hm = hc.copy(), hm.copy()

    is_volcano = fields["is_volcano"].copy()
    volcano_remaining = fields["volcano_active_years_remaining"].copy()
    crust_type_code = fields["crust_type_code"].copy()
    # Decompression melting (spec 2.3): a rift that just thinned past the critical
    # threshold erupts fresh crust in place -- same one-guaranteed-eruption convention
    # v1's stretch-volcano growth used. The erupted material's type depends on where it
    # surfaces: still standing above sea level (this node's *pre-melt* elevation, i.e.
    # its own elevation before today's deform() pass touched it) is
    # continental-type magmatism -- real continental rifts stay bimodal-volcanic land
    # for a long stretch before a true ocean opens (the East African Rift, well before
    # the Red Sea stage) -- while a node already at or below sea level (a drowned
    # margin, or an ordinary oceanic ridge) erupts ordinary mid-ocean-ridge oceanic
    # crust. See docs/simulation-model.md's "Magma-typed decompression melting".
    hc_before_melting = hc.copy()
    _erupt_melted_nodes(world, plate.plate_id, rng_index, hc, hm, crust_type_code, is_volcano, volcano_remaining, melting, elevation)
    continental_after_melting = effective_is_continental_from_codes(
        crust_type_code, plate.crust_type == "continental"
    )
    juvenile = np.where(
        melting & continental_after_melting,
        np.maximum(hc - hc_before_melting, 0.0),
        0.0,
    )
    continental_material += juvenile
    continental_ledger.record(
        world, "juvenile_additions_m3", float(np.dot(juvenile, budget_area_m2))
    )
    if np.any(melting):
        mobile_cover.record(world, "rift_reset_m3", float(np.dot(np.where(melting, cover, 0.0), budget_area_m2)))
        cover[melting] = 0.0
        cover_material[melting] = 0.0
    became_oceanic = melting & ~continental_after_melting
    if np.any(became_oceanic):
        removed = np.where(became_oceanic, continental_material, 0.0)
        continental_material[became_oceanic] = 0.0
        continental_ledger.record(
            world, "rift_thinned_m3", float(np.dot(removed, budget_area_m2))
        )
    phase_budget.record(world, plate, "decompression_melting", checkpoint_hc, checkpoint_hm, codes0, hc, hm, crust_type_code, **budget_areas)
    _ignite_early_rift_volcanoes(world, plate.plate_id, rng_index, is_volcano, volcano_remaining, newly_below_rift_onset)

    # Transform (strike-slip) pressure-ridge uplift: a modest, always-transpressional
    # bump on the transform band, kept as a direct elevation delta (like erosion's
    # own contributions) rather than an Hc change -- a strike-slip contact shoulders
    # up local relief without net crustal shortening. Gated by `fault_influence` in
    # "fault" mode so it tracks the boundary strike-slip fault families rather than
    # smearing along the whole polygon edge.
    transform_uplift = np.zeros(n)
    transform_uplift[transform] = TRANSFORM_UPLIFT_RATE_M_PER_MYR * years_myr * fault_influence[transform]

    elevation_after = lithosphere.isostatic_elevation(hc, hm, rho_c)

    # Issue #189 follow-up (see UNBACKED_RELIEF_DECAY_PER_MYR above): relax any
    # *existing* positive debt -- elevation this node already carries in excess of what
    # `elevation_before` says its own Hc/Hm column supports -- toward zero, before this
    # step's own fresh transform_uplift (still a bare delta by design) potentially
    # adds more. The column at the *start* of this step (`fields`, matching what
    # `elevation_before` was computed from), not `hc` (already mutated by the
    # convergent/divergent passes above) -- and v1 nodes with no Hc tracking at all
    # (all-zero) are left alone, same has_column gating lithosphere.back_elevation_gain uses.
    start_hc = fields["crustal_thickness_m"]
    existing_debt = np.where(start_hc > 0.0, np.clip(elevation - elevation_before, 0.0, None), 0.0)
    debt_relief = existing_debt * (1.0 - np.exp(-UNBACKED_RELIEF_DECAY_PER_MYR * years_myr))

    new_elevation = rheology.clip_elevation_bounds(elevation - debt_relief + (elevation_after - elevation_before) + transform_uplift)
    if np.any(melting):
        # The delta above used this plate's single nominal `rho_c` for both
        # `elevation_before`/`elevation_after` -- fine for every ordinary node, whose
        # crust_type_code is still CRUST_TYPE_INHERIT, but wrong for a node that just
        # melted into the *other* type (e.g. a continental plate's drowned margin
        # melting through to real oceanic crust): its fresh Hc/Hm reference column
        # should float at the density of what it actually is now, not the plate's own
        # nominal density. This is a brand-new column with no prior erosion history to
        # preserve (same "hard reset, not a delta" character the Hc/Hm reset above
        # already has), so read it exactly rather than folding it into the delta.
        melt_rho_c = lithosphere.node_crust_density(crust_type_code[melting], plate.crust_type)
        new_elevation[melting] = rheology.clip_elevation_bounds(lithosphere.isostatic_elevation(hc[melting], hm[melting], melt_rho_c))

    # Elevation-change provenance (diagnostic only -- see elevation_lines.ELEV_CHANGE_*
    # and render_image's "elevReason" view). Stamp whichever tectonic process moved a
    # node this step, gated on ELEV_CHANGE_MIN_DELTA_M so a node barely grazed by a
    # fading boundary force keeps its older provenance. The masks partition the
    # near-boundary band by motion (convergent / divergent / transform), so a plain
    # per-mask assignment needs no priority order. `faults._apply_plate_fault_relief`
    # runs after this pass and overwrites these with a FAULT_* code wherever a
    # boundary fault of the matching regime moved the node -- that is what paints the
    # fault families along every boundary in the elevReason view.
    reason = fields["elev_change_reason"].copy()
    moved = np.abs(new_elevation - elevation) >= ELEV_CHANGE_MIN_DELTA_M
    if plate.crust_type == "continental":
        reason[convergent & moved & ~neighbor_oceanic] = ELEV_CHANGE_COLLISION
        reason[convergent & moved & neighbor_oceanic] = ELEV_CHANGE_SUBDUCTION_ARC
        reason[arc_band & moved] = ELEV_CHANGE_SUBDUCTION_ARC
    else:
        reason[convergent & moved] = ELEV_CHANGE_TRENCH
    reason[divergent & moved] = ELEV_CHANGE_RIFT
    # Gated on transform_uplift itself, not just band membership -- issue #189
    # follow-up's debt-decay term can move a transform-band node's elevation on its own
    # (fault_influence == 0 in "fault" mode zeroes transform_uplift there, but leftover
    # debt still decays), which would otherwise mislabel a pure decay move as an active
    # transform pressure ridge.
    reason[(transform_uplift > 0.0) & moved] = ELEV_CHANGE_TRANSFORM
    reason[melting] = ELEV_CHANGE_VOLCANO

    # Cratonic crust thins with its column and is gone wherever the column melted through:
    # both are rifting (cratons.py). Thickening leaves it unchanged.
    new_craton = cratons.scale_with_column(craton, fields["crustal_thickness_m"], hc)
    new_craton[melting] = 0.0
    rifted = craton - new_craton
    if np.any(rifted > 0.0):
        cratons.record(world, "rifted_m3", float(np.sum(rifted * node_area_m2)))

    if strained is not None:
        strained.update(crustal_thickness_m=strained_hc, mantle_lithosphere_thickness_m=strained_hm, melted=melting)

    return {
        "elevation": new_elevation,
        "crustal_thickness_m": hc,
        "mantle_lithosphere_thickness_m": hm,
        "divergent_age_myr": new_age,
        "is_volcano": is_volcano,
        "volcano_active_years_remaining": volcano_remaining,
        "elev_change_reason": reason,
        "crust_type_code": crust_type_code,
        "craton_crust_m": new_craton,
        "continental_material_m": continental_material,
        "mobile_cover_m": cover,
        "mobile_cover_continental_m": np.minimum(cover_material, continental_material),
    }


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

# The `sketch` path's own continental hc_at (see generate_plates) -- replaces the noise-quantile
# `land_threshold` with a direct land/sea call from the drawn coastline (worldsketch.SketchMasks).
# `_SKETCH_LAND_OFFSET_UNITS`/`_SKETCH_SEA_OFFSET_UNITS` stand in for `(sample - land_threshold)`
# in the ordinary formula -- comfortably positive on drawn land, comfortably negative at sea (well
# outside `_UPLIFT_SEA_MARGIN` + `_UPLIFT_SEA_RAMP`'s ~0.40 band either way, so the uplift gate
# below reads land as land regardless of the noise texture riding on top of it), so the sketch --
# not the noise field -- decides the coastline. `_SKETCH_TEXTURE_WEIGHT` still lets
# `relief.sample()` show through at reduced strength for natural-looking local variation on each
# side. `_SKETCH_MOUNTAIN_BONUS_M`/`_SKETCH_RIVER_CARVE_M` are flat elevation-equivalent nudges
# (same "units = m / 2000" convention `_OROGENIC_RELIEF_M` etc. already use) applied only to
# painted-and-land nodes; the river carve is capped so it can never drop a node below the "land"
# side of the sketch (a drawn river is a hint for hydrology.py to find once the world is stepped,
# not a persisted river -- rivers/lakes have no such concept at generation time).
_SKETCH_LAND_OFFSET_UNITS = 0.55
_SKETCH_SEA_OFFSET_UNITS = -0.55
_SKETCH_TEXTURE_WEIGHT = 0.35
_SKETCH_MOUNTAIN_BONUS_M = 1800.0
_SKETCH_RIVER_CARVE_M = 250.0
_SKETCH_MOUNTAIN_BONUS_UNITS = _SKETCH_MOUNTAIN_BONUS_M / 2000.0
_SKETCH_RIVER_CARVE_UNITS = _SKETCH_RIVER_CARVE_M / 2000.0

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
# outline can get.
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


def build_plate_tiling(
    rng: np.random.Generator,
    num_plates: int,
    extra_sites_per_plate: int = EXTRA_SITES_PER_PLATE,
    primary_sites: np.ndarray | None = None,
) -> PlateTiling:
    """Place `num_plates` primary sites plus `num_plates * extra_sites_per_plate` extra sites
    uniformly on the sphere, then hand every extra site to a plate by region-growing: each
    round, the still-unassigned site closest (angularly) to any already-assigned site joins
    that site's plate. Prim-style growth keeps each plate's set of sites a compact cluster,
    so the union of their Voronoi cells stays a single lumpy blob rather than scattering
    disconnected islands across the sphere. Deterministic in `rng`.

    `primary_sites`, when given (unit vectors, shape `(num_plates, 3)`), replaces the usual
    random primary placement -- used by `generate_plates`'s `sketch` path
    (`worldsketch.sketch_plate_sites`) to seed plates from a drawn/loaded coastline's
    landmasses and open ocean instead of scattering them uniformly at random. `None` (every
    caller before this parameter existed) draws primaries the same single `rng.normal` call as
    always, so that path is unchanged."""
    num_extra = max(0, num_plates * extra_sites_per_plate)
    if primary_sites is None:
        site_xyz = rng.normal(size=(num_plates + num_extra, 3))
        site_xyz /= np.linalg.norm(site_xyz, axis=-1, keepdims=True)
    else:
        primaries = np.asarray(primary_sites, dtype=float)
        if len(primaries) != num_plates:
            raise ValueError(f"primary_sites must have exactly num_plates={num_plates} rows, got {len(primaries)}")
        primaries = primaries / np.linalg.norm(primaries, axis=-1, keepdims=True)
        if num_extra > 0:
            extra_xyz = rng.normal(size=(num_extra, 3))
            extra_xyz /= np.linalg.norm(extra_xyz, axis=-1, keepdims=True)
            site_xyz = np.concatenate([primaries, extra_xyz], axis=0)
        else:
            site_xyz = primaries

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
    voronoi_points: int | None = None,
    sketch: worldsketch.SketchMasks | None = None,
    premade_world_id: str | None = None,
) -> list[Plate]:
    """`plates.generate_plates`'s own seed-placement/Voronoi-tiling algorithm, extended so
    each plate owns the union of several adjacent Voronoi cells (see `build_plate_tiling` and
    `EXTRA_SITES_PER_PLATE`) rather than a single cell -- still deterministic per `seed`, still
    nearest-site-owns-the-node so the tiling has no gaps/overlaps by construction, just with
    lumpier, less convex plate outlines. Only the per-plate line-building step differs from
    v1: each node gets a reference Hc/Hm plus a composite relief field on Hc (see
    `terrain_noise.py` -- a low-frequency `sample()` that decides land/sea exactly as v1's
    single noise did, plus a land-gated non-negative `uplift()` carrying orogenic belts and
    plateaus; `_HC_NOISE_AMPLITUDE_*`/`_OROGENIC_*`/`_PLATEAU_*` above set the amplitudes),
    with `elevation` itself computed once via isostasy at the end.

    `sketch` (the "Human-made" Generate World tab, see `worldsketch.py`), when given, replaces
    two independent pieces of the usual random generation rather than running alongside it:
    plate *sites* come from `worldsketch.sketch_plate_sites` (landmasses -> continental sites,
    open ocean -> farthest-point-sampled oceanic sites) instead of `build_plate_tiling`'s own
    random placement, and each continental plate's `hc_at` reads land/sea (and, where painted,
    mountain/river) straight off the sketch instead of a noise-quantile threshold -- see the
    `_SKETCH_*` constants above. `land_fraction` is ignored in this case (the drawing decides
    land extent directly); `continental_fraction` keeps its meaning, now controlling how many
    plate slots the drawn landmasses are split across vs. pure ocean. Oceanic `hc_at` is
    untouched either way -- sketch masks only ever bias continental crust, the same way
    `land_fraction` already only touches the continental formula.

    `voronoi_points` (the UI's "Voronoi points" Advanced-settings slider), when given, is the
    *total* number of Voronoi seed points to scatter -- it overrides `extra_sites_per_plate`,
    which is derived once the final plate count is known as
    `round(voronoi_points / num_plates) - 1` (clamped at 0). More points -> lumpier, less
    convex plate outlines; `voronoi_points <= num_plates` recovers the one-cell-per-plate
    tiling. `None` keeps `extra_sites_per_plate` as passed.

    `premade_world_id` (the Generate World dialog's "Premade worlds" tab -- `"earth"`,
    `"pangaea"`, or `"got"`; `None` for every other tab) always swaps `ContinentalRelief`'s
    belt/plateau masks for real-geography/lore ones (see relief_regions.py) instead of the
    usual calibrated-random coverage -- see that module and `terrain_noise.ContinentalRelief`
    for what a belt/plateau mask actually changes (where one is allowed to appear, not the
    ridge/terrace texture inside it).

    `"earth"` additionally bypasses `sketch`/`num_plates`/`continental_fraction`/
    `voronoi_points` entirely for plate/coastline construction: `real_plates.
    build_exact_earth_plates` partitions the globe directly from the real 16-plate boundaries
    and a real coastline raster (see that function's own module comment for why -- the old
    Voronoi-approximated site placement below, kept for Pangaea/GoT, measurably merged real
    straits/erased real islands no site-count tuning could fix), so both plate count and land
    extent are simply whatever that partition contains, not a knob here. `sketch` is still
    parsed from whatever the frontend sent (every premade world still supplies one, even
    "earth") but only for mountain/river painting now -- see `real_plates.exact_earth_masks`,
    which is what actually reaches each continental plate's `hc_at` below.

    `"pangaea"` still replaces `sketch`'s own `sketch_plate_sites` for *site placement* with
    real (rigidly transformed) plate geometry (`real_plates.real_plate_sites`) the old way --
    `sketch` is still required alongside it and still decides land/sea/mountain/river in each
    continental plate's `hc_at` exactly as it does for "Human-made". `"got"` has no real-plate
    analog, so its site placement falls through to the ordinary sketch-driven path, same as
    "Human-made"."""
    rng = np.random.default_rng(seed)
    if num_plates is None:
        num_plates = int(rng.integers(MIN_AUTO_PLATES, MAX_AUTO_PLATES + 1))

    num_continents: int | None = None
    if continental_fraction is not None:
        continental_fraction = max(0.0, min(continental_fraction, 1.0))
        num_continents = round(continental_fraction * num_plates)
        num_plates = max(num_plates, num_continents + MIN_OCEANIC_PLATES)

    # Premade-worlds' real-geography relief regions (relief_regions.py) -- unlike the site-
    # placement override below, these apply for *every* premade world including "got" (Z&D
    # has no real-world data to ground plate placement/mantle convection in, but its named
    # lore regions still replace the usual random belt/plateau coverage).
    belt_mask = plateau_mask = None
    if premade_world_id is not None:
        # Local imports: real_plates.py/relief_regions.py are the newer, Premade-worlds-
        # specific modules; keeping the import here (rather than at module scope) avoids
        # paying for real_plates' own data-file parsing on every other generation path.
        from . import relief_regions

        if premade_world_id == "earth":
            belts, plateaus = relief_regions.EARTH_BELTS, relief_regions.EARTH_PLATEAUS
        elif premade_world_id == "pangaea":
            belts, plateaus = relief_regions.pangaea_belts(), relief_regions.pangaea_plateaus()
        elif premade_world_id == "got":
            belts, plateaus = relief_regions.GOT_BELTS, {}
        else:
            raise ValueError(f"unknown premade_world_id {premade_world_id!r}")
        belt_mask = relief_regions.build_belt_mask(belts)
        plateau_mask = relief_regions.build_plateau_mask(plateaus)

    # Non-None only for "earth" -- see its own module comment in real_plates.py. Its
    # `is_owned`/frame-seed come straight from a baked exact partition instead of the
    # Voronoi-tiling `tiling`/`owner_tree` every other path below builds, so most of the rest
    # of this function branches on whether this is set rather than touching `tiling` directly.
    exact_earth_plates = None

    if premade_world_id == "earth":
        from . import real_plates

        exact_earth_plates = real_plates.build_exact_earth_plates()
        num_plates = len(exact_earth_plates.crust_types)
        crust_types = exact_earth_plates.crust_types
        # The sketch's own (traced, resolution-limited) land grid is only still used for
        # mountain/river painting now -- land/sea itself comes from the exact partition's own
        # much finer raster instead (see exact_earth_masks' own docstring).
        sketch = real_plates.exact_earth_masks(sketch)
    elif premade_world_id == "pangaea":
        # "Predetermined plates": site placement from real plate geometry instead of the
        # sketch's own landmasses (see real_plates.py). Not (yet) rebuilt on the same exact
        # partition "earth" uses above -- Pangaea's plates/coastline are already an
        # approximation (a rigid per-continent transform, not a real reconstruction), so
        # Voronoi-approximated placement within that is a smaller loss than it is for Earth.
        from . import real_plates

        real_plate_list = real_plates.pangaea_real_plates()
        target_continents = num_continents if num_continents is not None else round(CONTINENTAL_FRACTION * num_plates)
        site_xyz, crust_types = real_plates.real_plate_sites(
            sketch, real_plate_list, num_plates, target_continents, rng, pooled_oceanic=True
        )
        num_plates = len(site_xyz)
        if voronoi_points is not None:
            extra_sites_per_plate = max(0, round(voronoi_points / max(num_plates, 1)) - 1)
        tiling = build_plate_tiling(rng, num_plates, extra_sites_per_plate, primary_sites=site_xyz)
    elif sketch is not None:
        target_continents = num_continents if num_continents is not None else round(CONTINENTAL_FRACTION * num_plates)
        site_xyz, crust_types = worldsketch.sketch_plate_sites(sketch, num_plates, target_continents, rng)
        num_plates = len(site_xyz)
        if voronoi_points is not None:
            extra_sites_per_plate = max(0, round(voronoi_points / max(num_plates, 1)) - 1)
        tiling = build_plate_tiling(rng, num_plates, extra_sites_per_plate, primary_sites=site_xyz)
    else:
        if voronoi_points is not None:
            extra_sites_per_plate = max(0, round(voronoi_points / max(num_plates, 1)) - 1)
        tiling = build_plate_tiling(rng, num_plates, extra_sites_per_plate)
        if num_continents is None:
            crust_types = ["continental" if rng.random() < CONTINENTAL_FRACTION else "oceanic" for _ in range(num_plates)]
        else:
            continental_indices = set(rng.choice(num_plates, size=num_continents, replace=False).tolist())
            crust_types = ["continental" if i in continental_indices else "oceanic" for i in range(num_plates)]

    # `exact_earth_plates` ("earth" only) has no `tiling`/`owner_tree` at all -- its `is_owned`
    # is a direct label-grid lookup (see the per-plate loop below) -- and `land_threshold`
    # below is never reached for it either (that branch requires `sketch is None`, never true
    # for a premade world), so neither `site_crust_types` nor `owner_tree` has a use to build
    # for it.
    if exact_earth_plates is None:
        seed_xyz = tiling.site_xyz

        # Per-site crust type (each site inherits its owning plate's) -- `_land_noise_threshold`
        # and the `is_owned` test below both index by nearest *site*, not nearest plate.
        site_crust_types = [crust_types[tiling.site_plate[s]] for s in range(len(seed_xyz))]

        owner_tree = cKDTree(seed_xyz)

    # Composite relief fields (see terrain_noise.py) -- the last consumers of `rng`, drawn in
    # a fixed order so a given seed reproduces the same terrain. `relief.sample()` stands in
    # for the old single `SphereNoise` (same std, same land/sea decision); `relief.uplift()`
    # adds the orogenic belts and plateaus, land-gated in `hc_at` below. Still drawn even in
    # the `sketch` path -- the sketch's own hc_at still layers this texture on top (see
    # `_SKETCH_TEXTURE_WEIGHT`), and skipping the draw would shift every later `rng` call.
    relief = terrain_noise.ContinentalRelief(
        rng,
        orogenic_units=_OROGENIC_RELIEF_UNITS,
        plateau_units=_PLATEAU_UPLIFT_UNITS,
        plateau_relief_units=_PLATEAU_RELIEF_UNITS,
        belt_mask=belt_mask,
        plateau_mask=plateau_mask,
    )
    ocean_relief = terrain_noise.OceanicRelief(rng)

    land_threshold = None
    if sketch is None and land_fraction is not None:
        land_fraction = max(0.0, min(land_fraction, 1.0))
        land_threshold = _land_noise_threshold(
            owner_tree, site_crust_types, relief, land_fraction, _continental_sealevel_noise_offset()
        )

    spacing_rad = line_spacing_rad(node_density)
    plates: list[Plate] = []
    for i in range(num_plates):
        crust_type = crust_types[i]
        hc0, hm0 = lithosphere.reference_thickness(crust_type)
        hc_amp = _HC_NOISE_AMPLITUDE_CONTINENTAL_M if crust_type == "continental" else _HC_NOISE_AMPLITUDE_OCEANIC_M

        if exact_earth_plates is None:
            frame = geometry.plate_frame_from_seed(tiling.primary_site(i))

            def is_owned(world_pts: np.ndarray, _i: int = i) -> np.ndarray:
                _, nearest_idx = owner_tree.query(world_pts)
                return tiling.site_plate[nearest_idx] == _i
        else:
            frame = geometry.plate_frame_from_seed(exact_earth_plates.frame_xyz[i])

            def is_owned(world_pts: np.ndarray, _i: int = i) -> np.ndarray:
                return exact_earth_plates.is_owned(_i, world_pts)

        if crust_type == "continental":
            if sketch is not None:

                def hc_at(world_pts: np.ndarray, _hc0: float = hc0, _amp: float = hc_amp, _sketch: worldsketch.SketchMasks = sketch) -> np.ndarray:
                    s = relief.sample(world_pts)
                    is_land = _sketch.sample_land(world_pts)
                    base_units = np.where(is_land, _SKETCH_LAND_OFFSET_UNITS, _SKETCH_SEA_OFFSET_UNITS)
                    combined = base_units + _SKETCH_TEXTURE_WEIGHT * s
                    hc = _hc0 + _amp * combined
                    gate = np.clip((combined - _UPLIFT_SEA_MARGIN) / _UPLIFT_SEA_RAMP, 0.0, 1.0)
                    hc = hc + _amp * gate * relief.uplift(world_pts)
                    mountain = _sketch.sample_mountain(world_pts) & is_land
                    if np.any(mountain):
                        hc = hc + np.where(mountain, _amp * _SKETCH_MOUNTAIN_BONUS_UNITS, 0.0)
                    river = _sketch.sample_river(world_pts) & is_land
                    if np.any(river):
                        min_hc = _hc0 + _amp * (_SKETCH_LAND_OFFSET_UNITS * 0.25)
                        hc = np.where(river, np.maximum(hc - _amp * _SKETCH_RIVER_CARVE_UNITS, min_hc), hc)
                    return hc
            else:
                _lt = 0.0 if land_threshold is None else land_threshold

                def hc_at(world_pts: np.ndarray, _hc0: float = hc0, _amp: float = hc_amp, _lt: float = _lt) -> np.ndarray:
                    s = relief.sample(world_pts)
                    gate = np.clip((s - _lt - _UPLIFT_SEA_MARGIN) / _UPLIFT_SEA_RAMP, 0.0, 1.0)
                    return _hc0 + _amp * (s - _lt) + _amp * gate * relief.uplift(world_pts)
        else:

            def hc_at(world_pts: np.ndarray, _hc0: float = hc0, _amp: float = hc_amp) -> np.ndarray:
                return _hc0 + _amp * ocean_relief.sample(world_pts)

        plate = PlateWithSparseQuadPatch.from_lattice(i, frame, crust_type, spacing_rad, is_owned)
        world_pts = plate.all_points_and_elevation()[0]
        # Also upper-clipped (issue #161): the noise term alone can occasionally seed a
        # continental node above MAX_CRUSTAL_THICKNESS_M (confirmed directly -- generation
        # noise landed at 85,349 m on one seed), which is generation-time noise, not real
        # tectonic mass, so it's simply clipped rather than conserved/redistributed.
        hc = np.clip(hc_at(world_pts), lithosphere.MIN_CRUSTAL_THICKNESS_M, lithosphere.MAX_CRUSTAL_THICKNESS_M)
        plate.set_fields_on_plate(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=np.full(len(hc), hm0))
        lithosphere.sync_plate_elevation(plate)
        plates.append(plate)

    # Each plate above was seeded purely from its own crust type: submerged continental crust
    # is a uniform bright shelf however far from land it sits, and every continent/ocean
    # boundary is a vertical cliff. Drown the offshore continental interiors and grade the
    # boundary steps into slopes -- see bathymetry.shape_initial_bathymetry.
    bathymetry.shape_initial_bathymetry(plates)
    return plates


def new_plate(
    plate_id: int,
    frame: np.ndarray,
    crust_type: str,
    spacing_rad: float,
    seed: int,
    is_owned=None,
    node_is_continental=None,
) -> Plate:
    """A brand-new `PlateWithSparseQuadPatch` seeded with reference Hc/Hm plus the same
    composite relief field `generate_plates` uses (see `terrain_noise.py`). Keyed off
    `(seed, plate_id, _TERRAIN_SEED_TAG)` so the crust stays attached to this plate.

    `is_owned` (default: every node in `frame`'s entire local lattice) restricts which nodes
    of that lattice actually become part of the plate -- see `gaps.py`'s use of this to carve
    out just one uncovered region rather than claiming the whole sphere.

    `node_is_continental` (default: every node matches `crust_type`) lets a caller seed a
    genuinely mixed-composition plate node-by-node -- see `gaps.py`'s land-adjacent gap-fill,
    which spawns continental nodes right at a coastline and oceanic ones everywhere else in
    the same gap. Each node gets the reference thickness/relief/density of its *own* decided
    type (`elevation_lines.CRUST_TYPE_OCEANIC`/`CONTINENTAL`, stamped explicitly since these
    are the authoritative source of that node's real composition, not an inherited label), and
    the plate's own `crust_type` is the majority, by cell area, of what actually got built
    rather than the `crust_type` argument verbatim --
    that argument is only the fallback/default when every node ends up the same type, which is
    every caller before this parameter existed."""
    if is_owned is None:

        def is_owned(world_pts: np.ndarray) -> np.ndarray:
            return np.ones(len(world_pts), dtype=bool)
    if node_is_continental is None:
        _default_continental = crust_type == "continental"

        def node_is_continental(world_pts: np.ndarray) -> np.ndarray:
            return np.full(len(world_pts), _default_continental)

    rng = np.random.default_rng((seed, plate_id, _TERRAIN_SEED_TAG))
    hc0_continental, hm0_continental = lithosphere.reference_thickness("continental")
    hc0_oceanic, hm0_oceanic = lithosphere.reference_thickness("oceanic")
    continental_relief = terrain_noise.ContinentalRelief(
        rng,
        orogenic_units=_OROGENIC_RELIEF_UNITS,
        plateau_units=_PLATEAU_UPLIFT_UNITS,
        plateau_relief_units=_PLATEAU_RELIEF_UNITS,
    )
    oceanic_relief = terrain_noise.OceanicRelief(rng)

    def hc_at(world_pts: np.ndarray, is_continental: np.ndarray) -> np.ndarray:
        hc = np.empty(len(world_pts))
        if np.any(is_continental):
            pts = world_pts[is_continental]
            s = continental_relief.sample(pts)
            gate = np.clip((s - _UPLIFT_SEA_MARGIN) / _UPLIFT_SEA_RAMP, 0.0, 1.0)
            hc[is_continental] = hc0_continental + _HC_NOISE_AMPLITUDE_CONTINENTAL_M * s + _HC_NOISE_AMPLITUDE_CONTINENTAL_M * gate * continental_relief.uplift(pts)
        if not np.all(is_continental):
            pts = world_pts[~is_continental]
            hc[~is_continental] = hc0_oceanic + _HC_NOISE_AMPLITUDE_OCEANIC_M * oceanic_relief.sample(pts)
        return hc

    cells = PlateWithSparseQuadPatch.from_lattice(plate_id, frame, crust_type, spacing_rad, is_owned)
    world_pts = cells.all_points_and_elevation()[0]
    is_continental = np.asarray(node_is_continental(world_pts), dtype=bool)
    areas = cells.node_areas_m2()
    continental_area, oceanic_area = float(areas[is_continental].sum()), float(areas[~is_continental].sum())
    majority = crust_type
    if not np.isclose(continental_area, oceanic_area):
        majority = "continental" if continental_area > oceanic_area else "oceanic"
    plate = PlateWithSparseQuadPatch(
        plate_id,
        frame,
        majority,
        cells.cells_per_edge,
        cells.cell_keys,
        fields={
            # Also upper-clipped (issue #161) -- see generate_plates' own note.
            "crustal_thickness_m": np.clip(hc_at(world_pts, is_continental), lithosphere.MIN_CRUSTAL_THICKNESS_M, lithosphere.MAX_CRUSTAL_THICKNESS_M),
            "mantle_lithosphere_thickness_m": np.where(is_continental, hm0_continental, hm0_oceanic),
            "crust_type_code": np.where(is_continental, CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC).astype(np.int8),
        },
    )
    lithosphere.sync_plate_elevation(plate)
    return plate

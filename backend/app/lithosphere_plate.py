"""`LithospherePlate`: the plate representation. Subclasses `PlateWithLines` (see
plates.py) rather than reinventing its plumbing -- everything not about *how a plate moves or
deforms* (outline tracing, containment, neighbour search, node iteration, the whole render/
erosion/hydrology/stats/persistence surface) is inherited unchanged, since every one of those
consumers only ever reads `Plate`'s abstract interface (`all_points_and_elevation`, `collect`,
`contains_batch`, ...), never `shift`/`deform`'s own internals. See the plan's "Why subclass
PlateWithLines" section.

Defined here: `shift` (torque integration, torque.py), `deform` (Mohr-Coulomb/isostasy,
rheology.py + lithosphere.py), `merge_with`/`_merge_nodes_with`/`split` (carrying
Hc/Hm through the representation's own resample/partition operations). `PlateWithLines`
itself now carries only the geometry/representation half of the interface -- the tectonic
deformation engine lives entirely on this subclass.
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
    CRUST_TYPE_CONTINENTAL,
    CRUST_TYPE_OCEANIC,
    COVERAGE_RADIUS_MULT,
    DEFRAG_CONNECT_RADIUS_MULT,
    ElevationLine,
    IRREGULARITY_TOLERANCE,
    PLANET_RADIUS_KM,
    build_lines_from_lattice,
    line_spacing_rad,
    majority_crust_type,
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
    query_workers,
)
from . import bathymetry, lithosphere, magma_transport, mantle, phase_budget, rheology, terrain_noise, torque, worldsketch

EXTEND_THRESHOLD_MULTIPLIER = 1.3  # same shape as v1's plates.EXTEND_THRESHOLD_RAD
MAX_EXTEND_NODES_PER_STEP = 400

# `_claim_adjacent_territory`'s own mass-conservation share (the phi/between-row half of
# rift-stretch closing -- see `_grow_or_shrink_line_for_deform`'s theta/within-row half for
# the other): rather than seed a brand-new row at the full oceanic reference column "for
# free," the volume it would otherwise get is drawn down across this many of the plate's own
# existing outermost rows on that side *plus* the new row itself, in proportion to how much of
# a full row-spacing's worth of genuine phi-direction separation is behind the claim (see
# rheology.stretch_components' own phi_gap). A small constant, not node-count-scaled, like
# CONTINENTAL_CONTESTED_RETREAT_MIN_RUN above -- this is about how many *rows* share a
# stretching event's mass deficit, not a distance or a per-node budget.
K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION = 2

# `_claim_adjacent_territory` used to claim exactly one new row per call, regardless of how
# much genuine phi-direction separation was actually behind it -- fine for an ordinary
# two-plate rift (next step's call claims the next row in turn), but at a triple junction,
# where three independently-oriented plate grids are all receding from a shared point, one row
# per plate per step structurally cannot keep pace with three-way divergence (confirmed on a
# real save: the void between three plates widened every step even with every other growth
# knob loosened). Looping lets one call claim as many rows as the actual gap calls for, capped
# so a single pathological step can't claim an unbounded amount -- same node-count-budget
# shape as MAX_EXTEND_NODES_PER_STEP, not a distance.
MAX_CLAIM_ROWS_PER_STEP = 4
MAX_CLAIM_NODES_PER_STEP = 400

# `_fill_corner_notch_frontier`'s own window and safety valve. Both `_stretch_end`
# (theta-aligned) and the generalized `_claim_adjacent_territory` (phi-aligned, whole rows) are
# still confined to this plate's own local (theta, phi) grid axes -- a genuinely diagonal
# sub-row notch, the shape independently-oriented plate grids leave open at a triple junction
# *or* a curved two-plate rift boundary, is invisible to both. `_fill_corner_notch_frontier` is
# meant to be the generalized closer for that whole family, so its window has to cover whatever
# a boundary could genuinely have opened up in one step, not a fixed, resolution-only row count -- a flat
# constant here structurally falls behind on a fast-diverging boundary regardless of how
# generous it is (confirmed on a real save, seed 430031492: a 2-plate rift's diagonal residual
# sat up to ~8x spacing deep with the old CORNER_NOTCH_WINDOW_ROWS=5 window, so it could never
# be fully reached no matter the per-step node budget). `mantle.MAX_PLATE_RATE * years` is the
# simulation's own hard per-plate speed ceiling applied to this step's real duration -- the
# worst-case distance *any* plate's edge could have advanced, not just this plate's actual
# current speed, so this stays sufficient even for a plate that speeds up later. Since this
# runs once per plate per step, both sides of a boundary reaching this far from their own edge
# already sums to the worst-case *combined* divergence without needing to double it here. The
# `CORNER_NOTCH_MIN_WINDOW_ROWS` floor keeps today's minimum footprint (a coarse world or a
# near-zero `years` step should never shrink the window below the old baseline).
CORNER_NOTCH_MIN_WINDOW_ROWS = 5
# How far a real neighbour may sit from a candidate notch point and still justify claiming it
# -- see the "requires a real neighbour" guard in _fill_corner_notch_frontier's own docstring
# for why this exists at all (without it, a lone plate's every open perimeter point looks
# claimable). Wider than the window itself (see `window_rad` in `_fill_corner_notch_frontier`)
# by one ordinary boundary-adjacency reach, since the neighbour that makes a deep interior notch-point genuine
# can itself be one of the *other* plates at the same junction, not necessarily the closest one
# to this exact point.
CORNER_NOTCH_NEIGHBOUR_REACH_MARGIN_MULT = torque.BOUNDARY_FORCE_REACH_MULTIPLIER
# A node count, not a distance -- comparable in scale to merge_split.DEFRAG_FRAGMENT_MIN_NODES
# at the old fixed window size, since this path has weaker geometric guarantees than the
# row-based claim above and was meant to stay a rare corner cleanup, never a bulk grower. Now
# that the window itself scales with a step's worst-case opening (see above), this is kept as
# a floor, not the cap: the real per-step budget scales with the window's own area so a
# genuinely large single-step opening isn't starved by a resolution-only constant, while this
# floor still bounds a call whose window collapses to the old minimum. The real limit stays
# geometric (only points within DEFRAG_CONNECT_RADIUS_MULT of an already-connected point --
# see the frontier-hop loop below -- are ever candidates at all); this is defense in depth,
# not the physical constraint.
MIN_CORNER_FILL_NODES_PER_STEP = 200

# A plate whose boundary is more deeply/widely overlapping a neighbour right now should
# crumple faster than one barely grazing -- on top of (not instead of) the existing
# distance-decay shape within the belt. `contested_all.mean()` (this plate's own fraction of
# near-boundary-band nodes currently classified contested) is a free, already-computed
# per-step severity signal; at severity 1.0 (the whole band contested -- a deep pile-up) the
# near-field contested-band strength is tripled.
OVERLAP_UPLIFT_SEVERITY_GAIN = 2.0

# Transform (strike-slip) boundary pressure-ridge uplift, applied as a direct elevation delta
# on the transform band in deform() (there is no net crustal shortening at a strike-slip
# contact, so it does not go through Hc/Hm). Gentler than either convergent case -- real
# transform relief is local pressure ridges / transtensional sags, not an orogen. Half v1's
# plates.TRANSFORM_UPLIFT_RATE_M_PER_MYR (200) since here it is not distance-tapered, only
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
# without dedicated erosion attention. So any *existing* positive debt on a line decays toward
# zero every step (this step's own fresh contribution is added after, so it isn't clawed back
# before it even shows up) -- same age/relax-toward-target idiom
# `rheology.relax_young_oceanic_mantle_lithosphere` already uses for a different field. Picked
# so a node sitting continuously in a fully-active transform band (worst case, no relief from
# ever leaving the band) settles at a steady-state debt on the order of a few thousand meters
# -- consistent with transform relief being "modest"/"local" by design, not a second orogeny
# -- rather than drifting arbitrarily far before the hard clip.
UNBACKED_RELIEF_DECAY_PER_MYR = 0.05

# A continental line's *contested* end is allowed to retreat -- one node per step -- whether
# the overriding neighbour is oceanic (a passive margin / accretion front: the ocean slab
# descends under it and the buried continental node cedes nothing the model should keep) or
# continental (a suture whose overlapping crust is consumed into the orogen -- its volume is
# not discarded but thrust back onto the plate's own surviving leading edge, see
# `_redistribute_accreted_column` / SUTURE_ACCRETION_SPREAD_NODES, so the belt builds real
# relief in proportion to the overlap it actually eats). Left un-retreatable, a contested end
# still grows at its *other* (divergent) side every step and never back -- the continental
# node ratchet that drives the unbounded node-count creep and the slow land-fraction decline
# (GitHub issue #119, "Node-count creep") -- and, for a continent-continent pile-up, a deep
# territory overlap that just sat there for tens of Myr until the forced-merge timer fused
# the pair (the `overlapAge` view's stalled multi-plate collisions). Retreat is gated:
#   - one node per step (the existing `n_distance_cap` / `max_extend_nodes` caps already do
#     this at continental drift rates), and
#   - only where the contested node is part of a run of at least this many consecutive
#     contested nodes, so a single stray node from bounding-polygon envelope fuzz can't
#     nibble a stable coastline or, worse, sever a lobe into a spurious defragmentation plate
#     (the failure the naive "retreat every continental contested node" experiment hit -- see
#     GitHub issue #119; the interior-subduction carve below also stays oceanic-only for the
#     same lobe-severing reason).
CONTINENTAL_CONTESTED_RETREAT_MIN_RUN = 3

# Whole-row retreat -- the reverse of `_claim_adjacent_territory`, and the *only* retreat op
# available in the "parallel suture" regime: when a neighbour overrides a continental plate's
# frontmost phi-row over its full theta width there is no uncontested end for
# `_grow_or_shrink_line_for_deform` to trim, and a continental row is never carved mid-span
# (that severs the landmass into a spurious defragmentation plate), so end-trim alone leaves
# that plate physically unable to give ground -- its trailing edge still grows every step, so
# the node pile ratchets outward regardless (GitHub issue #119, "Continental ratchet: solution
# design", mechanism 2). Once a plate's outermost row (either phi extreme) has been at least
# LEADING_ROW_CONTESTED_FRACTION contested for a cumulative LEADING_ROW_RETREAT_SUSTAINED_YEARS
# of deform time, the whole row is dropped. Whole-row removal keeps the plate contiguous -- the
# lobe-severing hazard is specific to *mid*-row carving -- so this is safe exactly where the
# interior-subduction carve is not.
#
# The dropped row's crustal volume is conserved (mass-conservation bug, confirmed 2026-09-12
# investigating a "land fraction rises then falls" report on the `two_colliding_pairs`
# Debugging World: a real, if quantitatively minor -- a few tenths of a percent of total
# continental crust volume over 30+ My in that scenario -- silent mass leak, since this used
# to just drop the row's Hc/Hm with nothing thrusting it anywhere, unlike every other
# continental retreat path in this file, which conserves via `_redistribute_accreted_column`).
# See `_accrete_dropped_row_volume`: it thrusts each dropped row's summed Hc onto the plate's
# new leading row at that same extreme, spread evenly across that whole row's nodes (a
# whole-row-wide suture, since the row being consumed spans the plate's full theta width --
# unlike `_redistribute_accreted_column`'s single-line-end case, there is no narrower "end"
# to concentrate it on), capped at the same `SUTURE_ACCRETION_MAX_HC_M` (the overflow
# delaminates, same as ordinary suture accretion).
LEADING_ROW_CONTESTED_FRACTION = 0.7
LEADING_ROW_RETREAT_SUSTAINED_YEARS = 5_000_000.0
# Never drop a row that would take the plate below this many rows -- a tiny plate has no
# "leading row" worth the name and the contiguity argument gets thin.
LEADING_ROW_DROP_MIN_ROWS = 4
# Volume-budget growth gate (GitHub issue #119, "Continental ratchet: solution design",
# mechanism 1). A lattice node's physical footprint is constant across the sphere by
# construction (`lithosphere.node_area_m2`), so a plate's total area is just its node count
# times that -- and its implied *mean* crustal thickness is `mean(crustal_thickness_m)`.
# The continental boundary ratchet dilutes this: `_grow_or_shrink_line_for_deform` and
# `_claim_adjacent_territory` seed every new margin node at the *oceanic* reference column
# (`growth_seed_thickness`), and nothing ever removes a whole leading row, so a
# shear-stretched continental plate tiles unbounded drowned passive-margin outward -- node
# count creeps ~+5-6% per 150 My and the plate interior isostatically oceanises into a
# "giant 80%-drowned continent" (GitHub issue #119 items 2 / 5, and issue #120's land-fraction decline).
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
# `deform()` computes that total once (mirroring the calculation `apply_arc_magmatic_
# thickening` performs per-line for real) into a shared budget, and `_grow_or_shrink_line_
# for_deform` spends it down end by end via `_budget_limited_removal`. Node area is constant
# per node, so summed Hc is directly comparable as "volume" on both sides of this cap, the
# same convention `_redistribute_accreted_column` already uses. A continent-continent suture
# retreat is untouched by this (it already fully conserves its own volume, so it isn't the
# uncapped channel this targets), and so is an oceanic self-plate's own ordinary subduction
# (not a passive margin at all -- that loss is expected, not a bug).
#
# 1.0 is a genuinely symmetric cap (this issue's own proposal: match the two channels' rates
# directly), kept as a named multiplier rather than folded into the budget calculation itself
# so it is the one place to loosen this if a strict 1:1 turns out to starve ordinary retreat
# unrelated to the runaway this exists to fix.
OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER = 1.0

# GitHub issue #216: an oceanic self-plate's own ordinary subduction (this constant's own
# neighbour above used to call it "not a passive margin at all -- that loss is expected, not a
# bug") is real geology -- subducted oceanic lithosphere genuinely leaves the surface reservoir
# -- but the issue's own phase-budget instrumentation (#222) found its *rate* has exactly the
# same uncapped-loss/capped-gain shape #177 already fixed for the continental case, just in the
# one retreat path that fix didn't reach: bounded only by the geometric/count caps
# (n_distance_cap/max_extend_nodes), with nothing tying how much Hc a step's subduction removes
# to how much this same plate actually created that step. Measured across two windows (issue
# comment), this uncompensated loss was the single largest net Hc/Hm driver once the
# self-canceling stretch/regularize bookkeeping pair (line_end_stretch vs line_regularization)
# was excluded.
#
# Arc magmatism doesn't apply here -- ARC_MARGIN_SEED_*/arc_band_all are continental-only, since
# arc volcanism happens on the *overriding* plate in this model (growth_seed_thickness's own
# land-runaway history is why only a continental leading edge ever seeds arc crust). The oceanic
# analog of "genuinely new mass added from the mantle this step" is decompression-melting
# eruption at this plate's own divergent/rift nodes (the `melting` mask
# `rheology.apply_divergent_deformation` already returns) -- the only mechanism that adds Hc to
# an oceanic self-plate at all. Same rationale and same 1:1 multiplier as
# OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER: cap the loss side by this step's real creation
# rather than invent an independent rate for it. Shares the same budget-spend plumbing
# (`_budget_limited_removal`, the `oceanic_override_retreat_budget_hc` array threaded through
# `_grow_or_shrink_line_for_deform`) -- only one of the two budget computations below ever
# applies to a given plate, so reusing one array/parameter name is safe.
OCEANIC_SELF_RETREAT_BUDGET_MULTIPLIER = 1.0


def _budget_limited_removal(hc: np.ndarray, n_remove: int, budget_hc: np.ndarray, from_high: bool) -> int:
    """How many of the candidate `n_remove` end nodes a budget-limited retreat may actually
    take this step, capped by `budget_hc[0]` (this plate's remaining same-step creation budget
    -- see deform()'s own `oceanic_override_retreat_budget_hc`, filled from arc-magmatic
    creation for a continental self-plate's oceanic-override retreat
    (OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER, issue #177) or from decompression-melting
    creation for an oceanic self-plate's own ordinary subduction
    (OCEANIC_SELF_RETREAT_BUDGET_MULTIPLIER, issue #216)). `budget_hc` is a shared, mutable
    single-element array -- every call across this plate's whole deform() pass (both ends,
    every line) spends down the same plate-wide budget in place, in call order; once it hits
    zero, every budget-limited retreat still to come this step is refused, not just throttled.

    Finds the largest prefix of the candidate window, counted in from the true (retreating)
    end, whose summed Hc stays within the remaining budget -- node area is constant per node,
    so that sum is directly comparable to the Hc-sum budget. Returns a possibly-smaller
    `n_remove` (down to 0) and debits whatever it actually spends from `budget_hc`."""
    if n_remove <= 0:
        return n_remove
    tail = hc[-n_remove:] if from_high else hc[:n_remove]
    ordered = tail[::-1] if from_high else tail  # index 0 = the true (retreating) end
    cum = np.cumsum(ordered)
    budget = max(float(budget_hc[0]), 0.0)
    within = cum <= budget
    k = n_remove if within.all() else int(np.argmax(~within))
    if k > 0:
        budget_hc[0] -= float(cum[k - 1])
    return k


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
# GitHub issue #120, "Land fraction slowly declines"). The seed lands as shallow forearc/shelf
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


def _erupt_melted_nodes(
    world: "World",  # noqa: F821
    plate_id: int,
    line_index: int,
    hc: np.ndarray,
    hm: np.ndarray,
    crust_type_code: np.ndarray,
    is_volcano: np.ndarray,
    volcano_remaining: np.ndarray,
    melting: np.ndarray,
    prior_elevation: np.ndarray,
) -> None:
    """In-place: any node flagged `melting` (its Hc just crossed below
    `RIFT_CRITICAL_THICKNESS_M` -- from ordinary divergent thinning
    (`rheology.apply_divergent_deformation`) or from rift-stretch thinning
    (`rheology.apply_stretch_thinning`) alike, thinning is thinning) erupts fresh crust in
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
    rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate_id, line_index))
    volcano_remaining[melting] = rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=int(melting.sum()))


def _ignite_early_rift_volcanoes(
    world: "World",  # noqa: F821
    plate_id: int,
    line_index: int,
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
    # event types don't share a random stream at the same (seed, step, plate, line).
    rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate_id, line_index, 1))
    volcano_remaining[ignite] = rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=int(ignite.sum()))


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


def _distance_to_mask_1d(mask: np.ndarray, width: int) -> np.ndarray:
    """Per-element distance (in node-index steps) to the nearest `True` in `mask`, within the
    1-D node order (used per line, so no wrap). Saturates at `width + 1` for anything `width`
    steps or farther away -- callers that only care about the ring within `width` steps of
    `mask` don't need an exact distance transform past its own edge. Backs the collision-uplift
    *reach* knob's near-field ring: a wider contested band -> the orogenic thickening spreads
    into a broader belt, the same "how far inland does a collision crumple crust" lever v1 had
    as COLLISION_RANGE_RAD, now with a smooth taper (COLLISION_NEAR_FIELD_INNER_FACTOR -> 0)
    across the ring rather than a flat rate -- which needs *how far into* the ring a node is,
    not just whether it's in the ring at all."""
    sentinel = width + 1
    dist = np.where(mask, 0, sentinel)
    if width <= 0 or not mask.any():
        return dist
    for shift in range(1, width + 1):
        dist[shift:] = np.minimum(dist[shift:], np.where(mask[:-shift], shift, sentinel))
        dist[:-shift] = np.minimum(dist[:-shift], np.where(mask[shift:], shift, sentinel))
    return dist


# The collision-uplift *reach* knob (World.collision_uplift_reach_multiplier) dilates the
# contested band feeding the orogenic thickening by this many physical km per unit of the
# knob (so reach 3x -> +1050 km each side of every contested stretch) -- linear in the knob,
# not "extra beyond 1.0": at the knob's own untuned value (1.0) this is already a real
# near-field belt (~350 km, a real orogen's crumple-zone width), not zero, so the near-field
# ring is part of the model's own baseline collision-uplift behaviour, not something that
# only exists once a user raises this knob above default (see GitHub issue #120, "Land
# fraction slowly declines" -- measured, this also modestly slows the land-fraction decline
# in its own right, since more of a collision's crust ends up thickened rather than left for
# erosion to plane down untouched).
#
# 2026-09-14 (GitHub issue #146, "Mountain ranges are too thin"): raised 200 -> 350 km. The
# real-world target is the topographic belt itself: the Himalaya proper (Main Frontal Thrust
# to the Indus-Tsangpo suture) runs
# ~150-350 km depending on strike segment, widening toward the Pakistan/Nanga Parbat and
# Arunachal syntaxes; the Andes run ~200-900 km along their length. 350 km sits at the wide
# end of the Himalaya range and mid-pack for the Andes, plus this ring is additive on top of
# the (much narrower) geometric contested band itself.
#
# Expressed in km, not a flat node count: `_distance_to_mask_1d` still operates on node
# indices (there is no cheaper way to widen a per-line band), but the index count converted to
# is divided by
# this step's *actual* line spacing (`spacing_rad`, which shrinks as `world.node_density`
# rises) so the belt's physical width stays ~350 km regardless of render/simulation
# resolution -- a flat node count would otherwise make mountains visibly narrower at higher
# node_density (confirmed: at density=4 a flat 2-node ring is only ~125 km, well under the
# real-world width it's meant to model).
COLLISION_NEAR_FIELD_REACH_KM_PER_UNIT = 350.0
# Near-field (dilated-but-not-contested) nodes thicken at a fraction of the contested rate
# that tapers (cell-centered linear ramp -- see the taper computation at its one call site)
# across the ring, from just under COLLISION_NEAR_FIELD_INNER_FACTOR right outside the
# contested band down to just above 0 at the ring's own outer edge (COLLISION_NEAR_FIELD_
# REACH_KM_PER_UNIT away) -- a collision belt's deformation fades outward from the suture, it
# doesn't step. 2026-09-16 (GitHub issue #146, "Mountain ranges are too thin", cause 4 of that
# investigation): this used to be a flat factor across the whole ring, then a hard drop to
# zero right past it -- a two-level "shelf" profile, not a falloff, which reads as an abrupt
# foothill line rather than a real orogen's gradual taper into its foreland. Doubling the old
# flat 0.4 into a ramp from ~0.8 down to ~0 keeps the ring's average thickening rate exactly
# the old flat value (a symmetric linear ramp always averages to its midpoint, so this knob's
# own existing tuning-knob tests, which only assert monotonicity, are unaffected) while making
# the taper itself continuous instead of a step.
COLLISION_NEAR_FIELD_INNER_FACTOR = 0.8


def _redistribute_accreted_column(
    persistent_fields: dict[str, np.ndarray],
    elevation: np.ndarray,
    rho_c: float,
    removed_hc: np.ndarray,
    removed_hm: np.ndarray,
    accrete_removed: np.ndarray,
    from_high: bool,
) -> None:
    """Conserve the crustal volume *and* attached mantle lithosphere of the continental-suture
    nodes just dropped from a line end (`removed_hc`/`removed_hm`, restricted to the
    `accrete_removed` subset) by thrusting both back onto the `SUTURE_ACCRETION_SPREAD_NODES`
    surviving nodes nearest that same end, each node's elevation bumped by the isostatic delta.
    Mutates `persistent_fields`' Hc/Hm arrays and `elevation` in place.

    No-op when nothing dropped was flagged for accretion -- a passive-margin retreat against
    an *oceanic* neighbour leaves `accrete_removed` all-False, and that column is genuinely
    subducted, not preserved. Node area is constant per node, so summed Hc/Hm *is* the conserved
    volume (see SUTURE_ACCRETION_SPREAD_NODES / SUTURE_ACCRETION_MAX_HC_M)."""
    if not np.any(accrete_removed):
        return
    add_hc = float(np.sum(removed_hc[accrete_removed]))
    add_hm = float(np.sum(removed_hm[accrete_removed]))
    hc = persistent_fields["crustal_thickness_m"]
    hm = persistent_fields["mantle_lithosphere_thickness_m"]
    n = len(hc)
    if n == 0 or add_hc <= 0.0:
        return
    k = min(SUTURE_ACCRETION_SPREAD_NODES, n)
    idx = np.arange(n - k, n) if from_high else np.arange(k)
    before = lithosphere.isostatic_elevation(hc[idx], hm[idx], rho_c)
    # Both reservoirs are thrust onto the survivors the same way -- the dropped nodes' own
    # summed Hc/Hm, spread evenly across the k survivors and added to what they already carry.
    # Previously Hm was only scaled by the survivor's own Hc growth ratio, which never
    # referenced the donor's actual removed_hm at all: the donor's mantle lithosphere was
    # silently discarded rather than conserved (GitHub issue #216's phase-budget instrumentation
    # flagged this while investigating the Hc/Hm decline). Each capped independently (issue
    # #161): a node that started thin can otherwise absorb most of a suture's whole volume in
    # one step.
    hc[idx] = np.minimum(hc[idx] + add_hc / k, SUTURE_ACCRETION_MAX_HC_M)
    hm[idx] = np.minimum(hm[idx] + add_hm / k, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)
    after = lithosphere.isostatic_elevation(hc[idx], hm[idx], rho_c)
    elevation[idx] = rheology.clip_elevation_bounds(elevation[idx] + (after - before))


class LithospherePlate(PlateWithLines):
    """A `PlateWithLines` whose per-node state is a lithospheric column (Hc/Hm) rather than
    an independently-set elevation -- see elevation_lines.py's own note on the two new
    OPTIONAL_FIELDS this relies on."""

    def crust_density(self) -> float:
        return lithosphere.crust_density(self.crust_type)

    def node_crust_density(self) -> np.ndarray:
        """Per-node counterpart of `crust_density()` -- see `lithosphere.node_crust_density`."""
        return lithosphere.node_crust_density(self.collect("crust_type_code"), self.crust_type)

    # -- Motion: torque.py's real implementation -----------------------------------------

    def shift(self, world: "World", years: float) -> float:  # noqa: F821 (World only for typing)
        # Debug-world-only escape hatch (see World.pinned_omegas): a scripted "Debugging
        # Worlds" scenario wants a plate to move exactly as configured every step, not
        # whatever the real torque balance (ridge-push/slab-pull/basal-drag against
        # world.mantle_centers) happens to settle it to -- empty for every ordinarily
        # generated world, so this is a no-op there.
        pinned = world.pinned_omegas.get(self.plate_id)
        if pinned is not None:
            old_points, _ = self.all_points_and_elevation()
            if len(old_points) == 0:
                return 0.0
            return torque.apply_omega_and_rotate(self, old_points, np.asarray(pinned, dtype=float), years)
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

        # GitHub issue #177 direction 1 -- see OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER's own
        # comment for the full rationale. This step's real arc-magmatic creation across the
        # whole plate, computed once here (mirrors the per-line calculation below) into a
        # shared budget the per-line loop's oceanic-override retreats spend down in place.
        #
        # GitHub issue #216 -- see OCEANIC_SELF_RETREAT_BUDGET_MULTIPLIER's own comment. An
        # oceanic self-plate has no arc-magmatic creation of its own (arc_band_all is always
        # empty for it), so it gets the same budget array filled from a different source: this
        # step's real decompression-melting eruption across the whole plate, previewed the same
        # way (mirrors the per-line divergent-deformation call below) -- only the *melted*
        # subset represents genuinely new mass; ordinary sub-critical thinning is a loss already
        # counted elsewhere (divergent_deformation), not creation, and including it here would
        # make the budget negative almost every step and block ordinary subduction outright.
        oceanic_override_retreat_budget_hc = np.zeros(1)
        if self.crust_type == "continental" and np.any(arc_band_all):
            hm_all = self.collect("mantle_lithosphere_thickness_m")
            grown_hc, _ = rheology.apply_arc_magmatic_thickening(
                hc_all[arc_band_all], hm_all[arc_band_all], closing_rate_all[arc_band_all],
                years_myr, arc_intensity_all[arc_band_all],
            )
            oceanic_override_retreat_budget_hc[0] = (
                OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER * float(np.sum(grown_hc - hc_all[arc_band_all]))
            )
        elif self.crust_type == "oceanic" and np.any(divergent_all):
            hc_div_all = self.collect("crustal_thickness_m")[divergent_all]
            hm_div_all = self.collect("mantle_lithosphere_thickness_m")[divergent_all]
            _, _, melting_all = rheology.apply_divergent_deformation(
                hc_div_all, hm_div_all, closing_rate_all[divergent_all], years_myr,
            )
            if np.any(melting_all):
                oceanic_override_retreat_budget_hc[0] = (
                    OCEANIC_SELF_RETREAT_BUDGET_MULTIPLIER
                    * float(np.sum(lithosphere.REFERENCE_HC_OCEANIC_M - hc_div_all[melting_all]))
                )

        # Collision-uplift tuning knobs (the "Controls" window, 1.0 == untuned -- see World).
        # `orogen_amount` scales the plastic thickening rate at contested nodes; `orogen_reach`
        # widens (>1) or narrows (<1) the belt it acts on -- see _distance_to_mask_1d /
        # COLLISION_NEAR_FIELD_*. At 1.0, `orogen_amount` leaves apply_convergent_deformation's
        # contested-band strength at exactly 1.0, same as ever -- but `orogen_reach` no longer
        # means "no near-field ring below/at 1.0, only above": the ring is linear in the knob
        # from 0 (see COLLISION_NEAR_FIELD_REACH_KM_PER_UNIT's own comment for why the model's
        # own baseline collision belt already carries one at the knob's untuned value).
        orogen_amount = world.collision_uplift_multiplier
        orogen_reach = world.collision_uplift_reach_multiplier
        # How much of this plate's own active margin is currently jammed in overlap right now
        # -- normalized against the near-boundary band, not the whole plate (a huge plate's
        # boundary is a small fraction of its own node count, which would dilute this to
        # near-zero for exactly the large-plate case that matters). A deeper/wider overlap
        # should crumple faster than a light graze, on top of the existing distance-decay
        # shape within the belt -- see OVERLAP_UPLIFT_SEVERITY_GAIN.
        band_all = inputs.dist_to_neighbor <= reach_rad
        overlap_severity = float(np.count_nonzero(contested_all)) / max(1, int(np.count_nonzero(band_all)))
        orogen_contested_strength = orogen_amount * min(orogen_reach, 1.0) * (1.0 + OVERLAP_UPLIFT_SEVERITY_GAIN * overlap_severity)
        orogen_dilation_nodes = (
            round(orogen_reach * COLLISION_NEAR_FIELD_REACH_KM_PER_UNIT / (spacing_rad * PLANET_RADIUS_KM))
            if orogen_reach > 0.0 and self.crust_type == "continental"
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
            # GitHub issue #216 Hc/Hm budget checkpoints -- see phase_budget.py. `codes0` is
            # this line's crust_type_code, unchanged until the decompression-melting checkpoint
            # below, so every intermediate checkpoint below reuses it for both before/after.
            codes0 = line.crust_type_code
            checkpoint_hc, checkpoint_hm = (hc.copy(), hm.copy()) if world.debug_diagnostics else (None, None)
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
            # `orogen_contested_strength`, plus (reach knob > 1) a dilated near-field ring
            # that tapers linearly from COLLISION_NEAR_FIELD_INNER_FACTOR right outside the
            # contested band down to 0 at the ring's own outer edge -- see that constant's own
            # comment (issue #146). `orogen_strength` is the per-node multiplier handed to
            # apply_convergent_deformation; > 0 exactly on the nodes that thicken.
            # `apply_convergent_deformation` still gates on each node's own closing rate
            # (below yield / not actually closing -> zero strain), so a node that is
            # `convergent` only via the `contested` deep-overlap fold and is no longer
            # actively closing simply thickens at zero.
            near_field_dist = _distance_to_mask_1d(convergent, orogen_dilation_nodes) if orogen_dilation_nodes > 0 else None
            near_field = (
                (near_field_dist <= orogen_dilation_nodes) & ~convergent & ~divergent
                if near_field_dist is not None
                else np.zeros(n, dtype=bool)
            )
            orogen_strength = np.where(convergent, orogen_contested_strength, 0.0)
            if np.any(near_field):
                # Cell-centered, not edge-to-edge: node `d` (1-indexed) is treated as sitting
                # at the middle of its own step, so the outermost ring node still gets a small
                # nonzero share (COLLISION_NEAR_FIELD_INNER_FACTOR / (2*orogen_dilation_nodes))
                # instead of tapering all the way to exactly 0 right at the last real node --
                # an edge-to-edge ramp would otherwise re-introduce a (smaller) hard step at a
                # small reach (e.g. a 1-2 node ring at coarse node_density), the same shelf-vs-
                # falloff problem this taper exists to fix. This also keeps the ring's mean
                # strength at exactly COLLISION_NEAR_FIELD_INNER_FACTOR / 2 regardless of
                # orogen_dilation_nodes (a symmetric linear ramp always averages to its
                # midpoint), matching the old flat factor's total contribution.
                taper = np.clip(1.0 - (near_field_dist - 0.5) / orogen_dilation_nodes, 0.0, 1.0)
                orogen_strength[near_field] = orogen_amount * COLLISION_NEAR_FIELD_INNER_FACTOR * taper[near_field]
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

                # Lateral magma export (GitHub issue #205, follow-up to #120's "Land fraction
                # slowly declines"): divert a fraction of the core convergent band's own strain
                # increment to a mobile magma parcel instead of thickening the node in place --
                # see rheology.magma_export_strength_and_volume's own docstring for the full
                # mechanism/reasoning. Masked to `convergent` only (never the near-field ring
                # below): the ring's own melt supply already comes from the delamination-
                # overflow path a few lines down, and skimming it here too would starve that
                # supply a second, independent way (see that function's own docstring on why
                # near-ceiling nodes are exempt for the same reason).
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
                        volumes_m3 = export_hc[exporting] * lithosphere.node_area_m2(spacing_rad)
                        origins = own_points[sl][node_idx]
                        world.pending_magma_parcels.extend(
                            magma_transport.MagmaParcel(
                                origin_xyz=origins[i], volume_m3=float(volumes_m3[i]), step_generated=world.steps_taken
                            )
                            for i in range(len(node_idx))
                        )

                new_hc, new_hm, overflow_hc = rheology.apply_convergent_deformation(
                    hc[thicken], hm[thicken], closing_rate[thicken], years_myr,
                    fault_factor[thicken], strength=used_strength,
                )
                hc[thicken] = new_hc
                hm[thicken] = new_hm

                # Hc that hit MAX_CRUSTAL_THICKNESS_M this step didn't just vanish (issue
                # #161) -- but it also doesn't reappear whole and instant on the foreland
                # either (issue #145's reopened investigation: that turned out to over-
                # thicken/elevate the majority of a run's continental land within tens of
                # Myr). It delaminates, partially remelts, and the buoyant melt fraction
                # intrudes the near-field ring at a bounded rate -- see
                # rheology.apply_delamination_melt_intrusion's own docstring for the full
                # reasoning. Only the core band's overflow is conserved (partially) this way;
                # the near-field ring's own overflow (rarer -- it thickens at a faded rate
                # already) has nowhere further out to spread to on this pass and delaminates
                # in full, same as suture accretion's own overflow past its cap. No-op when
                # there's no near-field ring to receive it (reach knob at 0, or an oceanic
                # plate, which never gets one).
                overflow_total = float(np.sum(overflow_hc[convergent[thicken]]))
                if overflow_total > 0.0 and np.any(near_field):
                    hc[near_field] = rheology.apply_delamination_melt_intrusion(hc[near_field], overflow_total, years_myr)

            if world.debug_diagnostics:
                phase_budget.record(world, self, "convergent_deformation", checkpoint_hc, checkpoint_hm, codes0, hc, hm, codes0)
                checkpoint_hc, checkpoint_hm = hc.copy(), hm.copy()

            # Continental arc magmatism: an oceanic slab subducting under this margin fluxes
            # the mantle wedge and underplates juvenile crust across the whole arc band --
            # extra Hc (added from the mantle, not conserved), the crust-building half of
            # "subduction under a continent makes more continent" (GitHub issue #120, "Land fraction
            # slowly declines"). Separate from the contested shortening above: the band is far
            # wider than the contact line. Bounded long-term by the CONTINENTAL_AREA_BUDGET_MULT
            # volume gate.
            if np.any(arc_band):
                hc[arc_band], hm[arc_band] = rheology.apply_arc_magmatic_thickening(
                    hc[arc_band], hm[arc_band], closing_rate[arc_band], years_myr, arc_intensity[arc_band]
                )

            if world.debug_diagnostics:
                phase_budget.record(world, self, "arc_magmatism", checkpoint_hc, checkpoint_hm, codes0, hc, hm, codes0)
                checkpoint_hc, checkpoint_hm = hc.copy(), hm.copy()

            prior_hc = hc.copy()
            melting = np.zeros(n, dtype=bool)
            newly_below_rift_onset = np.zeros(n, dtype=bool)
            if np.any(divergent):
                new_hc, new_hm, melt = rheology.apply_divergent_deformation(hc[divergent], hm[divergent], closing_rate[divergent], years_myr)
                # "fault" mode: scale the thinning delta by fault proximity (all-ones
                # otherwise). Melt (decompression volcanism) still fires on the geometric
                # rift threshold -- it's a discrete event, not a rate.
                infl = fault_influence[divergent]
                hc[divergent] = hc[divergent] + infl * (new_hc - hc[divergent])
                hm[divergent] = hm[divergent] + infl * (new_hm - hm[divergent])
                melting[divergent] = melt

                # Rift magmatic underplating (see rheology.apply_rift_magmatic_thickening): a
                # partial Hc offset for nodes that thinned past RIFT_VOLCANISM_ONSET_HC_M but
                # didn't melt all the way through this step -- nodes that did melt already got
                # the full reference-column reset below and don't need this on top of it.
                magmatic_band = divergent & ~melting
                if np.any(magmatic_band):
                    new_hc_mag, new_hm_mag = rheology.apply_rift_magmatic_thickening(
                        hc[magmatic_band], hm[magmatic_band], closing_rate[magmatic_band], years_myr
                    )
                    hc[magmatic_band] = new_hc_mag
                    hm[magmatic_band] = new_hm_mag
                    newly_below_rift_onset = (
                        magmatic_band & (prior_hc >= rheology.RIFT_VOLCANISM_ONSET_HC_M) & (hc < rheology.RIFT_VOLCANISM_ONSET_HC_M)
                    )

            if world.debug_diagnostics:
                phase_budget.record(world, self, "divergent_deformation", checkpoint_hc, checkpoint_hm, codes0, hc, hm, codes0)
                checkpoint_hc, checkpoint_hm = hc.copy(), hm.copy()

            prior_age = line.divergent_age_myr
            new_age = np.where(divergent, prior_age + years_myr, 0.0)
            if self.crust_type == "oceanic":
                hm = rheology.relax_young_oceanic_mantle_lithosphere(hm, new_age, years_myr)
                if world.debug_diagnostics:
                    phase_budget.record(world, self, "oceanic_cooling_relaxation", checkpoint_hc, checkpoint_hm, codes0, hc, hm, codes0)
                    checkpoint_hc, checkpoint_hm = hc.copy(), hm.copy()

            is_volcano = line.is_volcano.copy()
            volcano_remaining = line.volcano_active_years_remaining.copy()
            crust_type_code = line.crust_type_code.copy()
            # Decompression melting (spec 2.3): a rift that just thinned past the critical
            # threshold erupts fresh crust in place -- same one-guaranteed-eruption convention
            # v1's stretch-volcano growth used. The erupted material's type depends on where it
            # surfaces: still standing above sea level (this node's *pre-melt* elevation, i.e.
            # this line's own elevation before today's deform() pass touched it) is
            # continental-type magmatism -- real continental rifts stay bimodal-volcanic land
            # for a long stretch before a true ocean opens (the East African Rift, well before
            # the Red Sea stage) -- while a node already at or below sea level (a drowned
            # margin, or an ordinary oceanic ridge) erupts ordinary mid-ocean-ridge oceanic
            # crust. See docs/simulation-model.md's "Magma-typed decompression melting".
            _erupt_melted_nodes(world, self.plate_id, line_index, hc, hm, crust_type_code, is_volcano, volcano_remaining, melting, line.elevation)
            phase_budget.record(world, self, "decompression_melting", checkpoint_hc, checkpoint_hm, codes0, hc, hm, crust_type_code)
            _ignite_early_rift_volcanoes(world, self.plate_id, line_index, is_volcano, volcano_remaining, newly_below_rift_onset)

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

            # Issue #189 follow-up (see UNBACKED_RELIEF_DECAY_PER_MYR above): relax any
            # *existing* positive debt -- elevation this line already carries in excess of what
            # `elevation_before` says its own Hc/Hm column supports -- toward zero, before this
            # step's own fresh transform_uplift (still a bare delta by design) potentially
            # adds more. `line.crustal_thickness_m` (the column at the *start* of this step,
            # matching what `elevation_before` was computed from), not `hc` (already mutated
            # by the convergent/divergent passes above) -- and v1 lines with no Hc tracking at
            # all (all-zero) are left alone, same has_column gating
            # lithosphere.back_elevation_gain uses.
            existing_debt = np.where(line.crustal_thickness_m > 0.0, np.clip(line.elevation - elevation_before, 0.0, None), 0.0)
            debt_relief = existing_debt * (1.0 - np.exp(-UNBACKED_RELIEF_DECAY_PER_MYR * years_myr))

            new_elevation = rheology.clip_elevation_bounds(
                line.elevation - debt_relief + (elevation_after - elevation_before) + transform_uplift
            )
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
                melt_rho_c = lithosphere.node_crust_density(crust_type_code[melting], self.crust_type)
                new_elevation[melting] = rheology.clip_elevation_bounds(
                    lithosphere.isostatic_elevation(hc[melting], hm[melting], melt_rho_c)
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
            # Gated on transform_uplift itself, not just band membership -- issue #189
            # follow-up's debt-decay term can move a transform-band node's elevation on its own
            # (fault_influence == 0 in "fault" mode zeroes transform_uplift there, but leftover
            # debt still decays), which would otherwise mislabel a pure decay move as an active
            # transform pressure ridge.
            reason[(transform_uplift > 0.0) & moved] = ELEV_CHANGE_TRANSFORM
            reason[melting] = ELEV_CHANGE_VOLCANO

            updated_line = line.replace(
                elevation=new_elevation,
                crustal_thickness_m=hc,
                mantle_lithosphere_thickness_m=hm,
                divergent_age_myr=new_age,
                is_volcano=is_volcano,
                volcano_active_years_remaining=volcano_remaining,
                elev_change_reason=reason,
                crust_type_code=crust_type_code,
            )
            grown_lines = self._grow_or_shrink_line_for_deform(
                updated_line,
                inputs.dist_to_neighbor[sl],
                inputs.direction_to_neighbor[sl],
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
                oceanic_override_retreat_budget_hc,
                suppress_growth,
                arc_end_low,
                arc_end_high,
            )
            if world.debug_diagnostics:
                # Endpoint stretch-thinning, end growth, end/interior retreat, and accreted-
                # column redistribution all happen inside this one call (GitHub issue #216
                # items 1-3). Each sub-mechanism now has its own before/after snapshot taken
                # deep inside _grow_or_shrink_line_for_deform (line_end_stretch, line_end_
                # arc_grow, line_end_retreat, line_end_accretion, line_interior_carve) -- this
                # "line_growth_shrink" record is kept as the net line-boundary total, a
                # cross-check that the five sub-phases' combined delta matches the whole call's
                # actual effect (they should sum to it exactly, since together they cover every
                # mutation this function makes).
                grown_hc = np.concatenate([gl.crustal_thickness_m for gl in grown_lines]) if grown_lines else np.array([])
                grown_hm = np.concatenate([gl.mantle_lithosphere_thickness_m for gl in grown_lines]) if grown_lines else np.array([])
                grown_codes = (
                    np.concatenate([gl.crust_type_code for gl in grown_lines]) if grown_lines else np.array([], dtype=updated_line.crust_type_code.dtype)
                )
                phase_budget.record(
                    world, self, "line_growth_shrink",
                    updated_line.crustal_thickness_m, updated_line.mantle_lithosphere_thickness_m, updated_line.crust_type_code,
                    grown_hc, grown_hm, grown_codes,
                )
            new_lines.extend(gl for gl in grown_lines if len(gl) > 0)

        if self.crust_type == "continental":
            if world.debug_diagnostics:
                before_hc = np.concatenate([l.crustal_thickness_m for l in new_lines]) if new_lines else np.array([])
                before_hm = np.concatenate([l.mantle_lithosphere_thickness_m for l in new_lines]) if new_lines else np.array([])
                before_codes = np.concatenate([l.crust_type_code for l in new_lines]) if new_lines else np.array([])
            new_lines = self._retreat_contested_leading_rows(new_lines, contested_all, years)
            if world.debug_diagnostics:
                after_hc = np.concatenate([l.crustal_thickness_m for l in new_lines]) if new_lines else np.array([])
                after_hm = np.concatenate([l.mantle_lithosphere_thickness_m for l in new_lines]) if new_lines else np.array([])
                after_codes = np.concatenate([l.crust_type_code for l in new_lines]) if new_lines else np.array([])
                phase_budget.record(world, self, "contested_leading_row_retreat", before_hc, before_hm, before_codes, after_hc, after_hm, after_codes)

        self.set_lines(new_lines)
        if not suppress_growth:
            # Built once and shared between the two calls below (both need "distance to the
            # nearest neighbour node", nothing else about which neighbour) -- previously each
            # rebuilt its own identical cKDTree over the same concatenated neighbour-node
            # cloud from scratch (docs/profiling.md; neither of these two calls touches any
            # *neighbour's* own nodes, only self's, so the cloud can't have changed between
            # them). Not the torque.gather_boundary_force_inputs per-neighbour-tree pattern --
            # that trades one big batched query for several smaller ones, a good trade when
            # the query side is tiny (as in _claim_adjacent_territory's own per-row probe) but
            # a bad one for _fill_corner_notch_frontier's single large batched query, where a
            # single combined tree stays the cheaper query shape.
            neighbour_points = [p.all_points_and_elevation()[0] for p in neighbours if p.node_count() > 0]
            neighbour_tree = cKDTree(np.concatenate(neighbour_points, axis=0)) if neighbour_points else None
            # GitHub issue #216 items: adjacent-row claiming and corner-notch fill both grow
            # this plate's own node population at its own (plate-level, not per-line) expense
            # -- snapshotted via `self.collect(...)` around each call rather than per-line
            # since neither operates line-by-line.
            before_claim = phase_budget.snapshot(self) if world.debug_diagnostics else None
            self._claim_adjacent_territory(world, neighbours, spacing_rad, neighbour_tree=neighbour_tree)
            if world.debug_diagnostics:
                after_claim = phase_budget.snapshot(self)
                phase_budget.record(world, self, "adjacent_row_claim", *before_claim, *after_claim)
            self._fill_corner_notch_frontier(world, neighbours, spacing_rad, years, neighbour_tree=neighbour_tree)
            if world.debug_diagnostics:
                after_notch = phase_budget.snapshot(self)
                phase_budget.record(world, self, "corner_notch_fill", *after_claim, *after_notch)

        for line_index, line in enumerate(self.lines):
            if needs_regularizing(line, spacing_rad):
                if world.debug_diagnostics:
                    before = (line.crustal_thickness_m, line.mantle_lithosphere_thickness_m, line.crust_type_code)
                    regularized = regularize_line(line, spacing_rad)
                    phase_budget.record(
                        world, self, "line_regularization",
                        *before, regularized.crustal_thickness_m, regularized.mantle_lithosphere_thickness_m, regularized.crust_type_code,
                    )
                    self.replace_line(line_index, regularized)
                else:
                    self.replace_line(line_index, regularize_line(line, spacing_rad))

    def _count_open_prefix(self, theta_candidates: np.ndarray, phi: float, neighbours: list) -> int:
        if len(theta_candidates) == 0 or not neighbours:
            return len(theta_candidates)
        world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.full_like(theta_candidates, phi), theta_candidates))
        contested = _contested_by_any(world_pts, neighbours)
        first_contested = np.argmax(contested) if np.any(contested) else len(contested)
        return int(first_contested)

    def _separation_components(self, world: "World", phi: float, theta: float, direction_world: np.ndarray) -> tuple[float, float]:  # noqa: F821
        """(sep_theta, sep_phi): this node's own local rift-separation direction, decomposed
        along its row's own (theta, phi) tangent basis -- see rheology.stretch_components for
        what the caller does with it. "fault" mode prefers the nearest active fault's own
        tangent (faults.fault_tangent_components, swapped -- a fault opens *across* its own
        strike); every other mode, or a "fault"-mode plate with no active fault there yet,
        falls back to this node's own `direction_to_neighbor` projected into the same local
        basis (geometry.local_separation_components) -- the separation direction is normal to
        a rift's own trend, so this stands in for "the line that would pass through most of
        the empty space" without needing to fit one."""
        if getattr(world, "fault_deformation_mode", "fault") == "fault":
            from . import faults

            tangent = faults.fault_tangent_components(world, self, phi, theta)
            if tangent is not None:
                return tangent
        if float(np.linalg.norm(direction_world)) < 1e-9:
            # No neighbour close enough to have a meaningful direction toward at all (a
            # genuinely isolated plate/stretch of open water, `direction_world` left at zero
            # by torque.gather_boundary_force_inputs) -- fall back to this row's own theta
            # direction, matching ordinary end growth's long-standing behaviour when there is
            # nothing nearby to orient a stretch against.
            return 1.0, 0.0
        return geometry.local_separation_components(self.frame, phi, theta, direction_world)

    def _grow_or_shrink_line_for_deform(
        self,
        line: ElevationLine,
        dist: np.ndarray,
        direction: np.ndarray,
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
        oceanic_override_retreat_budget_hc: np.ndarray,
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
        margin against an oceanic slab) and is further rate-limited against
        `oceanic_override_retreat_budget_hc` (see `_budget_limited_removal` /
        OCEANIC_OVERRIDE_RETREAT_BUDGET_MULTIPLIER for a continental self-plate's oceanic-
        override case, issue #177 direction 1; OCEANIC_SELF_RETREAT_BUDGET_MULTIPLIER for an
        oceanic self-plate's own ordinary subduction, both end and interior-carve, issue #216 --
        only one of the two ever applies to a given plate, so both share this one array). That
        budget is shared and mutable across this plate's whole deform() pass (every line, both
        ends, and the interior carve), so it must be threaded through from the caller rather
        than recomputed here.

        `direction` (world-frame, this plate's own `torque.BoundaryForceInputs.
        direction_to_neighbor`, one per node) feeds the end-growth branch's own rift-stretch
        decomposition (`rheology.stretch_components` via `self._separation_components`) --
        see that branch's own comment for why growth now stretches an existing end node
        rather than always appending fresh ones."""
        theta = line.theta.copy()
        elevation = line.elevation.copy()
        contested = contested.copy()
        shrinkable = shrinkable.copy()
        accrete = accrete.copy()
        dist = dist.copy()
        direction = direction.copy()
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
            # Budget-limited whenever this retreat isn't conserved by accretion -- continental
            # vs. oceanic-override (issue #177) or, now, an oceanic self-plate's own ordinary
            # subduction (issue #216) -- `accrete[-1]` is always False for the latter (oceanic
            # self-plates never accrete, see `accrete_all`), so dropping the old
            # `self.crust_type == "continental"` restriction here just lets that case reach the
            # same budget machinery, spending from whichever of the two budgets `deform()`
            # actually populated for this plate's own crust_type.
            if n_remove > 0 and not accrete[-1]:
                n_remove = _budget_limited_removal(
                    persistent_fields["crustal_thickness_m"], n_remove, oceanic_override_retreat_budget_hc, from_high=True
                )
            if n_remove > 0:
                removed_hc = persistent_fields["crustal_thickness_m"][-n_remove:].copy()
                removed_hm = persistent_fields["mantle_lithosphere_thickness_m"][-n_remove:].copy()
                removed_codes = persistent_fields["crust_type_code"][-n_remove:].copy()
                accrete_removed = accrete[-n_remove:].copy()
                removed_world = geometry.to_world(self.frame, geometry.local_xyz(np.full(n_remove, line.phi), theta[-n_remove:]))
                world.record_removed_points(removed_world, self.plate_id)
                theta, elevation = theta[:-n_remove], elevation[:-n_remove]
                contested, shrinkable, accrete, dist = contested[:-n_remove], shrinkable[:-n_remove], accrete[:-n_remove], dist[:-n_remove]
                direction = direction[:-n_remove]
                persistent_fields = {name: values[:-n_remove] for name, values in persistent_fields.items()}
                if world.debug_diagnostics:
                    phase_budget.record(
                        world, self, "line_end_retreat",
                        removed_hc, removed_hm, removed_codes,
                        np.array([]), np.array([]), np.array([], dtype=removed_codes.dtype),
                    )
                accrete_pre_hc = accrete_pre_hm = accrete_pre_codes = None
                if world.debug_diagnostics and np.any(accrete_removed):
                    accrete_pre_hc = persistent_fields["crustal_thickness_m"].copy()
                    accrete_pre_hm = persistent_fields["mantle_lithosphere_thickness_m"].copy()
                    accrete_pre_codes = persistent_fields["crust_type_code"].copy()
                _redistribute_accreted_column(persistent_fields, elevation, rho_c, removed_hc, removed_hm, accrete_removed, from_high=True)
                if accrete_pre_hc is not None:
                    phase_budget.record(
                        world, self, "line_end_accretion",
                        accrete_pre_hc, accrete_pre_hm, accrete_pre_codes,
                        persistent_fields["crustal_thickness_m"], persistent_fields["mantle_lithosphere_thickness_m"], persistent_fields["crust_type_code"],
                    )

        if len(theta) == 0:
            return [ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)]

        if shrinkable[0]:
            n_remove = min(contested_run_from_end(shrinkable, from_high=False), n_distance_cap, max_extend_nodes, len(theta) - 1)
            # See the mirrored high-end block above for why the crust-type restriction is gone.
            if n_remove > 0 and not accrete[0]:
                n_remove = _budget_limited_removal(
                    persistent_fields["crustal_thickness_m"], n_remove, oceanic_override_retreat_budget_hc, from_high=False
                )
            if n_remove > 0:
                removed_hc = persistent_fields["crustal_thickness_m"][:n_remove].copy()
                removed_hm = persistent_fields["mantle_lithosphere_thickness_m"][:n_remove].copy()
                removed_codes = persistent_fields["crust_type_code"][:n_remove].copy()
                accrete_removed = accrete[:n_remove].copy()
                removed_world = geometry.to_world(self.frame, geometry.local_xyz(np.full(n_remove, line.phi), theta[:n_remove]))
                world.record_removed_points(removed_world, self.plate_id)
                theta, elevation = theta[n_remove:], elevation[n_remove:]
                contested, shrinkable, accrete, dist = contested[n_remove:], shrinkable[n_remove:], accrete[n_remove:], dist[n_remove:]
                direction = direction[n_remove:]
                persistent_fields = {name: values[n_remove:] for name, values in persistent_fields.items()}
                if world.debug_diagnostics:
                    phase_budget.record(
                        world, self, "line_end_retreat",
                        removed_hc, removed_hm, removed_codes,
                        np.array([]), np.array([]), np.array([], dtype=removed_codes.dtype),
                    )
                accrete_pre_hc = accrete_pre_hm = accrete_pre_codes = None
                if world.debug_diagnostics and np.any(accrete_removed):
                    accrete_pre_hc = persistent_fields["crustal_thickness_m"].copy()
                    accrete_pre_hm = persistent_fields["mantle_lithosphere_thickness_m"].copy()
                    accrete_pre_codes = persistent_fields["crust_type_code"].copy()
                _redistribute_accreted_column(persistent_fields, elevation, rho_c, removed_hc, removed_hm, accrete_removed, from_high=False)
                if accrete_pre_hc is not None:
                    phase_budget.record(
                        world, self, "line_end_accretion",
                        accrete_pre_hc, accrete_pre_hm, accrete_pre_codes,
                        persistent_fields["crustal_thickness_m"], persistent_fields["mantle_lithosphere_thickness_m"], persistent_fields["crust_type_code"],
                    )

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
                # GitHub issue #216: this is oceanic self-plate subduction exactly like the end
                # retreats above, just mid-row -- cap it by the same plate-wide Hc budget
                # (oceanic_override_retreat_budget_hc) those spend from, rather than leaving it
                # bounded only by node count. Counted in from the run's own start, mirroring
                # `_budget_limited_removal`'s own prefix-sum approach.
                run_hc = persistent_fields["crustal_thickness_m"][start : start + take]
                cum_hc = np.cumsum(run_hc)
                remaining_budget = max(float(oceanic_override_retreat_budget_hc[0]), 0.0)
                within_budget = cum_hc <= remaining_budget
                take = take if within_budget.all() else int(np.argmax(~within_budget))
                if take <= 0:
                    continue
                oceanic_override_retreat_budget_hc[0] -= float(cum_hc[take - 1])
                keep[start : start + take] = False
                budget -= take
            if not keep.all():
                removed_world = geometry.to_world(self.frame, geometry.local_xyz(np.full((~keep).sum(), line.phi), theta[~keep]))
                world.record_removed_points(removed_world, self.plate_id)
                if world.debug_diagnostics:
                    phase_budget.record(
                        world, self, "line_interior_carve",
                        persistent_fields["crustal_thickness_m"][~keep],
                        persistent_fields["mantle_lithosphere_thickness_m"][~keep],
                        persistent_fields["crust_type_code"][~keep],
                        np.array([]), np.array([]), np.array([], dtype=persistent_fields["crust_type_code"].dtype),
                    )
                theta, elevation = theta[keep], elevation[keep]
                contested, shrinkable, accrete, dist = contested[keep], shrinkable[keep], accrete[keep], dist[keep]
                direction = direction[keep]
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
                elif name == "node_created_years":
                    # These are genuinely brand-new arc-margin nodes (see the arc_end_high/
                    # arc_end_low append branches below) -- a zero-fill would misread as
                    # "created at year 0" rather than the -1.0 "unknown" sentinel, so stamp the
                    # real creation time explicitly, same as _seed_and_erupt_new_nodes does for
                    # its own two callers.
                    fill = np.full(n_new, world.elapsed_years, dtype=values.dtype)
                else:
                    fill = np.zeros(n_new, dtype=values.dtype)
                out[name] = fill
            return out

        def _stretch_end(index: int, sign: float, dist_end: float, direction_end: np.ndarray) -> None:
            """Close the ordinary (non-arc) share of a gap by stretching this row's own end
            node outward rather than appending brand-new full-thickness nodes -- see
            rheology.stretch_components/apply_stretch_thinning and this plate's own
            _separation_components. Only the theta-attributable share of the gap (how much of
            it runs along this row's own theta axis, vs. across to a neighbouring row -- see
            _claim_adjacent_territory for that other share) is closed here; capped per step at
            `(IRREGULARITY_TOLERANCE - 1) * dtheta`, *not* `max_extend_nodes`/`ring_room()`'s
            much larger per-step allowance (confirmed as a real bug: stretching the endpoint by
            several node-spacings in one jump leaves a gap past `CONTIGUOUS_RUN_GAP_MULT` to
            its inward neighbour, so `split_into_contiguous_runs` -- called on every result
            below -- reads the stretched node as a disconnected one-node fragment rather than
            part of a row that merely needs re-densifying). Staying inside
            `IRREGULARITY_TOLERANCE` instead means `needs_regularizing`/`regularize_line`
            (called every deform() pass, see below) always see this as "irregular, resample
            it," never "disconnected, split it" -- so a wide gap now closes the same way it did
            before this change: gradually, one bounded step at a time, just via repeated
            stretch-and-resample instead of repeated append. `regularize_line`'s own resample
            carries the now-thinned Hc/Hm through via interpolation exactly like every other
            persistent field already does -- so "no new points are created" here still ends up
            at target density, just one pass later, with the thinning already baked in."""
            # Clamp to a finite gap estimate *before* decomposing it -- `dist_end` is `inf`
            # whenever nothing is within `reach_rad` at all (an isolated end, or simply no
            # neighbour close enough yet), and `inf * 0.0` (a separation direction with an
            # exactly-zero component along one axis) is `nan`, not `0.0` -- the same
            # gap_estimate clamp the old node-appending growth always applied before this
            # change (see the arc branch just above, `gap_estimate = min(dist_end, ...)`).
            gap_estimate = min(dist_end, (n_distance_cap + 1) * spacing_rad)
            sep_theta, sep_phi = self._separation_components(world, line.phi, float(theta[index]), direction_end)
            theta_gap, _ = rheology.stretch_components(sep_theta, sep_phi, gap_estimate)
            theta_gap = min(theta_gap, (IRREGULARITY_TOLERANCE - 1.0) * dtheta, max(ring_room(), 0) * dtheta)
            if theta_gap <= 0.0:
                return
            candidate = np.array([theta[index] + sign * theta_gap])
            if self._count_open_prefix(candidate, line.phi, neighbours) == 0:
                return
            prior_elevation = float(elevation[index])
            if world.debug_diagnostics:
                stretch_before_hc = persistent_fields["crustal_thickness_m"][index, None].copy()
                stretch_before_hm = persistent_fields["mantle_lithosphere_thickness_m"][index, None].copy()
                stretch_before_codes = persistent_fields["crust_type_code"][index, None].copy()
            new_hc, new_hm, melt = rheology.apply_stretch_thinning(
                persistent_fields["crustal_thickness_m"][index, None],
                persistent_fields["mantle_lithosphere_thickness_m"][index, None],
                np.array([dtheta]),
                np.array([theta_gap]),
            )
            persistent_fields["crustal_thickness_m"][index] = new_hc[0]
            persistent_fields["mantle_lithosphere_thickness_m"][index] = new_hm[0]
            _erupt_melted_nodes(
                world,
                self.plate_id,
                line_index,
                persistent_fields["crustal_thickness_m"][index, None],
                persistent_fields["mantle_lithosphere_thickness_m"][index, None],
                persistent_fields["crust_type_code"][index, None],
                persistent_fields["is_volcano"][index, None],
                persistent_fields["volcano_active_years_remaining"][index, None],
                melt,
                np.array([prior_elevation]),
            )
            theta[index] = candidate[0]
            node_rho_c = lithosphere.node_crust_density(persistent_fields["crust_type_code"][index, None], self.crust_type)
            elevation[index] = lithosphere.isostatic_elevation(
                persistent_fields["crustal_thickness_m"][index, None], persistent_fields["mantle_lithosphere_thickness_m"][index, None], node_rho_c
            )[0]
            persistent_fields["elev_change_reason"][index] = ELEV_CHANGE_VOLCANO if melt[0] else ELEV_CHANGE_RIFT
            if world.debug_diagnostics:
                phase_budget.record(
                    world, self, "line_end_stretch",
                    stretch_before_hc, stretch_before_hm, stretch_before_codes,
                    persistent_fields["crustal_thickness_m"][index, None],
                    persistent_fields["mantle_lithosphere_thickness_m"][index, None],
                    persistent_fields["crust_type_code"][index, None],
                )

        # `suppress_growth` (the continental area-budget gate) only applies to the arc-seed
        # append branch below, which conjures brand-new full-thickness nodes for free. The
        # `_stretch_end` branch is mass-conserving (rheology.apply_stretch_thinning thins the
        # existing column by exact footprint conservation) until it crosses
        # RIFT_CRITICAL_THICKNESS_M, at which point it erupts via the same decompression-
        # melting/magma-upwelling path (_erupt_melted_nodes) ordinary divergent thinning
        # already uses -- not an unpaid-for area grab, so it shouldn't need this gate at all.
        if not contested[-1] and dist[-1] > extend_threshold_rad and ring_room() > 0:
            if arc_end_high:
                if not suppress_growth:
                    gap_estimate = min(dist[-1], (n_distance_cap + 1) * spacing_rad)
                    n_candidates = min(max(int(gap_estimate / spacing_rad), 1), n_distance_cap, max_extend_nodes, ring_room())
                    candidate_theta = theta[-1] + dtheta * np.arange(1, n_candidates + 1)
                    n_new = self._count_open_prefix(candidate_theta, line.phi, neighbours)
                    if n_new > 0:
                        hc_seed, hm_seed, elev_seed, reason_seed = _end_seed(True)
                        new_theta = candidate_theta[:n_new]
                        theta = np.append(theta, new_theta)
                        elevation = np.append(elevation, np.full(n_new, elev_seed))
                        for name, fill in _fill_new_nodes(n_new, hc_seed, hm_seed, reason_seed).items():
                            persistent_fields[name] = np.append(persistent_fields[name], fill)
                        if world.debug_diagnostics:
                            phase_budget.record(
                                world, self, "line_end_arc_grow",
                                np.array([]), np.array([]), np.array([], dtype=persistent_fields["crust_type_code"].dtype),
                                persistent_fields["crustal_thickness_m"][-n_new:],
                                persistent_fields["mantle_lithosphere_thickness_m"][-n_new:],
                                persistent_fields["crust_type_code"][-n_new:],
                            )
            else:
                _stretch_end(-1, 1.0, dist[-1], direction[-1])

        if not contested[0] and dist[0] > extend_threshold_rad and ring_room() > 0:
            if arc_end_low:
                if not suppress_growth:
                    gap_estimate = min(dist[0], (n_distance_cap + 1) * spacing_rad)
                    n_candidates = min(max(int(gap_estimate / spacing_rad), 1), n_distance_cap, max_extend_nodes, ring_room())
                    candidate_theta = theta[0] - dtheta * np.arange(1, n_candidates + 1)
                    n_new = self._count_open_prefix(candidate_theta, line.phi, neighbours)
                    if n_new > 0:
                        hc_seed, hm_seed, elev_seed, reason_seed = _end_seed(True)
                        new_theta = candidate_theta[:n_new][::-1]
                        theta = np.insert(theta, 0, new_theta)
                        elevation = np.insert(elevation, 0, np.full(n_new, elev_seed))
                        for name, fill in _fill_new_nodes(n_new, hc_seed, hm_seed, reason_seed).items():
                            persistent_fields[name] = np.insert(persistent_fields[name], 0, fill)
                        if world.debug_diagnostics:
                            phase_budget.record(
                                world, self, "line_end_arc_grow",
                                np.array([]), np.array([]), np.array([], dtype=persistent_fields["crust_type_code"].dtype),
                                persistent_fields["crustal_thickness_m"][:n_new],
                                persistent_fields["mantle_lithosphere_thickness_m"][:n_new],
                                persistent_fields["crust_type_code"][:n_new],
                            )
            else:
                _stretch_end(0, -1.0, dist[0], direction[0])

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

        drop_extremes: dict[float, str] = {}
        for extreme, phi_key in (("lo", rows_left[0]), ("hi", rows_left[-1])):
            n_nodes, n_contested = contested_by_phi.get(phi_key, (0.0, 0.0))
            fraction = n_contested / n_nodes if n_nodes else 0.0
            if fraction < LEADING_ROW_CONTESTED_FRACTION:
                tracker.pop(extreme, None)
                continue
            prev_phi, prev_years = tracker.get(extreme, (None, 0.0))
            accumulated = (prev_years if prev_phi == phi_key else 0.0) + years
            if accumulated >= LEADING_ROW_RETREAT_SUSTAINED_YEARS:
                drop_extremes[phi_key] = extreme
                tracker.pop(extreme, None)
            else:
                tracker[extreme] = (phi_key, accumulated)

        if not drop_extremes:
            return new_lines
        return self._accrete_dropped_row_volume(new_lines, drop_extremes)

    def _accrete_dropped_row_volume(
        self, lines: list[ElevationLine], drop_extremes: dict[float, str]
    ) -> list[ElevationLine]:
        """Conserve the crustal volume of the whole leading rows `_retreat_contested_leading_
        rows` just dropped (`drop_extremes`: this plate's own phi keys -> "lo"/"hi", which
        extreme each was) by thrusting each dropped row's summed Hc onto the plate's *new*
        leading row at that same extreme -- spread evenly across every node of that whole
        row, since the row being consumed spans the plate's full theta width, unlike
        `_redistribute_accreted_column`'s single line-end case, where there's a narrower
        "end" to concentrate the volume on. Node area is constant per node, so summed Hc *is*
        the conserved volume; capped at `SUTURE_ACCRETION_MAX_HC_M` same as ordinary suture
        accretion (the overflow delaminates -- see that constant's own comment)."""
        rho_c = self.crust_density()
        kept = [ln for ln in lines if round(float(ln.phi), 6) not in drop_extremes]
        if not kept:
            return kept
        kept_phis = sorted({round(float(ln.phi), 6) for ln in kept})

        for phi_key, extreme in drop_extremes.items():
            dropped_hc_sum = sum(
                float(ln.crustal_thickness_m.sum()) for ln in lines if round(float(ln.phi), 6) == phi_key
            )
            if dropped_hc_sum <= 0.0:
                continue
            target_phi = kept_phis[0] if extreme == "lo" else kept_phis[-1]
            target_idx = [i for i, ln in enumerate(kept) if round(float(ln.phi), 6) == target_phi]
            target_n = sum(len(kept[i]) for i in target_idx)
            if target_n == 0:
                continue
            add_hc_per_node = dropped_hc_sum / target_n
            for i in target_idx:
                ln = kept[i]
                hc = ln.crustal_thickness_m
                hm = ln.mantle_lithosphere_thickness_m
                before = lithosphere.isostatic_elevation(hc, hm, rho_c)
                new_hc = np.minimum(hc + add_hc_per_node, SUTURE_ACCRETION_MAX_HC_M)
                # Also capped (issue #161) -- see _redistribute_accreted_column's own note on
                # why a large new_hc/hc ratio here (a thin row absorbing a whole dropped row's
                # volume) can otherwise carry Hm past its own ceiling.
                new_hm = np.minimum(hm * (new_hc / hc), lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)
                after = lithosphere.isostatic_elevation(new_hc, new_hm, rho_c)
                new_elevation = rheology.clip_elevation_bounds(ln.elevation + (after - before))
                kept[i] = ln.replace(crustal_thickness_m=new_hc, mantle_lithosphere_thickness_m=new_hm, elevation=new_elevation)
        return kept

    def _seed_and_erupt_new_nodes(
        self, world: "World", line_index: int, world_pts: np.ndarray, thin_ratio: float,  # noqa: F821
        hc0: float, hm0: float, amp: float, texture: "terrain_noise.FractalTexture"
    ) -> dict[str, np.ndarray]:
        """Brand-new nodes carry no prior column to conserve, so -- exactly like ordinary
        divergent thinning and `_grow_or_shrink_line_for_deform`'s `_stretch_end` -- they are
        seeded thin (`thin_ratio` share of the oceanic reference column, `growth_seed_thickness`)
        and run straight through `_erupt_melted_nodes`, the same decompression-melting/magma-
        upwelling path every other new-crust event in `deform()` uses. This makes every node
        this produces a real eruption (typed, `is_volcano`-stamped, timed) rather than a
        distinct silent "spawn" concept -- shared by `_claim_adjacent_territory` and
        `_fill_corner_notch_frontier`, the two callers that ever originate genuinely new areal crust."""
        n = len(world_pts)
        hc = np.full(n, hc0 * thin_ratio) + amp * texture.sample(world_pts)
        hm = np.full(n, hm0 * thin_ratio)
        elevation = lithosphere.isostatic_elevation(hc, hm, self.crust_density())
        crust_type_code = np.zeros(n, dtype=np.int8)
        is_volcano = np.zeros(n, dtype=bool)
        volcano_remaining = np.zeros(n)
        melting = hc < rheology.RIFT_CRITICAL_THICKNESS_M
        _erupt_melted_nodes(
            world, self.plate_id, line_index,
            hc, hm, crust_type_code, is_volcano, volcano_remaining,
            melting, elevation,
        )
        elevation = lithosphere.isostatic_elevation(hc, hm, lithosphere.node_crust_density(crust_type_code, self.crust_type))
        return {
            "elevation": elevation,
            "crustal_thickness_m": hc,
            "mantle_lithosphere_thickness_m": hm,
            "crust_type_code": crust_type_code,
            "is_volcano": is_volcano,
            "volcano_active_years_remaining": volcano_remaining,
            "elev_change_reason": np.full(n, ELEV_CHANGE_NEW_CRUST, dtype=float),
            "node_created_years": np.full(n, world.elapsed_years, dtype=float),
        }

    def _claim_adjacent_territory(
        self, world: "World", neighbours: list, spacing_rad: float, neighbour_tree: cKDTree | None = None  # noqa: F821
    ) -> None:
        """Same shape as `PlateWithLines._claim_adjacent_territory` -- one or more brand-new
        phi rows just past this plate's own phi extremes, where open -- seeded with fresh Hc/Hm
        (oceanic reference, see `growth_seed_thickness`) plus `terrain_noise.FractalTexture`
        on Hc (an extension of an already-shaped plate, so texture rather than a fresh
        orogen), rather than a flat elevation baseline. Keyed off `(world.seed, plate_id,
        _TERRAIN_SEED_TAG)` so the texture stays attached to this plate as it grows.

        Gated on genuine phi-direction separation (`rheology.stretch_components`' `phi_gap`,
        evaluated at the candidate row's own midpoint via this plate's own
        `_separation_components`) rather than firing unconditionally whenever the ground ahead
        happens to be open -- a purely theta-aligned gap is `_grow_or_shrink_line_for_deform`'s
        own end-stretch to close, not this row-claim's. When it does fire, the new row is not
        seeded at the full reference column "for free": that volume is instead drawn down
        across the new row *and* `K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION` of this plate's own
        existing outermost rows on that side, in proportion to `phi_gap`'s own share of one
        full row spacing -- mass conservation, the phi-direction counterpart of the theta case's
        own stretch-thinning.

        Loops up to `MAX_CLAIM_ROWS_PER_STEP` rows per direction (bounded overall by
        `MAX_CLAIM_NODES_PER_STEP`), re-deriving the reference row each iteration, instead of
        claiming exactly one row and waiting for next step's call to claim the next: one row per
        plate per step structurally cannot keep pace with a fast-widening gap (see
        MAX_CLAIM_ROWS_PER_STEP's own comment). "Open" is a real coverage/proximity check
        (`COVERAGE_RADIUS_MULT * spacing_rad` against every neighbour's own node cloud, the
        same tolerance gaps.py's whole-sphere sweep uses) rather than `_contested_by_any`'s
        weaker "not inside a neighbour's polygon", and each candidate row is split into however
        many contiguous open runs it actually has (`split_into_contiguous_runs`) instead of
        always claiming the full old row's theta span with holes in it -- what lets a
        triple-junction wedge narrow or widen row by row instead of insisting on full-width
        rectangular strips.

        `neighbour_tree`, when passed (deform()'s own call does, sharing it with whichever
        corner-notch closer runs right after this in the same deform() call), is a cKDTree
        already built over every neighbour's concatenated node cloud -- built fresh here only
        for a caller without one (a direct unit-test call, or `neighbours=[]`)."""
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
        if neighbour_tree is None:
            neighbour_points = [p.all_points_and_elevation()[0] for p in neighbours if p.node_count() > 0]
            neighbour_tree = cKDTree(np.concatenate(neighbour_points, axis=0)) if neighbour_points else None
        coverage_radius_rad = COVERAGE_RADIUS_MULT * spacing_rad
        line_index_by_id = {id(row): i for i, row in enumerate(self.lines)}
        new_lines: list[ElevationLine] = []
        thinned: dict[int, ElevationLine] = {}

        def is_open(world_pts: np.ndarray) -> np.ndarray:
            if neighbour_tree is None:
                return np.ones(len(world_pts), dtype=bool)
            dist, _ = neighbour_tree.query(world_pts)
            return dist > coverage_radius_rad

        for direction in (-1, 1):
            reference = ordered[0] if direction < 0 else ordered[-1]
            # Rolling window of the K rows nearest the *current* frontier, for mass-
            # conservation draw-down -- starts as the plate's own original outermost K rows,
            # then shifts outward to include each newly claimed row in turn (dropping whichever
            # row is now farthest from the frontier). Recomputing this from the ORIGINAL
            # `ordered` on every loop iteration instead (a real bug caught by
            # test_continent_continent_suture_consumes_its_overlap_as_mass_conserving_accretion)
            # would draw the *same* original edge down once per row claimed this call --
            # correct for one row, but a multi-row claim would over-thin that edge N-fold for
            # no physical reason, since only the newest row's own gap should charge against
            # rows actually near it.
            recent_rows = ordered[:K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION] if direction < 0 else ordered[-K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION:]
            rows_claimed = 0
            nodes_claimed = 0
            while rows_claimed < MAX_CLAIM_ROWS_PER_STEP and nodes_claimed < MAX_CLAIM_NODES_PER_STEP:
                new_phi = reference.phi + direction * spacing_rad
                if abs(new_phi) > max_phi_limit:
                    break
                dtheta = spacing_rad / max(np.cos(new_phi), 1e-3)
                span = reference.theta[-1] - reference.theta[0]
                n_cols = max(int(round(span / dtheta)) + 1, 1)
                theta_candidates = reference.theta[0] + dtheta * np.arange(n_cols)
                world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.full(n_cols, new_phi), theta_candidates))

                open_mask = is_open(world_pts)
                if not np.any(open_mask):
                    break

                # Real coverage, not `_contested_by_any`'s polygon test, so a candidate row
                # against a diagonal neighbour boundary comes back with only its genuinely
                # uncovered stretch open -- probe_line lets split_into_contiguous_runs find
                # every contiguous open run instead of one all-or-nothing span.
                probe_line = ElevationLine(phi=new_phi, theta=theta_candidates, elevation=np.zeros(n_cols))
                probe_line = probe_line.masked(open_mask)
                runs = split_into_contiguous_runs(probe_line, dtheta)
                if not runs:
                    break

                # How far, and in what direction, does a genuine gap actually extend here? --
                # the same nearest-neighbour-node query torque.gather_boundary_force_inputs
                # uses for a live node, evaluated fresh at the largest run's own midpoint (one
                # representative point stands in for the whole candidate row, matching how
                # coarse a once-per-step per-plate-edge check this already was).
                largest_run = max(runs, key=len)
                rep_theta = float(largest_run.theta[len(largest_run) // 2])
                rep_world = geometry.to_world(self.frame, geometry.local_xyz(np.array([new_phi]), np.array([rep_theta])))[0]
                if neighbour_tree is not None:
                    best_dist, idx = neighbour_tree.query(rep_world)
                    direction_world = geometry.normalize(neighbour_tree.data[idx] - rep_world)
                else:
                    best_dist, direction_world = np.inf, np.zeros(3)
                sep_theta, sep_phi = self._separation_components(world, new_phi, rep_theta, direction_world)
                gap_estimate = min(float(best_dist), spacing_rad) if np.isfinite(best_dist) else spacing_rad
                _, phi_gap = rheology.stretch_components(sep_theta, sep_phi, gap_estimate)
                if phi_gap <= 0.0:
                    break

                share_count = len(recent_rows) + 1  # + the new row itself
                stretch_fraction = float(np.clip(phi_gap / spacing_rad, 0.0, 1.0))
                thin_ratio = 1.0 - stretch_fraction * (share_count - 1) / share_count

                for run in runs:
                    world_run = geometry.to_world(self.frame, geometry.local_xyz(np.full(len(run), new_phi), run.theta))
                    seeded = self._seed_and_erupt_new_nodes(
                        world, len(self.lines) + len(new_lines), world_run, thin_ratio, hc0, hm0, amp, texture
                    )
                    new_row = ElevationLine(phi=new_phi, theta=run.theta, **seeded)
                    new_lines.append(new_row)
                    # A row claimed this call can itself become a `recent_rows` mass-
                    # conservation source in a later iteration (see below) -- give it a stable
                    # index too, not just the rows that existed before this call.
                    line_index_by_id[id(new_row)] = len(self.lines) + len(new_lines) - 1
                    nodes_claimed += len(run)
                rows_claimed += 1

                # Pull the same fractional share out of the existing rows nearest this claim --
                # thinned in place (a row already thinned by an earlier iteration, or by the
                # *other* phi extreme's own claim this same call, keeps compounding, which is
                # fine: repeated genuine gaps opening on the same side really should thin it
                # further each time).
                for original_row in recent_rows:
                    row = thinned.get(id(original_row), original_row)
                    new_hc_row = row.crustal_thickness_m * thin_ratio
                    new_hm_row = row.mantle_lithosphere_thickness_m * thin_ratio
                    melting_row = (row.crustal_thickness_m >= rheology.RIFT_CRITICAL_THICKNESS_M) & (
                        new_hc_row < rheology.RIFT_CRITICAL_THICKNESS_M
                    )
                    crust_type_row = row.crust_type_code.copy()
                    is_volcano_row = row.is_volcano.copy()
                    volcano_remaining_row = row.volcano_active_years_remaining.copy()
                    _erupt_melted_nodes(
                        world, self.plate_id, line_index_by_id.get(id(original_row), -1),
                        new_hc_row, new_hm_row, crust_type_row, is_volcano_row, volcano_remaining_row,
                        melting_row, row.elevation,
                    )
                    new_elevation_row = lithosphere.isostatic_elevation(
                        new_hc_row, new_hm_row, lithosphere.node_crust_density(crust_type_row, self.crust_type)
                    )
                    thinned[id(original_row)] = row.replace(
                        crustal_thickness_m=new_hc_row,
                        mantle_lithosphere_thickness_m=new_hm_row,
                        crust_type_code=crust_type_row,
                        is_volcano=is_volcano_row,
                        volcano_active_years_remaining=volcano_remaining_row,
                        elevation=new_elevation_row,
                    )

                # Continue outward from the widest surviving run -- the next row's own span is
                # derived from its reference row's theta extent, and a fragmented claim narrows
                # naturally toward whichever piece the gap actually left open.
                reference = max(new_lines[-len(runs):], key=len)
                # Shift the mass-conservation window outward: the row just claimed becomes the
                # nearest neighbour to the *next* frontier, displacing whichever row is now
                # farthest from it.
                recent_rows = [reference] + recent_rows[: K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION - 1]

        if new_lines or thinned:
            # `recent_rows` can include a `new_lines` entry from an earlier iteration of the
            # same direction's loop (the row just claimed becomes a mass-conservation source
            # for the next one) -- so `thinned` may hold an updated version of a *new* line,
            # not just an original one. Apply it to both, or a later thinning pass on a
            # just-claimed row would silently apply to nothing.
            new_lines = [thinned.get(id(line), line) for line in new_lines]
            self.set_lines([thinned.get(id(line), line) for line in self.lines] + new_lines)

    def _fill_corner_notch_frontier(
        self, world: "World", neighbours: list, spacing_rad: float, years: float, neighbour_tree: cKDTree | None = None  # noqa: F821
    ) -> None:
        """Generalized fallback for the sub-row diagonal residual `_stretch_end` (theta-axis-only)
        and `_claim_adjacent_territory` (phi-axis-only, whole rows) structurally cannot reach:
        the shape independently-oriented plate grids leave open at a triple junction, or along a
        curved two-plate rift boundary that isn't a clean theta/phi-aligned edge. Called from the
        same place in `deform()`, right after `_claim_adjacent_territory` -- see
        `gap_fill_frontier.py`'s own module docstring for the growth algorithm. Once this
        method's own window/neighbour-reach gate finds candidate notch points,
        `gap_fill_frontier.fill_gap_by_growing_plates` grows *this plate's own existing lines*
        into them node by node instead of emitting disjoint new ones.

        `neighbour_tree`, when passed (deform()'s own call does, sharing it with
        `_claim_adjacent_territory`'s own call right before this in the same deform() call --
        neither call touches any neighbour's own nodes, so the same tree serves both), is a
        cKDTree already built over every neighbour's concatenated node cloud -- built fresh
        here only for a caller without one.

        Imports `gap_fill_frontier` locally (not at module scope) since that module itself
        imports back from here (`_TERRAIN_SEED_TAG`, `_erupt_melted_nodes`,
        `growth_seed_thickness`) -- same deferred-import shape `deform()`'s own `from . import
        faults` already uses to avoid a circular top-level import."""
        from . import gap_fill_frontier

        if neighbour_tree is None:
            neighbour_points = [p.all_points_and_elevation()[0] for p in neighbours if p.node_count() > 0]
            neighbour_tree = cKDTree(np.concatenate(neighbour_points, axis=0)) if neighbour_points else None
        if neighbour_tree is None:
            world.log_corner_notch({"plate_id": self.plate_id, "outcome": "no_neighbours", "nodes_added": 0, "algorithm": "frontier"})
            return
        own_points, _ = self.all_points_and_elevation()
        lines_with_nodes = [line for line in self.lines if len(line) > 0]
        if len(own_points) == 0 or not lines_with_nodes:
            world.log_corner_notch({"plate_id": self.plate_id, "outcome": "no_own_lines", "nodes_added": 0, "algorithm": "frontier"})
            return

        coverage_radius_rad = COVERAGE_RADIUS_MULT * spacing_rad
        window_rad = max(CORNER_NOTCH_MIN_WINDOW_ROWS * spacing_rad, mantle.MAX_PLATE_RATE * years)
        neighbour_reach_rad = window_rad + CORNER_NOTCH_NEIGHBOUR_REACH_MARGIN_MULT * spacing_rad

        max_abs_phi = np.pi / 2 - spacing_rad / 2  # matches iter_local_lattice's own bound
        max_phi_limit = np.pi / 2 - POLE_CAP_MARGIN_MULT * spacing_rad
        line_phis = np.array([line.phi for line in lines_with_nodes])
        phi_lo = max(line_phis.min() - window_rad, -max_abs_phi)
        phi_hi = min(line_phis.max() + window_rad, max_abs_phi)
        row_lo = int(np.floor((phi_lo + max_abs_phi) / spacing_rad))
        row_hi = int(np.ceil((phi_hi + max_abs_phi) / spacing_rad))
        phi_values = -max_abs_phi + spacing_rad * np.arange(row_lo, row_hi + 1)
        phi_values = phi_values[np.abs(phi_values) <= max_phi_limit]

        gap_chunks: list[np.ndarray] = []
        if len(phi_values) > 0:
            # Vectorized replacement for a per-row `min(lines_with_nodes, key=...)` scan:
            # a row split into multiple contiguous runs (e.g. one wrapping the antimeridian,
            # see split_into_contiguous_runs) can leave two-plus lines sharing the same phi,
            # and `min` breaks that tie by first list occurrence -- `np.unique`'s
            # `return_index` gives exactly that (each unique phi's *first* original index),
            # so searching the deduped, sorted phis and mapping back through it reproduces
            # `min`'s tie-break exactly (verified against a brute-force `min` scan across
            # thousands of randomized duplicate-phi cases) without ever re-scanning
            # `lines_with_nodes` per row. Batching every row's world points into one cKDTree
            # query below (instead of one query per row) is this method's other real cost --
            # together these were the two costs profiled here; detection logic/results are
            # otherwise unchanged.
            sorted_phis, first_idx = np.unique(line_phis, return_index=True)
            theta_first = np.array([lines_with_nodes[i].theta[0] for i in first_idx], dtype=float)
            theta_last = np.array([lines_with_nodes[i].theta[-1] for i in first_idx], dtype=float)
            if len(sorted_phis) == 1:
                nearest_idx = np.zeros(len(phi_values), dtype=int)
            else:
                insert_idx = np.clip(np.searchsorted(sorted_phis, phi_values), 1, len(sorted_phis) - 1)
                left_idx = insert_idx - 1
                choose_left = np.abs(sorted_phis[left_idx] - phi_values) <= np.abs(sorted_phis[insert_idx] - phi_values)
                nearest_idx = np.where(choose_left, left_idx, insert_idx)

            cos_phi = np.maximum(np.cos(phi_values), 1e-3)
            dtheta = spacing_rad / cos_phi
            margin = window_rad / cos_phi
            theta_lo = theta_first[nearest_idx] - margin
            theta_hi = theta_last[nearest_idx] + margin
            n_theta = np.maximum(np.round((theta_hi - theta_lo) / dtheta).astype(int) + 1, 1)
            n_theta = np.minimum(n_theta, np.maximum(np.round(2.0 * np.pi / dtheta).astype(int), 1))

            row_thetas = [lo + dth * np.arange(n) for lo, dth, n in zip(theta_lo, dtheta, n_theta)]
            row_phis = [np.full(n, phi) for phi, n in zip(phi_values, n_theta)]
            world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.concatenate(row_phis), np.concatenate(row_thetas)))

            neighbour_dist, _ = neighbour_tree.query(world_pts, workers=query_workers(len(world_pts)))
            near_neighbour = neighbour_dist <= neighbour_reach_rad
            covered_by_neighbour = neighbour_dist <= coverage_radius_rad
            base_mask = near_neighbour & ~covered_by_neighbour

            offset = 0
            for n in n_theta:
                row_mask = base_mask[offset : offset + n]
                if np.any(row_mask):
                    gap_chunks.append(world_pts[offset : offset + n][row_mask])
                offset += n

        if not gap_chunks:
            world.log_corner_notch({
                "plate_id": self.plate_id, "outcome": "no_candidate_rows", "nodes_added": 0, "algorithm": "frontier",
                "window_rad": float(window_rad), "phi_lo": float(phi_lo), "phi_hi": float(phi_hi),
            })
            return

        gap_points = np.concatenate(gap_chunks, axis=0)
        # Cap the frontier walk to this method's own `window_rad` reach --
        # `ceil(window_rad / connect_radius)` hops -- rather than letting
        # `fill_gap_by_growing_plates` auto-size from the candidate cloud's own (potentially
        # much wider) bounding sphere. Without this the notch-filler outgrows a continental
        # plate's leading edge faster than `_retreat_contested_leading_rows` can retreat it
        # against a parallel suture (regression caught by
        # test_lithosphere_contested_leading_row_is_dropped_after_sustained_override).
        connect_radius_rad = DEFRAG_CONNECT_RADIUS_MULT * spacing_rad
        max_hops = max(1, int(np.ceil(window_rad / connect_radius_rad)))
        added = gap_fill_frontier.fill_gap_by_growing_plates(world, gap_points, [self], spacing_rad, max_hops=max_hops)
        nodes_added = added.get(self.plate_id, 0)
        world.log_corner_notch({
            "plate_id": self.plate_id, "outcome": "claimed" if nodes_added else "no_claim", "nodes_added": nodes_added,
            "algorithm": "frontier", "window_rad": float(window_rad), "phi_lo": float(phi_lo), "phi_hi": float(phi_hi),
        })

    # -- Merge/split: carry Hc/Hm through, not just elevation -------------------------------

    def merge_with(self, other: "LithospherePlate", spacing_rad: float, coverage_radius_rad: float, other_points_xyz: np.ndarray) -> None:
        own_points, _ = self.all_points_and_elevation()
        other_points, _ = other.all_points_and_elevation()
        inertia_self = lithosphere.moment_of_inertia_tensor(
            own_points, self.collect("crustal_thickness_m"), self.collect("mantle_lithosphere_thickness_m"), self.node_crust_density(), spacing_rad
        )
        inertia_other = lithosphere.moment_of_inertia_tensor(
            other_points, other.collect("crustal_thickness_m"), other.collect("mantle_lithosphere_thickness_m"), other.node_crust_density(), spacing_rad
        )
        self._merge_nodes_with(other, spacing_rad, coverage_radius_rad, other_points_xyz)
        self.set_omega(torque.merge_omega(self, inertia_self, other, inertia_other))
        self.reset_age()

    def _merge_nodes_with(self, other: "LithospherePlate", spacing_rad: float, coverage_radius_rad: float, other_points_xyz: np.ndarray) -> None:
        keep_pts, _ = self.all_points_and_elevation()
        absorb_pts, _ = other.all_points_and_elevation()
        exclude_tree = cKDTree(other_points_xyz) if len(other_points_xyz) else None
        self.set_lines(
            _merge_lines_from_resample(
                self.frame,
                keep_pts, self.collect("crustal_thickness_m"), self.collect("mantle_lithosphere_thickness_m"),
                absorb_pts, other.collect("crustal_thickness_m"), other.collect("mantle_lithosphere_thickness_m"),
                coverage_radius_rad, spacing_rad, exclude_tree,
            )
        )
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

        # A daughter's own crust_type is the majority of what its nodes actually are, not a
        # blind copy of the parent's -- see elevation_lines.majority_crust_type. A no-op
        # (returns self.crust_type unchanged) unless this plate has ever had a magma-typing
        # event (rift decompression melting) whose composition ended up lopsided across the cut.
        crust_type_a = majority_crust_type(lines_a, self.crust_type)
        crust_type_b = majority_crust_type(lines_b, self.crust_type)
        plate_a = LithospherePlate(plate_id=self.plate_id, frame=self.frame.copy(), crust_type=crust_type_a, lines=lines_a)
        plate_b = LithospherePlate(plate_id=new_id, frame=self.frame.copy(), crust_type=crust_type_b, lines=lines_b)
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
        # pole-winding notes in GitHub issue #119) would put every node "next to the rift" and thin
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

    def relattice(self, spacing_rad: float) -> None:
        """Refit this plate's lattice to its own current outline and redistribute its
        existing total crustal volume onto the fresh node set (GitHub issue #119, "Continental
        ratchet: solution design," mechanism 4, "periodic conservative continental
        re-lattice") -- the 2-D generalisation of `elevation_lines.regularize_line`, which only
        ever re-evens spacing *within* one already-existing row, preserving its two endpoints
        exactly. That per-row guarantee is exactly what it can't fix: `_grow_or_shrink_line_
        for_deform` grows each row's own ends independently, one node at a time, at whatever
        rate that row's own local contact happens to demand, so nothing stops row-to-row phase
        drift (a diagonal, staircase boundary from repeated row-end extension -- issue #119's
        "streaking" symptom, a thin triangular tongue grown one row-end at a time) from
        compounding indefinitely even while every individual row stays evenly spaced.
        Continental crust only: an oceanic plate's footprint is already self-bounding via
        subduction, so there is no such drift here worth periodically re-fitting.

        Ownership of each fresh lattice site is `contains_batch` -- this plate's own polygon
        test, exact against its outline as of the last `deform()` -- not a coverage-radius
        dilation of the existing node cloud (`Plate.grow_into`, used for merges). A
        radius-based resample was already rejected for routine per-step use because its
        coverage radius around even a handful of points reconstructs far more lattice area
        than they actually cover (see docs/simulation-model.md's "Claiming adjacent
        territory"); testing the outline directly instead reproduces exactly this plate's
        existing footprint -- no smaller, no larger -- just resampled onto the canonical,
        evenly-phased lattice `iter_local_lattice` builds from scratch.

        Every `OPTIONAL_FIELDS` value (Hc/Hm included) is carried onto each new site from its
        nearest surviving node: a 2-D scatter resample has no single ordered axis to
        `np.interp` along the way `regularize_line` does, so nearest-neighbour is the natural
        generalisation (the same choice `regularize_line` itself already makes for its
        categorical fields). Because node area is constant, a nearest-neighbour carry alone
        doesn't exactly conserve total crustal volume -- a different-shaped node set
        overweights whichever few old nodes end up nearest the most new sites near the
        boundary -- so `crustal_thickness_m` (and `mantle_lithosphere_thickness_m`, scaled by
        the same ratio: a thicker resampled column carries a proportionally thicker attached
        mantle lid, the same convention `_redistribute_accreted_column` uses) is rescaled by
        one uniform factor afterward so the plate's `sum(Hc)` -- its total crustal volume,
        since per-node area is constant -- comes out exactly where it started."""
        if self.crust_type != "continental":
            return
        own_points, _ = self.all_points_and_elevation()
        if len(own_points) == 0:
            return
        field_names = ("elevation",) + ElevationLine.OPTIONAL_FIELDS
        own_fields = {name: self.collect(name) for name in field_names}
        total_hc_before = float(np.sum(own_fields["crustal_thickness_m"]))
        if total_hc_before <= 0.0:
            return

        from .elevation_lines import iter_local_lattice

        tree = cKDTree(own_points)
        raw_rows: list[tuple[float, np.ndarray, dict[str, np.ndarray]]] = []
        for phi, theta_candidates, world_pts in iter_local_lattice(self.frame, spacing_rad=spacing_rad):
            owned = self.contains_batch(world_pts)
            if not np.any(owned):
                continue
            _, idx = tree.query(world_pts[owned])
            row_fields = {name: values[idx] for name, values in own_fields.items()}
            raw_rows.append((phi, theta_candidates[owned], row_fields))

        if not raw_rows:
            return

        total_hc_after = float(sum(np.sum(fields["crustal_thickness_m"]) for _, _, fields in raw_rows))
        hc_scale = total_hc_before / total_hc_after if total_hc_after > 0.0 else 1.0

        # Conserve volume via the uniform `hc_scale` above, but never push a node past
        # `lithosphere.MAX_CRUSTAL_THICKNESS_M` doing it (issue #161) -- a node already near
        # the ceiling going into this resample would otherwise cross it under an ordinary >1
        # rescale. Clip first, then spread whatever got clipped off across every *other*
        # node's own remaining headroom below the ceiling, weighted by how much headroom each
        # has (one pass is enough in practice: relattice only ever corrects a small drift in
        # total volume, not redistributes a large fraction of it) -- the same conserve-then-
        # spread idiom `apply_convergent_deformation`'s own overflow uses, just plate-wide
        # instead of a local ring. A true residual (every node already pinned at the ceiling)
        # finally delaminates, same as everywhere else this ceiling applies.
        pre_scale_hc = np.concatenate([fields["crustal_thickness_m"] for _, _, fields in raw_rows])
        scaled_hc = pre_scale_hc * hc_scale
        capped_hc = np.minimum(scaled_hc, lithosphere.MAX_CRUSTAL_THICKNESS_M)
        overflow = float(np.sum(scaled_hc - capped_hc))
        if overflow > 0.0:
            headroom = lithosphere.MAX_CRUSTAL_THICKNESS_M - capped_hc
            total_headroom = float(np.sum(headroom))
            if total_headroom > 0.0:
                capped_hc = np.minimum(capped_hc + overflow * (headroom / total_headroom), lithosphere.MAX_CRUSTAL_THICKNESS_M)
        final_hc = np.maximum(capped_hc, lithosphere.MIN_CRUSTAL_THICKNESS_M)
        hm_ratio = final_hc / pre_scale_hc

        new_lines = []
        offset = 0
        for phi, theta_owned, fields in raw_rows:
            n = len(theta_owned)
            overrides = dict(fields)
            elevation = overrides.pop("elevation")
            overrides["crustal_thickness_m"] = final_hc[offset : offset + n]
            overrides["mantle_lithosphere_thickness_m"] = np.clip(
                overrides["mantle_lithosphere_thickness_m"] * hm_ratio[offset : offset + n],
                lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M,
                lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M,
            )
            new_lines.append(ElevationLine(phi=phi, theta=theta_owned, elevation=elevation, **overrides))
            offset += n

        self.set_lines(new_lines)
        lithosphere.sync_plate_elevation(self)


def _merge_lines_from_resample(
    frame: np.ndarray,
    keep_points: np.ndarray,
    keep_hc: np.ndarray,
    keep_hm: np.ndarray,
    absorb_points: np.ndarray,
    absorb_hc: np.ndarray,
    absorb_hm: np.ndarray,
    coverage_radius_rad: float,
    spacing_rad: float,
    exclude_tree: cKDTree | None = None,
) -> list[ElevationLine]:
    """Resample the fusing pair's own node clouds onto a fresh local lattice (see
    `_merge_nodes_with`'s own docstring) -- `keep`/`absorb` queried *separately*, not
    (as an earlier version of this did) via one nearest-neighbor lookup against their naive
    concatenation.

    That distinction matters exactly where a merge is most consequential: the deeply
    overlapping suture band a sustained collision leaves behind before the two plates
    actually fuse (see merge_split.py's own `FORCED_MERGE_OVERLAP_FRACTION` / ordinary
    closing-rate merge threshold, both tuned to only fire once real territorial overlap has
    built up). There, `keep` and `absorb` each still carry their own full-thickness column
    at roughly the same location -- two independent lattices interleaved at close to the
    *same* areal density the new, merged lattice also targets. A single nearest-of-the-union
    lookup at that density picks up only whichever one of the two happened to be closer to
    each new site and silently drops the other's entire column -- confirmed directly
    (2026-09-12, investigating a "land keeps declining across a continental collision"
    report): every one of 7 sampled merge events across 5 seeds of the `two_colliding_pairs`
    Debugging World lost 9-17% of the merging pair's total continental crustal volume at the
    moment of fusion, with no discrete accounting for where it went -- confirming "some
    volume gets consumed entirely by an early fusion" as a real bug, not just a hypothesis.

    Fix: query `keep`'s and `absorb`'s own clouds independently at every new lattice site.
    Where only one has a node within `coverage_radius_rad`, behavior is unchanged from
    before (an ordinary single-column resample). Where *both* do -- genuine suture overlap
    -- their Hc is summed (real collisions really do stack two independent columns into one
    thicker one; this is the literal orogenic-thrust-wedge physics
    `_redistribute_accreted_column` already applies to an ordinary contested-edge retreat,
    just applied here to the deep-overlap remainder that survives all the way to the final
    merge instead of retreating earlier), capped at the same `SUTURE_ACCRETION_MAX_HC_M` an
    ordinary suture uses (the overflow delaminates -- see that constant's own comment), with
    `keep`'s own Hm scaled by the same Hc growth ratio `_redistribute_accreted_column` uses
    (a doubled crustal column does not imply a doubled mantle-lithosphere lid) rather than
    also summed.

    `exclude_tree` (every *other*, uninvolved plate's own nodes) keeps its original meaning
    and does not change this plate's resulting footprint, only how Hc/Hm is sampled within
    it: a candidate site closer to some third plate than to either fusing plate is still
    never claimed here, the same exclusivity `plates.generate_plates`' Voronoi tiling
    guarantees at generation.

    `elevation` on the returned lines is a placeholder (zeros); the caller must run
    `lithosphere.sync_plate_elevation` right after to derive the real isostatic value."""
    from .elevation_lines import iter_local_lattice

    keep_tree = cKDTree(keep_points) if len(keep_points) else None
    absorb_tree = cKDTree(absorb_points) if len(absorb_points) else None

    def query(tree: cKDTree | None, world_pts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if tree is None:
            return np.full(len(world_pts), np.inf), np.zeros(len(world_pts), dtype=int)
        return tree.query(world_pts)

    lines: list[ElevationLine] = []
    for phi, theta_candidates, world_pts in iter_local_lattice(frame, spacing_rad=spacing_rad):
        keep_dist, keep_idx = query(keep_tree, world_pts)
        absorb_dist, absorb_idx = query(absorb_tree, world_pts)
        own_dist = np.minimum(keep_dist, absorb_dist)

        if exclude_tree is not None:
            exclude_dist, _ = exclude_tree.query(world_pts)
        else:
            exclude_dist = np.full(len(world_pts), np.inf)
        owned = (own_dist < coverage_radius_rad) & (own_dist < exclude_dist)
        if not np.any(owned):
            continue

        keep_close = keep_dist < coverage_radius_rad
        absorb_close = absorb_dist < coverage_radius_rad
        both = keep_close & absorb_close

        hc = np.zeros(len(world_pts))
        hm = np.zeros(len(world_pts))
        only_keep = keep_close & ~both
        hc[only_keep] = keep_hc[keep_idx[only_keep]]
        hm[only_keep] = keep_hm[keep_idx[only_keep]]
        only_absorb = absorb_close & ~both
        hc[only_absorb] = absorb_hc[absorb_idx[only_absorb]]
        hm[only_absorb] = absorb_hm[absorb_idx[only_absorb]]
        if np.any(both):
            base_hc = keep_hc[keep_idx[both]]
            new_hc = np.minimum(base_hc + absorb_hc[absorb_idx[both]], SUTURE_ACCRETION_MAX_HC_M)
            hc[both] = new_hc
            # Also capped (issue #161) -- see _redistribute_accreted_column's own note on why
            # a large new_hc/base_hc ratio here (a thin column absorbing a thick one) can
            # otherwise carry Hm past its own ceiling even though Hc's own is respected.
            hm[both] = np.minimum(keep_hm[keep_idx[both]] * (new_hc / base_hc), lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)

        theta_owned = theta_candidates[owned]
        lines.append(
            ElevationLine(
                phi=phi,
                theta=theta_owned,
                elevation=np.zeros(len(theta_owned)),
                crustal_thickness_m=hc[owned],
                mantle_lithosphere_thickness_m=hm[owned],
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
    plates: list[LithospherePlate] = []
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

        def elevation_at(world_pts: np.ndarray) -> np.ndarray:
            return np.zeros(len(world_pts))  # placeholder; synced from Hc/Hm below

        lines = build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)
        hc_lines = []
        for line in lines:
            world_pts = line.world_xyz(frame)
            # Also upper-clipped (issue #161): the noise term alone can occasionally seed a
            # continental node above MAX_CRUSTAL_THICKNESS_M (confirmed directly -- generation
            # noise landed at 85,349 m on one seed), which is generation-time noise, not real
            # tectonic mass, so it's simply clipped rather than conserved/redistributed.
            hc = np.clip(hc_at(world_pts), lithosphere.MIN_CRUSTAL_THICKNESS_M, lithosphere.MAX_CRUSTAL_THICKNESS_M)
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


def new_plate(
    plate_id: int,
    frame: np.ndarray,
    crust_type: str,
    spacing_rad: float,
    seed: int,
    is_owned=None,
    node_is_continental=None,
) -> LithospherePlate:
    """A brand-new `LithospherePlate` seeded with reference Hc/Hm plus the same composite
    relief field `generate_plates` uses (see `terrain_noise.py`) -- the v2 analogue of
    `plates.generate_plates`' own per-plate initial-line construction. Keyed off
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
    the plate's own `crust_type` is the majority of what actually got built
    (`elevation_lines.majority_crust_type`) rather than the `crust_type` argument verbatim --
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

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        return np.zeros(len(world_pts))

    lines = build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)
    hc_lines = []
    for line in lines:
        world_pts = line.world_xyz(frame)
        is_continental = node_is_continental(world_pts)
        # Also upper-clipped (issue #161) -- see generate_plates' own note above.
        hc = np.clip(hc_at(world_pts, is_continental), lithosphere.MIN_CRUSTAL_THICKNESS_M, lithosphere.MAX_CRUSTAL_THICKNESS_M)
        hm = np.where(is_continental, hm0_continental, hm0_oceanic)
        code = np.where(is_continental, CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC).astype(np.int8)
        hc_lines.append(line.replace(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm, crust_type_code=code))
    majority = majority_crust_type(hc_lines, crust_type)
    plate = LithospherePlate(plate_id=plate_id, frame=frame, crust_type=majority, lines=hc_lines)
    lithosphere.sync_plate_elevation(plate)
    return plate

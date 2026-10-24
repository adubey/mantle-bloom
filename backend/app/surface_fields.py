"""Representation-neutral metadata for every persistent plate-surface field."""

from dataclasses import dataclass
from enum import Enum

import numpy as np


class RemapClass(str, Enum):
    EXTENSIVE = "extensive"
    INTENSIVE = "intensive"
    CATEGORICAL = "categorical"
    BOOLEAN_PROVENANCE = "boolean_provenance"
    COUNTDOWN = "countdown"
    CLOCK = "clock"
    HISTORY = "history"
    WRITE_ONCE_HISTORY = "write_once_history"
    # A dated mark a merged or coarsened cell keeps only if most of its area carried one,
    # with the latest date among those (`majority_history`; RFC 0009 belts).
    MAJORITY_HISTORY = "majority_history"
    DERIVED = "derived"


@dataclass(frozen=True)
class SurfaceField:
    dtype: np.dtype
    default: float | int | bool
    remap_class: RemapClass
    sentinel: float | int | None = None
    coupled_to: tuple[str, ...] = ()


def _field(dtype, default, remap_class, *, sentinel=None, coupled_to=()):
    return SurfaceField(np.dtype(dtype), default, remap_class, sentinel, coupled_to)


# `craton_formed_years`' "no craton" sentinel. Not 0.0: a column quiet from year zero
# legitimately forms with that date. Finite, so the surface audits' finiteness checks hold, and
# far beyond any reachable simulation year.
CRATON_UNFORMED_YEARS = 1.0e18
# `belt_formed_years`' "never a belt" sentinel (RFC 0009). Below any real date, generated
# belts' negative ones included, and finite for the same reason as CRATON_UNFORMED_YEARS.
BELT_UNMARKED_YEARS = -1.0e18

# `channel_reference_elevation_m`'s "no record yet" value. It's finite for the same reason as
# CRATON_UNFORMED_YEARS. Erosion treats it, and any area-weighted blend that includes it, as no
# record (see CHANNEL_REFERENCE_VALID_BELOW_M).
CHANNEL_REFERENCE_UNSET_M = 1.0e18
# A reference counts as recorded only below this, far above any real elevation. An
# area-weighted blend that includes the unset sentinel is then also treated as unset, rather
# than as a hugely negative rise.
CHANNEL_REFERENCE_VALID_BELOW_M = 1.0e6


# Phase 3 owns the algorithms that consume these classes. Phase 1 centralizes the policy so
# adding a persistent field without declaring its transfer semantics fails a contract test.
SURFACE_FIELDS: dict[str, SurfaceField] = {
    "elevation": _field(float, 0.0, RemapClass.INTENSIVE, coupled_to=("crustal_thickness_m", "mantle_lithosphere_thickness_m")),
    "channel_depth": _field(float, 0.0, RemapClass.INTENSIVE),
    "channel_width": _field(float, 0.0, RemapClass.INTENSIVE, coupled_to=("channel_depth",)),
    "lake_depth": _field(float, 0.0, RemapClass.EXTENSIVE),
    "glacier_depth": _field(float, 0.0, RemapClass.EXTENSIVE),
    "silt_depth": _field(float, 0.0, RemapClass.EXTENSIVE),
    "is_volcano": _field(bool, False, RemapClass.BOOLEAN_PROVENANCE),
    "volcano_active_years_remaining": _field(float, 0.0, RemapClass.COUNTDOWN),
    "soil_depth": _field(float, 0.0, RemapClass.EXTENSIVE),
    "soil_mineral_content": _field(float, 0.0, RemapClass.INTENSIVE, coupled_to=("soil_depth",)),
    "soil_organic_content": _field(float, 0.0, RemapClass.INTENSIVE, coupled_to=("soil_depth",)),
    "coal_deposit_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    "oil_gas_deposit_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    "mineral_deposit_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    "divergent_age_myr": _field(float, 0.0, RemapClass.CLOCK),
    "elev_change_reason": _field(float, 0.0, RemapClass.CATEGORICAL),
    "overlap_onset_years": _field(float, 0.0, RemapClass.HISTORY, sentinel=0.0),
    "node_created_years": _field(float, -1.0, RemapClass.WRITE_ONCE_HISTORY, sentinel=-1.0),
    "crustal_thickness_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    "mantle_lithosphere_thickness_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    "crust_type_code": _field(np.int8, 0, RemapClass.CATEGORICAL),
    # Thickness-equivalent continental-derived material.  Extensive so quad refinement,
    # coarsening, merge and partition conserve its area-integrated volume independently of
    # the receiving cell's binary crust type.
    "continental_material_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    # Cratons (cratons.py). The cratonic share of Hc is extensive so every remap conserves its
    # volume; its formation date keeps the oldest contributor's (CRATON_UNFORMED_YEARS: no
    # craton); the formation clock blends.
    "craton_crust_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    "craton_formed_years": _field(float, CRATON_UNFORMED_YEARS, RemapClass.HISTORY, sentinel=CRATON_UNFORMED_YEARS),
    "stable_continental_myr": _field(float, 0.0, RemapClass.CLOCK),
    # When the column last became a belt -- welded in a collision or stretched by a rift --
    # so it never forms craton (cratons.mark_belt, RFC 0009). Majority by area, not HISTORY:
    # keeping any contributor's mark spread belts onto craton cells in G1 (9 M km^2 on seed 7
    # by 500 Myr).
    "belt_formed_years": _field(float, BELT_UNMARKED_YEARS, RemapClass.MAJORITY_HISTORY, sentinel=BELT_UNMARKED_YEARS),
    # The (<= 0) isostatic depression the current ice load has applied to `elevation` (see
    # lithosphere.ice_load_deflection). It records what is baked into `elevation`, so it
    # remaps the way `elevation` does rather than the way `glacier_depth` does.
    "ice_load_deflection_m": _field(float, 0.0, RemapClass.INTENSIVE),
    # RFC 0015 gate G0 prototype (#554; gate branch only): last step's basal frictional heat
    # (W m^-2), and the deflection arm (c) holds the mass-balance surface at.
    "g0_frictional_heat_w_m2": _field(float, 0.0, RemapClass.INTENSIVE),
    "g0_hold_deflection_m": _field(float, 0.0, RemapClass.INTENSIVE),
    # Each node's elevation right after last step's erosion (issue #297). Next step's erosion
    # reads any rise since then as uplift, which wears down channel_depth. Remaps like
    # elevation.
    "channel_reference_elevation_m": _field(float, CHANNEL_REFERENCE_UNSET_M, RemapClass.INTENSIVE),
    # Breach notch depth (breaching.py, issue #297): sub-cell relief cut through a pit's rim.
    # Only the passage elevation reads it. It fades like channel_depth, and coarsens and merges
    # by max like channel_depth (sparse_quad_patch, quad_merge), so a notch isn't averaged
    # shut.
    "breach_notch_depth_m": _field(float, 0.0, RemapClass.INTENSIVE),
    # Anatexis (orogeny.py, quad plates only): how far the Moho lags below its steady-state
    # temperature (C; 0, the default, is steady state), and the refractory melt residue at
    # the base of the crust (an extensive share of Hc, like the craton's).
    "moho_thermal_lag_c": _field(float, 0.0, RemapClass.INTENSIVE),
    "restite_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    # Mobile cover (erosion.py, issue #297 phase 2): the loose sediment/regolith at the top of
    # the crust column -- a share of Hc, so extensive like the craton's -- and the
    # continental-derived share of it (a share of `continental_material_m`).
    "mobile_cover_m": _field(float, 0.0, RemapClass.EXTENSIVE),
    "mobile_cover_continental_m": _field(float, 0.0, RemapClass.EXTENSIVE),
}


def majority_history(
    cell: np.ndarray,
    values: np.ndarray,
    weights: np.ndarray,
    count: int,
    sentinel: float,
    own: np.ndarray | None = None,
    own_weights: np.ndarray | None = None,
) -> np.ndarray:
    """The MAJORITY_HISTORY remap: per output cell (0..`count`-1), the latest of the marked
    `values` (those not `sentinel`) mapped to it by `cell`, if the marked contributors hold
    more than half its total weight (area); otherwise `sentinel`. `own` and `own_weights`
    are an existing cell's own value and weight, voting alongside what it takes in (a merge
    suture)."""
    cell = np.asarray(cell, dtype=np.int64)
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    marked = values != sentinel
    total = np.bincount(cell, weights=weights, minlength=count)
    marked_weight = np.bincount(cell[marked], weights=weights[marked], minlength=count)
    latest = np.full(count, -np.inf)
    np.maximum.at(latest, cell[marked], values[marked])
    if own is not None:
        own = np.asarray(own, dtype=float)
        own_weights = np.asarray(own_weights, dtype=float)
        own_marked = own != sentinel
        total = total + own_weights
        marked_weight = marked_weight + np.where(own_marked, own_weights, 0.0)
        latest = np.where(own_marked, np.maximum(latest, own), latest)
    # Relative tolerance: equal areas summed from sub-cells differ in the last ulp, and an
    # even split is not "most".
    return np.where(marked_weight > 0.5 * total * (1.0 + 1e-9), latest, sentinel)

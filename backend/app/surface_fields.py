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

# `channel_reference_elevation_m`'s "no record yet" value. It's finite for the same reason as
# CRATON_UNFORMED_YEARS. Because it is so high, `max(0, elevation - reference)` reads as no
# uplift, both for an unset node and for any area-weighted blend that includes one.
CHANNEL_REFERENCE_UNSET_M = 1.0e18


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
    # The (<= 0) isostatic depression the current ice load has applied to `elevation` (see
    # lithosphere.ice_load_deflection). It records what is baked into `elevation`, so it
    # remaps the way `elevation` does rather than the way `glacier_depth` does.
    "ice_load_deflection_m": _field(float, 0.0, RemapClass.INTENSIVE),
    # Each node's elevation right after last step's erosion (issue #297). Next step's erosion
    # reads any rise since then as uplift, which wears down channel_depth. Remaps like
    # elevation.
    "channel_reference_elevation_m": _field(float, CHANNEL_REFERENCE_UNSET_M, RemapClass.INTENSIVE),
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

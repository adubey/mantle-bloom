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
}

"""Persistent, area-weighted continental-material provenance accounting.

``continental_material_m`` is a thickness-equivalent tracer on every terrain node.  It is
independent of the host column's crust type: continental sediment can therefore remain
continental material after it is deposited on an oceanic plate or after a plate is retyped.
The world ledger stores only sources and terminal sinks; live surface inventories are always
recomputed from the node tracer, avoiding a second mutable copy of the same quantity.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, TypedDict

import numpy as np

from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad

if TYPE_CHECKING:
    from .world import World


LedgerAccount = Literal[
    "initial_continental_m3",
    "juvenile_additions_m3",
    "accreted_thickened_m3",
    "delaminated_lower_crust_m3",
    "deeply_subducted_m3",
    "remelted_relaminated_returns_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
]


class ContinentalMaterialLedger(TypedDict):
    initial_continental_m3: float
    juvenile_additions_m3: float
    accreted_thickened_m3: float
    delaminated_lower_crust_m3: float
    deeply_subducted_m3: float
    remelted_relaminated_returns_m3: float
    numerical_unplaced_m3: float
    # Continental sediment the `ocean_deposition_multiplier` knob (< 1) declines to settle --
    # a deliberate user-tuned shelf-starving sink, kept apart from numerical clipping.
    discarded_marine_sediment_m3: float


LEDGER_KEYS: tuple[LedgerAccount, ...] = (
    "initial_continental_m3",
    "juvenile_additions_m3",
    "accreted_thickened_m3",
    "delaminated_lower_crust_m3",
    "deeply_subducted_m3",
    "remelted_relaminated_returns_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
)


def empty_ledger() -> ContinentalMaterialLedger:
    return ContinentalMaterialLedger(**{key: 0.0 for key in LEDGER_KEYS})


def ensure_initialized(world: "World") -> None:
    """Backfill old saves and seed untracked line or quad surfaces."""
    if not hasattr(world, "continental_material_ledger"):
        world.continental_material_ledger = empty_ledger()
    for key in LEDGER_KEYS:
        world.continental_material_ledger.setdefault(key, 0.0)

    ledger_is_new = not any(world.continental_material_ledger.values())
    tracer_is_empty = not any(
        np.any(plate.collect("continental_material_m"))
        for plate in world.plates
    )
    if ledger_is_new and tracer_is_empty:
        for plate in world.plates:
            codes = plate.collect("crust_type_code")
            is_continental = effective_is_continental_from_codes(
                codes, plate.crust_type == "continental"
            )
            plate.set_fields_on_plate(
                continental_material_m=np.where(
                    is_continental, plate.collect("crustal_thickness_m"), 0.0
                )
            )
    if ledger_is_new and world.continental_material_ledger["initial_continental_m3"] == 0.0:
        world.continental_material_ledger["initial_continental_m3"] = surface_volume_m3(world)


def surface_volume_m3(world: "World") -> float:
    spacing = line_spacing_rad(world.node_density)
    return sum(
        float(np.sum(plate.collect("continental_material_m") * plate.accounting_areas_m2(spacing)))
        for plate in world.plates
    )


def record(world: "World", account: LedgerAccount, volume_m3: float) -> None:
    """Add a non-negative volume to one persisted account.

    Terrain processes should call this at the exact cap or transfer site.  Keeping the API in
    volume units prevents callers from accidentally summing thicknesses without cell area.
    """
    ensure_initialized(world)
    if not np.isfinite(volume_m3) or volume_m3 < 0.0:
        raise ValueError(f"ledger volume must be finite and non-negative, got {volume_m3!r}")
    if account not in LEDGER_KEYS:
        raise KeyError(f"unknown continental-material ledger account: {account}")
    world.continental_material_ledger[account] += float(volume_m3)


def inventories(world: "World") -> dict[str, float]:
    """Return current live inventories plus the persisted source/sink accounts."""
    ensure_initialized(world)
    spacing = line_spacing_rad(world.node_density)
    surface = 0.0
    sediment_on_ocean = 0.0
    for plate in world.plates:
        material = np.clip(
            plate.collect("continental_material_m"), 0.0, plate.collect("crustal_thickness_m")
        )
        areas = plate.accounting_areas_m2(spacing)
        surface += float(np.sum(material * areas))
        host_continental = effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        sediment_on_ocean += float(np.sum(material[~host_continental] * areas[~host_continental]))
    result = dict(world.continental_material_ledger)
    result["surface_continental_derived_m3"] = surface
    result["continental_sediment_on_oceanic_hosts_m3"] = sediment_on_ocean
    result["balance_error_m3"] = balance_error_m3(world, surface=surface)
    return result


def balance_error_m3(world: "World", *, surface: float | None = None) -> float:
    ledger = world.continental_material_ledger
    if surface is None:
        surface = surface_volume_m3(world)
    sources = ledger["initial_continental_m3"] + ledger["juvenile_additions_m3"]
    sinks_and_live = (
        surface
        + ledger["delaminated_lower_crust_m3"]
        + ledger["deeply_subducted_m3"]
        + ledger["numerical_unplaced_m3"]
        + ledger["discarded_marine_sediment_m3"]
    )
    return sinks_and_live - sources


def assert_closed(world: "World", *, relative_tolerance: float = 1e-10) -> None:
    """Fail a diagnostic check when instrumented sources, live material and sinks differ."""
    ensure_initialized(world)
    ledger = world.continental_material_ledger
    error = balance_error_m3(world)
    tolerance = max(1.0, relative_tolerance * max(ledger["initial_continental_m3"], 1.0))
    if abs(error) > tolerance:
        raise AssertionError(
            f"continental material ledger does not close: error={error:.6g} m^3, "
            f"tolerance={tolerance:.6g} m^3"
        )

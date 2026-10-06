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
    "collision_subducted_m3",
    "deeply_subducted_m3",
    "remelted_relaminated_returns_m3",
    "rift_thinned_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
    "overloaded_root_delaminated_m3",
    "topology_removed_m3",
]


class ContinentalMaterialLedger(TypedDict):
    initial_continental_m3: float
    juvenile_additions_m3: float
    accreted_thickened_m3: float
    delaminated_lower_crust_m3: float
    # Suture crust no receiver anywhere could hold (quad_tectonics._accrete_onto_survivors'
    # terminal remainder), which goes down with the consumed plate (issue #276). Kept apart
    # from trench subduction in `deeply_subducted_m3`.
    collision_subducted_m3: float
    deeply_subducted_m3: float
    remelted_relaminated_returns_m3: float
    # Rifting: stretch thinning a column's footprint can't hold on the fixed-area node, plus
    # continental columns that melt through and reset to oceanic ridge crust.
    rift_thinned_m3: float
    numerical_unplaced_m3: float
    # Continental sediment the `ocean_deposition_multiplier` knob (< 1) declines to settle --
    # a deliberate user-tuned shelf-starving sink, kept apart from numerical clipping.
    discarded_marine_sediment_m3: float
    # Sediment a column already at the Hc cap was loaded with after erosion's overflow carry
    # found no receiver with room in reach (issue #288): the overloaded root sheds it instead.
    # Kept apart from suture-accretion delamination and from numerical clipping.
    overloaded_root_delaminated_m3: float
    # Continental material on stranded fragments a plate's defragmentation drops, and on
    # plates removed with no territory left (merge_split.py) -- geometric cleanup, not
    # physics; kept apart so a save shows how much land it costs (issue #276).
    topology_removed_m3: float


LEDGER_KEYS: tuple[LedgerAccount, ...] = (
    "initial_continental_m3",
    "juvenile_additions_m3",
    "accreted_thickened_m3",
    "delaminated_lower_crust_m3",
    "collision_subducted_m3",
    "deeply_subducted_m3",
    "remelted_relaminated_returns_m3",
    "rift_thinned_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
    "overloaded_root_delaminated_m3",
    "topology_removed_m3",
)


def empty_ledger() -> ContinentalMaterialLedger:
    return ContinentalMaterialLedger(**{key: 0.0 for key in LEDGER_KEYS})


def ensure_initialized(world: "World") -> None:
    """Backfill old saves and seed untracked line or quad surfaces."""
    ledger_is_new = not hasattr(world, "continental_material_ledger") or not world.continental_material_ledger
    if not hasattr(world, "continental_material_ledger"):
        world.continental_material_ledger = {}
    for key in LEDGER_KEYS:
        world.continental_material_ledger.setdefault(key, 0.0)

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
    if ledger_is_new:
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


def add_material_thickness(
    world: "World",
    plate,
    delta_m: np.ndarray,
    account: LedgerAccount,
    *,
    eligible: np.ndarray | None = None,
) -> float:
    """Add provenance thickness to existing nodes and book its exact volume.

    ``delta_m`` is node-aligned and may contain zeros; negative values are rejected.  This is
    the common write path for juvenile volcanism and recycled magma returns, keeping the
    per-node tracer and the persisted volume account impossible to update separately.
    """
    ensure_initialized(world)
    delta = np.asarray(delta_m, dtype=float)
    if delta.shape != (plate.node_count(),):
        raise ValueError(f"material thickness must have shape ({plate.node_count()},), got {delta.shape}")
    if np.any(~np.isfinite(delta)) or np.any(delta < 0.0):
        raise ValueError("material thickness must be finite and non-negative")
    if eligible is not None:
        delta = np.where(np.asarray(eligible, dtype=bool), delta, 0.0)
    if not np.any(delta):
        return 0.0
    material = plate.collect("continental_material_m") + delta
    plate.set_fields_on_plate(continental_material_m=material)
    areas = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
    volume = float(np.dot(delta, areas))
    record(world, account, volume)
    return volume


def inventories(world: "World") -> dict[str, float]:
    """Return current live inventories plus the persisted source/sink accounts."""
    ensure_initialized(world)
    spacing = line_spacing_rad(world.node_density)
    surface = 0.0
    sediment_on_ocean = 0.0
    for plate in world.plates:
        material = plate.collect("continental_material_m")
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
    sources = (
        ledger["initial_continental_m3"]
        + ledger["juvenile_additions_m3"]
        + ledger["accreted_thickened_m3"]
        + ledger["remelted_relaminated_returns_m3"]
    )
    sinks_and_live = (
        surface
        + ledger["delaminated_lower_crust_m3"]
        + ledger["collision_subducted_m3"]
        + ledger["deeply_subducted_m3"]
        + ledger["rift_thinned_m3"]
        + ledger["numerical_unplaced_m3"]
        + ledger["discarded_marine_sediment_m3"]
        + ledger["overloaded_root_delaminated_m3"]
        + ledger["topology_removed_m3"]
    )
    return sinks_and_live - sources


def assert_closed(world: "World", *, relative_tolerance: float = 1e-10) -> None:
    """Fail a diagnostic check when instrumented sources, live material and sinks differ."""
    ensure_initialized(world)
    for plate in world.plates:
        material = plate.collect("continental_material_m")
        hc = plate.collect("crustal_thickness_m")
        if np.any(~np.isfinite(material)) or np.any(material < -1e-9) or np.any(material > hc + 1e-9):
            raise AssertionError("continental material must be finite and remain within [0, Hc]")
    ledger = world.continental_material_ledger
    error = balance_error_m3(world)
    tolerance = max(1.0, relative_tolerance * max(ledger["initial_continental_m3"], 1.0))
    if abs(error) > tolerance:
        raise AssertionError(
            f"continental material ledger does not close: error={error:.6g} m^3, "
            f"tolerance={tolerance:.6g} m^3"
        )

"""Persistent accounting for the mobile cover (issue #297 phase 2).

``mobile_cover_m`` is the loose sediment/regolith share of each column's Hc, at the top of the
column, and ``mobile_cover_continental_m`` the continental-derived share of that (see
erosion.MOBILE_COVER_*). Erosion fills and strips it. Tectonics moves or ends it:

- rift stretching thins it with its column, and spreads what it removes over the rifted cells
  the stretch covers (`quad_tectonics._open_rift`); a failed rift thins it in place
  (`thin_with_column`, ``rift_thinned_m3``); a column that melts through is reset to fresh
  magmatic crust, so its cover is gone (``rift_reset_m3``);
- orogenic collapse and ductile flow carry it with the crust they move, in proportion;
- cells consumed at a trench take it down with them (``subducted_m3``), and suture donors and
  a relocated terrane's displaced columns metamorphose it into the accreted crust
  (``accreted_m3``);
- erupted lava buries it, consolidating it into substrate (``volcanic_buried_m3``), in
  proportion to how much of the cell the lava covers (VOLCANIC_SEAL_THICKNESS_M);
- stranded fragments a plate's defragmentation drops take theirs with them (``stranded_m3``).

Shortening, underplating, anatexis, fault relief and the column caps change Hc at depth or
don't move material in this model, and leave it alone; what any of them leaves above its
column is clipped at the next erosion step (``clipped_m3``). Quad merge, refinement and
coarsening remap it as an extensive field, exactly. Strike-slip advection
(`faults._apply_plate_fault_shear`) copies every field from its nearest upstream node, which
conserves no extensive field -- Hc and the continental tracer included -- so its net change is
booked, signed, as ``fault_advection_m3``.

Like continental_ledger.py, the world stores only sources and sinks; the live inventory is
recomputed from the nodes, so ``balance_error_m3`` is zero when every change is booked. The
continental share is a share of ``continental_material_m``, whose own ledger covers it.
Booking is complete: every quad topology change books what it removes."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..elevation_lines import line_spacing_rad

if TYPE_CHECKING:
    from ..plates import Plate
    from ..world import World

FIELDS = ("mobile_cover_m", "mobile_cover_continental_m")
# Erupted lava seals the cover under it, but a cell-mean addition is not a uniform sheet: the
# thin tail of a volcanic-plain apron is lava over part of the cell's footprint. A cell whose
# surface rises this much (a few stacked flood-basalt flows, each typically 5-30 m) counts as
# fully sealed; less seals that fraction of its cover. Measured as the eruption's elevation
# gain. Under the model's isostasy, lava stacked on the surface raises it by only ~15% of its
# thickness (lithosphere.back_elevation_gain backs each metre of relief with ~6 m of
# continental Hc), so Hc gain bounds the lava thickness from above and elevation gain from
# below. Elevation gain is the better proxy because large igneous provinces emplace roughly
# 5-10x more intrusive than extrusive volume, which puts the surface lava near Hc gain / 6.
VOLCANIC_SEAL_THICKNESS_M = 20.0
SOURCES = ("initial_m3", "deposited_m3")
SINKS = (
    "entrained_m3",
    "consolidated_m3",
    "clipped_m3",
    "rift_thinned_m3",
    "rift_reset_m3",
    "subducted_m3",
    "accreted_m3",
    "volcanic_buried_m3",
    "fault_advection_m3",
    "stranded_m3",
)
# Accounts that may be negative: round-off clipped back up, and advection that gained volume.
SIGNED = ("clipped_m3", "fault_advection_m3")
ACCOUNTS = SOURCES + SINKS


def surface_volume_m3(world: "World") -> float:
    spacing = line_spacing_rad(world.node_density)
    return sum(float(np.dot(p.collect("mobile_cover_m"), p.accounting_areas_m2(spacing))) for p in world.plates)


def ensure_ledger(world: "World") -> dict[str, float]:
    """Backfill the accounts; the first time, seed ``initial_m3`` with whatever cover the
    world already holds (none, on any world from before the cover existed)."""
    ledger = getattr(world, "mobile_cover_ledger", None)
    if not isinstance(ledger, dict):
        ledger = world.mobile_cover_ledger = {}
    if "initial_m3" not in ledger:
        ledger["initial_m3"] = surface_volume_m3(world)
    for key in ACCOUNTS:
        ledger.setdefault(key, 0.0)
    return ledger


def record(world: "World | None", account: str, volume_m3: float) -> None:
    """Book `volume_m3` into one account. Only the SIGNED accounts may be negative."""
    if world is None:
        return
    if account not in ACCOUNTS or account == "initial_m3":
        raise KeyError(f"unknown mobile cover account: {account}")
    volume_m3 = float(volume_m3)
    if not np.isfinite(volume_m3) or (volume_m3 < 0.0 and account not in SIGNED):
        raise ValueError(f"mobile cover volume must be finite and non-negative, got {volume_m3!r}")
    ensure_ledger(world)[account] += volume_m3


def balance_error_m3(world: "World") -> float:
    ledger = ensure_ledger(world)
    expected = sum(ledger[key] for key in SOURCES) - sum(ledger[key] for key in SINKS)
    return surface_volume_m3(world) - expected


def book_removed(world: "World | None", plate: "Plate", mask: np.ndarray, account: str) -> None:
    """Book the cover on the `mask` cells, which the caller is about to remove."""
    if world is None or not np.any(mask):
        return
    areas = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
    record(world, account, float(np.dot(plate.collect("mobile_cover_m")[mask], areas[mask])))


def thin_with_column(world: "World | None", plate: "Plate", hc_before: np.ndarray, account: str) -> None:
    """After a whole-column thinning that kept node order (a failed rift), thin each node's
    cover by its Hc ratio, as `cratons.thin_with_column` does the craton, and book the loss
    to `account`. Thickening leaves it unchanged."""
    ratio = np.clip(
        np.divide(plate.collect("crustal_thickness_m"), hc_before, out=np.ones(len(hc_before)), where=hc_before > 0.0),
        0.0,
        1.0,
    )
    end(world, plate, 1.0 - ratio, account)


def end(world: "World | None", plate: "Plate", fraction: np.ndarray, account: str) -> None:
    """Book and clear `fraction` (per cell, in [0, 1]; a mask ends all of it) of each cell's
    cover, which stays in its column: it is buried, consolidated or reset into the substrate
    under it. Hc is unchanged."""
    fraction = np.clip(np.asarray(fraction, dtype=float), 0.0, 1.0)
    if not np.any(fraction > 0.0):
        return
    cover = plate.collect("mobile_cover_m")
    if world is not None:
        areas = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
        record(world, account, float(np.dot(cover * fraction, areas)))
    plate.set_fields_on_plate(**{name: plate.collect(name) * (1.0 - fraction) for name in FIELDS})

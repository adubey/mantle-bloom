"""Shared buoyant-crust terrane machinery (RFC 0007, issue #580).

Docking, buoyant-crust tracers (`plateau_crust_m`, RFC 0006's `arc_crust_m`) and
their age moment, and the buoyant-crust ledger -- built first by RFC 0007, consumed
as-is by RFC 0006. This module owns no plates or worlds: patch geometry arrives as
explicit index sets, tracer shares as arrays, and the ledger as a plain dict, so
every helper is a pure function of its inputs plus the one in-place audit.

Field ownership: `plateau_crust_m` and the shared age moment
`buoyant_crust_time_m_yr` are registered here (see `surface_fields.py`);
`arc_crust_m` is RFC 0006's registration. `World.buoyant_crust_ledger` is keyed by
source (`arc`, `plume_head`, `plume_tail`).
"""

from __future__ import annotations

import hashlib
import math

import numpy as np

SOURCES = ("arc", "plume_head", "plume_tail")
# Reservoirs that partition emplaced buoyant crust; the balance closes over these
# plus live tracers (see balance_residual). `docked_m3` counts docking events and
# sits outside the balance because docked crust stays live; `clipped_by_phase`
# records the per-phase share of `clipped_m3` outside it for the same reason.
BALANCE_ACCOUNTS = (
    "emplaced_m3",
    "thinned_m3",
    "reset_m3",
    "subducted_m3",
    "converted_m3",
    "foundered_m3",
    "clipped_m3",
    "unattributed_m3",
)
DOCKED_ACCOUNT = "docked_m3"
CLIPPED_BY_PHASE = "clipped_by_phase"

# Age-dependent docking threshold, in excess meters: a young Caribbean-thickness
# plateau docks while an old one subducts; Ontong Java (>= 23.5 km excess) resisted
# 100 Myr after forming (issue #580).
DOCK_E_YOUNG_M = 5_000.0
DOCK_E_OLD_M = 15_000.0
DOCK_COOLING_MYR = 20.0
DOCK_WIDTH_M = 3_000.0
# Ontong Java lost lower crust at the trench; no fraction is published, so this
# design choice is swept 0-0.4 in sensitivity runs (issue #580).
DOCK_LOWER_CRUST_LOSS = 0.2
# Budgets count accreted plateaus as continental growth though mafic; the
# sensitivity case halves it (issue #580).
PLATEAU_PRESERVED_FRACTION = 0.9


def plate_tracers(plate) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A plate's buoyant-crust tracers: (plateau, arc, age moment), float copies.

    `arc_crust_m` is RFC 0006's registration and may not exist yet; until it
    lands this returns zeros for the arc channel so every writer stays exact and
    generic over both sources (issue #581).
    """
    plateau = np.asarray(plate.collect("plateau_crust_m"), dtype=float)
    moment = np.asarray(plate.collect("buoyant_crust_time_m_yr"), dtype=float)
    try:
        arc = np.asarray(plate.collect("arc_crust_m"), dtype=float)
    except AttributeError:
        arc = np.zeros(plate.node_count(), dtype=float)
    return plateau, arc, moment


def live_volumes_m3(
    plateau_m: np.ndarray, arc_m: np.ndarray, areas_m2: np.ndarray
) -> dict[str, float]:
    """Live buoyant-crust volumes by ledger source, for `balance_residual`.

    The plateau channel is plume-head crust (tails emplace ordinary-thickness
    melt with no tracer, issue #581); the arc channel is RFC 0006's.
    """
    areas = np.asarray(areas_m2, dtype=float)
    return {
        "arc": float(np.sum(np.asarray(arc_m, dtype=float) * areas)),
        "plume_head": float(np.sum(np.asarray(plateau_m, dtype=float) * areas)),
        "plume_tail": 0.0,
    }


def world_live_volumes_m3(world) -> dict[str, float]:
    """Live buoyant-crust volumes by source, summed over every plate."""
    from .elevation_lines import line_spacing_rad

    spacing = line_spacing_rad(world.node_density)
    total = {"arc": 0.0, "plume_head": 0.0, "plume_tail": 0.0}
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        plateau, arc, _ = plate_tracers(plate)
        if not np.any(plateau) and not np.any(arc):
            continue
        for source, volume in live_volumes_m3(
            plateau, arc, plate.accounting_areas_m2(spacing)
        ).items():
            total[source] += volume
    return total


def audit_world(world, phase: str, reference_hc_m: float = 0.0) -> float:
    """Hold every plate's tracers to their column bounds, booking corrections.

    Phase-bracket backstop for every writer the lifecycle table routes through
    the audit (erosion, column caps, crust-type changes, anything a writer
    missed): `0 <= plateau + arc <= E` on oceanic-coded cells, 0 on continental
    ones. Returns the total booked volume. Idempotent; a no-op without tracers.
    """
    from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad

    ledger = ensure_ledger(world)
    spacing = line_spacing_rad(world.node_density)
    booked = 0.0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        plateau, arc, moment = plate_tracers(plate)
        if not np.any(plateau) and not np.any(arc):
            continue
        areas = plate.accounting_areas_m2(spacing)
        continental = effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        hc = plate.collect("crustal_thickness_m")
        # One call audits the plateau + arc total together (the bound is on the
        # sum). The arc channel is all zeros until RFC 0006 registers it, so the
        # whole correction is plateau crust today; when 0006 lands it owns the
        # per-source split of a mixed correction here.
        booked += apply_phase_audit(
            plateau,
            arc,
            moment,
            hc,
            areas,
            continental,
            ledger,
            "plume_head",
            phase,
            reference_hc_m,
        )
        plate.set_fields_on_plate(plateau_crust_m=plateau, buoyant_crust_time_m_yr=moment)
    return booked


def new_ledger() -> dict:
    """A fresh buoyant-crust ledger: zeroed balance accounts per source, with the
    outside-balance keys (`docked_m3`, per-phase clips) present so readers never
    depend on booking history (#580)."""
    return {
        source: {**{account: 0.0 for account in BALANCE_ACCOUNTS}, DOCKED_ACCOUNT: 0.0, CLIPPED_BY_PHASE: {}}
        for source in SOURCES
    }


def ensure_ledger(world) -> dict:
    """Create or repair the buoyant-crust ledger on a world (#580).

    Fresh and backfilled worlds start with an empty or partial dict; every booking
    path needs per-source accounts. Idempotent: existing values are preserved,
    only missing sources and accounts are filled.
    """
    ledger = getattr(world, "buoyant_crust_ledger", None)
    if not isinstance(ledger, dict):
        ledger = {}
        world.buoyant_crust_ledger = ledger
    for source in SOURCES:
        accounts = ledger.get(source)
        if not isinstance(accounts, dict):
            accounts = {}
            ledger[source] = accounts
        for account in BALANCE_ACCOUNTS:
            accounts.setdefault(account, 0.0)
        accounts.setdefault(DOCKED_ACCOUNT, 0.0)
        accounts.setdefault(CLIPPED_BY_PHASE, {})
    return ledger


def record(ledger: dict, source: str, account: str, volume_m3: float) -> float:
    """Add a non-negative volume to one source account, keeping the exact sum."""
    if source not in ledger or not isinstance(ledger[source], dict):
        raise KeyError(f"unknown buoyant-crust source: {source!r}")
    if account == DOCKED_ACCOUNT:
        raise ValueError("book dockings with record_docked, not record")
    if account == CLIPPED_BY_PHASE:
        raise ValueError("book per-phase clips with record_clipped_by_phase, not record")
    if account not in BALANCE_ACCOUNTS:
        raise KeyError(f"unknown buoyant-crust account: {account!r}")
    if not math.isfinite(volume_m3) or volume_m3 < 0.0:
        raise ValueError(f"ledger volume must be finite and non-negative, got {volume_m3!r}")
    ledger[source].setdefault(account, 0.0)
    ledger[source][account] += float(volume_m3)
    return float(volume_m3)


def record_docked(ledger: dict, source: str, volume_m3: float) -> float:
    """Count a docking event volume outside the balance (docked crust stays live)."""
    if source not in ledger:
        raise KeyError(f"unknown buoyant-crust source: {source!r}")
    if not math.isfinite(volume_m3) or volume_m3 < 0.0:
        raise ValueError(f"ledger volume must be finite and non-negative, got {volume_m3!r}")
    ledger[source][DOCKED_ACCOUNT] = ledger[source].get(DOCKED_ACCOUNT, 0.0) + float(volume_m3)
    return float(volume_m3)


def record_clipped_by_phase(ledger: dict, source: str, phase: str, volume_m3: float) -> float:
    """Book a phase audit's correction to `clipped_m3`, remembering its phase.

    Erosion's share is "eroded"; the per-phase breakdown sits outside the balance
    next to `docked_m3` while `clipped_m3` itself stays in it (#580).
    """
    if not math.isfinite(volume_m3) or volume_m3 < 0.0:
        raise ValueError(f"ledger volume must be finite and non-negative, got {volume_m3!r}")
    record(ledger, source, "clipped_m3", volume_m3)
    by_phase = ledger[source].setdefault(CLIPPED_BY_PHASE, {})
    by_phase[phase] = by_phase.get(phase, 0.0) + float(volume_m3)
    return float(volume_m3)


def threshold_excess_m(
    age_myr: float,
    e_young_m: float = DOCK_E_YOUNG_M,
    e_old_m: float = DOCK_E_OLD_M,
    cooling_myr: float = DOCK_COOLING_MYR,
) -> float:
    """Excess thickness a patch needs to resist at its age: young and thin may dock,
    old and thin subducts (#580)."""
    if not math.isfinite(age_myr):
        raise ValueError(f"age must be finite, got {age_myr!r}")
    age_myr = max(age_myr, 0.0)
    return e_young_m + (e_old_m - e_young_m) * (1.0 - math.exp(-age_myr / cooling_myr))


def docking_probability(
    excess_m: float,
    age_myr: float,
    e_young_m: float = DOCK_E_YOUNG_M,
    e_old_m: float = DOCK_E_OLD_M,
    cooling_myr: float = DOCK_COOLING_MYR,
    width_m: float = DOCK_WIDTH_M,
) -> float:
    """Docking probability: logistic in (excess - age threshold) / width (#580)."""
    if not math.isfinite(excess_m) or not math.isfinite(age_myr):
        raise ValueError(f"excess and age must be finite, got {excess_m!r}, {age_myr!r}")
    spread = max(width_m, 1e-9)
    odds = (excess_m - threshold_excess_m(age_myr, e_young_m, e_old_m, cooling_myr)) / spread
    odds = min(max(odds, -50.0), 50.0)
    return 1.0 / (1.0 + math.exp(-odds))


def docking_draw(seed: int, plate_id: int, cell_key: int, encounter_start_years: float, probability: float) -> bool:
    """One deterministic draw per patch encounter, held until it leaves the trench.

    Keyed on (seed, plate, the patch's lowest cell key, encounter start time), so
    the same encounter always draws the same outcome whatever the step partition.
    Pass the identical stored start value on re-entry: the key hashes its repr, so
    a bitwise-distinct float is a different encounter (#580).
    """
    if not 0.0 <= probability <= 1.0:
        raise ValueError(f"probability must be in [0, 1], got {probability!r}")
    digest = hashlib.sha256(f"{seed}/docking/{plate_id}/{cell_key}/{encounter_start_years!r}".encode()).digest()
    return int.from_bytes(digest[:8], "big") / 2**64 < probability


def group_patches(count: int, neighbor_pairs: np.ndarray, member_mask: np.ndarray) -> list[np.ndarray]:
    """Connected components of the member cells over explicit neighbor pairs.

    Pure union-find: no quad topology needed, so tests and the retreat wiring
    share the grouping (#580). Groups come out sorted by smallest member.
    """
    if len(np.asarray(member_mask).ravel()) != count:
        raise ValueError(f"mask length {len(np.asarray(member_mask).ravel())} != count {count}")
    member_mask = np.asarray(member_mask, dtype=bool)
    parent = list(range(count))

    def find(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    pairs = np.asarray(neighbor_pairs)
    if pairs.size and not np.all(np.equal(np.mod(pairs.ravel().astype(float), 1.0), 0.0)):
        raise ValueError("neighbor pairs must be integer node indices")
    pairs = pairs.astype(np.int64).reshape(-1, 2)
    for left, right in pairs:
        if 0 <= left < count and 0 <= right < count and member_mask[left] and member_mask[right]:
            left_root, right_root = find(int(left)), find(int(right))
            if left_root != right_root:
                parent[max(left_root, right_root)] = min(left_root, right_root)
    groups: dict[int, list[int]] = {}
    for node in range(count):
        if member_mask[node]:
            groups.setdefault(find(node), []).append(node)
    return [np.asarray(sorted(members), dtype=np.int64) for _, members in sorted(groups.items())]


def patch_excess_m(hc_m: np.ndarray, areas_m2: np.ndarray, mask: np.ndarray, reference_hc_m: float) -> float:
    """Area-weighted mean excess over the reference oceanic column (#580)."""
    from . import lithosphere

    reference = reference_hc_m if reference_hc_m > 0.0 else lithosphere.REFERENCE_HC_OCEANIC_M
    hc = np.asarray(hc_m, dtype=float)[mask]
    areas = np.asarray(areas_m2, dtype=float)[mask]
    total = float(np.sum(areas))
    if total <= 0.0:
        return 0.0
    return float(np.sum((hc - reference) * areas) / total)


def patch_age_years(
    plateau_m: np.ndarray,
    arc_m: np.ndarray,
    age_moment_m_yr: np.ndarray,
    mask: np.ndarray,
    now_years: float,
) -> float:
    """Volume-weighted buoyant-crust age from the age moment (#580)."""
    tracer = np.asarray(plateau_m, dtype=float)[mask] + np.asarray(arc_m, dtype=float)[mask]
    moment = np.asarray(age_moment_m_yr, dtype=float)[mask]
    total = float(np.sum(tracer))
    if total <= 0.0:
        return 0.0
    return float(now_years - np.sum(moment) / total)


def scale_tracers(
    plateau_m: np.ndarray, arc_m: np.ndarray, age_moment_m_yr: np.ndarray, ratio: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Order-preserving thinning: tracers and their age moment scale with Hc (#580)."""
    ratio = float(np.clip(ratio, 0.0, 1.0))
    return (
        np.asarray(plateau_m, dtype=float) * ratio,
        np.asarray(arc_m, dtype=float) * ratio,
        np.asarray(age_moment_m_yr, dtype=float) * ratio,
    )


def add_tracers(
    plateau_m: np.ndarray,
    arc_m: np.ndarray,
    age_moment_m_yr: np.ndarray,
    add_plateau_m: np.ndarray,
    add_arc_m: np.ndarray,
    add_age_moment_m_yr: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Carry tracer volumes with the Hc they ride on (transfers, accretion) (#580)."""
    return (
        np.asarray(plateau_m, dtype=float) + np.asarray(add_plateau_m, dtype=float),
        np.asarray(arc_m, dtype=float) + np.asarray(add_arc_m, dtype=float),
        np.asarray(age_moment_m_yr, dtype=float) + np.asarray(add_age_moment_m_yr, dtype=float),
    )


def clear_tracers(count: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Melt-through reset: buoyant shares vanish with their booking (#580)."""
    zeros = np.zeros(int(count), dtype=float)
    return zeros.copy(), zeros.copy(), zeros.copy()


def split_docking_loss(
    plateau_m: np.ndarray,
    arc_m: np.ndarray,
    age_moment_m_yr: np.ndarray,
    loss_fraction: float = DOCK_LOWER_CRUST_LOSS,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Partial docking loss: the lower-crust share goes down the trench (#580).

    Returns ((kept plateau, arc, moment), (lost plateau, arc, moment)); the caller
    books the lost share to `subducted_m3`.
    """
    kept = float(np.clip(1.0 - loss_fraction, 0.0, 1.0))
    plateau, arc, moment = (np.asarray(values, dtype=float) for values in (plateau_m, arc_m, age_moment_m_yr))
    return (plateau * kept, arc * kept, moment * kept), (plateau - plateau * kept, arc - arc * kept, moment - moment * kept)


def conversion_split(
    plateau_m: np.ndarray, arc_m: np.ndarray, preserved_fraction: float = PLATEAU_PRESERVED_FRACTION
) -> tuple[np.ndarray, np.ndarray]:
    """Docked-onto-continent conversion: preserved juvenile vs foundered mafic rock.

    Returns (juvenile thickness, foundered thickness); the caller books juvenile
    once, removes the foundered share from Hc, and clears the tracers (#580).
    """
    kept = float(np.clip(preserved_fraction, 0.0, 1.0))
    crust = np.asarray(plateau_m, dtype=float) + np.asarray(arc_m, dtype=float)
    return crust * kept, crust - crust * kept


def apply_phase_audit(
    plateau_m: np.ndarray,
    arc_m: np.ndarray,
    age_moment_m_yr: np.ndarray,
    hc_m: np.ndarray,
    areas_m2: np.ndarray,
    continental_host: np.ndarray,
    ledger: dict,
    source: str,
    phase: str,
    reference_hc_m: float = 0.0,
) -> float:
    """Hold tracers to `0 <= plateau + arc <= E` on oceanic hosts, 0 on continental
    ones, booking corrections to `clipped_m3` with the phase recorded (#580).

    Mutates the tracer and age-moment arrays in place (they must be float
    ndarrays, not lists: `asarray` would copy those and the write-back would be
    lost); returns the booked volume.
    """
    from . import lithosphere

    for name, values in (("plateau", plateau_m), ("arc", arc_m), ("age moment", age_moment_m_yr)):
        if not isinstance(values, np.ndarray) or values.dtype.kind != "f":
            raise TypeError(f"audit mutates float ndarrays in place, got {name}: {type(values)}")
    reference = reference_hc_m if reference_hc_m > 0.0 else lithosphere.REFERENCE_HC_OCEANIC_M
    hc = np.asarray(hc_m, dtype=float)
    areas = np.asarray(areas_m2, dtype=float)
    if not np.all(np.isfinite(hc)) or not np.all(np.isfinite(areas)):
        raise ValueError("audit inputs must be finite")
    excess = np.maximum(hc - reference, 0.0)
    host_continental = np.asarray(continental_host, dtype=bool)
    allowed = np.where(host_continental, 0.0, excess)
    total = plateau_m + arc_m
    over = np.maximum(total - allowed, 0.0)
    # Cap only: never inflate tracers up to the excess (#580).
    scale = np.minimum(np.divide(allowed, total, out=np.ones_like(total), where=total > 0.0), 1.0)
    plateau_m[:] = plateau_m * scale
    arc_m[:] = arc_m * scale
    age_moment_m_yr[:] = age_moment_m_yr * scale
    clipped = float(np.sum(over * areas))
    if clipped > 0.0:
        record_clipped_by_phase(ledger, source, phase, clipped)
    return clipped


def balance_residual(ledger: dict, live_by_source: dict[str, float]) -> dict[str, float]:
    """Per-source `emplaced - (live + sinks)`; closure means ~0 next to emplaced.

    `docked_m3` and per-phase clips sit outside the balance by design (#580).
    """
    residual = {}
    for source, accounts in ledger.items():
        if source not in SOURCES:
            raise KeyError(f"unknown buoyant-crust source in ledger: {source!r}")
        if not isinstance(accounts, dict):
            raise KeyError(f"buoyant-crust accounts for {source!r} must be a dict")
        sinks = sum(accounts.get(account, 0.0) for account in BALANCE_ACCOUNTS if account != "emplaced_m3")
        residual[source] = accounts.get("emplaced_m3", 0.0) - float(live_by_source.get(source, 0.0)) - sinks
    return residual

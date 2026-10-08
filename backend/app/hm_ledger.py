"""Debug-only source/sink accounting for mantle-lithosphere thickness (Hm).

Hm is a thermal/mechanical thickness, not conserved rock mass.  The useful identity over an
interval is therefore::

live Hm change = recorded sources - recorded sinks + type reclassification + signed residual

``phase_budget`` remains the per-phase diagnostic.  This module keeps a separate set of
physical source/sink accounts, populated from the same exact-area observations, so reports do
not add phase deltas and accounts together and double-book a change.  It is deliberately gated
by ``World.debug_diagnostics`` and reset with ``World.reset_phase_budget``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from . import elevation_lines

if TYPE_CHECKING:
    from .world import World


SCOPES = ("all", "continental_node", "oceanic_node")
_DIRECT_PHASES = {"boundary_retreat"}


def account_for_phase(phase: str) -> str:
    """Map an implementation phase to one physical Hm account.

    Accounts contain gross source and sink columns.  A phase can therefore contribute to
    both without claiming Hm itself is conserved.  Topology/remeshing is intentionally kept
    distinct from seeding and physical removal.
    """
    if phase == "oceanic_cooling_relaxation":
        return "relaxation"
    if phase == "convergent_deformation":
        return "convergent_strain"
    if phase in {"divergent_deformation", "failed_rift_thinning"}:
        return "extensional_thinning"
    if phase == "decompression_melting":
        return "melt_reset"
    if phase in {"boundary_advance", "adjacent_row_claim", "corner_notch_fill", "gap_fill"}:
        return "seeded_hm"
    if phase in {"boundary_retreat", "contested_leading_row_retreat"}:
        return "subduction_and_suture_transfer"
    if phase == "column_cap_clamp":
        return "floor_and_cap_clamps"
    if phase in {
        "plate_merge",
        "forced_plate_merge",
        "continental_relattice",
        "line_regularization",
        "plate_cleanup_removal",
        "line_growth_shrink",
    }:
        return "topology_and_regridding"
    if phase == "orogenic_relief":
        return "orogenic_relaxation"
    return "other_live_hm_writers"


def reset(world: "World") -> None:
    world.hm_source_sink_ledger = {}
    world.hm_suture_budget = {"fronts": 0, "donor_hm_m3": 0.0, "placed_hm_m3": 0.0, "unplaced_hm_m3": 0.0, "by_pair": {}}


def _new_account() -> dict:
    return {
        "calls": 0,
        "scopes": {
            scope: {"source_m3": 0.0, "sink_m3": 0.0, "reclassification_m3": 0.0}
            for scope in SCOPES
        },
    }


def _book(scope: dict, before_volume: float, after_volume: float) -> None:
    delta = after_volume - before_volume
    if delta >= 0.0:
        scope["source_m3"] += delta
    else:
        scope["sink_m3"] -= delta


def record_change(
    world: "World",
    plate,
    phase: str,
    hm_before: np.ndarray,
    codes_before: np.ndarray,
    hm_after: np.ndarray,
    codes_after: np.ndarray,
    area_before_m2: np.ndarray,
    area_after_m2: np.ndarray,
    node_ids_before: np.ndarray | None = None,
    node_ids_after: np.ndarray | None = None,
    plate_is_continental_before: bool | np.ndarray | None = None,
    plate_is_continental_after: bool | np.ndarray | None = None,
) -> None:
    """Book one measured live-Hm change into a physical account.

    Equal, aligned cells retain gross growth and thinning.  A topology or area change has no
    one-to-one column correspondence, so it is booked as the exact signed inventory change
    instead. For aligned cells, continental/oceanic scopes charge physical Hm change to the
    starting node type and book any type crossing as an explicit signed transfer.
    """
    if not getattr(world, "debug_diagnostics", False):
        return
    if phase in _DIRECT_PHASES:
        # Retreat is booked at its actual oceanic-subduction and suture-delamination sites;
        # its phase snapshot remains diagnostic only.
        return
    hm_before = np.asarray(hm_before, dtype=float)
    hm_after = np.asarray(hm_after, dtype=float)
    area_before = np.asarray(area_before_m2, dtype=float)
    area_after = np.asarray(area_after_m2, dtype=float)
    codes_before = np.asarray(codes_before)
    codes_after = np.asarray(codes_after)

    ledger = getattr(world, "hm_source_sink_ledger", None)
    if not isinstance(ledger, dict):
        ledger = world.hm_source_sink_ledger = {}
    entry = ledger.setdefault(account_for_phase(phase), _new_account())
    entry["calls"] += 1

    before_volume = hm_before * area_before
    after_volume = hm_after * area_after
    identity_aligned = (
        node_ids_before is None
        and node_ids_after is None
        or node_ids_before is not None
        and node_ids_after is not None
        and np.array_equal(node_ids_before, node_ids_after)
    )
    aligned = (
        identity_aligned
        and len(hm_before) == len(hm_after)
        and np.array_equal(area_before.shape, area_after.shape)
        and np.allclose(area_before, area_after, rtol=0.0, atol=0.0)
    )
    if aligned:
        delta = after_volume - before_volume
        entry["scopes"]["all"]["source_m3"] += float(np.maximum(delta, 0.0).sum())
        entry["scopes"]["all"]["sink_m3"] += float(np.maximum(-delta, 0.0).sum())
        default_before = plate.crust_type == "continental" if plate_is_continental_before is None else plate_is_continental_before
        default_after = plate.crust_type == "continental" if plate_is_continental_after is None else plate_is_continental_after
        cont_before = elevation_lines.effective_is_continental_from_codes(codes_before, default_before)
        cont_after = elevation_lines.effective_is_continental_from_codes(codes_after, default_after)
        # Charge Hm growth/thinning to the node's starting class. The after-volume crossing
        # the class boundary is then an explicit signed transfer between scopes.
        for scope, mask in (("continental_node", cont_before), ("oceanic_node", ~cont_before)):
            scope_delta = delta[mask]
            entry["scopes"][scope]["source_m3"] += float(np.maximum(scope_delta, 0.0).sum())
            entry["scopes"][scope]["sink_m3"] += float(np.maximum(-scope_delta, 0.0).sum())
        moved_to_cont = (~cont_before) & cont_after
        moved_to_ocean = cont_before & (~cont_after)
        transfer_to_cont = float(after_volume[moved_to_cont].sum())
        transfer_to_ocean = float(after_volume[moved_to_ocean].sum())
        entry["scopes"]["continental_node"]["reclassification_m3"] += transfer_to_cont - transfer_to_ocean
        entry["scopes"]["oceanic_node"]["reclassification_m3"] += transfer_to_ocean - transfer_to_cont
    else:
        _book(entry["scopes"]["all"], float(before_volume.sum()), float(after_volume.sum()))
        default_before = plate.crust_type == "continental" if plate_is_continental_before is None else plate_is_continental_before
        default_after = plate.crust_type == "continental" if plate_is_continental_after is None else plate_is_continental_after
        cont_before = elevation_lines.effective_is_continental_from_codes(codes_before, default_before)
        cont_after = elevation_lines.effective_is_continental_from_codes(codes_after, default_after)
        for name, before_mask, after_mask in (
            ("continental_node", cont_before, cont_after),
            ("oceanic_node", ~cont_before, ~cont_after),
        ):
            _book(
                entry["scopes"][name],
                float(before_volume[before_mask].sum()),
                float(after_volume[after_mask].sum()),
            )


def record_sink_by_mask(
    world: "World | None",
    plate,
    account: str,
    hm_m: np.ndarray,
    codes: np.ndarray,
    area_m2: np.ndarray,
    mask: np.ndarray,
) -> float:
    """Book columns about to leave the live surface to an exact-area sink account."""
    if world is None or not getattr(world, "debug_diagnostics", False) or not np.any(mask):
        return 0.0
    hm_m = np.asarray(hm_m, dtype=float)
    area_m2 = np.asarray(area_m2, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    volume = hm_m * area_m2
    cont = elevation_lines.effective_is_continental_from_codes(codes, plate.crust_type == "continental")
    entry = world.hm_source_sink_ledger.setdefault(account, _new_account())
    entry["calls"] += 1
    entry["scopes"]["all"]["sink_m3"] += float(volume[mask].sum())
    entry["scopes"]["continental_node"]["sink_m3"] += float(volume[mask & cont].sum())
    entry["scopes"]["oceanic_node"]["sink_m3"] += float(volume[mask & ~cont].sum())
    return float(volume[mask].sum())


def record_typed_sink(world: "World | None", account: str, volume_m3: float, *, continental: bool) -> None:
    """Book an already-measured terminal sink with its donor node type."""
    if world is None or not getattr(world, "debug_diagnostics", False):
        return
    volume = max(float(volume_m3), 0.0)
    entry = world.hm_source_sink_ledger.setdefault(account, _new_account())
    entry["calls"] += 1
    entry["scopes"]["all"]["sink_m3"] += volume
    scope = "continental_node" if continental else "oceanic_node"
    entry["scopes"][scope]["sink_m3"] += volume


def record_reclassification(
    world: "World | None",
    account: str,
    volume_m3: float,
    *,
    from_continental: bool,
    to_continental: bool,
) -> None:
    """Book live Hm transferred between typed-node scopes without changing total Hm."""
    if world is None or not getattr(world, "debug_diagnostics", False) or from_continental == to_continental:
        return
    volume = max(float(volume_m3), 0.0)
    entry = world.hm_source_sink_ledger.setdefault(account, _new_account())
    entry["calls"] += 1
    source_scope = "continental_node" if from_continental else "oceanic_node"
    destination_scope = "continental_node" if to_continental else "oceanic_node"
    entry["scopes"][source_scope]["reclassification_m3"] -= volume
    entry["scopes"][destination_scope]["reclassification_m3"] += volume


def record_suture_front(
    world: "World | None",
    donor_plate_id: int,
    neighbour_plate_ids: list[int],
    donor_hm_m3: float,
    placed_hm_m3: float,
    *,
    donor_is_continental: bool,
    placed_hm_continental_m3: float,
) -> None:
    """Record one connected donor front without feeding the closure accounts again.

    ``boundary_retreat`` already accounts for the resulting live change.  These counters are
    a transfer diagnostic: donor volume, survivor placement and the unplaced (delaminated)
    remainder, including bilateral consumption grouped by donor/neighbor pair.
    """
    if world is None or not getattr(world, "debug_diagnostics", False):
        return
    budget = getattr(world, "hm_suture_budget", None)
    if not isinstance(budget, dict) or "by_pair" not in budget:
        reset_suture = {"fronts": 0, "donor_hm_m3": 0.0, "placed_hm_m3": 0.0, "unplaced_hm_m3": 0.0, "by_pair": {}}
        world.hm_suture_budget = budget = reset_suture
    donor = max(float(donor_hm_m3), 0.0)
    placed = min(max(float(placed_hm_m3), 0.0), donor)
    unplaced = donor - placed
    record_typed_sink(world, "suture_delamination", unplaced, continental=donor_is_continental)
    placed_continental = min(max(float(placed_hm_continental_m3), 0.0), placed)
    placed_oceanic = placed - placed_continental
    cross_type_placed = placed_oceanic if donor_is_continental else placed_continental
    record_reclassification(
        world,
        "subduction_and_suture_transfer",
        cross_type_placed,
        from_continental=donor_is_continental,
        to_continental=not donor_is_continental,
    )
    budget["fronts"] += 1
    budget["donor_hm_m3"] += donor
    budget["placed_hm_m3"] += placed
    budget["unplaced_hm_m3"] += unplaced
    neighbours = sorted(set(int(i) for i in neighbour_plate_ids)) or [-1]
    # A front can touch several candidates.  Keep the full front under the stable donor ->
    # candidate-set key instead of fractionally inventing attribution among them.
    key = f"{int(donor_plate_id)}->{'|'.join(map(str, neighbours))}"
    row = budget["by_pair"].setdefault(
        key, {"fronts": 0, "donor_hm_m3": 0.0, "placed_hm_m3": 0.0, "unplaced_hm_m3": 0.0}
    )
    row["fronts"] += 1
    row["donor_hm_m3"] += donor
    row["placed_hm_m3"] += placed
    row["unplaced_hm_m3"] += unplaced


def cumulative_scopes(world: "World") -> dict[str, dict[str, float]]:
    """Sum every physical account by scope for interval-closure reporting."""
    out = {
        scope: {"source_m3": 0.0, "sink_m3": 0.0, "reclassification_m3": 0.0}
        for scope in SCOPES
    }
    for account in getattr(world, "hm_source_sink_ledger", {}).values():
        for scope in SCOPES:
            out[scope]["source_m3"] += account["scopes"][scope]["source_m3"]
            out[scope]["sink_m3"] += account["scopes"][scope]["sink_m3"]
            out[scope]["reclassification_m3"] += account["scopes"][scope].get("reclassification_m3", 0.0)
    return out


def inventory_scopes_m3(world: "World") -> dict[str, float]:
    """Current live Hm volume with each node's exact accounting area."""
    spacing = elevation_lines.line_spacing_rad(world.node_density)
    out = {scope: 0.0 for scope in SCOPES}
    for plate in world.plates:
        hm = np.asarray(plate.collect("mantle_lithosphere_thickness_m"), dtype=float)
        if len(hm) == 0:
            continue
        area = np.asarray(plate.accounting_areas_m2(spacing), dtype=float)
        volume = hm * area
        cont = elevation_lines.effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        out["all"] += float(volume.sum())
        out["continental_node"] += float(volume[cont].sum())
        out["oceanic_node"] += float(volume[~cont].sum())
    return out

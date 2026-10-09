"""Crust transfer at polarized collision fronts (issue #320): a consumed lower-plate column's
crust is partitioned into scraped, underthrust and lost shares and placed on the recorded
upper plate, with every share booked and every ledger closing on exact quad cell areas --
see crust_transfer.py."""

import numpy as np
import pytest

from app import (
    collision_polarity as cp,
    continental_ledger,
    cratons,
    crust_transfer,
    hm_ledger,
    lithosphere,
    mobile_cover,
    orogeny,
    quad_tectonics,
)
from app.elevation_lines import line_spacing_rad
from app.lithosphere_plate import CONTINENTAL_CONTESTED_RETREAT_MIN_RUN, SUTURE_ACCRETION_MAX_HC_M, boundary_context
from app.sparse_quad_patch import PlateWithSparseQuadPatch, cells_per_face_edge, pack_cell_keys, unpack_cell_keys
from app.world import World

DENSITY = 0.5
SPACING = line_spacing_rad(DENSITY)
N = cells_per_face_edge(SPACING)
STEP_YEARS = 100_000.0
COVER_M = 400.0
CRATON_M = 8_000.0
RESTITE_M = 2_000.0


def _block(i_range, j_range, face: int = 0) -> np.ndarray:
    jj, ii = np.meshgrid(np.arange(*j_range), np.arange(*i_range), indexing="ij")
    return pack_cell_keys(np.full(ii.size, face), ii.ravel(), jj.ravel())


def _plate(plate_id, keys, crust_type="continental", **fields) -> PlateWithSparseQuadPatch:
    count = len(keys)
    hc, hm = lithosphere.reference_thickness(crust_type)
    defaults = {"crustal_thickness_m": np.full(count, hc), "mantle_lithosphere_thickness_m": np.full(count, hm)}
    defaults.update(fields)
    plate = PlateWithSparseQuadPatch(plate_id, np.eye(3), crust_type, N, keys, fields=defaults)
    lithosphere.sync_plate_elevation(plate)
    return plate


def _lower(plate_id=1, i_range=(10, 20), j_range=(20, 32)) -> PlateWithSparseQuadPatch:
    """A continental plate whose columns carry every tracer the transfer must partition."""
    count = (i_range[1] - i_range[0]) * (j_range[1] - j_range[0])
    hc = np.full(count, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    return _plate(
        plate_id, _block(i_range, j_range), crustal_thickness_m=hc, continental_material_m=hc.copy(),
        craton_crust_m=np.full(count, CRATON_M), restite_m=np.full(count, RESTITE_M),
        mobile_cover_m=np.full(count, COVER_M), mobile_cover_continental_m=np.full(count, COVER_M),
    )


def _world(*plates) -> World:
    world = World(seed=0, plates=list(plates), next_plate_id=max(p.plate_id for p in plates) + 1, node_density=DENSITY)
    world.fault_deformation_mode = "smooth"
    continental_ledger.ensure_initialized(world)
    cratons.ensure_ledger(world)
    world.craton_ledger["initial_m3"] = cratons.live_volume_m3(world)
    mobile_cover.ensure_ledger(world)
    orogeny.ensure_budget(world)
    return world


def _columns(plate) -> tuple[np.ndarray, np.ndarray]:
    _, _, i, j = unpack_cell_keys(plate.cell_keys)
    return i, j


def _volume(plate, name, mask=None) -> float:
    values = plate.collect(name)
    areas = plate.node_areas_m2()
    if mask is None:
        return float(values @ areas)
    return float(values[mask] @ areas[mask])


def _consume(world, lower, donors, upper_ids, overriders=None, years=0.0) -> None:
    """What `_retreat` does to suture donors: book the cover of those with no upper plate
    (the transfer books the rest), accrete, remove them."""
    mobile_cover.book_removed(world, lower, donors & (upper_ids < 0), "accreted_m3")
    quad_tectonics._accrete_onto_survivors(
        lower, donors, ~donors, world, years=years, overriders=overriders, upper_plate_ids=upper_ids,
    )
    lower.remove_cells(donors)


def _assert_ledgers_close(world) -> None:
    continental_ledger.assert_closed(world)
    scale = max(world.craton_ledger["initial_m3"], 1.0)
    assert abs(cratons.balance_error_m3(world)) <= 1e-9 * scale
    assert abs(mobile_cover.balance_error_m3(world)) <= 1e-9 * max(mobile_cover.ensure_ledger(world)["initial_m3"], 1.0)


def _front(lower, upper_id, i=19):
    i_cols, _ = _columns(lower)
    donors = i_cols == i
    ids = np.where(donors, upper_id, -1)
    return donors, ids


def test_upper_plate_receives_the_configured_shares_and_the_rest_is_booked():
    lower = _lower()
    upper = _plate(2, _block((20, 40), (10, 40)))
    world = _world(lower, upper)
    donors, ids = _front(lower, 2)
    areas = lower.node_areas_m2()
    assert np.ptp(areas) > 0.0  # exact quad cells: unequal areas
    hc = _volume(lower, "crustal_thickness_m", donors)
    cover = _volume(lower, "mobile_cover_m", donors)
    material = _volume(lower, "continental_material_m", donors)
    body = hc - cover
    shares = crust_transfer.partition(world)
    upper_hc = _volume(upper, "crustal_thickness_m")
    upper_hm = upper.collect("mantle_lithosphere_thickness_m")
    upper_material = _volume(upper, "continental_material_m")
    lower_hm_survivors = lower.collect("mantle_lithosphere_thickness_m")[~donors]

    _consume(world, lower, donors, ids)

    budget = world.orogenic_relief_budget
    lost = shares.loss * body
    assert budget["suture_donated_m3"] == pytest.approx(hc)
    assert budget["suture_scraped_m3"] == pytest.approx(cover + shares.scrape * body)
    assert budget["suture_underthrust_m3"] == pytest.approx(shares.underthrust * body)
    assert budget["suture_underthrust_placed_m3"] == pytest.approx(shares.underthrust * body)
    assert budget["suture_lower_crust_subducted_m3"] == pytest.approx(lost)
    assert budget["no_outlet_subducted_m3"] == 0.0
    # The upper plate gains everything but the lost share; the lower plate keeps none of it.
    assert _volume(upper, "crustal_thickness_m") - upper_hc == pytest.approx(hc - lost, rel=1e-10)
    assert budget["upper_plate_placed_m3"] + budget["suture_underthrust_placed_m3"] == pytest.approx(hc - lost)
    assert np.allclose(lower.collect("crustal_thickness_m"), lithosphere.REFERENCE_HC_CONTINENTAL_M)
    # Continental material splits by the same source-volume fractions.
    assert _volume(upper, "continental_material_m") - upper_material == pytest.approx(material * (hc - lost) / hc, rel=1e-10)
    assert world.continental_material_ledger["collision_lower_crust_subducted_m3"] == pytest.approx(material * lost / hc)
    assert world.continental_material_ledger["collision_subducted_m3"] == 0.0
    # No Hm on either plate.
    assert np.array_equal(upper.collect("mantle_lithosphere_thickness_m"), upper_hm)
    assert np.array_equal(lower.collect("mantle_lithosphere_thickness_m"), lower_hm_survivors)
    _assert_ledgers_close(world)


def test_underthrust_stays_within_its_reach_and_buries_the_moho():
    lower = _lower()
    lower.set_fields_on_plate(mobile_cover_m=np.zeros(lower.node_count()), mobile_cover_continental_m=np.zeros(lower.node_count()))
    upper = _plate(2, _block((20, 40), (10, 40)))
    world = _world(lower, upper)
    world.suture_scrape_fraction, world.suture_underthrust_fraction, world.suture_lower_crust_loss_fraction = 0.0, 1.0, 0.0
    donors, ids = _front(lower, 2)
    hc_before = upper.collect("crustal_thickness_m")
    lag_before = upper.collect("moho_thermal_lag_c")
    elevation_before = upper.collect("elevation")

    _consume(world, lower, donors, ids)

    gained = upper.collect("crustal_thickness_m") - hc_before > 1e-6
    i_cols, _ = _columns(upper)
    hops = cratons._hops(crust_transfer.UNDERTHRUST_REACH_KM, SPACING)
    assert np.any(gained)
    assert i_cols[gained].max() <= 20 + hops
    assert world.orogenic_relief_budget["upper_plate_placed_m3"] == 0.0
    # Moho buried and elevation raised by the Hc each receiver actually gained.
    assert np.all(upper.collect("moho_thermal_lag_c")[gained] != lag_before[gained])
    hm = upper.collect("mantle_lithosphere_thickness_m")
    density = lithosphere.node_crust_density(upper.collect("crust_type_code"), upper.crust_type)
    expected = elevation_before + (
        lithosphere.isostatic_elevation(upper.collect("crustal_thickness_m"), hm, density)
        - lithosphere.isostatic_elevation(hc_before, hm, density)
    )
    assert np.allclose(upper.collect("elevation"), expected)
    _assert_ledgers_close(world)


def test_underthrust_the_reach_cannot_hold_is_thrust_up_with_the_scraped_share():
    lower = _lower()
    # The upper plate's front is at the cap; only its interior has room.
    upper_keys = _block((20, 40), (10, 40))
    upper_i = unpack_cell_keys(upper_keys)[2]
    hc = np.where(upper_i < 24, SUTURE_ACCRETION_MAX_HC_M, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    upper = _plate(2, upper_keys, crustal_thickness_m=hc)
    world = _world(lower, upper)
    donors, ids = _front(lower, 2)
    donated = _volume(lower, "crustal_thickness_m", donors)
    before = _volume(upper, "crustal_thickness_m")

    _consume(world, lower, donors, ids)

    budget = world.orogenic_relief_budget
    assert budget["suture_underthrust_m3"] > 0.0
    assert budget["suture_underthrust_placed_m3"] == 0.0
    lost = budget["suture_lower_crust_subducted_m3"]
    assert budget["upper_plate_placed_m3"] == pytest.approx(donated - lost)
    assert _volume(upper, "crustal_thickness_m") - before == pytest.approx(donated - lost, rel=1e-10)
    _assert_ledgers_close(world)


def test_a_full_upper_plate_hands_on_then_subducts_what_has_no_outlet():
    lower = _lower()
    full = SUTURE_ACCRETION_MAX_HC_M
    upper = _plate(2, _block((20, 24), (20, 32)), crustal_thickness_m=np.full(48, full))
    neighbour = _plate(3, _block((24, 26), (20, 32)), crustal_thickness_m=np.full(24, full - 1_000.0))
    world = _world(lower, upper, neighbour)
    donors, ids = _front(lower, 2)
    lower_areas = lower.node_areas_m2()[donors]
    donated = _volume(lower, "crustal_thickness_m", donors)
    room = float(np.full(24, 1_000.0) @ neighbour.node_areas_m2())
    material = _volume(lower, "continental_material_m", donors)

    _consume(world, lower, donors, ids, overriders=[upper, neighbour])

    budget = world.orogenic_relief_budget
    lost = budget["suture_lower_crust_subducted_m3"]
    assert budget["upper_plate_placed_m3"] == 0.0
    assert budget["overrider_placed_m3"] == pytest.approx(room, rel=1e-9)
    unplaced = donated - lost - room
    assert budget["no_outlet_subducted_m3"] == pytest.approx(unplaced, rel=1e-9)
    assert world.continental_material_ledger["collision_subducted_m3"] == pytest.approx(material * unplaced / donated)
    # The cover rides on the thrust-up crust, so it shares the no-outlet share's fate.
    cover = float(np.full(int(donors.sum()), COVER_M) @ lower_areas)
    covers = world.mobile_cover_ledger
    assert covers["subducted_m3"] == pytest.approx(cover * unplaced / (donated - lost), rel=1e-9)
    assert covers["accreted_m3"] == pytest.approx(cover * (1.0 - unplaced / (donated - lost)), rel=1e-9)
    _assert_ledgers_close(world)


def test_cratonic_donor_is_reworked_where_it_lands_and_subducted_with_the_lost_share():
    lower = _lower()
    upper = _plate(2, _block((20, 40), (10, 40)))
    world = _world(lower, upper)
    donors, ids = _front(lower, 2)
    craton = _volume(lower, "craton_crust_m", donors)
    restite = _volume(lower, "restite_m", donors)
    hc = _volume(lower, "crustal_thickness_m", donors)

    _consume(world, lower, donors, ids)

    lost_share = world.orogenic_relief_budget["suture_lower_crust_subducted_m3"] / hc
    assert world.craton_ledger["subducted_m3"] == pytest.approx(lost_share * craton)
    assert world.craton_ledger["collision_reworked_m3"] == pytest.approx((1.0 - lost_share) * craton)
    assert _volume(upper, "craton_crust_m") == 0.0
    assert world.orogenic_relief_budget["restite_subducted_m3"] == pytest.approx(lost_share * restite)
    assert _volume(upper, "restite_m") == pytest.approx((1.0 - lost_share) * restite, rel=1e-10)
    _assert_ledgers_close(world)


def test_separate_fronts_go_to_their_own_upper_plates():
    lower = _lower(i_range=(10, 20), j_range=(20, 32))
    east = _plate(2, _block((20, 30), (20, 32)))
    west = _plate(3, _block((0, 10), (20, 32)))
    world = _world(lower, east, west)
    i_cols, _ = _columns(lower)
    donors = (i_cols == 19) | (i_cols == 10)
    ids = np.where(i_cols == 19, 2, np.where(i_cols == 10, 3, -1))
    east_donated = _volume(lower, "crustal_thickness_m", i_cols == 19)
    west_donated = _volume(lower, "crustal_thickness_m", i_cols == 10)
    before = {p.plate_id: _volume(p, "crustal_thickness_m") for p in (east, west)}
    world.debug_diagnostics = True

    _consume(world, lower, donors, ids)

    keep = 1.0 - world.orogenic_relief_budget["suture_lower_crust_subducted_m3"] / (east_donated + west_donated)
    assert _volume(east, "crustal_thickness_m") - before[2] == pytest.approx(keep * east_donated, rel=1e-10)
    assert _volume(west, "crustal_thickness_m") - before[3] == pytest.approx(keep * west_donated, rel=1e-10)
    assert world.suture_transfer_stats["fronts"] == 2
    assert world.hm_suture_budget["placed_hm_m3"] == 0.0
    _assert_ledgers_close(world)


def test_a_vanished_upper_plate_falls_back_to_same_plate_accretion():
    lower = _lower()
    world = _world(lower)
    donors, ids = _front(lower, 7)
    donated = _volume(lower, "crustal_thickness_m", donors)

    _consume(world, lower, donors, ids)

    assert world.orogenic_relief_budget["suture_scraped_m3"] == 0.0
    assert _volume(lower, "crustal_thickness_m") - lithosphere.REFERENCE_HC_CONTINENTAL_M * float(
        lower.node_areas_m2().sum()
    ) == pytest.approx(donated, rel=1e-9)
    _assert_ledgers_close(world)


def test_polarized_retreat_moves_lower_plate_crust_to_the_upper_plate():
    """End to end through the prepass and `_retreat`: the lower plate retreats, the upper
    plate keeps its cells and takes the crust, and the Hm sink matches the donors' Hm."""
    from unit_tests.test_collision_polarity import _column, _front_with_lower_a, _move

    world, a, b, _ = _front_with_lower_a()
    count = a.node_count()
    hc = a.collect("crustal_thickness_m")
    a.set_fields_on_plate(
        continental_material_m=hc.copy(), restite_m=np.full(count, RESTITE_M), mobile_cover_m=np.full(count, COVER_M),
    )
    continental_ledger.ensure_initialized(world)
    cratons.ensure_ledger(world)
    world.craton_ledger["initial_m3"] = cratons.live_volume_m3(world)
    mobile_cover.ensure_ledger(world)
    _move(b, _column(20), _column(19))
    cp.observe_contacts(world, STEP_YEARS)
    world.debug_diagnostics = True
    ctx = boundary_context(
        world, a, [b], STEP_YEARS,
        lambda mask: quad_tectonics.components_of_at_least(a, mask, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN),
    )
    assert np.any(ctx.suture_hm_subduct)
    assert set(ctx.suture_upper_plate_id[ctx.suture_hm_subduct].tolist()) == {b.plate_id}
    assert np.all(ctx.suture_upper_plate_id[~ctx.suture_hm_subduct] == -1)
    hm_a = a.collect("mantle_lithosphere_thickness_m")
    areas = a.node_areas_m2()
    b_hc, b_count = _volume(b, "crustal_thickness_m"), b.node_count()

    survivors = quad_tectonics._retreat(a, world, ctx, SPACING * 2, a.node_count(), STEP_YEARS)

    removed = ~survivors
    assert np.any(removed) and np.all(ctx.suture_hm_subduct[removed])
    donated = float(hc[removed] @ areas[removed])
    lost = world.orogenic_relief_budget["suture_lower_crust_subducted_m3"]
    assert lost > 0.0
    assert b.node_count() == b_count
    assert _volume(b, "crustal_thickness_m") - b_hc == pytest.approx(donated - lost, rel=1e-10)
    sink = world.hm_source_sink_ledger["suture_hm_subducted_m3"]["scopes"]["all"]["sink_m3"]
    assert sink == pytest.approx(float(hm_a[removed] @ areas[removed]))
    assert world.hm_suture_budget["placed_hm_m3"] == 0.0
    _assert_ledgers_close(world)


def test_transfer_column_is_callable_without_the_source_plate():
    """The docking interface (issue #321): a column built from volumes alone."""
    upper = _plate(2, _block((20, 40), (10, 40)))
    world = _world(upper)
    world.debug_diagnostics = True
    front = upper.all_points_and_elevation()[0][_columns(upper)[0] == 20]
    column = crust_transfer.ConsumedColumn(
        plate_id=99, front_world=front * 1.0, approach_world=None, hc_m3=1e15, hm_m3=2e15,
        mobile_cover_m3=0.0, material_m3=1e15, craton_m3=0.0, restite_m3=0.0,
    )
    before = _volume(upper, "crustal_thickness_m")
    shares = crust_transfer.CrustPartition(scrape=0.5, underthrust=0.5, loss=0.0)

    result = crust_transfer.transfer_column(world, column, upper, shares=shares)

    assert result.scraped_m3 == pytest.approx(5e14) and result.underthrust_m3 == pytest.approx(5e14)
    assert _volume(upper, "crustal_thickness_m") - before == pytest.approx(1e15, rel=1e-10)
    assert world.hm_source_sink_ledger["suture_hm_subducted_m3"]["scopes"]["all"]["sink_m3"] == pytest.approx(2e15)
    stats = world.suture_transfer_stats
    assert stats["fronts"] == 1 and stats["steps_with_fronts"] == 1
    assert stats["max_step_seconds"] == pytest.approx(stats["seconds"])
    assert "by_step" not in stats


@pytest.mark.parametrize("shares", [(0.5, 0.5, 0.5), (1.2, -0.2, 0.0), (float("nan"), 0.5, 0.5)])
def test_partition_must_be_non_negative_and_sum_to_one(shares):
    world = World(seed=0)
    world.suture_scrape_fraction, world.suture_underthrust_fraction, world.suture_lower_crust_loss_fraction = shares
    with pytest.raises(ValueError):
        crust_transfer.partition(world)


def test_partition_is_reported_in_the_summary():
    world = World(seed=0)
    report = crust_transfer.summary(world)
    assert report["partition"]["scrape_fraction"] == world.suture_scrape_fraction
    assert report["partition"]["underthrust_reach_km"] == crust_transfer.UNDERTHRUST_REACH_KM
    assert sum(report["partition"][k] for k in ("scrape_fraction", "underthrust_fraction", "lower_crust_loss_fraction")) == pytest.approx(1.0)


def test_hm_ledger_closes_on_a_transfer():
    lower = _lower()
    upper = _plate(2, _block((20, 40), (10, 40)))
    world = _world(lower, upper)
    world.debug_diagnostics = True
    donors, ids = _front(lower, 2)
    inventory = hm_ledger.inventory_scopes_m3(world)["all"]
    donor_hm = _volume(lower, "mantle_lithosphere_thickness_m", donors)

    _consume(world, lower, donors, ids)

    assert inventory - hm_ledger.inventory_scopes_m3(world)["all"] == pytest.approx(donor_hm, rel=1e-10)
    sink = world.hm_source_sink_ledger["suture_hm_subducted_m3"]["scopes"]["all"]["sink_m3"]
    assert sink == pytest.approx(donor_hm)


def test_fronts_on_one_upper_plate_share_its_root_shedding_allowance():
    # Two separate fronts against one hot, nearly saturated upper plate in one step.
    lower = _lower(j_range=(10, 40))
    upper_hc = np.full(20 * 40, SUTURE_ACCRETION_MAX_HC_M - 300.0)
    # Thin mantle lid: a hot Moho, so its roots are eligible to founder.
    upper = _plate(
        2, _block((20, 40), (5, 45)), crustal_thickness_m=upper_hc,
        mantle_lithosphere_thickness_m=np.full(upper_hc.size, 20_000.0),
    )
    world = _world(lower, upper)
    world.suture_root_capacity = {}
    i_cols, j_cols = _columns(lower)
    donors = (i_cols == 19) & ((j_cols < 18) | (j_cols >= 32))
    ids = np.where(donors, 2, -1)
    continental = np.ones(upper.node_count(), dtype=bool)
    allowance = float(orogeny.plate_delamination_capacity_m3(upper, continental, 1_000_000.0).sum())
    assert allowance > 0.0

    _consume(world, lower, donors, ids, years=1_000_000.0)

    shed = world.orogenic_relief_budget["delamination_completed_m3"]
    assert shed > 0.0
    assert shed <= allowance * (1.0 + 1e-9)
    _, left = world.suture_root_capacity[2]
    assert float(left.sum()) == pytest.approx(allowance - shed, rel=1e-9)
    _assert_ledgers_close(world)


def test_overflow_is_never_handed_to_an_oceanic_neighbour():
    lower = _lower()
    upper = _plate(2, _block((20, 24), (20, 32)), crustal_thickness_m=np.full(48, SUTURE_ACCRETION_MAX_HC_M))
    ocean = _plate(3, _block((20, 30), (32, 40)), "oceanic")
    world = _world(lower, upper, ocean)
    donors, ids = _front(lower, 2)
    ocean_hc = ocean.collect("crustal_thickness_m")

    _consume(world, lower, donors, ids, overriders=[ocean, upper])

    assert np.array_equal(ocean.collect("crustal_thickness_m"), ocean_hc)
    assert world.orogenic_relief_budget["overrider_placed_m3"] == 0.0
    assert world.orogenic_relief_budget["no_outlet_subducted_m3"] > 0.0
    _assert_ledgers_close(world)


def test_an_upper_plate_with_no_continental_columns_takes_no_underthrust():
    upper = _plate(2, _block((20, 40), (10, 40)), "oceanic")
    world = _world(upper)
    front = upper.all_points_and_elevation()[0][_columns(upper)[0] == 20]
    column = crust_transfer.ConsumedColumn(
        plate_id=99, front_world=front * 1.0, approach_world=None, hc_m3=1e14, hm_m3=0.0,
        mobile_cover_m3=0.0, material_m3=1e14, craton_m3=0.0, restite_m3=0.0,
    )

    result = crust_transfer.transfer_column(world, column, upper)

    assert result.underthrust_m3 > 0.0
    assert result.underthrust_placed_m3 == 0.0
    assert result.upper_placed_m3 + result.unplaced_m3 == pytest.approx(result.scraped_m3 + result.underthrust_m3)


def test_the_approach_comes_from_the_nodes_facing_the_upper_plate(monkeypatch):
    lower = _lower()
    upper = _plate(2, _block((20, 40), (10, 40)))
    world = _world(lower, upper)
    donors, ids = _front(lower, 2)
    toward_upper, toward_third = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    convergence = np.tile(toward_third, (lower.node_count(), 1))
    neighbours = np.full(lower.node_count(), 3)
    facing = np.flatnonzero(donors)[::2]
    convergence[facing] = toward_upper
    neighbours[facing] = 2
    seen = []
    original = quad_tectonics._place_on_overrider

    def spy(world, over, front_world, volume, material, approach, *args, **kwargs):
        seen.append(approach)
        return original(world, over, front_world, volume, material, approach, *args, **kwargs)

    monkeypatch.setattr(quad_tectonics, "_place_on_overrider", spy)
    crust_transfer.transfer_fronts(world, lower, donors, ids, convergence, neighbours, 0.0, None)

    assert len(seen) == 1
    assert np.allclose(seen[0], -toward_upper)


def test_invalid_partition_fails_the_step_before_any_plate_deforms():
    from app import world as world_mod

    lower = _lower()
    upper = _plate(2, _block((20, 40), (10, 40)))
    world = _world(lower, upper)
    world.suture_scrape_fraction = 0.8
    deformed = []
    monkeypatch_target = quad_tectonics.deform

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(quad_tectonics, "deform", lambda *a, **k: deformed.append(a) or monkeypatch_target(*a, **k))
        with pytest.raises(ValueError):
            world_mod.step_world(world, STEP_YEARS)
    assert deformed == []
    assert "error" in crust_transfer.summary(world)["partition"]


def test_a_polarized_donor_nearest_an_oceanic_plate_still_transfers(monkeypatch):
    """At a triple junction a lower-plate front cell's nearest neighbour can be an oceanic
    third plate. It must still be a suture donor: not subducted, not charged to the #177
    oceanic-override budget, and its crust still goes to the frozen upper plate."""
    from app import torque
    from unit_tests.test_collision_polarity import _column, _front_with_lower_a, _move

    world, a, b, _ = _front_with_lower_a()
    _move(b, _column(20), _column(19))
    cp.observe_contacts(world, STEP_YEARS)
    original = torque.gather_boundary_force_inputs

    def oceanic_nearest(*args, **kwargs):
        inputs = original(*args, **kwargs)
        inputs.neighbor_is_oceanic[:] = True
        return inputs

    monkeypatch.setattr(torque, "gather_boundary_force_inputs", oceanic_nearest)
    ctx = boundary_context(
        world, a, [b], STEP_YEARS,
        lambda mask: quad_tectonics.components_of_at_least(a, mask, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN),
    )
    assert np.any(ctx.suture_hm_subduct)
    assert np.all(ctx.accrete[ctx.suture_hm_subduct])
    hc = a.collect("crustal_thickness_m")
    areas = a.node_areas_m2()
    b_hc = _volume(b, "crustal_thickness_m")
    deep_before = continental_ledger.inventories(world)["deeply_subducted_m3"]

    survivors = quad_tectonics._retreat(a, world, ctx, SPACING * 2, a.node_count(), STEP_YEARS)

    removed = ~survivors
    assert np.any(removed)
    donated = float(hc[removed] @ areas[removed])
    lost = world.orogenic_relief_budget["suture_lower_crust_subducted_m3"]
    assert _volume(b, "crustal_thickness_m") - b_hc == pytest.approx(donated - lost, rel=1e-10)
    assert continental_ledger.inventories(world)["deeply_subducted_m3"] == deep_before

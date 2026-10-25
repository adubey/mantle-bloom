import math

import numpy as np
import pytest

from app import persistence, terranes
from app.world import World, generate_world
from .quad_fixtures import quad_plate


def test_docking_probability_matches_rfc_examples():
    """Ontong-Java-like resists, old Caribbean-like subducts (#580)."""
    assert terranes.docking_probability(23_500.0, 100.0) == pytest.approx(0.95, abs=0.005)
    assert terranes.docking_probability(5_000.0, 100.0) == pytest.approx(0.04, abs=0.005)
    assert terranes.docking_probability(23_500.0, 100.0, width_m=100.0) > 0.99
    with pytest.raises(ValueError):
        terranes.docking_probability(math.nan, 100.0)
    with pytest.raises(ValueError):
        terranes.docking_probability(5_000.0, math.inf)


def test_docking_probability_stays_a_probability_and_orders_hosts():
    """In [0, 1], rising with excess and falling with age (#580)."""
    assert 0.0 <= terranes.docking_probability(-5_000.0, 1e9) <= 1.0
    assert terranes.docking_probability(20_000.0, 2.0) > terranes.docking_probability(5_000.0, 2.0)
    assert terranes.docking_probability(8_000.0, 2.0) > terranes.docking_probability(8_000.0, 100.0)


def test_threshold_excess_limits():
    """New crust faces the young threshold, ancient crust the old one (#580)."""
    assert terranes.threshold_excess_m(0.0) == pytest.approx(terranes.DOCK_E_YOUNG_M)
    assert terranes.threshold_excess_m(-5.0) == pytest.approx(terranes.DOCK_E_YOUNG_M)
    assert terranes.threshold_excess_m(1e9) == pytest.approx(terranes.DOCK_E_OLD_M)
    assert terranes.threshold_excess_m(100.0, e_young_m=1e4, e_old_m=2e4) == pytest.approx(
        1e4 + 1e4 * (1.0 - math.exp(-5.0))
    )
    with pytest.raises(ValueError):
        terranes.threshold_excess_m(math.nan)


def test_docking_draw_stable_per_encounter():
    """Same encounter always draws the same outcome; edges are exact (#580)."""
    args = (11, 3, 42, 12e6)
    first = terranes.docking_draw(*args, 0.5)
    assert isinstance(first, bool)
    assert terranes.docking_draw(*args, 0.5) == first
    assert terranes.docking_draw(11, 3, 42, 12e6, 1.0) is True
    assert terranes.docking_draw(11, 3, 42, 12e6, 0.0) is False
    with pytest.raises(ValueError):
        terranes.docking_draw(11, 3, 42, 12e6, math.nan)
    with pytest.raises(ValueError):
        terranes.docking_draw(11, 3, 42, 12e6, 1.5)


def test_group_patches_finds_components():
    """Union-find over explicit pairs, members only, sorted by minimum (#580)."""
    pairs = np.array([[0, 1], [1, 2], [0, 2], [4, 5], [6, 7]])
    mask = np.array([True, True, True, False, True, True, False, True])
    groups = terranes.group_patches(8, pairs, mask)
    assert [list(g) for g in groups] == [[0, 1, 2], [4, 5], [7]]
    assert terranes.group_patches(4, np.zeros((0, 2), dtype=int), np.zeros(4, dtype=bool)) == []
    with pytest.raises(ValueError):
        terranes.group_patches(4, np.zeros((0, 2), dtype=int), np.zeros(3, dtype=bool))
    with pytest.raises(ValueError):
        terranes.group_patches(4, np.array([[0.5, 1.5]]), np.ones(4, dtype=bool))


def test_patch_excess_and_age():
    """Area-weighted excess and volume-weighted age (#580)."""
    hc = np.array([26_000.0, 26_000.0, 11_000.0])
    areas = np.array([1e12, 3e12, 1e12])
    mask = np.array([True, True, False])
    assert terranes.patch_excess_m(hc, areas, mask, 0.0) == pytest.approx(15_000.0)
    plateau = np.array([15_000.0, 15_000.0, 0.0])
    arc = np.zeros(3)
    moment = np.array([15_000.0 * 10e6, 15_000.0 * 10e6, 0.0])
    assert terranes.patch_age_years(plateau, arc, moment, mask, 50e6) == pytest.approx(40e6)
    assert terranes.patch_age_years(np.zeros(3), np.zeros(3), np.zeros(3), mask, 50e6) == 0.0


def test_docking_loss_keeps_exact_lower_crust_share():
    """The trench takes exactly the lower-crust fraction (#580)."""
    plateau = np.array([15_000.0, 8_000.0])
    arc = np.array([0.0, 2_000.0])
    moment = np.array([15_000.0 * 10e6, 10_000.0 * 30e6])
    (kept_p, kept_a, kept_m), (lost_p, lost_a, lost_m) = terranes.split_docking_loss(
        plateau, arc, moment
    )
    assert kept_p == pytest.approx(plateau * 0.8)
    assert kept_a == pytest.approx(arc * 0.8)
    assert lost_a == pytest.approx(arc * 0.2)
    assert kept_p + lost_p == pytest.approx(plateau)
    assert kept_m + lost_m == pytest.approx(moment)


def test_add_tracers_carries_shares_with_hc():
    """Transfers and accretion sum the riding shares elementwise (#580)."""
    carried = terranes.add_tracers(
        np.array([1_000.0, 0.0]),
        np.array([0.0, 2_000.0]),
        np.array([1_000.0 * 5e6, 2_000.0 * 7e6]),
        np.array([500.0, 500.0]),
        np.array([500.0, 500.0]),
        np.array([500.0 * 9e6, 500.0 * 9e6]),
    )
    assert carried[0].tolist() == pytest.approx([1_500.0, 500.0])
    assert carried[1].tolist() == pytest.approx([500.0, 2_500.0])
    assert carried[2].tolist() == pytest.approx([1_000.0 * 5e6 + 500.0 * 9e6, 2_000.0 * 7e6 + 500.0 * 9e6])


def test_conversion_splits_preserved_and_foundered():
    """Preserved juvenile plus foundered mafic equals the docked crust (#580)."""
    plateau = np.array([12_000.0, 0.0])
    arc = np.array([0.0, 9_000.0])
    juvenile, foundered = terranes.conversion_split(plateau, arc)
    assert juvenile + foundered == pytest.approx(plateau + arc)
    assert juvenile == pytest.approx((plateau + arc) * terranes.PLATEAU_PRESERVED_FRACTION)
    custom, _ = terranes.conversion_split(plateau, arc, 0.5)
    assert custom == pytest.approx((plateau + arc) * 0.5)


def test_lifecycle_matrix_closes_ledger_every_stage():
    """Emplace, thin, reset, dock twice, convert, audit: balance holds (#580)."""
    areas = np.full(6, 1e12)
    ledger = terranes.new_ledger()

    plateau = np.full(6, 15_000.0)
    arc = np.zeros(6)
    moment = np.full(6, 15_000.0 * 10e6)
    hc = np.full(6, 26_000.0)
    oceanic = np.zeros(6, dtype=bool)
    live = {"plume_head": float(np.sum(plateau) * 1e12)}
    terranes.record(ledger, "plume_head", "emplaced_m3", live["plume_head"])
    assert terranes.balance_residual(ledger, live)["plume_head"] == pytest.approx(0.0, abs=1e6)

    plateau, arc, moment = terranes.scale_tracers(plateau, arc, moment, 0.8)
    thinned = live["plume_head"] * 0.2
    terranes.record(ledger, "plume_head", "thinned_m3", thinned)
    live["plume_head"] *= 0.8
    assert terranes.balance_residual(ledger, live)["plume_head"] == pytest.approx(0.0, abs=1e6)

    reset_cells = np.array([True, True, False, False, False, False])
    reset_vol = float(np.sum(plateau[reset_cells]) * 1e12)
    cleared = terranes.clear_tracers(int(np.sum(reset_cells)))
    plateau[reset_cells], arc[reset_cells], moment[reset_cells] = (column for column in cleared)
    terranes.record(ledger, "plume_head", "reset_m3", reset_vol)
    live["plume_head"] = float(np.sum(plateau) * 1e12)
    assert terranes.balance_residual(ledger, live)["plume_head"] == pytest.approx(0.0, abs=1e6)

    loss_cells = np.array([False, False, True, True, False, False])
    (kept_p, _, kept_m), (lost_p, _, _) = terranes.split_docking_loss(
        plateau[loss_cells], arc[loss_cells], moment[loss_cells]
    )
    lost_vol = float(np.sum(lost_p) * 1e12)
    terranes.record(ledger, "plume_head", "subducted_m3", lost_vol)
    plateau[loss_cells], moment[loss_cells] = kept_p, kept_m
    live["plume_head"] = float(np.sum(plateau) * 1e12)
    assert terranes.balance_residual(ledger, live)["plume_head"] == pytest.approx(0.0, abs=1e6)

    terranes.record_docked(ledger, "plume_head", 1e15)
    terranes.record_docked(ledger, "plume_head", 1e15)
    assert ledger["plume_head"][terranes.DOCKED_ACCOUNT] == pytest.approx(2e15)
    assert terranes.balance_residual(ledger, live)["plume_head"] == pytest.approx(0.0, abs=1e6)

    weathered_hc = hc.copy()
    weathered_hc[2] = 12_000.0
    clipped = terranes.apply_phase_audit(
        plateau, arc, moment, weathered_hc, areas, oceanic, ledger, "plume_head", "eroded"
    )
    assert clipped == pytest.approx(8.6e15)
    assert moment[2] == pytest.approx(9.6e10 / 9.6)
    assert ledger["plume_head"][terranes.CLIPPED_BY_PHASE]["eroded"] == pytest.approx(clipped)
    live["plume_head"] = float(np.sum(plateau) * 1e12)
    assert terranes.balance_residual(ledger, live)["plume_head"] == pytest.approx(0.0, abs=1e6)

    juvenile, foundered = terranes.conversion_split(plateau, arc)
    terranes.record(ledger, "plume_head", "converted_m3", float(np.sum(juvenile) * 1e12))
    terranes.record(ledger, "plume_head", "foundered_m3", float(np.sum(foundered) * 1e12))
    plateau[:], arc[:], moment[:] = 0.0, 0.0, 0.0
    live["plume_head"] = 0.0
    terranes.record(ledger, "plume_head", "unattributed_m3", 0.0)
    assert terranes.balance_residual(ledger, live)["plume_head"] == pytest.approx(0.0, abs=1e6)
    assert terranes.balance_residual(ledger, {"arc": 0.0, "plume_tail": 0.0})["arc"] == 0.0


def test_phase_audit_clips_continental_and_oceanic_overfill():
    """Tracers vanish on continental hosts and cap at excess on oceanic ones (#580)."""
    ledger = terranes.new_ledger()
    plateau = np.array([15_000.0, 15_000.0, 4_000.0])
    arc = np.zeros(3)
    moment = np.array([15_000.0 * 10e6, 15_000.0 * 10e6, 4_000.0 * 10e6])
    hc = np.array([26_000.0, 12_000.0, 26_000.0])
    areas = np.full(3, 1e12)
    host = np.array([False, False, True])
    terranes.record(ledger, "plume_head", "emplaced_m3", float(np.sum(plateau) * 1e12))
    clipped = terranes.apply_phase_audit(
        plateau, arc, moment, hc, areas, host, ledger, "plume_head", "rift"
    )
    assert clipped == pytest.approx((14_000.0 + 4_000.0) * 1e12)
    assert plateau.tolist() == pytest.approx([15_000.0, 1_000.0, 0.0])
    assert moment.tolist() == pytest.approx([15_000.0 * 10e6, 15_000.0 * 10e6 / 15.0, 0.0])
    assert ledger["plume_head"]["clipped_m3"] == pytest.approx(clipped)
    assert ledger["plume_head"][terranes.CLIPPED_BY_PHASE]["rift"] == pytest.approx(clipped)
    live = {"plume_head": float(np.sum(plateau) * 1e12)}
    assert terranes.balance_residual(ledger, live)["plume_head"] == pytest.approx(0.0, abs=1e6)
    with pytest.raises(TypeError):
        terranes.apply_phase_audit(
            [1_000.0], [0.0], [0.0], hc[:1], areas[:1], host[:1], ledger, "plume_head", "rift"
        )
    with pytest.raises(ValueError):
        terranes.apply_phase_audit(
            plateau, arc, moment, np.array([math.nan, 12_000.0, 26_000.0]), areas, host,
            ledger, "plume_head", "rift",
        )


def test_ledger_rejects_unknown_names_and_bad_volumes():
    """Unknown sources/accounts, negative and NaN volumes fail fast (#580)."""
    ledger = terranes.new_ledger()
    with pytest.raises(KeyError):
        terranes.record(ledger, "ocean", "emplaced_m3", 1.0)
    with pytest.raises(KeyError):
        terranes.record(ledger, "arc", "sunk_m3", 1.0)
    with pytest.raises(ValueError):
        terranes.record(ledger, "arc", "emplaced_m3", -1.0)
    with pytest.raises(ValueError):
        terranes.record(ledger, "arc", "emplaced_m3", math.nan)
    with pytest.raises(KeyError):
        terranes.record_docked(ledger, "ocean", 1.0)
    with pytest.raises(ValueError):
        terranes.record(ledger, "arc", terranes.CLIPPED_BY_PHASE, 1.0)
    with pytest.raises(ValueError):
        terranes.record(ledger, "arc", terranes.DOCKED_ACCOUNT, 1.0)
    with pytest.raises(ValueError):
        terranes.record_docked(ledger, "arc", -1.0)
    with pytest.raises(ValueError):
        terranes.record_docked(ledger, "arc", math.nan)
    with pytest.raises(KeyError):
        terranes.record_docked(ledger, "nope", 1.0)
    with pytest.raises(KeyError):
        terranes.balance_residual({"ocean": {}}, {})
    terranes.record(ledger, "arc", "unattributed_m3", 0.0)
    assert ledger["arc"]["unattributed_m3"] == 0.0


def test_world_ledger_field_and_backfill():
    """New saves carry a full ledger; old and partial saves backfill it (#580)."""
    world = World(seed=11, plates=[])
    assert world.buoyant_crust_ledger == terranes.new_ledger()
    terranes.record(world.buoyant_crust_ledger, "arc", "emplaced_m3", 1e12)
    del world.__dict__["buoyant_crust_ledger"]
    persistence._backfill_added_fields(world)
    assert world.buoyant_crust_ledger == terranes.new_ledger()
    world.buoyant_crust_ledger = {"plume_head": {"emplaced_m3": 5.0}}
    persistence._backfill_added_fields(world)
    repaired = world.buoyant_crust_ledger["plume_head"]
    assert repaired["emplaced_m3"] == 5.0
    assert repaired["docked_m3"] == 0.0
    assert repaired["clipped_by_phase"] == {}
    assert repaired["foundered_m3"] == 0.0


def test_round_trip_preserves_buoyant_ledger():
    """Booked ledger entries survive a save/load round trip (#580)."""
    world = generate_world(seed=11, node_density=0.5, num_plates=6)
    terranes.record(world.buoyant_crust_ledger, "plume_head", "emplaced_m3", 1e12)
    terranes.record_docked(world.buoyant_crust_ledger, "plume_head", 2e12)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.buoyant_crust_ledger == world.buoyant_crust_ledger


def test_record_clipped_by_phase_validates():
    """Per-phase clips need finite non-negative volumes (#580)."""
    ledger = terranes.new_ledger()
    with pytest.raises(ValueError):
        terranes.record_clipped_by_phase(ledger, "arc", "eroded", -1.0)
    with pytest.raises(ValueError):
        terranes.record_clipped_by_phase(ledger, "arc", "eroded", math.nan)
    terranes.record_clipped_by_phase(ledger, "arc", "eroded", 3.0)
    assert ledger["arc"]["clipped_m3"] == 3.0
    assert ledger["arc"][terranes.CLIPPED_BY_PHASE] == {"eroded": 3.0}


def test_old_plates_default_tracers_to_zero():
    """Plates made before the tracers read zero plateau and age moment (#580)."""
    plate = quad_plate(0, "oceanic")
    assert np.all(plate.collect("plateau_crust_m") == 0.0)
    assert np.all(plate.collect("buoyant_crust_time_m_yr") == 0.0)

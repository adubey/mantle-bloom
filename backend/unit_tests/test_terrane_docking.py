"""Terrane docking (issue #321): a continental terrane consumed with its oceanic carrier against
a continental overrider docks onto that plate through the crust-transfer partition, its Hm going
down with the carrier's slab. Without a continental overrider, or for what the overrider can't
hold, it stays on its carrier -- see crust_transfer.py's "Terranes"."""

import numpy as np
import pytest

from app import collision_polarity as cp, crust_transfer, lithosphere, quad_tectonics
from app.elevation_lines import CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC, line_spacing_rad
from app.lithosphere_plate import CONTINENTAL_CONTESTED_RETREAT_MIN_RUN, SUTURE_ACCRETION_MAX_HC_M, boundary_context
from unit_tests.test_crust_transfer import (
    COVER_M,
    CRATON_M,
    RESTITE_M,
    STEP_YEARS,
    _assert_ledgers_close,
    _block,
    _columns,
    _plate,
    _volume,
    _world,
)

TERRANE_HC_M = lithosphere.REFERENCE_HC_CONTINENTAL_M
TERRANE_HM_M = 90_000.0
J = (20, 30)


def _carrier(i_range=(10, 24), terrane_from=18, j_range=J, plate_id=1):
    """An oceanic plate whose cells from column `terrane_from` on are a continental terrane
    carrying every tracer a dock must move."""
    plate = _plate(plate_id, _block(i_range, j_range), "oceanic")
    i, _ = _columns(plate)
    terrane = i >= terrane_from
    count = plate.node_count()

    def on_terrane(value, base=0.0):
        return np.where(terrane, value, base)

    codes = plate.collect("crust_type_code")
    codes[terrane] = CRUST_TYPE_CONTINENTAL
    plate.set_fields_on_plate(
        crust_type_code=codes,
        crustal_thickness_m=on_terrane(TERRANE_HC_M, plate.collect("crustal_thickness_m")),
        mantle_lithosphere_thickness_m=on_terrane(TERRANE_HM_M, plate.collect("mantle_lithosphere_thickness_m")),
        continental_material_m=on_terrane(TERRANE_HC_M),
        craton_crust_m=on_terrane(CRATON_M),
        craton_formed_years=on_terrane(-2.0e9, plate.collect("craton_formed_years")),
        restite_m=on_terrane(RESTITE_M),
        # Overlapped long enough for the craton delay to release it (cratons.retreat_allowed).
        overlap_onset_years=on_terrane(-1.0e9),
        mobile_cover_m=np.full(count, COVER_M),
        mobile_cover_continental_m=on_terrane(COVER_M),
    )
    lithosphere.sync_plate_elevation(plate)
    return plate


def _setup(*others, carrier=None):
    a = carrier if carrier is not None else _carrier()
    world = _world(a, *others)
    world.debug_diagnostics = True
    world.reset_phase_budget()
    return world, a


def _retreat(world, a, others, layers=1.5):
    ctx = boundary_context(
        world, a, list(others), STEP_YEARS,
        lambda contested: quad_tectonics.components_of_at_least(a, contested, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN),
        node_weight=a.node_areas_m2() / lithosphere.node_area_m2(line_spacing_rad(world.node_density)),
    )
    terrane = a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    areas = a.node_areas_m2()
    before = {
        name: a.collect(name).copy()
        for name in ("crustal_thickness_m", "mantle_lithosphere_thickness_m", "continental_material_m", "craton_crust_m")
    }
    survivors = quad_tectonics._retreat(a, world, ctx, layers * line_spacing_rad(world.node_density), 1_000)
    removed = ~survivors & terrane
    assert np.any(removed)
    consumed = {name: float(values[removed] @ areas[removed]) for name, values in before.items()}
    return ctx, removed, consumed


def _terranes(world) -> dict:
    return world.suture_transfer_stats.get("terranes", {})


def _hm_sink(world, account, scope="continental_node") -> float:
    entry = world.hm_source_sink_ledger.get(account)
    return entry["scopes"][scope]["sink_m3"] if entry else 0.0


def test_cells_are_unequal():
    a = _carrier()
    areas = a.node_areas_m2()
    assert np.ptp(areas) > 1e-6 * areas.mean()


def test_a_terrane_docks_onto_the_continental_plate_it_runs_into():
    b = _plate(2, _block((20, 34), J))
    world, a = _setup(b)
    b_hc, b_hm, b_count = _volume(b, "crustal_thickness_m"), _volume(b, "mantle_lithosphere_thickness_m"), b.node_count()
    b_material = _volume(b, "continental_material_m")
    carrier_terrane_before = _volume(a, "crustal_thickness_m", a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL)

    _, removed, consumed = _retreat(world, a, [b])

    budget = world.orogenic_relief_budget
    lost = budget["suture_lower_crust_subducted_m3"]
    assert lost > 0.0
    # The upper plate takes the terrane's crust less the lost share, keeps its cells, and gains
    # no mantle lithosphere: the terrane's Hm goes down with its carrier's slab.
    assert b.node_count() == b_count
    assert _volume(b, "crustal_thickness_m") - b_hc == pytest.approx(consumed["crustal_thickness_m"] - lost, rel=1e-10)
    assert _volume(b, "mantle_lithosphere_thickness_m") == pytest.approx(b_hm, rel=1e-12)
    assert _hm_sink(world, "oceanic_and_deep_subduction") == pytest.approx(consumed["mantle_lithosphere_thickness_m"], rel=1e-12)
    assert world.hm_suture_budget["fronts"] == 0
    # Provenance follows the crust: material onto the upper plate, craton reworked there.
    lost_share = lost / consumed["crustal_thickness_m"]
    assert _volume(b, "continental_material_m") - b_material == pytest.approx(
        (1.0 - lost_share) * consumed["continental_material_m"], rel=1e-10
    )
    # Its craton lands as craton, with its date; only the lost share's sinks.
    assert world.craton_ledger["collision_reworked_m3"] == 0.0
    assert _volume(b, "craton_crust_m") == pytest.approx((1.0 - lost_share) * consumed["craton_crust_m"], rel=1e-10)
    assert world.craton_ledger["subducted_m3"] == pytest.approx(lost_share * consumed["craton_crust_m"], rel=1e-10)
    received = b.collect("craton_crust_m") > 0.0
    assert np.any(received) and np.all(b.collect("craton_formed_years")[received] == -2.0e9)
    assert world.suture_transfer_stats["terranes"]["docked_craton_m3"] == pytest.approx(consumed["craton_crust_m"], rel=1e-12)
    # Nothing stayed behind on the carrier.
    carrier_terrane_after = _volume(a, "crustal_thickness_m", a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL)
    assert carrier_terrane_after == pytest.approx(carrier_terrane_before - consumed["crustal_thickness_m"], rel=1e-12)
    stats = _terranes(world)
    assert stats["docked_fronts"] == 1 and stats["fallback_fronts"] == 0
    assert stats["contact_docking_cells"] == np.count_nonzero(removed)
    assert stats["docked_hc_m3"] == pytest.approx(consumed["crustal_thickness_m"], rel=1e-12)
    _assert_ledgers_close(world)


@pytest.mark.parametrize("neighbour", ["oceanic_plate", "oceanic_coded_margin"])
def test_a_noncontinental_neighbour_leaves_the_terrane_on_its_carrier(neighbour):
    if neighbour == "oceanic_plate":
        b = _plate(2, _block((20, 34), J), "oceanic")
    else:
        # A continental plate whose cells at the contact are oceanic-coded: no established
        # continental margin there.
        b = _plate(2, _block((20, 34), J))
        i, _ = _columns(b)
        codes = b.collect("crust_type_code")
        codes[i < 26] = CRUST_TYPE_OCEANIC
        b.set_fields_on_plate(crust_type_code=codes)
    world, a = _setup(b)
    b_hc = _volume(b, "crustal_thickness_m")
    terrane_before = _volume(a, "crustal_thickness_m", a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL)

    _retreat(world, a, [b])

    assert _volume(b, "crustal_thickness_m") == pytest.approx(b_hc, rel=1e-12)
    terrane_after = _volume(a, "crustal_thickness_m", a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL)
    assert terrane_after == pytest.approx(terrane_before, rel=1e-10)
    stats = _terranes(world)
    assert stats.get("docked_fronts", 0) == 0 and stats["fallback_fronts"] == 1
    _assert_ledgers_close(world)


def test_adjacent_fronts_dock_onto_their_own_overriders_and_not_onto_an_oceanic_one():
    # The terrane's leading edge meets three plates along strike: two continents and, between
    # them, an oceanic plate.
    a = _carrier(j_range=(10, 40))
    b = _plate(2, _block((20, 34), (10, 20)))
    o = _plate(3, _block((20, 34), (20, 30)), "oceanic")
    c = _plate(4, _block((20, 34), (30, 40)))
    world, a = _setup(b, o, c, carrier=a)
    hc = {p.plate_id: _volume(p, "crustal_thickness_m") for p in (b, o, c)}
    hc_a, areas = a.collect("crustal_thickness_m").copy(), a.node_areas_m2()

    ctx, removed, consumed = _retreat(world, a, [b, o, c])

    lost_fraction = orogeny_lost_fraction(world)
    assert np.any(removed & (ctx.inputs.neighbor_plate_id == o.plate_id))
    for plate in (b, c):
        front = removed & (ctx.inputs.neighbor_plate_id == plate.plate_id)
        expected = float(hc_a[front] @ areas[front])
        assert expected > 0.0
        gained = _volume(plate, "crustal_thickness_m") - hc[plate.plate_id]
        # Each continent takes its own stretch (less the lost share), and only that.
        assert gained == pytest.approx(expected * (1.0 - lost_fraction(front, hc_a, areas)), rel=1e-9)
    assert _volume(o, "crustal_thickness_m") == pytest.approx(hc[o.plate_id], rel=1e-12)
    stats = _terranes(world)
    assert stats["docked_fronts"] == 2 and stats["fallback_fronts"] == 1
    _assert_ledgers_close(world)


def orogeny_lost_fraction(world):
    """The lost share of a column's Hc, from the configured partition and its cover."""
    shares = crust_transfer.partition(world)

    def fraction(front, hc, areas):
        total = float(hc[front] @ areas[front])
        cover = COVER_M * float(areas[front].sum())
        return shares.loss * (total - cover) / total

    return fraction


def _full_upper(room_cells=0, room_m=0.0):
    """A continental plate at the suture cap everywhere, so a dock can only place what its
    far `room_cells` columns hold (`room_m` each)."""
    b = _plate(2, _block((20, 34), J))
    i, _ = _columns(b)
    hc = np.full(b.node_count(), SUTURE_ACCRETION_MAX_HC_M)
    if room_cells:
        hc[i >= 34 - room_cells] -= room_m
    b.set_fields_on_plate(crustal_thickness_m=hc, continental_material_m=hc.copy())
    lithosphere.sync_plate_elevation(b)
    return b


def test_partial_docking_returns_the_rest_to_the_carrier_with_its_share_of_the_column():
    b = _full_upper(room_cells=1, room_m=5_000.0)
    world, a = _setup(b)
    terrane = a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    totals = {name: _volume(a, name, terrane) + _volume(b, name) for name in ("crustal_thickness_m", "continental_material_m")}
    hm_a = _volume(a, "mantle_lithosphere_thickness_m")
    craton_before = _volume(a, "craton_crust_m") + _volume(b, "craton_crust_m")

    _, removed, consumed = _retreat(world, a, [b])

    budget = world.orogenic_relief_budget
    returned = budget["terrane_returned_m3"]
    assert 0.0 < returned < consumed["crustal_thickness_m"]
    assert budget["no_outlet_subducted_m3"] == 0.0
    # Continental crust survives except the lost share: what didn't dock is back on the carrier.
    lost = budget["suture_lower_crust_subducted_m3"]
    terrane = a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    after = _volume(a, "crustal_thickness_m", terrane) + _volume(b, "crustal_thickness_m")
    assert after == pytest.approx(totals["crustal_thickness_m"] - lost, rel=1e-10)
    # The donated, returned and partition accounts reconcile.
    hc = consumed["crustal_thickness_m"]
    assert budget["suture_scraped_m3"] + budget["suture_underthrust_m3"] + lost == pytest.approx(hc, rel=1e-12)
    # The transfer's donation excludes the returned share; the own-plate path donates it again.
    assert budget["suture_donated_m3"] == pytest.approx(hc, rel=1e-12)
    # Hm: the docked and the returned shares both go down with the carrier's slab, so none
    # piles onto the terrane's survivors.
    assert _hm_sink(world, "oceanic_and_deep_subduction") == pytest.approx(consumed["mantle_lithosphere_thickness_m"], rel=1e-10)
    assert _volume(a, "mantle_lithosphere_thickness_m") == pytest.approx(
        hm_a - consumed["mantle_lithosphere_thickness_m"], rel=1e-10
    )
    assert world.hm_suture_budget["subducted_hm_m3"] == 0.0
    stats = _terranes(world)
    # Counted once, as a partial dock, not again as a fallback.
    assert stats["partial_fronts"] == 1 and stats["fallback_fronts"] == 0 and stats["relocated_fronts"] == 0
    assert stats["returned_hc_m3"] == pytest.approx(returned)
    # The remainder's craton stays craton on the carrier: none of it is reworked.
    assert world.craton_ledger["collision_reworked_m3"] == 0.0
    assert _volume(a, "craton_crust_m") + _volume(b, "craton_crust_m") == pytest.approx(
        craton_before - world.craton_ledger["subducted_m3"], rel=1e-10
    )
    assert np.all(a.collect("mobile_cover_m") <= a.collect("crustal_thickness_m") + 1e-9)
    _assert_ledgers_close(world)


def test_a_terrane_no_overrider_can_hold_relocates_on_its_carrier_rather_than_subducting():
    # The whole terrane is the consumed column, so nothing of it survives on the carrier to
    # take the remainder: it relocates onto the carrier's oceanic footprint.
    a = _carrier(terrane_from=23)
    b = _full_upper()
    world, a = _setup(b, carrier=a)
    b_hc = _volume(b, "crustal_thickness_m")

    _, removed, consumed = _retreat(world, a, [b])

    assert _volume(b, "crustal_thickness_m") == pytest.approx(b_hc, rel=1e-12)
    budget = world.orogenic_relief_budget
    lost = budget["suture_lower_crust_subducted_m3"]
    assert budget["no_outlet_subducted_m3"] == 0.0
    terrane = a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    assert np.any(terrane)
    assert _volume(a, "crustal_thickness_m", terrane) == pytest.approx(consumed["crustal_thickness_m"] - lost, rel=1e-10)
    # It keeps the returned share of its mantle column: only the docked share's Hm sank.
    returned = budget["terrane_returned_m3"]
    hm_returned = consumed["mantle_lithosphere_thickness_m"] * returned / consumed["crustal_thickness_m"]
    assert _volume(a, "mantle_lithosphere_thickness_m", terrane) == pytest.approx(hm_returned, rel=1e-10)
    assert _hm_sink(world, "oceanic_and_deep_subduction") == pytest.approx(
        consumed["mantle_lithosphere_thickness_m"] - hm_returned, rel=1e-10
    )
    # Its craton and age come along whole.
    assert _volume(a, "craton_crust_m", terrane) == pytest.approx(consumed["craton_crust_m"] * (1.0 - lost / consumed["crustal_thickness_m"]), rel=1e-10)
    assert np.all(a.collect("craton_formed_years")[terrane] == -2.0e9)
    stats = _terranes(world)
    assert stats["partial_fronts"] == 1 and stats["relocated_fronts"] == 1
    _assert_ledgers_close(world)


def test_without_evidence_a_terrane_carrier_is_the_lower_plate():
    """The carrier tier of the polarity fallback: stationary plates leave no slab or arc
    evidence, and the oceanic carrier goes down even though it is the larger plate, which the
    size tier alone would have kept on top."""
    a = _carrier(i_range=(0, 21))
    b = _plate(2, _block((20, 26), J))
    world, a = _setup(b, carrier=a)
    assert a.node_count() > b.node_count()

    cp.observe_contacts(world, STEP_YEARS)

    (record,) = world.collision_fronts
    assert (record.source, record.fallback_basis, record.lower_plate_id) == ("fallback", "carrier", a.plate_id)
    assert world.collision_polarity_stats["fallback_carrier"] == 1


def test_a_polarized_terrane_front_docks_through_its_frozen_upper_plate():
    """End to end through the prepass: the carrier is the frozen lower plate, and its terrane
    docks onto the front's upper plate rather than by contact."""
    from unit_tests.test_collision_polarity import _drive, _plate_centroid

    b = _plate(2, _block((20, 34), J))
    a = _carrier(i_range=(10, 21))
    world, a = _setup(b, carrier=a)
    _drive(a, _plate_centroid(b), 5.0)
    cp.observe_contacts(world, STEP_YEARS)
    (record,) = world.collision_fronts
    assert record.lower_plate_id == a.plate_id
    b_hc = _volume(b, "crustal_thickness_m")

    ctx, removed, consumed = _retreat(world, a, [b])

    assert np.all(ctx.suture_upper_plate_id[removed] == b.plate_id)
    lost = world.orogenic_relief_budget["suture_lower_crust_subducted_m3"]
    assert _volume(b, "crustal_thickness_m") - b_hc == pytest.approx(consumed["crustal_thickness_m"] - lost, rel=1e-10)
    stats = _terranes(world)
    assert stats["polarized_cells"] == np.count_nonzero(removed) and stats["contact_docking_cells"] == 0
    assert stats["docked_fronts"] == 1
    _assert_ledgers_close(world)


def test_a_wholly_consumed_carrier_sends_what_no_overrider_holds_down_with_the_slab():
    """With no surviving carrier cell to return it to, the remainder is booked as subducted
    rather than dropped."""
    a = _carrier(i_range=(22, 24), terrane_from=22)
    b = _full_upper()
    world, a = _setup(b, carrier=a)
    donors = np.ones(a.node_count(), dtype=bool)
    ids = np.full(a.node_count(), b.plate_id)

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world, upper_plate_ids=ids, hm_subduct_mask=~donors)
    a.remove_cells(donors)

    budget = world.orogenic_relief_budget
    assert budget["terrane_returned_m3"] == 0.0 and budget["no_outlet_subducted_m3"] > 0.0
    _assert_ledgers_close(world)


@pytest.mark.parametrize(("upper_date", "expected"), [(-1.0e9, -2.0e9), (-3.0e9, -3.0e9)])
def test_a_receiving_craton_keeps_the_older_date(upper_date, expected):
    b = _plate(2, _block((20, 34), J))
    b.set_fields_on_plate(craton_crust_m=np.full(b.node_count(), 5_000.0), craton_formed_years=np.full(b.node_count(), upper_date))
    world, a = _setup(b)
    craton_before = b.collect("craton_crust_m").copy()

    _retreat(world, a, [b])

    received = b.collect("craton_crust_m") > craton_before
    assert np.any(received)
    assert np.all(b.collect("craton_formed_years")[received] == expected)
    assert np.all(b.collect("craton_formed_years")[~received] == upper_date)
    _assert_ledgers_close(world)


def test_a_terrane_docks_onto_a_terrane_riding_another_oceanic_plate():
    """Terranes amalgamate offshore: the overrider's contact is another carrier's terrane."""
    b = _carrier(i_range=(20, 34), terrane_from=20, plate_id=2)
    i, _ = _columns(b)
    codes = b.collect("crust_type_code")
    codes[i >= 28] = CRUST_TYPE_OCEANIC  # b is mostly terrane at the contact, ocean behind
    b.set_fields_on_plate(crust_type_code=codes)
    world, a = _setup(b)
    b_terrane = b.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    b_hc = _volume(b, "crustal_thickness_m", b_terrane)

    _, removed, consumed = _retreat(world, a, [b])

    lost = world.orogenic_relief_budget["suture_lower_crust_subducted_m3"]
    b_terrane = b.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    assert _volume(b, "crustal_thickness_m", b_terrane) - b_hc == pytest.approx(consumed["crustal_thickness_m"] - lost, rel=1e-10)
    stats = _terranes(world)
    assert stats["terrane_to_terrane_cells"] == np.count_nonzero(removed) and stats["docked_fronts"] == 1
    _assert_ledgers_close(world)

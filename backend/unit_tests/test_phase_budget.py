import numpy as np
import pytest
from app import elevation_lines, lithosphere, phase_budget
from app.world import generate_world, step_world


class _FakeWorld:
    def __init__(self, debug_diagnostics=True):
        self.debug_diagnostics = debug_diagnostics
        self.node_density = 1.0
        self.phase_budget = {}


class _FakePlate:
    def __init__(self, crust_type):
        self.crust_type = crust_type


def _codes(*values):
    return np.array(values, dtype=np.int64)

INHERIT = elevation_lines.CRUST_TYPE_INHERIT
OCEANIC = elevation_lines.CRUST_TYPE_OCEANIC
CONTINENTAL = elevation_lines.CRUST_TYPE_CONTINENTAL


def test_record_is_a_no_op_when_debug_diagnostics_is_off():
    world = _FakeWorld(debug_diagnostics=False)
    plate = _FakePlate("continental")
    phase_budget.record(
        world, plate, "some_phase",
        np.array([10.0]), np.array([20.0]), _codes(INHERIT),
        np.array([12.0]), np.array([22.0]), _codes(INHERIT),
    )
    assert world.phase_budget == {}


def test_record_is_a_no_op_for_an_empty_slice():
    world = _FakeWorld()
    plate = _FakePlate("oceanic")
    phase_budget.record(
        world, plate, "some_phase",
        np.array([]), np.array([]), _codes(),
        np.array([]), np.array([]), _codes(),
    )
    assert world.phase_budget == {}


def test_record_accumulates_the_all_scope_delta():
    world = _FakeWorld()
    plate = _FakePlate("continental")
    phase_budget.record(
        world, plate, "convergent_deformation",
        np.array([10.0, 20.0]), np.array([30.0, 40.0]), _codes(INHERIT, INHERIT),
        np.array([12.0, 21.0]), np.array([31.0, 41.0]), _codes(INHERIT, INHERIT),
    )
    totals = world.phase_budget["convergent_deformation"]
    assert totals["calls"] == 1
    all_scope = totals["scopes"]["all"]
    assert all_scope["count_before"] == 2
    assert all_scope["count_after"] == 2
    assert all_scope["sum_hc_before"] == pytest.approx(30.0)
    assert all_scope["sum_hc_after"] == pytest.approx(33.0)
    assert all_scope["sum_hm_before"] == pytest.approx(70.0)
    assert all_scope["sum_hm_after"] == pytest.approx(72.0)


def test_record_splits_by_owning_plate_crust_type():
    world = _FakeWorld()
    phase_budget.record(
        world, _FakePlate("continental"), "phase",
        np.array([10.0]), np.array([0.0]), _codes(INHERIT), np.array([10.0]), np.array([0.0]), _codes(INHERIT),
    )
    phase_budget.record(
        world, _FakePlate("oceanic"), "phase",
        np.array([5.0]), np.array([0.0]), _codes(INHERIT), np.array([5.0]), np.array([0.0]), _codes(INHERIT),
    )
    totals = world.phase_budget["phase"]
    assert totals["calls"] == 2
    assert totals["scopes"]["continental_plate"]["sum_hc_before"] == pytest.approx(10.0)
    assert totals["scopes"]["oceanic_plate"]["sum_hc_before"] == pytest.approx(5.0)
    assert totals["scopes"]["all"]["sum_hc_before"] == pytest.approx(15.0)


def test_record_splits_by_effective_node_type_not_just_plate_type():
    # A continental plate carrying one node explicitly stamped oceanic (e.g. a melted-through
    # rift node, GitHub issue #216 item 7) -- the node-type split should follow the node's own
    # code, not the plate's nominal crust_type, while the plate-type split still lumps it with
    # the rest of this continental plate's own call.
    world = _FakeWorld()
    plate = _FakePlate("continental")
    phase_budget.record(
        world, plate, "decompression_melting",
        np.array([40.0, 5.0]), np.array([80.0, 8.0]), _codes(INHERIT, INHERIT),
        np.array([40.0, 11.0]), np.array([80.0, 8.0]), _codes(INHERIT, OCEANIC),
    )
    totals = world.phase_budget["decompression_melting"]
    scopes = totals["scopes"]
    # plate-type bucket: both nodes belong to this one continental-plate call
    assert scopes["continental_plate"]["count_before"] == 2
    assert scopes["oceanic_plate"]["count_before"] == 0
    # node-type bucket: node 0 stays continental (INHERIT resolves against the plate), node 1
    # starts continental (INHERIT) but is explicitly reclassified oceanic on the "after" side
    assert scopes["continental_node"]["count_before"] == 2
    assert scopes["continental_node"]["count_after"] == 1
    assert scopes["oceanic_node"]["count_before"] == 0
    assert scopes["oceanic_node"]["count_after"] == 1
    assert scopes["continental_node"]["sum_hc_after"] == pytest.approx(40.0)
    assert scopes["oceanic_node"]["sum_hc_after"] == pytest.approx(11.0)


def test_record_handles_a_node_count_change_between_before_and_after():
    # A growth/shrink phase can add or remove nodes -- count_before/count_after (and the sums
    # behind them) must be tracked independently rather than assuming a fixed-length pair.
    world = _FakeWorld()
    plate = _FakePlate("oceanic")
    phase_budget.record(
        world, plate, "line_growth_shrink",
        np.array([10.0, 10.0, 10.0]), np.array([40.0, 40.0, 40.0]), _codes(INHERIT, INHERIT, INHERIT),
        np.array([10.0]), np.array([40.0]), _codes(INHERIT),
    )
    all_scope = world.phase_budget["line_growth_shrink"]["scopes"]["all"]
    assert all_scope["count_before"] == 3
    assert all_scope["count_after"] == 1
    assert all_scope["sum_hc_before"] == pytest.approx(30.0)
    assert all_scope["sum_hc_after"] == pytest.approx(10.0)


def test_record_accumulates_across_multiple_calls_to_the_same_phase():
    world = _FakeWorld()
    plate = _FakePlate("continental")
    for _ in range(3):
        phase_budget.record(
            world, plate, "phase",
            np.array([10.0]), np.array([0.0]), _codes(INHERIT), np.array([11.0]), np.array([0.0]), _codes(INHERIT),
        )
    totals = world.phase_budget["phase"]
    assert totals["calls"] == 3
    assert totals["scopes"]["all"]["sum_hc_before"] == pytest.approx(30.0)
    assert totals["scopes"]["all"]["sum_hc_after"] == pytest.approx(33.0)


def test_snapshot_reads_hc_hm_and_codes_from_a_real_plate():
    world = generate_world(seed=3, num_plates=4)
    plate = world.plates[0]
    spacing_rad = elevation_lines.line_spacing_rad(world.node_density)
    hc, hm, codes, area_m2 = phase_budget.snapshot(plate, spacing_rad)
    assert len(hc) == len(hm) == len(codes) == len(area_m2) == plate.node_count()
    assert hc.sum() == pytest.approx(plate.collect("crustal_thickness_m").sum())
    assert np.allclose(area_m2, lithosphere.node_area_m2(spacing_rad))


def test_record_weights_area_and_volume_by_each_nodes_own_area():
    # Issue #257: a phase that swaps one big cell for two small ones keeps the plain Hc sum
    # growing while the real volume is unchanged.
    world = _FakeWorld()
    plate = _FakePlate("continental")
    phase_budget.record(
        world, plate, "remap",
        np.array([10.0]), np.array([20.0]), _codes(INHERIT),
        np.array([10.0, 10.0]), np.array([20.0, 20.0]), _codes(INHERIT, INHERIT),
        area_before_m2=np.array([2.0]), area_after_m2=np.array([1.0, 1.0]),
    )
    all_scope = world.phase_budget["remap"]["scopes"]["all"]
    assert all_scope["sum_hc_after"] == pytest.approx(2 * all_scope["sum_hc_before"])
    assert all_scope["area_before_m2"] == pytest.approx(2.0)
    assert all_scope["area_after_m2"] == pytest.approx(2.0)
    assert all_scope["hc_volume_before_m3"] == pytest.approx(20.0)
    assert all_scope["hc_volume_after_m3"] == pytest.approx(20.0)
    assert all_scope["hm_volume_after_m3"] == pytest.approx(40.0)


def test_record_defaults_to_the_nominal_area_per_node():
    world = _FakeWorld()
    plate = _FakePlate("oceanic")
    phase_budget.record(
        world, plate, "phase",
        np.array([10.0, 30.0]), np.array([0.0, 0.0]), _codes(INHERIT, INHERIT),
        np.array([10.0]), np.array([0.0]), _codes(INHERIT),
    )
    nominal = lithosphere.node_area_m2(elevation_lines.line_spacing_rad(world.node_density))
    all_scope = world.phase_budget["phase"]["scopes"]["all"]
    assert all_scope["area_before_m2"] == pytest.approx(2 * nominal)
    assert all_scope["area_after_m2"] == pytest.approx(nominal)
    assert all_scope["hc_volume_before_m3"] == pytest.approx(40.0 * nominal)


def test_snapshot_of_a_quad_plate_carries_its_cells_own_areas():
    world = generate_world(seed=5, num_plates=5, surface="quad")
    plate = max(world.plates, key=lambda p: p.node_count())
    snap = phase_budget.snapshot(plate, elevation_lines.line_spacing_rad(world.node_density))
    assert np.array_equal(snap.area_m2, plate.surface_nodes().area_m2)
    world.debug_diagnostics = True
    phase_budget.record_snapshots(world, plate, "phase", snap, snap)
    all_scope = world.phase_budget["phase"]["scopes"]["all"]
    assert all_scope["area_before_m2"] == pytest.approx(float(snap.area_m2.sum()))
    assert all_scope["hc_volume_before_m3"] == pytest.approx(float(np.dot(snap.area_m2, snap.hc)))


def test_quad_deform_phases_book_each_cells_own_area():
    # Issue #257: quad deformation runs over the whole plate's columns at once, and the
    # budget weights them by the cells' areas, not the nominal one.
    world = generate_world(seed=5, num_plates=5, node_density=0.5, surface="quad")
    world.debug_diagnostics = True
    world.reset_phase_budget()
    step_world(world, years=100_000)
    scope = world.phase_budget["convergent_deformation"]["scopes"]["all"]
    nominal = lithosphere.node_area_m2(elevation_lines.line_spacing_rad(world.node_density))
    assert scope["count_before"] > 0
    assert scope["area_before_m2"] != pytest.approx(scope["count_before"] * nominal, rel=1e-3)

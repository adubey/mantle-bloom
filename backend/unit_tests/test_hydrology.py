import numpy as np
import pytest
from app import geometry
from app.elevation_lines import line_spacing_rad
from app.hydroclimate import breaching, erosion, hydrology, lakes
from app.plates import collect_all_accounting_areas_m2
from app.world import World, generate_world, step_world

from .quad_fixtures import quad_plate


def _row(count: int) -> np.ndarray:
    """Column offsets for a row of `count` cells centred on the face."""
    return np.arange(count) - count // 2


def _flow_line_plate(plate_id, theta, elevation):
    """A continental row of cells, one per entry of `theta` (only its length matters)."""
    return quad_plate(plate_id, "continental", columns=_row(len(theta)), elevation=np.asarray(elevation, dtype=float))


def _flow_line_plate_with_lake(plate_id, theta, elevation, lake_depth, glacier_depth=None):
    kw = {}
    if glacier_depth is not None:
        kw["glacier_depth"] = np.asarray(glacier_depth, dtype=float)
    return quad_plate(
        plate_id, "continental", columns=_row(len(theta)), elevation=np.asarray(elevation, dtype=float),
        lake_depth=np.asarray(lake_depth, dtype=float), **kw,
    )


def _normalize(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def _river_inspector_fields():
    # A hand-built 9-node HydrologyFields for the River Inspector's grouping logic, not run
    # through compute_hydrology (small enough to hit its FLOW_NEIGHBOR_COUNT bail-out path
    # anyway). Three separate is_river-classified networks, one per mouth_type:
    #   Network A -- a Y-confluence: 0 and 1 (two headwaters) both flow into 2, which flows
    #   into 3 (the mouth, highest fallback discharge in the group), which flows onward into ocean
    #   node 4. mouth_type "ocean".
    #   Network B -- a single lone node 5, flowing into node 7, a real lake (lake_depth above
    #   the visible threshold). Node 7 itself is deliberately *not* is_river -- exactly the real
    #   invariant compute_hydrology enforces (a lake node is never classified is_river, see
    #   hydrology.compute_hydrology's own `~is_lake` term), unlike an earlier version of this
    #   fixture that fabricated is_river=True directly on the flooded node -- a state
    #   compute_hydrology can never actually produce, which had been masking a real group_rivers
    #   bug (mouth_type read the mouth's own, structurally-always-dry lake_depth instead of the
    #   target it flows into). mouth_type "lake".
    #   Network C -- a single lone node 8, a genuine unresolved dry sink (flow_target -1, no
    #   lake there yet). mouth_type "other".
    #   Network D -- a single lone node 9, ending at dry node 10 (still-dry basin floor, well
    #   above the lake's current shoreline) which itself flows onward one more hop into the
    #   same lake node 7 -- reproducing a real basin sitting near its own rim, with only its
    #   lowest node or two actually underwater: mouth_target alone reads as neither ocean nor
    #   lake, but the chain still reaches real lake water. mouth_type "lake" (confirmed
    #   directly against a real save -- see group_rivers's own docstring).
    # Node 6 is ordinary land with flow but deliberately *not* in the is_river mask, to check
    # grouping doesn't pull in non-river neighbors.
    points = np.array(
        [
            _normalize([1.0, 0.05, 0.0]),  # 0: headwater
            _normalize([1.0, -0.05, 0.02]),  # 1: headwater
            _normalize([1.0, 0.0, 0.03]),  # 2: confluence
            _normalize([1.0, 0.0, 0.06]),  # 3: mouth
            _normalize([1.0, 0.0, 0.09]),  # 4: ocean
            _normalize([-1.0, 0.0, 0.0]),  # 5: lone river node, ends at a lake
            _normalize([0.0, 1.0, 0.0]),  # 6: non-river land
            _normalize([-1.0, 0.02, 0.0]),  # 7: the lake node 5 flows into (not itself a river)
            _normalize([0.0, -1.0, 0.0]),  # 8: lone river node, a still-dry interior sink
            _normalize([0.0, -1.0, 0.02]),  # 9: lone river node, mouth of network D
            _normalize([0.0, -1.0, 0.04]),  # 10: dry basin floor between node 9 and the lake
        ]
    )
    n = len(points)
    elevation = np.array([100.0, 110.0, 80.0, 50.0, -10.0, 60.0, 40.0, 55.0, 70.0, 65.0, 58.0])
    is_ocean = np.array([False, False, False, False, True, False, False, False, False, False, False])
    flow_target = np.array([2, 2, 3, 4, -1, 7, -1, -1, -1, 10, 7])
    flow_accum = np.array([10.0, 12.0, 25.0, 40.0, 0.0, 5.0, 0.0, 0.0, 3.0, 4.0, 0.0])
    is_river = np.array([True, True, True, True, False, True, False, False, True, True, False])
    lake_depth = np.zeros(n)
    lake_depth[7] = 5.0  # node 7 well above LAKE_MIN_VISIBLE_DEPTH_M
    zeros = np.zeros(n)
    return hydrology.HydrologyFields(
        points=points,
        elevation=elevation,
        is_ocean=is_ocean,
        neighbor_idx=np.zeros((n, 1), dtype=np.int64),
        flow_target=flow_target,
        flow_accum=flow_accum,
        water_deposited=zeros.copy(),
        filled_elevation=zeros.copy(),
        spill_target=np.full(n, -1, dtype=np.int64),
        is_river=is_river,
        lake_depth=lake_depth,
        glacier_depth=zeros.copy(),
        plates_in_order=[],
    )


def test_compute_flow_direction_sink_redirects_once_lake_reaches_spill_point():
    # Same 4-node setup as above. Node 0 has no lower neighbor (its only neighbor, node 1,
    # is higher) so it's a sink candidate -- unless its current water surface has already
    # reached the basin's true spill point (50, from the test above).
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    filled = np.array([50.0, 50.0, 20.0, -10.0])
    spill = np.array([2, 2, 3, -1])

    not_frozen = np.zeros(4, dtype=bool)
    not_yet_full, not_yet_spilling = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.zeros(4), filled, spill, np.zeros(4), not_frozen
    )
    assert not_yet_full[0] == -1  # water surface 30 < 50, still a plain sink
    assert not not_yet_spilling[0]

    now_full, now_spilling = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.array([25.0, 0.0, 0.0, 0.0]), filled, spill, np.zeros(4), not_frozen
    )
    assert now_full[0] == 2  # water surface 30+25=55 >= 50 -- redirects to spill_target
    assert now_spilling[0]


def test_compute_flow_direction_a_wet_sink_spills_over_a_rim_at_its_own_floor():
    # Node 0's rim is notched down to its own centre (filled == elevation). Holding a little
    # water, its hypsometric water surface sits 0.7 m below that centre. It still spills: a
    # floor-level rim holds no water at all (#503, a breached pit that held rivers forever).
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    filled = np.array([30.0, 50.0, 20.0, -10.0])
    spill = np.array([2, 2, 3, -1])

    flow_target, should_spill = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.array([-0.7, 0.0, 0.0, 0.0]), filled, spill, np.zeros(4), np.zeros(4, dtype=bool)
    )
    assert flow_target[0] == 2
    assert should_spill[0]


def _spill_basins(elevation, forest, governing, evaporation=0.0, frozen=None, neighbor_idx=None):
    n = len(elevation)
    if neighbor_idx is None:  # a chain: each node's neighbours are the nodes either side of it
        neighbor_idx = np.stack([np.maximum(np.arange(n) - 1, 0), np.minimum(np.arange(n) + 1, n - 1)], axis=1)
    return hydrology.SpillBasins(
        forest, governing, lakes.CellHypsometry(elevation, np.ones(n), np.zeros(n), 1.0), neighbor_idx,
        np.zeros(n), np.full(n, evaporation), np.zeros(n, dtype=bool) if frozen is None else frozen,
    )


def _route_for(elevation, is_ocean, source):
    return lambda targets: hydrology.route_downstream(elevation, is_ocean, targets, source)


def test_spill_within_step_carries_a_river_through_a_chain_of_dry_pits():
    # 5 -> pit B (4) -> 3 -> pit A (2) -> 1 -> ocean. Both pits start dry, so the lagged spill
    # test holds the river in B. This step's 1,000 m^3/yr fills each pit's few m^3 of room many
    # times over, so both open and the river reaches the sea within the step (#503).
    elevation = np.array([-10.0, 5.0, 3.0, 10.0, 6.0, 20.0])
    is_ocean = np.array([True, False, False, False, False, False])
    pit_a = lakes.Lake(0, np.array([2, 3]), 3.0, 3.0, 5.0, outlet_node_idx=2, sink_node_idx=2, outlet_target_idx=1)
    pit_b = lakes.Lake(1, np.array([4, 5]), 6.0, 6.0, 10.0, outlet_node_idx=4, sink_node_idx=4, outlet_target_idx=3)
    flow_target = np.array([-1, 0, -1, 2, -1, 4])
    source = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1000.0])
    no_spill = np.zeros(6, dtype=bool)

    targets, spilling, _, _, through, _ = hydrology._spill_within_step(
        _spill_basins(elevation, [pit_a, pit_b], [pit_a, pit_b]), flow_target, no_spill,
        elevation.copy(), np.full(6, -1), _route_for(elevation, is_ocean, source), years=1000.0,
    )
    assert targets[4] == 3 and targets[2] == 1
    assert spilling[4] and spilling[2]
    assert through[1] == pytest.approx(1000.0)

    # In a desert, a lake at B's rim would evaporate more than the river brings: B holds it.
    targets, spilling, _, _, through, _ = hydrology._spill_within_step(
        _spill_basins(elevation, [pit_a, pit_b], [pit_a, pit_b], evaporation=2000.0), flow_target, no_spill,
        elevation.copy(), np.full(6, -1), _route_for(elevation, is_ocean, source), years=1000.0,
    )
    assert targets[4] == -1 and not spilling[4]
    assert through[1] == 0.0


def test_spill_within_step_unites_sibling_basins_that_spill_into_each_other():
    # Sibling pits A (node 2) and B (node 3) both reached their shared saddle, so each spills
    # into the other and the water circles. United at their parent, the pair either spills once
    # at the parent's rim (to node 1, then the sea) or holds the water as one lake (#503).
    elevation = np.array([-10.0, 5.0, 2.0, 3.0])
    is_ocean = np.array([True, False, False, False])
    pit_a = lakes.Lake(0, np.array([2]), 2.0, 2.0, 4.0, outlet_node_idx=2, sink_node_idx=2, outlet_target_idx=3)
    pit_b = lakes.Lake(1, np.array([3]), 3.0, 3.0, 4.0, outlet_node_idx=3, sink_node_idx=3, outlet_target_idx=2)
    pair = lakes.Lake(2, np.array([2, 3]), 2.0, 4.0, 5.0, children=[pit_a, pit_b], outlet_node_idx=2, outlet_target_idx=1)
    circling = np.array([-1, 0, 3, 2])
    spilling = np.array([False, False, True, True])
    source = np.array([0.0, 0.0, 100.0, 100.0])

    def run(years):
        return hydrology._spill_within_step(
            _spill_basins(elevation, [pair], [pit_a, pit_b]), circling, spilling,
            np.array([-10.0, 5.0, 4.0, 4.0]), np.array([-1, 0, 3, 2]), _route_for(elevation, is_ocean, source), years=years,
        )

    targets, now_spilling, filled, spill_to, through, _ = run(years=1000.0)  # 200,000 m^3 >> 5 m^3 of room
    assert hydrology.flow_cycles(targets) == []
    assert targets[2] == 1 and targets[3] == 1
    assert list(filled[2:]) == [5.0, 5.0] and list(spill_to[2:]) == [1, 1]
    assert through[1] == pytest.approx(200.0)

    targets, now_spilling, _, _, through, settled = run(years=0.01)  # 2 m^3 < 5 m^3 of room
    assert hydrology.flow_cycles(targets) == []
    assert targets[2] == -1 and targets[3] == -1 and not now_spilling[2:].any()
    assert settled[2:].sum() == pytest.approx(200.0)


def test_spill_within_step_sends_a_channel_back_down_the_catchment_it_belongs_to():
    # Pit A (node 1) drains at its floor over to node 2, in basin B's catchment. Node 2's
    # steepest descent is B's outlet path (node 4, then the sea), but an old channel leads it to
    # node 3, which drains back into A: a cycle across two separate root basins, which share no
    # ancestor to unite in. Node 2 follows its steepest descent instead, and the water reaches
    # the sea (#503).
    elevation = np.array([-10.0, 4.0, 6.0, 5.0, 2.0])
    is_ocean = np.array([True, False, False, False, False])
    neighbor_idx = np.array([[1, 4], [2, 3], [3, 4], [1, 2], [2, 0]])
    pit_a = lakes.Lake(0, np.array([1, 3]), 4.0, 4.0, 4.0, outlet_node_idx=1, sink_node_idx=1, outlet_target_idx=2)
    basin_b = lakes.Lake(1, np.array([2, 4]), 2.0, 2.0, 6.0, outlet_node_idx=4, sink_node_idx=4, outlet_target_idx=0)
    circling = np.array([-1, 2, 3, 1, 0])
    source = np.array([0.0, 100.0, 0.0, 0.0, 0.0])

    targets, _, _, _, through, _ = hydrology._spill_within_step(
        _spill_basins(elevation, [pit_a, basin_b], [pit_a, basin_b], neighbor_idx=neighbor_idx), circling,
        np.array([False, True, False, False, False]), np.array([-10.0, 4.0, 6.0, 5.0, 2.0]), np.array([-1, 2, -1, -1, -1]),
        _route_for(elevation, is_ocean, source), years=1000.0,
    )
    assert hydrology.flow_cycles(targets) == []
    assert targets[2] == 4
    assert through[4] == pytest.approx(100.0)


def test_flow_cycles_lists_each_cycle_once():
    flow_target = np.array([1, 2, 0, 2, -1, 6, 5, 4])
    cycles = sorted(sorted(c.tolist()) for c in hydrology.flow_cycles(flow_target))
    assert cycles == [[0, 1, 2], [5, 6]]
    assert hydrology.flow_cycles(np.array([-1, 0, 1])) == []


def test_compute_flow_direction_frozen_node_blocks_flow_even_on_unobstructed_terrain():
    # Node 2 (elevation 20) has a real downhill neighbor (node 3, the ocean) with nothing
    # blocking it -- filled_elevation[2] is exactly its own elevation, so if the frozen
    # override were folded into the same is_sink/should_spill check the redirect above uses,
    # forcing node 2's flow_target to -1 would make it trivially satisfy
    # water_surface(20) >= filled_elevation(20) and get redirected right back onto its own
    # ordinary downhill neighbor via spill_target -- silently undoing the freeze. The override
    # must instead apply to the final flow_target, after that check has already run.
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    filled = np.array([50.0, 50.0, 20.0, -10.0])
    spill = np.array([2, 2, 3, -1])
    is_frozen = np.array([False, False, True, False])

    flow_target, should_spill = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.zeros(4), filled, spill, np.zeros(4), is_frozen
    )
    assert flow_target[2] == -1  # frozen -- blocked, not redirected to its own ordinary neighbor
    assert not should_spill[2]  # should_spill itself still reflects genuine basin-spill only


def test_route_downstream_accumulates_along_a_chain():
    # A straight 4-node downhill chain (0 -> 1 -> 2 -> ocean-adjacent 3), each node
    # contributing 1.0 of its own -- accumulated flux at each node should be its own
    # contribution plus everything upstream of it.
    elevation = np.array([40.0, 30.0, 20.0, -5.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 3, -1])
    source = np.array([1.0, 1.0, 1.0, 0.0])

    through_flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source)
    assert through_flux.tolist() == [1.0, 2.0, 3.0, 0.0]
    assert deposited[3] == 3.0  # everything reaching the ocean-adjacent node's target settles there
    assert deposited[0] == 0.0 and deposited[1] == 0.0 and deposited[2] == 0.0


def test_route_downstream_retain_fraction_deposits_partway():
    # Same chain, but node 1 retains half of whatever passes through it.
    elevation = np.array([40.0, 30.0, 20.0, -5.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 3, -1])
    source = np.array([1.0, 1.0, 1.0, 0.0])
    retain = np.array([0.0, 0.5, 0.0, 0.0])

    through_flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source, retain_fraction=retain)
    # Node 1 receives 1.0 (its own) + 1.0 (from node 0) = 2.0, retains half (1.0) locally,
    # passes 1.0 onward.
    assert deposited[1] == 1.0
    assert through_flux[1] == 1.0
    # Node 2 receives 1.0 (its own) + 1.0 (passed on from node 1) = 2.0, all continues on.
    assert through_flux[2] == 2.0
    assert deposited[3] == 2.0
    # Conservation: total deposited equals total source.
    assert np.isclose(deposited.sum(), source.sum())


def test_route_downstream_loss_fraction_evaporates_without_depositing():
    # Same chain, but node 1 evaporates half of whatever passes through it -- unlike
    # retain_fraction, that half must simply vanish (not show up in `deposited` anywhere),
    # since it left as atmospheric moisture, not as material staying in the channel.
    elevation = np.array([40.0, 30.0, 20.0, -5.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 3, -1])
    source = np.array([1.0, 1.0, 1.0, 0.0])
    loss = np.array([0.0, 0.5, 0.0, 0.0])

    through_flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source, loss_fraction=loss)
    # Node 1 receives 1.0 (its own) + 1.0 (from node 0) = 2.0, evaporates half (1.0), passes
    # 1.0 onward -- none of the evaporated half is deposited anywhere.
    assert through_flux[1] == 1.0
    assert deposited[1] == 0.0
    assert through_flux[2] == 2.0  # node 2's own 1.0 + the 1.0 that made it past node 1
    assert deposited[3] == 2.0
    # Conservation is broken on purpose here: 1.0 of the original 3.0 total source evaporated
    # away instead of ending up in `deposited`.
    assert np.isclose(deposited.sum(), source.sum() - 1.0)


def test_connected_ocean_mask_excludes_an_enclosed_interior_pit():
    # 40 nodes along an arc: a wide ocean (0-14), a thick land ridge (15-27), a small
    # sub-sea-level pocket (28-30), then more land. The ridge is far wider than the k-NN
    # neighbour count, so no below-sea-level edge ever bridges it.
    theta = 0.05 * np.arange(40)
    xyz = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0   # the world ocean
    elevation[28:31] = -60.0   # an enclosed interior depression

    mask = hydrology.connected_ocean_mask(xyz, elevation, sea_level_m=0.0)

    assert mask[0:15].all()          # the ocean is ocean
    assert not mask[28:31].any()     # the enclosed pit is not
    assert not mask[15:28].any() and not mask[31:].any()  # dry land never is


def test_connected_ocean_mask_leaves_a_single_connected_ocean_untouched():
    theta = 0.05 * np.arange(20)
    xyz = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)
    elevation = np.full(20, 100.0)
    elevation[0:12] = -200.0  # one contiguous below-sea-level body

    mask = hydrology.connected_ocean_mask(xyz, elevation, sea_level_m=0.0)
    assert mask.tolist() == [True] * 12 + [False] * 8


def test_compute_hydrology_treats_an_interior_pit_as_an_endorheic_basin_that_silts_in():
    d = 0.05
    theta = d * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    fields = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=1_000_000)

    assert fields.is_ocean[0:15].all()
    assert not fields.is_ocean[28:31].any()  # classified land, not a vanishing ocean sink
    # Rain routed into the closed pit pools there and drops sediment onto the wet floor.
    assert fields.lake_depth[28:31].max() > 0.0
    assert fields.silt_deposited[28:31].max() > 0.0


def test_compute_hydrology_breaches_a_shallow_pit_instead_of_ponding_it(monkeypatch):
    # A 20 m deep pit on a 100 m plateau, in a wet climate (2,000 mm/yr at 10 C), so its
    # runoff outpaces a rim-level lake's evaporation. One 100 kyr step carves 30 m of ordinary
    # rock, so with breaching on the rim is notched and the pit drains. With it off, rain
    # ponds there.
    d = 0.05
    theta = d * np.arange(40)
    elevation = np.full(40, 100.0)
    elevation[0:15] = -100.0
    elevation[28:31] = 80.0

    def run(precipitation=2000.0, temperature=10.0, lake_depth=None):
        plate = (
            _flow_line_plate(0, theta, elevation)
            if lake_depth is None
            else _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
        )
        world = World(seed=0, plates=[plate])
        return hydrology.compute_hydrology(world, np.full(40, precipitation), np.full(40, temperature), years=100_000)

    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", False)
    closed = run()
    assert closed.lake_depth[28:31].max() > 0.0
    assert closed.breach_notch_m.shape == (0,)

    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", True)
    breached = run()
    # The notch is cut to the pit's floor and leads lower, so even the low ground inside the
    # pit's cells drains through it (lakes.Lake.drains_at_floor).
    assert breached.lake_depth[28:31].max() == 0.0
    # Closed, the same water keeps rising toward the 100 m rim on the next step, once the
    # pit's cells have merged into one lake.
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", False)
    closed_again = run(lake_depth=closed.lake_depth)
    assert closed_again.lake_depth[28:31].min() > 10.0
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", True)
    assert breached.breach_notch_m.max() == pytest.approx(20.0)
    assert np.all(breached.filled_elevation[28:31] <= 80.0)

    # The same pit in a semi-arid climate (800 mm/yr at 15 C) would level off below its rim,
    # so it stays a closed endorheic basin.
    arid = run(precipitation=800.0, temperature=15.0)
    assert arid.lake_depth[28:31].max() > 0.0
    assert not arid.breach_notch_m.any()


def test_compute_hydrology_populates_is_sea_once_a_basin_stays_flooded():
    # Same interior-pit fixture as the test above. is_sea is derived from *last* step's own
    # flooded extent (lakes._classify_tier), so it can't promote on the very first pass -- but
    # this tiny synthetic world's implied node_area_km2 (whole-sphere area / 40 nodes) is
    # enormous, so a single wet node held over from a prior step already clears
    # SEA_MIN_FLOODED_AREA_KM2 on the next one.
    d = 0.05
    theta = d * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    first = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=1_000_000)
    assert first.lake_depth[28:31].max() > 0.0
    assert first.is_sea.shape == first.lake_depth.shape
    assert not first.is_sea.any()  # nothing was wet *last* step yet -- can't have promoted

    second_plate = _flow_line_plate_with_lake(0, theta, elevation, first.lake_depth)
    second_world = World(seed=0, plates=[second_plate])
    second = hydrology.compute_hydrology(second_world, np.full(40, 800.0), np.full(40, 15.0), years=1_000_000)
    assert second.lake_depth[28:31].max() > 0.0
    assert second.is_sea[28:31].any()
    assert not second.is_sea[0:15].any()  # ocean itself is never "sea tier" -- is_sea is lake-only


def test_compute_hydrology_end_to_end_on_a_small_synthetic_world():
    # 12 points along a monotonically descending line, the last three underwater --
    # exercises the whole pipeline (basin-spill, flow direction, accumulation, river
    # classification) on a real (if small) node cloud built the same way as the rest of the
    # test suite's synthetic-plate fixtures.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)  # last three points go underwater
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    precipitation_at_nodes = np.full(12, 800.0)
    temperature_at_nodes = np.full(12, 15.0)  # well above freezing -- no glaciation involved
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes, temperature_at_nodes, years=1_000_000)

    assert len(fields.points) == 12
    assert fields.is_ocean.tolist() == [False] * 9 + [True] * 3
    # Every land node should have a full path to the ocean on this simple monotonic line.
    assert np.all(np.isfinite(fields.filled_elevation[~fields.is_ocean]))
    # flow_accum should be non-decreasing heading downhill along this single chain (each
    # node picks up whatever accumulated above it, net of nothing lost here).
    land_accum = fields.flow_accum[~fields.is_ocean]
    assert np.all(np.diff(land_accum) >= -1e-6)


def test_lake_component_sizes_groups_connected_lake_nodes_and_ignores_dry_ones():
    # 0-1-2 chain all flooded, node 3 (also flooded) reachable only from node 1 via a
    # one-directional edge (1 lists 3 as a neighbor, 3 does not list 1 back) -- confirms
    # sizing treats neighbor_idx as undirected, not just whatever direction happens to be
    # listed. Node 4 is dry and must read 0 regardless of being graph-adjacent to lake nodes.
    is_lake = np.array([True, True, True, True, False])
    neighbor_idx = np.array([[1, 1], [2, 3], [1, 1], [1, 1], [2, 2]])

    sizes = hydrology._lake_component_sizes(is_lake, neighbor_idx)
    assert sizes.tolist() == [4, 4, 4, 4, 0]


def test_lake_component_sizes_splits_disconnected_lakes():
    # Two separate two-node ponds (0-1 and 2-3), never connected to each other.
    is_lake = np.array([True, True, True, True])
    neighbor_idx = np.array([[1], [0], [3], [2]])

    sizes = hydrology._lake_component_sizes(is_lake, neighbor_idx)
    assert sizes.tolist() == [2, 2, 2, 2]


def test_compute_hydrology_spilling_lake_erodes_its_outlet_from_its_own_surface_area():
    # A 22-node line: an upland ramp (0-9) feeding a sink (10) sitting behind a barrier
    # (11-13, peaking at 300) that descends to the ocean (19-21). filled_elevation[10] works
    # out to 200 (the true minimal-bottleneck crossing, confirmed directly against this exact
    # layout) -- lake_depth[10]=160 puts the sink's water surface at 50+160=210, over that
    # cap, so it's actively spilling this step.
    d = 0.002
    theta = d * np.arange(22)
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    lake_depth = np.zeros(22)
    lake_depth[10] = 160.0

    plate = _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
    world = World(seed=0, plates=[plate])

    precipitation_at_nodes = np.full(22, 10.0)  # small -- the breach term should dominate
    temperature_at_nodes = np.full(22, 15.0)
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes, temperature_at_nodes, years=1_000_000)

    # A single-node lake (only node 10 is above LAKE_MIN_VISIBLE_DEPTH_M), so the breach term
    # there is exactly one component-size's worth of LAKE_BREACH_EROSION_COEFFICIENT --
    # dwarfing the handful of small precip contributions reaching the same node.
    assert fields.flow_accum[10] >= hydrology.LAKE_BREACH_EROSION_COEFFICIENT
    assert fields.flow_accum[10] < hydrology.LAKE_BREACH_EROSION_COEFFICIENT * 1.1


def test_compute_hydrology_spilling_lake_also_erodes_its_own_rim_when_that_is_the_higher_side():
    # The multi-node counterpart to the test above: a two-node basin (floor 10, shore/rim 11)
    # behind a barrier, engineered so the basin's *own* rim member (11, elevation 210) is the
    # higher side of the governing boundary edge -- higher than the far side outside the lake
    # (15, elevation 200) that `max_depth` would otherwise appear to come from (confirmed
    # directly: `lakes.build_lake_hierarchy` on this exact layout gives max_depth=210 via
    # outlet_node_idx=11, not 200 via outlet_target_idx=15). Before the rim_node_idx fix, the
    # breach surge only ever rode `spill_target` from the sink (10) to the far side (15) --
    # bypassing 11 entirely -- so the *actual* dam (11) never got a drop of erosive flow no
    # matter how long the lake spilled. lake_depth[10]=170 puts the sink's water surface at
    # 50+170=220, over the 210 cap, so it's actively spilling this step.
    d = 0.002
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0, 210.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0, -150.0, -200.0, -250.0, -300.0]
    )
    n = len(elevation)
    theta = d * np.arange(n)
    lake_depth = np.zeros(n)
    lake_depth[10] = 170.0
    plate = _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
    world = World(seed=0, plates=[plate])

    precipitation_at_nodes = np.full(n, 10.0)  # small -- the breach term should dominate
    temperature_at_nodes = np.full(n, 15.0)
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes, temperature_at_nodes, years=1_000_000)

    # Both sides of the true bottleneck now carry roughly one component-size's worth of
    # LAKE_BREACH_EROSION_COEFFICIENT -- the lake's own rim (11, previously untouched) just as
    # much as the far side (15, which already worked before this fix).
    assert fields.flow_accum[11] >= hydrology.LAKE_BREACH_EROSION_COEFFICIENT
    assert fields.flow_accum[11] < hydrology.LAKE_BREACH_EROSION_COEFFICIENT * 1.1
    assert fields.flow_accum[15] >= hydrology.LAKE_BREACH_EROSION_COEFFICIENT


def test_breach_surge_drives_flow_accum_but_never_physical_discharge(monkeypatch):
    # The single-node spilling lake above. Its surge rides the legacy flow_accum (erosion's
    # input) but is kept out of discharge_m3_s (issue #355): turning the surge off leaves the
    # discharge exactly unchanged.
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    n = len(elevation)
    lake_depth = np.zeros(n)
    lake_depth[10] = 160.0

    def run():
        world = World(seed=0, plates=[_flow_line_plate_with_lake(0, 0.002 * np.arange(n), elevation, lake_depth)])
        return hydrology.compute_hydrology(world, np.full(n, 10.0), np.full(n, 15.0), years=1_000_000)

    surged = run()
    assert surged.breach_surge[10] == pytest.approx(hydrology.LAKE_BREACH_EROSION_COEFFICIENT)
    monkeypatch.setattr(hydrology, "LAKE_BREACH_EROSION_COEFFICIENT", 0.0)
    calm = run()
    assert np.all(calm.breach_surge == 0.0)
    assert surged.flow_accum[10] > calm.flow_accum[10] + 0.9 * surged.breach_surge[10]
    assert np.array_equal(surged.discharge_m3_s, calm.discharge_m3_s)
    assert surged.discharge_m3_s[10] > 0.0  # the lake's own runoff still passes through


def test_breach_surge_never_sizes_lake_silt(monkeypatch):
    # A spilling lake (sink 10, behind the 11-13 barrier) whose overflow runs into a pit (16)
    # before the ocean. The surge settles in the pit's legacy water_deposited, but silt is
    # sized from the physical inflow (issue #363), so turning the surge off leaves every
    # lake's silt unchanged. Breaching is off so the pit holds water instead of draining.
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0]
        + [40.0]
        + [120.0, 130.0, 140.0, 130.0, 60.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    n = len(elevation)
    lake_depth = np.zeros(n)
    lake_depth[10] = 160.0
    lake_depth[16] = 5.0
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", False)

    def run():
        world = World(seed=0, plates=[_flow_line_plate_with_lake(0, 0.002 * np.arange(n), elevation, lake_depth)])
        return hydrology.compute_hydrology(world, np.full(n, 2000.0), np.full(n, 5.0), years=1_000_000)

    surged = run()
    assert surged.breach_surge[10] > 0.0
    monkeypatch.setattr(hydrology, "LAKE_BREACH_EROSION_COEFFICIENT", 0.0)
    calm = run()
    assert surged.water_deposited[16] > calm.water_deposited[16]  # the surge reached the pit
    assert surged.silt_deposited[16] > 0.0 and surged.silt_deposited[10] > 0.0
    assert np.array_equal(surged.silt_deposited, calm.silt_deposited)


def test_a_spilling_sink_loses_the_same_transit_share_as_any_node():
    # The single-node spilling lake above. Physical water has one transit-loss rule (issue
    # #355): the spilling sink loses its share like any other node. Only the legacy routing
    # exempts it, for the breach surge's sake.
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    n = len(elevation)
    lake_depth = np.zeros(n)
    lake_depth[10] = 160.0
    world = World(seed=0, plates=[_flow_line_plate_with_lake(0, 0.002 * np.arange(n), elevation, lake_depth)])
    fields = hydrology.compute_hydrology(world, np.full(n, 10.0), np.full(n, 15.0), years=1_000_000)

    assert fields.breach_surge[10] > 0.0
    loss = hydrology.river_transit_loss_fraction(np.array([15.0]), np.array([False]))[0]
    assert loss > 0.0
    runoff_m3_s = breaching.runoff_mm(np.array([10.0]), breaching.open_water_evaporation_mm(np.array([15.0])))[0] / 1000.0
    runoff_m3_s *= _node_areas(world)[10] / hydrology.SECONDS_PER_YEAR
    inflow_m3_s = float(fields.discharge_m3_s[fields.flow_target == 10].sum())
    assert inflow_m3_s > 0.0
    assert fields.discharge_m3_s[10] == pytest.approx((inflow_m3_s + runoff_m3_s) * (1.0 - loss), rel=1e-12)


def test_frozen_river_ice_is_the_physical_water_settling_at_frozen_sinks():
    # Warm uplands drain into a frozen band: each frozen land node is a dead end, so the
    # physical water reaching it settles there and freezes at the glacier exchange rate
    # (WATER_TO_ICE, issue #355). The ice per year doesn't depend on step length, and
    # the breach surge never reaches it.
    theta = 0.002 * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    temperature = np.where(np.arange(12) < 6, 15.0, -10.0)
    per_year = []
    for years in (100_000, 1_000_000):
        world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
        fields = hydrology.compute_hydrology(world, np.full(12, 800.0), temperature, years=years)
        frozen_land = (temperature < 0.0) & ~fields.is_ocean
        assert np.all(fields.flow_target[frozen_land] == -1)
        settled_m3 = float(fields.discharge_m3_s[frozen_land].sum()) * hydrology.SECONDS_PER_YEAR * years
        assert settled_m3 > 0.0
        ledger = fields.ice.ledger
        assert ledger.frozen_river_capped_m3 == 0.0
        assert ledger.frozen_river_ice_m3 == pytest.approx(hydrology.WATER_TO_ICE * settled_m3, rel=1e-9)
        per_year.append(ledger.frozen_river_ice_m3 / years)
    assert per_year[1] == pytest.approx(per_year[0], rel=1e-6)  # snowmelt on the cold band varies a hair


def test_discharge_carries_the_ice_ledgers_meltwater_exactly_once_at_any_step_length():
    # An all-land chain under 200 m of melting ice: every node is frozen, so each is its own
    # sink and the discharge summed over nodes is each node's own meltwater. That's the
    # ledger's meltwater (surface and basal melt, overflow, capped ice, bare-ground snow; not
    # calving) spread over the step, in m^3/s (issue #355).
    theta = 0.002 * np.arange(12)
    elevation = 900.0 - 60.0 * np.arange(12)
    for years in (100_000, 1_000_000):
        world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, np.zeros(12), glacier_depth=np.full(12, 200.0))])
        fields = hydrology.compute_hydrology(
            world, np.full(12, 100.0), np.full(12, -6.0), years=years, seasonal_amplitude_at_nodes=np.full(12, 18.0)
        )
        ledger = fields.ice.ledger
        assert not np.any(fields.is_ocean) and np.all(fields.flow_target == -1)
        assert ledger.surface_melt_m3 > 0.0 and ledger.basal_melt_m3 > 0.0
        routed_m3 = float(fields.discharge_m3_s.sum()) * hydrology.SECONDS_PER_YEAR * years
        assert routed_m3 == pytest.approx(ledger.meltwater_m3, rel=1e-9)


def test_discharge_from_fixed_rainfall_does_not_depend_on_step_length():
    # A warm, ice-free chain to the sea: runoff and the transit loss are both per-year rates,
    # so the same rainfall gives the same discharge at a 100 kyr and a 1 Myr step (issue #355).
    # The legacy flow_accum still loses more water to evaporation on the longer step.
    theta = 0.002 * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)

    def run(years):
        world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
        return hydrology.compute_hydrology(world, np.full(12, 800.0), np.full(12, 20.0), years=years)

    short, long = run(100_000), run(1_000_000)
    land = ~short.is_ocean
    assert np.max(short.discharge_m3_s[land]) > 0.0
    assert long.discharge_m3_s == pytest.approx(short.discharge_m3_s, rel=1e-12)
    assert np.max(long.flow_accum[land]) < np.max(short.flow_accum[land])


def test_fixed_catchment_has_the_same_outlet_discharge_at_two_node_densities():
    # Split the same 10,000 km^2 rainy catchment into 4 or 64 accounting cells, the
    # fourfold linear-density change from density 1 to 4 (issue #355 step 5). Routing
    # area-weighted runoff must give the same outlet discharge and river forcing.
    catchment_area_m2 = 1e10
    runoff_m = breaching.runoff_mm(np.array([800.0]), breaching.open_water_evaporation_mm(np.array([15.0])))[0] / 1000.0
    outflows = []
    for density in (1, 4):
        source_count = 4 * density**2
        outlet, ocean = source_count, source_count + 1
        elevation = np.concatenate((np.full(source_count, 2.0), [1.0, 0.0]))
        is_ocean = np.arange(source_count + 2) == ocean
        target = np.concatenate((np.full(source_count, outlet), [ocean, -1]))
        source = np.concatenate((np.full(source_count, runoff_m * catchment_area_m2 / source_count), [0.0, 0.0]))
        through, settled = hydrology.route_downstream(elevation, is_ocean, target, source)
        assert settled[ocean] == pytest.approx(runoff_m * catchment_area_m2)
        outflows.append(through[outlet] / hydrology.SECONDS_PER_YEAR)

    assert outflows[0] == pytest.approx(outflows[1], rel=1e-12)
    assert outflows[0] == pytest.approx(runoff_m * catchment_area_m2 / hydrology.SECONDS_PER_YEAR)
    assert hydrology.discharge_m3_s_to_flow_accum_mm_yr(np.array(outflows))[0] == pytest.approx(
        hydrology.discharge_m3_s_to_flow_accum_mm_yr(np.array(outflows))[1]
    )
    assert hydrology.compute_river_speed(np.array([0.01, 0.01]), np.array(outflows))[0] == pytest.approx(
        hydrology.compute_river_speed(np.array([0.01, 0.01]), np.array(outflows))[1]
    )
    river = erosion._river_erosion_amount(
        np.array(outflows), np.zeros(2), np.array([0.01, 0.01]), np.zeros(2), 1.0, 1.0,
    )
    width = erosion._channel_width_growth(np.array(outflows), 1.0)
    assert river[0] > 0.0 and river[0] == pytest.approx(river[1])
    assert width[0] > 0.0 and width[0] == pytest.approx(width[1])


def test_compute_flow_direction_prefers_an_existing_channel_over_the_literal_steepest_neighbor():
    # Node 0's two downhill candidates: node 1 (elevation 20, steeper, no channel) and node 2
    # (elevation 25, less steep, but already has a real established channel).
    elevation = np.array([30.0, 20.0, 25.0])
    is_ocean = np.zeros(3, dtype=bool)
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])
    filled = np.array([30.0, 20.0, 25.0])
    spill = np.full(3, -1)

    not_frozen = np.zeros(3, dtype=bool)
    steepest, _ = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, np.zeros(3), not_frozen)
    assert steepest[0] == 1  # no channel anywhere -- plain steepest descent

    channel_depth = np.array([0.0, 0.0, 50.0])  # node 2's channel is well-established
    channelized, _ = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, channel_depth, not_frozen)
    assert channelized[0] == 2  # prefers the channelized neighbor even though it's less steep


def test_compute_flow_direction_ignores_a_barely_established_channel():
    # Same setup, but node 2's channel is below CHANNEL_PREFERENCE_THRESHOLD_M -- too new/
    # shallow to count as "an existing channel," so plain steepest descent still wins.
    elevation = np.array([30.0, 20.0, 25.0])
    is_ocean = np.zeros(3, dtype=bool)
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])
    filled = np.array([30.0, 20.0, 25.0])
    spill = np.full(3, -1)
    channel_depth = np.array([0.0, 0.0, hydrology.CHANNEL_PREFERENCE_THRESHOLD_M - 1.0])

    result, _ = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, channel_depth, np.zeros(3, dtype=bool))
    assert result[0] == 1


def _node_areas(world):
    return collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))


def test_seasonal_positive_degrees_matches_a_sampled_annual_cycle():
    mean = np.array([-20.0, -12.0, -5.0, 0.0, 3.0, 12.0, 25.0])
    amplitude = np.array([10.0, 15.0, 15.0, 8.0, 0.0, 10.0, 5.0])
    t = np.linspace(0.0, 2.0 * np.pi, 200_001)[:-1]
    sampled = np.mean(np.clip(mean[:, None] + amplitude[:, None] * np.cos(t)[None, :], 0.0, None), axis=1)
    assert hydrology.seasonal_positive_degrees(mean, amplitude) == pytest.approx(sampled, abs=1e-4)


def test_seasonal_positive_degrees_limits():
    # No swing: degrees above freezing of the annual mean itself. A swing that never reaches
    # the threshold: zero. A cycle centred on it: amplitude / pi. Fully above: the mean.
    assert np.array_equal(hydrology.seasonal_positive_degrees(np.array([-3.0, 0.0, 4.0]), np.zeros(3)), [0.0, 0.0, 4.0])
    assert hydrology.seasonal_positive_degrees(np.array([-15.0]), np.array([10.0]))[0] == 0.0
    assert hydrology.seasonal_positive_degrees(np.array([0.0]), np.array([9.0]))[0] == pytest.approx(9.0 / np.pi)
    assert hydrology.seasonal_positive_degrees(np.array([12.0]), np.array([10.0]))[0] == pytest.approx(12.0)
    # Increasing in both the mean and (once the swing reaches the threshold) the amplitude.
    means = np.linspace(-20.0, 5.0, 101)
    assert np.all(np.diff(hydrology.seasonal_positive_degrees(means, np.full_like(means, 12.0))[33:]) > 0.0)
    amps = np.linspace(0.0, 30.0, 101)
    assert np.all(np.diff(hydrology.seasonal_positive_degrees(np.full_like(amps, -6.0), amps)[21:]) > 0.0)


def test_compute_hydrology_threads_seasonal_amplitude_into_glacier_melt():
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    glacier = np.full(12, 200.0)
    precipitation = np.full(12, 100.0)
    temperature = np.full(12, -6.0)

    def run(amplitude):
        world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, np.zeros(12), glacier_depth=glacier)])
        amp = None if amplitude is None else np.full(12, amplitude)
        fields = hydrology.compute_hydrology(world, precipitation, temperature, years=1_000_000, seasonal_amplitude_at_nodes=amp)
        return fields.glacier_depth, fields.ice.ledger

    no_swing, still = run(None)
    assert np.array_equal(no_swing, run(0.0)[0])
    assert still.surface_melt_m3 == 0.0
    # A -6C mean with an 18C swing has summers well above freezing: the ice melts.
    _, melting = run(18.0)
    assert melting.surface_melt_m3 > 0.0
    assert melting.closing_m3 < still.closing_m3


def test_compute_hydrology_freezes_an_existing_lake_when_cold():
    # A single evenly-spaced, tightly-packed 16-point chain -- dense enough that every
    # point's FLOW_NEIGHBOR_COUNT=8 nearest neighbors really are just its ~4 closest
    # chain-neighbors on each side (keeping the whole chain one connected graph end to end,
    # so basin-spill can still find a real path to the ocean at the far end), with a genuine
    # local dip at index 7 -- lower than every one of its own 8 nearest neighbors, so it's a
    # real sink, not just a low point somewhere along a monotonic slope.
    d = 0.002
    theta = d * np.arange(16)
    elevation = np.array(
        [500, 450, 400, 350, 300, 250, 200, 50, 300, 250, 200, 150, 100, 50, 10, -50], dtype=float
    )
    lake_depth = np.zeros(16)
    lake_depth[7] = 30.0
    plate = _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
    world = World(seed=0, plates=[plate])

    precipitation = np.full(16, 500.0)
    cold_temperature = np.full(16, -20.0)  # everywhere well below freezing

    fields = hydrology.compute_hydrology(world, precipitation, cold_temperature, years=1_000_000)
    assert np.isfinite(fields.filled_elevation[7])  # confirms the sink really is graph-connected to the ocean

    assert fields.lake_depth[7] == 0.0  # the lake froze, not still sitting there as water
    # Its water joined the ice instead (issue #353: into its ice body, which may reshape it), at
    # the glacier model's water-to-ice rate, as frozen river water does -- not 1:1, which let
    # real-unit lake water outgrow all snowfall ~35x (#503).
    assert fields.ice.ledger.frozen_lake_m3 == pytest.approx(hydrology.WATER_TO_ICE * 30.0 * _node_areas(world)[7])
    assert fields.ice.ledger.closing_m3 >= fields.ice.ledger.frozen_lake_m3


def test_compute_hydrology_freezes_a_lake_below_freezing_even_above_glacier_accumulation_temp():
    # Same fixture as the -20C test above, but at a mild -3C, just below FREEZE_POINT_C (0,
    # the real freezing point). A lake here must still freeze solid: freezing is gated on the
    # annual mean, independent of the seasonal cycle glacier melt uses.
    d = 0.002
    theta = d * np.arange(16)
    elevation = np.array(
        [500, 450, 400, 350, 300, 250, 200, 50, 300, 250, 200, 150, 100, 50, 10, -50], dtype=float
    )
    lake_depth = np.zeros(16)
    lake_depth[7] = 30.0
    plate = _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
    world = World(seed=0, plates=[plate])

    precipitation = np.full(16, 500.0)
    mildly_cold_temperature = np.full(16, -3.0)  # just below FREEZE_POINT_C

    fields = hydrology.compute_hydrology(world, precipitation, mildly_cold_temperature, years=1_000_000)
    assert fields.lake_depth[7] == 0.0  # frozen solid despite being only mildly cold
    assert fields.ice.ledger.frozen_lake_m3 > 0.0


def test_compute_hydrology_freezing_river_blocks_flow_and_banks_ice_at_the_freeze_point():
    # Same 12-node monotonic-descent chain as the end-to-end test, but just cold enough (-3C)
    # to freeze. No liquid water should
    # reach the ocean at all this step -- it's blocked and frozen in place instead.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    precipitation_at_nodes = np.full(12, 800.0)
    temperature_at_nodes = np.full(12, -3.0)
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes, temperature_at_nodes, years=1_000_000)

    land = ~fields.is_ocean
    assert np.all(fields.flow_target[land] == -1)  # every land node frozen -- nothing routes onward
    assert np.all(fields.water_deposited[fields.is_ocean] == 0.0)  # so nothing reaches the ocean this step
    # The snow stays on the ground as ice instead (issue #353: gathered into its ice body).
    ledger = fields.ice.ledger
    assert ledger.accumulation_m3 > 0.0
    assert float(np.dot(fields.glacier_depth, _node_areas(world))) == pytest.approx(ledger.closing_m3 + ledger.frozen_river_ice_m3)


def test_compute_hydrology_snows_on_an_ice_cap_above_freezing_but_lets_its_water_drain():
    # Same 12-node descending chain as the freezing-river test, at a sampled temperature of
    # +2C -- above freezing. A bare chain routes rain to the sea. An otherwise-identical chain
    # that already carries a real ice cap (200 m of glacier) still takes its precipitation as
    # snow, but liquid water on it isn't frozen in place: it drains downhill, and none settles
    # as river ice. Freezing it there let a river rebuild a warm ice plug every step (#503).
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)  # last two nodes are below sea level -> ocean
    warm = np.full(12, 2.0)
    precip = np.full(12, 800.0)

    bare = _flow_line_plate(0, theta, elevation)
    bare_fields = hydrology.compute_hydrology(World(seed=0, plates=[bare]), precip, warm, years=1_000_000)
    bare_land = ~bare_fields.is_ocean
    assert np.any(bare_fields.flow_target[bare_land] >= 0)  # control: warm -> rivers route downhill
    assert bare_fields.water_deposited[bare_fields.is_ocean].sum() > 0.0  # control: rain reaches the sea

    iced = _flow_line_plate_with_lake(0, theta, elevation, np.zeros(12), glacier_depth=np.full(12, 200.0))
    iced_fields = hydrology.compute_hydrology(World(seed=0, plates=[iced]), precip, warm, years=1_000_000)
    iced_land = ~iced_fields.is_ocean
    assert np.all(iced_fields.flow_target[iced_land] >= 0)  # water on the warm cap drains downhill
    assert iced_fields.ice.ledger.accumulation_m3 > 0.0  # the precipitation still lands as snow
    assert iced_fields.ice.ledger.frozen_river_ice_m3 == 0.0  # and no water freezes where it settles


def test_group_rivers_groups_confluence_and_finds_max_flow_mouth():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    assert len(rivers) == 4

    network_a = next(r for r in rivers if len(r.member_idx) == 4)
    assert sorted(network_a.member_idx.tolist()) == [0, 1, 2, 3]
    assert network_a.mouth_idx == 3  # highest fallback discharge (converted from flow_accum)
    assert network_a.flow_rate == pytest.approx(
        float(hydrology.flow_accum_mm_yr_to_discharge_m3_s(40.0))
    )
    # Two headwaters (0 and 1) merging into one mainstem -- one tributary joins it.
    assert network_a.num_tributaries == 1


def test_river_selection_stays_connected_through_positive_downstream_discharge():
    # Transit evaporation can reduce discharge below the percentile cutoff along a real river.
    # Once selected upstream, its positive-flow path remains connected to the mouth.
    selected = hydrology._extend_river_mask_downstream(
        np.array([True, False, False, False, False]),
        np.array([1, 2, 3, 4, -1]),
        np.array([10.0, 0.5, 0.01, 0.0, 0.0]),
        np.zeros(5, dtype=bool),
        np.zeros(5, dtype=bool),
    )
    assert selected.tolist() == [True, True, True, False, False]

    # Multiple selected headwaters may share a mainstem; walking it once must still retain the
    # whole shared path without revisiting it for every tributary.
    confluence = hydrology._extend_river_mask_downstream(
        np.array([True, True, False, False, False]),
        np.array([2, 2, 3, 4, -1]),
        np.array([10.0, 5.0, 0.5, 0.01, 0.0]),
        np.zeros(5, dtype=bool),
        np.zeros(5, dtype=bool),
    )
    assert confluence.tolist() == [True, True, True, True, False]


def test_river_selection_stops_downstream_extension_at_lakes_and_oceans():
    flow_target = np.array([1, 2, 3, 4, -1])
    discharge = np.ones(5)
    is_ocean = np.zeros(5, dtype=bool)
    is_lake = np.zeros(5, dtype=bool)
    is_lake[2] = True
    lake_path = hydrology._extend_river_mask_downstream(
        np.array([True, False, False, False, False]), flow_target, discharge, is_ocean, is_lake
    )
    assert lake_path.tolist() == [True, True, False, False, False]

    ocean = np.zeros(5, dtype=bool)
    ocean[2] = True
    ocean_path = hydrology._extend_river_mask_downstream(
        np.array([True, False, False, False, False]), flow_target, discharge, ocean, np.zeros(5, dtype=bool)
    )
    assert ocean_path.tolist() == [True, True, False, False, False]


def test_group_rivers_classifies_mouth_type_ocean_lake_and_other():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)

    network_a = next(r for r in rivers if len(r.member_idx) == 4)
    assert network_a.mouth_type == "ocean"  # mouth's flow_target (4) is an ocean node

    network_b = next(r for r in rivers if 5 in r.member_idx.tolist())
    assert network_b.mouth_idx == 5
    # lake_depth is read off *node 7* (what node 5 flows into), not node 5's own -- node 5
    # itself, like every real river node, is never itself is_lake.
    assert network_b.mouth_type == "lake"
    assert network_b.num_tributaries == 0  # a single unbranched node has no tributaries

    network_c = next(r for r in rivers if 8 in r.member_idx.tolist())
    assert network_c.mouth_idx == 8
    assert network_c.mouth_type == "other"  # flow_target -1: a still-dry, unresolved sink

    network_d = next(r for r in rivers if 9 in r.member_idx.tolist())
    assert network_d.mouth_idx == 9
    # flow_target[9] (node 10) is itself still dry -- only the chain's *next* hop (node 7)
    # actually carries standing water. Stopping at node 10 alone would misclassify this as
    # "other", the exact bug a real save reproduced (see group_rivers's own docstring).
    assert network_d.mouth_type == "lake"


def test_group_rivers_excludes_non_river_nodes():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    all_members = {i for r in rivers for i in r.member_idx.tolist()}
    assert 4 not in all_members  # ocean node, never is_river
    assert 6 not in all_members  # land node deliberately not classified as a river
    assert 7 not in all_members  # the lake node itself -- is_river excludes it too
    assert 10 not in all_members  # dry basin floor between node 9 and the lake, not is_river


def test_group_rivers_on_no_rivers_returns_empty_list():
    fields = _river_inspector_fields()
    fields.is_river = np.zeros(len(fields.is_river), dtype=bool)
    assert hydrology.group_rivers(fields) == []


def test_river_at_finds_the_network_nearest_a_query_point():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    network_a_id = next(i for i, r in enumerate(rivers) if len(r.member_idx) == 4)
    network_b_id = next(i for i, r in enumerate(rivers) if 5 in r.member_idx.tolist())

    near_a = hydrology.river_at(fields, rivers, fields.points[0])
    assert near_a == network_a_id

    near_b = hydrology.river_at(fields, rivers, fields.points[5])
    assert near_b == network_b_id


def test_river_at_returns_none_when_there_are_no_rivers():
    fields = _river_inspector_fields()
    assert hydrology.river_at(fields, [], fields.points[0]) is None


def test_group_rivers_on_a_real_monotonic_chain_is_one_unbranched_network():
    # Reuses the same 12-node monotonic-descent fixture as the compute_hydrology end-to-end
    # test above -- a single chain has no confluence, so grouping it should yield exactly one
    # network with zero tributaries, mouthing into the ocean.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    fields = hydrology.compute_hydrology(world, np.full(12, 800.0), np.full(12, 15.0), years=1_000_000)
    rivers = hydrology.group_rivers(fields)

    assert len(rivers) == 1
    assert rivers[0].num_tributaries == 0
    assert rivers[0].mouth_type == "ocean"
    assert fields.is_ocean[fields.flow_target[rivers[0].mouth_idx]]


def test_compute_river_speed_increases_with_slope_and_flow():
    base = hydrology.compute_river_speed(np.array([0.01]), np.array([1.0]))[0]
    steeper = hydrology.compute_river_speed(np.array([0.1]), np.array([1.0]))[0]
    more_flow = hydrology.compute_river_speed(np.array([0.01]), np.array([100.0]))[0]
    assert steeper > base
    assert more_flow > base


# -- Issue #133 phase 2: node_cloud_resample_mode in sample_is_ocean/sample_is_sea -----------


def _grid_xyz(height, width):
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_grid = np.repeat(lon_deg[None, :], height, axis=0)
    return geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))


def test_sample_is_ocean_and_is_sea_healpix_mode_agree_closely_with_kdtree():
    world = generate_world(30, num_plates=8)
    step_world(world, 1_000_000)  # populates world.hydrology_cache
    world_xyz = _grid_xyz(31, 61)
    fallback = np.zeros((31, 61), dtype=bool)

    world.node_cloud_resample_mode = "kdtree"
    ocean_kd = hydrology.sample_is_ocean(world, world_xyz, fallback)
    sea_kd = hydrology.sample_is_sea(world, world_xyz, fallback)

    world.node_cloud_resample_mode = "healpix"
    ocean_hp = hydrology.sample_is_ocean(world, world_xyz, fallback)
    sea_hp = hydrology.sample_is_sea(world, world_xyz, fallback)

    assert (ocean_kd == ocean_hp).mean() > 0.95
    assert (sea_kd == sea_hp).mean() > 0.95


def test_nearest_hydro_node_idx_cache_does_not_cross_contaminate_between_modes():
    """The single-entry cache (`_LAST_NEAREST_HYDRO_NODE`) must key on `resample_mode` too --
    otherwise a "kdtree" call followed immediately by a "healpix" call against the identical
    `(hydro, query_xyz)` would incorrectly serve the first mode's stale result."""
    world = generate_world(31, num_plates=8)
    step_world(world, 1_000_000)
    hydro = world.hydrology_cache
    world_xyz = _grid_xyz(21, 41)

    idx_kd = hydrology._nearest_hydro_node_idx(hydro, world_xyz, "kdtree")
    hydrology._nearest_hydro_node_idx(hydro, world_xyz, "healpix")  # overwrites the single-entry cache
    idx_kd_again = hydrology._nearest_hydro_node_idx(hydro, world_xyz, "kdtree")

    assert np.array_equal(idx_kd, idx_kd_again)  # re-fetching "kdtree" wasn't clobbered by the "healpix" call
    tree = hydrology._ocean_node_tree(hydro)
    _, expected_kd = tree.query(world_xyz.reshape(-1, 3), workers=1)
    assert np.array_equal(idx_kd, expected_kd)


def test_sample_is_ocean_fallback_before_first_hydrology_pass_is_mode_independent():
    world = generate_world(32, num_plates=8)  # no step_world -- no hydrology_cache yet
    world_xyz = _grid_xyz(11, 21)
    fallback = np.zeros((11, 21), dtype=bool)
    fallback[0, 0] = True

    for mode in ("kdtree", "healpix"):
        world.node_cloud_resample_mode = mode
        result = hydrology.sample_is_ocean(world, world_xyz, fallback)
        assert np.array_equal(result, fallback)


def test_route_downstream_conserves_across_uphill_edges_and_cycles():
    # Issue #275: spill/ice routing can point a node at a *higher* neighbour. An elevation
    # sweep used to add that flux to a node it had already visited, so it was never passed on.
    elevation = np.array([100.0, 50.0, 80.0, 10.0, 0.0, 30.0, 31.0])
    is_ocean = np.array([False, False, False, False, True, False, False])
    # 0 -> 1 -> 2 (uphill spill) -> 3 -> ocean 4; 5 <-> 6 is a cycle.
    flow_target = np.array([1, 2, 3, 4, -1, 6, 5])
    source = np.array([1.0, 2.0, 3.0, 4.0, 0.0, 5.0, 6.0])
    retain = np.array([0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0])
    flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source, retain_fraction=retain)
    assert np.isclose(deposited.sum(), source.sum())
    assert np.isclose(flux[2], 6.0)  # nodes 0 and 1 both reached the spill node
    assert np.isclose(deposited[3], 5.0) and np.isclose(deposited[4], 5.0)
    assert np.isclose(deposited[5] + deposited[6], 11.0)  # a cycle keeps its own load


def test_route_downstream_conserves_on_an_all_downhill_graph():
    rng = np.random.default_rng(0)
    n = 200
    elevation = rng.uniform(0.0, 1000.0, n)
    is_ocean = elevation < 100.0
    order = np.argsort(elevation)
    flow_target = np.full(n, -1)
    for rank, i in enumerate(order[1:], start=1):
        flow_target[i] = order[rng.integers(0, rank)]  # always a strictly lower node
    source = rng.uniform(0.0, 1.0, n)
    flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source)
    assert np.isclose(deposited.sum(), source[~is_ocean].sum())


def test_route_downstream_applies_loss_and_retention_on_cycle_nodes():
    # Review on #278: a node on a flow cycle must still evaporate/retain like any other node
    # before keeping its remainder -- otherwise river water skips evaporation there.
    elevation = np.array([100.0, 30.0, 31.0, 0.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 1, -1])  # 0 feeds the 1 <-> 2 cycle
    source = np.array([4.0, 2.0, 2.0, 0.0])
    loss = np.array([0.0, 0.5, 0.25, 0.0])
    retain = np.array([0.0, 0.5, 0.0, 0.0])
    flux, deposited = hydrology.route_downstream(
        elevation, is_ocean, flow_target, source, retain_fraction=retain, loss_fraction=loss
    )
    # Node 1 receives 2 + 4 = 6, loses half (3), retains half the rest; node 2 loses a quarter.
    assert np.isclose(deposited[1], 3.0) and np.isclose(flux[1], 1.5)
    assert np.isclose(deposited[2], 1.5) and np.isclose(flux[2], 1.5)
    assert np.isclose(deposited.sum(), 4.5)  # 8 in, 3.5 evaporated

    # Without loss, the cycle still conserves everything it received.
    _, kept = hydrology.route_downstream(elevation, is_ocean, flow_target, source, retain_fraction=retain)
    assert np.isclose(kept.sum(), source.sum())


def test_physical_water_accounts_for_transit_loss_lake_retention_and_ocean_arrival():
    # Three land nodes drain through a lake (node 1) to the ocean. The yearly source volume
    # must end in atmospheric loss, lake retention or ocean arrival exactly once (#355).
    elevation = np.array([30.0, 20.0, 10.0, 0.0])
    is_ocean = np.array([False, False, False, True])
    target = np.array([1, 2, 3, -1])
    source_m3_per_yr = np.array([10.0, 5.0, 2.0, 0.0])
    loss = np.array([0.1, 0.2, 0.0, 0.0])
    lake_retention = np.array([0.0, 0.25, 0.0, 0.0])
    through, settled = hydrology.route_downstream(
        elevation, is_ocean, target, source_m3_per_yr,
        retain_fraction=lake_retention, loss_fraction=loss,
    )
    evaporated = through[0] * loss[0] / (1.0 - loss[0]) + (
        through[1] + settled[1]
    ) * loss[1] / (1.0 - loss[1])
    assert settled[1] == pytest.approx(2.8)
    assert settled[3] == pytest.approx(10.4)
    assert source_m3_per_yr.sum() == pytest.approx(evaporated + settled.sum())


def test_ice_ledger_closes_across_steps_with_frozen_river_ice_a_step_late():
    # Cold, snowy uplands above a band whose summers melt snow that refreezes in place as
    # river ice: that ice is added after the solve and joins next step's opening volume.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])
    precipitation = np.full(12, 800.0)
    temperature = np.where(np.arange(12) < 4, -20.0, -3.0)
    amplitude = np.where(np.arange(12) < 4, 0.0, 10.0)
    areas = _node_areas(world)
    ledgers = []
    for _ in range(2):
        fields = hydrology.compute_hydrology(world, precipitation, temperature, years=1_000_000, seasonal_amplitude_at_nodes=amplitude)
        ledgers.append(fields.ice.ledger)
        assert float(np.dot(fields.glacier_depth, areas)) == pytest.approx(fields.ice.ledger.closing_m3 + fields.ice.ledger.frozen_river_ice_m3)
        plate.set_fields_on_plate(glacier_depth=fields.glacier_depth)
    first, second = ledgers
    assert first.frozen_river_ice_m3 > 0.0
    assert second.opening_m3 == pytest.approx(first.closing_m3 + first.frozen_river_ice_m3)


def test_ice_flow_target_follows_the_ice_surface_on_the_ice_and_the_ground_off_it():
    # A dome on a flat bed: on the ice, material moves down the ice surface (away from the
    # divide), even where the bed under it doesn't slope; off the ice, the ground's own target.
    elevation = np.zeros(5)
    depth = np.array([0.0, 500.0, 1000.0, 500.0, 0.0])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 3]])
    ground = np.array([-1, 0, -1, 4, -1])
    target = hydrology._ice_surface_flow_target(elevation + depth, neighbor_idx, depth, ground)
    assert target[1] == 0 and target[3] == 4
    assert target[2] in (1, 3)
    assert target[0] == -1 and target[4] == -1


def test_a_lake_under_ice_on_warm_ground_stays_water():
    # Issue #353: only cold ground (annual mean below freezing) freezes a lake. A lake beside or
    # under ice on warm ground is refilled by runoff in real units every step; freezing it into
    # the ice 1:1 fed the conserved ice volume far faster than snowfall.
    d = 0.002
    theta = d * np.arange(16)
    elevation = np.array([500, 450, 400, 350, 300, 250, 200, 50, 300, 250, 200, 150, 100, 50, 10, -50], dtype=float)
    lake_depth = np.zeros(16)
    lake_depth[7] = 30.0
    world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, lake_depth, glacier_depth=np.full(16, 200.0))])
    fields = hydrology.compute_hydrology(world, np.full(16, 500.0), np.full(16, 2.0), years=1_000_000)
    ledger = fields.ice.ledger
    assert ledger.frozen_lake_m3 == 0.0
    assert ledger.closing_m3 <= ledger.opening_m3 + ledger.accumulation_m3  # the lake didn't become ice

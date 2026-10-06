import numpy as np
import pytest

from app import breaching, lakes


def _chain(n):
    """A 1-D chain, each node linked to its two neighbours (ends padded with a duplicate)."""
    left = np.maximum(np.arange(n) - 1, 0)
    right = np.minimum(np.arange(n) + 1, n - 1)
    left[0] = 1
    right[-1] = n - 2
    return np.stack([left, right], axis=1)


def _ordinary_rate(n):
    return np.full(n, breaching.BREACH_REFERENCE_CARVE_M_PER_MYR)


def test_least_climb_sums_every_barrier_on_the_path():
    # Ocean at node 4. Water from node 0 climbs 10 m onto node 1, drops, then climbs 20 m onto
    # node 3. Its total rise is 30 m, not the 20 m of the single highest barrier.
    elevation = np.array([0.0, 10.0, -5.0, 15.0, -100.0])
    is_ocean = np.array([False, False, False, False, True])
    cost, next_hop = breaching.least_climb_to_ocean(elevation, is_ocean, _chain(5), _ordinary_rate(5))
    assert cost[0] == pytest.approx(30.0, abs=1e-3)
    assert cost[2] == pytest.approx(20.0, abs=1e-3)
    assert cost[3] == pytest.approx(0.0, abs=1e-3)
    assert cost[4] == 0.0
    assert next_hop.tolist() == [1, 2, 3, 4, -1]


def test_unreachable_nodes_cost_infinity():
    elevation = np.array([0.0, 5.0, -10.0, 3.0])
    is_ocean = np.array([False, False, True, False])
    neighbor_idx = np.array([[1], [0], [2], [3]])  # 0/1 and 3 never reach node 2
    cost, next_hop = breaching.least_climb_to_ocean(elevation, is_ocean, neighbor_idx, _ordinary_rate(4))
    assert np.isinf(cost[[0, 1, 3]]).all()
    assert (next_hop[[0, 1, 3]] == -1).all()


def test_a_cheap_pit_is_breached_and_drains_at_its_own_floor():
    # A 20 m deep pit (node 1) behind a rim (node 2). One 100 kyr step carves 30 m of ordinary
    # rock, so the path is notched down to the pit floor and the lake hierarchy spills at 0 m.
    # Node 3 is past the rim but still 10 m above the floor, so it's cut too.
    elevation = np.array([30.0, 0.0, 20.0, 10.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    neighbor_idx = _chain(5)
    result = breaching.breach_depressions(elevation, is_ocean, neighbor_idx, np.zeros(5), _ordinary_rate(5), years=100_000)

    assert result.breached_pits.tolist() == [1]
    assert result.closed_pits.tolist() == []
    assert result.notch_m[2] == pytest.approx(20.0)
    assert result.notch_m[3] == pytest.approx(10.0)
    assert result.notch_m[[0, 1, 4]].tolist() == [0.0, 0.0, 0.0]

    passes = breaching.interface_pass_elevation(result.passage_m, neighbor_idx)
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx, interface_pass_elevation=passes)
    filled, spill, _ = lakes.compute_spill_routing(forest, elevation, np.zeros(5))
    assert filled[1] == pytest.approx(0.0)
    assert spill[1] == 3  # node 2 drains into the pit's own catchment, so the boundary is 2->3


def test_a_costly_pit_stays_closed_and_fills_to_its_rim():
    # The same pit 50 m deep exceeds one step's 30 m budget. It isn't notched, so lakes.py
    # still sees a 50 m rim.
    elevation = np.array([60.0, 0.0, 50.0, 10.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    neighbor_idx = _chain(5)
    result = breaching.breach_depressions(elevation, is_ocean, neighbor_idx, np.zeros(5), _ordinary_rate(5), years=100_000)

    assert result.breached_pits.tolist() == []
    assert result.closed_pits.tolist() == [1]
    assert not result.notch_m.any()
    passes = breaching.interface_pass_elevation(result.passage_m, neighbor_idx)
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx, interface_pass_elevation=passes)
    filled, _, _ = lakes.compute_spill_routing(forest, elevation, np.zeros(5))
    assert filled[1] == pytest.approx(50.0)


def test_a_longer_step_carves_a_deeper_breach():
    elevation = np.array([60.0, 0.0, 50.0, 10.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    result = breaching.breach_depressions(elevation, is_ocean, _chain(5), np.zeros(5), _ordinary_rate(5), years=1_000_000)
    assert result.breached_pits.tolist() == [1]
    assert result.notch_m[2] == pytest.approx(50.0)


def test_rock_strength_decides_whether_the_same_barrier_is_cut():
    elevation = np.array([30.0, 0.0, 20.0, 10.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    craton = breaching.bedrock_carve_rate_m_per_myr(np.ones(5))
    ordinary = breaching.bedrock_carve_rate_m_per_myr(np.zeros(5))
    assert craton[0] == pytest.approx(breaching.BREACH_STRONG_CARVE_M_PER_MYR)
    assert ordinary[0] == pytest.approx(breaching.BREACH_REFERENCE_CARVE_M_PER_MYR)

    strong = breaching.breach_depressions(elevation, is_ocean, _chain(5), np.zeros(5), craton, years=100_000)
    # Deep loose cover over ordinary rock: the whole 20 m climb is through cover.
    weak = breaching.breach_depressions(
        elevation, is_ocean, _chain(5), np.zeros(5), ordinary, years=100_000, mobile_cover_m=np.full(5, 100.0)
    )
    assert strong.breached_pits.tolist() == []  # 20 m of craton costs 60 m of ordinary rock
    assert weak.breached_pits.tolist() == [1]


def test_an_established_channel_is_a_notch_through_a_higher_centre():
    # Node 1 is a lake floor. Its downstream neighbour, node 2, has a centre 80 m higher, above
    # any one-step budget. But node 2 carries a 200 m channel, which notches it down to its
    # lowest neighbour (-1 m). The lake drains through the channel without inventing a notch,
    # and doesn't wait for node 2's whole centre to be submerged.
    elevation = np.array([90.0, 0.0, 80.0, -1.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    channel = np.array([0.0, 0.0, 200.0, 0.0, 0.0])
    neighbor_idx = _chain(5)
    result = breaching.breach_depressions(elevation, is_ocean, neighbor_idx, channel, _ordinary_rate(5), years=100_000)

    assert result.channel_passage_m[2] == pytest.approx(-1.0)
    assert result.cost_m[1] == pytest.approx(0.0, abs=1e-3)
    assert not result.notch_m.any()
    passes = breaching.interface_pass_elevation(result.passage_m, neighbor_idx)
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx, interface_pass_elevation=passes)
    filled, spill, _ = lakes.compute_spill_routing(forest, elevation, np.zeros(5))
    assert filled[1] == pytest.approx(0.0)
    assert spill[1] == 2


def test_a_channel_never_notches_below_the_cells_lowest_neighbour():
    elevation = np.array([100.0, 50.0, 40.0])
    neighbor_idx = _chain(3)
    passage = breaching.channel_passage_elevation(elevation, np.array([0.0, 500.0, 0.0]), neighbor_idx)
    assert passage.tolist() == [100.0, 40.0, 40.0]


def test_two_pits_sharing_a_path_cut_to_the_lower_arriving_level():
    # Pits 1 (at 0 m) and 3 (at 5 m) both drain through rim node 4 (20 m) to the ocean. The
    # rim is cut to the lower pit's level, and pit 3's own barrier (node 2) needs no cut
    # because node 3 is downstream of it.
    elevation = np.array([40.0, 0.0, 12.0, 5.0, 20.0, -50.0])
    is_ocean = np.array([False] * 5 + [True])
    result = breaching.breach_depressions(elevation, is_ocean, _chain(6), np.zeros(6), _ordinary_rate(6), years=100_000)
    assert sorted(result.breached_pits.tolist()) == [1, 3]
    assert result.notch_m[4] == pytest.approx(20.0)
    assert result.notch_m[2] == pytest.approx(12.0)
    assert result.passage_m[4] == pytest.approx(0.0)


def test_breaching_is_independent_of_node_order():
    rng = np.random.default_rng(3)
    n = 60
    elevation = rng.uniform(0.0, 40.0, n)
    elevation[:5] = -100.0
    is_ocean = elevation < 0.0
    points = rng.normal(size=(n, 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    dist = np.linalg.norm(points[:, None] - points[None], axis=2)
    neighbor_idx = np.argsort(dist, axis=1)[:, 1:7]
    rate = _ordinary_rate(n)

    base = breaching.breach_depressions(elevation, is_ocean, neighbor_idx, np.zeros(n), rate, years=100_000)
    perm = rng.permutation(n)
    inverse = np.argsort(perm)
    permuted = breaching.breach_depressions(
        elevation[perm], is_ocean[perm], inverse[neighbor_idx[perm]], np.zeros(n), rate, years=100_000
    )
    np.testing.assert_allclose(permuted.cost_m, base.cost_m[perm], atol=1e-4)
    np.testing.assert_allclose(permuted.notch_m, base.notch_m[perm], atol=1e-4)


def test_runoff_follows_the_budyko_limits():
    # Cold and wet: almost all rain runs off. Hot and dry: almost none does.
    wet = breaching.runoff_mm(np.array([2000.0]), np.array([100.0]))[0]
    dry = breaching.runoff_mm(np.array([200.0]), np.array([2000.0]))[0]
    assert wet > 0.9 * 2000.0
    assert dry < 0.01 * 200.0
    assert breaching.runoff_mm(np.array([0.0]), np.array([1000.0]))[0] == 0.0


def _pit_water(n, precipitation, temperature, area=1.0e10):
    return breaching.WaterBalance(np.full(n, precipitation), np.full(n, temperature), np.full(n, area))


def test_a_cheap_pit_with_no_water_stays_closed():
    elevation = np.array([30.0, 0.0, 20.0, 10.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    result = breaching.breach_depressions(
        elevation, is_ocean, _chain(5), np.zeros(5), _ordinary_rate(5), years=100_000,
        water=_pit_water(5, precipitation=0.0, temperature=20.0),
    )
    assert result.breached_pits.tolist() == []
    assert result.endorheic_pits.tolist() == [1]
    assert not result.notch_m.any()


def test_a_cheap_pit_breaches_only_when_its_catchment_outpaces_evaporation():
    # Pit 3 drains a long slope (nodes 0-2). Wet and cool, the slope's runoff swamps a lake at
    # the rim and the pit is breached. Hot and dry, the lake evaporates it all and stays closed.
    elevation = np.array([60.0, 40.0, 20.0, 0.0, 15.0, -50.0])
    is_ocean = np.array([False] * 5 + [True])
    args = (elevation, is_ocean, _chain(6), np.zeros(6), _ordinary_rate(6))
    wet = breaching.breach_depressions(*args, years=100_000, water=_pit_water(6, 1500.0, 5.0))
    dry = breaching.breach_depressions(*args, years=100_000, water=_pit_water(6, 250.0, 25.0))
    assert wet.breached_pits.tolist() == [3]
    assert wet.notch_m[4] == pytest.approx(15.0)
    assert dry.breached_pits.tolist() == []
    assert dry.endorheic_pits.tolist() == [3]


def test_endorheic_demand_factor_keeps_more_basins_closed(monkeypatch):
    elevation = np.array([60.0, 40.0, 20.0, 0.0, 15.0, -50.0])
    is_ocean = np.array([False] * 5 + [True])
    args = (elevation, is_ocean, _chain(6), np.zeros(6), _ordinary_rate(6))
    water = _pit_water(6, 1500.0, 5.0)
    assert breaching.breach_depressions(*args, years=100_000, water=water).breached_pits.tolist() == [3]
    monkeypatch.setattr(breaching, "ENDORHEIC_DEMAND_FACTOR", 100.0)
    assert breaching.breach_depressions(*args, years=100_000, water=water).breached_pits.tolist() == []


def test_a_pit_whose_path_crosses_rock_too_strong_to_drain_stays_closed():
    # Review repro (#301): a silt pit and silt rim, then a craton node 80 m above the floor.
    # The cost (27 m of reference rock) fits the 30 m budget, but the craton can only be cut
    # 10 m this step, so the path can't reach the pit floor. The pit stays closed and no
    # notches are banked.
    elevation = np.array([100.0, 0.0, 90.0, 80.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 3]])
    rate = breaching.bedrock_carve_rate_m_per_myr(np.array([0.0, 0.0, 0.0, 1.0, 0.0]))
    cover = np.array([100.0, 100.0, 100.0, 0.0, 0.0])
    result = breaching.breach_depressions(
        elevation, is_ocean, neighbor_idx, np.zeros(5), rate, years=100_000, mobile_cover_m=cover
    )

    assert result.cost_m[1] <= breaching.BREACH_REFERENCE_CARVE_M_PER_MYR * 0.1
    assert result.breached_pits.tolist() == []
    assert result.closed_pits.tolist() == [1]
    assert not result.notch_m.any()


def test_every_breached_pit_drains_to_its_floor():
    # Mixed rock where every cut fits: the pit's lake hierarchy spills at the floor.
    elevation = np.array([100.0, 0.0, 25.0, 8.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 3]])
    rate = breaching.bedrock_carve_rate_m_per_myr(np.array([0.0, 0.0, 0.0, 1.0, 0.0]))
    result = breaching.breach_depressions(elevation, is_ocean, neighbor_idx, np.zeros(5), rate, years=100_000)
    assert result.breached_pits.tolist() == [1]
    passes = breaching.interface_pass_elevation(result.passage_m, neighbor_idx)
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx, interface_pass_elevation=passes)
    filled, _, _ = lakes.compute_spill_routing(forest, elevation, np.zeros(5))
    assert filled[1] == pytest.approx(0.0)


def test_a_persisted_notch_lowers_the_passage_and_new_cuts_are_only_the_increment():
    # A 50 m pit. Last step already notched the rim 30 m, so this step only needs 20 m more.
    elevation = np.array([60.0, 0.0, 50.0, -50.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = _chain(4)
    prior = np.array([0.0, 0.0, 30.0, 0.0])
    result = breaching.breach_depressions(
        elevation, is_ocean, neighbor_idx, np.zeros(4), _ordinary_rate(4), years=100_000, prior_notch_m=prior
    )
    assert result.channel_passage_m[2] == pytest.approx(20.0)
    assert result.breached_pits.tolist() == [1]
    assert result.notch_m[2] == pytest.approx(20.0)


def test_the_water_gate_uses_the_rim_height_not_the_rock_weighted_cost():
    # Review repro (#301, second pass): a long wet slope drains into pit 6 behind a 90 m silt
    # rim. Silt's low resistance makes that rim cost only 27 m of reference rock, but a lake
    # at the real 90 m rim evaporates more than this catchment supplies, so the basin stays
    # closed.
    z = np.array([85.0, 70.0, 55.0, 40.0, 25.0, 10.0, 0.0, 90.0, -50.0])
    n = len(z)
    ocean = np.zeros(n, dtype=bool)
    ocean[-1] = True
    rate = breaching.bedrock_carve_rate_m_per_myr(np.zeros(n))
    cover = np.where(np.arange(n) >= 3, 100.0, 0.0)
    water = breaching.WaterBalance(np.full(n, 1600.0), np.full(n, 12.0), np.full(n, 1e9))
    result = breaching.breach_depressions(
        z, ocean, _chain(n), np.zeros(n), rate, 1_000_000, water=water, mobile_cover_m=cover
    )
    assert result.cost_m[6] == pytest.approx(27.0, abs=1e-3)
    assert result.breached_pits.tolist() == []
    assert result.endorheic_pits.tolist() == [6]
    assert not result.notch_m.any()


def test_the_gate_routes_runoff_with_the_same_river_losses_as_route_downstream():
    from app import hydrology

    elevation = np.array([60.0, 40.0, 20.0, 0.0, 15.0, -50.0])
    is_ocean = np.array([False] * 5 + [True])
    neighbor_idx = _chain(6)
    loss = np.array([0.1, 0.2, 0.3, 0.05, 0.0, 0.0])
    sink, carried = breaching._routed_to_sinks(elevation, neighbor_idx, is_ocean, 1.0 - loss)
    assert sink[:4].tolist() == [3, 3, 3, 3]
    amount = np.array([1.0, 2.0, 3.0, 4.0, 0.0, 0.0])
    flow_target = np.array([1, 2, 3, -1, 5, -1])
    _, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, amount, loss_fraction=loss)
    assert np.bincount(sink, weights=amount * carried, minlength=6)[3] == pytest.approx(deposited[3])


def test_only_the_share_of_a_cut_in_cover_carves_at_the_weak_rate():
    # Review (#301, fourth pass): 15 m of sediment over ordinary rock must not make a 60 m
    # climb weak all the way down.
    weak = breaching.BREACH_WEAK_CARVE_M_PER_MYR
    rock = breaching.BREACH_REFERENCE_CARVE_M_PER_MYR
    elevation = np.array([100.0, 0.0, 60.0, -50.0])
    is_ocean = np.array([False, False, False, True])
    cover = np.array([0.0, 0.0, 15.0, 0.0])
    cost, _ = breaching.least_climb_to_ocean(elevation, is_ocean, _chain(4), np.full(4, rock), cover)
    assert cost[1] == pytest.approx((15.0 / weak + 45.0 / rock) * rock, abs=1e-3)
    # In 100 kyr the river gets through the cover in 1.5 kyr, then cuts bedrock for the rest.
    limit = breaching._carve_limit_m(np.array([0.0, 15.0, 500.0]), np.full(3, rock), 0.1)
    np.testing.assert_allclose(limit, [rock * 0.1, 15.0 + (0.1 - 15.0 / weak) * rock, weak * 0.1])


def test_cover_above_an_existing_channel_floor_does_not_soften_the_cut():
    # Review repro (#301, fifth pass): the rim's 280 m channel leaves its passage 20 m above
    # the pit, far below its 15 m of surface cover. The remaining cut is all bedrock, so the
    # cover changes nothing: 20 m of ordinary rock exceeds a 50 kyr step's 15 m budget.
    z = np.array([400.0, 0.0, 300.0, -50.0])
    ocean = np.array([False, False, False, True])
    nb = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    channel = np.array([0.0, 0.0, 280.0, 0.0])
    rate = breaching.bedrock_carve_rate_m_per_myr(np.zeros(4))
    results = [
        breaching.breach_depressions(z, ocean, nb, channel, rate, 50_000, mobile_cover_m=np.array([0.0, 0.0, cover, 0.0]))
        for cover in (0.0, 15.0)
    ]
    for result in results:
        assert result.cost_m[1] == pytest.approx(20.0, abs=1e-3)
        assert result.breached_pits.tolist() == []
        assert not result.notch_m.any()

    # Cover thicker than the channel is deep still has some left below the passage: 300 m of
    # cover over a 280 m channel leaves 20 m, so the same cut is all cover and fits.
    thick = breaching.breach_depressions(z, ocean, nb, channel, rate, 50_000, mobile_cover_m=np.array([0.0, 0.0, 300.0, 0.0]))
    weak, rock = breaching.BREACH_WEAK_CARVE_M_PER_MYR, breaching.BREACH_REFERENCE_CARVE_M_PER_MYR
    assert thick.cost_m[1] == pytest.approx(20.0 / weak * rock, abs=1e-3)
    assert thick.breached_pits.tolist() == [1]

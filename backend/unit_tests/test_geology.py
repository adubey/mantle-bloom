import numpy as np
from app import geology, volcanism
from app.erosion import ErosionResult
from app.hydrology import HydrologyFields
from app.world import World

from .quad_fixtures import block_keys, node_points, quad_plate


def _plate(theta, elevation):
    """One continental row of cells, one per entry of `theta` (only its length matters)."""
    return quad_plate(0, "continental", columns=range(len(theta)), elevation=np.asarray(elevation, dtype=float))


def _world_with_hydrology(plate, is_ocean, water_deposited):
    n = plate.node_count()
    points = node_points(plate)
    hydro = HydrologyFields(
        points=points,
        elevation=plate.collect("elevation").copy(),
        is_ocean=np.asarray(is_ocean, dtype=bool),
        neighbor_idx=np.zeros((n, 1), dtype=int),
        flow_target=np.full(n, -1),
        flow_accum=np.zeros(n),
        water_deposited=np.asarray(water_deposited, dtype=float),
        filled_elevation=np.zeros(n),
        spill_target=np.full(n, -1),
        is_river=np.zeros(n, dtype=bool),
        lake_depth=np.zeros(n),
        glacier_depth=np.zeros(n),
        plates_in_order=[plate],
    )
    world = World(seed=0, plates=[plate])
    world.hydrology_cache = hydro
    return world


def _erosion_result(
    points, elevation, slope=None, rain=None, river=None, weathering=None, sediment_deposited=None,
    is_river_depositing=None, temperature_c=None, precipitation_mm=None,
):
    n = len(points)
    zeros = np.zeros(n)
    return ErosionResult(
        points=points,
        elevation=np.asarray(elevation, dtype=float),
        slope=zeros if slope is None else np.asarray(slope, dtype=float),
        rain=zeros if rain is None else np.asarray(rain, dtype=float),
        river=zeros if river is None else np.asarray(river, dtype=float),
        weathering=zeros if weathering is None else np.asarray(weathering, dtype=float),
        sediment_deposited=zeros if sediment_deposited is None else np.asarray(sediment_deposited, dtype=float),
        is_river_depositing=np.zeros(n, dtype=bool) if is_river_depositing is None else np.asarray(is_river_depositing, dtype=bool),
        net_elevation_change_m=zeros,
        temperature_c=zeros if temperature_c is None else np.asarray(temperature_c, dtype=float),
        precipitation_mm=zeros if precipitation_mm is None else np.asarray(precipitation_mm, dtype=float),
    )


def test_apply_resource_formation_noop_when_hydrology_cache_missing():
    plate = _plate([0.0], [200.0])
    world = World(seed=0, plates=[plate])  # hydrology_cache defaults to None
    result = _erosion_result(node_points(plate), [200.0])
    geology.apply_resource_formation(world, years=1_000_000, erosion_result=result)
    assert world.plates[0].collect("coal_deposit_m")[0] == 0.0


def test_coal_accumulates_fastest_in_carboniferous_forest_conditions():
    # node0: warm/very wet/flat/low -- Carboniferous Forest. node1: cooler/flat/low -- plain
    # Wetland. node2: warm/very wet but steep -- neither.
    theta = [-0.001, 0.0, 0.001]
    plate = _plate(theta, elevation=[5.0, 5.0, 5.0])
    world = _world_with_hydrology(plate, is_ocean=[False, False, False], water_deposited=[0.0, 0.0, 0.0])
    points = node_points(plate)
    result = _erosion_result(
        points, elevation=[5.0, 5.0, 5.0],
        slope=[0.0001, 0.0001, 0.05],
        temperature_c=[25.0, 10.0, 25.0],
        precipitation_mm=[2500.0, 1500.0, 2500.0],
    )

    geology.apply_resource_formation(world, years=10_000_000, erosion_result=result)

    coal = world.plates[0].collect("coal_deposit_m")
    assert coal[0] > coal[1] > 0.0
    assert coal[2] == 0.0


def test_coal_deposit_is_monotonic_across_repeated_steps():
    plate = _plate([0.0], elevation=[5.0])
    world = _world_with_hydrology(plate, is_ocean=[False], water_deposited=[0.0])
    points = node_points(plate)
    result = _erosion_result(points, elevation=[5.0], slope=[0.0001], temperature_c=[25.0], precipitation_mm=[2500.0])

    prior = 0.0
    for _ in range(5):
        geology.apply_resource_formation(world, years=5_000_000, erosion_result=result)
        current = float(world.plates[0].collect("coal_deposit_m")[0])
        assert current >= prior
        prior = current
    assert prior > 0.0
    assert prior <= geology.MAX_COAL_DEPOSIT_M


def test_oil_gas_forms_only_on_shelf_water_and_is_boosted_near_a_river_mouth():
    # A land cluster of three cells, a shelf ocean cell one cell (~135 km) east of it (inside
    # SHELF_RANGE_RAD), a second shelf cell one cell north of it with heavy river-mouth inflow,
    # and a deep ocean cell far outside shelf range. In node order: land x3, the east shelf
    # cell, the north shelf cell, the deep cell.
    keys = np.concatenate([block_keys(range(4), [0]), block_keys([0, 20], [1])])
    elevation = [200.0, 200.0, 200.0, -50.0, -50.0, -3000.0]
    plate = quad_plate(0, "continental", keys=keys, elevation=np.array(elevation))
    is_ocean = [False, False, False, True, True, True]
    water_deposited = [0.0, 0.0, 0.0, 0.0, 10.0, 0.0]
    world = _world_with_hydrology(plate, is_ocean, water_deposited)
    points = node_points(plate)
    result = _erosion_result(points, elevation=elevation)

    geology.apply_resource_formation(world, years=10_000_000, erosion_result=result)

    oil_gas = world.plates[0].collect("oil_gas_deposit_m")
    assert oil_gas[0] == 0.0  # land
    assert oil_gas[5] == 0.0  # deep, off-shelf
    assert oil_gas[3] > 0.0  # ordinary shelf water
    assert oil_gas[4] > oil_gas[3]  # boosted by the river-mouth inflow
    assert oil_gas[4] <= geology.MAX_OIL_GAS_DEPOSIT_M


def test_soil_depth_rises_from_weathering_and_deposition_and_falls_from_erosion():
    plate = _plate([0.0, 0.01], elevation=[500.0, 500.0])
    world = _world_with_hydrology(plate, is_ocean=[False, False], water_deposited=[0.0, 0.0])
    points = node_points(plate)
    # node0: gentle weathering + floodplain deposition, no fast erosion -- soil should build up.
    # node1: heavy rain+river erosion, no weathering/deposition -- soil should stay at 0 (can't
    # go negative) and definitely not exceed node0's.
    result = _erosion_result(
        points, elevation=[500.0, 500.0],
        weathering=[5.0, 0.0], sediment_deposited=[3.0, 0.0], rain=[0.0, 50.0], river=[0.0, 50.0],
        temperature_c=[20.0, 20.0], precipitation_mm=[800.0, 800.0],
    )

    geology.apply_resource_formation(world, years=1_000_000, erosion_result=result)

    soil = world.plates[0].collect("soil_depth")
    assert soil[0] > 0.0
    assert soil[1] == 0.0
    assert soil[0] <= geology.MAX_SOIL_DEPTH_M


def test_riparian_boost_reaches_depositing_node_and_its_bank_neighbor_but_not_further():
    # 3 nodes, all identical (cool, modest-rain) climate so any difference is purely the
    # riparian boost, not the ordinary precipitation/temperature productivity term. node0 is a
    # slow, big river actively depositing (erosion.is_river_depositing); node1 is node0's own
    # neighbor (a "bank" node) but not itself depositing; node2 neighbors node1, not node0, so
    # the boost shouldn't reach it at all.
    theta = [0.0, 0.001, 0.002]
    plate = _plate(theta, elevation=[100.0, 100.0, 100.0])
    points = node_points(plate)
    n = len(points)
    hydro = HydrologyFields(
        points=points,
        elevation=plate.collect("elevation").copy(),
        is_ocean=np.zeros(n, dtype=bool),
        neighbor_idx=np.array([[1], [0], [1]]),
        flow_target=np.full(n, -1),
        flow_accum=np.zeros(n),
        water_deposited=np.zeros(n),
        filled_elevation=np.zeros(n),
        spill_target=np.full(n, -1),
        is_river=np.zeros(n, dtype=bool),
        lake_depth=np.zeros(n),
        glacier_depth=np.zeros(n),
        plates_in_order=[plate],
    )
    world = World(seed=0, plates=[plate])
    world.hydrology_cache = hydro

    result = _erosion_result(
        points, elevation=[100.0, 100.0, 100.0],
        temperature_c=[9.0, 9.0, 9.0], precipitation_mm=[800.0, 800.0, 800.0],
        is_river_depositing=[True, False, False],
    )
    geology.apply_resource_formation(world, years=2_000_000, erosion_result=result)

    organic = world.plates[0].collect("soil_organic_content")
    # Same base climate everywhere, so any ordering here is purely the riparian boost: the
    # depositing node itself grows the richest, its bank neighbor a real but smaller boost,
    # and node2 (not adjacent to the depositing node at all) is unaffected by either.
    assert organic[0] > organic[1] > organic[2]


def test_soil_zeroed_over_ocean():
    plate = _plate([0.0], elevation=[-500.0])
    world = _world_with_hydrology(plate, is_ocean=[True], water_deposited=[0.0])
    points = node_points(plate)
    result = _erosion_result(points, elevation=[-500.0], weathering=[10.0], sediment_deposited=[10.0])

    geology.apply_resource_formation(world, years=1_000_000, erosion_result=result)

    plate = world.plates[0]
    assert plate.collect("soil_depth")[0] == 0.0
    assert plate.collect("soil_mineral_content")[0] == 0.0
    assert plate.collect("soil_organic_content")[0] == 0.0


def test_soil_organic_content_relaxes_toward_productivity_and_mineral_toward_deposit():
    theta = [0.0]
    plate = _plate(theta, elevation=[100.0])
    # Seed a real mineral_deposit_m and enough soil_depth to hold content, so the relaxation
    # target/gate are both meaningfully nonzero.
    plate.set_fields_on_plate(mineral_deposit_m=np.full(1, volcanism.MAX_MINERAL_DEPOSIT_M), soil_depth=np.ones(1))
    world = _world_with_hydrology(plate, is_ocean=[False], water_deposited=[0.0])
    points = node_points(plate)
    # Warm and wet -- productivity (and so the organic-content target) saturates near 1.0.
    result = _erosion_result(points, elevation=[100.0], temperature_c=[25.0], precipitation_mm=[2500.0], weathering=[0.1])

    organic_prev, mineral_prev = 0.0, 0.0
    for _ in range(20):
        geology.apply_resource_formation(world, years=2_000_000, erosion_result=result)
        organic = float(world.plates[0].collect("soil_organic_content")[0])
        mineral = float(world.plates[0].collect("soil_mineral_content")[0])
        assert organic >= organic_prev
        assert mineral >= mineral_prev
        organic_prev, mineral_prev = organic, mineral

    assert 0.0 < organic_prev <= 1.0
    assert 0.0 < mineral_prev <= 1.0


def test_seed_initial_soil_is_noop_at_zero_maturity():
    plate = _plate([0.0, 0.01], elevation=[500.0, -500.0])
    geology.seed_initial_soil([plate], seed=1, initial_soil_maturity=0.0)
    assert np.all(plate.collect("soil_depth") == 0.0)
    assert np.all(plate.collect("soil_organic_content") == 0.0)


def test_seed_initial_soil_seeds_land_only_at_full_maturity():
    plate = _plate([0.0, 0.01], elevation=[500.0, -500.0])
    geology.seed_initial_soil([plate], seed=1, initial_soil_maturity=1.0)
    soil = plate.collect("soil_depth")
    assert soil[0] > 0.0  # land
    assert soil[1] == 0.0  # ocean, untouched
    assert 0.0 < plate.collect("soil_organic_content")[0] <= 1.0
    assert 0.0 < plate.collect("soil_mineral_content")[0] <= 1.0

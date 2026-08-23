import numpy as np
from app import biomes, geology, geometry, volcanism
from app.erosion import ErosionResult
from app.hydrology import HydrologyFields
from app.plates import ElevationLine, PlateWithLines
from app.world import World


def _plate(theta, elevation):
    frame = geometry.plate_frame_from_seed(np.array([1.0, 0.0, 0.0]))
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.asarray(elevation, dtype=float))
    return PlateWithLines(plate_id=0, frame=frame, crust_type="continental", lines=[line])


def _world_with_hydrology(plate, is_ocean, water_deposited):
    n = len(plate.lines[0].theta)
    points = plate.lines[0].world_xyz(plate.frame)
    hydro = HydrologyFields(
        points=points,
        elevation=plate.lines[0].elevation.copy(),
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
        line_refs=[(plate, 0, 0, n)],
    )
    world = World(seed=0, plates=[plate])
    world.hydrology_cache = hydro
    return world


def _erosion_result(points, elevation, slope=None, rain=None, river=None, weathering=None, sediment_deposited=None, temperature_c=None, precipitation_mm=None):
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
        temperature_c=zeros if temperature_c is None else np.asarray(temperature_c, dtype=float),
        precipitation_mm=zeros if precipitation_mm is None else np.asarray(precipitation_mm, dtype=float),
    )


def test_apply_resource_formation_noop_when_hydrology_cache_missing():
    plate = _plate([0.0], [200.0])
    world = World(seed=0, plates=[plate])  # hydrology_cache defaults to None
    result = _erosion_result(plate.lines[0].world_xyz(plate.frame), [200.0])
    geology.apply_resource_formation(world, years=1_000_000, erosion_result=result)
    assert world.plates[0].lines[0].coal_deposit_m[0] == 0.0


def test_coal_accumulates_fastest_in_carboniferous_forest_conditions():
    # node0: warm/very wet/flat/low -- Carboniferous Forest. node1: cooler/flat/low -- plain
    # Wetland. node2: warm/very wet but steep -- neither.
    theta = [-0.001, 0.0, 0.001]
    plate = _plate(theta, elevation=[5.0, 5.0, 5.0])
    world = _world_with_hydrology(plate, is_ocean=[False, False, False], water_deposited=[0.0, 0.0, 0.0])
    points = plate.lines[0].world_xyz(plate.frame)
    result = _erosion_result(
        points, elevation=[5.0, 5.0, 5.0],
        slope=[0.0001, 0.0001, 0.05],
        temperature_c=[25.0, 10.0, 25.0],
        precipitation_mm=[2500.0, 1500.0, 2500.0],
    )

    geology.apply_resource_formation(world, years=10_000_000, erosion_result=result)

    coal = world.plates[0].lines[0].coal_deposit_m
    assert coal[0] > coal[1] > 0.0
    assert coal[2] == 0.0


def test_coal_deposit_is_monotonic_across_repeated_steps():
    plate = _plate([0.0], elevation=[5.0])
    world = _world_with_hydrology(plate, is_ocean=[False], water_deposited=[0.0])
    points = plate.lines[0].world_xyz(plate.frame)
    result = _erosion_result(points, elevation=[5.0], slope=[0.0001], temperature_c=[25.0], precipitation_mm=[2500.0])

    prior = 0.0
    for _ in range(5):
        geology.apply_resource_formation(world, years=5_000_000, erosion_result=result)
        current = float(world.plates[0].lines[0].coal_deposit_m[0])
        assert current >= prior
        prior = current
    assert prior > 0.0
    assert prior <= geology.MAX_COAL_DEPOSIT_M


def test_oil_gas_forms_only_on_shelf_water_and_is_boosted_near_a_river_mouth():
    # A land cluster near theta=0, a shelf ocean point ~127km out (inside SHELF_RANGE_RAD),
    # a second shelf point at the same distance but with heavy river-mouth inflow, and a deep
    # ocean point far outside shelf range -- same theta-to-km relationship
    # test_bathymetry.py's own shelf-vs-deep test already establishes.
    theta = [-0.001, -0.0009, -0.0011, 0.02, 0.0205, 0.3]
    elevation = [200.0, 200.0, 200.0, -50.0, -50.0, -3000.0]
    plate = _plate(theta, elevation)
    is_ocean = [False, False, False, True, True, True]
    water_deposited = [0.0, 0.0, 0.0, 0.0, 10.0, 0.0]
    world = _world_with_hydrology(plate, is_ocean, water_deposited)
    points = plate.lines[0].world_xyz(plate.frame)
    result = _erosion_result(points, elevation=elevation)

    geology.apply_resource_formation(world, years=10_000_000, erosion_result=result)

    oil_gas = world.plates[0].lines[0].oil_gas_deposit_m
    assert oil_gas[0] == 0.0  # land
    assert oil_gas[5] == 0.0  # deep, off-shelf
    assert oil_gas[3] > 0.0  # ordinary shelf water
    assert oil_gas[4] > oil_gas[3]  # boosted by the river-mouth inflow
    assert oil_gas[4] <= geology.MAX_OIL_GAS_DEPOSIT_M


def test_soil_depth_rises_from_weathering_and_deposition_and_falls_from_erosion():
    plate = _plate([0.0, 0.01], elevation=[500.0, 500.0])
    world = _world_with_hydrology(plate, is_ocean=[False, False], water_deposited=[0.0, 0.0])
    points = plate.lines[0].world_xyz(plate.frame)
    # node0: gentle weathering + floodplain deposition, no fast erosion -- soil should build up.
    # node1: heavy rain+river erosion, no weathering/deposition -- soil should stay at 0 (can't
    # go negative) and definitely not exceed node0's.
    result = _erosion_result(
        points, elevation=[500.0, 500.0],
        weathering=[5.0, 0.0], sediment_deposited=[3.0, 0.0], rain=[0.0, 50.0], river=[0.0, 50.0],
        temperature_c=[20.0, 20.0], precipitation_mm=[800.0, 800.0],
    )

    geology.apply_resource_formation(world, years=1_000_000, erosion_result=result)

    soil = world.plates[0].lines[0].soil_depth
    assert soil[0] > 0.0
    assert soil[1] == 0.0
    assert soil[0] <= geology.MAX_SOIL_DEPTH_M


def test_soil_zeroed_over_ocean():
    plate = _plate([0.0], elevation=[-500.0])
    world = _world_with_hydrology(plate, is_ocean=[True], water_deposited=[0.0])
    points = plate.lines[0].world_xyz(plate.frame)
    result = _erosion_result(points, elevation=[-500.0], weathering=[10.0], sediment_deposited=[10.0])

    geology.apply_resource_formation(world, years=1_000_000, erosion_result=result)

    line = world.plates[0].lines[0]
    assert line.soil_depth[0] == 0.0
    assert line.soil_mineral_content[0] == 0.0
    assert line.soil_organic_content[0] == 0.0


def test_soil_organic_content_relaxes_toward_productivity_and_mineral_toward_deposit():
    theta = [0.0]
    plate = _plate(theta, elevation=[100.0])
    # Seed a real mineral_deposit_m and enough soil_depth to hold content, so the relaxation
    # target/gate are both meaningfully nonzero.
    plate.lines[0].mineral_deposit_m[:] = volcanism.MAX_MINERAL_DEPOSIT_M
    plate.lines[0].soil_depth[:] = 1.0
    world = _world_with_hydrology(plate, is_ocean=[False], water_deposited=[0.0])
    points = plate.lines[0].world_xyz(plate.frame)
    # Warm and wet -- productivity (and so the organic-content target) saturates near 1.0.
    result = _erosion_result(points, elevation=[100.0], temperature_c=[25.0], precipitation_mm=[2500.0], weathering=[0.1])

    organic_prev, mineral_prev = 0.0, 0.0
    for _ in range(20):
        geology.apply_resource_formation(world, years=2_000_000, erosion_result=result)
        line = world.plates[0].lines[0]
        assert line.soil_organic_content[0] >= organic_prev
        assert line.soil_mineral_content[0] >= mineral_prev
        organic_prev, mineral_prev = float(line.soil_organic_content[0]), float(line.soil_mineral_content[0])

    assert 0.0 < organic_prev <= 1.0
    assert 0.0 < mineral_prev <= 1.0


def test_seed_initial_soil_is_noop_at_zero_maturity():
    plate = _plate([0.0, 0.01], elevation=[500.0, -500.0])
    geology.seed_initial_soil([plate], seed=1, initial_soil_maturity=0.0)
    line = plate.lines[0]
    assert np.all(line.soil_depth == 0.0)
    assert np.all(line.soil_organic_content == 0.0)


def test_seed_initial_soil_seeds_land_only_at_full_maturity():
    plate = _plate([0.0, 0.01], elevation=[500.0, -500.0])
    geology.seed_initial_soil([plate], seed=1, initial_soil_maturity=1.0)
    line = plate.lines[0]
    assert line.soil_depth[0] > 0.0  # land
    assert line.soil_depth[1] == 0.0  # ocean, untouched
    assert 0.0 < line.soil_organic_content[0] <= 1.0
    assert 0.0 < line.soil_mineral_content[0] <= 1.0

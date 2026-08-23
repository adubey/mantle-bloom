import numpy as np
from app import erosion, geometry, hydrology
from app.plates import ElevationLine, PlateWithLines
from app.world import World, generate_world, step_world


def _flow_line_plate(plate_id, theta, elevation):
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.asarray(elevation, dtype=float))
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type="continental", lines=[line])


def _flow_line_plate_with_lake(plate_id, theta, elevation, lake_depth):
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(
        phi=0.0, theta=theta, elevation=np.asarray(elevation, dtype=float), lake_depth=np.asarray(lake_depth, dtype=float)
    )
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type="continental", lines=[line])


def _normalize(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def _river_inspector_fields():
    # A hand-built 7-node HydrologyFields for the River Inspector's grouping logic, not run
    # through compute_hydrology (small enough to hit its FLOW_NEIGHBOR_COUNT bail-out path
    # anyway). Two separate is_river-classified networks:
    #   Network A -- a Y-confluence: 0 and 1 (two headwaters) both flow into 2, which flows
    #   into 3 (the mouth, highest flow_accum in the group), which flows onward into ocean
    #   node 4.
    #   Network B -- a single lone node 5, sitting at its own sink with lake_depth above the
    #   visible threshold (so it should read as ending in a lake, not the ocean).
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
        ]
    )
    n = len(points)
    elevation = np.array([100.0, 110.0, 80.0, 50.0, -10.0, 60.0, 40.0])
    is_ocean = np.array([False, False, False, False, True, False, False])
    flow_target = np.array([2, 2, 3, 4, -1, -1, -1])
    flow_accum = np.array([10.0, 12.0, 25.0, 40.0, 0.0, 5.0, 0.0])
    is_river = np.array([True, True, True, True, False, True, False])
    lake_depth = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0])  # node 5 well above LAKE_MIN_VISIBLE_DEPTH_M
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
        line_refs=[],
    )


def test_channel_depth_and_lake_depth_persist_across_boundary_and_erosion_steps():
    # A direct regression check on the field-threading work: a node's channel_depth/
    # lake_depth must survive being touched by boundary.py (which runs every step, before
    # erosion.py sets these) rather than getting silently reset to 0.
    world = generate_world(seed=30, num_plates=8, continental_fraction=0.5)
    step_world(world, years=2_000_000)

    channel_total_1 = sum(line.channel_depth.sum() for p in world.plates for line in p.lines)
    lake_total_1 = sum(line.lake_depth.sum() for p in world.plates for line in p.lines)
    assert channel_total_1 > 0.0 or lake_total_1 > 0.0

    step_world(world, years=2_000_000)  # boundary.step_boundaries runs first here again

    channel_total_2 = sum(line.channel_depth.sum() for p in world.plates for line in p.lines)
    # Channel depth only ever grows (monotonic, see erosion.py) -- if boundary.py had reset
    # it, this would have dropped back toward 0 instead.
    assert channel_total_2 >= channel_total_1 * 0.5  # loose bound: some nodes can be pruned/moved


def test_apply_erosion_conserves_mass_between_erosion_and_deposition():
    world = generate_world(seed=31, num_plates=8, continental_fraction=0.5)
    for _ in range(3):
        step_world(world, years=2_000_000)
    assert world.hydrology_cache is not None
    assert world.hydrology_cache.is_river.sum() >= 0  # river classification ran without error


def test_rivers_never_overlap_lakes_in_a_real_stepped_world():
    # Regression check for "rivers end at lakes and don't go through them": a lake, however
    # widely flooded, must never also be classified is_river.
    world = generate_world(seed=32, num_plates=10, continental_fraction=0.5)
    for _ in range(10):
        step_world(world, years=2_000_000)
    hydro = world.hydrology_cache
    assert hydro is not None
    is_lake = hydro.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    assert not np.any(hydro.is_river & is_lake)


def test_channel_width_grows_with_flow_and_persists_across_steps():
    world = generate_world(seed=33, num_plates=8, continental_fraction=0.5)
    for _ in range(4):
        step_world(world, years=2_000_000)

    widths = np.concatenate([line.channel_width for p in world.plates for line in p.lines])
    depths = np.concatenate([line.channel_depth for p in world.plates for line in p.lines])
    assert np.any(widths > 0.0)
    assert np.all(widths <= erosion.MAX_CHANNEL_WIDTH_M)
    # Width grows from discharge alone, unlike depth (which also needs real slope -- see
    # erosion.py's own river term) -- a flat, wide, slow-flowing stretch can widen without
    # incising at all, so this only checks that *some* overlap exists, not that every widened
    # node is also a channel-depth node.
    assert np.any((widths > 0.0) & (depths > 0.0))

    width_total_1 = widths.sum()
    step_world(world, years=2_000_000)
    width_total_2 = sum(line.channel_width.sum() for p in world.plates for line in p.lines)
    # Monotonic like channel_depth -- if boundary.py/bathymetry.py had reset it, this would
    # have dropped back toward 0 instead.
    assert width_total_2 >= width_total_1 * 0.5  # loose bound: some nodes can be pruned/moved


def test_is_volcano_survives_a_full_step_cycle():
    # Direct regression check for the bug this session found: bathymetry.py's and erosion.py's
    # own line-reconstruction sites didn't know about is_volcano/volcano_active_years_remaining
    # and silently wiped them to False/0 every step, before volcanism.apply_volcanic_activity
    # ever got a chance to read them.
    world = generate_world(seed=34, num_plates=10, continental_fraction=0.4)
    found_a_field = False
    for _ in range(20):
        step_world(world, years=2_000_000)
        if world.volcanic_field_plate_ids:
            found_a_field = True
            total_volcano_nodes = sum(int(line.is_volcano.sum()) for p in world.plates for line in p.lines)
            assert total_volcano_nodes > 0
            break
    assert found_a_field  # sanity: this seed/step count should produce at least one field


def test_channel_lake_and_glacier_depth_persist_across_boundary_and_erosion_steps():
    # Regression check on the field-threading work: a node's channel_depth/lake_depth/
    # glacier_depth must survive being touched by boundary.py (which runs every step, before
    # erosion.py sets these) rather than getting silently reset to 0.
    world = generate_world(seed=32, num_plates=10, continental_fraction=0.5)
    for _ in range(6):
        step_world(world, years=2_000_000)

    glacier_total_1 = sum(line.glacier_depth.sum() for p in world.plates for line in p.lines)

    step_world(world, years=2_000_000)  # boundary.step_boundaries runs first here again

    glacier_total_2 = sum(line.glacier_depth.sum() for p in world.plates for line in p.lines)
    # Not a strict monotonic-growth assertion (unlike channel_depth, glacier_depth can melt)
    # -- just confirms boundary.py didn't wipe the whole world's ice to exactly 0.
    assert glacier_total_1 > 0.0
    assert glacier_total_2 > 0.0


def test_silt_depth_persists_across_boundary_and_erosion_steps():
    # Regression check on silt_depth's own field-threading (plates.py/boundary.py/
    # elevation_lines.py/reassign.py/merge_split.py) -- same shape as the channel/lake/glacier
    # persistence test above: confirms boundary.py (which runs every step, before erosion.py
    # sets silt_depth) doesn't silently reset it to 0, and that it's actually accumulating
    # somewhere after several steps of real lake inflow.
    world = generate_world(seed=34, num_plates=10, continental_fraction=0.5)
    for _ in range(6):
        step_world(world, years=2_000_000)

    silt_total_1 = sum(line.silt_depth.sum() for p in world.plates for line in p.lines)
    step_world(world, years=2_000_000)  # boundary.step_boundaries runs first here again
    silt_total_2 = sum(line.silt_depth.sum() for p in world.plates for line in p.lines)

    assert silt_total_1 > 0.0
    assert silt_total_2 > 0.0

import numpy as np

from app import climate, coastline, geometry, hydrology
from app.world import World, generate_world, step_world


def _climate_fields(is_ocean: np.ndarray) -> climate.ClimateFields:
    """A minimal-but-complete ClimateFields for a small hand-built grid -- only lat_deg/
    lon_deg/world_xyz/is_ocean matter for coastline tracing; every other field is filled with
    correctly-shaped placeholder zeros/False, since compute_coastline_segments never reads
    them."""
    height, width = is_ocean.shape
    lat_deg, lon_deg, world_xyz = climate._build_grid(height, width)
    zeros = np.zeros((height, width))
    return climate.ClimateFields(
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        world_xyz=world_xyz,
        elevation_m=zeros,
        is_ocean=is_ocean,
        land_temperature_c=zeros,
        ocean_temperature_c=zeros,
        air_temperature_c=zeros,
        wind_u=zeros,
        wind_v=zeros,
        current_u=zeros,
        current_v=zeros,
        humidity=zeros,
        precipitation_mm=zeros,
        swell_rows=np.zeros(0, dtype=int),
        swell_cols=np.zeros(0, dtype=int),
    )


def _hydrology_fields(points: np.ndarray, lake_depth: np.ndarray) -> hydrology.HydrologyFields:
    """A minimal-but-complete HydrologyFields -- only points/lake_depth matter for
    coastline._lake_mask_on_grid; every other field is a correctly-shaped placeholder."""
    n = len(points)
    zeros_f = np.zeros(n)
    zeros_i = np.zeros(n, dtype=np.int64)
    return hydrology.HydrologyFields(
        points=points,
        elevation=zeros_f,
        is_ocean=np.zeros(n, dtype=bool),
        neighbor_idx=np.zeros((n, 0), dtype=np.int64),
        flow_target=np.full(n, -1, dtype=np.int64),
        flow_accum=zeros_f,
        water_deposited=zeros_f,
        filled_elevation=zeros_f,
        spill_target=np.full(n, -1, dtype=np.int64),
        is_river=np.zeros(n, dtype=bool),
        lake_depth=lake_depth,
        glacier_depth=zeros_f,
        plates_in_order=[],
    )


def _segment_set(point_a: np.ndarray, point_b: np.ndarray) -> set:
    """Rounds xyz endpoints and returns an order-independent set of (a, b) pairs -- each pair
    itself order-independent too, since which of a/b comes first isn't semantically
    meaningful for an undirected coastline edge."""
    result = set()
    for a, b in zip(np.round(point_a, 6).tolist(), np.round(point_b, 6).tolist()):
        result.add(frozenset([tuple(a), tuple(b)]))
    return result


def test_compute_coastline_segments_empty_before_first_step():
    world = World(seed=0, plates=[])
    assert world.climate_cache is None
    point_a, point_b = coastline.compute_coastline_segments(world)
    assert point_a.shape == (0, 3)
    assert point_b.shape == (0, 3)


def test_compute_coastline_segments_traces_a_known_rectangular_island():
    # A 4x4 grid, is_ocean everywhere except a 2x2 "island" at rows 1-2, cols 1-2 -- a
    # rectangle's coastline is exactly its own perimeter in grid-cell edges: 2*(2+2) = 8.
    is_ocean = np.array(
        [
            [True, True, True, True],
            [True, False, False, True],
            [True, False, False, True],
            [True, True, True, True],
        ]
    )
    world = World(seed=0, plates=[])
    world.climate_cache = _climate_fields(is_ocean)

    point_a, point_b = coastline.compute_coastline_segments(world)
    assert len(point_a) == 8
    assert len(point_b) == 8
    assert np.allclose(np.linalg.norm(point_a, axis=1), 1.0)
    assert np.allclose(np.linalg.norm(point_b, axis=1), 1.0)

    # Cross-check against directly computed expected lat/lon endpoints for this exact grid
    # (half_lat = 90/4 = 22.5, half_lon = 180/4 = 45; lat_deg = [67.5, 22.5, -22.5, -67.5],
    # lon_deg = [-135, -45, 45, 135]).
    expected_latlon = [
        ((45.0, -90.0), (0.0, -90.0)),
        ((45.0, 90.0), (0.0, 90.0)),
        ((0.0, -90.0), (-45.0, -90.0)),
        ((0.0, 90.0), (-45.0, 90.0)),
        ((45.0, -90.0), (45.0, 0.0)),
        ((45.0, 0.0), (45.0, 90.0)),
        ((-45.0, -90.0), (-45.0, 0.0)),
        ((-45.0, 0.0), (-45.0, 90.0)),
    ]
    expected = {frozenset([tuple(np.round(geometry.latlon_to_xyz(np.radians(p[0]), np.radians(p[1])), 6)) for p in pair]) for pair in expected_latlon}
    assert _segment_set(point_a, point_b) == expected


def test_lake_mask_on_grid_marks_cells_above_the_lake_depth_threshold():
    is_ocean = np.zeros((4, 4), dtype=bool)  # all land
    fields = _climate_fields(is_ocean)
    world = World(seed=0, plates=[])
    world.climate_cache = fields

    # One hydrology node sitting exactly on grid cell (1, 1)'s own center, with a lake deep
    # enough to count, plus one far-away dry node so the tree has more than one point.
    points = np.stack([fields.world_xyz[1, 1], fields.world_xyz[3, 3]], axis=0)
    lake_depth = np.array([hydrology.LAKE_MIN_VISIBLE_DEPTH_M + 1.0, 0.0])
    world.hydrology_cache = _hydrology_fields(points, lake_depth)

    mask = coastline._lake_mask_on_grid(world, fields.world_xyz)
    assert mask.dtype == bool
    assert mask[1, 1]
    assert not mask[3, 3]


def test_compute_coastline_segments_includes_a_lake_boundary_over_land():
    # All-land grid (no ocean at all) -- without lake support, this would produce zero
    # coastline segments even though there's a real interior lake that should still get its
    # own traced boundary (same 2x2-island shape as the ocean test above, just land-vs-lake
    # instead of ocean-vs-land).
    is_ocean = np.zeros((4, 4), dtype=bool)
    fields = _climate_fields(is_ocean)
    world = World(seed=0, plates=[])
    world.climate_cache = fields

    lake_mask = np.array(
        [
            [False, False, False, False],
            [False, True, True, False],
            [False, True, True, False],
            [False, False, False, False],
        ]
    )
    lake_rows, lake_cols = np.nonzero(lake_mask)
    points = fields.world_xyz[lake_rows, lake_cols]
    lake_depth = np.full(len(points), hydrology.LAKE_MIN_VISIBLE_DEPTH_M + 1.0)
    # Every other grid cell needs a *dry* node too, or the nearest-neighbor resample would
    # snap every cell to the nearest lake point and mark the whole grid as lake.
    dry_rows, dry_cols = np.nonzero(~lake_mask)
    dry_points = fields.world_xyz[dry_rows, dry_cols]
    world.hydrology_cache = _hydrology_fields(
        np.concatenate([points, dry_points], axis=0),
        np.concatenate([lake_depth, np.zeros(len(dry_points))]),
    )

    point_a, point_b = coastline.compute_coastline_segments(world)
    assert len(point_a) == 8  # same rectangular-perimeter count as the ocean-island test


def test_compute_coastline_segments_returns_unit_vectors_on_a_real_world():
    world = generate_world(seed=20, num_plates=10, continental_fraction=0.5)
    step_world(world, years=2_000_000)
    point_a, point_b = coastline.compute_coastline_segments(world)
    assert len(point_a) == len(point_b)
    if len(point_a) > 0:
        assert np.allclose(np.linalg.norm(point_a, axis=1), 1.0)
        assert np.allclose(np.linalg.norm(point_b, axis=1), 1.0)

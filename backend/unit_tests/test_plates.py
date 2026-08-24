import numpy as np
from app import geometry
from app.elevation_lines import ElevationLine
from app.plates import (
    ELLIPSE_OUTLINE_POINTS,
    MAX_AUTO_PLATES,
    MIN_AUTO_PLATES,
    MIN_OCEANIC_PLATES,
    NODE_DENSITY_CHOICES,
    PlateWithLines,
    collect_all_coal_deposit,
    collect_all_mineral_deposit,
    collect_all_oil_gas_deposit,
    collect_all_points,
    collect_all_soil_depth,
    collect_all_soil_mineral_content,
    collect_all_soil_organic_content,
    generate_plates,
    nearest_plate_id,
    plate_bounding_ellipse,
)
from app.world import generate_world


def _measured_land_fraction(plates_list) -> float:
    total = sum(p.node_count() for p in plates_list)
    land = sum(int(np.sum(line.elevation > 0)) for p in plates_list for line in p.lines)
    return land / total if total else 0.0


def _measured_continental_area_fraction(plates_list) -> float:
    total = sum(p.node_count() for p in plates_list)
    continental = sum(p.node_count() for p in plates_list if p.crust_type == "continental")
    return continental / total if total else 0.0


def test_generate_plates_node_density_quadruples_node_count():
    reference = generate_plates(seed=7, num_plates=8, node_density=1.0)
    denser = generate_plates(seed=7, num_plates=8, node_density=4.0)
    reference_nodes = sum(p.node_count() for p in reference)
    denser_nodes = sum(p.node_count() for p in denser)
    ratio = denser_nodes / reference_nodes
    assert 3.5 < ratio < 4.5  # not exact -- lattice row/column counts round to integers


def test_node_density_choices_all_produce_a_valid_world():
    for density in NODE_DENSITY_CHOICES:
        plates = generate_plates(seed=3, num_plates=6, node_density=density)
        assert all(p.node_count() > 0 for p in plates)


def test_generate_plates_count_and_crust_types():
    plates = generate_plates(seed=42, num_plates=10)
    assert len(plates) == 10
    assert all(p.crust_type in ("continental", "oceanic") for p in plates)


def test_every_plate_has_elevation_lines():
    plates = generate_plates(seed=1, num_plates=8)
    for p in plates:
        assert p.node_count() > 0, f"plate {p.plate_id} has no elevation nodes"


def test_frames_are_proper_rotations():
    plates = generate_plates(seed=2, num_plates=6)
    for p in plates:
        assert np.allclose(p.frame @ p.frame.T, np.eye(3), atol=1e-9)
        assert np.isclose(np.linalg.det(p.frame), 1.0)


def test_every_node_is_closest_to_its_own_plate_seed():
    """A node kept by a plate must actually be in that plate's Voronoi cell -- i.e. no
    other plate's seed is angularly closer to it."""
    plates = generate_plates(seed=3, num_plates=8)
    seeds = np.array([p.seed_world for p in plates])

    for p in plates:
        for line in p.lines:
            world_pts = line.world_xyz(p.frame)
            dists = geometry.angular_distance(world_pts[:, None, :], seeds[None, :, :])
            nearest = np.argmin(dists, axis=1)
            assert np.all(nearest == p.plate_id)


def test_lines_are_evenly_spaced_in_phi():
    # node_density pinned to 1.0 (not DEFAULT_NODE_DENSITY): the spacing check below compares
    # against the reference TARGET_LINE_SPACING_RAD, which only holds at density 1.0 -- at any
    # other density, actual line spacing is TARGET_LINE_SPACING_RAD / sqrt(node_density) (see
    # line_spacing_rad), and this test only cares that spacing is *even*, not what density
    # produced it.
    plates = generate_plates(seed=4, num_plates=6, node_density=1.0)
    for p in plates:
        phis = sorted(line.phi for line in p.lines)
        if len(phis) < 2:
            continue
        diffs = np.diff(phis)
        # All gaps should be an integer multiple of the target spacing (some plates
        # won't own every consecutive row near their boundary).
        from app.elevation_lines import TARGET_LINE_SPACING_RAD

        ratios = diffs / TARGET_LINE_SPACING_RAD
        assert np.allclose(ratios, np.round(ratios), atol=1e-6)


def test_generate_world_matches_plate_count():
    world = generate_world(seed=7, num_plates=9)
    assert len(world.plates) == 9
    assert world.elapsed_years == 0.0


def test_generation_is_deterministic_for_same_seed():
    w1 = generate_world(seed=123, num_plates=8)
    w2 = generate_world(seed=123, num_plates=8)
    assert len(w1.plates) == len(w2.plates)
    for p1, p2 in zip(w1.plates, w2.plates):
        assert p1.crust_type == p2.crust_type
        assert np.allclose(p1.frame, p2.frame)
        assert len(p1.lines) == len(p2.lines)


def test_generate_plates_auto_count_is_deterministic_for_same_seed():
    p1 = generate_plates(seed=99)
    p2 = generate_plates(seed=99)
    assert len(p1) == len(p2)
    for a, b in zip(p1, p2):
        assert a.crust_type == b.crust_type
        assert np.allclose(a.frame, b.frame)


def test_generate_plates_continental_fraction_gives_exact_continental_count():
    for n in range(1, 6):
        # n / 12 divides evenly, so round() introduces no rounding-error ambiguity here.
        plates = generate_plates(seed=3, num_plates=12, continental_fraction=n / 12)
        continental = [p for p in plates if p.crust_type == "continental"]
        assert len(continental) == n


def test_generate_plates_continental_fraction_bumps_up_total_plate_count_if_needed():
    plates = generate_plates(seed=4, num_plates=5, continental_fraction=1.0)
    assert len(plates) >= 5 + MIN_OCEANIC_PLATES
    assert sum(1 for p in plates if p.crust_type == "continental") == 5


def test_generate_plates_continental_fraction_is_clamped_to_one():
    plates = generate_plates(seed=5, num_plates=10, continental_fraction=999.0)
    continental = sum(1 for p in plates if p.crust_type == "continental")
    assert continental == 10  # clamped to 1.0 -> round(1.0 * 10), not literally 999 plates
    assert len(plates) == 10 + MIN_OCEANIC_PLATES


def test_generate_plates_land_fraction_matches_target_when_achievable():
    # 70% continental plates leaves comfortably more continental area than 29% land needs,
    # so the target should land almost exactly (bounded only by the sampling in
    # _land_noise_threshold, not by running out of continental crust to place land on).
    plates = generate_plates(seed=2, num_plates=14, continental_fraction=0.7, land_fraction=0.29)
    assert abs(_measured_land_fraction(plates) - 0.29) < 0.02


def test_generate_plates_land_fraction_is_capped_by_continental_area():
    # Only ~1/4 of plates (by count, roughly by area too) are continental, so 80% land is
    # not achievable -- every continental node should end up as land (elevation > 0) and no
    # more, capping measured land at roughly the continental area fraction itself.
    plates = generate_plates(seed=2, num_plates=14, continental_fraction=0.25, land_fraction=0.8)
    continental_area = _measured_continental_area_fraction(plates)
    assert abs(_measured_land_fraction(plates) - continental_area) < 0.02


def test_generate_plates_land_fraction_zero_gives_no_land():
    # Not an exact 0.0: the threshold is estimated from a coarser whole-sphere sample
    # (LAND_FRACTION_SAMPLE_SPACING_KM) than the actual plate lattice it's applied to, so a
    # handful of real nodes can have a noise value fractionally above that sample's max --
    # negligible (a few thousandths of a percent), not a sign the target was ignored.
    plates = generate_plates(seed=2, num_plates=14, continental_fraction=0.7, land_fraction=0.0)
    assert _measured_land_fraction(plates) < 0.001


def test_outline_world_traces_a_loop_covering_every_line():
    plates = generate_plates(seed=5, num_plates=8)
    for p in plates:
        lines_with_nodes = [line for line in p.lines if len(line.theta) > 0]
        outline = p.outline_world()
        assert len(outline) == 2 * len(lines_with_nodes)
        assert np.allclose(np.linalg.norm(outline, axis=-1), 1.0, atol=1e-9)


def test_outline_world_empty_for_plate_with_no_lines():
    plates = generate_plates(seed=6, num_plates=8)
    p = plates[0]
    p.set_lines([])
    assert len(p.outline_world()) == 0


def test_get_bounding_polygon_matches_outline_world():
    plates = generate_plates(seed=7, num_plates=8)
    for p in plates:
        assert np.array_equal(p.get_bounding_polygon(), p.outline_world())


def test_get_bounding_polygon_returns_the_same_cached_array_until_invalidated():
    p = generate_plates(seed=8, num_plates=8)[0]
    first = p.get_bounding_polygon()
    second = p.get_bounding_polygon()
    assert first is second  # same cached object, not recomputed


def test_get_bounding_polygon_cache_invalidated_by_rotate():
    p = generate_plates(seed=9, num_plates=8)[0]
    cached = p.get_bounding_polygon()
    p.rotate(geometry.plate_frame_from_seed(np.array([0.0, 1.0, 0.0])))
    rotated = p.get_bounding_polygon()
    assert rotated is not cached
    assert np.array_equal(rotated, p.outline_world())


def test_get_bounding_polygon_cache_invalidated_by_set_lines():
    p = generate_plates(seed=10, num_plates=8)[0]
    cached = p.get_bounding_polygon()
    p.set_lines(list(p.lines))
    refreshed = p.get_bounding_polygon()
    assert refreshed is not cached
    assert np.array_equal(refreshed, cached)  # same lines, so same outline -- just recomputed


def test_get_bounding_polygon_cache_invalidated_by_replace_line():
    p = generate_plates(seed=11, num_plates=8)[0]
    cached = p.get_bounding_polygon()
    p.replace_line(0, p.lines[0])
    refreshed = p.get_bounding_polygon()
    assert refreshed is not cached


def test_map_world_points_on_plate_fraction_spans_zero_to_one_along_each_line():
    frame = geometry.plate_frame_from_seed(np.array([1.0, 0.0, 0.0]))
    theta = np.array([0.0, 0.25, 1.0])
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.zeros(3))
    plate = PlateWithLines(plate_id=1, frame=frame, crust_type="continental", lines=[line])

    fractions = [fraction for _, _, fraction in plate.map_world_points_on_plate()]
    assert np.allclose(fractions, [0.0, 0.25, 1.0])


def test_map_world_points_on_plate_single_node_line_gives_midpoint_fraction():
    frame = geometry.plate_frame_from_seed(np.array([0.0, 1.0, 0.0]))
    line = ElevationLine(phi=0.1, theta=np.array([0.3]), elevation=np.zeros(1))
    plate = PlateWithLines(plate_id=2, frame=frame, crust_type="oceanic", lines=[line])

    ((_, _, fraction),) = list(plate.map_world_points_on_plate())
    assert fraction == 0.5


def test_map_world_points_on_plate_matches_map_world_points_order_and_positions():
    world_plates = generate_plates(seed=11, num_plates=6)
    for p in world_plates:
        if p.node_count() == 0:
            continue
        base = list(p.map_world_points())
        extended = list(p.map_world_points_on_plate())
        assert len(base) == len(extended)
        for (point_a, xyz_a), (point_b, xyz_b, fraction) in zip(base, extended):
            assert point_a.get_theta() == point_b.get_theta()
            assert np.allclose(xyz_a, xyz_b)
            assert 0.0 <= fraction <= 1.0


def test_plate_bounding_ellipse_empty_for_no_points():
    assert plate_bounding_ellipse(np.zeros((0, 3))) is None


def test_plate_bounding_ellipse_contains_a_clustered_point_cloud():
    rng = np.random.default_rng(7)
    centroid = geometry.normalize(np.array([1.0, 0.0, 0.0]))
    east, north = geometry.local_tangent_basis(centroid)
    # A compact cluster within ~10 degrees of centroid -- comfortably inside
    # azimuthal-equidistant's well-behaved regime (see its docstring).
    offsets = rng.normal(size=(200, 2)) * np.radians(4.0)
    points = geometry.azimuthal_equidistant_inverse(centroid, east, north, offsets)

    ellipse = plate_bounding_ellipse(points)
    assert ellipse is not None
    assert ellipse.diameter_a_km >= ellipse.diameter_b_km >= 0.0
    assert ellipse.outline_xyz.shape == (ELLIPSE_OUTLINE_POINTS, 3)
    assert np.allclose(np.linalg.norm(ellipse.outline_xyz, axis=-1), 1.0, atol=1e-9)
    assert np.isclose(np.linalg.norm(ellipse.center_xyz), 1.0, atol=1e-9)

    # Every input point should fall within the fitted ellipse, measured the same way the
    # ellipse itself was fit (azimuthal-equidistant km-plane around the point cloud's own
    # bounding_sphere centroid).
    from app.elevation_lines import PLANET_RADIUS_KM

    fit_centroid, _ = geometry.bounding_sphere(points)
    fit_east, fit_north = geometry.local_tangent_basis(fit_centroid)
    xy_km = geometry.azimuthal_equidistant_forward(fit_centroid, fit_east, fit_north, points) * PLANET_RADIUS_KM
    center_km = (
        geometry.azimuthal_equidistant_forward(fit_centroid, fit_east, fit_north, ellipse.center_xyz[None, :])[0]
        * PLANET_RADIUS_KM
    )
    rel = xy_km - center_km
    semi_a = max(ellipse.diameter_a_km / 2.0, 1e-9)
    semi_b = max(ellipse.diameter_b_km / 2.0, 1e-9)
    # Not axis-aligned in general, so bound by the enclosing circle of the larger semi-axis
    # rather than re-deriving the fitted rotation angle here (that's ellipse.py's own test).
    assert np.all(np.hypot(rel[:, 0], rel[:, 1]) <= max(semi_a, semi_b) + 1e-6)


def test_plate_bounding_ellipse_handles_more_than_a_hemisphere():
    """Adversarial case per azimuthal_equidistant_forward's documented antipodal-singularity
    limitation: points spread across more than a hemisphere from their own mean-direction
    centroid. Not asserting a tight/correct fit here (known limitation, not solved for v1)
    -- just that this doesn't crash or produce NaN/Inf, documenting current behavior."""
    rng = np.random.default_rng(8)
    lat = rng.uniform(-np.pi / 2, np.pi / 2, size=60)
    lon = rng.uniform(-np.pi, np.pi, size=60)  # spread across the *entire* sphere
    points = geometry.latlon_to_xyz(lat, lon)

    ellipse = plate_bounding_ellipse(points)
    assert ellipse is not None
    assert np.all(np.isfinite(ellipse.center_xyz))
    assert np.isfinite(ellipse.diameter_a_km)
    assert np.isfinite(ellipse.diameter_b_km)
    assert np.all(np.isfinite(ellipse.outline_xyz))


def test_collect_all_points_concatenates_across_plates():
    world_plates = generate_plates(seed=9, num_plates=6)
    collected = collect_all_points(world_plates)
    assert collected is not None
    points, elevation, owner = collected
    total = sum(p.node_count() for p in world_plates)
    assert len(points) == total
    assert len(elevation) == total
    assert len(owner) == total
    assert set(owner.tolist()) <= {p.plate_id for p in world_plates}


def test_collect_all_points_none_when_every_plate_is_empty():
    world_plates = generate_plates(seed=9, num_plates=3)
    for p in world_plates:
        p.set_lines([])
    assert collect_all_points(world_plates) is None


def test_nearest_plate_id_finds_the_owning_plate_at_its_own_seed():
    world_plates = generate_plates(seed=10, num_plates=8)
    for p in world_plates:
        if p.node_count() == 0:
            continue
        assert nearest_plate_id(world_plates, p.seed_world) == p.plate_id


def test_collect_all_soil_and_resource_fields_are_index_aligned_with_collect_all_points():
    world_plates = generate_plates(seed=9, num_plates=6)
    for p in world_plates:
        for i, line in enumerate(p.lines):
            n = len(line.theta)
            if n == 0:
                continue
            p.replace_line(i, line.replace(soil_depth=np.full(n, 2.5), coal_deposit_m=np.full(n, 1.5)))
    points, _, _ = collect_all_points(world_plates)
    soil_depth = collect_all_soil_depth(world_plates)
    coal = collect_all_coal_deposit(world_plates)
    mineral = collect_all_mineral_deposit(world_plates)
    oil_gas = collect_all_oil_gas_deposit(world_plates)
    soil_mineral = collect_all_soil_mineral_content(world_plates)
    soil_organic = collect_all_soil_organic_content(world_plates)
    for arr in (soil_depth, coal, mineral, oil_gas, soil_mineral, soil_organic):
        assert arr.shape == (len(points),)
    assert np.all(soil_depth == 2.5)
    assert np.all(coal == 1.5)
    assert np.all(mineral == 0.0)

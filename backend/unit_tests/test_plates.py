import numpy as np
import pytest
from app import geometry
from app.elevation_lines import ElevationLine, line_spacing_rad
from app.lithosphere_plate import build_plate_tiling, generate_plates
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
    nearest_plate_id,
    node_components,
    plate_bounding_ellipse,
)
from app.world import generate_world, step_world


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


def test_every_node_is_closest_to_a_site_of_its_own_plate():
    """A node kept by a plate must fall in one of that plate's own merged Voronoi cells --
    i.e. its angularly-nearest *site* is owned by that plate, no other plate's site is
    closer. Reconstructs the tiling the same way generate_plates does: default_rng(seed),
    then build_plate_tiling as its first draw (num_plates given, continental_fraction None)."""
    seed, num_plates = 3, 8
    tiling = build_plate_tiling(np.random.default_rng(seed), num_plates)
    plates = generate_plates(seed=seed, num_plates=num_plates)

    for p in plates:
        for line in p.lines:
            world_pts = line.world_xyz(p.frame)
            dists = geometry.angular_distance(world_pts[:, None, :], tiling.site_xyz[None, :, :])
            nearest_site = np.argmin(dists, axis=1)
            assert np.all(tiling.site_plate[nearest_site] == p.plate_id)


def test_build_plate_tiling_is_deterministic_and_covers_every_plate():
    a = build_plate_tiling(np.random.default_rng(11), num_plates=7)
    b = build_plate_tiling(np.random.default_rng(11), num_plates=7)
    assert np.array_equal(a.site_xyz, b.site_xyz)
    assert np.array_equal(a.site_plate, b.site_plate)
    # every plate owns at least its own primary cell
    assert set(a.site_plate.tolist()) == set(range(7))
    assert np.array_equal(a.site_plate[:7], np.arange(7))


def test_build_plate_tiling_extra_sites_zero_recovers_one_cell_per_plate():
    tiling = build_plate_tiling(np.random.default_rng(1), num_plates=6, extra_sites_per_plate=0)
    assert len(tiling.site_xyz) == 6
    assert np.array_equal(tiling.site_plate, np.arange(6))


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
    # 70% continental plates leaves comfortably more continental area than 29% land needs.
    # _land_noise_threshold corrects for the isostatic sea-level offset of the reference
    # continental column (see its `sealevel_noise_offset` param), so the measured land
    # fraction tracks the request closely rather than overshooting it.
    plates = generate_plates(seed=2, num_plates=14, continental_fraction=0.7, land_fraction=0.29)
    assert abs(_measured_land_fraction(plates) - 0.29) < 0.1


def test_generate_plates_land_fraction_is_capped_by_continental_area():
    # Only ~1/4 of plates (by count, roughly by area too) are continental, so 80% land is
    # not achievable -- every continental node should end up as land (elevation > 0) and no
    # more, capping measured land at roughly the continental area fraction itself.
    plates = generate_plates(seed=2, num_plates=14, continental_fraction=0.25, land_fraction=0.8)
    continental_area = _measured_continental_area_fraction(plates)
    assert abs(_measured_land_fraction(plates) - continental_area) < 0.02


def test_generate_plates_land_fraction_zero_gives_no_land():
    # Not guaranteed to be an exact 0.0: the threshold is estimated from a coarser
    # whole-sphere sample (LAND_FRACTION_SAMPLE_SPACING_KM) than the actual plate lattice it
    # is applied to, so a handful of real nodes can sit fractionally above that sample's max.
    plates = generate_plates(seed=2, num_plates=14, continental_fraction=0.7, land_fraction=0.0)
    assert _measured_land_fraction(plates) < 0.02


def _land_points_and_elevation(plates_list):
    from app.plates import gather_node_positions

    points, ordered = gather_node_positions(plates_list)
    elevation = np.concatenate([p.collect("elevation") for p in ordered])
    land = elevation > 0.0
    return points[land], elevation[land]


def test_generation_elevation_is_deterministic_for_same_seed():
    # Stronger than test_generation_is_deterministic_for_same_seed (which only checks
    # frame/crust type/line count): the composite relief field must reproduce the exact
    # per-node elevation for a given seed.
    e1 = np.concatenate([l.elevation for p in generate_plates(seed=321, num_plates=10) for l in p.lines])
    e2 = np.concatenate([l.elevation for p in generate_plates(seed=321, num_plates=10) for l in p.lines])
    assert np.array_equal(e1, e2)


def test_generation_relief_has_more_variety_than_a_smooth_base():
    # The orogenic/plateau uplift materially widens the spread of land elevation versus a
    # world seeded from the low-frequency sample() field alone.
    import app.terrain_noise as terrain_noise

    kw = dict(seed=6, num_plates=12, continental_fraction=0.55, land_fraction=0.29, node_density=1.0)
    _, varied = _land_points_and_elevation(generate_plates(**kw))

    original = terrain_noise.ContinentalRelief.uplift
    terrain_noise.ContinentalRelief.uplift = lambda self, xyz: np.zeros(np.shape(xyz)[:-1])
    try:
        _, smooth = _land_points_and_elevation(generate_plates(**kw))
    finally:
        terrain_noise.ContinentalRelief.uplift = original

    assert np.std(varied) > 1.5 * np.std(smooth)
    assert np.std(varied) > 900.0


def test_generation_has_clustered_mountain_ranges():
    # Nodes above 3 km are not scattered singletons -- they form connected belts.
    lp, le = _land_points_and_elevation(
        generate_plates(seed=6, num_plates=12, continental_fraction=0.55, land_fraction=0.29, node_density=1.0)
    )
    peaks = le > 3000.0
    assert peaks.sum() > 100
    labels = node_components(lp[peaks], 2.2 * line_spacing_rad(1.0))
    assert np.max(np.bincount(labels)) >= 50  # one contiguous range of >=50 peak nodes


def test_generation_has_elevated_flats():
    # A plateau reads as high ground that is also locally flat -- distinct from a peak,
    # which is high but locally rough.
    from scipy.spatial import cKDTree

    lp, le = _land_points_and_elevation(
        generate_plates(seed=3, num_plates=12, continental_fraction=0.55, land_fraction=0.29, node_density=1.0)
    )
    _, idx = cKDTree(lp).query(lp, k=10)
    local_std = le[idx].std(axis=1)
    elevated_flat = (le > 1600.0) & (local_std < 300.0)
    assert elevated_flat.sum() > 150


def test_outline_world_traces_a_loop_covering_every_line():
    # A staircase, not a smooth scanline (see PlateWithLines.outline_world's own docstring):
    # every line contributes at least one point per side (high/low theta), plus extra
    # "corner" points wherever two adjacent rows' theta extents actually differ -- which is
    # the normal case for a real, non-uniform plate, not the exception -- so the loop is at
    # least, not exactly, 2 points per line.
    plates = generate_plates(seed=5, num_plates=8)
    for p in plates:
        lines_with_nodes = [line for line in p.lines if len(line.theta) > 0]
        outline = p.outline_world()
        assert len(outline) >= 2 * len(lines_with_nodes)
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


def test_get_node_kdtree_is_cached_until_geometry_changes():
    p = generate_plates(seed=12, num_plates=8)[0]
    first = p.get_node_kdtree()
    assert first is p.get_node_kdtree()  # same cached tree
    p.rotate(geometry.plate_frame_from_seed(np.array([0.0, 1.0, 0.0])))
    rotated = p.get_node_kdtree()
    assert rotated is not first  # invalidated by rotate
    assert np.allclose(np.asarray(rotated.data), p.all_points_and_elevation()[0])
    p.set_lines(list(p.lines))
    assert p.get_node_kdtree() is not rotated  # invalidated by set_lines


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


def _sampled_overlap_fraction(plates_list, sample_per_plate: int = 20) -> float:
    """Fraction of sampled nodes (each plate's own nodes, thinned to at most
    `sample_per_plate`) found geometrically inside a *different* plate's current
    `get_bounding_polygon()` -- the closest testable proxy for "bounding polygons don't
    overlap" available with an envelope-based (not exact) outline. Not expected to be
    exactly zero: `PlateWithLines.deform`'s own docstring, and docs/simulation-model.md's
    account of this design's known limits, explain why -- one-turn processing lag (a
    neighbour that grows into adjacent space later in the same turn's randomized order
    isn't re-checked until this plate's own next turn) and residual envelope looseness for
    a genuinely non-convex, lateral-sheared plate shape. What *should* hold is that this
    stays bounded rather than climbing without limit turn over turn -- see
    stress_tests/test_world_stepping.py's own long-running version of this same check."""
    total = 0
    overlapping = 0
    for plate in plates_list:
        points, _ = plate.all_points_and_elevation()
        if len(points) == 0:
            continue
        sample = points[:: max(1, len(points) // sample_per_plate)][:sample_per_plate]
        for other in plates_list:
            if other.plate_id == plate.plate_id:
                continue
            polygon = other.get_bounding_polygon()
            if len(polygon) < 3:
                continue
            total += len(sample)
            overlapping += int(np.count_nonzero(geometry.points_in_spherical_polygon(sample, polygon)))
    return overlapping / total if total > 0 else 0.0


def test_deform_keeps_plate_overlap_bounded_not_runaway():
    world = generate_world(seed=3, num_plates=8, node_density=0.5)
    world.simulate_climate_biomes = False  # only plate geometry is checked here
    for _ in range(3):
        step_world(world, years=3_000_000)
    early = _sampled_overlap_fraction(world.plates)
    for _ in range(10):
        step_world(world, years=3_000_000)
    late = _sampled_overlap_fraction(world.plates)

    # Generous ceiling: this is a smoke check against a severe regression (e.g. a shrink/
    # grow bug letting contested territory balloon unchecked), not a precision bound --
    # confirmed empirically to sit in the 10-20% range for this seed, stable across many
    # steps, not climbing toward saturation.
    assert late < 0.35
    # And it shouldn't have grown much further from where it started -- a real runaway
    # would keep climbing step over step, not plateau.
    assert late < early + 0.15


# -- node_components / Plate.defragment / _plates_from_node_masks -----------------------
#
# Geometric plate cleanup -- see merge_split.defragment_plates and Plate.defragment.
# Ordinary deform() only shrinks a line from its ends and never deletes its last node, so
# subduction/transform can sever one Plate's node cloud into two disconnected landmasses or
# strand a comb of one-node rows; maybe_split_plate only cuts on mantle-flow disagreement,
# not geometry, so it never catches either. These exercise the pass that does.

_DEFRAG_SPACING_RAD = line_spacing_rad(1.0)
_DEFRAG_CONNECT_RAD = 2.5 * _DEFRAG_SPACING_RAD


def _lobed_plate(lobes, plate_id=0, rows=12, per_row=8, **plate_kwargs):
    """A `PlateWithLines` (frame = identity, so plate-local phi/theta are world lat/lon)
    whose nodes form one connected blob per entry in `lobes`. Each entry is either a theta
    centre (radians) or a `(centre, nodes_per_row)` tuple; centres must sit far enough apart
    to read as separate connected components at `_DEFRAG_CONNECT_RAD`. `rows` lines are
    stacked one spacing apart in phi, so a lobe contributes `rows * nodes_per_row` nodes."""
    lines = []
    for r in range(rows):
        chunks = []
        for lobe in lobes:
            centre, count = lobe if isinstance(lobe, tuple) else (lobe, per_row)
            chunks.append(centre + np.arange(count) * _DEFRAG_SPACING_RAD)
        theta = np.concatenate(chunks)
        lines.append(ElevationLine(phi=r * _DEFRAG_SPACING_RAD, theta=theta, elevation=np.zeros(len(theta))))
    return PlateWithLines(plate_id=plate_id, frame=np.eye(3), crust_type="oceanic", lines=lines, **plate_kwargs)


def test_node_components_labels_isolated_clusters_separately():
    points, _ = _lobed_plate([0.0, 0.6]).all_points_and_elevation()
    labels = node_components(points, _DEFRAG_CONNECT_RAD)
    _, counts = np.unique(labels, return_counts=True)
    assert sorted(counts.tolist()) == [96, 96]  # two equal lobes, 12 rows x 8 nodes each


def test_node_components_one_label_for_a_contiguous_blob():
    points, _ = _lobed_plate([0.0]).all_points_and_elevation()
    assert set(node_components(points, _DEFRAG_CONNECT_RAD).tolist()) == {0}


def test_node_components_empty_input():
    assert node_components(np.zeros((0, 3)), 0.1).shape == (0,)


def test_defragment_splits_a_severed_plate_and_keeps_identity_on_the_largest():
    plate = _lobed_plate([(0.0, 10), (0.6, 6)], plate_id=7, omega=np.array([0.1, 0.2, 0.3]), age_steps=9)
    before = plate.node_count()

    result = plate.defragment(next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50)
    assert result is not None
    replacements, consumed = result

    assert consumed == 1
    assert [p.plate_id for p in replacements] == [7, 20]
    assert sum(p.node_count() for p in replacements) == before
    assert replacements[0].node_count() > replacements[1].node_count()  # largest keeps the id
    # The identity-keeper carries this plate's own omega and age; the fresh fragment shares
    # the omega (it was co-moving, which is why nothing split it off) but resets age to 0.
    assert np.allclose(replacements[0].omega, [0.1, 0.2, 0.3])
    assert replacements[0].age_steps == 9
    assert np.allclose(replacements[1].omega, [0.1, 0.2, 0.3])
    assert replacements[1].age_steps == 0
    for p in replacements:
        pts, _ = p.all_points_and_elevation()
        assert len(np.unique(node_components(pts, _DEFRAG_CONNECT_RAD))) == 1


def test_defragment_sheds_stranded_nodes_without_splitting():
    # second lobe is 12 nodes (1 per row), well below min_fragment_nodes -- dropped, not
    # promoted to its own plate, and no new id is consumed.
    plate = _lobed_plate([(0.0, 10), (0.6, 1)], plate_id=3)
    before = plate.node_count()

    result = plate.defragment(next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50)
    assert result is not None
    replacements, consumed = result

    assert consumed == 0
    assert [p.plate_id for p in replacements] == [3]
    assert replacements[0].node_count() == before - 12


def test_defragment_leaves_a_contiguous_plate_alone():
    plate = _lobed_plate([0.0])
    assert plate.defragment(next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50) is None


def test_defragment_leaves_an_all_debris_plate_for_the_territory_check():
    # three lobes, none reaching min_fragment_nodes: defrag declines (returns None) rather
    # than deleting a whole plate itself -- has_negligible_territory / remove_defunct_plates
    # own that call.
    plate = _lobed_plate([(0.0, 2), (0.6, 2), (1.2, 2)], rows=10)
    assert plate.defragment(next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50) is None


def test_defragment_partition_carries_each_nodes_own_fields_to_the_right_fragment():
    plate = _lobed_plate([0.0, 0.6], plate_id=4)
    points, _ = plate.all_points_and_elevation()
    marker = np.arange(len(points), dtype=float)  # a distinct value per node
    offset = 0
    for i, line in enumerate(plate.lines):
        k = len(line)
        plate.replace_line(i, line.replace(channel_depth=marker[offset : offset + k]))
        offset += k

    replacements, _ = plate.defragment(
        next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50
    )
    recombined = np.concatenate([p.collect("channel_depth") for p in replacements])
    assert sorted(recombined.tolist()) == sorted(marker.tolist())


def test_has_negligible_territory_flags_a_comb_of_one_node_stubs():
    # Many lines but ~1 node each: deform() decayed a heavily-subducted plate into stranded
    # rows. High line count masks that there's no 2D patch left -- the original
    # len(lines) <= 1 test missed this.
    comb = PlateWithLines(
        plate_id=0,
        frame=np.eye(3),
        crust_type="oceanic",
        lines=[ElevationLine(phi=i * _DEFRAG_SPACING_RAD, theta=np.array([0.0]), elevation=np.zeros(1)) for i in range(40)],
    )
    assert comb.has_negligible_territory()


def test_has_negligible_territory_false_for_a_plate_with_real_rows():
    assert not _lobed_plate([0.0]).has_negligible_territory()


# -- pole cap / theta-winding guard ---------------------------------------------------
#
# A row is a circle of local latitude, so its theta extent can't physically exceed a full
# 2*pi revolution. Nothing treats theta as periodic, so before the guard in
# _grow_or_shrink_line_for_deform a plate that grew to encircle its own local pole -- where
# the "gap to nearest neighbour" is wide open forever, the cap belonging to nobody -- just
# kept winding the same ring, covering the same ground dozens of times (concentric-circle /
# hole artifacts in the Plate Inspector, unbounded overlap + node count on long runs).


def test_deform_never_winds_a_row_past_a_full_revolution():
    from app.world import World

    spacing = line_spacing_rad(1.0)
    # An isolated plate hugging its own local pole (frame = identity, so local phi is world
    # latitude): near-pole rows plus a mid-latitude row for contrast. No neighbours, so every
    # end is "wide open" every step.
    near_pole_phis = [np.pi / 2 - k * spacing for k in (8, 7, 6, 5, 4)]
    lines = [
        ElevationLine(phi=phi, theta=np.linspace(-0.5, 0.5, 6), elevation=np.zeros(6))
        for phi in [0.3] + near_pole_phis
    ]
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="oceanic", lines=lines)
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)

    for _ in range(25):
        plate.deform(world, [], years=1_000_000, max_distance=5 * spacing)

    for line in plate.lines:
        span = float(line.theta[-1] - line.theta[0])
        assert span <= 2.0 * np.pi + spacing, f"row at phi={line.phi:.3f} wound to {span:.2f} rad"
        assert np.all(np.diff(line.theta) > 0)  # still a sorted, contiguous span

    # A near-pole row *did* close its loop (the guard caps it, it isn't just slow growth).
    assert any(
        line.theta[-1] - line.theta[0] > 2.0 * np.pi - 3 * (spacing / max(np.cos(line.phi), 1e-3))
        for line in plate.lines
        if line.phi > 1.0
    )
    # ... and the mid-latitude row is still an ordinary partial arc, untouched by the cap.
    mid = min(plate.lines, key=lambda ln: ln.phi)
    assert mid.theta[-1] - mid.theta[0] < 2.0 * np.pi


def test_claim_adjacent_territory_keeps_a_margin_from_the_local_pole():
    from app.plates import POLE_CAP_MARGIN_MULT
    from app.world import World

    spacing = line_spacing_rad(1.0)
    # A plate whose poleward rows already sit just inside the pole-cap margin: claim should
    # refuse to add a new row past it, however open the space is.
    phis = [np.pi / 2 - k * spacing for k in (POLE_CAP_MARGIN_MULT + 2, POLE_CAP_MARGIN_MULT + 1)]
    lines = [
        ElevationLine(phi=phi, theta=np.linspace(-0.3, 0.3, 20), elevation=np.zeros(20))
        for phi in phis
    ]
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="oceanic", lines=lines)
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)

    for _ in range(10):
        plate.deform(world, [], years=1_000_000, max_distance=5 * spacing)

    assert max(ln.phi for ln in plate.lines) <= np.pi / 2 - POLE_CAP_MARGIN_MULT * spacing + 1e-9


# The v2 engine (LithospherePlate) overrides both _grow_or_shrink_line_for_deform and
# _claim_adjacent_territory, and originally shipped without either pole-winding guard -- it
# relied entirely on regularize_line's after-the-fact unwind, so rows still over-wound by up
# to a revolution every step. These mirror the two PlateWithLines tests above for the class
# that actually runs.


def _lithosphere_polar_plate(near_pole_phis, theta, crust_type="oceanic"):
    from app.lithosphere import reference_thickness
    from app.lithosphere_plate import LithospherePlate

    hc0, hm0 = reference_thickness(crust_type)
    lines = [
        ElevationLine(
            phi=phi,
            theta=theta.copy(),
            elevation=np.zeros(len(theta)),
            crustal_thickness_m=np.full(len(theta), hc0),
            mantle_lithosphere_thickness_m=np.full(len(theta), hm0),
        )
        for phi in near_pole_phis
    ]
    return LithospherePlate(plate_id=0, frame=np.eye(3), crust_type=crust_type, lines=lines)


def test_lithosphere_deform_never_winds_a_row_past_a_full_revolution():
    from app.world import World

    spacing = line_spacing_rad(1.0)
    near_pole_phis = [0.3] + [np.pi / 2 - k * spacing for k in (8, 7, 6, 5, 4)]
    plate = _lithosphere_polar_plate(near_pole_phis, np.linspace(-0.5, 0.5, 6))
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)

    from app.elevation_lines import needs_regularizing

    for _ in range(25):
        plate.deform(world, [], years=1_000_000, max_distance=5 * spacing)

    grew_a_full_ring = False
    for line in plate.lines:
        span = float(line.theta[-1] - line.theta[0])
        dtheta = spacing / max(np.cos(line.phi), 1e-3)
        assert span <= 2.0 * np.pi + spacing, f"row at phi={line.phi:.3f} wound to {span:.2f} rad"
        assert np.all(np.diff(line.theta) > 0)
        if line.phi > 1.0 and span > 2.0 * np.pi - 4 * dtheta:
            grew_a_full_ring = True

    # A near-pole row *did* grow all the way around to the ring cap (the cap stopped it there,
    # it isn't just slow growth).
    assert grew_a_full_ring
    # ... and the mid-latitude row is still an ordinary partial arc, untouched by the cap.
    mid = min(plate.lines, key=lambda ln: ln.phi)
    assert mid.theta[-1] - mid.theta[0] < 2.0 * np.pi

    # The point of the ring_room cap (vs. leaning on regularize_line's after-the-fact unwind):
    # once a near-pole ring has closed, end-growth stops there rather than over-winding past
    # 2*pi and being unwound again on the very next step, step after step. Spy on one more
    # step: no row should still be tripping regularize_line (pre-fix, the near-pole rings wind
    # past 2*pi and get unwound every single step, forever).
    assert not any(needs_regularizing(line, spacing) for line in plate.lines)
    import app.lithosphere_plate as _lp

    regularized_phis = []
    orig = _lp.regularize_line

    def _spy(line, *a, **k):
        regularized_phis.append(round(line.phi, 3))
        return orig(line, *a, **k)

    _lp.regularize_line = _spy
    try:
        plate.deform(world, [], years=1_000_000, max_distance=5 * spacing)
    finally:
        _lp.regularize_line = orig
    assert not regularized_phis, f"rows still churning through regularize every step: {regularized_phis}"


def test_runs_of_at_least_clears_short_true_runs():
    from app.lithosphere_plate import _runs_of_at_least

    mask = np.array([1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1], dtype=bool)
    assert list(_runs_of_at_least(mask, 3)) == [0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1]
    assert list(_runs_of_at_least(mask, 1)) == list(mask)
    assert not _runs_of_at_least(np.array([1, 1, 0, 0, 1], dtype=bool), 3).any()
    # a qualifying run flush against the end is kept
    assert list(_runs_of_at_least(np.array([0, 0, 1, 1, 1], dtype=bool), 3)) == [0, 0, 1, 1, 1]


def test_lithosphere_continental_contested_edge_retreats():
    """A continental line's contested end retreats one node per step whether the overriding
    neighbour is oceanic (a passive margin / accretion front, breaking the node ratchet) or
    continental (a suture whose territory overlap is consumed into the orogen rather than
    frozen for tens of Myr). Both cases lose real theta extent."""
    from app.lithosphere_plate import CONTINENTAL_CONTESTED_RETREAT_MIN_RUN, LithospherePlate
    from app.lithosphere import reference_thickness
    from app.world import World

    spacing = line_spacing_rad(1.0)

    def _plate(pid, crust_type, theta_lo, theta_hi, n):
        hc0, hm0 = reference_thickness(crust_type)
        theta = np.linspace(theta_lo, theta_hi, n)
        line = ElevationLine(
            phi=0.2,
            theta=theta,
            elevation=np.zeros(n),
            crustal_thickness_m=np.full(n, hc0),
            mantle_lithosphere_thickness_m=np.full(n, hm0),
        )
        filler = ElevationLine(
            phi=-0.6,
            theta=np.linspace(-0.2, 0.2, 8),
            elevation=np.zeros(8),
            crustal_thickness_m=np.full(8, hc0),
            mantle_lithosphere_thickness_m=np.full(8, hm0),
        )
        return LithospherePlate(plate_id=pid, frame=np.eye(3), crust_type=crust_type, lines=[line, filler])

    def _high_end_retreat(neighbour_crust: str) -> float:
        continent = _plate(0, "continental", -0.5, 0.5, 40)
        # Same frame, theta range overlapping the continent's high end -> the continent's
        # nodes with theta > ~0.15 fall inside this neighbour's polygon (contested). The
        # regularize pass keeps the node *count* ~constant, so it's the line's theta *extent*
        # (its actual territory) that retreats.
        neighbour = _plate(1, neighbour_crust, 0.15, 0.9, 40)
        world = World(seed=0, plates=[continent, neighbour], mantle_centers=[], node_density=1.0)

        def high_theta() -> float:
            return max(ln.theta[-1] for ln in continent.lines if abs(ln.phi - 0.2) < 1e-6)

        before = high_theta()
        for _ in range(6):
            continent.deform(world, [neighbour], years=200_000, max_distance=1.5 * spacing)
        return before - high_theta()

    retreat_against_ocean = _high_end_retreat("oceanic")
    retreat_against_continent = _high_end_retreat("continental")
    assert retreat_against_ocean > CONTINENTAL_CONTESTED_RETREAT_MIN_RUN * spacing
    assert retreat_against_continent > CONTINENTAL_CONTESTED_RETREAT_MIN_RUN * spacing


def test_redistribute_accreted_column_conserves_crustal_volume():
    """`_redistribute_accreted_column` moves the exact summed Hc/Hm of the dropped, accretion-
    flagged nodes onto the surviving edge nodes (node area is constant, so summed thickness is
    the conserved volume), and lifts their elevation by the matching isostatic delta. Dropped
    nodes not flagged (a passive margin against an oceanic slab) contribute nothing."""
    from app.lithosphere import crust_density, isostatic_elevation
    from app.lithosphere_plate import _redistribute_accreted_column, SUTURE_ACCRETION_SPREAD_NODES

    rho_c = crust_density("continental")
    hc = np.full(10, 35_000.0)
    hm = np.full(10, 100_000.0)
    fields = {"crustal_thickness_m": hc, "mantle_lithosphere_thickness_m": hm}
    elevation = isostatic_elevation(hc, hm, rho_c).copy()

    removed_hc = np.array([35_000.0, 35_000.0, 35_000.0])
    accrete_removed = np.array([True, True, False])  # last one was against ocean -> subducts

    total_hc_before = hc.sum()
    _redistribute_accreted_column(fields, elevation, rho_c, removed_hc, accrete_removed, from_high=True)

    assert fields["crustal_thickness_m"].sum() == pytest.approx(total_hc_before + 2 * 35_000.0)
    # spread over the last SUTURE_ACCRETION_SPREAD_NODES nodes, evenly
    k = SUTURE_ACCRETION_SPREAD_NODES
    assert np.allclose(fields["crustal_thickness_m"][-k:], 35_000.0 + 2 * 35_000.0 / k)
    assert np.allclose(fields["crustal_thickness_m"][:-k], 35_000.0)
    # thicker crust -> higher ground on exactly those nodes
    assert np.all(elevation[-k:] > elevation[:-k].max())

    # A suture that never heals would pile columns onto the same nodes forever -- Hc is capped
    # (the overflow delaminates), so it can't run away.
    from app.lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M

    hc2 = np.full(6, 60_000.0)
    fields2 = {"crustal_thickness_m": hc2, "mantle_lithosphere_thickness_m": np.full(6, 100_000.0)}
    elev2 = isostatic_elevation(hc2, fields2["mantle_lithosphere_thickness_m"], rho_c).copy()
    _redistribute_accreted_column(
        fields2, elev2, rho_c,
        np.full(5, 90_000.0), np.ones(5, dtype=bool), from_high=True,
    )
    assert np.all(fields2["crustal_thickness_m"] <= SUTURE_ACCRETION_MAX_HC_M + 1e-6)


def test_continent_continent_suture_consumes_its_overlap_as_mass_conserving_accretion():
    """The overlap a retreating continent-continent suture consumes is thrust onto the plate's
    own leading edge, not discarded: the retreated column's crustal volume reappears on the
    surviving edge nodes (so the belt thickens and the plate's total Hc is conserved). A
    retreat against an *oceanic* neighbour is a passive margin -- that column subducts, so the
    edge does not thicken and total Hc drops."""
    from app.lithosphere import reference_thickness
    from app.lithosphere_plate import LithospherePlate
    from app.world import World

    hc0, hm0 = reference_thickness("continental")
    spacing = line_spacing_rad(1.0)

    def _plate(pid, crust_type, theta_lo, theta_hi, n):
        theta = np.linspace(theta_lo, theta_hi, n)
        line = ElevationLine(
            phi=0.2,
            theta=theta,
            elevation=np.zeros(n),
            crustal_thickness_m=np.full(n, hc0 if crust_type == "continental" else reference_thickness("oceanic")[0]),
            mantle_lithosphere_thickness_m=np.full(n, hm0 if crust_type == "continental" else reference_thickness("oceanic")[1]),
        )
        filler = ElevationLine(
            phi=-0.6,
            theta=np.linspace(-0.2, 0.2, 8),
            elevation=np.zeros(8),
            crustal_thickness_m=np.full(8, hc0),
            mantle_lithosphere_thickness_m=np.full(8, hm0),
        )
        return LithospherePlate(plate_id=pid, frame=np.eye(3), crust_type=crust_type, lines=[line, filler])

    def _run(neighbour_crust: str) -> tuple[float, float]:
        continent = _plate(0, "continental", -0.5, 0.5, 40)
        neighbour = _plate(1, neighbour_crust, 0.15, 0.9, 40)  # overlaps the continent's high end
        world = World(seed=0, plates=[continent, neighbour], mantle_centers=[], node_density=1.0)

        def suture_line():
            return next(ln for ln in continent.lines if abs(ln.phi - 0.2) < 1e-6)

        edge_hc_before = float(suture_line().crustal_thickness_m[-4:].max())
        for _ in range(6):
            continent.deform(world, [neighbour], years=200_000, max_distance=1.5 * spacing)
        edge_hc_after = float(suture_line().crustal_thickness_m[-4:].max())
        return edge_hc_after - edge_hc_before, float(suture_line().crustal_thickness_m.sum())

    edge_gain_cc, suture_hc_cc = _run("continental")
    edge_gain_co, suture_hc_co = _run("oceanic")

    # Suture accretion piles the consumed columns onto the leading edge -- kilometres, not the
    # metres a bare sub-yield graze would add.
    assert edge_gain_cc > 3_000.0
    # Passive margin against an oceanic slab: the retreated column is subducted, so no edge
    # pile-up...
    assert edge_gain_co < edge_gain_cc / 3
    # ...and the suture line keeps several consumed continental columns' worth of extra crust
    # that the oceanic-neighbour run simply loses.
    assert suture_hc_cc > suture_hc_co + 3 * 30_000.0


def test_lithosphere_continental_volume_budget_suppresses_growth():
    """A continental plate whose node footprint has outrun its crustal volume -- most of its
    lattice diluted to the oceanic reference column by the boundary ratchet -- grows no new
    areal crust this step (neither end-growth nor a claimed new row), so it thins/drowns back
    toward budget instead of tiling drowned margin outward forever. A plate at genuine
    continental thickness everywhere is within budget and still grows normally."""
    from app.lithosphere import (
        REFERENCE_HC_CONTINENTAL_M,
        REFERENCE_HC_OCEANIC_M,
        REFERENCE_HM_CONTINENTAL_M,
    )
    from app.lithosphere_plate import LithospherePlate
    from app.world import World

    spacing = line_spacing_rad(1.0)

    def _continent(main_hc: float) -> LithospherePlate:
        n = 40
        main = ElevationLine(
            phi=0.2,
            theta=np.linspace(-0.5, 0.5, n),
            elevation=np.zeros(n),
            crustal_thickness_m=np.full(n, main_hc),
            mantle_lithosphere_thickness_m=np.full(n, REFERENCE_HM_CONTINENTAL_M),
        )
        # A genuinely-continental core so `n_continental` is never zero in either case.
        core = ElevationLine(
            phi=-0.6,
            theta=np.linspace(-0.2, 0.2, 8),
            elevation=np.zeros(8),
            crustal_thickness_m=np.full(8, REFERENCE_HC_CONTINENTAL_M),
            mantle_lithosphere_thickness_m=np.full(8, REFERENCE_HM_CONTINENTAL_M),
        )
        return LithospherePlate(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[main, core])

    def _main_span_growth(main_hc: float) -> float:
        plate = _continent(main_hc)
        world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)

        def main_span() -> float:
            line = next(ln for ln in plate.lines if abs(ln.phi - 0.2) < 1e-6)
            return float(line.theta[-1] - line.theta[0])

        before = main_span()
        # One step, with a wide max_distance so an unsuppressed end can take several nodes --
        # a later step would see `_claim_adjacent_territory` dilute the compact plate too (the
        # ratchet this gate exists to bound), so measure the first step alone.
        plate.deform(world, [], years=200_000, max_distance=5 * spacing)
        return main_span() - before

    # Diluted lattice (main row at the oceanic reference column): ~40 nodes vs ~8 genuine
    # continental -> well past CONTINENTAL_AREA_BUDGET_MULT -> the open end grows nothing.
    assert _main_span_growth(REFERENCE_HC_OCEANIC_M) < 0.5 * spacing
    # Genuine continental thickness everywhere -> within budget -> the open end still grows.
    assert _main_span_growth(REFERENCE_HC_CONTINENTAL_M) > 3 * spacing


def test_lithosphere_active_margin_grows_arc_crust_not_ocean_floor():
    """A continental plate's *leading* edge advancing into space a subducting oceanic slab is
    vacating grows juvenile arc / accreted-terrane crust (the thicker ARC_MARGIN_SEED_*
    column, stamped as a subduction arc), not the drowned oceanic reference column
    `growth_seed_thickness` seeds everywhere else. The active-margin signal here is the
    subduction-arc provenance stamp a recent convergent step left on the leading nodes."""
    from app.elevation_lines import ELEV_CHANGE_NEW_CRUST, ELEV_CHANGE_SUBDUCTION_ARC
    from app.lithosphere import REFERENCE_HC_CONTINENTAL_M, REFERENCE_HC_OCEANIC_M, REFERENCE_HM_CONTINENTAL_M
    from app.lithosphere_plate import ARC_MARGIN_SEED_HC_M, LithospherePlate
    from app.world import World

    spacing = line_spacing_rad(1.0)

    def _grown_high_end(mark_active_margin: bool):
        n = 40
        reason = np.zeros(n)
        if mark_active_margin:
            reason[-4:] = ELEV_CHANGE_SUBDUCTION_ARC
        main = ElevationLine(
            phi=0.2,
            theta=np.linspace(-0.5, 0.5, n),
            elevation=np.zeros(n),
            crustal_thickness_m=np.full(n, REFERENCE_HC_CONTINENTAL_M),
            mantle_lithosphere_thickness_m=np.full(n, REFERENCE_HM_CONTINENTAL_M),
            elev_change_reason=reason,
        )
        plate = LithospherePlate(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[main])
        world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)
        plate.deform(world, [], years=200_000, max_distance=5 * spacing)
        line = next(ln for ln in plate.lines if abs(ln.phi - 0.2) < 1e-6)
        new = line.theta > 0.5 + 1e-9
        assert new.any(), "expected the open high end to grow"
        return line.crustal_thickness_m[new], line.elev_change_reason[new]

    ocean_hc, ocean_reason = _grown_high_end(mark_active_margin=False)
    arc_hc, arc_reason = _grown_high_end(mark_active_margin=True)

    # Baseline: new margin nodes seeded at (or near, after regularize blends toward the old
    # continental interior) the oceanic reference column, stamped plain new crust.
    assert ocean_hc.mean() < 0.5 * (REFERENCE_HC_OCEANIC_M + REFERENCE_HC_CONTINENTAL_M)
    assert (ocean_reason == ELEV_CHANGE_NEW_CRUST).any()
    # Active margin: markedly thicker juvenile crust, stamped as a subduction arc.
    assert arc_hc.mean() > ocean_hc.mean() + 5_000.0
    assert arc_hc.max() >= ARC_MARGIN_SEED_HC_M - 1e-6 or arc_hc.mean() > 20_000.0
    assert (arc_reason == ELEV_CHANGE_SUBDUCTION_ARC).any()


def test_lithosphere_arc_magmatism_thickens_the_continental_margin_band(monkeypatch):
    """A converging oceanic neighbour underplates juvenile crust across the overriding
    continental plate's arc *band* (out to `reach_rad`, not just the contested contact line).
    Isolated by running the same setup with the magmatic rate zeroed -- the difference is the
    arc contribution alone, and it is a real multi-hundred-metre Hc gain over a swath many
    nodes wide, not confined to the two or three nodes that actually overlap."""
    from app import mantle, rheology
    from app.lithosphere import REFERENCE_HC_CONTINENTAL_M, REFERENCE_HC_OCEANIC_M, REFERENCE_HM_CONTINENTAL_M, REFERENCE_HM_OCEANIC_M
    from app.lithosphere_plate import LithospherePlate
    from app.world import World

    spacing = line_spacing_rad(1.0)

    def _line(phi, lo, hi, count, hc, hm, elev):
        theta = np.linspace(lo, hi, count)
        return ElevationLine(
            phi=phi, theta=theta, elevation=np.full(count, elev),
            crustal_thickness_m=np.full(count, hc), mantle_lithosphere_thickness_m=np.full(count, hm),
        )

    def _run(arc_rate: float):
        monkeypatch.setattr(rheology, "ARC_MAGMATIC_HC_RATE_M_PER_MYR", arc_rate)
        cont = LithospherePlate(
            plate_id=0, frame=np.eye(3), crust_type="continental",
            lines=[
                _line(0.15, -0.6, 0.6, 60, REFERENCE_HC_CONTINENTAL_M, REFERENCE_HM_CONTINENTAL_M, 0.0),
                _line(-0.05, -0.6, 0.6, 60, REFERENCE_HC_CONTINENTAL_M, REFERENCE_HM_CONTINENTAL_M, 0.0),
            ],
        )
        ocean = LithospherePlate(
            plate_id=1, frame=np.eye(3), crust_type="oceanic",
            lines=[
                _line(0.15, 0.55, 1.6, 50, REFERENCE_HC_OCEANIC_M, REFERENCE_HM_OCEANIC_M, -4000.0),
                _line(-0.05, 0.55, 1.6, 50, REFERENCE_HC_OCEANIC_M, REFERENCE_HM_OCEANIC_M, -4000.0),
            ],
        )
        cont.set_omega(np.zeros(3))
        ocean.set_omega(np.array([0.0, 0.0, -0.35 * mantle.MAX_PLATE_RATE]))  # drifts toward the continent's high end
        world = World(seed=1, plates=[cont, ocean], mantle_centers=[], node_density=1.0)
        for _ in range(6):
            cont.deform(world, [ocean], years=500_000, max_distance=1.5 * spacing)
        line = next(ln for ln in cont.lines if abs(ln.phi - 0.15) < 1e-6)
        return line.theta.copy(), line.crustal_thickness_m.copy()

    default_rate = rheology.ARC_MAGMATIC_HC_RATE_M_PER_MYR
    theta_off, hc_off = _run(0.0)
    theta_on, hc_on = _run(default_rate)

    # Compare on the shared theta support (regularize can shift node counts a hair). The
    # difference is the arc contribution alone -- convergent shortening at the contested
    # contact runs identically in both.
    hi = theta_on > 0.2
    gain = hc_on[hi] - np.interp(theta_on[hi], theta_off, hc_off)
    assert gain.max() > 400.0  # the margin band genuinely thickened over the 3 My run
    assert int((gain > 30.0).sum()) >= 3  # ... over a swath, not a single contact node
    # No effect on the trailing (non-margin) half.
    lo = theta_on < -0.3
    assert np.allclose(hc_on[lo], np.interp(theta_on[lo], theta_off, hc_off), atol=5.0)


def test_lithosphere_contested_leading_row_is_dropped_after_sustained_override():
    """The parallel-suture retreat op: a continental plate's outermost phi-row that a
    neighbour has overridden over its full theta width -- no uncontested end for end-trim,
    no mid-row carve allowed -- is dropped whole once the override has held for a cumulative
    LEADING_ROW_RETREAT_SUSTAINED_YEARS. Inner rows and a control plate with open ground are
    untouched."""
    from app.lithosphere import reference_thickness
    from app.lithosphere_plate import (
        LEADING_ROW_CONTESTED_FRACTION,
        LEADING_ROW_RETREAT_SUSTAINED_YEARS,
        LithospherePlate,
    )
    from app.world import World

    spacing = line_spacing_rad(1.0)
    hc0, hm0 = reference_thickness("continental")

    def _rows(phis, theta_lo, theta_hi, n):
        theta = np.linspace(theta_lo, theta_hi, n)
        return [
            ElevationLine(
                phi=phi,
                theta=theta.copy(),
                elevation=np.zeros(n),
                crustal_thickness_m=np.full(n, hc0),
                mantle_lithosphere_thickness_m=np.full(n, hm0),
            )
            for phi in phis
        ]

    front_phi = 0.30
    inner_phis = [0.10, 0.15, 0.20, 0.25]
    continent = LithospherePlate(
        plate_id=0, frame=np.eye(3), crust_type="continental",
        lines=_rows(inner_phis + [front_phi], -0.5, 0.5, 40),
    )
    # A continental neighbour whose (densely-spaced) rows straddle the continent's front row
    # and just outrun it in theta -> every node of the phi=0.30 row falls inside the
    # neighbour's polygon (a suture *parallel* to the continent's own rows: no uncontested
    # end, so end-trim can't touch it), while the lower inner rows stay clear of it.
    neighbour = LithospherePlate(
        plate_id=1, frame=np.eye(3), crust_type="continental",
        lines=_rows(list(front_phi - 0.5 * spacing + spacing * np.arange(10)), -0.53, 0.53, 60),
    )
    world = World(seed=0, plates=[continent, neighbour], mantle_centers=[], node_density=1.0)

    def front_row_present() -> bool:
        return any(abs(ln.phi - front_phi) < 1e-6 for ln in continent.lines)

    def inner_rows_present() -> bool:
        return all(any(abs(ln.phi - p) < 1e-6 for ln in continent.lines) for p in inner_phis)

    years_per_step = 1_000_000
    steps_to_drop = int(np.ceil(LEADING_ROW_RETREAT_SUSTAINED_YEARS / years_per_step))

    assert LEADING_ROW_CONTESTED_FRACTION <= 1.0
    for step in range(steps_to_drop):
        assert front_row_present(), f"front row gone early, on step {step}"
        continent.deform(world, [neighbour], years=years_per_step, max_distance=1.5 * spacing)

    assert not front_row_present(), "sustained-override front row should have been dropped"
    assert inner_rows_present(), "inner rows must survive a leading-row drop"

    # Control: same plate, no neighbour -- nothing is contested, nothing is dropped.
    lonely = LithospherePlate(
        plate_id=0, frame=np.eye(3), crust_type="continental",
        lines=_rows(inner_phis + [front_phi], -0.5, 0.5, 40),
    )
    lonely_world = World(seed=0, plates=[lonely], mantle_centers=[], node_density=1.0)
    for _ in range(steps_to_drop + 2):
        lonely.deform(lonely_world, [], years=years_per_step, max_distance=1.5 * spacing)
    assert any(abs(ln.phi - front_phi) < 1e-6 for ln in lonely.lines)


def test_lithosphere_claim_adjacent_territory_keeps_a_margin_from_the_local_pole():
    from app.plates import POLE_CAP_MARGIN_MULT
    from app.world import World

    spacing = line_spacing_rad(1.0)
    phis = [np.pi / 2 - k * spacing for k in (POLE_CAP_MARGIN_MULT + 2, POLE_CAP_MARGIN_MULT + 1)]
    plate = _lithosphere_polar_plate(phis, np.linspace(-0.3, 0.3, 20))
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)

    for _ in range(10):
        plate.deform(world, [], years=1_000_000, max_distance=5 * spacing)

    assert max(ln.phi for ln in plate.lines) <= np.pi / 2 - POLE_CAP_MARGIN_MULT * spacing + 1e-9


# -- split / defragment keep every row one contiguous arc -----------------------------
#
# A great circle can cut a row so one side's nodes land in its interior, leaving the other
# side holding two arcs with a gap. Carrying that row whole makes outline_world / the
# row-lookup fast path claim the gap -- i.e. the sibling's own territory -- which is the
# "split/defragmentation produces overlapping siblings" long-run degradation.


def _round_plate() -> PlateWithLines:
    spacing = line_spacing_rad(1.0)
    lines = []
    for r in range(-40, 41):
        phi = r * spacing
        if abs(phi) > 0.85:
            continue
        half_width = np.sqrt(max(0.85**2 - phi**2, 0.0))
        dtheta = spacing / max(np.cos(phi), 1e-3)
        theta = np.arange(-half_width, half_width, dtheta)
        if len(theta) < 3:
            continue
        lines.append(ElevationLine(phi=phi, theta=theta, elevation=np.zeros(len(theta))))
    return PlateWithLines(plate_id=1, frame=np.eye(3), crust_type="oceanic", lines=lines)


def test_split_yields_disjoint_daughters_when_a_cut_strands_a_row():
    plate = _round_plate()
    spacing = line_spacing_rad(1.0)
    # A cut plane tilted just off this plate's own pole: the mid rows keep both their ends on
    # the positive side with the negative side biting a lens out of their interior, so the
    # positive daughter is left holding two arcs per stranded row.
    cut_normal = geometry.normalize(np.array([-0.2, 0.0, 1.0]))

    result = plate.split(new_id=2, cut_normal=cut_normal, min_nodes=1)
    assert result is not None
    a, b = result

    # No row of either daughter carries an interior gap wider than a couple of node steps.
    for daughter in (a, b):
        for line in daughter.lines:
            if len(line) < 2:
                continue
            dtheta = spacing / max(np.cos(line.phi), 1e-3)
            assert np.all(np.diff(line.theta) < 4.0 * dtheta)

    # And neither daughter's envelope claims any of the other's nodes.
    pa, _ = a.all_points_and_elevation()
    pb, _ = b.all_points_and_elevation()
    assert not np.any(a.contains_batch(pb))
    assert not np.any(b.contains_batch(pa))


# -- interior subduction: an overridden mid-row patch is carved out and keyholed -----------
#
# deform()'s end-only shrink can't reach a run of overridden nodes stranded in the *middle*
# of an oceanic row (live nodes on both sides). _grow_or_shrink_line_for_deform carves it
# out and returns the row as two contiguous ElevationLines; outline_world / contains_batch
# then keyhole the gap out instead of claiming the neighbour's lobe. Fixes the frozen
# continental-over-oceanic overlap in the seed-888151728 world.


def _oceanic_slab(plate_id: int, frame: np.ndarray, half_theta: float = 0.6) -> PlateWithLines:
    spacing = line_spacing_rad(1.0)
    lines = []
    for r in range(-18, 19):
        phi = r * spacing
        dtheta = spacing / max(np.cos(phi), 1e-3)
        theta = np.arange(-half_theta, half_theta, dtheta)
        if len(theta) < 3:
            continue
        lines.append(ElevationLine(phi=phi, theta=theta, elevation=np.full(len(theta), -3800.0)))
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type="oceanic", lines=lines)


def test_deform_carves_interior_override_into_two_arcs_and_keyholes_it():
    from app.plates import _row_intervals
    from app.world import World

    # A continental plate whose territory sits squarely over the middle of an oceanic plate's
    # rows -- a neighbour's lobe planted in the interior, the frozen-overlap geometry.
    ocean = _oceanic_slab(1, np.eye(3))
    cont_frame = geometry.rotation_matrix(np.array([0.0, 1.0, 0.0]), 0.25)
    cont = _oceanic_slab(2, cont_frame, half_theta=0.18)
    cont._crust_type = "continental"
    world = World(seed=0, plates=[ocean, cont], mantle_centers=[], node_density=1.0)

    assert all(len(ivs) == 1 for _phi, ivs in _row_intervals(list(ocean.lines)))
    cont_pts, _ = cont.all_points_and_elevation()
    assert ocean.contains_batch(cont_pts).mean() > 0.5  # lobe starts squarely inside

    for _ in range(3):
        ocean.deform(world, [cont], years=1_000_000, max_distance=5 * line_spacing_rad(1.0))

    rows = _row_intervals([ln for ln in ocean.lines if len(ln) > 0])
    split_rows = [ivs for _phi, ivs in rows if len(ivs) > 1]
    assert split_rows, "expected some oceanic rows carved into two arcs"
    # every ElevationLine stays a single contiguous arc
    for line in ocean.lines:
        if len(line) >= 2:
            dtheta = line_spacing_rad(1.0) / max(np.cos(line.phi), 1e-3)
            assert np.all(np.diff(line.theta) < 4.0 * dtheta)

    # the continental lobe's own nodes are almost entirely out of the oceanic plate's
    # polygon now -- the frozen overlap is gone.
    assert ocean.contains_batch(cont_pts).mean() < 0.05


def test_outline_world_still_one_contiguous_array_and_excludes_the_hole():
    from app.plates import _plate_outline_loops, _row_intervals

    ocean = _oceanic_slab(1, np.eye(3))
    lines = list(ocean.lines)
    spacing = line_spacing_rad(1.0)
    holed = []
    for i, ln in enumerate(lines):
        if 12 <= i <= 24 and len(ln.theta) > 24:
            keep = np.ones(len(ln.theta), dtype=bool)
            mid = len(ln.theta) // 2
            keep[mid - 4 : mid + 4] = False
            from app.elevation_lines import split_into_contiguous_runs

            holed.extend(split_into_contiguous_runs(ln.masked(keep), spacing / max(np.cos(ln.phi), 1e-3)))
        else:
            holed.append(ln)
    ocean.set_lines(holed)

    rows = _row_intervals([ln for ln in ocean.lines if len(ln) > 0])
    loops = _plate_outline_loops(rows)
    assert len(loops) == 2  # outer boundary + one hole

    poly = ocean.get_bounding_polygon()
    assert poly.ndim == 2 and poly.shape[1] == 3
    assert np.allclose(np.linalg.norm(poly, axis=-1), 1.0, atol=1e-9)

    # a point in the hole is outside the polygon and outside contains_batch;
    # a point in a still-covered row is inside.
    hole_phi = ocean.lines[len(ocean.lines) // 2].phi
    gap = [ivs for phi, ivs in rows if abs(phi - hole_phi) < 1e-9][0]
    hole_theta = (gap[0][1] + gap[1][0]) / 2.0
    hole_pt = geometry.to_world(ocean.frame, geometry.local_xyz(np.array([hole_phi]), np.array([hole_theta])))
    assert not geometry.points_in_spherical_polygon(hole_pt, poly)[0]
    assert not ocean.contains_batch(hole_pt)[0]

    edge_row = ocean.lines[0]
    solid_pt = geometry.to_world(
        ocean.frame, geometry.local_xyz(np.array([edge_row.phi]), np.array([float(np.median(edge_row.theta))]))
    )
    assert ocean.contains_batch(solid_pt)[0]


def _melt_test_plates(elevation: float):
    """Two plates spreading apart at the equator (mirrors test_boundary.py's own "self spins
    -z away from an eastward neighbour" divergent geometry), with the near-boundary end
    already just above RIFT_CRITICAL_THICKNESS_M so one real divergent step thins it past the
    threshold and triggers decompression melting there. `elevation` is the boundary end's
    *pre-melt* elevation -- the thing that decides whether the erupted material comes back
    continental (still land) or oceanic (at/below sea level), see LithospherePlate.deform."""
    from app import mantle
    from app.lithosphere import reference_thickness
    from app.lithosphere_plate import LithospherePlate

    rate = mantle.cm_per_yr_to_rad_per_yr(5.0)

    def _plate(pid, crust_type, theta_lo, theta_hi, omega_z, hc):
        theta = np.linspace(theta_lo, theta_hi, 20)
        line = ElevationLine(
            phi=0.0,
            theta=theta,
            elevation=np.full(20, elevation if pid == 0 else -3000.0),
            crustal_thickness_m=np.full(20, hc),
            mantle_lithosphere_thickness_m=np.full(20, reference_thickness(crust_type)[1]),
        )
        filler = ElevationLine(
            phi=-0.6,
            theta=np.linspace(-0.2, 0.2, 8),
            elevation=np.zeros(8),
            crustal_thickness_m=np.full(8, reference_thickness(crust_type)[0]),
            mantle_lithosphere_thickness_m=np.full(8, reference_thickness(crust_type)[1]),
        )
        return LithospherePlate(plate_id=pid, frame=np.eye(3), crust_type=crust_type, omega=np.array([0.0, 0.0, omega_z]), lines=[line, filler])

    # Hc just above the melting threshold (rheology.RIFT_CRITICAL_THICKNESS_M, 5000m) -- a
    # single divergent step's thinning is enough to cross it.
    west = _plate(0, "continental", -0.5, -0.02, -rate, hc=5050.0)
    east = _plate(1, "oceanic", 0.02, 0.5, rate, hc=reference_thickness("oceanic")[0])
    return west, east


def test_decompression_melting_above_sea_level_erupts_continental_crust():
    from app.elevation_lines import CRUST_TYPE_CONTINENTAL, ELEV_CHANGE_VOLCANO
    from app.lithosphere import REFERENCE_HC_CONTINENTAL_M
    from app.world import World

    west, east = _melt_test_plates(elevation=500.0)  # still standing above sea level
    world = World(seed=0, plates=[west, east], mantle_centers=[], node_density=1.0)
    spacing = line_spacing_rad(1.0)
    west.deform(world, [east], years=300_000, max_distance=1.5 * spacing)

    line = next(ln for ln in west.lines if abs(ln.phi - 0.0) < 1e-6)
    assert line.crustal_thickness_m[-1] == REFERENCE_HC_CONTINENTAL_M
    assert line.crust_type_code[-1] == CRUST_TYPE_CONTINENTAL
    assert line.elev_change_reason[-1] == ELEV_CHANGE_VOLCANO
    # A fresh continental reference column floats near ordinary dry land, not the deep abyss.
    assert line.elevation[-1] > 0.0
    # Untouched, far-from-the-boundary nodes are unaffected.
    assert line.crustal_thickness_m[0] == 5050.0


def test_decompression_melting_at_or_below_sea_level_erupts_oceanic_crust():
    from app.elevation_lines import CRUST_TYPE_OCEANIC, ELEV_CHANGE_VOLCANO
    from app.lithosphere import REFERENCE_HC_OCEANIC_M
    from app.world import World

    west, east = _melt_test_plates(elevation=-2000.0)  # a drowned, already-submerged margin
    world = World(seed=0, plates=[west, east], mantle_centers=[], node_density=1.0)
    spacing = line_spacing_rad(1.0)
    west.deform(world, [east], years=300_000, max_distance=1.5 * spacing)

    line = next(ln for ln in west.lines if abs(ln.phi - 0.0) < 1e-6)
    assert line.crustal_thickness_m[-1] == REFERENCE_HC_OCEANIC_M
    assert line.crust_type_code[-1] == CRUST_TYPE_OCEANIC
    assert line.elev_change_reason[-1] == ELEV_CHANGE_VOLCANO
    # A fresh oceanic reference column floats at abyssal depth, not dry land.
    assert line.elevation[-1] < -3000.0

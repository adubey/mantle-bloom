import numpy as np
from app import geometry
from app.elevation_lines import (
    CRUST_TYPE_CONTINENTAL,
    CRUST_TYPE_INHERIT,
    CRUST_TYPE_OCEANIC,
    NODE_DENSITY_CHOICES,
    line_spacing_rad,
)
from app.lithosphere_plate import build_plate_tiling, generate_plates
from app import healpix_grid
from app.plates import (
    ELLIPSE_OUTLINE_POINTS,
    MAX_AUTO_PLATES,
    MIN_AUTO_PLATES,
    MIN_OCEANIC_PLATES,
    cached_node_healpix_index,
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
from app.sparse_quad_patch import cells_per_face_edge, unpack_cell_keys
from app.world import generate_world, step_world

from .quad_fixtures import lobed_plate


def _measured_land_fraction(plates_list) -> float:
    total = sum(float(p.node_areas_m2().sum()) for p in plates_list)
    land = sum(float(p.node_areas_m2()[p.collect("elevation") > 0].sum()) for p in plates_list)
    return land / total if total else 0.0


def _measured_continental_area_fraction(plates_list) -> float:
    total = sum(float(p.node_areas_m2().sum()) for p in plates_list)
    continental = sum(float(p.node_areas_m2().sum()) for p in plates_list if p.crust_type == "continental")
    return continental / total if total else 0.0


def _remove_every_cell(plate) -> None:
    plate.remove_cells(np.ones(plate.node_count(), dtype=bool))


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


def test_every_plate_has_nodes():
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
        world_pts = p.all_points_and_elevation()[0]
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


def test_voronoi_points_changes_the_tiling_but_not_the_plate_count():
    # `voronoi_points` is the total seed-point count; generate_plates turns it into a per-plate
    # extra-site count once the final plate count is known, then merges cells down to the
    # requested plate count. Different point counts -> different (lumpier vs. smoother) plate
    # outlines, same plate count.
    sparse = generate_plates(seed=1, num_plates=8, voronoi_points=8)
    dense = generate_plates(seed=1, num_plates=8, voronoi_points=120)
    assert len(sparse) == len(dense) == 8

    def outline_signature(plates):
        return [tuple(np.round(p.outline_world(), 4).ravel().tolist()) for p in plates]

    assert outline_signature(sparse) != outline_signature(dense)


def test_voronoi_points_is_deterministic():
    a = generate_plates(seed=5, num_plates=7, voronoi_points=60)
    b = generate_plates(seed=5, num_plates=7, voronoi_points=60)
    assert [p.outline_world().tolist() for p in a] == [p.outline_world().tolist() for p in b]


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
        assert np.array_equal(p1.cell_keys, p2.cell_keys)


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
    # frame/crust type/cells): the composite relief field must reproduce the exact per-node
    # elevation for a given seed.
    e1 = np.concatenate([p.collect("elevation") for p in generate_plates(seed=321, num_plates=10)])
    e2 = np.concatenate([p.collect("elevation") for p in generate_plates(seed=321, num_plates=10)])
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


def test_get_bounding_polygon_cache_invalidated_by_a_topology_change():
    p = generate_plates(seed=10, num_plates=8)[0]
    cached = p.get_bounding_polygon()
    drop = np.zeros(p.node_count(), dtype=bool)
    drop[0] = True
    p.remove_cells(drop)
    refreshed = p.get_bounding_polygon()
    assert refreshed is not cached
    assert np.array_equal(refreshed, p.outline_world())


def test_outline_world_empty_for_plate_with_no_cells():
    p = generate_plates(seed=6, num_plates=8)[0]
    _remove_every_cell(p)
    assert len(p.outline_world()) == 0


def test_get_node_kdtree_is_cached_until_geometry_changes():
    p = generate_plates(seed=12, num_plates=8)[0]
    first = p.get_node_kdtree()
    assert first is p.get_node_kdtree()  # same cached tree
    p.rotate(geometry.plate_frame_from_seed(np.array([0.0, 1.0, 0.0])))
    rotated = p.get_node_kdtree()
    assert rotated is not first  # invalidated by rotate
    assert np.allclose(np.asarray(rotated.data), p.all_points_and_elevation()[0])
    drop = np.zeros(p.node_count(), dtype=bool)
    drop[0] = True
    p.remove_cells(drop)
    assert p.get_node_kdtree() is not rotated  # invalidated by a topology change


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
        _remove_every_cell(p)
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
        n = p.node_count()
        p.set_fields_on_plate(soil_depth=np.full(n, 2.5), coal_deposit_m=np.full(n, 1.5))
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
    `get_bounding_polygon()`. Not expected to be exactly zero -- independently rotated cell
    lattices overlap a little where they meet (issue #255) -- but it should stay bounded
    rather than climbing without limit turn over turn -- see
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
# Boundary retreat can sever one Plate's cells into two disconnected landmasses or strand a
# few cells; maybe_split_plate only cuts on mantle-flow disagreement, not geometry, so it
# never catches either. These exercise the pass that does.

_DEFRAG_SPACING_RAD = line_spacing_rad(1.0)
_DEFRAG_CONNECT_RAD = 2.5 * _DEFRAG_SPACING_RAD
_DEFRAG_I0 = cells_per_face_edge(_DEFRAG_SPACING_RAD) // 2 - 20  # see quad_fixtures.lobed_plate


def _lobed_plate(lobes, plate_id=0, rows=12, per_row=8, crust_type="oceanic", **plate_kwargs):
    """See `quad_fixtures.lobed_plate`."""
    return lobed_plate(lobes, plate_id=plate_id, rows=rows, per_row=per_row, crust_type=crust_type, **plate_kwargs)


def test_node_components_labels_isolated_clusters_separately():
    points, _ = _lobed_plate([0, 30]).all_points_and_elevation()
    labels = node_components(points, _DEFRAG_CONNECT_RAD)
    _, counts = np.unique(labels, return_counts=True)
    assert sorted(counts.tolist()) == [96, 96]  # two equal lobes, 12 rows x 8 nodes each


def test_node_components_one_label_for_a_contiguous_blob():
    points, _ = _lobed_plate([0]).all_points_and_elevation()
    assert set(node_components(points, _DEFRAG_CONNECT_RAD).tolist()) == {0}


def test_node_components_empty_input():
    assert node_components(np.zeros((0, 3)), 0.1).shape == (0,)


def test_defragment_splits_a_severed_plate_and_keeps_identity_on_the_largest():
    from app.world import World

    plate = _lobed_plate([(0, 10), (30, 6)], plate_id=7, omega=np.array([0.1, 0.2, 0.3]), age_steps=9)
    before = plate.node_count()
    world = World(seed=0, plates=[plate], mantle_centers=[])

    result = plate.defragment(next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50, world=world)
    assert result is not None
    replacements, consumed, _ = result

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
    from app.world import World

    # second lobe is 12 nodes (1 per row), well below min_fragment_nodes -- dropped, not
    # promoted to its own plate, and no new id is consumed.
    plate = _lobed_plate([(0, 10), (30, 1)], plate_id=3)
    before = plate.node_count()
    world = World(seed=0, plates=[plate], mantle_centers=[])

    result = plate.defragment(next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50, world=world)
    assert result is not None
    replacements, consumed, _ = result

    assert consumed == 0
    assert [p.plate_id for p in replacements] == [3]
    assert replacements[0].node_count() == before - 12
    # The 12 shed nodes are recorded, not silently discarded -- see World.removed_points_log.
    assert len(world.removed_points_log) == 12
    assert all(plate_id == 3 for _, _, plate_id in world.removed_points_log)


def test_defragment_leaves_a_contiguous_plate_alone():
    from app.world import World

    plate = _lobed_plate([0])
    world = World(seed=0, plates=[plate], mantle_centers=[])
    assert plate.defragment(next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50, world=world) is None


def test_defragment_leaves_an_all_debris_plate_for_the_territory_check():
    from app.world import World

    # three lobes, none reaching min_fragment_nodes: defrag declines (returns None) rather
    # than deleting a whole plate itself -- has_negligible_territory / remove_defunct_plates
    # own that call.
    plate = _lobed_plate([(0, 2), (25, 2), (50, 2)], rows=10)
    world = World(seed=0, plates=[plate], mantle_centers=[])
    assert plate.defragment(next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50, world=world) is None


def test_defragment_partition_carries_each_nodes_own_fields_to_the_right_fragment():
    from app.world import World

    plate = _lobed_plate([0, 30], plate_id=4)
    marker = np.arange(plate.node_count(), dtype=float)  # a distinct value per node
    plate.set_fields_on_plate(channel_depth=marker)

    world = World(seed=0, plates=[plate], mantle_centers=[])
    replacements, _, _ = plate.defragment(
        next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50, world=world
    )
    recombined = np.concatenate([p.collect("channel_depth") for p in replacements])
    assert sorted(recombined.tolist()) == sorted(marker.tolist())


def test_defragment_freezes_inherited_crust_when_a_fragment_changes_type():
    # Issue #239: a continental plate whose second lobe was mostly re-marked oceanic by a
    # magma-typing event. That lobe's fragment comes out oceanic, but its one
    # still-CRUST_TYPE_INHERIT node per row was never re-marked -- it must stay continental
    # rather than silently turning oceanic with its new plate.
    from app.lithosphere import node_crust_density
    from app.world import World

    plate = _lobed_plate([0, 30], plate_id=4, crust_type="continental")
    column = unpack_cell_keys(plate.cell_keys)[2] - _DEFRAG_I0
    codes = np.full(plate.node_count(), CRUST_TYPE_INHERIT, dtype=np.int8)
    codes[(column >= 30) & (column < 37)] = CRUST_TYPE_OCEANIC  # second lobe: 7 of 8 cells per row oceanic, last inherits
    plate.set_fields_on_plate(crust_type_code=codes)
    points_before, _ = plate.all_points_and_elevation()
    density_before = node_crust_density(plate.collect("crust_type_code"), plate.crust_type)

    world = World(seed=0, plates=[plate], mantle_centers=[])
    replacements, _, _ = plate.defragment(
        next_id=20, connect_radius_rad=_DEFRAG_CONNECT_RAD, min_fragment_nodes=50, world=world
    )
    by_type = {p.crust_type: p for p in replacements}
    assert set(by_type) == {"continental", "oceanic"}
    assert set(by_type["continental"].collect("crust_type_code").tolist()) == {CRUST_TYPE_INHERIT}
    oceanic_codes = by_type["oceanic"].collect("crust_type_code")
    assert np.count_nonzero(oceanic_codes == CRUST_TYPE_CONTINENTAL) == 12
    assert not np.any(oceanic_codes == CRUST_TYPE_INHERIT)

    # Every node's own crust density is exactly what it was before the partition.
    before = {tuple(np.round(pt, 12)): rho for pt, rho in zip(points_before, density_before)}
    for p in replacements:
        pts, _ = p.all_points_and_elevation()
        rho = node_crust_density(p.collect("crust_type_code"), p.crust_type)
        assert [before[tuple(np.round(pt, 12))] for pt in pts] == rho.tolist()


def test_has_negligible_territory_false_for_a_plate_with_real_territory():
    assert not _lobed_plate([0]).has_negligible_territory()


# -- Issue #133 phase 2: cached_node_healpix_index -------------------------------------------


def test_cached_node_healpix_index_query_matches_cktree_for_nodes_that_own_their_pixel():
    """Same contract as healpix_grid's own `NodePixelIndex` regression test, but through the
    shared `plates.cached_node_healpix_index` entry point every "healpix"
    `node_cloud_resample_mode` caller now goes through (render_image._node_cloud_and_tree,
    climate._sample_elevation_and_crust)."""
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(2)
    xyz = rng.normal(size=(500, 3))
    xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)

    index = cached_node_healpix_index(None, xyz)
    assert not np.any(index.pixel_to_node == -1)

    pix = index.grid.ang2pix(np.arctan2(xyz[:, 1], xyz[:, 0]), np.arcsin(xyz[:, 2]))
    owns_own_pixel = index.pixel_to_node[pix] == np.arange(len(xyz))
    assert owns_own_pixel.sum() > len(xyz) // 2

    _, idx = index.query(xyz[owns_own_pixel])
    tree = cKDTree(xyz)
    _, tree_idx = tree.query(xyz[owns_own_pixel])
    assert np.array_equal(idx, tree_idx)


def test_cached_node_healpix_index_world_none_always_builds_fresh():
    rng = np.random.default_rng(3)
    xyz_a = rng.normal(size=(200, 3))
    xyz_a /= np.linalg.norm(xyz_a, axis=1, keepdims=True)
    xyz_b = rng.normal(size=(200, 3))
    xyz_b /= np.linalg.norm(xyz_b, axis=1, keepdims=True)

    index_a = cached_node_healpix_index(None, xyz_a)
    index_b = cached_node_healpix_index(None, xyz_b)
    assert not np.array_equal(index_a.pixel_to_node, index_b.pixel_to_node)


def test_cached_node_healpix_index_shares_grid_and_index_across_world_callers():
    """The whole point of moving this out of render_image.py's own private helper: a second
    caller against the same world this step (climate.py, per issue #133 phase 2) must reuse
    the first caller's scatter+fill, not rebuild -- the same "first caller wins" sharing
    `cached_node_position_tree` already provides for "kdtree" mode."""
    world = generate_world(5, num_plates=8)
    all_points, _all_elev, _all_owner = collect_all_points(world.plates)

    first = cached_node_healpix_index(world, all_points)
    assert world.node_healpix_grid_cache is not None
    assert world.node_healpix_index_cache is first

    second = cached_node_healpix_index(world, all_points)
    assert second is first  # reused, not rebuilt

    grid_cache_before = world.node_healpix_grid_cache
    step_world(world, 1_000_000)
    assert world.node_healpix_index_cache is None  # a node moved -- must rebuild
    assert world.node_healpix_grid_cache is grid_cache_before  # node count unchanged -- reused

    third = cached_node_healpix_index(world, all_points)
    assert third is not first

import numpy as np
from app import geometry, mantle
from app.elevation_lines import ElevationLine
from app.plates import Plate, PlateWithMesh, _contested_by_any, _max_boundary_effect_rad, line_spacing_rad
from app.world import World

_FRAME = geometry.plate_frame_from_seed(np.array([1.0, 0.0, 0.0]))


def _disk_cluster(rng, n, radius_rad):
    """A roughly disk-shaped cloud of (theta, phi) local coordinates around the origin --
    same helper test_plate_with_rtree.py's own _disk_cluster uses, kept separate per that
    file's own convention of not sharing test helpers across test modules."""
    r = radius_rad * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    angle = rng.uniform(0.0, 2 * np.pi, size=n)
    return r * np.cos(angle), r * np.sin(angle)


def test_plate_with_mesh_conforms_to_the_plate_interface():
    plate = PlateWithMesh(plate_id=1, frame=_FRAME, crust_type="continental")
    assert isinstance(plate, Plate)


def test_empty_plate_has_no_nodes():
    plate = PlateWithMesh(plate_id=1, frame=_FRAME, crust_type="oceanic")
    assert plate.node_count() == 0
    points, elevation = plate.all_points_and_elevation()
    assert points.shape == (0, 3)
    assert elevation.shape == (0,)
    assert plate.outline_world().shape == (0, 3)
    assert plate.collect("elevation").shape == (0,)
    assert plate.collect("is_volcano").dtype == bool


def test_all_points_and_elevation_matches_theta_phi_on_the_unit_sphere():
    rng = np.random.default_rng(0)
    n = 50
    theta, phi = _disk_cluster(rng, n, radius_rad=0.2)
    elevation = rng.uniform(-500.0, 500.0, size=n)
    plate = PlateWithMesh(plate_id=2, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=elevation)

    assert plate.node_count() == n
    points, got_elevation = plate.all_points_and_elevation()
    assert points.shape == (n, 3)
    assert np.allclose(np.linalg.norm(points, axis=-1), 1.0, atol=1e-9)
    assert np.array_equal(got_elevation, elevation)

    expected = geometry.to_world(_FRAME, geometry.local_xyz(phi, theta))
    assert np.allclose(points, expected)


def test_optional_fields_default_to_zero_and_is_volcano_defaults_to_false():
    theta = np.array([0.0, 0.1, 0.2])
    phi = np.zeros(3)
    plate = PlateWithMesh(plate_id=3, frame=_FRAME, crust_type="oceanic", theta=theta, phi=phi, elevation=np.zeros(3))
    for name in ElevationLine.OPTIONAL_FIELDS:
        values = plate.collect(name)
        assert values.shape == theta.shape
        if name == "is_volcano":
            assert values.dtype == bool
            assert not np.any(values)
        else:
            assert np.all(values == 0.0)


def test_collect_returns_the_field_passed_at_construction():
    theta = np.array([0.0, 0.1, 0.2])
    phi = np.zeros(3)
    plate = PlateWithMesh(
        plate_id=4,
        frame=_FRAME,
        crust_type="continental",
        theta=theta,
        phi=phi,
        elevation=np.zeros(3),
        soil_depth=np.array([1.0, 2.0, 3.0]),
        is_volcano=np.array([True, False, True]),
    )
    assert np.array_equal(plate.collect("soil_depth"), [1.0, 2.0, 3.0])
    assert np.array_equal(plate.collect("is_volcano"), [True, False, True])
    assert np.all(plate.collect("channel_depth") == 0.0)  # untouched field still defaults


def test_set_nodes_replaces_every_node_and_rebuilds_the_triangulation():
    plate = PlateWithMesh(
        plate_id=5, frame=_FRAME, crust_type="oceanic", theta=np.array([0.0]), phi=np.array([0.0]), elevation=np.array([1.0])
    )
    assert plate.node_count() == 1
    first_triangulation = plate._get_triangulation()

    new_theta = np.array([0.0, 0.1, 0.2, 0.3])
    new_phi = np.zeros(4)
    plate.set_nodes(new_theta, new_phi, np.array([1.0, 2.0, 3.0, 4.0]))
    assert plate.node_count() == 4
    assert plate._get_triangulation() is not first_triangulation


def test_map_world_points_on_plate_uses_plate_wide_theta_range():
    theta = np.array([0.0, 0.1, 0.2, 0.3])
    phi = np.array([0.0, 0.05, -0.05, 0.02])
    plate = PlateWithMesh(plate_id=9, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(4))

    results = list(plate.map_world_points_on_plate())
    assert len(results) == 4
    fractions = np.array([fraction for _, _, fraction in results])
    assert np.allclose(fractions, [0.0, 1 / 3, 2 / 3, 1.0])

    expected_world = geometry.to_world(_FRAME, geometry.local_xyz(phi, theta))
    got_world = np.array([xyz for _, xyz, _ in results])
    assert np.allclose(got_world, expected_world)


def test_map_world_points_on_plate_empty_for_no_nodes():
    plate = PlateWithMesh(plate_id=10, frame=_FRAME, crust_type="oceanic")
    assert list(plate.map_world_points_on_plate()) == []


def test_map_world_points_on_plate_single_node_gives_midpoint_fraction():
    plate = PlateWithMesh(
        plate_id=11, frame=_FRAME, crust_type="oceanic", theta=np.array([0.5]), phi=np.array([0.0]), elevation=np.array([0.0])
    )
    ((_, _, fraction),) = list(plate.map_world_points_on_plate())
    assert fraction == 0.5


def test_outline_world_is_a_finite_unit_vector_loop_for_a_disk_cluster():
    rng = np.random.default_rng(1)
    n = 400
    theta, phi = _disk_cluster(rng, n, radius_rad=0.3)
    elevation = np.zeros(n)
    plate = PlateWithMesh(plate_id=6, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=elevation)

    outline = plate.outline_world()
    assert outline.ndim == 2 and outline.shape[1] == 3
    assert len(outline) >= 3
    assert len(outline) < n  # the boundary loop of a disk is only its own rim, not every node
    assert np.all(np.isfinite(outline))
    assert np.allclose(np.linalg.norm(outline, axis=-1), 1.0, atol=1e-9)


def test_outline_world_falls_back_to_every_node_below_the_hull_threshold():
    theta = np.array([0.0, 0.1])
    phi = np.array([0.0, 0.0])
    plate = PlateWithMesh(plate_id=7, frame=_FRAME, crust_type="oceanic", theta=theta, phi=phi, elevation=np.zeros(2))
    outline = plate.outline_world()
    assert outline.shape == (2, 3)


def test_outline_world_does_not_crash_on_collinear_nodes():
    theta = np.linspace(0.0, 1.0, 20)
    phi = np.zeros(20)
    plate = PlateWithMesh(plate_id=8, frame=_FRAME, crust_type="oceanic", theta=theta, phi=phi, elevation=np.zeros(20))
    outline = plate.outline_world()
    assert np.all(np.isfinite(outline))
    assert np.allclose(np.linalg.norm(outline, axis=-1), 1.0, atol=1e-9)


def test_outline_world_traces_a_concave_notch():
    # A disk cluster with a wedge (theta in [-0.5, 0.5] rad, i.e. a ~60-degree slice) carved
    # out -- proves the alpha-shape boundary trace follows the real concave shape, unlike
    # PlateWithRTree.outline_world's convex-hull-only approximation, which would silently
    # claim the notch as still-owned territory.
    rng = np.random.default_rng(12)
    n = 1500
    theta, phi = _disk_cluster(rng, n, radius_rad=0.4)
    wedge = (theta > -0.5) & (theta < 0.5) & (phi > 0.0)
    keep = ~wedge
    theta, phi = theta[keep], phi[keep]
    plate = PlateWithMesh(plate_id=16, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(len(theta)))

    notch_lat, notch_lon = 0.2, 0.0  # inside the carved-out wedge (phi > 0, |theta| < 0.5)
    rim_lat, rim_lon = -0.2, 0.0  # inside the retained lower half, untouched by the wedge
    assert not plate.contains(notch_lat, notch_lon)
    assert plate.contains(rim_lat, rim_lon)


def test_get_bounding_polygon_matches_outline_world():
    rng = np.random.default_rng(2)
    theta, phi = _disk_cluster(rng, 100, radius_rad=0.3)
    plate = PlateWithMesh(plate_id=12, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(100))
    assert np.array_equal(plate.get_bounding_polygon(), plate.outline_world())


def test_get_bounding_polygon_returns_the_same_cached_array_until_invalidated():
    rng = np.random.default_rng(3)
    theta, phi = _disk_cluster(rng, 100, radius_rad=0.3)
    plate = PlateWithMesh(plate_id=13, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(100))
    first = plate.get_bounding_polygon()
    second = plate.get_bounding_polygon()
    assert first is second


def test_get_bounding_polygon_cache_invalidated_by_set_nodes():
    rng = np.random.default_rng(4)
    theta, phi = _disk_cluster(rng, 100, radius_rad=0.3)
    plate = PlateWithMesh(plate_id=14, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(100))
    cached = plate.get_bounding_polygon()
    new_theta, new_phi = _disk_cluster(rng, 100, radius_rad=0.3)
    plate.set_nodes(new_theta, new_phi, np.zeros(100))
    refreshed = plate.get_bounding_polygon()
    assert refreshed is not cached
    assert np.array_equal(refreshed, plate.outline_world())


def test_get_bounding_polygon_cache_invalidated_by_rotate():
    plate = PlateWithMesh(
        plate_id=15, frame=_FRAME, crust_type="oceanic", theta=np.array([0.0]), phi=np.array([0.0]), elevation=np.array([0.0])
    )
    cached = plate.get_bounding_polygon()
    plate.rotate(geometry.plate_frame_from_seed(np.array([0.0, 1.0, 0.0])))
    refreshed = plate.get_bounding_polygon()
    assert refreshed is not cached
    assert np.array_equal(refreshed, plate.outline_world())


def test_triangulation_cache_survives_rotate_and_is_only_invalidated_by_set_nodes():
    rng = np.random.default_rng(5)
    theta, phi = _disk_cluster(rng, 60, radius_rad=0.2)
    plate = PlateWithMesh(plate_id=17, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(60))
    triangulation = plate._get_triangulation()

    plate.rotate(geometry.plate_frame_from_seed(np.array([0.0, 1.0, 0.0])))
    assert plate._get_triangulation() is triangulation  # unchanged: local theta/phi didn't move

    plate.set_nodes(theta, phi, np.zeros(60))
    assert plate._get_triangulation() is not triangulation


def test_shift_rotates_rigidly_and_returns_max_displacement():
    rng = np.random.default_rng(6)
    theta, phi = _disk_cluster(rng, 30, radius_rad=0.1)
    plate = PlateWithMesh(plate_id=18, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(30))
    off_axis = geometry.to_world(_FRAME, geometry.local_xyz(np.array([0.3]), np.array([0.3]))[0])
    centers = [mantle.ConvectionCenter(position=off_axis, strength=mantle.MANTLE_FLOW_REFERENCE_RATE * 10, falloff=0.5)]
    world = World(seed=0, plates=[plate], mantle_centers=centers, node_density=1.0)

    before, _ = plate.all_points_and_elevation()
    d = plate.shift(world, years=1_000_000)
    after, _ = plate.all_points_and_elevation()

    assert d >= 0.0
    assert not np.allclose(before, after)
    assert np.isclose(d, float(geometry.angular_distance(before, after).max()))


def test_shift_on_an_empty_plate_returns_zero():
    plate = PlateWithMesh(plate_id=19, frame=_FRAME, crust_type="oceanic")
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)
    assert plate.shift(world, years=1_000_000) == 0.0


def test_split_partitions_by_cut_normal_and_preserves_every_field_exactly():
    theta = np.linspace(-0.3, 0.3, 40)
    phi = np.zeros(40)
    soil_depth = np.arange(40, dtype=float)
    plate = PlateWithMesh(
        plate_id=20, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(40), soil_depth=soil_depth
    )
    world_pts, _ = plate.all_points_and_elevation()
    centroid = geometry.normalize(world_pts.mean(axis=0))
    cut_normal = geometry.normalize(np.cross(centroid, np.array([0.0, 0.0, 1.0])))

    result = plate.split(new_id=99, cut_normal=cut_normal, min_nodes=1)
    assert result is not None
    plate_a, plate_b = result
    assert plate_a.plate_id == 20
    assert plate_b.plate_id == 99
    assert plate_a.node_count() + plate_b.node_count() == 40
    assert set(plate_a.collect("soil_depth").tolist()) | set(plate_b.collect("soil_depth").tolist()) == set(soil_depth.tolist())


def test_split_returns_none_below_min_nodes():
    theta = np.array([0.0, 0.001])
    phi = np.zeros(2)
    plate = PlateWithMesh(plate_id=21, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(2))
    world_pts, _ = plate.all_points_and_elevation()
    cut_normal = geometry.normalize(np.cross(world_pts.mean(axis=0), np.array([0.0, 1.0, 0.0])))
    assert plate.split(new_id=1, cut_normal=cut_normal, min_nodes=5) is None


def test_merge_with_unions_nodes_and_preserves_every_field_exactly():
    theta_a, phi_a = _disk_cluster(np.random.default_rng(7), 20, radius_rad=0.05)
    plate_a = PlateWithMesh(
        plate_id=22,
        frame=_FRAME,
        crust_type="continental",
        theta=theta_a,
        phi=phi_a,
        elevation=np.zeros(20),
        soil_depth=np.full(20, 5.0),
        is_volcano=np.ones(20, dtype=bool),
    )
    other_frame = geometry.plate_frame_from_seed(np.array([1.0, 0.0, 0.0]))
    theta_b, phi_b = _disk_cluster(np.random.default_rng(8), 15, radius_rad=0.05)
    theta_b += 0.5  # keep the two clusters from overlapping
    plate_b = PlateWithMesh(
        plate_id=23,
        frame=other_frame,
        crust_type="continental",
        theta=theta_b,
        phi=phi_b,
        elevation=np.zeros(15),
        soil_depth=np.full(15, 9.0),
    )

    plate_a.merge_with(plate_b, spacing_rad=0.01, coverage_radius_rad=0.05, other_points_xyz=np.zeros((0, 3)))
    assert plate_a.node_count() == 35
    soil = plate_a.collect("soil_depth")
    assert np.count_nonzero(soil == 5.0) == 20
    assert np.count_nonzero(soil == 9.0) == 15
    is_volcano = plate_a.collect("is_volcano")
    assert np.count_nonzero(is_volcano) == 20  # every original plate_a node's True survives exactly


def test_grow_into_appends_points_with_default_fields():
    theta, phi = _disk_cluster(np.random.default_rng(9), 20, radius_rad=0.05)
    plate = PlateWithMesh(plate_id=24, frame=_FRAME, crust_type="oceanic", theta=theta, phi=phi, elevation=np.zeros(20))
    new_local = geometry.local_xyz(np.array([0.0, 0.01]), np.array([0.5, 0.51]))
    new_world = geometry.to_world(_FRAME, new_local)
    new_elevation = np.array([100.0, 200.0])

    plate.grow_into(new_world, new_elevation, coverage_radius_rad=0.05, spacing_rad=0.01)
    assert plate.node_count() == 22
    assert np.count_nonzero(plate.collect("elevation") == 100.0) == 1
    assert np.count_nonzero(plate.collect("elevation") == 200.0) == 1
    assert not np.any(plate.collect("is_volcano")[-2:])


def _isolated_plate(plate_id, n=200, radius_rad=0.15):
    """A single compact plate with no neighbours -- deform()'s open-boundary growth path."""
    theta, phi = _disk_cluster(np.random.default_rng(plate_id), n, radius_rad)
    return PlateWithMesh(plate_id=plate_id, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=np.zeros(n))


def test_deform_isolated_plate_grows_its_boundary_outward():
    plate = _isolated_plate(30)
    spacing_rad = line_spacing_rad(1.0)
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)
    before_count = plate.node_count()

    plate.deform(world, other_plates=[], years=3_000_000, max_distance=spacing_rad * 5)
    assert plate.node_count() > before_count


def test_deform_never_moves_the_boundary_by_more_than_max_distance_worth_of_nodes():
    plate = _isolated_plate(31)
    spacing_rad = line_spacing_rad(1.0)
    world = World(seed=1, plates=[plate], mantle_centers=[], node_density=1.0)
    before_count = plate.node_count()
    max_distance = spacing_rad * 3  # allow only ~3 nodes' worth of growth this call

    plate.deform(world, other_plates=[], years=3_000_000, max_distance=max_distance)
    grown = plate.node_count() - before_count
    assert 0 <= grown  # some growth is expected...
    # ...but not unbounded: no single boundary vertex should have grown further than the cap.
    n_distance_cap = max(1, int(max_distance / spacing_rad))
    max_extend_nodes = 400  # MAX_EXTEND_NODES_PER_STEP at node_density=1.0
    assert grown <= before_count * min(n_distance_cap, max_extend_nodes)


def _hex_disk(radius_rad, spacing_rad):
    """An evenly-spaced hex-lattice disk of local (theta, phi) offsets, ~spacing_rad apart --
    unlike a uniformly random `_disk_cluster` (whose real point density has nothing to do
    with `spacing_rad` and so triggers heavy, dominating remesh-pass merging the moment
    deform() runs, swamping whatever contested/shrink/grow signal a test is actually after),
    this matches the density a plate generated through the normal lattice-sweep path
    (`iter_local_lattice`) would already have, so deform()'s own remesh pass stays a mild
    boundary correction rather than the dominant effect."""
    q_range = int(radius_rad / spacing_rad) + 2
    theta_list, phi_list = [], []
    for q in range(-q_range, q_range + 1):
        for r in range(-q_range, q_range + 1):
            if abs(q + r) > q_range:
                continue
            theta = spacing_rad * (q + r / 2.0)
            phi = spacing_rad * (r * (3**0.5) / 2.0)
            if theta**2 + phi**2 <= radius_rad**2:
                theta_list.append(theta)
                phi_list.append(phi)
    return np.array(theta_list), np.array(phi_list)


def _overlapping_disk_pair(id_a, crust_a, id_b, crust_b):
    """Two compact, genuinely 2D disk-shaped plates (not a single collinear row -- a real
    plate footprint is a 2D patch, and the Delaunay-on-a-sphere triangulation is degenerate
    for a perfectly collinear point set, same as any representation's own outline for a
    single-row plate) with substantially overlapping territory, for deform()'s contested-
    classification tests."""
    spacing_rad = line_spacing_rad(1.0)
    theta_a, phi_a = _hex_disk(0.08, spacing_rad)
    plate_a = PlateWithMesh(
        plate_id=id_a, frame=_FRAME, crust_type=crust_a, theta=theta_a, phi=phi_a,
        elevation=np.full(len(theta_a), -3000.0 if crust_a == "oceanic" else 200.0),
    )
    other_frame = geometry.plate_frame_from_seed(geometry.to_world(_FRAME, geometry.local_xyz(np.array([0.0]), np.array([0.06]))[0]))
    theta_b, phi_b = _hex_disk(0.08, spacing_rad)
    plate_b = PlateWithMesh(
        plate_id=id_b, frame=other_frame, crust_type=crust_b, theta=theta_b, phi=phi_b,
        elevation=np.full(len(theta_b), -3000.0 if crust_b == "oceanic" else 200.0),
    )
    return plate_a, plate_b


def test_deform_contested_oceanic_plate_shrinks():
    # Checks that at least one originally-contested node was actually removed -- not that
    # the plate's total node count went down, since this same deform() call also grows the
    # plate's own far, uncontested boundary (see test_deform_isolated_plate_grows_its_
    # boundary_outward) and can run its remesh pass (which can itself insert nodes by
    # splitting an overlong edge) -- for a small test plate, growth+remeshing elsewhere can
    # outweigh a handful of shrunk nodes on net, without the shrink itself having failed.
    ocean, continent = _overlapping_disk_pair(40, "oceanic", 41, "continental")
    own_points_before, _ = ocean.all_points_and_elevation()
    neighbours = ocean.get_neighbours([continent], threshold_rad=_max_boundary_effect_rad(line_spacing_rad(1.0)))
    contested_before = _contested_by_any(own_points_before, neighbours)
    assert np.any(contested_before)  # sanity: the two disks really do overlap
    contested_points_before = own_points_before[contested_before]

    world = World(seed=0, plates=[ocean, continent], mantle_centers=[], node_density=1.0)
    ocean.deform(world, other_plates=[continent], years=3_000_000, max_distance=line_spacing_rad(1.0) * 5)

    own_points_after, _ = ocean.all_points_and_elevation()
    still_present = np.array(
        [np.any(np.all(np.isclose(own_points_after, p), axis=-1)) for p in contested_points_before]
    )
    assert not np.all(still_present)  # at least one contested node was actually removed


def test_deform_contested_continental_plate_uplifts_without_shrinking():
    continent_a, continent_b = _overlapping_disk_pair(42, "continental", 43, "continental")
    world = World(seed=0, plates=[continent_a, continent_b], mantle_centers=[], node_density=1.0)
    before_count = continent_a.node_count()
    before_max_elevation = float(continent_a.collect("elevation").max())  # a uniform 200.0 before deform

    continent_a.deform(world, other_plates=[continent_b], years=3_000_000, max_distance=line_spacing_rad(1.0) * 5)
    # Continental crust crumples in place rather than subducting -- node_count can still grow
    # (this plate's own far/uncontested boundary is open territory and grows normally, same
    # as test_deform_isolated_plate_grows_its_boundary_outward), but it must never shrink.
    assert continent_a.node_count() >= before_count
    assert continent_a.collect("elevation").max() > before_max_elevation  # contested nodes uplifted


def test_deform_stretch_volcano_probability_spawns_a_volcano_when_forced(monkeypatch):
    plate = _isolated_plate(50, n=80, radius_rad=0.1)
    spacing_rad = line_spacing_rad(1.0)
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)

    class _AlwaysOverstretchedRNG:
        """Stands in for np.random.Generator everywhere deform() constructs one -- both the
        per-vertex growth-event RNG (only .random()/.uniform() with a plain float range) and
        (for a continental plate) SphereNoise's own fault-noise RNG (.normal()/.uniform() for
        direction/phase, see noise.py) -- so this needs to answer all three without crashing,
        not just the two the growth-event roll itself uses."""

        def random(self):
            return 0.0  # always below STRETCH_VOLCANO_PROBABILITY

        def uniform(self, low, high, size=None):
            return np.full(size, low) if size is not None else low

        def normal(self, size=None):
            return np.ones(size) if size is not None else 1.0

    monkeypatch.setattr(np.random, "default_rng", lambda *args, **kwargs: _AlwaysOverstretchedRNG())
    plate.deform(world, other_plates=[], years=3_000_000, max_distance=spacing_rad * 5)

    assert np.any(plate.collect("is_volcano"))


def test_deform_no_neighbours_does_not_crash_on_an_empty_plate():
    plate = PlateWithMesh(plate_id=60, frame=_FRAME, crust_type="continental")
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)
    plate.deform(world, other_plates=[], years=3_000_000, max_distance=0.01)
    assert plate.node_count() == 0

import numpy as np
from app import geometry
from app.elevation_lines import ElevationLine
from app.plates import Plate, PlateWithRTree

_FRAME = geometry.plate_frame_from_seed(np.array([1.0, 0.0, 0.0]))


def _disk_cluster(rng, n, radius_rad):
    """A roughly disk-shaped cloud of (theta, phi) local coordinates around the origin --
    stands in for "a plate's own territory" the way a compact PlateWithLines footprint
    would, without needing a real generate_plates call."""
    r = radius_rad * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    angle = rng.uniform(0.0, 2 * np.pi, size=n)
    return r * np.cos(angle), r * np.sin(angle)


def test_plate_with_rtree_conforms_to_the_plate_interface():
    plate = PlateWithRTree(plate_id=1, frame=_FRAME, crust_type="continental")
    assert isinstance(plate, Plate)


def test_empty_plate_has_no_nodes():
    plate = PlateWithRTree(plate_id=1, frame=_FRAME, crust_type="oceanic")
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
    plate = PlateWithRTree(plate_id=2, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=elevation)

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
    plate = PlateWithRTree(plate_id=3, frame=_FRAME, crust_type="oceanic", theta=theta, phi=phi, elevation=np.zeros(3))
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
    plate = PlateWithRTree(
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


def test_set_nodes_replaces_every_node_and_rebuilds_the_index():
    plate = PlateWithRTree(
        plate_id=5, frame=_FRAME, crust_type="oceanic", theta=np.array([0.0]), phi=np.array([0.0]), elevation=np.array([1.0])
    )
    assert plate.node_count() == 1

    new_theta = np.array([0.0, 0.1, 0.2, 0.3])
    new_phi = np.zeros(4)
    plate.set_nodes(new_theta, new_phi, np.array([1.0, 2.0, 3.0, 4.0]))
    assert plate.node_count() == 4
    assert plate.rtree.points.shape == (4, 2)
    nearest = plate.rtree.nearest_one(np.array([0.05, 0.0]))
    assert nearest is not None
    assert nearest[0] in (0, 1)  # closest to either theta=0.0 or theta=0.1


def test_outline_world_is_a_finite_unit_vector_hull_for_a_disk_cluster():
    rng = np.random.default_rng(1)
    n = 400
    theta, phi = _disk_cluster(rng, n, radius_rad=0.3)
    elevation = np.zeros(n)
    plate = PlateWithRTree(plate_id=6, frame=_FRAME, crust_type="continental", theta=theta, phi=phi, elevation=elevation)

    outline = plate.outline_world()
    assert outline.ndim == 2 and outline.shape[1] == 3
    assert len(outline) >= 3
    assert len(outline) < n  # the hull of a disk is only its own rim, not every node
    assert np.all(np.isfinite(outline))
    assert np.allclose(np.linalg.norm(outline, axis=-1), 1.0, atol=1e-9)


def test_outline_world_falls_back_to_every_node_below_the_hull_threshold():
    theta = np.array([0.0, 0.1])
    phi = np.array([0.0, 0.0])
    plate = PlateWithRTree(plate_id=7, frame=_FRAME, crust_type="oceanic", theta=theta, phi=phi, elevation=np.zeros(2))
    outline = plate.outline_world()
    assert outline.shape == (2, 3)


def test_outline_world_does_not_crash_on_collinear_nodes():
    theta = np.linspace(0.0, 1.0, 20)
    phi = np.zeros(20)
    plate = PlateWithRTree(plate_id=8, frame=_FRAME, crust_type="oceanic", theta=theta, phi=phi, elevation=np.zeros(20))
    outline = plate.outline_world()
    assert np.all(np.isfinite(outline))
    assert np.allclose(np.linalg.norm(outline, axis=-1), 1.0, atol=1e-9)

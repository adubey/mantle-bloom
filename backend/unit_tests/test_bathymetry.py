import numpy as np
from app import bathymetry, geometry
from app.plates import ElevationLine, PlateWithLines
from app.world import World, generate_world


def _plate(plate_id, crust_type, theta, base_elevation):
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(len(theta), base_elevation))
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=[line])


def test_submerged_continental_crust_relaxes_toward_shelf_or_deep_target():
    # Land cluster near theta=0; two submerged continental points at ~133km (inside the
    # shelf range) and ~324km (outside it, deep water).
    land = _plate(0, "continental", [-0.001, -0.0009, -0.0011], base_elevation=200.0)
    submerged = _plate(1, "continental", [0.02, 0.05], base_elevation=-10.0)
    world = World(seed=0, plates=[land, submerged])

    bathymetry.apply_bathymetry(world, years=5_000_000)

    updated = next(p for p in world.plates if p.plate_id == 1).lines[0].elevation
    # Shelf point relaxes toward SHELF_TARGET_M (-100), deep point toward
    # DEEP_CONTINENTAL_TARGET_M (-3000) -- both starting from -10.
    assert updated[0] < -10.0
    assert updated[0] > bathymetry.SHELF_TARGET_M  # relaxing toward it, not there yet
    assert updated[1] < updated[0]  # the deep point dropped much further
    assert updated[1] > bathymetry.DEEP_CONTINENTAL_TARGET_M


def test_land_elevation_untouched():
    land = _plate(0, "continental", [-0.001, -0.0009, -0.0011], base_elevation=200.0)
    submerged = _plate(1, "continental", [0.02], base_elevation=-10.0)
    world = World(seed=0, plates=[land, submerged])

    bathymetry.apply_bathymetry(world, years=5_000_000)

    assert np.array_equal(next(p for p in world.plates if p.plate_id == 0).lines[0].elevation, np.full(3, 200.0))


def test_oceanic_crust_untouched():
    land = _plate(0, "continental", [-0.001, -0.0009, -0.0011], base_elevation=200.0)
    oceanic = _plate(1, "oceanic", [0.02, 0.05], base_elevation=-3800.0)
    world = World(seed=0, plates=[land, oceanic])

    bathymetry.apply_bathymetry(world, years=5_000_000)

    assert np.array_equal(next(p for p in world.plates if p.plate_id == 1).lines[0].elevation, np.full(2, -3800.0))


def test_noop_with_no_land():
    all_submerged = _plate(0, "continental", [0.0, 0.1], base_elevation=-500.0)
    world = World(seed=0, plates=[all_submerged])

    bathymetry.apply_bathymetry(world, years=5_000_000)

    assert np.array_equal(all_submerged.lines[0].elevation, np.full(2, -500.0))


def test_noop_for_empty_world():
    world = World(seed=0, plates=[])
    bathymetry.apply_bathymetry(world, years=1_000_000)  # must not raise
    assert world.plates == []

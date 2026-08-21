import numpy as np
from app import bathymetry, geometry
from app.plates import ElevationLine, Plate
from app.world import World, generate_world


def _plate(plate_id, crust_type, theta, base_elevation):
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(len(theta), base_elevation))
    return Plate(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=[line])


def test_stepping_a_real_world_keeps_bathymetry_well_formed():
    world = generate_world(seed=15, num_plates=10, continental_fraction=0.5)
    from app.world import step_world

    for _ in range(10):
        step_world(world, years=2_000_000)

    for plate in world.plates:
        for line in plate.lines:
            assert np.all(np.isfinite(line.elevation))
            assert np.all(line.elevation >= bathymetry.MIN_ELEVATION_M)
            assert np.all(line.elevation <= bathymetry.MAX_ELEVATION_M)

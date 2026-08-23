import numpy as np
from app import climate, geometry
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def _world(seed=1, num_plates=12, continental_fraction=0.7, land_fraction=0.29, steps=0, years=5_000_000):
    world = generate_world(seed, num_plates=num_plates, continental_fraction=continental_fraction, land_fraction=land_fraction)
    for _ in range(steps):
        step_world(world, years)
    return world


def test_compute_climate_produces_finite_fields_for_a_real_world():
    # Threads the world's own climate_density through, same as every production caller
    # (erosion.py/main.py) -- world.climate_cache (populated by the step_world calls in
    # _world) is already shaped at the world's own density, so a bare compute_climate(world)
    # call here (defaulting to the *reference* GRID_HEIGHT/GRID_WIDTH) would mismatch against
    # it as soon as the world's density differs from the reference.
    world = _world(steps=3)
    height, width = climate.grid_dimensions(world.climate_density)
    fields = climate.compute_climate(world, height, width)
    for name in [
        "elevation_m", "land_temperature_c", "ocean_temperature_c", "air_temperature_c",
        "wind_u", "wind_v", "current_u", "current_v", "humidity", "precipitation_mm",
    ]:
        arr = getattr(fields, name)
        assert arr.shape == (height, width)
        assert np.all(np.isfinite(arr)), f"{name} has non-finite values"


def test_step_world_populates_climate_cache_via_erosion():
    world = _world(seed=4, num_plates=8, steps=0)
    assert world.climate_cache is None
    step_world(world, years=1_000_000)
    assert world.climate_cache is not None


def test_stats_and_render_reuse_the_same_step_cache():
    from app import render_image, stats

    world = _world(seed=5, num_plates=8, steps=1)
    cached = world.climate_cache
    assert cached is not None

    stats.compute_stats(world)
    assert world.climate_cache is cached  # compute_stats didn't replace it with a fresh one

    render_image.render_png(world, "eckert4", "temperature", 200, 120)
    assert world.climate_cache is cached  # neither did rendering a climate view

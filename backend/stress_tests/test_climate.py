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

    # precipitation/humidity/biome renders now pull the diagnostic fields via
    # compute_climate_cached (not the CFD state) -- still a same-turn cache reuse, not a
    # recompute.
    render_image.render_png(world, "eckert4", "precipitation", 200, 120)
    assert world.climate_cache is cached


def test_diagnostic_wind_model_tracks_the_cfd_biome_map():
    # The "ABL" wind model (World.wind_model == "diagnostic") is meant to reproduce *most* of
    # the shallow-water CFD's downstream climate for a fraction of the cost -- see
    # docs/simulation-model.md#wind-model and TODO.md. This locks in that it stays a decent
    # approximation, not a regression guard on an exact value: two copies of the same world
    # stepped the same number of times, one CFD, one diagnostic, should agree on the large
    # majority of the land biome map and keep precipitation close. A drop here most likely
    # means compute_wind, compute_air_temperature_diagnostic, or the biome thresholds moved.
    # Characterised at climate_density/fluid_density 2.0 (see the "~84-89%" figure in
    # docs/simulation-model.md#wind-model); the agreement is somewhat config-sensitive, so
    # the floor here is deliberately loose.
    kw = dict(num_plates=10, continental_fraction=0.7, land_fraction=0.29, climate_density=2.0, fluid_density=2.0)
    cfd_world = generate_world(7, **kw)
    abl_world = generate_world(7, **kw)
    abl_world.wind_model = "diagnostic"
    for _ in range(10):
        step_world(cfd_world, years=2_000_000)
        step_world(abl_world, years=2_000_000)

    cfd = climate.compute_climate_cached(cfd_world)
    abl = climate.compute_climate_cached(abl_world)
    land = ~cfd.is_ocean

    land_biome_agreement = np.mean(cfd.biome_ids[land] == abl.biome_ids[land])
    assert land_biome_agreement > 0.78, land_biome_agreement

    precip_rel_rms = np.sqrt(np.mean((cfd.precipitation_mm - abl.precipitation_mm) ** 2)) / np.sqrt(
        np.mean(cfd.precipitation_mm ** 2)
    )
    assert precip_rel_rms < 0.25, precip_rel_rms

    # The CFD state on the diagnostic world never advanced -- _advance_fluid_dynamics is a
    # no-op in that mode.
    assert abl_world.atmosphere_cfd_state.elapsed_seconds == 0.0

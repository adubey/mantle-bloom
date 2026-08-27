import numpy as np
from app import climate, geometry
from app.world import generate_world, step_world


def test_generate_world_starts_plates_at_rest():
    # Motion is torque-driven (see torque.py) rather than fit from the mantle flow field at
    # generation time -- a freshly generated world's plates haven't yet resolved any torque
    # balance, so they start at rest and only acquire an omega once step_world's own
    # Plate.shift runs.
    world = generate_world(seed=10, num_plates=8)
    rates = [np.linalg.norm(p.omega) for p in world.plates]
    assert all(r == 0.0 for r in rates)


def test_step_world_gives_plates_nonzero_omega():
    world = generate_world(seed=10, num_plates=8)
    step_world(world, years=1_000_000)
    rates = [np.linalg.norm(p.omega) for p in world.plates]
    assert any(r > 0 for r in rates)


def test_different_plates_generally_acquire_different_omegas_after_a_step():
    world = generate_world(seed=14, num_plates=10)
    step_world(world, years=1_000_000)
    omegas = np.array([p.omega for p in world.plates])
    # Not every plate should have picked up an identical rotation from a spatially-varying
    # mantle flow field.
    assert not np.allclose(omegas, omegas[0], atol=1e-12)


def test_generate_world_logs_a_generation_event():
    world = generate_world(seed=15, num_plates=8, continental_fraction=3 / 8)
    assert len(world.events) == 1
    elapsed, message = world.events[0]
    assert elapsed == 0.0
    assert "8 plates" in message and "3 continental" in message


def test_step_world_at_doubled_climate_density_does_not_crash_and_uses_the_finer_grid():
    # End-to-end: erosion.py (every step) and main.py's controls route both compute climate
    # against world.climate_density -- this confirms that actually happens, not just that
    # climate.compute_climate itself accepts a height/width override.
    world = generate_world(seed=16, num_plates=8, continental_fraction=0.5, land_fraction=0.35, climate_density=2.0)
    step_world(world, years=1_000_000)
    assert world.climate_cache is not None
    assert world.climate_cache.elevation_m.shape == climate.grid_dimensions(2.0)

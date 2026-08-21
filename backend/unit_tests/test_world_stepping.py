import numpy as np
from app import geometry
from app.world import generate_world, step_world


def test_generate_world_gives_plates_nonzero_omega():
    world = generate_world(seed=10, num_plates=8)
    rates = [np.linalg.norm(p.omega) for p in world.plates]
    assert any(r > 0 for r in rates)


def test_different_plates_generally_acquire_different_omegas():
    world = generate_world(seed=14, num_plates=10)
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

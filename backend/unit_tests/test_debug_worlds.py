"""debug_worlds.py -- the "Debugging Worlds" Generate World tab's scripted plate scenarios and
pinned-omega motion mechanism."""

import numpy as np
import pytest
from app import debug_worlds, persistence
from app.world import step_world


def test_generate_debug_world_unknown_scenario_raises():
    with pytest.raises(ValueError):
        debug_worlds.generate_debug_world("not_a_real_scenario")


@pytest.mark.parametrize("scenario", sorted(debug_worlds.DEBUG_SCENARIOS))
def test_every_scenario_builds_plates_with_nodes_and_pinned_omegas(scenario):
    world = debug_worlds.generate_debug_world(scenario, seed=1)
    assert len(world.plates) > 0
    assert all(p.node_count() > 0 for p in world.plates)
    assert world.debug_diagnostics is True
    assert set(world.pinned_omegas) == {p.plate_id for p in world.plates}
    for omega in world.pinned_omegas.values():
        assert np.linalg.norm(omega) > 0.0


@pytest.mark.parametrize("scenario", sorted(debug_worlds.DEBUG_SCENARIOS))
def test_every_scenario_steps_without_crashing(scenario):
    world = debug_worlds.generate_debug_world(scenario, seed=2)
    for _ in range(3):
        step_world(world, years=1_000_000)
    assert world.elapsed_years == 3_000_000.0


def test_pinned_omega_overrides_what_the_real_torque_balance_would_produce():
    """The core correctness property of the override: after a step, every plate's own omega
    is exactly its pinned value, not whatever torque.shift_plate's force balance would have
    computed against this scenario's own (irrelevant, empty) mantle_centers."""
    world = debug_worlds.generate_debug_world("two_plate_divergent", seed=1)
    step_world(world, years=1_000_000)
    for plate in world.plates:
        assert np.allclose(plate.omega, world.pinned_omegas[plate.plate_id])


def test_two_plate_divergent_scenario_actually_diverges():
    """Sanity check that the scenario's own pinned omegas do what their name says: the two
    plates' seed points move apart step over step, not together or in lockstep."""
    from app import geometry

    world = debug_worlds.generate_debug_world("two_plate_divergent", seed=1)
    plate_a, plate_b = world.plates
    seed_a = geometry.to_world(plate_a.frame, np.array([1.0, 0.0, 0.0]))
    seed_b = geometry.to_world(plate_b.frame, np.array([1.0, 0.0, 0.0]))
    before = geometry.angular_distance(seed_a, seed_b)

    for _ in range(3):
        step_world(world, years=1_000_000)

    seed_a_after = geometry.to_world(plate_a.frame, np.array([1.0, 0.0, 0.0]))
    seed_b_after = geometry.to_world(plate_b.frame, np.array([1.0, 0.0, 0.0]))
    after = geometry.angular_distance(seed_a_after, seed_b_after)
    assert after > before


def test_triple_junction_mixed_scenario_exercises_fill_corner_notch():
    """This is the exact scenario `_fill_corner_notch`'s own docstring calls out (a triple
    junction with mixed divergent/convergent legs) -- confirm the debug world actually drives
    real corner-notch activity, not just ordinary end-growth, so it's a genuine reproduction
    of the problem this whole diagnostic suite exists to investigate."""
    world = debug_worlds.generate_debug_world("triple_junction_mixed", seed=1)
    for _ in range(10):
        step_world(world, years=1_000_000)

    outcomes = {entry["outcome"] for entry in world.corner_notch_log}
    assert "claimed" in outcomes, "expected the notch-filler to have claimed real nodes at some point"


def test_pinned_omegas_backfilled_on_load_of_an_older_save():
    world = debug_worlds.generate_debug_world("two_plate_divergent", seed=1)
    del world.__dict__["pinned_omegas"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.pinned_omegas == {}

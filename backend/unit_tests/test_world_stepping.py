import numpy as np
from app import climate, geometry, magma_transport
from app.world import generate_world, step_world


def test_generate_world_starts_plates_at_rest():
    # Motion is torque-driven (see torque.py) rather than fit from the mantle flow field at
    # generation time -- a freshly generated world's plates haven't yet resolved any torque
    # balance, so they start at rest and only acquire an omega once step_world's own
    # Plate.shift runs.
    world = generate_world(seed=10, num_plates=8)
    rates = [np.linalg.norm(p.omega) for p in world.plates]
    assert all(r == 0.0 for r in rates)


def test_generate_world_records_baseline_stats_without_leaking_a_climate_cache():
    # finish_generation's own record_stats(force=True) computes stats.compute_stats once for
    # the baseline stats_history entry, but must not leave world.climate_cache populated as a
    # side effect of that -- other code (coastline.py, climate.py's vegetation-transpiration
    # source) reads climate_cache is None as "no step has run yet."
    world = generate_world(seed=10, num_plates=8)
    assert len(world.stats_history) == 1
    assert world.stats_history[0]["elapsed_years"] == 0.0
    assert world.climate_cache is None


def test_step_world_skips_recording_stats_when_climate_is_off():
    # World.record_stats(force=False) (step_world's own call) must not force a fresh climate
    # computation purely to snapshot stats -- that would defeat simulate_climate_biomes=False
    # existing to run tectonics-only steps fast. The elapsed_years=0 baseline from generation
    # is still there (finish_generation always forces it), just nothing gets appended on top.
    world = generate_world(seed=10, num_plates=8)
    world.simulate_climate_biomes = False
    assert len(world.stats_history) == 1
    step_world(world, years=1_000_000)
    assert world.climate_cache is None
    assert len(world.stats_history) == 1


def test_step_world_uses_selected_magma_transport_k(monkeypatch):
    world = generate_world(seed=10, num_plates=6, node_density=0.5, climate_density=0.5)
    world.simulate_climate_biomes = False
    world.steps_taken = magma_transport.MAGMA_TRANSPORT_INTERVAL_STEPS - 1
    world.magma_transport_k = 128
    seen = []

    def capture_transport(_world, _banked_myr, max_destinations_per_parcel=None):
        seen.append(max_destinations_per_parcel)
        return []

    monkeypatch.setattr(magma_transport, "run_magma_transport", capture_transport)
    step_world(world, years=100_000)
    assert seen == [128]


def test_step_world_gives_plates_nonzero_omega():
    world = generate_world(seed=10, num_plates=8)
    step_world(world, years=1_000_000)
    rates = [np.linalg.norm(p.omega) for p in world.plates]
    assert any(r > 0 for r in rates)


def test_hydrology_cache_step_freezes_when_climate_toggled_off():
    # World.hydrology_cache_step (see its own docstring, GitHub issue #121) is how
    # stats.compute_stats tells "hydrology_cache reflects this step" apart from "hydrology_
    # cache is frozen from whenever simulate_climate_biomes was last on" -- it must track
    # World.steps_taken exactly while climate is on, then stop advancing the instant it's
    # toggled off, even though tectonics (and steps_taken itself) keep moving.
    world = generate_world(seed=10, num_plates=8)
    step_world(world, years=1_000_000)
    assert world.hydrology_cache is not None
    assert world.hydrology_cache_step == world.steps_taken

    world.simulate_climate_biomes = False
    step_world(world, years=1_000_000)
    assert world.hydrology_cache is not None  # last-good value, not cleared
    assert world.hydrology_cache_step != world.steps_taken


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


def test_step_world_advances_the_step_counter():
    # steps_taken drives merge_split.py's defragmentation cadence (see DEFRAG_INTERVAL_STEPS)
    # -- one increment per step_world call regardless of years, since that pass tracks
    # accumulated topology drift, not elapsed time.
    world = generate_world(seed=10, num_plates=8)
    assert world.steps_taken == 0
    step_world(world, years=1_000_000)
    step_world(world, years=50_000)
    assert world.steps_taken == 2


def test_step_world_at_doubled_climate_density_does_not_crash_and_uses_the_finer_grid():
    # End-to-end: erosion.py (every step) and main.py's controls route both compute climate
    # against world.climate_density -- this confirms that actually happens, not just that
    # climate.compute_climate itself accepts a height/width override.
    world = generate_world(seed=16, num_plates=8, continental_fraction=0.5, land_fraction=0.35, climate_density=2.0)
    step_world(world, years=1_000_000)
    assert world.climate_cache is not None
    assert world.climate_cache.elevation_m.shape == climate.grid_dimensions(2.0)


def test_pinned_omega_overrides_the_real_torque_balance():
    """World.pinned_omegas (the "Debugging Worlds" tab's scripted-motion mechanism) should
    make LithospherePlate.shift use the pinned value verbatim, bypassing torque.shift_plate's
    own torque-balance recompute entirely -- not just happen to produce the same result."""
    world = generate_world(seed=10, num_plates=8)
    pinned_plate = world.plates[0]
    pinned_omega = np.array([0.0, 0.0, 0.05])
    world.pinned_omegas[pinned_plate.plate_id] = pinned_omega

    step_world(world, years=1_000_000)

    assert np.allclose(pinned_plate.omega, pinned_omega)
    # Every other plate still goes through the real torque balance and (almost certainly)
    # doesn't land on this exact hand-picked value.
    others = [p for p in world.plates if p.plate_id != pinned_plate.plate_id]
    assert any(not np.allclose(p.omega, pinned_omega) for p in others)


def test_pinned_omegas_backfilled_on_load_of_an_older_save():
    from app import persistence

    world = generate_world(seed=10, num_plates=4)
    del world.__dict__["pinned_omegas"]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.pinned_omegas == {}


def test_step_world_reconciles_gap_tracks_alongside_fill_gaps():
    from app import gaps

    world = generate_world(seed=10, num_plates=8)
    assert world.gap_tracks == []
    for _ in range(gaps.GAP_FILL_INTERVAL_STEPS):
        step_world(world, years=1_000_000)
    # Ordinary boundary-growth catch-up lag can leave small (sub-MIN_GAP_NODES) transient
    # slivers even on a healthy world -- GAP_AGE_MIN_CLUSTER_NODES is deliberately far below
    # MIN_GAP_NODES so the tracker surfaces those too (see gaps.py's own module docstring on
    # "ordinary catch-up lag" vs a genuine void). The real health signal is that none of them
    # are big enough for fill_gaps to react to.
    assert isinstance(world.gap_tracks, list)
    assert all(track.node_count < gaps.MIN_GAP_NODES for track in world.gap_tracks)


def test_record_removed_points_is_a_noop_for_an_empty_array():
    world = generate_world(seed=10, num_plates=4)
    world.record_removed_points(np.zeros((0, 3)), plate_id=0)
    assert world.removed_points_log == []


def test_record_removed_points_appends_with_current_elapsed_years():
    from app.world import World

    world = World(seed=0, plates=[], mantle_centers=[])
    world.elapsed_years = 5_000_000.0
    points = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    world.record_removed_points(points, plate_id=3)
    assert len(world.removed_points_log) == 2
    for point, removed_years, plate_id in world.removed_points_log:
        assert removed_years == 5_000_000.0
        assert plate_id == 3
    np.testing.assert_array_equal(np.array([p for p, _, _ in world.removed_points_log]), points)


def test_record_removed_points_caps_the_log_evicting_oldest_first():
    from app.world import MAX_REMOVED_POINTS_LOG, World

    world = World(seed=0, plates=[], mantle_centers=[])
    over_cap = MAX_REMOVED_POINTS_LOG + 50
    for i in range(over_cap):
        world.elapsed_years = float(i)
        world.record_removed_points(np.array([[float(i), 0.0, 0.0]]), plate_id=0)
    assert len(world.removed_points_log) == MAX_REMOVED_POINTS_LOG
    # Oldest entries (elapsed_years 0..49) were evicted; the most recent ones survive.
    surviving_years = {removed_years for _, removed_years, _ in world.removed_points_log}
    assert min(surviving_years) == float(over_cap - MAX_REMOVED_POINTS_LOG)
    assert max(surviving_years) == float(over_cap - 1)

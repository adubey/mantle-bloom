import numpy as np
from app.hydroclimate import climate
from app.torque import shift_plates_quasistatic
from app.world import PROGRESS_HEARTBEAT_INTERVAL_S, generate_world, generate_world_progress, step_world, step_world_progress


def test_generate_world_starts_plates_at_rest():
    # Motion is torque-driven (see torque.py) rather than fit from the mantle flow field at
    # generation time. Freshly generated plates acquire omega in the first coupled step.
    world = generate_world(seed=10, num_plates=8)
    rates = [np.linalg.norm(p.omega) for p in world.plates]
    assert all(r == 0.0 for r in rates)


def test_coupled_motion_returns_zero_displacement_for_empty_plates():
    world = generate_world(seed=10, num_plates=4)
    empty = world.plates[0]
    empty._keys = np.array([], dtype=np.int64)

    distances = shift_plates_quasistatic(world, years=1_000_000)

    assert distances[empty.plate_id] == 0.0


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


def test_step_world_progress_labels_phases_in_order():
    # The run overlay's phase text (#629): per-pass labels plus one label per tail
    # phase boundary, fractions non-decreasing to exactly 1.0.
    world = generate_world(seed=10, num_plates=8)
    progress = [m for m in step_world_progress(world, years=1_000_000) if m[0] == "progress"]
    fractions = [fraction for _, fraction, _ in progress]
    assert fractions[-1] == 1.0
    assert all(b >= a for a, b in zip(fractions, fractions[1:]))
    phases = [phase for _, _, phase in progress]
    assert phases[0] == "Moving plates"
    assert phases[-1] == "Balancing sea level"
    assert "Raising mantle plumes" in phases
    ordered = [
        "Moving plates",
        "Reshaping plates",
        "Merging and splitting plates",
        "Filling seafloor gaps",
        "Resolving glaciers",
        "Raising volcanoes",
        "Forming natural resources",
        "Balancing sea level",
    ]
    positions = [phases.index(phase) for phase in ordered]
    assert positions == sorted(positions)


def test_step_world_progress_with_everything_off_labels_only_the_final_pass():
    # Movement and climate off: no per-plate units and no tail phases, so the only
    # labels are the sea-level pass's own boundary and completion lines.
    world = generate_world(seed=10, num_plates=8)
    world.simulate_plate_movement = False
    world.simulate_climate_biomes = False
    assert list(step_world_progress(world, years=1_000_000)) == [
        ("progress", 0.0, "Balancing sea level"),
        ("progress", 1.0, "Balancing sea level"),
    ]


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
    """World.pinned_omegas prescribes one plate's rate while the other plates still solve."""
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


def test_quad_world_keeps_its_represented_area_between_gap_fill_intervals():
    """Issue #259: quad retreat removes cells every step, so gap repair has to keep up every
    step -- otherwise the area the surface represents sags ~0.5% of the sphere per step until
    the next interval and eustatic sea level jumps with it."""
    from app.sparse_quad_patch import PLANET_RADIUS_M

    world = generate_world(seed=3, node_density=0.5)
    sphere_m2 = 4.0 * np.pi * PLANET_RADIUS_M**2
    for _ in range(2):
        step_world(world, years=1_000_000)
        represented = sum(float(np.sum(p.surface_nodes().area_m2)) for p in world.plates if p.node_count())
        assert represented / sphere_m2 > 0.995


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


def test_step_world_ends_with_every_column_inside_its_caps(monkeypatch):
    """Issue #256: "Hc, Hm within caps" is a hard per-step invariant (docs/plate-surface-
    baseline.md section 3.1). A writer that leaves a column outside them -- simulated here by
    the last Hc/Hm writer in the step -- is caught by the end-of-step clamp, which books what
    it changed in `phase_budget`."""
    from app import lithosphere
    from app.hydroclimate import volcanism

    world = generate_world(seed=3, node_density=0.5)
    world.simulate_climate_biomes = False
    world.debug_diagnostics = True
    original = volcanism.apply_volcanic_activity

    def breach(world, years):
        original(world, years)
        plate = max(world.plates, key=lambda p: p.node_count())
        hm = plate.collect("mantle_lithosphere_thickness_m").copy()
        hc = plate.collect("crustal_thickness_m").copy()
        hm[0], hm[1], hc[2] = 1.0, 1.5 * lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M, 0.0
        plate.set_fields_on_plate(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm)

    monkeypatch.setattr(volcanism, "apply_volcanic_activity", breach)
    step_world(world, years=1_000_000)

    for plate in world.plates:
        hc = plate.collect("crustal_thickness_m")
        hm = plate.collect("mantle_lithosphere_thickness_m")
        assert np.all((hc >= lithosphere.MIN_CRUSTAL_THICKNESS_M) & (hc <= lithosphere.MAX_CRUSTAL_THICKNESS_M))
        assert np.all((hm >= lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M) & (hm <= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M))
    clamp = world.phase_budget["column_cap_clamp"]["scopes"]["all"]
    assert clamp["sum_hc_after"] > clamp["sum_hc_before"]
    assert clamp["sum_hm_after"] != clamp["sum_hm_before"]


def test_step_world_progress_sends_heartbeats_from_its_slow_tail(monkeypatch):
    # GitHub issue #611: the tail past the plate loops (topology, climate/erosion,
    # caps, sea level) has no per-plate loop of its own, so a slow step went silent.
    # An overdue tick must fire at each slow sub-boundary carrying the tail fraction --
    # while a fast step stays exactly as silent as before. Progress lines are tagged
    # ("progress", fraction, phase) (#629); heartbeats stay ("heartbeat", fraction).
    world = generate_world(seed=10, num_plates=8)
    now = [0.0]

    def _fast_clock():
        now[0] += PROGRESS_HEARTBEAT_INTERVAL_S + 1.0
        return now[0]

    monkeypatch.setattr("app.world._monotonic", _fast_clock)
    messages = list(step_world_progress(world, 1_000_000))
    progress = [m for m in messages if m[0] == "progress"]
    heartbeats = [m for m in messages if m[0] == "heartbeat"]
    assert len(progress) == 2 * 8 + 8  # per-plate units plus one labeled line per tail phase
    fractions = [m[1] for m in progress]
    assert fractions == sorted(fractions) and fractions[-1] == 1.0
    assert len(heartbeats) == 4  # topology, erosion, geology, caps
    assert all(m[1] == fractions[-2] for m in heartbeats)  # tail fraction, no progress implied


def test_step_world_progress_stays_silent_between_ticks_when_fast():
    # The mirror: a fast step's own per-plate ticks already prove life, so no
    # heartbeat should pad the stream (progress lines still carry their labels, #629).
    world = generate_world(seed=10, num_plates=8)
    messages = list(step_world_progress(world, 1_000_000))
    assert all(m[0] == "progress" for m in messages)
    assert len(messages) == 2 * 8 + 8
    assert "Moving magma" not in [m[2] for m in messages]  # interval-gated: step 1 isn't its turn


def test_step_world_progress_labels_interval_gated_tail_passes():
    # Interval-gated tail passes label only on their own steps (#629): with steps_taken
    # at 3, the increment to 4 fires the magma-transport interval.
    world = generate_world(seed=10, num_plates=8)
    world.steps_taken = 3
    progress = [m for m in step_world_progress(world, years=1_000_000) if m[0] == "progress"]
    assert "Moving magma" in [phase for _, _, phase in progress]


def test_generate_world_progress_sends_heartbeats_from_its_slow_tail(monkeypatch):
    # Same #611 reasoning for world creation: the relief/soil/finish-bootstrap
    # sub-calls subdivide the final third so a multi-minute generation reads as alive.
    now = [0.0]

    def _fast_clock():
        now[0] += PROGRESS_HEARTBEAT_INTERVAL_S + 1.0
        return now[0]

    monkeypatch.setattr("app.world._monotonic", _fast_clock)
    messages = list(generate_world_progress(seed=10, num_plates=8))
    assert [m[0] for m in messages] == ["progress", "progress", "heartbeat", "heartbeat", "heartbeat", "progress", "done"]
    assert [m[1] for m in messages if m[0] != "done"] == [1 / 3, 2 / 3, 3 / 4, 5 / 6, 11 / 12, 1.0]

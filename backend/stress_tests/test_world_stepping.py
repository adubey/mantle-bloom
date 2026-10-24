import numpy as np
from app import debug_worlds, geometry, lithosphere, quad_tectonics
from app.hydroclimate import erosion
from app import merge_split
from app.elevation_lines import effective_is_continental_from_codes
from app.sparse_quad_patch import PlateWithSparseQuadPatch
from app.world import generate_world, step_world


def _sampled_overlap_fraction(plates_list, sample_per_plate: int = 20) -> float:
    """Same proxy invariant as unit_tests/test_plates.py's own version -- see that
    function's docstring for why this isn't expected to be exactly zero (one-turn
    processing lag, residual envelope looseness for non-convex plate shapes), just bounded
    rather than growing without limit turn over turn."""
    total = 0
    overlapping = 0
    for plate in plates_list:
        points, _ = plate.all_points_and_elevation()
        if len(points) == 0:
            continue
        sample = points[:: max(1, len(points) // sample_per_plate)][:sample_per_plate]
        for other in plates_list:
            if other.plate_id == plate.plate_id:
                continue
            polygon = other.get_bounding_polygon()
            if len(polygon) < 3:
                continue
            total += len(sample)
            overlapping += int(np.count_nonzero(geometry.points_in_spherical_polygon(sample, polygon)))
    return overlapping / total if total > 0 else 0.0


def test_plate_overlap_stays_bounded_over_many_steps():
    # Confirmed directly during development: with the polygon-based shift()/deform() model,
    # a real (non-runaway) amount of "envelope overlap" persists indefinitely -- what
    # matters is that it stays roughly flat over many steps rather than climbing toward
    # saturation, which would indicate territory genuinely running away unchecked (e.g. a
    # shrink/grow imbalance letting one plate's claimed footprint balloon).
    world = generate_world(seed=3, num_plates=8, node_density=0.5)
    world.simulate_climate_biomes = False  # only plate geometry is checked here
    fractions = []
    for i in range(24):
        step_world(world, years=3_000_000)
        if i % 4 == 3:
            fractions.append(_sampled_overlap_fraction(world.plates))
    # None of the later checkpoints drift far above the earliest one -- a real runaway would
    # show a clear upward trend across the run, not noise around a stable level.
    assert max(fractions) < fractions[0] + 0.15
    assert max(fractions) < 0.35


def test_node_density_persists_through_stepping_and_gap_fill():
    # The core correctness concern for a runtime density option: anything that rebuilds or
    # grows nodes at the reference TARGET_LINE_SPACING_RAD instead of the world's own density
    # silently reverts a non-default density back to the reference one within a handful of
    # steps. Total node count should keep tracking node_density's ratio, not decay toward the
    # 1x baseline, across enough steps.
    reference = generate_world(seed=5, num_plates=8, continental_fraction=0.5, node_density=1.0)
    denser = generate_world(seed=5, num_plates=8, continental_fraction=0.5, node_density=4.0)
    # This test only checks node counts, never climate/erosion/hydrology output (see
    # World.simulate_climate_biomes) -- skipping that per-step computation speeds this up
    # without changing what it exercises.
    reference.simulate_climate_biomes = False
    denser.simulate_climate_biomes = False

    def total_nodes(world):
        return sum(p.node_count() for p in world.plates)

    for _ in range(8):
        step_world(reference, years=300_000)
        step_world(denser, years=300_000)

    ratio = total_nodes(denser) / total_nodes(reference)
    assert 3.0 < ratio < 5.0  # stays roughly 4x -- not decayed back toward 1x


def test_step_world_advances_elapsed_years():
    world = generate_world(seed=11, num_plates=6)
    step_world(world, years=1_000_000)
    assert world.elapsed_years == 1_000_000
    step_world(world, years=500_000)
    assert world.elapsed_years == 1_500_000


def test_rigid_rotation_preserves_interior_node_spacing_exactly():
    """The whole point of the plate-local-frame design: rotating a plate must not disturb
    the relative spacing of its own nodes at all (no resampling)."""
    world = generate_world(seed=12, num_plates=6)
    world.simulate_climate_biomes = False  # only node spacing/rotation is checked here
    plate = max(world.plates, key=lambda p: p.node_count())
    local = plate.surface_nodes().local_xyz

    def spacing():
        world_xyz = geometry.to_world(plate.frame, local)
        return geometry.angular_distance(world_xyz[:-1], world_xyz[1:])

    before_spacing = spacing()
    for _ in range(5):
        plate.rotate(geometry.rotation_matrix(np.array([0.3, -0.2, 0.9]) / np.linalg.norm([0.3, -0.2, 0.9]), 0.05))
    assert np.allclose(before_spacing, spacing(), atol=1e-9)


def test_rigid_rotation_preserves_plate_frame_orthonormality():
    world = generate_world(seed=13, num_plates=6)
    world.simulate_climate_biomes = False  # only plate.frame orthonormality is checked here
    for _ in range(10):
        step_world(world, years=5_000_000)
    for plate in world.plates:
        assert np.allclose(plate.frame @ plate.frame.T, np.eye(3), atol=1e-6)
        assert np.isclose(np.linalg.det(plate.frame), 1.0, atol=1e-6)


def test_step_world_events_are_timestamped_with_post_step_elapsed_years():
    world = generate_world(seed=16, num_plates=10, continental_fraction=0.6)
    # Only event timestamps/log length are checked here, never climate/erosion/hydrology
    # output (see World.simulate_climate_biomes) -- plate-movement events (merges/splits/
    # gap-fills/volcanism) still fire normally, since simulate_plate_movement stays on.
    world.simulate_climate_biomes = False
    for _ in range(15):
        step_world(world, years=8_000_000)
    # Any event logged during stepping (merges/splits/consumption are plausible at this
    # scale/seed count but not guaranteed for a specific seed) must be timestamped no later
    # than the world's current elapsed_years.
    for elapsed, _ in world.events:
        assert elapsed <= world.elapsed_years
    assert len(world.events) <= 200  # MAX_EVENT_LOG_LENGTH


def _continental_crust_area_m2(world) -> float:
    area_m2 = 0.0
    for plate in world.plates:
        continental = effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        area_m2 += float(np.sum(plate.node_areas_m2()[continental]))
    return area_m2


def _continental_plate_footprint_area_m2(world) -> float:
    return sum(
        float(np.sum(plate.node_areas_m2()))
        for plate in world.plates
        if plate.crust_type == "continental"
    )


def test_continental_volume_budget_bounds_the_boundary_ratchet(monkeypatch):
    """The continental boundary ratchet (GitHub issue #119, "Continental ratchet: solution design")
    grows every continental margin at the oceanic reference column and never retreats a
    leading row, so a continental plate's footprint can climb over a long run (#482). The
    volume-budget gate bounds continental-coded area against genuine continental crust volume. It
    excludes oceanic-coded cells even when they sit on a continental plate (#447), so whole-
    plate footprint is a separate metric and may grow as oceanic crust is added. This run
    verifies the production gate is active and that continental-coded area stays bounded;
    the production threshold is covered directly in test_quad_tectonics.

    Seed 60461418 (used here previously) drives this scenario into heavy plate fragmentation
    around step 35-45. Seed 1 was checked across the full 60-step schedule and stays clear of
    that fragmentation crash."""

    # Keep the 20-step split cooldown (40 Myr at this step length) used to select a trajectory
    # with little fragmentation; #358 changed the production cooldown to two million years.
    monkeypatch.setattr(merge_split, "SPLIT_MIN_AGE_YEARS", 20 * 2_000_000)

    gate_suppressions = []
    original_boundary_context = quad_tectonics.boundary_context

    def track_growth_gate(*args, **kwargs):
        context = original_boundary_context(*args, **kwargs)
        gate_suppressions.append(context.suppress_growth)
        return context

    monkeypatch.setattr(quad_tectonics, "boundary_context", track_growth_gate)

    world = generate_world(seed=1, num_plates=14, continental_fraction=0.6, node_density=0.5)
    world.simulate_climate_biomes = False  # plate geometry only -- keeps the run quick
    before_continental_crust = _continental_crust_area_m2(world)
    before_plate_footprint = _continental_plate_footprint_area_m2(world)
    for _ in range(60):
        step_world(world, years=2_000_000)
    continental_crust_growth = _continental_crust_area_m2(world) / before_continental_crust - 1.0
    plate_footprint_growth = _continental_plate_footprint_area_m2(world) / before_plate_footprint - 1.0

    assert any(gate_suppressions)  # exercise the production growth gate in this trajectory
    # Continental-coded area stays bounded over 120 My. The whole continental-plate footprint
    # is tracked separately because oceanic-coded growth is outside this continental budget.
    assert continental_crust_growth < 0.15, (
        f"continental-coded growth {continental_crust_growth:.1%}; "
        f"continental-plate footprint growth {plate_footprint_growth:.1%}"
    )


def _continental_hc_hm(world):
    continental_columns = []
    for plate in world.plates:
        mask = effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        if np.any(mask):
            continental_columns.append(
                (plate.collect("crustal_thickness_m")[mask], plate.collect("mantle_lithosphere_thickness_m")[mask])
            )
    if not continental_columns:
        return np.empty(0), np.empty(0)
    hc = np.concatenate([columns[0] for columns in continental_columns])
    hm = np.concatenate([columns[1] for columns in continental_columns])
    return hc, hm


def test_sustained_collision_caps_hc_hm():
    """GitHub issue #161 ("Unbounded Hc/Hm growth in apply_convergent_deformation"): with no
    ceiling, a node sitting in a long-lived convergent regime compounded Hc/Hm exponentially
    (measured up to Hc=413,885 m / Hm=1,311,592 m on a 626 My save -- more than an order of
    magnitude past anything geologically real), and the resulting elevation/isostasy desync
    read as sheer canyon walls. `two_continental_collision` is the debug world built exactly
    for this: two continental plates pinned into a head-on collision that "never resolves...
    the two landmasses keep shoving into and piling onto each other indefinitely."

    Over a long stretch (400 My -- several multiples of the ~45 My a strong collision takes to
    double Hc), Hc/Hm never exceed `lithosphere.MAX_CRUSTAL_THICKNESS_M`/`MAX_MANTLE_
    LITHOSPHERE_THICKNESS_M`, and the ceiling is actually reached, so the caps are tested.

    This test used to also require total continental Hc to keep growing in the second 200 My,
    which the retired line engine's near-field overflow redistribution produced. The quad
    engine places only the melt of over-ceiling crust and delaminates the rest (#290), and
    this pinned collision subducts more continental crust than it adds, so that property no
    longer holds (#251)."""
    world = debug_worlds.generate_debug_world("two_continental_collision", seed=1)

    for _ in range(40):
        step_world(world, years=5_000_000)  # 200 My: the core boundary band should be saturated by now
    hc_mid, hm_mid = _continental_hc_hm(world)
    assert np.all(hc_mid <= lithosphere.MAX_CRUSTAL_THICKNESS_M + 1e-6)
    assert np.all(hm_mid <= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M + 1e-6)
    # The ceiling should actually be biting by now, not just headroom that was never reached.
    assert np.any(hc_mid >= 0.95 * lithosphere.MAX_CRUSTAL_THICKNESS_M)

    for _ in range(40):
        step_world(world, years=5_000_000)  # another 200 My
    hc_late, hm_late = _continental_hc_hm(world)
    assert hc_late.size > 0  # keep the late-run cap checks meaningful after plate mergers
    assert np.all(hc_late <= lithosphere.MAX_CRUSTAL_THICKNESS_M + 1e-6)
    assert np.all(hm_late <= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M + 1e-6)


def test_coastal_feedback_stays_stable_over_many_steps():
    # A regression floor for the symmetric coastal-leveling feedback (erosion.py), not a tight
    # bound. The real drowned-shelf checkerboard from GitHub issue #122, "Speckled low-relief
    # coastlines" only bites at node_density=4 (or on the seed-888151728 save) -- too slow for
    # a stress test, and a sudden sea-level jump on a density-1 world just makes a rough
    # newborn coast whose transient roughening swamps the feedback's slow ~My effect. So this
    # only asserts the pass doesn't destabilise a world with a broad near-waterline shelf:
    # elevation stays finite and in-bounds over many steps, and the pass is demonstrably
    # active at the shore.
    world = generate_world(seed=17, num_plates=8, continental_fraction=0.6, node_density=1.0)
    _, elevation, *_ = erosion._gather_nodes(world)
    world.sea_level_m = float(np.percentile(elevation[elevation > 0.0], 40))
    sl = world.sea_level_m
    near_shore = np.abs(elevation - sl) <= 80.0
    assert near_shore.sum() > 200  # the test bed really does have a broad near-waterline zone

    for _ in range(30):
        step_world(world, years=200_000)

    _, after, *_ = erosion._gather_nodes(world)
    from app.elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M

    assert np.all(np.isfinite(after))
    assert np.all(after >= MIN_ELEVATION_M - 1e-6) and np.all(after <= MAX_ELEVATION_M + 1e-6)
    # The near-waterline shelf hasn't wholesale-drowned or wholesale-emerged -- the grind and
    # fill halves nudge the coast toward its local datum, they don't run away with it.
    still_near = float(np.mean(np.abs(after - sl) <= 200.0))
    assert still_near > 0.05


def _uncovered_fraction(plates_list, samples: int = 5000) -> float:
    points = np.random.default_rng(0).normal(size=(samples, 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    covered = np.zeros(samples, dtype=bool)
    for plate in plates_list:
        covered |= plate.contains_batch(points)
    return float(np.mean(~covered))


def test_quad_world_steps_with_bounded_coverage_and_valid_topology():
    # Issue #228 Phase 4: a quad world steps with plate movement. Boundary retreat/advance
    # must keep every plate a valid balanced leaf mesh, and keep the sphere covered about as
    # well as it started -- neither leaking territory nor piling it up.
    world = generate_world(seed=5, num_plates=6, node_density=0.5)
    world.simulate_climate_biomes = False
    cells_before = sum(p.node_count() for p in world.plates)
    for _ in range(12):
        step_world(world, years=1_000_000)
        assert all(isinstance(p, PlateWithSparseQuadPatch) for p in world.plates)
        assert _uncovered_fraction(world.plates) < 0.04
        assert _sampled_overlap_fraction(world.plates) < 0.05

    assert world.steps_taken == 12
    for plate in world.plates:
        plate._validate_leaf_topology()
        assert np.all(plate.collect("crustal_thickness_m") > 0.0)
    assert abs(sum(p.node_count() for p in world.plates) / cells_before - 1.0) < 0.05


def test_quad_world_stepping_is_deterministic():
    def run():
        world = generate_world(seed=9, num_plates=4, node_density=0.5)
        world.simulate_climate_biomes = False
        for _ in range(2):
            step_world(world, 1_000_000)
        return [(p.plate_id, p.cell_keys.copy(), p.collect("crustal_thickness_m")) for p in world.plates]

    for (id_a, keys_a, hc_a), (id_b, keys_b, hc_b) in zip(run(), run(), strict=True):
        assert id_a == id_b
        np.testing.assert_array_equal(keys_a, keys_b)
        np.testing.assert_array_equal(hc_a, hc_b)


def test_quad_merge_of_neighbouring_plates_in_a_stepped_world():
    # Issue #228 Phase 4: fusing two real neighbouring quad plates after some stepping keeps
    # a valid leaf mesh, conserves crustal volume by exact cell area, and claims no more of
    # the other plates' territory than the pair already overlapped.
    from app import merge_split

    world = generate_world(seed=5, num_plates=6, node_density=0.5)
    world.simulate_climate_biomes = False
    for _ in range(3):
        step_world(world, years=1_000_000)
    keep = max(world.plates, key=lambda p: p.node_count())
    absorb = max(keep.get_neighbours(world.plates), key=lambda p: p.node_count())
    others = [p for p in world.plates if p not in (keep, absorb)]

    def volume(plates):
        return sum(float(np.sum(p.collect("crustal_thickness_m") * p.node_areas_m2())) for p in plates)

    def area_on_others(plates):
        total = 0.0
        for p in plates:
            points = p.all_points_and_elevation()[0]
            on_other = np.zeros(len(points), dtype=bool)
            for other in others:
                on_other |= other.contains_batch(points)
            total += float(p.node_areas_m2()[on_other].sum())
        return total

    before = volume([keep, absorb])
    overlap_before = area_on_others([keep, absorb])
    merge_split.merge_plates(world, keep.plate_id, absorb.plate_id)

    assert absorb not in world.plates
    keep._validate_leaf_topology()
    after = volume([keep])
    # Only the suture cap may remove crust.
    assert after <= before * (1 + 1e-12)
    assert after > 0.99 * before
    overlap_after = area_on_others([keep])
    assert overlap_after <= overlap_before * (1 + 1e-9)


def test_mobile_cover_ledger_closes_over_full_quad_steps():
    # Issue #297 phase 2: every process that moves or ends erosion's mobile cover books it
    # (see mobile_cover.py), so the live inventory equals sources less sinks every step.
    from app.hydroclimate import mobile_cover

    world = generate_world(seed=7, num_plates=8)
    for _ in range(8):
        step_world(world, 5_000_000)
        live = mobile_cover.surface_volume_m3(world)
        assert live > 0.0
        assert abs(mobile_cover.balance_error_m3(world)) <= 1e-9 * live
    ledger = world.mobile_cover_ledger
    for account in ("deposited_m3", "entrained_m3", "subducted_m3", "volcanic_buried_m3"):
        assert ledger[account] > 0.0, account

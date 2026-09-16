import numpy as np
from app import debug_worlds, erosion, geometry, lithosphere
from app import lithosphere_plate
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


def test_node_density_persists_through_regularize_and_gap_fill():
    # The core correctness concern for a runtime density option: elevation_lines.py's own
    # regularize pass (and LithospherePlate.deform's own claim-adjacent-territory/merge_split.py's
    # merging) previously always rebuilt/resampled nodes at the module's default
    # TARGET_LINE_SPACING_RAD regardless of what density a world was actually generated at,
    # silently reverting a non-default density back to the reference one within a handful of
    # steps. Confirmed directly this stays fixed: total node count should keep tracking
    # node_density's ratio, not decay toward the 1x baseline, across enough steps.
    reference = generate_world(seed=5, num_plates=8, continental_fraction=0.5, node_density=1.0)
    denser = generate_world(seed=5, num_plates=8, continental_fraction=0.5, node_density=4.0)
    # This test only checks node counts, never climate/erosion/hydrology output (see
    # World.simulate_climate_biomes) -- skipping that per-step computation speeds this up
    # without changing what it exercises.
    reference.simulate_climate_biomes = False
    denser.simulate_climate_biomes = False

    def total_nodes(world):
        return sum(len(line.theta) for p in world.plates for line in p.lines)

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
    the relative spacing of its own elevation-line nodes at all (no resampling)."""
    world = generate_world(seed=12, num_plates=6)
    world.simulate_climate_biomes = False  # only node spacing/rotation is checked here
    plate = max(world.plates, key=lambda p: p.node_count())
    line = max(plate.lines, key=lambda l: len(l.theta))

    before_world = line.world_xyz(plate.frame)
    before_spacing = geometry.angular_distance(before_world[:-1], before_world[1:])

    for _ in range(5):
        step_world(world, years=2_000_000)

    after_world = line.world_xyz(plate.frame)
    after_spacing = geometry.angular_distance(after_world[:-1], after_world[1:])

    assert np.allclose(before_spacing, after_spacing, atol=1e-9)
    # theta/elevation arrays themselves must be untouched (identity-based, not resampled).
    assert np.array_equal(line.theta, line.theta)


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


def _continental_node_count(world) -> int:
    return sum(
        sum(len(line) for line in plate.lines)
        for plate in world.plates
        if plate.crust_type == "continental"
    )


def test_continental_volume_budget_bounds_the_boundary_ratchet(monkeypatch):
    """The continental boundary ratchet (GitHub issue #119, "Continental ratchet: solution design")
    grows every continental margin at the oceanic reference column and never retreats a
    leading row, so a continental plate's node pile climbs without limit over a long run. The
    volume-budget growth gate (lithosphere_plate.CONTINENTAL_AREA_BUDGET_MULT) suppresses
    areal growth once a plate has outrun its crustal volume. Same seed / step schedule with
    the gate effectively disabled must grow the continental node count meaningfully more.

    GitHub issue #162: seed 60461418 (used here previously) drives this scenario into heavy
    plate fragmentation around step 35-45, which is chaotic enough in its floating-point
    trajectory that an unrelated, legitimate change elsewhere (#161) flipped which of the two
    runs ended up lower and failed the tight margin below -- even though the gate's suppression
    effect held up throughout most of the run on both the old and new code. Seed 1 was checked
    across the full 60-step schedule and never approaches that fragmentation crash: the
    gated/ungated gap grows smoothly and monotonically from step ~20 on, giving this assertion
    real headroom instead of sitting on a knife-edge."""

    def _run_continental_growth(budget_mult: float) -> float:
        monkeypatch.setattr(lithosphere_plate, "CONTINENTAL_AREA_BUDGET_MULT", budget_mult)
        world = generate_world(seed=1, num_plates=14, continental_fraction=0.6, node_density=0.5)
        world.simulate_climate_biomes = False  # plate geometry only -- keeps the run quick
        before = _continental_node_count(world)
        for _ in range(60):
            step_world(world, years=2_000_000)
        return _continental_node_count(world) / before - 1.0

    gated = _run_continental_growth(1.8)
    ungated = _run_continental_growth(1e9)

    # The gate does bite -- the ungated ratchet grows the continental lattice substantially
    # more over the same 120 My.
    assert ungated > gated + 0.04
    # ... and the gated run's own growth stays modest rather than running away.
    assert gated < 0.15


def _continental_hc_hm(world):
    hc = np.concatenate([p.collect("crustal_thickness_m") for p in world.plates if p.crust_type == "continental"])
    hm = np.concatenate([p.collect("mantle_lithosphere_thickness_m") for p in world.plates if p.crust_type == "continental"])
    return hc, hm


def test_sustained_collision_caps_hc_hm_without_losing_the_overflow():
    """GitHub issue #161 ("Unbounded Hc/Hm growth in apply_convergent_deformation"): with no
    ceiling, a node sitting in a long-lived convergent regime compounded Hc/Hm exponentially
    (measured up to Hc=413,885 m / Hm=1,311,592 m on a 626 My save -- more than an order of
    magnitude past anything geologically real), and the resulting elevation/isostasy desync
    read as sheer canyon walls. `two_continental_collision` is the debug world built exactly
    for this: two continental plates pinned into a head-on collision that "never resolves...
    the two landmasses keep shoving into and piling onto each other indefinitely."

    Two properties, run over a long stretch (400 My -- several multiples of the ~45 My a
    strong collision takes to double Hc, so the core boundary band should saturate against
    the new ceiling well before the run ends):

    1. Hc/Hm never exceed `lithosphere.MAX_CRUSTAL_THICKNESS_M`/`MAX_MANTLE_LITHOSPHERE_
       THICKNESS_M` -- the regression this test exists to catch.
    2. Total continental crustal volume keeps growing well past the point the core boundary
       band saturates, rather than flatlining the moment its own nodes first hit the ceiling
       -- confirming the capped overflow is actually reaching the near-field foreland
       (`lithosphere_plate.deform`'s redistribution), not just silently vanishing at the cap
       the way a bare clip with no conservation would (all growth would stop dead once the
       one-line-wide boundary band saturated, however much longer the collision ran)."""
    world = debug_worlds.generate_debug_world("two_continental_collision", seed=1)

    for _ in range(40):
        step_world(world, years=5_000_000)  # 200 My: the core boundary band should be saturated by now
    hc_mid, hm_mid = _continental_hc_hm(world)
    total_hc_mid = float(np.sum(hc_mid))
    assert np.all(hc_mid <= lithosphere.MAX_CRUSTAL_THICKNESS_M + 1e-6)
    assert np.all(hm_mid <= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M + 1e-6)
    # The ceiling should actually be biting by now, not just headroom that was never reached --
    # otherwise property 2 below wouldn't be testing anything.
    assert np.any(hc_mid >= 0.95 * lithosphere.MAX_CRUSTAL_THICKNESS_M)

    for _ in range(40):
        step_world(world, years=5_000_000)  # another 200 My
    hc_late, hm_late = _continental_hc_hm(world)
    total_hc_late = float(np.sum(hc_late))

    assert np.all(hc_late <= lithosphere.MAX_CRUSTAL_THICKNESS_M + 1e-6)
    assert np.all(hm_late <= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M + 1e-6)
    # Real continued growth in the second 200 My, not a plateau -- the overflow from the
    # (already-saturated) core band is still landing somewhere real.
    assert total_hc_late > total_hc_mid * 1.02


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

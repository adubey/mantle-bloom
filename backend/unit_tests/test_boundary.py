import numpy as np
from app import boundary, geometry, mantle
from app.elevation_lines import ERUPTION_ELEVATION_M, line_spacing_rad
from app.plates import (
    DIVERGENT_RIFT_TARGET_M,
    ElevationLine,
    PlateWithLines,
    _band_intensity,
    _far_field_intensity,
)
from app.world import World


def _patch(plate_id, crust_type, seed_xyz, phi_offsets, theta_values, base_elevation, omega=(0.0, 0.0, 0.0)):
    """A small multi-row plate at `seed_xyz`'s own frame -- multiple lines (so it has a real
    2D bounding polygon, unlike a single-line plate, whose two-point "outline" never has
    enough vertices for point_in_spherical_polygon to register containment at all). Every
    row shares the same `theta_values` for simplicity."""
    frame = geometry.plate_frame_from_seed(np.asarray(seed_xyz, dtype=float))
    theta = np.asarray(theta_values, dtype=float)
    lines = [
        ElevationLine(phi=phi, theta=theta.copy(), elevation=np.full(len(theta), base_elevation))
        for phi in phi_offsets
    ]
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type=crust_type, omega=np.asarray(omega, dtype=float), lines=lines)


def test_closing_rate_positive_when_approaching():
    # Two points near the equator, close together in longitude.
    p_self = np.array([np.cos(0.0), np.sin(0.0), 0.0])
    p_neighbor = np.array([np.cos(0.01), np.sin(0.01), 0.0])

    # self plate spins +z (eastward motion at the equator) -> moves toward neighbor (east of it)
    approaching_omega = np.array([0.0, 0.0, 1.0])
    still_omega = np.zeros(3)
    closing = boundary.closing_rate(p_self[None, :], approaching_omega, still_omega, p_neighbor[None, :])
    assert closing[0] > 0

    # self plate spins -z (westward motion) -> moves away from neighbor
    receding_omega = np.array([0.0, 0.0, -1.0])
    closing = boundary.closing_rate(p_self[None, :], receding_omega, still_omega, p_neighbor[None, :])
    assert closing[0] < 0


def test_band_intensity_zero_outside_band_peaks_at_midpoint():
    dist = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    result = _band_intensity(dist, inner=0.1, outer=0.3)
    assert result[0] == 0.0  # below inner edge
    assert np.isclose(result[1], 0.0)  # exactly at inner edge
    assert result[2] == 1.0  # at the band's midpoint
    assert np.isclose(result[3], 0.0)  # exactly at outer edge
    assert result[4] == 0.0  # beyond outer edge
    assert np.all((result >= 0.0) & (result <= 1.0))


def test_far_field_intensity_zero_below_inner_ramps_to_zero_at_outer():
    dist = np.array([0.0, 0.05, 0.1, 0.2, 0.3])
    result = _far_field_intensity(dist, inner=0.1, outer=0.3)
    assert result[0] == 0.0  # well below inner
    assert result[1] == 0.0  # still below inner
    assert result[2] == 1.0  # right at inner edge -- full intensity
    assert np.isclose(result[3], 0.5)  # midway between inner and outer
    assert np.isclose(result[4], 0.0)  # at outer edge
    assert np.all((result >= 0.0) & (result <= 1.0))


def test_shift_rotates_and_reports_max_displacement():
    rate = mantle.cm_per_yr_to_rad_per_yr(5.0)
    spacing = line_spacing_rad(1.0)
    plate = _patch(
        0, "continental", [1.0, 0.0, 0.0],
        phi_offsets=[-spacing, 0.0, spacing],
        theta_values=[-2 * spacing, -spacing, 0.0, spacing, 2 * spacing],
        base_elevation=500.0,
        omega=[0.0, 0.0, rate],
    )
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)
    old_frame = plate.frame.copy()

    d = plate.shift(world, years=1_000_000)

    assert not np.allclose(plate.frame, old_frame)  # actually rotated
    assert d > 0.0
    # D is the greatest angular distance any node moved -- bounded above by the rotation
    # angle itself (achieved only by a node exactly 90 degrees from the axis).
    assert d <= mantle.MAX_PLATE_RATE * 1_000_000 * 1.001


def test_deform_uplifts_collision_zone_without_shrinking_continental_territory():
    # Continental crust doesn't subduct in reality -- a colliding continent crumples
    # (thickens/uplifts) in place rather than having its own territory deleted node by node,
    # so a continental self-plate's contested nodes must survive a collision untouched aside
    # from the elevation change (deform()'s own `shrinkable_all` is all-False for a
    # continental self-plate; only a genuinely subducting oceanic self-plate's own contested
    # nodes -- see test_deform_shrinks_oceanic_trench_on_subduction below -- still shrink).
    # merge_split.py's own slow (50-100 Myr) fusion is what actually resolves the overlap.
    spacing = line_spacing_rad(1.0)
    # Same frame (same seed) as plate_b -- so local theta directly overlaps. plate_a spans
    # theta in [-7*spacing, 0]; plate_b spans [-3*spacing, 4*spacing] -- a genuine overlap
    # in [-3*spacing, 0], not merely two patches touching at a single shared edge point
    # (which -- confirmed directly -- point_in_spherical_polygon's winding-number test
    # doesn't reliably register as "contained" at all).
    theta_a = np.arange(8) * spacing - 7 * spacing
    plate_a = _patch(0, "continental", [1.0, 0.0, 0.0], phi_offsets=[-spacing, 0.0, spacing], theta_values=theta_a, base_elevation=500.0)

    theta_b = np.arange(8) * spacing - 3 * spacing
    plate_b = _patch(1, "continental", [1.0, 0.0, 0.0], phi_offsets=[-spacing, 0.0, spacing], theta_values=theta_b, base_elevation=500.0)

    world = World(seed=0, plates=[plate_a, plate_b], mantle_centers=[], node_density=1.0)
    before_counts = [len(line.theta) for line in plate_a.lines[:3]]
    plate_a.deform(world, [plate_b], years=1_000_000, max_distance=10 * spacing)

    # Only the three original rows (phi = -spacing, 0, spacing) actually overlapped plate_b.
    # None of plate_a's own contested nodes were deleted -- every original node is still
    # there (deform() may still append brand-new nodes at each line's *other*, uncontested
    # end in this same call, so this checks "at least as many," not "exactly as many").
    for line, before in zip(plate_a.lines[:3], before_counts):
        assert len(line.theta) >= before
        assert np.all(np.diff(line.theta) > 0)  # stays a contiguous, sorted span

    # The contested nodes (originally at theta=0, deep inside plate_b's [-3s, 4s] span)
    # picked up collision uplift -- continental-continental convergence, not a trench --
    # rather than being removed.
    contested_node = plate_a.lines[1].elevation[plate_a.lines[1].theta == 0.0]
    assert len(contested_node) == 1
    assert contested_node[0] > 500.0


def test_deform_shrinks_oceanic_trench_on_subduction():
    # The crust-type asymmetry from the test above cuts the other way for a genuinely
    # subducting oceanic self-plate: its own contested end still retreats (the slab's own
    # trench self-destructing), exactly as before this fix.
    spacing = line_spacing_rad(1.0)
    theta_a = np.arange(8) * spacing - 7 * spacing
    plate_a = _patch(0, "oceanic", [1.0, 0.0, 0.0], phi_offsets=[-spacing, 0.0, spacing], theta_values=theta_a, base_elevation=-3000.0)

    theta_b = np.arange(8) * spacing - 3 * spacing
    plate_b = _patch(1, "continental", [1.0, 0.0, 0.0], phi_offsets=[-spacing, 0.0, spacing], theta_values=theta_b, base_elevation=500.0)

    world = World(seed=0, plates=[plate_a, plate_b], mantle_centers=[], node_density=1.0)
    plate_a.deform(world, [plate_b], years=1_000_000, max_distance=10 * spacing)

    for line in plate_a.lines[:3]:
        # The contested end (originally at theta=0, deep inside plate_b's [-3s, 4s] span)
        # retreated past the overlap entirely, back to plate_b's own low edge.
        assert line.theta.max() <= -3 * spacing + 1e-9
        assert np.all(np.diff(line.theta) > 0)  # stays a contiguous, sorted span


def test_deform_grows_isolated_plate_with_no_neighbours():
    spacing = line_spacing_rad(1.0)
    plate = _patch(0, "oceanic", [1.0, 0.0, 0.0], phi_offsets=[0.0], theta_values=[0.0, spacing, 2 * spacing], base_elevation=-3000.0)
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)
    before = len(plate.lines[0].theta)

    plate.deform(world, [], years=1_000_000, max_distance=10 * spacing)

    grown_line = plate.lines[0]
    assert len(grown_line.theta) > before  # nothing to contest -> both ends grow
    assert np.all(np.diff(grown_line.theta) > 0)


def test_deform_growth_occasionally_spawns_a_volcano():
    # Growth coming back volcanic is a per-growth-event probability roll
    # (STRETCH_VOLCANO_PROBABILITY), not a deterministic threshold on any single call (see
    # that constant's own comment for why a threshold-based design didn't work: sampled
    # against a real running simulation, growth essentially always adds exactly 1 node per
    # call, so no per-call node-count threshold above 1 is ever reachable at realistic step
    # sizes/plate rates). The RNG is keyed by (seed, elapsed_years, plate_id, line_index,
    # end_tag), so sweeping plate_id across many otherwise-identical isolated single-node
    # growth events is a deterministic way to sample many independent rolls -- with
    # STRETCH_VOLCANO_PROBABILITY == 0.02, at least one hit in 200 trials is a near-certainty
    # (1 - 0.98**200 ~= 98%), not a flaky assumption.
    spacing = line_spacing_rad(1.0)
    found_a_volcano = False
    for plate_id in range(200):
        plate = _patch(plate_id, "continental", [1.0, 0.0, 0.0], phi_offsets=[0.0], theta_values=[0.0, spacing], base_elevation=200.0)
        world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)
        plate.deform(world, [], years=1_000_000, max_distance=5 * spacing)

        grown_line = plate.lines[0]
        new_nodes = grown_line.theta > spacing
        assert np.count_nonzero(new_nodes) > 0  # ordinary growth always happens
        if np.any(grown_line.is_volcano[new_nodes]):
            found_a_volcano = True
            assert np.all(grown_line.elevation[new_nodes] == DIVERGENT_RIFT_TARGET_M + ERUPTION_ELEVATION_M)
            assert np.all(grown_line.volcano_active_years_remaining[new_nodes] > 0)
            break
    assert found_a_volcano
    assert np.all(grown_line.volcano_active_years_remaining[new_nodes] > 0)

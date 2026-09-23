"""gap_fill_frontier.fill_gap_by_growing_plates: grow the plate(s) adjacent to a gap into it
node by node, extending an existing line where one's close enough or opening a new one where
none is, instead of always spawning a new plate or always emitting brand-new lines."""

import numpy as np

from app import gap_fill_frontier as gff
from app import geometry
from app.elevation_lines import line_spacing_rad
from app.lithosphere_plate import new_plate
from app.world import World

SPACING = line_spacing_rad(1.0)
FRAME = np.eye(3)


def _strip_plate(pid, theta_lo, theta_hi, phi_lim=0.2, crust_type="oceanic"):
    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        local = geometry.to_local(FRAME, world_pts)
        phi, theta = geometry.xyz_to_latlon(local)
        return (np.abs(phi) < phi_lim) & (theta >= theta_lo) & (theta < theta_hi)

    return new_plate(pid, FRAME, crust_type, SPACING, seed=1, is_owned=is_owned)


def _world(plates):
    return World(seed=1, plates=list(plates), next_plate_id=len(plates), node_density=1.0, mantle_centers=[])


def test_empty_claimants_or_gap_points_is_a_no_op():
    plate = _strip_plate(0, -0.1, 0.0)
    world = _world([plate])
    assert gff.fill_gap_by_growing_plates(world, np.zeros((0, 3)), [plate], SPACING) == {}
    assert gff.fill_gap_by_growing_plates(world, np.array([[1.0, 0.0, 0.0]]), [], SPACING) == {}


def test_extends_an_existing_line_instead_of_opening_a_new_one():
    """A gap point sitting right past an existing line's own end should grow that same line
    (line count unchanged), not proliferate a disjoint new one at the same row."""
    plate = _strip_plate(0, -0.1, 0.0)
    n_lines_before = len(plate.lines)
    world = _world([plate])

    gap_pts = []
    for line in plate.lines:
        if len(line) == 0:
            continue
        dtheta = SPACING / max(np.cos(line.phi), 1e-3)
        theta = float(line.theta[-1]) + dtheta
        gap_pts.append(geometry.to_world(FRAME, geometry.local_xyz(np.array([line.phi]), np.array([theta])))[0])
    gap_pts = np.array(gap_pts)

    added = gff.fill_gap_by_growing_plates(world, gap_pts, [plate], SPACING)
    assert added == {0: len(gap_pts)}
    assert len(plate.lines) == n_lines_before
    for line in plate.lines:
        assert np.all(np.diff(line.theta) > 0), "lines must stay theta-sorted after appending"


def test_a_run_of_new_points_coalesces_into_one_line_not_many():
    """A cluster of gap points far from any existing line, but contiguous with each other,
    should end up as one new multi-node line -- not one one-node line per point."""
    plate = _strip_plate(0, -0.1, 0.0)
    n_lines_before = len(plate.lines)
    world = _world([plate])

    # A row this plate has no line at all (well outside its own phi_lim), so every point in
    # it is a genuinely new claim regardless of how close together they are -- exercised
    # directly against _claim_row_points (the real per-row claiming code) rather than through
    # the full multi-hop orchestrator, since within *one* row a same-plate candidate close
    # enough to be hop-reachable from this plate's own cloud is, by construction (gap_tol ==
    # CONTIGUOUS_RUN_GAP_MULT * dtheta > DEFRAG_CONNECT_RADIUS_MULT * spacing_rad), always
    # also within stretch-adjacency range of an existing line where one exists -- so a
    # same-row "brand new, multi-node" case only arises at a row with no existing line yet.
    from app.lithosphere_plate import growth_seed_thickness
    from app import terrain_noise

    wp = gff._WorkingPlate(plate, SPACING)
    phi = 0.9  # outside phi_lim=0.2, so wp.rows has nothing at this row's key
    row_key = gff._pole_aligned_row_key(phi, SPACING)
    assert row_key not in wp.rows
    dtheta = SPACING / max(np.cos(phi), 1e-3)
    thetas = dtheta * np.arange(6)
    hc0, hm0 = growth_seed_thickness()
    texture = terrain_noise.FractalTexture(np.random.default_rng(0))

    added, consumed = gff._claim_row_points(world, wp, row_key, phi, gff.CONTIGUOUS_RUN_GAP_MULT * dtheta, thetas, hc0, hm0, hc0 * 0.1, texture, budget=6)
    assert (added, consumed) == (6, 6)
    new_lines = [line for line in wp.lines if len(line) > 0 and abs(line.phi - phi) < 1e-9]
    assert len(new_lines) == 1
    assert len(new_lines[0]) == 6
    assert len(wp.lines) == n_lines_before + 1


def test_new_nodes_are_typed_oceanic_underwater_and_erupted_as_volcanoes_above():
    """Every claimed node is a real magma eruption (is_volcano set), typed by whether it was
    above or below sea level at the moment it erupted -- see _erupt_melted_nodes."""
    from app.lithosphere_plate import growth_seed_thickness
    from app import terrain_noise

    plate = _strip_plate(0, -0.1, 0.0)
    world = _world([plate])
    wp = gff._WorkingPlate(plate, SPACING)
    phi = 0.9
    row_key = gff._pole_aligned_row_key(phi, SPACING)
    dtheta = SPACING / max(np.cos(phi), 1e-3)
    thetas = dtheta * np.arange(4)
    hc0, hm0 = growth_seed_thickness()
    texture = terrain_noise.FractalTexture(np.random.default_rng(0))

    gff._claim_row_points(world, wp, row_key, phi, gff.CONTIGUOUS_RUN_GAP_MULT * dtheta, thetas, hc0, hm0, hc0 * 0.1, texture, budget=4)
    new_line = next(line for line in wp.lines if len(line) > 0 and abs(line.phi - phi) < 1e-9)
    assert np.all(new_line.is_volcano)
    # Seeded from the oceanic reference column and thin -- every brand-new node here should
    # land underwater (the ordinary "new sea floor" case), not above sea level.
    assert np.all(new_line.elevation < 0.0)


def test_stretching_an_existing_line_thins_its_own_nearby_source_nodes():
    """Extending a *pre-existing continental* line draws its material down from that line's own
    nearest K_STRETCH_SOURCE_NODES end nodes -- their crustal thickness should measurably drop,
    unlike a brand-new line's own first node, which has no source nodes to thin at all. Stretch
    is continental-only (GitHub issue #216) -- see
    test_oceanic_claim_erupts_fresh_crust_instead_of_stretching for the oceanic case."""
    plate = _strip_plate(0, -0.1, 0.0, crust_type="continental")
    world = _world([plate])
    line = max((ln for ln in plate.lines if len(ln) > 3), key=len)
    hc_before = line.crustal_thickness_m[-gff.K_STRETCH_SOURCE_NODES :].copy()
    dtheta = SPACING / max(np.cos(line.phi), 1e-3)
    theta = float(line.theta[-1]) + dtheta
    gap_pt = geometry.to_world(FRAME, geometry.local_xyz(np.array([line.phi]), np.array([theta])))

    gff.fill_gap_by_growing_plates(world, gap_pt, [plate], SPACING)
    grown = next(ln for ln in plate.lines if abs(ln.phi - line.phi) < 1e-9)
    # The source nodes are wherever they ended up after re-sorting -- locate by matching the
    # original thetas rather than assuming a fixed offset.
    idx = np.searchsorted(grown.theta, line.theta[-gff.K_STRETCH_SOURCE_NODES :])
    assert np.all(grown.crustal_thickness_m[idx] < hc_before)


def test_oceanic_claim_erupts_fresh_crust_instead_of_stretching():
    """GitHub issue #216: an *oceanic* claimant extending a pre-existing line represents
    mid-ocean-ridge seafloor spreading (continuous fresh magma injection), not a fixed
    reservoir thinning to cover more area -- unlike the continental case above, its own nearest
    end nodes must be left untouched, and the new node itself should land at (near) the full
    oceanic reference column, not the stretched-down fraction _STRETCH_THIN_RATIO would give."""
    from app.lithosphere import REFERENCE_HC_OCEANIC_M

    plate = _strip_plate(0, -0.1, 0.0)  # default crust_type="oceanic"
    world = _world([plate])
    line = max((ln for ln in plate.lines if len(ln) > 3), key=len)
    hc_before = line.crustal_thickness_m[-gff.K_STRETCH_SOURCE_NODES :].copy()
    dtheta = SPACING / max(np.cos(line.phi), 1e-3)
    theta = float(line.theta[-1]) + dtheta
    gap_pt = geometry.to_world(FRAME, geometry.local_xyz(np.array([line.phi]), np.array([theta])))

    gff.fill_gap_by_growing_plates(world, gap_pt, [plate], SPACING)
    grown = next(ln for ln in plate.lines if abs(ln.phi - line.phi) < 1e-9)
    idx = np.searchsorted(grown.theta, line.theta[-gff.K_STRETCH_SOURCE_NODES :])
    assert np.allclose(grown.crustal_thickness_m[idx], hc_before), "source nodes must stay untouched, not thinned"
    new_idx = int(np.searchsorted(grown.theta, theta))
    # Seeded thin (same convention as every other fresh eruption -- _NEW_MAGMA_SEED_THIN_RATIO
    # guarantees it melts straight through), so it lands exactly at the full oceanic reference
    # column via _erupt_melted_nodes -- not the fractional _STRETCH_THIN_RATIO value stretching
    # from a full-thickness neighbour would have given instead.
    assert grown.crustal_thickness_m[new_idx] == REFERENCE_HC_OCEANIC_M


def test_multiple_claimants_assign_each_point_to_the_nearest_one():
    # A real (non-degenerate) gap strip between the two plates: a ends at theta=-0.02, b
    # starts at theta=0.02, leaving a 0.04-rad gap neither owns.
    plate_a = _strip_plate(0, -0.2, -0.02)
    plate_b = _strip_plate(1, 0.02, 0.2)
    world = _world([plate_a, plate_b])

    line_a = min((ln for ln in plate_a.lines if len(ln) > 0), key=lambda ln: abs(ln.phi))
    # theta=-0.01 sits 0.01 rad from a's own edge (-0.02) and 0.03 rad from b's (0.02) --
    # clearly closer to a.
    gap_pt = geometry.to_world(FRAME, geometry.local_xyz(np.array([line_a.phi]), np.array([-0.01])))

    added = gff.fill_gap_by_growing_plates(world, gap_pt, [plate_a, plate_b], SPACING)
    assert added == {0: 1}


def test_respects_max_nodes_budget_leaving_the_rest_uncovered():
    """10 different rows' worth of one-node extensions, each individually reachable in hop 1
    -- a max_nodes=3 budget should still claim only 3 of them, not silently drop or over-claim."""
    plate = _strip_plate(0, -0.1, 0.0, phi_lim=0.4)
    world = _world([plate])
    gap_pts = []
    for line in plate.lines:
        if len(line) == 0:
            continue
        dtheta = SPACING / max(np.cos(line.phi), 1e-3)
        theta = float(line.theta[-1]) + dtheta
        gap_pts.append(geometry.to_world(FRAME, geometry.local_xyz(np.array([line.phi]), np.array([theta])))[0])
    gap_pts = np.array(gap_pts)
    assert len(gap_pts) >= 3

    added = gff.fill_gap_by_growing_plates(world, gap_pts, [plate], SPACING, max_hops=5, max_nodes=3)
    assert added == {0: 3}

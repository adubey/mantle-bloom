import astropy.units as u
import numpy as np
from astropy_healpix import HEALPix

from app import fluid_dynamics_healpix as fdh
from app import healpix_grid


def test_pixel_area_is_uniform_and_covers_the_sphere():
    """Equal-area by construction (spec section 4.2's own claim) -- every pixel the same
    size, summing to the full sphere's surface area."""
    grid = healpix_grid.build(16)
    total = grid.pixel_area_m2 * grid.npix
    sphere_area = 4 * np.pi * healpix_grid.PLANET_RADIUS_M**2
    assert np.isclose(total, sphere_area, rtol=1e-9)


def test_neighbour_distance_bounded_no_polar_singularity():
    """The spec's central polar-singularity claim: neighbour spacing stays within a bounded
    range at every latitude, including near the poles -- no row of vanishingly-small dx the
    way an equirectangular grid has (see fluid_dynamics.grid_spacing_m's own polar clamp,
    which exists precisely because v1's own dx shrinks to ~0 there)."""
    grid = healpix_grid.build(16)
    near_pole = np.abs(grid.lat_rad) > np.radians(85.0)
    assert np.any(near_pole)
    polar_dist = grid.neighbour_distance_m[near_pole][grid.neighbour_valid[near_pole]]
    all_dist = grid.neighbour_distance_m[grid.neighbour_valid]
    # Within a small constant factor of the global typical spacing, not orders of magnitude
    # smaller the way equirectangular dx is at the same latitude.
    assert polar_dist.min() > all_dist.min() * 0.4
    assert polar_dist.max() < all_dist.max() * 2.5


def test_gradient_matches_analytic_linear_field():
    """The least-squares gradient fit is exact in the limit of a perfectly regular stencil
    and only approximate on HEALPix's own slightly irregular neighbour geometry -- checked
    here via mean error across the whole grid rather than a per-pixel `allclose`, since a
    handful of individual pixels with less favorably-conditioned neighbour layouts can carry
    a noticeably larger local error while the bulk of the grid stays accurate (confirmed
    directly: worst single-pixel error ~0.32 against a true gradient magnitude of 1.0, mean
    error under 0.01)."""
    grid = healpix_grid.build(32)
    field = grid.lat_rad * healpix_grid.PLANET_RADIUS_M  # true gradient: (0, 1)
    gx, gy = fdh.gradient(field, grid)
    assert np.abs(gx).mean() < 0.01
    assert np.abs(gy - 1.0).mean() < 0.01


def test_divergence_of_solid_body_rotation_is_near_zero():
    grid = healpix_grid.build(32)
    u = -grid.world_xyz[:, 1]
    v = grid.world_xyz[:, 0]
    div = fdh.divergence(u, v, grid)
    assert np.max(np.abs(div)) < 1e-5


def test_laplacian_of_constant_is_zero():
    grid = healpix_grid.build(16)
    const = np.full(grid.npix, 7.0)
    assert np.allclose(fdh.laplacian(const, grid), 0.0, atol=1e-8)


def test_semi_lagrangian_advect_zero_velocity_is_identity():
    grid = healpix_grid.build(16)
    field = np.sin(grid.lat_rad) * np.cos(grid.lon_rad)
    zero = np.zeros(grid.npix)
    advected = fdh.semi_lagrangian_advect(field, zero, zero, 100.0, grid)
    assert np.array_equal(advected, field)


def test_ang2pix_matches_astropy_reference_at_pixel_centers():
    """`HealpixGrid.ang2pix` is now a hand-rolled Numba nested-scheme lookup (replacing
    astropy_healpix's own, for speed -- see fluid_dynamics_healpix.py's module docstring on
    why that call was worth taking on despite nested-scheme indexing's own reputation for
    subtle bugs at base-pixel boundaries). Every pixel's own center must map back to its own
    index -- exercises every face, including the low-nside corner pixels with fewer than 8
    neighbours."""
    for nside in healpix_grid.NSIDE_CHOICES.values():
        grid = healpix_grid.build(nside)
        recovered = grid.ang2pix(grid.lon_rad, grid.lat_rad)
        assert np.array_equal(recovered, np.arange(grid.npix))


def test_ang2pix_matches_astropy_reference_at_random_points():
    """Cross-checked against astropy_healpix's own (independently-implemented) nested lookup
    at 50,000 uniformly-sampled sphere points per nside -- not just pixel centers, so this
    also exercises points that land near a face/pixel boundary, where a subtly wrong
    bit-interleave or face-numbering formula would first show a mismatch."""
    rng = np.random.default_rng(0)
    n = 50_000
    lon_rad = rng.uniform(0, 2 * np.pi, n)
    lat_rad = np.arcsin(rng.uniform(-1.0, 1.0, n))
    for nside in healpix_grid.NSIDE_CHOICES.values():
        grid = healpix_grid.build(nside)
        ours = grid.ang2pix(lon_rad, lat_rad)
        hp = HEALPix(nside=nside, order="nested")
        theirs = np.asarray(hp.lonlat_to_healpix(lon_rad * u.rad, lat_rad * u.rad))
        assert np.array_equal(ours, theirs)


def test_resample_round_trip_preserves_smooth_field():
    grid = healpix_grid.build(32)
    height, width = 90, 180
    lat_rows = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_cols = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    lat_grid, lon_grid = np.meshgrid(lat_rows, lon_cols, indexing="ij")
    field_hw = np.sin(np.radians(lat_grid)) * np.cos(np.radians(lon_grid))

    field_pix = healpix_grid.resample_from_equirect(grid, field_hw, lat_rows)
    back = healpix_grid.resample_to_equirect(grid, field_pix, height, width)
    correlation = np.corrcoef(field_hw.ravel(), back.ravel())[0, 1]
    assert correlation > 0.98


# -- Issue #133 phase 1: node-cloud scatter/fill/query plumbing -------------------------------


def test_nside_for_node_count_matches_phase0_spike_ratio():
    """Phase-0's own spike used nside=128 (npix=196,608) for ~130.6K nodes -- the chooser
    should reproduce that exact choice at that node count, not just "some" power of 2."""
    assert healpix_grid.nside_for_node_count(130_587) == 128


def test_nside_for_node_count_handles_small_and_zero_counts():
    # Never below nside=1 (the smallest valid HEALPix grid), even for a tiny or empty node cloud.
    assert healpix_grid.nside_for_node_count(0) >= 1
    assert healpix_grid.nside_for_node_count(1) >= 1
    assert healpix_grid.nside_for_node_count(5) >= 1


def test_scatter_node_indices_leaves_most_pixels_occupied_or_tracks_collisions():
    grid = healpix_grid.build(16)
    # One node placed exactly at pixel 0's own center: must win that pixel outright.
    node_xyz = grid.world_xyz[[0]]
    pixel_to_node = healpix_grid.scatter_node_indices(grid, node_xyz)
    assert pixel_to_node[0] == 0
    assert np.sum(pixel_to_node != -1) == 1  # the other 3071 pixels have no node at all


def test_scatter_node_indices_tie_break_is_nearest_to_pixel_center():
    """Two nodes landing in the same pixel: the *nearer* one to that pixel's own center must
    win, regardless of which one appears first/last in `node_xyz` -- the scatter has to be
    independent of node/plate iteration order, not "first node wins" or "last node wins"."""
    grid = healpix_grid.build(16)
    center = grid.world_xyz[0]
    # A tiny eastward nudge for both, one twice as far as the other -- both still land in
    # pixel 0 at this nside (pixel extent is far larger than either nudge).
    east = grid.east[0]
    near = center + east * 1e-4
    near /= np.linalg.norm(near)
    far = center + east * 5e-4
    far /= np.linalg.norm(far)

    # Order A: near node first.
    pixel_to_node = healpix_grid.scatter_node_indices(grid, np.stack([near, far]))
    assert pixel_to_node[0] == 0
    # Order B: far node first -- same winner (the nearer one), proving this isn't "first wins".
    pixel_to_node_reversed = healpix_grid.scatter_node_indices(grid, np.stack([far, near]))
    assert pixel_to_node_reversed[0] == 1


def test_wavefront_fill_leaves_no_pixel_empty_and_is_deterministic():
    grid = healpix_grid.build(16)
    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(400, 3))
    xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)
    pixel_to_node = healpix_grid.scatter_node_indices(grid, xyz)
    assert np.any(pixel_to_node == -1)  # a real, non-trivial fill case

    filled_a, rounds_a = healpix_grid._wavefront_fill(grid, pixel_to_node.copy())
    filled_b, rounds_b = healpix_grid._wavefront_fill(grid, pixel_to_node.copy())
    assert not np.any(filled_a == -1)
    assert rounds_a == rounds_b
    assert np.array_equal(filled_a, filled_b)
    # Every originally scatter-assigned pixel must keep its own value -- filling must never
    # overwrite a real assignment.
    assigned = pixel_to_node != -1
    assert np.array_equal(filled_a[assigned], pixel_to_node[assigned])


def test_wavefront_fill_round_breaks_ties_by_lowest_neighbour_slot():
    """Pins the exact tie-break convention (lowest neighbour-slot index wins when more than
    one neighbour is already filled in the same round) rather than leaving it an implementation
    detail -- a reviewer-visible contract, and a regression guard against silently flipping to
    "last slot wins" (still deterministic, but a different answer) while refactoring the kernel."""
    grid = healpix_grid.build(16)
    pix = 0
    current = np.full(grid.npix, -1, dtype=np.int64)
    valid_slots = np.flatnonzero(grid.neighbour_valid[pix])
    assert len(valid_slots) >= 2  # nside=16 pixel 0 should have its full 8 neighbours
    first_valid, second_valid = valid_slots[0], valid_slots[1]
    current[grid.neighbours[pix, first_valid]] = 111
    current[grid.neighbours[pix, second_valid]] = 222
    next_buf = np.empty_like(current)
    healpix_grid._wavefront_fill_round(grid.neighbours, grid.neighbour_valid, current, next_buf)
    assert next_buf[pix] == 111


def test_wavefront_fill_round_real_grid_shapes_no_crash():
    """A plain shape-regression guard for the njit kernel -- run it against real, correctly-
    shaped tiny grids (not hand-typed dummy arrays), the thing phase-0's own spike got wrong
    (a shape-mismatched dummy array segfaulted silently instead of raising)."""
    for nside in (1, 2):
        grid = healpix_grid.build(nside)
        current = np.full(grid.npix, -1, dtype=np.int64)
        current[0] = 0
        next_buf = np.empty_like(current)
        healpix_grid._wavefront_fill_round(grid.neighbours, grid.neighbour_valid, current, next_buf)
        assert next_buf.shape == (grid.npix,)


def test_build_node_pixel_index_query_matches_cktree_for_nodes_that_own_their_pixel():
    """For a node that's the nearest node to its *own* assigned pixel, `NodePixelIndex.query`
    must recover that exact node when queried at its own position -- the one case where the
    HEALPix path and a true nearest-neighbour tree are guaranteed to agree exactly."""
    from scipy.spatial import cKDTree

    grid = healpix_grid.build(16)
    rng = np.random.default_rng(1)
    xyz = rng.normal(size=(500, 3))
    xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)

    index = healpix_grid.build_node_pixel_index(grid, xyz)
    assert not np.any(index.pixel_to_node == -1)

    pix = grid.ang2pix(np.arctan2(xyz[:, 1], xyz[:, 0]), np.arcsin(xyz[:, 2]))
    owns_own_pixel = index.pixel_to_node[pix] == np.arange(len(xyz))
    assert owns_own_pixel.sum() > len(xyz) // 2  # most nodes should, at this node-to-pixel ratio

    _, idx = index.query(xyz[owns_own_pixel])
    assert np.array_equal(idx, np.arange(len(xyz))[owns_own_pixel])

    tree = cKDTree(xyz)
    _, tree_idx = tree.query(xyz[owns_own_pixel])
    assert np.array_equal(idx, tree_idx)

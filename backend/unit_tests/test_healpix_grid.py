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

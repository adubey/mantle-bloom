"""HEALPix background grid (spec section 4.1/4.2): equal-area, isolatitude pixels, replacing
the equirectangular lat/lon grid `fluid_dynamics.py`/`atmosphere_cfd.py`/`ocean_cfd.py` run
on -- see fluid_dynamics_healpix.py for the stencils built on top of this module's neighbour
table, and this module's own docstring note on `resample_to_equirect` for how v2's CFD state
still hands wind/current vectors back to v1's `climate.py` (kept equirectangular; out of
scope to also port -- see the plan).

Uses `astropy_healpix` (pip-installable, pure Cython/numpy, no `healpy`/cfitsio system
dependency) rather than hand-rolling nested-scheme pixel indexing, which is notoriously easy
to get subtly wrong at base-pixel boundaries.
"""

from __future__ import annotations

from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy_healpix import HEALPix
from numba import njit, prange

from ..elevation_lines import PLANET_RADIUS_KM

PLANET_RADIUS_M = PLANET_RADIUS_KM * 1000.0

# Same JIT convention fluid_dynamics.py's own equirectangular stencils use (see that module's
# _NUMBA_JIT_KWARGS docstring) -- cache=True persists compiled machine code to disk so only
# the first call per process pays compilation cost.
_NUMBA_JIT_KWARGS = {"cache": True, "parallel": True, "fastmath": True}


@njit(cache=True, fastmath=True)
def ang2pix_nest_scalar(nside: int, order: int, lon_rad: float, lat_rad: float) -> int:
    """Standard z/phi HEALPix nested-scheme pixel lookup (Gorski et al. 2005's own
    `ang2pix_nest`) for one point -- a hand-rolled reimplementation of what
    `astropy_healpix`'s compiled `lonlat_to_healpix` already does, taken on despite this
    module's own prior docstring caution about nested-scheme indexing being "notoriously easy
    to get subtly wrong at base-pixel boundaries": validated bit-for-bit against
    `astropy_healpix` across every pixel center at every `NSIDE_CHOICES` value (an exact
    round-trip -- every pixel's own center must map back to itself) plus 200,000 uniformly-
    sampled random sphere points per nside, zero mismatches (see
    unit_tests/v2/test_healpix_grid.py's own regression test). Exists purely for speed: this
    was `semi_lagrangian_advect`'s single costliest call (profiling showed ~700us/call going
    through astropy's own Quantity-wrapped, single-threaded path). A scalar function (not
    array-in/array-out) so `fluid_dynamics_healpix._semi_lagrangian_advect_kernel` can inline
    it directly inside its own per-pixel loop -- fusing offset/normalize/lon-lat/pixel-lookup
    into one pass instead of allocating separate (npix,)/(npix, 3) intermediates and making a
    second jitted call for this step alone (`_ang2pix_nested_kernel` below still wraps this
    as a plain array lookup for `HealpixGrid.ang2pix`'s other caller, `resample_to_equirect`)."""
    two_pi = 2.0 * np.pi
    half_pi = np.pi / 2.0
    z = np.sin(lat_rad)
    phi = lon_rad % two_pi
    za = abs(z)
    tt = phi / half_pi  # in [0, 4)

    if za <= 2.0 / 3.0:  # equatorial belt
        temp1 = nside * (0.5 + tt)
        temp2 = nside * z * 0.75
        jp = int(np.floor(temp1 - temp2))  # ascending edge line index
        jm = int(np.floor(temp1 + temp2))  # descending edge line index
        ifp = jp // nside
        ifm = jm // nside
        if ifp == ifm:
            face_num = (ifp % 4) + 4
        elif ifp < ifm:
            face_num = ifp % 4
        else:
            face_num = (ifm % 4) + 8
        ix = jm % nside
        iy = nside - (jp % nside) - 1
    else:  # polar cap
        ntt = int(np.floor(tt))
        if ntt >= 4:
            ntt = 3
        tp = tt - ntt
        tmp = nside * np.sqrt(3.0 * (1.0 - za))
        jp = int(np.floor(tp * tmp))
        jm = int(np.floor((1.0 - tp) * tmp))
        if jp > nside - 1:
            jp = nside - 1
        if jm > nside - 1:
            jm = nside - 1
        if z >= 0.0:
            face_num = ntt
            ix = nside - jm - 1
            iy = nside - jp - 1
        else:
            face_num = ntt + 8
            ix = jp
            iy = jm

    ipf = 0  # bit-interleave ix (even bits) and iy (odd bits)
    for b in range(order):
        ipf |= ((ix >> b) & 1) << (2 * b)
        ipf |= ((iy >> b) & 1) << (2 * b + 1)
    return ipf + face_num * nside * nside


@njit(**_NUMBA_JIT_KWARGS)
def _ang2pix_nested_kernel(nside: int, order: int, lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    n = lon_rad.shape[0]
    out = np.empty(n, dtype=np.int64)
    for idx in prange(n):
        out[idx] = ang2pix_nest_scalar(nside, order, lon_rad[idx], lat_rad[idx])
    return out

# UI-facing choices for World.fluid_density_v2 -- npix ~= 12*nside^2, so this spans roughly
# 3,072 / 12,288 / 49,152 pixels, a comparable order of magnitude to v1's own
# climate.FLUID_DENSITY_CHOICES equirectangular grid sizes at the same density labels.
NSIDE_CHOICES = {0.5: 16, 1.0: 32, 2.0: 64}
DEFAULT_FLUID_DENSITY = 1.0


def nside_for_density(density: float) -> int:
    return NSIDE_CHOICES.get(density, NSIDE_CHOICES[DEFAULT_FLUID_DENSITY])


@dataclass
class HealpixGrid:
    nside: int
    order: int  # log2(nside) -- NSIDE_CHOICES are all powers of 2; used by ang2pix's bit-interleave
    npix: int
    lon_rad: np.ndarray  # (npix,)
    lat_rad: np.ndarray  # (npix,)
    world_xyz: np.ndarray  # (npix, 3)
    neighbours: np.ndarray  # (npix, 8) int, -1 where a pixel has fewer than 8 (rare, low-nside corners)
    neighbour_valid: np.ndarray  # (npix, 8) bool
    neighbour_distance_m: np.ndarray  # (npix, 8) great-circle distance to each neighbour, inf where invalid
    pixel_area_m2: float  # equal for every pixel, by construction
    east: np.ndarray  # (npix, 3) local tangent-plane "east" unit vector per pixel
    north: np.ndarray  # (npix, 3) local tangent-plane "north" unit vector per pixel
    neighbour_dx_m: np.ndarray  # (npix, 8) each valid neighbour's offset along `east`, meters
    neighbour_dy_m: np.ndarray  # (npix, 8) each valid neighbour's offset along `north`, meters
    gradient_inv: np.ndarray  # (npix, 2, 2) precomputed (A^T A)^-1 for the least-squares gradient fit -- see fluid_dynamics_healpix.gradient
    _healpix: HEALPix

    def ang2pix(self, lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
        lon_rad = np.ascontiguousarray(lon_rad, dtype=np.float64)
        lat_rad = np.ascontiguousarray(lat_rad, dtype=np.float64)
        return _ang2pix_nested_kernel(self.nside, self.order, lon_rad, lat_rad)


def _lonlat_to_xyz(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    coslat = np.cos(lat_rad)
    return np.stack([coslat * np.cos(lon_rad), coslat * np.sin(lon_rad), np.sin(lat_rad)], axis=-1)


def _great_circle_distance_m(a_lon: np.ndarray, a_lat: np.ndarray, b_lon: np.ndarray, b_lat: np.ndarray) -> np.ndarray:
    a_xyz = _lonlat_to_xyz(a_lon, a_lat)
    b_xyz = _lonlat_to_xyz(b_lon, b_lat)
    dot = np.clip(np.sum(a_xyz * b_xyz, axis=-1), -1.0, 1.0)
    return np.arccos(dot) * PLANET_RADIUS_M


def build(nside: int) -> HealpixGrid:
    hp = HEALPix(nside=nside, order="nested")
    npix = hp.npix
    pix = np.arange(npix)
    lon, lat = hp.healpix_to_lonlat(pix)
    lon_rad = lon.to_value(u.rad)
    lat_rad = lat.to_value(u.rad)
    world_xyz = _lonlat_to_xyz(lon_rad, lat_rad)

    # astropy_healpix.neighbours returns shape (8, npix), -1 for a missing slot (rare, only
    # ever near the 8 corners shared by 3 base pixels at very low nside) -- transposed here to
    # (npix, 8) to match this module's own per-pixel-row convention.
    raw_neighbours = np.asarray(hp.neighbours(pix)).T
    valid = raw_neighbours >= 0
    safe_neighbours = np.where(valid, raw_neighbours, 0)  # placeholder index, masked out via `valid` everywhere it's read
    neighbour_lon = lon_rad[safe_neighbours]
    neighbour_lat = lat_rad[safe_neighbours]
    distance = _great_circle_distance_m(lon_rad[:, None], lat_rad[:, None], neighbour_lon, neighbour_lat)
    distance = np.where(valid, distance, np.inf)

    pixel_area_sr = float(hp.pixel_area.to_value(u.steradian))
    pixel_area_m2 = pixel_area_sr * PLANET_RADIUS_M**2

    # Local tangent-plane basis per pixel -- "east" tangent to the local latitude circle,
    # "north" completing a right-handed (east, north, outward) frame. Falls back to a
    # different reference pole exactly at the poles themselves (where "east" from the
    # z-axis is undefined) -- HEALPix's own pixel centers never land exactly on a pole, but
    # this keeps the formula well-conditioned arbitrarily close to one regardless.
    z_axis = np.array([0.0, 0.0, 1.0])
    near_pole = np.abs(world_xyz[:, 2]) > 0.999
    pole_ref = np.where(near_pole[:, None], np.array([1.0, 0.0, 0.0]), z_axis)
    east = np.cross(pole_ref, world_xyz)
    east = east / np.linalg.norm(east, axis=-1, keepdims=True)
    north = np.cross(world_xyz, east)

    # Each valid neighbour's offset (meters) along this pixel's own (east, north) -- a linear
    # (chord-projection) approximation of the true geodesic offset, accurate to a small
    # fraction of a percent at real HEALPix neighbour separations (a few hundred km against a
    # 6371km radius -- see fluid_dynamics_healpix.py's own module docstring).
    delta = world_xyz[raw_neighbours.clip(min=0)] - world_xyz[:, None, :]
    neighbour_dx_m = np.einsum("nkj,nj->nk", delta, east) * PLANET_RADIUS_M
    neighbour_dy_m = np.einsum("nkj,nj->nk", delta, north) * PLANET_RADIUS_M
    neighbour_dx_m = np.where(valid, neighbour_dx_m, 0.0)
    neighbour_dy_m = np.where(valid, neighbour_dy_m, 0.0)

    # (A^T A)^-1 for the per-pixel weighted least-squares gradient fit (see
    # fluid_dynamics_healpix.gradient): A's rows are each valid neighbour's own (dx, dy), so
    # A^T A is the 2x2 matrix [[sum dx^2, sum dx*dy], [sum dx*dy, sum dy^2]] -- inverted in
    # closed form (a 2x2 inverse needs no iterative solve) once here, reused by every
    # gradient call for the life of this grid.
    sxx = np.sum(neighbour_dx_m**2, axis=1)
    syy = np.sum(neighbour_dy_m**2, axis=1)
    sxy = np.sum(neighbour_dx_m * neighbour_dy_m, axis=1)
    det = sxx * syy - sxy * sxy
    det = np.where(np.abs(det) < 1e-12, 1e-12, det)
    gradient_inv = np.empty((npix, 2, 2))
    gradient_inv[:, 0, 0] = syy / det
    gradient_inv[:, 0, 1] = -sxy / det
    gradient_inv[:, 1, 0] = -sxy / det
    gradient_inv[:, 1, 1] = sxx / det

    return HealpixGrid(
        nside=nside,
        order=int(round(np.log2(nside))),
        npix=npix,
        lon_rad=lon_rad,
        lat_rad=lat_rad,
        world_xyz=world_xyz,
        neighbours=raw_neighbours,
        neighbour_valid=valid,
        neighbour_distance_m=distance,
        pixel_area_m2=pixel_area_m2,
        east=east,
        north=north,
        # float32 (not float64, despite the float64 math above) -- these feed directly into
        # fluid_dynamics_healpix.gradient's per-substep hot loop against u/v/eta/... state
        # that ocean_cfd_v2.py/atmosphere_cfd_v2.py already keep float32 throughout (mirroring
        # v1's own fluid_dynamics.py discipline, see that module's dtype comment); leaving
        # these float64 silently promoted every gradient/laplacian/divergence call's output
        # back to float64 regardless of the input field's own dtype, doubling memory traffic
        # through the exact loop profiling identified as this solver's single biggest cost.
        neighbour_dx_m=neighbour_dx_m.astype(np.float32),
        neighbour_dy_m=neighbour_dy_m.astype(np.float32),
        gradient_inv=gradient_inv.astype(np.float32),
        _healpix=hp,
    )


def resample_from_equirect(grid: HealpixGrid, field_hw: np.ndarray, lat_deg_rows: np.ndarray) -> np.ndarray:
    """Nearest-cell lookup from an existing `(H, W)` equirectangular field (row 0 = north
    pole, `lat_deg_rows` its own row latitudes, columns evenly spaced in longitude -- the
    same convention `climate.py`/`fluid_dynamics.py` use) onto this HEALPix grid's `(npix,)`
    pixels -- used once, at `init_*_cfd`, to bootstrap v2's fluid state from v1's diagnostic
    climate snapshot the same way v1's own `init_atmosphere_cfd` bootstraps from
    `climate.compute_wind`."""
    height, width = field_hw.shape
    row = np.clip(np.searchsorted(-lat_deg_rows, -np.degrees(grid.lat_rad)), 0, height - 1)
    # Column 0 is longitude -180 degrees (this codebase's own convention -- see
    # climate.py's `_build_grid`/lat_long_grid.py), not 0 degrees -- normalize into
    # [-180, 180) first before mapping to a column index.
    lon_deg = ((np.degrees(grid.lon_rad) + 180.0) % 360.0) - 180.0
    col = np.clip(((lon_deg + 180.0) / 360.0 * width).astype(int), 0, width - 1)
    return field_hw[row, col]


def resample_to_equirect(grid: HealpixGrid, field_pix: np.ndarray, height: int, width: int) -> np.ndarray:
    """The reverse of `resample_from_equirect` -- HEALPix `(npix,)` values onto a fresh
    `(H, W)` equirectangular grid, nearest-pixel lookup via `ang2pix`. This is the seam
    `climate.py`'s `resample_uv_to_equirect` (added to both v1 and v2 CFD state classes)
    calls so `compute_climate`'s render/erosion grid stays equirectangular and oblivious to
    which grid actually produced the wind/current field it's reading."""
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_deg = (np.arange(width) + 0.5) * (360.0 / width) - 180.0
    lat_grid, lon_grid = np.meshgrid(np.radians(lat_deg), np.radians(lon_deg), indexing="ij")
    pix = grid.ang2pix(lon_grid.ravel(), lat_grid.ravel())
    return field_pix[pix].reshape(height, width)

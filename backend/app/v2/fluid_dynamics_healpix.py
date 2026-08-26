"""HEALPix-native twins of `fluid_dynamics.py`'s six stencil kernels (spec section 4.2) --
same signatures/semantics wherever possible, operating on flat `(npix,)` arrays over a
`healpix_grid.HealpixGrid`'s precomputed neighbour table instead of `(H, W)` row/col
indexing. `cfl_substeps`/`coriolis_parameter` are reused directly from `fluid_dynamics.py`
(confirmed grid-agnostic, see that module). The polar zonal filter/sponge-drag pair is simply
not ported here -- there is no polar singularity on an equal-area grid to compensate for.

**Gradient/Laplacian/divergence** use a per-pixel weighted least-squares fit against each
pixel's own (up to 8) neighbours' tangent-plane offsets (precomputed once at grid build time,
`HealpixGrid.gradient_inv`/`neighbour_dx_m`/`neighbour_dy_m`) rather than a fixed 5-point
stencil -- HEALPix neighbours aren't perfectly axis-aligned or equidistant the way a regular
lat/lon grid's row/col neighbours are, so a least-squares plane fit (the standard meshless-
gradient technique) is the appropriate generalization, not an approximation of convenience.

**Now `@njit`-compiled, not plain vectorized NumPy.** This module's own prior docstring
reasoned that HEALPix grids here (a few thousand to ~50k pixels, `NSIDE_CHOICES`) stay
"comfortably within NumPy's own vectorized performance envelope" -- true of any *one* call,
but profiling `step_world_v2` at production density (12 plates, node/climate/fluid density
1.0) found `gradient` alone consuming ~66% of total step time, because it's called thousands
of times per tectonic step (several times per fluid substep, `MAX_SUBSTEPS_PER_STEP` up to
2000 substeps per solver per step). Exactly the call-count-not-array-size mismatch
`fluid_dynamics.py`'s own equirectangular stencils hit before *their* JIT pass (see that
module's `_NUMBA_JIT_KWARGS` docstring) -- each plain-NumPy call here paid the same
allocate-a-temporary-per-op cost (`neighbour_field`, `b`, two `np.sum(..., axis=1)` reductions)
every single substep. The `@njit(parallel=True)` kernels below fuse each stencil into one pass
over pixels with no intermediate allocation, `prange` splitting the pixel loop across cores --
same technique, same `_NUMBA_JIT_KWARGS`, ported one level from a fixed 4-neighbour row/col
stencil to a variable-length (up to 8) neighbour-table gather.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from .healpix_grid import PLANET_RADIUS_M, HealpixGrid

# Same convention as fluid_dynamics.py's own _NUMBA_JIT_KWARGS (see that module's docstring
# for cache=True/parallel=True/fastmath=True's own rationale) -- kept as a separate constant
# here (not imported from there) since these kernels' signatures are HEALPix-specific.
_NUMBA_JIT_KWARGS = {"cache": True, "parallel": True, "fastmath": True}

# fluid_dynamics.cfl_substeps' own 0.4 safety factor was tuned against v1's fixed 5-point
# row/col stencil; this module's least-squares gradient (and the Laplacian composed from it)
# has a larger effective amplification factor on HEALPix's irregular neighbour spacing --
# confirmed directly: the shallow-water gravity-wave mode was still exponentially unstable at
# a CFL-computed dt with an extra safety factor of 2 (u/eta both diverging within ~50
# substeps), stable at an extra factor of 4, used here with headroom. Every `cfl_substeps`
# call in this module's callers (atmosphere_cfd_v2.py/ocean_cfd_v2.py) divides its own
# min-spacing input by this before passing it in, rather than each re-deriving/re-confirming
# the same factor independently.
CFL_STENCIL_SAFETY_DIVISOR = 6.0


def min_spacing_m(grid: HealpixGrid) -> float:
    return float(grid.neighbour_distance_m[grid.neighbour_valid].min()) / CFL_STENCIL_SAFETY_DIVISOR


@njit(**_NUMBA_JIT_KWARGS)
def _gradient_kernel(
    field: np.ndarray,
    neighbours: np.ndarray,
    neighbour_valid: np.ndarray,
    neighbour_dx_m: np.ndarray,
    neighbour_dy_m: np.ndarray,
    gradient_inv: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    npix = field.shape[0]
    n_neighbours = neighbours.shape[1]
    gx = np.empty(npix, dtype=np.float64)
    gy = np.empty(npix, dtype=np.float64)
    for i in prange(npix):
        fi = field[i]
        atb_x = 0.0
        atb_y = 0.0
        for k in range(n_neighbours):
            if neighbour_valid[i, k]:
                d = field[neighbours[i, k]] - fi
                atb_x += neighbour_dx_m[i, k] * d
                atb_y += neighbour_dy_m[i, k] * d
        gx[i] = gradient_inv[i, 0, 0] * atb_x + gradient_inv[i, 0, 1] * atb_y
        gy[i] = gradient_inv[i, 1, 0] * atb_x + gradient_inv[i, 1, 1] * atb_y
    return gx, gy


def gradient(field: np.ndarray, grid: HealpixGrid) -> tuple[np.ndarray, np.ndarray]:
    """(d/d-east, d/d-north) in real per-meter units -- the HEALPix analogue of
    `fluid_dynamics.gradient_m`. Solves, per pixel, the weighted least-squares system
    `A [gx, gy]^T = b` where `A`'s rows are each valid neighbour's own (dx, dy) offset and
    `b`'s entries are that neighbour's `field` value minus this pixel's own -- i.e. the best-
    fit local linear ramp explaining the observed differences to every neighbour at once,
    rather than a hand-picked pair of opposite neighbours the way a regular grid's 5-point
    stencil picks. See `_NUMBA_JIT_KWARGS`'s own comment for why the actual math lives in a
    jitted `_gradient_kernel` instead of here. `gx`/`gy` come back float64 regardless of
    `field`'s own dtype (Numba follows the same NumPy promotion rule the prior plain-NumPy
    form did: `gradient_inv`/`neighbour_dx_m`/`neighbour_dy_m` are float32, but the `atb_x`/
    `atb_y` accumulators are genuine Python floats -- float64 -- so every product involving
    them promotes back up) -- unchanged from this function's pre-JIT behavior, not a
    regression introduced by this pass."""
    return _gradient_kernel(field, grid.neighbours, grid.neighbour_valid, grid.neighbour_dx_m, grid.neighbour_dy_m, grid.gradient_inv)


def laplacian(field: np.ndarray, grid: HealpixGrid) -> np.ndarray:
    """div(grad(field)) -- applying the already-validated least-squares `gradient` operator
    twice, rather than a hand-derived single-pass mesh-Laplacian formula. An earlier version
    used the standard graph/mesh ("umbrella") formula `(2/n) * sum_k (field_k - field_center)
    / dist_k^2` directly; confirmed directly that its scaling doesn't match the stability
    assumptions `VISCOSITY_M2_S`/`cfl_substeps` were tuned against (ported unchanged from
    v1's own 5-point-stencil-calibrated equirectangular values) -- ocean_cfd_v2's own
    diffusion term blew up to >1e6 m/s within one step's substep loop. Composing the gradient
    operator with itself keeps the two dimensionally and numerically consistent by
    construction, at the cost of a second least-squares solve per call -- three total jitted
    `_gradient_kernel` calls, each now a single fused pass rather than several allocating
    NumPy ops, so this stays cheap even called `VISCOSITY_M2_S`-many times per substep."""
    gx, gy = gradient(field, grid)
    d2x, _ = gradient(gx, grid)
    _, d2y = gradient(gy, grid)
    return d2x + d2y


def divergence(u: np.ndarray, v: np.ndarray, grid: HealpixGrid) -> np.ndarray:
    """du/d-east + dv/d-north -- reuses `gradient`'s own least-squares fit for each
    component, same convention `fluid_dynamics.divergence_m` follows."""
    du_dx, _ = gradient(u, grid)
    _, dv_dy = gradient(v, grid)
    return du_dx + dv_dy


@njit(**_NUMBA_JIT_KWARGS)
def _grid_noise_filter_kernel(field: np.ndarray, neighbours: np.ndarray, neighbour_valid: np.ndarray, weight: float) -> np.ndarray:
    npix = field.shape[0]
    n_neighbours = neighbours.shape[1]
    out = np.empty_like(field)
    for i in prange(npix):
        total = 0.0
        count = 0
        for k in range(n_neighbours):
            if neighbour_valid[i, k]:
                total += field[neighbours[i, k]]
                count += 1
        if count == 0:
            count = 1
        neighbour_avg = total / count
        out[i] = field[i] + weight * (neighbour_avg - field[i])
    return out


def grid_noise_filter(field: np.ndarray, grid: HealpixGrid, weight: float = 0.05) -> np.ndarray:
    """`fluid_dynamics.grid_noise_filter`'s own light neighbour-average blend, over this
    pixel's real HEALPix neighbours instead of a 4-neighbour row/col stencil. See
    `_NUMBA_JIT_KWARGS`'s own comment for why the actual math lives in a jitted
    `_grid_noise_filter_kernel` instead of here."""
    return _grid_noise_filter_kernel(field, grid.neighbours, grid.neighbour_valid, weight)


def coastal_ocean_mask(is_ocean: np.ndarray, grid: HealpixGrid) -> np.ndarray:
    is_land = ~is_ocean
    neighbour_is_land = is_land[grid.neighbours.clip(min=0)] & grid.neighbour_valid
    return is_ocean & neighbour_is_land.any(axis=1)


def semi_lagrangian_advect(field: np.ndarray, u: np.ndarray, v: np.ndarray, dt_s: float, grid: HealpixGrid) -> np.ndarray:
    """Backward-trace each pixel along (u, v) by `dt_s` (real velocity in m/s along this
    pixel's own east/north tangent directions) and sample `field` there, nearest-pixel --
    same unconditionally-stable technique `fluid_dynamics.semi_lagrangian_advect` uses, with
    the row/col closed-form index lookup replaced by a real `ang2pix` query."""
    offset_m = (-u * dt_s)[:, None] * grid.east + (-v * dt_s)[:, None] * grid.north
    src_xyz = grid.world_xyz + offset_m / PLANET_RADIUS_M
    src_xyz = src_xyz / np.linalg.norm(src_xyz, axis=-1, keepdims=True)
    src_lon = np.arctan2(src_xyz[:, 1], src_xyz[:, 0])
    src_lat = np.arcsin(np.clip(src_xyz[:, 2], -1.0, 1.0))
    src_pix = grid.ang2pix(src_lon, src_lat)
    return field[src_pix]

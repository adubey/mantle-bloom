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
Implemented in plain vectorized NumPy rather than `@njit` kernels -- HEALPix grids here run at
a few thousand to ~50k pixels (`NSIDE_CHOICES`), comfortably within NumPy's own vectorized
performance envelope without needing per-cell JIT compilation the way `fluid_dynamics.py`'s
equirectangular kernels (running at similar or larger cell counts, but a much simpler,
JIT-friendly fixed 4-neighbour stencil) benefit from.
"""

from __future__ import annotations

import numpy as np

from .healpix_grid import PLANET_RADIUS_M, HealpixGrid

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


def gradient(field: np.ndarray, grid: HealpixGrid) -> tuple[np.ndarray, np.ndarray]:
    """(d/d-east, d/d-north) in real per-meter units -- the HEALPix analogue of
    `fluid_dynamics.gradient_m`. Solves, per pixel, the weighted least-squares system
    `A [gx, gy]^T = b` where `A`'s rows are each valid neighbour's own (dx, dy) offset and
    `b`'s entries are that neighbour's `field` value minus this pixel's own -- i.e. the best-
    fit local linear ramp explaining the observed differences to every neighbour at once,
    rather than a hand-picked pair of opposite neighbours the way a regular grid's 5-point
    stencil picks."""
    neighbour_field = field[grid.neighbours.clip(min=0)]
    b = np.where(grid.neighbour_valid, neighbour_field - field[:, None], 0.0)
    atb_x = np.sum(grid.neighbour_dx_m * b, axis=1)
    atb_y = np.sum(grid.neighbour_dy_m * b, axis=1)
    gx = grid.gradient_inv[:, 0, 0] * atb_x + grid.gradient_inv[:, 0, 1] * atb_y
    gy = grid.gradient_inv[:, 1, 0] * atb_x + grid.gradient_inv[:, 1, 1] * atb_y
    return gx, gy


def laplacian(field: np.ndarray, grid: HealpixGrid) -> np.ndarray:
    """div(grad(field)) -- applying the already-validated least-squares `gradient` operator
    twice, rather than a hand-derived single-pass mesh-Laplacian formula. An earlier version
    used the standard graph/mesh ("umbrella") formula `(2/n) * sum_k (field_k - field_center)
    / dist_k^2` directly; confirmed directly that its scaling doesn't match the stability
    assumptions `VISCOSITY_M2_S`/`cfl_substeps` were tuned against (ported unchanged from
    v1's own 5-point-stencil-calibrated equirectangular values) -- ocean_cfd_v2's own
    diffusion term blew up to >1e6 m/s within one step's substep loop. Composing the gradient
    operator with itself keeps the two dimensionally and numerically consistent by
    construction, at the cost of a second least-squares solve per call."""
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


def grid_noise_filter(field: np.ndarray, grid: HealpixGrid, weight: float = 0.05) -> np.ndarray:
    """`fluid_dynamics.grid_noise_filter`'s own light neighbour-average blend, over this
    pixel's real HEALPix neighbours instead of a 4-neighbour row/col stencil."""
    neighbour_field = field[grid.neighbours.clip(min=0)]
    n_valid = np.maximum(grid.neighbour_valid.sum(axis=1), 1)
    neighbour_avg = np.sum(np.where(grid.neighbour_valid, neighbour_field, 0.0), axis=1) / n_valid
    return field + weight * (neighbour_avg - field)


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

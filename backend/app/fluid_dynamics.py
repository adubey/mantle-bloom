"""Shared numerical primitives for ocean_cfd.py/atmosphere_cfd.py: real physical-unit
gradients/Laplacians/divergence on climate.py's own fixed equirectangular grid, real Coriolis
parameter, CFL-stable substep sizing, and semi-Lagrangian advection -- the low-level array
math both shallow-water solvers share, factored out once rather than duplicated (see
docs/simulation-model.md#ocean-atmospheric-fluid-dynamics for why both solvers exist and how
they use this).

**Why not climate.py's own `_centered_gradient`/`_smooth_field`.** Those compute *raw*
per-cell differences (fine for climate.py's own stylized, unitless heuristics -- a
"which-way-is-warmer" direction, not a real quantity in m/s^2), not gradients in real
physical units (1/m). A genuine momentum equation needs `d(eta)/dx` in real units, so every
helper here is explicit about, and divides through by, real cell spacing in meters -- see
`grid_spacing_m`.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from . import plates

# Every stencil below (gradient/laplacian/divergence/grid_noise_filter/polar_zonal_filter/
# semi_lagrangian_advect) is called several times per substep, up to MAX_SUBSTEPS_PER_STEP
# times per step_*_cfd call -- profiling a 1-week step at default resolution showed these
# dominating both solvers' wall time even after float32 (see ocean_cfd.init_ocean_cfd's own
# comment) and hoisting invariant work out of the substep loop. Their plain-NumPy form (each
# a handful of np.roll/np.clip/... calls) pays per-call Python/dispatch overhead and
# materializes a full-grid intermediate array for every one of those calls; the `_numba`
# functions below instead compile each into a single fused pass with no intermediate
# allocations and no per-call dispatch, `parallel=True`/`prange` splitting the outer
# (latitude) loop across cores. `cache=True` persists the compiled machine code to disk so
# only the first call per process pays compilation cost, not every one of this session's --
# let alone future sessions' -- step_*_cfd calls.
_NUMBA_JIT_KWARGS = {"cache": True, "parallel": True, "fastmath": True}

# Real Earth rotation rate (rad/s) -- unlike climate.py's own `coriolis_parameter` (a bare
# `sin(lat)` proxy, fine for that module's stylized, unitless deflection heuristics), a real
# momentum equation needs a dimensionally real Coriolis parameter `f = 2*OMEGA*sin(lat)`.
EARTH_ANGULAR_VELOCITY_RAD_S = 7.292e-5


def coriolis_parameter(lat_deg: np.ndarray) -> np.ndarray:
    """Real f = 2*Omega*sin(lat), (H,) broadcastable against an (H, W) field via `[:, None]`."""
    return 2.0 * EARTH_ANGULAR_VELOCITY_RAD_S * np.sin(np.radians(lat_deg))


def grid_spacing_m(lat_deg: np.ndarray, height: int, width: int) -> tuple[np.ndarray, float]:
    """Real east-west cell width per row (meters, shrinking by cos(lat) toward the poles --
    (H,) broadcastable via `[:, None]`) and the single real north-south cell height (meters,
    constant across rows -- a meridian step covers the same real distance everywhere). `dx_m`
    is clamped from below at its own value *at* POLAR_FILTER_START_LAT_DEG for every row
    poleward of it, matching `polar_zonal_filter`'s own suppression of finer east-west
    structure there (see that function's docstring): without this, `gradient_m`/`laplacian_m`
    -- especially the viscosity/diffusion terms' `1/dx^2` -- blow up at the pole-most rows
    even after the *state* has been smoothed nearly flat there each substep, since the raw,
    unclamped spacing keeps shrinking all the way to the pole regardless of what the filter
    just removed. Confirmed directly during development: with only a small safety-epsilon
    floor (not this latitude-based clamp), both solvers stayed numerically stable at a
    world's default (coarser) climate resolution but reliably blew up to `inf`/`NaN` within
    the first simulated day at the finest resolution -- exactly the rows this clamp now
    protects."""
    radius_m = plates.PLANET_RADIUS_KM * 1000.0
    dx_m = (2.0 * np.pi * radius_m / width) * np.cos(np.radians(lat_deg))
    # float(...) -- a genuine (dtype-"weak") Python float, not the strong np.float64 scalar
    # np.cos/np.radians would otherwise hand back. lat_deg (and so dx_m) is float32 in both
    # solvers (see ocean_cfd.init_ocean_cfd's own comment); np.maximum against a *strong*
    # float64 scalar here would silently promote dx_m -- and everything downstream of it in
    # gradient_m/laplacian_m/divergence_m -- back to float64 for the whole grid.
    dx_at_filter_start = float((2.0 * np.pi * radius_m / width) * np.cos(np.radians(POLAR_FILTER_START_LAT_DEG)))
    dx_m = np.maximum(dx_m, dx_at_filter_start)
    dy_m = (np.pi * radius_m) / height
    return dx_m, dy_m


# Latitude beyond which polar_zonal_filter progressively suppresses east-west structure, and
# grid_spacing_m clamps dx_m to avoid it shrinking further -- see both functions' own
# docstrings. The real, standard "pole problem" fix lat/lon-grid ocean/atmosphere models use
# (rather than a staggered/reduced/icosahedral grid, a much larger undertaking here).
POLAR_FILTER_START_LAT_DEG = 75.0


def stable_min_spacing_m(dx_m: np.ndarray, dy_m: float) -> float:
    """The real length scale CFL substep sizing should respect -- now that grid_spacing_m
    itself already clamps dx_m poleward of POLAR_FILTER_START_LAT_DEG, this is just the
    grid's own minimum spacing, kept as a named helper (rather than inlining `min(dx_m.min(),
    dy_m)` at each call site) so both solvers stay explicit about *why* they're taking this
    particular minimum for CFL purposes."""
    return float(min(dx_m.min(), dy_m))


def polar_sponge_drag_per_s(lat_deg: np.ndarray, max_extra_drag_per_s: float) -> np.ndarray:
    """An extra momentum-damping *rate* (per second, (H,) broadcastable via `[:, None]`),
    zero equatorward of POLAR_FILTER_START_LAT_DEG and ramping linearly to
    `max_extra_drag_per_s` at the pole-adjacent row -- a standard "sponge layer" real polar-
    cap treatments add alongside `polar_zonal_filter`: suppressing zonal structure there also
    suppresses the small-scale eddies that would normally dissipate momentum, so without
    extra damping specifically in that zone, ordinary forcing (wind stress, pressure
    gradient) can spin the whole zonally-averaged polar band up into a persistent, growing
    drift rather than the bounded circulation the rest of the grid settles into -- confirmed
    directly during development: without this, ocean_cfd.py's polar band kept accelerating
    past 10x a typical mid-latitude current speed over a 10-day run with no sign of
    leveling off, even though nothing elsewhere on the grid was unstable."""
    abs_lat = np.abs(lat_deg)
    span = 90.0 - POLAR_FILTER_START_LAT_DEG
    strength = np.clip((abs_lat - POLAR_FILTER_START_LAT_DEG) / span, 0.0, 1.0)
    return (max_extra_drag_per_s * strength)[:, None]


@njit(**_NUMBA_JIT_KWARGS)
def _polar_zonal_filter_kernel(field: np.ndarray, strength: np.ndarray) -> np.ndarray:
    height, width = field.shape
    out = np.empty_like(field)
    for i in prange(height):
        s = strength[i]
        if s == 0.0:
            for j in range(width):
                out[i, j] = field[i, j]
            continue
        total = 0.0
        for j in range(width):
            total += field[i, j]
        blended_mean = (total / width) * s
        one_minus_s = 1.0 - s
        for j in range(width):
            out[i, j] = field[i, j] * one_minus_s + blended_mean
    return out


def polar_zonal_filter(field: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
    """Blends each row within POLAR_FILTER_START_LAT_DEG of either pole toward its own zonal
    (east-west) mean, ramping linearly to a full zonal average exactly at the pole-adjacent
    row -- see stable_min_spacing_m's own docstring for why this is what lets CFL substep
    sizing ignore the poles' own shrinking spacing rather than being dictated by it. See
    `_NUMBA_JIT_KWARGS`'s own comment for why the actual math lives in a jitted
    `_polar_zonal_filter_kernel` instead of here -- `strength` itself (an (H,) row vector) is
    cheap enough to keep recomputing in plain NumPy every call."""
    abs_lat = np.abs(lat_deg)
    span = 90.0 - POLAR_FILTER_START_LAT_DEG
    strength = np.clip((abs_lat - POLAR_FILTER_START_LAT_DEG) / span, 0.0, 1.0)
    return _polar_zonal_filter_kernel(field, strength)


@njit(**_NUMBA_JIT_KWARGS)
def _gradient_m_kernel(field: np.ndarray, dx_m: np.ndarray, dy_m: float) -> tuple[np.ndarray, np.ndarray]:
    height, width = field.shape
    gx = np.empty_like(field)
    gy = np.empty_like(field)
    inv_2dy = 1.0 / (2.0 * dy_m)
    for i in prange(height):
        im1 = height - 1 if i == 0 else i - 1
        ip1 = 0 if i == height - 1 else i + 1
        inv_2dx = 1.0 / (2.0 * dx_m[i])
        for j in range(width):
            jm1 = width - 1 if j == 0 else j - 1
            jp1 = 0 if j == width - 1 else j + 1
            gx[i, j] = (field[i, jp1] - field[i, jm1]) * inv_2dx
            gy[i, j] = (field[im1, j] - field[ip1, j]) * inv_2dy
    return gx, gy


def gradient_m(field: np.ndarray, dx_m: np.ndarray, dy_m: float) -> tuple[np.ndarray, np.ndarray]:
    """(d/dx eastward, d/dy northward) in real per-meter units, longitude-wrapping. Row 0 =
    north pole (climate.py's own convention throughout), so "northward" is the *negative* row
    direction, matching climate.py's `_centered_gradient`. See `_NUMBA_JIT_KWARGS`'s own
    comment for why the actual math lives in a jitted `_gradient_m_kernel` instead of here."""
    return _gradient_m_kernel(field, dx_m, dy_m)


@njit(**_NUMBA_JIT_KWARGS)
def _laplacian_m_kernel(field: np.ndarray, dx_m: np.ndarray, dy_m: float) -> np.ndarray:
    height, width = field.shape
    out = np.empty_like(field)
    inv_dy2 = 1.0 / (dy_m * dy_m)
    for i in prange(height):
        im1 = height - 1 if i == 0 else i - 1
        ip1 = 0 if i == height - 1 else i + 1
        inv_dx2 = 1.0 / (dx_m[i] * dx_m[i])
        for j in range(width):
            jm1 = width - 1 if j == 0 else j - 1
            jp1 = 0 if j == width - 1 else j + 1
            d2x = (field[i, jp1] + field[i, jm1] - 2.0 * field[i, j]) * inv_dx2
            d2y = (field[ip1, j] + field[im1, j] - 2.0 * field[i, j]) * inv_dy2
            out[i, j] = d2x + d2y
    return out


def laplacian_m(field: np.ndarray, dx_m: np.ndarray, dy_m: float) -> np.ndarray:
    """Standard 5-point Laplacian in real per-square-meter units, longitude-wrapping. See
    `_NUMBA_JIT_KWARGS`'s own comment for why the actual math lives in a jitted
    `_laplacian_m_kernel` instead of here."""
    return _laplacian_m_kernel(field, dx_m, dy_m)


@njit(**_NUMBA_JIT_KWARGS)
def _divergence_m_kernel(u: np.ndarray, v: np.ndarray, dx_m: np.ndarray, dy_m: float) -> np.ndarray:
    height, width = u.shape
    out = np.empty_like(u)
    inv_2dy = 1.0 / (2.0 * dy_m)
    for i in prange(height):
        im1 = height - 1 if i == 0 else i - 1
        ip1 = 0 if i == height - 1 else i + 1
        inv_2dx = 1.0 / (2.0 * dx_m[i])
        for j in range(width):
            jm1 = width - 1 if j == 0 else j - 1
            jp1 = 0 if j == width - 1 else j + 1
            dudx = (u[i, jp1] - u[i, jm1]) * inv_2dx
            dvdy = (v[im1, j] - v[ip1, j]) * inv_2dy
            out[i, j] = dudx + dvdy
    return out


def divergence_m(u: np.ndarray, v: np.ndarray, dx_m: np.ndarray, dy_m: float) -> np.ndarray:
    """du/dx + dv/dy in real per-second units (given u, v in m/s), same gradient convention
    as gradient_m -- used for the shallow-water continuity equation's flux divergence. Computes
    just the two needed components (du/dx, dv/dy) directly, rather than calling gradient_m(u)/
    gradient_m(v) and discarding half of each. See `_NUMBA_JIT_KWARGS`'s own comment for why
    the actual math lives in a jitted `_divergence_m_kernel` instead of here."""
    return _divergence_m_kernel(u, v, dx_m, dy_m)


def advection_geometry(lat_deg: np.ndarray, width: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Precomputes the part of semi_lagrangian_advect's own backward-trace math that depends
    only on the grid (lat_deg, width) and never on the field/velocity being advected --
    lat_deg itself, lon_deg, cos_lat (both (H,)/(W,) row vectors: latitude-only quantities
    stay one row tall here rather than being broadcast out to the full (H, W) grid, since
    _semi_lagrangian_advect_kernel indexes them by row/column directly instead of relying on
    NumPy broadcasting), and meters_per_deg_lat. Both solvers' `state.lat_deg` is fixed for a
    whole session, so callers compute this once (e.g. per step_*_cfd call) and pass it into
    every semi_lagrangian_advect call that step takes, instead of each of those calls
    (temperature, sediment/humidity, ...) redundantly rebuilding the same arrays from
    scratch."""
    radius_m = plates.PLANET_RADIUS_KM * 1000.0
    # dtype=lat_deg.dtype (not a bare np.arange, which defaults to int64 and would promote
    # the float32 arithmetic below back to float64 -- both solvers pass a float32 lat_deg,
    # see ocean_cfd.init_ocean_cfd's own comment on why that matters here) so lon_deg lands
    # in the same dtype as everything else this feeds into semi_lagrangian_advect.
    lon_deg = -180.0 + (np.arange(width, dtype=lat_deg.dtype) + 0.5) * (360.0 / width)
    cos_lat = np.clip(np.cos(np.radians(lat_deg)), 0.15, 1.0)
    # meters -> degrees: a meridian degree is (pi*R/180) meters; a zonal degree shrinks by
    # cos(lat) the same way grid_spacing_m's own dx_m does.
    meters_per_deg_lat = (np.pi * radius_m) / 180.0
    return lat_deg, lon_deg, cos_lat, meters_per_deg_lat


@njit(**_NUMBA_JIT_KWARGS)
def _semi_lagrangian_advect_kernel(
    field: np.ndarray, u: np.ndarray, v: np.ndarray, dt_s: float, lat_deg: np.ndarray, lon_deg: np.ndarray, cos_lat: np.ndarray, meters_per_deg_lat: float
) -> np.ndarray:
    height, width = field.shape
    out = np.empty_like(field)
    row_per_deg = height / 180.0
    col_per_deg = width / 360.0
    for i in prange(height):
        lat_i = lat_deg[i]
        cos_lat_i = cos_lat[i]
        for j in range(width):
            src_lat = lat_i - (v[i, j] * dt_s) / meters_per_deg_lat
            src_lon = lon_deg[j] - (u[i, j] * dt_s) / (meters_per_deg_lat * cos_lat_i)

            src_row = int(round((90.0 - src_lat) * row_per_deg - 0.5))
            if src_row < 0:
                src_row = 0
            elif src_row > height - 1:
                src_row = height - 1

            src_col = int(round((src_lon + 180.0) * col_per_deg - 0.5)) % width
            out[i, j] = field[src_row, src_col]
    return out


def semi_lagrangian_advect(
    field: np.ndarray, u: np.ndarray, v: np.ndarray, dt_s: float, geometry: tuple[np.ndarray, np.ndarray, np.ndarray, float]
) -> np.ndarray:
    """Backward-trace each cell along (u, v) by dt_s and sample the field there, nearest-cell
    -- the same technique climate.py's own `advect_ocean_temperature`/`_sample_at_offset` use,
    chosen for the same reason: unconditionally stable regardless of how large `dt_s * speed`
    gets relative to one cell, so this never adds a second, stricter CFL constraint on top of
    the shallow-water gravity-wave one substep sizing is already built around. `geometry` is
    this grid's own advection_geometry(lat_deg, width) -- see its docstring for why callers
    precompute it once rather than this function rebuilding it every call. See
    `_NUMBA_JIT_KWARGS`'s own comment for why the actual backward-trace-and-gather math lives
    in a jitted `_semi_lagrangian_advect_kernel` instead of here -- this was the single
    costliest call in either solver's profile even after that hoisting, since (unlike the
    stencils above) it was never just avoidable per-call redundancy: the round/clip/mod/gather
    chain is inherent per-cell work that plain NumPy can't fuse into one pass on its own."""
    lat_deg, lon_deg, cos_lat, meters_per_deg_lat = geometry
    return _semi_lagrangian_advect_kernel(field, u, v, dt_s, lat_deg, lon_deg, cos_lat, meters_per_deg_lat)


def cfl_substeps(seconds: float, min_spacing_m: float, wave_speed_m_s: float, max_advect_speed_m_s: float, max_substeps: int) -> tuple[int, float]:
    """(n_substeps, dt_s) -- as many CFL-stable substeps as `min_spacing_m` (see
    stable_min_spacing_m -- deliberately *not* the raw grid's own shrinking-at-the-poles
    minimum, see its docstring) and gravity-wave speed (plus whatever advection speed is
    already present) demand to cover `seconds` of real time, capped at `max_substeps` (a hard
    safety ceiling on worst-case request time -- mirrors main.py's own MAX_ANIMATION_FRAMES
    precedent) rather than left unbounded. `dt_s` always evenly divides `seconds` (n_substeps
    computed first, dt derived from it) so repeated substepping lands on exactly `seconds`
    elapsed, not slightly over."""
    fastest = max(wave_speed_m_s, max_advect_speed_m_s, 1e-6)
    cfl_safety = 0.4
    dt_stable = cfl_safety * min_spacing_m / fastest
    n_substeps = max(1, int(np.ceil(seconds / dt_stable)))
    n_substeps = min(n_substeps, max_substeps)
    return n_substeps, seconds / n_substeps


@njit(**_NUMBA_JIT_KWARGS)
def _grid_noise_filter_kernel(field: np.ndarray, weight: float) -> np.ndarray:
    height, width = field.shape
    out = np.empty_like(field)
    for i in prange(height):
        im1 = height - 1 if i == 0 else i - 1
        ip1 = 0 if i == height - 1 else i + 1
        for j in range(width):
            jm1 = width - 1 if j == 0 else j - 1
            jp1 = 0 if j == width - 1 else j + 1
            neighbor_avg = 0.25 * (field[i, jp1] + field[i, jm1] + field[ip1, j] + field[im1, j])
            out[i, j] = field[i, j] + weight * (neighbor_avg - field[i, j])
    return out


def grid_noise_filter(field: np.ndarray, weight: float = 0.05) -> np.ndarray:
    """A light Jacobi neighbor-average blend, longitude-wrapping -- both shallow-water
    solvers here use a co-located (unstaggered) finite-difference grid, which is simpler than
    the staggered Arakawa-C grid real ocean/atmosphere models use for the momentum-continuity
    coupling, but supports a spurious grid-scale computational mode the C-grid avoids by
    construction: confirmed directly during development, a very small (~0.1%) per-substep
    amplification that's invisible over a single UI "Step" but compounds into an obviously
    unstable blow-up after enough steps. `weight` is deliberately small -- this needs to
    damp grid-scale (checkerboard) noise, not meaningfully diffuse genuine large-scale
    circulation, so it's applied every substep as a light standing filter rather than a
    one-off cleanup. See `_NUMBA_JIT_KWARGS`'s own comment for why the actual math lives in a
    jitted `_grid_noise_filter_kernel` instead of here."""
    return _grid_noise_filter_kernel(field, weight)


def coastal_ocean_mask(is_ocean: np.ndarray) -> np.ndarray:
    """Ocean cells with at least one land neighbor -- used by ocean_cfd.py's sediment pickup
    (erosion only happens right at the coast, not mid-ocean), same 4-neighbor longitude-
    wrapping adjacency check climate.py's own `_coastal_normal` uses for current deflection."""
    is_land = ~is_ocean
    has_land_neighbor = (
        np.roll(is_land, 1, axis=0) | np.roll(is_land, -1, axis=0) | np.roll(is_land, 1, axis=1) | np.roll(is_land, -1, axis=1)
    )
    return is_ocean & has_land_neighbor

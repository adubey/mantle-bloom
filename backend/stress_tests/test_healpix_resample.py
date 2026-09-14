"""Issue #133 phase 1: the real-scale (not unit-test-scale) `node_cloud_resample_mode` checks.

Phase-0's own spike (docs/profiling.md's "Item 4" section, landed via PR #135) measured its
throwaway scatter+fill script against `cKDTree` at seed=0, node_density=climate_density=4.0,
~130.6K nodes -- the config reproduced here against the real (not throwaway) implementation, so
these numbers are directly comparable to phase-0's own and replace them once Phase 1 lands.
Lives in stress_tests/, not unit_tests/, because building a node_density=4.0 world and stepping
it several times is exactly the "many-step, full-simulation" shape bin/stress_test.sh's own
module docstring describes -- too slow for the fast suite's ~1s-per-test budget.
"""

import time

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion

from app import healpix_grid, render_image
from app.world import generate_world, step_world

SEED = 0
NODE_DENSITY = 4.0
STEPS = 4
YEARS_PER_STEP = 1_000_000
# Matches render_image.GRID_SPACING_RAD's own effective resolution at this node_density --
# the same full-sphere render-grid shape phase-0's own write-up measured against.
GRID_H, GRID_W = 801, 1601


def _phase0_world():
    world = generate_world(SEED, node_density=NODE_DENSITY, climate_density=NODE_DENSITY)
    for _ in range(STEPS):
        step_world(world, YEARS_PER_STEP)
    return world


def _full_sphere_query_points():
    lat_deg = 90.0 - (np.arange(GRID_H) + 0.5) * (180.0 / GRID_H)
    lon_deg = (np.arange(GRID_W) + 0.5) * (360.0 / GRID_W) - 180.0
    lat_grid, lon_grid = np.meshgrid(np.radians(lat_deg), np.radians(lon_deg), indexing="ij")
    return np.stack(
        [np.cos(lat_grid) * np.cos(lon_grid), np.cos(lat_grid) * np.sin(lon_grid), np.sin(lat_grid)], axis=-1
    ).reshape(-1, 3)


def test_healpix_query_is_several_times_faster_than_kdtree_at_steady_state():
    """Repeated-trial timing (phase-0's own method: single-sample timing on this box was noisy
    enough to invert the result once) of the steady-state (already-built/cached) query cost
    only -- not the one-time tree/index build, which both paths pay once per step regardless."""
    world = _phase0_world()
    all_points, _all_elev, _all_owner, kdtree = render_image._node_cloud_and_tree(world)

    world.node_cloud_resample_mode = "healpix"
    world.node_kdtree_cache = None
    _hp_points, _hp_elev, _hp_owner, node_index = render_image._node_cloud_and_tree(world)

    query_points = _full_sphere_query_points()
    kdtree.query(query_points, workers=-1)  # warm-up: first call eats page-cache/JIT cost
    node_index.query(query_points, workers=-1)

    kd_times, hp_times = [], []
    for _ in range(5):
        t0 = time.perf_counter()
        kdtree.query(query_points, workers=-1)
        kd_times.append(time.perf_counter() - t0)
        t0 = time.perf_counter()
        node_index.query(query_points, workers=-1)
        hp_times.append(time.perf_counter() - t0)

    # Generous vs. the measured ~8-9x on this implementation (phase-0's own throwaway script
    # measured ~2.9-3.8x) -- a loose floor so this doesn't flake on a slower/shared box, not a
    # tight performance assertion.
    assert min(kd_times) > max(hp_times) * 2.0


def test_node_pixel_index_scatter_and_fill_cost_is_a_small_fraction_of_a_kdtree_build():
    """The fill step was phase-0's own open question ("the one number that makes or breaks the
    rest of the plan") -- confirms it's still cheap on the real implementation, not just the
    throwaway spike, and reports the round count actually used."""
    world = _phase0_world()
    all_points, _all_elev, _all_owner, _kdtree = render_image._node_cloud_and_tree(world)
    node_count = all_points.shape[0]
    grid = healpix_grid.build(healpix_grid.nside_for_node_count(node_count))

    healpix_grid.build_node_pixel_index(grid, all_points)  # warm-up
    fill_times = []
    rounds_used = None
    for _ in range(5):
        t0 = time.perf_counter()
        index = healpix_grid.build_node_pixel_index(grid, all_points)
        fill_times.append(time.perf_counter() - t0)
        rounds_used = index.fill_rounds_used

    assert rounds_used is not None and rounds_used < 10  # phase-0 found 2-4 rounds typical
    # cKDTree build over the same node cloud is documented elsewhere in this codebase
    # (render_image._node_cloud_and_tree's own docstring) as ~20ms at this node count --
    # scatter+fill should stay in the same ballpark, not balloon into a new bottleneck.
    assert max(fill_times) < 0.2


def test_healpix_and_kdtree_resample_agree_closely_with_outliers_explained_by_local_relief():
    """The pixel-for-pixel diff check issue #133's own phase-1 checklist calls for -- not
    `array_equal` (a HEALPix pixel-ownership lookup and a true nearest-neighbour tree query are
    genuinely different operations at pixel granularity, so exact agreement isn't the bar), but
    a real accuracy characterization: overall error bounds, *and* confirmation of where the
    large outliers actually land (phase-0 could only guess "coastline/fill-boundary,
    unconfirmed"). Measured directly here: outliers correlate with local elevation gradient
    (steep terrain -- coastlines are the steepest-and-most-common case of this, but not the
    only one) rather than being literally confined to the land/ocean boundary, so the
    assertions below check gradient correlation, not a coastline-only partition."""
    world = _phase0_world()
    all_points, all_elevation, _all_owner, kdtree = render_image._node_cloud_and_tree(world)

    world.node_cloud_resample_mode = "healpix"
    world.node_kdtree_cache = None
    _hp_points, hp_elevation, _hp_owner, node_index = render_image._node_cloud_and_tree(world)

    query_points = _full_sphere_query_points()
    _, idx_kd = kdtree.query(query_points, workers=-1)
    _, idx_hp = node_index.query(query_points, workers=-1)

    elev_kd = all_elevation[idx_kd].reshape(GRID_H, GRID_W)
    elev_hp = hp_elevation[idx_hp].reshape(GRID_H, GRID_W)
    diff = np.abs(elev_kd - elev_hp)

    mean_diff, p95_diff, max_diff = diff.mean(), np.percentile(diff, 95), diff.max()
    print(f"node_cloud_resample_mode diff: mean={mean_diff:.1f}m p95={p95_diff:.1f}m max={max_diff:.1f}m")

    # Generous vs. this implementation's own measured mean ~53m / p95 ~257m / max ~7.0km --
    # bounds this implementation's own numbers produce, not phase-0's unrelated spike numbers.
    assert mean_diff < 200.0
    assert p95_diff < 750.0
    assert max_diff < 12_000.0

    gy, gx = np.gradient(elev_kd)
    local_relief = np.hypot(gx, gy)
    correlation = np.corrcoef(diff.ravel(), local_relief.ravel())[0, 1]
    print(f"corr(diff, local relief) = {correlation:.2f}")
    assert correlation > 0.3  # outliers concentrate in high-relief terrain, not scattered randomly

    # Coastlines are the special case of high relief that matters most for is_ocean
    # correctness (see hydrology's own last-step/this-step semantics note in docs/
    # profiling.md) -- confirm disagreement there is measurably worse than the grid average,
    # even though (per the gradient-correlation check above) it isn't the *only* driver.
    is_ocean_kd = elev_kd <= world.sea_level_m
    coastal = binary_dilation(is_ocean_kd, iterations=2) & ~binary_erosion(is_ocean_kd, iterations=2)
    assert diff[coastal].mean() > diff[~coastal].mean()

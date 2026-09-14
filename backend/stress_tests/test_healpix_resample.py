"""Issue #133 phases 1-2: the real-scale (not unit-test-scale) `node_cloud_resample_mode` checks.

Phase-0's own spike (docs/profiling.md's "Item 4" section, landed via PR #135) measured its
throwaway scatter+fill script against `cKDTree` at seed=0, node_density=climate_density=4.0,
~130.6K nodes -- the config reproduced here against the real (not throwaway) implementation, so
these numbers are directly comparable to phase-0's own and replace them once Phase 1 lands.
Lives in stress_tests/, not unit_tests/, because building a node_density=4.0 world and stepping
it several times is exactly the "many-step, full-simulation" shape bin/stress_test.sh's own
module docstring describes -- too slow for the fast suite's ~1s-per-test budget.

Phase 2 extends the same checks to `climate._sample_elevation_and_crust` and
`hydrology.sample_is_ocean`/`sample_is_sea` -- the latter is the accuracy gate issue #133's own
text explicitly calls for before trusting the HEALPix path for hydrology's cross-step
(last-step-node-cloud-onto-this-step-query-grid) resample, since that's exactly where
`is_ocean` correctness matters most (coastlines)."""

import time

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion

from app import climate, healpix_grid, hydrology, plates, render_image
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


# -- Issue #133 phase 2: climate._sample_elevation_and_crust ----------------------------------


def test_climate_healpix_mode_reuses_the_render_paths_shared_index_population():
    """The point of factoring `plates.cached_node_healpix_index` out of `render_image.py`'s own
    private helper (issue #133 phase 2): a same-step `climate._sample_elevation_and_crust` call
    must reuse whichever `NodePixelIndex` `render_image._node_cloud_and_tree` already built for
    the identical node cloud, not scatter+fill a second one -- the same "first caller this step
    wins" sharing `World.node_position_tree_cache` already provides for "kdtree" mode."""
    world = _phase0_world()
    world.node_cloud_resample_mode = "healpix"

    all_points, _all_elev, _all_owner, render_index = render_image._node_cloud_and_tree(world)
    assert world.node_healpix_index_cache is render_index

    world_xyz = _full_sphere_query_points().reshape(GRID_H, GRID_W, 3)
    climate._sample_elevation_and_crust(world, world_xyz, node_cloud=(all_points, list(world.plates)))
    assert world.node_healpix_index_cache is render_index  # unchanged -- climate reused it, didn't rebuild


def test_climate_healpix_and_kdtree_elevation_agree_closely():
    world = _phase0_world()
    world_xyz = _full_sphere_query_points().reshape(GRID_H, GRID_W, 3)

    world.node_cloud_resample_mode = "kdtree"
    elev_kd, ocean_kd, _lake_kd, _chan_kd, sea_kd = climate._sample_elevation_and_crust(world, world_xyz)

    world.node_cloud_resample_mode = "healpix"
    world.node_kdtree_cache = None
    world.node_position_tree_cache = None
    elev_hp, ocean_hp, _lake_hp, _chan_hp, sea_hp = climate._sample_elevation_and_crust(world, world_xyz)

    diff = np.abs(elev_kd - elev_hp)
    mean_diff, p95_diff, max_diff = diff.mean(), np.percentile(diff, 95), diff.max()
    print(f"climate elevation diff: mean={mean_diff:.1f}m p95={p95_diff:.1f}m max={max_diff:.1f}m")
    # Same bounds as the render-path check above -- same node cloud, same scatter+fill, only the
    # query grid shape differs (climate's own H/W here happens to match the render grid's).
    assert mean_diff < 200.0
    assert p95_diff < 750.0
    assert max_diff < 12_000.0
    assert (ocean_kd == ocean_hp).mean() > 0.99
    assert (sea_kd == sea_hp).mean() > 0.99


# -- Issue #133 phase 2: hydrology.sample_is_ocean/sample_is_sea (the last-step/this-step gate) -


def test_hydrology_sample_is_ocean_healpix_and_kdtree_agree_with_coastal_band_measured():
    """The specific check issue #133's own text calls for before trusting this migration:
    `hydrology.sample_is_ocean`/`sample_is_sea` resample **last step's** node cloud, one step
    stale by design -- confirm "nearest node" vs. "nearest filled HEALPix pixel" don't disagree
    materially where it matters most, coastlines, not just assume phase-1's render-path numbers
    carry over unchanged to this different (cross-step) resample."""
    world = _phase0_world()
    world_xyz = _full_sphere_query_points().reshape(GRID_H, GRID_W, 3)
    fallback = np.zeros((GRID_H, GRID_W), dtype=bool)

    world.node_cloud_resample_mode = "kdtree"
    ocean_kd = hydrology.sample_is_ocean(world, world_xyz, fallback)
    sea_kd = hydrology.sample_is_sea(world, world_xyz, fallback)

    world.node_cloud_resample_mode = "healpix"
    ocean_hp = hydrology.sample_is_ocean(world, world_xyz, fallback)
    sea_hp = hydrology.sample_is_sea(world, world_xyz, fallback)

    mismatch = ocean_kd != ocean_hp
    mismatch_rate = mismatch.mean()
    print(f"hydrology is_ocean mismatch rate: {mismatch_rate:.4f}")
    # Measured on this implementation: ~0.4% overall -- near-exact away from coastlines (see
    # docs/profiling.md's "Phase 2 landed" section for the full writeup).
    assert mismatch_rate < 0.02

    coastal = binary_dilation(ocean_kd, iterations=2) & ~binary_erosion(ocean_kd, iterations=2)
    coastal_mismatch_rate = mismatch[coastal].mean()
    interior_mismatch_rate = mismatch[~coastal].mean()
    print(f"coastal mismatch rate: {coastal_mismatch_rate:.4f}, interior: {interior_mismatch_rate:.4f}")
    # This is the real number issue #133 asked for, and it's not small: measured ~16% coastal-
    # band mismatch (vs. ~0.03% interior) at this config -- confirming, not dispelling, the
    # issue's own worry that "nearest node" and "nearest filled HEALPix pixel" diverge most
    # exactly where is_ocean correctness matters most. Bound set with headroom above the
    # measured value (not tightened to "prove" accuracy) -- see docs/profiling.md's "Phase 2
    # landed" section for the caveat this measurement supports.
    assert coastal_mismatch_rate < 0.25
    assert (sea_kd != sea_hp).mean() < 0.02

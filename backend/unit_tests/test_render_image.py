import base64
import io
import threading

import av
import numpy as np
import pytest
from PIL import Image
from app import climate, geometry, healpix_grid, hydrology, render_image
from app.world import World, generate_world, step_world


def _world(seed=1, num_plates=10, continental_fraction=0.4):
    return generate_world(seed, num_plates=num_plates, continental_fraction=continental_fraction)


def _set_one_node_channel(world, depth_m, width_m):
    """Carves a single land node's channel_depth/channel_width in place (first node of the
    first plate's first line) -- everything else in the world is untouched, so any resulting
    hillshade difference at that one node can only be this incision, not some unrelated terrain
    feature."""
    plate = world.plates[0]
    line = plate.lines[0]
    depth = line.channel_depth.copy()
    width = line.channel_width.copy()
    depth[0] = depth_m
    width[0] = width_m
    plate.replace_line(0, line.replace(channel_depth=depth, channel_width=width))


def test_hillshade_shows_a_deep_wide_channel_as_real_relief():
    # Issue #190: a deep-and-wide-enough channel should show up as genuinely lit/shadowed
    # relief via the same directional hillshade mountains already use, not a flat darken.
    world = _world()
    baseline = render_image._hillshade_for_world(world)[0]

    world.node_hillshade_cache = None
    world.node_kdtree_cache = None
    _set_one_node_channel(world, depth_m=1200.0, width_m=3000.0)
    carved = render_image._hillshade_for_world(world)[0]

    assert carved != pytest.approx(baseline, abs=1e-9)


def test_hillshade_stays_flat_for_a_deep_but_skinny_channel():
    # The issue's own explicit requirement: even a deep channel may not show up if it's skinny.
    world = _world()
    baseline = render_image._hillshade_for_world(world)[0]

    world.node_hillshade_cache = None
    world.node_kdtree_cache = None
    _set_one_node_channel(world, depth_m=1200.0, width_m=20.0)
    carved = render_image._hillshade_for_world(world)[0]

    assert carved == pytest.approx(baseline, abs=1e-9)


def test_render_grid_arrays_cover_the_sphere_with_no_gaps():
    world = _world()
    xy, elevation, plate_id, lake_depth, glacier_depth, is_volcano, is_sea, half_w, half_h, _hillshade = render_image._render_grid_arrays(world, "behrmann", np.eye(3))

    n = len(xy)
    assert n > 1000  # a real full-sphere sweep, not a token few points
    assert elevation.shape == (n,)
    assert plate_id.shape == (n,)
    assert half_w.shape == (n,)
    assert half_h.shape == (n,)
    # Every cell must have a real, positive footprint -- a zero or negative half-extent
    # would mean a hole in the map regardless of how densely the grid was swept.
    assert np.all(half_w > 0)
    assert np.all(half_h > 0)
    # Sizes must actually vary row to row (projection distortion differs by latitude) -- not
    # a single fixed value applied everywhere, which would reintroduce the original gap.
    assert len(set(half_w.tolist())) > 1
    # Every grid plate_id must reference a plate that actually exists.
    live_ids = {p.plate_id for p in world.plates}
    assert set(plate_id.tolist()) <= live_ids


def test_render_grid_arrays_match_selected_projection_shape():
    world = _world(num_plates=8)
    behrmann_xy, *_ = render_image._render_grid_arrays(world, "behrmann", np.eye(3))
    eckert4_xy, *_ = render_image._render_grid_arrays(world, "eckert4", np.eye(3))
    # Same sample count (same underlying sweep), different projected coordinates.
    assert len(behrmann_xy) == len(eckert4_xy)
    assert not np.allclose(behrmann_xy[0], eckert4_xy[0])


def _close_points(n: int, step_km: float = 5.0):
    """`n` points along the equator, `step_km` apart -- close enough together that every one
    sits within TERRAIN_RELIEF_RADIUS_KM of every other for a small `n` (used by the
    _classify_terrain_relief tests below, which need a controlled local neighbourhood)."""
    from app.plates import PLANET_RADIUS_KM

    theta = np.arange(n) * (step_km / PLANET_RADIUS_KM)
    return geometry.local_xyz(np.zeros(n), theta)


def test_classify_terrain_relief_flat_land_is_plains_plateau_not_mountain():
    from scipy.spatial import cKDTree

    points = _close_points(10)
    elevation = np.full(10, 300.0)  # flat, well above sea level, no relief anywhere
    codes = render_image._classify_terrain_relief(points, elevation, cKDTree(points), sea_level_m=0.0)
    assert np.all(codes == render_image._TERRAIN_PLAINS_PLATEAU)


def test_classify_terrain_relief_sharp_step_is_mountain_not_plains_plateau():
    from scipy.spatial import cKDTree

    points = _close_points(10)
    # Alternating land/deep-ocean elevation: every land node's neighbourhood spans ~2300 m,
    # well past MOUNTAIN_RELIEF_THRESHOLD_M -- a real cliff/escarmpent, not a gentle plain.
    elevation = np.where(np.arange(10) % 2 == 0, 300.0, -2000.0)
    codes = render_image._classify_terrain_relief(points, elevation, cKDTree(points), sea_level_m=0.0)
    is_land = elevation > 0.0
    assert np.all(codes[is_land] == render_image._TERRAIN_MOUNTAIN)
    assert np.all(codes[~is_land] == render_image._TERRAIN_NONE)


def test_classify_terrain_relief_empty_world_returns_empty():
    from scipy.spatial import cKDTree

    empty = np.zeros((0, 3))
    codes = render_image._classify_terrain_relief(empty, np.zeros(0), cKDTree(_close_points(1)), sea_level_m=0.0)
    assert codes.shape == (0,)


def test_render_grid_arrays_terrain_relief_is_opt_in():
    world = _world(num_plates=8)
    default_grid = render_image._render_grid_arrays(world, "behrmann", np.eye(3))
    assert len(default_grid) == 10  # unaffected -- every existing caller's unpack still matches

    full_grid = render_image._render_grid_arrays(world, "behrmann", np.eye(3), include_terrain_relief=True)
    assert len(full_grid) == 11
    terrain = full_grid[-1]
    elevation = full_grid[1]
    assert terrain.shape == elevation.shape
    is_land = elevation > world.sea_level_m
    assert np.all(terrain[~is_land] == render_image._TERRAIN_NONE)
    assert np.all(np.isin(terrain[is_land], [render_image._TERRAIN_MOUNTAIN, render_image._TERRAIN_PLAINS_PLATEAU]))


def test_render_png_terrain_toggles_only_change_the_elevation_view():
    world = _world(num_plates=8)
    plain = render_image.render_png(world, "behrmann", "elevation", 200, 111)
    with_mountains = render_image.render_png(world, "behrmann", "elevation", 200, 111, show_mountains=True)
    with_plains = render_image.render_png(world, "behrmann", "elevation", 200, 111, show_plains_plateaus=True)
    assert with_mountains != plain
    assert with_plains != plain

    # A view with no relief information to tint shouldn't error, and the toggles shouldn't
    # perturb it -- render_png silently ignores them outside "elevation" (see main.py's own
    # render() docstring for the same contract at the API layer).
    plates_plain = render_image.render_png(world, "behrmann", "plates", 200, 111)
    plates_toggled = render_image.render_png(world, "behrmann", "plates", 200, 111, show_mountains=True)
    assert plates_plain == plates_toggled


def test_elevation_colors_matches_known_stops():
    # Exact stops from the hypsometric table should map to their exact color.
    colors = render_image.elevation_colors(np.array([-11000.0, 0.0, 9000.0]))
    assert tuple(colors[0]) == (10, 10, 40)
    assert tuple(colors[1]) == (150, 195, 222)  # waterline: pale light blue, not sandy tan
    assert tuple(colors[2]) == (222, 217, 210)
    # Never pure white -- reserved exclusively for ice cover (GLACIER_COLOR_RGB), see
    # elevation_colors' own docstring.
    assert tuple(colors[2]) != (255, 255, 255)


def test_elevation_colors_clamps_outside_the_stop_range():
    colors = render_image.elevation_colors(np.array([-999999.0, 999999.0]))
    assert tuple(colors[0]) == (10, 10, 40)
    assert tuple(colors[1]) == (222, 217, 210)


def test_elevation_colors_shifts_with_sea_level():
    # A cell right at the new sea level should get the waterline stop's color (0m in the
    # unshifted table), the same way elevation=0 does at the default sea_level_m=0.0.
    shifted = render_image.elevation_colors(np.array([500.0]), sea_level_m=500.0)
    baseline = render_image.elevation_colors(np.array([0.0]), sea_level_m=0.0)
    assert tuple(shifted[0]) == tuple(baseline[0])


def test_channel_incision_needs_both_deep_and_wide():
    depth = np.array([0.0, 500.0, 500.0, 500.0])
    width = np.array([500.0, 0.0, 500.0, 500.0])
    incision = render_image._channel_incision_m(depth, width)
    assert incision[0] == 0.0  # no depth at all -- untouched regardless of width
    assert incision[1] == 0.0  # wide but no real depth -- untouched regardless of width
    assert incision[2] == incision[3]  # same (deep, wide) input -> same incision
    assert incision[2] > 0.0  # deep AND wide -- carves real relief


def test_channel_incision_stays_near_zero_for_a_deep_but_skinny_channel():
    # The issue's own explicit requirement: even a deep channel may not show up if it's skinny
    # -- a channel narrower than CHANNEL_INCISION_MIN_WIDTH_M incises ~nothing regardless of how
    # deep it is, so it stays invisible to hillshade the same way a real relief map sampled this
    # coarsely wouldn't resolve it either.
    deep_and_skinny = render_image._channel_incision_m(np.array([1500.0]), np.array([20.0]))
    deep_and_wide = render_image._channel_incision_m(np.array([1500.0]), np.array([5000.0]))
    assert deep_and_skinny[0] == 0.0
    assert deep_and_wide[0] > 0.0


def test_channel_incision_grows_with_depth_past_the_floor():
    modest = render_image._channel_incision_m(np.array([200.0]), np.array([500.0]))
    extreme = render_image._channel_incision_m(np.array([1800.0]), np.array([500.0]))
    assert 0.0 < modest[0] < extreme[0]


def test_plate_colors_is_stable_and_wraps():
    ids = np.array([0, 1, len(render_image.PLATE_PALETTE)])
    colors = render_image.plate_colors(ids)
    assert tuple(colors[0]) == tuple(render_image.PLATE_PALETTE[0])
    # Wraps around the palette rather than indexing out of bounds.
    assert tuple(colors[2]) == tuple(colors[0])


def test_render_png_is_decodable_at_requested_size():
    world = _world()
    for view in render_image.VIEWS:
        png = render_image.render_png(world, "behrmann", view, 320, 180)
        image = Image.open(io.BytesIO(png))
        assert image.format == "PNG"
        assert image.size == (320, 180)


def test_stream_animation_mp4_stop_event_ends_the_run_early_but_still_yields_a_video():
    world = _world()
    stop_event = threading.Event()

    # Fire the stop signal after the 2nd frame's step -- stream_animation_mp4 checks
    # `stop_event` before every individual step_fn call, so with steps_per_frame=1 this
    # should stop the run after frame 3 (frame 1 is unstepped, frames 2 and 3 each call
    # step_fn once) rather than reaching the requested 10.
    steps_taken = []

    def _counting_step(w, years):
        steps_taken.append(years)
        step_world(w, years)
        if len(steps_taken) == 2:
            stop_event.set()

    messages = list(render_image.stream_animation_mp4(
        world, "behrmann", "elevation", 64, 64, None, 1_000_000, 1, 10,
        step_fn=_counting_step, stop_event=stop_event,
    ))

    progress = [m for m in messages if m[0] == "progress"]
    assert [m[1] for m in progress] == [1, 2, 3]  # frames 1-3 completed, then the loop broke
    assert all(m[2] == 10 for m in progress)  # `total` still reports the requested ceiling
    # elapsed_years climbs by one step per frame after the first (frame 1 is the unstepped
    # starting state), and never goes backward.
    assert [m[4] for m in progress] == sorted(m[4] for m in progress)
    assert progress[-1][4] > progress[0][4]

    done = messages[-1]
    assert done[0] == "done"
    assert done[2] is True  # stopped_early
    with av.open(io.BytesIO(done[1])) as container:
        assert sum(1 for _ in container.decode(video=0)) == 3


def test_stream_animation_mp4_runs_to_completion_when_never_stopped():
    world = _world()
    messages = list(render_image.stream_animation_mp4(
        world, "behrmann", "elevation", 64, 64, None, 1_000_000, 1, 3,
    ))
    done = messages[-1]
    assert done[0] == "done"
    assert done[2] is False  # stopped_early


def test_stream_animation_mp4_steps_per_frame_steps_every_time_but_renders_only_the_last():
    # steps_per_frame=5 real step_fn calls per frame after the first, but still only 3 frames
    # rendered/encoded -- not one bigger step_fn(world, 5 * step_years) call per frame (see
    # stream_animation_mp4's own docstring for why that distinction matters).
    world = _world()
    step_calls = []

    def _counting_step(w, years):
        step_calls.append(years)
        step_world(w, years)

    messages = list(render_image.stream_animation_mp4(
        world, "behrmann", "elevation", 64, 64, None, 200_000, 5, 3,
        step_fn=_counting_step,
    ))
    progress = [m for m in messages if m[0] == "progress"]
    assert [m[1] for m in progress] == [1, 2, 3]
    assert [m[4] for m in progress] == [0.0, 1_000_000.0, 2_000_000.0]  # 5 * 200_000 per frame

    assert step_calls == [200_000] * 10  # (3 - 1) * 5 real, individual step_fn calls
    assert world.steps_taken == 10

    done = messages[-1]
    with av.open(io.BytesIO(done[1])) as container:
        assert sum(1 for _ in container.decode(video=0)) == 3  # still only 3 frames encoded


def test_stream_animation_mp4_stop_event_can_interrupt_mid_frame_batch():
    # steps_per_frame=10 but the stop event fires after the 3rd sub-step of frame 2's batch --
    # that whole batch should be abandoned (frame 2 never renders) rather than either running
    # all 10 sub-steps or rendering a partially-stepped frame.
    world = _world()
    stop_event = threading.Event()
    step_calls = []

    def _counting_step(w, years):
        step_calls.append(years)
        step_world(w, years)
        if len(step_calls) == 3:
            stop_event.set()

    messages = list(render_image.stream_animation_mp4(
        world, "behrmann", "elevation", 64, 64, None, 100_000, 10, 5,
        step_fn=_counting_step, stop_event=stop_event,
    ))
    progress = [m for m in messages if m[0] == "progress"]
    assert [m[1] for m in progress] == [1]  # only frame 1 (unstepped) ever rendered
    assert len(step_calls) == 3  # abandoned mid-batch, not run to completion or rendered

    done = messages[-1]
    assert done[2] is True  # stopped_early


def test_geomorph_colors_diverge_around_zero():
    # Neutral grey at no net change, warm where the step net-lowered a node, cool where it
    # net-raised one -- and clamped to the end stops past the +-60 m band.
    neutral, erosion, deposition = render_image.geomorph_colors(np.array([0.0, -50.0, 50.0]))
    assert tuple(neutral) == (232, 232, 232)
    assert erosion[0] > erosion[2]  # more red than blue
    assert deposition[2] > deposition[0]  # more blue than red
    lo, hi = render_image.geomorph_colors(np.array([-9999.0, 9999.0]))
    assert tuple(lo) == tuple(render_image._GEOMORPH_STOP_RGB[0].astype(int))
    assert tuple(hi) == tuple(render_image._GEOMORPH_STOP_RGB[-1].astype(int))


def test_geomorph_view_renders_neutral_before_a_step_then_varies_after():
    world = _world(seed=7, num_plates=8)
    # erosion_cache is None until the first climate/erosion step -- the view falls back to a
    # flat neutral field (plus the coastline) rather than erroring.
    before = np.asarray(Image.open(io.BytesIO(render_image.render_png(world, "behrmann", "geomorph", 320, 180))).convert("RGB"))
    assert world.erosion_cache is None

    step_world(world, years=1_000_000)
    assert world.erosion_cache is not None
    assert world.erosion_cache.net_elevation_change_m.shape == world.erosion_cache.points.shape[:1]

    after = np.asarray(Image.open(io.BytesIO(render_image.render_png(world, "behrmann", "geomorph", 320, 180))).convert("RGB"))
    # A real geomorph field has erosion and deposition both -- more than one distinct color
    # away from the background, unlike the pre-step neutral fill.
    assert len(np.unique(after.reshape(-1, 3), axis=0)) > len(np.unique(before.reshape(-1, 3), axis=0))


def test_elev_reason_view_is_all_none_before_a_step_then_gains_process_codes():
    from app import plates
    from app.elevation_lines import ELEV_CHANGE_LABELS

    world = _world(seed=7, num_plates=10)
    # Freshly generated: nothing has moved elevation yet, so every node reads NONE.
    assert np.all(plates.collect_all_elev_change_reason(world.plates) == 0.0)
    before = render_image.render_png(world, "behrmann", "elevReason", 320, 180)

    for _ in range(6):
        step_world(world, years=1_000_000)

    codes = np.unique(np.round(plates.collect_all_elev_change_reason(world.plates)).astype(int))
    # More than one process code now (at minimum an erosion/deposition/marine geomorphic one),
    # and every code is a real ELEV_CHANGE_* index.
    assert len(codes) > 1
    assert codes.max() < len(ELEV_CHANGE_LABELS)
    after = render_image.render_png(world, "behrmann", "elevReason", 320, 180)
    assert after != before


def test_elev_reason_colors_are_flat_per_code_and_clamp_out_of_range():
    from app.elevation_lines import ELEV_CHANGE_LABELS

    codes = np.arange(len(ELEV_CHANGE_LABELS), dtype=float)
    colors = render_image.elev_reason_colors(codes)
    assert colors.shape == (len(ELEV_CHANGE_LABELS), 3)
    assert len({tuple(c) for c in colors}) == len(ELEV_CHANGE_LABELS)  # every code a distinct swatch
    # a code past the table (or a fractional one from a resample) clamps, never indexes OOB
    assert tuple(render_image.elev_reason_colors(np.array([999.0]))[0]) == tuple(colors[-1])


def test_node_added_and_removed_colors_are_distinct_ramps():
    # Deliberately disjoint from each other and from _OVERLAP_AGE_STOP_RGB -- warm for added,
    # cool for removed, so the three debug views are never confusable at a glance.
    fresh_added, old_added = render_image.node_added_colors(np.array([0.0, 999.0]))
    fresh_removed, old_removed = render_image.node_removed_colors(np.array([0.0, 999.0]))
    assert tuple(fresh_added) == tuple(render_image._NODE_ADDED_STOP_RGB[0].astype(int))
    assert tuple(old_added) == tuple(render_image._NODE_ADDED_STOP_RGB[-1].astype(int))
    assert tuple(fresh_removed) == tuple(render_image._NODE_REMOVED_STOP_RGB[0].astype(int))
    assert tuple(old_removed) == tuple(render_image._NODE_REMOVED_STOP_RGB[-1].astype(int))
    assert fresh_added[0] > fresh_added[2]  # added ramp reads warm (more red than blue)
    assert fresh_removed[2] > fresh_removed[0]  # removed ramp reads cool (more blue than red)


def test_overlap_age_view_renders_gap_age_dots_once_tracked():
    from app import gaps

    world = _world(seed=7, num_plates=8)
    world.plates.pop(1)  # opens a real gap
    gaps.reconcile_gap_tracks(world)
    assert len(world.gap_tracks) >= 1

    png = render_image.render_png(world, "behrmann", "overlapAge", 320, 180)
    pixels = np.asarray(Image.open(io.BytesIO(png)).convert("RGB")).reshape(-1, 3)
    # A freshly-tracked gap's age is exactly 0 -- the gap-age ramp's own first stop colour --
    # so it must appear verbatim among the rendered pixels, not just "some new colour".
    assert np.any(np.all(pixels == render_image._GAP_AGE_STOP_RGB[0].astype(int), axis=1))


def test_node_age_view_renders_more_once_points_have_been_added_and_removed():
    world = _world(seed=7, num_plates=8)
    # Freshly generated: every node predates tracking (-1.0 sentinel) and nothing has been
    # removed yet -- an all-backdrop render, same "healthy/quiet" convention overlapAge uses.
    before = np.asarray(Image.open(io.BytesIO(render_image.render_png(world, "behrmann", "nodeAge", 320, 180))).convert("RGB"))

    for _ in range(6):
        step_world(world, years=1_000_000)

    from app import plates as plates_module

    created = plates_module.collect_all_node_created_years(world.plates)
    assert np.any(created >= 0.0), "expected at least one node to have been created by boundary growth"

    after = np.asarray(Image.open(io.BytesIO(render_image.render_png(world, "behrmann", "nodeAge", 320, 180))).convert("RGB"))
    assert len(np.unique(after.reshape(-1, 3), axis=0)) > len(np.unique(before.reshape(-1, 3), axis=0))


def test_combined_view_encodes_biome_ids_in_the_alpha_channel():
    # Combined's per-pixel class id rides in alpha (see render_image.COMBINED_LAKE_ID_CODE's
    # comment): alpha = 255 - code, code 0 only for gaps between cells, biome_id + 1 for every
    # classified land (Köppen) or ocean (pelagic) cell, lake/glacier overlays above that.
    world = _world()
    png = render_image.render_png(world, "behrmann", "combined", 320, 180)
    image = Image.open(io.BytesIO(png))
    assert image.mode == "RGBA"

    alpha = np.asarray(image.convert("RGBA"))[:, :, 3]
    codes = 255 - alpha.astype(int)
    assert codes.min() >= 0
    assert codes.max() <= render_image.COMBINED_GLACIER_ID_CODE
    # Both land and ocean carry real codes now (ocean is no longer code 0).
    assert np.count_nonzero(codes) > 0.5 * codes.size

    # Other views stay plain RGB -- alpha is a Combined-only channel.
    elev = Image.open(io.BytesIO(render_image.render_png(world, "behrmann", "elevation", 320, 180)))
    assert elev.mode == "RGB"


def _pile_ice_everywhere(world, depth_m):
    for plate in world.plates:
        for i, line in enumerate(plate.lines):
            gd = line.glacier_depth.copy()
            gd[:] = depth_m
            plate.replace_line(i, line.replace(glacier_depth=gd))


@pytest.mark.parametrize("view", ["elevation", "combined"])
def test_glacier_overlay_is_shaded_by_hillshade_not_flat(view):
    # Regression for issue #169: ice caps read as flat, texture-less plateaus even over
    # genuinely varied terrain, because the glacier overlay painted a single flat
    # GLACIER_COLOR_RGB with no relation to the hillshade the land underneath it would show.
    # Glaciating the whole world (real, naturally rugged terrain from generation) means almost
    # every pixel is "under ice" with no glacier/non-glacier boundary to blur across, so any
    # color variation among those pixels can only come from hillshade, not edge antialiasing.
    world = _world(seed=7, num_plates=8, continental_fraction=0.5)
    _pile_ice_everywhere(world, hydrology.GLACIER_VISIBLE_DEPTH_M + 500.0)

    png = render_image.render_png(world, "behrmann", view, 320, 180)
    pixels = np.asarray(Image.open(io.BytesIO(png)).convert("RGB"))

    glacier_rgb = np.array(render_image.GLACIER_COLOR_RGB)
    near_glacier = np.all(np.abs(pixels.astype(int) - glacier_rgb) < 60, axis=-1)
    assert near_glacier.sum() > 500  # a real, sizeable ice-covered area to sample from

    # A handful of ice-adjacent pixels differ from the flat overlay color regardless -- blur
    # and the river-line overlay both blend a few pixels toward it without touching the code
    # under test (confirmed by running this same check against the pre-fix code: it produces
    # a couple hundred such incidental pixels out of ~40k). Real per-cell hillshade variation
    # is far larger -- roughly a third of the ice-covered area -- so the bar is set well above
    # that incidental-blend noise floor.
    not_flat = near_glacier & ~np.all(np.abs(pixels.astype(int) - glacier_rgb) < 3, axis=-1)
    assert not_flat.sum() > 2000


def test_land_only_glacier_shade_keeps_floating_ice_flat():
    hillshade = np.array([0.5, 1.5, 0.7, 1.2])
    is_ocean = np.array([True, False, False, False])
    is_lake = np.array([False, True, False, False])
    shade = render_image._land_only_glacier_shade(hillshade, is_ocean, is_lake)
    # Floating sea ice and frozen-lake ice: flat (1.0) regardless of the seafloor/lakebed
    # hillshade underneath -- the real values above (0.5, 1.5) would visibly darken/brighten
    # it if the guard were missing.
    assert shade[0] == 1.0
    assert shade[1] == 1.0
    # Land-based ice: real per-cell hillshade passes through unchanged.
    assert shade[2] == 0.7
    assert shade[3] == 1.2


def test_combined_view_ocean_ice_stays_flat_not_shaded_by_seafloor_relief():
    # Companion regression to the one above: a real orbital photo shows floating sea ice as a
    # flat white surface regardless of the seafloor relief beneath it -- the exact rule the
    # Combined view's own land_rgb hillshade already carves out for open ocean (see its
    # "land only" comment). Glaciating ocean nodes only (never land) means any color variation
    # among the resulting ice pixels could only come from that seafloor bathymetry leaking
    # through, which the fix explicitly guards against.
    world = _world(seed=11, num_plates=8, continental_fraction=0.5)
    depth = hydrology.GLACIER_VISIBLE_DEPTH_M + 500.0
    for plate in world.plates:
        for i, line in enumerate(plate.lines):
            gd = line.glacier_depth.copy()
            gd[line.elevation <= world.sea_level_m] = depth
            plate.replace_line(i, line.replace(glacier_depth=gd))

    png = render_image.render_png(world, "behrmann", "combined", 320, 180)
    pixels = np.asarray(Image.open(io.BytesIO(png)).convert("RGB"))

    glacier_rgb = np.array(render_image.GLACIER_COLOR_RGB)
    near_glacier = np.all(np.abs(pixels.astype(int) - glacier_rgb) < 60, axis=-1)
    assert near_glacier.sum() > 500  # a real, sizeable sea-ice-covered area to sample from

    not_flat = near_glacier & ~np.all(np.abs(pixels.astype(int) - glacier_rgb) < 3, axis=-1)
    # Unlike the sibling test above, land here is never glaciated, so ice-covered ocean cells
    # sit right next to varied-color land at every coastline -- the post-fill blur (see
    # CELL_BLUR_RADIUS_PX) blends a real, if modest, band of those coastline pixels toward
    # off-white regardless of this fix (confirmed directly: ~10% on this seed with the fix
    # applied, entirely from that blur). A real seafloor-hillshade leak would push this far
    # higher across the *whole* ocean interior, not just coastline blur -- confirmed directly
    # against the pre-fix code, which put it above 30% on this same seed -- so the bar sits
    # well clear of both.
    assert not_flat.sum() < 0.2 * near_glacier.sum()


def test_biome_view_smoothing_preserves_the_major_biomes_and_barely_moves_the_rest():
    # smooth_biome_field is a cleanup pass, not a reclassification: on the real biome render
    # grid it should change only a small slice of land and never erase a biome that has a
    # genuine regional presence.
    from app import biomes

    world = _world(seed=7, num_plates=12, continental_fraction=0.6)
    lat_deg, _lon, _xyz, elevation_m, is_ocean, air_temp, ocean_temp, precip, _lake, glacier_depth, _is_sea, _hillshade = render_image._biome_fields(
        world, *render_image.biome_grid_dimensions(world.climate_density)
    )
    display_temp = np.where(is_ocean, ocean_temp, air_temp)
    slope = biomes.grid_slope(elevation_m, lat_deg)
    lat_grid = np.broadcast_to(lat_deg[:, None], elevation_m.shape)
    # Same geometry inputs smooth_biome_field derives internally, so the diff isolates the
    # vote pass rather than the continentality/coast-distance inputs.
    raw = biomes.classify_biomes(
        display_temp, precip, elevation_m, slope, is_ocean, world.sea_level_m,
        lat_deg=lat_grid, axial_tilt_deg=world.axial_tilt_deg,
        continentality=biomes.grid_continentality(is_ocean),
        dist_to_land_rad=biomes.grid_dist_to_land_rad(is_ocean),
        has_sea_ice=is_ocean & (glacier_depth > 0.0),
    )
    smoothed = biomes.smooth_biome_field(
        display_temp, precip, elevation_m, slope, is_ocean, world.sea_level_m,
        lat_deg=lat_deg, axial_tilt_deg=world.axial_tilt_deg, glacier_depth_m=glacier_depth,
    )

    land = ~is_ocean
    n_land = int(land.sum())
    assert int(np.count_nonzero(raw[land] != smoothed[land])) < 0.06 * n_land

    for biome_id in np.unique(raw[land]):
        raw_share = np.count_nonzero(raw[land] == biome_id) / n_land
        if raw_share >= 0.02:
            assert np.count_nonzero(smoothed[land] == biome_id) / n_land >= 0.01


def test_biome_grid_dimensions_matches_reference_at_density_one():
    assert render_image.biome_grid_dimensions(1.0) == (render_image.BIOME_GRID_HEIGHT, render_image.BIOME_GRID_WIDTH)


def test_biome_grid_dimensions_doubles_each_dimension_at_density_two():
    # Within 1 of an exact double, not necessarily exact -- BIOME_GRID_HEIGHT/WIDTH are
    # themselves already round()-derived from an irrational spacing ratio, so a second
    # independent round() at half the spacing doesn't always commute perfectly with "times 2"
    # (confirmed directly: BIOME_GRID_WIDTH=400 but biome_grid_dimensions(2.0)'s width is 801,
    # not 800). Unlike climate.grid_dimensions, whose GRID_HEIGHT/WIDTH reference values are
    # already nice round numbers, so that doubling happens to land exactly.
    height, width = render_image.biome_grid_dimensions(2.0)
    assert abs(height - render_image.BIOME_GRID_HEIGHT * 2) <= 1
    assert abs(width - render_image.BIOME_GRID_WIDTH * 2) <= 1


def test_render_png_is_decodable_at_requested_size_with_doubled_climate_density():
    # Same sweep as test_render_png_is_decodable_at_requested_size, but at
    # World.climate_density=2.0 -- confirms every view (climate-grid-derived and
    # biome-grid-derived alike) still renders correctly at the finer resolution, not just the
    # default one.
    world = _world()
    world.climate_density = 2.0
    for view in render_image.VIEWS:
        png = render_image.render_png(world, "behrmann", view, 320, 180)
        image = Image.open(io.BytesIO(png))
        assert image.format == "PNG"
        assert image.size == (320, 180)


def test_render_png_empty_world_returns_background_only():
    world = World(seed=1)  # plates defaults to [] -- shouldn't happen via the API, but shouldn't crash
    png = render_image.render_png(world, "behrmann", "elevation", 100, 60)
    image = Image.open(io.BytesIO(png)).convert("RGB")
    pixels = np.asarray(image)
    assert np.all(pixels == np.array(render_image.BACKGROUND_RGB))


def test_render_png_base64_round_trips():
    world = _world()
    encoded = render_image.render_png_base64(world, "behrmann", "elevation", 100, 60)
    decoded = base64.b64decode(encoded)
    image = Image.open(io.BytesIO(decoded))
    assert image.size == (100, 60)


def test_draw_rivers_only_draws_segments_above_the_width_floor():
    # The main map views (unlike the River Inspector, which deliberately shows every
    # is_river-classified network regardless of size) additionally require channel_width to
    # clear the absolute RIVER_VISIBLE_MIN_WIDTH_M creek floor -- confirmed here by drawing the
    # exact same is_river/flow_target topology twice, once with channel_width below the floor
    # (no river pixels should appear at all) and once well above it (some must). Calls
    # _draw_rivers directly (rather than the full render_png) against a blank background so
    # any pixel that isn't BACKGROUND_RGB can only be the river line itself (now antialiased --
    # see RIVER_BLUR_RADIUS_PX -- so drawn pixels blend toward RIVER_COLOR_RGB rather than
    # matching it exactly).
    world = _world()
    world.node_density = 1.0
    world.climate_density = 1.0
    points = np.array([[1.0, 0.0, 0.0], [0.99, 0.05, 0.0], [0.0, 1.0, 0.0], [0.0, 0.95, 0.05]])
    points = points / np.linalg.norm(points, axis=1, keepdims=True)
    n = len(points)
    base_fields = dict(
        points=points,
        elevation=np.array([10.0, 5.0, 10.0, 5.0]),
        is_ocean=np.zeros(n, dtype=bool),
        neighbor_idx=np.zeros((n, 1), dtype=np.int64),
        flow_accum=np.zeros(n),
        water_deposited=np.zeros(n),
        filled_elevation=np.zeros(n),
        spill_target=np.full(n, -1, dtype=np.int64),
        is_river=np.array([True, False, True, False]),
        flow_target=np.array([1, -1, 3, -1]),
        lake_depth=np.zeros(n),
        glacier_depth=np.zeros(n),
        plates_in_order=[],
    )

    def river_pixel_count(channel_width):
        world.hydrology_cache = hydrology.HydrologyFields(channel_width=channel_width, **base_fields)
        image = Image.new("RGB", (300, 200), render_image.BACKGROUND_RGB)
        image = render_image._draw_rivers(image, world, "behrmann", 50.0, 150.0, 100.0, 1.0, np.eye(3))
        pixels = np.asarray(image)
        return int(np.any(pixels != np.array(render_image.BACKGROUND_RGB), axis=-1).sum())

    below_floor = river_pixel_count(np.full(n, render_image.RIVER_VISIBLE_MIN_WIDTH_M * 0.5))
    above_floor = river_pixel_count(np.full(n, render_image.RIVER_VISIBLE_MIN_WIDTH_M * 5.0))

    assert below_floor == 0
    assert above_floor > 0


def test_draw_rivers_wide_segment_has_a_dimmer_halo_than_its_own_core(monkeypatch):
    # A multi-pixel-wide river should show a real cross-sectional gradient -- a core close to
    # full RIVER_COLOR_RGB with a dimmer halo around it (RIVER_CORE_WIDTH_FRACTION/
    # RIVER_HALO_FILL_FRACTION) -- not one flat, uniformly-blended slab across its whole width.
    # _rivers_to_draw is monkeypatched to a single segment at a known, full color_alpha/
    # width_frac so this is a clean test of the halo/core drawing mechanism itself, not
    # entangled with the percentile-rank color scheme (which needs a large candidate pool to
    # behave sensibly and isn't what's under test here).
    world = _world()
    points = np.array([[1.0, 0.0, 0.0], [0.99, 0.05, 0.0]])
    points = points / np.linalg.norm(points, axis=1, keepdims=True)
    n = len(points)
    world.hydrology_cache = hydrology.HydrologyFields(
        points=points,
        elevation=np.array([10.0, 5.0]),
        is_ocean=np.zeros(n, dtype=bool),
        neighbor_idx=np.zeros((n, 1), dtype=np.int64),
        flow_accum=np.zeros(n),
        water_deposited=np.zeros(n),
        filled_elevation=np.zeros(n),
        spill_target=np.full(n, -1, dtype=np.int64),
        is_river=np.array([True, False]),
        flow_target=np.array([1, -1]),
        lake_depth=np.zeros(n),
        glacier_depth=np.zeros(n),
        plates_in_order=[],
        channel_width=np.array([5000.0, 0.0]),
    )
    monkeypatch.setattr(render_image, "_rivers_to_draw", lambda w: (np.array([0]), np.array([1.0]), np.array([1.0])))

    image = Image.new("RGB", (300, 200), render_image.BACKGROUND_RGB)
    image = render_image._draw_rivers(image, world, "behrmann", 50.0, 150.0, 100.0, 10.0, np.eye(3))
    pixels = np.asarray(image).astype(np.int64)
    background = np.array(render_image.BACKGROUND_RGB)
    river_color = np.array(render_image.RIVER_COLOR_RGB)

    is_river_px = np.any(pixels != background, axis=-1)
    assert is_river_px.sum() > 0
    near_full_color = is_river_px & np.all(np.abs(pixels - river_color) < 15, axis=-1)
    # A real halo ring: plenty of drawn pixels are neither background nor near-full color --
    # far more than the handful of pixels a purely edge-blur effect alone would leave, which is
    # what an all-one-fill single-pass line (see the width-floor test above) would show instead.
    halo_ring_px = int(is_river_px.sum()) - int(near_full_color.sum())
    assert halo_ring_px > int(near_full_color.sum()) * 0.2


def _linear_rivers_world(mouth_widths):
    """A synthetic world whose hydrology_cache holds one independent 3-node linear river per
    entry in `mouth_widths` (head -> mid -> mouth -> ocean sink). Node i of river r is index
    3*r + i; channel_width along each river is 20% / 50% / 100% of that river's own mouth
    width, so it tapers from head to mouth. flow_accum mirrors channel_width (a real world's
    two are correlated too, just channel_width also saturates at a hard cap over time in a way
    flow_accum doesn't -- irrelevant for these small, freshly-built fixtures) so tests can
    exercise _rivers_to_draw's flow-percentile width gate without a separate parameter. Ocean
    sink nodes come last."""
    world = _world()
    world.node_density = 1.0
    world.climate_density = 1.0
    r = len(mouth_widths)
    n = 3 * r + 1
    ocean = n - 1

    points = np.zeros((n, 3))
    points[:, 0] = 1.0
    points[:, 1] = np.linspace(-0.4, 0.4, n)  # spread them out so segments are real hops
    points = points / np.linalg.norm(points, axis=1, keepdims=True)

    flow_target = np.full(n, -1, dtype=np.int64)
    channel_width = np.zeros(n)
    is_river = np.zeros(n, dtype=bool)
    is_ocean = np.zeros(n, dtype=bool)
    is_ocean[ocean] = True
    for ri_, mw in enumerate(mouth_widths):
        head, mid, mouth = 3 * ri_, 3 * ri_ + 1, 3 * ri_ + 2
        flow_target[head], flow_target[mid], flow_target[mouth] = mid, mouth, ocean
        channel_width[head], channel_width[mid], channel_width[mouth] = 0.2 * mw, 0.5 * mw, mw
        is_river[head] = is_river[mid] = is_river[mouth] = True

    world.hydrology_cache = hydrology.HydrologyFields(
        points=points,
        elevation=np.linspace(50.0, 0.0, n),
        is_ocean=is_ocean,
        neighbor_idx=np.zeros((n, 1), dtype=np.int64),
        water_deposited=np.zeros(n),
        filled_elevation=np.zeros(n),
        spill_target=np.full(n, -1, dtype=np.int64),
        is_river=is_river,
        flow_target=flow_target,
        flow_accum=channel_width.copy(),
        lake_depth=np.zeros(n),
        glacier_depth=np.zeros(n),
        plates_in_order=[],
        channel_width=channel_width,
    )
    return world


def test_rivers_to_draw_treats_a_stale_pre_channel_width_cache_as_absent():
    # A hydrology_cache loaded from a save written before channel_width existed (or any other
    # hand-built fixture that doesn't set it) must be treated the same as "no cache at all" --
    # not crash, and not draw anything -- the same stale-cache guard convention
    # hydrology.sample_is_sea already uses for its own is_sea field.
    world = _linear_rivers_world([5000])

    # Simulate a pre-#190 cache: channel_width defaults to a shape-0 array, mismatched against
    # the real (nonzero) node count.
    stale = world.hydrology_cache
    stale.channel_width = np.zeros(0)
    src_idx, color_alpha, width_frac = render_image._rivers_to_draw(world)
    assert len(src_idx) == 0
    assert len(color_alpha) == 0
    assert len(width_frac) == 0

    # No hydrology_cache at all (before the world's first step) is the same case.
    world.hydrology_cache = None
    src_idx2, color_alpha2, width_frac2 = render_image._rivers_to_draw(world)
    assert len(src_idx2) == 0
    assert len(color_alpha2) == 0
    assert len(width_frac2) == 0


def test_rivers_to_draw_draws_every_wide_enough_network_regardless_of_count():
    # Unlike the old top-N-networks-by-flow scheme, there's no cap on how many distinct
    # drainage networks get drawn -- every one whose channel clears the width floor shows up,
    # including a 9th, 10th, ... river that an old fixed network-count cap would have dropped.
    mouth_widths = [4800, 4700, 4600, 4500, 4400, 4300, 4200, 4100, 4000]
    world = _linear_rivers_world(mouth_widths)
    src_idx, _, _ = render_image._rivers_to_draw(world)
    drawn_rivers = {int(i) // 3 for i in src_idx}
    assert drawn_rivers == set(range(len(mouth_widths)))

    # A genuinely narrow/creek-scale river (well under RIVER_VISIBLE_MIN_WIDTH_M) among the
    # same set is still excluded -- the width floor still does real work, it's just no longer a
    # fixed count of networks.
    world_with_creek = _linear_rivers_world(mouth_widths + [30])
    src_idx2, _, _ = render_image._rivers_to_draw(world_with_creek)
    assert len(mouth_widths) not in {int(i) // 3 for i in src_idx2}


def test_rivers_to_draw_color_is_a_percentile_rank_not_an_absolute_fraction():
    # Color blend fraction is a percentile rank of channel_width among this call's own visible
    # drawn segments, not a fixed fraction of the physical MAX_CHANNEL_WIDTH_M cap (which most
    # real rivers eventually saturate to, reading as uniformly bright) or of each segment's own
    # network's widest point (which, since is_river already only selects each network's own
    # highest-flow tail, also reads as mostly-saturated) -- confirmed here across three
    # independent rivers of very different absolute scale.
    world = _linear_rivers_world([5000, 2500, 500])
    src_idx, color_alpha, _ = render_image._rivers_to_draw(world)
    river_of = np.array([int(i) // 3 for i in src_idx])

    # This call's single widest segment overall (river 0's mouth, channel_width == 5000, the
    # unique global max among all 9 drawn nodes) reaches exactly full color; its single
    # narrowest (river 2's head, channel_width == 100, the unique global min) sits exactly at
    # the floor.
    assert color_alpha[river_of == 0].max() == pytest.approx(1.0)
    assert color_alpha[river_of == 2].min() == pytest.approx(render_image.RIVER_COLOR_BLEND_MIN_FRACTION)

    # Monotonic within each individual tapering river (head -> mid -> mouth) regardless of how
    # the three rivers' absolute scales compare to each other -- the "gets bluer toward a
    # confluence/mouth" behavior.
    for r in (0, 1, 2):
        segs = color_alpha[river_of == r]
        assert segs[0] < segs[1] < segs[2]

    # A much bigger river's segments generally rank higher than a much smaller river's -- the
    # scheme isn't simply "every river looks equally blue regardless of its own scale."
    assert color_alpha[river_of == 0].mean() > color_alpha[river_of == 2].mean()


def test_rivers_to_draw_width_only_widens_for_a_steep_top_flow_percentile():
    # Line width is gated on a *flow_accum* percentile rank (RIVER_WIDE_MIN_PERCENTILE), not
    # channel_width -- and deliberately a much steeper bar than color's own rank, so only the
    # single biggest river's own mouth (this fixture's global top-flow node) ever draws wider
    # than 1px, not every node along it.
    world = _linear_rivers_world([5000, 500])
    src_idx, _, width_frac = render_image._rivers_to_draw(world)
    river_of = np.array([int(i) // 3 for i in src_idx])

    # The big river's own mouth (this fixture's single highest-flow node) gets real extra
    # width...
    assert width_frac[river_of == 0].max() > 0.0
    # ...but nothing else does -- not even that same big river's own mid-course, and not the
    # small river anywhere along its length, despite both ranking reasonably on the (much more
    # generous) percentile color scale.
    assert width_frac[river_of == 0][:-1].max() == 0.0
    assert width_frac[river_of == 1].max() == 0.0


def test_render_png_scales_visual_constants_with_resolution():
    """Doubling the requested width (the "sharper, same displayed size" retina use case)
    should double pixel_scale, so a fixed-size feature like the pole marker should occupy
    roughly proportionally more pixels at 2x than at 1x -- i.e. the map doesn't get
    thinner-looking lines/markers just because more pixels were requested."""
    world = _world(num_plates=6, continental_fraction=1.0)
    small = render_image.render_png(world, "behrmann", "plates", 550, 306)
    large = render_image.render_png(world, "behrmann", "plates", 1100, 611)

    def non_background_fraction(png_bytes, size):
        image = Image.open(io.BytesIO(png_bytes)).convert("RGB").resize(size, Image.NEAREST)
        pixels = np.asarray(image)
        return np.mean(np.any(pixels != np.array(render_image.BACKGROUND_RGB), axis=-1))

    # Downscale both to the same size before comparing -- if line/marker widths scaled
    # correctly, the two should cover a similar fraction of non-background pixels.
    frac_small = non_background_fraction(small, (275, 153))
    frac_large = non_background_fraction(large, (275, 153))
    assert frac_small > 0
    assert frac_large > 0
    assert abs(frac_small - frac_large) < 0.05


def test_plate_tectonics_pole_is_the_true_euler_pole():
    """Euler poles can be anywhere on the map, not necessarily near the plate they belong
    to (this is physically normal for real plate tectonics too) -- pole_xyz should be
    exactly +omega/|omega|, with no adjustment toward the plate's own territory."""
    world = _world(seed=9, num_plates=12, continental_fraction=0.7)
    for plate in world.plates:
        speed = np.linalg.norm(plate.omega)
        if speed < 1e-15:
            continue
        info = render_image._plate_tectonics("eckert4", plate, np.eye(3))
        assert np.allclose(info["rotation_arc"]["pole_xyz"], plate.omega / speed)


def test_pole_marker_uses_plate_color_not_a_fixed_color():
    """Since a plate's true Euler pole can project anywhere on the map, color -- not
    position -- is what ties a pole marker back to its plate (see render_png); confirms the
    marker no longer uses a fixed color regardless of which plate it belongs to."""
    world = _world(seed=3, num_plates=6, continental_fraction=1.0)
    png = render_image.render_png(world, "eckert4", "plates", 600, 400)
    image = Image.open(io.BytesIO(png)).convert("RGB")
    pixels = np.asarray(image).reshape(-1, 3)
    old_fixed_pole_color = np.array([255, 45, 85])
    assert not np.any(np.all(pixels == old_fixed_pole_color, axis=-1))


def test_climate_views_render_as_distinct_images():
    world = _world(seed=4, num_plates=10, continental_fraction=0.5)
    pngs = {view: render_image.render_png(world, "behrmann", view, 320, 180) for view in render_image.CLIMATE_VIEWS}
    views = list(pngs)
    for i in range(len(views)):
        for j in range(i + 1, len(views)):
            assert pngs[views[i]] != pngs[views[j]], f"{views[i]} and {views[j]} rendered identically"


def test_ocean_currents_view_marks_swells_at_synthetic_convergence(monkeypatch):
    # The oceanCurrents view draws swells at climate.ClimateFields' own swell_rows/swell_cols
    # (climate.compute_ocean_swells' picked convergence cells) resampled to xyz. Monkeypatch
    # that pick to a known grid cell and confirm render_png actually draws a white marker
    # there -- the same "does the drawing step work" contract as before.
    from app import climate

    monkeypatch.setattr(climate, "compute_ocean_swells", lambda *a, **k: (np.array([20]), np.array([40])))

    png = render_image.render_png(_world(), "behrmann", "oceanCurrents", 320, 180)
    image = Image.open(io.BytesIO(png)).convert("RGB")
    pixels = np.asarray(image).reshape(-1, 3)
    swell_marker_white = np.array([255, 255, 255])
    assert np.any(np.all(pixels == swell_marker_white, axis=-1))


def _latlon_grid_points(n=12, spacing_deg=0.4):
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    pts = geometry.latlon_to_xyz(np.radians(ii.reshape(-1) * spacing_deg), np.radians(jj.reshape(-1) * spacing_deg))
    return pts, ii.reshape(-1), jj.reshape(-1)


def test_coastal_dither_fraction_flags_isolated_specks_but_not_a_coherent_coast():
    pts, ii, jj = _latlon_grid_points()

    # A mostly-ocean shelf with a scattering of lone, well-separated land nodes -- exactly the
    # single-pixel islands the investigation cared about. Each speck's whole neighbourhood is
    # the opposite class, so it pegs at 1.0 and clears the flag threshold.
    speck_elev = np.full(pts.shape[0], -20.0)
    speck_mask = (ii % 4 == 1) & (jj % 4 == 1)
    speck_elev[speck_mask] = 20.0
    speck_frac, speck_near = render_image.coastal_dither_fraction(pts, speck_elev, 0.0)
    assert speck_near.all()  # every |elev| = 20 < SPECKLE_NEAR_BAND_M
    assert np.all(speck_frac[speck_mask] >= 0.99)
    assert int((speck_frac >= render_image.SPECKLE_FLAG_FRACTION).sum()) == int(speck_mask.sum())

    # A gentle monotonic ramp across sea level is a coherent shoreline: only nodes straddling
    # the waterline see any disagreement at all, and none of it reaches the flag threshold.
    ramp_frac, ramp_near = render_image.coastal_dither_fraction(pts, (ii - 5.5) * 8.0, 0.0)
    assert not (ramp_frac >= render_image.SPECKLE_FLAG_FRACTION).any()
    assert speck_frac[speck_mask].max() > ramp_frac[ramp_near].max() + 0.3


def test_coastal_dither_fraction_is_zero_outside_the_near_band():
    pts, ii, jj = _latlon_grid_points()
    elev = np.where((ii + jj) % 2 == 0, 5000.0, -5000.0)  # a checkerboard, but nowhere near sea level
    frac, near = render_image.coastal_dither_fraction(pts, elev, 0.0)
    assert not near.any()
    assert np.all(frac == 0.0)


def test_speckle_view_differs_from_the_elevation_view():
    world = _world()
    assert render_image.render_png(world, "behrmann", "speckle", 320, 180) != render_image.render_png(
        world, "behrmann", "elevation", 320, 180
    )


def test_speckle_view_draws_flagged_nodes_in_the_flag_color(monkeypatch):
    # Feed the renderer a synthetic per-node fraction so a known slice of nodes clears
    # SPECKLE_FLAG_FRACTION -- those must show up as the oversized magenta flag marker.
    world = _world()
    all_points, _elev, _owner = render_image.plates.collect_all_points(world.plates)
    n = len(all_points)
    monkeypatch.setattr(
        render_image, "coastal_dither_fraction", lambda *a, **k: (np.linspace(0.0, 1.0, n), np.ones(n, dtype=bool))
    )
    png = render_image.render_png(world, "behrmann", "speckle", 500, 275)
    pixels = np.asarray(Image.open(io.BytesIO(png)).convert("RGB")).reshape(-1, 3)
    assert np.any(np.all(pixels == np.array(render_image.SPECKLE_FLAG_RGB), axis=-1))


def test_combined_view_draws_sea_color_for_a_sea_tier_lake_not_plain_lake_color(monkeypatch):
    # Force every node to read as "a lake" (plates.collect_all_lake_depth) and sea tier
    # (hydrology.sample_is_sea) -- easy to find on the rendered image, and isolates the
    # color-selection logic itself from needing a real, correctly-classified sea-scale basin.
    world = _world()
    all_points, _elev, _owner = render_image.plates.collect_all_points(world.plates)
    lake_depth = np.full(len(all_points), 50.0)
    monkeypatch.setattr(render_image.plates, "collect_all_lake_depth", lambda plates: lake_depth)
    monkeypatch.setattr(render_image.hydrology, "sample_is_sea", lambda world, xyz, fallback: np.ones(fallback.shape, dtype=bool))

    png = render_image.render_png(world, "behrmann", "combined", 320, 180)
    pixels = np.asarray(Image.open(io.BytesIO(png)).convert("RGB")).reshape(-1, 3)
    assert np.any(np.all(pixels == np.array(render_image.SEA_COLOR_RGB), axis=-1))
    assert not np.any(np.all(pixels == np.array(render_image.LAKE_COLOR_RGB), axis=-1))


def test_combined_view_draws_plain_lake_color_when_not_sea_tier(monkeypatch):
    world = _world()
    all_points, _elev, _owner = render_image.plates.collect_all_points(world.plates)
    lake_depth = np.full(len(all_points), 50.0)
    monkeypatch.setattr(render_image.plates, "collect_all_lake_depth", lambda plates: lake_depth)
    monkeypatch.setattr(render_image.hydrology, "sample_is_sea", lambda world, xyz, fallback: np.zeros(fallback.shape, dtype=bool))

    png = render_image.render_png(world, "behrmann", "combined", 320, 180)
    pixels = np.asarray(Image.open(io.BytesIO(png)).convert("RGB")).reshape(-1, 3)
    assert np.any(np.all(pixels == np.array(render_image.LAKE_COLOR_RGB), axis=-1))
    assert not np.any(np.all(pixels == np.array(render_image.SEA_COLOR_RGB), axis=-1))


def test_rotate_maps_a_known_point_to_its_expected_position():
    """The core operation the whole view-rotation feature rests on: a 180-degree rotation
    about the z-axis should send lat=0/lon=0 to lat=0/lon=180 (its antipode on the equator)."""
    rotation = geometry.rotation_matrix(np.array([0.0, 0.0, 1.0]), np.pi)
    origin = geometry.latlon_to_xyz(np.array([0.0]), np.array([0.0]))
    rotated = render_image._rotate(origin, rotation)
    lat, lon = geometry.xyz_to_latlon(rotated)
    assert np.isclose(lat[0], 0.0, atol=1e-9)
    assert np.isclose(abs(lon[0]), np.pi, atol=1e-9)


def test_render_grid_stays_gap_free_under_a_nontrivial_rotation():
    """The per-cell half-extent fix (see _render_grid_arrays) exists specifically so rotation
    doesn't reopen the gaps the render grid was built to avoid -- re-run the original
    no-gaps assertions at a rotation that mixes all three axes, not just identity/90/180."""
    world = _world(seed=8)
    rotation = geometry.rotation_matrix(np.array([0.4, -0.5, 0.7]), 1.3)
    xy, elevation, plate_id, lake_depth, glacier_depth, is_volcano, is_sea, half_w, half_h, _hillshade = render_image._render_grid_arrays(world, "eckert4", rotation)
    assert len(xy) > 1000
    assert np.all(half_w > 0)
    assert np.all(half_h > 0)
    assert len(set(half_w.round(6).tolist())) > 1  # genuinely per-cell now, not one value per row


def test_climate_grid_stays_gap_free_under_a_nontrivial_rotation():
    """A rotation reorients the oval map's *content*, not its outer shape (Eckert4's own
    coverage is unchanged), so identity and rotated renders should paint the same fraction of
    the canvas -- any gap the per-cell extent fix failed to close would show up as a lower
    non-background fraction specifically at the rotated orientation, not baked into the oval
    shape itself the way an absolute threshold would be."""
    world = _world(seed=8)
    rotation = geometry.rotation_matrix(np.array([0.4, -0.5, 0.7]), 1.3)

    def non_background_fraction(rotation_matrix):
        png = render_image.render_png(world, "eckert4", "temperature", 550, 306, rotation_matrix)
        pixels = np.asarray(Image.open(io.BytesIO(png)).convert("RGB"))
        return np.mean(np.any(pixels != np.array(render_image.BACKGROUND_RGB), axis=-1))

    frac_identity = non_background_fraction(np.eye(3))
    frac_rotated = non_background_fraction(rotation)
    assert abs(frac_identity - frac_rotated) < 0.01


def test_rotation_arc_direction_mirrors_when_omega_sign_flips():
    """Two plates at the same seed, rotating at the same rate but in opposite senses,
    should render different (mirror-image) arcs -- confirms the arc's sweep direction is
    actually sensitive to the sign of omega, not just its magnitude."""
    from app.plates import ElevationLine, PlateWithLines

    def make_plate(plate_id, omega_sign):
        seed_xyz = np.array([1.0, 0.0, 0.0])
        frame = geometry.plate_frame_from_seed(seed_xyz)
        lines = [
            ElevationLine(phi=float(phi), theta=np.linspace(-0.3, 0.3, 8), elevation=np.full(8, 100.0))
            for phi in np.linspace(-0.3, 0.3, 8)
        ]
        omega = omega_sign * 0.03 * seed_xyz  # pole exactly at the seed either way
        return PlateWithLines(plate_id=plate_id, frame=frame, crust_type="continental", omega=omega, lines=lines)

    world_pos = World(seed=1, plates=[make_plate(0, +1.0)])
    world_neg = World(seed=1, plates=[make_plate(0, -1.0)])
    png_pos = render_image.render_png(world_pos, "eckert4", "plates", 400, 400)
    png_neg = render_image.render_png(world_neg, "eckert4", "plates", 400, 400)
    assert png_pos != png_neg


# -- Issue #133 phase 1: node_cloud_resample_mode ---------------------------------------------


def test_node_cloud_resample_mode_default_is_kdtree_and_unaffected_by_the_new_flag_existing():
    """The default ("kdtree") path must render byte-identical output whether or not
    `node_cloud_resample_mode` is ever touched -- the flag is additive, opt-in only."""
    world_untouched = _world(seed=2, num_plates=8)
    world_explicit_default = _world(seed=2, num_plates=8)
    world_explicit_default.node_cloud_resample_mode = "kdtree"
    png_untouched = render_image.render_png(world_untouched, "eckert4", "elevation", 200, 100)
    png_explicit_default = render_image.render_png(world_explicit_default, "eckert4", "elevation", 200, 100)
    assert png_untouched == png_explicit_default


def test_node_cloud_and_tree_healpix_mode_returns_a_node_pixel_index():
    world = _world(seed=2, num_plates=8)
    world.node_cloud_resample_mode = "healpix"
    node_cloud = render_image._node_cloud_and_tree(world)
    assert node_cloud is not None
    all_points, all_elevation, all_owner, tree = node_cloud
    assert isinstance(tree, healpix_grid.NodePixelIndex)
    assert not np.any(tree.pixel_to_node == -1)
    # The 4th slot differs by mode, but the node cloud itself (points/elevation/owner) must
    # not -- both modes resample the exact same underlying data.
    world.node_cloud_resample_mode = "kdtree"
    world.node_kdtree_cache = None
    kd_points, kd_elevation, kd_owner, _kd_tree = render_image._node_cloud_and_tree(world)
    assert np.array_equal(all_points, kd_points)
    assert np.array_equal(all_elevation, kd_elevation)
    assert np.array_equal(all_owner, kd_owner)


def test_render_png_healpix_mode_renders_without_crashing_including_terrain_relief():
    """Exercises every `_node_cloud_and_tree` consumer end-to-end under "healpix" mode via a
    real `render_png` call, including the Elevation view's relief toggles -- the one path
    (`_classify_terrain_relief`'s radius search) that can't go through a `NodePixelIndex` at
    all and must fall back to `_relief_kdtree`'s own real `cKDTree`."""
    world = _world(seed=3, num_plates=8)
    world.node_cloud_resample_mode = "healpix"
    for view in ("elevation", "combined", "biome", "crustType"):
        png = render_image.render_png(world, "eckert4", view, 200, 100)
        assert len(png) > 0
    png_relief = render_image.render_png(
        world, "eckert4", "elevation", 200, 100, show_mountains=True, show_plains_plateaus=True
    )
    assert len(png_relief) > 0
    assert world.node_kdtree_relief_cache is not None


def test_step_world_rebuilds_healpix_index_cache_but_keeps_the_grid_cache():
    world = _world(seed=4, num_plates=8)
    world.node_cloud_resample_mode = "healpix"
    render_image.render_png(world, "eckert4", "elevation", 200, 100)
    assert world.node_healpix_index_cache is not None
    index_before = world.node_healpix_index_cache
    grid_cache_before = world.node_healpix_grid_cache
    assert grid_cache_before is not None

    step_world(world, 1_000_000)
    # A node moved, so the stale index must not survive -- but since issue #133 phase 2,
    # climate._sample_elevation_and_crust shares this same cache and runs every step_world, so
    # by the time step_world returns the cache has already been rebuilt in-step (the exact same
    # "climate rebuilds it before the next render even asks" behavior node_position_tree_cache
    # already had under "kdtree" mode) -- not left None until the next render call the way it
    # was in phase 1, before climate had a healpix branch of its own.
    assert world.node_healpix_index_cache is not None
    assert world.node_healpix_index_cache is not index_before
    assert world.node_healpix_grid_cache is grid_cache_before  # node count unchanged -- reused

    render_image.render_png(world, "eckert4", "elevation", 200, 100)
    assert world.node_healpix_index_cache is not None

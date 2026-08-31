import base64
import io
import numpy as np
from PIL import Image
from app import climate, geometry, hydrology, render_image
from app.world import World, generate_world


def _world(seed=1, num_plates=10, continental_fraction=0.4):
    return generate_world(seed, num_plates=num_plates, continental_fraction=continental_fraction)


def test_render_grid_arrays_cover_the_sphere_with_no_gaps():
    world = _world()
    xy, elevation, plate_id, lake_depth, glacier_depth, is_volcano, half_w, half_h = render_image._render_grid_arrays(world, "behrmann", np.eye(3))

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


def test_combined_view_encodes_biome_ids_in_the_alpha_channel():
    # Combined's per-pixel biome id rides in alpha (see render_image.COMBINED_LAKE_ID_CODE's
    # comment): alpha = 255 - code, code 0 for ocean/unclassified, biome_id + 1 for land.
    world = _world()
    png = render_image.render_png(world, "behrmann", "combined", 320, 180)
    image = Image.open(io.BytesIO(png))
    assert image.mode == "RGBA"

    alpha = np.asarray(image.convert("RGBA"))[:, :, 3]
    codes = 255 - alpha.astype(int)
    # Background/ocean stays code 0; land carries real biome codes within the encodable range.
    assert codes.min() == 0
    assert codes.max() > 0
    assert codes.max() <= render_image.COMBINED_GLACIER_ID_CODE
    assert np.count_nonzero(codes) > 0

    # Other views stay plain RGB -- alpha is a Combined-only channel.
    elev = Image.open(io.BytesIO(render_image.render_png(world, "behrmann", "elevation", 320, 180)))
    assert elev.mode == "RGB"


def test_biome_view_smoothing_preserves_the_major_biomes_and_barely_moves_the_rest():
    # smooth_biome_field is a cleanup pass, not a reclassification: on the real biome render
    # grid it should change only a small slice of land and never erase a biome that has a
    # genuine regional presence.
    from app import biomes

    world = _world(seed=7, num_plates=12, continental_fraction=0.6)
    lat_deg, _lon, _xyz, elevation_m, is_ocean, air_temp, ocean_temp, precip, _lake, _gl = render_image._biome_fields(
        world, *render_image.biome_grid_dimensions(world.climate_density)
    )
    display_temp = np.where(is_ocean, ocean_temp, air_temp)
    slope = biomes.grid_slope(elevation_m, lat_deg)
    raw = biomes.classify_biomes(
        display_temp.reshape(-1), precip.reshape(-1), elevation_m.reshape(-1), slope.reshape(-1), is_ocean.reshape(-1), world.sea_level_m
    ).reshape(elevation_m.shape)
    smoothed = biomes.smooth_biome_field(display_temp, precip, elevation_m, slope, is_ocean, world.sea_level_m)

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


def test_draw_rivers_only_draws_segments_above_the_flow_floor():
    # The main map views (unlike the River Inspector, which deliberately shows every
    # is_river-classified network regardless of flow_rate) additionally require flow_accum
    # to clear river_draw_min_flow(world) -- confirmed here by drawing the exact same
    # is_river/flow_target topology twice, once with flow_accum below the floor (no river
    # pixels should appear at all) and once above it (some must). Calls _draw_rivers directly
    # (rather than the full render_png) against a blank background so any pixel that isn't
    # BACKGROUND_RGB can only be the river line itself (now antialiased -- see
    # RIVER_BLUR_RADIUS_PX -- so drawn pixels blend toward RIVER_COLOR_RGB rather than
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
        water_deposited=np.zeros(n),
        filled_elevation=np.zeros(n),
        spill_target=np.full(n, -1, dtype=np.int64),
        is_river=np.array([True, False, True, False]),
        flow_target=np.array([1, -1, 3, -1]),
        lake_depth=np.zeros(n),
        glacier_depth=np.zeros(n),
        plates_in_order=[],
    )

    def river_pixel_count(flow_accum):
        world.hydrology_cache = hydrology.HydrologyFields(flow_accum=flow_accum, **base_fields)
        image = Image.new("RGB", (300, 200), render_image.BACKGROUND_RGB)
        image = render_image._draw_rivers(image, world, "behrmann", 50.0, 150.0, 100.0, 1.0, np.eye(3))
        pixels = np.asarray(image)
        return int(np.any(pixels != np.array(render_image.BACKGROUND_RGB), axis=-1).sum())

    floor = render_image.river_draw_min_flow(world)
    below_floor = river_pixel_count(np.full(n, floor * 0.5))
    above_floor = river_pixel_count(np.full(n, floor * 2.0))

    assert below_floor == 0
    assert above_floor > 0


def test_river_draw_min_flow_increases_with_generation_resolution():
    # A world generated at a higher elevation-point density and/or climate & biome resolution
    # should require a higher flow_accum to draw a river segment at all, since flow_accum is
    # being swept from a finer-grained network -- see river_draw_min_flow's own docstring.
    world = _world()
    world.node_density = 1.0
    world.climate_density = 1.0
    lowest = render_image.river_draw_min_flow(world)
    assert lowest == render_image.RIVER_DRAW_MIN_FLOW

    world.node_density = 4.0
    higher_node_density = render_image.river_draw_min_flow(world)
    assert higher_node_density > lowest

    world.node_density = 1.0
    world.climate_density = 2.0
    higher_climate_density = render_image.river_draw_min_flow(world)
    assert higher_climate_density > lowest

    world.node_density = 4.0
    world.climate_density = 2.0
    assert render_image.river_draw_min_flow(world) > max(higher_node_density, higher_climate_density)


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
    xy, elevation, plate_id, lake_depth, glacier_depth, is_volcano, half_w, half_h = render_image._render_grid_arrays(world, "eckert4", rotation)
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

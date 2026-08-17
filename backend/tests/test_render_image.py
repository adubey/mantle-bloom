import base64
import io

import numpy as np
from PIL import Image

from app import climate, geometry, render_image
from app.world import World, generate_world


def _world(seed=1, num_plates=10, continental_fraction=0.4):
    return generate_world(seed, num_plates=num_plates, continental_fraction=continental_fraction)


def test_render_grid_arrays_cover_the_sphere_with_no_gaps():
    world = _world()
    xy, elevation, plate_id, half_w, half_h = render_image._render_grid_arrays(world, "behrmann", np.eye(3))

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
    assert tuple(colors[1]) == (200, 210, 150)
    assert tuple(colors[2]) == (255, 255, 255)


def test_elevation_colors_clamps_outside_the_stop_range():
    colors = render_image.elevation_colors(np.array([-999999.0, 999999.0]))
    assert tuple(colors[0]) == (10, 10, 40)
    assert tuple(colors[1]) == (255, 255, 255)


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


def test_every_view_draws_a_legend_panel():
    world = _world(seed=5, num_plates=8, continental_fraction=0.6)
    for view in render_image.VIEWS:
        png = render_image.render_png(world, "behrmann", view, 550, 306)
        image = Image.open(io.BytesIO(png)).convert("RGB")
        pixels = np.asarray(image)
        # The legend panel is a solid-color rectangle anchored at the bottom-left corner --
        # its exact background color should appear somewhere in that corner for every view
        # now that all of them draw one.
        corner = pixels[-90:, :220]
        assert np.any(np.all(corner == np.array(render_image.LEGEND_BG_RGB), axis=-1)), f"{view} has no legend panel"


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


def _synthetic_converging_currents_fields(height=20, width=40) -> "climate.ClimateFields":
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_grid = np.repeat(lon_deg[None, :], height, axis=0)
    world_xyz = geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))

    is_ocean = np.ones((height, width), dtype=bool)
    current_u = np.zeros((height, width))
    current_u[:, : width // 2] = 3.0
    current_u[:, width // 2 :] = -3.0  # converges on the seam at width // 2
    current_v = np.zeros((height, width))

    rows, cols = climate.compute_ocean_swells(current_u, current_v, is_ocean, np.random.default_rng(0))
    zeros = np.zeros((height, width))
    return climate.ClimateFields(
        lat_deg=lat_deg, lon_deg=lon_deg, world_xyz=world_xyz,
        elevation_m=zeros, is_ocean=is_ocean,
        land_temperature_c=zeros, ocean_temperature_c=zeros, air_temperature_c=zeros,
        wind_u=zeros, wind_v=zeros, current_u=current_u, current_v=current_v,
        humidity=zeros, precipitation_mm=zeros, swell_rows=rows, swell_cols=cols,
    )


def test_ocean_currents_view_marks_swells_at_synthetic_convergence(monkeypatch):
    fields = _synthetic_converging_currents_fields()
    assert len(fields.swell_rows) > 0  # sanity: the synthetic setup actually produces swells
    monkeypatch.setattr(render_image.climate, "compute_climate", lambda world, **kwargs: fields)

    png = render_image.render_png(_world(), "behrmann", "oceanCurrents", 320, 180)
    image = Image.open(io.BytesIO(png)).convert("RGB")
    pixels = np.asarray(image).reshape(-1, 3)
    swell_marker_white = np.array([255, 255, 255])
    assert np.any(np.all(pixels == swell_marker_white, axis=-1))


def test_rotate_maps_a_known_point_to_its_expected_position():
    """The core operation the whole view-rotation feature rests on: a 180-degree rotation
    about the z-axis should send lat=0/lon=0 to lat=0/lon=180 (its antipode on the equator)."""
    rotation = geometry.rotation_matrix(np.array([0.0, 0.0, 1.0]), np.pi)
    origin = geometry.latlon_to_xyz(np.array([0.0]), np.array([0.0]))
    rotated = render_image._rotate(origin, rotation)
    lat, lon = geometry.xyz_to_latlon(rotated)
    assert np.isclose(lat[0], 0.0, atol=1e-9)
    assert np.isclose(abs(lon[0]), np.pi, atol=1e-9)


def test_identity_rotation_matches_omitting_rotation_entirely():
    world = _world(seed=6, num_plates=8, continental_fraction=0.5)
    for view in ("elevation", "plates", "temperature", "oceanCurrents"):
        default = render_image.render_png(world, "eckert4", view, 320, 180)
        explicit_identity = render_image.render_png(world, "eckert4", view, 320, 180, np.eye(3))
        assert default == explicit_identity, f"{view}: omitting rotation should match explicit identity"


def test_rotation_visibly_changes_every_view():
    world = _world(seed=7, num_plates=9, continental_fraction=0.5)
    rotation = geometry.rotation_matrix(np.array([0.3, 0.6, 0.1]), 2.1)  # an arbitrary, non-trivial rotation
    for view in ("elevation", "plates", "platesDetail", "temperature", "wind"):
        identity_png = render_image.render_png(world, "eckert4", view, 320, 180)
        rotated_png = render_image.render_png(world, "eckert4", view, 320, 180, rotation)
        assert identity_png != rotated_png, f"{view}: a non-trivial rotation should change the render"


def test_render_grid_stays_gap_free_under_a_nontrivial_rotation():
    """The per-cell half-extent fix (see _render_grid_arrays) exists specifically so rotation
    doesn't reopen the gaps the render grid was built to avoid -- re-run the original
    no-gaps assertions at a rotation that mixes all three axes, not just identity/90/180."""
    world = _world(seed=8)
    rotation = geometry.rotation_matrix(np.array([0.4, -0.5, 0.7]), 1.3)
    xy, elevation, plate_id, half_w, half_h = render_image._render_grid_arrays(world, "eckert4", rotation)
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
    from app.plates import ElevationLine, Plate

    def make_plate(plate_id, omega_sign):
        seed_xyz = np.array([1.0, 0.0, 0.0])
        frame = geometry.plate_frame_from_seed(seed_xyz)
        lines = [
            ElevationLine(phi=float(phi), theta=np.linspace(-0.3, 0.3, 8), elevation=np.full(8, 100.0))
            for phi in np.linspace(-0.3, 0.3, 8)
        ]
        omega = omega_sign * 0.03 * seed_xyz  # pole exactly at the seed either way
        return Plate(plate_id=plate_id, frame=frame, crust_type="continental", omega=omega, lines=lines)

    world_pos = World(seed=1, plates=[make_plate(0, +1.0)])
    world_neg = World(seed=1, plates=[make_plate(0, -1.0)])
    png_pos = render_image.render_png(world_pos, "eckert4", "plates", 400, 400)
    png_neg = render_image.render_png(world_neg, "eckert4", "plates", 400, 400)
    assert png_pos != png_neg

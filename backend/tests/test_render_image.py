import base64
import io

import numpy as np
from PIL import Image

from app import geometry, render_image
from app.world import World, generate_world


def _world(seed=1, num_plates=10, continental_fraction=0.4):
    return generate_world(seed, num_plates=num_plates, continental_fraction=continental_fraction)


def test_render_grid_arrays_cover_the_sphere_with_no_gaps():
    world = _world()
    xy, elevation, plate_id, half_w, half_h = render_image._render_grid_arrays(world, "behrmann")

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
    behrmann_xy, *_ = render_image._render_grid_arrays(world, "behrmann")
    eckert4_xy, *_ = render_image._render_grid_arrays(world, "eckert4")
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
        info = render_image._plate_tectonics("eckert4", plate)
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

import base64
import io
import numpy as np
from PIL import Image, ImageDraw
from app import climate, geometry, hydrology, render_image
from app.world import World, generate_world


def _world(seed=1, num_plates=10, continental_fraction=0.4):
    return generate_world(seed, num_plates=num_plates, continental_fraction=continental_fraction)


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

import base64
import io
from types import SimpleNamespace

import numpy as np
import pytest
from PIL import Image

from app import cross_section
from app.elevation_lines import CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC
from .quad_fixtures import node_points, row_plate


@pytest.fixture
def section_world():
    theta = np.linspace(-0.3, 0.3, 61)
    elevation = np.linspace(-8000.0, 3000.0, len(theta))
    plate = row_plate(
        1,
        "continental",
        theta,
        elevation=elevation,
        crustal_thickness_m=np.linspace(20_000.0, 50_000.0, len(theta)),
        mantle_lithosphere_thickness_m=np.linspace(80_000.0, 160_000.0, len(theta)),
        crust_type_code=np.where(np.arange(len(theta)) < len(theta) // 2, CRUST_TYPE_OCEANIC, CRUST_TYPE_CONTINENTAL),
        glacier_depth=np.where(elevation > 0.0, 150.0, 0.0),
        is_volcano=np.arange(len(theta)) % 12 == 0,
        volcano_active_years_remaining=np.where(np.arange(len(theta)) % 24 == 0, 1_000_000.0, 0.0),
        moho_thermal_lag_c=np.zeros(len(theta)),
    )
    return SimpleNamespace(
        plates=[plate],
        node_density=4.0,
        sea_level_m=0.0,
        solar_multiplier=1.0,
        elapsed_years=0.0,
        steps_taken=0,
        mantle_centers=[],
        pending_magma_parcels=[SimpleNamespace(origin_xyz=node_points(plate)[30], volume_m3=1e12)],
        node_position_tree_cache=None,
    )


def test_great_circle_samples_follow_rotated_plane_and_fixed_arc_length():
    points, tangent, distance_km = cross_section._great_circle_samples(90.0, 30, 60.0)

    assert np.allclose(np.linalg.norm(points, axis=1), 1.0)
    assert np.allclose(np.sum(points * tangent, axis=1), 0.0, atol=1e-12)
    assert distance_km[0] == pytest.approx(-cross_section.ARC_LENGTH_KM / 2.0)
    assert distance_km[-1] == pytest.approx(cross_section.ARC_LENGTH_KM / 2.0)
    assert points[len(points) // 2] == pytest.approx([0.0, 1.0, 0.0])


def test_section_sampling_keeps_categorical_crust_values_discrete(section_world):
    plate = section_world.plates[0]
    points = node_points(plate)[[0, -1]]
    sampled, continental, _ = cross_section._sample_fields(section_world, points)

    assert np.array_equal(sampled["crust_type_code"], plate.collect("crust_type_code")[[0, -1]])
    assert continental.tolist() == [False, True]
    assert np.array_equal(sampled["elevation"], plate.collect("elevation")[[0, -1]])


def test_render_returns_png_and_section_revision(section_world, monkeypatch):
    monkeypatch.setattr(
        cross_section.climate,
        "compute_climate_cached",
        lambda _world: SimpleNamespace(biome_ids=np.zeros((90, 180), dtype=int)),
    )

    image_base64, metadata = cross_section.render_png_base64(section_world, 0.0, 0, 320, 180)

    image = Image.open(io.BytesIO(base64.b64decode(image_base64)))
    assert image.format == "PNG"
    assert image.size == (320, 180)
    assert len(image.getcolors(maxcolors=320 * 180)) > 10
    assert metadata["image_base64"] == image_base64
    assert metadata["world_revision"].startswith(f"{id(section_world):x}:")
    assert metadata["arc_length_km"] == cross_section.ARC_LENGTH_KM
    assert metadata["section_depth_km"] == cross_section.SECTION_DEPTH_KM


def test_render_keeps_sampled_crust_and_lithosphere_above_schematic_asthenosphere(monkeypatch):
    sample_xyz = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
    tangent = np.zeros_like(sample_xyz)
    distances = np.array([-1500.0, 0.0, 1500.0])
    fields = {
        "elevation": np.zeros(3),
        "crustal_thickness_m": np.full(3, 35_000.0),
        "mantle_lithosphere_thickness_m": np.full(3, 120_000.0),
        "crust_type_code": np.full(3, CRUST_TYPE_OCEANIC),
        "glacier_depth": np.zeros(3),
        "is_volcano": np.zeros(3, dtype=bool),
        "volcano_active_years_remaining": np.zeros(3),
        "moho_thermal_lag_c": np.zeros(3),
    }
    monkeypatch.setattr(cross_section, "_great_circle_samples", lambda *_: (sample_xyz, tangent, distances))
    monkeypatch.setattr(cross_section, "_sample_fields", lambda *_: (fields, np.zeros(3, dtype=bool), sample_xyz))
    monkeypatch.setattr(cross_section, "_biome_colors", lambda *_: np.zeros((3, 3)))
    monkeypatch.setattr(cross_section.orogeny, "moho_temperature_c", lambda hc, hm, lag: np.zeros_like(hc))
    monkeypatch.setattr(cross_section, "_geotherm_temperature", lambda hc, hm, lag, depth: np.zeros_like(depth))
    monkeypatch.setattr(cross_section.mantle, "flow_at", lambda points, _centers: np.zeros_like(points))
    world = SimpleNamespace(
        node_density=4.0,
        sea_level_m=1000.0,
        elapsed_years=0.0,
        solar_multiplier=1.0,
        steps_taken=0,
        mantle_centers=[],
        pending_magma_parcels=[],
    )

    image_base64, _ = cross_section.render_png_base64(world, 0.0, 0, 800, 500)
    image = Image.open(io.BytesIO(base64.b64decode(image_base64))).convert("RGB")
    margin_top, margin_bottom = 38, 55
    vertical_scale = (500 - margin_top - margin_bottom) / (
        cross_section.SECTION_DEPTH_KM + cross_section.RELIEF_ALLOWANCE_KM + 10.0
    )
    lower_limit_y = 500 - margin_bottom - cross_section.RELIEF_ALLOWANCE_KM * 0.5 * vertical_scale
    surface_y = round(lower_limit_y - cross_section.SECTION_DEPTH_KM * vertical_scale)
    sample_x = 800 // 2 - 15  # Clear of the center depth-tick labels.

    assert image.getpixel((sample_x, surface_y + 10)) == cross_section.CRUST_OCEANIC
    assert image.getpixel((sample_x, surface_y + round(70 * vertical_scale))) == cross_section.MANTLE_LITHOSPHERE
    assert image.getpixel((sample_x, surface_y + round(250 * vertical_scale))) == cross_section.ASTHENOSPHERE


def test_sampling_empty_world_reports_no_active_surface_nodes():
    world = SimpleNamespace(plates=[], node_position_tree_cache=None)

    with pytest.raises(ValueError, match="no active surface nodes"):
        cross_section._sample_fields(world, np.array([[1.0, 0.0, 0.0]]))

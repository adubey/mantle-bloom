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
    assert np.array_equal(sampled["nearest_node"], [0, 60])


def test_section_sampling_blends_continuous_fields_along_the_arc(section_world):
    """#702: copying the nearest node drew each interface as a node-wide staircase."""
    plate = section_world.plates[0]
    nodes = node_points(plate)
    t = np.linspace(0.0, 1.0, 101)[:, None]
    points = nodes[29] * (1.0 - t) + nodes[31] * t
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    sampled, _, _ = cross_section._sample_fields(section_world, points)

    node_step = float(np.diff(plate.collect("elevation"))[0])
    steps = np.diff(sampled["elevation"])
    assert np.all(steps > 0.0)
    assert steps.max() < 0.1 * node_step


def test_sampling_handles_worlds_with_fewer_nodes_than_the_kernel(monkeypatch):
    nodes = np.array([[1.0, -0.005, 0.0], [1.0, 0.005, 0.0]])
    nodes /= np.linalg.norm(nodes, axis=1, keepdims=True)
    fields = {name: np.zeros(2) for name in cross_section.CONTINUOUS_FIELDS}
    fields["elevation"] = np.array([0.0, 10.0])
    monkeypatch.setattr(cross_section, "_world_surface", lambda _world: (nodes, fields, np.zeros(2, dtype=bool)))
    world = SimpleNamespace(node_density=4.0, node_position_tree_cache=None)

    sampled, _, _ = cross_section._sample_fields(world, np.array([[1.0, 0.0, 0.0]]))

    assert sampled["elevation"][0] == pytest.approx(5.0)


def test_biome_colors_interpolate_display_colors_between_cells(monkeypatch):
    """#702: blend colours, not ids -- the mean of two ids is an unrelated biome."""
    ids = np.zeros((2, 4), dtype=int)
    ids[:, 2:] = 1
    monkeypatch.setattr(
        cross_section.climate, "compute_climate_cached", lambda _world: SimpleNamespace(biome_ids=ids)
    )
    # Longitude 0 is the shared edge of columns 1 and 2, midway between their centres.
    colors = cross_section._biome_colors(None, np.array([[1.0, 0.0, 0.0]]), np.zeros(1))
    palette = cross_section.biomes.BIOME_COLORS.astype(float)
    assert colors[0] == pytest.approx((palette[0] + palette[1]) / 2.0)


def test_vent_samples_give_each_volcano_node_one_middle_sample():
    """#702: the fine arc gives a node many samples, and neighbouring vents stay separate."""
    node_ids = np.array([5, 5, 5, 8, 8, 9, 9, 9, 9, 2])
    volcano = np.array([False] * 3 + [True] * 2 + [True] * 4 + [False])
    assert cross_section._vent_samples(volcano, node_ids).tolist() == [3, 6]


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
    tangent = np.tile([0.0, 1.0, 0.0], (len(sample_xyz), 1))
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
    monkeypatch.setattr(
        cross_section.mantle,
        "flow_at",
        lambda points, _centers: np.tile([0.0, 1.0, 0.0], (len(points), 1)),
    )
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
    margin_x, margin_top = 20, 38
    arc_width = (800 - 2 * margin_x) * (690.0 / 820.0)
    horizontal_scale = arc_width / (
        2.0 * cross_section.PLANET_RADIUS_KM
        * np.sin(cross_section.ARC_LENGTH_KM / (2.0 * cross_section.PLANET_RADIUS_KM))
    )
    vertical_scale = cross_section.VERTICAL_EXAGGERATION * horizontal_scale
    surface_y = margin_top + round(500 * 0.08)
    sample_x = 800 // 2

    assert image.getpixel((sample_x, surface_y + 10)) == cross_section.CRUST_OCEANIC
    assert image.getpixel((sample_x, surface_y + round(70 * vertical_scale))) == cross_section.MANTLE_LITHOSPHERE
    assert image.getpixel((sample_x + 40, surface_y + round(300 * vertical_scale))) == cross_section.ASTHENOSPHERE
    flow_depth_km = (35_000.0 + 120_000.0) / 1000.0
    flow_depth_km += (cross_section.SECTION_DEPTH_KM - flow_depth_km) * 0.5
    flow_y = surface_y + round(flow_depth_km * vertical_scale)
    assert image.getpixel((sample_x + 1, flow_y)) == (71, 236, 210)
    assert image.getpixel((sample_x, surface_y - 2)) == (14, 20, 34)
    arc_half_angle = cross_section.ARC_LENGTH_KM / (2.0 * cross_section.PLANET_RADIUS_KM)
    edge_x = int(np.ceil(800 / 2 + arc_width / 2))
    edge_surface_y = surface_y + round(
        cross_section.VERTICAL_EXAGGERATION
        * cross_section.PLANET_RADIUS_KM
        * (1.0 - np.cos(arc_half_angle))
        * horizontal_scale
    )
    # The arc's last column is one drawn column wide, so the output filter blends it with
    # the background: it must still carry crust, just not at full coverage.
    edge = np.array(image.getpixel((edge_x, edge_surface_y + 2)), dtype=float)
    background = np.array((14.0, 20.0, 34.0))
    coverage = (edge - background) / (np.array(cross_section.CRUST_OCEANIC) - background)
    assert coverage == pytest.approx(np.full(3, coverage[0]), abs=0.05)
    assert coverage[0] > 0.3
    assert image.getpixel((40, surface_y + 20)) == (14, 20, 34)


def test_sampling_empty_world_reports_no_active_surface_nodes():
    world = SimpleNamespace(plates=[], node_position_tree_cache=None)

    with pytest.raises(ValueError, match="no active surface nodes"):
        cross_section._sample_fields(world, np.array([[1.0, 0.0, 0.0]]))


def _marker_world(plumes):
    return SimpleNamespace(mantle_plumes=plumes, elapsed_years=1_000_000.0)


def _live_head(position, head_end_years=3_000_000.0, tail_end_years=80_000_000.0):
    from app import plumes as plumes_mod

    return plumes_mod.MantlePlume(
        plume_id=0,
        position=tuple(position),
        born_years=0.0,
        head_volume_m3=8e15,
        head_remaining_m3=8e15,
        head_end_years=head_end_years,
        tail_flux_m3_per_yr=1e8,
        tail_end_years=tail_end_years,
    )


def _tail_only(position):
    import dataclasses

    return dataclasses.replace(
        _live_head(position),
        plume_id=-1,
        head_volume_m3=0.0,
        head_remaining_m3=0.0,
        head_end_years=0.0,
    )


def test_plume_markers_flag_head_footprint_and_tail_roles():
    """RFC 0007 (#669): a live head marks its nearest sample plus the
    intersected footprint extent; a live tail-only plume marks tail."""
    from app import plumes as plumes_mod

    sample_xyz, _, _ = cross_section._great_circle_samples(0.0, 0, 60.0)
    middle = len(sample_xyz) // 2
    world = _marker_world([_live_head(sample_xyz[middle]), _tail_only(sample_xyz[0])])

    markers = dict(cross_section._plume_markers(world, sample_xyz, 60.0))

    assert markers[middle] == "head"
    radius_km = (
        plumes_mod.head_kernel_radius_m(plumes_mod.head_footprint_area_m2(8e15)) / 1000.0
    )
    # Along-arc step between adjacent samples, from the sample geometry itself.
    step_km = float(
        np.arccos(np.clip(sample_xyz[0] @ sample_xyz[1], -1.0, 1.0))
        * cross_section.PLANET_RADIUS_KM
    )
    expected_footprint = {
        index
        for index in range(len(sample_xyz))
        if index != middle and abs(index - middle) * step_km <= radius_km
    }
    assert radius_km > 60.0  # the test only means something with real extent
    footprint = {index for index, role in markers.items() if role == "head_footprint"}
    assert footprint == expected_footprint
    tail_samples = [index for index, role in markers.items() if role == "tail"]
    assert tail_samples == [0]  # the tail-only plume's nearest sample


def test_render_draws_plume_conduit_and_glyphs_once(monkeypatch):
    """RFC 0007 (#669): the render loop visits every pixel column, so surface
    glyphs must paint once -- the head ellipse keeps its fill instead of being
    buried by repeated outlines. #702: the tail is a conduit of its schematic
    width, open at the section floor rather than a hairline."""
    sample_xyz = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
    tangent = np.tile([0.0, 1.0, 0.0], (len(sample_xyz), 1))
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
    monkeypatch.setattr(
        cross_section.mantle,
        "flow_at",
        lambda points, _centers: np.zeros((len(points), 3)),
    )
    world = SimpleNamespace(
        node_density=4.0,
        sea_level_m=1000.0,
        elapsed_years=0.0,
        solar_multiplier=1.0,
        steps_taken=0,
        mantle_centers=[],
        pending_magma_parcels=[],
        mantle_plumes=[_live_head(sample_xyz[1]), _tail_only(sample_xyz[0])],
    )

    image_base64, _ = cross_section.render_png_base64(world, 0.0, 0, 800, 500, show_plumes=True)
    image = Image.open(io.BytesIO(base64.b64decode(image_base64))).convert("RGB")
    pixels = np.asarray(image)
    conduit = np.all(pixels == cross_section.PLUME_TAIL_FILL, axis=-1)
    arc_width = (800 - 2 * 20) * (690.0 / 820.0)
    horizontal_scale = arc_width / (
        2.0 * cross_section.PLANET_RADIUS_KM
        * np.sin(cross_section.ARC_LENGTH_KM / (2.0 * cross_section.PLANET_RADIUS_KM))
    )
    # The tail sits on the arc's first sample, so the section edge clips its left half.
    half_width_px = cross_section.PLUME_TAIL_CONDUIT_RADIUS_KM * horizontal_scale
    assert conduit.sum(axis=1).max() == pytest.approx(half_width_px, abs=3.0)
    # Flush with the floor: the conduit's last row sits on the section background.
    column = int(np.argmax(conduit.sum(axis=0)))
    floor_row = int(np.flatnonzero(conduit[:, column])[-1])
    assert tuple(pixels[floor_row + 2, column]) == (14, 20, 34)
    assert not np.any(np.all(pixels[:, column] == cross_section.ASTHENOSPHERE, axis=-1)[floor_row:])
    # One head ellipse plus one tail ellipse worth of dark outline: the
    # per-column repaint drew one outline per owned column (~250 each here).
    # Found by colour so the pin survives layout changes.
    outlines = int(np.all(pixels == (45, 24, 24), axis=-1).sum())
    assert outlines <= 80


def _sized_ocean_fields(count):
    return {
        "elevation": np.zeros(count),
        "crustal_thickness_m": np.full(count, 35_000.0),
        "mantle_lithosphere_thickness_m": np.full(count, 120_000.0),
        "crust_type_code": np.full(count, CRUST_TYPE_OCEANIC),
        "glacier_depth": np.zeros(count),
        "is_volcano": np.zeros(count, dtype=bool),
        "volcano_active_years_remaining": np.zeros(count),
        "moho_thermal_lag_c": np.zeros(count),
    }


def test_render_draws_conduit_for_every_tail_marker_at_fractional_widths(monkeypatch):
    """RFC 0007 (#669): real sections have fractional sample x positions, and
    a marker must survive them -- over several render widths every tail marker
    gets its conduit (round() used to land on the neighbour's column for ~half
    the samples, dropping their glyphs)."""
    # Real geometry (fractional x_samples); everything else schematic.
    monkeypatch.setattr(
        cross_section, "_sample_fields",
        lambda _world, points: (
            _sized_ocean_fields(len(points)), np.zeros(len(points), dtype=bool), points,
        ),
    )
    monkeypatch.setattr(
        cross_section, "_biome_colors", lambda _w, pts, _e: np.zeros((len(pts), 3))
    )
    monkeypatch.setattr(cross_section.orogeny, "moho_temperature_c", lambda hc, hm, lag: np.zeros_like(hc))
    monkeypatch.setattr(cross_section, "_geotherm_temperature", lambda hc, hm, lag, depth: np.zeros_like(depth))
    monkeypatch.setattr(
        cross_section.mantle, "flow_at", lambda points, _c: np.zeros((len(points), 3))
    )
    sample_xyz, _, _ = cross_section._great_circle_samples(0.0, 0, 60.0)
    world = SimpleNamespace(
        node_density=4.0,
        sea_level_m=1000.0,
        elapsed_years=0.0,
        solar_multiplier=1.0,
        steps_taken=0,
        mantle_centers=[],
        pending_magma_parcels=[],
        mantle_plumes=[_tail_only(sample_xyz[i]) for i in (10, 40, 70)],
    )

    for width in (700, 1100, 1600):
        image_base64, _ = cross_section.render_png_base64(world, 0.0, 0, width, 500, show_plumes=True)
        image = Image.open(io.BytesIO(base64.b64decode(image_base64))).convert("RGB")
        pixels = np.asarray(image)
        conduit = np.all(pixels == cross_section.PLUME_TAIL_FILL, axis=-1)
        widest_row = conduit[int(np.argmax(conduit.sum(axis=1)))].astype(int)
        assert int(np.count_nonzero(np.diff(widest_row) == 1)) == 3


def test_plume_markers_omit_plumes_outside_corridor_and_dead_stages():
    """RFC 0007 (#669): a plume off the section plane is omitted, as is a plume
    whose head and tail stages both ended."""
    sample_xyz, _, _ = cross_section._great_circle_samples(0.0, 0, 60.0)
    middle = len(sample_xyz) // 2
    off_plane = _live_head(-sample_xyz[middle])
    assert dict(cross_section._plume_markers(_marker_world([off_plane]), sample_xyz, 60.0)) == {}
    spent = _live_head(sample_xyz[middle], head_end_years=100.0, tail_end_years=100.0)
    assert dict(cross_section._plume_markers(_marker_world([spent]), sample_xyz, 60.0)) == {}


def test_plume_markers_skip_nonfinite_positions():
    """RFC 0007 (#669): a plume with a NaN position must not mark sample 0 --
    argmax over NaNs returns 0 and a NaN distance passes a plain range check."""
    sample_xyz, _, _ = cross_section._great_circle_samples(0.0, 0, 60.0)
    bad = _live_head((float("nan"), float("nan"), float("nan")))
    assert dict(cross_section._plume_markers(_marker_world([bad]), sample_xyz, 60.0)) == {}


def test_plume_markers_draw_cut_footprint_without_ellipse_outside_corridor():
    """RFC 0007 (#669): a head whose centre is just outside the corridor still
    paints the footprint band where the section cuts it -- only the head
    ellipse requires the centre within the corridor. A tail outside the
    corridor stays fully omitted."""
    from app import plumes as plumes_mod

    sample_xyz, _, _ = cross_section._great_circle_samples(0.0, 0, 60.0)
    middle = len(sample_xyz) // 2
    # 200 km off the section plane: outside the 90 km corridor, well inside
    # the ~438 km footprint radius. (Plane 0 samples lie in the x-y plane.)
    off = sample_xyz[middle] + 0.0314 * np.array([0.0, 0.0, 1.0])
    off = tuple(off / np.linalg.norm(off))
    markers = dict(cross_section._plume_markers(_marker_world([_live_head(off)]), sample_xyz, 60.0))
    assert "head" not in markers.values()
    radius_km = (
        plumes_mod.head_kernel_radius_m(plumes_mod.head_footprint_area_m2(8e15)) / 1000.0
    )
    # Measured from the off-plane plume position, so the band is narrower than
    # the centred 29-sample extent -- it shrinks smoothly instead of vanishing.
    off_hat = np.asarray(off)
    expected = {
        index
        for index in range(len(sample_xyz))
        if np.arccos(np.clip(sample_xyz[index] @ off_hat, -1.0, 1.0))
        * cross_section.PLANET_RADIUS_KM
        <= radius_km
    }
    assert 0 < len(expected) < 29
    assert {index for index, role in markers.items() if role == "head_footprint"} == expected
    assert (
        dict(cross_section._plume_markers(_marker_world([_tail_only(off)]), sample_xyz, 60.0))
        == {}
    )


def test_merge_plume_markers_prefers_head_then_tail_over_footprint():
    """RFC 0007 (#669): at a shared sample the head glyph wins over another
    plume's footprint, and a tail wins over a footprint -- the render merge
    must not depend on plume order."""
    sample_xyz, _, _ = cross_section._great_circle_samples(0.0, 0, 60.0)
    middle = len(sample_xyz) // 2
    # B's footprint (radius ~437 km, samples ~30 km apart) covers the middle
    # and middle+5, but B sorts after the tail so a naive last-wins dict would
    # bury both the head and the tail under footprint paint.
    world = _marker_world(
        [
            _live_head(sample_xyz[middle]),
            _tail_only(sample_xyz[middle + 5]),
            _live_head(sample_xyz[middle + 10]),
        ]
    )
    merged = cross_section._merge_plume_markers(world, sample_xyz, 60.0)
    assert merged[middle] == "head"
    assert merged[middle + 5] == "tail"
    assert merged[middle + 10] == "head"


def test_parcel_markers_keep_only_parcels_inside_the_corridor():
    """#702: the bounded tree lookup must match the corridor rule and still scale
    strength by the largest parcel anywhere, including ones off the section."""
    sample_xyz, _, _ = cross_section._great_circle_samples(0.0, 0, 60.0)
    middle = len(sample_xyz) // 2
    off_plane = np.array([0.0, 0.0, 1.0])
    near = sample_xyz[middle] + (80.0 / cross_section.PLANET_RADIUS_KM) * off_plane
    far = sample_xyz[middle] + (120.0 / cross_section.PLANET_RADIUS_KM) * off_plane
    world = SimpleNamespace(pending_magma_parcels=[
        SimpleNamespace(origin_xyz=near / np.linalg.norm(near), volume_m3=1.0e9),
        SimpleNamespace(origin_xyz=far / np.linalg.norm(far), volume_m3=4.0e9),
    ])

    assert cross_section._parcel_markers(world, sample_xyz, 60.0) == [(middle, 0.25)]

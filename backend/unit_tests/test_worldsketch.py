import io

import numpy as np
import pytest
from app import worldsketch
from app.lithosphere_plate import build_plate_tiling, generate_plates
from app.worldsketch import (
    SKETCH_GRID_H,
    SKETCH_GRID_W,
    SketchMasks,
    parse_sketch_image,
    sketch_plate_sites,
)
from PIL import Image, ImageDraw


def _png_bytes(img: Image.Image) -> bytes:
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def _blank_canvas() -> Image.Image:
    return Image.new("RGB", (SKETCH_GRID_W, SKETCH_GRID_H), "white")


def _island_canvas(gap_px: int = 0) -> Image.Image:
    """A roughly-circular coastline outline in the middle of the map, drawn at working-grid
    resolution so pixel coordinates map 1:1 onto the resulting mask -- optionally with a small
    gap in the outline (to exercise gap-closing) of `gap_px` width."""
    img = _blank_canvas()
    draw = ImageDraw.Draw(img)
    cx, cy, r = SKETCH_GRID_W // 2, SKETCH_GRID_H // 2, 60
    draw.ellipse((cx - r, cy - r, cx + r, cy + r), outline=worldsketch.COAST_COLOR, width=3)
    if gap_px > 0:
        # Erase a small arc at the top of the circle by painting over it with white.
        draw.rectangle((cx - gap_px // 2, cy - r - 3, cx + gap_px // 2, cy - r + 3), fill="white")
    return img


def _measured_land_fraction(plates_list) -> float:
    total = sum(p.node_count() for p in plates_list)
    land = sum(int(np.sum(line.elevation > 0)) for p in plates_list for line in p.lines)
    return land / total if total else 0.0


def test_parse_sketch_image_flood_fills_closed_coastline():
    masks = parse_sketch_image(_png_bytes(_island_canvas()))
    cx, cy = SKETCH_GRID_W // 2, SKETCH_GRID_H // 2
    assert masks.land[cy, cx]  # island interior
    assert not masks.land[5, 5]  # a corner, presumed ocean
    assert not masks.land[cy, cx + 200]  # well outside the circle


def test_parse_sketch_image_closes_small_coastline_gaps():
    # A gap a couple of pixels wide should still resolve as a closed island (the tool's
    # documented tolerance for a rough, not-quite-closed sketch) -- see
    # worldsketch._COAST_GAP_CLOSE_ITERS.
    masks = parse_sketch_image(_png_bytes(_island_canvas(gap_px=2)))
    cx, cy = SKETCH_GRID_W // 2, SKETCH_GRID_H // 2
    assert masks.land[cy, cx]


def test_parse_sketch_image_all_white_is_all_ocean():
    masks = parse_sketch_image(_png_bytes(_blank_canvas()))
    assert not masks.land.any()


def test_parse_sketch_image_invalid_bytes_raises():
    with pytest.raises(ValueError):
        parse_sketch_image(b"not an image")


def test_parse_sketch_image_river_and_mountain_restricted_to_land():
    img = _island_canvas()
    draw = ImageDraw.Draw(img)
    cx, cy = SKETCH_GRID_W // 2, SKETCH_GRID_H // 2
    # A river stroke through the island's interior...
    draw.line((cx - 20, cy - 20, cx + 20, cy + 20), fill=worldsketch.RIVER_COLOR, width=3)
    # ...a mountain stroke also inside...
    draw.line((cx - 10, cy + 10, cx + 10, cy - 10), fill=worldsketch.MOUNTAIN_COLOR, width=3)
    # ...and a stray river stroke out in the open ocean, which should be dropped.
    draw.line((10, 10, 30, 10), fill=worldsketch.RIVER_COLOR, width=3)
    masks = parse_sketch_image(_png_bytes(img))
    # The two strokes cross exactly at (cy, cx) -- check a point on each stroke away from that
    # intersection, where only one ink color is actually present.
    assert masks.river[cy - 15, cx - 15]
    assert masks.mountain[cy - 5, cx + 5]
    assert not masks.river[10, 20]  # the ocean stroke never survives the land-only gate
    assert not (masks.river & ~masks.land).any()
    assert not (masks.mountain & ~masks.land).any()


def test_sample_land_matches_grid_at_cell_centers():
    from app import geometry

    masks = parse_sketch_image(_png_bytes(_island_canvas()))
    xyz = masks.cell_centers_xyz()
    sampled = masks.sample_land(xyz.reshape(-1, 3))
    assert np.array_equal(sampled.reshape(masks.land.shape), masks.land)


def _synthetic_two_landmass_masks() -> SketchMasks:
    """A hand-built SketchMasks (bypassing image parsing) with one big landmass spanning most
    of the northern hemisphere and one small island near the south pole -- for exercising
    sketch_plate_sites' area-proportional plate-count split without depending on flood-fill
    exactness."""
    land = np.zeros((SKETCH_GRID_H, SKETCH_GRID_W), dtype=bool)
    land[0:120, :] = True  # a broad northern continent
    land[340:355, 340:360] = True  # a small southern island
    empty = np.zeros_like(land)
    return SketchMasks(land=land, mountain=empty, river=empty)


def test_sketch_plate_sites_splits_large_landmass_and_keeps_small_one():
    rng = np.random.default_rng(0)
    masks = _synthetic_two_landmass_masks()
    site_xyz, crust_types = sketch_plate_sites(masks, num_plates=12, num_continents=4, rng=rng)
    assert len(site_xyz) == len(crust_types)
    n_continental = sum(1 for c in crust_types if c == "continental")
    n_oceanic = sum(1 for c in crust_types if c == "oceanic")
    # The big continent should have soaked up more than one of the 4 continental slots, and
    # the small island still got at least its guaranteed one.
    assert n_continental >= 2
    assert n_oceanic == 12 - n_continental
    # Every site is a unit vector.
    norms = np.linalg.norm(site_xyz, axis=-1)
    assert np.allclose(norms, 1.0, atol=1e-6)


def test_sketch_plate_sites_with_no_land_is_all_oceanic():
    rng = np.random.default_rng(1)
    masks = SketchMasks(
        land=np.zeros((SKETCH_GRID_H, SKETCH_GRID_W), dtype=bool),
        mountain=np.zeros((SKETCH_GRID_H, SKETCH_GRID_W), dtype=bool),
        river=np.zeros((SKETCH_GRID_H, SKETCH_GRID_W), dtype=bool),
    )
    site_xyz, crust_types = sketch_plate_sites(masks, num_plates=10, num_continents=3, rng=rng)
    assert len(site_xyz) == 10
    assert crust_types == ["oceanic"] * 10


def test_build_plate_tiling_primary_sites_are_used_and_normalized():
    rng = np.random.default_rng(2)
    raw = np.array([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]])
    tiling = build_plate_tiling(rng, num_plates=3, extra_sites_per_plate=0, primary_sites=raw)
    assert np.allclose(tiling.site_xyz, np.eye(3))


def test_build_plate_tiling_primary_sites_wrong_length_raises():
    rng = np.random.default_rng(3)
    with pytest.raises(ValueError):
        build_plate_tiling(rng, num_plates=3, primary_sites=np.eye(2))


def test_build_plate_tiling_none_primary_sites_unchanged():
    # Regression: the default (primary_sites=None) path must draw from `rng` exactly as it did
    # before this parameter existed -- a single combined rng.normal(num_plates + num_extra)
    # call, not two separate draws.
    seed = 17
    a = build_plate_tiling(np.random.default_rng(seed), num_plates=6, extra_sites_per_plate=2)
    expected_rng = np.random.default_rng(seed)
    expected = expected_rng.normal(size=(6 + 12, 3))
    expected /= np.linalg.norm(expected, axis=-1, keepdims=True)
    assert np.array_equal(a.site_xyz, expected)


def test_generate_plates_with_sketch_matches_drawn_land_roughly():
    masks = parse_sketch_image(_png_bytes(_island_canvas()))
    plates = generate_plates(seed=6, num_plates=10, continental_fraction=0.5, sketch=masks)
    # A drawn circle covering a modest slice of the map -- just confirm *some* land was
    # actually produced (not an all-ocean or all-land world) and it isn't wildly off from the
    # drawn fraction, without pinning down the terrain-noise texture's exact contribution.
    drawn_fraction = float(masks.land.mean())
    measured = _measured_land_fraction(plates)
    assert 0.0 < measured < 1.0
    assert abs(measured - drawn_fraction) < 0.25


def test_generate_plates_sketch_none_is_unaffected():
    # sketch=None must reproduce the exact same world generate_plates always has -- this is
    # mostly load-bearing on build_plate_tiling's own regression test above, but confirms it
    # end to end too.
    a = generate_plates(seed=41, num_plates=8, continental_fraction=0.6, land_fraction=0.3)
    b = generate_plates(seed=41, num_plates=8, continental_fraction=0.6, land_fraction=0.3, sketch=None)
    assert len(a) == len(b)
    for pa, pb in zip(a, b):
        assert pa.crust_type == pb.crust_type
        assert np.allclose(pa.frame, pb.frame)
        assert len(pa.lines) == len(pb.lines)

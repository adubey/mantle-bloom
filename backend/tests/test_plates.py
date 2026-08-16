import numpy as np

from app import geometry
from app.plates import (
    MAX_AUTO_PLATES,
    MIN_AUTO_PLATES,
    MIN_OCEANIC_PLATES,
    TARGET_LINE_SPACING_RAD,
    generate_plates,
    iter_local_lattice,
)
from app.world import generate_world


def test_iter_local_lattice_default_spacing_matches_target():
    rows = list(iter_local_lattice(np.eye(3)))
    assert len(rows) > 0
    phis = [phi for phi, _, _ in rows]
    assert np.allclose(np.diff(sorted(phis)), TARGET_LINE_SPACING_RAD, atol=1e-9)


def test_iter_local_lattice_custom_spacing_changes_row_count():
    coarse = list(iter_local_lattice(np.eye(3), spacing_rad=TARGET_LINE_SPACING_RAD * 4))
    fine = list(iter_local_lattice(np.eye(3), spacing_rad=TARGET_LINE_SPACING_RAD / 2))
    assert len(fine) > len(coarse)
    total_fine = sum(len(theta) for _, theta, _ in fine)
    total_coarse = sum(len(theta) for _, theta, _ in coarse)
    assert total_fine > total_coarse


def test_generate_plates_count_and_crust_types():
    plates = generate_plates(seed=42, num_plates=10)
    assert len(plates) == 10
    assert all(p.crust_type in ("continental", "oceanic") for p in plates)


def test_every_plate_has_elevation_lines():
    plates = generate_plates(seed=1, num_plates=8)
    for p in plates:
        assert p.node_count() > 0, f"plate {p.plate_id} has no elevation nodes"


def test_frames_are_proper_rotations():
    plates = generate_plates(seed=2, num_plates=6)
    for p in plates:
        assert np.allclose(p.frame @ p.frame.T, np.eye(3), atol=1e-9)
        assert np.isclose(np.linalg.det(p.frame), 1.0)


def test_every_node_is_closest_to_its_own_plate_seed():
    """A node kept by a plate must actually be in that plate's Voronoi cell -- i.e. no
    other plate's seed is angularly closer to it."""
    plates = generate_plates(seed=3, num_plates=8)
    seeds = np.array([p.seed_world for p in plates])

    for p in plates:
        for line in p.lines:
            world_pts = line.world_xyz(p.frame)
            dists = geometry.angular_distance(world_pts[:, None, :], seeds[None, :, :])
            nearest = np.argmin(dists, axis=1)
            assert np.all(nearest == p.plate_id)


def test_lines_are_evenly_spaced_in_phi():
    plates = generate_plates(seed=4, num_plates=6)
    for p in plates:
        phis = sorted(line.phi for line in p.lines)
        if len(phis) < 2:
            continue
        diffs = np.diff(phis)
        # All gaps should be an integer multiple of the target spacing (some plates
        # won't own every consecutive row near their boundary).
        from app.plates import TARGET_LINE_SPACING_RAD

        ratios = diffs / TARGET_LINE_SPACING_RAD
        assert np.allclose(ratios, np.round(ratios), atol=1e-6)


def test_generate_world_matches_plate_count():
    world = generate_world(seed=7, num_plates=9)
    assert len(world.plates) == 9
    assert world.elapsed_years == 0.0


def test_generation_is_deterministic_for_same_seed():
    w1 = generate_world(seed=123, num_plates=8)
    w2 = generate_world(seed=123, num_plates=8)
    assert len(w1.plates) == len(w2.plates)
    for p1, p2 in zip(w1.plates, w2.plates):
        assert p1.crust_type == p2.crust_type
        assert np.allclose(p1.frame, p2.frame)
        assert len(p1.lines) == len(p2.lines)


def test_generate_plates_without_num_plates_picks_a_plausible_count():
    for seed in range(20):
        plates = generate_plates(seed=seed)
        assert MIN_AUTO_PLATES <= len(plates) <= MAX_AUTO_PLATES


def test_generate_plates_auto_count_is_deterministic_for_same_seed():
    p1 = generate_plates(seed=99)
    p2 = generate_plates(seed=99)
    assert len(p1) == len(p2)
    for a, b in zip(p1, p2):
        assert a.crust_type == b.crust_type
        assert np.allclose(a.frame, b.frame)


def test_generate_plates_num_continents_gives_exact_continental_count():
    for n in range(1, 6):
        plates = generate_plates(seed=3, num_plates=12, num_continents=n)
        continental = [p for p in plates if p.crust_type == "continental"]
        assert len(continental) == n


def test_generate_plates_num_continents_bumps_up_total_plate_count_if_needed():
    plates = generate_plates(seed=4, num_plates=5, num_continents=5)
    assert len(plates) >= 5 + MIN_OCEANIC_PLATES
    assert sum(1 for p in plates if p.crust_type == "continental") == 5


def test_generate_plates_num_continents_is_clamped_to_max():
    from app.plates import MAX_CONTINENTS

    plates = generate_plates(seed=5, num_plates=10, num_continents=999)
    assert sum(1 for p in plates if p.crust_type == "continental") == MAX_CONTINENTS


def test_outline_world_traces_a_loop_covering_every_line():
    plates = generate_plates(seed=5, num_plates=8)
    for p in plates:
        lines_with_nodes = [line for line in p.lines if len(line.theta) > 0]
        outline = p.outline_world()
        assert len(outline) == 2 * len(lines_with_nodes)
        assert np.allclose(np.linalg.norm(outline, axis=-1), 1.0, atol=1e-9)


def test_outline_world_empty_for_plate_with_no_lines():
    plates = generate_plates(seed=6, num_plates=8)
    p = plates[0]
    p.lines = []
    assert len(p.outline_world()) == 0

"""gaps.fill_gaps -- whole-sphere coverage maintenance (spawning new oceanic crust into a
region no plate currently covers, e.g. one vacated by a fully-subducted plate)."""

import numpy as np
from scipy.spatial import cKDTree

from app import gaps, geometry
from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from app.lithosphere_plate import generate_plates, new_plate
from app.world import World


def _small_world(seed=1, num_plates=4, node_density=1.0):
    plates = generate_plates(seed=seed, num_plates=num_plates, node_density=node_density)
    return World(seed=seed, plates=plates, next_plate_id=len(plates), node_density=node_density, mantle_centers=[])


def test_fill_gaps_is_a_no_op_on_a_freshly_generated_world():
    """`generate_plates`' own Voronoi tiling already covers the whole sphere with no gaps by
    construction -- fill_gaps should find nothing to do."""
    world = _small_world()
    assert gaps.fill_gaps(world) == []
    assert len(world.plates) == 4


def test_fill_gaps_spawns_a_new_plate_covering_a_removed_plates_vacated_territory():
    world = _small_world()
    removed = world.plates.pop(1)
    removed_points, _ = removed.all_points_and_elevation()
    n_before = world.next_plate_id

    events = gaps.fill_gaps(world)

    assert len(events) == 1
    assert "New oceanic crust" in events[0]
    assert len(world.plates) == 4
    assert world.next_plate_id == n_before + 1

    spawned = world.plates[-1]
    assert spawned.crust_type == "oceanic"
    assert spawned.node_count() > 0.9 * removed.node_count()  # near-full replacement, not a sliver

    # The vacated territory is now actually covered: every point the removed plate used to
    # own has some live plate's node nearby.
    all_points = np.concatenate([p.all_points_and_elevation()[0] for p in world.plates], axis=0)
    tree = cKDTree(all_points)
    coverage_radius_rad = gaps.COVERAGE_RADIUS_MULT * line_spacing_rad(world.node_density)
    dist, _ = tree.query(removed_points)
    assert np.all(dist < coverage_radius_rad)


def test_fill_gaps_ignores_a_gap_smaller_than_the_minimum(monkeypatch):
    """A gap cluster below MIN_GAP_NODES is ordinary per-step catch-up lag, not a genuine
    void -- reacting to it would spawn a sliver plate at every busy divergent boundary every
    interval (see the module docstring)."""
    world = _small_world()
    removed = world.plates.pop(1)

    # Raise the floor above this removed plate's own node count -- from fill_gaps' own
    # perspective this is the same "a gap cluster this size is unremarkable" case, just
    # exercised without needing to fabricate a precisely-sized synthetic gap.
    monkeypatch.setattr(gaps, "MIN_GAP_NODES", removed.node_count() * 10)

    events = gaps.fill_gaps(world)
    assert events == []
    assert len(world.plates) == 3


def test_fill_gaps_handles_a_world_with_no_plates():
    world = World(seed=1, plates=[], next_plate_id=0, mantle_centers=[])
    assert gaps.fill_gaps(world) == []


def _lone_continental_island(seed=1, node_density=1.0, radius_rad=0.4):
    """A single small continental plate covering a disk around the local north pole
    (forced above sea level regardless of the noise field), with the rest of the sphere
    left as one big uncovered gap -- lets a gap-fill test control exactly where "land" is
    without depending on generate_plates' own random tiling."""
    frame = np.eye(3)
    spacing_rad = line_spacing_rad(node_density)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        return geometry.angular_distance(world_pts, frame[:, 2]) < radius_rad

    plate = new_plate(0, frame, "continental", spacing_rad, seed, is_owned=is_owned)
    n = plate.node_count()
    plate.set_fields_on_plate(elevation=np.full(n, 800.0))  # force unambiguous dry land
    world = World(seed=seed, plates=[plate], next_plate_id=1, node_density=node_density, mantle_centers=[])
    return world, spacing_rad


def test_fill_gaps_adopts_continental_type_only_right_at_a_real_coastline():
    world, spacing_rad = _lone_continental_island()
    events = gaps.fill_gaps(world)
    assert len(events) == 1

    spawned = world.plates[-1]
    assert spawned.node_count() > 0
    points, _ = spawned.all_points_and_elevation()
    is_continental = effective_is_continental_from_codes(spawned.collect("crust_type_code"), spawned.crust_type == "continental")

    # Nodes right at the island's own coastline (just past the parent plate's coverage
    # radius) adopt continental type; the antipodal deep-ocean interior of the gap does not
    # -- a gap virtually never turns land except where it's genuinely landlocked by a
    # still-standing coast (see gaps.GAP_LAND_ADOPTION_RADIUS_MULT).
    dist_from_pole = geometry.angular_distance(points, np.array([0.0, 0.0, 1.0]))
    near_coast = dist_from_pole < 0.4 + gaps.GAP_LAND_ADOPTION_RADIUS_MULT * spacing_rad
    far_interior = dist_from_pole > np.pi - 0.3  # near the antipode of the island

    assert np.any(is_continental[near_coast])
    assert not np.any(is_continental[far_interior])
    # The gap is overwhelmingly open ocean -- a small coastal fringe shouldn't flip the
    # spawned plate's own majority-vote label.
    assert spawned.crust_type == "oceanic"


def test_fill_gaps_mid_ocean_gap_stays_all_oceanic():
    """Regression guard: with no continental plate anywhere in the world, the new
    land-adoption rule has nothing to adopt from and every gap node comes back oceanic,
    exactly as before that rule existed."""
    frame = np.eye(3)
    node_density = 1.0
    spacing_rad = line_spacing_rad(node_density)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        return geometry.angular_distance(world_pts, frame[:, 2]) < 0.4

    plate = new_plate(0, frame, "oceanic", spacing_rad, seed=1, is_owned=is_owned)
    world = World(seed=1, plates=[plate], next_plate_id=1, node_density=node_density, mantle_centers=[])

    events = gaps.fill_gaps(world)
    assert len(events) == 1
    spawned = world.plates[-1]
    is_continental = effective_is_continental_from_codes(spawned.collect("crust_type_code"), spawned.crust_type == "continental")
    assert not np.any(is_continental)
    assert spawned.crust_type == "oceanic"

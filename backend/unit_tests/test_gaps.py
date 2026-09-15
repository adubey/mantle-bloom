"""gaps.py -- whole-sphere coverage maintenance: growing the plate(s) adjacent to a region no
plate currently covers into it (e.g. one vacated by a fully-subducted plate), falling back to
spawning a brand-new plate only when nothing is adjacent."""

import numpy as np
from scipy.spatial import cKDTree

from app import gaps, geometry
from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from app.lithosphere_plate import generate_plates, new_plate
from app.world import World


def _small_world(seed=1, num_plates=4, node_density=1.0):
    plates = generate_plates(seed=seed, num_plates=num_plates, node_density=node_density)
    return World(seed=seed, plates=plates, next_plate_id=len(plates), node_density=node_density, mantle_centers=[])


def test_fill_gaps_by_growing_neighbours_handles_a_world_with_no_plates():
    world = World(seed=1, plates=[], next_plate_id=0, mantle_centers=[])
    assert gaps.fill_gaps_by_growing_neighbours(world) == []


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


def _spawn_over_whole_gap(world: "World"):
    """Exercises `gaps._spawn_plate_from_gap` directly over this world's one big uncovered
    region -- the same crust-typing logic `fill_gaps_by_growing_neighbours` falls back to when
    a gap cluster has no adjacent plate, without depending on that whole detect/cluster/adjacency
    pipeline (these fixtures' single small plate leaves the rest of the sphere as one cluster
    that *is* adjacent to it, so the fallback path itself isn't what's under test here)."""
    existing_context = gaps._existing_node_tree(world)
    spacing_rad = line_spacing_rad(world.node_density)
    gap_points = gaps._find_gap_points(existing_context, spacing_rad)
    return gaps._spawn_plate_from_gap(world, gap_points, spacing_rad, existing_context)


def test_fill_gaps_adopts_continental_type_only_right_at_a_real_coastline():
    world, spacing_rad = _lone_continental_island()
    spawned = _spawn_over_whole_gap(world)
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

    spawned = _spawn_over_whole_gap(world)
    is_continental = effective_is_continental_from_codes(spawned.collect("crust_type_code"), spawned.crust_type == "continental")
    assert not np.any(is_continental)
    assert spawned.crust_type == "oceanic"


# -- gap-age tracking (world.gap_tracks) --------------------------------------------------


def test_reconcile_gap_tracks_is_a_no_op_on_a_freshly_generated_world():
    world = _small_world()
    gaps.reconcile_gap_tracks(world)
    assert world.gap_tracks == []


def test_reconcile_gap_tracks_persists_first_seen_years_and_accumulates_steps_seen():
    world = _small_world()
    world.plates.pop(1)  # opens a real, stable gap -- nothing refills it between calls
    world.elapsed_years = 1_000_000.0

    gaps.reconcile_gap_tracks(world)
    assert len(world.gap_tracks) == 1
    track = world.gap_tracks[0]
    assert track.first_seen_years == 1_000_000.0
    assert track.steps_seen == 1

    world.elapsed_years = 5_000_000.0
    gaps.reconcile_gap_tracks(world)
    assert len(world.gap_tracks) == 1
    track = world.gap_tracks[0]
    # Same gap, re-identified by centroid proximity -- first_seen_years is carried forward,
    # not reset, and steps_seen accumulates.
    assert track.first_seen_years == 1_000_000.0
    assert track.last_seen_years == 5_000_000.0
    assert track.steps_seen == 2


def test_reconcile_gap_tracks_drops_a_track_once_the_gap_heals():
    world = _small_world()
    removed = world.plates.pop(1)
    world.elapsed_years = 1_000_000.0
    gaps.reconcile_gap_tracks(world)
    assert len(world.gap_tracks) == 1

    # Heal the gap the same way fill_gaps_by_growing_neighbours would -- put a plate back over
    # the vacated ground.
    world.plates.append(removed)
    world.elapsed_years = 2_000_000.0
    gaps.reconcile_gap_tracks(world)
    assert world.gap_tracks == []


def test_reconcile_gap_tracks_surfaces_a_cluster_smaller_than_min_gap_nodes(monkeypatch):
    """The whole point of a separate, much lower GAP_AGE_MIN_CLUSTER_NODES floor: age-tracking
    must never hide a small persistent notch just because it's too small for
    fill_gaps_by_growing_neighbours to ever act on (MIN_GAP_NODES) -- same fixture as
    test_fill_gaps_by_growing_neighbours_ignores_a_gap_smaller_than_the_minimum, but checking
    the tracker instead."""
    world = _small_world()
    removed = world.plates.pop(1)
    monkeypatch.setattr(gaps, "MIN_GAP_NODES", removed.node_count() * 10)

    gaps.reconcile_gap_tracks(world)
    assert len(world.gap_tracks) >= 1
    # fill_gaps_by_growing_neighbours itself declines under the same raised floor -- confirms
    # the two floors are genuinely independent (GAP_AGE_MIN_CLUSTER_NODES doesn't move just
    # because MIN_GAP_NODES did).
    assert gaps.fill_gaps_by_growing_neighbours(world) == []


# -- fill_gaps_by_growing_neighbours -------------------------------------------------------


def test_fill_gaps_by_growing_neighbours_is_a_no_op_on_a_freshly_generated_world():
    """`generate_plates`' own Voronoi tiling already covers the whole sphere with no gaps by
    construction -- fill_gaps_by_growing_neighbours should find nothing to do."""
    world = _small_world()
    assert gaps.fill_gaps_by_growing_neighbours(world) == []
    assert len(world.plates) == 4


def test_fill_gaps_by_growing_neighbours_grows_existing_adjacent_plates_not_a_new_one():
    """The whole point of this algorithm: a removed plate's vacated territory is absorbed by
    its real, still-live neighbours -- no new plate is spawned."""
    world = _small_world()
    removed = world.plates.pop(1)
    removed_points, _ = removed.all_points_and_elevation()
    n_plates_before = len(world.plates)
    next_plate_id_before = world.next_plate_id

    events = gaps.fill_gaps_by_growing_neighbours(world)

    assert len(events) > 0
    assert all("grew by" in e for e in events)
    assert len(world.plates) == n_plates_before  # no new plate spawned
    assert world.next_plate_id == next_plate_id_before

    all_points = np.concatenate([p.all_points_and_elevation()[0] for p in world.plates], axis=0)
    tree = cKDTree(all_points)
    coverage_radius_rad = gaps.COVERAGE_RADIUS_MULT * line_spacing_rad(world.node_density)
    dist, _ = tree.query(removed_points)
    assert np.mean(dist < coverage_radius_rad) > 0.9  # near-full replacement, not a sliver


def test_fill_gaps_by_growing_neighbours_falls_back_to_spawning_when_nothing_is_adjacent(monkeypatch):
    """A gap cluster with no plate anywhere near it (e.g. every bordering plate that used to
    border it has itself fully vanished) has nothing to grow, so this falls back to spawning a
    new plate. On a real sphere a gap cluster's own boundary almost always touches *some*
    still-live plate within
    ADJACENT_PLATE_REACH_MULT (that's the common case the two tests above exercise) -- this
    directly forces the "found nothing nearby" branch instead of trying to construct that rare
    geometry from scratch."""
    world = _small_world()
    world.plates.pop(1)
    monkeypatch.setattr(gaps, "_adjacent_plates_to_cluster", lambda *a, **k: [])
    n_before = world.next_plate_id

    events = gaps.fill_gaps_by_growing_neighbours(world)

    assert len(events) == 1
    assert "no adjacent plate to grow" in events[0]
    assert len(world.plates) == 4
    assert world.next_plate_id == n_before + 1


def test_fill_gaps_by_growing_neighbours_ignores_a_gap_smaller_than_the_minimum(monkeypatch):
    world = _small_world()
    removed = world.plates.pop(1)
    monkeypatch.setattr(gaps, "MIN_GAP_NODES", removed.node_count() * 10)

    events = gaps.fill_gaps_by_growing_neighbours(world)
    assert events == []
    assert len(world.plates) == 3

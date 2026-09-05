"""gaps.fill_gaps -- whole-sphere coverage maintenance (spawning new oceanic crust into a
region no plate currently covers, e.g. one vacated by a fully-subducted plate)."""

import numpy as np
from scipy.spatial import cKDTree

from app import gaps
from app.elevation_lines import line_spacing_rad
from app.lithosphere_plate import generate_plates
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

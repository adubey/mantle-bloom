import numpy as np
from scipy.spatial import cKDTree
from app import gaps, geometry
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def _old_plate(plate_id: int) -> Plate:
    return Plate(plate_id=plate_id, frame=np.eye(3), crust_type="oceanic", age_steps=gaps.YOUNG_PLATE_AGE_STEPS + 1)


def _young_plate(plate_id: int) -> Plate:
    return Plate(plate_id=plate_id, frame=np.eye(3), crust_type="oceanic", age_steps=0)


def test_fill_gaps_leaves_no_large_gaps_after_a_long_simulation():
    # node_density=2.0 (half the default 4.0, see plates.NODE_DENSITY_CHOICES) -- this test
    # only checks that no large *gap* remains, a geometric property that doesn't need the
    # default's full node density to verify.
    world = generate_world(seed=55, num_plates=10, node_density=2.0)
    # gaps.py only reacts to plate/node geometry, never to climate/erosion/hydrology output
    # (see World.simulate_climate_biomes) -- skipping that per-step computation cuts this
    # long-running test's cost substantially without changing what it exercises.
    world.simulate_climate_biomes = False
    for _ in range(25):
        step_world(world, years=6_000_000)
    # A single fill_gaps call only claims one ring of a gap (see GROWTH_RING_RAD) --
    # gradual healing, not an instant full close -- so give it several passes, matching how
    # step_world actually calls it repeatedly over time.
    for _ in range(10):
        gaps.fill_gaps(world)

    existing_points, _ = gaps._all_existing_points(world)
    tree = cKDTree(existing_points)
    remaining_gap_points = gaps._find_gap_points(tree)
    if len(remaining_gap_points) == 0:
        return
    labels = gaps.cluster_points(remaining_gap_points, gaps.CLUSTER_RADIUS_RAD)
    counts = np.bincount(labels)
    assert counts.max() < gaps.MIN_GAP_POINTS


def test_gap_fill_is_deterministic_for_same_seed_and_step_sequence():
    def run():
        w = generate_world(seed=77, num_plates=10, node_density=2.0)  # see the test above for why this is safe
        w.simulate_climate_biomes = False  # see the test above for why this is safe here
        for _ in range(12):
            step_world(w, years=5_000_000)
        return w

    w1, w2 = run(), run()
    p1s = sorted(w1.plates, key=lambda p: p.plate_id)
    p2s = sorted(w2.plates, key=lambda p: p.plate_id)
    assert [p.plate_id for p in p1s] == [p.plate_id for p in p2s]
    for p1, p2 in zip(p1s, p2s):
        assert p1.crust_type == p2.crust_type
        assert np.allclose(p1.frame, p2.frame)
        assert p1.node_count() == p2.node_count()

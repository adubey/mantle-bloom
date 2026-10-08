"""merge_split.update_overlap_tracking / plates.compute_node_overlap / the
`overlap_onset_years` field."""

import numpy as np

from app import geometry, lithosphere, merge_split
from app.lithosphere_plate import growth_seed_thickness
from app.plates import compute_node_overlap
from app.world import World, generate_world, step_world

from .quad_fixtures import block_keys, quad_plate


def _plate(plate_id, filler_column, crust_type="continental", seed_xyz=(1.0, 0.0, 0.0)):
    """Eight cells on the face-centre row, plus two filler cells further north starting at
    `filler_column` -- so in node order the row comes first, then the filler."""
    keys = np.concatenate([block_keys(range(8), [0]), block_keys([filler_column, filler_column + 1], [20])])
    frame = geometry.plate_frame_from_seed(np.asarray(seed_xyz, dtype=float))
    return quad_plate(plate_id, crust_type, frame=frame, keys=keys)


def _overlapping_world(seed=1):
    """Plate B's row sits exactly on top of plate A's (same frame, same cells), so every row
    cell of each lies inside the other. Their filler cells are far apart, so those 2 cells
    each are NOT co-located."""
    a = _plate(0, filler_column=20)
    b = _plate(1, filler_column=-20)
    return World(seed=seed, plates=[a, b], next_plate_id=2)


def test_compute_node_overlap_flags_colocated_nodes_both_ways():
    world = _overlapping_world()
    overlap = compute_node_overlap(world.plates)

    for pid, other in ((0, 1), (1, 0)):
        info = overlap[pid]
        # The 8 row cells overlap; the 2 filler cells don't.
        assert info["overlap_mask"][:8].all()
        assert not info["overlap_mask"][8:].any()
        assert info["by_partner"] == {other: 8}
        np.testing.assert_array_equal(info["cover_count"], [1] * 8 + [0, 0])
        np.testing.assert_array_equal(info["continental_cover_count"], [1] * 8 + [0, 0])


def test_compute_node_overlap_cover_count_counts_every_plate_on_a_node():
    # A third plate stacked on the same row: each row node now sits on two other plates.
    world = _overlapping_world()
    world.plates.append(_plate(2, filler_column=0, crust_type="oceanic"))
    overlap = compute_node_overlap(world.plates)

    for pid in (0, 1, 2):
        np.testing.assert_array_equal(overlap[pid]["cover_count"], [2] * 8 + [0, 0])
    # Plate 2 is oceanic: the continental plates see one continental partner, plate 2 sees two.
    for pid, expected in ((0, 1), (1, 1), (2, 2)):
        np.testing.assert_array_equal(overlap[pid]["continental_cover_count"], [expected] * 8 + [0, 0])


def test_update_overlap_tracking_stamps_once_then_clears():
    world = _overlapping_world()
    world.elapsed_years = 5_000_000.0
    merge_split.update_overlap_tracking(world, 100_000.0)
    onset0 = world.plates[0].collect("overlap_onset_years")
    assert np.all(onset0[:8] == 5_000_000.0)
    assert np.all(onset0[8:] == 0.0)

    # Still overlapping a step later -> the onset year is not overwritten.
    world.elapsed_years = 9_000_000.0
    merge_split.update_overlap_tracking(world, 100_000.0)
    assert np.all(world.plates[0].collect("overlap_onset_years")[:8] == 5_000_000.0)

    # Move B off A entirely -> the stamp clears back to 0.
    world.plates[1] = _plate(1, filler_column=-20, seed_xyz=(-1.0, 0.0, 0.0))
    merge_split.update_overlap_tracking(world, 100_000.0)
    assert np.all(world.plates[0].collect("overlap_onset_years") == 0.0)


def test_growth_seed_thickness_is_oceanic_and_submerged():
    hc, hm = growth_seed_thickness()
    assert (hc, hm) == (lithosphere.REFERENCE_HC_OCEANIC_M, lithosphere.YOUNG_RIDGE_HM_M)
    # Even seeded onto a continental plate's crust density, a fresh grown node is deep ocean,
    # not +200 m land -- this is the land-area runaway fix.
    z = lithosphere.isostatic_elevation(np.array([hc]), np.array([hm]), lithosphere.RHO_CONTINENTAL_CRUST)
    assert z[0] < -1000.0


def test_stepping_a_world_keeps_land_bounded_and_populates_overlap_onset():
    world = generate_world(7, node_density=1.0)

    def land_fraction(w):
        elev = np.concatenate([p.all_points_and_elevation()[1] for p in w.plates])
        return float(np.mean(elev > w.sea_level_m))

    start = land_fraction(world)
    for _ in range(30):
        step_world(world, 1_000_000.0)

    # The pre-fix engine ran land fraction up ~0.07 in 30 Myr at this density; it must not
    # run away upward any more.
    assert land_fraction(world) < start + 0.04

    onset = np.concatenate([p.collect("overlap_onset_years") for p in world.plates])
    # Some transient envelope overlap always exists, and it must now carry an onset year.
    assert np.any(onset > 0.0)
    assert np.all(onset[onset > 0.0] <= world.elapsed_years)

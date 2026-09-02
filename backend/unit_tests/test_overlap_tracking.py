"""merge_split.update_overlap_tracking / plates.compute_node_overlap / the
ElevationLine.overlap_onset_years field and its pickle backfill."""

import pickle

import numpy as np

from app import geometry, lithosphere, merge_split
from app.elevation_lines import ElevationLine
from app.lithosphere_plate import growth_seed_thickness
from app.plates import PlateWithLines, compute_node_overlap, line_spacing_rad
from app.world import World, generate_world, step_world


def _plate(plate_id, seed_xyz, theta, phi=0.0, crust_type="continental", filler_phi=1.0):
    frame = geometry.plate_frame_from_seed(np.asarray(seed_xyz, dtype=float))
    line = ElevationLine(phi=phi, theta=np.asarray(theta, dtype=float), elevation=np.zeros(len(theta)))
    filler = ElevationLine(phi=filler_phi, theta=np.array([0.0, 0.1]), elevation=np.zeros(2))
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=[line, filler])


def _overlapping_world(seed=1):
    """Plate B's phi=0 row is placed right on top of plate A's phi=0 row (same frame, same
    thetas) so every one of B's row nodes reads as overlapping A and vice versa. Their filler
    rows sit at opposite latitudes so those 2 nodes each are NOT co-located."""
    theta = np.linspace(-0.05, 0.05, 8)
    a = _plate(0, [1.0, 0.0, 0.0], theta, filler_phi=1.0)
    b = _plate(1, [1.0, 0.0, 0.0], theta, filler_phi=-1.0)
    return World(seed=seed, plates=[a, b], next_plate_id=2)


def test_compute_node_overlap_flags_colocated_nodes_both_ways():
    world = _overlapping_world()
    tol = 0.5 * line_spacing_rad(world.node_density)
    overlap = compute_node_overlap(world.plates, tol)

    for pid, other in ((0, 1), (1, 0)):
        info = overlap[pid]
        # The 8 phi=0 nodes overlap; the 2 filler nodes at phi=1 don't.
        assert info["overlap_mask"][:8].all()
        assert not info["overlap_mask"][8:].any()
        assert info["by_partner"] == {other: 8}


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
    world.plates[1] = _plate(1, [-1.0, 0.0, 0.0], np.linspace(-0.05, 0.05, 8))
    merge_split.update_overlap_tracking(world, 100_000.0)
    assert np.all(world.plates[0].collect("overlap_onset_years") == 0.0)


def test_overlap_onset_years_survives_masked_and_pickle_backfill():
    line = ElevationLine(phi=0.0, theta=np.arange(5.0), elevation=np.zeros(5))
    line.set_fields(overlap_onset_years=np.array([0.0, 1.0, 2.0, 3.0, 4.0]))
    kept = line.masked(np.array([1, 3]))
    assert list(kept.overlap_onset_years) == [1.0, 3.0]

    # Simulate a save written before the field existed: drop the backing array, round-trip.
    stale = pickle.loads(pickle.dumps(line))
    del stale.__dict__["_overlap_onset_years"]
    restored = pickle.loads(pickle.dumps(stale))
    assert list(restored.overlap_onset_years) == [0.0] * 5


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

"""Lateral magma transport (GitHub issue #205) -- the transport/deposit half in
magma_transport.py (the source half, rheology.magma_export_strength_and_volume, has its own
tests in test_rheology.py)."""

import numpy as np
import pytest

from app import geometry, lithosphere, magma_transport
from app.elevation_lines import ELEV_CHANGE_LATERAL_MAGMA, line_spacing_rad
from app.lithosphere_plate import new_plate
from app.world import World


def _dest_index(xyz, hc, plate_id=None, line_index=None, node_index=None):
    n = len(xyz)
    return magma_transport._ContinentalNodeIndex(
        np.asarray(xyz, dtype=float),
        np.zeros(n, dtype=int) if plate_id is None else plate_id,
        np.zeros(n, dtype=int) if line_index is None else line_index,
        np.arange(n) if node_index is None else node_index,
        np.asarray(hc, dtype=float),
    )


def _range_rad():
    return magma_transport.MAGMA_TRANSPORT_RANGE_KM / lithosphere.PLANET_RADIUS_KM


def _two_adjacent_continental_plates(radius_rad=0.15, node_density=1.0, seed=1):
    """Two small continental plates sharing a seam near the north pole, split by x >= 0 vs
    x < 0 -- close enough (well under MAGMA_TRANSPORT_RANGE_KM) that a parcel generated on one
    plate can legitimately reach a destination on the other, the scenario cross-plate
    scatter-write needs to be exercised against. Every node starts at exactly
    REFERENCE_HC_CONTINENTAL_M (overriding new_plate's own relief-linked Hc noise) so a test
    can thin exactly one node and know, by construction, that every other node has zero
    destination weight (see magma_transport._weighted_destination_pairs)."""
    frame = np.eye(3)
    spacing_rad = line_spacing_rad(node_density)

    def in_disk(pts):
        return geometry.angular_distance(pts, frame[:, 2]) < radius_rad

    def is_owned_a(pts):
        return in_disk(pts) & (pts[:, 0] >= 0.0)

    def is_owned_b(pts):
        return in_disk(pts) & (pts[:, 0] < 0.0)

    plate_a = new_plate(0, frame, "continental", spacing_rad, seed, is_owned=is_owned_a)
    plate_b = new_plate(1, frame, "continental", spacing_rad, seed, is_owned=is_owned_b)
    plate_a.set_fields_on_plate(crustal_thickness_m=np.full(plate_a.node_count(), lithosphere.REFERENCE_HC_CONTINENTAL_M))
    plate_b.set_fields_on_plate(crustal_thickness_m=np.full(plate_b.node_count(), lithosphere.REFERENCE_HC_CONTINENTAL_M))

    world = World(seed=seed, plates=[plate_a, plate_b], next_plate_id=2, node_density=node_density, mantle_centers=[])
    return world, plate_a, plate_b


def _thin_nearest_node(plate, target_xyz, hc_m):
    """Sets the single node on `plate` nearest `target_xyz` to `hc_m`, returning
    `(line_index, node_index, actual_xyz)` -- the address a test needs to later check that
    exact node's own elev_change_reason/crustal_thickness_m."""
    points, _ = plate.all_points_and_elevation()
    global_idx = int(np.argmin(np.linalg.norm(points - target_xyz, axis=1)))
    offset = 0
    for line_index, line in enumerate(plate.lines):
        n = len(line)
        if n == 0:
            continue
        if offset <= global_idx < offset + n:
            node_index = global_idx - offset
            hc = line.crustal_thickness_m.copy()
            hc[node_index] = hc_m
            plate.replace_line(line_index, line.replace(crustal_thickness_m=hc))
            return line_index, node_index, points[global_idx]
        offset += n
    raise AssertionError("node not found")


# -- destination weighting -----------------------------------------------------------------


def test_weighted_destinations_favours_thinner_crust_at_equal_distance():
    origin = np.array([1.0, 0.0, 0.0])
    thin = geometry.normalize(origin + np.array([0.0, 0.002, 0.0]))
    thick = geometry.normalize(origin + np.array([0.0, -0.002, 0.0]))
    dest_index = _dest_index(np.array([thin, thick]), np.array([10_000.0, 34_000.0]))

    parcel_idx, dest_idx, weight = magma_transport._weighted_destination_pairs(np.array([origin]), dest_index, _range_rad())

    assert set(dest_idx.tolist()) == {0, 1}
    thin_weight = float(weight[dest_idx == 0][0])
    thick_weight = float(weight[dest_idx == 1][0])
    assert thin_weight > thick_weight


def test_weighted_destinations_ignores_a_candidate_at_or_above_reference_thickness():
    origin = np.array([1.0, 0.0, 0.0])
    at_reference = geometry.normalize(origin + np.array([0.0, 0.002, 0.0]))
    dest_index = _dest_index(np.array([at_reference]), np.array([lithosphere.REFERENCE_HC_CONTINENTAL_M]))

    _, dest_idx, _ = magma_transport._weighted_destination_pairs(np.array([origin]), dest_index, _range_rad())
    assert len(dest_idx) == 0


def test_weighted_destinations_ignores_a_candidate_past_the_transport_range():
    origin = np.array([1.0, 0.0, 0.0])
    far_rad = 1.5 * _range_rad()
    beyond_range = geometry.normalize(origin + np.array([0.0, far_rad, 0.0]))
    dest_index = _dest_index(np.array([beyond_range]), np.array([5_000.0]))

    _, dest_idx, _ = magma_transport._weighted_destination_pairs(np.array([origin]), dest_index, _range_rad())
    assert len(dest_idx) == 0


def test_fixed_k_matches_radius_search_when_all_candidates_fit():
    origin = np.array([1.0, 0.0, 0.0])
    destinations = geometry.normalize(np.array([
        [1.0, 0.001, 0.0], [1.0, -0.002, 0.0], [1.0, 0.003, 0.0],
    ]))
    dest_index = _dest_index(destinations, np.array([10_000.0, 20_000.0, 30_000.0]))
    exact = magma_transport._weighted_destination_pairs(np.array([origin]), dest_index, _range_rad())
    capped = magma_transport._weighted_destination_pairs(np.array([origin]), dest_index, _range_rad(), 3)
    for old, new in zip(exact, capped):
        np.testing.assert_array_equal(old, new)


def test_fixed_k_limits_destinations_per_parcel():
    origin = np.array([1.0, 0.0, 0.0])
    destinations = geometry.normalize(np.array([
        [1.0, 0.001, 0.0], [1.0, 0.002, 0.0], [1.0, 0.003, 0.0],
    ]))
    dest_index = _dest_index(destinations, np.array([10_000.0, 20_000.0, 30_000.0]))
    _, dest_idx, _ = magma_transport._weighted_destination_pairs(np.array([origin]), dest_index, _range_rad(), 2)
    assert dest_idx.tolist() == [0, 1]


def test_continental_node_index_excludes_oceanic_nodes():
    """GitHub issue #205's own destination filter: an oceanic destination is just ordinary
    seafloor volcanism, not the land-fraction fix this exists for -- a mixed-composition
    plate's oceanic-typed nodes must never appear in the whole-sphere destination index."""
    frame = np.eye(3)
    spacing_rad = line_spacing_rad(1.0)

    def is_owned(pts):
        return geometry.angular_distance(pts, frame[:, 2]) < 0.15

    def node_is_continental(pts):
        return pts[:, 0] >= 0.0  # half continental, half oceanic

    plate = new_plate(0, frame, "continental", spacing_rad, seed=1, is_owned=is_owned, node_is_continental=node_is_continental)
    world = World(seed=1, plates=[plate], next_plate_id=1, node_density=1.0, mantle_centers=[])
    _thin_nearest_node(plate, geometry.normalize(np.array([0.05, 0.0, 1.0])), hc_m=5_000.0)

    dest_index = magma_transport._build_continental_node_index(world)
    assert len(dest_index) > 0
    assert np.all(dest_index.xyz[:, 0] >= 0.0)


def test_continental_node_index_keeps_only_nodes_with_positive_destination_weight():
    world, _, plate_b = _two_adjacent_continental_plates()
    target = geometry.normalize(np.array([-0.05, 0.0, 1.0]))
    line_index, node_index, actual_xyz = _thin_nearest_node(plate_b, target, hc_m=5_000.0)

    dest_index = magma_transport._build_continental_node_index(world)

    assert len(dest_index) == 1
    assert int(dest_index.plate_id[0]) == plate_b.plate_id
    assert int(dest_index.line_index[0]) == line_index
    assert int(dest_index.node_index[0]) == node_index
    np.testing.assert_array_equal(dest_index.xyz[0], actual_xyz)


def test_run_magma_transport_is_a_no_op_with_no_pending_parcels():
    world, _, _ = _two_adjacent_continental_plates()
    assert world.pending_magma_parcels == []
    assert magma_transport.run_magma_transport(world, banked_myr=4.0) == []


# -- end-to-end transport + deposit ----------------------------------------------------------


def test_deposit_writes_into_a_different_plates_line_and_stamps_lateral_magma():
    world, plate_a, plate_b = _two_adjacent_continental_plates()
    target = geometry.normalize(np.array([-0.05, 0.0, 1.0]))  # plate_b's side, near the seam
    line_index, node_index, actual_xyz = _thin_nearest_node(plate_b, target, hc_m=5_000.0)

    origin = geometry.normalize(np.array([0.05, 0.0, 1.0]))  # plate_a's side, near the seam
    world.pending_magma_parcels.append(
        magma_transport.MagmaParcel(origin_xyz=origin, volume_m3=1e13, step_generated=0)
    )

    events = magma_transport.run_magma_transport(world, banked_myr=4.0)

    assert len(events) == 1
    updated_line = plate_b.lines[line_index]
    assert updated_line.crustal_thickness_m[node_index] > 5_000.0
    assert updated_line.elev_change_reason[node_index] == ELEV_CHANGE_LATERAL_MAGMA


@pytest.mark.parametrize("max_destinations_per_parcel", [None, 1])
def test_global_cap_bounds_a_single_destination_regardless_of_parcel_count(max_destinations_per_parcel):
    """GitHub issue #205's round-1 review point 4 / the #145-reopening regression this design
    specifically closes: several independent parcels all targeting the same thinnest node must
    not each get their own separate rate-capped share -- the node has exactly one ceiling,
    enforced once, however many parcels reach it."""
    world, plate_a, plate_b = _two_adjacent_continental_plates()
    target = geometry.normalize(np.array([-0.05, 0.0, 1.0]))
    line_index, node_index, _ = _thin_nearest_node(plate_b, target, hc_m=5_000.0)
    hc_before = plate_b.lines[line_index].crustal_thickness_m[node_index]

    origin = geometry.normalize(np.array([0.05, 0.0, 1.0]))
    # Five independent parcels, each individually large enough to blow past the rate cap on
    # its own -- if capping were per-parcel rather than global-per-destination, the node would
    # receive roughly 5x the intended ceiling.
    for _ in range(5):
        world.pending_magma_parcels.append(magma_transport.MagmaParcel(origin_xyz=origin, volume_m3=1e13, step_generated=0))

    banked_myr = 4.0
    magma_transport.run_magma_transport(world, banked_myr=banked_myr, max_destinations_per_parcel=max_destinations_per_parcel)

    hc_after = plate_b.lines[line_index].crustal_thickness_m[node_index]
    max_allowed = magma_transport.MAGMA_DEPOSIT_RATE_M_PER_MYR * banked_myr
    # The 5 parcels' combined request (well over max_allowed, by construction) must still be
    # clipped to a single ceiling -- not just bounded by it, but actually pinned at it, which
    # confirms the cap bound rather than merely never having been reached.
    assert np.isclose(hc_after - hc_before, max_allowed, rtol=1e-6)


def test_stale_parcel_with_no_headroom_is_dropped_after_max_age_cycles():
    world, plate_a, plate_b = _two_adjacent_continental_plates()
    target = geometry.normalize(np.array([-0.05, 0.0, 1.0]))
    _thin_nearest_node(plate_b, target, hc_m=5_000.0)

    origin = geometry.normalize(np.array([0.05, 0.0, 1.0]))
    world.pending_magma_parcels.append(magma_transport.MagmaParcel(origin_xyz=origin, volume_m3=1e10, step_generated=0))

    # banked_myr=0.0 -> a zero rate cap every firing, so this parcel can never find headroom
    # anywhere, however good its candidate destinations are.
    for _ in range(magma_transport.MAGMA_PARCEL_MAX_AGE_CYCLES):
        assert len(world.pending_magma_parcels) == 1
        magma_transport.run_magma_transport(world, banked_myr=0.0)

    assert world.pending_magma_parcels == []


def test_deposit_never_exceeds_a_destinations_own_headroom_to_the_ceiling():
    """A destination's headroom to MAX_CRUSTAL_THICKNESS_M -- not just the per-step rate cap
    (`cap_hc`) -- must bound both the actual Hc write (`_scatter_write_deposits` already clips
    there) and this pass's own accounting (`placed_per_parcel` / the logged deposited total),
    so an unusually large `banked_myr` (an outsized step's own `years`) can't make the pass
    credit more than what actually landed."""
    world, plate_a, plate_b = _two_adjacent_continental_plates()
    target = geometry.normalize(np.array([-0.05, 0.0, 1.0]))
    line_index, node_index, _ = _thin_nearest_node(plate_b, target, hc_m=34_000.0)

    origin = geometry.normalize(np.array([0.05, 0.0, 1.0]))
    world.pending_magma_parcels.append(magma_transport.MagmaParcel(origin_xyz=origin, volume_m3=1e16, step_generated=0))

    banked_myr = 1_000.0  # cap_hc = 300,000 m -- far past this node's own ~50,000 m headroom
    events = magma_transport.run_magma_transport(world, banked_myr=banked_myr)

    hc_after = plate_b.lines[line_index].crustal_thickness_m[node_index]
    assert hc_after == lithosphere.MAX_CRUSTAL_THICKNESS_M
    actual_delta = hc_after - 34_000.0

    assert len(events) == 1
    reported_m = float(events[0].split(" ")[4])
    assert np.isclose(reported_m, actual_delta, rtol=1e-3)


def test_fully_placed_parcel_is_removed_in_one_firing():
    world, plate_a, plate_b = _two_adjacent_continental_plates()
    target = geometry.normalize(np.array([-0.05, 0.0, 1.0]))
    _thin_nearest_node(plate_b, target, hc_m=5_000.0)

    origin = geometry.normalize(np.array([0.05, 0.0, 1.0]))
    # A tiny volume, well under one firing's rate-cap headroom -- placed in full immediately.
    world.pending_magma_parcels.append(magma_transport.MagmaParcel(origin_xyz=origin, volume_m3=1.0, step_generated=0))

    magma_transport.run_magma_transport(world, banked_myr=4.0)

    assert world.pending_magma_parcels == []

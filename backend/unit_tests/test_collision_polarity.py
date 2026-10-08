"""Collision polarity (issue #318): subduction evidence recorded on the plates where it
happens, persistent continental front records, frozen lower/upper decisions and masks, and
the boundary-search cache the prepass shares with deform() -- see collision_polarity.py."""

import copy

import numpy as np

from app import collision_polarity as cp
from app import geometry, lithosphere, mantle, persistence, quad_tectonics, torque
from app.elevation_lines import line_spacing_rad
from app.lithosphere_plate import CONTINENTAL_CONTESTED_RETREAT_MIN_RUN, boundary_context
from app.sparse_quad_patch import PlateWithSparseQuadPatch, cells_per_face_edge, pack_cell_keys, unpack_cell_keys
from app.world import World

DENSITY = 0.5
SPACING = line_spacing_rad(DENSITY)
N = cells_per_face_edge(SPACING)
J = (20, 32)
STEP_YEARS = 100_000.0


def _block(i_range, j_range=J) -> np.ndarray:
    jj, ii = np.meshgrid(np.arange(*j_range), np.arange(*i_range), indexing="ij")
    return pack_cell_keys(np.zeros(ii.size, dtype=int), ii.ravel(), jj.ravel())


def _plate(plate_id, keys, crust_type="continental") -> PlateWithSparseQuadPatch:
    hc, hm = lithosphere.reference_thickness(crust_type)
    count = len(keys)
    plate = PlateWithSparseQuadPatch(
        plate_id,
        np.eye(3),
        crust_type,
        N,
        keys,
        fields={"crustal_thickness_m": np.full(count, hc), "mantle_lithosphere_thickness_m": np.full(count, hm)},
    )
    lithosphere.sync_plate_elevation(plate)
    return plate


def _world(*plates) -> World:
    world = World(seed=0, plates=list(plates), next_plate_id=max(p.plate_id for p in plates) + 1, node_density=DENSITY)
    world.fault_deformation_mode = "smooth"
    return world


def _centroid(points: np.ndarray) -> np.ndarray:
    return geometry.normalize(np.mean(points, axis=0))


def _column(i: int, j_range=J) -> np.ndarray:
    """World centroid of lattice column `i` (identity frame)."""
    return _centroid(_plate(999, _block((i, i + 1), j_range)).all_points_and_elevation()[0])


def _plate_centroid(plate) -> np.ndarray:
    return _centroid(plate.all_points_and_elevation()[0])


def _drive(plate, target: np.ndarray, cm_per_yr: float) -> None:
    """Set `plate`'s omega so its centroid moves toward `target` at `cm_per_yr`."""
    c = _plate_centroid(plate)
    tangent = geometry.normalize(target - c * float(c @ target))
    plate.set_omega(np.cross(c, tangent) * mantle.cm_per_yr_to_rad_per_yr(cm_per_yr))


def _move(plate, src: np.ndarray, dst: np.ndarray) -> None:
    """Rotate `plate` rigidly so world point `src` lands on `dst` -- its evidence rides along."""
    axis = geometry.normalize(np.cross(src, dst))
    angle = float(np.arccos(np.clip(src @ dst, -1.0, 1.0)))
    rotation = geometry.rotation_matrix(axis, angle)
    if float((rotation @ src) @ dst) < 1.0 - 1e-9:
        rotation = geometry.rotation_matrix(axis, -angle)
    plate.rotate(rotation)


def _step(world: World) -> cp.PolarityFrame:
    frame = cp.observe_contacts(world, STEP_YEARS)
    world.elapsed_years += STEP_YEARS
    return frame


def _records(world: World, pair) -> list[cp.CollisionFront]:
    return [r for r in world.collision_fronts if r.plate_ids == pair]


def _near(plate, other, spacings: float = 1.5) -> np.ndarray:
    """`plate`'s node positions within `spacings` of `other`'s nodes."""
    points = plate.all_points_and_elevation()[0]
    dist, _ = other.get_node_kdtree().query(points)
    return points[dist <= spacings * SPACING]


def _remove_plate(world: World, plate_id: int) -> None:
    frames = cp.begin_topology(world)
    world.plates = [p for p in world.plates if p.plate_id != plate_id]
    cp.end_topology(world, frames)


# --- Evidence and decisions ------------------------------------------------------------------


def test_third_plate_ocean_closure_keeps_its_slab_history_on_the_upper_continent():
    # A (continental) overrides ocean O; continent B rides in behind O. O vanishes before A
    # and B touch, so the polarity must come from the history left on A.
    a = _plate(1, _block((10, 20)))
    o = _plate(2, _block((20, 30)), "oceanic")
    b = _plate(3, _block((30, 40)))
    world = _world(a, o, b)
    _drive(o, _plate_centroid(a), 5.0)
    b.set_omega(o.omega.copy())
    for _ in range(3):
        _step(world)
    assert set(world.collision_evidence[1].neighbour.tolist()) == {2}
    assert set(world.collision_evidence[1].source.tolist()) == {cp.SOURCE_SLAB, cp.SOURCE_ARC}
    assert np.all(world.collision_evidence[1].role == cp.ROLE_UPPER)

    _remove_plate(world, 2)
    assert 2 not in world.collision_evidence
    _move(b, _column(30), _column(20))
    # The fallback alone would put A down: A drives into the front far faster than B.
    _drive(a, _plate_centroid(b), 10.0)
    _drive(b, _plate_centroid(a), 1.0)
    frame = _step(world)

    (record,) = _records(world, (1, 3))
    assert record.lower_plate_id == 3 and record.upper_plate_id == 1
    assert record.source == "slab" and not record.ambiguous
    assert record.votes["slab"] > 0 and record.votes["arc"] > 0
    assert frame.polarity[record.front_id] == (3, 1)
    assert np.any(frame.masks[3].lower) and not np.any(frame.masks[3].upper)
    assert np.any(frame.masks[1].upper) and not np.any(frame.masks[1].lower)


def test_paired_arcs_are_ambiguous_and_fall_back():
    # Ocean O subducts beneath continents on both sides (Molucca Sea): when they meet, each
    # carries upper-plate evidence, so the evidence can't decide.
    a = _plate(1, _block((10, 20)))
    o = _plate(2, _block((20, 30)), "oceanic")
    c = _plate(3, _block((30, 40)))
    world = _world(a, o, c)
    _drive(a, _plate_centroid(o), 5.0)
    _drive(c, _plate_centroid(o), 5.0)
    for _ in range(3):
        _step(world)
    assert np.all(world.collision_evidence[3].role == cp.ROLE_UPPER)

    _remove_plate(world, 2)
    _move(c, _column(30), _column(20))
    _drive(a, _plate_centroid(c), 6.0)
    _drive(c, _plate_centroid(a), 2.0)
    _step(world)

    (record,) = _records(world, (1, 3))
    assert record.ambiguous and record.source == "fallback"
    assert record.fallback_basis == "motion" and record.lower_plate_id == 1
    assert world.collision_polarity_stats["ambiguous_physical"] == 1
    assert world.collision_polarity_stats["fallback_ambiguous"] == 1


def _contact(*, b_speed: float = 3.0) -> tuple[World, PlateWithSparseQuadPatch, PlateWithSparseQuadPatch]:
    a = _plate(1, _block((10, 20)))
    b = _plate(2, _block((20, 30)))
    world = _world(a, b)
    _drive(b, _plate_centroid(a), b_speed)
    return world, a, b


def test_arc_alone_decides_when_no_physical_evidence():
    world, a, b = _contact()
    cp.add_evidence(world, b, _near(b, a), cp.SOURCE_ARC, cp.ROLE_UPPER, 99)
    _step(world)
    (record,) = _records(world, (1, 2))
    assert record.source == "arc" and record.lower_plate_id == 1


def test_arc_cue_never_overrules_consumption():
    world, a, b = _contact()
    near = _near(b, a)
    for _ in range(10):
        cp.add_evidence(world, b, near, cp.SOURCE_ARC, cp.ROLE_UPPER, 99)
    cp.add_evidence(world, b, near[:3], cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    _step(world)
    (record,) = _records(world, (1, 2))
    assert record.source == "consumption" and record.lower_plate_id == 2
    assert world.collision_polarity_stats["arc_overruled"] == 1


def test_contradictory_consumption_is_ambiguous():
    world, a, b = _contact()
    cp.add_evidence(world, a, _near(a, b), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    cp.add_evidence(world, b, _near(b, a), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 98)
    _step(world)
    (record,) = _records(world, (1, 2))
    assert record.ambiguous and record.source == "fallback"
    # Only B moves, into the front, so the motion fallback puts B down.
    assert record.fallback_basis == "motion" and record.lower_plate_id == 2


def _long_contact():
    j_range = (0, 50)
    a = _plate(1, _block((10, 20), j_range))
    b = _plate(2, _block((20, 30), j_range))
    world = _world(a, b)
    _drive(b, _plate_centroid(a), 3.0)
    return world, a, b


def _far_end(plate, other, j_lo: int):
    """`plate`'s nodes facing `other` at lattice columns j >= `j_lo`."""
    points = plate.all_points_and_elevation()[0]
    dist, _ = other.get_node_kdtree().query(points)
    _, _, _, j = unpack_cell_keys(plate.cell_keys)
    return points[(dist <= 1.5 * SPACING) & (j >= j_lo)]


def test_pair_evidence_from_elsewhere_on_the_boundary_decides():
    # The pair's subduction was recorded far along the boundary from where the continents
    # meet: B overrode A's ocean floor there. The motion fallback alone would put B down.
    world, a, b = _long_contact()
    cp.add_evidence(world, b, _far_end(b, a, 40), cp.SOURCE_CONSUMPTION, cp.ROLE_UPPER, 1)
    # Evidence naming a third plate, however close, doesn't count at pair scope.
    cp.add_evidence(world, a, _far_end(a, b, 40), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    _, _, _, j = unpack_cell_keys(b.cell_keys)
    b.remove_cells(j >= 15)
    _step(world)
    (record,) = _records(world, (1, 2))
    assert (record.scope, record.source, record.lower_plate_id) == ("pair", "consumption", 1)
    assert world.collision_polarity_stats["scope_pair"] == 1


def test_front_evidence_outranks_pair_evidence():
    world, a, b = _long_contact()
    _, _, _, j = unpack_cell_keys(b.cell_keys)
    b.remove_cells(j >= 15)
    cp.add_evidence(world, a, _near(a, b), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    cp.add_evidence(world, b, _far_end(_plate(2, _block((20, 30), (0, 50))), a, 40), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 1)
    _step(world)
    (record,) = _records(world, (1, 2))
    assert (record.scope, record.lower_plate_id) == ("front", 1)


def test_new_front_without_evidence_inherits_the_pairs_existing_polarity():
    world, a, b = _long_contact()
    _, _, _, j = unpack_cell_keys(b.cell_keys)
    b.remove_cells((j >= 15) & (j < 35))
    # Two fronts, decided by front evidence at the j < 15 end only.
    cp.add_evidence(world, a, _near(a, b)[:3], cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    first_only = _plate(2, _block((20, 30), (0, 15)))
    first_only.set_omega(b.omega.copy())
    world.plates = [a, first_only]
    _step(world)
    (established,) = _records(world, (1, 2))
    assert established.lower_plate_id == 1
    # The far stretch of B now touches A too: no evidence there, so it takes the pair's
    # established polarity rather than the motion fallback (which would put B down).
    world.plates = [a, b]
    _step(world)
    (inherited,) = [r for r in _records(world, (1, 2)) if r.front_id != established.front_id]
    assert (inherited.source, inherited.scope, inherited.lower_plate_id) == ("inherited", "record", 1)
    assert inherited.parent_id == established.front_id
    assert world.collision_polarity_stats["decided_fallback"] == 0


def test_fragment_beside_an_established_front_copies_it_before_pair_evidence():
    # Pair evidence far along the boundary says A is lower; the established front beside the
    # new fragment says B is. The fragment is part of that collision.
    world, a, full = _long_contact()
    _, _, _, j = unpack_cell_keys(full.cell_keys)

    def b_with(mask):
        plate = _plate(2, full.cell_keys[mask])
        plate.set_omega(full.omega.copy())
        return plate

    first = b_with(j < 20)
    world.plates = [a, first]
    cp.add_evidence(world, first, _near(first, a), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    cp.add_evidence(world, a, _far_end(a, full, 45), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 2)
    _step(world)
    (established,) = _records(world, (1, 2))
    assert established.lower_plate_id == 2

    # A separate stretch eight cells along: past the link radius, inside FRONT_INHERIT_SPACINGS,
    # and well short of the far-end evidence.
    world.plates = [a, b_with((j < 20) | ((j >= 28) & (j < 36)))]
    _step(world)
    (fragment,) = [r for r in _records(world, (1, 2)) if r.front_id != established.front_id]
    assert (fragment.source, fragment.lower_plate_id, fragment.parent_id) == ("inherited", 2, established.front_id)


def test_disagreeing_fronts_that_grow_together_keep_their_own_stretches():
    # Two fronts of one pair with opposite polarities (as in the nearby-fronts test), then
    # B's gap fills in and they become one connected contact.
    a = _plate(1, _block((10, 20), (10, 44)))
    prongs = np.concatenate([_block((20, 26), (12, 17)), _block((20, 26), (24, 29)), _block((26, 32), (12, 29))])
    b = _plate(2, prongs)
    world = _world(a, b)
    _drive(b, _plate_centroid(a), 3.0)
    prong_1 = _column(19, (12, 17))
    prong_2 = _column(19, (24, 29))
    near_b, near_a = _near(b, a), _near(a, b)
    cp.add_evidence(world, b, near_b[near_b @ prong_1 > near_b @ prong_2], cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    cp.add_evidence(world, a, near_a[near_a @ prong_2 > near_a @ prong_1], cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    _step(world)
    lowers = {r.front_id: r.lower_plate_id for r in _records(world, (1, 2))}
    assert sorted(lowers.values()) == [1, 2]

    filled = _plate(2, np.concatenate([prongs, _block((20, 26), (17, 24))]))
    filled.set_omega(b.omega.copy())
    world.plates = [a, filled]
    for _ in range(3):
        frame = _step(world)
        assert {r.front_id: r.lower_plate_id for r in _records(world, (1, 2))} == lowers
        assert set(frame.polarity) == set(lowers)
    assert world.collision_polarity_stats["fronts_merged"] == 0
    assert world.collision_polarity_stats["conflicts_kept_apart"] == 3
    # Each record holds the stretch it started on: B is lower along prong 1, upper along prong 2.
    masks = world.collision_polarity_frame.masks[2]
    _, _, _, j = unpack_cell_keys(filled.cell_keys)
    assert np.any(masks.lower & (j < 17)) and not np.any(masks.lower & (j >= 24))
    assert np.any(masks.upper & (j >= 24)) and not np.any(masks.upper & (j < 17))


def test_agreeing_fronts_that_grow_together_merge_into_the_oldest():
    a = _plate(1, _block((10, 20), (10, 44)))
    prong_1 = _plate(2, np.concatenate([_block((20, 26), (24, 29)), _block((26, 32), (12, 29))]))
    world = _world(a, prong_1)
    _drive(prong_1, _plate_centroid(a), 3.0)
    cp.add_evidence(world, a, _near(a, prong_1), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    _step(world)
    (oldest,) = _records(world, (1, 2))
    world.elapsed_years += 1_000_000.0
    both = _plate(2, np.concatenate([prong_1.cell_keys, _block((20, 26), (12, 17))]))
    both.set_omega(prong_1.omega.copy())
    world.plates = [a, both]
    cp.add_evidence(world, a, _near(a, both), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    _step(world)
    assert len(_records(world, (1, 2))) == 2
    filled = _plate(2, np.concatenate([both.cell_keys, _block((20, 26), (17, 24))]))
    filled.set_omega(prong_1.omega.copy())
    world.plates = [a, filled]
    _step(world)
    (merged,) = _records(world, (1, 2))
    assert merged.front_id == oldest.front_id
    assert world.collision_polarity_stats["fronts_merged"] == 1


def test_a_contact_returning_after_expiry_is_decided_afresh():
    world, a, b = _contact()
    _step(world)
    (old,) = world.collision_fronts
    away = _column(26)
    _move(b, _column(20), away)
    b.set_omega(np.zeros(3))
    world.elapsed_years += cp.FRONT_EXPIRY_YEARS + 1.0
    # Back onto exactly its old footprint, driving in again.
    _move(b, away, _column(20))
    _drive(b, _plate_centroid(a), 3.0)
    _step(world)
    (fresh,) = world.collision_fronts
    assert fresh.front_id != old.front_id and fresh.parent_id is None
    assert world.collision_polarity_stats["fronts_expired"] == 1


# --- Fallback --------------------------------------------------------------------------------


def test_fallback_is_motion_then_size_then_plate_id():
    world, _, _ = _contact()
    _step(world)
    (record,) = _records(world, (1, 2))
    assert (record.source, record.fallback_basis, record.lower_plate_id) == ("fallback", "motion", 2)
    assert world.collision_polarity_stats["fallback_no_evidence"] == 1

    # Stationary and overlapping by one column: tied on motion, so the smaller plate goes down.
    world = _world(_plate(1, _block((10, 21))), _plate(2, _block((20, 26))))
    _step(world)
    (record,) = _records(world, (1, 2))
    assert (record.fallback_basis, record.lower_plate_id) == ("size", 2)

    # Mirror-image plates about the face centre: equal areas, so the lower id goes down --
    # whichever list order the plates are in.
    left = (N // 2 - 8, N // 2 + 1)
    right = (N - left[1], N - left[0])
    for order in (1, -1):
        world = _world(*[_plate(1, _block(left)), _plate(2, _block(right))][::order])
        _step(world)
        (record,) = _records(world, (1, 2))
        assert (record.fallback_basis, record.lower_plate_id) == ("plate_id", 1)


# --- Fronts ----------------------------------------------------------------------------------


def test_nearby_fronts_of_one_pair_keep_separate_records_and_polarities():
    # B touches A with two prongs four cells apart: two fronts, each with its own history.
    a = _plate(1, _block((10, 20), (10, 44)))
    b = _plate(2, np.concatenate([_block((20, 26), (12, 17)), _block((20, 26), (24, 29)), _block((26, 32), (12, 29))]))
    world = _world(a, b)
    _drive(b, _plate_centroid(a), 3.0)
    prong_1 = _centroid(_plate(9, _block((19, 21), (12, 17))).all_points_and_elevation()[0])
    prong_2 = _centroid(_plate(9, _block((19, 21), (24, 29))).all_points_and_elevation()[0])
    near_b = _near(b, a)
    near_a = _near(a, b)
    cp.add_evidence(world, b, near_b[near_b @ prong_1 > near_b @ prong_2], cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    cp.add_evidence(world, a, near_a[near_a @ prong_2 > near_a @ prong_1], cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)

    _step(world)
    records = _records(world, (1, 2))
    assert len(records) == 2
    lowers = {r.front_id: r.lower_plate_id for r in records}
    assert sorted(lowers.values()) == [1, 2]
    for _ in range(3):
        _step(world)
    assert {r.front_id: r.lower_plate_id for r in _records(world, (1, 2))} == lowers
    assert world.collision_polarity_stats["fronts_created"] == 2


def test_migrating_front_keeps_its_record():
    world, a, b = _contact()
    _step(world)
    (record,) = _records(world, (1, 2))
    for k in range(3):
        # B advances one more column over A each step: the front moves across A.
        _move(b, _column(21 + k), _column(20 + k))
        _step(world)
    assert _records(world, (1, 2)) == [record]
    assert record.contact_steps == 4
    assert world.collision_polarity_stats["fronts_created"] == 1


def test_front_split_inherits_polarity_and_merge_folds_back():
    j_range = (10, 40)
    a = _plate(1, _block((10, 20), j_range))
    b = _plate(2, _block((20, 30), j_range))
    world = _world(a, b)
    _drive(b, _plate_centroid(a), 3.0)
    cp.add_evidence(world, a, _near(a, b), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    _step(world)
    (parent,) = _records(world, (1, 2))
    assert parent.lower_plate_id == 1

    # A six-cell gap opens in B's margin: two fronts.
    _, _, _, j = unpack_cell_keys(b.cell_keys)
    gap = (j >= 22) & (j < 28)
    b.remove_cells(gap)
    _step(world)
    records = _records(world, (1, 2))
    assert len(records) == 2
    (child,) = [r for r in records if r.front_id != parent.front_id]
    assert child.parent_id == parent.front_id
    assert child.lower_plate_id == parent.lower_plate_id == 1
    assert world.collision_polarity_stats["fronts_split"] == 1

    # The gap closes again: one front, one record.
    world.plates = [a, b2 := _plate(2, _block((20, 30), j_range))]
    b2.set_omega(b.omega.copy())
    _step(world)
    (merged,) = _records(world, (1, 2))
    assert merged.lower_plate_id == 1 and merged.front_id == parent.front_id
    assert world.collision_polarity_stats["fronts_merged"] == 1
    assert world.collision_polarity_stats["conflicts_kept_apart"] == 0


def test_front_expires_after_sustained_loss_of_contact():
    world, a, b = _contact()
    _step(world)
    assert len(world.collision_fronts) == 1
    _move(b, _column(20), _column(26))
    b.set_omega(np.zeros(3))
    elapsed = 0.0
    while elapsed <= cp.FRONT_EXPIRY_YEARS:
        world.elapsed_years += 1_000_000.0
        elapsed += 1_000_000.0
        cp.observe_contacts(world, STEP_YEARS)
    assert world.collision_fronts == []
    assert world.collision_polarity_stats["fronts_expired"] == 1


def test_evidence_expires_after_the_lookback():
    world, a, b = _contact()
    _move(b, _column(20), _column(26))
    b.set_omega(np.zeros(3))
    cp.add_evidence(world, a, a.all_points_and_elevation()[0][:5], cp.SOURCE_SLAB, cp.ROLE_UPPER, 99)
    world.elapsed_years += cp.EVIDENCE_LOOKBACK_YEARS - 1.0
    cp.observe_contacts(world, STEP_YEARS)
    assert 1 in world.collision_evidence
    world.elapsed_years += 2.0
    cp.observe_contacts(world, STEP_YEARS)
    assert 1 not in world.collision_evidence


# --- Topology --------------------------------------------------------------------------------


def _front_with_history(j_range=(10, 40)):
    a = _plate(1, _block((10, 20), j_range))
    b = _plate(2, _block((20, 30), j_range))
    world = _world(a, b)
    world.next_plate_id = 10
    _drive(b, _plate_centroid(a), 3.0)
    cp.add_evidence(world, b, _near(b, a), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    cp.add_evidence(world, a, _near(a, b), cp.SOURCE_CONSUMPTION, cp.ROLE_UPPER, 99)
    _step(world)
    return world, a, b


def test_plate_split_hands_the_front_to_the_daughter_carrying_the_contact():
    world, a, b = _front_with_history()
    (record,) = world.collision_fronts
    frames = cp.begin_topology(world)
    interior = _plate(1, _block((10, 15), (10, 40)))
    margin = _plate(10, _block((15, 20), (10, 40)))
    world.plates = [interior, margin, b]
    cp.note_lineage(world, 1, 1)
    cp.note_lineage(world, 1, 10)
    cp.end_topology(world, frames)

    assert world.collision_fronts == [record]
    assert record.plate_ids == (2, 10) and record.lower_plate_id == 2
    assert 1 not in world.collision_evidence
    assert np.all(world.collision_evidence[10].role == cp.ROLE_UPPER)
    _step(world)
    assert _records(world, (2, 10)) == [record]
    assert world.collision_polarity_stats["fronts_created"] == 1


def _front_with_lower_a(j_range=(0, 50)):
    """A long A-B front where evidence puts A down, though B drives into A -- so the motion
    fallback alone would decide the opposite."""
    a = _plate(1, _block((10, 20), j_range))
    b = _plate(2, _block((20, 30), j_range))
    world = _world(a, b)
    world.next_plate_id = 10
    _drive(b, _plate_centroid(a), 3.0)
    cp.add_evidence(world, a, _near(a, b), cp.SOURCE_CONSUMPTION, cp.ROLE_LOWER, 99)
    _step(world)
    (record,) = world.collision_fronts
    assert record.lower_plate_id == 1
    return world, a, b, record


def _split_a(world: World, pieces: dict[int, tuple[int, int]], j_range=(0, 50)) -> None:
    """Replace plate 1 by daughters, each a j-range of its cells, all descending from 1."""
    frames = cp.begin_topology(world)
    b = next(p for p in world.plates if p.plate_id == 2)
    daughters = [_plate(pid, _block((10, 20), span)) for pid, span in pieces.items()]
    world.plates = [*daughters, b]
    for pid in pieces:
        cp.note_lineage(world, 1, pid)
    cp.end_topology(world, frames)


def test_plate_split_across_a_front_keeps_the_decision_on_every_daughter():
    world, a, b, record = _front_with_lower_a()
    _split_a(world, {1: (0, 25), 10: (25, 50)})

    lowers = {r.plate_ids: r.lower_plate_id for r in world.collision_fronts}
    assert lowers == {(1, 2): 1, (2, 10): 10}
    (child,) = [r for r in world.collision_fronts if r.front_id != record.front_id]
    assert child.parent_id == record.front_id and child.source == record.source
    assert world.collision_polarity_stats["fronts_split_topology"] == 1

    # Both stretches match their records next step: nothing is decided afresh.
    frame = _step(world)
    assert world.collision_polarity_stats["fronts_created"] == 1
    assert sorted(frame.polarity.values()) == [(1, 2), (10, 2)]


def test_narrow_stretch_left_on_a_daughter_keeps_the_decision():
    # A long front split so one daughter holds a one-column stretch of the contact: two front
    # nodes on each side, a valid front (FRONT_MIN_NODES counts both sides). Through a
    # 64-point sample of the front (over 120 nodes a side) it would hold fewer anchors than
    # that and be dropped.
    world, a, b, record = _front_with_lower_a()
    assert len(record.side_points[1]) > 64
    _split_a(world, {1: (0, 30), 10: (30, 31), 11: (31, 50)})

    lowers = {r.plate_ids: r.lower_plate_id for r in world.collision_fronts}
    assert lowers == {(1, 2): 1, (2, 10): 10, (2, 11): 11}
    narrow = next(r for r in world.collision_fronts if r.plate_ids == (2, 10))
    assert {pid: len(points) for pid, points in narrow.side_points.items()} == {10: 2, 2: 2}
    _step(world)
    assert world.collision_polarity_stats["fronts_created"] == 1
    assert {r.plate_ids: r.lower_plate_id for r in world.collision_fronts} == lowers


def test_a_small_front_survives_its_plate_being_rehomed():
    # A three-node front -- the minimum, split two nodes to one across the contact. Absorbing
    # one of its plates into another is re-homing, not a split, so the record must survive
    # however few nodes it has on that side.
    world, a, b, record = _front_with_lower_a()
    record.side_points = {pid: points[:2] if pid == 1 else points[:1] for pid, points in record.side_points.items()}
    frames = cp.begin_topology(world)
    world.plates = [_plate(5, a.cell_keys), b]
    cp.note_lineage(world, 1, 5)
    cp.end_topology(world, frames)
    assert world.collision_fronts == [record]
    assert record.plate_ids == (2, 5) and record.lower_plate_id == 5
    assert world.collision_polarity_stats["fronts_dropped_topology"] == 0


def test_merging_the_pair_drops_the_front_and_merging_a_third_plate_rehomes_it():
    world, a, b = _front_with_history()
    frames = cp.begin_topology(world)
    world.plates = [_plate(1, np.concatenate([a.cell_keys, b.cell_keys]))]
    cp.note_lineage(world, 2, 1)
    cp.end_topology(world, frames)
    assert world.collision_fronts == []
    assert world.collision_polarity_stats["fronts_fused"] == 1
    assert np.any(world.collision_evidence[1].role == cp.ROLE_LOWER)

    world, a, b = _front_with_history()
    (record,) = world.collision_fronts
    frames = cp.begin_topology(world)
    world.plates = [_plate(5, a.cell_keys), b]
    cp.note_lineage(world, 1, 5)
    cp.end_topology(world, frames)
    assert record.plate_ids == (2, 5) and record.lower_plate_id == 2


def test_vanished_plate_takes_its_fronts_and_evidence():
    world, a, b = _front_with_history()
    _remove_plate(world, 2)
    assert world.collision_fronts == []
    assert 2 not in world.collision_evidence
    assert world.collision_polarity_stats["fronts_dropped_topology"] == 1


# --- Persistence, determinism, cost ----------------------------------------------------------


def test_records_and_evidence_survive_save_and_load():
    world, a, b = _front_with_history()
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    (before,) = world.collision_fronts
    (after,) = loaded.collision_fronts
    assert (after.front_id, after.plate_ids, after.lower_plate_id, after.source) == (
        before.front_id, before.plate_ids, before.lower_plate_id, before.source
    )
    for pid, store in world.collision_evidence.items():
        np.testing.assert_array_equal(loaded.collision_evidence[pid].keys, store.keys)
        np.testing.assert_array_equal(loaded.collision_evidence[pid].count, store.count)
    _step(loaded)
    assert loaded.collision_fronts[0].front_id == before.front_id
    assert loaded.collision_polarity_stats["fronts_created"] == 1


def test_old_save_without_polarity_state_loads_empty():
    world, _, _ = _contact()
    for name in ("collision_evidence", "collision_fronts", "collision_polarity_stats"):
        del world.__dict__[name]
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.collision_evidence == {} and loaded.collision_fronts == [] and loaded.collision_polarity_stats == {}
    _step(loaded)
    assert len(loaded.collision_fronts) == 1


def _busy_world() -> World:
    """Several contacts at once: a continental front with history, an active margin, and a
    stationary overlap the fallback decides."""
    a = _plate(1, _block((10, 20), (10, 30)))
    b = _plate(2, _block((20, 30), (10, 30)))
    o = _plate(3, _block((30, 40), (10, 30)), "oceanic")
    c = _plate(4, _block((10, 20), (29, 36)))
    world = _world(a, b, o, c)
    _drive(b, _plate_centroid(a), 3.0)
    _drive(o, _plate_centroid(b), 4.0)
    for _ in range(2):
        _step(world)
    return world


def _frame_signature(world: World):
    frame = world.collision_polarity_frame
    masks = {
        pid: tuple(np.flatnonzero(getattr(m, name)).tolist() for name in ("lower", "upper", "retreat_eligible", "override"))
        + (m.front_id.tolist(),)
        for pid, m in frame.masks.items()
    }
    records = sorted((r.front_id, r.plate_ids, r.lower_plate_id, r.source) for r in world.collision_fronts)
    evidence = {pid: (s.keys.tolist(), s.source.tolist(), s.role.tolist(), s.neighbour.tolist(), s.count.tolist()) for pid, s in world.collision_evidence.items()}
    return frame.polarity, masks, records, evidence


def test_plate_order_and_deform_order_do_not_change_polarity_or_masks():
    world = _busy_world()
    reversed_world = copy.deepcopy(world)
    reversed_world.plates = reversed_world.plates[::-1]

    cp.observe_contacts(world, STEP_YEARS)
    cp.observe_contacts(reversed_world, STEP_YEARS)
    signature = _frame_signature(world)
    assert signature == _frame_signature(reversed_world)
    assert len(signature[0]) >= 2

    # The decisions are frozen before deform: deforming in opposite orders leaves them alone.
    for w, order in ((world, 1), (reversed_world, -1)):
        for plate in w.plates[::order]:
            others = [p for p in w.plates if p.plate_id != plate.plate_id]
            plate.deform(w, others, STEP_YEARS, SPACING)
    assert _frame_signature(world)[:3] == signature[:3]
    assert _frame_signature(reversed_world)[:3] == signature[:3]


def test_search_cache_answers_unchanged_neighbours_exactly():
    world = _busy_world()
    a, b = world.plates[0], world.plates[1]

    def context(cache):
        world.boundary_search_cache = cache
        try:
            return boundary_context(world, a, world.plates[1:], STEP_YEARS, lambda c: quad_tectonics.components_of_at_least(a, c, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN))
        finally:
            world.boundary_search_cache = None

    def same(x, y):
        for name in ("convergent", "divergent", "transform", "contested", "closing_rate"):
            np.testing.assert_array_equal(getattr(x, name), getattr(y, name))
        for name in ("dist_to_neighbor", "direction_to_neighbor", "neighbor_plate_id", "neighbor_node_index"):
            np.testing.assert_array_equal(getattr(x.inputs, name), getattr(y.inputs, name))

    cache = torque.BoundarySearchCache()
    fresh = context(None)
    same(context(cache), fresh)
    searches = cache.searches
    same(context(cache), fresh)
    assert cache.searches == searches and cache.reused > 0

    # A neighbour that changed is searched again.
    b.remove_cells(np.arange(b.node_count()) == 0)
    fresh = context(None)
    same(context(cache), fresh)
    assert cache.searches > searches

    # The new neighbour identity: every node's nearest other-plate node, by plate and index.
    inputs = fresh.inputs
    ids = {p.plate_id: p for p in world.plates}
    for k in np.flatnonzero(np.isfinite(inputs.dist_to_neighbor))[:20]:
        neighbour = ids[int(inputs.neighbor_plate_id[k])]
        point = neighbour.all_points_and_elevation()[0][inputs.neighbor_node_index[k]]
        best = min(float(np.min(np.linalg.norm(p.all_points_and_elevation()[0] - inputs.own_points[k], axis=1))) for p in world.plates[1:])
        assert np.isclose(np.linalg.norm(point - inputs.own_points[k]), best)


def test_retreat_records_consumed_ocean_floor_on_both_sides():
    a = _plate(1, _block((10, 21)))
    o = _plate(2, _block((20, 30)), "oceanic")
    world = _world(a, o)
    quad_tectonics.deform(o, world, [a], STEP_YEARS, SPACING)

    lower = world.collision_evidence[2]
    upper = world.collision_evidence[1]
    assert set(lower.source.tolist()) == {cp.SOURCE_CONSUMPTION} and np.all(lower.role == cp.ROLE_LOWER)
    assert set(lower.neighbour.tolist()) == {1}
    assert set(upper.source.tolist()) == {cp.SOURCE_CONSUMPTION} and np.all(upper.role == cp.ROLE_UPPER)
    assert set(upper.neighbour.tolist()) == {2}


def test_step_world_runs_the_prepass_and_reports_its_cost():
    from app import world as world_mod

    world = _busy_world()
    world_mod.step_world(world, STEP_YEARS)
    report = cp.summary(world)
    assert world.boundary_search_cache is None
    assert report["counters"]["prepass_calls"] == 3
    assert report["counters"]["deform_searches_reused"] > 0
    assert report["prepass_seconds_per_step"] > 0.0
    assert world.collision_polarity_frame is not None

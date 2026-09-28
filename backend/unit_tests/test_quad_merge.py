"""Cross-plate merge on sparse quad plates (issue #228 Phase 4): the frame remap of an absorbed
plate onto the surviving plate's lattice, checked by exact cell area -- see quad_merge.py."""

import copy

import numpy as np
import pytest

from app import geometry, lithosphere, merge_split, quad_merge
from app.elevation_lines import CRUST_TYPE_INHERIT, CRUST_TYPE_OCEANIC, line_spacing_rad
from app.lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M, new_plate
from app.sparse_quad_patch import PLANET_RADIUS_M, PlateWithSparseQuadPatch, unpack_cell_keys
from app.world import World

DENSITY = 0.5
SPACING = line_spacing_rad(DENSITY)
RADIUS = 8 * SPACING

# An arbitrary, fixed rotation, so the absorbed plate's lattice lines up with the survivor's
# nowhere -- the case a frame remap exists for.
_Q, _ = np.linalg.qr(np.random.default_rng(3).normal(size=(3, 3)))
ROTATED = _Q * np.sign(np.linalg.det(_Q))


def _direction(angle: float) -> np.ndarray:
    """World unit vector `angle` radians from +x toward +y."""
    return np.array([np.cos(angle), np.sin(angle), 0.0])


def _cap(plate_id, frame, centre, radius=RADIUS, crust_type="continental") -> PlateWithSparseQuadPatch:
    plate = PlateWithSparseQuadPatch.from_lattice(plate_id, frame, crust_type, SPACING, lambda pts: pts @ centre > np.cos(radius))
    world = plate.all_points_and_elevation()[0]
    hc, hm = lithosphere.reference_thickness(crust_type)
    plate.set_fields_on_plate(
        crustal_thickness_m=hc + 4000.0 * world[:, 1],
        mantle_lithosphere_thickness_m=np.full(len(world), hm),
    )
    lithosphere.sync_plate_elevation(plate)
    return plate


def _volume(plates, name) -> float:
    return float(sum(np.sum(p.collect(name) * p.node_areas_m2()) for p in plates))


def _area(plates) -> float:
    return float(sum(np.sum(p.node_areas_m2()) for p in plates))


def _world(*plates) -> World:
    return World(seed=0, plates=list(plates), next_plate_id=max(p.plate_id for p in plates) + 1, node_density=DENSITY)


def test_subsamples_partition_each_cell_by_exact_area():
    plate = _cap(1, ROTATED, _direction(0.0))

    source, local, area = quad_merge.subsample_cells(plate, 4)

    assert np.allclose(np.bincount(source, weights=area), plate.node_areas_m2(), rtol=1e-14, atol=0.0)
    # Every sub-cell centre lies inside the cell it came from.
    assert np.array_equal(plate._leaf_key_at(local), plate.cell_keys[source])


def test_merge_across_rotated_frames_conserves_every_extensive_field_exactly():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    rng = np.random.default_rng(0)
    for plate in (keep, absorb):
        plate.set_fields_on_plate(
            silt_depth=rng.uniform(0.0, 5.0, plate.node_count()),
            coal_deposit_m=rng.uniform(0.0, 2.0, plate.node_count()),
        )
    names = ("crustal_thickness_m", "mantle_lithosphere_thickness_m", "silt_depth", "coal_deposit_m")
    before = {name: _volume([keep, absorb], name) for name in names}
    area_before = _area([keep, absorb])

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    for name in names:
        assert _volume([keep], name) == pytest.approx(before[name], rel=1e-12), name
    # The footprint comes out within a boundary ring of the two plates' combined area.
    ring = 2 * np.pi * RADIUS * SPACING * PLANET_RADIUS_M**2
    assert abs(_area([keep]) - area_before) < 0.5 * ring
    keep._validate_leaf_topology()


def test_merge_keeps_the_survivors_cells_ids_and_values():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    before = dict(zip(map(int, keep.cell_keys), keep.collect("elevation")))

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    after = dict(zip(map(int, keep.cell_keys), keep.collect("elevation")))
    assert set(before) <= set(after)
    untouched = [key for key in before if after[key] == before[key]]
    # Only the survivor's cells along the seam can take in stray sub-cells.
    assert len(untouched) > 0.9 * len(before)
    assert keep.frame is not None and np.array_equal(keep.frame, np.eye(3))


def test_remapped_territory_covers_the_absorbed_plate():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    absorbed_points = absorb.all_points_and_elevation()[0].copy()

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    assert keep.contains_batch(absorbed_points).mean() > 0.97


def test_frame_remap_puts_values_where_they_were_in_world_space():
    # A field that is a smooth function of world position must come out as the same function
    # after the remap -- the check that the absorbed cells land in the right place.
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    absorb.set_fields_on_plate(divergent_age_myr=100.0 * absorb.all_points_and_elevation()[0][:, 1])
    old_keys = set(map(int, keep.cell_keys))

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    new = np.array([int(k) not in old_keys for k in keep.cell_keys])
    world = keep.all_points_and_elevation()[0][new]
    error = keep.collect("divergent_age_myr")[new] - 100.0 * world[:, 1]
    # Interior cells average a linear field exactly up to sub-cell resolution; boundary
    # cells, which are only partly covered, are off by at most about a cell.
    assert np.median(np.abs(error)) < 100.0 * 0.2 * SPACING
    assert np.max(np.abs(error)) < 100.0 * 1.5 * SPACING


def test_partly_covered_rim_cells_are_neither_thin_nor_thick():
    # Remapping volume cell by cell leaves the rim of the absorbed plate as a ring of thin
    # (partly covered) and thick (fed stray sub-cells) columns, i.e. false relief. The remap
    # instead keeps every new cell within the single conserving ratio of the source.
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    old_keys = set(map(int, keep.cell_keys))

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    new = np.array([int(k) not in old_keys for k in keep.cell_keys])
    world = keep.all_points_and_elevation()[0][new]
    expected = lithosphere.reference_thickness("continental")[0] + 4000.0 * world[:, 1]
    assert np.max(np.abs(keep.collect("crustal_thickness_m")[new] / expected - 1.0)) < 0.03


def test_suture_overlap_stacks_crust_and_conserves_it_below_the_cap():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS - 3 * SPACING))
    overlap = keep.contains_batch(absorb.all_points_and_elevation()[0])
    assert overlap.sum() > 20
    before = _volume([keep, absorb], "crustal_thickness_m")
    keep_hc = dict(zip(map(int, keep.cell_keys), keep.collect("crustal_thickness_m")))

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    assert _volume([keep], "crustal_thickness_m") == pytest.approx(before, rel=1e-12)
    hc = dict(zip(map(int, keep.cell_keys), keep.collect("crustal_thickness_m")))
    thickened = [key for key in keep_hc if hc[key] > keep_hc[key] + 1000.0]
    assert len(thickened) > 20
    assert max(hc.values()) <= SUTURE_ACCRETION_MAX_HC_M


def test_suture_cap_is_the_only_volume_that_leaves(monkeypatch):
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS - 3 * SPACING))
    keep.set_fields_on_plate(crustal_thickness_m=np.full(keep.node_count(), 0.8 * SUTURE_ACCRETION_MAX_HC_M))
    before = _volume([keep, absorb], "crustal_thickness_m")
    old = dict(zip(map(int, keep.cell_keys), keep.collect("crustal_thickness_m")))

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    hc = keep.collect("crustal_thickness_m")
    capped = np.isclose(hc, SUTURE_ACCRETION_MAX_HC_M)
    assert np.any(capped)
    assert np.all([int(k) in old for k in keep.cell_keys[capped]])
    lost = before - _volume([keep], "crustal_thickness_m")
    assert lost > 0.0
    # Exactly the stacked volume above the cap: the uncapped remap of this pair minus what
    # the capped columns hold.
    keep2 = _cap(1, np.eye(3), _direction(0.0))
    keep2.set_fields_on_plate(crustal_thickness_m=np.full(keep2.node_count(), 0.8 * SUTURE_ACCRETION_MAX_HC_M))
    absorb2 = _cap(2, ROTATED, _direction(2 * RADIUS - 3 * SPACING))
    cap = SUTURE_ACCRETION_MAX_HC_M
    monkeypatch.setattr(quad_merge, "SUTURE_ACCRETION_MAX_HC_M", np.inf)
    quad_merge.merge(keep2, absorb2, np.zeros((0, 3)))
    uncapped = keep2.collect("crustal_thickness_m")
    assert np.array_equal(keep2.cell_keys, keep.cell_keys)
    excess = np.sum(np.maximum(uncapped - cap, 0.0) * keep.node_areas_m2())
    assert lost == pytest.approx(excess, rel=1e-9)


def test_merge_does_not_grow_over_a_third_plate_and_still_conserves():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    # A third plate sitting on part of the absorbed plate's territory.
    third = _cap(3, np.eye(3), _direction(3 * RADIUS), radius=0.6 * RADIUS, crust_type="oceanic")
    assert np.any(third.contains_batch(absorb.all_points_and_elevation()[0]))
    before = _volume([keep, absorb], "crustal_thickness_m")

    quad_merge.merge(keep, absorb, third.all_points_and_elevation()[0])

    assert _volume([keep], "crustal_thickness_m") == pytest.approx(before, rel=1e-12)
    # No new cell sits deep inside the third plate: at most a boundary sliver shared by the
    # Voronoi split.
    deep = third.contains_batch(keep.all_points_and_elevation()[0])
    shrunk = _cap(4, np.eye(3), _direction(3 * RADIUS), radius=0.6 * RADIUS - 1.5 * SPACING, crust_type="oceanic")
    assert not np.any(shrunk.contains_batch(keep.all_points_and_elevation()[0]))
    assert deep.sum() < 0.1 * third.node_count()


def test_merge_carries_composition_provenance_and_history():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    world = absorb.all_points_and_elevation()[0]
    oceanic = world[:, 2] > 0.0
    codes = np.where(oceanic, CRUST_TYPE_OCEANIC, CRUST_TYPE_INHERIT).astype(np.int8)
    volcano = np.zeros(absorb.node_count(), dtype=bool)
    volcano[np.argmax(world[:, 1])] = True
    absorb.set_fields_on_plate(
        crust_type_code=codes,
        is_volcano=volcano,
        volcano_active_years_remaining=np.where(volcano, 5e5, 0.0),
        node_created_years=np.where(oceanic, 2e6, 1e6),
    )
    old_keys = set(map(int, keep.cell_keys))

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    new = np.array([int(k) not in old_keys for k in keep.cell_keys])
    world_new = keep.all_points_and_elevation()[0][new]
    new_codes = keep.collect("crust_type_code")[new]
    clear = np.abs(world_new[:, 2]) > 1.5 * SPACING
    assert np.all(new_codes[clear & (world_new[:, 2] > 0)] == CRUST_TYPE_OCEANIC)
    # Continental cells stay inherited on a continental survivor.
    assert np.all(new_codes[clear & (world_new[:, 2] < 0)] == CRUST_TYPE_INHERIT)
    assert keep.collect("is_volcano")[new].sum() >= 1
    assert keep.collect("volcano_active_years_remaining")[new].max() == 5e5
    created = keep.collect("node_created_years")[new]
    assert set(np.unique(created)) <= {1e6, 2e6}


def test_merge_of_an_oceanic_plate_freezes_its_inherited_composition():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING), crust_type="oceanic")
    old_keys = set(map(int, keep.cell_keys))

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    new = np.array([int(k) not in old_keys for k in keep.cell_keys])
    assert np.all(keep.collect("crust_type_code")[new] == CRUST_TYPE_OCEANIC)


def test_merged_elevation_is_isostatic_plus_the_carried_residual():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    absorb.set_fields_on_plate(elevation=absorb.collect("elevation") - 250.0)
    old_keys = set(map(int, keep.cell_keys))

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    new = np.array([int(k) not in old_keys for k in keep.cell_keys])
    density = lithosphere.node_crust_density(keep.collect("crust_type_code"), keep.crust_type)
    isostatic = lithosphere.isostatic_elevation(keep.collect("crustal_thickness_m"), keep.collect("mantle_lithosphere_thickness_m"), density)
    assert np.allclose((keep.collect("elevation") - isostatic)[new], -250.0)


def _momentum(*plates) -> np.ndarray:
    return sum(quad_merge._inertia(p) @ p.omega for p in plates)


def test_merge_conserves_angular_momentum_against_the_merged_plates_own_inertia():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    keep.set_omega(np.array([0.0, 0.0, 1e-9]))
    absorb.set_omega(np.array([2e-9, 0.0, 0.0]))
    momentum = _momentum(keep, absorb)

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    assert np.allclose(quad_merge._inertia(keep) @ keep.omega, momentum, rtol=1e-12, atol=0.0)
    assert keep.age_steps == 0


def test_crust_removed_at_the_suture_cap_leaves_with_the_absorbed_plates_momentum():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS - 3 * SPACING))
    keep.set_fields_on_plate(crustal_thickness_m=np.full(keep.node_count(), 0.95 * SUTURE_ACCRETION_MAX_HC_M))
    keep.set_omega(np.array([0.0, 0.0, 1e-9]))
    absorb.set_omega(np.array([2e-9, 0.0, 0.0]))
    momentum = _momentum(keep, absorb)
    probe_keep, probe_absorb = copy.deepcopy(keep), copy.deepcopy(absorb)
    lost = quad_merge._transfer(probe_keep, probe_absorb, np.zeros((0, 3)))
    assert np.linalg.norm(lost) > 0.0

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    retained = momentum - lost @ absorb.omega
    assert np.allclose(quad_merge._inertia(keep) @ keep.omega, retained, rtol=1e-12, atol=0.0)
    assert not np.allclose(retained, momentum, rtol=1e-6, atol=0.0)


def _overlapping_pair():
    """A survivor and an absorbed plate on the same lattice, the absorbed one wholly on top
    of the survivor: every absorbed cell lands on exactly one survivor cell (the suture)."""
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, np.eye(3), _direction(0.0), radius=0.5 * RADIUS)
    hc, _ = lithosphere.reference_thickness("continental")
    keep.set_fields_on_plate(crustal_thickness_m=np.full(keep.node_count(), 0.5 * hc))
    absorb.set_fields_on_plate(crustal_thickness_m=np.full(absorb.node_count(), 0.5 * hc))
    lithosphere.sync_plate_elevation(keep)
    lithosphere.sync_plate_elevation(absorb)
    under = keep.contains_batch(absorb.all_points_and_elevation()[0])
    assert under.all()
    stacked = np.isin(keep.cell_keys, absorb.cell_keys)
    return keep, absorb, stacked


def test_a_wholly_overlapping_merge_keeps_the_absorbed_volcanoes():
    keep, absorb, stacked = _overlapping_pair()
    absorb.set_fields_on_plate(
        is_volcano=np.ones(absorb.node_count(), dtype=bool),
        volcano_active_years_remaining=np.full(absorb.node_count(), 5e5),
    )
    count = keep.node_count()

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    assert keep.node_count() == count
    assert np.all(keep.collect("is_volcano")[stacked])
    assert np.all(keep.collect("volcano_active_years_remaining")[stacked] == 5e5)
    assert not np.any(keep.collect("is_volcano")[~stacked])


def test_suture_cells_combine_both_plates_by_remap_class():
    keep, absorb, stacked = _overlapping_pair()
    keep.set_fields_on_plate(
        divergent_age_myr=np.full(keep.node_count(), 10.0),
        node_created_years=np.full(keep.node_count(), 3e6),
        overlap_onset_years=np.zeros(keep.node_count()),
        elev_change_reason=np.full(keep.node_count(), 3.0),
    )
    absorb.set_fields_on_plate(
        divergent_age_myr=np.full(absorb.node_count(), 30.0),
        node_created_years=np.full(absorb.node_count(), 1e6),
        overlap_onset_years=np.full(absorb.node_count(), 2e6),
        elev_change_reason=np.full(absorb.node_count(), 5.0),
        # Thinner explicitly oceanic crust: outvoted by the survivor's by volume.
        crust_type_code=np.full(absorb.node_count(), CRUST_TYPE_OCEANIC, dtype=np.int8),
        crustal_thickness_m=np.full(absorb.node_count(), 7000.0),
    )

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    # Intensive: area-weighted mean of two equal-area contributions.
    assert np.allclose(keep.collect("divergent_age_myr")[stacked], 20.0)
    assert np.allclose(keep.collect("divergent_age_myr")[~stacked], 10.0)
    # History: the earliest valid time, ignoring the sentinel.
    assert np.all(keep.collect("node_created_years")[stacked] == 1e6)
    assert np.all(keep.collect("overlap_onset_years")[stacked] == 2e6)
    # Categorical: an equal-area tie keeps the survivor's; composition votes by volume.
    assert np.all(keep.collect("elev_change_reason")[stacked] == 3.0)
    assert np.all(keep.collect("crust_type_code")[stacked] == CRUST_TYPE_INHERIT)


def test_suture_elevation_blends_both_residuals():
    keep, absorb, stacked = _overlapping_pair()
    absorb.set_fields_on_plate(elevation=absorb.collect("elevation") - 200.0)

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    density = lithosphere.node_crust_density(keep.collect("crust_type_code"), keep.crust_type)
    isostatic = lithosphere.isostatic_elevation(keep.collect("crustal_thickness_m"), keep.collect("mantle_lithosphere_thickness_m"), density)
    assert np.allclose((keep.collect("elevation") - isostatic)[stacked], -100.0)


def test_quad_pairs_are_offered_to_merge_and_fuse_through_merge_plates():
    keep = _cap(1, np.eye(3), _direction(0.0))
    absorb = _cap(2, ROTATED, _direction(2 * RADIUS + 0.5 * SPACING))
    bystander = _cap(3, np.eye(3), _direction(-2.5 * RADIUS), crust_type="oceanic")
    world = _world(keep, absorb, bystander)
    world.debug_diagnostics = True
    before = _volume([keep, absorb], "crustal_thickness_m")

    assert merge_split._supports_merge(world, 1, 2)
    merge_split.merge_plates(world, 1, 2)

    assert [p.plate_id for p in world.plates] == [1, 3]
    assert _volume([keep], "crustal_thickness_m") == pytest.approx(before, rel=1e-12)
    # The phase budget books the merge by exact area: no spurious gain or loss from the
    # remap onto differently sized cells.
    scope = world.phase_budget["plate_merge"]["scopes"]["all"]
    assert scope["sum_hc_after"] == pytest.approx(scope["sum_hc_before"], rel=1e-12)


def test_a_line_plate_and_a_quad_plate_are_not_offered_to_merge():
    quad = _cap(1, np.eye(3), _direction(0.0))
    line = new_plate(2, np.eye(3), "continental", SPACING, seed=0, is_owned=lambda pts: pts @ _direction(2 * RADIUS) > np.cos(RADIUS))
    world = _world(quad, line)
    world.overlap_progress[(1, 2)] = merge_split.FORCED_MERGE_SUSTAINED_YEARS * 2

    assert not merge_split._supports_merge(world, 1, 2)
    forced = merge_split.pop_ready_forced_merge(world, can_merge=lambda x, y: merge_split._supports_merge(world, x, y))
    assert forced is None
    assert world.overlap_progress[(1, 2)] == merge_split.FORCED_MERGE_SUSTAINED_YEARS * 2


def test_merge_across_a_cube_face_seam_of_the_survivors_lattice():
    # The absorbed plate sits 45+ degrees from the survivor's seed, so it lands on the
    # survivor's lattice across a cube-face edge.
    keep = _cap(1, np.eye(3), geometry.normalize(np.array([1.0, 0.9, 0.0])))
    absorb = _cap(2, ROTATED, geometry.normalize(geometry.normalize(np.array([1.0, 0.9, 0.0])) + np.array([0.0, 0.0, 2.05 * RADIUS])))
    before = _volume([keep, absorb], "crustal_thickness_m")

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    assert _volume([keep], "crustal_thickness_m") == pytest.approx(before, rel=1e-12)
    faces = np.unique(unpack_cell_keys(keep.cell_keys)[0])
    assert len(faces) > 1
    keep._validate_leaf_topology()

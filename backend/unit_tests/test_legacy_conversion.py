"""Issue #248: one-way conversion of line-backed worlds to sparse quads (legacy_conversion.py)
and the compatibility policy in docs/save-compatibility.md."""

import copy
import pickle

import numpy as np
import pytest
from app import geometry, legacy_conversion, persistence, surface_parity
from app.elevation_lines import line_spacing_rad
from app.sparse_quad_patch import PlateWithSparseQuadPatch
from app.surface_fields import SURFACE_FIELDS, RemapClass
from app.world import World, generate_world, step_world

from .legacy_lines import legacy_line_plate, line_world_like, retired_class_names

DENSITY = 1.0
SPACING = line_spacing_rad(DENSITY)


def _xyz(lat_deg: float, lon_deg: float) -> np.ndarray:
    return geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))


def _angle(points: np.ndarray, centre: np.ndarray) -> np.ndarray:
    return np.arccos(np.clip(points @ centre, -1.0, 1.0))


def _in_slot(local: np.ndarray) -> np.ndarray:
    return (np.abs(local[:, 2]) < 0.06) & (local[:, 1] > 0.03)


def _shape_world() -> World:
    """Five isolated plates, each one awkward for a row-based or chart-based representation:
    an annulus (a hole), a C (a concave notch), a disc on the antimeridian, a polar cap, and
    one plate made of two disconnected discs. Everything else is open sphere."""
    annulus = _xyz(30, 0)
    c_shape = _xyz(-30, 90)
    c_frame = geometry.plate_frame_from_seed(c_shape)
    antimeridian = _xyz(0, 180)
    pole = _xyz(90, 0)
    fragment_a, fragment_b = _xyz(-50, -110), _xyz(-50, -60)
    shapes = {
        0: (annulus, "continental", lambda p: (_angle(p, annulus) > 0.12) & (_angle(p, annulus) < 0.30)),
        # A straight-sided slot six spacings wide, cut from the centre out through the rim.
        1: (c_shape, "oceanic", lambda p: (_angle(p, c_shape) < 0.25) & ~_in_slot(geometry.to_local(c_frame, p))),
        2: (antimeridian, "continental", lambda p: _angle(p, antimeridian) < 0.2),
        3: (pole, "oceanic", lambda p: _angle(p, pole) < 0.25),
        4: (fragment_a, "continental", lambda p: (_angle(p, fragment_a) < 0.12) | (_angle(p, fragment_b) < 0.12)),
    }
    frames = {1: c_frame}
    plates = [
        legacy_line_plate(plate_id, frames.get(plate_id, geometry.plate_frame_from_seed(centre)), crust_type, SPACING, owned)
        for plate_id, (centre, crust_type, owned) in shapes.items()
    ]
    world = World(seed=11, plates=plates, next_plate_id=len(plates), node_density=DENSITY)
    from app.hydroclimate import eustasy

    eustasy.initialize_water_budget(world)
    return world


@pytest.fixture(scope="module")
def shape_world() -> World:
    return _shape_world()


@pytest.fixture(scope="module")
def converted_shapes(shape_world):
    world = copy.deepcopy(shape_world)
    report = legacy_conversion.convert_world_to_quad(world)
    return world, report


def _line_nodes(world: World) -> tuple[np.ndarray, np.ndarray]:
    points, owners = [], []
    for plate in world.plates:
        xyz = plate.all_points_and_elevation()[0]
        points.append(xyz)
        owners.append(np.full(len(xyz), plate.plate_id))
    return np.concatenate(points), np.concatenate(owners)


def test_territory_keeps_every_shape_and_every_void(shape_world, converted_shapes):
    world, _ = converted_shapes
    sample = surface_parity.fibonacci_sphere(60_000)
    nodes, owners = _line_nodes(shape_world)
    from scipy.spatial import cKDTree

    dist, nearest = cKDTree(nodes).query(sample)
    dist_rad = 2.0 * np.arcsin(np.clip(dist / 2.0, 0.0, 1.0))
    claims = np.stack([plate.contains_batch(sample) for plate in world.plates])
    ids = np.array([plate.plate_id for plate in world.plates])

    # The annulus' hole, the C's notch and all the open sphere stay open.
    far = dist_rad > 2.5 * SPACING
    assert far.sum() > 10_000
    assert not np.any(claims[:, far])
    # Crust stays where it was, and with its own plate. An edge moves by up to half a cell
    # either way, so only points well clear of the open sphere (as resolved by a sample about
    # 0.75 spacings apart) are held to that exactly.
    near = dist_rad < 0.5 * SPACING
    open_dist, _ = cKDTree(sample[dist_rad > 1.5 * SPACING]).query(sample)
    deep = near & (2.0 * np.arcsin(np.clip(open_dist / 2.0, 0.0, 1.0)) > 3.0 * SPACING)
    assert deep.sum() > 1_000
    assert np.all(claims[:, deep].sum(axis=0) == 1)
    assert np.all(ids[np.argmax(claims[:, deep], axis=0)] == owners[nearest[deep]])
    assert np.mean(claims[:, near].sum(axis=0) == 1) > 0.97


def test_holes_fragments_poles_and_the_antimeridian_survive(converted_shapes):
    world, report = converted_shapes
    by_id = {plate.plate_id: plate for plate in world.plates}
    assert all(isinstance(plate, PlateWithSparseQuadPatch) and plate.node_count() > 0 for plate in world.plates)

    def holes(plate):
        return sum(surface_parity._loop_orientation(loop) < 0 for loop in plate.boundary_loops_world() if len(loop) >= 3)

    assert holes(by_id[0]) == 1
    assert [holes(by_id[k]) for k in (1, 2, 3, 4)] == [0, 0, 0, 0]
    components = {entry["plate_id"]: (entry["line_components"], entry["quad_components"]) for entry in report.plates}
    assert components == {0: (1, 1), 1: (1, 1), 2: (1, 1), 3: (1, 1), 4: (2, 2)}
    # The polar cap really contains the pole; the antimeridian disc spans both sides of it.
    assert by_id[3].contains_batch(_xyz(90, 0)[None])[0]
    assert by_id[2].contains_batch(np.stack([_xyz(0, 179), _xyz(0, -179)])).all()


def test_isolated_shapes_convert_without_stacking_and_conserve_crust(converted_shapes):
    _, report = converted_shapes
    assert report.targets["stacked"] == 0
    assert report.clamped["crustal_thickness_m"] == pytest.approx(0.0, abs=1e-12)
    for name in ("crustal_thickness_m", "mantle_lithosphere_thickness_m"):
        assert abs(report.extensive[name]["delta_vs_deduplicated"]) < 0.02
    # Small plates have long edges for their area; each edge moves by up to half a cell.
    assert all(abs(entry["area_delta"]) < 0.05 and abs(entry["crust_volume_delta"]) < 0.05 for entry in report.plates)
    assert report.preserved["provenance_violations"] == {"is_volcano": 0, "volcano_active_years_remaining": 0, "node_created_years": 0, "untargeted": 0}


def test_converted_world_passes_the_hard_invariant_audit_and_round_trips(converted_shapes):
    world, _ = converted_shapes
    violations, _ = surface_parity.audit_world(world, {})
    assert violations == []
    data = persistence.save_world_bytes(world)
    assert pickle.loads(data)["surface"] == "quad"
    assert surface_parity.state_hash(persistence.load_world_bytes(data)) == surface_parity.state_hash(world)


def test_plate_identity_and_motion_carry_over(shape_world, converted_shapes):
    world, _ = converted_shapes
    for line, quad in zip(shape_world.plates, world.plates):
        assert quad.plate_id == line.plate_id
        assert quad.crust_type == line.crust_type
        np.testing.assert_array_equal(quad.frame, line.frame)
        np.testing.assert_array_equal(quad.omega, line.omega)
        assert quad.age_steps == line.age_steps


def test_conversion_is_deterministic(shape_world):
    first, second = copy.deepcopy(shape_world), copy.deepcopy(shape_world)
    legacy_conversion.convert_world_to_quad(first)
    legacy_conversion.convert_world_to_quad(second)
    assert surface_parity.state_hash(first) == surface_parity.state_hash(second)


def test_conversion_is_one_way(converted_shapes):
    world, _ = converted_shapes
    with pytest.raises(ValueError, match="no line-backed plates"):
        legacy_conversion.convert_world_to_quad(copy.deepcopy(world))


# --- Field semantics ------------------------------------------------------------------------


def _stamped_world() -> World:
    world = _shape_world()
    plate = world.plates[2]
    count = plate.node_count()
    rng = np.random.default_rng(0)
    volcano = np.zeros(count, dtype=bool)
    volcano[count // 2] = True
    countdown = np.where(volcano, 2.5e6, 0.0)
    created = np.full(count, -1.0)
    created[: count // 3] = rng.uniform(1e6, 5e6, count // 3)
    plate.set_fields_on_plate(
        is_volcano=volcano,
        volcano_active_years_remaining=countdown,
        node_created_years=created,
        soil_depth=rng.uniform(0.0, 2.0, count),
        soil_mineral_content=rng.uniform(0.0, 1.0, count),
        elev_change_reason=rng.integers(0, 20, count).astype(float),
    )
    return world


def test_point_features_and_history_are_neither_lost_nor_multiplied():
    world = _stamped_world()
    plate = world.plates[2]
    earliest = plate.collect("node_created_years")
    earliest = earliest[earliest >= 0].min()

    report = legacy_conversion.convert_world_to_quad(world)
    quad = world.plates[2]
    assert int(quad.collect("is_volcano").sum()) == 1
    assert quad.collect("volcano_active_years_remaining")[quad.collect("is_volcano")][0] == 2.5e6
    created = quad.collect("node_created_years")
    assert created[created >= 0].min() == earliest
    assert not any(report.preserved["provenance_violations"].values())
    assert set(np.unique(quad.collect("elev_change_reason"))) <= set(range(20))
    soil = quad.collect("soil_mineral_content")
    assert np.all((soil >= 0.0) & (soil <= 1.0))


def test_every_surface_field_is_handled_by_the_transfer():
    # A field added to the registry must get a transfer rule here too: the conversion would
    # otherwise give it a default silently.
    world = _stamped_world()
    legacy_conversion.convert_world_to_quad(world)
    quad = world.plates[2]
    assert set(quad._fields) == set(SURFACE_FIELDS)
    assert {spec.remap_class for spec in SURFACE_FIELDS.values()} <= {
        RemapClass.EXTENSIVE,
        RemapClass.INTENSIVE,
        RemapClass.CATEGORICAL,
        RemapClass.BOOLEAN_PROVENANCE,
        RemapClass.COUNTDOWN,
        RemapClass.CLOCK,
        RemapClass.HISTORY,
        RemapClass.WRITE_ONCE_HISTORY,
    }


def test_lines_that_predate_optional_fields_convert_with_their_defaults():
    world = _shape_world()
    for plate in world.plates:
        for line in plate.lines:
            for name in ("node_created_years", "crust_type_code", "elev_change_reason", "overlap_onset_years"):
                line.__dict__.pop(f"_{name}", None)
    legacy_conversion.convert_world_to_quad(world)
    for plate in world.plates:
        assert np.all(plate.collect("node_created_years") == -1.0)
        assert np.all(plate.collect("overlap_onset_years") == 0.0)


def test_state_from_a_newer_build_is_refused_rather_than_dropped():
    world = _shape_world()
    world.plates[1].lines[0].__dict__["_future_field"] = np.zeros(len(world.plates[1].lines[0].theta))
    with pytest.raises(ValueError, match="_future_field"):
        legacy_conversion.convert_world_to_quad(world)

    world = _shape_world()
    world.plates[1].__dict__["_future_plate_state"] = 1
    with pytest.raises(ValueError, match="_future_plate_state"):
        legacy_conversion.convert_world_to_quad(world)


def test_corrupt_line_fields_are_refused():
    world = _shape_world()
    line = world.plates[0].lines[0]
    line.__dict__["_soil_depth"] = np.full(len(line.theta), np.nan)
    with pytest.raises(ValueError, match="non-finite"):
        legacy_conversion.convert_world_to_quad(world)

    world = _shape_world()
    line = world.plates[0].lines[0]
    line.__dict__["_soil_depth"] = np.zeros(len(line.theta) + 1)
    with pytest.raises(ValueError, match="soil_depth"):
        legacy_conversion.convert_world_to_quad(world)


def test_a_line_save_converts_on_load_exactly_as_its_world_converts_in_memory(shape_world):
    # A save written before #251 pickled the retired line classes. Loading it reads their
    # state through inert stand-ins (legacy_conversion.legacy_unpickler) and converts it: the
    # same world converting the in-memory line plates gives.
    with retired_class_names():
        data = persistence.save_world_bytes(shape_world)
    assert legacy_conversion.legacy_unpickler(data).load()["surface"] == "lines"

    loaded = persistence.load_world_bytes(data)
    direct = copy.deepcopy(shape_world)
    legacy_conversion.convert_world_to_quad(direct)

    assert all(isinstance(plate, PlateWithSparseQuadPatch) for plate in loaded.plates)
    assert loaded.surface_conversion["from"] == "lines"
    assert surface_parity.state_hash(loaded) == surface_parity.state_hash(direct)


# --- Stepped worlds -------------------------------------------------------------------------


def _line_world_like(world: World) -> World:
    return line_world_like(world, SPACING)


@pytest.fixture(scope="module")
def stepped_quad_world() -> World:
    world = generate_world(seed=7, num_plates=8, node_density=DENSITY)
    for _ in range(4):
        step_world(world, 1_000_000)
    return world


@pytest.fixture(scope="module")
def stepped_line_world(stepped_quad_world) -> World:
    return _line_world_like(stepped_quad_world)


def _land_fraction(world: World) -> float:
    """A line world's land fraction: every node stands for the same nominal footprint."""
    elevation = np.concatenate([plate.collect("elevation") for plate in world.plates])
    return float(np.mean(elevation > world.sea_level_m))


def test_sea_level_and_land_carry_over_and_the_world_keeps_stepping(stepped_quad_world, stepped_line_world):
    source = copy.deepcopy(stepped_quad_world)
    lines = copy.deepcopy(stepped_line_world)
    quad = copy.deepcopy(stepped_line_world)
    report = legacy_conversion.convert_world_to_quad(quad)

    assert quad.sea_level_m == lines.sea_level_m
    assert report.sea_level["after_m"] == report.sea_level["before_m"]
    assert abs(surface_parity.totals(quad)["land_fraction"] - _land_fraction(lines)) < 0.01
    assert abs(report.extensive["crustal_thickness_m"]["delta_vs_deduplicated"]) < 0.01
    assert report.coverage["before"] is None
    assert report.coverage["after"]["uncovered"] < 0.01
    assert quad.surface_conversion["from"] == "lines"
    assert "Converted line-backed world" in quad.events[-1][1]

    # The converted world keeps stepping cleanly and tracks the quad world its lines were
    # sampled from.
    signatures = {}
    for _ in range(2):
        step_world(source, 1_000_000)
        step_world(quad, 1_000_000)
        violations, signatures = surface_parity.audit_world(quad, signatures)
        assert violations == []
    assert abs(quad.sea_level_m - source.sea_level_m) < 50.0


def test_overlapping_plates_stack_their_crust_instead_of_dropping_it(stepped_line_world):
    # Push one plate's nodes half a spacing into its neighbours so territories overlap: the
    # nodes that lose their cell must deposit their crust on whoever won it.
    world = copy.deepcopy(stepped_line_world)
    plate = world.plates[0]
    centre = plate.all_points_and_elevation()[0].mean(axis=0)
    axis = np.cross(centre, [0.0, 0.0, 1.0])
    axis /= np.linalg.norm(axis)
    plate.rotate(geometry.rotation_matrix(axis, 3.0 * SPACING))
    report = legacy_conversion.convert_world_to_quad(world)
    assert report.targets["stacked"] > 0
    lost = report.clamped["crustal_thickness_m"]
    assert abs(report.extensive["crustal_thickness_m"]["delta_vs_deduplicated"] + lost) < 0.01
    assert not any(report.preserved["provenance_violations"].values())


def test_stacking_overflow_spreads_across_the_plate_and_conserves_volume():
    world = generate_world(seed=3, num_plates=4, node_density=0.5)
    plate = world.plates[0]
    cap = 84_000.0
    before = np.full(plate.node_count(), 30_000.0)
    stacked = before.copy()
    stacked[plate.node_count() // 2] = 250_000.0
    spread = legacy_conversion._spread_overflow(plate, stacked, before, cap)
    area = plate.node_areas_m2()
    assert np.sum(spread * area) == pytest.approx(np.sum(stacked * area), rel=1e-12)
    assert spread.max() <= cap * 1.01
    # A column a legacy save already held past the cap is not spread: the line engine clamps it.
    legacy_over = legacy_conversion._spread_overflow(plate, stacked, stacked, cap)
    np.testing.assert_array_equal(legacy_over, stacked)


def test_overlap_nodes_bring_every_field_onto_the_cell_they_stack_on():
    # Plate 0 is pushed into its neighbours, and every field on it is stamped differently from
    # theirs. Where its nodes stack onto a neighbour's cell, each field must combine by its
    # remap class (quad_merge's suture rules), not keep only the receiving cell's value.
    from app import lithosphere

    world = _line_world_like(generate_world(seed=7, num_plates=8, node_density=DENSITY))
    for plate in world.plates:
        count = plate.node_count()
        mine = plate is world.plates[0]
        plate.set_fields_on_plate(
            soil_depth=np.ones(count),
            soil_mineral_content=np.full(count, 0.9 if mine else 0.1),
            divergent_age_myr=np.full(count, 7.0 if mine else 1.0),
            channel_depth=np.ones(count),
            channel_width=np.full(count, 50.0 if mine else 10.0),
            elev_change_reason=np.full(count, 5.0 if mine else 11.0),
        )
        lithosphere.sync_plate_elevation(plate)
        if mine:
            plate.set_fields_on_plate(elevation=plate.collect("elevation") + 500.0)
    plate = world.plates[0]
    centre = plate.all_points_and_elevation()[0].mean(axis=0)
    axis = np.cross(centre, [0.0, 0.0, 1.0])
    plate.rotate(geometry.rotation_matrix(axis / np.linalg.norm(axis), 3.0 * SPACING))

    report = legacy_conversion.convert_world_to_quad(world)
    assert report.targets["stacked"] > 0
    blended, residuals, reasons = 0, [], set()
    for quad in world.plates[1:]:
        if not quad.node_count():
            continue
        soil = quad.collect("soil_mineral_content")
        receiving = (soil > 0.1 + 1e-9) & (soil < 0.9 - 1e-9)
        blended += int(receiving.sum())
        for name, low, high in (("divergent_age_myr", 1.0, 7.0), ("channel_width", 10.0, 50.0)):
            values = quad.collect(name)[receiving]
            assert np.all((values > low) & (values < high)), name
        density = lithosphere.node_crust_density(quad.collect("crust_type_code"), quad.crust_type)
        residual = quad.collect("elevation") - lithosphere.isostatic_elevation(
            quad.collect("crustal_thickness_m"), quad.collect("mantle_lithosphere_thickness_m"), density
        )
        residuals.append(residual[receiving])
        reasons |= set(np.unique(quad.collect("elev_change_reason")[receiving]))
    assert blended > 0
    residuals = np.concatenate(residuals)
    # Elevation keeps a share of the incoming crust's +500 m residual (clipping aside).
    assert np.median(residuals) > 1.0 and np.all(residuals < 500.0 + 1e-6)
    assert reasons <= {5.0, 11.0}
    assert 5.0 in {float(r) for quad in world.plates[1:] if quad.node_count() for r in np.unique(quad.collect("elev_change_reason"))}


@pytest.mark.parametrize(
    "corrupt, message",
    [
        (lambda p: p.__dict__.__setitem__("_frame", np.full((3, 3), np.nan)), "frame must be a finite"),
        (lambda p: p.__dict__.__setitem__("_frame", np.eye(3)[:2]), "frame must be a finite"),
        (lambda p: p.__dict__.__setitem__("_frame", 2.0 * np.asarray(p.__dict__["_frame"])), "proper rotation"),
        (lambda p: p.__dict__.__setitem__("_frame", -np.asarray(p.__dict__["_frame"])), "proper rotation"),
        (lambda p: p.__dict__.__setitem__("_omega", np.array([np.nan, 0.0, 0.0])), "omega"),
        (lambda p: p.__dict__.__setitem__("_omega", np.zeros(2)), "omega"),
        (lambda p: p.__dict__.__setitem__("_internal_stress", np.inf), "internal stress"),
        (lambda p: p.__dict__.__setitem__("_age_steps", -1), "age"),
        (lambda p: p.__dict__.__setitem__("_age_steps", 2.5), "age"),
        (lambda p: p.__dict__.__setitem__("_crust_type", "basaltic"), "crust type"),
        (lambda p: p.lines[0].__dict__["_theta"].__setitem__(0, np.nan), "coordinates"),
        (lambda p: p.lines[0].__dict__.__setitem__("_phi", np.inf), "coordinates"),
    ],
)
def test_malformed_plate_state_is_refused(corrupt, message):
    world = _shape_world()
    corrupt(world.plates[2])
    with pytest.raises(ValueError, match=message):
        legacy_conversion.convert_world_to_quad(world)

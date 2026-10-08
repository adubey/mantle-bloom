"""The issue #247 parity harness: its audits catch the faults they claim to, every topology and
geometry mutation reaches the right cache-invalidation path (the derived-index audit asked for
on #247), and its gates judge results as documented in docs/surface-parity.md."""

import pickle

import numpy as np
import pytest

from app import eustasy, lithosphere, persistence, plates as plates_mod, surface_parity as sp, surface_parity_gates as gates
from app.elevation_lines import line_spacing_rad
from app.plates import gather_node_positions
from app.sparse_quad_patch import PlateWithSparseQuadPatch
from app.world import World

DENSITY = 0.5
SPACING = line_spacing_rad(DENSITY)
_Q, _ = np.linalg.qr(np.random.default_rng(3).normal(size=(3, 3)))
ROTATED = _Q * np.sign(np.linalg.det(_Q))


def _cap(plate_id=1, frame=None, centre=(1.0, 0.0, 0.0), radius=6 * SPACING) -> PlateWithSparseQuadPatch:
    centre = np.asarray(centre, dtype=float)
    plate = PlateWithSparseQuadPatch.from_lattice(
        plate_id, np.eye(3) if frame is None else frame, "continental", SPACING, lambda pts: pts @ centre > np.cos(radius)
    )
    hc, hm = lithosphere.reference_thickness("continental")
    n = plate.node_count()
    plate.set_fields_on_plate(crustal_thickness_m=np.full(n, hc), mantle_lithosphere_thickness_m=np.full(n, hm))
    lithosphere.sync_plate_elevation(plate)
    return plate


def _prime(plate) -> None:
    """Populate every derived cache the audit knows about."""
    for builder in sp.cache_builders(plate).values():
        getattr(plate, builder)()


def _populated(plate) -> list[str]:
    return [name for name in sp.cache_builders(plate) if getattr(plate, name, None) is not None]


# --- Audit sensitivity ----------------------------------------------------------------------


def test_a_healthy_quad_plate_has_no_audit_findings():
    plate = _cap()
    _prime(plate)
    assert sp.stale_plate_caches(plate) == []
    assert sp.quad_topology_violations(plate) == []
    assert sp.quad_neighbour_violations(plate) == []
    assert not sp.quad_cell_shape(plate)["folded"].any()
    assert sp.field_violations(plate) == [] and sp.field_bound_violations(plate) == []
    assert sp.frame_violation(plate) is None


def test_audit_catches_a_topology_edit_that_bypasses_invalidation():
    plate = _cap()
    _prime(plate)
    before = sp.plate_signature(plate)
    # A buggy edit: drop a cell without _replace_topology (no revision bump, no cache reset).
    plate._keys = plate._keys[:-1].copy()
    plate._fields = {name: values[:-1].copy() for name, values in plate._fields.items()}

    assert {"_node_kdtree_cache", "_adjacency_cache", "_world_points_cache"} <= set(sp.stale_plate_caches(plate))
    assert sp.revision_violations(before, sp.plate_signature(plate)) == ["local node set changed without a topology revision bump"]


def test_audit_catches_a_frame_change_without_a_geometry_revision():
    plate = _cap()
    before = sp.plate_signature(plate)
    plate._frame = ROTATED.copy()
    assert sp.revision_violations(before, sp.plate_signature(plate)) == ["frame changed without a geometry revision bump"]


def test_audit_catches_broken_neighbours():
    plate = _cap()
    adjacency = plate._neighbour_indices().copy()
    far = int(np.argmax(np.linalg.norm(plate._local_centres() - plate._local_centres()[0], axis=1)))
    adjacency[0, np.flatnonzero(adjacency[0] < 0)[0]] = far  # one-way, non-adjacent entry
    plate._adjacency_cache = adjacency

    problems = sp.quad_neighbour_violations(plate)
    assert "adjacency is not symmetric" in problems
    assert any("do not share an edge" in p for p in problems)


def test_neighbour_check_accepts_2_to_1_level_jumps_and_cube_face_seams():
    plate = _cap(radius=30 * SPACING, centre=(1.0, 1.0, 0.0) / np.sqrt(2))  # straddles faces
    plate.refine_cells(plate.cell_keys[:5])
    assert sp.quad_neighbour_violations(plate) == []
    assert sp.quad_topology_violations(plate) == []


def test_folded_cell_detection(monkeypatch):
    plate = _cap()
    corners = sp.quad_cell_corners_local(plate)
    monkeypatch.setattr(sp, "quad_cell_corners_local", lambda _: corners[:, ::-1])
    assert sp.quad_cell_shape(plate)["folded"].all()


def test_field_audit_catches_non_finite_and_out_of_range_values():
    plate = _cap()
    elevation = plate.collect("elevation")
    elevation[0] = np.nan
    hc = plate.collect("crustal_thickness_m")
    hc[1] = 10.0
    plate.set_fields_on_plate(elevation=elevation, crustal_thickness_m=hc)
    assert any(p.startswith("elevation: 1 non-finite") for p in sp.field_violations(plate))
    assert any(p.startswith("crustal_thickness_m: 1 values outside") for p in sp.field_bound_violations(plate))


def test_a_replacement_plate_reusing_an_id_is_not_a_revision_regression():
    plate = _cap()
    plate.rotate(ROTATED)
    before = sp.plate_signature(plate)
    centre = plate._get_world_points().mean(axis=0)
    left, _ = plate.split(5, np.cross(centre, [0.0, 0.0, 1.0]), min_nodes=3)  # keeps plate_id, fresh counters
    assert left.plate_id == plate.plate_id and left.geometry_revision < before["geometry_revision"]
    assert sp.revision_violations(before, sp.plate_signature(left)) == []


def test_world_cache_audit_catches_an_index_built_before_a_plate_moved():
    world = World(seed=0, plates=[_cap()], next_plate_id=2, node_density=DENSITY)
    points, _ = gather_node_positions(world.plates)
    plates_mod.cached_node_position_tree(world, points)
    assert sp.stale_world_caches(world) == []

    world.plates[0].rotate(ROTATED)
    assert sp.stale_world_caches(world) == ["node_position_tree_cache"]


# --- Every mutation reaches the right invalidation path --------------------------------------


def _insert(plate):
    plate.insert_cells(np.asarray(plate.empty_neighbour_keys(np.arange(plate.node_count()))[0][:3]))
    return plate


def _merge(plate):
    other = _cap(9, ROTATED, centre=plate._get_world_points()[0], radius=3 * SPACING)
    plate.merge_with(other, other.all_points_and_elevation()[0])
    return plate


QUAD_TOPOLOGY_EDITS = {
    "insert_cells": _insert,
    "remove_cells": lambda p: (p.remove_cells(np.arange(p.node_count()) < 3), p)[1],
    "refine_cells": lambda p: (p.refine_cells(p.cell_keys[:2]), p)[1],
    "merge_with": _merge,
}


@pytest.mark.parametrize("edit", sorted(QUAD_TOPOLOGY_EDITS))
def test_quad_topology_edits_bump_both_revisions_and_leave_no_stale_cache(edit):
    plate = _cap()
    _prime(plate)
    before = sp.plate_signature(plate)
    plate = QUAD_TOPOLOGY_EDITS[edit](plate)
    after = sp.plate_signature(plate)

    assert after["local"] != before["local"]
    assert after["topology_revision"] > before["topology_revision"]
    assert after["geometry_revision"] > before["geometry_revision"]
    assert sp.stale_plate_caches(plate) == []
    assert sp.quad_topology_violations(plate) == [] and sp.quad_neighbour_violations(plate) == []


def test_quad_coarsen_invalidates_like_refine():
    plate = _cap()
    lineage = plate.refine_cells(plate.cell_keys[:1])
    _prime(plate)
    before = sp.plate_signature(plate)
    plate.coarsen_cells(np.array(list(lineage)[:1]))
    assert sp.plate_signature(plate)["topology_revision"] > before["topology_revision"]
    assert sp.stale_plate_caches(plate) == []


def test_quad_rotation_bumps_geometry_only_and_keeps_local_caches():
    plate = _cap()
    _prime(plate)
    adjacency, loops = plate._adjacency_cache, plate._local_loops_cache
    before = sp.plate_signature(plate)
    plate.rotate(ROTATED)
    after = sp.plate_signature(plate)

    assert after["topology_revision"] == before["topology_revision"]
    assert after["geometry_revision"] == before["geometry_revision"] + 1
    assert plate._adjacency_cache is adjacency and plate._local_loops_cache is loops
    assert plate._node_kdtree_cache is None and plate._world_points_cache is None
    _prime(plate)
    assert sp.stale_plate_caches(plate) == []


def test_field_only_writes_keep_position_indexes():
    plate = _cap()
    tree, outline = plate.get_node_kdtree(), plate.get_bounding_polygon()
    revisions = (plate.topology_revision, plate.geometry_revision)
    plate.set_fields_on_plate(elevation=plate.collect("elevation") + 1.0)

    assert (plate.topology_revision, plate.geometry_revision) == revisions
    assert plate.get_node_kdtree() is tree and plate.get_bounding_polygon() is outline


def test_quad_split_and_failed_rift_leave_consistent_caches():
    plate = _cap(radius=10 * SPACING)
    _prime(plate)
    halves = plate.split(5, np.array([0.0, 1.0, 0.0]), min_nodes=3)
    assert halves is not None
    for half in halves:
        assert sp.stale_plate_caches(half) == [] and sp.quad_neighbour_violations(half) == []

    _prime(plate)
    before = sp.plate_signature(plate)
    plate.apply_failed_rift(np.array([0.0, 1.0, 0.0]), SPACING)
    assert sp.plate_signature(plate)["local"] == before["local"]  # topology unchanged
    assert sp.stale_plate_caches(plate) == []


def test_pickling_drops_every_derived_cache():
    plate = _cap()
    _prime(plate)
    assert _populated(plate)
    loaded = pickle.loads(pickle.dumps(plate))
    assert _populated(loaded) == []
    assert sp.stale_plate_caches(loaded) == []


def test_world_load_drops_every_derived_index_and_keeps_authoritative_state():
    world = World(seed=0, plates=[_cap()], next_plate_id=3, node_density=DENSITY)
    # A generated world always has its water budget; without one, load's backfill would
    # rebuild the area caches it reads.
    eustasy.initialize_water_budget(world)
    for plate in world.plates:
        _prime(plate)
    points, _ = gather_node_positions(world.plates)
    plates_mod.cached_node_position_tree(world, points)
    plates_mod.cached_node_healpix_index(world, points)

    result, loaded = sp.load_check(world)
    assert result == {"state_identical": True, "derived_state_after_load": []}
    assert sp.state_hash(loaded) == sp.state_hash(world)


def test_state_hash_ignores_caches_but_not_fields():
    world = World(seed=0, plates=[_cap()], next_plate_id=2, node_density=DENSITY)
    digest = sp.state_hash(world)
    _prime(world.plates[0])
    assert sp.state_hash(world) == digest
    world.plates[0].set_fields_on_plate(soil_depth=np.ones(world.plates[0].node_count()))
    assert sp.state_hash(world) != digest


# --- Instrumentation ------------------------------------------------------------------------


def test_instrumentation_restores_everything_it_patched():
    originals = {(id(owner), name): (owner.__dict__[name] if isinstance(owner, type) else getattr(owner, name)) for owner, name, _ in sp._PHASES}
    instrumentation = sp.Instrumentation()
    with instrumentation.installed():
        assert plates_mod.cached_node_position_tree is not originals.get((id(plates_mod), "cached_node_position_tree"))
        plate = _cap()
        plate.get_node_kdtree()
        plate.get_node_kdtree()
    for owner, name, _ in sp._PHASES:
        current = owner.__dict__[name] if isinstance(owner, type) else getattr(owner, name)
        assert current is originals[(id(owner), name)]
    assert instrumentation.index_builds["plate_node_kdtree"] == 1
    assert instrumentation.index_hits["plate_node_kdtree"] == 1


def test_instrumentation_times_generator_phases_and_does_not_double_count_nesting():
    instrumentation = sp.Instrumentation()

    def inner():
        return 1

    def outer():
        yield wrapped_inner()
        yield wrapped_inner()

    wrapped_inner = instrumentation._timed(inner, "inner")
    wrapped_outer = instrumentation._timed(outer, "outer")
    assert list(wrapped_outer()) == [1, 1]
    assert "outer" in instrumentation.phase_seconds and "inner" not in instrumentation.phase_seconds


# --- Gates ----------------------------------------------------------------------------------


def _checkpoint(age, **overrides):
    checkpoint = {
        "elapsed_myr": age,
        "step": int(age),
        "totals": {
            "plates": 10,
            "nodes": 1000,
            "area_is_exact": False,
            "area_over_sphere": 1.0,
            "hc_volume_km3": 100.0,
            "hm_volume_km3": 200.0,
            "continental_hc_volume_km3": 50.0,
            "land_fraction": 0.3,
            "sea_level_m": 0.0,
            "elevation": {"p05": -5000.0, "p50": -3000.0, "p95": 1000.0},
        },
        "coverage": {"uncovered": 0.02, "multiply_covered": 0.01, "void": 0.0, "nodes_inside_other_plate": 0.01},
        "sample_cloud": {"stacked": 0.0, "anisotropic": 0.0, "anisotropic_row_alignment": None, "thin": 0.001, "extra_components": 0, "isolated_nodes": 0},
        "quad_lattice": None,
        "climate_hydrology": {"stats": {"air_temperature_mean_c": 10.0}, "hydrology": {"finite": True}},
        "index_parity": {"cached_kdtree_mismatches": 0},
    }
    for path, value in overrides.items():
        target = checkpoint
        *parents, leaf = path.split(".")
        for key in parents:
            target = target[key]
        target[leaf] = value
    return checkpoint


def _run(seed=1, checkpoints=None, violations=(), series=None):
    return {
        "config": {"name": "test"},
        "seed": seed,
        "checkpoints": checkpoints or [_checkpoint(0.0), _checkpoint(10.0)],
        "violations": list(violations),
        "audits": 3,
        "series": series or [],
        "load_checks": [{"step": 10, "state_identical": True, "derived_state_after_load": [], "continuation_identical": True}],
    }


def _gate(evaluation, key):
    return evaluation["gates"][key]["status"]


def test_a_clean_run_passes():
    evaluation = gates.evaluate([(_run(), None)])
    assert evaluation["verdict"] == gates.PASS


def test_a_hard_invariant_violation_fails():
    violation = {"kind": "quad_folded", "detail": "1 folded cells", "step": 4}
    evaluation = gates.evaluate([(_run(violations=[violation]), None)])
    assert evaluation["verdict"] == gates.FAIL and _gate(evaluation, "H3:quad_folded") == gates.FAIL

    stale = {"kind": "stale_plate_cache", "detail": "_node_kdtree_cache", "step": 4}
    assert _gate(gates.evaluate([(_run(violations=[stale]), None)]), "H6:derived_caches") == gates.FAIL


def test_a_cap_breach_warns():
    breach = {"kind": "field_bounds", "detail": "mantle_lithosphere_thickness_m: 1 values outside", "step": 4}
    evaluation = gates.evaluate([(_run(violations=[breach]), None)])
    assert _gate(evaluation, "H11:field_caps") == gates.WARN and evaluation["verdict"] == gates.WARN


def test_a_void_fails():
    run = _run(checkpoints=[_checkpoint(0.0), _checkpoint(10.0, **{"coverage.void": 0.001})])
    assert _gate(gates.evaluate([(run, None)]), "C2:void") == gates.FAIL


def test_multiply_covered_uses_its_exact_cell_tolerance():
    def status(value):
        run = _run(checkpoints=[_checkpoint(0.0), _checkpoint(10.0, **{"coverage.multiply_covered": value})])
        return gates.evaluate([(run, None)])

    assert _gate(status(0.014), "C3:multiply_covered") == gates.PASS
    result = status(0.019)
    assert _gate(result, "C3:multiply_covered") == gates.WARN
    assert result["gates"]["C3:multiply_covered"]["worst"]["fail_above"] == gates.QUAD_MULTIPLY_COVERED_FAIL
    assert _gate(status(0.021), "C3:multiply_covered") == gates.FAIL


def test_area_accounting_must_match_the_covered_sphere():
    exact = {"totals.area_is_exact": True, "totals.area_over_sphere": 0.99}
    run = _run(checkpoints=[_checkpoint(0.0), _checkpoint(10.0, **exact)])
    assert _gate(gates.evaluate([(run, None)]), "K3:quad_area_accounting") == gates.PASS
    run = _run(checkpoints=[_checkpoint(0.0), _checkpoint(10.0, **{**exact, "totals.area_over_sphere": 0.9})])
    assert _gate(gates.evaluate([(run, None)]), "K3:quad_area_accounting") == gates.WARN


def test_report_renders_every_group_present():
    timings = {"steps": [{"wall_s": 1.2, "phases": {"deform": 0.2, "climate_erosion_hydrology": 1.0}, "index_builds": {}}]}
    evaluation = gates.evaluate([(_run(), timings)])
    report = gates.render_report(evaluation, [(_run(), timings)], ["cmd"])
    assert "**Verdict: PASS**" in report and "## Hard invariants" in report and "## Checkpoints, seed 1" in report
    assert "## Performance" in report


def test_presets_have_consistent_checkpoints():
    for config in sp.PRESETS.values():
        assert config.checkpoint_steps[0] == 0 and config.steps == config.checkpoint_steps[-1]
    with pytest.raises(ValueError):
        sp.RunConfig("bad", 1.0, 1e6, (2.5,), 1, 10).checkpoint_steps


def test_run_can_continue_a_saved_world(tmp_path):
    from app.world import generate_world

    world = generate_world(seed=4, node_density=0.25)
    path = tmp_path / "w.mbworld"
    path.write_bytes(persistence.save_world_bytes(world))
    config = sp.RunConfig("resume", node_density=0.25, step_years=1e6, checkpoints_myr=(), audit_every=1, samples=500, load_checks=False)

    document, _ = sp.run_world(config, 4, initial_world=path, log=lambda _: None)
    assert document["initial_world"]["name"] == "w.mbworld"
    assert document["checkpoints"][0]["totals"]["nodes"] == sum(p.node_count() for p in world.plates)
    with pytest.raises(ValueError):
        sp.run_world(config, 5, initial_world=path, log=lambda _: None)


def test_hydrology_finite_check_allows_the_no_rim_sentinel_only():
    from types import SimpleNamespace

    n = 4
    fields = SimpleNamespace(
        points=np.zeros((n, 3)),
        is_ocean=np.zeros(n, dtype=bool),
        is_river=np.zeros(n, dtype=bool),
        lake_depth=np.zeros(n),
        glacier_depth=np.zeros(n),
        flow_target=np.zeros(n, dtype=int),
        flow_accum=np.ones(n),
        filled_elevation=np.array([0.0, 1.0, np.inf, 2.0]),
    )
    world = SimpleNamespace(stats_history=[], hydrology_cache=fields)
    assert sp.climate_hydrology(world)["hydrology"]["finite"]
    fields.filled_elevation[0] = np.nan
    assert not sp.climate_hydrology(world)["hydrology"]["finite"]

import numpy as np
import pytest
from app import faults, persistence
from app.elevation_lines import (
    ELEV_CHANGE_FAULT_NORMAL,
    ELEV_CHANGE_FAULT_REVERSE,
    ELEV_CHANGE_FAULT_STRIKE_SLIP,
)
from app.faults import (
    Fault,
    FaultSystem,
    _KIND_NORMAL,
    _KIND_REVERSE,
    _KIND_STRIKE_SLIP,
    _cull_inactive_systems,
    _cull_scars,
    _regime_from_closing,
    reconcile_faults,
    update_faults,
)
from app.world import World, generate_world, step_world

# A high BASE_SPAWN_RATE keeps the step-count (and so the runtime) of the integration tests
# low while still exercising spawn/age/retire/reconcile -- the default rate is tuned for
# plausibility over a long run, not for a 15-step test.
pytestmark = pytest.mark.usefixtures("_fast_faults")


@pytest.fixture
def _fast_faults(monkeypatch):
    monkeypatch.setattr(faults, "BASE_SPAWN_RATE_PER_MYR", 40.0)
    monkeypatch.setattr(faults, "LIFESPAN_MIN_MYR", 2.0)
    monkeypatch.setattr(faults, "LIFESPAN_MAX_MYR", 6.0)
    # Systems too: spawn them often, and let them retire inside a ~15-step test.
    monkeypatch.setattr(faults, "SYSTEM_SPAWN_FRACTION", 0.5)
    monkeypatch.setattr(faults, "SYSTEM_LIFESPAN_MIN_MYR", 4.0)
    monkeypatch.setattr(faults, "SYSTEM_LIFESPAN_MAX_MYR", 10.0)


def _run(seed: int, *, plates: int = 8, steps: int = 15, dt: float = 1_000_000) -> World:
    world = generate_world(seed=seed, num_plates=plates)
    for _ in range(steps):
        step_world(world, dt)
    return world


def _fault(**kw) -> Fault:
    """A minimal Fault for the pure-helper tests -- geometry is a short arc near (1, 0, 0) in
    the identity frame; only the fields a given test reads need to be meaningful."""
    n = kw.pop("n_nodes", 5)
    ang = np.linspace(-0.02, 0.02, n)
    phi = kw.pop("local_phi", ang)
    theta = kw.pop("local_theta", np.zeros(n))
    defaults = dict(
        fault_id=0,
        plate_id=0,
        kind=_KIND_REVERSE,
        local_phi=phi,
        local_theta=theta,
        slip_rate_m_per_myr=3000.0,
        dip_deg=30.0,
        strike_sense=1,
        dip_dir_local=np.array([0.0, 1.0, 0.0]),
        lifespan_myr=10.0,
        birth_years=0.0,
        birth_distance_from_boundary_km=100.0,
    )
    defaults.update(kw)
    return Fault(**defaults)


# --------------------------------------------------------------------------- pure helpers


def test_regime_from_closing_follows_andersonian_sign():
    from app import boundary

    thr = boundary.TRANSFORM_RATE_THRESHOLD
    assert _regime_from_closing(2.0 * thr) == _KIND_REVERSE  # shortening -> thrust
    assert _regime_from_closing(-2.0 * thr) == _KIND_NORMAL  # extension -> normal
    assert _regime_from_closing(0.0) == _KIND_STRIKE_SLIP  # neither -> wrench


def test_length_km_is_the_summed_great_circle_arc():
    n = 6
    f = _fault(local_phi=np.zeros(n), local_theta=np.linspace(0.0, 0.1, n))
    assert f.length_km() == pytest.approx(0.1 * faults.PLANET_RADIUS_KM, rel=1e-6)


def test_cull_scars_caps_inactive_scars_per_plate_keeping_the_youngest():
    world = World(seed=0, plates=[])
    n = faults.MAX_SCARS_PER_PLATE + 5
    for i in range(n):
        world.faults.append(_fault(fault_id=i, active=False, age_myr=float(i)))  # older = higher age
    world.faults.append(_fault(fault_id=999, active=True, age_myr=500.0))

    _cull_scars(world)

    kept = world.faults
    assert sum(not f.active for f in kept) == faults.MAX_SCARS_PER_PLATE
    assert any(f.fault_id == 999 for f in kept)  # active never culled
    assert min(f.age_myr for f in kept if not f.active) == 5.0  # the 5 oldest scars dropped


def test_cull_scars_is_per_plate_not_global():
    world = World(seed=0, plates=[])
    for pid in (0, 1):
        for i in range(faults.MAX_SCARS_PER_PLATE - 1):
            world.faults.append(_fault(fault_id=pid * 1000 + i, plate_id=pid, active=False, age_myr=float(i)))
    _cull_scars(world)
    assert len(world.faults) == 2 * (faults.MAX_SCARS_PER_PLATE - 1)  # neither plate over the cap


def test_update_faults_noop_for_empty_world():
    world = World(seed=0, plates=[])
    assert update_faults(world, years=1_000_000) is None
    assert world.faults == []


# --------------------------------------------------------------------------- spawn / age / retire


def test_faults_accumulate_with_a_mix_of_regimes_and_lock_up_as_scars():
    world = _run(seed=3, steps=10)

    assert len(world.faults) > 5
    kinds = {f.kind for f in world.faults}
    assert kinds <= {_KIND_NORMAL, _KIND_REVERSE, _KIND_STRIKE_SLIP}
    assert len(kinds) >= 2

    matured = [f for f in world.faults if f.age_myr > 0.0]
    assert matured and all(f.cumulative_offset_m > 0.0 for f in matured if f.active)
    assert any(not f.active for f in world.faults)  # some have locked up
    for f in world.faults:
        if not f.active:
            assert f.age_myr >= f.lifespan_myr  # only retires past its drawn lifespan

    messages = " ".join(m for _, m in world.events)
    assert "locked up" in messages


def test_spawning_is_deterministic_per_seed():
    key = lambda w: sorted(
        (f.fault_id, f.kind, round(f.cumulative_offset_m, 3), round(f.age_myr, 6)) for f in w.faults
    )
    assert key(_run(seed=11, steps=10)) == key(_run(seed=11, steps=10))
    assert key(_run(seed=11, steps=10)) != key(_run(seed=12, steps=10))


# --------------------------------------------------------------------------- relief


def test_fault_relief_stamps_a_fault_reason_code_on_nearby_crust():
    world = _run(seed=7, steps=10)
    assert any(f.active for f in world.faults)

    fault_reasons = {ELEV_CHANGE_FAULT_NORMAL, ELEV_CHANGE_FAULT_REVERSE, ELEV_CHANGE_FAULT_STRIKE_SLIP}
    stamped: set[int] = set()
    for plate in world.plates:
        for line in plate.lines:
            stamped |= set(np.unique(line.elev_change_reason).astype(int))
    assert stamped & fault_reasons


def test_reverse_fault_uplifts_and_normal_fault_drops_its_hanging_wall():
    # Drive one plate's crust past a single hand-placed fault of each kind and check the
    # sign of the relief it applies, isolated from the spawn model.
    world = generate_world(seed=2, num_plates=6)
    step_world(world, 1_000_000)
    plate = max(world.plates, key=lambda p: p.node_count())

    def relief_delta(kind: str) -> np.ndarray:
        pts = plate.all_points_and_elevation()[0]
        mid = np.asarray(pts[len(pts) // 2])
        local_mid = mid @ plate.frame  # to_local for an orthonormal frame
        phi0 = float(np.arcsin(np.clip(local_mid[2], -1, 1)))
        theta0 = float(np.arctan2(local_mid[1], local_mid[0]))
        span = np.linspace(-0.03, 0.03, 6)
        f = _fault(
            kind=kind,
            local_phi=np.full(6, phi0),
            local_theta=theta0 + span,
            dip_dir_local=np.array([np.cos(phi0 + 0.5), 0.0, np.sin(phi0 + 0.5)]),
            slip_rate_m_per_myr=faults.SLIP_RATE_REF_M_PER_MYR,
            lifespan_myr=1e9,
            plate_id=plate.plate_id,
        )
        world.faults = [f]
        before = np.concatenate([ln.elevation.copy() for ln in plate.lines])
        faults._apply_plate_fault_relief(world, plate, years_myr=1.0)
        after = np.concatenate([ln.elevation for ln in plate.lines])
        return after - before

    rev = relief_delta(_KIND_REVERSE)
    assert np.max(rev) > 0.0 and np.min(rev) >= -1e-9  # thrust only pushes up

    nrm = relief_delta(_KIND_NORMAL)
    assert np.min(nrm) < 0.0 and np.max(nrm) > 0.0  # graben down, footwall shoulder up


# --------------------------------------------------------------------------- topology reconciliation


def test_reconcile_faults_drops_faults_whose_plate_is_gone():
    world = generate_world(seed=4, num_plates=6)
    step_world(world, 1_000_000)
    ghost = _fault(fault_id=42, plate_id=99_999)  # a plate id that does not exist
    ghost.world_polyline = None
    world.faults = [ghost]
    reconcile_faults(world)
    assert world.faults == []


def test_reconcile_faults_rehomes_a_fault_onto_the_plate_now_under_its_trace():
    world = generate_world(seed=4, num_plates=6)
    step_world(world, 1_000_000)
    donor = max(world.plates, key=lambda p: p.node_count())
    other = min((p for p in world.plates if p.plate_id != donor.plate_id), key=lambda p: p.node_count())

    other_pts = other.all_points_and_elevation()[0]
    trace_world = np.asarray(other_pts[: min(6, len(other_pts))], dtype=float)
    trace_world = trace_world / np.linalg.norm(trace_world, axis=-1, keepdims=True)
    local = trace_world @ donor.frame
    phi = np.arcsin(np.clip(local[:, 2], -1, 1))
    theta = np.arctan2(local[:, 1], local[:, 0])
    f = _fault(fault_id=7, plate_id=donor.plate_id, local_phi=phi, local_theta=theta)
    f.world_polyline = trace_world
    world.faults = [f]

    reconcile_faults(world)

    assert len(world.faults) == 1
    assert world.faults[0].plate_id == other.plate_id  # re-homed
    assert world.faults[0].age_myr == f.age_myr  # history preserved


def test_faults_never_dangle_on_a_vanished_plate_over_a_run():
    world = _run(seed=9, steps=10)
    live = {p.plate_id for p in world.plates}
    assert all(f.plate_id in live for f in world.faults)


# --------------------------------------------------------------------------- fault systems


def test_fault_systems_spawn_with_long_master_traces_and_strand_families():
    world = _run(seed=3, steps=12)
    assert world.fault_systems, "expected some fault systems"

    # Master lineaments are an order of magnitude longer than a lone fault's ~200 km cap.
    assert max(s.master_length_km() for s in world.fault_systems) > 600.0
    for s in world.fault_systems:
        assert s.kind in {_KIND_NORMAL, _KIND_REVERSE, _KIND_STRIKE_SLIP}
        assert len(s.master_local_phi) >= 6

    # Strands carry a valid system_id, and can run past the lone-fault LENGTH_MAX_KM (200)
    # -- the widened distribution.
    live_system_ids = {s.system_id for s in world.fault_systems}
    strand_system_ids = {f.system_id for f in world.faults if f.system_id is not None}
    assert strand_system_ids
    assert strand_system_ids & live_system_ids  # strands point at systems that still exist
    assert max(f.length_km() for f in world.faults if f.system_id is not None) > faults.LENGTH_MAX_KM

    # Lone faults (no system) still obey the original tight length cap.
    lone = [f.length_km() for f in world.faults if f.system_id is None]
    if lone:
        assert max(lone) <= faults.LENGTH_MAX_KM + 1.0


def test_fault_systems_age_and_go_inactive():
    world = _run(seed=8, steps=15)
    assert any(not s.active for s in world.fault_systems), "some systems should have timed out"
    for s in world.fault_systems:
        if not s.active:
            assert s.age_myr >= s.lifespan_myr
    assert "fault system" in " ".join(m for _, m in world.events)


def test_cull_inactive_systems_caps_scars_per_plate_keeping_youngest():
    world = World(seed=0, plates=[])
    n = faults.MAX_INACTIVE_SYSTEMS_PER_PLATE + 4
    for i in range(n):
        world.fault_systems.append(
            FaultSystem(
                system_id=i, plate_id=0, kind=_KIND_REVERSE,
                master_local_phi=np.zeros(3), master_local_theta=np.zeros(3),
                length_km=2000.0, birth_years=0.0, lifespan_myr=1.0,
                age_myr=float(i), active=False,
            )
        )
    world.fault_systems.append(
        FaultSystem(
            system_id=999, plate_id=0, kind=_KIND_REVERSE,
            master_local_phi=np.zeros(3), master_local_theta=np.zeros(3),
            length_km=2000.0, birth_years=0.0, lifespan_myr=1.0, age_myr=500.0, active=True,
        )
    )
    _cull_inactive_systems(world)
    assert sum(not s.active for s in world.fault_systems) == faults.MAX_INACTIVE_SYSTEMS_PER_PLATE
    assert any(s.system_id == 999 for s in world.fault_systems)


def test_fault_systems_round_trip_through_a_save(tmp_path):
    world = _run(seed=7, steps=10)
    assert world.fault_systems
    save = tmp_path / "sys7.mbworld"
    save.write_bytes(persistence.save_world_bytes(world))
    loaded = persistence.load_world_bytes(save.read_bytes())
    assert len(loaded.fault_systems) == len(world.fault_systems)
    assert loaded.next_fault_system_id == world.next_fault_system_id
    step_world(loaded, 1_000_000)  # keeps stepping cleanly


def test_old_save_without_fault_systems_field_still_loads(tmp_path):
    world = generate_world(seed=7, num_plates=6)
    del world.fault_systems
    save = tmp_path / "old_sys.mbworld"
    save.write_bytes(persistence.save_world_bytes(world))
    loaded = persistence.load_world_bytes(save.read_bytes())
    assert loaded.fault_systems == []


# --------------------------------------------------------------------------- persistence


def test_old_save_without_faults_field_still_loads(tmp_path):
    world = generate_world(seed=7, num_plates=6)
    del world.faults  # simulate a pickle written before the field existed
    save = tmp_path / "old.mbworld"
    save.write_bytes(persistence.save_world_bytes(world))
    loaded = persistence.load_world_bytes(save.read_bytes())
    assert loaded.faults == []


def test_faults_round_trip_through_a_save(tmp_path):
    world = _run(seed=7, steps=8)
    assert world.faults
    save = tmp_path / "seed7.mbworld"
    save.write_bytes(persistence.save_world_bytes(world))
    loaded = persistence.load_world_bytes(save.read_bytes())
    assert len(loaded.faults) == len(world.faults)
    assert loaded.next_fault_id == world.next_fault_id
    step_world(loaded, 1_000_000)
    assert len({f.fault_id for f in loaded.faults}) == len(loaded.faults)  # no id collisions

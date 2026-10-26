import json

import numpy as np
import pytest
from app import cratons, persistence, replay_log
from app.replay_log import Health, ReplayLog, StepOutcome


@pytest.fixture(scope="module")
def saved_world() -> bytes:
    from app.world import generate_world

    world = generate_world(seed=5, num_plates=2, node_density=0.5)
    # Tectonics only: a quarter of the step time, and these tests don't need climate.
    world.simulate_climate_biomes = False
    return persistence.save_world_bytes(world)


@pytest.fixture
def world(saved_world):
    return persistence.load_world_bytes(saved_world)


def _fake_step(mutate):
    """A step that advances the clock and applies `mutate(world)`, for forcing a trigger."""
    def step(world, years):
        world.elapsed_years += years
        mutate(world)
        yield 1.0
    return step


def _run(log, world, years=1_000.0):
    for _ in log.step(world, years):
        pass


def _drop_plate(world):
    world.plates = world.plates[:1]


def _break_continental_ledger(world):
    world.continental_material_ledger["initial_continental_m3"] *= 2.0


def _book_craton_loss(world):
    world.craton_ledger["subducted_m3"] += 0.05 * cratons.live_volume_m3(world)


def _raise(world):
    raise RuntimeError("boom")


@pytest.mark.parametrize(
    "mutate, expected",
    [
        (_book_craton_loss, "craton_step_loss"),
        # The same booking without the craton leaving the plates breaks the balance.
        (_book_craton_loss, "craton_ledger_balance"),
        (_break_continental_ledger, "continental_ledger_open"),
        (_drop_plate, "plate_count_jump"),
        (lambda world: setattr(world, "sea_level_m", world.sea_level_m + 5_000.0), "land_fraction_jump"),
    ],
)
def test_each_trigger_writes_a_record_with_the_sidecar(tmp_path, world, mutate, expected):
    log = ReplayLog(tmp_path, step_progress=_fake_step(mutate))
    elapsed = world.elapsed_years
    _run(log, world)
    log.flush()

    [record] = log.records()
    assert expected in record["triggers"]
    sidecar = json.loads((tmp_path / record["name"] / replay_log.RECORD_FILE).read_text())
    assert expected in sidecar["triggers"]
    assert sidecar["step_years"] == 1_000.0
    assert sidecar["elapsed_years_before"] == elapsed
    assert sidecar["seed"] == world.seed
    assert sidecar["node_density"] == world.node_density
    assert sidecar["controls"]["simulate_climate_biomes"] is False
    assert set(sidecar) >= {"git", "step_seconds", "recorded_at", "health_before", "health_after", "fingerprint_after"}
    before = persistence.load_world_bytes((tmp_path / record["name"] / replay_log.BEFORE_FILE).read_bytes())
    assert before.elapsed_years == elapsed


def test_continental_area_jump_fires(tmp_path, world, monkeypatch):
    # Raising the host threshold mid-step leaves no genuine continent behind.
    log = ReplayLog(tmp_path, step_progress=_fake_step(lambda w: monkeypatch.setattr(cratons, "CRATON_HOST_MIN_HC_M", 1e9)))
    _run(log, world)
    log.flush()
    [record] = log.records()
    assert "continental_area_jump" in record["triggers"]


@pytest.mark.parametrize("trigger, field", [
    ("hc_at_max", "crustal_thickness_m"),
    ("hm_at_max", "mantle_lithosphere_thickness_m"),
    ("glacier_at_max", "glacier_depth"),
])
def test_a_step_that_puts_nodes_at_a_thickness_cap_fires(tmp_path, world, trigger, field):
    cap = replay_log._CAPPED_FIELDS[f"{trigger}_fraction"][1]

    def to_cap(world):
        plate = world.plates[0]
        plate.set_fields_on_plate(**{field: np.full(plate.node_count(), cap)})

    log = ReplayLog(tmp_path, step_progress=_fake_step(to_cap))
    _run(log, world)
    log.flush()
    [record] = log.records()
    assert trigger in record["triggers"]
    sidecar = json.loads((tmp_path / record["name"] / replay_log.RECORD_FILE).read_text())
    assert sidecar["health_before"][f"{trigger}_fraction"] == 0.0
    assert sidecar["health_after"][f"{trigger}_fraction"] > replay_log.AT_CAP_FRACTION


def test_a_step_that_raises_saves_the_world_before_it(tmp_path, world):
    log = ReplayLog(tmp_path, step_progress=_fake_step(_raise))
    elapsed = world.elapsed_years
    with pytest.raises(RuntimeError, match="boom"):
        _run(log, world)
    log.flush()

    [record] = log.records()
    assert record["triggers"] == ["step_exception"]
    sidecar = json.loads((tmp_path / record["name"] / replay_log.RECORD_FILE).read_text())
    assert sidecar["triggers"]["step_exception"]["message"] == "boom"
    assert sidecar["fingerprint_after"] is None
    before = persistence.load_world_bytes((tmp_path / record["name"] / replay_log.BEFORE_FILE).read_bytes())
    assert before.elapsed_years == elapsed


def test_a_healthy_step_writes_nothing(tmp_path, world):
    log = ReplayLog(tmp_path, step_progress=_fake_step(lambda w: None))
    _run(log, world)
    log.flush()
    assert log.records() == []
    assert log.recent_events() == []


def test_a_replay_reproduces_the_logged_step(tmp_path, world):
    log = ReplayLog(tmp_path, triggers={"always": lambda outcome: {"forced": True}})
    _run(log, world, years=1_000_000.0)
    log.flush()

    [record] = log.records()
    result = replay_log.replay(record["path"])
    assert result.matches, result.mismatched
    assert len(result.world.plates) == len(world.plates)
    assert result.world.craton_ledger == world.craton_ledger


def test_a_replay_that_differs_is_reported(tmp_path, world):
    log = ReplayLog(tmp_path, triggers={"always": lambda outcome: {"forced": True}})
    _run(log, world, years=1_000_000.0)
    log.flush()
    [record] = log.records()
    path = tmp_path / record["name"] / replay_log.RECORD_FILE
    sidecar = json.loads(path.read_text())
    sidecar["fingerprint_after"]["sea_level_m"] = "0" * 16
    sidecar["health_after"]["plate_count"] += 1
    path.write_text(json.dumps(sidecar))

    result = replay_log.replay(record["path"])
    assert not result.matches
    assert result.mismatched == ["sea_level_m"]
    assert not result.plate_count_matches


def test_rate_limit_allows_one_record_per_trigger_per_window(tmp_path, world):
    log = ReplayLog(tmp_path, min_steps_between=3, triggers={"always": lambda outcome: {}}, step_progress=_fake_step(lambda w: None))
    for _ in range(7):
        _run(log, world)
    log.flush()
    # Steps 1, 4 and 7.
    assert len(log.records()) == 3


def test_a_rate_limited_trigger_does_not_hold_back_another(tmp_path, world):
    fired = {"noisy": lambda outcome: {}, "rare": lambda outcome: {} if outcome.before.plate_count == 1 else None}
    log = ReplayLog(tmp_path, min_steps_between=10, triggers=fired, step_progress=_fake_step(lambda w: None))
    _run(log, world)
    world.plates = world.plates[:1]
    _run(log, world)
    log.flush()
    # "noisy" is still cooling down on the second step, so it stays out of that record.
    assert [r["triggers"] for r in log.records()] == [["rare"], ["noisy"]]


def test_another_trigger_firing_does_not_reset_a_cooldown(tmp_path, world):
    # "noisy" fires every step and "other" every other step: "noisy" must still wait out its
    # own window rather than ride along with "other" (PR #498 review).
    log = ReplayLog(
        tmp_path, min_steps_between=3,
        triggers={"noisy": lambda outcome: {}, "other": lambda outcome: {} if log._steps % 2 == 0 else None},
        step_progress=_fake_step(lambda w: None),
    )
    for _ in range(7):
        _run(log, world)
    log.flush()
    noisy = [r["elapsed_years_before"] for r in log.records() if "noisy" in r["triggers"]]
    gaps = [a - b for a, b in zip(noisy, noisy[1:])]
    assert gaps and all(gap >= 3 * 1_000.0 for gap in gaps)


def test_disk_cap_removes_the_oldest_records(tmp_path, world, saved_world):
    log = ReplayLog(
        tmp_path, max_bytes=int(2.5 * len(saved_world)), min_steps_between=1,
        triggers={"always": lambda outcome: {}}, step_progress=_fake_step(lambda w: None),
    )
    for _ in range(4):
        _run(log, world)
    log.flush()
    records = log.records()
    assert len(records) == 2
    # The newest two survive.
    assert [r["elapsed_years_before"] for r in records] == [world.elapsed_years - 1_000.0, world.elapsed_years - 2_000.0]


def test_the_newest_record_stays_even_over_the_cap(tmp_path, world):
    log = ReplayLog(tmp_path, max_bytes=1, triggers={"always": lambda outcome: {}}, step_progress=_fake_step(lambda w: None))
    _run(log, world)
    log.flush()
    assert len(log.records()) == 1


def test_a_broken_trigger_is_recorded_and_does_not_stop_the_step(tmp_path, world):
    def broken(outcome):
        raise ValueError("bad trigger")

    log = ReplayLog(tmp_path, triggers={"broken": broken}, step_progress=_fake_step(lambda w: None))
    elapsed = world.elapsed_years
    _run(log, world)
    log.flush()
    assert world.elapsed_years == elapsed + 1_000.0
    [record] = log.records()
    assert record["triggers"] == ["broken"]


def test_annotations_from_inside_the_step_reach_the_record(tmp_path, world):
    log = ReplayLog(
        tmp_path, triggers={"always": lambda outcome: {}},
        step_progress=_fake_step(lambda w: replay_log.annotate("merge", {"keep": 3, "absorb": 7})),
    )
    _run(log, world)
    log.flush()
    [record] = log.records()
    sidecar = json.loads((tmp_path / record["name"] / replay_log.RECORD_FILE).read_text())
    assert sidecar["annotations"] == {"merge": {"keep": 3, "absorb": 7}}
    # Outside a logged step it's a no-op.
    replay_log.annotate("ignored", 1)


def test_leftover_partial_writes_are_removed(tmp_path):
    partial = tmp_path / ".replay-log-abc"
    partial.mkdir()
    (partial / replay_log.BEFORE_FILE).write_bytes(b"cut off")
    log = ReplayLog(tmp_path)
    assert not partial.exists()
    assert log.records() == []


def test_from_env(monkeypatch, tmp_path):
    monkeypatch.delenv(replay_log.REPLAY_LOG_DIR_ENV, raising=False)
    assert ReplayLog.from_env() is None
    monkeypatch.setenv(replay_log.REPLAY_LOG_DIR_ENV, str(tmp_path / "logs"))
    log = ReplayLog.from_env()
    assert log.directory == tmp_path / "logs"
    assert log.directory.is_dir()


def test_fingerprint_ignores_wall_clock_timings_but_not_state(world):
    base = replay_log.fingerprint(world)
    world.collision_polarity_stats["prepass_seconds"] = 123.0
    assert replay_log.fingerprint(world) == base
    world.sea_level_m += 1.0
    changed = replay_log.fingerprint(world)
    assert [name for name in base if base[name] != changed[name]] == ["sea_level_m"]


def test_fingerprint_survives_a_save_and_load(world):
    data = persistence.save_record_bytes(world)
    assert replay_log.fingerprint(persistence.load_world_bytes(data)) == replay_log.fingerprint(world)


# --- Trigger thresholds on synthetic metrics -----------------------------------------------


def _health(**overrides) -> Health:
    values = dict(
        plate_count=10,
        craton_live_m3=1e18,
        craton_ledger={key: 0.0 for key in cratons.CRATON_SINK_ACCOUNTS},
        craton_balance_error_m3=0.0,
        continental_ledger_error=None,
        genuine_continental_area_m2=1e14,
        land_fraction=0.3,
    )
    values.update(overrides)
    return Health(**values)


def _fires(name, before, after):
    return replay_log.TRIGGERS[name](StepOutcome(1e5, before, after)) is not None


def test_craton_step_loss_threshold():
    ledger = {key: 0.0 for key in cratons.CRATON_SINK_ACCOUNTS}
    assert not _fires("craton_step_loss", _health(), _health(craton_ledger={**ledger, "eroded_m3": 0.009e18}))
    assert _fires("craton_step_loss", _health(), _health(craton_ledger={**ledger, "eroded_m3": 0.011e18}))


def test_craton_balance_fires_on_the_change_not_the_level():
    assert not _fires("craton_ledger_balance", _health(craton_balance_error_m3=5e12), _health(craton_balance_error_m3=5e12))
    assert _fires("craton_ledger_balance", _health(), _health(craton_balance_error_m3=2e9))


def test_continental_ledger_fires_only_when_it_opens():
    assert _fires("continental_ledger_open", _health(), _health(continental_ledger_error="open"))
    assert not _fires("continental_ledger_open", _health(continental_ledger_error="open"), _health(continental_ledger_error="open"))


@pytest.mark.parametrize("trigger", ["hc_at_max", "hm_at_max", "glacier_at_max"])
def test_at_cap_triggers_fire_on_the_crossing_not_the_level(trigger):
    metric = f"{trigger}_fraction"
    assert _fires(trigger, _health(**{metric: 0.0}), _health(**{metric: 0.0006}))
    assert _fires(trigger, _health(**{metric: 0.0005}), _health(**{metric: 0.0006}))
    assert not _fires(trigger, _health(**{metric: 0.0}), _health(**{metric: 0.0005}))
    # Already over the threshold before the step: not this step's doing.
    assert not _fires(trigger, _health(**{metric: 0.0006}), _health(**{metric: 0.01}))


def test_jump_thresholds():
    assert not _fires("continental_area_jump", _health(), _health(genuine_continental_area_m2=1.04e14))
    assert _fires("continental_area_jump", _health(), _health(genuine_continental_area_m2=0.9e14))
    assert not _fires("land_fraction_jump", _health(), _health(land_fraction=0.34))
    assert _fires("land_fraction_jump", _health(), _health(land_fraction=0.36))
    assert not _fires("plate_count_jump", _health(), _health(plate_count=8))
    assert _fires("plate_count_jump", _health(), _health(plate_count=13))


# --- Review fixes (#498) -------------------------------------------------------------------


def test_pruning_never_touches_directories_the_log_did_not_write(tmp_path, world):
    # The log may share a directory with the user's own captures.
    mine = tmp_path / "issue452"
    mine.mkdir()
    (mine / "capture.pkl").write_bytes(b"x" * 1000)
    log = ReplayLog(tmp_path, max_bytes=1, min_steps_between=1, triggers={"always": lambda outcome: {}}, step_progress=_fake_step(lambda w: None))
    _run(log, world)
    _run(log, world)
    log.flush()
    assert (mine / "capture.pkl").exists()
    assert len(log.records()) == 1


def test_annotations_survive_a_step_driven_from_several_threads(tmp_path, world):
    # /world/step's StreamingResponse runs each next() on whichever threadpool worker is free.
    import threading

    def two_part_step(world, years):
        replay_log.annotate("first", 1)
        yield 0.5
        replay_log.annotate("second", 2)
        world.elapsed_years += years
        yield 1.0

    log = ReplayLog(tmp_path, triggers={"always": lambda outcome: {}}, step_progress=two_part_step)
    progress = log.step(world, 1_000.0)
    for _ in range(3):
        thread = threading.Thread(target=lambda: next(progress, None))
        thread.start()
        thread.join()
    log.flush()
    [record] = log.records()
    sidecar = json.loads((tmp_path / record["name"] / replay_log.RECORD_FILE).read_text())
    assert sidecar["annotations"] == {"first": 1, "second": 2}
    assert replay_log._active_annotations is None


def test_a_logged_crash_reproduces_only_with_the_same_exception(tmp_path, world, monkeypatch):
    log = ReplayLog(tmp_path, step_progress=_fake_step(_raise))
    with pytest.raises(RuntimeError):
        _run(log, world)
    log.flush()
    [record] = log.records()

    def raises(exc):
        def step(world, years):
            yield from ()
            raise exc
        return step

    monkeypatch.setattr(replay_log, "step_world_progress", raises(RuntimeError("boom")))
    assert replay_log.replay(record["path"]).matches
    monkeypatch.setattr(replay_log, "step_world_progress", raises(KeyError("other")))
    assert not replay_log.replay(record["path"]).matches
    # Same type, different failure.
    monkeypatch.setattr(replay_log, "step_world_progress", raises(RuntimeError("elsewhere")))
    assert not replay_log.replay(record["path"]).matches


def test_records_past_the_queue_bound_are_dropped(tmp_path, world, monkeypatch):
    log = ReplayLog(tmp_path, min_steps_between=1, triggers={"always": lambda outcome: {}}, step_progress=_fake_step(lambda w: None))
    # No writer: everything queued stays queued.
    monkeypatch.setattr(log, "_ensure_writer", lambda: None)
    for _ in range(replay_log.MAX_QUEUED_RECORDS + 2):
        _run(log, world)
    assert log._queue.qsize() == replay_log.MAX_QUEUED_RECORDS
    assert len(log.recent_events()) == replay_log.MAX_QUEUED_RECORDS
    assert not log.flush(timeout_s=0.01)


def test_fingerprint_of_a_set_of_objects_survives_a_save_and_load(world):
    vars(world)["token_set"] = frozenset({(1, 2), (3, 4), "five"})
    reloaded = persistence.load_world_bytes(persistence.save_record_bytes(world))
    assert replay_log.fingerprint(reloaded)["token_set"] == replay_log.fingerprint(world)["token_set"]

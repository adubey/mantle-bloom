"""The multi-world registry, sessions and compute gate (#563)."""

import re
import threading
import time

import pytest

from app.world_registry import FairComputeGate, UnknownWorldError, WorldLimitReached, WorldRegistry


def test_create_list_get_remove():
    registry = WorldRegistry(max_worlds=4)
    a = registry.create(name="A")
    b = registry.create()
    assert a.world_id != b.world_id
    # Unnamed worlds default to an exoplanet-style name (#597), not "Untitled world".
    assert re.fullmatch(r".+ [b-g]", b.name), b.name
    assert registry.get(a.world_id) is a
    assert [s.world_id for s in registry.list()] == [a.world_id, b.world_id]
    assert registry.remove(a.world_id) is a
    with pytest.raises(UnknownWorldError):
        registry.get(a.world_id)
    with pytest.raises(UnknownWorldError):
        registry.remove(a.world_id)


def test_cancelled_ticket_is_skipped_not_wedged():
    # A waiter that times out cancels its ticket; the next release must skip it
    # so later acquires proceed (#563 review F10: verified no stall, kept as a guard).
    gate = FairComputeGate()
    assert gate.acquire()
    outcome: list[bool] = []

    def waiter() -> None:
        outcome.append(gate.acquire(timeout=0.05))

    thread = threading.Thread(target=waiter)
    thread.start()
    thread.join(timeout=10)
    assert not thread.is_alive()
    assert outcome == [False]
    gate.release()
    assert gate.acquire(timeout=0)
    gate.release()


def test_cap_blocks_new_records_but_adopt_restores_past_it():
    registry = WorldRegistry(max_worlds=1)
    registry.create()
    with pytest.raises(WorldLimitReached):
        registry.create()
    from app.world_registry import WorldSession

    registry.adopt(WorldSession(world_id="restored"))
    assert len(registry) == 2


def test_record_round_trip_fields():
    registry = WorldRegistry(max_worlds=2)
    session = registry.create(name="Kepler", is_public=True)
    record = session.record()
    assert record == {
        "world_id": session.world_id,
        "name": "Kepler",
        "is_public": True,
        "created_at": session.created_at,
        "updated_at": session.updated_at,
    }


def test_session_operation_busy_flag():
    registry = WorldRegistry(max_worlds=1)
    session = registry.create()
    assert session.begin_operation() is True
    assert session.begin_operation() is False
    session.end_operation()
    assert session.begin_operation() is True


def test_gate_serves_in_fifo_order():
    gate = FairComputeGate()
    order = []

    def worker(name, delay):
        time.sleep(delay)
        with gate.hold():
            order.append(name)

    # The later arrivals start first at the gate; FIFO tickets must still
    # serve them in arrival (ticket) order, not wakeup order.
    threads = [threading.Thread(target=worker, args=(name, delay)) for name, delay in [("a", 0.0), ("b", 0.05), ("c", 0.1)]]
    with gate.hold():
        for t in threads:
            t.start()
        time.sleep(0.3)
    for t in threads:
        t.join(5)
    assert order == ["a", "b", "c"]


def test_gate_timeout_cancels_ticket():
    gate = FairComputeGate()
    assert gate.acquire(timeout=0) is True
    assert gate.acquire(timeout=0.05) is False
    gate.release()
    # The cancelled ticket is skipped: the next acquire succeeds at once.
    assert gate.acquire(timeout=0) is True
    gate.release()

"""Plume step-partition replay agreement (RFC 0007, #624).

Slow, full-simulation checks live here rather than in backend/unit_tests/ --
see bin/stress_test.sh. The 100 Myr pinned replay below is a proxy for the
RFC's full-replay box: it exercises the same schedule/emplacement code path
at both partitions, and the 2% bar absorbs `deform` path-dependence across
partitions (which is why it is 2% here but 1e-9 on the shared-schedule unit
test). Single seed 7: the schedule is a pure function of seed.
"""

import pytest

from app import debug_worlds, plumes
from app import world as world_mod


def test_replay_totals_agree_across_partitions():
    world_fine = debug_worlds.scenario_two_plate_convergent(7)
    world_coarse = debug_worlds.scenario_two_plate_convergent(7)
    fine_head = fine_tail = coarse_head = coarse_tail = 0.0
    original = plumes.apply

    def apply(world, years, *args, **kwargs):
        nonlocal fine_head, fine_tail, coarse_head, coarse_tail
        out = original(world, years, *args, **kwargs)
        if world is world_fine:
            fine_head += out["head_volume_m3"]
            fine_tail += out["tail_volume_m3"]
        else:
            coarse_head += out["head_volume_m3"]
            coarse_tail += out["tail_volume_m3"]
        return out

    plumes.apply = apply
    try:
        t = 0.0
        while t < 100e6 - 1e-9:
            world_mod.step_world(world_fine, 1e6)
            if int(t / 1e6) % 10 == 0:
                world_mod.step_world(world_coarse, 10e6)
            t += 1e6
    finally:
        plumes.apply = original

    assert fine_head == pytest.approx(coarse_head, rel=0.02)
    assert fine_tail == pytest.approx(coarse_tail, rel=0.02)

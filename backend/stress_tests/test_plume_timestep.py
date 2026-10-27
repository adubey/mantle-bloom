"""Plume step-partition replay agreement (RFC 0007, #624).

Slow, full-simulation checks live here rather than in backend/unit_tests/ --
see bin/stress_test.sh. The 100 Myr pinned replay below is a proxy for the
RFC's full-replay box: it exercises the same schedule/emplacement code path
at both partitions, and the 2% bar absorbs `deform` path-dependence across
partitions (which is why it is 2% here but 1e-9 on the shared-schedule unit
test). Single seed 7: the schedule is a pure function of seed.
"""

import math

import numpy as np
import pytest

from app import debug_worlds, geometry, plumes
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


def _track_geometry(seed, monkeypatch):
    """(angle_deg, age_corr, n) of the injected tail's track on the host plate.

    Per-step deposition maxima are un-rotated into the t=0 frame with the
    host plate's own pinned omega; the track direction is the youngest-third
    centroid minus the oldest-third centroid. Observed on seed 7: 7.0
    degrees, correlation 0.99."""
    world = debug_worlds.scenario_two_plate_convergent(seed)
    # Same silenced-stream convention as unit_tests.test_plumes (#579): the
    # only tail is the injected one. 1e18 keeps the state JSON-safe.
    world.plume_stream = {"head": [10**9, 1e18], "tail": [10**9, 1e18]}
    position = geometry.latlon_to_xyz(math.radians(10.0), math.radians(-30.0))
    plates = [plate for plate in world.plates if plate.node_count() > 0]
    continental, host = plumes._center_host(plates, np.asarray(position))
    assert not continental
    omega = np.asarray(world.pinned_omegas[host], dtype=float)
    tail = plumes.MantlePlume(
        plume_id=-(10**9),
        position=tuple(float(v) for v in position),
        born_years=0.0,
        head_volume_m3=0.0,
        head_remaining_m3=0.0,
        head_end_years=0.0,
        tail_flux_m3_per_yr=2e8,
        tail_end_years=80e6,
    )
    world.mantle_plumes.append(tail)

    steps = []
    original = plumes._deposit_on_plate

    def tap(world_, plate, fields, mask, weights, norm, volume_m3, stage):
        out = original(world_, plate, fields, mask, weights, norm, volume_m3, stage)
        if stage == "tail" and out > 0.0 and plate.plate_id == host:
            xyz = np.asarray(plumes._unit_vectors(plate), dtype=float)[mask]
            share = np.asarray(weights, dtype=float)
            centroid = (xyz * (share / share.sum())[:, None]).sum(axis=0)
            steps.append((float(world_.elapsed_years), centroid / np.linalg.norm(centroid)))
        return out

    monkeypatch.setattr(plumes, "_deposit_on_plate", tap)
    for _ in range(30):
        world_mod.step_world(world, 1e6)
    # Slack is for mid-run topology changes: if the host plate splits or
    # merges, its plate_id can vanish and later steps record elsewhere.
    assert len(steps) >= 20

    plume = np.asarray(position, dtype=float)
    plume = plume / np.linalg.norm(plume)
    axis = omega / np.linalg.norm(omega)
    motion = np.cross(axis, plume)
    motion = motion / np.linalg.norm(motion)
    start_frame = [
        (year, geometry.rotation_matrix_from_omega(-omega, year) @ centroid)
        for year, centroid in steps
    ]
    ordered = sorted(start_frame, key=lambda entry: entry[0])
    third = max(1, len(ordered) // 3)
    old = np.mean([point for _, point in ordered[:third]], axis=0)
    young = np.mean([point for _, point in ordered[-third:]], axis=0)
    direction = young - old
    direction = direction / np.linalg.norm(direction)
    angle = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(direction, -motion))))))
    years = np.array([year for year, _ in ordered])
    upstream = np.array([(point - plume) @ -motion for _, point in ordered])
    return angle, float(np.corrcoef(years, upstream)[0, 1]), len(ordered)


def test_hotspot_track_aligns_with_plate_motion(monkeypatch):
    """RFC 0007 (#666): on a pinned-motion world, a tail's per-step deposition
    maxima form a track within 20 degrees of the plate's motion direction at
    the plume and get younger toward it. Two seeds so the cutoffs are not
    tuned to one noise realization."""
    for seed in (7, 11):
        angle, age_corr, count = _track_geometry(seed, monkeypatch)
        assert count >= 20
        assert angle < 20.0, f"seed {seed}: track {angle:.1f} deg off motion"
        assert age_corr > 0.9, f"seed {seed}: age correlation {age_corr:.3f}"

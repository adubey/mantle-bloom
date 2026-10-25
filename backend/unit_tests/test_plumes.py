import math

import numpy as np
import pytest

from app import lithosphere, mantle, persistence, plumes
from app.world import World
from .quad_fixtures import globe_plate, node_points, quad_plate


def _median_head_m3() -> float:
    return (
        math.sqrt(
            plumes.PLUME_HEAD_VOLUME_MIN_KM3 * plumes.PLUME_HEAD_VOLUME_MAX_KM3
        )
        * plumes.KM3_TO_M3
    )


def _peak_relief_m(volume_m3: float) -> float:
    """Relief of the footprint peak over reference oceanic lithosphere (#577)."""
    hc0 = lithosphere.REFERENCE_HC_OCEANIC_M
    hm = lithosphere.REFERENCE_HM_OCEANIC_M
    rho = lithosphere.crust_density("oceanic")
    before = float(lithosphere.isostatic_elevation(hc0, hm, rho))
    peak_excess = volume_m3 / plumes.head_footprint_area_m2(volume_m3)
    after = float(lithosphere.isostatic_elevation(hc0 + peak_excess, hm, rho))
    return after - before


def test_median_head_gives_plateau_relief():
    """Isolated-head scale check: a median head makes 2.0-2.6 km of relief (#577)."""
    volume = _median_head_m3()
    area = plumes.head_footprint_area_m2(volume)
    assert (
        plumes.PLATEAU_AREA_MIN_KM2 * plumes.KM2_TO_M2
        < area
        < plumes.PLATEAU_AREA_MAX_KM2 * plumes.KM2_TO_M2
    )
    assert 2_000.0 <= _peak_relief_m(volume) <= 2_600.0


def test_min_and_max_heads_stay_within_relief_bounds():
    """Small heads are area-clipped, Ontong-scale heads capped, both sane (#577)."""
    for volume_km3 in (
        plumes.PLUME_HEAD_VOLUME_MIN_KM3,
        plumes.PLUME_HEAD_VOLUME_MAX_KM3,
    ):
        assert 1_200.0 <= _peak_relief_m(volume_km3 * plumes.KM3_TO_M3) <= 4_500.0


def test_kernel_places_scheduled_volume():
    """The flat-topped kernel integrates to the head volume (#577)."""
    volume = _median_head_m3()
    area = plumes.head_footprint_area_m2(volume)
    radius = plumes.head_kernel_radius_m(area)
    edges = np.linspace(0.0, 4.0 * radius, 4001)
    mid = 0.5 * (edges[:-1] + edges[1:])
    ring_area = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    placed = float(np.sum(plumes.head_excess_thickness_m(mid, volume) * ring_area))
    assert placed == pytest.approx(volume, rel=5e-3)
    assert float(plumes.head_excess_thickness_m(np.zeros(1), volume)[0]) == volume / area


def _test_centers() -> list:
    return [
        mantle.ConvectionCenter(position=np.array([0.0, 0.0, 1.0]), strength=2.0, falloff=0.6),
        mantle.ConvectionCenter(position=np.array([0.0, 0.0, -1.0]), strength=1.0, falloff=1.0),
    ]


def _event_key(plume: plumes.MantlePlume) -> tuple:
    return (
        plume.born_years,
        tuple(plume.position),
        plume.head_volume_m3,
        plume.tail_flux_m3_per_yr,
    )


def test_event_stream_identical_across_partitions():
    """Same seed gives the same events under any step partition (#579)."""
    centers = _test_centers()
    whole = plumes.consume_arrivals(11, 1.0, False, centers, {}, 0.0, 10e6)
    split_state: dict = {}
    first = plumes.consume_arrivals(11, 1.0, False, centers, split_state, 0.0, 5e6)
    second = plumes.consume_arrivals(11, 1.0, False, centers, split_state, 5e6, 10e6)
    # Same multiset (streams are consumed head-then-tail per chunk, so chunk order
    # differs; the scheduled events do not).
    assert sorted(_event_key(p) for p in whole) == sorted(_event_key(p) for p in first + second)


def test_lower_multiplier_keeps_subset():
    """Kept events at 0.5 are a subset of those at 1.0 (#579)."""
    centers = _test_centers()
    kept_half = {_event_key(p) for p in plumes.consume_arrivals(11, 0.5, False, centers, {}, 0.0, 50e6)}
    kept_full = {_event_key(p) for p in plumes.consume_arrivals(11, 1.0, False, centers, {}, 0.0, 50e6)}
    assert kept_half
    assert kept_half <= kept_full


def test_cox_modulation_bounded_and_deterministic():
    """Cox variant stays in [0, 3] and is a pure function of time (#579)."""
    assert plumes._modulation(11, "head", 12e6, False) == 1.0
    first = plumes._modulation(11, "head", 12e6, True)
    assert 0.0 <= first <= plumes.PLUME_RATE_MODULATION_MAX
    assert plumes._modulation(11, "head", 12e6, True) == first


def _stream_silenced() -> dict:
    """Stream state whose next arrivals sit past any test window (#579)."""
    return {"head": [10**9, 1e18], "tail": [10**9, 1e18]}


def _head_plume(position: tuple[float, float, float]) -> plumes.MantlePlume:
    volume = 7e15
    return plumes.MantlePlume(
        plume_id=0,
        position=position,
        born_years=2e6,
        head_volume_m3=volume,
        head_remaining_m3=volume,
        head_end_years=5e6,
        tail_flux_m3_per_yr=2e7,
        tail_end_years=85e6,
    )


def _tail_plume(position: tuple[float, float, float]) -> plumes.MantlePlume:
    return plumes.MantlePlume(
        plume_id=1,
        position=position,
        born_years=-10e6,
        tail_flux_m3_per_yr=5e7,
        tail_end_years=70e6,
    )


def _plume_world(position: tuple[float, float, float]) -> World:
    world = World(seed=11, plates=[globe_plate("oceanic")])
    world.plume_stream = _stream_silenced()
    world.mantle_plumes = [_head_plume(position), _tail_plume(position)]
    return world


def _placed_volume(world: World, hc_before: np.ndarray) -> float:
    from app.elevation_lines import line_spacing_rad

    plate = world.plates[0]
    areas = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
    return float(np.sum((plate.collect("crustal_thickness_m") - hc_before) * areas))


def test_delivery_exact_across_partitions():
    """Total delivery matches the schedule within 1e-9 for 1 and 10 Myr steps (#579)."""
    head_volume = 7e15
    tail_volume = 2e7 * 5e6 + 5e7 * 10e6
    position = (1.0, 0.0, 0.0)

    world_a = _plume_world(position)
    hc0_a = world_a.plates[0].collect("crustal_thickness_m").copy()
    result_a = {"head_volume_m3": 0.0, "tail_volume_m3": 0.0}
    for step in range(10):
        world_a.elapsed_years = (step + 1) * 1e6
        out = plumes.apply(world_a, 1e6)
        result_a["head_volume_m3"] += out["head_volume_m3"]
        result_a["tail_volume_m3"] += out["tail_volume_m3"]

    world_b = _plume_world(position)
    hc0_b = world_b.plates[0].collect("crustal_thickness_m").copy()
    world_b.elapsed_years = 10e6
    out_b = plumes.apply(world_b, 10e6)

    assert result_a["head_volume_m3"] == pytest.approx(head_volume, rel=1e-9)
    assert out_b["head_volume_m3"] == pytest.approx(head_volume, rel=1e-9)
    assert result_a["tail_volume_m3"] == pytest.approx(tail_volume, rel=1e-9)
    assert out_b["tail_volume_m3"] == pytest.approx(tail_volume, rel=1e-9)
    assert _placed_volume(world_a, hc0_a) == pytest.approx(head_volume + tail_volume, rel=1e-9)
    assert _placed_volume(world_b, hc0_b) == pytest.approx(head_volume + tail_volume, rel=1e-9)
    assert all(
        p.head_remaining_m3 == pytest.approx(0.0, abs=1e6) for p in world_a.mantle_plumes if not p.is_tail_only
    )
    assert not np.any(world_a.plates[0].collect("is_volcano"))
    assert not np.any(world_b.plates[0].collect("is_volcano"))


def test_hotspot_track_follows_motion():
    """A tail emplaced at successive positions leaves a line along motion (#579)."""
    plate = quad_plate(0, "oceanic", columns=range(16), rows=range(8), density=2.0)
    points = node_points(plate)
    center = points[len(points) // 2]
    center = center / np.linalg.norm(center)
    delta_rad = 400e3 / 6371e3
    cos_d, sin_d = math.cos(delta_rad), math.sin(delta_rad)
    rotated = np.array(
        [
            [cos_d, -sin_d, 0.0],
            [sin_d, cos_d, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    p0 = tuple(float(v) for v in center)
    p1 = tuple(float(v) for v in rotated @ center)

    world = World(seed=11, plates=[plate])
    world.plume_stream = _stream_silenced()
    hc0 = plate.collect("crustal_thickness_m").copy()
    world.mantle_plumes = [_tail_plume(p0)]
    world.elapsed_years = 1e6
    plumes.apply(world, 1e6)
    bump_a = node_points(plate)[int(np.argmax(plate.collect("crustal_thickness_m") - hc0))]
    hc1 = plate.collect("crustal_thickness_m").copy()
    world.mantle_plumes = [_tail_plume(p1)]
    world.elapsed_years = 2e6
    plumes.apply(world, 1e6)
    bump_b = node_points(plate)[int(np.argmax(plate.collect("crustal_thickness_m") - hc1))]

    motion = np.asarray(p1) - np.asarray(p0)
    separation = bump_b - bump_a
    cosang = float(np.dot(separation, motion) / (np.linalg.norm(separation) * np.linalg.norm(motion)))
    assert math.degrees(math.acos(np.clip(cosang, -1.0, 1.0))) < 20.0
    assert np.linalg.norm(bump_b - np.asarray(p1)) < np.linalg.norm(bump_b - np.asarray(p0))


def test_multiplier_zero_is_noop():
    """`plume_multiplier = 0` changes nothing, like a matched plume-off run (#579)."""
    world = World(seed=11, plates=[quad_plate(0, "oceanic")])
    world.plume_multiplier = 0.0
    plate = world.plates[0]
    before = {
        name: plate.collect(name).copy()
        for name in ("crustal_thickness_m", "elevation", "plume_heat_myr", "elev_change_reason")
    }
    ledger_before = dict(world.continental_material_ledger)
    world.elapsed_years = 5e6
    out = plumes.apply(world, 5e6)
    assert out == {"new_heads": 0, "new_tails": 0, "head_volume_m3": 0.0, "tail_volume_m3": 0.0}
    for name, values in before.items():
        assert np.array_equal(plate.collect(name), values)
    assert world.continental_material_ledger == ledger_before
    assert world.mantle_plumes == []
    assert world.plume_stream == {}


def test_old_save_backfills_plume_fields():
    """A save predating plumes loads with an empty population and stream (#579)."""
    world = World(seed=11, plates=[])
    del world.__dict__["mantle_plumes"]
    del world.__dict__["plume_stream"]
    persistence._backfill_added_fields(world)
    assert world.mantle_plumes == []
    assert world.plume_stream == {}

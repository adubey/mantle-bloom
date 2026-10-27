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
    # plume_id first: same-time/same-volume events on different streams would
    # otherwise collide (#624 made the key collision-free on purpose).
    return (
        plume.plume_id,
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
    assert out == {"new_heads": 0, "new_tails": 0, "head_volume_m3": 0.0, "tail_volume_m3": 0.0, "credited_heads": 0}
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


def test_stream_identical_over_many_partitions():
    """Ten chunks and random boundaries schedule the same events (#579)."""
    centers = _test_centers()
    whole = plumes.consume_arrivals(11, 1.0, False, centers, {}, 0.0, 50e6)
    state: dict = {}
    chunks = []
    for stop in range(5, 55, 5):
        chunks.extend(plumes.consume_arrivals(11, 1.0, False, centers, state, (stop - 5) * 1e6, stop * 1e6))
    assert sorted(_event_key(p) for p in whole) == sorted(_event_key(p) for p in chunks)
    rng = np.random.default_rng(0)
    bounds = sorted(rng.uniform(0.0, 50e6, 7).tolist())
    state = {}
    jittered = []
    prev = 0.0
    for bound in bounds + [50e6]:
        jittered.extend(plumes.consume_arrivals(11, 1.0, False, centers, state, prev, bound))
        prev = bound
    assert sorted(_event_key(p) for p in whole) == sorted(_event_key(p) for p in jittered)


def test_cox_stream_identical_across_partitions():
    """The Cox variant is partition-invariant too (#579)."""
    centers = _test_centers()
    whole = plumes.consume_arrivals(11, 1.0, True, centers, {}, 0.0, 30e6)
    state: dict = {}
    first = plumes.consume_arrivals(11, 1.0, True, centers, state, 0.0, 17e6)
    second = plumes.consume_arrivals(11, 1.0, True, centers, state, 17e6, 30e6)
    assert sorted(_event_key(p) for p in whole) == sorted(_event_key(p) for p in first + second)


def test_plume_ids_unique_across_streams():
    """Head ids count up from zero; tail-only ids mirror negative (#579)."""
    centers = _test_centers()
    events = plumes.consume_arrivals(11, 1.0, False, centers, {}, 0.0, 200e6)
    heads = [p.plume_id for p in events if not p.is_tail_only]
    tails = [p.plume_id for p in events if p.is_tail_only]
    assert heads and tails
    assert all(i >= 0 for i in heads)
    assert all(i < 0 for i in tails)
    assert len(set(heads)) == len(heads) and len(set(tails)) == len(tails)


def _straddling_world() -> tuple[World, tuple[float, float, float]]:
    left = quad_plate(0, "oceanic", columns=range(8), rows=range(8), density=2.0)
    right = quad_plate(1, "oceanic", columns=range(8, 16), rows=range(8), density=2.0)
    points = node_points(left)
    position = tuple(float(v) for v in points[3 * 8 + 7] / np.linalg.norm(points[3 * 8 + 7]))
    world = World(seed=11, plates=[left, right])
    world.plume_stream = _stream_silenced()
    world.mantle_plumes = [_head_plume(position)]
    return world, position


def test_cross_plate_footprint_places_exact_total():
    """A boundary-straddling head still places its volume exactly once (#579)."""
    world, _ = _straddling_world()
    before = [p.collect("crustal_thickness_m").copy() for p in world.plates]
    world.elapsed_years = 10e6
    out = plumes.apply(world, 10e6)
    from app.elevation_lines import line_spacing_rad

    placed = 0.0
    got_right = 0.0
    for plate, hc0 in zip(world.plates, before):
        areas = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
        delta = float(np.sum((plate.collect("crustal_thickness_m") - hc0) * areas))
        placed += delta
        if plate.plate_id == 1:
            got_right = delta
    assert out["head_volume_m3"] == pytest.approx(7e15, rel=1e-9)
    # The head plume's tail also emplaces over the window (tail [5, 85] Myr).
    assert placed == pytest.approx(7e15 + out["tail_volume_m3"], rel=1e-9)
    assert got_right > 0.0


def test_heat_refreshes_and_counts_down():
    """Head-active cells refresh to 20 Myr; heat decays; tails stay cold (#579)."""
    world, position = _straddling_world()
    plate = world.plates[0]
    world.elapsed_years = 3e6
    plumes.apply(world, 1e6)
    heat = plate.collect("plume_heat_myr")
    # Refresh sets 20 Myr, then the step's own 1 Myr passes: end-of-step state is
    # the true remaining time.
    assert heat.max() == pytest.approx(plumes.PLUME_HEAT_MYR - 1.0)
    globe_world = World(seed=11, plates=[globe_plate("oceanic")])
    globe_world.plume_stream = _stream_silenced()
    globe_world.mantle_plumes = [_head_plume((1.0, 0.0, 0.0))]
    globe_world.elapsed_years = 3e6
    plumes.apply(globe_world, 1e6)
    globe_heat = globe_world.plates[0].collect("plume_heat_myr")
    assert globe_heat.max() == pytest.approx(plumes.PLUME_HEAT_MYR - 1.0)
    assert np.any(globe_heat == 0.0)
    world.mantle_plumes = []
    world.elapsed_years = 8e6
    plumes.apply(world, 5e6)
    assert plate.collect("plume_heat_myr").max() == pytest.approx(plumes.PLUME_HEAT_MYR - 6.0)
    tail_world = World(seed=11, plates=[quad_plate(0, "oceanic")])
    tail_world.plume_stream = _stream_silenced()
    tail_world.mantle_plumes = [_tail_plume(position)]
    tail_world.elapsed_years = 1e6
    plumes.apply(tail_world, 1e6)
    assert np.all(tail_world.plates[0].collect("plume_heat_myr") == 0.0)


def test_mid_step_birth_keeps_full_post_birth_halo():
    """A head born partway through the step ends it with PLUME_HEAT_MYR minus
    only the post-birth time, not minus the whole step (RFC 0007, #582)."""
    world = World(seed=11, plates=[globe_plate("oceanic")])
    world.plume_stream = _stream_silenced()
    volume = 7e15
    world.mantle_plumes = [plumes.MantlePlume(
        plume_id=0,
        position=(1.0, 0.0, 0.0),
        born_years=2.5e6,
        head_volume_m3=volume,
        head_remaining_m3=volume,
        head_end_years=5.5e6,
        tail_flux_m3_per_yr=2e7,
        tail_end_years=85.5e6,
    )]
    world.elapsed_years = 3e6
    plumes.apply(world, 1e6)
    assert world.plates[0].collect("plume_heat_myr").max() == pytest.approx(plumes.PLUME_HEAT_MYR - 0.5)


def test_dead_plumes_drop_at_step_end():
    """A plume whose tail ended is gone after the step (#579)."""
    world = World(seed=11, plates=[quad_plate(0, "oceanic")])
    world.plume_stream = _stream_silenced()
    spent = _tail_plume((1.0, 0.0, 0.0))
    spent.tail_end_years = -20e6
    world.mantle_plumes = [spent]
    world.elapsed_years = 0.0
    out = plumes.apply(world, 1e6)
    assert out["tail_volume_m3"] == 0.0
    assert world.mantle_plumes == []


def test_continental_head_books_juvenile():
    """Flood-basalt Hc on continental hosts books juvenile volume (#579)."""
    plate = quad_plate(0, "continental", columns=range(4), rows=range(4), density=2.0)
    points = node_points(plate)
    position = tuple(float(v) for v in points[len(points) // 2] / np.linalg.norm(points[len(points) // 2]))
    world = World(seed=11, plates=[plate])
    world.plume_stream = _stream_silenced()
    world.mantle_plumes = [_head_plume(position)]
    before = world.continental_material_ledger.get("juvenile_additions_m3", 0.0)
    world.elapsed_years = 10e6
    out = plumes.apply(world, 10e6)
    booked = world.continental_material_ledger.get("juvenile_additions_m3", 0.0) - before
    # Head flood basalt plus the head plume's tail, both on continental hosts.
    assert booked == pytest.approx(out["head_volume_m3"] + out["tail_volume_m3"], rel=1e-9)
    assert out["head_volume_m3"] == pytest.approx(7e15, rel=1e-9)


def test_oceanic_head_stamps_plateau_and_books_emplaced():
    """A head on oceanic hosts stamps plateau crust plus its age moment and books
    EMPLACED to the plume-head buoyant ledger; tails stamp nothing (#581)."""
    from app import terranes

    world = World(seed=11, plates=[globe_plate("oceanic")])
    world.plume_stream = _stream_silenced()
    world.mantle_plumes = [_head_plume((1.0, 0.0, 0.0))]
    world.elapsed_years = 10e6
    out = plumes.apply(world, 10e6)
    plate = world.plates[0]
    from app.elevation_lines import line_spacing_rad

    areas = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
    plateau = plate.collect("plateau_crust_m")
    moment = plate.collect("buoyant_crust_time_m_yr")
    live = float(np.sum(plateau * areas))
    assert live == pytest.approx(out["head_volume_m3"], rel=1e-9)
    assert float(np.sum(moment * areas)) == pytest.approx(live * 10e6, rel=1e-9)
    ledger = terranes.ensure_ledger(world)
    assert ledger["plume_head"]["emplaced_m3"] == pytest.approx(out["head_volume_m3"], rel=1e-9)
    # The head plume's own tail emplaces ordinary melt: no tracer, no ledger.
    assert not np.any(plate.collect("plateau_crust_m") < 0.0)
    assert out["tail_volume_m3"] > 0.0
    assert terranes.balance_residual(ledger, terranes.world_live_volumes_m3(world))["plume_head"] == pytest.approx(
        0.0, abs=1.0
    )


def test_oceanic_tail_stamps_no_tracer():
    """Tail-only plumes emplace ordinary-thickness melt with no buoyant tracer (#581)."""
    from app import terranes

    world = World(seed=11, plates=[globe_plate("oceanic")])
    world.plume_stream = _stream_silenced()
    world.mantle_plumes = [_tail_plume((1.0, 0.0, 0.0))]
    world.elapsed_years = 10e6
    out = plumes.apply(world, 10e6)
    assert out["tail_volume_m3"] > 0.0
    assert not np.any(world.plates[0].collect("plateau_crust_m"))
    ledger = terranes.ensure_ledger(world)
    assert ledger["plume_head"]["emplaced_m3"] == 0.0
    assert ledger["plume_tail"]["emplaced_m3"] == 0.0


def test_continental_plume_books_juvenile_by_source():
    """Continental heads book flood-basalt juvenile, tails hotspot juvenile, each
    under the `juvenile_additions_m3` total (#581)."""
    plate = quad_plate(0, "continental", columns=range(4), rows=range(4), density=2.0)
    points = node_points(plate)
    position = tuple(float(v) for v in points[len(points) // 2] / np.linalg.norm(points[len(points) // 2]))
    world = World(seed=11, plates=[plate])
    world.plume_stream = _stream_silenced()
    world.mantle_plumes = [_head_plume(position)]
    world.elapsed_years = 10e6
    out = plumes.apply(world, 10e6)
    ledger = world.continental_material_ledger
    assert ledger["juvenile_flood_basalt_m3"] == pytest.approx(out["head_volume_m3"], rel=1e-9)
    assert ledger["juvenile_hotspot_m3"] == pytest.approx(out["tail_volume_m3"], rel=1e-9)
    assert ledger["juvenile_additions_m3"] == pytest.approx(
        out["head_volume_m3"] + out["tail_volume_m3"], rel=1e-9
    )


def test_head_logs_event_with_id_host_and_volume():
    """Each stream-born head logs id, volume, and host plate (#579)."""
    world = World(seed=11, plates=[globe_plate("oceanic")])
    world.elapsed_years = 200e6
    plumes.apply(world, 200e6)
    heads = [message for _, message in world.events if message.startswith("Mantle plume head")]
    assert heads
    assert all("M km3" in message and "plate" in message for message in heads)


def test_old_plates_default_heat_to_zero():
    """Plates made before the heat field read zero heat (#579)."""
    plate = quad_plate(0, "oceanic")
    assert np.all(plate.collect("plume_heat_myr") == 0.0)


def _head_stress_world(crust_type):
    """Globe world with a live head stream forced to arrive in-window (#582).

    Returns (world, arrival); the globe plate hosts every sampled head, so all
    arrivals are `crust_type` heads."""
    from .quad_fixtures import globe_plate

    world = World(seed=11, plates=[globe_plate(crust_type)])
    world.plume_multiplier = plumes.PLUME_MULTIPLIER_MAX
    arrival = 3e6
    kept = next(
        k
        for k in range(200)
        if plumes._sample_event(seed=11, stream="head", index=k, born_years=arrival,
                                multiplier=world.plume_multiplier, cox=False,
                                centers=world.mantle_centers) is not None
    )
    world.plume_stream = {"head": [kept, arrival], "tail": [10**9, 1e18]}
    world.mantle_plumes = []
    world.elapsed_years = arrival
    return world, arrival


def test_continental_head_adds_plume_stress_to_host():
    """RFC 0007 (#582): every continental head arrival adds PLUME_HEAD_STRESS to
    its host plate's internal_stress; oceanic heads add nothing."""
    world, _ = _head_stress_world("continental")
    out = plumes.apply(world, 1e6)
    assert out["new_heads"] >= 1
    assert world.plates[0].internal_stress == pytest.approx(out["new_heads"] * plumes.PLUME_HEAD_STRESS)
    assert out["credited_heads"] == out["new_heads"]


def test_oceanic_head_adds_no_plume_stress():
    """RFC 0007 (#582): oceanic heads leave internal_stress untouched."""
    world, _ = _head_stress_world("oceanic")
    out = plumes.apply(world, 1e6)
    assert out["new_heads"] >= 1
    assert world.plates[0].internal_stress == pytest.approx(0.0)
    assert out["credited_heads"] == 0


def test_head_stress_credit_uses_plate_ids():
    """RFC 0007 (#582): the credit follows the host plate id, not the
    plates_data position; oceanic heads and unknown hosts credit nothing."""
    from .quad_fixtures import quad_plate

    oceanic = quad_plate(3, "oceanic")
    continental = quad_plate(7, "continental")
    entries = [{"plate": oceanic}, {"plate": continental}]
    assert plumes._credit_head_stress(entries, True, 7) is True
    assert continental.internal_stress == pytest.approx(plumes.PLUME_HEAD_STRESS)
    assert oceanic.internal_stress == pytest.approx(0.0)
    assert plumes._credit_head_stress(entries, False, 7) is False
    assert plumes._credit_head_stress(entries, True, 999) is False
    assert continental.internal_stress == pytest.approx(plumes.PLUME_HEAD_STRESS)


def test_head_volume_distribution_default_matches_log_uniform_draw_for_draw():
    """RFC 0007 (#584): the default distribution is exactly the #579 draw, so
    scheduled volumes are unchanged unless the knob moves."""
    from app import plumes

    assert plumes.PLUME_VOLUME_DISTRIBUTION == "log_uniform"
    for seed in range(5):
        rng_a = np.random.default_rng(seed)
        rng_b = np.random.default_rng(seed)
        assert plumes._sample_head_volume_km3(rng_a) == plumes._log_uniform(
            rng_b, plumes.PLUME_HEAD_VOLUME_MIN_KM3, plumes.PLUME_HEAD_VOLUME_MAX_KM3
        )


def test_head_volume_distribution_variants_stay_in_range_and_reject_unknown(monkeypatch):
    """RFC 0007 (#584): lognormal and power-law samplers stay inside the shared
    range; the truncated lognormal leaves no pile-up on the clamps; an unknown
    name and a degenerate power-law alpha fail loudly."""
    from app import plumes

    rng = np.random.default_rng(0)
    for name in ("log_uniform", "lognormal", "power_law"):
        monkeypatch.setattr(plumes, "PLUME_VOLUME_DISTRIBUTION", name)
        for _ in range(2000):
            volume = plumes._sample_head_volume_km3(rng)
            assert plumes.PLUME_HEAD_VOLUME_MIN_KM3 < volume < plumes.PLUME_HEAD_VOLUME_MAX_KM3
    monkeypatch.setattr(plumes, "PLUME_VOLUME_DISTRIBUTION", "gaussian")
    with pytest.raises(ValueError, match="unknown plume volume distribution"):
        plumes._sample_head_volume_km3(rng)
    monkeypatch.setattr(plumes, "PLUME_VOLUME_DISTRIBUTION", "power_law")
    monkeypatch.setattr(plumes, "PLUME_VOLUME_POWER_LAW_ALPHA", 1.0)
    with pytest.raises(ValueError, match="degenerate"):
        plumes._sample_head_volume_km3(rng)


def test_diagnostics_counts_population():
    """diagnostics() reports the live census (#579)."""
    world = World(seed=11, plates=[])
    world.mantle_plumes = [_head_plume((1.0, 0.0, 0.0)), _tail_plume((0.0, 1.0, 0.0))]
    census = plumes.diagnostics(world)
    assert census["active_plumes"] == 2
    assert census["heads"] == 1
    assert census["tails"] == 1


def _consume_schedule(seed, step_years, total_years, multiplier=1.0, cox=False):
    """Fresh plumes consumed over [0, total] in `step_years` increments."""
    state: dict = {}
    fresh: list = []
    t = 0.0
    while t < total_years - 1e-9:
        fresh.extend(
            plumes.consume_arrivals(seed, multiplier, cox, [], state, t, t + step_years)
        )
        t += step_years
    return fresh


def test_event_stream_multiset_identical_across_step_partitions():
    """RFC 0007 (#624): the same schedule consumed as 10x1 Myr or 1x10 Myr gives
    the same event multiset. Order legitimately groups per step (stream-major
    per consume call), so compare the file's exact-key convention sorted.
    Single seed 7: the schedule is a pure function of seed, so one seed pins it."""
    fine = _consume_schedule(7, 1e6, 10e6)
    coarse = _consume_schedule(7, 10e6, 10e6)
    assert len(fine) == len(coarse) > 0
    assert sorted(_event_key(p) for p in fine) == sorted(_event_key(p) for p in coarse)


def _pinned_branches(seed, step_years, total_years):
    """plume_id -> oceanic_host per head on a pinned-motion convergent world.

    Pinned means rotation only (`pinned_omegas`): per-step `deform` stays
    path-dependent, which is exactly the branch-flip mechanism this measures.
    Single seed 7: the schedule is a pure function of seed, and the 95% bar
    absorbs one-off flips."""
    from app import debug_worlds
    from app import world as world_mod

    world = debug_worlds.scenario_two_plate_convergent(seed)
    t = 0.0
    while t < total_years - 1e-9:
        world_mod.step_world(world, min(step_years, total_years - t))
        t += step_years
    return {
        plume.plume_id: plume.oceanic_host
        for plume in world.mantle_plumes
        if not plume.is_tail_only
    }


def test_pinned_branch_agrees_across_partitions():
    """RFC 0007 (#624): on pinned motion, the oceanic/continental branch agrees
    on over 95% of events between 1 Myr and 10 Myr partitions."""
    fine = _pinned_branches(7, 1e6, 20e6)
    coarse = _pinned_branches(7, 10e6, 20e6)
    # Drops are partition-independent (same schedule, same drop rule), so the
    # id sets must match exactly -- intersection would hide a divergence.
    assert set(fine) == set(coarse), f"only in fine: {set(fine) - set(coarse)}, only in coarse: {set(coarse) - set(fine)}"
    assert len(fine) > 0
    assert all(host is not None for host in fine.values())
    disagree = [pid for pid in fine if fine[pid] != coarse[pid]]
    agreement = 1.0 - len(disagree) / len(fine)
    assert agreement > 0.95, f"flipped branches: {disagree}"


def _world_fingerprint(world):
    """Per-node max-abs state plus leak-prone scalars.

    Sums can mask compensating writes, so fields compare by max absolute node
    value; heat, stress, ledgers, counters, stream state, and event count are
    included because a heat- or stress-only leak would otherwise be invisible.
    """
    from app import terranes

    material_ledger = getattr(world, "continental_material_ledger", None) or {}
    fingerprint = [world.elapsed_years, len(world.mantle_plumes)]
    fingerprint.append(float(np.nansum(material_ledger.get("juvenile_flood_basalt_m3", 0.0))))
    fingerprint.append(float(np.nansum(material_ledger.get("juvenile_hotspot_m3", 0.0))))
    ledger = terranes.ensure_ledger(world)
    for source in terranes.SOURCES:
        for account in terranes.BALANCE_ACCOUNTS:
            fingerprint.append(float(ledger[source][account]))
    fingerprint.append(world.plume_heads_emplaced)
    fingerprint.append(repr(world.plume_stream))
    fingerprint.append(len(world.events))
    nodes = []
    for plate in world.plates:
        for field in (
            "crustal_thickness_m",
            "elevation",
            "continental_material_m",
            "plume_heat_myr",
            "internal_stress",
        ):
            try:
                nodes.append(float(np.nanmax(np.abs(plate.collect(field)))))
            except Exception:
                continue
    fingerprint.extend(nodes)
    return fingerprint


def test_multiplier_zero_matches_plume_off_replay(monkeypatch):
    """RFC 0007 (#624): `plume_multiplier = 0` is a strict no-op, matching a
    10-step replay with the plume pass patched out."""
    from app import debug_worlds
    from app import world as world_mod

    zero_result = {"new_heads": 0, "new_tails": 0, "head_volume_m3": 0.0,
                   "tail_volume_m3": 0.0, "credited_heads": 0}

    def replay(patch):
        world = debug_worlds.scenario_two_plate_convergent(7)
        world.plume_multiplier = 0.0
        if patch:
            # years=0 takes the real early-return path without touching the
            # world, so this pins the patched dict to the real contract.
            assert set(zero_result) == set(plumes.apply(world, 0.0))
            monkeypatch.setattr(
                plumes, "apply", lambda wd, yrs, *a, **k: dict(zero_result)
            )
        for _ in range(10):
            world_mod.step_world(world, 1e6)
        return _world_fingerprint(world)

    unpatched, patched = replay(False), replay(True)
    assert len(unpatched) == len(patched)
    for got, want in zip(unpatched, patched):
        if isinstance(got, float):
            assert got == pytest.approx(want)
        else:
            assert got == want




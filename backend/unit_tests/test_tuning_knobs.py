"""The live geomorphic-budget tuning knobs (World.*_multiplier, set via /world/controls).

Every knob is a dimensionless multiplier, 1.0 == the model's untuned behaviour. These tests
pin the two properties that matter: 1.0 is a genuine no-op (bit-identical to the knob not
existing), and turning a knob up/down actually moves the corresponding process.
"""

import numpy as np
import pytest

from app import erosion, volcanism
from app.plates import ElevationLine, PlateWithLines
from app.world import TUNING_MULTIPLIER_FIELDS, World, generate_world, step_world


def _elevation_snapshot(world: World) -> np.ndarray:
    _, elevation, _, _, _, _ = erosion._gather_nodes(world)
    return elevation.copy()


def test_all_multipliers_default_to_one():
    world = World(seed=0)
    for name in TUNING_MULTIPLIER_FIELDS:
        assert getattr(world, name) == 1.0


def test_old_pickle_without_the_fields_reads_one():
    # Plain-scalar dataclass defaults are class attributes, so a World pickled before these
    # fields existed (simulated by deleting them from __dict__) still reads 1.0 on load.
    from app import persistence

    world = generate_world(seed=5, num_plates=8)
    for name in TUNING_MULTIPLIER_FIELDS:
        world.__dict__.pop(name, None)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    for name in TUNING_MULTIPLIER_FIELDS:
        assert getattr(loaded, name) == 1.0


def test_erosion_multipliers_at_one_are_bit_identical_to_untuned():
    a = generate_world(seed=21, num_plates=8)
    b = generate_world(seed=21, num_plates=8)
    for name in TUNING_MULTIPLIER_FIELDS:
        setattr(b, name, 1.0)  # explicitly, vs. `a` which never touches them
    erosion.apply_erosion(a, years=5_000_000)
    erosion.apply_erosion(b, years=5_000_000)
    np.testing.assert_array_equal(_elevation_snapshot(a), _elevation_snapshot(b))


def test_zeroing_every_erosion_knob_nearly_freezes_the_surface():
    frozen = generate_world(seed=21, num_plates=8)
    normal = generate_world(seed=21, num_plates=8)
    before = _elevation_snapshot(frozen)
    for name in (
        "rain_erosion_multiplier", "river_erosion_multiplier", "wind_erosion_multiplier",
        "ocean_erosion_multiplier", "coastal_leveling_multiplier", "glacier_erosion_multiplier",
        "seismic_erosion_multiplier", "river_deposition_multiplier", "ocean_deposition_multiplier",
    ):
        setattr(frozen, name, 0.0)

    erosion.apply_erosion(frozen, years=5_000_000)
    erosion.apply_erosion(normal, years=5_000_000)

    frozen_motion = np.abs(_elevation_snapshot(frozen) - before).mean()
    normal_motion = np.abs(_elevation_snapshot(normal) - before).mean()
    assert frozen_motion < 0.2 * normal_motion


@pytest.mark.parametrize("knob", ["rain_erosion_multiplier", "river_erosion_multiplier", "ocean_erosion_multiplier"])
def test_turning_an_erosion_knob_up_erodes_more(knob):
    before = _elevation_snapshot(generate_world(seed=21, num_plates=8))
    low = generate_world(seed=21, num_plates=8)
    high = generate_world(seed=21, num_plates=8)
    setattr(low, knob, 0.0)
    setattr(high, knob, 3.0)
    erosion.apply_erosion(low, years=5_000_000)
    erosion.apply_erosion(high, years=5_000_000)
    low_motion = np.abs(_elevation_snapshot(low) - before).mean()
    high_motion = np.abs(_elevation_snapshot(high) - before).mean()
    assert high_motion > low_motion


def test_ocean_deposition_knob_builds_or_starves_the_shelf():
    before = _elevation_snapshot(generate_world(seed=21, num_plates=8))
    starved = generate_world(seed=21, num_plates=8)
    fed = generate_world(seed=21, num_plates=8)
    starved.ocean_deposition_multiplier = 0.0
    fed.ocean_deposition_multiplier = 3.0
    erosion.apply_erosion(starved, years=5_000_000)
    erosion.apply_erosion(fed, years=5_000_000)
    starved_elev = _elevation_snapshot(starved)
    fed_elev = _elevation_snapshot(fed)
    ocean = before <= 0.0
    # More marine sediment settling -> ocean-side nodes end up higher on average.
    assert fed_elev[ocean].mean() > starved_elev[ocean].mean()


def _all_active_volcano_world(n: int = 200) -> World:
    line = ElevationLine(
        phi=0.0, theta=np.arange(n) * 0.001, elevation=np.full(n, 200.0),
        is_volcano=np.ones(n, dtype=bool), volcano_active_years_remaining=np.full(n, 5_000_000.0),
    )
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    return World(seed=0, plates=[plate])


def test_volcanism_multiplier_zero_stops_all_eruptions():
    world = _all_active_volcano_world()
    world.volcanism_multiplier = 0.0
    before = world.plates[0].lines[0].elevation.copy()
    for _ in range(10):
        volcanism.apply_volcanic_activity(world, years=200_000)
    np.testing.assert_array_equal(world.plates[0].lines[0].elevation, before)


def test_volcanism_multiplier_scales_total_land_built():
    def built(multiplier: float) -> float:
        world = _all_active_volcano_world()
        world.volcanism_multiplier = multiplier
        before = world.plates[0].lines[0].elevation.copy()
        for _ in range(10):
            volcanism.apply_volcanic_activity(world, years=200_000)
        return float((world.plates[0].lines[0].elevation - before).sum())

    assert built(3.0) > built(1.0) > 0.0


def _stepped_crust_and_relief(seed: int, steps: int, **knobs) -> tuple[float, float]:
    from app import plates as plates_mod

    world = generate_world(seed=seed, num_plates=8)
    for name, value in knobs.items():
        setattr(world, name, value)
    for _ in range(steps):
        step_world(world, years=2_000_000)
    hc = plates_mod.collect_all_crustal_thickness(world.plates)
    _, elevation, _, _, _, _ = erosion._gather_nodes(world)
    land = elevation[elevation > 0.0]
    return float(hc.sum()), float(np.percentile(land, 90)) if len(land) else 0.0


def test_collision_uplift_amount_scales_crustal_thickening():
    # A colliding-plate seed, few enough steps that no peak has railed at MAX_ELEVATION_M.
    strong_hc, _ = _stepped_crust_and_relief(559394024, 6, collision_uplift_multiplier=3.0)
    weak_hc, _ = _stepped_crust_and_relief(559394024, 6, collision_uplift_multiplier=0.2)
    assert strong_hc > weak_hc


def test_collision_uplift_amount_at_one_is_untuned():
    tuned_hc, _ = _stepped_crust_and_relief(559394024, 6, collision_uplift_multiplier=1.0,
                                            collision_uplift_reach_multiplier=1.0)
    untuned_hc, _ = _stepped_crust_and_relief(559394024, 6)
    assert tuned_hc == untuned_hc


def test_collision_uplift_reach_widens_the_thickened_belt():
    # Wider reach -> more crust participates in the orogenic thickening -> more total Hc.
    wide_hc, _ = _stepped_crust_and_relief(559394024, 6, collision_uplift_reach_multiplier=3.0)
    narrow_hc, _ = _stepped_crust_and_relief(559394024, 6, collision_uplift_reach_multiplier=1.0)
    assert wide_hc > narrow_hc

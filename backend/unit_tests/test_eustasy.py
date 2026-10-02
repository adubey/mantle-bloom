"""Eustatic sea level (`app.eustasy`): `World.sea_level_m` re-solved each step from a
conserved ocean water volume against the world's hypsometry."""

import numpy as np
import pytest

from app import eustasy
from app.world import generate_world as _generate_world, step_world


def generate_world(*args, **kwargs):
    kwargs.setdefault("surface", "lines")
    return _generate_world(*args, **kwargs)


def test_total_water_volume_is_monotonic_and_area_weighted():
    z = np.array([-4000.0, -1000.0, -200.0, 50.0, 800.0])
    area = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    # Only the nodes below the datum contribute, each its own depth times its own area.
    assert eustasy.total_water_volume_m3(z, area, 0.0) == 4000.0 + 2 * 1000.0 + 3 * 200.0
    assert eustasy.total_water_volume_m3(z, area, -300.0) == 3700.0 + 2 * 700.0
    assert eustasy.total_water_volume_m3(z, area, -5000.0) == 0.0
    # Strictly increasing in sea level.
    levels = np.linspace(-4000.0, 1000.0, 50)
    w = [eustasy.total_water_volume_m3(z, area, h) for h in levels]
    assert np.all(np.diff(w) >= 0.0)


def test_solve_sea_level_inverts_total_water_volume():
    rng = np.random.default_rng(0)
    z = np.sort(rng.uniform(-6000.0, 3000.0, 500))
    area = rng.uniform(0.8e9, 1.2e9, 500)
    for target_level in (-3000.0, -500.0, 0.0, 250.0):
        w = eustasy.total_water_volume_m3(z, area, target_level)
        solved = eustasy.solve_sea_level(z, area, w)
        assert abs(solved - target_level) < 1.0


def test_moving_water_onto_smaller_cells_keeps_its_volume():
    """Issue #257: the budget is a volume, so the same water over cells of a different size
    stands at the level that volume actually reaches -- the node count alone doesn't decide it."""
    z = np.full(4, -100.0)
    big = np.full(4, 2.0e9)
    budget = eustasy.total_water_volume_m3(z, big, 0.0)
    # The same basin retiled into cells half the area: the same water now stands twice as deep.
    assert eustasy.solve_sea_level(z, big / 2.0, budget) == pytest.approx(100.0, abs=1.0)


def test_deepening_the_ocean_basins_drops_sea_level():
    """The whole point: move crust down (a spreading basin) and the same water volume no
    longer reaches as high -- sea level falls, handing dry land back as freeboard."""
    world = generate_world(seed=555, num_plates=8, continental_fraction=0.5, node_density=0.3)
    eustasy.update_sea_level(world)
    sea_before = world.sea_level_m

    for plate in world.plates:
        if plate.crust_type != "oceanic":
            continue
        for i, line in enumerate(plate.lines):
            plate.replace_line(i, line.replace(elevation=line.elevation - 800.0))
    eustasy.update_sea_level(world)

    assert world.sea_level_m < sea_before - 50.0
    # Water volume itself is unchanged -- only where it sits.
    assert eustasy.water_volume_for_sea_level(world, world.sea_level_m) == \
        pytest.approx(world.ocean_water_volume_m3, rel=1e-6)


def test_generate_world_snapshots_the_budget_and_keeps_sea_level_at_zero():
    world = generate_world(seed=7, num_plates=8, continental_fraction=0.5, node_density=0.3)
    assert world.ocean_water_volume_m3 is not None
    assert world.sea_level_m == 0.0  # the flat starting datum
    # Re-solving an unchanged world reproduces the same level (the solve is the exact inverse).
    eustasy.update_sea_level(world)
    assert abs(world.sea_level_m) < 1.0


def test_set_sea_level_via_water_budget_is_conserved_across_a_step():
    world = generate_world(seed=9, num_plates=8, continental_fraction=0.5, node_density=0.3)
    world.simulate_plate_movement = False
    world.simulate_climate_biomes = False

    eustasy.set_sea_level_via_water_budget(world, 120.0)
    assert world.sea_level_m == 120.0
    budget = world.ocean_water_volume_m3

    step_world(world, years=1_000_000)
    # Nothing moved (both sims off), so the level and the budget both hold.
    assert abs(world.sea_level_m - 120.0) < 1.0
    assert world.ocean_water_volume_m3 == budget


class _FakeHydrologyCache:
    """Just enough of `HydrologyFields` for `eustasy._ocean_connected_mask` to work with:
    `points` (only its length is consulted), `neighbor_idx` (the actual connectivity graph),
    and `is_ocean` (the seed `_seeded_connected_mask` anchors to -- see that function's own
    docstring for why a fixed seed, not a fresh "largest wins" pass, is what the eustasy solve
    needs). Lets the connectivity-aware path be exercised without spinning up a real `World`/
    plate tectonics run."""

    def __init__(self, n: int, neighbor_idx: np.ndarray, is_ocean: np.ndarray):
        self.points = np.zeros((n, 3))
        self.neighbor_idx = neighbor_idx
        self.is_ocean = is_ocean


def _ring_neighbors(indices: list[int]) -> list[tuple[int, int]]:
    """Every consecutive (wraparound) pair in `indices` -- builds one connected ring
    component, isolated from any other index range that gets its own separate ring."""
    return [(indices[i], indices[(i + 1) % len(indices)]) for i in range(len(indices))]


def _make_fake_world_with_two_below_sea_level_rings():
    """20 nodes, three disjoint rings, no edges crossing between them: a 10-node "ocean" ring
    at -100 m, a 5-node "closed pit" ring at -500 m (below sea level, but never connected to
    the ocean ring), and a 5-node "land" ring at +50 m -- a minimal fixture for
    `connected_ocean_mask` to draw exactly the distinction this module's whole fix is about:
    the pit is below sea level but isn't the ocean, so its own volume shouldn't count as free
    ocean depth (see eustasy.py's own module docstring). Every node has unit area."""
    n = 20
    elevations = np.concatenate([np.full(10, -100.0), np.full(5, -500.0), np.full(5, 50.0)])
    edges = _ring_neighbors(list(range(0, 10))) + _ring_neighbors(list(range(10, 15))) + _ring_neighbors(list(range(15, 20)))
    neighbor_idx = np.full((n, 2), -1, dtype=np.int64)
    fill = np.zeros(n, dtype=np.int64)
    for a, b in edges:
        neighbor_idx[a, fill[a]] = b
        neighbor_idx[b, fill[b]] = a
        fill[a] += 1
        fill[b] += 1
    # The seed is last step's own real `is_ocean` -- the 10-node ring is the only component
    # ever actually classified ocean; the pit is deliberately left out of it, exactly like a
    # real save's `hydrology_cache.is_ocean` would leave out a genuine closed basin.
    is_ocean = np.zeros(n, dtype=bool)
    is_ocean[:10] = True
    world = type("FakeWorld", (), {})()
    world.hydrology_cache = _FakeHydrologyCache(n, neighbor_idx, is_ocean)
    return world, elevations, np.ones(n)


def test_ocean_water_volume_excludes_a_disconnected_below_sea_level_pit():
    world, elevations, areas = _make_fake_world_with_two_below_sea_level_rings()
    sea_level_m = 0.0

    # The bare formula (what the old code used) double-books the closed pit's own volume.
    assert eustasy.total_water_volume_m3(elevations, areas, sea_level_m) == 10 * 100.0 + 5 * 500.0
    # The connectivity-aware one counts only the actual (connected) ocean ring.
    assert eustasy.ocean_water_volume_m3(world, elevations, areas, sea_level_m) == 10 * 100.0


def test_solve_sea_level_connected_ignores_the_pits_volume():
    world, elevations, areas = _make_fake_world_with_two_below_sea_level_rings()
    # A water budget that would exactly float the ocean ring at -20 m, ignoring the pit.
    target_level = -20.0
    ocean_only_budget = 10 * (target_level - (-100.0))
    solved = eustasy._solve_sea_level_connected(world, elevations, areas, ocean_only_budget)
    assert abs(solved - target_level) < 1.0
    # The same budget fed to the bare (connectivity-oblivious) solver reads as a much lower
    # level, since it also expects to have filled the pit "for free" along the way.
    naive = eustasy.solve_sea_level(elevations, areas, ocean_only_budget)
    assert naive < solved - 1.0


def _pile_ice_everywhere(world, depth_m):
    for plate in world.plates:
        for i, line in enumerate(plate.lines):
            gd = line.glacier_depth.copy()
            gd[:] = depth_m
            plate.replace_line(i, line.replace(glacier_depth=gd))


def test_trapped_ice_lowers_sea_level_and_total_budget_is_conserved():
    """Water frozen into ice caps/glaciers is debited from the ocean's share of the conserved
    total budget -- so the shoreline falls as the ice grows (glacio-eustasy), while the total
    surface-water budget itself doesn't change."""
    world = generate_world(seed=321, num_plates=8, continental_fraction=0.5, node_density=0.3)
    eustasy.update_sea_level(world)
    sea_before = world.sea_level_m
    total_before = world.ocean_water_volume_m3
    assert eustasy.trapped_water_volume_m3(world) == 0.0

    _pile_ice_everywhere(world, 400.0)
    eustasy.update_sea_level(world)

    assert eustasy.trapped_water_volume_m3(world) > 0.0
    assert world.sea_level_m < sea_before - 50.0  # ice locked up ocean water -> lower stand
    assert world.ocean_water_volume_m3 == total_before  # the *total* budget is untouched
    # The ocean's own volume plus what's trapped still sums back to the conserved total.
    assert eustasy.water_volume_for_sea_level(world, world.sea_level_m) + eustasy.trapped_water_volume_m3(world) == \
        pytest.approx(world.ocean_water_volume_m3, rel=1e-6)

    # Melt it all back off -> sea level returns to where it started.
    _pile_ice_everywhere(world, 0.0)
    eustasy.update_sea_level(world)
    assert abs(world.sea_level_m - sea_before) < 1.0


def test_set_sea_level_via_water_budget_folds_in_trapped_water():
    world = generate_world(seed=654, num_plates=8, continental_fraction=0.5, node_density=0.3)
    _pile_ice_everywhere(world, 300.0)

    eustasy.set_sea_level_via_water_budget(world, 0.0)
    trapped = eustasy.trapped_water_volume_m3(world)
    assert trapped > 0.0
    # The stored budget is ocean-at-0 *plus* the currently-trapped water.
    assert world.ocean_water_volume_m3 == pytest.approx(
        eustasy.water_volume_for_sea_level(world, 0.0) + trapped, rel=1e-9
    )
    # Re-solving against that budget reproduces the requested level (ice unchanged).
    eustasy.update_sea_level(world)
    assert abs(world.sea_level_m) < 1.0

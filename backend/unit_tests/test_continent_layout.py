import numpy as np
import pytest
from scipy.sparse.csgraph import dijkstra

from app import continent_layout as cl
from app import plates as plates_mod
from app.elevation_lines import PLANET_RADIUS_KM, line_spacing_rad
from app.noise import SphereNoise
from app.lithosphere_plate import cratons_first_layout

SPACING = line_spacing_rad(1.0)


@pytest.fixture(scope="module")
def layout() -> cl.ContinentLayout:
    return cratons_first_layout(7)[0]


def _links(lat: cl.Lattice) -> tuple[np.ndarray, np.ndarray]:
    return np.repeat(np.arange(lat.size), np.diff(lat.adj_offsets)), lat.adj_neighbours


def test_every_point_has_one_owner_and_every_plate_is_connected(layout):
    assert layout.plate.shape == (layout.lattice.size,)
    assert set(np.unique(layout.plate)) == set(range(layout.num_plates))
    assert cl._label_components(layout.lattice, layout.plate) == layout.num_plates


def test_cores_are_indivisible(layout):
    for core_id in np.unique(layout.core[layout.is_core]):
        assert len(np.unique(layout.plate[layout.core == core_id])) == 1


def test_path_distances_through_a_core_stay_geodesic(layout):
    """Dijkstra crosses a core as ordinary lattice, not as one contracted node, so path length
    over the 8-connected lattice stays close to the great-circle distance: the octile metric's
    8% at worst on a square grid, a little more where cube-sphere cells are skewed."""
    lat = layout.lattice
    core_id = int(np.bincount(layout.core[layout.is_core]).argmax())
    members = np.flatnonzero(layout.core == core_id)
    start = int(lat.tree.query(layout.nuclei[core_id])[1])
    assert layout.core[start] == core_id
    path_km = dijkstra(lat.graph(ring=True), directed=False, indices=start)[members]
    great_km = np.arccos(np.clip(lat.xyz[members] @ lat.xyz[start], -1.0, 1.0)) * PLANET_RADIUS_KM
    far = great_km > 500.0
    assert np.any(far)
    ratio = path_km[far] / great_km[far]
    assert ratio.min() >= 0.99
    assert np.median(ratio) <= 1.06
    assert ratio.max() <= 1.15


def test_each_plate_stays_in_its_area_band_and_above_the_node_floor(layout):
    area = layout.lattice.area_km2
    start = np.bincount(layout.initial_plate, weights=area, minlength=layout.num_plates)
    final = np.bincount(layout.plate, weights=area, minlength=layout.num_plates)
    assert np.all(final >= cl.PLATE_AREA_BAND[0] * start - 1e-6)
    assert np.all(final <= cl.PLATE_AREA_BAND[1] * start + 1e-6)
    assert np.bincount(layout.plate).min() >= cl._node_floor(1.0)


def test_refinement_never_raises_the_boundary_cost_and_stops_within_its_sweeps(layout):
    costs = np.asarray(layout.refine_costs)
    assert len(costs) - 1 <= cl.REFINE_MAX_SWEEPS
    assert np.all(np.diff(costs) <= 1e-6)
    assert costs[-1] < costs[0]
    final = cl.boundary_cost(layout.lattice, layout.plate, cl.point_costs(layout.continent, layout.core, layout.belt))
    assert final == pytest.approx(costs[-1], rel=1e-9)
    # Cores never move.
    assert np.array_equal(layout.plate[layout.is_core], layout.initial_plate[layout.is_core])


def _cross_at(lat: cl.Lattice, p: int) -> tuple[list[int], list[int]]:
    """p's edge neighbours, and its diagonal neighbours (the rest of its 1-ring)."""
    edge = lat.adj_neighbours[lat.adj_offsets[p] : lat.adj_offsets[p + 1]].tolist()
    ring = lat.ring_neighbours[lat.ring_offsets[p] : lat.ring_offsets[p + 1]].tolist()
    return edge, [q for q in ring if q not in edge]


def test_simple_point_check_refuses_a_move_that_would_cut_a_plate():
    lat = cl.lattice(SPACING)
    # A face-centre cell: four edge neighbours, four diagonals.
    p = int(lat.tree.query(np.array([1.0, 0.0, 0.0]))[1])
    edge, diagonal = _cross_at(lat, p)
    assert len(edge) == 4 and len(diagonal) == 4
    labels = np.ones(lat.size, dtype=np.int64)
    labels[p] = 0
    # Two opposite edge neighbours of p in plate 0, joined only through p.
    def around(q: int) -> set[int]:
        return set(lat.adj_neighbours[lat.adj_offsets[q] : lat.adj_offsets[q + 1]].tolist()) - {p}

    # The edge neighbour across from edge[0] shares no corner cell with it.
    opposite = [q for q in edge if q != edge[0] and not (around(q) & around(edge[0]))]
    labels[[edge[0], opposite[0]]] = 0
    args = (lat.adj_offsets, lat.adj_neighbours, lat.ring_offsets, lat.ring_neighbours)
    assert not cl._is_simple(p, 0, labels, *args)
    # Join them around the ring through a diagonal and the edge neighbour between: now p is simple.
    side = [q for q in edge if q not in (edge[0], opposite[0])][0]
    corners = [d for d in diagonal if np.any(lat.adj_neighbours[lat.adj_offsets[d] : lat.adj_offsets[d + 1]] == side)]
    labels[side] = 0
    labels[corners] = 0
    assert cl._is_simple(p, 0, labels, *args)


def test_belts_separate_the_cores_with_no_unmarked_gap(layout):
    lat = layout.lattice
    rows, cols = _links(lat)
    a, b = layout.core[rows], layout.core[cols]
    # No link joins two different cores.
    assert not np.any((a >= 0) & (b >= 0) & (a != b))
    # Along the arc between neighbouring nuclei of one continent, everything is core or belt.
    spacing_rad = cl.CRATON_CORE_SPACING_KM / PLANET_RADIUS_KM
    checked = 0
    for i in range(len(layout.nuclei)):
        for j in range(i + 1, len(layout.nuclei)):
            if layout.nucleus_continent[i] != layout.nucleus_continent[j]:
                continue
            if np.arccos(np.clip(layout.nuclei[i] @ layout.nuclei[j], -1, 1)) > 1.3 * spacing_rad:
                continue
            if not (np.any(layout.core == i) and np.any(layout.core == j)):
                continue
            t = np.linspace(0.0, 1.0, 50)[:, None]
            arc = (1 - t) * layout.nuclei[i] + t * layout.nuclei[j]
            arc /= np.linalg.norm(arc, axis=1, keepdims=True)
            idx = layout.nearest(arc)
            assert np.all(layout.is_core[idx] | layout.belt[idx])
            checked += 1
    assert checked > 0
    assert not np.any(layout.belt & ~layout.is_continental)
    assert not np.any(layout.belt & layout.is_core)
    assert np.array_equal(layout.margin, layout.is_continental & ~layout.is_core & ~layout.belt)


def test_land_threshold_reads_the_per_point_continental_mask():
    noise = SphereNoise(np.random.default_rng(3))
    north = plates_mod.land_noise_threshold_for_mask(lambda pts: pts[:, 2] > 0.0, noise, 0.25)
    south = plates_mod.land_noise_threshold_for_mask(lambda pts: pts[:, 2] < 0.0, noise, 0.25)
    assert north != south
    samples = np.random.default_rng(4).normal(size=(20000, 3))
    samples /= np.linalg.norm(samples, axis=1, keepdims=True)
    land = (samples[:, 2] > 0.0) & (noise.sample(samples) > north)
    assert land.mean() == pytest.approx(0.25, abs=0.02)
    assert plates_mod.land_noise_threshold_for_mask(lambda pts: np.zeros(len(pts), dtype=bool), noise, 0.25) is None


def test_continental_and_craton_areas_hit_their_targets(layout):
    area = layout.lattice.area_km2
    sphere = float(area.sum())
    continental = float(area[layout.is_continental].sum())
    assert continental / sphere == pytest.approx(0.5, rel=0.05)
    np.testing.assert_allclose([c["area_km2"] for c in layout.diagnostics()["continents"]], layout.continent_area_km2, rtol=0.05)
    cored = np.isin(layout.continent, np.unique(layout.continent[layout.is_core]))
    assert float(area[layout.is_core].sum()) / float(area[cored].sum()) == pytest.approx(cl.CRATON_SHARE_OF_CONTINENT, rel=0.05)


def test_same_seed_gives_the_same_layout(layout):
    again = cratons_first_layout(7)[0]
    for name in ("continent", "core", "belt", "plate"):
        assert np.array_equal(getattr(again, name), getattr(layout, name))
    assert not np.array_equal(cratons_first_layout(8)[0].plate, layout.plate)


def test_a_continent_too_small_for_a_core_gets_none():
    small = cratons_first_layout(4, continental_fraction=0.05)[0]
    disk = cl._cap_area_km2(cl.CRATON_CORE_RADIUS_KM)
    too_small = np.flatnonzero(small.continent_area_km2 < disk)
    assert len(too_small) > 0
    for k in too_small:
        assert not np.any(small.is_core & (small.continent == k))
        assert np.any(small.continent == k)

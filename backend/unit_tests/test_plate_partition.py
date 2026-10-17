import numpy as np
import pytest
from scipy.sparse.csgraph import dijkstra

from app import plate_partition as pp
from app.continent_layout import ContinentLayout
from app.elevation_lines import PLANET_RADIUS_KM, line_spacing_rad
from app.layout_lattice import Lattice, lattice
from app.lithosphere_plate import cratons_first_layout

SPACING = line_spacing_rad(1.0)


@pytest.fixture(scope="module")
def layout() -> ContinentLayout:
    return cratons_first_layout(7, land_fraction=0.29).layout


def test_every_point_has_one_owner_and_every_plate_is_connected(layout):
    assert layout.plate.shape == (layout.lattice.size,)
    assert set(np.unique(layout.plate)) == set(range(layout.num_plates))
    assert pp.label_components(layout.lattice, layout.plate) == layout.num_plates


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
    assert np.all(final >= pp.PLATE_AREA_BAND[0] * start - 1e-6)
    assert np.all(final <= pp.PLATE_AREA_BAND[1] * start + 1e-6)
    assert np.bincount(layout.plate).min() >= pp.node_floor(1.0)


def test_refinement_never_raises_the_boundary_cost_and_stops_within_its_sweeps(layout):
    costs = np.asarray(layout.refine_costs)
    assert len(costs) - 1 <= pp.REFINE_MAX_SWEEPS
    assert np.all(np.diff(costs) <= 1e-6)
    assert costs[-1] < costs[0]
    final = pp.boundary_cost(layout.lattice, layout.plate, pp.point_costs(layout.continent, layout.core, layout.belt), layout.land)
    assert final == pytest.approx(costs[-1], rel=1e-9)
    # Cores never move.
    assert np.array_equal(layout.plate[layout.is_core], layout.initial_plate[layout.is_core])


def _cross_at(lat: Lattice, p: int) -> tuple[list[int], list[int]]:
    """p's edge neighbours, and its diagonal neighbours (the rest of its 1-ring)."""
    edge = lat.adj_neighbours[lat.adj_offsets[p] : lat.adj_offsets[p + 1]].tolist()
    ring = lat.ring_neighbours[lat.ring_offsets[p] : lat.ring_offsets[p + 1]].tolist()
    return edge, [q for q in ring if q not in edge]


def test_simple_point_check_refuses_a_move_that_would_cut_a_plate():
    lat = lattice(SPACING)
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
    assert not pp._is_simple(p, 0, labels, *args)
    # Join them around the ring through a diagonal and the edge neighbour between: now p is simple.
    side = [q for q in edge if q not in (edge[0], opposite[0])][0]
    corners = [d for d in diagonal if np.any(lat.adj_neighbours[lat.adj_offsets[d] : lat.adj_offsets[d + 1]] == side)]
    labels[side] = 0
    labels[corners] = 0
    assert pp._is_simple(p, 0, labels, *args)


def test_a_land_link_costs_more_only_where_land_is_on_both_sides(layout):
    lat = layout.lattice
    cost = pp.point_costs(layout.continent, layout.core, layout.belt)
    rows, cols = lat.adj_rows, lat.adj_neighbours
    crossing = layout.plate[rows] != layout.plate[cols]
    land_both = crossing & layout.land[rows] & layout.land[cols]
    assert np.any(land_both)
    extra = pp.boundary_cost(lat, layout.plate, cost, layout.land) - pp.boundary_cost(lat, layout.plate, cost)
    # Each link appears once from each end.
    assert extra == pytest.approx(0.5 * pp.PLATE_COST_LAND_LINK * lat.adj_length_km[land_both].sum(), rel=1e-9)


def test_the_kept_partition_is_the_cheapest_candidate(layout, monkeypatch):
    """Drawing more candidates can only lower the starting J: the first draw is the same."""

    def start_cost(candidates: int) -> float:
        monkeypatch.setattr(pp, "PARTITION_CANDIDATES", candidates)
        rng = np.random.default_rng(1)
        partition = pp.partition_plates(rng, layout.lattice, 12, layout.continent, layout.core, layout.belt, 1.0, layout.land)
        return partition.refine_costs[0]

    assert start_cost(4) <= start_cost(1)

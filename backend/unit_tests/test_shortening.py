"""Collisional shortening carried into plate interiors (issue #314): the ring cascade on
synthetic rectangular cell grids."""

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from app import lithosphere, rheology, shortening
from app.mantle import PLANET_RADIUS_KM

DRIVE_PA = 1e8  # a ~3 cm/yr collision's normal-stress proxy


def _grid(rows: int, cols: int, cell_km: float, areas_km2: np.ndarray | None = None):
    """A rows x cols 4-connected grid, column 0 the collision front."""
    index = np.arange(rows * cols).reshape(rows, cols)
    pairs = [(index[:, :-1].ravel(), index[:, 1:].ravel()), (index[:-1, :].ravel(), index[1:, :].ravel())]
    i = np.concatenate([a for a, b in pairs] + [b for a, b in pairs])
    j = np.concatenate([b for a, b in pairs] + [a for a, b in pairs])
    n = rows * cols
    adjacency = csr_matrix((np.ones(len(i), dtype=np.int8), (i, j)), shape=(n, n))
    areas = np.full(n, cell_km**2 * 1e6) if areas_km2 is None else areas_km2.ravel() * 1e6
    col = np.tile(np.arange(cols), rows)
    row = np.repeat(np.arange(rows), cols)
    return adjacency, areas, row, col


def _problem(rows, cols, cell_km=62.5, demand_per_km=1_000.0, craton=None, relief=None, hc=None, hm=None, areas_km2=None):
    adjacency, areas, row, col = _grid(rows, cols, cell_km, areas_km2)
    n = rows * cols
    front = col == 0
    # Metres of convergence per metre of boundary times each front cell's own boundary length.
    demand = np.where(front, demand_per_km * np.sqrt(areas), 0.0)
    problem = shortening.ShorteningProblem(
        adjacency=adjacency,
        areas_m2=areas,
        front=front,
        demand_m2=demand,
        hc_m=np.full(n, lithosphere.REFERENCE_HC_CONTINENTAL_M) if hc is None else hc,
        hm_m=np.full(n, lithosphere.REFERENCE_HM_CONTINENTAL_M) if hm is None else hm,
        craton_strength=np.zeros(n) if craton is None else craton.astype(float),
        relief_m=np.zeros(n) if relief is None else relief,
        crust_density=np.full(n, lithosphere.RHO_CONTINENTAL_CRUST),
        drive_stress_pa=DRIVE_PA,
        spacing_rad=cell_km / PLANET_RADIUS_KM,
    )
    return problem, row, col


def _assert_conserved(problem, result):
    total = float(result.absorbed_m2.sum() + result.returned_m2.sum())
    assert total == pytest.approx(float(problem.demand_m2.sum()), rel=1e-12)
    assert np.all(result.absorbed_m2 >= 0.0) and np.all(result.returned_m2 >= 0.0)


def _by_column(values, col):
    return np.bincount(col, weights=values)


def test_uniform_strip_conserves_shortening_and_fades_inland():
    problem, _, col = _problem(8, 40)
    result = shortening.ring_cascade(problem)

    _assert_conserved(problem, result)
    per_column = _by_column(result.absorbed_m2, col)
    assert np.all(np.diff(per_column[:20]) < 0.0)
    # Weak ground takes up most of it; basal drag returns the rest to the boundary.
    absorbed = result.absorbed_m2.sum() / problem.demand_m2.sum()
    assert 0.7 < absorbed < 0.95
    assert result.returned_decay_m2.sum() > 0.0
    # Routing adds no along-strike artifacts: every row of a uniform strip is the same.
    per_row = result.absorbed_m2.reshape(8, 40)
    np.testing.assert_allclose(per_row, np.broadcast_to(per_row[0], per_row.shape), rtol=1e-9)


def test_a_plate_with_no_collision_is_a_no_op():
    problem, _, _ = _problem(5, 10, demand_per_km=0.0)
    result = shortening.ring_cascade(problem)
    assert not np.any(result.absorbed_m2) and not np.any(result.returned_m2)


def test_a_craton_across_the_whole_front_passes_shortening_through():
    craton_cols = (3, 9)
    problem, _, col = _problem(8, 40)
    craton = (col >= craton_cols[0]) & (col < craton_cols[1])
    blocked, _, _ = _problem(8, 40, craton=craton)

    open_result = shortening.ring_cascade(problem)
    result = shortening.ring_cascade(blocked)

    _assert_conserved(blocked, result)
    assert result.absorbed_m2[craton].sum() < 1e-6 * blocked.demand_m2.sum()
    behind = col >= craton_cols[1]
    # The rigid indenter: the weak ground behind it takes up what the craton didn't.
    assert result.absorbed_m2[behind].sum() > 2.0 * open_result.absorbed_m2[behind].sum()


def test_shortening_escapes_around_a_craton_when_a_weak_route_exists():
    rows, cols = 15, 30
    problem, row, col = _problem(rows, cols)
    craton = (row >= 6) & (row <= 8) & (col >= 1) & (col <= 12)
    blocked, _, _ = _problem(rows, cols, craton=craton)

    open_result = shortening.ring_cascade(problem)
    result = shortening.ring_cascade(blocked)

    _assert_conserved(blocked, result)
    assert result.absorbed_m2[craton].sum() < 1e-6 * blocked.demand_m2.sum()
    flanks = ((row == 4) | (row == 5) | (row == 9) | (row == 10)) & (col >= 1) & (col <= 12)
    assert result.absorbed_m2[flanks].sum() > open_result.absorbed_m2[flanks].sum()
    # Going around costs little: most of what met the craton is still taken up.
    assert result.absorbed_m2.sum() > 0.95 * open_result.absorbed_m2.sum()


def test_a_capped_column_passes_its_shortening_on():
    problem, _, col = _problem(6, 30, demand_per_km=20_000.0)
    hc = np.full(len(col), lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hc[col <= 1] = lithosphere.MAX_CRUSTAL_THICKNESS_M
    capped, _, _ = _problem(6, 30, demand_per_km=20_000.0, hc=hc)

    open_result = shortening.ring_cascade(problem)
    result = shortening.ring_cascade(capped)

    _assert_conserved(capped, result)
    assert not np.any(result.absorbed_m2[col <= 1])
    assert result.absorbed_m2[col == 2].sum() > open_result.absorbed_m2[col == 2].sum()
    new_hc, new_hm = shortening.apply_strain(capped.hc_m, capped.hm_m, result.absorbed_m2 / capped.areas_m2)
    assert np.all(new_hc <= lithosphere.MAX_CRUSTAL_THICKNESS_M)
    assert np.all(new_hm <= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)


def test_a_column_nearing_its_cap_passes_shortening_on_before_it_fills():
    problem, _, col = _problem(6, 30, demand_per_km=50_000.0)
    hc = np.where(col <= 3, 0.95 * lithosphere.MAX_CRUSTAL_THICKNESS_M, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    full, _, _ = _problem(6, 30, demand_per_km=50_000.0, hc=hc)

    fresh = shortening.ring_cascade(problem)
    result = shortening.ring_cascade(full)

    _assert_conserved(full, result)
    room = shortening.room_m2(full.hc_m, full.hm_m, full.areas_m2)
    assert np.all(result.absorbed_m2 <= room * (1.0 + 1e-12))
    near_cap = col <= 3
    assert result.absorbed_m2[near_cap].sum() < 0.5 * fresh.absorbed_m2[near_cap].sum()
    assert result.absorbed_m2[~near_cap].sum() > fresh.absorbed_m2[~near_cap].sum()
    new_hc, _ = shortening.apply_strain(full.hc_m, full.hm_m, result.absorbed_m2 / full.areas_m2)
    assert new_hc.max() < lithosphere.MAX_CRUSTAL_THICKNESS_M


def test_cells_outside_the_host_take_their_own_demand_in_place():
    # A continental terrane (columns >= 3) on an oceanic plate: only it routes shortening on.
    problem, _, col = _problem(4, 20)
    host = col >= 3
    problem.demand_m2 = np.where((col == 0) | (col == 3), 1e9, 0.0)
    problem.front = (col == 0) | (col == 3)
    problem.host = host

    result = shortening.ring_cascade(problem)

    _assert_conserved(problem, result)
    np.testing.assert_allclose(result.absorbed_m2[col == 0], 1e9)
    assert not np.any(result.absorbed_m2[(col == 1) | (col == 2)])
    assert np.all(result.absorbed_m2[col == 4] > 0.0)


def test_relief_pushes_the_belt_outward():
    problem, _, col = _problem(6, 40)
    relief = np.where(col < 4, 6_000.0, 0.0)
    high, _, _ = _problem(6, 40, relief=relief)

    flat = shortening.ring_cascade(problem)
    result = shortening.ring_cascade(high)

    _assert_conserved(high, result)
    centroid = lambda r: float(np.dot(r.absorbed_m2, col) / r.absorbed_m2.sum())  # noqa: E731
    assert centroid(result) > centroid(flat) + 1.0
    assert result.absorbed_m2[col < 4].sum() < 0.5 * flat.absorbed_m2[col < 4].sum()


def test_unequal_cells_conserve_shortening():
    rng = np.random.default_rng(3)
    areas_km2 = 62.5**2 * rng.uniform(0.6, 1.4, size=(8, 30))
    problem, _, _ = _problem(8, 30, areas_km2=areas_km2)
    result = shortening.ring_cascade(problem)

    _assert_conserved(problem, result)
    # Each cell's strain is its own absorbed shortening over its own area.
    strain = result.absorbed_m2 / problem.areas_m2
    new_hc, _ = shortening.apply_strain(problem.hc_m, problem.hm_m, strain)
    volume_added = float(np.dot(new_hc - problem.hc_m, problem.areas_m2))
    assert volume_added == pytest.approx(float(np.dot(result.absorbed_m2, problem.hc_m)), rel=1e-9)


@pytest.mark.parametrize("craton", [False, True])
def test_node_density_does_not_change_the_result(craton):
    width_km, length_km = 1_000.0, 2_500.0
    profiles = []
    for cell_km in (62.5, 31.25):
        rows, cols = round(width_km / cell_km), round(length_km / cell_km)
        problem, _, col = _problem(rows, cols, cell_km=cell_km)
        if craton:
            x_km = (col + 0.5) * cell_km
            problem.craton_strength = ((x_km > 250.0) & (x_km < 625.0)).astype(float)
        result = shortening.ring_cascade(problem)
        _assert_conserved(problem, result)
        x_km = (col + 0.5) * cell_km
        bins = np.digitize(x_km, np.arange(0.0, length_km + 1.0, 250.0))
        profiles.append(np.bincount(bins, weights=result.absorbed_m2, minlength=12) / problem.demand_m2.sum())
    coarse, fine = profiles
    assert fine.sum() == pytest.approx(coarse.sum(), abs=0.02)
    np.testing.assert_allclose(fine, coarse, atol=0.03)


def test_demand_is_independent_of_the_band_s_width_in_cells():
    # The convergent band is a fixed number of cells wide, so a finer lattice has a thinner band.
    strain = np.full(3, 0.01)
    for spacing_km in (62.5, 31.25, 15.625):
        areas = np.full(3, (spacing_km * 1e3) ** 2)
        per_km_boundary = shortening.demand_m2(strain, areas, spacing_km / PLANET_RADIUS_KM).sum() / spacing_km
        assert per_km_boundary == pytest.approx(0.01 * 3 * shortening.REFERENCE_SPACING_KM * 1e6, rel=1e-12)


def test_convergent_strain_matches_the_band_local_deformation():
    rate = np.array([0.03, 0.05]) / rheology.SECONDS_PER_YEAR
    strain = rheology.convergent_strain(rate, 0.1, np.ones(2), 1.5)
    hc, hm = np.full(2, 30_000.0), np.full(2, 90_000.0)
    new_hc, new_hm, _ = rheology.apply_convergent_deformation(hc, hm, rate, 0.1, np.ones(2), strength=1.5)
    np.testing.assert_allclose(new_hc, hc * (1.0 + strain))
    np.testing.assert_allclose(new_hm, hm * (1.0 + strain))

import numpy as np
import pytest

from app import cratons, lithosphere, orogeny, synthetic_relief
from app.plates import PLANET_RADIUS_KM


def _arc(lon_deg: np.ndarray) -> np.ndarray:
    """Points along the equator at the given longitudes."""
    lon = np.radians(lon_deg)
    return np.stack([np.cos(lon), np.sin(lon), np.zeros_like(lon)], axis=1)


def _two_cratons(continent_of_core: dict[int, int]):
    # Craton 0 at -15..-5 deg, craton 1 at +5..+15 deg, open ground between them. Their
    # distances sum to ~1,110 km at the bisector, short of the far-pair fade.
    lon = np.linspace(-15.0, 15.0, 301)
    core = np.where(lon <= -5.0, 0, np.where(lon >= 5.0, 1, -1))
    return lon, _arc(lon), core, continent_of_core


@pytest.fixture
def one_middle_class(monkeypatch):
    # Every pair survives and draws the middle class, so the geometry is what's under test.
    monkeypatch.setattr(synthetic_relief, "SUTURE_DROP_PROBABILITY", 0.0)
    monkeypatch.setattr(
        synthetic_relief, "SUTURE_CLASSES", ((0.0, 4_500.0, 180.0, 3_500.0), (1.0, 2_500.0, 120.0, 0.0), (0.0, 900.0, 150.0, 0.0))
    )


@pytest.mark.usefixtures("one_middle_class")
def test_suture_belt_peaks_on_the_bisector_between_neighbouring_cratons():
    lon, points, core, continent_of_core = _two_cratons({0: 0, 1: 0})

    relief, young = synthetic_relief.suture_relief(points, core, continent_of_core, np.zeros(len(lon)), seed=5)

    assert lon[np.argmax(relief)] == pytest.approx(0.0, abs=0.2)
    assert relief.max() == pytest.approx(2_500.0, rel=1e-6)
    # A craton's edge is ~560 km from the bisector, far outside a 120 km sigma.
    assert relief[np.abs(lon) >= 5.0].max() < 1.0
    assert not young.any()  # The middle class is not young.


@pytest.mark.usefixtures("one_middle_class")
def test_suture_warp_moves_the_belt_off_the_bisector_without_splitting_it():
    lon, points, core, continent_of_core = _two_cratons({0: 0, 1: 0})

    for warp, direction in ((0.5, 1.0), (-0.5, -1.0)):
        relief, _ = synthetic_relief.suture_relief(points, core, continent_of_core, np.full(len(lon), warp), seed=5)

        # A positive warp moves the belt 150 km toward the higher-id craton, a negative one
        # away from it; either way one belt at full height, not two lowered ones.
        shift_km = lon[np.argmax(relief)] * np.pi / 180.0 * PLANET_RADIUS_KM
        assert shift_km == pytest.approx(direction * 0.5 * synthetic_relief.SUTURE_WARP_KM, abs=15.0)
        assert relief.max() == pytest.approx(2_500.0, rel=5e-3)  # peak between 11 km samples
        peaks = np.flatnonzero((relief[1:-1] > relief[:-2]) & (relief[1:-1] > relief[2:]))
        assert len(peaks) == 1


@pytest.mark.usefixtures("one_middle_class")
def test_cratons_on_different_continents_draw_no_suture():
    lon, points, core, continent_of_core = _two_cratons({0: 0, 1: 1})

    relief, young = synthetic_relief.suture_relief(points, core, continent_of_core, np.zeros(len(lon)), seed=5)

    assert not relief.any() and not young.any()


def test_suture_draws_are_keyed_by_seed_and_pair_not_node_order():
    lon, points, core, continent_of_core = _two_cratons({0: 0, 1: 0})
    order = np.random.default_rng(0).permutation(len(lon))

    relief, _ = synthetic_relief.suture_relief(points, core, continent_of_core, np.zeros(len(lon)), seed=11)
    shuffled, _ = synthetic_relief.suture_relief(points[order], core[order], continent_of_core, np.zeros(len(lon)), seed=11)

    assert np.allclose(shuffled, relief[order])


def test_land_mask_takes_the_highest_potential_continental_nodes_only():
    potential = np.array([5.0, 9.0, 1.0, 7.0, 8.0])
    continental = np.array([True, False, True, True, True])

    land = synthetic_relief.land_mask(potential, continental, 2)
    all_of_it = synthetic_relief.land_mask(potential, continental, 10)

    assert land.tolist() == [False, False, False, True, True]  # Node 1 is oceanic.
    assert all_of_it.tolist() == continental.tolist()


def test_back_columns_floats_each_column_at_its_target_within_the_caps():
    target = np.array([-2_000.0, -150.0, 5.0, 400.0, 4_000.0, np.nan])
    hc0 = np.full(6, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hm = np.full(6, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    rho = np.full(6, lithosphere.RHO_CONTINENTAL_CRUST)

    hc, lag = synthetic_relief.back_columns(target, hc0, hm, rho, np.zeros(6, dtype=bool), np.zeros(6))

    floated = lithosphere.isostatic_elevation(hc, hm, rho)
    assert np.allclose(floated[:5], target[:5], atol=1e-6)
    assert hc[5] == hc0[5]  # A NaN target (oceanic crust) keeps its column.
    assert np.all((hc >= lithosphere.MIN_CRUSTAL_THICKNESS_M) & (hc <= lithosphere.MAX_CRUSTAL_THICKNESS_M))
    assert not lag.any()  # Nothing young: every column starts at steady state.


def test_back_columns_keeps_cratonic_cores_at_host_thickness():
    target = np.array([-2_000.0, -2_000.0])
    hc0 = np.full(2, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hm = np.full(2, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    rho = np.full(2, lithosphere.RHO_CONTINENTAL_CRUST)

    hc, _ = synthetic_relief.back_columns(target, hc0, hm, rho, np.array([True, False]), np.zeros(2))

    assert hc[0] == pytest.approx(cratons.CRATON_HOST_MIN_HC_M)
    assert hc[1] < cratons.CRATON_HOST_MIN_HC_M


def test_young_belts_start_with_the_lag_of_freshly_buried_crust():
    target = np.array([4_000.0, 4_000.0, 4_000.0])
    hc0 = np.full(3, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hm = np.full(3, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    rho = np.full(3, lithosphere.RHO_CONTINENTAL_CRUST)

    hc, lag = synthetic_relief.back_columns(target, hc0, hm, rho, np.zeros(3, dtype=bool), np.array([1.0, 0.5, 0.0]))

    steady_ref = orogeny.steady_moho_temperature_c(lithosphere.REFERENCE_HC_CONTINENTAL_M, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    # Fully young: the thickened Moho keeps the reference column's temperature.
    assert orogeny.moho_temperature_c(hc[0], hm[0], lag[0]) == pytest.approx(steady_ref, abs=1e-6)
    assert lag[1] == pytest.approx(0.5 * lag[0])
    assert lag[2] == 0.0

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


@pytest.fixture
def young_only(monkeypatch):
    # Every pair survives as a young belt with a plateau, the case that can reach a craton.
    monkeypatch.setattr(synthetic_relief, "SUTURE_DROP_PROBABILITY", 0.0)
    monkeypatch.setattr(
        synthetic_relief, "SUTURE_CLASSES", ((1.0, 4_500.0, 180.0, 3_500.0), (0.0, 2_500.0, 120.0, 0.0), (0.0, 900.0, 150.0, 0.0))
    )


def test_a_young_plateau_stops_before_the_craton_it_faces(monkeypatch):
    # Cratons ~450 km apart, as the layout's 400 km belt width allows.
    lon = np.linspace(-15.0, 15.0, 301)
    core = np.where(lon <= -2.0, 0, np.where(lon >= 2.0, 1, -1))
    points = _arc(lon)
    monkeypatch.setattr(synthetic_relief, "SUTURE_DROP_PROBABILITY", 0.0)

    def relief_with(plateau_m: float, points, core) -> np.ndarray:
        classes = ((1.0, 4_500.0, 180.0, plateau_m), (0.0, 2_500.0, 120.0, 0.0), (0.0, 900.0, 150.0, 0.0))
        monkeypatch.setattr(synthetic_relief, "SUTURE_CLASSES", classes)
        return synthetic_relief.suture_relief(points, core, {0: 0, 1: 0}, np.zeros(len(points)), seed=5)[0]

    # #218 review: the plateau used to run 450 km to one side, straight onto a shield that
    # then held the thick crust as craton. It now fades out before the craton it faces, so
    # on craton nodes only the belt's own flank remains.
    with_plateau, belt_only = relief_with(3_500.0, points, core), relief_with(0.0, points, core)
    assert np.allclose(with_plateau[core >= 0], belt_only[core >= 0], atol=1.0)

    # Where the gap is wide enough, the plateau still stands beside the belt.
    _, wide_points, wide_core, _ = _two_cratons({0: 0, 1: 0})
    wide, wide_belt = relief_with(3_500.0, wide_points, wide_core), relief_with(0.0, wide_points, wide_core)
    assert np.any(wide[wide_core < 0] > wide_belt[wide_core < 0] + 500.0)


def _grid(half_deg: float = 12.0, step_deg: float = 0.18):
    lat, lon = np.meshgrid(np.arange(-half_deg, half_deg, step_deg), np.arange(-half_deg, half_deg, step_deg), indexing="ij")
    la, lo = np.radians(lat), np.radians(lon)
    points = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)
    return lat, lon, points


@pytest.mark.usefixtures("young_only")
def test_three_cratons_meet_without_cliffs():
    lat, lon, points = _grid()
    centres = ((-6.0, -6.0), (-6.0, 6.0), (6.0, 0.0))
    core = np.full(lat.shape, -1)
    for i, (clat, clon) in enumerate(centres):
        core[np.hypot(lat - clat, lon - clon) < 2.5] = i

    for seed in range(4):
        warp = np.zeros(lat.size)
        relief, _ = synthetic_relief.suture_relief(points.reshape(-1, 3), core.ravel(), {0: 0, 1: 0, 2: 0}, warp, seed)
        relief = relief.reshape(lat.shape)
        # #218 review: a node took only its two nearest cratons' belt, so belts stopped dead
        # (up to 3.6 km between 20 km neighbours) where the second-nearest craton changed.
        # The steepest designed edge, a plateau's 150 km smoothstep, is ~35 m/km.
        jump = max(np.abs(np.diff(relief, axis=0)).max(), np.abs(np.diff(relief, axis=1)).max())
        assert jump < 1_000.0


def test_leading_edge_scores_convergence_relative_to_the_ocean_plate():
    point = _arc(np.array([0.0]))
    ocean = _arc(np.array([2.0]))  # to the east
    eastward = np.array([[0.0, 0.0, 1e-8]])  # rad/yr about +z: ~6.4 cm/yr east at the equator
    still = np.zeros((1, 3))

    toward, passive = synthetic_relief.leading_edge(point, eastward, np.array([0]), ocean, still, np.array([1]))
    away, _ = synthetic_relief.leading_edge(point, still, np.array([0]), ocean, eastward, np.array([1]))
    own_ocean, own_passive = synthetic_relief.leading_edge(point, eastward, np.array([0]), ocean, eastward, np.array([0]))

    assert toward[0] == pytest.approx(1.0, abs=1e-6) and not passive[0]
    assert away[0] == pytest.approx(-1.0, abs=1e-6)
    # Ocean on the continent's own plate: nothing converges, whatever the plate's motion.
    assert own_ocean[0] == 0.0 and own_passive[0]


def test_back_columns_caps_cratonic_cores_at_the_host_maximum():
    target = np.array([4_000.0, 4_000.0])
    hc0 = np.full(2, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hm = np.full(2, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    rho = np.full(2, lithosphere.RHO_CONTINENTAL_CRUST)

    hc, _ = synthetic_relief.back_columns(target, hc0, hm, rho, np.array([True, False]), np.zeros(2))

    assert hc[0] == pytest.approx(cratons.CRATON_HOST_MAX_HC_M)
    assert hc[1] > cratons.CRATON_HOST_MAX_HC_M

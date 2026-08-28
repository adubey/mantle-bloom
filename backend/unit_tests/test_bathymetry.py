import numpy as np
from scipy.spatial import cKDTree

from app import bathymetry, lithosphere
from app.elevation_lines import line_spacing_rad
from app.lithosphere_plate import generate_plates

_GEN_KWARGS = dict(seed=7, num_plates=12, continental_fraction=0.5, land_fraction=0.2, node_density=1.0)


def _node_arrays(plates):
    """(points, elevation, is_continental, dist_to_land_km) concatenated across every plate."""
    points = np.concatenate([p.all_points_and_elevation()[0] for p in plates])
    elevation = np.concatenate([np.asarray(p.collect("elevation"), dtype=float) for p in plates])
    is_continental = np.concatenate(
        [np.full(p.node_count(), p.crust_type == "continental", dtype=bool) for p in plates]
    )
    land = elevation > 0.0
    dist_chord, _ = cKDTree(points[land]).query(points)
    dist_km = 2.0 * np.arcsin(np.clip(dist_chord / 2.0, 0.0, 1.0)) * 6371.0
    return points, elevation, is_continental, dist_km


def _generate(monkeypatch, subside=True, margins=True):
    if not subside:
        monkeypatch.setattr(bathymetry, "_subside_offshore_continental_crust", lambda plates: None)
    if not margins:
        monkeypatch.setattr(bathymetry, "_smooth_continental_margins", lambda plates: None)
    plates = generate_plates(**_GEN_KWARGS)
    monkeypatch.undo()
    return plates


def test_offshore_continental_crust_subsides_with_distance_from_land(monkeypatch):
    """Submerged continental crust should be drawn down toward abyssal depth as it gets
    farther from land -- flat and shelf-shallow at every distance without the pass."""
    _, raw_elev, raw_isc, raw_dist = _node_arrays(_generate(monkeypatch, subside=False))
    _, sm_elev, sm_isc, sm_dist = _node_arrays(_generate(monkeypatch))

    def median_depth(elev, isc, dist, lo, hi):
        m = isc & (elev < 0.0) & (dist >= lo) & (dist < hi)
        return np.median(elev[m]) if np.any(m) else None

    # Without the pass, depth barely varies with distance from land.
    raw_near = median_depth(raw_elev, raw_isc, raw_dist, 0.0, 200.0)
    raw_far = median_depth(raw_elev, raw_isc, raw_dist, bathymetry.OFFSHORE_ABYSSAL_KM, 1e9)
    assert raw_near is not None and raw_far is not None
    assert abs(raw_far - raw_near) < 2500.0

    # With it, near-shore holds while far-offshore reaches basin depth.
    sm_near = median_depth(sm_elev, sm_isc, sm_dist, 0.0, 200.0)
    sm_mid = median_depth(sm_elev, sm_isc, sm_dist, 400.0, 1000.0)
    sm_far = median_depth(sm_elev, sm_isc, sm_dist, bathymetry.OFFSHORE_ABYSSAL_KM, 1e9)
    assert sm_near > -800.0
    assert sm_far < bathymetry.ABYSSAL_REFERENCE_DEPTH_M + 800.0
    assert sm_near > sm_mid > sm_far  # monotonic ramp


def test_near_shore_shelf_is_left_shallow(monkeypatch):
    """Submerged continental crust within the shelf range keeps its shallow depth -- the
    subsidence is a distance-from-land ramp, not a blanket deepening."""
    _, elev, isc, dist = _node_arrays(_generate(monkeypatch))
    shelf = isc & (elev < 0.0) & (dist < bathymetry.OFFSHORE_SHELF_KM)
    assert np.any(shelf)
    assert np.median(elev[shelf]) > -1500.0


def test_subsidence_never_lifts_crust(monkeypatch):
    """The pass only drowns crust -- no submerged continental node ends up shallower than the
    generation noise left it (margins pass disabled in both, to isolate the subsidence)."""
    _, raw_elev, raw_isc, _ = _node_arrays(_generate(monkeypatch, subside=False, margins=False))
    _, sm_elev, _, _ = _node_arrays(_generate(monkeypatch, subside=True, margins=False))
    moved = raw_isc & (raw_elev < 0.0)
    assert np.all(sm_elev[moved] <= raw_elev[moved] + 1e-6)


def test_margin_smoothing_softens_the_continent_ocean_cliff(monkeypatch):
    """Grading the plate-boundary columns should cut the elevation step across a
    continent/ocean contact substantially."""
    spacing_rad = line_spacing_rad(_GEN_KWARGS["node_density"])
    neighbour_chord = 2.0 * np.sin(1.5 * spacing_rad / 2.0)

    def cliff(plates):
        points, elevation, is_continental, _ = _node_arrays(plates)
        cont = is_continental & (elevation < 0.0)
        ocean = ~is_continental
        dist, idx = cKDTree(points[ocean]).query(points[cont])
        close = dist < neighbour_chord
        return np.abs(elevation[cont][close] - elevation[ocean][idx][close])

    raw = cliff(_generate(monkeypatch, subside=False, margins=False))
    smoothed = cliff(_generate(monkeypatch, subside=False, margins=True))
    assert len(raw) > 20 and len(smoothed) > 20
    assert np.median(smoothed) < 0.65 * np.median(raw)


def test_deep_ocean_interior_stays_deep():
    """Oceanic crust far from any continent keeps full abyssal depth -- neither pass touches
    it."""
    plates = generate_plates(**_GEN_KWARGS)
    points, elevation, is_continental, _ = _node_arrays(plates)
    dist_to_continent, _ = cKDTree(points[is_continental]).query(points)
    far = ~is_continental & (dist_to_continent > 2.0 * np.sin(bathymetry.MARGIN_TRANSITION_RAD / 2.0))
    assert elevation[far].min() < -4000.0


def test_bathymetry_shaping_keeps_elevation_consistent_with_isostasy():
    """Both passes mutate Hc/Hm and must re-sync `elevation`, not leave it stale."""
    for plate in generate_plates(**_GEN_KWARGS):
        rho_c = lithosphere.crust_density(plate.crust_type)
        for line in plate.lines:
            if len(line) == 0:
                continue
            expected = lithosphere.isostatic_elevation(
                line.crustal_thickness_m, line.mantle_lithosphere_thickness_m, rho_c
            )
            assert np.allclose(line.elevation, expected, atol=1e-6)


def test_shape_initial_bathymetry_no_op_on_single_crust_world():
    """A world with no continents has neither margins nor drowned continental crust -- must
    not raise or change anything."""
    plates = generate_plates(seed=3, num_plates=6, continental_fraction=0.0, node_density=0.5)
    before = np.concatenate([np.asarray(p.collect("elevation"), dtype=float) for p in plates])
    bathymetry.shape_initial_bathymetry(plates)
    after = np.concatenate([np.asarray(p.collect("elevation"), dtype=float) for p in plates])
    assert np.array_equal(before, after)


def test_crustal_thickness_inverse_round_trips():
    """`crustal_thickness_for_submerged_elevation` is the exact inverse of
    `isostatic_elevation`'s water-loaded branch."""
    hm = np.full(6, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    target = np.array([-500.0, -1500.0, -3000.0, -4500.0, -5222.0, -6000.0])
    hc = lithosphere.crustal_thickness_for_submerged_elevation(target, hm, lithosphere.RHO_CONTINENTAL_CRUST)
    back = lithosphere.isostatic_elevation(hc, hm, lithosphere.RHO_CONTINENTAL_CRUST)
    assert np.allclose(back, target, atol=1e-6)

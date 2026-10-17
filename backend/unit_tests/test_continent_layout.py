import numpy as np
import pytest

from app import continent_layout as cl
from app import plates as plates_mod
from app.continent_layout import ContinentLayout
from app.elevation_lines import PLANET_RADIUS_KM
from app.layout_lattice import cap_area_km2
from app.lithosphere_plate import cratons_first_layout
from app.noise import SphereNoise


@pytest.fixture(scope="module")
def layout() -> ContinentLayout:
    return cratons_first_layout(7, land_fraction=0.29).layout


def test_belts_separate_the_cores_with_no_unmarked_gap(layout):
    lat = layout.lattice
    rows, cols = lat.adj_rows, lat.adj_neighbours
    a, b = layout.core[rows], layout.core[cols]
    # No link joins two different cores.
    assert not np.any((a >= 0) & (b >= 0) & (a != b))
    # Along the arc between neighbouring nuclei of two cratons of one continent, everything
    # is craton or belt.
    spacing_rad = cl.CRATON_CORE_SPACING_KM / PLANET_RADIUS_KM
    checked = 0
    for i in range(len(layout.nuclei)):
        for j in range(i + 1, len(layout.nuclei)):
            ci, cj = layout.nucleus_craton[i], layout.nucleus_craton[j]
            if ci == cj or layout.nucleus_continent[i] != layout.nucleus_continent[j]:
                continue
            if np.arccos(np.clip(layout.nuclei[i] @ layout.nuclei[j], -1, 1)) > 1.3 * spacing_rad:
                continue
            if not (np.any(layout.core == ci) and np.any(layout.core == cj)):
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
    again = cratons_first_layout(7, land_fraction=0.29).layout
    for name in ("continent", "core", "belt", "plate"):
        assert np.array_equal(getattr(again, name), getattr(layout, name))
    assert not np.array_equal(cratons_first_layout(8, land_fraction=0.29).layout.plate, layout.plate)


def test_a_continent_too_small_for_a_core_gets_none():
    small = cratons_first_layout(4, continental_fraction=0.05).layout
    disk = cap_area_km2(cl.CRATON_CORE_RADIUS_KM)
    too_small = np.flatnonzero(small.continent_area_km2 < disk)
    assert len(too_small) > 0
    for k in too_small:
        assert not np.any(small.is_core & (small.continent == k))
        assert np.any(small.continent == k)



def test_some_nuclei_join_an_existing_craton(layout):
    """Cratons are composites (#470): some cratons have several overlapping nuclei, and no
    craton id is shared across continents."""
    counts = np.bincount(layout.nucleus_craton)
    assert counts.max() > 1
    assert np.sum(counts == 1) > 0
    for craton in np.unique(layout.nucleus_craton):
        assert len(np.unique(layout.nucleus_continent[layout.nucleus_craton == craton])) == 1
    # A joining nucleus overlaps its craton, so its nucleus point lies inside it.
    joined = np.flatnonzero(counts[layout.nucleus_craton] > 1)
    inside = layout.core[layout.nearest(layout.nuclei[joined])] == layout.nucleus_craton[joined]
    assert inside.mean() > 0.8


def test_without_clustering_every_nucleus_is_its_own_craton(monkeypatch):
    monkeypatch.setattr(cl, "CRATON_CLUSTER_CONCENTRATION", 1e12)
    layout = cratons_first_layout(7, land_fraction=0.29).layout
    assert len(np.unique(layout.nucleus_craton)) == len(layout.nucleus_craton)


def test_some_cratons_reach_a_continental_edge(layout, monkeypatch):
    """Variable growth and the centre pull leave some cratons near their continent's edge
    (#470); with neither, the even margin keeps every craton inside."""
    assert layout.craton_diagnostics()["cratons_near_edge_share"] > 0.1
    monkeypatch.setattr(cl, "CONTINENT_CENTRE_PULL", 0.0)
    monkeypatch.setattr(cl, "CONTINENT_GROWTH_SPEED_SIGMA", 0.0)
    even = cratons_first_layout(7, land_fraction=0.29).layout
    assert even.craton_diagnostics()["cratons_near_edge_share"] < layout.craton_diagnostics()["cratons_near_edge_share"]


def test_craton_sizes_spread_below_one_disk_and_stop_near_the_cap(layout):
    """Size scales let a craton no nucleus joins stay well under one core disk, and joining
    stops at CRATON_MAX_AREA_KM2 (#470 review); the final trim can leave a craton a little
    over it."""
    area = layout.lattice.area_km2
    sizes = np.array([area[layout.core == c].sum() for c in np.unique(layout.core[layout.is_core])])
    disk = cap_area_km2(cl.CRATON_CORE_RADIUS_KM)
    assert sizes.min() < 0.5 * disk
    assert sizes.max() < 1.5 * cl.CRATON_MAX_AREA_KM2

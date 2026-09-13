import numpy as np
from app import mantle, real_plates, worldsketch


def test_load_major_plates_returns_sane_unit_geometry():
    plates = real_plates.load_major_plates()
    assert len(plates) == 16
    names = {p.name for p in plates}
    assert {"Pacific", "North America", "Eurasia", "Africa", "Antarctica"} <= names
    for p in plates:
        assert len(p.boundary_xyz) >= 3
        norms = np.linalg.norm(p.boundary_xyz, axis=-1)
        assert np.allclose(norms, 1.0, atol=1e-6)
        assert p.omega.shape == (3,)
        # A real absolute plate rate is on the order of a few cm/yr, not zero and not wildly
        # larger than mantle.py's own plausible-rate bounds.
        rate = np.linalg.norm(p.omega)
        assert 0 < rate < mantle.MAX_PLATE_RATE


def test_load_ridge_trench_candidates_have_both_polarities():
    candidates = real_plates.load_ridge_trench_candidates()
    assert len(candidates) > 20
    signs = {c.sign for c in candidates}
    assert signs == {1.0, -1.0}
    for c in candidates:
        assert np.isclose(np.linalg.norm(c.position), 1.0, atol=1e-6)


def _earth_sketch_masks():
    # A synthetic half-land/half-sea checkerboard is enough to exercise real_plate_sites'
    # land/sea split without depending on the frontend's baked-in Earth PNG -- every real
    # plate's own polygon still needs the actual boundary geometry, but `sketch.land`/
    # `sketch.cell_weight` only need to exist with the right shape and a mix of both.
    h, w = worldsketch.SKETCH_GRID_H, worldsketch.SKETCH_GRID_W
    rows, cols = np.indices((h, w))
    land = ((rows // 20 + cols // 20) % 2 == 0)
    mountain = np.zeros_like(land)
    river = np.zeros_like(land)
    return worldsketch.SketchMasks(land=land, mountain=mountain, river=river)


def test_real_plate_sites_distributes_across_land_and_sea():
    masks = _earth_sketch_masks()
    plates = real_plates.load_major_plates()
    rng = np.random.default_rng(0)

    site_xyz, crust_types = real_plates.real_plate_sites(masks, plates, num_plates=32, num_continents=16, rng=rng)

    assert len(site_xyz) == len(crust_types) == 32
    assert set(crust_types) <= {"continental", "oceanic"}
    assert crust_types.count("continental") + crust_types.count("oceanic") == 32
    norms = np.linalg.norm(site_xyz, axis=-1)
    assert np.allclose(norms, 1.0, atol=1e-6)


def test_fit_mantle_centers_beats_a_naive_sign_only_baseline():
    plates = real_plates.load_major_plates()
    candidates = real_plates.load_ridge_trench_candidates()
    rng = np.random.default_rng(0)

    fitted = real_plates.fit_mantle_centers(plates, candidates, rng)
    assert len(fitted) == len(candidates)

    points, targets = real_plates._target_samples(plates, np.random.default_rng(0))
    fitted_flow = mantle.flow_at(points, fitted)
    fitted_err = np.linalg.norm(fitted_flow - targets, axis=-1).mean()

    naive = [
        mantle.ConvectionCenter(position=c.position, strength=c.sign * real_plates._INITIAL_STRENGTH, falloff=real_plates._CANDIDATE_FALLOFF_RAD)
        for c in candidates
    ]
    naive_flow = mantle.flow_at(points, naive)
    naive_err = np.linalg.norm(naive_flow - targets, axis=-1).mean()

    assert fitted_err < naive_err

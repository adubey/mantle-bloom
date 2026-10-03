import numpy as np

from app import lithosphere


def test_isostatic_elevation_reference_columns_land_near_v1_baselines():
    """Reference Hc/Hm for each crust type should land close to v1's own
    BASE_CONTINENTAL_M(200)/BASE_OCEANIC_M(-3800) baselines, and continents should sit well
    above oceans -- the calibration `ISOSTATIC_REFERENCE_OFFSET_M` exists for (see
    lithosphere.py's own docstring)."""
    hc_c, hm_c = lithosphere.reference_thickness("continental")
    z_c = lithosphere.isostatic_elevation(np.array([hc_c]), np.array([hm_c]), lithosphere.RHO_CONTINENTAL_CRUST)[0]
    assert 0.0 <= z_c <= 1000.0

    hc_o, hm_o = lithosphere.REFERENCE_HC_OCEANIC_M, lithosphere.REFERENCE_HM_OCEANIC_M
    z_o = lithosphere.isostatic_elevation(np.array([hc_o]), np.array([hm_o]), lithosphere.RHO_OCEANIC_CRUST)[0]
    assert -7000.0 <= z_o <= -2000.0
    assert z_c - z_o > 3000.0  # a real continent/ocean contrast, not a flat planet


def test_isostatic_elevation_monotonic_in_crustal_thickness():
    """Thicker crust should always float higher (Airy isostasy's core physical claim) --
    holding Hm fixed, z must be strictly increasing in Hc."""
    hm = np.full(20, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    hc = np.linspace(5_000.0, 60_000.0, 20)
    z = lithosphere.isostatic_elevation(hc, hm, lithosphere.RHO_CONTINENTAL_CRUST)
    assert np.all(np.diff(z) > 0)


def test_isostatic_elevation_thicker_mantle_lithosphere_sinks_the_column():
    """Hm is *denser* than the asthenosphere it displaces (rho_m > rho_a) -- more of it
    should pull the column down, not lift it."""
    hc = np.full(20, lithosphere.REFERENCE_HC_OCEANIC_M)
    hm = np.linspace(10_000.0, 150_000.0, 20)
    z = lithosphere.isostatic_elevation(hc, hm, lithosphere.RHO_OCEANIC_CRUST)
    assert np.all(np.diff(z) < 0)


def test_isostatic_elevation_broadcasts_a_per_node_crust_density():
    """erosion.apply_erosion spans several plates at once, so it passes rho_c as a per-node
    array -- the result must match calling the scalar form plate by plate, and a small Hc
    change must move the surface by well under that change (Airy rebound eats most of it)."""
    hc = np.array([40_000.0, 40_000.0, 18_000.0, 18_000.0])
    hm = np.full(4, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    rho_c = np.array(
        [lithosphere.RHO_CONTINENTAL_CRUST, lithosphere.RHO_OCEANIC_CRUST] * 2
    )

    vectorized = lithosphere.isostatic_elevation(hc, hm, rho_c)
    per_node = np.array(
        [lithosphere.isostatic_elevation(hc[i : i + 1], hm[i : i + 1], rho_c[i])[0] for i in range(4)]
    )
    assert np.allclose(vectorized, per_node)

    # Remove 100 m of crustal column; the subaerial columns (indices 0-1) are still well
    # above sea level and drop far less than 100 m.
    dropped = lithosphere.isostatic_elevation(hc - 100.0, hm, rho_c)
    surface_drop = vectorized[:2] - dropped[:2]
    assert np.all(surface_drop > 0.0) and np.all(surface_drop < 25.0)


def test_crustal_thickness_for_elevation_round_trips_both_branches():
    """`crustal_thickness_for_elevation` (issue #173) is the general inverse of
    `isostatic_elevation`, covering both the dry (above-datum) and water-loaded branches --
    unlike `crustal_thickness_for_submerged_elevation`, which only covers `z < 0`."""
    hm = np.full(6, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    target = np.array([-6000.0, -3000.0, -500.0, 200.0, 1500.0, 4000.0])
    hc = lithosphere.crustal_thickness_for_elevation(target, hm, lithosphere.RHO_CONTINENTAL_CRUST)
    back = lithosphere.isostatic_elevation(hc, hm, lithosphere.RHO_CONTINENTAL_CRUST)
    assert np.allclose(back, target, atol=1e-6)


def test_isostatic_elevation_clips_to_world_bounds():
    from app.elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M

    z = lithosphere.isostatic_elevation(np.array([1e9, -1e9]), np.array([0.0, 0.0]), lithosphere.RHO_CONTINENTAL_CRUST)
    assert z[0] == MAX_ELEVATION_M
    assert z[1] >= MIN_ELEVATION_M


def test_moment_of_inertia_tensor_symmetric_and_positive_definite():
    rng = np.random.default_rng(0)
    points = rng.normal(size=(200, 3))
    points /= np.linalg.norm(points, axis=-1, keepdims=True)
    hc = np.full(200, 35_000.0)
    hm = np.full(200, 100_000.0)
    inertia = lithosphere.moment_of_inertia_tensor(points, hc, hm, lithosphere.RHO_CONTINENTAL_CRUST, spacing_rad=0.02)
    assert np.allclose(inertia, inertia.T)
    eigenvalues = np.linalg.eigvalsh(inertia)
    assert np.all(eigenvalues > 0)


def test_omega_from_angular_momentum_inverts_angular_momentum():
    rng = np.random.default_rng(1)
    points = rng.normal(size=(100, 3))
    points /= np.linalg.norm(points, axis=-1, keepdims=True)
    hc = np.full(100, 35_000.0)
    hm = np.full(100, 100_000.0)
    inertia = lithosphere.moment_of_inertia_tensor(points, hc, hm, lithosphere.RHO_CONTINENTAL_CRUST, spacing_rad=0.02)
    omega = np.array([1e-9, -2e-9, 3e-9])
    l = lithosphere.angular_momentum(inertia, omega)
    recovered = lithosphere.omega_from_angular_momentum(inertia, l)
    assert np.allclose(recovered, omega, rtol=1e-6)


def test_clamp_column_caps_clips_both_bounds_and_keeps_the_elevation_residual():
    from app.elevation_lines import line_spacing_rad
    from app.lithosphere_plate import new_plate

    plate = new_plate(0, np.eye(3), "continental", line_spacing_rad(1.0), seed=1, is_owned=lambda pts: pts[:, 2] > 0.95)
    n = plate.node_count()
    assert n >= 4
    hc = np.full(n, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hm = np.full(n, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    hc[0], hm[1] = 1.5 * lithosphere.MAX_CRUSTAL_THICKNESS_M, 0.5 * lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
    hm[2] = 1.01 * lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M
    residual = np.linspace(-50.0, 50.0, n)
    rho_c = lithosphere.crust_density("continental")
    plate.set_fields_on_plate(
        crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm, elevation=lithosphere.isostatic_elevation(hc, hm, rho_c) + residual
    )

    assert lithosphere.clamp_column_caps(plate)

    new_hc = plate.collect("crustal_thickness_m")
    new_hm = plate.collect("mantle_lithosphere_thickness_m")
    assert new_hc[0] == lithosphere.MAX_CRUSTAL_THICKNESS_M
    assert new_hm[1] == lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
    assert new_hm[2] == lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M
    np.testing.assert_array_equal(new_hc[1:], hc[1:])
    np.testing.assert_allclose(plate.collect("elevation"), lithosphere.isostatic_elevation(new_hc, new_hm, rho_c) + residual)
    assert not lithosphere.clamp_column_caps(plate)


# --- Issue #275 phase 3: reversible ice loading ----------------------------------------------


def test_grounded_ice_load_counts_only_ice_the_water_column_cannot_float():
    sea = 0.0
    depth = np.array([1000.0, 1000.0, 1000.0, 1000.0, 0.0])
    bed = np.array([500.0, -500.0, -917.0, -2000.0, 500.0])
    load = lithosphere.grounded_ice_load_kg_m2(depth, bed, sea)
    assert load[0] == lithosphere.RHO_ICE * 1000.0  # on land: the whole ice column
    assert np.isclose(load[1], lithosphere.RHO_ICE * 1000.0 - lithosphere.RHO_WATER * 500.0)
    assert np.isclose(load[2], 0.0)  # exactly at flotation
    assert load[3] == 0.0 and load[4] == 0.0  # floating shelf / sea ice, and no ice


def test_ice_load_deflection_sinks_dry_land_by_rho_ice_over_rho_a():
    hc, hm = lithosphere.reference_thickness("continental")
    hc = np.full(3, hc + 5_000.0)  # comfortably above sea level, so the dry branch holds
    hm = np.full(3, hm)
    ice = np.array([0.0, 1000.0, 2000.0])
    load = lithosphere.RHO_ICE * ice
    w = lithosphere.ice_load_deflection(hc, hm, lithosphere.RHO_CONTINENTAL_CRUST, load)
    assert w[0] == 0.0
    assert np.allclose(w[1:], -ice[1:] * lithosphere.RHO_ICE / lithosphere.RHO_ASTHENOSPHERE)
    # No Hc tracking (v1 lines): the same dry-land response.
    v1 = lithosphere.ice_load_deflection(np.zeros(3), np.zeros(3), lithosphere.RHO_CONTINENTAL_CRUST, load)
    assert np.allclose(v1, w)


def test_ice_load_deflection_is_deeper_once_the_depression_floods():
    hm = np.full(1, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    hc = lithosphere.crustal_thickness_for_elevation(np.array([-1000.0]), hm, lithosphere.RHO_CONTINENTAL_CRUST)
    load = np.array([lithosphere.RHO_ICE * 3000.0])
    w = lithosphere.ice_load_deflection(hc, hm, lithosphere.RHO_CONTINENTAL_CRUST, load)
    dry = -load / lithosphere.RHO_ASTHENOSPHERE
    assert np.isclose(w[0], dry[0] * lithosphere.RHO_ASTHENOSPHERE / (lithosphere.RHO_ASTHENOSPHERE - lithosphere.RHO_WATER))


def test_ice_load_pa_is_weight_per_area():
    assert np.isclose(lithosphere.ice_load_pa(np.array([1000.0]))[0], 1000.0 * lithosphere.GRAVITY_M_S2)


def test_sync_line_elevation_keeps_the_ice_load_deflection():
    from app.elevation_lines import ElevationLine

    hc, hm = lithosphere.reference_thickness("continental")
    line = ElevationLine(
        phi=0.0,
        theta=np.array([0.1, 0.2]),
        elevation=np.zeros(2),
        crustal_thickness_m=np.full(2, hc),
        mantle_lithosphere_thickness_m=np.full(2, hm),
        ice_load_deflection_m=np.array([0.0, -300.0]),
    )
    synced = lithosphere.sync_line_elevation(line, lithosphere.RHO_CONTINENTAL_CRUST)
    bare = lithosphere.isostatic_elevation(line.crustal_thickness_m, line.mantle_lithosphere_thickness_m, lithosphere.RHO_CONTINENTAL_CRUST)
    assert np.allclose(synced.elevation, bare + np.array([0.0, -300.0]))


def test_sync_line_elevation_stores_only_the_deflection_the_floor_lets_through():
    from app.elevation_lines import MIN_ELEVATION_M, ElevationLine

    hm = np.full(1, lithosphere.REFERENCE_HM_CONTINENTAL_M)
    hc = lithosphere.crustal_thickness_for_elevation(np.array([MIN_ELEVATION_M + 50.0]), hm, lithosphere.RHO_CONTINENTAL_CRUST)
    line = ElevationLine(
        phi=0.0,
        theta=np.array([0.1]),
        elevation=np.zeros(1),
        crustal_thickness_m=hc,
        mantle_lithosphere_thickness_m=hm,
        ice_load_deflection_m=np.array([-300.0]),
    )
    synced = lithosphere.sync_line_elevation(line, lithosphere.RHO_CONTINENTAL_CRUST)
    assert synced.elevation[0] == MIN_ELEVATION_M
    assert np.isclose(synced.ice_load_deflection_m[0], -50.0)

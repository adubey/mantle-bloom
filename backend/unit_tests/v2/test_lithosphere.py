import numpy as np

from app.v2 import lithosphere


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

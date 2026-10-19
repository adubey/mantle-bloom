import numpy as np

from app.quasistatic import column_gpe_anomaly_n_per_m, solve_quasistatic_contacts


def test_column_gpe_anomaly_is_zero_for_the_reference_column():
    result = column_gpe_anomaly_n_per_m(
        np.array([20_000.0, 30_000.0]),
        np.array([3_250.0, 3_250.0]),
        300_000.0,
        asthenosphere_density_kg_per_m3=3_250.0,
    )

    assert result == 0.0


def test_column_gpe_anomaly_uses_one_shared_compensation_plane():
    upper = column_gpe_anomaly_n_per_m(
        np.array([10_000.0]),
        np.array([3_300.0]),
        300_000.0,
        asthenosphere_density_kg_per_m3=3_250.0,
    )
    deeper = column_gpe_anomaly_n_per_m(
        np.array([10_000.0, 10_000.0]),
        np.array([3_250.0, 3_300.0]),
        300_000.0,
        asthenosphere_density_kg_per_m3=3_250.0,
    )

    assert upper > deeper > 0.0


def test_zero_load_and_no_resistance_selects_zero_rate():
    result = solve_quasistatic_contacts(np.zeros((3, 3)), np.zeros(3))

    np.testing.assert_array_equal(result.omega, np.zeros(3))
    np.testing.assert_array_equal(result.residual_torque_nm, np.zeros(3))


def test_unbalanced_torque_is_reported_instead_of_becoming_motion():
    torque = np.array([0.0, 0.0, 12.0])
    result = solve_quasistatic_contacts(np.zeros((3, 3)), torque)

    np.testing.assert_array_equal(result.omega, np.zeros(3))
    np.testing.assert_array_equal(result.residual_torque_nm, -torque)


def test_contact_force_is_equal_opposite_and_changes_both_plate_rates():
    drag = np.eye(2)
    torque = np.array([1.0, -1.0])
    contact_torque = np.array([[-1.0], [1.0]])
    closing_map = np.array([[1.0, -1.0]])

    result = solve_quasistatic_contacts(drag, torque, contact_torque, closing_map, np.array([5.0]))

    assert result.converged
    assert result.contact_force_n_per_m[0] == 1.0
    np.testing.assert_allclose(result.omega, [0.0, 0.0], atol=1e-12)
    assert result.closing_rate[0] <= 1e-12
    np.testing.assert_allclose(contact_torque[:, 0].sum(), 0.0, atol=0.0)
    np.testing.assert_allclose(result.residual_torque_nm, [0.0, 0.0], atol=1e-12)


def test_contact_at_capacity_yields_and_leaves_excess_closing_motion():
    result = solve_quasistatic_contacts(
        np.eye(2),
        np.array([1.0, -1.0]),
        np.array([[-1.0], [1.0]]),
        np.array([[1.0, -1.0]]),
        np.array([0.5]),
    )

    assert result.converged
    assert result.contact_force_n_per_m[0] == 0.5
    np.testing.assert_allclose(result.omega, [0.5, -0.5], atol=1e-12)
    assert result.closing_rate[0] > 0.0
    np.testing.assert_allclose(result.residual_torque_nm, [0.0, 0.0], atol=1e-12)


def test_si_scale_contact_compliance_is_not_lost_to_absolute_epsilon():
    result = solve_quasistatic_contacts(
        np.eye(2) * 1e27,
        np.array([1e21, -1e21]),
        np.array([[-1e13], [1e13]]),
        np.array([[2e-4, -2e-4]]),
        np.array([5e12]),
        rate_tolerance=1e-18,
    )

    assert result.converged
    assert 0.0 < result.contact_force_n_per_m[0] < 5e12
    assert abs(result.closing_rate[0]) <= 1e-18

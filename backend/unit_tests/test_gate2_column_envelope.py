"""Analytical checks for the offline RFC 0010 Gate 2 column diagnostic."""

import math
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin" / "debug"))

sys.path.insert(0, str(ROOT / "backend"))
from diagnose_gate2_column_envelope import (  # noqa: E402
    FLOW_LAWS,
    GAS_CONSTANT_J_MOL_K,
    Scenario,
    _additional_brittle_extension_pa,
    _pore_pressure_profile,
    creep_stress_pa,
)


def _scenario(**updates):
    values = {
        "base_temperature_c": 1330.0,
        "strain_rate_s": 1e-15,
        "pore_pressure_ratio": 0.0,
        "k0": 1.0,
        "cohesion_pa": 0.0,
        "friction_angle_deg": 30.0,
        "lower_crust_law": "plagioclase_an75",
        "mantle_law": "dry_olivine",
    }
    values.update(updates)
    return Scenario(**values)


def test_hydrostatic_pore_pressure_uses_water_density_and_depth():
    depth = np.array([0.0, 10_000.0, 100_000.0])
    lithostatic = np.array([0.0, 2700.0, 3200.0]) * 9.81 * depth
    result = _pore_pressure_profile(
        depth, lithostatic, _scenario(pore_pressure_mode="hydrostatic_water")
    )
    expected = 1000.0 * 9.81 * depth
    np.testing.assert_allclose(result, expected)
    assert result[0] == 0.0
    assert result[1] / lithostatic[1] == pytest.approx(1000.0 / 2700.0)
    assert result[2] / lithostatic[2] < result[1] / lithostatic[1]


def test_fixed_lithostatic_ratio_is_separate_from_hydrostatic_water():
    depth = np.array([10_000.0, 100_000.0])
    lithostatic = np.array([2700.0, 3200.0]) * 9.81 * depth
    result = _pore_pressure_profile(
        depth, lithostatic,
        _scenario(pore_pressure_ratio=0.42, pore_pressure_mode="fixed_lithostatic_ratio"),
    )
    np.testing.assert_allclose(result, 0.42 * lithostatic)


def test_crust_hydrostatic_pressure_is_held_at_moho_into_mantle():
    depth = np.array([0.0, 10_000.0, 35_000.0, 100_000.0])
    lithostatic = np.array([0.0, 2700.0, 2800.0, 3200.0]) * 9.81 * depth
    scenario = _scenario(pore_pressure_mode="hydrostatic_crust_moho_hold")
    result = _pore_pressure_profile(depth, lithostatic, scenario, crust_thickness_m=35_000.0)
    expected = 1000.0 * 9.81 * np.array([0.0, 10_000.0, 35_000.0, 35_000.0])
    np.testing.assert_allclose(result, expected)


def test_normal_fault_increment_matches_closed_form_for_k0_one():
    vertical_effective = np.array([100e6])
    result = _additional_brittle_extension_pa(
        vertical_effective, k0=1.0, cohesion_pa=0.0, friction_angle_deg=30.0
    )
    expected = 2.0 * math.sin(math.radians(30.0)) / (1.0 + math.sin(math.radians(30.0))) * 100e6
    np.testing.assert_allclose(result, [expected])


def test_normal_fault_increment_includes_cohesion_and_never_goes_negative():
    result = _additional_brittle_extension_pa(
        np.array([0.0, 100e6]), k0=0.6, cohesion_pa=20e6, friction_angle_deg=30.0
    )
    assert np.all(result >= 0.0)
    expected_at_depth = (
        0.6 * 100e6
        - (100e6 * (1.0 - math.sin(math.radians(30.0)))
           - 2.0 * 20e6 * math.cos(math.radians(30.0)))
        / (1.0 + math.sin(math.radians(30.0)))
    )
    assert result[1] == pytest.approx(expected_at_depth)


@pytest.mark.parametrize("law_name", ["dry_olivine", "wet_olivine"])
def test_creep_law_inversion_recovers_stress_that_generates_rate(law_name):
    law = FLOW_LAWS[law_name]
    pressure_pa = np.array([1.2e9])
    temperature_k = np.array([1100.0])
    stress_mpa = 75.0
    log_rate = (
        math.log(law.a_mpa_n_s)
        + law.n * math.log(stress_mpa)
        + law.water_exponent * math.log(max(law.water_ppm_h_si, 1.0))
        - (law.activation_j_mol + pressure_pa[0] * law.activation_volume_m3_mol)
        / (GAS_CONSTANT_J_MOL_K * temperature_k[0])
    )
    recovered = creep_stress_pa(
        law, math.exp(log_rate), pressure_pa, temperature_k
    )
    np.testing.assert_allclose(recovered, [stress_mpa * 1e6], rtol=1e-12)

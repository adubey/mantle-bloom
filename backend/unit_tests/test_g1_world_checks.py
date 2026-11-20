"""Tests for bin/debug/g1_world_checks.py: alternative D's static bathymetry and the ledger read-out (#980)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from app import lithosphere
from app.elevation_lines import CRUST_TYPE_INHERIT

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "debug" / "g1_world_checks.py"
spec = importlib.util.spec_from_file_location("g1_world_checks", SCRIPT)
assert spec is not None and spec.loader is not None
checks = importlib.util.module_from_spec(spec)
sys.modules["g1_world_checks"] = checks
spec.loader.exec_module(checks)

MYR = 1e6
AREA_M2 = 1e10


class _Plate:
    def __init__(self, crust_type: str, hc, hm, created_years, elevation):
        self.crust_type = crust_type
        self.fields = {"crustal_thickness_m": np.asarray(hc, float), "mantle_lithosphere_thickness_m": np.asarray(hm, float),
                       "node_created_years": np.asarray(created_years, float),
                       "crust_type_code": np.full(len(hc), CRUST_TYPE_INHERIT, dtype=np.int8)}
        self.elevation = np.asarray(elevation, float)

    def node_count(self) -> int:
        return len(self.elevation)

    def collect(self, name: str) -> np.ndarray:
        return self.fields[name]

    def accounting_areas_m2(self, spacing_rad: float) -> np.ndarray:
        return np.full(len(self.elevation), AREA_M2)

    def all_points_and_elevation(self):
        return np.zeros((len(self.elevation), 3)), self.elevation


def _world(plates, **fields):
    return SimpleNamespace(plates=plates, node_density=1.0, sea_level_m=0.0, elapsed_years=200 * MYR, age_force_column=False,
                           rift_line_params={}, **fields)


def _ocean(hm_old: float = 40e3):
    # A 1 Myr-old ridge cell (its age column, 13 km less 11 km of crust, sits below its own 8 km) and a 100 Myr cell.
    hc, hm = [11e3, 11e3], [8e3, hm_old]
    elevation = lithosphere.isostatic_elevation(np.array(hc), np.array(hm), lithosphere.RHO_OCEANIC_CRUST)
    return _Plate("oceanic", hc, hm, [199 * MYR, 100 * MYR], elevation)


def test_sea_level_for_volume_inverts_the_water_volume():
    elevation = np.array([-5000.0, -3000.0, -1000.0, 500.0])
    area = np.array([1.0, 2.0, 1.0, 3.0]) * AREA_M2
    water = float(np.dot(area, np.clip(-200.0 - elevation, 0.0, None)))
    assert checks.sea_level_for_volume(elevation, area, water) == pytest.approx(-200.0, abs=1e-6)


def test_alternative_d_deepens_old_floor_by_its_age_column_and_keeps_the_water():
    land = _Plate("continental", [35e3], [100e3], [-1.0], [800.0])
    result = checks.alternative_d(_world([_ocean(), land]))
    # Only the 100 Myr cell moves: Hm 40 -> 84 km (min(13 sqrt(100), 95) - 11), submerged, so it drops by the
    # water-loaded Airy factor.
    drop = 44e3 * (lithosphere.RHO_LITHOSPHERE_MANTLE / lithosphere.RHO_ASTHENOSPHERE - 1.0) * (
        lithosphere.RHO_ASTHENOSPHERE / (lithosphere.RHO_ASTHENOSPHERE - lithosphere.RHO_WATER))
    assert result["oceanic_area_share_deepened"] == pytest.approx(0.5)
    assert result["oceanic_hm_mean_d_km"] == pytest.approx((8.0 + 84.0) / 2)
    # Two ocean cells share the deeper basin's extra drop, so sea level falls by half of it.
    assert result["sea_level_change_m"] == pytest.approx(-drop / 2, rel=1e-6)
    old = result["depth_by_age"][-1]
    assert old["mean_age_myr"] == pytest.approx(100.0)
    assert old["depth_d_m"] - old["depth_now_m"] == pytest.approx(drop / 2, rel=1e-6)
    assert result["land_fraction_d"] == result["land_fraction_now"] == pytest.approx(1 / 3)


def test_alternative_d_leaves_floor_already_thicker_than_its_age_column():
    result = checks.alternative_d(_world([_ocean(hm_old=90e3)]))
    assert result["oceanic_area_share_deepened"] == 0.0
    assert result["sea_level_change_m"] == pytest.approx(0.0, abs=1e-6)


def test_gdh1_is_continuous_at_20_myr():
    assert checks.gdh1_depth_m(0.0) == 2600.0
    assert checks.gdh1_depth_m(19.999) == pytest.approx(checks.gdh1_depth_m(20.0), abs=15.0)


def test_ledger_gaps_report_missing_ledgers_without_creating_them():
    world = _world([_ocean()])
    assert checks.ledger_gaps(world) == {"continental_material": None, "craton": None, "mobile_cover": None, "hm": None}
    assert not hasattr(world, "craton_ledger") and not hasattr(world, "mobile_cover_ledger")


def test_ledger_gaps_measure_the_craton_error_against_its_sources():
    world = _world([_ocean()], craton_ledger={"initial_m3": 100.0, "formed_m3": 0.0})
    # The ocean plate holds no craton, so all 100 m3 are unbooked.
    world.plates[0].fields["craton_crust_m"] = np.zeros(2)
    assert checks.ledger_gaps(world)["craton"] == pytest.approx(1.0)


def test_pool_takes_each_ledgers_worst_gap_and_pools_the_last_samples_bins():
    world = _world([_ocean()])
    one = checks.sample(world)
    one["ledgers"].update(craton=1e-15, continental_material=2e-15)
    two = checks.sample(world)
    two["ledgers"].update(craton=3e-15, continental_material=None)
    runs = [{"samples": [{"myr": 10.0, "world_checks": one}, {"myr": 20.0, "world_checks": two}]},
            {"samples": [{"myr": 10.0}]}]
    pooled = checks.pool(runs)
    assert pooled["ledger_gap_max"]["craton"] == 3e-15
    assert pooled["ledger_gap_max"]["continental_material"] == 2e-15
    assert pooled["ledger_gap_max"]["hm"] is None
    final = pooled["alternative_d_final"]
    assert final["worlds"] == 1 and final["myr"] == [20.0]
    for pooled_bin, own in zip(final["depth_by_age"], two["alternative_d"]["depth_by_age"], strict=True):
        assert pooled_bin.pop("age_myr") == own["age_myr"]
        assert pooled_bin == pytest.approx({k: v for k, v in own.items() if k != "age_myr"})
    assert checks.pool([{"samples": [{"myr": 10.0}]}]) is None

"""Saved-history summaries must include every measured loss and its physical interval."""

import importlib.util
from pathlib import Path

import pytest


_SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "debug" / "ledger_history.py"
_SPEC = importlib.util.spec_from_file_location("ledger_history", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
ledger_history = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ledger_history)


def _row(myr, volume=100.0, craton=None, material=None):
    continental = {
        "surface_continental_derived_m3": 1e18,
        "continental_sediment_on_oceanic_hosts_m3": 0.0,
        "juvenile_additions_m3": 0.0,
        "remelted_relaminated_returns_m3": 0.0,
        "deeply_subducted_m3": 0.0,
        "rift_thinned_m3": 0.0,
        "collision_subducted_m3": 0.0,
        "collision_lower_crust_subducted_m3": 0.0,
        "delaminated_lower_crust_m3": 0.0,
        "overloaded_root_delaminated_m3": 0.0,
        "topology_removed_m3": 0.0,
        "numerical_unplaced_m3": 0.0,
        "discarded_marine_sediment_m3": 0.0,
    }
    continental.update(material or {})
    return {
        "elapsed_years": myr * 1e6,
        "craton_volume_km3": volume,
        "craton_area_km2": 10.0,
        "craton_continental_fraction": 0.5,
        "land_fraction": 0.3,
        "craton_ledger_km3": craton or {},
        "continental_material_ledger": continental,
    }


def test_gyr_windows_include_endpoint_losses_without_counting_them_twice():
    history = [
        _row(0),
        _row(500, 95.0, {"delaminated_m3": 5.0}),
        _row(1000, 80.0, {"delaminated_m3": 20.0}),
        _row(1500, 70.0, {"delaminated_m3": 30.0}),
        _row(2000, 50.0, {"delaminated_m3": 50.0}),
    ]

    windows = ledger_history.summarize(history)["craton_loss_per_gyr"]

    assert [row["window_myr"] for row in windows] == [[0.0, 1000.0], [1000.0, 2000.0]]
    assert [row["delaminated_m3"] for row in windows] == pytest.approx([20.0 / 100.0, 30.0 / 80.0])
    assert [row["total"] for row in windows] == pytest.approx([0.2, 0.375])


def test_partial_final_gyr_uses_measured_duration_for_sinks_and_rift_paths():
    history = [
        _row(0),
        _row(1000, 80.0, {"rifted_m3": 20.0, "rifted_failed_m3": 20.0}),
        _row(1500, 70.0, {"rifted_m3": 30.0, "rifted_failed_m3": 30.0}),
    ]

    first, partial = ledger_history.summarize(history)["craton_loss_per_gyr"]

    assert first["rifted_m3"] == pytest.approx(0.2)
    assert partial["window_myr"] == [1000.0, 1500.0]
    assert partial["rifted_m3"] == pytest.approx(10.0 / 80.0 / 0.5)
    assert partial["rifted_failed_m3"] == pytest.approx(10.0 / 80.0 / 0.5)
    assert partial["total"] == pytest.approx(0.25)


def test_gyr_rate_uses_physical_samples_when_history_starts_after_zero():
    history = [_row(100), _row(600, 90.0, {"subducted_m3": 10.0})]

    window, = ledger_history.summarize(history)["craton_loss_per_gyr"]

    assert window["window_myr"] == [100.0, 600.0]
    assert window["subducted_m3"] == pytest.approx(10.0 / 100.0 / 0.5)


def test_gyr_total_includes_audit_topology_and_signed_residual_accounts():
    history = [
        _row(0),
        _row(1000, 76.0, {
            "rifted_m3": 10.0,
            "topology_removed_m3": 5.0,
            "audit_clipped_m3": 10.0,
            "unattributed_m3": -1.0,
        }),
    ]

    window, = ledger_history.summarize(history)["craton_loss_per_gyr"]

    assert window["topology_removed_m3"] == pytest.approx(0.05)
    assert window["audit_clipped_m3"] == pytest.approx(0.1)
    assert window["unattributed_m3"] == pytest.approx(-0.01)
    assert window["total"] == pytest.approx(0.24)


def test_material_flux_counts_numerical_and_marine_terminal_sinks():
    history = [
        _row(0, material={"numerical_unplaced_m3": 5e15, "discarded_marine_sediment_m3": 7e15}),
        _row(100, material={
            "deeply_subducted_m3": 1e16,
            "numerical_unplaced_m3": 2.5e16,
            "discarded_marine_sediment_m3": 3.7e16,
        }),
    ]

    flux, = ledger_history.summarize(history)["continental_fluxes"]

    # Cumulative starting balances are excluded: 10 + 20 + 30 million km^3 / 100 Myr.
    assert flux["sinks_km3_per_yr"] == pytest.approx(0.6)
    assert flux["deeply_subducted_km3_per_yr"] == pytest.approx(0.1)


def test_old_save_without_new_terminal_accounts_treats_them_as_zero():
    history = [_row(0), _row(100, material={"deeply_subducted_m3": 1e16})]
    for row in history:
        del row["continental_material_ledger"]["numerical_unplaced_m3"]
        del row["continental_material_ledger"]["discarded_marine_sediment_m3"]

    flux, = ledger_history.summarize(history)["continental_fluxes"]

    assert flux["sinks_km3_per_yr"] == pytest.approx(0.1)

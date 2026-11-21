"""Regression tests for the D4 lag diagnostic (#1005)."""

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from app import world as world_module


def _geometry(x, glacier_depth):
    points = np.column_stack((x, np.zeros((len(x), 2))))
    return SimpleNamespace(
        node_count=lambda: len(x),
        all_points_and_elevation=lambda: (points, np.zeros(len(x))),
        collect=lambda _: np.asarray(glacier_depth),
    )


def _hydrology(glacier_depth, lake_depth):
    n = len(glacier_depth)
    return SimpleNamespace(
        glacier_depth=np.asarray(glacier_depth), lake_depth=np.asarray(lake_depth),
        hypsometry=SimpleNamespace(area_m2=np.ones(n)), neighbor_idx=np.arange(n)[:, None],
        lake_water_ledger=SimpleNamespace(glacier_transfer_m3=0.0, residual_m3=0.0),
        land_water_ledger=SimpleNamespace(residual_m3=0.0),
        ice=SimpleNamespace(ledger=SimpleNamespace(residual_m3=0.0)),
    )


@pytest.mark.parametrize("cached", [False, True], ids=["compact-save", "cached-save"])
def test_d4_lag_uses_previous_step_ice_and_positions(monkeypatch, tmp_path, cached):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / "bin" / "debug"))
    diagnostic = importlib.import_module("measure_d4_lake_lag")
    world = SimpleNamespace(
        plates=[_geometry([0, 10, 30], [20, 0, 0])],
        hydrology_cache=_hydrology([20, 0, 0], [1, 0, 1000]) if cached else None,
        seed=1, ice_age_period_years=0, elapsed_years=0,
    )
    steps = iter([
        (_geometry([0, 10, 15, 30], [20, 20, 0, 0]), _hydrology([20, 20, 0, 0], [2, 3, 0, 1000])),
        (_geometry([30, 15, 10, 0], [0, 20, 20, 20]), _hydrology([0, 20, 20, 20], [1000, 5, 3, 2])),
    ])

    def step_world(world, years):
        plate, world.hydrology_cache = next(steps)
        world.plates = [plate]

    result = {"backend": "", "checkpoints": []}
    monkeypatch.setattr(world_module, "step_world", step_world)
    monkeypatch.setattr(diagnostic, "load_replay", lambda *_: (world, result))
    monkeypatch.setattr(sys, "argv", [
        "measure_d4_lake_lag.py", "test.bloom", "--myr", "0.5", "--step-years", "250000",
        "--out", str(tmp_path / "replay.json"),
    ])

    diagnostic.main()

    assert [record["newly_covered_liquid_m3"] for record in result["checkpoints"]] == [3.0, 5.0]

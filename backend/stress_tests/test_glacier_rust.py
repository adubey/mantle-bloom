"""Density-1 same-opening-state parity checks for `mantle_bloom_rust` (#1036)."""

from __future__ import annotations

import copy

import numpy as np

from app import world as world_mod
from app.hydroclimate import _glacier_rust, ice_sheet


def _step_from_world(world, monkeypatch, *, use_rust: bool):
    monkeypatch.setattr(_glacier_rust, "RUST_AVAILABLE", use_rust)
    melt_total = 0.0
    last_closing = 0.0
    orig_step_ice = ice_sheet.step_ice

    def _capture(*args, **kwargs):
        nonlocal melt_total, last_closing
        depth, melt_m, ice = orig_step_ice(*args, **kwargs)
        melt_total += float(np.dot(melt_m, args[2]))
        last_closing = float(ice.ledger.closing_m3)
        return depth, melt_m, ice

    monkeypatch.setattr(ice_sheet, "step_ice", _capture)
    w = copy.deepcopy(world)
    world_mod.step_world(w, 100_000.0)
    march_log = dict(ice_sheet.MARCH_LOG[-1])
    return {
        "substeps": int(march_log["substeps"]),
        "sea_level": float(w.sea_level_m),
        "closing_m3": last_closing,
        "melt_m3": melt_total,
    }


def test_density_one_same_opening_state_parity_seed1(monkeypatch):
    monkeypatch.setattr(_glacier_rust, "RUST_AVAILABLE", False)
    opening = world_mod.generate_world(seed=1, node_density=1.0, climate_density=1.0, fluid_density=1.0)
    opening.wind_model = "diagnostic"

    py = _step_from_world(opening, monkeypatch, use_rust=False)
    ru = _step_from_world(opening, monkeypatch, use_rust=True)

    assert ru["substeps"] == py["substeps"]
    assert abs(ru["sea_level"] - py["sea_level"]) < 1e-2
    assert abs(ru["closing_m3"] - py["closing_m3"]) / max(py["closing_m3"], 1.0) < 1e-3
    assert abs(ru["melt_m3"] - py["melt_m3"]) / max(py["melt_m3"], 1.0) < 1e-5


def test_density_one_seed0_matches_python(monkeypatch):
    monkeypatch.setattr(_glacier_rust, "RUST_AVAILABLE", False)
    opening = world_mod.generate_world(seed=0, node_density=1.0, climate_density=1.0, fluid_density=1.0)
    opening.wind_model = "diagnostic"

    py = _step_from_world(opening, monkeypatch, use_rust=False)
    ru = _step_from_world(opening, monkeypatch, use_rust=True)

    assert ru["substeps"] == py["substeps"]
    assert abs(ru["sea_level"] - py["sea_level"]) < 1e-5
    assert abs(ru["closing_m3"] - py["closing_m3"]) / max(py["closing_m3"], 1.0) < 1e-8
    assert abs(ru["melt_m3"] - py["melt_m3"]) / max(py["melt_m3"], 1.0) < 1e-8

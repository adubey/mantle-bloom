"""Density-1 same-opening-state parity checks for `mantle_bloom_rust` (#1036)."""

from __future__ import annotations

import copy
import math

import numpy as np

from app import world as world_mod
from app.hydroclimate import _glacier_rust, ice_sheet


class _FsumNpProxy:
    """Proxy for `ice_sheet.np` replacing 1-D `np.dot` with exact IEEE `math.fsum`."""

    @staticmethod
    def dot(a, b):
        if np.ndim(a) == 1 and np.ndim(b) == 1:
            return math.fsum((np.asarray(a, dtype=float) * np.asarray(b, dtype=float)).tolist())
        return np.dot(a, b)

    def __getattr__(self, name: str):
        return getattr(np, name)


def _step_from_world(world, monkeypatch, *, use_rust: bool, use_fsum: bool = False):
    monkeypatch.setattr(_glacier_rust, "RUST_AVAILABLE", use_rust)
    if use_fsum:
        monkeypatch.setattr(ice_sheet, "np", _FsumNpProxy())
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


def test_density_one_seed0_matches_exact_fsum_python(monkeypatch):
    monkeypatch.setattr(_glacier_rust, "RUST_AVAILABLE", False)
    opening = world_mod.generate_world(seed=0, node_density=1.0, climate_density=1.0, fluid_density=1.0)
    opening.wind_model = "diagnostic"

    py_fsum = _step_from_world(opening, monkeypatch, use_rust=False, use_fsum=True)
    ru = _step_from_world(opening, monkeypatch, use_rust=True)

    assert ru["substeps"] == py_fsum["substeps"]
    assert abs(ru["sea_level"] - py_fsum["sea_level"]) < 1e-5
    assert abs(ru["closing_m3"] - py_fsum["closing_m3"]) / max(py_fsum["closing_m3"], 1.0) < 1e-8
    assert abs(ru["melt_m3"] - py_fsum["melt_m3"]) / max(py_fsum["melt_m3"], 1.0) < 1e-8

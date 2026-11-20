"""Tests for bin/debug/calibrate_plate_forces.py's setting overrides and its resume check (#975, #979 review)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from app import torque

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "debug" / "calibrate_plate_forces.py"
spec = importlib.util.spec_from_file_location("calibrate_plate_forces", SCRIPT)
assert spec is not None and spec.loader is not None
calibrate = importlib.util.module_from_spec(spec)
sys.modules["calibrate_plate_forces"] = calibrate
spec.loader.exec_module(calibrate)


def test_module_knobs_set_known_constants_and_refuse_bad_values(monkeypatch):
    monkeypatch.setattr(torque, "SUBDUCTION_LSINK_M", torque.SUBDUCTION_LSINK_M)  # restores the constant afterwards
    calibrate.apply_module_knobs({"torque.SUBDUCTION_LSINK_M": 150e3})
    assert torque.SUBDUCTION_LSINK_M == 150e3
    for bad in ({"torque.SLAB_MANTLE_VISCOSITY_PA_S": -1.0}, {"torque.SLAB_MANTLE_VISCOSITY_PA_S": 0.0},
                {"torque.SUBDUCTION_LSINK_M": float("nan")}, {"torque.SUBDUCTION_LSINK_M": float("inf")},
                {"torque.GRAVITY": 1.0}, {"nowhere.SUBDUCTION_LSINK_M": 1.0}):
        with pytest.raises(SystemExit):
            calibrate.apply_module_knobs(bad)
    assert torque.SUBDUCTION_LSINK_M == 150e3


def _settings_file(tmp_path, module_knobs=None):
    path = tmp_path / "settings.json"
    setting = {"params": {}, "world_knobs": {"slab_drag_form": "trench_normal"}}
    if module_knobs:
        setting["module_knobs"] = module_knobs
    path.write_text(json.dumps({"common_params": {}, "settings": {"s": setting}}))
    return path


def _saved_run(tmp_path, settings, **override):
    saved = {"seed": 1, "step_myr": 1.0, "every": 10.0, "years_myr": 200.0, "setting": calibrate.load_setting(settings, "s")}
    saved.update(override)
    path = tmp_path / "run.json"
    path.write_text(json.dumps(saved))
    return path


def _args(settings, **override):
    values = {"settings": settings, "setting": "s", "step_myr": 1.0, "every": 10.0, "years": 200.0, "no_contacts": False}
    return SimpleNamespace(**{**values, **override})


def test_resume_refuses_a_finished_run_made_for_another_experiment(tmp_path):
    settings = _settings_file(tmp_path)
    args = _args(settings)
    path = _saved_run(tmp_path, settings)
    calibrate.check_cached_run(path, args)  # the same experiment: reused
    with pytest.raises(SystemExit, match="step_myr"):
        calibrate.check_cached_run(path, _args(settings, step_myr=0.5))
    other = _settings_file(tmp_path, module_knobs={"torque.SUBDUCTION_LSINK_M": 150e3})
    with pytest.raises(SystemExit, match="setting"):
        calibrate.check_cached_run(path, _args(other))
    with pytest.raises(SystemExit, match="years_myr"):
        calibrate.check_cached_run(path, _args(settings, years=100.0))
    with pytest.raises(SystemExit, match="contacts"):
        calibrate.check_cached_run(path, _args(settings, no_contacts=True))
    calibrate.check_cached_run(_saved_run(tmp_path, settings, contacts=False), _args(settings, no_contacts=True))

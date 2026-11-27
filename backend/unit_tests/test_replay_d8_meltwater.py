"""Tests for the G6 ice-margin network diagnostic."""

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


def test_only_margin_and_adjacent_headwaters_count(monkeypatch):
    debug_dir = Path(__file__).resolve().parents[2] / "bin" / "debug"
    monkeypatch.syspath_prepend(str(debug_dir))
    spec = importlib.util.spec_from_file_location("replay_d8_meltwater", debug_dir / "replay_d8_meltwater.py")
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    spec.loader.exec_module(module)
    hydro = SimpleNamespace(
        glacier_depth=np.array([10.0, 10.0, 0.0, 0.0]),
        neighbor_idx=np.array([[1, 1], [0, 2], [1, 3], [2, 2]]),
        flow_target=np.full(4, -1),
    )
    rivers = [SimpleNamespace(member_idx=np.array([i]), flow_rate=float(i + 1)) for i in range(4)]
    hydrology = SimpleNamespace(GLACIER_VISIBLE_DEPTH_M=1.0, group_rivers=lambda _: rivers)

    assert module._glacier_rivers(hydro, hydrology) == [3.0, 2.0]

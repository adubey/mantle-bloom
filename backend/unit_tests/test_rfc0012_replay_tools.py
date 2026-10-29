import argparse
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_script(name: str):
    path = REPO_ROOT / "bin" / "debug" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_checkpoint_trajectory_handles_all_undefined_quantile_deltas():
    summarize = _load_script("summarize_rfc0012_merge_replay.py")
    row = {
        "myr": 10.0,
        "plate_count": 2,
        "land_fraction": 0.5,
        "hc_area_weighted_km": {name: 0.0 for name in summarize.QUANTILES},
        "hm_area_weighted_km": {name: 0.0 for name in summarize.QUANTILES},
    }
    changed = {
        **row,
        "hc_area_weighted_km": {name: 1.0 for name in summarize.QUANTILES},
        "hm_area_weighted_km": {name: 1.0 for name in summarize.QUANTILES},
    }

    result = summarize._checkpoint_trajectory({"checkpoints": [row]}, {"checkpoints": [changed]})

    assert result["max_abs_hc_hm_quantile_delta_percent"] == [None]


def test_resume_rolls_metrics_back_to_saved_checkpoint_and_enables_diagnostics(tmp_path, monkeypatch):
    replay = _load_script("replay_rfc0012_merge_rate.py")
    world_path = tmp_path / "world.mbworld"
    world_path.write_bytes(b"test")
    world = SimpleNamespace(elapsed_years=30e6, phase_budget={})
    monkeypatch.setattr(replay.persistence, "load_world_bytes", lambda _: world)
    checkpoint = {
        "seed": 7,
        "step_myr": 1.0,
        "world_save_path": str(world_path),
        "world_save_myr": 30.0,
        "checkpoints": [
            {"myr": 0.0, "hm_inventory_m3": {"all": 100.0}},
            {
                "myr": 30.0,
                "deformation_since_checkpoint": {"merge": {"hm_m3": 3.0}},
                "motion_cumulative": {"rate_samples": 30},
                "runtime_seconds": 12.0,
                "peak_rss_mib": 200.0,
            },
            {
                "myr": 40.0,
                "deformation_since_checkpoint": {"merge": {"hm_m3": 4.0}},
                "motion_cumulative": {"rate_samples": 40},
            },
        ],
        "topology_events": [[29e6, "merge"], [35e6, "split"]],
        "motion": {"rate_samples": 40},
        "deformation": {"merge": {"hm_m3": 7.0}},
    }
    checkpoint_path = tmp_path / "resume.json"
    checkpoint_path.write_text(json.dumps(checkpoint))
    args = SimpleNamespace(seed=7, step_myr=1.0)

    resumed_world, result = replay._load_resume(checkpoint_path, args)

    assert resumed_world.debug_diagnostics is True
    assert [row["myr"] for row in result["checkpoints"]] == [0.0, 30.0]
    assert result["motion"]["rate_samples"] == 30
    assert result["deformation"]["merge"]["hm_m3"] == 3.0
    assert result["topology_events"] == [[29e6, "merge"]]
    assert result["runtime_seconds"] == 12.0


def test_world_save_must_align_to_checkpoint_boundary():
    replay = _load_script("replay_rfc0012_merge_rate.py")
    parser = argparse.ArgumentParser()
    args = SimpleNamespace(world_save_myr=35.0, checkpoint_myr=10.0, step_myr=1.0, myr=50.0)

    with pytest.raises(SystemExit):
        replay._validate_intervals(parser, args)


def test_solve_contact_metrics_are_paired_with_prior_step_topology_events():
    replay = _load_script("replay_rfc0012_merge_rate.py")
    events = [[10e6, "merge"], [20e6, "split"]]

    result = replay._topology_events_for_previous_step(events, elapsed_years=20e6, step_years=10e6)

    assert result == ["merge"]


@pytest.mark.parametrize(("platform", "raw", "expected"), [("linux", 1024.0, 1.0), ("darwin", 2**20, 1.0)])
def test_peak_rss_conversion_is_platform_specific(monkeypatch, platform, raw, expected):
    replay = _load_script("replay_rfc0012_merge_rate.py")
    monkeypatch.setattr(replay.sys, "platform", platform)
    monkeypatch.setattr(replay.resource, "getrusage", lambda _: SimpleNamespace(ru_maxrss=raw))

    assert replay._peak_rss_mib() == expected

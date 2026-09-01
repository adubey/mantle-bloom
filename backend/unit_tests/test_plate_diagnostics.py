import json

import numpy as np
import pytest
from app import persistence
from app.plate_diagnostics import build_report, clean_tiling_node_estimate, format_report, main
from app.world import generate_world, step_world


@pytest.fixture
def stepped_world():
    world = generate_world(seed=7, num_plates=6)
    step_world(world, 5_000_000)
    step_world(world, 5_000_000)
    return world


def test_clean_tiling_estimate_matches_the_known_density_4_figure():
    # docs/TODO.md records ~130k for a clean tiling at node_density 4, ~32.6k at 1x.
    assert clean_tiling_node_estimate(4.0) == pytest.approx(130_000, rel=0.02)
    assert clean_tiling_node_estimate(1.0) == pytest.approx(32_600, rel=0.02)


def test_build_report_has_one_row_per_plate_and_a_node_budget(stepped_world):
    report = build_report(stepped_world)
    assert report["seed"] == 7
    assert len(report["plates"]) == len(stepped_world.plates)
    assert report["approx_steps"] == round(stepped_world.elapsed_years / 100_000)

    budget = report["node_budget"]
    expected_total = sum(p.node_count() for p in stepped_world.plates)
    assert budget["total_nodes"] == expected_total
    assert budget["clean_tiling_estimate"] == round(clean_tiling_node_estimate(stepped_world.node_density))
    assert budget["ratio"] == pytest.approx(expected_total / clean_tiling_node_estimate(stepped_world.node_density), abs=5e-4)

    for row in report["plates"]:
        assert row["crust_type"] in ("continental", "oceanic")
        assert row["num_points"] >= 0
        assert 0.0 <= row["submerged_fraction"] <= 1.0


def test_report_is_json_serializable(stepped_world):
    json.dumps(build_report(stepped_world))  # no numpy scalars / tuples left in


def test_format_report_renders_every_section(stepped_world):
    text = format_report(build_report(stepped_world))
    assert "mantle-bloom plate diagnostics" in text
    assert "per-plate" in text
    assert "territory overlaps" in text
    assert "sustained-collision timers" in text
    assert "node budget" in text
    # one table row per plate
    for plate in stepped_world.plates:
        assert f"\n  {plate.plate_id:>3}  " in text


def test_main_round_trips_a_saved_file(tmp_path, stepped_world, capsys):
    save = tmp_path / "seed7.mbworld"
    save.write_bytes(persistence.save_world_bytes(stepped_world))

    assert main([str(save)]) == 0
    assert "mantle-bloom plate diagnostics" in capsys.readouterr().out

    assert main([str(save), "--json"]) == 0
    parsed = json.loads(capsys.readouterr().out)
    assert parsed["seed"] == 7


def test_main_errors_on_a_missing_file(tmp_path):
    with pytest.raises(SystemExit):
        main([str(tmp_path / "nope.mbworld")])

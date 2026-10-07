import json

import pytest
from app import persistence
from app.phase_budget import SCOPES
from app.phase_budget_diagnostics import CONVENTIONAL_YEARS_PER_STEP, build_report, format_report, main
from app.world import generate_world


@pytest.fixture
def fresh_world():
    return generate_world(seed=11, num_plates=6, node_density=1.0)


def test_build_report_steps_forward_and_turns_on_diagnostics(fresh_world):
    start_years = fresh_world.elapsed_years
    report = build_report(fresh_world, total_years=300_000)

    assert fresh_world.debug_diagnostics is True
    assert report["seed"] == 11
    assert report["start_elapsed_years"] == start_years
    assert report["steps"] == 3
    assert report["end_elapsed_years"] == pytest.approx(start_years + 3 * CONVENTIONAL_YEARS_PER_STEP)
    assert fresh_world.elapsed_years == pytest.approx(report["end_elapsed_years"])


def test_build_report_rows_have_every_scope_and_derived_field(fresh_world):
    report = build_report(fresh_world, total_years=200_000)
    assert report["phases"]  # a two-step run on any seed touches at least the deform phases

    for row in report["phases"]:
        assert row["calls"] > 0
        assert set(row["scopes"]) == set(SCOPES)
        for scope in SCOPES:
            fields = row["scopes"][scope]
            assert fields["delta_count"] == fields["count_after"] - fields["count_before"]
            assert fields["delta_sum_hc"] == pytest.approx(fields["sum_hc_after"] - fields["sum_hc_before"])
            assert fields["delta_sum_hm"] == pytest.approx(fields["sum_hm_after"] - fields["sum_hm_before"])
            assert fields["area_before_m2"] >= 0.0
            assert fields["area_after_m2"] >= 0.0
            assert fields["delta_area_m2"] == pytest.approx(fields["area_after_m2"] - fields["area_before_m2"])
            assert fields["delta_hc_volume_m3"] == pytest.approx(fields["hc_volume_after_m3"] - fields["hc_volume_before_m3"])
            assert fields["delta_hm_volume_m3"] == pytest.approx(fields["hm_volume_after_m3"] - fields["hm_volume_before_m3"])


def test_build_report_sorts_phases_by_descending_hc_impact(fresh_world):
    report = build_report(fresh_world, total_years=500_000)
    deltas = [abs(row["scopes"]["all"]["delta_hc_volume_m3"]) for row in report["phases"]]
    assert deltas == sorted(deltas, reverse=True)


def test_build_report_all_scope_is_the_sum_of_the_two_plate_type_scopes(fresh_world):
    report = build_report(fresh_world, total_years=300_000)
    for row in report["phases"]:
        scopes = row["scopes"]
        for field in ("count_before", "count_after"):
            assert scopes["all"][field] == scopes["continental_plate"][field] + scopes["oceanic_plate"][field]
        for field in ("count_before", "count_after"):
            assert scopes["all"][field] == scopes["continental_node"][field] + scopes["oceanic_node"][field]


def test_report_is_json_serializable(fresh_world):
    json.dumps(build_report(fresh_world, total_years=200_000))


def test_hm_source_sink_identity_closes_and_stays_separate_from_phases(fresh_world):
    report = build_report(fresh_world, total_years=200_000)
    assert report["hm_source_sink_accounts"]
    assert all("hm_cap_transitions" in row for row in report["phases"])
    for row in report["hm_closure"].values():
        assert row["live_change_m3"] == pytest.approx(
            row["source_m3"] - row["sink_m3"] + row["reclassification_m3"] + row["signed_residual_m3"],
            abs=1.0,
        )


def test_format_report_renders_both_tables(fresh_world):
    text = format_report(build_report(fresh_world, total_years=200_000))
    assert "mantle-bloom Hc/Hm phase budget" in text
    assert "continental/oceanic node-type split" in text
    assert "convergent_deformation" in text


def test_main_round_trips_a_saved_file(tmp_path, fresh_world, capsys):
    save = tmp_path / "seed11.mbworld"
    save.write_bytes(persistence.save_world_bytes(fresh_world))

    assert main([str(save), "--years", "200000"]) == 0
    assert "mantle-bloom Hc/Hm phase budget" in capsys.readouterr().out

    assert main([str(save), "--years", "200000", "--json"]) == 0
    parsed = json.loads(capsys.readouterr().out)
    assert parsed["seed"] == 11


def test_main_generates_a_fresh_world_when_no_save_is_given(capsys):
    assert main(["--seed", "99", "--node-density", "1", "--years", "200000"]) == 0
    assert "seed:          99" in capsys.readouterr().out


def test_main_errors_on_a_missing_file(tmp_path):
    with pytest.raises(SystemExit):
        main([str(tmp_path / "nope.mbworld"), "--years", "100000"])


def test_shortening_ledger_is_reported_apart_from_the_phases(fresh_world):
    # Issue #314: the collisional-shortening entry has no Hc/Hm scopes of its own.
    from app.phase_budget import SHORTENING_PHASE

    report = build_report(fresh_world, CONVENTIONAL_YEARS_PER_STEP)
    assert all(row["phase"] != SHORTENING_PHASE for row in report["phases"])
    book = report["shortening"]
    assert book["demand_m2"] > 0.0
    returned = book["returned_decay_m2"] + book["returned_unrouted_m2"]
    assert book["absorbed_m2"] + returned == pytest.approx(book["demand_m2"], rel=1e-12)
    assert "collisional shortening" in format_report(report)

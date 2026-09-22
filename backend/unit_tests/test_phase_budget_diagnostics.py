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


def test_build_report_sorts_phases_by_descending_hc_impact(fresh_world):
    report = build_report(fresh_world, total_years=500_000)
    deltas = [abs(row["scopes"]["all"]["delta_sum_hc"]) for row in report["phases"]]
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


def test_line_growth_shrink_sub_phases_sum_to_the_aggregate():
    # GitHub issue #216: line_growth_shrink bundles five mechanistically distinct sub-events
    # (endpoint stretch-thinning, brand-new arc-margin nodes, plain end/interior deletion, and
    # accreted-column redistribution) inside one call to _grow_or_shrink_line_for_deform. Each
    # now gets its own phase_budget.record() call in addition to the whole call's own net
    # total -- this asserts the five sub-phases actually cover every mutation that total makes.
    #
    # The *net* deltas must telescope to the aggregate's own net delta (each sub-phase records
    # only the specific slice it touches, while the aggregate call records the whole line's
    # count/sum on every call regardless of whether that line changed at all -- so raw
    # count_before/count_after and sum_before/sum_after aren't comparable directly, only the
    # before-after deltas are).
    world = generate_world(seed=23, num_plates=8, node_density=1.0)
    report = build_report(world, total_years=3_000_000)
    rows = {row["phase"]: row for row in report["phases"]}
    assert "line_growth_shrink" in rows

    sub_phases = ["line_end_stretch", "line_end_arc_grow", "line_end_retreat", "line_end_accretion", "line_interior_carve"]
    present = [name for name in sub_phases if name in rows]
    assert present, "expected at least one line_growth_shrink sub-phase to fire over this replay"

    for scope in SCOPES:
        aggregate = rows["line_growth_shrink"]["scopes"][scope]
        aggregate_delta_count = aggregate["count_after"] - aggregate["count_before"]
        aggregate_delta_sum_hc = aggregate["sum_hc_after"] - aggregate["sum_hc_before"]
        aggregate_delta_sum_hm = aggregate["sum_hm_after"] - aggregate["sum_hm_before"]

        combined_delta_count = sum(rows[name]["scopes"][scope]["count_after"] - rows[name]["scopes"][scope]["count_before"] for name in present)
        combined_delta_sum_hc = sum(rows[name]["scopes"][scope]["sum_hc_after"] - rows[name]["scopes"][scope]["sum_hc_before"] for name in present)
        combined_delta_sum_hm = sum(rows[name]["scopes"][scope]["sum_hm_after"] - rows[name]["scopes"][scope]["sum_hm_before"] for name in present)

        # rel tolerance, not abs -- both sides accumulate the same underlying per-node deltas
        # across thousands of calls, just summed in a different grouping/order, so float64
        # accumulation noise scales with the totals themselves (~1e6 m here), not a fixed floor.
        assert combined_delta_count == aggregate_delta_count, scope
        assert combined_delta_sum_hc == pytest.approx(aggregate_delta_sum_hc, rel=1e-6, abs=1.0), scope
        assert combined_delta_sum_hm == pytest.approx(aggregate_delta_sum_hm, rel=1e-6, abs=1.0), scope

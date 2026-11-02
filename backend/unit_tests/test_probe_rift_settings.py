"""probe_rift_settings.py -- Gate D1's fragment merging, persistence and event accounting (#737)."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin" / "debug"))

import probe_rift_settings as probe  # noqa: E402


def _arc(start_deg, stop_deg, count):
    angles = np.radians(np.linspace(start_deg, stop_deg, count))
    return np.column_stack((np.cos(angles), np.sin(angles), np.zeros(count)))


def _tracker():
    return probe.Tracker(probe.ONSET_MYR, probe.LAPSE_MYR, np.radians(2.0))


def _event(onset, plate=1, setting="transtension", site=(1.0, 0.0, 0.0), growth=None, last=None):
    row = {"onset_myr": onset, "plate_id": plate, "route_status": "separates", "setting": setting, "site_local": list(site),
           "last_myr": onset + 100.0 if last is None else last}
    if growth is not None:
        row["growth_myr"] = growth
    return row


def test_fragments_facing_one_neighbour_merge_across_a_small_gap():
    points = np.vstack((_arc(0.0, 2.0, 5), _arc(3.0, 5.0, 5)))
    neighbours = np.full(10, 7)
    fronts = probe.merged_fronts(points, neighbours, np.arange(10), np.radians(1.5))
    assert [(neighbour, len(nodes)) for neighbour, nodes in fronts] == [(7, 10)]


def test_fragments_stay_separate_across_a_wide_gap_or_another_neighbour():
    points = np.vstack((_arc(0.0, 2.0, 5), _arc(30.0, 32.0, 5)))
    assert len(probe.merged_fronts(points, np.full(10, 7), np.arange(10), np.radians(1.5))) == 2
    points = _arc(0.0, 4.5, 10)
    neighbours = np.array([7] * 5 + [8] * 5)
    fronts = probe.merged_fronts(points, neighbours, np.arange(10), np.radians(1.5))
    assert {neighbour for neighbour, _nodes in fronts} == {7, 8}


def test_front_length_counts_faced_nodes_so_band_width_does_not_inflate_it():
    # Three rows of own nodes facing the same ten neighbour nodes are one 10-node front.
    state = {"inputs": SimpleNamespace(neighbor_node_index=np.tile(np.arange(10), 3)), "spacing": 1e-3}
    expected = 10 * 1e-3 * probe.lithosphere.PLANET_RADIUS_M / 1e3
    assert probe.front_km(state, np.arange(30)) == expected


def test_front_only_becomes_an_event_after_persisting_for_onset():
    tracker = _tracker()
    local = _arc(0.0, 3.0, 4)
    detection = [(1, "transtension", 2, local)]
    assert tracker.update(1.0, detection, {1})[0] == []
    assert tracker.update(2.0, detection, {1})[0] == []
    started, _landed = tracker.update(3.0, detection, {1})
    assert [track.onset_myr for _index, track in started] == [3.0]
    assert tracker.update(4.0, detection, {1})[0] == []


def test_a_missed_step_before_onset_resets_the_clock():
    tracker = _tracker()
    tracker.update(1.0, [(1, "passive", -1, None)], {1})
    tracker.update(2.0, [], {1})
    started, _landed = tracker.update(3.0, [(1, "passive", -1, None)], {1})
    assert started == []
    assert tracker.closed[0][1] == "lapsed"


def test_passive_persists_per_plate_though_its_site_moves():
    tracker = _tracker()
    for now in (1.0, 2.0):
        tracker.update(now, [(1, "passive", -1, None)], {1})
    started, _landed = tracker.update(3.0, [(1, "passive", -1, None)], {1})
    assert len(started) == 1
    assert len(tracker.open) == 1


def test_a_front_elsewhere_on_the_boundary_is_a_new_track():
    tracker = _tracker()
    tracker.update(1.0, [(1, "rollback", 2, _arc(0.0, 3.0, 4))], {1})
    _started, landed = tracker.update(2.0, [(1, "rollback", 2, _arc(60.0, 63.0, 4))], {1})
    assert landed[0].first_myr == 2.0


def test_a_running_event_bridges_gaps_up_to_the_lapse_and_then_ends():
    tracker = _tracker()
    for now in (1.0, 2.0, 3.0):
        tracker.update(now, [(1, "passive", -1, None)], {1})
    tracker.update(8.0, [], {1})
    assert len(tracker.open) == 1
    tracker.update(9.0, [], {1})
    assert not tracker.open
    track, ended = tracker.closed[0]
    assert (track.last_myr, ended) == (3.0, "lapsed")


def test_a_gone_plate_ends_its_tracks_with_the_reported_cause():
    tracker = _tracker()
    tracker.update(1.0, [(1, "passive", -1, None)], {1})
    tracker.update(2.0, [], {2})
    assert tracker.closed[0][1] == "plate_gone"


def test_a_split_ends_tracks_on_a_parent_that_keeps_its_id():
    tracker = _tracker()
    census = SimpleNamespace(tracker=tracker, production=probe.Counter())
    tracker.update(1.0, [(4, "passive", -1, None)], {4})
    world = SimpleNamespace(events=[(1.0, "Plate 4 split into plates 4 and 9.")])
    probe.Census.note_topology(census, world, 0.0)
    assert tracker.closed[0][1] == "production_split"
    assert census.production["splits"] == 1


def test_a_merge_ends_tracks_on_both_plates():
    tracker = _tracker()
    census = SimpleNamespace(tracker=tracker, production=probe.Counter())
    tracker.update(1.0, [(4, "passive", -1, None), (5, "passive", -1, None)], {4, 5})
    world = SimpleNamespace(events=[(1.0, "Plates 4 and 5 collided and merged into plate 4 (x).")])
    probe.Census.note_topology(census, world, 0.0)
    assert [cause for _track, cause in tracker.closed] == ["plate_merged", "plate_merged"]


def test_gating_allows_one_rift_per_plate_at_a_time():
    kept = probe.gated_events([_event(10.0, setting="passive", growth=15.0),
                               _event(12.0, setting="rollback", growth=15.0),
                               _event(12.0, plate=2, growth=15.0)])
    assert [(row["plate_id"], row["setting"]) for row in kept] == [(1, "passive"), (2, "transtension")]


def test_gating_holds_the_plate_for_the_shorter_of_growth_and_driver_life():
    kept = probe.gated_events([_event(10.0, setting="passive", growth=40.0, last=11.0),
                               _event(20.0, setting="rollback", growth=15.0)])
    assert len(kept) == 2


def test_gating_refuses_a_repeat_near_an_old_site_until_the_retry_passes():
    near = (np.cos(0.01), np.sin(0.01), 0.0)
    far = (0.0, 1.0, 0.0)
    kept = probe.gated_events([_event(10.0, growth=1.0), _event(20.0, growth=1.0, site=near),
                               _event(30.0, growth=1.0, site=far), _event(61.0, growth=1.0, site=near)])
    assert [row["onset_myr"] for row in kept] == [10.0, 30.0, 61.0]


def _report(seed, events, final=20, splits=4):
    return {"seed": seed, "myr": 200.0, "plates": [[0.0, 18], [200.0, final]], "plates_max": final,
            "production": {"splits": splits, "merges": 3}, "pre_onset_tracks": {},
            "polarity_front_steps": {"ocean_ocean:upper": 9, "ocean_ocean:ambiguous": 1,
                                     "continent_ocean:upper": 8, "continent_ocean:lower": 2},
            "polarity_fronts": {"ocean_ocean": 10}, "polarity_flips": {"ocean_ocean": 1}, "events": events}


def _full_event(onset, setting, **extra):
    return {**_event(onset, plate=int(onset), setting=setting, growth=10.0), "first_myr": onset - 2.0,
            "ended": "lapsed", "duration_myr": 12.0, "after_head": False, **extra}


def test_summary_counts_settings_on_two_of_three_seeds_and_splits_rollback_by_upper_plate():
    reports = [_report(1, [_full_event(10.0, "passive"), _full_event(20.0, "rollback", upper_type="oceanic")]),
               _report(2, [_full_event(10.0, "passive"), _full_event(20.0, "rollback", upper_type="continental")]),
               _report(3, [])]
    summary = probe.summarize(reports)
    assert sorted(summary["gate_settings_passing"]) == ["passive", "rollback"]
    assert summary["settings"]["passive"]["events_per_100myr"] == round(2 / 6, 2)
    assert summary["rollback_by_upper_plate"]["oceanic"]["events"] == 1
    assert summary["rollback_by_upper_plate"]["continental"]["events"] == 1
    assert summary["settings"]["passive"]["outlast_growth"] == 2


def test_oceanic_rollback_cannot_carry_the_setting_when_polarity_fails():
    reports = [_report(seed, [_full_event(20.0, "rollback", upper_type="oceanic")]) for seed in (1, 2)]
    for report in reports:
        report["polarity_front_steps"]["ocean_ocean:ambiguous"] = 90
    assert "rollback" not in probe.summarize(reports)["gate_settings_passing"]


def test_summary_applies_the_polarity_thresholds_and_truth_test():
    polarity = probe.polarity_summary([_report(1, [])])
    assert polarity["ocean_ocean"]["decided_share"] == 0.9
    assert polarity["ocean_ocean"]["works"] is True
    assert polarity["continent_ocean"]["ocean_lower_share"] == 0.8


def test_older_ocean_is_lower_only_past_the_age_gap():
    young, old = np.full(3, 150e6), np.full(3, 20e6)
    assert probe.age_verdict(old, young, 200.0) == "lower"
    assert probe.age_verdict(young, old, 200.0) == "upper"
    assert probe.age_verdict(np.full(3, 150e6), np.full(3, 145e6), 200.0) == "tied"


def test_negative_starting_dates_decide_and_sentinels_read_as_rock_from_the_start():
    assert probe.age_verdict(np.full(3, -80e6), np.full(3, -20e6), 200.0) == "lower"
    assert probe.age_verdict(np.full(3, -20e6), np.full(3, -80e6), 200.0) == "upper"
    assert probe.age_verdict(np.zeros(3), np.full(3, -1.0), 200.0) == "tied"
    assert probe.age_verdict(np.zeros(3), np.full(3, 100e6), 200.0) == "lower"


def _census_row(age):
    return {"kind": "ocean_ocean", "neighbour": 2, "nodes": np.arange(4), "counted": True, "age": age, "verdict": "ambiguous"}


def _census_plate():
    plate = SimpleNamespace(plate_id=1, frame=np.eye(3))
    return plate, {"points": _arc(0.0, 3.0, 4)}


def test_a_fronts_age_verdict_freezes_and_a_reversal_is_counted_once():
    census, (plate, state) = probe.Census(np.radians(1.0)), _census_plate()
    for now, age in ((1.0, "tied"), (2.0, "upper"), (3.0, "lower"), (4.0, "lower")):
        census.observe_polarity([(plate, state, _census_row(age))], now, {1: plate, 2: plate})
    assert census.frozen_polarity(plate, state, 2, np.arange(4)) == "upper"
    assert census.age_steps["ocean_ocean", "reversed_fronts"] == 1
    assert census.age_steps["ocean_ocean", "undecided"] == 1


def test_a_front_with_no_record_reads_as_tied():
    census, (plate, state) = probe.Census(np.radians(1.0)), _census_plate()
    assert census.frozen_polarity(plate, state, 2, np.arange(4)) == "tied"


def test_age_read_replaces_consumption_for_the_rollback_rule():
    reports = [_report(seed, [_full_event(20.0, "rollback", upper_type="oceanic")]) for seed in (1, 2)]
    for report in reports:
        report["polarity_front_steps"]["ocean_ocean:ambiguous"] = 90
        report["age_polarity_steps"] = {"ocean_ocean:upper": 9, "ocean_ocean:lower": 9, "ocean_ocean:undecided": 1}
    summary = probe.summarize(reports)
    assert summary["age_polarity"]["works"] is True
    assert "rollback" in summary["gate_settings_passing"]


def test_worst_case_plate_count_replaces_production_splits_with_gated_events():
    events = [_full_event(float(5 + 20 * i), "passive") for i in range(9)]
    row = probe.summarize([_report(1, events, final=36, splits=4)])["plate_count"][0]
    assert row["worst_case_plates"] == 41.0
    assert row["compatible"] is False


class _GridPlate:
    """A square of cells on the sphere, 4-connected, for routing tests."""

    def __init__(self, size, step_deg, fields):
        lat, lon = np.meshgrid(np.radians((np.arange(size) - size // 2) * step_deg),
                               np.radians((np.arange(size) - size // 2) * step_deg), indexing="ij")
        self.points = np.column_stack(((np.cos(lat) * np.cos(lon)).ravel(), (np.cos(lat) * np.sin(lon)).ravel(), np.sin(lat).ravel()))
        index = np.arange(size * size).reshape(size, size)
        links = [[] for _ in range(size * size)]
        for a, b in ((index[:-1], index[1:]), (index[:, :-1], index[:, 1:])):
            for i, j in zip(a.ravel(), b.ravel()):
                links[i].append(j)
                links[j].append(i)
        self._adjacency = SimpleNamespace(offsets=np.r_[0, np.cumsum([len(row) for row in links])],
                                          neighbours=np.array([j for row in links for j in row]))
        self.fields, self.size, self.index = fields, size, index

    def adjacency(self):
        return self._adjacency

    def node_count(self):
        return len(self.points)

    def collect(self, name):
        return self.fields[name]

    def node_areas_m2(self):
        return np.ones(len(self.points))


def _craton_ring_plate(heated, keel_m):
    """Craton everywhere but a small hole around the centre, so every path out crosses its interior."""
    size, step = 41, 0.5
    centre = size // 2
    rows, cols = np.indices((size, size))
    hole = ((rows - centre) ** 2 + (cols - centre) ** 2 <= 9).ravel()
    craton = np.where(hole, 0.0, 2.0 * probe.cratons.CRATON_FULL_STRENGTH_HC_M)
    fields = {"craton_crust_m": craton, "plume_heat_myr": np.full(size * size, 5.0 if heated else 0.0),
              "mantle_lithosphere_thickness_m": np.full(size * size, keel_m),
              "belt_formed_years": np.zeros(size * size)}
    return _GridPlate(size, step, fields), int(centre * size + centre), np.radians(step)


def _route_ring(heated, keel_m):
    plate, seed, spacing = _craton_ring_plate(heated, keel_m)
    cratons_mask = plate.collect("craton_crust_m") > 0.0
    forbidden, _interior = probe.craton_barrier(plate, cratons_mask, spacing)
    resistance = probe.route_resistance(plate, np.zeros(plate.node_count(), dtype=bool), forbidden,
                                        probe.craton_margin(plate, cratons_mask, spacing))
    found, status, _cheapest = probe.route_by_daughters(plate, plate.points, seed, np.array([0.0, 1.0, 0.0]),
                                                        forbidden, resistance, plate.node_areas_m2(), np.radians(1.0))
    return found, status, forbidden


def test_plume_heat_alone_leaves_an_intact_craton_a_barrier():
    found, status, _forbidden = _route_ring(heated=True, keel_m=100e3)
    assert (found, status) == (None, "barrier")
    assert _route_ring(heated=False, keel_m=40e3)[1] == "barrier"


def test_heat_with_an_eroded_keel_opens_the_craton_and_the_route_enters_no_blocked_cell():
    found, status, forbidden = _route_ring(heated=True, keel_m=40e3)
    assert status == "separates"
    route, mask = found
    assert not forbidden[route].any()
    assert min(mask.sum(), (~mask).sum()) >= probe.MIN_DAUGHTER_SHARE * len(mask)


def test_router_picks_endpoints_that_leave_both_daughters_large_enough():
    size = 41
    fields = {"craton_crust_m": np.zeros(size * size), "plume_heat_myr": np.zeros(size * size),
              "mantle_lithosphere_thickness_m": np.full(size * size, 100e3), "belt_formed_years": np.zeros(size * size)}
    plate = _GridPlate(size, 0.5, fields)
    # A site 3 cells from the south edge, opening north-south: the nearest boundary either side
    # of the east-west strike is the south edge, which would cut a sliver.
    seed = int(plate.index[3, size // 2])
    none = np.zeros(plate.node_count(), dtype=bool)
    found, status, _cheapest = probe.route_by_daughters(plate, plate.points, seed, np.array([0.0, 0.0, 1.0]), none,
                                                        np.ones(plate.node_count()), plate.node_areas_m2(), np.radians(1.0))
    assert status == "separates"
    assert min(found[1].sum(), (~found[1]).sum()) >= probe.MIN_DAUGHTER_SHARE * plate.node_count()


def test_margin_share_ignores_an_arms_approach_but_counts_a_hook_along_the_edge():
    plate = _GridPlate(21, 0.5, {})
    edge = probe._edge_cells(plate)
    band = np.radians(1.0)
    straight = list(plate.index[:, 10])
    assert probe.margin_share(plate.points, straight, edge, band) == 0.0
    hook = list(plate.index[:12, 10]) + list(plate.index[11, 11:20]) + list(plate.index[11:, 20][1:])
    assert probe.margin_share(plate.points, hook, edge, band) > 0.0


def _routed(onset, setting, status, **extra):
    return _full_event(onset, setting, route_status=status, world_craton_km2=1e6, **extra)


def test_route_summary_excludes_barriers_and_applies_the_outcome_and_shape_limits():
    events = ([_routed(float(i), "passive", "separates", margin_share=0.0) for i in range(9)]
              + [_routed(50.0, "passive", "barrier"), _routed(60.0, "passive", "daughter_too_small")])
    row = probe.route_summary([_report(1, events)])["settings"]["passive"]
    assert (row["non_barrier"], row["breakup_share"], row["outcome_passes"]) == (10, 0.9, True)
    hooked = events[:8] + [_routed(9.0, "passive", "separates", margin_share=0.5)] * 2
    assert probe.route_summary([_report(1, hooked)])["settings"]["passive"]["shape_passes"] is False


def test_craton_budget_scales_with_the_run_and_counts_only_gated_breakups():
    # 2%/Gyr over 200 Myr of a 1e6 km2 craton is 4,000 km2.
    events = [_routed(10.0, "passive", "separates", craton_detached_km2=3000.0),
              _routed(80.0, "transtension", "daughter_too_small", craton_detached_km2=9000.0)]
    row = probe.craton_budget_row(_report(1, events))
    assert (row["budget_km2"], row["detached_km2"], row["passes"]) == (4000.0, 3000.0, True)

"""Candidate recorder and realism-target scorer for the rift-line split method (RFC 0018 D1, #888), reading the
production method's per-step records (`app.rift_damage.RiftRuntime`). Diagnostic only: nothing here is saved state.

Ported from the D1 prototype's `rift_damage_recorder.py` (PR #890) so tuning under RFC 0021 can score the method that
is in the model, not the prototype. The rules are the same.

``Recorder.observe(world, years)`` is called after every ``step_world``. It follows RFC 0018 step 5: a *physical
candidate* is a connected component of cells with ``rift_strain >= eps_cand``. Each step it samples every component's
area, length (extent along its principal strike), width (area / length), mean strain rate, and the share of its cells
on inherited weakness. Identity follows overlap with the previous step by area:

* a new component takes the identity it overlaps most;
* when one identity is claimed by several components, the largest keeps it and the rest start new identities with the
  parent's onset time; when a component overlaps several identities it keeps the oldest and the others are
  recorded as merged;
* a candidate is *active* while any of its cells strains, *failed* once none has for ``T_FAIL_YR``, and *successful*
  when a split path runs through it;
* candidates still active at the end of a run are censored, and candidates already present at the first observation
  have unknown onset: both are counted but kept out of the duration, onset and failure statistics.

``score(...)`` turns the pooled candidates and split events of one or more runs into the RFC's realism-target table;
``pool(directory, ...)`` does that for every world run of one arm.
"""

import json
import pickle
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components, dijkstra

# This checkout's backend, whatever PYTHONPATH says (docs/debugging.md).
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import rift_damage, rift_lines  # noqa: E402
from app.rift_strength import RADIUS_M  # noqa: E402

T_FAIL_YR = 5e6  # RFC 0018's t_fail: a candidate none of whose cells has strained for this long has failed
MIN_CELLS_TO_COUNT = 3  # a candidate counts once it has reached this many cells (predeclared in the #888 stage 1 post)
SUPER_SEA_WIDTH_M = 375e3
LOCALIZATION_COMPONENTS, LOCALIZATION_SHARE = 3, 0.30
SCORED = ("onset_width_km", "width_growth", "duration_myr", "transition_ratio", "inherited_weakness", "failed_rifts",
          "localization", "no_microplates", "line_shape", "splits_happen")
RANGES = {"onset_width_km": (100.0, 150.0), "duration_myr": (10.0, 50.0), "failed_rifts": (0.2, 0.8)}


@dataclass
class Candidate:
    ident: int
    plate_id: int
    onset_years: float
    onset_width_m: float
    from_start: bool
    samples: list = field(default_factory=list)  # (years, cells, area_m2, length_m, width_m, mean_rate_per_yr, weak_share)
    last_strain_years: float = 0.0
    max_cells: int = 0
    status: str = "active"  # active | failed | successful | merged
    end_years: float | None = None
    ever_failed: bool = False
    parent: int | None = None


def random_path_weak_share(step, plate_id: int, elapsed_years: float, trials: int = 20) -> float:
    """Mean inherited-weakness share along matched random edge-to-edge paths: shortest paths between random pairs of
    boundary cells at least 60% of the plate's longest boundary-to-boundary span apart (the falsification test's
    baseline for the inherited-weakness target). Deterministic in (plate, time)."""
    rng = np.random.default_rng((plate_id, round(elapsed_years)))
    cells = np.flatnonzero(step.boundary)
    if len(cells) < 4:
        return float("nan")
    span = 0.6 * float(np.max(np.linalg.norm(step.points[cells][:, None] - step.points[cells][None, :], axis=2)))
    shares = []
    for _ in range(trials * 5):
        a, b = rng.choice(cells, 2, replace=False)
        if np.linalg.norm(step.points[a] - step.points[b]) < span:
            continue
        _, pred = dijkstra(step.graph, directed=False, indices=a, return_predecessors=True)
        shares.append(float(step.weak[rift_lines._trace(pred, b)].mean()))
        if len(shares) == trials:
            break
    return float(np.mean(shares)) if shares else float("nan")


@dataclass
class Recorder:
    eps_cand: float
    candidates: dict = field(default_factory=dict)
    key_owner: dict = field(default_factory=dict)  # cell key -> identity at the previous observation
    next_ident: int = 0
    first_step: bool = True
    split_events_seen: int = 0
    super_sea_max_width_m: float = 0.0
    super_sea_violations: int = 0
    start_years: float = 0.0
    last_years: float = 0.0

    def _components(self, step):
        hit = step.eps >= self.eps_cand
        if not hit.any():
            return []
        cells = np.flatnonzero(hit)
        count, labels = connected_components(step.graph[hit][:, hit], directed=False)
        return [cells[labels == k] for k in range(count)]

    @staticmethod
    def _geometry(step, cells):
        pts = step.points[cells]
        area = float(step.areas[cells].sum())
        spacing = float(np.sqrt(np.median(step.areas[cells])))
        if len(cells) < 2:
            return area, spacing, area / spacing
        _, _, vt = np.linalg.svd(pts - pts.mean(axis=0), full_matrices=False)
        coord = (pts @ vt[0]) * RADIUS_M
        length = float(coord.max() - coord.min()) + spacing
        return area, length, area / length

    def _match(self, state) -> list:
        """Every strained component of this step, with its overlap (by area) with last step's identities."""
        pending = []
        for step in state.steps.values():
            for cells in self._components(step):
                keys = step.keys[cells]
                overlap: dict[int, float] = {}
                for key, area in zip(keys.tolist(), step.areas[cells].tolist()):
                    owner = self.key_owner.get(key)
                    if owner is not None:
                        overlap[owner] = overlap.get(owner, 0.0) + area
                pending.append((step, cells, keys, overlap))
        return pending

    def _identity(self, pending, index, claimed, now) -> int:
        step, cells, _, overlap = pending[index]
        if not overlap:
            _, _, width = self._geometry(step, cells)
            return self._new(now, width, step.plate_id, self.first_step, now)
        best = max(overlap, key=overlap.get)
        if index != max(claimed[best], key=lambda i: pending[i][3].get(best, 0.0)):
            parent = self.candidates[best]  # a piece split off a larger component: a new identity with the parent's onset
            return self._new(parent.onset_years, parent.onset_width_m, step.plate_id, parent.from_start, now, parent=best)
        ident = min(overlap, key=lambda i: self.candidates[i].onset_years)  # the oldest identity wins a merge
        for other in overlap:
            if other != ident and self.candidates[other].status == "active":
                self.candidates[other].status, self.candidates[other].end_years = "merged", now
        return ident

    def observe(self, world, years: float) -> None:
        state = rift_damage.runtime(world)
        now = world.elapsed_years
        if self.first_step:
            self.start_years = now
        self.last_years = now
        pending = self._match(state)
        claimed: dict[int, list] = {}
        for index, (_, _, _, overlap) in enumerate(pending):
            if overlap:
                claimed.setdefault(max(overlap, key=overlap.get), []).append(index)
        owners_now: dict[int, int] = {}
        for index, (step, cells, keys, _) in enumerate(pending):
            ident = self._identity(pending, index, claimed, now)
            self._sample(self.candidates[ident], step, cells, now, years)
            for key in keys.tolist():
                owners_now[key] = ident
        self.key_owner = owners_now
        self._close_out(state, now)
        self._check_super_seas(world, now)
        self.first_step = False

    def _new(self, onset_years, onset_width, plate_id, from_start, now, parent=None) -> int:
        ident = self.next_ident
        self.next_ident += 1
        self.candidates[ident] = Candidate(ident, plate_id, onset_years, onset_width, from_start, last_strain_years=now, parent=parent)
        return ident

    def _sample(self, cand: Candidate, step, cells, now, years) -> None:
        if cand.status in ("successful", "merged"):
            return  # the series ends at breakup (or at the merge)
        area, length, width = self._geometry(step, cells)
        rate = float(np.average(step.d_eps[cells], weights=step.areas[cells]) / max(years, 1.0))
        if (step.d_eps[cells] > 0).any():
            cand.last_strain_years = now
        cand.max_cells = max(cand.max_cells, len(cells))
        cand.plate_id = step.plate_id
        cand.samples.append((now, len(cells), area, length, width, rate, float(step.weak[cells].mean())))
        if cand.status == "failed":
            cand.status, cand.end_years = "active", None  # a scar reopened

    def _close_out(self, state, now) -> None:
        """Mark the candidates this step's splits ran through as successful, add the random-path weakness baseline to
        each new split event (read from the parent's record of this step's pass), and fail the candidates that have
        stopped straining."""
        for event in state.events[self.split_events_seen:]:
            parent = state.steps.get(event["plate_id"])
            event["random_weak_share"] = (random_path_weak_share(parent, event["plate_id"], event["elapsed_years"])
                                          if parent is not None else float("nan"))
            for ident in {self.key_owner.get(k) for k in event["path_keys"]} - {None}:
                cand = self.candidates[ident]
                if cand.status in ("active", "failed"):
                    cand.status, cand.end_years = "successful", event["elapsed_years"]
        self.split_events_seen = len(state.events)
        for cand in self.candidates.values():
            if cand.status == "active" and now - cand.last_strain_years >= T_FAIL_YR:
                cand.status, cand.end_years, cand.ever_failed = "failed", cand.last_strain_years, True

    def _check_super_seas(self, world, now) -> None:
        """Below-sea-level basins made by interior thinning wider than three onset widths across the narrow axis, outside
        crust thicker than the plate median (a proxy for thick, hot crust), checked every 5 Myr."""
        if round((now - self.start_years) / 1e6) % 5:
            return
        for plate in world.plates:
            wet = (plate.collect("rift_strain") >= self.eps_cand) & (plate.collect("elevation") < 0.0)
            if wet.sum() < 4:
                continue
            hc = plate.collect("crustal_thickness_m")
            wet &= hc <= np.median(hc)
            if wet.sum() < 4:
                continue
            points, _ = plate.all_points_and_elevation()
            count, labels = connected_components(rift_damage.plate_graph(plate, points)[wet][:, wet], directed=False)
            cells = np.flatnonzero(wet)
            for k in range(count):
                members = cells[labels == k]
                if len(members) < 4:
                    continue
                pts = points[members]
                _, _, vt = np.linalg.svd(pts - pts.mean(axis=0), full_matrices=False)
                narrow = float(np.ptp((pts @ vt[1]) * RADIUS_M))
                self.super_sea_max_width_m = max(self.super_sea_max_width_m, narrow)
                self.super_sea_violations += narrow > SUPER_SEA_WIDTH_M

    def finish(self, world) -> None:
        """Mark what is still active at the end as censored (status stays 'active')."""
        self.last_years = world.elapsed_years


# --- scoring -------------------------------------------------------------------------------------


def kaplan_meier_median(durations, observed) -> float | None:
    """Median of the survival curve for time-to-breakup, with failed/active candidates censored at their last
    observed time. None if the curve never reaches 0.5."""
    durations, observed = np.asarray(durations, float), np.asarray(observed, bool)
    if not len(durations):
        return None
    survival = 1.0
    for t in np.unique(durations[observed]):
        at_risk = (durations >= t).sum()
        survival *= 1.0 - (observed & (durations == t)).sum() / at_risk
        if survival <= 0.5:
            return float(t)
    return None


def _transition(c: Candidate) -> float | None:
    """Strain rate over the last quarter of a successful rift's life over that over its first quarter."""
    ordered = sorted(c.samples)
    times = np.array([s[0] for s in ordered])
    rates = np.array([s[5] for s in ordered])
    if len(times) < 4 or times[-1] <= times[0]:
        return None
    span = times[-1] - times[0]
    first = rates[times <= times[0] + span / 4].mean()
    last = rates[times >= times[-1] - span / 4].mean()
    return last / first if first > 0 else np.inf


def candidate_stats(cands: list[Candidate], last_years: float, drop_fragments: bool = True) -> dict:
    """Per-candidate quantities behind the targets, for candidates with known onset that reached the size filter.
    ``drop_fragments`` (the stage 1 correction, #888) leaves out split children seen at a single observation: pieces of
    an older cluster that were re-labelled when the cluster split, not rifts that started and stalled."""
    counted = [c for c in cands if not c.from_start and c.max_cells >= MIN_CELLS_TO_COUNT and c.samples
               and not (drop_fragments and c.parent is not None and len(c.samples) < 2)]
    by_status = {s: [c for c in counted if c.status == s] for s in ("successful", "failed", "active", "merged")}
    durations, observed = [], []
    for c in counted:
        if c.status == "merged":
            continue
        end = c.end_years if c.end_years is not None else last_years
        durations.append((end - c.onset_years) / 1e6)
        observed.append(c.status == "successful")
    transition = [r for r in map(_transition, by_status["successful"]) if r is not None]
    return {
        "counted": len(counted), "successful": len(by_status["successful"]), "failed": len(by_status["failed"]),
        "censored_active": len(by_status["active"]), "merged": len(by_status["merged"]),
        "from_start_excluded": sum(c.from_start for c in cands),
        "below_size_filter": sum((not c.from_start) and c.max_cells < MIN_CELLS_TO_COUNT for c in cands),
        "onset_width_km": [c.onset_width_m / 1e3 for c in counted if c.status != "merged"],
        "width_grew": [sorted(c.samples)[-1][4] > c.onset_width_m for c in by_status["successful"]],
        "duration_median_myr": kaplan_meier_median(durations, observed),
        "durations_completed_myr": [d for d, o in zip(durations, observed) if o],
        "transition_ratio": transition,
    }


def _rift_targets(stats: dict) -> dict:
    """The targets read from the candidates: onset width, width growth, duration, transition and failed share."""
    widths, grew, completed, ratios = (stats[k] for k in ("onset_width_km", "width_grew", "durations_completed_myr", "transition_ratio"))
    median_width = float(np.median(widths)) if widths else None
    # Scored on completed rifts (a correction to the predeclared survival median, #888): with most candidates failing the
    # survival curve never reaches 0.5, so the median is undefined whatever the breakups look like. Kept as an observation.
    median_completed = float(np.median(completed)) if completed else None
    completed_n = stats["successful"] + stats["failed"]
    fail_share = stats["failed"] / completed_n if completed_n else None
    return {
        "onset_width_km": {"value": median_width, "range": [100, 150], "ok": median_width is not None and 100 <= median_width <= 150},
        "width_growth": {"value": float(np.mean(grew)) if grew else None, "rule": "median successful band wider at breakup than onset",
                         "ok": bool(grew) and np.mean(grew) > 0.5},
        "duration_myr": {"value": median_completed, "survival_median": stats["duration_median_myr"], "range": [10, 50],
                         "ok": median_completed is not None and 10 <= median_completed <= 50},
        "transition_ratio": {"value": float(np.median(ratios)) if ratios else None, "minimum": 3.0, "ok": bool(ratios) and np.median(ratios) >= 3.0},
        "failed_rifts": {"value": fail_share, "range": [0.2, 0.8], "ok": fail_share is not None and 0.2 <= fail_share <= 0.8},
    }


def _split_targets(splits: list[dict], a_min: float) -> dict:
    """The targets read from the split events: inherited weakness, localization, daughter size and line shape."""
    out = {}
    if splits:
        share = np.mean([e["weak_share"] > 0.5 for e in splits])
        out["inherited_weakness"] = {"value": float(share), "minimum": 2 / 3, "mean_path_weak_share": float(np.mean([e["weak_share"] for e in splits])),
                                     "random_path_weak_share": float(np.nanmean([e["random_weak_share"] for e in splits])), "ok": share >= 2 / 3}
    else:
        out["inherited_weakness"] = {"value": None, "minimum": 2 / 3, "ok": False}
    bad = [e for e in splits if e["damaged_components"] > LOCALIZATION_COMPONENTS or e["damaged_share"] >= LOCALIZATION_SHARE]
    out["localization"] = {"splits": len(splits), "violations": len(bad), "max_damaged_share": max((e["damaged_share"] for e in splits), default=None),
                           "max_damaged_components": max((e["damaged_components"] for e in splits), default=None), "ok": not bad}
    small = [e["small_share"] for e in splits]
    out["no_microplates"] = {"splits": len(splits), "violations": int(sum(x < a_min - 1e-9 for x in small)), "min_small_share": min(small, default=None),
                             "median_small_share": float(np.median(small)) if small else None, "ok": bool(splits) and all(x >= a_min - 1e-9 for x in small)}
    ratios = [e["length_ratio"] for e in splits]
    out["line_shape"] = {"median_length_ratio": float(np.median(ratios)) if ratios else None, "minimum": 1.1, "ok": bool(ratios) and np.median(ratios) >= 1.1}
    out["plume_link"] = {"share_with_plume_on_path": float(np.mean([e["plume_share"] > 0 for e in splits])) if splits else None, "descriptive": True}
    out["oceanic_share_of_successful"] = {"value": float(np.mean([e["oceanic_share"] > 0.5 for e in splits])) if splits else None, "fallback_above": 1 / 3}
    return out


def score(candidates: list[Candidate], events: list[dict], last_years: float, super_sea_max_m: float, super_sea_violations: int,
          a_min: float, production_splits_per_world_100myr: float | None = None, world_myr: float | None = None,
          drop_fragments: bool = True) -> dict:
    """The RFC's realism-target table for the pooled runs. Each entry has its value, its pass range, and ``ok``.
    ``production_splits_per_world_100myr`` is the great-circle splitter's rate on the same worlds, against which the
    rift-line rate is judged at half or more; ``world_myr`` is the total simulated time over all pooled worlds in Myr."""
    stats = candidate_stats(candidates, last_years, drop_fragments)
    splits = [e for e in events if e["kind"] == "split"]
    out = {"candidates": {k: v for k, v in stats.items() if not isinstance(v, list)}, **_rift_targets(stats), **_split_targets(splits, a_min)}
    out["no_super_seas"] = {"max_narrow_width_km": super_sea_max_m / 1e3, "violations": super_sea_violations, "limit_km": SUPER_SEA_WIDTH_M / 1e3,
                            "ok": super_sea_violations == 0, "descriptive": True}
    rate = len(splits) / world_myr * 100.0 if world_myr else None
    out["splits_happen"] = {"rate_per_world_100myr": rate, "production_rate": production_splits_per_world_100myr,
                            "ok": bool(rate is not None and production_splits_per_world_100myr is not None and rate >= 0.5 * production_splits_per_world_100myr),
                            "note": "judged on the rift-line arm alone against one great-circle run per development world"}
    out["splits"] = {"value": len(splits)}
    out["targets_passed"] = int(sum(bool(out[k]["ok"]) for k in SCORED))
    out["targets_scored"] = len(SCORED)
    return out


# --- pooling the runs of one arm -----------------------------------------------------------------


def violation(name: str, entry: dict) -> float:
    """How far a target is from passing, as a fraction of its range (0 when met)."""
    if entry.get("ok"):
        return 0.0
    value = entry.get("value")
    if name in RANGES and value is not None:
        lo, hi = RANGES[name]
        return max(lo - value, value - hi, 0.0) / (hi - lo)
    floors = {"transition_ratio": ("value", 3.0), "inherited_weakness": ("value", 2 / 3), "line_shape": ("median_length_ratio", 1.1),
              "no_microplates": ("min_small_share", 0.15)}
    if name in floors and entry.get(floors[name][0]) is not None:
        key, floor = floors[name]
        return max(floor - entry[key], 0.0) / floor
    if name == "splits_happen" and entry.get("rate_per_world_100myr") is not None and entry.get("production_rate"):
        need = 0.5 * entry["production_rate"]
        return max(need - entry["rate_per_world_100myr"], 0.0) / max(need, 1e-9)
    return 1.0


def events_behind(name: str, result: dict) -> int:
    """How many pooled events a target's value rests on."""
    stats = result["candidates"]
    if name in ("width_growth", "transition_ratio", "duration_myr"):
        return stats["successful"]
    if name in ("onset_width_km", "failed_rifts"):
        return stats["counted"]
    return result["splits"]["value"] if name in ("inherited_weakness", "no_microplates", "line_shape", "localization") else 10**9


def pool(directory: Path, production_rate: float | None, a_min: float, min_events: int = 5, drop_fragments: bool = True) -> dict:
    """Score every world run of one arm together. Targets resting on fewer than ``min_events`` events count as not met."""
    cands, events, super_max, super_viol, myr, ledger, errors, last = [], [], 0.0, 0, 0.0, 0.0, [], 0.0
    for run in sorted(p for p in directory.iterdir() if p.is_dir()):
        pickled, summary = run / "recorder.pkl", run / "run.json"
        if not (pickled.exists() and summary.exists()):
            errors.append(f"{run.name}: no result")
            continue
        with open(pickled, "rb") as handle:
            data = pickle.load(handle)
        offset = len(cands)
        for c in data["candidates"]:
            c.ident += offset
            cands.append(c)
        events += data["events"]
        super_max, super_viol, myr = max(super_max, data["super_sea_max_m"]), super_viol + data["super_sea_violations"], myr + data["years_myr"]
        last = max(last, data["last_years"])
        ledger = max(ledger, max(json.loads(summary.read_text())["ledger_relative_gap_max"].values(), default=0.0))
    result = score(cands, events, last, super_max, super_viol, a_min, production_rate, myr, drop_fragments)
    for name in SCORED:
        if events_behind(name, result) < min_events:
            result[name]["ok"], result[name]["not_evaluable"] = False, True
    result["targets_passed"] = int(sum(bool(result[n]["ok"]) for n in SCORED))
    result["violation_total"] = float(sum(violation(n, result[n]) for n in SCORED))
    result["ledger_gap"], result["errors"] = ledger, errors
    return result

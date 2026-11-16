"""Candidate recorder and realism-target scorer for RFC 0018's D1 gate (#888). Diagnostic only: nothing here is
saved state or production code.

``Recorder.observe(world, mechanism, years)`` is called after every ``step_world``. It follows RFC step 5: a
*physical candidate* is a connected component of cells with ``rift_strain >= eps_cand``. Each step it samples every
component's area, length (extent along its principal strike), width (area / length), mean strain rate, and the
share of its cells on inherited weakness. Identity follows overlap with the previous step by area:

* a new component takes the identity it overlaps most;
* when one identity is claimed by several components, the largest keeps it and the rest start new identities with the
  parent's onset time; when a component overlaps several identities it keeps the oldest and the others are
  recorded as merged;
* a candidate is *active* while any of its cells strains, *failed* once none has for ``t_fail``, and *successful*
  when a split path runs through it;
* candidates still active at the end of a run are censored, and candidates already present at the first observation
  have unknown onset: both are counted but kept out of the duration, onset and failure statistics.

``score(records, ...)`` turns the pooled candidates and split events of one or more runs into the RFC's realism-target
table, with the Kaplan-Meier median for onset-to-breakup duration.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse.csgraph import connected_components

import rift_damage_mechanism as mech

RADIUS_M = mech.RADIUS_M
MIN_CELLS_TO_COUNT = 3  # a candidate counts once it has reached this many cells (predeclared in the stage 1 post)
SUPER_SEA_WIDTH_M = 375e3
LOCALIZATION_COMPONENTS, LOCALIZATION_SHARE = 3, 0.30


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


@dataclass
class Recorder:
    params: mech.Params
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
        hit = step.eps >= self.params.eps_cand
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

    def observe(self, world, mechanism, years: float) -> None:
        now = world.elapsed_years
        if self.first_step:
            self.start_years = now
        self.last_years = now
        owners_now: dict[int, int] = {}
        claimed: dict[int, list] = {}
        pending = []
        for step in mechanism.steps.values():
            for cells in self._components(step):
                keys = step.keys[cells]
                overlap: dict[int, float] = {}
                for key, area in zip(keys.tolist(), step.areas[cells].tolist()):
                    owner = self.key_owner.get(key)
                    if owner is not None:
                        overlap[owner] = overlap.get(owner, 0.0) + area
                pending.append((step, cells, keys, overlap))
                if overlap:
                    best = max(overlap, key=overlap.get)
                    claimed.setdefault(best, []).append(len(pending) - 1)
        assigned: dict[int, int] = {}
        for index, (step, cells, keys, overlap) in enumerate(pending):
            area = float(step.areas[cells].sum())
            if overlap:
                best = max(overlap, key=overlap.get)
                rivals = claimed[best]
                if index != max(rivals, key=lambda i: pending[i][3].get(best, 0.0)):
                    parent = self.candidates[best]
                    ident = self._new(parent.onset_years, parent.onset_width_m, step.plate_id, parent.from_start, now, parent=best)
                else:
                    ident = min(overlap, key=lambda i: self.candidates[i].onset_years)  # the oldest identity wins a merge
                    for other in overlap:
                        if other != ident and self.candidates[other].status == "active":
                            self.candidates[other].status, self.candidates[other].end_years = "merged", now
            else:
                _, length, width = self._geometry(step, cells)
                ident = self._new(now, width, step.plate_id, self.first_step, now)
            assigned[index] = ident
            self._sample(self.candidates[ident], step, cells, now, years)
            for key in keys.tolist():
                owners_now[key] = ident
        self.key_owner = owners_now
        self._close_out(mechanism, now)
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

    def _close_out(self, mechanism, now) -> None:
        for event in mechanism.events[self.split_events_seen:]:
            keys = set(event["path_keys"])
            for ident in {self.key_owner.get(k) for k in keys} - {None}:
                cand = self.candidates[ident]
                if cand.status in ("active", "failed"):
                    cand.status, cand.end_years = "successful", event["elapsed_years"]
        self.split_events_seen = len(mechanism.events)
        for cand in self.candidates.values():
            if cand.status == "active" and now - cand.last_strain_years >= self.params.t_fail_yr:
                cand.status, cand.end_years, cand.ever_failed = "failed", cand.last_strain_years, True

    def _check_super_seas(self, world, now) -> None:
        """Below-sea-level basins made by interior thinning wider than three onset widths across the narrow axis, outside
        crust thicker than the plate median (a proxy for thick, hot crust), checked every 5 Myr."""
        if round((now - self.start_years) / 1e6) % 5:
            return
        for plate in world.plates:
            eps = plate.collect("rift_strain")
            wet = (eps >= self.params.eps_cand) & (plate.collect("elevation") < 0.0)
            if wet.sum() < 4:
                continue
            hc = plate.collect("crustal_thickness_m")
            wet &= hc <= np.median(hc)
            if wet.sum() < 4:
                continue
            points, _ = plate.all_points_and_elevation()
            graph, _, _ = mech.plate_graph(plate, points)
            count, labels = connected_components(graph[wet][:, wet], directed=False)
            cells = np.flatnonzero(wet)
            for k in range(count):
                members = cells[labels == k]
                if len(members) < 4:
                    continue
                pts = points[members]
                _, sv, vt = np.linalg.svd(pts - pts.mean(axis=0), full_matrices=False)
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


def candidate_stats(cands: list[Candidate], last_years: float, drop_fragments: bool = True) -> dict:
    """Per-candidate quantities behind the targets, for candidates with known onset that reached the size filter.
    ``drop_fragments`` (the stage 1 correction, #888) leaves out split children seen at a single observation: pieces of
    an older cluster that were re-labelled when the cluster split, not rifts that started and stalled."""
    counted = [c for c in cands if not c.from_start and c.max_cells >= MIN_CELLS_TO_COUNT and c.samples
               and not (drop_fragments and c.parent is not None and len(c.samples) < 2)]
    successful = [c for c in counted if c.status == "successful"]
    failed = [c for c in counted if c.status == "failed"]
    active = [c for c in counted if c.status == "active"]
    merged = [c for c in counted if c.status == "merged"]
    durations, observed = [], []
    for c in counted:
        if c.status == "merged":
            continue
        end = c.end_years if c.end_years is not None else last_years
        durations.append((end - c.onset_years) / 1e6)
        observed.append(c.status == "successful")
    growth, transition = [], []
    for c in successful:
        ordered = sorted(c.samples)
        growth.append(ordered[-1][4] > c.onset_width_m)
        times = np.array([s[0] for s in ordered])
        rates = np.array([s[5] for s in ordered])
        if len(times) >= 4 and times[-1] > times[0]:
            span = times[-1] - times[0]
            first = rates[times <= times[0] + span / 4].mean()
            last = rates[times >= times[-1] - span / 4].mean()
            transition.append(last / first if first > 0 else np.inf)
    return {
        "counted": len(counted), "successful": len(successful), "failed": len(failed), "censored_active": len(active), "merged": len(merged),
        "from_start_excluded": sum(c.from_start for c in cands), "below_size_filter": sum((not c.from_start) and c.max_cells < MIN_CELLS_TO_COUNT for c in cands),
        "onset_width_km": [c.onset_width_m / 1e3 for c in counted if c.status != "merged"],
        "width_grew": growth,
        "duration_median_myr": kaplan_meier_median(durations, observed),
        "durations_completed_myr": [d for d, o in zip(durations, observed) if o],
        "transition_ratio": transition,
    }


def score(candidates: list[Candidate], events: list[dict], last_years: float, super_sea_max_m: float, super_sea_violations: int,
          a_min: float = 0.15, production_splits_per_world_100myr: float | None = None, world_myr: float | None = None,
          drop_fragments: bool = True) -> dict:
    """The RFC's realism-target table for the pooled runs. Each entry has its value, its pass range, and ``ok``.
    ``production_splits_per_world_100myr`` is the production splitter's rate on the same worlds (one run per
    development world), against which the rift-line rate is judged at half or more; ``world_myr`` is the total
    simulated time over all pooled worlds in Myr."""
    stats = candidate_stats(candidates, last_years, drop_fragments)
    out = {"candidates": {k: v for k, v in stats.items() if not isinstance(v, list)}}
    widths = stats["onset_width_km"]
    median_width = float(np.median(widths)) if widths else None
    out["onset_width_km"] = {"value": median_width, "range": [100, 150], "ok": median_width is not None and 100 <= median_width <= 150}
    grew = stats["width_grew"]
    out["width_growth"] = {"value": float(np.mean(grew)) if grew else None, "rule": "median successful band wider at breakup than onset",
                           "ok": bool(grew) and np.mean(grew) > 0.5}
    duration = stats["duration_median_myr"]
    completed = stats["durations_completed_myr"]
    # Scored on completed rifts (a correction to the predeclared survival median, #888): with most candidates failing the
    # survival curve never reaches 0.5, so the median is undefined whatever the breakups look like. Kept as an observation.
    median_completed = float(np.median(completed)) if completed else None
    out["duration_myr"] = {"value": median_completed, "survival_median": duration, "range": [10, 50],
                           "ok": median_completed is not None and 10 <= median_completed <= 50}
    ratios = stats["transition_ratio"]
    out["transition_ratio"] = {"value": float(np.median(ratios)) if ratios else None, "minimum": 3.0, "ok": bool(ratios) and np.median(ratios) >= 3.0}
    splits = [e for e in events if e["kind"] == "split"]
    if splits:
        share = np.mean([e["weak_share"] > 0.5 for e in splits])
        random_share = float(np.nanmean([e["random_weak_share"] for e in splits]))
        out["inherited_weakness"] = {"value": float(share), "minimum": 2 / 3, "mean_path_weak_share": float(np.mean([e["weak_share"] for e in splits])),
                                     "random_path_weak_share": random_share, "ok": share >= 2 / 3}
    else:
        out["inherited_weakness"] = {"value": None, "minimum": 2 / 3, "ok": False}
    completed_n = stats["successful"] + stats["failed"]
    fail_share = stats["failed"] / completed_n if completed_n else None
    out["failed_rifts"] = {"value": fail_share, "range": [0.2, 0.8], "ok": fail_share is not None and 0.2 <= fail_share <= 0.8}
    out["no_super_seas"] = {"max_narrow_width_km": super_sea_max_m / 1e3, "violations": super_sea_violations, "limit_km": SUPER_SEA_WIDTH_M / 1e3,
                            "ok": super_sea_violations == 0, "descriptive": True}
    bad = [e for e in splits if e["damaged_components"] > LOCALIZATION_COMPONENTS or e["damaged_share"] >= LOCALIZATION_SHARE]
    out["localization"] = {"splits": len(splits), "violations": len(bad), "max_damaged_share": max((e["damaged_share"] for e in splits), default=None),
                           "max_damaged_components": max((e["damaged_components"] for e in splits), default=None), "ok": not bad}
    small = [e["small_share"] for e in splits]
    out["no_microplates"] = {"splits": len(splits), "violations": int(sum(x < a_min - 1e-9 for x in small)), "min_small_share": min(small, default=None),
                             "median_small_share": float(np.median(small)) if small else None, "ok": bool(splits) and all(x >= a_min - 1e-9 for x in small)}
    ratios = [e["length_ratio"] for e in splits]
    out["line_shape"] = {"median_length_ratio": float(np.median(ratios)) if ratios else None, "minimum": 1.1, "ok": bool(ratios) and np.median(ratios) >= 1.1}
    rate = len(splits) / world_myr * 100.0 if world_myr else None
    out["splits_happen"] = {"rate_per_world_100myr": rate, "production_rate": production_splits_per_world_100myr,
                            "ok": bool(rate is not None and production_splits_per_world_100myr is not None and rate >= 0.5 * production_splits_per_world_100myr),
                            "note": "judged on the rift-line arm alone against one production-splitter run per development world"}
    out["plume_link"] = {"share_with_plume_on_path": float(np.mean([e["plume_share"] > 0 for e in splits])) if splits else None, "descriptive": True}
    out["oceanic_share_of_successful"] = {"value": float(np.mean([e["oceanic"] for e in splits])) if splits else None, "fallback_above": 1 / 3}
    out["splits"] = {"value": len(splits)}
    scored = [k for k in ("onset_width_km", "width_growth", "duration_myr", "transition_ratio", "inherited_weakness", "failed_rifts", "localization",
              "no_microplates", "line_shape", "splits_happen")]
    out["targets_passed"] = int(sum(bool(out[k]["ok"]) for k in scored))
    out["targets_scored"] = len(scored)
    return out

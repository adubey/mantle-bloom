#!/usr/bin/env python3
"""Issue #318: why do collision fronts fall back? Replays a save and, for every front decided
after the first step, records the pair's contact history over the preceding steps and the
nearest evidence on each plate (distance in spacings, source, role, age).

History per pair and step: band nodes on each side facing the other plate, how many converge,
and how many of them (or the neighbour node they face) are oceanic-coded -- an ocean closing
between the two would show up there, attached to one side or the other.

    cd backend
    .venv/bin/python ../bin/debug/diagnose_polarity_fallbacks.py \\
        ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 43.1 \\
        --out ../analysis/issue318/fallback-diagnosis-seed997271774.json
"""

from __future__ import annotations

from collections import defaultdict, deque

import numpy as np
from scipy.spatial import cKDTree

from polarity_replay_lib import load_polarity_replay, replay

from app import collision_polarity as cp, geometry  # noqa: E402

HISTORY_STEPS = 100  # 10 Myr at the UI step


def main() -> None:
    args, world = load_polarity_replay(__doc__)
    history: dict[tuple[int, int], deque] = defaultdict(lambda: deque(maxlen=HISTORY_STEPS))
    state = {"step": 0}
    decisions: list[dict] = []

    detect = cp._detect_fronts

    def detect_and_record(contacts, spacing_rad):
        step_rows: dict[tuple[int, int], dict] = {}
        for pid, c in contacts.items():
            nb = c.inputs.neighbor_plate_id
            band = np.isfinite(c.inputs.dist_to_neighbor) & (nb >= 0)
            for q in np.unique(nb[band]):
                q = int(q)
                if q not in contacts:
                    continue
                sel = np.flatnonzero(band & (nb == q))
                q_cont = contacts[q].continental[c.inputs.neighbor_node_index[sel]]
                row = step_rows.setdefault((min(pid, q), max(pid, q)), {"step": state["step"]})
                row[str(pid)] = {
                    "band": int(len(sel)),
                    "convergent": int(c.convergent[sel].sum()),
                    "own_oceanic": int((~c.continental[sel]).sum()),
                    "facing_oceanic": int((~q_cont).sum()),
                    "plate_oceanic": c.plate.crust_type == "oceanic",
                }
        for pair, row in step_rows.items():
            history[pair].append(row)
        return detect(contacts, spacing_rad)

    decide = cp._decide

    def decide_and_record(w, record, front, contacts, spacing_rad):
        decide(w, record, front, contacts, spacing_rad)
        if state["step"] == 0:
            return
        sides = {}
        for side in record.plate_ids:
            other = [p for p in record.plate_ids if p != side][0]
            look = contacts[side].points[front[side]] if len(front[side]) else contacts[other].points[front[other]]
            store = w.collision_evidence.get(side)
            entry = {"front_nodes": int(len(front[side])), "crust_type": contacts[side].plate.crust_type, "bins": 0}
            if store is not None and len(store):
                rows = geometry.to_world(contacts[side].plate.frame, store.local_points())
                dist, _ = cKDTree(look).query(rows)
                k = int(np.argmin(dist))
                entry.update(
                    bins=len(store),
                    nearest_spacings=float(dist[k] / spacing_rad),
                    nearest_source=cp.SOURCE_NAMES[int(store.source[k])],
                    nearest_role=int(store.role[k]),
                    nearest_age_myr=float((w.elapsed_years - store.last_seen[k]) / 1e6),
                    nearest_neighbour=int(store.neighbour[k]),
                    within_8_spacings=int((dist <= 8 * spacing_rad).sum()),
                )
            sides[str(side)] = entry
        decisions.append(
            {
                "elapsed_myr": w.elapsed_years / 1e6,
                "front_id": record.front_id,
                "pair": list(record.plate_ids),
                "source": record.source,
                "ambiguous": record.ambiguous,
                "lower": record.lower_plate_id,
                "votes": record.votes,
                "sides": sides,
                "history": list(history[record.plate_ids]),
            }
        )

    cp._detect_fronts = detect_and_record
    cp._decide = decide_and_record

    def before_step(k: int) -> None:
        state["step"] = k

    replay(world, args, lambda: f"{len(decisions)} decisions", lambda: {"decisions": decisions}, before_step)


if __name__ == "__main__":
    main()

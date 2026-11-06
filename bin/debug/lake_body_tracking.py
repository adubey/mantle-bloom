"""Lake-body tracking shared by the lake-flicker replay drivers (#588, #754, #811).

A body is a connected component of wet nodes. Bodies are matched between steps by node overlap, and
a large-body flicker reversal is #754's definition: the body is over 10,000 km^2 at both steps, and
two consecutive volume changes have opposite signs and each exceeds half its earlier volume.
"""

from __future__ import annotations

import numpy as np

LARGE_BODY_AREA_KM2 = 10_000.0
REVERSAL_FRACTION = 0.5


def components(fields) -> list[dict]:
    from app.hydroclimate import hydrology, lakes  # noqa: E402

    wet = fields.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    members = hydrology.lake_components(wet, fields.neighbor_idx)
    hypsometry = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    area = hypsometry.area_m2
    records = []
    for nodes in members:
        node_points = np.asarray(fields.points)[nodes]
        centre = node_points.mean(axis=0)
        records.append({
            "nodes": nodes,
            # The wet node nearest the normalized centroid, #680's anchor rule.
            "anchor": node_points[int(np.argmax(node_points @ (centre / np.linalg.norm(centre))))].tolist(),
            "area_km2": float(area[nodes].sum() / 1e6),
            "volume_m3": float(np.dot(fields.lake_depth[nodes], area[nodes])),
        })
    return records


def match_bodies(previous: list[dict], current: list[dict]) -> dict[int, int]:
    candidates = []
    for old_i, old in enumerate(previous):
        for new_i, new in enumerate(current):
            overlap = np.intersect1d(old["nodes"], new["nodes"], assume_unique=True).size
            if overlap:
                candidates.append((overlap, old_i, new_i))
    candidates.sort(reverse=True)
    old_for_new: dict[int, int] = {}
    used_old: set[int] = set()
    for overlap, old_i, new_i in candidates:
        if old_i in used_old or new_i in old_for_new:
            continue
        smaller = min(len(previous[old_i]["nodes"]), len(current[new_i]["nodes"]))
        if overlap >= max(1, smaller // 2):
            old_for_new[new_i] = old_i
            used_old.add(old_i)
    return old_for_new


def count_reversals(previous: list[dict], current: list[dict]) -> tuple[int, int, list[dict]]:
    """Return (reversals, comparable large bodies, events); stores each body's `delta_m3` on `current`."""
    reversals = eligible = 0
    events = []
    for new_i, old_i in match_bodies(previous, current).items():
        old, new = previous[old_i], current[new_i]
        if old["area_km2"] <= LARGE_BODY_AREA_KM2 or new["area_km2"] <= LARGE_BODY_AREA_KM2:
            continue
        eligible += 1
        prior_delta = old.get("delta_m3")
        delta = new["volume_m3"] - old["volume_m3"]
        if prior_delta is not None and prior_delta * delta < 0.0:
            baseline = max(old["volume_m3"], 1.0)
            if abs(prior_delta) > REVERSAL_FRACTION * baseline and abs(delta) > REVERSAL_FRACTION * baseline:
                reversals += 1
                # RFC 0016's window scan needs to know which lake reversed, not just how many.
                events.append({
                    "prior_delta_km3": prior_delta / 1e9,
                    "delta_km3": delta / 1e9,
                    "volume_before_km3": old["volume_m3"] / 1e9,
                    "volume_after_km3": new["volume_m3"] / 1e9,
                    "area_before_km2": old["area_km2"],
                    "area_after_km2": new["area_km2"],
                    "anchor": new["anchor"],
                })
        new["delta_m3"] = delta
    return reversals, eligible, events

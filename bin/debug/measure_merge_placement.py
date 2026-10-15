#!/usr/bin/env python3
"""Replay a captured collision merge and attribute its craton loss (RFC 0001 gate, #440).

`probe_craton_clips.py --capture-dir` saves the world just before each merge with RFC 0001's
signature. This loads one, runs the real `merge_split._merge_plates_with_craton_audit` on it
with read-only taps, checks that the replay books the same craton as the run did, and
reports what decides between RFC 0001's hypotheses:

- **room (H1):** the Hc stacked past the suture cap (`overflow`), the eligible cells
  `quad_merge._place_stacking_overflow` hands `_place_suture_crust`, their area and the room
  under the cap among them, all of it and in the continental components the front touches;
- **mask (H2):** the eligible area against the two plates' continental area before the merge,
  and thick crust (Hc >= the craton host floor) left out of the mask;
- **placement:** the volume each stage of `_place_suture_crust` placed, and the share of the
  overflow placed;
- **carry (H4):** the craton moving with the overflow, and the craton still above its column
  or on a non-continental cell right after placement;
- **later writers (H3):** that excess after the whole merge, split by where it sits (overflow
  front, a cell the placement changed, elsewhere; above Hc or off-continent), against what
  the audit books.

With `--contacts`, the pair's RFC 0008 contact rows over `--history-myr` before the merge
(closing rate, jam, brake, overlap fraction and craton-on-craton overlap volume).

    backend/.venv/bin/python bin/debug/measure_merge_placement.py /tmp/rfc0001/seed7-*.pkl \\
        --contacts /tmp/rfc0001/contacts-7.jsonl --json /tmp/rfc0001/placement-7.json
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import continental_ledger, cratons, merge_split, persistence, quad_merge, quad_tectonics  # noqa: E402
from app.elevation_lines import effective_is_continental_from_codes  # noqa: E402
from app.lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M  # noqa: E402

MKM3 = 1e15
KM2 = 1e6


def _continental(codes: np.ndarray, plate) -> np.ndarray:
    return effective_is_continental_from_codes(codes, plate.crust_type == "continental")


def _excess(craton: np.ndarray, hc: np.ndarray, continental: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(craton above a continental column, craton on a non-continental cell), per cell."""
    above = np.where(continental, np.maximum(craton - np.maximum(hc, 0.0), 0.0), 0.0)
    off = np.where(continental, 0.0, np.maximum(craton, 0.0))
    return above, off


def _room(keep, hc, eligible, front, areas) -> dict:
    """Room under the cap among eligible cells, and in the eligible components the front touches."""
    room = np.where(eligible, np.maximum(SUTURE_ACCRETION_MAX_HC_M - hc, 0.0), 0.0) * areas
    adjacency = quad_tectonics._adjacency_matrix(keep)
    at = np.flatnonzero(eligible)
    _, labels = connected_components(adjacency[at][:, at], directed=False)
    component = np.full(len(hc), -1)
    component[at] = labels
    touched = np.unique(component[front][component[front] >= 0])
    connected = np.isin(component, touched) & eligible
    return {
        "eligible_room_mkm3": float(room.sum()) / MKM3,
        "front_component_room_mkm3": float(room[connected].sum()) / MKM3,
        "eligible_components": int(labels.max() + 1) if len(at) else 0,
        "front_components": int(len(touched)),
        "front_cells_off_mask": int(np.count_nonzero(component[front] < 0)),
    }


class Taps:
    """Read-only taps on one merge's placement; `record` collects what they see."""

    def __init__(self) -> None:
        self.record: dict = {}
        self.stages: dict | None = None
        self.state: dict = {}

    def install(self) -> None:
        place_overflow = quad_merge._place_stacking_overflow
        place_suture = quad_tectonics._place_suture_crust
        merge_plates = merge_split.merge_plates

        def suture(*args, **kwargs):
            changed, stages = place_suture(*args, **kwargs)
            self.stages = dict(stages)
            return changed, stages

        def overflow_tap(keep, out, overflow, areas):
            self._before_placement(keep, out, overflow, areas)
            changed, placed = place_overflow(keep, out, overflow, areas)
            self._after_placement(keep, out, overflow, areas, changed, placed)
            return changed, placed

        def merge(world, id_keep, id_absorb):
            merge_plates(world, id_keep, id_absorb)
            self._after_merge(next(p for p in world.plates if p.plate_id == id_keep))

        quad_tectonics._place_suture_crust = suture
        quad_merge._place_stacking_overflow = overflow_tap
        merge_split.merge_plates = merge

    def _before_placement(self, keep, out, overflow, areas) -> None:
        hc = out["crustal_thickness_m"]
        continental = _continental(out["crust_type_code"], keep)
        eligible = continental if np.any(continental) else np.ones(len(hc), dtype=bool)
        front = np.flatnonzero(overflow > 0.0)
        craton = out.get("craton_crust_m", np.zeros(len(hc)))
        thick = hc >= cratons.CRATON_HOST_MIN_HC_M
        new = ~np.isin(keep.cell_keys, self.state["keys_before"])
        self.record["placement"] = {
            "overflow_stacked_mkm3": float(overflow[~new] @ areas[~new]) / MKM3,
            "overflow_new_mkm3": float(overflow[new] @ areas[new]) / MKM3,
            "front_stacked_area_km2": float(areas[(overflow > 0.0) & ~new].sum()) / KM2,
            "front_new_area_km2": float(areas[(overflow > 0.0) & new].sum()) / KM2,
            "new_area_km2": float(areas[new].sum()) / KM2,
            "new_hc_capped_mean_km": float(hc[new] @ areas[new]) / max(float(areas[new].sum()), 1.0) / 1e3,
            "new_hc_uncapped_mean_km": float((hc + overflow)[new] @ areas[new]) / max(float(areas[new].sum()), 1.0) / 1e3,
            "overflow_mkm3": float(overflow @ areas) / MKM3,
            "front_cells": int(len(front)),
            "front_area_km2": float(areas[front].sum()) / KM2,
            "merged_cells": int(len(hc)),
            "merged_area_km2": float(areas.sum()) / KM2,
            "eligible_cells": int(eligible.sum()),
            "eligible_area_km2": float(areas[eligible].sum()) / KM2,
            "thick_off_mask_area_km2": float(areas[thick & ~eligible].sum()) / KM2,
            "thick_off_mask_room_mkm3": float(np.maximum(SUTURE_ACCRETION_MAX_HC_M - hc, 0.0)[thick & ~eligible] @ areas[thick & ~eligible]) / MKM3,
            "craton_moving_mkm3": float(np.maximum(craton - hc, 0.0) @ areas) / MKM3,
            "craton_moving_new_mkm3": float(np.maximum(craton - hc, 0.0)[new] @ areas[new]) / MKM3,
            # Craton above even the uncapped column: the remap gave a cell more craton than
            # crust (each extensive field is restored by its own ratio), so no share of the
            # overflow can carry it.
            "craton_above_uncapped_mkm3": float(np.maximum(craton - hc - overflow, 0.0) @ areas) / MKM3,
            "craton_on_front_mkm3": float(craton[front] @ areas[front]) / MKM3,
            **_room(keep, hc, eligible, front, areas),
        }
        self.state.update(front=front, eligible=eligible)

    def _after_placement(self, keep, out, overflow, areas, changed, placed) -> None:
        above, off = _excess(out.get("craton_crust_m", np.zeros(len(areas))), out["crustal_thickness_m"],
                             _continental(out["crust_type_code"], keep))
        self.record["placement"].update({
            "placed_share": float(placed),
            "changed_cells": int(changed.sum()),
            "stages_mkm3": {k: v / MKM3 for k, v in (self.stages or {}).items() if v},
            "craton_above_after_placement_mkm3": float(above @ areas) / MKM3,
            "craton_off_continent_after_placement_mkm3": float(off @ areas) / MKM3,
        })
        self.state.update(changed=changed.copy(), keys=keep.cell_keys.copy())

    def before_merge(self, keep, absorb, others) -> None:
        """Where the absorbed plate's area goes: onto `keep` (stacked), onto a third plate's
        territory (its cells are never activated, `quad_merge._on_other_plates`), or open
        lattice; and the Hc volume on each, which the remap conserves."""
        self.state = {"keys_before": keep.cell_keys.copy()}
        source, local, area = quad_merge.subsample_cells(absorb, quad_merge.MERGE_REMAP_SUBSAMPLES)
        world_xyz = quad_merge.geometry.to_world(absorb.frame, local)
        stacked = keep._leaf_key_at(quad_merge.geometry.to_local(keep.frame, world_xyz)) >= 0
        on_other = np.zeros(len(area), dtype=bool)
        for other in others:
            on_other[~stacked] |= other.contains_batch(world_xyz[~stacked])
        hc = absorb.collect("crustal_thickness_m")[source] * area
        open_ = ~stacked & ~on_other
        self.record["absorbed"] = {
            "stacked_area_km2": float(area[stacked].sum()) / KM2,
            "on_third_plates_area_km2": float(area[on_other].sum()) / KM2,
            "open_area_km2": float(area[open_].sum()) / KM2,
            "stacked_hc_mkm3": float(hc[stacked].sum()) / MKM3,
            "on_third_plates_hc_mkm3": float(hc[on_other].sum()) / MKM3,
            "open_hc_mkm3": float(hc[open_].sum()) / MKM3,
        }

    def _after_merge(self, keep) -> None:
        """Where the craton the audit will clip sits, on the merged plate before the audit."""
        areas = keep.node_areas_m2()
        above, off = _excess(keep.collect("craton_crust_m"), keep.collect("crustal_thickness_m"),
                             _continental(keep.collect("crust_type_code"), keep))
        excess = above + off
        front = np.zeros(len(areas), dtype=bool)
        changed = np.zeros(len(areas), dtype=bool)
        if "keys" in self.state and np.array_equal(self.state["keys"], keep.cell_keys):
            front[self.state["front"]] = True
            changed = self.state["changed"] & ~front
        elsewhere = ~front & ~changed
        self.record["after_merge"] = {
            "excess_mkm3": float(excess @ areas) / MKM3,
            "above_hc_mkm3": float(above @ areas) / MKM3,
            "off_continent_mkm3": float(off @ areas) / MKM3,
            "on_front_mkm3": float(excess[front] @ areas[front]) / MKM3,
            "on_changed_mkm3": float(excess[changed] @ areas[changed]) / MKM3,
            "elsewhere_mkm3": float(excess[elsewhere] @ areas[elsewhere]) / MKM3,
            "cells_with_excess": int(np.count_nonzero(excess > 1.0)),
        }


def _pre_merge(world, keep, absorb) -> dict:
    out = {}
    for name, plate in (("keep", keep), ("absorb", absorb)):
        areas = plate.node_areas_m2()
        continental = _continental(plate.collect("crust_type_code"), plate)
        hc = plate.collect("crustal_thickness_m")
        above, off = _excess(plate.collect("craton_crust_m"), hc, continental)
        out[name] = {
            "cells": plate.node_count(),
            "area_km2": float(areas.sum()) / KM2,
            "continental_area_km2": float(areas[continental].sum()) / KM2,
            "hc_volume_mkm3": float(hc @ areas) / MKM3,
            "continental_hc_mean_km": float(hc[continental] @ areas[continental]) / max(float(areas[continental].sum()), 1.0) / 1e3,
            "area_at_70km_plus_km2": float(areas[hc >= 70_000.0].sum()) / KM2,
            "craton_mkm3": cratons.plate_volume_m3(world, plate) / MKM3,
            "craton_outside_column_mkm3": float((above + off) @ areas) / MKM3,
        }
    return out


def _history(path: Path | None, myr: float, pair: tuple[int, int], span: float) -> list[dict]:
    if path is None:
        return []
    key = f"{min(pair)}-{max(pair)}"
    rows = []
    for line in path.read_text().splitlines():
        row = json.loads(line)
        at = row.get("myr")
        if at is not None and myr - span - 1e-6 <= at <= myr + 1e-6 and key in row.get("pairs", {}):
            rows.append({"myr": at, **row["pairs"][key]})
    return rows


def measure(path: Path, contacts: Path | None, span: float) -> dict:
    capture = pickle.loads(path.read_bytes())
    world = persistence.load_world_bytes(capture["world"])
    cratons.ensure_ledger(world)
    continental_ledger.ensure_initialized(world)
    keep = next(p for p in world.plates if p.plate_id == capture["keep"])
    absorb = next(p for p in world.plates if p.plate_id == capture["absorb"])
    result = {k: v for k, v in capture.items() if k != "world"}
    result["pre_merge"] = _pre_merge(world, keep, absorb)
    taps = Taps()
    taps.install()
    taps.before_merge(keep, absorb, [p for p in world.plates if p is not keep and p is not absorb])
    sinks = dict(world.craton_ledger)
    material = world.continental_material_ledger.get("collision_subducted_m3", 0.0)
    started = time.time()
    merge_split._merge_plates_with_craton_audit(world, capture["keep"], capture["absorb"])
    result["replay_wall_s"] = time.time() - started
    booked = {k: (world.craton_ledger[k] - sinks.get(k, 0.0)) / MKM3 for k in cratons.CRATON_SINK_ACCOUNTS
              if abs(world.craton_ledger[k] - sinks.get(k, 0.0)) > 1e9}
    result.update(taps.record)
    result["booked_mkm3"] = booked
    result["material_clipped_mkm3"] = (world.continental_material_ledger["collision_subducted_m3"] - material) / MKM3
    result["reproduced"] = bool(np.isclose(sum(booked.values()) * MKM3, capture["booked_m3"], rtol=1e-6))
    result["history"] = _history(contacts, capture["myr"], (capture["keep"], capture["absorb"]), span)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("captures", type=Path, nargs="+")
    parser.add_argument("--contacts", type=Path, help="probe_craton_clips.py --contacts output for the same run")
    parser.add_argument("--history-myr", type=float, default=30.0)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    results = []
    for path in args.captures:
        # One process per capture keeps the taps from stacking; run them in turn.
        result = measure(path, args.contacts, args.history_myr)
        print(json.dumps({k: v for k, v in result.items() if k != "history"}, indent=1), flush=True)
        results.append(result)
        if len(args.captures) > 1:
            _untap()
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(results, indent=1) + "\n")


_ORIGINALS = (quad_merge._place_stacking_overflow, quad_tectonics._place_suture_crust, merge_split.merge_plates)


def _untap() -> None:
    quad_merge._place_stacking_overflow, quad_tectonics._place_suture_crust, merge_split.merge_plates = _ORIGINALS


if __name__ == "__main__":
    main()

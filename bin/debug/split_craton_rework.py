#!/usr/bin/env python3
"""Replay a saved world and split the craton ledger's sinks by the site that booked them
(issue #332).

`collision_reworked_m3` is booked at several sites; this taps `cratons.record` and keys
each booking by its caller, so the share from each mechanism can be read off:

- **polarized transfer** (`crust_transfer._book`, #320), split into crust kept on the upper
  plate (underthrust and thrust-up placement) and crust handed to the front's other
  overriders, in proportion to the Hc each took.
- **own-plate suture accretion** (`quad_tectonics` suture retreat, #290), split into crust
  placed on the consumed plate's own survivors and crust handed to its overrider, by the
  path's `placed_share` / `handed_share`, and by the donor plate's nominal crust type.
- **anatexis relamination** and **orogen relaxation** (`orogeny`).
- **phase-audit clips and residuals** (`cratons.PhaseAudit.settle`, the faults phase).

Every other sink account is split by caller too, so a phase audit's clips and residuals
(`cratons.PhaseAudit.settle`) read apart from the sink's own site (e.g. real rift thinning
in `rifted_m3`, #340 item 6). Checkpoints
also record the craton inventory, the cratonic area (cells with `craton_crust_m > 0`) and
the mean relief of the non-cratonic continental crust within CRATON_EDGE_RELIEF_WINDOW_KM
of a craton (every craton edge, not only those facing convergent fronts). For #340's
side effects they also record land and emergent continental area, land hypsometry, the share
of cratonic area under orogenic crust (Hc above the craton host window), the orogenic
`ceiling_overflow_m3` (ponding against cratons sends more crust to the Hc cap), the
area-weighted mean plate surface speed (collision braking moves it), and the plate merges
(collision or forced fusion) so far.

Given `--baseline` (this script's JSON from `main`, same save and span), it checks issue
#340's acceptance targets, both uncalibrated (tracked in #313):

- **Non-attrition loss rate.** The craton volume booked to every sink account
  (`cratons.CRATON_SINK_ACCOUNTS`: rifting, subduction, delamination, collisional reworking,
  erosion, topology removal, audit clips and unattributed loss) over the replay, less what
  trench attrition booked (the `cratons.attrite` site of `subducted_m3`), divided by the
  replay's span in Myr. Sources (`formed_m3`) are not netted off. It must be at most
  NON_ATTRITION_LOSS_MAX_FRACTION of the baseline's rate on the same save and span.
- **Craton-edge relief** above the baseline's at the last checkpoint.

`--rescore RUN_JSON --baseline MAIN_JSON` recomputes the targets from finished replays
without stepping anything.

Runs on any branch: sites a branch doesn't have simply don't appear.

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/split_craton_rework.py \\
        ~/Downloads/mantle-bloom-seed910211954-220200000y.mbworld --myr 20.2 --out split.json \\
        [--baseline split-main.json]
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/split_craton_rework.py \
        --rescore split.json --baseline split-main.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import continental_ledger, cratons, merge_split, persistence, quad_tectonics, shortening, world as world_mod  # noqa: E402
from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad  # noqa: E402
from app.mantle import PLANET_RADIUS_KM  # noqa: E402

KM2 = 1e6
KM3 = 1e9
SPLIT_ACCOUNTS = cratons.CRATON_SINK_ACCOUNTS
# Issue #340's acceptance targets, both uncalibrated (#313). The non-attrition craton loss
# rate may be at most this fraction of `main`'s on the same save and span.
NON_ATTRITION_LOSS_MAX_FRACTION = 0.2
# Craton-edge relief is the mean relief of non-cratonic continental crust within this
# distance of a craton: roughly the width of a foreland thrust belt against a cratonic
# backstop.
CRATON_EDGE_RELIEF_WINDOW_KM = 300.0
# Margin attrition is the one craton sink #340 allows; it books to `subducted_m3` from this
# site, alongside no-outlet suture subduction, which still counts as loss.
ATTRITION_SITE = ("subducted_m3", "cratons.attrite")
MERGES = [0]


def _count_merges() -> None:
    original = merge_split._merge_plates_with_craton_audit

    def counted(*args, **kwargs):
        MERGES[0] += 1
        return original(*args, **kwargs)

    merge_split._merge_plates_with_craton_audit = counted


def _site(frame) -> str:
    module = Path(frame.f_code.co_filename).stem
    return f"{module}.{frame.f_code.co_name}"


def _split(account: str, volume: float, frame) -> dict[str, float]:
    """`volume` booked from `frame`, keyed by sub-source where the caller's locals say."""
    site = _site(frame)
    local = frame.f_locals
    if site == "crust_transfer._book" and account == "collision_reworked_m3":
        result = local["result"]
        upper = result.underthrust_placed_m3 + result.upper_placed_m3
        total = upper + result.handed_m3
        if total <= 0.0:
            return {f"{site}:upper": volume}
        return {f"{site}:upper": volume * upper / total, f"{site}:handed": volume * result.handed_m3 / total}
    if site.startswith("quad_tectonics.") and "lost_share" in local:
        kind = getattr(local.get("plate"), "crust_type", "?")
        site = f"{site}[{kind}]"
        if account == "collision_reworked_m3":
            placed, handed = float(local["placed_share"]), float(local["handed_share"])
            if placed + handed <= 0.0:
                return {f"{site}:placed": volume}
            return {
                f"{site}:placed": volume * placed / (placed + handed),
                f"{site}:handed": volume * handed / (placed + handed),
            }
    if site in ("cratons.clip_to_column", "cratons.thin_with_column"):
        site = f"{site}<-{_site(frame.f_back)}"
    return {site: volume}


class RecordTap:
    """Sums each split account's bookings by `_split` key."""

    def __init__(self):
        self.by_account: dict[str, dict[str, float]] = {account: defaultdict(float) for account in SPLIT_ACCOUNTS}

    def __enter__(self):
        self.original = cratons.record

        def wrapped(world, account, volume_m3):
            if account in self.by_account and volume_m3:
                for key, part in _split(account, float(volume_m3), sys._getframe(1)).items():
                    self.by_account[account][key] += part
            return self.original(world, account, volume_m3)

        cratons.record = wrapped
        return self

    def __exit__(self, *exc):
        cratons.record = self.original

    def snapshot(self) -> dict:
        return {
            account: {key: value / KM3 for key, value in sorted(sites.items(), key=lambda kv: -kv[1])}
            for account, sites in self.by_account.items()
        }


def _edge_window(plate, cratonic: np.ndarray, hops: int) -> np.ndarray:
    """Non-cratonic continental cells within `hops` cell hops of a cratonic one."""
    continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
    ring = shortening.rings(quad_tectonics._adjacency_matrix(plate), cratonic, hops)
    return continental & ~cratonic & (ring >= 1) & (ring <= hops)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    if values.size == 0:
        return float("nan")
    order = np.argsort(values)
    cum = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cum, q * cum[-1])])


def metrics(world) -> dict:
    out = {"elapsed_myr": world.elapsed_years / 1e6, "craton_mkm3": 0.0, "cratonic_area_mkm2": 0.0, "by_plate_type": {}}
    hops = max(1, round(CRATON_EDGE_RELIEF_WINDOW_KM / (line_spacing_rad(world.node_density) * PLANET_RADIUS_KM)))
    edge_relief_m3 = edge_area_m2 = 0.0
    land_elev, land_area, emergent_m2, buried_craton_m2, speed_m2 = [], [], 0.0, 0.0, 0.0
    total_m2 = 0.0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        areas = plate.node_areas_m2()
        craton = plate.collect("craton_crust_m")
        elevation = plate.collect("elevation")
        land = elevation > world.sea_level_m
        continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        land_elev.append(elevation[land])
        land_area.append(areas[land])
        emergent_m2 += float(areas[land & continental].sum())
        buried_craton_m2 += float(areas[(craton > 0.0) & (plate.collect("crustal_thickness_m") > cratons.CRATON_HOST_MAX_HC_M)].sum())
        points, _ = plate.all_points_and_elevation()
        speed_m2 += float(np.linalg.norm(np.cross(plate.omega, points), axis=1) @ areas)
        total_m2 += float(areas.sum())
        if np.any(craton > 0.0):
            window = _edge_window(plate, craton > 0.0, hops)
            relief = np.maximum(plate.collect("elevation") - world.sea_level_m, 0.0)
            edge_relief_m3 += float(relief[window] @ areas[window])
            edge_area_m2 += float(areas[window].sum())
        volume = float(craton @ areas) / KM3 / 1e6
        area = float(areas[craton > 0.0].sum()) / KM2 / 1e6
        out["craton_mkm3"] += volume
        out["cratonic_area_mkm2"] += area
        row = out["by_plate_type"].setdefault(plate.crust_type, {"craton_mkm3": 0.0, "cratonic_area_mkm2": 0.0})
        row["craton_mkm3"] += volume
        row["cratonic_area_mkm2"] += area
    out["craton_edge_relief_m"] = edge_relief_m3 / edge_area_m2 if edge_area_m2 > 0.0 else 0.0
    elev, area = np.concatenate(land_elev), np.concatenate(land_area)
    out["land_area_mkm2"] = float(area.sum()) / KM2 / 1e6
    out["emergent_continental_area_mkm2"] = emergent_m2 / KM2 / 1e6
    out["land_elevation_m"] = {f"p{q}": _weighted_quantile(elev, area, q / 100.0) for q in (10, 50, 90, 99)}
    out["cratonic_area_under_orogen_mkm2"] = buried_craton_m2 / KM2 / 1e6
    # omega is rad/yr; x R (m) x 100 -> cm/yr.
    out["mean_surface_speed_cm_yr"] = speed_m2 / total_m2 * PLANET_RADIUS_KM * 1e5 if total_m2 > 0.0 else 0.0
    out["ceiling_overflow_km3"] = float((getattr(world, "orogenic_relief_budget", {}) or {}).get("ceiling_overflow_m3", 0.0)) / KM3
    out["merges"] = MERGES[0]
    out["craton_ledger_km3"] = {k: v / KM3 for k, v in (getattr(world, "craton_ledger", {}) or {}).items()}
    out["craton_balance_error_km3"] = cratons.balance_error_m3(world) / KM3
    return out


def non_attrition_loss_km3_per_myr(run: dict) -> float:
    """Every sink's booked craton over the replay, less trench attrition, per Myr."""
    first, last = run["checkpoints"][0], run["checkpoints"][-1]
    lost = sum(last["craton_ledger_km3"].get(k, 0.0) - first["craton_ledger_km3"].get(k, 0.0) for k in cratons.CRATON_SINK_ACCOUNTS)
    account, site = ATTRITION_SITE
    lost -= last.get("by_site_km3", {}).get(account, {}).get(site, 0.0)
    return lost / (last["elapsed_myr"] - first["elapsed_myr"])


def acceptance(run: dict, baseline: dict) -> dict:
    """Issue #340's targets for `run` against `baseline` (same save and span)."""
    loss, base_loss = non_attrition_loss_km3_per_myr(run), non_attrition_loss_km3_per_myr(baseline)
    ratio = loss / base_loss if base_loss > 0.0 else float("inf")
    relief = run["checkpoints"][-1].get("craton_edge_relief_m", 0.0)
    base_relief = baseline["checkpoints"][-1].get("craton_edge_relief_m", 0.0)
    return {
        "non_attrition_loss_km3_per_myr": loss,
        "baseline_non_attrition_loss_km3_per_myr": base_loss,
        "non_attrition_loss_ratio": ratio,
        "non_attrition_loss_pass": ratio <= NON_ATTRITION_LOSS_MAX_FRACTION,
        "craton_edge_relief_m": relief,
        "baseline_craton_edge_relief_m": base_relief,
        "craton_edge_relief_pass": relief > base_relief,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("save", type=Path, nargs="?")
    parser.add_argument("--myr", type=float, default=20.0)
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--checkpoint-myr", type=float, default=2.0)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--baseline", type=Path, help="this script's JSON from main, same save and span")
    parser.add_argument("--rescore", type=Path, help="a finished replay's JSON: recompute the targets against --baseline")
    args = parser.parse_args()
    if args.rescore is not None:
        if args.baseline is None:
            parser.error("--rescore needs --baseline")
        print(json.dumps(acceptance(json.loads(args.rescore.read_text()), json.loads(args.baseline.read_text())), indent=1))
        return
    if args.save is None or args.out is None:
        parser.error("a save and --out are required unless --rescore is given")

    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=False).stdout.strip()
    _count_merges()
    world = persistence.load_world_bytes(args.save.read_bytes())
    world.reset_phase_budget()
    continental_ledger.ensure_initialized(world)
    cratons.ensure_initialized(world)
    steps = int(round(args.myr * 1e6 / args.step_years))
    every = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
    result = {"save": str(args.save), "commit": commit, "myr": args.myr, "checkpoints": [metrics(world)]}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with RecordTap() as tap:
        for i in range(1, steps + 1):
            world_mod.step_world(world, args.step_years)
            if i % every == 0 or i == steps:
                row = metrics(world)
                row["by_site_km3"] = tap.snapshot()
                result["checkpoints"].append(row)
                reworked = row["by_site_km3"]["collision_reworked_m3"]
                print(
                    f"{i}/{steps} {row['elapsed_myr']:.1f} Myr  craton {row['craton_mkm3']:.1f} M km3  "
                    f"area {row['cratonic_area_mkm2']:.2f} M km2  reworked {sum(reworked.values()) / 1e6:.1f} M km3",
                    flush=True,
                )
                args.out.write_text(json.dumps(result, indent=1))
    result["runtime_s"] = time.time() - started
    if args.baseline is not None:
        result["acceptance"] = acceptance(result, json.loads(args.baseline.read_text()))
        print(json.dumps(result["acceptance"], indent=1))
    args.out.write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()

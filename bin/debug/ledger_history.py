#!/usr/bin/env python3
"""Craton and continental-material budgets over a saved world's whole run (rfc/0001-0007).

Reads `World.stats_history` (one row per stats sample, ~10 Myr apart on long runs), so it
needs no replay: the save alone says when cratons and continental crust were lost, through
which ledger account, and whether a sink was steady or came in a few large steps.

Prints three tables:

- checkpoints: craton volume and area, genuine continental area, land fraction, and each
  craton sink's running total;
- per-Gyr craton destruction by sink, as a share of the craton volume at the window start
  (rifting also by path, RFC 0003), plus each sink's largest single-sample jumps (lumpy sinks
  are one-off events);
- continental-material fluxes per window, in km^3/yr: juvenile additions, remelted returns,
  and the terminal sinks, against the mean live inventory.

Usage:
    backend/.venv/bin/python bin/debug/ledger_history.py ~/Downloads/world.mbworld [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import persistence  # noqa: E402

CRATON_SINKS = (
    "rifted_m3", "subducted_m3", "delaminated_m3", "collision_reworked_m3", "eroded_m3",
    "topology_removed_m3", "audit_clipped_m3", "unattributed_m3",
)
# rifted_m3 by path (RFC 0003); zero in saves from before the split.
RIFT_PATHS = ("rifted_divergent_m3", "rifted_stretch_m3", "rifted_failed_m3")
CONTINENTAL_SINKS = (
    "deeply_subducted_m3",
    "rift_thinned_m3",
    "collision_subducted_m3",
    "collision_lower_crust_subducted_m3",
    "delaminated_lower_crust_m3",
    "overloaded_root_delaminated_m3",
    "topology_removed_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
)
CHECKPOINTS_MYR = (0, 100, 200, 500, 1000, 2000, 3000, 4000, 6000, 8000, 10000, 15000, 20000)
FLUX_WINDOWS_MYR = ((0, 100), (100, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 8000), (8000, 12000), (12000, 20000))
MKM3 = 1e6  # stats rows hold km^3; tables print millions of km^3


def _index(t_myr: np.ndarray, at_myr: float) -> int:
    return int(min(np.searchsorted(t_myr, at_myr), len(t_myr) - 1))


def summarize(history: list[dict]) -> dict:
    t = np.array([row["elapsed_years"] for row in history]) / 1e6
    craton = {key: np.array([row["craton_ledger_km3"].get(key, 0.0) for row in history]) for key in (*CRATON_SINKS, *RIFT_PATHS, "formed_m3")}
    volume = np.array([row["craton_volume_km3"] for row in history])
    material = [row["continental_material_ledger"] for row in history]

    checkpoints = []
    for at in sorted({_index(t, m) for m in CHECKPOINTS_MYR}):
        row = history[at]
        fraction = row["craton_continental_fraction"]
        checkpoints.append({
            "myr": float(t[at]),
            "craton_volume_mkm3": volume[at] / MKM3,
            "craton_area_mkm2": row["craton_area_km2"] / 1e6,
            "genuine_continental_area_mkm2": row["craton_area_km2"] / fraction / 1e6 if fraction > 0.0 else 0.0,
            "land_fraction": row["land_fraction"],
            **{f"craton_{key}_mkm3": craton[key][at] / MKM3 for key in (*CRATON_SINKS, "formed_m3")},
        })

    per_gyr = []
    for start in range(0, int(t[-1] // 1000) + 1):
        # Include the endpoint: subtracting one silently dropped the last sample's
        # loss, including potentially a large merge, from every Gyr (#403).
        i, j = _index(t, start * 1000), _index(t, (start + 1) * 1000)
        # Once the cratons are nearly gone, per-window shares are noise.
        if j <= i or volume[i] <= 1e-3 * volume[0]:
            continue
        duration_gyr = (t[j] - t[i]) / 1000
        shares = {key: (craton[key][j] - craton[key][i]) / volume[i] / duration_gyr for key in CRATON_SINKS}
        paths = {key: (craton[key][j] - craton[key][i]) / volume[i] / duration_gyr for key in RIFT_PATHS}
        per_gyr.append({"gyr": start, "window_myr": [float(t[i]), float(t[j])], "start_volume_mkm3": volume[i] / MKM3,
                        **shares, "total": sum(shares.values()), **paths})

    jumps = {}
    for key in (*CRATON_SINKS, "formed_m3"):
        step = np.diff(craton[key])
        order = np.argsort(step)[::-1][:5]
        jumps[key] = [
            {"myr": float(t[k + 1]), "mkm3": step[k] / MKM3, "share_of_live": step[k] / volume[k] if volume[k] > 0.0 else 0.0}
            for k in order if step[k] > 0.0
        ]

    fluxes = []
    for lo, hi in FLUX_WINDOWS_MYR:
        i, j = _index(t, lo), _index(t, hi)
        years = (t[j] - t[i]) * 1e6
        if years <= 0.0:
            continue

        def rate(key: str, i: int = i, j: int = j, years: float = years) -> float:
            return (material[j].get(key, 0.0) - material[i].get(key, 0.0)) / 1e9 / years

        fluxes.append({
            "window_myr": [float(t[i]), float(t[j])],
            "juvenile_km3_per_yr": rate("juvenile_additions_m3"),
            "returns_km3_per_yr": rate("remelted_relaminated_returns_m3"),
            "sinks_km3_per_yr": sum(rate(key) for key in CONTINENTAL_SINKS),
            "deeply_subducted_km3_per_yr": rate("deeply_subducted_m3"),
            "mean_inventory_mkm3": float(np.mean([m["surface_continental_derived_m3"] for m in material[i : j + 1]])) / 1e15,
            "on_oceanic_hosts_fraction": material[j]["continental_sediment_on_oceanic_hosts_m3"]
            / max(material[j]["surface_continental_derived_m3"], 1.0),
        })
    return {"checkpoints": checkpoints, "craton_loss_per_gyr": per_gyr, "craton_largest_jumps": jumps, "continental_fluxes": fluxes}


def _print(summary: dict) -> None:
    print(f"{'Myr':>7} {'craton':>8} {'cr.area':>8} {'cont.area':>9} {'land':>6} " + " ".join(f"{k.removesuffix('_m3')[:9]:>9}" for k in (*CRATON_SINKS, "formed_m3")))
    print(f"{'':>7} {'Mkm3':>8} {'Mkm2':>8} {'Mkm2':>9} {'':>6}   (craton ledger running totals, Mkm3)")
    for row in summary["checkpoints"]:
        print(
            f"{row['myr']:7.0f} {row['craton_volume_mkm3']:8.0f} {row['craton_area_mkm2']:8.1f} {row['genuine_continental_area_mkm2']:9.1f} "
            f"{row['land_fraction']:6.3f} " + " ".join(f"{row[f'craton_{k}_mkm3']:9.0f}" for k in (*CRATON_SINKS, "formed_m3"))
        )
    print("\ncraton destruction per Gyr, share of the volume at the measured window start (unattributed is signed)")
    for row in summary["craton_loss_per_gyr"]:
        lo, hi = row["window_myr"]
        print(f"{lo / 1000:5g}-{hi / 1000:<5g} Gyr V0={row['start_volume_mkm3']:7.0f} Mkm3 " + " ".join(f"{k.removesuffix('_m3')[:9]}={row[k]:4.0%}" for k in CRATON_SINKS) + f"  total={row['total']:.0%}"
              + "  rifted by path: " + " ".join(f"{k.removeprefix('rifted_').removesuffix('_m3')}={row[k]:4.1%}" for k in RIFT_PATHS))
    print("\nlargest single-sample craton ledger jumps")
    for key, rows in summary["craton_largest_jumps"].items():
        print(f"  {key:24s} " + ", ".join(f"{r['myr']:.0f} My: {r['mkm3']:.0f} Mkm3 ({r['share_of_live']:.0%})" for r in rows))
    print("\ncontinental material fluxes (km^3/yr)")
    print(f"{'window Myr':>14} {'juvenile':>9} {'returns':>8} {'sinks':>7} {'deep.sub':>8} {'inventory Mkm3':>14} {'on ocean':>8}")
    for row in summary["continental_fluxes"]:
        lo, hi = row["window_myr"]
        print(
            f"{lo:6.0f}-{hi:7.0f} {row['juvenile_km3_per_yr']:9.2f} {row['returns_km3_per_yr']:8.2f} {row['sinks_km3_per_yr']:7.2f} "
            f"{row['deeply_subducted_km3_per_yr']:8.2f} {row['mean_inventory_mkm3']:14.0f} {row['on_oceanic_hosts_fraction']:8.0%}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("save", type=Path)
    parser.add_argument("--json", type=Path, help="also write the summary here")
    args = parser.parse_args()
    world = persistence.load_world_bytes(args.save.read_bytes())
    summary = summarize(world.stats_history)
    _print(summary)
    if args.json:
        args.json.write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()

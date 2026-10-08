#!/usr/bin/env python3
"""Issue #297 part 1: measure least-cost depression breaching (app/breaching.py) on real runs.

For each seed this generates the world twice, once with `hydrology.BREACH_DEPRESSIONS_ENABLED`
off and once with it on, and steps both identically. Every hydrology call records:

- pits (land nodes with no lower neighbour), breached pits, pits too costly to carve, and
  pits cheap enough to carve but kept closed by their water balance (endorheic);
- invented notch depth (sum and max) and how many nodes were notched;
- boundary passes lowered below the higher centre: hierarchy edges where the passage
  elevation lets water cross below `max(z_i, z_j)`, from an established channel or a new notch;
- standing lakes: flooded node count, flooded area (whole cells with > 1 m of water, and the
  partial-cell wet area from lakes.CellHypsometry), the visible area in lakes filled to their
  rim and spilling, and water volume;
- wall time of the breaching pass and of the whole hydrology call.

The breaching statistics are computed in the "off" run too, as a diagnostic that doesn't feed
back. The "off" run shows how many pits *would* be breached on unchanged terrain.

    cd backend
    .venv/bin/python ../bin/debug/measure_breaching.py --seeds 1 2 --steps 50 --years 100000
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--steps", type=int, default=50)
    parser.add_argument("--years", type=float, default=100_000.0)
    parser.add_argument("--every", type=int, default=10, help="print every Nth step")
    parser.add_argument("--modes", default="off,on", help="breaching modes to run, e.g. 'on' for an A/B against another --backend")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    sys.path.insert(0, str(args.backend.resolve()))
    from app import breaching, cratons, erosion, hydrology, lakes  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    original_compute = hydrology.compute_hydrology
    original_breach = breaching.breach_depressions
    current: dict = {}

    def timed_breach(*a, **kw):
        start = time.perf_counter()
        result = original_breach(*a, **kw)
        current["breach_s"] = current.get("breach_s", 0.0) + time.perf_counter() - start
        current["breach"] = result
        return result

    def wrapped(world, precipitation, temperature, years, **kw):
        current.clear()
        start = time.perf_counter()
        fields = original_compute(world, precipitation, temperature, years, **kw)
        current["hydrology_s"] = time.perf_counter() - start
        breach = current.get("breach")
        if breach is None:
            # Breaching was off: compute it as a pure diagnostic on the same inputs.
            plates = fields.plates_in_order
            strength = cratons.strength(np.concatenate([p.collect("craton_crust_m") for p in plates]))
            cover = np.concatenate([p.collect("mobile_cover_m") for p in plates])
            channel = np.concatenate([p.collect("channel_depth") for p in plates])
            prior_notch = np.concatenate([p.collect("breach_notch_depth_m") for p in plates])
            areas = erosion._gather_areas(world, plates)
            breach = original_breach(
                fields.elevation, fields.is_ocean, fields.neighbor_idx, channel,
                breaching.bedrock_carve_rate_m_per_myr(strength), years,
                water=breaching.WaterBalance(
                    np.where(temperature < hydrology.FREEZE_POINT_C, 0.0, precipitation), temperature, areas,
                    loss_fraction=hydrology.river_transit_loss_fraction(
                        temperature, fields.is_ocean | (temperature < hydrology.FREEZE_POINT_C)
                    ),
                ),
                mobile_cover_m=cover,
                prior_notch_m=prior_notch,
            )
        current["record"] = _record(fields, breach, erosion._gather_areas(world, fields.plates_in_order), current, hydrology, lakes)
        return fields

    breaching.breach_depressions = timed_breach
    hydrology.compute_hydrology = wrapped

    results: dict = {}
    for seed in args.seeds:
        for enabled in [mode.strip() == "on" for mode in args.modes.split(",")]:
            hydrology.BREACH_DEPRESSIONS_ENABLED = enabled
            label = f"seed={seed} breaching={'on' if enabled else 'off'}"
            world = generate_world(seed=seed)
            rows = []
            for step in range(1, args.steps + 1):
                step_world(world, args.years)
                record = dict(current.get("record", {}), step=step)
                rows.append(record)
                if step % args.every == 0 or step == args.steps:
                    print(f"{label} step {step}: " + ", ".join(f"{k}={_fmt(v)}" for k, v in record.items() if k != "step"), flush=True)
            results[label] = rows
            print(f"{label} mean: " + ", ".join(f"{k}={_fmt(v)}" for k, v in _means(rows).items()), flush=True)

    if args.out:
        args.out.write_text(json.dumps(results, indent=1))


def _record(fields, breach, area_m2, current, hydrology, lakes) -> dict:
    land = ~fields.is_ocean
    z = fields.elevation
    nb = fields.neighbor_idx
    passes = np.maximum(breach.passage_m[:, None], breach.passage_m[nb])
    centre_max = np.maximum(z[:, None], z[nb])
    lowered = land[:, None] & (passes < centre_max - 1.0)
    lake = fields.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    notched = breach.notch_m > 0.0
    return {
        "land_nodes": int(land.sum()),
        "pits": int(len(breach.breached_pits) + len(breach.closed_pits) + len(breach.endorheic_pits)),
        "breached_pits": int(len(breach.breached_pits)),
        "closed_pits": int(len(breach.closed_pits)),
        "endorheic_pits": int(len(breach.endorheic_pits)),
        "notched_nodes": int(notched.sum()),
        "notch_sum_m": float(breach.notch_m.sum()),
        "notch_max_m": float(breach.notch_m.max(initial=0.0)),
        "lowered_pass_edges": int(lowered.sum()),
        "lowered_pass_share": float(lowered.sum() / max(int(land.sum()) * nb.shape[1], 1)),
        "lake_nodes": int(lake.sum()),
        "lake_area_km2": float(area_m2[lake].sum() / 1e6),
        "lake_wet_area_km2": _wet_area_km2(fields, area_m2, lakes),
        "lake_at_rim_area_km2": _at_rim_area_km2(fields, area_m2, hydrology, lakes),
        "lake_volume_km3": float((fields.lake_depth * area_m2).sum() / 1e9),
        "breach_s": float(current.get("breach_s", 0.0)),
        "hydrology_s": float(current["hydrology_s"]),
    }


def _at_rim_area_km2(fields, area_m2, hydrology, lakes) -> float:
    """Visible lake area (cells > 1 m) in lakes filled to their rim and spilling this step."""
    at_rim = np.zeros(len(area_m2), dtype=bool)
    for lake in lakes.iter_all_lakes(fields.lake_forest):
        if lake.is_spilling:
            at_rim[lake.members] = True
    visible = fields.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    return float(area_m2[at_rim & visible].sum() / 1e6)


def _wet_area_km2(fields, area_m2, lakes) -> float:
    """Flooded area counting partly flooded cells by their wet fraction. Falls back to whole
    cells on a backend without lakes.CellHypsometry."""
    wet = np.flatnonzero(fields.lake_depth > 0.0)
    if not hasattr(lakes, "CellHypsometry"):
        return float(area_m2[wet].sum() / 1e6)
    hyps = lakes.CellHypsometry.from_terrain(fields.elevation, area_m2, fields.neighbor_idx)
    levels = hyps.level_from_depth(fields.lake_depth[wet], wet)
    fraction = np.array([hyps.wet_fraction(h, np.array([i]))[0] for h, i in zip(levels, wet)])
    return float((fraction * area_m2[wet]).sum() / 1e6)


def _means(rows: list[dict]) -> dict:
    keys = [k for k in rows[0] if k != "step"]
    return {k: float(np.mean([r[k] for r in rows])) for k in keys}


def _fmt(v) -> str:
    return f"{v:.3g}" if isinstance(v, float) else str(v)


if __name__ == "__main__":
    main()

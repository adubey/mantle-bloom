#!/usr/bin/env python3
"""#822: compare the glacier-march cost options on one world, one preset per run.

Steps a world N steps of `--myr` Myr with `step_world` under one `ice_sheet.MARCH_COST`
preset (`d2` is the production march) and records, per step: whole-step and march wall
time, the march log entry (substeps, capped, largest substep, sea-level re-solves),
the ice ledger (closing volume, accumulation, melt, residual), the land ledger closure
and sea level. The closing glacier depth of every step goes to an `.npz` beside the JSON
so `summarize_822_march_cost.py` can diff each preset against the `d2` arm.

Make a shared start state from a generated world first (an exact save, so every arm
opens from identical bytes), then run each preset as its own process, one at a time:

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_822_march_cost.py \\
        --generate '{"seed": 5, "node_density": 0.25}' --write-save <scratch>/gen5.mbworld
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_822_march_cost.py \\
        <save> --preset cap-hold --steps 1 --myr 0.25 --out <scratch>/gen5-cap-hold.json

`--set key=value` overrides one MarchCost field on top of the preset (repeatable), e.g.
`--set max_substeps=400`. `--substep-years 125` reruns an arm at a shorter base substep.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import persistence  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, ice_sheet  # noqa: E402


def _parse_override(text: str) -> tuple[str, object]:
    key, _, raw = text.partition("=")
    fields = {f.name: f for f in dataclasses.fields(ice_sheet.MarchCost)}
    if key not in fields:
        raise SystemExit(f"unknown MarchCost field {key!r}; known: {sorted(fields)}")
    value = json.loads(raw) if raw not in ("hold", "endpoint") else raw
    return key, value


def _generate(spec: str, out: Path) -> None:
    options = json.loads(spec)
    seed = options.pop("seed")
    world = world_mod.generate_world(seed, **options)
    out.write_bytes(persistence.save_exact_bytes(world))
    print(f"wrote {out} sha256={hashlib.sha256(out.read_bytes()).hexdigest()}")


def _timed_march():
    """Wrap hydrology._march_glacier to sum its wall time per step."""
    clock = {"s": 0.0}
    real = hydrology._march_glacier

    def timed(*args, **kwargs):
        start = time.perf_counter()
        try:
            return real(*args, **kwargs)
        finally:
            clock["s"] += time.perf_counter() - start

    hydrology._march_glacier = timed
    return clock


def _step_record(i, world, wall, march_s):
    land = world.hydrology_cache.land_water_ledger
    ice = world.hydrology_cache.ice.ledger
    external = land.rain_m3 + land.meltwater_m3
    march = dict(ice_sheet.MARCH_LOG[-1]) if ice_sheet.MARCH_LOG else None
    if march is not None:
        march.pop("volume_path_km3", None)
    loss = ice.surface_melt_m3 + ice.basal_melt_m3 + ice.calved_m3 + ice.overflow_m3 + ice.capped_m3
    return {
        "step": i,
        "wall_s": round(wall, 2),
        "march_s": round(march_s, 2),
        "march": march,
        "ice_opening_km3": ice.opening_m3 / 1e9,
        "ice_closing_km3": ice.closing_m3 / 1e9,
        "ice_accumulation_km3": ice.accumulation_m3 / 1e9,
        "ice_loss_km3": loss / 1e9,
        "ice_residual_km3": ice.residual_m3 / 1e9,
        "ice_closure": abs(ice.residual_m3) / max(ice.accumulation_m3, 1.0),
        "land_closure": abs(land.residual_m3) / max(external, 1.0),
        "sea_level_m": float(world.sea_level_m),
    }


def main() -> None:
    args = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    args.add_argument("save", type=Path, nargs="?")
    args.add_argument("--generate", help="JSON generate_world options incl. seed; writes --write-save")
    args.add_argument("--write-save", type=Path)
    args.add_argument("--preset", default="d2", choices=sorted(ice_sheet.MARCH_COST_PRESETS))
    args.add_argument("--set", action="append", default=[], metavar="FIELD=VALUE")
    args.add_argument("--steps", type=int, default=1)
    args.add_argument("--myr", type=float, default=0.25)
    args.add_argument("--out", type=Path)
    args.add_argument("--substep-years", type=float, help="override MARCH_SUBSTEP_YEARS (a d2 convergence control)")
    ns = args.parse_args()

    if ns.generate:
        if ns.write_save is None:
            raise SystemExit("--generate needs --write-save")
        _generate(ns.generate, ns.write_save)
        return
    if ns.save is None or ns.out is None:
        raise SystemExit("a save and --out are required")

    cost = dataclasses.replace(ice_sheet.march_cost(ns.preset), **dict(_parse_override(t) for t in ns.set))
    ice_sheet.MARCH_COST = cost
    if ns.substep_years:
        ice_sheet.MARCH_SUBSTEP_YEARS = ns.substep_years
    raw = ns.save.read_bytes()
    world = persistence.load_world_bytes(raw)
    ice_sheet.MARCH_LOG.clear()
    clock = _timed_march()

    steps, depths = [], {}
    for i in range(ns.steps):
        clock["s"] = 0.0
        start = time.perf_counter()
        world_mod.step_world(world, ns.myr * 1e6)
        wall = time.perf_counter() - start
        steps.append(_step_record(i, world, wall, clock["s"]))
        depths[f"step{i}"] = np.asarray(world.hydrology_cache.glacier_depth, dtype=float)
        rec = steps[-1]
        m = rec["march"] or {}
        print(
            f"[{ns.preset}] step {i}: {wall:.1f}s (march {rec['march_s']:.1f}s) substeps={m.get('substeps')} "
            f"settled={m.get('settled')} capped={m.get('capped')} ice={rec['ice_closing_km3']:.1f} km3 "
            f"ice_closure={rec['ice_closure']:.2e} sea={rec['sea_level_m']:.2f}",
            flush=True,
        )

    np.savez_compressed(ns.out.with_suffix(".npz"), **depths)
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    ns.out.write_text(json.dumps({
        "save": ns.save.name,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "preset": ns.preset,
        "cost": dataclasses.asdict(cost),
        "steps": ns.steps,
        "myr_per_step": ns.myr,
        "peak_rss_gb": round(peak_gb, 2),
        "per_step": steps,
    }, indent=1) + "\n")
    print(f"wrote {ns.out} peak_rss_gb={peak_gb:.2f}")


if __name__ == "__main__":
    main()

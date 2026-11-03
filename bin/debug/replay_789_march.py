#!/usr/bin/env python3
"""RFC 0015 #789: replay a saved world step by step, recording the march and ledgers.

Steps a `.mbworld` save forward N steps of `--myr` Myr each with `step_world`, recording
per step: wall time, the land water ledger's residual against its external inputs, the ice
ledger's residual against its accumulation, the march log entry (substeps, marched/held
years, settlement, volume path), and the closing ice volume and sea level. Used for the
D1-D2 acceptance replay: every step must close the land ledger to 0.1% of its external
inputs, and the march must settle on the diagnosis replays.

    cd backend
    PYTHONPATH=backend .venv/bin/python ../bin/debug/replay_789_march.py \
        <save> --steps 10 --myr 0.25 --out <scratch-output>.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import persistence  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.hydroclimate import ice_sheet  # noqa: E402


def main() -> None:
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("save", type=Path)
    args.add_argument("--steps", type=int, default=10)
    args.add_argument("--myr", type=float, default=0.25)
    args.add_argument("--out", type=Path, required=True)
    ns = args.parse_args()

    raw = ns.save.read_bytes()
    world = persistence.load_world_bytes(raw)
    ice_sheet.MARCH_LOG.clear()

    steps = []
    for i in range(ns.steps):
        start = time.perf_counter()
        world_mod.step_world(world, ns.myr * 1e6)
        elapsed = time.perf_counter() - start
        land = world.hydrology_cache.land_water_ledger
        ice = world.hydrology_cache.ice.ledger
        external = land.rain_m3 + land.meltwater_m3
        march = dict(ice_sheet.MARCH_LOG[-1])
        steps.append({
            "step": i,
            "wall_s": round(elapsed, 1),
            "land_residual_m3": land.residual_m3,
            "land_external_m3": external,
            "land_closure": abs(land.residual_m3) / max(external, 1.0),
            "ice_residual_m3": ice.residual_m3,
            "ice_accumulation_m3": ice.accumulation_m3,
            "ice_closure": abs(ice.residual_m3) / max(ice.accumulation_m3, 1.0),
            "ice_closing_km3": round(ice.closing_m3 / 1e9, 3),
            "sea_level_m": round(float(world.sea_level_m), 2),
            "march": march,
        })
        print(
            f"step {i}: {elapsed:.0f}s settled={march['settled']} "
            f"substeps={march['substeps']} held={march['held_years']:.0f} "
            f"land_closure={steps[-1]['land_closure']:.2e} "
            f"ice_closure={steps[-1]['ice_closure']:.2e}",
            flush=True,
        )

    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    ns.out.write_text(json.dumps({
        "save": ns.save.name,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "steps": ns.steps,
        "myr_per_step": ns.myr,
        "peak_rss_gb": round(peak_gb, 2),
        "per_step": steps,
    }, indent=1) + "\n")
    print(f"wrote {ns.out} peak_rss_gb={peak_gb:.2f}")


if __name__ == "__main__":
    main()

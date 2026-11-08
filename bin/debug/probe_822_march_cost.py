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
`--set max_substeps=400`. `--substep-years 125` reruns an arm at a shorter base substep, and
`--ice-age-period-myr 8` turns ice ages on (bin/debug/march_cli.py). `--ensemble` (#842 E3)
re-marches the step after the `--steps` warm-up steps once per #765 opening-ice member
(depths scaled by +/-0.01-0.5%), recording each member's closing ice and largest
single-substep volume move.
"""

from __future__ import annotations

import argparse
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from march_cli import add_march_args, apply_march_args, apply_world_args  # noqa: E402

# #765's opening-ice ensemble: an unperturbed member 0, then 12 members scaled by
# +/-0.01-0.5% (6 magnitudes, both signs).
ENSEMBLE_SCALES = [0.0] + [sign * m for m in (1e-4, 3e-4, 1e-3, 2e-3, 3.5e-3, 5e-3) for sign in (1.0, -1.0)]


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
    add_march_args(args)
    args.add_argument("--ensemble", action="store_true",
                      help="re-march the step after --steps warm-up steps once per ENSEMBLE_SCALES member")
    args.add_argument("--steps", type=int, default=1)
    args.add_argument("--myr", type=float, default=0.25)
    args.add_argument("--out", type=Path)
    ns = args.parse_args()

    if ns.generate:
        if ns.write_save is None:
            raise SystemExit("--generate needs --write-save")
        _generate(ns.generate, ns.write_save)
        return
    if ns.save is None or ns.out is None:
        raise SystemExit("a save and --out are required")

    choice = apply_march_args(ns, ice_sheet)
    raw = ns.save.read_bytes()
    world = persistence.load_world_bytes(raw)
    apply_world_args(ns, world)
    ice_sheet.MARCH_LOG.clear()
    clock = _timed_march()

    steps, depths = _run_steps(ns, world, clock)
    members = _run_ensemble(ns, world, clock) if ns.ensemble else None

    np.savez_compressed(ns.out.with_suffix(".npz"), **depths)
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    ns.out.write_text(json.dumps({
        "save": ns.save.name,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "preset": ns.preset,
        "march_choice": choice,
        "cost": choice["cost"],
        "steps": ns.steps,
        "myr_per_step": ns.myr,
        "peak_rss_gb": round(peak_gb, 2),
        "per_step": steps,
        "ensemble": members,
    }, indent=1) + "\n")
    print(f"wrote {ns.out} peak_rss_gb={peak_gb:.2f}")


def _run_steps(ns, world, clock):
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
    return steps, depths


def _run_ensemble(ns, world, clock):
    """#765's opening-ice ensemble: from the state after the warm-up steps, scale every
    node's glacier depth by (1 + s) for each s in ENSEMBLE_SCALES and step once."""
    snapshot = persistence.save_exact_bytes(world)
    members = []
    for k, scale in enumerate(ENSEMBLE_SCALES):
        member = persistence.load_world_bytes(snapshot)
        for plate in member.plates:
            if plate.node_count() > 0:
                plate.set_fields_on_plate(glacier_depth=plate.collect("glacier_depth") * (1.0 + scale))
        clock["s"] = 0.0
        start = time.perf_counter()
        world_mod.step_world(member, ns.myr * 1e6)
        rec = _step_record(ns.steps, member, time.perf_counter() - start, clock["s"])
        rec["member"], rec["scale"] = k, scale
        members.append(rec)
        m = rec["march"] or {}
        print(
            f"[{ns.preset}] member {k} ({scale:+.2%}): ice={rec['ice_closing_km3']:.1f} km3 "
            f"substeps={m.get('substeps')} max_substep_change={m.get('max_substep_change')} "
            f"settle={m.get('settle_test')}",
            flush=True,
        )
    return members

if __name__ == "__main__":
    main()

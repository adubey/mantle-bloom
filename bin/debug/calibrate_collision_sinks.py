#!/usr/bin/env python3
"""Issue #276: how much continental crust collisions lose, per step and compounded, under
the current constants or any override of them.

Steps a generated world (`--seeds`) or a save (`--save`) and appends one JSON line per
snapshot to `--out`. Every snapshot reports both rates the issue asks to keep apart:

- **retreat-processed loss** -- of the Hc volume suture retreat donated (`suture_donated_m3`),
  the share lost: eligible-root delamination plus the terminal no-outlet remainder, which
  subducts. This is the issue's "4.15% of retreat-processed material".
- **whole-inventory loss** -- the continental-material ledger's collision sinks
  (`delaminated_lower_crust_m3`, `collision_subducted_m3`, `overloaded_root_delaminated_m3`)
  as a share of the live
  inventory at the start of each step, compounded step by step (`collision_survival`), and
  the same for every sink together (`all_sink_survival`, which adds trench subduction,
  rifting, discarded sediment, numerical clipping and dropped plate fragments). `*_500myr_projection` extrapolates the run's
  mean per-Myr rate over 500 Myr, as #272's "~72% over 500 Myr" did.

Alongside those: land and continental-host area, live material against the start, the
continental Hc distribution and area at the Hc cap (runaway thickening), the orogenic
relief budget, and the convergent ceiling overflow's residue (Hc the shortening would have
made, which never enters the material tracer).

`--set NAME=VALUE` overrides a module constant before the run. It is applied to every
`app.*` module holding that name, so a constant imported by value (e.g.
`SUTURE_ACCRETION_MAX_HC_M`, which `quad_tectonics` and `quad_merge` import from
`lithosphere_plate`) changes everywhere it is read:

    cd backend
    .venv/bin/python ../bin/debug/calibrate_collision_sinks.py --seeds 1 2 3 --steps 100 --out ../results/276.jsonl
    .venv/bin/python ../bin/debug/calibrate_collision_sinks.py --save ~/Downloads/x.mbworld --steps 50 \\
        --set SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION=0.05 --label frac05 --out ../results/276.jsonl
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

COLLISION_SINKS = ("delaminated_lower_crust_m3", "collision_subducted_m3", "overloaded_root_delaminated_m3")
ALL_SINKS = COLLISION_SINKS + (
    "deeply_subducted_m3",
    "rift_thinned_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
    "topology_removed_m3",
)


def apply_overrides(settings: list[str]) -> dict[str, float]:
    import app

    modules = [m for name, m in sys.modules.items() if name == "app" or name.startswith("app.")]
    applied = {}
    for setting in settings:
        name, _, raw = setting.partition("=")
        value = int(raw) if raw.lstrip("-").isdigit() else float(raw)
        holders = [m for m in modules if hasattr(m, name)]
        if not holders:
            raise SystemExit(f"no app module defines {name}")
        for module in holders:
            setattr(module, name, value)
        applied[name] = value
    del app
    return applied


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--seeds", type=int, nargs="+")
    source.add_argument("--save", type=Path)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--every", type=int, default=10, help="snapshot every Nth step")
    parser.add_argument("--set", dest="settings", action="append", default=[], metavar="NAME=VALUE")
    parser.add_argument("--label", default="baseline")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    sys.path.insert(0, str(args.backend.resolve()))
    from app import continental_ledger, lithosphere, orogeny, persistence, quad_tectonics  # noqa: E402,F401
    from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    overrides = apply_overrides(args.settings)
    from app import lithosphere_plate  # noqa: E402

    cap = lithosphere_plate.SUTURE_ACCRETION_MAX_HC_M
    runs = [("save", args.save)] if args.save else [("seed", seed) for seed in args.seeds]
    for kind, key in runs:
        if kind == "save":
            world = persistence.load_world_bytes(key.expanduser().read_bytes())
            run_id = key.name
        else:
            world = generate_world(seed=key)
            run_id = f"seed{key}"
        continental_ledger.ensure_initialized(world)
        budget = orogeny.ensure_budget(world)
        ledger = world.continental_material_ledger
        start_inventory = continental_ledger.surface_volume_m3(world)
        start_budget = dict(budget)
        start_ledger = dict(ledger)
        collision_survival = all_survival = 1.0
        first_hc_volume = None
        start = time.perf_counter()
        for step in range(1, args.steps + 1):
            inventory = continental_ledger.surface_volume_m3(world)
            before = {k: ledger[k] for k in ALL_SINKS}
            step_world(world, args.years)
            if inventory > 0.0:
                collision_survival *= 1.0 - sum(ledger[k] - before[k] for k in COLLISION_SINKS) / inventory
                all_survival *= 1.0 - sum(ledger[k] - before[k] for k in ALL_SINKS) / inventory
            if step % args.every and step != args.steps:
                continue

            spacing = line_spacing_rad(world.node_density)
            hc_parts, area_parts, cont_parts = [], [], []
            for plate in world.plates:
                if plate.node_count() == 0:
                    continue
                hc_parts.append(plate.collect("crustal_thickness_m"))
                area_parts.append(plate.accounting_areas_m2(spacing))
                cont_parts.append(
                    effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
                )
            hc, areas, continental = (np.concatenate(p) for p in (hc_parts, area_parts, cont_parts))
            total_area = float(areas.sum())
            hydro = world.hydrology_cache
            land_area = (
                float(np.sum(np.concatenate([p.accounting_areas_m2(spacing) for p in hydro.plates_in_order])[~hydro.is_ocean]))
                if hydro is not None
                else float("nan")
            )
            cont_hc = hc[continental]
            cont_hc_volume = float(cont_hc @ areas[continental])
            first_hc_volume = first_hc_volume or cont_hc_volume
            d_budget = {k: budget[k] - start_budget.get(k, 0.0) for k in budget}
            d_ledger = {k: ledger[k] - start_ledger.get(k, 0.0) for k in ledger}
            donated = d_budget["suture_donated_m3"]
            suture_lost = d_budget["delamination_completed_m3"] + d_budget["no_outlet_subducted_m3"]
            myr = step * args.years / 1e6
            collision_rate = 1.0 - collision_survival ** (1.0 / myr)
            all_rate = 1.0 - all_survival ** (1.0 / myr)
            live = continental_ledger.surface_volume_m3(world)
            row = {
                "run": run_id,
                "label": args.label,
                "overrides": overrides,
                "step": step,
                "seconds": round(time.perf_counter() - start, 1),
                "balance_error_vs_inventory": continental_ledger.balance_error_m3(world) / start_inventory,
                # Retreat-processed (Hc volume through suture placement).
                "suture_donated_vs_inventory": donated / start_inventory,
                "retreat_processed_loss_share": suture_lost / donated if donated > 0.0 else 0.0,
                "root_delamination_share": d_budget["delamination_completed_m3"] / donated if donated > 0.0 else 0.0,
                "no_outlet_share": d_budget["no_outlet_subducted_m3"] / donated if donated > 0.0 else 0.0,
                # Whole-inventory (continental material).
                "collision_sink_vs_start": sum(d_ledger[k] for k in COLLISION_SINKS) / start_inventory,
                "collision_survival": collision_survival,
                "collision_rate_per_myr": collision_rate,
                "collision_500myr_projection": 1.0 - (1.0 - collision_rate) ** 500,
                "all_sink_survival": all_survival,
                "all_sink_rate_per_myr": all_rate,
                "all_sink_500myr_projection": 1.0 - (1.0 - all_rate) ** 500,
                **{f"{k[:-3]}_vs_start": d_ledger[k] / start_inventory for k in ALL_SINKS if k not in COLLISION_SINKS},
                "juvenile_vs_start": d_ledger["juvenile_additions_m3"] / start_inventory,
                "live_material_vs_start": live / start_inventory,
                # Hc volume the convergent ceiling dropped (never in the material tracer).
                "ceiling_residue_vs_inventory": d_budget["ceiling_overflow_residue_m3"] / start_inventory,
                # Area and thickness.
                "land_fraction": land_area / total_area,
                "continental_area_fraction": float(areas[continental].sum()) / total_area,
                "continental_hc_mean_m": float(np.average(cont_hc, weights=areas[continental])) if cont_hc.size else 0.0,
                "continental_hc_p99_m": float(np.percentile(cont_hc, 99)) if cont_hc.size else 0.0,
                "hc_at_cap_area_fraction": float(areas[hc >= cap - 1.0].sum()) / total_area,
                "continental_hc_volume_vs_first": cont_hc_volume / first_hc_volume,
                **{f"budget_{k}": v for k, v in d_budget.items()},
            }
            print(json.dumps({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
            if args.out:
                args.out.parent.mkdir(parents=True, exist_ok=True)
                with args.out.open("a") as handle:
                    handle.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()

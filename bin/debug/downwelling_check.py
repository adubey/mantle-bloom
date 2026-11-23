"""Does subduction land on the model's downwelling centres? (#998, exploration)

The static mantle-flow field is a set of signed Gaussian centres (`mantle.ConvectionCenter`: strength > 0 upwelling,
< 0 downwelling) that never see the plates. Nothing in the model makes a trench form where a downwelling is. This
measures whether it does anyway. For each sample time of one seed it records, over the plates' boundary band:
  * `D` at every node, the downwelling potential  D(x) = -sum_c strength_c * exp(-(angle(x, c) / falloff_c)^2)
    (positive over downwellings, negative over upwellings), in units of the reference flow rate;
  * the mean D at subducting nodes (`torque.subducting_boundary_mask`), at ridge nodes (divergent) and at all band nodes;
  * a null: the same mean D at subducting nodes, with the centres redrawn from `generate_convection_centers`
    (200 random sets, same count and strength distribution), so the statistic is a percentile among chance placements.
Run per seed, then aggregate:

    intervention: `--flow-scale 0` (flow off) or `--flow-scale -4` (every centre's sign flipped) with the default
    `--reference-scale 4`; trench placement is still scored against the original centres.
    run        backend/.venv/bin/python bin/debug/downwelling_check.py run --settings analysis/issue998/settings.json \\
                   --setting b123-brake-0.5-for-flow4 --seed 591453549 --at 0 5 20 50 --flow-scale 4 --out <captures>/dw/591453549.json
    aggregate  backend/.venv/bin/python bin/debug/downwelling_check.py aggregate --runs <captures>/dw --out analysis/issue998/downwelling-check.json

Exploration only; no `backend/app` change. About 1 GB per seed (compute-budget first).
"""

import argparse
import json

import numpy as np
from sampler_common import boundary_context, build_parser, build_world, plate_boundary, reference_centers, scale_flow, step_and_sample

NULL_SETS = 200


def downwelling_potential(points: np.ndarray, centers) -> np.ndarray:
    total = np.zeros(len(points))
    for c in centers:
        angle = np.arccos(np.clip(points @ c.position, -1.0, 1.0))
        total -= c.strength * np.exp(-((angle / c.falloff) ** 2))
    return total


def band_nodes(world) -> dict:
    from app import torque

    ctx = boundary_context(world)
    out = {"all": [], "subducting": [], "ridge": []}
    for plate in ctx.plates:
        neighbours, inputs, band = plate_boundary(world, ctx, plate)
        if not band.any():
            continue
        _convergent, divergent, _transform, _contested = torque.classify_boundary_nodes(plate, neighbours, inputs, ctx.reach, ctx.cache)
        subducting = torque.subducting_boundary_mask(plate, inputs, ctx.reach)
        out["all"].append(inputs.own_points[band])
        out["subducting"].append(inputs.own_points[subducting & band])
        out["ridge"].append(inputs.own_points[divergent & band])
    return {k: (np.concatenate(v) if v else np.zeros((0, 3))) for k, v in out.items()}


def run_one(args: argparse.Namespace) -> None:
    from app import mantle

    setting, world = build_world(args)
    # D is always measured against the ORIGINAL centres at `--reference-scale` (the baseline arm's flow x4), whatever flow
    # the run actually uses, so an intervention that turns the flow off or flips it is scored against where the
    # downwellings were: if trenches follow the flow, the tendency should vanish (off) or reverse (flipped).
    centers = reference_centers(world, args.reference_scale)
    scale_flow(world, args.flow_scale)
    rng = np.random.default_rng(args.seed)
    reference = mantle.MANTLE_FLOW_REFERENCE_RATE * args.reference_scale
    null_sets = [mantle.generate_convection_centers(rng, len(centers), reference) for _ in range(NULL_SETS)]
    samples = []

    def take(myr: float) -> None:
        nodes = band_nodes(world)
        row = {"myr": myr, "centres": [{"strength": float(c.strength / reference), "falloff": float(c.falloff)} for c in centers],
               "nodes": {k: len(v) for k, v in nodes.items()}}
        for k, pts in nodes.items():
            row[f"mean_D_{k}"] = float(downwelling_potential(pts, centers).mean() / reference) if len(pts) else None
        if len(nodes["subducting"]):
            null = np.array([downwelling_potential(nodes["subducting"], s).mean() / reference for s in null_sets])
            row["null_mean_D_subducting"] = {"mean": float(null.mean()), "sd": float(null.std())}
            row["percentile_subducting"] = float((null < row["mean_D_subducting"]).mean())
        if len(nodes["ridge"]):
            null = np.array([downwelling_potential(nodes["ridge"], s).mean() / reference for s in null_sets])
            row["percentile_ridge_in_upwelling"] = float((null > row["mean_D_ridge"]).mean())
        samples.append(row)
        print(json.dumps({"myr": myr, "subducting_nodes": row["nodes"]["subducting"]}), flush=True)

    step_and_sample(world, setting, args, take)
    args.out.write_text(json.dumps({"seed": args.seed, "setting": args.setting, "flow_scale": args.flow_scale, "reference_scale": args.reference_scale, "samples": samples}) + "\n")


def aggregate(args: argparse.Namespace) -> None:
    rows = [json.loads(p.read_text()) for p in sorted(args.runs.glob("*.json"))]
    result = {"seeds": len(rows), "by_time": {}}
    for t in sorted({s["myr"] for r in rows for s in r["samples"]}):
        sel = [s for r in rows for s in r["samples"] if s["myr"] == t]
        pct = [s["percentile_subducting"] for s in sel if "percentile_subducting" in s]
        ridge = [s["percentile_ridge_in_upwelling"] for s in sel if "percentile_ridge_in_upwelling" in s]
        z = [(s["mean_D_subducting"] - s["null_mean_D_subducting"]["mean"]) / s["null_mean_D_subducting"]["sd"] for s in sel
             if "null_mean_D_subducting" in s and s["null_mean_D_subducting"]["sd"] > 0]
        result["by_time"][str(t)] = {
            "seeds_with_subduction": len(pct),
            "subducting_mean_D": float(np.mean([s["mean_D_subducting"] for s in sel if s.get("mean_D_subducting") is not None])) if pct else None,
            "all_band_mean_D": float(np.mean([s["mean_D_all"] for s in sel if s.get("mean_D_all") is not None])),
            "percentile_subducting_median": float(np.median(pct)) if pct else None,
            "percentile_subducting_share_above_0.9": float(np.mean(np.array(pct) > 0.9)) if pct else None,
            "z_subducting_mean": float(np.mean(z)) if z else None,
            "percentile_ridge_in_upwelling_median": float(np.median(ridge)) if ridge else None,
        }
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    for t, b in result["by_time"].items():
        print(t, b)


def main() -> None:
    args = build_parser(__doc__, __file__, run_one, forward=("flow_scale", "reference_scale"), aggregate=aggregate, reference_scale=True).parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()

"""RFC 0018 D1's line check on the development worlds (#888), before any stepping search.

Usage (repo root, this worktree's venv; follow the compute-budget skill, one heavy run at a time, about 40 min):
    backend/.venv/bin/python bin/debug/check_rift_lines_world_d1.py --out ~/mantle-bloom-captures/issue888/line-check \
        [--seeds 106533960 ...] [--half-step-seeds 106533960 1097624328 1097710202] [--years 50]

For each development world (production ``generate_world`` defaults, plates-only), computes the rift line of every
plate with at least ``SPLIT_MIN_NODES`` cells (scaled by node density) at 0 Myr and after ``--years`` Myr stepped
by the production splitter at 1 Myr, without straining anything (step 3 only; the age gate is not applied because
the line is only a placement). The same 50 Myr is repeated from 0.5 Myr steps on three seeds. Per plate it records
the smaller side's area share, the line length over the straight distance between its ends, the craton and oceanic
shares of the line, ``Λ``, and whether cutting along it gives a valid two-daughter split (``a_min`` and
``a_stray`` rules). Writes ``<out>/<seed>-<step>.json`` per world (existing files are kept, so a run resumes) and
``summary.json``.

Pass: at least 90% of lines give a valid split, and the line crosses less craton than the straight chord between its
ends on at least two thirds of the plates that have craton on the line or on the chord.
A held-out or reserve seed from ``analysis/issue797/g2-seed-preregistration.json`` is refused.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import rift_damage_mechanism as mech  # noqa: E402
import rift_damage_step as st  # noqa: E402
from run_rift_damage_d1 import forbidden_seeds  # noqa: E402

DEV_SEEDS = [106533960, 1097624328, 1097710202, 1397817354, 1439649077, 146046339, 1629738962, 1630166320, 1643146570, 1854897345,
             487506346, 494568917, 591453549, 600370317, 604997962, 704718994, 763144949, 870996356, 954941394, 969845080]
HALF_STEP_SEEDS = DEV_SEEDS[:3]


def eligible_rows(world, params) -> list[dict]:
    from app import merge_split

    min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
    rows = []
    for plate in world.plates:
        if plate.node_count() >= min_nodes:
            row = st.plate_line(world, plate, params)
            row["age_myr"] = float(plate.age_years / 1e6)
            rows.append(row)
    return rows


def run_world(seed: int, years: float, step_myr: float, out: Path, params) -> dict:
    target = out / f"{seed}-{step_myr:g}.json"
    if target.exists():
        return json.loads(target.read_text())
    from app import world as world_module

    started = time.perf_counter()
    world = world_module.generate_world(seed)
    world.simulate_climate_biomes = False
    result = {"seed": seed, "step_myr": step_myr, "t0": eligible_rows(world, params), "plates_t0": len(world.plates)}
    for _ in range(round(years / step_myr)):
        world_module.step_world(world, step_myr * 1e6)
    result["t_end"] = eligible_rows(world, params)
    result.update(plates_end=len(world.plates), years=years, wall_s=time.perf_counter() - started)
    target.write_text(json.dumps(result, indent=1, default=float) + "\n")
    print(seed, step_myr, {k: len(v) for k, v in result.items() if isinstance(v, list)}, round(result["wall_s"]), flush=True)
    return result


def aggregate(rows: list[dict]) -> dict:
    with_line = [r for r in rows if r["line"]]
    valid = [r for r in with_line if r["valid"]]
    craton = [r for r in with_line if r["craton_share"] > 0 or r["chord_craton_share"] > 0]
    less = [r for r in craton if r["craton_share"] < r["chord_craton_share"]]
    small = [r["smaller_share"] for r in valid]
    ratio = [r["length_ratio"] for r in with_line]
    reasons: dict = {}
    for r in with_line:
        if not r["valid"]:
            reasons[r["reject"]] = reasons.get(r["reject"], 0) + 1
    return {
        "eligible_plates": len(rows), "with_line": len(with_line), "no_line": len(rows) - len(with_line), "valid": len(valid),
        "valid_share_of_eligible": len(valid) / max(len(rows), 1), "valid_share_of_lines": len(valid) / max(len(with_line), 1),
        "reject_reasons": reasons, "craton_plates": len(craton), "line_less_craton": len(less),
        "line_less_craton_share": len(less) / max(len(craton), 1) if craton else None,
        "smaller_share_median": float(np.median(small)) if small else None,
        "smaller_share_quartiles": [float(x) for x in np.percentile(small, [25, 75])] if small else None,
        "length_ratio_median": float(np.median(ratio)) if ratio else None,
        "oceanic_share_median": float(np.median([r["oceanic_share"] for r in with_line])) if with_line else None,
        "lam_median": float(np.median([r["lam"] for r in with_line])) if with_line else None,
        "open_lines": sum(r["open"] for r in with_line),
        "weak_share_median": float(np.median([r["weak_share"] for r in with_line])) if with_line else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=DEV_SEEDS)
    parser.add_argument("--half-step-seeds", type=int, nargs="*", default=HALF_STEP_SEEDS)
    parser.add_argument("--years", type=float, default=50.0)
    args = parser.parse_args()
    banned = forbidden_seeds()
    bad = [s for s in args.seeds + args.half_step_seeds if s in banned]
    if bad:
        sys.exit(f"seeds {bad} are held-out or reserve seeds; D1 must not use them")
    args.out.mkdir(parents=True, exist_ok=True)
    mech.register_fields()
    params = mech.Params()
    runs = [(s, 1.0) for s in args.seeds] + [(s, 0.5) for s in args.half_step_seeds]
    results = [run_world(seed, args.years, step, args.out, params) for seed, step in runs]
    summary = {"params": params.__dict__, "years": args.years}
    for key in ("t0", "t_end"):
        for step in (1.0, 0.5):
            rows = [r for res in results if res["step_myr"] == step for r in res[key]]
            if rows:
                summary[f"{key}_step{step:g}"] = aggregate(rows)
    both = [r for res in results if res["step_myr"] == 1.0 for key in ("t0", "t_end") for r in res[key]]
    summary["pooled_step1"] = aggregate(both)
    pooled = summary["pooled_step1"]
    summary["pass"] = {"valid_share": pooled["valid_share_of_eligible"] >= 0.9,
                       "craton_share": pooled["line_less_craton_share"] is not None and pooled["line_less_craton_share"] >= 2 / 3}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    print(json.dumps(summary["pooled_step1"], indent=1, default=float))
    print("pass", summary["pass"])


if __name__ == "__main__":
    main()

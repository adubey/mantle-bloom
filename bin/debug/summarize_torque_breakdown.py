"""Pool `plate_torque_breakdown.py` outputs into the tables of analysis/issue942 (#942).

Usage (repo root):
    backend/.venv/bin/python bin/debug/summarize_torque_breakdown.py ~/mantle-bloom-captures/issue942/breakdown \
        [--out analysis/issue942/torque-breakdown.json]

For each world and for both pooled, over every plate-sample: area-weighted crust speed; each driving term's share of
the summed driving torque and each resisting term's share of the summed resisting torque (components along each
plate's rotation axis, so a term that opposes motion on one plate and helps it on another is counted on the side it
acts on); the share of plate area on which each resisting term is the largest; and the share of closing down-going
and of opening boundary length that passes the 1 cm/yr motion gate. Each is given for all plates and split at an
`|omega| R` of 1 cm/yr. Where the runs carry them, it also gives the area-weighted speeds each plate would have
alone (mantle flow only; all own loads with no contact), the mean prescribed mantle-flow speed under the plates,
and the slab terms' work budget (`slab`).
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

SLOW_CM_PER_YR = 1.0


def _summarize(rows: list[dict]) -> dict:
    drive: dict[str, float] = defaultdict(float)
    resist: dict[str, float] = defaultdict(float)
    largest: dict[str, float] = defaultdict(float)
    km: dict[str, float] = defaultdict(float)
    moved = {"continental": [0.0, 0.0], "oceanic": [0.0, 0.0]}
    alone: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0])
    slab: dict[str, float] = defaultdict(float)
    area = pull = pull_ungated = 0.0
    for row in rows:
        area += row["area_km2"]
        cont = row["continental_area_km2"]
        for kind, kind_area in (("continental", cont), ("oceanic", row["area_km2"] - cont)):
            speed = row["speed_cm_per_yr"][kind]
            if speed is not None:
                moved[kind][0] += speed * kind_area
                moved[kind][1] += kind_area
        for name, speeds in row.get("counterfactual_speed_cm_per_yr", {}).items():
            for kind, kind_area in (("all", row["area_km2"]), ("continental", cont), ("oceanic", row["area_km2"] - cont)):
                if speeds[kind] is not None:
                    alone[f"{name}:{kind}"][0] += speeds[kind] * kind_area
                    alone[f"{name}:{kind}"][1] += kind_area
        if "mantle_flow_under_plate_cm_per_yr" in row:
            alone["mantle_flow_under_plate:all"][0] += row["mantle_flow_under_plate_cm_per_yr"] * row["area_km2"]
            alone["mantle_flow_under_plate:all"][1] += row["area_km2"]
        for name, value in row.get("slab", {}).items():
            if name == "hm_km_length_weighted":
                slab["hm_km_times_km"] += value * row["slab"]["subducting_km"]
            elif name == "pull_torque_surviving_cancellation":
                slab["pull_torque_surviving_nm"] += (value or 0.0) * row["slab"]["pull_moment_sum_nm"]
            else:
                slab[name] += value
        for name, value in row["boundary_km"].items():
            km[name] += value
        terms = row.get("torque_along_axis_nm")
        if not terms:
            continue
        for name, value in terms.items():
            (drive if value > 0 else resist)[name] += abs(value)
        resisting = {name: -value for name, value in terms.items() if value < 0}
        if resisting:
            largest[max(resisting, key=resisting.get)] += row["area_km2"]
        pull += row["torque_norm_nm"]["slab_pull"]
        pull_ungated += row["torque_norm_nm"]["slab_pull_ungated"]

    def shares(values: dict[str, float]) -> dict[str, float]:
        total = sum(values.values())
        return {name: value / total for name, value in sorted(values.items(), key=lambda kv: -kv[1])} if total else {}

    return {
        "plate_samples": len(rows),
        "plate_area_share": None,
        "speed_cm_per_yr": {kind: v[0] / v[1] if v[1] else None for kind, v in moved.items()},
        "speed_alone_cm_per_yr": {name: v[0] / v[1] for name, v in sorted(alone.items()) if v[1]},
        "slab": _slab_summary(slab),
        "driving_share": shares(drive),
        "resisting_share": shares(resist),
        "largest_resisting_term_area_share": shares(largest),
        "closing_down_going_length_past_gate": km["down_going_closing_past_gate"] / km["down_going_closing_any_rate"] if km["down_going_closing_any_rate"] else None,
        "opening_length_past_gate": km["opening_past_gate"] / km["opening_any_rate"] if km["opening_any_rate"] else None,
        "closing_down_going_share_of_band": km["down_going_closing_any_rate"] / km["band_total"] if km["band_total"] else None,
        "slab_pull_ungated_over_applied": pull_ungated / pull if pull else None,
        "_area": area,
    }


def _slab_summary(slab: dict[str, float]) -> dict:
    """Pooled slab terms: powers as shares of the work slab pull does where plates move toward their trench."""
    km = slab["subducting_km"]
    work = slab["pull_power_toward_trench_w"]
    if not km or not work:
        return {}
    return {
        "subducting_km_per_sample_sum": km,
        "pull_n_per_m": slab["pull_n"] / (km * 1e3),
        "hm_km": slab["hm_km_times_km"] / km,
        "length_share_no_pull": slab["no_pull_km"] / km,
        "length_share_moving_away_from_trench": slab["moving_away_from_trench_km"] / km,
        "length_share_past_terminal_speed": slab["past_terminal_speed_km"] / km,
        "pull_torque_surviving_cancellation": slab["pull_torque_surviving_nm"] / slab["pull_moment_sum_nm"],
        "power_over_pull_toward_trench": {
            "pull_away_from_trench": slab["pull_power_away_from_trench_w"] / work,
            "drag_trench_normal": -slab["drag_power_trench_normal_w"] / work,
            "drag_trench_parallel": -slab["drag_power_trench_parallel_w"] / work,
            "net": (work + slab["pull_power_away_from_trench_w"] - slab["drag_power_trench_normal_w"] - slab["drag_power_trench_parallel_w"]) / work,
        },
    }


def summarize(rows: list[dict]) -> dict:
    groups = {
        "all": rows,
        "slow": [r for r in rows if r["rate_cm_per_yr"] < SLOW_CM_PER_YR],
        "fast": [r for r in rows if r["rate_cm_per_yr"] >= SLOW_CM_PER_YR],
    }
    out = {name: _summarize(group) for name, group in groups.items()}
    total = out["all"]["_area"]
    for item in out.values():
        item["plate_area_share"] = item.pop("_area") / total if total else None
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("runs", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result: dict = {"slow_threshold_cm_per_yr": SLOW_CM_PER_YR, "worlds": {}}
    pooled: list[dict] = []
    for path in sorted(args.runs.glob("*.json")):
        run = json.loads(path.read_text())
        rows = [row for sample in run["samples"] for row in sample["plates"]]
        pooled += rows
        result["worlds"][str(run["seed"])] = {
            "years": run["years"], "every": run["every"], "world_knobs": run["world_knobs"],
            "solver_converged_all": all(s["solver"]["solver_converged"] for s in run["samples"]),
            "blocked_contacts_max": max(s["solver"]["blocked_contacts"] for s in run["samples"]),
            **summarize(rows),
        }
    result["pooled"] = summarize(pooled)
    text = json.dumps(result, indent=1) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text)


if __name__ == "__main__":
    main()

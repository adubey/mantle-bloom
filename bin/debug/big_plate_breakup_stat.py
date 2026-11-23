"""The share of large plates that break up, per arm, from `big_plate_rift_probe.py` outputs (#995).

For every probe sample (every 5 Myr) and every plate holding at least `--big` of the planet's area, the plate "breaks up" if
the run records a split of that plate id within the next `--window` Myr. Reports, per arm (a probe `--setting`):

- the share of big plate-samples that break up within the window (pooled), and the same per seed with a bootstrap
  interval over seeds;
- how many distinct big plates there were and how many ever split while big;
- splits per run per 100 Myr, and the mean parent share of the planet.

    python3 -I bin/debug/big_plate_breakup_stat.py --dir ~/mantle-bloom-captures/issue995/runs --out analysis/issue995/breakup-stat.json \\
        brake-o0.25-c1 brake-o0.25-c1-heal250 brake-o0.25-c1-kcap07
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np


def arm_runs(directory: Path, setting: str) -> list[dict]:
    pattern = re.compile(rf"rift-state-{re.escape(setting)}-(\d+)\.json$")
    return [json.loads(f.read_text()) for f in sorted(directory.iterdir()) if pattern.match(f.name)]


def run_stats(run: dict, big: float, window: float) -> dict:
    split_times: dict[int, list[float]] = {}
    for s in run["splits"]:
        split_times.setdefault(s["plate"], []).append(s["myr"])
    hit = total = 0
    seen: set[int] = set()
    broke: set[int] = set()
    for sample in run["samples"]:
        t = sample["myr"]
        if t + window > run["samples"][-1]["myr"]:
            continue  # the window must fit in the run
        for p in sample["plates"]:
            if p["area_share"] < big:
                continue
            total += 1
            seen.add(p["plate"])
            if any(t < x <= t + window for x in split_times.get(p["plate"], [])):
                hit += 1
                broke.add(p["plate"])
    parents = [s["parent_area_m2"] / s["total_area_m2"] for s in run["splits"]]
    return {"seed": run["seed"], "big_samples": total, "big_breakups": hit, "big_plates": len(seen), "big_plates_split": len(broke),
            "splits": len(run["splits"]), "mean_parent_share": float(np.mean(parents)) if parents else 0.0}


def big_samples(run: dict, big: float, window: float) -> list[dict]:
    """One record per big plate-sample: craton share, Myr the plate has been continuously big, and whether it split in the window.

    Supercontinent-like plates (high craton share, long spell above `big`) are the Earth analogue that does break up (Pangea
    after ~100-150 Myr, Rodinia and Nuna longer), so the pooled share mixes plates that should and should not split (#995).
    """
    split_times: dict[int, list[float]] = {}
    for s in run["splits"]:
        split_times.setdefault(s["plate"], []).append(s["myr"])
    spell_start: dict[int, float] = {}
    rows = []
    end = run["samples"][-1]["myr"]
    for sample in run["samples"]:
        t = sample["myr"]
        shares = {p["plate"]: p for p in sample["plates"]}
        for pid in list(spell_start):
            if pid not in shares or shares[pid]["area_share"] < big:
                del spell_start[pid]
        for pid, p in shares.items():
            if p["area_share"] >= big:
                spell_start.setdefault(pid, t)
                if t + window <= end:
                    hit = any(t < x <= t + window for x in split_times.get(pid, []))
                    rows.append({"seed": run["seed"], "craton": p["craton_share"], "spell": t - spell_start[pid], "hit": hit})
    return rows


def stratify(rows: list[dict], key: str, edges: list[float], rng: np.random.Generator) -> list[dict]:
    seeds = sorted({r["seed"] for r in rows})
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = [r for r in rows if lo <= r[key] < hi]
        if not sel:
            continue
        per_seed = {s: [r["hit"] for r in sel if r["seed"] == s] for s in seeds}
        boots = []
        for _ in range(2000):
            pick = rng.choice(seeds, len(seeds))
            h = [x for s in pick for x in per_seed[s]]
            boots.append(np.mean(h) if h else np.nan)
        ci = np.nanquantile(boots, [0.025, 0.975])
        p = float(np.mean([r["hit"] for r in sel]))
        out.append({"range": [lo, hi], "samples": len(sel), "break_up_pct": round(100 * p, 1),
                    "bootstrap_95_pct": [round(100 * float(ci[0]), 1), round(100 * float(ci[1]), 1)],
                    "hazard_per_myr_pct": round(100 * -np.log(max(1 - p, 1e-9)) / 20.0, 2)})
    return out


def parse_arm_args() -> argparse.Namespace:
    """The arguments of the arm scripts: `--dir` of probe outputs, `--out`, the big-plate share, the window and the arms."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--big", type=float, default=0.10)
    ap.add_argument("--window", type=float, default=20.0)
    ap.add_argument("settings", nargs="+")
    return ap.parse_args()


def main() -> None:
    args = parse_arm_args()
    rng = np.random.default_rng(995)
    out = {}
    for setting in args.settings:
        rows = [run_stats(r, args.big, args.window) for r in arm_runs(args.dir, setting)]
        if not rows:
            continue
        hit = sum(r["big_breakups"] for r in rows)
        total = sum(r["big_samples"] for r in rows)
        boots = []
        for _ in range(2000):
            pick = rng.integers(0, len(rows), len(rows))
            h = sum(rows[i]["big_breakups"] for i in pick)
            n = sum(rows[i]["big_samples"] for i in pick)
            boots.append(h / n if n else np.nan)
        lo, hi = np.nanquantile(boots, [0.025, 0.975])
        out[setting] = {
            "runs": len(rows), "big_plate_samples": total, "break_up_within_window_pct": round(100 * hit / total, 1) if total else None,
            "bootstrap_95_pct": [round(100 * float(lo), 1), round(100 * float(hi), 1)],
            "big_plates": sum(r["big_plates"] for r in rows), "big_plates_that_split": sum(r["big_plates_split"] for r in rows),
            "splits_per_run_per_100myr": round(float(np.mean([r["splits"] for r in rows])) / 2.0, 1),
            "mean_parent_share": round(float(np.mean([r["mean_parent_share"] for r in rows])), 3),
        }
        recs = [x for r in arm_runs(args.dir, setting) for x in big_samples(r, args.big, args.window)]
        cr = np.quantile([x["craton"] for x in recs], [1 / 3, 2 / 3]) if recs else []
        out[setting]["by_craton_share_tercile"] = stratify(recs, "craton", [-1, float(cr[0]), float(cr[1]), 2], rng) if recs else []
        out[setting]["by_myr_spent_big"] = stratify(recs, "spell", [0, 25, 50, 100, 1e9], rng)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

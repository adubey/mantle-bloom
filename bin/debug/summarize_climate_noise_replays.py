#!/usr/bin/env python3
"""Issue #741: compare drifting and redrawn climate noise across `measure_climate_fields.py` runs.

Reads `<dir>/<save>-<arm><draw>.json` (arm D drift, R redraw) and applies the decision rule in
`analysis/issue741/README.md`. Each run gives one number per metric, its mean over checkpoints
after the opening, so the draw is the statistical unit. It reports each metric's gap between
the arms with a 95% Welch interval, and each drift draw's gap in three time windows.

    backend/.venv/bin/python bin/debug/summarize_climate_noise_replays.py <dir> --out analysis/issue741/climate.json
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

GLOBAL_METRICS = ("precip_land_mm", "precip_ocean_mm", "humidity_land_g_kg", "humidity_ocean_g_kg")
BAND_METRICS = ("precip_land_mm", "precip_ocean_mm", "humidity_land_g_kg", "humidity_ocean_g_kg")
BAND_QUANTILE_METRICS = ("precip_land_quantiles", "precip_ocean_quantiles", "humidity_land_quantiles", "humidity_ocean_quantiles")
EXTRA_METRICS = ("ice_cap_pct_of_land", "ice_volume_km3", "moisture_source_vegetation_g_kg", "moisture_source_water_g_kg")
QUANTILE_NAMES = ("p10", "p50", "p90")
WINDOWS_MYR = ((0.0, 16.0), (16.0, 32.0), (32.0, 48.0))
MIN_LAND_SHARE_PCT = 5.0
UNCHANGED_LIMIT_PCT = 2.0


def checkpoint_metrics(checkpoint: dict) -> dict[str, float | None]:
    """Every metric of one checkpoint, flat."""
    out = {name: checkpoint.get(name) for name in (*GLOBAL_METRICS, *EXTRA_METRICS)}
    for band in checkpoint.get("bands", []):
        for name in BAND_METRICS:
            out[f"band {band['band']} {name}"] = band.get(name)
        for name in BAND_QUANTILE_METRICS:
            values = band.get(name) or [None] * len(QUANTILE_NAMES)
            for label, value in zip(QUANTILE_NAMES, values):
                out[f"band {band['band']} {name.removesuffix('_quantiles')} {label}"] = value
    return out


def run_series(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, float]]:
    """A run's times since its opening, its metric series, and each band's mean land share,
    all without the opening checkpoint."""
    checkpoints = json.loads(path.read_text())["checkpoints"][1:]
    times = np.array([c["elapsed_myr"] for c in checkpoints])
    rows = [checkpoint_metrics(c) for c in checkpoints]
    series = {k: np.array([np.nan if r.get(k) is None else r[k] for r in rows], dtype=float) for k in rows[0]}
    shares: dict[str, list[float]] = defaultdict(list)
    for c in checkpoints:
        for band in c.get("bands", []):
            if "land_area_share_pct" in band:
                shares[band["band"]].append(band["land_area_share_pct"])
    return times - (times[0] - 1.0), series, {k: float(np.mean(v)) for k, v in shares.items()}


def welch_gap(drift: np.ndarray, redraw: np.ndarray) -> dict:
    """The arms' gap as a percentage of the redraw mean, with a 95% Welch interval."""
    n_d, n_r = len(drift), len(redraw)
    mean_r = float(redraw.mean())
    difference = float(drift.mean()) - mean_r
    var_d, var_r = drift.var(ddof=1) / n_d, redraw.var(ddof=1) / n_r
    se = float(np.sqrt(var_d + var_r))
    df = (var_d + var_r) ** 2 / (var_d**2 / (n_d - 1) + var_r**2 / (n_r - 1)) if se > 0.0 else 1.0
    half = float(stats.t.ppf(0.975, df)) * se
    scale = 100.0 / abs(mean_r)
    return {
        "gap_pct": round(difference * scale, 2),
        "low_pct": round((difference - half) * scale, 2),
        "high_pct": round((difference + half) * scale, 2),
        "half_width_pct": round(half * scale, 2),
        "df": round(float(df), 1),
        "drift_draws_pct": [round((float(v) - mean_r) * scale, 2) for v in drift],
        "redraw_draws_pct": [round((float(v) - mean_r) * scale, 2) for v in redraw],
    }


def save_verdict(gap: dict) -> str:
    if gap["low_pct"] >= -UNCHANGED_LIMIT_PCT and gap["high_pct"] <= UNCHANGED_LIMIT_PCT:
        return "unchanged"
    return "excludes 0" if gap["low_pct"] > 0.0 or gap["high_pct"] < 0.0 else "unresolved"


def overall_verdict(gaps: dict[str, dict]) -> str:
    """Unchanged needs every save's interval inside the limit; shifted needs every save's
    interval off zero on one side."""
    verdicts = [save_verdict(g) for g in gaps.values()]
    if all(v == "unchanged" for v in verdicts):
        return "unchanged"
    signs = {np.sign(g["gap_pct"]) for g in gaps.values()}
    return "shifted" if len(gaps) > 1 and len(signs) == 1 and all(g["low_pct"] > 0.0 or g["high_pct"] < 0.0 for g in gaps.values()) else "unresolved"


def window_gaps(times: dict[str, np.ndarray], runs: dict[str, dict[str, np.ndarray]], metric: str) -> dict:
    """Each drift draw's gap to the redraw mean, in each window."""
    out = {}
    redraw = [r for r in runs if r[0] == "R"]
    for lo, hi in WINDOWS_MYR:
        def window_mean(run: str) -> float:
            in_window = (times[run] > lo) & (times[run] <= hi)
            return float(np.nanmean(runs[run][metric][in_window]))

        base = float(np.mean([window_mean(r) for r in redraw]))
        out[f"{int(lo)}-{int(hi)}"] = {
            "drift": {r: round(100.0 * (window_mean(r) / base - 1.0), 2) for r in sorted(runs) if r[0] == "D"},
            "redraw": {r: round(100.0 * (window_mean(r) / base - 1.0), 2) for r in redraw},
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    saves: dict[str, dict] = defaultdict(lambda: {"times": {}, "runs": {}, "shares": {}})
    for path in sorted(args.runs.glob("*.json")):
        match = re.fullmatch(r"(.+)-([DR]\d+)", path.stem)
        if match is None or not path.with_suffix(".done").exists():
            continue
        save, run = match.groups()
        times, series, shares = run_series(path)
        saves[save]["times"][run] = times
        saves[save]["runs"][run] = series
        if shares:
            saves[save]["shares"] = shares

    by_metric: dict[str, dict[str, dict]] = defaultdict(dict)
    report: dict = {"saves": {}, "metrics": {}}
    for save, data in saves.items():
        report["saves"][save] = {
            "runs": sorted(data["runs"]),
            "band_land_share_pct": {k: round(v, 1) for k, v in data["shares"].items()},
            "windows": {m: window_gaps(data["times"], data["runs"], m) for m in ("precip_land_mm", "humidity_land_g_kg", "band 15-30 precip_land_mm")},
        }
        names = next(iter(data["runs"].values())).keys()
        for metric in names:
            drift = np.array([np.nanmean(s[metric]) for r, s in sorted(data["runs"].items()) if r[0] == "D" and metric in s])
            redraw = np.array([np.nanmean(s[metric]) for r, s in sorted(data["runs"].items()) if r[0] == "R" and metric in s])
            drift, redraw = drift[np.isfinite(drift)], redraw[np.isfinite(redraw)]
            if len(drift) >= 2 and len(redraw) >= 2:
                by_metric[metric][save] = welch_gap(drift, redraw)
    for metric, gaps in by_metric.items():
        report["metrics"][metric] = {"verdict": overall_verdict(gaps), **gaps}
    args.out.write_text(json.dumps(report, indent=1) + "\n")
    print(f"{'metric':46s} {'verdict':11s} " + " ".join(f"{s:>30s}" for s in saves))
    for metric, row in report["metrics"].items():
        cells = " ".join(
            f"{row[s]['gap_pct']:+7.2f} [{row[s]['low_pct']:+7.2f},{row[s]['high_pct']:+7.2f}]    " if s in row else " " * 30 for s in saves
        )
        print(f"{metric:46s} {row['verdict']:11s} {cells}")


if __name__ == "__main__":
    main()

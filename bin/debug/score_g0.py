#!/usr/bin/env python3
"""RFC 0015 gate G0 (#554): score `replay_g0.py` outputs against the gate's metrics.

For each arm, given its 0.25 Myr run and its 1 Myr run on one save:
- *swing*: the largest step-to-step change in global ice volume, |dV| / V, after the first
  `--spin-up-myr` (2 by default), on each run;
- *step gap*: |V(0.25 Myr steps) - V(1 Myr steps)| / V(0.25 Myr steps) at the shared end.

With ice ages on, it also reports, on the 0.25 Myr run:
- the second forcing cycle's range of V over the first's (pass < 1.25);
- V's lag behind the intensity curve at its maximum, per cycle (pass <= 1 Myr).

Without ice ages an arm passes at swing < 25% (both runs) and step gap < 10%; with ice ages,
at step gap < 10% and the two forcing conditions.

    python bin/debug/score_g0.py --fine <g0-a-0.25.json> --coarse <g0-a-1.json> [--out score.json]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

SWING_LIMIT = 0.25
GAP_LIMIT = 0.10
CYCLE_GROWTH_LIMIT = 1.25
LAG_LIMIT_MYR = 1.0


def _series(path: Path) -> tuple[list[float], list[float], list[float], dict]:
    run = json.loads(path.read_text())
    cps = run["checkpoints"]
    return [c["elapsed_myr"] for c in cps], [c["ice_volume_km3"] for c in cps], [c.get("ice_age_intensity", 0.0) for c in cps], run


def swing(t: list[float], v: list[float], spin_up_myr: float) -> float | None:
    start = t[0] + spin_up_myr
    worst = None
    for i in range(1, len(v)):
        if t[i - 1] >= start - 1e-9 and v[i - 1] > 0.0:
            change = abs(v[i] - v[i - 1]) / v[i - 1]
            worst = change if worst is None else max(worst, change)
    return worst


def forcing(t: list[float], v: list[float], intensity: list[float], period_myr: float, spin_up_myr: float) -> dict:
    """Range of V and its lag behind the intensity maximum, per forcing period counted from the
    end of spin-up (the RFC's two periods plus 2 Myr)."""
    start = t[0] + spin_up_myr
    cycles = []
    lo = start
    while True:
        hi = lo + period_myr
        idx = [i for i, ti in enumerate(t) if lo - 1e-9 <= ti <= hi + 1e-9]
        if not idx or t[idx[-1]] < hi - 1e-9:
            break
        vs = [v[i] for i in idx]
        peak_v = idx[vs.index(max(vs))]
        peak_f = idx[[intensity[i] for i in idx].index(max(intensity[i] for i in idx))]
        cycles.append({"from_myr": round(lo, 3), "range_km3": max(vs) - min(vs), "lag_myr": round(t[peak_v] - t[peak_f], 3)})
        lo = hi
    out = {"cycles": cycles}
    if len(cycles) >= 2 and cycles[0]["range_km3"] > 0.0:
        out["second_over_first"] = cycles[1]["range_km3"] / cycles[0]["range_km3"]
    return out


def score(fine: Path, coarse: Path | None, spin_up_myr: float) -> dict:
    t, v, f, run = _series(fine)
    period = run["config"].get("ice_age_period_years", 0.0) / 1e6
    out = {"arm": run["config"]["arm"], "save": Path(run["save"]).name, "swing_fine": swing(t, v, spin_up_myr)}
    if coarse is not None:
        tc, vc, _, _ = _series(coarse)
        out["swing_coarse"] = swing(tc, vc, spin_up_myr)
        end = min(t[-1], tc[-1])
        vf = next(v[i] for i in range(len(t) - 1, -1, -1) if t[i] <= end + 1e-9)
        vcoarse = next(vc[i] for i in range(len(tc) - 1, -1, -1) if tc[i] <= end + 1e-9)
        out["step_gap"] = abs(vf - vcoarse) / vf if vf > 0.0 else None
        out["end_myr"] = end
    if period > 0.0:
        out["forcing"] = forcing(t, v, f, period, spin_up_myr)
    gap_ok = out.get("step_gap") is not None and out["step_gap"] < GAP_LIMIT
    if period > 0.0:
        fc = out["forcing"]
        growth_ok = fc.get("second_over_first") is not None and fc["second_over_first"] < CYCLE_GROWTH_LIMIT
        lag_ok = bool(fc["cycles"]) and max(abs(c["lag_myr"]) for c in fc["cycles"]) <= LAG_LIMIT_MYR
        out["passes"] = gap_ok and growth_ok and lag_ok
    else:
        swings = [s for s in (out["swing_fine"], out.get("swing_coarse")) if s is not None]
        out["passes"] = gap_ok and bool(swings) and max(swings) < SWING_LIMIT
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fine", type=Path, required=True)
    ap.add_argument("--coarse", type=Path)
    ap.add_argument("--spin-up-myr", type=float, default=2.0)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    result = score(args.fine, args.coarse, args.spin_up_myr)
    text = json.dumps(result, indent=1)
    if args.out:
        args.out.write_text(text)
    print(text)


if __name__ == "__main__":
    main()

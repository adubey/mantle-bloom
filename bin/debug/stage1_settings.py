"""The 15 stage 1 settings of RFC 0018's D1 search (#888): a Latin hypercube over six knobs, fixed by its seed.

Usage:
    backend/.venv/bin/python bin/debug/stage1_settings.py [--out analysis/issue888/stage1-settings.json]

Knobs and ranges (RFC 0018 starting-values table, with ``k_cap`` extended to 0.1-3): ``eps_dot0`` log-uniform
2e-8 to 3e-7 /yr, ``s_warm`` uniform 0 to 0.8, ``k_cap`` log-uniform 0.1 to 3, ``n_plume`` uniform 0 to 2 TN/m,
``lam_grow`` uniform 0.1 to 0.6, ``l_tip_m`` log-uniform 125 to 500 km. The scrambled Latin hypercube uses seed 888 and
15 points, so the list is reproducible. Everything else is fixed at ``rift_damage_mechanism.Params``.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

SEED, POINTS = 888, 15
KNOBS = (("eps_dot0", 2e-8, 3e-7, True), ("s_warm", 0.0, 0.8, False), ("k_cap", 0.1, 3.0, True), ("n_plume", 0.0, 2e12, False),
         ("lam_grow", 0.1, 0.6, False), ("l_tip_m", 125e3, 500e3, True))


def settings() -> list[dict]:
    unit = qmc.LatinHypercube(d=len(KNOBS), seed=SEED).random(POINTS)
    out = []
    for row in unit:
        out.append({name: float(np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))) if log else lo + u * (hi - lo))
                    for u, (name, lo, hi, log) in zip(row, KNOBS)})
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    table = settings()
    for i, row in enumerate(table):
        print(i, {k: float(f"{v:.3g}") for k, v in row.items()})
    if args.out:
        args.out.write_text(json.dumps({"seed": SEED, "points": POINTS, "knobs": [list(k) for k in KNOBS], "settings": table}, indent=1) + "\n")


if __name__ == "__main__":
    main()

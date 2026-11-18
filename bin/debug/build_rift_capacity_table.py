"""Build or verify the capacity table RFC 0018's rift lines interpolate (`backend/app/data/rift_capacity_table.npz`, #846).

Usage (repo root, this worktree's venv):
    backend/.venv/bin/python bin/debug/build_rift_capacity_table.py --verify [--samples 200]
    backend/.venv/bin/python bin/debug/build_rift_capacity_table.py --write

The table holds `N_cap` (N/m), a column's integrated extensional capacity to normal-fault yield, on a grid of crust
kind (continental, oceanic at a ridge, oceanic mature), crust thickness, mantle-lithosphere thickness and Moho
thermal lag. Each entry is one evaluation of the RFC 0010 Gate 2 yield envelope
(`replay_gate2_capacity._column_capacity`: hydrostatic pore pressure, K0 = 0.8, dry olivine mantle; see
`analysis/issue487/`). `--verify` recomputes a random sample of entries and fails if any differs from the shipped
table by more than 1e-9 relative. `--write` rebuilds every entry (about 23,000 columns, several minutes) and
overwrites the file; do that only when the envelope or the grid changes, and say so in the PR, because every
rift-line world's strength moves with it.

The mantle-lithosphere grid is 2 km steps below 24 km and 8 km above: capacity bends sharply where the lithosphere
is thin, and a uniform 8 km grid missed a young-ridge column by 17% (PR #890 review).
"""

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parent), str(HERE.parents[2] / "backend")]

import replay_gate2_capacity as gate2  # noqa: E402

from app import rift_strength  # noqa: E402

HC_GRID = np.arange(2_000.0, 90_001.0, 2_000.0)
HM_GRID = np.concatenate([np.arange(2_000.0, 24_000.0, 2_000.0), np.arange(24_000.0, 200_001.0, 8_000.0)])
LAG_GRID = np.array([0.0, 50.0, 100.0, 200.0, 400.0])
TABLE_PATH = HERE.parents[2] / "backend" / "app" / "data" / "rift_capacity_table.npz"


def entry(kind: int, hc: int, hm: int, lag: int) -> float:
    return gate2._column_capacity(rift_strength.KINDS[kind], float(HC_GRID[hc]), float(HM_GRID[hm]), float(LAG_GRID[lag]))


def verify(samples: int) -> None:
    table, hc, hm, lag = rift_strength.capacity_table()
    for name, shipped, grid in (("hc", hc, HC_GRID), ("hm", hm, HM_GRID), ("lag", lag, LAG_GRID)):
        if not np.array_equal(shipped, grid):
            sys.exit(f"the shipped table's {name} grid differs from this script's")
    rng = np.random.default_rng(846)
    worst = 0.0
    for _ in range(samples):
        index = tuple(int(rng.integers(size)) for size in table.shape)
        expected = entry(*index)
        worst = max(worst, abs(table[index] - expected) / max(abs(expected), 1.0))
    print(f"{samples} entries recomputed; worst relative difference {worst:.2e}")
    if worst > 1e-9:
        sys.exit("the shipped table does not match the Gate 2 envelope")


def write() -> None:
    table = np.zeros((len(rift_strength.KINDS), len(HC_GRID), len(HM_GRID), len(LAG_GRID)))
    for index in np.ndindex(table.shape):
        table[index] = entry(*index)
    np.savez(TABLE_PATH, table=table, hc=HC_GRID, hm=HM_GRID, lag=LAG_GRID)
    print(f"wrote {TABLE_PATH}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--write", action="store_true")
    parser.add_argument("--samples", type=int, default=200)
    args = parser.parse_args()
    if args.write:
        write()
    else:
        verify(args.samples)


if __name__ == "__main__":
    main()

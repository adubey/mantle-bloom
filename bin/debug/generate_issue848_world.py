"""Generate and age the issue #848 profiling world, saving checkpoints on the way.

Usage (repo root, this worktree's venv; follow the compute-budget skill first):
    backend/.venv/bin/python bin/debug/generate_issue848_world.py \
        --out ~/mantle-bloom-captures/issue848 [--seed N] [--years-myr 200] \
        [--step-myr 1] [--save-at 50,100,200]

The world is for profiling RFC 0018's candidate search on mature plates (#848), not for
evaluating the method. It uses a seed outside every registry: derived the registry way as
the first four SHA-256 bytes of ``mantle-bloom/rfc0018/profile/v1/00``, big-endian, masked to
31 bits. It uses the cratons-first layout with relief from a synthetic tectonic history
(#218), both set explicitly because ``synthetic_relief`` defaults to off. Steps default to
1 Myr to match the G2 protocol, so plate sizes are comparable. The production (great-circle)
splitter runs here, so this is a baseline-arm world.

Stepping is plates-only: ``World.simulate_climate_biomes`` is turned off after generation,
which skips climate, erosion, hydrology and resource formation (the first full-default run
took ~385 s per 1 Myr step, ~21 h for 200 Myr). Plate movement, merges and splits still run,
which is what the breakup profile needs. Surface relief therefore stays as plates and
isostasy leave it, with no erosion, so it is a tectonics-only history.

Writes ``world_<myr>myr.mbworld`` per checkpoint, ``steps.json`` (per-step wall seconds,
plate counts) and ``SHA256SUMS`` into ``--out``. Saves stay out of Git.
"""

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import persistence  # noqa: E402
from app import world as world_mod  # noqa: E402


def derived_seed(label: str = "profile", index: int = 0) -> int:
    digest = hashlib.sha256(f"mantle-bloom/rfc0018/{label}/v1/{index:02d}".encode()).digest()
    return int.from_bytes(digest[:4], "big") & 0x7FFFFFFF


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=derived_seed())
    parser.add_argument("--years-myr", type=float, default=200.0)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--save-at", default="50,100,200")
    args = parser.parse_args()
    out = args.out.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    save_at = {float(x) for x in args.save_at.split(",") if x}

    world = world_mod.generate_world(args.seed, cratons_first=True, synthetic_relief=True)
    world.simulate_climate_biomes = False
    steps = []
    n_steps = round(args.years_myr / args.step_myr)
    for i in range(1, n_steps + 1):
        t0 = time.perf_counter()
        world_mod.step_world(world, args.step_myr * 1e6)
        myr = i * args.step_myr
        steps.append({"myr": myr, "wall_s": round(time.perf_counter() - t0, 3), "plates": len(world.plates)})
        if myr in save_at:
            (out / f"world_{round(myr)}myr.mbworld").write_bytes(persistence.save_exact_bytes(world))
        print(json.dumps(steps[-1]), flush=True)
    (out / "steps.json").write_text(json.dumps({"seed": args.seed, "step_myr": args.step_myr, "steps": steps}))
    sums = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}" for p in sorted(out.glob("*.mbworld"))]
    (out / "SHA256SUMS").write_text("\n".join(sums) + "\n")


if __name__ == "__main__":
    main()

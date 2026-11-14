"""Log why the damage-band prototype rejects ligament-passing paths (#894, from #888).

For each rejected path the sorted component sizes (or daughter sizes), plate cells, the minimum daughter size and
the path length are recorded, so a reader can tell stray pieces from edge slivers.

Usage (run from backend/, plates-only, one heavy run at a time):
    python ../bin/debug/measure_split_rejections.py SEED MYR OUT.json
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rift_damage_step as st  # noqa: E402
from app import world as W  # noqa: E402


def main(seed: int, steps: int, out: str) -> None:
    st.REJECT_SIZES = []
    mechanism = st.install()
    world = W.generate_world(seed)
    world.simulate_climate_biomes = False
    for _ in range(steps):
        W.step_world(world, 1e6)
    json.dump({"counters": dict(mechanism.counters), "rejects": st.REJECT_SIZES}, open(out, "w"))
    print(mechanism.counters)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])

"""Freeze RFC 0021's 20 held-out seeds, before any G2 run (#958).

Usage (repo root):
    python3 bin/debug/freeze_rfc0021_heldout.py [--check]

Derivation, as RFC 0021's G0 states it and as RFC 0018's G2 preregistration did: the SHA-256 of the UTF-8 label
``mantle-bloom/rfc0021/v1/heldout/{index:02d}``, first four bytes as a big-endian uint32 masked with ``0x7fffffff``,
for indices 0, 1, ... A value is skipped if it is zero, repeats an earlier one, is any seed in
``analysis/issue797/g2-seed-preregistration.json`` (its held-out, reserve and excluded sets) or is one of the 20
development worlds (``analysis/issue958/development-seeds.json``). Stops at 20 seeds and asserts the held-out and
development sets are disjoint. Writes ``analysis/issue958/heldout-seeds.json``; with ``--check``, verifies that file
instead of writing it.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "analysis/issue958/heldout-seeds.json"
TEMPLATE = "mantle-bloom/rfc0021/v1/heldout/{index:02d}"
COUNT = 20


def seed_for(label: str) -> int:
    return int.from_bytes(hashlib.sha256(label.encode("utf-8")).digest()[:4], "big") & 0x7FFFFFFF


def _ids(value) -> set[int]:
    if isinstance(value, dict):
        return set().union(*(_ids(v) for v in value.values())) if value else set()
    if isinstance(value, list):
        return {int(v) for v in value if isinstance(v, int)}
    return set()


def derive() -> dict:
    prereg = json.loads((REPO / "analysis/issue797/g2-seed-preregistration.json").read_text())
    rfc0018 = _ids(prereg["excluded_sets"]) | _ids(prereg["heldout_seeds"]) | _ids(prereg["reserve_seeds"])
    development = {int(s) for s in json.loads((REPO / "analysis/issue958/development-seeds.json").read_text())["seeds"]}
    seeds, labels, skipped = [], [], []
    index = 0
    while len(seeds) < COUNT:
        label = TEMPLATE.format(index=index)
        seed = seed_for(label)
        reason = ("zero" if seed == 0 else "repeat" if seed in seeds else "rfc0018_preregistration" if seed in rfc0018
                  else "development" if seed in development else None)
        if reason:
            skipped.append({"label": label, "seed": seed, "reason": reason})
        else:
            seeds.append(seed)
            labels.append(label)
        index += 1
    assert not set(seeds) & development, "held-out and development sets overlap"
    assert not set(seeds) & rfc0018, "held-out set overlaps RFC 0018's preregistration"
    assert len(set(seeds)) == COUNT
    return {
        "protocol": "RFC0021-heldout-v1", "status": "frozen_before_any_G2_run",
        "derivation": "SHA-256(UTF-8 label), first 4 bytes big-endian uint32 masked with 0x7fffffff; skip zero, repeats, "
                      "every seed in analysis/issue797/g2-seed-preregistration.json and the development seeds",
        "label_template": TEMPLATE, "heldout_seeds": seeds, "heldout_labels": labels, "skipped": skipped,
        "development_overlap_count": 0, "rfc0018_overlap_count": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = derive()
    if args.check:
        frozen = json.loads(OUT.read_text())
        if frozen["heldout_seeds"] != result["heldout_seeds"]:
            sys.exit("the frozen held-out set does not match its derivation")
        print("held-out set matches its derivation")
        return
    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result["heldout_seeds"]))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify the Gate 3 finite-strain constitutive mechanics without running its sweep.

Usage:
    backend/.venv/bin/python bin/debug/verify_gate3_finite_strain.py \
        --output analysis/issue487/gate3-finite-strain-verification-review.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from prototype_gate3_finite_strain import verify_finite_strain_mechanics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "title": "RFC 0010 Gate 3 finite-strain mechanics verification",
        "scope": (
            "Elastic stress objectivity is checked with plastic deformation held at identity. "
            "Plastic-flow update objectivity under finite rotations is not checked."
        ),
        "checks": verify_finite_strain_mechanics(),
    }
    encoded = json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")
    if not result["checks"]["passed"]:
        raise SystemExit("finite-strain mechanics verification failed")


if __name__ == "__main__":
    main()

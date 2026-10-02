#!/usr/bin/env python3
"""Issue #248: inventory the save-format versions and persistent state of `.mbworld` files.

For each save, records what the compatibility path has to handle, without changing anything:

- the envelope (`persistence.SAVE_FORMAT_VERSION`; 1 = bare pickled `World`) and whether the
  current build loads it, with the error if not;
- which plate/surface classes it pickles (via a recording `Unpickler`);
- `World` dataclass fields absent from the pickle (backfilled on load) and any extra
  attributes the current `World` no longer declares;
- per surface field, how many line plates carried a stored array versus fell back to
  `ElevationLine.__getattr__`'s lazy default (a field added after the save was written);
- line topology counts (plates, lines, nodes, one-node lines) and quad plates, if any.

Writes one JSON document (keys sorted, floats rounded) so reruns diff cleanly.

Usage (from repo root, this repo's venv):
    backend/.venv/bin/python bin/debug/inventory_legacy_saves.py ~/Downloads/*.mbworld \\
        --out analysis/issue248/save_inventory.json
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import io
import json
import pickle
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

from app import persistence  # noqa: E402
from app.elevation_lines import ElevationLine  # noqa: E402
from app.plates import PlateWithLines  # noqa: E402
from app.sparse_quad_patch import PlateWithSparseQuadPatch  # noqa: E402
from app.surface_fields import SURFACE_FIELDS  # noqa: E402
from app.world import World  # noqa: E402


class _RecordingUnpickler(pickle.Unpickler):
    def __init__(self, data: bytes) -> None:
        super().__init__(io.BytesIO(data))
        self.classes: Counter[str] = Counter()

    def find_class(self, module: str, name: str):
        if module.startswith("app."):
            self.classes[f"{module}.{name}"] += 1
        return super().find_class(module, name)


def _raw_world(data: bytes) -> tuple[object, int | None, list[str]]:
    unpickler = _RecordingUnpickler(data)
    payload = unpickler.load()
    version = None
    world = payload
    if isinstance(payload, dict) and payload.get("format") == persistence.SAVE_FORMAT:
        version = payload.get("version")
        world = payload.get("world")
    else:
        version = 1
    return world, version, sorted(unpickler.classes)


def inventory(path: Path) -> dict:
    data = path.read_bytes()
    entry: dict = {"name": path.name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    try:
        raw, version, classes = _raw_world(data)
    except Exception as exc:  # noqa: BLE001 - recorded, not raised: this is an inventory
        entry["raw_unpickle_error"] = f"{type(exc).__name__}: {exc}"
        return entry
    entry["envelope_version"] = version
    entry["classes"] = classes
    if isinstance(raw, World):
        declared = {f.name for f in dataclasses.fields(World)}
        present = set(raw.__dict__)
        entry["world_fields_missing"] = sorted(declared - present)
        entry["world_attrs_undeclared"] = sorted(present - declared)
        entry["elapsed_years"] = float(raw.__dict__.get("elapsed_years", 0.0))
        entry["seed"] = raw.__dict__.get("seed")
        entry["node_density"] = raw.__dict__.get("node_density")
        line_plates = [p for p in raw.plates if isinstance(p, PlateWithLines)]
        quad_plates = [p for p in raw.plates if isinstance(p, PlateWithSparseQuadPatch)]
        entry["plate_classes"] = dict(sorted(Counter(type(p).__name__ for p in raw.plates).items()))
        lines = [line for p in line_plates for line in p.__dict__.get("_lines", [])]
        stored = Counter()
        for line in lines:
            for name in SURFACE_FIELDS:
                if f"_{name}" in line.__dict__:
                    stored[name] += 1
        entry["line_plates"] = len(line_plates)
        entry["quad_plates"] = len(quad_plates)
        entry["lines"] = len(lines)
        entry["line_nodes"] = int(sum(len(line.__dict__["_theta"]) for line in lines))
        entry["one_node_lines"] = int(sum(len(line.__dict__["_theta"]) == 1 for line in lines))
        entry["line_fields_lazily_defaulted"] = {
            name: len(lines) - stored[name] for name in sorted(SURFACE_FIELDS) if stored[name] < len(lines)
        }
        entry["line_attrs_unknown"] = sorted(
            {
                key
                for line in lines
                for key in line.__dict__
                if key not in ("_phi", "_theta") and key[1:] not in SURFACE_FIELDS
            }
        )
        plate_attrs = Counter(key for p in line_plates for key in p.__dict__)
        entry["line_plate_attrs"] = dict(sorted(plate_attrs.items()))
        del raw
    try:
        world = persistence.load_world_bytes(data)
        entry["loads"] = True
        entry["loaded_plates"] = len(world.plates)
    except Exception as exc:  # noqa: BLE001
        entry["loads"] = False
        entry["load_error"] = f"{type(exc).__name__}: {exc}"
    return entry


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("saves", nargs="+", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    entries = []
    for path in sorted(args.saves):
        print(f"{path.name} ...", file=sys.stderr, flush=True)
        entries.append(inventory(path))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"saves": entries}, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

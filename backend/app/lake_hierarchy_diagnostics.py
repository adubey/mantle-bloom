"""Lake-hierarchy depth / catchment-size instrumentation for GitHub issue #144.

Issue #117 measured seed 23097282 @ 79.2 My producing a `lakes.build_lake_hierarchy` merge
forest whose deepest subtree was ~3,500 levels -- a near-linear chain of thousands of tiny
sub-resolution catchments each spilling into the next, rather than siltation collapsing them
into a handful of real basins (see lakes.py's own module docstring for the two-phase
catchment+Kruskal algorithm that produces this tree). #143 closed one of the two literal
"never holds water, gets nothing" silt gaps the issue named (a chronically-frozen catchment),
but left open whether that's enough to actually shrink the cascade, and #144 -- the direct
follow-up -- found no existing tool measures catchment-size or hierarchy-depth distribution at
all, so any future call on `lakes.SILT_ACCUMULATION_COEFFICIENT` or a depression pre-fill pass
would be tuned by feel rather than against real numbers. This module is that instrumentation.

Like `stranded_basins.py`/`plate_diagnostics.py`, it's a pure read of one already-resolved
`lakes.Lake` forest (the same `HydrologyFields.lake_forest` those two modules read), exposed
both as library functions and as the same `python -m app.<name> <save.mbworld> [--json]`
offline-dump CLI shape those two use -- see docs/debugging.md.

Two numbers this reports that neither existing tool does:

- **Hierarchy depth** -- the length of the longest root-to-leaf chain in the merge forest (1
  for a bare leaf with no children). This is exactly the "~3,500 levels" figure issue #117
  measured, and the number any future fix needs to move. Computed iteratively (post-order,
  explicit stack) rather than recursively, matching lakes.py's own `_resolve`/
  `build_lake_hierarchy`/`_catchment_roots` style for the same reason: a real long spill
  cascade runs thousands of levels deep, well past Python's recursion limit.
- **Catchment size** -- `len(Lake.members)` counted only for *leaves* (`children == []`). A
  leaf is one genuine phase-1 catchment; a parent's `members` is just the union of its
  children's, so counting parents too would double-count every node once per merge level it
  sits under. The "thousands of tiny sub-resolution catchments" pathology issue #117 names is
  a statement about leaf sizes specifically, not the merged tree's internal nodes.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from . import lakes

if TYPE_CHECKING:
    from .world import World

# Convention only -- nothing in the engine hardcodes a step size -- but every long-run save
# examined in GitHub issues #119, #120, and #126, and the UI's Step/Play buttons, uses 100 ky,
# so reporting an approximate step count next to raw elapsed years is what makes
# "79.2 My (~792 steps)" legible. Mirrors plate_diagnostics.CONVENTIONAL_YEARS_PER_STEP /
# stranded_basins.CONVENTIONAL_YEARS_PER_STEP.
CONVENTIONAL_YEARS_PER_STEP = 100_000.0

# Leaf catchment-size bucket edges (node count, inclusive lower bound; the last edge is an
# open-ended "this or more" bucket). Log-ish spacing, not linear: a real forest is dominated
# by 1-5 node leaves with a long tail out to a few hundred, so linear buckets would dump
# almost everything into one bin and say nothing about the tail that actually matters for
# "is siltation collapsing the small ones."
_CATCHMENT_SIZE_BUCKET_EDGES = [1, 2, 6, 21, 101, 501]

# Hierarchy-depth bucket edges (chain length, same inclusive/open-ended convention). Issue
# #117's own number (~3,500) sets the top of this range -- 1000+ is deliberately one bucket
# rather than resolving finer, since at that depth the exact count matters far less than
# "is there any chain even close to this."
_DEPTH_BUCKET_EDGES = [1, 2, 6, 21, 101, 501, 1001]


def _bucket_labels(edges: list[int]) -> list[str]:
    labels = [f"{edges[i]}-{edges[i + 1] - 1}" for i in range(len(edges) - 1)]
    labels.append(f"{edges[-1]}+")
    return labels


def _bucket_counts(values: list[int], edges: list[int]) -> list[int]:
    counts = [0] * len(edges)
    for v in values:
        idx = 0
        for i in range(len(edges) - 1, -1, -1):
            if v >= edges[i]:
                idx = i
                break
        counts[idx] += 1
    return counts


def leaf_catchment_sizes(roots: list["lakes.Lake"]) -> list[int]:
    """`len(members)` for every leaf (`children == []`) in the forest, roots and descendants
    alike -- see this module's own docstring for why leaves only. Plain iterative DFS (no
    depth tracking needed here, unlike `hierarchy_depths`), matching `lakes.iter_all_lakes`'s
    own traversal but filtered to leaves inline rather than materializing the full walk first."""
    sizes: list[int] = []
    stack = list(roots)
    while stack:
        node = stack.pop()
        if node.children:
            stack.extend(node.children)
        else:
            sizes.append(int(len(node.members)))
    return sizes


def hierarchy_depths(roots: list["lakes.Lake"]) -> list[int]:
    """One entry per root: the length of the longest root-to-leaf chain under it (1 for a
    bare leaf root, i.e. a catchment that never merged with anything). Iterative post-order --
    push `(node, False)`, then `(node, True)` ahead of its children so children resolve
    first -- the same explicit-stack shape `lakes._resolve`/`build_lake_hierarchy` already use
    for the same reason: a long spill cascade nests thousands of levels deep, well past
    Python's recursion limit."""
    depth: dict[int, int] = {}
    stack: list[tuple[lakes.Lake, bool]] = [(r, False) for r in roots]
    while stack:
        node, ready = stack.pop()
        if ready:
            depth[id(node)] = 1 if not node.children else 1 + max(depth[id(c)] for c in node.children)
        else:
            stack.append((node, True))
            stack.extend((c, False) for c in node.children)
    return [depth[id(r)] for r in roots]


def build_report(world: "World") -> dict:
    """Structured lake-hierarchy diagnostics for `world` -- the payload behind both the text
    dump and `--json`. Reads `world.hydrology_cache.lake_forest` (persisted in the save, same
    field `stranded_basins.py` reads); a world never stepped with climate on has no hydrology
    snapshot and reports an empty forest."""
    fields = getattr(world, "hydrology_cache", None)
    have_hydrology = fields is not None and len(fields.points) > 0
    forest = fields.lake_forest if have_hydrology else []

    sizes = leaf_catchment_sizes(forest)
    depths = hierarchy_depths(forest)

    return {
        "seed": world.seed,
        "elapsed_years": world.elapsed_years,
        "approx_steps": round(world.elapsed_years / CONVENTIONAL_YEARS_PER_STEP),
        "node_density": world.node_density,
        "have_hydrology_snapshot": bool(have_hydrology),
        "num_roots": len(forest),
        "num_leaf_catchments": len(sizes),
        "max_hierarchy_depth": max(depths, default=0),
        "mean_hierarchy_depth": round(sum(depths) / len(depths), 2) if depths else 0.0,
        "leaf_catchment_size_histogram": {
            "labels": _bucket_labels(_CATCHMENT_SIZE_BUCKET_EDGES),
            "counts": _bucket_counts(sizes, _CATCHMENT_SIZE_BUCKET_EDGES),
        },
        "root_depth_histogram": {
            "labels": _bucket_labels(_DEPTH_BUCKET_EDGES),
            "counts": _bucket_counts(depths, _DEPTH_BUCKET_EDGES),
        },
    }


def format_report(report: dict) -> str:
    lines: list[str] = []
    lines.append("mantle-bloom lake-hierarchy diagnostics")
    lines.append(f"  seed:          {report['seed']}")
    lines.append(
        f"  elapsed:       {report['elapsed_years']:,.0f} yr"
        f"  (~{report['approx_steps']:,} steps @ 100 ky)"
    )
    lines.append(f"  node_density:  {report['node_density']}")
    if not report["have_hydrology_snapshot"]:
        lines.append("")
        lines.append("  no hydrology snapshot in this save (never stepped with climate on) -- nothing to report")
        return "\n".join(line.rstrip() for line in lines)

    lines.append(f"  roots:         {report['num_roots']}")
    lines.append(f"  leaf catchments: {report['num_leaf_catchments']}")
    lines.append(
        f"  hierarchy depth: max {report['max_hierarchy_depth']}"
        f"   mean {report['mean_hierarchy_depth']:.2f}"
    )
    lines.append("")

    lines.append("leaf catchment size (node count)")
    hist = report["leaf_catchment_size_histogram"]
    for label, count in zip(hist["labels"], hist["counts"]):
        lines.append(f"  {label:>10}: {count}")
    lines.append("")

    lines.append("root-to-leaf hierarchy depth (chain length)")
    hist = report["root_depth_histogram"]
    for label, count in zip(hist["labels"], hist["counts"]):
        lines.append(f"  {label:>10}: {count}")
    return "\n".join(line.rstrip() for line in lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m app.lake_hierarchy_diagnostics",
        description="Report lake-hierarchy depth / leaf catchment-size distribution for a saved world.",
    )
    parser.add_argument("save", type=Path, help="path to a .mbworld save file")
    parser.add_argument("--json", action="store_true", help="emit the structured report as JSON")
    args = parser.parse_args(argv)

    if not args.save.is_file():
        parser.error(f"no such file: {args.save}")
    from . import persistence  # local: keeps world.py -> lake_hierarchy_diagnostics -> persistence -> world out of the import cycle

    world = persistence.load_world_bytes(args.save.read_bytes())
    report = build_report(world)

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(format_report(report))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

#!/usr/bin/env python3
"""Issue #520: count river lines drawn across lakes, on `main` or on the branch.

A spilling sink's `flow_target` jumps to the far side of its lake's rim, which can be many cells
away. At the loaded save and after each checkpoint, it measures:

- spill hops: flow edges whose endpoints aren't k-NN neighbours, all and from river nodes, with
  their longest length (km);
- the median flow-edge length (km) over every edge;
- River Inspector segments (`group_rivers` networks) and drawn map segments
  (`render_image._rivers_to_draw`, elevation views): how many are longer than 3x that median and
  how many join non-neighbour nodes.

Segments come from `hydrology.flow_segments` when the backend has it, else from each river
node's own `flow_target` edge, as before #520. With `--render-dir` it also writes an elevation
globe centred on `--center LAT LON` at each checkpoint. `--backend` points at the backend to
import, so the same script measures both trees:

    cd backend
    .venv/bin/python ../bin/debug/measure_river_hops.py <local-save> --myr 0.1 --checkpoint-myr 0.1 \\
        --center 43 153 --render-dir <scratch>/branch --out <scratch>/branch.json
    .venv/bin/python ../bin/debug/measure_river_hops.py <local-save> --myr 0.1 --checkpoint-myr 0.1 \\
        --backend <baseline-worktree>/backend --center 43 153 --render-dir <scratch>/main \\
        --out <scratch>/main.json
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from replay_lib import load_replay, parse_replay_args, run_replay


def _add_arguments(ap) -> None:
    ap.add_argument("--center", type=float, nargs=2, default=(0.0, 0.0), metavar=("LAT", "LON"))
    ap.add_argument("--render-dir", type=Path)


def _km(points, a, b, radius_km):
    from app import geometry  # noqa: E402

    return geometry.angular_distance(points[a], points[b]) * radius_km


def _segments(hydrology, fields, sources):
    if hasattr(hydrology, "flow_segments"):
        from_idx, to_idx, _, _ = hydrology.flow_segments(fields, sources)
        return from_idx, to_idx
    sources = sources[fields.flow_target[sources] >= 0]
    return sources, fields.flow_target[sources]


def _segment_counts(fields, from_idx, to_idx, median_km, radius_km) -> dict:
    neighbour = (fields.neighbor_idx[from_idx] == to_idx[:, None]).any(axis=1)
    neighbour |= (fields.neighbor_idx[to_idx] == from_idx[:, None]).any(axis=1)
    km = _km(fields.points, from_idx, to_idx, radius_km)
    return {
        "segments": int(len(from_idx)),
        "over_3x_median": int(np.sum(km > 3 * median_km)),
        "non_neighbour": int(np.sum(~neighbour)),
        "longest_km": float(km.max()) if len(km) else 0.0,
    }


def _view_rotation(lat_deg: float, lon_deg: float) -> np.ndarray:
    """A rotation taking (lat, lon) to the view centre (0, 0)."""
    from app import geometry  # noqa: E402

    center = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    front = geometry.latlon_to_xyz(0.0, 0.0)
    axis = np.cross(center, front)
    if np.linalg.norm(axis) < 1e-12:
        return np.eye(3)
    rotation = geometry.rotation_matrix(axis / np.linalg.norm(axis), float(np.arccos(np.clip(center @ front, -1, 1))))
    assert np.allclose(rotation @ center, front, atol=1e-9)
    return rotation


def main() -> None:
    args = parse_replay_args(__doc__, _add_arguments)
    world, result = load_replay(args, issue=520)
    from app import render_image  # noqa: E402
    from app.elevation_lines import PLANET_RADIUS_KM  # noqa: E402
    from app.hydroclimate import hydrology  # noqa: E402

    result["center"] = list(args.center)
    rotation = _view_rotation(*args.center)
    snapshots = [0]

    def snapshot() -> dict:
        fields = world.hydrology_cache
        snap = {"elapsed_myr": world.elapsed_years / 1e6, "nodes": int(len(fields.points))}
        sources = np.flatnonzero(fields.flow_target >= 0)
        targets = fields.flow_target[sources]
        median_km = float(np.median(_km(fields.points, sources, targets, PLANET_RADIUS_KM)))
        neighbour = (fields.neighbor_idx[sources] == targets[:, None]).any(axis=1)
        hops = sources[~neighbour]
        hop_km = _km(fields.points, hops, fields.flow_target[hops], PLANET_RADIUS_KM)
        snap["median_edge_km"] = median_km
        snap["spill_hops"] = int(len(hops))
        snap["river_spill_hops"] = int(fields.is_river[hops].sum())
        snap["longest_river_hop_km"] = float(hop_km[fields.is_river[hops]].max()) if fields.is_river[hops].any() else 0.0

        rivers = hydrology.group_rivers(fields)
        snap["rivers"] = len(rivers)
        members = np.concatenate([r.member_idx for r in rivers]) if rivers else np.empty(0, dtype=np.int64)
        snap["inspector"] = _segment_counts(fields, *_segments(hydrology, fields, members), median_km, PLANET_RADIUS_KM)
        drawn, _, _ = render_image._rivers_to_draw(world, hide_short_networks=True)
        snap["map"] = _segment_counts(fields, *_segments(hydrology, fields, drawn), median_km, PLANET_RADIUS_KM)

        if args.render_dir is not None:
            args.render_dir.mkdir(parents=True, exist_ok=True)
            png = render_image.render_png(world, "globe", "elevation", 1000, 1000, rotation)
            path = args.render_dir / f"elevation-{snapshots[0]:02d}.png"
            path.write_bytes(png)
            snap["render"] = path.name
        snapshots[0] += 1
        return snap

    run_replay(world, args, result, snapshot, lambda _seconds: None, lambda _steps: {})


if __name__ == "__main__":
    main()

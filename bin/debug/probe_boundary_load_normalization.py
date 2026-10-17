#!/usr/bin/env python3
"""Boundary-load normalization (RFC 0012 G1/#465 and follow-up #468): are slab pull, ridge push, slab
drag and contact drag the same force per metre of boundary whatever the node density, the
boundary band's width, and how many rows sample it? Each case reports load per represented
front-node length and the fraction of the geometric test span represented by that length.

Builds two oceanic plates split along the 0/180 meridian great circle (`debug_worlds`'
nearest-seed partition, seeds at (0, -10) and (0, 10)). The west plate rotates east about
the pole at `--rate-cm-yr`; the east plate is still. So the 0 meridian is a trench the west
plate closes on, and the 180 meridian is a ridge it opens from, each with a normal rate of
v cos(lat). Only nodes within `--lat-max-deg` of the equator count, which keeps the poles,
where the two boundaries meet, out of it. The west plate's Hm is the oceanic reference
except near the ridge (`--ridge-hm-km`), so ridge push has relief to act on.

For each case the west plate's own band (`torque.gather_boundary_force_inputs`,
`classify_boundary_nodes`, `subducting_boundary_mask`, exactly as `shift_plate` reads them)
gives:

- slab pull: summed |force| of `slab_pull_force_n` / trench length, against the column's
  buoyancy x `SUBDUCTION_LSINK_M` (N/m);
- ridge push: the same for `ridge_push_torque`'s per-node force (recomputed here from its
  own formula, which the torque sums), against 0.5 g (rho_a - rho_w) e_r^2 (N/m);
- slab drag: the viscous resistance per metre of trench at the plate's own rate, from
  `slab_drag_coefficient_matrix`'s per-node coefficient, against (mu / d_s) Lsink v (N/m);
- contact drag: the resistance per metre at `CONTACT_RESISTANCE_REF_CM_PER_YR` with jam = 1
  on the trench band (`contact_drag_pa_s_per_m` x node area x v), against
  `CONTACT_RESISTANCE_N_PER_M`.

The geometric comparison span is R x the counted latitude span. Each load is normalized by
its represented front length (the sum of its shared per-front-node weights), and the output
reports what fraction of that geometric span the sampled band represents. Sweeps: node
density (`--densities`), `BOUNDARY_FORCE_REACH_MULTIPLIER` (`--reach-mults`), and repeated
rows (every band node counted `--repeat` times, as an over-sampled boundary would be). Each
density is built at `--offsets` east-seed offsets spread over one spacing, since
how the two lattices line up sets how many rows fall in the band; a level's value is the mean
over those alignments, and the alignment spread is reported alongside. Gate G1 passes a load
if its spread (max/min - 1) across a sweep's levels is <= 10%.

Usage:
    backend/.venv/bin/python bin/debug/probe_boundary_load_normalization.py \\
        --out analysis/rfc0012/load_normalization.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import debug_worlds, lithosphere, mantle, torque  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402

R = lithosphere.PLANET_RADIUS_M


def build_plates(density: float, rate_cm_yr: float, ridge_hm_m: float, east_offset_deg: float = 0.0, seed: int = 0):
    """The two plates. `east_offset_deg` moves only the east seed, so the east lattice shifts
    against the west one and the boundary: it sets how far apart the two plates' nearest
    nodes sit (the "gap"), which a real world's gap filling leaves anywhere up to a spacing."""
    seeds = [(0.0, -10.0), (0.0, 10.0 + east_offset_deg)]
    plates, _ = debug_worlds._build_plates(seed, seeds, ["oceanic", "oceanic"], density)
    west, east = plates
    west.set_omega(np.array([0.0, 0.0, mantle.cm_per_yr_to_rad_per_yr(rate_cm_yr)]))
    east.set_omega(np.zeros(3))
    for plate in plates:
        pts, _ = plate.all_points_and_elevation()
        n = len(pts)
        lon = np.arctan2(pts[:, 1], pts[:, 0])
        # Uniform columns, so each load's ideal per-metre value is one number. Thin mantle
        # lithosphere within 15 deg of the ridge gives ridge push its relief.
        near_ridge = np.abs(np.abs(lon) - np.pi) < np.radians(15.0)
        hm = np.where(near_ridge, ridge_hm_m, lithosphere.REFERENCE_HM_OCEANIC_M)
        plate.set_fields_on_plate(
            crustal_thickness_m=np.full(n, lithosphere.REFERENCE_HC_OCEANIC_M),
            mantle_lithosphere_thickness_m=hm,
            craton_crust_m=np.zeros(n),
        )
    return west, east


def repeat_inputs(inputs: torque.BoundaryForceInputs, k: int) -> torque.BoundaryForceInputs:
    if k == 1:
        return inputs
    fields = {}
    for name in inputs.__dataclass_fields__:
        value = getattr(inputs, name)
        fields[name] = np.concatenate([value] * k, axis=0) if isinstance(value, np.ndarray) else value
    return torque.BoundaryForceInputs(**fields)


class RepeatedPlate:
    """The west plate seen through `repeat_inputs`: the reads the force functions make,
    repeated to match."""

    def __init__(self, plate, k: int):
        self._plate, self._k = plate, k
        self.omega = plate.omega
        self.crust_type = plate.crust_type

    def accounting_areas_m2(self, spacing_rad):
        return np.concatenate([self._plate.accounting_areas_m2(spacing_rad)] * self._k)


def measure(west, east, density: float, reach_mult: float, repeat: int, args) -> dict:
    torque.BOUNDARY_FORCE_REACH_MULTIPLIER = reach_mult
    spacing = line_spacing_rad(density)
    reach = torque.boundary_force_reach_rad(spacing)
    inputs = torque.gather_boundary_force_inputs(west, [east], spacing, reach)
    _conv, divergent, _tr, _contested = torque.classify_boundary_nodes(west, [east], inputs, reach)
    subducting = torque.subducting_boundary_mask(west, inputs, reach)
    lat = np.arcsin(np.clip(inputs.own_points[:, 2], -1.0, 1.0))
    window = np.abs(lat) <= np.radians(args.lat_max_deg)
    length_m = R * 2.0 * np.radians(args.lat_max_deg)
    lon = np.arctan2(inputs.own_points[:, 1], inputs.own_points[:, 0])
    trench_side = np.abs(lon) < np.pi / 2

    inputs_k = repeat_inputs(inputs, repeat)
    plate_k = RepeatedPlate(west, repeat)
    tile = lambda m: np.concatenate([m] * repeat)  # noqa: E731
    sub_k = tile(subducting & window)
    div_k = tile(divergent & window & ~trench_side)
    band_k = tile((inputs.dist_to_neighbor <= reach) & window & trench_side)
    ds = spacing * R
    boundary_lengths = torque.boundary_node_length_m(inputs_k, spacing)

    pull_length_m = boundary_lengths[sub_k].sum()
    pull = torque.slab_pull_force_n(plate_k, inputs_k, sub_k, spacing).sum() / pull_length_m
    buoyancy = torque.column_buoyancy_pa(np.zeros(1), np.array([lithosphere.REFERENCE_HM_OCEANIC_M]), np.zeros(1))[0]
    pull_ideal = buoyancy * torque.SUBDUCTION_LSINK_M

    # ridge_push_torque's own per-node magnitude (it sums these against direction x lever arm).
    rho_c = lithosphere.node_crust_density(inputs_k.own_crust_type_codes[div_k], west.crust_type)
    z = lithosphere.isostatic_elevation(inputs_k.own_hc[div_k], inputs_k.own_hm[div_k], rho_c)
    e_r = np.clip(z - torque.ABYSSAL_PLAIN_REFERENCE_ELEVATION_M, 0.0, None)
    push_scale = 0.5 * lithosphere.GRAVITY_M_S2 * (lithosphere.RHO_ASTHENOSPHERE - lithosphere.RHO_WATER)
    ridge_length_m = boundary_lengths[div_k].sum()
    push = (push_scale * e_r**2 * boundary_lengths[div_k]).sum() / ridge_length_m
    push_ideal = push_scale * float(np.average(e_r**2, weights=boundary_lengths[div_k])) if len(e_r) else 0.0

    # Slab drag per node: (mu / d_s) Lsink ds, times the node's speed, summed per metre.
    v_m_s = np.linalg.norm(np.cross(west.omega, inputs_k.own_points[sub_k]), axis=1) * R / torque.SECONDS_PER_YEAR
    eta_slab = torque.SLAB_MANTLE_VISCOSITY_PA_S / torque.SHEAR_ZONE_THICKNESS_M
    slab_drag = (eta_slab * torque.SUBDUCTION_LSINK_M * boundary_lengths[sub_k] * v_m_s).sum() / pull_length_m
    slab_drag_ideal = (
        eta_slab * torque.SUBDUCTION_LSINK_M * float(np.average(v_m_s, weights=boundary_lengths[sub_k])) if len(v_m_s) else 0.0
    )

    v_ref = torque.CONTACT_RESISTANCE_REF_CM_PER_YR * 0.01 / torque.SECONDS_PER_YEAR
    jam = band_k.astype(float)
    coeff = torque.contact_drag_pa_s_per_m(plate_k, inputs_k, np.zeros(len(jam), dtype=bool), spacing, jam=jam)
    contact = (coeff * plate_k.accounting_areas_m2(spacing) * v_ref).sum() / boundary_lengths[band_k].sum()

    band_rows = (inputs.dist_to_neighbor[window & trench_side] <= reach).sum() * ds / length_m
    near = torque.gather_boundary_force_inputs(west, [east], spacing, np.inf).dist_to_neighbor
    # Edge gap: in each one-spacing latitude bin, the side's nearest node to the other plate;
    # its median over the bins. It varies along the boundary, so a minimum would hide it.
    bins = np.floor(lat / spacing).astype(int)
    gaps = []
    for side in (trench_side, ~trench_side):
        sel = window & side
        edge = [near[sel][bins[sel] == b].min() for b in np.unique(bins[sel])]
        gaps.append(float(np.median(edge) / spacing))
    return {
        "density": density,
        "reach_mult": reach_mult,
        "repeat": repeat,
        "nodes": int(len(inputs.own_points)),
        "band_rows_per_spacing": float(band_rows),
        "trench_gap_spacings": gaps[0],
        "ridge_gap_spacings": gaps[1],
        "slab_pull_ratio": float(pull / pull_ideal),
        "ridge_push_ratio": float(push / push_ideal) if push_ideal else None,
        "slab_drag_ratio": float(slab_drag / slab_drag_ideal) if slab_drag_ideal else None,
        "contact_drag_ratio": float(contact / torque.CONTACT_RESISTANCE_N_PER_M),
        "slab_pull_n_per_m": float(pull),
        "ridge_push_n_per_m": float(push),
        "represented_length_fraction": {
            "slab_pull": float(pull_length_m / length_m),
            "ridge_push": float(ridge_length_m / length_m),
            "contact_drag": float(boundary_lengths[band_k].sum() / length_m),
        },
    }


LOADS = ("slab_pull_ratio", "ridge_push_ratio", "slab_drag_ratio", "contact_drag_ratio")


def spread(values: list[float]) -> float | None:
    values = [v for v in values if v]
    return float(max(values) / min(values) - 1.0) if values else None


def summarize(rows: list[dict]) -> dict:
    """Per load, the mean over alignment offsets, and the spread alignment alone causes."""
    return {
        key: {"mean": float(np.mean([r[key] for r in rows])), "alignment_spread": spread([r[key] for r in rows])}
        for key in LOADS
        if all(r[key] for r in rows)
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--densities", type=float, nargs="+", default=[1.0, 2.0, 4.0])
    parser.add_argument("--reach-mults", type=float, nargs="+", default=[1.5, 3.0, 6.0])
    parser.add_argument("--repeats", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--offsets", type=int, default=6, help="east-seed offsets per density, spread over one spacing")
    parser.add_argument("--rate-cm-yr", type=float, default=5.0)
    parser.add_argument("--ridge-hm-km", type=float, default=20.0)
    parser.add_argument("--lat-max-deg", type=float, default=60.0)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    default_reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER
    cases = {"density": [], "reach": [], "repeat": []}
    for density in args.densities:
        step_deg = np.degrees(line_spacing_rad(density)) / args.offsets
        for i in range(args.offsets):
            west, east = build_plates(density, args.rate_cm_yr, args.ridge_hm_km * 1000.0, i * step_deg)
            cases["density"].append(measure(west, east, density, default_reach, 1, args))
            if density == 1.0:
                cases["reach"] += [measure(west, east, density, m, 1, args) for m in args.reach_mults]
                cases["repeat"] += [measure(west, east, density, default_reach, k, args) for k in args.repeats]
    torque.BOUNDARY_FORCE_REACH_MULTIPLIER = default_reach
    axis = {"density": "density", "reach": "reach_mult", "repeat": "repeat"}

    result = {
        "commit": subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip(),
        "args": {k: v for k, v in vars(args).items() if k != "out"},
        "sweeps": {},
    }
    for name, rows in cases.items():
        levels = sorted({r[axis[name]] for r in rows})
        by_level = {level: summarize([r for r in rows if r[axis[name]] == level]) for level in levels}
        sweep_spread = {
            key: spread([by_level[level][key]["mean"] for level in levels]) for key in LOADS if all(key in by_level[lv] for lv in levels)
        }
        result["sweeps"][name] = {"levels": {str(k): v for k, v in by_level.items()}, "spread": sweep_spread, "rows": rows}
        print(f"\n{name} sweep (ratio to the ideal per-metre load; mean over {args.offsets} alignments, alignment spread)")
        for level in levels:
            cells = "  ".join(
                f"{key.removesuffix('_ratio')}={v['mean']:.2f} ({v['alignment_spread']:.0%})" for key, v in by_level[level].items()
            )
            gaps = [r["trench_gap_spacings"] for r in rows if r[axis[name]] == level]
            print(f"  {axis[name]}={level:g}  median trench gap {min(gaps):.2f}-{max(gaps):.2f} sp  {cells}")
        print("  spread across levels: " + "  ".join(f"{k.removesuffix('_ratio')}={v:.0%}" for k, v in sweep_spread.items()))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()

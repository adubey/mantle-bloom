"""Issue #213 comparison only: the additive fault-relief pass before conservation.

Keep this out of app/. It is the pre-issue implementation, including its summed
elevation requests and Hc backing, so a paired sweep can reproduce that behavior.
"""

from __future__ import annotations

import numpy as np

from app import faults, geometry, lithosphere
from app.elevation_lines import ELEV_CHANGE_MIN_DELTA_M, PLANET_RADIUS_KM


def apply_plate_fault_relief(world, plate, years_myr, _cache=None):
    if not hasattr(plate, "lines"):
        return
    active = [f for f in faults._all_faults(world) if f.plate_id == plate.plate_id and f.active]
    if not active:
        return
    own_points, tree = faults._own_points_and_tree(plate, _cache)
    if len(own_points) == 0:
        return
    rate_scale, reach_scale = faults._relief_mode_scales(world)
    reach_rad = reach_scale * faults.MAX_FAULT_REACH_KM / PLANET_RADIUS_KM

    delta = np.zeros(len(own_points))
    reason = np.zeros(len(own_points), dtype=float)
    for fault in active:
        trace = faults.fault_world_points(fault, plate)
        neighbours = tree.query_ball_point(trace, reach_rad)
        affected = sorted({i for sub in neighbours for i in sub})
        if not affected:
            continue
        affected = np.array(affected)
        pts = own_points[affected]
        d = np.min(np.linalg.norm(pts[:, None, :] - trace[None, :, :], axis=-1), axis=1)
        taper = np.clip(1.0 - d / reach_rad, 0.0, 1.0)
        slip_norm = float(np.clip(fault.slip_rate_m_per_myr / faults.SLIP_RATE_REF_M_PER_MYR, 0.2, 3.0))
        fault_scale = faults.BOUNDARY_FAULT_RELIEF_SCALE if fault.boundary else 1.0
        mag = taper * slip_norm * years_myr * rate_scale * fault_scale

        if fault.kind == faults._KIND_REVERSE:
            contrib = faults.REVERSE_UPLIFT_M_PER_MYR * mag
        elif fault.kind == faults._KIND_NORMAL:
            dip_dir_world = geometry.to_world(plate.frame, fault.dip_dir_local)
            mid = trace[len(trace) // 2]
            hanging = (pts - mid) @ dip_dir_world > 0.0
            contrib = np.where(hanging, -faults.NORMAL_THROW_M_PER_MYR, faults.NORMAL_SHOULDER_UPLIFT_M_PER_MYR) * mag
        else:
            contrib = (faults.STRIKE_SLIP_RIDGE_M_PER_MYR + fault.strike_sense * faults.STRIKE_SLIP_BEND_M_PER_MYR) * mag

        delta[affected] += contrib
        if not (fault.boundary and fault.kind == faults._KIND_REVERSE):
            reason[affected] = faults._KIND_REASON[fault.kind]

    if not np.any(delta):
        return

    new_lines = []
    offset = 0
    changed = False
    for line in plate.lines:
        n = len(line)
        if n == 0:
            new_lines.append(line)
            continue
        seg_delta = delta[offset : offset + n]
        seg_reason = reason[offset : offset + n]
        offset += n
        if not np.any(seg_delta):
            new_lines.append(line)
            continue
        new_hc, new_elev = lithosphere.back_elevation_gain(line, plate, seg_delta, seg_delta != 0.0)
        moved = np.abs(new_elev - line.elevation) >= ELEV_CHANGE_MIN_DELTA_M
        new_reason = np.where(moved & (seg_reason > 0), seg_reason, line.elev_change_reason)
        new_lines.append(line.replace(elevation=new_elev, crustal_thickness_m=new_hc, elev_change_reason=new_reason))
        changed = True
    if changed:
        plate.set_lines(new_lines)

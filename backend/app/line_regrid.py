"""Periodic line regularization.

Per-step boundary evolution (boundary.py) only ever touches the two ends of a line --
inserting at target spacing when growing, deleting when shrinking -- so interior spacing
stays regular on its own. What it can't fix is spacing that's drifted at a *transform*
boundary (nodes sheared along the line without insertion/deletion) or after several steps'
worth of end-growth at a slightly different rate than the line's original spacing. This
module re-derives a fresh evenly-spaced node set spanning each line's *existing* extent
(the two endpoints are preserved exactly -- regularizing never changes where a line's
physical edge is, only how regularly it's sampled) and interpolates elevation onto it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .plates import TARGET_LINE_SPACING_RAD, ElevationLine, Plate

if TYPE_CHECKING:
    from .world import World

REGULARIZE_INTERVAL_STEPS = 5
IRREGULARITY_TOLERANCE = 1.5  # regularize a line if any gap exceeds this multiple of target


def needs_regularizing(line: ElevationLine) -> bool:
    if len(line.theta) < 3:
        return False
    dtheta_target = TARGET_LINE_SPACING_RAD / max(np.cos(line.phi), 1e-3)
    gaps = np.diff(line.theta)
    ratio = gaps / dtheta_target
    return bool(np.any(ratio > IRREGULARITY_TOLERANCE) or np.any(ratio < 1.0 / IRREGULARITY_TOLERANCE))


def regularize_line(line: ElevationLine) -> ElevationLine:
    if len(line.theta) < 3:
        return line

    dtheta_target = TARGET_LINE_SPACING_RAD / max(np.cos(line.phi), 1e-3)
    theta_min, theta_max = line.theta[0], line.theta[-1]
    span = theta_max - theta_min
    n = max(int(round(span / dtheta_target)) + 1, 2)

    new_theta = np.linspace(theta_min, theta_max, n)
    new_elevation = np.interp(new_theta, line.theta, line.elevation)
    # channel_depth/lake_depth/glacier_depth interpolated the same way -- a plain reset to 0
    # here would wipe out a river's carved channel (or a glacier) every time this line's
    # spacing drifts enough to trigger regularizing, which runs periodically throughout the
    # simulation (see REGULARIZE_INTERVAL_STEPS), not as a rare one-off event like a
    # merge/split resample.
    new_channel_depth = np.interp(new_theta, line.theta, line.channel_depth)
    new_lake_depth = np.interp(new_theta, line.theta, line.lake_depth)
    new_glacier_depth = np.interp(new_theta, line.theta, line.glacier_depth)
    # volcano_active_years_remaining interpolates the same way; is_volcano is interpolated as
    # a float (blending a volcano node's 1.0 against a non-volcano neighbor's 0.0) then
    # thresholded back to bool, same spirit as the others -- a resampled node keeps "was this
    # near a volcano" rather than silently losing volcanic provenance every regularize pass.
    new_volcano_active_years_remaining = np.interp(new_theta, line.theta, line.volcano_active_years_remaining)
    new_is_volcano = np.interp(new_theta, line.theta, line.is_volcano.astype(float)) > 0.5
    return ElevationLine(
        phi=line.phi,
        theta=new_theta,
        elevation=new_elevation,
        channel_depth=new_channel_depth,
        lake_depth=new_lake_depth,
        glacier_depth=new_glacier_depth,
        is_volcano=new_is_volcano,
        volcano_active_years_remaining=new_volcano_active_years_remaining,
    )


def regularize_plate_lines(plate: Plate) -> None:
    plate.lines = [
        regularize_line(line) if needs_regularizing(line) else line for line in plate.lines
    ]


def regularize_world_lines(world: "World") -> None:
    for plate in world.plates:
        regularize_plate_lines(plate)

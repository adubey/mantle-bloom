"""Periodic line regularization ("garbage collection").

Per-step boundary evolution (boundary.py) only ever touches the two ends of a line --
inserting at target spacing when growing, deleting when shrinking -- so interior spacing
stays regular on its own. What it can't fix is spacing that's drifted at a *transform*
boundary (nodes sheared along the line without insertion/deletion) or after several steps'
worth of end-growth at a slightly different rate than the line's original spacing. This
module re-derives a fresh evenly-spaced node set spanning each line's *existing* extent
(the two endpoints are preserved exactly -- GC never changes where a line's physical edge
is, only how regularly it's sampled) and interpolates elevation onto it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .plates import TARGET_LINE_SPACING_RAD, ElevationLine, Plate

if TYPE_CHECKING:
    from .world import World

GC_INTERVAL_STEPS = 5
IRREGULARITY_TOLERANCE = 1.5  # trigger GC on a line if any gap exceeds this multiple of target


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
    return ElevationLine(phi=line.phi, theta=new_theta, elevation=new_elevation)


def garbage_collect_plate(plate: Plate) -> None:
    plate.lines = [
        regularize_line(line) if needs_regularizing(line) else line for line in plate.lines
    ]


def garbage_collect_world(world: "World") -> None:
    for plate in world.plates:
        garbage_collect_plate(plate)

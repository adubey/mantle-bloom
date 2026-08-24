"""Shared boundary-motion primitives that survive the move to polygon-based deformation.

Per-step boundary evolution itself now lives on `PlateWithLines.deform` (see plates.py) --
classification is by *geometry* (did a plate's rotated territory end up overlapping a
neighbor's, or open up unclaimed space) rather than by the *velocity* decomposition this
module used to compute directly. What's left here is `closing_rate` (still needed by
merge_split.py to detect a genuine, actively-converging continental collision, as opposed to
two plates that merely tile-neighbor) and the two threshold constants that both
merge_split.py and plates.py's deform() share.
"""

from __future__ import annotations

import numpy as np

from .elevation_lines import TARGET_LINE_SPACING_RAD
from . import mantle

TRANSFORM_RATE_THRESHOLD = mantle.cm_per_yr_to_rad_per_yr(1.0)

# A node this close to (and, in plates.py's deform(), also geometrically overlapping) a
# neighbor's territory counts as consumed/collided rather than merely adjacent. Shared with
# merge_split.py, which reuses it as its own continental-contact distance -- see that
# module's MERGE_CONTACT_DISTANCE_RAD.
MERGE_THRESHOLD_RAD = 0.4 * TARGET_LINE_SPACING_RAD


def closing_rate(
    points: np.ndarray, self_omega: np.ndarray, neighbor_omega: np.ndarray, neighbor_points: np.ndarray
) -> np.ndarray:
    """Positive = this plate's material is moving toward the neighbor's (convergent) at
    this point; negative = moving apart (divergent). Used by merge_split.py to confirm a
    pair of continental plates is actively colliding, not just touching -- two plates are
    already touching along their entire shared boundary by construction (plates.py's tiling
    has no gaps), so proximity alone can't distinguish an actively-colliding pair from any
    other pair of neighbors; it has to check motion."""
    v_self = np.cross(self_omega, points)
    v_neighbor = np.cross(neighbor_omega, points)
    normal_dir = neighbor_points - points
    norm = np.linalg.norm(normal_dir, axis=-1, keepdims=True)
    safe_norm = np.where(norm < 1e-12, 1.0, norm)
    normal_dir = normal_dir / safe_norm
    return np.sum((v_self - v_neighbor) * normal_dir, axis=-1)

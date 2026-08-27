"""Continental shelf range: how far from a coastline "shallow shelf water" extends.

Isostasy (see lithosphere.py) determines submerged continental crust's own depth directly, so
nothing here relaxes elevation toward a shelf/deep-water target any more -- this module now
just holds the shared shelf-range constant geology.py's own oil & gas deposit formation uses
to tell shelf water (shallow, oil/gas-favorable) from open ocean.
"""

from __future__ import annotations

from .elevation_lines import PLANET_RADIUS_KM

# Real continental shelves are shallow and comparatively narrow before the "shelf break"
# drops off toward deep water -- SHELF_RANGE_KM has no exact figure to port (none was given),
# picked as a reasonable shelf width.
SHELF_RANGE_KM = 200.0
SHELF_RANGE_RAD = SHELF_RANGE_KM / PLANET_RADIUS_KM

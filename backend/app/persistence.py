"""Whole-`World` save/load to a single opaque file -- the "File > Save/Load" feature.

Deliberately just `pickle`, not a hand-written interchange format: every field on `World`
(see world.py) -- plates, mantle centers, the collision-progress dict, the climate/hydrology
caches -- is already a plain dataclass or numpy array with no
open handles or unpicklable state, so pickling the object graph directly round-trips it
exactly with no bespoke (de)serialization code to keep in sync as `World`'s own fields
change. This deliberately makes no promise of compatibility across app versions (pickling by
class identity means a later renamed/restructured field breaks old files) -- an accepted
trade for "just get me back what I had," not a stable export format (see geodesic.py's
`export_hexgrid` for that end of the spectrum instead).

Security note: unpickling is equivalent to running arbitrary code from the file's bytes.
Acceptable here because this server is a single-user localhost dev tool already (see
main.py's CORS allowlist) -- the same trust boundary every other route already assumes --
but worth keeping in mind before ever exposing this beyond localhost.
"""

from __future__ import annotations

import pickle

from .world import World


def save_world_bytes(world: World) -> bytes:
    return pickle.dumps(world)


def load_world_bytes(data: bytes) -> World:
    """Raises whatever pickle itself raises on malformed/foreign input (UnpicklingError,
    EOFError, AttributeError for an unknown class, etc.) -- the caller (main.py) is
    responsible for catching broadly and mapping to a 400, not this function."""
    world = pickle.loads(data)
    if not isinstance(world, World):
        raise TypeError(f"expected a World, got {type(world).__name__}")
    _backfill_added_fields(world)
    _drop_derived_caches(world)
    return world


def _backfill_added_fields(world: World) -> None:
    """Give a `default_factory` field added to `World` after this save was written its empty
    default. Unlike a plain-default field (`steps_taken: int = 0`, a class attribute an old
    pickle falls through to), a `field(default_factory=...)` sets no class attribute, so an
    old save's `__dict__` simply won't have the key and the first access `AttributeError`s.
    Only mutable-default fields need listing here."""
    if not hasattr(world, "stranded_basin_tracks"):
        world.stranded_basin_tracks = []
    if not hasattr(world, "overlap_progress"):
        world.overlap_progress = {}
    if not hasattr(world, "faults"):
        world.faults = []
    if not hasattr(world, "boundary_faults"):
        world.boundary_faults = []
    if not hasattr(world, "fault_systems"):
        world.fault_systems = []
    if not hasattr(world, "earthquakes"):
        world.earthquakes = []
    # Eustatic sea level (eustasy.py): a save written before this existed has a fixed
    # sea_level_m and no water budget -- snapshot the budget from that save's own hypsometry
    # + sea level so loading it doesn't jump the shoreline, then let it be conserved onward.
    if getattr(world, "ocean_water_column_m", None) is None:
        from . import eustasy

        eustasy.initialize_water_budget(world)


def _drop_derived_caches(world: World) -> None:
    """Clear every plate's lazily-rebuilt geometry cache (bounding polygon, its k-d tree,
    the contains_batch row lookup). These are pure functions of a plate's current lines +
    frame, so a stale one from an older app version -- e.g. a pre-keyhole `outline_world`
    result, or a `_RowLookup` from before it grew per-arc intervals -- would otherwise be
    trusted as-is on load. Cheap: each rebuilds on first use after this."""
    for plate in world.plates:
        invalidate = getattr(plate, "_invalidate_bounding_polygon", None)
        if callable(invalidate):
            invalidate()
    # The render path's cached node-cloud k-d tree (see World.node_kdtree_cache) -- a pure
    # function of the just-invalidated plate geometry, rebuilt on the first render after load.
    world.node_kdtree_cache = None

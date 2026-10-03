"""Whole-`World` save/load to a single opaque file -- the "File > Save/Load" feature.

Deliberately just `pickle`, not a hand-written interchange format: every field on `World`
(see world.py) -- plates, mantle centers, the collision-progress dict, the climate/hydrology
caches -- is already a plain dataclass or numpy array with no
open handles or unpicklable state, so pickling the object graph directly round-trips it
exactly with no bespoke (de)serialization code to keep in sync as `World`'s own fields
change. Pickling by class identity ties a save to the code, so this is "just get me back what
I had," not a stable export format (see geodesic.py's `export_hexgrid` for that end of the
spectrum instead). What it does promise is in docs/save-compatibility.md: which older save
versions load, how fields added since are backfilled, how a line-backed save moves to sparse
quads (legacy_conversion.py), and that a save this build can't read fails with
`CorruptSaveError` or `UnsupportedSaveVersionError` rather than loading partially.

Security note: unpickling is equivalent to running arbitrary code from the file's bytes.
Acceptable here because this server is a single-user localhost dev tool already (see
main.py's CORS allowlist) -- the same trust boundary every other route already assumes --
but worth keeping in mind before ever exposing this beyond localhost.
"""

from __future__ import annotations

import pickle

from .world import World


# The saved file is a pickled `{"format": SAVE_FORMAT, "version": N, ...}` envelope:
#
# - Version 1 is the original bare pickled `World` (no envelope), still loadable.
# - Version 2 added the envelope `{"format", "version", "world"}` so a save can carry, and a
#   loader can check, which surface representations it may contain -- see
#   sparse_quad_patch.py, whose plates also version their own pickled state.
# - Version 3 adds `"surface"`: `"lines"`, `"quad"`, or `"empty"` (no plates), checked
#   against the plates actually unpickled. A world never mixes the two surfaces.
#
# Bump on any change an older build can't read. docs/save-compatibility.md is the policy for
# line-backed (legacy) saves and how they move to quads (legacy_conversion.py).
SAVE_FORMAT = "mantle-bloom-world"
SAVE_FORMAT_VERSION = 3
SURFACES = ("lines", "quad", "empty")


class SaveFormatError(ValueError):
    """A save this build can't load. The message is meant for the user."""


class CorruptSaveError(SaveFormatError):
    """The bytes aren't a readable mantle-bloom save, or the save contradicts itself."""


class UnsupportedSaveVersionError(SaveFormatError):
    """A save (or one of its plates) in a format version this build doesn't read."""


def world_surface(world: World) -> str:
    """`"lines"`, `"quad"`, or `"empty"`. Raises `CorruptSaveError` on a world that mixes
    line-backed plates with others, which no code path produces."""
    from .legacy_conversion import is_line_plate

    if not world.plates:
        return "empty"
    lines = [is_line_plate(p) for p in world.plates]
    if all(lines):
        return "lines"
    if not any(lines):
        return "quad"
    raise CorruptSaveError("save mixes line-backed and sparse-quad plates")


def save_world_bytes(world: World) -> bytes:
    return pickle.dumps({"format": SAVE_FORMAT, "version": SAVE_FORMAT_VERSION, "surface": world_surface(world), "world": world})


def load_world_bytes(data: bytes, *, convert_lines: bool = False) -> World:
    """The `World` in a save of any supported version. Raises `CorruptSaveError` for bytes
    that aren't a readable save (pickle's own errors are wrapped, never passed through) or a
    save that contradicts itself, and `UnsupportedSaveVersionError` for a save, or a sparse
    quad plate inside one, written in a newer format than this build reads. Nothing falls back
    silently to a partial load.

    A line-backed save loads as a line world unless `convert_lines`, which converts it to
    sparse quads one way (legacy_conversion.convert_world_to_quad); the conversion summary is
    then on `World.surface_conversion`."""
    try:
        payload = pickle.loads(data)
    except SaveFormatError:
        raise
    except Exception as exc:  # noqa: BLE001 - pickle raises many types on foreign input
        raise CorruptSaveError(f"not a readable mantle-bloom save ({type(exc).__name__}: {exc})") from exc
    declared = None
    if isinstance(payload, dict) and "format" in payload:
        if payload.get("format") != SAVE_FORMAT:
            raise CorruptSaveError(f"not a mantle-bloom save (format {payload.get('format')!r})")
        version = payload.get("version")
        if not isinstance(version, int) or isinstance(version, bool) or not 2 <= version <= SAVE_FORMAT_VERSION:
            raise UnsupportedSaveVersionError(f"unsupported save format version {version!r}; this build reads 1-{SAVE_FORMAT_VERSION}")
        if version >= 3:
            declared = payload.get("surface")
            if declared not in SURFACES:
                raise CorruptSaveError(f"save declares unknown surface {declared!r}")
        world = payload.get("world")
    else:
        world = payload  # version 1: a bare pickled World
    if not isinstance(world, World):
        raise CorruptSaveError(f"expected a World, got {type(world).__name__}")
    surface = world_surface(world)
    if declared is not None and declared != surface:
        raise CorruptSaveError(f"save declares a {declared!r} world but holds a {surface!r} one")
    _backfill_added_fields(world)
    _drop_derived_caches(world)
    if convert_lines and surface == "lines":
        from .legacy_conversion import convert_world_to_quad

        # Snapshots the water budget against the quad hypsometry itself.
        convert_world_to_quad(world)
    _backfill_water_budget(world)
    return world


def _backfill_added_fields(world: World) -> None:
    """Give a `default_factory` field added to `World` after this save was written its empty
    default. Unlike a plain-default field (`steps_taken: int = 0`, a class attribute an old
    pickle falls through to), a `field(default_factory=...)` sets no class attribute, so an
    old save's `__dict__` simply won't have the key and the first access `AttributeError`s.
    Only mutable-default fields need listing here."""
    if not hasattr(world, "stranded_basin_tracks"):
        world.stranded_basin_tracks = []
    if not hasattr(world, "gap_tracks"):
        world.gap_tracks = []
    if not hasattr(world, "pending_magma_parcels"):
        world.pending_magma_parcels = []
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
    if not hasattr(world, "removed_points_log"):
        world.removed_points_log = []
    if not hasattr(world, "corner_notch_log"):
        world.corner_notch_log = []
    if not hasattr(world, "pinned_omegas"):
        world.pinned_omegas = {}
    if not hasattr(world, "stats_history"):
        world.stats_history = []
    if not hasattr(world, "phase_budget"):
        world.phase_budget = {}
    ledger_missing = not hasattr(world, "continental_material_ledger")
    if ledger_missing:
        world.continental_material_ledger = {}
        from . import continental_ledger

        # Infer the tracer and opening balance only for saves that predate the ledger.
        continental_ledger.ensure_initialized(world)
    from . import cratons

    # Accounts only: a save from before cratons existed seeds them on its first step.
    cratons.ensure_ledger(world)


def _backfill_water_budget(world: World) -> None:
    """Eustatic sea level (eustasy.py): a save written before this existed has a fixed
    sea_level_m and no water budget -- snapshot the budget from that save's own hypsometry +
    sea level so loading it doesn't jump the shoreline, then let it be conserved onward.
    Saves from before issue #257 kept the budget as a summed water column over nominal-area
    nodes (`ocean_water_column_m`); that can't be converted exactly, so it is dropped and
    re-snapshotted in m^3 the same way. Runs after any conversion, which snapshots the budget
    against the converted hypsometry itself."""
    world.__dict__.pop("ocean_water_column_m", None)
    if getattr(world, "ocean_water_volume_m3", None) is None:
        from . import eustasy

        eustasy.initialize_water_budget(world)


def _drop_derived_caches(world: World) -> None:
    """Clear every plate's lazily-rebuilt geometry cache (bounding polygon, its k-d tree,
    the contains_batch row lookup). These are pure functions of a plate's current lines +
    frame, so a stale one from an older app version -- e.g. a pre-keyhole `outline_world`
    result, or a `_RowLookup` from before it grew per-arc intervals -- would otherwise be
    trusted as-is on load. Cheap: each rebuilds on first use after this."""
    for plate in world.plates:
        reset = getattr(plate, "_reset_caches", None)
        if callable(reset):
            reset()
            continue
        invalidate = getattr(plate, "_invalidate_bounding_polygon", None)
        if callable(invalidate):
            invalidate()
    # The render path's cached node-cloud k-d tree and its positions-only sibling shared with
    # climate.py (see World.node_kdtree_cache/node_position_tree_cache) -- both pure functions
    # of the just-invalidated plate geometry, rebuilt on first use after load.
    world.node_kdtree_cache = None
    world.node_position_tree_cache = None
    # "healpix" node_cloud_resample_mode's own caches -- same rationale, see their own
    # docstrings on World. node_healpix_grid_cache is kept (it's a pure function of node
    # *count*, not the just-invalidated geometry) rather than dropped alongside the others.
    world.node_healpix_index_cache = None
    world.node_kdtree_relief_cache = None
    world.node_hillshade_cache = None

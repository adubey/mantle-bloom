"""Whole-`World` save/load to a single opaque file -- the "File > Save/Load" feature.

New saves are v5: a zip of `.npy` arrays plus a JSON manifest (see `npz_save.py`),
loaded with `allow_pickle=False`, so a v5 save carries no code. Every field on `World`
(see world.py) -- plates, mantle centers, the collision-progress dict -- is already a
plain dataclass, scalar, container or numpy array, so the zip round-trips it exactly
with no bespoke per-field (de)serialization code to keep in sync as `World`'s own
fields change. What v5 promises is in docs/save-compatibility.md: which older save
versions load, how fields added since are backfilled, and that a save this build can't
read fails with `CorruptSaveError` or `UnsupportedSaveVersionError` rather than
loading partially.

Security note: unpickling is equivalent to running arbitrary code from the file's
bytes. v5 saves avoid that entirely, but saves from before v5 (versions 1-4) still
load through the old pickle path, so the #555 F1 trust decision -- never load another
user's old save -- stands regardless of this module.
"""

from __future__ import annotations

import numpy as np

from . import npz_save
from .world import World


# The saved file is a versioned envelope:
#
# - Version 1 is the original bare pickled `World` (no envelope), still loadable.
# - Version 2 added the envelope `{"format", "version", "world"}` so a save can carry, and a
#   loader can check, which surface representations it may contain -- see
#   sparse_quad_patch.py, whose plates also version their own pickled state.
# - Version 3 adds `"surface"`: `"lines"`, `"quad"`, or `"empty"` (no plates), checked
#   against the plates actually unpickled. A world never mixes the two surfaces, and this
#   build only wrote `"quad"` or `"empty"`: a `"lines"` save converts on load.
# - Version 4 moves climate, hydrology, and surface-process modules under
#   `app.hydroclimate`. Pickle records class module paths, so earlier builds cannot read
#   the new paths.
# - Version 5 replaces pickle with a zip of `.npy` arrays plus a JSON manifest (see
#   `npz_save.py`), dropping only the erosion cache and derived scratch: every cache
#   the next step reads stays, so a post-load step lands on identical state (H9).
#   Versions 1-4 still load through the pickle path; only v5 is written.
#
# Bump on any change an older build can't read. docs/save-compatibility.md is the policy for
# line-backed (legacy) saves and how they move to quads (legacy_conversion.py).
SAVE_FORMAT = "mantle-bloom-world"
SAVE_FORMAT_VERSION = npz_save.VERSION
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
    """The default v5 save of `world`: a zip of `.npy` arrays plus a JSON manifest.
    Keeps every cache the next step reads (so post-load evolution matches), drops only
    the erosion cache and derived scratch (see `npz_save.py` for the drop sets). What
    user-facing Save, checkpoints, and diagnostics share; `save_exact_bytes` is the
    bit-identical variant for replay-log records and anywhere exactness is needed,
    and `save_compact_bytes` the smaller approximate one."""
    return npz_save.save_world_v5(world, profile="default")


def save_exact_bytes(world: World) -> bytes:
    """The exact v5 save of `world`, keeping the erosion cache too: a replay re-steps
    the record and must reproduce its fingerprint exactly, including the debug views
    that read the erosion cache."""
    return npz_save.save_world_v5(world, profile="exact")


def save_compact_bytes(world: World) -> bytes:
    """The compact v5 save of `world`: additionally drops the step-read caches
    (climate, hydrology, CFD, polarity frame) for a smaller file. Loading one
    rebuilds them on the first post-load step, which lands close to but not on the
    uninterrupted state -- check continuation with tolerances, never `state_hash`
    equality (see `test_compact_profile_step_lands_within_tolerance`). Nothing
    user-facing or replay-critical uses this profile; it is for cheap snapshots and
    transfer where approximate continuation is acceptable."""
    return npz_save.save_world_v5(world, profile="compact")


def save_kind(data: bytes) -> str:
    """Which save encoding `data` holds, sniffed from its content (never the file
    extension): `"v5"` for a zip-manifest save, `"pickle"` for a legacy pickled
    (v1-v4) save. Garbage also sniffs as `"pickle"` -- kind is a hint for debug
    scripts choosing how to report a file, not a validity check; `load_world_bytes`
    is what decides whether bytes actually load."""
    return "v5" if npz_save.is_v5(data) else "pickle"


def load_world_bytes(data: bytes) -> World:
    """The `World` in a save of any supported version. Raises `CorruptSaveError` for bytes
    that aren't a readable save (pickle's own errors are wrapped, never passed through) or a
    save that contradicts itself, and `UnsupportedSaveVersionError` for a save, or a sparse
    quad plate inside one, written in a newer format than this build reads. Nothing falls back
    silently to a partial load.

    A line-backed save is unpickled with inert stand-ins for the retired line classes
    (`legacy_conversion.legacy_unpickler`) and converted to sparse quads one way
    (`legacy_conversion.convert_world_to_quad`); the conversion summary is then on
    `World.surface_conversion`."""
    if npz_save.is_v5(data):
        world = npz_save.load_world_v5(data)
    else:
        world = _load_pickled_world(data)
    _migrate_plate_ages(world)
    _backfill_added_fields(world)
    _drop_retired_state(world)
    _drop_derived_caches(world)
    _migrate_climate_humidity_units(world)
    _backfill_water_budget(world)
    return world


def _load_pickled_world(data: bytes) -> World:
    """The `World` in a pre-v5 (pickled) save, converted to sparse quads when it is
    line-backed. Carries the envelope and surface checks the pickle format needs; v5
    validates its own manifest inside `npz_save.load_world_v5`."""
    from .legacy_conversion import convert_world_to_quad, legacy_unpickler

    try:
        payload = legacy_unpickler(data).load()
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
    if surface == "lines":
        # First: the stand-in line plates have no behaviour, and the backfills below read
        # plate fields. Snapshots the water budget against the quad hypsometry itself.
        try:
            convert_world_to_quad(world)
        except ValueError as exc:
            raise CorruptSaveError(f"line-backed save can't be converted: {exc}") from exc
    return world


def _migrate_climate_humidity_units(world: World) -> None:
    """Convert a saved `climate_cache` from the unitless humidity index to g/kg (#522). The
    cache is kept rather than dropped because the next step's transpiration source reads its
    precipitation. Its humidity is only rendered until that step recomputes it, so the index
    is scaled by the same factor its constants were (`climate.HUMIDITY_INDEX_G_KG`), not
    re-derived."""
    cache = world.climate_cache
    if cache is None or "humidity" not in vars(cache):
        return
    from .hydroclimate import climate

    cache.humidity_g_kg = vars(cache).pop("humidity") * climate.HUMIDITY_INDEX_G_KG


def _migrate_plate_ages(world: World) -> None:
    """A quad plate saved before #358 (surface format version 1) counted its age in steps.
    Convert it to years with the world's own mean step length (see
    `legacy_conversion.legacy_age_years`)."""
    from .legacy_conversion import legacy_age_years

    for plate in world.plates:
        state = plate.__dict__
        if "_age_years" in state:
            continue
        # Without either age the plate would load and then fail on its first `age_years` read.
        if "_age_steps" not in state:
            raise CorruptSaveError(f"plate {plate.plate_id!r}: no age in years or steps")
        age_steps = state.pop("_age_steps")
        if isinstance(age_steps, bool) or not isinstance(age_steps, (int, np.integer)) or age_steps < 0:
            raise CorruptSaveError(f"plate {plate.plate_id!r}: age must be a non-negative step count")
        plate.set_age_years(legacy_age_years(int(age_steps), world))


def _backfill_added_fields(world: World) -> None:
    """Give a `default_factory` field added to `World` after this save was written its empty
    default. Unlike a plain-default field (`steps_taken: int = 0`, a class attribute an old
    pickle falls through to), a `field(default_factory=...)` sets no class attribute, so an
    old save's `__dict__` simply won't have the key and the first access `AttributeError`s.
    Only mutable-default fields need listing here."""
    for name, make in (
        ("stranded_basin_tracks", list),
        ("gap_tracks", list),
        ("pending_magma_parcels", list),
        ("mantle_plumes", list),
        ("plume_stream", dict),
        ("overlap_progress", dict),
        ("faults", list),
        ("boundary_faults", list),
        ("fault_systems", list),
        ("earthquakes", list),
        ("removed_points_log", list),
        ("pinned_omegas", dict),
        ("stats_history", list),
        ("phase_budget", dict),
        ("collision_evidence", dict),
        ("collision_fronts", list),
        ("collision_polarity_stats", dict),
        ("hm_source_sink_ledger", dict),
        ("hm_suture_budget", dict),
        ("suture_transfer_stats", dict),
    ):
        if not hasattr(world, name):
            setattr(world, name, make())
    # Buoyant-crust ledger needs more than a missing-key default: backfill must
    # also repair a bare {} or partial dict, so ensure_ledger runs unconditionally.
    from . import terranes

    terranes.ensure_ledger(world)
    ledger_missing = not hasattr(world, "continental_material_ledger")
    if ledger_missing:
        world.continental_material_ledger = {}
        from . import continental_ledger

        # Infer the tracer and opening balance only for saves that predate the ledger.
        continental_ledger.ensure_initialized(world)
    from . import cratons, orogeny
    from .hydroclimate import mobile_cover

    # Accounts only: a save from before cratons existed seeds them on its first step.
    cratons.ensure_ledger(world)
    orogeny.ensure_budget(world)
    # A save from before the mobile cover holds none, so this seeds an empty opening balance.
    mobile_cover.ensure_ledger(world)


def _drop_retired_state(world: World) -> None:
    """Drop `World` attributes the line surface (retired in #251) kept, so a loaded world --
    and anything it is saved as -- carries no line state. `corner_notch_log` recorded the line
    engine's corner-notch filler; every save before #251, quad ones included, pickled it.
    `gap_fill_algorithm` chose the line engine's gap filler, in saves older still."""
    world.__dict__.pop("corner_notch_log", None)
    world.__dict__.pop("gap_fill_algorithm", None)


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
        from .hydroclimate import eustasy

        eustasy.initialize_water_budget(world)


def _drop_derived_caches(world: World) -> None:
    """Clear every plate's lazily-rebuilt geometry caches (bounding polygon, k-d trees,
    adjacency). These are pure functions of a plate's current cells + frame, so a stale one
    from an older app version would otherwise be trusted as-is on load. Cheap: each rebuilds
    on first use after this."""
    for plate in world.plates:
        plate._reset_caches()
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

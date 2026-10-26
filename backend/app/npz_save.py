"""v5 world saves: a zip of `.npy` arrays plus a JSON manifest (#596).

Replaces pickle for newly written saves. Every value the save holds is `None`, a
scalar, a string, a dict/list/tuple/set, an `np.ndarray`, or a dataclass composed of
those (see the Stage 1 inventory in `analysis/issue594/inventory.json`) -- anything
else fails closed with `TypeError` at save time rather than being silently dropped.
The two exceptions are explicit encodings: live plate references re-link by id, and the
HEALPix grid (a pure scaffold of `healpix_grid.build`) rebuilds from `nside`.
Arrays load with `allow_pickle=False`, so unlike a pickle a v5 save carries no code.

Layout: `manifest.json` holds the envelope (`format`, `version`, `surface`), the
world's scalars and small containers, one entry per plate, and references to array
entries; each array lives at `a/<path>.npy`. Array-valued fields of dataclass lists
(fault traces, magma parcels, quakes, ...) are packed into one entry per
(list, field) as a ragged concatenation plus offsets -- one entry per array would
make thousands of zip entries and dominate write time (see #594).

Two drop profiles share the codec. User saves drop the erosion cache (overwritten
whole by every erosion call, read only by a debug view); both profiles drop the
unserializable derived scratch (k-d trees, the HEALPix grid, the boundary-search
cache), which rebuilds lazily; replay-log records keep the erosion cache too. The
climate, hydrology and CFD caches and the polarity frame are always kept: the H9 hard gate
(`surface_parity_gates.load_continuation`) requires a post-load step to land on
identical authoritative state, and the next step reads all four before refreshing
them. The loader recomputes a missing `climate_cache` for a save that was stepped
(its `hydrology_cache_step` is kept) as a safety net for saves that predate it; a
never-stepped world keeps the empty-cache signal its readers rely on.
"""

from __future__ import annotations

import dataclasses
import importlib
import io
import json
import math
import zipfile

import numpy as np

from . import healpix_grid
from .sparse_quad_patch import PlateWithSparseQuadPatch
from .world import World

FORMAT = "mantle-bloom-world"
VERSION = 5
SURFACES = ("quad", "empty")

# The H9 hard gate (`surface_parity_gates.load_continuation`) requires a step from a
# loaded world to land on identical authoritative state to a step from the in-memory
# world -- so every cache the next step reads before recomputing it must survive the
# save. Measured during #596: dropping the climate, hydrology or CFD caches or the
# polarity frame diverges the next step (transpiration reads the climate cache, and the
# hydro/CFD/polarity state feeds the step before it is refreshed), while dropping the
# erosion cache (overwritten whole by every erosion call, read only by a debug view)
# steps identically. Replay-log records additionally keep the erosion cache so a replay
# re-steps exactly including its debug views.
_SAVE_DROP = frozenset({"erosion_cache"})
_DERIVED_SCRATCH = frozenset(
    {
        "land_kdtree_cache",
        "node_kdtree_cache",
        "node_position_tree_cache",
        "node_healpix_grid_cache",
        "node_healpix_index_cache",
        "node_kdtree_relief_cache",
        "node_hillshade_cache",
        "boundary_search_cache",
    }
)


def _drop_for(*, keep_caches: bool) -> frozenset:
    # Exact records keep every data cache; only unserializable scratch is dropped.
    if keep_caches:
        return _DERIVED_SCRATCH
    return _DERIVED_SCRATCH | _SAVE_DROP


def _freeze_scalar(value, path: str):
    """`None`, an exact-typed scalar, or a tagged non-finite float. Exact-type matching
    matters: `np.float64` subclasses `float`, and silently keeping one as a JSON float
    loses the type (`replay_log._feed` hashes type names)."""
    if value is None or type(value) in (bool, int, str):
        return value
    if type(value) is float:
        if math.isfinite(value):
            return value
        return {"__nonfinite__": "nan" if math.isnan(value) else ("inf" if value > 0 else "-inf")}
    if isinstance(value, bytes):
        # The manifest is JSON, so raw bytes have no encoding: fail here with the path
        # instead of dying later in json.dumps with no location (#596 review).
        raise TypeError(f"unencodable bytes at {path}; tag-encode the field or drop it from the save")
    if isinstance(value, np.generic):
        item = value.item()
        if isinstance(item, bytes):
            raise TypeError(f"unencodable numpy bytes scalar ({value.dtype}) at {path}")
        return {"__np_scalar__": {"dtype": str(value.dtype), "value": item}}
    raise TypeError(f"unencodable scalar {type(value).__name__} at {path}")


def _pack_arrays(field: str, path: str, vals: list, arrays: dict) -> dict:
    """Pack one array-valued column of a dataclass list into a single zip entry: the
    concatenation of every present array plus per-item lengths, shapes, and an explicit
    null bitmap (an empty array and a missing one both have length 0, so lengths alone
    cannot tell them apart)."""
    nulls = [i for i, v in enumerate(vals) if v is None]
    present = [(i, v) for i, v in enumerate(vals) if v is not None]
    if not present:
        return {"__packed__": f"{path}.{field}", "dtype": None, "lengths": [0] * len(vals), "shapes": [], "nulls": nulls}
    dtype = present[0][1].dtype
    for i, v in present:
        if v.dtype != dtype:
            raise TypeError(f"mixed dtypes in {path}.{field}: {dtype} vs {v.dtype} (item {i})")
    flat = np.concatenate([np.ravel(np.ascontiguousarray(v)) for _, v in present])
    entry = f"{path}.{field}"
    arrays[entry] = flat
    present_at = dict(present)
    lengths, shapes = [], []
    for i in range(len(vals)):
        v = present_at.get(i)
        lengths.append(0 if v is None else v.size)
        shapes.append([] if v is None else list(v.shape))
    assert sum(lengths) == flat.size
    return {"__packed__": entry, "dtype": str(dtype), "lengths": lengths, "shapes": shapes, "nulls": nulls}


def _freeze_dataclass_list(items: list, path: str, arrays: dict) -> dict:
    """Encode a same-typed dataclass list column-wise: scalars inline in the manifest,
    each array-valued field as one packed entry (see `_pack_arrays`)."""
    cls = type(items[0])
    if any(type(item) is not cls for item in items):
        raise TypeError(f"mixed dataclass list at {path}")
    cols = {}
    for f in dataclasses.fields(cls):
        vals = [getattr(item, f.name) for item in items]
        if all(v is None or isinstance(v, np.ndarray) for v in vals):
            cols[f.name] = _pack_arrays(f.name, path, vals, arrays)
        elif all(isinstance(v, np.generic) for v in vals):
            # Route through _freeze_scalar so numpy bytes scalars fail with a path (#596).
            cols[f.name] = [_freeze_scalar(v, f"{path}.{f.name}[{i}]") for i, v in enumerate(vals)]
        else:
            # Scalars, tuples, nested dataclasses -- encoded per item; correctness over
            # packing for the rare non-array column.
            cols[f.name] = [freeze(v, f"{path}.{f.name}[{i}]", arrays) for i, v in enumerate(vals)]
    return {"__dclist__": f"{cls.__module__}.{cls.__qualname__}", "cols": cols}


def freeze(value, path: str, arrays: dict):
    """Encode `value` into manifest-safe structure, stashing arrays in `arrays`."""
    if isinstance(value, PlateWithSparseQuadPatch):
        # A live plate reference (see `HydrologyFields.plates_in_order`): the loader
        # re-links it to the loaded plate with the same id (see `thaw`).
        return {"__plate_ref__": value.plate_id}
    if isinstance(value, healpix_grid.HealpixGrid):
        # The grid is a pure scaffold: `healpix_grid.build` is its only constructor, so
        # the loader rebuilds it from `nside` instead of carrying every derived array.
        # `NodePixelIndex` (the query cache built over a grid) is dropped, never saved.
        if not isinstance(value.nside, int) or value.nside <= 0:
            raise TypeError(f"bad HEALPix nside {value.nside!r} at {path}")
        return {"__healpix_grid__": value.nside}
    if value is None or type(value) in (bool, int, float, str, bytes) or isinstance(value, np.generic):
        return _freeze_scalar(value, path)
    if isinstance(value, np.ndarray):
        arrays[path] = np.ascontiguousarray(value)
        return {"__npy__": path}
    if isinstance(value, (set, frozenset)):
        tag = "__frozenset__" if isinstance(value, frozenset) else "__set__"
        items = sorted((freeze(v, f"{path}#{i}", arrays) for i, v in enumerate(value)), key=repr)
        return {tag: items}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        cls = type(value)
        return {
            "__dc__": f"{cls.__module__}.{cls.__qualname__}",
            "fields": {f.name: freeze(getattr(value, f.name), f"{path}.{f.name}", arrays) for f in dataclasses.fields(value)},
        }
    if isinstance(value, (list, tuple)):
        if value and all(dataclasses.is_dataclass(v) and not isinstance(v, type) for v in value):
            if len({type(v) for v in value}) == 1:
                return _freeze_dataclass_list(list(value), path, arrays)
        tag = "__tuple__" if isinstance(value, tuple) else "__list__"
        return {tag: [freeze(v, f"{path}[{i}]", arrays) for i, v in enumerate(value)]}
    if isinstance(value, dict):
        return {"__dict__": [[freeze(k, f"{path}<k{i}>", arrays), freeze(v, f"{path}<{i}>", arrays)] for i, (k, v) in enumerate(value.items())]}
    raise TypeError(f"unencodable {type(value).__name__} at {path}")


def _resolve(dotted: str):
    """Import `a.b.C` by trying the longest importable module prefix first."""
    parts = dotted.split(".")
    for i in range(len(parts) - 1, 0, -1):
        try:
            obj = importlib.import_module(".".join(parts[:i]))
        except ImportError:
            continue
        for part in parts[i:]:
            obj = getattr(obj, part)
        return obj
    raise ImportError(dotted)


def _thaw_packed(node: dict, zf: zipfile.ZipFile) -> list:
    """Split one packed column back into per-item arrays (`None` exactly where the null
    bitmap says so -- a present-but-empty array keeps its shape)."""
    if node["dtype"] is None:
        return [None] * len(node["lengths"])
    nulls = set(node["nulls"])
    try:
        with zf.open(f"a/{node['__packed__']}.npy") as fh:
            flat = np.load(fh, allow_pickle=False)
    except KeyError as exc:
        raise _corrupt(f"v5 save is missing array {node['__packed__']!r}") from exc
    except ValueError as exc:
        raise _corrupt(f"v5 save holds an unreadable array {node['__packed__']!r}: {exc}") from exc
    dtype = np.dtype(node["dtype"])
    if flat.dtype != dtype:
        raise _corrupt(f"v5 array {node['__packed__']!r} changed dtype {flat.dtype} to {dtype}")
    items, cursor = [], 0
    for i, (length, shape) in enumerate(zip(node["lengths"], node["shapes"])):
        try:
            items.append(None if i in nulls else flat[cursor : cursor + length].reshape(shape))
        except ValueError as exc:
            raise _corrupt(f"v5 array {node['__packed__']!r} item {i} does not fit shape {shape}: {exc}") from exc
        cursor += length
    if cursor != flat.size:
        raise _corrupt(f"v5 array {node['__packed__']!r} has inconsistent offsets")
    return items


def thaw(node, zf: zipfile.ZipFile, plates_by_id: dict | None = None):
    """Rebuild a `freeze`d value, reading arrays from the open zip. Plate references
    (`__plate_ref__`) resolve against `plates_by_id`, so the world dict must thaw after
    the plates are built; a reference with no table (or to an unknown id) is corrupt."""
    if isinstance(node, dict) and set(node) == {"__dc__", "fields"}:
        cls = _resolve(node["__dc__"])
        try:
            return cls(**{k: thaw(v, zf, plates_by_id) for k, v in node["fields"].items()})
        except TypeError as exc:
            raise _corrupt(f"v5 dataclass {node['__dc__']!r} changed shape: {exc}") from exc
    if isinstance(node, dict) and set(node) == {"__dclist__", "cols"}:
        cls = _resolve(node["__dclist__"])
        names = [f.name for f in dataclasses.fields(cls)]
        if set(node["cols"]) != set(names):
            raise _corrupt(f"v5 dataclass list {node['__dclist__']!r} changed fields")
        cols = {k: _thaw_column(v, zf, plates_by_id) for k, v in node["cols"].items()}
        try:
            return [cls(**{k: cols[k][i] for k in names}) for i in range(len(next(iter(cols.values()))))]
        except (TypeError, IndexError) as exc:
            raise _corrupt(f"v5 dataclass list {node['__dclist__']!r} changed shape: {exc}") from exc
    if isinstance(node, dict) and len(node) == 1:
        [[tag, payload]] = node.items()
        return _thaw_tagged(tag, payload, node, zf, plates_by_id)
    if isinstance(node, dict):
        raise _corrupt(f"v5 manifest holds a bare dict where a tagged value belongs: {str(node)[:120]}")
    return node


def _thaw_tagged(tag: str, payload, node: dict, zf: zipfile.ZipFile, plates_by_id: dict | None):
    """Rebuild one tagged manifest value. Kept separate so `thaw` stays under the
    complexity cap."""
    if tag == "__plate_ref__":
        if plates_by_id is None or payload not in plates_by_id:
            raise _corrupt(f"v5 save references unknown plate {payload!r}")
        return plates_by_id[payload]
    if tag == "__healpix_grid__":
        try:
            return healpix_grid.build(int(payload))
        except Exception as exc:
            raise _corrupt(f"v5 save holds a bad HEALPix nside {payload!r}: {exc}") from exc
    if tag == "__npy__":
        try:
            with zf.open(f"a/{payload}.npy") as fh:
                return np.load(fh, allow_pickle=False)
        except KeyError as exc:
            raise _corrupt(f"v5 save is missing array {payload!r}") from exc
        except ValueError as exc:
            raise _corrupt(f"v5 save holds an unreadable array {payload!r}: {exc}") from exc
    if tag == "__packed__":
        return _thaw_packed(node, zf)
    if tag == "__nonfinite__":
        return {"nan": math.nan, "inf": math.inf, "-inf": -math.inf}[payload]
    if tag == "__np_scalar__":
        return np.dtype(payload["dtype"]).type(payload["value"])
    if tag in ("__list__", "__tuple__"):
        items = [thaw(v, zf, plates_by_id) for v in payload]
        return tuple(items) if tag == "__tuple__" else items
    if tag in ("__set__", "__frozenset__"):
        items = [thaw(v, zf, plates_by_id) for v in payload]
        try:
            return frozenset(items) if tag == "__frozenset__" else set(items)
        except TypeError as exc:
            raise _corrupt(f"v5 set holds an unhashable item: {exc}") from exc
    if tag == "__dict__":
        try:
            return {thaw(k, zf, plates_by_id): thaw(v, zf, plates_by_id) for k, v in payload}
        except TypeError as exc:
            raise _corrupt(f"v5 dict holds an unhashable key: {exc}") from exc
    raise _corrupt(f"v5 manifest holds an unknown tag {tag!r}")


def _thaw_column(col, zf: zipfile.ZipFile, plates_by_id: dict | None = None) -> list:
    if isinstance(col, dict) and set(col) == {"__packed__", "dtype", "lengths", "shapes", "nulls"}:
        return _thaw_packed(col, zf)
    if isinstance(col, list):
        return [thaw(v, zf, plates_by_id) for v in col]
    raise _corrupt(f"v5 dataclass-list column changed shape: {str(col)[:120]}")


def _corrupt(reason: str) -> Exception:
    """`persistence.CorruptSaveError`, imported lazily: this module is imported by
    `persistence`, so a module-level import would be circular."""
    from .persistence import CorruptSaveError

    return CorruptSaveError(reason)


def _unsupported(reason: str) -> Exception:
    from .persistence import UnsupportedSaveVersionError

    return UnsupportedSaveVersionError(reason)


def is_v5(data: bytes) -> bool:
    """Whether `data` is a v5 zip save (rather than a pickled v1-v4 save)."""
    return data[:4] == b"PK\x03\x04"


def save_world_v5(world: World, *, keep_caches: bool) -> bytes:
    """Serialize `world` as a v5 zip (`manifest.json` plus one `.npy` per array).
    `keep_caches=False` (user saves, checkpoints) drops the erosion cache;
    `keep_caches=True` (replay-log records) keeps it too, so a replay re-steps exactly
    including its debug views. Unserializable derived scratch is dropped by both
    profiles; every other cache is kept so a post-load step lands on identical state (H9)."""
    from . import sparse_quad_patch
    from .persistence import world_surface
    from .sparse_quad_patch import PlateWithSparseQuadPatch

    # v5 never holds line-backed plates: they only exist transiently while an old save
    # converts on load. `world_surface` refuses a mixed world with the same message the
    # pickle path always raised.
    surface = world_surface(world)
    if surface == "lines":
        raise _corrupt("v5 saves cannot hold line-backed plates")
    drop = _drop_for(keep_caches=keep_caches)
    skip = set(PlateWithSparseQuadPatch._TOPOLOGY_CACHES) | set(PlateWithSparseQuadPatch._GEOMETRY_CACHES)
    arrays: dict[str, np.ndarray] = {}
    plates = []
    for i, plate in enumerate(world.plates):
        state = {k: freeze(v, f"plates/{i}/{k}", arrays) for k, v in vars(plate).items() if k not in skip}
        state["_plate_surface_version"] = sparse_quad_patch.QUAD_SURFACE_FORMAT_VERSION
        plates.append(state)
    manifest = {
        "format": FORMAT,
        "version": VERSION,
        "surface": surface,
        "plates": plates,
        # "plates" itself is excluded here: it rides in the per-plate list above.
        "world": {k: freeze(v, f"world/{k}", arrays) for k, v in vars(world).items() if k not in drop and k != "plates"},
    }
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_STORED) as zf:
        zf.writestr("manifest.json", json.dumps(manifest))
        for path in sorted(arrays):
            # Write straight into the zip: `np.save` streams header-then-data, so the
            # unseekable entry works and each array is copied once, not twice (#596).
            with zf.open(f"a/{path}.npy", "w") as entry:
                np.save(entry, arrays[path], allow_pickle=False)
    return buf.getvalue()


def load_world_v5(data: bytes) -> World:
    """The `World` in v5 bytes. Raises `CorruptSaveError` for anything that is not a
    readable v5 save and `UnsupportedSaveVersionError` for a newer v5 or a plate
    surface format this build does not read. Never loads partially."""
    from .sparse_quad_patch import (
        LEGACY_QUAD_SURFACE_FORMAT_VERSIONS,
        QUAD_SURFACE_FORMAT_VERSION,
        PlateWithSparseQuadPatch,
    )

    from .persistence import SaveFormatError

    try:
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            try:
                manifest = json.loads(zf.read("manifest.json"))
            except Exception as exc:
                raise _corrupt(f"not a readable mantle-bloom save ({type(exc).__name__}: {exc})") from exc
            if not isinstance(manifest, dict) or manifest.get("format") != FORMAT:
                raise _corrupt(f"not a mantle-bloom save (format {manifest.get('format')!r})" if isinstance(manifest, dict) else "not a mantle-bloom save")
            version = manifest.get("version")
            if version != VERSION:
                raise _unsupported(f"unsupported save format version {version!r}; this build reads v5 plates in 1-{QUAD_SURFACE_FORMAT_VERSION}")
            if manifest.get("surface") not in SURFACES:
                raise _corrupt(f"save declares unknown surface {manifest.get('surface')!r}")
            try:
                plates = []
                for state in manifest["plates"]:
                    plate_version = state.get("_plate_surface_version")
                    if plate_version != QUAD_SURFACE_FORMAT_VERSION and plate_version not in LEGACY_QUAD_SURFACE_FORMAT_VERSIONS:
                        raise _unsupported(
                            f"unsupported sparse quad surface format version {plate_version!r}; "
                            f"this build reads versions 1-{QUAD_SURFACE_FORMAT_VERSION}"
                        )
                    plate = PlateWithSparseQuadPatch.__new__(PlateWithSparseQuadPatch)
                    plate.__dict__.update({k: thaw(v, zf) for k, v in state.items() if k != "_plate_surface_version"})
                    plate._reset_caches()
                    plate._validate_leaf_topology()
                    plates.append(plate)
                # Plate references elsewhere (see `HydrologyFields.plates_in_order`) re-link to
                # the loaded plates by id, so the cache reads the world's own plates.
                plates_by_id = {plate.plate_id: plate for plate in plates}
                if len(plates_by_id) != len(plates):
                    raise _corrupt("v5 save holds two plates with one id")
                world_dict = {k: thaw(v, zf, plates_by_id) for k, v in manifest["world"].items()}
                actual = "empty" if not plates else "quad"
                if manifest.get("surface") != actual:
                    raise _corrupt(f"save declares a {manifest.get('surface')!r} world but holds a {actual!r} one")
            except Exception as exc:  # noqa: BLE001 - every unreadable shape maps to CorruptSaveError
                if isinstance(exc, SaveFormatError):
                    raise
                raise _corrupt(f"v5 save contradicts itself ({type(exc).__name__}: {exc})") from exc
    except Exception as exc:
        if isinstance(exc, SaveFormatError):
            raise
        raise _corrupt(f"not a readable mantle-bloom save ({type(exc).__name__}: {exc})") from exc
    world = World.__new__(World)
    world.__dict__.update(world_dict)
    world.plates = plates
    for name in _SAVE_DROP | _DERIVED_SCRATCH:
        world.__dict__.setdefault(name, None)
    _rebuild_climate(world)
    return world


def _rebuild_climate(world: World) -> None:
    """Recompute `climate_cache` for a save that was stepped, so the first post-load
    step's transpiration reads a real field instead of the zeros path. A never-stepped
    world (`hydrology_cache_step is None`) keeps the empty-cache signal its readers
    rely on, as does a stepped world the record profile saved with its cache."""
    if world.climate_cache is not None or world.hydrology_cache_step is None or not world.plates:
        return
    from .hydroclimate import climate

    world.climate_cache = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density))


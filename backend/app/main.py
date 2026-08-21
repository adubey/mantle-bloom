"""FastAPI routes. Single in-memory world for v1 (see docs/architecture.md) -- no world id,
no persistence, one world at a time, matching the "elevation view only" v1 scope.
"""

from __future__ import annotations

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from scipy.spatial import cKDTree

from . import climate, coastline, geometry, hydrology, lakes, plates, projections, render_image, stats
from .world import DEFAULT_MANTLE_CENTERS, World, generate_world, step_world

# A generous ceiling on requested image dimensions -- width/height come straight from the
# client's query string, and PIL will happily try to allocate whatever it's told, so an
# unbounded value here would let a malformed or malicious request force a huge allocation.
MAX_RENDER_DIMENSION_PX = 4000

app = FastAPI(title="mantle-bloom")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://localhost:5174"],
    allow_methods=["*"],
    allow_headers=["*"],
)

_state: dict[str, World | None] = {"world": None}


class GenerateRequest(BaseModel):
    seed: int = 0
    # Optional: the world tiles itself into a plausible plate count when omitted (see
    # plates.generate_plates) -- the frontend doesn't ask for one.
    num_plates: int | None = None
    # The UI's "continental plates" and "initial land" sliders (0 to 1) -- optional (each
    # falls back to its own default behavior, see plates.generate_plates) but the frontend
    # always sends both.
    continental_fraction: float | None = None
    land_fraction: float | None = None
    # The UI's "axial tilt" slider (degrees) -- optional, falls back to
    # world.DEFAULT_AXIAL_TILT_DEG, but the frontend always sends it.
    axial_tilt_deg: float | None = None
    num_mantle_centers: int = DEFAULT_MANTLE_CENTERS
    # The UI's "point density" choice -- how many elevation-line nodes each plate starts
    # with, relative to plates.DEFAULT_NODE_DENSITY (1.0 = plates.TARGET_LINE_SPACING_KM's
    # own default spacing). Validated against plates.NODE_DENSITY_CHOICES below rather than
    # accepted as an arbitrary float -- there's no continuous "in-between" density the UI
    # offers, only a fixed set of multipliers.
    node_density: float = plates.DEFAULT_NODE_DENSITY
    # The UI's "initial soil maturity" slider (0 to 1) -- optional, falls back to
    # world.generate_world's own default (0.0, a fully barren starting world -- see
    # geology.seed_initial_soil), but the frontend always sends it.
    initial_soil_maturity: float | None = None


class StepRequest(BaseModel):
    years: float


class ControlsRequest(BaseModel):
    # Both optional and independently settable -- the UI's "Controls" window (see
    # StatsModal.tsx's "Controls" button) sends only whichever slider the user actually
    # moved. Omitted fields leave that World attribute untouched. See World.sea_level_m/
    # World.solar_multiplier for what each one does physically.
    sea_level_m: float | None = None
    solar_multiplier: float | None = None


def _parse_view_rotation(rotation: str | None) -> np.ndarray:
    """The map's current view orientation (see docs/simulation-model.md#rotating-the-view),
    a row-major 3x3 rotation matrix as 9 comma-separated floats -- default identity (today's
    behavior, center at lat=0/lon=0) when omitted. This is purely a render-time transform;
    it's never stored on `World` (see world.py) since it's client-local view state, not
    simulation state. Only shape/finiteness is validated here, not strict orthonormality --
    the frontend only ever sends matrices built by composing exact rotation matrices, so any
    drift is negligible, and this guard exists to reject malformed/malicious input before it
    reaches numpy, not to re-derive a "correct" rotation from an arbitrary one."""
    if rotation is None:
        return np.eye(3)
    parts = rotation.split(",")
    if len(parts) != 9:
        raise HTTPException(status_code=400, detail=f"rotation must be 9 comma-separated floats, got {len(parts)}")
    try:
        values = [float(p) for p in parts]
    except ValueError as exc:
        raise HTTPException(status_code=400, detail="rotation must be 9 comma-separated floats") from exc
    if not all(np.isfinite(values)):
        raise HTTPException(status_code=400, detail="rotation values must be finite")
    return np.array(values, dtype=float).reshape(3, 3)


def _require_world() -> World:
    world = _state["world"]
    if world is None:
        raise HTTPException(status_code=404, detail="no world generated yet")
    return world


def _summary(world: World) -> dict:
    return {
        "seed": world.seed,
        "elapsed_years": world.elapsed_years,
        "num_plates": len(world.plates),
        # Full current log (bounded, see world.MAX_EVENT_LOG_LENGTH) every time -- simpler
        # for the frontend than diffing "new since last call", and small enough not to
        # matter on the wire.
        "events": [{"elapsed_years": e, "message": m} for e, m in world.events],
    }


#  Rounding every coordinate sent to the client to this many decimal places -- plenty for a
# unit sphere (1e-6 rad is sub-meter at planet scale, far finer than anything visible on
# screen) -- meaningfully shrinks the /world/plates payload now that it includes every node
# point per plate (`_plate_summary`'s "points" field), not just the much shorter outline/
# ellipse loops.
_COORD_DECIMALS = 6


def _round_coords(arr: np.ndarray) -> list:
    return np.round(arr, _COORD_DECIMALS).tolist()


def _plate_summary(plate: plates.Plate) -> dict:
    outline = plate.outline_world()
    node_points, _ = plate.all_points_and_elevation()
    ellipse = plates.plate_bounding_ellipse(node_points)
    return {
        "plate_id": plate.plate_id,
        "crust_type": plate.crust_type,
        # Lines with zero nodes are a real state a plate can be in (see merge_split.py's
        # own "no_land"/consumption checks) -- excluded here to match outline_world()'s own
        # filtering, so this doesn't look inconsistent next to num_points.
        "num_rows": sum(1 for line in plate.lines if len(line.theta) > 0),
        "num_points": plate.node_count(),
        "outline": _round_coords(outline),
        # Every node's own position (not just the outline loop) -- lets the client plot each
        # point individually (see docs/simulation-model.md#plate-inspector), highlighted for
        # the selected plate and dimmed for every other one.
        "points": _round_coords(node_points),
        "bounding_ellipse": None
        if ellipse is None
        else {
            "center_xyz": _round_coords(ellipse.center_xyz),
            "diameter_a_km": ellipse.diameter_a_km,
            "diameter_b_km": ellipse.diameter_b_km,
            "outline": _round_coords(ellipse.outline_xyz),
        },
    }


def _river_summary(fields: hydrology.HydrologyFields, river: hydrology.RiverInfo, river_id: int) -> dict:
    targets = fields.flow_target[river.member_idx]
    has_target = targets >= 0
    from_points = fields.points[river.member_idx[has_target]]
    to_points = fields.points[targets[has_target]]
    segments = np.round(np.stack([from_points, to_points], axis=1), _COORD_DECIMALS).tolist()
    return {
        "river_id": river_id,
        "num_nodes": int(len(river.member_idx)),
        "segments": segments,
        "mouth_xyz": _round_coords(fields.points[river.mouth_idx]),
        "mouth_type": river.mouth_type,
        "flow_rate": river.flow_rate,
        "speed": river.speed,
        "num_tributaries": river.num_tributaries,
    }


def _coastline_segments_json(world: World) -> list:
    point_a, point_b = coastline.compute_coastline_segments(world)
    if len(point_a) == 0:
        return []
    return np.round(np.stack([point_a, point_b], axis=1), _COORD_DECIMALS).tolist()


def _leaf_lakes_by_node(forest: list[lakes.Lake]) -> dict[int, lakes.Lake]:
    """Every land node that belongs to some real catchment (i.e. isn't on a chain that drains
    straight to the ocean without ever passing a local minimum -- see lakes.py's own
    `_OCEAN_CATCHMENT`), mapped to its own *leaf* `Lake` -- the basin the Lake Inspector shows
    for a click that lands on dry, never-merged land."""
    by_node: dict[int, lakes.Lake] = {}
    for lake in lakes.iter_all_lakes(forest):
        if lake.children:
            continue
        for member in lake.members.tolist():
            by_node[member] = lake
    return by_node


def _smallest_lake_containing(forest: list[lakes.Lake], wet_set: set[int]) -> lakes.Lake | None:
    """The smallest `Lake` tree node (leaf or merged parent) whose own `members` fully contains
    `wet_set` -- i.e. the basin exactly as currently manifested by a connected body of standing
    water that may already span more than one original leaf catchment. Searches top-down rather
    than needing parent pointers on `Lake` (which nothing else in lakes.py needs): a level whose
    `members` doesn't yet contain the whole wet set can't be the answer, so every other branch at
    that level is skipped."""

    def search(lake: lakes.Lake) -> lakes.Lake | None:
        if not wet_set.issubset(lake.members.tolist()):
            return None
        for child in lake.children:
            found = search(child)
            if found is not None:
                return found
        return lake

    for root in forest:
        found = search(root)
        if found is not None:
            return found
    return None


def _lake_components_sorted(fields: hydrology.HydrologyFields) -> list[np.ndarray]:
    """Every currently-visible lake's node group (`hydrology.lake_components`), in a fixed,
    deterministic order (by each component's own lowest node index) -- so the index into this
    list (`lake_id`) means the same thing across a `/world/lakes` call and a same-step
    `/world/lake_at` call, the same "meaningful only against this step" contract `river_id`
    already has (see hydrology.group_rivers)."""
    components = hydrology.lake_components(fields.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M, fields.neighbor_idx)
    components.sort(key=lambda members: int(members.min()))
    return components


def _lake_basin_summary(
    fields: hydrology.HydrologyFields, lake: lakes.Lake, lake_id: int | None, rivers: list[hydrology.RiverInfo], is_lake: bool
) -> dict:
    """One basin's full inspector payload, for both a currently-wet lake (`is_lake=True`, `lake`
    is whatever `_smallest_lake_containing` resolved) and a dry, never-flooded basin
    (`is_lake=False`, `lake` is the leaf catchment straight from `_leaf_lakes_by_node`) -- same
    shape either way, since "a basin nonetheless" (the user's own framing) should read as the
    same kind of information, just without live water. `floor_*` is this basin's own lowest
    point (`lake.floor_elevation`'s own node); `outlet_*` is "the lowest point of the edge of the
    basin" -- `lake.max_depth`'s own saddle node, `None` for an unresolved closed/endorheic basin
    with no known spill (a legitimate case, see lakes.py's own docstring), not a bug. Rivers
    "feeding into" this basin are simply every `RiverInfo` whose own mouth lands on one of this
    basin's members, regardless of that river's own `mouth_type` label ("lake", or "other" for a
    river dead-ending in a currently-dry sink) -- from this basin's own point of view both are
    equally "a river that ends here"."""
    members = lake.members
    effective_elev = fields.elevation + fields.silt_depth

    # `lake.members` is the *geometric* catchment -- every node on the way down to this basin's
    # own floor, not just whichever of them are currently underwater (a leaf's own catchment
    # routinely includes dry higher ground, see lakes.py's own docstring). For a currently-wet
    # lake, the *displayed* extent (what's drawn as water, and what a click needs to land on to
    # read back as "lake") should be the actual flooded subset (`fields.lake_depth`, this step's
    # real, continuity-resolved water state), not the full basin outline -- a dry basin has no
    # such distinction to make, so it keeps showing its whole catchment.
    wet_mask = fields.lake_depth[members] > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    wet_members = members[wet_mask]
    display_members = wet_members if is_lake and len(wet_members) > 0 else members

    floor_node = int(display_members[np.argmin(effective_elev[display_members])])
    outlet_elev = lake.max_depth
    outlet_node = lake.outlet_node_idx

    # `lake_depth` is depth above the *silt-adjusted* bed (lakes.py's own `effective_elevation`
    # basis, see step_lakes), so the actual water surface needs `silt_depth` added back in too,
    # not just bare `elevation` -- otherwise a silted-in lake reads as underwater at its own bed.
    water_elevation = (
        float(np.max(fields.elevation[wet_members] + fields.silt_depth[wet_members] + fields.lake_depth[wet_members]))
        if len(wet_members) > 0
        else None
    )
    is_spilling = outlet_elev is not None and water_elevation is not None and water_elevation >= outlet_elev

    # Rivers "feeding into" this basin are matched against the full geometric catchment, not
    # just `display_members` -- a river can plausibly end anywhere in the basin (e.g. a dry
    # higher-elevation sink not yet part of the standing water), and that's still meaningfully
    # "a river that ends in this basin" from the user's point of view.
    member_set = set(members.tolist())
    inflow_rivers = [
        {
            "mouth_xyz": _round_coords(fields.points[river.mouth_idx]),
            "flow_rate": river.flow_rate,
            "num_nodes": int(len(river.member_idx)),
        }
        for river in rivers
        if river.mouth_idx in member_set
    ]

    return {
        "lake_id": lake_id,
        "is_lake": is_lake,
        "member_count": int(len(display_members)),
        "member_xyz": _round_coords(fields.points[display_members]),
        "floor_xyz": _round_coords(fields.points[floor_node]),
        "floor_elevation_m": float(effective_elev[floor_node]),
        "outlet_xyz": None if outlet_elev is None else _round_coords(fields.points[outlet_node]),
        "outlet_elevation_m": outlet_elev,
        "water_elevation_m": water_elevation,
        "is_spilling": bool(is_spilling),
        "inflow_rivers": inflow_rivers,
    }


@app.post("/world/generate")
def generate(req: GenerateRequest) -> dict:
    if req.node_density not in plates.NODE_DENSITY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown node_density {req.node_density!r}; choices are {plates.NODE_DENSITY_CHOICES}")
    world = generate_world(
        req.seed,
        num_plates=req.num_plates,
        continental_fraction=req.continental_fraction,
        land_fraction=req.land_fraction,
        num_mantle_centers=req.num_mantle_centers,
        axial_tilt_deg=req.axial_tilt_deg,
        node_density=req.node_density,
        initial_soil_maturity=req.initial_soil_maturity,
    )
    _state["world"] = world
    return _summary(world)


@app.post("/world/step")
def step(req: StepRequest) -> dict:
    world = _require_world()
    step_world(world, req.years)
    return _summary(world)


@app.get("/world/render")
def render(
    projection: str = "behrmann", view: str = "elevation", width: int = 1100, height: int = 611, rotation: str | None = None
) -> dict:
    """Renders `view` of the current world, in `projection`, as a `width`x`height` PNG,
    returned base64-encoded in `image_base64` -- see render_image.py for the actual
    rasterization. `rotation` is the map's current view orientation (see
    _parse_view_rotation), default identity. `400` for an unrecognized projection/view, an
    out-of-range width/height, or a malformed rotation, `404` if no world has been generated
    yet."""
    world = _require_world()
    if projection not in projections.PROJECTIONS:
        raise HTTPException(status_code=400, detail=f"unknown projection {projection!r}")
    if view not in render_image.VIEWS:
        raise HTTPException(status_code=400, detail=f"unknown view {view!r}; choices are {render_image.VIEWS}")
    if not (1 <= width <= MAX_RENDER_DIMENSION_PX and 1 <= height <= MAX_RENDER_DIMENSION_PX):
        raise HTTPException(status_code=400, detail=f"width/height must be in [1, {MAX_RENDER_DIMENSION_PX}]")
    view_rotation = _parse_view_rotation(rotation)

    return {
        "projection": projection,
        "elapsed_years": world.elapsed_years,
        "image_base64": render_image.render_png_base64(world, projection, view, width, height, view_rotation),
    }


@app.post("/world/controls")
def set_controls(req: ControlsRequest) -> dict:
    """Live-adjusts sea level and/or solar heat on the current world (the "Controls" window
    -- see World.sea_level_m/World.solar_multiplier) and immediately recomputes
    world.climate_cache so the very next /world/render or /world/stats call reflects the
    change without waiting for a step -- climate_cache is otherwise only refreshed by
    erosion.py once per step (see World.climate_cache), which would make a real-time control
    feel laggy or unresponsive between steps. `404` if no world has been generated yet."""
    world = _require_world()
    if req.sea_level_m is not None:
        world.sea_level_m = req.sea_level_m
    if req.solar_multiplier is not None:
        world.solar_multiplier = req.solar_multiplier
    world.climate_cache = climate.compute_climate(world)
    return {"sea_level_m": world.sea_level_m, "solar_multiplier": world.solar_multiplier}


@app.get("/world/plates")
def list_plates() -> dict:
    """Every plate's outline + metadata (row/point counts, bounding ellipse) as JSON, for the
    "Plate Inspector" map mode -- unlike /world/render, the client renders this itself
    interactively rather than receiving a baked PNG. Un-rotated/true-frame throughout (no
    `rotation` param): the client applies its current view rotation only at draw time, same
    philosophy as climate/render-grid geometry (see docs/simulation-model.md#rotating-the-view).
    `404` if no world has been generated yet."""
    world = _require_world()
    return {"elapsed_years": world.elapsed_years, "plates": [_plate_summary(p) for p in world.plates]}


@app.get("/world/stats")
def get_stats() -> dict:
    """Aggregate physical/temperature/precipitation statistics for the Stats panel -- see
    stats.py. Stateless: recomputed fresh from the current world state every call (the
    frontend is what accumulates a history over time, by calling this after every
    generate/step, same as it already does for /world/render). `404` if no world has been
    generated yet."""
    world = _require_world()
    return stats.compute_stats(world)


@app.get("/world/plate_at")
def plate_at(lat_deg: float, lon_deg: float) -> dict:
    """Which plate owns the node nearest (lat_deg, lon_deg) -- the Plate Inspector's click
    hit-test. The client unprojects its click (through whatever view rotation is currently
    active) to a true lat/lon before calling this, so this endpoint itself never needs to
    know about rotation. `400` for non-finite input, `404` if no world has been generated
    yet."""
    world = _require_world()
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    return {"plate_id": plates.nearest_plate_id(world.plates, query_xyz)}


@app.get("/world/rivers")
def list_rivers() -> dict:
    """Every distinct river network's path + endpoint metadata, plus the current coastline, as
    JSON, for the "River Inspector" map mode -- same philosophy as /world/plates: the client
    renders and drives this itself rather than receiving a baked PNG (see
    docs/simulation-model.md#hydrology and #river-inspector). Rivers are grouped fresh from
    the current world.hydrology_cache (see hydrology.group_rivers) on every call rather than
    persisted -- hydrology_cache itself is already at most one step stale by design (see
    World.hydrology_cache), and grouping it is cheap. `coastline_segments` is computed the
    same way and drawn server-side in the temperature/humidity/precipitation views (see
    coastline.py) -- included here too since this view has no other land/ocean cue at all
    (unlike those, it draws no filled backdrop). Both are `[]` before the first step
    (climate_cache/hydrology_cache are None until erosion.py runs once). Un-rotated/
    true-frame throughout, same as /world/plates. `404` if no world has been generated yet."""
    world = _require_world()
    if world.hydrology_cache is None:
        return {"elapsed_years": world.elapsed_years, "rivers": [], "coastline_segments": []}
    rivers = hydrology.group_rivers(world.hydrology_cache)
    return {
        "elapsed_years": world.elapsed_years,
        "rivers": [_river_summary(world.hydrology_cache, river, i) for i, river in enumerate(rivers)],
        "coastline_segments": _coastline_segments_json(world),
    }


@app.get("/world/river_at")
def river_at(lat_deg: float, lon_deg: float) -> dict:
    """Which river network owns the node nearest (lat_deg, lon_deg) -- the River Inspector's
    click hit-test, same pattern as /world/plate_at. `river_id` is an index into whatever
    /world/rivers returns for *this same step* (rivers are regrouped fresh every call, not
    given a persistent identity -- see hydrology.group_rivers), so the client should treat it
    as only meaningful against its most recent /world/rivers response. `400` for non-finite
    input, `404` if no world has been generated yet."""
    world = _require_world()
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    if world.hydrology_cache is None:
        return {"river_id": None}
    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    rivers = hydrology.group_rivers(world.hydrology_cache)
    return {"river_id": hydrology.river_at(world.hydrology_cache, rivers, query_xyz)}


@app.get("/world/lakes")
def list_lakes() -> dict:
    """Every currently-visible lake's basin info, for the "Lake Inspector" map mode -- same
    "raw JSON, client renders it" philosophy as `/world/rivers`. A "lake" here means one
    connected component of `lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M` nodes (see
    `hydrology.lake_components`) -- exactly what's drawn as standing water everywhere else in
    this codebase -- resolved back up to its own basin in this step's freshly-rebuilt
    depression hierarchy (`lakes.build_lake_hierarchy`) so a lake that's currently the result of
    several original catchments merging together reports its own *current* outlet, not some
    interior saddle that's already flooded. Regrouped fresh from `world.hydrology_cache` on
    every call, same as `/world/rivers` -- `lake_id` is only meaningful against this same
    response, not a stable identity across steps. `coastline_segments` included for the same
    reason `/world/rivers` includes it: this view draws no filled backdrop of its own. Both are
    `[]` before the first step. `404` if no world has been generated yet."""
    world = _require_world()
    fields = world.hydrology_cache
    if fields is None:
        return {"elapsed_years": world.elapsed_years, "lakes": [], "coastline_segments": []}
    effective_elev = fields.elevation + fields.silt_depth
    forest = lakes.build_lake_hierarchy(effective_elev, fields.is_ocean, fields.neighbor_idx)
    components = _lake_components_sorted(fields)
    rivers = hydrology.group_rivers(fields)
    result = []
    for lake_id, members in enumerate(components):
        lake_obj = _smallest_lake_containing(forest, set(members.tolist()))
        if lake_obj is None:
            continue  # defensive only -- every wet node belongs to some catchment by construction
        result.append(_lake_basin_summary(fields, lake_obj, lake_id, rivers, is_lake=True))
    return {"elapsed_years": world.elapsed_years, "lakes": result, "coastline_segments": _coastline_segments_json(world)}


@app.get("/world/lake_at")
def lake_at(lat_deg: float, lon_deg: float) -> dict:
    """The Lake Inspector's click hit-test, resolving to whichever *basin* owns the node
    nearest `(lat_deg, lon_deg)` -- not just a currently-wet lake, so a click on ordinary dry
    land still returns "the basin nonetheless" (the user's own spec), including its lowest
    point, as long as that land actually drains into some enclosed catchment at all (see
    `kind` below). `400` for non-finite input, `404` if no world has been generated yet.

    `kind` is one of:
    - `"lake"` -- the node is currently flooded; `basin` is that connected lake's own info
      (`lake_id` set, matching the same node's `/world/lakes` entry for this same step).
    - `"basin"` -- dry land, but part of a real (possibly still-unresolved/endorheic) basin;
      `basin` is that basin's own leaf catchment info (`lake_id` is `None` -- a dry basin has
      no visible-lake identity to be stable against).
    - `"no_basin"` -- dry land whose own steepest-descent chain drains straight to the ocean
      without ever passing through a local minimum first -- an ordinary hillslope, never part
      of any basin (see lakes.py's own `_OCEAN_CATCHMENT`). `basin` is `null`.
    - `"ocean"` -- the nearest node is open ocean. `basin` is `null`.
    """
    world = _require_world()
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    fields = world.hydrology_cache
    if fields is None or len(fields.points) == 0:
        return {"kind": "no_basin", "basin": None}

    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    tree = cKDTree(fields.points)
    _, nearest = tree.query(query_xyz.reshape(1, 3), k=1)
    node_idx = int(nearest[0])

    if fields.is_ocean[node_idx]:
        return {"kind": "ocean", "basin": None}

    effective_elev = fields.elevation + fields.silt_depth
    forest = lakes.build_lake_hierarchy(effective_elev, fields.is_ocean, fields.neighbor_idx)
    leaf = _leaf_lakes_by_node(forest).get(node_idx)
    if leaf is None:
        return {"kind": "no_basin", "basin": None}

    rivers = hydrology.group_rivers(fields)
    if fields.lake_depth[node_idx] > hydrology.LAKE_MIN_VISIBLE_DEPTH_M:
        components = _lake_components_sorted(fields)
        lake_id, members = next((i, m) for i, m in enumerate(components) if node_idx in m.tolist())
        lake_obj = _smallest_lake_containing(forest, set(members.tolist())) or leaf
        basin = _lake_basin_summary(fields, lake_obj, lake_id, rivers, is_lake=True)
        return {"kind": "lake", "basin": basin}

    basin = _lake_basin_summary(fields, leaf, None, rivers, is_lake=False)
    return {"kind": "basin", "basin": basin}

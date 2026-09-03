"""FastAPI routes -- a single in-memory world at a time, no world id (see docs/architecture.md)."""

from __future__ import annotations

import base64
import json
import os
import threading
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from scipy.spatial import cKDTree

from . import (
    biomes,
    climate,
    coastline,
    eustasy,
    geodesic,
    geometry,
    hydrology,
    lakes,
    mantle,
    persistence,
    plates,
    projections,
    render_image,
    stats,
    stranded_basins,
)
from .world import DEFAULT_MANTLE_CENTERS, World, generate_world, step_world

# A generous ceiling on requested image dimensions -- width/height come straight from the
# client's query string, and PIL will happily try to allocate whatever it's told, so an
# unbounded value here would let a malformed or malicious request force a huge allocation.
MAX_RENDER_DIMENSION_PX = 4000
# Each animation frame costs a full step_world + render, so this bounds worst-case request
# time for POST /world/animate the same way MAX_RENDER_DIMENSION_PX bounds /world/render's
# own worst case.
MAX_ANIMATION_FRAMES = 240

app = FastAPI(title="mantle-bloom")

# FRONTEND_PORT lets bin/restart.sh --frontend-port point CORS at a non-default frontend
# port; 5174 stays allowed alongside it since that's Vite dev's own fallback when its primary
# port is taken.
_frontend_port = os.environ.get("FRONTEND_PORT", "5173")
app.add_middleware(
    CORSMiddleware,
    allow_origins=list({f"http://localhost:{_frontend_port}", "http://localhost:5174"}),
    allow_methods=["*"],
    allow_headers=["*"],
)

_state: dict[str, World | None] = {"world": None}

# FastAPI runs sync `def` routes in a thread pool, so two requests can otherwise race on the
# same World's mutable numpy arrays -- and, worse, launch Numba's `parallel=True` kernels
# (fluid_dynamics.py, fluid_dynamics_healpix.py, healpix_grid.py) concurrently from two
# threads, which the default `workqueue` threading layer is not safe against and crashes the
# interpreter outright. So every route that runs the simulation compute path or mutates
# shared World state holds this one lock: step / animate / generate / controls / stats /
# render. The read-mostly inspector routes (plates, rivers, lakes, *_at) touch neither a
# Numba kernel nor a shared-array write and are left lock-free.
#
# Long-running writes where an overlap means a genuine conflict (step, animate) reject with a
# 503 rather than blocking -- the second request would otherwise silently hang for as long as
# the first takes (seconds on a large world); see api.ts's stepWorld for the client-side
# retry. The quick or read-side holders (generate, controls, stats, render) just block and
# wait the brief window out -- a render especially is called far too often, and is far too
# cheap, to fail every time it lands during a step (see the render route's own docstring).
_world_lock = threading.Lock()


@contextmanager
def _reject_if_busy(detail: str):
    """Hold `_world_lock` for the block, or 503 immediately if another request already has it
    -- for the long-running writes (step, animate) where an overlap is a real conflict the
    caller should retry, not something to silently queue behind."""
    if not _world_lock.acquire(blocking=False):
        raise HTTPException(status_code=503, detail=detail)
    try:
        yield
    finally:
        _world_lock.release()


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
    # with, relative to elevation_lines.TARGET_LINE_SPACING_KM's own default spacing.
    # Validated against plates.NODE_DENSITY_CHOICES below rather than accepted as an
    # arbitrary float -- there's no continuous "in-between" density the UI offers, only a
    # fixed set of multipliers.
    node_density: float = 1.0
    # The UI's "initial soil maturity" slider (0 to 1) -- optional, falls back to
    # world.generate_world's own default (no seeding, a fully barren starting world -- see
    # geology.seed_initial_soil), but the frontend always sends it.
    initial_soil_maturity: float | None = None
    # The UI's "climate & biome resolution" choice -- how finely climate.py's own grid
    # resolves temperature/wind/humidity/precipitation, and (scaled the same way) how finely
    # the Biome/Combined/Resources/Soil-Quality views' own render grid resolves them.
    # Validated against climate.CLIMATE_DENSITY_CHOICES below, same reasoning node_density's
    # own validation gives.
    climate_density: float = 1.0
    # The UI's "Fluid dynamics resolution" Advanced-settings choice -- how finely Ocean/
    # Atmospheric Fluid Dynamics resolves currents/wind, independent of climate_density (see
    # World.fluid_density's own comment for why). Validated against
    # climate.FLUID_DENSITY_CHOICES below -- a smaller, lower-capped set than
    # climate_density's own, see that constant's own comment for why.
    fluid_density: float = 1.0


class StepRequest(BaseModel):
    years: float


class AnimateRequest(BaseModel):
    # Same shape/defaults as /world/render's own query params (see `render`), plus the two
    # animation-specific knobs below. `rotation` is a 9-comma-separated-float string, same
    # convention _parse_view_rotation already expects.
    projection: str = "behrmann"
    view: str = "elevation"
    width: int = 1100
    height: int = 611
    rotation: str | None = None
    # The UI's "years per frame" choice -- the frontend defaults this to 1_000_000, matching
    # the feature's own framing ("a new frame every million years"), but it's sent explicitly
    # rather than defaulted here so a future UI change doesn't need a matching backend change.
    years_per_frame: float
    # Total frames, including frame 0 (the world's current, unstepped state) -- so this
    # permanently advances the world by (num_frames - 1) * years_per_frame years (see
    # render_image.stream_animation_mp4's own docstring).
    num_frames: int


class ExportHexGridRequest(BaseModel):
    frequency: int = geodesic.DEFAULT_FREQUENCY


class ControlsRequest(BaseModel):
    # All optional and independently settable -- the UI's "Controls" window (see
    # StatsModal.tsx's "Controls" button) sends only whichever control the user actually
    # touched. Omitted fields leave that World attribute untouched. See World.sea_level_m/
    # World.solar_multiplier/World.simulate_plate_movement/World.simulate_climate_biomes/
    # World.wind_model for what each one does physically.
    sea_level_m: float | None = None
    solar_multiplier: float | None = None
    simulate_plate_movement: bool | None = None
    simulate_climate_biomes: bool | None = None
    wind_model: str | None = None


WIND_MODEL_CHOICES = ("cfd", "diagnostic")


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


def _plate_overlaps(world: World) -> dict[int, list[dict]]:
    """For every plate, which other plates its territory currently overlaps and by how much
    (`fraction` = share of *this* plate's own nodes that sit on top of that other plate),
    plus `since_years` -- the earliest `world.elapsed_years` at which any of this plate's
    still-overlapping nodes first went over that partner (`ElevationLine.overlap_onset_years`,
    stamped by merge_split.update_overlap_tracking; None if the save predates the field or the
    overlap only appeared this step). The node-cloud overlap itself is
    `plates.compute_node_overlap`, shared with that tracker so the two can't drift."""
    tol = plates.OVERLAP_TOLERANCE_MULT * plates.line_spacing_rad(world.node_density)
    overlap = plates.compute_node_overlap(world.plates, tol)
    active = [p for p in world.plates if p.node_count() > 0]
    result: dict[int, list[dict]] = {p.plate_id: [] for p in active}
    for src_plate in active:
        info = overlap.get(src_plate.plate_id)
        if info is None or not info["by_partner"]:
            continue
        count = src_plate.node_count()
        mask = info["overlap_mask"]
        onset = src_plate.collect("overlap_onset_years")
        # A per-node partner id isn't tracked, so `since_years` is the earliest non-zero
        # onset over *all* of this plate's overlapping nodes -- close enough for "how long
        # has this plate been stuck on top of something", which is what the number is for.
        overlapping_onsets = onset[mask & (onset > 0.0)]
        since = float(overlapping_onsets.min()) if len(overlapping_onsets) else None
        for partner_id, n in info["by_partner"].items():
            result[src_plate.plate_id].append(
                {
                    "plate_id": partner_id,
                    "fraction": round(n / count, 4),
                    "since_years": since,
                }
            )
    for plate_id in result:
        result[plate_id].sort(key=lambda entry: -entry["fraction"])
    return result


def _plate_summary(plate: plates.Plate, world: World, overlaps: dict[int, list[dict]]) -> dict:
    outline = plate.get_bounding_polygon()
    node_points, node_elevation = plate.all_points_and_elevation()
    ellipse = plates.plate_bounding_ellipse(node_points)
    omega = np.asarray(plate.omega, dtype=float)
    rate_rad = float(np.linalg.norm(omega))
    euler_pole = None
    if rate_rad > 0.0:
        pole_lat, pole_lon = geometry.xyz_to_latlon((omega / rate_rad)[None, :])
        euler_pole = {"lat_deg": round(float(np.degrees(pole_lat[0])), 2), "lon_deg": round(float(np.degrees(pole_lon[0])), 2)}
    submerged = (
        float(np.mean(node_elevation <= world.sea_level_m)) if len(node_elevation) else 0.0
    )
    collisions = [
        {"plate_id": other, "years": years}
        for (a, b), years in world.collision_progress.items()
        for other in ((b,) if a == plate.plate_id else (a,) if b == plate.plate_id else ())
    ]
    return {
        "plate_id": plate.plate_id,
        "crust_type": plate.crust_type,
        # Plate motion (torque.py's dynamic omega): surface speed at the planet's radius, the
        # Euler-rotation pole, and whether the plate is pinned at mantle.MAX_PLATE_RATE (a
        # long-run pathology when it's true for most plates -- see docs/TODO.md).
        "speed_cm_per_yr": round(mantle.rad_per_yr_to_cm_per_yr(rate_rad), 3),
        "at_max_rate": bool(rate_rad >= mantle.MAX_PLATE_RATE - 1e-12),
        "euler_pole": euler_pole,
        "age_steps": plate.age_steps,
        # Median height of the plate's own nodes and the share of them at/below sea level --
        # a continental plate reading mostly-submerged is a red flag for over-stretching.
        "median_elevation_m": round(float(np.median(node_elevation)), 1) if len(node_elevation) else None,
        "submerged_fraction": round(submerged, 3),
        # Other plates this plate's territory currently sits on top of (see _plate_overlaps),
        # and any sustained-collision timers involving it (merge_split.update_collision_progress).
        "overlaps": overlaps.get(plate.plate_id, []),
        "collisions": collisions,
        # A PlateWithLines-only concept (how many of its own rows still have at least one
        # node -- see merge_split.py's own "no_land"/consumption checks, which mirror this
        # same "zero nodes is a real, filtered-out state" reasoning); None for any other
        # representation, which has no notion of "rows" at all.
        "num_rows": (
            sum(1 for line in plate.lines if len(line) > 0) if isinstance(plate, plates.PlateWithLines) else None
        ),
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
    same kind of information, just without live water. `floor_*`/`outlet_*`/`water_elevation_m`/
    `is_spilling` are read directly off `lake` itself (`floor_elevation`/`max_depth`/
    `current_water_elevation`/`is_spilling`, already resolved by `lakes.step_lakes` -- see
    `HydrologyFields.lake_forest`'s own docstring for why this must be the *same* `Lake` object
    that produced `fields.lake_depth`, not a hierarchy rebuilt from a different elevation/silt
    snapshot: the two can genuinely disagree on basin topology, which used to surface as a
    lake's own `outlet_elevation_m` reading *below* its `floor_elevation_m`) -- `outlet_*` is
    `None` for an unresolved closed/endorheic basin with no known spill (a legitimate case, see
    lakes.py's own docstring), not a bug. Rivers "feeding into" this basin are simply every
    `RiverInfo` whose own mouth lands on one of this basin's members, regardless of that
    river's own `mouth_type` label ("lake", or "other" for a river dead-ending in a currently-dry
    sink) -- from this basin's own point of view both are equally "a river that ends here"."""
    members = lake.members

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

    # Bare elevation only picks *which* node represents the floor for the map marker's
    # position -- the reported elevation *value* comes from `lake.floor_elevation` itself
    # (below), not a separate recomputation, so it's guaranteed consistent with `outlet_elev`.
    floor_node = int(display_members[np.argmin(fields.elevation[display_members])])
    outlet_elev = lake.max_depth
    outlet_node = lake.outlet_node_idx

    water_elevation = lake.current_water_elevation if len(wet_members) > 0 else None
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
        "floor_elevation_m": float(lake.floor_elevation),
        "outlet_xyz": None if outlet_elev is None else _round_coords(fields.points[outlet_node]),
        "outlet_elevation_m": outlet_elev,
        "water_elevation_m": None if water_elevation is None else float(water_elevation),
        "is_spilling": bool(is_spilling),
        "inflow_rivers": inflow_rivers,
    }


@app.post("/world/generate")
def generate(req: GenerateRequest) -> dict:
    if req.node_density not in plates.NODE_DENSITY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown node_density {req.node_density!r}; choices are {plates.NODE_DENSITY_CHOICES}")
    if req.climate_density not in climate.CLIMATE_DENSITY_CHOICES:
        raise HTTPException(
            status_code=400, detail=f"unknown climate_density {req.climate_density!r}; choices are {climate.CLIMATE_DENSITY_CHOICES}"
        )
    if req.fluid_density not in climate.FLUID_DENSITY_CHOICES:
        raise HTTPException(
            status_code=400, detail=f"unknown fluid_density {req.fluid_density!r}; choices are {climate.FLUID_DENSITY_CHOICES}"
        )
    with _world_lock:
        world = generate_world(
            req.seed,
            num_plates=req.num_plates,
            continental_fraction=req.continental_fraction,
            land_fraction=req.land_fraction,
            num_mantle_centers=req.num_mantle_centers,
            axial_tilt_deg=req.axial_tilt_deg,
            node_density=req.node_density,
            initial_soil_maturity=req.initial_soil_maturity,
            climate_density=req.climate_density,
            fluid_density=req.fluid_density,
        )
        _state["world"] = world
    return _summary(world)


@app.get("/world/summary")
def get_summary() -> dict:
    """The current world's summary (same shape /world/generate and /world/step return) --
    lets the client re-sync with whatever world is already sitting in server memory (e.g.
    after a browser refresh, see App.tsx's restore-on-mount effect and its view-state cookie)
    without generating a new one. `404` if no world has been generated yet, same as every
    other /world/* endpoint."""
    world = _require_world()
    return _summary(world)


@app.post("/world/step")
def step(req: StepRequest) -> dict:
    """Advances plate tectonics/climate/erosion by `req.years` (see step_world) -- and, gated
    on World.simulate_climate_biomes the same way erosion/hydrology already are, Ocean/
    Atmospheric Fluid Dynamics by their own fixed real-time increment regardless of `req.years`
    (see world.py's `_advance_fluid_dynamics`)."""
    world = _require_world()
    with _reject_if_busy("a step is already in progress"):
        step_world(world, req.years)
        return _summary(world)


@app.get("/world/save")
def save_world() -> Response:
    """The "File > Save World" download -- the *entire* current world (every plate/line,
    mantle field, caches, event log -- see World) pickled as a single opaque file (see
    persistence.py for why pickle, deliberately with no cross-version compatibility
    promise). `404` if no world has been generated yet."""
    world = _require_world()
    body = persistence.save_world_bytes(world)
    filename = f"mantle-bloom-seed{world.seed}-{int(world.elapsed_years)}y.mbworld"
    return Response(
        content=body, media_type="application/octet-stream", headers={"Content-Disposition": f'attachment; filename="{filename}"'}
    )


@app.post("/world/load")
async def load_world(request: Request) -> dict:
    """The "File > Load World" upload -- the raw bytes of a file /world/save previously
    produced, as the request body (not JSON -- see persistence.py). Replaces whatever world
    previously existed, same as /world/generate. `400` if the bytes aren't a valid
    mantle-bloom world file (pickle can raise many different exception types on malformed
    or foreign input, so this is caught broadly)."""
    body = await request.body()
    try:
        world = persistence.load_world_bytes(body)
    except Exception as exc:
        raise HTTPException(status_code=400, detail="invalid or incompatible world file") from exc
    _state["world"] = world
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
    yet.

    Takes `_world_lock` (blocking) around the actual read, unlike step/animate which 503 on a
    busy lock and let the caller retry. A render is a read, called far more often (every view
    switch, every post-step refresh, each animation frame) and each mutation-collecting pass
    it makes (plates.collect_all_points, collect_all_lake_depth, ...) walks world.plates
    independently -- concurrently with an in-flight step actually mutating those same plates
    (deform() grows/shrinks nodes), those passes can see different node counts and desync,
    which crashes with an IndexError rather than just rendering a stale frame. Blocking here
    (confirmed to actually reproduce that IndexError without it, via a deliberately raced
    step+render) waits out that brief window instead of failing every render that happens to
    land during one -- a hard failure on nearly every routine step->refresh overlap would be
    a far worse trade for a request this frequent and this cheap to just wait a moment for."""
    world = _require_world()
    if projection not in projections.PROJECTIONS:
        raise HTTPException(status_code=400, detail=f"unknown projection {projection!r}")
    if view not in render_image.VIEWS:
        raise HTTPException(status_code=400, detail=f"unknown view {view!r}; choices are {render_image.VIEWS}")
    if not (1 <= width <= MAX_RENDER_DIMENSION_PX and 1 <= height <= MAX_RENDER_DIMENSION_PX):
        raise HTTPException(status_code=400, detail=f"width/height must be in [1, {MAX_RENDER_DIMENSION_PX}]")
    view_rotation = _parse_view_rotation(rotation)

    with _world_lock:
        image_base64 = render_image.render_png_base64(world, projection, view, width, height, view_rotation)
    return {
        "projection": projection,
        "elapsed_years": world.elapsed_years,
        "image_base64": image_base64,
    }


@app.post("/world/animate")
def animate(req: AnimateRequest) -> StreamingResponse:
    """The "File > Make Animation" action: streams newline-delimited JSON progress while
    rendering an H.264/MP4 video of `view`/`projection`'s progress -- one frame for the
    world's current state plus `req.num_frames - 1` more, each `req.years_per_frame` further
    along (see render_image.stream_animation_mp4 for the encoding details). Each response
    line is one JSON object: `{"type": "progress", "frame": n, "total": N}` as each frame
    finishes, then a final `{"type": "done", "video_base64": ..., "mime": "video/mp4",
    ...world summary fields}`. If rendering raises partway through, a `{"type": "error",
    "detail": ...}` line is emitted instead -- the HTTP status is already 200 by then, since
    the stream has started.

    **This permanently advances the world** by `(req.num_frames - 1) * req.years_per_frame`
    years, the same as calling /world/step that many times -- not a side-effect-free preview.
    Same validation as /world/render for projection/view/width/height/rotation, plus
    `num_frames` bounded to `[1, MAX_ANIMATION_FRAMES]` (each frame costs a full step_world +
    render). `404` if no world has been generated yet, `503` if a step or another animation
    is already in progress."""
    world = _require_world()
    if req.projection not in projections.PROJECTIONS:
        raise HTTPException(status_code=400, detail=f"unknown projection {req.projection!r}")
    if req.view not in render_image.VIEWS:
        raise HTTPException(status_code=400, detail=f"unknown view {req.view!r}; choices are {render_image.VIEWS}")
    if not (1 <= req.width <= MAX_RENDER_DIMENSION_PX and 1 <= req.height <= MAX_RENDER_DIMENSION_PX):
        raise HTTPException(status_code=400, detail=f"width/height must be in [1, {MAX_RENDER_DIMENSION_PX}]")
    if not (1 <= req.num_frames <= MAX_ANIMATION_FRAMES):
        raise HTTPException(status_code=400, detail=f"num_frames must be in [1, {MAX_ANIMATION_FRAMES}]")
    view_rotation = _parse_view_rotation(req.rotation)

    # Acquire `_world_lock` synchronously (like _reject_if_busy, but the response streams so
    # the release has to happen in the generator's `finally`, not a `with` here) -- an
    # overlapping step/animation still 503s before the stream starts.
    if not _world_lock.acquire(blocking=False):
        raise HTTPException(status_code=503, detail="a step or animation is already in progress")

    def _stream():
        try:
            for message in render_image.stream_animation_mp4(
                world, req.projection, req.view, req.width, req.height, view_rotation, req.years_per_frame, req.num_frames
            ):
                if message[0] == "progress":
                    _, frame, total = message
                    yield json.dumps({"type": "progress", "frame": frame, "total": total}) + "\n"
                else:
                    _, mp4_bytes = message
                    yield json.dumps({
                        "type": "done",
                        "mime": "video/mp4",
                        "video_base64": base64.b64encode(mp4_bytes).decode("ascii"),
                        **_summary(world),
                    }) + "\n"
        except Exception as exc:  # noqa: BLE001 -- status is already 200, so surface it as a data line
            yield json.dumps({"type": "error", "detail": f"{type(exc).__name__}: {exc}"}) + "\n"
        finally:
            _world_lock.release()

    return StreamingResponse(_stream(), media_type="application/x-ndjson")


@app.post("/world/controls")
def set_controls(req: ControlsRequest) -> dict:
    """Live-adjusts sea level, solar heat, the plate-movement/climate-biomes step toggles,
    and/or the wind model on the current world (the "Controls" window -- see World.sea_level_m/
    World.solar_multiplier/World.simulate_plate_movement/World.simulate_climate_biomes/
    World.wind_model) and immediately recomputes world.climate_cache so the very next /world/render or
    /world/stats call reflects the change without waiting for a step -- climate_cache is
    otherwise only refreshed by erosion.py once per step (see World.climate_cache), which
    would make a real-time control feel laggy or unresponsive between steps. This one-time
    recompute happens regardless of the current simulate_climate_biomes setting -- it's a
    single call triggered by a UI action, not a per-step cost, so it's cheap even when the
    user has turned off per-step climate computation to speed up stepping. `404` if no world
    has been generated yet."""
    world = _require_world()
    if req.wind_model is not None and req.wind_model not in WIND_MODEL_CHOICES:
        raise HTTPException(
            status_code=400, detail=f"unknown wind_model {req.wind_model!r}; choices are {WIND_MODEL_CHOICES}"
        )
    with _world_lock:
        if req.sea_level_m is not None:
            # Sea level is eustatic now (see eustasy.py) -- the slider sets the conserved
            # ocean water volume to whatever floats the current hypsometry at this level,
            # rather than pinning sea_level_m directly (which the next step would re-solve away).
            eustasy.set_sea_level_via_water_budget(world, req.sea_level_m)
        if req.solar_multiplier is not None:
            world.solar_multiplier = req.solar_multiplier
        if req.simulate_plate_movement is not None:
            world.simulate_plate_movement = req.simulate_plate_movement
        if req.simulate_climate_biomes is not None:
            world.simulate_climate_biomes = req.simulate_climate_biomes
        if req.wind_model is not None:
            world.wind_model = req.wind_model
        world.climate_cache = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density))
    return {
        "sea_level_m": world.sea_level_m,
        "solar_multiplier": world.solar_multiplier,
        "simulate_plate_movement": world.simulate_plate_movement,
        "simulate_climate_biomes": world.simulate_climate_biomes,
        "wind_model": world.wind_model,
    }


@app.get("/world/plates")
def list_plates() -> dict:
    """Every plate's outline + metadata (row/point counts, bounding ellipse) as JSON, for the
    "Plate Inspector" map mode -- unlike /world/render, the client renders this itself
    interactively rather than receiving a baked PNG. Un-rotated/true-frame throughout (no
    `rotation` param): the client applies its current view rotation only at draw time, same
    philosophy as climate/render-grid geometry (see docs/simulation-model.md#rotating-the-view).
    `404` if no world has been generated yet."""
    world = _require_world()
    overlaps = _plate_overlaps(world)
    return {"elapsed_years": world.elapsed_years, "plates": [_plate_summary(p, world, overlaps) for p in world.plates]}


@app.get("/world/stats")
def get_stats() -> dict:
    """Aggregate physical/temperature/precipitation statistics for the Stats panel -- see
    stats.py. Stateless: recomputed fresh from the current world state every call (the
    frontend is what accumulates a history over time, by calling this after every
    generate/step, same as it already does for /world/render). `404` if no world has been
    generated yet."""
    world = _require_world()
    with _world_lock:
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


@app.get("/world/sample_at")
def sample_at(lat_deg: float, lon_deg: float) -> dict:
    """Everything the "Elevation & Biome" view's click-to-inspect popup shows for one point
    (also offered on the Elevation and Biome views): elevation, biome, precipitation,
    temperature, and which plate owns it. Elevation/biome/precipitation/temperature are read
    from the same climate grid `/world/render` and `/world/stats` already use
    (`climate.compute_climate_cached` -- the world's own `climate_density` resolution),
    picked by nearest grid cell; `temperature_c` is the same land-air / ocean-surface
    composite the temperature view and biome classification use (see
    `climate.compute_climate`). `biome` may differ from the rendered pixel right at a
    coastline -- the Combined/Biome *renders* re-run `biomes.smooth_biome_field` on a finer
    grid (see `render_image._biome_fields`) -- but this is the canonical per-cell field every
    other consumer reads, and it carries the same boundary-cleanup pass. `plate_id` is the same nearest-node hit-test as `/world/plate_at`. The
    client unprojects its click through the active view rotation to a true lat/lon first,
    same as `/world/plate_at`, so this endpoint never needs to know about rotation. `400`
    for non-finite input, `404` if no world has been generated yet."""
    world = _require_world()
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    with _world_lock:
        fields = climate.compute_climate_cached(world)
        query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
        plate_id = plates.nearest_plate_id(world.plates, query_xyz)
    height = fields.lat_deg.shape[0]
    width = fields.lon_deg.shape[0]
    # Invert climate._build_grid's cell centering: row 0 spans the north pole, column 0 the
    # antimeridian, columns wrap.
    row = int(np.clip(round((90.0 - lat_deg) * height / 180.0 - 0.5), 0, height - 1))
    col = int(round((lon_deg + 180.0) * width / 360.0 - 0.5)) % width
    biome_id = int(fields.biome_ids[row, col])
    is_ocean = bool(fields.is_ocean[row, col])
    temperature_c = float(
        fields.ocean_temperature_c[row, col] if is_ocean else fields.air_temperature_c[row, col]
    )
    return {
        "lat_deg": lat_deg,
        "lon_deg": lon_deg,
        "elevation_m": float(fields.elevation_m[row, col]),
        "is_ocean": is_ocean,
        "biome_id": biome_id,
        "biome": biomes.BIOME_NAMES[biome_id],
        "temperature_c": temperature_c,
        "precipitation_mm": float(fields.precipitation_mm[row, col]),
        "plate_id": plate_id,
    }


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
    this codebase -- resolved back up to its own basin in `fields.lake_forest`, this step's
    already-resolved depression hierarchy (see `HydrologyFields.lake_forest`'s own docstring
    for why this must be the *same* tree `lakes.step_lakes` used, not a separately rebuilt
    one), so a lake that's currently the result of several original catchments merging
    together reports its own *current* outlet, not some interior saddle that's already
    flooded. Regrouped fresh from `world.hydrology_cache` on
    every call, same as `/world/rivers` -- `lake_id` is only meaningful against this same
    response, not a stable identity across steps. `coastline_segments` included for the same
    reason `/world/rivers` includes it: this view draws no filled backdrop of its own. Both are
    `[]` before the first step. `404` if no world has been generated yet."""
    world = _require_world()
    fields = world.hydrology_cache
    if fields is None:
        return {"elapsed_years": world.elapsed_years, "lakes": [], "coastline_segments": []}
    forest = fields.lake_forest
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

    forest = fields.lake_forest
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


@app.get("/world/stranded_basins")
def list_stranded_basins() -> dict:
    """Endorheic depressions whose floor sits below sea level and that have no drainage path
    to the ocean at all -- the "land-locked coastal pit" pathology from docs/TODO.md's
    coastal-speckle section, which the event log records but drowns in near-sea-level
    transient-pond churn (see docs/debugging.md). Read straight off this step's already-
    resolved depression hierarchy (`world.hydrology_cache.lake_forest`): a top-level basin
    with no known spill (`Lake.max_depth is None`) whose floor is below `world.sea_level_m`.
    `persisted_years`/`steps_seen` come from `world.stranded_basin_tracks`, the cross-step
    centroid-matched tracker `world.step_world` maintains (at most one step stale, same as
    `hydrology_cache` itself). `stranded_basins` is `[]` before the first step. `404` if no
    world has been generated yet."""
    world = _require_world()
    fields = world.hydrology_cache
    if fields is None or len(fields.points) == 0:
        return {"elapsed_years": world.elapsed_years, "sea_level_m": world.sea_level_m, "stranded_basins": []}
    basins = stranded_basins.find_stranded_basins(
        fields.lake_forest, fields.elevation, fields.points, fields.lake_depth, world.sea_level_m
    )
    stranded_basins.enrich_with_persistence(basins, world.stranded_basin_tracks, world.elapsed_years)
    return {
        "elapsed_years": world.elapsed_years,
        "sea_level_m": world.sea_level_m,
        "stranded_basins": [_stranded_basin_summary(b) for b in basins],
    }


def _stranded_basin_summary(basin: stranded_basins.StrandedBasin) -> dict:
    return {
        "floor_elevation_m": round(basin.floor_elevation_m, 1),
        "depth_below_sea_level_m": round(basin.depth_below_sea_level_m, 1),
        "catchment_node_count": basin.catchment_node_count,
        "flooded_node_count": basin.flooded_node_count,
        "water_elevation_m": None if basin.water_elevation_m is None else round(basin.water_elevation_m, 1),
        "centroid_xyz": [round(c, _COORD_DECIMALS) for c in basin.centroid_xyz],
        "centroid_lat_deg": basin.centroid_lat_deg,
        "centroid_lon_deg": basin.centroid_lon_deg,
        "floor_xyz": [round(c, _COORD_DECIMALS) for c in basin.floor_xyz],
        "first_seen_years": basin.first_seen_years,
        "persisted_years": basin.persisted_years,
        "steps_seen": basin.steps_seen,
    }


@app.post("/world/export_hexgrid")
def export_hexgrid(req: ExportHexGridRequest) -> dict:
    """The "File > Export Hex Grid" action: tiles the sphere into a geodesic-icosahedron
    hex/pentagon dome (independent of the plate simulation's own node cloud -- see
    geodesic.py) at `req.frequency`, samples the current world's elevation/ocean/biome onto
    each tile, and returns the whole tiling as JSON for an external application (see
    docs/hex-export-format.md for the file format and the neighbor-finding pseudocode its
    own `neighbor_ids`/`corner_vertex_ids` fields back). `400` if `frequency` isn't one of
    `geodesic.FREQUENCY_CHOICES`, `404` if no world has been generated yet."""
    world = _require_world()
    if req.frequency not in geodesic.FREQUENCY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown frequency {req.frequency!r}; choices are {geodesic.FREQUENCY_CHOICES}")
    return geodesic.export_hexgrid(world, req.frequency)


# --- Static frontend -----------------------------------------------------------------------
# In the packaged single-process app there is no separate Vite server: FastAPI serves the
# built React bundle itself, from the same origin as the API (so api.ts's API_BASE is ""). The
# dev workflow (bin/restart.sh) is unaffected -- it runs its own frontend server and never
# sets MANTLE_BLOOM_FRONTEND_DIST.


def mount_frontend(app: FastAPI, dist_dir: str | os.PathLike[str]) -> None:
    """Serve a `vite build` output directory (its `index.html` + `assets/`) at `/`.

    Mounted last, after every API route, so `/world/*` and `/docs` keep priority; anything
    else falls through to the bundle, with unknown non-API paths served `index.html` so a
    hard refresh on a client-side view doesn't 404."""
    dist = Path(dist_dir)
    index = dist / "index.html"
    if not index.is_file():
        raise RuntimeError(f"no index.html under {dist!r} -- run `npm run build` in frontend/ first")

    app.mount("/assets", StaticFiles(directory=dist / "assets"), name="assets")

    @app.get("/", include_in_schema=False)
    def _index() -> FileResponse:
        return FileResponse(index)

    @app.get("/{path:path}", include_in_schema=False)
    def _spa_fallback(path: str) -> FileResponse:
        candidate = (dist / path).resolve()
        if dist.resolve() in candidate.parents and candidate.is_file():
            return FileResponse(candidate)
        return FileResponse(index)


_frontend_dist = os.environ.get("MANTLE_BLOOM_FRONTEND_DIST")
if _frontend_dist:
    mount_frontend(app, _frontend_dist)

"""V2 FastAPI app -- same route shapes as `main.py`, mounted at `/v2` (see `main.py`'s
`app.mount("/v2", v2_app)`). Reuses `main.py`'s own request/response helper functions
directly (`_summary`, `_plate_summary`, `_river_summary`, `_coastline_segments_json`,
`_lake_basin_summary`, `_parse_view_rotation`, ...) since every one of them is already
generic over `World`/`Plate` -- see world_v2.py's own docstring for why v2 populates the same
`World` dataclass v1 does. A separate module-level `_state` dict keeps a v2 world entirely
independent of whatever v1 world is loaded in the same process (see the plan's frontend
toggle section).

Not supported here: the two ocean-sediment map views (`render_image.FLUID_VIEWS` --
oceanCfdSediment/oceanCfdDeposition) -- v2's CFD state is a flat `(npix,)` HEALPix array, not
the `(H, W)` grid that pair's rendering code assumes, and sediment has no HEALPix port yet.
Every other view (elevation/plates/platesDetail/biome/combined/climate/resources/soilQuality)
works unchanged; `CLIMATE_VIEWS` (temperature/wind/oceanCurrents/humidity/precipitation/biome)
in particular render natively off the HEALPix grid for a v2 world (see
render_image._render_climate_view_healpix), not through `resample_uv_to_equirect`'s (H, W)
bridge -- that seam still exists (see atmosphere_cfd_v2.py) but is now only used by other,
still-equirectangular callers (e.g. erosion.py/hydrology.py's own grids).
"""

from __future__ import annotations

import os
import threading

import numpy as np
from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from scipy.spatial import cKDTree

from .. import climate, geodesic, geometry, hydrology, persistence, plates, projections, render_image, stats
from ..elevation_lines import NODE_DENSITY_CHOICES
from ..main import (
    MAX_ANIMATION_FRAMES,
    MAX_RENDER_DIMENSION_PX,
    _coastline_segments_json,
    _lake_basin_summary,
    _lake_components_sorted,
    _leaf_lakes_by_node,
    _parse_view_rotation,
    _plate_summary,
    _river_summary,
    _smallest_lake_containing,
    _summary,
)
from ..world import DEFAULT_MANTLE_CENTERS, World
from .world_v2 import generate_world_v2, step_world_v2

app = FastAPI(title="mantle-bloom-v2")

_frontend_port = os.environ.get("FRONTEND_PORT", "5173")
app.add_middleware(
    CORSMiddleware,
    allow_origins=list({f"http://localhost:{_frontend_port}", "http://localhost:5174"}),
    allow_methods=["*"],
    allow_headers=["*"],
)

_state: dict[str, World | None] = {"world": None}
_step_lock = threading.Lock()

FLUID_DENSITY_CHOICES = (0.5, 1.0, 2.0)  # healpix_grid.NSIDE_CHOICES' own keys


class GenerateRequest(BaseModel):
    seed: int = 0
    num_plates: int | None = None
    continental_fraction: float | None = None
    land_fraction: float | None = None
    axial_tilt_deg: float | None = None
    num_mantle_centers: int = DEFAULT_MANTLE_CENTERS
    node_density: float = 1.0
    initial_soil_maturity: float | None = None
    climate_density: float = 1.0
    fluid_density: float = 1.0


class StepRequest(BaseModel):
    years: float


class AnimateRequest(BaseModel):
    projection: str = "behrmann"
    view: str = "elevation"
    width: int = 1100
    height: int = 611
    rotation: str | None = None
    years_per_frame: float
    num_frames: int


class ExportHexGridRequest(BaseModel):
    frequency: int = geodesic.DEFAULT_FREQUENCY


class ControlsRequest(BaseModel):
    sea_level_m: float | None = None
    solar_multiplier: float | None = None
    simulate_plate_movement: bool | None = None
    simulate_climate_biomes: bool | None = None


def _require_world() -> World:
    world = _state["world"]
    if world is None:
        raise HTTPException(status_code=404, detail="no world generated yet")
    return world


@app.post("/world/generate")
def generate(req: GenerateRequest) -> dict:
    if req.node_density not in NODE_DENSITY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown node_density {req.node_density!r}; choices are {NODE_DENSITY_CHOICES}")
    if req.climate_density not in NODE_DENSITY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown climate_density {req.climate_density!r}; choices are {NODE_DENSITY_CHOICES}")
    if req.fluid_density not in FLUID_DENSITY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown fluid_density {req.fluid_density!r}; choices are {FLUID_DENSITY_CHOICES}")
    world = generate_world_v2(
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
    return _summary(_require_world())


@app.post("/world/step")
def step(req: StepRequest) -> dict:
    world = _require_world()
    if not _step_lock.acquire(blocking=False):
        raise HTTPException(status_code=503, detail="a step is already in progress")
    try:
        step_world_v2(world, req.years)
        return _summary(world)
    finally:
        _step_lock.release()


@app.get("/world/save")
def save_world() -> Response:
    world = _require_world()
    body = persistence.save_world_bytes(world)
    filename = f"mantle-bloom-v2-seed{world.seed}-{int(world.elapsed_years)}y.mbworld"
    return Response(
        content=body, media_type="application/octet-stream", headers={"Content-Disposition": f'attachment; filename="{filename}"'}
    )


@app.post("/world/load")
async def load_world(request: Request) -> dict:
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
    world = _require_world()
    if projection not in projections.PROJECTIONS:
        raise HTTPException(status_code=400, detail=f"unknown projection {projection!r}")
    if view not in render_image.VIEWS:
        raise HTTPException(status_code=400, detail=f"unknown view {view!r}; choices are {render_image.VIEWS}")
    if view in render_image.FLUID_VIEWS:
        raise HTTPException(status_code=400, detail=f"view {view!r} is not available on v2 -- its CFD state is HEALPix-native, not a (H, W) grid")
    if not (1 <= width <= MAX_RENDER_DIMENSION_PX and 1 <= height <= MAX_RENDER_DIMENSION_PX):
        raise HTTPException(status_code=400, detail=f"width/height must be in [1, {MAX_RENDER_DIMENSION_PX}]")
    view_rotation = _parse_view_rotation(rotation)

    return {
        "projection": projection,
        "elapsed_years": world.elapsed_years,
        "image_base64": render_image.render_png_base64(world, projection, view, width, height, view_rotation),
    }


@app.post("/world/animate")
def animate(req: AnimateRequest) -> dict:
    world = _require_world()
    if req.projection not in projections.PROJECTIONS:
        raise HTTPException(status_code=400, detail=f"unknown projection {req.projection!r}")
    if req.view not in render_image.VIEWS:
        raise HTTPException(status_code=400, detail=f"unknown view {req.view!r}; choices are {render_image.VIEWS}")
    if req.view in render_image.FLUID_VIEWS:
        raise HTTPException(status_code=400, detail=f"view {req.view!r} is not available on v2")
    if not (1 <= req.width <= MAX_RENDER_DIMENSION_PX and 1 <= req.height <= MAX_RENDER_DIMENSION_PX):
        raise HTTPException(status_code=400, detail=f"width/height must be in [1, {MAX_RENDER_DIMENSION_PX}]")
    if not (1 <= req.num_frames <= MAX_ANIMATION_FRAMES):
        raise HTTPException(status_code=400, detail=f"num_frames must be in [1, {MAX_ANIMATION_FRAMES}]")
    view_rotation = _parse_view_rotation(req.rotation)

    def _step(w: World, years: float) -> None:
        step_world_v2(w, years)

    image_base64 = render_image.render_animation_gif_base64(
        world, req.projection, req.view, req.width, req.height, view_rotation, req.years_per_frame, req.num_frames, step_fn=_step
    )
    return {**_summary(world), "image_base64": image_base64}


@app.post("/world/controls")
def set_controls(req: ControlsRequest) -> dict:
    world = _require_world()
    if req.sea_level_m is not None:
        world.sea_level_m = req.sea_level_m
    if req.solar_multiplier is not None:
        world.solar_multiplier = req.solar_multiplier
    if req.simulate_plate_movement is not None:
        world.simulate_plate_movement = req.simulate_plate_movement
    if req.simulate_climate_biomes is not None:
        world.simulate_climate_biomes = req.simulate_climate_biomes
    world.climate_cache = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density))
    return {
        "sea_level_m": world.sea_level_m,
        "solar_multiplier": world.solar_multiplier,
        "simulate_plate_movement": world.simulate_plate_movement,
        "simulate_climate_biomes": world.simulate_climate_biomes,
    }


@app.get("/world/plates")
def list_plates() -> dict:
    world = _require_world()
    return {"elapsed_years": world.elapsed_years, "plates": [_plate_summary(p) for p in world.plates]}


@app.get("/world/stats")
def get_stats() -> dict:
    return stats.compute_stats(_require_world())


@app.get("/world/plate_at")
def plate_at(lat_deg: float, lon_deg: float) -> dict:
    world = _require_world()
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    return {"plate_id": plates.nearest_plate_id(world.plates, query_xyz)}


@app.get("/world/rivers")
def list_rivers() -> dict:
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
            continue
        result.append(_lake_basin_summary(fields, lake_obj, lake_id, rivers, is_lake=True))
    return {"elapsed_years": world.elapsed_years, "lakes": result, "coastline_segments": _coastline_segments_json(world)}


@app.get("/world/lake_at")
def lake_at(lat_deg: float, lon_deg: float) -> dict:
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


@app.post("/world/export_hexgrid")
def export_hexgrid(req: ExportHexGridRequest) -> dict:
    world = _require_world()
    if req.frequency not in geodesic.FREQUENCY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown frequency {req.frequency!r}; choices are {geodesic.FREQUENCY_CHOICES}")
    return geodesic.export_hexgrid(world, req.frequency)

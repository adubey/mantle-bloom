"""FastAPI routes. Single in-memory world for v1 (see docs/architecture.md) -- no world id,
no persistence, one world at a time, matching the "elevation view only" v1 scope.
"""

from __future__ import annotations

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from . import geometry, plates, projections, render_image
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


class StepRequest(BaseModel):
    years: float


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
        "outline": outline.tolist(),
        "bounding_ellipse": None
        if ellipse is None
        else {
            "center_xyz": ellipse.center_xyz.tolist(),
            "diameter_a_km": ellipse.diameter_a_km,
            "diameter_b_km": ellipse.diameter_b_km,
            "outline": ellipse.outline_xyz.tolist(),
        },
    }


@app.post("/world/generate")
def generate(req: GenerateRequest) -> dict:
    world = generate_world(
        req.seed,
        num_plates=req.num_plates,
        continental_fraction=req.continental_fraction,
        land_fraction=req.land_fraction,
        num_mantle_centers=req.num_mantle_centers,
        axial_tilt_deg=req.axial_tilt_deg,
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

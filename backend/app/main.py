"""FastAPI routes. Single in-memory world for v1 (see docs/architecture.md) -- no world id,
no persistence, one world at a time, matching the "elevation view only" v1 scope.
"""

from __future__ import annotations

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from . import geometry, mantle, projections
from .world import DEFAULT_MANTLE_CENTERS, World, generate_world, step_world

# Baseline angular length (radians) of a plate's velocity arrow on the map, before scaling
# by how fast that plate is actually moving relative to the fastest allowed rate.
ARROW_BASE_ANGULAR_LENGTH_RAD = 0.15

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
    num_mantle_centers: int = DEFAULT_MANTLE_CENTERS


class StepRequest(BaseModel):
    years: float


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
    }


@app.post("/world/generate")
def generate(req: GenerateRequest) -> dict:
    world = generate_world(req.seed, num_plates=req.num_plates, num_mantle_centers=req.num_mantle_centers)
    _state["world"] = world
    return _summary(world)


@app.post("/world/step")
def step(req: StepRequest) -> dict:
    world = _require_world()
    step_world(world, req.years)
    return _summary(world)


def _project_points(projection: str, world_pts: np.ndarray) -> list[list[float]]:
    lat, lon = geometry.xyz_to_latlon(world_pts)
    x, y = projections.project(projection, lat, lon)
    return np.stack([x, y], axis=-1).tolist()


def _plate_tectonics(projection: str, plate) -> dict:
    """Pole marker, velocity arrow, and boundary outline for the "Plates" map view --
    everything except the elevation lines themselves."""
    seed_xyz = plate.seed_world
    speed = float(np.linalg.norm(plate.omega))

    pole = None
    arrow = None
    if speed > 1e-15:
        pole_xyz = plate.omega / speed
        pole = _project_points(projection, pole_xyz[None, :])[0]

        direction = np.cross(plate.omega, seed_xyz) / speed
        intensity = np.clip(speed / mantle.MAX_PLATE_RATE, 0.3, 1.0)
        arrow_len = ARROW_BASE_ANGULAR_LENGTH_RAD * intensity
        end_xyz = np.cos(arrow_len) * seed_xyz + np.sin(arrow_len) * direction
        end_xyz = end_xyz / np.linalg.norm(end_xyz)
        start, end = _project_points(projection, np.stack([seed_xyz, end_xyz]))
        arrow = {"start": start, "end": end}

    outline_world = plate.outline_world()
    boundary = _project_points(projection, outline_world) if len(outline_world) > 0 else []

    return {
        "pole": pole,
        "rotation_rate_deg_per_myr": np.degrees(speed) * 1e6,
        "velocity_arrow": arrow,
        "boundary": boundary,
    }


@app.get("/world/render")
def render(projection: str = "behrmann") -> dict:
    world = _require_world()
    if projection not in projections.PROJECTIONS:
        raise HTTPException(status_code=400, detail=f"unknown projection {projection!r}")

    plates_out = []
    for plate in world.plates:
        lines_out = []
        for line in plate.lines:
            if len(line.theta) == 0:
                continue
            world_pts = line.world_xyz(plate.frame)
            lines_out.append(
                {
                    "points": _project_points(projection, world_pts),
                    "elevation": line.elevation.tolist(),
                }
            )
        plates_out.append(
            {
                "plate_id": plate.plate_id,
                "crust_type": plate.crust_type,
                "lines": lines_out,
                **_plate_tectonics(projection, plate),
            }
        )

    return {"projection": projection, "elapsed_years": world.elapsed_years, "plates": plates_out}

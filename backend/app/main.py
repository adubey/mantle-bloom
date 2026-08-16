"""FastAPI routes. Single in-memory world for v1 (see docs/architecture.md) -- no world id,
no persistence, one world at a time, matching the "elevation view only" v1 scope.
"""

from __future__ import annotations

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from . import geometry, projections
from .world import DEFAULT_MANTLE_CENTERS, DEFAULT_NUM_PLATES, World, generate_world, step_world

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
    num_plates: int = DEFAULT_NUM_PLATES
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
            lat, lon = geometry.xyz_to_latlon(world_pts)
            x, y = projections.project(projection, lat, lon)
            lines_out.append(
                {
                    "points": np.stack([x, y], axis=-1).tolist(),
                    "elevation": line.elevation.tolist(),
                }
            )
        plates_out.append(
            {
                "plate_id": plate.plate_id,
                "crust_type": plate.crust_type,
                "lines": lines_out,
            }
        )

    return {"projection": projection, "elapsed_years": world.elapsed_years, "plates": plates_out}

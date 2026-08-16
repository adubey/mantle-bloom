"""FastAPI routes. Single in-memory world for v1 (see docs/architecture.md) -- no world id,
no persistence, one world at a time, matching the "elevation view only" v1 scope.
"""

from __future__ import annotations

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from scipy.spatial import cKDTree

from . import geometry, mantle, plates, projections
from .world import DEFAULT_MANTLE_CENTERS, World, generate_world, step_world

# Baseline angular length (radians) of a plate's velocity arrow on the map, before scaling
# by how fast that plate is actually moving relative to the fastest allowed rate.
ARROW_BASE_ANGULAR_LENGTH_RAD = 0.15
# Resolution of the render grid (see _render_grid) -- the same physical spacing the
# simulation itself uses, but swept on a plate-independent global grid so the rendered map's
# coverage never depends on how sparse any one plate's own line data looks once projected.
GRID_SPACING_RAD = plates.TARGET_LINE_SPACING_RAD

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
    # The UI's continents slider -- optional (falls back to an independent per-plate coin
    # flip, see plates.generate_plates) but the frontend always sends one.
    num_continents: int | None = None
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
        # Full current log (bounded, see world.MAX_EVENT_LOG_LENGTH) every time -- simpler
        # for the frontend than diffing "new since last call", and small enough not to
        # matter on the wire.
        "events": [{"elapsed_years": e, "message": m} for e, m in world.events],
    }


@app.post("/world/generate")
def generate(req: GenerateRequest) -> dict:
    world = generate_world(
        req.seed,
        num_plates=req.num_plates,
        num_continents=req.num_continents,
        num_mantle_centers=req.num_mantle_centers,
    )
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


def _collect_all_points(world: World) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Every plate's current elevation-node positions, elevations, and owning plate_id,
    concatenated -- the source data _render_grid resamples from."""
    points_list, elevation_list, owner_list = [], [], []
    for plate in world.plates:
        pts, elev = plate.all_points_and_elevation()
        if len(pts) == 0:
            continue
        points_list.append(pts)
        elevation_list.append(elev)
        owner_list.append(np.full(len(pts), plate.plate_id))
    if not points_list:
        return None
    return (
        np.concatenate(points_list, axis=0),
        np.concatenate(elevation_list, axis=0),
        np.concatenate(owner_list, axis=0),
    )


def _row_cell_half_extent(projection: str, phi: float, dtheta: float) -> tuple[float, float]:
    """How far apart (in *projected* units) this row's grid points end up, in each
    direction -- half the theta-neighbor distance and half the phi-neighbor distance, both
    measured directly by projecting nearby samples rather than assumed. A fixed on-screen
    dot size can't work here: projected spacing isn't uniform (Behrmann, for instance,
    stretches longitude spacing near the poles while compressing latitude spacing), so a
    dot sized to look right at the equator leaves gaps near the poles, and one sized to
    close those gaps grossly overlaps everywhere else. Sizing each row's cells from the
    projection's actual local behavior is what makes the render grid tile with no gaps
    anywhere on the map. theta=0 is used for the phi-neighbor sample so only lat moves (in
    both current projections, longitude-only movement at fixed latitude doesn't perturb the
    projected y coordinate, and vice versa) -- an exact, not approximate, isolation of each
    derivative."""
    origin = geometry.local_xyz(np.array([phi]), np.array([0.0]))
    theta_neighbor = geometry.local_xyz(np.array([phi]), np.array([dtheta]))
    phi_neighbor = geometry.local_xyz(np.array([phi + GRID_SPACING_RAD]), np.array([0.0]))

    (ox, oy), (tx, ty), (px, py) = (
        _project_points(projection, origin)[0],
        _project_points(projection, theta_neighbor)[0],
        _project_points(projection, phi_neighbor)[0],
    )
    half_width = abs(tx - ox) / 2
    half_height = abs(py - oy) / 2
    return half_width, half_height


def _render_grid(world: World, projection: str) -> dict:
    """A uniform lat/lon grid covering the whole sphere (GRID_SPACING_RAD, independent of
    any plate's own line spacing), each cell assigned its nearest elevation node's
    elevation and owning plate. The raw per-plate line data (plates_out below) is spaced at
    the same physical resolution, but only within each plate's own footprint and along its
    own local frame -- once projected, that leaves visible gaps wherever a plate's lines
    don't happen to align densely with the screen, especially since projected spacing isn't
    uniform (e.g. Behrmann stretches it near the poles). A grid swept independently of any
    single plate's data, with every cell given its closest sample regardless of distance,
    covers the whole map with no such gaps -- as long as each cell is actually drawn at the
    right size (cell_half_width/cell_half_height, see _row_cell_half_extent) rather than a
    fixed dot size."""
    collected = _collect_all_points(world)
    if collected is None:
        return {
            "points": [],
            "elevation": [],
            "plate_id": [],
            "crust_type": [],
            "cell_half_width": [],
            "cell_half_height": [],
        }
    all_points, all_elevation, all_owner = collected
    tree = cKDTree(all_points)
    crust_type_by_id = {p.plate_id: p.crust_type for p in world.plates}

    point_chunks = []
    elevation_chunks = []
    owner_chunks = []
    half_width_chunks = []
    half_height_chunks = []
    for phi, theta_candidates, world_pts in plates.iter_local_lattice(np.eye(3), spacing_rad=GRID_SPACING_RAD):
        _, idx = tree.query(world_pts)
        point_chunks.append(world_pts)
        elevation_chunks.append(all_elevation[idx])
        owner_chunks.append(all_owner[idx])

        dtheta = GRID_SPACING_RAD / max(np.cos(phi), 1e-3)
        half_w, half_h = _row_cell_half_extent(projection, phi, dtheta)
        half_width_chunks.append(np.full(len(theta_candidates), half_w))
        half_height_chunks.append(np.full(len(theta_candidates), half_h))

    grid_points = np.concatenate(point_chunks, axis=0)
    grid_elevation = np.concatenate(elevation_chunks, axis=0)
    grid_plate_id = np.concatenate(owner_chunks, axis=0)
    grid_crust_type = [crust_type_by_id[int(i)] for i in grid_plate_id]

    return {
        "points": _project_points(projection, grid_points),
        "elevation": grid_elevation.tolist(),
        "plate_id": grid_plate_id.tolist(),
        "crust_type": grid_crust_type,
        "cell_half_width": np.concatenate(half_width_chunks).tolist(),
        "cell_half_height": np.concatenate(half_height_chunks).tolist(),
    }


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

    return {
        "projection": projection,
        "elapsed_years": world.elapsed_years,
        "plates": plates_out,
        "grid": _render_grid(world, projection),
    }

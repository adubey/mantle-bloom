"""FastAPI routes -- several worlds in one process, each addressed by an
explicit world ID (see docs/architecture.md, #563)."""

from __future__ import annotations

import base64
import json
import logging
import os
import threading
from contextlib import asynccontextmanager
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, RedirectResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from scipy.spatial import cKDTree

from . import (
    auth,
    checkpoint,
    cross_section as cross_section_view,
    debug_worlds,
    faults,
    geodesic,
    geometry,
    healpix_grid,
    mantle,
    persistence,
    plates,
    projections,
    public_view,
    render_image,
    replay_log,
    run_control,
    stats,
    world_names,
    worldsketch,
)
from .world_registry import FairComputeGate, WorldLimitReached, WorldRegistry, WorldSession
from .hydroclimate import biomes, climate, coastline, eustasy, hydrology, lakes, stranded_basins
from .elevation_lines import NODE_DENSITY_CHOICES
from .world import (
    DEFAULT_MANTLE_CENTERS,
    TUNING_MULTIPLIER_FIELDS,
    World,
    generate_world_progress,
    step_world_progress,
)

# A generous ceiling on requested image dimensions -- width/height come straight from the
# client's query string, and PIL will happily try to allocate whatever it's told, so an
# unbounded value here would let a malformed or malicious request force a huge allocation.
MAX_RENDER_DIMENSION_PX = 4000
# Each animation frame costs a full step_world + render, so this bounds worst-case request
# time for POST /worlds/{world_id}/animate the same way MAX_RENDER_DIMENSION_PX bounds /worlds/{world_id}/render's
# own worst case. Also doubles as the safety ceiling for "keep going until Stop is pressed"
# (see _animation_stop_event below) -- that mode is meant to be ended by the user, not by
# this cap, but the cap still applies underneath so a forgotten run can't loop forever.
MAX_ANIMATION_FRAMES = 480
# Ceilings on one step's years for /worlds/{world_id}/step, /worlds/{world_id}/animate and /worlds/{world_id}/run
# (#557): the UI's largest "Years per step" choice is 10 Myr and debug replays
# step 10 Myr at a time, so 100 Myr leaves 10x headroom while stopping a single
# request from asking for an unbounded advance. Positivity/finiteness reuse
# checkpoint.valid_step_years, the same rule /worlds/{world_id}/run already applies.
MAX_STEP_YEARS = 100_000_000.0
# The animation dialog offers steps-per-frame 1/10/100, so this matches the largest
# offered choice exactly: a single frame can hold the compute gate for at most 100
# real steps (down from 1000), bounding how long one animation frame can delay every
# other world -- and the whole run to the dialog's own maximum of 479 x 100 = 47,900
# steps (down from 479,000). Any lower total would reject combinations the dialog
# itself offers, so the per-frame bound plus the per-frame gate release carry the
# fairness guarantee instead of a separate total.
MAX_STEPS_PER_FRAME = 100



@asynccontextmanager
async def _lifespan(_app: FastAPI):
    _restore_checkpoint()
    yield
    for session in _registry.list():
        if session.runner is not None:
            session.runner.close()
        # The writer is a daemon thread: give a record of the last rare event time to land.
        if session.replay_log is not None and not session.replay_log.flush(timeout_s=REPLAY_LOG_SHUTDOWN_FLUSH_S):
            _logger.warning("replay log: shutting down with records still unwritten")


app = FastAPI(title="mantle-bloom", lifespan=_lifespan)

# FRONTEND_PORT lets bin/restart.sh --frontend-port point CORS at a non-default frontend
# port; 5174 stays allowed alongside it since that's Vite dev's own fallback when its primary
# port is taken.
_frontend_port = os.environ.get("FRONTEND_PORT", "5173")
_frontend_origins = frozenset({f"http://localhost:{_frontend_port}", "http://localhost:5174"})
# Added before CORSMiddleware so CORS stays outermost and its headers reach the auth guard's
# 401/403 refusals too (#372). The Vite dev origin counts as same-site for the CSRF check.
app.add_middleware(auth.AuthMiddleware, extra_origins=_frontend_origins)
app.include_router(auth.router)


@app.get("/demo/", include_in_schema=False)
@app.get("/demo", include_in_schema=False)
def demo_login_page() -> RedirectResponse:
    """Show the login page in a local demo, even though local use skips authentication."""
    return RedirectResponse("/login?demo=true", status_code=302)


app.add_middleware(
    CORSMiddleware,
    allow_origins=list(_frontend_origins),
    allow_methods=["*"],
    allow_headers=["*"],
)

_logger = logging.getLogger(__name__)


def _max_worlds() -> int:
    raw = os.environ.get("MANTLE_BLOOM_MAX_WORLDS", "").strip()
    if not raw:
        return 4
    try:
        value = int(raw)
    except ValueError:
        _logger.warning("MANTLE_BLOOM_MAX_WORLDS=%r is not an integer; using 4", raw)
        return 4
    if value < 1:
        _logger.warning("MANTLE_BLOOM_MAX_WORLDS=%r must be at least 1; using 4", raw)
        return 4
    return value


# Every world in this process (#563). See world_registry.py.
_registry = WorldRegistry(max_worlds=_max_worlds())

# FastAPI runs sync `def` routes in a thread pool, so two requests can otherwise race on the
# same World's mutable numpy arrays -- and, worse, launch Numba's `parallel=True` kernels
# (fluid_dynamics.py, fluid_dynamics_healpix.py, healpix_grid.py) concurrently from two
# threads, which the default `workqueue` threading layer is not safe against and crashes the
# interpreter outright. Two locks now divide that old `_world_lock`'s job:
#
# - Each WorldSession's own lock protects that world's mutable state: replacement, mutation
#   and stable reads. The /worlds/{id}/plates inspector holds it because overlap masks and
#   per-node fields must come from one stable topology. The other read-mostly inspector
#   routes (rivers, lakes, *_at) touch neither a Numba kernel nor a shared-array write and
#   only snapshot the session's world reference.
# - `_compute_gate` serializes the simulation compute path (generate / step / animate, and
#   kernel-backed reads such as render / controls / stats) across *all* worlds, FIFO, so
#   independently running worlds take turns instead of entering Numba concurrently.
#
# Long-running writes where an overlap on the *same* world is a genuine conflict (step,
# animate) still reject with a 503 rather than blocking (session.begin_operation); an
# operation on a different world waits at the compute gate instead. Quick read-side holders
# just wait their brief window out.
_compute_gate = FairComputeGate()


class SketchRequest(BaseModel):
    # A PNG (see worldsketch.py's module docstring for the ink convention: near-black strokes
    # = coastline, blue-ish = rivers, brown/orange-ish = mountains, everything else unmarked),
    # base64-encoded so it rides along in the same JSON body as the rest of GenerateRequest --
    # produced either by the "Draw a map" in-app tool or a "Load an image" file picker.
    image_base64: str


class GenerateRequest(BaseModel):
    seed: int = 0
    # The Random tab's optional craton-first terrain/plate layout (RFC 0009). None keeps
    # the configured default; sketches and premade worlds ignore this mode.
    cratons_first: bool | None = None
    # With cratons_first, relief from a synthetic tectonic history (#218, synthetic_relief.py).
    # None keeps the configured default; ignored when the cratons-first layout doesn't apply.
    synthetic_relief: bool | None = None
    # Optional: the world tiles itself into a plausible plate count when omitted (see
    # plates.generate_plates) -- the frontend doesn't ask for one.
    num_plates: int | None = None
    # The UI's "Voronoi points" Advanced-settings slider -- the total number of Voronoi seed
    # points the plate tiling scatters before merging cells down to `num_plates` plates.
    # Optional; when omitted the tiling uses its fixed per-plate default (see
    # lithosphere_plate.EXTRA_SITES_PER_PLATE). More points -> lumpier plate outlines.
    voronoi_points: int | None = None
    # The UI's "continental plates" and "initial land" sliders (0 to 1) -- optional (each
    # falls back to its own default behavior, see plates.generate_plates) but the frontend
    # always sends both.
    continental_fraction: float | None = None
    land_fraction: float | None = None
    # The UI's "axial tilt" slider (degrees) -- optional, falls back to
    # world.DEFAULT_AXIAL_TILT_DEG, but the frontend always sends it.
    axial_tilt_deg: float | None = None
    num_mantle_centers: int = DEFAULT_MANTLE_CENTERS
    # The UI's "point density" choice -- how many surface nodes each plate starts
    # with, relative to elevation_lines.TARGET_LINE_SPACING_KM's own default spacing.
    # Validated against NODE_DENSITY_CHOICES below rather than accepted as an
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
    # The "Human-made" Generate World tab's drawn-or-loaded coastline (see worldsketch.py) --
    # when present, land/sea (and, where painted, mountains/rivers) come from this image
    # instead of the usual noise-driven generation, and plate sites/boundaries are fit to its
    # landmasses rather than scattered uniformly at random (see
    # lithosphere_plate.generate_plates's own `sketch` param). `land_fraction` above is
    # ignored when this is set -- the drawing decides land extent directly.
    sketch: SketchRequest | None = None
    # The "Premade worlds" tab's pick -- `"earth"`, `"pangaea"`, or `"got"` (see
    # premadeWorlds.ts's `backendId`) -- `sketch` above still supplies that world's coastline,
    # but its mountain belts/plateaus come from real-geography/lore regions instead of random
    # coverage, and for `"earth"`/`"pangaea"` specifically plate sites/mantle convection
    # additionally come from real plate geometry/motion (see world.generate_world's own
    # `premade_world_id` param). `None` (every other tab) is unaffected.
    premade_world_id: str | None = None
    # The Generate World dialog's name field and Public checkbox: applied to the
    # record being generated into (a fresh world keeps its automatic name when
    # this is omitted). Same bounds as PATCH /worlds/{world_id}'s renames.
    name: str | None = Field(default=None, min_length=1, max_length=100)
    is_public: bool | None = None


class StepRequest(BaseModel):
    years: float


class AnimateRequest(BaseModel):
    # Same shape/defaults as /worlds/{world_id}/render's own query params (see `render`), plus the two
    # animation-specific knobs below. `rotation` is a 9-comma-separated-float string, same
    # convention _parse_view_rotation already expects.
    projection: str = "behrmann"
    view: str = "elevation"
    width: int = 1100
    height: int = 611
    rotation: str | None = None
    # The UI's "Years per step" choice (same value /worlds/{world_id}/step's own `years` would use) --
    # sent explicitly rather than defaulted here so a future UI change doesn't need a matching
    # backend change.
    step_years: float
    # The UI's "Steps per frame" choice -- each frame after the first costs this many real
    # step_world calls of `step_years` each, only the last of which gets rendered (see
    # render_image.stream_animation_mp4's own docstring for why that's not the same thing as
    # one bigger step_world(world, steps_per_frame * step_years) call).
    steps_per_frame: int = 1
    # Total frames, including frame 0 (the world's current, unstepped state) -- so this
    # permanently advances the world by (num_frames - 1) * steps_per_frame * step_years years.
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
    # Full period (years) of the glacial<->interglacial cycle; 0 disables it. See
    # World.ice_age_period_years / climate.ice_age_cooling_c.
    ice_age_period_years: float | None = None
    simulate_plate_movement: bool | None = None
    simulate_climate_biomes: bool | None = None
    wind_model: str | None = None
    fault_deformation_mode: str | None = None
    # Issue #133 phase-1 proving-out flag -- see World.node_cloud_resample_mode. "kdtree"
    # (default) or "healpix". Backend/API-only for now, no Controls-panel entry.
    node_cloud_resample_mode: str | None = None
    # Geomorphic-budget tuning knobs -- dimensionless multipliers, 1.0 == untuned (see
    # world.TUNING_MULTIPLIER_FIELDS and World's field group). Same "only the touched one is
    # sent" convention as the fields above; each is rejected below if negative.
    rain_erosion_multiplier: float | None = None
    river_erosion_multiplier: float | None = None
    wind_erosion_multiplier: float | None = None
    ocean_erosion_multiplier: float | None = None
    coastal_leveling_multiplier: float | None = None
    glacier_erosion_multiplier: float | None = None
    seismic_erosion_multiplier: float | None = None
    river_deposition_multiplier: float | None = None
    ocean_deposition_multiplier: float | None = None
    collision_uplift_multiplier: float | None = None
    collision_uplift_reach_multiplier: float | None = None
    volcanism_multiplier: float | None = None
    fault_relief_multiplier: float | None = None
    plume_multiplier: float | None = None
    # Gate for the debug-only diagnostics (see World.debug_diagnostics) -- on by default for a
    # "Debugging Worlds" tab world, off/toggleable here for any other loaded save.
    debug_diagnostics: bool | None = None


WIND_MODEL_CHOICES = ("cfd", "diagnostic")
FAULT_DEFORMATION_MODE_CHOICES = faults.FAULT_DEFORMATION_MODES
NODE_CLOUD_RESAMPLE_MODE_CHOICES = healpix_grid.NODE_CLOUD_RESAMPLE_MODE_CHOICES


def _parse_view_rotation(rotation: str | None) -> np.ndarray:
    """The map's current view orientation (see docs/simulation-model/rendering-and-inspectors.md#rotating-the-view),
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


def _session_or_404(world_id: str) -> WorldSession:
    session = _registry.find(world_id)
    if session is None:
        raise HTTPException(status_code=404, detail=f"unknown world {world_id!r}")
    return session


def _session_running(session: WorldSession) -> bool:
    """Whether this world's run loop is on: what /view calls live (#597)."""
    runner = session.runner
    return runner is not None and runner.status().running


def _apply_generate_identity(session: WorldSession, name: str | None, is_public: bool | None) -> None:
    """Name (and publish) the record a generate/debug-generate is about to fill (#597)."""
    if name is None and is_public is None:
        return
    with session.lock:
        changed = False
        if name is not None and name != session.name:
            session.name = name
            changed = True
        if is_public is not None and is_public != session.is_public:
            session.is_public = is_public
            changed = True
    if changed:
        _save_index()


def _public_session_or_404(world_id: str) -> WorldSession:
    """A session the public viewer may address. Private idle worlds and unknown
    IDs are indistinguishable (both 404) so the public API cannot reveal private
    worlds -- but a world whose run loop is on answers while it runs (#597), so
    anyone can watch a planet that is currently stepping."""
    session = _registry.find(world_id)
    if session is None or session.world is None:
        raise HTTPException(status_code=404, detail="no such public world")
    if not session.is_public and not _session_running(session):
        raise HTTPException(status_code=404, detail="no such public world")
    return session


def _world_changed(session: WorldSession, completed: bool = True) -> None:
    """Publish a new revision of `session`'s world to its viewers. Call it after every change
    to the world, before releasing the session lock: a viewer's render takes the lock and
    caches what it draws under the revision it reads there, so a change published after the
    release could leave the new world cached under the old revision.

    A completed change is also checkpointed (#375), under the same lock so nothing mutates
    the world mid-pickle. `completed=False` is for a step that raised partway: viewers still
    see what it left, but a restart goes back to the last whole step rather than resuming
    from a half-stepped world."""
    session.touch()
    # Persist the touched timestamp, not just the world: the frontend auto-selects
    # the most recently updated world, and without this a restart would fall back
    # to the stale updated_at from the record's creation. A tiny atomic write.
    _save_index()
    session.public_renders.clear()
    session.feed.publish(session.world)
    if completed:
        _checkpoint_world(session)


def _set_world(session: WorldSession, world: World) -> None:
    """Replace `session`'s world; the caller holds the session lock."""
    session.world = world
    session.load_error = None
    if session.replay_log is not None:
        session.replay_log.clear_recent_events()
    _world_changed(session)


def _require_world(session: WorldSession) -> World:
    world = session.world
    if world is None:
        raise HTTPException(status_code=404, detail="no world generated yet")
    return world


def _record_json(session: WorldSession) -> dict:
    """One registry record for the admin API: metadata plus, when populated,
    the world summary and run state."""
    record = session.record()
    record["load_error"] = session.load_error
    runner = session.runner
    record["run"] = None if runner is None else {
        "running": runner.status().running,
        "step_years": runner.status().step_years,
        "stepping": runner.stepping(),
        "error": runner.error(),
    }
    record["summary"] = None if session.world is None else _summary(session, session.world)
    return record


def _require_valid_step_years(years: float, field: str = "years") -> None:
    """400 unless `years` is a positive, finite step of at most MAX_STEP_YEARS (#557)."""
    if not checkpoint.valid_step_years(years) or years > MAX_STEP_YEARS:
        raise HTTPException(
            status_code=400,
            detail=f"{field} must be positive, finite, and at most {MAX_STEP_YEARS:,.0f}",
        )


def _summary(session: WorldSession, world: World) -> dict:
    return {
        "world_id": session.world_id,
        "revision": session.feed.revision(),
        "seed": world.seed,
        "elapsed_years": world.elapsed_years,
        "num_plates": len(world.plates),
        "surface": persistence.world_surface(world),
        "surface_conversion": world.surface_conversion,
        # Full current log (bounded, see world.MAX_EVENT_LOG_LENGTH) every time -- simpler
        # for the frontend than diffing "new since last call", and small enough not to
        # matter on the wire.
        "events": [{"elapsed_years": e, "message": m} for e, m in world.events],
        # Replay-log records this process wrote (#442), for the ?deb Event Console. Kept out of
        # `events`, which is pickled World state.
        "replay_logs": session.replay_log.recent_events() if session.replay_log is not None else [],
    }


#  Rounding every coordinate sent to the client to this many decimal places -- plenty for a
# unit sphere (1e-6 rad is sub-meter at planet scale, far finer than anything visible on
# screen) -- meaningfully shrinks the /worlds/{world_id}/plates payload now that it includes every node
# point per plate (`_plate_summary`'s "points" field), not just the much shorter outline/
# ellipse loops.
_COORD_DECIMALS = 6


def _round_coords(arr: np.ndarray) -> list:
    return np.round(arr, _COORD_DECIMALS).tolist()


def _plate_overlaps(world: World) -> dict[int, list[dict]]:
    """For every plate, which other plates its territory currently overlaps and by how much
    (`fraction` = share of *this* plate's own nodes that sit on top of that other plate),
    plus `since_years` -- the earliest `world.elapsed_years` at which any of this plate's
    still-overlapping nodes first went over that partner (the `overlap_onset_years` field,
    stamped by merge_split.update_overlap_tracking; None if the save predates the field or the
    overlap only appeared this step). The node-cloud overlap itself is
    `plates.compute_node_overlap`, shared with that tracker so the two can't drift."""
    overlap = plates.compute_node_overlap(world.plates)
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
        # long-run pathology when it's true for most plates -- see GitHub issue #119).
        "speed_cm_per_yr": round(mantle.rad_per_yr_to_cm_per_yr(rate_rad), 3),
        "at_max_rate": bool(rate_rad >= mantle.MAX_PLATE_RATE - 1e-12),
        "euler_pole": euler_pole,
        "age_years": plate.age_years,
        # Median height of the plate's own nodes and the share of them at/below sea level --
        # a continental plate reading mostly-submerged is a red flag for over-stretching.
        "median_elevation_m": round(float(np.median(node_elevation)), 1) if len(node_elevation) else None,
        "submerged_fraction": round(submerged, 3),
        # Other plates this plate's territory currently sits on top of (see _plate_overlaps),
        # and any sustained-collision timers involving it (merge_split.update_collision_progress).
        "overlaps": overlaps.get(plate.plate_id, []),
        "collisions": collisions,
        "num_points": plate.node_count(),
        "outline": _round_coords(outline),
        # Every node's own position (not just the outline loop) -- lets the client plot each
        # point individually (see docs/simulation-model/rendering-and-inspectors.md#plate-inspector), highlighted for
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
    from_idx, to_idx, _, _ = hydrology.flow_segments(fields, river.member_idx)
    from_points = fields.points[from_idx]
    to_points = fields.points[to_idx]
    segments = np.round(np.stack([from_points, to_points], axis=1), _COORD_DECIMALS).tolist()
    return {
        "river_id": river_id,
        "num_nodes": int(len(river.member_idx)),
        "length_m": river.length_m,
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


def _river_by_node(rivers: list[hydrology.RiverInfo]) -> dict[int, tuple[int, hydrology.RiverInfo]]:
    """Every river-network member node mapped to `(river_id, river)` -- built once per request
    (same "compute the list, then look members up in it" shape as `_leaf_lakes_by_node`) so
    `_lake_basin_summary`'s outflow-river lookup (see its own docstring) doesn't re-scan every
    river's `member_idx` for every basin in a `/worlds/{world_id}/lakes` call."""
    by_node: dict[int, tuple[int, hydrology.RiverInfo]] = {}
    for river_id, river in enumerate(rivers):
        for member in river.member_idx.tolist():
            by_node[member] = (river_id, river)
    return by_node


def _lake_components_sorted(fields: hydrology.HydrologyFields) -> list[np.ndarray]:
    """Every currently-visible lake's node group (`hydrology.lake_components`), in a fixed,
    deterministic order (by each component's own lowest node index) -- so the index into this
    list (`lake_id`) means the same thing across a `/worlds/{world_id}/lakes` call and a same-step
    `/worlds/{world_id}/lake_at` call, the same "meaningful only against this step" contract `river_id`
    already has (see hydrology.group_rivers)."""
    components = hydrology.lake_components(fields.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M, fields.neighbor_idx)
    components.sort(key=lambda members: int(members.min()))
    return components


def _lake_basin_summary(
    fields: hydrology.HydrologyFields,
    lake: lakes.Lake,
    lake_id: int | None,
    rivers: list[hydrology.RiverInfo],
    river_by_node: dict[int, tuple[int, hydrology.RiverInfo]],
    is_lake: bool,
) -> dict:
    """One basin's full inspector payload, for both a currently-wet lake (`is_lake=True`, `lake`
    is whatever `_smallest_lake_containing` resolved) and a dry, never-flooded basin
    (`is_lake=False`, `lake` is the leaf catchment straight from `_leaf_lakes_by_node`) -- same
    shape either way, since "a basin nonetheless" (the user's own framing) should read as the
    same kind of information, just without live water. `outlet_*`/`water_elevation_m`/
    `is_spilling` are read directly off `lake` itself (`max_depth`/
    `current_water_elevation`/`is_spilling`, already resolved by `lakes.step_lakes` -- see
    `HydrologyFields.lake_forest`'s own docstring for why this must be the *same* `Lake` object
    that produced `fields.lake_depth`, not a hierarchy rebuilt from a different elevation/silt
    snapshot: the two can genuinely disagree on basin topology, which used to surface as a
    lake's own `outlet_elevation_m` reading *below* its `floor_elevation_m`) -- `outlet_*` is
    `None` for an unresolved closed/endorheic basin with no known spill (a legitimate case, see
    lakes.py's own docstring), not a bug. Rivers "feeding into" this basin are simply every
    `RiverInfo` whose own mouth lands on one of this basin's members, regardless of that
    river's own `mouth_type` label ("lake", or "other" for a river dead-ending in a currently-dry
    sink) -- from this basin's own point of view both are equally "a river that ends here".

    The river actually carrying this basin's *outflow* is different in kind, not just
    direction, from an inflow: an inflow is identified by its mouth landing somewhere inside
    `members` (this basin is wherever it ends), but a basin's own outflow starts just *outside*
    `members` -- at `lake.outlet_target_idx`, the node on the far side of the governing rim
    (see lakes.py's own `compute_spill_routing` docstring) -- so it's found by which river
    network (if any) that specific node belongs to, via `river_by_node`. `None` whenever
    there's no live outflow to show at all: a closed/endorheic basin with no known spill
    (`outlet_target_idx == -1`), a resolved rim that isn't currently `is_spilling` (nothing is
    actually flowing out this step), or -- a real, if narrow, gap right at the moment a lake
    first breaches -- a spilling rim whose own flow hasn't yet cleared
    `hydrology.RIVER_FLOW_PERCENTILE` to register as its own river network at all."""
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
    # position. The reported value is the basin's lowest ground under the same hypsometry
    # `lakes.resolve_lakes` used: `lake.floor_elevation` is the lowest cell *mean*, and a lake
    # wet only in its floor cell's low ground has its surface below that (issue #344). That
    # ground is never above `lake.floor_elevation`, so it stays at or below `outlet_elev`.
    floor_node = int(display_members[np.argmin(fields.elevation[display_members])])
    floor_elevation = lake.floor_elevation if fields.hypsometry is None else fields.hypsometry.dry_level(members)
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
            "length_m": river.length_m,
        }
        for river in rivers
        if river.mouth_idx in member_set
    ]

    # See this function's own docstring for why an outflow is looked up by river membership at
    # `outlet_target_idx`, not by mouth like an inflow -- and why `None` covers three distinct,
    # all-legitimate cases (no known spill at all, a resolved-but-dry rim, and a spilling rim
    # too fresh to have grown its own river network yet).
    outflow_river = None
    if is_spilling and lake.outlet_target_idx >= 0:
        found = river_by_node.get(lake.outlet_target_idx)
        if found is not None:
            river_id, river = found
            outflow_river = _river_summary(fields, river, river_id)

    return {
        "lake_id": lake_id,
        "is_lake": is_lake,
        "member_count": int(len(display_members)),
        "member_xyz": _round_coords(fields.points[display_members]),
        "floor_xyz": _round_coords(fields.points[floor_node]),
        "floor_elevation_m": float(floor_elevation),
        "outlet_xyz": None if outlet_elev is None else _round_coords(fields.points[outlet_node]),
        "outlet_elevation_m": outlet_elev,
        "water_elevation_m": None if water_elevation is None else float(water_elevation),
        "is_spilling": bool(is_spilling),
        "inflow_rivers": inflow_rivers,
        "outflow_river": outflow_river,
    }


def _drain_progress(progress) -> World | None:
    """Finishes driving a generate_world_progress/step_world_progress generator to
    completion with nobody left listening -- shared by generate()'s and step()'s own
    `except GeneratorExit` handling below (see either's comment on it for why a disconnected
    client can't just abandon one of these generators wherever it happened to be suspended).
    Returns the World from generate_world_progress's final `("done", world)` message, or
    `None` for step_world_progress (which mutates its own `world` argument in place and
    yields plain floats, never a `("done", ...)` tuple -- there's nothing further for a
    caller to do with its return value).

    If generation/stepping itself raises while finishing up here, there's no request left to
    report it to (the client that would have received a `{"type": "error"}` line is already
    gone) -- letting that exception propagate out of this function would instead replace the
    caller's pending `GeneratorExit` mid-`.close()`, which Python can only report as an
    unreachable "exception ignored in generator" at GC time. Logging it here instead is the
    difference between a genuine failure (e.g. a degenerate-geometry edge case) being visible
    at all versus vanishing silently, at the cost of the session world being
    left however far the drain got before the error -- already true of a mid-request failure
    on a normal (non-disconnected) call, so no worse a guarantee than that."""
    world = None
    try:
        for message in progress:
            if isinstance(message, tuple) and message[0] == "done":
                world = message[1]
    except Exception:
        _logger.exception("world mutation failed while finishing up after a client disconnect")
    return world


class CreateWorldRequest(BaseModel):
    # Bounded: an empty name renders as a blank option (and restores as
    # "Untitled world", so live and persisted names would differ), while a
    # huge one is rewritten to worlds.json on every record save.
    name: str | None = Field(default=None, min_length=1, max_length=100)
    is_public: bool = False


class PatchWorldRequest(BaseModel):
    name: str | None = Field(default=None, min_length=1, max_length=100)
    is_public: bool | None = None


@app.get("/worlds")
def list_worlds() -> dict:
    """Every world record in this process (#563)."""
    worlds = []
    for session in _registry.list():
        # The summary reads world arrays a concurrent step mutates in place;
        # hold the session lock so a record can't tear mid-step.
        with session.lock:
            worlds.append(_record_json(session))
    return {"worlds": worlds}


@app.post("/worlds", status_code=201)
def create_world(req: CreateWorldRequest) -> dict:
    """Create an empty world record. Generate or Load populates it afterwards."""
    try:
        session = _registry.create(name=req.name, is_public=req.is_public)
    except WorldLimitReached as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    _wire_session(session)
    _save_index()
    return _record_json(session)


@app.get("/worlds/{world_id}")
def get_world(world_id: str) -> dict:
    session = _session_or_404(world_id)
    # Like list_worlds and get_summary: the record carries a world summary over
    # mutable arrays, so hold the lock against a mid-step torn read.
    with session.lock:
        return _record_json(session)


@app.patch("/worlds/{world_id}")
def patch_world(world_id: str, req: PatchWorldRequest) -> dict:
    session = _session_or_404(world_id)
    with session.lock:
        if req.name is not None:
            session.name = req.name
        if req.is_public is not None:
            session.is_public = req.is_public
        session.touch()
        _save_index()
        return _record_json(session)


@app.delete("/worlds/{world_id}")
def delete_world(world_id: str) -> dict:
    session = _session_or_404(world_id)
    if session.operation_busy():
        raise HTTPException(status_code=409, detail="a generation, step, or animation is in progress on this world")
    if session.guest_steps.busy():
        raise HTTPException(status_code=409, detail="a guest step is in progress on this world; try again in a moment")
    if session.runner is not None and session.runner.stepping():
        raise HTTPException(status_code=409, detail="a run step is under way on this world; try again in a moment")
    # Mark first: streams that queued at the compute gate before this delete re-check
    # `deleted` under the session lock when their turn arrives and bail out instead
    # of working on a detached session (see each _stream's check). A runner waking
    # between here and close() sees the same flag in _run_step_once and stands down.
    session.deleted = True
    if session.runner is not None:
        session.runner.set_status(checkpoint.RunStatus(running=False, step_years=session.runner.status().step_years), persist=False)
        session.runner.close()
    _registry.remove(world_id)
    if _checkpoint_root is not None and session.checkpoint_store is not None:
        import shutil

        # Persist the removal before deleting the files: a crash between the two
        # must leave an orphaned directory (restore ignores it -- no record
        # references it), never an index entry pointing at files that no longer
        # exist. A checkpoint can't be mid-save here: writers either hold the
        # operation flag (refused above), are the runner (closed above), or are a
        # guest step (refused above).
        _save_index()
        shutil.rmtree(session.checkpoint_store.directory, ignore_errors=True)
    return {"deleted": world_id}


@app.post("/worlds/{world_id}/generate")
def generate(world_id: str, req: GenerateRequest) -> StreamingResponse:
    """Streams newline-delimited JSON progress while generating a fresh world, replacing
    whatever world previously existed -- same NDJSON progress/done/error contract as
    /worlds/{world_id}/animate (see that endpoint's own docstring): a `{"type": "progress", "fraction":
    f}` line (`f` from 0 to 1) after each of world.generate_world_progress's three phase
    boundaries (plate/site generation, mantle-center fitting, the finish_generation
    bootstrap -- coarser-grained than /worlds/{world_id}/step's per-plate progress, see that generator's
    own docstring for why), then a final `{"type": "done", ...summary fields}` once the world
    is fully built, or a `{"type": "error", "detail": ...}` line if generation raises
    partway through (the HTTP status is already 200 by then, same as animate). `404` for an
    unknown world id (there's no existing world *content* required); `503` if a
    generation/step/animation is already in progress on this world (see
    `session.begin_operation` below). Every other validation
    (density/voronoi_points/premade_world_id/sketch) still happens synchronously before the
    stream starts, so a bad request still gets an ordinary 400 response rather than an error
    line."""
    session = _session_or_404(world_id)
    _apply_generate_identity(session, req.name, req.is_public)
    if req.node_density not in NODE_DENSITY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown node_density {req.node_density!r}; choices are {NODE_DENSITY_CHOICES}")
    if req.climate_density not in climate.CLIMATE_DENSITY_CHOICES:
        raise HTTPException(
            status_code=400, detail=f"unknown climate_density {req.climate_density!r}; choices are {climate.CLIMATE_DENSITY_CHOICES}"
        )
    if req.fluid_density not in climate.FLUID_DENSITY_CHOICES:
        raise HTTPException(
            status_code=400, detail=f"unknown fluid_density {req.fluid_density!r}; choices are {climate.FLUID_DENSITY_CHOICES}"
        )
    if req.voronoi_points is not None and not 1 <= req.voronoi_points <= 2000:
        raise HTTPException(status_code=400, detail="voronoi_points must be between 1 and 2000")
    if req.premade_world_id is not None and req.premade_world_id not in ("earth", "pangaea", "got"):
        raise HTTPException(status_code=400, detail=f"unknown premade_world_id {req.premade_world_id!r}")
    if req.premade_world_id is not None and req.sketch is None:
        raise HTTPException(status_code=400, detail="premade_world_id requires sketch")
    sketch_masks = None
    if req.sketch is not None:
        try:
            sketch_masks = worldsketch.parse_sketch_image(base64.b64decode(req.sketch.image_base64))
        except Exception as exc:
            raise HTTPException(status_code=400, detail=f"invalid sketch image: {exc}") from exc

    # Mark synchronously (like animate's own, see its comment) -- an overlapping
    # generate/step/animation on this world still 503s before the stream starts.
    # The compute gate itself is taken inside the stream, so work on another
    # world queues FIFO instead of failing.
    if not session.begin_operation():
        raise HTTPException(status_code=503, detail="a generation, step, or animation is already in progress")

    def _stream():
        # Taken here, not synchronously above: another world's compute queues
        # at the gate FIFO instead of 503ing. Guard the release: a generator
        # closed before its first step holds nothing yet.
        acquired = False
        _compute_gate.acquire()
        try:
            session.lock.acquire()
        except BaseException:
            # Lock.acquire() raises only on interpreter-level events, but if it
            # did, `acquired` would stay False while the gate is held -- a wedge
            # every world shares and no flag tracks. Release it first.
            _compute_gate.release()
            raise
        acquired = True
        try:
            # The world may have been deleted while this stream queued at the gate.
            # Bail out quietly rather than generating into a detached session. This
            # check sits inside the try on purpose: returning before the try would
            # skip the finally below and leak the compute gate (wedging every
            # world), the session lock and the operation flag.
            if session.deleted:
                return
            progress = generate_world_progress(
                req.seed,
                cratons_first=req.cratons_first,
                synthetic_relief=req.synthetic_relief,
                num_plates=req.num_plates,
                voronoi_points=req.voronoi_points,
                continental_fraction=req.continental_fraction,
                land_fraction=req.land_fraction,
                num_mantle_centers=req.num_mantle_centers,
                axial_tilt_deg=req.axial_tilt_deg,
                node_density=req.node_density,
                initial_soil_maturity=req.initial_soil_maturity,
                climate_density=req.climate_density,
                fluid_density=req.fluid_density,
                sketch=sketch_masks,
                premade_world_id=req.premade_world_id,
            )
            world = None
            try:
                for message in progress:
                    if message[0] == "progress":
                        yield json.dumps({"type": "progress", "fraction": message[1]}) + "\n"
                    else:
                        world = message[1]
                _set_world(session, world)
                yield json.dumps({"type": "done", **_summary(session, world)}) + "\n"
            except GeneratorExit:
                # The client disconnected mid-stream -- Starlette closes a sync generator behind
                # StreamingResponse by throwing this in at whichever `yield` above is suspended
                # (see iterate_in_threadpool), same as it would for a `for` loop that `break`s
                # early. Generation was never a resumable/abortable transaction (unlike
                # /worlds/{world_id}/animate's per-frame stopEvent, which is deliberate and cooperative, and
                # unlike animate's own `_stream` -- see its comment -- which doesn't need this at
                # all) -- world.py's generate_world_progress mutates its own local `world` across
                # multiple phases with no valid "half-generated" state to leave sitting in
                # `_state`, so _drain_progress finishes driving it to completion here regardless
                # of whether anyone is still listening, the same as this endpoint's single
                # synchronous response body used to unconditionally do before it streamed
                # progress. Only `_state` needs updating at this point -- there's no client left
                # to send progress/done lines to.
                world = _drain_progress(progress) or world
                if world is not None:
                    _set_world(session, world)
                raise
            except Exception as exc:  # noqa: BLE001 -- status is already 200, so surface it as a data line
                yield json.dumps({"type": "error", "detail": f"{type(exc).__name__}: {exc}"}) + "\n"
        finally:
            if acquired:
                session.lock.release()
                _compute_gate.release()
            session.end_operation()

    return StreamingResponse(_stream(), media_type="application/x-ndjson")


class GenerateDebugRequest(BaseModel):
    scenario: str
    seed: int = 0
    # Same dialog name/public choice as GenerateRequest above.
    name: str | None = Field(default=None, min_length=1, max_length=100)
    is_public: bool | None = None


@app.post("/worlds/{world_id}/generate_debug")
def generate_debug(world_id: str, req: GenerateDebugRequest) -> dict:
    """The "Debugging Worlds" Generate World tab -- tiny, hand-scripted plate configurations
    with pinned motion (see debug_worlds.py), for fast iteration on the gap-filling problem
    independent of any real save. Unlike `/worlds/{world_id}/generate`, there's no Voronoi tiling, plate
    count, or land/continental fraction to choose -- each scenario name is its own complete,
    fixed configuration; only `seed` (for the per-node terrain texture, see
    `lithosphere_plate.new_plate`) varies. `400` for an unknown scenario name, `404`
    for an unknown world id, `503` while a generation, step, or animation is in
    progress."""
    session = _session_or_404(world_id)
    _apply_generate_identity(session, req.name, req.is_public)
    if req.scenario not in debug_worlds.DEBUG_SCENARIOS:
        raise HTTPException(
            status_code=400,
            detail=f"unknown scenario {req.scenario!r}; choices are {sorted(debug_worlds.DEBUG_SCENARIOS)}",
        )
    # Joins the operation-flag protocol like every other world replacer: step
    # captures its world object at request time and then queues at the compute
    # gate, so a debug generation slipping in between would swap session.world
    # out from under it -- the step would mutate a detached world whose work
    # never reaches the published state.
    if not session.begin_operation():
        raise HTTPException(status_code=503, detail="a generation, step, or animation is already in progress")
    try:
        with _compute_gate.hold(), session.lock:
            if session.deleted:
                raise HTTPException(status_code=404, detail="no such world")
            world = debug_worlds.generate_debug_world(req.scenario, req.seed)
            _set_world(session, world)
            # Summarized under the lock: the arrays a concurrent step mutates in
            # place must not be read mid-step.
            return _summary(session, world)
    finally:
        session.end_operation()


@app.get("/debug/scenarios")
def debug_scenarios() -> dict:
    """The "Debugging Worlds" tab's scenario picker options -- name plus a short human label
    (debug_worlds.DEBUG_SCENARIO_LABELS), so the frontend doesn't hardcode the list."""
    return {
        "scenarios": [
            {"name": name, "label": debug_worlds.DEBUG_SCENARIO_LABELS[name]}
            for name in debug_worlds.DEBUG_SCENARIOS
        ]
    }


@app.get("/worlds/{world_id}/summary")
def get_summary(world_id: str) -> dict:
    """The addressed world's summary (same shape its generate and step return) --
    lets the client re-sync with whatever world is already sitting in server memory (e.g.
    after a browser refresh, see App.tsx's restore-on-mount effect and its view-state cookie)
    without generating a new one. `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    # Under the lock: a concurrent step mutates the world's numpy arrays in place,
    # and reading them mid-step would report a torn, half-stepped world.
    with session.lock:
        world = _require_world(session)
        return _summary(session, world)


@app.post("/worlds/{world_id}/step")
def step(world_id: str, req: StepRequest) -> StreamingResponse:
    """Advances plate tectonics/climate/erosion by `req.years` (see world.step_world_progress)
    -- and, gated on World.simulate_climate_biomes the same way erosion/hydrology already
    are, Ocean/Atmospheric Fluid Dynamics by their own fixed real-time increment regardless
    of `req.years` (see world.py's `_advance_fluid_dynamics`). Streams newline-delimited
    JSON progress the same shape /worlds/{world_id}/generate and /worlds/{world_id}/animate do: a `{"type":
    "progress", "fraction": f}` line after each per-plate shift()/deform() call (see
    step_world_progress's own docstring for why that's a real, if coarse, measure of
    completion), then a final `{"type": "done", ...summary fields}`, or a `{"type": "error",
    ...}` line if the step raises partway through. `404` if no world has been generated yet,
    `503` (before the stream starts, same as /worlds/{world_id}/animate) if a generation, step, or
    animation is already in progress."""
    session = _session_or_404(world_id)
    _require_valid_step_years(req.years)
    if not session.begin_operation():
        raise HTTPException(status_code=503, detail="a generation, step, or animation is already in progress")
    # Captured under the flag: a replacer swapping session.world between this read
    # and the stream's turn at the gate would leave the step mutating a detached
    # world whose work never reaches the published state.
    world = _require_world(session)

    def _stream():
        acquired = False
        _compute_gate.acquire()
        try:
            session.lock.acquire()
        except BaseException:
            # Lock.acquire() raises only on interpreter-level events, but if it
            # did, `acquired` would stay False while the gate is held -- a wedge
            # every world shares and no flag tracks. Release it first.
            _compute_gate.release()
            raise
        acquired = True
        try:
            # The world may have been deleted while this stream queued at the gate.
            # Bail out quietly rather than stepping a detached session. Inside the
            # try on purpose: returning before it would skip the finally and leak
            # the gate, the session lock and the operation flag (see generate).
            # Initialized before the bail: the finally below reads it.
            completed = False
            if session.deleted:
                return
            progress = _step_progress(session, world, req.years)
            try:
                for fraction in progress:
                    yield json.dumps({"type": "progress", "fraction": fraction}) + "\n"
                completed = True
                yield json.dumps({"type": "done", **_summary(session, world)}) + "\n"
            except GeneratorExit:
                # The client disconnected mid-stream (see /worlds/{world_id}/generate's own comment on this
                # exact exception for why Starlette raises it here). A step mutates `world` in
                # place across many per-plate checkpoints with no valid "half-stepped" state to
                # leave it in -- world.elapsed_years not yet advanced but some plates already
                # shifted/deformed, no topology/climate/erosion/sea-level pass -- so finish
                # driving step_world_progress to completion here regardless, the same as this
                # endpoint's single synchronous response body used to unconditionally do before
                # it streamed progress. `world` is mutated in place (unlike generate_world_progress,
                # which builds a fresh World to hand back), so there's nothing further to write to
                # `_state` here -- just let the step finish.
                _drain_progress(progress)
                completed = True
                raise
            except Exception as exc:  # noqa: BLE001 -- status is already 200, so surface it as a data line
                yield json.dumps({"type": "error", "detail": f"{type(exc).__name__}: {exc}"}) + "\n"
        finally:
            if acquired:
                # Even a failed step has changed the world partway.
                _world_changed(session, completed)
                session.lock.release()
                _compute_gate.release()
            session.end_operation()

    return StreamingResponse(_stream(), media_type="application/x-ndjson")


@app.get("/worlds/{world_id}/save")
def save_world(world_id: str) -> Response:
    """The "File > Save World" download -- the *entire* current world (every plate/line,
    mantle field, caches, event log -- see World) pickled as a single opaque file (see
    persistence.py for why pickle, and docs/save-compatibility.md for which older saves load).
    `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    # Serialized under the lock: a concurrent step mutates the world's numpy arrays in
    # place, and pickling them mid-step would write a torn, unloadable file.
    with session.lock:
        world = _require_world(session)
        body = persistence.save_world_bytes(world)
        filename = f"mantle-bloom-seed{world.seed}-{int(world.elapsed_years)}y.mbworld"
    return Response(
        content=body, media_type="application/octet-stream", headers={"Content-Disposition": f'attachment; filename="{filename}"'}
    )


@app.post("/worlds/{world_id}/load")
async def load_world(world_id: str, request: Request) -> dict:
    """The "File > Load World" upload -- the raw bytes of a file /worlds/{world_id}/save previously
    produced, as the request body (not JSON -- see persistence.py). Replaces whatever world
    previously existed, same as /worlds/{world_id}/generate. A line-backed (legacy) save is converted to
    sparse quads on load, one way -- see docs/save-compatibility.md; the summary's
    `surface_conversion` then reports what changed. `400` with the reason if the
    bytes are corrupt or in a save format this build doesn't read; any other failure is
    still a `400`, without detail. `503` while a generation, step, or animation is in
    progress."""
    body = await request.body()
    try:
        world = persistence.load_world_bytes(body)
    except persistence.SaveFormatError as exc:
        raise HTTPException(status_code=400, detail=f"invalid or incompatible world file: {exc}") from exc
    except Exception as exc:
        raise HTTPException(status_code=400, detail="invalid or incompatible world file") from exc
    # Physical discharge is derived routing state. Saves from before #355 may carry a cached
    # HydrologyFields without it; discard that cache so River/Lake Inspector APIs don't label
    # legacy stacked-mm/yr values as m^3/s. The next simulated step rebuilds the cache.
    if world.hydrology_cache is not None and getattr(world.hydrology_cache, "discharge_m3_s", None) is None:
        world.hydrology_cache = None
        world.hydrology_cache_step = None
    session = _session_or_404(world_id)
    # Under the lock, so a step can't still be mutating the world viewers are told replaced it.
    if not session.begin_operation():
        raise HTTPException(status_code=503, detail="a generation, step, or animation is already in progress")
    try:
        # Non-blocking: this route is `async`, so parking on the session lock here
        # would freeze the event loop -- every long-poll and async route on every
        # world -- for as long as the step or animation holding it runs. The runner
        # and guest steps hold the lock without setting the busy flag, so refuse
        # instead, like an overlapping manual operation would.
        if not session.lock.acquire(blocking=False):
            raise HTTPException(status_code=503, detail="a generation, step, or animation is already in progress")
        try:
            if session.deleted:
                raise HTTPException(status_code=404, detail="no such world")
            _set_world(session, world)
            return _summary(session, world)
        finally:
            session.lock.release()
    finally:
        session.end_operation()


@app.get("/worlds/{world_id}/render")
def render(
    world_id: str,
    projection: str = "behrmann",
    view: str = "elevation",
    width: int = 1100,
    height: int = 611,
    rotation: str | None = None,
    show_mountains: bool = False,
    show_plains_plateaus: bool = False,
) -> dict:
    """Renders `view` of the current world, in `projection`, as a `width`x`height` PNG,
    returned base64-encoded in `image_base64` -- see render_image.py for the actual
    rasterization. `rotation` is the map's current view orientation (see
    _parse_view_rotation), default identity. `400` for an unrecognized projection/view, an
    out-of-range width/height, or a malformed rotation, `404` if no world has been generated
    yet.

    `show_mountains`/`show_plains_plateaus` are the Elevation view's own legend toggles (see
    render_image._classify_terrain_relief) -- both default `False`, and are silently ignored
    (no validation error) on every `view` but `"elevation"`, same as any other view-specific
    param in this codebase.

    Takes the compute gate and session lock (blocking) around the actual read, unlike
    step/animate which 503 on a busy world and let the caller retry. A render is a read, called far more often (every view
    switch, every post-step refresh, each animation frame) and each mutation-collecting pass
    it makes (plates.collect_all_points, collect_all_lake_depth, ...) walks world.plates
    independently -- concurrently with an in-flight step actually mutating those same plates
    (deform() grows/shrinks nodes), those passes can see different node counts and desync,
    which crashes with an IndexError rather than just rendering a stale frame. Blocking here
    (confirmed to actually reproduce that IndexError without it, via a deliberately raced
    step+render) waits out that brief window instead of failing every render that happens to
    land during one -- a hard failure on nearly every routine step->refresh overlap would be
    a far worse trade for a request this frequent and this cheap to just wait a moment for."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    if projection not in projections.PROJECTIONS:
        raise HTTPException(status_code=400, detail=f"unknown projection {projection!r}")
    if view not in render_image.VIEWS + render_image.LAYER_VIEWS:
        raise HTTPException(
            status_code=400, detail=f"unknown view {view!r}; choices are {render_image.VIEWS + render_image.LAYER_VIEWS}"
        )
    if not (1 <= width <= MAX_RENDER_DIMENSION_PX and 1 <= height <= MAX_RENDER_DIMENSION_PX):
        raise HTTPException(status_code=400, detail=f"width/height must be in [1, {MAX_RENDER_DIMENSION_PX}]")
    view_rotation = _parse_view_rotation(rotation)

    with _compute_gate.hold(), session.lock:
        image_base64 = render_image.render_png_base64(
            world, projection, view, width, height, view_rotation, show_mountains, show_plains_plateaus
        )
    return {
        "projection": projection,
        "elapsed_years": world.elapsed_years,
        "image_base64": image_base64,
    }


@app.get("/worlds/{world_id}/cross-section")
def cross_section(
    world_id: str,
    plane_angle_deg: float = 0.0,
    position_index: int = 0,
    width: int = 1100,
    height: int = 611,
    show_boundaries: bool = False,
    show_faults: bool = False,
) -> dict:
    """Render one surface-sampled section for the shallow mantle view (#360).

    The image uses existing per-column crust/lithosphere thicknesses and an illustrative
    geotherm; the mantle below modeled lithosphere is schematic. `404` if no world exists.
    """
    session = _session_or_404(world_id)
    with _compute_gate.hold(), session.lock:
        world = _require_world(session)
        if not -90.0 <= plane_angle_deg <= 90.0:
            raise HTTPException(status_code=400, detail="plane_angle_deg must be in [-90, 90]")
        if not 0 <= position_index < 120:
            raise HTTPException(status_code=400, detail="position_index must be in [0, 119]")
        if not (1 <= width <= MAX_RENDER_DIMENSION_PX and 1 <= height <= MAX_RENDER_DIMENSION_PX):
            raise HTTPException(status_code=400, detail=f"width/height must be in [1, {MAX_RENDER_DIMENSION_PX}]")
        try:
            _, result = cross_section_view.render_png_base64(
                world, plane_angle_deg, position_index, width, height,
                show_boundaries=show_boundaries, show_faults=show_faults,
            )
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
    return result


@app.post("/worlds/{world_id}/animate")
def animate(world_id: str, req: AnimateRequest) -> StreamingResponse:
    """The "File > Make Animation" action: streams newline-delimited JSON progress while
    rendering an H.264/MP4 video of `view`/`projection`'s progress -- one frame for the
    world's current state plus `req.num_frames - 1` more, each `req.steps_per_frame` real
    steps of `req.step_years` further along, with only the last of those steps rendered (see
    render_image.stream_animation_mp4 for why that's not the same as one bigger step). Each
    response line is one JSON object: `{"type": "progress", "frame": n, "total": N,
    "image_base64": <that frame's PNG>, "elapsed_years": <world.elapsed_years after this
    frame>, "stats": <latest World.stats_history snapshot, or null>}` as each frame finishes,
    then a final `{"type": "done",
    "video_base64": ..., "mime": "video/mp4", ...world summary fields}`. If rendering raises partway through, a `{"type": "error",
    "detail": ...}` line is emitted instead -- the HTTP status is already 200 by then, since
    the stream has started.

    **This permanently advances the world** by
    `(req.num_frames - 1) * req.steps_per_frame * req.step_years` years, the same as calling
    /worlds/{world_id}/step that many times -- not a side-effect-free preview. Same validation as
    /worlds/{world_id}/render for projection/view/width/height/rotation, plus `num_frames` bounded to
    `[1, MAX_ANIMATION_FRAMES]` and `steps_per_frame` required to be `>= 1` (each frame costs
    up to `steps_per_frame` full step_world calls + one render). `404` if no world has been
    generated yet, `503` if a generation, step, or animation is already in progress."""
    session = _session_or_404(world_id)
    if req.projection not in projections.PROJECTIONS:
        raise HTTPException(status_code=400, detail=f"unknown projection {req.projection!r}")
    if req.view not in render_image.VIEWS:
        raise HTTPException(status_code=400, detail=f"unknown view {req.view!r}; choices are {render_image.VIEWS}")
    if not (1 <= req.width <= MAX_RENDER_DIMENSION_PX and 1 <= req.height <= MAX_RENDER_DIMENSION_PX):
        raise HTTPException(status_code=400, detail=f"width/height must be in [1, {MAX_RENDER_DIMENSION_PX}]")
    if not (1 <= req.num_frames <= MAX_ANIMATION_FRAMES):
        raise HTTPException(status_code=400, detail=f"num_frames must be in [1, {MAX_ANIMATION_FRAMES}]")
    if not 1 <= req.steps_per_frame <= MAX_STEPS_PER_FRAME:
        raise HTTPException(
            status_code=400, detail=f"steps_per_frame must be in [1, {MAX_STEPS_PER_FRAME}]"
        )
    _require_valid_step_years(req.step_years, field="step_years")
    view_rotation = _parse_view_rotation(req.rotation)

    # Mark busy synchronously, not via a `with` block -- the response streams, so
    # the release has to happen in the generator's own `finally` instead. An overlapping
    # step/animation on this world still 503s before the stream starts.
    if not session.begin_operation():
        raise HTTPException(status_code=503, detail="a generation, step, or animation is already in progress")
    # Captured under the flag (see step): no replacer can swap session.world between
    # this read and the stream's turn at the gate.
    world = _require_world(session)

    # Clear any stale signal from a previous run (e.g. a stop request that arrived after that
    # run had already finished on its own) so this fresh run starts unsignaled.
    session.animation_stop_event.clear()

    def _stream():
        acquired = False
        _compute_gate.acquire()
        try:
            session.lock.acquire()
        except BaseException:
            # Lock.acquire() raises only on interpreter-level events, but if it
            # did, `acquired` would stay False while the gate is held -- a wedge
            # every world shares and no flag tracks. Release it first.
            _compute_gate.release()
            raise
        acquired = True
        try:
            # The world may have been deleted while this stream queued at the gate.
            # Bail out quietly rather than animating a detached session. Inside the
            # try on purpose (see generate): a pre-try return would skip the finally
            # and leak the gate, the session lock and the operation flag. Both flags
            # are initialized before the bail: the finally below reads them.
            last_frame = 0
            completed = False
            if session.deleted:
                return
            animation = render_image.stream_animation_mp4(
                world, req.projection, req.view, req.width, req.height, view_rotation,
                req.step_years, req.steps_per_frame, req.num_frames,
                step_fn=lambda w, y: _step_world(session, w, y), stop_event=session.animation_stop_event,
            )
            try:
                while True:
                    # Advance one frame while holding the locks, then let them go
                    # while the line goes out to the client: holding the global
                    # gate across up to 480 frames would starve every other world
                    # at the FIFO gate (and pin a worker thread for the whole
                    # run). The operation flag stays held throughout, so no other
                    # world replacer can swap session.world underneath; each frame
                    # ends on a whole step, so anything running between frames
                    # sees a consistent world.
                    try:
                        message = next(animation)
                    except StopIteration:
                        break
                    if message[0] == "progress":
                        _, frame, total, frame_png, elapsed_years = message
                        last_frame = frame
                        # Each frame ends on a whole step. Saved now, a run that ends badly later
                        # (a disconnect, a failed step or encoder) still keeps the frames before it.
                        _checkpoint_world(session)
                        line = json.dumps({
                            "type": "progress",
                            "frame": frame,
                            "total": total,
                            "image_base64": base64.b64encode(frame_png).decode("ascii"),
                            "elapsed_years": elapsed_years,
                            # Each frame's step_fn call already records a stats.py snapshot onto
                            # World.stats_history (see World.record_stats, called from step_world)
                            # -- ride the latest one along here so the Stats panel can update live
                            # instead of sitting frozen until the run finishes (GET stats
                            # would just wait at the compute gate for the whole run, see its own
                            # docstring). `None` if record_stats skipped recording this step (e.g.
                            # simulate_climate_biomes off leaves climate_cache empty) -- the
                            # frontend just keeps showing its last-known snapshot in that case.
                            "stats": world.stats_history[-1] if world.stats_history else None,
                        }) + "\n"
                    else:
                        completed = True
                        _, mp4_bytes, stopped_early = message
                        line = json.dumps({
                            "type": "done",
                            "mime": "video/mp4",
                            "video_base64": base64.b64encode(mp4_bytes).decode("ascii"),
                            "stopped_early": stopped_early,
                            **_summary(session, world),
                        }) + "\n"
                    session.lock.release()
                    _compute_gate.release()
                    acquired = False
                    try:
                        yield line
                    except GeneratorExit:
                        # The client disconnected mid-stream -- close the inner generator
                        # so its try/finally flushes and closes the encoder instead of
                        # leaking the native libx264 context (see stream_animation_mp4).
                        animation.close()
                        raise
                    _compute_gate.acquire()
                    try:
                        session.lock.acquire()
                    except BaseException:
                        _compute_gate.release()
                        raise
                    acquired = True
            except GeneratorExit:
                # A failure line (below) is yielded while holding the locks, so a
                # disconnect there lands here rather than in the loop's handler.
                animation.close()
                raise
            except Exception as exc:  # noqa: BLE001 -- status is already 200, so surface it as a data line
                _logger.exception(
                    "animation stream failed at frame %s/%s (projection=%s view=%s size=%sx%s steps_per_frame=%s)",
                    last_frame,
                    req.num_frames,
                    req.projection,
                    req.view,
                    req.width,
                    req.height,
                    req.steps_per_frame,
                )
                yield json.dumps({"type": "error", "detail": f"{type(exc).__name__}: {exc}"}) + "\n"
        finally:
            if acquired:
                # Published once at the end; between frames other worlds (and this
                # world's own renders) get their turns at the gate instead of
                # waiting out the whole run.
                _world_changed(session, completed)
                session.lock.release()
                _compute_gate.release()
            elif last_frame:
                # Disconnected mid-run: the completed frames were checkpointed
                # individually, but without this the feed revision, render cache
                # and index would still describe the pre-animation world. Publish
                # them without re-checkpointing (completed=False skips the
                # redundant pickle). Non-blocking: disconnect finalization must
                # never park on the session lock.
                if session.lock.acquire(blocking=False):
                    try:
                        _world_changed(session, False)
                    finally:
                        session.lock.release()
            session.end_operation()

    return StreamingResponse(_stream(), media_type="application/x-ndjson")


@app.post("/worlds/{world_id}/animate/stop")
def animate_stop(world_id: str) -> dict:
    """Signals an in-progress POST /worlds/{world_id}/animate stream (the "Stop" toolbar button while
    recording -- see App.tsx) to end cleanly after its current frame: the stream still
    flushes the encoder and returns the video made so far, unlike the client just aborting
    the connection, which drops the video entirely. This is what makes "keep going until
    Stop is pressed" mode useful -- without it, stopping would lose the recording. A no-op
    (still 200) if no animation is currently running or it already finished; the signal is
    cleared at the start of each new animate() run, not here, so it can't leak into one."""
    _session_or_404(world_id).animation_stop_event.set()
    return {"status": "ok"}


@app.post("/worlds/{world_id}/controls")
def set_controls(world_id: str, req: ControlsRequest) -> dict:
    """Live-adjusts sea level, solar heat, ice-age frequency, the plate-movement/climate-biomes
    step toggles, the wind model, and/or the geomorphic-budget tuning multipliers
    (world.TUNING_MULTIPLIER_FIELDS -- rain/river/wind/ocean/glacier/seismic erosion,
    coastal leveling, river/ocean deposition, collision uplift rate & reach, volcanism,
    fault relief strength, plume production; all dimensionless, 1.0 == untuned) on the current world (the "Controls" window -- see
    World.sea_level_m/World.solar_multiplier/World.simulate_plate_movement/
    World.simulate_climate_biomes/World.wind_model and World's tuning-knob field group) and
    immediately recomputes world.climate_cache so the very next /worlds/{world_id}/render or
    /worlds/{world_id}/stats call reflects the change without waiting for a step -- climate_cache is
    otherwise only refreshed by erosion.py once per step (see World.climate_cache), which
    would make a real-time control feel laggy or unresponsive between steps. This one-time
    recompute happens regardless of the current simulate_climate_biomes setting -- it's a
    single call triggered by a UI action, not a per-step cost, so it's cheap even when the
    user has turned off per-step climate computation to speed up stepping. `400` for an
    unknown wind model or a negative multiplier; `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    if req.wind_model is not None and req.wind_model not in WIND_MODEL_CHOICES:
        raise HTTPException(
            status_code=400, detail=f"unknown wind_model {req.wind_model!r}; choices are {WIND_MODEL_CHOICES}"
        )
    if req.fault_deformation_mode is not None and req.fault_deformation_mode not in FAULT_DEFORMATION_MODE_CHOICES:
        raise HTTPException(
            status_code=400,
            detail=f"unknown fault_deformation_mode {req.fault_deformation_mode!r}; choices are {FAULT_DEFORMATION_MODE_CHOICES}",
        )
    if req.node_cloud_resample_mode is not None and req.node_cloud_resample_mode not in NODE_CLOUD_RESAMPLE_MODE_CHOICES:
        raise HTTPException(
            status_code=400,
            detail=(
                f"unknown node_cloud_resample_mode {req.node_cloud_resample_mode!r}; "
                f"choices are {NODE_CLOUD_RESAMPLE_MODE_CHOICES}"
            ),
        )
    tuning_updates = {name: getattr(req, name) for name in TUNING_MULTIPLIER_FIELDS if getattr(req, name) is not None}
    for name, value in tuning_updates.items():
        if not np.isfinite(value) or value < 0.0:
            raise HTTPException(status_code=400, detail=f"{name} must be finite and >= 0 (got {value})")
    with _compute_gate.hold(), session.lock:
        if req.sea_level_m is not None:
            # Sea level is eustatic now (see eustasy.py) -- the slider sets the conserved
            # ocean water volume to whatever floats the current hypsometry at this level,
            # rather than pinning sea_level_m directly (which the next step would re-solve away).
            eustasy.set_sea_level_via_water_budget(world, req.sea_level_m)
        if req.solar_multiplier is not None:
            world.solar_multiplier = req.solar_multiplier
        if req.ice_age_period_years is not None:
            world.ice_age_period_years = max(0.0, float(req.ice_age_period_years))
        if req.simulate_plate_movement is not None:
            world.simulate_plate_movement = req.simulate_plate_movement
        if req.simulate_climate_biomes is not None:
            world.simulate_climate_biomes = req.simulate_climate_biomes
        if req.wind_model is not None:
            world.wind_model = req.wind_model
        if req.fault_deformation_mode is not None:
            world.fault_deformation_mode = req.fault_deformation_mode
        if req.node_cloud_resample_mode is not None:
            world.node_cloud_resample_mode = req.node_cloud_resample_mode
        if req.debug_diagnostics is not None:
            world.debug_diagnostics = req.debug_diagnostics
        for name, value in tuning_updates.items():
            setattr(world, name, float(value))
        world.climate_cache = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density))
        # Sea level and climate change what the viewer's renders show.
        _world_changed(session)
    return {
        "sea_level_m": world.sea_level_m,
        "solar_multiplier": world.solar_multiplier,
        "ice_age_period_years": world.ice_age_period_years,
        "simulate_plate_movement": world.simulate_plate_movement,
        "simulate_climate_biomes": world.simulate_climate_biomes,
        "wind_model": world.wind_model,
        "fault_deformation_mode": world.fault_deformation_mode,
        "node_cloud_resample_mode": world.node_cloud_resample_mode,
        "debug_diagnostics": world.debug_diagnostics,
        **{name: getattr(world, name) for name in TUNING_MULTIPLIER_FIELDS},
    }


@app.get("/worlds/{world_id}/plates")
def list_plates(world_id: str) -> dict:
    """Every plate's outline + metadata (point count, bounding ellipse) as JSON, for the
    "Plate Inspector" map mode -- unlike /worlds/{world_id}/render, the client renders this itself
    interactively rather than receiving a baked PNG. Un-rotated/true-frame throughout (no
    `rotation` param): the client applies its current view rotation only at draw time, same
    philosophy as climate/render-grid geometry (see docs/simulation-model/rendering-and-inspectors.md#rotating-the-view).
    `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    with session.lock:
        world = _require_world(session)
        overlaps = _plate_overlaps(world)
        return {
            "elapsed_years": world.elapsed_years,
            "plates": [_plate_summary(p, world, overlaps) for p in world.plates],
        }


@app.get("/worlds/{world_id}/stats")
def get_stats(world_id: str) -> dict:
    """Aggregate physical/temperature/precipitation statistics for the Stats panel -- see
    stats.py. Stateless: recomputed fresh from the current world state every call. Separate
    from the actual recorded time series (see /worlds/{world_id}/stats_history / World.stats_history) --
    this is just the current instant, for a caller that only wants "right now" without paying
    for a full-history response. `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    with _compute_gate.hold(), session.lock:
        return stats.compute_stats(world)


@app.get("/worlds/{world_id}/stats_history")
def get_stats_history(world_id: str) -> dict:
    """The full recorded stats time series (see World.stats_history) -- one compute_stats
    snapshot per generate/step, persisted on `World` itself (and so in every saved .mbworld,
    see persistence.py) rather than only kept in the frontend's own React state. The frontend
    calls this once after a Load (see App.tsx's handleWorldReplaced) to restore the Stats
    panel's history charts from a save, instead of the pre-history-persistence behaviour of
    just resetting them to empty. `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    # Under the lock: record_stats appends to stats_history at the end of every step,
    # so transposing it mid-step would read a torn list.
    with session.lock:
        world = _require_world(session)
        history = list(world.stats_history)
    return _columnar_stats_history(history)


def _columnar_stats_history(history: list[dict]) -> dict:
    """Transpose recorded stats snapshots into per-field arrays for compact transfer."""
    fields = list(dict.fromkeys(key for snapshot in history for key in snapshot))
    columns: dict[str, object] = {}
    for field in fields:
        values = [snapshot.get(field) for snapshot in history]
        # The stats breakdowns are maps whose keys repeat in every snapshot. Transpose those
        # too, so biome names and ledger categories are sent once for the whole time series.
        if any(isinstance(value, dict) for value in values) and all(
            value is None or isinstance(value, dict) for value in values
        ):
            map_keys = list(dict.fromkeys(key for value in values if isinstance(value, dict) for key in value))
            columns[field] = {
                "keys": map_keys,
                "columns": {
                    key: [value.get(key) if isinstance(value, dict) and key in value else None for value in values]
                    for key in map_keys
                },
                # 0 = null/missing, 1 = empty map, 2 = non-empty map.
                "states": [0 if value is None else 2 if value else 1 for value in values],
            }
        else:
            columns[field] = values
    return {"format": "columnar-v1", "length": len(history), "columns": columns}


@app.get("/worlds/{world_id}/plate_at")
def plate_at(world_id: str, lat_deg: float, lon_deg: float) -> dict:
    """Which plate owns the node nearest (lat_deg, lon_deg) -- the Plate Inspector's click
    hit-test. The client unprojects its click (through whatever view rotation is currently
    active) to a true lat/lon before calling this, so this endpoint itself never needs to
    know about rotation. `400` for non-finite input, `404` if no world has been generated
    yet."""
    world = _require_world(_session_or_404(world_id))
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    return {"plate_id": plates.nearest_plate_id(world.plates, query_xyz)}


# A click near a genuinely removed point should register as "there was something here" even
# though the click itself is unlikely to land exactly on the old node's position -- this is
# the same order of magnitude as gaps.CLUSTER_RADIUS_MULT's own line-spacing multiple, just
# expressed as a flat angular radius since /worlds/{world_id}/node_at has no single plate's spacing to
# scale against (a removed point's own owning plate may no longer even exist).
_REMOVED_POINT_MATCH_RAD = 0.02  # ~127 km at PLANET_RADIUS_KM


@app.get("/worlds/{world_id}/node_at")
def node_at(world_id: str, lat_deg: float, lon_deg: float) -> dict:
    """The "Added/Removed Points" (`nodeAge`) debug view's click-to-inspect popup: the single
    live surface node nearest (lat_deg, lon_deg) -- its owning plate, plate-local
    `phi`/`theta`, elevation, and `node_created_years` -- plus, independently, the nearest
    `World.removed_points_log` entry if one sits within `_REMOVED_POINT_MATCH_RAD` of the
    click (a removed point has no live node to report instead, so this is reported alongside,
    not merged with, whatever live node was also found).

    Unlike `/worlds/{world_id}/sample_at` (which reads the climate grid), this reports the actual nearest
    node's own plate-local coordinates -- recovered by inverse-transforming its world position
    through its owning plate's frame, which round-trips exactly to the value the node was
    created with (see geometry.py's module docstring for the phi/theta convention). `node` is `null` if no world has any live nodes yet (shouldn't happen via the
    API, but see `plates.nearest_node_index`'s own docstring). `400` for non-finite input,
    `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    with session.lock:
        node_info = None
        idx = plates.nearest_node_index(world.plates, query_xyz)
        if idx is not None:
            points, elevation, owner = plates.collect_all_points(world.plates)
            plate_id = int(owner[idx])
            plate = next((p for p in world.plates if p.plate_id == plate_id), None)
            if plate is not None:
                phi, theta = geometry.xyz_to_latlon(geometry.to_local(plate.frame, points[idx]))
                created = plates.collect_all_node_created_years(world.plates)
                node_info = {
                    "plate_id": plate_id,
                    "phi": float(phi),
                    "theta": float(theta),
                    "elevation_m": float(elevation[idx]),
                    "node_created_years": float(created[idx]),
                }

        removed_info = None
        if world.removed_points_log:
            removed_xyz = np.array([point for point, _, _ in world.removed_points_log])
            dist_rad = geometry.angular_distance(removed_xyz, query_xyz)
            best = int(np.argmin(dist_rad))
            if dist_rad[best] < _REMOVED_POINT_MATCH_RAD:
                _, removed_years, removed_plate_id = world.removed_points_log[best]
                removed_info = {
                    "removed_years": float(removed_years),
                    "plate_id": int(removed_plate_id),
                    "distance_rad": float(dist_rad[best]),
                }

    return {"lat_deg": lat_deg, "lon_deg": lon_deg, "node": node_info, "removed": removed_info}


def _surface_point_summary(plate: plates.Plate, point_index: int) -> dict:
    """The Points-view payload for one node: its real plate-local coordinates, authoritative
    surface fields, and world position (`point_xyz`, for the map highlight)."""
    nodes = plate.surface_nodes("elevation", "node_created_years")
    point_index = max(0, min(point_index, len(nodes.world_xyz) - 1))
    phi, theta = geometry.xyz_to_latlon(nodes.local_xyz[point_index])
    return {
        "plate_id": plate.plate_id,
        "point": {
            "phi": float(phi),
            "theta": float(theta),
            "elevation_m": float(nodes.fields["elevation"][point_index]),
            "node_created_years": float(nodes.fields["node_created_years"][point_index]),
        },
        "point_xyz": [float(value) for value in nodes.world_xyz[point_index]],
    }


def _surface_node_summary(plate, flat_index: int) -> dict:
    """Representation-neutral inspector payload for one current surface node."""
    nodes = plate.surface_nodes("elevation", "node_created_years")
    if not 0 <= flat_index < len(nodes.world_xyz):
        raise IndexError("surface node index out of range")
    return {
        "plate_id": plate.plate_id,
        # Decimal strings, not JSON numbers: the hash word spans uint64 and would lose
        # precision in a browser's IEEE-754 Number before being sent back to the lookup API.
        "node_id": [str(int(word)) for word in nodes.node_ids[flat_index]],
        "topology_revision": plate.topology_revision,
        "geometry_revision": plate.geometry_revision,
        "local_xyz": [float(value) for value in nodes.local_xyz[flat_index]],
        "world_xyz": [float(value) for value in nodes.world_xyz[flat_index]],
        "area_m2": float(nodes.area_m2[flat_index]),
        "area_is_exact": nodes.area_is_exact,
        "elevation_m": float(nodes.fields["elevation"][flat_index]),
        "node_created_years": float(nodes.fields["node_created_years"][flat_index]),
    }


def _parse_node_id_word(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail="node ID words must be decimal uint64 strings") from exc
    if not 0 <= parsed < 2**64:
        raise HTTPException(status_code=400, detail="node ID words must be decimal uint64 strings")
    return parsed


@app.get("/worlds/{world_id}/surface_node")
def surface_node(world_id: str, plate_id: int, node_id_hi: str, node_id_lo: str) -> dict:
    """Look up a live node by its storage-neutral current surface address."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    with session.lock:
        plate = next((p for p in world.plates if p.plate_id == plate_id), None)
        if plate is None:
            raise HTTPException(status_code=404, detail=f"no plate {plate_id}")
        flat_index = plate.node_index_for_id(
            (_parse_node_id_word(node_id_hi), _parse_node_id_word(node_id_lo))
        )
        if flat_index is None:
            raise HTTPException(status_code=404, detail="node is not live on this plate")
        return _surface_node_summary(plate, flat_index)


@app.get("/worlds/{world_id}/elevation_point_at")
def elevation_point_at(world_id: str, lat_deg: float, lon_deg: float) -> dict:
    """The "Points" (`platesDetail`) debug view's click-to-inspect: the nearest terrain node
    (see `_surface_point_summary`). `400` for non-finite input, `404` if no world has been
    generated yet or no plate has any live nodes."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    with session.lock:
        idx = plates.nearest_node_index(world.plates, query_xyz)
        if idx is None:
            raise HTTPException(status_code=404, detail="world has no live nodes")
        points, _elevation, owner = plates.collect_all_points(world.plates)
        plate_id = int(owner[idx])
        plate = next((p for p in world.plates if p.plate_id == plate_id), None)
        if plate is None:
            raise HTTPException(status_code=404, detail="owning plate not found")
        tree = plate.get_node_kdtree()
        if tree is None:
            raise HTTPException(status_code=404, detail="owning plate has no live nodes")
        _distance, point_index = tree.query(points[idx])
        return _surface_point_summary(plate, int(point_index))


def _fault_summary(fault, plate, other_plate_tree) -> dict:
    """One fault's inspector payload -- trace in true world coords (rounded like
    `_plate_summary`'s points), plus its type/motion/age metadata and how far the trace
    midpoint currently sits from the nearest cross-plate boundary node."""
    trace = faults.fault_world_points(fault, plate)
    mid = trace[len(trace) // 2]
    if other_plate_tree is not None:
        dist_rad = float(other_plate_tree.query(mid)[0])
    else:
        dist_rad = 0.0
    lifespan = float(fault.lifespan_myr)
    return {
        "fault_id": fault.fault_id,
        "plate_id": fault.plate_id,
        "kind": fault.kind,
        "trace": _round_coords(trace),
        "active": bool(fault.active),
        # A boundary fault (see faults.generate_boundary_faults) is regenerated every step
        # tracking its plate-boundary segment -- no age / lifespan / scar.
        "boundary": bool(getattr(fault, "boundary", False)),
        "slip_rate_m_per_myr": round(float(fault.slip_rate_m_per_myr), 1),
        "cumulative_offset_m": round(float(fault.cumulative_offset_m), 1),
        "age_myr": round(float(fault.age_myr), 2),
        "lifespan_myr": None if not np.isfinite(lifespan) else round(lifespan, 2),
        "dip_deg": round(float(fault.dip_deg), 1),
        "length_km": round(float(fault.length_km()), 1),
        "set_id": fault.set_id,
        "system_id": fault.system_id,
        "birth_distance_from_boundary_km": round(float(fault.birth_distance_from_boundary_km), 1),
        "distance_from_boundary_km": round(float(dist_rad) * mantle.PLANET_RADIUS_KM, 1),
    }


def _fault_system_summary(system, plate) -> dict:
    """One fault system's inspector payload -- its master lineament in true world coords,
    plus the belt-scale metadata. The strand `Fault`s carry `system_id` and are returned in
    the same `faults` list."""
    master = faults.system_world_points(system, plate)
    return {
        "system_id": system.system_id,
        "plate_id": system.plate_id,
        "kind": system.kind,
        "trace": _round_coords(master),
        "active": bool(system.active),
        "length_km": round(float(system.length_km), 1),
        "age_myr": round(float(system.age_myr), 2),
        "lifespan_myr": round(float(system.lifespan_myr), 2),
    }


def _other_plate_trees(world: World) -> dict:
    """plate_id -> cKDTree over every *other* plate's nodes, so a fault's distance to the
    nearest cross-plate boundary is one query. Built once per /worlds/{world_id}/faults call."""
    collected = plates.collect_all_points(world.plates)
    if collected is None:
        return {}
    points, _, owner = collected
    trees: dict[int, cKDTree] = {}
    for pid in {p.plate_id for p in world.plates}:
        mask = owner != pid
        trees[pid] = cKDTree(points[mask]) if np.any(mask) else None
    return trees


@app.get("/worlds/{world_id}/faults")
def list_faults(world_id: str) -> dict:
    """Every intraplate fault line (see faults.py) as plain JSON, for the "Fault Line
    Inspector" map mode -- same philosophy as /worlds/{world_id}/plates: un-rotated/true-frame, the
    client renders and drives interaction itself. `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    with session.lock:
        by_id = {p.plate_id: p for p in world.plates}
        trees = _other_plate_trees(world)
        summaries = [
            _fault_summary(fault, by_id[fault.plate_id], trees.get(fault.plate_id))
            for fault in (*world.faults, *world.boundary_faults)
            if fault.plate_id in by_id
        ]
        systems = [
            _fault_system_summary(system, by_id[system.plate_id])
            for system in world.fault_systems
            if system.plate_id in by_id
        ]
        return {"elapsed_years": world.elapsed_years, "faults": summaries, "fault_systems": systems}


@app.get("/worlds/{world_id}/earthquakes")
def list_earthquakes(world_id: str) -> dict:
    """Recent earthquakes (see faults.Earthquake / faults._generate_earthquakes) as plain
    JSON, for the fading epicentre overlay in the "Fault lines" view. The backend prunes this
    list after faults.EARTHQUAKE_RETAIN_MYR, so it only ever holds the last few Myr. `404` if
    no world has been generated yet."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    with session.lock:
        elapsed = world.elapsed_years
        quakes = [
            {
                "earthquake_id": q.earthquake_id,
                "fault_id": q.fault_id,
                "plate_id": q.plate_id,
                "kind": q.kind,
                "epicenter": _round_coords(np.asarray(q.epicenter_world)),
                "magnitude": round(float(q.magnitude), 2),
                "age_myr": round((elapsed - q.birth_years) / 1e6, 3),
                "birth_years": q.birth_years,
                "trigger": getattr(q, "trigger", faults.TRIGGER_TECTONIC),
            }
            for q in world.earthquakes
        ]
        return {"elapsed_years": elapsed, "earthquakes": quakes}


@app.get("/worlds/{world_id}/volcanoes")
def list_volcanoes(world_id: str) -> dict:
    """Every current volcano node (see volcanism.py / elevation_lines.is_volcano) as plain
    JSON, for the "Plates & Faults" view's overlay -- un-rotated/true-frame world positions,
    the client projects and draws them itself. `active` is true while the node still has
    eruption potential (`volcano_active_years_remaining > 0`). `404` if no world yet."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    with session.lock:
        collected = plates.collect_all_points(world.plates)
        if collected is None:
            return {"elapsed_years": world.elapsed_years, "volcanoes": []}
        points, _elev, _owner = collected
        is_volcano = plates.collect_all_is_volcano(world.plates)
        active_remaining = plates.collect_all_volcano_active_years_remaining(world.plates)
        idx = np.flatnonzero(is_volcano)
        volcanoes = [
            {"position": _round_coords(points[i]), "active": bool(active_remaining[i] > 0)}
            for i in idx
        ]
        return {"elapsed_years": world.elapsed_years, "volcanoes": volcanoes}


@app.get("/worlds/{world_id}/fault_at")
def fault_at(world_id: str, lat_deg: float, lon_deg: float) -> dict:
    """Which fault trace is nearest (lat_deg, lon_deg), or `null` if none is within a small
    threshold -- the Fault Line Inspector's click hit-test. Same true-frame / finite-input
    contract as `/worlds/{world_id}/plate_at`."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    with session.lock:
        by_id = {p.plate_id: p for p in world.plates}
        threshold_rad = 120.0 / mantle.PLANET_RADIUS_KM
        best_id = None
        best_dist = threshold_rad
        for fault in world.faults:
            plate = by_id.get(fault.plate_id)
            if plate is None:
                continue
            trace = faults.fault_world_points(fault, plate)
            d = float(np.min(np.arccos(np.clip(trace @ query_xyz, -1.0, 1.0))))
            if d < best_dist:
                best_dist = d
                best_id = fault.fault_id
    return {"fault_id": best_id}


@app.get("/worlds/{world_id}/sample_at")
def sample_at(world_id: str, lat_deg: float, lon_deg: float) -> dict:
    """Everything the "Elevation & Biome" view's click-to-inspect popup shows for one point
    (also offered on the Elevation and Biome views): elevation, biome, precipitation,
    temperature, and which plate owns it. Elevation/biome/precipitation/temperature are read
    from the same climate grid `/worlds/{world_id}/render` and `/worlds/{world_id}/stats` already use
    (`climate.compute_climate_cached` -- the world's own `climate_density` resolution),
    picked by nearest grid cell; `temperature_c` is the same land-air / ocean-surface
    composite the temperature view and biome classification use (see
    `climate.compute_climate`). `biome` may differ from the rendered pixel right at a
    coastline -- the Combined/Biome *renders* re-run `biomes.smooth_biome_field` on a finer
    grid (see `render_image._biome_fields`) -- but this is the canonical per-cell field every
    other consumer reads, and it carries the same boundary-cleanup pass. `plate_id` is the same nearest-node hit-test as `/worlds/{world_id}/plate_at`. The
    client unprojects its click through the active view rotation to a true lat/lon first,
    same as `/worlds/{world_id}/plate_at`, so this endpoint never needs to know about rotation. `400`
    for non-finite input, `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    world = _require_world(session)
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    with _compute_gate.hold(), session.lock:
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


@app.get("/worlds/{world_id}/rivers")
def list_rivers(world_id: str) -> dict:
    """Every distinct river network's path + endpoint metadata, plus the current coastline, as
    JSON, for the "River Inspector" map mode -- same philosophy as /worlds/{world_id}/plates: the client
    renders and drives this itself rather than receiving a baked PNG (see
    docs/simulation-model/hydrology.md#hydrology and docs/simulation-model/hydrology.md#river-inspector). Rivers are grouped fresh from
    the current world.hydrology_cache (see hydrology.group_rivers) on every call rather than
    persisted -- hydrology_cache itself is already at most one step stale by design (see
    World.hydrology_cache), and grouping it is cheap. `coastline_segments` is computed the
    same way and drawn server-side in the temperature/humidity/precipitation views (see
    coastline.py) -- included here too since this view has no other land/ocean cue at all
    (unlike those, it draws no filled backdrop). Both are `[]` before the first step
    (climate_cache/hydrology_cache are None until erosion.py runs once). Un-rotated/
    true-frame throughout, same as /worlds/{world_id}/plates. `404` if no world has been generated yet."""
    world = _require_world(_session_or_404(world_id))
    if world.hydrology_cache is None:
        return {"elapsed_years": world.elapsed_years, "rivers": [], "coastline_segments": []}
    rivers = hydrology.group_rivers(world.hydrology_cache)
    return {
        "elapsed_years": world.elapsed_years,
        "rivers": [_river_summary(world.hydrology_cache, river, i) for i, river in enumerate(rivers)],
        "coastline_segments": _coastline_segments_json(world),
    }


@app.get("/worlds/{world_id}/river_at")
def river_at(world_id: str, lat_deg: float, lon_deg: float) -> dict:
    """Which river network owns the node nearest (lat_deg, lon_deg) -- the River Inspector's
    click hit-test, same pattern as /worlds/{world_id}/plate_at. `river_id` is an index into whatever
    /worlds/{world_id}/rivers returns for *this same step* (rivers are regrouped fresh every call, not
    given a persistent identity -- see hydrology.group_rivers), so the client should treat it
    as only meaningful against its most recent /worlds/{world_id}/rivers response. `400` for non-finite
    input, `404` if no world has been generated yet."""
    world = _require_world(_session_or_404(world_id))
    if not (np.isfinite(lat_deg) and np.isfinite(lon_deg)):
        raise HTTPException(status_code=400, detail="lat_deg/lon_deg must be finite")
    if world.hydrology_cache is None:
        return {"river_id": None}
    query_xyz = geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg))
    rivers = hydrology.group_rivers(world.hydrology_cache)
    return {"river_id": hydrology.river_at(world.hydrology_cache, rivers, query_xyz)}


@app.get("/worlds/{world_id}/lakes")
def list_lakes(world_id: str) -> dict:
    """Every currently-visible lake's basin info, for the "Lake Inspector" map mode -- same
    "raw JSON, client renders it" philosophy as `/worlds/{world_id}/rivers`. A "lake" here means one
    connected component of `lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M` nodes (see
    `hydrology.lake_components`) -- exactly what's drawn as standing water everywhere else in
    this codebase -- resolved back up to its own basin in `fields.lake_forest`, this step's
    already-resolved depression hierarchy (see `HydrologyFields.lake_forest`'s own docstring
    for why this must be the *same* tree `lakes.step_lakes` used, not a separately rebuilt
    one), so a lake that's currently the result of several original catchments merging
    together reports its own *current* outlet, not some interior saddle that's already
    flooded. Regrouped fresh from `world.hydrology_cache` on
    every call, same as `/worlds/{world_id}/rivers` -- `lake_id` is only meaningful against this same
    response, not a stable identity across steps. `coastline_segments` included for the same
    reason `/worlds/{world_id}/rivers` includes it: this view draws no filled backdrop of its own. Both are
    `[]` before the first step. `404` if no world has been generated yet."""
    world = _require_world(_session_or_404(world_id))
    fields = world.hydrology_cache
    if fields is None:
        return {"elapsed_years": world.elapsed_years, "lakes": [], "coastline_segments": []}
    forest = fields.lake_forest
    components = _lake_components_sorted(fields)
    rivers = hydrology.group_rivers(fields)
    river_by_node = _river_by_node(rivers)
    result = []
    for lake_id, members in enumerate(components):
        lake_obj = _smallest_lake_containing(forest, set(members.tolist()))
        if lake_obj is None:
            continue  # defensive only -- every wet node belongs to some catchment by construction
        result.append(_lake_basin_summary(fields, lake_obj, lake_id, rivers, river_by_node, is_lake=True))
    return {"elapsed_years": world.elapsed_years, "lakes": result, "coastline_segments": _coastline_segments_json(world)}


@app.get("/worlds/{world_id}/lake_at")
def lake_at(world_id: str, lat_deg: float, lon_deg: float) -> dict:
    """The Lake Inspector's click hit-test, resolving to whichever *basin* owns the node
    nearest `(lat_deg, lon_deg)` -- not just a currently-wet lake, so a click on ordinary dry
    land still returns "the basin nonetheless" (the user's own spec), including its lowest
    point, as long as that land actually drains into some enclosed catchment at all (see
    `kind` below). `400` for non-finite input, `404` if no world has been generated yet.

    `kind` is one of:
    - `"lake"` -- the node is currently flooded; `basin` is that connected lake's own info
      (`lake_id` set, matching the same node's `/worlds/{world_id}/lakes` entry for this same step).
    - `"basin"` -- dry land, but part of a real (possibly still-unresolved/endorheic) basin;
      `basin` is that basin's own leaf catchment info (`lake_id` is `None` -- a dry basin has
      no visible-lake identity to be stable against).
    - `"no_basin"` -- dry land whose own steepest-descent chain drains straight to the ocean
      without ever passing through a local minimum first -- an ordinary hillslope, never part
      of any basin (see lakes.py's own `_OCEAN_CATCHMENT`). `basin` is `null`.
    - `"ocean"` -- the nearest node is open ocean. `basin` is `null`.
    """
    world = _require_world(_session_or_404(world_id))
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
    river_by_node = _river_by_node(rivers)
    if fields.lake_depth[node_idx] > hydrology.LAKE_MIN_VISIBLE_DEPTH_M:
        components = _lake_components_sorted(fields)
        lake_id, members = next((i, m) for i, m in enumerate(components) if node_idx in m.tolist())
        lake_obj = _smallest_lake_containing(forest, set(members.tolist())) or leaf
        basin = _lake_basin_summary(fields, lake_obj, lake_id, rivers, river_by_node, is_lake=True)
        return {"kind": "lake", "basin": basin}

    basin = _lake_basin_summary(fields, leaf, None, rivers, river_by_node, is_lake=False)
    return {"kind": "basin", "basin": basin}


@app.get("/worlds/{world_id}/stranded_basins")
def list_stranded_basins(world_id: str) -> dict:
    """Endorheic depressions whose floor sits below sea level and that have no drainage path
    to the ocean at all -- the "land-locked coastal pit" pathology from GitHub issue #122's
    coastal-speckle section, which the event log records but drowns in near-sea-level
    transient-pond churn (see docs/debugging/debug-ui.md). Read straight off this step's already-
    resolved depression hierarchy (`world.hydrology_cache.lake_forest`): a top-level basin
    with no known spill (`Lake.max_depth is None`) whose floor is below `world.sea_level_m`.
    `persisted_years`/`steps_seen` come from `world.stranded_basin_tracks`, the cross-step
    centroid-matched tracker `world.step_world` maintains (at most one step stale, same as
    `hydrology_cache` itself). `stranded_basins` is `[]` before the first step. `404` if no
    world has been generated yet."""
    world = _require_world(_session_or_404(world_id))
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


@app.post("/worlds/{world_id}/export_hexgrid")
def export_hexgrid(world_id: str, req: ExportHexGridRequest) -> dict:
    """The "File > Export Hex Grid" action: tiles the sphere into a geodesic-icosahedron
    hex/pentagon dome (independent of the plate simulation's own node cloud -- see
    geodesic.py) at `req.frequency`, samples the current world's elevation/ocean/biome onto
    each tile, and returns the whole tiling as JSON for an external application (see
    docs/hex-export-format.md for the file format and the neighbor-finding pseudocode its
    own `neighbor_ids`/`corner_vertex_ids` fields back). `400` if `frequency` isn't one of
    `geodesic.FREQUENCY_CHOICES`, `404` if no world has been generated yet."""
    session = _session_or_404(world_id)
    if req.frequency not in geodesic.FREQUENCY_CHOICES:
        raise HTTPException(status_code=400, detail=f"unknown frequency {req.frequency!r}; choices are {geodesic.FREQUENCY_CHOICES}")
    # Sampling reads elevation/ocean/biome arrays a step mutates in place, so hold the
    # session lock (and the gate, like the other kernel-backed reads); deleted worlds
    # 404 rather than exporting a detached session's leftovers.
    with _compute_gate.hold(), session.lock:
        if session.deleted:
            raise HTTPException(status_code=404, detail="no such world")
        return geodesic.export_hexgrid(_require_world(session), req.frequency)


# --- Public viewer --------------------------------------------------------------------------
# The logged-out viewer's routes (#373): read-only renders plus an allowance-limited step,
# answering to anyone. AuthMiddleware only guards /worlds/*, so these are the whole of what
# an anonymous client can do to a world. Each public route addresses one world that the
# owner has marked public (#563); private and unknown IDs both 404. See public_view.py.


class _SessionComputeLock:
    """Lock adapter for GuestStepCoordinator: takes the compute gate before
    the session lock (the global lock order) as one acquire/release pair."""

    def __init__(self, session: WorldSession) -> None:
        self._session = session

    def acquire(self, blocking: bool = True, timeout: float = -1) -> bool:
        gate_timeout = None if blocking and timeout == -1 else max(timeout, 0.0)
        if not _compute_gate.acquire(timeout=gate_timeout):
            return False
        # Non-blocking here even when the caller would wait: this runs inside the
        # coordinator's mutex, where parking would hold the compute gate (stalling
        # every other world) and pile up guest threads behind it.
        if not self._session.lock.acquire(blocking=False):
            _compute_gate.release()
            return False
        return True

    def release(self) -> None:
        self._session.lock.release()
        _compute_gate.release()


class PublicStepRequest(BaseModel):
    # The revision (`/public/worlds/{world_id}/world`'s `world.revision`) the viewer is looking at.
    revision: str


def _guest_allowance(request: Request, guest_demo: bool = False) -> tuple[auth.AuthConfig, int | None]:
    """(config, steps left). None means unlimited: the owner viewing /view isn't charged.

    The local viewer demo opts into the guest allowance so the guest UI can be previewed."""
    config = auth.load_config()
    if auth.is_admin(request, config) and not guest_demo:
        return config, None
    return config, max(0, public_view.GUEST_STEP_LIMIT - public_view.steps_used(request, config))


@app.get("/public/worlds")
def public_worlds() -> dict:
    """The public worlds a logged-out viewer may open (#563)."""
    worlds = []
    for session in _registry.list():
        # name/is_public are mutated under the session lock by patch_world; the
        # feed summary has its own mutex.
        with session.lock:
            summary = session.feed.summary()
            running = _session_running(session)
            # A generated world answers while it is public -- and while its run
            # loop is on even when private, so /view can watch a planet that is
            # currently stepping (#597). It drops back out when the run stops.
            if summary is not None and (session.is_public or running):
                worlds.append({"world_id": session.world_id, "name": session.name, "running": running, **summary})
    return {"worlds": worlds}


@app.get("/public/worlds/{world_id}/world")
async def public_world(world_id: str, request: Request, after: str | None = None, guest_demo: bool = False) -> dict:
    """What the viewer shows for one public world: its revision and age, and this
    browser's remaining step allowance.

    With `after`, a long poll (#374): while `after` is still the current revision, wait up to
    public_view.GUEST_POLL_WAIT_S for the world to change before answering. This is how every
    viewer follows steps that others (or the owner) take without computing anything itself.
    Async, so a waiting viewer holds no worker thread."""
    session = _public_session_or_404(world_id)
    if after is not None:
        await session.feed.wait_past(after, public_view.GUEST_POLL_WAIT_S)
    _, remaining = _guest_allowance(request, guest_demo)
    return {
        "world": session.feed.summary(),
        "running": _session_running(session),
        "steps_remaining": remaining,
        "step_limit": public_view.GUEST_STEP_LIMIT,
        "step_years": public_view.GUEST_STEP_YEARS,
    }


@app.get("/public/worlds/{world_id}/render")
def public_render(world_id: str, view: str = "elevation", width: int = 1100, height: int = 611) -> dict:
    """The admin render restricted to the viewer's fixed choices: Eckert IV, the default
    orientation, and the viewer's supported views. Each (view, size) is rendered once per
    world revision and then served from that world's RenderCache without touching the
    world. `503` if the owner holds the world for longer than GUEST_LOCK_WAIT_S."""
    session = _public_session_or_404(world_id)
    if view not in public_view.GUEST_VIEWS:
        raise HTTPException(status_code=400, detail=f"unknown view {view!r}; choices are {public_view.GUEST_VIEWS}")
    limit = public_view.GUEST_MAX_RENDER_PX
    if not (1 <= width <= limit and 1 <= height <= limit):
        raise HTTPException(status_code=400, detail=f"width/height must be in [1, {limit}]")
    key = (view, width, height)
    cached = session.public_renders.get(session.feed.revision(), key)
    if cached is not None:
        return cached
    if not _compute_gate.acquire(timeout=public_view.GUEST_LOCK_WAIT_S):
        raise HTTPException(status_code=503, detail="the world is busy; try again shortly")
    try:
        # Timed, not blocking: the owner can hold the session lock for minutes, and
        # this runs with the gate held, so parking here would stall every world.
        if not session.lock.acquire(timeout=public_view.GUEST_LOCK_WAIT_S):
            raise HTTPException(status_code=503, detail="the world is busy; try again shortly")
        try:
            # Look again: viewers that arrive together after a step queue here, and only
            # the first of them needs to render.
            revision = session.feed.revision()
            entry = session.public_renders.get(revision, key)
            if entry is None:
                image_base64 = render_image.render_png_base64(
                    _require_world(session), public_view.GUEST_PROJECTION, view, width, height, _parse_view_rotation(None), False, False
                )
                entry = {**session.feed.summary(), "image_base64": image_base64}
                session.public_renders.put(revision, key, entry)
        finally:
            session.lock.release()
    finally:
        _compute_gate.release()
    return entry


def _run_guest_step(session: WorldSession) -> None:
    # Runs with the coordinator's lock pair held. A delete that slipped in between the
    # coordinator's acquire and here stands down instead of stepping a detached world.
    completed = False
    try:
        if session.deleted:
            return
        for _ in _step_progress(session, _require_world(session), public_view.GUEST_STEP_YEARS):
            pass
        completed = True
    finally:
        _world_changed(session, completed)


@app.post("/public/worlds/{world_id}/step")
def public_step(world_id: str, req: PublicStepRequest, request: Request, response: Response, guest_demo: bool = False) -> dict:
    """Advances one public world by public_view.GUEST_STEP_YEARS, charging this browser's
    allowance only if this request is the one that computes it. The allowance is per
    browser across all public worlds, not per world. `computed` says whether it did; a
    stale `revision`, or a step another viewer already has running, is answered with the
    resulting world for free. `429` when a computation is needed and the allowance is
    spent, `503` if the owner is using the world, `404` for a private or unknown world."""
    session = _public_session_or_404(world_id)
    config, remaining = _guest_allowance(request, guest_demo)
    outcome = session.guest_steps.request_step(
        req.revision,
        may_compute=remaining is None or remaining > 0,
        current_revision=session.feed.revision,
        world_lock=_SessionComputeLock(session),
        compute=lambda: _run_guest_step(session),
    )
    if outcome.status == "exhausted":
        raise HTTPException(status_code=429, detail="step allowance used up")
    if outcome.status == "busy":
        raise HTTPException(status_code=503, detail="the world is busy; try again shortly")
    if outcome.status == "computed" and remaining is not None:
        remaining -= 1
        public_view.set_steps_used(response, public_view.GUEST_STEP_LIMIT - remaining, config)
    return {
        "world": session.feed.summary(),
        "computed": outcome.status == "computed",
        "steps_remaining": remaining,
    }


# --- Run status and checkpoint (#375, #563) -------------------------------------------------
# The owner's Play/Stop is a server-side status per world, stepped by that world's
# background thread (run_control.py), and each world and its status are checkpointed to
# persistent storage when one is configured (checkpoint.py), so a hosted instance that
# sleeps or restarts resumes them.
_checkpoint_root: checkpoint.CheckpointStore | None = None
# Why a configured checkpoint directory couldn't be used at startup, for /healthz to report.
_checkpoint_error: str | None = None


# Serializes _save_index's read-all-records + atomic write. Without it, a step's
# per-change save on one world interleaves with a rename/create/delete save on
# another and the last writer wins with a stale snapshot (lost update), and a
# save can catch a half-mutated record. Every mutation is followed by its own
# save, so serialized saves always converge on current memory. The registry lock
# inside list() stays a leaf (never held across a session lock), so this adds no
# lock-order edge. Deadlock audit: index -> registry only.
_index_lock = threading.Lock()


def _save_index() -> None:
    # Best-effort, like _checkpoint_world below: a failing index write (a full or
    # unmounted volume) is logged rather than raised, so it can never skip the
    # lock/gate releases in a stream's finally or wedge every world.
    if _checkpoint_root is None:
        return
    try:
        with _index_lock:
            _checkpoint_root.save_index([session.record() for session in _registry.list()])
    except Exception:  # noqa: BLE001 -- e.g. a full volume mustn't fail the request
        _logger.exception("checkpoint: couldn't save the registry index")


def _checkpoint_world(session: WorldSession) -> None:
    """Write `session`'s world to its checkpoint; the caller holds the session
    lock. A failed write is logged rather than raised: the world in memory is
    fine, and the next change retries."""
    world = session.world
    if session.checkpoint_store is None or world is None or session.deleted:
        return
    try:
        session.checkpoint_store.save_world(world)
    except Exception:  # noqa: BLE001 -- e.g. a full or unmounted volume mustn't fail the step
        _logger.exception("checkpoint: couldn't save the world")
        return
    if session.runner is not None:
        session.runner.retry_unsaved_status()


def _save_run_status(session: WorldSession, status: checkpoint.RunStatus) -> None:
    if session.checkpoint_store is not None and not session.deleted:
        session.checkpoint_store.save_status(status)


def _run_step_once(session: WorldSession, step_years: float, should_step) -> bool:
    """One step of `session`'s run loop. Waits at the compute gate rather than
    refusing like POST /worlds/{id}/step: the owner's own step or animation on
    any world just delays the run."""
    with _compute_gate.hold(), session.lock:
        if session.deleted:
            return True
        world = session.world
        if world is None:
            return False
        if not should_step():
            return True
        completed = False
        try:
            for _ in _step_progress(session, world, step_years):
                pass
            completed = True
        finally:
            _world_changed(session, completed)
    return True


def _wire_session(session: WorldSession) -> WorldSession:
    """Attach the per-session services that need main.py's helpers: the runner,
    the replay log, and (when checkpointing is on) the per-world store."""
    if session.runner is None:
        session.runner = run_control.WorldRunner(
            lambda years, should_step: _run_step_once(session, years, should_step),
            lambda status: _save_run_status(session, status),
            on_change=session.feed.wake,
        )
    if session.replay_log is None:
        session.replay_log = _make_replay_log(session.world_id)
    if _checkpoint_root is not None and session.checkpoint_store is None:
        session.checkpoint_store = _checkpoint_root.for_world(session.world_id)
    return session


def _restore_checkpoint() -> None:
    """At startup: load every checkpointed world and run status, if checkpointing
    is configured, and resume running where it was. Worlds are restored exactly
    as saved; the time the server spent asleep isn't simulated."""
    global _checkpoint_root, _checkpoint_error
    store = checkpoint.CheckpointStore.from_env()
    if store is None:
        _logger.info("checkpoint: %s isn't set; worlds won't survive a restart", checkpoint.CHECKPOINT_DIR_ENV)
        return
    _checkpoint_root = store
    try:
        records = store.load_index()
    except checkpoint.RegistryIndexError as exc:
        _logger.exception("checkpoint: can't use %s; running without checkpoints", store.directory)
        _checkpoint_error = str(exc)
        # Detach: later creates/renames/deletes must not _save_index() over the corrupt
        # worlds.json and destroy the evidence the fail-closed path preserved.
        _checkpoint_root = None
        return
    except OSError as exc:
        _logger.exception("checkpoint: can't use %s; running without checkpoints", store.directory)
        _checkpoint_error = f"can't use checkpoint directory {store.directory}: {exc}"
        _checkpoint_root = None
        return
    _checkpoint_error = None
    if records is None:
        _import_legacy_checkpoint(store)
        return
    # Most recently updated first: when the index holds more worlds than the cap
    # allows (a lowered cap, a tampered index), the least relevant ones are the
    # ones left unloaded below.
    by_recency = sorted(records, key=lambda r: float(r.get("updated_at") or 0.0), reverse=True)
    restored = 0
    for record in by_recency:
        # Reuse a live session with this ID when there is one (a restart in the
        # same process, as in tests): update it in place so references held by
        # callers stay valid.
        session = _registry.find(record["world_id"])
        if session is None:
            session = WorldSession(
                world_id=record["world_id"],
                name=record.get("name") or "Untitled world",
                is_public=bool(record.get("is_public")),
                created_at=float(record.get("created_at") or 0.0),
                updated_at=float(record.get("updated_at") or 0.0),
            )
            _registry.adopt(session)
        else:
            session.name = record.get("name") or session.name
            session.is_public = bool(record.get("is_public"))
        _wire_session(session)
        if restored >= _registry.max_worlds:
            # Memory bound: adopt() deliberately bypasses the create cap so saved
            # work is never lost, but eagerly loading every indexed world would let
            # an over-cap index (or a lowered cap) exhaust RAM and spawn unbounded
            # runners. The record stays visible and deletable, its files stay on
            # disk, and raising the cap + restarting resumes it.
            session.load_error = (
                f"not restored: the checkpoint holds more worlds than "
                f"MANTLE_BLOOM_MAX_WORLDS allows ({_registry.max_worlds}); "
                f"raise the cap and restart, or delete this record"
            )
            if session.world is not None:
                # Same-process restart (tests) reusing a live populated session:
                # don't keep serving a world that is over the cap.
                with session.lock:
                    session.world = None
                    session.feed.publish(None)
            continue
        restored += 1
        # Read before load(): an unreadable file makes load() return None *and* renames
        # the file aside, so checking existence only afterwards can't tell "corrupt"
        # from "never saved".
        had_world_file = session.checkpoint_store.world_path.exists()
        try:
            world, status = session.checkpoint_store.load()
        except OSError as exc:
            _logger.exception("checkpoint: couldn't load world %s", session.world_id)
            session.load_error = f"couldn't load the checkpointed world: {exc}"
            continue
        if world is not None:
            with session.lock:
                session.world = world
                session.feed.publish(world)
        elif had_world_file:
            # The record says a checkpoint should be here, but nothing usable came back:
            # checkpoint.py renamed the unreadable file aside with an `.invalid-` suffix
            # rather than overwriting it. Say so, instead of serving a silent empty world.
            session.load_error = (
                "the checkpointed world couldn't be loaded; the unreadable file "
                "was moved aside with an `.invalid-` suffix for inspection"
            )
        session.runner.set_status(status, persist=False)


def _import_legacy_checkpoint(store: checkpoint.CheckpointStore) -> None:
    """First startup after multi-world support: import the single root
    checkpoint as one public world, leaving the legacy files in place. The
    registry index takes precedence on later starts."""
    global _checkpoint_error, _checkpoint_root
    try:
        world, status = store.load()
    except OSError as exc:
        _logger.exception("checkpoint: can't use %s; running without checkpoints", store.directory)
        _checkpoint_error = f"can't use checkpoint directory {store.directory}: {exc}"
        _checkpoint_root = None
        return
    if world is None:
        return
    import uuid

    session = WorldSession(
        world_id=uuid.uuid4().hex,
        name=world_names.default_world_name({s.name for s in _registry.list()}),
        is_public=True,
    )
    _wire_session(session)
    _registry.adopt(session)
    with session.lock:
        session.world = world
        session.feed.publish(world)
    session.runner.set_status(status, persist=False)
    _save_index()
    _checkpoint_world(session)
    session.checkpoint_store.save_status(status)


# --- Replay log (#442, #563) ------------------------------------------------------------------
# Every step of a world goes through `_step_progress`, which, when
# MANTLE_BLOOM_REPLAY_LOG_DIR is set, saves the pre-step world of any step whose health
# metrics go askew (replay_log.py). Each world session gets its own log under
# `<dir>/<world_id>` so records and recent events never mix across worlds.
# How long shutdown waits for queued replay-log records (a 60 MB write each).
REPLAY_LOG_SHUTDOWN_FLUSH_S = 10.0


def _make_replay_log(world_id: str) -> replay_log.ReplayLog | None:
    directory = os.environ.get(replay_log.REPLAY_LOG_DIR_ENV, "").strip()
    if not directory:
        return None
    try:
        return replay_log.ReplayLog(Path(directory).expanduser() / world_id)
    except OSError:
        _logger.exception("replay log: can't use %s; running without it", directory)
        return None


def _step_progress(session: WorldSession, world: World, years: float):
    if session.replay_log is None:
        return step_world_progress(world, years)
    return session.replay_log.step(world, years)


def _step_world(session: WorldSession, world: World, years: float) -> None:
    for _ in _step_progress(session, world, years):
        pass


@app.get("/worlds/{world_id}/replay-logs")
def replay_logs(world_id: str) -> dict:
    """The replay-log records on disk for this world, newest first (`enabled:
    false` and none when MANTLE_BLOOM_REPLAY_LOG_DIR isn't set). Each names the
    triggers that fired and the record's directory, which
    `bin/debug/replay_logged_event.py` replays."""
    session = _session_or_404(world_id)
    if session.replay_log is None:
        return {"enabled": False, "directory": None, "records": []}
    return {"enabled": True, "directory": str(session.replay_log.directory), "records": session.replay_log.records()}


class RunRequest(BaseModel):
    running: bool
    # The UI's "Years per step". Omitted: keep the current one.
    step_years: float | None = None


def _run_summary(session: WorldSession) -> dict:
    runner = session.runner
    # The version first: anything that changes after it's read moves it again, so a poller
    # can't miss a change by having seen a newer status under an older version.
    run_version = runner.version()
    status = runner.status()
    return {
        "running": status.running,
        "step_years": status.step_years,
        "stepping": runner.stepping(),
        "error": runner.error(),
        "run_version": run_version,
        "world": session.feed.summary(),
    }


@app.get("/worlds/{world_id}/run")
async def get_run(world_id: str, after: str | None = None, run_version: int | None = None) -> dict:
    """The run status (`running`, `step_years`, `stepping`: a step is under way, which can
    still be true just after Stop, and `error`: why the run last stopped on its own, or null),
    its `run_version`, and the world's current revision. With `after`, a long poll like
    /public/worlds/{world_id}/world's that also answers when the run state moves past
    `run_version`: the full app follows the running world with it, and sees Stop pressed in
    another tab."""
    session = _session_or_404(world_id)
    if after is not None:
        runner = session.runner
        changed = None if run_version is None else (lambda: runner.version() != run_version)
        await session.feed.wait_past(after, public_view.GUEST_POLL_WAIT_S, changed)
    return _run_summary(session)


@app.post("/worlds/{world_id}/run")
def set_run(world_id: str, req: RunRequest) -> dict:
    """Play (`running: true`) or Stop. While running, the server steps the world by `step_years`
    after each run_control.RUN_INTERVAL_S pause, whether or not anyone is watching, and the
    status survives a restart when checkpointing is configured. Stop takes effect before the
    next step; a step already under way finishes. `400` for a non-positive `step_years`, `404`
    to start without a world, `500` if the status took effect but couldn't be saved (it's
    retried with the next world checkpoint)."""
    session = _session_or_404(world_id)
    runner = session.runner
    current = runner.status()
    step_years = current.step_years if req.step_years is None else req.step_years
    _require_valid_step_years(step_years, field="step_years")
    if req.running:
        _require_world(session)
    try:
        runner.set_status(checkpoint.RunStatus(running=req.running, step_years=float(step_years)))
    except run_control.StatusNotSaved as exc:
        # In effect now, but a restart would go back to the saved status (#375 review).
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    return _run_summary(session)


@app.get("/healthz", include_in_schema=False)
def healthz() -> dict:
    """Railway's deploy health check (.railway/railway.ts, #376). Public and lock-free: a step
    or render holding the world lock mustn't fail a deploy. It can't answer before the
    lifespan's checkpoint restore finishes. `503` when a checkpoint directory is configured
    but couldn't be used, so a deploy whose volume can't be written fails instead of running
    a world that the next sleep would lose."""
    if _checkpoint_error is not None:
        raise HTTPException(status_code=503, detail=_checkpoint_error)
    return {"status": "ok"}


# In the packaged single-process app and the Railway image (Dockerfile, #376) there is no
# separate Vite server: FastAPI serves the built React bundle itself, from the same origin as
# the API (so api.ts's API_BASE is ""). The dev workflow (bin/restart.sh) is unaffected -- it
# runs its own frontend server and never sets MANTLE_BLOOM_FRONTEND_DIST.


def mount_frontend(app: FastAPI, dist_dir: str | os.PathLike[str]) -> None:
    """Serve a `vite build` output directory (its `index.html` + `assets/`) at `/`.

    Mounted last, after every API route, so `/worlds/{world_id}/*` and `/docs` keep priority; anything
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

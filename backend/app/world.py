"""Top-level world state: a collection of plates, a mantle-flow field, and elapsed time."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import bathymetry, boundary, climate, erosion, gaps, geometry, hydrology, line_regrid, mantle, merge_split, reassign, volcanism
from .plates import DEFAULT_NODE_DENSITY, Plate, generate_plates

DEFAULT_MANTLE_CENTERS = 8
DEFAULT_AXIAL_TILT_DEG = 23.5
# Bounds how large World.events can grow over a long play session -- the UI's console only
# ever needs recent history, not an unbounded transcript.
MAX_EVENT_LOG_LENGTH = 200


@dataclass
class World:
    seed: int
    plates: list[Plate] = field(default_factory=list)
    mantle_centers: list[mantle.ConvectionCenter] = field(default_factory=list)
    elapsed_years: float = 0.0
    steps_since_regularize: int = 0
    # Separate from steps_since_regularize (and deliberately never triggered on the same
    # step -- see step_world) so reassign.py's whole-sphere neighbor query never runs against
    # a plate layout gaps.fill_gaps just changed out from under it, or vice versa.
    steps_since_reassign: int = 0
    next_plate_id: int = 0
    # A fixed per-world property, like `seed` -- set once at generation and read again on
    # every future climate render (see climate.py's compute_insolation), not rendering/cache
    # state. The one deliberate exception to climate being otherwise fully stateless.
    axial_tilt_deg: float = DEFAULT_AXIAL_TILT_DEG
    # Another fixed per-world property, set once at generation (see plates.generate_plates'
    # own node_density parameter) and read for the rest of this world's life by every module
    # that builds new elevation-line nodes or derives a distance/count threshold from
    # plates.TARGET_LINE_SPACING_RAD (line_regrid.py, boundary.py, gaps.py, merge_split.py,
    # volcanism.py -- see plates.line_spacing_rad's own docstring for why each of those needs
    # this rather than reading the bare module constant), so a world generated at a
    # non-default density stays self-consistent through regularizing/gap-filling/merging/
    # splitting/volcanism, not just at the moment it's generated.
    node_density: float = DEFAULT_NODE_DENSITY
    # Number of times gaps.fill_gaps has actually run -- part of the deterministic RNG seed
    # for gap-fill's new-crust noise texture (see gaps.py), not just a counter.
    gap_fill_calls: int = 0
    # Sustained-collision tracking for merge_split.py: (plate_id, plate_id) -> accumulated
    # convergent years. See merge_split.update_collision_progress.
    collision_progress: dict[tuple[int, int], float] = field(default_factory=dict)
    # plate_ids currently tracked as an active volcanic field -- removed (and relabeled as an
    # ordinary continental plate) once fewer than volcanism.VOLCANO_FRACTION_DORMANT_
    # THRESHOLD of the plate's own nodes are still is_volcano. See volcanism.py.
    volcanic_field_plate_ids: set[int] = field(default_factory=set)
    # Number of times volcanism.detect_and_spawn_volcanic_fields has actually run -- part of
    # the deterministic RNG seed for that pass (see volcanism.py), same role
    # World.gap_fill_calls plays for gaps.py's own new-crust noise texture.
    volcanic_field_calls: int = 0
    # Human-readable log for the UI's event console, each entry (elapsed_years, message).
    events: list[tuple[float, str]] = field(default_factory=list)
    # This step's climate snapshot (see climate.py), populated by erosion.py -- which needs
    # a fresh one every step regardless -- and reused by /world/stats and a climate map
    # render so they don't each trigger their own (~50ms) recomputation the same turn. See
    # climate.compute_climate_cached and climate.py's own module docstring for why reusing
    # a value that's up to one step stale is an accepted simplification here, not a bug.
    climate_cache: climate.ClimateFields | None = None
    # This step's flow-routing snapshot (see hydrology.py), populated by erosion.py
    # alongside climate_cache -- same reuse pattern, same one-step-stale simplification.
    hydrology_cache: hydrology.HydrologyFields | None = None

    def log_event(self, message: str) -> None:
        self.events.append((self.elapsed_years, message))
        if len(self.events) > MAX_EVENT_LOG_LENGTH:
            del self.events[: len(self.events) - MAX_EVENT_LOG_LENGTH]


def _plate_sample_points(plate: Plate) -> np.ndarray:
    """World positions used to sample the mantle flow field for fitting this plate's Euler
    pole: every elevation-line node -- its current footprint, for free (no separate
    sampling grid needed). This already includes each line's two endpoints, i.e. the
    plate's actual edge, so no separate boundary sample is needed."""
    points, _ = plate.all_points_and_elevation()
    return points


def _update_plate_omega(
    plate: Plate, centers: list[mantle.ConvectionCenter], damping: float | None
) -> None:
    sample_points = _plate_sample_points(plate)
    if len(sample_points) == 0:
        return
    velocities = mantle.flow_at(sample_points, centers)
    target_omega = mantle.fit_euler_pole(sample_points, velocities)
    if damping is None:
        new_omega = target_omega
    else:
        new_omega = plate.omega + damping * (target_omega - plate.omega)
    plate.omega = mantle.clamp_rate(new_omega)


def generate_world(
    seed: int,
    num_plates: int | None = None,
    continental_fraction: float | None = None,
    land_fraction: float | None = None,
    num_mantle_centers: int = DEFAULT_MANTLE_CENTERS,
    axial_tilt_deg: float | None = None,
    node_density: float = DEFAULT_NODE_DENSITY,
) -> World:
    """`num_plates` is optional -- see plates.generate_plates for why: the world tiles
    itself into a plausible number of plates rather than requiring the caller to pick one.
    `continental_fraction`/`land_fraction` are the UI's generation sliders -- see
    plates.generate_plates. `axial_tilt_deg` is the UI's third generation slider (degrees,
    defaults to DEFAULT_AXIAL_TILT_DEG = Earth's real tilt) -- doesn't affect plate
    generation at all, only climate.py's insolation, read at render time long after
    generation, which is why it's stored on World rather than consumed once here.
    `node_density` (the UI's "point density" choice) similarly affects only generation
    itself directly, but -- unlike axial_tilt_deg -- is stored on World because every later
    step also needs it (see World.node_density's own comment)."""
    plates = generate_plates(
        seed, num_plates=num_plates, continental_fraction=continental_fraction, land_fraction=land_fraction, node_density=node_density
    )
    # Separate RNG stream so changing num_mantle_centers doesn't reshuffle plate layout.
    mantle_rng = np.random.default_rng(seed + 1)
    mantle_centers = mantle.generate_convection_centers(mantle_rng, n_centers=num_mantle_centers)

    world = World(
        seed=seed,
        plates=plates,
        mantle_centers=mantle_centers,
        elapsed_years=0.0,
        next_plate_id=len(plates),
        axial_tilt_deg=DEFAULT_AXIAL_TILT_DEG if axial_tilt_deg is None else axial_tilt_deg,
        node_density=node_density,
    )
    for plate in world.plates:
        _update_plate_omega(plate, world.mantle_centers, damping=None)

    n_continents = sum(1 for p in plates if p.crust_type == "continental")
    world.log_event(f"World generated with {len(plates)} plates ({n_continents} continental).")
    return world


def step_world(world: World, years: float) -> None:
    """Advance the world by `years`: refit each plate's Euler pole from the mantle flow
    field (damped toward the new target), rotate the plate rigidly by that pole for
    `years` (exact for every carried point, no resampling), then let boundaries evolve --
    uplift/trench/ridge/rift elevation deltas and line growth/shrinkage where plates are
    now close to each other. Then topology changes: fully-subducted plates disappear,
    colliding continental plates merge (at most one per step, only after a sustained
    50-100 Myr collision -- see merge_split.py), and plates whose flow field no longer fits
    one rigid rotation well can split; any resulting events are logged to world.events for
    the UI's console. Every `line_regrid.REGULARIZE_INTERVAL_STEPS` calls, also fills any
    sphere-coverage gaps (a plate growing toward its own pole, or territory a subducted
    plate left unclaimed -- see gaps.py, which now logs each newly spawned plate) and
    regularizes any line whose interior spacing has drifted (line_regrid.py). On steps that
    aren't a regularize/gap-fill step, also periodically reassigns any node that's ended up
    geometrically embedded in a neighboring plate's own territory (see reassign.py) --
    deliberately never the same step as the gap-fill pass above, so each one's picture of
    "which plate owns what" stays current when it acts. Every step also erodes elevation
    based on the world's current climate (see erosion.py) -- rain/sheet erosion and
    weathering, the other half of the weather<->geology coupling from climate.py's own
    terrain-influences-weather mechanics (lapse rate, mountain wind deflection, orographic
    rain shadow) -- relaxes submerged continental crust toward a shelf-or-deep-water target
    based on distance to the nearest land (see bathymetry.py), and rolls each active
    volcano's own eruption chance (see volcanism.py). On the same cadence as the gap-fill/
    regularize pass, also scans for divergent gaps between plates and spawns any new
    volcanic fields they warrant (see volcanism.py)."""
    for plate in world.plates:
        _update_plate_omega(plate, world.mantle_centers, damping=mantle.VELOCITY_DAMPING)
        increment = geometry.rotation_matrix_from_omega(plate.omega, years)
        plate.frame = increment @ plate.frame
    boundary.step_boundaries(world, years)
    world.elapsed_years += years
    for message in merge_split.apply_topology_changes(world, years):
        world.log_event(message)

    erosion.apply_erosion(world, years)
    bathymetry.apply_bathymetry(world, years)
    for message in volcanism.apply_volcanic_activity(world, years):
        world.log_event(message)

    world.steps_since_regularize += 1
    run_regularize_this_step = world.steps_since_regularize >= line_regrid.REGULARIZE_INTERVAL_STEPS
    if run_regularize_this_step:
        for message in volcanism.detect_and_spawn_volcanic_fields(world):
            world.log_event(message)
        for message in gaps.fill_gaps(world):
            world.log_event(message)
        line_regrid.regularize_world_lines(world)
        world.steps_since_regularize = 0

    world.steps_since_reassign += 1
    if not run_regularize_this_step and world.steps_since_reassign >= reassign.REASSIGN_INTERVAL_STEPS:
        reassign.reassign_misplaced_points(world)
        world.steps_since_reassign = 0

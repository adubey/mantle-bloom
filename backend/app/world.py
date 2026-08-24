"""Top-level world state: a collection of plates, a mantle-flow field, and elapsed time."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from . import bathymetry, boundary, climate, elevation_lines, erosion, gaps, geology, geometry, hydrology, mantle, merge_split, reassign, volcanism
from .elevation_lines import DEFAULT_NODE_DENSITY
from .plates import Plate, gather_node_positions, generate_plates, query_workers

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
    # plates.TARGET_LINE_SPACING_RAD (elevation_lines.py, boundary.py, gaps.py, merge_split.py,
    # volcanism.py -- see plates.line_spacing_rad's own docstring for why each of those needs
    # this rather than reading the bare module constant), so a world generated at a
    # non-default density stays self-consistent through regularizing/gap-filling/merging/
    # splitting/volcanism, not just at the moment it's generated.
    node_density: float = DEFAULT_NODE_DENSITY
    # Another fixed per-world property, set once at generation (the UI's "climate & biome
    # resolution" choice, see climate.CLIMATE_DENSITY_CHOICES) and read for the rest of this
    # world's life by every caller that computes climate (erosion.py every step,
    # climate.compute_climate_cached for render_image.py/stats.py, main.py's /world/controls)
    # or the Biome/Combined/Resources/Soil-Quality views' own finer render grid
    # (render_image.biome_grid_dimensions) -- see climate.grid_dimensions' own docstring. Not
    # purely cosmetic: erosion.py samples precipitation/wind/humidity/temperature from this
    # same grid every step, so a finer grid resolves orographic rain shadow and mountain wind
    # deflection more precisely, which can subtly change simulated erosion outcomes too, not
    # just how smoothly the climate/biome maps render -- still fundamentally a resolution
    # choice, not a different climate model, but not entirely inconsequential to physics
    # either the way a pure render-quality setting would be.
    climate_density: float = climate.DEFAULT_CLIMATE_DENSITY
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
    # FIFO round-robin queue for reassign.run_round_robin_check: one plate's own nodes get
    # audited for territory outliers per step_world call, cycling continuously rather than in
    # a periodic whole-world burst (unlike reassign_misplaced_points below, which this
    # complements, not replaces). Reconciled against world.plates (new ids appended, dead ids
    # dropped) at the top of every run_round_robin_check call, so nothing else needs to keep
    # this in sync at every plate-creation/removal site. merge_split.merge_plates jumps the
    # merge survivor to the front on creation -- see that function's own comment.
    plate_check_queue: list[int] = field(default_factory=list)
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
    # Nearest-land spatial index backing distance_from_land_approx (below) -- reset to None
    # once per step (step_world, right where node_cloud is gathered, since land nodes' own
    # world positions and elevations can both have changed since the last build) and rebuilt
    # lazily on first use after that. None means "not built yet this step," not "no land" --
    # distance_from_land_approx itself is what tells those two cases apart.
    land_kdtree_cache: cKDTree | None = None
    # Live-adjustable via POST /world/controls (see main.py) for the UI's "Controls" window
    # -- unlike axial_tilt_deg/node_density (fixed at generation), these are meant to be
    # tweaked mid-simulation. sea_level_m replaces the bare `elevation <= 0.0` convention
    # every is_ocean check in this codebase used to hardcode (climate.py, hydrology.py,
    # bathymetry.py); solar_multiplier scales climate.py's own SUNLIGHT constant. Changing
    # either forces an immediate climate_cache recompute (see main.py's controls route) so
    # /world/render and /world/stats reflect it right away, without waiting for a step.
    sea_level_m: float = 0.0
    solar_multiplier: float = 1.0
    # Live-adjustable via POST /world/controls, same pattern as sea_level_m/solar_multiplier
    # above -- the UI's "Controls" window lets the user run *just* plate tectonics or *just*
    # climate & biomes. When False, step_world skips plate rotation, boundary evolution
    # (uplift/trench/ridge/rift), topology changes (merge/split), volcanism, and the periodic
    # gap-fill/regularize/reassign passes entirely -- world.plates is otherwise frozen, so
    # climate & biomes (if simulate_climate_biomes is still True) can be watched evolving on a
    # static landscape. elapsed_years still advances either way -- these two flags gate *what*
    # a step computes, not whether time passes.
    simulate_plate_movement: bool = True
    # Companion to simulate_plate_movement above, same live-adjustable pattern. When False,
    # step_world skips climate/erosion/hydrology/bathymetry/resource-formation entirely -- by
    # far the most expensive part of a step (climate.py's grid computation and hydrology.py's
    # flow routing) -- leaving climate_cache/hydrology_cache at whatever they were last
    # computed to (stale, the same one-step-behind tolerance World.climate_cache already
    # documents) rather than None, so a render/stats call right after toggling this off still
    # shows the last real climate snapshot instead of going blank.
    simulate_climate_biomes: bool = True

    def log_event(self, message: str) -> None:
        self.events.append((self.elapsed_years, message))
        if len(self.events) > MAX_EVENT_LOG_LENGTH:
            del self.events[: len(self.events) - MAX_EVENT_LOG_LENGTH]

    def distance_from_land_approx(self, points: np.ndarray) -> np.ndarray:
        """Approximate distance from each given world-xyz point (shape (n, 3)) to the
        nearest land node (elevation > sea_level_m) anywhere in this world -- lazily builds
        and reuses land_kdtree_cache (see its own docstring) off every plate's own public
        Plate.map_world_points()/ElevationPoint.get_elevation() interface, rather than
        reaching into any one representation's own storage the way gather_node_positions
        does (PlateWithLines-only). bathymetry.py uses this in place of building its own
        land-only tree from scratch every call. np.inf for every point if this world has no
        land anywhere yet."""
        if self.land_kdtree_cache is None:
            self.land_kdtree_cache = _build_land_kdtree(self)
        if self.land_kdtree_cache is None:
            return np.full(len(points), np.inf)
        if len(points) == 0:
            return np.zeros(0)
        dist, _ = self.land_kdtree_cache.query(points, workers=query_workers(len(points)))
        return dist


def _build_land_kdtree(world: World) -> cKDTree | None:
    """Every land node's (elevation > world.sea_level_m) world position, across every
    plate, gathered via Plate's own public map_world_points()/ElevationPoint.get_elevation()
    -- representation-agnostic, unlike plates.gather_node_positions -- and indexed for
    World.distance_from_land_approx. None if this world has no land nodes at all."""
    land_points = [
        world_xyz
        for plate in world.plates
        for point, world_xyz in plate.map_world_points()
        if point.get_elevation() > world.sea_level_m
    ]
    if not land_points:
        return None
    return cKDTree(np.array(land_points))


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
    plate.set_omega(mantle.clamp_rate(new_omega))


def generate_world(
    seed: int,
    num_plates: int | None = None,
    continental_fraction: float | None = None,
    land_fraction: float | None = None,
    num_mantle_centers: int = DEFAULT_MANTLE_CENTERS,
    axial_tilt_deg: float | None = None,
    node_density: float = DEFAULT_NODE_DENSITY,
    initial_soil_maturity: float | None = None,
    climate_density: float = climate.DEFAULT_CLIMATE_DENSITY,
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
    step also needs it (see World.node_density's own comment). `initial_soil_maturity` (the
    UI's "initial soil maturity" slider, 0 to 1, default None -> 0.0/barren) is a one-time
    generation-time seed -- like continental_fraction/land_fraction, not stored on World,
    since nothing later needs to know what it was (see geology.seed_initial_soil).
    `climate_density` (the UI's "climate & biome resolution" choice, see
    climate.CLIMATE_DENSITY_CHOICES) doesn't affect plate generation at all, only how finely
    climate.py's own grid resolves the world's climate every future step/render -- stored on
    World for the same reason node_density is (see World.climate_density's own comment)."""
    plates = generate_plates(
        seed, num_plates=num_plates, continental_fraction=continental_fraction, land_fraction=land_fraction, node_density=node_density
    )
    geology.seed_initial_soil(plates, seed, initial_soil_maturity or 0.0)
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
        climate_density=climate_density,
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
    the UI's console. Every `elevation_lines.REGULARIZE_INTERVAL_STEPS` calls, also fills any
    sphere-coverage gaps (a plate growing toward its own pole, or territory a subducted
    plate left unclaimed -- see gaps.py, which now logs each newly spawned plate) and
    regularizes any line whose interior spacing has drifted (elevation_lines.py). On steps that
    aren't a regularize/gap-fill step, also periodically reassigns any node that's ended up
    geometrically embedded in a neighboring plate's own territory (see reassign.py) --
    deliberately never the same step as the gap-fill pass above, so each one's picture of
    "which plate owns what" stays current when it acts. Every step, regardless of that
    cadence, also pops one plate off World.plate_check_queue and self-audits just its own
    nodes for the same kind of territory outlier -- a continuous, one-plate-at-a-time
    complement to the periodic whole-world reassign pass, giving fast, spread-out coverage
    rather than a periodic burst (see reassign.run_round_robin_check). Every step also erodes
    elevation based on the world's current climate (see erosion.py) -- rain/sheet erosion and
    weathering, the other half of the weather<->geology coupling from climate.py's own
    terrain-influences-weather mechanics (lapse rate, mountain wind deflection, orographic
    rain shadow) -- relaxes submerged continental crust toward a shelf-or-deep-water target
    based on distance to the nearest land (see bathymetry.py), and rolls each active
    volcano's own eruption chance, growing mineral_deposit_m wherever one erupts (see
    volcanism.py). Right after that, grows/relaxes soil and coal/oil-gas deposits from this
    same step's erosion/flow-routing results (see geology.py). On the same cadence as the
    gap-fill/regularize pass, also scans for divergent gaps between plates and spawns any new
    volcanic fields they warrant, fuses any currently-tracked volcanic fields that have ended
    up close to each other (bridging the gap between them with new interpolated land/
    volcanoes), and grows any volcanic field whose own edge has no other plate nearby to bound
    it (see volcanism.py).

    Both halves above are individually skippable, live, via World.simulate_plate_movement/
    World.simulate_climate_biomes (see their own docstrings and main.py's /world/controls) --
    "plate movement" below means rotation, boundary evolution, topology changes, volcanism,
    and the periodic gap-fill/regularize/reassign passes; "climate & biomes" means erosion
    (which itself computes this step's climate.py fields), hydrology, bathymetry, and
    geology's resource formation. elapsed_years always advances regardless of either flag."""
    if world.simulate_plate_movement:
        for plate in world.plates:
            _update_plate_omega(plate, world.mantle_centers, damping=mantle.VELOCITY_DAMPING)
            increment = geometry.rotation_matrix_from_omega(plate.omega, years)
            plate.rotate(increment)
        boundary.step_boundaries(world, years)
    world.elapsed_years += years
    if world.simulate_plate_movement:
        for message in merge_split.apply_topology_changes(world, years):
            world.log_event(message)

    erosion_result = None
    if world.simulate_climate_biomes:
        # Gathered once here (rather than independently inside erosion.apply_erosion/
        # bathymetry.apply_bathymetry) since node positions are fixed for the rest of this
        # step -- nothing between here and the next step's rotation moves a node or changes
        # line topology (only elevation and other per-node fields still change, which each of
        # climate.py/erosion.py/hydrology.py/bathymetry.py still reads fresh off world.plates
        # itself) -- see plates.gather_node_positions's own docstring for why this was worth
        # factoring out. Skipped entirely, alongside erosion/bathymetry below, when
        # simulate_climate_biomes is off -- no plate-movement module reads climate/hydrology
        # output (see World.simulate_climate_biomes), so there's nothing here for them to miss.
        node_cloud = gather_node_positions(world.plates)
        # Same "fixed for the rest of this step" reasoning as node_cloud above applies to
        # land_kdtree_cache -- reset here, then lazily rebuilt on whichever of erosion (which
        # can move the coastline) or bathymetry (see World.distance_from_land_approx) reads
        # it first.
        world.land_kdtree_cache = None
        erosion_result = erosion.apply_erosion(world, years, node_cloud=node_cloud)
        bathymetry.apply_bathymetry(world, years, node_cloud=node_cloud)
    if world.simulate_plate_movement:
        for message in volcanism.apply_volcanic_activity(world, years):
            world.log_event(message)
    if erosion_result is not None:
        geology.apply_resource_formation(world, years, erosion_result)

    world.steps_since_regularize += 1
    run_regularize_this_step = (
        world.simulate_plate_movement and world.steps_since_regularize >= elevation_lines.REGULARIZE_INTERVAL_STEPS
    )
    if run_regularize_this_step:
        for message in volcanism.detect_and_spawn_volcanic_fields(world):
            world.log_event(message)
        for message in volcanism.merge_close_volcanic_fields(world):
            world.log_event(message)
        volcanism.grow_isolated_volcanic_fields(world)
        for message in gaps.fill_gaps(world):
            world.log_event(message)
        elevation_lines.regularize_world_lines(world)
        world.steps_since_regularize = 0

    world.steps_since_reassign += 1
    if (
        world.simulate_plate_movement
        and not run_regularize_this_step
        and world.steps_since_reassign >= reassign.REASSIGN_INTERVAL_STEPS
    ):
        reassign.reassign_misplaced_points(world)
        world.steps_since_reassign = 0

    # Continuous, one-plate-per-step complement to the periodic whole-world pass just above --
    # see reassign.run_round_robin_check's own docstring for why this doesn't need the same
    # "never the same step as gap-fill" mutual exclusion (it only ever touches one plate
    # against a fresh snapshot taken at call time, not a whole-world batch).
    if world.simulate_plate_movement:
        reassign.run_round_robin_check(world)

"""Top-level world state: a collection of plates, a mantle-flow field, and elapsed time."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from . import atmosphere_cfd, bathymetry, climate, erosion, geology, hydrology, mantle, merge_split, ocean_cfd, volcanism
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
    next_plate_id: int = 0
    # A fixed per-world property, like `seed` -- set once at generation and read again on
    # every future climate render (see climate.py's compute_insolation), not rendering/cache
    # state. The one deliberate exception to climate being otherwise fully stateless.
    axial_tilt_deg: float = DEFAULT_AXIAL_TILT_DEG
    # Another fixed per-world property, set once at generation (see plates.generate_plates'
    # own node_density parameter) and read for the rest of this world's life by every module
    # that builds new elevation-line nodes or derives a distance/count threshold from
    # elevation_lines.TARGET_LINE_SPACING_RAD (elevation_lines.py, plates.py's deform(),
    # merge_split.py -- see plates.line_spacing_rad's own docstring for why each of those
    # needs this rather than reading the bare module constant), so a world generated at a
    # non-default density stays self-consistent through deforming/merging/splitting, not
    # just at the moment it's generated.
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
    # Another fixed per-world property, set once at generation (the UI's "Fluid dynamics
    # resolution" Advanced-settings choice, same climate.CLIMATE_DENSITY_CHOICES set
    # climate_density itself uses) -- independent of climate_density, read only by
    # ocean_cfd.init_ocean_cfd/atmosphere_cfd.init_atmosphere_cfd to size *their own* grid.
    # Unlike climate_density, lowering this doesn't touch the Biome/Combined/Resources/
    # Soil-Quality render grid or erosion.py's own climate sampling at all -- it only trades
    # off how finely Ocean/Atmospheric Fluid Dynamics mode resolves currents/wind against how
    # many substeps (and so how long) each step_ocean_cfd/step_atmosphere_cfd call takes,
    # since CFL substep count scales with grid spacing (see fluid_dynamics.cfl_substeps) and
    # per-substep cost scales with cell count. Defaults to climate_density's own default so a
    # world generated without touching this new setting behaves exactly as before.
    fluid_density: float = climate.DEFAULT_CLIMATE_DENSITY
    # Sustained-collision tracking for merge_split.py: (plate_id, plate_id) -> accumulated
    # convergent years. See merge_split.update_collision_progress.
    collision_progress: dict[tuple[int, int], float] = field(default_factory=dict)
    # plate_ids currently tracked as an active volcanic field -- removed (and relabeled as an
    # ordinary continental plate) once fewer than volcanism.VOLCANO_FRACTION_DORMANT_
    # THRESHOLD of the plate's own nodes are still is_volcano. Nothing populates this set any
    # more (PlateWithLines.deform spawns overstretch volcanoes as new nodes on the same
    # plate's own existing line, not as a separate tracked plate) -- kept for now since
    # apply_volcanic_activity's own per-node eruption rolling (which reads `is_volcano`
    # directly off every plate's lines, independent of this set) still needs somewhere to
    # report a field "cooling" if that tracking comes back later. See volcanism.py.
    volcanic_field_plate_ids: set[int] = field(default_factory=set)
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
    # Which of the three mutually-exclusive top-level simulation modes is currently active --
    # the UI's "Mode" toggle (see main.py's POST /world/mode and
    # docs/simulation-model.md#ocean-atmospheric-fluid-dynamics). "tectonics_climate" (the
    # default) is everything above: step_world runs exactly as it always has, completely
    # unaware this field even exists. Switching to "ocean_cfd"/"atmosphere_cfd" freezes plate
    # tectonics and the climate/erosion model (real ocean/atmosphere fluid dynamics needs
    # steps of hours-to-days; tectonics needs steps of thousands-to-millions of years -- no
    # single "Step" can mean both), and hands stepping over to ocean_cfd.step_ocean_cfd/
    # atmosphere_cfd.step_atmosphere_cfd instead (see main.py's POST /world/step_fluid) against
    # a frozen snapshot of whatever elevation/climate this world had the moment the mode was
    # entered -- world.plates itself is never touched while either FD mode is active.
    fluid_mode: str = "tectonics_climate"
    # Populated by main.py's POST /world/mode when fluid_mode switches to "ocean_cfd"/
    # "atmosphere_cfd" (init_ocean_cfd/init_atmosphere_cfd -- a fresh climate/elevation
    # snapshot every time, but seeded with whatever wind/ocean-current state
    # remembered_wind_u/remembered_ocean_u (etc., above) still hold from an earlier FD
    # session, if any -- see those fields' own docstrings), and cleared again when it switches
    # away. None whenever the matching mode isn't the active one.
    ocean_cfd_state: ocean_cfd.OceanCFDState | None = None
    atmosphere_cfd_state: atmosphere_cfd.AtmosphereCFDState | None = None
    # Wind (eastward/northward, m/s) carried across a fluid-dynamics mode switch -- written
    # by atmosphere_cfd.remember_atmosphere_state whenever "atmosphere_cfd" (the only mode
    # that actually evolves wind) is left, read by ocean_cfd.init_ocean_cfd/
    # atmosphere_cfd.init_atmosphere_cfd as their own starting wind instead of a fresh
    # climate.compute_climate diagnostic. Cleared by step_world (see its own docstring)
    # whenever a "tectonics_climate" step actually moves plates or recomputes climate, since
    # a wind baseline from before that step no longer matches the terrain it came from. None
    # means there's nothing to resume yet.
    remembered_wind_u: np.ndarray | None = None
    remembered_wind_v: np.ndarray | None = None
    # Ocean current/sea-surface state carried the same way -- written by
    # ocean_cfd.remember_ocean_state whenever "ocean_cfd" is left, read by
    # ocean_cfd.init_ocean_cfd as its own starting current instead of starting at rest.
    # Cleared alongside remembered_wind_u/v above.
    remembered_ocean_u: np.ndarray | None = None
    remembered_ocean_v: np.ndarray | None = None
    remembered_ocean_eta: np.ndarray | None = None
    remembered_ocean_temperature_c: np.ndarray | None = None
    remembered_ocean_sediment_concentration: np.ndarray | None = None
    remembered_ocean_sediment_deposited_m: np.ndarray | None = None

    def log_event(self, message: str) -> None:
        self.events.append((self.elapsed_years, message))
        if len(self.events) > MAX_EVENT_LOG_LENGTH:
            del self.events[: len(self.events) - MAX_EVENT_LOG_LENGTH]

    def distance_from_land_approx(self, points: np.ndarray) -> np.ndarray:
        """Approximate distance from each given world-xyz point (shape (n, 3)) to the
        nearest land node (elevation > sea_level_m) anywhere in this world -- lazily builds
        and reuses land_kdtree_cache (see its own docstring) off every plate's own public
        Plate.map_world_points()/ElevationPoint.get_elevation() interface. bathymetry.py uses
        this in place of building its own land-only tree from scratch every call. np.inf for
        every point if this world has no land anywhere yet."""
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
    plate, gathered via Plate's own public all_points_and_elevation() (a bulk per-plate
    array read, not a per-node point object) and indexed for World.distance_from_land_approx.
    None if this world has no land nodes at all."""
    chunks = []
    for plate in world.plates:
        points, elevation = plate.all_points_and_elevation()
        if len(points) == 0:
            continue
        chunks.append(points[elevation > world.sea_level_m])
    land_points = np.concatenate(chunks, axis=0) if chunks else np.zeros((0, 3))
    if len(land_points) == 0:
        return None
    # balanced_tree=False/compact_nodes=False -- same build-time speedup hydrology.py's/
    # plates.py's own per-plate cKDTrees use: land_points forms contiguous coastal blobs
    # rather than uniformly-scattered points, and the default (True/True) construction
    # degrades badly on that kind of clustered data (benchmarked ~8x slower query for a
    # ~45k-point land tree of this shape).
    return cKDTree(land_points, balanced_tree=False, compact_nodes=False)


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
    fluid_density: float = climate.DEFAULT_CLIMATE_DENSITY,
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
    World for the same reason node_density is (see World.climate_density's own comment).
    `fluid_density` (the UI's "Fluid dynamics resolution" Advanced-settings choice) similarly
    doesn't affect generation, only how finely Ocean/Atmospheric Fluid Dynamics mode's own
    grid resolves currents/wind -- see World.fluid_density's own comment for why it's a
    separate knob from climate_density rather than reusing it."""
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
        fluid_density=fluid_density,
    )
    for plate in world.plates:
        _update_plate_omega(plate, world.mantle_centers, damping=None)

    n_continents = sum(1 for p in plates if p.crust_type == "continental")
    world.log_event(f"World generated with {len(plates)} plates ({n_continents} continental).")
    return world


def step_world(world: World, years: float) -> None:
    """Advance the world by `years`.

    Plate movement (skippable via World.simulate_plate_movement) is now two per-plate
    passes rather than the old rotate-then-classify-by-velocity pipeline:

    1. `Plate.shift(world, years)`, for every plate: refit its Euler pole from the mantle
       flow field (damped toward the new target), rotate rigidly by `years` (exact, no
       resampling), and return `D` -- the greatest angular distance any of that plate's own
       nodes actually moved this step. Order doesn't matter here; rotation only touches the
       rotating plate's own frame.
    2. `Plate.deform(world, other_plates, years, D)`, for every plate, in a freshly
       randomized order each turn: reconcile this plate's actual post-shift footprint
       against the footprint it's entitled to occupy -- the sphere minus every *other*
       currently-live plate's own bounding polygon. Territory now overlapping a neighbor is
       collision/subducted (uplift/trench elevation, then the affected line end shrinks, up
       to `D` worth of nodes); territory nobody else claims is a rift (grows, or -- if
       already stretched thin -- spawns a fresh volcano instead); everything else close to a
       neighbor but not overlapping it is transform. Regularizes any line whose spacing has
       drifted, and claims any adjacent unclaimed territory a line's own end-growth can't
       reach (a plate growing toward its own pole, or reclaiming a subducted neighbor's
       vacated ground), all as part of this same call -- see PlateWithLines.deform's own
       docstring. Randomizing the processing order each turn is what keeps two neighbors
       from both claiming the same contested/unclaimed space in the same turn: each plate's
       "what am I entitled to" check runs against whatever state its neighbors are
       *currently* in -- already-deformed neighbors reflect this turn's change, not-yet-
       deformed ones don't -- and average fairness comes from the order changing every turn.

    Then topology changes: fully-subducted plates disappear, colliding continental plates
    merge (at most one per step, only after a sustained 50-100 Myr collision -- see
    merge_split.py), and plates whose flow field no longer fits one rigid rotation well can
    split; any resulting events are logged to world.events for the UI's console. Every step
    also erodes elevation based on the world's current climate (see erosion.py) --
    rain/sheet erosion and weathering, the other half of the weather<->geology coupling from
    climate.py's own terrain-influences-weather mechanics (lapse rate, mountain wind
    deflection, orographic rain shadow) -- relaxes submerged continental crust toward a
    shelf-or-deep-water target based on distance to the nearest land (see bathymetry.py),
    and rolls each active volcano's own eruption chance, growing mineral_deposit_m wherever
    one erupts (see volcanism.py). Right after that, grows/relaxes soil and coal/oil-gas
    deposits from this same step's erosion/flow-routing results (see geology.py).

    Both halves above are individually skippable, live, via World.simulate_plate_movement/
    World.simulate_climate_biomes (see their own docstrings and main.py's /world/controls) --
    "plate movement" means shift/deform, topology changes, and volcanism; "climate & biomes"
    means erosion (which itself computes this step's climate.py fields), hydrology,
    bathymetry, and geology's resource formation. elapsed_years always advances regardless
    of either flag.

    Only ever called while world.fluid_mode == "tectonics_climate" (see main.py's /world/
    step), so a call reaching here is itself "the tectonics & climate mode" acting -- if
    either half above is live, plates and/or climate are about to change, which invalidates
    any wind/ocean-current baseline an earlier Ocean/Atmospheric FD session left behind (see
    World.remembered_wind_u's own docstring), so both are dropped up front rather than left
    to seed a future FD session against terrain they no longer match.
    """
    if world.simulate_plate_movement or world.simulate_climate_biomes:
        world.remembered_wind_u = None
        world.remembered_wind_v = None
        world.remembered_ocean_u = None
        world.remembered_ocean_v = None
        world.remembered_ocean_eta = None
        world.remembered_ocean_temperature_c = None
        world.remembered_ocean_sediment_concentration = None
        world.remembered_ocean_sediment_deposited_m = None
    if world.simulate_plate_movement:
        distances = {plate.plate_id: plate.shift(world, years) for plate in world.plates}
        order = list(world.plates)
        # Deterministic per (seed, elapsed_years) so a replayed session still deforms plates
        # in the same order -- not the same order every turn, which is the whole point (see
        # docstring above), but reproducible given the same seed and step history.
        np.random.default_rng((world.seed, round(world.elapsed_years))).shuffle(order)
        for plate in order:
            others = [p for p in world.plates if p.plate_id != plate.plate_id]
            plate.deform(world, others, years, distances[plate.plate_id])
    world.elapsed_years += years
    if world.simulate_plate_movement:
        for message in merge_split.apply_topology_changes(world, years):
            world.log_event(message)

    erosion_result = None
    if world.simulate_climate_biomes:
        # Gathered once here (rather than independently inside erosion.apply_erosion/
        # bathymetry.apply_bathymetry) since node positions are fixed for the rest of this
        # step -- nothing between here and the next step's shift() moves a node or changes
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

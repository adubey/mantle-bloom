"""Top-level world state: a collection of plates, a mantle-flow field, and elapsed time."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from . import atmosphere_cfd, climate, erosion, geology, hydrology, mantle, merge_split, stranded_basins, volcanism
from .elevation_lines import DEFAULT_NODE_DENSITY
from . import lithosphere_plate
from .lithosphere_plate import generate_plates
from .plates import Plate, gather_node_positions, query_workers

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
    # Count of step_world calls on this world. Drives the cadence of merge_split.py's
    # geometric defragmentation pass (see DEFRAG_INTERVAL_STEPS) -- a step counter, not a
    # year counter, since that pass is about accumulated topology drift, not elapsed time,
    # and step sizes vary. A plain-int default is a class attribute, so worlds pickled
    # before this field existed still load (reading 0) -- see persistence.py.
    steps_taken: int = 0
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
    # resolution" Advanced-settings choice, climate.FLUID_DENSITY_CHOICES -- capped lower than
    # climate_density's own choices, see that constant's own comment for why) -- independent
    # of climate_density, read only by atmosphere_cfd.init_atmosphere_cfd and world.py's own
    # _advance_fluid_dynamics to size the atmospheric wind-solver grid. Unlike climate_density,
    # lowering this doesn't touch the Biome/Combined/Resources/Soil-Quality render grid or
    # erosion.py's own climate sampling at all -- it only trades off how finely the wind solve
    # resolves the flow against how long each tectonics step's own fixed-real-time CFD
    # advancement takes, since CFL substep count scales with grid spacing (see
    # fluid_dynamics.cfl_substeps) and per-substep cost scales with cell count -- a real
    # per-step cost, since the wind solve runs continuously alongside tectonics.
    fluid_density: float = climate.DEFAULT_FLUID_DENSITY
    # Sustained-collision tracking for merge_split.py: (plate_id, plate_id) -> accumulated
    # convergent years. See merge_split.update_collision_progress.
    collision_progress: dict[tuple[int, int], float] = field(default_factory=dict)
    # The territory-overlap sibling of collision_progress: (plate_id, plate_id) -> accumulated
    # years the two continental plates have sat deeply superimposed. Feeds the forced
    # continental merge that resolves a stuck multi-plate pile-up the closing-rate timer above
    # can't (they overlap too completely to register a closing rate). See
    # merge_split.update_overlap_progress. A default_factory field -> backfilled on load of an
    # older save (persistence._backfill_added_fields).
    overlap_progress: dict[tuple[int, int], float] = field(default_factory=dict)
    # Cross-step memory for the stranded-basin diagnostic (docs/debugging.md): one
    # `stranded_basins.StrandedBasinTrack` per endorheic, below-sea-level, ocean-disconnected
    # basin currently present, reconciled by centroid proximity every hydrology step (see
    # stranded_basins.reconcile_world_tracks) so the report can say how long each pit has
    # persisted. Same "lightweight per-key first-seen tracker" role `collision_progress` plays
    # for plate pairs -- diagnostic only, nothing in the physics reads it back. A
    # `default_factory` field, so an older save without it is backfilled on load (see
    # persistence._backfill_added_fields).
    stranded_basin_tracks: list = field(default_factory=list)
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
    # Last step's erosion breakdown (see erosion.ErosionResult), retained here purely so the
    # Geomorph Rate debug view (render_image._render_geomorph_view) can colour every node by
    # its net elevation change this step -- geology.py still receives its own copy as a direct
    # step_world argument (that's its only same-turn consumer; see ErosionResult's docstring).
    # None until the first climate/erosion step runs, and left at last-good (never reset to
    # None) if simulate_climate_biomes is later toggled off, same tolerance as the two caches
    # above. Not persisted -- a freshly loaded save shows the neutral field until it's stepped.
    erosion_cache: erosion.ErosionResult | None = None
    # Nearest-land spatial index backing distance_from_land_approx (below) -- reset to None
    # once per step (step_world, right where node_cloud is gathered, since land nodes' own
    # world positions and elevations can both have changed since the last build) and rebuilt
    # lazily on first use after that. None means "not built yet this step," not "no land" --
    # distance_from_land_approx itself is what tells those two cases apart.
    land_kdtree_cache: cKDTree | None = None
    # The render path's full node-cloud k-d tree (a cKDTree over plates.collect_all_points'
    # concatenated node positions -- ~131 K at node_density 4), paired with the concatenated
    # (points, elevation, owner) arrays it indexes into -- see
    # render_image._node_cloud_and_tree. The tree build is ~20 ms at that size (docs/profiling.md
    # #6 -- the query over the render grid is the larger cost and is already workers=parallel);
    # the node cloud is fixed between steps (only elevation/other per-node fields still move
    # mid-step, and the render path never runs mid-step), so this is built once by the first
    # render after a step and reused by every subsequent render -- and by the several separate
    # resamples within a single combined/elevation render -- until step_world resets it.
    # Persisted like the other caches but dropped on load (persistence._drop_derived_caches).
    node_kdtree_cache: tuple[np.ndarray, np.ndarray, np.ndarray, cKDTree] | None = None
    # Live-adjustable via POST /world/controls (see main.py) for the UI's "Controls" window
    # -- unlike axial_tilt_deg/node_density (fixed at generation), these are meant to be
    # tweaked mid-simulation. sea_level_m replaces the bare `elevation <= 0.0` convention
    # every is_ocean check in this codebase used to hardcode (climate.py, hydrology.py);
    # solar_multiplier scales climate.py's own SUNLIGHT constant. Changing
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
    # step_world skips climate/erosion/hydrology/resource-formation entirely -- by
    # far the most expensive part of a step (climate.py's grid computation and hydrology.py's
    # flow routing) -- leaving climate_cache/hydrology_cache at whatever they were last
    # computed to (stale, the same one-step-behind tolerance World.climate_cache already
    # documents) rather than None, so a render/stats call right after toggling this off still
    # shows the last real climate snapshot instead of going blank.
    simulate_climate_biomes: bool = True
    # Which wind field feeds climate.py (live-adjustable via POST /world/controls, the
    # "Controls" window). "cfd": the genuine time-integrated shallow-water solve in
    # atmosphere_cfd.py, advanced once per step by _advance_fluid_dynamics -- the most
    # expensive single piece of a step at the default fluid_density (see
    # docs/simulation-model.md#fd-performance). "diagnostic" (default): skip that solve entirely
    # and let climate.compute_climate rebuild wind/air-temperature from its own closed-form
    # ABL-style formulas (compute_wind + compute_air_temperature_diagnostic) every call, the same
    # way it already does during the pre-CFD cold-start bootstrap. Reproduces ~85-90% of the land
    # biome map and precipitation within ~10% for a fraction of the cost -- see
    # docs/simulation-model.md#wind-model and docs/TODO.md -- so it's the default; opt into "cfd"
    # via Controls for the full solve. atmosphere_cfd_state is still kept in sync (init'd at
    # generation, never cleared) so switching to "cfd" mid-session resumes from a real, if
    # now-stale, state rather than a cold start.
    wind_model: str = "diagnostic"
    # Atmospheric wind-solver state -- see docs/simulation-model.md#atmospheric-fluid-dynamics.
    # Always on, not a mode: generate_world populates it immediately after constructing this
    # World (atmosphere_cfd.init_atmosphere_cfd) and it's never re-initialized or cleared again
    # for the rest of this world's life -- step_world's own _advance_fluid_dynamics just keeps
    # advancing it by a fixed SECONDS_PER_TECTONIC_STEP (one simulated day) every tectonics
    # step, gated on simulate_climate_biomes the same way erosion/hydrology already are.
    # climate.py's compute_climate sources wind_u/wind_v (and air_temperature_c) from this
    # state's continuously-evolving fields rather than its own diagnostic compute_wind, which
    # now exists only as the one-time cold-start bootstrap init_atmosphere_cfd falls back to
    # before this field is ever populated. Ocean currents and precipitation are *not* CFD-
    # solved -- they're diagnostic in climate.py every step (the ocean shallow-water solver was
    # retired for lack of a stable operating point; see climate.py's module docstring). Typed
    # `| None` only because a dataclass field default can't construct it inline -- never
    # actually `None` once generate_world has returned.
    atmosphere_cfd_state: atmosphere_cfd.AtmosphereCFDState | None = None

    def log_event(self, message: str) -> None:
        self.events.append((self.elapsed_years, message))
        if len(self.events) > MAX_EVENT_LOG_LENGTH:
            del self.events[: len(self.events) - MAX_EVENT_LOG_LENGTH]

    def distance_from_land_approx(self, points: np.ndarray) -> np.ndarray:
        """Approximate distance from each given world-xyz point (shape (n, 3)) to the
        nearest land node (elevation > sea_level_m) anywhere in this world -- lazily builds
        and reuses land_kdtree_cache (see its own docstring) off every plate's own public
        Plate.map_world_points()/ElevationPoint.get_elevation() interface. geology.py uses
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


def generate_world(
    seed: int,
    num_plates: int | None = None,
    continental_fraction: float | None = None,
    land_fraction: float | None = None,
    num_mantle_centers: int = DEFAULT_MANTLE_CENTERS,
    axial_tilt_deg: float | None = None,
    node_density: float = 1.0,
    initial_soil_maturity: float | None = None,
    climate_density: float = climate.DEFAULT_CLIMATE_DENSITY,
    fluid_density: float = 1.0,
    extra_sites_per_plate: int = lithosphere_plate.EXTRA_SITES_PER_PLATE,
) -> World:
    """`num_plates` is optional -- see lithosphere_plate.generate_plates for why: the world
    tiles itself into a plausible number of plates rather than requiring the caller to pick
    one. `continental_fraction`/`land_fraction` are the UI's generation sliders -- see
    lithosphere_plate.generate_plates. `axial_tilt_deg` is the UI's third generation slider
    (degrees, defaults to DEFAULT_AXIAL_TILT_DEG = Earth's real tilt) -- doesn't affect plate
    generation at all, only climate.py's insolation, read at render time long after
    generation, which is why it's stored on World rather than consumed once here.
    `node_density` (the UI's "point density" choice) similarly affects only generation
    itself directly, but -- unlike axial_tilt_deg -- is stored on World because every later
    step also needs it (see World.node_density's own comment). `initial_soil_maturity` (the
    UI's "initial soil maturity" slider, 0 to 1, default None -> no seeding) is a one-time
    generation-time seed -- like continental_fraction/land_fraction, not stored on World,
    since nothing later needs to know what it was (see geology.seed_initial_soil).
    `climate_density` (the UI's "climate & biome resolution" choice, see
    climate.CLIMATE_DENSITY_CHOICES) doesn't affect plate generation at all, only how finely
    climate.py's own grid resolves the world's climate every future step/render -- stored on
    World for the same reason node_density is (see World.climate_density's own comment).
    `fluid_density` (the UI's "Fluid dynamics resolution" Advanced-settings choice) similarly
    doesn't affect the plates/mantle generated above, but *does* immediately seed this world's
    permanent atmosphere_cfd_state (see below) at its own resolution -- see World.fluid_density's
    own comment for why it's a separate knob from climate_density rather than reusing it.
    `extra_sites_per_plate` controls how many adjacent Voronoi cells each plate fuses at
    generation (see lithosphere_plate.build_plate_tiling) -- higher means lumpier, less
    convex initial plate outlines; 0 is the old one-cell-per-plate tiling."""
    plates = generate_plates(
        seed, num_plates, continental_fraction, land_fraction, node_density, extra_sites_per_plate=extra_sites_per_plate
    )
    rng = np.random.default_rng(seed)
    mantle_centers = mantle.generate_convection_centers(rng, n_centers=num_mantle_centers)

    world = World(
        seed=seed,
        plates=plates,
        mantle_centers=mantle_centers,
        next_plate_id=len(plates),
        axial_tilt_deg=axial_tilt_deg if axial_tilt_deg is not None else DEFAULT_AXIAL_TILT_DEG,
        node_density=node_density,
        climate_density=climate_density,
        fluid_density=fluid_density,
    )
    if initial_soil_maturity is not None:
        geology.seed_initial_soil(world.plates, seed, initial_soil_maturity)

    # See World.atmosphere_cfd_state's own comment for why it's populated here,
    # unconditionally, rather than lazily. `terrain` is the diagnostic bootstrap snapshot the
    # wind solver is seeded from (compute_wind's own latitude-banded field), before the state
    # it will read from every step after exists.
    height, width = climate.grid_dimensions(world.climate_density)
    terrain = climate.compute_climate(world, height, width)
    world.atmosphere_cfd_state = atmosphere_cfd.init_atmosphere_cfd(world, terrain)

    n_continents = sum(1 for p in plates if p.crust_type == "continental")
    world.log_event(f"World generated with {len(plates)} plates ({n_continents} continental).")
    return world


def _advance_fluid_dynamics(world: World, node_cloud: tuple[np.ndarray, list[Plate]]) -> None:
    """Advances World.atmosphere_cfd_state by its fixed SECONDS_PER_TECTONIC_STEP (one
    simulated day) once per tectonics step, *before* erosion/hydrology each step (unlike a
    naive post-erosion ordering) -- so erosion/hydrology read post-substep, not pre-substep,
    wind. No-op when `world.wind_model != "cfd"` (see World.wind_model): the diagnostic wind
    model rebuilds wind from climate.py's own formulas every call and never reads the CFD
    state, so advancing it would be wasted work -- the single biggest saving of that mode.
    `world.climate_cache` still holds last step's snapshot at this point, so a fresh
    climate snapshot is always computed here regardless of whether fluid_density matches
    climate_density; erosion.apply_erosion computes its own snapshot right after this returns
    (on the same, still-unchanged post-tectonics world.plates), so this is one extra
    compute_climate call per step in exchange for correct wind forcing -- passed
    skip_moisture=True since refresh_forcing consumes only elevation/is_ocean/the temperature
    baseline, so this call doesn't pay for the humidity/precipitation sweep."""
    if world.wind_model != "cfd":
        return
    terrain = climate.compute_climate(
        world, *climate.grid_dimensions(world.fluid_density), node_cloud=node_cloud, skip_moisture=True
    )
    atmosphere_cfd.refresh_forcing(world, world.atmosphere_cfd_state, terrain)
    atmosphere_cfd.step_atmosphere_cfd(world, world.atmosphere_cfd_state, atmosphere_cfd.SECONDS_PER_TECTONIC_STEP)


def step_world(world: World, years: float) -> None:
    """Advance the world by `years`.

    Plate movement (skippable via World.simulate_plate_movement) is two per-plate passes:
    `Plate.shift(world, years)` for every plate (refit Euler pole from torque balance, rotate
    rigidly), then `Plate.deform(world, other_plates, years, D)` for every plate in a freshly
    randomized order each turn (Mohr-Coulomb yield/isostasy -- see lithosphere_plate.py).
    Randomizing the processing order each turn is what keeps two neighbors from both claiming
    the same contested/unclaimed space in the same turn.

    Then topology changes: fully-subducted plates disappear, colliding continental plates
    merge, and plates whose flow field no longer fits one rigid rotation well can split; any
    resulting events are logged to world.events for the UI's console (see merge_split.py).

    The atmospheric wind solve (see _advance_fluid_dynamics) advances *before* erosion/
    hydrology each step, so erosion/hydrology read post-substep wind. Every step also
    erodes elevation based on the world's current climate (see erosion.py), and rolls each
    active volcano's own eruption chance (see volcanism.py). Right after that, grows/relaxes
    soil and coal/oil-gas deposits from this same step's erosion/flow-routing results (see
    geology.py). Isostasy (lithosphere_plate.py/lithosphere.py) supersedes what a separate
    bathymetry-relaxation pass used to do for submerged continental crust.

    Both plate movement and climate/biomes are individually skippable, live, via
    World.simulate_plate_movement/World.simulate_climate_biomes (see their own docstrings and
    main.py's /world/controls) -- elapsed_years always advances regardless of either flag.
    """
    world.steps_taken += 1
    # The render path's cached node-cloud k-d tree (see World.node_kdtree_cache) is a pure
    # function of node positions, which shift()/deform()/topology changes below are about to
    # move -- drop it now so the first render after this step rebuilds it. (land_kdtree_cache
    # is reset separately, inside the simulate_climate_biomes block, since only that path
    # reads it.)
    world.node_kdtree_cache = None
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
        # Stamp/clear ElevationLine.overlap_onset_years and advance World.overlap_progress
        # against this step's final geometry (see merge_split.update_overlap_tracking /
        # docs/debugging.md).
        merge_split.update_overlap_tracking(world, years)

    erosion_result = None
    if world.simulate_climate_biomes:
        # Gathered once here since node positions are fixed for the rest of this step --
        # nothing between here and the next step's shift() moves a node or changes line
        # topology (only elevation and other per-node fields still change, which each of
        # climate.py/erosion.py/hydrology.py still reads fresh off world.plates itself) --
        # see plates.gather_node_positions's own docstring for why this was worth factoring
        # out. Skipped entirely, alongside erosion below, when simulate_climate_biomes is off.
        node_cloud = gather_node_positions(world.plates)
        # Same "fixed for the rest of this step" reasoning as node_cloud above applies to
        # land_kdtree_cache -- reset here, then lazily rebuilt on whichever future caller
        # this step reads it first.
        world.land_kdtree_cache = None
        _advance_fluid_dynamics(world, node_cloud)
        erosion_result = erosion.apply_erosion(world, years, node_cloud=node_cloud)
        world.erosion_cache = erosion_result
    if world.simulate_plate_movement:
        volcanism.apply_volcanic_activity(world, years)
    if erosion_result is not None:
        geology.apply_resource_formation(world, years, erosion_result)
        # Reconcile the stranded-basin tracker against this step's freshly-rebuilt depression
        # hierarchy (world.hydrology_cache, just set by erosion) -- only on a step that
        # actually recomputed hydrology, so persistence timers count simulated hydrology steps.
        stranded_basins.reconcile_world_tracks(world)



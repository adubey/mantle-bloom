"""Top-level world state: a collection of plates, a mantle-flow field, and elapsed time."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from . import atmosphere_cfd, climate, erosion, eustasy, faults, gaps, geology, healpix_grid, hydrology, magma_transport, mantle, merge_split, stranded_basins, volcanism, worldsketch
from .elevation_lines import DEFAULT_NODE_DENSITY
from . import lithosphere_plate
from .lithosphere_plate import generate_plates
from .plates import Plate, gather_node_positions, query_workers

DEFAULT_MANTLE_CENTERS = 8
DEFAULT_AXIAL_TILT_DEG = 23.5
# Bounds how large World.events can grow over a long play session -- the UI's console only
# ever needs recent history, not an unbounded transcript.
MAX_EVENT_LOG_LENGTH = 200
# Bounds World.removed_points_log the same way -- a count cap, not an age cap, since
# subduction/defragmentation/merge remove nodes essentially every step on a full-size save
# and an age-based retention window would grow unboundedly at high node_density. Sized well
# past MAX_EVENT_LOG_LENGTH since removal is per-node, not per-event -- see
# World.record_removed_points and the "Added/Removed Points" (`nodeAge`) debug render view.
MAX_REMOVED_POINTS_LOG = 20_000
# Bounds World.corner_notch_log -- one entry per plate per step while World.debug_diagnostics
# is on, so a longer debugging session can still build up a lot of entries. More generous than
# MAX_EVENT_LOG_LENGTH since this is meant for an active debugging session (see
# World.log_corner_notch), not indefinite retention.
MAX_CORNER_NOTCH_LOG_LENGTH = 2_000

# The dimensionless geomorphic-budget tuning knobs on World (see the field group below),
# named once here so main.py's /world/controls route can validate/apply/echo them without
# repeating the list. All default 1.0; all must stay >= 0.
TUNING_MULTIPLIER_FIELDS = (
    "rain_erosion_multiplier",
    "river_erosion_multiplier",
    "wind_erosion_multiplier",
    "ocean_erosion_multiplier",
    "coastal_leveling_multiplier",
    "glacier_erosion_multiplier",
    "seismic_erosion_multiplier",
    "river_deposition_multiplier",
    "ocean_deposition_multiplier",
    "collision_uplift_multiplier",
    "collision_uplift_reach_multiplier",
    "volcanism_multiplier",
    "fault_relief_multiplier",
)


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
    # Another fixed per-world property: which Premade-worlds tab entry (if any) this world
    # was generated from ("earth"/"pangaea"/"got"), `None` for every other Generate World tab.
    # Stored (not just consumed once at generation) so a later step, long after generation,
    # can still tell -- hydrology.py's own Gibraltar note is the one thing that currently
    # reads this, to apply an "earth"-only override; a plain default means a world pickled
    # before this field existed still loads (reading None, no override) same as any other
    # field added here, see persistence.py.
    premade_world_id: str | None = None
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
    # Cross-step memory for the gap-age diagnostic (docs/debugging.md's overlapAge section):
    # one `gaps.GapTrack` per currently-uncovered lattice cluster, reconciled by centroid
    # proximity at the same cadence as gaps.fill_gaps_by_growing_neighbours (see
    # gaps.reconcile_gap_tracks) -- the same "lightweight per-key first-seen tracker" role
    # stranded_basin_tracks plays for
    # basins, since a gap cluster has no persistent identity across steps any more than a
    # basin does. Diagnostic only, nothing in the physics reads it back. A `default_factory`
    # field -> backfilled on load (see persistence._backfill_added_fields).
    gap_tracks: list = field(default_factory=list)
    # Cross-step memory for lateral magma transport (GitHub issue #205, magma_transport.py):
    # exported convergent-boundary melt banked here every step (by LithospherePlate.deform())
    # and drained -- fully or partially, see MagmaParcel.unplaced_cycles -- every
    # magma_transport.MAGMA_TRANSPORT_INTERVAL_STEPS'th step's magma_transport.
    # run_magma_transport call below. A `default_factory` field -> backfilled on load (see
    # persistence._backfill_added_fields).
    pending_magma_parcels: list = field(default_factory=list)
    # Years accumulated since the last run_magma_transport firing -- that pass's own per-
    # destination-node rate cap needs the *banked* elapsed time (parcels queue across
    # MAGMA_TRANSPORT_INTERVAL_STEPS steps), not just the current step's own `years`. Reset to
    # 0.0 after each firing. A plain-float default, so an old pickle falls through to 0.0 with
    # no persistence backfill needed (same as steps_taken).
    magma_transport_banked_years: float = 0.0
    # Cross-step memory for the "Added/Removed Points" debug view's removed-node half (the
    # added half needs no cross-step state -- it reads straight off each live node's own
    # ElevationLine.node_created_years). A node vanishes from every plate's own node cloud the
    # instant it's removed (retreat, interior-subduction carve, merge absorption, defrag
    # stripping, whole-plate subduction -- see World.record_removed_points' own call sites),
    # so unlike every other per-node diagnostic field this can't ride along on the node itself
    # -- it has to be its own short-lived history buffer here. Each entry is
    # (world_xyz, removed_years, plate_id); capped by count, not age, see
    # MAX_REMOVED_POINTS_LOG's own comment. `default_factory` -> backfilled on load (see
    # persistence._backfill_added_fields).
    removed_points_log: list[tuple[np.ndarray, float, int]] = field(default_factory=list)
    # Gate for the verbose, structured `_fill_corner_notch_frontier` decision log below -- off by
    # default (a plain-scalar field, so an old pickle falls through to False with no
    # persistence backfill needed), on by default for a "Debugging Worlds" tab world, and
    # toggleable for any other loaded save via POST /world/controls. Kept as an explicit flag
    # rather than always logging: this runs once per plate per step, and the codebase's own
    # precedent (lakes.summarize_lake_events, see docs/debugging.md) is that this volume of
    # per-step diagnostic detail does not belong in the always-on Event Console -- see
    # corner_notch_log's own comment for where it goes instead.
    debug_diagnostics: bool = False
    # Verbose, structured decision log for `LithospherePlate._fill_corner_notch_frontier` (see
    # World.log_corner_notch) -- populated only while `debug_diagnostics` is True. Deliberately
    # separate from `events` (the always-on Event Console): this can fire once per plate per
    # step, far higher volume than that log is meant to carry, so it gets its own capped buffer
    # and its own `GET /world/corner_notch_log` endpoint / panel instead of ever going through
    # `log_event`. Capped by count like `events`, just a more generous ceiling since it's
    # meant for an active debugging session, not indefinite retention. A `default_factory`
    # field -> backfilled on load (see persistence._backfill_added_fields).
    corner_notch_log: list[dict] = field(default_factory=list)
    # Debug-world-only: plate_id -> a fixed world-frame angular velocity (rad/s) that
    # LithospherePlate.shift uses verbatim every step, bypassing torque.shift_plate's own
    # torque-balance recompute entirely for that plate (see torque.py) -- the "Debugging
    # Worlds" tab's scripted-motion mechanism (see debug_worlds.py), so a tiny hand-built
    # scenario moves exactly as scripted every step rather than however real ridge-push/
    # slab-pull/basal-drag against world.mantle_centers happens to settle it. Empty for every
    # ordinarily-generated world. A `default_factory` field -> backfilled on load (see
    # persistence._backfill_added_fields).
    pinned_omegas: dict[int, np.ndarray] = field(default_factory=dict)
    # Intraplate fault lines (see faults.py) -- a `faults.Fault` per trace, geometry stored
    # in the owning plate's local frame so it rotates with the crust. Grown/retired/applied
    # by faults.update_faults every step and re-homed across topology changes by
    # faults.reconcile_faults. A `default_factory` field -> backfilled on load of an older
    # save (persistence._backfill_added_fields). `next_fault_id` is the monotonic id source,
    # a plain-int default so an old pickle falls through to 0.
    faults: list = field(default_factory=list)
    next_fault_id: int = 0
    # Boundary faults (see faults.generate_boundary_faults) -- `faults.Fault` traces laid
    # along every live plate-boundary segment, classified by that segment's motion
    # (convergent -> reverse, divergent -> normal, transform -> strike-slip). Rebuilt from
    # scratch every step against the current geometry, so unlike `faults` they are never
    # aged, retired, culled, or re-homed. They feed the same relief / earthquake / rendering
    # / fault_influence path as intraplate faults. `default_factory` -> backfilled on load.
    boundary_faults: list = field(default_factory=list)
    # Fault systems (see faults.FaultSystem) -- a first-class zone above the individual
    # trace: one long curving master lineament plus the strand family (`Fault`s carrying its
    # `system_id`) scattered along its belt. Same local-frame storage, same
    # update_faults / reconcile_faults lifecycle. `default_factory` -> backfilled on load;
    # `next_fault_system_id` is a plain-int default so an old pickle falls through to 0.
    fault_systems: list = field(default_factory=list)
    next_fault_system_id: int = 0
    # Recent earthquakes (see faults.Earthquake) -- transient located events produced when an
    # active fault ruptures (faults.update_faults), each with a true-world-frame epicenter, a
    # magnitude, and the id of the fault that slipped. Not persistent geology: pruned once
    # older than faults.EARTHQUAKE_RETAIN_MYR every step, and never re-homed across topology
    # changes (an epicenter is a fixed point in space, and the window is short). Read by
    # erosion.py (a local seismic-erosion burst around each epicenter) and exposed via
    # GET /world/earthquakes for the "Fault lines" view's fading epicenter overlay.
    # `default_factory` -> backfilled on load; `next_earthquake_id` is a plain-int default an
    # old pickle falls through to 0.
    earthquakes: list = field(default_factory=list)
    next_earthquake_id: int = 0
    # Which model applies plate-boundary deformation (uplift / rift / transform). Live-
    # adjustable via POST /world/controls, same pattern as `wind_model`:
    #   "fault" (default) -- that boundary thickening is gated to proximity to an active fault
    #     trace (faults.fault_influence), and faults.py's own relief layer is scaled up to
    #     carry the deformation the bands give up, so transformation localises onto fault
    #     lines. Faults spawn boundary-hugging (see faults.SPAWN_PLACE_*), so the collision
    #     zone still deforms -- just as a family of fault-tracking ridges rather than one
    #     smooth swell.
    #   "boundary" -- LithospherePlate.deform's smooth distance-band thickening at the polygon
    #     edge, exactly as before the faults rework (bit-identical to a pre-field pickle
    #     stepped in "boundary" mode).
    #   "both" -- boundary bands at full strength *and* the scaled-up fault relief layer.
    # See faults.FAULT_DEFORMATION_MODES and LithospherePlate.deform.
    fault_deformation_mode: str = "fault"
    # Which structure `render_image._node_cloud_and_tree` resamples the node cloud through --
    # live-adjustable via POST /world/controls, but backend/API-only for now (no Controls-panel
    # entry): this is issue #133's phase-1 proving-out flag, not yet a user-facing tuning knob.
    # See healpix_grid.NODE_CLOUD_RESAMPLE_MODE_CHOICES:
    #   "kdtree" (default) -- today's `cKDTree(all_points).query(...)`, unchanged.
    #   "healpix" -- scatter the node cloud onto a `HealpixGrid` sized to the node count
    #     (healpix_grid.nside_for_node_count), wavefront-fill the empty pixels, and resolve
    #     every render query through an `ang2pix` lookup instead of a tree query
    #     (healpix_grid.build_node_pixel_index / NodePixelIndex). `_classify_terrain_relief`'s
    #     radius search is excluded (see its own call site in render_image.py) and still uses a
    #     real `cKDTree` regardless of this flag -- out of scope for issue #133's phase 1.
    node_cloud_resample_mode: str = "kdtree"
    # Human-readable log for the UI's event console, each entry (elapsed_years, message).
    events: list[tuple[float, str]] = field(default_factory=list)
    # One stats.compute_stats(self) snapshot per real advance (generate_world, then every
    # step_world -- see record_stats), for the Stats panel's history charts. Previously kept
    # only in the frontend's own React state (built one fetch at a time from the stateless
    # GET /world/stats -- see stats.py's own module docstring), which meant a Save/Load
    # round-trip silently dropped the whole run's history even though every other per-step
    # record (events, corner_notch_log) already survives one. Recorded here instead so it
    # rides along in the pickle like everything else on World, and GET /world/stats_history
    # (main.py) hands the *full* series back after a load. Deliberately uncapped, unlike
    # events/corner_notch_log -- a history chart needs the whole run, not just recent
    # activity, and mirrors what the frontend already kept unbounded client-side before this.
    # `default_factory` field -> backfilled on load (see persistence._backfill_added_fields).
    stats_history: list[dict] = field(default_factory=list)
    # This step's climate snapshot (see climate.py), populated by erosion.py -- which needs
    # a fresh one every step regardless -- and reused by /world/stats and a climate map
    # render so they don't each trigger their own (~50ms) recomputation the same turn. See
    # climate.compute_climate_cached and climate.py's own module docstring for why reusing
    # a value that's up to one step stale is an accepted simplification here, not a bug.
    climate_cache: climate.ClimateFields | None = None
    # This step's flow-routing snapshot (see hydrology.py), populated by erosion.py
    # alongside climate_cache -- same reuse pattern, same one-step-stale simplification.
    hydrology_cache: hydrology.HydrologyFields | None = None
    # steps_taken at the moment hydrology_cache was last (re)populated by erosion.py, or None
    # while hydrology_cache itself is None. Lets a reader (stats.compute_stats) tell "one step
    # stale" (the tolerated case above -- this step just hasn't run erosion *yet*) apart from
    # "frozen for N steps" (simulate_climate_biomes toggled off, so hydrology_cache.is_ocean
    # keeps answering for a coastline tectonics has long since moved on from -- see GitHub
    # issue #121): compare against World.steps_taken at read time.
    hydrology_cache_step: int | None = None
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
    # (points, elevation, owner) arrays it indexes into -- see render_image._node_cloud_and_tree.
    # The tree build is ~20 ms at that size (docs/profiling.md #6 -- the query over the render
    # grid is the larger cost and is already workers=parallel); the node cloud is fixed between
    # steps (only elevation/other per-node fields still move mid-step, and the render path
    # never runs mid-step), so this is built once by the first render after a step and reused
    # by every subsequent render -- and by the several separate resamples within a single
    # combined/elevation render -- until step_world resets it. Elevation is deliberately NOT
    # shared any earlier than this (e.g. with climate.py's own per-step resample, which runs
    # before erosion has finished mutating it that same step) -- see
    # node_position_tree_cache's own docstring for the piece that *is* safe to share mid-step.
    # Persisted like the other caches but dropped on load (persistence._drop_derived_caches).
    # The 4th slot is a `cKDTree` under the default "kdtree" node_cloud_resample_mode, or a
    # `healpix_grid.NodePixelIndex` under "healpix" -- both expose the same `.query()` shape
    # every caller here actually uses (see render_image._node_cloud_and_tree).
    node_kdtree_cache: tuple[np.ndarray, np.ndarray, np.ndarray, cKDTree | healpix_grid.NodePixelIndex] | None = None
    # The *positions-only* half of the cache above -- (all_points, cKDTree over them), no
    # elevation/owner -- shared between climate.py's own per-step resample
    # (`_sample_elevation_and_crust`, which runs early in a step, before erosion has finished
    # mutating elevation this same step) and render_image._node_cloud_and_tree (which runs
    # only after a step fully settles). Node *positions* really are fixed for the rest of a
    # step the instant shift()/deform()/topology changes finish (nothing moves a node again
    # until the next step's shift()), so this narrower cache is safe to write and read at any
    # point in that window -- unlike node_kdtree_cache's own elevation, which erosion.py is
    # still about to change when climate.py builds this. Reset alongside node_kdtree_cache
    # (same invalidation event: a node moved). Persisted like the other caches but dropped on
    # load (persistence._drop_derived_caches).
    node_position_tree_cache: tuple[np.ndarray, cKDTree] | None = None
    # "healpix" node_cloud_resample_mode's own caches, same invalidation event (a node moved)
    # as node_kdtree_cache/node_position_tree_cache above, reset alongside them. Split into two
    # because they vary on different things: the `HealpixGrid` itself (keyed by nside, which
    # only depends on node *count*) is stable for a world's whole life and only ever rebuilt if
    # that count actually changes class, while the node-pixel assignment (`NodePixelIndex`, the
    # scatter+fill result) depends on node *positions* and is rebuilt every step like the
    # cKDTree caches above. Built/reused via plates.cached_node_healpix_index -- the "healpix"
    # counterpart to node_position_tree_cache/cached_node_position_tree, shared the same way
    # across render_image._node_cloud_and_tree and climate._sample_elevation_and_crust (issue
    # #133 phase 2). Persisted like the other caches but dropped on load
    # (persistence._drop_derived_caches).
    node_healpix_grid_cache: tuple[int, healpix_grid.HealpixGrid] | None = None
    node_healpix_index_cache: healpix_grid.NodePixelIndex | None = None
    # `_classify_terrain_relief`'s own radius-search tree (see its call site in
    # render_image._render_grid_arrays) -- only ever built under "healpix" node_cloud_resample_
    # mode, since the default "kdtree" mode already has a real cKDTree in node_kdtree_cache's
    # 4th slot. `query_ball_point` has no HEALPix-pixel-index equivalent (issue #133 excludes
    # this call from scope), so this is a real `cKDTree`, built lazily only when the Elevation
    # view's relief toggles are on. Reset alongside the caches above.
    node_kdtree_relief_cache: cKDTree | None = None
    # Per-node hillshade brightness multiplier (see render_image._hillshade_for_world), a pure
    # function of node positions/elevation like the caches above -- computed once (a local
    # plane-fit gradient over each node's own k nearest neighbours, see
    # render_image._compute_hillshade) and reused by every render this step (Elevation and
    # Combined views both resample it the same way they already resample channel_depth/
    # lake_depth). Reset alongside node_kdtree_cache -- same invalidation event (a node moved
    # or its elevation changed), since it's derived directly from that cache's own tuple.
    node_hillshade_cache: np.ndarray | None = None
    # Live-adjustable via POST /world/controls (see main.py) for the UI's "Controls" window
    # -- unlike axial_tilt_deg/node_density (fixed at generation), these are meant to be
    # tweaked mid-simulation. sea_level_m replaces the bare `elevation <= 0.0` convention
    # every is_ocean check in this codebase used to hardcode (climate.py, hydrology.py);
    # solar_multiplier scales climate.py's own SUNLIGHT constant. Changing
    # either forces an immediate climate_cache recompute (see main.py's controls route) so
    # /world/render and /world/stats reflect it right away, without waiting for a step.
    sea_level_m: float = 0.0
    # Eustatic sea level (see eustasy.py). `sea_level_m` above is no longer a fixed input --
    # `step_world` re-solves it every step so it tracks the ocean volume this conserved
    # water-column budget represents against the world's changing hypsometry (deeper basins /
    # drowned continents -> lower stand). `None` until first initialized (a freshly built
    # World, or a save written before eustasy existed); `eustasy.initialize_water_budget`
    # snapshots it from the flat starting sea level at generation. The `/world/controls`
    # slider sets this budget rather than `sea_level_m` directly (adds/removes ocean water).
    ocean_water_column_m: float | None = None
    solar_multiplier: float = 1.0
    # The "Ice Age Frequency" Controls slider (Climate tab): the full period, in years, of a
    # slow glacial<->interglacial temperature oscillation driven purely by `elapsed_years`.
    # 0.0 (the default) disables it entirely -- no ice-age cycle, bit-identical to before this
    # field existed and to an older save (a plain-scalar default an old pickle falls through
    # to, no persistence backfill needed). See climate.ice_age_glacial_intensity /
    # climate.ice_age_cooling_c: a raised-cosine intensity in [0, 1] (0 at generation and at
    # every whole multiple of the period, 1 at the half-period glacial maximum) scales a
    # uniform cooling applied to climate.py's land+ocean temperature baselines, so every
    # downstream consumer (air temperature, biomes, hydrology's freeze/glacier growth,
    # erosion) feels the ice age through fields it already reads. During a deep enough glacial
    # phase hydrology also lets ice caps spread out over freezing polar ocean (see
    # hydrology.compute_hydrology's `polar_ocean_ice` mask), and eustasy.py debits the water
    # locked into that ice (and into glaciers/large lakes) from the ocean budget, so sea level
    # falls as the ice sheets grow.
    ice_age_period_years: float = 0.0
    # Geomorphic-budget tuning knobs -- live-adjustable via POST /world/controls, same
    # "Controls" window / same immediate-climate-recompute pattern as sea_level_m/
    # solar_multiplier. Every one is a dimensionless multiplier, 1.0 == the model's
    # untuned behaviour exactly (so an old save, or a user who never opens the panel, is
    # bit-identical to before these existed). They exist because mantle-bloom's long-run
    # land budget is a balance between a handful of erosion, deposition, uplift and
    # volcanism processes whose individual hard-coded rates can't be re-tuned for one seed
    # without regressing another -- these let the user rebalance a *live* world and watch
    # the result. Read directly off `world` by erosion.apply_erosion (the *_erosion_/
    # *_deposition_ knobs), lithosphere_plate.LithospherePlate.deform (the collision_uplift_*
    # knobs, via rheology.apply_convergent_deformation's `strength`) and
    # volcanism.apply_volcanic_activity (volcanism_multiplier) -- none of those need a
    # signature change since they already take `world`. Plain-scalar defaults, so an older
    # pickle falls through to 1.0 with no persistence backfill needed (see persistence.py).
    rain_erosion_multiplier: float = 1.0
    river_erosion_multiplier: float = 1.0
    wind_erosion_multiplier: float = 1.0  # wind-driven weathering (erosion.WEATHERING_COEFFICIENT term)
    ocean_erosion_multiplier: float = 1.0  # submarine + coastal (wave/frost) sea-side erosion
    coastal_leveling_multiplier: float = 1.0  # the symmetric near-shore planation grind -- a prime long-run land drain
    glacier_erosion_multiplier: float = 1.0  # glacial abrasion + the ice-flattening blur
    seismic_erosion_multiplier: float = 1.0  # earthquake-triggered landsliding in active ranges
    river_deposition_multiplier: float = 1.0  # floodplain/delta retain fraction (erosion.DEPOSITION_FRACTION)
    ocean_deposition_multiplier: float = 1.0  # marine + beach sediment settling onto the shelf
    collision_uplift_multiplier: float = 1.0  # convergent mountain-building rate (rate, not reach)
    collision_uplift_reach_multiplier: float = 1.0  # horizontal reach of the collision/subduction uplift bands
    volcanism_multiplier: float = 1.0  # per-node eruption frequency *and* per-eruption elevation added
    # Fault-based-only: scales the relief individual fault traces stamp onto the terrain (see
    # faults._apply_plate_fault_relief's `rate_scale`). A true no-op in "boundary"
    # fault_deformation_mode -- faults._relief_mode_scales hard-returns 1.0 there regardless of
    # this value, since boundary mode's smooth polygon-edge bands (not fault traces) carry the
    # deformation. Meaningful only in "fault"/"both" mode, where it directly controls how sharp
    # fault-line ridges/scarps read relative to the surrounding boundary swell.
    fault_relief_multiplier: float = 1.0
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
    # docs/simulation-model.md#wind-model and GitHub issue #118 -- so it's the default; opt into "cfd"
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

    def record_stats(self, force: bool = False) -> None:
        """Append this instant's `stats.compute_stats(self)` snapshot to `stats_history` --
        called once at the end of `finish_generation` (`force=True`) and once at the end of
        every `step_world` (see `stats_history`'s own comment for why this lives on World at
        all now). Deduped by `elapsed_years`, the same guard the frontend's own client-side
        accumulation used to apply, so a `years=0` step (or any other call that leaves
        elapsed_years unchanged) doesn't pile up a redundant entry.

        `stats.compute_stats` reads `climate.compute_climate_cached`, which computes a fresh
        (~50ms) climate snapshot on demand whenever `climate_cache` is empty rather than ever
        skipping it -- fine as a one-time generation cost (`force=True`, unconditionally), but
        *not* something an ordinary step should trigger purely to record a stats snapshot: a
        `simulate_climate_biomes=False` run deliberately skips climate every step to run fast,
        and `climate_cache is None` the whole time it does, so `step_world`'s own call
        (`force=False`) simply skips recording rather than undoing that skip's whole point.

        A `force=True` call that finds `climate_cache` still empty (generation) restores it to
        `None` afterward -- `compute_climate_cached` itself stashes whatever it computes back
        onto `world.climate_cache` as a side effect, and several call sites elsewhere
        (coastline.py, climate.py's vegetation-transpiration source) read `climate_cache is
        None` as "no step has run yet"; letting this method's own one-off snapshot leak into
        that field would make a freshly generated, never-stepped world look like it had.

        A local import: stats.py imports World for its own type hint, so importing it back at
        module scope here would be circular."""
        if not force and self.climate_cache is None:
            return
        from . import stats

        if self.stats_history and self.stats_history[-1]["elapsed_years"] == self.elapsed_years:
            return
        had_no_cache = self.climate_cache is None
        snapshot = stats.compute_stats(self)
        if had_no_cache:
            self.climate_cache = None
        self.stats_history.append(snapshot)

    def record_removed_points(self, points_xyz: np.ndarray, plate_id: int) -> None:
        """Append one `removed_points_log` entry per point in `points_xyz` (world-frame, shape
        (n, 3)) at the current `elapsed_years`, against `plate_id` -- the plate that owned them
        just before they were dropped (retreat, interior-subduction carve, merge absorption,
        defrag stripping, or the whole plate being removed as defunct). A no-op for an empty
        array, so every call site can call this unconditionally rather than guarding on "did
        anything actually get removed this call." See removed_points_log's own comment for why
        this can't just be another ElevationLine field."""
        if len(points_xyz) == 0:
            return
        for point in points_xyz:
            self.removed_points_log.append((point, self.elapsed_years, plate_id))
        overflow = len(self.removed_points_log) - MAX_REMOVED_POINTS_LOG
        if overflow > 0:
            del self.removed_points_log[:overflow]

    def log_corner_notch(self, entry: dict) -> None:
        """Append one structured decision record from `_fill_corner_notch_frontier` -- a no-op unless
        `debug_diagnostics` is on, so a caller can build `entry` unconditionally without
        worrying about cost on an ordinary (non-debugging) world; see this method's own
        callers in lithosphere_plate.py for the exact guard-then-build pattern that keeps this
        genuinely free when diagnostics are off. Stamps `elapsed_years` onto the entry itself
        so a caller doesn't need to repeat it."""
        if not self.debug_diagnostics:
            return
        entry["elapsed_years"] = self.elapsed_years
        self.corner_notch_log.append(entry)
        overflow = len(self.corner_notch_log) - MAX_CORNER_NOTCH_LOG_LENGTH
        if overflow > 0:
            del self.corner_notch_log[:overflow]

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
    voronoi_points: int | None = None,
    sketch: worldsketch.SketchMasks | None = None,
    premade_world_id: str | None = None,
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
    convex initial plate outlines; 0 is the old one-cell-per-plate tiling. `voronoi_points`
    (the UI's "Voronoi points" Advanced-settings slider) is the total seed-point count and,
    when given, overrides `extra_sites_per_plate` -- see lithosphere_plate.generate_plates.
    `sketch` (the
    "Human-made" Generate World tab, see worldsketch.py) replaces noise-driven land/sea and
    random plate-site placement with a drawn/loaded coastline -- see
    lithosphere_plate.generate_plates's own `sketch` param for what it changes; `None` (every
    caller before this parameter existed) is unaffected. `premade_world_id` (`"earth"`,
    `"pangaea"`, or `"got"` -- the "Premade worlds" tab; see lithosphere_plate.generate_plates's
    own param of the same name for what it changes there, real-geography relief regions for
    every one of the three) additionally swaps the random mantle convection centers below for
    ones fit to reproduce known (or, for Pangaea, directionally-approximated) real plate
    motion, for `"earth"`/`"pangaea"` specifically -- see real_plates.py's
    `fit_mantle_centers`. `"got"` has no real-world motion to fit to and keeps the ordinary
    random centers."""
    plates = generate_plates(
        seed,
        num_plates,
        continental_fraction,
        land_fraction,
        node_density,
        extra_sites_per_plate=extra_sites_per_plate,
        voronoi_points=voronoi_points,
        sketch=sketch,
        premade_world_id=premade_world_id,
    )
    rng = np.random.default_rng(seed)
    if premade_world_id in ("earth", "pangaea"):
        # Local import, same reasoning as generate_plates' own (avoid paying for
        # real_plates.py's data-file parsing on every other generation path). "got" has no
        # real-world motion to fit to, so it keeps the ordinary random centers below.
        from . import real_plates

        if premade_world_id == "earth":
            real_plate_list = real_plates.load_major_plates()
            candidates = real_plates.load_ridge_trench_candidates()
        else:
            real_plate_list = real_plates.pangaea_real_plates()
            candidates = real_plates.pangaea_ridge_trench_candidates(sketch, rng)
        mantle_centers = real_plates.fit_mantle_centers(real_plate_list, candidates, rng)
    else:
        mantle_centers = mantle.generate_convection_centers(rng, n_centers=num_mantle_centers)

    world = World(
        seed=seed,
        plates=plates,
        mantle_centers=mantle_centers,
        next_plate_id=len(plates),
        axial_tilt_deg=axial_tilt_deg if axial_tilt_deg is not None else DEFAULT_AXIAL_TILT_DEG,
        premade_world_id=premade_world_id,
        node_density=node_density,
        climate_density=climate_density,
        fluid_density=fluid_density,
    )
    if initial_soil_maturity is not None:
        geology.seed_initial_soil(world.plates, seed, initial_soil_maturity)

    n_continents = sum(1 for p in plates if p.crust_type == "continental")
    finish_generation(world, f"World generated with {len(plates)} plates ({n_continents} continental).")
    return world


def finish_generation(world: World, log_message: str) -> None:
    """Shared tail every world-generation path (`generate_world`, `debug_worlds.
    generate_debug_world`) runs once `world.plates` is populated: bootstrap the permanent
    atmosphere wind-solver state, snapshot the eustatic water budget from the flat starting
    sea level, and log the one generation event. Pulled out of `generate_world` so the
    "Debugging Worlds" tab's scripted scenarios get the exact same bootstrap without
    duplicating it."""
    # See World.atmosphere_cfd_state's own comment for why it's populated here,
    # unconditionally, rather than lazily. `terrain` is the diagnostic bootstrap snapshot the
    # wind solver is seeded from (compute_wind's own latitude-banded field), before the state
    # it will read from every step after exists.
    height, width = climate.grid_dimensions(world.climate_density)
    terrain = climate.compute_climate(world, height, width)
    world.atmosphere_cfd_state = atmosphere_cfd.init_atmosphere_cfd(world, terrain)

    # Snapshot the ocean water volume from the flat starting sea level -- from here on
    # step_world re-solves sea_level_m against this fixed budget every step (see eustasy.py).
    eustasy.initialize_water_budget(world)

    world.log_event(log_message)
    world.record_stats(force=True)  # the elapsed_years=0 baseline entry -- see World.stats_history


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
    `LithospherePlate.shift(world, years)` for every plate (refit Euler pole from torque
    balance, rotate rigidly), then `LithospherePlate.deform(world, other_plates, years, D)` for every plate in a freshly
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
    # The render path's cached node-cloud k-d tree (see World.node_kdtree_cache) and its
    # positions-only sibling shared with climate.py (World.node_position_tree_cache) are both
    # a pure function of node positions, which shift()/deform()/topology changes below are
    # about to move -- drop both now so the first render after this step rebuilds
    # node_kdtree_cache fresh (reusing node_position_tree_cache's tree if climate.py already
    # rebuilt that one this step -- see its own docstring). (land_kdtree_cache is reset
    # separately, inside the simulate_climate_biomes block, since only that path reads it.)
    world.node_kdtree_cache = None
    world.node_position_tree_cache = None
    world.node_healpix_index_cache = None
    world.node_kdtree_relief_cache = None
    world.node_hillshade_cache = None
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
        # Intraplate faults: age/spawn/retire and apply their own relief, on top of (never
        # replacing) deform()'s boundary classification -- see faults.py. Before topology
        # changes so a fresh fault's relief is in place when merge/split geometry is judged.
        faults.update_faults(world, years)
        # Lateral magma export (magma_transport.py) just banked this step's own share of
        # exported convergent-boundary melt into world.pending_magma_parcels above (inside
        # deform()) -- bank the elapsed time alongside it so the eventual transport-pass firing
        # below can rate-cap its deposit against the *banked* interval, not just its own step.
        world.magma_transport_banked_years += years
    world.elapsed_years += years
    if world.simulate_plate_movement:
        for message in merge_split.apply_topology_changes(world, years):
            world.log_event(message)
        # Re-home faults onto surviving plates after any merge/split, drop subducted ones.
        faults.reconcile_faults(world)
        # Stamp/clear ElevationLine.overlap_onset_years and advance World.overlap_progress
        # against this step's final geometry (see merge_split.update_overlap_tracking /
        # docs/debugging.md).
        merge_split.update_overlap_tracking(world, years)
        # Whole-sphere coverage maintenance: spawn new oceanic crust into any region no plate
        # has reached in a long time (see gaps.py) -- e.g. ocean floor a fully-subducted
        # plate vacated with no neighbour left nearby to grow into it. Gated to the same
        # cadence as defragment_plates above (a whole-world pass, not needed every step).
        if world.steps_taken % gaps.GAP_FILL_INTERVAL_STEPS == 0:
            for message in gaps.fill_gaps_by_growing_neighbours(world):
                world.log_event(message)
            # Gap-age diagnostic (see docs/debugging.md's overlapAge section): reconciles
            # world.gap_tracks against this step's uncovered-lattice clusters at the same
            # cadence as fill_gaps_by_growing_neighbours above, since both are the same
            # whole-sphere sweep -- see gaps.reconcile_gap_tracks.
            gaps.reconcile_gap_tracks(world)
        # Lateral magma transport (GitHub issue #205, magma_transport.py): another whole-sphere
        # pass, so gated the same way and (deliberately) placed after topology has fully
        # settled for this step -- same reasoning as gaps.py's own placement here, since this
        # pass writes into destination lines by (plate_id, line_index) and running it earlier
        # could target a line that subducts, splits, or merges away in this same step.
        if world.steps_taken % magma_transport.MAGMA_TRANSPORT_INTERVAL_STEPS == 0:
            for message in magma_transport.run_magma_transport(world, world.magma_transport_banked_years / 1_000_000.0):
                world.log_event(message)
            world.magma_transport_banked_years = 0.0

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

    # Eustatic sea level: re-solve world.sea_level_m against this step's final hypsometry,
    # holding the conserved ocean water volume fixed (see eustasy.py). Unconditional -- both
    # tectonics and erosion reshape the basins, and even a movement-and-climate-off step
    # should keep sea level self-consistent if a control just changed the water budget.
    eustasy.update_sea_level(world)

    world.record_stats()



"""V2 world orchestration: `generate_world_v2`/`step_world_v2`, the direct analogues of
`world.py`'s `generate_world`/`step_world`.

Reuses v1's own `World` dataclass unchanged (not a separate `WorldV2` class) -- `World` is a
plain dataclass with no runtime type enforcement on its `plates`/`ocean_cfd_state`/
`atmosphere_cfd_state` fields, and every module this needs (`erosion.py`, `hydrology.py`,
`lakes.py`, `coastline.py`, `biomes.py`, `geology.py`, `stats.py`, `render_image.py`'s
non-CFD-native views, `persistence.py`) only ever reads it through `Plate`'s own abstract
interface or plain fields (`seed`, `sea_level_m`, `node_density`, ...) -- never a v1-specific
type. Populating `world.plates` with `LithospherePlate`s and `world.ocean_cfd_state`/
`atmosphere_cfd_state` with the `*_v2` HEALPix state classes is exactly as valid as v1's own
population of those same fields; this is the mechanism the plan's "why subclass
PlateWithLines" section describes, extended one level up to `World` itself.

**Multi-rate coupling, scoped down from the spec's full time-averaging design** (see the
plan's own implementation-order note): fluid dynamics is advanced *before* erosion each step
(so erosion/hydrology read post-substep, not pre-substep, wind/precipitation -- the ordering
half of the spec's Section 5.2 loop), and updated topography feeds back into the fluid grid's
boundary conditions for the *next* step (`refresh_forcing`, same pattern v1 already uses) --
but rather than explicitly time-averaging precipitation/wind-stress/sediment-capacity across
the fluid sub-stepping cycle, erosion/hydrology read the *end-of-cycle* snapshot
(`climate.compute_climate`, the same instantaneous-snapshot mechanism v1 already treats as
"up to one step stale, an accepted simplification" -- see `World.climate_cache`'s own
docstring). A real time-average would need `erosion.py` itself (a large, intricate v1 module)
to accept a pre-averaged field set instead of computing its own snapshot -- out of scope for
this pass; flagged here rather than silently simplified.

`bathymetry.apply_bathymetry` is dropped entirely (Section 2.2: isostasy supersedes it).
`merge_split.py` (merge/consumption/split) and `volcanism.py`/`geology.py` (resource
formation) are reused completely unchanged -- see plates_v2.py's own docstring for why split's
trigger (still `mantle`-flow-residual-based, not yet a stress-rupture trigger) is a further
scope simplification, not an oversight.
"""

from __future__ import annotations

import numpy as np

from .. import climate, erosion, geology, mantle, merge_split, volcanism
from ..plates import gather_node_positions
from ..world import DEFAULT_AXIAL_TILT_DEG, DEFAULT_MANTLE_CENTERS, World
from . import atmosphere_cfd_v2, ocean_cfd_v2, plates_v2


def generate_world_v2(
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
) -> World:
    plates = plates_v2.generate_plates_v2(seed, num_plates, continental_fraction, land_fraction, node_density)
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

    height, width = climate.grid_dimensions(world.climate_density)
    terrain = climate.compute_climate(world, height, width)
    world.atmosphere_cfd_state = atmosphere_cfd_v2.init_atmosphere_cfd(world, terrain)
    world.ocean_cfd_state = ocean_cfd_v2.init_ocean_cfd(world, terrain, world.atmosphere_cfd_state)

    world.log_event(
        f"World generated with {len(plates)} plates "
        f"({sum(1 for p in plates if p.crust_type == 'continental')} continental)."
    )
    return world


def _advance_fluid_dynamics_v2(world: World, node_cloud) -> None:
    """Unlike v1's own `_advance_fluid_dynamics` (which runs *after* erosion and so can
    reuse erosion's own freshly-populated `world.climate_cache` when the two share a
    resolution), this runs *before* erosion each step (see step_world_v2's own docstring
    note on reordering) -- `world.climate_cache` would still hold last step's snapshot at
    this point, so a fresh one is always computed here regardless of whether fluid_density
    matches climate_density. `erosion.apply_erosion` computes its own snapshot right after
    this returns (on the same, still-unchanged post-tectonics `world.plates`), so this is one
    redundant `compute_climate` call per step in exchange for correct fluid forcing --
    acceptable at this pass's scope, not a correctness bug."""
    terrain = climate.compute_climate(world, *climate.grid_dimensions(world.fluid_density), node_cloud=node_cloud)
    atmosphere_cfd_v2.refresh_forcing(world, world.atmosphere_cfd_state, terrain)
    atmosphere_cfd_v2.step_atmosphere_cfd(world, world.atmosphere_cfd_state, atmosphere_cfd_v2.SECONDS_PER_TECTONIC_STEP)
    ocean_cfd_v2.refresh_forcing(world, world.ocean_cfd_state, terrain)
    ocean_cfd_v2.step_ocean_cfd(world, world.ocean_cfd_state, ocean_cfd_v2.SECONDS_PER_TECTONIC_STEP)


def step_world_v2(world: World, years: float) -> None:
    if world.simulate_plate_movement:
        distances = {plate.plate_id: plate.shift(world, years) for plate in world.plates}
        order = list(world.plates)
        np.random.default_rng((world.seed, round(world.elapsed_years))).shuffle(order)
        for plate in order:
            others = [p for p in world.plates if p.plate_id != plate.plate_id]
            plate.deform(world, others, years, distances[plate.plate_id])
    world.elapsed_years += years
    if world.simulate_plate_movement:
        # TODO(v2): merge_split.maybe_split_plate's own trigger (does one rigid rotation fit
        # this plate's *static mantle-flow field* sample poorly?) is a v1-kinematic criterion
        # that doesn't really map onto v2, where plate motion is no longer fit to that field
        # at all -- it's reused unchanged here for now (see plates_v2.py's module docstring).
        # A torque-driven model has a more natural, physically-motivated trigger available:
        # spatial heterogeneity in local strain rate/stress (rheology.py's own per-node
        # closing-rate-driven yield state) across a plate's own footprint -- rework the
        # trigger to key off that instead of mantle.fit_euler_pole's residual.
        for message in merge_split.apply_topology_changes(world, years):
            world.log_event(message)

    erosion_result = None
    if world.simulate_climate_biomes:
        node_cloud = gather_node_positions(world.plates)
        world.land_kdtree_cache = None
        # Fluid dynamics advances *before* erosion/hydrology here (unlike v1's own ordering,
        # which runs erosion first) -- the multi-rate loop's own "topographic feedback ->
        # fluid sub-stepping -> [time-averaged] fields -> geological update" sequencing (see
        # module docstring for the one piece of that loop this pass doesn't fully implement).
        _advance_fluid_dynamics_v2(world, node_cloud)
        erosion_result = erosion.apply_erosion(world, years, node_cloud=node_cloud)
    if world.simulate_plate_movement:
        for message in volcanism.apply_volcanic_activity(world, years):
            world.log_event(message)
    if erosion_result is not None:
        geology.apply_resource_formation(world, years, erosion_result)

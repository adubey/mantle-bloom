"""Per-node terrain constants and the node-density/spacing choices that govern a plate's
surface.

Node spacing (`TARGET_LINE_SPACING_RAD`, `line_spacing_rad`, `NODE_DENSITY_CHOICES`),
elevation bounds, the per-node elevation-change provenance and crust-type codes, volcano
constants, and the `ElevationPoint` per-node view protocol live here, shared by every module
that reads or writes plate surface nodes. `iter_local_lattice` sweeps a plate-local sampling
lattice for whole-sphere coverage tests and rendering; it is a derived sampling grid, not plate
topology. See docs/simulation-model.md for the full design writeup."""

from __future__ import annotations

from typing import Protocol

import numpy as np

from . import geometry
from .surface_fields import SURFACE_FIELDS

PLANET_RADIUS_KM = 6371.0

# Halving this doubles resolution in each dimension, i.e. ~4x the nodes per plate. Several other modules define *absolute node-count*
# thresholds (not distances, which already scale automatically as multiples of
# TARGET_LINE_SPACING_RAD) that represent a physical area or distance in terms of the *old*
# density -- those were rescaled alongside this (merge_split.SPLIT_MIN_NODES,
# gaps.MIN_GAP_POINTS/MAX_ABSORB_NODES_PER_PLATE_PER_CALL by ~4x for area,
# boundary.MAX_EXTEND_NODES_PER_STEP by ~2x for a 1D distance) -- see each for the reasoning.
# This is the reference value for the default node_density=1.0 -- see line_spacing_rad below
# for how a world's own chosen density (World.node_density, set once at generation and read
# by every module in this same list for the rest of that world's life) scales it at runtime,
# now that density is a per-world user choice rather than a hardcoded, one-off code change.
TARGET_LINE_SPACING_KM = 125.0
TARGET_LINE_SPACING_RAD = TARGET_LINE_SPACING_KM / PLANET_RADIUS_KM

# UI-facing choices for World.node_density -- a discrete set (not a free-form slider) since
# there's no natural continuous unit for "how many points," only "how many times as many."
# 6.0 (the UI's "High" choice, 1.5x the default multiplier) sits above the default -- more
# boundary geometry than 4.0 for the rare step where even that isn't enough, at a real per-step
# cost (see line_spacing_rad: node count scales with the *square* of this multiplier, so 1.5x
# the density is already 2.25x the nodes). 2.0 (half the default multiplier) is a
# lower-resolution middle ground -- fewer nodes than the default, so plate-movement-only
# stepping (World.simulate_plate_movement, World.simulate_climate_biomes off) runs faster,
# without dropping all the way to 1.0's much coarser boundary geometry. 0.5 (an eighth of the
# default) is coarser still -- the fastest, lowest-fidelity option, useful where even 1.0's
# geometry is more than a given step needs.
NODE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0, 6.0)
DEFAULT_NODE_DENSITY = 4.0

# Physical elevation bounds every module that modifies elevation clips against (boundary.py,
# erosion.py, volcanism.py) -- kept in one place so
# they can't drift out of sync between call sites.
MIN_ELEVATION_M = -11000.0
MAX_ELEVATION_M = 9000.0

# --- Elevation-change provenance ("why did this node's elevation last move") -------------
#
# The `elev_change_reason` surface field (see surface_fields.SURFACE_FIELDS) holds
# one of these integer codes per node: the dominant process that last moved that node's
# `elevation` by a non-trivial amount. The tectonic engine stamps the tectonic codes,
# volcanism.py the eruption code, erosion.py the geomorphic codes -- each only where its own
# per-step delta clears `ELEV_CHANGE_MIN_DELTA_M`, so a quiescent low-relief node keeps
# whatever last genuinely shaped it (often NONE -- untouched since generation) rather than
# being relabelled every step by sub-metre erosion noise. Diagnostic only; nothing in the
# physics reads it back. Surfaced by render_image.py's "elevReason" debug view -- built to
# answer "why is so much of this world flat: never uplifted, or actively planed down?".
ELEV_CHANGE_MIN_DELTA_M = 2.0

# Erosion runs every step and brushes almost every land node a little, so without a guard it
# would relabel an actively-rising mountain belt "erosion" purely because a few m/step of
# rain wash also happened there -- burying the tectonic signal the view exists to show. A
# structural code (ELEV_CHANGE_COLLISION..ELEV_CHANGE_VOLCANO, re-stamped by deform() every
# step the belt is still active) is therefore only overwritten by a geomorphic code when this
# step's net geomorphic change is itself large -- at least this rate, comparable to a real
# uplift increment, not ordinary background wash. Non-structural prior codes (NONE, or an
# earlier geomorphic one) are overwritten on any move past ELEV_CHANGE_MIN_DELTA_M.
ELEV_CHANGE_STRUCTURAL_OVERRIDE_M_PER_MYR = 100.0

ELEV_CHANGE_NONE = 0  # untouched since initial generation (or below the per-step threshold)
ELEV_CHANGE_COLLISION = 1  # continent-continent near-field collision uplift
ELEV_CHANGE_COLLISION_FAR_FIELD = 2  # legacy broad far-field collision uplift, retired in #206
ELEV_CHANGE_SUBDUCTION_ARC = 3  # oceanic-under-continental volcanic-arc uplift
ELEV_CHANGE_TRENCH = 4  # subducting oceanic plate's own trench subsidence
ELEV_CHANGE_TRANSFORM = 5  # transform-boundary pressure-ridge uplift
ELEV_CHANGE_RIFT = 6  # divergent ridge/rift relaxation toward the spreading target
ELEV_CHANGE_NEW_CRUST = 7  # brand-new crust inserted at a growing spreading edge
ELEV_CHANGE_VOLCANO = 8  # volcanic eruption (deform-spawned or ongoing lifecycle)
ELEV_CHANGE_EROSION = 9  # subaerial erosion (rain/river/weathering/glacier/seismic)
ELEV_CHANGE_DEPOSITION = 10  # fluvial / wind / glacial-transport sediment deposition
ELEV_CHANGE_COASTAL_LEVELING = 11  # wave-cut planation + sheltered-shelf infill toward sea level
ELEV_CHANGE_MARINE = 12  # submarine erosion / marine sediment spread on the sea floor
ELEV_CHANGE_GLACIAL_FLATTEN = 13  # glacier flattening (broad sub-ice smoothing)
ELEV_CHANGE_LAKE_SILT = 14  # lake / endorheic-basin siltation raising a basin floor
# Intraplate fault relief (see faults.py) -- a fault line that is *not* a plate boundary,
# so these are distinct from ELEV_CHANGE_TRANSFORM (which is boundary-only). Structural, so
# they get the same erosion-override protection as the boundary codes above.
ELEV_CHANGE_FAULT_NORMAL = 15  # extensional graben subsidence / footwall-shoulder uplift
ELEV_CHANGE_FAULT_REVERSE = 16  # intraplate thrust / fold-belt uplift away from a boundary
ELEV_CHANGE_FAULT_STRIKE_SLIP = 17  # strike-slip transpressional ridge / transtensional sag
# A broad, low-relief apron an eruption spreads around itself (see volcanism.py's
# VOLCANIC_PLAIN_*) -- distinct from ELEV_CHANGE_VOLCANO, which is reserved for the sharp
# point bump at the vent itself.
ELEV_CHANGE_VOLCANIC_PLAIN = 18
# Lateral magma transport deposit (GitHub issue #205, magma_transport.py): a cross-plate melt
# deposit landing far from the collision boundary that generated it, so it can't ride along
# inside the tectonic engine's own reason-stamping block the way delamination
# melt does (that one shares ELEV_CHANGE_COLLISION only because it's applied *inside* the same
# deform() call). Structural, so it gets the same erosion-override protection as the other codes.
ELEV_CHANGE_LATERAL_MAGMA = 19

# Human-readable label per code, index == code -- kept here (not in render_image.py or the
# frontend) as the single source both sync against, same precedent as biomes.BIOME_NAMES.
ELEV_CHANGE_LABELS = (
    "Unchanged since generation",
    "Continental collision uplift",
    "Legacy far-field collision uplift",
    "Subduction-arc uplift",
    "Oceanic trench subsidence",
    "Transform pressure ridge",
    "Divergent rift / ridge",
    "New crust at spreading edge",
    "Volcanic eruption",
    "Erosion (worn down)",
    "Sediment deposition",
    "Coastal planation / infill",
    "Submarine erosion / sediment",
    "Glacial flattening",
    "Lake / basin siltation",
    "Fault: normal (graben)",
    "Fault: reverse (thrust)",
    "Fault: strike-slip",
    "Volcanic plain",
    "Lateral magma transport deposit",
)

# --- Per-node crust type (the `crust_type_code` surface field) ---------------------------
#
# A plate's `crust_type` ("oceanic"/"continental") is the *usual* case, but real crust is
# genuinely composite: a rift can erupt continental-type magma while sitting on an oceanic
# plate (a volcanic island breaching the surface) or oceanic-type magma while sitting on a
# continental plate (a drowned continental margin finally thinning through to true seafloor),
# and a plate spawned to fill a whole-sphere gap can straddle both if it borders a continent.
# CRUST_TYPE_INHERIT (0, the field's registry default) means "same as the owning plate," which is both the safe
# backward-compatible reading for a pre-existing save and the exactly-correct reading for
# every node created by ordinary generation/growth/merge/split -- those are never stamped
# otherwise, so nothing about their existing (calibrated) physics changes. Only decompression
# melting (lithosphere_plate.py) and gap-fill (gaps.py) ever stamp an explicit value.
CRUST_TYPE_INHERIT = 0
CRUST_TYPE_OCEANIC = 1
CRUST_TYPE_CONTINENTAL = 2


def effective_is_continental_from_codes(codes: np.ndarray, plate_is_continental: bool) -> np.ndarray:
    """Per-node bool: is this node's crust actually continental, resolving
    CRUST_TYPE_INHERIT against the owning plate's own `crust_type` and taking an explicit
    CRUST_TYPE_OCEANIC/CONTINENTAL code at face value. Takes a raw code array directly (e.g.
    `Plate.collect("crust_type_code")`)."""
    return np.where(codes == CRUST_TYPE_INHERIT, plate_is_continental, codes == CRUST_TYPE_CONTINENTAL)


# Shared between volcanism.py (per-step eruption rolling for every existing volcano node)
# and the tectonic engine (spawning a brand-new volcano when a rift has stretched too thin to
# keep filling with plain ridge/rift crust) -- kept here, rather than in volcanism.py, so the
# engine can use them without importing volcanism.py (which itself imports from plates.py).
VOLCANO_ACTIVE_MIN_YEARS = 100_000
VOLCANO_ACTIVE_MAX_YEARS = 1_000_000
# A single eruption's land contribution. Volcano nodes only ever spawn where a rift has
# stretched too thin (plates.py's own deform() -- there is no arc/hotspot spawning yet, see
# GitHub issue #120, "Land fraction slowly declines"), so a volcano's whole visible-land contribution
# is however many of these `ERUPTION_RATE_PER_MYR`-paced eruptions it manages during its short
# (VOLCANO_ACTIVE_MIN/MAX_YEARS) active life -- at the old 100 m this topped out well under a
# real shield/stratovolcano's actual relief and mostly read as erosion noise by the next few
# steps (measured: across a 188 My save, volcano nodes averaged *below* sea level). Raised
# 2026-09-04 so a volcano that does erupt a few times actually builds a lasting peak.
ERUPTION_ELEVATION_M = 300.0

# A shield-volcano/flood-basalt apron: real effusive eruptions don't just build a point peak,
# they spread lava broadly around the vent (Hawaiian shield flanks, the Deccan/Siberian Traps
# at the extreme end) -- a real, physically distinct landform from the sharp vent bump above,
# so it gets its own reason code (ELEV_CHANGE_VOLCANIC_PLAIN) and shows up separately in the
# "Last elevation change" view. Applied by volcanism.py alongside the existing point bump,
# same eruption roll, tapering linearly from the vent to zero at VOLCANIC_PLAIN_REACH_KM --
# see `_spread_volcanic_plains`. Deliberately much lower than ERUPTION_ELEVATION_M (a plain is
# low-relief by definition) and reaches roughly one neighbour ring at default node density, so
# it reads as a genuine broad apron rather than an isolated spike.
VOLCANIC_PLAIN_REACH_KM = 120.0
VOLCANIC_PLAIN_ELEVATION_M = 60.0

def line_spacing_rad(node_density: float) -> float:
    """The line spacing (radians) that gives a plate ~node_density times as many nodes as
    the default TARGET_LINE_SPACING_RAD would. Node count for a fixed physical area scales
    with the *square* of resolution (see TARGET_LINE_SPACING_KM's own comment -- halving
    spacing quadruples node count), so this divides by sqrt(node_density), not
    node_density itself. Every module that derives a distance threshold or an absolute
    node-count cap from TARGET_LINE_SPACING_RAD calls this (with the world's own
    node_density) instead of reading the bare module constant directly, so that a world
    generated at a non-default density stays self-consistent for its entire life -- not just
    at generation, but through every later gap-fill/merge/split/volcanism pass
    too (each of those modules' own docstrings/comments explain why its own particular
    thresholds need this)."""
    return TARGET_LINE_SPACING_RAD / np.sqrt(node_density)


# Shared geometry-tolerance constants for "is this lattice point already covered by real
# crust" / "is this node still part of the same contiguous patch" checks -- kept in one place
# (rather than each caller defining its own copy) so gaps.py's whole-sphere sweep,
# merge_split.py's defragmentation pass, and quad_tectonics.py's boundary growth all agree on
# the same tolerances.

# A lattice point counts as "covered" if some real node sits within this multiple of node
# spacing of it -- comfortably more than one spacing so ordinary per-step catch-up growth
# isn't mistaken for a genuine void, but tight enough that a real neighbouring patch is never
# missed. See gaps.py's own module docstring for the whole-sphere case this was written for.
COVERAGE_RADIUS_MULT = 1.5

# Two nodes count as connected (the same contiguous patch) if they're within this multiple of
# node spacing of each other -- an edge neighbour is ~1x spacing away and a diagonal
# neighbour ~1.4x, so 2.5x comfortably links a genuinely contiguous patch while still
# separating two lobes across a real (>~300km) subduction gap. Validated against real saved
# worlds by merge_split.py's defragmentation pass: every healthy plate comes back as a single
# component at this radius.
DEFRAG_CONNECT_RADIUS_MULT = 2.5


class ElevationPoint(Protocol):
    """A single surface node's data -- structural (not a base class), so a
    representation-specific backing (`sparse_quad_patch.ElevationPointInPatch`) can satisfy it
    without sharing a base.

    `phi`/`get_theta()` are this point's fixed plate-local position, get-only: no code in this
    simulation ever moves a single node in place -- position changes always go through a
    topology change on the whole plate, so a per-point position setter would just invite a
    caller to silently desync the surface's own ordering. `elevation` and every other
    `surface_fields.SURFACE_FIELDS` name get real setters, since per-node *value* mutation
    (an eroded elevation, a grown channel, a newly lit volcano) is exactly what per-step
    simulation passes do."""

    @property
    def phi(self) -> float: ...
    def get_theta(self) -> float: ...

    def get_elevation(self) -> float: ...
    def set_elevation(self, value: float) -> None: ...

    def get_channel_depth(self) -> float: ...
    def set_channel_depth(self, value: float) -> None: ...

    def get_channel_width(self) -> float: ...
    def set_channel_width(self, value: float) -> None: ...

    def get_lake_depth(self) -> float: ...
    def set_lake_depth(self, value: float) -> None: ...

    def get_glacier_depth(self) -> float: ...
    def set_glacier_depth(self, value: float) -> None: ...

    def get_silt_depth(self) -> float: ...
    def set_silt_depth(self, value: float) -> None: ...

    def get_is_volcano(self) -> bool: ...
    def set_is_volcano(self, value: bool) -> None: ...

    def get_volcano_active_years_remaining(self) -> float: ...
    def set_volcano_active_years_remaining(self, value: float) -> None: ...

    def get_soil_depth(self) -> float: ...
    def set_soil_depth(self, value: float) -> None: ...

    def get_soil_mineral_content(self) -> float: ...
    def set_soil_mineral_content(self, value: float) -> None: ...

    def get_soil_organic_content(self) -> float: ...
    def set_soil_organic_content(self, value: float) -> None: ...

    def get_coal_deposit_m(self) -> float: ...
    def set_coal_deposit_m(self, value: float) -> None: ...

    def get_oil_gas_deposit_m(self) -> float: ...
    def set_oil_gas_deposit_m(self, value: float) -> None: ...

    def get_mineral_deposit_m(self) -> float: ...
    def set_mineral_deposit_m(self, value: float) -> None: ...

    def get_elev_change_reason(self) -> float: ...
    def set_elev_change_reason(self, value: float) -> None: ...

    def get_overlap_onset_years(self) -> float: ...
    def set_overlap_onset_years(self, value: float) -> None: ...

    def get_node_created_years(self) -> float: ...
    def set_node_created_years(self, value: float) -> None: ...


def _point_field_getter(name: str):
    def getter(self) -> float:
        return self._field_array(name)[self._index]

    getter.__name__ = f"get_{name}"
    return getter


def _point_field_setter(name: str):
    def setter(self, value) -> None:
        self._field_array(name)[self._index] = value

    setter.__name__ = f"set_{name}"
    return setter


def install_point_field_accessors(cls: type) -> type:
    """Class decorator attaching `get_theta` plus a `get_<name>`/`set_<name>` pair for every
    `surface_fields.SURFACE_FIELDS` name to `cls`, which need only provide
    `_field_array(self, name) -> np.ndarray` and an `_index` attribute -- so an
    `ElevationPoint` implementation stays wired to the one field registry rather than
    hand-writing (and risking silently forgetting) a method per field."""
    setattr(cls, "get_theta", _point_field_getter("theta"))
    for _name in SURFACE_FIELDS:
        setattr(cls, f"get_{_name}", _point_field_getter(_name))
        setattr(cls, f"set_{_name}", _point_field_setter(_name))
    return cls


def iter_local_lattice(frame: np.ndarray, spacing_rad: float = TARGET_LINE_SPACING_RAD):
    """Sweep a full plate-local (phi, theta) lattice at `spacing_rad` resolution, yielding
    (phi, theta_candidates, world_pts) per row. A derived sampling grid: used by whole-sphere
    coverage sweeps (gaps.py) and, at a resolution independent of the physical node spacing,
    by the render-grid sweep (see render_image.py's
    _render_grid_arrays) that gives the rendered map full coverage regardless of how sparse
    the underlying physical data is once projected."""
    max_abs_phi = np.pi / 2 - spacing_rad / 2
    phi_values = np.arange(-max_abs_phi, max_abs_phi, spacing_rad)
    for phi in phi_values:
        dtheta = spacing_rad / max(np.cos(phi), 1e-3)
        n_theta = max(int(np.round(2 * np.pi / dtheta)), 1)
        theta_candidates = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)

        local_pts = geometry.local_xyz(np.full_like(theta_candidates, phi), theta_candidates)
        world_pts = geometry.to_world(frame, local_pts)
        yield float(phi), theta_candidates, world_pts

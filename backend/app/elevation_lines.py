"""Elevation lines: a plate's terrain nodes, and the node-density/spacing choices that
govern them.

An `ElevationLine` sits at a fixed plate-local latitude `phi`, holding elevation (and other
persistent, land-only or lake/volcano/soil/resource) samples at plate-local longitude nodes
`theta`. Because a plate's local (phi, theta) coordinates never change -- only its `frame`
(rotation matrix, local -> world) does, see `Plate` in plates.py -- rotating a plate rigidly
never needs resampling: it's exact for every carried point. See docs/simulation-model.md for
the full design writeup.

`TARGET_LINE_SPACING_RAD`/`line_spacing_rad`/`NODE_DENSITY_CHOICES` (how densely a plate's
lattice is sampled) and `iter_local_lattice`/`build_lines_from_lattice` (sweeping that
lattice to build a fresh set of lines) live here too, rather than in plates.py, since they
only ever produce or describe `ElevationLine`s -- nothing about them depends on how those
lines get bundled into a `Plate`.

This module also owns periodic line regularization. Per-step boundary evolution
(boundary.py) only ever touches the two ends of a line -- inserting at target spacing when
growing, deleting when shrinking -- so interior spacing stays regular on its own. What it
can't fix is spacing that's drifted at a *transform* boundary (nodes sheared along the line
without insertion/deletion) or after several steps' worth of end-growth at a slightly
different rate than the line's original spacing. `regularize_line` re-derives a fresh
evenly-spaced node set spanning each line's *existing* extent (the two endpoints are
preserved exactly -- regularizing never changes where a line's physical edge is, only how
regularly it's sampled) and interpolates elevation onto it.

`spacing_rad` (default `TARGET_LINE_SPACING_RAD`, the reference density) should always be
`line_spacing_rad(world.node_density)` in practice -- `regularize_world_lines` computes it
once per call and threads it down. Without this, a world generated at a non-default density
would regularize itself right back down to the reference density the first time any line's
spacing drifted enough to trigger this pass (every REGULARIZE_INTERVAL_STEPS steps) --
confirmed directly as the failure mode that made a "just build denser lines at generation"
version of a density option pointless within a handful of steps."""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterator, Protocol

import numpy as np

from . import geometry

if TYPE_CHECKING:
    from .plates import PlateWithLines
    from .world import World

PLANET_RADIUS_KM = 6371.0

# Halving this doubles resolution in each dimension (phi rows and theta samples per row),
# i.e. ~4x the nodes per plate. Several other modules define *absolute node-count*
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
# 2.0 (half the default multiplier) is a lower-resolution middle ground -- fewer nodes than the
# default, so plate-movement-only stepping (World.simulate_plate_movement,
# World.simulate_climate_biomes off) runs faster, without dropping all the way to 1.0's much
# coarser boundary geometry. 0.5 (an eighth of the default) is coarser still -- the fastest,
# lowest-fidelity option, useful where even 1.0's geometry is more than a given step needs.
NODE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0)
DEFAULT_NODE_DENSITY = 4.0

# Physical elevation bounds every module that modifies elevation clips against (boundary.py,
# erosion.py, volcanism.py, and this module's own crumpling below) -- kept in one place so
# they can't drift out of sync between call sites.
MIN_ELEVATION_M = -11000.0
MAX_ELEVATION_M = 9000.0

REGULARIZE_INTERVAL_STEPS = 5
IRREGULARITY_TOLERANCE = 1.5  # regularize a line if any gap exceeds this multiple of target

# Shared between volcanism.py (per-step eruption rolling for every existing volcano node)
# and plates.py (PlateWithLines.deform spawning a brand-new volcano when a rift has
# stretched too thin to keep filling with plain ridge/rift crust) -- kept here, rather than
# in volcanism.py, so plates.py can use them without importing volcanism.py (which itself
# imports from plates.py).
VOLCANO_ACTIVE_MIN_YEARS = 100_000
VOLCANO_ACTIVE_MAX_YEARS = 1_000_000
# A single eruption's land contribution -- comparable order of magnitude to
# plates.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR (800 m/Myr) applied over a fraction of a Myr,
# consistent with "a discrete volcanic event" rather than a smooth continuous uplift rate.
ERUPTION_ELEVATION_M = 100.0

# Self-affine scaling exponent used by _crumple_elevation below: real terrain roughened by
# compressing a profile horizontally by k doesn't just get resampled at the new spacing, its
# vertical amplitude grows by roughly k**-CRUMPLE_HURST_EXPONENT (a Hurst exponent -- 0.5 is
# the standard "random walk" / Brownian terrain default used when no better estimate of a
# specific landscape's roughness is available). This is what makes the vulcanism-driven
# density increase that triggers crumpling look like real compression -- ridges pushed
# together get taller, not just thinned out -- rather than plain decimation, which would
# leave peak/valley heights untouched and only make the line coarser.
CRUMPLE_HURST_EXPONENT = 0.5


def line_spacing_rad(node_density: float) -> float:
    """The line spacing (radians) that gives a plate ~node_density times as many nodes as
    the default TARGET_LINE_SPACING_RAD would. Node count for a fixed physical area scales
    with the *square* of resolution (see TARGET_LINE_SPACING_KM's own comment -- halving
    spacing quadruples node count), so this divides by sqrt(node_density), not
    node_density itself. Every module that derives a distance threshold or an absolute
    node-count cap from TARGET_LINE_SPACING_RAD calls this (with the world's own
    node_density) instead of reading the bare module constant directly, so that a world
    generated at a non-default density stays self-consistent for its entire life -- not just
    at generation, but through every later regularize/gap-fill/merge/split/volcanism pass
    too (each of those modules' own docstrings/comments explain why its own particular
    thresholds need this)."""
    return TARGET_LINE_SPACING_RAD / np.sqrt(node_density)


class ElevationLine:
    """A fixed plate-local latitude `phi` holding elevation (and other persistent, land-only
    or lake/volcano/soil/resource) samples at plate-local longitude nodes `theta`.

    Iterating a line (`for point in line`, `line[i]`, `len(line)`) yields `ElevationPointOnLine`
    instances -- one per node, each a live view onto this line's own arrays (see that class
    below). A point's `set_*` methods mutate this line's data in place; that's the one way an
    existing line's *data* changes without going through the whole-array methods below. What
    stays off-limits to per-point mutation is the node set itself (`theta`'s shape/order) --
    for that, use:
    - `replace(...)` swaps in a subset of fields, keeping `theta`'s shape/order untouched
      (the common case -- most steps only ever change `elevation` or one or two persistent
      fields for the same set of nodes, in bulk).
    - `masked(mask)` filters and/or reorders every field together by a boolean mask or
      fancy index (plate split, node removal, node reassignment reordering).
    - `with_new_nodes(theta, elevation)` appends brand-new nodes (zero/False for every
      persistent field -- no history to carry) to the end, unsorted.

    Threading every persistent field through a single generic method here (keyed off
    `OPTIONAL_FIELDS`, one list) is what earlier avoided a real, previously confirmed bug:
    erosion.py's and bathymetry.py's own hand-written reconstruction sites were both written
    before is_volcano/volcano_active_years_remaining existed, so neither passed them through
    -- silently wiping every node's volcanic status to False every single step. A call site
    that constructs a new field-by-field `ElevationLine(...)` directly (bypassing these
    methods) reintroduces exactly that risk."""

    # All persistent, land-only, meters (unless noted), same shape as theta -- see
    # hydrology.py/lakes.py/volcanism.py/geology.py. Because the grid is plate-local and
    # rotates with a plate's `frame` rather than sitting fixed in world space, these ride
    # along for free just by being an ordinary parallel array on this same line, exactly
    # like elevation itself -- no explicit semi-Lagrangian advection needed every step:
    # rotating a plate only ever touches `frame`, never these arrays.
    OPTIONAL_FIELDS = (
        "channel_depth",  # river channel incision, self-reinforcing
        "channel_width",  # river channel width, grows with flow -- see erosion.py
        "lake_depth",  # standing lake water depth
        "glacier_depth",  # accumulated ice, meters ice-equivalent
        # Sediment settled on a lake's own bed, monotonically increasing (never erodes back
        # away, same self-reinforcing character as channel_depth) -- raises the *effective*
        # floor a lake's own depth is measured against without touching real terrain
        # `elevation` itself, see lakes.py's own module docstring for why. Always 0 outside
        # an active lake.
        "silt_depth",
        # Two more of the same "rides along for free" persistent fields, see volcanism.py.
        # is_volcano never reverts to False once set (permanent provenance -- a dormant
        # volcano is still excluded from being redetected as a fresh rift gap);
        # volcano_active_years_remaining is a countdown, 0 once dormant (whether or not
        # is_volcano is set).
        "is_volcano",  # bool
        "volcano_active_years_remaining",  # years
        # Soil, land-only -- see geology.py. Unlike every other field here, these three can
        # both rise *and* fall (real soil forms and erodes), not just accumulate.
        "soil_depth",  # meters, regolith/soil thickness
        "soil_mineral_content",  # [0, 1], weathered/hydrothermal richness
        "soil_organic_content",  # [0, 1], accumulated organic matter
        # Resource deposits -- see geology.py/volcanism.py. All monotonically non-decreasing,
        # the same self-reinforcing "once formed, never erodes back away" convention
        # silt_depth already uses (buried peat/hydrocarbons/ore aren't un-buried by a later
        # climate shift).
        "coal_deposit_m",  # land-only
        "oil_gas_deposit_m",  # ocean-only
        "mineral_deposit_m",  # either -- grown by volcanism.py's own eruptions
    )

    def __init__(
        self,
        phi: float,
        theta: np.ndarray,
        elevation: np.ndarray,
        channel_depth: np.ndarray | None = None,
        channel_width: np.ndarray | None = None,
        lake_depth: np.ndarray | None = None,
        glacier_depth: np.ndarray | None = None,
        silt_depth: np.ndarray | None = None,
        is_volcano: np.ndarray | None = None,
        volcano_active_years_remaining: np.ndarray | None = None,
        soil_depth: np.ndarray | None = None,
        soil_mineral_content: np.ndarray | None = None,
        soil_organic_content: np.ndarray | None = None,
        coal_deposit_m: np.ndarray | None = None,
        oil_gas_deposit_m: np.ndarray | None = None,
        mineral_deposit_m: np.ndarray | None = None,
    ) -> None:
        self._phi = phi
        self._theta = theta
        self._elevation = elevation
        self._channel_depth = channel_depth if channel_depth is not None else np.zeros_like(theta)
        self._channel_width = channel_width if channel_width is not None else np.zeros_like(theta)
        self._lake_depth = lake_depth if lake_depth is not None else np.zeros_like(theta)
        self._glacier_depth = glacier_depth if glacier_depth is not None else np.zeros_like(theta)
        self._silt_depth = silt_depth if silt_depth is not None else np.zeros_like(theta)
        self._is_volcano = is_volcano if is_volcano is not None else np.zeros_like(theta, dtype=bool)
        self._volcano_active_years_remaining = (
            volcano_active_years_remaining if volcano_active_years_remaining is not None else np.zeros_like(theta)
        )
        self._soil_depth = soil_depth if soil_depth is not None else np.zeros_like(theta)
        self._soil_mineral_content = soil_mineral_content if soil_mineral_content is not None else np.zeros_like(theta)
        self._soil_organic_content = soil_organic_content if soil_organic_content is not None else np.zeros_like(theta)
        self._coal_deposit_m = coal_deposit_m if coal_deposit_m is not None else np.zeros_like(theta)
        self._oil_gas_deposit_m = oil_gas_deposit_m if oil_gas_deposit_m is not None else np.zeros_like(theta)
        self._mineral_deposit_m = mineral_deposit_m if mineral_deposit_m is not None else np.zeros_like(theta)

    @property
    def phi(self) -> float:
        return self._phi

    @property
    def theta(self) -> np.ndarray:
        return self._theta

    @property
    def elevation(self) -> np.ndarray:
        return self._elevation

    @property
    def channel_depth(self) -> np.ndarray:
        return self._channel_depth

    @property
    def channel_width(self) -> np.ndarray:
        return self._channel_width

    @property
    def lake_depth(self) -> np.ndarray:
        return self._lake_depth

    @property
    def glacier_depth(self) -> np.ndarray:
        return self._glacier_depth

    @property
    def silt_depth(self) -> np.ndarray:
        return self._silt_depth

    @property
    def is_volcano(self) -> np.ndarray:
        return self._is_volcano

    @property
    def volcano_active_years_remaining(self) -> np.ndarray:
        return self._volcano_active_years_remaining

    @property
    def soil_depth(self) -> np.ndarray:
        return self._soil_depth

    @property
    def soil_mineral_content(self) -> np.ndarray:
        return self._soil_mineral_content

    @property
    def soil_organic_content(self) -> np.ndarray:
        return self._soil_organic_content

    @property
    def coal_deposit_m(self) -> np.ndarray:
        return self._coal_deposit_m

    @property
    def oil_gas_deposit_m(self) -> np.ndarray:
        return self._oil_gas_deposit_m

    @property
    def mineral_deposit_m(self) -> np.ndarray:
        return self._mineral_deposit_m

    def world_xyz(self, frame: np.ndarray) -> np.ndarray:
        phi_arr = np.full_like(self.theta, self.phi)
        local = geometry.local_xyz(phi_arr, self.theta)
        return geometry.to_world(frame, local)

    def __len__(self) -> int:
        return len(self._theta)

    def __iter__(self) -> Iterator["ElevationPointOnLine"]:
        for i in range(len(self._theta)):
            yield ElevationPointOnLine(self, i)

    def __getitem__(self, index: int) -> "ElevationPointOnLine":
        return ElevationPointOnLine(self, index)

    def replace(self, **overrides: np.ndarray) -> "ElevationLine":
        """A new line with the given fields (elevation and/or any of OPTIONAL_FIELDS)
        swapped in and every other field copied from this one unchanged -- `theta`/`phi`
        are never touched here, so only use this when the node set itself isn't changing."""
        kwargs: dict[str, np.ndarray] = {name: getattr(self, name) for name in self.OPTIONAL_FIELDS}
        kwargs["elevation"] = self.elevation
        kwargs.update(overrides)
        return ElevationLine(phi=self.phi, theta=self.theta, **kwargs)

    def masked(self, mask) -> "ElevationLine":
        """A new line with `theta`, `elevation`, and every OPTIONAL_FIELDS array filtered
        and/or reordered together by a boolean mask or fancy index -- for removing nodes
        (plate split, node reassignment) or reordering them (after concatenating in new
        nodes at the end)."""
        kwargs = {name: getattr(self, name)[mask] for name in self.OPTIONAL_FIELDS}
        return ElevationLine(phi=self.phi, theta=self.theta[mask], elevation=self.elevation[mask], **kwargs)

    def set_fields(self, **fields: np.ndarray) -> None:
        """In-place bulk write for `elevation` and/or any `OPTIONAL_FIELDS` name, straight
        into this line's own backing arrays -- unlike `replace`, no new `ElevationLine` comes
        back. For a caller (`Plate.set_fields_on_plate`) writing values already aligned 1:1
        with this line's existing node order (theta's shape/order unchanged), this is the
        vectorized counterpart to looping `ElevationPointOnLine.set_*` one node at a time."""
        for name, values in fields.items():
            getattr(self, f"_{name}")[:] = values

    def with_new_nodes(self, theta: np.ndarray, elevation: np.ndarray) -> "ElevationLine":
        """A new line with `theta`/`elevation` nodes appended at the end -- every
        OPTIONAL_FIELDS value for the new nodes starts at zero/False, no history to carry.
        The result is unsorted by theta; follow with `.masked(np.argsort(new_line.theta))`
        if ascending order matters to the caller."""
        n = len(theta)
        kwargs = {
            name: np.concatenate([getattr(self, name), np.zeros(n, dtype=getattr(self, name).dtype)])
            for name in self.OPTIONAL_FIELDS
        }
        return ElevationLine(
            phi=self.phi,
            theta=np.concatenate([self.theta, theta]),
            elevation=np.concatenate([self.elevation, elevation]),
            **kwargs,
        )


class ElevationPoint(Protocol):
    """A single node's worth of `ElevationLine` data -- structural (not a base class), so any
    representation-specific backing (`ElevationPointOnLine`, plates.py's own
    `ElevationPointInCloud` for `PlateWithRTree`) can satisfy it without sharing a base.

    `phi`/`get_theta()` are this point's fixed position, get-only: no code in this simulation
    ever moves a single node in place -- position changes always go through a whole-line or
    whole-plate rebuild (`ElevationLine.replace`/`masked`/`with_new_nodes`, boundary growth,
    `regularize_line`, `PlateWithRTree.set_nodes`), so a per-point position setter would just
    invite a caller to silently desync a line's ordering or an R-tree's index instead of going
    through one of those. `elevation` and every `ElevationLine.OPTIONAL_FIELDS` name get real
    setters, since per-node *value* mutation (an eroded elevation, a grown channel, a newly lit
    volcano) is exactly what per-step simulation passes do."""

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
    """Class decorator attaching `get_theta` plus a `get_<name>`/`set_<name>` pair for
    `elevation` and every `ElevationLine.OPTIONAL_FIELDS` name to `cls`, which need only
    provide `_field_array(self, name) -> np.ndarray` and an `_index` attribute -- shared by
    `ElevationPointOnLine` below and plates.py's `ElevationPointInCloud`, so both
    `ElevationPoint` implementations stay wired to the same one field list `ElevationLine`'s
    own `replace`/`masked`/`with_new_nodes` are keyed off, rather than each hand-writing (and
    risking silently forgetting) a method per field -- see `ElevationLine`'s own docstring for
    the bug class that's avoided by never doing this field-by-field by hand."""
    setattr(cls, "get_theta", _point_field_getter("theta"))
    for _name in ("elevation",) + ElevationLine.OPTIONAL_FIELDS:
        setattr(cls, f"get_{_name}", _point_field_getter(_name))
        setattr(cls, f"set_{_name}", _point_field_setter(_name))
    return cls


@install_point_field_accessors
class ElevationPointOnLine:
    """One node of an `ElevationLine`: a pointer to the line plus its index within it. A live
    view, not a snapshot -- `get_*` reads the line's own arrays and `set_*` mutates them in
    place at `index`, so a point handed out by iterating a line (or plate) stays valid and
    stays wired to that same underlying data for as long as the line's node set itself doesn't
    change shape (a `replace`/`masked`/`with_new_nodes` call, or `regularize_line`, produces a
    *new* `ElevationLine` -- any point held from before that call is now stale, the same way a
    Python list index would be after the list it was taken from got reassigned elsewhere)."""

    def __init__(self, line: ElevationLine, index: int) -> None:
        n = len(line)
        if not -n <= index < n:
            raise IndexError(f"ElevationLine point index {index} out of range for length {n}")
        self._line = line
        self._index = index % n

    @property
    def line(self) -> ElevationLine:
        return self._line

    @property
    def index(self) -> int:
        """This point's (always non-negative) position within `line`."""
        return self._index

    @property
    def phi(self) -> float:
        return self._line.phi

    def _field_array(self, name: str) -> np.ndarray:
        return getattr(self._line, f"_{name}")


def iter_local_lattice(frame: np.ndarray, spacing_rad: float = TARGET_LINE_SPACING_RAD):
    """Sweep a full plate-local (phi, theta) lattice at `spacing_rad` resolution, yielding
    (phi, theta_candidates, world_pts) per row. Shared by initial generation and by
    plate-merge resampling (see merge_split.py), and, at a resolution independent of the
    physical line spacing, by the render-grid sweep (see render_image.py's
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


def build_lines_from_lattice(frame: np.ndarray, is_owned, elevation_at, spacing_rad: float = TARGET_LINE_SPACING_RAD) -> list[ElevationLine]:
    """Build a plate's elevation lines by sweeping its local lattice and keeping whichever
    nodes `is_owned(world_pts) -> bool array` selects, with elevation from
    `elevation_at(owned_world_pts) -> array`. `spacing_rad` defaults to the reference
    density (1.0) -- every caller that has a `World` in hand should instead pass
    `line_spacing_rad(world.node_density)`, so newly-built lines (initial generation, gap
    absorption/spawning, plate merges, volcanic fields) match whatever density that world was
    actually generated at, not silently fall back to the default."""
    lines: list[ElevationLine] = []
    for phi, theta_candidates, world_pts in iter_local_lattice(frame, spacing_rad=spacing_rad):
        owned = is_owned(world_pts)
        if not np.any(owned):
            continue
        theta_owned = theta_candidates[owned]
        elevation = elevation_at(world_pts[owned])
        lines.append(ElevationLine(phi=phi, theta=theta_owned, elevation=elevation))
    return lines


def needs_regularizing(line: ElevationLine, spacing_rad: float = TARGET_LINE_SPACING_RAD) -> bool:
    if len(line) < 3:
        return False
    dtheta_target = spacing_rad / max(np.cos(line.phi), 1e-3)
    gaps = np.diff(line.theta)
    ratio = gaps / dtheta_target
    return bool(np.any(ratio > IRREGULARITY_TOLERANCE) or np.any(ratio < 1.0 / IRREGULARITY_TOLERANCE))


def _crumple_elevation(elevation: np.ndarray, m: int, hurst: float = CRUMPLE_HURST_EXPONENT) -> np.ndarray:
    """Replace n points' worth of elevation with m < n points' worth by "crumpling": fit a
    smooth curve e = f(x) to the n original points (x = 0..n-1, plain sample index -- the
    fit doesn't need to know about theta/phi, just the shape), then read m new values off a
    horizontally squashed version of that same curve, e' = f(x / k), where k = m/n < 1 is how
    squashed the m points are relative to the n they replace. Dividing by k (rather than
    multiplying) is what makes k < 1 actually squash the domain: as the m new sample indices
    range over [0, m-1], x/k ranges over [0, (m-1)/k] = [0, n-1] -- i.e. the same few new
    points now have to cover the *entire* original curve's span, packing all of its shape
    into fewer samples, exactly like real crumpling packs the same strip of material into
    less room.

    Squashing alone (no amplitude change) would keep every new sample within the original
    curve's min/max -- steeper-looking between points, but never actually taller. Real
    crumpled terrain isn't just steeper, it's taller: compressing a self-affine profile
    horizontally by k grows its vertical amplitude by k**-hurst (see CRUMPLE_HURST_EXPONENT),
    so peaks get pushed higher and valleys pulled lower in proportion to how aggressively
    this call is squashing, not by some unrelated fixed multiplier.

    The fit itself is a truncated cosine series (a real, non-periodic basis -- unlike a raw
    FFT, it has no wraparound artifact at the two ends of what is an open curve, never a
    periodic one) with only m+1 terms, not n -- deliberately under-resolved relative to the n
    input points, so the fit is a smoothing regression through them rather than an exact
    interpolation. That smoothing is what discards the sub-target-spacing detail crumpling is
    supposed to be discarding in the first place; fitting all n harmonics would just
    reconstruct every original point exactly and defeat the point of thinning them out.
    """
    n = len(elevation)
    x = np.arange(n, dtype=float)
    num_harmonics = min(n - 1, max(2, m))
    denom = max(n - 1, 1)
    basis = np.stack([np.cos(np.pi * p * x / denom) for p in range(num_harmonics + 1)], axis=1)
    coeffs, *_ = np.linalg.lstsq(basis, elevation, rcond=None)

    k = m / n
    x_new = np.clip(np.arange(m, dtype=float) / k, 0.0, n - 1)
    new_basis = np.stack([np.cos(np.pi * p * x_new / denom) for p in range(num_harmonics + 1)], axis=1)
    fitted = new_basis @ coeffs

    amplitude = k**-hurst
    mean_e = elevation.mean()
    crumpled = mean_e + amplitude * (fitted - mean_e)
    # amplitude > 1 (k < 1) means crumpling can push a peak/valley past what the original n
    # points ever reached -- clamp back into the world's elevation bounds the same way every
    # other module that modifies elevation does (boundary.py, bathymetry.py, erosion.py,
    # volcanism.py), since nothing downstream of regularizing re-checks this.
    crumpled = np.clip(crumpled, MIN_ELEVATION_M, MAX_ELEVATION_M)
    # The fit is a smoothing regression, not an exact interpolant, so it can drift slightly
    # from the original data even at x=0/x=n-1 -- force the two ends back to the real
    # original values so a crumpled line still butts up exactly against its neighbors'
    # elevation at the endpoints regularize_line preserves the position of.
    crumpled[0] = elevation[0]
    crumpled[-1] = elevation[-1]
    return crumpled


def regularize_line(line: ElevationLine, spacing_rad: float = TARGET_LINE_SPACING_RAD) -> ElevationLine:
    if len(line) < 3:
        return line

    dtheta_target = spacing_rad / max(np.cos(line.phi), 1e-3)
    theta_min, theta_max = line.theta[0], line.theta[-1]
    span = theta_max - theta_min
    n = max(int(round(span / dtheta_target)) + 1, 2)

    new_theta = np.linspace(theta_min, theta_max, n)
    # Fewer new nodes than the line already has -- vulcanism-driven density increases (fresh
    # volcano nodes inserted mid-line, see volcanism.py) can push points closer together than
    # target spacing without ever widening a gap, so this is the "too close" direction
    # needs_regularizing also fires on. Crumple instead of linearly resampling here: a plain
    # np.interp thin-out can smooth away or altogether skip a narrow peak that happens to fall
    # between two kept sample points, where crumpling fits the whole n-point shape first and
    # only then reads fewer values off it, so a peak influences every new sample near it
    # rather than being invisible to all but its two immediate neighbors.
    if n < len(line):
        new_elevation = _crumple_elevation(line.elevation, n)
    else:
        new_elevation = np.interp(new_theta, line.theta, line.elevation)
    # channel_depth/channel_width/lake_depth/glacier_depth interpolated the same way -- a
    # plain reset to 0 here would wipe out a river's carved channel (or a glacier) every time
    # this line's spacing drifts enough to trigger regularizing, which runs periodically
    # throughout the simulation (see REGULARIZE_INTERVAL_STEPS), not as a rare one-off event
    # like a merge/split resample.
    new_channel_depth = np.interp(new_theta, line.theta, line.channel_depth)
    new_channel_width = np.interp(new_theta, line.theta, line.channel_width)
    new_lake_depth = np.interp(new_theta, line.theta, line.lake_depth)
    new_glacier_depth = np.interp(new_theta, line.theta, line.glacier_depth)
    new_silt_depth = np.interp(new_theta, line.theta, line.silt_depth)
    # volcano_active_years_remaining interpolates the same way; is_volcano is interpolated as
    # a float (blending a volcano node's 1.0 against a non-volcano neighbor's 0.0) then
    # thresholded back to bool, same spirit as the others -- a resampled node keeps "was this
    # near a volcano" rather than silently losing volcanic provenance every regularize pass.
    new_volcano_active_years_remaining = np.interp(new_theta, line.theta, line.volcano_active_years_remaining)
    new_is_volcano = np.interp(new_theta, line.theta, line.is_volcano.astype(float)) > 0.5
    # Soil/resource fields (see geology.py/volcanism.py) interpolated the same way as the
    # rest -- a plain reset to 0 here would wipe out accumulated soil/coal/oil-gas/mineral
    # deposits every time a line's spacing drifts enough to trigger regularizing.
    new_soil_depth = np.interp(new_theta, line.theta, line.soil_depth)
    new_soil_mineral_content = np.interp(new_theta, line.theta, line.soil_mineral_content)
    new_soil_organic_content = np.interp(new_theta, line.theta, line.soil_organic_content)
    new_coal_deposit_m = np.interp(new_theta, line.theta, line.coal_deposit_m)
    new_oil_gas_deposit_m = np.interp(new_theta, line.theta, line.oil_gas_deposit_m)
    new_mineral_deposit_m = np.interp(new_theta, line.theta, line.mineral_deposit_m)
    return ElevationLine(
        phi=line.phi,
        theta=new_theta,
        elevation=new_elevation,
        channel_depth=new_channel_depth,
        channel_width=new_channel_width,
        lake_depth=new_lake_depth,
        glacier_depth=new_glacier_depth,
        silt_depth=new_silt_depth,
        is_volcano=new_is_volcano,
        volcano_active_years_remaining=new_volcano_active_years_remaining,
        soil_depth=new_soil_depth,
        soil_mineral_content=new_soil_mineral_content,
        soil_organic_content=new_soil_organic_content,
        coal_deposit_m=new_coal_deposit_m,
        oil_gas_deposit_m=new_oil_gas_deposit_m,
        mineral_deposit_m=new_mineral_deposit_m,
    )


def regularize_plate_lines(plate: "PlateWithLines", spacing_rad: float = TARGET_LINE_SPACING_RAD) -> None:
    plate.set_lines(
        [regularize_line(line, spacing_rad) if needs_regularizing(line, spacing_rad) else line for line in plate.lines]
    )


def regularize_world_lines(world: "World") -> None:
    spacing_rad = line_spacing_rad(world.node_density)
    for plate in world.plates:
        regularize_plate_lines(plate, spacing_rad)

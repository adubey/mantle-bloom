"""Intraplate fault lines -- a first-class tectonic feature that is *not* a plate boundary.

Plate boundaries carry their own deformation in `PlateWithLines.deform` (plates.py):
classification there is geometric (contested territory -> convergent, uncontested-but-near
-> transform, wider -> divergent). That model has no notion of a fault line sitting *inside*
a plate, away from any edge -- yet in reality faults nucleate at a wide range of distances
from plate boundaries (most within ~200 km, but stable-continental-region faults sit well
over 1000 km away -- the New Madrid seismic zone is ~1500 km from the nearest boundary),
come in sub-parallel families (Basin-and-Range horst/graben trains, en echelon step-overs),
and stay individually active for a few to a few tens of Myr before locking up and surviving
as inert scars.

By default this module is an **additive** layer: it never touches the live deform()
classification (`LithospherePlate.deform`). Each step it

1. ages every existing fault, accumulating slip on the active ones and retiring those past
   their drawn lifespan (kept forever after as an inactive scar, like `is_volcano`);
2. rolls a stress-weighted Poisson spawn per plate -- probability high in a band near the
   boundary, decaying exponentially into the interior with a small nonzero floor everywhere,
   **and lifted to `OVERLAP_STRESS_WEIGHT` wherever this plate's nodes sit on top of a
   neighbour's** (a point overlap -- a stalled collision, a plate drifting bodily over one it
   can't merge with -- is locally super-stressed crust), the regime (normal / reverse /
   strike-slip) picked from the local closing rate per Andersonian faulting theory, and with
   `SET_PROBABILITY` chance a whole sub-parallel family rather than a lone trace;
3. applies each active fault's own relief to the nearby crust -- reverse: an uplift ridge;
   normal: a hanging-wall graben with a footwall shoulder; strike-slip: a modest
   transpressional ridge or transtensional sag (relief only -- the node field is *not*
   physically sheared across the trace, see docs/TODO.md);
4. rolls each active fault's **earthquakes** for the step (`_generate_earthquakes`): a
   Poisson count from `slip_rate * dt / CHARACTERISTIC_SLIP_PER_QUAKE_M`, each a transient
   located `Earthquake` (magnitude from trace length + slip rate) appended to
   `World.earthquakes`, pruned after `EARTHQUAKE_RETAIN_MYR`. `erosion.py` reads them for a
   local seismic-erosion burst; the "Fault lines" view draws them as a fading overlay.

When `World.fault_deformation_mode` is `"fault"` or `"both"` the layer stops being purely
additive: `_apply_plate_fault_relief`'s rates/reach scale up (`FAULT_RELIEF_MODE_*`) and, in
`"fault"` mode, `LithospherePlate.deform` gates its own boundary thickening by
`fault_influence()` so plate-boundary transformation localises onto fault lines rather than a
smooth band at the polygon edge. `"boundary"` (the default) is bit-identical to before.

Geometry is stored in the owning plate's **local frame** (`local_phi` / `local_theta`), so
a fault rides along with the crust as the plate rotates for free -- the same "attached to
the crust, not the world" property every persistent `ElevationLine` field already has.
`reconcile_faults` re-homes faults across merges/splits and drops those whose plate
subducted (see world.step_world).

**Fault systems** (`FaultSystem`) sit one level above the individual trace: with
`SYSTEM_SPAWN_FRACTION` of spawns, instead of a lone fault / tight set we lay a long,
gently curving *master lineament* (up to ~5500 km -- East African Rift / Anatolian /
Sunda scale) and scatter a wide sub-parallel family of strands along its belt, each strand
an ordinary `Fault` (with a widened length distribution, tail to ~1300 km) carrying the
system's `system_id`. The master trace applies no relief of its own -- it is an organising
scaffold; the strands carry the relief exactly as lone faults do.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import boundary, geometry
from .elevation_lines import (
    ELEV_CHANGE_FAULT_NORMAL,
    ELEV_CHANGE_FAULT_REVERSE,
    ELEV_CHANGE_FAULT_STRIKE_SLIP,
    ELEV_CHANGE_MIN_DELTA_M,
    MAX_ELEVATION_M,
    MIN_ELEVATION_M,
    PLANET_RADIUS_KM,
    line_spacing_rad,
)
from .plates import OVERLAP_TOLERANCE_MULT, Plate, collect_all_points, query_workers

if TYPE_CHECKING:
    from .world import World

# RNG stream tag, keyed alongside (world.seed, round(world.elapsed_years), plate_id) -- the
# same "deterministic per (seed, elapsed_years, plate)" convention volcanism.py uses so a
# replayed session spawns the identical faults. Arbitrary int, distinct from every other tag.
_FAULT_SEED_TAG = 7331

# --- Spawn model (see the real-world table in docs/simulation-model.md#faults) ---
# Expected new fault *systems* per Myr over the whole sphere at full boundary stress; the
# per-plate roll scales this by the plate's area fraction and its mean stress weight, so the
# realised rate is far lower (most of a plate is low-stress interior). Tunable -- raise it to
# see faults accumulate faster.
BASE_SPAWN_RATE_PER_MYR = 3.0
# Stress weight vs distance from the nearest cross-plate boundary: exp(-d / DECAY_LEN) with a
# floor, so faulting concentrates near boundaries but never drops to zero in the deep
# interior (real intraplate seismicity).
SPAWN_DECAY_LEN_KM = 500.0
SPAWN_INTERIOR_FLOOR = 0.03

# Segment length: lognormal, a few km to ~200 km. NOTE: this (and the fault-set spread
# below) is known to be far too short -- real faults run to ~1300 km and fault *systems* to
# ~5500 km. See docs/TODO.md "Intraplate faults: follow-ups" item 1.
LENGTH_MEDIAN_KM = 45.0
LENGTH_SIGMA = 0.6
LENGTH_MIN_KM = 12.0
LENGTH_MAX_KM = 200.0
# Nodes sampled along a fault trace -- enough to render a gently curved line, capped so a
# long fault doesn't dominate the relief query.
FAULT_NODES_MIN = 4
FAULT_NODES_MAX = 14
# Along-strike curvature: a fraction of the fault length, perpendicular to strike -- gives
# strike-slip faults the restraining/releasing character their relief keys off, and keeps
# fault sets from looking like a ruled grid.
BEND_MAX_FRACTION = 0.12

# Slip rate: intraplate faults ~0.1-5 mm/yr (100-5000 m/Myr); a major boundary-adjacent
# transform up to ~35 mm/yr. Scaled between these by the seed node's stress weight.
SLIP_RATE_MIN_M_PER_MYR = 150.0
SLIP_RATE_MAX_M_PER_MYR = 30000.0
SLIP_RATE_REF_M_PER_MYR = 3000.0  # relief magnitude is normalised against this

# Active lifespan (Myr), then a permanent inactive scar. Higher-slip faults live longer.
LIFESPAN_MIN_MYR = 2.0
LIFESPAN_MAX_MYR = 25.0

# Andersonian dips.
DIP_NORMAL_DEG = 60.0
DIP_REVERSE_DEG = 30.0
DIP_STRIKE_SLIP_DEG = 90.0
# Strike-slip faults strike obliquely to the local shortening direction.
STRIKE_SLIP_OBLIQUITY_DEG = 30.0

# Fault sets -- tight sub-parallel families born of a single lone spawn. Basin-and-Range
# major normal faults sit ~15-30 km apart; en echelon step-overs ~1-5 km. (For the larger,
# spatially-extended sub-parallel family that shares one belt over 1000s of km, see fault
# systems below -- a set is the local cluster, a system is the whole zone.)
SET_PROBABILITY = 0.4
SET_MIN_MEMBERS = 2
SET_MAX_MEMBERS = 5
SET_SPACING_KM = 20.0
SET_ECHELON_STEP_KM = 3.0

# --- Fault systems: a first-class structure one level above the individual trace ---
# A real fault zone / system is many sub-parallel and en echelon strands acting together
# along one gently curving belt: the East African Rift, the Anatolian system and the
# Sunda/Sumatran system all reach ~5500 km, and individual continuous strands within them
# run to ~1300 km (San Andreas ~1200, North Anatolian ~1500, Great Sumatran ~1900). The
# lone-fault/set path above tops out near ~200 km, so systems are spawned as their own
# thing: with SYSTEM_SPAWN_FRACTION of Poisson spawns, instead of a lone trace/set we lay a
# long master lineament and scatter a wide strand family along it. Everything else about a
# strand (relief, aging, scars, reconcile) is identical to a lone fault -- a system is an
# organising scaffold, not a new relief mechanism.
SYSTEM_SPAWN_FRACTION = 0.18

# Master lineament length -- lognormal, long right tail, no hard clamp near the median.
SYSTEM_LENGTH_MEDIAN_KM = 2200.0
SYSTEM_LENGTH_SIGMA = 0.5
SYSTEM_LENGTH_MIN_KM = 600.0
SYSTEM_LENGTH_MAX_KM = 5500.0
# The belt curves: a few low-frequency lobes of lateral wander, as a fraction of length.
SYSTEM_BEND_LOBES = 2
SYSTEM_BEND_MAX_FRACTION = 0.06
SYSTEM_MASTER_NODE_KM = 70.0
SYSTEM_MASTER_NODES_MAX = 80

# Strand family: sub-parallel traces scattered along and across the belt.
SYSTEM_STRAND_SPACING_KM = 65.0  # mean along-belt gap between strand seed points
SYSTEM_BELT_HALF_WIDTH_KM = 130.0  # strands scatter this far either side of the master trace
SYSTEM_STRAND_COUNT_MIN = 5
SYSTEM_STRAND_COUNT_MAX = 16
SYSTEM_STRAND_STRIKE_JITTER_DEG = 12.0
SYSTEM_OFFREGIME_FRACTION = 0.15  # strands that buck the system's dominant regime
# Strand length -- widened vs the lone-fault LENGTH_* (median ~45, max 200): a system strand
# is a major fault in its own right.
SYSTEM_STRAND_LENGTH_MEDIAN_KM = 150.0
SYSTEM_STRAND_LENGTH_SIGMA = 0.7
SYSTEM_STRAND_LENGTH_MIN_KM = 20.0
SYSTEM_STRAND_LENGTH_MAX_KM = 1300.0
SYSTEM_STRAND_NODE_KM = 35.0
SYSTEM_STRAND_NODES_MAX = 40

# A system outlives its individual strands (2-25 Myr) by an order of magnitude -- the belt
# stays a locus of faulting long after any one strand locks up.
SYSTEM_LIFESPAN_MIN_MYR = 25.0
SYSTEM_LIFESPAN_MAX_MYR = 140.0
MAX_INACTIVE_SYSTEMS_PER_PLATE = 12

_FAULT_SYSTEM_SEED_TAG = 7332

# Relief (per Myr, at the trace, tapering linearly to zero at MAX_FAULT_REACH_KM). Kept well
# below the boundary rates in plates.py (CONVERGENT_MOUNTAIN_RATE_M_PER_MYR = 800) so this
# additive layer doesn't disturb long-run hypsometry tuning.
MAX_FAULT_REACH_KM = 45.0
REVERSE_UPLIFT_M_PER_MYR = 220.0
NORMAL_THROW_M_PER_MYR = 180.0  # hanging-wall down
NORMAL_SHOULDER_UPLIFT_M_PER_MYR = 55.0  # footwall up
STRIKE_SLIP_RIDGE_M_PER_MYR = 70.0  # small always-transpressional component
STRIKE_SLIP_BEND_M_PER_MYR = 130.0  # +restraining (uplift) / -releasing (sag), per strike_sense

# Bound total memory: a plate keeps at most this many inactive scars (oldest culled first);
# active faults are never culled.
MAX_SCARS_PER_PLATE = 40

# --- Point-overlap spawning (see _maybe_spawn_faults / World docstring) ---
# A node sitting within this many line-spacings of another plate's node is a genuine point
# overlap (two plates tiling the same patch of sphere), not a normal shared boundary -- the
# same OVERLAP_TOLERANCE_MULT the Plate Inspector / merge_split.update_overlap_tracking use.
# Such a node is treated as locally super-stressed: its spawn weight is lifted to
# OVERLAP_STRESS_WEIGHT (> 1.0, i.e. more stressed than a clean convergent edge), which both
# raises the plate's mean weight (more fault systems spawn) and pulls seeds into the overlap.
OVERLAP_STRESS_WEIGHT = 1.5

# --- Fault-deformation mode (World.fault_deformation_mode; see LithospherePlate.deform) ---
FAULT_DEFORMATION_MODES = ("boundary", "fault", "both")
# In "fault" mode, boundary thickening in LithospherePlate.deform is multiplied by
# fault_influence(): 1.0 within FAULT_DEFORM_REACH_KM of an active fault trace, tapering to
# FAULT_DEFORM_FLOOR far from one (never 0 -- a contested zone with no fault yet still
# deforms while Piece-1 spawning fills it in).
FAULT_DEFORM_REACH_KM = 120.0
FAULT_DEFORM_FLOOR = 0.15
# In "fault"/"both" mode faults.py's own relief layer is scaled up to carry the deformation
# the smooth boundary bands give up -- rates toward plates.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR
# (800), reach wider. Exactly 1.0 in "boundary" mode (bit-identical to before this existed).
FAULT_RELIEF_MODE_RATE_SCALE = 3.0
FAULT_RELIEF_MODE_REACH_SCALE = 1.6

# --- Earthquakes (see Earthquake / _generate_earthquakes) ---
_QUAKE_SEED_TAG = 7333
# A fault that accumulated at least this much slip this step gets one `Earthquake` record --
# the characteristic (roughly largest) rupture of that interval. Over a Myr a real active
# fault ruptures thousands of times; storing each is pointless, so one representative event
# per active fault per step is the model. Below this the fault is treated as aseismic for the
# step (a barely-creeping trace).
MIN_STEP_SLIP_FOR_QUAKE_M = 25.0
# A fault born in a point-overlap zone (a stalled collision -- a seismic hotspot) adds this to
# its characteristic magnitude.
OVERLAP_QUAKE_MW_BONUS = 0.4
# Moment magnitude from trace length and slip rate: Mw ~ base + len_coeff*log10(km) +
# slip_coeff*log10(m/Myr), clamped, plus a small per-event jitter. Eyeballed to put a
# ~200 km / 3000 m-per-Myr fault near Mw 6.5 and a ~1300 km system strand near Mw 8.
QUAKE_MW_BASE = 3.4
QUAKE_MW_LENGTH_COEFF = 1.15
QUAKE_MW_SLIP_COEFF = 0.35
QUAKE_MW_JITTER = 0.3
QUAKE_MW_MIN = 4.0
QUAKE_MW_MAX = 9.3
# At most one earthquake per step is written to the UI event console -- the step's largest,
# and only if it clears this magnitude (a genuinely major event). The console is for
# structural/topology history; routine seismicity lives on the "Fault lines" overlay.
EARTHQUAKE_LOG_MIN_MW = 7.5
# An earthquake is dropped from World.earthquakes once older than this -- bounds memory and
# matches the fading-overlay window in the "Fault lines" view.
EARTHQUAKE_RETAIN_MYR = 5.0

_KIND_NORMAL = "normal"
_KIND_REVERSE = "reverse"
_KIND_STRIKE_SLIP = "strike_slip"
_KIND_LABEL = {_KIND_NORMAL: "normal", _KIND_REVERSE: "reverse", _KIND_STRIKE_SLIP: "strike-slip"}
_KIND_REASON = {
    _KIND_NORMAL: ELEV_CHANGE_FAULT_NORMAL,
    _KIND_REVERSE: ELEV_CHANGE_FAULT_REVERSE,
    _KIND_STRIKE_SLIP: ELEV_CHANGE_FAULT_STRIKE_SLIP,
}


@dataclass
class Fault:
    """One fault trace, geometry in the owning plate's local frame (see module docstring).

    `dip_dir_local` is a plate-local unit-ish vector in the tangent plane at the fault
    midpoint pointing toward the hanging-wall / downthrown side (normal/reverse only;
    unused for strike-slip). `strike_sense` is +1 (restraining / transpressional) or -1
    (releasing / transtensional) for strike-slip; for a normal fault it flips which side
    of a set member is downthrown.
    """

    fault_id: int
    plate_id: int
    kind: str
    local_phi: np.ndarray
    local_theta: np.ndarray
    slip_rate_m_per_myr: float
    dip_deg: float
    strike_sense: int
    dip_dir_local: np.ndarray
    lifespan_myr: float
    birth_years: float
    birth_distance_from_boundary_km: float
    set_id: int | None = None
    # The fault system (see FaultSystem) this trace is a strand of, or None for a lone
    # fault / tight set. A plain-default field, so an older pickle without it reads None.
    system_id: int | None = None
    # True if the seed node lay in a point overlap with another plate (see
    # OVERLAP_STRESS_WEIGHT) -- a plain-default field an older pickle reads as False. Drives
    # a higher rupture rate in _generate_earthquakes.
    born_in_overlap: bool = False
    age_myr: float = 0.0
    cumulative_offset_m: float = 0.0
    active: bool = True
    # Refreshed at the end of update_faults every step (world-space polyline, true frame) so
    # reconcile_faults can re-home a fault even after its plate was absorbed by a merge.
    # Not authoritative geometry -- local_phi/local_theta are.
    world_polyline: np.ndarray | None = field(default=None, repr=False, compare=False)

    def length_km(self) -> float:
        return _polyline_length_km(self.local_phi, self.local_theta)


@dataclass
class FaultSystem:
    """A fault zone / system: one long, gently curving master lineament plus the family of
    sub-parallel strands (`Fault`s carrying this `system_id`) scattered along its belt. The
    master trace is an organising scaffold -- it applies no relief of its own; each strand
    does, exactly as a lone fault would. Geometry is stored in the owning plate's local
    frame, same as a `Fault`, so the whole belt rides with the crust.

    A system outlives its strands: it stays `active` (a locus that keeps spawning fresh
    strands is a future refinement -- for now `active` just drives rendering and culling)
    for `lifespan_myr`, then becomes an inert scar bundle like a locked-up fault.
    """

    system_id: int
    plate_id: int
    kind: str  # the belt's dominant regime
    master_local_phi: np.ndarray
    master_local_theta: np.ndarray
    length_km: float
    birth_years: float
    lifespan_myr: float
    age_myr: float = 0.0
    active: bool = True
    world_polyline: np.ndarray | None = field(default=None, repr=False, compare=False)

    def master_length_km(self) -> float:
        return _polyline_length_km(self.master_local_phi, self.master_local_theta)


@dataclass
class Earthquake:
    """One rupture on an active fault this-or-a-recent step. `epicenter_world` is a fixed
    unit vector in the true (un-rotated) world frame -- an earthquake is an event at a place,
    not a persistent crustal feature, so it is never re-homed across a merge/split; it just
    ages out of `World.earthquakes` after `EARTHQUAKE_RETAIN_MYR`. `magnitude` is a moment
    magnitude (see `_earthquake_magnitude`)."""

    earthquake_id: int
    fault_id: int
    plate_id: int
    kind: str
    epicenter_world: np.ndarray
    magnitude: float
    slip_m: float
    birth_years: float


def _polyline_length_km(local_phi: np.ndarray, local_theta: np.ndarray) -> float:
    pts = geometry.local_xyz(local_phi, local_theta)
    seg = np.arccos(np.clip(np.sum(pts[:-1] * pts[1:], axis=-1), -1.0, 1.0))
    return float(np.sum(seg) * PLANET_RADIUS_KM)


def plate_by_id(world: "World") -> dict[int, Plate]:
    return {p.plate_id: p for p in world.plates}


def fault_world_points(fault: Fault, plate: Plate) -> np.ndarray:
    """The fault's trace in true world coordinates, via the owning plate's current frame."""
    return geometry.to_world(plate.frame, geometry.local_xyz(fault.local_phi, fault.local_theta))


def system_world_points(system: FaultSystem, plate: Plate) -> np.ndarray:
    """The system's master lineament in true world coordinates."""
    return geometry.to_world(plate.frame, geometry.local_xyz(system.master_local_phi, system.master_local_theta))


# --------------------------------------------------------------------------- step entry point


def update_faults(world: "World", years: float) -> None:
    """Age / spawn / retire faults, apply their relief, and roll each active fault's
    earthquakes. Called once per step from world.step_world, inside the
    simulate_plate_movement block, right after the deform loop and before
    merge_split.apply_topology_changes."""
    years_myr = years / 1_000_000.0

    # Drop earthquakes that have aged out of the retention window (memory bound + the
    # fading-overlay window). Done first so a rupture logged this step isn't immediately
    # eligible for pruning.
    if world.earthquakes:
        cutoff = world.elapsed_years - EARTHQUAKE_RETAIN_MYR * 1_000_000.0
        world.earthquakes = [q for q in world.earthquakes if q.birth_years >= cutoff]

    for fault in world.faults:
        fault.age_myr += years_myr
        if fault.active:
            fault.cumulative_offset_m += fault.slip_rate_m_per_myr * years_myr
            if fault.age_myr >= fault.lifespan_myr:
                fault.active = False
                world.log_event(
                    f"fault #{fault.fault_id} ({_KIND_LABEL[fault.kind]}) on plate {fault.plate_id} "
                    f"locked up after {fault.age_myr:.0f} Myr"
                )

    for system in world.fault_systems:
        system.age_myr += years_myr
        if system.active and system.age_myr >= system.lifespan_myr:
            system.active = False
            world.log_event(
                f"fault system #{system.system_id} ({_KIND_LABEL[system.kind]}, ~{system.length_km:.0f} km) "
                f"on plate {system.plate_id} went inactive after {system.age_myr:.0f} Myr"
            )

    total_nodes = max(1, sum(p.node_count() for p in world.plates))
    for plate in world.plates:
        _maybe_spawn_faults(world, plate, total_nodes, years_myr)

    _cull_scars(world)
    _cull_inactive_systems(world)

    for plate in world.plates:
        _apply_plate_fault_relief(world, plate, years_myr)

    _generate_earthquakes(world, years_myr)

    by_id = plate_by_id(world)
    for fault in world.faults:
        plate = by_id.get(fault.plate_id)
        if plate is not None:
            fault.world_polyline = fault_world_points(fault, plate)
    for system in world.fault_systems:
        plate = by_id.get(system.plate_id)
        if plate is not None:
            system.world_polyline = system_world_points(system, plate)


def reconcile_faults(world: "World") -> None:
    """After a topology change: drop faults whose plate subducted, and re-home faults whose
    trace midpoint now lies in a different surviving plate's territory (covers both a merge
    -- the absorbed plate's id vanishes -- and a split -- a new id appears near the cut).
    Recomputes local coordinates in the new plate's frame; preserves age / offset / active.
    Mirrors stranded_basins.reconcile_world_tracks. Fault systems' master lineaments are
    re-homed by the same midpoint rule."""
    if not world.faults and not world.fault_systems:
        return
    live = plate_by_id(world)
    if not live:
        world.faults = []
        world.fault_systems = []
        return
    collected = collect_all_points(world.plates)
    if collected is None:
        world.faults = []
        world.fault_systems = []
        return
    points, _, owner = collected
    tree = cKDTree(points)

    kept: list[Fault] = []
    for fault in world.faults:
        world_poly = fault.world_polyline
        if world_poly is None:
            plate = live.get(fault.plate_id)
            if plate is None:
                continue  # no geometry to re-home with -- drop it
            world_poly = fault_world_points(fault, plate)
        mid = world_poly[len(world_poly) // 2]
        _, idx = tree.query(mid)
        new_pid = int(owner[idx])
        if new_pid == fault.plate_id and fault.plate_id in live:
            kept.append(fault)
            continue
        new_plate = live.get(new_pid)
        if new_plate is None:
            continue
        local = geometry.to_local(new_plate.frame, world_poly)
        phi, theta = geometry.xyz_to_latlon(local)
        fault.plate_id = new_pid
        fault.local_phi = phi
        fault.local_theta = theta
        fault.dip_dir_local = _rehome_dip_dir(world_poly, new_plate)
        fault.world_polyline = fault_world_points(fault, new_plate)
        kept.append(fault)
    world.faults = kept

    # Re-home the systems' master lineaments the same way (by their midpoint). A strand and
    # its system can land on different plates after a split cuts the belt -- that's fine,
    # each carries its own plate_id; `system_id` just links them for display.
    kept_systems: list[FaultSystem] = []
    for system in world.fault_systems:
        world_poly = system.world_polyline
        if world_poly is None:
            plate = live.get(system.plate_id)
            if plate is None:
                continue
            world_poly = system_world_points(system, plate)
        mid = world_poly[len(world_poly) // 2]
        _, idx = tree.query(mid)
        new_pid = int(owner[idx])
        if new_pid != system.plate_id or system.plate_id not in live:
            new_plate = live.get(new_pid)
            if new_plate is None:
                continue
            local = geometry.to_local(new_plate.frame, world_poly)
            phi, theta = geometry.xyz_to_latlon(local)
            system.plate_id = new_pid
            system.master_local_phi = phi
            system.master_local_theta = theta
            system.world_polyline = system_world_points(system, new_plate)
        kept_systems.append(system)
    world.fault_systems = kept_systems


def _rehome_dip_dir(world_poly: np.ndarray, new_plate: Plate) -> np.ndarray:
    """Rebuild dip_dir in the new plate's local frame: take the fault's own strike from its
    (already re-homed) local polyline and a tangent-plane perpendicular at the midpoint, then
    keep whichever perpendicular sense the old world-space dip direction pointed."""
    mid = world_poly[len(world_poly) // 2]
    east, _ = geometry.local_tangent_basis(mid)
    a = world_poly[0]
    b = world_poly[-1]
    strike = b - a
    strike = strike - np.dot(strike, mid) * mid
    n = np.linalg.norm(strike)
    strike = strike / n if n > 1e-9 else east
    perp = np.cross(mid, strike)
    pn = np.linalg.norm(perp)
    perp = perp / pn if pn > 1e-9 else east
    return geometry.to_local(new_plate.frame, perp)


# --------------------------------------------------------------------------- spawning


def _plate_stress(world: "World", plate: Plate):
    """Per own-node: world position, distance (rad) to the nearest cross-plate boundary
    node, that neighbour's omega, and that neighbour's position. None if this plate has no
    neighbours or no nodes."""
    own_points = plate.all_points_and_elevation()[0]
    if len(own_points) == 0:
        return None
    neighbours = plate.get_neighbours(world.plates)
    if not neighbours:
        return None
    nb_pts = np.concatenate([p.all_points_and_elevation()[0] for p in neighbours], axis=0)
    nb_omega = np.concatenate(
        [np.tile(np.asarray(p.omega, dtype=float), (p.node_count(), 1)) for p in neighbours], axis=0
    )
    if len(nb_pts) == 0:
        return None
    tree = cKDTree(nb_pts, balanced_tree=False, compact_nodes=False)
    dist, idx = tree.query(own_points, workers=query_workers(len(own_points)))
    return own_points, dist, nb_omega[idx], nb_pts[idx]


def _maybe_spawn_faults(world: "World", plate: Plate, total_nodes: int, years_myr: float) -> None:
    stress = _plate_stress(world, plate)
    if stress is None:
        return
    own_points, dist, nn_omega, nn_pts = stress
    rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate.plate_id, _FAULT_SEED_TAG))

    dist_km = dist * PLANET_RADIUS_KM
    weight = np.exp(-dist_km / SPAWN_DECAY_LEN_KM)
    weight = SPAWN_INTERIOR_FLOOR + (1.0 - SPAWN_INTERIOR_FLOOR) * np.clip(weight, 0.0, 1.0)

    # Point overlap -> locally super-stressed crust (see OVERLAP_STRESS_WEIGHT). A node this
    # close to a neighbour's node isn't a normal shared boundary (those sit ~1 spacing apart);
    # it's two plates tiling the same ground -- a stalled collision, a plate drifting bodily
    # over one it can't merge with. Lifting the weight here both raises the plate's mean
    # weight (`expected` below -> more fault systems) and pulls seeds into the overlap.
    overlap_tol_rad = OVERLAP_TOLERANCE_MULT * line_spacing_rad(world.node_density)
    overlap_mask = dist < overlap_tol_rad
    if np.any(overlap_mask):
        weight[overlap_mask] = np.maximum(weight[overlap_mask], OVERLAP_STRESS_WEIGHT)

    area_frac = len(own_points) / total_nodes
    expected = BASE_SPAWN_RATE_PER_MYR * area_frac * float(np.mean(weight)) * years_myr
    n_spawn = int(rng.poisson(max(expected, 0.0)))
    if n_spawn == 0:
        return

    probs = weight / weight.sum()
    seeds = rng.choice(len(own_points), size=n_spawn, p=probs)
    for seed_idx in seeds:
        spawn = _spawn_fault_system if rng.random() < SYSTEM_SPAWN_FRACTION else _spawn_one_or_set
        before = len(world.faults)
        spawn(
            world, plate, rng,
            seed_world=own_points[seed_idx],
            seed_dist_rad=float(dist[seed_idx]),
            seed_weight=float(weight[seed_idx]),
            nn_omega=nn_omega[seed_idx],
            nn_point=nn_pts[seed_idx],
        )
        # Tag the traces this spawn just appended (a lone fault, a whole set, or a system's
        # strand family) so _generate_earthquakes can rupture overlap-born faults harder.
        if bool(overlap_mask[seed_idx]):
            for fault in world.faults[before:]:
                fault.born_in_overlap = True


def _regime_from_closing(closing_rate_rad_per_yr: float) -> str:
    thr = boundary.TRANSFORM_RATE_THRESHOLD
    if closing_rate_rad_per_yr > thr:
        return _KIND_REVERSE
    if closing_rate_rad_per_yr < -thr:
        return _KIND_NORMAL
    return _KIND_STRIKE_SLIP


def _spawn_one_or_set(
    world: "World", plate: Plate, rng: np.random.Generator,
    seed_world: np.ndarray, seed_dist_rad: float, seed_weight: float,
    nn_omega: np.ndarray, nn_point: np.ndarray,
) -> None:
    seed_world = seed_world / np.linalg.norm(seed_world)
    closing = float(
        boundary.closing_rate(seed_world[None], np.asarray(plate.omega, dtype=float), nn_omega[None], nn_point[None])[0]
    )
    kind = _regime_from_closing(closing)

    # Local tangent-plane axes at the seed: t_perp points toward the boundary (the local
    # shortening / extension direction), t_par runs along it.
    toward = nn_point - seed_world
    t_perp = toward - np.dot(toward, seed_world) * seed_world
    n = np.linalg.norm(t_perp)
    if n < 1e-9:
        east, _ = geometry.local_tangent_basis(seed_world)
        t_perp = east
    else:
        t_perp = t_perp / n
    t_par = np.cross(seed_world, t_perp)
    t_par = t_par / max(np.linalg.norm(t_par), 1e-12)

    if kind == _KIND_STRIKE_SLIP:
        sign = 1.0 if rng.random() < 0.5 else -1.0
        ang = np.radians(STRIKE_SLIP_OBLIQUITY_DEG) * sign
        strike = np.cos(ang) * t_par + np.sin(ang) * t_perp
        strike = strike / max(np.linalg.norm(strike), 1e-12)
    else:
        strike = t_par

    n_members = 1
    set_id: int | None = None
    if rng.random() < SET_PROBABILITY:
        n_members = int(rng.integers(SET_MIN_MEMBERS, SET_MAX_MEMBERS + 1))
        set_id = world.next_fault_id  # the first member's id doubles as the family id

    length_km = float(np.clip(rng.lognormal(np.log(LENGTH_MEDIAN_KM), LENGTH_SIGMA), LENGTH_MIN_KM, LENGTH_MAX_KM))
    slip_rate = SLIP_RATE_MIN_M_PER_MYR + (SLIP_RATE_MAX_M_PER_MYR - SLIP_RATE_MIN_M_PER_MYR) * (seed_weight**2)
    slip_rate *= float(rng.uniform(0.6, 1.4))
    lifespan = LIFESPAN_MIN_MYR + (LIFESPAN_MAX_MYR - LIFESPAN_MIN_MYR) * (0.3 + 0.7 * seed_weight)
    lifespan *= float(rng.uniform(0.7, 1.3))
    dip_base = {_KIND_NORMAL: DIP_NORMAL_DEG, _KIND_REVERSE: DIP_REVERSE_DEG, _KIND_STRIKE_SLIP: DIP_STRIKE_SLIP_DEG}[kind]

    spacing_rad = SET_SPACING_KM / PLANET_RADIUS_KM
    echelon_rad = SET_ECHELON_STEP_KM / PLANET_RADIUS_KM
    for member in range(n_members):
        centre = seed_world
        if member > 0:
            offset = (member - (n_members - 1) / 2.0) * spacing_rad
            along = (member - (n_members - 1) / 2.0) * echelon_rad
            centre = seed_world + offset * t_perp + along * strike
            centre = centre / np.linalg.norm(centre)
        member_sense = 1 if rng.random() < 0.5 else -1
        if kind == _KIND_NORMAL and n_members > 1:
            member_sense = 1 if member % 2 == 0 else -1  # alternating horst / graben polarity
        dip_dir_world = member_sense * t_perp
        fault = _build_fault(
            world, plate, rng, centre, strike, t_perp, kind,
            length_km=length_km,
            slip_rate=float(np.clip(slip_rate, SLIP_RATE_MIN_M_PER_MYR, SLIP_RATE_MAX_M_PER_MYR)),
            dip_deg=float(dip_base + rng.uniform(-5.0, 5.0)),
            strike_sense=member_sense,
            dip_dir_world=dip_dir_world,
            lifespan_myr=float(np.clip(lifespan, LIFESPAN_MIN_MYR, LIFESPAN_MAX_MYR + 15.0)),
            seed_dist_km=seed_dist_rad * PLANET_RADIUS_KM,
            set_id=set_id,
        )
        world.faults.append(fault)
        world.next_fault_id += 1

    label = f"set of {n_members}" if n_members > 1 else "single"
    world.log_event(
        f"fault system born on plate {plate.plate_id} ({_KIND_LABEL[kind]}, {label}, ~{length_km:.0f} km, "
        f"{seed_dist_rad * PLANET_RADIUS_KM:.0f} km from boundary)"
    )


def _seed_frame(plate: Plate, seed_world: np.ndarray, nn_omega: np.ndarray, nn_point: np.ndarray):
    """(seed_world, closing_rate, kind, t_perp, t_par) at a spawn seed. `t_perp` points
    toward the nearest boundary (local shortening/extension direction), `t_par` runs along
    it. Shared by the lone-fault and the fault-system spawn paths."""
    seed_world = seed_world / np.linalg.norm(seed_world)
    closing = float(
        boundary.closing_rate(seed_world[None], np.asarray(plate.omega, dtype=float), nn_omega[None], nn_point[None])[0]
    )
    kind = _regime_from_closing(closing)
    toward = nn_point - seed_world
    t_perp = toward - np.dot(toward, seed_world) * seed_world
    n = np.linalg.norm(t_perp)
    if n < 1e-9:
        east, _ = geometry.local_tangent_basis(seed_world)
        t_perp = east
    else:
        t_perp = t_perp / n
    t_par = np.cross(seed_world, t_perp)
    t_par = t_par / max(np.linalg.norm(t_par), 1e-12)
    return seed_world, closing, kind, t_perp, t_par


def _spawn_fault_system(
    world: "World", plate: Plate, rng: np.random.Generator,
    seed_world: np.ndarray, seed_dist_rad: float, seed_weight: float,
    nn_omega: np.ndarray, nn_point: np.ndarray,
) -> None:
    """A whole fault zone: one long, gently curving master lineament plus a family of
    sub-parallel strands scattered along its belt. Signature matches `_spawn_one_or_set` so
    `_maybe_spawn_faults` can pick either."""
    seed_world, _closing, kind, t_perp, t_par = _seed_frame(plate, seed_world, nn_omega, nn_point)

    if kind == _KIND_STRIKE_SLIP:
        sign = 1.0 if rng.random() < 0.5 else -1.0
        ang = np.radians(STRIKE_SLIP_OBLIQUITY_DEG) * sign
        master_strike = np.cos(ang) * t_par + np.sin(ang) * t_perp
        master_strike = master_strike / max(np.linalg.norm(master_strike), 1e-12)
    else:
        master_strike = t_par

    length_km = float(
        np.clip(
            rng.lognormal(np.log(SYSTEM_LENGTH_MEDIAN_KM), SYSTEM_LENGTH_SIGMA),
            SYSTEM_LENGTH_MIN_KM,
            SYSTEM_LENGTH_MAX_KM,
        )
    )
    length_rad = length_km / PLANET_RADIUS_KM

    # Master lineament: a great-circle arc through the seed along master_strike, warped by a
    # few low-frequency lateral lobes so the belt gently curves rather than running straight.
    n_master = int(np.clip(round(length_km / SYSTEM_MASTER_NODE_KM) + 1, 6, SYSTEM_MASTER_NODES_MAX))
    u = np.linspace(-length_rad / 2.0, length_rad / 2.0, n_master)
    s01 = (u - u[0]) / max(u[-1] - u[0], 1e-12)
    bend = np.zeros(n_master)
    for lobe in range(1, SYSTEM_BEND_LOBES + 1):
        amp = rng.uniform(-1.0, 1.0) * SYSTEM_BEND_MAX_FRACTION * length_rad / lobe
        bend += amp * np.sin(lobe * np.pi * s01 + rng.uniform(0.0, np.pi))
    master = (
        np.cos(u)[:, None] * seed_world[None, :]
        + np.sin(u)[:, None] * master_strike[None, :]
        + bend[:, None] * t_perp[None, :]
    )
    master = master / np.linalg.norm(master, axis=-1, keepdims=True)

    system_id = world.next_fault_system_id
    world.next_fault_system_id += 1
    system = FaultSystem(
        system_id=system_id,
        plate_id=plate.plate_id,
        kind=kind,
        master_local_phi=np.zeros(0),
        master_local_theta=np.zeros(0),
        length_km=length_km,
        birth_years=world.elapsed_years,
        lifespan_myr=float(rng.uniform(SYSTEM_LIFESPAN_MIN_MYR, SYSTEM_LIFESPAN_MAX_MYR)),
    )
    local_master = geometry.to_local(plate.frame, master)
    system.master_local_phi, system.master_local_theta = geometry.xyz_to_latlon(local_master)
    system.world_polyline = master
    world.fault_systems.append(system)

    # Strands: seeded at points stepped along the master, offset across the belt, striking
    # along the local master tangent with a little jitter.
    n_strands = int(
        np.clip(
            round(length_km / SYSTEM_STRAND_SPACING_KM * rng.uniform(0.7, 1.3)),
            SYSTEM_STRAND_COUNT_MIN,
            SYSTEM_STRAND_COUNT_MAX,
        )
    )
    tangents = np.gradient(master, axis=0)
    tangents = tangents / np.clip(np.linalg.norm(tangents, axis=-1, keepdims=True), 1e-12, None)
    dip_base = {_KIND_NORMAL: DIP_NORMAL_DEG, _KIND_REVERSE: DIP_REVERSE_DEG, _KIND_STRIKE_SLIP: DIP_STRIKE_SLIP_DEG}
    half_width_rad = SYSTEM_BELT_HALF_WIDTH_KM / PLANET_RADIUS_KM

    for i in range(n_strands):
        frac = (i + rng.uniform(-0.35, 0.35)) / max(n_strands - 1, 1)
        frac = float(np.clip(frac, 0.0, 1.0))
        m_idx = int(round(frac * (n_master - 1)))
        c0 = master[m_idx]
        tang = tangents[m_idx]
        perp = np.cross(c0, tang)
        perp = perp / max(np.linalg.norm(perp), 1e-12)
        centre = c0 + rng.uniform(-1.0, 1.0) * half_width_rad * perp
        centre = centre / np.linalg.norm(centre)

        jitter = np.radians(rng.uniform(-SYSTEM_STRAND_STRIKE_JITTER_DEG, SYSTEM_STRAND_STRIKE_JITTER_DEG))
        strike = np.cos(jitter) * tang + np.sin(jitter) * perp
        strike = strike / max(np.linalg.norm(strike), 1e-12)

        strand_kind = kind
        if rng.random() < SYSTEM_OFFREGIME_FRACTION:
            strand_kind = rng.choice([_KIND_NORMAL, _KIND_REVERSE, _KIND_STRIKE_SLIP])

        strand_len = float(
            np.clip(
                rng.lognormal(np.log(SYSTEM_STRAND_LENGTH_MEDIAN_KM), SYSTEM_STRAND_LENGTH_SIGMA),
                SYSTEM_STRAND_LENGTH_MIN_KM,
                SYSTEM_STRAND_LENGTH_MAX_KM,
            )
        )
        slip_rate = SLIP_RATE_MIN_M_PER_MYR + (SLIP_RATE_MAX_M_PER_MYR - SLIP_RATE_MIN_M_PER_MYR) * (seed_weight**2)
        slip_rate *= float(rng.uniform(0.6, 1.4))
        lifespan = LIFESPAN_MIN_MYR + (LIFESPAN_MAX_MYR - LIFESPAN_MIN_MYR) * (0.3 + 0.7 * seed_weight)
        lifespan *= float(rng.uniform(0.7, 1.3))
        member_sense = 1 if rng.random() < 0.5 else -1

        fault = _build_fault(
            world, plate, rng, centre, strike, perp, strand_kind,
            length_km=strand_len,
            slip_rate=float(np.clip(slip_rate, SLIP_RATE_MIN_M_PER_MYR, SLIP_RATE_MAX_M_PER_MYR)),
            dip_deg=float(dip_base[strand_kind] + rng.uniform(-5.0, 5.0)),
            strike_sense=member_sense,
            dip_dir_world=member_sense * perp,
            lifespan_myr=float(np.clip(lifespan, LIFESPAN_MIN_MYR, LIFESPAN_MAX_MYR + 15.0)),
            seed_dist_km=seed_dist_rad * PLANET_RADIUS_KM,
            set_id=None,
            system_id=system_id,
            node_km=SYSTEM_STRAND_NODE_KM,
            max_nodes=SYSTEM_STRAND_NODES_MAX,
        )
        world.faults.append(fault)
        world.next_fault_id += 1

    world.log_event(
        f"fault system #{system_id} born on plate {plate.plate_id} ({_KIND_LABEL[kind]}, "
        f"~{length_km:.0f} km master, {n_strands} strands, "
        f"{seed_dist_rad * PLANET_RADIUS_KM:.0f} km from boundary)"
    )


def _build_fault(
    world: "World", plate: Plate, rng: np.random.Generator,
    centre: np.ndarray, strike: np.ndarray, t_perp: np.ndarray, kind: str,
    length_km: float, slip_rate: float, dip_deg: float, strike_sense: int,
    dip_dir_world: np.ndarray, lifespan_myr: float, seed_dist_km: float, set_id: int | None,
    system_id: int | None = None, node_km: float = 15.0, max_nodes: int = FAULT_NODES_MAX,
) -> Fault:
    length_rad = length_km / PLANET_RADIUS_KM
    # ~1 node per `node_km` of trace, clamped -- enough to render a gently curved line and
    # give the relief query even coverage along a long strand.
    n_nodes = int(np.clip(round(length_km / node_km) + 1, FAULT_NODES_MIN, max_nodes))
    u = np.linspace(-length_rad / 2.0, length_rad / 2.0, n_nodes)
    base = np.cos(u)[:, None] * centre[None, :] + np.sin(u)[:, None] * strike[None, :]
    bend_amp = rng.uniform(-BEND_MAX_FRACTION, BEND_MAX_FRACTION) * length_rad
    bend = bend_amp * np.sin(np.pi * (u - u[0]) / max(length_rad, 1e-9))
    pts = base + bend[:, None] * t_perp[None, :]
    pts = pts / np.linalg.norm(pts, axis=-1, keepdims=True)

    local = geometry.to_local(plate.frame, pts)
    phi, theta = geometry.xyz_to_latlon(local)
    dip_dir_local = geometry.to_local(plate.frame, dip_dir_world / max(np.linalg.norm(dip_dir_world), 1e-12))

    return Fault(
        fault_id=world.next_fault_id,
        plate_id=plate.plate_id,
        kind=kind,
        local_phi=phi,
        local_theta=theta,
        slip_rate_m_per_myr=slip_rate,
        dip_deg=dip_deg,
        strike_sense=int(strike_sense),
        dip_dir_local=dip_dir_local,
        lifespan_myr=lifespan_myr,
        birth_years=world.elapsed_years,
        birth_distance_from_boundary_km=seed_dist_km,
        set_id=set_id,
        system_id=system_id,
    )


def _cull_scars(world: "World") -> None:
    by_plate: dict[int, list[Fault]] = {}
    for fault in world.faults:
        if not fault.active:
            by_plate.setdefault(fault.plate_id, []).append(fault)
    drop: set[int] = set()
    for scars in by_plate.values():
        if len(scars) <= MAX_SCARS_PER_PLATE:
            continue
        scars.sort(key=lambda f: f.age_myr, reverse=True)  # oldest first
        for fault in scars[MAX_SCARS_PER_PLATE:]:
            drop.add(id(fault))
    if drop:
        world.faults = [f for f in world.faults if id(f) not in drop]


def _cull_inactive_systems(world: "World") -> None:
    """Keep at most MAX_INACTIVE_SYSTEMS_PER_PLATE inert system scars per plate (oldest
    first); active systems are never culled. A dropped system's strands are untouched --
    they age out and get culled on their own by `_cull_scars`."""
    by_plate: dict[int, list[FaultSystem]] = {}
    for system in world.fault_systems:
        if not system.active:
            by_plate.setdefault(system.plate_id, []).append(system)
    drop: set[int] = set()
    for scars in by_plate.values():
        if len(scars) <= MAX_INACTIVE_SYSTEMS_PER_PLATE:
            continue
        scars.sort(key=lambda s: s.age_myr, reverse=True)
        for system in scars[MAX_INACTIVE_SYSTEMS_PER_PLATE:]:
            drop.add(id(system))
    if drop:
        world.fault_systems = [s for s in world.fault_systems if id(s) not in drop]


# --------------------------------------------------------------------------- relief


def _relief_mode_scales(world: "World") -> tuple[float, float]:
    """(rate_scale, reach_scale) for `_apply_plate_fault_relief`, keyed off
    `world.fault_deformation_mode`. Both 1.0 in the default "boundary" mode -- bit-identical
    to before the mode existed. In "fault"/"both" mode the fault-relief layer is scaled up to
    carry the deformation the smooth boundary bands give up (see FAULT_RELIEF_MODE_*)."""
    if getattr(world, "fault_deformation_mode", "boundary") in ("fault", "both"):
        return FAULT_RELIEF_MODE_RATE_SCALE, FAULT_RELIEF_MODE_REACH_SCALE
    return 1.0, 1.0


def fault_influence(
    world: "World", plate: Plate, own_points: np.ndarray,
    reach_km: float = FAULT_DEFORM_REACH_KM, floor: float = FAULT_DEFORM_FLOOR,
) -> np.ndarray:
    """Per own-node, 1.0 within `reach_km` of one of this plate's active fault traces,
    tapering linearly to `floor` beyond -- never 0, so a contested zone that has no fault yet
    still deforms while Piece-1 spawning fills it in. All-ones if the plate has no active
    fault. Used by LithospherePlate.deform in "fault" mode to localise boundary thickening
    onto fault lines (see World.fault_deformation_mode). Uses each fault's `world_polyline`
    (refreshed at the end of update_faults last step -- deform runs before update_faults, so
    "last step's faults" is the right, and only available, set)."""
    if len(own_points) == 0:
        return np.ones(0)
    traces = [
        f.world_polyline
        for f in world.faults
        if f.plate_id == plate.plate_id and f.active and f.world_polyline is not None
    ]
    if not traces:
        return np.ones(len(own_points))
    reach_rad = reach_km / PLANET_RADIUS_KM
    trace_points = np.concatenate(traces, axis=0)
    d, _ = cKDTree(trace_points).query(own_points, workers=query_workers(len(own_points)))
    return np.clip(1.0 - d / reach_rad, floor, 1.0)


def _apply_plate_fault_relief(world: "World", plate: Plate, years_myr: float) -> None:
    if not hasattr(plate, "lines"):
        return
    active = [f for f in world.faults if f.plate_id == plate.plate_id and f.active]
    if not active:
        return
    own_points = plate.all_points_and_elevation()[0]
    if len(own_points) == 0:
        return
    rate_scale, reach_scale = _relief_mode_scales(world)
    tree = cKDTree(own_points, balanced_tree=False, compact_nodes=False)
    reach_rad = reach_scale * MAX_FAULT_REACH_KM / PLANET_RADIUS_KM

    delta = np.zeros(len(own_points))
    reason = np.zeros(len(own_points), dtype=float)
    for fault in active:
        trace = fault_world_points(fault, plate)
        neighbours = tree.query_ball_point(trace, reach_rad)
        affected = sorted({i for sub in neighbours for i in sub})
        if not affected:
            continue
        affected = np.array(affected)
        pts = own_points[affected]
        d, _ = cKDTree(trace).query(pts)
        taper = np.clip(1.0 - d / reach_rad, 0.0, 1.0)
        slip_norm = float(np.clip(fault.slip_rate_m_per_myr / SLIP_RATE_REF_M_PER_MYR, 0.2, 3.0))
        mag = taper * slip_norm * years_myr * rate_scale

        if fault.kind == _KIND_REVERSE:
            contrib = REVERSE_UPLIFT_M_PER_MYR * mag
        elif fault.kind == _KIND_NORMAL:
            dip_dir_world = geometry.to_world(plate.frame, fault.dip_dir_local)
            mid = trace[len(trace) // 2]
            hanging = (pts - mid) @ dip_dir_world > 0.0
            contrib = np.where(hanging, -NORMAL_THROW_M_PER_MYR, NORMAL_SHOULDER_UPLIFT_M_PER_MYR) * mag
        else:  # strike-slip
            contrib = (STRIKE_SLIP_RIDGE_M_PER_MYR + fault.strike_sense * STRIKE_SLIP_BEND_M_PER_MYR) * mag

        delta[affected] += contrib
        reason[affected] = _KIND_REASON[fault.kind]

    if not np.any(delta):
        return

    new_lines = []
    offset = 0
    changed = False
    for line in plate.lines:
        n = len(line)
        if n == 0:
            new_lines.append(line)
            continue
        seg_delta = delta[offset : offset + n]
        seg_reason = reason[offset : offset + n]
        offset += n
        if not np.any(seg_delta):
            new_lines.append(line)
            continue
        new_elev = np.clip(line.elevation + seg_delta, MIN_ELEVATION_M, MAX_ELEVATION_M)
        moved = np.abs(new_elev - line.elevation) >= ELEV_CHANGE_MIN_DELTA_M
        new_reason = np.where(moved & (seg_reason > 0), seg_reason, line.elev_change_reason)
        new_lines.append(line.replace(elevation=new_elev, elev_change_reason=new_reason))
        changed = True
    if changed:
        plate.set_lines(new_lines)


# --------------------------------------------------------------------------- earthquakes


def _earthquake_magnitude(length_km: float, slip_rate_m_per_myr: float, born_in_overlap: bool, jitter: float) -> float:
    """Moment magnitude from trace length and slip rate (see QUAKE_MW_* constants), plus an
    overlap bonus and a per-event jitter."""
    mw = (
        QUAKE_MW_BASE
        + QUAKE_MW_LENGTH_COEFF * np.log10(max(length_km, 1.0))
        + QUAKE_MW_SLIP_COEFF * np.log10(max(slip_rate_m_per_myr, 1.0))
        + (OVERLAP_QUAKE_MW_BONUS if born_in_overlap else 0.0)
        + jitter
    )
    return float(np.clip(mw, QUAKE_MW_MIN, QUAKE_MW_MAX))


def _generate_earthquakes(world: "World", years_myr: float) -> None:
    """One characteristic `Earthquake` per active fault that accumulated at least
    `MIN_STEP_SLIP_FOR_QUAKE_M` of slip this step (a real active fault ruptures thousands of
    times per Myr -- one representative event per step is all we keep). Epicentre is a node
    drawn along the trace; magnitude from length + slip + an overlap bonus + a small jitter.
    Deterministic per `(seed, round(elapsed_years), fault_id, _QUAKE_SEED_TAG)` so a replayed
    session produces identical earthquakes."""
    by_id = plate_by_id(world)
    step_quakes: list[Earthquake] = []
    for fault in world.faults:
        if not fault.active or fault.slip_rate_m_per_myr * years_myr < MIN_STEP_SLIP_FOR_QUAKE_M:
            continue
        plate = by_id.get(fault.plate_id)
        if plate is None:
            continue
        rng = np.random.default_rng((world.seed, round(world.elapsed_years), fault.fault_id, _QUAKE_SEED_TAG))
        trace = fault.world_polyline
        if trace is None or len(trace) == 0:
            trace = fault_world_points(fault, plate)
        magnitude = _earthquake_magnitude(
            fault.length_km(), fault.slip_rate_m_per_myr, fault.born_in_overlap,
            float(rng.uniform(-QUAKE_MW_JITTER, QUAKE_MW_JITTER)),
        )
        epicenter = trace[int(rng.integers(len(trace)))]
        epicenter = epicenter / max(float(np.linalg.norm(epicenter)), 1e-12)
        quake = Earthquake(
            earthquake_id=world.next_earthquake_id,
            fault_id=fault.fault_id,
            plate_id=fault.plate_id,
            kind=fault.kind,
            epicenter_world=epicenter,
            magnitude=magnitude,
            slip_m=fault.slip_rate_m_per_myr * years_myr,
            birth_years=world.elapsed_years,
        )
        world.earthquakes.append(quake)
        world.next_earthquake_id += 1
        step_quakes.append(quake)

    # One console line per step at most: the step's largest, if it's a major event.
    if step_quakes:
        biggest = max(step_quakes, key=lambda q: q.magnitude)
        if biggest.magnitude >= EARTHQUAKE_LOG_MIN_MW:
            world.log_event(
                f"M{biggest.magnitude:.1f} earthquake on fault #{biggest.fault_id} "
                f"({_KIND_LABEL[biggest.kind]}) on plate {biggest.plate_id}"
            )

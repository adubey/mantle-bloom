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

This module is an **additive** layer: it never touches deform()'s own boundary
classification. Each step it

1. ages every existing fault, accumulating slip on the active ones and retiring those past
   their drawn lifespan (kept forever after as an inactive scar, like `is_volcano`);
2. rolls a stress-weighted Poisson spawn per plate -- probability high in a band near the
   boundary, decaying exponentially into the interior with a small nonzero floor everywhere,
   the regime (normal / reverse / strike-slip) picked from the local closing rate per
   Andersonian faulting theory, and with `SET_PROBABILITY` chance a whole sub-parallel
   family rather than a lone trace;
3. applies each active fault's own relief to the nearby crust -- reverse: an uplift ridge;
   normal: a hanging-wall graben with a footwall shoulder; strike-slip: a modest
   transpressional ridge or transtensional sag (relief only -- the node field is *not*
   physically sheared across the trace, see docs/TODO.md).

Geometry is stored in the owning plate's **local frame** (`local_phi` / `local_theta`), so
a fault rides along with the crust as the plate rotates for free -- the same "attached to
the crust, not the world" property every persistent `ElevationLine` field already has.
`reconcile_faults` re-homes faults across merges/splits and drops those whose plate
subducted (see world.step_world).
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
)
from .plates import Plate, collect_all_points, query_workers

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

# Fault sets -- sub-parallel families. Basin-and-Range major normal faults sit ~15-30 km
# apart; en echelon step-overs ~1-5 km.
SET_PROBABILITY = 0.4
SET_MIN_MEMBERS = 2
SET_MAX_MEMBERS = 5
SET_SPACING_KM = 20.0
SET_ECHELON_STEP_KM = 3.0

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
    age_myr: float = 0.0
    cumulative_offset_m: float = 0.0
    active: bool = True
    # Refreshed at the end of update_faults every step (world-space polyline, true frame) so
    # reconcile_faults can re-home a fault even after its plate was absorbed by a merge.
    # Not authoritative geometry -- local_phi/local_theta are.
    world_polyline: np.ndarray | None = field(default=None, repr=False, compare=False)

    def length_km(self) -> float:
        pts = geometry.local_xyz(self.local_phi, self.local_theta)
        seg = np.arccos(np.clip(np.sum(pts[:-1] * pts[1:], axis=-1), -1.0, 1.0))
        return float(np.sum(seg) * PLANET_RADIUS_KM)


def plate_by_id(world: "World") -> dict[int, Plate]:
    return {p.plate_id: p for p in world.plates}


def fault_world_points(fault: Fault, plate: Plate) -> np.ndarray:
    """The fault's trace in true world coordinates, via the owning plate's current frame."""
    return geometry.to_world(plate.frame, geometry.local_xyz(fault.local_phi, fault.local_theta))


# --------------------------------------------------------------------------- step entry point


def update_faults(world: "World", years: float) -> None:
    """Age / spawn / retire faults and apply their relief. Called once per step from
    world.step_world, inside the simulate_plate_movement block, right after the deform loop
    and before merge_split.apply_topology_changes."""
    years_myr = years / 1_000_000.0

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

    total_nodes = max(1, sum(p.node_count() for p in world.plates))
    for plate in world.plates:
        _maybe_spawn_faults(world, plate, total_nodes, years_myr)

    _cull_scars(world)

    for plate in world.plates:
        _apply_plate_fault_relief(world, plate, years_myr)

    by_id = plate_by_id(world)
    for fault in world.faults:
        plate = by_id.get(fault.plate_id)
        if plate is not None:
            fault.world_polyline = fault_world_points(fault, plate)


def reconcile_faults(world: "World") -> None:
    """After a topology change: drop faults whose plate subducted, and re-home faults whose
    trace midpoint now lies in a different surviving plate's territory (covers both a merge
    -- the absorbed plate's id vanishes -- and a split -- a new id appears near the cut).
    Recomputes local coordinates in the new plate's frame; preserves age / offset / active.
    Mirrors stranded_basins.reconcile_world_tracks."""
    if not world.faults:
        return
    live = plate_by_id(world)
    if not live:
        world.faults = []
        return
    collected = collect_all_points(world.plates)
    if collected is None:
        world.faults = []
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

    area_frac = len(own_points) / total_nodes
    expected = BASE_SPAWN_RATE_PER_MYR * area_frac * float(np.mean(weight)) * years_myr
    n_spawn = int(rng.poisson(max(expected, 0.0)))
    if n_spawn == 0:
        return

    probs = weight / weight.sum()
    seeds = rng.choice(len(own_points), size=n_spawn, p=probs)
    for seed_idx in seeds:
        _spawn_one_or_set(
            world, plate, rng,
            seed_world=own_points[seed_idx],
            seed_dist_rad=float(dist[seed_idx]),
            seed_weight=float(weight[seed_idx]),
            nn_omega=nn_omega[seed_idx],
            nn_point=nn_pts[seed_idx],
        )


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


def _build_fault(
    world: "World", plate: Plate, rng: np.random.Generator,
    centre: np.ndarray, strike: np.ndarray, t_perp: np.ndarray, kind: str,
    length_km: float, slip_rate: float, dip_deg: float, strike_sense: int,
    dip_dir_world: np.ndarray, lifespan_myr: float, seed_dist_km: float, set_id: int | None,
) -> Fault:
    length_rad = length_km / PLANET_RADIUS_KM
    # ~1 node per 15 km of trace, clamped -- enough to render a gently curved line.
    n_nodes = int(np.clip(round(length_km / 15.0) + 1, FAULT_NODES_MIN, FAULT_NODES_MAX))
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


# --------------------------------------------------------------------------- relief


def _apply_plate_fault_relief(world: "World", plate: Plate, years_myr: float) -> None:
    if not hasattr(plate, "lines"):
        return
    active = [f for f in world.faults if f.plate_id == plate.plate_id and f.active]
    if not active:
        return
    own_points = plate.all_points_and_elevation()[0]
    if len(own_points) == 0:
        return
    tree = cKDTree(own_points, balanced_tree=False, compact_nodes=False)
    reach_rad = MAX_FAULT_REACH_KM / PLANET_RADIUS_KM

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
        mag = taper * slip_norm * years_myr

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

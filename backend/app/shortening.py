"""Collisional shortening carried into plate interiors (GitHub issue #314).

A plate's convergent band demands shortening every step: the plastic strain the band would
take on locally (`rheology.convergent_strain`), times each cell's area. That demand is
*shortening area* -- metres of convergence times metres of boundary, m^2 -- and it is what
this module conserves, not volume. A cell that takes up `A` m^2 of it thickens its whole column
by the factor `1 + A / area` (pure shear, Hc and Hm together), the same multiplier the
band-local `rheology.apply_convergent_deformation` applies, so unequal cells and node densities
are handled by each cell's own area.

The accommodation interface is `accommodate(problem) -> ShorteningResult`. Option 1 of the
issue, `ring_cascade`, is the only implementation: it walks rings of cells outward from the
collision front, each cell keeping a share of what reaches it by its strength, its relief and
the room left under the Hc/Hm caps, and passing the rest on. An intraplate stress solve (the
issue's option 2) would plug in behind the same signature.

Whatever the plate doesn't absorb -- lost to basal drag along the way (`returned_decay_m2`) or
with nowhere left to go (`returned_unrouted_m2`: past `CASCADE_MAX_REACH_KM`, a dead end, or a
convergent cell the front can't reach) -- stays in the boundary overlap, where boundary retreat
and suture accretion consume it, as convergence the band never took up always has. It is
booked, never clipped silently (`phase_budget.record_shortening`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.sparse import csr_matrix, diags

from . import lithosphere, rheology
from .elevation_lines import DEFAULT_NODE_DENSITY, line_spacing_rad
from .mantle import PLANET_RADIUS_KM

# --- Knobs (issue #314) -------------------------------------------------------------------

# How far from the front shortening may travel before whatever is left returns to the
# boundary: far enough for Tien Shan-style reactivation 1,500-2,000 km behind a collision.
CASCADE_MAX_REACH_KM = 1500.0
# Stress lost to basal drag: the travelling remainder decays by exp(-h / this) per cell of
# width h it crosses.
CASCADE_DECAY_LENGTH_KM = 1500.0
# The e-folding width over which weak, unthickened ground takes up shortening: a cell of width
# h and acceptance factor w keeps 1 - exp(-h * w / this) of what reaches it. Set in km rather
# than as a per-cell fraction so the belt's width doesn't depend on node density.
CASCADE_ABSORPTION_LENGTH_KM = 250.0
# Yield stress of ordinary continental lithosphere, Y_0, and how much stronger a full-strength
# craton is (Y = Y_0 * (1 + gain * cratons.strength)). At 8 a full craton's yield stress is
# well above any collision's driving stress, so it passes shortening through undeformed.
CASCADE_BASE_YIELD_STRESS_PA = rheology.COHESION_PA
CASCADE_CRATON_YIELD_GAIN = 8.0
# The share of a column's relief load, rho_c * g * relief, that pushes back on further
# thickening. At 0.5 a typical collision stops thickening a column at ~7-8 km of relief and
# pushes the deformation outward instead -- how a belt widens rather than piling up.
CASCADE_RELIEF_PUSHBACK = 0.5
# A column takes up less and less of what reaches it once either Hc or Hm passes this fraction
# of its cap, reaching nothing at the cap itself, and routing steers away from it the same way
# -- so it hands shortening outward well before it fills up, instead of saturating at the cap.
CASCADE_ROOM_TAPER_ONSET = 0.75
# Routing weight of a fully resistant target, relative to a fully weak one of the same area:
# shortening goes around a craton when it can, and through it when every route forward is
# craton.
CASCADE_ROUTING_FLOOR = 1e-3
# Demand is normalized to the default density's band width: the convergent band is a fixed
# number of cells wide (torque.BOUNDARY_FORCE_REACH_MULTIPLIER), so without this a finer
# lattice would collide a thinner band and demand less shortening per km of boundary.
REFERENCE_SPACING_KM = line_spacing_rad(DEFAULT_NODE_DENSITY) * PLANET_RADIUS_KM


@dataclass
class ShorteningProblem:
    """One plate's shortening to accommodate this step. Every array is per cell."""

    adjacency: csr_matrix  # symmetric cell adjacency
    areas_m2: np.ndarray
    front: np.ndarray  # bool: the cells the rings count out from
    demand_m2: np.ndarray  # shortening area each cell's own convergence demands
    hc_m: np.ndarray
    hm_m: np.ndarray
    craton_strength: np.ndarray  # cratons.strength, 0..1
    relief_m: np.ndarray  # isostatic elevation above the continental reference, >= 0
    crust_density: np.ndarray
    drive_stress_pa: float
    spacing_rad: float
    reach_scale: float = 1.0  # World.collision_uplift_reach_multiplier
    # Cells the cascade may route through; the rest take up their own demand in place, as far
    # as `room_m2` allows. None: every cell.
    host: np.ndarray | None = None


@dataclass
class ShorteningResult:
    absorbed_m2: np.ndarray  # per cell
    returned_decay_m2: np.ndarray  # per cell: lost to basal drag on leaving it
    returned_unrouted_m2: np.ndarray  # per cell: left over here with no route onward

    @property
    def returned_m2(self) -> np.ndarray:
        return self.returned_decay_m2 + self.returned_unrouted_m2


Accommodation = Callable[[ShorteningProblem], ShorteningResult]


def demand_m2(strain: np.ndarray, areas_m2: np.ndarray, spacing_rad: float) -> np.ndarray:
    """Shortening area a convergent band's fractional `strain` demands, per cell."""
    return np.asarray(strain) * areas_m2 * (REFERENCE_SPACING_KM / (spacing_rad * PLANET_RADIUS_KM))


def drive_stress_pa(closing_rate_m_per_s: np.ndarray, demand: np.ndarray) -> float:
    """The collision's driving stress: the band's normal-stress proxy
    (`rheology.EFFECTIVE_LITHOSPHERE_VISCOSITY_PA_S_PER_M` times closing rate), averaged over
    the cells by their demand."""
    total = float(np.sum(demand))
    if total <= 0.0:
        return 0.0
    sigma = np.abs(closing_rate_m_per_s) * rheology.EFFECTIVE_LITHOSPHERE_VISCOSITY_PA_S_PER_M
    return float(np.dot(sigma, demand) / total)


def acceptance_factor(problem: ShorteningProblem) -> np.ndarray:
    """min(a_strength, a_gpe, a_room) per cell, each in [0, 1]: how readily a cell deforms
    under the collision's driving stress, against its own yield strength, its relief's pushback
    and how close it already is to the caps (CASCADE_ROOM_TAPER_ONSET)."""
    sigma = problem.drive_stress_pa
    if sigma <= 0.0:
        return np.zeros(len(problem.areas_m2))
    yield_pa = CASCADE_BASE_YIELD_STRESS_PA * (1.0 + CASCADE_CRATON_YIELD_GAIN * problem.craton_strength)
    a_strength = np.clip((sigma - yield_pa) / sigma, 0.0, 1.0)
    relief_pa = CASCADE_RELIEF_PUSHBACK * problem.crust_density * lithosphere.GRAVITY_M_S2 * problem.relief_m
    a_gpe = np.clip(1.0 - relief_pa / sigma, 0.0, 1.0)
    fullness = np.maximum(problem.hc_m / lithosphere.MAX_CRUSTAL_THICKNESS_M, problem.hm_m / lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)
    a_room = np.clip((1.0 - fullness) / (1.0 - CASCADE_ROOM_TAPER_ONSET), 0.0, 1.0)
    return np.minimum(np.minimum(a_strength, a_gpe), a_room)


def room_m2(hc_m: np.ndarray, hm_m: np.ndarray, areas_m2: np.ndarray) -> np.ndarray:
    """Shortening area a column can take up before either Hc or Hm reaches its cap."""
    hc_room = np.divide(lithosphere.MAX_CRUSTAL_THICKNESS_M, hc_m, out=np.full(len(hc_m), np.inf), where=hc_m > 0.0) - 1.0
    hm_room = np.divide(lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M, hm_m, out=np.full(len(hm_m), np.inf), where=hm_m > 0.0) - 1.0
    return areas_m2 * np.clip(np.minimum(hc_room, hm_room), 0.0, None)


def rings(adjacency: csr_matrix, front: np.ndarray, max_rings: int) -> np.ndarray:
    """Per cell, edge hops to the nearest `front` cell, saturating at `max_rings + 1`."""
    ring = np.where(front, 0, max_rings + 1)
    reached = np.asarray(front, dtype=bool).copy()
    frontier = reached.copy()
    for hop in range(1, max_rings + 1):
        frontier = (adjacency @ frontier.astype(np.int8) > 0) & ~reached
        if not np.any(frontier):
            break
        ring[frontier] = hop
        reached |= frontier
    return ring


def _routing(adjacency: csr_matrix, ring: np.ndarray, max_rings: int, areas_m2: np.ndarray, w: np.ndarray) -> csr_matrix:
    """Row-normalized routing from each ring-r cell to ring r + 1. Straight ahead (the cell's
    own ring r + 1 neighbours) by target area times `w + CASCADE_ROUTING_FLOOR`, so a craton
    ahead with nothing better beside it passes shortening through. One cell sideways (ring
    r + 1 neighbours of its ring-r neighbours) only toward targets weaker than straight ahead,
    by area times how much weaker -- so it finds its way around a strong patch, while uniform
    ground routes straight on and stays uniform."""
    n = len(ring)
    coo = adjacency.tocoo()
    i, j = coo.row, coo.col
    inside = (ring[i] <= max_rings) & (ring[j] <= max_rings)
    forward = inside & (ring[j] == ring[i] + 1)
    lateral = inside & (ring[j] == ring[i])
    step = csr_matrix((np.ones(int(forward.sum())), (i[forward], j[forward])), shape=(n, n))
    side = csr_matrix((np.ones(int(lateral.sum())), (i[lateral], j[lateral])), shape=(n, n))
    ahead_count = np.asarray(step.sum(axis=1)).ravel()
    ahead_w = np.divide(step @ w, ahead_count, out=np.zeros(n), where=ahead_count > 0.0)
    around = (side @ step).tocoo()
    gain = areas_m2[around.col] * np.clip(w[around.col] - ahead_w[around.row], 0.0, None)
    weighted = step @ diags(areas_m2 * (w + CASCADE_ROUTING_FLOOR)) + csr_matrix((gain, (around.row, around.col)), shape=(n, n))
    row_sum = np.asarray(weighted.sum(axis=1)).ravel()
    return (diags(np.divide(1.0, row_sum, out=np.zeros(n), where=row_sum > 0.0)) @ weighted).tocsr()


def ring_cascade(problem: ShorteningProblem) -> ShorteningResult:
    """Option 1: pass shortening outward from the front ring by ring. A ring-r cell takes in
    what reaches it plus its own demand, keeps `1 - exp(-h * w / L)` of it (w from
    `acceptance_factor`, L = CASCADE_ABSORPTION_LENGTH_KM) up to its `room_m2`, and passes the
    rest -- less exp(-h / CASCADE_DECAY_LENGTH_KM) to basal drag -- on to ring r + 1, split by
    target area times acceptance. Conserves demand = absorbed + returned exactly."""
    n = len(problem.areas_m2)
    demand = np.asarray(problem.demand_m2, dtype=float)
    absorbed = np.zeros(n)
    returned_decay = np.zeros(n)
    returned_unrouted = np.zeros(n)
    if n == 0 or not np.any(demand > 0.0) or problem.reach_scale <= 0.0:
        return ShorteningResult(absorbed, returned_decay, demand.copy())
    room = room_m2(problem.hc_m, problem.hm_m, problem.areas_m2)
    host = np.ones(n, dtype=bool) if problem.host is None else np.asarray(problem.host, dtype=bool)
    if not np.all(host):
        absorbed[~host] = np.minimum(demand[~host], room[~host])
        returned_unrouted[~host] = demand[~host] - absorbed[~host]
        demand = np.where(host, demand, 0.0)
    reach_scale = problem.reach_scale
    max_rings = max(1, round(CASCADE_MAX_REACH_KM * reach_scale / (problem.spacing_rad * PLANET_RADIUS_KM)))
    adjacency = problem.adjacency if problem.host is None else (diags(host.astype(float)) @ problem.adjacency @ diags(host.astype(float))).tocsr()
    adjacency.eliminate_zeros()
    ring = rings(adjacency, problem.front & host, max_rings)
    width_km = np.sqrt(problem.areas_m2) / 1e3
    w = acceptance_factor(problem)
    keep = 1.0 - np.exp(-width_km * w / CASCADE_ABSORPTION_LENGTH_KM)
    carry = np.exp(-width_km / (CASCADE_DECAY_LENGTH_KM * reach_scale))
    routing = _routing(adjacency, ring, max_rings, problem.areas_m2, w)
    routed = np.asarray(routing.sum(axis=1)).ravel() > 0.0
    inbound = routing.T.tocsr()

    reachable = host & (ring <= max_rings)
    returned_unrouted[host & ~reachable] = demand[host & ~reachable]
    incoming = np.where(reachable, demand, 0.0)
    order = np.argsort(ring, kind="stable")
    bounds = np.searchsorted(ring[order], np.arange(max_rings + 2))
    for r in range(max_rings + 1):
        cells = order[bounds[r]:bounds[r + 1]]
        if not len(cells):
            break
        arriving = incoming[cells]
        taken = np.minimum(keep[cells] * arriving, room[cells])
        absorbed[cells] = taken
        passing = (arriving - taken) * carry[cells]
        returned_decay[cells] = arriving - taken - passing
        onward = routed[cells]
        returned_unrouted[cells[~onward]] += passing[~onward]
        send = np.zeros(n)
        send[cells[onward]] = passing[onward]
        incoming += inbound @ send
    return ShorteningResult(absorbed, returned_decay, returned_unrouted)


def apply_strain(hc_m: np.ndarray, hm_m: np.ndarray, strain: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Thicken each column by its absorbed strain (absorbed shortening / area), pure shear.
    `room_m2` already kept every cell under both caps; the clip only trims rounding."""
    return (
        np.minimum(hc_m * (1.0 + strain), np.maximum(hc_m, lithosphere.MAX_CRUSTAL_THICKNESS_M)),
        np.minimum(hm_m * (1.0 + strain), np.maximum(hm_m, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)),
    )

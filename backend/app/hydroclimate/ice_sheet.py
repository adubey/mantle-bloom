"""Glacier ice: volume conserved per ice body, shaped by a plastic ice-sheet approximation
(issue #353).

Explicit shallow-ice flow is far too stiff for 100 kyr-Myr steps, so ice does not *flow* here
in the time-stepping sense. Each step instead:

1. **Groups ice into bodies**: connected components, over the symmetrised hydrology graph, of
   nodes that hold ice or have positive surface balance at the reference surface (bed + last
   step's ice).
2. **Conserves each body's volume** with a backward-Euler update, V(d) = V_open + dt * B(d),
   solved for `d`, the position of the body's margin in a nested family of extents (signed
   graph distance from the positive-balance area). B is evaluated on the solved surface, with
   temperature lapse-corrected up to the ice surface, so accumulation magnitude matters (not
   just its sign), and the higher-is-colder feedback gives Weertman (1961) hysteresis. The
   root nearest last step's extent is kept.
3. **Shapes the volume** with Nye's (1952) perfectly plastic profile, |grad s| = c / H with
   c = tau_b / (rho_i g), solved by a monotone, causal ordered-upwind (fast-marching-style)
   sweep on the irregular graph. The profile is a geometry approximation only, as ICESHEET
   (Gowan et al. 2016) uses it for reconstruction; how much ice exists comes from step 2.
4. **Diagnoses flux and sliding**: balance flux routed down the ice surface with
   multiple-direction weights (Le Brocq et al. 2006), turned into a depth-averaged velocity
   and a basal sliding speed on warm-based beds, which erosion.py's glacial abrasion reads.

Glacier mass balance is still in this model's own units (`hydrology.GLACIER_ACCUMULATION_RATE`
is ~50,000x slower than real snowfall; see issue #353's open decision), so ice is *not* assumed
to reach balance within a step, and the velocities here are a calibrated proxy, not m/yr of
real ice. Nothing here is persisted: the opening volume comes from `glacier_depth`, and the
previous extent from where that ice sits.
"""

from __future__ import annotations

import heapq
import logging
from dataclasses import dataclass, field

import numpy as np
from numba import njit
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra

from .. import geometry
from ..elevation_lines import PLANET_RADIUS_KM
from ..lithosphere import GRAVITY_M_S2, RHO_ICE, RHO_WATER
from . import hydrology
from .climate import LAPSE_RATE_C_PER_KM

_logger = logging.getLogger(__name__)

# Basal yield stress (Pa) of the plastic profile. A frozen bed, or a thawed one on bare rock,
# holds about 100 kPa (Nye's classic value, and the upper range ice-sheet reconstructions use);
# a thawed bed on thick soft sediment yields far sooner, around 20-50 kPa. Between the two the
# yield stress is interpolated by how much mobile cover the bed carries, saturating at
# ICE_SOFT_BED_COVER_M.
ICE_YIELD_STRESS_HARD_PA = 100_000.0
ICE_YIELD_STRESS_SOFT_PA = 30_000.0
ICE_SOFT_BED_COVER_M = 20.0

# Thermal state of the bed: warm-based once steady conduction brings it to the pressure-melting
# point, T_s + G * H / k >= -ICE_PRESSURE_MELT_C_PER_M * H. Ignoring advection makes this
# generous; at T_s = -20 C it thaws under about 700 m of ice. G is a constant for now; coupling
# it to Hm is a follow-up (issue #353, out of scope).
ICE_GEOTHERMAL_FLUX_W_M2 = 0.06
ICE_THERMAL_CONDUCTIVITY_W_M_K = 2.1
ICE_PRESSURE_MELT_C_PER_M = 8.7e-4

# Geothermal basal melt on a warm bed, G / (rho_i L_f), about 6 mm/yr of real ice. It is scaled
# into the same model units as accumulation: GLACIER_ACCUMULATION_RATE (m of ice per Myr per
# mm/yr of snowfall) over the physical rate, 1 mm/yr of water = 1000/917 mm/yr of ice =
# ~1090 m/Myr. Replaces the old depth-squared basal melt, whose only job was to stop ice
# piling up on flat beds where bed-slope flow left it -- the plastic geometry spreads it instead.
ICE_LATENT_HEAT_J_KG = 334_000.0
_SECONDS_PER_MYR = 3.15576e13
PHYSICAL_ACCUMULATION_M_PER_MM_MYR = 1000.0 / RHO_ICE * 1.0e-3 * 1.0e6
MODEL_UNIT_SCALE = hydrology.GLACIER_ACCUMULATION_RATE / PHYSICAL_ACCUMULATION_M_PER_MM_MYR
GEOTHERMAL_BASAL_MELT_M_PER_MYR = ICE_GEOTHERMAL_FLUX_W_M2 / (RHO_ICE * ICE_LATENT_HEAT_J_KG) * _SECONDS_PER_MYR * MODEL_UNIT_SCALE

# Temperature ramps (C) that keep the mass balance continuous in ice thickness, so the
# backward-Euler root find never lands on a jump (a jump in B is multiplied by dt, and a long
# step turns a small one into a large volume error): snow starts falling on a thickening sheet
# over the ICE_SNOW_RAMP_C below freezing, and geothermal melt starts over the
# ICE_THAW_RAMP_C of bed warming past the pressure-melting point.
ICE_SNOW_RAMP_C = 1.0
ICE_THAW_RAMP_C = 1.0

# A plastic profile on one or two cells is meaningless: a body with room for fewer cells than
# this (its domain: its own nodes plus the ground it could advance onto) spreads its volume
# evenly over its own nodes, under the same volume contract. The test is the room, not the
# body's current size: a one- or two-node body can hold kilometres of ice (lake and river ice
# froze in place like that until #539), and held to its own nodes it saturates and sheds the
# excess as meltwater every step.
ICE_BODY_MIN_NODES = 3

# Basal sliding: a share of the depth-averaged balance velocity on a warm bed, none on a frozen
# one. Velocities divide by thickness, floored here so a thin margin doesn't blow up.
ICE_SLIDING_FRACTION = 0.7
ICE_VELOCITY_MIN_THICKNESS_M = 10.0

# Root find for the margin position d (metres of graph distance). The bracket is grown outward
# from last step's extent, starting at about one edge length; the margin is then located to
# ICE_ROOT_TOLERANCE_M, a fraction of a percent of a 60 km cell.
ICE_ROOT_TOLERANCE_M = 100.0
ICE_ROOT_MAX_ITERATIONS = 40
# How many times an evaluation re-solves after dropping marine nodes the ice can't ground.
ICE_FLOTATION_PASSES = 6

_FLOTATION_RATIO = RHO_WATER / RHO_ICE


# --- RFC 0015 gate G0 prototype (#554; gate branch only, never merged) -------------------
# Every switch off reproduces origin/main. `bin/debug/replay_g0.py` sets them per arm.
@dataclass
class G0Config:
    physical: bool = False  # D1's physical mass balance (m/yr of ice, converted to m/Myr here)
    f_snow: bool = False  # arm (e): seasonal snow fraction instead of the is_frozen switch
    root_rule: str = "nearest"  # "nearest" (today), "largest" (h1) or "smallest" (h2) stable root
    enumerate_roots: bool = False  # scan every sign change on [d_min, d_max] and log it
    hold_deflection: bool = False  # arm (c): mass-balance surface on the run's opening deflection
    hold_sea_level: bool = False  # arm (b): read by the driver
    reevaluate_membership: bool = False  # arm (f): membership and yield stress at the solved thickness
    couple_sea_level: bool = False  # arm (b'): sea level solved within the step (D2's bracket)
    couple_deflection: bool = False  # arm (c'): mass-balance surface on the closing deflection
    current_friction: bool = False  # arm (j): frictional heat from the current step, by iteration
    fixed_yield: bool = False  # diagnostic arm (y): yield stress always the hard-bed value
    no_friction: bool = False  # diagnostic arm (q): no frictional heat
    march: bool = False  # arm (m): march to equilibrium in substeps, then hold
    march_substep_years: float = 1_000.0
    march_tolerance: float = 1.0e-3  # |dV|/V between two-substep means counted as settled
    march_patience: int = 3  # settled substeps in a row before holding
    march_max_years: float = 100_000.0  # longest march before holding regardless (logged)
    ddf_snow_mm: float = 4.0
    ddf_ice_mm: float = 8.0
    rain_snow_c: float = 1.0
    # Per-node frictional heat (W m^-2) from last step, aligned with step_ice's nodes; set by
    # hydrology before each solve.
    frictional_heat_w_m2: np.ndarray | None = None


G0 = G0Config()
G0_ROOT_LOG: list[dict] = []
G0_MEMBERSHIP_LOG: list[dict] = []
G0_OUTER_LOG: list[dict] = []
G0_MARCH_LOG: list[dict] = []
G0_OUTER_ITERATIONS = 5
G0_SEA_LEVEL_EVALUATIONS = 8
G0_SEA_LEVEL_TOLERANCE_M = 1.0
G0_SEA_LEVEL_BRACKET_M = 500.0
SECONDS_PER_YEAR = 3.15576e7
W_M2_TO_M_ICE_PER_YR = SECONDS_PER_YEAR / (RHO_ICE * ICE_LATENT_HEAT_J_KG)
G0_MEMBERSHIP_ITERATIONS = 5


def snow_fraction(temperature_c: np.ndarray, amplitude_c: np.ndarray, threshold_c: float) -> np.ndarray:
    """D1: share of a year-round sinusoid below the rain-snow threshold, with the continuous
    ICE_SNOW_RAMP_C ramp where the amplitude vanishes."""
    amp = np.asarray(amplitude_c, dtype=float)
    ratio = np.clip((threshold_c - temperature_c) / np.maximum(amp, 1e-6), -1.0, 1.0)
    seasonal = 0.5 + np.arcsin(ratio) / np.pi
    ramp = np.clip(0.5 + (threshold_c - temperature_c) / ICE_SNOW_RAMP_C, 0.0, 1.0)
    return np.where(amp > 1e-3, seasonal, ramp)


def basal_heat_w_m2(thickness_m: np.ndarray, surface_c: np.ndarray, frictional_w_m2) -> np.ndarray:
    """D1: net basal heat G + Q_f - k max(0, T_pmp - T_surf) / H (W m^-2); -inf without ice."""
    h = np.asarray(thickness_m, dtype=float)
    t_pmp = -ICE_PRESSURE_MELT_C_PER_M * h
    conducted = ICE_THERMAL_CONDUCTIVITY_W_M_K * np.clip(t_pmp - surface_c, 0.0, None) / np.maximum(h, 1e-9)
    q = ICE_GEOTHERMAL_FLUX_W_M2 + frictional_w_m2 - conducted
    return np.where(h > 0.0, q, -np.inf)


@dataclass
class SurfaceClimate:
    """Per-node climate the mass balance needs. `temperature_c` is the annual mean at the
    node's current (ice-free, already deflected) surface, so the ice surface sits
    LAPSE_RATE_C_PER_KM colder per km of ice."""

    precipitation_mm: np.ndarray
    temperature_c: np.ndarray
    seasonal_amplitude_c: np.ndarray
    is_frozen: np.ndarray
    frictional_heat_w_m2: np.ndarray | None = None  # G0 prototype: last step's tau_b u_b / Y

    def subset(self, idx: np.ndarray) -> "SurfaceClimate":
        return SurfaceClimate(
            self.precipitation_mm[idx], self.temperature_c[idx], self.seasonal_amplitude_c[idx], self.is_frozen[idx],
            None if self.frictional_heat_w_m2 is None else self.frictional_heat_w_m2[idx],
        )

    def _physical_rates(self, thickness: np.ndarray, surface_c: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """G0 prototype of RFC 0015 D1, in m/Myr of ice: snowfall by mass, degree-day melt of
        the year's snow then bare ice, and the basal energy balance."""
        to_ice = RHO_WATER / RHO_ICE
        if G0.f_snow:
            snows = snow_fraction(surface_c, self.seasonal_amplitude_c, G0.rain_snow_c)
        else:
            snows = np.where(self.is_frozen, 1.0, np.clip((hydrology.FREEZE_POINT_C - surface_c) / ICE_SNOW_RAMP_C, 0.0, 1.0))
        snow_we = snows * self.precipitation_mm * 1.0e-3  # m w.e. per yr
        pdd = 365.25 * hydrology.seasonal_positive_degrees(surface_c, self.seasonal_amplitude_c, 0.0)
        snow_capacity = G0.ddf_snow_mm * 1.0e-3 * pdd
        snow_melt = np.minimum(snow_we, snow_capacity)
        pdd_left = pdd - snow_melt / (G0.ddf_snow_mm * 1.0e-3)
        ice_melt = np.where(thickness > 0.0, G0.ddf_ice_mm * 1.0e-3 * pdd_left, 0.0)
        frictional = 0.0 if self.frictional_heat_w_m2 is None else self.frictional_heat_w_m2
        basal = W_M2_TO_M_ICE_PER_YR * np.clip(basal_heat_w_m2(thickness, surface_c, frictional), 0.0, None)
        return snow_we * to_ice * 1.0e6, (snow_melt + ice_melt) * to_ice * 1.0e6, basal * 1.0e6

    def rates(self, thickness_m: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(accumulation, surface melt, basal melt), metres of ice per Myr in model units, under
        `thickness_m` of ice. Snow falls where the node already freezes or where the ice surface
        is cold enough; surface melt is the same positive-degree law as before, on the
        lapse-corrected surface; basal melt is geothermal, on warm beds only. Both switch on
        over a ramp (ICE_SNOW_RAMP_C, ICE_THAW_RAMP_C) so the rates are continuous in thickness."""
        thickness = np.clip(thickness_m, 0.0, None)
        surface_c = self.temperature_c - LAPSE_RATE_C_PER_KM * thickness / 1000.0
        if G0.physical:
            return self._physical_rates(thickness, surface_c)
        snows = np.where(self.is_frozen, 1.0, np.clip((hydrology.FREEZE_POINT_C - surface_c) / ICE_SNOW_RAMP_C, 0.0, 1.0))
        accumulation = snows * hydrology.GLACIER_ACCUMULATION_RATE * self.precipitation_mm
        melt_factor = np.clip(
            hydrology.seasonal_positive_degrees(surface_c, self.seasonal_amplitude_c) / hydrology.GLACIER_MELT_REFERENCE_DEGREES_C,
            0.0,
            hydrology.GLACIER_MELT_MAX_FACTOR,
        )
        surface_melt = hydrology.GLACIER_MELT_RATE_M_PER_MYR * melt_factor
        thawed = np.where(thickness > 0.0, np.clip(_bed_excess_c(thickness, surface_c) / ICE_THAW_RAMP_C, 0.0, 1.0), 0.0)
        return accumulation, surface_melt, thawed * GEOTHERMAL_BASAL_MELT_M_PER_MYR


def _bed_excess_c(thickness_m: np.ndarray, surface_temperature_c: np.ndarray) -> np.ndarray:
    """How far (C) steady conduction through the ice puts the bed above its pressure-melting
    point (negative: frozen)."""
    gradient = ICE_GEOTHERMAL_FLUX_W_M2 / ICE_THERMAL_CONDUCTIVITY_W_M_K
    thickness = np.asarray(thickness_m, dtype=float)
    return surface_temperature_c + gradient * thickness + ICE_PRESSURE_MELT_C_PER_M * thickness


def warm_based(thickness_m: np.ndarray, surface_temperature_c: np.ndarray) -> np.ndarray:
    """True where steady conduction through `thickness_m` of ice brings the bed to the
    pressure-melting point."""
    if G0.physical:
        # D1's single thermal state: thawed where the net basal heat, with last step's
        # frictional heat, is >= 0. Called by step_ice on whole-cloud arrays only.
        frictional = G0.frictional_heat_w_m2
        if frictional is None or len(frictional) != len(np.atleast_1d(thickness_m)):
            frictional = 0.0
        return basal_heat_w_m2(thickness_m, surface_temperature_c, frictional) >= 0.0
    return (np.asarray(thickness_m) > 0.0) & (_bed_excess_c(thickness_m, surface_temperature_c) >= 0.0)


def yield_stress_pa(mobile_cover_m: np.ndarray, is_warm: np.ndarray) -> np.ndarray:
    """Basal yield stress: soft only where a thawed bed has sediment to deform."""
    if G0.fixed_yield:
        return np.full(np.shape(mobile_cover_m), ICE_YIELD_STRESS_HARD_PA)
    softness = np.where(is_warm, np.clip(np.asarray(mobile_cover_m, dtype=float) / ICE_SOFT_BED_COVER_M, 0.0, 1.0), 0.0)
    return ICE_YIELD_STRESS_HARD_PA + (ICE_YIELD_STRESS_SOFT_PA - ICE_YIELD_STRESS_HARD_PA) * softness


def flotation_thickness_m(bed_m: np.ndarray, sea_level_m: float) -> np.ndarray:
    """Ice thickness at which a column on `bed_m` just floats."""
    return _FLOTATION_RATIO * np.clip(sea_level_m - np.asarray(bed_m, dtype=float), 0.0, None)


def is_grounded(depth_m: np.ndarray, unloaded_bed_m: np.ndarray, sea_level_m: float) -> np.ndarray:
    """Grounded ice, by exactly the test `lithosphere.grounded_ice_load_kg_m2` loads with: the
    ice outweighs the water it displaces on the unloaded bed."""
    water_depth = np.clip(sea_level_m - unloaded_bed_m, 0.0, None)
    return RHO_ICE * np.clip(depth_m, 0.0, None) - RHO_WATER * water_depth > 0.0


def symmetric_graph(points: np.ndarray, neighbor_idx: np.ndarray) -> csr_matrix:
    """The k-NN graph made symmetric (k-NN isn't), weighted by great-circle edge length in
    metres."""
    n, k = neighbor_idx.shape
    rows = np.repeat(np.arange(n), k)
    cols = neighbor_idx.ravel()
    keep = rows != cols
    rows, cols = rows[keep], cols[keep]
    lo, hi = np.minimum(rows, cols), np.maximum(rows, cols)
    key = np.unique(lo.astype(np.int64) * n + hi)
    lo, hi = key // n, key % n
    length = np.maximum(geometry.angular_distance(points[lo], points[hi]) * PLANET_RADIUS_KM * 1000.0, 1.0)
    graph = coo_matrix(
        (np.concatenate([length, length]), (np.concatenate([lo, hi]), np.concatenate([hi, lo]))), shape=(n, n)
    )
    return graph.tocsr()


def _restrict(graph: csr_matrix, keep: np.ndarray) -> csr_matrix:
    """`graph` with every edge touching a node outside `keep` removed (same shape)."""
    coo = graph.tocoo()
    ok = keep[coo.row] & keep[coo.col]
    return csr_matrix((coo.data[ok], (coo.row[ok], coo.col[ok])), shape=graph.shape)


@njit(cache=True)
def _upwind(s_up: float, b_up: float, b_i: float, c_i: float, step_m: float) -> float:
    """Surface at node i from an accepted upwind point at distance `step_m`. Nye's relation
    with the edge's mean thickness, ((H_i + H_up) / 2) * (s_i - s_up) = c * step, which is
    exact along a straight flat-bed path (s^2 grows by 2 c dx). The one-point root of
    issue #353 uses H_i alone and under-reads the first cells off a margin by up to ~30%.

    Clamping the discriminant's linear term at zero keeps the update monotone (ds_i/ds_up in
    [0, 1)) and causal (s_i > s_up for step > 0), and it never divides by H, so an H = 0 margin
    is regular. The result is never below the bed: a bump that pokes through stays ice-free."""
    x = max(0.0, 2.0 * s_up - b_i - b_up)
    s = 0.5 * ((b_i + b_up) + np.sqrt(x * x + 8.0 * c_i * step_m))
    return max(s, b_i)


@njit(cache=True)
def _margin_point(i, e, j, phi, nbr_phi, nbr_bed, nbr_cap, bed, sea_level, d):
    """Where frontier edge e (active i, inactive j) puts the margin: fraction t of the edge
    from i, found by interpolating phi linearly to `d`, so the margin moves continuously with
    d. An edge toward a node the body can never occupy (a wall, or a domain node dropped as
    floating, marked by phi = -inf) stops at mid-edge. Returns (t, bed there, thickness
    there): zero, or the flotation thickness where the bed is below sea level (a grounding
    line)."""
    if j >= 0 and phi[j] <= d:
        t = 0.5
    elif nbr_phi[e] > phi[i]:
        t = (d - phi[i]) / (nbr_phi[e] - phi[i])
    else:
        t = 1.0
    t = min(max(t, 0.0), nbr_cap[e])
    b_m = bed[i] + t * (nbr_bed[e] - bed[i])
    return t, b_m, _FLOTATION_RATIO * max(0.0, sea_level - b_m)


@njit(cache=True)
def _segment_mean(h0, h1):
    """Mean thickness along a segment whose H^2 runs linearly from h0^2 to h1^2 -- exactly
    Nye's flat-bed profile, so margin cells hold the right volume at 60-180 km spacing."""
    a = h0 * h0
    b = h1 * h1
    if abs(b - a) <= 1e-9 * (a + b) or a + b == 0.0:
        return 0.5 * (h0 + h1)
    return (2.0 / 3.0) * (b * np.sqrt(b) - a * np.sqrt(a)) / (b - a)


@njit(cache=True)
def _plastic_sweep(indptr, nbr, length, nbr_phi, nbr_bed, nbr_cap, phi, bed, c, sea_level, d, active, max_thickness):
    """Ordered upwind sweep for the plastic surface over the `active` nodes of one domain.

    Frontier edges seed virtual margin points (`_margin_point`); nodes are then accepted in
    increasing surface order, each taking the lowest (earliest-arrival) update from its
    accepted neighbours.

    Cell contents come from the edges, not the node alone: each edge is split at its midpoint
    and each half belongs to the cell at that end, including the stretch of a frontier edge
    the margin has pushed past mid-edge, which lies over the outside node's cell. Along an
    edge H^2 is taken as linear. A cell's mean depth is the average over its edges of the ice
    on its halves, and its cover the average covered share. A node entering the extent takes
    over exactly the ice already on its half-edges, so the volume moves continuously with d.

    Thickness is capped at `max_thickness` (GLACIER_MAX_DEPTH_M) after the sweep: a sheet too
    big for the cap gets a flat top, and the volume it can't stack spreads its margin outward
    through the root find instead of being shed. The cap is a min of a monotone field, so V(d)
    keeps its order and continuity.

    Returns (surface, mean depth, cover)."""
    n = len(phi)
    surface = np.full(n, np.nan)
    tentative = np.full(n, np.inf)
    done = np.zeros(n, dtype=np.bool_)
    heap = [(0.0, np.int64(0))]
    heap.pop()
    for i in range(n):
        if not active[i]:
            continue
        best = np.inf
        for e in range(indptr[i], indptr[i + 1]):
            j = nbr[e]
            if j >= 0 and active[j]:
                continue
            t, b_m, h_m = _margin_point(i, e, j, phi, nbr_phi, nbr_bed, nbr_cap, bed, sea_level, d)
            best = min(best, _upwind(b_m + h_m, b_m, bed[i], c[i], t * length[e]))
        if best < np.inf:
            tentative[i] = best
            heapq.heappush(heap, (best, np.int64(i)))
    while len(heap) > 0:
        s, i = heapq.heappop(heap)
        if done[i] or s > tentative[i]:
            continue
        done[i] = True
        surface[i] = s
        for e in range(indptr[i], indptr[i + 1]):
            j = nbr[e]
            if j < 0 or not active[j] or done[j]:
                continue
            candidate = _upwind(s, bed[i], bed[j], c[j], length[e])
            if candidate < tentative[j]:
                tentative[j] = candidate
                heapq.heappush(heap, (candidate, np.int64(j)))
    for i in range(n):
        if active[i] and not done[i]:
            surface[i] = bed[i]  # an active patch with no margin at all: no profile to build
        elif active[i]:
            surface[i] = min(surface[i], bed[i] + max_thickness)
    depth = np.zeros(n)
    cover = np.zeros(n)
    for i in range(n):
        if not active[i]:
            continue
        h_i = max(surface[i] - bed[i], 0.0)
        degree = indptr[i + 1] - indptr[i]
        if degree == 0:
            depth[i] += h_i
            cover[i] += 1.0
            continue
        for e in range(indptr[i], indptr[i + 1]):
            j = nbr[e]
            if j >= 0 and active[j]:
                h_j = max(surface[j] - bed[j], 0.0)
                depth[i] += _segment_mean(h_i, np.sqrt(0.5 * (h_i * h_i + h_j * h_j))) / degree
                cover[i] += 1.0 / degree
                continue
            t, _, h_m = _margin_point(i, e, j, phi, nbr_phi, nbr_bed, nbr_cap, bed, sea_level, d)
            if t <= 0.5:
                depth[i] += 2.0 * t * _segment_mean(h_i, h_m) / degree
                cover[i] += 2.0 * t / degree
                continue
            # Past mid-edge: the far half (out to the margin) is over j's cell.
            h_half = np.sqrt(max(h_i * h_i + (0.5 / t) * (h_m * h_m - h_i * h_i), 0.0))
            depth[i] += _segment_mean(h_i, h_half) / degree
            cover[i] += 1.0 / degree
            degree_j = indptr[j + 1] - indptr[j]
            depth[j] += 2.0 * (t - 0.5) * _segment_mean(h_half, h_m) / degree_j
            cover[j] += 2.0 * (t - 0.5) / degree_j
    return surface, depth, cover


@njit(cache=True)
def _balance_flux(order, indptr, indices, length, surface, is_ice, source):
    """Accumulate `source` down the ice surface with multiple-direction weights proportional
    to each lower neighbour's surface slope (Le Brocq et al. 2006: single-direction routing
    streaks balance flux into lines). `order` lists ice nodes highest surface first, so each
    node is complete before it passes flux on. Negative accumulated flux is clipped to zero;
    flux leaving the ice (to a node not in `is_ice`) leaves the sum."""
    flux = np.zeros(len(surface))
    for k in range(len(order)):
        i = order[k]
        q = max(flux[i] + source[i], 0.0)
        flux[i] = q
        total = 0.0
        for e in range(indptr[i], indptr[i + 1]):
            j = indices[e]
            if is_ice[j] and surface[j] < surface[i]:
                total += (surface[i] - surface[j]) / length[e]
        if total <= 0.0 or q <= 0.0:
            continue
        for e in range(indptr[i], indptr[i + 1]):
            j = indices[e]
            if is_ice[j] and surface[j] < surface[i]:
                flux[j] += q * ((surface[i] - surface[j]) / length[e]) / total
    return flux


@dataclass
class BodyBalance:
    """One ice body's volume budget for the step (m^3). The contract is
    closing = opening + accumulation - surface_melt - basal_melt - calved - overflow - capped.
    `overflow_m3` is what a saturated body can't hold; `capped_m3` is ice shed at the
    GLACIER_MAX_DEPTH_M guard (see `IceLedger`)."""

    nodes: np.ndarray
    opening_m3: float
    accumulation_m3: float
    surface_melt_m3: float
    basal_melt_m3: float
    calved_m3: float
    overflow_m3: float
    closing_m3: float
    capped_m3: float = 0.0
    margin_m: float = 0.0  # the solved d (signed graph distance of the margin from S), metres
    saturated: bool = False
    # Every node whose ice this body sets: its own nodes, plus the ground it can advance onto.
    domain_nodes: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))

    @property
    def residual_m3(self) -> float:
        return (
            self.opening_m3 + self.accumulation_m3 - self.surface_melt_m3 - self.basal_melt_m3
            - self.calved_m3 - self.overflow_m3 - self.capped_m3 - self.closing_m3
        )


@dataclass
class IceLedger:
    """The world's glacier-ice budget for one step (m^3 of ice). `opening_m3` is last step's
    glacier_depth and `closing_m3` the ice after the solve, which is next step's opening.
    `snow_melted_m3` is snowfall on ice-free ground that melts within the step and never joins
    the ice.

    `frozen_lake_m3`, `frozen_river_ice_m3` and `frozen_river_capped_m3` are always 0: liquid
    water no longer freezes into the ice (RFC 0014 D2(b), #539). They stay so older code and
    caches reading them still work, and RFC 0015 D5 deletes them.

    `capped_m3` is ice shed as meltwater at the GLACIER_MAX_DEPTH_M guard, and
    `capped_nodes` how many nodes reached it (including sheets whose plastic profile was
    flattened there, which sheds nothing). Real ice spreads and melts at its base long before
    5 km, so either being nonzero means some sink or flow is missing for that ice, and it is
    logged as a warning (`warn_if_capped`)."""

    opening_m3: float = 0.0
    frozen_lake_m3: float = 0.0
    accumulation_m3: float = 0.0
    surface_melt_m3: float = 0.0
    basal_melt_m3: float = 0.0
    calved_m3: float = 0.0
    overflow_m3: float = 0.0
    closing_m3: float = 0.0
    frozen_river_ice_m3: float = 0.0
    snow_melted_m3: float = 0.0
    capped_m3: float = 0.0
    capped_nodes: int = 0
    frozen_river_capped_m3: float = 0.0

    @property
    def meltwater_m3(self) -> float:
        """Ice that left as liquid water this step, the volume `step_ice`'s meltwater carries:
        surface and basal melt, overflow, ice shed at the depth guard, and snow that melted on
        bare ground. Calved ice isn't included because it goes to the ocean unrouted. See
        hydrology's "Water units" (issue #355)."""
        return self.surface_melt_m3 + self.basal_melt_m3 + self.overflow_m3 + self.capped_m3 + self.snow_melted_m3

    @property
    def residual_m3(self) -> float:
        return (
            self.opening_m3 + self.frozen_lake_m3 + self.accumulation_m3 - self.surface_melt_m3 - self.basal_melt_m3
            - self.calved_m3 - self.overflow_m3 - self.capped_m3 - self.closing_m3
        )


def warn_if_capped(ledger: IceLedger) -> None:
    """Log a warning when any ice reached GLACIER_MAX_DEPTH_M in the solve (see `IceLedger`)."""
    if ledger.capped_nodes == 0 and ledger.capped_m3 <= 0.0:
        return
    _logger.warning(
        "glacier ice reached GLACIER_MAX_DEPTH_M (%.0f m) at %d nodes, %.3g km^3 shed as meltwater: "
        "a missing sink or flow for ice this thick",
        hydrology.GLACIER_MAX_DEPTH_M, ledger.capped_nodes, ledger.capped_m3 / 1e9,
    )


@dataclass
class IceStep:
    """One glacier step's diagnostics, aligned with hydrology's node arrays and kept on
    `HydrologyFields.ice` (so in saves: the new depth and meltwater, which hydrology already
    carries elsewhere, aren't repeated here).

    `thickness` is the plastic profile's ice thickness at each node and `covered` the cell's
    ice-covered fraction. `calved_m3` is ice that left for the ocean, unrouted. The sliding
    speed is m/yr in model units (a calibrated proxy while mass balance stays in model units)."""

    thickness: np.ndarray
    covered: np.ndarray
    warm_based: np.ndarray
    balance_flux_m3_per_yr: np.ndarray
    sliding_speed_m_per_yr: np.ndarray
    calved_m3: np.ndarray
    ledger: IceLedger
    bodies: list[BodyBalance] = field(default_factory=list)


@dataclass
class _Domain:
    """One body's candidate nodes, as a local graph. Edges to nodes outside the domain are
    kept as walls (`nbr` = -1)."""

    nodes: np.ndarray
    indptr: np.ndarray
    nbr: np.ndarray
    length: np.ndarray
    nbr_phi: np.ndarray
    nbr_bed: np.ndarray
    nbr_cap: np.ndarray
    nbr_ocean: np.ndarray
    phi: np.ndarray
    bed: np.ndarray
    c: np.ndarray
    area: np.ndarray
    water_depth: np.ndarray
    climate: SurfaceClimate
    bare_rates: tuple[np.ndarray, np.ndarray, np.ndarray]


def _build_domain(graph, nodes, phi, bed, c, area, sea_level_m, is_ocean, climate) -> _Domain:
    local = np.full(graph.shape[0], -1, dtype=np.int64)
    local[nodes] = np.arange(len(nodes))
    sub = graph[nodes]
    cols = sub.indices
    nbr = local[cols]
    rows = np.repeat(np.arange(len(nodes)), np.diff(sub.indptr))
    wall = nbr < 0
    length = sub.data.astype(float)
    nbr_phi = np.where(wall, phi[nodes][rows] + length, phi[cols])
    sub_climate = climate.subset(nodes)
    return _Domain(
        nodes=nodes,
        indptr=sub.indptr.astype(np.int64),
        nbr=nbr,
        length=length,
        nbr_phi=nbr_phi,
        nbr_bed=bed[cols].astype(float),
        nbr_cap=np.where(wall, 0.5, 1.0),
        nbr_ocean=is_ocean[cols],
        phi=phi[nodes],
        bed=bed[nodes],
        c=c[nodes],
        area=area[nodes],
        water_depth=np.clip(sea_level_m - bed[nodes], 0.0, None),
        climate=sub_climate,
        bare_rates=sub_climate.rates(np.zeros(len(nodes))),
    )


@dataclass
class _Shape:
    """One candidate sheet: the extent phi <= d, its surface and node thickness, and each
    cell's mean depth (what glacier_depth stores) and covered fraction."""

    d: float
    active: np.ndarray
    surface: np.ndarray
    thickness: np.ndarray
    depth: np.ndarray
    cover: np.ndarray
    volume: float
    rates: tuple[np.ndarray, np.ndarray, np.ndarray]

    def scaled(self, factor: float) -> "_Shape":
        return _Shape(
            self.d, self.active, self.surface, self.thickness * factor, self.depth * factor, self.cover,
            self.volume * factor, self.rates,
        )

    def budget_rates(self, dom: _Domain) -> tuple[float, float, float]:
        """Area-weighted (accumulation, surface melt, basal melt), m^3/Myr: the covered area at
        the solved thickness, plus positive-balance bare ground, whose snow feeds the body."""
        acc, sm, bm = self.rates
        acc0, sm0, _ = dom.bare_rates
        covered = self.cover * dom.area
        bare = (1.0 - np.minimum(self.cover, 1.0)) * dom.area * ((acc0 - sm0) > 0.0)
        return (
            float(np.dot(covered, acc) + np.dot(bare, acc0)),
            float(np.dot(covered, sm) + np.dot(bare, sm0)),
            float(np.dot(covered, bm)),
        )


def _shape(dom: _Domain, d: float, sea_level_m: float) -> _Shape:
    """The plastic sheet filling the extent phi <= d, less marine cells it can't ground. A
    dropped cell becomes a wall, which can thin its neighbours enough to float them too, so
    this repeats; anything still afloat after the last pass is left bare."""
    active = dom.phi <= d
    phi = dom.phi.copy()
    for attempt in range(ICE_FLOTATION_PASSES + 1):
        surface, depth, cover = _plastic_sweep(
            dom.indptr, dom.nbr, dom.length, dom.nbr_phi, dom.nbr_bed, dom.nbr_cap, phi, dom.bed, dom.c, sea_level_m, d, active,
            hydrology.GLACIER_MAX_DEPTH_M,
        )
        # Any cell holding ice it can't ground -- in the extent, or outside it under the margin
        # past mid-edge -- is dropped: it becomes a wall the ice stops half-way to.
        floating = (depth > 0.0) & (dom.water_depth > 0.0) & ~is_grounded(depth, dom.bed, sea_level_m)
        if not np.any(floating):
            break
        active = active & ~floating
        if attempt < ICE_FLOTATION_PASSES:
            phi[floating] = -np.inf
        else:
            depth = np.where(floating, 0.0, depth)
            cover = np.where(floating, 0.0, cover)
    # Outside nodes can hold the margin ice past mid-edge; dropped (floating) ones hold none.
    holds = active | ((cover > 0.0) & (phi > -np.inf))
    depth = np.where(holds, depth, 0.0)
    cover = np.where(holds, cover, 0.0)
    # The mass balance reads each cell's mean thickness where it's covered, not the node's: a
    # node joining the extent has zero thickness at its own margin, but the ice already over
    # its cell doesn't change, so the balance stays continuous in d.
    mean_thickness = np.divide(depth, cover, out=np.zeros_like(depth), where=cover > 0.0)
    thickness = np.where(active, np.clip(surface - dom.bed, 0.0, None), mean_thickness)
    surface = np.where(active, surface, np.nan)
    volume = float(np.dot(depth, dom.area))
    return _Shape(d, active, surface, thickness, depth, cover, volume, dom.climate.rates(mean_thickness))


def _bracketed_root(f, lo: float, hi: float, f_lo: float, f_hi: float, tolerance: float):
    """Illinois false position on a bracket with f(lo) <= 0 < f(hi). Returns the final
    (lo, hi): f may jump (a node freezing over, a cell grounding), so the caller picks a side
    rather than trusting a point root."""
    side = 0
    for _ in range(ICE_ROOT_MAX_ITERATIONS):
        if hi - lo <= tolerance:
            break
        x = hi - f_hi * (hi - lo) / (f_hi - f_lo) if f_hi != f_lo else 0.5 * (lo + hi)
        if not lo < x < hi:
            x = 0.5 * (lo + hi)
        fx = f(x)
        if fx <= 0.0:
            lo, f_lo = x, fx
            if side == -1:
                f_hi *= 0.5
            side = -1
        else:
            hi, f_hi = x, fx
            if side == 1:
                f_lo *= 0.5
            side = 1
    return lo, hi


def _solve_plastic_body(dom: _Domain, opening_m3: float, dt_myr: float, sea_level_m: float, step_m: float):
    """Backward-Euler volume update for one body. Returns (shape, accumulation, surface melt,
    basal melt, saturation excess, saturated) with volumes in m^3, the shape already scaled to
    the closing volume.

    Last step's extent is the one whose profile holds the opening volume, V(d_prev) = V_open
    (V is nondecreasing in d, so bisection finds it). Where the ice actually sits is a poor
    guide: snow builds ice on any cold node, often far beyond a body's sheet (and before #539,
    frozen river and lake ice joined a body wherever they froze)."""
    cache: dict[float, _Shape] = {}

    def shape(d: float) -> _Shape:
        if d not in cache:
            cache[d] = _shape(dom, d, sea_level_m)
        return cache[d]

    def f(d: float) -> float:
        s = shape(d)
        acc, sm, bm = s.budget_rates(dom)
        return s.volume - opening_m3 - dt_myr * (acc - sm - bm)

    # The extent runs from empty (d_min) to the whole domain (d_max) -- less its farthest node
    # when the domain has no wall, since with every node in nothing would seed a margin.
    d_min = float(dom.phi.min()) - 1.0
    d_max = max(float(dom.phi.max()) - (0.0 if np.any(dom.nbr < 0) else 1.0), d_min)
    d0 = _volume_extent(shape, opening_m3, d_min, d_max, step_m)
    f0 = f(d0)
    saturated = False
    # Walk outward from last step's extent one edge at a time, so the bracket holds the root
    # nearest it: B(d) can change sign more than once (a big sheet's cold, high surface), and
    # a growing step could jump a pair of roots.
    if f0 <= 0.0:
        lo, f_lo = d0, f0
        while True:
            hi = min(lo + step_m, d_max)
            f_hi = f(hi)
            if f_hi > 0.0:
                break
            lo, f_lo = hi, f_hi
            if hi >= d_max:
                saturated = True
                break
    else:
        hi, f_hi = d0, f0
        while True:
            lo = max(hi - step_m, d_min)
            f_lo = f(lo)
            if f_lo <= 0.0:
                break
            hi, f_hi = lo, f_lo
    if G0.enumerate_roots or G0.root_rule != "nearest":
        lo, hi, f_lo, f_hi, saturated = _g0_choose_root(
            shape, f, opening_m3, dt_myr, d_min, d_max, step_m, (lo, hi, f_lo, f_hi, saturated), len(dom.nodes)
        )
    if saturated:
        final = shape(d_max)
    else:
        lo, hi = _bracketed_root(f, lo, hi, f_lo, f_hi, ICE_ROOT_TOLERANCE_M)
        # The closing volume is reached by scaling the nearer side's profile (f can still jump
        # a little, where a cell grounds or a node joins with a different area per edge).
        final = shape(hi) if shape(lo).volume <= 0.0 or abs(f(hi)) < abs(f(lo)) else shape(lo)
    acc, sm, bm = (dt_myr * r for r in final.budget_rates(dom))
    target = opening_m3 + acc - sm - bm
    excess = 0.0
    if saturated:
        excess = max(target - final.volume, 0.0)
        target -= excess
    if target < 0.0:
        # The scaled-to side melted more than existed: melt only what there was.
        shortfall = -target
        sm_share = sm / max(sm + bm, 1e-30)
        sm -= shortfall * sm_share
        bm -= shortfall * (1.0 - sm_share)
        target = 0.0
    if final.volume > 0.0:
        final = final.scaled(target / final.volume)
    return final, acc, sm, bm, excess, saturated


def _g0_choose_root(shape, f, opening_m3, dt_myr, d_min, d_max, step_m, walked, n_nodes):
    """G0 prototype (D2's root choice): every sign change of the residual on [d_min, d_max]
    at one-edge spacing, logged with which are stable (an upward crossing where the balance
    falls as the volume grows), and the bracket the root rule keeps. `walked` is today's
    walk's (lo, hi, f_lo, f_hi, saturated)."""
    grid = np.unique(np.append(np.arange(d_min, d_max, max(step_m, 1.0)), d_max))
    fs = np.array([f(d) for d in grid])
    vols = np.array([shape(d).volume for d in grid])
    balance = (vols - opening_m3 - fs) / dt_myr
    roots = []
    for k in range(len(grid) - 1):
        a, b = fs[k], fs[k + 1]
        if a <= 0.0 < b:
            roots.append((grid[k], grid[k + 1], a, b, bool(balance[k + 1] <= balance[k])))
        elif a > 0.0 >= b:
            roots.append((grid[k], grid[k + 1], b, a, False))
    stable = [r for r in roots if r[4]]
    lo, hi, f_lo, f_hi, saturated = walked
    walk_rank = next((i for i, r in enumerate(stable) if r[0] <= hi and lo <= r[1]), -1)
    kept = "walk"
    if G0.root_rule != "nearest" and stable and not saturated:
        r = stable[-1] if G0.root_rule == "largest" else stable[0]
        lo, hi, f_lo, f_hi = r[0], r[1], r[2], r[3]
        kept = G0.root_rule
    G0_ROOT_LOG.append({
        "nodes": n_nodes, "opening_km3": opening_m3 / 1e9, "roots": len(roots), "stable": len(stable),
        "walk_rank": walk_rank, "kept": kept, "saturated": bool(saturated),
        "stable_volumes_km3": [round(float(shape(r[1]).volume) / 1e9, 3) for r in stable][:8],
    })
    return lo, hi, f_lo, f_hi, saturated


def _volume_extent(shape, volume_m3: float, d_min: float, d_max: float, tolerance: float) -> float:
    """The smallest extent d whose profile holds `volume_m3` (or d_max if none does), to
    within `tolerance`."""
    if volume_m3 <= 0.0:
        return d_min
    if shape(d_max).volume < volume_m3:
        return d_max
    lo, hi = d_min, d_max
    while hi - lo > tolerance:
        mid = 0.5 * (lo + hi)
        if shape(mid).volume < volume_m3:
            lo = mid
        else:
            hi = mid
    return hi


def _solve_uniform_body(nodes, opening_m3, area, climate: SurfaceClimate, dt_myr):
    """A body too small for a plastic profile: one thickness over all its nodes, with the
    same backward-Euler volume update. Returns (thickness, accumulation, surface melt, basal
    melt, saturation excess), volumes in m^3."""
    total_area = float(area.sum())

    def budget(h: float):
        acc, sm, bm = climate.rates(np.full(len(nodes), h))
        return dt_myr * float(np.dot(area, acc)), dt_myr * float(np.dot(area, sm)), dt_myr * float(np.dot(area, bm))

    def g(h: float) -> float:
        acc, sm, bm = budget(h)
        return h * total_area - opening_m3 - (acc - sm - bm)

    g0 = g(0.0)
    if g0 >= 0.0:
        acc, sm, bm = budget(0.0)
        # Everything melts: the melt is whatever existed plus what fell.
        return 0.0, acc, opening_m3 + acc, 0.0, 0.0
    h_max = hydrology.GLACIER_MAX_DEPTH_M
    if g(h_max) < 0.0:
        acc, sm, bm = budget(h_max)
        excess = opening_m3 + acc - sm - bm - h_max * total_area
        return h_max, acc, sm, bm, excess
    lo, hi = _bracketed_root(g, 0.0, h_max, g0, g(h_max), 1e-3)
    acc, sm, bm = budget(hi)
    target = opening_m3 + acc - sm - bm
    if target < 0.0:
        acc, sm, bm = budget(lo)
        target = max(opening_m3 + acc - sm - bm, 0.0)
    return target / total_area, acc, sm, bm, 0.0


def _signed_distance(graph, eligible_graph, positive: np.ndarray):
    """phi = (distance to S - distance to not-S) / 2 over graph edges: negative inside the
    positive-balance set S, positive outside, zero half-way across each boundary edge. The
    outside distance runs over the eligible graph only. Also returns, per node, the S node
    nearest to it there (-1 if unreachable)."""
    n = graph.shape[0]
    s_idx = np.flatnonzero(positive)
    d_out, _, sources = dijkstra(eligible_graph, indices=s_idx, min_only=True, return_predecessors=True)
    not_s = np.flatnonzero(~positive)
    if len(not_s):
        d_in = dijkstra(graph, indices=not_s, min_only=True)
    else:
        d_in = np.full(n, np.inf)
    d_in = np.where(positive, d_in, 0.0)
    d_out = np.where(positive, 0.0, d_out)
    return 0.5 * (d_out - d_in), sources


def step_ice(
    points: np.ndarray,
    neighbor_idx: np.ndarray,
    area_m2: np.ndarray,
    elevation: np.ndarray,
    prior_deflection_m: np.ndarray,
    sea_level_m: float,
    prev_depth: np.ndarray,
    climate: SurfaceClimate,
    mobile_cover_m: np.ndarray,
    is_ocean: np.ndarray,
    allow_ocean_ice: np.ndarray,
    years: float,
) -> tuple[np.ndarray, np.ndarray, IceStep]:
    """One glacier step over the whole node cloud. See the module docstring for the method.

    Returns (new glacier_depth, meltwater, diagnostics). glacier_depth is the mean ice depth
    over each cell, so depth * area is its volume. Meltwater -- surface and basal melt, snow
    that melts on bare ground, and overflow -- is in metres of ice over the cell, for
    hydrology's water_source; calved ice is not in it.

    Bed convention (issue #353 section 4): the thickness solve and the flotation test both
    use the unloaded bed, elevation minus last step's `ice_load_deflection_m` -- the bed
    `lithosphere.grounded_ice_load_kg_m2` loads -- and stored ice is grounded exactly when
    that function gives it a load. The mass-balance surface is the current (already
    depressed) elevation plus the ice, a one-step lag on the deflection.

    The ice-age polar sea-ice cap (`allow_ocean_ice`) keeps its own per-node rule (accumulate,
    melt, never flow) and sits outside the bodies. Other ice on the ocean that isn't grounded
    calves."""
    n = len(points)
    area = np.asarray(area_m2, dtype=float)
    dt_myr = years / 1_000_000.0
    bed = elevation - prior_deflection_m
    opening_depth = np.clip(prev_depth, 0.0, None)
    ledger = IceLedger(opening_m3=float(np.dot(opening_depth, area)))
    h0 = opening_depth
    out = _Output(n)

    # The sea-ice cap: the old per-node rule, minus flow.
    cap = allow_ocean_ice & is_ocean
    if np.any(cap):
        acc, sm, _ = climate.subset(cap).rates(np.zeros(int(cap.sum())))
        gained = acc * dt_myr
        lost = np.minimum(sm * dt_myr, h0[cap] + gained)
        depth = h0[cap] + gained - lost
        # The same depth guard as the bodies (review of #359): a sustained cold spell with
        # little summer melt would otherwise stack the floating cap past it.
        over = np.clip(depth - hydrology.GLACIER_MAX_DEPTH_M, 0.0, None)
        out.depth[cap] = depth - over
        out.melt_m3[cap] += (lost + over) * area[cap]
        ledger.accumulation_m3 += float(np.dot(gained, area[cap]))
        ledger.surface_melt_m3 += float(np.dot(lost, area[cap]))
        ledger.capped_m3 += float(np.dot(over, area[cap]))
        ledger.capped_nodes += int(np.count_nonzero(over > 0.0))

    # Nodes a body may occupy: never the sea-ice cap, never ocean too deep to ground the
    # thickest ice allowed. Ice already afloat anywhere else calves.
    eligible = ~cap & (flotation_thickness_m(bed, sea_level_m) <= hydrology.GLACIER_MAX_DEPTH_M)
    stray = ~cap & (h0 > 0.0) & (~eligible | ~is_grounded(h0, bed, sea_level_m))
    out.calved[stray] = h0[stray] * area[stray]
    h0 = np.where(stray, 0.0, h0)

    bare_acc = climate.rates(np.zeros(n))[0]
    base_out = out

    def solve_bodies(ref_h: np.ndarray):
        """Membership and yield stress read at `ref_h` (the opening ice; with G0's arm (f),
        the last solve's thickness). Returns the filled output, bodies, domain mask, a scratch
        ledger and the membership used."""
        trial = _Output.copy_of(base_out)
        scratch = IceLedger()
        ref_acc, ref_sm, ref_bm = climate.rates(ref_h)
        reference_balance = ref_acc - ref_sm - ref_bm
        positive = eligible & (reference_balance > 0.0) & (~is_ocean | (h0 > 0.0))
        members = eligible & ((h0 > 0.0) | positive)

        bodies: list[BodyBalance] = []
        in_domain = np.zeros(n, dtype=bool)
        if np.any(members):
            graph = symmetric_graph(points, neighbor_idx)
            _, comp = connected_components(_restrict(graph, members), directed=False)
            comp = np.where(members, comp, -1)
            body_ids = np.unique(comp[members])
            # A body with no positive-balance node (a dying one) shrinks toward its
            # least-negative node instead. One that is positive everywhere and borders nothing
            # else (it fills its whole patch of the graph) has no margin to measure from, so its
            # least-positive node, lowest first, stands in for one.
            seeds = positive.copy()
            touches_outside = graph @ (~positive).astype(float) > 0.0
            for b in body_ids:
                idx = np.flatnonzero(comp == b)
                if not np.any(positive[idx]):
                    seeds[idx[np.argmax(reference_balance[idx])]] = True
                elif np.all(positive[idx]) and not np.any(touches_outside[idx]):
                    seeds[idx[np.lexsort((bed[idx], reference_balance[idx]))[0]]] = False
            phi, sources = _signed_distance(graph, _restrict(graph, eligible), seeds)
            # Each eligible node belongs to the body whose positive-balance area is nearest; a
            # body's own nodes always belong to it.
            label = np.where(eligible & (sources >= 0), comp[np.clip(sources, 0, None)], -1)
            label[members] = comp[members]
            # Yield stress from last step's thermal state (a one-step lag: the new thickness
            # isn't known until the profile it shapes is solved).
            is_warm0 = warm_based(ref_h, climate.temperature_c - LAPSE_RATE_C_PER_KM * ref_h / 1000.0)
            c = yield_stress_pa(mobile_cover_m, is_warm0) / (RHO_ICE * GRAVITY_M_S2)
            step_m = float(np.median(graph.data)) if graph.nnz else 1.0

            order = np.argsort(label, kind="stable")
            sorted_labels = label[order]
            starts = np.searchsorted(sorted_labels, body_ids, side="left")
            ends = np.searchsorted(sorted_labels, body_ids, side="right")
            for b, start, end in zip(body_ids, starts, ends):
                dom_nodes = order[start:end]
                body_nodes = np.flatnonzero(comp == b)
                in_domain[dom_nodes] = True
                opening = float(np.dot(h0[body_nodes], area[body_nodes]))
                if len(dom_nodes) < ICE_BODY_MIN_NODES:
                    body = _uniform_body(trial, body_nodes, opening, area, bed, climate, is_ocean, dt_myr)
                    rest = np.setdiff1d(dom_nodes, body_nodes)
                    scratch.snow_melted_m3 += _melt_bare_snow(trial, rest, bare_acc, area, dt_myr)
                else:
                    dom = _build_domain(graph, dom_nodes, phi, bed, c, area, sea_level_m, is_ocean, climate)
                    body = _plastic_body(trial, dom, body_nodes, opening, area, dt_myr, sea_level_m, step_m, scratch)
                _guard_body(trial, body, bed, area, sea_level_m)
                bodies.append(body)
        return trial, bodies, in_domain, scratch, (members, is_warm0 if np.any(members) else None)

    out, bodies, in_domain, scratch, used = solve_bodies(h0)
    if G0.reevaluate_membership:
        # Arm (f): re-read membership and yield stress at the solved thickness, at most
        # G0_MEMBERSHIP_ITERATIONS times; keep the last iterate if it hasn't settled.
        converged = False
        for iteration in range(G0_MEMBERSHIP_ITERATIONS):
            trial = solve_bodies(out.thickness)
            same = np.array_equal(trial[4][0], used[0]) and (
                (trial[4][1] is None and used[1] is None)
                or (trial[4][1] is not None and used[1] is not None and np.array_equal(trial[4][1], used[1]))
            )
            out, bodies, in_domain, scratch, used = trial
            if same:
                converged = True
                break
        G0_MEMBERSHIP_LOG.append({"iterations": iteration + 1, "converged": converged})
    ledger.snow_melted_m3 += scratch.snow_melted_m3

    # Snow on frozen ground no body claims melts in place (it can't be positive balance there,
    # or the node would be a body member).
    loose = np.flatnonzero(~in_domain & ~cap & ~is_ocean)
    ledger.snow_melted_m3 += _melt_bare_snow(out, loose, bare_acc, area, dt_myr)

    for body in bodies:
        ledger.accumulation_m3 += body.accumulation_m3
        ledger.surface_melt_m3 += body.surface_melt_m3
        ledger.basal_melt_m3 += body.basal_melt_m3
        ledger.overflow_m3 += body.overflow_m3
        ledger.capped_m3 += body.capped_m3
    ledger.calved_m3 = float(out.calved.sum())
    ledger.closing_m3 = float(np.dot(out.depth, area))
    # Body nodes at the depth guard, whether shed there or flattened by the profile's cap.
    at_cap = ~cap & (out.depth > 0.0) & (out.thickness >= hydrology.GLACIER_MAX_DEPTH_M * (1.0 - 1e-9))
    ledger.capped_nodes += int(np.count_nonzero(at_cap))
    warn_if_capped(ledger)

    thickness = out.thickness
    is_warm = warm_based(thickness, climate.temperature_c - LAPSE_RATE_C_PER_KM * thickness / 1000.0)
    flux, _, sliding = _diagnose_flow(
        points, neighbor_idx, area, out.surface, out.depth, opening_depth, thickness, out.covered, climate, is_warm, years,
    )
    melt_m = out.melt_m3 / np.where(area > 0.0, area, 1.0)
    return out.depth, melt_m, IceStep(
        thickness=thickness,
        covered=out.covered,
        warm_based=is_warm,
        balance_flux_m3_per_yr=flux,
        sliding_speed_m_per_yr=sliding,
        calved_m3=out.calved,
        ledger=ledger,
        bodies=bodies,
    )


class _Output:
    """The per-node arrays `step_ice` fills body by body."""

    def __init__(self, n: int):
        self.depth = np.zeros(n)
        self.melt_m3 = np.zeros(n)
        self.calved = np.zeros(n)
        self.thickness = np.zeros(n)
        self.covered = np.zeros(n)
        self.surface = np.full(n, np.nan)

    @classmethod
    def copy_of(cls, other: "_Output") -> "_Output":
        copy = cls(len(other.depth))
        for name in ("depth", "melt_m3", "calved", "thickness", "covered", "surface"):
            setattr(copy, name, getattr(other, name).copy())
        return copy


def _uniform_body(out: _Output, nodes, opening, area, bed, climate, is_ocean, dt_myr) -> BodyBalance:
    h, acc, sm, bm, excess = _solve_uniform_body(nodes, opening, area[nodes], climate.subset(nodes), dt_myr)
    out.thickness[nodes] = h
    out.covered[nodes] = 1.0 if h > 0.0 else 0.0
    out.depth[nodes] = h
    out.surface[nodes] = bed[nodes] + h if h > 0.0 else np.nan
    out.melt_m3[nodes] += (sm + bm) * area[nodes] / area[nodes].sum()
    body = BodyBalance(nodes, opening, acc, sm, bm, 0.0, 0.0, 0.0, saturated=excess > 0.0, domain_nodes=nodes)
    _shed_excess(out, body, excess, nodes, bool(np.any(is_ocean[nodes])), area)
    return body


def _plastic_body(out: _Output, dom: _Domain, body_nodes, opening, area, dt_myr, sea_level_m, step_m, ledger) -> BodyBalance:
    nodes = dom.nodes
    shape, acc, sm, bm, excess, saturated = _solve_plastic_body(dom, opening, dt_myr, sea_level_m, step_m)
    cov = np.minimum(shape.cover, 1.0)
    out.thickness[nodes] = shape.thickness
    out.covered[nodes] = cov
    out.depth[nodes] = shape.depth
    out.surface[nodes] = np.where(shape.active, dom.bed + shape.thickness, np.nan)
    _, r_sm, r_bm = shape.rates
    acc0, sm0, _ = dom.bare_rates
    feeds = (acc0 - sm0) > 0.0
    # Booked melt goes to the nodes in proportion to where it happens.
    where_melt = shape.cover * dom.area * (r_sm + r_bm) + (1.0 - cov) * dom.area * np.where(feeds, sm0, 0.0)
    if where_melt.sum() > 0.0:
        out.melt_m3[nodes] += (sm + bm) * where_melt / where_melt.sum()
    elif sm + bm > 0.0:
        out.melt_m3[body_nodes] += (sm + bm) * area[body_nodes] / area[body_nodes].sum()
    # Snow on bare ground that can't feed the body melts where it fell.
    bare_snow = (1.0 - cov) * dom.area * np.where(feeds, 0.0, acc0) * dt_myr
    out.melt_m3[nodes] += bare_snow
    ledger.snow_melted_m3 += float(bare_snow.sum())
    body = BodyBalance(body_nodes, opening, acc, sm, bm, 0.0, 0.0, 0.0, margin_m=shape.d, saturated=saturated, domain_nodes=nodes)
    frontier = shape.active & (cov < 1.0)
    margin = nodes[frontier] if np.any(frontier) else nodes[shape.active] if np.any(shape.active) else body_nodes
    marine = bool(np.any(shape.active & (dom.water_depth > 0.0))) or _touches_ocean(dom, shape.active)
    _shed_excess(out, body, excess, margin, marine, area)
    return body


def _shed_excess(out: _Output, body: BodyBalance, excess: float, margin, marine: bool, area) -> None:
    """Volume a saturated body can't hold calves from a marine margin, or else leaves the
    margin as meltwater (as GLACIER_MAX_DEPTH_M's overflow always has)."""
    if excess <= 0.0:
        return
    share = area[margin] / area[margin].sum()
    if marine:
        out.calved[margin] += excess * share
        body.calved_m3 += excess
    else:
        out.melt_m3[margin] += excess * share
        body.overflow_m3 += excess


def _guard_body(out: _Output, body: BodyBalance, bed, area, sea_level_m) -> None:
    """Per-node guards on a solved body, booked in its own budget: ice over
    GLACIER_MAX_DEPTH_M is shed as meltwater (`capped_m3`), and any marine cell the final (scaled) ice
    no longer grounds calves -- so stored ice is grounded exactly where
    `grounded_ice_load_kg_m2` loads it. Then the closing volume."""
    nodes = body.domain_nodes
    over = np.clip(out.depth[nodes] - hydrology.GLACIER_MAX_DEPTH_M, 0.0, None)
    if np.any(over > 0.0):
        kept = out.depth[nodes] - over
        out.thickness[nodes] *= np.divide(kept, out.depth[nodes], out=np.ones_like(kept), where=out.depth[nodes] > 0.0)
        out.depth[nodes] = kept
        out.melt_m3[nodes] += over * area[nodes]
        body.capped_m3 += float(np.dot(over, area[nodes]))
    afloat = (out.depth[nodes] > 0.0) & ~is_grounded(out.depth[nodes], bed[nodes], sea_level_m)
    if np.any(afloat):
        idx = nodes[afloat]
        out.calved[idx] += out.depth[idx] * area[idx]
        body.calved_m3 += float(np.dot(out.depth[idx], area[idx]))
        out.depth[idx] = 0.0
        out.thickness[idx] = 0.0
        out.covered[idx] = 0.0
        out.surface[idx] = np.nan
    body.closing_m3 = float(np.dot(out.depth[nodes], area[nodes]))


def _melt_bare_snow(out: _Output, idx: np.ndarray, bare_acc: np.ndarray, area: np.ndarray, dt_myr: float) -> float:
    """Snowfall on ice-free cells `idx` over the step melts where it fell. Returns m^3."""
    water = bare_acc[idx] * area[idx] * dt_myr
    out.melt_m3[idx] += water
    return float(water.sum())


def _touches_ocean(dom: _Domain, active: np.ndarray) -> bool:
    """Whether any active node of the domain borders the ocean (a marine margin)."""
    for i in np.flatnonzero(active):
        if np.any(dom.nbr_ocean[dom.indptr[i] : dom.indptr[i + 1]]):
            return True
    return False


def _diagnose_flow(points, neighbor_idx, area, surface, depth, prev_depth, thickness, covered, climate, is_warm, years):
    """Balance flux, depth-averaged velocity and basal sliding speed (issue #353 section 6).

    The source per node is (b - d(depth)/dt) * area in m^3/yr: what the ice gains there and
    doesn't keep. It is routed down the plastic surface; the velocity is flux over the cell
    width (sqrt(area)) and thickness. Diagnostic only: volume moves through the body solve, so
    clipping flux at zero conserves nothing and loses nothing."""
    n = len(surface)
    flux = np.zeros(n)
    velocity = np.zeros(n)
    sliding = np.zeros(n)
    is_ice = np.isfinite(surface) & (thickness > 0.0)
    if not np.any(is_ice) or years <= 0.0:
        return flux, velocity, sliding
    acc, sm, bm = climate.rates(thickness)
    balance_m_per_yr = covered * (acc - sm - bm) / 1.0e6
    source = np.where(is_ice, (balance_m_per_yr - (depth - prev_depth) / years) * area, 0.0)
    graph = symmetric_graph(points, neighbor_idx)
    ice_idx = np.flatnonzero(is_ice)
    order = ice_idx[np.argsort(-surface[ice_idx], kind="stable")]
    flux = _balance_flux(
        order.astype(np.int64), graph.indptr.astype(np.int64), graph.indices.astype(np.int64), graph.data.astype(float),
        np.where(is_ice, surface, -np.inf), is_ice, source,
    )
    width = np.sqrt(np.clip(area, 1.0, None))
    velocity = np.where(is_ice, flux / (width * np.maximum(thickness, ICE_VELOCITY_MIN_THICKNESS_M)), 0.0)
    sliding = np.where(is_ice & is_warm, ICE_SLIDING_FRACTION * velocity, 0.0)
    return flux, velocity, sliding

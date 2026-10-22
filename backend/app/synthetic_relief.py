"""Generation-time continental relief from a synthetic tectonic history (#218).

A freshly generated world has no tectonic history, so relief drawn from noise alone looks
unlike Earth: no mountain belts where continents once collided, no cordilleras on leading
margins, coastlines at noise-quantile or plate-boundary edges. This pass, run once on a
cratons-first world (RFC 0009) after its plates are built, draws a plausible history from
the layout and the first predicted plate motion and backs the result with crustal columns:

- **Sutures.** Neighbouring cratons on one continent are joined by a belt along their
  (noise-warped) bisector. Each pair draws an age class: young (high, with a plateau on one
  side), middle, or old and eroded; some pairs are welded with no surviving belt.
- **Active margins.** Continental crust moving toward its nearest ocean relative to the plate
  that owns it (the leading edge) gets an Andean cordillera ~200 km inland and a narrow shelf;
  trailing and passive margins (ocean on the continent's own plate) get a low rift-shoulder
  escarpment and wide shelves. Predicted continent-continent convergence adds a collision
  belt. Belt relief is halved on cratons, rising to full strength 150 km outside them.
- **Coastlines.** Land is the requested share of nodes with the highest land potential:
  distance inward from the continental-crust edge, noise-scaled, plus the mountain relief.
- **Interiors.** Cratons split into exposed shields (rolling hills) and covered platforms
  (gentle, with broad basins); other interior land is hill country or plains; high-latitude
  shields and lowlands get glacial-scour pits.
- **Submerged continental crust.** A shelf near the coast, deepening to marginal-plateau
  depth far from it, with a slope at the crust edge.

Every changed column is backed: Hc comes from `lithosphere.crustal_thickness_for_elevation`
on the column's existing Hm, and elevation is then re-derived from the column, so the
relief is a faithful readout of Hc/Hm like everything else. Hm is left at the generated
reference lid: it moves elevation ~11x less per km than Hc, and the isostasy has no
depleted craton keel to justify a thicker lid under shields (#218). Young belts start with
the Moho thermal lag of crust that was just buried (`orogeny.bury_moho`), so they heat and
evolve over the next tens of Myr instead of starting at steady state.

Oceanic crust is not changed. Age-depth bathymetry, ridges and trenches need thermal
buoyancy the column model does not represent: an 11 km oceanic column floats no shallower
than about -4.4 km, so a -2.6 km ridge can't be backed (#218).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import cratons, lithosphere, orogeny, rheology, torque
from .elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from .noise import IsotropicSphereNoise
from .plates import PLANET_RADIUS_KM

if TYPE_CHECKING:
    from .continent_layout import ContinentLayout
    from .world import World

# Off until an RFC settles its open questions (ocean relief, lake share, early evolution;
# see #218's exploration summary). Cratons-first only.
SYNTHETIC_RELIEF_DEFAULT = False
FIRST_MOTION_STEP_YEARS = 5_000_000.0
# Boundary nodes within this distance of another plate read that contact's regime; at least
# REGIME_REACH_SPACINGS node spacings, so a coarse world still finds its contacts (principle 4).
REGIME_REACH_KM = 150.0
REGIME_REACH_SPACINGS = 1.5

# Sutures. Pairs of cratons whose summed distance exceeds this are not neighbours.
SUTURE_PAIR_MAX_KM = 2_400.0
# Bisector warp amplitude: enough to bend a belt off the straight bisector between two
# cratons, which otherwise draws a ladder of straight belts across a many-craton continent.
SUTURE_WARP_KM = 300.0
# Share of neighbouring craton pairs welded long ago with no surviving belt.
SUTURE_DROP_PROBABILITY = 0.3
# (probability, peak m, Gaussian sigma km, plateau m). Young belts are Himalaya-like with a
# Tibet-like plateau on one side, middle ones Alps/Rockies-like, old ones Appalachian-like.
SUTURE_YOUNG = 0  # index of the young class in SUTURE_CLASSES
SUTURE_CLASSES = ((0.25, 4_500.0, 180.0, 3_500.0), (0.35, 2_500.0, 120.0, 0.0), (0.40, 900.0, 150.0, 0.0))
SUTURE_PLATEAU_KM = 300.0
SUTURE_PLATEAU_EDGE_KM = 150.0
# A pair's belt gives way where a third craton is nearer than either of the pair: full
# strength while the third is this much farther than the pair's nearer craton, none when it
# is this much nearer, so belts meet without a cliff where the nearest pair changes.
SUTURE_THIRD_CRATON_TAPER_KM = 150.0

# Active margins: Andean cordillera on leading edges, peak ~200 km inland (the Andes' crest
# sits 150-250 km from the trench), scaled by how directly and fast the plate advances.
CORDILLERA_PEAK_M = 4_000.0
CORDILLERA_INLAND_KM = 200.0
CORDILLERA_SIGMA_KM = 90.0
LEADING_EDGE_FULL_SPEED_CM_YR = 3.0
COLLISION_PEAK_M = 4_200.0
COLLISION_SIGMA_KM = 160.0
ESCARPMENT_PEAK_M = 500.0
ESCARPMENT_INLAND_KM = 250.0
ESCARPMENT_SIGMA_KM = 100.0
# Cratons resist shortening (principle 7): belt relief on a craton is scaled by this, rising
# smoothly to full strength over CRATON_EDGE_FADE_KM outside it (a step would leave a cliff
# at the craton edge). A plateau fades out over the same distance before the craton it faces:
# Tibet sits on accreted terranes, not on a shield.
CRATON_OROGEN_FACTOR = 0.5
CRATON_EDGE_FADE_KM = 150.0
# A passive margin (ocean on the continent's own plate) is trailing by construction; its rift
# shoulder gets this share of the full escarpment.
PASSIVE_MARGIN_ESCARPMENT = 0.5

# Land potential: km of inland distance, noise-scaled by exp(this x noise) so shelf widths
# vary, plus coastline noise in km and orogen relief at this many km per metre.
LAND_REACH_NOISE = 0.6
LAND_COAST_NOISE_KM = 250.0
LAND_OROGEN_KM_PER_M = 0.12

# Land surface. A coastal plain rises inland toward a regional upland level.
COAST_PLAIN_M = 30.0
INLAND_RISE_M = 400.0
INLAND_RISE_KM = 600.0
TEXTURE_M = 60.0
SHIELD_HILLS_M = 250.0
PLATFORM_HILLS_M = 60.0
PLATFORM_BASIN_M = 150.0
HILL_COUNTRY_SHARE = 0.35
HILL_COUNTRY_M = 400.0
HILL_COUNTRY_LIFT_M = 200.0
PLAINS_HILLS_M = 100.0
GLACIAL_LATITUDE_DEG = 50.0
GLACIAL_PIT_M = 60.0
INTERIOR_COAST_FADE_KM = 150.0
OROGEN_INTERIOR_FADE_M = 1_500.0
MIN_LAND_M = 5.0

# Submerged continental crust: shelf to ~150 m, then marginal-plateau depths further out,
# and a slope down to the oceanic contact at the crust edge.
SHELF_MIN_KM, SHELF_SPAN_KM = 80.0, 220.0
SHELF_LEADING_FACTOR = 0.25
SHELF_BREAK_M = 150.0
MARGINAL_PLATEAU_DROP_M = 1_800.0
MARGINAL_PLATEAU_KM = 400.0
CRUST_EDGE_SLOPE_KM = 100.0
CRUST_EDGE_DEPTH_M = -3_200.0
MAX_SUBMERGED_M = -5.0


@dataclass
class _Cloud:
    """The world's nodes in plate order, with per-node fields this pass reads."""

    points: np.ndarray
    elevation: np.ndarray
    continental: np.ndarray
    craton: np.ndarray
    plate_id: np.ndarray
    omega: np.ndarray  # (N, 3) each node's plate's predicted angular velocity, rad/yr


def _smoothstep(t: np.ndarray) -> np.ndarray:
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _distance_km(sources: np.ndarray, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Great-circle distance from each of `points` to the nearest of `sources`, and its index."""
    if len(sources) == 0:
        return np.full(len(points), np.inf), np.zeros(len(points), dtype=int)
    chord, index = cKDTree(sources).query(points, workers=-1)
    return 2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0)) * PLANET_RADIUS_KM, index


def _gather(world: "World", omegas: dict[int, np.ndarray]) -> _Cloud:
    parts: dict[str, list[np.ndarray]] = {k: [] for k in ("points", "elevation", "continental", "craton", "plate_id", "omega")}
    for plate in world.plates:
        points, elevation = plate.all_points_and_elevation()
        parts["points"].append(points)
        parts["elevation"].append(elevation)
        parts["continental"].append(
            effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        )
        parts["craton"].append(plate.collect("craton_crust_m") > 0.0)
        parts["plate_id"].append(np.full(len(points), plate.plate_id))
        parts["omega"].append(np.tile(omegas[plate.plate_id], (len(points), 1)))
    joined = {k: np.concatenate(v) for k, v in parts.items()}
    return _Cloud(joined["points"], joined["elevation"], joined["continental"], joined["craton"], joined["plate_id"], joined["omega"])


def _collision_contacts(world: "World", omegas: dict[int, np.ndarray]) -> np.ndarray:
    """Per node in plate order: continental crust within the regime reach of another plate's
    continental crust that the first predicted motion closes on."""
    spacing = line_spacing_rad(world.node_density)
    reach = max(REGIME_REACH_KM / PLANET_RADIUS_KM, REGIME_REACH_SPACINGS * spacing)
    masks = []
    for plate in world.plates:
        n = len(plate.cell_keys)
        neighbours = plate.get_neighbours(world.plates, threshold_rad=reach)
        if not neighbours:
            masks.append(np.zeros(n, dtype=bool))
            continue
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach)
        neighbour_omega = np.zeros_like(inputs.neighbor_omega)
        for neighbour in neighbours:
            neighbour_omega[inputs.neighbor_plate_id == neighbour.plate_id] = omegas[neighbour.plate_id]
        closing = rheology.normal_closing_rate_m_per_s(omegas[plate.plate_id], neighbour_omega, inputs.own_points, inputs.direction_to_neighbor)
        own_continental = effective_is_continental_from_codes(inputs.own_crust_type_codes, plate.crust_type == "continental")
        masks.append(
            (inputs.dist_to_neighbor <= reach)
            & (inputs.neighbor_plate_id >= 0)
            & (closing > torque.BOUNDARY_TRANSFORM_RATE_M_PER_S)
            & own_continental
            & ~np.asarray(inputs.neighbor_is_oceanic, dtype=bool)
        )
    return np.concatenate(masks)


def _third_craton_km(order: np.ndarray, distances: np.ndarray, lo_col: int, hi_col: int) -> np.ndarray:
    """Per node, the distance to its nearest craton other than columns `lo_col`/`hi_col`, from
    each node's three nearest (`order`, columns of `distances`); inf with only two cratons."""
    if order.shape[1] < 3:
        return np.full(len(distances), np.inf)
    rows = np.arange(len(distances))
    first_other = np.where((order[:, 0] != lo_col) & (order[:, 0] != hi_col), 0, np.where((order[:, 1] != lo_col) & (order[:, 1] != hi_col), 1, 2))
    return distances[rows, order[rows, first_other]]


def _neighbour_pairs(order: np.ndarray, distances: np.ndarray, ids: list[int], continent_of_core: dict[int, int]) -> list[tuple[int, int]]:
    """(lo, hi) column pairs of cratons that are some node's two nearest, on one continent and
    within SUTURE_PAIR_MAX_KM of each other there."""
    rows = np.arange(len(distances))
    near = distances[rows, order[:, 0]] + distances[rows, order[:, 1]] <= SUTURE_PAIR_MAX_KM
    pairs = {(min(a, b), max(a, b)) for a, b in zip(order[near, 0].tolist(), order[near, 1].tolist())}
    return sorted((a, b) for a, b in pairs if continent_of_core[ids[a]] == continent_of_core[ids[b]])


def suture_relief(
    points: np.ndarray, core_id: np.ndarray, continent_of_core: dict[int, int], warp: np.ndarray, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    """(relief m, young-class relief m) per point from belts between neighbouring cratons.

    `core_id` is each point's craton (-1 outside cratons); belts join only cratons on the same
    continent that are some point's two nearest. Each surviving pair's belt is evaluated on
    every point at signed distance (d_hi - d_lo) / 2 from the pair's bisector (positive on the
    lower-id craton's side) plus `SUTURE_WARP_KM * warp`; orienting by the pair lets the warp
    move a belt sideways instead of splitting it. A belt fades where the pair is far apart
    and where a third craton is nearer than both, and a point takes the highest belt, so
    belts meet without cliffs. A young belt's plateau lies on one side and fades out before
    reaching that side's craton. Each pair's class, plateau side and survival come from an rng
    keyed by the seed and the pair, so they don't depend on node order or density."""
    relief = np.zeros(len(points))
    young = np.zeros(len(points))
    ids = sorted(continent_of_core)
    if len(ids) < 2:
        return relief, young
    distances = np.stack([_distance_km(points[core_id == c], points)[0] for c in ids], axis=1)
    order = np.argsort(distances, axis=1)[:, :3]
    for lo_col, hi_col in _neighbour_pairs(order, distances, ids, continent_of_core):
        lo, hi = ids[lo_col], ids[hi_col]
        rng = np.random.default_rng((seed, lo, hi, 0x5B2))
        if rng.random() < SUTURE_DROP_PROBABILITY:
            continue
        cls = int(rng.choice(len(SUTURE_CLASSES), p=[c[0] for c in SUTURE_CLASSES]))
        _, peak, sigma, plateau = SUTURE_CLASSES[cls]
        plateau_on_lo = rng.random() < 0.5
        d_lo, d_hi = distances[:, lo_col], distances[:, hi_col]
        s = (d_hi - d_lo) / 2.0 + SUTURE_WARP_KM * warp
        r = peak * np.exp(-0.5 * (s / sigma) ** 2)
        if plateau:
            # s >= 0 is the lower-id craton's side of the warped bisector.
            on_side = (s >= 0.0) == plateau_on_lo
            reach = 1.0 - _smoothstep((np.abs(s) - SUTURE_PLATEAU_KM) / SUTURE_PLATEAU_EDGE_KM)
            before_craton = _smoothstep(np.where(plateau_on_lo, d_lo, d_hi) / CRATON_EDGE_FADE_KM)
            r = np.maximum(r, np.where(on_side, plateau * reach * before_craton, 0.0))
        pair_km = d_lo + d_hi
        fade = 1.0 - _smoothstep((pair_km - 0.75 * SUTURE_PAIR_MAX_KM) / (0.25 * SUTURE_PAIR_MAX_KM))
        third = _third_craton_km(order, distances, lo_col, hi_col) - np.minimum(d_lo, d_hi)
        r *= fade * _smoothstep(0.5 + third / (2.0 * SUTURE_THIRD_CRATON_TAPER_KM))
        relief = np.maximum(relief, r)
        if cls == SUTURE_YOUNG:
            young = np.maximum(young, r)
    return relief, young


def land_mask(potential: np.ndarray, continental: np.ndarray, n_land: int) -> np.ndarray:
    """The `n_land` continental nodes with the highest land potential (all of them if fewer).
    Oceanic crust never becomes land, and ties break by node order."""
    ranked = np.where(continental, potential, -np.inf)
    count = min(int(n_land), int(continental.sum()))
    land = np.zeros(len(potential), dtype=bool)
    land[np.argsort(-ranked, kind="stable")[:count]] = True
    return land


@dataclass
class _Noise:
    regional: IsotropicSphereNoise
    along: IsotropicSphereNoise
    ridged: IsotropicSphereNoise
    shelf_width: IsotropicSphereNoise
    coast: IsotropicSphereNoise
    reach: IsotropicSphereNoise
    texture: IsotropicSphereNoise
    suture_warp: IsotropicSphereNoise
    shield: IsotropicSphereNoise
    hill_mask: IsotropicSphereNoise
    hills: IsotropicSphereNoise
    basin: IsotropicSphereNoise
    pits: IsotropicSphereNoise

    @classmethod
    def for_seed(cls, seed: int) -> "_Noise":
        # Fixed draw order; a separate stream from every other generation rng.
        rng = np.random.default_rng((seed, 0x218))
        return cls(
            regional=IsotropicSphereNoise(rng, octaves=2, base_freq=3.0),
            along=IsotropicSphereNoise(rng, octaves=3, base_freq=8.0),
            ridged=IsotropicSphereNoise(rng, octaves=4, base_freq=40.0, gain=0.55),
            shelf_width=IsotropicSphereNoise(rng, octaves=2, base_freq=6.0),
            coast=IsotropicSphereNoise(rng, octaves=6, base_freq=10.0, gain=0.6),
            reach=IsotropicSphereNoise(rng, octaves=2, base_freq=4.0),
            texture=IsotropicSphereNoise(rng, octaves=6, base_freq=12.0),
            suture_warp=IsotropicSphereNoise(rng, octaves=2, base_freq=5.0),
            shield=IsotropicSphereNoise(rng, octaves=2, base_freq=8.0),
            hill_mask=IsotropicSphereNoise(rng, octaves=2, base_freq=7.0),
            hills=IsotropicSphereNoise(rng, octaves=3, base_freq=50.0, gain=0.6),
            basin=IsotropicSphereNoise(rng, octaves=2, base_freq=14.0),
            pits=IsotropicSphereNoise(rng, octaves=2, base_freq=160.0),
        )


def leading_edge(
    points: np.ndarray, omega: np.ndarray, plate_id: np.ndarray, edge_point: np.ndarray, edge_omega: np.ndarray, edge_plate_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """(leading score in [-1, 1], passive mask) per continental node.

    The score is the node's velocity relative to the plate that owns its nearest oceanic
    node (`edge_*`), along the outward tangent toward that node, over the relative speed and
    scaled down below `LEADING_EDGE_FULL_SPEED_CM_YR`: positive where the continent overrides
    that ocean. Where the ocean is on the node's own plate nothing converges, so the score is
    0 and the margin is passive. omega is rad/yr, so omega x point is in planet radii per year."""
    outward = edge_point - points
    outward -= np.sum(outward * points, axis=1, keepdims=True) * points
    outward /= np.maximum(np.linalg.norm(outward, axis=1, keepdims=True), 1e-12)
    relative = np.cross(omega - edge_omega, points) * PLANET_RADIUS_KM * 1e5  # cm/yr
    speed = np.linalg.norm(relative, axis=1)
    along = np.sum(relative * outward, axis=1) / np.maximum(speed, 1e-9)
    passive = edge_plate_id == plate_id
    score = along * np.clip(speed / LEADING_EDGE_FULL_SPEED_CM_YR, 0.0, 1.0)
    return np.where(passive, 0.0, score), passive


def _continent_of_core(layout: "ContinentLayout") -> dict[int, int]:
    """Each craton's continent (its most common one); cratons with no continental point are left out."""
    out = {}
    for c in np.unique(layout.core[layout.core >= 0]):
        continents = layout.continent[(layout.core == c) & (layout.continent >= 0)]
        if len(continents):
            out[int(c)] = int(np.bincount(continents).argmax())
    return out


def _orogens(cloud, noise, layout, seed, d_edge, lead, passive, collision_contact):
    """(orogen relief m, young share in [0, 1], escarpment m) on continental crust."""
    p = cloud.points
    core_id = np.asarray(layout.core[layout.nearest(p)])
    suture, suture_young = suture_relief(p, core_id, _continent_of_core(layout), np.tanh(noise.suture_warp.sample(p)), seed)
    cordillera = CORDILLERA_PEAK_M * np.clip(lead, 0.0, 1.0) * np.exp(-0.5 * ((d_edge - CORDILLERA_INLAND_KM) / CORDILLERA_SIGMA_KM) ** 2)
    d_collision, _ = _distance_km(p[collision_contact], p)
    collision = COLLISION_PEAK_M * np.exp(-0.5 * (d_collision / COLLISION_SIGMA_KM) ** 2)
    margin = np.maximum(cordillera, collision)
    d_craton, _ = _distance_km(p[cloud.craton], p)
    craton_damping = CRATON_OROGEN_FACTOR + (1.0 - CRATON_OROGEN_FACTOR) * _smoothstep(d_craton / CRATON_EDGE_FADE_KM)
    along = np.clip(0.75 + 0.35 * noise.along.sample(p), 0.3, 1.3)
    ridged = np.clip(1.0 - np.abs(noise.ridged.sample(p)) / 1.8, 0.0, 1.0) ** 1.5
    shape = along * (0.55 + 0.45 * ridged) * craton_damping
    orogen = np.where(cloud.continental, (suture + margin) * shape, 0.0)
    young = np.where(cloud.continental, (suture_young + margin) * shape, 0.0)
    young_share = np.where(orogen > 0.0, np.clip(young / np.maximum(orogen, 1e-9), 0.0, 1.0), 0.0)
    trailing = np.where(passive, PASSIVE_MARGIN_ESCARPMENT, np.clip(-lead, 0.0, 1.0))
    escarpment = ESCARPMENT_PEAK_M * trailing * np.exp(-0.5 * ((d_edge - ESCARPMENT_INLAND_KM) / ESCARPMENT_SIGMA_KM) ** 2)
    return orogen, young_share, escarpment


def _interior_relief(cloud: _Cloud, noise: _Noise, base: np.ndarray, d_coast: np.ndarray) -> np.ndarray:
    """Interior terrain types (m) on land, before fading inside mountain belts."""
    p = cloud.points
    hills = noise.hills.sample(p)
    shield = cloud.craton & (noise.shield.sample(p) > 0.0)
    platform = cloud.craton & ~shield
    hill_field = noise.hill_mask.sample(p)
    threshold = np.quantile(hill_field[cloud.continental], 1.0 - HILL_COUNTRY_SHARE) if cloud.continental.any() else np.inf
    hilly = ~cloud.craton & (hill_field > threshold)
    relief = np.where(shield, SHIELD_HILLS_M * hills, 0.0)
    relief += np.where(platform, PLATFORM_HILLS_M * hills - PLATFORM_BASIN_M * _smoothstep(-noise.basin.sample(p)), 0.0)
    relief += np.where(hilly, HILL_COUNTRY_M * hills + HILL_COUNTRY_LIFT_M, 0.0)
    relief += np.where(~cloud.craton & ~hilly, PLAINS_HILLS_M * hills, 0.0)
    latitude = np.degrees(np.arcsin(np.clip(p[:, 2], -1.0, 1.0)))
    glacial = (np.abs(latitude) > GLACIAL_LATITUDE_DEG) & (shield | (base < INLAND_RISE_M))
    relief -= np.where(glacial, GLACIAL_PIT_M * _smoothstep(noise.pits.sample(p) - 0.6), 0.0)
    return relief * np.minimum(1.0, d_coast / INTERIOR_COAST_FADE_KM)


def _target_elevation(cloud, noise, layout, seed, n_land, collision_contact):
    """(target elevation relative to sea level, land mask, young share) per node; oceanic
    nodes get NaN (left as generated)."""
    p = cloud.points
    oceanic = ~cloud.continental
    d_edge, edge_index = _distance_km(p[oceanic], p)
    if oceanic.any():
        lead, passive = leading_edge(
            p, cloud.omega, cloud.plate_id, p[oceanic][edge_index], cloud.omega[oceanic][edge_index], cloud.plate_id[oceanic][edge_index]
        )
    else:  # an all-continental world has no margins
        lead, passive = np.zeros(len(p)), np.zeros(len(p), dtype=bool)
    orogen, young_share, escarpment = _orogens(cloud, noise, layout, seed, d_edge, lead, passive, collision_contact)

    potential = d_edge * np.exp(LAND_REACH_NOISE * noise.reach.sample(p)) + LAND_COAST_NOISE_KM * noise.coast.sample(p)
    land = land_mask(potential + LAND_OROGEN_KM_PER_M * orogen, cloud.continental, n_land)
    d_coast, _ = _distance_km(p[~land], p)  # land nodes: distance to the sea
    d_land, _ = _distance_km(p[land], p)  # sea nodes: distance to land

    regional = 0.5 + 0.5 * np.tanh(noise.regional.sample(p))
    base = COAST_PLAIN_M + INLAND_RISE_M * (1.0 - np.exp(-d_coast / INLAND_RISE_KM)) * (0.6 + 0.8 * regional)
    texture = TEXTURE_M * noise.texture.sample(p) * (1.0 + base / 500.0)
    interior = _interior_relief(cloud, noise, base, d_coast) * (1.0 - _smoothstep(orogen / OROGEN_INTERIOR_FADE_M))
    z_land = np.maximum(base + texture + interior + orogen + escarpment, MIN_LAND_M)

    width = (SHELF_MIN_KM + SHELF_SPAN_KM * (0.5 + 0.5 * np.tanh(noise.shelf_width.sample(p)))) * np.where(lead > 0.3, SHELF_LEADING_FACTOR, 1.0)
    z_shelf = -(20.0 + (SHELF_BREAK_M - 20.0) * _smoothstep(d_land / width))
    z_shelf -= MARGINAL_PLATEAU_DROP_M * _smoothstep((d_land - width) / MARGINAL_PLATEAU_KM)
    t = _smoothstep(1.0 - d_edge / CRUST_EDGE_SLOPE_KM)
    z_sea = np.minimum((1.0 - t) * z_shelf + t * CRUST_EDGE_DEPTH_M, MAX_SUBMERGED_M)

    target = np.where(land, z_land, z_sea)
    return np.where(cloud.continental, target, np.nan), land, young_share


def back_columns(
    target_m: np.ndarray, hc_m: np.ndarray, hm_m: np.ndarray, rho_c: np.ndarray, core: np.ndarray, young_share: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """(Hc, Moho thermal lag C) that float each column at `target_m` (absolute elevation) on
    its own Hm. NaN targets keep their Hc. Cratonic cores stay within the craton host window
    (`cratons.CRATON_HOST_MIN_HC_M`-`CRATON_HOST_MAX_HC_M`): neither a drowned shelf nor an
    orogen, which a core's crust couldn't become later either (#218 review). Young belts take the lag of crust buried from the reference column, scaled by
    their young share; others start at steady state."""
    changed = np.isfinite(target_m)
    hc = np.where(changed, lithosphere.crustal_thickness_for_elevation(np.where(changed, target_m, 0.0), hm_m, rho_c), hc_m)
    hc = np.clip(hc, lithosphere.MIN_CRUSTAL_THICKNESS_M, lithosphere.MAX_CRUSTAL_THICKNESS_M)
    hc = np.where(changed & core, np.clip(hc, cratons.CRATON_HOST_MIN_HC_M, cratons.CRATON_HOST_MAX_HC_M), hc)
    hc0, hm0 = lithosphere.reference_thickness("continental")
    buried = orogeny.bury_moho(np.zeros(len(hc)), np.full(len(hc), hc0), np.full(len(hc), hm0), hc, hm_m)
    lag = np.where(changed, np.clip(young_share, 0.0, 1.0) * np.maximum(buried, 0.0), 0.0)
    return hc, lag


def apply_synthetic_relief(world: "World", layout: "ContinentLayout") -> None:
    """Replace continental relief on a freshly generated cratons-first world with relief from
    a synthetic tectonic history, backed by Hc (see the module docstring). Keeps the land
    node count; leaves oceanic crust, plates and elapsed time untouched."""
    omegas = {p.plate_id: torque.predict_plate_motion(p, world, world.plates, FIRST_MOTION_STEP_YEARS) for p in world.plates}
    cloud = _gather(world, omegas)
    if not cloud.continental.any():
        return
    n_land = int((cloud.elevation > world.sea_level_m).sum())
    target, _, young_share = _target_elevation(
        cloud, _Noise.for_seed(world.seed), layout, world.seed, n_land, _collision_contacts(world, omegas)
    )
    offset = 0
    for plate in world.plates:
        n = len(plate.cell_keys)
        sl = slice(offset, offset + n)
        offset += n
        craton_m = plate.collect("craton_crust_m")
        core = craton_m > 0.0
        hc, lag = back_columns(
            target[sl] + world.sea_level_m,
            plate.collect("crustal_thickness_m"),
            plate.collect("mantle_lithosphere_thickness_m"),
            lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type),
            core,
            young_share[sl],
        )
        # A generated core is cratonic through its whole crust (cratons.seed_cores).
        plate.set_fields_on_plate(crustal_thickness_m=hc, moho_thermal_lag_c=lag, craton_crust_m=np.where(core, hc, craton_m))
        lithosphere.sync_plate_elevation(plate)

"""Cratons-first world layout: cores, then continents, then plates (RFC 0009, #462).

Random worlds used to be built plates first: a Voronoi tiling, each plate wholly continental
or oceanic, and a craton in the middle of every continental plate. Every continental edge was
then a plate boundary, and the tiling alone set each craton's place, size and shape. On
Earth the order runs the other way: cratons are the oldest parts of the continents, the
continents were assembled around them along orogenic belts, and today's plates are younger
and cut across that assembly. This module lays a world out in that order, before any plate
exists, on a whole-sphere lattice at the world's node spacing (`layout_lattice.Lattice`):

1. **Continents.** `CONTINENT_COUNT_RANGE` centres, Poisson-disc sampled
   `CONTINENT_SPACING_KM` apart, split the continental area (`continental_fraction` of the
   sphere) between them by a Dirichlet draw (`CONTINENT_SIZE_CONCENTRATION`).
2. **Cratons.** Each continent's core nuclei chain out from its centre. Each new nucleus
   starts a new craton, at least `CRATON_CORE_SPACING_KM` from every other craton's nuclei,
   with probability `CRATON_CLUSTER_CONCENTRATION / (that + n)`; otherwise it joins an
   existing craton, chosen in proportion to its nuclei, overlapping one of its nuclei (a
   Chinese restaurant process, #470). A craton is the points within a noise-perturbed
   `CRATON_CORE_RADIUS_KM` of its nuclei, overlapping claims going to the nearest nucleus,
   less a band either side of the line halfway between two cratons, so cratons are at least
   `CRATON_BELT_WIDTH_KM` apart. Nuclei are added until the cratons can cover
   `CRATON_SHARE_OF_CONTINENT` of the continent; the radius is then trimmed to hit that
   share. A continent smaller than one core disk gets none.
3. **Belts and margins.** Each continent grows outward from its cratons (from its centre if
   it has none) until it reaches its area, in order of growth time: the distance to the
   cratons plus `CONTINENT_CENTRE_PULL` times the distance to the centre, over a smoothly
   varying growth speed, plus outline noise. So cratons on the far side of a continent, or
   where it grows slowly, come near its edge. The closure of its cores
   (dilated, then eroded, by `BELT_CLOSURE_KM`) is the core-and-belt interior: everything in
   it that isn't core is belt, so the whole gap between neighbouring cores is belt and no
   unmarked gap is left to form craton later. Continental points outside it are margin.
4. **Plates** are laid over the result by `plate_partition`, on the same
   `layout_lattice.Lattice`. `ContinentLayout` holds both stages' per-point output.

`lithosphere_plate.generate_plates` builds mixed-crust plates from the result, behind
`CRATONS_FIRST_DEFAULT` until RFC 0009's gate G3 decides whether it becomes the default.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass, fields

import numpy as np
from scipy.spatial import cKDTree

from .elevation_lines import PLANET_RADIUS_KM
from .layout_lattice import Lattice, arc_km, cap_area_km2, lattice, poisson_disc
from .noise import SphereNoise
from .plate_partition import PlatePartition

# Off until RFC 0009's gate G3 picks the layout: main's random worlds stay plates-first.
CRATONS_FIRST_DEFAULT = False

# --- Continents --------------------------------------------------------------------------

# Design choices (RFC 0009 knobs), tuned to the layout measures over 20 seeds.
CONTINENT_COUNT_RANGE = (4, 7)
CONTINENT_SPACING_KM = 4000.0
# Dirichlet concentration: 2 gives clearly unequal continents without slivers.
CONTINENT_SIZE_CONCENTRATION = 2.0
# Continental area share of the sphere is clamped here so there is always some of each.
CONTINENTAL_FRACTION_BOUNDS = (0.05, 0.9)
# Amplitude of the noise on the distance a continent grows outward from its cores, so
# coastlines and margins aren't parallel to the core outlines.
CONTINENT_OUTLINE_NOISE_KM = 400.0
# A continent grows outward from its cores at a speed that varies smoothly with position,
# exp(CONTINENT_GROWTH_SPEED_SIGMA * noise), so some sides get wide margins and others leave
# a craton at the coast, as West Africa's, Australia's and Baltica's are (#470). The noise's
# lowest frequency is CONTINENT_GROWTH_NOISE_BASE_FREQ cycles per radian: a continent-scale
# pattern, ~13,000 km down to ~3,300 km. Design choice; tuned to the share of cratons at an
# edge.
CONTINENT_GROWTH_SPEED_SIGMA = 1.0
# Weight of a point's distance to its continent's centre in the growth order, beside its
# distance to the cores: the outer side of a craton far from the centre is reached late, so
# its margin is thin. With the speed field, a median 35% (11-49%) of cratons come within
# CRATON_EDGE_NEAR_KM of the continental edge over 40 seeds, against 0% with neither (#470).
CONTINENT_CENTRE_PULL = 1.0
CONTINENT_GROWTH_NOISE_BASE_FREQ = 3.0

# --- Cores and belts ---------------------------------------------------------------------

# Fixed to today's mean starting craton share (37.6-50.7% of genuine continental area on
# seeds 7, 11 and 2024, #403) so craton calibrations carry over. Capped so every continent
# keeps some non-cratonic crust.
CRATON_SHARE_OF_CONTINENT = 0.4
CRATON_SHARE_MAX = 0.7
# Minimum nucleus separation between two full-size cratons, scaled by their mean size scale
# for smaller ones -- not a count of geological cratons.
CRATON_CORE_SPACING_KM = 1500.0
# About 3.1 M km^2 per full disk.
CRATON_CORE_RADIUS_KM = 1000.0
# Orogens welded onto shields are a few hundred km wide (Cawood et al., 2009). 400 km, in
# the RFC's 200-600 km sweep, rather than its starting 300: with composite cratons (#470)
# there are fewer, larger cores and so less gap between them, and at 300-350 km belts fell
# to 9-10% of continental area on some seeds, under the RFC's 10-30%; at 400 km, 11-20%.
CRATON_BELT_WIDTH_KM = 400.0
# Noise amplitude on a point's distance to its nucleus, for irregular core outlines.
CRATON_CORE_OUTLINE_NOISE_KM = 250.0
# A nucleus that starts a new craton is placed this far (in multiples of the core spacing,
# scaled by the two cratons' mean size scale) from an existing nucleus of the same continent,
# so a continent's cratons stay a connected assembly.
NUCLEUS_STEP_RANGE = (1.0, 1.25)
NUCLEUS_MAX_TRIES = 1000
# Cratons are composites (#470): within a continent each new nucleus starts a new craton
# with probability alpha / (alpha + n), n the continent's nuclei so far, and otherwise joins
# an existing craton chosen in proportion to its nuclei (a Chinese restaurant process).
# Design choice; tuned to the craton count and size spread over the layout seeds.
# At 8, with the size scales below, a median 2 cratons per world reach 10 M km^2 (0 at 20,
# #470's first value) and 20-55 cratons form over 40 seeds; at 5-6, belts fell to 5-7% of
# continental area on some seeds. At 2 with every scale 1, 89% of nuclei joined.
CRATON_CLUSTER_CONCENTRATION = 8.0
# A joining nucleus is placed this far (in multiples of its craton's reach, the size scale
# times the core radius) from a nucleus of its craton, so the two disks overlap and merge into
# one irregular core.
NUCLEUS_JOIN_STEP_RANGE = (0.8, 1.3)
# Each craton draws a size scale from this range when it starts, and reaches its scale times
# CRATON_CORE_RADIUS_KM from its nuclei, so a craton that no nucleus joins can be much smaller
# than one disk. The top stays at 1: the largest cratons are composites of several nuclei.
# 0.5-1 puts 30% of cratons under 1 M km^2 (5th percentile 0.51 M km^2, against 1.33 with
# every scale 1) over 40 seeds (#470 review); at 0.3 with spacing scaled to size, half were
# under 1 M km^2 and the smallest 0.09.
CRATON_SIZE_SCALE_RANGE = (0.5, 1.0)
# No nucleus joins a craton that already covers this much, so composites stop growing at
# about the largest of #470's first layouts (20-23 M km^2 over 40 seeds); the trim to the
# craton share can leave one a little larger (22 M km^2).
CRATON_MAX_AREA_KM2 = 18e6
# Closing radius for the core-and-belt interior: bridges gaps between cores up to twice
# this wide, more than the widest gap the nucleus step and a trimmed core radius leave. At
# 600 km belts were 9.7-16% of continental area over 20 seeds, under the RFC's 10-30%; at
# 800 km, 10.8-18% (#462).
BELT_CLOSURE_KM = 800.0
# A point whose distances to its two nearest cratons sum to at most this many belt widths lies between them
# and is belt: the gap a closing can't bridge between two small cratons far apart (#470
# review). Twice the belt width: with size-aware spacing, neighbouring cratons' gap is about
# a belt wide.
BELT_LENS_BELT_WIDTHS = 2.0
# The sinusoid noise's lowest frequency (cycles per radian): wavelengths of ~5,000 km down
# to ~1,250 km over its three octaves, the scale of a core's or continent's outline.
OUTLINE_NOISE_BASE_FREQ = 8.0

# A craton is "at an edge" when it comes within this distance of the continental-crust
# edge (diagnostics only, #470).
CRATON_EDGE_NEAR_KM = 300.0
# Continent centres are drawn by dart throwing; this bounds the draws.
CONTINENT_CENTRE_MAX_TRIES = 5000


@dataclass
class ContinentGeography:
    """Per lattice point: `continent` (-1 for ocean), `core` (craton id, -1 outside cratons),
    and `belt` and `margin` masks. `nuclei` are the core nuclei: `nucleus_craton` names each
    one's craton (several overlapping nuclei can make one), `nucleus_continent` its
    continent. `continent_area_km2` is each continent's target area."""

    lattice: Lattice
    continent: np.ndarray
    core: np.ndarray
    belt: np.ndarray
    margin: np.ndarray
    nuclei: np.ndarray
    nucleus_craton: np.ndarray
    nucleus_continent: np.ndarray
    continent_area_km2: np.ndarray

    @property
    def is_continental(self) -> np.ndarray:
        return self.continent >= 0

    @property
    def is_core(self) -> np.ndarray:
        return self.core >= 0

    def nearest(self, points_xyz: np.ndarray) -> np.ndarray:
        """Index of the lattice point nearest each of `points_xyz`."""
        return self.lattice.tree.query(points_xyz, workers=-1)[1]


@dataclass
class ContinentLayout(ContinentGeography):
    """A `ContinentGeography` with the plates laid over it (`plate_partition.PlatePartition`'s
    fields) and `land`, the generator's land mask per lattice point."""

    plate: np.ndarray
    initial_plate: np.ndarray
    sites: np.ndarray
    refine_costs: list[float]
    partition_attempts: int
    land: np.ndarray

    @classmethod
    def from_parts(cls, geography: ContinentGeography, partition: PlatePartition, land: np.ndarray) -> "ContinentLayout":
        own = {f.name: getattr(geography, f.name) for f in fields(ContinentGeography)}
        return cls(
            **own,
            plate=partition.plate,
            initial_plate=partition.initial_plate,
            sites=partition.sites,
            refine_costs=partition.refine_costs,
            partition_attempts=partition.attempts,
            land=land,
        )

    def craton_diagnostics(self, edge_km: float = CRATON_EDGE_NEAR_KM) -> dict:
        """Craton count, the share of nuclei that joined an existing craton, the craton-area
        coefficient of variation, and the share of cratons with some point within `edge_km`
        of the continental-crust edge (#470)."""
        area = self.lattice.area_km2
        ids = np.unique(self.core[self.is_core])
        sizes = np.array([area[self.core == c].sum() for c in ids])
        ocean = ~self.is_continental
        near = 0
        if len(ids) and np.any(ocean):
            gap = arc_km(cKDTree(self.lattice.xyz[ocean]).query(self.lattice.xyz[self.is_core], workers=-1)[0])
            nearest_gap = np.full(int(ids.max()) + 1, np.inf)
            np.minimum.at(nearest_gap, self.core[self.is_core], gap)
            near = int(np.sum(nearest_gap[ids] <= edge_km))
        return {
            "cratons": int(len(ids)),
            "joined_nucleus_share": 1.0 - len(np.unique(self.nucleus_craton)) / max(len(self.nucleus_craton), 1),
            "craton_area_cv": float(sizes.std() / sizes.mean()) if len(sizes) else 0.0,
            "cratons_near_edge_share": near / max(len(ids), 1),
        }

    @property
    def num_plates(self) -> int:
        return len(self.sites)

    def plate_mask(self, plate_id: int, points_xyz: np.ndarray) -> np.ndarray:
        """Whether each of `points_xyz` falls in plate `plate_id`: its nearest lattice point
        is labelled with it. Points outside the plate's bounding cap around its site (padded
        by a cell) are rejected before the nearest-point query."""
        owned = np.zeros(len(points_xyz), dtype=bool)
        near = np.flatnonzero(points_xyz @ self.sites[plate_id] >= self._cap_cos[plate_id])
        if len(near):
            owned[near] = self.plate[self.nearest(points_xyz[near])] == plate_id
        return owned

    @functools.cached_property
    def _cap_cos(self) -> np.ndarray:
        pad = 2.0 * float(np.sqrt(self.lattice.area_km2.max())) / PLANET_RADIUS_KM
        angle = np.arccos(np.clip(np.einsum("ij,ij->i", self.lattice.xyz, self.sites[self.plate]), -1.0, 1.0))
        reach = np.zeros(self.num_plates)
        np.maximum.at(reach, self.plate, angle)
        return np.cos(np.minimum(reach + pad, np.pi))

    def diagnostics(self) -> dict:
        area = self.lattice.area_km2
        continental = float(area[self.is_continental].sum())
        per_continent = []
        disk = cap_area_km2(CRATON_CORE_RADIUS_KM)
        for k, target in enumerate(self.continent_area_km2):
            in_k = self.continent == k
            per_continent.append(
                {
                    "area_km2": float(area[in_k].sum()),
                    "target_km2": float(target),
                    "cores": int(len(np.unique(self.core[in_k & self.is_core]))),
                    "room_for_two_cores": bool(CRATON_SHARE_OF_CONTINENT * target > disk),
                }
            )
        return {
            **self.craton_diagnostics(),
            "continental_area_fraction": continental / float(area.sum()),
            "core_share_of_continent": float(area[self.is_core].sum()) / max(continental, 1.0),
            "belt_share_of_continent": float(area[self.belt].sum()) / max(continental, 1.0),
            "margin_share_of_continent": float(area[self.margin].sum()) / max(continental, 1.0),
            "continents": per_continent,
            "refine_sweeps": len(self.refine_costs) - 1,
            "partition_attempts": self.partition_attempts,
        }


# --- Sampling helpers --------------------------------------------------------------------


def _step_from(point: np.ndarray, distance_rad: float, bearing: float) -> np.ndarray:
    """The point `distance_rad` from `point` along `bearing` (any fixed tangent basis)."""
    helper = np.array([0.0, 0.0, 1.0]) if abs(point[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(point, helper)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(point, e1)
    direction = np.cos(bearing) * e1 + np.sin(bearing) * e2
    return np.cos(distance_rad) * point + np.sin(distance_rad) * direction


def _fill_to_area(order_values: np.ndarray, area: np.ndarray, target: float) -> np.ndarray:
    """Boolean mask of the points with the smallest `order_values` whose area first reaches
    `target` (all of them when they can't)."""
    order = np.argsort(order_values, kind="stable")
    cumulative = np.cumsum(area[order])
    count = int(np.searchsorted(cumulative, target)) + 1
    mask = np.zeros(len(order_values), dtype=bool)
    mask[order[: min(count, len(order))]] = True
    return mask


# --- 1-3: continents, cores, belts -------------------------------------------------------


class _NucleusDistances:
    """Per lattice point: the distance in km to the nearest nucleus (`d1`), which craton that
    nucleus belongs to (`craton`), and the distance to the nearest nucleus of any other craton
    (`d2`), updated one nucleus at a time."""

    def __init__(self, lat: Lattice, core_noise_km: np.ndarray) -> None:
        self.lat = lat
        self.noise = core_noise_km
        self.d1 = np.full(lat.size, np.inf)
        self.d2 = np.full(lat.size, np.inf)
        self.craton = np.zeros(lat.size, dtype=np.int64)
        # Each craton's size scale (CRATON_SIZE_SCALE_RANGE), by craton id.
        self.scale: list[float] = []

    def add(self, nucleus: np.ndarray, craton: int) -> None:
        d = np.arccos(np.clip(self.lat.xyz @ nucleus, -1.0, 1.0)) * PLANET_RADIUS_KM
        same = self.craton == craton
        # The same craton: only its nearest distance can fall. Another craton: if it is now
        # nearest, the old nearest becomes the nearest other craton.
        closer = ~same & (d < self.d1)
        self.d2 = np.where(closer, self.d1, np.where(same, self.d2, np.minimum(self.d2, d)))
        self.d1 = np.where(same | closer, np.minimum(self.d1, d), self.d1)
        self.craton = np.where(closer, craton, self.craton)

    def eligible(self, craton_continent: np.ndarray, k: int) -> np.ndarray:
        """Points continent `k`'s cratons can claim: nearest to one of its nuclei, within the
        perturbed core radius, and clear of the belt band. Half the difference of the
        distances to the two nearest cratons is at most the distance to their bisector, so
        cratons stay at least the belt width apart."""
        clears_belt = (self.d2 - self.d1) / 2.0 >= CRATON_BELT_WIDTH_KM / 2.0
        return (craton_continent[self.craton] == k) & clears_belt & (self.perturbed <= CRATON_CORE_RADIUS_KM)

    def craton_areas_km2(self, craton_continent: np.ndarray, k: int) -> np.ndarray:
        """Area each craton could claim now (its eligible points), by craton id."""
        eligible = self.eligible(craton_continent, k)
        return np.bincount(self.craton[eligible], weights=self.lat.area_km2[eligible], minlength=len(craton_continent))

    @property
    def perturbed(self) -> np.ndarray:
        """Distance to the nearest nucleus over its craton's size scale, plus the outline
        noise: a craton claims points where this is within CRATON_CORE_RADIUS_KM."""
        scale = np.asarray(self.scale, dtype=float)[self.craton] if self.scale else 1.0
        return self.d1 / scale + self.noise


def _place_nuclei(rng: np.random.Generator, field: _NucleusDistances, centres: np.ndarray, core_targets: np.ndarray):
    """Nuclei for every continent with a core target, chained out from its centre until its
    cratons could cover the target. Each new nucleus starts a new craton or joins an existing
    one (`CRATON_CLUSTER_CONCENTRATION`); a target larger than one core disk always gets at
    least two cratons. Returns (nuclei, each nucleus's craton, each craton's continent)."""
    nuclei: list[np.ndarray] = []
    nucleus_craton: list[int] = []
    craton_continent: list[int] = []
    spacing_rad = CRATON_CORE_SPACING_KM / PLANET_RADIUS_KM
    disk = cap_area_km2(CRATON_CORE_RADIUS_KM)

    def add_and_check_cover(candidate: np.ndarray, craton: int, k: int, scale: float) -> bool:
        """Add `candidate` as a nucleus of `craton` (a new one, of size scale `scale`, when it
        equals the craton count) in continent `k` -- always, whatever this returns -- and say
        whether `k`'s cratons can now cover its target, i.e. whether to stop adding."""
        if craton == len(craton_continent):
            craton_continent.append(k)
            field.scale.append(scale)
        nuclei.append(candidate)
        nucleus_craton.append(craton)
        field.add(candidate, craton)
        if core_targets[k] > disk and craton_continent.count(k) < 2:
            return False
        return float(field.lat.area_km2[field.eligible(np.asarray(craton_continent), k)].sum()) >= core_targets[k]

    def far_from_other_cratons(candidate: np.ndarray, craton: int, scale: float) -> bool:
        """Whether `candidate`, of a craton with size scale `scale`, is at least the craton
        spacing times the two cratons' mean scale from every other craton's nuclei: two
        full-size cratons stay CRATON_CORE_SPACING_KM apart, and small ones pack closer, so
        the gap between neighbours stays about a belt wide."""
        others = [i for i, c in enumerate(nucleus_craton) if c != craton]
        if not others:
            return True
        angles = np.arccos(np.clip(np.asarray([nuclei[i] for i in others]) @ candidate, -1.0, 1.0))
        need = spacing_rad * 0.5 * (scale + np.asarray([field.scale[nucleus_craton[i]] for i in others]))
        return bool(np.all(angles >= need))

    for k, target in enumerate(core_targets):
        if target <= 0.0:
            continue
        if add_and_check_cover(centres[k], len(craton_continent), k, float(rng.uniform(*CRATON_SIZE_SCALE_RANGE))):
            continue
        for _ in range(NUCLEUS_MAX_TRIES):
            mine = [i for i, c in enumerate(nucleus_craton) if craton_continent[c] == k]
            cratons_k = sorted({nucleus_craton[i] for i in mine})
            sizes = np.array([sum(1 for i in mine if nucleus_craton[i] == c) for c in cratons_k], dtype=float)
            sizes[field.craton_areas_km2(np.asarray(craton_continent), k)[cratons_k] >= CRATON_MAX_AREA_KM2] = 0.0
            # A new craton when the draw says so, when the target needs a second craton, or
            # when every craton is at the size cap.
            force_new = (target > disk and len(cratons_k) < 2) or not np.any(sizes > 0.0)
            if force_new or rng.random() < CRATON_CLUSTER_CONCENTRATION / (CRATON_CLUSTER_CONCENTRATION + len(mine)):
                craton = len(craton_continent)
                scale = float(rng.uniform(*CRATON_SIZE_SCALE_RANGE))
                # From the assembly's edge outward: a nucleus deep inside a large composite
                # has no room for a new craton beside it.
                reach = np.arccos(np.clip(np.asarray([nuclei[i] for i in mine]) @ centres[k], -1.0, 1.0)) + 1e-3
                pick = mine[int(rng.choice(len(mine), p=reach**2 / np.sum(reach**2)))]
                parent = nuclei[pick]
                mean_scale = 0.5 * (scale + field.scale[nucleus_craton[pick]])
                step = rng.uniform(*NUCLEUS_STEP_RANGE) * spacing_rad * mean_scale
            else:
                craton = cratons_k[int(rng.choice(len(cratons_k), p=sizes / sizes.sum()))]
                scale = field.scale[craton]
                members = [i for i in mine if nucleus_craton[i] == craton]
                parent = nuclei[members[int(rng.integers(len(members)))]]
                # Within the craton's own reach, so its disks overlap however small it is.
                step = rng.uniform(*NUCLEUS_JOIN_STEP_RANGE) * field.scale[craton] * CRATON_CORE_RADIUS_KM / PLANET_RADIUS_KM
            candidate = _step_from(parent, step, rng.uniform(0.0, 2.0 * np.pi))
            if not far_from_other_cratons(candidate, craton, scale):
                continue
            if int(np.argmax(centres @ candidate)) != k:
                continue
            if add_and_check_cover(candidate, craton, k, scale):
                break
        _inflate_to_cover(field, np.asarray(craton_continent), k, target)
    return np.asarray(nuclei).reshape(-1, 3), np.asarray(nucleus_craton, dtype=np.int64), np.asarray(craton_continent, dtype=np.int64)


def _inflate_to_cover(field: _NucleusDistances, craton_continent: np.ndarray, k: int, target: float) -> None:
    """When continent `k`'s nuclei ran out of room before its cratons could cover `target`
    (a crowded continent of small cratons), grow its cratons' size scales together, each
    capped at 1, until they can."""
    mine = np.flatnonzero(craton_continent == k)
    while float(field.lat.area_km2[field.eligible(craton_continent, k)].sum()) < target:
        scales = np.asarray(field.scale)
        if np.all(scales[mine] >= 1.0):
            return
        scales[mine] = np.minimum(scales[mine] * 1.1, 1.0)
        field.scale = scales.tolist()


def _cores(field: _NucleusDistances, craton_continent: np.ndarray, core_targets: np.ndarray) -> np.ndarray:
    """Craton id per lattice point, -1 outside every craton: each continent's eligible
    points nearest their nuclei, up to its core target."""
    core = np.full(field.lat.size, -1, dtype=np.int64)
    perturbed = field.perturbed
    for k in np.unique(craton_continent):
        eligible = np.flatnonzero(field.eligible(craton_continent, k))
        chosen = eligible[_fill_to_area(perturbed[eligible], field.lat.area_km2[eligible], core_targets[k])]
        core[chosen] = field.craton[chosen]
    return core


def _continents(
    lat: Lattice,
    centres: np.ndarray,
    areas: np.ndarray,
    core: np.ndarray,
    belt: np.ndarray,
    owners: np.ndarray,
    outline_noise_km: np.ndarray,
    growth_speed: np.ndarray,
):
    """Continent id per point, -1 for ocean. A continent takes the points with the
    smallest growth time, distance over `growth_speed`, plus the outline noise, after its
    cratons and belts (`belt`, continent id per point). `owners` is each craton's continent."""
    n_c = len(centres)
    distance = np.empty((lat.size, n_c))
    for k in range(n_c):
        in_core = (core >= 0) & (owners[np.maximum(core, 0)] == k)
        if np.any(in_core):
            distance[:, k] = arc_km(cKDTree(lat.xyz[in_core]).query(lat.xyz, workers=-1)[0])
        else:
            distance[:, k] = np.arccos(np.clip(lat.xyz @ centres[k], -1.0, 1.0)) * PLANET_RADIUS_KM
    to_centre = np.arccos(np.clip(lat.xyz @ centres.T, -1.0, 1.0)) * PLANET_RADIUS_KM
    perturbed = (distance + CONTINENT_CENTRE_PULL * to_centre) / growth_speed[:, None] + outline_noise_km[:, None]
    owner = np.argmin(perturbed, axis=1)
    continent = np.full(lat.size, -1, dtype=np.int64)
    for k in range(n_c):
        in_core = ((core >= 0) & (owners[np.maximum(core, 0)] == k)) | (belt == k)
        continent[in_core] = k
        rest = np.flatnonzero((owner == k) & ~in_core & (continent < 0))
        remaining = areas[k] - float(lat.area_km2[in_core].sum())
        if remaining > 0.0 and len(rest):
            continent[rest[_fill_to_area(perturbed[rest, k], lat.area_km2[rest], remaining)]] = k
    # A continent whose own region was too small for its area (the centre pull can shrink
    # it) takes the shortfall from points no continent claimed, in its own growth order.
    for k in range(n_c):
        shortfall = areas[k] - float(lat.area_km2[continent == k].sum())
        free = np.flatnonzero(continent < 0)
        if shortfall > lat.area_km2.max() and len(free):
            continent[free[_fill_to_area(perturbed[free, k], lat.area_km2[free], shortfall)]] = k
    return continent


def _belts(lat: Lattice, core: np.ndarray, craton_continent: np.ndarray) -> np.ndarray:
    """Continent id of each belt point (-1 elsewhere): the non-craton points in the closure of
    a continent's cratons, or in the lens between two of them. Computed before the continent
    grows, which then always contains them, so no sea opens between two of its cratons (#470
    review). The closing runs over the whole sphere; the lens -- points whose distances to the
    two nearest cratons sum to at most BELT_LENS_BELT_WIDTHS belt widths -- covers the gap between two small
    cratons far apart, which a closing can't bridge."""
    belt = np.full(lat.size, -1, dtype=np.int64)
    has_core = core >= 0
    for k in np.unique(craton_continent):
        mine = has_core & (craton_continent[np.maximum(core, 0)] == k)
        if not np.any(mine):
            continue
        dilated = arc_km(cKDTree(lat.xyz[mine]).query(lat.xyz, workers=-1)[0]) <= BELT_CLOSURE_KM
        candidates = np.flatnonzero(dilated & ~has_core & (belt < 0))
        if len(candidates) == 0:
            continue
        inside = _lens(lat, core, np.unique(core[mine]), candidates)
        if not np.all(dilated):
            inside |= arc_km(cKDTree(lat.xyz[~dilated]).query(lat.xyz[candidates], workers=-1)[0]) > BELT_CLOSURE_KM
        belt[candidates[inside]] = k
    return belt


def _lens(lat: Lattice, core: np.ndarray, cratons: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Whether each of `points` lies between two of `cratons`: its distances to the nearest
    two sum to at most BELT_LENS_BELT_WIDTHS belt widths."""
    if len(cratons) < 2:
        return np.zeros(len(points), dtype=bool)
    lens_km = BELT_LENS_BELT_WIDTHS * CRATON_BELT_WIDTH_KM
    bound = 2.0 * np.sin(lens_km / PLANET_RADIUS_KM / 2.0)
    d1 = np.full(len(points), np.inf)
    d2 = np.full(len(points), np.inf)
    for craton in cratons:
        chord = cKDTree(lat.xyz[core == craton]).query(lat.xyz[points], distance_upper_bound=bound, workers=-1)[0]
        d = np.where(np.isfinite(chord), arc_km(np.where(np.isfinite(chord), chord, 0.0)), np.inf)
        d2 = np.where(d < d1, d1, np.minimum(d2, d))
        d1 = np.minimum(d1, d)
    return d1 + d2 <= lens_km


def build_geography(rng: np.random.Generator, spacing_rad: float, continental_fraction: float) -> ContinentGeography:
    """Continents, cores, belts and margins (see the module docstring), deterministic in `rng`."""
    lat = lattice(spacing_rad)
    continental_fraction = float(np.clip(continental_fraction, *CONTINENTAL_FRACTION_BOUNDS))
    share = min(CRATON_SHARE_OF_CONTINENT, CRATON_SHARE_MAX)

    core_noise_km = CRATON_CORE_OUTLINE_NOISE_KM * SphereNoise(rng, base_freq=OUTLINE_NOISE_BASE_FREQ).sample(lat.xyz)
    outline_noise_km = CONTINENT_OUTLINE_NOISE_KM * SphereNoise(rng, base_freq=OUTLINE_NOISE_BASE_FREQ).sample(lat.xyz)
    growth_speed = np.exp(CONTINENT_GROWTH_SPEED_SIGMA * SphereNoise(rng, base_freq=CONTINENT_GROWTH_NOISE_BASE_FREQ).sample(lat.xyz))

    n_c = int(rng.integers(CONTINENT_COUNT_RANGE[0], CONTINENT_COUNT_RANGE[1] + 1))
    centres = poisson_disc(rng, n_c, CONTINENT_SPACING_KM / PLANET_RADIUS_KM, lat.xyz, lat.area_km2, CONTINENT_CENTRE_MAX_TRIES)
    total_area = continental_fraction * float(lat.area_km2.sum())
    areas = total_area * rng.dirichlet(np.full(len(centres), CONTINENT_SIZE_CONCENTRATION))
    disk = cap_area_km2(CRATON_CORE_RADIUS_KM)
    core_targets = np.where(areas >= disk, share * areas, 0.0)

    field = _NucleusDistances(lat, core_noise_km)
    nuclei, nucleus_craton, craton_continent = _place_nuclei(rng, field, centres, core_targets)
    core = _cores(field, craton_continent, core_targets)
    belt_of = _belts(lat, core, craton_continent)
    continent = _continents(lat, centres, areas, core, belt_of, craton_continent, outline_noise_km, growth_speed)
    belt = (belt_of >= 0) & (belt_of == continent)
    return ContinentGeography(
        lattice=lat,
        continent=continent,
        core=core,
        belt=belt,
        margin=(continent >= 0) & (core < 0) & ~belt,
        nuclei=nuclei,
        nucleus_craton=nucleus_craton,
        nucleus_continent=craton_continent[nucleus_craton] if len(nucleus_craton) else nucleus_craton,
        continent_area_km2=areas,
    )

"""Composite relief fields for seeding believable initial terrain.

`noise.SphereNoise` on its own is a short sum of random sinusoids -- continuous and cheap,
but band-limited and blobby: it gives a plate interior gentle rolling undulation and nothing
else. A freshly generated world built on it alone has no orogenic belts, no valleys between
ridges, and no plateaus -- those only ever appear later, once tectonics has run at a plate
boundary for tens of millions of years.

This module composes several `SphereNoise` fields into a richer relief signal so *newly
seeded crust already looks like a real landmass*. Two generators:

- `ContinentalRelief` -- the field continental generation seeds `Hc` from. It deliberately
  splits into two parts (see the class docstring): a low-frequency `sample()` that behaves
  just like the old single `SphereNoise` (and is the only thing that decides land vs sea,
  so `plates._land_noise_threshold`'s coarse whole-sphere quantile stays a good estimator),
  plus a non-negative `uplift()` carrying the orogenic belts and plateaus, added only on
  land so it can never move a coastline.
- `OceanicRelief` -- abyssal-hill texture plus sparse seamount chains, kept gentle so deep
  ocean interiors stay deep.

Both duck-type `SphereNoise` (`.sample(xyz) -> array`), are fully deterministic in the
`np.random.Generator` handed to `__init__`, and express their output in the same
"1 unit ~= 2000 m of elevation through Airy isostasy" convention
`lithosphere_plate._HC_NOISE_AMPLITUDE_*` already uses, so a caller keeps one Hc-noise
amplitude and multiplies.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np

from .noise import SphereNoise


class ReliefField(Protocol):
    """The structural type `plates._land_noise_threshold` and the generation `hc_at`
    closures depend on -- `SphereNoise`, `ContinentalRelief` and `OceanicRelief` all
    satisfy it."""

    def sample(self, xyz: np.ndarray) -> np.ndarray: ...


# --- tuning constants -------------------------------------------------------------------
# Every octave/frequency knob is a module constant so the whole feel of generated terrain
# can be tuned in one place (and monkeypatched in a tuning script) without touching logic.

# `sample()` is renormalized so its standard deviation over the sphere matches this -- the
# measured std of the old `SphereNoise(octaves=4, base_freq=2.5)` continental field, so land
# extent and coastline character are unchanged from before this module existed.
_TARGET_SAMPLE_STD = 0.30

# Domain-warp displacement, as a fraction of the unit radius. Small on purpose: enough to
# bend features off the sinusoidal grid SphereNoise is built on, not so much that the warp
# dominates the low-frequency signal it is warping.
_WARP_STRENGTH = 0.10
_WARP_OCTAVES, _WARP_FREQ = 2, 0.8

# `sample()` = fractal base (domain-warped) + a very-large-scale regional term (so different
# continents differ in overall character). Base frequency stays near the old field's 2.5 so
# the undulation is not domier than before; the extra octaves add finer texture on top.
_BASE_OCTAVES, _BASE_FREQ = 6, 2.6
_REGIONAL_OCTAVES, _REGIONAL_FREQ = 2, 0.8
_BASE_WEIGHT = 0.80
_REGIONAL_WEIGHT = 0.38

# Orogenic belts. `_BELT_COVERAGE` is the fraction of the sphere the belt mask covers at
# full strength before its smoothstep edge; `_BELT_EDGE` is that edge width as a fraction of
# the mask field's sampled range. A moderate frequency gives a handful of distinct belts
# rather than one continent-spanning smear.
_BELT_OCTAVES, _BELT_FREQ = 3, 1.6
_BELT_COVERAGE = 0.40
_BELT_EDGE = 0.18
# A floor under the belt mask: even outside a named belt, crust gets this fraction of the
# ridged relief, so every continent has some hill country / low ranges (old eroded orogens,
# hotspot swells) rather than a coin-flip between "alps everywhere" and "billiard table".
_BELT_FLOOR = 0.16
# Ridged noise inside a belt: `r = (1 - |fbm|) ** _RIDGE_SHARPEN` in [0, 1] (lower sharpen ->
# broader massifs, not needles), scaled by `_RIDGE_GAIN`. A belt interior between ridges
# gets ~0 uplift, so it sits at the local `sample()` baseline while the crests tower over
# it -- that gap is the intermontane valley. `_RIDGE_FREQ` is high on purpose: it sets the
# ridge-to-ridge spacing (~a few hundred km), so a belt reads as a train of sub-parallel
# ranges rather than one dome.
_RIDGE_OCTAVES, _RIDGE_FREQ = 5, 8.0
_RIDGE_SHARPEN = 1.3
_RIDGE_GAIN = 1.4
# The ridge field is sampled through a heavier domain warp than the rest of the relief, so
# its otherwise near-parallel crest lines bend into curved, branching ranges.
_RIDGE_WARP_STRENGTH = 0.28

# Plateaus: mask coverage / edge as for belts, plus the number of discrete terrace levels
# the plateau surface is quantized toward and the fraction of one level's height the riser
# between two levels occupies (small -> abrupt, cliff-like edges).
_PLATEAU_MASK_OCTAVES, _PLATEAU_MASK_FREQ = 2, 1.0
_PLATEAU_FIELD_OCTAVES, _PLATEAU_FIELD_FREQ = 3, 1.3
_PLATEAU_COVERAGE = 0.14
_PLATEAU_EDGE = 0.22
_PLATEAU_LEVELS = 3
_PLATEAU_RISER = 0.20

# Oceanic abyssal hills (std, in units) and seamount chains (mask coverage / edge, and the
# crest gain in units above the surrounding hills). Kept sparse and low so
# `test_deep_ocean_interior_stays_deep` and the bathymetry margin passes are unaffected
# across the open ocean.
_ABYSSAL_HILL_OCTAVES, _ABYSSAL_HILL_FREQ = 3, 12.0
_ABYSSAL_HILL_AMPLITUDE = 0.18
_SEAMOUNT_MASK_OCTAVES, _SEAMOUNT_MASK_FREQ = 2, 2.2
_SEAMOUNT_OCTAVES, _SEAMOUNT_FREQ = 3, 13.0
_SEAMOUNT_COVERAGE = 0.03
_SEAMOUNT_EDGE = 0.30
_SEAMOUNT_GAIN = 0.3

# Number of points on the internal Fibonacci-sphere lattice used once per generator to
# calibrate every field's arbitrary output scale into the concrete coverage fractions /
# target std above. Coarse on purpose -- it only needs to be statistically representative.
_CALIBRATION_POINTS = 4096


def _fibonacci_sphere(n: int) -> np.ndarray:
    """`n` roughly equally spaced unit vectors on the sphere (the golden-spiral lattice) --
    a fixed, deterministic sample set for calibrating a noise field's scale."""
    i = np.arange(n) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    phi = np.pi * (1.0 + np.sqrt(5.0)) * i
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=-1)


def _smoothstep(lo: float, hi: float, x: np.ndarray) -> np.ndarray:
    """Hermite smoothstep: 0 at/below `lo`, 1 at/above `hi`, an S-curve with zero slope at
    both ends between (so a mask built from it ramps in with no kink)."""
    t = np.clip((x - lo) / (hi - lo), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _terrace(x: np.ndarray, n_levels: int, riser: float) -> np.ndarray:
    """Quantize `x` (values in roughly [-1, 1]) toward `n_levels` flat levels, connected by
    risers that occupy only `riser` (0..1) of each level's height -- so most of the range is
    dead-flat plateau surface and the transitions between levels are abrupt steps rather
    than a linear ramp. Returns values back in roughly [-1, 1]."""
    f = (np.clip(x, -1.0, 1.0) * 0.5 + 0.5) * n_levels
    level = np.floor(np.minimum(f, n_levels - 1e-9))
    frac = f - level
    step = _smoothstep(0.5 - riser / 2.0, 0.5 + riser / 2.0, frac)
    return ((level + step) / n_levels) * 2.0 - 1.0


class _CalibratedMask:
    """A `SphereNoise` field turned into a [0, 1] coverage mask: at construction it samples
    the field on the fixed `_fibonacci_sphere` lattice and picks the threshold whose
    upper tail covers exactly `coverage` of the sphere, so the mask's area fraction is
    controlled regardless of the underlying field's arbitrary output scale. `edge` is the
    smoothstep width past that threshold, as a fraction of the field's own sampled range."""

    def __init__(self, noise: SphereNoise, coverage: float, edge: float, lattice: np.ndarray):
        self._noise = noise
        vals = noise.sample(lattice)
        self._lo = float(np.quantile(vals, 1.0 - coverage))
        spread = float(vals.max() - vals.min()) or 1.0
        self._hi = self._lo + edge * spread

    def __call__(self, xyz: np.ndarray) -> np.ndarray:
        return _smoothstep(self._lo, self._hi, self._noise.sample(xyz))


class FractalTexture:
    """A cheap domain-warped multi-octave fractal in roughly [-1, 1], with no per-instance
    calibration sweep. For adding believable texture to crust that is grown or claimed
    mid-simulation (`LithospherePlate._claim_adjacent_territory`), where a full
    `ContinentalRelief` -- orogenic belts, plateaus, calibration -- is neither wanted (that
    crust is an extension of an already-shaped plate, not a fresh orogen) nor worth the
    cost. Duck-types `SphereNoise` via `sample()`. Deterministic in `rng`."""

    def __init__(self, rng: np.random.Generator):
        self._warp = [SphereNoise(rng, octaves=_WARP_OCTAVES, base_freq=_WARP_FREQ) for _ in range(3)]
        self._base = SphereNoise(rng, octaves=5, base_freq=2.4)

    def sample(self, xyz: np.ndarray) -> np.ndarray:
        w = np.stack([n.sample(xyz) for n in self._warp], axis=-1)
        return self._base.sample(np.asarray(xyz, dtype=float) + _WARP_STRENGTH * w)


class ContinentalRelief:
    """Relief field for continental crust, duck-typing `SphereNoise` via `sample()`.
    Deterministic in `rng`.

    `sample(xyz)` -- a low-frequency, domain-warped fractal field renormalized to the same
    standard deviation as the old `SphereNoise(octaves=4, base_freq=2.5)`. This is the ONLY
    part `plates._land_noise_threshold` looks at and the only part that decides which nodes
    end up above sea level, so land extent and the coarse-sample quantile behave exactly as
    before.

    `uplift(xyz)` -- a NON-NEGATIVE field, in the same units as `sample()`, carrying the
    orogenic belts (a low-frequency belt mask times sharpened ridged noise) and the
    plateaus (a low-frequency plateau mask times a terraced surface). Generation adds this
    only where `sample()` is already comfortably above the land threshold (see
    `lithosphere_plate.generate_plates`), so it can raise mountains and plateaus on land but
    never push a marine node up into land or a coastal node down into the sea -- the land
    set is identical with or without it.
    """

    def __init__(
        self,
        rng: np.random.Generator,
        *,
        orogenic_units: float,
        plateau_units: float,
        plateau_relief_units: float,
    ):
        self._orogenic_units = orogenic_units
        self._plateau_units = plateau_units
        self._plateau_relief_units = plateau_relief_units

        # Fixed draw order -- determinism depends on nothing between here reordering.
        self._warp = [SphereNoise(rng, octaves=_WARP_OCTAVES, base_freq=_WARP_FREQ) for _ in range(3)]
        self._base = SphereNoise(rng, octaves=_BASE_OCTAVES, base_freq=_BASE_FREQ)
        self._regional = SphereNoise(rng, octaves=_REGIONAL_OCTAVES, base_freq=_REGIONAL_FREQ)
        belt_noise = SphereNoise(rng, octaves=_BELT_OCTAVES, base_freq=_BELT_FREQ)
        self._ridge = SphereNoise(rng, octaves=_RIDGE_OCTAVES, base_freq=_RIDGE_FREQ)
        plateau_mask_noise = SphereNoise(rng, octaves=_PLATEAU_MASK_OCTAVES, base_freq=_PLATEAU_MASK_FREQ)
        self._plateau_field = SphereNoise(rng, octaves=_PLATEAU_FIELD_OCTAVES, base_freq=_PLATEAU_FIELD_FREQ)

        lattice = _fibonacci_sphere(_CALIBRATION_POINTS)
        self._belt_mask = _CalibratedMask(belt_noise, _BELT_COVERAGE, _BELT_EDGE, lattice)
        self._plateau_mask = _CalibratedMask(plateau_mask_noise, _PLATEAU_COVERAGE, _PLATEAU_EDGE, lattice)

        raw_std = float(np.std(self._raw_sample(lattice))) or 1.0
        self._sample_renorm = _TARGET_SAMPLE_STD / raw_std
        ridge_std = float(np.std(self._ridge.sample(lattice))) or 1.0
        self._ridge_scale = 1.0 / (2.2 * ridge_std)  # maps the ridge field to roughly [-1, 1]

    def _warp_vec(self, xyz: np.ndarray) -> np.ndarray:
        return np.stack([n.sample(xyz) for n in self._warp], axis=-1)

    def _warped(self, xyz: np.ndarray) -> np.ndarray:
        return np.asarray(xyz, dtype=float) + _WARP_STRENGTH * self._warp_vec(xyz)

    def _raw_sample(self, xyz: np.ndarray) -> np.ndarray:
        return _BASE_WEIGHT * self._base.sample(self._warped(xyz)) + _REGIONAL_WEIGHT * self._regional.sample(xyz)

    def sample(self, xyz: np.ndarray) -> np.ndarray:
        return self._raw_sample(xyz) * self._sample_renorm

    def uplift(self, xyz: np.ndarray) -> np.ndarray:
        xyz = np.asarray(xyz, dtype=float)
        warp = self._warp_vec(xyz)
        warped = xyz + _WARP_STRENGTH * warp
        ridge_warped = xyz + _RIDGE_WARP_STRENGTH * warp

        belt = _BELT_FLOOR + (1.0 - _BELT_FLOOR) * self._belt_mask(warped)
        ridged = np.clip(1.0 - np.abs(self._ridge.sample(ridge_warped) * self._ridge_scale), 0.0, 1.0) ** _RIDGE_SHARPEN
        mountains = belt * _RIDGE_GAIN * ridged * self._orogenic_units

        pmask = self._plateau_mask(xyz) * (1.0 - 0.5 * belt)  # let a belt win where the two overlap
        surface = 0.5 * (_terrace(self._plateau_field.sample(warped), _PLATEAU_LEVELS, _PLATEAU_RISER) + 1.0)
        plateau = pmask * (self._plateau_units + self._plateau_relief_units * surface)

        return np.clip(mountains + plateau, 0.0, None)


class OceanicRelief:
    """Relief field for oceanic crust, duck-typing `SphereNoise` via `sample()`.
    Deterministic in `rng`. Fine abyssal-hill texture everywhere, plus sparse seamount
    chains -- gentle enough that the open abyssal plain is unchanged at render resolution
    and only the chains themselves stand out."""

    def __init__(self, rng: np.random.Generator):
        self._hills = SphereNoise(rng, octaves=_ABYSSAL_HILL_OCTAVES, base_freq=_ABYSSAL_HILL_FREQ)
        seamount_mask_noise = SphereNoise(rng, octaves=_SEAMOUNT_MASK_OCTAVES, base_freq=_SEAMOUNT_MASK_FREQ)
        self._seamount = SphereNoise(rng, octaves=_SEAMOUNT_OCTAVES, base_freq=_SEAMOUNT_FREQ)

        lattice = _fibonacci_sphere(_CALIBRATION_POINTS)
        self._seamount_mask = _CalibratedMask(seamount_mask_noise, _SEAMOUNT_COVERAGE, _SEAMOUNT_EDGE, lattice)
        hill_std = float(np.std(self._hills.sample(lattice))) or 1.0
        self._hill_renorm = _ABYSSAL_HILL_AMPLITUDE / hill_std
        peak_std = float(np.std(self._seamount.sample(lattice))) or 1.0
        self._peak_scale = 1.0 / (2.2 * peak_std)

    def sample(self, xyz: np.ndarray) -> np.ndarray:
        hills = self._hills.sample(xyz) * self._hill_renorm
        peaks = np.clip(1.0 - np.abs(self._seamount.sample(xyz) * self._peak_scale), 0.0, 1.0) ** 2
        return hills + _SEAMOUNT_GAIN * self._seamount_mask(xyz) * peaks

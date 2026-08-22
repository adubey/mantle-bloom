"""Periodic line regularization.

Per-step boundary evolution (boundary.py) only ever touches the two ends of a line --
inserting at target spacing when growing, deleting when shrinking -- so interior spacing
stays regular on its own. What it can't fix is spacing that's drifted at a *transform*
boundary (nodes sheared along the line without insertion/deletion) or after several steps'
worth of end-growth at a slightly different rate than the line's original spacing. This
module re-derives a fresh evenly-spaced node set spanning each line's *existing* extent
(the two endpoints are preserved exactly -- regularizing never changes where a line's
physical edge is, only how regularly it's sampled) and interpolates elevation onto it.

`spacing_rad` (default plates.TARGET_LINE_SPACING_RAD, the reference density) should always
be `plates.line_spacing_rad(world.node_density)` in practice -- regularize_world_lines
computes it once per call and threads it down. Without this, a world generated at a
non-default density would regularize itself right back down to the reference density the
first time any line's spacing drifted enough to trigger this pass (every
REGULARIZE_INTERVAL_STEPS steps) -- confirmed directly as the failure mode that made a
"just build denser lines at generation" version of a density option pointless within a
handful of steps."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .plates import TARGET_LINE_SPACING_RAD, ElevationLine, Plate, line_spacing_rad

if TYPE_CHECKING:
    from .world import World

REGULARIZE_INTERVAL_STEPS = 5
IRREGULARITY_TOLERANCE = 1.5  # regularize a line if any gap exceeds this multiple of target

# Self-affine scaling exponent used by _crumple_elevation below: real terrain roughened by
# compressing a profile horizontally by k doesn't just get resampled at the new spacing, its
# vertical amplitude grows by roughly k**-CRUMPLE_HURST_EXPONENT (a Hurst exponent -- 0.5 is
# the standard "random walk" / Brownian terrain default used when no better estimate of a
# specific landscape's roughness is available). This is what makes the vulcanism-driven
# density increase that triggers crumpling look like real compression -- ridges pushed
# together get taller, not just thinned out -- rather than plain decimation, which would
# leave peak/valley heights untouched and only make the line coarser.
CRUMPLE_HURST_EXPONENT = 0.5


def needs_regularizing(line: ElevationLine, spacing_rad: float = TARGET_LINE_SPACING_RAD) -> bool:
    if len(line.theta) < 3:
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
    # The fit is a smoothing regression, not an exact interpolant, so it can drift slightly
    # from the original data even at x=0/x=n-1 -- force the two ends back to the real
    # original values so a crumpled line still butts up exactly against its neighbors'
    # elevation at the endpoints regularize_line preserves the position of.
    crumpled[0] = elevation[0]
    crumpled[-1] = elevation[-1]
    return crumpled


def regularize_line(line: ElevationLine, spacing_rad: float = TARGET_LINE_SPACING_RAD) -> ElevationLine:
    if len(line.theta) < 3:
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
    if n < len(line.theta):
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


def regularize_plate_lines(plate: Plate, spacing_rad: float = TARGET_LINE_SPACING_RAD) -> None:
    plate.lines = [
        regularize_line(line, spacing_rad) if needs_regularizing(line, spacing_rad) else line for line in plate.lines
    ]


def regularize_world_lines(world: "World") -> None:
    spacing_rad = line_spacing_rad(world.node_density)
    for plate in world.plates:
        regularize_plate_lines(plate, spacing_rad)

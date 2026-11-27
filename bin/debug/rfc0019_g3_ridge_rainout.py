"""Test a conservative windward rainout rule on captured G3 tuning fields.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3_ridge_rainout.py CAPTURE [--maps DIR]
Compares the half-strength broad transport and overflow spread with two
interception strengths. Use --weak-sensitivity for the single weaker arm,
--mid-sensitivity for the intermediate 5% arm, or --dry-sensitivity to limit
the 10% arm by the control donor moisture. Use --strength-sensitivity to
compare all four strengths against the same control. Use
--coast-rainout-sensitivity for the preregistered, tuning-only coast sink screen.
Use --upslope-sensitivity to compare the original ridge rule with rainout on
every uphill land face using a moisture-flux upslope approximation.
No saved world is loaded or stepped.

Each arm also reports native-grid rainfall along fixed, flow-aligned transects centered
on the highest native land cell inside each eligible coarse crest cell. Offsets span
three native cells upwind through three downwind; these are diagnostic samples, not
independent mountain ranges or an acceptance threshold.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import scipy
import PIL

import rfc0019_g3_compare as g3
from rfc0019_coarse_moisture import (
    FaceInterception, TransportFaces, aggregate_inputs, build_transport_rates, distribute_rain,
    solve_bounded_budget, solve_budget,
)
from rfc0019_g3_archive_cap import (
    _belt_contrasts,
    _belt_row_shares,
    _broad_meridional_capacity,
    _convergence_correction,
    _distribute_overflow,
    _signed_faces,
)
from rfc0019_g3_bounded_column import COLD_COAST_DISTANCE_KM
from rfc0019_g3_mixing import GRID_SHAPE, RAIN_RATE, TUNING_SAVES
from rfc0019_g3_mixing import _regional_metrics

RIDGE_RISE_M = 250.0
RISE_SCALE_M = 1000.0
WATER_VAPOR_GAS_CONSTANT_J_KG_K = 461.5
LATENT_HEAT_VAPORIZATION_J_KG = 2.5e6
MOIST_LAPSE_RATE_K_PER_KM = 6.0
CONDENSATE_FALLOUT_TIME_S = 6.0 * 3600.0
DRY_SCALE_KG_M2 = 10.0
MIN_REGION_AREA_FRACTION = 1e-3
MIN_RIDGE_CROSSINGS = 100
MIN_DISTINCT_CRESTS = 25
PROTOCOL_ID = "g3-three-round-plan-2026-11-18-v2"
SOURCE_FILES = (
    "backend/app/__init__.py",
    "backend/app/hydroclimate/__init__.py",
    "backend/app/hydroclimate/climate.py",
    "backend/app/hydroclimate/hydrology.py",
    "backend/constraints.txt",
    "bin/debug/rfc0019_coarse_moisture.py",
    "bin/debug/rfc0019_g3_archive_cap.py",
    "bin/debug/rfc0019_g3_bounded_column.py",
    "bin/debug/rfc0019_g3_compare.py",
    "bin/debug/rfc0019_g3_mixing.py",
    "bin/debug/rfc0019_g3_ridge_rainout.py",
    "bin/debug/score_rfc0019_g3_tuning_round.py",
    "bin/debug/rfc0019_transport_prototypes.py",
)
ARMS = {"control": 0.0, "ridge_010": 0.10, "ridge_025": 0.25}
WEAK_ARMS = {"control": 0.0, "ridge_0025": 0.025}
MID_ARMS = {"control": 0.0, "ridge_005": 0.05}
DRY_ARMS = {"control": 0.0, "ridge_010_dry": 0.10}
STRENGTH_ARMS = {
    "control": 0.0, "ridge_0025": 0.025, "ridge_005": 0.05,
    "ridge_010": 0.10, "ridge_025": 0.25,
}


@dataclass(frozen=True)
class Crossings:
    face: np.ndarray
    toward_second: np.ndarray
    donor: np.ndarray
    crest: np.ndarray
    lee: np.ndarray
    rise_m: np.ndarray


def _directed_pair_flow(flow_by_pair, origin, destination):
    """Return signed flow in origin→destination order from a canonical pair map."""
    key = (min(origin, destination), max(origin, destination))
    canonical_flow = flow_by_pair.get(key, 0.0)
    return canonical_flow if origin < destination else -canonical_flow


def _assert_pair_flow_orientation():
    flow = {(2, 7): 3.5}
    if _directed_pair_flow(flow, 2, 7) != 3.5 or _directed_pair_flow(flow, 7, 2) != -3.5:
        raise AssertionError("canonical ridge-face flow lost its direction")


def _ridge_crossings(grid, faces, elevation):
    """Follow each directed face across a crest and onward to the next lower land cell."""
    first, second = faces.first, faces.second
    nlat, nlon = grid.coarse_shape
    row, col = np.divmod(grid.coarse_cell, nlon)
    east_axis = row[first] == row[second]
    signed = faces.forward - faces.reverse
    toward_second = signed > 0.0
    step = np.sign(signed).astype(int)
    entry = np.where(toward_second, first, second)
    crest = np.where(toward_second, second, first)
    lee_row = row[crest] + step * ~east_axis
    lee_col = (col[crest] + step * east_axis) % nlon
    land_reservoir = np.full(nlat * nlon, -1, dtype=np.intp)
    land_reservoir[grid.coarse_cell[~grid.is_ocean]] = np.flatnonzero(~grid.is_ocean)
    lee_cell = np.clip(lee_row, 0, nlat - 1) * nlon + lee_col
    lee = land_reservoir[lee_cell]

    # Check the exact crest-to-lee face, summing its native boundary segments.
    face_pairs = np.column_stack((first, second))
    canonical_pairs = np.sort(face_pairs, axis=1)
    unique_pairs, pair_index = np.unique(canonical_pairs, axis=0, return_inverse=True)
    orientation = np.where(first < second, 1.0, -1.0)
    pair_flow = np.bincount(pair_index, weights=signed * orientation)
    flow_by_pair = {tuple(pair): value for pair, value in zip(unique_pairs, pair_flow, strict=True)}
    exit_first = np.where(step > 0, crest, lee)
    exit_second = np.where(step > 0, lee, crest)
    exit_flow = np.fromiter(
        (_directed_pair_flow(flow_by_pair, int(a), int(b))
         for a, b in zip(exit_first, exit_second, strict=True)),
        dtype=float, count=step.size,
    ) * step
    rise = elevation[crest] - elevation[entry]
    fall = elevation[crest] - elevation[np.maximum(lee, 0)]
    land = ~grid.is_ocean
    valid = ((step != 0) & (lee_row >= 0) & (lee_row < nlat) & (lee >= 0) & (exit_flow > 0)
             & land[entry] & land[crest] & land[np.maximum(lee, 0)]
             & (rise >= RIDGE_RISE_M) & (fall >= RIDGE_RISE_M))
    return Crossings(np.flatnonzero(valid), toward_second[valid], entry[valid], crest[valid], lee[valid], rise[valid])


def _interception(crossings, face_count, strength):
    forward = np.zeros(face_count)
    reverse = np.zeros(face_count)
    fraction = -strength * np.expm1(-crossings.rise_m / RISE_SCALE_M)
    east_or_north = crossings.toward_second
    forward[crossings.face[east_or_north]] = fraction[east_or_north]
    reverse[crossings.face[~east_or_north]] = fraction[~east_or_north]
    return FaceInterception(forward, reverse)


def _saturation_scale_height_m(temperature_k, moist_lapse_rate_k_per_km):
    """Approximate saturation-density scale from Clausius-Clapeyron."""
    if moist_lapse_rate_k_per_km <= 0.0:
        raise ValueError("moist lapse rate must be positive")
    lapse_rate_k_per_m = moist_lapse_rate_k_per_km / 1000.0
    return WATER_VAPOR_GAS_CONSTANT_J_KG_K * temperature_k**2 / (
        LATENT_HEAT_VAPORIZATION_J_KG * lapse_rate_k_per_m
    )


def _upslope_interception(grid, faces, elevation, temperature_k, moist_lapse_rate_k_per_km):
    """Apply the basic upslope rainout fraction using transported moisture flux.

    Assumes saturated ascent and immediate fallout. The transported face flux
    supplies available water; the thermodynamic scale varies with face temperature.
    """
    if not np.all(np.isfinite(temperature_k) & (temperature_k > 0)):
        raise ValueError("expected positive finite reservoir temperatures in kelvin")
    first, second = faces.first, faces.second
    both_land = ~grid.is_ocean[first] & ~grid.is_ocean[second]
    rise = elevation[second] - elevation[first]
    forward_rise = np.where(both_land, np.maximum(rise, 0.0), 0.0)
    reverse_rise = np.where(both_land, np.maximum(-rise, 0.0), 0.0)
    face_temperature = 0.5 * (temperature_k[first] + temperature_k[second])
    scale_height = _saturation_scale_height_m(face_temperature, moist_lapse_rate_k_per_km)
    return FaceInterception(
        -np.expm1(-forward_rise / scale_height),
        -np.expm1(-reverse_rise / scale_height),
    )


def _upslope_receiver_mask(grid, faces, elevation):
    """Mark land reservoirs receiving any directed moisture over rising terrain."""
    first, second = faces.first, faces.second
    land = ~grid.is_ocean
    mask = np.zeros(grid.area.size, dtype=bool)
    forward = land[first] & land[second] & (elevation[second] > elevation[first]) & (faces.forward > 0)
    reverse = land[first] & land[second] & (elevation[first] > elevation[second]) & (faces.reverse > 0)
    mask[second[forward]] = True
    mask[first[reverse]] = True
    return mask


def _upslope_rain_metrics(grid, rain, control, receiver_mask):
    """Compare native rain over all land reservoirs receiving upslope flux."""
    native_mask = receiver_mask[grid.native_ids][::-1]
    native_land = ~grid.is_ocean[grid.native_ids][::-1]
    mask = native_mask & native_land
    area = grid.native_area[::-1]
    weights = area[mask]
    if weights.size == 0:
        return {"receiver_area_fraction": 0.0, "candidate_mean_mm_yr": None,
                "control_mean_mm_yr": None, "change_mm_yr": None}
    candidate_mean = float(np.average(rain[mask], weights=weights))
    control_mean = float(np.average(control[mask], weights=weights))
    land_area = float(area[native_land].sum())
    return {
        "receiver_area_fraction": float(weights.sum() / land_area) if land_area else 0.0,
        "candidate_mean_mm_yr": candidate_mean,
        "control_mean_mm_yr": control_mean,
        "change_mm_yr": candidate_mean - control_mean,
        "receiver_cells": int(mask.sum()),
    }


def _paired_rain(crossings, rain, control, grid):
    if not crossings.face.size:
        return {"crossings": 0}
    crest, lee = crossings.crest, crossings.lee
    native_area = grid.native_area.ravel()
    on_path = np.bincount(
        grid.native_ids.ravel(), weights=rain[::-1].ravel() * native_area, minlength=grid.area.size,
    ) / grid.area
    baseline = np.bincount(
        grid.native_ids.ravel(), weights=control[::-1].ravel() * native_area, minlength=grid.area.size,
    ) / grid.area
    crest_change = on_path[crest] - baseline[crest]
    lee_change = on_path[lee] - baseline[lee]
    return {
        "crossings": int(crest.size),
        "distinct_crests": int(np.unique(crest).size),
        "median_crest_rain_mm_yr": float(np.median(on_path[crest])),
        "median_lee_rain_mm_yr": float(np.median(on_path[lee])),
        "median_crest_minus_lee_mm_yr": float(np.median(on_path[crest] - on_path[lee])),
        "median_contrast_gain_mm_yr": float(np.median(crest_change - lee_change)),
        "median_lee_change_mm_yr": float(np.median(lee_change)),
        "crest_wetter_and_lee_drier_share": float(np.mean((crest_change > 0) & (lee_change < 0))),
    }


def _native_crossing_profiles(crossings, rain, control, grid, fields, max_offset=3):
    """Summarize native-grid rain along unique flow-aligned crest transects."""
    if not crossings.face.size:
        return {"offsets_native_cells": [], "profile_count": 0, "profiles": [], "summary": []}

    nlat, nlon = grid.coarse_shape
    native_height, native_width = rain.shape
    factor = grid.factor
    if (native_height, native_width) != (nlat * factor, nlon * factor):
        raise ValueError("native rain grid does not match the coarse crossing grid")

    reservoir_cell = grid.coarse_cell

    native_ids = grid.native_ids.ravel()
    member_order = np.argsort(native_ids, kind="stable")
    member_counts = np.bincount(native_ids, minlength=grid.area.size)
    member_starts = np.concatenate(([0], np.cumsum(member_counts)))
    elevation = fields.elevation_m[::-1].ravel()
    native_ocean = fields.is_ocean[::-1].ravel()
    native_rain = rain[::-1].ravel()
    native_control = control[::-1].ravel()
    native_area = grid.native_area.ravel()

    profiles = []
    seen = set()
    for entry, crest, lee in zip(crossings.donor, crossings.crest, crossings.lee, strict=True):
        entry_cell = reservoir_cell[entry]
        crest_cell = reservoir_cell[crest]
        lee_cell = reservoir_cell[lee]
        erow, _ = divmod(int(entry_cell), nlon)
        crow, ccol = divmod(int(crest_cell), nlon)
        lrow, lcol = divmod(int(lee_cell), nlon)
        east_axis = erow == crow
        if east_axis:
            step = 1 if (lcol - ccol) % nlon == 1 else -1
            row_step, col_step = 0, step
            axis = "east_west"
        else:
            step = 1 if lrow > crow else -1
            row_step, col_step = step, 0
            axis = "north_south"
        key = (int(crest_cell), axis, step)
        if key in seen:
            continue
        seen.add(key)

        reservoir_start = member_starts[crest]
        reservoir_end = member_starts[crest + 1]
        members = member_order[reservoir_start:reservoir_end]
        if members.size == 0 or np.any(native_ocean[members]):
            raise ValueError("crossing crest reservoir has no all-land native footprint")
        peak = int(members[np.argmax(elevation[members])])
        peak_row, peak_col = divmod(peak, native_width)
        samples = {}
        for offset in range(-max_offset, max_offset + 1):
            row = peak_row + offset * row_step
            col = (peak_col + offset * col_step) % native_width
            if row < 0 or row >= native_height:
                samples[str(offset)] = None
                continue
            native_index = row * native_width + col
            if native_ocean[native_index]:
                samples[str(offset)] = None
                continue
            samples[str(offset)] = {
                "rain_mm_yr": float(native_rain[native_index]),
                "control_rain_mm_yr": float(native_control[native_index]),
                "candidate_minus_control_mm_yr": float(
                    native_rain[native_index] - native_control[native_index]
                ),
                "cell_area_m2": float(native_area[native_index]),
            }
        profiles.append({
            "crest_coarse_cell": int(crest_cell),
            "crest_reservoir": int(crest),
            "axis": axis,
            "downwind_step": int(step),
            "native_peak_row_solver_order": int(peak_row),
            "native_peak_col": int(peak_col),
            "samples_by_offset_native_cells": samples,
        })

    summary = []
    for offset in range(-max_offset, max_offset + 1):
        rows = [p["samples_by_offset_native_cells"][str(offset)] for p in profiles]
        valid = [row for row in rows if row is not None]
        rain_values = np.asarray([row["rain_mm_yr"] for row in valid], dtype=float)
        control_values = np.asarray([row["control_rain_mm_yr"] for row in valid], dtype=float)
        changes = np.asarray([row["candidate_minus_control_mm_yr"] for row in valid], dtype=float)
        weights = np.asarray([row["cell_area_m2"] for row in valid], dtype=float)
        summary.append({
            "offset_native_cells": offset,
            "valid_profiles": len(valid),
            "candidate_median_rain_mm_yr": (
                float(np.median(rain_values)) if valid else None
            ),
            "candidate_p10_rain_mm_yr": (
                float(np.percentile(rain_values, 10)) if valid else None
            ),
            "candidate_p90_rain_mm_yr": (
                float(np.percentile(rain_values, 90)) if valid else None
            ),
            "candidate_area_weighted_mean_rain_mm_yr": (
                float(np.average(rain_values, weights=weights)) if valid else None
            ),
            "control_median_rain_mm_yr": (
                float(np.median(control_values)) if valid else None
            ),
            "control_p10_rain_mm_yr": (
                float(np.percentile(control_values, 10)) if valid else None
            ),
            "control_p90_rain_mm_yr": (
                float(np.percentile(control_values, 90)) if valid else None
            ),
            "median_candidate_minus_control_mm_yr": (
                float(np.median(changes)) if valid else None
            ),
        })
    windward_contrasts = []
    control_contrasts = []
    contrast_gains = []
    for profile in profiles:
        samples = profile["samples_by_offset_native_cells"]
        windward = [samples[str(offset)] for offset in range(-max_offset, 0)]
        leeward = [samples[str(offset)] for offset in range(1, max_offset + 1)]
        windward = [sample for sample in windward if sample is not None]
        leeward = [sample for sample in leeward if sample is not None]
        if min(len(windward), len(leeward)) < 2:
            continue
        candidate_contrast = (
            np.mean([sample["rain_mm_yr"] for sample in windward])
            - np.mean([sample["rain_mm_yr"] for sample in leeward])
        )
        control_contrast = (
            np.mean([sample["control_rain_mm_yr"] for sample in windward])
            - np.mean([sample["control_rain_mm_yr"] for sample in leeward])
        )
        windward_contrasts.append(float(candidate_contrast))
        control_contrasts.append(float(control_contrast))
        contrast_gains.append(float(candidate_contrast - control_contrast))

    def summarize_contrasts(values):
        if not values:
            return {"median_mm_yr": None, "p10_mm_yr": None, "p90_mm_yr": None}
        return {
            "median_mm_yr": float(np.median(values)),
            "p10_mm_yr": float(np.percentile(values, 10)),
            "p90_mm_yr": float(np.percentile(values, 90)),
        }

    return {
        "offsets_native_cells": list(range(-max_offset, max_offset + 1)),
        "profile_count": len(profiles),
        "crest_selection": "highest native land cell in each coarse transport crest cell",
        "positive_offsets_follow_transport_downwind": True,
        "north_south_rows_use_solver_order": True,
        "profiles": profiles,
        "summary": summary,
        "multi_cell_windward_minus_lee": {
            "side_offsets_native_cells": {
                "windward": list(range(-max_offset, 0)),
                "lee": list(range(1, max_offset + 1)),
            },
            "minimum_valid_cells_per_side": 2,
            "eligible_profiles": len(windward_contrasts),
            "candidate": summarize_contrasts(windward_contrasts),
            "control": summarize_contrasts(control_contrasts),
            "candidate_minus_control": summarize_contrasts(contrast_gains),
            "candidate_fraction_windward_wetter": (
                float(np.mean(np.asarray(windward_contrasts) > 0))
                if windward_contrasts else None
            ),
            "fraction_with_positive_contrast_gain": (
                float(np.mean(np.asarray(contrast_gains) > 0)) if contrast_gains else None
            ),
        },
    }


def _rain_metrics(fields, rain, area, distance, freeze):
    cold_land = ~fields.is_ocean & (fields.air_temperature_c < freeze)
    coastal = distance <= COLD_COAST_DISTANCE_KM
    metrics = _regional_metrics(fields, rain, area, distance, freeze)
    for name, mask in (("cold_coast", cold_land & coastal), ("cold_inland", cold_land & ~coastal)):
        metrics[name + "_rain_mm_yr"] = g3._weighted_mean(rain[mask], area[mask])
    return metrics | {
        "ocean_belt_contrasts": _belt_contrasts(fields, rain, area),
        "ocean_belt_wettest_row_shares": _belt_row_shares(fields, rain, area),
    }


def _tuning_screen(candidate, control, paired_control):
    """Score one paired tuning case against the preregistered pattern limits."""
    support = candidate["geographic_support"]
    distances = candidate["land_rain_by_distance_to_ocean_km"]
    near = distances["0_to_500"]["rain_mm_yr"]
    inland = distances["1500_plus"]["rain_mm_yr"]
    inland_ratio = None if near is None or near <= 0 or inland is None else inland / near
    belt = candidate["ocean_belt_contrasts"]["0"]["ratio"]
    ridge = candidate["ridge_transects"]
    ridge_share = ridge.get("crest_wetter_and_lee_drier_share")
    ridge_gain = ridge.get("median_contrast_gain_mm_yr")
    cold_coast = candidate.get("cold_coast_rain_mm_yr")
    control_inland = control["land_rain_by_distance_to_ocean_km"]["1500_plus"]["rain_mm_yr"]
    inland_vs_control = None if control_inland in (None, 0) or inland is None else inland / control_inland
    land_share_delta = candidate["land_rain_share"] - control["land_rain_share"]

    shortfall = {
        "equatorial_belt": (
            None if not support["equatorial_belt"]["evaluable"]
            else 1.0 if belt is None else max(0.0, (1.5 - belt) / 1.5)
        ),
        "inland_coast": (
            None if not support["inland_coast"]["evaluable"]
            else 1.0 if inland_ratio is None else max(0.0, (inland_ratio - 0.6) / 0.6)
        ),
        "ridge_direction_share": (
            None if not support["ridge_shadow"]["evaluable"]
            else 1.0 if ridge_share is None else max(0.0, (0.6 - ridge_share) / 0.6)
        ),
        "ridge_contrast_gain": (
            None if not support["ridge_shadow"]["evaluable"]
            else 1.0 if ridge_gain is None or ridge_gain <= 0 else 0.0
        ),
        "cold_coast_positive": (
            None if not support["cold_coast"]["evaluable"]
            else 1.0 if cold_coast is None or cold_coast <= 0 else 0.0
        ),
    }
    evaluable = [value for value in shortfall.values() if value is not None]
    safeguards = {
        "positive_control_inland_rain": control_inland is not None and control_inland > 0,
        "inland_at_least_80_percent_of_control": inland_vs_control is not None and inland_vs_control >= 0.8,
        "land_share_within_2pp_of_control": abs(land_share_delta) <= 0.02,
    }
    hard_gates = {
        **{f"candidate_{name}": passed for name, passed in candidate["hard_gates"].items()},
        **{f"paired_control_{name}": passed for name, passed in paired_control["hard_gates"].items()},
        **{f"frozen_baseline_{name}": passed for name, passed in control["hard_gates"].items()},
    }
    eligible = all(safeguards.values()) and all(hard_gates.values()) and bool(evaluable)
    pattern_pass = bool(evaluable) and all(value == 0.0 for value in evaluable)
    return {
        "inland_to_coast_ratio": inland_ratio,
        "inland_rain_vs_control": inland_vs_control,
        "land_rain_share_delta_vs_control": land_share_delta,
        "geographic_support": support,
        "equatorial_belt_ratio": belt,
        "ridge_direction_share": ridge_share,
        "ridge_median_contrast_gain_mm_yr": ridge_gain,
        "cold_coast_rain_mm_yr": cold_coast,
        "normalized_shortfalls": shortfall,
        "minimax_score": max(evaluable) if evaluable else None,
        "hard_gates": hard_gates,
        "all_hard_gates_pass": all(hard_gates.values()),
        "safeguards": safeguards,
        "eligible": eligible,
        "pattern_pass": pattern_pass,
    }


def _geographic_support(fields, area, distance, freeze, crossings):
    """Determine metric eligibility from fixed geography, never candidate rainfall."""
    total_area = float(area.sum())
    abs_lat = np.abs(fields.lat_deg[:, None])
    land = ~fields.is_ocean
    masks = {
        "belt_center": fields.is_ocean & (abs_lat < 8.0),
        "belt_flanks": fields.is_ocean & (abs_lat >= 8.0) & (abs_lat < 20.0),
        "near_coast": land & (distance < 500.0),
        "inland": land & (distance > 1500.0),
        "cold_coast": land & (fields.air_temperature_c < freeze) & (distance <= 500.0),
    }
    fractions = {name: float(area[mask].sum() / total_area) for name, mask in masks.items()}
    ridge_count = int(crossings.face.size)
    distinct_crests = int(np.unique(crossings.crest).size)
    return {
        "area_fraction": fractions,
        "minimum_region_area_fraction": MIN_REGION_AREA_FRACTION,
        "ridge_crossings": ridge_count,
        "distinct_crests": distinct_crests,
        "equatorial_belt": {
            "evaluable": min(fractions["belt_center"], fractions["belt_flanks"])
            >= MIN_REGION_AREA_FRACTION,
        },
        "inland_coast": {
            "evaluable": min(fractions["near_coast"], fractions["inland"])
            >= MIN_REGION_AREA_FRACTION,
        },
        "ridge_shadow": {
            "evaluable": ridge_count >= MIN_RIDGE_CROSSINGS
            and distinct_crests >= MIN_DISTINCT_CRESTS,
            "minimum_crossings": MIN_RIDGE_CROSSINGS,
            "minimum_distinct_crests": MIN_DISTINCT_CRESTS,
        },
        "cold_coast": {"evaluable": fractions["cold_coast"] >= MIN_REGION_AREA_FRACTION},
    }


def _distribution_metrics(fields, rain, area):
    """Summarize global and folded-latitude rain for the G3 density check."""
    mass_rate = rain * area / g3.SECONDS_PER_YEAR
    land = ~fields.is_ocean
    total = float(mass_rate.sum())
    land_total = float(mass_rate[land].sum())
    ocean_total = float(mass_rate[~land].sum())
    bands = []
    absolute_latitude = np.abs(fields.lat_deg)
    for lower in range(0, 90, 15):
        mask = (absolute_latitude >= lower) & (absolute_latitude < lower + 15)
        band_area = area[mask]
        band_land = land[mask]
        bands.append({
            "absolute_latitude_deg": [lower, lower + 15],
            "land_mean_mm_yr": g3._weighted_mean(rain[mask][band_land], band_area[band_land]),
            "ocean_mean_mm_yr": g3._weighted_mean(rain[mask][~band_land], band_area[~band_land]),
        })
    return {
        "total_rain_kg_s": total,
        "global_mean_rain_mm_yr": float(np.sum(rain * area) / area.sum()),
        "land_rain_share": land_total / max(total, 1e-300),
        "ocean_rain_share": ocean_total / max(total, 1e-300),
        "folded_latitude_bands": bands,
    }


def _place_surface_rain(grid, budget, overflow, area, fields, freeze):
    """Place reservoir rain on native cells, spreading only saturation overflow."""
    weights = np.ones(grid.native_ids.shape)
    ordinary_rain = budget.rain - overflow
    tolerance = 100 * np.finfo(float).eps * np.maximum(budget.rain, 1.0)
    if np.any(ordinary_rain < -tolerance):
        raise AssertionError("overflow exceeds the reservoir's total rain")
    ordinary = distribute_rain(grid, np.maximum(ordinary_rain, 0.0), weights)[::-1]
    excess = distribute_rain(grid, overflow, weights)[::-1]
    surface_class = np.where(fields.is_ocean, 2, fields.air_temperature_c < freeze)
    deposition = _distribute_overflow(excess, area, surface_class)
    before = np.bincount(surface_class.ravel(), weights=excess.ravel(), minlength=3)
    after = np.bincount(surface_class.ravel(), weights=deposition.ravel(), minlength=3)
    scale = max(float(budget.rain.sum()), 1e-300)
    class_residual = float(np.abs(after - before).sum() / scale)
    placed = float((ordinary + deposition).sum())
    mass_residual = abs(placed - scale) / scale
    residual = max(class_residual, mass_residual)
    if residual > 1e-12:
        raise AssertionError(f"overflow placement did not preserve surface-class budgets: {residual:g}")
    return (ordinary + deposition) / area * g3.SECONDS_PER_YEAR, residual


def _coast_sink_profiles(grid, distance_km, length_km, alphas=(1.0, 3.0)):
    """Return coast-focused and equal-integral uniform land sink arrays."""
    if distance_km.shape != grid.native_ids.shape:
        raise ValueError("coast distance must match the archived native grid")
    if not np.isfinite(length_km) or length_km <= 0:
        raise ValueError("coast distance scale must be finite and positive")
    distance = distance_km[::-1]
    if not np.all(np.isfinite(distance) & (distance >= 0)):
        raise ValueError("coast distances must be finite and nonnegative")

    native_profile = np.exp(-distance / length_km)
    profile = np.bincount(
        grid.native_ids.ravel(),
        weights=native_profile.ravel() * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    land = ~grid.is_ocean
    profile = np.where(land, profile, 0.0)
    if not np.all(np.isfinite(profile) & (profile >= 0) & (profile <= 1)):
        raise AssertionError("coast profile must be finite and in [0, 1]")

    sinks = {"control": np.full(grid.area.shape, RAIN_RATE)}
    for alpha in alphas:
        if not np.isfinite(alpha) or alpha <= 0:
            raise ValueError("coast sink alpha must be finite and positive")
        name = f"{alpha:g}"
        coast = np.where(land, RAIN_RATE * (1.0 + alpha * profile), RAIN_RATE)
        mean_profile = float(np.average(profile[land], weights=grid.area[land]))
        uniform = np.where(land, RAIN_RATE * (1.0 + alpha * mean_profile), RAIN_RATE)
        if not np.isclose(np.dot(grid.area, coast), np.dot(grid.area, uniform), rtol=1e-14, atol=0):
            raise AssertionError("matched uniform sink must have the same integrated land sink")
        if np.any(coast[grid.is_ocean] != RAIN_RATE) or np.any(coast[land] > (1.0 + alpha) * RAIN_RATE):
            raise AssertionError("coast sink is outside its preregistered bounds")
        sinks[f"coast_{name}"] = coast
        sinks[f"uniform_{name}"] = uniform
    return profile, sinks


def _source_provenance(root):
    """Fingerprint every local source file that can affect the diagnostic results."""
    hashes = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in SOURCE_FILES
    }
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, check=True, capture_output=True, text=True,
    ).stdout.strip()
    dirty = bool(subprocess.run(
        ["git", "status", "--porcelain"], cwd=root, check=True, capture_output=True, text=True,
    ).stdout.strip())
    runtime = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pillow": PIL.__version__,
    }
    fingerprint = hashlib.sha256(json.dumps({
        "commit": commit, "dirty": dirty, "files": hashes, "runtime": runtime,
        "protocol_id": PROTOCOL_ID,
    }, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {
        "commit": commit,
        "protocol_id": PROTOCOL_ID,
        "worktree_dirty": dirty,
        "files": hashes,
        "fingerprint": fingerprint,
        "runtime": runtime,
    }


def _assert_coast_profile_orientation():
    """Catch a missing archive-to-reservoir row flip on a synthetic coast pair."""
    field_ocean = np.array([[False, False, False, False, True, True],
                            [True, True, False, False, True, True]])
    grid = aggregate_inputs(field_ocean[::-1], {}, factor=2)
    archived_distance = np.array([[0.0] * 6, [3000.0] * 6])
    profile, _ = _coast_sink_profiles(grid, archived_distance, 750.0)
    reservoir = np.full((1, 3, 2), -1, dtype=int)
    for idx, cell in enumerate(grid.coarse_cell):
        row, col = divmod(int(cell), grid.coarse_shape[1])
        reservoir[row, col, int(grid.is_ocean[idx])] = idx
    coastal_land = reservoir[0, 0, 0]
    inland_land = reservoir[0, 1, 0]
    if coastal_land < 0 or inland_land < 0 or profile[coastal_land] <= profile[inland_land]:
        raise AssertionError("coast profile row orientation is reversed")


def _validate_capture_identity(provenance, profile, parser):
    """Reject unregistered or spent-as-tuning captures before opening any NPZ arrays."""
    save = Path(provenance["save"])
    if provenance.get("heldout"):
        digest = provenance.get("save_sha1")
        active = g3.HELDOUT_SHA1.get(save.name)
        spent = g3.SPENT_HELDOUT_SHA1.get(save.name)
        if digest not in (active, spent) or digest is None:
            parser.error("held-out capture does not match a registered active or spent SHA-1")
        if profile not in (None, "upslope") and (active is None or digest != active):
            parser.error("spent held-out captures cannot be used for tuning")
        if profile not in (None, "upslope"):
            parser.error("held-out captures cannot be used for any tuning sensitivity")
    else:
        expected = TUNING_SAVES.get(save.name)
        if expected is None or provenance.get("save_sha256") != expected:
            parser.error("capture does not match an exact registered tuning save")


def _validate_upslope_parameters(args, heldout, parser):
    if not 4.0 <= args.moist_lapse_rate_k_per_km <= 8.0:
        parser.error("moist lapse rate must be between 4 and 8 K/km")
    if heldout and args.profile == "upslope" and args.moist_lapse_rate_k_per_km != MOIST_LAPSE_RATE_K_PER_KM:
        parser.error("held-out upslope comparisons use the fixed 6 K/km moist-lapse assumption")
    if args.condensate_fallout_hours <= 0:
        parser.error("condensate fallout time must be positive")
    if heldout and args.profile == "upslope" and args.condensate_fallout_hours not in (6.0, 12.0):
        parser.error("held-out delayed-fallout comparisons use fixed 6- or 12-hour travel-scale tests")


def _cold_inland_delivery(fields, distance, freeze, grid, faces, source, budget, surface_rain):
    """Separate atmospheric delivery from any spread across the coastal-distance mask."""
    cold_inland = ~fields.is_ocean & (fields.air_temperature_c < freeze) & (distance > COLD_COAST_DISTANCE_KM)
    inland = np.bincount(
        grid.native_ids.ravel(),
        weights=cold_inland[::-1].ravel() * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    inward = (inland[faces.second] - inland[faces.first]) * budget.face_flux
    atmospheric_rain = float(np.dot(budget.rain, inland))
    native_surface = surface_rain[::-1].ravel() * grid.native_area.ravel() / g3.SECONDS_PER_YEAR
    local_source = float(np.dot(source, inland))
    if abs(local_source + inward.sum() - atmospheric_rain) > 1e-12 * max(1.0, source.sum()):
        raise AssertionError("cold-inland regional budget did not close")
    inland_native = cold_inland[::-1].ravel()
    return {
        "inflow_kg_s": float(np.maximum(inward, 0).sum()),
        "outflow_kg_s": float(np.maximum(-inward, 0).sum()),
        "net_inflow_kg_s": float(inward.sum()),
        "local_source_kg_s": local_source,
        "atmospheric_rain_kg_s": atmospheric_rain,
        "spread_contribution_kg_s": float(native_surface[inland_native].sum()) - atmospheric_rain,
    }


def _wind_arrow_faces(grid, native_u, native_v):
    """Make coarse, locally directed neighbor transfers with shared coast exchange."""
    nlat, nlon = grid.coarse_shape
    res_area = grid.area
    ids = grid.native_ids.ravel()
    native_area = grid.native_area.ravel()
    u = np.bincount(ids, weights=native_u.ravel() * native_area, minlength=res_area.size) / res_area
    v = np.bincount(ids, weights=native_v.ravel() * native_area, minlength=res_area.size) / res_area
    row, col = np.divmod(grid.coarse_cell, nlon)
    class_lookup = np.full((nlat * nlon, 2), -1, dtype=np.intp)
    class_lookup[grid.coarse_cell, grid.is_ocean.astype(np.intp)] = np.arange(res_area.size)
    edge = np.linspace(-np.pi / 2, np.pi / 2, nlat + 1)
    latitude = 0.5 * (edge[:-1] + edge[1:])
    longitude_step = 2 * np.pi / nlon
    arcs = []

    def add_arc(origin, row_step, col_step, speed, distance):
        if speed == 0.0:
            return
        dest_row = row[origin] + row_step
        if not 0 <= dest_row < nlat:
            return
        dest_col = (col[origin] + col_step) % nlon
        cell = dest_row * nlon + dest_col
        surface = int(grid.is_ocean[origin])
        dest = class_lookup[cell, surface]
        if dest < 0:
            dest = class_lookup[cell, 1 - surface]
        if dest >= 0:
            arcs.append((int(origin), int(dest), abs(float(speed)) / distance))

    for origin in range(res_area.size):
        if u[origin] != 0.0:
            step = 1 if u[origin] > 0.0 else -1
            east_distance = 6_371_000.0 * np.cos(latitude[row[origin]]) * longitude_step
            add_arc(origin, 0, step, u[origin], east_distance)
        if v[origin] != 0.0:
            step = 1 if v[origin] > 0.0 else -1
            neighbor = row[origin] + step
            if 0 <= neighbor < nlat:
                north_distance = 6_371_000.0 * abs(latitude[neighbor] - latitude[row[origin]])
                add_arc(origin, step, 0, v[origin], north_distance)

    flux_by_pair = {}
    for origin, destination, rate in arcs:
        if row[origin] == row[destination]:
            eastward = (col[destination] - col[origin]) % nlon == 1
            first, second = (origin, destination) if eastward else (destination, origin)
        else:
            northward = row[destination] == row[origin] + 1
            first, second = (origin, destination) if northward else (destination, origin)
        pair = flux_by_pair.setdefault((first, second), [0.0, 0.0])
        direction = 0 if (first, second) == (origin, destination) else 1
        pair[direction] += rate * res_area[origin]

    if flux_by_pair:
        pairs = sorted(flux_by_pair)
        first = np.fromiter((pair[0] for pair in pairs), dtype=np.intp)
        second = np.fromiter((pair[1] for pair in pairs), dtype=np.intp)
        forward = np.fromiter((flux_by_pair[pair][0] for pair in pairs), dtype=float)
        reverse = np.fromiter((flux_by_pair[pair][1] for pair in pairs), dtype=float)
    else:
        first = second = np.array([], dtype=np.intp)
        forward = reverse = np.array([], dtype=float)
    coast = build_transport_rates(
        grid, np.zeros(grid.native_ids.shape), np.zeros(grid.native_ids.shape), 0.0,
    )
    has_exchange = coast.exchange > 0.0
    return TransportFaces(
        np.concatenate((first, coast.first[has_exchange])),
        np.concatenate((second, coast.second[has_exchange])),
        np.concatenate((forward, coast.forward[has_exchange])),
        np.concatenate((reverse, coast.reverse[has_exchange])),
        np.concatenate((np.zeros(forward.size), coast.exchange[has_exchange])),
    )


def _candidate_faces(grid, wind_u, wind_v, solver):
    if solver == "finite_volume":
        original = build_transport_rates(grid, wind_u, wind_v, 0.0)
        signed = original.forward - original.reverse
        correction = _convergence_correction(grid, original.first, original.second, signed)
        capacity = signed - correction + 0.5 * _broad_meridional_capacity(grid, original)
        return _signed_faces(original.first, original.second, capacity)
    if solver == "wind_arrows":
        return _wind_arrow_faces(grid, wind_u, wind_v)
    raise ValueError(f"unsupported transport solver: {solver}")


def _solve_arms(
    fields, archive, freeze, arms, *, solver="finite_volume", dry_scale=None, sink_profiles=None,
    upslope_scale_height_k_per_km=MOIST_LAPSE_RATE_K_PER_KM,
    condensate_fallout_time_s=CONDENSATE_FALLOUT_TIME_S,
):
    density_factor = fields.is_ocean.shape[0] // GRID_SHAPE[0]
    if fields.is_ocean.shape != (GRID_SHAPE[0] * density_factor, GRID_SHAPE[1] * density_factor):
        raise ValueError("capture grid must be density 1 or 4")
    grid = aggregate_inputs(fields.is_ocean[::-1], {}, factor=density_factor)
    faces = _candidate_faces(grid, fields.wind_u[::-1], fields.wind_v[::-1], solver)
    native_elevation = fields.elevation_m[::-1].ravel()
    elevation_by_reservoir = np.bincount(
        grid.native_ids.ravel(),
        weights=native_elevation * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    temperature_k_by_reservoir = np.bincount(
        grid.native_ids.ravel(),
        weights=(fields.air_temperature_c[::-1].ravel() + 273.15) * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    native_land = ~fields.is_ocean[::-1]
    native_wind_speed = np.hypot(fields.wind_u[::-1], fields.wind_v[::-1])
    mean_land_wind_m_s = float(np.average(
        native_wind_speed[native_land], weights=grid.native_area[native_land],
    ))
    coarse_north_south_spacing_m = float(np.pi * g3.EARTH_RADIUS_M / grid.coarse_shape[0])
    crossings = _ridge_crossings(grid, faces, elevation_by_reservoir)
    upslope_receivers = _upslope_receiver_mask(grid, faces, elevation_by_reservoir)
    geographic_support = _geographic_support(fields, grid.native_area[::-1], archive["coast_distance_km"],
                                              freeze, crossings)
    source = archive["source_kg_s"]
    cap = archive["cap_kg_m2"] * grid.area
    base_sink = np.full_like(source, RAIN_RATE)
    area = grid.native_area[::-1]
    distance = archive["coast_distance_km"]
    reference = archive["reference_precipitation_mm_yr"]
    maps, results = {}, {}
    moisture_factor = np.ones_like(crossings.rise_m)
    for name, arm in arms.items():
        if isinstance(arm, tuple):
            strength, sink_name, paired_control = arm
        else:
            strength, sink_name, paired_control = arm, name, "control"
        if name in ("upslope_basic", "upslope_delayed"):
            interception = _upslope_interception(
                grid, faces, elevation_by_reservoir, temperature_k_by_reservoir,
                upslope_scale_height_k_per_km,
            )
        else:
            interception = _interception(crossings, faces.first.size, strength * moisture_factor)
        sink = base_sink if sink_profiles is None else sink_profiles[sink_name]
        budget, excess, _ = solve_bounded_budget(grid.area, source, sink, faces, cap, interception=interception)
        if dry_scale is not None and name == "control":
            donor = crossings.donor
            donor_column = budget.mass[donor] / grid.area[donor]
            moisture_factor = np.minimum(donor_column / dry_scale, 1.0)**2
        rain, placement_residual = _place_surface_rain(grid, budget, excess, area, fields, freeze)
        condensate_budget = None
        if name == "upslope_delayed":
            generated_condensate = budget.intercepted_rain_kg_s
            if generated_condensate is None:
                raise AssertionError("upslope arm did not produce a condensate source")
            vapor_rain = budget.rain - generated_condensate
            tolerance = 100 * np.finfo(float).eps * np.maximum(budget.rain, 1.0)
            if np.any(vapor_rain < -tolerance):
                raise AssertionError("intercepted condensate exceeds vapor rain budget")
            vapor_only = SimpleNamespace(rain=np.maximum(vapor_rain, 0.0))
            rain, placement_residual = _place_surface_rain(
                grid, vapor_only, excess, area, fields, freeze,
            )
            condensate_budget = solve_budget(
                grid.area, generated_condensate,
                np.full_like(generated_condensate, 1.0 / condensate_fallout_time_s), faces,
            )
            condensate_surface, condensate_residual = _place_surface_rain(
                grid, condensate_budget, np.zeros_like(generated_condensate),
                area, fields, freeze,
            )
            rain = rain + condensate_surface
            placement_residual = max(placement_residual, condensate_residual)
        maps[name] = rain
        paired_baseline = maps[paired_control] if paired_control in maps else rain
        metrics = _rain_metrics(fields, rain, area, distance, freeze) | _distribution_metrics(fields, rain, area)
        mass_tolerance = 100 * np.finfo(float).eps * cap
        geographic_metrics_present = all((
            not geographic_support[metric]["evaluable"] or value is not None
        ) for metric, value in (
            ("equatorial_belt", metrics["ocean_belt_contrasts"]["0"]["ratio"]),
            ("inland_coast", metrics["land_rain_by_distance_to_ocean_km"]["0_to_500"]["rain_mm_yr"]),
            ("inland_coast", metrics["land_rain_by_distance_to_ocean_km"]["1500_plus"]["rain_mm_yr"]),
            ("ridge_shadow", _paired_rain(crossings, rain, paired_baseline, grid).get(
                "crest_wetter_and_lee_drier_share",
            )),
            ("cold_coast", metrics["cold_coast_rain_mm_yr"]),
        ))
        hard_gates = {
            "finite_nonnegative_mass": bool(np.all(np.isfinite(budget.mass) & (budget.mass >= 0))),
            "within_mass_cap": bool(np.all(budget.mass <= cap + mass_tolerance)),
            "finite_nonnegative_overflow": bool(np.all(np.isfinite(excess) & (excess >= 0))),
            "finite_nonnegative_source": bool(np.all(np.isfinite(source) & (source >= 0))),
            "finite_nonnegative_budget_rain": bool(
                np.all(np.isfinite(budget.rain) & (budget.rain >= 0))
            ),
            "finite_nonnegative_rain": bool(np.all(np.isfinite(rain) & (rain >= 0))),
            "global_closure": budget.global_relative <= 1e-8,
            "summed_cell_closure": budget.summed_cell_relative <= 1e-8,
            "source_equals_budget_rain": bool(
                abs(budget.rain.sum() - source.sum())
                <= 1e-8 * max(float(source.sum()), 1e-300)
            ),
            "placement_closure": placement_residual <= 1e-12,
            "condensate_closure": condensate_budget is None or condensate_budget.global_relative <= 1e-8,
        }
        results[name] = metrics | {
            "all_upslope_receivers": _upslope_rain_metrics(
                grid, rain, maps.get("control", rain), upslope_receivers,
            ) if name in ("upslope_basic", "upslope_delayed") else None,
            "condensate_transport": ({
                "fallout_time_s": condensate_fallout_time_s,
                "mean_land_wind_m_s": mean_land_wind_m_s,
                "estimated_advection_distance_m": mean_land_wind_m_s * condensate_fallout_time_s,
                "coarse_meridional_cell_spacing_m": coarse_north_south_spacing_m,
                "estimated_cells_traveled": (
                    mean_land_wind_m_s * condensate_fallout_time_s / coarse_north_south_spacing_m
                ),
                "generated_kg_s": float(budget.intercepted_rain_kg_s.sum()),
                "fell_out_kg_s": float(condensate_budget.rain.sum()),
                "fallout_closure_relative": condensate_budget.global_relative,
            } if condensate_budget is not None else None),
            "ridge_transects": _paired_rain(crossings, rain, paired_baseline, grid),
            "ridge_native_cell_profiles": _native_crossing_profiles(
                crossings, rain, paired_baseline, grid, fields,
            ),
            "geographic_metrics_present": geographic_metrics_present,
            "cold_inland_budget": _cold_inland_delivery(fields, distance, freeze, grid, faces, source, budget, rain),
            "raw_rms_to_legacy": g3._rms_pair(reference, rain, area)["raw_rel_rms"],
            "summed_cell_closure": budget.summed_cell_relative,
            "global_closure": budget.global_relative,
            "source_total_kg_s": float(source.sum()),
            "budget_rain_total_kg_s": float(budget.rain.sum()),
            "source_rain_relative_residual": float(
                abs(budget.rain.sum() - source.sum()) / max(float(source.sum()), 1e-300)
            ),
            "placed_rain_total_kg_s": float(
                np.sum(rain * area / g3.SECONDS_PER_YEAR)
            ),
            "minimum_column_water_kg_m2": float(np.min(budget.mass / grid.area)),
            "maximum_mass_cap_violation_kg": float(np.maximum(budget.mass - cap, 0.0).max()),
            "overflow_total_kg_s": float(excess.sum()),
            "overflow_max_kg_s": float(excess.max(initial=0.0)),
            "condensate_placement_relative_residual": placement_residual,
            "intercepted_rain_fraction": float(budget.intercepted_rain_kg_s.sum() / budget.rain.sum()),
            "sink_rate_min_s_inv": float(sink.min()),
            "sink_rate_max_s_inv": float(sink.max()),
            "geographic_support": geographic_support,
            "hard_gates": hard_gates,
        }
    return results, maps, area


def _compare_capture(
    capture, map_dir, latitude, freeze, arms, *, solver="finite_volume", dry_scale=None,
    sink_profiles=None, coast_profile_signature=None, coast_length_km=None, coast_alphas=(),
    upslope_scale_height_k_per_km=MOIST_LAPSE_RATE_K_PER_KM,
    condensate_fallout_time_s=CONDENSATE_FALLOUT_TIME_S,
):
    provenance = json.loads((capture / "results.json").read_text())
    save = Path(provenance["save"])
    if provenance.get("heldout"):
        expected = g3.HELDOUT_SHA1.get(save.name) or g3.SPENT_HELDOUT_SHA1.get(save.name)
        if expected != provenance.get("save_sha1"):
            raise ValueError("held-out capture does not match an active or spent registered SHA-1")
    elif TUNING_SAVES.get(save.name) != provenance.get("save_sha256"):
        raise ValueError("expected a pinned tuning capture")
    density = int(provenance.get("density", provenance["grid_shape"][0] // GRID_SHAPE[0]))
    previous_path = capture / "overflow-placement" / "results.json"
    previous = json.loads(previous_path.read_text()) if previous_path.exists() else None
    result = {}
    for wind in ("diagnostic", "cfd"):
        with np.load(capture / f"{wind}-columns.npz", allow_pickle=False) as archive:
            fields = SimpleNamespace(
                is_ocean=archive["is_ocean"], lat_deg=latitude,
                air_temperature_c=archive["air_temperature_c"], elevation_m=archive["elevation_m"],
                wind_u=archive["wind_u_m_s"], wind_v=archive["wind_v_m_s"],
            )
            if coast_profile_signature is not None:
                grid = aggregate_inputs(fields.is_ocean[::-1], {}, factor=density)
                if (not np.array_equal(grid.native_ids, coast_profile_signature["native_ids"])
                        or not np.array_equal(grid.is_ocean, coast_profile_signature["is_ocean"])
                        or not np.array_equal(grid.area, coast_profile_signature["area"])):
                    raise ValueError("diagnostic and CFD captures have different reservoir grids")
                profile, _ = _coast_sink_profiles(
                    grid, archive["coast_distance_km"], coast_length_km, alphas=coast_alphas,
                )
                if not np.allclose(profile, coast_profile_signature["profile"], rtol=0, atol=1e-15):
                    raise ValueError("diagnostic and CFD coast-distance profiles do not align")
            result[wind], maps, area = _solve_arms(
                fields, archive, freeze, arms, solver=solver, dry_scale=dry_scale,
                sink_profiles=sink_profiles,
                upslope_scale_height_k_per_km=upslope_scale_height_k_per_km,
                condensate_fallout_time_s=condensate_fallout_time_s,
            )
        if sink_profiles is None and previous is not None and solver == "finite_volume":
            expected = previous[wind]["broad_transport_half_spread"]["raw_rms_to_legacy"]
            np.testing.assert_allclose(result[wind]["control"]["raw_rms_to_legacy"], expected, rtol=1e-10, atol=0)
        elif sink_profiles is not None:
            for name, values in result[wind].items():
                arm = arms[name]
                strength = arm[0] if isinstance(arm, tuple) else arm
                if strength > 0:
                    paired_control = result[wind][arm[2]] if isinstance(arm, tuple) else result[wind]["control"]
                    values["tuning_screen"] = _tuning_screen(
                        values, result[wind]["ridge_025"], paired_control,
                    )
        if map_dir is not None:
            for name in arms.keys() - {"control"}:
                if name not in ("ridge_025", "upslope_basic", "upslope_delayed") and not name.startswith(("coast_", "uniform_")):
                    continue
                candidate = "ridge_025" if name == "ridge_025" else name
                baseline = maps["control"] if name == "ridge_025" else maps["ridge_025"]
                g3._save_common_scale_rain_map(
                    map_dir, save, density, wind, solver, candidate,
                    baseline, maps[name], area,
                    labels=("Broad transport and spread" if name == "ridge_025" else "Frozen ridge_025",
                            f"With {name}", "Change"),
                )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path)
    parser.add_argument("--maps", type=Path)
    parser.add_argument("--out", type=Path, help="write compact result JSON to this path")
    parser.add_argument("--solver", choices=("finite_volume", "wind_arrows", "both"), default="both")
    parser.add_argument(
        "--include-profile-samples", action="store_true",
        help="include every crossing profile instead of only the compact summaries",
    )
    sensitivity = parser.add_mutually_exclusive_group()
    sensitivity.add_argument("--weak-sensitivity", action="store_const", dest="profile", const="weak")
    sensitivity.add_argument("--mid-sensitivity", action="store_const", dest="profile", const="mid")
    sensitivity.add_argument("--dry-sensitivity", action="store_const", dest="profile", const="dry")
    sensitivity.add_argument("--strength-sensitivity", action="store_const", dest="profile", const="strengths")
    sensitivity.add_argument("--upslope-sensitivity", action="store_const", dest="profile", const="upslope")
    sensitivity.add_argument(
        "--coast-rainout-sensitivity", action="store_const", dest="profile", const="coast",
    )
    parser.add_argument("--coast-length-km", type=float, choices=(500.0, 750.0, 1500.0), default=750.0)
    parser.add_argument("--coast-alpha", type=float, choices=(1.0, 2.0, 3.0))
    parser.add_argument("--tuning-round", type=int, choices=(1, 2, 3))
    parser.add_argument("--moist-lapse-rate-k-per-km", type=float,
                        default=MOIST_LAPSE_RATE_K_PER_KM)
    parser.add_argument("--condensate-fallout-hours", type=float,
                        default=CONDENSATE_FALLOUT_TIME_S / 3600.0)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
    from app.hydroclimate import climate, hydrology

    provenance = json.loads((args.capture / "results.json").read_text())
    heldout = bool(provenance.get("heldout"))
    _validate_capture_identity(provenance, args.profile, parser)
    _validate_upslope_parameters(args, heldout, parser)
    height, width = provenance["grid_shape"]
    latitude, _, _ = climate._build_grid(height, width)
    if args.profile == "coast":
        if args.tuning_round is None:
            parser.error("coast sensitivity requires --tuning-round")
        if args.tuning_round == 1 and (args.coast_alpha is not None or args.coast_length_km != 750.0):
            parser.error("round 1 is fixed to alpha=1,3 at L=750 km")
        if args.tuning_round == 2 and (args.coast_alpha not in (1.0, 3.0)
                                       or args.coast_length_km == 750.0):
            parser.error("round 2 requires one selected alpha and L=500 or 1500 km")
        if args.tuning_round == 3 and args.coast_alpha != 2.0:
            parser.error("round 3 is fixed to alpha=2 at the selected coast length")
    elif args.tuning_round is not None or args.coast_alpha is not None or args.coast_length_km != 750.0:
        parser.error("coast round options require --coast-rainout-sensitivity")
    coast_alphas = (1.0, 3.0) if args.coast_alpha is None else (args.coast_alpha,)
    if args.profile == "coast":
        arms = {
            "control": (0.0, "control", "control"),
            "ridge_025": (0.25, "control", "control"),
        }
        for alpha in coast_alphas:
            suffix = f"{alpha:g}"
            arms[f"control_coast_{suffix}"] = (0.0, f"coast_{suffix}", f"control_coast_{suffix}")
            arms[f"coast_{suffix}"] = (0.25, f"coast_{suffix}", f"control_coast_{suffix}")
            arms[f"control_uniform_{suffix}"] = (0.0, f"uniform_{suffix}", f"control_uniform_{suffix}")
            arms[f"uniform_{suffix}"] = (0.25, f"uniform_{suffix}", f"control_uniform_{suffix}")
    else:
        arms = ({"control": 0.0, "ridge_025": 0.25,
                 "upslope_basic": 1.0, "upslope_delayed": 1.0} if heldout and args.profile == "upslope" else
                {"control": 0.0, "ridge_025": 0.25} if heldout else
                {"weak": WEAK_ARMS, "mid": MID_ARMS, "dry": DRY_ARMS,
                 "strengths": STRENGTH_ARMS,
                 "upslope": {"control": 0.0, "ridge_025": 0.25,
                             "upslope_basic": 1.0, "upslope_delayed": 1.0}}
                .get(args.profile, ARMS))
    dry_scale = DRY_SCALE_KG_M2 if args.profile == "dry" else None
    if args.profile == "coast" and args.solver != "finite_volume":
        parser.error("coast rainout tuning is preregistered for the finite-volume solver only")
    if args.profile == "coast":
        _assert_pair_flow_orientation()
        _assert_coast_profile_orientation()
    solvers = ("finite_volume", "wind_arrows") if args.solver == "both" else (args.solver,)
    results = {}
    sink_parameters = None
    for solver in solvers:
        sink_profiles = None
        coast_profile_signature = None
        if args.profile == "coast":
            first_archive = np.load(args.capture / "diagnostic-columns.npz", allow_pickle=False)
            first_grid = aggregate_inputs(first_archive["is_ocean"][::-1], {}, factor=int(
                first_archive["is_ocean"].shape[0] // GRID_SHAPE[0]
            ))
            profile, sink_profiles = _coast_sink_profiles(
                first_grid, first_archive["coast_distance_km"], args.coast_length_km,
                alphas=coast_alphas,
            )
            np.testing.assert_array_equal(sink_profiles["control"], np.full(first_grid.area.shape, RAIN_RATE))
            sink_parameters = {
                "length_km": args.coast_length_km,
                "alphas": {f"coast_{alpha:g}": alpha for alpha in coast_alphas},
                "uniform_land_multipliers": {
                    f"uniform_{alpha:g}": float(
                        sink_profiles[f"uniform_{alpha:g}"][~first_grid.is_ocean][0] / RAIN_RATE
                    ) for alpha in coast_alphas
                },
                "profile_min": float(profile[~first_grid.is_ocean].min()),
                "profile_max": float(profile[~first_grid.is_ocean].max()),
            }
            coast_profile_signature = {
                "native_ids": first_grid.native_ids.copy(),
                "is_ocean": first_grid.is_ocean.copy(),
                "area": first_grid.area.copy(),
                "profile": profile.copy(),
            }
            first_archive.close()
        results[solver] = _compare_capture(
            args.capture, args.maps, latitude, hydrology.FREEZE_POINT_C, arms,
            solver=solver, dry_scale=dry_scale, sink_profiles=sink_profiles,
            coast_profile_signature=coast_profile_signature,
            coast_length_km=args.coast_length_km, coast_alphas=coast_alphas,
            upslope_scale_height_k_per_km=args.moist_lapse_rate_k_per_km,
            condensate_fallout_time_s=args.condensate_fallout_hours * 3600.0,
        )
    if not args.include_profile_samples:
        for solver_results in results.values():
            for wind_results in solver_results.values():
                for arm_results in wind_results.values():
                    arm_results["ridge_native_cell_profiles"].pop("profiles", None)
    source_provenance = _source_provenance(Path(__file__).resolve().parents[2])
    report = {
        "save": provenance["save"],
        "heldout": heldout,
        "save_digest": provenance.get("save_sha1", provenance.get("save_sha256")),
        "density": provenance.get("density", height // GRID_SHAPE[0]),
        "candidate": (
            "upslope_basic_and_delayed_fallout" if args.profile == "upslope" else
            "ridge_025" if heldout else args.profile or "default"
        ),
        "baseline_reference": (
            "same-run ridge_025 arm; within-run comparisons only" if args.profile == "coast" else None
        ),
        "solvers": list(results),
        "code_commit": source_provenance["commit"],
        "source_provenance": source_provenance,
        "tuning_round": args.tuning_round,
        "parameters": {
            "arm_strength_by_name": {
                name: arm[0] if isinstance(arm, tuple) else arm for name, arm in arms.items()
            },
            "baseline_rainout_s_inv": RAIN_RATE,
            "upslope_rule": ({
                "formula": "1 - exp(-positive_rise_m / Hsat)",
                "Hsat_formula": "Rv * face_temperature_K^2 / (Lv * moist_lapse_rate_K_m)",
                "moist_lapse_rate_k_per_km": args.moist_lapse_rate_k_per_km,
                "temperature_source": "mean captured surface air temperature at face endpoints",
                "moisture_flux": "model transported face flux; no additional humidity multiplier",
                "scope": "directed land-to-land faces; local rain at receiving reservoir",
                "assumptions": "saturated ascent; instant fallout in basic arm; delayed transported fallout in paired arm",
                "delayed_arm": ({
                    "fallout_time_s": args.condensate_fallout_hours * 3600.0,
                    "transport": "same directed face network as vapor",
                } if any(name == "upslope_delayed" for name in arms) else None),
            } if args.profile == "upslope" else None),
            "coast_sink": sink_parameters,
        },
        "capture_provenance": {
            "results_sha256": hashlib.sha256((args.capture / "results.json").read_bytes()).hexdigest(),
            "capture_commit": provenance.get("base_commit"),
            "capture_scripts": provenance.get("script_sha256", provenance.get("scripts")),
        },
        "results": results,
    }
    payload = json.dumps(report, indent=2, allow_nan=False)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload + "\n")
    print(payload)


if __name__ == "__main__":
    main()

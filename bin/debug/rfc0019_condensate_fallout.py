"""Screen #978's six-hour ocean condensate perturbation on registered G3 captures.

Usage:
  backend/.venv/bin/python bin/debug/rfc0019_condensate_fallout.py CAPTURE ... --out PRIVATE.json
  backend/.venv/bin/python bin/debug/rfc0019_condensate_fallout.py \
      --live-save SAVE --g4-runner PINNED_RUNNER.py --out PRIVATE.json
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
from time import perf_counter
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import rfc0019_g3_ridge_rainout as ridge
import rfc0019_g3_compare as g3
from rfc0019_coarse_moisture import aggregate_inputs, distribute_rain, solve_bounded_budget, solve_budget
from rfc0019_g3_compare import SECONDS_PER_YEAR

FALLOUT_SECONDS = 6 * 3600
G4_SHA256 = "09bc7b43eeba05c440734e2ba74fb7b5cd2b6fb14458bbf73e230089e6496527"


def load_g4(path):
    if g3._sha256(path) != G4_SHA256:
        raise ValueError("G4 runner does not match the published #977 checksum")
    spec = importlib.util.spec_from_file_location("rfc0019_g4_coupling", path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    return runner


def install_fallout(runner):
    """Use the pinned runner's face construction, vapor solve and final biome update."""
    candidate = runner._candidate_precipitation

    def update(world, fields, climate, hydrology, biomes):
        ledger = {}
        faces = None

        def make_faces(*args):
            nonlocal faces
            faces = ridge._candidate_faces(*args)
            return faces

        def place(grid, budget, overflow, area, fields, freeze):
            rain, checks = place_condensate_rain(grid, faces, budget, overflow, fields, freeze)
            ledger.update(checks)
            return rain, checks["combined_placement_relative_residual"]

        placement = SimpleNamespace(
            _candidate_faces=make_faces, _place_surface_rain=place,
            _ridge_crossings=ridge._ridge_crossings, _interception=ridge._interception,
        )
        with patch.object(runner, "ridge", placement):
            humidity, rain, biome_ids, diagnostics = candidate(world, fields, climate, hydrology, biomes)
        return humidity, rain, biome_ids, diagnostics | ledger

    runner._candidate_precipitation = update


def time_update(update, world, shape, climate, hydrology, biomes):
    started = perf_counter()
    fields, discarded = g3._compute_candidate_fields(world, *shape, climate, biomes)
    _, rain, _, ledger = update(world, fields, climate, hydrology, biomes)
    return perf_counter() - started - discarded, fields, rain, ledger


def check_live(save, parser, runner, immediate_update, out):
    results = {}
    for density in (1, 4):
        _, digest, world, biomes, climate, hydrology, _, shape = g3._load_tuning_world(
            save, parser, ridge.TUNING_SAVES, density,
        )
        area = g3._grid_area(*shape)[2]
        for wind in ("diagnostic", "cfd"):
            world.wind_model = wind
            current, immediate, variant = [], [], []
            for _ in range(4 if density == 4 else 1):
                started = perf_counter()
                climate.compute_climate(world, *shape)
                current.append(perf_counter() - started)
                if density == 4:
                    elapsed, *_ = time_update(immediate_update, world, shape, climate, hydrology, biomes)
                    immediate.append(elapsed)
                elapsed, fields, rain, ledger = time_update(
                    runner._candidate_precipitation, world, shape, climate, hydrology, biomes,
                )
                variant.append(elapsed)
            total = np.sum(rain * area)
            results[f"{density}-{wind}"] = {
                "global_rain_mm_yr": float(total / area.sum()),
                "land_share": float(np.sum((rain * area)[~fields.is_ocean]) / total),
                "ledger": ledger,
                "runtime_ratio": float(np.median(variant[1:]) / np.median(current[1:])) if density == 4 else None,
                "legacy_seconds": current[1:], "variant_seconds": variant[1:],
                "immediate_seconds": immediate[1:],
            }
    report = {
        "save_sha256": digest, "results": results,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


def place_condensate_rain(grid, faces, vapor, overflow, fields, freeze):
    """Transfer ocean cap overflow to condensate, then deposit at its fallout location."""
    if not np.all(np.isfinite(overflow) & (overflow >= 0)):
        raise AssertionError("invalid vapor cap overflow")
    to_condensate = np.where(grid.is_ocean, overflow, 0.0)
    immediate = overflow - to_condensate
    direct_vapor = replace(vapor, rain=vapor.rain - to_condensate)
    area = grid.native_area[::-1]
    if np.any(direct_vapor.rain):
        direct_rain, placement = ridge._place_surface_rain(grid, direct_vapor, immediate, area, fields, freeze)
    else:
        direct_rain, placement = np.zeros_like(area), 0.0
    condensate = solve_budget(grid.area, to_condensate, np.full_like(overflow, 1 / FALLOUT_SECONDS), faces)
    fallout = distribute_rain(grid, condensate.rain, np.ones_like(grid.native_area))[::-1]
    rain = direct_rain + fallout / area * SECONDS_PER_YEAR
    residual = abs(float(np.sum(rain * area / SECONDS_PER_YEAR) - vapor.rain.sum()))
    residual /= max(float(vapor.rain.sum()), 1.0)
    if residual > 1e-12 or not np.all(np.isfinite(rain) & (rain >= 0)):
        raise AssertionError("combined vapor and condensate precipitation failed placement")
    return rain, {
        "vapor_relative_residual": max(vapor.global_relative, vapor.summed_cell_relative),
        "condensate_relative_residual": max(condensate.global_relative, condensate.summed_cell_relative),
        "combined_placement_relative_residual": max(placement, residual),
        "ocean_cap_source_kg_s": float(to_condensate.sum()),
        "land_fallout_kg_s": float(fallout[~fields.is_ocean].sum()),
    }


def screen_wind(fields, archive, freeze):
    height = fields.is_ocean.shape[0]
    grid = aggregate_inputs(fields.is_ocean[::-1], {}, factor=height // ridge.GRID_SHAPE[0])
    area = grid.native_area[::-1]
    distance = archive["coast_distance_km"]
    _, sinks = ridge._coast_sink_profiles(grid, distance, 1500.0, alphas=(1.0,))
    faces = ridge._candidate_faces(grid, fields.wind_u[::-1], fields.wind_v[::-1], "finite_volume")
    elevation = np.bincount(
        grid.native_ids.ravel(), weights=(fields.elevation_m[::-1] * grid.native_area).ravel(),
        minlength=grid.area.size,
    ) / grid.area
    crossings = ridge._ridge_crossings(grid, faces, elevation)
    support = ridge._geographic_support(fields, area, distance, freeze, crossings)
    source = archive["source_kg_s"]
    cap = archive["cap_kg_m2"] * grid.area

    def solve(sink, strength):
        interception = ridge._interception(crossings, faces.first.size, strength)
        return solve_bounded_budget(grid.area, source, sinks[sink], faces, cap, interception=interception)[:2]

    def immediate_rain(budget, excess):
        return ridge._place_surface_rain(grid, budget, excess, area, fields, freeze)[0]

    ridge_only = immediate_rain(*solve("control", 0.25))
    coast_without_ridge = immediate_rain(*solve("coast_1", 0.0))
    budget, excess = solve("coast_1", 0.25)
    control = immediate_rain(budget, excess)
    variant, ledger = place_condensate_rain(grid, faces, budget, excess, fields, freeze)
    return score_pair(
        fields, area, distance, freeze, grid, crossings, support,
        ridge_only, coast_without_ridge, control, variant,
    ) | {"budget": ledger, "global_source_kg_s": float(source.sum())}


def score_pair(fields, area, distance, freeze, grid, crossings, support,
               ridge_only, coast_without_ridge, control, variant):
    def metrics(rain):
        region = ridge._rain_metrics(fields, rain, area, distance, freeze)
        return region | ridge._distribution_metrics(fields, rain, area) | {
            "geographic_support": support,
            "ridge_transects": ridge._paired_rain(crossings, rain, coast_without_ridge, grid),
            "hard_gates": {},
        }

    reference = metrics(ridge_only)
    baseline = ridge._tuning_screen(metrics(control), reference, reference)
    trial = ridge._tuning_screen(metrics(variant), reference, reference)
    baseline_ratio = baseline["inland_to_coast_ratio"]
    ratio = trial["inland_to_coast_ratio"]
    safeguards = trial["safeguards"] | {
        name: value == 0.0 for name, value in trial["normalized_shortfalls"].items()
        if name != "inland_coast" and value is not None
    }
    if support["inland_coast"]["evaluable"]:
        safeguards["inland_coast_no_worse"] = (
            ratio is not None and baseline_ratio is not None and ratio <= baseline_ratio + 1e-12
        )
    return {"baseline": baseline, "variant": trial, "safeguards": safeguards}


def screen_capture(folder, parser):
    provenance = json.loads((folder / "results.json").read_text())
    ridge._validate_capture_identity(provenance, "coast", parser)
    from app.hydroclimate import climate, hydrology

    lat, _, _ = climate._build_grid(*provenance["grid_shape"])
    results = {}
    for wind in ("diagnostic", "cfd"):
        path = folder / f"{wind}-columns.npz"
        with path.open("rb") as stream:
            sha256 = hashlib.file_digest(stream, "sha256").hexdigest()
        with np.load(path, allow_pickle=False) as archive:
            fields = SimpleNamespace(
                lat_deg=lat, is_ocean=archive["is_ocean"], air_temperature_c=archive["air_temperature_c"],
                elevation_m=archive["elevation_m"], wind_u=archive["wind_u_m_s"], wind_v=archive["wind_v_m_s"],
            )
            results[wind] = screen_wind(fields, archive, hydrology.FREEZE_POINT_C) | {"archive_sha256": sha256}
    return {"save_sha256": provenance["save_sha256"], "winds": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("captures", nargs="*", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--live-save", type=Path)
    parser.add_argument("--g4-runner", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
    if args.live_save:
        if not args.g4_runner or args.captures:
            parser.error("--live-save requires --g4-runner and no capture directories")
        runner = load_g4(args.g4_runner)
        immediate = runner._candidate_precipitation
        install_fallout(runner)
        check_live(args.live_save, parser, runner, immediate, args.out)
        return
    if not args.captures:
        parser.error("supply capture directories or --live-save")
    started = perf_counter()
    results = {capture.name: screen_capture(capture, parser) for capture in args.captures}
    report = {
        "fallout_seconds": FALLOUT_SECONDS, "captures": results,
        "wall_seconds": perf_counter() - started, "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()

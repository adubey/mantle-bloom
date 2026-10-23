#!/usr/bin/env python3
"""RFC 0015: measure how a saved world represents frozen water, without stepping it.

Reads one `.mbworld` save and reports, as JSON:

- **Land ice on the map** (the Biome/Combined views' fine grid, area-weighted by cos(lat)):
  land under visible ice (>= GLACIER_VISIBLE_DEPTH_M), the share of it the Biome view
  classifies as Köppen EF (Ice Cap) or ET (Tundra), the EF land that carries no ice, the
  mean ice thickness on visibly iced land, and the share of iced land whose warmest month is
  above 0 C at the bed but would fall below it at the ice surface (LAPSE_RATE_C_PER_KM per km
  of ice).
- **Snowfall on frozen land** (same grid): hydrology treats all precipitation on land whose
  annual mean is below FREEZE_POINT_C as snow. Under the model's own sinusoidal annual
  cycle (`biomes.seasonal_temp_amplitude`), the precipitation-weighted share of the year
  above 0 C on that land is the share it counts as snow but that falls in thawing months.
- **Ice and sea level** (plate nodes, exact areas): the ice volume `eustasy` debits from the
  ocean (all `glacier_depth`, x 0.917), split into floating ice, grounded marine ice and land
  ice, against the water equivalent of the ice's volume above flotation.

Usage:

    cd backend
    .venv/bin/python ../bin/debug/measure_frozen_water.py <save>.mbworld --out <scratch>/frozen.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


def area_weights(lat_deg: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    return np.broadcast_to(np.cos(np.radians(lat_deg)).reshape(-1, 1), shape)


def share(mask: np.ndarray, within: np.ndarray, weights: np.ndarray) -> float:
    total = float(weights[within].sum())
    return float(weights[within & mask].sum()) / total if total > 0 else 0.0


def map_measures(world) -> dict:
    """Land ice and snowfall on the Biome view's own grid and classification."""
    from app import render_image
    from app.hydroclimate import biomes, hydrology
    from app.hydroclimate.climate import LAPSE_RATE_C_PER_KM

    fields, winner_ids, _share, _runner = render_image._classify_biome_grid(world)
    lat_deg, _lon, _xyz, elevation, is_ocean, air_temp, _ocean_temp, precip, lake_depth, glacier_depth, _is_sea, _hs = fields
    shape = elevation.shape
    w = area_weights(lat_deg, shape)
    biome = winner_ids.reshape(shape)
    land = ~is_ocean & (lake_depth <= hydrology.LAKE_MIN_VISIBLE_DEPTH_M)
    iced = land & (glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M)
    lat_grid = np.broadcast_to(lat_deg.reshape(-1, 1), shape)
    amplitude = biomes.seasonal_temp_amplitude(lat_grid, biomes.grid_continentality(is_ocean), world.axial_tilt_deg)
    warm_bed = air_temp + amplitude
    warm_ice = warm_bed - LAPSE_RATE_C_PER_KM * glacier_depth / 1000.0
    frozen = land & (air_temp < hydrology.FREEZE_POINT_C)
    d = np.clip(air_temp / np.maximum(amplitude, 1e-9), -1.0, 1.0)
    thaw_share = np.arccos(-d) / np.pi  # share of the year above 0 C under the sinusoid
    pw = w * precip
    return {
        "land_area_share_of_globe": float(w[land].sum() / w.sum()),
        "iced_share_of_land": share(iced, land, w),
        "iced_land_classified_EF": share(biome == biomes.ICE_CAP, iced, w),
        "iced_land_classified_ET": share(biome == biomes.koppen_index("ET"), iced, w),
        "EF_land_without_ice": share(~iced, land & (biome == biomes.ICE_CAP), w),
        "EF_share_of_land": share(biome == biomes.ICE_CAP, land, w),
        "iced_land_mean_thickness_m": float(np.average(glacier_depth[iced], weights=w[iced])) if iced.any() else 0.0,
        "iced_land_warm_month_above_0_at_bed_below_0_at_surface": share((warm_bed >= 0.0) & (warm_ice < 0.0), iced, w),
        "frozen_share_of_land": share(frozen, land, w),
        "frozen_land_precip_in_thawing_months": float((pw[frozen] * thaw_share[frozen]).sum() / max(pw[frozen].sum(), 1e-9)),
    }


def sea_level_measures(world) -> dict:
    """The ice eustasy debits from the ocean, by where it sits, against its volume above flotation."""
    from app import plates as plates_mod
    from app.elevation_lines import line_spacing_rad
    from app.hydroclimate import ice_sheet
    from app.lithosphere import RHO_ICE, RHO_WATER

    area = plates_mod.collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))
    depth = np.clip(plates_mod.collect_all_glacier_depth(world.plates), 0.0, None)
    bed = plates_mod.collect_all_elevation(world.plates) - plates_mod.collect_all_ice_load_deflection(world.plates)
    sea = float(world.sea_level_m)
    water_depth = np.clip(sea - bed, 0.0, None)
    grounded = ice_sheet.is_grounded(depth, bed, sea)
    marine = water_depth > 0.0
    above_flotation = np.where(grounded, depth - ice_sheet.flotation_thickness_m(bed, sea), 0.0)
    ratio = RHO_ICE / RHO_WATER

    def km3_we(mask: np.ndarray, column: np.ndarray = depth) -> float:
        return ratio * float(np.dot(np.where(mask, column, 0.0), area)) / 1e9

    return {
        "sea_level_m": sea,
        "ice_debited_km3_we": km3_we(np.ones_like(grounded)),
        "floating_ice_km3_we": km3_we(~grounded & (depth > 0.0)),
        "grounded_marine_ice_km3_we": km3_we(grounded & marine),
        "land_ice_km3_we": km3_we(grounded & ~marine),
        "above_flotation_km3_we": km3_we(grounded, np.clip(above_flotation, 0.0, None)),
        "ocean_area_km2": float(area[marine & ~grounded].sum()) / 1e6,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    sys.path.insert(0, str(args.backend.resolve()))
    from app import persistence  # noqa: E402

    data = args.save.read_bytes()
    world = persistence.load_world_bytes(data)
    result = {
        "save": args.save.name,
        "save_sha256": hashlib.sha256(data).hexdigest(),
        "elapsed_myr": world.elapsed_years / 1e6,
        "map": map_measures(world),
        "sea_level": sea_level_measures(world),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()

"""Shared lake-state snapshot for the RFC 0014 G5 production replay drivers."""

from __future__ import annotations

def lake_ice_snapshot(world, plates, line_spacing_rad):
    """Return aligned lake-state arrays and masks from a world's current hydro cache."""
    hydro = world.hydrology_cache
    plate_order = hydro.plates_in_order if hydro is not None else world.plates
    liquid_m = plates.collect_all_lake_depth(plate_order)
    ice_we_m = plates.collect_all_lake_ice(plate_order)
    elevations = plates.collect_all_elevation(plate_order)
    areas_m2 = plates.collect_all_accounting_areas_m2(
        plate_order, line_spacing_rad(world.node_density)
    )
    if hydro is None or len(hydro.is_ocean) != len(liquid_m):
        is_ocean = elevations <= world.sea_level_m
    else:
        is_ocean = hydro.is_ocean
    total_water_m = liquid_m + ice_we_m
    wet = (total_water_m > 1e-6) & ~is_ocean & (areas_m2 > 0.0)
    shallow = wet & (total_water_m <= 2.0)
    deep = wet & (total_water_m > 2.0)
    bedfast = wet & (liquid_m <= 1e-6) & (ice_we_m > 1e-6)
    return {
        "hydro": hydro,
        "liquid_m": liquid_m,
        "ice_we_m": ice_we_m,
        "total_water_m": total_water_m,
        "areas_m2": areas_m2,
        "wet": wet,
        "shallow": shallow,
        "deep": deep,
        "bedfast": bedfast,
    }

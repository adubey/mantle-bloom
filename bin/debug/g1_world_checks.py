"""World-state checks for RFC 0021's plates-only probes (#980): ledger closure and alternative D's static effect.

`calibrate_plate_forces.py run --world-checks` calls these at each breakdown sample and stores the result in the
sample. Both only read the world, so a run with them steps exactly as one without.

- `ledger_gaps`: each world ledger's closure error over its own scale. The continental-material, craton and
  mobile-cover ledgers are always on. The Hm ledger books only with `World.debug_diagnostics`, so it is
  reported only when the run turned that on (`--hm-ledger`), against the Hm inventory at the start.
  RFC 0021 G1 adopts a correction only if "ledgers close", and #979's review found the harness recorded none of
  them (run.json's `ledger_relative_gap_max` is the rift-line splitter's own interior-thinning closure).
- `alternative_d`: what RFC 0021's alternative D (oceanic Hm itself following age) would do to bathymetry and land
  if switched on with nothing else changed. Every oceanic cell's Hm becomes B.3's force column
  (`torque.force_column_hm`'s law: `max(Hm, min(h1 sqrt(age), H_max) - Hc)`). Its elevation moves by the
  isostatic difference, and sea level is re-solved so the ocean holds the same water. The re-solve is geometric:
  it doesn't separate basins cut off from the ocean, which is enough for a first-order answer. The depth-age
  bins are compared with GDH1 (Stein & Stein 1992).

Usage: imported, not run; see `calibrate_plate_forces.py`.
"""

from __future__ import annotations

import numpy as np

from app import continental_ledger, cratons, hm_ledger, lithosphere, rift_strength, torque
from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from app.hydroclimate import mobile_cover

AGE_BINS_MYR = (0.0, 10.0, 30.0, 60.0, 100.0, 200.0, np.inf)


def hm_inventory_m3(world) -> float | None:
    """The Hm inventory to measure the Hm ledger against, when the world books it."""
    return hm_ledger.inventory_scopes_m3(world)["all"] if getattr(world, "debug_diagnostics", False) else None


def _relative(error: float, scale: float) -> float:
    return abs(float(error)) / max(abs(float(scale)), 1.0)


def ledger_gaps(world, initial_hm_m3: float | None = None) -> dict[str, float | None]:
    """Closure error of each world ledger, relative to its own scale; None where the world doesn't keep it. Reads
    the ledgers without initializing them, so a check never seeds one the run itself hadn't."""
    gaps: dict[str, float | None] = {"continental_material": None, "craton": None, "mobile_cover": None, "hm": None}
    material = getattr(world, "continental_material_ledger", None)
    if material:
        gaps["continental_material"] = _relative(continental_ledger.balance_error_m3(world), material["initial_continental_m3"])
    craton = getattr(world, "craton_ledger", None)
    if craton and "initial_m3" in craton:
        sources = sum(float(craton.get(key, 0.0)) for key in cratons.CRATON_SOURCE_ACCOUNTS)
        gaps["craton"] = _relative(cratons.balance_error_m3(world), sources)
    cover = getattr(world, "mobile_cover_ledger", None)
    if isinstance(cover, dict) and "initial_m3" in cover and all(key in cover for key in mobile_cover.ACCOUNTS):
        sources = sum(cover[key] for key in mobile_cover.SOURCES)
        gaps["mobile_cover"] = _relative(mobile_cover.balance_error_m3(world), sources)
    if initial_hm_m3 is not None and getattr(world, "debug_diagnostics", False):
        accounts = hm_ledger.cumulative_scopes(world)["all"]
        booked = accounts["source_m3"] - accounts["sink_m3"] + accounts["reclassification_m3"]
        gaps["hm"] = _relative(hm_ledger.inventory_scopes_m3(world)["all"] - initial_hm_m3 - booked, initial_hm_m3)
    return gaps


def _columns(world):
    """Per node over every plate: area, elevation, oceanic mask, age (yr), and Hm now and under alternative D."""
    # The world's own cooling law, which B.3 (when on) and the rift-line ocean-floor strength both read.
    law = rift_strength.RiftLineParams.with_overrides(getattr(world, "rift_line_params", None))
    spacing = line_spacing_rad(world.node_density)
    parts = []
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        hc = plate.collect("crustal_thickness_m")
        hm = plate.collect("mantle_lithosphere_thickness_m")
        age = torque.node_age_yr(world, plate)
        oceanic = ~effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        aged = np.where(oceanic, rift_strength.oceanic_strength_hm(hc, hm, age, 0.0, law), hm)
        shift = (lithosphere.isostatic_elevation(hc, aged, lithosphere.RHO_OCEANIC_CRUST)
                 - lithosphere.isostatic_elevation(hc, hm, lithosphere.RHO_OCEANIC_CRUST))
        _, elevation = plate.all_points_and_elevation()
        parts.append((plate.accounting_areas_m2(spacing), elevation, elevation + np.where(oceanic, shift, 0.0), oceanic, age, hm, aged))
    return [np.concatenate(column) for column in zip(*parts)]


def _water_m3(elevation, area, sea_level_m: float) -> float:
    return float(np.dot(area, np.clip(sea_level_m - elevation, 0.0, None)))


def sea_level_for_volume(elevation, area, water_m3: float) -> float:
    """The level at which `water_m3` stands over the surface `elevation`, by bisection (geometric: every cell below
    it counts as ocean)."""
    lo, hi = float(elevation.min()), float(elevation.max()) + water_m3 / float(area.sum()) + 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if _water_m3(elevation, area, mid) < water_m3 else (lo, mid)
    return 0.5 * (lo + hi)


def gdh1_depth_m(age_myr: float) -> float:
    """GDH1 ocean-floor depth below sea level (Stein & Stein 1992): 2600 + 365 sqrt(t) to 20 Myr, then
    5651 - 2473 exp(-0.0278 t)."""
    return float(2600.0 + 365.0 * np.sqrt(age_myr) if age_myr < 20.0 else 5651.0 - 2473.0 * np.exp(-0.0278 * age_myr))


def _depth_bins(area, depth_now, depth_d, age_myr) -> list[dict]:
    bins = []
    for lo, hi in zip(AGE_BINS_MYR[:-1], AGE_BINS_MYR[1:]):
        at = (age_myr >= lo) & (age_myr < hi)
        if not np.any(at):
            continue
        w = area[at]
        mean_age = float(np.average(age_myr[at], weights=w))
        bins.append({"age_myr": [lo, None if np.isinf(hi) else hi], "area_km2": float(w.sum() / 1e6), "mean_age_myr": mean_age,
                     "depth_now_m": float(np.average(depth_now[at], weights=w)), "depth_d_m": float(np.average(depth_d[at], weights=w)),
                     "depth_gdh1_m": gdh1_depth_m(mean_age)})
    return bins


def alternative_d(world) -> dict:
    """Land fraction, sea level, mean ocean depth and the oceanic depth-age curve now and with alternative D."""
    area, z_now, z_d, oceanic, age_yr, hm, aged = _columns(world)
    sea_now = float(world.sea_level_m)
    water = _water_m3(z_now, area, sea_now)
    sea_d = sea_level_for_volume(z_d, area, water)
    total = float(area.sum())
    wet_now, wet_d = z_now < sea_now, z_d < sea_d
    dated = oceanic & np.isfinite(age_yr)
    age_myr = np.where(dated, age_yr, 0.0) / 1e6
    ocean_area = float(area[oceanic].sum())
    return {
        "land_fraction_now": float(area[~wet_now].sum() / total),
        "land_fraction_d": float(area[~wet_d].sum() / total),
        "sea_level_change_m": sea_d - sea_now,
        "mean_ocean_depth_now_m": float(np.average(sea_now - z_now[wet_now], weights=area[wet_now])) if np.any(wet_now) else None,
        "mean_ocean_depth_d_m": float(np.average(sea_d - z_d[wet_d], weights=area[wet_d])) if np.any(wet_d) else None,
        "oceanic_area_share_deepened": float(area[oceanic & (aged > hm)].sum() / ocean_area) if ocean_area else None,
        "oceanic_hm_mean_now_km": float(np.average(hm[oceanic], weights=area[oceanic]) / 1e3) if ocean_area else None,
        "oceanic_hm_mean_d_km": float(np.average(aged[oceanic], weights=area[oceanic]) / 1e3) if ocean_area else None,
        "depth_by_age": _depth_bins(area[dated], (sea_now - z_now)[dated], (sea_d - z_d)[dated], age_myr[dated]) if np.any(dated) else [],
    }


def sample(world, initial_hm_m3: float | None = None) -> dict:
    """What `calibrate_plate_forces.py` stores under a breakdown sample's ``world_checks``."""
    return {"ledgers": ledger_gaps(world, initial_hm_m3), "alternative_d": alternative_d(world)}


def _pool_bins(finals: list[dict]) -> list[dict]:
    rows = {}
    for final in finals:
        for row in final["depth_by_age"]:
            rows.setdefault(tuple(row["age_myr"]), []).append(row)
    pooled = []
    for key, group in sorted(rows.items(), key=lambda item: item[0][0]):
        w = np.array([row["area_km2"] for row in group])
        mean = lambda name: float(np.average([row[name] for row in group], weights=w))  # noqa: E731
        age = mean("mean_age_myr")
        pooled.append({"age_myr": list(key), "area_km2": float(w.sum()), "mean_age_myr": age, "depth_now_m": mean("depth_now_m"),
                       "depth_d_m": mean("depth_d_m"), "depth_gdh1_m": gdh1_depth_m(age)})
    return pooled


def pool(runs: list[dict]) -> dict | None:
    """A setting's checks over its worlds: each ledger's worst gap over every sample, and alternative D at each world's
    last sample (the mature ocean floor), averaged over worlds with the depth-age bins pooled by area."""
    checked = [[s for s in run["samples"] if "world_checks" in s] for run in runs]
    if not any(checked):
        return None
    names = checked[next(i for i, samples in enumerate(checked) if samples)][0]["world_checks"]["ledgers"]
    found = {name: [s["world_checks"]["ledgers"][name] for samples in checked for s in samples] for name in names}
    gaps = {name: max((g for g in values if g is not None), default=None) for name, values in found.items()}
    finals = [samples[-1]["world_checks"]["alternative_d"] for samples in checked if samples]
    keys = ("land_fraction_now", "land_fraction_d", "sea_level_change_m", "mean_ocean_depth_now_m", "mean_ocean_depth_d_m",
            "oceanic_area_share_deepened", "oceanic_hm_mean_now_km", "oceanic_hm_mean_d_km")
    final = {key: float(np.mean([f[key] for f in finals if f[key] is not None])) for key in keys}
    return {"ledger_gap_max": gaps, "alternative_d_final": {"worlds": len(finals), "myr": [samples[-1]["myr"] for samples in checked if samples],
                                                            **final, "depth_by_age": _pool_bins(finals)}}

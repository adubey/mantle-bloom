#!/usr/bin/env python3
"""RFC 0016 G1b's treatment: the G1a land law, installed on a replay by wrapping hydrology.

Not a production path. `install` replaces `_partition_precipitation_sources` and, unless
`humidity` is off, `_vegetation_transpiration_source` for one replay, so no model file changes.
`replay_rfc0016_river_water.py --treatment ...` calls it. See `analysis/rfc0016/README.md`,
"G1a preregistration", for the equations; the step numbers below refer to that list.

The partition's `dry_et` return carries ET plus the store's change, so the land ledger's rain,
which sums the partition, stays the arm's real rain. `pop_step` reports the two apart. The
store remaps between steps by nearest node, scaled to keep its volume, and any volume clipped
at capacity is reported as `store_remap_m3`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numba
import numpy as np
from scipy.spatial import cKDTree

MONTHS = 12
LATENT_HEAT_MJ_MM = 2.45
W_M2_TO_MJ_M2_DAY = 0.0864
DAYS_PER_MONTH = 365.25 / MONTHS
REPEAT_TOLERANCE = 1e-6
MAX_SPIN_YEARS = 500
IDENTITY_LIMIT_MM = 1e-3
# Wei et al. (2017) table 1, T/(E+T+I), by Köppen group. BW is a design choice; EF has no plants.
SPLIT_BY_GROUP = {"A": 0.61, "C": 0.55, "D": 0.55, "BS": 0.65, "ET": 0.65, "BW": 0.25, "EF": 0.0}
SPLIT_OVER_OCEAN_CELL = 0.55


@dataclass(frozen=True)
class Arm:
    capacity_mm: float = 300.0
    depletion_p: float = 0.5
    demand_factor: float = 1.0
    split_shift: float = 0.0
    humidity: bool = True


def transpiration_split_table(biomes, shift: float = 0.0) -> np.ndarray:
    """The potential transpiration share by biome id, shifted and clipped."""
    shares = []
    for code in biomes.KOPPEN_CODES:
        group = next(key for key in (code, code[:2], code[0]) if key in SPLIT_BY_GROUP)
        shares.append(SPLIT_BY_GROUP[group])
    table = np.array(shares + [SPLIT_OVER_OCEAN_CELL] * len(biomes.PELAGIC_NAMES))
    plantless = table == 0.0
    return np.where(plantless, 0.0, np.clip(table + shift, 0.0, 1.0))


def monthly_insolation_w_m2(lat_deg, tilt_deg: float, solar_constant_w_m2: float) -> np.ndarray:
    """Daily-mean top-of-atmosphere flux, (months, nodes), by the formula of
    `climate.compute_insolation_w_m2`. Month 0 starts at the northern spring equinox."""
    lat = np.radians(lat_deg)[None, :]
    longitude = (np.arange(MONTHS) + 0.5) * 2.0 * np.pi / MONTHS
    declination = np.arcsin(np.sin(np.radians(tilt_deg)) * np.sin(longitude))[:, None]
    hour = np.arccos(np.clip(-np.tan(lat) * np.tan(declination), -1.0, 1.0))
    daily = hour * np.sin(lat) * np.sin(declination) + np.cos(lat) * np.cos(declination) * np.sin(hour)
    return np.clip(daily, 0.0, None) * solar_constant_w_m2 / np.pi


def monthly_forcing(lat_deg, temperature_c, amplitude_c, rain_mm, tilt_deg, solar_constant_w_m2, biomes):
    """Monthly rain and Oudin demand in mm, each (months, nodes): G1a steps 1-3."""
    warm_peak = np.where(lat_deg >= 0.0, 3.5, 9.5)
    phase = 2.0 * np.pi * ((np.arange(MONTHS) + 0.5)[:, None] - warm_peak) / MONTHS
    monthly_temperature = temperature_c + amplitude_c * np.cos(phase)
    strength = np.clip(tilt_deg / biomes._SEASON_TILT_REFERENCE_DEG, 0.0, 1.3)
    summer_share, _ = biomes._precip_season(lat_deg, strength)
    # Exactly six months: a test on `cos(phase)` lands on its zeros and can take a seventh (#751).
    warm_half = np.mod((np.arange(MONTHS) + 0.5)[:, None] - (warm_peak - 3.0), MONTHS) < 6.0
    monthly_rain = rain_mm * np.where(warm_half, summer_share, 1.0 - summer_share) / 6.0
    radiation = monthly_insolation_w_m2(lat_deg, tilt_deg, solar_constant_w_m2) * W_M2_TO_MJ_M2_DAY
    demand = radiation / LATENT_HEAT_MJ_MM * np.clip(monthly_temperature + 5.0, 0.0, None) / 100.0 * DAYS_PER_MONTH
    return monthly_rain, demand


@numba.njit(cache=True)
def month_update(store, rain, demand, capacity, critical):
    """One month of `dS/dt = rain - demand * min(1, S/critical)`, capped at `capacity`, in closed
    form: G1a step 5. Returns (store, et, excess)."""
    remaining, et, excess = 1.0, 0.0, 0.0
    while remaining > 0.0:
        if store >= critical:
            net = rain - demand
            if net > 0.0:
                to_full = (capacity - store) / net
                if to_full >= remaining:
                    return store + net * remaining, et + demand * remaining, excess
                remaining -= to_full
                et += demand * remaining + demand * to_full
                return capacity, et, excess + net * remaining
            if net == 0.0 or store + net * remaining >= critical:
                return store + net * remaining, et + demand * remaining, excess
            to_critical = (store - critical) / -net
            et += demand * to_critical
            store = critical
            remaining -= to_critical
            # At the critical store, a falling store enters the limited regime below.
        rate = demand / critical
        if rate == 0.0:
            to_critical = (critical - store) / rain if rain > 0.0 else np.inf
            if to_critical >= remaining:
                return store + rain * remaining, et, excess
            store = critical
            remaining -= to_critical
            continue
        balance = rain / rate
        if balance > critical and store < critical:
            to_critical = np.log((balance - store) / (balance - critical)) / rate
            if to_critical < remaining:
                et += rain * to_critical - (critical - store)
                store = critical
                remaining -= to_critical
                continue
        final = balance + (store - balance) * np.exp(-rate * remaining)
        return final, et + rain * remaining - (final - store), excess
    return store, et, excess


@numba.njit(cache=True)
def year_update(store, rain, demand, capacity, critical):
    et_year, excess_year = 0.0, 0.0
    for month in range(rain.shape[0]):
        store, et, excess = month_update(store, rain[month], demand[month], capacity, critical)
        et_year += et
        excess_year += excess
    return store, et_year, excess_year


@numba.njit(cache=True)
def interval_update(store, rain, demand, capacity, critical, years, tolerance, max_spin):
    """`years` of the bucket from `store`, per node: G1a step 8. Returns the closing store and the
    interval's ET and excess, in mm, and the most years any node took to repeat."""
    nodes = store.shape[0]
    closing, et_total, excess_total = np.zeros(nodes), np.zeros(nodes), np.zeros(nodes)
    longest = 0
    for node in range(nodes):
        current = store[node]
        whole_years = int(years)
        for year in range(1, min(whole_years, max_spin) + 1):
            previous = current
            current, et, excess = year_update(current, rain[:, node], demand[:, node], capacity, critical)
            et_total[node] += et
            excess_total[node] += excess
            if abs(current - previous) <= tolerance * capacity or year == min(whole_years, max_spin):
                # A year within tolerance still drifts by its store change, so the repeated year
                # is closed on rain, or the drift times 1e5 years would be unbooked water.
                drift = rain[:, node].sum() - et - excess
                if excess > 0.0:
                    excess += drift
                else:
                    et += drift
                et_total[node] += et * (years - year)
                excess_total[node] += excess * (years - year)
                longest = max(longest, year)
                break
        closing[node] = current
    return closing, et_total, excess_total, longest


class LandStore:
    """One replay's store, written by the partition wrapper and read by the humidity wrapper."""

    def __init__(self, world, arm: Arm, biomes, solar_constant_w_m2: float) -> None:
        self.world, self.arm, self.biomes, self.solar = world, arm, biomes, solar_constant_w_m2
        self.split_table = transpiration_split_table(biomes, arm.split_shift)
        self.points = self.depth_mm = self.volume_weight = None
        self.call = {}
        self.transpiration = None
        self.record = {}

    def node_split(self, points) -> np.ndarray:
        fields = self.world.climate_cache
        _, cell = cKDTree(fields.world_xyz.reshape(-1, 3)).query(points, k=1)
        return self.split_table[fields.biome_ids.ravel()[cell]]

    def opening(self, points, weight, rain, demand, eligible, capacity, critical) -> np.ndarray:
        """This step's opening store. The first opens on the repeating year (G1a step 9)."""
        if self.points is None:
            spun, _, _, _ = interval_update(
                np.full(len(points), 0.5 * capacity), rain, demand, capacity, critical,
                float(MAX_SPIN_YEARS), REPEAT_TOLERANCE, MAX_SPIN_YEARS,
            )
            self.record["store_remap_m3"] = 0.0
            return np.where(eligible, spun, 0.0)
        _, nearest = cKDTree(self.points).query(points, k=1)
        remapped = np.where(eligible, self.depth_mm[nearest], 0.0)
        before = float(np.dot(self.depth_mm, self.volume_weight))
        after = float(np.dot(remapped, weight))
        scaled = remapped * (before / after) if after > 0.0 else remapped
        opening = np.minimum(scaled, capacity)
        self.record["store_remap_m3"] = (before - float(np.dot(opening, weight))) / 1000.0
        return opening

    def partition(self, liquid_mm, areas_m2, wet_fraction, budyko_runoff_m3):
        call = self.call
        points, temperature, amplitude, years = call["points"], call["temperature"], call["amplitude"], call["years"]
        arm = self.arm
        liquid = np.clip(np.asarray(liquid_mm, dtype=float), 0.0, None)
        dry = 1.0 - np.clip(np.asarray(wet_fraction, dtype=float), 0.0, 1.0)
        lat_deg = np.degrees(np.arcsin(np.clip(points[:, 2] / np.linalg.norm(points, axis=1), -1.0, 1.0)))
        rain, demand = monthly_forcing(
            lat_deg, temperature, amplitude, liquid, self.world.axial_tilt_deg, self.solar, self.biomes
        )
        demand = demand * arm.demand_factor
        eligible = (liquid > 0.0) & (dry > 0.0)
        capacity, critical = arm.capacity_mm, (1.0 - arm.depletion_p) * arm.capacity_mm
        weight = areas_m2 * dry
        opening = self.opening(points, weight, rain, demand, eligible, capacity, critical)
        closing, et_mm, excess_mm, longest = interval_update(
            opening, rain, demand, capacity, critical, float(years), REPEAT_TOLERANCE, MAX_SPIN_YEARS
        )
        et_mm, excess_mm = np.where(eligible, et_mm, 0.0), np.where(eligible, excess_mm, 0.0)
        closing = np.where(eligible, closing, 0.0)
        identity = np.where(eligible, liquid * years - et_mm - excess_mm - (closing - opening), 0.0)
        # The land ledger sums this partition, so it can't see rain lost here; only this can.
        if np.abs(identity).max() > IDENTITY_LIMIT_MM:
            raise RuntimeError(f"store water balance off by {np.abs(identity).max():.3g} mm")
        split = self.node_split(points)
        to_m3_per_yr = areas_m2 * dry / 1000.0 / years
        runoff = excess_mm * to_m3_per_yr
        et = et_mm * to_m3_per_yr
        change = (closing - opening) * to_m3_per_yr
        self.points, self.depth_mm, self.volume_weight = points, closing, weight
        self.transpiration = (points, split * et_mm / years * dry, split * et * years, areas_m2, eligible)
        total_budyko = float(budyko_runoff_m3.sum())
        self.record.update(
            excess_ratio=float(runoff.sum()) / total_budyko,
            et_m3=float(et.sum()) * years,
            transpiration_m3=float(np.dot(split, et)) * years,
            store_m3=float(np.dot(closing, weight)) / 1000.0,
            store_change_m3=float(change.sum()) * years,
            max_identity_mm=float(np.abs(identity).max()),
            max_repeat_years=int(longest),
        )
        return runoff, et + change

    def humidity_source(self, control_source, is_ocean, precip_mm_per_g_kg):
        """Last step's transpiration in the climate's units, on cells with store nodes (step 10)."""
        fields = self.world.climate_cache
        if self.transpiration is None or fields is None or fields.is_ocean.shape != is_ocean.shape:
            return control_source
        points, transpiration_mm, _, areas, eligible = self.transpiration
        _, cell = cKDTree(fields.world_xyz.reshape(-1, 3)).query(points, k=1)
        size = is_ocean.size
        volume = np.bincount(cell, weights=transpiration_mm * areas, minlength=size)
        area = np.bincount(cell, weights=areas, minlength=size)
        held = np.bincount(cell, weights=eligible.astype(float), minlength=size) > 0.0
        depth = np.divide(volume, area, out=np.zeros(size), where=area > 0.0).reshape(is_ocean.shape)
        source = np.where(held.reshape(is_ocean.shape) & ~is_ocean, depth / precip_mm_per_g_kg, control_source)
        land = ~is_ocean
        self.record["humidity_source_ratio"] = float(source[land].sum() / max(control_source[land].sum(), 1e-12))
        return source

    def pop_step(self) -> dict:
        record, self.record = self.record, {}
        return record


def install(world, arm: Arm, hydrology, climate, biomes, plates_mod) -> LandStore:
    store = LandStore(world, arm, biomes, climate.SOLAR_CONSTANT_W_M2)
    control_hydrology = hydrology.compute_hydrology
    control_partition = hydrology._partition_precipitation_sources
    control_source = climate._vegetation_transpiration_source

    def compute(world_, precipitation, temperature, years, *rest, **kw):
        cloud = kw.get("node_cloud") or (rest[0] if rest else None)
        points = cloud[0] if cloud is not None else plates_mod.gather_node_positions(world_.plates)[0]
        store.call = {
            "points": np.asarray(points), "temperature": np.asarray(temperature),
            "amplitude": kw["seasonal_amplitude_at_nodes"], "years": years,
        }
        return control_hydrology(world_, precipitation, temperature, years, *rest, **kw)

    def partition(liquid, evaporation, areas, wet_fraction):
        wet, runoff, et = control_partition(liquid, evaporation, areas, wet_fraction)
        if store.call.get("years", 0.0) <= 0.0:
            return wet, runoff, et
        new_runoff, new_et = store.partition(liquid, areas, wet_fraction, runoff)
        return wet, new_runoff, new_et

    def source(world_, elevation, is_ocean):
        return store.humidity_source(control_source(world_, elevation, is_ocean), is_ocean, climate.PRECIP_MM_PER_G_KG)

    hydrology.compute_hydrology = compute
    hydrology._partition_precipitation_sources = partition
    if arm.humidity:
        climate._vegetation_transpiration_source = source
    return store


class RunoffScale:
    """The partition-match and sham arms: control hydrology with runoff scaled for each step."""

    def __init__(self, factors) -> None:
        self.factors = list(factors)
        self.record = {}

    def next_factor(self) -> float:
        return self.factors[min(self.record.get("calls", 0), len(self.factors) - 1)]

    def pop_step(self) -> dict:
        calls = self.record.get("calls", 0)
        record, self.record = dict(self.record), {"calls": calls}
        return record


def install_runoff_scale(factors, hydrology) -> RunoffScale:
    scale = RunoffScale(factors)
    control_partition = hydrology._partition_precipitation_sources

    def partition(liquid, evaporation, areas, wet_fraction):
        wet, runoff, et = control_partition(liquid, evaporation, areas, wet_fraction)
        factor = scale.next_factor()
        available = runoff + et
        scaled = np.minimum(runoff * factor, available)
        scale.record.update(calls=scale.record.get("calls", 0) + 1, runoff_factor=factor)
        return wet, scaled, available - scaled

    hydrology._partition_precipitation_sources = partition
    return scale

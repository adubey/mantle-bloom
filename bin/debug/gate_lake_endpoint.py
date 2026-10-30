#!/usr/bin/env python3
"""RFC 0014 G1 refresh for the endpoint candidate and C fallback (#589).

Load one preregistered diagnosis save, run one 10 kyr world step, and capture the current
production C forcing and D2(c) phase inputs. On isolated leaf basins that neither merge nor
spill, compare equilibrium endpoint volume with an independently substepped C integration.
The default is 1,000 substeps; use `--transient-substeps 4000` for a convergence check. Merged,
incoming-transfer, capped-overflow, invalid, and unsampled bodies are counted as C fallback.
This is an offline qualification screen; its per-body transient measurements are not available
to production code and must not be treated as a production eligibility cache.

Run one save at a time, never alongside another replay, probe, or test:

    backend/.venv/bin/python bin/debug/gate_lake_endpoint.py \\
        ~/Downloads/mantle-bloom-seed910211954-220200000y.mbworld \\
        --out analysis/issue589/seed910211954.json

Repeat for seeds 635056485 and 317739385. The JSON records the save checksum, backend commit,
measurement coverage, volume-weighted pass/fallback shares, and body-level measurements. Keep
the saves and any full transient traces out of the repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "backend"))
from app.hydroclimate import hydrology, lakes

DIAGNOSIS_SAVES = {
    "bd4874092efc805e8b174e115fe567c82d35fc2ee606686938b5b6eb96b060d2": 910211954,
    "ec27d5077bdf43b3a6f0339b5726fa4d110da9eb1f143f96a3004d79f309b045": 635056485,
    "6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd": 317739385,
}
MIN_DIRECT_AREA_KM2 = 10_000.0
TRANSIENT_SUBSTEPS = 1_000
ENDPOINT_TOLERANCE = 0.01


def _args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("save", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--small-sample-fraction", type=float, default=0.05)
    parser.add_argument("--sample-seed", type=int, default=589)
    parser.add_argument("--transient-substeps", type=int, default=TRANSIENT_SUBSTEPS)
    return parser.parse_args()


def _endpoint_level(body, climate, hypsometry, cap_level: float, dry_level: float, supply: float) -> float:
    low, high = dry_level, cap_level
    if supply <= 0.0:
        return low
    if lakes.lake_evaporation_m3_per_yr(high, body.members, hypsometry, climate) <= supply:
        return high
    for _ in range(64):
        mid = 0.5 * (low + high)
        if lakes.lake_evaporation_m3_per_yr(mid, body.members, hypsometry, climate) < supply:
            low = mid
        else:
            high = mid
        if high - low <= 1e-6:
            break
    return 0.5 * (low + high)


def _time_constants(body, climate, hypsometry, level: float) -> tuple[float, float]:
    members = body.members
    span = max(float(hypsometry.elevation[members].max() - hypsometry.elevation[members].min()), 1.0)
    delta = max(1e-3, span * 1e-6)
    lo = max(hypsometry.dry_level(members), level - delta)
    hi = level + delta
    dvolume_dlevel = (hypsometry.volume(hi, members) - hypsometry.volume(lo, members)) / max(hi - lo, 1e-12)
    dloss_dlevel = (
        lakes.lake_evaporation_m3_per_yr(hi, members, hypsometry, climate)
        - lakes.lake_evaporation_m3_per_yr(lo, members, hypsometry, climate)
    ) / max(hi - lo, 1e-12)
    tau = dvolume_dlevel / dloss_dlevel if dloss_dlevel > 0.0 else float("inf")
    return tau, dvolume_dlevel


def _integrate_body_c(body, stored_m3, supply_m3_yr, years, climate, hypsometry, cap_level, dry_level, substeps):
    dt = years / substeps
    volume = stored_m3
    level = dry_level
    for _ in range(substeps):
        volume, _, overflow, level = lakes._integrate_body_volume(
            volume + supply_m3_yr * dt,
            body.members,
            cap_level,
            dry_level,
            dt,
            climate.evaporation_m_per_yr,
            climate.sublimation_m_per_yr,
            hypsometry,
        )
        if overflow > max(1e-6, 1e-12 * max(stored_m3, supply_m3_yr * years, 1.0)):
            return None, "transient_spill"
    return hypsometry.mean_depth(level, body.members), None


def _phase_state(members, liquid_depth, phase_inputs):
    if phase_inputs is None:
        return None
    ice, liquid, _ = hydrology.advance_lake_ice(
        phase_inputs["prior_ice"][members], liquid_depth, phase_inputs["area"][members],
        phase_inputs["mean_c"][members], phase_inputs["amplitude_c"][members],
        phase_inputs["ice_cover"][members], phase_inputs["years"],
    )
    return ice, liquid


def _routing_fallback_reason(body, leaf_of_node, incoming_leaf_ids) -> str | None:
    if body.children:
        return "merged_body"
    members_leafs = {int(item) for item in leaf_of_node[body.members] if item >= 0}
    if len(members_leafs) != 1:
        return "multi_leaf_body"
    if next(iter(members_leafs)) in incoming_leaf_ids:
        return "incoming_transfer_route"
    return None


def _valid_forcing_and_storage(climate, members, stored: float) -> bool:
    source_rates = np.asarray(climate.external_supply_m3_per_yr[members], dtype=float)
    evaporation_rates = np.asarray(climate.evaporation_m_per_yr[members], dtype=float)
    if (
        not np.all(np.isfinite(source_rates))
        or np.any(source_rates < 0.0)
        or not np.all(np.isfinite(evaporation_rates))
        or np.any(evaporation_rates < 0.0)
        or not np.isfinite(stored)
        or stored < 0.0
    ):
        return False
    if climate.sublimation_m_per_yr is not None:
        sublimation_rates = np.asarray(climate.sublimation_m_per_yr[members], dtype=float)
        if not np.all(np.isfinite(sublimation_rates)) or np.any(sublimation_rates < 0.0):
            return False
    return True


def _candidate_records(capture, sample_fraction: float, sample_seed: int, substeps: int) -> list[dict]:
    forest, elevation, previous_depth, years, climate, hypsometry, is_ocean = capture["solve"]
    phase_inputs = capture.get("phase")
    node_count = len(elevation)
    node_area_km2 = 4.0 * np.pi * lakes.PLANET_RADIUS_KM**2 / max(node_count, 1)
    active = lakes.governing_bodies(forest, elevation, previous_depth, hypsometry)
    all_leaves = [lake for lake in lakes.iter_all_lakes(forest) if not lake.children]
    leaf_of_node = np.full(node_count, -1, dtype=np.int64)
    for leaf in all_leaves:
        leaf_of_node[leaf.members] = leaf.lake_id
    incoming_leaf_ids = set()
    for body in active:
        target = body.outlet_target_idx
        if 0 <= target < node_count:
            target_leaf = int(leaf_of_node[target])
            if target_leaf >= 0:
                incoming_leaf_ids.add(target_leaf)
    areas = [float(hypsometry.area_m2[body.members].sum()) for body in active]
    area_km2 = [value / 1e6 for value in areas]
    direct = [index for index, area in enumerate(area_km2) if area >= MIN_DIRECT_AREA_KM2]
    small = [index for index, area in enumerate(area_km2) if area < MIN_DIRECT_AREA_KM2]
    rng = np.random.default_rng(sample_seed)
    take = min(len(small), int(np.ceil(len(small) * np.clip(sample_fraction, 0.0, 1.0))))
    sampled = set(direct) | set(rng.choice(small, size=take, replace=False).tolist())

    records = []
    for index, body in enumerate(active):
        members = body.members
        stored = float(np.dot(previous_depth[members], hypsometry.area_m2[members]))
        record = {
            "lake_id": int(body.lake_id),
            "area_km2": area_km2[index],
            "stored_m3": stored,
            "sampled": index in sampled,
            "eligible": False,
            "fallback_reason": None,
        }
        if index not in sampled:
            record["fallback_reason"] = "small_body_not_sampled"
            records.append(record)
            continue
        routing_reason = _routing_fallback_reason(body, leaf_of_node, incoming_leaf_ids)
        if routing_reason:
            record["fallback_reason"] = routing_reason
            records.append(record)
            continue

        if not _valid_forcing_and_storage(climate, members, stored):
            record["fallback_reason"] = "invalid_forcing_or_storage"
            records.append(record)
            continue
        source_rates = np.asarray(climate.external_supply_m3_per_yr[members], dtype=float)
        supply = float(source_rates.sum())
        _, tier_depth = lakes._classify_tier(members, previous_depth, node_area_km2)
        rim = body.max_depth if body.max_depth is not None else float("inf")
        cap_level = min(rim, body.floor_elevation + tier_depth)
        dry_level = hypsometry.dry_level(members)
        if body.drains_at_floor:
            cap_level = dry_level
        end_evap = lakes.lake_evaporation_m3_per_yr(cap_level, members, hypsometry, climate)
        dry_evap = lakes.lake_evaporation_m3_per_yr(dry_level, members, hypsometry, climate)
        if not np.isfinite(end_evap) or not np.isfinite(dry_evap) or end_evap + 1e-9 < dry_evap:
            record["fallback_reason"] = "non_monotone_or_invalid_loss"
            records.append(record)
            continue
        if supply > end_evap + 1e-9 * max(supply, end_evap, 1.0):
            record["fallback_reason"] = "endpoint_spill"
            records.append(record)
            continue

        endpoint_start = time.perf_counter()
        endpoint_level = _endpoint_level(body, climate, hypsometry, cap_level, dry_level, supply)
        endpoint_depth = hypsometry.mean_depth(endpoint_level, members)
        endpoint_volume = float(np.dot(endpoint_depth, hypsometry.area_m2[members]))
        tau_relax, dvolume_dlevel = _time_constants(body, climate, hypsometry, endpoint_level)
        tau_fill = endpoint_volume / supply if supply > 0.0 else (0.0 if endpoint_volume == 0.0 else float("inf"))
        tau_limit = 0.1 * years
        time_pass = max(tau_relax, tau_fill) <= tau_limit
        record.update(
            supply_m3_yr=supply,
            endpoint_level_m=endpoint_level,
            endpoint_volume_m3=endpoint_volume,
            tau_relax_years=tau_relax,
            tau_fill_years=tau_fill,
            dvolume_dlevel_m2=dvolume_dlevel,
            time_bound_pass=bool(time_pass),
        )
        if not time_pass:
            record["fallback_reason"] = "time_bound"
            records.append(record)
            continue

        previous_wet = members[previous_depth[members] > 0.0]
        old_level = (
            float(np.max(hypsometry.level_from_depth(previous_depth[previous_wet], previous_wet)))
            if len(previous_wet) else dry_level
        )
        old_evap = lakes.lake_evaporation_m3_per_yr(old_level, members, hypsometry, climate)
        inferred_evap = supply - (endpoint_volume - stored) / years
        if (
            not np.isfinite(inferred_evap)
            or inferred_evap < -1e-9 * max(supply, 1.0)
            or inferred_evap > max(old_evap, lakes.lake_evaporation_m3_per_yr(endpoint_level, members, hypsometry, climate))
            + 1e-9 * max(supply, old_evap, 1.0)
        ):
            record["fallback_reason"] = "endpoint_budget_bound"
            records.append(record)
            continue
        record["inferred_step_evaporation_m3_yr"] = max(inferred_evap, 0.0)
        endpoint_eval_runtime_s = time.perf_counter() - endpoint_start

        direct_depth, integration_reason = _integrate_body_c(
            body, stored, supply, years, climate, hypsometry, cap_level, dry_level, substeps
        )
        if integration_reason:
            record["fallback_reason"] = integration_reason
            records.append(record)
            continue
        c_slice_start = time.perf_counter()
        lakes._integrate_body_volume(
            stored + supply * years, members, cap_level, dry_level, years,
            climate.evaporation_m_per_yr, climate.sublimation_m_per_yr, hypsometry,
        )
        c_slice_runtime_s = time.perf_counter() - c_slice_start
        c_volume = float(np.dot(direct_depth, hypsometry.area_m2[members]))
        liquid_diff = abs(c_volume - endpoint_volume) / max(c_volume, endpoint_volume, 1.0)
        record.update(
            c_transient_volume_m3=c_volume,
            liquid_transient_difference_fraction=liquid_diff,
            endpoint_eval_runtime_s=endpoint_eval_runtime_s,
            c_slice_runtime_s=c_slice_runtime_s,
        )
        if phase_inputs is not None:
            endpoint_phase = _phase_state(members, endpoint_depth, phase_inputs)
            c_phase = _phase_state(members, direct_depth, phase_inputs)
            if endpoint_phase is None or c_phase is None:
                record["fallback_reason"] = "missing_phase_measurement"
                records.append(record)
                continue
            endpoint_ice, endpoint_liquid = endpoint_phase
            c_ice, c_liquid = c_phase
            endpoint_total = float(np.dot(endpoint_ice + endpoint_liquid, hypsometry.area_m2[members]))
            c_total = float(np.dot(c_ice + c_liquid, hypsometry.area_m2[members]))
            record["endpoint_combined_volume_m3"] = endpoint_total
            record["c_transient_combined_volume_m3"] = c_total
            record["combined_transient_difference_fraction"] = abs(c_total - endpoint_total) / max(c_total, endpoint_total, 1.0)

        transient_pass = liquid_diff <= ENDPOINT_TOLERANCE
        record["transient_bound_pass"] = bool(transient_pass)
        record["eligible"] = bool(transient_pass)
        if not record["eligible"]:
            reasons = []
            if not transient_pass:
                reasons.append("transient_bound")
            record["fallback_reason"] = "+".join(reasons)
        records.append(record)
    return records


def main() -> None:
    args = _args()
    raw = args.save.read_bytes()
    checksum = hashlib.sha256(raw).hexdigest()
    if checksum not in DIAGNOSIS_SAVES:
        raise SystemExit(f"save checksum is not registered for RFC 0014 G1: {checksum}")

    from app import persistence, world as world_mod
    world = persistence.load_world_bytes(raw)
    world.ice_age_period_years = 0.0
    capture: dict = {}
    original_solve = lakes.solve_lakes_prognostic
    original_ice = hydrology.advance_lake_ice

    def capture_solve(forest, elevation, previous_depth, years, climate, hypsometry, **kwargs):
        frozen_elevation = np.asarray(elevation).copy()
        frozen_hypsometry = lakes.CellHypsometry(
            frozen_elevation,
            np.asarray(hypsometry.area_m2).copy(),
            np.asarray(hypsometry.relief_m).copy(),
            hypsometry.reference_area_m2,
        )
        frozen_climate = lakes.LakeClimate(**{
            name: np.asarray(value).copy() if isinstance(value, np.ndarray) else value
            for name, value in vars(climate).items()
        })
        capture["solve"] = (
            forest, frozen_elevation, np.asarray(previous_depth).copy(),
            years, frozen_climate, frozen_hypsometry,
            None if kwargs.get("is_ocean") is None else np.asarray(kwargs["is_ocean"]).copy(),
        )
        return original_solve(forest, elevation, previous_depth, years, climate, hypsometry, **kwargs)

    def capture_ice(prior_ice, liquid_depth, area, mean_c, amplitude_c, ice_cover_fraction, years):
        capture["phase"] = {
            "prior_ice": np.asarray(prior_ice).copy(),
            "area": np.asarray(area).copy(),
            "mean_c": np.asarray(mean_c).copy(),
            "amplitude_c": np.asarray(amplitude_c).copy(),
            "ice_cover": np.asarray(ice_cover_fraction).copy(),
            "years": years,
        }
        return original_ice(prior_ice, liquid_depth, area, mean_c, amplitude_c, ice_cover_fraction, years)

    lakes.solve_lakes_prognostic = capture_solve
    hydrology.advance_lake_ice = capture_ice
    start = time.perf_counter()
    world_mod.step_world(world, 10_000.0)
    step_runtime = time.perf_counter() - start
    lakes.solve_lakes_prognostic = original_solve
    hydrology.advance_lake_ice = original_ice
    if "solve" not in capture:
        raise RuntimeError("did not capture the production lake solve")

    if args.transient_substeps < 1:
        raise SystemExit("--transient-substeps must be at least one")
    records = _candidate_records(capture, args.small_sample_fraction, args.sample_seed, args.transient_substeps)
    measured = [row for row in records if row["sampled"]]
    total_volume = sum(row["stored_m3"] for row in records)
    measured_volume = sum(row["stored_m3"] for row in measured)
    eligible_volume = sum(row["stored_m3"] for row in measured if row["eligible"])
    time_pass_volume = sum(row["stored_m3"] for row in measured if row.get("time_bound_pass"))
    transient_pass_volume = sum(row["stored_m3"] for row in measured if row.get("transient_bound_pass"))
    eligible_rows = [row for row in measured if row["eligible"]]
    reason_counts: dict[str, int] = {}
    for row in records:
        reason = row["fallback_reason"] or "eligible"
        reason_counts[reason] = reason_counts.get(reason, 0) + 1

    result = {
        "issue": 589,
        "rfc": "0014",
        "seed": DIAGNOSIS_SAVES[checksum],
        "save_sha256": checksum,
        "backend_commit": __import__("subprocess").run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout.strip(),
        "ice_age_control": "never",
        "step_years": 10_000.0,
        "step_runtime_s": step_runtime,
        "transient_substeps": args.transient_substeps,
        "small_body_sample_fraction": args.small_sample_fraction,
        "sample_seed": args.sample_seed,
        "active_body_count": len(records),
        "measured_body_count": len(measured),
        "total_starting_liquid_volume_m3": total_volume,
        "measured_starting_liquid_volume_m3": measured_volume,
        "measured_volume_coverage_fraction": measured_volume / max(total_volume, 1.0),
        "eligible_measured_volume_share": eligible_volume / max(measured_volume, 1.0),
        "time_bound_pass_measured_volume_share": time_pass_volume / max(measured_volume, 1.0),
        "transient_bound_pass_measured_volume_share": transient_pass_volume / max(measured_volume, 1.0),
        "eligible_endpoint_eval_runtime_s": sum(row["endpoint_eval_runtime_s"] for row in eligible_rows),
        "eligible_c_slice_runtime_s": sum(row["c_slice_runtime_s"] for row in eligible_rows),
        "fallback_reason_counts": reason_counts,
        "qualification_note": "Offline direct C integration is not a production eligibility measurement; unmeasured or unsupported bodies use C.",
        "bodies": records,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps({key: value for key, value in result.items() if key != "bodies"}, indent=1))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

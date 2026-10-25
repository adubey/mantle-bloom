"""Mantle plume head footprints (RFC 0007; issue #577 scale-check slice).

Pure geometry only: no event stream, no `World` integration, no tracers or ledger.
The head footprint derives the emplacement area from the head volume through one mean
excess thickness, so the peak excess is exactly `V / A` and typical heads reproduce
the 2-3 km relief plateaus stand above the abyssal floor
[[Kerr & Mahoney, 2007]](#ref-kerr2007). Full design and knob table live in
`rfc/0007-mantle-plumes-and-plateaus.md`; the event stream, tails, weakening,
docking, and ledger arrive with the implementation issues tracked in #573.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass

import numpy as np

KM3_TO_M3 = 1e9
KM2_TO_M2 = 1e6

# Mean excess crust thickness tying head volume to footprint area: with the model's
# water-loaded oceanic isostasy each km of excess raises the surface by ~155 m, so
# 15 km gives ~2.3 km of relief, inside the observed 2-3 km plateau range (#577).
PLATEAU_MEAN_EXCESS_M = 15_000.0
# LIPs cover more than 0.1 M km²; plateaus reach 2 M km² (#577).
PLATEAU_AREA_MIN_KM2 = 1e5
PLATEAU_AREA_MAX_KM2 = 2e6
# Calibration assumption: the range spans typical ≥1 M km³ LIPs to Ontong Java's
# 44-57 M km³; the distribution shape is not established (#577).
PLUME_HEAD_VOLUME_MIN_KM3 = 1e6
PLUME_HEAD_VOLUME_MAX_KM3 = 5e7


def head_footprint_area_m2(volume_m3: float) -> float:
    """Oceanic head footprint area: volume tied to area by one mean excess thickness.

    Very small heads would otherwise vanish below the grid; very large ones are
    capped at the largest observed plateaus, accepting the overshoot the RFC records
    for Ontong-Java-scale heads (issue #577).
    """
    area = volume_m3 / PLATEAU_MEAN_EXCESS_M
    return float(
        np.clip(
            area,
            PLATEAU_AREA_MIN_KM2 * KM2_TO_M2,
            PLATEAU_AREA_MAX_KM2 * KM2_TO_M2,
        )
    )


def head_kernel_radius_m(area_m2: float) -> float:
    """Flat-topped kernel radius whose own integral over the plane equals the area.

    With `w = exp(-(d/r)^4)`, the plane integral is `pi * r^2 * Gamma(3/2)`, so this
    `r` makes the kernel place exactly the head volume at peak `V / A` (#577).
    """
    return math.sqrt(area_m2 / (math.pi * math.gamma(1.5)))


def head_excess_thickness_m(distances_m: np.ndarray, volume_m3: float) -> np.ndarray:
    """Excess crust thickness over the footprint: peak `V / A` with flat-topped falloff.

    The peak is exact by construction; the integral over cells equals the head volume
    up to the grid discretization, which the emplacement pass normalizes away (#577).
    """
    distances = np.asarray(distances_m, dtype=float)
    area = head_footprint_area_m2(volume_m3)
    radius = head_kernel_radius_m(area)
    peak = volume_m3 / area
    return peak * np.exp(-((distances / radius) ** 4))


# ---- Event stream, emplacement, and heat (issue #579) ----
#
# Pure-index event scheduling (every draw keyed by event index, never by step) plus
# the per-step emplacement pass. No tracers, ledger sub-accounts, docking,
# weakening, or views yet -- those arrive with their own issues (#580-#583).

PLUME_HEAD_RATE_PER_MYR = 0.1
PLUME_TAIL_ONLY_RATE_PER_MYR = 0.3
PLUME_HEAD_DURATION_MYR = 3.0
PLUME_TAIL_LIFETIME_MYR = 80.0
PLUME_TAIL_FLUX_MIN_KM3_PER_YR = 0.01
PLUME_TAIL_FLUX_MAX_KM3_PER_YR = 0.25
PLUME_TAIL_RADIUS_KM = 150.0
PLUME_HEAD_RADIUS_KM = 1_100.0
PLUME_UNIFORM_SHARE = 0.3
PLUME_RATE_CORRELATION_MYR = 100.0
PLUME_RATE_MODULATION_MAX = 3.0
PLUME_HEAT_MYR = 20.0
# Lithospheric weakening (#582): while `plume_heat_myr` covers emplaced cells,
# craton-rooted rift resistance is scaled by (1 - PLUME_CRATON_WEAKENING), with
# the product clamped like `cratons.rift_resistance`. Default 0.75. The weakened
# state decays back automatically as the heat countdown expires; decay is gradual
# (interval-overlap scaling), not an abrupt return to full strength.
PLUME_CRATON_WEAKENING = 0.75
# Stream ceiling: the knob's upper bound for thinning. Changing `plume_multiplier`
# never reshuffles kept events because every draw below is index-keyed (#579).
PLUME_MULTIPLIER_MAX = 4.0
# Breakup stress (RFC 0007, #582): every continental head arrival adds this much to
# its host plate's `internal_stress` (merge_split.py), pushing plume-affected
# continents toward the split gates. Default 1.0 against STRESS_CERTAIN_RIFT 3.0;
# provisional -- tuned in #584. Oceanic heads add nothing (no craton to weaken, and
# oceanic lithosphere already rifts freely).
PLUME_HEAD_STRESS = 1.0
COX_LOG_SIGMA = 1.0

HEAD_STREAM = "head"
TAIL_STREAM = "tail"


@dataclass
class MantlePlume:
    """A discrete long-lived plume fixed in the world (mantle) frame (#579).

    All defaults so old saves load (see docs/save-compatibility.md): a plume only
    gains a tail after its head stage, and a tail-only plume carries no head.
    """

    plume_id: int = -1
    position: tuple[float, float, float] = (0.0, 0.0, 1.0)
    born_years: float = 0.0
    head_volume_m3: float = 0.0
    head_remaining_m3: float = 0.0
    head_end_years: float = 0.0
    tail_flux_m3_per_yr: float = 0.0
    tail_end_years: float = 0.0
    # Host type judged once at first emplacement (RFC: at born_years). Plain default
    # so old saves load; None means "not judged yet" (#579).
    oceanic_host: bool | None = None

    @property
    def is_tail_only(self) -> bool:
        return self.head_volume_m3 <= 0.0

    @property
    def tail_start_years(self) -> float:
        return self.head_end_years if not self.is_tail_only else self.born_years


def _base_rate_per_myr(stream: str) -> float:
    return PLUME_HEAD_RATE_PER_MYR if stream == HEAD_STREAM else PLUME_TAIL_ONLY_RATE_PER_MYR


def _max_rate_per_myr(stream: str) -> float:
    return _base_rate_per_myr(stream) * PLUME_MULTIPLIER_MAX * PLUME_RATE_MODULATION_MAX


def _event_rng(seed: int, stream: str, index: int) -> np.random.Generator:
    """The event's private generator: every draw below is a pure function of the
    event index, so step partitions never change events (#579)."""
    digest = hashlib.sha256(f"{seed}/plume/{stream}/{index}".encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "big"))


def _sample_gap_years(seed: int, stream: str, index: int) -> float:
    rng = _event_rng(seed, stream, index)
    return float(rng.exponential(1.0 / _max_rate_per_myr(stream))) * 1e6


def _modulation(seed: int, stream: str, t_years: float, cox: bool) -> float:
    """Arrival-rate modulation at time t: 1 stationary, or a seed-keyed log-normal
    series (correlation time PLUME_RATE_CORRELATION_MYR, clipped) for the Cox
    sensitivity variant (#579)."""
    if not cox:
        return 1.0
    slot = int(math.floor(t_years / (PLUME_RATE_CORRELATION_MYR * 1e6)))
    draw = _event_rng(seed, f"cox/{stream}", slot).standard_normal()
    return min(math.exp(COX_LOG_SIGMA * draw - 0.5 * COX_LOG_SIGMA**2), PLUME_RATE_MODULATION_MAX)


def _uniform_sphere(rng: np.random.Generator) -> np.ndarray:
    vec = rng.standard_normal(3)
    norm = float(np.linalg.norm(vec))
    if norm == 0.0:
        return np.array([0.0, 0.0, 1.0])
    return vec / norm


def _sample_position(rng: np.random.Generator, centers: list) -> tuple[float, float, float]:
    """Plume source: upwelling-biased by analogy with LIPs favouring deep-mantle
    margins, plus a uniform floor so not every plume sits on an upwelling (#579)."""
    up = [(np.asarray(c.position, dtype=float), c.strength, c.falloff) for c in centers if c.strength > 0.0]
    total = float(sum(s for _, s, _ in up))
    floor = total * PLUME_UNIFORM_SHARE / (1.0 - PLUME_UNIFORM_SHARE) if total > 0.0 else 1.0
    if total <= 0.0 or rng.random() < floor / (total + floor):
        return tuple(float(v) for v in _uniform_sphere(rng))
    pick = rng.random() * total
    for pos, strength, falloff in up:
        pick -= strength
        if pick <= 0.0:
            axis = pos / np.linalg.norm(pos)
            ref = np.array([0.0, 0.0, 1.0]) if abs(axis[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
            t1 = np.cross(axis, ref)
            t1 /= np.linalg.norm(t1)
            t2 = np.cross(axis, t1)
            offset = falloff * (rng.standard_normal() * t1 + rng.standard_normal() * t2)
            placed = axis + offset
            return tuple(float(v) for v in placed / np.linalg.norm(placed))
    pos, _, _ = up[-1]
    return tuple(float(v) for v in pos / np.linalg.norm(pos))


def _log_uniform(rng: np.random.Generator, low: float, high: float) -> float:
    return float(math.exp(rng.uniform(math.log(low), math.log(high))))


def _sample_event(
    seed: int,
    stream: str,
    index: int,
    born_years: float,
    multiplier: float,
    cox: bool,
    centers: list,
) -> MantlePlume | None:
    """The event at index, or None when thinned out.

    Replays the gap draw first so the private generator stays aligned with
    `_sample_gap_years`; kept events at a lower multiplier are then a strict subset
    of those at a higher one (#579).
    """
    rng = _event_rng(seed, stream, index)
    rng.exponential(1.0 / _max_rate_per_myr(stream))
    ceiling = PLUME_MULTIPLIER_MAX * PLUME_RATE_MODULATION_MAX
    if rng.random() >= multiplier * _modulation(seed, stream, born_years, cox) / ceiling:
        return None
    position = _sample_position(rng, centers)
    if stream == HEAD_STREAM:
        volume = _log_uniform(rng, PLUME_HEAD_VOLUME_MIN_KM3, PLUME_HEAD_VOLUME_MAX_KM3) * KM3_TO_M3
        head_end = born_years + PLUME_HEAD_DURATION_MYR * 1e6
        flux = _log_uniform(rng, PLUME_TAIL_FLUX_MIN_KM3_PER_YR, PLUME_TAIL_FLUX_MAX_KM3_PER_YR) * KM3_TO_M3
        return MantlePlume(
            plume_id=index,
            position=position,
            born_years=born_years,
            head_volume_m3=volume,
            head_remaining_m3=volume,
            head_end_years=head_end,
            tail_flux_m3_per_yr=flux,
            tail_end_years=head_end + PLUME_TAIL_LIFETIME_MYR * 1e6,
        )
    flux = _log_uniform(rng, PLUME_TAIL_FLUX_MIN_KM3_PER_YR, PLUME_TAIL_FLUX_MAX_KM3_PER_YR) * KM3_TO_M3
    # Negative ids: each stream numbers from zero, so tail-only plumes take the
    # mirrored range to keep event ids unambiguous (#579).
    return MantlePlume(
        plume_id=-(index + 1),
        position=position,
        born_years=born_years,
        tail_flux_m3_per_yr=flux,
        tail_end_years=born_years + PLUME_TAIL_LIFETIME_MYR * 1e6,
    )


def consume_arrivals(
    seed: int,
    multiplier: float,
    cox: bool,
    centers: list,
    state: dict,
    _t0_years: float,
    t1_years: float,
) -> list[MantlePlume]:
    """Consume every event with arrival time in (t0, t1], advancing `state` in place.

    `state` maps each stream to `[next_index, next_arrival_years]`. Arrivals at or
    before t0 were consumed by earlier steps, so only t1 is read here; t0 is
    documented for the interval contract and ignored (#579).
    """
    fresh: list[MantlePlume] = []
    for stream in (HEAD_STREAM, TAIL_STREAM):
        if stream not in state:
            state[stream] = [0, _sample_gap_years(seed, stream, 0)]
        index, arrival = state[stream]
        while arrival <= t1_years:
            plume = _sample_event(seed, stream, index, arrival, multiplier, cox, centers)
            if plume is not None:
                fresh.append(plume)
            index += 1
            arrival = arrival + _sample_gap_years(seed, stream, index)
        state[stream] = [index, arrival]
    return fresh


def _stage_volumes(plume: MantlePlume, t0_years: float, t1_years: float) -> tuple[float, float, bool]:
    """(head m³, tail m³, head active) delivered in [t0, t1] by interval overlap.

    Delivery never exceeds the scheduled volume whatever the step lengths (#579).
    """
    head_vol = 0.0
    head_active = False
    if not plume.is_tail_only and plume.head_remaining_m3 > 0.0:
        overlap = max(0.0, min(t1_years, plume.head_end_years) - max(t0_years, plume.born_years))
        duration = plume.head_end_years - plume.born_years
        if overlap > 0.0 and duration > 0.0:
            head_vol = min(plume.head_remaining_m3, plume.head_volume_m3 * overlap / duration)
            head_active = True
    tail_vol = 0.0
    if plume.tail_flux_m3_per_yr > 0.0:
        overlap = max(0.0, min(t1_years, plume.tail_end_years) - max(t0_years, plume.tail_start_years))
        if overlap > 0.0:
            tail_vol = plume.tail_flux_m3_per_yr * overlap
    return head_vol, tail_vol, head_active


def _kernel_weights(distances_m: np.ndarray, radius_m: float, gaussian: bool) -> np.ndarray:
    if gaussian:
        return np.exp(-((distances_m / radius_m) ** 2))
    return np.exp(-((distances_m / radius_m) ** 4))


def _unit_vectors(plate) -> np.ndarray:
    """Unit world-frame node vectors in the plate's node order."""
    points = plate.all_points_and_elevation()[0]
    norms = np.linalg.norm(points, axis=1, keepdims=True)
    return points / np.maximum(norms, 1e-9)


def _heat_top_up_myr(plume_born_years: float, t0_years: float) -> float:
    """Pre-decay heat top-up for one active head (RFC 0007, #582).

    A head born mid-step has heated for only part of it, so topping to
    `PLUME_HEAT_MYR` before this step's decay would shorten its halo by the
    pre-birth time (partition dependence). Topping to
    `PLUME_HEAT_MYR + (born - t0)` makes the post-decay halo exactly
    `PLUME_HEAT_MYR - (t1 - born)`; older heads top up to the full countdown
    exactly as before.
    """
    if plume_born_years <= t0_years:
        return PLUME_HEAT_MYR
    return PLUME_HEAT_MYR + (plume_born_years - t0_years) / 1e6


def _credit_head_stress(plates_data: list, continental: bool, host_plate: int) -> bool:
    """Add PLUME_HEAD_STRESS to the continental head's host plate (RFC 0007, #582).

    `host_plate` is a plate id, not a `plates_data` position (ids go sparse after
    merges). Oceanic heads and unknown hosts credit nothing. Returns True when a
    plate was credited.
    """
    if not continental or host_plate < 0:
        return False
    for entry in plates_data:
        if entry["plate"].plate_id == host_plate:
            host = entry["plate"]
            host.set_internal_stress(host.internal_stress + PLUME_HEAD_STRESS)
            return True
    return False


def _credit_head_arrival(plates_data: list, continental: bool, host_plate: int, result: dict) -> None:
    """Credit one continental head arrival (RFC 0007, #582): host stress plus the
    `credited_heads` step-diagnostic count. Split out of `apply` for complexity."""
    if _credit_head_stress(plates_data, continental, host_plate):
        result["credited_heads"] += 1


def _center_host(plates: list, position: np.ndarray) -> tuple[bool, int]:
    """Host at the plume centre: (is continental, plate id) of the nearest node."""
    from .elevation_lines import effective_is_continental_from_codes

    best_cos = -2.0
    host_continental = False
    host_plate = -1
    for plate in plates:
        if plate.node_count() == 0:
            continue
        cosang = _unit_vectors(plate) @ position
        top = int(np.argmax(cosang))
        if float(cosang[top]) > best_cos:
            best_cos = float(cosang[top])
            codes = plate.collect("crust_type_code")
            host_continental = bool(
                effective_is_continental_from_codes(codes[top : top + 1], plate.crust_type == "continental")[0]
            )
            host_plate = plate.plate_id
    return host_continental, host_plate


def _footprint_share(
    unit: np.ndarray,
    areas: np.ndarray,
    position: np.ndarray,
    radius_m: float,
    gaussian: bool,
    truncate_radii: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """One plate's share of a stage footprint: (cell mask, kernel weights, norm).

    The caller partitions the stage volume across plates by norm, so a footprint
    straddling a boundary still places its volume exactly once (#579).
    """
    from .elevation_lines import PLANET_RADIUS_KM

    earth_m = PLANET_RADIUS_KM * 1000.0
    dist = np.arccos(np.clip(unit @ position, -1.0, 1.0)) * earth_m
    mask = dist <= truncate_radii * radius_m
    if not np.any(mask):
        # Nowhere in reach (only on degenerate cover): park the whole stage on the
        # nearest cell with a uniform weight so scheduled volume is never lost (#579).
        mask = np.zeros(len(dist), dtype=bool)
        mask[int(np.argmin(dist))] = True
        return mask, np.ones(1), float(areas[mask].sum())
    weights = _kernel_weights(dist[mask], radius_m, gaussian)
    return mask, weights, float(np.sum(weights * areas[mask]))


def _deposit_on_plate(world, plate, fields: dict, mask: np.ndarray, weights: np.ndarray, norm: float, volume_m3: float) -> float:
    """Add this plate's partitioned share; returns the volume placed here.

    A zero norm (footprint cells whose weights all underflowed far from the
    plume) places nothing -- the caller partitions around such shares (#579).

    Oceanic-coded cells gain Hc; continental-coded cells additionally book the
    added crust as juvenile. Never sets `is_volcano`: plume relief comes from the
    plume's own budget alone (see `volcanism.apply_volcanic_activity`) (#579).
    """
    from . import continental_ledger, lithosphere, rheology
    from .elevation_lines import ELEV_CHANGE_VOLCANO, effective_is_continental_from_codes

    if norm <= 0.0:
        return 0.0
    thick = volume_m3 * weights / norm
    hc = fields["hc"]
    hm = fields["hm"]
    density = lithosphere.node_crust_density(fields["codes"], plate.crust_type)
    # Incremental shift from THIS deposit only: isostasy is pointwise, so later
    # deposits in the same step each add their own marginal relief (#579).
    hc_pre = hc[mask].copy()
    hc[mask] = hc_pre + thick
    shift = lithosphere.isostatic_elevation(hc[mask], hm[mask], density[mask]) - lithosphere.isostatic_elevation(
        hc_pre, hm[mask], density[mask]
    )
    elevation = fields["elevation"]
    elevation[mask] = rheology.clip_elevation_bounds(elevation[mask] + shift)
    fields["reason"][mask] = ELEV_CHANGE_VOLCANO
    continental = effective_is_continental_from_codes(fields["codes"], plate.crust_type == "continental")
    juvenile = np.zeros(len(hc))
    juvenile[mask & continental] = thick[continental[mask]]
    if np.any(juvenile):
        continental_ledger.add_material_thickness(
            world, plate, juvenile, "juvenile_additions_m3"
        )
    return float(np.sum(thick * fields["areas"][mask]))


def _stage_kernel(plume: MantlePlume, stage: str, oceanic_host: bool) -> tuple[float, bool, float]:
    """(emplacement radius, gaussian, truncation) for one stage volume (#579)."""
    if stage == "head":
        if oceanic_host:
            return head_kernel_radius_m(head_footprint_area_m2(plume.head_volume_m3)), False, 4.0
        return PLUME_HEAD_RADIUS_KM * 1000.0, True, 4.0
    return PLUME_TAIL_RADIUS_KM * 1000.0 / 2.0, False, 2.0


def apply(world, years: float) -> dict:
    """Emplace one step of plume head/tail/heat; returns step diagnostics (#579).

    Interval-overlap delivery makes any step length exact; `plume_multiplier = 0`
    is a strict no-op reproducing a matched plume-off run. Stage volumes partition
    across plates by footprint weight, so a boundary-straddling footprint still
    places its volume exactly once. Plumes whose tail has ended are dropped at the
    end of the step.

    Heat/stress take effect the step *after* arrival: deform and topology run
    before this pass, so a fresh head's halo and stress credit land on next
    step's rifting. Harmless at 5 Myr steps against a 20 Myr halo.
    """
    from .elevation_lines import PLANET_RADIUS_KM, line_spacing_rad

    result = {"new_heads": 0, "new_tails": 0, "head_volume_m3": 0.0, "tail_volume_m3": 0.0, "credited_heads": 0}
    if years <= 0.0 or not (world.plume_multiplier > 0.0) or not world.simulate_plate_movement:
        # `not (... > 0.0)` rather than `<= 0.0` so a NaN multiplier no-ops instead
        # of producing at ceiling rate; set_controls also rejects non-finite (#579).
        return result
    t1_years = world.elapsed_years
    t0_years = t1_years - years
    fresh = consume_arrivals(
        world.seed,
        world.plume_multiplier,
        world.plume_cox_modulation,
        world.mantle_centers,
        world.plume_stream,
        t0_years,
        t1_years,
    )
    world.mantle_plumes.extend(fresh)
    spacing = line_spacing_rad(world.node_density)
    plates_data = []
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        plates_data.append(
            {
                "plate": plate,
                "unit": _unit_vectors(plate),
                "areas": plate.accounting_areas_m2(spacing),
                "hc": plate.collect("crustal_thickness_m"),
                "hm": plate.collect("mantle_lithosphere_thickness_m"),
                "codes": plate.collect("crust_type_code"),
                "elevation": plate.collect("elevation"),
                "reason": plate.collect("elev_change_reason"),
                "heat": plate.collect("plume_heat_myr"),
            }
        )
    for plume in fresh:
        if plume.is_tail_only:
            result["new_tails"] += 1
        else:
            result["new_heads"] += 1
            continental, host_plate = _center_host([e["plate"] for e in plates_data], np.asarray(plume.position))
            if plates_data:
                plume.oceanic_host = not continental
            if host_plate >= 0:
                world.log_event(
                    f"Mantle plume head {plume.plume_id}: {plume.head_volume_m3 / 1e15:.1f} M km3 "
                    f"on {'oceanic' if plume.oceanic_host is not False else 'continental'} plate {host_plate}."
                )
            # A plate taking a continental head carries that stress toward breakup
            # (RFC 0007, #582). Per arrival, not per step: the head's thermal halo
            # persists separately as `plume_heat_myr`.
            _credit_head_arrival(plates_data, continental, host_plate, result)
    earth_m = PLANET_RADIUS_KM * 1000.0
    for plume in world.mantle_plumes:
        head_vol, tail_vol, active = _stage_volumes(plume, t0_years, t1_years)
        if head_vol <= 0.0 and tail_vol <= 0.0 and not active:
            continue
        position = np.asarray(plume.position, dtype=float)
        if plume.oceanic_host is None and plates_data:
            continental, _ = _center_host([e["plate"] for e in plates_data], position)
            plume.oceanic_host = not continental
        # No plates: nothing to emplace on, so the kernel choice is dead.
        oceanic = plume.oceanic_host is not False
        # Tail-only plumes never refresh heat: weakening follows the head stage's
        # broad radius (RFC 0007), so an 80 Myr tail outlives its thermal halo (#579).
        for stage, volume in (("head", head_vol), ("tail", tail_vol)):
            if volume <= 0.0:
                continue
            radius_m, gaussian, truncate = _stage_kernel(plume, stage, oceanic)
            shares = []
            total = 0.0
            for entry in plates_data:
                mask, weights, norm = _footprint_share(
                    entry["unit"], entry["areas"], position, radius_m, gaussian, truncate
                )
                shares.append((entry, mask, weights, norm))
                total += norm
            if total <= 0.0:
                continue
            for entry, mask, weights, norm in shares:
                got = _deposit_on_plate(
                    world, entry["plate"], entry, mask, weights, norm, volume * norm / total
                )
                result[f"{stage}_volume_m3"] += got
                if stage == "head":
                    plume.head_remaining_m3 = max(0.0, plume.head_remaining_m3 - got)
        if active:
            top_up = _heat_top_up_myr(plume.born_years, t0_years)
            for entry in plates_data:
                dist = np.arccos(np.clip(entry["unit"] @ position, -1.0, 1.0)) * earth_m
                refresh = dist <= 2.0 * PLUME_HEAD_RADIUS_KM * 1000.0
                entry["heat"][refresh] = np.maximum(entry["heat"][refresh], top_up)
    for entry in plates_data:
        heat = np.maximum(0.0, entry["heat"] - years / 1e6)
        entry["plate"].set_fields_on_plate(
            crustal_thickness_m=entry["hc"],
            elevation=entry["elevation"],
            elev_change_reason=entry["reason"],
            plume_heat_myr=heat,
        )
    world.mantle_plumes = [p for p in world.mantle_plumes if p.tail_end_years > t1_years]
    return result


def diagnostics(world) -> dict:
    """Live plume census for stats and tests (#579)."""
    return {
        "active_plumes": len(world.mantle_plumes),
        "heads": sum(1 for p in world.mantle_plumes if not p.is_tail_only),
        "tails": sum(1 for p in world.mantle_plumes if p.is_tail_only),
        "stream": {name: list(state) for name, state in world.plume_stream.items()},
    }

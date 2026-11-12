"""Experimental river-substep arms for the #867 exploration. Every arm is off by default.

Issue #867 asks whether rivers need intermediate chances to incise and reroute inside a long
tectonic step. The #866 diagnostic found that when a basin's outlet changed between 100 kyr
steps, the old outlet path needed a median cut of 81 m, and that modeled channel growth on
that path was 0 m at the median. The owner's plan is an assumption ladder, run from the same
checkpoint with the same tectonic forcing:

- **Outlet-preservation bound** (`outlet_bound`): a deliberately extreme upper bound. A pit
  whose water left for the ocean last step (a *threatened outlet*) may cut its way out at
  `outlet_rate_m_per_yr` (1 m/yr by default) times the step, stopping at the cut its route
  needs. Everything else is the ordinary breaching rule.
- **Uplift-limited bound** (`uplift_bound`): the same trigger, but each cell on the route may
  be cut only by the ordinary carving limit plus the uplift that cell accumulated since the
  last erosion pass (`channel_reference_elevation_m`). A river that keeps pace with uplift,
  and no more.
- **Simulation-driven march** (`march`): the river incision and routing of one step run in
  `substeps` equal slices, like the glacier march (RFC 0015 D2): each slice re-breaches,
  rebuilds the depression hierarchy and flow directions on the terrain the previous slices
  cut, re-routes the step's steady water sources, and applies stream-power incision for its
  share of the step. Climate, the lake solve and the ice solve stay at their opening,
  full-step values; the march only lets the rivers see their own cuts.

The arms act through two hooks, both no-ops when `ARM.mode == "off"`:
`hydrology.compute_hydrology` passes `outlet_relief` to `breaching.breach_depressions`
(the two bounds), and `erosion.apply_erosion` replaces its one-shot river term with
`march_rivers`'s summed demand (the march). With `substeps == 1` the march is skipped, so a
one-slice march is bit-identical to the baseline.

The knobs live on this module (like `hydrology.MARCH_GLACIERS` and `ice_sheet.MARCH_COST`),
not on `World`: saves stay unchanged and the UI never runs an arm. The driver
`bin/debug/river_substep_arms.py` sets `ARM`, records every step's `LOG` entry, and writes
the checkpoint saves and images. No model change is proposed here.
"""

from __future__ import annotations

import time
from collections import deque
from dataclasses import dataclass

import numpy as np

from ..plates import (
    collect_all_breach_notch_depth,
    collect_all_channel_reference_elevation,
    collect_all_lake_ice,
)
from ..surface_fields import CHANNEL_REFERENCE_VALID_BELOW_M

MODES = ("off", "outlet_bound", "uplift_bound", "march")
TRIGGERS = ("threatened", "all")


@dataclass(frozen=True)
class RiverSubstepArm:
    """One arm of the #867 comparison. `mode` picks the arm; the other fields only matter to
    the arms named in their comments."""

    mode: str = "off"
    # outlet_bound / uplift_bound: which pits get the extra carving allowance. "threatened"
    # is a pit whose water reached the ocean last step; "all" is every land pit, which bounds
    # what any incision rule could do at all.
    trigger: str = "threatened"
    # outlet_bound: the allowance, metres per year of step (1 m/yr is the owner's extreme bound).
    outlet_rate_m_per_yr: float = 1.0
    # march: slices per step. 1 is the baseline.
    substeps: int = 1

    def __post_init__(self) -> None:
        if self.mode not in MODES:
            raise ValueError(f"unknown river-substep mode {self.mode!r}; known: {MODES}")
        if self.trigger not in TRIGGERS:
            raise ValueError(f"unknown trigger {self.trigger!r}; known: {TRIGGERS}")
        if self.substeps < 1:
            raise ValueError("substeps must be at least 1")

    @property
    def relieves_outlets(self) -> bool:
        return self.mode in ("outlet_bound", "uplift_bound")

    @property
    def label(self) -> str:
        """The command-line spelling `parse_arm` reads back."""
        if self.mode == "march":
            return f"march:{self.substeps}"
        if self.relieves_outlets:
            name = self.mode.replace("_", "-")
            return name if self.trigger == "threatened" else f"{name}:{self.trigger}"
        return "off"


ARM = RiverSubstepArm()

# One entry per hydrology/erosion step while an arm is on (the driver drains it). Bounded,
# like `ice_sheet.MARCH_LOG`, so a long interactive session can't grow it without limit.
LOG: deque[dict] = deque(maxlen=64)


def parse_arm(text: str) -> RiverSubstepArm:
    """`off`, `outlet-bound[:all]`, `uplift-bound[:all]` or `march:<substeps>`."""
    name, has_option, option = text.strip().partition(":")
    name = name.replace("_", "-")
    usage = f"use off, outlet-bound[:all], uplift-bound[:all] or march:<n>, not {text!r}"
    # A mistyped spelling must fail loudly: a preregistered run can't quietly become `off`.
    if name == "off" and not has_option:
        return RiverSubstepArm()
    if name == "march" and option.isdigit():
        return RiverSubstepArm(mode="march", substeps=int(option))
    if name in ("outlet-bound", "uplift-bound") and (not has_option or option in TRIGGERS):
        return RiverSubstepArm(mode=name.replace("-", "_"), trigger=option or "threatened")
    raise ValueError(usage)


def set_arm(arm: RiverSubstepArm) -> None:
    global ARM
    ARM = arm
    LOG.clear()


# --- Threatened outlets ------------------------------------------------------------------------


def exported_to_ocean(flow_target: np.ndarray, is_ocean: np.ndarray) -> np.ndarray:
    """Whether each node's routed water reached the ocean: the end of its `flow_target` chain
    is an ocean node. A chain that ends at a land sink, or runs in a cycle, did not. Pointer
    jumping, so it takes about log2(longest path) passes; the pass count is capped so a cycle
    can't spin it."""
    flow_target = np.asarray(flow_target, dtype=np.int64)
    is_ocean = np.asarray(is_ocean, dtype=bool)
    n = len(flow_target)
    if n == 0:
        return np.zeros(0, dtype=bool)
    receiver = np.where((flow_target >= 0) & (flow_target < n), flow_target, np.arange(n))
    receiver = np.where(is_ocean, np.arange(n), receiver)
    for _ in range(int(np.ceil(np.log2(n))) + 2):
        following = receiver[receiver]
        if np.array_equal(following, receiver):
            break
        receiver = following
    settled = receiver[receiver] == receiver
    return settled & is_ocean[receiver]


def threatened_mask(world, points: np.ndarray) -> np.ndarray | None:
    """Nodes whose nearest node in last step's routing (`world.hydrology_cache`) sent its
    water to the ocean, resampled onto this step's `points` the same way `sample_is_ocean`
    resamples the ocean mask. None before the first step, when there is no previous routing.
    A pit among these nodes is a threatened outlet: last step the water here got out."""
    from . import hydrology

    hydro = getattr(world, "hydrology_cache", None)
    flow_target = getattr(hydro, "flow_target", None)
    if hydro is None or flow_target is None or len(hydro.points) == 0 or len(flow_target) != len(hydro.points):
        return None
    exported = exported_to_ocean(flow_target, hydro.is_ocean)
    idx = hydrology._nearest_hydro_node_idx(hydro, points, getattr(world, "node_cloud_resample_mode", "kdtree"))
    return exported[idx]


def uplift_since_reference_m(elevation: np.ndarray, plates_in_order) -> np.ndarray:
    """How much each node rose since the last erosion pass recorded its
    `channel_reference_elevation_m` (tectonics, volcanism and isostasy this step), clipped at
    zero. A node with no valid reference counts as no uplift, as erosion.py's channel fade
    does."""
    reference = collect_all_channel_reference_elevation(plates_in_order)
    rise = np.where(reference < CHANNEL_REFERENCE_VALID_BELOW_M, np.asarray(elevation, dtype=float) - reference, 0.0)
    return np.clip(rise, 0.0, None)


def outlet_relief(world, points: np.ndarray, elevation: np.ndarray, plates_in_order, is_ocean: np.ndarray, years: float):
    """The `breaching.OutletRelief` for the two bound arms, or None when neither is on (or
    there is no previous routing to name a threatened outlet)."""
    from . import breaching

    arm = ARM
    if not arm.relieves_outlets:
        return None
    n = len(points)
    if arm.trigger == "all":
        eligible = ~np.asarray(is_ocean, dtype=bool)
    else:
        eligible = threatened_mask(world, points)
        if eligible is None:
            return None
        eligible = eligible & ~np.asarray(is_ocean, dtype=bool)
    if arm.mode == "outlet_bound":
        extra = np.full(n, arm.outlet_rate_m_per_yr * max(float(years), 0.0))
    else:
        extra = uplift_since_reference_m(elevation, plates_in_order)
    return breaching.OutletRelief(eligible=eligible, extra_cut_m=extra)


def note_breach(relief, breach, elevation: np.ndarray, is_ocean: np.ndarray, neighbor_idx: np.ndarray) -> None:
    """Record what the bound arms did this step: how many pits were eligible, how many were
    breached only because of the allowance, and the notch those breaches needed."""
    from . import breaching

    if relief is None or breach is None:
        return
    pits = breaching._pits(np.asarray(elevation, dtype=float), np.asarray(is_ocean, dtype=bool), neighbor_idx)
    LOG.append({
        "mode": ARM.mode,
        "trigger": ARM.trigger,
        "pits": int(len(pits)),
        "eligible_pits": int(np.count_nonzero(relief.eligible[pits])) if len(pits) else 0,
        "breached_pits": int(len(breach.breached_pits)),
        "relieved_pits": int(len(breach.relieved_pits)),
        "notch_sum_m": float(breach.notch_m.sum()),
        "notch_max_m": float(breach.notch_m.max()) if len(breach.notch_m) else 0.0,
        "extra_cut_max_m": float(np.max(relief.extra_cut_m)) if len(relief.extra_cut_m) else 0.0,
    })


# --- The river march ----------------------------------------------------------------------------


@dataclass
class RiverMarch:
    """What `march_rivers` hands back to `erosion.apply_erosion`, all per node unless noted."""

    # The stream-power demand summed over the slices. Replaces the one-shot `river` term and
    # goes through erosion.py's usual cover-first, headroom-capped, craton-resisted pipeline.
    demand_m: np.ndarray
    # The cut the march applied to its working terrain, slice by slice, each slice capped at
    # that slice's drop to the lowest neighbour. erosion.py lets the step's cap rise to this
    # where it exceeds the opening drop, since the march already honoured the drop as the
    # valley deepened.
    working_cut_m: np.ndarray
    # Routing on the terrain the last slice left: what deposition and channel width read.
    flow_target: np.ndarray
    discharge_m3_s: np.ndarray
    # Breach notches the later slices invented (the opening slice's are already on `hydro`).
    extra_notch_m: np.ndarray
    substeps: int
    seconds: float


def _reroute(world, hydro, elevation, channel, cover, notch, remaining_years, solved_water_m, areas):
    """One slice's routing on the march's working terrain: breach, rebuild the hierarchy,
    choose flow directions against the held lake surface, and re-route the step's steady
    sources. Returns `(flow_target, discharge_m3_s, breach)`."""
    from . import breaching, hydrology, lakes

    is_ocean = hydro.is_ocean
    neighbor_idx = hydro.neighbor_idx
    breach = breaching.breach_depressions(
        elevation, is_ocean, neighbor_idx, channel, getattr(hydro, "breach_carve_rate", None), remaining_years,
        water=getattr(hydro, "breach_water_balance", None), prior_notch_m=notch, mobile_cover_m=cover,
    )
    interface = breaching.interface_pass_elevation(breach.passage_m, neighbor_idx)
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx, interface_pass_elevation=interface)
    hypsometry = lakes.CellHypsometry.from_terrain(elevation, areas, neighbor_idx)
    filled, spill_target, _ = lakes.compute_spill_routing(forest, elevation, solved_water_m, hypsometry)
    offset = np.zeros(len(elevation))
    wet = np.flatnonzero(solved_water_m > 0.0)
    if len(wet):
        offset[wet] = hypsometry.level_from_depth(solved_water_m[wet], wet) - elevation[wet]
    flow_target, _ = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, offset, filled, spill_target, channel)
    routed = hydrology.route_downstream(
        elevation, is_ocean, flow_target, hydro.river_source_m3_per_yr, loss_fraction=getattr(hydro, "river_transit_loss", None)
    )
    discharge = np.where(is_ocean, routed.deposited, routed.through_flux) / hydrology.SECONDS_PER_YEAR
    return flow_target, discharge, breach


def march_rivers(
    world,
    hydro,
    points: np.ndarray,
    elevation: np.ndarray,
    prior_channel_depth: np.ndarray,
    prior_hc: np.ndarray,
    prior_hm: np.ndarray,
    cover_m: np.ndarray,
    craton_keep: np.ndarray,
    removable_m: np.ndarray,
    years: float,
    multiplier: float,
) -> RiverMarch | None:
    """March the step's river incision in `ARM.substeps` slices (see the module docstring).
    None when the march is off or has one slice, or when `hydro` carries no routable sources
    (a hand-built cache).

    Slice 0 reads the opening routing exactly as the baseline does, over `years / substeps`.
    Every later slice re-routes on the terrain the earlier slices cut (`_reroute`), with the
    time left in the step as its breaching budget, so a pit the opening routing could not open
    is judged again once the rivers have lowered its route. Each slice applies the capped cut
    to the working terrain: the cover first, then substrate under craton resistance, never
    past the slice's drop to the lowest neighbour or the column's headroom. The surface drops
    by the Airy response to the removed rock, as the end-of-step write-back will, while the
    channel deepens by the full cut, as `breaching.channel_passage_elevation` expects.

    What the march holds fixed, and why: climate, the lake solve (its levels and spills) and
    the ice solve are full-step quantities already computed on the opening routing; the water
    ledger closes on that routing. Re-solving lakes each slice is the broad step-integrated
    landscape solve of #297, not this issue's targeted question. A basin the march opens
    mid-step keeps its opening water until the next step's solve drains it, the same one-step
    lag a breach has today."""
    from .. import lithosphere
    from . import erosion

    arm = ARM
    if arm.mode != "march" or arm.substeps <= 1:
        return None
    n = len(points)
    sources = getattr(hydro, "river_source_m3_per_yr", None)
    if sources is None or len(sources) != n or hydro.hypsometry is None or getattr(hydro, "breach_carve_rate", None) is None:
        LOG.append({"mode": "march", "substeps": arm.substeps, "skipped": "no routable sources on the hydrology cache"})
        return None
    started = time.perf_counter()
    plates_in_order = hydro.plates_in_order
    is_ocean = np.asarray(hydro.is_ocean, dtype=bool)
    areas = hydro.hypsometry.area_m2
    rho_c = np.concatenate([lithosphere.node_crust_density(p.collect("crust_type_code"), p.crust_type) for p in plates_in_order])
    has_column = prior_hc > 0.0
    lake_ice = collect_all_lake_ice(plates_in_order)
    solved_water = np.asarray(hydro.lake_depth, dtype=float) + (lake_ice if len(lake_ice) == n else 0.0)
    opening_notch = hydro.breach_notch_m if len(hydro.breach_notch_m) == n else np.zeros(n)
    notch = collect_all_breach_notch_depth(plates_in_order) + opening_notch

    elev = np.asarray(elevation, dtype=float).copy()
    channel = np.asarray(prior_channel_depth, dtype=float).copy()
    cover = np.clip(np.asarray(cover_m, dtype=float), 0.0, None).copy()
    room = np.asarray(removable_m, dtype=float).copy()
    hc = np.asarray(prior_hc, dtype=float).copy()
    flow_target = np.asarray(hydro.flow_target).copy()
    discharge = np.asarray(hydro.discharge_m3_s, dtype=float).copy()
    surge = hydro.breach_surge_accum if hydro.breach_surge_accum is not None and len(hydro.breach_surge_accum) == n else np.zeros(n)
    extra_notch = np.zeros(n)
    demand_total = np.zeros(n)
    cut_total = np.zeros(n)
    later_breached = 0
    dt = years / arm.substeps
    dt_myr = dt / 1_000_000.0
    for k in range(arm.substeps):
        if k > 0:
            flow_target, discharge, breach = _reroute(
                world, hydro, elev, channel, cover, notch + extra_notch, years - k * dt, solved_water, areas,
            )
            extra_notch = extra_notch + breach.notch_m
            later_breached += int(len(breach.breached_pits))
        slope, drop = erosion.compute_slope(points, elev, world=world)
        demand = np.where(is_ocean, 0.0, erosion._river_erosion_amount(discharge, surge, slope, channel, multiplier, dt_myr))
        limit = np.minimum(drop, room)
        cover_cut, substrate = erosion._strip_mobile_cover(cover, demand)
        cover_cut = np.minimum(cover_cut, limit)
        cut = np.clip(cover_cut + np.minimum(substrate * craton_keep, limit - cover_cut), 0.0, None)
        cover = cover - cover_cut
        room = room - cut
        new_hc = np.where(has_column, hc - cut, hc)
        delta = np.where(
            has_column,
            lithosphere.isostatic_elevation(new_hc, prior_hm, rho_c) - lithosphere.isostatic_elevation(hc, prior_hm, rho_c),
            -cut,
        )
        elev = elev + delta
        hc = new_hc
        channel = np.where(is_ocean, 0.0, np.clip(channel + cut, 0.0, erosion.MAX_CHANNEL_DEPTH_M))
        demand_total += demand
        cut_total += cut
    seconds = time.perf_counter() - started
    rerouted = int(np.count_nonzero(flow_target != hydro.flow_target))
    LOG.append({
        "mode": "march",
        "substeps": arm.substeps,
        "seconds": seconds,
        "rerouted_nodes": rerouted,
        "later_breached_pits": later_breached,
        "extra_notch_sum_m": float(extra_notch.sum()),
        "demand_sum_m": float(demand_total.sum()),
        "working_cut_sum_m": float(cut_total.sum()),
    })
    return RiverMarch(
        demand_m=demand_total,
        working_cut_m=cut_total,
        flow_target=flow_target,
        discharge_m3_s=discharge,
        extra_notch_m=extra_notch,
        substeps=arm.substeps,
        seconds=seconds,
    )

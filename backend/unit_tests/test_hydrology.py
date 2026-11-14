import pickle

import numpy as np
import pytest
from app import geometry, lithosphere, trace
from app.elevation_lines import line_spacing_rad
from app.hydroclimate import breaching, erosion, eustasy, hydrology, ice_sheet, lakes
from app.plates import collect_all_accounting_areas_m2
from app.world import World, generate_world, step_world

from .quad_fixtures import quad_plate
from .river_fixtures import DRY_SPILL_PATH, segment_pairs, spill_hop_fields


def _row(count: int) -> np.ndarray:
    """Column offsets for a row of `count` cells centred on the face."""
    return np.arange(count) - count // 2


def _flow_line_plate(plate_id, theta, elevation):
    """A continental row of cells, one per entry of `theta` (only its length matters)."""
    return quad_plate(plate_id, "continental", columns=_row(len(theta)), elevation=np.asarray(elevation, dtype=float))


def _flow_line_plate_with_lake(plate_id, theta, elevation, lake_depth, glacier_depth=None):
    kw = {}
    if glacier_depth is not None:
        kw["glacier_depth"] = np.asarray(glacier_depth, dtype=float)
    return quad_plate(
        plate_id, "continental", columns=_row(len(theta)), elevation=np.asarray(elevation, dtype=float),
        lake_depth=np.asarray(lake_depth, dtype=float), **kw,
    )


def _normalize(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def test_precipitation_source_partition_accounts_for_dry_wet_and_partial_cells():
    precipitation_mm = np.array([800.0, 800.0, 800.0])
    evaporation_mm = np.array([600.0, 600.0, 600.0])
    areas_m2 = np.array([1.0e6, 2.0e6, 3.0e6])
    wet_fraction = np.array([0.0, 1.0, 0.25])

    wet_rain, dry_runoff, dry_et = hydrology._partition_precipitation_sources(
        precipitation_mm, evaporation_mm, areas_m2, wet_fraction
    )

    liquid_input = precipitation_mm / 1000.0 * areas_m2
    assert wet_rain[0] == 0.0 and dry_runoff[0] > 0.0
    assert wet_rain[1] == pytest.approx(liquid_input[1])
    assert dry_runoff[1] == 0.0 and dry_et[1] == 0.0
    assert wet_rain[2] == pytest.approx(liquid_input[2] * 0.25)
    assert dry_runoff[2] == pytest.approx(
        breaching.runoff_mm(precipitation_mm[2:3], evaporation_mm[2:3])[0] / 1000.0 * areas_m2[2] * 0.75
    )
    assert wet_rain + dry_runoff + dry_et == pytest.approx(liquid_input)


def test_external_water_sources_keep_their_own_transit_loss_rules():
    elevation = np.array([2.0, 1.0, 0.0])
    is_ocean = np.zeros(3, dtype=bool)
    targets = np.array([1, 2, -1])
    dry_loss = np.array([0.1, 0.2, 0.3])
    transit_loss = np.full(3, 0.2)
    routed = hydrology._route_external_water_sources(
        elevation, is_ocean, targets, np.array([100.0, 0.0, 0.0]),
        np.array([100.0, 0.0, 0.0]), np.array([100.0, 0.0, 0.0]), dry_loss, transit_loss,
    )

    assert routed["dry_runoff"][1][2] == pytest.approx(100.0 * 0.9 * 0.8 * 0.7)
    assert routed["meltwater"][1][2] == pytest.approx(100.0 * 0.8**3)
    assert routed["wet_precip"][1][2] == pytest.approx(100.0)


def test_compute_hydrology_exposes_wet_precipitation_once_in_external_supply():
    theta = 0.05 * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    prev_lake_depth = np.zeros(40)
    prev_lake_depth[28:31] = 35.0
    world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, prev_lake_depth)])
    precipitation = np.full(40, 800.0)
    fields = hydrology.compute_hydrology(world, precipitation, np.full(40, 15.0), years=100_000)

    wet_fraction = hydrology._lagged_lake_wet_fraction(prev_lake_depth, fields.hypsometry)
    wet_fraction[fields.is_ocean] = 0.0
    areas = fields.hypsometry.area_m2
    expected_direct_rain = float((precipitation / 1000.0 * areas * wet_fraction).sum())
    assert fields.lake_external_wet_precip_m3_per_yr.sum() == pytest.approx(expected_direct_rain)
    assert fields.lake_external_supply_m3_per_yr == pytest.approx(
        fields.lake_external_wet_precip_m3_per_yr
        + fields.lake_external_dry_runoff_m3_per_yr
        + fields.lake_external_meltwater_m3_per_yr
    )


@pytest.mark.usefixtures("substep_lakes")
def test_compute_hydrology_passes_capped_substeps_to_lake_solve(monkeypatch):
    # #875 production wiring: with the knob on, compute_hydrology must hand the
    # lake solver lake_substeps(years), not the old single slice.
    seen = {}
    original = lakes.solve_lakes_prognostic

    def spy(*args, **kwargs):
        seen["substeps"] = kwargs.get("substeps", 1)
        return original(*args, **kwargs)

    monkeypatch.setattr(lakes, "solve_lakes_prognostic", spy)
    theta = 0.05 * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    prev_lake_depth = np.zeros(40)
    prev_lake_depth[28:31] = 35.0
    world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, prev_lake_depth)])
    hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=100_000)
    assert seen["substeps"] == hydrology.lake_substeps(100_000.0) == 100


@pytest.mark.usefixtures("march_glaciers")
def test_lake_takes_full_wet_share_including_snow():
    # The lake's direct supply is its FULL wet share -- rain plus snow -- exactly, not a
    # re-weighted slice of dry rain (#806 review). Wet snow has no glacier to join and
    # no D4 cap hand-off yet, so routing only the rain share deletes it; and it cannot
    # go through the meltwater stream instead, which drops source-pass sink settlement
    # -- exactly where lake snow starts. The snow share is still booked, as in-step
    # snowmelt in the ice ledger, so the conservation checks count it.
    theta = 0.05 * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    prev_lake_depth = np.zeros(40)
    prev_lake_depth[28:31] = 35.0
    world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, prev_lake_depth)])
    precipitation = np.full(40, 800.0)
    years = 100_000.0
    fields = hydrology.compute_hydrology(
        world, precipitation, np.full(40, -3.0), years=years,
        seasonal_amplitude_at_nodes=np.full(40, 10.0),
    )
    wet_fraction = hydrology._lagged_lake_wet_fraction(prev_lake_depth, fields.hypsometry)
    wet_fraction[fields.is_ocean] = 0.0
    areas = fields.hypsometry.area_m2
    snow = ice_sheet.snow_fraction(np.full(40, -3.0), np.full(40, 10.0))
    assert float(snow.sum()) > 0.0  # the fixture must actually have snow to split
    expected_direct = float((precipitation / 1000.0 * areas * wet_fraction).sum())
    assert expected_direct > 0.0
    assert fields.lake_external_wet_precip_m3_per_yr.sum() == pytest.approx(expected_direct)
    # The snow part of that direct supply is booked as in-step snowmelt (ice m^3), on
    # top of whatever bare-ground snow the march itself melted.
    expected_snow_ice_m3 = float((snow * precipitation / 1000.0 * areas * wet_fraction).sum()) * years
    expected_snow_ice_m3 *= ice_sheet.RHO_WATER / ice_sheet.RHO_ICE
    assert fields.ice.ledger.snow_melted_m3 >= expected_snow_ice_m3


def test_seasonal_lake_inputs_leave_dry_land_budyko_pet_unchanged():
    # At 25 C every drop is rain at any seasonal swing (RFC 0015 D1's snow share is 0), so
    # the amplitude must not move dry-land runoff: Budyko PET reads the annual mean.
    theta = 0.05 * np.arange(40)
    elevation = np.linspace(100.0, -100.0, 40)
    world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
    precipitation = np.full(40, 800.0)
    temperature = np.full(40, 25.0)

    steady = hydrology.compute_hydrology(
        world, precipitation, temperature, years=100_000, seasonal_amplitude_at_nodes=np.zeros(40)
    )
    seasonal = hydrology.compute_hydrology(
        world, precipitation, temperature, years=100_000, seasonal_amplitude_at_nodes=np.full(40, 20.0)
    )

    assert seasonal.lake_external_dry_runoff_m3_per_yr == pytest.approx(
        steady.lake_external_dry_runoff_m3_per_yr
    )


def test_lake_temperature_uses_lagged_water_surface_and_connected_ocean_mask():
    elevation = np.array([-100.0, -100.0, 200.0])
    previous_depth = np.array([64.0, 64.0, 0.0])
    hypsometry = lakes.CellHypsometry(
        elevation=elevation,
        area_m2=np.ones(3),
        relief_m=np.array([100.0, 100.0, 0.0]),
        reference_area_m2=1.0,
    )

    water_temperature = hydrology._lake_water_surface_temperature_c(
        air_temperature_at_bed_c=np.full(3, 10.0),
        ocean_temperature_c=np.full(3, 5.0),
        is_ocean=np.array([False, True, False]),
        elevation_m=elevation,
        previous_lake_depth_m=previous_depth,
        hypsometry=hypsometry,
    )

    # A mean depth of 64 m in a cell with 100 m sub-cell relief puts its water surface at
    # -40 m, not bed elevation + mean depth (RFC 0014's shared hypsometry contract).
    # The enclosed below-sea-level lake uses air temperature there; connected ocean uses ocean.
    assert water_temperature[0] == pytest.approx(10.0 - 6.5 * 0.06)
    assert water_temperature[1] == pytest.approx(5.0)
    # A newly formed lake has no previous water level, so this step uses the bed-level air temp.
    assert water_temperature[2] == pytest.approx(10.0)


def _river_inspector_fields():
    # A hand-built 9-node HydrologyFields for the River Inspector's grouping logic, not run
    # through compute_hydrology (small enough to hit its FLOW_NEIGHBOR_COUNT bail-out path
    # anyway). Three separate is_river-classified networks, one per mouth_type:
    #   Network A -- a Y-confluence: 0 and 1 (two headwaters) both flow into 2, which flows
    #   into 3 (the mouth, highest fallback discharge in the group), which flows onward into ocean
    #   node 4. mouth_type "ocean".
    #   Network B -- a single lone node 5, flowing into node 7, a real lake (lake_depth above
    #   the visible threshold). Node 7 itself is deliberately *not* is_river -- exactly the real
    #   invariant compute_hydrology enforces (a lake node is never classified is_river, see
    #   hydrology.compute_hydrology's own `~is_lake` term), unlike an earlier version of this
    #   fixture that fabricated is_river=True directly on the flooded node -- a state
    #   compute_hydrology can never actually produce, which had been masking a real group_rivers
    #   bug (mouth_type read the mouth's own, structurally-always-dry lake_depth instead of the
    #   target it flows into). mouth_type "lake".
    #   Network C -- a single lone node 8, a genuine unresolved dry sink (flow_target -1, no
    #   lake there yet). mouth_type "other".
    #   Network D -- a single lone node 9, ending at dry node 10 (still-dry basin floor, well
    #   above the lake's current shoreline) which itself flows onward one more hop into the
    #   same lake node 7 -- reproducing a real basin sitting near its own rim, with only its
    #   lowest node or two actually underwater: mouth_target alone reads as neither ocean nor
    #   lake, but the chain still reaches real lake water. mouth_type "lake" (confirmed
    #   directly against a real save -- see group_rivers's own docstring).
    # Node 6 is ordinary land with flow but deliberately *not* in the is_river mask, to check
    # grouping doesn't pull in non-river neighbors.
    points = np.array(
        [
            _normalize([1.0, 0.05, 0.0]),  # 0: headwater
            _normalize([1.0, -0.05, 0.02]),  # 1: headwater
            _normalize([1.0, 0.0, 0.03]),  # 2: confluence
            _normalize([1.0, 0.0, 0.06]),  # 3: mouth
            _normalize([1.0, 0.0, 0.09]),  # 4: ocean
            _normalize([-1.0, 0.0, 0.0]),  # 5: lone river node, ends at a lake
            _normalize([0.0, 1.0, 0.0]),  # 6: non-river land
            _normalize([-1.0, 0.02, 0.0]),  # 7: the lake node 5 flows into (not itself a river)
            _normalize([0.0, -1.0, 0.0]),  # 8: lone river node, a still-dry interior sink
            _normalize([0.0, -1.0, 0.02]),  # 9: lone river node, mouth of network D
            _normalize([0.0, -1.0, 0.04]),  # 10: dry basin floor between node 9 and the lake
        ]
    )
    n = len(points)
    elevation = np.array([100.0, 110.0, 80.0, 50.0, -10.0, 60.0, 40.0, 55.0, 70.0, 65.0, 58.0])
    is_ocean = np.array([False, False, False, False, True, False, False, False, False, False, False])
    flow_target = np.array([2, 2, 3, 4, -1, 7, -1, -1, -1, 10, 7])
    flow_accum = np.array([10.0, 12.0, 25.0, 40.0, 0.0, 5.0, 0.0, 0.0, 3.0, 4.0, 0.0])
    is_river = np.array([True, True, True, True, False, True, False, False, True, True, False])
    lake_depth = np.zeros(n)
    lake_depth[7] = 5.0  # node 7 well above LAKE_MIN_VISIBLE_DEPTH_M
    zeros = np.zeros(n)
    return hydrology.HydrologyFields(
        points=points,
        elevation=elevation,
        is_ocean=is_ocean,
        neighbor_idx=np.zeros((n, 1), dtype=np.int64),
        flow_target=flow_target,
        flow_accum=flow_accum,
        water_deposited=zeros.copy(),
        filled_elevation=zeros.copy(),
        spill_target=np.full(n, -1, dtype=np.int64),
        is_river=is_river,
        lake_depth=lake_depth,
        glacier_depth=zeros.copy(),
        plates_in_order=[],
    )


def test_compute_flow_direction_sink_redirects_once_lake_reaches_spill_point():
    # Same 4-node setup as above. Node 0 has no lower neighbor (its only neighbor, node 1,
    # is higher) so it's a sink candidate -- unless its current water surface has already
    # reached the basin's true spill point (50, from the test above).
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    filled = np.array([50.0, 50.0, 20.0, -10.0])
    spill = np.array([2, 2, 3, -1])

    not_yet_full, not_yet_spilling = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.zeros(4), filled, spill, np.zeros(4)
    )
    assert not_yet_full[0] == -1  # water surface 30 < 50, still a plain sink
    assert not not_yet_spilling[0]

    now_full, now_spilling = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.array([25.0, 0.0, 0.0, 0.0]), filled, spill, np.zeros(4)
    )
    assert now_full[0] == 2  # water surface 30+25=55 >= 50 -- redirects to spill_target
    assert now_spilling[0]


def test_compute_flow_direction_a_wet_sink_spills_over_a_rim_at_its_own_floor():
    # Node 0's rim is notched down to its own centre (filled == elevation). Holding a little
    # water, its hypsometric water surface sits 0.7 m below that centre. It still spills: a
    # floor-level rim holds no water at all (#503, a breached pit that held rivers forever).
    elevation = np.array([30.0, 50.0, 20.0, -10.0])
    is_ocean = np.array([False, False, False, True])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    filled = np.array([30.0, 50.0, 20.0, -10.0])
    spill = np.array([2, 2, 3, -1])

    flow_target, should_spill = hydrology._compute_flow_direction(
        elevation, is_ocean, neighbor_idx, np.array([-0.7, 0.0, 0.0, 0.0]), filled, spill, np.zeros(4)
    )
    assert flow_target[0] == 2
    assert should_spill[0]


def test_flow_cycles_lists_each_cycle_once():
    flow_target = np.array([1, 2, 0, 2, -1, 6, 5, 4])
    cycles = sorted(sorted(c.tolist()) for c in hydrology.flow_cycles(flow_target))
    assert cycles == [[0, 1, 2], [5, 6]]
    assert hydrology.flow_cycles(np.array([-1, 0, 1])) == []


def test_route_downstream_accumulates_along_a_chain():
    # A straight 4-node downhill chain (0 -> 1 -> 2 -> ocean-adjacent 3), each node
    # contributing 1.0 of its own -- accumulated flux at each node should be its own
    # contribution plus everything upstream of it.
    elevation = np.array([40.0, 30.0, 20.0, -5.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 3, -1])
    source = np.array([1.0, 1.0, 1.0, 0.0])

    through_flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source)
    assert through_flux.tolist() == [1.0, 2.0, 3.0, 0.0]
    assert deposited[3] == 3.0  # everything reaching the ocean-adjacent node's target settles there
    assert deposited[0] == 0.0 and deposited[1] == 0.0 and deposited[2] == 0.0


def test_route_downstream_retain_fraction_deposits_partway():
    # Same chain, but node 1 retains half of whatever passes through it.
    elevation = np.array([40.0, 30.0, 20.0, -5.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 3, -1])
    source = np.array([1.0, 1.0, 1.0, 0.0])
    retain = np.array([0.0, 0.5, 0.0, 0.0])

    through_flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source, retain_fraction=retain)
    # Node 1 receives 1.0 (its own) + 1.0 (from node 0) = 2.0, retains half (1.0) locally,
    # passes 1.0 onward.
    assert deposited[1] == 1.0
    assert through_flux[1] == 1.0
    # Node 2 receives 1.0 (its own) + 1.0 (passed on from node 1) = 2.0, all continues on.
    assert through_flux[2] == 2.0
    assert deposited[3] == 2.0
    # Conservation: total deposited equals total source.
    assert np.isclose(deposited.sum(), source.sum())


def test_route_downstream_loss_fraction_evaporates_without_depositing():
    # Same chain, but node 1 evaporates half of whatever passes through it -- unlike
    # retain_fraction, that half must simply vanish (not show up in `deposited` anywhere),
    # since it left as atmospheric moisture, not as material staying in the channel.
    elevation = np.array([40.0, 30.0, 20.0, -5.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 3, -1])
    source = np.array([1.0, 1.0, 1.0, 0.0])
    loss = np.array([0.0, 0.5, 0.0, 0.0])

    through_flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source, loss_fraction=loss)
    # Node 1 receives 1.0 (its own) + 1.0 (from node 0) = 2.0, evaporates half (1.0), passes
    # 1.0 onward -- none of the evaporated half is deposited anywhere.
    assert through_flux[1] == 1.0
    assert deposited[1] == 0.0
    assert through_flux[2] == 2.0  # node 2's own 1.0 + the 1.0 that made it past node 1
    assert deposited[3] == 2.0
    # Conservation is broken on purpose here: 1.0 of the original 3.0 total source evaporated
    # away instead of ending up in `deposited`.
    assert np.isclose(deposited.sum(), source.sum() - 1.0)


def test_connected_ocean_mask_excludes_an_enclosed_interior_pit():
    # 40 nodes along an arc: a wide ocean (0-14), a thick land ridge (15-27), a small
    # sub-sea-level pocket (28-30), then more land. The ridge is far wider than the k-NN
    # neighbour count, so no below-sea-level edge ever bridges it.
    theta = 0.05 * np.arange(40)
    xyz = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0   # the world ocean
    elevation[28:31] = -60.0   # an enclosed interior depression

    # Real cells, not 1/40 of the sphere each: at that size the pit would be an ocean basin.
    mask = hydrology.connected_ocean_mask(xyz, elevation, sea_level_m=0.0, areas_m2=np.full(40, 4.0e9))

    assert mask[0:15].all()          # the ocean is ocean
    assert not mask[28:31].any()     # the enclosed pit is not
    assert not mask[15:28].any() and not mask[31:].any()  # dry land never is


def test_connected_ocean_mask_leaves_a_single_connected_ocean_untouched():
    theta = 0.05 * np.arange(20)
    xyz = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)
    elevation = np.full(20, 100.0)
    elevation[0:12] = -200.0  # one contiguous below-sea-level body

    mask = hydrology.connected_ocean_mask(xyz, elevation, sea_level_m=0.0)
    assert mask.tolist() == [True] * 12 + [False] * 8


def _enclosed_basin_arc(basin_nodes: int, cell_km2: float = 1.0e6):
    """60 nodes along an arc, each `cell_km2`: a 20-node ocean (0-19), a land ridge far wider
    than the k-NN neighbour count (20-34), then an enclosed basin of `basin_nodes` nodes at
    -3,000 m starting at node 35. Returns (xyz, elevation, areas_m2, basin slice)."""
    theta = 0.05 * np.arange(60)
    xyz = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)
    elevation = np.full(60, 500.0)
    elevation[0:20] = -4_000.0
    basin = slice(35, 35 + basin_nodes)
    elevation[basin] = -3_000.0
    return xyz, elevation, np.full(60, cell_km2 * 1e6), basin


def test_connected_ocean_mask_counts_an_ocean_sized_enclosed_basin_as_ocean():
    # #758: a body cut off from the main ocean is ocean when it is ocean-sized, by its area in
    # km^2: 9 million km^2 here. A small enclosed pit is still a lake (3 million km^2).
    xyz, elevation, areas, basin = _enclosed_basin_arc(9)
    mask = hydrology.connected_ocean_mask(xyz, elevation, 0.0, areas_m2=areas)
    assert mask[basin].all() and mask[0:20].all()

    xyz, elevation, areas, pit = _enclosed_basin_arc(3)
    mask = hydrology.connected_ocean_mask(xyz, elevation, 0.0, areas_m2=areas)
    assert not mask[pit].any() and mask[0:20].all()


def test_connected_ocean_mask_keeps_a_basin_near_the_threshold_in_its_class():
    # An 8.1 million km^2 basin crosses the threshold when sea level, or a few metres of uplift
    # on its rim cell, moves one cell's worth of area. Its class follows last step's ocean,
    # so it doesn't flip with each small elevation change (design principle 10).
    xyz, elevation, areas, basin = _enclosed_basin_arc(5, cell_km2=1.62e6)
    rim = basin.stop - 1

    def classify(rim_elevation_m, prior):
        z = elevation.copy()
        z[rim] = rim_elevation_m
        return hydrology.connected_ocean_mask(xyz, z, 0.0, areas_m2=areas, prior_ocean=prior)

    entered = classify(-5.0, None)
    assert entered[basin].all()
    # The rim cell rises 10 m above the sea: 6.48 million km^2 left, but it was ocean.
    kept = classify(5.0, entered)
    assert kept[basin.start : rim].all() and not kept[rim]
    assert classify(-5.0, kept)[basin].all()
    # Without that history the same 6.48 million km^2 basin is a lake, and it stays one.
    assert not classify(5.0, None)[basin].any()
    assert not classify(5.0, np.zeros(60, dtype=bool))[basin].any()
    # The hysteresis has a floor: a former ocean cut down to 3.24 million km^2 is a lake.
    z = elevation.copy()
    z[basin.start + 2 : basin.stop] = 5.0
    shrunk = hydrology.connected_ocean_mask(xyz, z, 0.0, areas_m2=areas, prior_ocean=entered)
    assert not shrunk[basin].any()


def test_connected_ocean_mask_one_prior_ocean_cell_does_not_latch_hysteresis():
    # #758 review: `prior_ocean` is a nearest-neighbour resample, so a stray cell can alias
    # across a one-cell ridge. One aliased cell must not drop a lake between the keep and
    # entry areas (6.5 million km^2 here) to the keep threshold -- latched, its own cells
    # would become next step's prior and hold it ocean. Most of the body having been ocean
    # does keep it.
    xyz, elevation, areas, basin = _enclosed_basin_arc(5, cell_km2=1.3e6)

    stray = np.zeros(60, dtype=bool)
    stray[basin.start] = True
    assert not hydrology.connected_ocean_mask(xyz, elevation, 0.0, areas_m2=areas, prior_ocean=stray)[basin].any()

    majority = np.zeros(60, dtype=bool)
    majority[basin.start : basin.start + 3] = True
    assert hydrology.connected_ocean_mask(xyz, elevation, 0.0, areas_m2=areas, prior_ocean=majority)[basin].all()


def test_connected_ocean_mask_keeps_a_sea_that_splits_off_the_main_ocean():
    # A 5 million km^2 sea joined to the main ocean across a shallow sill (cells 20-25) is part
    # of the largest component. When the sill rises 20 m above the sea, the sea stays ocean: it
    # held last step's ocean, so it doesn't become a dry lake that eustasy must refill when the
    # sill next drops (the flicker test_sea_level_and_land_carry_over_and_the_world_keeps_
    # stepping caught, with a 41-47 million km^2 shelf sea and the earlier crust rule).
    xyz, elevation, areas, _ = _enclosed_basin_arc(0)
    sea = slice(26, 31)
    elevation[sea] = -3_000.0

    def classify(sill_elevation_m, prior):
        z = elevation.copy()
        z[20:26] = sill_elevation_m
        return hydrology.connected_ocean_mask(xyz, z, 0.0, areas_m2=areas, prior_ocean=prior)

    joined = classify(-10.0, None)
    assert joined[sea].all()
    assert classify(10.0, joined)[sea].all()
    assert not classify(10.0, None)[sea].any()


def _enclosed_ocean_basin_world(lake_depth_m: float, lake_ice_m: float):
    """The 40-node interior-pit row below, with its 3-node basin (28-30) at -3,000 m holding
    `lake_depth_m` of water under `lake_ice_m` of ice."""
    elevation = np.full(40, 300.0)
    elevation[0:15] = -4_000.0
    elevation[28:31] = -3_000.0
    depth = np.zeros(40)
    depth[28:31] = lake_depth_m
    ice = np.zeros(40)
    ice[28:31] = lake_ice_m
    plate = quad_plate(0, "continental", columns=_row(40), elevation=elevation, lake_depth=depth, lake_ice_m=ice)
    return World(seed=0, plates=[plate]), plate


def test_reclassified_basin_hands_its_lake_water_to_the_ocean_and_closes_the_ledgers(monkeypatch):
    # #758: the moment an enclosed basin becomes ocean, the lake water and lake ice it holds
    # leave the lake reservoirs for the ocean, booked in the lake and land ledgers, and
    # eustasy fills it to sea level from the conserved total budget. The 3-node basin is about
    # 51,000 km^2, so the thresholds are scaled to this fixture once it has stepped as a lake.
    world, plate = _enclosed_ocean_basin_world(lake_depth_m=1_000.0, lake_ice_m=1.5)
    areas = _node_areas(world)
    basin = slice(28, 31)
    # 5 C keeps some lake ice through the step; colder, sublimation in the solve strips it.
    precipitation, temperature, years = np.full(40, 800.0), np.full(40, 5.0), 10_000.0

    def step():
        fields = hydrology.compute_hydrology(world, precipitation, temperature, years=years)
        world.hydrology_cache = fields
        plate.set_fields_on_plate(
            lake_depth=fields.lake_depth, lake_ice_m=fields.lake_ice_m, river_water_m=fields.river_water_m
        )
        return fields

    as_lake = step()
    assert not as_lake.is_ocean[basin].any()
    eustasy.initialize_water_budget(world)
    total_m3 = world.ocean_water_volume_m3
    liquid_m3 = float(np.dot(areas[basin], as_lake.lake_depth[basin]))
    ice_m3 = float(np.dot(areas[basin], as_lake.lake_ice_m[basin]))
    assert liquid_m3 > 0.0 and ice_m3 > 0.0

    monkeypatch.setattr(hydrology, "ENCLOSED_OCEAN_MIN_AREA_KM2", 40_000.0)
    monkeypatch.setattr(hydrology, "ENCLOSED_OCEAN_KEEP_AREA_KM2", 20_000.0)
    fields = step()

    assert fields.is_ocean[basin].all()
    assert not fields.lake_depth[basin].any() and not fields.lake_ice_m[basin].any()
    lake, land = fields.lake_water_ledger, fields.land_water_ledger
    assert lake.melting_m3 >= ice_m3 * (1 - 1e-12)
    assert lake.overflow_out_m3 >= (liquid_m3 + ice_m3) * (1 - 1e-12)
    assert abs(lake.residual_m3) <= 1e-9 * liquid_m3
    assert abs(land.residual_m3) <= 1e-9 * liquid_m3

    # Eustasy: the basin's water is no longer trapped, the ocean fills it to sea level, and
    # the total surface-water budget is unchanged.
    eustasy.update_sea_level(world)
    assert world.ocean_water_volume_m3 == total_m3
    ocean_m3 = eustasy.water_volume_for_sea_level(world, world.sea_level_m)
    solve_tolerance_m3 = eustasy._SOLVE_TOLERANCE_M * float(areas.mean())
    assert ocean_m3 + eustasy.trapped_water_volume_m3(world) == pytest.approx(total_m3, abs=solve_tolerance_m3)
    assert eustasy._ocean_connected_mask(world, plate.collect("elevation"), world.sea_level_m)[basin].all()

    # On the next step it stays ocean, and the ledgers still close.
    again = step()
    assert again.is_ocean[basin].all()
    assert abs(again.lake_water_ledger.residual_m3) <= 1e-9 * liquid_m3
    assert abs(again.land_water_ledger.residual_m3) <= 1e-9 * liquid_m3


def test_reclassified_basin_gets_ocean_currents_on_the_next_climate_step(monkeypatch):
    # #758 review: once hydrology classes the basin as ocean, climate must give it real
    # currents, not just the `is_ocean` flag. Climate runs before hydrology in a world step,
    # so the currents arrive on the step after reclassification, read off the new
    # `hydrology_cache` (`climate._sample_elevation_and_crust` -> `sample_is_ocean` ->
    # `compute_ocean_currents`, which has no connectivity gate of its own).
    from app.hydroclimate import climate

    world, plate = _enclosed_ocean_basin_world(lake_depth_m=1_000.0, lake_ice_m=1.5)
    basin = slice(28, 31)
    precipitation, temperature, years = np.full(40, 800.0), np.full(40, 5.0), 10_000.0

    def step():
        fields = hydrology.compute_hydrology(world, precipitation, temperature, years=years)
        world.hydrology_cache = fields
        plate.set_fields_on_plate(
            lake_depth=fields.lake_depth, lake_ice_m=fields.lake_ice_m, river_water_m=fields.river_water_m
        )
        return fields

    def basin_cells(fields):
        cells = climate.compute_climate(world)
        points = fields.points[basin]
        lat = np.degrees(np.arcsin(points[:, 2] / np.linalg.norm(points, axis=1)))
        lon = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
        rows = np.argmin(np.abs(cells.lat_deg[:, None] - lat[None, :]), axis=0)
        cols = np.argmin(np.abs(cells.lon_deg[:, None] - lon[None, :]), axis=0)
        return cells, rows, cols

    as_lake = step()
    cells, rows, cols = basin_cells(as_lake)
    assert not cells.is_ocean[rows, cols].any()
    assert not np.hypot(cells.current_u[rows, cols], cells.current_v[rows, cols]).any()

    monkeypatch.setattr(hydrology, "ENCLOSED_OCEAN_MIN_AREA_KM2", 40_000.0)
    monkeypatch.setattr(hydrology, "ENCLOSED_OCEAN_KEEP_AREA_KM2", 20_000.0)
    as_ocean = step()
    assert as_ocean.is_ocean[basin].all()
    cells, rows, cols = basin_cells(as_ocean)
    assert cells.is_ocean[rows, cols].all()
    assert np.hypot(cells.current_u[rows, cols], cells.current_v[rows, cols]).max() > 0.0


def test_compute_hydrology_treats_an_interior_pit_as_an_endorheic_basin_that_silts_in():
    d = 0.05
    theta = d * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    fields = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=1_000_000)

    assert fields.is_ocean[0:15].all()
    assert not fields.is_ocean[28:31].any()  # classified land, not a vanishing ocean sink
    # Rain routed into the closed pit pools there and drops sediment onto the wet floor.
    assert fields.lake_depth[28:31].max() > 0.0
    assert fields.silt_deposited[28:31].max() > 0.0


def test_compute_hydrology_breaches_a_shallow_pit_instead_of_ponding_it(monkeypatch):
    # A 20 m deep pit on a 100 m plateau, in a wet climate (2,000 mm/yr at 10 C), so its
    # runoff outpaces a rim-level lake's evaporation. One 100 kyr step carves 30 m of ordinary
    # rock, so with breaching on the rim is notched and the pit drains. With it off, rain
    # ponds there.
    d = 0.05
    theta = d * np.arange(40)
    elevation = np.full(40, 100.0)
    elevation[0:15] = -100.0
    elevation[28:31] = 80.0

    def run(precipitation=2000.0, temperature=10.0, lake_depth=None):
        plate = (
            _flow_line_plate(0, theta, elevation)
            if lake_depth is None
            else _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
        )
        world = World(seed=0, plates=[plate])
        return hydrology.compute_hydrology(world, np.full(40, precipitation), np.full(40, temperature), years=100_000)

    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", False)
    closed = run()
    assert closed.lake_depth[28:31].max() > 0.0
    assert closed.breach_notch_m.shape == (0,)

    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", True)
    breached = run()
    # The notch is cut to the pit's floor and leads lower, so even the low ground inside the
    # pit's cells drains through it (lakes.Lake.drains_at_floor).
    assert breached.lake_depth[28:31].max() == 0.0
    # Closed, the same water keeps rising toward the 100 m rim on the next step, once the
    # pit's cells have merged into one lake.
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", False)
    closed_again = run(lake_depth=closed.lake_depth)
    assert closed_again.lake_depth[28:31].min() > 10.0
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", True)
    assert breached.breach_notch_m.max() == pytest.approx(20.0)
    assert np.all(breached.filled_elevation[28:31] <= 80.0)

    # The same pit in a semi-arid climate (800 mm/yr at 15 C) would level off below its rim,
    # so it stays a closed endorheic basin.
    arid = run(precipitation=800.0, temperature=15.0)
    assert arid.lake_depth[28:31].max() > 0.0
    assert not arid.breach_notch_m.any()


def _plateau_pit_run(precipitation_mm: float, temperature_c: float, years: float = 0.0, **kw):
    """The 20 m pit on a 100 m plateau above, by default in a zero-duration pass like
    generation's river bootstrap (#218), with `kw` passed on to compute_hydrology."""
    theta = 0.05 * np.arange(40)
    elevation = np.full(40, 100.0)
    elevation[0:15] = -100.0
    elevation[28:31] = 80.0
    world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
    return hydrology.compute_hydrology(world, np.full(40, precipitation_mm), np.full(40, temperature_c), years=years, **kw)


def _reaches_ocean(fields, node: int) -> bool:
    for _ in range(len(fields.flow_target)):
        if fields.is_ocean[node]:
            return True
        if fields.flow_target[node] < 0:
            return False
        node = int(fields.flow_target[node])
    return False


def test_compute_hydrology_history_years_gives_a_zero_duration_pass_a_breach_budget(monkeypatch):
    # With no time to carve, the wet pit ends its river; `history_years` gives the breach
    # budget a span of its own.
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", True)
    instant = _plateau_pit_run(2000.0, 10.0)
    assert not instant.breach_notch_m.any()
    assert not _reaches_ocean(instant, 29)

    history = _plateau_pit_run(2000.0, 10.0, history_years=100_000)
    assert history.breach_notch_m.max() == pytest.approx(20.0)
    assert history.lake_depth[28:31].max() == 0.0
    assert _reaches_ocean(history, 29)


def test_compute_hydrology_history_years_does_not_advance_lake_storage(monkeypatch):
    # D1b integrates lake storage over `years`; `history_years` supplies a carving history only.
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", False)
    instant = _plateau_pit_run(2000.0, 10.0)
    assert not _reaches_ocean(instant, 29)

    history = _plateau_pit_run(2000.0, 10.0, history_years=100_000)
    assert not _reaches_ocean(history, 29)
    assert not history.lake_depth.any()
    assert history.discharge_m3_s[29] == 0.0


@pytest.mark.usefixtures("march_glaciers")
def test_compute_hydrology_history_years_leaves_ice_on_the_zero_length_step(monkeypatch):
    # Cold and snowy: snow makes ice over the step's own `years` only, so a
    # zero-duration pass grows none, whatever span basins get to fill.
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", True)
    instant = _plateau_pit_run(1000.0, -10.0)
    history = _plateau_pit_run(1000.0, -10.0, history_years=1_000_000)
    np.testing.assert_array_equal(history.glacier_depth, instant.glacier_depth)
    assert not history.glacier_depth.any()
    # Over a real step the same climate does grow ice, so the check above has teeth.
    stepped = _plateau_pit_run(1000.0, -10.0, years=100_000)
    assert stepped.glacier_depth.any()


def test_compute_hydrology_populates_is_sea_once_a_basin_stays_flooded():
    # Same interior-pit fixture as the test above. is_sea is derived from *last* step's own
    # flooded extent (lakes._classify_tier), so it can't promote on the very first pass -- but
    # this tiny synthetic world's implied node_area_km2 (whole-sphere area / 40 nodes) is
    # enormous, so a single wet node held over from a prior step already clears
    # SEA_MIN_FLOODED_AREA_KM2 on the next one.
    d = 0.05
    theta = d * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    first = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=1_000_000)
    assert first.lake_depth[28:31].max() > 0.0
    assert first.is_sea.shape == first.lake_depth.shape
    assert not first.is_sea.any()  # nothing was wet *last* step yet -- can't have promoted

    second_plate = _flow_line_plate_with_lake(0, theta, elevation, first.lake_depth)
    second_world = World(seed=0, plates=[second_plate])
    second = hydrology.compute_hydrology(second_world, np.full(40, 800.0), np.full(40, 15.0), years=1_000_000)
    assert second.lake_depth[28:31].max() > 0.0
    assert second.is_sea[28:31].any()
    assert not second.is_sea[0:15].any()  # ocean itself is never "sea tier" -- is_sea is lake-only


def test_compute_hydrology_end_to_end_on_a_small_synthetic_world():
    # 12 points along a monotonically descending line, the last three underwater --
    # exercises the whole pipeline (basin-spill, flow direction, accumulation, river
    # classification) on a real (if small) node cloud built the same way as the rest of the
    # test suite's synthetic-plate fixtures.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)  # last three points go underwater
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    precipitation_at_nodes = np.full(12, 800.0)
    temperature_at_nodes = np.full(12, 15.0)  # well above freezing -- no glaciation involved
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes, temperature_at_nodes, years=1)

    assert len(fields.points) == 12
    assert fields.is_ocean.tolist() == [False] * 9 + [True] * 3
    # Every land node should have a full path to the ocean on this simple monotonic line.
    assert np.all(np.isfinite(fields.filled_elevation[~fields.is_ocean]))
    # flow_accum should be non-decreasing heading downhill along this single chain (each
    # node picks up whatever accumulated above it, net of nothing lost here).
    land_accum = fields.flow_accum[~fields.is_ocean]
    assert np.all(np.diff(land_accum) >= -1e-6)


def test_lake_component_sizes_groups_connected_lake_nodes_and_ignores_dry_ones():
    # 0-1-2 chain all flooded, node 3 (also flooded) reachable only from node 1 via a
    # one-directional edge (1 lists 3 as a neighbor, 3 does not list 1 back) -- confirms
    # sizing treats neighbor_idx as undirected, not just whatever direction happens to be
    # listed. Node 4 is dry and must read 0 regardless of being graph-adjacent to lake nodes.
    is_lake = np.array([True, True, True, True, False])
    neighbor_idx = np.array([[1, 1], [2, 3], [1, 1], [1, 1], [2, 2]])

    sizes = hydrology._lake_component_sizes(is_lake, neighbor_idx)
    assert sizes.tolist() == [4, 4, 4, 4, 0]


def test_lake_component_sizes_splits_disconnected_lakes():
    # Two separate two-node ponds (0-1 and 2-3), never connected to each other.
    is_lake = np.array([True, True, True, True])
    neighbor_idx = np.array([[1], [0], [3], [2]])

    sizes = hydrology._lake_component_sizes(is_lake, neighbor_idx)
    assert sizes.tolist() == [2, 2, 2, 2]


def test_compute_hydrology_spilling_lake_erodes_its_outlet_from_its_own_surface_area():
    # A 22-node line: an upland ramp (0-9) feeding a sink (10) sitting behind a barrier
    # (11-13, peaking at 300) that descends to the ocean (19-21). filled_elevation[10] works
    # out to 200 (the true minimal-bottleneck crossing, confirmed directly against this exact
    # layout) -- lake_depth[10]=160 puts the sink's water surface at 50+160=210, over that
    # cap, so it's actively spilling this step.
    d = 0.002
    theta = d * np.arange(22)
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    lake_depth = np.zeros(22)
    lake_depth[10] = 160.0

    plate = _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
    world = World(seed=0, plates=[plate])

    precipitation_at_nodes = np.full(22, 10.0)  # small -- the breach term should dominate
    temperature_at_nodes = np.full(22, 15.0)
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes, temperature_at_nodes, years=1)

    # A single-node lake (only node 10 is above LAKE_MIN_VISIBLE_DEPTH_M), so the breach term
    # there is exactly one component-size's worth of LAKE_BREACH_EROSION_COEFFICIENT --
    # dwarfing the handful of small precip contributions reaching the same node.
    assert fields.flow_accum[10] >= hydrology.LAKE_BREACH_EROSION_COEFFICIENT
    assert fields.flow_accum[10] < hydrology.LAKE_BREACH_EROSION_COEFFICIENT * 1.1


def test_compute_hydrology_spilling_lake_also_erodes_its_own_rim_when_that_is_the_higher_side():
    # The multi-node counterpart to the test above: a two-node basin (floor 10, shore/rim 11)
    # behind a barrier, engineered so the basin's *own* rim member (11, elevation 210) is the
    # higher side of the governing boundary edge -- higher than the far side outside the lake
    # (15, elevation 200) that `max_depth` would otherwise appear to come from (confirmed
    # directly: `lakes.build_lake_hierarchy` on this exact layout gives max_depth=210 via
    # outlet_node_idx=11, not 200 via outlet_target_idx=15). The one-year fixture starts with
    # enough stored water to export surplus through this rim and exercise both erosion sites.
    d = 0.002
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0, 210.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0, -150.0, -200.0, -250.0, -300.0]
    )
    n = len(elevation)
    theta = d * np.arange(n)
    lake_depth = np.zeros(n)
    lake_depth[10] = 300.0
    plate = _flow_line_plate_with_lake(0, theta, elevation, lake_depth)
    world = World(seed=0, plates=[plate])

    precipitation_at_nodes = np.full(n, 10.0)  # small -- the breach term should dominate
    temperature_at_nodes = np.full(n, 15.0)
    fields = hydrology.compute_hydrology(world, precipitation_at_nodes, temperature_at_nodes, years=1)

    # Both nodes of the two-node basin (floor 10 and partly flooded rim member 11) exceed
    # LAKE_MIN_VISIBLE_DEPTH_M after the prognostic solve, so both sides of the bottleneck
    # carry two nodes' worth of LAKE_BREACH_EROSION_COEFFICIENT (#506).
    assert fields.flow_accum[11] >= 2.0 * hydrology.LAKE_BREACH_EROSION_COEFFICIENT
    assert fields.flow_accum[11] < 2.1 * hydrology.LAKE_BREACH_EROSION_COEFFICIENT
    assert fields.flow_accum[15] >= 2.0 * hydrology.LAKE_BREACH_EROSION_COEFFICIENT
    lake = fields.lake_forest[0]
    assert fields.should_spill[lake.sink_node_idx]
    assert lake.is_spilling
    assert fields.flow_target[lake.sink_node_idx] == lake.outlet_target_idx
    assert fields.discharge_m3_s[lake.outlet_target_idx] > 0.0
    from app import main

    inspector = main._lake_basin_summary(fields, lake, 0, [], {}, is_lake=True)
    assert inspector["is_spilling"]


def test_breach_surge_drives_flow_accum_but_never_physical_discharge(monkeypatch):
    # The single-node spilling lake above. Its surge rides the legacy flow_accum (erosion's
    # input) but is kept out of discharge_m3_s (issue #355): turning the surge off leaves the
    # discharge exactly unchanged.
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    n = len(elevation)
    lake_depth = np.zeros(n)
    lake_depth[10] = 160.0

    def run():
        world = World(seed=0, plates=[_flow_line_plate_with_lake(0, 0.002 * np.arange(n), elevation, lake_depth)])
        return hydrology.compute_hydrology(world, np.full(n, 10.0), np.full(n, 15.0), years=1)

    surged = run()
    assert surged.breach_surge[10] == pytest.approx(hydrology.LAKE_BREACH_EROSION_COEFFICIENT)
    monkeypatch.setattr(hydrology, "LAKE_BREACH_EROSION_COEFFICIENT", 0.0)
    calm = run()
    assert np.all(calm.breach_surge == 0.0)
    assert surged.flow_accum[10] > calm.flow_accum[10] + 0.9 * surged.breach_surge[10]
    assert np.array_equal(surged.discharge_m3_s, calm.discharge_m3_s)
    # External supply now enters storage; only surplus leaves at the hierarchy outlet (#588).
    assert surged.discharge_m3_s[10] == 0.0
    assert surged.lake_water_ledger.residual_m3 == pytest.approx(0.0, abs=1.0e6)


def test_breach_surge_never_sizes_lake_silt(monkeypatch):
    # A spilling lake (sink 10, behind the 11-13 barrier) whose overflow runs into a pit (16)
    # before the ocean. The surge settles in the pit's legacy water_deposited, but silt is
    # sized from the physical inflow (issue #363), so turning the surge off leaves every
    # lake's silt unchanged. Breaching is off so the pit holds water instead of draining.
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0]
        + [40.0]
        + [120.0, 130.0, 140.0, 130.0, 60.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    n = len(elevation)
    lake_depth = np.zeros(n)
    lake_depth[10] = 160.0
    lake_depth[16] = 5.0
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", False)

    def run():
        world = World(seed=0, plates=[_flow_line_plate_with_lake(0, 0.002 * np.arange(n), elevation, lake_depth)])
        return hydrology.compute_hydrology(world, np.full(n, 2000.0), np.full(n, 5.0), years=1)

    surged = run()
    assert surged.breach_surge[10] > 0.0
    monkeypatch.setattr(hydrology, "LAKE_BREACH_EROSION_COEFFICIENT", 0.0)
    calm = run()
    assert surged.water_deposited[16] > calm.water_deposited[16]  # the surge reached the pit
    assert surged.silt_deposited[16] > 0.0 and surged.silt_deposited[10] > 0.0
    assert np.array_equal(surged.silt_deposited, calm.silt_deposited)


def test_a_lake_sink_stores_external_supply_instead_of_passing_full_inflow():
    # The single-node lake above. External sources reach its storage, while only overflow leaves
    # the lake and enters downstream discharge (issue #588).
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    n = len(elevation)
    lake_depth = np.zeros(n)
    lake_depth[10] = 160.0
    world = World(seed=0, plates=[_flow_line_plate_with_lake(0, 0.002 * np.arange(n), elevation, lake_depth)])
    fields = hydrology.compute_hydrology(world, np.full(n, 10.0), np.full(n, 15.0), years=1)

    assert fields.breach_surge[10] > 0.0
    inflow_m3_s = float(fields.discharge_m3_s[fields.flow_target == 10].sum())
    assert inflow_m3_s > 0.0
    # The sink's external supply is now balanced into storage, evaporation and overflow;
    # discharge at the lake itself contains no full-inflow pass-on (#588).
    assert fields.discharge_m3_s[10] == 0.0
    assert fields.lake_water_ledger.residual_m3 == pytest.approx(0.0, abs=1.0e6)


def test_rivers_flow_on_through_a_cold_band_to_the_sea():
    # Warm uplands drain through a band below 0 C on the annual mean to the sea. The band's
    # nodes used to be dead ends where the water froze at WATER_TO_ICE, 1/50,000 of its
    # volume, so nearly all of it was deleted (#355, #519). Now it flows on and reaches the
    # ocean less only the transit loss (RFC 0014 D2(b), #539).
    theta = 0.002 * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    temperature = np.where(np.arange(12) < 6, 15.0, -10.0)
    world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
    fields = hydrology.compute_hydrology(world, np.full(12, 800.0), temperature, years=100_000)
    cold_land = (temperature < 0.0) & ~fields.is_ocean
    assert np.all(fields.flow_target[cold_land] >= 0)
    upland_m3_s = float(fields.discharge_m3_s[5])
    keep = np.prod(1.0 - hydrology.river_transit_loss_fraction(temperature[cold_land], np.zeros(int(cold_land.sum()), dtype=bool)))
    assert upland_m3_s > 0.0
    assert float(fields.discharge_m3_s[fields.is_ocean].sum()) >= upland_m3_s * keep * (1.0 - 1e-9)


@pytest.mark.usefixtures("march_glaciers")
def test_discharge_carries_the_ice_ledgers_meltwater_exactly_once_at_any_step_length():
    # An all-land chain under 200 m of melting ice, draining to its lowest node. The discharge
    # arriving at that sink is the ledger's meltwater (surface and basal melt, overflow, capped ice, bare-ground snow; not
    # calving) spread over the step, in m^3/s (issue #355).
    theta = 0.002 * np.arange(12)
    elevation = 900.0 - 60.0 * np.arange(12)
    for years in (100_000, 1_000_000):
        world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, np.zeros(12), glacier_depth=np.full(12, 200.0))])
        fields = hydrology.compute_hydrology(
            world, np.full(12, 100.0), np.full(12, -6.0), years=years, seasonal_amplitude_at_nodes=np.full(12, 18.0)
        )
        ledger = fields.ice.ledger
        sink = fields.flow_target == -1
        assert not np.any(fields.is_ocean) and np.count_nonzero(sink) == 1
        assert ledger.surface_melt_m3 > 0.0 and ledger.basal_melt_m3 > 0.0
        assert fields.lake_external_meltwater_m3_per_yr.sum() > 0.0
        assert np.all(fields.discharge_m3_s >= 0.0)
        # Physical fluxes are ~70x the old model-unit ones, so the lake's absolute
        # residual scales with its throughput (RFC 0015 D1, #789).
        supply_m3 = (
            ledger.meltwater_m3 * lithosphere.RHO_ICE / lithosphere.RHO_WATER
            + fields.lake_external_dry_runoff_m3_per_yr.sum() * years
        )
        assert fields.lake_water_ledger.residual_m3 == pytest.approx(0.0, abs=1e-9 * supply_m3)
        # A reach starts at `Q tau`, so what the river held back of this melt is matched by
        # what its inventory released (#725). The match is to a part per million: the march
        # (RFC 0015 D2, #789) sums hundreds of substep ledgers, and the initial-vs-final
        # reach inventory differs by the transit transient, so single-solve 1e-9 no longer
        # holds.
        # Meltwater converts by mass (RFC 0015 D1, #789): the routed water is the ice
        # ledger's meltwater times RHO_ICE / RHO_WATER.
        assert (
            fields.lake_external_meltwater_m3_per_yr.sum()
            + fields.lake_external_river_storage_m3_per_yr.sum()
        ) == pytest.approx(
            ledger.meltwater_m3 * lithosphere.RHO_ICE / lithosphere.RHO_WATER / years, rel=1e-6
        )
        assert fields.lake_external_wet_precip_m3_per_yr.sum() == 0.0
        # Rain on the ice (a third of the precipitation at -6 C with an 18 C swing) drains
        # as ordinary runoff for now; D8 will route it on the ice instead.
        assert fields.lake_external_dry_runoff_m3_per_yr.sum() > 0.0


def test_discharge_from_fixed_rainfall_does_not_depend_on_step_length():
    # A warm, ice-free chain to the sea: runoff and the transit loss are both per-year rates,
    # so the same rainfall gives the same discharge at a 100 kyr and a 1 Myr step (issue #355).
    # The legacy flow_accum still loses more water to evaporation on the longer step.
    theta = 0.002 * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)

    def run(years):
        world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
        return hydrology.compute_hydrology(world, np.full(12, 800.0), np.full(12, 20.0), years=years)

    short, long = run(100_000), run(1_000_000)
    land = ~short.is_ocean
    assert np.max(short.discharge_m3_s[land]) > 0.0
    # Rivers hold days of a 100 kyr interval, so only their tiny transient differs (#725).
    assert long.discharge_m3_s == pytest.approx(short.discharge_m3_s, rel=1e-7)
    assert np.max(long.flow_accum[land]) < np.max(short.flow_accum[land])


def test_fixed_catchment_has_the_same_outlet_discharge_at_two_node_densities():
    # Split the same 10,000 km^2 rainy catchment into 4 or 64 accounting cells, the
    # fourfold linear-density change from density 1 to 4 (issue #355 step 5). Routing
    # area-weighted runoff must give the same outlet discharge and river forcing.
    catchment_area_m2 = 1e10
    runoff_m = breaching.runoff_mm(np.array([800.0]), breaching.open_water_evaporation_mm(np.array([15.0])))[0] / 1000.0
    outflows = []
    for density in (1, 4):
        source_count = 4 * density**2
        outlet, ocean = source_count, source_count + 1
        elevation = np.concatenate((np.full(source_count, 2.0), [1.0, 0.0]))
        is_ocean = np.arange(source_count + 2) == ocean
        target = np.concatenate((np.full(source_count, outlet), [ocean, -1]))
        source = np.concatenate((np.full(source_count, runoff_m * catchment_area_m2 / source_count), [0.0, 0.0]))
        through, settled = hydrology.route_downstream(elevation, is_ocean, target, source)
        assert settled[ocean] == pytest.approx(runoff_m * catchment_area_m2)
        outflows.append(through[outlet] / hydrology.SECONDS_PER_YEAR)

    assert outflows[0] == pytest.approx(outflows[1], rel=1e-12)
    assert outflows[0] == pytest.approx(runoff_m * catchment_area_m2 / hydrology.SECONDS_PER_YEAR)
    assert hydrology.discharge_m3_s_to_flow_accum_mm_yr(np.array(outflows))[0] == pytest.approx(
        hydrology.discharge_m3_s_to_flow_accum_mm_yr(np.array(outflows))[1]
    )
    assert hydrology.compute_river_speed(np.array([0.01, 0.01]), np.array(outflows))[0] == pytest.approx(
        hydrology.compute_river_speed(np.array([0.01, 0.01]), np.array(outflows))[1]
    )
    river = erosion._river_erosion_amount(
        np.array(outflows), np.zeros(2), np.array([0.01, 0.01]), np.zeros(2), 1.0, 1.0,
    )
    width = erosion._channel_width_growth(np.array(outflows), 1.0)
    assert river[0] > 0.0 and river[0] == pytest.approx(river[1])
    assert width[0] > 0.0 and width[0] == pytest.approx(width[1])


def test_compute_flow_direction_prefers_an_existing_channel_over_the_literal_steepest_neighbor():
    # Node 0's two downhill candidates: node 1 (elevation 20, steeper, no channel) and node 2
    # (elevation 25, less steep, but already has a real established channel).
    elevation = np.array([30.0, 20.0, 25.0])
    is_ocean = np.zeros(3, dtype=bool)
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])
    filled = np.array([30.0, 20.0, 25.0])
    spill = np.full(3, -1)

    steepest, _ = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, np.zeros(3))
    assert steepest[0] == 1  # no channel anywhere -- plain steepest descent

    channel_depth = np.array([0.0, 0.0, 50.0])  # node 2's channel is well-established
    channelized, _ = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, channel_depth)
    assert channelized[0] == 2  # prefers the channelized neighbor even though it's less steep


def test_compute_flow_direction_ignores_a_barely_established_channel():
    # Same setup, but node 2's channel is below CHANNEL_PREFERENCE_THRESHOLD_M -- too new/
    # shallow to count as "an existing channel," so plain steepest descent still wins.
    elevation = np.array([30.0, 20.0, 25.0])
    is_ocean = np.zeros(3, dtype=bool)
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])
    filled = np.array([30.0, 20.0, 25.0])
    spill = np.full(3, -1)
    channel_depth = np.array([0.0, 0.0, hydrology.CHANNEL_PREFERENCE_THRESHOLD_M - 1.0])

    result, _ = hydrology._compute_flow_direction(elevation, is_ocean, neighbor_idx, np.zeros(3), filled, spill, channel_depth)
    assert result[0] == 1


def _node_areas(world):
    return collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))


def test_seasonal_positive_degrees_matches_a_sampled_annual_cycle():
    mean = np.array([-20.0, -12.0, -5.0, 0.0, 3.0, 12.0, 25.0])
    amplitude = np.array([10.0, 15.0, 15.0, 8.0, 0.0, 10.0, 5.0])
    t = np.linspace(0.0, 2.0 * np.pi, 200_001)[:-1]
    sampled = np.mean(np.clip(mean[:, None] + amplitude[:, None] * np.cos(t)[None, :], 0.0, None), axis=1)
    assert hydrology.seasonal_positive_degrees(mean, amplitude) == pytest.approx(sampled, abs=1e-4)


def test_seasonal_positive_degrees_limits():
    # No swing: degrees above freezing of the annual mean itself. A swing that never reaches
    # the threshold: zero. A cycle centred on it: amplitude / pi. Fully above: the mean.
    assert np.array_equal(hydrology.seasonal_positive_degrees(np.array([-3.0, 0.0, 4.0]), np.zeros(3)), [0.0, 0.0, 4.0])
    assert hydrology.seasonal_positive_degrees(np.array([-15.0]), np.array([10.0]))[0] == 0.0
    assert hydrology.seasonal_positive_degrees(np.array([0.0]), np.array([9.0]))[0] == pytest.approx(9.0 / np.pi)
    assert hydrology.seasonal_positive_degrees(np.array([12.0]), np.array([10.0]))[0] == pytest.approx(12.0)
    # Increasing in both the mean and (once the swing reaches the threshold) the amplitude.
    means = np.linspace(-20.0, 5.0, 101)
    assert np.all(np.diff(hydrology.seasonal_positive_degrees(means, np.full_like(means, 12.0))[33:]) > 0.0)
    amps = np.linspace(0.0, 30.0, 101)
    assert np.all(np.diff(hydrology.seasonal_positive_degrees(np.full_like(amps, -6.0), amps)[21:]) > 0.0)


def test_seasonal_lake_ice_helpers_match_analytic_cycle_and_are_continuous():
    means = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    amplitude = np.full_like(means, 2.0)
    fraction = hydrology.seasonal_fraction_above(means, amplitude)
    assert fraction == pytest.approx([0.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 1.0])
    assert hydrology.seasonal_open_mean_c(np.array([-1.0]), np.array([2.0]))[0] == pytest.approx(
        -1.0 + 2.0 * np.sin(np.pi / 3.0) / (np.pi / 3.0)
    )
    assert hydrology.seasonal_freezing_degrees(np.array([-1.0]), np.array([2.0]))[0] == pytest.approx(
        2.0 / 3.0 + np.sqrt(3.0) / np.pi
    )
    # Zero amplitude and threshold-touching limits are defined; the open-season mean tends
    # to the threshold as the ice-free fraction tends to zero.
    assert hydrology.seasonal_fraction_above(np.array([-1.0, 0.0, 1.0]), np.zeros(3)) == pytest.approx(
        [0.0, 0.5, 1.0]
    )
    epsilon = 1e-8
    near = hydrology.seasonal_open_mean_c(np.array([-2.0 - epsilon, -2.0, -2.0 + epsilon]), np.full(3, 2.0))
    assert np.all(np.isfinite(near))
    assert near[1] == pytest.approx(0.0)
    assert near[2] == pytest.approx(0.0, abs=1e-5)


@pytest.mark.usefixtures("march_glaciers")
def test_compute_hydrology_threads_seasonal_amplitude_into_glacier_melt():
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    glacier = np.full(12, 200.0)
    precipitation = np.full(12, 100.0)
    temperature = np.full(12, -6.0)

    def run(amplitude):
        world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, np.zeros(12), glacier_depth=glacier)])
        amp = None if amplitude is None else np.full(12, amplitude)
        fields = hydrology.compute_hydrology(world, precipitation, temperature, years=1_000_000, seasonal_amplitude_at_nodes=amp)
        return fields.glacier_depth, fields.ice.ledger

    no_swing, still = run(None)
    assert np.array_equal(no_swing, run(0.0)[0])
    # No degree-days at a -6 C mean with no swing: only the depression feedback's fringe
    # melt, far below a percent of accumulation (RFC 0015 D1-D2, #789).
    assert still.surface_melt_m3 < 0.01 * still.accumulation_m3
    # A -6C mean with an 18C swing has summers well above freezing: the ice melts.
    _, melting = run(18.0)
    assert melting.surface_melt_m3 > 0.0
    assert melting.closing_m3 < still.closing_m3


@pytest.mark.parametrize("cold", [slice(28, 31), slice(33, 34)], ids=["cold-pit", "one-cold-member"])
def test_compute_hydrology_keeps_a_lake_where_its_water_is_cold(cold, monkeypatch):
    # A closed 360 m pit (nodes 28-30) on a 300 m plateau, fed by warm rain on the plateau.
    # Liquid water isn't frozen by the annual mean (RFC 0014 D2(b), #539): water reaching a pit
    # below 0 C pools as a lake there, and a lake with one cold member in its catchment keeps
    # its balance. Both used to dry: the pit's water froze into the ice at WATER_TO_ICE, and
    # `lakes._resolve` dried a whole lake when any member was frozen, so the water its inflow
    # brought was booked nowhere (#519).
    theta = 0.05 * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    temperature = np.full(40, 15.0)
    temperature[cold] = -10.0
    world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])

    # Isolate D2(b): D2(c) surface losses can exhaust an unsupplied reservoir over a million years.
    monkeypatch.setattr(hydrology, "LAKE_ICE_SUBLIMATION_M_PER_YR", 0.0)
    fields = hydrology.compute_hydrology(world, np.full(40, 800.0), temperature, years=1_000_000)
    assert not fields.is_ocean[28:31].any()
    assert fields.lake_depth[28:31].max() > 0.0
    ledger = fields.ice.ledger
    assert ledger.closing_m3 <= ledger.opening_m3 + ledger.accumulation_m3  # no water became ice

def test_compute_hydrology_sublimates_covered_lake_ice_and_closes_its_ledger(monkeypatch):
    theta = 0.05 * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    temperature = np.full(40, 15.0)
    temperature[28:31] = -10.0
    precipitation = np.full(40, 100.0)
    world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
    years = 1_000.0
    monkeypatch.setattr(hydrology, "LAKE_ICE_SUBLIMATION_M_PER_YR", 1e-4)

    with_sublimation = hydrology.compute_hydrology(
        world, precipitation, temperature, years=years, seasonal_amplitude_at_nodes=np.zeros(40)
    )
    monkeypatch.setattr(hydrology, "LAKE_ICE_SUBLIMATION_M_PER_YR", 0.0)
    without_sublimation = hydrology.compute_hydrology(
        world, precipitation, temperature, years=years, seasonal_amplitude_at_nodes=np.zeros(40)
    )

    assert with_sublimation.lake_depth[28:31].max() > 0.0
    assert np.any(with_sublimation.lake_ice_m[28:31] > 0.0)
    assert with_sublimation.lake_water_ledger.sublimation_m3 > 0.0
    assert without_sublimation.lake_water_ledger.sublimation_m3 == 0.0
    assert abs(with_sublimation.lake_water_ledger.residual_m3) <= 1e-9 * (
        with_sublimation.lake_water_ledger.external_supply_m3
    )


def test_lake_ice_thickness_follows_its_stefan_trajectory():
    cap = hydrology.LAKE_ICE_MAX_THICKNESS_M
    thickness = hydrology._lake_ice_thickness_m
    # No loss: pure Stefan growth, sqrt(h0^2 + 2 q t), held at the cap.
    assert thickness(np.array([0.5, 0.5]), np.array([0.1, 10.0]), np.zeros(2), cap, 3.0) == pytest.approx(
        [np.sqrt(0.25 + 0.6), cap]
    )
    # No growth: a linear decline that stops at zero.
    assert thickness(np.array([1.0, 1.0]), np.zeros(2), np.array([0.2, 2.0]), cap, 2.0) == pytest.approx([0.6, 0.0])
    # A long step ends at equilibrium q/loss from either side, or at the cap below it.
    q, loss = np.array([0.5, 0.5, 4.0]), np.array([0.8, 0.8, 0.8])
    assert thickness(np.array([0.0, cap, 0.0]), q, loss, cap, 1_000.0) == pytest.approx([0.625, 0.625, cap])
    # A short step lands on the integrated trajectory of dh/dt = q/h - loss.
    for h0 in (0.1, 1.9):
        h, dt = h0, 1e-4
        for _ in range(int(round(1.5 / dt))):
            h += dt * (0.5 / h - 0.8)
        assert thickness(np.array([h0]), np.array([0.5]), np.array([0.8]), cap, 1.5)[0] == pytest.approx(h, rel=1e-3)


def test_advance_lake_ice_splits_the_solved_column_without_losing_water():
    rng = np.random.default_rng(757)
    n = 200
    prior = rng.uniform(0.0, 1.834, n)
    column = np.where(rng.random(n) < 0.2, 0.0, rng.uniform(0.0, 5.0, n))
    area = rng.uniform(1e8, 1e10, n)
    mean, amplitude, cover = rng.uniform(-20.0, 10.0, n), rng.uniform(0.0, 15.0, n), rng.random(n)
    ice, liquid, budget = hydrology.advance_lake_ice(prior, column, area, mean, amplitude, cover, 100_000.0)
    assert np.allclose(ice + liquid, column, rtol=0.0, atol=1e-12)
    assert np.all(ice >= 0.0) and np.all(liquid >= 0.0)
    assert np.all(ice <= np.minimum(hydrology.LAKE_ICE_MAX_THICKNESS_M * hydrology.LAKE_ICE_DENSITY_RATIO, column) + 1e-12)
    assert np.all(ice[column == 0.0] == 0.0)  # a cell the solve left dry keeps no ice
    assert "sublimation_m3" not in budget  # the solve books it
    assert budget["freezing_m3"] - budget["melting_m3"] == pytest.approx(budget["ice_storage_change_m3"], rel=1e-12)
    assert budget["ice_storage_change_m3"] == pytest.approx(float(np.dot(ice - prior, area)), rel=1e-12)
    # A zero-length step changes no phase.
    held, kept, _ = hydrology.advance_lake_ice(prior, prior + 1.0, area, mean, amplitude, cover, 0.0)
    assert np.array_equal(held, prior) and np.allclose(kept, 1.0)


def _cold_pit(prior_ice_m=0.0):
    # A one-cell closed pit (node 29, 360 m deep) in a cold, dry band on a 300 m plateau. Warm
    # rain beyond the band runs off into it. With no seasonal swing its water is ice-covered all
    # year, so it loses water only by sublimation.
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[29] = -60.0
    temperature = np.full(40, 15.0)
    temperature[22:37] = -3.0
    precipitation = np.where(temperature < 0.0, 0.0, 800.0)
    lake_ice = np.zeros(40)
    lake_ice[29] = prior_ice_m
    plate = quad_plate(0, "continental", columns=_row(40), elevation=elevation, lake_ice_m=lake_ice)
    return World(seed=0, plates=[plate]), plate, precipitation, temperature


def _step_cold_pit(world, plate, precipitation, temperature, years, amplitude_c=0.0):
    fields = hydrology.compute_hydrology(
        world, precipitation, temperature, years=years, seasonal_amplitude_at_nodes=np.full(40, amplitude_c)
    )
    plate.set_fields_on_plate(
        lake_depth=fields.lake_depth, lake_ice_m=fields.lake_ice_m,
        glacier_depth=fields.glacier_depth, river_water_m=fields.river_water_m,
    )
    return fields


def _hold_lake_surface_at_air_temperature(monkeypatch):
    """Take lake water temperature from the air at the bed, as #754's level-lag knockout did,
    so last step's water level can't change this step's climate."""
    monkeypatch.setattr(
        hydrology, "_lake_water_surface_temperature_c",
        lambda air_c, ocean_c, is_ocean, *_: np.where(is_ocean, ocean_c, air_c),
    )


def test_lake_ice_sublimation_grows_with_step_length_at_its_closed_form_rate(monkeypatch):
    # Sublimation used to be charged after the lake solve, through an integral that saturated,
    # so a cold lake lost about 42 m whatever the step length (#754). It's now a loss in the
    # solve's implicit balance: the rate times the solved wet area times the step. From one
    # settled state, the lake loses 1,000 times as much in 1 Myr as in 1 kyr (#757).
    _hold_lake_surface_at_air_temperature(monkeypatch)
    world, plate, precipitation, temperature = _cold_pit()
    for _ in range(3):
        _step_cold_pit(world, plate, precipitation, temperature, 100_000.0)
    settled = {
        name: plate.collect(name).copy()
        for name in ("lake_depth", "lake_ice_m", "glacier_depth", "river_water_m")
    }
    sublimated_m3 = {}
    for years in (1_000.0, 100_000.0, 1_000_000.0):
        plate.set_fields_on_plate(**settled)
        fields = _step_cold_pit(world, plate, precipitation, temperature, years)
        column = fields.lake_depth + fields.lake_ice_m
        assert np.flatnonzero(column > 0.0).tolist() == [29]
        level = fields.hypsometry.level_from_depth(column[[29]], np.array([29]))[0]
        wet_area_m2 = fields.hypsometry.wet_area(level, np.array([29]))
        ledger = fields.lake_water_ledger
        assert ledger.sublimation_m3 == pytest.approx(
            hydrology.LAKE_ICE_SUBLIMATION_M_PER_YR * wet_area_m2 * years, rel=1e-9
        )
        assert ledger.evaporation_m3 == 0.0  # no open-water season
        assert abs(ledger.residual_m3) <= 1e-9 * ledger.external_supply_m3
        sublimated_m3[years] = ledger.sublimation_m3
    assert sublimated_m3[100_000.0] / sublimated_m3[1_000.0] == pytest.approx(100.0, rel=1e-6)
    assert sublimated_m3[1_000_000.0] / sublimated_m3[1_000.0] == pytest.approx(1_000.0, rel=1e-6)


def test_a_lake_frozen_to_its_bed_steps_like_the_same_water_left_liquid():
    # Last step's lake is its liquid and ice together, since floating ice displaces its own
    # mass of water. Every lagged reader (wet area for rain, water-surface temperature, spill
    # routing, river reaches) must see that column, as the solve does, or rain on a lake frozen
    # to its bed is routed as dry-land runoff (#757 review).
    def step(lake_depth, lake_ice):
        elevation = np.full(40, 300.0)
        elevation[0:15] = -100.0
        elevation[29] = 100.0
        temperature = np.full(40, 15.0)
        temperature[22:37] = 2.0
        depth, ice = np.zeros(40), np.zeros(40)
        depth[29], ice[29] = lake_depth, lake_ice
        plate = quad_plate(
            0, "continental", columns=_row(40), elevation=elevation, lake_depth=depth, lake_ice_m=ice
        )
        return hydrology.compute_hydrology(
            World(seed=0, plates=[plate]), np.full(40, 800.0), temperature, years=100_000.0,
            seasonal_amplitude_at_nodes=np.full(40, 15.0),
        )

    liquid, frozen = step(1.0, 0.0), step(0.0, 1.0)
    assert frozen.lake_external_wet_precip_m3_per_yr[29] > 0.0
    for name in ("lake_external_wet_precip_m3_per_yr", "lake_external_supply_m3_per_yr", "lake_depth", "lake_ice_m"):
        np.testing.assert_allclose(getattr(frozen, name), getattr(liquid, name), rtol=1e-12, atol=1e-12)
    assert frozen.lake_water_ledger.sublimation_m3 == pytest.approx(liquid.lake_water_ledger.sublimation_m3, rel=1e-12)


@pytest.mark.parametrize("lagged_surface", [False, True], ids=["fixed-surface", "lagged-surface"])
def test_cold_closed_lake_reaches_one_level_whatever_its_prior_ice(lagged_surface, monkeypatch):
    # With sublimation charged after the solve, this pit filled to its rim and was stripped
    # whenever its stored ice sat on equilibrium, so at fixed climate it went from full to
    # empty every step (#754). In the solve, sublimation balances supply at one level, and the
    # ice layer is a split of that column, so prior ice changes neither (#757). The lagged run
    # adds a seasonal swing, so the ice-free share, and with it evaporation against
    # sublimation, follows the water-surface temperature that last step's level sets.
    amplitude_c = 6.0 if lagged_surface else 0.0
    if not lagged_surface:
        _hold_lake_surface_at_air_temperature(monkeypatch)
    runs = []
    for prior_ice_m in (0.0, hydrology.LAKE_ICE_MAX_THICKNESS_M * hydrology.LAKE_ICE_DENSITY_RATIO):
        world, plate, precipitation, temperature = _cold_pit(prior_ice_m)
        runs.append([
            _step_cold_pit(world, plate, precipitation, temperature, 100_000.0, amplitude_c) for _ in range(6)
        ])
    for steps in runs:
        column = np.array([fields.lake_depth[29] + fields.lake_ice_m[29] for fields in steps])
        assert column.min() > 0.0  # never stripped
        # The implicit solve approaches its level from one side, so it never swings back.
        changes = np.diff(column)
        significant = np.abs(changes) > 1e-12 * column[-1]
        assert np.all(np.sign(changes[significant]) == np.sign(changes[0]))
        # Each step closes most of the remaining gap, with or without the lagged feedback.
        later = np.abs(changes[1:])[significant[1:]]
        assert np.all(later[1:] <= 0.1 * later[:-1])
        assert abs(changes[-1]) <= 1e-9 * column[-1]
        assert abs(steps[-1].lake_depth[29] - steps[-2].lake_depth[29]) <= 1e-9 * column[-1]
        final = steps[-1].lake_water_ledger
        assert final.overflow_out_m3 == 0.0
        assert final.evaporation_m3 + final.sublimation_m3 == pytest.approx(final.external_supply_m3, rel=1e-6)
        assert final.sublimation_m3 > 0.0
        assert (final.evaporation_m3 > 0.0) == lagged_surface  # only the swing opens an ice-free season
        for fields in steps:
            ledger, land = fields.lake_water_ledger, fields.land_water_ledger
            assert abs(ledger.residual_m3) <= 1e-9 * ledger.external_supply_m3
            assert abs(land.residual_m3) <= 1e-9 * (land.rain_m3 + land.meltwater_m3)
    # With the lag, the first step's surface temperature reads the opening column, prior ice
    # included, so that step's level differs by a few percent; one step later it doesn't.
    first, rel = (1, 1e-2) if lagged_surface else (0, 1e-3)
    for no_ice, full_ice in zip(runs[0][first:], runs[1][first:]):
        for name in ("lake_depth", "lake_ice_m"):
            assert getattr(full_ice, name)[29] == pytest.approx(getattr(no_ice, name)[29], rel=rel)
    for name in ("lake_depth", "lake_ice_m"):
        assert getattr(runs[1][-1], name)[29] == pytest.approx(getattr(runs[0][-1], name)[29], rel=1e-9)


@pytest.mark.parametrize(
    "cold", [slice(28, 31), slice(33, 34), slice(20, 34)], ids=["cold-pit", "one-cold-member", "cold-band"]
)
def test_compute_hydrology_delivers_every_drop_of_liquid_water_past_cold_nodes(cold, monkeypatch):
    # The closed pit of the test above, with evaporation off and every warm node at 0 C, so no
    # water is lost from lakes or in transit (river_transit_loss_fraction is 0 at or below
    # 0 C). The pit fills to its rim in the first step; from then on every m^3 of liquid
    # runoff and meltwater must reach the ocean, wherever the cold nodes are (RFC 0014 D2(b),
    # #539). On main the frozen nodes deleted 58% (cold pit), 4% (one cold member) and 53%
    # (cold band) of it.
    monkeypatch.setattr(breaching, "open_water_evaporation_mm", lambda t: np.zeros_like(np.asarray(t, dtype=float)))
    monkeypatch.setattr(hydrology, "LAKE_ICE_SUBLIMATION_M_PER_YR", 0.0)
    theta = 0.05 * np.arange(40)
    elevation = np.full(40, 300.0)
    elevation[0:15] = -100.0
    elevation[28:31] = -60.0
    temperature = np.zeros(40)
    # Below freezing on the annual mean (the D2(b) point) but mild enough with the 10 C
    # swing that its snow melts in the step: at -10 C physical ice (RFC 0015 D1, #789)
    # would sheet over in kyr and eat this bathtub sea, which is G4's problem, not this
    # test's.
    temperature[cold] = -2.0
    precipitation = np.full(40, 800.0)
    # A 10 C seasonal swing so the 0 C plateau gets rain as well as snow (RFC 0015 D1,
    # #789): at no swing a 0 C mean sits below the 1 C rain-snow threshold and nothing
    # is liquid.
    amplitude = np.full(40, 10.0)
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])
    areas = _node_areas(world)
    years = 10_000.0

    for _ in range(2):
        fields = hydrology.compute_hydrology(
            world, precipitation, temperature, years=years, seasonal_amplitude_at_nodes=amplitude
        )
        plate.set_fields_on_plate(lake_depth=fields.lake_depth, glacier_depth=fields.glacier_depth)
    assert fields.lake_depth[28:31].min() > 0.0  # the pit holds a lake
    rain_mm = (1.0 - ice_sheet.snow_fraction(temperature, amplitude)) * precipitation
    rain_mm = np.where(fields.is_ocean, 0.0, rain_mm)
    liquid_m3 = float(np.dot(rain_mm / 1000.0, areas))
    # Meltwater converts by mass (RFC 0015 D1, #789). Budyko runoff is all of the rain
    # with no evaporation at or below 0 C.
    source_m3 = (
        liquid_m3
        + fields.ice.ledger.meltwater_m3 * lithosphere.RHO_ICE / lithosphere.RHO_WATER / years
    )
    ocean_m3 = float(fields.discharge_m3_s[fields.is_ocean].sum()) * hydrology.SECONDS_PER_YEAR
    ledger = fields.lake_water_ledger
    river = fields.river_water_ledger
    # Physical ice (RFC 0015 D1, #789) stores kyr of snowfall within a step, so the glacier's
    # own storage change -- grown ice plus calved export, in water -- belongs in the balance;
    # under model units it was small enough to leave out.
    ice = fields.ice.ledger
    ice_storage_m3_per_year = (
        (ice.closing_m3 - ice.opening_m3 + ice.calved_m3)
        * lithosphere.RHO_ICE / lithosphere.RHO_WATER / years
    )
    stored_or_lost_m3_per_year = (
        ledger.storage_change_m3 + ledger.ice_storage_change_m3 + ledger.evaporation_m3
        + ledger.sublimation_m3 + river.final_storage_m3 - river.initial_storage_m3
    ) / years + ice_storage_m3_per_year
    assert ocean_m3 == pytest.approx(source_m3 - stored_or_lost_m3_per_year, rel=1e-9)
    assert abs(ledger.residual_m3) <= 1e-9 * ledger.external_supply_m3


@pytest.mark.usefixtures("march_glaciers")
def test_compute_hydrology_banks_a_cold_chains_snow_as_ice():
    # Same 12-node monotonic-descent chain as the end-to-end test, just below freezing (-3C)
    # with no seasonal swing: every drop falls as snow and none melts, so no liquid water
    # reaches the ocean and the snow stays on the ground as ice (issue #353: gathered into its
    # ice body).
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    fields = hydrology.compute_hydrology(world, np.full(12, 800.0), np.full(12, -3.0), years=1_000_000)
    ledger = fields.ice.ledger
    # Only a trace of the ice's own meltwater moves; no rain does.
    ocean_m3 = float(fields.discharge_m3_s[fields.is_ocean].sum()) * hydrology.SECONDS_PER_YEAR * 1_000_000
    assert ocean_m3 <= ledger.meltwater_m3 * (1.0 + 1e-9)
    assert ledger.accumulation_m3 > 0.0
    assert float(np.dot(fields.glacier_depth, _node_areas(world))) == pytest.approx(ledger.closing_m3)


@pytest.mark.usefixtures("march_glaciers")
def test_compute_hydrology_snows_on_an_ice_cap_above_freezing_but_lets_its_water_drain():
    # Same 12-node descending chain as the cold-chain test, at a sampled temperature of
    # +2C -- above freezing. A bare chain routes rain to the sea. An otherwise-identical chain
    # that already carries a real ice cap (200 m of glacier) still takes its precipitation as
    # snow, but liquid water on it isn't frozen in place: it drains downhill. Freezing it there let a river rebuild a warm ice plug every step (#503).
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)  # last two nodes are below sea level -> ocean
    warm = np.full(12, 2.0)
    precip = np.full(12, 800.0)

    bare = _flow_line_plate(0, theta, elevation)
    bare_fields = hydrology.compute_hydrology(World(seed=0, plates=[bare]), precip, warm, years=1_000_000)
    bare_land = ~bare_fields.is_ocean
    assert np.any(bare_fields.flow_target[bare_land] >= 0)  # control: warm -> rivers route downhill
    assert bare_fields.water_deposited[bare_fields.is_ocean].sum() > 0.0  # control: rain reaches the sea

    iced = _flow_line_plate_with_lake(0, theta, elevation, np.zeros(12), glacier_depth=np.full(12, 200.0))
    iced_fields = hydrology.compute_hydrology(World(seed=0, plates=[iced]), precip, warm, years=1_000_000)
    iced_land = ~iced_fields.is_ocean
    assert np.all(iced_fields.flow_target[iced_land] >= 0)  # water on the warm cap drains downhill
    assert iced_fields.ice.ledger.accumulation_m3 > 0.0  # the precipitation still lands as snow


def test_group_rivers_groups_confluence_and_finds_max_flow_mouth():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    assert len(rivers) == 4

    network_a = next(r for r in rivers if len(r.member_idx) == 4)
    assert sorted(network_a.member_idx.tolist()) == [0, 1, 2, 3]
    assert network_a.mouth_idx == 3  # highest fallback discharge (converted from flow_accum)
    assert network_a.flow_rate == pytest.approx(
        float(hydrology.flow_accum_mm_yr_to_discharge_m3_s(40.0))
    )
    # Two headwaters (0 and 1) merging into one mainstem -- one tributary joins it.
    assert network_a.num_tributaries == 1


def test_river_selection_stays_connected_through_positive_downstream_discharge():
    # Transit evaporation can reduce discharge below the percentile cutoff along a real river.
    # Once selected upstream, its positive-flow path remains connected to the mouth.
    selected = hydrology._extend_river_mask_downstream(
        np.array([True, False, False, False, False]),
        np.array([1, 2, 3, 4, -1]),
        np.array([10.0, 0.5, 0.01, 0.0, 0.0]),
        np.zeros(5, dtype=bool),
        np.zeros(5, dtype=bool),
    )
    assert selected.tolist() == [True, True, True, False, False]

    # Multiple selected headwaters may share a mainstem; walking it once must still retain the
    # whole shared path without revisiting it for every tributary.
    confluence = hydrology._extend_river_mask_downstream(
        np.array([True, True, False, False, False]),
        np.array([2, 2, 3, 4, -1]),
        np.array([10.0, 5.0, 0.5, 0.01, 0.0]),
        np.zeros(5, dtype=bool),
        np.zeros(5, dtype=bool),
    )
    assert confluence.tolist() == [True, True, True, True, False]


def test_river_selection_stops_downstream_extension_at_lakes_and_oceans():
    flow_target = np.array([1, 2, 3, 4, -1])
    discharge = np.ones(5)
    is_ocean = np.zeros(5, dtype=bool)
    is_lake = np.zeros(5, dtype=bool)
    is_lake[2] = True
    lake_path = hydrology._extend_river_mask_downstream(
        np.array([True, False, False, False, False]), flow_target, discharge, is_ocean, is_lake
    )
    assert lake_path.tolist() == [True, True, False, False, False]

    ocean = np.zeros(5, dtype=bool)
    ocean[2] = True
    ocean_path = hydrology._extend_river_mask_downstream(
        np.array([True, False, False, False, False]), flow_target, discharge, ocean, np.zeros(5, dtype=bool)
    )
    assert ocean_path.tolist() == [True, True, False, False, False]


def test_group_rivers_classifies_mouth_type_ocean_lake_and_other():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)

    network_a = next(r for r in rivers if len(r.member_idx) == 4)
    assert network_a.mouth_type == "ocean"  # mouth's flow_target (4) is an ocean node

    network_b = next(r for r in rivers if 5 in r.member_idx.tolist())
    assert network_b.mouth_idx == 5
    # lake_depth is read off *node 7* (what node 5 flows into), not node 5's own -- node 5
    # itself, like every real river node, is never itself is_lake.
    assert network_b.mouth_type == "lake"
    assert network_b.num_tributaries == 0  # a single unbranched node has no tributaries

    network_c = next(r for r in rivers if 8 in r.member_idx.tolist())
    assert network_c.mouth_idx == 8
    assert network_c.mouth_type == "other"  # flow_target -1: a still-dry, unresolved sink

    network_d = next(r for r in rivers if 9 in r.member_idx.tolist())
    assert network_d.mouth_idx == 9
    # flow_target[9] (node 10) is itself still dry -- only the chain's *next* hop (node 7)
    # actually carries standing water. Stopping at node 10 alone would misclassify this as
    # "other", the exact bug a real save reproduced (see group_rivers's own docstring).
    assert network_d.mouth_type == "lake"


def test_group_rivers_excludes_non_river_nodes():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    all_members = {i for r in rivers for i in r.member_idx.tolist()}
    assert 4 not in all_members  # ocean node, never is_river
    assert 6 not in all_members  # land node deliberately not classified as a river
    assert 7 not in all_members  # the lake node itself -- is_river excludes it too
    assert 10 not in all_members  # dry basin floor between node 9 and the lake, not is_river


def test_group_rivers_on_no_rivers_returns_empty_list():
    fields = _river_inspector_fields()
    fields.is_river = np.zeros(len(fields.is_river), dtype=bool)
    assert hydrology.group_rivers(fields) == []


def test_river_at_finds_the_network_nearest_a_query_point():
    fields = _river_inspector_fields()
    rivers = hydrology.group_rivers(fields)
    network_a_id = next(i for i, r in enumerate(rivers) if len(r.member_idx) == 4)
    network_b_id = next(i for i, r in enumerate(rivers) if 5 in r.member_idx.tolist())

    near_a = hydrology.river_at(fields, rivers, fields.points[0])
    assert near_a == network_a_id

    near_b = hydrology.river_at(fields, rivers, fields.points[5])
    assert near_b == network_b_id


def test_river_at_returns_none_when_there_are_no_rivers():
    fields = _river_inspector_fields()
    assert hydrology.river_at(fields, [], fields.points[0]) is None


def test_group_rivers_on_a_real_monotonic_chain_is_one_unbranched_network():
    # Reuses the same 12-node monotonic-descent fixture as the compute_hydrology end-to-end
    # test above -- a single chain has no confluence, so grouping it should yield exactly one
    # network with zero tributaries, mouthing into the ocean.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])

    fields = hydrology.compute_hydrology(world, np.full(12, 800.0), np.full(12, 15.0), years=1_000_000)
    rivers = hydrology.group_rivers(fields)

    assert len(rivers) == 1
    assert rivers[0].num_tributaries == 0
    assert rivers[0].mouth_type == "ocean"
    assert fields.is_ocean[fields.flow_target[rivers[0].mouth_idx]]


def _is_neighbour_pair(fields, a: int, b: int) -> bool:
    return b in fields.neighbor_idx[a] or a in fields.neighbor_idx[b]


def test_a_spill_hop_across_a_lake_ends_the_river_at_the_lake():
    fields = spill_hop_fields(lake_m=15.0)

    rivers = sorted(hydrology.group_rivers(fields), key=lambda r: r.mouth_idx)
    inflow, outflow = rivers
    assert sorted(inflow.member_idx.tolist()) == [0, 1, 2]
    assert (inflow.mouth_idx, inflow.mouth_type) == (2, "lake")
    assert outflow.member_idx.tolist() == [6]
    assert outflow.mouth_type == "ocean"

    _, _, _, crosses_lake = hydrology.flow_segments(fields, inflow.member_idx)
    pairs = segment_pairs(fields, inflow.member_idx)
    # Drawn into the lake's edge and no further: nothing across the water or beyond the rim.
    assert pairs == {(0, 1), (1, 2), (2, 3)}
    assert all(_is_neighbour_pair(fields, a, b) for a, b in pairs)
    assert crosses_lake.tolist() == [m == 2 for m in inflow.member_idx.tolist()]


def test_a_spill_hop_over_a_dry_basin_is_drawn_along_its_spill_path():
    # A pre-#520 cache has no should_spill; the hop is still found from spill_target.
    fields = spill_hop_fields(lake_m=0.0, store_should_spill=False)

    rivers = hydrology.group_rivers(fields)
    assert len(rivers) == 1
    assert sorted(rivers[0].member_idx.tolist()) == [0, 1, 2, 6]

    assert segment_pairs(fields, rivers[0].member_idx) == {(0, 1), (1, 2), (6, 7)} | DRY_SPILL_PATH
    step_m = geometry.angular_distance(fields.points[0], fields.points[1]) * hydrology.PLANET_RADIUS_KM * 1000.0
    assert rivers[0].length_m == pytest.approx(7 * step_m, rel=1e-3)


def test_a_united_lake_pooling_hop_is_found_from_should_spill():
    # A united lake's sinks are routed to its lowest sink, not to their spill_target
    # after same-step hierarchy resolution, so only should_spill marks the hop.
    fields = spill_hop_fields(lake_m=0.0)
    fields.spill_target[2] = 5

    assert segment_pairs(fields, [2]) == DRY_SPILL_PATH


def test_a_spill_hop_over_a_one_way_rim_edge_still_follows_its_spill_path():
    # k-NN is directed: here only node 6 lists rim 5 as a neighbour, not the other way round.
    fields = spill_hop_fields(lake_m=0.0)
    fields.neighbor_idx[5] = [4, 4]

    rivers = hydrology.group_rivers(fields)
    assert len(rivers) == 1
    assert rivers[0].mouth_type == "ocean"
    assert segment_pairs(fields, [2]) == DRY_SPILL_PATH


def test_a_spill_hop_with_no_spill_path_keeps_its_straight_segment():
    # Node 6 is reachable only through the ocean, so there's no spill path over land.
    fields = spill_hop_fields(lake_m=0.0)
    fields.neighbor_idx[5] = [4, 4]
    fields.neighbor_idx[6] = [7, 7]

    _, _, _, crosses_lake = hydrology.flow_segments(fields, np.array([2]))
    assert segment_pairs(fields, [2]) == {(2, 6)}
    assert not crosses_lake[0]
    rivers = hydrology.group_rivers(fields)
    assert len(rivers) == 1
    assert rivers[0].mouth_type == "ocean"


def test_a_spill_path_never_climbs_past_the_deepest_notch_above_the_spill_level():
    # The routed spill level is 40 m; a rim more than BREACH_MAX_NOTCH_M above that can't be
    # the pass routing used, so the search gives up and the straight segment stays.
    fields = spill_hop_fields(lake_m=0.0)
    fields.filled_elevation[2] = 40.0
    fields.elevation[5] = 40.0 + breaching.BREACH_MAX_NOTCH_M + 1.0

    assert segment_pairs(fields, [2]) == {(2, 6)}


def test_spill_paths_are_cached_per_fields_and_left_out_of_saves():
    fields = spill_hop_fields(lake_m=0.0)
    hydrology.group_rivers(fields)
    assert set(fields._spill_path_cache.paths) == {2}

    restored = pickle.loads(pickle.dumps(fields))
    assert "_spill_path_cache" not in restored.__dict__
    assert segment_pairs(restored, [2]) == DRY_SPILL_PATH


def test_a_spill_hop_from_a_sink_under_visible_water_draws_nothing():
    fields = spill_hop_fields(lake_m=15.0, sink_m=15.0)

    _, _, _, crosses_lake = hydrology.flow_segments(fields, np.array([2]))
    assert segment_pairs(fields, [2]) == set()
    assert crosses_lake[0]


def test_compute_river_speed_increases_with_slope_and_flow():
    base = hydrology.compute_river_speed(np.array([0.01]), np.array([1.0]))[0]
    steeper = hydrology.compute_river_speed(np.array([0.1]), np.array([1.0]))[0]
    more_flow = hydrology.compute_river_speed(np.array([0.01]), np.array([100.0]))[0]
    assert steeper > base
    assert more_flow > base


# -- Issue #133 phase 2: node_cloud_resample_mode in sample_is_ocean/sample_is_sea -----------


def _grid_xyz(height, width):
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_grid = np.repeat(lon_deg[None, :], height, axis=0)
    return geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))


def test_sample_is_ocean_and_is_sea_healpix_mode_agree_closely_with_kdtree():
    world = generate_world(30, num_plates=8)
    # A short step: the assertions only need a populated hydrology cache and compare the
    # two modes against the same world, so a full 1 Myr march (thousands of substeps)
    # buys nothing (#806: this test alone cost ~20 min).
    step_world(world, 10_000)  # populates world.hydrology_cache
    world_xyz = _grid_xyz(31, 61)
    fallback = np.zeros((31, 61), dtype=bool)

    world.node_cloud_resample_mode = "kdtree"
    ocean_kd = hydrology.sample_is_ocean(world, world_xyz, fallback)
    sea_kd = hydrology.sample_is_sea(world, world_xyz, fallback)

    world.node_cloud_resample_mode = "healpix"
    ocean_hp = hydrology.sample_is_ocean(world, world_xyz, fallback)
    sea_hp = hydrology.sample_is_sea(world, world_xyz, fallback)

    assert (ocean_kd == ocean_hp).mean() > 0.95
    assert (sea_kd == sea_hp).mean() > 0.95


def test_nearest_hydro_node_idx_cache_does_not_cross_contaminate_between_modes():
    """The single-entry cache (`_LAST_NEAREST_HYDRO_NODE`) must key on `resample_mode` too --
    otherwise a "kdtree" call followed immediately by a "healpix" call against the identical
    `(hydro, query_xyz)` would incorrectly serve the first mode's stale result."""
    world = generate_world(31, num_plates=8)
    # A short step: only a populated cache is needed, not a marched world (#806).
    step_world(world, 10_000)
    hydro = world.hydrology_cache
    world_xyz = _grid_xyz(21, 41)

    idx_kd = hydrology._nearest_hydro_node_idx(hydro, world_xyz, "kdtree")
    hydrology._nearest_hydro_node_idx(hydro, world_xyz, "healpix")  # overwrites the single-entry cache
    idx_kd_again = hydrology._nearest_hydro_node_idx(hydro, world_xyz, "kdtree")

    assert np.array_equal(idx_kd, idx_kd_again)  # re-fetching "kdtree" wasn't clobbered by the "healpix" call
    tree = hydrology._ocean_node_tree(hydro)
    _, expected_kd = tree.query(world_xyz.reshape(-1, 3), workers=1)
    assert np.array_equal(idx_kd, expected_kd)


def test_sample_is_ocean_fallback_before_first_hydrology_pass_is_mode_independent():
    world = generate_world(32, num_plates=8, cratons_first=False)  # no step_world -- no hydrology_cache yet
    world_xyz = _grid_xyz(11, 21)
    fallback = np.zeros((11, 21), dtype=bool)
    fallback[0, 0] = True

    for mode in ("kdtree", "healpix"):
        world.node_cloud_resample_mode = mode
        result = hydrology.sample_is_ocean(world, world_xyz, fallback)
        assert np.array_equal(result, fallback)


def test_route_downstream_conserves_across_uphill_edges_and_cycles():
    # Issue #275: spill/ice routing can point a node at a *higher* neighbour. An elevation
    # sweep used to add that flux to a node it had already visited, so it was never passed on.
    elevation = np.array([100.0, 50.0, 80.0, 10.0, 0.0, 30.0, 31.0])
    is_ocean = np.array([False, False, False, False, True, False, False])
    # 0 -> 1 -> 2 (uphill spill) -> 3 -> ocean 4; 5 <-> 6 is a cycle.
    flow_target = np.array([1, 2, 3, 4, -1, 6, 5])
    source = np.array([1.0, 2.0, 3.0, 4.0, 0.0, 5.0, 6.0])
    retain = np.array([0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0])
    routed = hydrology.route_downstream(elevation, is_ocean, flow_target, source, retain_fraction=retain)
    flux, deposited = routed
    assert np.isclose(deposited.sum(), source.sum())
    assert np.isclose(flux[2], 6.0)  # nodes 0 and 1 both reached the spill node
    assert np.isclose(deposited[3], 5.0) and np.isclose(deposited[4], 5.0)
    assert np.isclose(deposited[5] + deposited[6], 11.0)  # a cycle keeps its own load
    assert routed.cycle_settlement[5:].sum() == pytest.approx(11.0)
    assert routed.sink_settlement.sum() == 0.0
    assert routed.residual == pytest.approx(0.0)
    assert (
        routed.retained.sum() + routed.sink_settlement.sum() + routed.ocean_export.sum()
        + routed.cycle_settlement.sum() + routed.evaporated.sum()
    ) == pytest.approx(source.sum())


def test_route_downstream_names_overflow_stopping_at_a_spilling_leaf_sink():
    # A leaf sink's spill edge can send exported lake water through an uphill cycle. The route
    # result names the exact per-node remainder so the caller can transfer it through the lake
    # hierarchy instead of silently treating it as an ordinary terminal deposit (#731).
    elevation = np.array([10.0, 20.0, 0.0])
    is_ocean = np.array([False, False, True])
    target = np.array([1, 0, -1])  # node 1 is the spilling leaf sink; 1 -> 0 is its spill edge
    overflow = np.array([0.0, 7.0, 0.0])

    routed = hydrology.route_downstream(elevation, is_ocean, target, overflow)

    assert routed.cycle_settlement[[0, 1]].sum() == pytest.approx(7.0)
    assert routed.sink_settlement.sum() == 0.0
    assert routed.ocean_export.sum() == 0.0
    assert routed.residual == pytest.approx(0.0)


def test_route_downstream_names_water_reaching_a_spilling_leaf_outlet():
    # When the spill edge has an ocean receiver, the same exported volume is credited to that
    # reservoir and not left on the leaf sink (#731).
    elevation = np.array([20.0, 10.0, 0.0])
    is_ocean = np.array([False, False, True])
    target = np.array([1, 2, -1])
    overflow = np.array([0.0, 5.0, 0.0])

    routed = hydrology.route_downstream(elevation, is_ocean, target, overflow)

    assert routed.ocean_export[2] == pytest.approx(5.0)
    assert routed.cycle_settlement.sum() == 0.0
    assert routed.sink_settlement.sum() == 0.0
    assert routed.residual == pytest.approx(0.0)


def test_route_downstream_conserves_on_an_all_downhill_graph():
    rng = np.random.default_rng(0)
    n = 200
    elevation = rng.uniform(0.0, 1000.0, n)
    is_ocean = elevation < 100.0
    order = np.argsort(elevation)
    flow_target = np.full(n, -1)
    for rank, i in enumerate(order[1:], start=1):
        flow_target[i] = order[rng.integers(0, rank)]  # always a strictly lower node
    source = rng.uniform(0.0, 1.0, n)
    flux, deposited = hydrology.route_downstream(elevation, is_ocean, flow_target, source)
    assert np.isclose(deposited.sum(), source[~is_ocean].sum())


def test_route_downstream_applies_loss_and_retention_on_cycle_nodes():
    # Review on #278: a node on a flow cycle must still evaporate/retain like any other node
    # before keeping its remainder -- otherwise river water skips evaporation there.
    elevation = np.array([100.0, 30.0, 31.0, 0.0])
    is_ocean = np.array([False, False, False, True])
    flow_target = np.array([1, 2, 1, -1])  # 0 feeds the 1 <-> 2 cycle
    source = np.array([4.0, 2.0, 2.0, 0.0])
    loss = np.array([0.0, 0.5, 0.25, 0.0])
    retain = np.array([0.0, 0.5, 0.0, 0.0])
    flux, deposited = hydrology.route_downstream(
        elevation, is_ocean, flow_target, source, retain_fraction=retain, loss_fraction=loss
    )
    # Node 1 receives 2 + 4 = 6, loses half (3), retains half the rest; node 2 loses a quarter.
    assert np.isclose(deposited[1], 3.0) and np.isclose(flux[1], 1.5)
    assert np.isclose(deposited[2], 1.5) and np.isclose(flux[2], 1.5)
    assert np.isclose(deposited.sum(), 4.5)  # 8 in, 3.5 evaporated

    # Without loss, the cycle still conserves everything it received.
    _, kept = hydrology.route_downstream(elevation, is_ocean, flow_target, source, retain_fraction=retain)
    assert np.isclose(kept.sum(), source.sum())


def test_physical_water_accounts_for_transit_loss_lake_retention_and_ocean_arrival():
    # Three land nodes drain through a lake (node 1) to the ocean. The yearly source volume
    # must end in atmospheric loss, lake retention or ocean arrival exactly once (#355).
    elevation = np.array([30.0, 20.0, 10.0, 0.0])
    is_ocean = np.array([False, False, False, True])
    target = np.array([1, 2, 3, -1])
    source_m3_per_yr = np.array([10.0, 5.0, 2.0, 0.0])
    loss = np.array([0.1, 0.2, 0.0, 0.0])
    lake_retention = np.array([0.0, 0.25, 0.0, 0.0])
    through, settled = hydrology.route_downstream(
        elevation, is_ocean, target, source_m3_per_yr,
        retain_fraction=lake_retention, loss_fraction=loss,
    )
    evaporated = through[0] * loss[0] / (1.0 - loss[0]) + (
        through[1] + settled[1]
    ) * loss[1] / (1.0 - loss[1])
    assert settled[1] == pytest.approx(2.8)
    assert settled[3] == pytest.approx(10.4)
    assert source_m3_per_yr.sum() == pytest.approx(evaporated + settled.sum())


@pytest.mark.usefixtures("march_glaciers")
def test_ice_ledger_closes_across_steps_and_meltwater_never_refreezes_as_river_ice():
    # Cold, snowy uplands above a band whose summers melt snow. That meltwater used to settle on
    # the band and refreeze as river ice, added after the solve; now it flows on as water (RFC
    # 0014 D2(b), #539), and each step's closing ice is the next step's opening.
    d = 0.002
    theta = d * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    plate = _flow_line_plate(0, theta, elevation)
    world = World(seed=0, plates=[plate])
    precipitation = np.full(12, 800.0)
    temperature = np.where(np.arange(12) < 4, -20.0, -3.0)
    amplitude = np.where(np.arange(12) < 4, 0.0, 10.0)
    areas = _node_areas(world)
    ledgers = []
    for _ in range(2):
        fields = hydrology.compute_hydrology(world, precipitation, temperature, years=1_000_000, seasonal_amplitude_at_nodes=amplitude)
        ledgers.append(fields.ice.ledger)
        assert float(np.dot(fields.glacier_depth, areas)) == pytest.approx(fields.ice.ledger.closing_m3)
        plate.set_fields_on_plate(glacier_depth=fields.glacier_depth)
    first, second = ledgers
    assert first.meltwater_m3 > 0.0
    assert second.opening_m3 == pytest.approx(first.closing_m3)


def test_ice_flow_target_follows_the_ice_surface_on_the_ice_and_the_ground_off_it():
    # A dome on a flat bed: on the ice, material moves down the ice surface (away from the
    # divide), even where the bed under it doesn't slope; off the ice, the ground's own target.
    elevation = np.zeros(5)
    depth = np.array([0.0, 500.0, 1000.0, 500.0, 0.0])
    neighbor_idx = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 3]])
    ground = np.array([-1, 0, -1, 4, -1])
    target = hydrology._ice_surface_flow_target(elevation + depth, neighbor_idx, depth, ground)
    assert target[1] == 0 and target[3] == 4
    assert target[2] in (1, 3)
    assert target[0] == -1 and target[4] == -1


def test_a_lake_under_ice_on_warm_ground_stays_water():
    # Issue #353: only cold ground (annual mean below freezing) freezes a lake. A lake beside or
    # under ice on warm ground is refilled by runoff in real units every step; freezing it into
    # the ice 1:1 fed the conserved ice volume far faster than snowfall.
    d = 0.002
    theta = d * np.arange(16)
    elevation = np.array([500, 450, 400, 350, 300, 250, 200, 50, 300, 250, 200, 150, 100, 50, 10, -50], dtype=float)
    lake_depth = np.zeros(16)
    lake_depth[7] = 30.0
    world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, lake_depth, glacier_depth=np.full(16, 200.0))])
    fields = hydrology.compute_hydrology(world, np.full(16, 500.0), np.full(16, 2.0), years=1_000_000)
    ledger = fields.ice.ledger
    assert ledger.frozen_lake_m3 == 0.0
    assert ledger.closing_m3 <= ledger.opening_m3 + ledger.accumulation_m3  # the lake didn't become ice


def _river_network():
    """Twelve nodes in a line on the unit sphere's equator, draining east to an ocean at node
    11, with a lake-covered reach at node 5 and a tributary (6) joining node 7."""
    theta = np.arange(12) * 0.002
    points = np.column_stack([np.cos(theta), np.sin(theta), np.zeros(12)])
    elevation = 2000.0 - 150.0 * np.arange(12.0)
    is_ocean = np.zeros(12, dtype=bool)
    is_ocean[11] = True
    target = np.array([1, 2, 3, 4, 5, 7, 7, 8, 9, 10, 11, -1], dtype=np.int64)
    stores = ~is_ocean & (target >= 0)
    stores[5] = False
    return points, elevation, is_ocean, target, stores


def _route_rivers(runoff, melt, prior, years=100_000.0, loss=0.0, melt_loss=None, **kw):
    points, elevation, is_ocean, target, stores = _river_network()
    q = np.full(12, 50.0)
    runoff_loss = np.full(12, loss)
    meltwater_loss = runoff_loss if melt_loss is None else np.full(12, melt_loss)
    return hydrology.route_river_reservoirs(
        points, elevation, is_ocean, target, stores, runoff, melt, prior, q, years,
        runoff_loss, meltwater_loss, **kw,
    )


def test_steady_river_reservoirs_deliver_what_route_downstream_does():
    # The reservoir only delays water: in steady state each source reaches the same place, with
    # the same source-specific transit loss, as the control's sweep (RFC 0016 G0, #725). A long
    # interval makes the stored volume a negligible share of what passes.
    points, elevation, is_ocean, target, _ = _river_network()
    runoff = np.array([9.0, 1, 0, 4, 0, 2, 7, 0, 3, 0, 1, 5]) * 1e9
    melt = np.array([0.0, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 0]) * 1e8
    runoff_loss = np.linspace(0.0, 0.02, 12)
    melt_loss = np.linspace(0.03, 0.0, 12)
    through_runoff, deposited_runoff = hydrology.route_downstream(
        elevation, is_ocean, target, runoff, loss_fraction=runoff_loss
    )
    through_melt, deposited_melt = hydrology.route_downstream(
        elevation, is_ocean, target, melt, loss_fraction=melt_loss
    )
    route = _route_rivers(runoff, melt, np.zeros(12), years=1e7)
    route = hydrology.route_river_reservoirs(
        points, elevation, is_ocean, target, np.where(is_ocean, False, target >= 0), runoff, melt,
        route.storage_m3, np.full(12, 50.0), 1e7, runoff_loss, melt_loss,
    )
    assert route.routed_m3[:11] / 1e7 == pytest.approx((through_runoff + through_melt)[:11], rel=1e-6)
    assert route.settled_runoff_m3[11] / 1e7 == pytest.approx(deposited_runoff[11], rel=1e-6)
    assert route.settled_meltwater_m3[11] / 1e7 == pytest.approx(deposited_melt[11], rel=1e-6)


def test_river_reach_update_is_exact_and_each_reach_closes():
    runoff = np.full(12, 3e8)
    prior = np.linspace(1e12, 2e12, 12)
    prior[[5, 11]] = 0.0
    route = _route_rivers(runoff, np.zeros(12), prior, years=2e-5, loss=0.01)
    _, _, _, target, stores = _river_network()
    points = _river_network()[0]
    tau = hydrology.river_residence_time_s(hydrology._flow_edge_length_m(points, target), np.full(12, 50.0))
    dt = 2e-5 * hydrology.SECONDS_PER_YEAR
    decay = np.exp(-dt / tau)
    inflow_rate = route.inflow_m3 / dt
    expected = prior * decay + inflow_rate * tau * (1.0 - decay)
    assert route.storage_m3[stores] == pytest.approx(expected[stores], rel=1e-12)
    # V1 - V0 = I dt - O at every reach, and the transit loss is charged after the outflow.
    assert np.allclose(route.storage_m3 - prior, route.inflow_m3 - route.outflow_m3, atol=1e-2, rtol=1e-12)
    assert route.evaporation_m3 == pytest.approx(route.outflow_m3 * 0.01, rel=1e-12)
    ledger = route.ledger
    assert abs(ledger.residual_m3) <= 1e-10 * (ledger.initial_storage_m3 + ledger.input_m3)


def test_river_reservoir_is_independent_of_step_partition_for_steady_inflow():
    runoff = np.zeros(12)
    runoff[0] = 4e9
    prior = np.zeros(12)
    prior[3] = 7e11
    whole = _route_rivers(runoff, np.zeros(12), prior, years=4e-5)
    half = _route_rivers(runoff, np.zeros(12), prior, years=2e-5)
    second = _route_rivers(runoff, np.zeros(12), half.storage_m3, years=2e-5)
    # Upstream delay within an interval is first order in dt/tau, so compare a reach whose
    # own inflow is the steady source: the head of the river.
    assert second.storage_m3[0] == pytest.approx(whole.storage_m3[0], rel=1e-12)
    # Every drop is accounted for either way.
    assert whole.ledger.final_storage_m3 + whole.ledger.ocean_export_m3 == pytest.approx(
        second.ledger.final_storage_m3 + half.ledger.ocean_export_m3 + second.ledger.ocean_export_m3, rel=1e-12
    )


def test_river_water_never_rests_on_lake_or_ocean_nodes():
    # Inventory stranded by a rising lake or sea moves into that reservoir as an interval
    # transfer; a volume has one owner (RFC 0016 G0, #725).
    prior = np.zeros(12)
    prior[5] = 6e10
    prior[11] = 2e10
    route = _route_rivers(np.zeros(12), np.zeros(12), prior, years=1e5)
    assert route.storage_m3[5] == 0.0
    assert route.storage_m3[11] == 0.0
    # The lake node's water carries on downstream; only the reaches below it hold any now.
    assert route.storage_m3[:5].sum() == 0.0
    assert route.ledger.ocean_export_m3 + route.ledger.final_storage_m3 == pytest.approx(8e10, rel=1e-12)
    assert route.outflow_m3[11] == pytest.approx(2e10, rel=1e-12)
    assert route.ledger.residual_m3 == pytest.approx(0.0, abs=1e-3)


def test_river_water_on_a_sink_or_flow_cycle_settles_there_without_being_lost():
    points, elevation, is_ocean, target, stores = _river_network()
    target = target.copy()
    target[3], target[4] = 4, 3  # a two-node cycle
    target[9] = -1  # a land sink
    route = hydrology.route_river_reservoirs(
        points, elevation, is_ocean, target, stores, np.full(12, 1e9), np.zeros(12),
        np.full(12, 1e10), np.full(12, 50.0), 1.0, np.zeros(12), np.zeros(12),
    )
    assert route.storage_m3[[3, 4, 9]].sum() == 0.0
    assert route.settled_m3[[3, 4, 9]].min() > 0.0
    assert route.cycle_settlement_m3[[3, 4]].sum() == pytest.approx(route.settled_m3[[3, 4]].sum())
    assert route.cycle_settlement_m3[9] == 0.0
    assert route.ledger.residual_m3 == pytest.approx(0.0, abs=1e-2)
    assert (
        route.storage_m3.sum() + route.ledger.evaporation_m3 + route.ledger.ocean_export_m3
        + route.ledger.land_settlement_m3
    ) == pytest.approx(route.ledger.input_m3 + route.ledger.initial_storage_m3, abs=1e-2)


def test_exported_river_overflow_cycle_and_sink_settlement_transfer_to_ocean():
    points, elevation, is_ocean, target, stores = _river_network()
    target = target.copy()
    target[3], target[4] = 4, 3
    target[9] = -1
    route = hydrology.route_river_reservoirs(
        points, elevation, is_ocean, target, stores, np.full(12, 1e9), np.zeros(12),
        np.full(12, 1e10), np.full(12, 50.0), 1.0, np.zeros(12), np.zeros(12),
    )
    ocean_settlement_before = route.settled_m3[is_ocean].copy()
    ocean_export_before = route.ocean_export_m3.copy()
    stranded_m3, stranded_node = hydrology._deliver_unowned_river_overflow_to_ocean(route, is_ocean)

    assert stranded_m3 == 0.0 and stranded_node == -1
    assert route.ledger.unlocated_ocean_export_m3 > 0.0
    assert np.array_equal(route.settled_m3[is_ocean], ocean_settlement_before)
    assert np.array_equal(route.ocean_export_m3, ocean_export_before)
    assert route.ledger.ocean_export_m3 == pytest.approx(
        route.ocean_export_m3.sum() + route.ledger.unlocated_ocean_export_m3
    )
    assert route.ledger.land_settlement_m3 == 0.0
    assert route.ledger.ocean_export_m3 > 0.0
    assert route.cycle_settlement_m3.sum() == 0.0
    assert route.settled_m3[~is_ocean].sum() == 0.0
    assert route.ledger.residual_m3 == pytest.approx(0.0, abs=1e-2)


def test_exported_river_settlement_without_ocean_is_reported_at_largest_node():
    points, elevation, is_ocean, target, stores = _river_network()
    is_ocean = np.zeros_like(is_ocean)
    target = target.copy()
    target[3], target[4] = 4, 3
    route = hydrology.route_river_reservoirs(
        points, elevation, is_ocean, target, stores, np.full(12, 1e9), np.zeros(12),
        np.full(12, 1e10), np.full(12, 50.0), 1.0, np.zeros(12), np.zeros(12),
    )

    stranded_m3, stranded_node = hydrology._deliver_unowned_river_overflow_to_ocean(route, is_ocean)

    assert stranded_m3 == pytest.approx(route.ledger.land_settlement_m3)
    assert stranded_node in np.flatnonzero(route.settled_m3 > 0.0)


def test_stranded_overflow_warning_names_volume_and_node(caplog):
    hydrology._warn_if_stranded_overflow(12.5, 7)

    assert "12.5 m^3" in caplog.text
    assert "total; largest settlement at node 7" in caplog.text


def test_dry_reach_drains_its_stored_water_at_the_residence_cap():
    assert hydrology.river_residence_time_s(np.array([20_000.0]), np.array([0.0]))[0] == hydrology.RIVER_MAX_RESIDENCE_SECONDS
    points, elevation, is_ocean, target, stores = _river_network()
    prior = np.zeros(12)
    prior[2] = 1e9
    route = hydrology.route_river_reservoirs(
        points, elevation, is_ocean, target, stores, np.zeros(12), np.zeros(12), prior,
        np.zeros(12), 1.0, np.zeros(12), np.zeros(12),
    )
    assert route.storage_m3[2] == pytest.approx(1e9 * np.exp(-1.0), rel=1e-3)


def test_river_initial_state_is_q_tau_and_not_an_interval_input():
    points, _, _, target, stores = _river_network()
    q = np.full(12, 50.0)
    tau = hydrology.river_residence_time_s(hydrology._flow_edge_length_m(points, target), q)
    route = _route_rivers(np.zeros(12), np.zeros(12), np.zeros(12), initial_discharge_m3_s=q)
    assert route.ledger.initial_storage_m3 == pytest.approx((q * tau)[stores].sum(), rel=1e-12)
    assert route.ledger.input_m3 == 0.0
    assert route.ledger.residual_m3 == pytest.approx(0.0, abs=1e-3)
    # The per-reach readout starts from that state, so each reach's identity holds on step one.
    assert np.allclose(
        route.storage_m3 - route.initial_storage_m3, route.inflow_m3 - route.outflow_m3,
        atol=1e-2, rtol=1e-12,
    )


def _river_world():
    theta = np.arange(40) * 0.01
    elevation = np.concatenate([np.linspace(1500.0, 100.0, 30), np.full(10, -500.0)])
    plate = _flow_line_plate(1, theta, elevation)
    return World(plates=[plate], seed=0, sea_level_m=0.0)


def test_hydrology_persists_river_water_and_closes_the_rain_only_land_ledger():
    world = _river_world()
    fields = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=100_000)
    assert fields.river_water_m.sum() > 0.0
    world.plates[0].set_fields_on_plate(river_water_m=fields.river_water_m)
    again = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=100_000)
    ledger = again.land_water_ledger
    assert ledger.rain_m3 > 0.0
    assert abs(ledger.residual_m3) <= 1e-9 * ledger.rain_m3
    # Steady forcing leaves the inventory where it was.
    assert again.river_water_m.sum() == pytest.approx(fields.river_water_m.sum(), rel=1e-7)
    assert abs(again.river_water_ledger.residual_m3) <= 1e-9 * ledger.rain_m3
    assert again.river_water_ledger.max_reach_residual_m3 <= 1e-6 * again.river_water_m.max() * 1e12


def test_precipitation_partition_sums_to_diagnostic_precipitation_per_node():
    precipitation = np.array([0.0, 300.0, 700.0, 1200.0])
    is_frozen = np.array([True, True, False, False])
    area = np.array([1e9, 2e9, 3e9, 4e9])
    liquid, frozen = hydrology.partition_precipitation_m3(precipitation, is_frozen, area, 1000.0)
    assert np.array_equal(frozen > 0.0, np.array([False, True, False, False]))
    np.testing.assert_allclose(liquid + frozen, precipitation / 1000.0 * area * 1000.0, rtol=1e-15)


def test_precipitation_partition_clips_negative_precipitation_like_the_liquid_sources():
    liquid, frozen = hydrology.partition_precipitation_m3(np.array([-5.0, -5.0]), np.array([True, False]), np.ones(2), 1.0)
    assert np.all(frozen == 0.0) and np.all(liquid == 0.0)


def _glacier_ledger_run(precipitation, temperature, amplitude, steps=3, glacier_depth=None):
    theta = 0.002 * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    plate = _flow_line_plate(0, theta, elevation)
    if glacier_depth is not None:
        plate.set_fields_on_plate(glacier_depth=glacier_depth)
    world = World(seed=0, plates=[plate])
    runs = []
    for _ in range(steps):
        fields = hydrology.compute_hydrology(
            world, precipitation, temperature, years=100_000, seasonal_amplitude_at_nodes=amplitude
        )
        plate.set_fields_on_plate(glacier_depth=fields.glacier_depth, river_water_m=fields.river_water_m)
        world.glacier_water_account_m3 = fields.glacier_water_ledger.closing_m3
        runs.append(fields)
    return world, runs


@pytest.mark.usefixtures("march_glaciers")
def test_glacier_account_chains_steps_and_closes_the_full_precipitation_ledger():
    temperature = np.where(np.arange(12) < 4, -20.0, -3.0)
    amplitude = np.where(np.arange(12) < 4, 0.0, 10.0)
    world, runs = _glacier_ledger_run(np.full(12, 800.0), temperature, amplitude)
    assert runs[0].glacier_water_ledger.initialized_from_ice
    assert not runs[1].glacier_water_ledger.initialized_from_ice
    assert runs[1].glacier_water_ledger.opening_m3 == runs[0].glacier_water_ledger.closing_m3
    areas = _node_areas(world)
    opening = np.zeros(12)
    for fields in runs:
        land, glacier = fields.land_water_ledger, fields.glacier_water_ledger
        external = land.rain_m3 + glacier.external_input_m3 + land.meltwater_m3
        assert abs(glacier.residual_m3) <= 1e-12 * max(glacier.closing_m3, external)
        frozen = glacier.frozen
        assert glacier.external_input_m3 == pytest.approx(frozen.retained_m3 + frozen.bare_melted_m3, rel=1e-12)
        # The frozen terms cancel by construction: this ledger names them, it doesn't test them.
        assert hydrology.full_land_water_residual_m3(land, glacier) == pytest.approx(land.residual_m3, abs=1e-9 * external)
        # D1's partition reaches the ledgers (RFC 0015 D1, #789): this chain has no pits so
        # no lake holds water and the wet share is 0 throughout.
        assert fields.lake_depth.max() == 0.0
        cover = (np.clip(opening, 0.0, None) >= hydrology.GLACIER_VISIBLE_DEPTH_M).astype(float)
        snow_ground = ice_sheet.snow_fraction(temperature, amplitude)
        snow_ice = ice_sheet.snow_fraction(
            temperature - ice_sheet.LAPSE_RATE_C_PER_KM * np.clip(opening, 0.0, None) / 1000.0, amplitude
        )
        parts = hydrology.glacier_rain_partition_mm(
            np.full(12, 800.0), np.zeros(12), cover, snow_ground, snow_ice
        )
        _, _, ice_rain, _, bare_rain = (np.asarray(p, dtype=float) for p in parts)
        land_mask = ~fields.is_ocean
        # The land ledger's rain is exactly the rain shares over land.
        assert land.rain_m3 == pytest.approx(
            float(np.dot((ice_rain + bare_rain)[land_mask] / 1000.0, areas[land_mask])) * 100_000
        )
        # Every land snowflake is kept as ice or melts in the step (land never vanishes);
        # ocean snow can additionally feed the ice, and ocean bare snow outside any domain
        # is the pre-existing gap, so the kept-plus-melted snow sits between the land snow
        # and the all-node snow, in water.
        snow_water = (fields.ice.ledger.accumulation_m3 + fields.ice.ledger.snow_melted_m3) * (
            lithosphere.RHO_ICE / lithosphere.RHO_WATER
        )
        snow_land = float(np.dot(((800.0 - (ice_rain + bare_rain)) * land_mask) / 1000.0, areas)) * 100_000
        # The snow share is at most the whole precipitation (f_snow <= 1 everywhere), so the
        # kept-plus-melted snow is bounded above by all of P, in water. The cover only grows
        # past the opening one here, so the opening land snow bounds it below.
        snow_all = float(np.dot(np.full(12, 800.0) / 1000.0, areas)) * 100_000
        assert snow_water >= snow_land * (1.0 - 1e-9)
        assert snow_water <= snow_all * (1.0 + 1e-9)
        opening = fields.glacier_depth


@pytest.mark.usefixtures("march_glaciers")
def test_melting_ice_releases_its_account_outward_and_leaves_control_ice_to_the_solve():
    temperature = np.full(12, 8.0)
    amplitude = np.full(12, 15.0)
    _, runs = _glacier_ledger_run(np.full(12, 800.0), temperature, amplitude, steps=1, glacier_depth=np.full(12, 40.0))
    glacier, land = runs[0].glacier_water_ledger, runs[0].land_water_ledger
    ice = runs[0].ice.ledger
    assert glacier.melt_release_m3 > 0.0
    assert glacier.closing_m3 < glacier.opening_m3
    loss = (ice.surface_melt_m3 + ice.basal_melt_m3 + ice.overflow_m3 + ice.capped_m3 + ice.calved_m3) / ice.opening_m3
    assert glacier.release_m3 == pytest.approx(min(loss, 1.0) * (glacier.opening_m3 + glacier.frozen.retained_m3))
    assert abs(hydrology.full_land_water_residual_m3(land, glacier)) <= 1e-9 * (
        land.rain_m3 + land.meltwater_m3 + glacier.external_input_m3 + glacier.outward_m3
    )


@pytest.mark.usefixtures("march_glaciers")
def test_bare_ground_snow_melts_within_the_step_and_never_enters_the_account():
    # Frozen mask, summers warm enough to clear it: the snow passes straight through.
    world, runs = _glacier_ledger_run(np.full(12, 800.0), np.full(12, -1.0), np.full(12, 15.0), steps=1)
    glacier = runs[0].glacier_water_ledger
    assert glacier.frozen.land_m3 > 0.0
    assert glacier.frozen.bare_melted_m3 == pytest.approx(glacier.frozen.land_m3)
    assert glacier.closing_m3 == 0.0
    assert world.glacier_water_account_m3 == 0.0


def test_each_world_step_commits_the_glacier_account_for_the_next():
    world = generate_world(seed=42, num_plates=2, node_density=0.25, climate_density=0.25, fluid_density=0.25)
    assert world.glacier_water_account_m3 is None
    step_world(world, 100_000)
    first = world.hydrology_cache.glacier_water_ledger
    assert first.initialized_from_ice
    assert world.glacier_water_account_m3 == first.closing_m3
    assert type(world.glacier_water_account_m3) is float
    step_world(world, 100_000)
    assert not world.hydrology_cache.glacier_water_ledger.initialized_from_ice
    assert world.hydrology_cache.glacier_water_ledger.opening_m3 == first.closing_m3


def test_ice_free_world_has_an_empty_glacier_account():
    _, runs = _glacier_ledger_run(np.full(12, 800.0), np.full(12, 15.0), np.zeros(12), steps=2)
    glacier = runs[-1].glacier_water_ledger
    assert glacier.opening_m3 == glacier.closing_m3 == glacier.outward_m3 == glacier.external_input_m3 == 0.0


def test_zero_length_hydrology_leaves_the_glacier_account_unopened():
    world = _river_world()
    assert hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=0.0).glacier_water_ledger is None
    assert world.glacier_water_account_m3 is None


def test_glacier_account_stays_out_of_trapped_water():
    world = _river_world()
    before = eustasy.trapped_water_volume_m3(world)
    world.glacier_water_account_m3 = 1.0e20
    assert eustasy.trapped_water_volume_m3(world) == before


def test_save_without_river_state_starts_from_its_saved_routed_flow():
    # RFC 0016: a checkpoint's initial inventory is its saved `Q tau`, with no reach eligible
    # unless it stores water, and it is recorded as initial state rather than input.
    world = _river_world()
    fields = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=100_000)
    no_cache = fields.river_water_ledger.initial_storage_m3
    fields.discharge_m3_s = fields.discharge_m3_s * 8.0
    world.hydrology_cache = fields
    resumed = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=100_000)
    assert resumed.river_water_ledger.initial_storage_m3 > 3.0 * no_cache
    assert resumed.river_water_ledger.input_m3 == pytest.approx(fields.river_water_ledger.input_m3)
    assert abs(resumed.land_water_ledger.residual_m3) <= 1e-9 * resumed.land_water_ledger.rain_m3


def test_lake_outlet_river_holds_exported_overflow_as_inventory():
    # Reach `Q` below a lake is mostly overflow. Its water is river inventory too, and the
    # overflow pass books it, with the lake ledger's export as its input (#725).
    theta = 0.002 * np.arange(12)
    elevation = np.array([900.0, 800, 500, 450, 500, 400, 300, 200, 100, -100, -200, -300])
    depth = np.zeros(12)
    depth[3] = 60.0
    world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, depth)])
    fields = hydrology.compute_hydrology(world, np.full(12, 900.0), np.full(12, 12.0), years=100_000)
    ledger = fields.river_water_ledger
    assert ledger.lake_overflow_m3 > 0.0
    assert ledger.lake_overflow_m3 == pytest.approx(fields.lake_water_ledger.overflow_out_m3, rel=1e-9)
    assert abs(fields.land_water_ledger.residual_m3) <= 1e-9 * fields.land_water_ledger.rain_m3
    assert fields.land_water_ledger.stranded_overflow_m3 == 0.0
    assert fields.river_water_m[5:9].min() > 0.0


def test_zero_length_hydrology_reports_steady_discharge_and_leaves_river_state():
    world = _river_world()
    stepped = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=100_000)
    bootstrap = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=0.0)
    assert np.all(bootstrap.river_water_m == 0.0)
    assert bootstrap.land_water_ledger is None
    assert bootstrap.discharge_m3_s == pytest.approx(stepped.discharge_m3_s, rel=1e-6)


def test_saved_hydrology_cache_drops_per_reach_readouts_but_keeps_river_state():
    world = _river_world()
    fields = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 15.0), years=100_000)
    assert fields.river_water_outflow_m3 is not None
    state = fields.__getstate__()
    assert "river_water_outflow_m3" not in state
    assert "river_water_m" in state


def test_river_water_adds_to_trapped_surface_water():
    world = _river_world()
    before = eustasy.trapped_water_volume_m3(world)
    plate = world.plates[0]
    plate.set_fields_on_plate(river_water_m=np.full(40, 0.01))
    areas = collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))
    assert eustasy.trapped_water_volume_m3(world) - before == pytest.approx(0.01 * areas.sum(), rel=1e-12)


def test_routing_stuck_nodes_are_cycle_members_not_separate_terminal_sinks():
    target = np.array([1, 0, 0, -1, 3])
    ocean = np.zeros(5, dtype=bool)
    elevation = np.arange(5, dtype=float)
    _, stuck = hydrology._routing_order(elevation, ocean, target)
    cycle_members = np.concatenate(hydrology.flow_cycles(target))
    assert sorted(stuck) == sorted(cycle_members) == [0, 1]
    result = hydrology.route_downstream(elevation, ocean, target, np.ones(5))
    assert result.cycle_settlement.sum() == pytest.approx(3.0)
    assert result.sink_settlement[3] == pytest.approx(2.0)


def test_missing_owned_lake_supply_stays_in_land_residual(monkeypatch):
    from dataclasses import replace

    original_solve = lakes.solve_lakes_prognostic

    def drop_supply(forest, elevation, previous, years, climate, hypsometry, **kwargs):
        # Simulate a solver regression: owned supply never enters a lake reservoir.
        climate = replace(climate, inflow_m3_per_yr=np.zeros_like(elevation),
                          external_supply_m3_per_yr=np.zeros_like(elevation))
        return original_solve(forest, elevation, previous, years, climate, hypsometry, **kwargs)

    elevation = np.concatenate([np.linspace(900., 100., 15), np.linspace(100., 900., 15),
                                np.linspace(900., -300., 10)])
    world = World(seed=0, plates=[_flow_line_plate(0, np.arange(40), elevation)])
    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", False)
    rain = np.full(40, 900.)
    temperature = np.full(40, 12.)
    control = hydrology.compute_hydrology(world, rain, temperature, years=100_000)
    monkeypatch.setattr(lakes, "solve_lakes_prognostic", drop_supply)
    broken = hydrology.compute_hydrology(world, rain, temperature, years=100_000)
    assert control.lake_water_ledger.external_supply_m3 > 0
    assert broken.land_water_ledger.land_sink_settlement_m3 == control.land_water_ledger.land_sink_settlement_m3
    assert broken.land_water_ledger.residual_m3 == pytest.approx(
        control.land_water_ledger.residual_m3 + control.lake_water_ledger.external_supply_m3, rel=1e-8
    )


def test_d1_precipitation_partition_sums_to_p():
    # Each cell's precipitation splits into disjoint parts that sum to it (RFC 0015 D1,
    # #789): the lake's wet share, and on the dry share the ice-covered part's snow and
    # rain plus the bare part's snow and rain. No liquid water is deleted.
    rng = np.random.default_rng(789)
    precipitation = rng.uniform(0.0, 2000.0, size=50)
    wet = rng.uniform(0.0, 1.0, size=50)
    cover = rng.uniform(0.0, 1.0, size=50)
    snow_ground = rng.uniform(0.0, 1.0, size=50)
    snow_ice = rng.uniform(0.0, 1.0, size=50)
    parts = hydrology.glacier_rain_partition_mm(precipitation, wet, cover, snow_ground, snow_ice)
    assert len(parts) == 5
    for part in parts:
        assert np.all(np.asarray(part) >= 0.0)
        assert np.all(np.asarray(part) <= precipitation)
    assert sum(np.asarray(p) for p in parts) == pytest.approx(precipitation)
    lake, ice_snow, ice_rain, bare_snow, bare_rain = (np.asarray(p) for p in parts)
    assert lake == pytest.approx(wet * precipitation)
    dry = (1.0 - wet) * precipitation
    assert ice_snow == pytest.approx(dry * cover * snow_ice)
    assert ice_rain == pytest.approx(dry * cover * (1.0 - snow_ice))
    assert bare_snow == pytest.approx(dry * (1.0 - cover) * snow_ground)
    assert bare_rain == pytest.approx(dry * (1.0 - cover) * (1.0 - snow_ground))


@pytest.mark.usefixtures("march_glaciers")
def test_d2_march_hold_never_books_more_loss_than_exists():
    # A sheet that melts away early must not keep booking its transient melt over the
    # held remainder: held losses are clamped to the available ice, so the residual
    # stays ~zero instead of going deeply negative (#806 review).
    ice_sheet.MARCH_LOG.clear()
    theta = 0.002 * np.arange(12)
    elevation = 500.0 - 30.0 * np.arange(12)
    world = World(seed=0, plates=[_flow_line_plate_with_lake(
        0, theta, elevation, np.zeros(12), glacier_depth=np.full(12, 5.0),
    )])
    fields = hydrology.compute_hydrology(world, np.full(12, 500.0), np.full(12, 20.0), years=250_000)
    ledger = fields.ice.ledger
    available = ledger.opening_m3 + ledger.accumulation_m3
    assert available > 0.0
    assert ledger.closing_m3 >= 0.0
    assert ledger.closing_m3 < 0.01 * ledger.opening_m3
    assert abs(ledger.residual_m3) <= 0.01 * available
    entry = ice_sheet.MARCH_LOG[-1]
    assert entry["settled"] is True
    assert entry["held_years"] > 0.0


@pytest.mark.usefixtures("march_glaciers")
def test_d2_march_fast_path_on_an_ice_free_world():
    # With no ice and no snowfall the march takes the fast path -- one endpoint solve
    # over the whole step instead of settling a few substeps in and holding the rest
    # (review of #806: same answer, what main ran). The ledger still closes exactly
    # (RFC 0015 D2, #789).
    ice_sheet.MARCH_LOG.clear()
    theta = 0.05 * np.arange(40)
    elevation = np.linspace(100.0, -100.0, 40)
    world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
    fields = hydrology.compute_hydrology(world, np.full(40, 800.0), np.full(40, 25.0), years=100_000)
    entry = ice_sheet.MARCH_LOG[-1]
    assert entry["substeps"] == 0
    assert entry["settled"] is True
    assert entry["marched_years"] == pytest.approx(100_000.0)
    assert entry["held_years"] == 0.0
    assert entry["tipping_growth"] is None
    assert len(entry["volume_path_km3"]) > 0
    assert fields.ice.ledger.opening_m3 == 0.0
    assert fields.ice.ledger.closing_m3 == 0.0
    assert fields.ice.ledger.residual_m3 == 0.0
    assert np.all(fields.glacier_depth == 0.0)


@pytest.mark.usefixtures("march_glaciers")
def test_d2_march_ledger_closes_on_growing_ice():
    # The march's summed substep ledgers close against the closing ice to 0.1% of the
    # step's accumulation, and every step logs its substeps, held years, settlement and
    # volume path (RFC 0015 D2, #789).
    ice_sheet.MARCH_LOG.clear()
    theta = 0.002 * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    world = World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])
    fields = hydrology.compute_hydrology(world, np.full(12, 800.0), np.full(12, -3.0), years=100_000)
    ledger = fields.ice.ledger
    assert ledger.accumulation_m3 > 0.0
    assert abs(ledger.residual_m3) <= 1e-3 * ledger.accumulation_m3
    entry = ice_sheet.MARCH_LOG[-1]
    assert entry["substeps"] >= 1
    assert entry["marched_years"] + entry["held_years"] == pytest.approx(100_000.0)
    assert len(entry["volume_path_km3"]) > 0
    assert "tipping_growth" in entry
    assert float(np.dot(fields.glacier_depth, _node_areas(world))) == pytest.approx(ledger.closing_m3)


def _flow_line_march_world():
    theta = 0.002 * np.arange(12)
    elevation = 500.0 - 60.0 * np.arange(12)
    return World(seed=0, plates=[_flow_line_plate(0, theta, elevation)])


@pytest.mark.usefixtures("march_glaciers")
def test_march_carries_state_across_three_substeps():
    # Review of #806: the bounded march integration case -- exactly three 250 yr
    # substeps (750 yr). Each substep opens with the previous one's closing ice, and
    # the summed ledger and meltwater stay consistent. Full-length marches belong in
    # stress_tests; this is the unit gate's march coverage.
    ice_sheet.MARCH_LOG.clear()
    world = _flow_line_march_world()
    fields = hydrology.compute_hydrology(world, np.full(12, 800.0), np.full(12, -3.0), years=750)
    entry = ice_sheet.MARCH_LOG[-1]
    assert entry["substeps"] == 3
    assert entry["marched_years"] == pytest.approx(750.0)
    assert entry["held_years"] == 0.0
    assert entry["settled"] is False
    assert len(entry["volume_path_km3"]) == 4
    volumes = entry["volume_path_km3"]
    assert volumes[0] == 0.0
    assert volumes[-1] > 0.0
    ledger = fields.ice.ledger
    assert ledger.accumulation_m3 > 0.0
    assert abs(ledger.residual_m3) <= 1e-3 * ledger.accumulation_m3
    assert float(np.dot(fields.glacier_depth, _node_areas(world))) == pytest.approx(ledger.closing_m3)


@pytest.mark.usefixtures("march_glaciers")
def test_world_flag_off_takes_the_same_path_as_the_test_knob_off(monkeypatch):
    # World.simulate_glacier_march=False gates the same zero-length path as the
    # MARCH_GLACIERS test knob off (#806): one opening-state solve, ice frozen.
    precip, temp, years = np.full(12, 800.0), np.full(12, -3.0), 10_000
    ice_sheet.MARCH_LOG.clear()
    hydrology.compute_hydrology(_flow_line_march_world(), precip, temp, years=years)
    # Positive control: the default world really marches, so the equality below
    # proves the flag gates the path rather than both arms being dead code.
    assert ice_sheet.MARCH_LOG[-1]["substeps"] >= 1
    flagged = _flow_line_march_world()
    flagged.simulate_glacier_march = False
    off_by_flag = hydrology.compute_hydrology(flagged, precip, temp, years=years)
    assert ice_sheet.MARCH_LOG[-1]["substeps"] == 0
    knobbed = _flow_line_march_world()
    monkeypatch.setattr(hydrology, "MARCH_GLACIERS", False)
    off_by_knob = hydrology.compute_hydrology(knobbed, precip, temp, years=years)
    assert ice_sheet.MARCH_LOG[-1]["substeps"] == 0
    np.testing.assert_array_equal(off_by_flag.glacier_depth, off_by_knob.glacier_depth)


@pytest.mark.usefixtures("march_glaciers")
def test_march_takes_the_fast_path_with_no_ice_and_no_snow():
    # Review of #806: with no ice and nothing that could start any, the march would run
    # its minimum 8 substeps to settle on an unchanged volume; one endpoint solve is
    # the same answer. Both a dry world and warm rain (far above freezing, so no
    # snowfall whatever the seasonal amplitude) must take it.
    for precip, temp in ((np.zeros(12), 15.0), (np.full(12, 800.0), 50.0)):
        ice_sheet.MARCH_LOG.clear()
        fields = hydrology.compute_hydrology(
            _flow_line_march_world(), precip, np.full(12, temp), years=10_000
        )
        entry = ice_sheet.MARCH_LOG[-1]
        assert entry["substeps"] == 0
        assert entry["marched_years"] == pytest.approx(10_000.0)
        assert entry["settled"] is True
        assert np.all(fields.glacier_depth == 0.0)
        assert fields.ice.ledger.accumulation_m3 == 0.0


def test_march_sea_level_resolve_uses_the_seeded_mask(monkeypatch):
    # Review of #806 item 6a: the march's per-substep re-solve must take the seeded
    # path (march's own current-step mask, #806 review), never silently fall through
    # to the cache path on a misaligned node cloud. Needs a world with a real ocean
    # budget: on the 12-node flow-line fixture the ice debit exceeds its puddle
    # ocean, the budget clamps to 0, and the bisection early-returns before any mask
    # call at all.
    calls = []
    real = eustasy._seeded_connected_mask

    def spy(elevations, sea_level_m, neighbor_idx, seed):
        calls.append((len(elevations), len(seed)))
        return real(elevations, sea_level_m, neighbor_idx, seed)

    monkeypatch.setattr(eustasy, "_seeded_connected_mask", spy)
    world = generate_world(seed=42, num_plates=2, node_density=0.25, climate_density=0.25, fluid_density=0.25)
    step_world(world, 5_000)
    assert calls, "no re-solve reached the seeded mask"
    assert all(n == s for n, s in calls)


@pytest.mark.usefixtures("march_glaciers")
def test_march_skips_the_final_sea_level_resolve(monkeypatch):
    # Review of #806: the last substep's re-solve feeds no further solve (each solve
    # reads the previous level) and nothing downstream reads the march's own level,
    # so a two-substep march re-solves exactly once.
    calls = []
    real = eustasy._solve_sea_level_connected

    def spy(world, elevations, areas, budget, **kwargs):
        if kwargs.get("seed_is_ocean") is not None:
            calls.append(1)
        return real(world, elevations, areas, budget, **kwargs)

    monkeypatch.setattr(eustasy, "_solve_sea_level_connected", spy)
    hydrology.compute_hydrology(_flow_line_march_world(), np.full(12, 800.0), np.full(12, -3.0), years=500)
    assert len(calls) == 1


@pytest.mark.usefixtures("march_glaciers")
def test_step_world_honors_the_march_flag():
    # The flag reaches the march through a real step, not just compute_hydrology.
    kw = dict(seed=42, num_plates=2, node_density=0.25, climate_density=0.25, fluid_density=0.25)
    ice_sheet.MARCH_LOG.clear()
    step_world(generate_world(**kw), 10_000)
    assert ice_sheet.MARCH_LOG[-1]["substeps"] >= 1
    ice_sheet.MARCH_LOG.clear()
    world = generate_world(**kw)
    world.simulate_glacier_march = False
    step_world(world, 10_000)
    assert ice_sheet.MARCH_LOG[-1]["substeps"] == 0


@pytest.mark.usefixtures("march_glaciers")
def test_march_profiler_spans_fire_only_above_the_cost_gate():
    # The march's discrete phases (graph build, ice solves, sea-level re-solves) earn
    # profiler spans only when the step may cost 500ms+ (MARCH_PROFILE_MIN_YEARS, #801);
    # short marches stay out of recorded trees.
    with trace.recording() as rec:
        hydrology.compute_hydrology(_flow_line_march_world(), np.full(12, 800.0), np.full(12, -3.0), years=5_000)
    labels = [" > ".join(path) for path in rec.snapshot()]
    assert any("Building the glacier-march graph" in p for p in labels)
    assert any("Solving march ice" in p for p in labels)
    assert any("Re-solving march sea level" in p for p in labels)

    # Just under the gate: a full march runs but records no phase spans.
    with trace.recording() as rec:
        hydrology.compute_hydrology(_flow_line_march_world(), np.full(12, 800.0), np.full(12, -3.0), years=4_999)
    labels = [" > ".join(path) for path in rec.snapshot()]
    assert not any("glacier-march" in p or "march ice" in p or "march sea level" in p for p in labels)

    with trace.recording() as rec:
        hydrology.compute_hydrology(_flow_line_march_world(), np.full(12, 800.0), np.full(12, -3.0), years=100)
    labels = [" > ".join(path) for path in rec.snapshot()]
    assert not any("glacier-march" in p or "march ice" in p or "march sea level" in p for p in labels)


def _march_with(monkeypatch, cost, years, world=None):
    """One march of the flow-line fixture under a #822 cost option; returns
    (fields, march log entry)."""
    monkeypatch.setattr(ice_sheet, "MARCH_COST", cost)
    ice_sheet.MARCH_LOG.clear()
    fields = hydrology.compute_hydrology(
        world or _flow_line_march_world(), np.full(12, 800.0), np.full(12, -3.0), years=years
    )
    return fields, ice_sheet.MARCH_LOG[-1]


@pytest.mark.usefixtures("march_glaciers")
def test_march_substep_cap_holds_the_remainder(monkeypatch):
    # #822 option 1: an unsettled march stops at the cap and holds the rest of the step
    # on its window's mean flows, like a settled hold.
    fields, entry = _march_with(monkeypatch, ice_sheet.MarchCost(max_substeps=6), 10_000)
    assert entry["capped"] is True
    assert entry["settled"] is False
    assert entry["substeps"] == 6
    assert entry["marched_years"] == pytest.approx(1_500.0)
    assert entry["held_years"] == pytest.approx(8_500.0)
    # Held ice doesn't move: the closing ice is the capped substep's.
    assert float(np.dot(fields.glacier_depth, _node_areas(_flow_line_march_world()))) == pytest.approx(
        fields.ice.ledger.closing_m3
    )
    assert entry["volume_path_km3"][-1] == pytest.approx(fields.ice.ledger.closing_m3 / 1e9, rel=1e-5)


@pytest.mark.usefixtures("march_glaciers")
def test_march_substep_cap_endpoint_closes_its_ledger(monkeypatch):
    # #822 option 1, endpoint fill: one backward-Euler substep over the remainder after
    # the cap, so the step is fully marched and the ice ledger closes like D2's.
    fields, entry = _march_with(monkeypatch, ice_sheet.MarchCost(max_substeps=6, cap_fill="endpoint"), 10_000)
    assert entry["capped"] is True
    assert entry["substeps"] == 7
    assert entry["marched_years"] == pytest.approx(10_000.0)
    assert entry["held_years"] == 0.0
    assert entry["max_substep_years"] == pytest.approx(8_500.0)
    ledger = fields.ice.ledger
    assert ledger.accumulation_m3 > 0.0
    assert abs(ledger.residual_m3) <= 1e-3 * ledger.accumulation_m3


@pytest.mark.usefixtures("march_glaciers")
def test_adaptive_march_grows_substeps_in_pairs(monkeypatch):
    # #822 option 2: with growth forced on every pair, substeps double per pair up to
    # the cap, the step still sums to its length, and the ledger closes.
    cost = ice_sheet.MarchCost(adaptive=True, adaptive_grow_below=10.0, adaptive_shrink_above=20.0,
                               adaptive_max_substep_years=1_000.0)
    fields, entry = _march_with(monkeypatch, cost, 10_000)
    # The first settle check needs five volumes: four 250 yr substeps, a 500 yr pair,
    # then 1000 yr substeps for the last 8000 yr.
    assert entry["substeps"] == 14
    assert entry["max_substep_years"] == pytest.approx(1_000.0)
    assert entry["marched_years"] + entry["held_years"] == pytest.approx(10_000.0)
    ledger = fields.ice.ledger
    assert abs(ledger.residual_m3) <= 1e-3 * ledger.accumulation_m3


def test_adaptive_substep_shrinks_back_on_a_jump():
    cost = ice_sheet.MarchCost(adaptive=True)
    assert hydrology._adapt_substep_years(cost, 250.0, 250.0, None) == 250.0
    assert hydrology._adapt_substep_years(cost, 250.0, 250.0, 1e-3) == 500.0
    assert hydrology._adapt_substep_years(cost, 4_000.0, 250.0, 1e-3) == 4_000.0
    assert hydrology._adapt_substep_years(cost, 2_000.0, 250.0, 2e-2) == 2_000.0
    assert hydrology._adapt_substep_years(cost, 2_000.0, 250.0, 0.1) == 250.0


@pytest.mark.usefixtures("march_glaciers")
def test_march_sea_level_skip_counts_resolves(monkeypatch):
    # #822 option 3: the D2 march re-solves after every substep but the last; a skip
    # threshold no debit can reach re-solves none, and the substeps are unchanged.
    _, d2 = _march_with(monkeypatch, ice_sheet.MarchCost(), 1_000)
    assert d2["substeps"] == 4
    assert (d2["sea_level_resolves"], d2["sea_level_skipped"]) == (3, 0)
    _, skip = _march_with(monkeypatch, ice_sheet.MarchCost(sea_level_skip_m=1e12), 1_000)
    assert (skip["sea_level_resolves"], skip["sea_level_skipped"]) == (0, 3)


@pytest.mark.usefixtures("march_glaciers")
def test_local_dijkstra_march_is_bit_identical_to_d2(monkeypatch):
    # #822 option 3: the local inside-distance search is exact, so a march with it is
    # the D2 march bit for bit.
    d2, _ = _march_with(monkeypatch, ice_sheet.MarchCost(), 2_000)
    local, _ = _march_with(monkeypatch, ice_sheet.MarchCost(local_dijkstra=True), 2_000)
    np.testing.assert_array_equal(local.glacier_depth, d2.glacier_depth)
    assert local.ice.ledger.closing_m3 == d2.ice.ledger.closing_m3


def _watch(cost):
    return hydrology._SettleWatch(cost, np.ones(4), 250.0, ice_sheet.MARCH_HOLD_WINDOW, ice_sheet.MARCH_SETTLE_PATIENCE)


def _feed(watch, states, sub=250.0, remaining=1e5):
    volumes = [float(np.sum(states[0]))]
    for i, depth in enumerate(states[1:], start=1):
        volumes.append(float(np.sum(depth)))
        settled, next_sub = watch.observe(np.asarray(depth, dtype=float), volumes, sub, remaining)
        if settled:
            return i, next_sub
        sub = next_sub or sub
    return None, sub


def test_settle_watch_is_a_no_op_at_the_d2_default():
    watch = _watch(ice_sheet.MarchCost())
    assert watch.window == ice_sheet.MARCH_HOLD_WINDOW
    cycle = [np.array([1.0, 2.0, 3.0, 4.0]) * (1 + (i % 3)) for i in range(40)]
    assert _feed(watch, cycle) == (None, 250.0)
    assert watch.settled_by is None and watch.hold_window == ice_sheet.MARCH_HOLD_WINDOW


def test_recurrence_settles_an_exact_cycle_on_whole_periods():
    # #822: a period-4 cycle (pair means never coincide on it either) settles once it has
    # repeated MARCH_SETTLE_PATIENCE substeps running; the hold covers whole periods.
    base = np.array([10.0, 20.0, 30.0, 40.0])
    cycle = [base * (1.0, 1.01, 0.98, 1.03)[i % 4] for i in range(40)]
    watch = _watch(ice_sheet.MarchCost(recurrence=True))
    settled_at, _ = _feed(watch, cycle)
    # Substep 5 is the first to see its state 4 substeps back (the opening isn't kept).
    assert settled_at == 5 + ice_sheet.MARCH_SETTLE_PATIENCE - 1
    assert watch.settled_by == "recurrence" and watch.period == 4
    assert watch.hold_window == 4  # only one whole period has run when it settles


def test_recurrence_never_settles_a_drift():
    states = [np.array([10.0, 20.0, 30.0, 40.0]) * (1.0 + 1e-4 * i) for i in range(60)]
    watch = _watch(ice_sheet.MarchCost(recurrence=True))
    assert _feed(watch, states)[0] is None


def test_trend_holds_a_decaying_change_and_grows_while_it_shrinks():
    # Window means approach a fixed state geometrically (ratio 0.5): the first ratio
    # predicts ~D of change still to come, which is too much at first (grow), then holds.
    target = np.array([10.0, 20.0, 30.0, 40.0])
    states = [target * (1.0 + 0.2 * 0.5 ** (i // 2)) for i in range(40)]
    watch = _watch(ice_sheet.MarchCost(trend=True, trend_window=2))
    settled_at, _ = _feed(watch, states)
    assert watch.settled_by == "trend" and watch.hold_window == 2
    assert settled_at is not None and settled_at % 2 == 0
    grown = _watch(ice_sheet.MarchCost(trend=True, trend_window=2, trend_tol=1e-12))
    assert _feed(grown, states[:12])[1] > 250.0


def test_trend_drops_back_on_a_growing_change():
    states = [np.array([10.0, 20.0, 30.0, 40.0]) * (1.0 + 1e-3 * 3 ** (i // 2)) for i in range(12)]
    watch = _watch(ice_sheet.MarchCost(trend=True, trend_window=2))
    assert _feed(watch, states, sub=1_000.0) == (None, 250.0)


def test_march_cost_rejects_trend_with_adaptive():
    with pytest.raises(ValueError):
        ice_sheet.MarchCost(trend=True, adaptive=True)
    with pytest.raises(ValueError):
        ice_sheet.MarchCost(trend=True, trend_window=3)


@pytest.mark.usefixtures("march_glaciers")
@pytest.mark.parametrize("preset", ["cap-hold-recur", "trend-adapt", "trend-hold"])
def test_new_settle_presets_cover_the_step(monkeypatch, preset):
    fields, entry = _march_with(monkeypatch, ice_sheet.march_cost(preset), 10_000)
    assert entry["marched_years"] + entry["held_years"] == pytest.approx(10_000.0)
    assert entry["settle_test"] in (None, "pair", "recurrence", "trend")


def test_breach_surge_fires_when_a_dry_basin_opens_within_step_and_stops_when_reheld():
    # Issue #506: a dry basin (prev_lake_depth == 0) that fills and spills within the step
    # gets its breach surge this step, while a previously full basin whose level falls below
    # its rim under evaporation gets none.
    elevation = np.array(
        [500.0, 480.0, 460.0, 440.0, 420.0, 400.0, 380.0, 360.0, 340.0, 320.0]
        + [50.0]
        + [250.0, 300.0, 260.0, 200.0, 150.0, 100.0, 50.0, 20.0]
        + [-10.0, -50.0, -100.0]
    )
    n = len(elevation)
    theta = 0.002 * np.arange(n)

    dry_world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, np.zeros(n))])
    opened = hydrology.compute_hydrology(dry_world, np.full(n, 2000.0), np.full(n, 5.0), years=10_000)
    assert opened.should_spill[10]
    assert opened.breach_surge[10] == pytest.approx(hydrology.LAKE_BREACH_EROSION_COEFFICIENT)
    assert opened.breach_surge_accum[14] > 0.0

    wet_depth = np.zeros(n)
    wet_depth[10] = 160.0
    wet_world = World(seed=0, plates=[_flow_line_plate_with_lake(0, theta, elevation, wet_depth)])
    reheld = hydrology.compute_hydrology(wet_world, np.zeros(n), np.full(n, 30.0), years=1_000)
    assert not reheld.should_spill[10]
    assert np.all(reheld.breach_surge == 0.0)
    assert np.all(reheld.breach_surge_accum == 0.0)


def test_breach_surge_accum_leaves_no_flux_stranded_on_flow_cycles():
    # Issue #506: if breach_surge enters a flow cycle, route_downstream clears erosion_only_flux
    # on the cycle nodes alongside settling through_flux into cycle_settlement.
    elevation = np.array([40.0, 30.0, 20.0, 10.0, -10.0])
    is_ocean = np.array([False, False, False, False, True])
    # Node 0 drains into the 1 <-> 2 cycle; node 3 drains to ocean 4.
    flow_target = np.array([1, 2, 1, 4, -1])
    source = np.array([100.0, 10.0, 10.0, 50.0, 0.0])
    surge = np.array([100.0, 0.0, 10.0, 50.0, 0.0])

    routed = hydrology.route_downstream(
        elevation, is_ocean, flow_target, source, erosion_only_source_amount=surge
    )
    _, _, breach_surge_accum = routed
    assert routed.cycle_settlement[[1, 2]].sum() == pytest.approx(120.0)
    assert breach_surge_accum[0] == pytest.approx(100.0)
    assert np.all(breach_surge_accum[[1, 2]] == 0.0)
    assert breach_surge_accum[3] == pytest.approx(50.0)


@pytest.mark.usefixtures("march_glaciers")
def test_real_march_settles_by_recurrence_on_the_flow_line_cycle(monkeypatch):
    # With the pair test off (tolerance 0), the fixture falls into an exact period-3 cycle
    # that only the recurrence test sees; the rest of the step is held (#839 review).
    monkeypatch.setattr(ice_sheet, "MARCH_SETTLE_TOLERANCE", 0.0)
    _, entry = _march_with(monkeypatch, ice_sheet.MarchCost(recurrence=True), 200_000)
    assert entry["settled"] is True and entry["settle_test"] == "recurrence"
    assert entry["period"] == 3
    assert entry["held_years"] > 0.0
    assert entry["marched_years"] + entry["held_years"] == pytest.approx(200_000.0)


@pytest.mark.usefixtures("march_glaciers")
def test_real_march_holds_on_the_trend_test(monkeypatch):
    monkeypatch.setattr(ice_sheet, "MARCH_SETTLE_TOLERANCE", 0.0)
    fields, entry = _march_with(monkeypatch, ice_sheet.MarchCost(trend=True, trend_grow=False), 200_000)
    assert entry["settled"] is True and entry["settle_test"] == "trend"
    assert entry["period"] is None and entry["trend_changes"]
    assert entry["marched_years"] + entry["held_years"] == pytest.approx(200_000.0)
    # The held ice is the last substep's, so the ledger's closing ice is the march's final ice.
    assert fields.ice.ledger.closing_m3 == pytest.approx(entry["volume_path_km3"][-1] * 1e9, rel=1e-2)


@pytest.mark.usefixtures("march_glaciers")
def test_default_march_logs_pair_settles_only(monkeypatch):
    _, entry = _march_with(monkeypatch, ice_sheet.MarchCost(), 10_000)
    assert entry["settle_test"] in ("pair", None)
    assert entry["period"] is None and entry["trend_changes"] == []


def test_recurrence_keeps_enough_substeps_for_its_longest_period():
    assert _watch(ice_sheet.MarchCost(recurrence=True)).window == 2 * 64
    assert _watch(ice_sheet.MarchCost(trend=True)).window == 30


def test_trend_settles_a_still_growing_sheet_near_the_end_of_the_step():
    # A steadily growing sheet has no shrinking ratio, so the projection is the change
    # times the windows left: with a hundredth of a window left, it is below tolerance.
    states = [np.array([10.0, 20.0, 30.0, 40.0]) * (1.0 + 0.01 * i) for i in range(12)]
    cost = ice_sheet.MarchCost(trend=True, trend_window=2)
    assert _feed(_watch(cost), states, remaining=1e5)[0] is None
    near_end = _watch(cost)
    assert _feed(near_end, states, remaining=5.0)[0] is not None and near_end.settled_by == "trend"


def test_march_cost_rejects_a_recurrence_period_below_one():
    with pytest.raises(ValueError):
        ice_sheet.MarchCost(recurrence=True, recurrence_max_period=0)

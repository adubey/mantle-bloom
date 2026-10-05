import numpy as np

from app import erosion, faults, geometry, plates
from app.elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M
from app.world import World, generate_world as _generate_world


def generate_world(*args, **kwargs):
    kwargs.setdefault("surface", "lines")
    return _generate_world(*args, **kwargs)


def _unloaded_elevation(world) -> np.ndarray:
    """Elevation with the ice-load depression (issue #275 phase 3) taken back out, so a test
    about erosion isn't fooled by a glaciated node sinking under its ice."""
    _, elevation, _, _, _, plates_in_order = erosion._gather_nodes(world)
    return elevation - plates.collect_all_ice_load_deflection(plates_in_order)


def test_climate_grid_indices_matches_build_grid_convention():
    # Row 0 = north pole, row increases southward; column increases eastward from lon=-180.
    world_xyz = np.array(
        [
            geometry.latlon_to_xyz(np.radians(89.0), np.radians(-179.0)),
            geometry.latlon_to_xyz(np.radians(-89.0), np.radians(179.0)),
            geometry.latlon_to_xyz(0.0, 0.0),
        ]
    )
    row, col = erosion.climate_grid_indices(world_xyz, height=90, width=180)
    assert row[0] == 0
    assert row[1] == 89
    assert col[0] == 0
    assert col[2] == 90  # lon=0 falls in the middle column


def test_compute_slope_zero_for_flat_cluster():
    # A tight cluster of same-elevation points -- nobody is lower than anybody else, so
    # both the slope and its cap (drop to lowest neighbor) must be exactly 0.
    theta = np.linspace(0.0, 0.01, 8)
    points = geometry.local_xyz(np.zeros_like(theta), theta)
    elevation = np.full(8, 500.0)
    slope, drop_m = erosion.compute_slope(points, elevation)
    assert np.allclose(slope, 0.0)
    assert np.allclose(drop_m, 0.0)


def test_compute_slope_matches_known_gradient():
    # 6 points along a line with elevation decreasing monotonically by 100m per step --
    # point 0's SLOPE_NEIGHBOR_COUNT=4 nearest neighbors are points 1-4, and since
    # elevation decreases monotonically, the *lowest* of those 4 is point 4 (farthest in
    # this small set, but lowest).
    d = 0.01
    theta = d * np.arange(6)
    points = geometry.local_xyz(np.zeros_like(theta), theta)
    elevation = 500.0 - 100.0 * np.arange(6)

    slope, drop_m = erosion.compute_slope(points, elevation)

    expected_drop = elevation[0] - elevation[4]
    expected_run_m = geometry.angular_distance(points[0], points[4]) * plates.PLANET_RADIUS_KM * 1000.0
    assert np.isclose(drop_m[0], expected_drop)
    assert np.isclose(slope[0], expected_drop / expected_run_m)


def test_weathering_relief_factor_suppresses_flat_terrain_and_saturates_on_steep_terrain():
    # Regression check for the "flat coastal plain weathers exactly as fast as a mountain"
    # bug: weathering must scale down toward 0 as slope -> 0 (same clip-and-normalize idiom
    # humidity_norm already uses), and saturate to full strength on genuinely steep terrain,
    # rather than being a flat function of wind/humidity alone.
    d = 0.01
    theta = d * np.arange(6)
    points = geometry.local_xyz(np.zeros_like(theta), theta)

    flat_elevation = np.full(6, 30.0)  # uniform -- see test_compute_slope_zero_for_flat_cluster
    flat_slope, _ = erosion.compute_slope(points, flat_elevation)
    flat_relief_factor = np.clip(flat_slope / erosion.WEATHERING_RELIEF_REFERENCE_SLOPE, 0.0, 1.0)
    assert np.allclose(flat_relief_factor, 0.0)

    steep_elevation = 30.0 + 2000.0 * np.arange(6)  # steep monotonic gradient
    steep_slope, _ = erosion.compute_slope(points, steep_elevation)
    steep_relief_factor = np.clip(steep_slope / erosion.WEATHERING_RELIEF_REFERENCE_SLOPE, 0.0, 1.0)
    # Every point except the global minimum (point 0, which has no lower neighbor at all)
    # should be steep enough to saturate weathering to full strength.
    assert np.allclose(steep_relief_factor[1:], 1.0)


def test_submarine_erosion_scales_with_slope_and_depth_only_over_ocean():
    # Flat sea floor (slope 0) erodes nothing; a steeper, deeper submerged scarp erodes more
    # than a steeper, shallower one; subaerial nodes are untouched by this term.
    is_ocean = np.array([True, True, True, False])
    elevation = np.array([-100.0, -100.0, -5000.0, 500.0])
    slope = np.array([0.0, 0.03, 0.03, 0.03])
    amt = erosion.submarine_erosion_amount(elevation, slope, is_ocean, dt_myr=5.0)
    assert amt[0] == 0.0  # flat
    assert amt[2] > amt[1] > 0.0  # deeper (more pressure) erodes faster at the same slope
    assert amt[3] == 0.0  # subaerial


def test_coastal_erosion_confined_to_band_and_peaks_near_freezing():
    elevation = np.array([0.0, 0.0, 150.0, 5000.0, -5000.0])
    freezing = np.full(5, erosion.COASTAL_FROST_PEAK_C)
    warm = np.full(5, 30.0)
    at_freezing = erosion.coastal_erosion_amount(elevation, freezing, dt_myr=1.0)
    when_warm = erosion.coastal_erosion_amount(elevation, warm, dt_myr=1.0)
    assert at_freezing[0] > when_warm[0] > 0.0  # frost adds on top of wave attack at the shore
    assert at_freezing[2] > 0.0 and at_freezing[2] < at_freezing[0]  # in-band but tapering
    assert at_freezing[3] == 0.0 and at_freezing[4] == 0.0  # far above / far below sea level


def test_apply_erosion_can_lower_a_submerged_range():
    # End to end: over a generated world, some ocean node loses net elevation once submarine +
    # coastal erosion are in play (they were entirely erosion-exempt before).
    world = generate_world(seed=20, num_plates=8)
    _, elevation_before, _, _, _, _ = erosion._gather_nodes(world)
    is_ocean = elevation_before <= 0.0
    erosion.apply_erosion(world, years=5_000_000)
    _, elevation_after, _, _, _, _ = erosion._gather_nodes(world)
    assert len(elevation_after) == len(elevation_before)
    assert np.any(elevation_after[is_ocean] < elevation_before[is_ocean] - 1e-6)


def test_apply_erosion_thins_crust_where_it_erodes_and_isostasy_compensates():
    # Erosion now books its rock removal against Hc and only lets the isostatically
    # compensated fraction reach the surface -- so a node that erodes drops its Hc by more
    # than its surface, and `elevation` stays a faithful readout of isostatic_elevation(Hc).
    from app import lithosphere

    world = generate_world(seed=24, num_plates=8)
    _, elev_before, _, _, _, plates_in_order = erosion._gather_nodes(world)
    hc_before = plates.collect_all_crustal_thickness(plates_in_order)
    assert np.all(hc_before > 0.0)  # v2 world

    erosion.apply_erosion(world, years=5_000_000)

    _, _, _, _, _, plates_after = erosion._gather_nodes(world)
    elev_after = _unloaded_elevation(world)
    hc_after = plates.collect_all_crustal_thickness(plates_after)

    eroded = elev_after < elev_before - 5.0  # nodes that lost real height
    assert np.any(eroded)
    # The crustal column fell by more than the surface did (the rest is isostatic rebound).
    assert np.all((hc_before - hc_after)[eroded] > (elev_before - elev_after)[eroded])

    # elevation still tracks the Airy readout, per-plate (rho_c is per plate).
    resid = []
    for plate in plates_after:
        z = lithosphere.isostatic_elevation(
            plate.collect("crustal_thickness_m"),
            plate.collect("mantle_lithosphere_thickness_m"),
            lithosphere.crust_density(plate.crust_type),
        )
        resid.append(np.abs(z - plate.collect("elevation")))
    resid = np.concatenate(resid)
    assert np.median(resid) < 1.0
    assert np.percentile(resid, 99) < 25.0


def test_apply_erosion_keeps_elevation_finite_and_changing():
    world = generate_world(seed=21, num_plates=8)
    _, elevation_before, _, _, _, _ = erosion._gather_nodes(world)
    erosion.apply_erosion(world, years=5_000_000)
    _, elevation_after, _, _, _, _ = erosion._gather_nodes(world)
    # Erosion/deposition/weathering can now raise *or* lower any given node (deposition can
    # outweigh local erosion at a floodplain or delta) -- just confirm the pass actually did
    # something and didn't produce garbage.
    assert np.all(np.isfinite(elevation_after))
    assert not np.allclose(elevation_after, elevation_before)


def test_apply_erosion_stamps_geomorphic_provenance_but_leaves_a_sticky_structural_code():
    from app.elevation_lines import ELEV_CHANGE_COLLISION, ELEV_CHANGE_NONE

    world = generate_world(seed=21, num_plates=8)
    # Pre-stamp every land node's provenance as a structural (collision) code -- ordinary
    # background erosion, one step, must leave the overwhelming majority of them alone (only a
    # node with a genuinely large net geomorphic step overrides a structural code).
    for p in world.plates:
        for l in p.lines:
            if len(l):
                l.set_fields(elev_change_reason=np.full(len(l), ELEV_CHANGE_COLLISION, dtype=float))

    erosion.apply_erosion(world, years=1_000_000)

    after = plates.collect_all_elev_change_reason(world.plates)
    kept = np.mean(after == ELEV_CHANGE_COLLISION)
    assert kept > 0.8  # structural code is sticky against background wash
    # ...but some nodes did flip to a geomorphic code (the pass really ran on real terrain).
    assert np.any((after != ELEV_CHANGE_COLLISION) & (after != ELEV_CHANGE_NONE))


def test_apply_erosion_treats_lateral_magma_as_a_sticky_structural_code():
    """GitHub issue #205: ELEV_CHANGE_LATERAL_MAGMA is documented (elevation_lines.py,
    magma_transport.py) as getting the same erosion-override protection as every other
    structural code -- confirm `apply_erosion`'s own `prior_structural` mask actually includes
    it, not just the comment claiming so."""
    from app.elevation_lines import ELEV_CHANGE_LATERAL_MAGMA, ELEV_CHANGE_NONE

    world = generate_world(seed=21, num_plates=8)
    for p in world.plates:
        for l in p.lines:
            if len(l):
                l.set_fields(elev_change_reason=np.full(len(l), ELEV_CHANGE_LATERAL_MAGMA, dtype=float))

    erosion.apply_erosion(world, years=1_000_000)

    after = plates.collect_all_elev_change_reason(world.plates)
    kept = np.mean(after == ELEV_CHANGE_LATERAL_MAGMA)
    assert kept > 0.8
    assert np.any((after != ELEV_CHANGE_LATERAL_MAGMA) & (after != ELEV_CHANGE_NONE))


def test_earthquake_erosion_multiplier_bumps_near_the_epicentre_only():
    world = World(seed=0, plates=[])
    world.elapsed_years = 1_000_000
    ang = np.linspace(0.0, 0.5, 40)
    points = np.stack([np.cos(ang), np.sin(ang), np.zeros_like(ang)], axis=1)
    assert np.all(erosion._earthquake_erosion_multiplier(world, points) == 1.0)  # nothing yet

    world.earthquakes = [
        faults.Earthquake(
            earthquake_id=0, fault_id=0, plate_id=0, kind="reverse",
            epicenter_world=points[0], magnitude=7.5, slip_m=100.0, birth_years=1_000_000,
        )
    ]
    mult = erosion._earthquake_erosion_multiplier(world, points)
    assert mult[0] > 1.5  # right at a fresh M7.5 epicentre
    assert mult[-1] == 1.0  # far side of the arc: untouched
    assert np.all(np.diff(mult) <= 1e-9)


def test_earthquakes_increase_seismic_erosion_on_a_generated_world():
    base = generate_world(seed=21, num_plates=8)
    pts0, elev0, _, _, _, _ = erosion._gather_nodes(base)
    erosion.apply_erosion(base, years=1_000_000)
    elev_no_quake = _unloaded_elevation(base)
    # Highest land node that actually eroded on the plain run -- it has the height + slope the
    # seismic term keys off (and isn't submerged, where subaerial erosion is zeroed).
    eroded_land = np.where((elev0 > 800.0) & (elev0 - elev_no_quake > 1.0), elev0, -np.inf)
    target = int(np.argmax(eroded_land))
    assert np.isfinite(eroded_land[target])

    quaked = generate_world(seed=21, num_plates=8)
    quaked.elapsed_years = 0.0
    quaked.earthquakes = [
        faults.Earthquake(
            earthquake_id=0, fault_id=0, plate_id=0, kind="reverse",
            epicenter_world=pts0[target], magnitude=8.5, slip_m=200.0, birth_years=0.0,
        )
    ]
    erosion.apply_erosion(quaked, years=1_000_000)
    elev_quake = _unloaded_elevation(quaked)
    assert elev_quake[target] < elev_no_quake[target] - 1.0


def test_apply_erosion_respects_elevation_bounds():
    world = generate_world(seed=22, num_plates=8)
    erosion.apply_erosion(world, years=5_000_000)
    for plate in world.plates:
        for line in plate.lines:
            assert np.all(line.elevation >= MIN_ELEVATION_M - 1e-6)
            assert np.all(line.elevation <= MAX_ELEVATION_M + 1e-6)


def test_apply_erosion_noop_for_empty_world():
    world = World(seed=0, plates=[])
    erosion.apply_erosion(world, years=1_000_000)  # must not raise
    assert world.plates == []


# --- Symmetric coastal-leveling feedback (GitHub issue #122, "Speckled low-relief coastlines") ---


def _equator_lattice(half_span_deg: float, step_deg: float) -> np.ndarray:
    """A square lattice of unit vectors straddling (lat, lon) = (0, 0), for the coastal-
    feedback helpers (which only care about relative node geometry, not which plate owns
    what)."""
    grid = np.radians(np.arange(-half_span_deg, half_span_deg + 1e-9, step_deg))
    lat, lon = np.meshgrid(grid, grid, indexing="ij")
    return geometry.latlon_to_xyz(lat.ravel(), lon.ravel())


def test_coastal_openness_tracks_open_water_fraction():
    points = _equator_lattice(half_span_deg=3.0, step_deg=0.4)
    lon_deg = np.degrees(geometry.xyz_to_latlon(points)[1])
    is_ocean = lon_deg > 0.0  # a straight N-S coastline down the prime meridian

    openness = erosion._coastal_openness(points, is_ocean)

    deep_land = np.argmin(lon_deg)  # far to the west, ringed by land
    deep_ocean = np.argmax(lon_deg)  # far to the east, ringed by ocean
    coast = np.argmin(np.abs(lon_deg) + np.abs(np.degrees(geometry.xyz_to_latlon(points)[0])))
    assert openness[deep_land] < 0.05
    assert openness[deep_ocean] > 0.95
    assert 0.3 < openness[coast] < 0.7


def test_coastal_openness_zero_when_no_open_ocean():
    points = _equator_lattice(half_span_deg=2.0, step_deg=0.4)
    assert np.all(erosion._coastal_openness(points, np.zeros(len(points), dtype=bool)) == 0.0)


def test_leveling_datum_is_a_continuous_function_of_exposure():
    # Exposure pulls the datum below the waterline (wave-cut platform); shelter lifts it above
    # (marsh crest); a straight open coast lands near sea level.
    openness = np.array([0.6, 0.0, 0.3])  # exposed / fully enclosed / straight coast
    dist_to_land = np.full(3, 1.0)  # nowhere near a barrier
    datum = erosion.leveling_datum_m(0.0, openness, dist_to_land)
    assert np.isclose(datum[0], -erosion.LEVELING_PLATFORM_UNDERCUT_M)  # exposed -> undercut
    assert np.isclose(datum[1], erosion.LEVELING_MARSH_CREST_M)  # enclosed -> marsh crest
    assert datum[0] < datum[2] < datum[1]  # monotone decreasing in openness

    # A barrier candidate (hugs land, still faces open water) instead targets a bar crest.
    barrier = erosion.leveling_datum_m(
        0.0, np.array([0.4]), np.array([erosion.BARRIER_LANDWARD_RAD * 0.5])
    )
    assert np.isclose(barrier[0], erosion.BARRIER_CREST_M)


def test_coastal_leveling_grind_planes_anything_above_its_datum_including_submerged():
    # Per-node datum; nodes straddle sea level. Only the ones standing above their datum grind
    # -- a just-submerged shoal included (the old planation gate ignored everything below 0).
    elevation = np.array([10.0, 30.0, 100.0, -8.0, -30.0])
    datum = np.array([0.0, 0.0, 0.0, -20.0, 0.0])
    openness = np.full(5, 0.5)  # well above LEVELING_EXPOSURE_REF -> full "coastal" factor
    amt = erosion.coastal_leveling_grind(elevation, 0.0, openness, datum, dt_myr=0.05)
    assert amt[0] > amt[1] > 0.0  # both in band, the one closer to sea level planes faster
    assert amt[2] == 0.0  # above the band
    assert amt[3] > 0.0  # submerged shoal standing 12 m above its own datum still grinds
    assert amt[4] == 0.0  # below datum -- fill's job, not grind's

    landlocked = erosion.coastal_leveling_grind(elevation, 0.0, np.zeros(5), datum, dt_myr=0.05)
    assert np.all(landlocked == 0.0)  # openness 0 -> untouched

    # Never past its datum, however long the step.
    huge = erosion.coastal_leveling_grind(np.array([12.0]), 0.0, np.array([1.0]), np.array([0.0]), dt_myr=99.0)
    assert huge[0] == 12.0


def test_spread_coastal_leveling_conserves_mass_including_fallback():
    rng = np.random.default_rng(0)
    points = _equator_lattice(half_span_deg=3.0, step_deg=0.5)
    n = len(points)
    elevation = rng.uniform(-500.0, 200.0, n)
    openness = rng.uniform(0.0, 1.0, n)
    dist_to_land = rng.uniform(0.0, 0.02, n)
    datum = erosion.leveling_datum_m(0.0, openness, dist_to_land)
    source = np.zeros(n)
    source[rng.choice(n, 5, replace=False)] = rng.uniform(1.0, 10.0, 5)

    out, _ = erosion._spread_coastal_leveling(points, elevation, openness, dist_to_land, 0.0, datum, source, dt_myr=5.0)
    assert np.isclose(out.sum(), source.sum())

    # No band node below its datum anywhere -> every source keeps its own amount in place.
    no_sink, _ = erosion._spread_coastal_leveling(
        points, elevation, openness, dist_to_land, 0.0, elevation - 100.0, source, dt_myr=5.0
    )
    assert np.allclose(no_sink, source)


def test_spread_coastal_leveling_carries_a_tagged_share_and_conserves_volume_on_unequal_areas():
    # Issue #275: with per-node areas the pool is a volume, sinks hold their thickness capacity
    # times their own area, and a tagged (continental) sub-share follows the same transfers.
    rng = np.random.default_rng(1)
    points = _equator_lattice(half_span_deg=3.0, step_deg=0.5)
    n = len(points)
    elevation = rng.uniform(-30.0, 30.0, n)
    openness = rng.uniform(0.0, 1.0, n)
    dist_to_land = rng.uniform(0.0, 0.02, n)
    datum = erosion.leveling_datum_m(0.0, openness, dist_to_land)
    area = rng.uniform(1.0e9, 4.0e9, n)
    source = np.zeros(n)
    picked = rng.choice(n, 8, replace=False)
    source[picked] = rng.uniform(1.0, 50.0, 8) * area[picked]
    tagged = source * rng.uniform(0.0, 1.0, n)

    out, out_tagged = erosion._spread_coastal_leveling(
        points, elevation, openness, dist_to_land, 0.0, datum, source, dt_myr=5.0, area_m2=area, tagged_amount=tagged
    )
    assert np.isclose(out.sum(), source.sum())
    assert np.isclose(out_tagged.sum(), tagged.sum())
    assert np.all(out_tagged <= out * (1.0 + 1e-12))
    # No sink was filled past its own datum (capacity is thickness x area).
    received = np.where(source > 0.0, 0.0, out) / area
    assert np.all(received <= np.clip(datum - elevation, 0.0, None) + 1e-9)


def test_spread_coastal_leveling_prefers_sheltered_hollow_and_barrier_sinks():
    # A land source just east of three candidate sinks; only geometry/openness/depth differ.
    lat = np.radians(np.array([0.0, 0.0, 0.0, 0.0]))
    lon = np.radians(np.array([0.30, 0.10, -0.10, -0.30]))  # source, sheltered, deep-out-of-band, barrier
    points = geometry.latlon_to_xyz(lat, lon)
    elevation = np.array([20.0, -10.0, -2000.0, -1.0])
    openness = np.array([0.0, 0.1, 0.9, 0.4])  # sheltered bay / open abyss / barrier edge
    dist_to_land = np.array([1.0, 0.05, 0.05, 0.001])  # barrier node hugs the coast
    datum = erosion.leveling_datum_m(0.0, openness, dist_to_land)
    source = np.array([100.0, 0.0, 0.0, 0.0])

    out, _ = erosion._spread_coastal_leveling(points, elevation, openness, dist_to_land, 0.0, datum, source, dt_myr=5.0)
    assert out[1] > 0.0  # sheltered shallow water silts up
    assert out[3] > 0.0  # barrier candidate accretes despite facing open water
    assert out[2] < out[1] and out[2] < out[3]  # deep node is out of band -> gets ~nothing
    assert np.isclose(out.sum(), 100.0)


def test_spread_coastal_leveling_declumps_a_single_source_across_neighbours():
    # The whole point of the pass: one lump on one band node ends up spread across many.
    points = _equator_lattice(half_span_deg=2.0, step_deg=0.3)
    n = len(points)
    elevation = np.full(n, -5.0)  # a flat, uniformly just-submerged shelf
    openness = np.full(n, 0.2)
    dist_to_land = np.full(n, 1.0)
    datum = np.zeros(n)  # every node sits 5 m below its datum -> every node is a sink
    source = np.zeros(n)
    source[n // 2] = 500.0

    out, _ = erosion._spread_coastal_leveling(points, elevation, openness, dist_to_land, 0.0, datum, source, dt_myr=5.0)
    assert np.isclose(out.sum(), 500.0)
    assert np.count_nonzero(out > 1.0) >= 6  # the lump reached many neighbours
    assert out[n // 2] < 500.0  # and did not all stay on the source node


def test_spread_lake_sediment_leaves_dry_sinks_untouched():
    # A closed basin with no standing water yet (lake_depth all zero) keeps route_downstream's
    # old single-point-pile behavior -- nothing here is a real lake to spread across.
    lake_depth = np.zeros(4)
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1], [0, 0]])
    source = np.array([50.0, 0.0, 0.0, 0.0])

    out = erosion._spread_lake_sediment(lake_depth, neighbor_idx, source)
    assert out.tolist() == [50.0, 0.0, 0.0, 0.0]


def test_spread_lake_sediment_conserves_mass_and_reaches_every_member():
    # A 4-node flooded lake, all connected; sediment lands only at node 0 (its own catchment
    # sink) but should spread across all four members instead of piling entirely on node 0.
    lake_depth = np.array([5.0, 5.0, 5.0, 5.0])
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1], [0, 1]])
    source = np.zeros(4)
    source[0] = 40.0

    out = erosion._spread_lake_sediment(lake_depth, neighbor_idx, source)
    assert np.isclose(out.sum(), 40.0)
    assert np.all(out > 0.0)  # reached every member, not just the sink
    assert out[0] < 40.0  # and did not all stay piled on the sink


def test_spread_lake_sediment_fills_the_deepest_member_more_than_a_shallow_one():
    # Same lake, but node 2 is much deeper (a real valley within the basin) than the others --
    # it should receive more of the redistributed sediment, the "fill valleys more than hills"
    # behavior that flattens a basin's floor over repeated steps.
    lake_depth = np.array([2.0, 2.0, 40.0, 2.0])
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1], [0, 1]])
    source = np.zeros(4)
    source[0] = 40.0

    out = erosion._spread_lake_sediment(lake_depth, neighbor_idx, source)
    assert out[2] > out[1]
    assert out[2] > out[3]
    assert np.isclose(out.sum(), 40.0)


def test_spread_lake_sediment_keeps_two_separate_lakes_disjoint():
    # Two disconnected two-node ponds -- sediment landing in one must never leak into the other.
    lake_depth = np.array([3.0, 3.0, 3.0, 3.0])
    neighbor_idx = np.array([[1], [0], [3], [2]])
    source = np.array([10.0, 0.0, 20.0, 0.0])

    out = erosion._spread_lake_sediment(lake_depth, neighbor_idx, source)
    assert np.isclose(out[0] + out[1], 10.0)
    assert np.isclose(out[2] + out[3], 20.0)
    assert out[0] > 0.0 and out[1] > 0.0
    assert out[2] > 0.0 and out[3] > 0.0


def test_apply_erosion_coastal_feedback_keeps_a_generated_world_sane():
    world = generate_world(seed=23, num_plates=8)
    _, before, _, _, _, _ = erosion._gather_nodes(world)
    erosion.apply_erosion(world, years=2_000_000)
    _, after, _, _, _, _ = erosion._gather_nodes(world)
    assert np.all(np.isfinite(after))
    assert np.all(after >= MIN_ELEVATION_M - 1e-6) and np.all(after <= MAX_ELEVATION_M + 1e-6)
    # The feedback nudges the coast, it doesn't rewrite the whole map in one step.
    assert np.median(np.abs(after - before)) < 50.0


# --- Issue #275: conservative, provenance-aware redistribution -------------------------------

_LEDGER_SINKS = (
    "delaminated_lower_crust_m3",
    "deeply_subducted_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
    "overloaded_root_delaminated_m3",
)


def _tracked_continental_m3(world) -> float:
    """Live continental-derived material plus every declared sink -- constant across a step
    that only moves material."""
    from app import continental_ledger

    inventory = continental_ledger.inventories(world)
    return inventory["surface_continental_derived_m3"] + sum(inventory[k] for k in _LEDGER_SINKS)


def test_apply_erosion_closes_volume_and_continental_ledger_on_a_quad_world():
    from app import continental_ledger

    world = _generate_world(seed=3, num_plates=8, surface="quad")
    assert continental_ledger.inventories(world)["continental_sediment_on_oceanic_hosts_m3"] == 0.0
    # Two steps: the second has glaciers and spill/ice edges that point uphill.
    for _ in range(2):
        tracked_before = _tracked_continental_m3(world)
        result = erosion.apply_erosion(world, years=5_000_000)
        budget = result.budget
        assert budget["removed_m3"] > 0.0
        assert np.isclose(budget["deposited_m3"], budget["removed_m3"], rtol=1e-9)
        assert np.isclose(budget["continental_deposited_m3"] + budget["continental_discarded_m3"], budget["continental_removed_m3"], rtol=1e-9)
        assert np.isclose(_tracked_continental_m3(world), tracked_before, rtol=1e-12)
    # Continental-derived sediment shed off the coast stays identifiable on oceanic hosts.
    assert continental_ledger.inventories(world)["continental_sediment_on_oceanic_hosts_m3"] > 0.0


def test_apply_erosion_never_routes_more_than_a_column_can_give_up():
    from app import lithosphere

    world = _generate_world(seed=3, num_plates=8, surface="quad")
    headroom_m = 5.0
    for plate in world.plates:
        hc = np.full(plate.node_count(), lithosphere.MIN_CRUSTAL_THICKNESS_M + headroom_m)
        plate.set_fields_on_plate(crustal_thickness_m=hc, continental_material_m=hc)
    tracked_before = _tracked_continental_m3(world)
    result = erosion.apply_erosion(world, years=5_000_000)
    assert np.isclose(result.budget["deposited_m3"], result.budget["removed_m3"], rtol=1e-9)
    assert np.isclose(_tracked_continental_m3(world), tracked_before, rtol=1e-12)
    hc_after = np.concatenate([p.collect("crustal_thickness_m") for p in world.plates])
    assert hc_after.min() >= lithosphere.MIN_CRUSTAL_THICKNESS_M - 1e-6


def _profile_neighbors(n: int) -> np.ndarray:
    """Left/right neighbours along a 1-D profile (an end repeats itself)."""
    idx = np.arange(n)
    return np.stack([np.clip(idx - 1, 0, n - 1), np.clip(idx + 1, 0, n - 1)], axis=1)


def test_mass_wasting_runs_off_the_flank_and_settles_on_the_foreland():
    # Summit, two steep flank nodes, three flat foreland nodes, then the sea.
    elevation = np.array([6000.0, 4000.0, 2000.0, 300.0, 200.0, 100.0, -500.0])
    slope = np.array([0.05, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0])
    is_ocean = elevation < 0.0
    source = np.array([10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    tagged = 0.4 * source
    capacity = np.array([np.inf, np.inf, np.inf, 3.0, 2.0, 1.0, np.inf])
    headroom = np.full(len(elevation), 1.0e4)

    land, arrival, _, land_tagged, arrival_tagged, _ = erosion._route_mass_wasting(
        elevation, is_ocean, _profile_neighbors(len(elevation)), slope, source, tagged, capacity, headroom
    )
    # Nothing settles on the source or the steep flank; the flat foreland fills to capacity
    # and the rest reaches the sea.
    assert np.all(land[:3] == 0.0)
    assert np.allclose(land[3:6], [3.0, 2.0, 1.0])
    assert np.isclose(arrival[6], 4.0)
    assert np.isclose(land.sum() + arrival.sum(), source.sum())
    # The continental share travels with the debris.
    assert np.allclose(land_tagged, 0.4 * land)
    assert np.allclose(arrival_tagged, 0.4 * arrival)


def test_mass_wasting_runs_past_a_column_at_the_hc_cap():
    elevation = np.array([6000.0, 5000.0, 400.0, 300.0])
    slope = np.zeros(4)
    is_ocean = np.zeros(4, dtype=bool)
    source = np.array([5.0, 0.0, 0.0, 0.0])
    headroom = np.array([0.0, 0.0, 1.0e4, 1.0e4])  # summit and next node sit at the cap

    land, arrival, *_ = erosion._route_mass_wasting(
        elevation, is_ocean, _profile_neighbors(4), slope, source, source, np.full(4, np.inf), headroom
    )
    assert land[1] == 0.0  # flat, but no headroom: debris runs on
    assert np.isclose(land[2], 5.0)
    assert arrival.sum() == 0.0


def test_mass_wasting_debris_on_a_pit_or_from_a_low_source_stays_put():
    # Node 0 has no lower neighbour, so its own debris stays; node 1's debris reaches the pit at 2.
    elevation = np.array([100.0, 900.0, 50.0, 800.0])
    neighbors = np.array([[1, 1], [0, 2], [1, 3], [2, 2]])
    source = np.array([1.0, 2.0, 0.0, 0.0])
    land, arrival, *_ = erosion._route_mass_wasting(
        elevation, np.zeros(4, dtype=bool), neighbors, np.full(4, 0.05), source, source, np.zeros(4), np.full(4, 1.0e4)
    )
    assert np.allclose(land, [1.0, 0.0, 2.0, 0.0])
    assert arrival.sum() == 0.0


def test_mass_wasting_fills_a_pit_to_its_room_and_spills_the_rest_over_the_rim():
    # Summit 0 drains into pit 1, which has 2 units of room; its basin spills to 2, which
    # drains to the sea at 3.
    elevation = np.array([3000.0, 100.0, 500.0, -200.0])
    neighbors = np.array([[1, 1], [0, 2], [3, 1], [2, 2]])
    is_ocean = elevation < 0.0
    source = np.array([5.0, 0.0, 0.0, 0.0])
    capacity = np.array([np.inf, 2.0, 0.0, np.inf])
    spill = np.array([-1, 2, -1, -1])
    land, arrival, _, land_tagged, _, _ = erosion._route_mass_wasting(
        elevation, is_ocean, neighbors, np.full(4, 0.05), source, 0.5 * source, capacity, np.full(4, 1.0e4), spill_target=spill
    )
    assert np.allclose(land, [0.0, 2.0, 0.0, 0.0])
    assert np.isclose(arrival[3], 3.0)
    assert np.allclose(land_tagged, 0.5 * land)


def test_mass_wasting_hands_debris_that_reaches_the_ice_to_the_glacier():
    # A bare summit sheds onto a glaciated node (2); a glaciated summit's own debris (4) rides
    # its ice from the start. Neither settles on the ice-free flat in between.
    elevation = np.array([4000.0, 3000.0, 2000.0, 300.0, 3500.0])
    neighbors = np.array([[1, 1], [0, 2], [1, 3], [2, 2], [3, 3]])
    on_ice = np.array([False, False, True, False, True])
    source = np.array([6.0, 0.0, 0.0, 0.0, 4.0])
    land, arrival, ice, land_tagged, _, ice_tagged = erosion._route_mass_wasting(
        elevation, np.zeros(5, dtype=bool), neighbors, np.full(5, 0.05), source, 0.5 * source,
        np.full(5, np.inf), np.full(5, 1.0e4), on_ice=on_ice,
    )
    assert np.allclose(ice, [0.0, 0.0, 6.0, 0.0, 4.0])
    assert land.sum() == 0.0 and arrival.sum() == 0.0
    assert np.allclose(ice_tagged, 0.5 * ice)


def test_mass_wasting_inflow_settles_where_the_ice_drops_it_up_to_its_room():
    # The ice put 5 units down on flat margin node 1, which has room for 2; the rest runs on.
    elevation = np.array([900.0, 400.0, 300.0, -100.0])
    neighbors = _profile_neighbors(4)
    inflow = np.array([0.0, 5.0, 0.0, 0.0])
    land, arrival, ice, land_tagged, *_ = erosion._route_mass_wasting(
        elevation, elevation < 0.0, neighbors, np.zeros(4), np.zeros(4), np.zeros(4),
        np.array([np.inf, 2.0, 1.0, np.inf]), np.full(4, 1.0e4), inflow_vol=inflow, inflow_tagged=0.2 * inflow,
    )
    assert np.allclose(land, [0.0, 2.0, 1.0, 0.0])
    assert np.isclose(arrival[3], 2.0)
    assert ice.sum() == 0.0
    assert np.allclose(land_tagged, 0.2 * land)


def test_capped_lake_spread_keeps_debris_off_a_member_at_the_hc_cap():
    # One four-member lake; member 1 is at the cap. 10 units landed on member 0.
    lake_depth = np.array([50.0, 50.0, 50.0, 50.0])
    neighbor_idx = np.array([[1, 3], [0, 2], [1, 3], [2, 0]])
    source = np.array([10.0, 0.0, 0.0, 0.0])
    capacity = np.array([20.0, 0.0, 3.0, 20.0])
    out, tagged = erosion._spread_lake_sediment_capped(lake_depth, neighbor_idx, source, 0.3 * source, capacity)
    assert out[1] == 0.0
    assert np.all(out <= capacity + 1e-12)
    assert np.isclose(out[2], 3.0)  # filled to its room, the rest went to members with room
    assert np.isclose(out.sum(), 10.0)
    assert np.allclose(tagged, 0.3 * out)


def test_capped_lake_spread_leaves_what_no_member_can_hold_where_it_landed():
    lake_depth = np.array([50.0, 50.0, 0.0])
    neighbor_idx = np.array([[1, 2], [0, 2], [0, 1]])
    source = np.array([6.0, 0.0, 2.0])  # node 2 is dry land: untouched
    capacity = np.array([1.0, 2.0, 0.0])
    out, tagged = erosion._spread_lake_sediment_capped(lake_depth, neighbor_idx, source, source, capacity)
    assert np.allclose(out, [4.0, 2.0, 2.0])
    assert np.allclose(tagged, out)


def test_apply_erosion_routes_landslide_debris_into_basins_and_closes_the_ledger():
    world = _generate_world(seed=3, num_plates=8, surface="quad")
    world.seismic_erosion_multiplier = 4.0
    tracked_before = _tracked_continental_m3(world)
    budget = erosion.apply_erosion(world, years=5_000_000).budget
    assert budget["mass_wasting_removed_m3"] > 0.0
    assert np.isclose(
        budget["mass_wasting_foreland_m3"] + budget["mass_wasting_marine_m3"], budget["mass_wasting_removed_m3"], rtol=1e-9
    )
    assert np.isclose(budget["deposited_m3"], budget["removed_m3"], rtol=1e-9)
    assert np.isclose(_tracked_continental_m3(world), tracked_before, rtol=1e-12)


def test_landslides_lower_a_frozen_summit_that_water_cannot_drain(monkeypatch):
    # Hydrology gives frozen land no flow target, so water-routed sediment eroded there settles
    # straight back. Landslide debris moves by gravity and must still leave the summit.
    import dataclasses

    from app import hydrology

    real = hydrology.compute_hydrology

    def all_frozen(*args, **kwargs):
        hydro = real(*args, **kwargs)
        return dataclasses.replace(hydro, flow_target=np.where(hydro.is_ocean, hydro.flow_target, -1))

    monkeypatch.setattr(hydrology, "compute_hydrology", all_frozen)
    world = _generate_world(seed=3, num_plates=8, surface="quad")
    world.glacier_erosion_multiplier = 0.0
    world.wind_erosion_multiplier = 0.0
    _, elevation, _, _, _, order = erosion._gather_nodes(world)
    summit = int(np.argmax(elevation))
    hc_before = plates.collect_all_crustal_thickness(order)[summit]
    erosion.apply_erosion(world, years=5_000_000)
    assert plates.collect_all_crustal_thickness(order)[summit] < hc_before - 1.0


def test_ocean_deposition_knob_below_one_declares_the_continental_sediment_it_withholds():
    world = _generate_world(seed=3, num_plates=8, surface="quad")
    world.ocean_deposition_multiplier = 0.5
    tracked_before = _tracked_continental_m3(world)
    result = erosion.apply_erosion(world, years=5_000_000)
    assert result.budget["ocean_deposition_knob_m3"] < 0.0
    assert world.continental_material_ledger["discarded_marine_sediment_m3"] == result.budget["continental_discarded_m3"] > 0.0
    assert np.isclose(_tracked_continental_m3(world), tracked_before, rtol=1e-12)


def test_flatten_is_a_volume_conserving_downhill_exchange():
    from types import SimpleNamespace

    rng = np.random.default_rng(2)
    n, k = 40, 6
    neighbor_idx = np.array([rng.choice(np.delete(np.arange(n), i), k, replace=False) for i in range(n)])
    hydro = SimpleNamespace(elevation=rng.uniform(0.0, 3000.0, n), neighbor_idx=neighbor_idx)
    ice_factor = rng.uniform(0.0, 2.0, n)
    area = rng.uniform(1.0e9, 4.0e9, n)
    removable = rng.uniform(0.0, 50.0, n)

    send = erosion._flatten(hydro, ice_factor, years=5_000_000, removable_m=removable)
    assert np.all(send >= 0.0)
    # Only ever sends downhill, and never more than the column's headroom.
    assert np.all(send[hydro.elevation[:, None] <= hydro.elevation[neighbor_idx]] == 0.0)
    assert np.all(send.sum(axis=1) <= removable + 1e-9)
    received = np.zeros(n)
    np.add.at(received, neighbor_idx.ravel(), (send * area[:, None]).ravel())
    assert np.isclose(received.sum(), (send.sum(axis=1) * area).sum())


def test_apply_erosion_passes_seasonal_amplitude_to_hydrology(monkeypatch):
    from app import hydrology

    world = generate_world(seed=20, num_plates=8)
    seen = {}
    real = hydrology.compute_hydrology

    def spy(*args, **kwargs):
        seen["amplitude"] = kwargs.get("seasonal_amplitude_at_nodes")
        return real(*args, **kwargs)

    monkeypatch.setattr(hydrology, "compute_hydrology", spy)
    erosion.apply_erosion(world, years=1_000_000)
    amplitude = seen["amplitude"]
    assert amplitude is not None and np.all(np.isfinite(amplitude)) and np.all(amplitude >= 0.0)
    assert amplitude.max() > 10.0  # a real high-latitude / interior swing somewhere


# --- Issue #275 phase 3: reversible ice loading ----------------------------------------------


def _force_glacier_depth(monkeypatch, depth_for_step):
    """Make compute_hydrology return its real fields but with `glacier_depth` replaced by
    `depth_for_step(n)` for the current call."""
    import dataclasses

    from app import hydrology

    real = hydrology.compute_hydrology

    def forced(*args, **kwargs):
        hydro = real(*args, **kwargs)
        return dataclasses.replace(hydro, glacier_depth=depth_for_step(len(hydro.glacier_depth)))

    monkeypatch.setattr(hydrology, "compute_hydrology", forced)


def test_ice_load_depresses_the_surface_and_meltback_rebounds_it_exactly(monkeypatch):
    from app import lithosphere

    world = _generate_world(seed=3, num_plates=8, surface="quad")
    _, elevation, _, _, _, plates_in_order = erosion._gather_nodes(world)
    hc = np.concatenate([p.collect("crustal_thickness_m") for p in plates_in_order])
    # Thick ice on high, comfortably-dry continental columns only.
    loaded = (elevation > world.sea_level_m + 1500.0) & (hc > 0.0)
    assert loaded.sum() > 10
    ice_m = 2000.0
    schedule = iter([np.where(loaded, ice_m, 0.0), np.zeros(len(loaded))])
    _force_glacier_depth(monkeypatch, lambda n: next(schedule))

    def deflection():
        return np.concatenate([p.collect("ice_load_deflection_m") for p in world.plates])

    # Loading step: the surface sinks by the Airy response to the ice, on top of erosion.
    before = erosion._gather_nodes(world)[1]
    loading = erosion.apply_erosion(world, years=1_000_000)
    after_load = erosion._gather_nodes(world)[1]
    expected_w = -ice_m * lithosphere.RHO_ICE / lithosphere.RHO_ASTHENOSPHERE
    assert np.allclose(deflection()[loaded], expected_w, rtol=0.02)
    assert np.all(deflection()[~loaded] == 0.0)
    assert np.allclose(after_load - before, loading.net_elevation_change_m + loading.ice_deflection_change_m, atol=1e-6)
    assert np.allclose(loading.ice_load_pa[loaded], ice_m * lithosphere.RHO_ICE * lithosphere.GRAVITY_M_S2)
    assert np.all(loading.ice_load_change_pa[loaded] > 0.0)

    # Meltback step: the deflection comes off exactly, through the same explicit path.
    applied = deflection().copy()
    unloading = erosion.apply_erosion(world, years=1_000_000)
    after_melt = erosion._gather_nodes(world)[1]
    assert np.all(deflection() == 0.0)
    assert np.allclose(unloading.ice_deflection_change_m, -applied)
    assert np.allclose(after_melt - after_load, unloading.net_elevation_change_m - applied, atol=1e-6)
    assert np.all(unloading.ice_load_change_pa[loaded] < 0.0)


def test_ice_load_near_the_elevation_floor_stores_only_the_applied_deflection():
    # A legacy (no-Hc) column 100 m above MIN_ELEVATION_M under 2 km of grounded ice: the raw
    # dry-land response (~-564 m) would push it through the floor. The stored deflection must
    # be the 100 m the clip actually let through, so meltback lands back on the original bed
    # rather than ~464 m above it.
    from app import lithosphere

    bed = np.array([MIN_ELEVATION_M + 100.0, 500.0])
    zeros = np.zeros(2)
    rho = np.full(2, lithosphere.RHO_CONTINENTAL_CRUST)
    sea_level = MIN_ELEVATION_M - 1000.0  # both columns dry, so the full ice column loads
    ice = np.full(2, 2000.0)

    loaded, deflection, load = erosion._apply_ice_load(bed, zeros, ice, zeros, zeros, rho, sea_level)
    raw = -ice * lithosphere.RHO_ICE / lithosphere.RHO_ASTHENOSPHERE
    assert np.all(load > 0.0)
    assert loaded[0] == MIN_ELEVATION_M and np.isclose(deflection[0], -100.0)
    assert np.isclose(deflection[1], raw[1]) and np.isclose(loaded[1], 500.0 + raw[1])
    assert np.allclose(loaded - bed, deflection)

    melted, after, _ = erosion._apply_ice_load(loaded, deflection, zeros, zeros, zeros, rho, sea_level)
    assert np.all(after == 0.0)
    assert np.allclose(melted, bed)


# --- Issue #288: depositional Hc-cap overflow is carried on, not clipped ----------------------


def _carry(elevation, neighbors, room, excess, *, area=None, is_ocean=None, on_ice=None, lake_depth=None,
           spill=None, ice_target=None, scour_limit=None, scour_material=None, points=None, **cover):
    """`_carry_overflow` on a hand-built profile, with an all-river excess matrix built from
    `excess` (volume per node) and a 0.4 continental share."""
    n = len(elevation)
    columns = len(erosion.OVERFLOW_PATHWAYS)
    load = np.zeros((n, columns + 1))
    load[:, 0] = excess
    load[:, -1] = 0.4 * np.asarray(excess, dtype=float)
    if points is None:
        # Nodes ~1 km apart along the equator, so every marine/search range reaches them all.
        lon = np.arange(n) * (1.0 / 6371.0)
        points = np.stack([np.cos(lon), np.sin(lon), np.zeros(n)], axis=1)
    return erosion._carry_overflow(
        points,
        np.asarray(elevation, dtype=float),
        np.zeros(n, dtype=bool) if is_ocean is None else is_ocean,
        np.zeros(n, dtype=bool) if on_ice is None else on_ice,
        np.zeros(n) if lake_depth is None else lake_depth,
        neighbors,
        np.full(n, -1) if spill is None else spill,
        np.full(n, -1) if ice_target is None else ice_target,
        np.ones(n) if area is None else area,
        np.asarray(room, dtype=float),
        load,
        np.zeros(n) if scour_limit is None else scour_limit,
        np.zeros(n) if scour_material is None else scour_material,
        **cover,
    ), load


def test_overflow_runs_on_through_several_saturated_receivers():
    # Node 0 overflows by 6; nodes 1 and 2 are already full; 3 and 4 have 2 and 1 of room;
    # the rest reaches the sea at 5.
    elevation = np.array([900.0, 800.0, 700.0, 600.0, 500.0, -100.0])
    room = np.array([0.0, 0.0, 0.0, 2.0, 1.0, 1.0e9])
    carry, load = _carry(elevation, _profile_neighbors(6), room, [6.0, 0, 0, 0, 0, 0], is_ocean=elevation < 0)
    volume = carry.placed[:, :-1].sum(axis=1)
    assert np.allclose(volume, [0.0, 0.0, 0.0, 2.0, 1.0, 3.0])
    assert carry.terminal.sum() == 0.0
    # The continental share and the pathway column travel with the load.
    assert np.allclose(carry.placed[:, -1], 0.4 * volume)
    assert np.allclose(carry.placed[:, 0], volume)


def test_overflow_conserves_volume_on_unequal_area_cells():
    # Thickness room is the same everywhere (1 m), so the large cell takes 4x the volume.
    elevation = np.array([900.0, 800.0, 700.0, 600.0])
    area = np.array([1.0e9, 4.0e9, 0.5e9, 2.0e9])
    room = 1.0 * area * np.array([0.0, 1.0, 1.0, 1.0])
    excess = np.array([5.0e9, 0.0, 0.0, 0.0])
    carry, load = _carry(elevation, _profile_neighbors(4), room, excess, area=area)
    volume = carry.placed[:, :-1].sum(axis=1)
    assert np.allclose(volume / area, [0.0, 1.0, 1.0, 0.25])
    assert np.isclose(volume.sum() + carry.terminal[:-1].sum(), excess.sum())
    assert np.isclose(carry.placed[:, -1].sum(), 0.4 * excess.sum())


def test_overflow_spreads_across_a_lake_and_leaves_over_its_spill_point():
    # Pit 1 is a three-member lake (1, 2, 3) with 1 unit of room per member; its outlet spills
    # to 4, which drains to the sea at 5.
    elevation = np.array([900.0, 100.0, 120.0, 130.0, 400.0, -100.0])
    neighbors = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [5, 5], [4, 4]])
    lake_depth = np.array([0.0, 30.0, 20.0, 10.0, 0.0, 0.0])
    spill = np.array([-1, 4, 4, 4, -1, -1])
    room = np.array([0.0, 1.0, 1.0, 1.0, 0.0, 1.0e9])
    carry, _ = _carry(elevation, neighbors, room, [5.0, 0, 0, 0, 0, 0], is_ocean=elevation < 0, lake_depth=lake_depth, spill=spill)
    volume = carry.placed[:, :-1].sum(axis=1)
    assert np.allclose(volume[1:4], 1.0)
    assert np.isclose(volume[5], 2.0)
    assert np.isclose(volume.sum(), 5.0)


def test_overflow_at_sea_spreads_onto_lower_ocean_nodes_with_room():
    elevation = np.array([-100.0, -200.0, -300.0, -400.0])
    room = np.array([0.0, 1.0, 1.0, 10.0])
    carry, _ = _carry(elevation, _profile_neighbors(4), room, [4.0, 0, 0, 0], is_ocean=np.ones(4, dtype=bool))
    volume = carry.placed[:, :-1].sum(axis=1)
    assert volume[0] == 0.0
    assert np.all(volume <= room + 1e-9)
    assert np.isclose(volume.sum(), 4.0)


def test_overflow_rides_the_ice_and_its_scour_deepens_the_bed():
    # Node 0's overflow lands on the ice (1, 2), which carries it to its margin at 3. Each ice
    # node crossed gives up rock, up to its scour limit, all continental here.
    elevation = np.array([3000.0, 2000.0, 1500.0, 500.0, 400.0])
    on_ice = np.array([False, True, True, False, False])
    ice_target = np.array([-1, 2, 3, -1, -1])
    room = np.array([0.0, 50.0, 50.0, 100.0, 100.0])
    scour_limit = np.array([0.0, 0.2, 0.3, 0.0, 0.0])
    carry, _ = _carry(
        elevation, _profile_neighbors(5), room, [10.0, 0, 0, 0, 0], on_ice=on_ice, ice_target=ice_target,
        scour_limit=scour_limit, scour_material=np.full(5, 1.0e3),
    )
    volume = carry.placed[:, :-1].sum(axis=1)
    assert volume[1] == 0.0 and volume[2] == 0.0  # nothing settles under the ice
    assert np.allclose(carry.scour_m, [0.0, 0.2, 0.3, 0.0, 0.0])
    assert np.isclose(volume[3], 10.5)
    # The scoured rock rides as glacial load, carrying its own continental share.
    assert np.isclose(carry.placed[:, erosion.OVERFLOW_GLACIAL].sum(), 0.5)
    assert np.isclose(carry.placed[:, -1].sum(), 4.0 + 0.5)


def test_overflow_ice_scour_takes_mobile_cover_before_substrate():
    # As above, but node 1 has 0.15 m of cover (a third continental) over all-continental rock,
    # and node 2 none: node 1's scour is mostly cover, tagged at the cover's own fraction.
    elevation = np.array([3000.0, 2000.0, 1500.0, 500.0, 400.0])
    on_ice = np.array([False, True, True, False, False])
    ice_target = np.array([-1, 2, 3, -1, -1])
    room = np.array([0.0, 50.0, 50.0, 100.0, 100.0])
    scour_limit = np.array([0.0, 0.2, 0.3, 0.0, 0.0])
    cover = np.array([0.0, 0.15, 0.0, 0.0, 0.0])
    carry, _ = _carry(
        elevation, _profile_neighbors(5), room, [10.0, 0, 0, 0, 0], on_ice=on_ice, ice_target=ice_target,
        scour_limit=scour_limit, scour_material=np.full(5, 1.0e3), scour_cover_m=cover, scour_cover_material_m=cover / 3.0,
    )
    assert np.allclose(carry.scour_m, [0.0, 0.2, 0.3, 0.0, 0.0])
    assert np.allclose(carry.scour_cover_m, [0.0, 0.15, 0.0, 0.0, 0.0])
    assert np.allclose(carry.scour_tagged_m, [0.0, 0.05 + 0.05, 0.3, 0.0, 0.0])
    assert np.isclose(carry.placed[:, -1].sum(), 4.0 + 0.4)


def test_overflow_caught_in_an_ice_loop_drains_off_by_gravity():
    # The ice on 1 and 2 flows in a loop (as ice_flow_target can on a capped ice cap); the load
    # still leaves it downhill, to the margin at 3.
    elevation = np.array([3000.0, 2000.0, 1900.0, 500.0])
    on_ice = np.array([False, True, True, False])
    ice_target = np.array([-1, 2, 1, -1])
    room = np.array([0.0, 0.0, 0.0, 100.0])
    carry, _ = _carry(elevation, _profile_neighbors(4), room, [4.0, 0, 0, 0], on_ice=on_ice, ice_target=ice_target)
    assert np.isclose(carry.placed[3, :-1].sum(), 4.0)
    assert carry.stuck_m3["hop_limit"] == 0.0 and carry.terminal.sum() == 0.0


def test_overflow_in_a_full_closed_basin_spills_over_its_lowest_saddle():
    # Node 0's overflow drains into full pit 1, which has no spill target. The basin fills and
    # overflows its saddle (2, also full) onto 3 below it -- not onto 4, higher up the far side,
    # nor back up the flank it came down.
    elevation = np.array([900.0, 100.0, 500.0, 200.0, 800.0])
    room = np.array([0.0, 0.0, 0.0, 10.0, 10.0])
    carry, _ = _carry(elevation, _profile_neighbors(5), room, [3.0, 0, 0, 0, 0])
    assert np.isclose(carry.placed[3, :-1].sum(), 3.0)
    assert carry.stuck_m3["pit"] == 3.0
    assert np.isclose(carry.basin_fill_m3, 3.0) and carry.search_m3 == 0.0


def test_overflow_with_no_receiver_in_reach_becomes_the_terminal_remainder():
    # Every node full, a saturated pit with no outlet: nothing can take it.
    elevation = np.array([900.0, 100.0, 500.0])
    neighbors = np.array([[1, 1], [0, 2], [1, 1]])
    carry, load = _carry(elevation, neighbors, np.zeros(3), [3.0, 0.0, 0.0])
    assert carry.placed.sum() == 0.0
    assert np.allclose(carry.terminal, load.sum(axis=0))


def test_overflow_no_basin_fill_can_reach_goes_to_the_nearest_node_with_room():
    # Same saturated pit, but the only room is on node 2, which no neighbour edge reaches.
    elevation = np.array([900.0, 100.0, 500.0])
    neighbors = np.array([[1, 1], [0, 0], [2, 2]])
    carry, _ = _carry(elevation, neighbors, np.array([0.0, 0.0, 5.0]), [3.0, 0.0, 0.0])
    assert np.isclose(carry.placed[2, :-1].sum(), 3.0)
    assert carry.basin_fill_m3 == 0.0 and np.isclose(carry.search_m3, 3.0)
    assert carry.terminal.sum() == 0.0


def _saturate_lowlands(world, quantile=0.5):
    """Push every land column below the `quantile` land elevation to the Hc cap -- the
    foreland/basin receivers erosion deposits on."""
    from app import lithosphere

    sea = world.sea_level_m
    land_elev = np.concatenate([p.collect("elevation") for p in world.plates])
    cutoff = np.quantile(land_elev[land_elev > sea], quantile)
    for plate in world.plates:
        elev = plate.collect("elevation")
        hc = plate.collect("crustal_thickness_m")
        material = plate.collect("continental_material_m")
        full = (elev > sea) & (elev < cutoff)
        new_hc = np.where(full, lithosphere.MAX_CRUSTAL_THICKNESS_M, hc)
        plate.set_fields_on_plate(crustal_thickness_m=new_hc, continental_material_m=np.where(full, new_hc, material))


def test_apply_erosion_carries_cap_overflow_on_instead_of_booking_numerical_loss():
    from app import continental_ledger, lithosphere

    world = _generate_world(seed=3, num_plates=8, surface="quad")
    continental_ledger.ensure_initialized(world)
    _saturate_lowlands(world)
    for _ in range(2):
        unplaced_before = world.continental_material_ledger["numerical_unplaced_m3"]
        tracked_before = _tracked_continental_m3(world)
        budget = erosion.apply_erosion(world, years=5_000_000).budget
        assert budget["hc_cap_overflow_m3"] > 0.0 and budget["hc_cap_overflow_nodes"] > 0
        assert budget["hc_cap_overflow_placed_m3"] > 0.0
        # Round-off at most (columns land on the cap to ~1e-11 m), never the overflow itself.
        assert world.continental_material_ledger["numerical_unplaced_m3"] - unplaced_before < 1e-9 * budget["continental_overflow_m3"]
        assert np.isclose(_tracked_continental_m3(world), tracked_before, rtol=1e-12)
        # Volume closes once the terminal reservoir is counted.
        terminal_sediment = budget["hc_cap_overflow_terminal_m3"] - budget["hc_cap_overflow_lake_silt_terminal_m3"]
        assert np.isclose(budget["removed_m3"], budget["deposited_m3"] + terminal_sediment, rtol=1e-9)
        assert np.isclose(
            budget["continental_removed_m3"],
            budget["continental_deposited_m3"] + budget["continental_discarded_m3"] + budget["continental_overflow_terminal_m3"],
            rtol=1e-9,
        )
        # The per-pathway diagnostics partition the totals.
        pathways = erosion.OVERFLOW_PATHWAYS
        assert np.isclose(sum(budget[f"hc_cap_overflow_{p}_m3"] for p in pathways), budget["hc_cap_overflow_m3"], rtol=1e-9)
        assert np.isclose(
            sum(budget[f"hc_cap_overflow_{p}_placed_m3"] + budget[f"hc_cap_overflow_{p}_terminal_m3"] for p in pathways),
            budget["hc_cap_overflow_m3"] + budget["hc_cap_overflow_ice_scour_m3"],
            rtol=1e-9,
        )
        hc = np.concatenate([p.collect("crustal_thickness_m") for p in world.plates])
        assert hc.max() <= lithosphere.MAX_CRUSTAL_THICKNESS_M + 1e-6


def test_apply_erosion_books_unplaceable_overflow_as_overloaded_root_delamination(monkeypatch):
    from app import continental_ledger

    def nowhere(*args, **kwargs):
        excess = args[10]
        return erosion.OverflowCarry(np.zeros_like(excess), excess.sum(axis=0), np.zeros(len(excess)), np.zeros(len(excess)))

    monkeypatch.setattr(erosion, "_carry_overflow", nowhere)
    world = _generate_world(seed=3, num_plates=8, surface="quad")
    continental_ledger.ensure_initialized(world)
    _saturate_lowlands(world)
    unplaced_before = world.continental_material_ledger["numerical_unplaced_m3"]
    tracked_before = _tracked_continental_m3(world)
    budget = erosion.apply_erosion(world, years=5_000_000).budget
    assert budget["continental_overflow_terminal_m3"] > 0.0
    assert world.continental_material_ledger["overloaded_root_delaminated_m3"] == budget["continental_overflow_terminal_m3"]
    # Round-off at most (columns land on the cap to ~1e-11 m).
    assert world.continental_material_ledger["numerical_unplaced_m3"] - unplaced_before < 1e-9 * budget["continental_overflow_m3"]
    assert np.isclose(_tracked_continental_m3(world), tracked_before, rtol=1e-12)


def test_glacial_erosion_carves_channels_rivers_can_inherit():
    world = _generate_world(seed=3, num_plates=8, surface="quad")
    world.rain_erosion_multiplier = 0.0
    world.river_erosion_multiplier = 0.0
    for plate in world.plates:
        plate.set_fields_on_plate(glacier_depth=np.where(plate.collect("elevation") > world.sea_level_m, 500.0, 0.0))
    _, _, before, _, _, _ = erosion._gather_nodes(world)
    erosion.apply_erosion(world, years=5_000_000)
    _, _, after, _, _, _ = erosion._gather_nodes(world)
    # No river erosion at all, so every metre of new channel is a glacial trough.
    assert np.max(after - before) > 1.0


def _set_channel_state(world, depth_m, uplift_m):
    """Every node gets `depth_m` of channel and a reference elevation `uplift_m` below its
    current elevation, as if it had risen that much since last step's erosion."""
    for p in world.plates:
        for line in p.lines:
            if len(line):
                line.set_fields(
                    channel_depth=np.full(len(line), depth_m),
                    channel_reference_elevation_m=line.elevation - uplift_m,
                )


def test_apply_erosion_records_the_reference_elevation_for_next_steps_uplift():
    world = generate_world(seed=21, num_plates=8)
    erosion.apply_erosion(world, years=1_000_000)
    _, elevation, _, _, _, plates_in_order = erosion._gather_nodes(world)
    np.testing.assert_allclose(plates.collect_all_channel_reference_elevation(plates_in_order), elevation)


def test_uplift_since_last_step_wears_channels_down_by_the_rise():
    # Two identical worlds, one of which rose 300 m since its last erosion pass. Erosion
    # itself doesn't read the reference, so the only difference is the fade.
    still = generate_world(seed=21, num_plates=8)
    risen = generate_world(seed=21, num_plates=8)
    _set_channel_state(still, 500.0, 0.0)
    _set_channel_state(risen, 500.0, 300.0)
    erosion.apply_erosion(still, years=1_000_000)
    erosion.apply_erosion(risen, years=1_000_000)

    still_depth = plates.collect_all_channel_depth(still.plates)
    risen_depth = plates.collect_all_channel_depth(risen.plates)
    land = still_depth > 0.0
    assert land.any()
    np.testing.assert_allclose(risen_depth[land], np.clip(still_depth[land] - 300.0, 0.0, None))


def test_deposition_fills_channels_back_in():
    # A fresh world has no ice, so wherever no river erodes nothing carves the channel. There
    # the channel loses exactly what settled into it: deposited sediment plus lake silt.
    world = generate_world(seed=21, num_plates=8)
    _set_channel_state(world, 500.0, 0.0)
    result = erosion.apply_erosion(world, years=1_000_000)

    depth = plates.collect_all_channel_depth(world.plates)
    filled = np.clip(result.sediment_deposited, 0.0, None) + world.hydrology_cache.silt_deposited
    quiet = ~world.hydrology_cache.is_ocean & (result.river == 0.0) & (filled > 0.0)
    assert quiet.any()
    np.testing.assert_allclose(depth[quiet], np.clip(500.0 - filled[quiet], 0.0, None))

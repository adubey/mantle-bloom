import numpy as np

from app import erosion, faults, geometry, plates
from app.elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M
from app.world import World, generate_world


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

    _, elev_after, _, _, _, plates_after = erosion._gather_nodes(world)
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
    _, elev_no_quake, _, _, _, _ = erosion._gather_nodes(base)
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
    _, elev_quake, _, _, _, _ = erosion._gather_nodes(quaked)
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


# --- Symmetric coastal-leveling feedback (docs/TODO.md "Speckled low-relief coastlines") ---


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

    out = erosion._spread_coastal_leveling(points, elevation, openness, dist_to_land, 0.0, datum, source, dt_myr=5.0)
    assert np.isclose(out.sum(), source.sum())

    # No band node below its datum anywhere -> every source keeps its own amount in place.
    no_sink = erosion._spread_coastal_leveling(
        points, elevation, openness, dist_to_land, 0.0, elevation - 100.0, source, dt_myr=5.0
    )
    assert np.allclose(no_sink, source)


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

    out = erosion._spread_coastal_leveling(points, elevation, openness, dist_to_land, 0.0, datum, source, dt_myr=5.0)
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

    out = erosion._spread_coastal_leveling(points, elevation, openness, dist_to_land, 0.0, datum, source, dt_myr=5.0)
    assert np.isclose(out.sum(), 500.0)
    assert np.count_nonzero(out > 1.0) >= 6  # the lump reached many neighbours
    assert out[n // 2] < 500.0  # and did not all stay on the source node


def test_apply_erosion_coastal_feedback_keeps_a_generated_world_sane():
    world = generate_world(seed=23, num_plates=8)
    _, before, _, _, _, _ = erosion._gather_nodes(world)
    erosion.apply_erosion(world, years=2_000_000)
    _, after, _, _, _, _ = erosion._gather_nodes(world)
    assert np.all(np.isfinite(after))
    assert np.all(after >= MIN_ELEVATION_M - 1e-6) and np.all(after <= MAX_ELEVATION_M + 1e-6)
    # The feedback nudges the coast, it doesn't rewrite the whole map in one step.
    assert np.median(np.abs(after - before)) < 50.0

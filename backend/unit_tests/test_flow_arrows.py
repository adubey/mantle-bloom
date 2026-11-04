import numpy as np
import pytest
from app import flow_arrows, geometry

# A coarse 2-degree grid keeps each trace well under a second; the module itself is
# resolution-independent (spacing rules are in degrees).
H, W = 90, 180
LAT_DEG = 90.0 - (np.arange(H) + 0.5) * (180.0 / H)
LON_DEG = -180.0 + (np.arange(W) + 0.5) * (360.0 / W)
LAT, LON = np.meshgrid(np.radians(LAT_DEG), np.radians(LON_DEG), indexing="ij")


def _eastward():
    return np.ones((H, W)), np.zeros((H, W))


def _vortex_about_x():
    """Solid-body rotation about the x axis, as (east, north) components: every streamline is
    a circle around (1, 0, 0), so a trace keeps turning."""
    xyz = np.stack([np.cos(LAT) * np.cos(LON), np.cos(LAT) * np.sin(LON), np.sin(LAT)], axis=-1)
    velocity = np.cross(np.array([1.0, 0.0, 0.0]), xyz)
    east = np.stack([-np.sin(LON), np.cos(LON), np.zeros_like(LON)], axis=-1)
    north = np.cross(xyz, east)
    return np.sum(velocity * east, axis=-1), np.sum(velocity * north, axis=-1)


def _arc_length_deg(points):
    return float(np.degrees(np.sum(geometry.angular_distance(points[:-1], points[1:]))))


def _turning_deg(points):
    """(net signed, unsigned) turning along a path, measured about the local vertical so a
    great circle turns 0 -- the same quantity the tracer limits."""
    d = np.diff(points, axis=0)
    d -= np.sum(d * points[:-1], axis=1, keepdims=True) * points[:-1]
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    angle = np.arctan2(np.sum(np.cross(d[:-1], d[1:]) * points[1:-1], axis=1), np.sum(d[:-1] * d[1:], axis=1))
    return float(np.degrees(angle.sum())), float(np.degrees(np.abs(angle).sum()))


def test_no_flow_draws_no_arrows():
    zeros = np.zeros((H, W))
    assert flow_arrows.trace_flow_arrows(zeros, zeros, np.ones((H, W), dtype=bool)) == []


def test_uniform_eastward_flow_draws_separated_eastward_arrows():
    u, v = _eastward()
    params = flow_arrows.FlowArrowParams()
    arrows = flow_arrows.trace_flow_arrows(u, v, np.ones((H, W), dtype=bool), params)
    assert len(arrows) > 10
    for arrow in arrows:
        lat, lon = geometry.xyz_to_latlon(arrow.points_xyz)
        assert np.ptp(lat) < np.radians(1.0)
        assert np.all(np.diff(np.unwrap(lon)) > 0), "arrow should run tail to tip with the flow"
        assert params.min_length_deg - 1e-6 <= _arc_length_deg(arrow.points_xyz) <= params.max_length_deg + 1.0
    # A growing arrow stops `stop_deg` short of any earlier arrow; the one-degree occupancy
    # grid and the spline's sampling allow about a degree of slack.
    for i, a in enumerate(arrows):
        for b in arrows[i + 1:]:
            closest = np.min(geometry.angular_distance(a.points_xyz[:, None, :], b.points_xyz[None, :, :]))
            assert np.degrees(closest) > params.stop_deg - 1.5


def test_arrows_stay_inside_the_valid_region():
    # Currents trace only over ocean: a half-planet "ocean" must hold every arrow point.
    u, v = _vortex_about_x()
    valid = np.abs(LON) < np.radians(90.0)
    arrows = flow_arrows.trace_flow_arrows(u, v, valid)
    assert arrows
    for arrow in arrows:
        _, lon = geometry.xyz_to_latlon(arrow.points_xyz)
        assert np.all(np.abs(lon) < np.radians(91.0))


def test_turn_limit_splits_a_vortex_into_arcs():
    # A closed circulation must not become one ring-shaped arrow: the turn limit cuts it into
    # arcs, which is how atlas maps draw a gyre.
    u, v = _vortex_about_x()
    # Lengths well past the real cap (42 degrees) so turning, not length, ends the arcs. Under
    # 180 degrees, since the vortex's widest streamline is a great circle that never turns.
    params = flow_arrows.FlowArrowParams(max_length_deg=150.0, max_turn_deg=120.0)
    arrows = flow_arrows.trace_flow_arrows(u, v, np.ones((H, W), dtype=bool), params)
    assert arrows
    for arrow in arrows:
        net, total = _turning_deg(arrow.points_xyz)
        assert abs(net) < params.max_turn_deg + 15.0
        assert total < params.max_total_turn_deg + 15.0
        assert _arc_length_deg(arrow.points_xyz) < params.max_length_deg + 1.0
    assert max(abs(_turning_deg(a.points_xyz)[0]) for a in arrows) > 90.0, "the turn limit, not length, should end some arcs"


def test_tracing_is_deterministic():
    u, v = _vortex_about_x()
    valid = np.ones((H, W), dtype=bool)
    a = flow_arrows.trace_flow_arrows(u, v, valid)
    b = flow_arrows.trace_flow_arrows(u, v, valid)
    assert len(a) == len(b)
    for x, y in zip(a, b):
        np.testing.assert_array_equal(x.points_xyz, y.points_xyz)


def test_faster_flow_reports_a_higher_mean_speed():
    u, v = _eastward()
    u = np.where(LAT > 0, 3.0, 1.0)
    arrows = flow_arrows.trace_flow_arrows(u, v, np.ones((H, W), dtype=bool))
    north = [a.mean_speed for a in arrows if a.points_xyz[0, 2] > 0.2]
    south = [a.mean_speed for a in arrows if a.points_xyz[0, 2] < -0.2]
    assert north and south
    assert min(north) > max(south)


def test_smooth_spline_keeps_its_ends_and_rounds_corners():
    # A zigzag control polygon: the spline must start and end on the end points, stay on the
    # sphere, be sampled about every half degree, and round each 60-degree corner into a bend
    # with no sharp turn between neighbouring samples.
    lat = np.radians([0.0, 3.0, 0.0, 3.0, 0.0, 3.0])
    lon = np.radians([0.0, 3.0, 6.0, 9.0, 12.0, 15.0])
    control = geometry.latlon_to_xyz(lat, lon)
    curve = flow_arrows._smooth_spline(control)
    np.testing.assert_allclose(curve[[0, -1]], control[[0, -1]], atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(curve, axis=1), 1.0)
    assert np.degrees(geometry.angular_distance(curve[:-1], curve[1:])).max() < 0.6
    _, total = _turning_deg(curve)
    _, polygon_total = _turning_deg(control)
    assert total < 0.6 * polygon_total


def test_one_non_finite_cell_does_not_remove_the_arrows():
    u, v = _eastward()
    valid = np.ones((H, W), dtype=bool)
    expected = len(flow_arrows.trace_flow_arrows(u, v, valid))
    u = u.copy()
    u[H // 2, W // 2] = np.nan
    assert len(flow_arrows.trace_flow_arrows(u, v, valid)) >= expected - 2


@pytest.mark.parametrize("lat_deg", [80.0, 60.0])
def test_spacing_grid_marks_every_cell_within_the_radius_at_high_latitude(lat_deg):
    # The cap around a path point reaches furthest in longitude poleward of its centre, at
    # sin(lat) = sin(lat0) / cos(r), not on the point's own row, so the column range must come
    # from the cap and not from each row's latitude. At 80 N the 18-degree cap covers the pole
    # (the full-width branch); at 60 N it does not (the arcsin branch).
    radius = np.radians(18.0)
    path = geometry.latlon_to_xyz(np.radians([lat_deg, lat_deg]), np.radians([0.0, 2.0]))
    distance = np.full((180, 360), np.inf)
    flow_arrows._stamp_distance(distance, path, radius)
    cell_lat = np.radians(90.0 - (np.arange(180) + 0.5))
    cell_lon = np.radians(-180.0 + (np.arange(360) + 0.5))
    cells = geometry.latlon_to_xyz(*np.meshgrid(cell_lat, cell_lon, indexing="ij"))
    true = np.min(geometry.angular_distance(cells[:, :, None, :], path[None, None, :, :]), axis=-1)
    inside = true < radius
    assert np.all(np.isfinite(distance[inside]))
    np.testing.assert_allclose(distance[inside], true[inside], atol=1e-9)


# A 1-degree grid for the coast tests: the coastal band is 2 degrees, one cell on the coarse
# grid above.
H1, W1 = 180, 360
LON1 = np.broadcast_to(np.radians(-180.0 + (np.arange(W1) + 0.5)), (H1, W1))


def _eastward_into_a_continent():
    """Eastward flow everywhere, with a continent from 20 to 60 degrees east, pole to pole."""
    valid = ~((LON1 > np.radians(20.0)) & (LON1 < np.radians(60.0)))
    return np.ones((H1, W1)), np.zeros((H1, W1)), valid


def test_coast_weight_is_one_at_the_coast_and_zero_beyond_the_band():
    _, _, valid = _eastward_into_a_continent()
    weight = flow_arrows._coast_weight(valid, 2.0)
    row = H1 // 2 - 30  # 30 N
    lon_deg = np.degrees(LON1[row])
    # 1 at the shore itself; a coastal cell's centre is half a cell out.
    coastal, next_out = weight[row, np.argmin(np.abs(lon_deg - 19.5))], weight[row, np.argmin(np.abs(lon_deg - 18.5))]
    assert 0.75 < coastal <= 1.0
    assert 0.0 < next_out < coastal
    assert np.all(weight[row, (lon_deg < 17.0) & (lon_deg > -170.0)] == 0.0)
    assert np.all(weight[~valid] == 0.0)
    assert not flow_arrows._coast_weight(np.ones((H1, W1), dtype=bool), 2.0).any()


def test_coast_weight_reaches_cells_next_to_land_on_a_coarse_grid():
    # On a 2-degree grid (climate density 1.0) cell centres sit a degree from the shore. At
    # 29 N the second cell out is about 3.5 degrees from the nearest land centre, 2.5 from the
    # shore, so a 3-degree band must still reach it; measuring the band between centres cut
    # it off.
    lon = np.radians(-180.0 + (np.arange(180) + 0.5) * 2.0)
    valid = np.broadcast_to(~((lon > np.radians(20.0)) & (lon < np.radians(60.0))), (90, 180))
    weight = flow_arrows._coast_weight(valid, 3.0)
    lon_deg = np.degrees(lon)
    row = 30  # 29 N
    assert weight[row, np.argmin(np.abs(lon_deg - 19.0))] > 0.5
    assert weight[row, np.argmin(np.abs(lon_deg - 17.0))] > 0.05
    assert weight[row, np.argmin(np.abs(lon_deg - 13.0))] == 0.0


def test_currents_turn_along_a_coast_in_the_models_hemisphere_sense():
    # climate._deflect_into_coast sends flow that meets land along the coast: clockwise around
    # the land in the north, anticlockwise in the south. Flow striking a coast that faces west
    # therefore turns south in the northern hemisphere and north in the southern one. The
    # arrows must show that turn rather than stop at the coast.
    u, v, valid = _eastward_into_a_continent()
    arrows = flow_arrows.trace_flow_arrows(u, v, valid, flow_arrows.CURRENT_PARAMS)
    turned = {"north": 0.0, "south": 0.0}
    for arrow in arrows:
        lat, lon = geometry.xyz_to_latlon(arrow.points_xyz)
        near_coast = (np.degrees(lon) > 15.0) & (np.degrees(lon) < 20.0)
        dlat = np.degrees(np.diff(lat))[near_coast[1:] & near_coast[:-1]]
        hemisphere = "north" if np.mean(lat) > 0 else "south"
        expected = -1.0 if hemisphere == "north" else 1.0
        turned[hemisphere] += float(np.sum(dlat[np.sign(dlat) == expected]) * expected)
        assert not np.any(np.sum(dlat[np.sign(dlat) == -expected]) * -expected > 2.0), "ran along the coast the wrong way"
    assert turned["north"] > 5.0 and turned["south"] > 5.0


def test_arrows_turning_at_a_coast_stay_off_the_land():
    u, v, valid = _eastward_into_a_continent()
    for arrow in flow_arrows.trace_flow_arrows(u, v, valid, flow_arrows.CURRENT_PARAMS):
        _, lon = geometry.xyz_to_latlon(arrow.points_xyz)
        assert not np.any((np.degrees(lon) > 20.5) & (np.degrees(lon) < 59.5))



def test_a_non_finite_patch_in_an_all_valid_field_gets_no_coastal_turning():
    # Wind is passed as all valid. A patch of NaN is a hole arrows stop at, not a coast: the
    # ocean's coastal rule must not curve uniform eastward flow around it. Counting the NaN
    # cells as land turned the flow into impact seeds that bent arrows south and north.
    u, v, _ = _eastward_into_a_continent()
    u = u.copy()
    u[86:94, 20 + 180 : 28 + 180] = np.nan
    arrows = flow_arrows.trace_flow_arrows(u, v, np.ones((H1, W1), dtype=bool), flow_arrows.CURRENT_PARAMS)
    assert arrows
    for arrow in arrows:
        lat, _ = geometry.xyz_to_latlon(arrow.points_xyz)
        assert np.ptp(lat) < np.radians(1.0), "an arrow bent around a NaN patch"


def test_the_spline_falls_back_to_the_traced_path_when_every_spacing_crosses_land():
    # Land fills lon > 0 and 0 < lat < 10. A trace runs north along its west side and east
    # along its north side, 0.1 degrees off the shore, turning the corner at (0, 10). Any
    # smoothing rounds that corner by a quarter of the control spacing, into the land, even at
    # the tightest spacing: the arrow must come back as the traced path itself.
    lat_c = 90.0 - (np.arange(180) + 0.5)
    lon_c = -180.0 + (np.arange(360) + 0.5)
    lat_g, lon_g = np.meshgrid(lat_c, lon_c, indexing="ij")
    mask = np.where((lon_g > 0) & (lat_g > 0) & (lat_g < 10), 0.0, 1.0)
    west = np.stack([np.full(50, -0.1), np.linspace(0.5, 10.1, 50)], axis=-1)  # (lon, lat)
    north = np.stack([np.linspace(-0.1, 10.5, 53)[1:], np.full(52, 10.1)], axis=-1)
    lonlat = np.concatenate([west, north])
    path = geometry.latlon_to_xyz(np.radians(lonlat[:, 1]), np.radians(lonlat[:, 0]))
    coast = np.ones(len(path))
    params = flow_arrows.FlowArrowParams()

    assert flow_arrows._off_water_count(mask, path) == 0
    tightest = params.min_control_spacing_deg
    smoothest = flow_arrows._smooth_spline(path[flow_arrows._control_indices(path, coast, tightest, tightest)])
    assert flow_arrows._off_water_count(mask, smoothest) > 0, "the scenario must defeat the tightest smoothing"

    curve = flow_arrows._spline_off_land(path, coast, mask, params)
    assert flow_arrows._off_water_count(mask, curve) == 0
    on_path = np.min(geometry.angular_distance(curve[:, None, :], path[None, :, :]), axis=1)
    np.testing.assert_allclose(on_path, 0.0, atol=1e-6)  # arccos near 1 is only good to about 1e-8


def test_spline_off_land_keeps_the_smoothest_arc_that_stays_on_water():
    # Open water: the widest spacing already stays on water, so the arrow is the smooth
    # spline at that spacing, not the traced path.
    lat = np.radians(np.linspace(0.0, 4.0, 40))
    lon = np.radians(np.linspace(0.0, 30.0, 40) + 0.8 * np.sin(np.linspace(0.0, 6.0, 40)))
    path = geometry.latlon_to_xyz(lat, lon)
    params = flow_arrows.FlowArrowParams()
    curve = flow_arrows._spline_off_land(path, np.zeros(len(path)), np.ones((180, 360)), params)
    expected = flow_arrows._smooth_spline(path[flow_arrows._control_indices(path, np.zeros(len(path)), params.control_spacing_deg, params.coast_control_spacing_deg)])
    np.testing.assert_array_equal(curve, expected)


def test_spline_off_land_uses_the_tracers_bilinear_land_test():
    # `_is_valid` treats a point as water when the bilinear mask is at least 0.5. A point in a
    # one-cell hole but near its corner is water by that test (the hole's weight there is only
    # 0.36), though its containing cell is land; the check must agree with the tracer, not
    # with the cell, or it tightens arrows for nothing.
    mask = np.ones((180, 360))
    mask[90, 180] = 0.0  # lat -1..0, lon 0..1
    near_corner = geometry.latlon_to_xyz(np.radians([-0.9]), np.radians([0.1]))
    assert flow_arrows._off_water_count(mask, near_corner) == 0
    assert mask[flow_arrows._cells(mask.shape, near_corner)][0] == 0.0, "the point must be inside the hole cell"
    centre = geometry.latlon_to_xyz(np.radians([-0.5]), np.radians([0.5]))
    assert flow_arrows._off_water_count(mask, centre) == 1


def test_min_control_spacing_must_be_positive():
    with pytest.raises(ValueError):
        flow_arrows.FlowArrowParams(min_control_spacing_deg=0.0)

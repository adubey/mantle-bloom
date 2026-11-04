import numpy as np
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
    params = flow_arrows.FlowArrowParams(max_length_deg=400.0, max_turn_deg=120.0)
    arrows = flow_arrows.trace_flow_arrows(u, v, np.ones((H, W), dtype=bool), params)
    assert arrows
    for arrow in arrows:
        pts = arrow.points_xyz
        # Total turning of the chord directions in each point's tangent plane.
        d = np.diff(pts, axis=0)
        d -= np.sum(d * pts[:-1], axis=1, keepdims=True) * pts[:-1]
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        turning = np.degrees(np.sum(np.arccos(np.clip(np.sum(d[1:] * d[:-1], axis=1), -1, 1))))
        assert turning < params.max_turn_deg + 15.0
        assert _arc_length_deg(pts) < 300.0


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


def test_catmull_rom_passes_through_its_control_points():
    lat = np.radians([0.0, 5.0, 8.0, 6.0, 10.0])
    lon = np.radians([0.0, 4.0, 10.0, 16.0, 20.0])
    control = geometry.latlon_to_xyz(lat, lon)
    curve = flow_arrows._catmull_rom(control, samples_per_segment=6)
    np.testing.assert_allclose(curve[::6], control, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(curve, axis=1), 1.0)

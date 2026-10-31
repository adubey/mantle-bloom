import numpy as np
import pytest

from app import slab_body
from app.sparse_quad_patch import PlateWithSparseQuadPatch, pack_cell_keys
from app.lithosphere import PLANET_RADIUS_M


def _ribbon() -> slab_body.SlabBody:
    return slab_body.create_slab_body(
        body_id=4,
        plate_id=2,
        plate_frame=np.eye(3),
        trench_world_xyz=np.array([1.0, 0.0, 0.0]),
        down_dip_world_xyz=np.array([0.0, 1.0, 0.0]),
        source_cell_keys=np.array([10, 11, 12]),
        width_m=200_000.0,
        thickness_m=100_000.0,
        max_depth_m=660_000.0,
        station_count=8,
        dip_degrees=45.0,
        density_anomaly_kg_m3=80.0,
        surface_temperature_k=273.0,
        geothermal_gradient_k_per_m=3e-4,
        bending_modulus_pa=70e9,
    )


def test_ribbon_has_measured_3d_geometry_and_volume_scaled_attachment_force():
    body = _ribbon()

    metrics = body.metrics()
    force_on_plate, reaction_on_slab = slab_body.attachment_force_local_n(body)

    assert metrics["station_count"] == 8
    assert metrics["length_m"] > metrics["max_depth_m"]
    assert metrics["max_depth_m"] == pytest.approx(660_000.0)
    assert metrics["dip_degrees"] == pytest.approx(45.0, abs=0.5)
    assert metrics["width_m"] == 200_000.0
    assert metrics["thickness_m"] == 100_000.0
    assert body.consumed_volume_m3 == pytest.approx(body.width_m * body.thickness_m * metrics["length_m"])
    np.testing.assert_allclose(force_on_plate + reaction_on_slab, np.zeros(3), atol=0.0)
    assert force_on_plate[1] > 0.0
    assert np.linalg.norm(force_on_plate) == pytest.approx(
        np.dot(body.station_volume_m3, body.density_anomaly_kg_m3) * slab_body.GRAVITY_M_S2
    )


def test_implicit_bending_and_drag_advance_is_finite_and_accounts_material():
    body = _ribbon()
    volume_before = body.consumed_volume_m3
    path_before = body.centerline_local_xyz_m.copy()

    budget = slab_body.advance_slab_body(
        body,
        2_000_000.0,
        mantle_viscosity_pa_s=1e21,
        shear_zone_thickness_m=20_000.0,
    )

    assert not np.allclose(body.centerline_local_xyz_m, path_before)
    assert np.all(np.isfinite(body.centerline_local_xyz_m))
    assert body.consumed_volume_m3 == pytest.approx(volume_before)
    assert 0.0 <= body.deep_reservoir_volume_m3 <= volume_before
    assert budget["work_residual_fraction"] < 0.05
    assert budget["max_station_displacement_m"] > 0.0
    assert body.metrics()["max_depth_m"] < PLANET_RADIUS_M


def test_attached_ribbons_drive_shell_strain_and_endogenous_patch_geometry():
    n = 12
    jj, ii = np.meshgrid(np.arange(3, 9), np.arange(3, 9), indexing="ij")
    keys = pack_cell_keys(np.zeros(ii.size, dtype=int), ii.ravel(), jj.ravel())
    plate = PlateWithSparseQuadPatch(
        2,
        np.eye(3),
        "oceanic",
        n,
        keys,
        fields={
            "crustal_thickness_m": np.full(len(keys), 12_000.0),
            "mantle_lithosphere_thickness_m": np.full(len(keys), 80_000.0),
        },
    )
    body = slab_body.create_slab_body(
        body_id=9,
        plate_id=plate.plate_id,
        plate_frame=plate.frame,
        trench_world_xyz=np.array([1.0, 0.0, 0.0]),
        down_dip_world_xyz=np.array([0.0, 1.0, 0.0]),
        source_cell_keys=np.array([10, 11]),
        width_m=300_000.0,
        thickness_m=100_000.0,
        max_depth_m=300_000.0,
        station_count=6,
        dip_degrees=45.0,
        density_anomaly_kg_m3=80.0,
        surface_temperature_k=273.0,
        geothermal_gradient_k_per_m=3e-4,
        bending_modulus_pa=70e9,
    )
    opposite_body = slab_body.create_slab_body(
        body_id=10,
        plate_id=plate.plate_id,
        plate_frame=plate.frame,
        trench_world_xyz=np.array([np.cos(0.25), -np.sin(0.25), 0.0]),
        down_dip_world_xyz=np.array([-np.sin(0.25), -np.cos(0.25), 0.0]),
        source_cell_keys=np.array([12, 13]),
        width_m=300_000.0,
        thickness_m=100_000.0,
        max_depth_m=300_000.0,
        station_count=6,
        dip_degrees=45.0,
        density_anomaly_kg_m3=80.0,
        surface_temperature_k=273.0,
        geothermal_gradient_k_per_m=3e-4,
        bending_modulus_pa=70e9,
    )
    before_geometry = plate.geometry_revision

    response = slab_body.apply_slab_load_to_plate_geometry(
        plate,
        [body, opposite_body],
        2_000_000.0,
        effective_extensional_viscosity_pa_s=1e24,
        load_spread_degrees=15.0,
        extensional_capacity_resultant_n_per_m=1e13,
    )

    first_force = slab_body.attachment_force_local_n(body)[0]
    second_force = slab_body.attachment_force_local_n(opposite_body)[0]
    assert first_force[1] * second_force[1] < 0.0
    assert np.linalg.norm(response["force_transfer_error_n"]) < 1e-12 * np.linalg.norm(first_force)
    assert np.linalg.norm(response["strain_increment"]) > 0.0
    assert response["equilibrium_residual_n"] < 1e-2 * np.linalg.norm(response["applied_force_n"])
    assert plate.geometry_revision > before_geometry
    assert np.all(np.isfinite(response["failure_ratio"]))

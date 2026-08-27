import numpy as np
import pytest

from app.world import generate_world, step_world

# Coarse settings throughout -- this is a smoke test (nothing crashes, state stays finite and
# in a sane range across several real steps), not a physics-precision test; see
# unit_tests/test_lithosphere.py/test_torque.py/test_healpix_grid.py for those.
_COARSE_KWARGS = dict(node_density=0.5, climate_density=0.5, fluid_density=0.5, num_plates=6)


@pytest.fixture(scope="module")
def stepped_world():
    world = generate_world(seed=11, **_COARSE_KWARGS)
    for _ in range(3):
        step_world(world, 3_000_000.0)
    return world


def test_generate_world_produces_plates_with_finite_elevation():
    world = generate_world(seed=11, **_COARSE_KWARGS)
    assert len(world.plates) == 6
    all_elevation = np.concatenate([p.collect("elevation") for p in world.plates])
    assert len(all_elevation) > 0
    assert np.all(np.isfinite(all_elevation))
    assert world.atmosphere_cfd_state is not None
    assert world.ocean_cfd_state is not None


def test_step_world_advances_time_and_stays_finite(stepped_world):
    world = stepped_world
    assert world.elapsed_years == pytest.approx(9_000_000.0)
    all_elevation = np.concatenate([p.collect("elevation") for p in world.plates])
    assert np.all(np.isfinite(all_elevation))
    assert all_elevation.min() >= -11000.0
    assert all_elevation.max() <= 9000.0


def test_step_world_keeps_fluid_state_bounded(stepped_world):
    world = stepped_world
    assert np.all(np.isfinite(world.atmosphere_cfd_state.u))
    assert np.all(np.isfinite(world.atmosphere_cfd_state.v))
    assert np.all(np.isfinite(world.ocean_cfd_state.u))
    assert np.all(np.isfinite(world.ocean_cfd_state.eta))
    # A real stability check, not just finiteness -- current/wind speeds should stay in a
    # physically plausible range, not have blown up (see fluid_dynamics_healpix.py's own
    # CFL_STENCIL_SAFETY_DIVISOR note for the instability this guards against regressing).
    assert np.abs(world.ocean_cfd_state.u).max() < 50.0
    assert np.abs(world.atmosphere_cfd_state.u).max() < 300.0


def test_step_world_populates_hydrology_and_climate_caches(stepped_world):
    world = stepped_world
    assert world.climate_cache is not None
    assert world.hydrology_cache is not None


def test_isostasy_derived_elevation_matches_hc_hm_state_at_generation():
    """Immediately after generation (no deform/erosion has run yet), every plate's
    `elevation` is exactly `isostatic_elevation(Hc, Hm)` -- the derived-field contract
    lithosphere_plate.py/lithosphere.py establish at construction time."""
    from app import lithosphere

    world = generate_world(seed=11, **_COARSE_KWARGS)
    for plate in world.plates:
        rho_c = lithosphere.crust_density(plate.crust_type)
        for line in plate.lines:
            if len(line) == 0:
                continue
            expected = lithosphere.isostatic_elevation(line.crustal_thickness_m, line.mantle_lithosphere_thickness_m, rho_c)
            assert np.allclose(line.elevation, expected, atol=1e-6)


def test_deform_applies_isostasy_as_a_delta_not_an_overwrite(stepped_world):
    """After stepping (which runs erosion every step, mutating `elevation` directly with no
    notion of Hc/Hm -- see lithosphere_plate.LithospherePlate.deform's own comment on why this
    must be a delta), `elevation` should generally have drifted *away* from a bare
    `isostatic_elevation(Hc, Hm)` readout on at least some nodes -- confirming erosion's
    contribution survived the subsequent deform() call rather than being silently erased by
    an unconditional overwrite back to pure isostasy."""
    from app import lithosphere

    world = stepped_world
    any_diverged = False
    for plate in world.plates:
        rho_c = lithosphere.crust_density(plate.crust_type)
        for line in plate.lines:
            if len(line) == 0:
                continue
            pure_isostasy = lithosphere.isostatic_elevation(line.crustal_thickness_m, line.mantle_lithosphere_thickness_m, rho_c)
            if not np.allclose(line.elevation, pure_isostasy, atol=1e-6):
                any_diverged = True
    assert any_diverged

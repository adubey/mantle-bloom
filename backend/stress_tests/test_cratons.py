"""Stepped-world craton checks too slow for unit_tests (issue #340)."""

import numpy as np

from app import cratons
from app import world as world_mod


def test_a_column_left_below_its_craton_is_booked_to_audit_clipped_not_a_mechanism(monkeypatch):
    world = world_mod.generate_world(5, num_plates=4, node_density=0.1)
    world.simulate_climate_biomes = False
    plate = max(world.plates, key=lambda p: cratons.plate_volume_m3(world, p))
    cell = int(np.argmax(plate.collect("craton_crust_m")))
    clamp = world_mod._clamp_column_caps

    def unbooked_thinning(w):
        # A writer that thins a cratonic column without booking it.
        hc = plate.collect("crustal_thickness_m")
        hc[cell] = 0.5 * plate.collect("craton_crust_m")[cell]
        plate.set_fields_on_plate(crustal_thickness_m=hc)
        clamp(w)

    monkeypatch.setattr(world_mod, "_clamp_column_caps", unbooked_thinning)
    before = dict(world.craton_ledger)

    world_mod.step_world(world, 100_000)

    assert world.craton_ledger["audit_clipped_m3"] > 0.0
    assert world.craton_ledger["delaminated_m3"] == before["delaminated_m3"]
    assert abs(cratons.balance_error_m3(world)) <= 1e-9 * world.craton_ledger["initial_m3"]

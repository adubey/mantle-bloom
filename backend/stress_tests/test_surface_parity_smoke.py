"""The issue #247/#249 plate-surface audit harness's CI-sized smoke configuration: one seed,
4 Myr at density 0.5. Hard invariants must hold, and the metrics document must be reproducible
byte for byte. Coverage and mesh gates are evaluated but not asserted here -- they track the
engine's current quality, which the long campaign (#249) judges."""

from app import surface_parity as sp
from app import surface_parity_gates as gates

HARD = ("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "S1")


def test_smoke_run_holds_every_hard_invariant_and_is_reproducible(tmp_path):
    config = sp.PRESETS["smoke"]
    run = sp.run_world(config, 3, tmp_path, log=lambda _: None)
    evaluation = gates.evaluate([run])

    hard = {key: gate for key, gate in evaluation["gates"].items() if gate["gate"] in HARD}
    assert {gate["gate"] for gate in hard.values()} == set(HARD)
    assert [key for key, gate in hard.items() if gate["status"] == gates.FAIL] == []
    document, _ = run
    assert document["audits"] == config.steps + 1
    assert len(document["checkpoints"]) == len(config.checkpoint_steps)
    assert all(check["state_identical"] for check in document["load_checks"])

    first = (tmp_path / "seed3.json").read_text()
    sp.run_world(config, 3, tmp_path, log=lambda _: None)
    assert (tmp_path / "seed3.json").read_text() == first

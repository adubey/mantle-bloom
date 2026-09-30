"""The issue #247 parity harness's CI-sized smoke configuration: one seed, both surfaces, 4 Myr at
density 0.5 (about 20 s). Hard invariants must hold on both surfaces, and the metrics document
must be reproducible byte for byte. Coverage/conservation/parity gates are evaluated but not
asserted here -- they track the engine's current quality, which the long campaign (#249)
judges."""

from app import surface_parity as sp
from app import surface_parity_gates as gates

HARD = ("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "S1")


def test_smoke_parity_run_holds_every_hard_invariant_and_is_reproducible(tmp_path):
    config = sp.PRESETS["smoke"]
    runs = [sp.run_surface(config, 3, surface, tmp_path, log=lambda _: None) for surface in sp.SURFACES]
    evaluation = gates.evaluate(runs)

    hard = {key: gate for key, gate in evaluation["gates"].items() if gate["gate"] in HARD}
    assert {gate["gate"] for gate in hard.values()} == set(HARD)
    # Line-surface findings are `info`; any quad or neutral failure is a regression.
    assert [key for key, gate in hard.items() if gate["status"] == gates.FAIL] == []
    for run, _ in runs:
        assert run["audits"] == config.steps + 1
        assert len(run["checkpoints"]) == len(config.checkpoint_steps)
        assert all(check["state_identical"] for check in run["load_checks"])

    first = (tmp_path / "seed3-quad.json").read_text()
    sp.run_surface(config, 3, "quad", tmp_path, log=lambda _: None)
    assert (tmp_path / "seed3-quad.json").read_text() == first

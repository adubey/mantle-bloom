"""Unit-suite default: the glacier march is opt-in (see `hydrology.MARCH_GLACIERS`).

The march is production behavior (`World.simulate_glacier_march`, default on) and costs
minutes per long step on generated worlds (#806, #822): the unit suite cannot fit CI
with it on everywhere, so it is off by default here and tests that assert marched ice
request the `march_glaciers` fixture below. Everything else runs the same zero-length
path as the world flag off (one opening-state solve, ice frozen for the step).
"""

import pytest

from app.hydroclimate import hydrology


@pytest.fixture(autouse=True, scope="session")
def _unit_tests_do_not_march_glaciers():
    # Session scope: module-scoped world fixtures (test_world_smoke's stepped_world) are set up
    # before any function-scoped autouse fixture, so a function-scoped patch left them marching.
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(hydrology, "MARCH_GLACIERS", False)
        yield


@pytest.fixture
def march_glaciers(monkeypatch):
    """Opt back into the glacier march, for the few unit tests that assert marched ice
    (short, bounded marches only -- full-length marches belong in stress_tests)."""
    monkeypatch.setattr(hydrology, "MARCH_GLACIERS", True)

import numpy as np
from app import geometry, volcanism
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def _line_plate(plate_id: int, theta: np.ndarray) -> Plate:
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(len(theta), 200.0))
    return Plate(plate_id=plate_id, frame=np.eye(3), crust_type="continental", lines=[line])


def _two_plate_world(gap: float, d: float = 0.01, n: int = 10) -> World:
    """Two plates, each a single dense line (spacing d) at phi=0 -- plate A's line ends at
    theta=(n-1)*d, plate B's line starts at that same value plus `gap`. Both share the
    identity frame purely so their local theta coordinates line up directly with world-space
    angular distance, for an exact, hand-verifiable test geometry. An *open* chain (plate A's
    own start and plate B's own end are dead ends, near nothing) -- fine for a "does this
    fire" check, but see _two_plate_ring_world for why a real "should NOT fire" control needs
    a closed ring instead."""
    theta_a = d * np.arange(n)
    plate_a = _line_plate(0, theta_a)
    theta_b = theta_a[-1] + gap + d * np.arange(n)
    plate_b = _line_plate(1, theta_b)
    return World(seed=0, plates=[plate_a, plate_b], next_plate_id=2)


def _two_plate_ring_world(gap_multiplier: float, n: int = 10) -> World:
    """Two plates whose lines close into a full ring around the sphere (both junctions --
    A's end to B's start, and B's end wrapping back to A's start -- the same size,
    `gap_multiplier * d`), so neither plate has an open "dead end" boundary point isolated
    from everything. Real generated plates always fully tile the sphere this way (every
    boundary point has a genuine neighbor on some side); a plain open two-plate chain
    (_two_plate_world) doesn't, and its own dead ends read as spuriously "infinitely far from
    any other plate" -- which is real outlier detection working correctly on an unrealistic
    shape, not a bug, but it means a "does this normal boundary NOT fire" control needs this
    closed-ring construction instead to mean what it says."""
    d = np.pi / ((n - 1) + gap_multiplier)
    gap = gap_multiplier * d
    theta_a = d * np.arange(n)
    plate_a = _line_plate(0, theta_a)
    theta_b = theta_a[-1] + gap + d * np.arange(n)
    plate_b = _line_plate(1, theta_b)
    return World(seed=0, plates=[plate_a, plate_b], next_plate_id=2)


def test_volcanism_integration_is_deterministic_for_the_same_seed():
    def run():
        # node_density=0.5 (the coarsest choice, an eighth of the default 4.0, see
        # plates.NODE_DENSITY_CHOICES) -- this test only checks determinism-with-itself, not a
        # specific outcome, so it doesn't need the default's full node density.
        world = generate_world(seed=17, num_plates=10, continental_fraction=0.4, node_density=0.5)
        # volcanism.py only reacts to plate/boundary geometry, never to climate/erosion/
        # hydrology output (see World.simulate_climate_biomes) -- skipping that per-step
        # computation speeds this up without changing what it exercises.
        world.simulate_climate_biomes = False
        for _ in range(10):
            step_world(world, years=2_000_000)
        return sorted(world.volcanic_field_plate_ids), len(world.plates)

    assert run() == run()

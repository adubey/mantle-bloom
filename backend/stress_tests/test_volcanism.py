from app.world import generate_world, step_world


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
        return (
            len(world.plates),
            sorted((p.plate_id, int(p.collect("is_volcano").sum())) for p in world.plates),
        )

    assert run() == run()

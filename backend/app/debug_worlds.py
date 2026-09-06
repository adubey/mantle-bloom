"""Small, hand-scripted plate configurations for the "Debugging Worlds" Generate World tab
(see main.py's `POST /world/generate_debug`) -- fast-iteration fixtures for the gap-filling
problem (`lithosphere_plate._fill_corner_notch`, `gaps.py`) that don't depend on any real
save's history.

Each scenario places a handful of explicit seed points directly (no Voronoi tiling -- see
`_build_plates`'s nearest-seed-wins ownership predicate), assigns every plate an explicit
`World.pinned_omegas` entry so it moves exactly as scripted every step regardless of what the
real torque balance would have done (see `LithospherePlate.shift`), and runs at the coarsest
real `node_density` choice so a step is fast enough to iterate on quickly. `debug_diagnostics`
is on by default so the corner-notch decision log (docs/debugging.md) is populated from the
very first step; both are ordinary Controls the user can change afterward.

Every plate's own omega is derived from a set of pairwise "this boundary is
convergent/divergent" relationships (see `_omegas_from_relationships`) rather than hand-picked
per plate, so a plate at a junction of several boundaries (the triple-junction scenario) gets
a single rotation that's the best simultaneous fit to all of them, the same way a real
mantle-driven plate's motion is a net result of every boundary force acting on it at once --
just picked directly here instead of emerging from torque.py's own force balance."""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, mantle
from .elevation_lines import line_spacing_rad
from .lithosphere_plate import new_plate
from .world import World, finish_generation

# Coarsest real node_density choice (elevation_lines.NODE_DENSITY_CHOICES) -- these scenarios
# are meant for fast iteration, not visual fidelity, so the fewer nodes per plate the better.
DEBUG_WORLD_NODE_DENSITY = 0.5

# A fixed fraction of MAX_PLATE_RATE, not the rail itself -- keeps every scenario's motion
# representative of an ordinary plate (so `_fill_corner_notch`'s own window/budget math, which
# scales with `mantle.MAX_PLATE_RATE * years`, sees realistic per-step gaps to close) without
# the corner-notch fallback's reach constant alone swallowing whatever gap opens each step.
_DEBUG_WORLD_RATE = 0.3 * mantle.MAX_PLATE_RATE


def _tangent_toward(p_i: np.ndarray, p_j: np.ndarray) -> np.ndarray:
    """Unit tangent vector at `p_i` pointing along the great circle toward `p_j`."""
    return geometry.normalize(p_j - p_i * np.dot(p_i, p_j))


def _omegas_from_relationships(
    seeds_xyz: list[np.ndarray], relationships: list[tuple[int, int, str]], rate: float
) -> list[np.ndarray]:
    """One omega per seed in `seeds_xyz`, each the rotation that best (in a simple summed-
    tangent sense) satisfies every relationship it appears in. `relationships` is a list of
    `(i, j, "divergent" | "convergent")` -- symmetric in effect (both plates move apart, or
    both move together), so the same sign contributes to both `i`'s and `j`'s own tangent sum,
    each using its own direction toward the other. A plate combines its every relationship's
    tangent contribution by summing (so a triple-junction plate's own single rotation reflects
    every one of its boundaries at once, not just its first) before converting the net
    direction to a rotation axis: for unit `p` and tangent `v` (`v` orthogonal to `p`),
    `omega = cross(p, v)` is exactly the rotation with `cross(omega, p) = v` at rate `|v|`, so
    scaling `v` to unit length before this cross product gives a rotation of exactly `rate`."""
    sign = {"divergent": -1.0, "convergent": 1.0}
    tangent_sums = [np.zeros(3) for _ in seeds_xyz]
    for i, j, kind in relationships:
        s = sign[kind]
        tangent_sums[i] += s * _tangent_toward(seeds_xyz[i], seeds_xyz[j])
        tangent_sums[j] += s * _tangent_toward(seeds_xyz[j], seeds_xyz[i])

    omegas = []
    for p, v in zip(seeds_xyz, tangent_sums):
        if np.linalg.norm(v) < 1e-12:
            omegas.append(np.zeros(3))
            continue
        omegas.append(np.cross(p, geometry.normalize(v)) * rate)
    return omegas


def _build_plates(
    seed: int, seeds_latlon_deg: list[tuple[float, float]], crust_types: list[str], node_density: float
) -> tuple[list, list[np.ndarray]]:
    """One `LithospherePlate` per seed, partitioning the *entire* sphere between them by
    nearest-seed-wins (a manual, exact-position stand-in for `generate_plates`' own random
    Voronoi tiling -- see `new_plate`'s `is_owned` parameter). Returns the plates alongside
    their world-xyz seed positions, which callers need again to compute pinned omegas."""
    spacing_rad = line_spacing_rad(node_density)
    seeds_xyz = [
        geometry.latlon_to_xyz(np.radians(lat_deg), np.radians(lon_deg)) for lat_deg, lon_deg in seeds_latlon_deg
    ]
    tree = cKDTree(np.array(seeds_xyz))

    plates = []
    for i, (xyz, crust_type) in enumerate(zip(seeds_xyz, crust_types)):
        frame = geometry.plate_frame_from_seed(xyz)

        def is_owned(world_pts: np.ndarray, _i: int = i) -> np.ndarray:
            _, idx = tree.query(world_pts)
            return idx == _i

        plates.append(new_plate(i, frame, crust_type, spacing_rad, seed, is_owned=is_owned))
    return plates, seeds_xyz


def _build_debug_world(
    seed: int,
    seeds_latlon_deg: list[tuple[float, float]],
    crust_types: list[str],
    relationships: list[tuple[int, int, str]],
    scenario_message: str,
) -> World:
    plates, seeds_xyz = _build_plates(seed, seeds_latlon_deg, crust_types, DEBUG_WORLD_NODE_DENSITY)
    omegas = _omegas_from_relationships(seeds_xyz, relationships, _DEBUG_WORLD_RATE)

    world = World(
        seed=seed,
        plates=plates,
        mantle_centers=[],  # irrelevant -- every plate's motion is pinned, never torque-driven
        next_plate_id=len(plates),
        node_density=DEBUG_WORLD_NODE_DENSITY,
        pinned_omegas={plate.plate_id: omega for plate, omega in zip(plates, omegas)},
        debug_diagnostics=True,
    )
    finish_generation(world, scenario_message)
    return world


def scenario_two_plate_divergent(seed: int) -> World:
    """Two oceanic plates straddling a shared meridian, spinning directly apart -- the
    simplest possible rift-boundary gap-filling fixture."""
    return _build_debug_world(
        seed,
        seeds_latlon_deg=[(0.0, -10.0), (0.0, 10.0)],
        crust_types=["oceanic", "oceanic"],
        relationships=[(0, 1, "divergent")],
        scenario_message="Debugging world generated: two-plate divergent boundary.",
    )


def scenario_two_plate_convergent(seed: int) -> World:
    """Same layout as `scenario_two_plate_divergent`, spinning together instead -- a closing
    boundary (subduction/collision test), one plate continental so both an oceanic-subduction
    and a continental-suture margin can be exercised by switching which plate the boundary
    favours."""
    return _build_debug_world(
        seed,
        seeds_latlon_deg=[(0.0, -10.0), (0.0, 10.0)],
        crust_types=["oceanic", "continental"],
        relationships=[(0, 1, "convergent")],
        scenario_message="Debugging world generated: two-plate convergent boundary.",
    )


def scenario_triple_junction_mixed(seed: int) -> World:
    """Three plates meeting at one point, two legs divergent and one convergent -- the exact
    shape `_fill_corner_notch`'s own docstring calls out as the case `_stretch_end` (theta-
    axis-only) and `_claim_adjacent_territory` (phi-axis-only, whole rows) structurally can't
    reach on their own (confirmed on a real save, seed 430031492)."""
    return _build_debug_world(
        seed,
        seeds_latlon_deg=[(10.0, 0.0), (-10.0, 10.0), (-10.0, -10.0)],
        crust_types=["continental", "oceanic", "oceanic"],
        relationships=[(0, 1, "divergent"), (0, 2, "divergent"), (1, 2, "convergent")],
        scenario_message="Debugging world generated: triple junction (mixed divergent/convergent).",
    )


def scenario_four_plate_grid(seed: int) -> World:
    """Four oceanic plates in a 2x2 arrangement, every edge of the square divergent -- four
    simultaneous triple-junction-like corners opening at once, around one shared center point."""
    return _build_debug_world(
        seed,
        seeds_latlon_deg=[(10.0, -10.0), (10.0, 10.0), (-10.0, 10.0), (-10.0, -10.0)],  # NW, NE, SE, SW
        crust_types=["oceanic", "oceanic", "oceanic", "oceanic"],
        relationships=[(0, 1, "divergent"), (1, 2, "divergent"), (2, 3, "divergent"), (3, 0, "divergent")],
        scenario_message="Debugging world generated: four-plate grid, all edges divergent.",
    )


def scenario_five_plate_irregular(seed: int) -> World:
    """Five plates at irregular (non-grid) seed positions with mixed divergent/convergent
    boundaries in a closed ring -- closest single scenario to a real save's messiness while
    staying small enough to iterate on quickly."""
    return _build_debug_world(
        seed,
        seeds_latlon_deg=[(15.0, 0.0), (-5.0, 20.0), (-15.0, -10.0), (5.0, -25.0), (0.0, 5.0)],
        crust_types=["continental", "continental", "oceanic", "oceanic", "oceanic"],
        relationships=[
            (0, 1, "divergent"),
            (1, 2, "convergent"),
            (2, 3, "divergent"),
            (3, 4, "convergent"),
            (4, 0, "divergent"),
        ],
        scenario_message="Debugging world generated: five-plate irregular ring.",
    )


DEBUG_SCENARIOS = {
    "two_plate_divergent": scenario_two_plate_divergent,
    "two_plate_convergent": scenario_two_plate_convergent,
    "triple_junction_mixed": scenario_triple_junction_mixed,
    "four_plate_grid": scenario_four_plate_grid,
    "five_plate_irregular": scenario_five_plate_irregular,
}

# Short human-readable labels for the Generate World "Debugging Worlds" tab's scenario picker
# (see main.py's GET /world/debug_scenarios) -- pulled from each scenario function's own
# one-line intent rather than duplicating the docstring's full explanation.
DEBUG_SCENARIO_LABELS = {
    "two_plate_divergent": "Two plates, divergent boundary",
    "two_plate_convergent": "Two plates, convergent boundary",
    "triple_junction_mixed": "Triple junction (2 divergent, 1 convergent)",
    "four_plate_grid": "Four-plate grid, all edges divergent",
    "five_plate_irregular": "Five plates, irregular ring",
}


def generate_debug_world(scenario: str, seed: int = 0) -> World:
    if scenario not in DEBUG_SCENARIOS:
        raise ValueError(f"unknown debug scenario {scenario!r}; choices are {sorted(DEBUG_SCENARIOS)}")
    return DEBUG_SCENARIOS[scenario](seed)

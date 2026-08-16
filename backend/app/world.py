"""Top-level world state: a collection of plates plus elapsed simulation time."""

from __future__ import annotations

from dataclasses import dataclass, field

from .plates import Plate, generate_plates

DEFAULT_NUM_PLATES = 12


@dataclass
class World:
    seed: int
    plates: list[Plate] = field(default_factory=list)
    elapsed_years: float = 0.0


def generate_world(seed: int, num_plates: int = DEFAULT_NUM_PLATES) -> World:
    plates = generate_plates(seed, num_plates=num_plates)
    return World(seed=seed, plates=plates, elapsed_years=0.0)

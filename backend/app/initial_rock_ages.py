"""Initial per-plate rock ages, assigned without changing world-generation randomness."""

from __future__ import annotations

import numpy as np

from . import elevation_lines
from .plates import Plate

OCEANIC_MEAN_AGE_YEARS = 64.2e6
OCEANIC_MAX_AGE_YEARS = 3.0 * OCEANIC_MEAN_AGE_YEARS

# Geological boundaries are ages before present. These uniform within-era ranges are a
# deliberately simple proxy; the era shares are the regional Laurentia split in #744.
CONTINENTAL_AGE_RANGES_YEARS = {
    "Archean": (2.5e9, 4.0e9),
    "Proterozoic": (541e6, 2.5e9),
    "Phanerozoic": (0.0, 541e6),
}
CONTINENTAL_ERA_NAMES = ("Archean", "Proterozoic", "Phanerozoic")
CONTINENTAL_ERA_PROBABILITIES = (0.33, 0.43, 0.24)
_AGE_STREAM_TAG = 744


def assign_initial_rock_ages(plates: list[Plate], seed: int) -> dict[int, dict[str, float | str]]:
    """Stamp each starting cell with its plate's age for that cell's crust type.

    Age draws use a private seed stream so adding this initialization does not change any
    existing plate, mantle, or climate random choices. Stored ages are negative simulation
    years: ``-1`` remains the distinct pre-tracking sentinel for legacy worlds.
    """
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), _AGE_STREAM_TAG]))
    initial_ages: dict[int, dict[str, float | str]] = {}
    for plate in plates:
        oceanic_age_years, continental_age_years = sample_plate_ages(rng)

        codes = plate.collect("crust_type_code")
        is_continental = elevation_lines.effective_is_continental_from_codes(
            codes, plate.crust_type == "continental"
        )
        cell_birth_years = np.where(
            is_continental, -continental_age_years, -oceanic_age_years
        )
        plate.set_fields_on_plate(node_created_years=cell_birth_years)
        initial_ages[plate.plate_id] = {
            "primary_crust_type": plate.crust_type,
            "oceanic_age_years": oceanic_age_years,
            "continental_age_years": continental_age_years,
        }
    return initial_ages


def sample_plate_ages(rng: np.random.Generator) -> tuple[float, float]:
    """Draw one oceanic and continental age in positive years before t=0."""
    oceanic_age_years = float(rng.triangular(0.0, 0.0, OCEANIC_MAX_AGE_YEARS))
    era = str(rng.choice(CONTINENTAL_ERA_NAMES, p=CONTINENTAL_ERA_PROBABILITIES))
    era_min, era_max = CONTINENTAL_AGE_RANGES_YEARS[era]
    continental_age_years = float(rng.uniform(era_min, era_max))
    return oceanic_age_years, continental_age_years

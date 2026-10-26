"""Default world names in exoplanet style: a star name plus a lowercase discovery letter (#597)."""

import random
from collections.abc import Collection

# Real star names: known exoplanet hosts first, then bright proper-name stars.
STAR_NAMES: tuple[str, ...] = (
    "Kepler-22",
    "Kepler-62",
    "Kepler-69",
    "Kepler-186",
    "Kepler-442",
    "TRAPPIST-1",
    "Proxima Centauri",
    "Tau Ceti",
    "Epsilon Eridani",
    "51 Pegasi",
    "HD 209458",
    "Gliese 581",
    "Gliese 667 C",
    "Gliese 876",
    "LHS 1140",
    "K2-18",
    "TOI-700",
    "Ross 128",
    "Wolf 1061",
    "55 Cancri",
    "Upsilon Andromedae",
    "Mu Arae",
    "HD 40307",
    "47 Ursae Majoris",
    "Pollux",
    "Aldebaran",
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Achernar",
    "Betelgeuse",
    "Altair",
    "Antares",
    "Spica",
    "Fomalhaut",
    "Deneb",
    "Regulus",
    "Castor",
    "Bellatrix",
    "Alnilam",
    "Polaris",
    "Algol",
    "Alpheratz",
    "Diphda",
    "Mirfak",
    "Mizar",
)

# Discovery letters start at "b" ("a" is the star itself); c/d/e come up more
# often so a fresh world sounds like a mid-system planet (#597).
LETTER_WEIGHTS: tuple[tuple[str, int], ...] = (
    ("b", 1),
    ("c", 3),
    ("d", 3),
    ("e", 3),
    ("f", 1),
    ("g", 1),
)


def default_world_name(used: Collection[str] = (), rng: random.Random | None = None) -> str:
    """A name like "Kepler-442 d", avoiding every name in `used` (#597)."""
    draw = rng if rng is not None else random.Random()
    letters = [letter for letter, _ in LETTER_WEIGHTS]
    weights = [weight for _, weight in LETTER_WEIGHTS]
    for _ in range(100):
        name = f"{draw.choice(STAR_NAMES)} {draw.choices(letters, weights=weights)[0]}"
        if name not in used:
            return name
    # Unreachable in practice (hundreds of combinations): suffix a numeral
    # rather than looping forever once every combination is taken.
    n = 2
    while True:
        name = f"{draw.choice(STAR_NAMES)} {draw.choices(letters, weights=weights)[0]} {n}"
        if name not in used:
            return name
        n += 1

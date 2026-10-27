"""Default world names in exoplanet style: a star name plus a lowercase discovery letter (#597)."""

import random
from collections.abc import Collection

# Sun-like stars with no known (confirmed) planets, so a fresh world never borrows
# the name of a real planetary system (#640). Every entry below was checked against
# the NASA Exoplanet Archive's confirmed-planets table (`pscomppars`, distinct
# hostnames, October 2026) and came back absent; solar-twin/analog status follows
# the usual literature lists. That check also ruled out lookalikes that do host
# planets despite their status -- 16 Cygni B, HD 42618, HD 141004 and HD 45184
# among them -- plus debated cases (59 Virginis) and unconfirmed candidates
# (Alpha Centauri A), which are left out. The claim can go stale with new
# discoveries; re-run the archive check before extending this list.
STAR_NAMES: tuple[str, ...] = (
    "18 Scorpii",
    "HD 98618",
    "HIP 56948",
    "HIP 79672",
    "HIP 100963",
    "HIP 73815",
    "HD 143436",
    "HD 44594",
    "HD 162826",
    "HD 10307",
    "HD 34411",
    "HD 138573",
    "HD 197027",
    "HD 76151",
    "Sigma Draconis",
    "HD 142373",
    "HD 110897",
    "HD 207129",
    "HD 30495",
    "HD 158633",
    "61 Ursae Majoris",
    "36 Ursae Majoris",
    "HD 72905",
    "Kappa Ceti",
    "HD 39587",
    "HD 56124",
    "HD 76440",
    "HD 9986",
    "HD 73350",
    "HD 187237",
    "HD 191785",
    "HD 19632",
    "37 Geminorum",
    "Beta Canum Venaticorum",
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
    """A name like "18 Scorpii d", avoiding every name in `used` (#597)."""
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

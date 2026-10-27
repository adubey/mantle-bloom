"""world_names.py -- exoplanet-style default world names (#597)."""

import random
import re

from app.world_names import LETTER_WEIGHTS, STAR_NAMES, default_world_name

NAME_RE = re.compile(r".+ [b-g]$")


def test_generated_name_is_a_star_plus_a_letter():
    name = default_world_name(rng=random.Random(7))
    assert NAME_RE.fullmatch(name), name
    star, _, _letter = name.rpartition(" ")
    assert star in STAR_NAMES


def test_same_seed_gives_same_name():
    assert default_world_name(rng=random.Random(42)) == default_world_name(rng=random.Random(42))


def test_middle_letters_come_up_more_often():
    counts: dict[str, int] = {}
    rng = random.Random(1234)
    for _ in range(1200):
        letter = default_world_name(rng=rng).rpartition(" ")[2]
        counts[letter] = counts.get(letter, 0) + 1
    middle = sum(counts.get(letter, 0) for letter in "cde")
    assert middle > sum(counts.values()) - middle


def test_letter_weights_cover_b_through_g():
    assert [letter for letter, _ in LETTER_WEIGHTS] == ["b", "c", "d", "e", "f", "g"]


def test_avoids_names_already_taken():
    first = default_world_name(rng=random.Random(9))
    second = default_world_name(used={first}, rng=random.Random(9))
    assert second != first
    assert NAME_RE.fullmatch(second), second


def test_numeral_fallback_when_every_combination_is_taken():
    used = {f"{star} {letter}" for star in STAR_NAMES for letter, _ in LETTER_WEIGHTS}
    name = default_world_name(used=used, rng=random.Random(3))
    assert re.fullmatch(r".+ [b-g] \d+", name), name


def test_star_names_are_not_known_exoplanet_hosts():
    # #640: a fresh world must never borrow a real planetary system's name. The full
    # planet-free verification is recorded in world_names.py's STAR_NAMES comment
    # (NASA Exoplanet Archive check, October 2026 -- it needs the network, so it
    # can't run here); this pins the previously-used hosts as gone for good.
    known_hosts = {
        "Kepler-22", "Kepler-62", "Kepler-69", "Kepler-186", "Kepler-442",
        "TRAPPIST-1", "Proxima Centauri", "Gliese 581", "Gliese 667 C", "Gliese 876",
        "LHS 1140", "K2-18", "TOI-700", "Ross 128", "Wolf 1061", "55 Cancri",
        "Upsilon Andromedae", "Mu Arae", "HD 40307", "HD 209458",
        "51 Pegasi", "47 Ursae Majoris", "Pollux",
        # Removed by the alias-aware re-check, pinned against re-adding: HD 191785
        # hosts TOI-5789 b-e, and HIP 79672 duplicates 18 Scorpii.
        "HD 191785", "HIP 79672",
    }
    assert not (set(STAR_NAMES) & known_hosts)

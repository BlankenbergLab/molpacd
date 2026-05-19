from __future__ import annotations

import pytest

from molpacd.names import choose_unused_residue_name, normalize_residue_name


def _all_candidate_names() -> set[str]:
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    return {
        f"{first}{second}{third}"
        for first in alphabet[:26]
        for second in alphabet
        for third in alphabet
    }


def test_choose_unused_residue_name_is_three_characters_and_unused() -> None:
    name = choose_unused_residue_name({"ALA", "HOH", "ABC"}, seed=1)

    assert len(name) == 3
    assert name not in {"ALA", "HOH", "ABC"}


def test_choose_unused_residue_name_falls_back_to_lexical_search() -> None:
    blocked = _all_candidate_names() - {"A00"}

    assert choose_unused_residue_name(blocked, seed=1) == "A00"


def test_choose_unused_residue_name_raises_when_exhausted() -> None:
    with pytest.raises(RuntimeError, match="unused"):
        choose_unused_residue_name(_all_candidate_names())


def test_normalize_residue_name_uppercases_valid_names() -> None:
    assert normalize_residue_name(" d1 ") == "D1"


def test_normalize_residue_name_rejects_invalid_names() -> None:
    with pytest.raises(ValueError):
        normalize_residue_name("TOOLONG")

from __future__ import annotations

import pytest

from molpacd.names import choose_unused_residue_name, normalize_residue_name


def test_choose_unused_residue_name_is_three_characters_and_unused() -> None:
    name = choose_unused_residue_name({"ALA", "HOH", "ABC"}, seed=1)

    assert len(name) == 3
    assert name not in {"ALA", "HOH", "ABC"}


def test_normalize_residue_name_rejects_invalid_names() -> None:
    with pytest.raises(ValueError):
        normalize_residue_name("TOOLONG")

from __future__ import annotations

import json
from pathlib import Path

import pytest

from molpacd.lipids import (
    AMINO_ACID_RESIDUE_NAMES,
    classify_lipid_atoms,
    load_lipid_fragments,
    load_lipid_residue_names,
)
from molpacd.models import AtomRecord


def _atom(**overrides: object) -> AtomRecord:
    values = {
        "record": "ATOM",
        "serial": 1,
        "name": "C1",
        "resname": "PA",
        "chain_id": "A",
        "res_seq": 1,
        "x": 0.0,
        "y": 0.0,
        "z": 0.0,
    }
    values.update(overrides)
    return AtomRecord(**values)


def test_load_lipid_fragments_returns_bundled_mapping() -> None:
    data = load_lipid_fragments()

    assert set(data["POPC"]) == {"OL", "PA", "PC"}
    assert "DOPE" in data


def test_load_lipid_residue_names_includes_fragments_and_full_names() -> None:
    names = load_lipid_residue_names()

    # Full lipid name and its head/tail fragment codes should both resolve.
    assert "POPC" in names
    assert {"OL", "PA", "PC"} <= names
    # Fragments only reachable through the bundled JSON (not the old
    # hardcoded script list) must also be present.
    assert {"AR", "DHA", "LAL", "PGR", "PH-", "SA", "SPM", "CHL"} <= names


def test_load_lipid_residue_names_accepts_explicit_override(tmp_path: Path) -> None:
    override = tmp_path / "custom.json"
    override.write_text(json.dumps({"XYZ": ["FOO", "BAR"]}), encoding="utf-8")

    names = load_lipid_residue_names(override)

    assert {"XYZ", "FOO", "BAR"} <= names
    # Fallback names are still merged in for robustness.
    assert "POPC" in names


def test_load_lipid_fragments_rejects_non_object_json(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(["not", "an", "object"]), encoding="utf-8")

    with pytest.raises(ValueError, match="JSON object"):
        load_lipid_fragments(bad)


def test_classify_lipid_atoms_groups_by_chain_and_residue_number() -> None:
    lipid_names = load_lipid_residue_names()
    atoms = [
        _atom(name="CA", resname="ALA", res_seq=1),
        _atom(name="C1", resname="PA", chain_id="A", res_seq=10),
        _atom(name="C2", resname="PA", chain_id="A", res_seq=10),
        _atom(name="C1", resname="OL", chain_id="A", res_seq=11),
        _atom(name="O", resname="WAT", chain_id="W", res_seq=1),
    ]

    protein, lipids, other = classify_lipid_atoms(atoms, lipid_names)

    assert [atom.resname for atom in protein] == ["ALA"]
    assert set(lipids.keys()) == {("A", 10), ("A", 11)}
    assert len(lipids[("A", 10)]) == 2
    assert [atom.resname for atom in other] == ["WAT"]


def test_amino_acid_residue_names_covers_common_variants() -> None:
    assert {"ALA", "HIE", "HID", "HIP", "CYX"} <= AMINO_ACID_RESIDUE_NAMES

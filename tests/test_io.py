from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from molpacd.io import _cif_value, atom_coordinates, read_structure, write_structure
from molpacd.models import AtomRecord, StructureData


def _atom(**overrides: object) -> AtomRecord:
    values = {
        "record": "ATOM",
        "serial": 1,
        "name": "CA",
        "resname": "GLY",
        "chain_id": "A",
        "res_seq": 1,
        "x": 1.0,
        "y": 2.0,
        "z": 3.0,
        "occupancy": 0.5,
        "bfactor": 10.0,
        "element": "C",
        "model": 1,
    }
    values.update(overrides)
    return AtomRecord(**values)


def test_read_and_write_reject_unsupported_formats(tmp_path: Path) -> None:
    unsupported_input = tmp_path / "structure.xyz"
    unsupported_input.write_text("", encoding="utf-8")
    structure = StructureData(atoms=[_atom()], source_format="pdb")

    with pytest.raises(ValueError, match="unsupported structure format"):
        read_structure(unsupported_input)

    with pytest.raises(ValueError, match="unsupported output format"):
        write_structure(structure, tmp_path / "structure.xyz")


def test_pdb_writer_preserves_header_metadata_and_models(tmp_path: Path) -> None:
    atom = _atom()
    structure = StructureData(
        atoms=[atom, replace(atom, serial=2, model=2, x=4.0)],
        source_format="pdb",
        header_lines=["HEADER    TEST STRUCTURE\n"],
        metadata={"resname": "DUM", "chain": "Z"},
    )
    path = tmp_path / "multi-model.pdb"

    write_structure(structure, path)
    text = path.read_text(encoding="utf-8")

    assert text.startswith("HEADER    TEST STRUCTURE\n")
    assert "REMARK MOLPACD CHAIN Z\n" in text
    assert "REMARK MOLPACD RESNAME DUM\n" in text
    assert text.count("MODEL") == 2
    assert text.count("ENDMDL") == 2

    reread = read_structure(path)
    assert reread.metadata == {"chain": "Z", "resname": "DUM"}
    assert reread.model_ids() == {1, 2}


def test_cif_writer_quotes_special_values_and_atom_coordinates(tmp_path: Path) -> None:
    atom = _atom(
        record="HETATM",
        name="_O",
        resname="D'M",
        chain_id="",
        element="",
        x=1.25,
        y=2.5,
        z=3.75,
    )
    structure = StructureData(atoms=[atom], source_format="cif")
    path = tmp_path / "quoted.cif"

    write_structure(structure, path)
    text = path.read_text(encoding="utf-8")

    assert "'_O'" in text
    assert "'D''M'" in text
    assert atom_coordinates([atom]) == [(1.25, 2.5, 3.75)]
    assert _cif_value("") == "."
    assert _cif_value("has space") == "'has space'"
    assert _cif_value("plain") == "plain"

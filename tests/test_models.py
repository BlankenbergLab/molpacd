from __future__ import annotations

import pytest

from molpacd.models import (
    AtomRecord,
    CapOptions,
    StructureData,
    normalize_atom_name,
    normalize_chain_id,
    normalize_element,
    validate_cap_options,
)


def _atom(**overrides: object) -> AtomRecord:
    values = {
        "record": "ATOM",
        "serial": 1,
        "name": "CA",
        "resname": "ala",
        "chain_id": "A",
        "res_seq": 1,
        "x": 1.0,
        "y": 2.0,
        "z": 3.0,
        "element": "C",
        "model": 2,
    }
    values.update(overrides)
    return AtomRecord(**values)


def test_atom_and_structure_helpers_copy_metadata() -> None:
    atom = _atom()
    structure = StructureData(
        atoms=[atom],
        source_format="pdb",
        header_lines=["HEADER    TEST\n"],
        metadata={"resname": "DUM"},
    )

    assert atom.coord.tolist() == [1.0, 2.0, 3.0]
    assert structure.residue_names() == {"ALA"}
    assert structure.model_ids() == {2}

    copied = structure.with_atoms([])
    assert copied.atoms == []
    assert copied.header_lines == ["HEADER    TEST\n"]
    assert copied.metadata == {"resname": "DUM"}
    assert copied.metadata is not structure.metadata

    overridden = structure.with_atoms([atom], metadata={"count": "0"})
    assert overridden.metadata == {"count": "0"}


@pytest.mark.parametrize(
    ("options", "match"),
    [
        (CapOptions(selection="water"), "selection"),
        (CapOptions(sides="middle"), "sides"),
        (CapOptions(spacing=0.0), "spacing"),
        (CapOptions(radius_scale=0.0), "radius scale"),
        (CapOptions(window=0.0), "window"),
        (CapOptions(min_atoms=2), "min_atoms"),
        (CapOptions(min_clearance=-0.1), "min_clearance"),
        (CapOptions(shared_radius="yes"), "shared_radius"),
        (CapOptions(atom_name="A B"), "atom name"),
        (CapOptions(element="C1"), "element"),
    ],
)
def test_validate_cap_options_rejects_invalid_values(options: CapOptions, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        validate_cap_options(options)


def test_normalizers_apply_defaults_and_reject_invalid_values() -> None:
    assert normalize_chain_id(" A ") == "A"
    assert normalize_chain_id("", default=None) is None
    assert normalize_chain_id(" ", default="Z") == "Z"
    assert normalize_atom_name(" ca ") == "CA"
    assert normalize_atom_name("", default=None) is None
    assert normalize_element("", "CA") == "C"
    assert normalize_element("cl", "CA") == "CL"

    with pytest.raises(ValueError, match="chain identifier"):
        normalize_chain_id("AB")
    with pytest.raises(ValueError, match="atom name"):
        normalize_atom_name("ABCDE")
    with pytest.raises(ValueError, match="atom name"):
        normalize_atom_name("A B")
    with pytest.raises(ValueError, match="element"):
        normalize_element("", "")
    with pytest.raises(ValueError, match="element"):
        normalize_element("XYZ", "CA")

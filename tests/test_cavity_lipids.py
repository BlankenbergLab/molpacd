from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path
from typing import List

import numpy as np
import pytest

from molpacd.cavity_lipids import remove_cavity_lipids
from molpacd.geometry import principal_axis_pca
from molpacd.models import AtomRecord, CavityLipidOptions, StructureData


def _barrel_ca_atoms(radius: float = 8.0) -> List[AtomRecord]:
    """A ring of CA atoms stacked along z, approximating a beta barrel."""
    atoms: List[AtomRecord] = []
    serial = 1
    for z in range(-15, 16, 5):
        for i in range(12):
            angle = 2 * math.pi * i / 12
            atoms.append(
                AtomRecord(
                    record="ATOM",
                    serial=serial,
                    name="CA",
                    resname="ALA",
                    chain_id="P",
                    res_seq=serial,
                    x=radius * math.cos(angle),
                    y=radius * math.sin(angle),
                    z=float(z),
                )
            )
            serial += 1
    return atoms


def _lipid_atom(res_seq: int, x: float, y: float, z: float, resname: str = "PC") -> AtomRecord:
    return AtomRecord(
        record="ATOM",
        serial=1000 + res_seq,
        name="C1",
        resname=resname,
        chain_id="L",
        res_seq=res_seq,
        x=x,
        y=y,
        z=z,
    )


def _structure(lipids: List[AtomRecord]) -> StructureData:
    return StructureData(atoms=[*_barrel_ca_atoms(), *lipids], source_format="pdb")


def test_removes_only_lipids_whose_center_falls_inside_the_cavity() -> None:
    inside = _lipid_atom(1, 0.0, 0.0, 0.0)  # dead center of the barrel
    outside_radially = _lipid_atom(2, 20.0, 0.0, 0.0)  # in the bilayer, not the pore
    outside_axially = _lipid_atom(3, 0.0, 0.0, 50.0)  # far past the barrel's z-range
    structure = _structure([inside, outside_radially, outside_axially])

    cleaned, result = remove_cavity_lipids(structure)

    assert result.total_lipid_residues == 3
    assert result.removed_lipid_residues == 1
    assert result.kept_lipid_residues == 2
    kept_keys = {(atom.chain_id, atom.res_seq) for atom in cleaned.atoms if atom.chain_id == "L"}
    assert kept_keys == {("L", 2), ("L", 3)}
    assert len(cleaned.atoms) == len(structure.atoms) - 1


def test_cavity_radius_override_widens_the_removed_region() -> None:
    outside_radially = _lipid_atom(2, 20.0, 0.0, 0.0)
    structure = _structure([outside_radially])

    _, default_result = remove_cavity_lipids(structure)
    assert default_result.removed_lipid_residues == 0

    _, widened_result = remove_cavity_lipids(structure, CavityLipidOptions(cavity_radius=25.0))
    assert widened_result.removed_lipid_residues == 1


def test_margin_option_is_recorded_and_expands_the_boundary() -> None:
    # Just past the default inner radius + margin (8 * 0.6 + 2.0 = 6.8).
    borderline = _lipid_atom(1, 7.5, 0.0, 0.0)
    structure = _structure([borderline])

    _, tight_result = remove_cavity_lipids(structure, CavityLipidOptions(margin=0.0))
    assert tight_result.removed_lipid_residues == 0

    _, loose_result = remove_cavity_lipids(structure, CavityLipidOptions(margin=5.0))
    assert loose_result.removed_lipid_residues == 1
    assert loose_result.margin == 5.0


def test_lipid_fragments_json_override_recognizes_additional_residue_names(
    tmp_path: Path,
) -> None:
    novel = _lipid_atom(1, 0.0, 0.0, 0.0, resname="ZZZ")
    structure = _structure([novel])

    _, default_result = remove_cavity_lipids(structure)
    assert default_result.total_lipid_residues == 0  # "ZZZ" isn't a known lipid name

    override = tmp_path / "custom.json"
    override.write_text(json.dumps({"WEIRD": ["ZZZ"]}), encoding="utf-8")

    _, override_result = remove_cavity_lipids(
        structure, CavityLipidOptions(lipid_fragments_json=str(override))
    )
    assert override_result.total_lipid_residues == 1
    assert override_result.removed_lipid_residues == 1


def test_invalid_options_are_rejected() -> None:
    structure = _structure([])

    with pytest.raises(ValueError, match="margin"):
        remove_cavity_lipids(structure, CavityLipidOptions(margin=-1.0))

    with pytest.raises(ValueError, match="cavity_radius"):
        remove_cavity_lipids(structure, CavityLipidOptions(cavity_radius=0.0))


def test_requires_protein_atoms() -> None:
    structure = StructureData(atoms=[_lipid_atom(1, 0.0, 0.0, 0.0)], source_format="pdb")

    with pytest.raises(ValueError, match="no protein atoms"):
        remove_cavity_lipids(structure)


def test_rejects_multi_model_structures() -> None:
    structure = _structure([_lipid_atom(1, 0.0, 0.0, 0.0)])
    multimodel = structure.with_atoms([replace(structure.atoms[0], model=2), *structure.atoms[1:]])

    with pytest.raises(ValueError, match="multi-model"):
        remove_cavity_lipids(multimodel)


def test_preserves_existing_structure_metadata() -> None:
    structure = _structure([_lipid_atom(1, 0.0, 0.0, 0.0)])
    tagged = structure.with_atoms(structure.atoms, metadata={"resname": "DUM"})

    cleaned, _ = remove_cavity_lipids(tagged)

    assert cleaned.metadata == {"resname": "DUM"}


def test_principal_axis_flips_a_raw_negative_z_eigenvector() -> None:
    # These points' raw covariance eigenvector (before the z>=0 convention
    # is enforced) points toward negative z; verifies the flip is applied.
    coords = np.array(
        [
            [0.33962000824864264, 0.42377135285334727, 3.7122741773625885],
            [0.3827571602707609, 0.3194142202523809, -3.5891330853862047],
            [-1.9016352983759945, -0.10891472790742324, -8.037318485206766],
            [1.0801634125378852, -0.2887665059953775, 0.8347535610700986],
            [-0.8496059556101431, -0.5106224678981267, -0.11533061686586776],
        ]
    )

    _, axis = principal_axis_pca(coords)

    assert axis[2] >= 0


def test_falls_back_to_all_protein_atoms_when_no_ca_atoms_present() -> None:
    backbone_only = [replace(atom, name="N") for atom in _barrel_ca_atoms()]
    structure = StructureData(
        atoms=[*backbone_only, _lipid_atom(1, 0.0, 0.0, 0.0)], source_format="pdb"
    )

    _, result = remove_cavity_lipids(structure)

    assert result.removed_lipid_residues == 1

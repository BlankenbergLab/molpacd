from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from molpacd import (
    CapOptions,
    add_caps,
    analyze_structure,
    read_structure,
    remove_caps,
    write_structure,
)

FIXTURE = Path(__file__).parent / "data" / "2iww.pdb"


def test_read_and_analyze_pdb_fixture() -> None:
    structure = read_structure(FIXTURE)

    assert len(structure.atoms) > 2000
    analysis = analyze_structure(structure, CapOptions(axis="z"))
    assert analysis.negative.atom_count >= 12
    assert analysis.positive.atom_count >= 12
    assert analysis.negative.radius > 5.0
    assert analysis.positive.radius > 5.0


def test_add_and_remove_caps_round_trip(tmp_path: Path) -> None:
    structure = read_structure(FIXTURE)
    capped, result = add_caps(structure, CapOptions(axis="z", seed=7))

    assert result.added_count > 0
    assert result.resname not in structure.residue_names()
    assert capped.metadata["resname"] == result.resname

    capped_path = tmp_path / "capped.pdb"
    write_structure(capped, capped_path)
    reread = read_structure(capped_path)

    decapped, remove_result = remove_caps(reread)
    assert remove_result.removed_count == result.added_count
    assert len(decapped.atoms) == len(structure.atoms)


def test_remove_uses_metadata_ranges_without_removing_matching_decoys() -> None:
    structure = read_structure(FIXTURE)
    capped, result = add_caps(structure, CapOptions(axis="z", seed=7))
    assert result.serial_end is not None
    assert result.res_seq_end is not None

    decoy = replace(
        capped.atoms[0],
        record="HETATM",
        serial=result.serial_end + 100,
        name=result.atom_name,
        resname=result.resname,
        chain_id=result.chain,
        res_seq=result.res_seq_end + 100,
    )
    with_decoy = capped.with_atoms([*capped.atoms, decoy], metadata=capped.metadata)

    decapped, remove_result = remove_caps(with_decoy)

    assert remove_result.removed_count == result.added_count
    assert any(atom.serial == decoy.serial for atom in decapped.atoms)


def test_remove_metadata_override_requires_force() -> None:
    structure = read_structure(FIXTURE)
    capped, _ = add_caps(structure, CapOptions(axis="z", seed=7))

    with pytest.raises(ValueError, match="requires --force"):
        remove_caps(capped, chain="Q")

    with pytest.raises(ValueError, match="requires --force"):
        remove_caps(capped, chain="")

    decapped, remove_result = remove_caps(capped, chain="Q", force=True)
    assert remove_result.removed_count == 0
    assert len(decapped.atoms) == len(capped.atoms)


def test_write_and_read_mmcif_round_trip(tmp_path: Path) -> None:
    structure = read_structure(FIXTURE)
    cif_path = tmp_path / "fixture.cif"

    write_structure(structure, cif_path)
    reread = read_structure(cif_path)

    assert len(reread.atoms) == len(structure.atoms)
    assert reread.atoms[0].resname == structure.atoms[0].resname


def test_mmcif_preserves_molpacd_metadata_for_removal(tmp_path: Path) -> None:
    structure = read_structure(FIXTURE)
    capped, result = add_caps(structure, CapOptions(axis="z", seed=7))
    cif_path = tmp_path / "capped.cif"

    write_structure(capped, cif_path)
    reread = read_structure(cif_path)
    decapped, remove_result = remove_caps(reread)

    assert reread.metadata["count"] == str(result.added_count)
    assert remove_result.removed_count == result.added_count
    assert len(decapped.atoms) == len(structure.atoms)


def test_remove_requires_force_without_metadata() -> None:
    structure = read_structure(FIXTURE)

    with pytest.raises(ValueError, match="requires --force"):
        remove_caps(structure, resname="DUM")


def test_invalid_options_are_rejected() -> None:
    structure = read_structure(FIXTURE)

    with pytest.raises(ValueError, match="chain identifier"):
        add_caps(structure, CapOptions(axis="z", chain="ZZ"))

    with pytest.raises(ValueError, match="min_atoms"):
        analyze_structure(structure, CapOptions(axis="z", min_atoms=2))

    with pytest.raises(ValueError, match="min_clearance"):
        add_caps(structure, CapOptions(axis="z", min_clearance=-1.0))


def test_analyze_and_add_reject_multi_model_structures() -> None:
    structure = read_structure(FIXTURE)
    multimodel = structure.with_atoms(
        [replace(structure.atoms[0], model=2), *structure.atoms[1:]],
        metadata=structure.metadata,
    )

    with pytest.raises(ValueError, match="multi-model"):
        analyze_structure(multimodel, CapOptions(axis="z"))

    with pytest.raises(ValueError, match="multi-model"):
        add_caps(multimodel, CapOptions(axis="z"))


def test_cli_add_json_dry_run(tmp_path: Path) -> None:
    capped = tmp_path / "capped.pdb"

    add_result = subprocess.run(
        [
            sys.executable,
            "-m",
            "molpacd",
            "add",
            str(FIXTURE),
            "-o",
            str(capped),
            "--axis",
            "z",
            "--seed",
            "9",
            "--dry-run",
            "--json",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(add_result.stdout)
    assert payload["added_count"] > 0
    assert payload["dry_run"] is True
    assert not capped.exists()


def test_cli_add_and_remove(tmp_path: Path) -> None:
    capped = tmp_path / "capped.pdb"
    decapped = tmp_path / "decapped.pdb"

    add_result = subprocess.run(
        [
            sys.executable,
            "-m",
            "molpacd",
            "add",
            str(FIXTURE),
            "-o",
            str(capped),
            "--axis",
            "z",
            "--seed",
            "9",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "total added atoms" in add_result.stdout
    assert capped.exists()

    remove_result = subprocess.run(
        [
            sys.executable,
            "-m",
            "molpacd",
            "remove",
            str(capped),
            "-o",
            str(decapped),
            "--json",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(remove_result.stdout)
    assert payload["removed_count"] > 0
    assert decapped.exists()

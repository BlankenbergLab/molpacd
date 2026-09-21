from __future__ import annotations

import math
from pathlib import Path
from typing import List

import pytest

from molpacd.cavity_lipids_pdb import remove_cavity_lipids_pdb
from molpacd.models import CavityLipidOptions


def _pdb_line(
    record: str,
    serial: int,
    name: str,
    resname: str,
    chain: str,
    res_seq: int,
    x: float,
    y: float,
    z: float,
    tail: str = "",
) -> str:
    return (
        f"{record:<6}"
        f"{serial:>5d} "
        f"{name:<4}"
        f" "
        f"{resname:>3} "
        f"{chain}"
        f"{res_seq:>4d}"
        f"    "
        f"{x:>8.3f}"
        f"{y:>8.3f}"
        f"{z:>8.3f}"
        f"{1.00:>6.2f}"
        f"{0.00:>6.2f}"
        f"{tail}\n"
    )


def _barrel_and_lipid_lines() -> List[str]:
    lines: List[str] = ["HEADER    SYNTHETIC BARREL\n", "CRYST1    1    1    1  90  90  90 P 1\n"]
    serial = 1
    for z in range(-15, 16, 5):
        for i in range(12):
            angle = 2 * math.pi * i / 12
            lines.append(
                _pdb_line(
                    "ATOM",
                    serial,
                    "CA",
                    "ALA",
                    "P",
                    serial,
                    8.0 * math.cos(angle),
                    8.0 * math.sin(angle),
                    float(z),
                )
            )
            serial += 1
    # A lipid dead center of the barrel (should be removed), with a
    # packmol-memgen-style segment-id tail field past the standard columns.
    lines.append(_pdb_line("ATOM", 9001, "C1", "PC", "L", 1, 0.0, 0.0, 0.0, tail="      MEMB     "))
    # A lipid out in the bilayer (should be kept), same tail styling.
    lines.append(
        _pdb_line("ATOM", 9002, "C1", "PC", "L", 2, 20.0, 0.0, 0.0, tail="      MEMB     ")
    )
    lines.append("END\n")
    return lines


def _write(path: Path, lines: List[str]) -> None:
    path.write_text("".join(lines))


def test_preserves_kept_lines_byte_for_byte(tmp_path: Path) -> None:
    input_path = tmp_path / "packed.pdb"
    output_path = tmp_path / "cleaned.pdb"
    lines = _barrel_and_lipid_lines()
    _write(input_path, lines)

    result = remove_cavity_lipids_pdb(input_path, output_path)

    assert result.total_lipid_residues == 2
    assert result.removed_lipid_residues == 1

    output_text = output_path.read_text()
    # Every kept original line (header, CRYST1, protein, the kept lipid, END)
    # must appear byte-for-byte, including the packmol-style segid tail.
    for line in lines:
        if "L   1" in line and "PC " in line:
            assert line not in output_text
        else:
            assert line in output_text


def test_inserts_remark_before_first_coordinate_line(tmp_path: Path) -> None:
    input_path = tmp_path / "packed.pdb"
    output_path = tmp_path / "cleaned.pdb"
    _write(input_path, _barrel_and_lipid_lines())

    remove_cavity_lipids_pdb(input_path, output_path)

    output_lines = output_path.read_text().splitlines()
    header_index = output_lines.index("HEADER    SYNTHETIC BARREL")
    remark_index = next(i for i, line in enumerate(output_lines) if "REMARK" in line)
    first_atom_index = next(i for i, line in enumerate(output_lines) if line.startswith("ATOM"))

    assert header_index < remark_index < first_atom_index
    assert "Removed 1 of 2 lipid residues" in output_lines[remark_index + 1]


def test_reused_chain_letters_across_ter_separated_segments_are_not_dropped(
    tmp_path: Path,
) -> None:
    # Mirrors real packmol-memgen output: many independent lipid molecules,
    # each terminated with TER, reusing a small set of chain letters and
    # residue numbers. Biopython's structure builder can silently drop such
    # atoms (see the io.py limitation this module works around); this raw
    # parser must account for every ATOM/HETATM line regardless.
    lines = _barrel_and_lipid_lines()[:-1]  # drop the trailing END for now
    extra: List[str] = []
    for i in range(3, 8):
        extra.append(_pdb_line("ATOM", 9000 + i, "C1", "PC", "L", i, 30.0 + i, 0.0, 0.0))
        extra.append("TER\n")
    lines = lines + extra + ["END\n"]
    input_atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    input_path = tmp_path / "packed.pdb"
    output_path = tmp_path / "cleaned.pdb"
    _write(input_path, lines)

    result = remove_cavity_lipids_pdb(input_path, output_path)

    # 2 original lipid residues + 5 more distinct ones (unique resseq each).
    assert result.total_lipid_residues == 7
    assert result.removed_lipid_residues == 1

    output_atom_lines = [
        line
        for line in output_path.read_text().splitlines(keepends=True)
        if line.startswith(("ATOM", "HETATM"))
    ]
    # Every input atom line is accounted for: either still present in the
    # output, or explicitly removed as part of the one cavity-bound residue.
    assert len(output_atom_lines) == len(input_atom_lines) - 1
    assert set(output_atom_lines) == set(input_atom_lines) - {
        line for line in input_atom_lines if "L   1" in line and " PC " in line
    }


def test_rejects_non_pdb_extensions(tmp_path: Path) -> None:
    input_path = tmp_path / "packed.cif"
    input_path.write_text("data_x\n")
    output_path = tmp_path / "cleaned.pdb"

    with pytest.raises(ValueError, match="expected a PDB file"):
        remove_cavity_lipids_pdb(input_path, output_path)

    pdb_input = tmp_path / "packed.pdb"
    _write(pdb_input, _barrel_and_lipid_lines())
    with pytest.raises(ValueError, match="expected a PDB file"):
        remove_cavity_lipids_pdb(pdb_input, tmp_path / "cleaned.cif")


def test_rejects_malformed_atom_records(tmp_path: Path) -> None:
    input_path = tmp_path / "bad.pdb"
    input_path.write_text("ATOM      1  CA  ALA A   1     not-a-number\n")

    with pytest.raises(ValueError, match="malformed PDB"):
        remove_cavity_lipids_pdb(input_path, tmp_path / "out.pdb")


def test_rejects_files_with_no_atom_records(tmp_path: Path) -> None:
    input_path = tmp_path / "empty.pdb"
    input_path.write_text("HEADER    EMPTY\nEND\n")

    with pytest.raises(ValueError, match="no ATOM/HETATM records"):
        remove_cavity_lipids_pdb(input_path, tmp_path / "out.pdb")


def test_requires_protein_atoms(tmp_path: Path) -> None:
    input_path = tmp_path / "lipids-only.pdb"
    _write(input_path, [_pdb_line("ATOM", 1, "C1", "PC", "L", 1, 0.0, 0.0, 0.0)])

    with pytest.raises(ValueError, match="no protein atoms"):
        remove_cavity_lipids_pdb(input_path, tmp_path / "out.pdb")


def test_rejects_multi_model_structures(tmp_path: Path) -> None:
    lines = ["MODEL        1\n", *_barrel_and_lipid_lines()[:-1], "ENDMDL\n"]
    lines += ["MODEL        2\n", *_barrel_and_lipid_lines()[2:-1], "ENDMDL\n", "END\n"]
    input_path = tmp_path / "multimodel.pdb"
    _write(input_path, lines)

    with pytest.raises(ValueError, match="multi-model"):
        remove_cavity_lipids_pdb(input_path, tmp_path / "out.pdb")


def test_ter_records_are_never_touched_even_for_removed_residues(tmp_path: Path) -> None:
    # A TER directly after the removed lipid's atoms, and one after the kept
    # lipid's atoms: both must survive untouched, exactly where they were in
    # the input, regardless of which residue's atoms were removed.
    base_lines = _barrel_and_lipid_lines()
    removed_atom_index = next(
        i for i, line in enumerate(base_lines) if "L   1" in line and " PC " in line
    )
    kept_atom_index = next(
        i for i, line in enumerate(base_lines) if "L   2" in line and " PC " in line
    )
    lines = [*base_lines[: removed_atom_index + 1], "TER\n", *base_lines[removed_atom_index + 1 :]]
    kept_atom_index += 1  # shifted by the TER inserted above
    lines = [*lines[: kept_atom_index + 1], "TER\n", *lines[kept_atom_index + 1 :]]
    input_path = tmp_path / "packed.pdb"
    output_path = tmp_path / "cleaned.pdb"
    _write(input_path, lines)
    ter_count_in = sum(1 for line in lines if line.startswith("TER"))
    assert ter_count_in == 2

    result = remove_cavity_lipids_pdb(input_path, output_path)

    assert result.removed_lipid_residues == 1
    output_lines = output_path.read_text().splitlines(keepends=True)
    ter_count_out = sum(1 for line in output_lines if line.startswith("TER"))
    assert ter_count_out == 2
    for line in lines:
        if line.startswith("TER"):
            assert line in output_lines


def test_cavity_radius_and_margin_options_are_applied(tmp_path: Path) -> None:
    input_path = tmp_path / "packed.pdb"
    output_path = tmp_path / "cleaned.pdb"
    _write(input_path, _barrel_and_lipid_lines())

    result = remove_cavity_lipids_pdb(
        input_path, output_path, CavityLipidOptions(cavity_radius=25.0, margin=0.0)
    )

    # With a widened cavity radius, the bilayer lipid at radius 20 is also removed.
    assert result.removed_lipid_residues == 2

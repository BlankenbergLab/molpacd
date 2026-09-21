from __future__ import annotations

import json
import math
import runpy
import sys
from pathlib import Path
from typing import List

import pytest

from molpacd.cli import main
from molpacd.io import write_structure
from molpacd.models import AtomRecord, StructureData

FIXTURE = Path(__file__).parent / "data" / "2iww.pdb"


def _cavity_lipid_fixture(path: Path) -> None:
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
                    x=8.0 * math.cos(angle),
                    y=8.0 * math.sin(angle),
                    z=float(z),
                )
            )
            serial += 1
    atoms.append(
        AtomRecord(
            record="ATOM",
            serial=9001,
            name="C1",
            resname="PC",
            chain_id="L",
            res_seq=1,
            x=0.0,
            y=0.0,
            z=0.0,
        )
    )
    atoms.append(
        AtomRecord(
            record="ATOM",
            serial=9002,
            name="C1",
            resname="PC",
            chain_id="L",
            res_seq=2,
            x=20.0,
            y=0.0,
            z=0.0,
        )
    )
    write_structure(StructureData(atoms=atoms, source_format="pdb"), path)


def test_cli_analyze_json(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["analyze", str(FIXTURE), "--axis", "z", "--json"]) == 0

    payload = json.loads(capsys.readouterr().out)
    assert payload["axis"] == [0.0, 0.0, 1.0]
    assert payload["negative"]["atom_count"] >= 12
    assert payload["positive"]["atom_count"] >= 12


def test_cli_analyze_text(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["analyze", str(FIXTURE), "--axis", "z"]) == 0

    output = capsys.readouterr().out
    assert "axis: (0.000, 0.000, 1.000)" in output
    assert "negative: atoms=" in output
    assert "positive: atoms=" in output


def test_cli_add_text_dry_run(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    capped = tmp_path / "capped.pdb"

    assert (
        main(
            [
                "add",
                str(FIXTURE),
                "-o",
                str(capped),
                "--axis",
                "z",
                "--seed",
                "9",
                "--dry-run",
            ]
        )
        == 0
    )

    output = capsys.readouterr().out
    assert "total added atoms" in output
    assert "wrote" not in output
    assert not capped.exists()


def test_cli_add_json_and_remove_outputs(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    capped = tmp_path / "capped.pdb"
    decapped_text = tmp_path / "decapped-text.pdb"
    decapped_json = tmp_path / "decapped-json.pdb"

    assert (
        main(
            [
                "add",
                str(FIXTURE),
                "-o",
                str(capped),
                "--axis",
                "z",
                "--seed",
                "9",
                "--json",
            ]
        )
        == 0
    )
    add_payload = json.loads(capsys.readouterr().out)
    assert add_payload["added_count"] > 0
    assert add_payload["wrote"] == str(capped)
    assert capped.exists()

    assert main(["remove", str(capped), "-o", str(decapped_text)]) == 0
    remove_text = capsys.readouterr().out
    assert "removed " in remove_text
    assert "wrote" in remove_text
    assert decapped_text.exists()

    assert main(["remove", str(capped), "-o", str(decapped_json), "--json"]) == 0
    remove_payload = json.loads(capsys.readouterr().out)
    assert remove_payload["removed_count"] == add_payload["added_count"]
    assert remove_payload["wrote"] == str(decapped_json)
    assert decapped_json.exists()


def test_cli_remove_lipids_text_and_json(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fixture = tmp_path / "packed.pdb"
    _cavity_lipid_fixture(fixture)
    cleaned_text = tmp_path / "cleaned-text.pdb"
    cleaned_json = tmp_path / "cleaned-json.pdb"

    assert main(["remove-lipids", str(fixture), "-o", str(cleaned_text)]) == 0
    output = capsys.readouterr().out
    assert "removed 1 of 2 lipid residues" in output
    assert "wrote" in output
    assert cleaned_text.exists()

    assert main(["remove-lipids", str(fixture), "-o", str(cleaned_json), "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["total_lipid_residues"] == 2
    assert payload["removed_lipid_residues"] == 1
    assert payload["kept_lipid_residues"] == 1
    assert payload["wrote"] == str(cleaned_json)
    assert cleaned_json.exists()


def test_cli_remove_lipids_cavity_radius_and_margin_overrides(tmp_path: Path) -> None:
    fixture = tmp_path / "packed.pdb"
    _cavity_lipid_fixture(fixture)
    cleaned = tmp_path / "cleaned.pdb"

    assert (
        main(
            [
                "remove-lipids",
                str(fixture),
                "-o",
                str(cleaned),
                "--cavity-radius",
                "25",
                "--margin",
                "0",
                "--json",
            ]
        )
        == 0
    )


def test_cli_reports_errors_and_debug_reraises(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    missing = tmp_path / "missing.pdb"

    assert main(["analyze", str(missing)]) == 1
    assert "molpacd: error:" in capsys.readouterr().err

    with pytest.raises(FileNotFoundError):
        main(["--debug", "analyze", str(missing)])


def test_cli_version_and_module_entrypoint(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as direct_exit:
        main(["--version"])
    assert direct_exit.value.code == 0
    assert "molpacd " in capsys.readouterr().out

    monkeypatch.setattr(sys, "argv", ["molpacd", "--version"])
    with pytest.raises(SystemExit) as module_exit:
        runpy.run_module("molpacd.__main__", run_name="__main__")
    assert module_exit.value.code == 0
    assert "molpacd " in capsys.readouterr().out

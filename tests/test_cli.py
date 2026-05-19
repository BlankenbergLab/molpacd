from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path

import pytest

from molpacd.cli import main

FIXTURE = Path(__file__).parent / "data" / "2iww.pdb"


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

from __future__ import annotations

import numpy as np
import pytest

from molpacd.geometry import (
    coords_from_atoms,
    filter_collisions,
    generate_cap_disk,
    parse_axis,
    select_atoms,
)
from molpacd.models import AtomRecord, CapOptions


def _atom(name: str, record: str = "ATOM") -> AtomRecord:
    return AtomRecord(
        record=record,
        serial=1,
        name=name,
        resname="GLY",
        chain_id="A",
        res_seq=1,
        x=0.0,
        y=0.0,
        z=0.0,
    )


def test_coords_from_atoms_and_selection_filters() -> None:
    atoms = [
        _atom("CA"),
        _atom("CA"),
        _atom("CA"),
        _atom("N"),
        _atom("C"),
        _atom("O"),
        _atom("CB"),
        _atom("CA", record="HETATM"),
    ]

    assert coords_from_atoms(atoms).shape == (8, 3)
    assert [atom.name for atom in select_atoms(atoms, CapOptions(selection="ca"))] == [
        "CA",
        "CA",
        "CA",
    ]
    assert [atom.name for atom in select_atoms(atoms, CapOptions(selection="backbone"))] == [
        "CA",
        "CA",
        "CA",
        "N",
        "C",
        "O",
    ]
    assert len(select_atoms(atoms, CapOptions(selection="all"))) == 8


def test_coords_and_selection_reject_empty_results() -> None:
    with pytest.raises(ValueError, match="no atoms selected"):
        coords_from_atoms([])

    with pytest.raises(ValueError, match="at least 3 required"):
        select_atoms([_atom("CB")], CapOptions(selection="ca"))


def test_parse_axis_named_auto_vector_and_errors() -> None:
    coords = np.array(
        [
            [-2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ],
        dtype=float,
    )

    assert np.allclose(parse_axis("auto", coords), np.array([1.0, 0.0, 0.0]))
    assert np.allclose(parse_axis("x", coords), np.array([1.0, 0.0, 0.0]))
    assert np.allclose(parse_axis("y", coords), np.array([0.0, 1.0, 0.0]))
    assert np.allclose(parse_axis("0,0,2", coords), np.array([0.0, 0.0, 1.0]))

    with pytest.raises(ValueError, match="3-value vector"):
        parse_axis("1,2", coords)
    with pytest.raises(ValueError, match="numeric"):
        parse_axis("a,b,c", coords)
    with pytest.raises(ValueError, match="must not be zero"):
        parse_axis("0,0,0", coords)


def test_generate_cap_disk_is_planar_and_centered() -> None:
    positions = generate_cap_disk(
        centroid=(1.0, 2.0, 3.0),
        axis=(0.0, 0.0, 1.0),
        radius=2.0,
        spacing=1.0,
        side_sign=1.0,
    )

    assert positions.shape[1] == 3
    assert np.allclose(positions[:, 2], 3.0)
    assert np.allclose(np.mean(positions, axis=0), np.array([1.0, 2.0, 3.0]))


def test_generate_cap_disk_handles_non_z_axis_and_rejects_invalid_inputs() -> None:
    positions = generate_cap_disk(
        centroid=(0.0, 0.0, 0.0),
        axis=(1.0, 0.0, 0.0),
        radius=2.0,
        spacing=1.0,
        side_sign=1.0,
    )

    assert np.allclose(positions[:, 0], 0.0)

    with pytest.raises(ValueError, match="spacing"):
        generate_cap_disk((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 2.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="radius"):
        generate_cap_disk((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.0, 1.0, 1.0)


def test_filter_collisions_skips_positions_too_close_to_existing_atoms() -> None:
    positions = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]], dtype=float)
    existing = np.array([[0.5, 0.0, 0.0]], dtype=float)

    accepted, skipped = filter_collisions(positions, existing, min_clearance=1.0)

    assert skipped == 1
    assert accepted.tolist() == [[5.0, 0.0, 0.0]]


def test_filter_collisions_handles_disabled_and_all_colliding_cases() -> None:
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]], dtype=float)
    existing = np.array([[0.0, 0.0, 0.0]], dtype=float)

    accepted, skipped = filter_collisions(positions, existing, min_clearance=0.0)
    assert skipped == 0
    assert accepted is positions

    accepted, skipped = filter_collisions(positions, np.empty((0, 3)), min_clearance=1.0)
    assert skipped == 0
    assert accepted is positions

    accepted, skipped = filter_collisions(positions, existing, min_clearance=1.0)
    assert skipped == 2
    assert accepted.shape == (0, 3)

from __future__ import annotations

import contextlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
from numpy.typing import NDArray

from molpacd.geometry import barrel_radii, point_inside_cylinder, principal_axis_pca
from molpacd.lipids import AMINO_ACID_RESIDUE_NAMES, load_lipid_residue_names
from molpacd.models import CavityLipidOptions, CavityLipidResult, validate_cavity_lipid_options

_PDB_EXTENSIONS = {".pdb", ".ent"}


@dataclass(frozen=True)
class _RawAtom:
    line_index: int
    model: int
    name: str
    resname: str
    chain_id: str
    res_seq: int
    x: float
    y: float
    z: float


def remove_cavity_lipids_pdb(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    options: Optional[CavityLipidOptions] = None,
) -> CavityLipidResult:
    """Remove cavity lipids from a PDB file, preserving its original formatting.

    Unlike `molpacd.cavity_lipids.remove_cavity_lipids`, this reads and
    writes raw PDB text directly instead of going through molpacd's
    Biopython-based structure I/O. Every line that is not part of a removed
    lipid residue is copied through byte-for-byte from the input file, so
    column widths, TER records, segment IDs, and any other fields molpacd's
    generic writer does not reproduce are preserved exactly. This also
    avoids residue-numbering edge cases in large multi-chain packmol-memgen
    output where Biopython's structure builder can silently drop atoms.
    """
    options = options or CavityLipidOptions()
    validate_cavity_lipid_options(options)

    input_path = Path(input_path)
    output_path = Path(output_path)
    if input_path.suffix.lower() not in _PDB_EXTENSIONS:
        raise ValueError(f"unsupported input format for {input_path}; expected a PDB file")
    if output_path.suffix.lower() not in _PDB_EXTENSIONS:
        raise ValueError(f"unsupported output format for {output_path}; expected a PDB file")

    lines = input_path.read_text().splitlines(keepends=True)
    atoms = _parse_pdb_atoms(lines)
    if not atoms:
        raise ValueError(f"no ATOM/HETATM records found in {input_path}")
    if len({atom.model for atom in atoms}) > 1:
        raise ValueError(
            "multi-model structures are not supported for cavity lipid removal; "
            "split the structure into a single model first"
        )

    lipid_residue_names = load_lipid_residue_names(options.lipid_fragments_json)
    protein_atoms, lipid_residues = _classify_raw_atoms(atoms, lipid_residue_names)
    if not protein_atoms:
        raise ValueError("no protein atoms found; cannot determine a cavity")

    ca_atoms = [atom for atom in protein_atoms if atom.name.strip().upper() == "CA"]
    axis_atoms = ca_atoms if ca_atoms else protein_atoms
    axis_coords = _coords(axis_atoms)

    center, axis = principal_axis_pca(axis_coords)
    inner_radius, outer_radius, z_min, z_max = barrel_radii(axis_coords, center, axis)
    if options.cavity_radius is not None:
        inner_radius = options.cavity_radius

    removed_keys: Set[Tuple[str, int]] = set()
    removed_line_indices: Set[int] = set()
    for key, residue_atoms in lipid_residues.items():
        residue_center = np.mean(_coords(residue_atoms), axis=0)
        if point_inside_cylinder(
            residue_center, center, axis, inner_radius, z_min, z_max, options.margin
        ):
            removed_keys.add(key)
            removed_line_indices.update(atom.line_index for atom in residue_atoms)

    result = CavityLipidResult(
        total_lipid_residues=len(lipid_residues),
        removed_lipid_residues=len(removed_keys),
        kept_lipid_residues=len(lipid_residues) - len(removed_keys),
        center=(float(center[0]), float(center[1]), float(center[2])),
        axis=(float(axis[0]), float(axis[1]), float(axis[2])),
        inner_radius=float(inner_radius),
        outer_radius=float(outer_radius),
        z_min=float(z_min),
        z_max=float(z_max),
        margin=options.margin,
    )
    _write_pdb(lines, removed_line_indices, result, output_path)
    return result


def _parse_pdb_atoms(lines: Sequence[str]) -> List[_RawAtom]:
    atoms: List[_RawAtom] = []
    model = 1
    for index, line in enumerate(lines):
        if line.startswith("MODEL"):
            with contextlib.suppress(ValueError):
                model = int(line[10:14])
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain_id = line[21].strip() or " "
            res_seq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"malformed PDB ATOM/HETATM record at line {index + 1}: {line!r}"
            ) from exc
        atoms.append(
            _RawAtom(
                line_index=index,
                model=model,
                name=name,
                resname=resname,
                chain_id=chain_id,
                res_seq=res_seq,
                x=x,
                y=y,
                z=z,
            )
        )
    return atoms


def _classify_raw_atoms(
    atoms: Sequence[_RawAtom], lipid_residue_names: Set[str]
) -> Tuple[List[_RawAtom], Dict[Tuple[str, int], List[_RawAtom]]]:
    protein_atoms: List[_RawAtom] = []
    lipid_residues: Dict[Tuple[str, int], List[_RawAtom]] = defaultdict(list)
    for atom in atoms:
        resname = atom.resname.strip().upper()
        if resname in AMINO_ACID_RESIDUE_NAMES:
            protein_atoms.append(atom)
        elif resname in lipid_residue_names:
            lipid_residues[(atom.chain_id, atom.res_seq)].append(atom)
    return protein_atoms, dict(lipid_residues)


def _coords(atoms: Sequence[_RawAtom]) -> NDArray[np.float64]:
    return np.array([[atom.x, atom.y, atom.z] for atom in atoms], dtype=float)


def _write_pdb(
    lines: Sequence[str],
    removed_line_indices: Set[int],
    result: CavityLipidResult,
    output_path: Path,
) -> None:
    remark_lines = [
        "REMARK   Processed by molpacd remove-lipids\n",
        f"REMARK   Removed {result.removed_lipid_residues} of "
        f"{result.total_lipid_residues} lipid residues from the cavity\n",
    ]
    output_lines: List[str] = []
    inserted = False
    for index, line in enumerate(lines):
        # `atoms` is guaranteed non-empty by the caller, so some line here
        # always starts with ATOM/HETATM/MODEL and this always fires once.
        if not inserted and line.startswith(("ATOM", "HETATM", "MODEL")):
            output_lines.extend(remark_lines)
            inserted = True
        if index in removed_line_indices:
            continue
        output_lines.append(line)
    _write_text_atomic(output_path, "".join(output_lines))


def _write_text_atomic(path: Path, text: str) -> None:
    temp_path: Optional[Path] = None
    try:
        with NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temp_path = Path(handle.name)
            handle.write(text)
        temp_path.replace(path)
    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()


__all__ = ["remove_cavity_lipids_pdb"]

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from molpacd.geometry import (
    barrel_radii,
    coord_tuple,
    coords_from_atoms,
    point_inside_cylinder,
    principal_axis_pca,
)
from molpacd.lipids import classify_lipid_atoms, load_lipid_residue_names
from molpacd.models import (
    CavityLipidOptions,
    CavityLipidResult,
    StructureData,
    validate_cavity_lipid_options,
)


def remove_cavity_lipids(
    structure: StructureData,
    options: Optional[CavityLipidOptions] = None,
) -> Tuple[StructureData, CavityLipidResult]:
    """Remove lipids packed inside a protein cavity (e.g. a beta barrel).

    Packmol-memgen can pack lipid tails into a protein's interior cavity
    instead of only around its exterior. This estimates the cavity from the
    protein's CA atoms and removes lipid residues whose center of mass falls
    inside it, leaving lipids in the surrounding bilayer untouched.

    This operates on an already-parsed `StructureData` (PDB or mmCIF, via
    `molpacd.io.read_structure`), so writing the result back out goes through
    molpacd's generic structure writer and does not reproduce the exact
    formatting of the original input file. To preserve a PDB file's original
    formatting exactly, use `remove_cavity_lipids_pdb` instead.
    """
    options = options or CavityLipidOptions()
    validate_cavity_lipid_options(options)
    _ensure_single_model(structure)

    lipid_residue_names = load_lipid_residue_names(options.lipid_fragments_json)
    protein_atoms, lipid_residues, _other_atoms = classify_lipid_atoms(
        structure.atoms, lipid_residue_names
    )
    if not protein_atoms:
        raise ValueError("no protein atoms found; cannot determine a cavity")

    ca_atoms = [atom for atom in protein_atoms if atom.name.strip().upper() == "CA"]
    axis_atoms = ca_atoms if ca_atoms else protein_atoms
    axis_coords = coords_from_atoms(axis_atoms)

    center, axis = principal_axis_pca(axis_coords)
    inner_radius, outer_radius, z_min, z_max = barrel_radii(axis_coords, center, axis)
    if options.cavity_radius is not None:
        inner_radius = options.cavity_radius

    removed_keys = {
        key
        for key, atoms in lipid_residues.items()
        if point_inside_cylinder(
            np.mean(coords_from_atoms(atoms), axis=0),
            center,
            axis,
            inner_radius,
            z_min,
            z_max,
            options.margin,
        )
    }

    kept_atoms = [
        atom
        for atom in structure.atoms
        if not (
            atom.resname.strip().upper() in lipid_residue_names
            and (atom.chain_id, atom.res_seq) in removed_keys
        )
    ]

    result = CavityLipidResult(
        total_lipid_residues=len(lipid_residues),
        removed_lipid_residues=len(removed_keys),
        kept_lipid_residues=len(lipid_residues) - len(removed_keys),
        center=coord_tuple(center),
        axis=coord_tuple(axis),
        inner_radius=float(inner_radius),
        outer_radius=float(outer_radius),
        z_min=float(z_min),
        z_max=float(z_max),
        margin=options.margin,
    )
    cleaned = structure.with_atoms(kept_atoms)
    return cleaned, result


def _ensure_single_model(structure: StructureData) -> None:
    model_ids = structure.model_ids()
    if len(model_ids) > 1:
        raise ValueError(
            "multi-model structures are not supported for cavity lipid removal; "
            "split the structure into a single model first"
        )


__all__ = ["remove_cavity_lipids"]

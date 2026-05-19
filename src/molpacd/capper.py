from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from molpacd._version import __version__
from molpacd.geometry import (
    analyze_openings,
    coords_from_atoms,
    filter_collisions,
    generate_cap_disk,
    select_atoms,
)
from molpacd.models import (
    AnalysisResult,
    AtomRecord,
    CapOptions,
    CapResult,
    CapSideResult,
    Opening,
    RemoveResult,
    StructureData,
    normalize_atom_name,
    normalize_chain_id,
    normalize_element,
    validate_cap_options,
)
from molpacd.names import choose_unused_residue_name, normalize_residue_name


@dataclass(frozen=True)
class _CapProvenance:
    count: int
    serial_start: Optional[int] = None
    serial_end: Optional[int] = None
    res_seq_start: Optional[int] = None
    res_seq_end: Optional[int] = None


def analyze_structure(
    structure: StructureData, options: Optional[CapOptions] = None
) -> AnalysisResult:
    options = options or CapOptions()
    validate_cap_options(options)
    _ensure_single_model(structure)
    selected = select_atoms(structure.atoms, options)
    axis, center, negative, positive = analyze_openings(selected, options)
    return AnalysisResult(
        axis=_tuple(axis),
        center=_tuple(center),
        negative=negative,
        positive=positive,
    )


def add_caps(
    structure: StructureData,
    options: Optional[CapOptions] = None,
) -> Tuple[StructureData, CapResult]:
    options = options or CapOptions()
    validate_cap_options(options)
    _ensure_single_model(structure)

    selected = select_atoms(structure.atoms, options)
    axis, center, negative, positive = analyze_openings(selected, options)
    analysis = AnalysisResult(
        axis=_tuple(axis),
        center=_tuple(center),
        negative=negative,
        positive=positive,
    )

    resname = (
        normalize_residue_name(options.resname)
        if options.resname
        else choose_unused_residue_name(structure.residue_names(), options.seed)
    )
    chain = normalize_chain_id(options.chain) or "Z"
    atom_name = normalize_atom_name(options.atom_name) or "O"
    element = normalize_element(options.element, atom_name)
    serial = max((atom.serial for atom in structure.atoms), default=0) + 1
    res_seq = _next_residue_number(structure.atoms, chain)
    existing_coords = coords_from_atoms(structure.atoms)

    cap_atoms: List[AtomRecord] = []
    side_results: List[CapSideResult] = []
    target_sides = _target_openings(options.sides, negative, positive)
    shared_radius = (
        max(opening.radius for opening in target_sides) if options.shared_radius else None
    )

    for opening in target_sides:
        side_sign = -1.0 if opening.side == "negative" else 1.0
        opening_radius = shared_radius if shared_radius is not None else opening.radius
        radius = opening_radius * options.radius_scale
        requested = generate_cap_disk(
            opening.centroid, _tuple(axis), radius, options.spacing, side_sign
        )
        accepted, skipped = filter_collisions(requested, existing_coords, options.min_clearance)

        for position in accepted:
            cap_atoms.append(
                AtomRecord(
                    record="HETATM",
                    serial=serial,
                    name=atom_name,
                    resname=resname,
                    chain_id=chain,
                    res_seq=res_seq,
                    x=float(position[0]),
                    y=float(position[1]),
                    z=float(position[2]),
                    occupancy=1.0,
                    bfactor=30.0,
                    element=element,
                )
            )
            serial += 1
            res_seq += 1

        if accepted.size:
            existing_coords = np.vstack([existing_coords, accepted])

        side_results.append(
            CapSideResult(
                side=opening.side,
                requested_count=int(len(requested)),
                added_count=int(len(accepted)),
                skipped_collision_count=skipped,
                centroid=opening.centroid,
                radius=float(radius),
            )
        )

    metadata = _cap_metadata(
        resname=resname,
        chain=chain,
        atom_name=atom_name,
        axis=axis,
        spacing=options.spacing,
        sides=options.sides,
        added_count=len(cap_atoms),
        serial_start=cap_atoms[0].serial if cap_atoms else None,
        serial_end=cap_atoms[-1].serial if cap_atoms else None,
        res_seq_start=cap_atoms[0].res_seq if cap_atoms else None,
        res_seq_end=cap_atoms[-1].res_seq if cap_atoms else None,
        side_results=side_results,
    )
    capped = structure.with_atoms([*structure.atoms, *cap_atoms], metadata=metadata)
    result = CapResult(
        resname=resname,
        chain=chain,
        atom_name=atom_name,
        added_count=len(cap_atoms),
        analysis=analysis,
        sides=side_results,
        serial_start=cap_atoms[0].serial if cap_atoms else None,
        serial_end=cap_atoms[-1].serial if cap_atoms else None,
        res_seq_start=cap_atoms[0].res_seq if cap_atoms else None,
        res_seq_end=cap_atoms[-1].res_seq if cap_atoms else None,
    )
    return capped, result


def remove_caps(
    structure: StructureData,
    resname: Optional[str] = None,
    chain: Optional[str] = None,
    atom_name: Optional[str] = None,
    force: bool = False,
) -> Tuple[StructureData, RemoveResult]:
    metadata_resname = structure.metadata.get("resname")
    metadata_present = bool(metadata_resname)

    if not metadata_present and not force:
        raise ValueError("cap removal without MolPACD metadata requires --force")

    metadata_effective_resname = (
        normalize_residue_name(metadata_resname) if metadata_resname else None
    )
    override_resname = normalize_residue_name(resname) if resname else None
    resname_mismatch = _metadata_override_mismatch(
        "resname",
        override_resname,
        metadata_effective_resname,
        override_resname is not None,
        metadata_present,
        force,
    )
    effective_resname = override_resname or metadata_effective_resname
    if not effective_resname:
        raise ValueError("cap removal requires a residue name when MolPACD metadata is absent")

    metadata_chain = normalize_chain_id(structure.metadata.get("chain"), default=None)
    override_chain = normalize_chain_id(chain, default=None) if chain is not None else None
    chain_mismatch = _metadata_override_mismatch(
        "chain",
        override_chain,
        metadata_chain,
        chain is not None,
        metadata_present,
        force,
    )
    effective_chain = override_chain if chain is not None else metadata_chain

    metadata_atom_name = normalize_atom_name(structure.metadata.get("atom_name"), default=None)
    override_atom_name = (
        normalize_atom_name(atom_name, default=None) if atom_name is not None else None
    )
    atom_name_mismatch = _metadata_override_mismatch(
        "atom name",
        override_atom_name,
        metadata_atom_name,
        atom_name is not None,
        metadata_present,
        force,
    )
    effective_atom_name = override_atom_name if atom_name is not None else metadata_atom_name

    use_metadata_provenance = not (
        force and (resname_mismatch or chain_mismatch or atom_name_mismatch)
    )
    provenance = (
        _cap_provenance_from_metadata(structure.metadata)
        if metadata_present and use_metadata_provenance
        else None
    )

    kept: List[AtomRecord] = []
    removed = 0
    for atom in structure.atoms:
        if _matches_cap(
            atom,
            effective_resname,
            effective_chain,
            effective_atom_name,
            provenance,
        ):
            removed += 1
        else:
            kept.append(atom)

    metadata: Dict[str, str] = {}
    decapped = structure.with_atoms(kept, metadata=metadata)
    return (
        decapped,
        RemoveResult(
            removed_count=removed,
            resname=effective_resname,
            chain=effective_chain,
            atom_name=effective_atom_name,
        ),
    )


def _matches_cap(
    atom: AtomRecord,
    resname: str,
    chain: Optional[str],
    atom_name: Optional[str],
    provenance: Optional[_CapProvenance] = None,
) -> bool:
    if atom.resname.upper() != resname:
        return False
    if chain is not None and atom.chain_id != chain:
        return False
    if atom_name is not None and atom.name.strip().upper() != atom_name:
        return False
    if provenance is None:
        return True
    if provenance.count == 0:
        return False
    if provenance.serial_start is not None and atom.serial < provenance.serial_start:
        return False
    if provenance.serial_end is not None and atom.serial > provenance.serial_end:
        return False
    if provenance.res_seq_start is not None and atom.res_seq < provenance.res_seq_start:
        return False
    return not (provenance.res_seq_end is not None and atom.res_seq > provenance.res_seq_end)


def _next_residue_number(atoms: Sequence[AtomRecord], chain: str) -> int:
    chain_numbers = [atom.res_seq for atom in atoms if atom.chain_id == chain]
    return max(chain_numbers, default=0) + 1


def _target_openings(sides: str, negative: Opening, positive: Opening) -> List[Opening]:
    if sides == "both":
        return [negative, positive]
    if sides == "negative":
        return [negative]
    if sides == "positive":
        return [positive]
    raise ValueError("sides must be one of both, negative, or positive")


def _cap_metadata(
    resname: str,
    chain: str,
    atom_name: str,
    axis: np.ndarray,
    spacing: float,
    sides: str,
    added_count: int,
    serial_start: Optional[int],
    serial_end: Optional[int],
    res_seq_start: Optional[int],
    res_seq_end: Optional[int],
    side_results: Sequence[CapSideResult],
) -> Dict[str, str]:
    metadata = {
        "version": __version__,
        "resname": resname,
        "chain": chain,
        "atom_name": atom_name,
        "axis": " ".join(f"{value:.8f}" for value in axis),
        "spacing": f"{spacing:.3f}",
        "sides": sides,
        "count": str(added_count),
        "side_counts": ",".join(f"{side.side}:{side.added_count}" for side in side_results),
    }
    if serial_start is not None:
        metadata["serial_start"] = str(serial_start)
    if serial_end is not None:
        metadata["serial_end"] = str(serial_end)
    if res_seq_start is not None:
        metadata["res_seq_start"] = str(res_seq_start)
    if res_seq_end is not None:
        metadata["res_seq_end"] = str(res_seq_end)
    return metadata


def _tuple(values: np.ndarray) -> Tuple[float, float, float]:
    return (float(values[0]), float(values[1]), float(values[2]))


def _ensure_single_model(structure: StructureData) -> None:
    model_ids = structure.model_ids()
    if len(model_ids) > 1:
        raise ValueError(
            "multi-model structures are not supported for analysis or cap addition; "
            "split the structure into a single model first"
        )


def _metadata_override_mismatch(
    label: str,
    override: Optional[str],
    metadata_value: Optional[str],
    override_provided: bool,
    metadata_present: bool,
    force: bool,
) -> bool:
    mismatch = override_provided and metadata_value is not None and override != metadata_value
    if metadata_present and mismatch and not force:
        raise ValueError(f"overriding MolPACD metadata {label} requires --force")
    return mismatch


def _cap_provenance_from_metadata(metadata: Dict[str, str]) -> Optional[_CapProvenance]:
    if "count" not in metadata:
        return None

    count = _metadata_int(metadata, "count")
    if count < 0:
        raise ValueError("invalid MolPACD metadata: count must be zero or greater")
    if count == 0:
        return _CapProvenance(count=0)

    provenance = _CapProvenance(
        count=count,
        serial_start=_metadata_int(metadata, "serial_start"),
        serial_end=_metadata_int(metadata, "serial_end"),
        res_seq_start=_metadata_int(metadata, "res_seq_start"),
        res_seq_end=_metadata_int(metadata, "res_seq_end"),
    )
    if (
        provenance.serial_start is not None
        and provenance.serial_end is not None
        and provenance.serial_start > provenance.serial_end
    ):
        raise ValueError("invalid MolPACD metadata: serial range is reversed")
    if (
        provenance.res_seq_start is not None
        and provenance.res_seq_end is not None
        and provenance.res_seq_start > provenance.res_seq_end
    ):
        raise ValueError("invalid MolPACD metadata: residue sequence range is reversed")
    return provenance


def _metadata_int(metadata: Dict[str, str], key: str) -> int:
    value = metadata.get(key)
    if value is None:
        raise ValueError(f"invalid MolPACD metadata: {key} is missing")
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"invalid MolPACD metadata: {key} must be an integer") from exc

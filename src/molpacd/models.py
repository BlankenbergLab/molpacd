from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Literal, Optional, Tuple

import numpy as np

Selection = Literal["ca", "backbone", "all"]
Sides = Literal["both", "positive", "negative"]


@dataclass(frozen=True)
class AtomRecord:
    record: str
    serial: int
    name: str
    resname: str
    chain_id: str
    res_seq: int
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    bfactor: float = 0.0
    element: str = ""
    altloc: str = ""
    insertion_code: str = ""
    model: int = 1

    @property
    def coord(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)


@dataclass(frozen=True)
class StructureData:
    atoms: List[AtomRecord]
    source_format: Literal["pdb", "cif"]
    header_lines: List[str] = field(default_factory=list)
    metadata: Dict[str, str] = field(default_factory=dict)

    def residue_names(self) -> set[str]:
        return {atom.resname.upper() for atom in self.atoms}

    def model_ids(self) -> set[int]:
        return {atom.model for atom in self.atoms}

    def with_atoms(
        self,
        atoms: List[AtomRecord],
        metadata: Optional[Dict[str, str]] = None,
    ) -> "StructureData":
        return StructureData(
            atoms=atoms,
            source_format=self.source_format,
            header_lines=list(self.header_lines),
            metadata=dict(self.metadata if metadata is None else metadata),
        )


@dataclass(frozen=True)
class CapOptions:
    axis: str = "auto"
    selection: Selection = "ca"
    sides: Sides = "both"
    spacing: float = 2.8
    radius_scale: float = 1.0
    window: float = 6.0
    min_atoms: int = 12
    chain: str = "Z"
    resname: Optional[str] = None
    atom_name: str = "O"
    element: str = "O"
    seed: Optional[int] = None
    min_clearance: float = 1.4
    shared_radius: bool = True


@dataclass(frozen=True)
class Opening:
    side: Literal["negative", "positive"]
    atom_count: int
    centroid: Tuple[float, float, float]
    radius: float
    projection: float


@dataclass(frozen=True)
class AnalysisResult:
    axis: Tuple[float, float, float]
    center: Tuple[float, float, float]
    negative: Opening
    positive: Opening


@dataclass(frozen=True)
class CapSideResult:
    side: Literal["negative", "positive"]
    requested_count: int
    added_count: int
    skipped_collision_count: int
    centroid: Tuple[float, float, float]
    radius: float


@dataclass(frozen=True)
class CapResult:
    resname: str
    chain: str
    atom_name: str
    added_count: int
    analysis: AnalysisResult
    sides: List[CapSideResult]
    serial_start: Optional[int] = None
    serial_end: Optional[int] = None
    res_seq_start: Optional[int] = None
    res_seq_end: Optional[int] = None


@dataclass(frozen=True)
class RemoveResult:
    removed_count: int
    resname: str
    chain: Optional[str]
    atom_name: Optional[str]


def validate_cap_options(options: CapOptions) -> None:
    if options.selection not in {"ca", "backbone", "all"}:
        raise ValueError("selection must be one of ca, backbone, or all")
    if options.sides not in {"both", "negative", "positive"}:
        raise ValueError("sides must be one of both, negative, or positive")
    if options.spacing <= 0:
        raise ValueError("spacing must be greater than zero")
    if options.radius_scale <= 0:
        raise ValueError("radius scale must be greater than zero")
    if options.window <= 0:
        raise ValueError("window must be greater than zero")
    if options.min_atoms < 3:
        raise ValueError("min_atoms must be at least 3")
    if options.min_clearance < 0:
        raise ValueError("min_clearance must be zero or greater")
    if not isinstance(options.shared_radius, bool):
        raise ValueError("shared_radius must be a boolean")
    normalize_chain_id(options.chain)
    atom_name = normalize_atom_name(options.atom_name)
    normalize_element(options.element, atom_name or "")


def normalize_chain_id(chain: Optional[str], default: Optional[str] = "Z") -> Optional[str]:
    cleaned = (chain or "").strip()
    if not cleaned:
        return default
    if len(cleaned) != 1:
        raise ValueError("chain identifier must be exactly one non-whitespace character")
    return cleaned


def normalize_atom_name(atom_name: Optional[str], default: Optional[str] = "O") -> Optional[str]:
    cleaned = (atom_name or "").strip().upper()
    if not cleaned:
        return default
    if len(cleaned) > 4 or any(char.isspace() for char in cleaned):
        raise ValueError("atom name must be 1 to 4 non-whitespace characters")
    return cleaned


def normalize_element(element: str, atom_name: str) -> str:
    cleaned = element.strip().upper()
    if not cleaned:
        cleaned = atom_name.strip().upper()[:1]
    if not 1 <= len(cleaned) <= 2 or not cleaned.isalpha():
        raise ValueError("element must be a 1 or 2 letter chemical symbol")
    return cleaned

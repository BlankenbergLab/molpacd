from __future__ import annotations

from io import StringIO
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Dict, Iterable, List, Optional, Sequence

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser

from molpacd.models import AtomRecord, StructureData


def read_structure(path: Path) -> StructureData:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".pdb", ".ent"}:
        return _read_pdb(path)
    if suffix in {".cif", ".mmcif"}:
        return _read_cif(path)
    raise ValueError(f"unsupported structure format for {path}; expected PDB or mmCIF")


def write_structure(structure: StructureData, path: Path) -> None:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".pdb", ".ent"}:
        _write_pdb(structure, path)
        return
    if suffix in {".cif", ".mmcif"}:
        _write_cif(structure, path)
        return
    raise ValueError(f"unsupported output format for {path}; expected PDB or mmCIF")


def _read_pdb(path: Path) -> StructureData:
    header_lines, metadata = _scan_pdb_metadata(path)
    parser: Any = PDBParser(QUIET=True)  # type: ignore[no-untyped-call]
    structure = parser.get_structure("molpacd", StringIO(_pdb_coordinate_text(path)))
    atoms = _records_from_biopython(structure)
    return StructureData(
        atoms=atoms,
        source_format="pdb",
        header_lines=header_lines,
        metadata=metadata,
    )


def _read_cif(path: Path) -> StructureData:
    metadata = _scan_cif_metadata(path)
    parser: Any = MMCIFParser(QUIET=True)  # type: ignore[no-untyped-call]
    structure = parser.get_structure("molpacd", str(path))
    atoms = _records_from_biopython(structure)
    return StructureData(atoms=atoms, source_format="cif", metadata=metadata)


def _pdb_coordinate_text(path: Path) -> str:
    prefixes = ("ATOM", "HETATM", "MODEL", "ENDMDL", "TER", "END")
    return "".join(
        line for line in path.read_text().splitlines(keepends=True) if line.startswith(prefixes)
    )


def _records_from_biopython(structure: Any) -> List[AtomRecord]:
    records: List[AtomRecord] = []
    fallback_serial = 1
    for model in structure:
        model_id = int(model.id) + 1 if isinstance(model.id, int) else 1
        for chain in model:
            chain_id = str(chain.id).strip() or " "
            for residue in chain:
                hetflag, resseq, insertion_code = residue.id
                record = "ATOM" if str(hetflag).strip() == "" else "HETATM"
                resname = residue.get_resname().strip()
                for atom in residue:
                    serial = int(atom.get_serial_number() or fallback_serial)
                    fallback_serial = max(fallback_serial + 1, serial + 1)
                    coord = atom.get_coord()
                    occupancy = atom.get_occupancy()
                    bfactor = atom.get_bfactor()
                    element = str(getattr(atom, "element", "") or "").strip()
                    records.append(
                        AtomRecord(
                            record=record,
                            serial=serial,
                            name=atom.get_name().strip(),
                            resname=resname,
                            chain_id=chain_id,
                            res_seq=int(resseq),
                            insertion_code=str(insertion_code).strip(),
                            x=float(coord[0]),
                            y=float(coord[1]),
                            z=float(coord[2]),
                            occupancy=float(occupancy if occupancy is not None else 1.0),
                            bfactor=float(bfactor if bfactor is not None else 0.0),
                            element=element,
                            altloc=str(atom.get_altloc()).strip(),
                            model=model_id,
                        )
                    )
    return records


def _scan_pdb_metadata(path: Path) -> tuple[List[str], Dict[str, str]]:
    header_lines: List[str] = []
    metadata: Dict[str, str] = {}
    in_atoms = False
    for line in path.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM", "MODEL")):
            in_atoms = True
        if line.startswith("REMARK MOLPACD"):
            _parse_metadata_tokens(line.split()[2:], metadata)
            continue
        if not in_atoms:
            header_lines.append(line)
    return header_lines, metadata


def _scan_cif_metadata(path: Path) -> Dict[str, str]:
    metadata: Dict[str, str] = {}
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("# MOLPACD"):
            _parse_metadata_tokens(stripped.split()[2:], metadata)
    return metadata


def _parse_metadata_tokens(tokens: Sequence[str], metadata: Dict[str, str]) -> None:
    if not tokens:
        return
    key = tokens[0].lower()
    value = " ".join(tokens[1:]).strip()
    if key:
        metadata[key] = value


def _write_pdb(structure: StructureData, path: Path) -> None:
    lines: List[str] = []
    lines.extend(structure.header_lines)
    lines.extend(_pdb_metadata_lines(structure.metadata))
    current_model = None
    has_multiple_models = len({atom.model for atom in structure.atoms}) > 1
    for atom in structure.atoms:
        if has_multiple_models and atom.model != current_model:
            if current_model is not None:
                lines.append("ENDMDL\n")
            current_model = atom.model
            lines.append(f"MODEL     {atom.model:>4d}\n")
        lines.append(_format_pdb_atom(atom))
    if has_multiple_models:
        lines.append("ENDMDL\n")
    lines.append("END\n")
    _write_text_atomic(path, "".join(lines))


def _pdb_metadata_lines(metadata: Dict[str, str]) -> List[str]:
    return [f"REMARK MOLPACD {key.upper()} {value}\n" for key, value in sorted(metadata.items())]


def _format_pdb_atom(atom: AtomRecord) -> str:
    atom_name = atom.name[:4]
    altloc = (atom.altloc or " ")[:1]
    insertion_code = (atom.insertion_code or " ")[:1]
    chain = (atom.chain_id or " ")[:1]
    element = (atom.element or atom.name[:1]).strip().upper()[:2]
    return (
        f"{atom.record:<6}"
        f"{atom.serial:>5d} "
        f"{atom_name:<4}"
        f"{altloc}"
        f"{atom.resname[:3]:>3} "
        f"{chain}"
        f"{atom.res_seq:>4d}"
        f"{insertion_code}   "
        f"{atom.x:>8.3f}"
        f"{atom.y:>8.3f}"
        f"{atom.z:>8.3f}"
        f"{atom.occupancy:>6.2f}"
        f"{atom.bfactor:>6.2f}"
        f"          "
        f"{element:>2}\n"
    )


def _write_cif(structure: StructureData, path: Path) -> None:
    lines = [
        "data_molpacd\n",
        "#\n",
    ]
    for key, value in sorted(structure.metadata.items()):
        lines.append(f"# MOLPACD {key.upper()} {value}\n")
    lines.extend(
        [
            "loop_\n",
            "_atom_site.group_PDB\n",
            "_atom_site.id\n",
            "_atom_site.type_symbol\n",
            "_atom_site.label_atom_id\n",
            "_atom_site.label_alt_id\n",
            "_atom_site.label_comp_id\n",
            "_atom_site.label_asym_id\n",
            "_atom_site.label_entity_id\n",
            "_atom_site.label_seq_id\n",
            "_atom_site.pdbx_PDB_ins_code\n",
            "_atom_site.Cartn_x\n",
            "_atom_site.Cartn_y\n",
            "_atom_site.Cartn_z\n",
            "_atom_site.occupancy\n",
            "_atom_site.B_iso_or_equiv\n",
            "_atom_site.auth_seq_id\n",
            "_atom_site.auth_comp_id\n",
            "_atom_site.auth_asym_id\n",
            "_atom_site.auth_atom_id\n",
            "_atom_site.pdbx_PDB_model_num\n",
        ]
    )
    for index, atom in enumerate(structure.atoms, start=1):
        element = (atom.element or atom.name[:1]).strip().upper()
        altloc = atom.altloc or "."
        insertion_code = atom.insertion_code or "?"
        values = [
            atom.record,
            str(index),
            element,
            atom.name,
            altloc,
            atom.resname,
            atom.chain_id or ".",
            "1",
            str(atom.res_seq),
            insertion_code,
            f"{atom.x:.3f}",
            f"{atom.y:.3f}",
            f"{atom.z:.3f}",
            f"{atom.occupancy:.2f}",
            f"{atom.bfactor:.2f}",
            str(atom.res_seq),
            atom.resname,
            atom.chain_id or ".",
            atom.name,
            str(atom.model),
        ]
        lines.append(" ".join(_cif_value(value) for value in values) + "\n")
    lines.append("#\n")
    _write_text_atomic(path, "".join(lines))


def _cif_value(value: str) -> str:
    if value == "":
        return "."
    if any(char.isspace() for char in value) or value[0] in {"_", "#", ";"} or "'" in value:
        escaped = value.replace("'", "''")
        return f"'{escaped}'"
    return value


def atom_coordinates(atoms: Iterable[AtomRecord]) -> List[tuple[float, float, float]]:
    return [(atom.x, atom.y, atom.z) for atom in atoms]


def _write_text_atomic(path: Path, text: str) -> None:
    path = Path(path)
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

from __future__ import annotations

import json
from collections import defaultdict
from importlib import resources
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple, Union

from molpacd.models import AtomRecord

# Standard amino acid residue names, including common AMBER/CHARMM
# protonation-state variants, used to identify protein atoms.
AMINO_ACID_RESIDUE_NAMES: Set[str] = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "HID",
    "HIE",
    "HIP",
    "HSD",
    "HSE",
    "HSP",  # histidine protonation variants
    "CYX",
    "CYM",
    "ASH",
    "GLH",
    "LYN",  # other protonation/disulfide variants
}

# Packmol-memgen (lipid21 / AMBER lipid builder) writes each lipid into the
# structure as separate head/tail *fragment* residues (e.g. POPC becomes
# "OL", "PA", "PC" residues), not as one residue named after the full lipid.
# These fallback sets mirror data/lipid_fragments.json and are merged with it
# so lipid detection keeps working even if that file cannot be read.
_FALLBACK_LIPID_FRAGMENT_CODES: Set[str] = {
    "AR",
    "CHL",
    "DHA",
    "LAL",
    "MY",
    "OL",
    "PA",
    "PC",
    "PE",
    "PGR",
    "PH-",
    "PS",
    "SA",
    "SPM",
    "ST",
}

_FALLBACK_FULL_LIPID_NAMES: Set[str] = {
    "AHPA",
    "AHPC",
    "AHPE",
    "AHPG",
    "AHPS",
    "ALPA",
    "ALPC",
    "ALPE",
    "ALPG",
    "ALPS",
    "AMPA",
    "AMPC",
    "AMPE",
    "AMPG",
    "AMPS",
    "AOPA",
    "AOPC",
    "AOPE",
    "AOPG",
    "AOPS",
    "APPA",
    "APPC",
    "APPE",
    "APPG",
    "APPS",
    "ASM",
    "ASPA",
    "ASPC",
    "ASPE",
    "ASPG",
    "ASPS",
    "CHL1",
    "DAPA",
    "DAPC",
    "DAPE",
    "DAPG",
    "DAPS",
    "DHPA",
    "DHPC",
    "DHPE",
    "DHPG",
    "DHPS",
    "DLPA",
    "DLPC",
    "DLPE",
    "DLPG",
    "DLPS",
    "DMPA",
    "DMPC",
    "DMPE",
    "DMPG",
    "DMPS",
    "DOPA",
    "DOPC",
    "DOPE",
    "DOPG",
    "DOPS",
    "DPPA",
    "DPPC",
    "DPPE",
    "DPPG",
    "DPPS",
    "DSPA",
    "DSPC",
    "DSPE",
    "DSPG",
    "DSPS",
    "HAPA",
    "HAPC",
    "HAPE",
    "HAPG",
    "HAPS",
    "HLPA",
    "HLPC",
    "HLPE",
    "HLPG",
    "HLPS",
    "HMPA",
    "HMPC",
    "HMPE",
    "HMPG",
    "HMPS",
    "HOPA",
    "HOPC",
    "HOPE",
    "HOPG",
    "HOPS",
    "HPPA",
    "HPPC",
    "HPPE",
    "HPPG",
    "HPPS",
    "HSM",
    "HSPA",
    "HSPC",
    "HSPE",
    "HSPG",
    "HSPS",
    "LAPA",
    "LAPC",
    "LAPE",
    "LAPG",
    "LAPS",
    "LHPA",
    "LHPC",
    "LHPE",
    "LHPG",
    "LHPS",
    "LMPA",
    "LMPC",
    "LMPE",
    "LMPG",
    "LMPS",
    "LOPA",
    "LOPC",
    "LOPE",
    "LOPG",
    "LOPS",
    "LPPA",
    "LPPC",
    "LPPE",
    "LPPG",
    "LPPS",
    "LSM",
    "LSPA",
    "LSPC",
    "LSPE",
    "LSPG",
    "LSPS",
    "MAPA",
    "MAPC",
    "MAPE",
    "MAPG",
    "MAPS",
    "MHPA",
    "MHPC",
    "MHPE",
    "MHPG",
    "MHPS",
    "MLPA",
    "MLPC",
    "MLPE",
    "MLPG",
    "MLPS",
    "MOPA",
    "MOPC",
    "MOPE",
    "MOPG",
    "MOPS",
    "MPPA",
    "MPPC",
    "MPPE",
    "MPPG",
    "MPPS",
    "MSM",
    "MSPA",
    "MSPC",
    "MSPE",
    "MSPG",
    "MSPS",
    "OAPA",
    "OAPC",
    "OAPE",
    "OAPG",
    "OAPS",
    "OHPA",
    "OHPC",
    "OHPE",
    "OHPG",
    "OHPS",
    "OLPA",
    "OLPC",
    "OLPE",
    "OLPG",
    "OLPS",
    "OMPA",
    "OMPC",
    "OMPE",
    "OMPG",
    "OMPS",
    "OPPA",
    "OPPC",
    "OPPE",
    "OPPG",
    "OPPS",
    "OSM",
    "OSPA",
    "OSPC",
    "OSPE",
    "OSPG",
    "OSPS",
    "PAPA",
    "PAPC",
    "PAPE",
    "PAPG",
    "PAPS",
    "PHPA",
    "PHPC",
    "PHPE",
    "PHPG",
    "PHPS",
    "PLPA",
    "PLPC",
    "PLPE",
    "PLPG",
    "PLPS",
    "PMPA",
    "PMPC",
    "PMPE",
    "PMPG",
    "PMPS",
    "POPA",
    "POPC",
    "POPE",
    "POPG",
    "POPS",
    "PSM",
    "PSPA",
    "PSPC",
    "PSPE",
    "PSPG",
    "PSPS",
    "SAPA",
    "SAPC",
    "SAPE",
    "SAPG",
    "SAPS",
    "SDPA",
    "SDPC",
    "SDPE",
    "SDPG",
    "SDPS",
    "SHPA",
    "SHPC",
    "SHPE",
    "SHPG",
    "SHPS",
    "SLPA",
    "SLPC",
    "SLPE",
    "SLPG",
    "SLPS",
    "SMPA",
    "SMPC",
    "SMPE",
    "SMPG",
    "SMPS",
    "SOPA",
    "SOPC",
    "SOPE",
    "SOPG",
    "SOPS",
    "SPPA",
    "SPPC",
    "SPPE",
    "SPPG",
    "SPPS",
    "SSM",
}

# A few common alternate/legacy names not present in data/lipid_fragments.json
# that still show up in some packmol-memgen / CHARMM-GUI style structures.
_EXTRA_LIPID_RESIDUE_NAMES: Set[str] = {"PG", "PI", "POPI", "DOPI", "DPPI", "CHOL", "LA"}

_FALLBACK_LIPID_RESIDUE_NAMES: Set[str] = (
    _FALLBACK_LIPID_FRAGMENT_CODES | _FALLBACK_FULL_LIPID_NAMES | _EXTRA_LIPID_RESIDUE_NAMES
)


def load_lipid_fragments(
    lipid_fragments_json: Optional[Union[str, Path]] = None,
) -> Dict[str, List[str]]:
    """Load the lipid-name to head/tail-fragment mapping.

    With no path given, loads the copy of ``lipid_fragments.json`` bundled
    with the package. Pass an explicit path to load a different or updated
    mapping instead (for example, a newer copy from a MemGen checkout).
    """
    if lipid_fragments_json is not None:
        with Path(lipid_fragments_json).open(encoding="utf-8") as handle:
            data = json.load(handle)
    else:
        with (
            resources.files("molpacd")
            .joinpath("data/lipid_fragments.json")
            .open("r", encoding="utf-8") as handle
        ):
            data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("lipid_fragments.json must contain a JSON object")
    return data


def load_lipid_residue_names(lipid_fragments_json: Optional[Union[str, Path]] = None) -> Set[str]:
    """Build the set of residue names that should be treated as lipid atoms.

    Combines every full lipid name and every head/tail fragment code from
    ``load_lipid_fragments`` with a hardcoded fallback set mirroring the
    bundled data, so lipid detection keeps working even in the unlikely case
    the bundled or overriding JSON file cannot be fully read.
    """
    data = load_lipid_fragments(lipid_fragments_json)
    names = set(_FALLBACK_LIPID_RESIDUE_NAMES)
    names.update(name.upper() for name in data)
    for fragments in data.values():
        names.update(str(fragment).upper() for fragment in fragments)
    return names


def classify_lipid_atoms(
    atoms: Sequence[AtomRecord], lipid_residue_names: Set[str]
) -> Tuple[List[AtomRecord], Dict[Tuple[str, int], List[AtomRecord]], List[AtomRecord]]:
    """Split atoms into protein, lipid residues, and everything else.

    Lipid residues are grouped by ``(chain_id, res_seq)`` since a single
    lipid is typically stored as several same-numbered head/tail fragment
    residues.
    """
    protein_atoms: List[AtomRecord] = []
    lipid_residues: Dict[Tuple[str, int], List[AtomRecord]] = defaultdict(list)
    other_atoms: List[AtomRecord] = []

    for atom in atoms:
        resname = atom.resname.strip().upper()
        if resname in AMINO_ACID_RESIDUE_NAMES:
            protein_atoms.append(atom)
        elif resname in lipid_residue_names:
            lipid_residues[(atom.chain_id, atom.res_seq)].append(atom)
        else:
            other_atoms.append(atom)

    return protein_atoms, dict(lipid_residues), other_atoms

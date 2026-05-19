from __future__ import annotations

import random
import string
from typing import Iterable, Optional

RESERVED_RESIDUE_NAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
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
    "A",
    "C",
    "G",
    "U",
    "DA",
    "DC",
    "DG",
    "DT",
    "HOH",
    "WAT",
    "TIP",
    "SOL",
    "NA",
    "CL",
    "K",
    "CA",
    "MG",
    "ZN",
    "MN",
    "FE",
    "CU",
    "CO",
    "NI",
}


def normalize_residue_name(resname: str) -> str:
    cleaned = resname.strip().upper()
    if not cleaned.isalnum() or not 1 <= len(cleaned) <= 3:
        raise ValueError("cap residue name must be 1 to 3 alphanumeric characters")
    return cleaned


def choose_unused_residue_name(existing: Iterable[str], seed: Optional[int] = None) -> str:
    used = {name.strip().upper() for name in existing}
    blocked = used | RESERVED_RESIDUE_NAMES
    rng = random.Random(seed)
    alphabet = string.ascii_uppercase + string.digits

    for _ in range(10000):
        candidate = "".join(rng.choice(alphabet) for _ in range(3))
        if candidate not in blocked and not candidate[0].isdigit():
            return candidate

    for first in string.ascii_uppercase:
        for second in alphabet:
            for third in alphabet:
                candidate = f"{first}{second}{third}"
                if candidate not in blocked:
                    return candidate

    raise RuntimeError("could not find an unused 3-character residue name")

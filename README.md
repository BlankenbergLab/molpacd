# MolPACD

**MolPACD** is the **Molecular Pore Aperture Cap Designer**. It designs, adds,
analyzes, and removes molecular aperture caps in protein structures.

The first implementation is derived from the MemGen beta-barrel water-cap script,
but generalizes the workflow for PDB and mmCIF structures with configurable axes,
cap identity, and cap removal metadata.

## Install

```bash
python -m pip install .
```

For development:

```bash
python -m pip install -e ".[dev]"
```

With Conda, keep the environment local to the checkout:

```bash
conda env create -p .conda/molpacd-dev -f environment-dev.yml
conda run -p .conda/molpacd-dev python -m nox -s lint format mypy tests-3.13 build
```

Common development tasks are also available through `make`:

```bash
make env
make check
```

## CLI

Analyze candidate apertures:

```bash
molpacd analyze input.pdb --json
```

Add caps to both apertures:

```bash
molpacd add input.pdb -o capped.pdb
```

Use a membrane-oriented z-axis and a fixed cap residue name:

```bash
molpacd add input.pdb -o capped.pdb --axis z --resname DUM
```

Use each aperture's inferred radius instead of sharing the larger radius:

```bash
molpacd add input.pdb -o capped.pdb --axis z --independent-radius
```

Remove MolPACD-generated caps:

```bash
molpacd remove capped.pdb -o decapped.pdb
```

Remove matching caps from a file without MolPACD metadata:

```bash
molpacd remove capped.pdb -o decapped.pdb --resname DUM --chain Z --force
```

MolPACD writes cap provenance metadata and uses it during removal so matching
non-cap atoms are not removed accidentally. If metadata is absent, or if you
override metadata values such as residue name, chain, or atom name, removal
requires `--force`.

## Format Notes

MolPACD reads PDB and mmCIF files and writes minimal PDB/mmCIF outputs focused
on atom records plus MolPACD metadata. It preserves PDB header lines before the
coordinate section, but it does not attempt to reproduce every original record
or mmCIF category. Multi-model structures are not supported for analysis or cap
addition; split those structures into single-model inputs before running
MolPACD.

## Library

```python
from pathlib import Path

from molpacd import CapOptions, add_caps, read_structure, write_structure

structure = read_structure(Path("input.pdb"))
capped, result = add_caps(structure, CapOptions(axis="z", resname="DUM"))
write_structure(capped, Path("capped.pdb"))

print(result.resname, result.added_count)
```

## Development Checks

```bash
python -m nox
```

The nox sessions run tests, linting, formatting checks, type checks, and package
build validation. GitHub Actions runs the same checks on pushes and pull requests.

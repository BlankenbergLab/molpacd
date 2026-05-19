# Quick Start

## Install

```bash
python -m pip install .
```

For a development checkout:

```bash
python -m pip install -e ".[dev]"
```

## Analyze A Structure

```bash
molpacd analyze input.pdb --axis z --json
```

The JSON output reports the inferred analysis axis, structure center, and the
negative and positive aperture summaries.

## Add Caps

```bash
molpacd add input.pdb -o capped.pdb --axis z
```

Use a fixed residue name when downstream tools expect a known cap identity:

```bash
molpacd add input.pdb -o capped.pdb --axis z --resname DUM
```

Use each aperture's inferred radius instead of sharing the larger radius:

```bash
molpacd add input.pdb -o capped.pdb --axis z --independent-radius
```

Preview the design without writing an output file:

```bash
molpacd add input.pdb -o capped.pdb --axis z --dry-run --json
```

## Remove Caps

```bash
molpacd remove capped.pdb -o decapped.pdb --json
```

MolPACD writes cap provenance metadata and uses it during removal so matching
non-cap atoms are not removed accidentally. If metadata is absent, or if you
override metadata values such as residue name, chain, or atom name, removal
requires `--force`.

```bash
molpacd remove capped.pdb -o decapped.pdb --resname DUM --chain Z --force
```

## Python API

```python
from pathlib import Path

from molpacd import CapOptions, add_caps, read_structure, write_structure

structure = read_structure(Path("input.pdb"))
capped, result = add_caps(structure, CapOptions(axis="z", resname="DUM"))
write_structure(capped, Path("capped.pdb"))

print(result.resname, result.added_count)
```

## Format Notes

Multi-model structures are not supported for analysis or cap addition. Split
those structures into single-model inputs before running MolPACD.

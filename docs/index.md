# MolPACD

**MolPACD** is the **Molecular Pore Aperture Cap Designer**. It provides a
Python API and command line interface for analyzing pore openings, adding
molecular aperture caps, and removing MolPACD-generated caps from protein
structure files.

MolPACD works with PDB and mmCIF inputs and writes minimal PDB/mmCIF outputs
focused on atom records plus MolPACD provenance metadata.

## Installation

For a minimal install from a local checkout:

```bash
python -m pip install .
```

For development:

```bash
python -m pip install -e ".[dev]"
```

For documentation work:

```bash
make docs-deps
make docs-serve
```

## Start Here

- Follow the [Quick Start](quickstart.md) for common CLI and Python examples.
- Use the [Command Line](cli.md) reference for `analyze`, `add`, and `remove`.
- Browse the [Basic API](api.md) for the import-root interface.
- Use the [Advanced API](advanced-api.md) for module-level details.

## Blankenberg Lab

MolPACD is developed by the
[Blankenberg Lab](https://www.blankenberglab.org/). Related tools are available
through [tools.blankenberglab.org](https://tools.blankenberglab.org/).

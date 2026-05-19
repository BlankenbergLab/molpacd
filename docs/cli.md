# Command Line

MolPACD installs the `molpacd` command.

```bash
molpacd --help
molpacd --version
```

Use `--debug` before the subcommand when you need a traceback for unexpected
errors:

```bash
molpacd --debug analyze input.pdb
```

## Analyze

Analyze candidate apertures without writing a structure file.

```bash
molpacd analyze input.pdb --axis auto
molpacd analyze input.pdb --axis z --selection backbone --json
```

Shared analysis options:

- `--axis`: `auto`, `x`, `y`, `z`, or a 3-value vector such as `0,0,1`.
- `--selection`: `ca`, `backbone`, or `all`.
- `--sides`: `both`, `negative`, or `positive`.
- `--spacing`: cap lattice spacing in Angstrom.
- `--radius-scale`: scale inferred opening radii.
- `--window`: opening slice width in Angstrom.
- `--min-atoms`: minimum atoms per opening slice.

## Add

Add molecular aperture caps and write a capped structure.

```bash
molpacd add input.pdb -o capped.pdb --axis z
molpacd add input.cif -o capped.cif --axis 0,0,1 --resname DUM
```

Additional options:

- `--chain`: one-character chain identifier for cap atoms.
- `--resname`: 1-3 character cap residue name.
- `--atom-name`: cap atom name.
- `--element`: cap atom element.
- `--seed`: seed for automatic residue-name generation.
- `--min-clearance`: minimum distance from existing atoms before skipping a cap atom.
- `--shared-radius`: use the larger opening radius for all requested sides.
- `--independent-radius`: use each opening's own inferred radius.
- `--dry-run`: report the cap design without writing output.
- `--json`: write machine-readable JSON.

## Remove

Remove MolPACD-generated caps.

```bash
molpacd remove capped.pdb -o decapped.pdb
molpacd remove capped.cif -o decapped.cif --json
```

When metadata is absent, or when metadata values are intentionally overridden,
pass `--force` with explicit matching criteria:

```bash
molpacd remove capped.pdb -o decapped.pdb --resname DUM --chain Z --force
```

Additional options:

- `--resname`: cap residue name to remove.
- `--chain`: cap chain identifier to remove.
- `--atom-name`: cap atom name to remove.
- `--force`: allow removal by residue/chain/atom match when metadata is absent
  or overridden.
- `--json`: write machine-readable JSON.

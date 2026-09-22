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

## Remove Cavity Lipids

Remove lipids that packing tools (such as packmol-memgen) placed inside a
protein cavity, e.g. the interior of a beta barrel, instead of only in the
surrounding bilayer.

```bash
molpacd remove-lipids packed.pdb -o cleaned.pdb
molpacd remove-lipids packed.pdb -o cleaned.pdb --cavity-radius 8.0 --margin 2.0 --json
```

The cavity is estimated from the protein's CA atoms (falling back to all
protein atoms if none are found): a principal axis and center from PCA, and
an inner radius from the average backbone radius. Lipid residues (grouped by
chain and residue number, since a lipid is often stored as separate head/tail
fragment residues) whose center of mass falls inside that cavity, plus
`--margin`, are removed.

Lipid residue names are recognized from a bundled `lipid_fragments.json`
mapping full lipid names (e.g. `POPC`) to their head/tail fragment residue
names (e.g. `OL`, `PA`, `PC`). Additional options:

- `--cavity-radius`: override the automatically estimated inner radius, in
  Angstrom.
- `--margin`: additional margin added around the estimated cavity boundary
  (default: 2.0 Angstrom).
- `--lipid-fragments-json`: path to a `lipid_fragments.json` file overriding
  the bundled lipid fragment list.
- `--json`: write machine-readable JSON.

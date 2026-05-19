# MolPACD Examples

## Analyze A Structure

```bash
molpacd analyze input.pdb --axis z --json
```

## Add Caps

```bash
molpacd add input.pdb -o capped.pdb --axis z --resname DUM
```

Use independent radii when the two openings should be capped with their own
inferred sizes:

```bash
molpacd add input.pdb -o capped.pdb --axis z --independent-radius
```

## Remove MolPACD Caps

```bash
molpacd remove capped.pdb -o decapped.pdb --json
```

Removal uses MolPACD metadata when present. If metadata is absent, or if you
intentionally override metadata values, pass `--force` with explicit matching
criteria.

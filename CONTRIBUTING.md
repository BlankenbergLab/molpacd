# Contributing

## Development

Create a local Conda environment and run the full check suite:

```bash
make env
make check
```

Run a narrower check while iterating:

```bash
make test
make lint
make mypy
```

Install and serve the documentation locally:

```bash
make docs-deps
make docs-serve
```

## Expectations

- Keep generated environments, caches, and build outputs out of version control.
- Add tests for behavior changes, especially structure parsing, cap metadata, CLI output, and removal safety.
- Prefer small deterministic fixtures for geometry edge cases. Real structure fixtures should include provenance in `tests/data/README.md`.
- Keep documentation updates in sync with CLI and API behavior.
- Run formatting and type checks before opening a pull request.

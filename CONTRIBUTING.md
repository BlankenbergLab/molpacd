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

## Release

Releases are published from GitHub Actions using PyPI Trusted Publishing. Before
the first release, configure trusted publishers for the `testpypi` and `pypi`
GitHub environments in the `BlankenbergLab/molpacd` project on TestPyPI and
PyPI.

For each release:

1. Update the version in `pyproject.toml`, `src/molpacd/_version.py`, and
   `CITATION.cff`.
2. Move relevant `CHANGELOG.md` entries from `Unreleased` to the release
   version.
3. Run `python -m nox -r` and `make docs`.
4. Commit the release changes and push them to the GitHub default branch.
5. Run the `Publish` workflow manually to publish the distribution to TestPyPI.
6. Smoke test the TestPyPI package in a fresh environment:

```bash
python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  molpacd==0.1.0
```

7. Create and push an annotated version tag to publish to PyPI:

```bash
git tag -a v0.1.0 -m "molpacd 0.1.0"
git push origin v0.1.0
```

## Expectations

- Keep generated environments, caches, and build outputs out of version control.
- Add tests for behavior changes, especially structure parsing, cap metadata, CLI output, and removal safety.
- Prefer small deterministic fixtures for geometry edge cases. Real structure fixtures should include provenance in `tests/data/README.md`.
- Keep documentation updates in sync with CLI and API behavior.
- Run formatting and type checks before opening a pull request.

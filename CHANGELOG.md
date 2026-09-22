# Changelog

## Unreleased

- Added `remove-lipids` (`remove_cavity_lipids`) to detect and remove lipids
  packed inside a protein cavity (e.g. a beta barrel interior), ported from
  the standalone MemGen `remove_cavity_lipids.py` script. Lipid residues are
  recognized from a bundled `lipid_fragments.json`, optionally overridden
  with `--lipid-fragments-json`.

## 0.2.0 - 2026-07-29

- Added configurable cap inversion through the Python API and command-line interface.
- Updated the GitHub checkout action to its Node.js 24-based release.

## 0.1.0 - 2026-06-26

- Added safer cap-removal metadata with generated atom and residue ranges.
- Added stricter option validation for analysis and cap generation.
- Added explicit rejection of multi-model structures for analysis and cap addition.
- Added JSON output options for cap addition and removal.
- Added PyPI and TestPyPI publishing workflow using GitHub Trusted Publishing.
- Added PyPI installation instructions to the README and documentation.
- Expanded automated test coverage for release readiness.
- Added developer workflow, citation, contribution, and fixture provenance docs.

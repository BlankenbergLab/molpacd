from __future__ import annotations

from pathlib import Path

import nox

nox.options.sessions = ["lint", "format", "mypy", "tests", "build"]
nox.options.error_on_missing_interpreters = False
nox.options.envdir = ".nox"

PYTHONS = ["3.9", "3.10", "3.11", "3.12", "3.13"]
PIP_CACHE_DIR = Path(".cache/pip").resolve()


def _use_local_pip_cache(session: nox.Session) -> None:
    session.env["PIP_CACHE_DIR"] = str(PIP_CACHE_DIR)


@nox.session(python=PYTHONS)
def tests(session: nox.Session) -> None:
    _use_local_pip_cache(session)
    session.install("-e", ".[dev]")
    session.run("pytest", "--cov=molpacd", "--cov-report=term-missing", *session.posargs)


@nox.session(python="3.13")
def lint(session: nox.Session) -> None:
    _use_local_pip_cache(session)
    session.install("ruff")
    session.run("ruff", "check", ".")


@nox.session(python="3.13")
def format(session: nox.Session) -> None:
    _use_local_pip_cache(session)
    session.install("ruff")
    session.run("ruff", "format", "--check", ".")


@nox.session(python="3.13")
def mypy(session: nox.Session) -> None:
    _use_local_pip_cache(session)
    session.install("-e", ".[dev]")
    session.run("mypy")


@nox.session(python="3.13")
def build(session: nox.Session) -> None:
    _use_local_pip_cache(session)
    session.install("build", "setuptools", "twine", "wheel")
    session.run("python", "-m", "build")
    dists = sorted(str(path) for path in Path("dist").glob("*"))
    session.run("twine", "check", *dists)

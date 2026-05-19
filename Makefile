CONDA ?= conda
ENV_PREFIX ?= .conda/molpacd-dev
PY ?= 3.13
NOX = $(CONDA) run -p $(ENV_PREFIX) python -m nox

.PHONY: help env env-update install test test-all lint format-check mypy build check clean

help:
	@printf '%s\n' \
		'Available targets:' \
		'  env           Create the local development conda environment' \
		'  env-update    Update and prune the local development environment' \
		'  install       Install the package with development dependencies' \
		'  test          Run tests for PY, defaulting to 3.13' \
		'  test-all      Run tests across all configured nox Python versions' \
		'  lint          Run ruff lint checks' \
		'  format-check  Run ruff formatting checks' \
		'  mypy          Run static type checks' \
		'  build         Build and validate distributions' \
		'  check         Run lint, format, mypy, tests, and build' \
		'  clean         Remove local build and test artifacts'

env:
	$(CONDA) env create -p $(ENV_PREFIX) -f environment-dev.yml

env-update:
	$(CONDA) env update -p $(ENV_PREFIX) -f environment-dev.yml --prune

install:
	$(CONDA) run -p $(ENV_PREFIX) python -m pip install -e ".[dev]"

test:
	$(NOX) -s tests-$(PY)

test-all:
	$(NOX) -s tests

lint:
	$(NOX) -s lint

format-check:
	$(NOX) -s format

mypy:
	$(NOX) -s mypy

build:
	$(NOX) -s build

check:
	$(NOX) -s lint format mypy tests-$(PY) build

clean:
	rm -rf .cache/pip .nox build dist src/*.egg-info

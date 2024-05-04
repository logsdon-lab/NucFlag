.PHONY: build test venv clean dev install upload

BIN=venv/bin/
PROJECT_NAME=nucflag

test:
	$(BIN)python3 -m pip install pytest
	$(BIN)python3 -m pytest -vv

build:
	$(BIN)python3 -m pip install --upgrade build
	$(BIN)python3 -m build

install:
	$(BIN)python3 -m pip uninstall -y $(PROJECT_NAME)
	$(BIN)python3 -m pip install $(shell find dist -name "*.whl" | sort -r | head -1)

dev:
	$(BIN)python3 -m pip install -r requirements-dev.txt

venv:
	python3 -m virtualenv venv

clean:
	rm -rf dist/ venv/ .*cache/ *.egg-info/ benchmarks/ logs/ output/ .snakemake/

upload:
	$(BIN)python3 -m pip install --upgrade twine
	$(BIN)python3 -m twine upload --repository pypi dist/*

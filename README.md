# Computationally-assisted synthesis planning for recyclable-by-design polymers

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/p7k/bottle/GH6_binder_deploy_poetry)

## Description
Runs network expansions using Pickaxe, finds paths, post-processes, and renders spreadsheet + pdf deliverables.

## Installation. 

### Use [poetry](https://python-poetry.org/docs/) to install dependencies

It's better to install `poetry` globally, using something like `brew` or `mise`.

```shell
poetry install
```

Installing poetry without admin priveleges. Assumes correct version of python install, e.g., 3.12

```shell
curl https://bootstrap.pypa.io/get-pip.py | python3.12
python3.12 -m pip install --user pipx
python3.12 -m pipx ensurepath
pipx install --python /usr/bin/python3.12 poetry
```

### Configure filepaths

This script will create a yaml file with filepaths importable through src.config. Note: should either place data and results directories outside project directory or gitignore them to avoid issues with Github file size limit.

```shell
python intialize_configs.py /path/to/data /path/to/results
```

### Download equilibrator database

```python
from equilibrator_assets.local_compound_cache import LocalCompoundCache
lc = LocalCompoundCache()
lc.generate_local_cache_from_default_zenodo('compounds.sqlite')
```

### Convert equilibrator database from sqlite to postgres by following sqlite_to_postgres.ipynb

### Acquire Chemaxon license 

License must cover pKa calculator and ```cxcalc``` must exist in PATH. See more from [equilibrator docs](https://equilibrator.readthedocs.io/en/latest/local_cache.html)

## Directories
1. artifacts (ancillary small data files)
2. logs
3. scripts
4. notebooks
5. src (packages and modules)
6. tests

## Usage
### Running a network expansion with Pickaxe

```
python scripts/expand.py rules=mechinformed_rules_w_coreactants starters=ccm_aa targets=bottle_targets_24
```

See conf/expand.yaml for more options

### Processing network expansion

Expansion processing has three main parts

1. Path finding

TODO

2. Mechanism processing

Analyzing the likelihood that an existing or engineered enzyme can significantly catalyze the predicted reactions. The main thrust of this is to compare reactant structures to those of characterized enzymatic reactions from [Rhea](https://www.rhea-db.org/) and [UniProt](https://www.uniprot.org/), extracted and formatted according to schemas in [this repository](https://github.com/stefanpate/enz-rxn-data).

3. Thermodynamic calculations

TODO: calculates thermodynamic values for predicted pathways using eQuilibrator (4) stores processed paths, predicted reactions, and known reaction information to json

### Visualizing processed synthesis paths

## Dev Notes

1. you will need the development requirements

```shell
poetry install --with dev
```

2. running unit tests with coverage

```shell
poe test
```

## Building the docker of viola PathViewer to Railway

### Pre-requisites
- Docker installed and running on your machine
- Basic knowledge of Docker
  - https://docs.docker.com/get-started/
  - https://docs.docker.com/get-started/docker-overview/
- poetry installed
  - https://python-poetry.org/docs/
- poe-the-poet plugin for poetry installed
  - e.g. `poetry self add 'poethepoet[poetry_plugin]'`
- Familiarity with environment variables
  - `BOTTLE_EXPANSION_ASSETS` set and pointing to the directory containing the expansion files
    - e.g. `export BOTTLE_EXPANSION_ASSETS=/path/to/expansion/assets/2024/11/15`
    - or defined and activated in your `.env` file (but that's just your convenience)

```shell
poetry poe build-image
```
if want to pass an assets directory bespoke to this particular build, use this:
```shell
BOTTLE_EXPANSION_ASSETS=/path/to/expansion/assets/2024/11/15 poetry poe build-image
```
## Deploying the docker image to Docker Hub

### Pre-requisites
- Login to Docker Hub (https://hub.docker.com/u/synbiorox) ( you can guess the password - ask anyone at TyoLab)
```shell
docker login -u synbiorox -p <password>
```

- Poetry, poe (see Docker Build Above)

```shell
poetry poe push-image
```
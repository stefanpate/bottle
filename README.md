# Computationally-assisted synthesis planning for recyclable-by-design polymers

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/p7k/bottle/GH6_binder_deploy_poetry)

## Description
Runs network expansions using Pickaxe, finds paths, post-processes, and renders spreadsheet + pdf deliverables.

## Installation. 

### Using `poetry` https://python-poetry.org/docs/

It's better to install `poetry` globally, using something like `brew` or `mise`.

```shell
poetry install
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
1. artifacts (stores results)
2. data (input data: known reactions, mappings, starters, targets, rules)
3. logs
4. scripts
5. notebooks
6. src (packages and modules)
7. tests

## Usage
- run_pickaxe.py does a network expansion from starters to targets.
- process_expansion.py from Pickaxe output, (1) finds paths between starter and target molecules (2) assigns known reaction analogues to the predicted reactions in found paths using a reaction-based similarity score (3) calculates thermodynamic values for predicted pathways using eQuilibrator (4) stores processed paths, predicted reactions, and known reaction information to json
- map_operators.py maps operators to reaction with option to return the reaction center
- sprhea_extract.py does most of the processing of known protein-reaction data from SwissProt + Rhea flat files
- download_rhea_smi_name_pairs.py pulls molecule SMILES-name pairs from Rhea

## Dev Notes

1. you will need the development requirements

```shell
poetry install --with dev
```

2. running unit tests with coverage

```shell
poe test
```

## Scripts (in order)

1. `bottle` a master command, supports `--help`
2. `bottle filter-smiles` allows filtering of input streams of smiles lines (files or stdin) by smarts patterns. Example
   usage:

```shell
bottle filter-smiles --mol-smarts=C4_FG123 ~/Downloads/known_compounds_240310_v2.txt -
```

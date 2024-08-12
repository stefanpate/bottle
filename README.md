# Computationally-assisted synthesis planning for recyclable-by-design polymers

## Description
Runs network expansions using Pickaxe, finds paths, post-processes, and renders spreadsheet + pdf deliverables.

## Installation. 

### Using `poetry` https://python-poetry.org/docs/

It's better to install `poetry` globally, using something like `brew` or `mise`.

```shell
poetry install
```

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
- prune_pickaxe.py takes a pickaxe output, finds paths from starters to targets, saves a pruned pickaxe object and a ProcessedExpansion object
- [Placeholder for thermodynamics. Currently calculations done outside this directory]
- analyze_paths.py does reaction-center-maximum-common-substructure calculations and enters thermo info into ProcessedExpansion object
- vis_pathways.ipynb filters & sorts paths in ProcessedExpansion object and renders pdfs of reaction formula and spreadsheets w/ further info

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

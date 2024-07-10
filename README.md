# Computationally-assisted synthesis planning for recyclable-by-design polymers

## Installation.

### Using `poetry` https://python-poetry.org/docs/

It's better to install `poetry` globally, using something like `brew` or `mise`.

```shell
poetry install
```

## Directories

2. data
3. scripts
4. src (packages and modules)
5. tests
6. artifacts (stores results)

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

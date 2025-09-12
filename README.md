# Computationally-assisted synthesis planning for recyclable-by-design polymers

## Description
Runs network expansions using Pickaxe, finds paths, post-processes, and renders spreadsheet + pdf deliverables.

## Installation

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

### Acquire Chemaxon license 

License must cover pKa calculator and ```cxcalc``` must exist in PATH. See more from [equilibrator docs](https://equilibrator.readthedocs.io/en/latest/local_cache.html)

## Directories
1. artifacts (ancillary small data files)
2. logs
3. scripts
4. notebooks
5. src (packages and modules)
6. pathway_viewer
7. tests

## Usage
### Running a network expansion with Pickaxe

```
python scripts/expand.py rules=mechinformed_rules_w_coreactants starters=ccm_aa targets=bottle_targets_24
```

See conf/expand.yaml for more options

### Processing network expansion

0. Parse expansion

Outputs atom-mapped smarts in the `src/schemas.expansion_reactions_schema` schema. The script `parse_expansion.py` can handle forward expansions, retro expansion, or a combination of the two.

Example usage:
```
python parse_expansion.py fwd_exp=my_forward_expansion rev_expansion=my_reverse_expansion
```

1. Path finding

There are two separate files available, ```linear_pathfinding.py```and ```retrosynthesis.py```. Both require the location of the parsed expansion, the output location, and maximum depth desired. Each has several more optional arguments.

Example usage:
```
python linear_pathfinding.py expansion=my/expansion/location casp_study=my/casp/study max_depth=5
```

```retrosynthesis.py``` utilizes a network of known enzymatic reactions in addition to predicted reactions from the expansion to find synthetic paths.

2. Reaction structural analysis

Looks up known reactions similar to the predicted reactions and collects their enzymes. Feasibility of the predicted reactions is scored with a [machine learning model](https://github.com/tyo-nu/DORA_XGB).

Example usage:
```
python analyze_structures.py casp_study=my/casp/study
```

3. Thermodynamic calculations (Optional)

Calculates thermodynamic values (e.g., $\Delta G$) for predicted pathways using [eQuilibrator](https://equilibrator.readthedocs.io/en/latest/local_cache.html).

Example usage:
```
python analyze_thermo.py casp_study=my/casp/study
```

Example usage:
```
python analyze_thermo.py casp_study=my/casp/study
```

4. Reaction drawing

Draws reactions to svg

Example usage:
```
python draw_reactions.py casp_study=my/casp/study
```

### Visualizing processed synthesis paths

Processed paths can be visualized in a path viewer app. Locally, this can be launched running

```
streamlit run path_viewer/list_view.py
```
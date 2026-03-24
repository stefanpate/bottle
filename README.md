# Biological computer aided synthesis planning

## Overview

This is a computational pipeline designed for biological synthesis planning. The system performs network expansions using Pickaxe, discovers synthetic pathways, processes results, and visualizes results in an interactive app.

### Key Features

- 🧪 **Network Expansion**: Automated chemical reaction network generation using Pickaxe
- 🔍 **Pathway Discovery**: Advanced algorithms for finding optimal synthetic routes
- 📊 **Structural Analysis**: Machine learning-based feasibility scoring for predicted reactions
- ⚡ **Thermodynamic Calculations**: Integration with eQuilibrator for thermodynamic analysis
- 📈 **Visualization**: Interactive pathway viewer and reaction drawing tools
- 📋 **Comprehensive Reporting**: Automated generation of analysis reports and deliverables

## 🚀 Installation

### Prerequisites

- Python 3.12 or higher
- [Poetry](https://python-poetry.org/docs/) for dependency management

### Using Poetry (Recommended)

It's recommended to install Poetry globally using a package manager like `brew` or `mise`.

```bash
poetry install
```

### Installing Poetry Without Admin Privileges

If you don't have administrative privileges, you can install Poetry using the following method. This assumes Python 3.12 is correctly installed:

```bash
curl https://bootstrap.pypa.io/get-pip.py | python3.12
python3.12 -m pip install --user pipx
python3.12 -m pipx ensurepath
pipx install --python /usr/bin/python3.12 poetry
```

### ⚙️ Configuration Setup

#### 1. Configure File Paths

Create a config file at `conf/filepaths/filepaths.yaml` following the template in `FILEPATHS_TEMPLATE.yaml`.

#### 2. Download Equilibrator Database

```python
from equilibrator_assets.local_compound_cache import LocalCompoundCache
lc = LocalCompoundCache()
lc.generate_local_cache_from_default_zenodo('compounds.sqlite')
```

## 📁 Project Structure

1. **artifacts** - Ancillary small data files
2. **logs** - Log files and error tracking
3. **scripts** - Main executable scripts
4. **notebooks** - Jupyter notebooks for analysis
5. **src** - Source packages and modules
6. **path_viewer** - Visualization tools
7. **tests** - Test suite

## 🔧 Usage

### Running a Network Expansion with Pickaxe

```bash
python scripts/expand.py rules=mechinformed_rules_w_coreactants starters=ccm_aa targets=bottle_targets_24
```

See `conf/expand.yaml` for more configuration options.

### Processing Network Expansion

The pipeline consists of several sequential steps for processing network expansions:

**1. Path finding**

There are two separate files available: `linear_pathfinding.py` and `retrosynthesis.py`. Both require the location of the parsed expansion, the output location, and maximum depth desired. Each has several additional optional arguments.

**Example usage:**
```bash
python linear_pathfinding.py expansion=my/expansion/location casp_study=my/casp/study max_depth=5
```

> **Note:** `retrosynthesis.py` utilizes a network of known enzymatic reactions in addition to predicted reactions from the expansion to find synthetic paths.

**2. Reaction structural analysis**

Looks up known reactions similar to the predicted reactions and collects their enzymes. Feasibility of the predicted reactions is scored with a [machine learning model](https://github.com/tyo-nu/DORA_XGB).

**Example usage:**
```bash
python analyze_structures.py casp_study=my/casp/study
```

**3. Thermodynamic calculations (Optional)**

Calculates thermodynamic values (e.g., ΔG) for predicted pathways using [eQuilibrator](https://equilibrator.readthedocs.io/en/latest/local_cache.html).

**Example usage:**
```bash
python analyze_thermo.py casp_study=my/casp/study
```

**4. Reaction drawing**

Draws reactions to SVG format.

**Example usage:**
```bash
python draw_reactions.py casp_study=my/casp/study
```

### 🖥️ Visualizing Processed Synthesis Paths

Processed paths can be visualized using the interactive path viewer app. To launch it locally:

```bash
streamlit run path_viewer/list_view.py -- --casp-study my_casp_study
```

Where `my_casp_study` is the name of the subdirectory in your processed data location that contains the processed parquet and SVG files.

## 🔧 Using Other Network Expansion Software

`parse_expansion.py` is written to extract from a Pickaxe-generated network expansion. Users may use different network expansion software upstream of the BOTTLE post-processing pipeline, provided they implement a script analogous to `parse_expansion.py` that outputs data conforming to the three interim data schemas defined in `src/schemas.py`.

Each schema specifies required columns that must be present in the output files. The expected columns for each schema are:

### 📋 Data Schema Requirements

#### 1. Expansion Reactions Schema

| Column Name       | Description                                                |
|-------------------|------------------------------------------------------------|
| `smarts`          | Reaction SMARTS representation                             |
| `am_smarts`       | Atom-mapped reaction SMARTS                                |
| `rules`           | List of rule names/IDs that generated the reaction        |
| `half_expansion`  | Type of expansion: `'forward'` or `'retro'`               |
| `size`            | Number of atoms specified by the rule(s)                  |

#### 2. Compounds Schema

| Column Name | Description                                                     |
|-------------|-----------------------------------------------------------------|
| `id`        | Unique identifier for the compound                              |
| `smiles`    | SMILES string representation of the compound                   |
| `type`      | Compound type: `'source'`, `'target'`, `'known'`, `'helper'`, or `'checkpoint'` |
| `name`      | Compound name (optional)                                       |

#### 3. Generations Schema

| Column Name       | Description                                                      |
|-------------------|------------------------------------------------------------------|
| `half_expansion`  | Type of expansion: `'forward'` or `'retro'`                     |
| `generation`      | Number of generations/steps the half expansion was applied for  |

---